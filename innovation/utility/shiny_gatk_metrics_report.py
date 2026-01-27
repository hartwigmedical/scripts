#!/usr/bin/env python3
"""
Shiny dashboard for WGS metrics (Picard + runlevel TSV + flagstat)

- Bucket selection + mode chooser (SBX batches vs Manual batch creation)
- Manual batch builder (create multiple batches; selected samples disappear from pool)
- Unified, normalized meta across both paths (same BatchKey/BatchLabel)
- Stable batch tiles, filters, plots (no list-type crashes)
- Yield (Gb) computed from Total_Reads × mean_read_length / 1e9
"""

from __future__ import annotations

# ========= Standard lib =========
import ast
import datetime as dt
import hashlib
import re
import sys
from io import StringIO
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

# ========= Third party =========
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import yaml
from gcsfs import GCSFileSystem
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_widget

# ========= Config =========
DEFAULT_BUCKET = "gs://gatk-wgs-metrics-clean"

METRICS_TO_PLOT: List[str] = [
    "GENOME_TERRITORY", "MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE", "MAD_COVERAGE",
    "PCT_EXC_ADAPTER", "PCT_EXC_MAPQ", "PCT_EXC_DUPE", "PCT_EXC_UNPAIRED", "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP", "PCT_EXC_CAPPED", "PCT_EXC_TOTAL", "PCT_1X", "PCT_5X", "PCT_10X",
    "PCT_15X", "PCT_20X", "PCT_25X", "PCT_30X", "PCT_40X", "PCT_50X", "PCT_60X", "PCT_70X",
    "PCT_80X", "PCT_90X", "PCT_100X", "FOLD_80_BASE_PENALTY", "FOLD_90_BASE_PENALTY",
    "FOLD_95_BASE_PENALTY", "HET_SNP_SENSITIVITY", "HET_SNP_Q", "AT_DROPOUT", "GC_DROPOUT",
    "Phred_All", "Phred_Substitutions", "Phred_Insertions", "Phred_Deletions", "Phred_Indel",
    "median_read_length", "mean_read_length", "ratio_90_to_10percentile",
    "MEDIAN_INSERT_SIZE", "MEDIAN_ABSOLUTE_DEVIATION",
    "Num_Full_Length_Reads", "Num_One_Plus_Reads", "READ_LENGTH_N50", "TOTAL_READS",
    "mmetrics_median_read_length", "PF_ALIGNED_BASES", "FLAGSTAT_DUPLICATES", "Yield"
]

BASE_PALETTE = [
    "#e41a1c", "#ff7f00", "#377eb8", "#4daf4a", "#984ea3",
    "#a65628", "#f781bf", "#999999", "#dede00", "#66c2a5",
    "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
    "#e5c494", "#b3b3b3",
]

THRESHOLDS: Dict[str, Dict[str, Any]] = {
    "MEAN_COVERAGE": {"lines": [
        {"value": 30, "color": "orange", "label": "Reference Min (30x)"},
        {"value": 90, "color": "red",    "label": "Tumor Min (90x)"},
    ]},
    "MEDIAN_COVERAGE": {"lines": [
        {"value": 30, "color": "orange", "label": "Reference Min (30x)"},
        {"value": 90, "color": "red",    "label": "Tumor Min (90x)"},
    ]},
    "PCT_30X": {"lines": [{"value": 0.95, "color": "red"}], "format": ".0%"},
    "AT_DROPOUT": {"lines": [], "format": ".1%"},
    "GC_DROPOUT": {"lines": [], "format": ".1%"},
    "percent_duplex":  {"lines": [], "format": ".1%"},
    "percent_oneplus": {"lines": [], "format": ".1%"},
}

PCT_LIKE_COLUMNS = ["AT_DROPOUT", "GC_DROPOUT"]

DOWNLOAD_FILENAME_FMT = "table_overview_{:%Y%m%d_%H%M%S}.csv"

# ========= Debug =========
def dbg(tag: str, **kv: Any) -> None:
    ts = dt.datetime.now().strftime("%H:%M:%S")
    parts = []
    for k, v in kv.items():
        if isinstance(v, (list, tuple)) and len(v) > 10:
            parts.append(f"{k}=[{', '.join(map(repr, v[:10]))}, … +{len(v)-10}]")
        else:
            parts.append(f"{k}={repr(v)}")
    print(f"[{ts}] {tag} " + " ".join(parts), file=sys.stderr, flush=True)

def _head(seq, n=20):
    try:
        lst = list(seq)
        return lst[:n]
    except Exception:
        return []

# ========= Helpers: scalarization / normalization =========
def _scalarize(x):
    if isinstance(x, (list, tuple)):
        return x[0] if x else None
    if isinstance(x, set):
        for it in x:
            return it
        return None
    if isinstance(x, dict):
        return None
    return x

def _maybe_unlist_string(s):
    if not isinstance(s, str):
        return s
    t = s.strip()
    if t.startswith("[") and t.endswith("]"):
        try:
            val = ast.literal_eval(t)
            if isinstance(val, (list, tuple)) and val:
                return str(val[0])
        except Exception:
            pass
        t2 = t.strip("[]").strip().strip("'").strip('"')
        return t2
    return s

def _clean_str(x):
    x = _scalarize(x)
    x = _maybe_unlist_string(x) if isinstance(x, str) else x
    if x is None:
        return ""
    return str(x).strip()

def fmt_date_ddmmyyyy(iso_date: Optional[str]) -> str:
    if not iso_date:
        return "unknown"
    try:
        return dt.date.fromisoformat(iso_date).strftime("%d-%m-%Y")
    except Exception:
        return "unknown"

def run_from_cycle(cycle: Optional[str]) -> str:
    if not cycle:
        return "Run unknown"
    m = re.search(r"(\d+)", str(cycle))
    return f"Run {m.group(1)}" if m else "Run unknown"

def is_illumina_sample_id(sample_id: str) -> bool:
    return isinstance(sample_id, str) and ("nseq6" in sample_id.lower())

def platform_from_sid(sid: str) -> str:
    return "Illumina" if is_illumina_sample_id(sid) else "SBX"

def batch_key(platform: str, date_iso: Optional[str], cycle: Optional[str]) -> str:
    return f"{platform}|{date_iso or 'unknown'}|{cycle or 'unknown'}"

def batch_label(platform: str, date_iso: Optional[str], cycle: Optional[str]) -> str:
    return f"{platform} | Batch {fmt_date_ddmmyyyy(date_iso)} | {run_from_cycle(cycle)}"

# ========= GCS I/O =========
def read_text_if_exists(gcs: GCSFileSystem, path: str) -> Optional[str]:
    try:
        if gcs.exists(path):
            with gcs.open(path, "r") as f:
                return f.read()
    except Exception:
        pass
    return None

def find_sample_folders(gcs: GCSFileSystem, bucket_path: str) -> List[str]:
    try:
        all_paths = gcs.ls(bucket_path, detail=False)
        return [p.rsplit("/", 1)[-1] for p in all_paths if gcs.isdir(p)]
    except Exception:
        return []

def find_flagstat_txt(gcs: GCSFileSystem, base_bucket: str, sample_id: str) -> Optional[str]:
    try:
        paths = gcs.glob(f"{base_bucket.rstrip('/')}/{sample_id}/run-flagstat/*.txt")
        return paths[0] if paths else None
    except Exception:
        return None

# ========= YAML + parsing =========
def parse_yaml_field(text: Optional[str], dotted_key: str, fallback: Optional[str] = None) -> Optional[str]:
    if text is None:
        return fallback
    try:
        cur: Any = yaml.safe_load(text) or {}
        for part in dotted_key.split("."):
            if not isinstance(cur, dict) or part not in cur:
                return fallback
            cur = cur[part]
        return str(cur)
    except Exception:
        pass
    tail = dotted_key.split(".")[-1]
    for line in text.splitlines():
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        if k.strip() == tail:
            return v.strip().strip('"').strip("'")
    return fallback

def get_any_seq_uri(text: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    if text is None:
        return None, None
    candidate_keys = [
        "params.bam_uri", "params.cram_uri", "params.bam_file", "params.cram_file",
        "inputs.bam_uri", "inputs.cram_uri", "inputs.bam_file", "inputs.cram_file",
        "bam_uri", "cram_uri", "bam_file", "cram_file",
    ]
    for key in candidate_keys:
        v = parse_yaml_field(text, key, fallback=None)
        if isinstance(v, str) and v.strip():
            return v.strip(), key
    return None, None

def extract_date_cycle_from_bam_uri(bam_uri: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    if not bam_uri:
        return None, None
    s = str(bam_uri); sl = s.lower()
    iso = None
    m = re.search(r"\b(\d{4})-(\d{2})-(\d{2})\b", sl)
    if m:
        y, mo, d = map(int, m.groups())
        try: iso = dt.date(y, mo, d).isoformat()
        except ValueError: pass
    if not iso:
        m = re.search(r"\b(\d{2})-(\d{2})-(\d{4})\b", sl)
        if m:
            d, mo, y = map(int, m.groups())
            try: iso = dt.date(y, mo, d).isoformat()
            except ValueError: pass
    if not iso:
        for seg in s.split("/"):
            m = re.search(r"(\d{6})", seg)
            if not m: continue
            yy, mo, d = int(m.group(1)[:2]), int(m.group(1)[2:4]), int(m.group(1)[4:6])
            try:
                iso = dt.date(2000 + yy, mo, d).isoformat()
                break
            except ValueError:
                continue
    m_cycle = re.search(r"(cycle\s*\d{1,3})", sl, flags=re.I)
    cycle = re.sub(r"\s+", "", m_cycle.group(1)) if m_cycle else None
    return iso, cycle

# ========= Picard / Runlevel parsing =========
def parse_picard_single_metrics(gcs: GCSFileSystem, file_path: str) -> Optional[pd.Series]:
    try:
        with gcs.open(file_path, "r") as f:
            content = f.read()
        m = content.find("## METRICS CLASS")
        if m == -1:
            return None
        block = content[m:]
        df = pd.read_csv(StringIO(block), sep="\t", comment="#", nrows=1)
        if df.empty:
            return None
        df.columns = df.columns.astype(str).str.strip()
        s = df.iloc[0]; s.index = s.index.astype(str).str.strip()
        return s
    except Exception:
        return None

def parse_alignment_summary_histogram(gcs: GCSFileSystem, file_path: str) -> Optional[pd.DataFrame]:
    try:
        with gcs.open(file_path, "r") as f:
            content = f.read()
        m = re.search(r"^##\s*HISTOGRAM\s+java\.lang\.(?:Integer|Double)\s*$", content, flags=re.M)
        if not m:
            return None
        tail = content[m.end():].lstrip("\n")
        next_section = re.search(r"^\#\#\s", tail, flags=re.M)
        hist_text = tail[: next_section.start()] if next_section else tail
        df = pd.read_csv(StringIO(hist_text), sep="\t", comment="#", header=0, engine="python")
        if df is None or df.empty:
            return None
        # detect cols
        read_len_col = next((c for c in df.columns if c.strip().upper() == "READ_LENGTH"), None)
        if not read_len_col:
            read_len_col = next((c for c in df.columns if ("read" in c.lower() and "length" in c.lower())), None)
        if not read_len_col:
            return None
        def find_col(*need: str) -> Optional[str]:
            for c in df.columns:
                lc = c.lower()
                if all(s in lc for s in need):
                    return c
            return None
        total_col   = find_col("total", "length", "count") or find_col("total", "reads") or find_col("total")
        aligned_col = find_col("aligned", "length", "count") or find_col("aligned", "reads") or find_col("aligned")
        out = pd.DataFrame()
        out["READ_LENGTH"] = pd.to_numeric(df[read_len_col], errors="coerce")
        if total_col and total_col in df:   out["TOTAL"]   = pd.to_numeric(df[total_col], errors="coerce")
        if aligned_col and aligned_col in df: out["ALIGNED"] = pd.to_numeric(df[aligned_col], errors="coerce")
        out = out.dropna(subset=["READ_LENGTH"])
        return out if not out.empty else None
    except Exception:
        return None

def parse_alignment_summary_table(gcs: GCSFileSystem, file_path: str) -> Optional[pd.Series]:
    try:
        with gcs.open(file_path, "r") as f:
            content = f.read()
        m = re.search(r"^##\s+METRICS\s+CLASS\s+picard\.analysis\.AlignmentSummaryMetrics\s*$", content, flags=re.M)
        if not m:
            return None
        tail = content[m.end():]
        next_section = re.search(r"^\#\#\s", tail, flags=re.M)
        block = tail[: next_section.start()] if next_section else tail
        df = pd.read_csv(StringIO(block), sep="\t", comment="#", engine="python")
        if df is None or df.empty:
            return None
        df.columns = df.columns.astype(str).str.strip()
        row = None
        if "CATEGORY" in df.columns:
            if (df["CATEGORY"] == "PAIR").any():
                row = df.loc[df["CATEGORY"] == "PAIR"].iloc[0]
            elif (df["CATEGORY"] == "UNPAIRED").any():
                row = df.loc[df["CATEGORY"] == "UNPAIRED"].iloc[0]
            else:
                row = df.iloc[0]
        else:
            row = df.iloc[0]
        s = row.copy(); s.index = s.index.astype(str).str.strip()
        return s
    except Exception:
        return None

def parse_runlevel_for_sample(gcs: GCSFileSystem, tsv_path: str, sample_id: str) -> Optional[Dict[str, float]]:
    try:
        with gcs.open(tsv_path, "r") as f:
            text = f.read()
    except Exception as e:
        dbg("runlevel_open_error", path=tsv_path, error=repr(e)); return None
    lines = text.splitlines()
    header_idx = next((i for i, ln in enumerate(lines) if "sample_name" in ln.split("\t")), None)
    if header_idx is None:
        return None
    try:
        df = pd.read_csv(StringIO("\n".join(lines[header_idx:])), sep="\t")
    except Exception as e:
        dbg("runlevel_read_error", path=tsv_path, error=repr(e)); return None
    if df.empty or "sample_name" not in df.columns:
        return None
    df = df.drop(columns=[c for c in [
        "run_id","Analysis Start","Analysis End","ppa2bam Start","ppa2bam End",
        "Sort Start","Sort End","BAM Start","BAM End","Diablo Start","Diablo End",
        "Time to Analyze (m)","pct_files_analyzed","SID_Number",
    ] if c in df.columns], errors="ignore")
    df = df[~df["sample_name"].isna()].copy()
    if df.empty:
        return None
    def root_and_is_ref(name: str) -> Tuple[str, bool]:
        if not name or str(name).lower() == "nan":
            return ("", False)
        s = re.sub(r"-(sbx|nseq6)$", "", str(name).strip(), flags=re.I)
        tokens = [t for t in re.split(r"[-_]", s) if t]
        root = next((t.upper() for t in tokens if re.match(r"(?i)^h\d{5,}$", t)), "")
        if not root:
            m = re.search(r"[A-Za-z0-9]+", s)
            root = m.group(0).upper() if m else ""
        is_ref = any((t.lower() == "ref") or re.match(r"(?i)^R(?:\d+|I+)?$", t) for t in tokens)
        return (root, is_ref)
    df["__name__"] = df["sample_name"].astype(str)
    root_flag = df["__name__"].map(root_and_is_ref)
    df["__root__"] = root_flag.map(lambda x: x[0])
    df["__is_ref__"] = root_flag.map(lambda x: x[1])
    target_root, target_is_ref = root_and_is_ref(sample_id)
    sub = df[(df["__root__"] == target_root) & (df["__is_ref__"] == target_is_ref)]
    if sub.empty and target_root:
        base = re.escape(target_root)
        if target_is_ref:
            pat = rf"(?i)^{base}(?:[-_]ref|[-_]R(?:\d+|I+)?)(?:[-_].*)?$"
        else:
            pat = rf"(?i)^{base}(?:(?:[-_]T(?:\d+|I+)?)?(?:[-_].*)?)$"
        df_pat = df[df["__name__"].str.match(pat, na=False)]
        sub = df_pat[df_pat["__is_ref__"] == target_is_ref]
    if sub.empty:
        sub = df[df["__name__"] == str(sample_id)]
    if sub.empty:
        return None
    row = sub[sub["__name__"] == str(sample_id)].iloc[0] if (sub["__name__"] == str(sample_id)).any() else \
          sub.sort_values(["__root__", "__is_ref__", "__name__"]).iloc[0]
    keep = {
        "Phred_All", "Phred_Substitutions", "Phred_Insertions", "Phred_Deletions", "Phred_Indel",
        "median_read_length", "mean_read_length", "ratio_90_to_10percentile",
        "Num_Full_Length_Reads", "Num_One_Plus_Reads", "Total_Reads",
    }
    out: Dict[str, float] = {}
    for k in keep:
        if k in row:
            out[k] = pd.to_numeric(row[k], errors="coerce")
    return out or None

def _tsv_contains_required_field(gcs: GCSFileSystem, tsv_path: str) -> bool:
    """Check if a TSV file contains the Num_Full_Length_Reads field (case-insensitive)."""
    try:
        with gcs.open(tsv_path, "r") as f:
            # Read only the first few KB to find headers
            content = f.read(8192)
        # Look for the field in headers (first line or any line with 'sample_name')
        for line in content.splitlines():
            lower_line = line.lower()
            if "num_full_length_reads" in lower_line:
                return True
        return False
    except Exception as e:
        dbg("tsv_field_check_error", path=tsv_path, error=repr(e))
        return False


def runlevel_tsv_from_bam_uri(gcs: GCSFileSystem, bam_uri: Optional[str]) -> Optional[str]:
    """
    Find the correct runlevel TSV file based on bam/cram URI path.

    Search strategy:
    1. Look at the same directory level as the bam/cram file
    2. Look one level up from the bam/cram file
    3. Look two levels up from the bam/cram file

    For each level, we:
    - Skip files named 'metrics-report.tsv' (these lack required fields)
    - Validate that the TSV contains 'Num_Full_Length_Reads' field
    """
    if not bam_uri:
        return None

    s = str(bam_uri).rstrip("/")

    # Get the directory containing the bam/cram file
    # The bam_uri typically points to a file, so we need its parent directory
    if s.endswith((".bam", ".cram")):
        base_dir = "/".join(s.split("/")[:-1])
    else:
        # If it doesn't end with .bam/.cram, assume it's already a directory or use as-is
        base_dir = s

    # Build list of directories to search (same level, one up, two up)
    dirs_to_search: List[str] = []

    # Same level as bam/cram
    dirs_to_search.append(base_dir)

    # One level up
    parts = base_dir.rstrip("/").split("/")
    if len(parts) > 1:
        one_up = "/".join(parts[:-1])
        dirs_to_search.append(one_up)

    # Two levels up
    if len(parts) > 2:
        two_up = "/".join(parts[:-2])
        dirs_to_search.append(two_up)

    dbg("runlevel_tsv_search", bam_uri=bam_uri, dirs_to_search=dirs_to_search)

    # Search each directory level
    for search_dir in dirs_to_search:
        try:
            candidates = gcs.glob(f"{search_dir}/*.tsv")
            dbg("runlevel_tsv_candidates", dir=search_dir, candidates=candidates)

            # Filter out metrics-report.tsv and validate remaining files
            for tsv_path in sorted(candidates):
                filename = tsv_path.split("/")[-1].lower()

                # Skip metrics-report.tsv (it lacks required fields)
                if filename == "metrics-report.tsv":
                    dbg("runlevel_tsv_skip_metrics_report", path=tsv_path)
                    continue

                # Validate this TSV has the required Num_Full_Length_Reads field
                if _tsv_contains_required_field(gcs, tsv_path):
                    dbg("runlevel_tsv_found_valid", path=tsv_path)
                    return tsv_path
                else:
                    dbg("runlevel_tsv_missing_field", path=tsv_path)
        except Exception as e:
            dbg("runlevel_tsv_search_error", dir=search_dir, error=repr(e))
            continue

    dbg("runlevel_tsv_not_found", bam_uri=bam_uri)
    return None

# ========= Meta normalization =========
def normalize_meta(df_meta_in: pd.DataFrame) -> pd.DataFrame:
    if df_meta_in is None or df_meta_in.empty:
        return pd.DataFrame(columns=["SampleID","Platform","Date","Cycle","bam_uri","yaml_found",
                                     "DateDisplay","Run","BatchKey","BatchLabel"])
    df = df_meta_in.copy()
    for c in ["SampleID","Platform","Date","Cycle","bam_uri","yaml_found"]:
        if c not in df.columns:
            df[c] = ""
    df["SampleID"]  = df["SampleID"].map(_clean_str)
    df["Platform"]  = df["Platform"].map(_clean_str)
    df["Date"]      = df["Date"].map(_clean_str)
    df["Cycle"]     = df["Cycle"].map(_clean_str)
    df["bam_uri"]   = df["bam_uri"].map(_clean_str)
    df["yaml_found"]= df["yaml_found"].apply(lambda x: bool(x))

    df["Platform"] = df["Platform"].apply(lambda p: "Illumina" if str(p).lower()=="illumina" or "nseq6" in str(p).lower() else "SBX")

    def _to_iso(s):
        s = _clean_str(s)
        if not s:
            return ""
        try:
            return dt.date.fromisoformat(s).isoformat()
        except Exception:
            m = re.match(r"^(\d{2})-(\d{2})-(\d{4})$", s)
            if m:
                d, mo, y = map(int, m.groups())
                try:
                    return dt.date(y, mo, d).isoformat()
                except Exception:
                    return ""
            return s
    df["Date"] = df["Date"].map(_to_iso)

    df["DateDisplay"] = df["Date"].apply(fmt_date_ddmmyyyy)
    df["Run"]         = df["Cycle"].apply(run_from_cycle)

    def _key(r):
        return batch_key(r["Platform"], (r["Date"] if r["Date"] else None), (r["Cycle"] if r["Cycle"] else None))
    def _lab(r):
        return batch_label(r["Platform"], (r["Date"] if r["Date"] else None), (r["Cycle"] if r["Cycle"] else None))
    df["BatchKey"]   = df.apply(_key, axis=1)
    df["BatchLabel"] = df.apply(_lab, axis=1)

    return df[["SampleID","Platform","Date","Cycle","bam_uri","yaml_found",
               "DateDisplay","Run","BatchKey","BatchLabel"]]

def manual_batches_to_meta(batches: List[Dict[str, Any]]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for b in (batches or []):
        platform = str(b.get("platform") or "SBX")
        date_iso = b.get("date_iso") or ""
        comment  = b.get("cycle") or ""     # comment stored in Cycle
        for sid in (b.get("samples") or []):
            rows.append({
                "SampleID":  str(sid),
                "Platform":  platform,
                "Date":      date_iso,
                "Cycle":     comment,
                "bam_uri":   "",
                "yaml_found": False,
            })
    return pd.DataFrame(rows, columns=["SampleID","Platform","Date","Cycle","bam_uri","yaml_found"])

def build_batch_choices(df_meta: pd.DataFrame) -> Dict[str, str]:
    if df_meta is None or df_meta.empty:
        return {}
    uniq = df_meta.drop_duplicates(subset=["BatchKey"])[["BatchKey","BatchLabel","Date"]].copy()
    uniq["_Date"] = pd.to_datetime(uniq["Date"], errors="coerce")
    uniq = uniq.sort_values("_Date", ascending=False, na_position="last").drop(columns=["_Date","Date"])
    return dict(zip(uniq["BatchKey"], uniq["BatchLabel"]))

def build_batch_index(df_meta: pd.DataFrame) -> Dict[str, List[str]]:
    if df_meta is None or df_meta.empty:
        return {}
    out: Dict[str, List[str]] = {}
    for _, r in df_meta.iterrows():
        out.setdefault(r["BatchKey"], []).append(r["SampleID"])
    return out

# ========= Colors / plots =========

def lighten_hex(hex_color: str, factor: float = 0.45) -> str:
    hc = hex_color.lstrip("#")
    if len(hc) != 6:
        return hex_color
    r, g, b = int(hc[0:2], 16), int(hc[2:4], 16), int(hc[4:6], 16)
    r = int(r + (255 - r) * factor)
    g = int(g + (255 - g) * factor)
    b = int(b + (255 - b) * factor)
    return f"#{r:02x}{g:02x}{b:02x}"

def build_two_tone_color_map(df: pd.DataFrame) -> Tuple[Dict[str, str], List[str]]:
    df = df.copy()
    for col in ("Group","Date","Cycle"):
        if col in df.columns:
            df[col] = df[col].map(_clean_str)
    tmp = df[["Group","Date","Cycle"]].drop_duplicates().copy()
    tmp["_Date"] = pd.to_datetime(tmp["Date"], errors="coerce")
    def _cycle_num(c: str) -> int:
        m = re.search(r"(\d+)", str(c) if c else "")
        return int(m.group(1)) if m else -1
    tmp["_CycleNum"] = tmp["Cycle"].map(_cycle_num)
    ordered_groups = tmp.sort_values(["_Date","_CycleNum","Group"], ascending=[True, True, True])["Group"].tolist()
    if "Illumina" in ordered_groups:
        ordered_groups = ["Illumina"] + [g for g in ordered_groups if g != "Illumina"]
    base_color_map = {grp: BASE_PALETTE[i % len(BASE_PALETTE)] for i, grp in enumerate(ordered_groups)}
    two_tone = {f"{grp} — Tumor": base for grp, base in base_color_map.items()}
    two_tone.update({f"{grp} — Reference": lighten_hex(base, 0.45) for grp, base in base_color_map.items()})
    two_tone.update({grp: base for grp, base in base_color_map.items()})
    legend_order: List[str] = []
    for grp in ordered_groups:
        legend_order += [f"{grp} — Tumor", f"{grp} — Reference"]
    return two_tone, legend_order

def make_metric_bar_figure(plot_df: pd.DataFrame, metric: str,
                           color_map: Mapping[str, str],
                           legend_order: Sequence[str]) -> go.Figure:
    # Batch order from legend_order → extract unique Group sequence
    group_order: List[str] = []
    for ck in legend_order:
        grp = ck.split(" — ")[0] if " — " in ck else ck
        if grp not in group_order:
            group_order.append(grp)

    # X order = for each batch (Group) in order, list its SampleIDs alphabetically
    x_order: List[str] = []
    for grp in group_order:
        sids = (plot_df.loc[plot_df["Group"] == grp, "SampleID"]
                        .astype(str)
                        .dropna()
                        .unique()
                        .tolist())
        x_order.extend(sorted(sids, key=str.lower))

    axis_fmt = THRESHOLDS.get(metric, {}).get("format")
    hover_fmt = ":.2%" if (axis_fmt or "").endswith("%") else ":.4f"

    fig = px.bar(
        plot_df,  # no value-based sort
        x="SampleID", y=metric,
        color="ColorKey",
        category_orders={
            "SampleID": x_order,
            "ColorKey": list(legend_order),
        },
        color_discrete_map=dict(color_map),
        hover_data={"SampleID": True, metric: hover_fmt, "Date": True, "Cycle": True,
                    "Group": True, "Role": True, "ColorKey": False},
        title=f"{metric} per Sample",
        labels={"SampleID": "Sample", metric: metric, "ColorKey": "Date/Cycle/Platform — Role"},
    )

    if metric in THRESHOLDS:
        for line in THRESHOLDS[metric].get("lines", []):
            fig.add_hline(
                y=line["value"], line_width=2, line_dash="dash", line_color=line["color"],
                annotation_text=line.get("label"), annotation_position="top left",
            )
        if axis_fmt:
            fig.update_yaxes(tickformat=axis_fmt)

    fig.update_layout(margin=dict(l=40, r=40, t=50, b=80), bargap=0.25, legend_title_text="Date / Cycle")
    fig.update_xaxes(tickangle=-45)
    return fig

def make_duplex_oneplus_figure(plot_df: pd.DataFrame) -> go.Figure:
    label_map = {"percent_duplex": "% Duplex", "percent_oneplus": "% OnePlus"}
    plot_df["MetricLabel"] = plot_df["Metric"].map(label_map).fillna(plot_df["Metric"])

    # Derive batch order from the same Group logic as metric plots
    # Approximate legend group order by Date/Cycle chronology already encoded in 'Group'
    # Keep the per-batch grouping and alpha order within batch
    # Build group order by the appearance order in the data grouped by Date/Cycle
    # (If you want to enforce the exact same order as metric tabs, pass it in similarly.)
    groups_in_order: List[str] = []
    for grp in plot_df["Group"].astype(str):
        if grp not in groups_in_order:
            groups_in_order.append(grp)

    x_order: List[str] = []
    for grp in groups_in_order:
        sids = (plot_df.loc[plot_df["Group"] == grp, "SampleID"]
                        .astype(str)
                        .dropna()
                        .unique()
                        .tolist())
        x_order.extend(sorted(sids, key=str.lower))

    fig = px.bar(
        plot_df,  # no sort_values
        x="SampleID", y="Value", color="MetricLabel",
        barmode="group",
        category_orders={"SampleID": x_order},
        hover_data={"SampleID": True, "MetricLabel": True, "Value": ":.1%", "Date": True, "Cycle": True,
                    "Group": True, "Role": True},
        title="Duplex vs OnePlus (% of Total Reads)",
        labels={"SampleID": "Sample", "Value": "Percentage", "MetricLabel": "Type"},
    )
    fig.update_yaxes(tickformat=".0%")
    fig.update_layout(margin=dict(l=40, r=40, t=50, b=80), bargap=0.25, legend_title_text="Type")
    fig.update_xaxes(tickangle=-45)
    return fig

# ========= UI =========
START_SCREEN = ui.nav_panel(
    "screen_start",
    ui.card(
        ui.h3("WGS Metrics – Bucket & Mode"),
        ui.input_text("bucket_input", "GCS bucket", value=DEFAULT_BUCKET, width="100%"),
        ui.input_radio_buttons(
            "mode", "Mode",
            choices={"sbx": "SBX Batches mode (auto from execution-definition.yaml)",
                     "manual": "Batch creation mode (manual)"},
            selected="sbx", inline=False
        ),
        ui.input_action_button("scan_bucket", "Scan bucket", class_="btn-primary"),
        ui.output_text("scan_summary"),
    )
)

MANUAL_SCREEN = ui.nav_panel(
    "screen_manual",
    ui.card(
        ui.h3("Manual batch creation"),
        ui.layout_columns(
            ui.card(
                ui.h4("Selectable samples"),
                ui.input_checkbox_group("manual_samples", None, choices={}, selected=[], width="100%", inline=False),
                ui.p("Selected samples disappear after you add them to a batch."),
            ),
            ui.card(
                ui.h4("Batch details"),
                ui.input_date("manual_date", "Date", value=None),
                ui.input_select("manual_platform", "Technology", choices=["SBX", "Illumina"], selected="SBX"),
                ui.input_text("manual_comment", "Extra comment (e.g., cycle01, downsampled-illumina)", value=""),
                ui.input_action_button("manual_create_batch", "Create batch from selected", class_="btn-primary"),
                ui.hr(),
                ui.h5("Created batches"),
                ui.output_ui("manual_batches_table"),
                ui.hr(),
                ui.input_action_button("manual_go_to_select", "Go to batch selection", class_="btn-success"),
            ),
            col_widths=(7,5)
        )
    )
)

SAMPLE_SCREEN = ui.nav_panel(
    "screen_samples",
    ui.div(
        ui.h2("Select batches to visualize"),
        ui.tags.style("""
        .batch-tiles .form-check{margin:0;padding:0}
        .batch-tiles .form-check-label{display:none}
        .batch-tiles .batch-item{cursor:pointer;}
        .batch-tiles .form-check-input{
            position:absolute;width:0;height:0;opacity:0;pointer-events:none;margin:0;
        }
        .batch-tiles .batch-item{
          display:grid; grid-template-columns: 40px 1fr 40px;
          align-items:center; justify-items:center;
          gap:10px; margin:6px auto; padding:10px 14px;
          min-height:44px; min-width:380px; max-width:640px; width:80%;
          background:rgba(255,255,255,.60); border:1px solid rgba(17,24,39,.10);
          border-radius:12px; backdrop-filter:blur(6px);
          transition:transform .12s, box-shadow .12s, background .12s, border-color .12s;
        }
        .batch-tiles .batch-item:hover{transform:scale(1.015);box-shadow:0 8px 18px rgba(0,0,0,.10)}
        .batch-tiles .batch-text{font-weight:600; letter-spacing:.2px; color:#111827; text-align:center; line-height:1.25;}
        .batch-tiles .batch-item:has(input[type="checkbox"]:checked){
           background:#f3efff; border-color:#7b68ee; box-shadow:0 8px 20px rgba(123,104,238,.25)
        }"""),
        ui.div(ui.output_ui("batch_feed_items"), class_="batch-tiles"),
        ui.layout_columns(
            ui.input_action_button("batch_select_all", "Select all", class_="btn-sm"),
            ui.input_action_button("batch_clear", "Clear", class_="btn-sm"),
            ui.input_action_button("load_data", "Load Data & View Dashboard", class_="btn-primary"),
            col_widths=(2,2,8)
        ),
        ui.div(
            ui.input_radio_buttons(
                "coverage_type", "Coverage data",
                choices={"cov": "General (MAPQ ≥ 20)",
                         "cov-hc":"High confidence regions statistics",
                         "cov-mq0":"All coverage statistics (MAPQ ≥ 0)"},
                selected="cov", inline=True,
            )
        ),
    ),
)

DASHBOARD_SCREEN = ui.nav_panel(
    "screen_dashboard",
    ui.page_sidebar(
        ui.sidebar(
            ui.card(
                ui.h4("Batch filters"),
                ui.input_date_range("batch_date_range", "Date range (batch)", start=None, end=None),
                ui.input_selectize("batch_cycles", "Runs (Cycle)", choices=[], multiple=True,
                                   options={"placeholder": "All runs"}),
                ui.input_radio_buttons(
                    "batch_role", "Sample type",
                    choices={"Any": "Any", "Tumor": "Tumor", "Normal": "Normal"},
                    selected="Any", inline=True,
                ),
                ui.layout_columns(
                    ui.input_action_button("batch_apply", "Apply",  class_="btn-primary",  style="width:100%"),
                    ui.input_action_button("batch_reset", "Reset",  class_="btn-secondary", style="width:100%"),
                    col_widths=(6, 6),
                ),
            ),
            ui.card(
                ui.h4("Displayed Samples"),
                ui.input_checkbox_group("filter_samples", None, choices=[]),
                ui.layout_columns(
                    ui.input_action_button("filter_select_all",   "Select All",   class_="btn-sm", style="width:100%"),
                    ui.input_action_button("filter_deselect_all", "Deselect All", class_="btn-sm", style="width:100%"),
                    col_widths=(6,6),
                ),
            ),
        ),
        ui.output_ui("dashboard_content"),
    ),
)

APP_UI = ui.page_fluid(
    ui.navset_hidden(
        START_SCREEN,
        MANUAL_SCREEN,
        SAMPLE_SCREEN,
        DASHBOARD_SCREEN,
        id="wizard",
    ),
)

# ========= Server =========
def server(input, output, session):
    # State
    gcs = reactive.Value[Optional[GCSFileSystem]](None)
    current_bucket = reactive.Value[str](DEFAULT_BUCKET)
    mode = reactive.Value[str]("sbx")

    # From scanning
    scanned_samples = reactive.Value[List[str]]([])   # for manual pool
    sample_meta = reactive.Value[Optional[pd.DataFrame]](None)  # normalized meta for tile screen (both modes)

    # Manual batches (list of dicts)
    manual_batches = reactive.Value[List[Dict[str, Any]]]([])

    # Tile bookkeeping
    batch_id_map = reactive.Value[Dict[str, str]]({})
    preselect_batch_keys = reactive.Value[set](set())

    # Loaded data
    full_data = reactive.Value[Optional[pd.DataFrame]](None)
    readlen_data = reactive.Value[Dict[str, pd.DataFrame]]({})
    role_color_map = reactive.Value[Dict[str, str]]({})
    color_key_order = reactive.Value[List[str]]([])

    # Connect GCS
    @reactive.Effect
    def _connect_gcs():
        if gcs.get() is not None:
            return
        try:
            fs = GCSFileSystem()
            gcs.set(fs)
        except Exception as e:
            ui.notification_show(f"GCS connection failed: {e}", duration=8, type="error")

    # Scan bucket
    @reactive.Effect
    @reactive.event(input.scan_bucket)
    def _scan_bucket():
        bkt = input.bucket_input().strip()
        if not bkt:
            ui.notification_show("Enter a bucket.", duration=5, type="warning"); return
        if not bkt.startswith("gs://"):
            ui.notification_show("Bucket must start with gs://", duration=5, type="warning"); return

        mode.set(input.mode())
        current_bucket.set(bkt)
        fs = gcs.get()
        if fs is None or not fs.exists(bkt):
            ui.notification_show(f"Bucket not found: {bkt}", duration=7, type="error"); return

        samples = find_sample_folders(fs, bkt)
        scanned_samples.set(sorted(samples))
        dbg("scan_bucket_done", bucket=bkt, mode=mode.get(), n_samples=len(samples), samples_head=_head(samples, 15))

        # SBX mode builds meta now; Manual goes to builder
        if mode.get() == "sbx":
            rows: List[Dict[str, Any]] = []
            with ui.Progress(min=0, max=max(1, len(samples))) as p:
                p.set(0, message="Reading execution-definition.yaml")
                for i, sid in enumerate(samples, start=1):
                    ypath = f"{bkt.rstrip('/')}/{sid}/execution-definition.yaml"
                    text = read_text_if_exists(fs, ypath)
                    uri, key_used = get_any_seq_uri(text)
                    date_iso, cycle = extract_date_cycle_from_bam_uri(uri)
                    dbg("meta_probe_sbx", bucket=bkt, sample_id=sid, yaml_found=bool(text),
                        key_used=key_used, uri=uri, date=date_iso, cycle=cycle)
                    rows.append({
                        "SampleID": sid,
                        "Platform": platform_from_sid(sid),
                        "Date": date_iso or "",
                        "Cycle": cycle or "",
                        "bam_uri": uri or "",
                        "yaml_found": bool(text),
                    })
                    p.set(i, detail=sid)
            meta_df = normalize_meta(pd.DataFrame(rows))
            sample_meta.set(meta_df)
            preselect_batch_keys.set(set())  # no auto preselect in SBX mode
            ui.update_navs("wizard", selected="screen_samples")
        else:
            # Manual mode → go to builder
            # populate choices
            ui.update_checkbox_group("manual_samples", choices={s: s for s in samples}, selected=[])
            ui.update_navs("wizard", selected="screen_manual")

    @render.text
    def scan_summary():
        return f"Bucket: {current_bucket.get()} | Mode: {mode.get()} | Found samples: {len(scanned_samples.get() or [])}"

    # Manual: create batch
    @reactive.Effect
    @reactive.event(input.manual_create_batch)
    def _manual_add_batch():
        pool = scanned_samples.get() or []
        selected = input.manual_samples() or []
        if not selected:
            ui.notification_show("Select at least one sample.", duration=5, type="warning"); return
        date_val = input.manual_date()
        date_iso = (date_val.isoformat() if date_val else "")
        platform = input.manual_platform() or "SBX"
        comment  = input.manual_comment().strip()

        batches = manual_batches.get()
        batches.append({
            "platform": platform,
            "date_iso": date_iso,
            "cycle": comment,           # we store comment in Cycle
            "samples": list(selected),
        })
        manual_batches.set(batches)

        # Remove from selectable pool
        remaining = [s for s in pool if s not in set(selected)]
        scanned_samples.set(remaining)
        ui.update_checkbox_group("manual_samples", choices={s: s for s in remaining}, selected=[])

        dbg("manual_batch_created", platform=platform, date_iso=date_iso, cycle=comment,
            n_added=len(selected), n_remaining=len(remaining))

    # Manual: render batches list
    @render.ui
    def manual_batches_table():
        batches = manual_batches.get() or []
        if not batches:
            return ui.p("No batches created yet.")
        rows = []
        for b in batches:
            rows.append(ui.tr(
                ui.td(b.get("date_iso") or ""),
                ui.td(b.get("platform") or ""),
                ui.td(b.get("cycle") or ""),
                ui.td(str(len(b.get("samples") or []))),
            ))
        return ui.table(
            {"class": "table table-sm"},
            ui.thead(ui.tr(ui.th("Date"), ui.th("Technology"), ui.th("Comment"), ui.th("# Samples"))),
            ui.tbody(*rows),
        )

    # Manual: proceed to selection
    @reactive.Effect
    @reactive.event(input.manual_go_to_select)
    def _manual_go_to_select_click():
        batches = manual_batches.get() or []
        dbg("manual_go_to_select_click",
            n_batches=len(batches),
            batches=[{
                "platform": b.get("platform"),
                "date_iso": b.get("date_iso"),
                "cycle": b.get("cycle"),
                "n_samples": len(b.get("samples") or []),
                "samples_head": (b.get("samples") or [])[:10],
            } for b in batches])

        if not batches:
            ui.notification_show("Create at least one batch first.", duration=6, type="warning"); return

        manual_df_raw = manual_batches_to_meta(batches)
        manual_df     = normalize_meta(manual_df_raw)
        sample_meta.set(manual_df)
        keys = {bk for bk in manual_df["BatchKey"].unique().tolist()}
        preselect_batch_keys.set(keys)
        dbg("manual_go_to_select_built", meta_rows=int(manual_df.shape[0]),
            meta_cols=list(manual_df.columns), preselect_keys=list(keys)[:20])
        ui.update_navs("wizard", selected="screen_samples")

    # Batch tiles
    @reactive.Calc
    def _batch_choices() -> Dict[str, str]:
        df = sample_meta.get()
        return build_batch_choices(df) if isinstance(df, pd.DataFrame) else {}

    @render.ui
    def batch_feed_items():
        choices = _batch_choices()
        preselect = preselect_batch_keys.get() or set()
        dbg("batch_feed_items_render", n_choices=len(choices), choice_keys=list(choices.keys()),
            preselect_keys=list(preselect))
        id_map: Dict[str, str] = {}
        items: List[Any] = []
        for key, lab in choices.items():
            sid = hashlib.md5(key.encode("utf-8")).hexdigest()[:8]
            safe_key = re.sub(r"\W+", "_", key)    
            safe_id = f"bb_{safe_key}_{sid}"        
            id_map[safe_id] = key
            default_checked = (key in preselect)
            dbg("batch_tile", key=key, input_id=safe_id, default_checked=default_checked, label=lab)
            items.append(
                ui.tags.label(
                    ui.input_checkbox(id=safe_id, label="", value=default_checked),
                    ui.tags.span(lab, class_="batch-text"),
                    ui.tags.span("", class_="tile-spacer"),
                    class_="batch-item",
                    **{"for": safe_id}
                )
            )
        batch_id_map.set(id_map)
        dbg("batch_feed_items_render_done", n_inputs=len(id_map), input_ids=list(id_map.keys()))
        return ui.tags.div(*items)

    @reactive.Effect 
    @reactive.event(input.batch_select_all)
    def _select_all_batches():
        for sid in batch_id_map.get().keys():
            ui.update_checkbox(sid, value=True)

    @reactive.Effect 
    @reactive.event(input.batch_clear)
    def _clear_batches():
        for sid in batch_id_map.get().keys():
            ui.update_checkbox(sid, value=False)

    # Selected samples from tiles
    def _selected_samples_from_batches() -> List[str]:
        df = sample_meta.get()
        if df is None or df.empty:
            return []
        idmap = batch_id_map.get()
        tile_states = []
        for safe_id, orig_key in idmap.items():
            try:
                tile_states.append({"input_id": safe_id, "key": orig_key, "val": bool(getattr(input, safe_id)())})
            except Exception:
                tile_states.append({"input_id": safe_id, "key": orig_key, "val": False})
        dbg("selected_samples_scan_start", n_inputs=len(tile_states))
        trues = [t for t in tile_states if t["val"]]
        dbg("selected_samples_tiles", n_true=len(trues), examples_true=_head(trues, 5), examples_false=_head([t for t in tile_states if not t["val"]], 5))
        idx = build_batch_index(df)
        dbg("selected_samples_index_built", n_batches=len(idx), batch_keys=_head(idx.keys(), 10))
        out: List[str] = []
        for t in trues:
            key = t["key"]
            sids = idx.get(key, [])
            dbg("selected_samples_match_key_indexed", key=key, n_matched=len(sids), sids_head=_head(sids, 10))
            out.extend(sids)
        out = sorted(set(out))
        dbg("selected_samples_from_batches_done", n=len(out), sample_ids_head=_head(out, 15))
        return out

    # Populate cycles in sidebar
    @reactive.Effect
    def _populate_cycles():
        df_loaded = full_data.get()
        if isinstance(df_loaded, pd.DataFrame) and not df_loaded.empty and "Cycle" in df_loaded.columns:
            cycles = sorted({c.strip() for c in df_loaded["Cycle"].astype(str) if c and c.strip() and c.lower() != "nan"})
            dbg("populate_cycles", source="full_data", n=len(cycles), cycles=_head(cycles, 30))
            current = input.batch_cycles() or []
            new_sel = [c for c in current if c in cycles] or cycles
            ui.update_selectize("batch_cycles", choices=cycles, selected=new_sel)
            return
        df_meta = sample_meta.get()
        if not isinstance(df_meta, pd.DataFrame) or df_meta.empty or "Cycle" not in df_meta.columns:
            dbg("populate_cycles", source="meta", status="empty_meta")
            ui.update_selectize("batch_cycles", choices=[], selected=[])
            return
        cyc = df_meta["Cycle"].astype(str)
        cyc = cyc[cyc.str.strip().ne("") & cyc.str.lower().ne("nan")]
        cycles = sorted(set(cyc.tolist()))
        dbg("populate_cycles", source="meta", n=len(cycles), cycles=_head(cycles, 30))
        current = input.batch_cycles() or []
        new_sel = [c for c in current if c in cycles] or cycles
        ui.update_selectize("batch_cycles", choices=cycles, selected=new_sel)

    # Enter sample tile screen debug
    @reactive.Effect
    def _on_enter_screen_samples():
        try:
            if input.wizard() == "screen_samples":
                dfm = sample_meta.get()
                choices = build_batch_choices(dfm) if isinstance(dfm, pd.DataFrame) else {}
                idx = build_batch_index(dfm) if isinstance(dfm, pd.DataFrame) else {}
                dbg("enter_screen_samples",
                    n_meta_rows=(0 if not isinstance(dfm, pd.DataFrame) else int(dfm.shape[0])),
                    sample_ids_head=(_head(dfm["SampleID"].astype(str).tolist(), 20) if isinstance(dfm, pd.DataFrame) else []),
                    batch_choices=_head(list(choices.keys()), 20),
                    index_keys=_head(list(idx.keys()), 20),
                    preselect_keys=_head(list(preselect_batch_keys.get() or set()), 20),
                    n_input_ids=len(batch_id_map.get() or {}),
                )
        except Exception as e:
            dbg("enter_screen_samples_error", error=repr(e))

    # Load data
    @reactive.Effect
    @reactive.event(input.load_data)
    def _load_data():
        bucket = current_bucket.get()
        dbg("load_click", bucket=bucket, mode=mode.get(),
            n_choices=len(_batch_choices()), choice_keys=list(_batch_choices().keys()),
            n_inputs=len(batch_id_map.get()), input_ids=list(batch_id_map.get().keys()),
            preselect_keys=list(preselect_batch_keys.get()))
        selected_samples = _selected_samples_from_batches()
        if not selected_samples:
            dbg("load_click_no_selection",
                tile_states=[{"input_id": sid, "key": key, "val": bool(getattr(input, sid)())} for sid, key in batch_id_map.get().items()],
                note="UI says none checked; see tile_states to confirm.")
            ui.notification_show("Pick at least one batch.", duration=5, type="warning")
            return

        fs = gcs.get()
        if fs is None:
            ui.notification_show("No GCS connection.", duration=5, type="error"); return
        meta_local = sample_meta.get()
        meta_local = meta_local if isinstance(meta_local, pd.DataFrame) else pd.DataFrame(columns=["SampleID","Date","Cycle","bam_uri","yaml_found","Platform"])

        series_rows: List[pd.Series] = []
        histograms: Dict[str, pd.DataFrame] = {}

        with ui.Progress(min=0, max=max(1, len(selected_samples))) as p:
            p.set(message="Loading metrics files...", detail="Starting...")
            for i, sample_id in enumerate(selected_samples, start=1):
                p.set(i-1, detail=f"Processing {sample_id}")
                cov_folder = input.coverage_type()
                cov_files = fs.glob(f"{bucket.rstrip('/')}/{sample_id}/{cov_folder}/*.wgs.metrics.txt")
                gc_files  = fs.glob(f"{bucket.rstrip('/')}/{sample_id}/mmetrics/*gc_bias.summary_metrics*")
                ins_files = fs.glob(f"{bucket.rstrip('/')}/{sample_id}/mmetrics/*insert_size_metrics*")
                aln_files = fs.glob(f"{bucket.rstrip('/')}/{sample_id}/mmetrics/*alignment_summary_metrics*")

                cov_row = parse_picard_single_metrics(fs, cov_files[0]) if cov_files else None
                gc_row  = parse_picard_single_metrics(fs, gc_files[0])  if gc_files  else None
                ins_row = parse_picard_single_metrics(fs, ins_files[0]) if ins_files else None

                aln_table = None
                if aln_files:
                    hist_df = parse_alignment_summary_histogram(fs, aln_files[0])
                    if hist_df is not None and not hist_df.empty:
                        histograms[sample_id] = hist_df
                    aln_table = parse_alignment_summary_table(fs, aln_files[0])

                merged: Dict[str, Any] = {}
                for s in (cov_row, gc_row, ins_row):
                    if s is not None:
                        merged.update({str(k).strip(): v for k, v in s.to_dict().items()})
                if aln_table is not None:
                    merged.update({str(k).strip(): v for k, v in aln_table.to_dict().items()})

                # runlevel TSV attach if available
                row_meta = meta_local.loc[meta_local["SampleID"] == sample_id] if not meta_local.empty else pd.DataFrame()
                bam_uri = str(row_meta["bam_uri"].iat[0]).strip() if (not row_meta.empty and "bam_uri" in row_meta) else ""
                if bam_uri:
                    tsv = runlevel_tsv_from_bam_uri(fs, bam_uri)
                    if tsv:
                        rv = parse_runlevel_for_sample(fs, tsv, sample_id)
                        if rv:
                            merged.update(rv)

                # ensure expected fields
                needed = [
                    "Phred_All","Phred_Substitutions","Phred_Insertions","Phred_Deletions","Phred_Indel",
                    "median_read_length","mean_read_length","ratio_90_to_10percentile",
                    "Num_Full_Length_Reads","Num_One_Plus_Reads","Total_Reads",
                ]
                for k in needed:
                    merged.setdefault(k, None)

                # derived duplex/oneplus %
                try:
                    fl = float(merged.get("Num_Full_Length_Reads")) if merged.get("Num_Full_Length_Reads") is not None else float("nan")
                    op = float(merged.get("Num_One_Plus_Reads"))  if merged.get("Num_One_Plus_Reads")  is not None else float("nan")
                    tot = float(merged.get("Total_Reads"))        if merged.get("Total_Reads")        is not None else float("nan")
                    merged["percent_duplex"]  = (fl / tot) if (pd.notna(fl) and pd.notna(tot) and tot > 0) else None
                    merged["percent_oneplus"] = (op / tot) if (pd.notna(op) and pd.notna(tot) and tot > 0) else None
                except Exception:
                    merged["percent_duplex"]  = None
                    merged["percent_oneplus"] = None

                # histogram-derived stats
                hist_df = histograms.get(sample_id)
                def read_length_n50_from_hist(hd: Optional[pd.DataFrame]) -> Optional[float]:
                    if hd is None or hd.empty or "READ_LENGTH" not in hd.columns:
                        return None
                    count_col = "TOTAL" if "TOTAL" in hd.columns else ("ALIGNED" if "ALIGNED" in hd.columns else None)
                    if not count_col:
                        return None
                    df = hd[["READ_LENGTH", count_col]].dropna().copy()
                    df["READ_LENGTH"] = pd.to_numeric(df["READ_LENGTH"], errors="coerce")
                    df[count_col]     = pd.to_numeric(df[count_col], errors="coerce")
                    df = df.dropna()
                    if df.empty:
                        return None
                    df = df.sort_values("READ_LENGTH", ascending=False)
                    df["bases"] = df["READ_LENGTH"] * df[count_col]
                    total_bases = df["bases"].sum()
                    if total_bases <= 0:
                        return None
                    df["cum_bases"] = df["bases"].cumsum()
                    idx = df.index[df["cum_bases"] >= (0.5 * total_bases)]
                    return float(df.loc[idx[0], "READ_LENGTH"]) if len(idx) else None
                def median_read_length_from_hist(hd: Optional[pd.DataFrame]) -> Optional[float]:
                    if hd is None or hd.empty or "READ_LENGTH" not in hd.columns:
                        return None
                    count_col = "TOTAL" if "TOTAL" in hd.columns else ("ALIGNED" if "ALIGNED" in hd.columns else None)
                    if not count_col:
                        return None
                    df = hd[["READ_LENGTH", count_col]].dropna().copy()
                    df["READ_LENGTH"] = pd.to_numeric(df["READ_LENGTH"], errors="coerce")
                    df[count_col]     = pd.to_numeric(df[count_col], errors="coerce")
                    df = df[(df["READ_LENGTH"] > 0) & df[count_col].notna()]
                    if df.empty:
                        return None
                    df = df.sort_values("READ_LENGTH", ascending=True)
                    total_reads = df[count_col].sum()
                    if total_reads <= 0:
                        return None
                    df["cum_reads"] = df[count_col].cumsum()
                    idx = df.index[df["cum_reads"] >= (0.5 * total_reads)]
                    return float(df.loc[idx[0], "READ_LENGTH"]) if len(idx) else None

                merged["READ_LENGTH_N50"] = read_length_n50_from_hist(hist_df)

                # Prefer Alignment Summary's MEDIAN_READ_LENGTH; fallback to histogram
                if "MEDIAN_READ_LENGTH" in merged and pd.notna(merged["MEDIAN_READ_LENGTH"]):
                    merged["mmetrics_median_read_length"] = merged["MEDIAN_READ_LENGTH"]
                else:
                    merged["mmetrics_median_read_length"] = median_read_length_from_hist(hist_df)

                # flagstat: duplicates + total reads (passed + failed)
                try:
                    flag_path = find_flagstat_txt(fs, bucket, sample_id)
                    if flag_path:
                        txt = read_text_if_exists(fs, flag_path) or ""
                        # duplicates
                        m_dup = re.search(r"^\s*(\d+)\s+\+\s+\d+\s+duplicates\s*$", txt, flags=re.M | re.I)
                        if m_dup:
                            merged["FLAGSTAT_DUPLICATES"] = float(m_dup.group(1))

                        # total reads: "<passed> + <failed> in total (QC-passed reads + QC-failed reads)"
                        m_tot = re.search(r"^\s*(\d+)\s+\+\s+(\d+)\s+in\s+total\b", txt, flags=re.M | re.I)
                        if m_tot:
                            passed = int(m_tot.group(1))
                            failed = int(m_tot.group(2))
                            merged["FLAGSTAT_TOTAL_READS"] = float(passed + failed)
                            dbg("flagstat_total_reads", sample_id=sample_id, passed=passed, failed=failed, total=int(passed+failed))
                        else:
                            # Fallback: grab the first integer on the "in total" line if format varies
                            m_line = re.search(r"^(.*in\s+total.*)$", txt, flags=re.M | re.I)
                            if m_line:
                                nums = re.findall(r"\d+", m_line.group(1))
                                if len(nums) >= 2:
                                    merged["FLAGSTAT_TOTAL_READS"] = float(int(nums[0]) + int(nums[1]))
                                    dbg("flagstat_total_reads_fallback_sum", sample_id=sample_id, nums=nums[:4])
                                elif len(nums) >= 1:
                                    merged["FLAGSTAT_TOTAL_READS"] = float(int(nums[0]))
                                    dbg("flagstat_total_reads_fallback_single", sample_id=sample_id, nums=nums[:4])
                except Exception as e:
                    dbg("flagstat_error", sample_id=sample_id, error=repr(e))

                if merged:
                    merged["SampleID"] = sample_id
                    series_rows.append(pd.Series(merged))

        if not series_rows:
            ui.notification_show("Failed to load data.", duration=5, type="error"); return

        df = pd.DataFrame(series_rows)
        df.columns = df.columns.astype(str).str.strip()

        # NEW (fixed):
        def _get_numeric_series(dataframe: pd.DataFrame, col_name: str) -> pd.Series:
            """Safely get a numeric Series from a DataFrame column, returning NA series if missing."""
            if col_name in dataframe.columns:
                return pd.to_numeric(dataframe[col_name], errors="coerce")
            return pd.Series([pd.NA] * len(dataframe), index=dataframe.index)

        flag_tot = _get_numeric_series(df, "FLAGSTAT_TOTAL_READS")
        tsv_tot = _get_numeric_series(df, "Total_Reads") if "Total_Reads" in df.columns else _get_numeric_series(df, "TOTAL_READS")
        total_reads = flag_tot.where(flag_tot.notna(), tsv_tot)

        mean_len_as = _get_numeric_series(df, "MEAN_READ_LENGTH")
        mean_len_fb = _get_numeric_series(df, "mean_read_length")
        mean_len = mean_len_as.where(mean_len_as.notna(), mean_len_fb)

        df["Yield"] = (total_reads * mean_len) / 1e9

        # ===== Debug logging for Yield =====
        try:
            dbg(
                "yield_sources_overview",
                n_rows=int(df.shape[0]),
                n_flagstat_total= int(flag_tot.notna().sum()) if flag_tot is not None else 0,
                n_tsv_total=      int(tsv_tot.notna().sum())  if tsv_tot  is not None else 0,
                n_total_reads=    int(total_reads.notna().sum()),
                n_mean_as=        int(mean_len_as.notna().sum()) if mean_len_as is not None else 0,
                n_mean_fb=        int(mean_len_fb.notna().sum()) if mean_len_fb is not None else 0,
                n_mean_final=     int(mean_len.notna().sum()),
                n_yield=          int(df["Yield"].notna().sum()),
                n_yield_missing=  int(df["Yield"].isna().sum())
            )

            # Per-sample trace (first 50)
            cols = [c for c in ["SampleID","FLAGSTAT_TOTAL_READS","Total_Reads","TOTAL_READS","MEAN_READ_LENGTH","mean_read_length","Yield"] if c in df.columns]
            peek = df[cols].head(50).copy()
            for _, row in peek.iterrows():
                sid = str(row.get("SampleID"))
                fr  = pd.to_numeric(row.get("FLAGSTAT_TOTAL_READS"), errors="coerce")
                trA = pd.to_numeric(row.get("Total_Reads"), errors="coerce")
                trB = pd.to_numeric(row.get("TOTAL_READS"), errors="coerce")
                tr  = fr if pd.notna(fr) else (trA if pd.notna(trA) else trB)
                mlA = pd.to_numeric(row.get("MEAN_READ_LENGTH"), errors="coerce")
                mlB = pd.to_numeric(row.get("mean_read_length"), errors="coerce")
                ml  = mlA if pd.notna(mlA) else mlB
                src_tr = "FLAGSTAT_TOTAL_READS" if pd.notna(fr) else ("Total_Reads" if pd.notna(trA) else ("TOTAL_READS" if pd.notna(trB) else None))
                src_ml = "MEAN_READ_LENGTH" if pd.notna(mlA) else ("mean_read_length" if pd.notna(mlB) else None)
                yv  = pd.to_numeric(row.get("Yield"), errors="coerce")
                dbg("yield_calc_sample", sample_id=sid, total_reads_source=src_tr, total_reads=(float(tr) if pd.notna(tr) else None),
                    mean_length_source=src_ml, mean_length=(float(ml) if pd.notna(ml) else None),
                    yield_gb=(float(yv) if pd.notna(yv) else None))
        except Exception as e:
            dbg("yield_logging_error", error=repr(e))

        # numeric coercion
        for col in METRICS_TO_PLOT:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        extras = [
            "Phred_All","Phred_Substitutions","Phred_Insertions","Phred_Deletions","Phred_Indel",
            "median_read_length","mean_read_length","ratio_90_to_10percentile",
            "Num_Full_Length_Reads","Num_One_Plus_Reads","Total_Reads",
            "percent_duplex","percent_oneplus","READ_LENGTH_N50","TOTAL_READS",
            "mmetrics_median_read_length","PF_ALIGNED_BASES","FLAGSTAT_DUPLICATES",
            "MEAN_READ_LENGTH","MEDIAN_READ_LENGTH", "FLAGSTAT_TOTAL_READS",
        ]
        for c in extras:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        for c in PCT_LIKE_COLUMNS:
            if c in df.columns:
                s = pd.to_numeric(df[c], errors="coerce")
                if s.dropna().any() and 1.0 < float(s.max()) <= 100.0:
                    df[c] = s / 100.0

        # attach meta for display
        meta_df = sample_meta.get()
        meta_df = meta_df if isinstance(meta_df, pd.DataFrame) else pd.DataFrame(columns=["SampleID","Date","Cycle","Platform"])
        meta = meta_df[["SampleID","Date","Cycle","Platform","DateDisplay","Run"]].copy()
        df = df.merge(meta, on="SampleID", how="left")

        df["Group"] = df.apply(lambda r: "Illumina" if r.get("Platform") == "Illumina"
                               else f"{r.get('DateDisplay','unknown')} • {r.get('Run','Run unknown')}", axis=1)
        # Role inference
        def classify_role(sample_id: str) -> str:
            sid = re.sub(r"-(sbx|nseq6)$", "", str(sample_id).strip(), flags=re.I)
            if re.search(r"(?i)(?:^|[-_])ref(?:$|[-_])", sid): return "Reference"
            if re.match(r"(?i)^h\d{5,}", sid): return "Tumor"
            if re.search(r"(?i)T(?:\d+|I+)?$", sid): return "Tumor"
            if re.search(r"(?i)R(?:\d+|I+)?$", sid): return "Reference"
            return "Unknown"
        df["Role"] = df["SampleID"].astype(str).map(classify_role)
        df["ColorKey"] = df.apply(lambda r: f"{r['Group']} — {r['Role']}"
                                  if r["Role"] in ("Tumor","Reference") else r["Group"], axis=1)

        for c in ("Group","Date","Cycle"):
            if c in df.columns:
                df[c] = df[c].map(_clean_str)

        cmap, order = build_two_tone_color_map(df)
        role_color_map.set(cmap)
        existing = set(df["ColorKey"].dropna().astype(str).unique())
        color_key_order.set([k for k in order if k in existing])

        full_data.set(df); readlen_data.set(histograms)

        # init sidebar filters
        try:
            dd = pd.to_datetime(df["Date"], errors="coerce")
            if dd.notna().any():
                dmin = dd.min().date(); dmax = dd.max().date()
                dbg("init_date_range", min=str(dmin), max=str(dmax))
                ui.update_date_range("batch_date_range", start=dmin, end=dmax)
            cycles = sorted({c.strip() for c in df["Cycle"].astype(str) if c and c.strip() and c.lower() != "nan"})
            dbg("init_cycles", cycles=cycles)
            ui.update_selectize("batch_cycles", choices=cycles, selected=cycles)
        except Exception as e:
            dbg("init_filters_error", error=repr(e))

        # initialize "Displayed Samples"
        sids = df["SampleID"].dropna().astype(str).tolist()
        ui.update_checkbox_group("filter_samples", choices={sid: sid for sid in sids}, selected=sids)
        ui.update_navs("wizard", selected="screen_dashboard")

        dbg("loaded_summary",
            n_rows=int(df.shape[0]),
            n_samples=int(df["SampleID"].nunique()),
            has_columns=list(df.columns)[:50],
            role_counts=df["Role"].astype(str).value_counts(dropna=False).to_dict(),
            date_min=str(pd.to_datetime(df["Date"], errors="coerce").min().date()) if pd.to_datetime(df["Date"], errors="coerce").notna().any() else None,
            date_max=str(pd.to_datetime(df["Date"], errors="coerce").max().date()) if pd.to_datetime(df["Date"], errors="coerce").notna().any() else None,
            n_unique_cycles=len(set(df["Cycle"].astype(str))),
            unique_cycles=_head(sorted(set(df["Cycle"].astype(str)))),
            preview=[{"SampleID": r["SampleID"], "Date": r["Date"], "Cycle": r["Cycle"], "Role": r["Role"]}
                     for _, r in df[["SampleID","Date","Cycle","Role"]].head(8).iterrows()]
        )

    # Sidebar: select all/none
    @reactive.Effect 
    @reactive.event(input.filter_select_all)
    def _filter_all():
        df = full_data.get()
        if df is None: return
        sids = df["SampleID"].astype(str).tolist()
        ui.update_checkbox_group("filter_samples", selected=sids)

    @reactive.Effect 
    @reactive.event(input.filter_deselect_all)
    def _filter_none():
        ui.update_checkbox_group("filter_samples", selected=[])

    # Batch filters (sidebar)
    def ids_matching_batch() -> List[str]:
        df = full_data.get()
        if df is None or df.empty:
            return []
        use = df.loc[:, ["SampleID","Date","Cycle","Role"]].copy()
        use["_Date"] = pd.to_datetime(use["Date"], errors="coerce")
        start, end = input.batch_date_range()
        if start is not None and end is not None:
            m = use["_Date"].notna() & (use["_Date"].dt.date >= start) & (use["_Date"].dt.date <= end)
            use = use.loc[m]
            if use.empty:
                return []
        sel_cycles = [c for c in (input.batch_cycles() or []) if c]
        if sel_cycles:
            use = use[use["Cycle"].astype(str).str.strip().isin(sel_cycles)]
        want_role = input.batch_role()
        if want_role and want_role != "Any":
            map_role = use["Role"].map({"Reference": "Normal", "Tumor": "Tumor"}).fillna("Unknown")
            use = use.loc[map_role == want_role]
        return sorted(use["SampleID"].astype(str).unique().tolist())

    @reactive.Effect 
    @reactive.event(input.batch_apply)
    def _batch_apply():
        ids = ids_matching_batch()
        if not ids:
            ui.notification_show("No samples match the chosen batch filters.", duration=5, type="warning"); return
        df = full_data.get()
        if df is not None and not df.empty:
            existing = set(df["SampleID"].astype(str))
            ids = [s for s in ids if s in existing]
        ui.update_checkbox_group("filter_samples", selected=ids)

    @reactive.Effect 
    @reactive.event(input.batch_reset)
    def _batch_reset():
        ui.update_selectize("batch_cycles", selected=[])
        ui.update_radio_buttons("batch_role", selected="Any")
        ui.update_date_range("batch_date_range", start=None, end=None)
        df = full_data.get()
        all_ids = df["SampleID"].astype(str).tolist() if (df is not None and not df.empty) else []
        if all_ids:
            ui.update_checkbox_group("filter_samples", selected=all_ids)

    # Dashboard content
    @render.ui
    def dashboard_content():
        if full_data.get() is None:
            return ui.p("Data is loading or no data was loaded.")
        tabs = [ui.nav_panel(m, output_widget(f"plot_{m}")) for m in METRICS_TO_PLOT]
        tabs.append(ui.nav_panel(
            "Table Overview",
            ui.div(ui.download_button("download_table_csv", "Download CSV"), style="margin-bottom:8px"),
            ui.output_data_frame("table_overview"),
        ))
        tabs.append(ui.nav_panel("Duplex vs OnePlus (%)", output_widget("plot_duplex_oneplus")))
        return ui.navset_card_tab(*tabs)

    # Download CSV
    @render.download(filename=lambda: DOWNLOAD_FILENAME_FMT.format(dt.datetime.now()))
    def download_table_csv():
        df = full_data.get()
        if df is None or df.empty:
            yield "No data loaded\n".encode("utf-8"); return
        selected = input.filter_samples() or []
        if not selected:
            yield "No samples selected\n".encode("utf-8"); return
        view = df[df["SampleID"].astype(str).isin(selected)].copy()
        role_series = view["Role"].fillna("Unknown").astype(str)
        view["SampleType"] = role_series.map({"Reference": "Normal", "Tumor": "Tumor"}).fillna("Unknown")
        cols_present = ["SampleID", "Date", "Cycle", "SampleType", "percent_duplex", "percent_oneplus"]
        for col in METRICS_TO_PLOT:
            if col in view.columns:
                view[col] = pd.to_numeric(view[col], errors="coerce")
                cols_present.append(col)
        view["_Date"] = pd.to_datetime(view["Date"], errors="coerce")
        view = view.sort_values(by=["_Date", "Cycle", "SampleID"], ascending=[False, True, True]).drop(columns=["_Date"], errors="ignore")
        view = view[[c for c in cols_present if c in view.columns]].reset_index(drop=True)
        numeric_cols = [c for c in view.columns if c not in ["SampleID","Date","Cycle","SampleType"]]
        for c in numeric_cols: view[c] = pd.to_numeric(view[c], errors="coerce")
        avg_vals = view[numeric_cols].mean(skipna=True)
        sum_vals = view[numeric_cols].sum(skipna=True)
        avg_row = {c: "" for c in view.columns}; avg_row["SampleID"]="AVERAGE (selected)"; avg_row["SampleType"]="—"
        for c in numeric_cols: avg_row[c] = avg_vals.get(c, pd.NA)
        sum_row = {c: "" for c in view.columns}; sum_row["SampleID"]="SUM (selected)"; sum_row["SampleType"]="—"
        for c in numeric_cols: sum_row[c] = sum_vals.get(c, pd.NA)
        out_df = pd.concat([pd.DataFrame([avg_row, sum_row])[view.columns], view], ignore_index=True)
        yield out_df.to_csv(index=False).encode("utf-8")

    # Table render
    @render.data_frame
    def table_overview():
        df = full_data.get()
        if df is None or df.empty:
            return render.DataTable(pd.DataFrame({"Info": ["No data loaded"]}), selection_mode="none")
        selected = input.filter_samples() or []
        if not selected:
            return render.DataTable(pd.DataFrame({"Info": ["No samples selected"]}), selection_mode="none")
        view = df[df["SampleID"].astype(str).isin(selected)].copy()
        role_series = view["Role"].fillna("Unknown").astype(str)
        view["SampleType"] = role_series.map({"Reference": "Normal", "Tumor": "Tumor"}).fillna("Unknown")
        cols_present = ["SampleID", "Date", "Cycle", "SampleType", "percent_duplex", "percent_oneplus"]
        for col in METRICS_TO_PLOT:
            if col in view.columns:
                view[col] = pd.to_numeric(view[col], errors="coerce")
                cols_present.append(col)
        view["_Date"] = pd.to_datetime(view["Date"], errors="coerce")
        view = view.sort_values(by=["_Date", "Cycle", "SampleID"], ascending=[False, True, True]).drop(columns=["_Date"], errors="ignore")
        numeric_cols: List[str] = []
        for c in view.columns:
            if c not in ["SampleID","Date","Cycle","SampleType"]:
                view[c] = pd.to_numeric(view[c], errors="coerce")
                numeric_cols.append(c)
        avg_vals = view[numeric_cols].mean(skipna=True)
        sum_vals = view[numeric_cols].sum(skipna=True)
        avg_row = {c: "" for c in view.columns}; avg_row["SampleID"]="AVERAGE (selected)"; avg_row["SampleType"]="—"
        for c in numeric_cols: avg_row[c] = avg_vals.get(c, pd.NA)
        sum_row = {c: "" for c in view.columns}; sum_row["SampleID"]="SUM (selected)"; sum_row["SampleType"]="—"
        for c in numeric_cols: sum_row[c] = sum_vals.get(c, pd.NA)
        final_cols = [c for c in cols_present if c in view.columns]
        view = view[final_cols].reset_index(drop=True)
        view = pd.concat([pd.DataFrame([avg_row, sum_row])[final_cols], view], ignore_index=True)
        return render.DataTable(view, selection_mode="none")

    # Register metric plots
    def _make_metric_output(metric_name: str):
        @output(id=f"plot_{metric_name}")
        @render_widget
        def _plot():
            df = full_data.get()
            selected = input.filter_samples() or []
            if df is None:
                return go.Figure(layout={"title": f"No data loaded for {metric_name}"})
            if not selected:
                return go.Figure(layout={"title": f"No samples selected for {metric_name}"})
            if metric_name not in df.columns:
                return go.Figure(layout={"title": f"{metric_name} not found in data"})
            plot_df = df[df["SampleID"].isin(selected)].copy()
            cmap = role_color_map.get() or {}
            order = color_key_order.get() or []
            return make_metric_bar_figure(plot_df, metric_name, cmap, order)
        return _plot
    for _metric in METRICS_TO_PLOT:
        _make_metric_output(_metric)

    # Duplex vs OnePlus plot
    @output(id="plot_duplex_oneplus")
    @render_widget
    def _plot_duplex_oneplus():
        df = full_data.get()
        selected = input.filter_samples() or []
        if df is None or not selected:
            return go.Figure(layout={"title": "No samples selected for Duplex vs OnePlus"})
        use = df[df["SampleID"].astype(str).isin(selected)].copy()
        if use.empty:
            return go.Figure(layout={"title": "No data"})
        cols = [c for c in ("percent_duplex", "percent_oneplus") if c in use.columns]
        if not cols:
            return go.Figure(layout={"title": "Required fields not found"})
        plot_df = use[["SampleID","Date","Cycle","Group","Role","ColorKey"] + cols].melt(
            id_vars=["SampleID","Date","Cycle","Group","Role","ColorKey"],
            value_vars=cols, var_name="Metric", value_name="Value"
        ).dropna(subset=["Value"])
        return make_duplex_oneplus_figure(plot_df)

# ========= Run =========
app = App(APP_UI, server)