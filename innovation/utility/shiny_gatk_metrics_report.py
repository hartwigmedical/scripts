#!/usr/bin/env python3
"""
Shiny dashboard for WGS metrics (Picard + runlevel TSV + flagstat).

Refactor highlights:
- Strong typing & docstrings
- Centralized constants/config
- Consolidated helpers (ID/labels, parsing, normalization)
- Robust matching for sample ↔ runlevel TSV
- Unified plot builders
- Reduced repetition in table/download logic
"""

from __future__ import annotations

# ========= Standard lib =========
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
    "mmetrics_median_read_length", "PF_ALIGNED_BASES", "FLAGSTAT_DUPLICATES",
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


# ========= Debug / Logging =========

def _log_filter_inputs(action: str, input_obj) -> None:   # <-- pass input explicitly
    try:
        start, end = input_obj.batch_date_range()         # <-- use input_obj
        sel_cycles = input_obj.batch_cycles() or []
        want_role = input_obj.batch_role()
        dbg("filter_inputs",
            action=action,
            want_role=want_role,
            start_type=type(start).__name__ if start is not None else None,
            end_type=type(end).__name__ if end is not None else None,
            start=str(start) if start is not None else None,
            end=str(end) if end is not None else None,
            n_sel_cycles=len(sel_cycles),
            sel_cycles=sel_cycles,
        )
    except Exception as e:
        dbg("filter_inputs_error", action=action, error=repr(e))

def _log_loaded_dataset_summary(df: pd.DataFrame) -> None:
    if df is None or df.empty:
        dbg("loaded_summary", status="empty_or_none"); return
    role_counts = (df["Role"].astype(str).value_counts(dropna=False).to_dict()
                   if "Role" in df.columns else {})
    # Parsed date range
    dd = pd.to_datetime(df["Date"], errors="coerce")
    date_min = str(pd.to_datetime(dd.min()).date()) if dd.notna().any() else None
    date_max = str(pd.to_datetime(dd.max()).date()) if dd.notna().any() else None
    # Unique cycles (stringified)
    uniq_cycles = sorted(set(df["Cycle"].astype(str))) if "Cycle" in df.columns else []
    # Preview the exact fields we filter on
    preview = df.loc[:, ["SampleID","Date","Cycle","Role"]].head(30).to_dict("records") \
        if all(c in df.columns for c in ["SampleID","Date","Cycle","Role"]) else []
    dbg("loaded_summary",
        n_rows=int(df.shape[0]),
        n_samples=int(df["SampleID"].nunique()) if "SampleID" in df.columns else None,
        has_columns=list(df.columns)[:50],
        role_counts=role_counts,
        date_min=date_min, date_max=date_max,
        n_unique_cycles=len(uniq_cycles), unique_cycles=uniq_cycles[:50],
        preview_count=len(preview), preview=preview
    )

def dbg(tag: str, **kv: Any) -> None:
    """Compact stderr logger with timestamp and truncated sequences."""
    ts = dt.datetime.now().strftime("%H:%M:%S")
    parts = []
    for k, v in kv.items():
        if isinstance(v, (list, tuple)) and len(v) > 10:
            parts.append(f"{k}=[{', '.join(map(repr, v[:10]))}, … +{len(v)-10}]")
        else:
            parts.append(f"{k}={repr(v)}")
    print(f"[{ts}] {tag} " + " ".join(parts), file=sys.stderr, flush=True)


# ========= Helpers: IDs / labels / colors =========

def safe_input_id(prefix: str, key: str) -> str:
    slug = re.sub(r"\W+", "_", key)
    h = hashlib.md5(key.encode("utf-8")).hexdigest()[:8]
    return f"{prefix}{slug}_{h}"

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

def lighten_hex(hex_color: str, factor: float = 0.45) -> str:
    """Lighten HEX by mixing with white."""
    hc = hex_color.lstrip("#")
    if len(hc) != 6:
        return hex_color
    r, g, b = int(hc[0:2], 16), int(hc[2:4], 16), int(hc[4:6], 16)
    r = int(r + (255 - r) * factor)
    g = int(g + (255 - g) * factor)
    b = int(b + (255 - b) * factor)
    return f"#{r:02x}{g:02x}{b:02x}"


# ========= Helpers: sample classification / parsing utilities =========

def clean_sample_id_remove_sbx(sample_id: str) -> str:
    return re.sub(r"-sbx$", "", str(sample_id), flags=re.I) if sample_id else sample_id

def classify_role(sample_id: str) -> str:
    """
    Return 'Tumor', 'Reference', or 'Unknown' from SampleID.
    - '-ref' token anywhere → Reference
    - 'h\\d+' start → Tumor (Illumina style)
    - trailing T/T2/TI… → Tumor; trailing R/R2/RI… → Reference
    """
    if not sample_id:
        return "Unknown"
    sid = re.sub(r"-(sbx|nseq6)$", "", str(sample_id).strip(), flags=re.I)
    if re.search(r"(?i)(?:^|[-_])ref(?:$|[-_])", sid):
        return "Reference"
    if re.match(r"(?i)^h\d{5,}", sid):
        return "Tumor"
    if re.search(r"(?i)T(?:\d+|I+)?$", sid):
        return "Tumor"
    if re.search(r"(?i)R(?:\d+|I+)?$", sid):
        return "Reference"
    return "Unknown"

def root_and_is_ref(name: str) -> Tuple[str, bool]:
    """Extract root H-number and reference flag from varied SampleID forms."""
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

def normalize_selected(values: Optional[Sequence[str]]) -> List[str]:
    out: List[str] = []
    for s in values or []:
        out.append(s.split(" — ", 1)[0] if " — " in s else s)
    return out

def extract_date_cycle_from_bam_uri(bam_uri: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    """Find YYYY-MM-DD (preferred) or YYMMDD, and 'cycleNN' tokens in a GCS path."""
    if not bam_uri:
        return None, None
    path_lower = bam_uri.lower()

    # ISO date
    iso = None
    m_iso = re.search(r"(\d{4})-(\d{2})-(\d{2})", path_lower)
    if m_iso:
        y, m, d = map(int, m_iso.groups())
        try:
            iso = dt.date(y, m, d).isoformat()
        except ValueError:
            iso = None

    # YYMMDD at segment start
    yymmdd = None
    if not iso:
        for seg in bam_uri.split("/"):
            m = re.match(r"(\d{6})", seg)
            if not m:
                continue
            yy, mm, dd = int(m.group(1)[:2]), int(m.group(1)[2:4]), int(m.group(1)[4:6])
            try:
                yymmdd = dt.date(2000 + yy, mm, dd).isoformat()
                break
            except ValueError:
                continue

    m_cycle = re.search(r"(cycle\d{1,3})", path_lower, flags=re.I)
    cycle = m_cycle.group(1) if m_cycle else None
    return (iso or yymmdd), cycle

def normalize_pct_columns_inplace(df: pd.DataFrame, cols: Sequence[str]) -> None:
    """Ensure %-like cols are on 0–1 scale."""
    for c in cols:
        if c not in df.columns:
            continue
        s = pd.to_numeric(df[c], errors="coerce")
        if s.dropna().empty:
            continue
        df[c] = s / 100.0 if 1.0 < float(s.max()) <= 100.0 else s


# ========= Helpers: GCS I/O =========

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


# ========= Helpers: parsing =========

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
    # simple fallback
    tail = dotted_key.split(".")[-1]
    for line in text.splitlines():
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        if k.strip() == tail:
            return v.strip().strip('"').strip("'")
    return fallback

def parse_picard_single_metrics(gcs: GCSFileSystem, file_path: str) -> Optional[pd.Series]:
    """Read first METRICS row under '## METRICS CLASS ...' block."""
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
        s = df.iloc[0]
        s.index = s.index.astype(str).str.strip()
        return s
    except Exception:
        return None

def parse_alignment_summary_histogram(gcs: GCSFileSystem, file_path: str) -> Optional[pd.DataFrame]:
    """
    Parse '## HISTOGRAM java.lang.Integer|Double' to standardized DF:
      columns: READ_LENGTH (numeric), TOTAL (optional numeric), ALIGNED (optional numeric)
    """
    try:
        with gcs.open(file_path, "r") as f:
            content = f.read()
        m = re.search(r"^##\s*HISTOGRAM\s+java\.lang\.(?:Integer|Double)\s*$", content, flags=re.M)
        if not m:
            dbg("aln_hist_header_not_found", file=file_path)
            return None
        tail = content[m.end():].lstrip("\n")
        next_section = re.search(r"^\#\#\s", tail, flags=re.M)
        hist_text = tail[: next_section.start()] if next_section else tail

        df = pd.read_csv(StringIO(hist_text), sep="\t", comment="#", header=0, engine="python")
        if df is None or df.empty:
            dbg("aln_hist_df_empty", file=file_path); return None

        # Column discovery
        def find_col(*need: str) -> Optional[str]:
            for c in df.columns:
                lc = c.lower()
                if all(s in lc for s in need):
                    return c
            return None

        read_len_col = next((c for c in df.columns if c.strip().upper() == "READ_LENGTH"), None)
        if not read_len_col:
            read_len_col = next((c for c in df.columns if ("read" in c.lower() and "length" in c.lower())), None)
        if not read_len_col:
            dbg("aln_hist_no_read_length_col", file=file_path, cols=list(df.columns)); return None

        total_col = find_col("total", "length", "count") or find_col("total", "reads") or find_col("total")
        aligned_col = find_col("aligned", "length", "count") or find_col("aligned", "reads") or find_col("aligned")

        out = pd.DataFrame()
        out["READ_LENGTH"] = pd.to_numeric(df[read_len_col], errors="coerce")
        if total_col and total_col in df:
            out["TOTAL"] = pd.to_numeric(df[total_col], errors="coerce")
        if aligned_col and aligned_col in df:
            out["ALIGNED"] = pd.to_numeric(df[aligned_col], errors="coerce")
        out = out.dropna(subset=["READ_LENGTH"])
        return out if not out.empty else None
    except Exception as e:
        dbg("aln_hist_parse_error", file=file_path, error=repr(e))
        return None

def parse_alignment_summary_table(gcs: GCSFileSystem, file_path: str) -> Optional[pd.Series]:
    """
    Parse '## METRICS CLASS picard.analysis.AlignmentSummaryMetrics' table.
    Prefer CATEGORY == 'PAIR', else 'UNPAIRED', else first row.
    """
    try:
        with gcs.open(file_path, "r") as f:
            content = f.read()
        m = re.search(r"^##\s+METRICS\s+CLASS\s+picard\.analysis\.AlignmentSummaryMetrics\s*$",
                      content, flags=re.M)
        if not m:
            dbg("aln_table_no_metrics_header", file=file_path); return None
        tail = content[m.end():]
        next_section = re.search(r"^\#\#\s", tail, flags=re.M)
        block = tail[: next_section.start()] if next_section else tail

        df = pd.read_csv(StringIO(block), sep="\t", comment="#", engine="python")
        if df is None or df.empty:
            dbg("aln_table_empty_df", file=file_path); return None

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

        s = row.copy()
        s.index = s.index.astype(str).str.strip()
        return s
    except Exception as e:
        dbg("aln_table_parse_error", file=file_path, error=repr(e))
        return None

def parse_runlevel_for_sample(gcs: GCSFileSystem, tsv_path: str, sample_id: str) -> Optional[Dict[str, float]]:
    """Extract relevant runlevel TSV metrics for the row that truly corresponds to sample_id."""
    try:
        with gcs.open(tsv_path, "r") as f:
            text = f.read()
    except Exception as e:
        dbg("runlevel_open_error", path=tsv_path, error=repr(e)); return None

    lines = text.splitlines()
    header_idx = next((i for i, ln in enumerate(lines) if "sample_name" in ln.split("\t")), None)
    if header_idx is None:
        dbg("runlevel_header_not_found", path=tsv_path); return None

    try:
        df = pd.read_csv(StringIO("\n".join(lines[header_idx:])), sep="\t")
    except Exception as e:
        dbg("runlevel_read_error", path=tsv_path, error=repr(e)); return None
    if df.empty or "sample_name" not in df.columns:
        dbg("runlevel_df_invalid", rows=int(df.shape[0])); return None

    drop_cols = [
        "run_id","Analysis Start","Analysis End","ppa2bam Start","ppa2bam End",
        "Sort Start","Sort End","BAM Start","BAM End","Diablo Start","Diablo End",
        "Time to Analyze (m)","pct_files_analyzed","SID_Number",
    ]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")
    df = df[~df["sample_name"].isna()].copy()
    if df.empty:
        return None

    df["__name__"] = df["sample_name"].astype(str)
    root_flag = df["__name__"].map(root_and_is_ref)
    df["__root__"] = root_flag.map(lambda x: x[0])
    df["__is_ref__"] = root_flag.map(lambda x: x[1])

    target_root, target_is_ref = root_and_is_ref(sample_id)

    # Step 1: exact (root, ref)
    sub = df[(df["__root__"] == target_root) & (df["__is_ref__"] == target_is_ref)]

    # Step 2: token-aware anchored regex
    if sub.empty and target_root:
        base = re.escape(target_root)
        if target_is_ref:
            pat = rf"(?i)^{base}(?:[-_]ref|[-_]R(?:\d+|I+)?)(?:[-_].*)?$"
        else:
            pat = rf"(?i)^{base}(?:(?:[-_]T(?:\d+|I+)?)?(?:[-_].*)?)$"
        df_pat = df[df["__name__"].str.match(pat, na=False)]
        sub = df_pat[df_pat["__is_ref__"] == target_is_ref]

    # Step 3: literal fallback
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

def runlevel_tsv_from_bam_uri(gcs: GCSFileSystem, bam_uri: Optional[str]) -> Optional[str]:
    """
    From gs://.../<run>/cycleXX/output/SAMPLE/... → gs://.../<run>/cycleXX/*.tsv
    """
    if not bam_uri:
        return None
    s = str(bam_uri).rstrip("/")
    m = re.search(r"/output(?:/|$)", s, flags=re.I)
    if not m:
        return None
    cycle_dir = s[:m.start()].rstrip("/")
    try:
        candidates = gcs.glob(f"{cycle_dir}/*.tsv")
        return sorted(candidates)[0] if candidates else None
    except Exception:
        return None

def parse_flagstat_duplicates(text: str) -> Optional[int]:
    m = re.search(r"^\s*(\d+)\s+\+\s+\d+\s+duplicates\s*$", text, flags=re.M | re.I)
    return int(m.group(1)) if m else None


# ========= Helpers: histogram stats =========

def read_length_n50_from_hist(hist_df: Optional[pd.DataFrame]) -> Optional[float]:
    if hist_df is None or hist_df.empty or "READ_LENGTH" not in hist_df.columns:
        return None
    count_col = "TOTAL" if "TOTAL" in hist_df.columns else ("ALIGNED" if "ALIGNED" in hist_df.columns else None)
    if not count_col:
        return None
    df = hist_df[["READ_LENGTH", count_col]].dropna().copy()
    df["READ_LENGTH"] = pd.to_numeric(df["READ_LENGTH"], errors="coerce")
    df[count_col] = pd.to_numeric(df[count_col], errors="coerce")
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

def median_read_length_from_hist(hist_df: Optional[pd.DataFrame]) -> Optional[float]:
    if hist_df is None or hist_df.empty or "READ_LENGTH" not in hist_df.columns:
        return None
    count_col = "TOTAL" if "TOTAL" in hist_df.columns else ("ALIGNED" if "ALIGNED" in hist_df.columns else None)
    if not count_col:
        return None
    df = hist_df[["READ_LENGTH", count_col]].dropna().copy()
    df["READ_LENGTH"] = pd.to_numeric(df["READ_LENGTH"], errors="coerce")
    df[count_col] = pd.to_numeric(df[count_col], errors="coerce")
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


# ========= UI =========

SAMPLE_SCREEN = ui.nav_panel(
    "screen_samples",
    ui.div(
        ui.h2("Metrics visualizer", style="text-align:center; font-weight:700; margin-bottom:8px;"),
        ui.p("Pick one or more batches below, then load data.",
             style="text-align:center; margin-top:0; margin-bottom:10px;"),
        ui.tags.style("""
        .feed-center{max-width:720px;margin:8px auto 0}
        .batch-feed-border{background:linear-gradient(135deg,#6a5acd,#7b68ee,#8a2be2);padding:2px;border-radius:14px}
        .batch-feed-wrapper{background:#fff;border-radius:12px;max-height:420px;overflow-y:auto;padding:8px 10px;box-shadow:0 4px 10px rgba(0,0,0,.07)}
        .batch-feed-wrapper::-webkit-scrollbar{width:8px}
        .batch-feed-wrapper::-webkit-scrollbar-thumb{background:rgba(0,0,0,.2);border-radius:6px}
        .batch-tiles .form-check{margin:0;padding:0}
        .batch-tiles .form-check-label{display:none}
        .batch-tiles .batch-item{cursor:pointer;}
        .batch-tiles .batch-item:has(input[type="checkbox"]:focus-visible){
        outline:2px solid #7b68ee; outline-offset:2px;
        }
        .batch-tiles .form-check-input{
        position:absolute;
        width:0; height:0;
        opacity:0;
        pointer-events:none;
        margin:0;
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
        }
        @media (max-width:520px){ .batch-tiles .batch-item{min-width:280px;width:95%} }
        """),
        ui.div(
            ui.div(ui.output_ui("batch_feed_items"), class_="batch-feed-wrapper batch-tiles"),
            class_="batch-feed-border feed-center"
        ),
        ui.div(
            ui.layout_columns(
                ui.input_action_button("batch_select_all", "Select all", class_="btn-sm"),
                ui.input_action_button("batch_clear", "Clear", class_="btn-sm"),
                ui.input_action_button("load_data", "Load Data & View Dashboard", class_="btn-primary"),
                col_widths=(2,2,8)
            ),
            style="max-width: 1100px; margin: 8px auto;"
        ),
        ui.div(
            ui.input_radio_buttons(
                "coverage_type", "Coverage data",
                choices={"cov": "General (MAPQ ≥ 20)",
                         "cov-hc":"High confidence regions statistics",
                         "cov-mq0":"All coverage statistics (MAPQ ≥ 0)"},
                selected="cov", inline=True,
            ),
            style="max-width: 1100px; margin: 6px auto;"
        ),
        style="padding: 20px 8px;"
    ),
)

DASHBOARD_SCREEN = ui.nav_panel(
    "screen_dashboard",
    ui.page_sidebar(
        ui.sidebar(
            ui.card(
                ui.h4("Batch selection"),
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
                ui.h4("Filter Samples"),
                ui.input_checkbox_group("filter_samples", "Displayed Samples", choices=[]),
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
        SAMPLE_SCREEN,
        DASHBOARD_SCREEN,
        id="wizard",
    ),
)


# ========= Data assembly helpers =========

def meta_df_or_empty(sample_meta_val: reactive.Value) -> pd.DataFrame:
    df = sample_meta_val.get()
    if df is None or not isinstance(df, pd.DataFrame):
        return pd.DataFrame(columns=["SampleID","Date","Cycle","bam_uri","yaml_found"])
    for c in ["SampleID","Date","Cycle","bam_uri","yaml_found"]:
        if c not in df.columns:
            df[c] = ""
    return df

def build_batch_choices(df_meta: pd.DataFrame) -> Dict[str, str]:
    if df_meta.empty:
        return {}
    df2 = df_meta.copy()
    df2["DateN"] = df2["Date"].apply(lambda x: x if str(x).strip() else None)
    df2["CycleN"] = df2["Cycle"].apply(lambda x: x if str(x).strip() else None)
    df2["_Date"] = pd.to_datetime(df2["DateN"], errors="coerce")

    def cycnum(c: Any) -> int:
        m = re.search(r"(\d+)", str(c) if c else "")
        return int(m.group(1)) if m else 10**9

    df2["_C"] = df2["CycleN"].map(cycnum)
    df2 = df2.sort_values(["_Date","_C","Platform"], ascending=[False, True, True], na_position="last")
    uniq = df2.drop_duplicates(subset=["Platform","DateN","CycleN"])[["Platform","DateN","CycleN"]]
    return {batch_key(r["Platform"], r["DateN"], r["CycleN"]): batch_label(r["Platform"], r["DateN"], r["CycleN"])
            for _, r in uniq.iterrows()}

def build_two_tone_color_map(df: pd.DataFrame) -> Tuple[Dict[str, str], List[str]]:
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


# ========= Plot builders =========

def make_metric_bar_figure(plot_df: pd.DataFrame, metric: str,
                           color_map: Mapping[str, str],
                           legend_order: Sequence[str]) -> go.Figure:
    axis_fmt = THRESHOLDS.get(metric, {}).get("format")
    hover_fmt = ":.2%" if (axis_fmt or "").endswith("%") else ":.4f"

    fig = px.bar(
        plot_df.sort_values(metric, na_position="last"),
        x="SampleID", y=metric,
        color="ColorKey",
        category_orders={"ColorKey": list(legend_order)},
        color_discrete_map=dict(color_map),
        hover_data={
            "SampleID": True, metric: hover_fmt, "Date": True, "Cycle": True,
            "Group": True, "Role": True, "ColorKey": False,
        },
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
    fig = px.bar(
        plot_df.sort_values(["SampleID","MetricLabel"]),
        x="SampleID", y="Value", color="MetricLabel",
        barmode="group",
        hover_data={"SampleID": True, "MetricLabel": True, "Value": ":.1%", "Date": True, "Cycle": True,
                    "Group": True, "Role": True},
        title="Duplex vs OnePlus (% of Total Reads)",
        labels={"SampleID": "Sample", "Value": "Percentage", "MetricLabel": "Type"},
    )
    fig.update_yaxes(tickformat=".0%")
    fig.update_layout(margin=dict(l=40, r=40, t=50, b=80), bargap=0.25, legend_title_text="Type")
    fig.update_xaxes(tickangle=-45)
    return fig

def make_read_length_figure(hist_by_sample: Dict[str, pd.DataFrame],
                            df_metrics: pd.DataFrame,
                            selected: Sequence[str],
                            color_map: Mapping[str, str]) -> go.Figure:
    key_map = {}
    if not df_metrics.empty:
        key_map = (df_metrics.drop_duplicates("SampleID")
                   .set_index("SampleID")["ColorKey"].to_dict())

    fig = go.Figure()
    any_curve = False
    seen_keys: set[str] = set()

    for sid in selected:
        hist = hist_by_sample.get(sid)
        if hist is None or hist.empty or "READ_LENGTH" not in hist.columns:
            continue
        color_key = key_map.get(sid, "Unknown")
        color = color_map.get(color_key)

        if color_key not in seen_keys:
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", name=color_key,
                                     line={"color": color} if color else None,
                                     hoverinfo="skip", showlegend=True))
            seen_keys.add(color_key)

        for series in ("TOTAL", "ALIGNED"):
            if series in hist.columns:
                fig.add_trace(go.Scatter(
                    x=hist["READ_LENGTH"], y=hist[series], mode="lines", name=f"{sid} — {series}",
                    showlegend=False, line={"color": color} if color else None,
                    hovertemplate=("Sample: %{customdata[0]}<br>Series: " + series +
                                   "<br>Read length: %{x}<br>Count: %{y}<extra></extra>"),
                    customdata=[[sid]] * len(hist),
                ))
                any_curve = True

    if not any_curve:
        return go.Figure(layout={"title": "No histogram data available for selected samples"})

    fig.update_layout(
        title="Read Length Distribution (AlignmentSummaryMetrics histogram)",
        xaxis_title="READ_LENGTH (bp)", yaxis_title="Count",
        margin=dict(l=50, r=30, t=60, b=50), legend_title="Date/Cycle/Platform — Role",
    )
    return fig


# ========= Server =========

def server(input, output, session):
    # --- Reactive state ---
    gcs = reactive.Value[Optional[GCSFileSystem]](None)
    sample_meta = reactive.Value[Optional[pd.DataFrame]](None)  # SampleID, Date, Cycle, bam_uri, yaml_found, Platform
    full_data = reactive.Value[Optional[pd.DataFrame]](None)
    readlen_data = reactive.Value[Dict[str, pd.DataFrame]]({})
    batch_id_map = reactive.Value[Dict[str, str]]({})

    role_color_map = reactive.Value[Dict[str, str]]({})
    color_key_order = reactive.Value[List[str]]([])

    # --- Bootstrap: connect GCS & scan bucket ---
    @reactive.Effect
    def _bootstrap():
        if gcs.get() is not None and sample_meta.get() is not None:
            return
        try:
            fs = GCSFileSystem()
            if not fs.exists(DEFAULT_BUCKET):
                ui.notification_show(f"Default bucket not found: {DEFAULT_BUCKET}", duration=7, type="error")
                return
            gcs.set(fs)
        except Exception as e:
            ui.notification_show(f"GCS connection failed: {e}", duration=8, type="error")
            return

        with ui.Progress(min=0, max=1) as p:
            p.set(0.1, message="Searching for samples...")
            samples = find_sample_folders(gcs.get(), DEFAULT_BUCKET)
            p.set(1)
        if not samples:
            ui.notification_show("No sample folders found in default bucket.", duration=6, type="warning")
            return

        rows: List[Dict[str, Any]] = []
        fs = gcs.get()
        with ui.Progress(min=0, max=max(1, len(samples))) as p:
            p.set(0, message="Reading sample metadata...", detail="execution-definition.yaml")
            for i, sid in enumerate(samples, start=1):
                ypath = f"{DEFAULT_BUCKET}/{sid}/execution-definition.yaml"
                text = read_text_if_exists(fs, ypath)
                bam_uri = parse_yaml_field(text, "params.bam_uri", fallback=None)
                date_iso, cycle = extract_date_cycle_from_bam_uri(bam_uri)
                rows.append({
                    "SampleID": sid,
                    "Date": date_iso or "",
                    "Cycle": cycle or "",
                    "bam_uri": bam_uri or "",
                    "yaml_found": bool(text),
                    "Platform": platform_from_sid(sid),
                })
                p.set(i, detail=sid)

        meta_df = pd.DataFrame(rows, columns=["SampleID","Date","Cycle","bam_uri","yaml_found","Platform"])
        meta_df["_Date"] = pd.to_datetime(meta_df["Date"].replace("", pd.NA), errors="coerce")
        meta_df = meta_df.sort_values("_Date", ascending=False, na_position="last").drop(columns=["_Date"])
        sample_meta.set(meta_df)
        ui.update_navset("wizard", selected="screen_samples")

    # --- Batch feed tiles ---
    @reactive.Calc
    def _batch_choices() -> Dict[str, str]:
        df = sample_meta.get()
        if df is None or df.empty:
            return {}
        return build_batch_choices(df)

    @render.ui
    def batch_feed_items():
        choices = _batch_choices()
        id_map: Dict[str, str] = {}
        items: List[Any] = []
        for key, lab in choices.items():
            sid = safe_input_id("bb_", key)
            id_map[sid] = key
            items.append(
                ui.tags.label(
                    ui.input_checkbox(id=sid, label="", value=False),
                    ui.tags.span(lab, class_="batch-text"),
                    ui.tags.span("", class_="tile-spacer"),
                    class_="batch-item",
                    **{"for": sid}  # not strictly required because input is nested, but harmless
                )
            )
        batch_id_map.set(id_map)
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

    # --- Populate cycles for sidebar ---
    @reactive.Effect
    def _populate_cycles():
        df_loaded = full_data.get()
        if df_loaded is not None and not df_loaded.empty and "Cycle" in df_loaded.columns:
            cycles = sorted({c.strip() for c in df_loaded["Cycle"].astype(str) if c and c.strip() and c.lower() != "nan"})
            dbg("populate_cycles", source="full_data", n=len(cycles), cycles=cycles[:50])
            # Keep prior selection if still valid; otherwise select all
            current = input.batch_cycles() or []
            new_sel = [c for c in current if c in cycles] or cycles
            ui.update_selectize("batch_cycles", choices=cycles, selected=new_sel)
            return

        df_meta = meta_df_or_empty(sample_meta)
        if df_meta.empty or "Cycle" not in df_meta.columns:
            dbg("populate_cycles", source="meta", status="no_cycles")
            ui.update_selectize("batch_cycles", choices=[], selected=[])
            return
        cycles = sorted({c.strip() for c in df_meta["Cycle"].astype(str) if c and c.strip() and c.lower() != "nan"})
        dbg("populate_cycles", source="meta", n=len(cycles), cycles=cycles[:50])
        current = input.batch_cycles() or []
        new_sel = [c for c in current if c in cycles] or cycles
        ui.update_selectize("batch_cycles", choices=cycles, selected=new_sel)

    # --- Selected samples derived from chosen batch tiles ---
    @reactive.Calc
    def _selected_samples_from_batches() -> List[str]:
        df = sample_meta.get()
        if df is None or df.empty:
            return []
        picked_keys: List[str] = []
        for safe_id, orig_key in batch_id_map.get().items():
            try:
                if getattr(input, safe_id)():
                    picked_keys.append(orig_key)
            except Exception:
                pass
        if not picked_keys:
            return []

        out: List[str] = []
        for key in picked_keys:
            platform, date_iso, cycle = key.split("|", 2)
            date_iso = None if date_iso == "unknown" else date_iso
            cycle    = None if cycle == "unknown" else cycle
            m = (
                (df["Platform"] == platform) &
                (df["Date"].apply(lambda x: x if str(x).strip() else None) == date_iso) &
                (df["Cycle"].apply(lambda x: x if str(x).strip() else None) == cycle)
            )
            out.extend(df.loc[m, "SampleID"].astype(str).tolist())
        return sorted(set(out))

    # --- Load data (main) ---
    @reactive.Effect
    @reactive.event(input.load_data)
    def _load_data():
        selected_samples = _selected_samples_from_batches()
        if not selected_samples:
            ui.notification_show("Pick at least one batch.", duration=5, type="warning")
            return

        meta_local = meta_df_or_empty(sample_meta)
        fs = gcs.get()
        bucket = DEFAULT_BUCKET

        series_rows: List[pd.Series] = []
        histograms: Dict[str, pd.DataFrame] = {}

        with ui.Progress(min=0, max=max(1, len(selected_samples))) as p:
            p.set(message="Loading metrics files...", detail="Starting...")
            for i, sample_id in enumerate(selected_samples, start=1):
                p.set(i-1, detail=f"Processing {sample_id}")

                # Picard (coverage/gc/insert/alignment)
                cov_folder = input.coverage_type()
                cov_files = fs.glob(f"{bucket}/{sample_id}/{cov_folder}/*.wgs.metrics.txt")
                gc_files  = fs.glob(f"{bucket}/{sample_id}/mmetrics/*gc_bias.summary_metrics*")
                ins_files = fs.glob(f"{bucket}/{sample_id}/mmetrics/*insert_size_metrics*")
                aln_files = fs.glob(f"{bucket}/{sample_id}/mmetrics/*alignment_summary_metrics*")

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

                # runlevel TSV attach
                row_meta = meta_local.loc[meta_local["SampleID"] == sample_id]
                bam_uri = str(row_meta["bam_uri"].iat[0]).strip() if not row_meta.empty else ""
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
                merged["READ_LENGTH_N50"] = read_length_n50_from_hist(hist_df)
                merged["mmetrics_median_read_length"] = median_read_length_from_hist(hist_df)

                # flagstat duplicates
                try:
                    flag_path = find_flagstat_txt(fs, bucket, sample_id)
                    if flag_path:
                        txt = read_text_if_exists(fs, flag_path) or ""
                        dups = parse_flagstat_duplicates(txt)
                        if dups is not None:
                            merged["FLAGSTAT_DUPLICATES"] = float(dups)
                except Exception as e:
                    dbg("flagstat_error", sample_id=sample_id, error=repr(e))

                if merged:
                    merged["SampleID"] = sample_id
                    series_rows.append(pd.Series(merged))

            p.set(max(1, len(selected_samples)), detail="Done.")

        if not series_rows:
            ui.notification_show("Failed to load data.", duration=5, type="error")
            return

        df = pd.DataFrame(series_rows)
        df.columns = df.columns.astype(str).str.strip()

        # numeric coercion for metrics we plot
        for col in METRICS_TO_PLOT:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        # extras
        extras = [
            "Phred_All","Phred_Substitutions","Phred_Insertions","Phred_Deletions","Phred_Indel",
            "median_read_length","mean_read_length","ratio_90_to_10percentile",
            "Num_Full_Length_Reads","Num_One_Plus_Reads","Total_Reads",
            "percent_duplex","percent_oneplus","READ_LENGTH_N50","TOTAL_READS",
            "mmetrics_median_read_length","PF_ALIGNED_BASES","FLAGSTAT_DUPLICATES",
        ]
        for c in extras:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        normalize_pct_columns_inplace(df, PCT_LIKE_COLUMNS)

        # add meta & display fields
        meta = meta_df_or_empty(sample_meta)[["SampleID","Date","Cycle","Platform"]].copy()
        df = df.merge(meta, on="SampleID", how="left")
        df["DateDisplay"] = df["Date"].apply(fmt_date_ddmmyyyy)
        df["Run"] = df["Cycle"].apply(run_from_cycle)
        df["Group"] = df.apply(lambda r: "Illumina" if r.get("Platform") == "Illumina"
                               else f"{r['DateDisplay']} • {r['Run']}", axis=1)
        df["Role"] = df["SampleID"].astype(str).map(lambda sid: classify_role(re.sub(r"-sbx\d*$", "", sid, flags=re.I)))
        df["ColorKey"] = df.apply(lambda r: f"{r['Group']} — {r['Role']}"
                                  if r["Role"] in ("Tumor","Reference") else r["Group"], axis=1)

        cmap, order = build_two_tone_color_map(df)
        role_color_map.set(cmap)
        existing = set(df["ColorKey"].dropna().astype(str).unique())
        color_key_order.set([k for k in order if k in existing])

        full_data.set(df)
        readlen_data.set(histograms)

        # DEBUG: summarize loaded dataset and filter fields
        _log_loaded_dataset_summary(df)

        try:
            dd = pd.to_datetime(df["Date"], errors="coerce")
            if dd.notna().any():
                dmin = dd.min().date()
                dmax = dd.max().date()
                dbg("init_date_range", min=str(dmin), max=str(dmax))
                ui.update_date_range("batch_date_range", start=dmin, end=dmax)
            # Default cycles to all present in full_data
            if "Cycle" in df.columns:
                cycles = sorted({c.strip() for c in df["Cycle"].astype(str) if c and c.strip() and c.lower() != "nan"})
                dbg("init_cycles", cycles=cycles)
                ui.update_selectize("batch_cycles", choices=cycles, selected=cycles)
        except Exception as e:
            dbg("init_filters_error", error=repr(e))

        # initialize filter choices
        sids = df["SampleID"].dropna().astype(str).tolist()
        ui.update_checkbox_group("filter_samples", choices={sid: sid for sid in sids}, selected=sids)
        ui.update_navset("wizard", selected="screen_dashboard")

    # --- Sidebar actions (select all/none) ---
    @reactive.Effect
    @reactive.event(input.filter_select_all)
    def _filter_all():
        df = full_data.get()
        if df is None:
            return
        sids = df["SampleID"].astype(str).tolist()
        ui.update_checkbox_group("filter_samples", selected=sids)

    @reactive.Effect
    @reactive.event(input.filter_deselect_all)
    def _filter_none():
        ui.update_checkbox_group("filter_samples", selected=[])

    # --- Batch filters (sidebar) ---
    def ids_matching_batch() -> List[str]:
        df = full_data.get()
        dbg("ids_matching_batch_called", have_full_data=(df is not None), rows=(0 if df is None else int(df.shape[0])))
        if df is None or df.empty:
            return []

        if not all(c in df.columns for c in ["SampleID","Date","Cycle","Role"]):
            dbg("ids_matching_batch_missing_cols", columns_present=list(df.columns))
            # Fall back to whatever we can
            needed = [c for c in ["SampleID","Date","Cycle","Role"] if c not in df.columns]
            dbg("ids_matching_batch_needed_missing", missing=needed)
            return []

        use = df.loc[:, ["SampleID","Date","Cycle","Role"]].copy()
        use["_Date"] = pd.to_datetime(use["Date"], errors="coerce")
        dbg("filter_step_start", rows=int(use.shape[0]),
            n_na_dates=int(use["_Date"].isna().sum()),
            unique_roles=sorted(use["Role"].astype(str).unique().tolist()),
            unique_cycles=sorted(use["Cycle"].astype(str).unique().tolist())[:50])

        # Log current inputs
        _log_filter_inputs("apply_pre", input)

        # Date range (inclusive). If outside loaded span → warn & return [].
        start, end = input.batch_date_range()
        if start is not None and end is not None:
            m = use["_Date"].notna() & (use["_Date"].dt.date >= start) & (use["_Date"].dt.date <= end)
            before = int(use.shape[0]); after = int(m.sum())
            # Compute loaded span for diagnostics
            dd_all = use["_Date"]
            have_dates = dd_all.notna().any()
            loaded_min = str(dd_all.min().date()) if have_dates else None
            loaded_max = str(dd_all.max().date()) if have_dates else None
            dbg("filter_date", start=str(start), end=str(end), before=before, after=after,
                loaded_min=loaded_min, loaded_max=loaded_max)
            use = use.loc[m]
            if after == 0:
                ui.notification_show(
                    f"No samples in date range {start} → {end}. "
                    f"Loaded data window is {loaded_min} → {loaded_max}.",
                    duration=6, type="warning"
                )
                # Early exit so Apply handler can show 0 cleanly
                return []

        # Cycle filter
        sel_cycles = [c for c in (input.batch_cycles() or []) if c]
        if sel_cycles:
            before = int(use.shape[0])
            use = use[use["Cycle"].astype(str).str.strip().isin(sel_cycles)]
            after = int(use.shape[0])
            dbg("filter_cycle", sel_cycles=sel_cycles, before=before, after=after)

        # Role filter (UI: Any/Tumor/Normal; data: Tumor/Reference)
        want_role = input.batch_role()
        if want_role and want_role != "Any":
            map_role = use["Role"].map({"Reference": "Normal", "Tumor": "Tumor"}).fillna("Unknown")
            before = int(use.shape[0])
            use = use.loc[map_role == want_role]
            after = int(use.shape[0])
            dbg("filter_role", want_role=want_role, before=before, after=after,
                role_counts_before=map_role.value_counts(dropna=False).to_dict())

        result = sorted(use["SampleID"].astype(str).unique().tolist())
        dbg("filter_result", n=len(result), sample_ids=result[:20])
        return result

    @reactive.Effect
    @reactive.event(input.batch_apply)
    def _batch_apply():
        _log_filter_inputs("apply_click", input)
        ids = ids_matching_batch()
        dbg("apply_click_result", n_ids=len(ids), ids_preview=ids[:20])
        if not ids:
            ui.notification_show("No samples match the chosen batch filters.", duration=5, type="warning")
            return
        df = full_data.get()
        if df is not None and not df.empty:
            existing = set(df["SampleID"].astype(str))
            ids = [s for s in ids if s in existing]
        ui.update_checkbox_group("filter_samples", selected=ids)

    @reactive.Effect
    @reactive.event(input.batch_reset)
    def _batch_reset():
        _log_filter_inputs("reset_click", input)
        ui.update_selectize("batch_cycles", selected=[])
        ui.update_radio_buttons("batch_role", selected="Any")
        ui.update_date_range("batch_date_range", start=None, end=None)
        df = full_data.get()
        all_ids = df["SampleID"].astype(str).tolist() if (df is not None and not df.empty) else []
        dbg("reset_click_result", n_ids=len(all_ids), ids_preview=all_ids[:20])
        if all_ids:
            ui.update_checkbox_group("filter_samples", selected=all_ids)

    # --- Dashboard content (tabs) ---
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
        tabs.append(ui.nav_panel("Read Length Distribution", output_widget("plot_read_length")))
        return ui.navset_card_tab(*tabs)

    # --- Download CSV (overview with average row) ---
    @render.download(filename=lambda: DOWNLOAD_FILENAME_FMT.format(dt.datetime.now()))
    def download_table_csv():
        df = full_data.get()
        if df is None or df.empty:
            yield "No data loaded\n".encode("utf-8"); return
        selected = normalize_selected(input.filter_samples())
        if not selected:
            yield "No samples selected\n".encode("utf-8"); return

        view = df[df["SampleID"].astype(str).isin(selected)].copy()
        role_series = view["Role"].fillna("Unknown").astype(str) if "Role" in view.columns else \
                      view["SampleID"].astype(str).map(lambda sid: classify_role(re.sub(r"-sbx\d*$", "", sid, flags=re.I)))
        map_role = {"Reference": "Normal", "Tumor": "Tumor"}
        view["SampleType"] = role_series.map(map_role).fillna("Unknown")

        if "percent_duplex" not in view.columns and {"Num_Full_Length_Reads","Total_Reads"}.issubset(view.columns):
            view["percent_duplex"] = pd.to_numeric(view["Num_Full_Length_Reads"], errors="coerce") / pd.to_numeric(view["Total_Reads"], errors="coerce")
        if "percent_oneplus" not in view.columns and {"Num_One_Plus_Reads","Total_Reads"}.issubset(view.columns):
            view["percent_oneplus"] = pd.to_numeric(view["Num_One_Plus_Reads"], errors="coerce") / pd.to_numeric(view["Total_Reads"], errors="coerce")

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
        avg_row = {c: "" for c in view.columns}
        avg_row["SampleID"] = "AVERAGE (selected)"
        avg_row["SampleType"] = "—"
        for c in numeric_cols: avg_row[c] = avg_vals.get(c, pd.NA)

        out_df = pd.concat([pd.DataFrame([avg_row])[view.columns], view], ignore_index=True)
        yield out_df.to_csv(index=False).encode("utf-8")

    # --- Table (same as download, rendered) ---
    @render.data_frame
    def table_overview():
        df = full_data.get()
        if df is None or df.empty:
            return render.DataTable(pd.DataFrame({"Info": ["No data loaded"]}), selection_mode="none")
        selected = normalize_selected(input.filter_samples())
        if not selected:
            return render.DataTable(pd.DataFrame({"Info": ["No samples selected"]}), selection_mode="none")

        view = df[df["SampleID"].astype(str).isin(selected)].copy()
        role_series = view["Role"].fillna("Unknown").astype(str) if "Role" in view.columns else \
                      view["SampleID"].astype(str).map(lambda sid: classify_role(re.sub(r"-sbx\d*$", "", sid, flags=re.I)))
        view["SampleType"] = role_series.map({"Reference": "Normal", "Tumor": "Tumor"}).fillna("Unknown")

        if "percent_duplex" not in view.columns and {"Num_Full_Length_Reads","Total_Reads"}.issubset(view.columns):
            view["percent_duplex"] = pd.to_numeric(view["Num_Full_Length_Reads"], errors="coerce") / pd.to_numeric(view["Total_Reads"], errors="coerce")
        if "percent_oneplus" not in view.columns and {"Num_One_Plus_Reads","Total_Reads"}.issubset(view.columns):
            view["percent_oneplus"] = pd.to_numeric(view["Num_One_Plus_Reads"], errors="coerce") / pd.to_numeric(view["Total_Reads"], errors="coerce")

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
        avg_row = {c: "" for c in view.columns}
        avg_row["SampleID"] = "AVERAGE (selected)"
        avg_row["SampleType"] = "—"
        for c in numeric_cols: avg_row[c] = avg_vals.get(c, pd.NA)

        final_cols = [c for c in cols_present if c in view.columns]
        view = view[final_cols].reset_index(drop=True)
        view = pd.concat([pd.DataFrame([avg_row])[final_cols], view], ignore_index=True)
        return render.DataTable(view, selection_mode="none")

    # --- Metric plots: register one output per metric using a factory to avoid closure gotchas ---
    def _make_metric_output(metric_name: str):
        @output(id=f"plot_{metric_name}")
        @render_widget
        def _plot():
            df = full_data.get()
            selected = normalize_selected(input.filter_samples())
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
        _make_metric_output(_metric)  # registers the output

    # Duplex vs OnePlus plot
    @output(id="plot_duplex_oneplus")
    @render_widget
    def _plot_duplex_oneplus():
        df = full_data.get()
        selected = normalize_selected(input.filter_samples())
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

    # Read length histogram plot
    @output(id="plot_read_length")
    @render_widget
    def _plot_read_length():
        selected = normalize_selected(input.filter_samples())
        hist_by_sample = readlen_data.get()
        if not hist_by_sample or not selected:
            return go.Figure(layout={"title": "No samples selected for Read Length Distribution"})
        df_loaded = full_data.get()
        if isinstance(df_loaded, pd.DataFrame) and not df_loaded.empty:
            df_metrics = df_loaded
        else:
            df_metrics = pd.DataFrame(columns=["SampleID","ColorKey"])
        return make_read_length_figure(hist_by_sample, df_metrics, selected, role_color_map.get() or {})


# ========= Run =========
app = App(APP_UI, server)