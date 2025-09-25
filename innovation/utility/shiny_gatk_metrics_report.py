#!/usr/bin/env python3
import re
from io import StringIO
import datetime as dt
import sys

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from gcsfs import GCSFileSystem
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_widget
import yaml

# ==============================
# Config
# ==============================
METRICS_TO_PLOT = [
    "GENOME_TERRITORY","MEAN_COVERAGE","SD_COVERAGE","MEDIAN_COVERAGE","MAD_COVERAGE",
    "PCT_EXC_ADAPTER","PCT_EXC_MAPQ","PCT_EXC_DUPE","PCT_EXC_UNPAIRED","PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP","PCT_EXC_CAPPED","PCT_EXC_TOTAL","PCT_1X","PCT_5X","PCT_10X","PCT_15X",
    "PCT_20X","PCT_25X","PCT_30X","PCT_40X","PCT_50X","PCT_60X","PCT_70X","PCT_80X","PCT_90X",
    "PCT_100X","FOLD_80_BASE_PENALTY","FOLD_90_BASE_PENALTY","FOLD_95_BASE_PENALTY",
    "HET_SNP_SENSITIVITY","HET_SNP_Q","AT_DROPOUT","GC_DROPOUT",
    "Phred_All","Phred_Substitutions","Phred_Insertions","Phred_Deletions","Phred_Indel",
    "median_read_length","mean_read_length","ratio_90_to_10percentile"
]

PALETTE = [
    "#e41a1c", "#ff7f00", "#377eb8", "#4daf4a", "#984ea3",
    "#a65628", "#f781bf", "#999999", "#dede00", "#66c2a5",
    "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
    "#e5c494", "#b3b3b3",
]

THRESHOLDS = {
    "MEAN_COVERAGE": {"lines": [
        {"value": 30, "color": "orange", "label": "Reference Min (30x)"},
        {"value": 90, "color": "red", "label": "Tumor Min (90x)"},
    ]},
    "MEDIAN_COVERAGE": {"lines": [
        {"value": 30, "color": "orange", "label": "Reference Min (30x)"},
        {"value": 90, "color": "red", "label": "Tumor Min (90x)"},
    ]},
    "PCT_30X": {"lines": [{"value": 0.95, "color": "red"}], "format": ".0%"},
    "AT_DROPOUT": {"lines": [], "format": ".1%"},
    "GC_DROPOUT": {"lines": [], "format": ".1%"},
}

def dbg(tag: str, **kv):
    ts = dt.datetime.now().strftime("%H:%M:%S")
    parts = []
    for k, v in kv.items():
        # keep lists short in logs
        if isinstance(v, (list, tuple)) and len(v) > 10:
            parts.append(f"{k}=[{', '.join(map(repr, v[:10]))}, … +{len(v)-10}]")
        else:
            parts.append(f"{k}={repr(v)}")
    print(f"[{ts}] {tag} " + " ".join(parts), file=sys.stderr, flush=True)

# ==============================
# Helpers
# ==============================

def _get_meta_df(sample_meta) -> pd.DataFrame:
    """Return a real DataFrame with the expected columns, never None."""
    df = sample_meta.get()
    if df is None or not isinstance(df, pd.DataFrame):
        return pd.DataFrame(columns=["SampleID","Date","Cycle","bam_uri","yaml_found"])
    # Ensure required columns exist
    for c in ["SampleID","Date","Cycle","bam_uri","yaml_found"]:
        if c not in df.columns:
            df[c] = ""
    return df

def _clean_sample_id_remove_sbx(sample_id: str) -> str:
    """Strip a trailing '-sbx' (case-insensitive). Do NOT touch '-ref'."""
    if not sample_id:
        return sample_id
    return re.sub(r"-sbx$", "", str(sample_id), flags=re.I)

def _metrics_report_path_from_bam_uri(bam_uri: str | None) -> str | None:
    """
    bam_uri example:
      gs://.../cycle01/output/SAMPLE001/
    report path:
      gs://.../cycle01/output/metrics-report.tsv
    """
    if not bam_uri:
        dbg("metrics_report_path", bam_uri=None)
        return None
    s = bam_uri.rstrip("/")
    parent1 = s.rsplit("/", 1)[0] if "/" in s else s          # up from SAMPLEID
    parent2 = parent1.rsplit("/", 1)[0] if "/" in parent1 else parent1  # maybe up to /output
    base = parent1 if parent1.endswith("/output") else parent2
    report = f"{base}/metrics-report.tsv" if base else None
    dbg("metrics_report_path", bam_uri=bam_uri, parent1=parent1, parent2=parent2, chosen_base=base, report=report)
    return report

def _parse_metrics_report_for_sample(gcs: GCSFileSystem, report_path: str, sample_id: str) -> dict | None:
    """
    Read metrics-report.tsv and return selected fields for the matching row.
    Matching uses 'sample_name' CONTAINS cleaned SampleID (strip trailing -sbx).
    """
    if not report_path:
        dbg("metrics_report_skip", reason="no_report_path", sample_id=sample_id)
        return None

    try:
        if not gcs.exists(report_path):
            dbg("metrics_report_missing", report_path=report_path, sample_id=sample_id)
            return None

        with gcs.open(report_path, "r") as f:
            df = pd.read_csv(f, sep="\t")

        dbg("metrics_report_loaded",
            report_path=report_path,
            rows=int(df.shape[0]),
            cols=list(df.columns)[:20])

    except Exception as e:
        dbg("metrics_report_read_error", report_path=report_path, error=repr(e))
        return None

    if df.empty or "sample_name" not in df.columns:
        dbg("metrics_report_bad_df", empty=df.empty, has_sample_name=("sample_name" in df.columns))
        return None

    target = _clean_sample_id_remove_sbx(sample_id)
    # primary: case-sensitive contains
    mask_cs = df["sample_name"].astype(str).str.contains(re.escape(target), regex=True)
    n_cs = int(mask_cs.sum())

    # fallback: case-insensitive if no case-sensitive match
    mask_ci = None
    n_ci = 0
    if n_cs == 0:
        mask_ci = df["sample_name"].astype(str).str.contains(re.escape(target), regex=True, case=False)
        n_ci = int(mask_ci.sum())

    dbg("metrics_report_match_stats",
        sample_id=sample_id, cleaned_target=target, case_sensitive_hits=n_cs, case_insensitive_hits=n_ci)

    if n_cs > 0:
        sub = df.loc[mask_cs].copy()
    elif n_ci > 0:
        dbg("metrics_report_using_case_insensitive", sample_id=sample_id, target=target)
        sub = df.loc[mask_ci].copy()
    else:
        # show a few sample_name values for quick eyeballing
        dbg("metrics_report_no_match",
            sample_id=sample_id,
            target=target,
            sample_name_head=df["sample_name"].astype(str).head(5).tolist())
        return None

    if len(sub) > 1:
        dbg("metrics_report_multimatch",
            sample_id=sample_id,
            target=target,
            first_5=sub["sample_name"].astype(str).head(5).tolist(),
            total=int(len(sub)))

    row = sub.iloc[0]

    keep = [
        "Phred_All","Phred_Substitutions","Phred_Insertions","Phred_Deletions","Phred_Indel",
        "median_read_length","mean_read_length","ratio_90_to_10percentile",
    ]
    out = {}
    for k in keep:
        if k in row.index:
            out[k] = pd.to_numeric(row[k], errors="coerce")

    dbg("metrics_report_selected_values",
        sample_id=sample_id,
        selected={k: out.get(k, None) for k in keep})

    return out if out else None

def _normalize_pct_columns_inplace(df: pd.DataFrame, cols: list[str]) -> None:
    """
    Ensure percent-like columns are on 0–1 scale for plotting with % tickformat.
    If values look like 0–100, divide by 100 (but do nothing if already 0–1).
    """
    for c in cols:
        if c not in df.columns:
            continue
        s = pd.to_numeric(df[c], errors="coerce")
        if s.dropna().empty:
            continue
        maxv = float(s.max())
        if 1.0 < maxv <= 100.0:
            df[c] = s / 100.0
        else:
            df[c] = s  # already 0–1 (or weird outlier >100, leave as-is)


def _normalize_selected(sel: list[str]) -> list[str]:
    """Convert checkbox selections to pure SampleIDs.
    Works for both correct mapping (values are SampleIDs) and old/bad mapping (labels with ' — ')."""
    out = []
    for s in sel or []:
        if " — " in s:
            out.append(s.split(" — ", 1)[0])  # take leading SampleID
        else:
            out.append(s)
    return out

def _extract_date_cycle_from_bam_uri(bam_uri: str | None) -> tuple[str | None, str | None]:
    """
    Detect:
      - YYYY-MM-DD anywhere (preferred)
      - YYMMDD at start of any segment (assume 20YY)
      - 'cycleNN' anywhere
    Return (YYYY-MM-DD or None, cycleNN or None)
    """
    if not bam_uri:
        return None, None

    path_lower = bam_uri.lower()

    # ISO date YYYY-MM-DD
    iso = None
    m_iso = re.search(r"(\d{4})-(\d{2})-(\d{2})", path_lower)
    if m_iso:
        y, m, d = map(int, m_iso.groups())
        try:
            iso = dt.date(y, m, d).isoformat()
        except ValueError:
            iso = None

    # Fallback YYMMDD at start of segment
    yymmdd = None
    if not iso:
        for seg in bam_uri.split("/"):
            m = re.match(r"(\d{6})", seg)
            if m:
                raw = m.group(1)
                yy, mm, dd = int(raw[:2]), int(raw[2:4]), int(raw[4:6])
                yyyy = 2000 + yy
                try:
                    yymmdd = dt.date(yyyy, mm, dd).isoformat()
                    break
                except ValueError:
                    continue

    # cycleNN
    m_cycle = re.search(r"(cycle\d{1,3})", path_lower, flags=re.I)
    cycle = m_cycle.group(1) if m_cycle else None

    return (iso or yymmdd), cycle


def _parse_yaml_field(text: str | None, dotted_key: str, fallback: str | None = None) -> str | None:
    if text is None:
        return fallback
    try:
        data = yaml.safe_load(text) or {}
        cur = data
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


def _read_text_if_exists(gcs: GCSFileSystem, path: str) -> str | None:
    try:
        if gcs.exists(path):
            with gcs.open(path, "r") as f:
                return f.read()
    except Exception:
        pass
    return None


def find_sample_folders(gcs: GCSFileSystem, bucket_path: str) -> list[str]:
    try:
        all_paths = gcs.ls(bucket_path, detail=False)
        return [p.rsplit("/", 1)[-1] for p in all_paths if gcs.isdir(p)]
    except Exception:
        return []


def parse_picard_single_metrics(gcs: GCSFileSystem, file_path: str) -> pd.Series | None:
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


def parse_alignment_summary_histogram(gcs: GCSFileSystem, file_path: str) -> pd.DataFrame | None:
    try:
        with gcs.open(file_path, "r") as f:
            content = f.read()
        m = re.search(r"^## HISTOGRAM\s+java\.lang\.Integer\s*$", content, flags=re.M)
        if not m:
            return None
        start = m.end()
        tail = content[start:].lstrip("\n")
        next_section = re.search(r"^\#\#\s", tail, flags=re.M)
        hist_text = tail[: next_section.start()] if next_section else tail

        df = pd.read_csv(StringIO(hist_text), sep="\t", comment="#", header=0, engine="python")
        df.columns = df.columns.str.strip()
        for c in ["READ_LENGTH","UNPAIRED_TOTAL_LENGTH_COUNT","UNPAIRED_ALIGNED_LENGTH_COUNT"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=["READ_LENGTH"])
        return df
    except Exception:
        return None

# ==============================
# UI
# ==============================
screen1_bucket_input = ui.nav_panel(
    "screen_1_bucket",
    ui.div(
        ui.h2("GATK WGS Metrics Dashboard"),
        ui.p("Enter the path to your GCS bucket containing the sample folders."),
        ui.input_text("bucket_path", "GCS Bucket Path", "gs://gatk-wgs-metrics"),
        ui.input_action_button("submit_bucket", "Find Samples", class_="btn-primary"),
        style="max-width: 600px; margin: auto; padding-top: 50px;",
    ),
)

screen2_sample_selection = ui.nav_panel(
    "screen_2_samples",
    ui.div(
        ui.h2("Select Samples"),
        ui.p("Pick a date range and cycles; then exclude specific rows if needed."),
        ui.layout_columns(
            ui.card(
                ui.h5("Criteria"),
                ui.input_date_range("date_range", "Date range", start=None, end=None),
                ui.input_selectize("cycle_multi", "Cycles", choices=[], multiple=True, options={"placeholder": "All cycles"}),
                ui.output_ui("selection_summary"),
            ),
            ui.card(
                ui.h5("Per-sample overrides"),
                ui.p("Select rows below, then Exclude / Re-include."),
                ui.input_action_button("exclude_mark", "Exclude selected rows", class_="btn-warning"),
                ui.input_action_button("exclude_unmark", "Re-include selected rows", class_="btn-secondary", style="margin-left:8px"),
                ui.input_action_button("exclude_clear", "Clear all exclusions", class_="btn-link", style="margin-left:8px"),
                ui.output_text("excluded_list"),
            ),
            col_widths=(6, 6),
        ),
        ui.hr(),
        ui.h4("Samples"),
        ui.output_data_frame("sample_table"),
        ui.div(
            ui.input_action_button("load_data", "Load Data & View Dashboard", class_="btn-primary"),
            style="margin-top: 16px;"
        ),
        style="max-width: 1100px; margin: auto; padding-top: 20px;",
    ),
)

screen3_dashboard = ui.nav_panel(
    "screen_3_dashboard",
    ui.page_sidebar(
        ui.sidebar(
            ui.h4("Filter Samples"),
            ui.input_checkbox_group("filter_samples", "Displayed Samples", choices=[]),
            ui.input_action_button("filter_select_all", "Select All", class_="btn-sm"),
            ui.input_action_button("filter_deselect_all", "Deselect All", class_="btn-sm"),
            ui.hr(),
            ui.h5("Displayed metadata"),
            ui.output_data_frame("filter_meta_table"),
        ),
        ui.output_ui("dashboard_content"),
    ),
)

app_ui = ui.page_fluid(
    ui.navset_hidden(
        screen1_bucket_input,
        screen2_sample_selection,
        screen3_dashboard,
        id="wizard",
    ),
)

# ==============================
# Server
# ==============================
def server(input, output, session):
    # --- reactive state ---
    gcs = reactive.Value[GCSFileSystem | None](None)
    sample_meta = reactive.Value[pd.DataFrame | None](None)  # SampleID, Date, Cycle, bam_uri, yaml_found
    full_data = reactive.Value[pd.DataFrame | None](None)    # metrics joined per SampleID
    readlen_data = reactive.Value[dict[str, pd.DataFrame]]({})  # hist df per SampleID
    excluded_ids = reactive.Value(set())                     # exclusions in Screen 2
    selected_initialized = reactive.Value(False)             # one-shot select-all on Screen 3

    # ---------- Screen 1: scan bucket ----------
    @reactive.Effect
    @reactive.event(input.submit_bucket)
    def handle_bucket_submission():
        bucket_path = input.bucket_path().strip()
        if not bucket_path.startswith("gs://"):
            ui.notification_show("Invalid path. Must start with 'gs://'.", duration=5, type="error")
            return
        try:
            fs = GCSFileSystem()
            if not fs.exists(bucket_path):
                ui.notification_show("Bucket path not found.", duration=5, type="error")
                return
            gcs.set(fs)
        except Exception as e:
            ui.notification_show(f"GCS connection failed: {e}", duration=8, type="error")
            return

        with ui.Progress(min=0, max=1) as p:
            p.set(0.1, message="Searching for samples...")
            samples = find_sample_folders(gcs.get(), bucket_path)
            p.set(1)

        if not samples:
            ui.notification_show("No sample folders found.", duration=5, type="warning")
            return

        # Gather metadata
        with ui.Progress(min=0, max=max(1, len(samples))) as p:
            p.set(0, message="Reading sample metadata...", detail="execution-definition.yaml")
            rows = []
            fs = gcs.get()
            for i, sid in enumerate(samples, start=1):
                ypath = f"{bucket_path}/{sid}/execution-definition.yaml"
                text = _read_text_if_exists(fs, ypath)
                bam_uri = _parse_yaml_field(text, "params.bam_uri", fallback=None)
                date_str, cycle = _extract_date_cycle_from_bam_uri(bam_uri)
                rows.append({
                    "SampleID": sid,
                    "Date": date_str or "",
                    "Cycle": cycle or "",
                    "bam_uri": bam_uri or "",
                    "yaml_found": bool(text),
                })
                p.set(i, detail=sid)

        meta_df = pd.DataFrame(rows, columns=["SampleID","Date","Cycle","bam_uri","yaml_found"])
        meta_df["_sort_date"] = pd.to_datetime(meta_df["Date"], errors="coerce")
        meta_df = meta_df.sort_values("_sort_date", ascending=False, na_position="last").drop(columns=["_sort_date"])
        sample_meta.set(meta_df)

        # Initialize criteria widgets
        dt_col = pd.to_datetime(meta_df["Date"].replace("", pd.NA), errors="coerce")
        min_dt = pd.to_datetime(dt_col.min()) if not dt_col.isna().all() else None
        max_dt = pd.to_datetime(dt_col.max()) if not dt_col.isna().all() else None
        cycles = sorted([c for c in meta_df["Cycle"].unique().tolist() if c])

        ui.update_date_range("date_range", start=(min_dt.date() if min_dt is not None else None),
                                             end=(max_dt.date() if max_dt is not None else None))
        ui.update_selectize("cycle_multi", choices=cycles, selected=[])  # empty = all

        # reset exclusions
        excluded_ids.set(set())

        ui.update_navs("wizard", selected="screen_2_samples")
        ui.notification_show(f"Found {len(samples)} potential samples.", duration=5)

    # ---------- Screen 2: criteria + exclusions ----------
    @reactive.Calc
    def meta_table() -> pd.DataFrame:
        df = sample_meta.get()
        return df if df is not None else pd.DataFrame(columns=["SampleID","Date","Cycle","bam_uri","yaml_found"])

    @reactive.Calc
    def criteria_filtered_ids() -> list[str]:
        df = meta_table().copy()
        if df.empty:
            return []

        # Parse date column once
        df["_Date"] = pd.to_datetime(df["Date"].replace("", pd.NA), errors="coerce")

        # --- Date filter: include unknown dates as well ---
        start, end = input.date_range()
        if start is not None and end is not None:
            in_range = (df["_Date"].notna()) & (df["_Date"].dt.date >= start) & (df["_Date"].dt.date <= end)
            unknown_date = df["_Date"].isna()
            df = df[in_range | unknown_date]   # <— keep unknown dates

        # --- Cycle filter: include unknown cycles as well ---
        sel_cycles = input.cycle_multi()
        if sel_cycles:
            # cycle is stored as empty string "" when unknown
            cycle_str = df["Cycle"].astype(str).str.strip()
            in_cycle = cycle_str.isin(sel_cycles)
            unknown_cycle = (cycle_str == "") | df["Cycle"].isna()
            df = df[in_cycle | unknown_cycle]  # <— keep unknown cycles

        return df["SampleID"].astype(str).tolist()

    @reactive.Calc
    def final_selected_ids() -> list[str]:
        base = set(criteria_filtered_ids())
        excl = excluded_ids.get()
        return sorted(list(base - excl))

    def _selected_rows_to_ids() -> list[str]:
        rows = input.sample_table_selected_rows()
        df = meta_table()
        if not rows or df.empty:
            return []
        return df.iloc[list(rows)]["SampleID"].astype(str).tolist()

    @reactive.Effect
    @reactive.event(input.exclude_mark)
    def _exclude_mark():
        ids = set(_selected_rows_to_ids())
        if not ids:
            return
        new = set(excluded_ids.get())
        new.update(ids)
        excluded_ids.set(new)

    @reactive.Effect
    @reactive.event(input.exclude_unmark)
    def _exclude_unmark():
        ids = set(_selected_rows_to_ids())
        if not ids:
            return
        new = set(excluded_ids.get())
        new.difference_update(ids)
        excluded_ids.set(new)

    @reactive.Effect
    @reactive.event(input.exclude_clear)
    def _exclude_clear():
        excluded_ids.set(set())

    @render.data_frame
    def sample_table():
        df = meta_table()
        if df.empty:
            return pd.DataFrame({"Info": ["No metadata found"]})
        # User selects rows to exclude/include
        return render.DataTable(df, selection_mode="rows")

    @render.ui
    def selection_summary():
        base = criteria_filtered_ids()
        excl = excluded_ids.get()
        final_n = len(set(base) - excl)
        return ui.div(
            ui.tags.div(f"Matching criteria: {len(base)} samples"),
            ui.tags.div(f"Excluded: {len(excl)} samples"),
            ui.tags.div(ui.tags.b(f"Will load: {final_n} samples")),
            class_="text-muted", style="margin-top:8px"
        )

    @render.text
    def excluded_list():
        ids = sorted(list(excluded_ids.get()))
        if not ids:
            return "No samples excluded."
        short = ", ".join(ids[:12])
        more = f" … (+{len(ids)-12} more)" if len(ids) > 12 else ""
        return f"Excluded: {short}{more}"

    # ---------- Load data and move to Screen 3 ----------
    @reactive.Effect
    @reactive.event(input.load_data)
    def load_metrics_data():
        selected_samples = final_selected_ids()
        if not selected_samples:
            ui.notification_show("No samples match the criteria after exclusions.", duration=5, type="warning")
            return

        fs = gcs.get()
        bucket = input.bucket_path().strip()
        all_metrics: list[pd.Series] = []
        histograms: dict[str, pd.DataFrame] = {}

        with ui.Progress(min=0, max=max(1, len(selected_samples))) as p:
            p.set(message="Loading metrics files...", detail="Starting...")
            for i, sample_id in enumerate(selected_samples, start=1):
                p.set(i-1, detail=f"Processing {sample_id}")

                # 1) WGS metrics
                cov_glob = f"{bucket}/{sample_id}/cov/*.wgs.metrics.txt"
                cov_files = fs.glob(cov_glob)
                cov_row = parse_picard_single_metrics(fs, cov_files[0]) if cov_files else None

                # 2) GC bias summary
                gc_glob = f"{bucket}/{sample_id}/mmetrics/*gc_bias.summary_metrics*"
                gc_files = fs.glob(gc_glob)
                gc_row = parse_picard_single_metrics(fs, gc_files[0]) if gc_files else None

                merged = {}
                if cov_row is not None:
                    merged.update({str(k).strip(): v for k, v in cov_row.to_dict().items()})
                if gc_row is not None:
                    merged.update({str(k).strip(): v for k, v in gc_row.to_dict().items()})
                if merged:
                    merged["SampleID"] = sample_id
                    all_metrics.append(pd.Series(merged))

                # 2.5) metrics-report.tsv next to bam_uri (if available)
                #     Use the cleaned SampleID to match 'sample_name'
                # Get the latest stored metadata DataFrame
                meta = _get_meta_df(sample_meta)

                row = meta.loc[meta["SampleID"] == sample_id]
                bam_uri = None
                if not row.empty:
                    val = str(row["bam_uri"].iat[0]).strip()
                    bam_uri = val if val else None

                dbg("metrics_report_probe_start", sample_id=sample_id, bam_uri=bam_uri)

                report_path = _metrics_report_path_from_bam_uri(bam_uri)

                # quick visibility into the /output directory contents
                if report_path:
                    base_dir = report_path.rsplit("/", 1)[0]
                    try:
                        exists = fs.exists(base_dir)
                    except Exception:
                        exists = False
                    dbg("metrics_report_base_dir", base_dir=base_dir, exists=exists)

                    if exists:
                        try:
                            # list up to 12 entries for debugging
                            entries = fs.ls(base_dir)
                            dbg("metrics_report_dir_list",
                                base_dir=base_dir,
                                n_entries=len(entries),
                                head=entries[:12])
                        except Exception as e:
                            dbg("metrics_report_dir_list_error", base_dir=base_dir, error=repr(e))

                if report_path:
                    report_vals = _parse_metrics_report_for_sample(fs, report_path, sample_id)

                if report_vals:
                    dbg("metrics_report_merge", sample_id=sample_id, keys=list(report_vals.keys()))
                    merged.update(report_vals)
                else:
                    dbg("metrics_report_skip_merge", sample_id=sample_id)

                # 3) Alignment summary histogram
                aln_glob = f"{bucket}/{sample_id}/mmetrics/*alignment_summary_metrics*"
                aln_files = fs.glob(aln_glob)
                if aln_files:
                    hist_df = parse_alignment_summary_histogram(fs, aln_files[0])
                    if hist_df is not None and not hist_df.empty:
                        histograms[sample_id] = hist_df

            p.set(max(1, len(selected_samples)), detail="Done.")

        if not all_metrics:
            ui.notification_show("Failed to load data.", duration=5, type="error")
            return

        df = pd.DataFrame(all_metrics)
        df.columns = df.columns.astype(str).str.strip()
        for col in METRICS_TO_PLOT:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        
        dbg("metrics_df_built", shape=df.shape, cols=list(df.columns))

        # Normalize percent-like GC-bias metrics to fractions for proper % axes
        _normalize_pct_columns_inplace(df, ["AT_DROPOUT", "GC_DROPOUT"])

        full_data.set(df)
        readlen_data.set(histograms)

        # Build sidebar labels with Date/Cycle
        sample_ids = df["SampleID"].dropna().astype(str).tolist()
        dbg("loaded_metrics_df", df_shape=df.shape, df_cols=list(df.columns), n_samples=len(sample_ids))

        meta = sample_meta.get()
        if meta is None:
            meta = pd.DataFrame(columns=["SampleID", "Date", "Cycle"])

        df = df.merge(meta[["SampleID", "Date", "Cycle"]], on="SampleID", how="left")
        df["Date"] = df["Date"].fillna("")
        df["Cycle"] = df["Cycle"].fillna("")
        df["Group"] = (df["Date"].replace("", "?") + " • " + df["Cycle"].replace("", "?"))

        # ---- Build deterministic order and color map for groups ----
        tmp = df[["Date", "Cycle", "Group"]].drop_duplicates().copy()
        tmp["_Date"] = pd.to_datetime(tmp["Date"], errors="coerce")

        def _cycle_num(c: str) -> int:
            m = re.search(r"(\d+)", str(c) if c else "")
            return int(m.group(1)) if m else -1  # unknown cycles sort first

        tmp["_CycleNum"] = tmp["Cycle"].map(_cycle_num)

        ordered_groups = (
            tmp.sort_values(["_Date", "_CycleNum", "Group"], ascending=[True, True, True])
            ["Group"].dropna().tolist()
        )
        color_map = {grp: PALETTE[i % len(PALETTE)] for i, grp in enumerate(ordered_groups)}

        group_order.set(ordered_groups)
        group_color_map.set(color_map)

        # Persist full dataset with metadata for plotting
        full_data.set(df)
        readlen_data.set(histograms)

        # values (SampleIDs) -> labels (pretty text)
        choices = {}
        for sid in sample_ids:
            row = meta.loc[meta["SampleID"] == sid]
            date = row["Date"].iat[0] if not row.empty else ""
            cycle = row["Cycle"].iat[0] if not row.empty else ""
            choices[sid] = f"{sid} — {date or '?'} — {cycle or '?'}"

        dbg("update_filter_choices", n_choices=len(choices))
        ui.update_checkbox_group("filter_samples", choices=choices, selected=sample_ids)
        dbg("update_filter_selected", selected_count=len(sample_ids), selected_head=sample_ids[:10])

        ui.update_navs("wizard", selected="screen_3_dashboard")

    # ---------- Screen 3: sidebar behaviors ----------

    group_order = reactive.Value[list[str]]([])
    group_color_map = reactive.Value[dict[str, str]]({})

    @reactive.Effect
    def _debug_filter_selection_changes():
        sel_raw = input.filter_samples()
        sel = _normalize_selected(sel_raw)
        dbg("filter_samples_changed", n=len(sel), sel_head=tuple(sel[:10]))

    @reactive.Effect
    @reactive.event(input.filter_select_all)
    def _filter_all():
        dbg("filter_select_all_clicked")
        df = full_data.get()
        if df is not None:
            sids = df["SampleID"].astype(str).tolist()
            dbg("filter_select_all_set", count=len(sids), head=sids[:10])
            ui.update_checkbox_group("filter_samples", selected=sids)
        else:
            dbg("filter_select_all_skipped", reason="full_data None")

    @reactive.Effect
    @reactive.event(input.filter_deselect_all)
    def _filter_none():
        dbg("filter_deselect_all_clicked")
        ui.update_checkbox_group("filter_samples", selected=[])

    @render.data_frame
    def filter_meta_table():
        df_all = sample_meta.get()
        if df_all is None or df_all.empty:
            dbg("filter_meta_table", note="no metadata")
            return pd.DataFrame({"Info": ["No metadata available"]})

        shown_raw = input.filter_samples()
        shown = _normalize_selected(shown_raw)
        dbg("filter_meta_table", shown_n=len(shown), shown_head=shown[:10])

        if not shown:
            return render.DataTable(pd.DataFrame(columns=["SampleID","Date","Cycle"]), selection_mode="none")

        view = df_all.loc[df_all["SampleID"].isin(shown), ["SampleID","Date","Cycle"]].copy()
        view["_Date"] = pd.to_datetime(view["Date"], errors="coerce")
        view = (view.sort_values(by=["_Date","Cycle","SampleID"], ascending=[False, True, True])
                .drop(columns=["_Date"]).reset_index(drop=True))
        dbg("filter_meta_table_view", rows=len(view))
        return render.DataTable(view, selection_mode="none")

    @render.data_frame
    def table_overview():
        # Full metrics with metadata
        df = full_data.get()
        if df is None or df.empty:
            return render.DataTable(pd.DataFrame({"Info": ["No data loaded"]}), selection_mode="none")

        # Respect sidebar selection
        selected_raw = input.filter_samples()
        selected = _normalize_selected(selected_raw)
        if not selected:
            return render.DataTable(pd.DataFrame({"Info": ["No samples selected"]}), selection_mode="none")

        # Keep only selected samples
        view = df[df["SampleID"].astype(str).isin(selected)].copy()

        # Ensure the columns we want exist and are numeric where applicable
        cols_present = ["SampleID", "Date", "Cycle"]
        for col in METRICS_TO_PLOT:
            if col in view.columns:
                view[col] = pd.to_numeric(view[col], errors="coerce")
                cols_present.append(col)

        # Order: Date desc, Cycle asc, then SampleID
        view["_Date"] = pd.to_datetime(view["Date"], errors="coerce")
        view = (view.sort_values(by=["_Date", "Cycle", "SampleID"], ascending=[False, True, True])
                    .drop(columns=["_Date"], errors="ignore"))

        # Final column order (only existing ones)
        view = view[[c for c in cols_present if c in view.columns]].reset_index(drop=True)

        # Render as a sortable data table, no row selection needed
        return render.DataTable(view, selection_mode="none")

    # ---------- Plots ----------
    @render.ui
    def dashboard_content():
        if full_data.get() is None:
            return ui.p("Data is loading or no data was loaded.")
        tabs = [ui.nav_panel(m, output_widget(f"plot_{m}")) for m in METRICS_TO_PLOT]
        tabs.append(ui.nav_panel("Table Overview", ui.output_data_frame("table_overview")))
        tabs.append(ui.nav_panel("Read Length Distribution", output_widget("plot_read_length")))
        return ui.navset_card_tab(*tabs)

    # Metric plots (react to sidebar selection)
    for metric in METRICS_TO_PLOT:
        @output(id=f"plot_{metric}")
        @render_widget
        def _(metric_name: str = metric):
            df = full_data.get()
            selected_raw = input.filter_samples()
            selected = _normalize_selected(selected_raw)

            if df is None or not selected:
                return go.Figure(layout={"title": f"No samples selected for {metric_name}"})

            plot_df = df[df["SampleID"].isin(selected)].copy()
            if metric_name not in plot_df.columns:
                return go.Figure(layout={"title": f"{metric_name} not found in data"})

            hover_fmt = ":.2%" if THRESHOLDS.get(metric_name, {}).get("format", "").endswith("%") else ":.4f"

            fig = px.bar(
                plot_df.sort_values(metric_name, na_position="last"),
                x="SampleID",
                y=metric_name,
                color="Group",
                category_orders={"Group": group_order.get() or []},
                color_discrete_map=group_color_map.get() or {},
                hover_data={"SampleID": True, metric_name: hover_fmt, "Date": True, "Cycle": True, "Group": True},
                title=f"{metric_name} per Sample",
                labels={"SampleID": "Sample", metric_name: metric_name, "Group": "Date / Cycle"},
            )

            # Threshold lines (unchanged)
            if metric_name in THRESHOLDS:
                for line in THRESHOLDS[metric_name].get("lines", []):
                    fig.add_hline(
                        y=line["value"],
                        line_width=2,
                        line_dash="dash",
                        line_color=line["color"],
                        annotation_text=line.get("label"),
                        annotation_position="top left",
                    )
                axis_format = THRESHOLDS[metric_name].get("format")
                if axis_format:
                    fig.update_yaxes(tickformat=axis_format)

            fig.update_layout(margin=dict(l=40, r=40, t=50, b=80), bargap=0.25, legend_title_text="Date / Cycle")
            fig.update_xaxes(tickangle=-45)
            return fig

    @output(id="plot_read_length")
    @render_widget
    def _read_length_plot():
        df_all = readlen_data.get()
        selected_raw = input.filter_samples()
        selected = _normalize_selected(selected_raw)

        if not df_all or not selected:
            return go.Figure(layout={"title": "No samples selected for Read Length Distribution"})

        # Get metrics (and ensure Group exists)
        df_metrics = full_data.get()
        if df_metrics is None:
            df_metrics = pd.DataFrame(columns=["SampleID", "Date", "Cycle", "Group"])

        if "Group" not in df_metrics.columns:
            # Fallback: rebuild Group from metadata if needed
            _meta = sample_meta.get() or pd.DataFrame(columns=["SampleID", "Date", "Cycle"])
            df_metrics = df_metrics.merge(_meta[["SampleID", "Date", "Cycle"]], on="SampleID", how="left")
            df_metrics["Date"] = df_metrics["Date"].fillna("")
            df_metrics["Cycle"] = df_metrics["Cycle"].fillna("")
            df_metrics["Group"] = (df_metrics["Date"].replace("", "?") + " • " + df_metrics["Cycle"].replace("", "?"))

        # Map sample -> group
        sample_to_group = {}
        if not df_metrics.empty:
            sample_to_group = (
                df_metrics.drop_duplicates("SampleID")
                .set_index("SampleID")["Group"]
                .to_dict()
            )

        # Color map for groups
        cmap = group_color_map.get() or {}

        fig = go.Figure()
        any_curve = False
        seen_groups = set()

        for sid in selected:
            hist = df_all.get(sid)
            if hist is None or hist.empty:
                continue

            grp = sample_to_group.get(sid, "? • ?")
            color = cmap.get(grp)

            # One legend entry per group (clean, not per sample)
            if grp not in seen_groups:
                legend_kwargs = {"line": {"color": color}} if color else {}
                fig.add_trace(
                    go.Scatter(
                        x=[None], y=[None], mode="lines", name=grp,
                        hoverinfo="skip", showlegend=True, **legend_kwargs
                    )
                )
                seen_groups.add(grp)

            # Plot TOTAL
            if "UNPAIRED_TOTAL_LENGTH_COUNT" in hist.columns:
                trace_kwargs = {"line": {"color": color}} if color else {}
                fig.add_trace(
                    go.Scatter(
                        x=hist["READ_LENGTH"],
                        y=hist["UNPAIRED_TOTAL_LENGTH_COUNT"],
                        mode="lines",
                        name=f"{sid} — TOTAL",
                        showlegend=False,  # legend handled by group item above
                        hovertemplate="Sample: %{customdata[0]}<br>Series: TOTAL<br>Read length: %{x}<br>Count: %{y}<extra></extra>",
                        customdata=[[sid]] * len(hist),
                        **trace_kwargs,
                    )
                )
                any_curve = True

            # Plot ALIGNED
            if "UNPAIRED_ALIGNED_LENGTH_COUNT" in hist.columns:
                trace_kwargs = {"line": {"color": color}} if color else {}
                fig.add_trace(
                    go.Scatter(
                        x=hist["READ_LENGTH"],
                        y=hist["UNPAIRED_ALIGNED_LENGTH_COUNT"],
                        mode="lines",
                        name=f"{sid} — ALIGNED",
                        showlegend=False,  # legend handled by group item above
                        hovertemplate="Sample: %{customdata[0]}<br>Series: ALIGNED<br>Read length: %{x}<br>Count: %{y}<extra></extra>",
                        customdata=[[sid]] * len(hist),
                        **trace_kwargs,
                    )
                )
                any_curve = True

        if not any_curve:
            return go.Figure(layout={"title": "No histogram data available for selected samples"})

        fig.update_layout(
            title="Read Length Distribution (mmetrics alignment_summary_metrics)",
            xaxis_title="READ_LENGTH (bp)",
            yaxis_title="Count",
            margin=dict(l=50, r=30, t=60, b=50),
            legend_title="Date / Cycle",
        )
        return fig

# --- Run ---
app = App(app_ui, server)
