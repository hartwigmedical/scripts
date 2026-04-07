#!/usr/bin/env python3
"""
Compare QSEE multisample output against cohort percentile distributions.

Outputs:
  - TSV: Summary table metrics with percentile ranks
  - PDF: Distribution plots with cohort percentile bands

Usage:
    python utility/qsee_percentile_report.py \
      --sample-data multisample.qsee.vis.data.tsv.gz \
      --percentiles qsee.cohort.percentiles.tsv.gz \
      --output-prefix qsee_report
"""

import argparse
import sys
from bisect import bisect_left

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


# ── Style constants (consistent with generate_qc_report.py) ─────────────────
HEADER_BG = "#2c3e50"
HEADER_FG = "white"
BAND_LIGHT = "#aec7e8"
BAND_DARK = "#1f77b4"
AVG_COLOR = "#d62728"
MEDIAN_COLOR = "#2c3e50"
BORDER_COLOR = "#bdc3c7"
ROW_ALT = "#f7f9fb"

SKIP_FEATURE_TYPES = {"BQR_PER_SNV96_CONTEXT", "MISSED_VARIANT_LIKELIHOOD"}

FEATURE_X_FIELD = {
    "COVERAGE_DISTRIBUTION": "ReadDepth",
    "GC_BIAS": "GCBucket",
    "DUPLICATE_FREQ": "ReadCount",
    "MS_INDEL_ERROR_RATES": "RefNumUnits",
    "MS_INDEL_ERROR_BIAS": "RefNumUnits",
    "BQR_PER_ORIG_QUAL": "StandardMutation",
    "DISCORDANT_FRAG_FREQ": "DiscordantFragType",
}

PCT_COLS = [f"Pct_{i}" for i in range(101)]


def load_data(path):
    return pd.read_csv(path, sep="\t", compression="infer")


def compute_percentile_rank(value, pct_values):
    idx = bisect_left(pct_values, value)
    return min(idx, 100)


def parse_feature_name(feature_name, x_field):
    """Parse FeatureName into (subgroup, x_value).

    Splits semicolon-separated key:value pairs. The pair matching x_field
    becomes the x-axis value; remaining pairs form the subgroup label.
    """
    parts = feature_name.split(";")
    kvs = {}
    for part in parts:
        key, _, val = part.partition(":")
        kvs[key] = val

    x_val = kvs.pop(x_field, feature_name)

    try:
        x_num = float(x_val)
        x_val = int(x_num) if x_num == int(x_num) else x_num
    except (ValueError, TypeError):
        pass

    subgroup = " / ".join(kvs.values()) if kvs else ""
    return subgroup, x_val


AXIS_LABELS = {
    "COVERAGE_DISTRIBUTION": ("Read depth", "Proportion of bases"),
    "GC_BIAS": ("GC content (%)", "Normalized coverage"),
    "DUPLICATE_FREQ": ("Duplicate read count", "Frequency"),
    "MS_INDEL_ERROR_RATES": ("Repeat unit count", "Error rate (%)"),
    "MS_INDEL_ERROR_BIAS": ("Repeat unit count", "Error bias"),
    "BQR_PER_ORIG_QUAL": ("Mutation type", "BQ adjustment"),
    "DISCORDANT_FRAG_FREQ": ("SV type", "Frequency"),
    "FRAGMENT_LENGTH": ("Fragment length (bp)", "Frequency"),
}


def plot_distribution(ax, x_values, avg_values, pct_data, title, is_numeric,
                      xlabel=None, ylabel=None):
    """Plot a single distribution with percentile bands.

    For numeric x: line plot with fill_between bands.
    For categorical x: markers with shaded rectangles.
    """
    if is_numeric:
        order = sorted(range(len(x_values)), key=lambda i: x_values[i])
        xs = [x_values[i] for i in order]
        avgs = [avg_values[i] for i in order]

        p10 = [pct_data[x_values[i]][10] for i in order]
        p25 = [pct_data[x_values[i]][25] for i in order]
        p50 = [pct_data[x_values[i]][50] for i in order]
        p75 = [pct_data[x_values[i]][75] for i in order]
        p90 = [pct_data[x_values[i]][90] for i in order]

        ax.fill_between(xs, p10, p90, alpha=0.15, color=BAND_DARK, label="10th-90th pct")
        ax.fill_between(xs, p25, p75, alpha=0.25, color=BAND_DARK, label="25th-75th pct")
        ax.plot(xs, p50, "--", color=MEDIAN_COLOR, linewidth=1, alpha=0.7, label="Cohort median")
        ax.plot(xs, avgs, "-", color=AVG_COLOR, linewidth=1.5, label="Input average")
    else:
        x_pos = list(range(len(x_values)))

        for i, x in enumerate(x_values):
            p = pct_data[x]
            ax.fill_between([i - 0.3, i + 0.3], [p[10]] * 2, [p[90]] * 2,
                            alpha=0.15, color=BAND_DARK)
            ax.fill_between([i - 0.3, i + 0.3], [p[25]] * 2, [p[75]] * 2,
                            alpha=0.25, color=BAND_DARK)

        p50 = [pct_data[x][50] for x in x_values]
        ax.plot(x_pos, p50, "s", color=MEDIAN_COLOR, markersize=5, alpha=0.7, label="Cohort median")
        ax.plot(x_pos, avg_values, "o", color=AVG_COLOR, markersize=6, label="Input average", zorder=5)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_values, rotation=45, ha="right", fontsize=7)

    ax.set_title(title, fontsize=9, fontweight="bold", color=HEADER_BG)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=8)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.3)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def collect_subplots(tool_data, sample_type, source_tool, pct_df):
    """Collect all subplot specs for a given SourceTool."""
    subplots = []

    for feature_type in sorted(tool_data["FeatureType"].unique()):
        x_field = FEATURE_X_FIELD.get(feature_type)
        if x_field is None:
            continue

        ft_data = tool_data[tool_data["FeatureType"] == feature_type]

        # Parse feature names into subgroups
        parsed = {}
        for _, row in ft_data.iterrows():
            subgroup, x_val = parse_feature_name(row["FeatureName"], x_field)
            parsed.setdefault(subgroup, []).append((x_val, row["FeatureName"], row["FeatureValue"]))

        for subgroup in sorted(parsed.keys()):
            items = parsed[subgroup]

            # Average across samples for each x_value
            x_val_groups = {}
            x_val_fname = {}
            for x_val, fname, fval in items:
                x_val_groups.setdefault(x_val, []).append(fval)
                x_val_fname[x_val] = fname

            x_values = list(x_val_groups.keys())
            avg_values = [np.mean(vs) for vs in x_val_groups.values()]
            is_numeric = all(isinstance(x, (int, float)) for x in x_values)

            # Look up percentile data
            pct_data = {}
            for x_val in x_values:
                pct_row = pct_df[
                    (pct_df["SampleType"] == sample_type)
                    & (pct_df["SourceTool"] == source_tool)
                    & (pct_df["FeatureType"] == feature_type)
                    & (pct_df["FeatureName"] == x_val_fname[x_val])
                ]
                if not pct_row.empty:
                    pct_data[x_val] = pct_row[PCT_COLS].values[0].astype(float)

            if not pct_data:
                continue

            valid_x = [x for x in x_values if x in pct_data]
            valid_avg = [avg_values[x_values.index(x)] for x in valid_x]

            title_parts = [feature_type]
            if subgroup:
                title_parts.append(subgroup)
            title = " \u2014 ".join(title_parts)

            xlabel, ylabel = AXIS_LABELS.get(feature_type, (None, None))
            subplots.append((title, valid_x, valid_avg, pct_data, is_numeric, xlabel, ylabel))

    return subplots


def _parse_plot_metadata(meta_str):
    """Extract DisplayName, NumberFormat, and FeatureGroup from PlotMetadata."""
    display_name, num_format, feature_group = None, None, None
    if pd.isna(meta_str) or not meta_str:
        return display_name, num_format, feature_group
    for part in str(meta_str).split(";"):
        key, _, val = part.partition(":")
        if key == "DisplayName":
            display_name = val
        elif key == "NumberFormat":
            num_format = val
        elif key == "FeatureGroup":
            feature_group = val
    return display_name, num_format, feature_group


def _format_value(value, num_format):
    """Format a numeric value based on its NumberFormat."""
    if isinstance(value, str):
        return value
    if num_format == "PERCENT":
        return f"{value * 100:.2f}%"
    elif num_format == "LOG10":
        return f"{value:.2f}"
    else:
        return f"{value:.2f}"



# ── Summary page: grouped sections with shared axes ─────────────────────────

# Which FeatureNames belong in which section, and the axis config for each.
# Order within each section defines the row order on the page.
SUMMARY_SECTIONS = [
    {
        "title": "Proportions (0\u2013100%)",
        "xlim": (0, 100),
        "transform": lambda v: v * 100,   # 0-1 fraction → 0-100%
        "fmt": lambda v: f"{v:.1f}%",
        "tick_values": [0, 25, 50, 75, 100],
        "features": [
            "MAPPED_PROPORTION", "COVERAGE_ABOVE_10", "COVERAGE_ABOVE_20",
            "COVERAGE_ABOVE_30", "COVERAGE_ABOVE_60", "COVERAGE_ABOVE_100",
            "COVERAGE_ABOVE_250", "LOW_BASE_QUAL", "LOW_MAP_QUAL",
            "DUPLICATE_READS", "DUAL_STRAND_READS", "PURITY", "LOH_PERCENT",
            "TINC", "CONTAMINATION",
        ],
    },
    {
        "title": "Counts (0\u20131000)",
        "xlim": (0, 1000),
        "transform": lambda v: v,
        "fmt": lambda v: f"{v:.0f}",
        "tick_values": [0, 250, 500, 750, 1000],
        "features": ["DELETED_GENES", "UNSUPPORTED_CN_SEGMENTS"],
    },
    {
        "title": "Mutational burden (log scale)",
        "xlim": (0.01, 1000),
        "log": True,
        "transform": lambda v: v,
        "fmt": lambda v: f"{v:.2f}",
        "tick_values": [0.01, 0.1, 1, 10, 100, 1000],
        "features": ["TMB_SMALL_VARIANTS", "TMB_MS_INDELS", "TMB_STRUCTURAL_VARIANTS"],
    },
    {
        "title": "Mean coverage",
        "xlim": (0, 180),
        "transform": lambda v: v,
        "fmt": lambda v: f"{v:.1f}",
        "tick_values": [0, 45, 90, 135, 180],
        "features": ["MEAN_COVERAGE"],
    },
    {
        "title": "Ploidy",
        "xlim": (0, 5),
        "transform": lambda v: v,
        "fmt": lambda v: f"{v:.2f}",
        "tick_values": [0, 1, 2, 3, 4, 5],
        "features": ["PLOIDY"],
    },
]


def _draw_section(ax, section, section_metrics):
    """Draw one section of the summary page on the given axes."""
    xlim = section["xlim"]
    is_log = section.get("log", False)
    tr = section["transform"]
    fmt = section["fmt"]

    n = len(section_metrics)
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()

    if is_log:
        ax.set_xscale("log")
    ax.set_xlim(xlim)
    ax.set_xticks(section["tick_values"])
    ax.set_xticklabels([str(t) for t in section["tick_values"]], fontsize=7)
    ax.xaxis.set_tick_params(length=3)

    ax.set_yticks(range(n))
    ax.set_yticklabels([m["display_name"] for m in section_metrics], fontsize=8)

    ax.set_title(section["title"], fontsize=9, fontweight="bold", color=HEADER_BG,
                 loc="left")
    ax.grid(axis="x", linestyle="--", alpha=0.3)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    bar_h = 0.55
    for i, m in enumerate(section_metrics):
        # Apply transform so raw values map to the section's x-axis
        avg = tr(m["avg"])
        p10, p25, p50, p75, p90 = (tr(m["p10"]), tr(m["p25"]), tr(m["p50"]),
                                    tr(m["p75"]), tr(m["p90"]))

        # Bands
        ax.barh(i, p90 - p10, left=p10, height=bar_h,
                color=BAND_DARK, alpha=0.15, edgecolor="none")
        ax.barh(i, p75 - p25, left=p25, height=bar_h,
                color=BAND_DARK, alpha=0.25, edgecolor="none")

        # Median tick
        ax.plot([p50, p50], [i - bar_h / 2, i + bar_h / 2],
                color=MEDIAN_COLOR, linewidth=1.5, alpha=0.7)

        # Input average dot
        ax.plot(avg, i, "o", color=AVG_COLOR, markersize=7, zorder=5)

        # Value annotations on the right (outside plot area)
        avg_str = fmt(avg)
        med_str = fmt(p50)
        ax.annotate(avg_str, xy=(1.02, i), xycoords=("axes fraction", "data"),
                    fontsize=7, color=AVG_COLOR, fontweight="bold", va="center")
        ax.annotate(f"({med_str})", xy=(1.18, i), xycoords=("axes fraction", "data"),
                    fontsize=7, color="#7f8c8d", va="center")

        # Alternating row background
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=ROW_ALT, zorder=0)


def build_summary_page(pdf, metrics, sample_type, n_samples):
    """Add a visual dot-range summary page with grouped shared axes."""
    # Index metrics by FeatureName for lookup
    metrics_by_name = {m["feature_name"]: m for m in metrics}

    # Build sections with available data
    active_sections = []
    for section in SUMMARY_SECTIONS:
        section_metrics = [metrics_by_name[f] for f in section["features"]
                           if f in metrics_by_name]
        if section_metrics:
            active_sections.append((section, section_metrics))

    if not active_sections:
        return

    # Calculate height ratios: proportional to number of metrics + space for title
    height_ratios = [0.4]  # title
    for section, section_metrics in active_sections:
        height_ratios.append(len(section_metrics) + 1.2)  # +1.2 for section title & padding
    height_ratios.append(0.3)  # legend at bottom

    fig, all_axes = plt.subplots(len(height_ratios), 1,
                                  figsize=(11.69, 8.27),
                                  gridspec_kw={"height_ratios": height_ratios})

    # Title
    ax_title = all_axes[0]
    ax_title.text(0.0, 0.8, "QSEE Percentile Report", fontsize=18, fontweight="bold",
                  color=HEADER_BG, transform=ax_title.transAxes, va="top")
    ax_title.text(0.5, 0.8, f"Sample type: {sample_type}   |   Samples: {n_samples}",
                  fontsize=10, color="#7f8c8d", transform=ax_title.transAxes, va="top")
    ax_title.axis("off")

    # Sections
    for idx, (section, section_metrics) in enumerate(active_sections):
        ax = all_axes[idx + 1]
        _draw_section(ax, section, section_metrics)

    # Legend at bottom
    ax_legend = all_axes[-1]
    ax_legend.axis("off")
    ax_legend.barh([], 0, color=BAND_DARK, alpha=0.15, label="10th\u201390th pct")
    ax_legend.barh([], 0, color=BAND_DARK, alpha=0.25, label="25th\u201375th pct")
    ax_legend.plot([], [], "|", color=MEDIAN_COLOR, markersize=10, label="Cohort median")
    ax_legend.plot([], [], "o", color=AVG_COLOR, markersize=7, label="Input average")
    ax_legend.legend(fontsize=8, loc="center", framealpha=0.9, ncol=4)

    fig.tight_layout(rect=[0.01, 0, 1, 1], h_pad=0.8)
    pdf.savefig(fig, facecolor="white")
    plt.close(fig)


def generate_pdf(sample_df, pct_df, output_path):
    dist_df = sample_df[
        (sample_df["FeatureType"] != "SUMMARY_TABLE")
        & (~sample_df["FeatureType"].isin(SKIP_FEATURE_TYPES))
    ].copy()
    dist_df["FeatureValue"] = pd.to_numeric(dist_df["FeatureValue"], errors="coerce")

    # Build summary page data
    summary = sample_df[sample_df["FeatureType"] == "SUMMARY_TABLE"].copy()
    summary["FeatureValue"] = pd.to_numeric(summary["FeatureValue"], errors="coerce")
    n_samples = sample_df["SampleId"].nunique()

    with PdfPages(output_path) as pdf:
        # Summary dot-range page(s) first
        for sample_type in sorted(summary["SampleType"].unique()):
            st_summary = summary[summary["SampleType"] == sample_type]
            metrics = []
            for feat_name in st_summary["FeatureName"].unique():
                feat_data = st_summary[st_summary["FeatureName"] == feat_name]
                avg_val = feat_data["FeatureValue"].mean()

                meta = feat_data["PlotMetadata"].iloc[0]
                display_name, num_format, feature_group = _parse_plot_metadata(meta)
                if not display_name:
                    display_name = feat_name

                pct_row = pct_df[
                    (pct_df["SampleType"] == sample_type)
                    & (pct_df["FeatureType"] == "SUMMARY_TABLE")
                    & (pct_df["FeatureName"] == feat_name)
                ]

                if pct_row.empty:
                    continue

                pct_values = pct_row[PCT_COLS].values[0].astype(float)
                metrics.append({
                    "feature_name": feat_name,
                    "display_name": display_name,
                    "num_format": num_format,
                    "feature_group": feature_group,
                    "avg": avg_val,
                    "p10": pct_values[10],
                    "p25": pct_values[25],
                    "p50": pct_values[50],
                    "p75": pct_values[75],
                    "p90": pct_values[90],
                })

            if metrics:
                build_summary_page(pdf, metrics, sample_type, n_samples)

        # Distribution plots
        OVERVIEW_TOOLS = {"BAM_METRICS", "COBALT", "ESVEE"}

        for sample_type in sorted(dist_df["SampleType"].unique()):
            st_data = dist_df[dist_df["SampleType"] == sample_type]

            # ── Overview page: coverage, GC bias, discordant frag, fragment length
            overview_subplots = []
            for tool in ["BAM_METRICS", "COBALT", "ESVEE"]:
                tool_data = st_data[st_data["SourceTool"] == tool]
                if not tool_data.empty:
                    overview_subplots.extend(
                        collect_subplots(tool_data, sample_type, tool, pct_df))

            fig, axes = plt.subplots(2, 2, figsize=(11.69, 8.27))
            fig.suptitle(f"{sample_type} \u2014 Overview", fontsize=14,
                         fontweight="bold", color=HEADER_BG, y=0.98)

            for i, (title, x_vals, avg_vals, pct_d, is_num, xl, yl) in enumerate(overview_subplots[:4]):
                r, c = divmod(i, 2)
                plot_distribution(axes[r, c], x_vals, avg_vals, pct_d, title,
                                  is_num, xl, yl)

            # Fill remaining slots (fragment length placeholder if < 4 plots)
            for i in range(len(overview_subplots[:4]), 4):
                r, c = divmod(i, 2)
                ax = axes[r, c]
                ax.set_title("FRAGMENT_LENGTH", fontsize=9, fontweight="bold",
                             color=HEADER_BG)
                xl, yl = AXIS_LABELS.get("FRAGMENT_LENGTH", (None, None))
                if xl:
                    ax.set_xlabel(xl, fontsize=8)
                if yl:
                    ax.set_ylabel(yl, fontsize=8)
                ax.text(0.5, 0.5, "No data available", fontsize=11, color="#bdc3c7",
                        ha="center", va="center", transform=ax.transAxes)
                ax.grid(True, linestyle="--", alpha=0.3)
                for spine in ["top", "right"]:
                    ax.spines[spine].set_visible(False)

            axes[0, 0].legend(fontsize=6, loc="best")
            fig.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)

            # ── Remaining tools (REDUX etc.) on their own pages
            for source_tool in sorted(st_data["SourceTool"].unique()):
                if source_tool in OVERVIEW_TOOLS:
                    continue

                tool_data = st_data[st_data["SourceTool"] == source_tool]
                subplots = collect_subplots(tool_data, sample_type, source_tool, pct_df)

                if not subplots:
                    continue

                plots_per_page = 6
                for page_start in range(0, len(subplots), plots_per_page):
                    page_plots = subplots[page_start:page_start + plots_per_page]
                    n = len(page_plots)
                    nrows = min(3, (n + 1) // 2)
                    ncols = 2 if n > 1 else 1

                    fig, axes_grid = plt.subplots(nrows, ncols, figsize=(11.69, 8.27))

                    page_title = f"{sample_type} \u2014 {source_tool}"
                    fig.suptitle(page_title, fontsize=14, fontweight="bold",
                                 color=HEADER_BG, y=0.98)

                    if nrows == 1 and ncols == 1:
                        axes_grid = np.array([[axes_grid]])
                    elif nrows == 1:
                        axes_grid = axes_grid.reshape(1, -1)
                    elif ncols == 1:
                        axes_grid = axes_grid.reshape(-1, 1)

                    for i, (title, x_vals, avg_vals, pct_d, is_num, xl, yl) in enumerate(page_plots):
                        row, col = divmod(i, ncols)
                        plot_distribution(axes_grid[row, col], x_vals, avg_vals,
                                          pct_d, title, is_num, xl, yl)

                    axes_grid[0, 0].legend(fontsize=6, loc="best")

                    for i in range(n, nrows * ncols):
                        row, col = divmod(i, ncols)
                        axes_grid[row, col].set_visible(False)

                    fig.tight_layout(rect=[0, 0, 1, 0.95])
                    pdf.savefig(fig, facecolor="white")
                    plt.close(fig)

    print(f"PDF report written to {output_path}")


def generate_summary_tsv(sample_df, pct_df, output_path):
    summary = sample_df[sample_df["FeatureType"] == "SUMMARY_TABLE"].copy()
    summary["FeatureValue"] = pd.to_numeric(summary["FeatureValue"], errors="coerce")

    rows = []
    for sample_type in sorted(summary["SampleType"].unique()):
        st_data = summary[summary["SampleType"] == sample_type]
        for feat_name in st_data["FeatureName"].unique():
            feat_data = st_data[st_data["FeatureName"] == feat_name]
            avg_val = feat_data["FeatureValue"].mean()

            pct_row = pct_df[
                (pct_df["SampleType"] == sample_type)
                & (pct_df["FeatureType"] == "SUMMARY_TABLE")
                & (pct_df["FeatureName"] == feat_name)
            ]

            if pct_row.empty:
                median = np.nan
                pct_rank = "N/A"
            else:
                pct_values = pct_row[PCT_COLS].values[0].astype(float)
                median = pct_values[50]
                pct_rank = compute_percentile_rank(avg_val, list(pct_values))

            rows.append({
                "SampleType": sample_type,
                "FeatureName": feat_name,
                "InputSetAvg": round(avg_val, 6),
                "CohortMedian": round(median, 6) if not np.isnan(median) else "N/A",
                "PercentileRank": pct_rank,
            })

    result = pd.DataFrame(rows)
    result.to_csv(output_path, sep="\t", index=False)
    print(f"Summary TSV written to {output_path} ({len(rows)} features)")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Compare QSEE sample data against cohort percentiles."
    )
    parser.add_argument("--sample-data", required=True,
                        help="Multisample qsee vis data TSV (.tsv or .tsv.gz)")
    parser.add_argument("--percentiles", required=True,
                        help="Cohort percentiles TSV (.tsv or .tsv.gz)")
    parser.add_argument("--output-prefix", default="qsee_percentile_report",
                        help="Output prefix (produces <prefix>.tsv and <prefix>.pdf)")
    args = parser.parse_args()

    print("Loading data...")
    sample_df = load_data(args.sample_data)
    pct_df = load_data(args.percentiles)

    sample_types = sorted(sample_df["SampleType"].unique())
    n_samples = sample_df["SampleId"].nunique()
    print(f"Found {n_samples} samples, types: {', '.join(sample_types)}")

    tsv_path = f"{args.output_prefix}.tsv"
    generate_summary_tsv(sample_df, pct_df, tsv_path)

    pdf_path = f"{args.output_prefix}.pdf"
    generate_pdf(sample_df, pct_df, pdf_path)

    print("Done!")


if __name__ == "__main__":
    main()
