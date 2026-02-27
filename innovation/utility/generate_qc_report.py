#!/usr/bin/env python3
"""
Generate a PDF QC report from a metrics TSV (e.g. WGS/HS-metrics output).

Usage:
    python utility/generate_qc_report.py --input metrics.tsv --project "FFPE Sandbox" --output qc_report.pdf
"""

import argparse
import os
import sys
from datetime import date

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


# ── Style constants ──────────────────────────────────────────────────────────
COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]
HEADER_BG = "#2c3e50"
HEADER_FG = "white"
ROW_ALT = "#f7f9fb"
BORDER_COLOR = "#bdc3c7"

DEPTH_LEVELS = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]


def load_data(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    return df


def compute_yield_gb(row):
    return (row["TotalReads"] * row["MeanReadLength"]) / 1e9


def page_header(ax, title, subtitle=""):
    """Draw a consistent page header."""
    ax.text(0.0, 0.85, title, fontsize=16, fontweight="bold", color=HEADER_BG,
            transform=ax.transAxes, va="top")
    if subtitle:
        ax.text(0.0, 0.55, subtitle, fontsize=9, color="#7f8c8d",
                transform=ax.transAxes, va="top")
    ax.axis("off")


def draw_table(ax, col_labels, row_data, col_widths=None):
    """Draw a styled table on the given axes."""
    ax.axis("off")
    n_cols = len(col_labels)
    n_rows = len(row_data)

    table = ax.table(
        cellText=row_data,
        colLabels=col_labels,
        loc="upper center",
        cellLoc="center",
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.4)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(BORDER_COLOR)
        cell.set_linewidth(0.5)
        if row == 0:
            cell.set_facecolor(HEADER_BG)
            cell.set_text_props(color=HEADER_FG, fontweight="bold")
        elif row % 2 == 0:
            cell.set_facecolor(ROW_ALT)
        else:
            cell.set_facecolor("white")


# ── Page builders ────────────────────────────────────────────────────────────

def build_summary_page(fig, df, project):
    """Page 1: summary header + AT/GC dropout mirrored bar plot."""
    gs = fig.add_gridspec(3, 1, height_ratios=[0.8, 0.4, 4], hspace=0.15,
                          left=0.10, right=0.92, top=0.94, bottom=0.06)

    # ── Title bar ────────────────────────────────────────────────────────
    ax_title = fig.add_subplot(gs[0])
    ax_title.text(0.0, 0.9, "QC Metrics Report", fontsize=20, fontweight="bold",
                  color=HEADER_BG, transform=ax_title.transAxes, va="top")
    ax_title.text(0.0, 0.25, f"Project: {project}   |   Samples: {len(df)}   |   "
                  f"Generated: {date.today().isoformat()}",
                  fontsize=9, color="#7f8c8d", transform=ax_title.transAxes, va="top")
    ax_title.axis("off")

    # ── Separator ────────────────────────────────────────────────────────
    ax_sep = fig.add_subplot(gs[1])
    ax_sep.axhline(0.5, color=BORDER_COLOR, linewidth=0.8)
    ax_sep.set_xlim(0, 1)
    ax_sep.set_ylim(0, 1)
    ax_sep.axis("off")

    # ── AT / GC dropout mirrored bar chart ───────────────────────────────
    ax = fig.add_subplot(gs[2])
    samples = df["SampleID"].values
    at_vals = pd.to_numeric(df["AT_DROPOUT"], errors="coerce").fillna(0).values
    gc_vals = pd.to_numeric(df["GC_DROPOUT"], errors="coerce").fillna(0).values

    y_pos = np.arange(len(samples))
    bar_h = 0.55

    ax.barh(y_pos, at_vals, height=bar_h, color="#e74c3c", label="AT Dropout", zorder=3)
    ax.barh(y_pos, -gc_vals, height=bar_h, color="#3498db", label="GC Dropout", zorder=3)

    # Zero line
    ax.axvline(0, color="#2c3e50", linewidth=0.8, zorder=4)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples, fontsize=8)
    ax.invert_yaxis()

    # Symmetric x limits
    max_val = max(at_vals.max(), gc_vals.max(), 1) * 1.15
    ax.set_xlim(-max_val, max_val)

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{abs(v):.1f}"))
    ax.set_xlabel("Dropout (%)", fontsize=9)
    ax.set_title("AT / GC Dropout", fontsize=12, fontweight="bold", pad=10, color=HEADER_BG)

    # Value labels on bars
    for i, (at, gc) in enumerate(zip(at_vals, gc_vals)):
        if at > 0:
            ax.text(at + max_val * 0.02, i, f"{at:.2f}", va="center", fontsize=7, color="#c0392b")
        if gc > 0:
            ax.text(-gc - max_val * 0.02, i, f"{gc:.2f}", va="center", fontsize=7,
                    ha="right", color="#2980b9")

    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(axis="x", linestyle="--", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def build_coverage_page(fig, df):
    """Page 2: depth-of-coverage line plot."""
    gs = fig.add_gridspec(2, 1, height_ratios=[0.6, 5], hspace=0.08,
                          left=0.10, right=0.92, top=0.94, bottom=0.08)

    ax_hdr = fig.add_subplot(gs[0])
    page_header(ax_hdr, "Depth-of-Coverage Curve",
                "Fraction of target bases covered at each minimum depth threshold")

    ax = fig.add_subplot(gs[1])

    depth_cols = [f"DepthCoverage_{d}" for d in DEPTH_LEVELS]
    present = [c for c in depth_cols if c in df.columns]
    present_levels = [int(c.split("_")[1]) for c in present]

    for idx, (_, row) in enumerate(df.iterrows()):
        vals = [float(row[c]) * 100 if pd.notna(row[c]) else np.nan for c in present]
        color = COLORS[idx % len(COLORS)]
        ax.plot(present_levels, vals, marker="o", markersize=4, linewidth=1.8,
                color=color, label=row["SampleID"], zorder=3)

    ax.set_xlabel("Minimum Depth (x)", fontsize=10)
    ax.set_ylabel("Bases Covered (%)", fontsize=10)
    ax.set_title("Coverage Distribution", fontsize=12, fontweight="bold", pad=10, color=HEADER_BG)
    ax.set_ylim(-2, 105)
    ax.set_xticks(present_levels)
    ax.legend(fontsize=7, loc="lower left", framealpha=0.9)
    ax.grid(True, linestyle="--", alpha=0.3, zorder=0)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def build_stats_page(fig, df):
    """Page 3: aggregate stats + per-sample table."""
    gs = fig.add_gridspec(3, 1, height_ratios=[0.6, 2.2, 4], hspace=0.25,
                          left=0.05, right=0.95, top=0.94, bottom=0.04)

    # ── Header ───────────────────────────────────────────────────────────
    ax_hdr = fig.add_subplot(gs[0])
    page_header(ax_hdr, "General Statistics Overview",
                "Cohort averages and per-sample breakdown")

    # ── Aggregate summary boxes ──────────────────────────────────────────
    ax_agg = fig.add_subplot(gs[1])
    ax_agg.axis("off")

    df["Yield_Gb"] = df.apply(compute_yield_gb, axis=1)

    metrics = [
        ("Median Coverage", "MedianCoverage", "x", ".1f"),
        ("MAD Coverage", "MadCoverage", "", ".1f"),
        ("Mean Read\nLength", "MeanReadLength", "bp", ".0f"),
        ("Phred score\nsubstitutions", "PhredSubstitutions", "", ".0f"),
        ("Yield", "Yield_Gb", "Gb", ".1f"),
        ("Duplicate %", "DuplicatePercent", "%", ".2f"),
        ("Low BaseQ %", "LowBaseQualPercent", "%", ".4f"),
        ("Low MapQ %", "LowMapQualPercent", "%", ".4f"),
    ]

    n = len(metrics)
    box_w = 1.0 / n
    for i, (label, col, unit, fmt) in enumerate(metrics):
        vals = pd.to_numeric(df[col], errors="coerce")
        avg = vals.mean()
        x_center = box_w * i + box_w / 2

        # Background box
        rect = plt.Rectangle((box_w * i + 0.005, 0.10), box_w - 0.01, 0.80,
                              transform=ax_agg.transAxes, facecolor="#eaf2f8",
                              edgecolor=BORDER_COLOR, linewidth=0.5, clip_on=False)
        ax_agg.add_patch(rect)

        ax_agg.text(x_center, 0.72, label, transform=ax_agg.transAxes,
                    ha="center", va="center", fontsize=7, color="#7f8c8d", fontweight="bold")
        display = f"{avg:{fmt}}"
        if unit:
            display += f" {unit}"
        ax_agg.text(x_center, 0.38, display, transform=ax_agg.transAxes,
                    ha="center", va="center", fontsize=12, color=HEADER_BG, fontweight="bold")

    # ── Per-sample table ─────────────────────────────────────────────────
    ax_tbl = fig.add_subplot(gs[2])

    col_labels = ["Sample", "Median Cov", "MAD Cov", "Mean RL",
                  "Yield (Gb)", "Phred score\n substitutions", "Dup %", "LowBQ %", "LowMQ %"]
    rows = []
    for _, row in df.iterrows():
        rows.append([
            str(row["SampleID"]),
            f"{float(row['MedianCoverage']):.1f}",
            f"{float(row['MadCoverage']):.1f}",
            f"{float(row['MeanReadLength']):.0f}",
            f"{compute_yield_gb(row):.1f}",
            f"{float(row['PhredSubstitutions']):.2f}",
            f"{float(row['DuplicatePercent']):.2f}",
            f"{float(row['LowBaseQualPercent']):.4f}",
            f"{float(row['LowMapQualPercent']):.4f}",
        ])
    # Append average row
    rows.append([
        "AVERAGE",
        f"{df['MedianCoverage'].astype(float).mean():.1f}",
        f"{df['MadCoverage'].astype(float).mean():.1f}",
        f"{df['MeanReadLength'].astype(float).mean():.0f}",
        f"{df['Yield_Gb'].mean():.1f}",
        f"{df['PhredSubstitutions'].astype(float).mean():.2f}",
        f"{df['DuplicatePercent'].astype(float).mean():.2f}",
        f"{df['LowBaseQualPercent'].astype(float).mean():.4f}",
        f"{df['LowMapQualPercent'].astype(float).mean():.4f}",
    ])

    draw_table(ax_tbl, col_labels, rows)

    # Bold the average row
    table = ax_tbl.tables[0]
    avg_row = len(rows)
    for col in range(len(col_labels)):
        cell = table[avg_row, col]
        cell.set_text_props(fontweight="bold")
        cell.set_facecolor("#d5e8d4")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate a PDF QC report from a metrics TSV."
    )
    parser.add_argument("--input", required=True, help="Input metrics TSV file")
    parser.add_argument("--project", default="", help="Project description for the report header")
    parser.add_argument("--output", default=None, help="Output PDF path (default: <input_stem>_qc_report.pdf)")
    args = parser.parse_args()

    if args.output is None:
        stem = os.path.splitext(args.input)[0]
        args.output = f"{stem}_qc_report.pdf"

    if not args.project:
        args.project = os.path.splitext(os.path.basename(args.input))[0]

    df = load_data(args.input)
    print(f"Loaded {len(df)} samples from {args.input}")

    with PdfPages(args.output) as pdf:
        # Page 1 — Summary + AT/GC Dropout
        fig1 = plt.figure(figsize=(8.27, 11.69))  # A4
        build_summary_page(fig1, df, args.project)
        pdf.savefig(fig1, facecolor="white")
        plt.close(fig1)

        # Page 2 — Coverage curve
        fig2 = plt.figure(figsize=(8.27, 11.69))
        build_coverage_page(fig2, df)
        pdf.savefig(fig2, facecolor="white")
        plt.close(fig2)

        # Page 3 — General stats + table
        fig3 = plt.figure(figsize=(8.27, 11.69))
        build_stats_page(fig3, df)
        pdf.savefig(fig3, facecolor="white")
        plt.close(fig3)

    print(f"QC report saved to {args.output}")


if __name__ == "__main__":
    main()
