#!/usr/bin/env python3
"""
Generate an ID SNP fingerprint image from AMBER SNP VCF files.

Input:  TSV with columns 'sampleId' and 'amberVcfPath' (GCS gs:// paths)
Output: PNG image showing genotypes and allele frequencies across all samples
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# ── Genotype codes (same as match_amber_sample.py) ──────────────────────────
GENOTYPE_MAP = {
    "DO_NOT_MATCH": 0,
    "HOM_REF": 1,
    "HET": 2,
    "HOM_ALT": 3,
}

GENOTYPE_COLORS = {
    0: "#f0f0f0",   # light grey — no call
    1: "#2ca02c",   # green     — HOM_REF
    2: "#ff7f0e",   # orange    — HET
    3: "#1f77b4",   # blue      — HOM_ALT
}

GENOTYPE_LABELS = {
    0: "No call",
    1: "Hom Ref",
    2: "Het",
    3: "Hom Alt",
}


def derive_genotype(ref_count, alt_count, total_coverage):
    """
    Derives a genotype code (0, 1, 2, 3) from AD/DP values.

    Uses lenient thresholds to maximise sensitivity for contamination
    detection: homozygous calls tolerate up to 10 % noise and the HET
    band spans everything in between (AF 0.10 – 0.90).  Contaminated
    samples will therefore show extra HET calls and shifted AF dots
    rather than being hidden as DO_NOT_MATCH.
    """
    if total_coverage < 10:
        return GENOTYPE_MAP["DO_NOT_MATCH"]

    allele_fraction = alt_count / total_coverage

    if allele_fraction <= 0.10:
        return GENOTYPE_MAP["HOM_REF"]

    if allele_fraction >= 0.90:
        return GENOTYPE_MAP["HOM_ALT"]

    return GENOTYPE_MAP["HET"]


def derive_noise_genotype(genotype):
    """
    Flip a genotype to represent the expected noise allele from contamination.

    HOM_REF → HOM_ALT (noise would carry alt alleles)
    HOM_ALT → HOM_REF (noise would carry ref alleles)
    HET / DO_NOT_MATCH → DO_NOT_MATCH (uninformative for contamination)
    """
    if genotype == GENOTYPE_MAP["HOM_REF"]:
        return GENOTYPE_MAP["HOM_ALT"]
    if genotype == GENOTYPE_MAP["HOM_ALT"]:
        return GENOTYPE_MAP["HOM_REF"]
    return GENOTYPE_MAP["DO_NOT_MATCH"]


def build_noise_profile(genotypes):
    """
    Build a noise genotype profile from a sample's majority genotype vector.

    Applies derive_noise_genotype element-wise.
    """
    return [derive_noise_genotype(gt) for gt in genotypes]


def stream_vcf_from_gcs(gcs_path):
    """
    Stream a VCF from GCS using 'gcloud storage cat' and decompress if gzipped.
    Returns the VCF content as a string.
    """
    result = subprocess.run(
        ["gcloud", "storage", "cat", gcs_path],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to stream {gcs_path}: {result.stderr.decode().strip()}"
        )
    raw = result.stdout
    # Decompress if gzipped
    if gcs_path.endswith(".gz"):
        raw = gzip.decompress(raw)
    return raw.decode("utf-8")


def parse_vcf_text(vcf_text, ordered_loci=None):
    """
    Parse an AMBER SNP VCF text and extract loci, genotypes and allele
    frequencies.  The loci order is taken from the VCF rows themselves.

    The AMBER VCF FORMAT is GT:AD:DP, sample column is e.g. 0/1:48,0:48.
    We ignore the GT field and derive genotype from AD/DP.

    If *ordered_loci* is provided the results are aligned to that order;
    otherwise the order is derived from the VCF rows and returned as well.

    Returns:
        loci:      list of (chrom, pos, ref, alt) tuples
        genotypes: list of genotype codes (0-3), one per locus
        afs:       list of allele frequencies (float, NaN for no-call)
    """
    # Parse every data row in file order
    vcf_data = {}   # (chrom, pos) -> {genotype, af}
    vcf_loci = []   # keeps insertion order

    for line in vcf_text.splitlines():
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 10:
            continue

        chrom = parts[0]
        pos = int(parts[1])
        ref = parts[3]
        alt = parts[4]
        fmt = parts[8]
        sample_data = parts[9]

        format_fields = fmt.split(":")
        sample_values = sample_data.split(":")
        fields = dict(zip(format_fields, sample_values))

        try:
            ad_parts = fields["AD"].split(",")
            ref_count = int(ad_parts[0])
            alt_count = int(ad_parts[1])
            dp = int(fields["DP"])
        except (KeyError, ValueError, IndexError):
            vcf_data[(chrom, pos)] = {
                'genotype': GENOTYPE_MAP["DO_NOT_MATCH"],
                'af': float('nan'),
            }
            vcf_loci.append((chrom, pos, ref, alt))
            continue

        genotype = derive_genotype(ref_count, alt_count, dp)
        af = alt_count / dp if dp > 0 else float('nan')

        vcf_data[(chrom, pos)] = {'genotype': genotype, 'af': af}
        vcf_loci.append((chrom, pos, ref, alt))

    # Use provided loci order, or the order from this VCF
    loci = ordered_loci if ordered_loci is not None else vcf_loci

    genotypes = []
    afs = []
    for chrom, pos, _ref, _alt in loci:
        entry = vcf_data.get((chrom, pos))
        if entry:
            genotypes.append(entry['genotype'])
            afs.append(entry['af'])
        else:
            genotypes.append(GENOTYPE_MAP["DO_NOT_MATCH"])
            afs.append(float('nan'))

    return loci, genotypes, afs


def compute_genotype_distance(gt_a, gt_b):
    """
    Compute distance between two genotype vectors.
    Only sites where both have a call (non-zero) are compared.
    Distance = 1 - (matching_sites / comparable_sites).
    """
    comparable = 0
    matches = 0
    for g1, g2 in zip(gt_a, gt_b):
        if g1 != 0 and g2 != 0:
            comparable += 1
            if g1 == g2:
                matches += 1
    if comparable == 0:
        return 1.0
    return 1.0 - (matches / comparable)


def compute_contamination_distance(noise_gt, source_gt):
    """
    Compute compatibility distance between a noise profile and a candidate
    source's majority genotype profile.

    At each informative site (neither side is DO_NOT_MATCH):
      - noise=HOM_ALT + source=HOM_ALT or HET → compatible (source has alt)
      - noise=HOM_REF + source=HOM_REF or HET → compatible (source has ref)
      - Opposite directions → incompatible

    Returns 1.0 - (compatible / informative_sites), or 1.0 if none.
    """
    informative = 0
    compatible = 0
    for n, s in zip(noise_gt, source_gt):
        if n == GENOTYPE_MAP["DO_NOT_MATCH"] or s == GENOTYPE_MAP["DO_NOT_MATCH"]:
            continue
        informative += 1
        if n == GENOTYPE_MAP["HOM_ALT"] and s in (GENOTYPE_MAP["HOM_ALT"], GENOTYPE_MAP["HET"]):
            compatible += 1
        elif n == GENOTYPE_MAP["HOM_REF"] and s in (GENOTYPE_MAP["HOM_REF"], GENOTYPE_MAP["HET"]):
            compatible += 1
    if informative == 0:
        return 1.0
    return 1.0 - (compatible / informative)


def generate_image(sample_ids, genotype_matrix, af_matrix, snp_labels, output_path):
    """Generate the ID SNP fingerprint image."""
    n_samples = len(sample_ids)
    n_snps = len(snp_labels)

    # ── Compute dendrogram ───────────────────────────────────────────────
    if n_samples >= 2:
        dist_matrix = np.zeros((n_samples, n_samples))
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                d = compute_genotype_distance(genotype_matrix[i], genotype_matrix[j])
                dist_matrix[i, j] = d
                dist_matrix[j, i] = d
        condensed = squareform(dist_matrix)
        Z = linkage(condensed, method='average')
        dendro_result = dendrogram(Z, no_plot=True)
        leaf_order = dendro_result['leaves']
    else:
        leaf_order = list(range(n_samples))
        Z = None

    # Reorder data by dendrogram leaf order
    sample_ids = [sample_ids[i] for i in leaf_order]
    genotype_matrix = [genotype_matrix[i] for i in leaf_order]
    af_matrix = [af_matrix[i] for i in leaf_order]

    # ── Figure layout ────────────────────────────────────────────────────
    row_height = 0.6
    fig_width = max(30, n_snps * 0.35)
    fig_height = max(6, n_samples * row_height + 4)

    if n_samples >= 2:
        fig = plt.figure(figsize=(fig_width + 3, fig_height))
        gs = GridSpec(1, 2, width_ratios=[20, 1], wspace=0.02, figure=fig)
        ax_main = fig.add_subplot(gs[0, 0])
        ax_dendro = fig.add_subplot(gs[0, 1])
    else:
        fig, ax_main = plt.subplots(figsize=(fig_width, fig_height))
        ax_dendro = None

    # ── Main heatmap + scatter ───────────────────────────────────────────
    # With invert_yaxis, row_idx is the TOP and row_idx+1 is the BOTTOM
    # of each row.  We want AF 0.0 at bottom and 1.0 at top, so place
    # dots at  row_idx + (1 - af).
    af_pad = 0.08  # inset so 0.0/1.0 labels don't overlap between rows

    for row_idx in range(n_samples):
        for col_idx in range(n_snps):
            gt = genotype_matrix[row_idx][col_idx]
            color = GENOTYPE_COLORS[gt]

            rect = plt.Rectangle(
                (col_idx, row_idx), 1, 1,
                facecolor=color, edgecolor='white', linewidth=0.5,
            )
            ax_main.add_patch(rect)

            af = af_matrix[row_idx][col_idx]
            if not np.isnan(af) and gt != 0:
                # Map af into the padded range within the row
                dot_y = row_idx + af_pad + (1 - af) * (1 - 2 * af_pad)
                ax_main.plot(
                    col_idx + 0.5, dot_y, 'o',
                    color='black', markersize=3, zorder=5,
                )

    ax_main.set_xlim(0, n_snps)
    ax_main.set_ylim(0, n_samples)
    ax_main.invert_yaxis()

    # Y-axis: sample IDs
    ax_main.set_yticks([i + 0.5 for i in range(n_samples)])
    ax_main.set_yticklabels(sample_ids, fontsize=8)

    # Mini AF scale ticks on the left (0.0 at bottom, 1.0 at top of each row)
    for row_idx in range(n_samples):
        top_y = row_idx + af_pad          # 1.0 position (top of row)
        bot_y = row_idx + 1 - af_pad      # 0.0 position (bottom of row)
        mid_y = row_idx + 0.5             # 0.5 position

        for y_pos, af_label in [(top_y, "1.0"), (mid_y, ""), (bot_y, "0.0")]:
            ax_main.plot(
                [-0.3, 0], [y_pos, y_pos],
                color='grey', linewidth=0.5, clip_on=False,
            )
            if af_label:
                ax_main.text(
                    -0.5, y_pos, af_label,
                    fontsize=5, ha='right', va='center', color='grey',
                )

    # X-axis: SNP labels
    ax_main.set_xticks([i + 0.5 for i in range(n_snps)])
    ax_main.set_xticklabels(snp_labels, rotation=90, fontsize=5, ha='center')
    ax_main.tick_params(axis='x', which='both', length=2)

    # Grid lines between rows
    for i in range(n_samples + 1):
        ax_main.axhline(y=i, color='grey', linewidth=0.3, zorder=1)

    ax_main.set_title("ID SNP Fingerprint", fontsize=14, fontweight='bold', pad=15)

    # ── Legend ───────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(facecolor=GENOTYPE_COLORS[gt], edgecolor='black', label=GENOTYPE_LABELS[gt])
        for gt in [1, 2, 3, 0]
    ]
    legend_patches.append(
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                   markersize=5, label='Allele Freq')
    )
    ax_main.legend(
        handles=legend_patches,
        loc='lower left', bbox_to_anchor=(0, 1.02),
        ncol=5, fontsize=8, frameon=True,
    )

    # ── Dendrogram ───────────────────────────────────────────────────────
    if ax_dendro is not None and Z is not None:
        dendro_plot = dendrogram(
            Z, orientation='right', ax=ax_dendro,
            leaf_rotation=0, leaf_font_size=0,
            color_threshold=0, above_threshold_color='#333333',
        )
        # Annotate each merge node with its similarity score (1 - distance)
        for i, (x_coords, y_coords) in enumerate(
            zip(dendro_plot['dcoord'], dendro_plot['icoord'])
        ):
            # The merge point is the rightmost x value; y is the midpoint
            x_merge = max(x_coords)
            y_merge = (y_coords[1] + y_coords[2]) / 2
            similarity = 1.0 - x_merge
            ax_dendro.text(
                x_merge, y_merge, f"{similarity:.2f}",
                fontsize=6, va='bottom', ha='center', color='#555555',
            )
        ax_dendro.set_xticks([])
        ax_dendro.set_yticks([])
        for spine in ax_dendro.spines.values():
            spine.set_visible(False)
        ax_dendro.set_title("Similarity", fontsize=9)

    # ── Save ─────────────────────────────────────────────────────────────
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Image saved to {output_path}")


def generate_contamination_image(
    target_id, noise_profile, other_ids, other_genotypes, other_afs,
    compatibility_scores, snp_labels, output_path,
):
    """
    Generate a contamination source analysis PNG for a single target sample.

    Top panel:  noise genotype profile of the target (colored rectangles only).
    Bottom panel: majority genotype profiles of all other samples with AF dots,
                  ranked by compatibility score (best match at top).
    """
    n_others = len(other_ids)
    n_snps = len(snp_labels)

    # ── Rank candidates by compatibility score (highest first) ───────────
    rank_order = sorted(
        range(n_others), key=lambda i: compatibility_scores[i], reverse=True,
    )
    other_ids = [other_ids[i] for i in rank_order]
    other_genotypes = [other_genotypes[i] for i in rank_order]
    other_afs = [other_afs[i] for i in rank_order]
    compatibility_scores = [compatibility_scores[i] for i in rank_order]

    # ── Figure layout ────────────────────────────────────────────────────
    row_height = 0.6
    fig_width = max(30, n_snps * 0.35)
    noise_height = row_height + 1.5  # space for title + legend + one row
    main_height = max(4, n_others * row_height + 2)

    fig = plt.figure(figsize=(fig_width, noise_height + main_height))
    gs = GridSpec(
        2, 1,
        height_ratios=[noise_height, main_height],
        hspace=0.3,
        figure=fig,
    )
    ax_noise = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[1, 0])

    af_pad = 0.08

    # ── Top panel: noise genotype profile (colored rectangles only) ──────
    for col_idx in range(n_snps):
        gt = noise_profile[col_idx]
        color = GENOTYPE_COLORS[gt]
        rect = plt.Rectangle(
            (col_idx, 0), 1, 1,
            facecolor=color, edgecolor='white', linewidth=0.5,
        )
        ax_noise.add_patch(rect)

    ax_noise.set_xlim(0, n_snps)
    ax_noise.set_ylim(0, 1)
    ax_noise.invert_yaxis()
    ax_noise.set_yticks([0.5])
    ax_noise.set_yticklabels([f"{target_id} (noise)"], fontsize=8)
    ax_noise.set_xticks([])
    ax_noise.axhline(y=0, color='grey', linewidth=0.3, zorder=1)
    ax_noise.axhline(y=1, color='grey', linewidth=0.3, zorder=1)

    ax_noise.set_title(
        f"Contamination Source Analysis: {target_id}",
        fontsize=14, fontweight='bold', pad=15,
    )

    # Legend on noise panel
    legend_patches = [
        mpatches.Patch(facecolor=GENOTYPE_COLORS[gt], edgecolor='black', label=GENOTYPE_LABELS[gt])
        for gt in [1, 2, 3, 0]
    ]
    legend_patches.append(
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                   markersize=5, label='Allele Freq')
    )
    ax_noise.legend(
        handles=legend_patches,
        loc='lower left', bbox_to_anchor=(0, 1.02),
        ncol=5, fontsize=8, frameon=True,
    )

    # ── Bottom panel: other samples with AF dots ─────────────────────────
    y_labels = [
        f"{sid} ({score:.2f})"
        for sid, score in zip(other_ids, compatibility_scores)
    ]

    for row_idx in range(n_others):
        for col_idx in range(n_snps):
            gt = other_genotypes[row_idx][col_idx]
            color = GENOTYPE_COLORS[gt]
            rect = plt.Rectangle(
                (col_idx, row_idx), 1, 1,
                facecolor=color, edgecolor='white', linewidth=0.5,
            )
            ax_main.add_patch(rect)

            af = other_afs[row_idx][col_idx]
            if not np.isnan(af) and gt != 0:
                dot_y = row_idx + af_pad + (1 - af) * (1 - 2 * af_pad)
                ax_main.plot(
                    col_idx + 0.5, dot_y, 'o',
                    color='black', markersize=3, zorder=5,
                )

    ax_main.set_xlim(0, n_snps)
    ax_main.set_ylim(0, n_others)
    ax_main.invert_yaxis()

    ax_main.set_yticks([i + 0.5 for i in range(n_others)])
    ax_main.set_yticklabels(y_labels, fontsize=8)

    # Mini AF scale ticks
    for row_idx in range(n_others):
        top_y = row_idx + af_pad
        bot_y = row_idx + 1 - af_pad
        mid_y = row_idx + 0.5
        for y_pos, af_label in [(top_y, "1.0"), (mid_y, ""), (bot_y, "0.0")]:
            ax_main.plot(
                [-0.3, 0], [y_pos, y_pos],
                color='grey', linewidth=0.5, clip_on=False,
            )
            if af_label:
                ax_main.text(
                    -0.5, y_pos, af_label,
                    fontsize=5, ha='right', va='center', color='grey',
                )

    # X-axis: SNP labels
    ax_main.set_xticks([i + 0.5 for i in range(n_snps)])
    ax_main.set_xticklabels(snp_labels, rotation=90, fontsize=5, ha='center')
    ax_main.tick_params(axis='x', which='both', length=2)

    # Grid lines
    for i in range(n_others + 1):
        ax_main.axhline(y=i, color='grey', linewidth=0.3, zorder=1)

    ax_main.set_title("Candidate Sources (compatibility score)", fontsize=10, pad=8)

    # ── Save ─────────────────────────────────────────────────────────────
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Contamination image saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate an ID SNP fingerprint image from AMBER SNP VCF files."
    )
    parser.add_argument(
        "--input", required=True,
        help="TSV file with columns 'sampleId' and 'amberVcfPath'"
    )
    parser.add_argument(
        "--output", default="idsnp_image.png",
        help="Output PNG path (default: idsnp_image.png)"
    )
    args = parser.parse_args()

    # Read input TSV
    samples = []
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            samples.append({
                'sampleId': row['sampleId'],
                'amberVcfPath': row['amberVcfPath'],
            })

    if not samples:
        print("Error: No samples found in input TSV.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(samples)} samples in input TSV")

    # Stream and parse each VCF directly from GCS.
    # Loci order is derived from the first VCF; subsequent samples are aligned to it.
    sample_ids = []
    all_genotypes = []
    all_afs = []
    ordered_loci = None

    for i, sample in enumerate(samples):
        sid = sample['sampleId']
        gcs_path = sample['amberVcfPath']
        print(f"  [{i+1}/{len(samples)}] Streaming {sid}...")

        vcf_text = stream_vcf_from_gcs(gcs_path)
        loci, genotypes, afs = parse_vcf_text(vcf_text, ordered_loci)

        if ordered_loci is None:
            ordered_loci = loci
            print(f"  Derived {len(ordered_loci)} SNP loci from first sample")

        sample_ids.append(sid)
        all_genotypes.append(genotypes)
        all_afs.append(afs)

    # Build SNP labels from the loci
    snp_labels = [
        f"{chrom}:{pos} {ref}>{alt}"
        for chrom, pos, ref, alt in ordered_loci
    ]

    # Generate the image
    generate_image(sample_ids, all_genotypes, all_afs, snp_labels, args.output)

    # ── Contamination source analysis ────────────────────────────────────
    if len(sample_ids) >= 2:
        output_stem, _ext = os.path.splitext(args.output)
        print(f"\nRunning contamination source analysis for {len(sample_ids)} samples...")

        for idx in range(len(sample_ids)):
            target_id = sample_ids[idx]
            noise_profile = build_noise_profile(all_genotypes[idx])

            other_ids = [sample_ids[j] for j in range(len(sample_ids)) if j != idx]
            other_genotypes = [all_genotypes[j] for j in range(len(sample_ids)) if j != idx]
            other_afs = [all_afs[j] for j in range(len(sample_ids)) if j != idx]

            compatibility_scores = [
                1.0 - compute_contamination_distance(noise_profile, gt)
                for gt in other_genotypes
            ]

            cont_output = f"{output_stem}_contamination_{target_id}.png"
            generate_contamination_image(
                target_id, noise_profile,
                other_ids, other_genotypes, other_afs,
                compatibility_scores, snp_labels, cont_output,
            )
    else:
        print("\nSkipping contamination analysis (need >= 2 samples)")


if __name__ == "__main__":
    main()
