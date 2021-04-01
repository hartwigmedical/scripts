#!/usr/bin/env bash

pilot_dir=/data/common/dbs/serve/pilot_output
prod_dir=/data/common/dbs/serve/prod_output

echo "[INFO] Verfiy actionable charateristics"
diff ${prod_dir}/ActionableSignatures.37.tsv ${pilot_dir}/ActionableCharacteristics.37.tsv

echo "[INFO] Verfiy actionable fusions"
diff ${prod_dir}/ActionableFusions.37.tsv ${pilot_dir}/ActionableFusions.37.tsv

echo "[INFO] Verfiy actionable genes"
diff ${prod_dir}/ActionableGenes.37.tsv ${pilot_dir}/ActionableGenes.37.tsv

echo "[INFO] Verfiy actionable hotspots"
diff ${prod_dir}/ActionableHotspots.37.tsv ${pilot_dir}/ActionableHotspots.37.tsv

echo "[INFO] Verfiy actionable ranges"
diff ${prod_dir}/ActionableRanges.37.tsv ${pilot_dir}/ActionableRanges.37.tsv

echo "[INFO] Verfiy known codons"
diff ${prod_dir}/KnownCodons.SERVE.37.tsv ${pilot_dir}/KnownCodons.SERVE.37.tsv

echo "[INFO] Verfiy known copy numbers"
diff ${prod_dir}/KnownCopyNumbers.SERVE.37.tsv ${pilot_dir}/KnownCopyNumbers.SERVE.37.tsv

echo "[INFO] Verfiy known exons"
diff ${prod_dir}/KnownExons.SERVE.37.tsv ${pilot_dir}/KnownExons.SERVE.37.tsv

echo "[INFO] Verfiy known fusion pairs"
diff ${prod_dir}/KnownFusionPairs.SERVE.37.tsv ${pilot_dir}/KnownFusionPairs.SERVE.37.tsv

echo "[INFO] Verfiy known hotspots"
zdiff ${prod_dir}/KnownHotspots.SERVE.37.vcf.gz ${pilot_dir}/KnownHotspots.SERVE.37.vcf.gz