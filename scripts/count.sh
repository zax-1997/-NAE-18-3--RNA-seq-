#!/usr/bin/env bash
set -euo pipefail

BAM=~/AA/bam
GTF=~/AA/ref/Arabidopsis_thaliana.TAIR10.60.gtf
OUT=~/AA/counts
LOG=~/AA/logs
THREADS=8

mkdir -p "$OUT" "$LOG"

featureCounts \
    -T $THREADS \
    -p \
    --countReadPairs \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o "$OUT/gene_counts.txt" \
    "$BAM"/*.sorted.bam \
    2> "$LOG/featureCounts.log"

echo "定量完成"
