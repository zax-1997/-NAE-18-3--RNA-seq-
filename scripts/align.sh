#!/usr/bin/env bash
set -euo pipefail

INDEX=~/AA/ref/ath_tair10
CLEAN=~/AA/clean_data
BAM=~/AA/bam
LOG=~/AA/logs
THREADS=8

mkdir -p "$BAM" "$LOG"

for r1 in "$CLEAN"/*_1.clean.fastq.gz; do
    sample=$(basename "$r1" _1.clean.fastq.gz)
    r2="$CLEAN/${sample}_2.clean.fastq.gz"
    echo "=== 比对 $sample ==="

    hisat2 -x "$INDEX" \
        -1 "$r1" -2 "$r2" \
        -p $THREADS \
        --dta \
        2> "$LOG/${sample}.hisat2.log" \
    | samtools sort -@ $THREADS -o "$BAM/${sample}.sorted.bam" -

    samtools index "$BAM/${sample}.sorted.bam"
    samtools flagstat "$BAM/${sample}.sorted.bam" > "$LOG/${sample}.flagstat.txt"

    echo "=== $sample 完成 ==="
done

echo "比对全部完成"
