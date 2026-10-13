#!/usr/bin/env bash
set -euo pipefail

RAW=~/AA/fastq
CLEAN=~/AA/clean_data
FQC_RAW=~/AA/qc/fastqc_raw
FQC_CLEAN=~/AA/qc/fastqc_clean
FPR=~/AA/qc/fastp_report
THREADS=8

mkdir -p "$CLEAN" "$FQC_RAW" "$FQC_CLEAN" "$FPR"

echo "=== 1. 原始数据 FastQC ==="
fastqc -t $THREADS -o "$FQC_RAW" "$RAW"/*.fastq.gz

echo "=== 2. fastp 去接头 + 过滤 ==="
for r1 in "$RAW"/*_1.fastq.gz; do
    sample=$(basename "$r1" _1.fastq.gz)
    r2="$RAW/${sample}_2.fastq.gz"
    echo "[fastp] $sample"
    fastp \
        -i "$r1" -I "$r2" \
        -o "$CLEAN/${sample}_1.clean.fastq.gz" \
        -O "$CLEAN/${sample}_2.clean.fastq.gz" \
        -w $THREADS \
        -j "$FPR/${sample}.fastp.json" \
        -h "$FPR/${sample}.fastp.html" \
        -l 36 -q 20 -u 40 \
        2> "$FPR/${sample}.fastp.log"
done

echo "=== 3. 清洁数据 FastQC ==="
fastqc -t $THREADS -o "$FQC_CLEAN" "$CLEAN"/*.clean.fastq.gz

echo "=== 4. MultiQC 汇总 ==="
multiqc -f -o ~/AA/qc \
    "$FQC_RAW" "$FQC_CLEAN" "$FPR" \
    --title "AA Arabidopsis RNA-seq QC"

echo "QC 完成"
