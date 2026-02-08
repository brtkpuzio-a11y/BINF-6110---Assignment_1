#!/usr/bin/env bash
set -euo pipefail

ACC=SRR32410565
REF="GCF_002507875.2_ASM250787v2_genomic.fna"
BAM="mapping_${ACC}/${ACC}.ref.sorted.bam"
OUT="bcftools_${ACC}.vcf.gz"

module load StdEnv/2023 samtools bcftools 2>/dev/null || true
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"
[[ -f "${BAM}.bai" ]] || samtools index "$BAM"

bcftools mpileup -Ou -f "$REF" -q 10 -Q 10 --threads 8 "$BAM" \
| bcftools call -mv -Oz --threads 8 -o "$OUT"

bcftools index -t "$OUT"
echo "VCF: $OUT"
