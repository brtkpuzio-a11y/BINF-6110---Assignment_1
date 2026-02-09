#!/usr/bin/env bash
set -euo pipefail

ACC="SRR32410565"
REF="GCF_002507875.2_ASM250787v2_genomic.fna"
READS="fastq/${ACC}.nanofilt.q10.l1000.fastq.gz"

OUTDIR="vcf_${ACC}_vs_asm"
BAM="${OUTDIR}/${ACC}.vs_asm.sorted.bam"
VCF_RAW="${OUTDIR}/${ACC}.raw.vcf.gz"
VCF_IGV="${OUTDIR}/${ACC}.igv.vcf.gz"

module load StdEnv/2023 minimap2 samtools bcftools 2>/dev/null || true
mkdir -p "$OUTDIR"

# index reference for mpileup + IGV
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

# align reads to your assembly, sort + index
minimap2 -t 8 -ax map-ont "$REF" "$READS" \
| samtools sort -o "$BAM" -
samtools index "$BAM"

# call variants
bcftools mpileup -Ou -f "$REF" -q 30 -Q 15 -d 10000 -a FORMAT/DP --threads 8 "$BAM" \
| bcftools call -mv --ploidy 1 -Oz --threads 8 -o "$VCF_RAW"
bcftools index -t "$VCF_RAW"

# IGV-friendly filter (QUAL>=50 and DP>=10), keep SNPs+indels
bcftools view -v snps,indels "$VCF_RAW" \
| bcftools filter -i 'QUAL>=50 && FORMAT/DP>=10' -Oz -o "$VCF_IGV"
bcftools index -t "$VCF_IGV"

# quick stats
bcftools stats "$VCF_IGV" | grep -E '^SN|TSTV'

echo "REF: $REF"
echo "BAM: $BAM"
echo "VCF raw: $VCF_RAW"
echo "VCF IGV: $VCF_IGV"
