#!/usr/bin/env bash
#Module load samtools
#Module load bcftools
set -euo pipefail

ACC="${1:-SRR32410565}"
REF="GCF_000006945.2_ASM694v2_genomic.fna"
READS="fastq/${ACC}.nanofilt.q10.l1000.fastq.gz"
OUTDIR="variants_${ACC}_vs_ref"

BAM="${OUTDIR}/${ACC}.vs_ref.sorted.bam"
VCF_RAW="${OUTDIR}/${ACC}.raw.vcf.gz"
VCF_FILT="${OUTDIR}/${ACC}.filtered.vcf.gz"

mkdir -p "$OUTDIR"

# Index reference
if [[ ! -f "${REF}.fai" ]]; then
   samtools faidx "$REF"
fi

# Align reads to reference genome
 minimap2 -t 8 -ax map-ont "$REF" "$READS" | \
 samtools sort -@ 8 -o "$BAM" -

 samtools index "$BAM"

# Alignment statistics
echo "Alignment Statistics"
 samtools flagstat "$BAM"

# Call variants
 bcftools mpileup \
    -Ou -f "$REF" \
    -q 10 -Q 15 \
    -d 10000 \
    -a FORMAT/DP,FORMAT/AD \
    --threads 8 "$BAM" | \
 bcftools call \
    -mv --ploidy 1 \
    -Oz --threads 8 \
    -o "$VCF_RAW"

 bcftools index -t "$VCF_RAW"

# Filter variants
 bcftools view -v snps,indels "$VCF_RAW" | \
 bcftools filter \
    -i 'QUAL>=50 && FORMAT/DP>=20' \
    -Oz -o "$VCF_FILT"

 bcftools index -t "$VCF_FILT"

 # Variant statistics
echo "Variant Statistics"
 bcftools stats "$VCF_FILT" | grep -E '^SN|TSTV'

echo "Reference: $REF"
echo "BAM: $BAM"
echo "VCF raw: $VCF_RAW"
echo "VCF filtered: $VCF_FILT"
echo "View variants with: bcftools view $VCF_FILT | less"
