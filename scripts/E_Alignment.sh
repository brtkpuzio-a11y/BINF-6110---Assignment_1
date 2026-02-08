#!/usr/bin/env bash
#module load StdEnv2023
#module load minimap2
#module load samtools
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/507/875/GCF_00250787>
#gunzip GCF_002507875.2_ASM250787v2_genomic.fna.gz

set -euo pipefail

ACC=SRR32410565

# refined reads
READS=medaka_SRR32410565/consensus.fasta
REF=GCF_002507875.2_ASM250787v2_genomic.fna
OUTDIR="mapping_${ACC}"

mkdir -p "$OUTDIR"

# Align ONT reads to reference
minimap2 -t 8 -ax map-ont "${REF}" "$READS" \
| samtools view -b - \
| samtools sort -o "${OUTDIR}/${ACC}.ref.sorted.bam" -

samtools index "${OUTDIR}/${ACC}.ref.sorted.bam"

echo "Done."
echo "BAM: ${OUTDIR}/${ACC}.ref.sorted.bam"
echo "BAI: ${OUTDIR}/${ACC}.ref.sorted.bam.bai"
