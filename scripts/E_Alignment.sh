#!/usr/bin/env bash
#module load StdEnv/2023
#module load minimap2
#module load samtools
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/507/875/GCF_002507875.2_ASM250787v2/GCF_002507875.2_ASM250787v2_genomic.fna.gz
#gunzip GCF_002507875.2_ASM250787v2_genomic.fna.gz

set -euo pipefail

ACC=SRR32410565

# refined reads
ASM="medaka_${ACC}/consensus.fasta"
REF="GCF_002507875.2_ASM250787v2_genomic.fna"
OUTDIR="mapping_${ACC}"

mkdir -p "$OUTDIR"

# Align assembly to reference using the asm5 preset, output BAM, then sort + index
minimap2 -t 8 -ax asm5 "${REF}" "$ASM" \
| samtools view -b - \
| samtools sort -o "${OUTDIR}/${ACC}.asm_vs_ref.sorted.bam" -

samtools index "${OUTDIR}/${ACC}.asm_vs_ref.sorted.bam"

echo "Done."
echo "BAM: ${OUTDIR}/${ACC}.asm_vs_ref.sorted.bam"
echo "BAI: ${OUTDIR}/${ACC}.asm_vs_ref.sorted.bam.bai"
