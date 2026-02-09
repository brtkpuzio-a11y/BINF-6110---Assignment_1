#!/usr/bin/env bash
#module load racon
#module load minimap2
#apptainer pull docker://quay.io/biocontainers/medaka:2.2.0--py312ha65e1bd_0
set -euo pipefail

ACC=SRR32410565
ASM="flye_${ACC}/assembly.fasta"
READS="fastq/${ACC}.nanofilt.q10.l1000.fastq.gz"
OUTDIR="polished_${ACC}"

mkdir -p "$OUTDIR"

# Racon Round 1
 minimap2 \
  -t 8 -ax map-ont "$ASM" "$READS" > "$OUTDIR/aln1.sam"

 racon \
  -t 8 "$READS" "$OUTDIR/aln1.sam" "$ASM" > "$OUTDIR/racon1.fasta"

# Racon Round 2
 minimap2 \
  -t 8 -ax map-ont "$OUTDIR/racon1.fasta" "$READS" > "$OUTDIR/aln2.sam"

 racon \
  -t 8 "$READS" "$OUTDIR/aln2.sam" "$OUTDIR/racon1.fasta" > "$OUTDIR/racon2.fasta"

# Medaka final polish
apptainer exec medaka_2.2.0--py312ha65e1bd_0.sif medaka_consensus \
  -i "$READS" \
  -d "$OUTDIR/racon2.fasta" \
  -o "$OUTDIR/medaka" \
  -t 8

echo "Final polished assembly: ${OUTDIR}/medaka/consensus.fasta"
