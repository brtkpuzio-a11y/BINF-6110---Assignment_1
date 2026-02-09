#!/usr/bin/env bash
#apptainer pull docker://quay.io/biocontainers/medaka:2.2.0--py312ha65e1bd_0

set -euo pipefail

ACC=SRR32410565
ASM="flye_${ACC}/assembly.fasta"
READS="fastq/SRR32410565.nanofilt.q15.l2000.fastq.gz"
OUTDIR="medaka_${ACC}"

mkdir -p "$OUTDIR"

apptainer exec medaka_2.2.0--py312ha65e1bd_0.sif medaka_consensus \
  -i "$READS" \
  -d "$ASM" \
  -o "$OUTDIR" \
  -t 8

echo "Polished assembly: ${OUTDIR}/consensus.fasta"
