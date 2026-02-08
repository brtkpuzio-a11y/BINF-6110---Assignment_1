#!/usr/bin/env bash
#apptainer pull docker://quay.io/biocontainers/medaka:2.2.0--py312ha65e1bd_0

set -euo pipefail

ACC=SRR32410565
ASM="flye_${ACC}/assembly.fasta"
READS="flye_${ACC}/10-consensus/consensus.fasta"
OUTDIR="medaka_${ACC}"

mkdir -p "$OUTDIR"

apptainer exec medaka_2.2.0--py312ha65e1bd_0.sif medaka_consensus \
  -i "$READS" \
  -d "$ASM" \
  -o "$OUTDIR" \
  -t 8

echo "Polished assembly: ${OUTDIR}/consensus.fasta"
