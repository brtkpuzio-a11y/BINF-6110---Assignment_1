#!/usr/bin/env bash
#load module apptainer
#apptainer pull docker://quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0

set -euo pipefail

ACC=SRR32410565
READS="fastq/${ACC}.nanofilt.q10.l1000.fastq.gz"
OUTDIR="flye_${ACC}"

apptainer exec flye_2.9.6--py311h2de2dd3_0.sif flye --nano-hq "$READS" --out-dir "$OUTDIR" --threads 8
echo "Assembly: ${OUTDIR}/assembly.fasta"
