#!/usr/bin/env bash
#load module apptainer
#apptainer pull docker://quay.io/biocontainers/nanofilt:2.8.0--py_0

set -euo pipefail

ACC="SRR32410565"
IN="fastq/${ACC}.fastq.gz"
OUT="fastq/${ACC}.nanofilt.q10.l1000.fastq.gz"

MINLEN=1000
MINQ=10
THREADS=8

zcat "$IN" | apptainer exec nanofilt_2.8.0--py_0.sif NanoFilt -l 1000 -q 10 | gzip -c > "$OUT"
echo "Wrote: $OUT"
