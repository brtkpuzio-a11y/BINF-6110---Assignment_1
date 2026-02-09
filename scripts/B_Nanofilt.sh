#!/usr/bin/env bash
#load module apptainer
#apptainer pull docker://quay.io/biocontainers/nanofilt:2.8.0--py_0

set -euo pipefail

ACC="SRR32410565"
IN="fastq/${ACC}.fastq.gz"
OUT="fastq/${ACC}.nanofilt.q15.l2000.fastq.gz"

zcat "$IN" | apptainer exec nanofilt_2.8.0--py_0.sif NanoFilt -l 2000 -q 15 | gzip -c > "$OUT"
echo "Wrote: $OUT"
