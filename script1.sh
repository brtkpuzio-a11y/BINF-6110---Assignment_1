#!/usr/bin/env bash
#load module apptainer
#apptainer pull docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
#apptainer pull docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1

set -euo pipefail

ACC=SRR32410565
SRA="${ACC}.sra"
FASTQ="${SRA%.sra}.fastq"

# Organization directories
mkdir -p sra fastq qc

# Download SRA file into sra directory
apptainer exec sra-tools_3.2.1--h4304569_1.sif prefetch -O sra "$ACC"

# Convert SRA -> FASTQ
apptainer exec sra-tools_3.2.1--h4304569_1.sif fasterq-dump \
  --threads 8 \
  -O fastq "sra/${ACC}/${ACC}.sra"

# Zip for efficiency
gzip -f fastq/${ACC}*.fastq

# QC
apptainer exec fastqc_0.12.1--hdfd78af_0.sif fastqc \
 --extract -o qc fastq/${ACC}*.fastq.gz

#Print QC
head -n 40 "qc/${ACC}_fastqc/fastqc_data.txt"
head -n 10 "qc/${ACC}_fastqc/summary.txt"
