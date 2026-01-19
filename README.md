# BINF*6110 Assignment 1: Assembly and Reference-Based Comparison of a *Salmonella enterica* Genome from Oxford Nanopore R10 Reads

## Introduction 

   Pathogenic bacteria employ diverse strategies to invade hosts, evade immune responses, and persist in changing environments. In this context, high-throughput genome sequencing has become central to modern microbiology and epidemiology, enabling detailed investigations of bacterial genomes at a scale and speed that was previously impractical [1]. However, inaccurate or fragmented assemblies can disrupt coding sequences, and inflate apparent variation. As a result, genome assembly remains a foundational step in bacterial genomics and downstream comparative analyses [2].
  *Salmonella enterica* is a gram-negative, facultative bacterium [3]. It is highly adaptable and persists across a wide range of environmental reservoirs and host species, making it a major pathogen of both humans and animals. The S. enterica genome is approximately 4.8 Mb in length, exhibits high gene density, and often contains mobile genetic elements such as plasmids that influence adaptation and pathogenicity [5]. These biological features, combined with the extensive availability of high-quality reference genomes and annotations, make S. enterica a strong model for evaluating bacterial genome assembly and reference-based comparative genomics [5].
  Reconstructing bacterial genomes from sequencing reads requires minimal fragmentation and high base-level accuracy. Oxford Nanopore Technologies (ONT) long-read sequencing is useful in this setting because it can generate long reads that span repetitive regions and structural variants [2]. This improves assembly continuity relative to other methods such as short-read approaches [6]. ONT sequencing supports real-time base-calling, which can accelerate pathogen characterization in urgent public health situations [12]. Despite these advantages, ONT reads can contain base-calling error and these errors may propagate into draft assemblies if not addressed [7]. Because of this, post-assembly refinement is commonly applied to improve consensus accuracy and to ensure that downstream analyses such as variant calling are not dominated by sequencing induced errors [8].
  Long-read bacterial assembly has also improved due to the development of autonomous assembler programs.  Flye is a long-read assembler designed for error-prone single-molecule data and uses repeat-graph methods to resolve complex repeats and generate contiguous assemblies [7]. Even when Flye produces highly contiguous assemblies, residual read-level errors can remain in the consensus sequence requiring additional refinement. Among refinement tools, Medaka has been shown to produce low-error genomes suitable for downstream comparative analyses [8]. Once a refined consensus genome is obtained, comparison to a high-quality reference enables identification of single nucleotide polymorphisms (SNPs) and other mutations such as insertions and deletions that characterize strain-level differences [9]. For mapping long reads to a reference, minimap2 is widely used due to its accuracy and computational efficiency across large genomic datasets [10]. Finally, visualization tools such as Integrative Genomics Viewer (IGV) enable an inspection of alignments and candidate variants, providing an important step in quality control [13].

## Proposed methods 

### Sequencing data acquisition and preprocessing

  Oxford Nanopore long reads generated using R10 chemistry (expected accuracy Q20+, N50: 5–15 kb) were obtained in FASTQ format from NCBI. Read quality and length distributions were summarized using NanoPlot to confirm that the dataset contained sufficient long-read content for assembly. To reduce noise and improve assembly performance, ultra-short reads were removed using NanoFilt, because such reads contribute limited long-range structural information and can increase spurious overlaps in long-read assembly workflows.

### Genome assembly

  De novo assembly was performed using Flye (v2.9.2), an autonomous long-read assembler optimized for error-prone single-molecule sequencing reads [7]. Assembly was conducted using the --nano-hq option to specify high-quality Nanopore reads, while other parameters were kept at default settings. Flye’s repeat graph–based approach is intended to improve resolution of repetitive regions and support generation of highly contiguous bacterial genome assemblies [7].

### Assembly polishing

  The draft assembly was polished using Medaka (v1.7.3) [8]. Medaka applies neural network models trained on ONT datasets to correct error patterns and improve accuracy [8]. The model corresponding to R10 chemistry was selected, and polishing was run with default settings. Polishing was performed prior to comparative analysis to reduce the likelihood that sequencing errors would appear as false-positive SNPs or indels during reference-based variant detection.

### Reference genome selection

  A high-quality, well-annotated S. enterica reference genome (NCBI Assembly ASM250787v2) was retrieved from the NCBI Assembly database [5]. This reference was chosen due to its completeness and annotation quality, providing a reliable framework for alignment-based comparison.

### Read alignment to the reference genome

  Refined sequencing reads were aligned to the reference genome using minimap2 (v2.26), a commonly used and benchmarked long-read aligner [10,12]. Mapping was performed using the -ax map-ont preset, which is designed for ONT read characteristics and balances computational speed with alignment accuracy. The resulting alignments were stored in BAM format for downstream analyses.

### Variant identification and visualization

  Aligned reads were used to identify genomic differences between the sequenced sample and the reference genome, with a focus on SNPs and indels as key indicators of strain-level divergence [9]. Alignments and variants were visualized in IGV (v2.16.2), which supports the inspection of mapping quality and sequence information around candidate variants [13].

## References








