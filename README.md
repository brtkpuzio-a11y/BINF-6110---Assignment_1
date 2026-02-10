# BINF*6110 Assignment 1: Assembly and Reference-Based Comparison of a *Salmonella enterica* Genome from Oxford Nanopore R10 Reads

## Introduction 

   Pathogenic bacteria employ diverse strategies to invade hosts, evade immune responses, and persist in changing environments. In this context, high-throughput genome sequencing has become central to modern epidemiology, enabling detailed investigations of bacterial genomes at a scale that was previously impractical [1]. However, inaccurate or fragmented assemblies can disrupt coding sequences and inflate apparent variation. As a result, genome assembly remains a foundational step in bacterial genomics and downstream comparative analyses [2].
  *Salmonella enterica* is a gram-negative, facultative bacterium [3]. It is highly adaptable and persists across a wide range of environmental reservoirs and host species, making it a major pathogen of both humans and animals. The *S. enterica* genome is approximately 4.8 Mb in length, exhibits high gene density, and often contains mobile genetic elements such as plasmids that influence adaptation and pathogenicity [5]. These features, combined with the availability of high-quality reference genomes and annotations, make *S. enterica* a strong model organism for evaluating bacterial genome assembly and reference-based comparative genomics [5].
  Reconstructing bacterial genomes from sequencing reads requires minimal fragmentation and high base-level accuracy. Oxford Nanopore Technologies (ONT) long-read sequencing is useful in this setting because it can generate long reads that span repetitive regions and structural variants [2]. This improves assembly continuity relative to other methods such as short-read approaches [6]. ONT sequencing supports real-time base-calling, which can accelerate pathogen characterization in urgent public health situations [11]. Despite these advantages, ONT reads can contain base-calling error and these errors may propagate into draft assemblies if not addressed [7]. Because of this, post-assembly refinement is commonly applied to improve accuracy and to ensure that downstream analyses such as variant calling are not dominated by sequencing induced errors [8].
  Long-read bacterial assembly has also improved due to the development of autonomous assembler programs.  Flye is a long-read assembler designed for error-prone single-molecule data and uses repeat-graph methods to resolve complex repeats and generate contiguous assemblies [7]. Even when Flye produces highly contiguous assemblies, residual read-level errors can remain in the consensus sequence requiring additional refinement. Among refinement tools, Medaka has been shown to produce low-error genomes suitable for downstream comparative analyses [8]. Once a refined consensus genome is obtained, comparison to a high-quality reference enables identification of single nucleotide polymorphisms (SNPs) and other mutations such as insertions and deletions that characterize strain-level differences [9]. For mapping long reads to a reference, minimap2 is widely used due to its accuracy and computational efficiency across large genomic datasets [10]. Finally, visualization tools such as Integrative Genomics Viewer (IGV) enable an inspection of alignments and candidate variants, providing an important step in quality control [12].

## Proposed methods 

### Sequencing data acquisition and preprocessing

  Oxford Nanopore long reads generated using R10 chemistry (expected accuracy Q20+, N50: 5–15 kb) were obtained in FASTQ format from NCBI. Read quality and length distributions were summarized using NanoPlot to confirm that the dataset contained sufficient long-read content for assembly. To improve assembly performance, ultra-short reads were removed using NanoFilt, because such reads contribute limited long-range structural information and can increase overlaps in long-read assembly workflows.

### Genome assembly

  De novo assembly was performed using Flye (v2.9.2), an autonomous long-read assembler optimized for error-prone single-molecule sequencing reads [7]. Assembly was conducted using the --nano-hq option to specify high-quality Nanopore reads, while other parameters were kept at default settings. Flye’s repeat graph–based approach is intended to improve resolution of repetitive regions and support generation of highly contiguous bacterial genome assemblies [7].

### Assembly polishing

  The draft assembly was polished using Medaka (v1.7.3) [8]. Medaka applies neural network models trained on ONT datasets to correct error patterns and improve accuracy [8]. The model corresponding to R10 chemistry was selected, and polishing was run with default settings. Polishing was performed prior to comparative analysis to reduce the likelihood that sequencing errors would appear as false-positive SNPs or indels during reference-based variant detection.

### Reference genome selection

  A high-quality, well-annotated *S. enterica* reference genome (NCBI Assembly ASM694v2) was retrieved from the NCBI Assembly database [5]. This reference was chosen due to its completeness and annotation quality, providing a reliable framework for alignment-based comparison.

### Read alignment to the reference genome

  Refined sequencing reads were aligned to the reference genome using minimap2 (v2.26), a commonly used and benchmarked long-read aligner [10,11]. Mapping was performed using the -ax map-ont preset, which is designed for ONT read characteristics and balances computational speed with alignment accuracy. The resulting alignments were stored in BAM format for downstream analyses.

### Variant identification and visualization

  Aligned reads were used to identify genomic differences between the sequenced sample and the reference genome, with a focus on SNPs and indels as key indicators of strain-level divergence [9]. Alignments and variants were visualized in IGV (v2.16.2), which supports the inspection of mapping quality and sequence information around candidate variants [12].

## References

1.	Donkor E. S. (2013). Sequencing of bacterial genomes: principles and insights into pathogenesis and development of antibiotics. Genes, 4(4), 556–572. https://doi.org/10.3390/genes4040556
2.	Kumar, M. S. et al. (2025). Benchmarking long-read assembly tools and preprocessing strategies for bacterial genomes: A case study on E. coli DH5α. Biotechnology reports (Amsterdam, Netherlands), 48, e00931. https://doi.org/10.1016/j.btre.2025.e00931
3.	Andino, A., Hanning, I., Salmonella enterica: Survival, Colonization, and Virulence Differences among Serovars, The Scientific World Journal, 2015, 520179, 16 pages, 2015. https://doi.org/10.1155/2015/520179
4.	Brown EW, Bell R, Zhang G, Timme R, Zheng J, Hammack TS, Allard MW. 2021.Salmonella Genomics in Public Health and Food Safety. 9: eESP-0008-2020. https://doi.org/10.1128/ecosalplus.ESP-0008-2020
5.	 NCBI. Salmonella enterica genome assembly ASM250787v2. NCBI Assembly database. https://www.ncbi.nlm.nih.gov/assembly/ASM250787v2
6.	Wick, R. R., Judd, L. M., & Holt, K. E. (2023). Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing. PLoS computational biology, 19(3), e1010905. https://doi.org/10.1371/journal.pcbi.1010905
7.	Kolmogorov, M., Yuan, J., Lin, Y. et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol 37, 540–546 (2019). https://doi.org/10.1038/s41587-019-0072-8
8.	Lee, J.Y., Kong, M., Oh, J. et al. Comparative evaluation of Nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis. Sci Rep 11, 20740 (2021). https://doi.org/10.1038/s41598-021-00178-w
9.	Olson, N. D. et al. (2015). Best practices for evaluating single nucleotide variant calling methods for Microbial Genomics. Frontiers in Genetics, 6. https://doi.org/10.3389/fgene.2015.00235
10.	Liyanage, K., Samarakoon, H., Parameswaran, S. et al. Efficient end-to-end long-read sequence mapping using minimap2-fpga integrated with hardware accelerated chaining. Sci Rep 13, 20174 (2023). https://doi.org/10.1038/s41598-023-47354-8
11.	Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191
12.	Robinson, J. T. et al. (2011). Integrative genomics viewer. Nature biotechnology, 29(1), 24–26. https://doi.org/10.1038/nbt.1754


   







