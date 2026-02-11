# BINF*6110 Assignment 1: Assembly and Reference-Based Comparison of a *Salmonella enterica* Genome from Oxford Nanopore R10 Reads

## Introduction 

   Pathogenic bacteria employ diverse strategies to invade hosts, evade immune responses, and persist in changing environments. In this context, high-throughput genome sequencing has become central to modern epidemiology, enabling detailed investigations of bacterial genomes at a scale that was previously impractical [1]. However, inaccurate or fragmented assemblies can disrupt coding sequences and inflate apparent variation. As a result, genome assembly remains a foundational step in bacterial genomics and downstream comparative analyses [2].
  *Salmonella enterica* is a gram-negative, facultative bacterium [3]. It is highly adaptable and persists across a wide range of environmental reservoirs and host species, making it a major pathogen of both humans and animals. The *S. enterica* genome is approximately 4.8 Mb in length, exhibits high gene density, and often contains mobile genetic elements such as plasmids that influence adaptation and pathogenicity [5]. These features, combined with the availability of high-quality reference genomes and annotations, make *S. enterica* a strong model organism for evaluating bacterial genome assembly and reference-based comparative genomics [5].
  Reconstructing bacterial genomes from sequencing reads requires minimal fragmentation and high base-level accuracy. Oxford Nanopore Technologies (ONT) long-read sequencing is useful in this setting because it can generate long reads that span repetitive regions and structural variants [2]. This improves assembly continuity relative to other methods such as short-read approaches [6]. ONT sequencing supports real-time base-calling, which can accelerate pathogen characterization in urgent public health situations [11]. Despite these advantages, ONT reads can contain base-calling error and these errors may propagate into draft assemblies if not addressed [7]. Because of this, post-assembly refinement is commonly applied to improve accuracy and to ensure that downstream analyses such as variant calling are not dominated by sequencing induced errors [8].
  Long-read bacterial assembly has also improved due to the development of autonomous assembler programs.  Flye is a long-read assembler designed for error-prone single-molecule data and uses repeat-graph methods to resolve complex repeats and generate contiguous assemblies [7]. Even when Flye produces highly contiguous assemblies, residual read-level errors can remain in the consensus sequence requiring additional refinement. Among refinement tools, Medaka has been shown to produce low-error genomes suitable for downstream comparative analyses [8]. Once a refined consensus genome is obtained, comparison to a high-quality reference enables identification of single nucleotide polymorphisms (SNPs) and other mutations such as insertions and deletions that characterize strain-level differences [9]. For mapping long reads to a reference, minimap2 is widely used due to its accuracy and computational efficiency across large genomic datasets [10]. Finally, visualization tools such as Integrative Genomics Viewer (IGV) enable an inspection of alignments and candidate variants, providing an important step in quality control [12].

## Proposed methods 

### Sequencing data acquisition and preprocessing

  Oxford Nanopore long reads generated using R10 chemistry (expected accuracy Q20+, N50: 5–15 kb) were obtained in FASTQ format from NCBI. To improve assembly performance, ultra-short reads were removed using NanoFilt, because such reads contribute limited long-range structural information and can increase overlaps in long-read assembly workflows.

### Genome assembly

  De novo assembly was performed using Flye (v2.9.2), an autonomous long-read assembler optimized for error-prone single-molecule sequencing reads [7]. Assembly was conducted using the --nano-hq option to specify high-quality Nanopore reads, while other parameters were kept at default settings. Flye’s repeat graph–based approach is intended to improve resolution of repetitive regions and support generation of highly contiguous bacterial genome assemblies [7].

### Assembly polishing

  The draft assembly was polished using racon (v1.5.0) and Medaka (v1.7.3) [8]. Medaka applies neural network models trained on ONT datasets to correct error patterns and improve accuracy [8]. The model corresponding to R10 chemistry was selected, and polishing was run with default settings. Polishing was performed prior to comparative analysis to reduce the likelihood that sequencing errors would appear as false-positive SNPs or indels during reference-based variant detection.

### Reference genome selection

  A high-quality, well-annotated *S. enterica* reference genome (NCBI Assembly ASM694v2) was retrieved from the NCBI Assembly database [5]. This reference was chosen due to its completeness and annotation quality, providing a reliable framework for alignment-based comparison.

### Read alignment to the reference genome

  Refined sequencing reads were aligned to the reference genome using minimap2 (v2.26), a commonly used and benchmarked long-read aligner [10,11]. Mapping was performed using the -ax map-ont preset, which is designed for ONT read characteristics and balances computational speed with alignment accuracy. The resulting alignments were stored in BAM format for downstream analyses.

### Variant identification and visualization

  Aligned reads were used to identify genomic differences between the sequenced sample and the reference genome, with a focus on SNPs and indels as key indicators of strain-level divergence [9]. Alignments and variants were visualized in IGV (v2.16.2), which supports the inspection of mapping quality and sequence information around candidate variants [12].

## Discussion
This study generated a highly contiguous ONT-based assembly for the sequenced Salmonella enterica isolate and used reference-based comparison to quantify divergence from the NCBI ASM694v2 reference. Read preprocessing improved input data suitability for long-read assembly: filtering reduced the dataset from 196,031 reads (809 Mb, min 9 bp, max 58,041 bp) to 164,065 reads (715 Mb, min 2,000 bp, max 17,593 bp) while increasing mean read length (4.1 kb to 4.3 kb). This retained high overall depth and removed ultra-short fragments that contribute little to repeat spanning but can increase overlaps during assembly.

The de novo assembly showed strong contiguity, producing three contigs totaling 5.10 Mb (N50:3.3 MB), with one contig flagged as circular and higher-coverage (215x) than the two larger contigs (139–150x). The circular, high-coverage contig is consistent with a plasmid-like replicon and suggests either increased plasmid copy number relative to the chromosome or preferential sequencing of that element. The assembly graph also showed a short connecting segment between the two large chromosomal components, consistent with an unresolved repeat junction that prevents closure into a single circular chromosome. Polishing with racon/Medaka had a modest effect on total length (−912 bp), implying that polishing corrected local consensus errors rather than restructuring the assembly.

Read mapping statistics further support both overall comparability and the presence of reference-missing sequence. Approximately 94.2% of filtered primary reads mapped to the reference, leaving 5.7% unmapped. Consistent with this, the assembly length (5.1 Mb) exceeds the combined length of the LT2 chromosome and plasmid (4.9 Mb), suggesting that the sample may carry additional or expanded mobile elements relative to the reference and/or that some inserted sequence is present that cannot be represented purely as SNP/short indel differences in a reference-based variant callset. While the VCF summarizes base-level differences where mapping is feasible, it will under-capture novel insertions or accessory regions absent from the reference.

<img width="4764" height="2949" alt="variant_analysis" src="https://github.com/user-attachments/assets/7f860277-0e3f-4d6a-a16d-359ed37bab13" />
Figure 1. Variant calling summary and quality metrics across the Salmonella enterica genome. Long-read reads generated by Oxford Nanopore Technologies were aligned to the S. enterica reference assembly NCBI ASM694v2, and variants were filtered at QUAL ≥ 50. A total of 8,838 variants passed filtering, comprising 8,770 SNPs (99.2%) and 68 indels (0.8%), with no multiallelic sites detected. Variant quality was uniformly high (QUAL min 51.0, max 228.4, mean 224.2), and variant positions show a non-uniform distribution with localized peaks in variant density along the genome.

When comparing the sample to the reference via read mapping and variant calling, the genome showed substantial sequence divergence. After filtering (QUAL ≥ 50), 8,838 total variants were detected, dominated by SNPs (8,770) with few indels (68). The low indel fraction is notable for ONT data, because raw ONT errors are often enriched for indels in low-complexity contexts; here, the small indel count suggests polishing and quality filtering were effective at suppressing common ONT-induced artifacts. Variant QUAL values were uniformly high (mean 224, min 51), and called sites were supported by strong read depth (mean depth at variants 170x versus genome-wide mean depth 137x), increasing confidence that most reported differences reflect genuine strain-level variation.


<img width="4162" height="3816" alt="circos_plot" src="https://github.com/user-attachments/assets/3592b809-aff0-49a4-9ff1-988dad809808" />
Figure 2. Circos view of genome-wide variant distribution across reference replicons. Circos plot summarizing the spatial distribution of filtered variants relative to reference replicons (outer ring; chromosome NC_003197.2 and plasmid NC_003277.2). The middle ring shows binned variant density, while the inner ring plots individual events (SNPs in red, indels in blue). Clustering of points highlights genomic intervals with elevated divergence from the reference, consistent with localized strain differences and/or mobile element-associated variation.


Variant distribution was non-uniform across the genome, with localized peaks in variant density. This pattern is consistent with *S. enterica* genomes where divergence can be concentrated in horizontally acquired regions, as well as in loci under different selective pressures. The Circos and density plots emphasize that divergence is not evenly spread across the chromosome and plasmid replicons, supporting the interpretation that the sequenced isolate is not genetically very similar to the reference ASM694v2 but instead represents a distinct strain/lineage. The transition/transversion ratio (Ts/Tv = 1.1) shows an enrichment of transitions, which is consistent with biological mutation patterns. The Ts/Tv ratio of 1.1 is lower than typical biological expectations of 2-3, which may reflect a combination of the ONT error profile or technical influences from the variant calling process.


<img width="3000" height="3000" alt="bandage_labeled" src="https://github.com/user-attachments/assets/5d9be72b-de1b-42b1-8057-ada3380996d8" />
Figure 3. Assembly graph topology indicates a highly contiguous genome with one circular element and an unresolved repeat junction. Assembly graph visualization (Bandage) of the Flye de novo assembly shows two large chromosomal components connected by a short segment, plus a separate circular element consistent with a plasmid-like replicon.


<img width="1915" height="1107" alt="igv 1" src="https://github.com/user-attachments/assets/6480e696-a14e-48d3-a98c-35213c953e7d" />
Figure 4. IGV validation of a representative high-confidence variant cluster within the dcoB locus. Snapshot from IGV showing aligned ONT reads and the corresponding filtered VCF track. Multiple SNPs and small indel signals are supported by consistent read evidence across the pileup, and the gene model track indicates the variants fall within the coding region of dcoB, illustrating how sequence differences relative to the reference remain after polishing and filtering.

Finally, inspection in IGV provides locus-level validation that called variants have coherent read support and occur within annotated genes. For example, the dcoB region shows clustered mismatches, illustrating how reference divergence can translate into coding changes that may affect protein sequence and function. However, functional interpretation requires downstream annotation of variant effects, which is not captured by counts alone.


<img width="1918" height="1110" alt="igc as" src="https://github.com/user-attachments/assets/816955fd-598a-42c4-b0a5-b17d364ee9d1" />
Figure 5. Broad IGV view of alignment/variant landscape across a representative 582 kb chromosomal interval. Genome browser view showing read alignment depth and dense SNV/indel signals relative to the reference. Mean sequencing depth across the alignment was 137×, with called variant sites supported by higher average depth (170× at variant positions), indicating that most reported variants occur in well-covered regions. This panel provides a quality-control overview linking coverage structure to variant calls over a large genomic span.

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


   







