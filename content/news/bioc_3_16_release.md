November 2, 2022

Bioconductors:

We are pleased to announce Bioconductor 3.16, consisting of
2183 software packages, 416 experiment data packages, 909 annotation
packages, 28 workflows and 8 books.

There are 71 new software packages, 9 new data experiment packages,
1 new annotation packages, no new workflows, no new books, and many updates and
improvements to existing packages.

Bioconductor 3.16 is compatible with R 4.2, and is supported on Linux,
64-bit Windows, and Intel 64-bit macOS 10.13 (High Sierra) or higher.
Bioconductor is excited to start supporting arm64 with this release. This
release will also include updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.16](#getting-started-with-bioconductor-316)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [NEWS from existing workflows](#news-from-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.16
======================================

To update to or install Bioconductor 3.16:

1. Install R 4.2. Bioconductor 3.16 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 71 new software packages in this release of Bioconductor.

- [ATACCoGAPS](/packages/ATACCoGAPS) Provides tools for running the
  CoGAPS algorithm (Fertig et al, 2010) on single-cell ATAC
  sequencing data and analysis of the results. Can be used to perform
  analyses at the level of genes, motifs, TFs, or pathways.
  Additionally provides tools for transfer learning and data
  integration with single-cell RNA sequencing data.

- [ATACseqTFEA](/packages/ATACseqTFEA) Assay for Transpose-Accessible
  Chromatin using sequencing (ATAC-seq) is a technique to assess
  genome-wide chromatin accessibility by probing open chromatin with
  hyperactive mutant Tn5 Transposase that inserts sequencing adapters
  into open regions of the genome. ATACseqTFEA is an improvement of
  the current computational method that detects differential activity
  of transcription factors (TFs). ATACseqTFEA not only uses the
  difference of open region information, but also (or emphasizes) the
  difference of TFs footprints (cutting sites or insertion sites).
  ATACseqTFEA provides an easy, rigorous way to broadly assess TF
  activity changes between two conditions.

- [BASiCStan](/packages/BASiCStan) Provides an interface to infer the
  parameters of BASiCS using the variational inference (ADVI), Markov
  chain Monte Carlo (NUTS), and maximum a posteriori (BFGS) inference
  engines in the Stan programming language. BASiCS is a Bayesian
  hierarchical model that uses an adaptive Metropolis within Gibbs
  sampling scheme. Alternative inference methods provided by Stan may
  be preferable in some situations, for example for particularly
  large data or posterior distributions with difficult geometries.

- [BiocBaseUtils](/packages/BiocBaseUtils) The package provides
  utility functions related to package development. These include
  functions that replace slots, and selectors for show methods. It
  aims to coalesce the various helper functions often re-used
  throughout the Bioconductor ecosystem.

- [BiocFHIR](/packages/BiocFHIR) FHIR R4 bundles in JSON format are
  derived from https://synthea.mitre.org/downloads. Transformation
  inspired by a kaggle notebook published by Dr Alexander Scarlat,
  https://www.kaggle.com/code/drscarlat/fhir-starter-parse-healthcare-bundles-into-tables.
  This is a very limited illustration of some basic parsing and
  reorganization processes.  Additional tooling will be required to
  move beyond the Synthea data illustrations.

- [BioNAR](/packages/BioNAR) the R package BioNAR, developed to step
  by step analysis of PPI network. The aim is to quantify and rank
  each protein’s simultaneous impact into multiple complexes based on
  network topology and clustering. Package also enables estimating of
  co-occurrence of diseases across the network and specific clusters
  pointing towards shared/common mechanisms.

- [cardelino](/packages/cardelino) Methods to infer clonal tree
  configuration for a population of cells using single-cell RNA-seq
  data (scRNA-seq), and possibly other data modalities. Methods are
  also provided to assign cells to inferred clones and explore
  differences in gene expression between clones. These methods can
  flexibly integrate information from imperfect clonal trees inferred
  based on bulk exome-seq data, and sparse variant alleles expressed
  in scRNA-seq data. A flexible beta-binomial error model that
  accounts for stochastic dropout events as well as systematic
  allelic imbalance is used.

- [ccImpute](/packages/ccImpute) Dropout events make the lowly
  expressed genes indistinguishable from true zero expression and
  different than the low expression present in cells of the same
  type. This issue makes any subsequent downstream analysis
  difficult. ccImpute is an imputation algorithm that uses cell
  similarity established by consensus clustering to impute the most
  probable dropout events in the scRNA-seq datasets. ccImpute
  demonstrated performance which exceeds the performance of existing
  imputation approaches while introducing the least amount of new
  noise as measured by clustering performance characteristics on
  datasets with known cell identities.

- [CircSeqAlignTk](/packages/CircSeqAlignTk) CircSeqAlignTk is
  designed for end-to-end RNA-Seq data analysis of circular genome
  sequences, from alignment to visualization. It mainly targets
  viroids which are composed of 246-401 nt circular RNAs. In
  addition, CircSeqAlignTk implements a tidy interface to generate
  synthetic sequencing data that mimic real RNA-Seq data, allowing
  developers to evaluate the performance of alignment tools and
  workflows.

- [consICA](/packages/consICA) consICA implements a data-driven
  deconvolution method – consensus independent component analysis
  (ICA) to decompose heterogeneous omics data and extract features
  suitable for patient diagnostics and prognostics. The method
  separates biologically relevant transcriptional signals from
  technical effects and provides information about the cellular
  composition and biological processes. The implementation of
  parallel computing in the package ensures efficient analysis of
  modern multicore systems.

- [crisprDesign](/packages/crisprDesign) Provides a comprehensive
  suite of functions to design and annotate CRISPR guide RNA (gRNAs)
  sequences. This includes on- and off-target search, on-target
  efficiency scoring, off-target scoring, full gene and TSS
  contextual annotations, and SNP annotation (human only). It
  currently support five types of CRISPR modalities (modes of
  perturbations): CRISPR knockout, CRISPR activation, CRISPR
  inhibition, CRISPR base editing, and CRISPR knockdown. All types of
  CRISPR nucleases are supported, including DNA- and RNA-target
  nucleases such as Cas9, Cas12a, and Cas13d. All types of base
  editors are also supported. gRNA design can be performed on
  reference genomes, transcriptomes, and custom DNA and RNA
  sequences. Both unpaired and paired gRNA designs are enabled.

- [crisprVerse](/packages/crisprVerse) The crisprVerse is a modular
  ecosystem of R packages developed for the design and manipulation
  of CRISPR guide RNAs (gRNAs). All packages share a common language
  and design principles. This package is designed to make it easy to
  install and load the crisprVerse packages in a single step. To
  learn more about the crisprVerse, visit
  <https://www.github.com/crisprVerse>.

- [crisprViz](/packages/crisprViz) Provides functionalities to
  visualize and contextualize CRISPR guide RNAs (gRNAs) on genomic
  tracks across nucleases and applications. Works in conjunction with
  the crisprBase and crisprDesign Bioconductor packages. Plots are
  produced using the Gviz framework.

- [CTSV](/packages/CTSV) The R package CTSV implements the CTSV
  approach developed by Jinge Yu and Xiangyu Luo that detects
  cell-type-specific spatially variable genes accounting for excess
  zeros. CTSV directly models sparse raw count data through a
  zero-inflated negative binomial regression model, incorporates
  cell-type proportions, and performs hypothesis testing based on R
  package pscl. The package outputs p-values and q-values for genes
  in each cell type, and CTSV is scalable to datasets with tens of
  thousands of genes measured on hundreds of spots. CTSV can be
  installed in Windows, Linux, and Mac OS.

- [demuxmix](/packages/demuxmix) A package for demultiplexing
  single-cell sequencing experiments of pooled cells labeled with
  barcode oligonucleotides. The package implements methods to fit
  regression mixture models for a probabilistic classification of
  cells, including multiplet detection. Demultiplexing error rates
  can be estimated, and methods for quality control are provided.

- [DNAfusion](/packages/DNAfusion) Paired-end sequencing of cfDNA
  generated BAM files can be used as input to discover EML4-ALK
  variants. This package was developed using position deduplicated
  BAM files generated with the AVENIO Oncology Analysis Software.
  These files are made using the AVENIO ctDNA surveillance kit and
  Illumina Nextseq 500 sequencing. This is a targeted hybridization
  NGS approach and includes ALK-specific but not EML4-specific
  probes.

- [EasyCellType](/packages/EasyCellType) We developed EasyCellType
  which can automatically examine the input marker lists obtained
  from existing software such as Seurat over the cell
  markerdatabases. Two quantification approaches to annotate cell
  types are provided: Gene set enrichment analysis (GSEA) and a
  modified versio of Fisher's exact test. The function presents
  annotation recommendations in graphical outcomes: bar plots for
  each cluster showing candidate cell types, as well as a dot plot
  summarizing the top 5 significant annotations for each cluster.

- [eds](/packages/eds) This packages provides a single function,
  readEDS. This is a low-level utility for reading in Alevin EDS
  format into R. This function is not designed for end-users but
  instead the package is predominantly for simplifying package
  dependency graph for other Bioconductor packages.

- [EpiMix](/packages/EpiMix) EpiMix is a comprehensive tool for the
  integrative analysis of high-throughput DNA methylation data and
  gene expression data. EpiMix enables automated data downloading
  (from TCGA or GEO), preprocessing, methylation modeling,
  interactive visualization and functional annotation.To identify
  hypo- or hypermethylated CpG sites across physiological or
  pathological conditions, EpiMix uses a beta mixture modeling to
  identify the methylation states of each CpG probe and compares the
  methylation of the experimental group to the control group.The
  output from EpiMix is the functional DNA methylation that is
  predictive of gene expression. EpiMix incorporates specialized
  algorithms to identify functional DNA methylation at various
  genetic elements, including proximal cis-regulatory elements of
  protein-coding genes, distal enhancers, and genes encoding
  microRNAs and lncRNAs.

- [epistasisGA](/packages/epistasisGA) This package runs the GADGETS
  method to identify epistatic effects in nuclear family studies. It
  also provides functions for permutation-based inference and
  graphical visualization of the results.

- [factR](/packages/factR) factR contain tools to process and
  interact with custom-assembled transcriptomes (GTF). At its core,
  factR constructs CDS information on custom transcripts and
  subsequently predicts its functional output. In addition, factR has
  tools capable of plotting transcripts, correcting chromosome and
  gene information and shortlisting new transcripts.

- [FuseSOM](/packages/FuseSOM) A correlation-based multiview
  self-organizing map for the characterization of cell types in
  highly multiplexed in situ imaging cytometry assays (`FuseSOM`) is
  a tool for unsupervised clustering. `FuseSOM` is robust and
  achieves high accuracy by combining a `Self Organizing Map`
  architecture and a `Multiview` integration of correlation based
  metrics. This allows FuseSOM to cluster highly multiplexed in situ
  imaging cytometry assays.

- [gemma.R](/packages/gemma.R) Low- and high-level wrappers for
  Gemma's RESTful API. They enable access to curated expression and
  differential expression data from over 10,000 published studies.
  Gemma is a web site, database and a set of tools for the
  meta-analysis, re-use and sharing of genomics data, currently
  primarily targeted at the analysis of gene expression profiles.

- [GenomAutomorphism](/packages/GenomAutomorphism) This is a R
  package to compute the automorphisms between pairwise aligned DNA
  sequences represented as elements from a Genomic Abelian group. In
  a general scenario, from genomic regions till the whole genomes
  from a given population (from any species or close related species)
  can be algebraically represented as a direct sum of cyclic groups
  or more specifically Abelian p-groups. Basically, we propose the
  representation of multiple sequence alignments of length N bp as
  element of a finite Abelian group created by the direct sum of
  homocyclic Abelian group of prime-power order.

- [ggtreeDendro](/packages/ggtreeDendro) Offers a set of 'autoplot'
  methods to visualize tree-like structures (e.g., hierarchical
  clustering and classification/regression trees) using 'ggtree'. You
  can adjust graphical parameters using grammar of graphic syntax and
  integrate external data to the tree.

- [goSorensen](/packages/goSorensen) This package implements
  inferential methods to compare gene lists (in this first release,
  to prove equivalence) in terms of their biological meaning as
  expressed in the GO. The compared gene lists are characterized by
  cross-tabulation frequency tables of enriched GO items.
  Dissimilarity between gene lists is evaluated using the
  Sorensen-Dice index. The fundamental guiding principle is that two
  gene lists are taken as similar if they share a great proportion of
  common enriched GO items.

- [HiCDOC](/packages/HiCDOC) HiCDOC normalizes intrachromosomal Hi-C
  matrices, uses unsupervised learning to predict A/B compartments
  from multiple replicates, and detects significant compartment
  changes between experiment conditions. It provides a collection of
  functions assembled into a pipeline to filter and normalize the
  data, predict the compartments and visualize the results. It
  accepts several type of data: tabular `.tsv` files, Cooler `.cool`
  or `.mcool` files, Juicer `.hic` files or HiC-Pro `.matrix` and
  `.bed` files.

- [HiContacts](/packages/HiContacts) HiContacts: R interface to
  (m)cool files and other Hi-C processed file formats. HiContacts
  provides a collection of tools to analyse and visualize Hi-C
  datasets. It can import data from pairs or (m)cool files.

- [iSEEhex](/packages/iSEEhex) This package provides panels
  summarising data points in hexagonal bins for `iSEE`. It is part of
  `iSEEu`, the iSEE universe of panels that extend the `iSEE`
  package.

- [iSEEhub](/packages/iSEEhub) This package defines a custom landing
  page for an iSEE app interfacing with the Bioconductor
  ExperimentHub. The landing page allows users to browse the
  ExperimentHub, select a data set, download and cache it, and import
  it directly into a Bioconductor iSEE app.

- [ISLET](/packages/ISLET) ISLET is a method to conduct signal
  deconvolution for general -omics data. It can estimate the
  individual-specific and cell-type-specific reference panels, when
  there are multiple samples observed from each subject. It takes the
  input of the observed mixture data (feature by sample matrix), and
  the cell type mixture proportions (sample by cell type matrix), and
  the sample-to-subject information. It can solve for the reference
  panel on the individual-basis. It can also conduct test to identify
  cell-type-specific differential expression (csDE) genes.

- [katdetectr](/packages/katdetectr) Kataegis refers to the
  occurrence of regional hypermutation and is a phenomenon observed
  in a wide range of malignancies. Using changepoint detection
  katdetectr aims to identify putative kataegis foci from common
  data-formats housing genomic variants. Katdetectr has shown to be a
  robust package for the detection, characterization and
  visualization of kataegis.

- [magrene](/packages/magrene) magrene allows the identification and
  analysis of graph motifs in (duplicated) gene regulatory networks
  (GRNs), including lambda, V, PPI V, delta, and bifan motifs. GRNs
  can be tested for motif enrichment by comparing motif frequencies
  to a null distribution generated from degree-preserving simulated
  GRNs. Motif frequencies can be analyzed in the context of gene
  duplications to explore the impact of small-scale and whole-genome
  duplications on gene regulatory networks. Finally, users can
  calculate interaction similarity for gene pairs based on the
  Sorensen-Dice similarity index.

- [metabinR](/packages/metabinR) Provide functions for performing
  abundance and compositional based binning on metagenomic samples,
  directly from FASTA or FASTQ files. Functions are implemented in
  Java and called via rJava. Parallel implementation that operates
  directly on input FASTA/FASTQ files for fast execution.

- [MetaPhOR](/packages/MetaPhOR) MetaPhOR was developed to enable
  users to assess metabolic dysregulation using transcriptomic-level
  data (RNA-sequencing and Microarray data) and produce
  publication-quality figures. A list of differentially expressed
  genes (DEGs), which includes fold change and p value, from DESeq2
  or limma, can be used as input, with sample size for MetaPhOR, and
  will produce a data frame of scores for each KEGG pathway. These
  scores represent the magnitude and direction of transcriptional
  change within the pathway, along with estimated p-values.MetaPhOR
  then uses these scores to visualize metabolic profiles within and
  between samples through a variety of mechanisms, including: bubble
  plots, heatmaps, and pathway models.

- [MsExperiment](/packages/MsExperiment) Infrastructure to store and
  manage all aspects related to a complete proteomics or metabolomics
  mass spectrometry (MS) experiment. The MsExperiment package
  provides light-weight and flexible containers for MS experiments
  building on the new MS infrastructure provided by the Spectra,
  QFeatures and related packages. Along with raw data
  representations, links to original data files and sample
  annotations, additional metadata or annotations can also be stored
  within the MsExperiment container. To guarantee maximum flexibility
  only minimal constraints are put on the type and content of the
  data within the containers.

- [mslp](/packages/mslp) An integrated pipeline to predict the
  potential synthetic lethality partners (SLPs) of tumour mutations,
  based on gene expression, mutation profiling and cell line genetic
  screens data. It has builtd-in support for data from cBioPortal.
  The primary SLPs correlating with muations in WT and compensating
  for the loss of function of mutations are predicted by random
  forest based methods (GENIE3) and Rank Products, respectively.
  Genetic screens are employed to identfy consensus SLPs leads to
  reduced cell viability when perturbed.

- [MSstatsShiny](/packages/MSstatsShiny) MSstatsShiny is an R-Shiny
  graphical user interface (GUI) integrated with the R packages
  MSstats, MSstatsTMT, and MSstatsPTM. It provides a point and click
  end-to-end analysis pipeline applicable to a wide variety of
  experimental designs. These include data-dependedent acquisitions
  (DDA) which are label-free or tandem mass tag (TMT)-based, as well
  as DIA, SRM, and PRM acquisitions and those targeting
  post-translational modifications (PTMs). The application
  automatically saves users selections and builds an R script that
  recreates their analysis, supporting reproducible data analysis.

- [NetActivity](/packages/NetActivity) #' NetActivity enables to
  compute gene set scores from previously trained sparsely-connected
  autoencoders. The package contains a function to prepare the data
  (`prepareSummarizedExperiment`) and a function to compute the gene
  set scores (`computeGeneSetScores`). The package `NetActivityData`
  contains different pre-trained models to be directly applied to the
  data. Alternatively, the users might use the package to compute
  gene set scores using custom models.

- [octad](/packages/octad) OCTAD provides a platform for virtually
  screening compounds targeting precise cancer patient groups. The
  essential idea is to identify drugs that reverse the gene
  expression signature of disease by tamping down over-expressed
  genes and stimulating weakly expressed ones. The package offers
  deep-learning based reference tissue selection, disease gene
  expression signature creation, pathway enrichment analysis, drug
  reversal potency scoring, cancer cell line selection, drug
  enrichment analysis and in silico hit validation. It currently
  covers ~20,000 patient tissue samples covering 50 cancer types, and
  expression profiles for ~12,000 distinct compounds.

- [omada](/packages/omada) Symptomatic heterogeneity in complex
  diseases reveals differences in molecular states that need to be
  investigated. However, selecting the numerous parameters of an
  exploratory clustering analysis in RNA profiling studies requires
  deep understanding of machine learning and extensive computational
  experimentation. Tools that assist with such decisions without
  prior field knowledge are nonexistent and further gene association
  analyses need to be performed independently. We have developed a
  suite of tools to automate these processes and make robust
  unsupervised clustering of transcriptomic data more accessible
  through automated machine learning based functions. The efficiency
  of each tool was tested with four datasets characterised by
  different expression signal strengths. Our toolkit’s decisions
  reflected the real number of stable partitions in datasets where
  the subgroups are discernible. Even in datasets with less clear
  biological distinctions, stable subgroups with different expression
  profiles and clinical associations were found.

- [oncoscanR](/packages/oncoscanR) The software uses the copy number
  segments from a text file and identifies all chromosome arms that
  are globally altered and computes various genome-wide scores. The
  following HRD scores (characteristic of BRCA-mutated cancers) are
  included: LST, HR-LOH, nLST and gLOH. the package is tailored for
  the ThermoFisher Oncoscan assay analyzed with their Chromosome
  Alteration Suite (ChAS) but can be adapted to any input.

- [PanViz](/packages/PanViz) This pacakge integrates data from the
  Kyoto Encyclopedia of Genes and Genomes (KEGG) with summary-level
  genome-wide association (GWAS) data, such as that provided by the
  GWAS Catalog or GWAS Central databases, or a user's own study or
  dataset, in order to produce biological networks, termed IMONs
  (Integrated Multi-Omic Networks). IMONs can be used to analyse
  trait-specific polymorphic data within the context of biochemical
  and metabolic reaction networks, providing greater biological
  interpretability for GWAS data.

- [phenomis](/packages/phenomis) The 'phenomis' package provides
  methods to perform post-processing (i.e. quality control and
  normalization) as well as univariate statistical analysis of single
  and multi-omics data sets. These methods include quality control
  metrics, signal drift and batch effect correction, intensity
  transformation, univariate hypothesis testing, but also clustering
  (as well as annotation of metabolomics data). The data are handled
  in the standard Bioconductor formats (i.e. SummarizedExperiment and
  MultiAssayExperiment for single and multi-omics datasets,
  respectively; the alternative ExpressionSet and MultiDataSet
  formats are also supported for convenience). As a result, all
  methods can be readily chained as workflows. The pipeline can be
  further enriched by multivariate analysis and feature selection, by
  using the 'ropls' and 'biosigner' packages, which support the same
  formats. Data can be conveniently imported from and exported to
  text files. Although the methods were initially targeted to
  metabolomics data, most of the methods can be applied to other
  types of omics data (e.g., transcriptomics, proteomics).

- [proteasy](/packages/proteasy) Retrieval of experimentally derived
  protease- and cleavage data derived from the MEROPS database.
  Proteasy contains functions for mapping peptide termini to known
  sites where a protease cleaves. This package also makes it possible
  to quickly look up known substrates based on a list of (potential)
  proteases, or vice versa - look up proteases based on a list of
  substrates.

- [RedisParam](/packages/RedisParam) This package provides a
  Redis-based back-end for BiocParallel, enabling an alternative
  mechanism for distributed computation. The The 'manager'
  distributes tasks to a 'worker' pool through a central Redis
  server, rather than directly to workers as with other BiocParallel
  implementations. This means that the worker pool can change
  dynamically during job evaluation. All features of BiocParallel are
  supported, including reproducible random number streams, logging to
  the manager, and alternative 'load balancing' task distributions.

- [regioneReloaded](/packages/regioneReloaded) RegioneReloaded is a
  package that allows simultaneous analysis of associations between
  genomic region sets, enabling clustering of data and the creation
  of ready-to-publish graphs. It takes over and expands on all the
  features of its predecessor regioneR. It also incorporates a
  strategy to improve p-value calculations and normalize z-scores
  coming from multiple analysis to allow for their direct comparison.
  RegioneReloaded builds upon regioneR by adding new plotting
  functions for obtaining publication-ready graphs.

- [RESOLVE](/packages/RESOLVE) Cancer is a genetic disease caused by
  somatic mutations in genes controlling key biological functions
  such as cellular growth and division. Such mutations may arise both
  through cell-intrinsic and exogenous processes, generating
  characteristic mutational patterns over the genome named mutational
  signatures. The study of mutational signatures have become a
  standard component of modern genomics studies, since it can reveal
  which (environmental and endogenous) mutagenic processes are active
  in a tumor, and may highlight markers for therapeutic response.
  Mutational signatures computational analysis presents many
  pitfalls. First, the task of determining the number of signatures
  is very complex and depends on heuristics. Second, several
  signatures have no clear etiology, casting doubt on them being
  computational artifacts rather than due to mutagenic processes.
  Last, approaches for signatures assignment are greatly influenced
  by the set of signatures used for the analysis. To overcome these
  limitations, we developed RESOLVE (Robust EStimation Of mutationaL
  signatures Via rEgularization), a framework that allows the
  efficient extraction and assignment of mutational signatures.
  RESOLVE implements a novel algorithm that enables (i) the efficient
  extraction, (ii) exposure estimation, and (iii) confidence
  assessment during the computational inference of mutational
  signatures.

- [RgnTX](/packages/RgnTX) RgnTX allows the integration of
  transcriptome annotations so as to model the complex alternative
  splicing patterns. It supports the testing of transcriptome
  elements without clear isoform association, which is often the real
  scenario due to technical limitations. It involves functions that
  do permutaion test for evaluating association between features and
  transcriptome regions.

- [scBubbletree](/packages/scBubbletree) scBubbletree is a
  quantitative method for visual exploration of scRNA-seq data. It
  preserves biologically meaningful properties of scRNA-seq data,
  such as local and global cell distances, as well as the density
  distribution of cells across the sample. scBubbletree is scalable
  and avoids the overplotting problem, and is able to visualize
  diverse cell attributes derived from multiomic single-cell
  experiments. Importantly, Importantly, scBubbletree is easy to use
  and to integrate with popular approaches for scRNA-seq data
  analysis.

- [scDDboost](/packages/scDDboost) scDDboost is an R package to
  analyze changes in the distribution of single-cell expression data
  between two experimental conditions. Compared to other methods that
  assess differential expression, scDDboost benefits uniquely from
  information conveyed by the clustering of cells into cellular
  subtypes. Through a novel empirical Bayesian formulation it
  calculates gene-specific posterior probabilities that the marginal
  expression distribution is the same (or different) between the two
  conditions. The implementation in scDDboost treats gene-level
  expression data within each condition as a mixture of negative
  binomial distributions.

- [scifer](/packages/scifer) Have you ever index sorted cells in a 96
  or 384-well plate and then sequenced using Sanger sequencing? If
  so, you probably had some struggles to either check the
  electropherogram of each cell sequenced manually, or when you tried
  to identify which cell was sorted where after sequencing the plate.
  Scifer was developed to solve this issue by performing basic
  quality control of Sanger sequences and merging flow cytometry data
  from probed single-cell sorted B cells with sequencing data. scifer
  can export summary tables, 'fasta' files, electropherograms for
  visual inspection, and generate reports.

- [scMET](/packages/scMET) High-throughput single-cell measurements
  of DNA methylomes can quantify methylation heterogeneity and
  uncover its role in gene regulation. However, technical limitations
  and sparse coverage can preclude this task. scMET is a hierarchical
  Bayesian model which overcomes sparsity, sharing information across
  cells and genomic features to robustly quantify genuine biological
  heterogeneity. scMET can identify highly variable features that
  drive epigenetic heterogeneity, and perform differential
  methylation and variability analyses. We illustrate how scMET
  facilitates the characterization of epigenetically distinct cell
  populations and how it enables the formulation of novel hypotheses
  on the epigenetic regulation of gene expression.

- [ScreenR](/packages/ScreenR) ScreenR is a package suitable to
  perform hit identification in loss of function High Throughput
  Biological Screenings performed using barcoded shRNA-based
  libraries. ScreenR combines the computing power of software such as
  edgeR with the simplicity of use of the Tidyverse metapackage.
  ScreenR executes a pipeline able to find candidate hits from
  barcode counts, and integrates a wide range of visualization modes
  for each step of the analysis.

- [signifinder](/packages/signifinder) signifinder is an R package
  for computing and exploring a compendium of tumor signatures. It
  allows computing signatures scores providing the only gene
  expression values and returns a single-sample score. Currently,
  signifinder contains 46 distinct signatures collected from the
  literature.

- [SimBu](/packages/SimBu) SimBu can be used to simulate bulk RNA-seq
  datasets with known cell type fractions. You can either use your
  own single-cell study for the simulation or the sfaira database.
  Different pre-defined simulation scenarios exist, as are options to
  run custom simulations. Additionally, expression values can be
  adapted by adding an mRNA bias, which produces more biologically
  relevant simulations.

- [simpleSeg](/packages/simpleSeg) Image segmentation is the process
  of identifying the borders of individual objects (in this case
  cells) within an image. This allows for the features of cells such
  as marker expression and morphology to be extracted, stored and
  analysed. simpleSeg provides functionality for user friendly,
  watershed based segmentation on multiplexed cellular images in R
  based on the intensity of user specified protein marker channels.
  simpleSeg can also be used for the normalization of single cell
  data obtained from multiple images.

- [spaSim](/packages/spaSim) A suite of functions for simulating
  spatial patterns of cells in tissue images. Output images are
  multitype point data in SingleCellExperiment format. Each point
  represents a cell, with its 2D locations and cell type. Potential
  cell patterns include background cells, tumour/immune cell
  clusters, immune rings, and blood/lymphatic vessels.

- [SpatialFeatureExperiment](/packages/SpatialFeatureExperiment) A
  new S4 class integrating Simple Features with the R package sf to
  bring geospatial data analysis methods based on vector data to
  spatial transcriptomics. Also implements management of spatial
  neighborhood graphs and geometric operations. This pakage builds
  upon SpatialExperiment and SingleCellExperiment, hence methods for
  these parent classes can still be used.

- [SPIAT](/packages/SPIAT) SPIAT (**Sp**atial **I**mage **A**nalysis
  of **T**issues) is an R package with a suite of data processing,
  quality control, visualization and data analysis tools. SPIAT is
  compatible with data generated from single-cell spatial proteomics
  platforms (e.g. OPAL, CODEX, MIBI, cellprofiler). SPIAT reads
  spatial data in the form of X and Y coordinates of cells, marker
  intensities and cell phenotypes. SPIAT includes six analysis
  modules that allow visualization, calculation of cell
  colocalization, categorization of the immune microenvironment
  relative to tumor areas, analysis of cellular neighborhoods, and
  the quantification of spatial heterogeneity, providing a
  comprehensive toolkit for spatial data analysis.

- [SpliceWiz](/packages/SpliceWiz) Reads and fragments aligned to
  splice junctions can be used to quantify alternative splicing
  events (ASE). However, overlapping ASEs can confound their
  quantification. SpliceWiz quantifies ASEs, calculating
  percent-spliced-in (PSI) using junction reads, and intron retention
  using IRFinder-based quantitation. Novel filters identify ASEs that
  are relatively less confounded by overlapping events, whereby PSIs
  can be calculated with higher confidence. SpliceWiz is ultra-fast,
  using multi-threaded processing of BAM files. It can be run using a
  graphical user or command line interfaces. GUI-based interactive
  visualization of differential ASEs, including novel group-based
  RNA-seq coverage visualization, simplifies short-read RNA-seq
  analysis in R.

- [SpotClean](/packages/SpotClean) SpotClean is a computational
  method to adjust for spot swapping in spatial transcriptomics data.
  Recent spatial transcriptomics experiments utilize slides
  containing thousands of spots with spot-specific barcodes that bind
  mRNA. Ideally, unique molecular identifiers at a spot measure
  spot-specific expression, but this is often not the case due to
  bleed from nearby spots, an artifact we refer to as spot swapping.
  SpotClean is able to estimate the contamination rate in observed
  data and decontaminate the spot swapping effect, thus increase the
  sensitivity and precision of downstream analyses.

- [Statial](/packages/Statial) Statial is a suite of functions for
  identifying changes in cell state. The functionality provided by
  Statial provides robust quantification of cell type localisation
  which are invariant to changes in tissue structure. In addition to
  this Statial uncovers changes in marker expression associated with
  varying levels of localisation. These features can be used to
  explore how the structure and function of different cell types may
  be altered by the agents they are surrounded with.

- [stJoincount](/packages/stJoincount) stJoincount, the application
  of join count analysis to the spatial transcriptomics dataset. This
  tool converts the spatial map into a raster object (a
  two-dimensional image as a rectangular matrix or grid of square
  pixels), where clusters labelled spots are converted to adjacent
  pixels with a calculated resolution. A neighbors' list was created
  based on the rasterized sample, which identifies adjacent and
  diagonal neighbors for each pixel. After adding binary spatial
  weights to the neighbors' list, a multi-categorical join count
  analysis is then performed, allowing all possible combinations of
  cluster pairings to be tabulated. The function returns the observed
  join counts, the expected count under conditions of spatial
  randomness, and the variance of observed to expected calculated
  under non-free sampling. The z-score is then calculated as the
  difference between observed and expected counts, divided by the
  square root of the variance.

- [SUITOR](/packages/SUITOR) An unsupervised cross-validation method
  to select the optimal number of mutational signatures. A data set
  of mutational counts is split into training and validation
  data.Signatures are estimated in the training data and then used to
  predict the mutations in the validation data.

- [syntenet](/packages/syntenet) syntenet can be used to infer
  synteny networks from whole-genome protein sequences and analyze
  them. Anchor pairs are detected with the MCScanX algorithm, which
  was ported to this package with the Rcpp framework for R and C++
  integration. Anchor pairs from synteny analyses are treated as an
  undirected unweighted graph (i.e., a synteny network), and users
  can perform: i. network clustering; ii. phylogenomic profiling (by
  identifying which species contain which clusters) and; iii.
  microsynteny-based phylogeny reconstruction with maximum
  likelihood.

- [TENxIO](/packages/TENxIO) Provides a structured S4 approach to
  importing data files from the 10X pipelines. It mainly supports
  Single Cell Multiome ATAC + Gene Expression data among other data
  types. The main Bioconductor data representations used are
  SingleCellExperiment and RaggedExperiment.

- [VDJdive](/packages/VDJdive) This package provides functions for
  handling and analyzing immune receptor repertoire data, such as
  produced by the CellRanger V(D)J pipeline. This includes reading
  the data into R, merging it with paired single-cell data,
  quantifying clonotype abundances, calculating diversity metrics,
  and producing common plots. It implements the E-M Algorithm for
  clonotype assignment, along with other methods, which makes use of
  ambiguous cells for improved quantification.

- [Voyager](/packages/Voyager) Voyager to SpatialFeatureExperiment
  (SFE) is just like scater to SingleCellExperiment. While SFE is a
  new S4 class, Voyager implements basic exploratory spatial data
  analysis (ESDA) methods for SFE. This first version supports
  univariate global spatial ESDA methods such as Moran's I,
  permutation testing for Moran's I, and correlograms. Voyager also
  implements plotting functions to plot SFE data and ESDA results.
  Multivariate ESDA and univariate local metrics will be added in
  later versions.

- [vsclust](/packages/vsclust) Feature-based variance-sensitive
  clustering of omics data. Optimizes cluster assignment by taking
  into account individual feature variance. Includes several modules
  for statistical testing, clustering and enrichment analysis.

- [zenith](/packages/zenith) Zenith performs gene set analysis on the
  result of differential expression using linear (mixed) modeling
  with dream by considering the correlation between gene expression
  traits.  This package implements the camera method from the limma
  package proposed by Wu and Smyth (2012).  Zenith is a simple
  extension of camera to be compatible with linear mixed models
  implemented in variancePartition::dream().

New Data Experiment Packages
=====================

There are 9 new data experiment packages in this release of Bioconductor.


- [BioPlex](/packages/BioPlex) The BioPlex package implements access
  to the BioPlex protein-protein interaction networks and related
  resources from within R. Besides protein-protein interaction
  networks for HEK293 and HCT116 cells, this includes access to CORUM
  protein complex data, and transcriptome and proteome data for the
  two cell lines. Functionality focuses on importing the various data
  resources and storing them in dedicated Bioconductor data
  structures, as a foundation for integrative downstream analysis of
  the data.

- [EpiMix.data](/packages/EpiMix.data) Supporting data for the EpiMix
  R package. It include: - HM450_lncRNA_probes.rda -
  HM450_miRNA_probes.rda - EPIC_lncRNA_probes.rda -
  EPIC_miRNA_probes.rda - EpigenomeMap.rda - LUAD.sample.annotation -
  TCGA_BatchData - MET.data - mRNA.data - microRNA.data - lncRNA.data
  - Sample_EpiMixResults_lncRNA - Sample_EpiMixResults_miRNA -
  Sample_EpiMixResults_Regular - lncRNA expression data of tumors
  from TCGA that are stored in the ExperimentHub.

- [HiContactsData](/packages/HiContactsData) Provides a collection of
  Hi-C files (pairs and (m)cool). These datasets can be read into R
  and further investigated and visualized with the HiContacts
  package. Data includes yeast Hi-C data generated by the Koszul lab
  from the Pasteur Institute.

- [MerfishData](/packages/MerfishData) MerfishData is an
  ExperimentHub package that serves publicly available datasets
  obtained with Multiplexed Error-Robust Fluorescence in situ
  Hybridization (MERFISH). MERFISH is a massively multiplexed
  single-molecule imaging technology capable of simultaneously
  measuring the copy number and spatial distribution of hundreds to
  tens of thousands of RNA species in individual cells. The scope of
  the package is to provide MERFISH data for benchmarking and
  analysis.

- [MicrobiomeBenchmarkData](/packages/MicrobiomeBenchmarkData) The
  MicrobiomeBenchmarkData package provides functionality to access
  microbiome datasets suitable for benchmarking. These datasets have
  some biological truth, which allows to have expected results for
  comparison. The datasets come from various published sources and
  are provided as TreeSummarizedExperiment objects. Currently, only
  datasets suitable for benchmarking differential abundance methods
  are available.

- [NetActivityData](/packages/NetActivityData) This package contains
  the weights from pre-trained shallow sparsely-connected
  autoencoders. This data is required for getting the gene set scores
  with NetActivity package.

- [octad.db](/packages/octad.db) Open Cancer TherApeutic Discovery
  (OCTAD) package implies sRGES approach for the drug discovery. The
  essential idea is to identify drugs that reverse the gene
  expression signature of a disease by tamping down over-expressed
  genes and stimulating weakly expressed ones. The following package
  contains all required precomputed data for whole OCTAD pipeline
  computation.

- [SFEData](/packages/SFEData) Example spatial transcriptomics
  datasets with Simple Feature annotations as
  SpatialFeatureExperiment objects. Technologies include Visium,
  slide-seq, Nanostring CoxMX, Vizgen MERFISH, and 10X Xenium.
  Tissues include mouse skeletal muscle, human melanoma metastasis,
  human lung, breast cancer, and mouse liver.

- [WeberDivechaLCdata](/packages/WeberDivechaLCdata)
  Spatially-resolved transcriptomics (SRT) and single-nucleus
  RNA-sequencing (snRNA-seq) data from the locus coeruleus (LC) in
  postmortem human brain samples. Data were generated with the 10x
  Genomics Visium SRT and 10x Genomics Chromium snRNA-seq platforms.
  Datasets are stored in SpatialExperiment and SingleCellExperiment
  formats.

New Annotation Packages
=====================

There are 2 new annotation packages in this release of Bioconductor.

- [HDO.db](/packages/HDO.db)
  A set of annotation maps describing the entire Human Disease Ontology
  assembled using data from DO. Its annotation data comes from
  https://github.com/DiseaseOntology/HumanDiseaseOntology/tree/main/src/ontology.

- [UniProtKeywords](/packages/UniProtKeywords)
  UniProt database provides a list of controlled vocabulary represented as
  keywords for genes or proteins. This is useful for summarizing gene functions
  in a compact way. This package provides data of keywords hierarchy and
  gene-keyword relations.

New Workflow Packages
=====================

There are no new workflow packages.


New Books
=====================

There are no new online books.


NEWS from existing Software Packages
===================================


[affxparser](/packages/affxparser)
----------

                 Changes in version 1.69.1 (2022-04-28)                 

BUG FIXES

- Ported bug fix from affxparser 1.68.1.

                 Changes in version 1.69.0 (2022-04-26)                 

- The version number was bumped for the Bioconductor devel version,
  which is
  now BioC 3.16 for R-devel.

                 Changes in version 1.68.1 (2022-04-28)                 

BUG FIXES

- affxparser (>= 1.67.1) failed to install with R built with '-fpic'
  flag. The symptom was a linking error 'ld: 000.init.o: relocation
  R_X86_64_32 against `.rodata' can not be used when making a shared
  object; recompile with -fPIC collect2: error: ld returned 1 exit
  status'.

[alevinQC](/packages/alevinQC)
--------

                       Changes in version 1.13.2                        

- Allow processing of alevin-fry data with unfiltered permitlist

                       Changes in version 1.13.1                        

- Bug fix; use nbr mapped reads to calculate fraction of reads in
barcode list for alevin-fry

[AlpsNMR](/packages/AlpsNMR)
-------

                 Changes in version 4.1.1 (2022-11-02)                  

- Remove archive dependency
- Try fixing build on palomino4, due to race condition in R CMD check

                 Changes in version 3.99.7 (2022-10-27)                 

Minor changes

- Fix build issue on palomino4, simplifying helper function

                 Changes in version 3.99.6 (2022-10-26)                 

Minor changes

- When saving, normalize extra information is saved as well.
- Updated downsampled demo data for examples.
- More robust nmr_baseline_threshold()
- Faster examples

                 Changes in version 3.99.5 (2022-10-26)                 

Minor changes

- Remove call to deprecated ggplot2::qplot()

                 Changes in version 3.99.4 (2022-10-19)                 

Major changes

- Improved the download_MTBLS242() function, allowing to either
download the parts of MTBLS242 needed for the tutorial or the whole
dataset, which may be nice to have if you want to play beyond the
tutorial.

- When reading a Bruker sample from a zip file, you now can specify
in
the file name the zip subdirectory. For instance,
"/path/to/sample.zip!/sample/3", when sample.zip contains a folder
named sample with a subfolder named 3 that includes the sample data
you want to actually read.

Minor changes

- Remove Bioconductor Build System workaround, since
https://github.com/Bioconductor/BBS/issues/220 was fixed.

                 Changes in version 3.99.3 (2022-10-17)                 

- Add libarchive as a SystemRequirement to workaround a limitation of
the Bioconductor build system (BBS), that can't pick system
requirements recursively. Thanks to Jennifer Wokaty for checking the
BBS and providing this suggestion.

                 Changes in version 3.99.2 (2022-10-14)                 

Breaking changes

- Set fix_baseline = FALSE in nmr_integrate_regions() as default. The
former TRUE approach here did not make much sense if peak boundaries
were not perfectly established.

Major changes

- Baseline estimation: We now offer nmr_baseline_estimation() besides
nmr_baseline_removal(). The estimation function computes the
baseline and saves it instead of subtracting it from the signal.
This is a better approach because it lets each step of the pipeline
decide whether it makes sense to subtract the baseline or not. The
nmr_baseline_removal() is for now still available, but it will be
deprecated in a future version.

- For the baselineThresh argument in nmr_detect_peaks() we now
suggest
using nmr_baseline_threshold(dataset, method = "median3mad"). This
is more robust than the former (but still the default) method.

- Peak detection and integration: We want to approach the peak
detection, clustering an integration in a different way. While the
old pipeline still works as expected, we have introduced new
arguments to peak detection, with backwards compatible defaults and
a peak clustering function. We still provide the vignette with the
former workflow, because it is still relevant but we may deprecate
it in a future version, once we are confident the changes we are
making are robust across several datasets.

- Parallellization: We are switching from the future package to
BiocParallel, to better integrate in the Bioconductor ecosystem. In
this version, if you use a different future plan you may get a
warning to switch to BiocParallel. In a future version we will
remove our dependency with the (awesome) future package.

Minor changes

- You can now set experiment names (NMRExperiment) with
names(dataset)
<- c("Sample1", "Sample2").
- You can now pass a named vector with the sample names to the
nmr_read_samples function. The names will be used as the sample
names.
- Peak detection has a more robust baseline threshold estimation
- Peak detection estimates the baseline threshold on each sample
individually. The threshold is calculated using only the sample
where we are currently detecting the peaks.
- Peak detection includes a simple but effective lorentzian fitting
(for area and width estimation)
- Add functions to evaluate the quality of the peak detection using
plots
- More fine grained interpolation axis if axis = NULL is given in
nmr_interpolate_1D()
- Save list of excluded regions in the nmr_dataset object.
- Drop MassSpecWavelet workaround on partial argument matching since
it was fixed upstream
- Documentation: Start providing verbose messages with tips in
functions
- Remove unused deprecated imports from the future package (#65,
thanks to @HenrikBengtsson)
- Add URL and BugReports to the DESCRIPTION (#64, thanks to
@HenrikBengtsson)
- Reading bruker samples is now a bit more robust and gives detailed
tracebacks in case of error.

[ANCOMBC](/packages/ANCOMBC)
-------

                  Changes in version 2.0 (2022-10-19)                   

- add ancombc2 function

[AnnotationForge](/packages/AnnotationForge)
---------------

                       Changes in version 1.39.0                        

MODIFICATIONS

- 1.39.1 Update to REST query for uniprot

- 1.39.2 Update viableID.rda for upcoming release

[AnnotationHub](/packages/AnnotationHub)
-------------

                        Changes in version 3.5.0                        

NEW FEATURES

- (3.5.2) Add dcf dispatchclass

- (3.5.1) Add CompDb dispatchclass

[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.10.0                        

NEW FEATURES

- (v 1.9.1) add drs_access_url() to returned signed https:// URLs
from
drs:// URIs. Enhance drs_cp().

- (v 1.9.4) add auto_unbox= argument to Service class, allowing other
developers flexibility in unboxing values passed to REST APIs.

- (v 1.9.7) add developer facilities for tracking API changes in
Rawls, Terra, and Leonardo services

USER VISIBLE CHANGES

- (v 1.9.2) Deprecate AnVIL::install() & friends in favor of
BiocManager::install(), which now knows about container binary
repositories.

- (v 1.9.8) Update Rawls, Terra, and Leonardo services. Changed
endpoints include:

## Rawls
$removed
[1] admin_delete_refresh_token admin_statistics_get
[3] refreshToken refreshTokenDate

$updated
[1] listUserBillingAccounts createWorkspace getTags
[4] clone entity_type_metadata get_entity
[7] entityQuery createSubmission validateSubmission

## Terra
$removed
[1] userTrial listImportPFBJobs importPFBStatus

$updated
[1] deleteBillingProject billingAccounts
[3] createWorkspace cloneWorkspace
[5] entityQuery flexibleImportEntities
[7] importEntities createSubmission
[9] validateSubmission browserDownloadEntitiesTSV
[11] setProfile

## Leonardo
$removed
[1] batchNodepoolCreate

$updated
[1] listApp listAppByProject deleteApp
[4] createApp listDisks listDisksByProject
[7] createDisk updateRuntime createRuntime
[10] setCookie proxyClusterJupyter proxyClusterJupyterLab
[13] proxyClusterRStudio

- (v 1.9.9) add 'gadgets' (simple graphical interfaces) to key
functions, avworkspace_gadget(), avtable_gadget(),
avworkflow_gadget(). Also browse_workspace() for opening a terra
workspace in the browser.

BUG FIXES

- (v 1.9.3 / 1.8.2) avworkflow_localize() looks for submissionId
files
correctly.

- (v 1.9.5 / 1.8.3) drs_stat() works when accessUrl is included in
response.

- (v 1.9.6 / 1.8.5) gsutil_cp() and gsutil_rsync() use
normalizePath()
on source and destination arguments to avoid creating directories in
unexpected locations when provided with paths containing ~, . or ...

- (v 19.10 / v 1.8.6) gcloud_account("<new account>") did not
invalidate cached access tokens.
https://github.com/Bioconductor/AnVIL/issues/66

- (v 1.9.11 / v 1.8.7) avoid changing status of 'Done' workflows to
'Aborted' https://github.com/Bioconductor/AnVIL/issues/64

- (v 1.9.11 / v 1.8.7) allow 'NULL' for entity arguments of
avworkflow_run() https://github.com/Bioconductor/AnVIL/issues/65

[AnVILPublish](/packages/AnVILPublish)
------------

                        Changes in version 1.8.0                        

New Features

- (v. 1.7.1) Support publication output '.Rmd'

- (v. 1.7.2) Support 'Quarto' for Rmd to ipynb conversion

- (v. 1.7.3) Support 'Qmd' files.

User Visible Changes

- (v. 1.7.4) Update workspace publication -- put 'Description'
directly under title so it appears in WORKSPACE summary; add current
'Date' to citation if not already defined.

Bug Fixes

- (v. 1.7.4) Workspace dashboard and startup link to notebooks in
'analysis' rather than 'notebook' location; clean up vignette author
list.

[aroma.light](/packages/aroma.light)
-----------

                 Changes in version 3.27.0 (2022-04-26)                 

- The version number was bumped for the Bioconductor devel version,
  which is now BioC 3.16 for R-devel.

[ATACCoGAPS](/packages/ATACCoGAPS)
----------

                 Changes in version 0.99.0 (2022-03-04)                 

- Submitted to Bioconductor

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.21.2                        

- Add error message for NA seqlengths for factorFootprints.R

                       Changes in version 1.21.1                        

- remove the limits of BSgenome object for vPlot

[AUCell](/packages/AUCell)
------

                        Changes in version 1.19                         

- New function: AUCell_assignCells()

[bandle](/packages/bandle)
------

                         Changes in version 1.1                         

bandle 1.1.3

- add more details in the documentation on how to fitGPs and pass
hyppar to fitGPmaternPC
- updated typos in vignettes and removed calls to unecessary
libraries
- added installation troubleshooting section to README to advise
users
of potential missing binary libraries

bandle 1.1.2

- fix bug in plotConvergence function

bandle 1.1.1

- use gh actions to build site

version 1.1.0

- fix typos
- upstream changes from bioconductor

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.9.7 (2022-09-06)                  

- Bugfix in BASiCS_DetectHVG/BASiCS_DetectLVG.

                 Changes in version 2.9.6 (2022-09-05)                  

- Add TransLogit parameter to BASiCS_PlotDE to have volcano plots on a
  logit scale.
  As the posterior probabilities are in the interval [0, 1] this may
  make it
  easier to visualise difference in "interesting" regions (e.g., 0.6-1)

                 Changes in version 2.9.5 (2022-07-29)                  

- Further bugfix in subset method for BASiCS_Chain objects. subset now
  respects the ordering implied by character or numeric indices, where
  it
  previously did not.

                 Changes in version 2.9.4 (2022-07-28)                  

- Update roxygen version

                 Changes in version 2.9.3 (2022-07-18)                  

- Bugfix for divide and conquer code - handle arguments like
  StoreChains and
  RunName correctly in this context.

                 Changes in version 2.9.2 (2022-07-18)                  

- Implement parallel tests

- Remove variance decomp plotting code from an internal code block to
  the public BASiCS_PlotVarianceDecomp.

                 Changes in version 2.9.1 (2022-05-04)                  

- switch from using DOUBLE_EPS to a function call to
  std::numeric_limits<double>::epsilon() - credit Tomas Kalibera for
  informing by private correspondence.

[BASiCStan](/packages/BASiCStan)
---------

                         Changes in version 1.0                         

- Initial version with support for WithSpikes and Regression modes.
No-spike version is supported with fixed scaling normalisation
factors (scran size factors are used by default).

[beachmat](/packages/beachmat)
--------

                       Changes in version 2.14.0                        

- Added the tatami C++ library for LinkingTo from other packages.

[BEclear](/packages/BEclear)
-------

                 Changes in version 2.13.2 (2022-09-21)                 

- fixed issue with tests

- fixed issue with folders

[beer](/packages/beer)
----

                 Changes in version 1.1.3 (2022-06-26)                  

- Added helper for plotting observed versus expected read counts
  (`getExpected()`)

- Added helper for plotting Bayes factors (`getBF()`)

                 Changes in version 1.1.2 (2022-06-06)                  

- Added `codetools` to `Suggests` for `BiocParallel` dependency.

                 Changes in version 1.1.1 (2022-06-06)                  

- Added support for running `edgeR` with `glmQLFTest` instead of
  `exactTest`.

[benchdamic](/packages/benchdamic)
----------

                 Changes in version 1.3.11 (2022-10-31)                 

- Updated DA_ANCOM with the new ancombc2 method

                 Changes in version 1.3.10 (2022-10-20)                 

- Re-adapted p-values extraction from DA_ALDEx2() with test = "glm"

                 Changes in version 1.3.9 (2022-09-18)                  

- Re-added corncob after temporary removal in v1.3.5

                 Changes in version 1.3.8 (2022-09-18)                  

- Bug-fix for table headers in chapter 3 of the vignette

                 Changes in version 1.3.7 (2022-09-18)                  

- Added more informations about the methods in chapter 3 of the
  vignette

                 Changes in version 1.3.6 (2022-09-10)                  

- Simplified the vignette

                 Changes in version 1.3.5 (2022-09-07)                  

- Temporarily removed corncob. New version of detectseparation v0.3
  broke it.

- Updated intro vignette accordingly.

                 Changes in version 1.3.4 (2022-09-07)                  

- Updated intro vignette.

- Removed HMP16SData chunks.

- Added conditional usage of kableExtra

- Removed dependency from HMP16SData and ffpe

- Implemented CAT function instead of using ffpe::CATplot

                 Changes in version 1.3.3 (2022-09-06)                  

- In DA_basic changed logFC calculation when paired = TRUE

- In intro vignette corrected the number of comparisons from 10 to 2

                 Changes in version 1.3.2 (2022-09-03)                  

- Changes after the first revision from the "Bioinformatics" journal

- Improved vignette with more descriptions, new methods, and more
  chapters

- Improved warning and error messages

- Support for TreeSummarizedExperiment object

- New function get_count_metadata() to access counts and metadata of a
  phyloseq or TreeSummarizedExperiment object.

- Added the "fitFeatureModel" implementation in DA_metagenomeSeq()

- "CSS" normalization factors are now returned as they are, not divided
  by a scaling factor

- Replaced "CSS_median" and "CSS_default" normalizations with "CSS"

- Added the "kw" and "glm" implementation in DA_ALDEx2()

- Added the "paired.test" option in DA_ALDEx2() for "t" and "wilcox"
  tests

- Added support for "CLR", "RC", and "none" normalization methods in
  DA_Seurat()

- Added support for "scale.factor" in DA_Seurat()

- Added more controls and verbosity regarding normalizations and DA
  methods

- Added the plotLogP() function for plottin the negative log p-values
  distribution

- Added DA_NOISeq(), DA_dearseq(), DA_ANCOM(), and DA_basic() to the
  pool of methods

- added detectseparation (<= 0.2.0) to the "Imports", new version 0.3
  not compatible with corncob 0.2

                 Changes in version 1.3.1 (2022-05-16)                  

- Updated intro vignette

- Added more citations and authors

- Added figure captions and references

                 Changes in version 1.3.0 (2022-04-26)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_15 branch

                 Changes in version 1.2.5 (2022-09-10)                  

- Porting the changes of devel version 1.3.6 to release

                 Changes in version 1.2.4 (2022-09-07)                  

- Porting the changes of devel version 1.3.5 to release

                 Changes in version 1.2.3 (2022-09-07)                  

- Porting the changes of devel version 1.3.4 to release

                 Changes in version 1.2.2 (2022-09-06)                  

- Porting the changes of devel version 1.3.3 to release

                 Changes in version 1.2.1 (2022-09-03)                  

- Porting the changes of devel version 1.3.2 to release

[BindingSiteFinder](/packages/BindingSiteFinder)
-----------------

                         Changes in version 1.4                         

- major re-work of the vignette

- implemented bindingSiteCoveragePlot() function

- added subsetting by index functionalities

- parameters minCrosslinks and minClSites can be deactivated when set
  to 0

- coverageOverRanges() function fails with error message on subscript
  out of
  bounds

- fixed bug in ReproducibilityFilter() function

- fix silent option in BSFDataSet constructor

[BiocBaseUtils](/packages/BiocBaseUtils)
-------------

                        Changes in version 1.0.0                        

- Added a NEWS.md file to track changes to the package.

[BiocCheck](/packages/BiocCheck)
---------

                        Changes in version 1.34                         

NEW FEATURES

- Redundant package dependencies checks between `DESCRIPTION` and
  `NAMESPACE` have been removed. These are already present in `R CMD
  check`
  as "checking package dependencies".

- Use `callr` to run `BiocCheck` in a separate process, this avoids
  interference with loaded packages (@vjcitn, #158)

BUG FIXES AND MINOR IMPROVEMENTS

- Update `checkVigInstalls` and `checkVigBiocInst` to avoid false
  positives
  (@almeidasilvaf, #170).

- Only count non evaluated chunks when there are any present in the
  vignette

- Fix false positive WARNING "Import {pkg} in NAMESPACE as well as
  DESCRIPTION." where pkg was not in NAMESPACE but it was used using
  double
  colons pkg::function inside an S4 method. (@zeehio, #166)

- Fix bug where inputs to `getDirFile` were vectors in
  `checkForValueSection` (@zeehio, #163)

- Allow lookback for matching T/F and exclude list elements, e.g.,
  `list$F`
  (@lshep, #161)

- Fix indentation count by excluding yaml front matter from vignettes
  (@harpomaxx, #100)

- Update internal documentation of the `BiocCheck-class`

- Fix bug where line numbers were off due to removal of empty lines at
  parsing (@lshep, #159)

- Slightly improve sentence counter for Description field check
  (@lshep, #160)

- Update documentation links to point to contributions.bioconductor.org
  (@LiNk-NY, #157)

[BiocFHIR](/packages/BiocFHIR)
--------

                        Changes in version 1.0.0                        

- Specifically works with JSON data from
  https://synthea.mitre.org/downloads

[BiocFileCache](/packages/BiocFileCache)
-------------

                         Changes in version 2.5                         

ENHANCEMENT

- (2.5.2) Add option to override unique identifer when adding files to
  the
  cache. This will allow exact match of original file name.

BUG FIX

- (2.5.1) ERROR if interactively decide not to removebfc

[BiocParallel](/packages/BiocParallel)
------------

                        Changes in version 1.32                         

NEW FEATURES

- (v 1.31.10) bpiterate() when ITER is not a function will use
  bpiterateAlong() to attempt to iterate over elements ITER[[1]],
  ITER[[2]], etc.
  https://stat.ethz.ch/pipermail/bioc-devel/2022-July/019075.html

USER VISIBLE CHANGES

- (v 1.31.3) Deprecate BatchJobsParam in favor of BatchtoolsParam

- (v 1.31.11) Replace Random Number .Rnw vignette with Rmd (html)
  version (thanks Madelyn Carlson!)
  https://github.com/Bioconductor/BiocParallel/pull/215

- (v 1.31.12) clarify default number of cores, and use on shared
  clusters (thanks Dario Strbenac)
  https://github.com/Bioconductor/BiocParallel/pull/218
  https://github.com/Bioconductor/BiocParallel/issues/217

- (v 1.31.15) Replace Introduction to BiocParallel .Rnw vignette
  with Rmd (html) version (thanks Phylis Atieno!)
  https://github.com/Bioconductor/BiocParallel/pull/226

BUG FIXES

- (v 1.31.1) suppress package startup messages on workers
  https://github.com/Bioconductor/BiocParallel/issues/198

- (v 1.31.1) coerce timeout to integer (typically from numeric)
  https://github.com/Bioconductor/BiocParallel/issues/200

- (v 1.31.2) avoid segfault when ipcmutex() functions generate C++
  errors. This occurs very rarely, for instance when the directory
  used by boost for file locking (under /tmp) was created by another
  user.
  https://github.com/Bioconductor/BiocParallel/pull/202

- (v 1.31.2) resetting bpRNGseed() after bpstart() is reproducible
  https://github.com/Bioconductor/BiocParallel/pull/204

- (v 1.31.5) enable logs for multiple managers sharing the same
  workers.
  https://github.com/Bioconductor/BiocParallel/pull/207

- (v 1.31.13 / v 1.30.4) only export variables in `.GlobalEnv` or
  `package:`
  <https://github.com/Bioconductor/BiocParallel/issues/223>
  <https://github.com/Bioconductor/BiocParallel/pull/224>

- (v 1.31.14) Reduce bpmapply memory usage. Thanks Sergio Oller.
  <https://github.com/Bioconductor/BiocParallel/pull/227>

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.16.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- `get_cre_orcids` and `orcid_table` allow the querying of the ORCID
  API
  to obtain information about the maintainer whose ORCID is documented
  in the `DESCRIPTION` file (@vjcitn)

- `biocBuildEmail` now allows a `to` maintainer email override along
  with
  `cc` and `bcc` inputs.

- `biocLastBuildDate` gives the date of the last release build for a
  particular release version as shown in the bioconductor.org website

- The argument defaults for `buildPkgDependencyDataFrame` have changed
  to
  match that of `tools::package_dependencies` which is more appropriate
  for
  identifying "strong" dependencies (@vjcitn, #55).

- `latestPkgStats`, `activitySince`, and `pkgDownloadRank` are designed
  to
  pull package statistics from the GitHub API using `gh` packages or
  from
  Bioconductor. They are esp. useful for reporting to funding agencies.

[biocthis](/packages/biocthis)
--------

                        Changes in version 1.7.0                        

NEW FEATURES

- use_bioc_pkgdown_css(): helps you style your pkgdown website with
Bioconductor colors. See
https://github.com/lcolladotor/biocthis/issues/34 for details.
- use_bioc_badges(): helps you list all the Bioconductor badges (for
software packages). See
https://github.com/lcolladotor/biocthis/pull/35 for details.
- use_bioc_feature_request_template(): creates a feature request
template for your GitHub repository. See
https://github.com/lcolladotor/biocthis/pull/33 for details.
use_bioc_issue_template() and use_bioc_support() were also updated
to be more Bioconductor-centric.
- use_tinytest() adds support for tinytest. See
https://github.com/lcolladotor/biocthis/pull/32 for details.

BUG FIXES

- Fixed pkgdown website creation. See
https://github.com/lcolladotor/biocthis/issues/29 for details. This
is also tangentially related to
https://github.com/lcolladotor/biocthis/issues/31.
- Fixed usage of r-lib/actions. We now use v2. See
https://github.com/lcolladotor/biocthis/issues/36,
https://github.com/lcolladotor/biocthis/pull/37, and
https://github.com/r-lib/actions/issues/639 for more.
- Thanks to everyone who contributed pull requests and GitHub issues!

[biocViews](/packages/biocViews)
---------

                       Changes in version 1.65.0                        

ENHANCEMENT

- (1.65.4) Change to allow Mac ARM64 binaries

- (1.65.2) Convert Sweave vignettes to Rmd

- (1.65.2) Add .gitignore

[BioNAR](/packages/BioNAR)
------

                      Changes in version 0.99.001                       

- Added a NEWS.md file to track changes to the package.

[biosigner](/packages/biosigner)
---------

                       Changes in version 1.25.6                        

MINOR MODIFICATION

- minor modification update

                       Changes in version 1.25.4                        

MINOR MODIFICATION

- biosigner.Rmd vignette: update

                       Changes in version 1.25.2                        

NEW FEATURE

- biosign: now handles SummarizedExperiment and MultiAssayExperiment
  instances

[BridgeDbR](/packages/BridgeDbR)
---------

                        Changes in version 2.7.5                        

BUG FIXES

- Fixed a method definition

                        Changes in version 2.7.4                        

NEW FEATURES

- improved code

                        Changes in version 2.7.3                        

NEW FEATURES

- Updated to BridgeDb 3.0.14

- Support primary/secondary identifier annotation

- Added mapping Bioregistry.io compact identifiers

- Added getting a DataSource by their Bioregistry.io prefix

BUG FIXES

- Fixed some metadata in the package Rmd

[cageminer](/packages/cageminer)
---------

                        Changes in version 1.3.3                        

NEW FEATURES

- Added CITATION file linking to published paper.
- Increased unit test coverage to 100%.
- Changed pkgdown website template to match new standard of all of my
packages.


[Cardinal](/packages/Cardinal)
--------

                 Changes in version 2.99.1 (2022-10-31)                 

SIGNIFICANT USER-VISIBLE CHANGES

- Updated vignettes for Cardinal 3

- Widened default m/z 'tolerance' for sparse spectra

- Switched to linear interpolation for sparse spectra

                 Changes in version 2.99.0 (2022-10-26)                 

SIGNIFICANT USER-VISIBLE CHANGES

- Updated out-of-memory backend to Matter 2.0

- Removed support for legacy classes and methods

[cbaf](/packages/cbaf)
----

                 Changes in version 1.20.0 (2022-10-24)                 

New Features

- Introducing oneOfEach! using this option would force the
  availableData() to return one study for each cancer type that
  contain a technique of interest, instead of storing all the data in
  an excel
  file.

- It is now possible to generate Heatmap for a single gene.

[cBioPortalData](/packages/cBioPortalData)
--------------

                       Changes in version 2.10.0                        

New features

- Add cgdsr to cBioPortalData migration vignette (@kmezhoud, #54)
- Unmapped experiments are now added to the metadata in cBioDataPack
- Set default api. = api/v2/api-docs in cBioPortal to access the API
protocol's new location
- The fetchData developer function added to handle both molecularData
and mutationData requests as deduced from the molecularProfileIds
vector

Bug fixes and minor improvements

- Check for valid studyIds with getStudies in cBioDataPack
- ask argument correctly passed down to caching mechanism in
cBioDataPack
- check_build option available in cBioDataPack particularly for new
studies that have not been checked against.

[celda](/packages/celda)
-----

                 Changes in version 1.13.0 (2022-10-20)                 

- Bug fixes related to cluster labels stored as factors and plotting
- Updated sparse matrix conversion to work with Matrix v1.4-2

[CellBarcode](/packages/CellBarcode)
-----------

                        Changes in version 1.2.2                        

Enhancements

- Reduce memory usage by removing the original sequence in
barcodeObj.

- Impromve fastq reading efficiency by using Rcpp.

New features

- bc_splitVDJ split VDJ barcode (experimental).

- bc_extract_10X_scSeq function used to extract transcripted barcode
in 10X Genomics scRNASeq data.

Bug fix

- fix the error of bc_cure_cluster.

[cellxgenedp](/packages/cellxgenedp)
-----------

                         Changes in version 1.2                         

SIGNIFICANT USER-VISIBLE CHANGES

- (v. 1.1.4) allow custom files_download() cache. Thanks
@stemangiola,
https://github.com/mtmorgan/cellxgenedp/pull/9

- (v. 1.1.6) datasets ethnicity field renamed to
self_reported_ethnicity

- (v. 1.1.7) use zellkonverter's basilisk-based Python parser to read
H5AD files in the vignette, see
https://github.com/theislab/zellkonverter/issues/78

OTHER

- (v. 1.1.2) reset cache on build machines weekly

- (v. 1.1.6) use {rjsoncons} CRAN package for queries, rather than
local implementation. Thanks @LiNk-NY,
https://github.com/mtmorgan/cellxgenedp/pull/12

[ChIPanalyser](/packages/ChIPanalyser)
------------

                       Changes in version 1.19.0                        

Package Updates

- New model - integrating chromatin affinity and chromatin states

- Genetic algorithm to infer optimal paramters

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.31.4                        

- update the oligoFrequency function to remove the limitation of
  sequence
  length

                       Changes in version 3.31.3                        

- support `CSV` format for `toGRanges`

                       Changes in version 3.31.2                        

- use 'entrez_id' instead of 'EntrezID' for getEnrichedPATH() output

                       Changes in version 3.31.1                        

- add parameter check of proximal.promoter.cutoff and
  immediate.downstream.cutoff for assignChromosomeRegion.

- fix the bug in `plotBinOverRegions`

[ChIPQC](/packages/ChIPQC)
------

                       Changes in version 1.33.1                        

- Make use of BiocParallel more robust - now a Depends

- Bugfix for(is.na(peaks)) issue

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.33.4                        

- add citation Q. Wang (2022) (2022-10-29, Sat)

                       Changes in version 1.33.3                        

- allows passing user defined color to vennpie() (2022-10-20, Thu,
#202, #207)
- add columns paramter to annotatePeak() to better support passing
EnsDb to annoDb (#193, #205)
- export getAnnoStat() (#200, #204)

                       Changes in version 1.33.2                        

- supports by = "ggVennDiagram" in vennplot function (2022-09-13,
Tue)

                       Changes in version 1.33.1                        

- plotPeakProf() allows passing GRanges object or a list of GRanges
objects to TxDb parameter (2022-06-04, Sat)
- add test files for getTagMatrix() and plotTagMatrix()
- getBioRegion() supports UTR regions (3'UTR + 5'UTR)
- makeBioRegionFromGranges() supports generating windoes from
self-made GRanges object
- allow specify colors in covplot() (2022-05-09, Mon, #185, #188)

[CHRONOS](/packages/CHRONOS)
-------

                       Changes in version 1.25.1                        

- Updated broken link in downloadKEGGPathwayList().

[CircSeqAlignTk](/packages/CircSeqAlignTk)
--------------

                 Changes in version 0.99.0 (2022-08-09)                 

- Released CircSeqAlignTk

[ClassifyR](/packages/ClassifyR)
---------

                        Changes in version 3.2.0                        

- 
  Fast Cox survival analysis.

- 
  Simple parameter sets, as used by crossVaildate, now come with
  tuning parameter grid as standard.

- 
  Wrappers are greatly simplified. Now, there is only one method
  for a data frame and they are not exported because they are not
  used directly by the end-user anyway.

- 
  prepareData function to filter and subset input data using
  common ways, such as missingness and variability.

- 
  Invalid column names of data (e.g. spaces, hyphens) are
  converted into safe names before modelling but converted back
  into original names for tracking ranked and selected features.

- 
  available function shows the keywords corresponding to
  transformation, selection, classifier functions.

- 
  More functions have automatically-selected parameters based on
  input data, reducing required user-specified parameters.

- 
  New classifiers added for random survival forests and extreme
  gradient boosting.

- 
  Adaptive sampling for modelling with uncertainty of class
  labels can be enabled with adaptiveResamplingDelta.

- 
  Parameter tuning fixed to only use samples from the training
  set.

[cleaver](/packages/cleaver)
-------

                 Changes in version 1.35.1 (2022-08-31)                 

- Fix: adapt vignette code to new UniProt.ws.

[clusterExperiment](/packages/clusterExperiment)
-----------------

                       Changes in version 2.17.3                        

- Fixed error in `plotBarplot` in handling non-assigned (-1) when
  comparing two clusterings (colors to clusters would be incorrect)

- Fixed error in subsetting of `ClusterExperiment` object if the
  unassigned cluster category (-1) was given a label in `clusterLegend`
  different than "-1" (would make those samples regular cluster).

[clusterProfiler](/packages/clusterProfiler)
---------------

                        Changes in version 4.5.3                        

- GSEA() supports GSONList object (2022-09-21, Wed)
- enricher() supports GSONList object (2022-09-06, Tue)

                        Changes in version 4.5.2                        

- support passing a GSON object to enricher(USER_DATA) and
GSEA(USER_DATA) (2022-8-01, Mon)
- gson_kegg_mapper() allows building a gson object from outputs of
KEGG Mapper service (2022-07-29, Fri, #492)
- fix show method for compareClusterResult (2022-06-21, Tue, #473)
- gson_KEGG() download latest KEGG and output a GSON object
(2022-06-08, Wed)
- support passing a GSON object to gseKEGG(organism)
- support passing a GSON object to enrichKEGG(organism) (2022-06-06,
Mon)

                        Changes in version 4.5.1                        

- follow KEGG api upgrade that change from http to https (2022-06-06,
Mon)
- use 'wininet' to download KEGG data when .Platform$OS.type =
"windows" (2022-06-03, Fri)
- mv read.gmt and read.gmt.wp to the 'gson' package and reexport
these
two functions from 'gson' (2022-04-28, Thu)
- fix compareCluster when fun = enrichPathway(2022-4-28, Thu)

[CNVfilteR](/packages/CNVfilteR)
---------

                       Changes in version 1.11.1                        

BUG FIXES

- Minor bug fixed: now loadVCFs() properly overwrites vcf and generates
  tabix files if temp dir is used. It prevents an unlikely but unwanted
  situation: when, in a same R session, a user updates VCF content
  without renaming it. In this case, file already existed and the new
  VCF/GZ file was not being used by CNVfilteR.

[cogeqc](/packages/cogeqc)
------

                        Changes in version 1.1.8                        

CHANGES

- Synteny assessment formula now also considers scale-free topology
fit.

BUG FIXES

- Reference-based orthogroup inference does not require the exact
same
set of species anymore.

                        Changes in version 1.1.7                        

BUG FIXES

- Variable Duplications_50 of the duplications data frame was not
matching variable Dups of the stats data frame in the output of
read_orthofinder_stats()

NEW FEATURES

- Replaced dispersal formula with a more meaningful and interpretable
one.
- Added a max_size param to plot_og_sizes() to ignore OGs larger than
a specific size.

                        Changes in version 1.1.4                        

NEW FEATURES

- Added option to scale scores by the maximum value

                        Changes in version 1.1.3                        

NEW FEATURES

- Added a correction for overclustering in calculate_H() that
penalizes protein domains in multiple orthogroups.
- Updated vignette to provide a detailed explanation of how
homogeneity scores are calculated.

[coMET](/packages/coMET)
-----

                 Changes in version 1.29.2 (2022-10-26)                 

- Remove dependancy colortools

- Add functions to replace the package colortools

[compEpiTools](/packages/compEpiTools)
------------

                       Changes in version 1.31.1                        

- The following function is updated
  + GRbaseCoverage/GRcoverage/GRenrichment: The maxDepth setting of
  ApplyPileupsParam was raised to 1000000

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

                       Changes in version 2.13.4                        

- `anno_barplot()`: fixed a bug when split is set, the bars are wrongly
  plotted under besides = TRUE.

- `anno_boxplot()`: add two new argumetn: `add_points` and `pt_gp`.

- fixed a bug of size of column title wrongly calculated.

                       Changes in version 2.13.2                        

- `HeatmapAnnotation()`: fixed a bug where annotation legends are not
  all generated when `df` is set.

- `UpSet()`: now `bg_col` can be a vector of length more than two.

- `oncoPrint()`: Add `pct_include` argument.

- `anno_density()`: fixed a bug where `xlim` is ignored for "heatmap".

                       Changes in version 2.13.1                        

- `column_title_rot` can be set with any degree value.

- automatically recognize Jupyter environment.

- `UpSet()`: `comb_col` now is correctly assigned when the combination
  matrix is transposed.

[CompoundDb](/packages/CompoundDb)
----------

                         Changes in version 1.1                         

Changes in version 1.1.6

- CompDb tests also for NA input.

Changes in version 1.1.5

- MsBackendCompDb always returns collisionEnergy as numeric.

Changes in version 1.1.4

- Add script to create a CompDb from a MassBank database.

Changes in version 1.1.3

- Expand vignette with examples to create CompDb databases from
scratch.

Changes in version 1.1.2

- Add insertCompound and deleteCompound functions to add or remove
compounds from a CompDb or IonDb.

Changes in version 1.1.1

- Fix wrong warning message in deleteIon.
- Change database data type for internal ion_id from character to
integer.

[CoreGx](/packages/CoreGx)
------

                        Changes in version 2.1.7                        

- Fixed a bug when deleting a TreatmentResponseExperiment assay via
NULL assignment
- Added a names S4 method fo TreatmentResponseExperiment to enable
tab
autocomplete with $ access in interative sessions
- Added additional methods for drug combination modelling; these
changes will be documented in a new vignette once we have thoroughly
tested the new functions
- For now these are experimental and should not be considered a
stable API

                        Changes in version 2.1.6                        

- Changed default parallelization strategy inside aggregate2 (and
therefore inside aggregate,TreatmentResponseExperiment-method and
endoaggregate) to split the table into nthread tables instead of
using by
- Result should be (1) parallelization is now always faster than
serial computations, which was not true previously
- Memory usage of parallelization should be much smaller, since we
aren't splitting into a very long list of tables
- Optimized the internal representation of the
TreatmentResponseExperiment assay index to remove storage of NA for
rowKey-colKey combinations with no observations in any assay
- This was causing memory usage to baloon if both rowKey and
colKey were a large sequence
- Prepended a "." to the internal assay index column names for each
assay
- This should reduce name clashes between internal TRE metadata
and the column names of an assay (specifically, you can now have
a column with the same name as the assay)

                        Changes in version 2.1.5                        

- Add error message to CoreSet,show-method which lets users know to
use updateObject if the slot names are not valid

                        Changes in version 2.1.4                        

- Add endoaggregate method to compute TreatmentResponseExperiment
assay aggregations within the object
- Add mergeAssays method to allow joining assaying within a
TreatmentResponseExperiment

                        Changes in version 2.1.3                        

- Updated CoreSet vignette to reflect recent changes to the object
structure
- Renamed the LongTable vignette to TreatmentResponseExperiment and
updated the content to reflect the changes in class structure
from 2.1.1
- Generated new TreatmentResponseExperiment class and structure
diagrams and inlcuded them in the TreamentResponseExperiment
vignette
- Added new example TreatmentResponseExperiment object to package
data
- Added various unit tests for the LongTable (and therefore also the
TreatmentResponseExperiment)
- Added proper documentation object for the TREDataMapper-accessors
- Added aggregate methods for data.table and LongTable
- Added endoaggregate method for LongTable, which uses aggregate
internally but assigns the result back to the object and returns the
updated object. Thus this method is an endomorphic version of
aggregate.
- Added new argument summarize to assay,LongTable-method which only
attaches columns which have been summarized over if FALSE
- Added assayCols and assayKeys helper methods to retrive valid assay
column names or the key columns for an assay, respectively.

                        Changes in version 2.1.2                        

- Fix bug in logLogisticRegression causing tests to fail in
Bioconductor 3.16 daily builds

                        Changes in version 2.1.1                        

- First update since Bioconductor 3.15 release
- Merged rework of the LongTable class back into main branch
- The object has now been updated to

[cTRAP](/packages/cTRAP)
-----

                       Changes in version 1.14.1                        

- Add cTRAP version to welcome modal
- Fix crash when submitting Celery job (included missing floweRy
function)
- Fix large JSON responses from DT blocked by reverse proxy default
settings

[cytomapper](/packages/cytomapper)
----------

                 Changes in version 1.9.2 (2022-10-24)                  

- measureObjects accounts for single objects in images

                 Changes in version 1.9.1 (2022-07-24)                  

- measureObjects function transfers mcols of images and masks to
  colData(sce)

- allow returned object to be of SpatialExperiment

[CytoML](/packages/CytoML)
------

                        Changes in version 3.11                         

API Changes

- Rename argument sampNLoc -> sample_names_from in open_flowjo_xml
- All parsers (flowjo/cytobank/diva_to_gatingset) now return
GatingSet
based on cytoset rather than ncdfFlowSet
- Add trans argument to cytobank_to_gatingset to allow overriding of
transformations from gatingML file (#76)
- gatingset_to_flowjo now uses a docker image with a compiled
converter: hub.docker.com/r/wjiang2/gs-to-flowjo
- Some updates to how flowjo_to_gatingset searches for FCS files
(#77)
- Add include_empty_tree option to flowjo_to_gatingset to include
samples without gates
- Allow gatingset_to_flowjo to take a path to a GatingSet archive
directory
- Add gating_graphGML to replace gating.graphGML method for
openCyto::gating generic
- Filter samples by panel when parsing cytobank experiment and add
ce_get_samples, ce_get_panels

Fixes/internal changes

- Automatic time scaling of samples from FlowJo workspaces now
handled
by flowjo_to_gatingset RGLab/cytolib#33
- Handle change to default stringsAsFactors=FALSE in R 4.0
- Eliminated extra intermediate files left in temp directory during
workspace parsing
- Switch usage of GatingSetList to merge_gs_list
- Solve some Windows build issues
- Switch from experimental::filesystem to boost::filesystem in C++
FlowJo parser
- Add CytoML XSD to installation

                        Changes in version 3.10                         

API Changes

- Change handling of quad gates according to RGLab/cytolib#16

- Renaming of methods:

- openWorkspace -> open_diva_xml, open_flowjo_xml
- cytobankExperiment -> open_cytobank_experiment
- cytobank2GatingSet -> cytobank_to_gatingset
- parseWorkspace -> flowjo_to_gatingset, diva_to_gatingset
- getSampleGroups -> fj_ws_get_sample_groups,
diva_get_sample_groups
- getSamples -> fj_ws_get_samples, diva_get_samples
- getKeywords -> fj_ws_get_keywords
- getCompensationMatrices -> ce_get_compensations
- getTransformation -> ce_get_transformations
- compare.counts -> gs_compare_cytobank_counts

- Renaming of classes:

- divaWorkspace -> diva_workspace
- flowJoWorkspace -> flowjo_workspace

- Add CytoML.par.set, CytoML.par.get for setting parameters in CytoML
namespace

Fixes/internal changes

- Make gatingset_to_cytobank export cytobank ML with attribute
namespaces
- Allow diva_to_gatingset to use compensation matrix from xml
- Pass ... args from cytobank_to_gatingset appropriately down to FCS
parser
- Fix some issues with scaling of gates parsed from Diva workspace
(#64)
- Guard against unsupported transformations being added to GatingSet
during Diva parsing
- Switch diva_to_gatingset to using flowjo_log_trans instead of
logtGml2_trans
- Fix ported flowUtils::xmlTag to enable self-closing tags
- Make gating.graphGML lookup tailored gates by FCS name as well as
file id
- Add some flexibility to getSpilloverMat used in gatingset_to_flowjo

[dagLogo](/packages/dagLogo)
-------

                       Changes in version 1.35.1                        

- Update uniprot rest api.

[dcanr](/packages/dcanr)
-----

                 Changes in version 1.14.0 (2022-10-11)                 

- Removed implementation for the ECF method as the package for
computing the statistic, COSINE, is now deprecated
- Added a method to run the z-score method with greater flexibility

[dearseq](/packages/dearseq)
-------

                 Changes in version 1.9.4 (2022-07-20)                  

- use scattermore in plot_weights()

                 Changes in version 1.8.3 (2022-07-13)                  

- fixed a bug in NA handling with which_weights != "none"

                 Changes in version 1.8.1 (2022-04-28)                  

- fixed a small bug in permutation p-values introduced with the 1.8.0
release

[decompTumor2Sig](/packages/decompTumor2Sig)
---------------

                 Changes in version 2.13.1 (2022-05-09)                 

- Fix: corrected recognition of EstimatedParameters (from package
  pmsignature)
  objects obtained with the use of a background signature.

- Fix: corrected internal scripts for production of external example
  data.

[deepSNV](/packages/deepSNV)
-------

                 Changes in version 1.99.3 (2013-07-25)                 

Updates

- A few changes to shearwater vignette

- Renamed arguments pi.gene and pi.backgr in makePrior()

Bugfixes

- Fixed bug in bf2Vcf() when no variant is called

                 Changes in version 1.99.2 (2013-07-11)                 

Updates

- Updated CITATION

- Added verbose option to bam2R to suppress output

- Changed mode() to "integer" for value of loadAllData()

Bugfixes

- Fixed bug when only one variant is called in bf2Vcf()

                 Changes in version 1.99.1 (2013-06-25)                 

Updates

- Using knitr for prettier vignettes

- Including shearwater vignette

Bugfixes

- fixed issues with deletions in bf2Vcf()

- makePrior() adds background on all sites

                 Changes in version 1.99.0 (2013-04-30)                 

Updates

- New shearwater algorithm

- Including VCF output through summary(deepSNV, value="VCF")

                 Changes in version 1.42.1 (2022-05-31)                 

Bugfixes

- Fully adjust VCF creation to use VariantAnnotation 1.27.6 and later's
  API
  (fixes broken VCF output like
  `##fileformat=<ID=fileformat,Value=VCFv4.1>`)

- Fix crash when encountering reads that don't have an NM field

- Fix crash when compiled with CentOS 8's GCC 8.3.1-4

- Correct the template VCF files' `#CHROM` header lines

[DegNorm](/packages/DegNorm)
-------

                        Changes in version 1.7.1                        

- .gtf_parse function has been updated so overlapping exons are merged
  in total transcript

[DEGreport](/packages/DEGreport)
---------

                       Changes in version 1.33.2                        

- Remove lasso2 from dependencies

                       Changes in version 1.33.1                        

- Remove Nozzle.R1 from dependencies

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.24.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Move the aperm() S4 generic to BiocGenerics.

[DExMA](/packages/DExMA)
-----

                        Changes in version 1.5.1                        

- Improvement of some functions

- Remove seeCov function

- ACAT p-values combination method added

- Imputation indicators added to the missGenesImput fucntion

[DiffBind](/packages/DiffBind)
--------

                         Changes in version 3.8                         

- Rollover bugfixes from Release branch

- Re-normalize after applying black/greylists

- Apply Rhtslib patch

- Update vignette

- Remove contact info for Gord

[DifferentialRegulation](/packages/DifferentialRegulation)
----------------------

                        Changes in version 1.0.7                        

- bug (prompting when analyzing 1 cluster only) in
  DifferentialRegulation function fixed

                        Changes in version 1.0.6                        

- bug giving error (in plot_pi function) fixed

                        Changes in version 1.0.5                        

- information about convergence added (in DifferentialRegulation
  function)

                        Changes in version 1.0.4                        

- bug giving error (in DifferentialRegulation function) fixed

                        Changes in version 1.0.3                        

- vignettes extended

                        Changes in version 1.0.2                        

- function 'plot_pi' added

                        Changes in version 1.0.1                        

- memory usage (of DifferentialRegulation function) reduced

[Dino](/packages/Dino)
----

                        Changes in version 1.3.1                        

- Bug fix for input data formats: matrices with named rownames and
DelayedMatrix
- Vignette updates
- Citation updates
- Contact updates

[dittoSeq](/packages/dittoSeq)
--------

                        Changes in version 1.10                         

- Added ability to plot multiple 'var' in a single 'dittoPlot()',
'dittoDimPlot()', 'dittoScatterPlot()', 'dittoDimHex()', and
'dittoScatterHex()' call by giving a vector of genes or continuous
metadata to the 'var' or 'color.var' input. Customization of how the
"multivar" data is displayed can be controlled with: 1-
'multivar.aes' (context: 'dittoPlot()' only) - which plot aesthetic
is utilized for displaying var-values. 2- 'multivar.split.dir' -
faceting direction to use for var-data when combining with an
additional 'split.by' variable.
- Improved the compatibility with 'split.by'/faceting customizations,
specifically with 'split.adjust = list(scales = "free")', by making
implementations of 'min'/'max' inputs less intrusive. Note: This
change very minorly alters the default output of some plotters.
- Improved error messaging for cases where 'object' does not have
cell/column names.

[DOSE](/packages/DOSE)
----

                       Changes in version 3.23.3                        

- replace DO.db to HDO.db (2022-10-7, Fri)
- add values of organism, keytype and setType for GSEA_internal()
(2022-09-21, Wed)
- add values of organism, keytype and ontology for
enricher_internal()
(2022-09-21, Wed)
- move inst/extdata/parse-obo.R to HDO.db package (2022-08-29, Mon)
- rename qvalues to qvalue in gseaResult object (2022-08-29, Mon)

                       Changes in version 3.23.2                        

- Support GSON object in GSEA_internal() (2022-06-08, Wed)

                       Changes in version 3.23.1                        

- Support GSON object in enricher_internal() (2022-06-06, Mon)

[DropletUtils](/packages/DropletUtils)
------------

                       Changes in version 1.18.0                        

- Added an intersect.genes= option to read10xCounts() for samples
  with inconsistent gene information. Automatically fix empty
  chromosome names for mitochondrial genes in certain Cellranger
  outputs.

[edgeR](/packages/edgeR)
-----

                       Changes in version 3.40.0                        

- 
  New argument 'hairpinBeforeBarcode' for processAmplicons(). The
  revised function can process reads where the
  hairpins/sgRNAs/sample index sequences are in variable
  positions within each read. When 'plotPositions=TRUE' a density
  plot of the match positions is created to allow the user to
  assess whether they occur in the expected positions.

- 
  Update C++ BLAS calls to account for USE_FC_LEN_T setting in R
  4.3.0.

- 
  Bug fix to R_compute_apl.cpp to make sure GLM working weights
  are zero when fitted mu=0.

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.17.4                        

- rename parameters of emapplot(), centplot() and treeplot()
(2022-09-11, Sun)

                       Changes in version 1.17.3                        

- align the dots in treeplot() (2022-10-1, Sat)
- fix a bug in color legend of treeplot() (2022-10-1, Sat)

                       Changes in version 1.17.2                        

- autofacet to automatically split barplot and dotplot into several
facets (2022-09-06, Tue)
- dotplot method for enrichResultList object
- add parameters hilight_category, alpha_hilight, alpha_nohilight for
cnetplot() and emapplot (2022-09-4, Sun)
- change round digits of cnetplot scatterpie legend to 1 (2022_8_29,
Mon).
- gsearank() can export result as a table when output = "table"
(2022-08-29, Mon, #184)
- fix a bug in fc_readable() (2022-08-29, Mon, #189)
- allows passing color="NES" to dotplot() for gseaResult object
(2022-08-29, Mon, #14)

                       Changes in version 1.17.1                        

- fix a bug in
https://github.com/YuLab-SMU/clusterProfiler/issues/488
(2022-08-25, Thu)
- support multiple gene sets in geom_gsea_gene layer (2022-08-25,
Thu)
- geom_gsea_gene layer (2022-08-24, Wed)
- add parameters symbol and pvalue for heatplot.enrichResult()
(2022-08-20, Sat)
- change default values of group_category and node_label in ssplot()
(2022-07-04, Mon)
- update document of ssplot() (2022-07-04, Mon)
- gseaplot() and gseaplot2() return gglist object instead of plotting
the figure (2022-05-05, Thu)
- fix ridgeplot when x@readable = TRUE (2022-04-30, Sat)

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.21.5                        

- Fix spell mistakes (thanks to Mike Love).

                       Changes in version 2.21.3                        

- Add minimal required Ensembl API version for certain data fields
  (stable_id_version). Issue
  https://github.com/jorainer/ensembldb/issues/139

                       Changes in version 2.21.1                        

- Since *is circular* flag for chromosomes extracted using the Ensembl
  Perl API
  is always `FALSE`: manually set `isCircular` to `TRUE` for
  chromosome(s)
  named `"MT"` (issue
  https://github.com/jorainer/ensembldb/issues/133).

[EpiCompare](/packages/EpiCompare)
----------

                        Changes in version 1.1.2                        

New features

- rebin_peaks:
- Added arg drop_empty_chr to automatically drop chroms that
aren't in any of the peakfiles.
- Added "score" as one of the default intensity_cols in all
relevant functions.
- Make examples use 5000bp bins to speed up.
- translate_genome:
- Add default_genome arg to handle genome=NULL.
- bpplapply:
- New exported function to automate handling of known issues with
BiocParallel across OS platforms.
- Enable users to specify their own apply function.
- get_bpparam: Add args to allow users to choose which BiocParallel
func to use.
- checkCache: Make default arg cache=BiocFileCache::BiocFileCache(ask
= FALSE) to skip user input during runtime.
- precision_recall:
- Change increment_threshold arg to n_threshold arg, using the
seq(length.out=) feature to avoid accidentally choosing an
inappropriately large increment_threshold.
- gather_files:
- Replace iterator with bpplapply.
- Pass up args from bpplapply.
- Provide warning message, not error, when 0 files found. Returns
NULL.
- Add "multiqc" as a search option.
- Add dedicated subfunctions for reading in a variety of
nf-core/cutandrun outputs files:
read_picard,read_multiqc,read_bowtie,
read_trimgalore,read_bam,read_peaks
- Add file paths to each object.
- Add new arg rbind_list.
- rebin_peaks/compute_corr: -Change defaultbin_size from 100 --> 5kb
to improve efficiency and align with other defaults of other
packages (e.g Signac).
- tss_plot:
- Pass up more arg for specifying upstream/downstream.
- EpiCompare: Pass up new args:
- bin_size
- n_threshold
- workers

Bug fixes

- Fix rebin_peaks unit tests.
- Fix pkg size issue by adding inst/report to .Rbuildignore.
- EpiCompare wasn't being run when reference was a single unlisted
GRanges object because it was indeed length>1, but the names were
all NULL. Now fixed.
- plot_precision_recall: Set default initial_threshold= to 0.
- Switch from BiocParallel to parallel, as the former is extremely
buggy and inconsistent.

                        Changes in version 1.1.1                        

New features

- check_genome_build: Add translate_genome as prestep.
- rebin_peaks:
1.  Move all steps that could be done just once (e.g. creating the
genome-wide tiles object) outside of the BiocParallel::bpmapply
iterator.
2.  Ensure all outputs of BiocParallel::bpmapply are of the same
length, within the exact same bins, so that we can return just
the bare minimum data needed to create the matrix (1 numeric
vector/sample).
3.  Instead of rbinding the results and then casting them back into
a matrix (which is safer bc it can handle vectors of different
lengths), simply cbind all vectors into one matrix directly and
name the rows using the predefined genome-wide tiles.
4.  Because we are no longer rbinding a series of very long tables,
this avoids the issue encountered here #103. This means this
function is now much more scalable to many hundreds/thousands of
samples (cells) even at very small bin sizes (e.g. 100bp).
5.  A new argument keep_chr allows users to specify whether they
want to restrict which chromosomes are used during binning. By
default, all chromosomes in the reference genome are used
(keep_chr=NULL), but specifying a subset of chromosomes (e.g.
paste0("chr",seq_len(12))) can drastically speed up compute time
and reduce memory usage. It can also be useful for removing
non-standard chromosomes (e.g. "chr21_gl383579_alt",
"chrUns...", "chrRand...").
6.  As a bonus, rebin_peaks now reports the final binned matrix
dimensions and a sparsity metric.

- compute_corr:
- Added unit tests at different bin sizes.
- Allow reference to be NULL.
- Updated README to reflect latest vesion of EpiCompare with
gather_files.

Bug fixes

- Bumped version to align with Bioc devel (currently 1.1.0).
- compute_percentiles:
- Making default initial_threshold=0, so as not to assume any
particular threshold.
- rebin_peaks:
- Addressed error that occurs when there's many samples/cells with
small bins.
- plot_precision_recall: Don't plot the reference as part of the PR
curve.

[epimutacions](/packages/epimutacions)
------------

                        Changes in version 1.1.2                        

- Get plots fixed
- Get annotation fixed
- Removed extra slashes from annotations

[EWCE](/packages/EWCE)
----

                        Changes in version 1.5.8                        

Bug fixes

- GHA fix.

                        Changes in version 1.5.7                        

Bug fixes

- orthogene dependency has been replacing user entered background
gene
list with one generated from all known genes when species across
gene lists and reference dataset are the same. This has now been
fixed.

                        Changes in version 1.5.5                        

Bug fixes

- standardise_ctd:
- Always force "specificity_quantiles" to be one of the matrices
in each level.

                        Changes in version 1.5.4                        

New features

- filter_ctd_genes
- Now exported.
- Can handle standardized CTD format.
- get_ctd_matrix_names: New function to get a list of all data
matrices in CTD.

Bug fixes

- check_ewce_genelist_inputs:
- User reported potential bug in
code:https://github.com/NathanSkene/EWCE/issues/71
- Fixed by removing conditional and instead always filtering out
genes not present in CTD/SCT.
- standardise_ctd:
- Add check_species()
- Ensure all matrices become sparse when as_sparse=TRUE.
- Generalize to matrices of any name.
- fix_celltype_names:
- Ensure all celltype names are unique after standardization.

                        Changes in version 1.5.3                        

New features

- genelistSpecies now passed to prepare_genesize_control_network in
bootstrap_enrichment_test meaning gene list species will be inferred
from user input.

                        Changes in version 1.5.2                        

New features

- drop_uninformative_genes:
- Expose new args: dge_method, dge_test, min_variance_decile
- merged_ctd: Actually merge the CTDs into one when as_SCE=FALSE.

Bug fixes

- Remove hard-coded file path separators (e.g.
sprintf("%s/MRK_List2.rpt", tempdir())) to be more compatible with
Windows.

                        Changes in version 1.5.1                        

New features

- Made substantial updates to orthogene, so going through and making
sure everything still works / is able to take advantage of new
features (e.g. separation of non121_strategy and agg_func args,
many:many mapping):
- filter_nonorthologs: Pass up args from
orthogene::convert_orthologs.
- generate_celltype_data: @inheritDotParams
- Update GHA.
- Bump to R (>= 4.2) now that we're developing on Bioc 3.16.

Bug fixes

- Avoid downloading large "MRK_List2.rpt" file any more than is
necessary for testing.

[exomePeak2](/packages/exomePeak2)
----------

                 Changes in version 1.8.1 (2022-05-16)                  

- Add parameters to initiate mode of motif based analysis.

- Add parameters to configure saving directories.

- Fix saving issue when -log pvalues have NA/NaN (no change in default
  method).

[ExploreModelMatrix](/packages/ExploreModelMatrix)
------------------

                        Changes in version 1.9.2                        

- Fix inconsistency in cooccurrence plots obtained when providing the
design matrix directly, if the sample table contained additional
columns not used for the design matrix.

[extraChIPs](/packages/extraChIPs)
----------

                 Changes in version 1.1.5 (2022-10-11)                  

- Added makeConsensus() and updated vignette

                 Changes in version 1.1.3 (2022-08-01)                  

- Added plotOverlaps() for generation of Venn Diagrams and
ComplexUpset plots

                 Changes in version 1.1.2 (2022-07-22)                  

- Added collapseTranscripts = "auto" as the default for plotHFGC()
- getProfileData() now returns log2 transformed data by default
- getProfileData() now uses bplapply() internally
- Bugfixes for plotPie(), distinctMC() colToRanges() and
stitchRanges()

                 Changes in version 1.1.1 (2022-06-01)                  

- Bugfix so as_tibble() respects original column names

[fastreeR](/packages/fastreeR)
--------

                 Changes in version 1.3.0 (2022-11-01)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_16
  branch.

                 Changes in version 1.2.0 (2022-11-01)                  

- Bump x.y.z version to even y prior to creation of RELEASE_3_16
  branch.

                 Changes in version 1.1.6 (2022-09-26)                  

- Update java backend to BioInfoJavaUtils-1.2.4 (revamped
  FastaManager).

                 Changes in version 1.1.5 (2022-06-23)                  

- Update vignette minor bug.

                 Changes in version 1.1.4 (2022-06-16)                  

- Update vignette to handle getting sample files through https
  failure.

                 Changes in version 1.1.3 (2022-06-06)                  

- Update vignette to handle getting sample files through https
  failure.

                 Changes in version 1.1.2 (2022-05-15)                  

- Update possibly corrupted samples.vcf.gz.

                 Changes in version 1.1.1 (2022-05-08)                  

- Update tests to improve coverage.

                 Changes in version 1.1.0 (2022-05-08)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_15
  branch.

[fgga](/packages/fgga)
----

                 Changes in version 1.5.0 (2022-07-13)                  

- Made the following significant changes
  o changed the names of the varianceSVM and svmSVM functions to
  varianceOnto and svmOnto respectively
  o added the PO, ZFA and HPO ontologies to the preCoreFG function
  o added functions to evaluate the performance of the classification
  process

[fgsea](/packages/fgsea)
-----

                       Changes in version 1.23.1                        

- Introduced GESECA method for multi-conditional gene set enrichment
  (see geseca-tutorial vignette for details)

- Enrichment table plots now use cowplot (issue #101)

- Enrichment table plots are more easily fine-tuned (issue #29)

[FindIT2](/packages/FindIT2)
-------

                        Changes in version 1.2.3                        

- fix a bug when Txdb have no gene scaffold

                        Changes in version 1.2.1                        

- add FindIT2 publication info so people can cite

[fishpond](/packages/fishpond)
--------

                        Changes in version 2.4.0                        

- For simple paired swish analysis, adding a fast=1 method which uses
a one-sample z-score on the paired LFCs (averaged over samples, then
median over inferential replicates). The permutation is computed by
changing the signs of the LFC matrix and recomputing z-scores.
Testing on the vignette example, but using all the transcripts, the
one-sample z-score method takes < 20 seconds while the signed rank
method takes > 200 seconds (12x speedup), while they have a high
rate of agreement on the detected set (30:1 in common vs
discordant).
- Removed the fast=0 methods that were previously implemented where
ranks could optionally be recomputed for every permutation. This was
much slower and didn't have any appreciable benefit.
- readEDS() has moved to the eds package, such that fishpond no
longer
requires Rcpp and C++ code compilation.
- Fix bug identified by GitHub user @JosephLalli, where
importAllelicCounts would find a1 and a2 strings internal to gene
IDs, instead of at the suffix.

                       Changes in version 2.3.22                        

- Fix bug identified by GitHub user @JosephLalli, where
importAllelicCounts would find a1 and a2 strings internal to gene
IDs, instead of at the suffix.

                       Changes in version 2.3.14                        

- For simple paired swish analysis, adding a fast=1 method which uses
a one-sample z-score on the paired LFCs (averaged over samples, then
median over inferential replicates). The permutation is computed by
changing the signs of the LFC matrix and recomputing z-scores.
Testing on the vignette example, but using all the transcripts, the
one-sample z-score method takes <20 seconds while the signed rank
method takes >200 seconds (12x speedup), while they have a high rate
of agreement on the detected set ( 30:1 in common vs discordant).
- Removed the fast=0 methods that were previously implemented where
ranks could optionally be recomputed for every permutation. This was
much slower and didn't have any appreciable benefit.

                        Changes in version 2.3.7                        

- readEDS() has moved to the eds package.

[FlowSOM](/packages/FlowSOM)
-------

                        Changes in version 2.5.6                        

- Reversed PlotDimred back to scattermore instead of geom_hex. Uses
  geom_point
  if scattermore is not available.

- Added support for custom colors and limits in PlotDimRed

                        Changes in version 2.5.5                        

- AddAnnotation now works correctly after using UpdateMetaclusters. The
  function
  also works slightly different. See examples.

- Fixed issue with wrong coloring in PlotDimRed

                        Changes in version 2.5.4                        

- percentage_positives in GetFeatures now works with a FlowSOM with
  updated metacluster labels

                        Changes in version 2.5.3                        

- Bugfix in GetChannels

                        Changes in version 2.5.2                        

- Bugfixes to pass BioConductor check

- Bugfixes in Plot2DScatters: when ggpointdensity is not installed,
  normal
  ggplot colors are used

- ParseNodeSize is now exported

                        Changes in version 2.5.1                        

- Version bump to get the same version as on BioConductor

- Changed yMargin parameter in PlotFileScatters to yLim for consistency

- Added possibility to use abbrevations in GetFeatures, GetCounts,
  GetPercentages,
  PlotFileScatters, PlotDimRed, ParseQuery, Plot2DScatters and
  PlotNumbers

- Moved UpdateFlowSOM to 0_FlowSOM.R

- Changed examples in PlotFlowSOM, PlotVariable and PlotStars

- Fixed issues from check: deleted examples with GetFlowJoLabels and
  used a csv
  instead, changed geom_scattermore to geom_hex for less dependencies.
  Moved packages to Suggest instead of Import.

- ParseNodeSize is now exported

[flowSpecs](/packages/flowSpecs)
---------

                  Changes in version 1.11 (2022-07-23)                  

- Updating the specMatCalc function, to allow for inclusion of samples
  with an
  internal negative control.

- Removing the dependency hexbin, as it is not in use (ggplot2 since
  1.9.3)

[fmrs](/packages/fmrs)
----

                        Changes in version 2.0.0                        

IMPROVEMENTS SINCE LAST RELEASE

- The package is rewritten using .Call function.
- The codes for Weibull distribution are improved.

BUG FIXES

- Several bugs are fixed which caused the results to be different for
the same analysis.

[FRASER](/packages/FRASER)
------

                        Changes in version 1.8.1                        

- Bugfix in merging splicing counts (#41)

[gage](/packages/gage)
----

                       Changes in version 2.47.1                        

- korg now include 8282 KEGG species or 1449 new species beyond 2020.

- updated khier to included newly added reference pathways. kegg.gsets
  can work with 477 pathways now.

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.34.0                        

UTILITIES

- update the web links

[GeneNetworkBuilder](/packages/GeneNetworkBuilder)
------------------

                       Changes in version 1.39.2                        

- Add subsetNetwork function.

                       Changes in version 1.39.1                        

- Add unrooted network.

- Add sample code for subset graph.

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 2.2.0                        

New features

- gs_heatmap gains the winsorize_threshold parameter, to control the
behavior of the geneset heatmap in presence of extreme values,
either negative or positive ones. If not specified, the heatmap is
not introducing any winsorization.
- map2color() has a behavior that better accounts for asymmetric
ranges of values. This propagates to some of the functions that use
it for mapping to colors, such as enrichment_map(), or
ggs_backbone().

Other notes

- Fixed the behavior of the reactive elements after uploading the
GeneTonicList object at runtime.
- Fixed the label namings for the gs_heatmap function
- The enhance_table() function can handle the case where a gene is in
the enrichment results table but not present in the annotation (e.g.
annotations are updated, so some correspondences might get lost). It
also presents an informative message on which genesets/genes are
potentially responsible for the behavior.
- Some additional checks are in place for controlling the cases where
the z_score of a geneset is detected as NA (e.g. because there was a
mismatch between gene names and identifiers in the annotation).

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.34.0                        

NEW FEATURES

- Add update() method for Seqinfo() objects. See '?Seqinfo'.

- Implement getChromInfoFromUCSC() "offline mode" for a selection of
  genomes. See '?getChromInfoFromUCSC'.

- Register the following NCBI assemblies:
  - a few Triticum aestivum assemblies (bread wheat)
  - a Pteropus alecto assembly
  - an Eucalyptus grandis assembly
  - a few Plasmodium falciparum assemblies (malaria parasite)
  - the Dog10K_Boxer_Tasha assembly
  - the Felis_catus_9.0 assembly
  - the UCB_Xtro_10.0 assembly

- Register the following UCSC genomes: equCab1, equCab2, equCab3,
  mpxvRivers, hs1, canFam6, felCat9, xenTro10.

- Export and document find_NCBI_assembly_ftp_dir().

DEPRECATED AND DEFUNCT

- Remove releaseName() method for GenomeDescription objects.
  The releaseName() method for GenomeDescription objects was deprecated
  in Bioconductor 3.12 and defunct in Bioconductor 3.15.
  Also move the releaseName() generic function to the BSgenome package.

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.34.0                        

BUG FIXES

- Fix buffer overflow in cigarNarrow()/cigarQNarrow()
  See commit 5d0a29b2795a95dfb006e8e54c7532b12bd5e79c.

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.50.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Use useEnsemblGenomes() instead of useMart() wherever possible.

- Improve pcoverageByTranscript() implementation. The function now uses
  a chunking strategy to handle bigger CompressedGRangesList objects
  and
  to reduce memory footprint. Note that even with this improvement, the
  function is still not very efficient.

- Rename proteinToGenome() argument 'txdb' -> 'db'.

BUG FIXES

- Change test for 'circ_seqs' argument in makeTxDbPackageFromUCSC()
  from
  an error to a warning if argument set to NULL.

[GenomicFiles](/packages/GenomicFiles)
------------

                        Changes in version 1.34                         

BUG FIXES

- (v 1.33.1) reduceRanges() supports `files` as character vectors.

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.50.0                        

- No changes in this version.

[GenomicScores](/packages/GenomicScores)
-------------

                       Changes in version 2.10.0                        

USER VISIBLE CHANGES

- Added two new score sets and corresponding metadata for mm39:
  phastCons35way.UCSC.mm39 and phyloP35way.UCSC.mm39.

- Added new method 'wgscores()' to find which genomic scores are
  present in a given set of genomic ranges. Feature request by
  https://github.com/rcastelo/GenomicScores/issues/20

- Improved error messages about optional arguments without name or with
  a mispecified name in the call to 'gscores()'.

BUG FIXES

- Bugfix for parameters 'ref' and 'alt' when used with HDF5 backends.

[GeomxTools](/packages/GeomxTools)
----------

                        Changes in version 3.1.1                        

- Documentation correction

                        Changes in version 3.1.0                        

- No changes from 2.99.2

[ggmsa](/packages/ggmsa)
-----

                        Changes in version 1.3.3                        

- calling \dontrun{} for examples on ggmsa()

                        Changes in version 1.3.2                        

- bugfix: geom_msaBar conservation layer incorrectly aligned
issues#34(2022-5-13, Fri)

                        Changes in version 1.3.1                        

- A new feature--selects ancestral sequence on Tree-MSA plot
treeMSA_plot (2022-4-14, Thu)
- A new feature--visualization of genome alignment ggmaf (2022-4-14,
Thu)
- A test feature--visualization protein-protein interactive
(2022-4-14, Thu)
- updated the way smooth is invoked on simplot(2022-01-03, Mon)

[ggspavis](/packages/ggspavis)
--------

                 Changes in version 1.3.1 (2022-10-11)                  

- support for multiple-sample datasets

[ggtree](/packages/ggtree)
------

                        Changes in version 3.6.0                        

- Bioconductor RELEASE_3_16 (2022-11-02, Wed)

                        Changes in version 3.5.3                        

- add new citation (the iMeta 2022 paper) (2022-09-26, Mon)
- move scale_color_subtree() to the 'ggtreeDendro' package
(2022-09-23, Fri)
- update fortify method for pvclust object (2022-08-15, Mon)
- add citation of the tree data book (2022-08-13, Sat)

                        Changes in version 3.5.2                        

- scale_color_subtree() now supports passing a numeric value and
internally it will call cutree(tree, k) (2022-08-11, Thu)
- support 'linkage' class defined in the 'mdendro' package
(2022-08-11, Thu)
- clone the plot environment before assigning layout (2022-07-19,
Tue,
#516)
- bug fixed in 'equal_angle' layout (2022-07-08, Fri, #514)
- optimize geom_tiplab to better compatible with dendrogram layout
(2022-06-23, Thu, #508)

                        Changes in version 3.5.1                        

- as.phylo.hclust2 to correct edge length as displayed in
stats:::plot.hclust (2022-06-21, Tue)
- add outline to nodepies (2022-06-20, Mon, #506)
- new 'slanted' layout for branch.length = 'none' (2022-04-29, Fri,
#497)
- only works for Cartesian coordination, that means it will not
work for layout = 'radial'

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                        Changes in version 1.7.1                        

- fix a bug of subset data (when the value mapped to x aesthetic
contains negative value.)
- fix a bug when x of mapping in data is only one unique.

- https://github.com/YuLab-SMU/ggtreeExtra/issues/24

[glmGamPoi](/packages/glmGamPoi)
---------

                         Changes in version 1.9                         

- Breaking change to the way that non-standard evaluation parameters
  are handled.
  Variables in arguments such as 'pseudobulk_by' or 'subset_to' which
  evaluate
  to a single string are no longer interpreted as referring to a
  column.
  This change makes the handling of NSE more consistent.

- Add new function 'pseudobulk_sce' to easily form pseudobulk samples

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.23.1                        

- Replacing DO.db with HDO.db (2022-07-29, Mon)

[GRaNIE](/packages/GRaNIE)
------

                 Changes in version 1.1.21 (2022-12-13)                 

Major changes

- major object changes and optimizations, particularly related to
storing the count matrices in an optimized and simpler format. In
short, the count matrices are now stored either as normal or sparse
matrices, depending on the amount of zeros present. In addition,
only the counts after normalization are saved, the raw counts before
applying normalization are not stored anymore. If no normalization
is wished by the user, as before, the "normalized" counts are equal
to the raw counts. GRaNIE is now more readily applicable for larger
analyses and single-cell analysis even though we just started
actively optimizing for it, so we cannot yet recommend applying our
framework in a single-cell manner. Older GRN objects are
automatically changed internally when executing the major functions
upon the first invocation.
- various Documentation and R help updates
- the function generateStatsSummary now doesnt alter the stored
filtered connections in the object anymore. This makes its usage
more intuitive and it can be used anywhere in the workflow.
- removed redundant biomaRt calls in the code. This saves time and
makes the code less vulnerable to timeout issues caused by remote
services
- due to the changes described above, the function plotPCA_all now
can
only plot the normalized counts and not the raw counts anymore
(except when no normalization is wanted)
- the GO enrichments are now also storing, for each GO term, the
ENSEMBL IDs of the genes that were found in the foreground. This
facilitates further exploration of the enrichment results.

Minor changes

- many small changes in the code

                 Changes in version 1.1.13 (2022-09-13)                 

Major changes

- many Documentation and R help updates, the Package Details Vignette
is online
- The workflow vignette is now improved: better figure resolution,
figure aspect ratios are optimized, and a few other changes
- the eGRN graph structure as built by build_eGRN_graph() in the
GRaNIE object is now reset whenever the function
filterGRNAndConnectGenes() is successfully executed to make sure
that enrichment functions etc are not using an outdated graph
structure.
- the landing page of the website has been extended and overhauled
- removed some dependency packages and moved others into Suggests to
lower the installation burden of the package. In addition, removed
topGO from the Depends section (now in Suggests) and removed
tidyverse altogether (before in Depends). Detailed explanations when
and how the packages listed under Suggests are needed can now be
found in the new Package Details Vignette and are clearly given to
the user when executing the respective functions
- major updates to the function getGRNConnections, which now has more
arguments allowing a more fine-tuned and rich retrieval of eGRN
connections, features and feature metadata
- a new function add_featureVariation to quantify and interpret
multiple sources of biological and technical variation for features
(TFs, peaks, and genes) in a GRN object, see the R help for more
information
- filterGRNAndConnectGenes now doesnt include feature metadata
columns
to save space in the result data frame that is created. The help has
been updated to make clear that getGRNConnections includes these
features now.

Minor changes

- small changes in the GRN object structure, moved
GRN@data$TFs@translationTable to GRN@annotation@TFs. All exported
functions run automatically a small helper function to make this
change for any GRN object automatically to adapt to the new
structure
- many small changes in the code, updated argument checking, and
preparing rigorous unit test inclusion
- internally renaming the (recently changed / renamed) gene type
lncRNA from biomaRt to lincRNA to be compatible with older versions
of GRaNIE

                  Changes in version 1.1 (2022-05-31)                   

Minor changes

- added the argument maxWidth_nchar_plot to all functions that plot
enrichments, and changed the default from 100 to 50.

Bug fixes

- fixed a small bug that resulted in the enrichment plots to ignore
the value of maxWidth_nchar_plot

[graphite](/packages/graphite)
--------

                 Changes in version 1.43.2 (2022-10-31)                 

- Updated all pathway data.

[GreyListChIP](/packages/GreyListChIP)
------------

                       Changes in version 1.29.1                        

- Change maintainer email from Gord to Rory.

[GSVA](/packages/GSVA)
----

                        Changes in version 1.46                         

BUG FIXES

- Bugfix for https://github.com/rcastelo/GSVA/issues/61 to enable using
  the ssgsea method with one single column (sample) in the input data
  container.

- Bugfix when input is a SummarizedExperiment and assays contain a
  data.frame instead of a matrix.

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.26.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Try harder to find and load the matrix rownames of a 10x Genomics
  dataset.
  See commit abafbb9e99ad54a64e5013305486b97daa9442bc.

BUG FIXES

- Handle HDF5 sparse matrices where shape is not an integer vector.
  When the shape returned by internal helper .read_h5sparse_dim() is a
  double vector it is now coerced to an integer vector. Integer
  overflows
  resulting from this coercion trigger an error with an informative
  error
  message.
  See GitHub issue #48.

[hermes](/packages/hermes)
------

                        Changes in version 1.0.1                        

- Additional version argument for connect_biomart to specify an
Ensembl version.
- Fixed tests.

[Herper](/packages/Herper)
------

                        Changes in version 1.6.1                        

- Added vignette describing how to work with Herper on machines with
  old GCC compilers

- Added additonal messaging and control of messaging

- Fixed bug that causes error on new miniconda install when path is
  abbreviated with '~'

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.34.0                        

- fix the compiler issue on Mac M1 chip

- fix the C stack error in RStudio

[HiCDCPlus](/packages/HiCDCPlus)
---------

                 Changes in version 1.5.2 (2022-06-06)                  

- Fixed bug that prevented juicer to setup

                 Changes in version 1.5.1 (2022-05-08)                  

- Updated references to juicer

[hipathia](/packages/hipathia)
--------

                 Changes in version 2.13.1 (2022-07-27)                 

- Fixing bug in nodes DE with limma with high rates of 0-variance
  genes.

[HPAanalyze](/packages/HPAanalyze)
----------

                        Changes in version 1.15                         

- Changes in version 1.15.1
  + Update built-in data to HPA version 21.1.

- Changes in version 1.15.0
  + Starting devel for Bioconductor 3.16

[IgGeneUsage](/packages/IgGeneUsage)
-----------

             Changes in version 1.11.2-1.11.11 (2022-10-24)             

- man files updated

- vignette updated

- d_zibb data and inst/scripts/d_zibb.R added: simulated data from
  zero-inflated beta-binomial (ZIBB) distribution

- minor fixes, function clarity/readibility improved

- probe compilation errors on Windows

- stanmodels used properly in v. 1.11.10, no extra compilation

- stanmodels in tests/ fixed, no extra compilation

[illuminaio](/packages/illuminaio)
----------

                 Changes in version 0.41.0 (2022-11-01)                 

- The version number was bumped for the Bioconductor devel version,
which is now Bioconductor 3.17 for R-devel.

                 Changes in version 0.40.0 (2022-11-01)                 

- The version number was bumped for the Bioconductor release version,
which is now Bioconductor 3.16 for R (>= 4.2.2).

                 Changes in version 0.39.2 (2022-10-31)                 

Miscellaneous

- Fix NEWS.md version formatting.

                 Changes in version 0.39.1 (2022-09-23)                 

Bug Fixes

- readIDAT() would only return the first five Unknown.N fields.
Additional ones would be dropped.

- readIDAT() could produce Warning message: In readChar(con, nchars =
n) : truncating string with embedded nuls if the IDAT file had an
Unknown.6 field. Until we know what that field represents, it is
parsed as an (nbytes, <byte sequence>) integer vector.

                 Changes in version 0.39.0 (2022-04-26)                 

- The version number was bumped for the Bioconductor devel version,
which is now Bioconductor 3.16 for R-devel.

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.3.10 (2022-10-22)                 

- fix axis.ratio to 1 in plotSpatial function and set scales = "fixed"

                 Changes in version 1.3.9 (2022-10-14)                  

- new function for spatial community detection

                 Changes in version 1.3.8 (2022-10-13)                  

- exclude a single validity check to to class changes in the ggraph
  package

                 Changes in version 1.3.7 (2022-08-11)                  

- new function to measure minimal distances to cells of interest

                 Changes in version 1.3.6 (2022-07-28)                  

- new function to measure size of patch

                 Changes in version 1.3.5 (2022-06-15)                  

- Added three spatial context functions

                 Changes in version 1.3.4 (2022-06-14)                  

- Avoid duplicated columns when generating nodes for plotSpatial

                 Changes in version 1.3.3 (2022-06-03)                  

- Bug fix: specify the number of nodes in patchDetection function

                 Changes in version 1.3.2 (2022-05-30)                  

- Bug fix: avoid duplication of metadata entries

                 Changes in version 1.3.1 (2022-04-28)                  

- Bug fix: avoid using factors in testInteraction function for group_by

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.29.2                        

- CHANGES
  * code cleaning

                       Changes in version 1.29.1                        

- CHANGES
  * minor code improvements

[infercnv](/packages/infercnv)
--------

                 Changes in version 1.13.1 (2022-10-17)                 

- Updated how the Leiden subclustering is run to use igraph instead of
  the Python implementation called through reticulate, and to run a PCA
  before building neighboring graph.

- Added option to run a second round of subclustering that is per
  chromosome (for the HMM predictions).

- Added method that runs the i6 HMM predictions based on per chromosome
  subclustering.

- Add arguments/options for Leiden settings and change defaults.

- Change default analysis mode to subclusters.

- Fix plotting when no hclust is already stored in the object and have
  clusters with only 1 or 2 cells.

- Fix to add_to_seurat so that unsorted chromosomes that are only named
  by their number are read as character and not integers, which
  resulted in the named features using the unsorted index as a name.

- Allow adding a custom column prefix when using add_to_seurat (by
  @matt-sd-watson )

- Add an "assay_name" option to add_to_seurat in case the assay is not
  named "RNA" like expected.

- More helpful message displayed before erroring when trying to run
  add_to_seurat without having run the HMM.

- Add option to filter genes used in subclustering based on zscore.

- Fix i3 settings calculation.

- Change to not output "scaled" CNV proportion values when running the
  HMM in i3 mode as we don't know if a gain/loss is of only 1 copy or
  more.

- Update actually used "levels" in the factor used to store gene
  chromosome information so that contigs entirely filtered out don't
  mess up the color bars on the heatmap.

- Linked write_expr_matrix option in plot_cnv to the methods that plot
  the observations and references so that the split matrix also don't
  get output when set to false.

- Fix for plot_cnv(plot_chr_scale=T) to handle chromosomes with a
  single gene on them, though such chromosomes should be dropped
  entirely.

[InteractiveComplexHeatmap](/packages/InteractiveComplexHeatmap)
-------------------------

                        Changes in version 1.5.1                        

- `do_default_click_action` and `do_default_bruch_action` are reset to
  `FALSE`
  if `subHeatmapOutput()` is not used.

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.32.0                        

NEW FEATURES

- splitAsList() can now perform a "dumb split", that is, when
  no split factor is supplied, 'splitAsList(x)' is equivalent
  to 'unname(splitAsList(x, seq_along(x)))' but is slightly more
  efficient.

SIGNIFICANT USER-VISIBLE CHANGES

- Add ellipsis argument (...) to the gaps() generic function.

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.7.6 (2022-10-25)                  

- Fixed minor bugs and typos

                 Changes in version 1.7.5 (2022-10-05)                  

- Fixed build issues

                 Changes in version 1.7.4 (2022-10-04)                  

VISIBLE USER CHANGES

- Progress bars for long processing functions are now implemented via
the package progressr, added a wrapper function for fast enabling
progress bars, enable_progress_bars()
- Introduced logging for issues in HSC_population_size_estimate() -
signals eventual problems in computing estimates and why

BUG FIXES AND MINOR CHANGES

- Fixed minor bugs and typos

                 Changes in version 1.7.3 (2022-06-17)                  

BUG FIXES AND MINOR CHANGES

- All functions that check for options now have a default value if
option is not set
- CIS_grubbs function is now faster (removed dependency from
psych::describe)

NEW

- New functions CIS_grubbs_overtime() and associated plotting
function
top_cis_overtime_heatmap() to compute CIS_grubbs test over time

                 Changes in version 1.7.2 (2022-05-23)                  

BUG FIXES AND MINOR CHANGES

- Fixed minor issues in import_association_file() - function had
minor
issues when importing *.xlsx files and missing optional columns
threw errors
- Fixed bug in as_sparse_matrix() - function failed when trying to
process an aggregated matrix

NEW

- Added 2 new utility functions export_ISA_settings() and
import_ISA_settings() that allow a faster workflow setup

                 Changes in version 1.7.1 (2022-05-04)                  

BUG FIXES AND MINOR CHANGES

- Fixed minor issue in compute_near_integrations() - function errored
when report_path argument was set to NULL
- Fixed dplyr warning in integration_alluvial_plot() internals
- Fixed issue with report of VISPA2 stats - report failed due to
minor
error in rmd fragment
- Internals of remove_collisions() use again dplyr internally for
joining and grouping operations - needed because of performance
issues with data.table
- fisher_scatterplot() has 2 new arguments that allow the disabling
of
highlighting for some genes even if their p-value is under the
threshold

[iSEE](/packages/iSEE)
----

                       Changes in version 2.9.12                        

- Export constants used in iSEEu.

                       Changes in version 2.9.11                        

- Fix R CMD check warnings about missing documentation.

                       Changes in version 2.9.10                        

- Enable customisation of hovering tooltip in DotPlot panels
including
colData or rowData information.

                        Changes in version 2.9.9                        

- Allow screenshots in vignettes to use full width of pkgdown site.

                        Changes in version 2.9.8                        

- Enable autocompletion for feature names in the heatmap feature
selection modal.

                        Changes in version 2.9.7                        

- Update FontAwesome icon question-circle to circle-question (v6).

                        Changes in version 2.9.6                        

- Bugfix related to https://github.com/rstudio/shiny/issues/3125,
mainly applicable to custom landing pages that render a
DT::datatable() prior to selectInput().

                        Changes in version 2.9.5                        

- Bugfix reverting a change in 2.9.3 breaking re-rendering of
reactivated panels.
- Complete bugfix to prevent unnecessary re-rendering of
ComplexHeatmapPlot panel when dimension of an incoming multiple
selection is dismissed by the options of the child panel.

                        Changes in version 2.9.4                        

- Bugfix setting the active multi-selection info of Table panels to a
fixed message, as the panel is not re-rendered when search boxes are
used.
- Bugfix re-rendering ComplexHeatmapPlot panels when displaying
incoming column selection.

                        Changes in version 2.9.3                        

- Partial bugfix avoiding re-rendering of ComplexHeatmapPlot panel
when an incoming row selection changes if custom rows are in use.
The partial bugfix only applies if the ComplexHeatmapPlot also
disables the restriction on any incoming column selection.

                        Changes in version 2.9.2                        

- Document the existing panel modification modes.

                        Changes in version 2.9.1                        

- Add spinner to ComplexHeatmapPlot.

[iSEEu](/packages/iSEEu)
-----

                        Changes in version 1.9.3                        

- Fix numbering issue affecting the coloring of LogFCLogFCPlot
panels.

                        Changes in version 1.9.2                        

- Depend on the iSEEhex package to initiate the "iSEEverse".

[isomiRs](/packages/isomiRs)
-------

                       Changes in version 1.25.1                        

FIX

- remove DiscriMiner dependency and functions since it is not available
  for 3.16

[kebabs](/packages/kebabs)
------

                       Changes in version 1.31.2                        

- changed dependency to 'Matrix' package (now requires at least version
  1.5-0)

                       Changes in version 1.31.1                        

- minor changes in coercions in accordance with latest version of
  'Matrix' package

                       Changes in version 1.31.0                        

- new branch for Bioconductor 3.16 devel

[KEGGREST](/packages/KEGGREST)
--------

                       Changes in version 1.37.0                        

BUG CORRECTION

- 1.37.1 Fixes new endpoint

- 1.37.2 http to https fixes windows error

[limma](/packages/limma)
-----

                       Changes in version 3.54.0                        

- 
  New function goanaTrend() to estimate a covariate-dependent
  trend in the probability of differential expression for use in
  goana() and kegga() with `trend=TRUE`. The trend estimated by
  goanaTrend() is squeezed slightly towards constancy to provide
  stability when the number of DE genes is small. The amount of
  squeezing decreases with the number of DE genes.
  
  goana() and kegga() now call goanaTrend() when the `trend`
  argument is used. When `plot=TRUE` the plot is now created by
  goanaTrend() instead of by barcodeplot().

- 
  Rename argument `prior.prob` to `null.prob` in goana() and
  kegga().

- 
  Update goana() and kegga() code to catch NA gene IDs or
  covariate values.

- 
  goana() and kegga() no longer give an error when a trend is
  estimated but there are no DE genes.

- 
  kegga() now uses https instead of http links when reading from
  rest.kegg.jp.

- 
  Add new argument `fc` to topTable() and topTableF().  Change
  default for `lfc` to NULL instead of 0.

- 
  lmFit() now checks explicitly for NAs in the design matrix.

- 
  Update references listed on help pages to use DOIs instead of
  URLs.

[MassSpecWavelet](/packages/MassSpecWavelet)
---------------

                 Changes in version 1.63.6 (2022-10-15)                 

- Fix regression in cwt() when scales has length 1.

                 Changes in version 1.63.5 (2022-10-11)                 

- getRidge() supports a scaleToWinSize parameter. This argument
controls how scales get mapped to window sizes. These windows are
used to track the local maxima into ridges. MassSpecWavelet had a
criteria of winsize <- 2*scale+1, while xcms modified it to winsize
<- floor(scale/2). This new argument enables xcms maintainers to
call MassSpecWavelet's getRidge (if they want to) using their
criteria, while it still lets us preserve backwards compatibility in
our results. See ?getRidge for further details.

- The getLocalMaximumCWT() is_amp_thres_relative parameter is now
isAmpThreshRelative, for consistency with other parameter
capitalization in the package. Since it was introduced 10 days ago,
I don't think there will be more than one user using it.

- getLocalMaximumCWT() and peakDetectionCWT have a
exclude0scaleAmpThresh parameter. When computing the relative
amp.Th, if this parameter is set to TRUE, the amp.Th will exclude
the zero-th scale from the max(wCoefs). The zero-th scale
corresponds to the original signal, that may have a much larger
baseline than the wavelet coefficients and can distort the threshold
calculation. The default value is FALSE to preserve backwards
compatibility.

- peakDetectionCWT lets the user pass custom arguments to getRidge().

                 Changes in version 1.63.4 (2022-10-10)                 

- The improvements in localMaxima() and cwt() provide significant
speed-ups to peakDetectionCWT() as well as better scalability.

- A prepareWavelets() function lets the user pre-compute the daughter
wavelets for more efficient cwt() calculations when applied on
multiple spectra. When used transforming 1000 spectra, of 2000
points long each, using 25 different scales, cwt() is twice as fast
as in previous versions. Further improvements to avoid some memory
allocations are still feasible in future versions.

- Through the prepareWavelets() function, we provide the
extendLengthScales argument, that provides the same functionality
than the extendLengthMSW argument in xcms:::MSW.cwt().

- The peakDetectionCWT() function accepts a prepared_wavelets object
in the scales argument for better efficiency.

                 Changes in version 1.63.3 (2022-10-10)                 

- localMaxima() has a more efficient implementation of the algorithm,
now being 10x faster than before, while giving the same results.

- Experimentally, localMaxima() can use a new and different algorithm
for detecting local maxima. See the new "Finding local maxima"
vignette for further details.

                 Changes in version 1.63.2 (2022-09-29)                 

- Let getLocalMaximumCWT() have a relative amp.Th. Related to #4.
- Fix bug on identifyMajorPeaks() where nearbyWinSize was forced
to 150 if nearbyPeak was set to TRUE. Related to #4.
- Added excludeBoundariesSize argument to identifyMajorPeaks().
Before, nearbyWinSize was used for two different but related
criteria: the range for including peaks close to a large peak AND
the range to exclude peaks close to the beginning and end of the
signal. Now, we have two independent arguments for each setting. The
current behaviour does not change, but it is now more flexible.
Related to #4.

                 Changes in version 1.63.1 (2022-07-08)                 

- Drop sav.gol alias from the documentation. Related to #2
- Explain how the scales argument is defined in relation to the
mother
wavelet.
- Add sessionInfo() to the vignette

[MatrixGenerics](/packages/MatrixGenerics)
--------------

                        Changes in version 1.9.1                        

- Fix for functions whose first argument is not x
  (<https://github.com/Bioconductor/MatrixGenerics/issues/28> and
  <https://github.com/Bioconductor/MatrixGenerics/pull/29>).

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.5.9 (2022-10-19)                  

- add ShinyApps to biocViews

                 Changes in version 1.5.8 (2022-09-30)                  

- replace vegan::metaMDS by MASS::isoMDS for NMDS dimension reduction

                 Changes in version 1.5.7 (2022-09-29)                  

- change renderUI functions to update...Input whenever possible

- remove the dual interface for measured/missing values help pages in
  the
  values tab

- improve the formula/expression checks in the DE tab using the
  model.matrix and makeContrasts functions

                 Changes in version 1.5.6 (2022-09-23)                  

- speed up non-server functions

                 Changes in version 1.5.5 (2022-08-23)                  

- add unit tests for imputation of missing values

- bug fix (ncol instead of nrow) for calculation of BPCA imputation

                 Changes in version 1.5.4 (2022-08-18)                  

- improve the imputation of missing values

                 Changes in version 1.5.3 (2022-08-01)                  

- only use a subset of maximum 5000 features in the calculation of
  Hoeffding's
  D values. In case there are less than 10000 features in the
  SummarizedExperiment object, all features of the SummarizedExperiment
  are
  taken

- createBoxplot precalculates the values of the boxplot instead on
  relying
  on geom_boxplot for calculation of the statistics

                 Changes in version 1.5.2 (2022-07-20)                  

- remove the functions biocrates, maxQuant, and spectronaut from
  MatrixQCvis
  and move to the MatrixQCUtils package

                 Changes in version 1.5.1 (2022-07-18)                  

- load package MatrixQCvis in report_qc.Rmd

[matter](/packages/matter)
------

                 Changes in version 1.99.2 (2022-10-31)                 

BUG FIXES

- Fix 'chunkApply()' behavior for when nchunks == 1

- Fix 'atoms' failure when 'source' is a compressed 'drle'

- Properly ignore missing groups in 'rowsweep()' and 'colsweep()'

                 Changes in version 1.99.1 (2022-10-30)                 

BUG FIXES

- Minor bug fixes to 'matter_list' and 'struct()' behavior

- Minor bug fix to endomorphic subsetting of sparse matrices

- Minor bug fix to 'binvec()' behavior at vector endpoints

                 Changes in version 1.99.0 (2022-10-23)                 

NEW FEATURES

- Complete re-implementation of sparse array code in C++

- Complete re-implementation of matter object C++ backends

- New 'sparse_arr' class with 'sparse_mat' and 'sparse_vec'

- New 'matter_arr' class with 'matter_mat' and 'matter_vec'

- New 'asearch()' function for approximate search with interpolation

- New 'chunkApply()' function replacing 'chunk_apply()', etc.

- New 'colsweep()' and 'rowsweep()' functions w/ group parameter

- New 'colscale()' and 'rowscale()' functions w/ group parameter

- New 'drle_fct' class for delta-run-length encoded factors

- New interpolation options for 'sparse_arr' objects

- Simplified deferred arithmetic operations interface

SIGNIFICANT USER-VISIBLE CHANGES

- New simplified constructors for most matter objects

- Updated S4 internals for most matter objects

- Updated S4 internals for sparse array objects

- Deprecated some rarely-used functions and classes

- Deprecated some S4 generic conflicts with core BioC packages

BUG FIXES

- Fix to 'stream_var' and 'stream_sd' behavior when n=0 or n=1

[memes](/packages/memes)
-----

                        Changes in version 1.5.2                        

- added a warning in runFimo in case bfile argument matches an
existing local file.

                        Changes in version 1.4.1                        

- added an informative error message during importMeme when
parse_genomic_coords fails with custom fasta import.

[metabCombiner](/packages/metabCombiner)
-------------

                        Changes in version 1.7.1                        

- Addition of filtered features to metabData objects and filtered
  method

- Bug Fix in metabData():
  + zero values now treated as missing in Q calculations (zero = TRUE)

- Bug Fix in calcScores()/ evaluateParams():
  + rtrange calculations now account for missing values

- Bug Fix in reduceTable()/ reduceTableParam():
  + method argument included with "mzrt" option

- small logic fix in labelRows() / reduceTable()
  + previous version caused duplicate feature matches in rare cases

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                         Changes in version 1.1                         

Changes in 1.1.6

- scoreVariables function to return the names of the score variables
in a Matched object.

Changes in 1.1.5

- Fix issues on BioC build machines.

Changes in 1.1.4

- matchSpectra: support a CompDb with parameter target.
- Add CompAnnotionSource classes to support definition of references
to annotation resources.
- Add CompDbSource class defining a reference to a CompDb database.
- matchSpectra: support for CompDbSource with parameter target.

Changes in 1.1.3

- Extend filterMatches framework (issue #86). ScoreThresholdParam
added to perform filtering the matches based on a threshold for the
scores.

Changes in 1.1.2

- lapply and endoapply methods (issue #84). lapply allows to apply
any
function to each subset of matches for each query element and
returns the corresponding list of results. endoapply is similar but
applies a function returning a Matched and returns a Matched
representing updated matches.

Changes in 1.1.1

- Extend filterMatches framework (issue #81). SelectMatchesParam and
TopRankedMatchesParam added to perform respectively manual filtering
and keeping only the best ranked matches for each query element.

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                         Changes in version 1.5                         

MetaboCoreUtils 1.5.2

- substractElements drops elements with zero counts (issue #57).

MetaboCoreUtils 1.5.1

- Add functions formula2mz, adductFormula and multiplyElements (issue
#55).

[MetCirc](/packages/MetCirc)
-------

                 Changes in version 1.27.2 (2022-10-19)                 

- use Spectra from Spectra for spectral representation instead of
  MSpectra
  from MSnbase

- remove compare_Spectra and normalizeddotproduct functions and use
  compareSpectra from Spectra package

- change vignette from Rnw to Rmd

                 Changes in version 1.27.1 (2022-10-12)                 

- add ShinyApps in biocViews

[MetNet](/packages/MetNet)
------

                 Changes in version 1.15.2 (2022-05-18)                 

- improve documentation for function correlation and statistical and
  include information that also ggm can be used

                 Changes in version 1.15.1 (2022-04-27)                 

- set ci in corr.test to FALSE to speed up calculation of correlation
  values
  (function correlation)

[mia](/packages/mia)
---

                         Changes in version 1.5                         

- Added HintikkaXOData

- Added sample metadata option to getExperimentCrossAssociation

- estimateFaith: add support for multiple rowTrees

- calculateRDA/CCA: added variable argument & replaced altexp argument
  with altExp

- getExpCrossCorr: bugfix; samples should match when correlations
  between features are calculated

- getExpCrossCorr: Kendall's tau is the default method

- mergeSEs: bugfix; links between trees and rows/cols were wrong &
  rowData did not include all info

- calculateDPCoA, calculateUnifrac & merge: add support for multiple
  trees

- altExp parameter to altexp

- agglomerateByRank: make rownames unique by default

- removed calculateDistance and calculateUniFrac alias

[miaViz](/packages/miaViz)
------

                         Changes in version 1.5                         

- plot*Tree & *TreeData: Add support for multiple trees

- plot*Tree layout bugfix

[microbiomeExplorer](/packages/microbiomeExplorer)
------------------

                        Changes in version 1.7.1                        

- removed ZILN method as metagenomeSeq's fitFeatureModel runs into
  issue with changed approach to compare ncol and nrow in limma's lmFit

- cast matrixes into data.frame before calling lmFit to circumvent
  identical check for ncol & nrow in lmFit (if dimensions are named,
  results are equal, but not identical)

[microbiomeMarker](/packages/microbiomeMarker)
----------------

                        Changes in version 1.3.2                        

- fix error on subgroup in lefse, #62, #55

                 Changes in version 1.3.1 (2022-05-26)                  

- Development version on Bioconductor.

                 Changes in version 1.2.1 (2022-05-26)                  

- Confounder analysis.
- Comparison of different methods.

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                        Changes in version 1.9.5                        

- fix the color of mp_plot_diff_boxplot and update mp_plot_abundance.
(2022-10-27, Thu)

                        Changes in version 1.9.4                        

- keep the consistent color between the panels of
mp_plot_diff_boxplot. (2022-10-18, Tue)
- optimizing the import of dtplyr. (2022-10-14, Fri)
- update mp_plot_diff_res to support the custom DAA results.
(2022-09-27, Tue)
- update show method of MPSE to avoid the colname advice.
(2022-09-27,
Tue)
- add mp_dmn, mp_dmngroup, mp_divergence and update
mp_plot_diff_boxplot, mp_plot_diff_cladogram to support the custom
style. (2022-09-25, Sun)
- update mp_cal_dist to support specifying distmethod to a function.
(2022-09-20, Tue)
- update left_join to support joining the dist class. (2022-09-20,
Tue)

                        Changes in version 1.9.3                        

- fix a bug of mp_plot_diff_boxplot when taxatree slot is NULL.
(2022-08-22, Mon)

                        Changes in version 1.9.2                        

- fixed the local vignettes. (2022-07-06, Wed)
- add mp_import_biom to build MPSE class from biom-format file.
(2022-07-13, Wed)
- add 'mp_plot_diff_boxplot' to replace ggdiffbox. (2022-07-29, Fri)

                        Changes in version 1.9.1                        

- fix the color of legend in mp_plot_diff_cladogram. (2022-05-25)
- add mp_plot_diff_cladogram in vignetters. (2022-05-14)
- fixed a bug when the total counts of sample is less than chunks in
mp_plot_rarecurve. (2022-05-14)
- add rm.zero argument in mp_plot_abundance to control whether mask
the zero abundance of species. (2022-05-06)
- add taxa.class argument in mp_diff_analysis to test the specified
taxa level. (2022-04-28)

[midasHLA](/packages/midasHLA)
--------

                        Changes in version 1.5.2                        

- alleles and KIR frequenceis updated.

- failed test fixed after update in tide dependence.

[minfi](/packages/minfi)
-----

                        Changes in version 1.43                         

- v1.43.1 Tightened the recognition code for 27k data to avoid
  conflict with the next release of the Allergy and Asthma array.

[miRspongeR](/packages/miRspongeR)
----------

                     Changes in version 2.1.1-2.1.3                     

- Add new reference and fix bug <2022-09-22, Thur>.

[mirTarRnaSeq](/packages/mirTarRnaSeq)
------------

                 Changes in version 1.5.1 (2022-07-21)                  

- added stickers and README.md file

[mistyR](/packages/mistyR)
------

                         Changes in version 1.5                         

- Added different modeling functions. Might not be completely
backwards compatible!

[Modstrings](/packages/Modstrings)
----------

                 Changes in version 1.13.1 (2022-08-13)                 

- fixed bug in combineIntoModstrings

[monaLisa](/packages/monaLisa)
--------

                        Changes in version 1.3.1                        

- update citation information

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.41.1                        

- Fix the issue of importMatrix for mixed cases of alphabet in MEME
  format.

[msa](/packages/msa)
---

                       Changes in version 1.29.3                        

- fix for possibly malformed inputs: all sequences are forced to
  uppercase characters (previously, ClustalW and ClustalOmega produced
  wrong results when called with lowercase sequences)

                       Changes in version 1.29.2                        

- fix in texshade.sty as suggested on TeXshade homepage at CTAN

                       Changes in version 1.29.1                        

- fix in argtable library (ClustalOmega) to avoid compilation errors on
  newest Mac OS

                       Changes in version 1.29.0                        

- new branch for Bioconductor 3.16 devel

[MSA2dist](/packages/MSA2dist)
--------

                 Changes in version 1.1.7 (2022-10-05)                  

Major changes

- changed licence into GPL-3 to account for KaKs Calculator 2.0
licence

Minor improvements and bug fixes

                 Changes in version 1.1.6 (2022-09-22)                  

Major changes

Minor improvements and bug fixes

- changed test-rcpp_KaKs

                 Changes in version 1.1.5 (2022-09-19)                  

Major changes

Minor improvements and bug fixes

- changed rcpp_KaKs data access from data.frame to list

                 Changes in version 1.1.4 (2022-07-16)                  

Major changes

Minor improvements and bug fixes

- fixed rcpp Warnings

                 Changes in version 1.1.3 (2022-07-11)                  

Major changes

Minor improvements and bug fixes

- fixed rcpp indentation Warnings

                 Changes in version 1.1.2 (2022-07-05)                  

Major changes

- Added KaKs Calculator 2.0 models as rcpp implementation

Minor improvements and bug fixes

- fixed some typos

                 Changes in version 1.1.1 (2022-06-30)                  

Major changes

Minor improvements and bug fixes

- Added cds2codonaln

[MsBackendMgf](/packages/MsBackendMgf)
------------

                         Changes in version 1.5                         

Changes in 1.5.1

- Fix export method to fail if one or more columns contain either S4
classes or list-like structures.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                         Changes in version 1.9                         

MsCoreUtils 1.9.2

- feat: imputation is compatible with HDF5Matrix objects
- feat: normalization is compatible with HDF5Matrix objects
- feat: matrix aggregation is compatible with HDF5Matrix objects
- fix+feat: aggregate_by_matrix now correctly handles missing data
and
implements 'na.rm'
- Fix rla/rowRla man page.

MsCoreUtils 1.9.1

- Random forest imputation (using missForest) is now available
(`method = "RF")

MsCoreUtils 1.9.0

- New Bioc devel version

[msImpute](/packages/msImpute)
--------

                        Changes in version 1.7.1                        

- add width and shift params of the gaussian distribution in the
  down-shift mode of msImpute to provide more control over the
  shape of the distribution when observed peptides have high
  average expression and missing is group-specific.

[mslp](/packages/mslp)
----

                       Changes in version 0.99.6                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

FIXES

- Added unit test for comp_slp and corr_slp.

- Updated vignette.

- Various fixes to reviewers questions.

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.23                         

Changes in 2.23.2

- Fix robust aggregation for MsCoreUtils 1.9.2.

Changes in 2.23.1

- Import spectrapply generic from ProtGenerics.

Changes in 2.23.0

- New release (Bioc 3.16)

[msqrob2](/packages/msqrob2)
-------

                        Changes in version 1.5.3                        

- Added the option to use mixed models without ridge regression
- Added the option to use rowdata variables in the mixed models

                        Changes in version 1.5.1                        

- Fix weighted variance covariance matrix and QR decomposition in the
msqrobLmer function

[MSstatsTMT](/packages/MSstatsTMT)
----------

                 Changes in version 2.4.1 (2022-09-12)                  

- Minor change: fix the bug of data table in performing local
  normalization

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.24.0                        

New features

- replicates provides the actual colnames identified as replicate
observations for a particular biological unit in the sampleMap

Bug fixes and minor improvements

- Added an assay<- replacement method for robustifying
saveHDF5MultiAssayExperiment with plain matrices
- Use BiocBaseUtils::setSlots and avoid warnings of triple colon use.
- Resolve issue when colData has one column when merging two
MultiAssayExperiment objects, i.e., using the c method (@cvanderaa,
#315)
- Increase efficiency in colnames and rownames methods (@cvanderaa,
#314)
- Make 'prefix' inputs consistent in saveHDF5MultiAssayExperiment and
loadHDF5MultiAssayExperiment (@asiyeka, #313)
- Improve performance for replicated method
- Update wideFormat documentation, when replicates present additional
sets of columns will be appended to the produced DataFrame (@DarioS,
#312)

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.5.18                        

Bug fix

- GHA fix.

                       Changes in version 1.5.17                        

New features

- By default ES taken as BETA new parameter added so users can
specify
if this isn't the case (es_is_beta). If set to FALSE, mapping
removed.
- Imputing BETA ordering has been changed so log(OR) will be sued
before calculating from Z, SE.

                       Changes in version 1.5.16                        

New features

- A new method for computing the Z-score of a sumstats (compute_z
input) has been added: BETA/SE. To use it set compute_z = 'BETA' to
continue to use the P-value calculation use compute_z = 'P'. Note
the default is stil compute_z = FALSE.

Bug fix

- Remove erroneous print statement.

                       Changes in version 1.5.15                        

Bug fixes

- Fix NA representation for tabular outputs - By default,
data.table::fread() leaves NAs blank instead of including a literal
NA. That's fine for CSVs and if the output is read in by fread, but
it breaks other tools for TSVs and is hard to read. Updated that and
added a message when the table is switched to uncompressed for
indexing.

                       Changes in version 1.5.14                        

New features

- read_header:
- Can now read entire files by setting n=NULL.
- Improved reading in of VCF files (can read .vcf.bgz now).
- Now exported.
- Added unit tests.
- Remove seqminer from all code (too buggy).
- Automatically remove residual .tsv files after tabix indexing.
- import_sumstats:
- Use @inheritDotParams format_sumstats for better documentation.
- parse_logs: Added new fields.
- format_sumstats: Added time report at the end (minutes taken
total).
Since this is a message, will be included in the logs, and is now
parsed by parse_logs and put into the column "time".

Bug fixes

- index_tabular: Fixed by replacing seqminer with Rsamtools.
- When SNP ID's passed with format 1:123456789, it will now be dealt
with appropriately.
- compute_n can't handle SNP level N values for imputation only
population level. An explanatory error message has now been added.

                       Changes in version 1.5.13                        

Bug fixes

- Special characters causing issues with find empty columns function.
Now fixed.

                       Changes in version 1.5.12                        

Bug fixes

- Mitchondrial (MT) SNPs' chromosome value were being forced to NA by
sort_coords function. This has been fixed.

                       Changes in version 1.5.11                        

Bug fixes

- Had to pass check_dups to other checks so they also wouldn't be
run.
Now independent of non-biallelic check.

                       Changes in version 1.5.10                        

New features

- check_dups parameter added so duplicates won't be removed if
formatting QTL datasets

                        Changes in version 1.5.9                        

Bug fixes

- validate_parameters checks for incorrect version of dbSNP package,
corrected.

                        Changes in version 1.5.6                        

Bug fixes

- MSS can now impute CHR, BP at a SNP level. For cases where CHR
and/or BP are NA but the RS ID is present, these will now be imputed
fromt he reference genome. Note previously, this imputation was done
when the chr and/or bp column was missing.
- Print statement from liftover silenced when no liftover required
- check missing data function will no longer remove cases with NA's
in
SNP_INFO column. The SNP_INFO column is created by MSS for cases
with RS ID and some other information in the same SNP column (like
rs1234:.....). Rather than throw out this info, it is stored in a
new column - SNP_INFO. However, the remove missing data function was
also looking in this column to remove SNPs. This has been corrected.
- find_sumstats():
- Fix N column in metadata.

                        Changes in version 1.5.5                        

New features

- save_format parameter created for format_sumstats. This will
replace
ldsc_format which is now deprecated. Use save_format="LDSC" instead.
Other options for save_format are generic standardised (NULL) and
IEU Open GWAS VCF format ("openGWAS").
- dbSNP version 155 has now been added. Users can now control the
version of dbSNP to be used for imputation (144 or 155). Note that
with the 9x more SNPs in dbSNP 155 vs 144, run times will increase.

Bug fixes

- Change where sex chromosomes were made lower case removed to match
UCSC

                        Changes in version 1.5.4                        

New features

- Further mappings added

Bug fixes

- Duplication of non-bi-allelic and indels fixed
- Correct compute_nsize documentation

                        Changes in version 1.5.1                        

New features

- Export vcf2df.
- Move some post-processing function inside this function (e.g.
drop duplicate cols/rows).
- read_vcf can now be parallised: splits query into chunks, imports
them, and (optionally) converts them to data.table before rbinding
them back into one object.
- Added report of VCF size (variants x samples) before processing
to give user an idea of long it will take to process.
- Added arg mt_thresh to avoid using parallelisation when VCFs are
small, due to the overhead outweighing the benefits in these
cases.
- Added Linux installation instructions for axel downloader.
- Added 2nd tryCatch to downloader with different download.file
parameters that may work better on certain machines.
- Avoid using file.path to specify URL in:
- get_chain_file
- import_sumstats
- Allow download_vcf to pass URLs directly (without downloading the
files) when vcf_download=FALSE.
- download_vcf:
- Make timeout 10min instead of 30min.
- Make axel verbose.
- load_ref_genome_data:
- Give more informative messages that let user know which steps
take a long time.
- Speed up substring preprocessing.
- read_vcf_genome: more robust way to get genome build from VCF.
- read_sumstats: Speed up by using remove_empty_cols(sampled_rows=),
and only run for tabular file (read_vcf already does this
internally).

Bug fixes

- select_vcf_field: Got rid of "REF col doesn't exists" warning by
omitting rowRanges.
- Ensured several unevaluated code chunks in
vignettes/MungeSumstats.Rmd were surrounding by ticks.
- vcf2df: Accounted for scenarios where writeVcf accidentally
converts
geno data into redundant 3D matrices.
- Use data.table::rbindlist(fill=TRUE) to bind chunks back
together.
- Remove unused functions after read_vcf upgrades:
- infer_vcf_sample_ids
- is_vcf_parsed
- check_tab_delimited
- read_vcf_data
- remove_nonstandard_vcf_cols
- Remove redundant dt_to_granges by merging functionality into
to_granges.
- Adjusted liftover to accommodate the slight change.
- Fix is_tabix (I had incorrectly made path all lowercase).
- Let index_vcf recognize all compressed vcf suffixes.
- Add extra error handling when .gz is not actually
bgz-compressed.
- Set BiocParallel registered threads back to 1 after
read_vcf_parallel finishes, to avoid potential conflicts with
downstream steps.

                        Changes in version 1.5.0                        

New features

- Added "query" column to find_sumstats output to keep track of
search
parameters.
- import_sumstats:
- Check if formatted file (save_path) exists before downloading to
save time.
- Pass up force_new in additional to force_new_vcf.
- Updated Description tag in DESCRIPTION file to better reflect the
scope of MungeSumstats.
- Upgraded read_vcf to be more robust.
- Edited Deps/Suggests
- Elevate IRanges to Imports.
- Remove stringr (no longer used)
- Add new internal function is_tabix to check whether a file is
already tabix-indexed.
- read_sumstats:
- now takes samples as an arg.
- Parallises reading VCF using GenomicFiles.
- read_sumstats: now takes samples as an arg.
By default, only uses first sample (if multiple are present in
file).
- Remove INFO_filter= from ALS VCF examples in vignettes (no longer
necessary now that INFO parsing has been corrected).
- download_vcf can now handle situations with vcf_url= is actually a
local file (not remote).

Bug fixes

- AF (allele frequency) was accidentally being assigned as INFO
column
in VCFs where the INFO rows started with "AF". This caused a large
number of SNPs to be incorrectly dropped during the check_info_score
step.
- If INFO score is not available, INFO column is now dropped entirely
(rather than assigning all 1s).
- Adjusted test-vcf_formatting to reflect this. This avoids
ambiguity about whether the INFO score is real or not.
- check_info_score:
- Added extra messages in various conditions where INFO is not
used for filtering, and don't add log_files$info_filter in these
instances.
- Added unit tests.
- check_empty_cols was accidentally dropping more columns than it
should have.
- Fix GHA pkgdown building:
- The newest version of git introduced bugs when building pkgdown
sites from within Docker containers (e.g. via my Linux GHA
workflow). Adjusting GHA to fix this.
- Fix write_sumstats when indexing VCF.
- Ensure read_sumstats can read in any VCF files (local/remote,
indexed/non-indexed).
- Fix test-vcf_formatting.R
- line 51: had wrong AF value in string
- line 109: encountering error? due to duplicate SNPs?
- Fix test-check_impute_se_beta
- lines 51/52: setkey on SNP (now automatically renamed from ID by
read_vcf).
- Fix test-read_sumstats:
- standardising of headers is now handled internally by
read_sumstats.
- Ensure CHR is a character vector when being read in.
- line 44: Ensure extra cols in vcf_ss are dropped.
- parse_logs: Add lines to parsing subfunctions to allow handling of
logs that don't contain certain info (thus avoid warnings when
creating the final data.table).
- 'Avoid the use of 'paste' in condition signals' fixed:
- check_pos_se
- check_signed_col
- Used to rely on gunzip to read bgz files, but apparently this
functionality is no longer supported (possibly due to changes to how
Rsamtools::bgzip does compression in Bioc 3.15. Switched to using
fread + readLines in:
- read_header
- read_sumstats
- read_header: wasn't reading in enough lines to get past the VCF
header. Increase to readLines(n=1000).
- read_vcf: Would sometimes induce duplicate rows. Now only unique
rows are used (after sample and columns filtering).
- Issue with mix of chr:bp:a1:a2 and chr:bp and rs id resolved

[muscat](/packages/muscat)
------

                       Changes in version 1.11.1                        

- bug fix in pbHeatmap(): previously failed for results from 'mmDS()'

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 2.4.0                        

- Fixed plot_region_heatmap() producing the wrong plot when a factor
is used for the chromosome.
- Fixed nanopolish and f5c import positions being off by 1.
- Fixed broken samples() setter for NanoMethResults.
- Added plot_agg_genes() function as a shorthand for
plot_agg_regions(x, exons_to_genes(exons(x))).
- Added the ability to interrupt methy_to_bsseq() calls.
- Added handling for NanoMethResults objects in filter_methy(). If
NanoMethResult is used as input, then NanoMethResult is invisibly
returned as output.
- Added black outlines to exons in annotation to distinguish
contiguous segments for features like tandem repeats.
- Added line_size argument to plot_gene(), plot_region() and
plot_granges() plots for adjusting line size.
- Added subsample argument to heatmap plots, default 50. This reduces
the number of rows shown the plot to the specified amount.
- Added get_exons_mm10(), get_exons_hg19(), and get_exons_hg38() as
replacements for get_exons_mus_musculus() and
get_exons_homo_sapiens().
- Changed heatmaps to no longer plot samples that are absent from
sample annotations.
- Changed heatmap labels to appear on the right rather than on top.
- Changed heatmap alpha from 0.33 to 0.5.
- Changed arrows in exon connectors to appear in the middle as open
arrow instead of at the end as closed arrow.
- Changed default X axis labels to be rescaled to appropriate
SI-style. e.g. Kb, Mb, Gb.

[NanoTube](/packages/NanoTube)
--------

                        Changes in version 1.3.7                        

- In addition to ruv::RUVIII, the RUVSeq::RUVg method can now be used
for data normalization. Options have also been added to allow tuning
of the parameters for these methods. More details are provided in
the vignette.
- A csv or txt file containing a design matrix can be input as a
'sampleTab' in processNanostringData(). This facilitates easier
differential expression analysis with more complex models, and an
example has been added to the vignette.
- Various other improvements to vignette.

                        Changes in version 1.3.6                        

- NanoTube can now process zipped and tarred (.zip or .tar)
directories, as well as gzipped (.gz) RCC files, such as those
downloaded from GEO in many cases.

                        Changes in version 1.3.5                        

- Corrected a bug that caused NanoTube not to recognize reporters
labeled as "Endogenous1", "Endogenous2", etc. as Endogenous.

[ndexr](/packages/ndexr)
-----

                       Changes in version 1.19.1                        

- **UPDATE: Using the RCX package for working with networks.**
  **Defunct Functions:**

- *rcx_fromJSON:* `RCX::readJSON()`

- *rcx_toJSON:* `RCX::toCX()`

- *rcx_aspect_toJSON:* `rcx_aspect_toJSON`

- *rcx_new:* `RCX::createRCX()`

- *rcx_asNewNetwork:* `RCX::createRCX()`

- *rcx_updateMetaData:* `RCX::updateMetaData()`

- *print.RCX:* `RCX::print.RCX()`

- *rcx_toRCXgraph:* `RCX::toIgraph()`

- *rcxgraph_toRCX* `RCX::fromIgraph()`

[netZooR](/packages/netZooR)
-------

                       Changes in version 1.1.12                        

- Reactivated unit tests for Ubuntu GitHub actions.
- LIONESS can now build single-sample coexpression networks using
@kshutta's implementation
- Fix for ALPACA singleton community case (detected by @talkhanz)
- Fix for CRANE significance test on constant modularity scores
(detected by @talkhanz)
- Improved method description by @kshutta
- Fix for PANDA edge case when only expression is provided

[ngsReports](/packages/ngsReports)
----------

                        Changes in version 2.0.0                        

- Added status bar to all plots from FastQC reports using pachwork

                       Changes in version 1.13.3                        

- Bug fix for importing macs2 logs

                       Changes in version 1.13.2                        

- Bug fix for importing DuplicationMetrics

                       Changes in version 1.13.1                        

- Bug fixes for importing macs2 logs

- Bug fixes for importing bowtie2 logs

[nnSVG](/packages/nnSVG)
-----

                 Changes in version 1.1.2 (2022-05-17)                  

- enable non-SpatialExperiment inputs

[NxtIRFcore](/packages/NxtIRFcore)
----------

                 Changes in version 1.3.2 (2022-10-25)                  

- Announcing that NxtIRFcore will be superceded by SpliceWiz from Bioc
  3.16 onwards.

                 Changes in version 1.2.1 (2022-06-22)                  

- Bugfix for BuildReference: ignores "N" bases when translating codons
  (instead of returning an error)

- Bugfix for BuildReference: more versatile handling of cache locations
  for BiocFileCache

[oncoscanR](/packages/oncoscanR)
---------

                        Changes in version 1.1.0                        

- To simplify the workflow, the gender of the patient is not taken
into account anymore. That implies that in a male sample, a gain
of 3 extra copies on the X or Y chromosome is considered as a gain
and not an amplification anymore. For female samples, nothing
changes.
- The oncoscan coverage has been corrected to reflect only areas
where
there are groups of probes. Isolated probes where causing issues to
identify arm-level alterations as ChAS segments where never extended
to these probes and the 90% threshold could never be met
(particularly on chromosomal arms 9p and Yq).
- Minor corrections in vignette

                        Changes in version 1.0.0                        

- The nLST test has been clinically validated on 384 patients from
the
PAOLA-1 trial and the recommended threshold is now >=15.
- The default value for arm-level alterations has been set to 90% as
mentioned in the publication [Christinat et al., J Mol Diagn 2021;
PMID: 34454110].
- The genomic LOH score (percent of LOH bases) has been added;
score_gloh.
- Adds a flag "no tumor?" if the percentage of altered bases is less
than 1%.
- Package to be released on Bioconductor

                        Changes in version 0.2.0                        

- Novel HRD score (nLST: number of LSTs, normalized by ploidy):
score_nlst
- Change in Oncoscan workflow to use the nLST score and thresholds.
- New function to compute the number of Mb altered (with or without
LOH): score_mbalt
- New function to compute the average copy number: score_avgcn
- New function to estimate the number of whole-genome duplication
events (based on the average copy number and the thresholds defined
by Carter et al.): score_estwgd

[OncoSimulR](/packages/OncoSimulR)
----------

                         Changes in version 4.0                         

- MAJOR change: users can simulate interventions and adaptive therapy,
  as well as define, track, and use in fitness-dependent specifications
  user-defined variables.

- MAJOR change: users can specify birth and death rates as they want,
  including making both or just one frequency-dependent.

- onlyCancer = FALSE by default in all calls to oncoSimul*.
  This is a possible BRAEKING CHANGE.

- Plot fails with a meaningful error if simulations
  had unrecoverable exception or hit max wall time
  or max num tries.

- Improvements in vignette speed.

- Improvements in compilation speed.

                Changes in version 3.99.12 (2022-10-19)                 

- Try to allow vignette to build successfully in the ARM64 Mac build:
  use a try(plot) in chunk fdfmutex2.

- Use apa style citation and a few fixes to bib.

- Minor additions to vignette.

                Changes in version 3.99.11 (2022-10-13)                 

- Increase max.wall.time of two tests, as in very slow Windoze machines
  they can hit max wall time.

                Changes in version 3.99.10 (2022-10-13)                 

- onlyCancer = FALSE by default in all calls to oncoSimul*.
  This is a possible BRAEKING CHANGE.

- Plot fails with a meaningful error if simulations
  had unrecoverable exception or hit max wall time
  or max num tries.

- Faster vignette.

                 Changes in version 3.99.9 (2022-09-30)                 

- Minor vignette split of code, to try to detect
  location of problem with ARM64 Mac build.

                 Changes in version 3.99.8 (2022-09-15)                 

- Using "Unity builds" to decrease build (compilation) time;
  see file inst/miscell/README_Unity_compilation.

                 Changes in version 3.99.7 (2022-09-13)                 

- Use c++14 and, in windows, do not use -O3: try to minimize
  build time.

                 Changes in version 3.99.6 (2022-09-12)                 

- Udpated exprtk.

                 Changes in version 3.99.5 (2022-07-19)                 

- Funding logo added to README.md and vignette.

                 Changes in version 3.99.4 (2022-07-04)                 

- Failing a test in test.Z-intervention.R in Mac OS; did not fix
  correctly previous time.

                 Changes in version 3.99.3 (2022-06-30)                 

- Failing a test in test.Z-intervention.R in Mac OS

                 Changes in version 3.99.2 (2022-06-29)                 

- Vignette: author list and references.

                 Changes in version 3.99.1 (2022-06-25)                 

- Users can specify birth and death rates as they want, including
  making both or just one frequency-dependent (thanks to Alberto
  González Klein).

- Interventions (thanks to Javier Muñoz Haro).

- User-defined variables and adaptive therapy (thanks to Javier López
  Cano).

[orthogene](/packages/orthogene)
---------

                        Changes in version 1.3.4                        

Bug fixes

- Remove Matrix.utils since it's now deprecated.
- Reimplement the dMcast function as a new internal function
within orthogene, since that's the only function I use from
Matrix.utils.
- Fix GHA workflow now that r-lib/actions@master has been removed.

                        Changes in version 1.3.3                        

Bug fixes

- Make test-map_orthologs_babelgene less stringent with the number of
genes expected.

                        Changes in version 1.3.2                        

New features

- Add inst/grofiler_namespace.csv.gz for documentation purposes.

Bug fixes

- create_background:
- Ensure user-supplied bg gets used:
https://github.com/neurogenomics/orthogene/issues/22
- Properly document internal data so that devtools::document doesn't
expect them to be exported objects.

                        Changes in version 1.3.1                        

New features

- plot_orthotree: Pass up tree_source arg.

                        Changes in version 1.3.0                        

New features

- aggregate_mapped_genes:
- Pass up additional args from map_genes.
- Add map_orthologs as a way to create gene_map automatically,
when gene_map=NULL and input_species!=output_species.
- Split species into input_species and output_species args.
- Change method --> agg_method, and use method to pass to
map_orthologs instead.
- Pass up additional args from map_orthologs.
- Add link to detailed explanation of matrix aggregation/expansion
in many2many_rows docs.
- Automatically pick best method for many:1 or many:many mapping.
- as_integers: new arg that uses floor.
- Rename FUN to agg_fun.
- Add new data gprofiler_namespace. Used to validate target= arg in
gprofiler2 functions.
- Upgraded aggregate_mapped_genes:
- Can now used gene_map made by map_orthologs or map_genes.
- Can now handle many:many relationships.
- Will automatically pick the best method to perform aggregation
and/or expansion.
- Removed internal function aggregate_mapped_genes_twice
- Extracted aggregation args from non121_strategy and placed them in
their own new own (agg_fun) since these options are no longer
mutually exclusive due to many:many expansion/aggregation.
- Pass up as_DelayedArray
- Bump to v1.3.0 and R >=4.2 now that we're developing on Bioc 3.16.
- Add ISSUE_TEMPLATE.
- prepare_tree:
- Add tree_source options: path / URL / OmaDB / UCSC / timetree

Bug fixes

- map_genes: Fix report at completion.
- Add safeguards against using aggregation when gene_df isn't a
matrix.
- Removed DelayedMatrixStats Import (no longer needed).
- Fix all unit tests and examples after making all updates.
- Recognize sparse/dense matrix or delayedarray in check_agg_args.

[pathview](/packages/pathview)
--------

                       Changes in version 1.37.1                        

- fix bug in download.kegg caused by recent change in kegg rest API
  urls (http to https).

- korg now include 8282 KEGG species or 1449 new species beyond 2020.

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.11.1 (2022-08-05)                 

- Bugfix: manage mxML files without timestamp tag

[pengls](/packages/pengls)
------

                        Changes in version 1.3.2                        

- With predict function for cv.pengls objects

                        Changes in version 1.3.1                        

- Extend to MSE loss function
- Random fold splitting results in more even folds

[phantasus](/packages/phantasus)
---------

                       Changes in version 1.17.5                        

- All GEO datasets with <100K rows can be loaded

                       Changes in version 1.17.4                        

- Updated DESeq2 usage

[PharmacoGx](/packages/PharmacoGx)
----------

                        Changes in version 3.1.4                        

- Modified downloadPSet function to automatically update the
PharmacoSet class structure and resave the updated object after
download
- This work around is necessary until we can rerun our data
engineering pipelines to regenerate all of our PharmacoSet using
the 3.1.0 package updates
- Added a number of additional methods for computing drug synergy
metrics

                        Changes in version 3.1.1                        

3.1.0

- Update to slot names "cell" -> "sample" and "drug" -> "treatment"
- Update standardized identifier column names to match the above slot
nomenclature: "cellid" -> "sampleid", "drugid" -> "treatmentid"

[PhyloProfile](/packages/PhyloProfile)
------------

                       Changes in version 1.10.3                        

- fixed bug parsing taxon names for input with more than 9999
  taxa

                       Changes in version 1.10.1                        

- fixed bug button for tree uploading disappear #119

- improved parsing seq ID to create links for external DB #117

- fixed bug grepping domains for group comparison fn

- fixed bug parsing data for a list of input gene IDs #122

[plotgardener](/packages/plotgardener)
------------

                        Changes in version 1.3.9                        

BUG FIXES

- Fixed vignette links in "Introduction to plotgardener" vignette.

                        Changes in version 1.3.7                        

BUG FIXES

- plotGenes and related functions will appropriately check for and
handle custom OrgDbs.
- getExons will double-check for appropriate chromosome data to avoid
incorrect plotting based on related chromosome contigs.

NEW FEATURES

- plotManhattan p-value data can be scaled according to a custom
transformation, rather than being limited to -log10.
- plotRanges elements can be ordered randomly or by decreasing width
before plotted row assignment.

                        Changes in version 1.3.6                        

BUG FIXES

- mapColors can appropriately map colors to a numeric vector with the
same values, so long as a range is provided.

                        Changes in version 1.3.5                        

NEW FEATURES

- plotMultiSignal function can plot multiple signal track data sets
in
line with each other.
- calcSignalRange helper function will calculate an appropriate range
for multiple signal data sets.
- pageLayoutRow and pageLayoutCol functions for generating row and
column positions for a number of plot elements.
- Gene transcripts can be highlighted by gene name or transcript name
with the parameter transcriptHighlights in plotTranscripts.

                        Changes in version 1.3.4                        

BUG FIXES

- plotPairsArches Bezier curve height calculations were fixed for
pairs with different sized anchors.

                        Changes in version 1.3.3                        

NEW FEATURES

- plotSignal can now plot negative signal data alone or listed as a
second file.
- A label parameter has been added for plotSignal for convenient
labeling.

                        Changes in version 1.3.2                        

BUG FIXES

- plotSignal range parsing bug fixes were resolved.

- Note about double page rendering has been added to pageGuideHide()
documentation.

                        Changes in version 1.3.1                        

- plotSignal bug fixes related to function not finding posSignal2 and
negSignal2 variables with insufficient data.

- Documentation to introduction vignette has been added to explain
double page rendering when using any removal function, particularly
pageGuideHide().

                        Changes in version 1.3.0                        

Version bump for Bioconductor 3.15 release.

[polyester](/packages/polyester)
---------

                       Changes in version 1.99.3                        

- NB function now exported

- note that version 1.99.3 on GitHub was version 1.1.0 on Bioconductor.

                       Changes in version 1.99.2                        

- bug fix in fragment generation (last 2 bases of transcript were never
  sequenced)

[POMA](/packages/POMA)
----

                       Changes in version 1.7.19                        

- New biocViews and Description
- Call external packages within each POMA function for consistency
- New methods: UMAP, PCR (Principal Components Regression), RNA-seq
analysis (via DESeq2 package)
- Add different statistical methods in PomaVolcano()
- Estimate relative quality weights in PomaLimma()
- Update tests
- Update vignettes
- Several bug fixes

[proDA](/packages/proDA)
-----

                        Changes in version 1.11                         

- Fix copy and paste bug that meant that the wrong value was used in
  the unregularized
  feature estimation procedure (see #18. Thanks to @amackey for
  spotting the bug).

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.37                         

Changes in version 1.37.1

- Fix bug in PerTubro classifiction function (see #146 and #147),
contributed by mgerault.

[pRolocGUI](/packages/pRolocGUI)
---------

                         Changes in version 2.7                         

CHANGES IN VERSION 2.7.2

- Explicitly add remap function for aggregation app

CHANGES IN VERSION 2.7.1

- Move LMB as maintainer
- Update URLs for Github repo
- Removed deprecated cogs icon to gears
- Fix bug in shinyhelper

CHANGES IN VERSION 2.7.0

- New version for Bioc 3.15

[psichomics](/packages/psichomics)
----------

                       Changes in version 1.24.0                        

- Bug fix: psichomics app now opens as expected instead of crashing

                       Changes in version 1.22.1                        

- Bug fix: GTEx data v8 can now be downloaded using newer versions of
R
- Bug fix: UniProt REST API call now uses the updated API
- Bug fix: Updated deprecated setting to remove ggplot2 guide
- Documentation: Clarify how to get survival plots in CLI tutorial

[QDNAseq](/packages/QDNAseq)
-------

                 Changes in version 1.33.1 (2022-04-27)                 

MISCELLANEOUS

- Now the package gives an informative error message when an outdated
  version of the 'future' package is used.  It requires future (>=
  1.22.1).

BUG FIXES

- A few functions used class(x) == "data.frame" rather than
  inherits(x, "data.frame").

                 Changes in version 1.33.0 (2022-04-26)                 

- The version number was bumped for the Bioconductor devel version,
  which is now BioC 3.16 for R-devel.

[QFeatures](/packages/QFeatures)
---------

                         Changes in version 1.7                         

QFeatures 1.7.3

- fix: fixed filterFeatures when selection contains environment
variables

QFeatures 1.7.2

- feat: added 'c' methods to combine QFeatures objects.
- feat: added nrows and ncols methods. Also added use.names argument
(cf ?BiocGenerics::dims)
- docs: improved docs for filterFeatures()
- tests: improved unit tests for filterFeatures()
- feat: added a keep argument in filterFeatures() to control whether
to keep or remove features for assays that do not contain the filter
variable. Also added message printing for a better overview of which
variable were found.
- fix: fixed addAssay() to solve issue #104.
- refactor: refactored addAssay() and dramatically improved the usage
of computational resources.
- feat: colData is automatically transferred from the assay to the
QFeatures object.
- feat: implemented removeAssay() and replaceAssay(). Together with
addAssay(), these functions are used to implement the replacement
method [[<- required to solve issue #57.
- Add CC-BY-SA license for vignettes.

QFeatures 1.7.1

- refactor: imputation now adds a new assay instead of replacing
values.

QFeatures 1.7.0

- New Bioc devel version.

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

                       Changes in version 1.15.3                        

- Fix for error in devtools::check tests caused by differing locales

[quantiseqr](/packages/quantiseqr)
----------

                        Changes in version 1.6.0                        

Other notes

- Fixed the behavior of the mapGenes function when no duplicate names
are found, avoiding unnecessary remapping and having a consistent
behavior with the original implementation

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.22.0                        

- Add as.list and as.data.frame methods to RaggedExperiment to
facilitate extraction of mcols and conversion to table format,
respectively.

[rawrr](/packages/rawrr)
-----

                  Changes in version 1.5 (2022-08-22)                   

- Add centroid.PreferredNoises.

[Rcpi](/packages/Rcpi)
----

                 Changes in version 1.33.2 (2022-07-17)                 

Improvements

- Updated the endpoint URL of UniProt API to fix access issues (#14).

                 Changes in version 1.33.1 (2022-05-07)                 

Improvements

- Remove the Enhances field in DESCRIPTION to improve clarity.

                 Changes in version 1.33.0 (2022-04-25)                 

Improvements

- Fixed a build error on macOS in the devel branch due to
dependencies
not available.

[RCy3](/packages/RCy3)
----

                       Changes in version 2.18.0                        

- New functions:
  - exportPNG
  - exportJPG
  - exportPDF
  - exportPS
  - exportSVG
  - importFileFromUrl
  - selectEdgesAdjacentToNodes

[ReactomePA](/packages/ReactomePA)
----------

                       Changes in version 1.41.1                        

- add function gson_Reactome (2022-7-13, Wed)

[reconsi](/packages/reconsi)
-------

                        Changes in version 1.9.1                        

- Explicitly load Vandeputte data

[rGREAT](/packages/rGREAT)
------

                       Changes in version 1.99.7                        

- Calculation has been speeded up.

                       Changes in version 1.99.6                        

- set `getGeneSetsFromBioMart()` public

                       Changes in version 1.99.2                        

- by default exclude gap regions

- chr prefix in `gr` is removed when `biomart_dataset` is set.

                       Changes in version 1.99.0                        

- add `great()` and related functions to support local GREAT.

- add `shinyReport()` to export results into a shiny application.

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.42.0                        

CHANGES

- Function H5Ocopy() has been included.

- UTF-8 encoded character datsets will be marked as having the
  same encoding when read into an R session.

BUG FIXES

- h5write() no longer truncates multibyte UTF-8 strings (Thanks
  to Aaron Lun @LTLA for reporting this and providing a fix,
  https://github.com/grimbough/rhdf5/issues/111).

[Rhdf5lib](/packages/Rhdf5lib)
--------

                        Changes in version 1.20                         

Bug fixes

- Installation on a system with an existing SZIP library would lead to
  an
  empty linking location when linking against the package. This could
  cause
  downstream compilation problems on some systems.
  (thanks to @brgew for reporting and Philippe Bordron @bordron for
  helpful
  diagnostics https://github.com/grimbough/rhdf5/issues/109)

- Fixed problem where, on some systems, the AEC library would be
  installed
  to a directory named 'lib64' rather than 'lib' as expected by the
  package
  Makevars
  (thanks to Robby Engelmann @robby81 for reporting this,
  https://github.com/grimbough/Rhdf5lib/issues/44)

[Rhisat2](/packages/Rhisat2)
-------

                       Changes in version 1.13.1                        

- Update hisat2 to v2.2.1

- Include SIMD Everywhere library
  (https://github.com/simd-everywhere/simde) to support arm64

[Rhtslib](/packages/Rhtslib)
-------

                        Changes in version 2.0.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Update htslib to 1.15.1 (was 1.7)

[RLSeq](/packages/RLSeq)
-----

                        Changes in version 1.3.0                        

Bioc 3.16 release with the addition of a noise analysis function.

[RMassBank](/packages/RMassBank)
---------

                        Changes in version 3.7.1                        

- Remove getCompTox() because EPA AcTORWS web services have been
  retired

- Fix issue #309 when a FileList had no "mode" column

[rmelting](/packages/rmelting)
--------

                       Changes in version 1.13.3                        

- Revert the Java (>= 11) as SystemRequirements. Removed access to
private KeySet in the code.

                       Changes in version 1.13.2                        

- Added Java (>= 11) as SystemRequirements.

                       Changes in version 1.13.1                        

- Update from MELTING 5.1.1 to MELTING 5.2.0
- Added options for consecutive locked nucleic acids
(method.consecutive.locked and method.consecutive.locked.singleMM).

[RnBeads](/packages/RnBeads)
-------

                       Changes in version 2.15.1                        

- Added custom reimplementations of dye-bias and background correction
  methods from SeSaMe

[rols](/packages/rols)
----

                        Changes in version 2.25                         

CHANGES IN VERSION 2.25.5

- Fail with useful message when ols query can't be converted to
terms.

CHANGES IN VERSION 2.25.4

- Fix error due to failing term().

CHANGES IN VERSION 2.25.3

- Fix failing unit test since change in GO.

CHANGES IN VERSION 2.25.2

- Fix new error when converting an Ontology to data.frame.

CHANGES IN VERSION 2.25.1

- Fix error when converting term with missing slots to data.frame.

CHANGES IN VERSION 2.25.0

- New devel version (Bioc 3.16)

[ropls](/packages/ropls)
-----

                       Changes in version 1.29.32                       

- view: minor modification

                       Changes in version 1.29.30                       

- view: minor modification

                       Changes in version 1.29.28                       

- view: minor modification

                       Changes in version 1.29.26                       

- gg_scoreplot: minor modification

                       Changes in version 1.29.24                       

- gg_scoreplot: minor modification

                       Changes in version 1.29.22                       

- gg_scoreplot: generic method to be applied either to a
  SummarizedExperiment or an opls object

                       Changes in version 1.29.20                       

- gg_scoreplot: modification of the info to display with plotly

                       Changes in version 1.29.18                       

MINOR MODIFICATION

- minor documentation update

                       Changes in version 1.29.16                       

MINOR MODIFICATION

- minor documentation update

                       Changes in version 1.29.14                       

NEW FEATURE

- NCI60 dataset from omicade4 available in MultiAssayExperiment and
  MultiDataSet

                       Changes in version 1.29.12                       

MINOR MODIFICATION

- minor documentation update

                       Changes in version 1.29.10                       

MINOR MODIFICATION

- minor documentation update

                       Changes in version 1.29.8                        

NEW FEATURE

- 'opls' method now handles MultiAssayExperiment

                       Changes in version 1.29.6                        

NEW FEATURE

- 'opls' method now handles SummarizedExperiment

                       Changes in version 1.29.4                        

NEW FEATURE

- 'view' method now handles SummarizedExperiment

- 'sacurine' dataset includes the ExpressionSet and
  SummarizedExperiment formats

                       Changes in version 1.29.2                        

MINOR MODIFICATION

- new bioc release

[rpx](/packages/rpx)
---

                         Changes in version 2.5                         

rpx 2.5.1

- Address failing unit test since PRIDE URL has changed back.

rpx 2.5.0

- New bioconductor devel version

[rrvgo](/packages/rrvgo)
-----

                        Changes in version 1.8.1                        

- New features in scatterPlot():
  - dimensionality reduction with UMAP
  - plot only parents

- Fix picking the right point size in scatterPlot()

- Improve vignette

[Rsubread](/packages/Rsubread)
--------

                       Changes in version 2.12.0                        

- 
  Improve the data structure used by cellCounts to further
  improve its speed.

- 
  More checks for input parameters to prevent cellCounts from
  crashing.

- 
  Improve the screen output of cellCounts.

[RTCGA](/packages/RTCGA)
-----

                       Changes in version 1.27.2                        

- moved createTCGA() out of the package

                       Changes in version 1.27.1                        

- Fixed http and https calls to https://gdac.broadinstitute.org/

[RTCGAToolbox](/packages/RTCGAToolbox)
------------

                       Changes in version 2.28.0                        

New features

- Resolve disparate columns in mutation files from GBMLGG
(@andreagrioni, #45)
- Update codebase to download https resources from newer layout in
gdac.broadinstitute.org (@biostars-nyc, #44).
- Update makeSummarizedExperimentFromGISTIC interface with rownameCol
input, type checking, and set row names when the are not duplicated.

Bug fixes and minor improvements

- Use cache directory obtained from tools::R_user_dir
- Clean up class membership checks, e.g., with is(x, "classname")
- Set missing rows to "" for downstream compatibility with
SummarizedExperiment.

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.36.0                        

NEW FEATURES

- Implement DataFrameFactor objects. See '?DataFrameFactor'.

- Add droplevels() method for Factor objects.

- Factor objects now use a raw vector instead of an integer vector to
  store
  their internal index when they have 255 levels or less. This greatly
  reduces their memory footprint.

SIGNIFICANT USER-VISIBLE CHANGES

- Subsetting a DataFrame object no longer tries to keep its colnames
  unique.

[sangeranalyseR](/packages/sangeranalyseR)
--------------

                        Changes in version 99.1                         

- Base class: SangerReads is designed to store each forward/revers
  reads.

                        Changes in version 1.6.1                        

- Fix chromatogram color issue.

[satuRn](/packages/satuRn)
------

                        Changes in version 1.4.1                        

We report a bug in satuRn 1.4.0. (Bioconductor release 3.15). The bug
was inadvertently introduced in satuRn 1.3.1 (from the former
Bioconductor devel). Note that the bug was not thus present in any of
the older Bioconductor releases 3.13 and 3.14 (satuRn 1.0.x, 1.1.x
and
1.2.x).

Bug details:

Imagine a gene with three isoforms and two cell types. The goal is to
assess DTU between cell types. All isoforms are expressed in all
cells
of cell type 1. However, none of the isoforms are expressed in any of
the cells in cell type 2 (i.e., the gene is not expressed in cell
type 2).

satuRn computes the log-odds of picking a certain isoform from the
pool
of isoforms in each cell type, and then compares these log-odds
estimates between the cell types. However, in this example, the
log-odds
of picking a certain isoform from the pool of isoforms in cell type 2
cannot be computed, as there is no data. Hence, the DTU test
statistic
should be NA. However, due to erroneous handling of NA estimates,
which
was inadvertently introduced in satuRn 1.3.1. while aiming to resolve
github issue 16, the log-odds in cell type 1 will be compared to
zero.
Hence, (erroneous) results can be obtained for this contrast, even
when
there are no data in cell type 2.

Note that in many cases such isoforms may not pass filtering and
would
not get evaluated altogether. However, when analyzing sprase
scRNA-Seq
datasets with a lenient filtering criterium, this problem will apply,
and will result in mistakes in the inference.

[scater](/packages/scater)
------

                       Changes in version 1.26.0                        

- Add projectReducedDim function to project points into an
  existing reduced dimensionality embedding.

- Support "color" and "colour" spellings in all plotting
  functions.

- Add order_by argument to cellwise plot functions.

- Add rasterise argument to plotReducedDim using rasterise.

[SCFA](/packages/SCFA)
----

                 Changes in version 1.7.1 (2022-07-29)                  

- Changed neural network computation from Tensorflow to Torch.

[scp](/packages/scp)
---

                         Changes in version 1.7                         

scp 1.7.4

- Updated CITATION

scp 1.7.3

- refactor: package complies with BiocCheck
- docs: fixed bug in vignette

scp 1.7.2

- Add CC-BY-SA license for vignettes.

scp 1.7.1

- refactor: removed deprecated function rowDataToDF()
- tests: fixed some tests failing because of SCE version differences.
- feat: users can now specify sep when sample names are automatically
generated.

scp 1.7.0

- New devel (Bioc 3.16)

[scRepertoire](/packages/scRepertoire)
------------

                        Changes in version 1.7.2                        

- Rebumping the version change with new release of Bioconductor

- Added mean call to the heatmap of vizGenes()

- To combineTCR, filteringMulti now checks to remove list elements with
  0 cells.

- Removed top_n() call as it is now deprecated, using slice_max()
  without ties.

- Add arrange() call during parseTCR() to organize the chains

- Correct the gd flip in the combineContig and subsequent functions

- Removed viridis call in the clonalNetwork() function that was leading
  to errors

- Matched syntax for strict clonotype in combineBCR()

- Added group.by variable to all applicable visualizations

- Added return.boots to clonalDiversity(), allow for export of all
  bootstrapped values

[scTGIF](/packages/scTGIF)
------

                       Changes in version 1.11.1                        

- Some bugs are fixed in reportTGIF()

[scuttle](/packages/scuttle)
-------

                        Changes in version 1.8.0                        

- Removed support for use.altexps= in aggregateAcrossCells() and
  logNormCounts().

[seqArchR](/packages/seqArchR)
--------

                        Changes in version 1.1.1                        

- (User-facing) Fixed: Too many warnings from Python's scikit-learn
were being printed to console (warning about an expected change in
future scikit-learn version)

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.38.0                        

UTILITIES

- new option 'ext_nbyte' in seqGet2bGeno()

- `seqAlleleCount()` and `seqGetAF_AC_Missing()` return NA instead of
  zero
  when all genotypes are missing at a site

- `seqGDS2VCF()` does not output the FORMAT column if there is no
  selected
  sample (e.g., site-only VCF files)

- `seqGetData(, "$chrom_pos2")` is similar to `seqGetData(,
  "$chrom_pos")`
  except the duplicates with the suffix ("_1", "_2" or >2)

NEW FEATURES

- `seqGDS2BED()` can convert to PLINK BED files with the best-guess
  genotypes when there are only numeric dosages in the GDS file

- `seqEmptyFile()` outputs an empty GDS file

                       Changes in version 1.36.2                        

BUG FIXES

- fix the bug at multi-allelic sites with more than 15 different
  alleles,
  see https://github.com/zhengxwen/SeqArray/issues/78

                       Changes in version 1.36.1                        

BUG FIXES

- `seqExport()` failed when there is no variant

- `seqSetFilter(, ret.idx=TRUE)`, see
  https://github.com/zhengxwen/SeqArray/issues/80

[shinyepico](/packages/shinyepico)
----------

                 Changes in version 1.4.2 (2022-05-18)                  

- Bug fixed: dplyr problem with bisulfite plot and contrasts table

                 Changes in version 1.4.1 (2022-05-10)                  

- Bug fixed: dplyr problem with bisulfite plot and contrasts table

[signeR](/packages/signeR)
------

                        Changes in version 2.0.0                        

- interactive shiny interface: signeRFlow

- Option to work with previously defined signatures, just estimating
  exposures.

- Parallel processing,

- New method to find the number of signatures (or groups) present in
  data, based on statistical testes of BIC values,

- Tests of signature association with survival data or continuous
  variables,

- Exposure-based clustering methods

[simplifyEnrichment](/packages/simplifyEnrichment)
------------------

                        Changes in version 1.7.2                        

- `anno_word_cloud()`: fixed a bug where `exclude_words` is not
  properly
  passed to internal functions.

- `simplifyGOFromMultipleLists()`: fixed a bug where `control` is not
  properly passed to internal functions.

- `simplifyGOFromMultipleLists()`: draw barplots which show numbers of
  signifiant
  GO terms in clusters.

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.7.3 (2022-10-25)                  

- Fixed bugs related to dependency updates

                 Changes in version 2.7.2 (2022-10-19)                  

- Deprecated findMarkerDiffExp(), findMarkerTopTable() and
plotMarkerDiffExp(), which are replaced by runFindMarker(),
getFindMarkerTopTable() and plotFindMarkerHeatmap(), respectively
- Added useReducedDim, detectThresh arguments for find marker
functions
- Deprecated getUMAP() and getTSNE(), which are replaced by runUMAP()
and runTSNE(), respectively
- Added runQuickUMAP() and runQuickTSNE() functions which directly
compute the proper embedding from raw counts matrices with a
simplified argument set
- Added arguments aggregateRow and aggregateCol to plotSCEHeatmap()
- Updated output metadata structure of QC functions, as well as
combineSCE() which merges the new structure properly
- Refined batch correction function set
- Fixed bugs related to UI and console functions

                 Changes in version 2.7.1 (2022-06-29)                  

- Refactored scaling related parts of the workflow, including
variable
feature detection, dimension reduction and 2D embedding
- Redesigned UI landing page and UI running prompt
- Added marker table module across UI
- Added more unit tests
- Fixed bugs in TSCAN UI
- Other minor bug fixes

[SingleR](/packages/SingleR)
-------

                        Changes in version 2.0.0                        

- The format of the output of trainSingleR() has changed and is
  no longer back-compatible.

- recompute=FALSE in trainSingleR() does nothing; all integrated
  analyses are now done with recompute=TRUE. To that end,
  combineCommonResults() is also deprecated.

- genes = "sd" and its associated options in trainSingleR() are
  no longer supported.

- first.labels is no longer reported in classifySingleR().

- Added another parallelization mechanism via num.threads= and
  C++11 threads. This should be much more memory efficient than
  using BiocParallel.

- combineRecomputedScores() will automatically handle mismatches
  in the input references by default.

[SNPRelate](/packages/SNPRelate)
---------

                       Changes in version 1.32.0                        

- update the web links

                       Changes in version 1.30.1                        

- fix a portable issue with Fortran character strings via
  USE_FC_LEN_T & FCONE

[SpatialDecon](/packages/SpatialDecon)
------------

                 Changes in version 1.8.0 (2022-08-09)                  

- updated package to logNormReg 0.4.0+

- weights were not used in logNormReg 0.3.0,
  weights are now being used as originally intended

- Bioconductor release 3.16 version

[SpatialExperiment](/packages/SpatialExperiment)
-----------------

                 Changes in version 1.7.2 (2022-10-07)                  

- support for seeing colData names with $ in RStudio

                 Changes in version 1.7.1 (2022-07-29)                  

- support for reading Space Ranger V2 outputs with read10xVisium()

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 2.3.1 (2022-10-22)                  

- In co-visualization, a new feature of bulk-to-cell mapping was
  developed, co-clustering method was simplified, and Shiny app was
  updated with options of Annotation/manual, Automatic, cell2bulk,
  bulk2cell.

- An S4 class of "SVG" was developed to store aSVG instances.

- Vignettes were updated.

[Spectra](/packages/Spectra)
-------

                         Changes in version 1.7                         

Changes in 1.7.5

- Force serial processing in some unit tests to avoid potential
failures on some Bioconductor build and check servers (under some
circumstances).

Changes in 1.7.4

- Add MsBackendMemory backend class providing a more efficient
in-memory data representation than MsBackendDataFrame.

Changes in 1.7.3

- Import spectrapply from ProtGenerics.

Changes in 1.7.2

- Fix setBackend if provided Spectra is empty.
- backendInitialize,Spectra,MsBackendDataFrame returns a Spectra
object with the full provided spectra data.

Changes in 1.7.1

- Add uniqueMsLevels function to allow more efficient,
backend-specific, implementations for retrieving unique MS levels
from a data set.

[splatter](/packages/splatter)
--------

                 Changes in version 1.22.0 (2022-10-31)                 

- 
  Fixed a bgg in BASiSSimulate() when spike.means is resampled

- 
  Fixed bugs in splatPopSimulate() with non-matching rownames and
  when sampling batches

[sSNAPPY](/packages/sSNAPPY)
-------

                 Changes in version 1.1.1 (2022-08-03)                  

- Fixed a missing column name problems in weight_ss_fc

                 Changes in version 1.1.0 (2022-06-27)                  

- Fixed a few typos

[statTarget](/packages/statTarget)
----------

                         Changes in version 2.0                         

NEW FEATURES

- New GUI
  o Mouse Hover for help information
  o .log file

- New Signal correction
  o Combat for QC-free Signal correction
  o QC-RFSC methods for metabolomics and proteomics data

- New feature slection
  o Random Forest and the Permutation based variable importance
  measures
  o new MDSplot for Random Forest
  o P-value based importance plot

- New data preprocessing
  o PQN/SUM/none normalization
  o center/none Scaling method

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.28.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- SummarizedExperiment objects now accept NAs in their rownames.
  Important notes:
  - NAs in the **rownames** are now tolerated but will cause problems
  downstream e.g. they break the rowData() getter unless
  'use.names=FALSE'
  is used.
  - NAs in the **colnames** are not and cannot be supported at the
  moment!
  Right now they break the SummarizedExperiment() constructor in an
  ugly
  way (error message not super helpful):
  > SummarizedExperiment(m)
  Error in DataFrame(x = seq_len(ncol(a1)), row.names = nms) :
  missing values in 'row.names'
  This will need to be improved.
  - At the root of these problems is the fact that at the moment
  DataFrame
  objects do NOT support NAs in their rownames.
  Bottom line: NAs in the dimnames of a SummarizedExperiment object
  should
  be avoided at all cost. One way to deal with them is to replace them
  with
  empty strings ("").
  See commit 71872cc03b7c0195fb80d1d09409243f049ebb3f.

- Small tweak to combineRows/combineCols: combineRows() and
  combineCols()
  no longer need to "fix" the dimnames that end up on the combined
  assays
  of the returned SummarizedExperiment object. So the assay dimnames
  are
  now returned as-is.
  See commit 1d6610eb168330f32433273e4fe49da173dcd33b.

[synapter](/packages/synapter)
--------

                        Changes in version 2.21                         

Changes in version 2.21.1

- Render fragmentmatching.Rmd with rmarkdown::html_document() (see
#139)

[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.9.20                        

- Adds new GeneralizedRF function to calculated information-theoretic
Generalized Robinson-Foulds distance between two dendrograms.
- Documentation for new function
- New ProtWeaver predictor based off of GeneralizedRF metric
- New internal C source code for GeneralizedRF

                       Changes in version 1.9.19                        

- Adds new DPhyloStatistic function to calculate the D-statistic for
a
binary state against a phylogeny following Fritz and Purvis (2009).
- Documentation for new function
- new internal C source code for DPhyloStatistic
- new internal C source code for random utility functions, currently
only has functions to generate random numbers

                       Changes in version 1.9.18                        

- Various internal improvements to presence/absence profile methods

                       Changes in version 1.9.17                        

- Adds new prediction algorithm GainLoss
- Adds new internal C implementation of dendrograms, significantly
faster than R dendrograms
- ProtWeaver methods Behdenna and GainLoss can now infer a species
tree when possible
- Updates Jaccard and Hamming methods to use C implementations for
distance calculation
- Adds HammingGL method to calculate Hamming distance of gain/loss
events
- Minor bugfixes to ProtWeaver methods relating to subsetting
- Updates to various man pages

                       Changes in version 1.9.16                        

- Removes flatdendrapply, function was already included in SynExtend
- minor bugfixes to ProtWeaver

                       Changes in version 1.9.15                        

- Edits to SelectByK, function can work as intended, but is still too
conservative at false positive removal.

                       Changes in version 1.9.14                        

- Adds new function flatdendrapply for more options to apply
functions
to dendrograms. Function is used in SuperTree.
- Adds new function SuperTree to construct a species tree from a set
of gene trees.
- Adds new dataset SuperTreeEx for SuperTree and flatdendrapply
examples.

                       Changes in version 1.9.13                        

- SelectByK function argument ClusterSelect switched to
ClusterScalar.
Cluster number selection now performed by fitting sum of total
within cluster sum of squares to a right hyperbola and taking the
ceiling of the half-max. Scalar allows a user to pick different
tolerances to select more, or less clusters. Plotting behavior
updated.

                       Changes in version 1.9.12                        

- simMat class now supports empty indexing (s[])
- simMat class now supports logical accession (s[c(T,F,T),])

                       Changes in version 1.9.11                        

- Added the function SelectByK that allows for quick removal of false
positive predicted pairs based on a relatively simple k-means
approach. Function is currently designed for use on the single
genome-to-genome pairwise comparison, and not on an all-vs-all many
genomes scale, though it may provide acceptable results on that
scale.

                       Changes in version 1.9.10                        

- New simMat class for dist-like similarity matrices that can be
manipulated like base matrices
- Major update to ProtWeaver internals
- All internal calls use simMat objects whenever possible to
decrease memory footprint
- Note ContextTree and ProfDCA require matrices internally
- ProtWeb objects now inherit from simMat
- ProtWeb.show and ProtWeb.print now display predictions in a more
natural way
- GetProtWebData() deprecated; ProtWeb now inherits
as.matrix.simMat and as.data.frame.simMat
- New documentation pages for simMat class
- GetProtWebData documentation page reworked into ProtWeb
documentation file.
- Fixes new bug in Method='Hamming' introduced in SynExtend 1.9.9

                        Changes in version 1.9.9                        

- Fixes minor bug in Method='Hamming'
- Moves some code around

                        Changes in version 1.9.8                        

- Major refactor to file structure of ProtWeaver to make individual
files more manageable
- Adds new documentation files for individual prediction streams of
predict.ProtWeaver

                        Changes in version 1.9.7                        

- BlockReconciliation now returns a an object of class PairSummaries.

                        Changes in version 1.9.6                        

- Fixes an error where warnings were mistakenly output to the user

                        Changes in version 1.9.5                        

- Moves platform-specific files in src/ (originally added by mistake)

                        Changes in version 1.9.4                        

- Lots of bugfixes to ResidueMI.ProtWeaver
- predict.ProtWeaver now correctly labels rows/columns with gene
names, not numbers
- predict.ProtWeaver now correctly handles Subset arguments
- predict.ProtWeaver(..., Subset=3) will correctly predict for all
pairs involving gene 3 (or for any gene x, as long as Subset is a
length 1 character or integer vector).

                        Changes in version 1.9.2                        

- Adds residue MI method to ProtWeaver
- Various bugfixes for ProtWeaver

[systemPipeShiny](/packages/systemPipeShiny)
---------------

                       Changes in version 1.7.03                        

General:

- bump up the version requirement for SPR to 2.2.0, SPRdata to 2.0.0,
spsComps to 0.3.2, drawer to 0.2.0.

Minor change:

- Since xl modal is not supported in bs3, change all modal in WF
module to l size.
- Fixed some links
- add BiocView keywords to description.

Bug fix:

- As the new version of SPR, when adding a Linewise step R code, if
the code is stored inside a variable, this is not supported.
However, SPS has to capture user input code as variable. A
workaround is used for now to fix this problem, waiting SPR to
support a better way of code in variable or quoted code.

[TargetDecoy](/packages/TargetDecoy)
-----------

                 Changes in version 1.3.4 (2022-10-21)                  

- [Feature]: add maxPoints argument to plotting function to allow
limiting the number of dots drawn in the PP-plots (#8,
@lievenclement)
- [Fix]: change y-limit in zoomed plots to the ECP of the target with
largest decoy score (#8, @lievenclement)
- Vignette: generate Figure 1 from code instead of image file (#8,
@lievenclement)
- Update gadget screenshots with new color code (#8, @lievenclement)

                 Changes in version 1.3.3 (2022-09-12)                  

- Shiny gadget: allow non-numerical variables for Score input

                 Changes in version 1.3.2 (2022-06-17)                  

- [Fix]: Pass initial argument choices to the gadget
- [Fix]: Made argument name for the log10-transformation consistent
(log10)
- evalTargetDecoysHist(): Updated default colors of targets and
decoys

                 Changes in version 1.3.1 (2022-05-20)                  

- Added a NEWS.md file to track changes to the package.
- Added Shiny gadget to interactively select variable names (#7)
- Moved decoyScoreTable() to its own file: R/decoyScoreTable.R

[TargetSearch](/packages/TargetSearch)
------------

                        Changes in version 2.0.0                        

NEW FEATURES

- Despite being a 2.0 release, the new features are minimal as
  this release contains mostly documentation improvements,
  code refactoring and clean-ups, and bug fixes.

- New methods for combining `tsLib`, `tsRim` and `tsSample` objects.
  These objects can be combined with the `c` operator. The method
  `length` reports the number of markers for objects of class `tsRim`.

- The function `quantMatrix` has been fixed and now it accepts
  three methods to generate its output. `quantmass` uses the library's
  quantification mass (QM), `maxint` takes the most abundant mass as
  QM, and `maxobs` takes the QM with the most observations. In
  addition the parameter `selmass` allows for selection of selective
  masses only if turned on.

- The function `Write.Results` has a new parameter `selmass` that is
  passed to `quantMatrix`, and the argument `quantMatrix` accepts
  the value `quantmass`, also passed to `quantMatrix`. The argument
  `prefix` changed its default value from `NA` to `NULL`.

BUG FIXES

- C code: use variables of same size (size_t vs int)

- C code: replace and update user-controlled memory interface macros.

- Baseline quantiles. Update quantiles computation with partially
  sorted
  data instead of sorted data. This yields a three-fold speed increase.

- Function `Profile`: set correct rownames. The slots of the
  `msProfile`
  object did not have correct rownames.

- Example data: fix object inconsistencies due to incorrect rownames.

- Plotting functions: they should return `invisible()`.

- Modernize the DESCRIPTION file.

- Import the whole `stats` package instead of listing each function.

- Replace instances of T/F with TRUE/FALSE in functions and vignettes.

- Add examples to functions that lacked them.

- Fix `quantMatrix`. The function has been broken for years as it did
  not work as intended. Nevertheless, the default option `quantmass`
  means that the output is now correct (formely, the options `maxint`
  and `maxobs` were simply ignored).

[TBSignatureProfiler](/packages/TBSignatureProfiler)
-------------------

                       Changes in version 1.92.0                        

Bug Fixes

- Fixed the tableAUC bootstrapped confidence interval to be the 2.5
and 97.5 percentiles instead of the 5 and 95 percentiles
- Fixed the upper CI value for the pROC/DeLong AUC CI method in the
bootstrapAUC function

Major Changes

- Changed tableAUC confidence interval default to bootstrapped
instead
of DeLong (pROC argument)

Minor Changes

- Updated the github and website introductions.
- Added Natarajan_7, Kaul_3 signatures
- Added Francisco_OD_2, Kwan_186 signatures
- Shortened some example run times

[tomoseqr](/packages/tomoseqr)
--------

                 Changes in version 1.11.1 (2022-10-25)                 

- Made the following changes:
  o Added 3D visualize function to imageViewer
  o Changed `normCount` and `normMask` options in
  `estimate3dExpressions()`
  to `normalize` option. When it is `TRUE` (default), the function
  works as if
  `normCount = "count", normMask = TRUE`. When it is `FALSE`, the
  function works as if `normCount = "none", normMask = FALSE`.

[topdownr](/packages/topdownr)
--------

                        Changes in version 1.19                         

Changes in version 1.19.2

- Avoid coercion methods that will be deprecated in Matrix 1.4-2.
- Replace stopifnot calls with if (...) stop calls.

Changes in version 1.19.1

- New version for Bioc 3.17 (devel)
- Rerun latest roxygen2.
- Remove if (getGeneric(...)) statements in AllGenerics.R do avoid
problems while reloading the package.

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.33.8                        

- Fix a bug if there is no name for the freatures.

                       Changes in version 1.33.7                        

- print the more pseudocode in browseTracks.

                       Changes in version 1.33.6                        

- Fix the bug for 'invalid color name' in browseTracks

- print the pseudocode in browseTracks

                       Changes in version 1.33.5                        

- Change the cool file open flag to "H5F_ACC_RDONLY"

                       Changes in version 1.33.4                        

- Fix a issue for plot tracks by insert 0 to the first signals.

                       Changes in version 1.33.3                        

- Add support for motif logo for lollipop plot.

                       Changes in version 1.33.2                        

- Add support for multiple operators for viewTrack.

                       Changes in version 1.33.1                        

- Fix the bug when locateScale return NA values.

[transomics2cytoscape](/packages/transomics2cytoscape)
--------------------

                        Changes in version 1.7.1                        

DATA or DOCUMENT CHANGES

- Changed the name of the network layer index in the extdata/usecase1
  and 2.

- Changed the Z height value in the extdata/usecase1 and 2.

[treeio](/packages/treeio)
------

                       Changes in version 1.21.3                        

- update as.phylo and as.treedata for data.frame object (2022-10-31,
Mon, #88)
- as.phylo() method for list (2022-09-14, Wed, #86)

                       Changes in version 1.21.2                        

- update as.treedata.pvclust method (2022-08-15, Mon, #82)
- add citation of the tree data book (2022-08-13, Sat)

                       Changes in version 1.21.1                        

- read.nextstrain.json() function to parse json tree file from
NextStrain (2022-08-03, Wed, #81)

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.15.3                        

- Up to GENCODE 42 (H.s.), M31 (M.m), and Ensembl 108

                       Changes in version 1.15.2                        

- Up to GENCODE 41 (H.s.), M30 (M.m), and Ensembl 107

[tximport](/packages/tximport)
--------

                       Changes in version 1.25.1                        

- Now the 'eds' package exports readEDS() instead of 'fishpond'.

[UCell](/packages/UCell)
-----

                        Changes in version 2.1.2                        

- New function SmoothKNN() for k-nearest neighbor smoothing of
  UCell scores. It can be applied both on SingleCellExperiment
  and Seurat objects (S3 method).

- Add two new vignettes: along with basic usage (vignette 1),
  there are now dedicated vignettes for running UCell with
  SingleCellExperiment objects (vignette 2) and Seurat objects
  (vignette 3). kNN smoothing is illustrated for both object
  types.

- Fixing a bug that prevented storing of feature ranks.

[UniProt.ws](/packages/UniProt.ws)
----------

                       Changes in version 2.38.0                        

USER VISIBLE CHANGES

- `UniProt.ws` uses the https://rest.uniprot.org/ API interface for
queries.
- `mapUniprot` is an exported function to directly map identifiers
via
UniProt and is used by the `select` method.
- `allToKeys` and `allFromKeys` provide all the available `to` and
`from`
keys for mapping identifiers
- `returnFields` provides all the possible inputs to the `columns`
argument
in the `select` and `mapUniProt` functions (ids in the `name` column)

[Uniquorn](/packages/Uniquorn)
--------

                 Changes in version 2.17.3 (2022-10-13)                 

Fixed parsing problem

- Fixd critial install error in R dpendency "Runit""

                 Changes in version 2.17.1 (2022-06-28)                 

Fixed parsing problem

- minor bugfixes

                 Changes in version 2.17.0 (2022-03-22)                 

Fixed parsing problem

- Fixed a parsing problem regarding the Cosmic and CCLE.
  databases

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.16.0                        

NEW FEATURES

- write_transfac(name.tag, altname.tag): New arguments to manually set
  the
  name and altname tags in the final TRANSFAC motifs.

MINOR CHANGES

- write_matrix(positions): Partial argument matching allowed.

- motif_tree(): Silence messages from ggtree.

BUG FIXES

- create_motif(): Don't ignore the nsites argument when generating
  random
  motifs.

- read_meme(): Correctly parse background if lines are prepended with a
  space.

- convert_type(pseudocount): Change the pseudocount within the motif
  when
  performing a type conversion and the option is set.

                       Changes in version 1.14.1                        

BUG FIXES

- read_meme(): Handle motif files with custom alphabets and no 'END'
  line in alphabet definition. Thanks to @manoloff (#24) for the bug
  report, and Spencer Nystrom for the fix (#25).

[updateObject](/packages/updateObject)
------------

                        Changes in version 1.2.0                        

- No changes in this version.

[variancePartition](/packages/variancePartition)
-----------------

                       Changes in version 1.27.6                        

- remove FMT functions

                       Changes in version 1.27.5                        

- update dependencies
- make topTable() generic to work with R 4.2.1 and Bioc 3.16

                       Changes in version 1.27.4                        

- update dependencies

                       Changes in version 1.27.3                        

- update filtering of covariates, especially for when many samples
are
dropped

                       Changes in version 1.27.1                        

- update plotPercentBars arguments

                       Changes in version 1.26.1                        

- patch to update docs

[vissE](/packages/vissE)
-----

                        Changes in version 1.6.0                        

- Added informative error messages. Improved edge case handling.

[wppi](/packages/wppi)
----

                 Changes in version 1.5.1 (2022-06-07)                  

- The workflow calculates Protein-Protein Interaction weights and
scores genes
- Database knowledge is automatically fetched from OmniPath, Gene
Ontology and Human Phenotype Ontology
- Submitted to Bioconductor

[zellkonverter](/packages/zellkonverter)
-------------

                        Changes in version 1.8.0                        

Major changes

- 
  Improve compatibility with the R *anndata* package. This
  required modifying conversion functions so that Python objects
  are explicitly converted rather than relying on automatic
  conversion.

- 
  Added support for *numpy* recarrays. This solves a
  long-standing issue and allows results from *scanpy*'s
  rank_genes_groups() function to be read.

Minor changes

- 
  The Python version is now pinned in the *anndata* v0.7.6
  environment for compatibility with changes in *basilisk*

- 
  Instatiate Python environments so they can be properly picked
  up by basilisk::configureBasiliskEnv()

- 
  Allow missing obs/var names when use_hdf5 = TRUE

- 
  Minor changes to the UI functions for compatibility with *cli*
  v3.4.0

- 
  Minor changes for compatibility with *Matrix* v1.4-2

- 
  Improvements to the UI for warnings

- 
  Updates and improvments to tests


NEWS from existing Data Experiment Packages
===================================


[curatedMetagenomicData](/packages/curatedMetagenomicData)
----------------------

                        Changes in version 3.6.0                        

- curatedMetagenomicData now contains 22,588 samples from 93 studies

- A total of 2,055 samples added since Bioconductor 3.15 (April 2022)

- Studies added since Bioconductor 3.15 (April 2022):

- BedarfJR_2017 (59 samples)

- IaniroG_2022 (165 samples)

- MetaCardis_2020_a (1,831 samples)

- Both "short" and "NCBI" row names were re-validated against NCBI
  Taxonomy

[depmap](/packages/depmap)
------

                      Changes in version 1.11.2                        

- 22Q1 data added for crispr, copyNumber, TPM, mutationCalls, metadata and 
  achilles datasets. Note: 22Q2 is the last release to follow a quarterly 
  release schedule. Future Depmap releases will follow a bi-annual release
  schedule with dataset updates every 6 months.

                      Changes in version 1.11.1                        

- EH numbers have been pdated in the latest dataset (see #78)

[imcdatasets](/packages/imcdatasets)
-----------

                 Changes in version 1.5.4 (2022-10-21)                  

- Possibility to retrieve single cell data in the SpatialExperiment
  format.

                 Changes in version 1.5.3 (2022-10-19)                  

- Updated all datasets to have consistent object formatting.

- Added dataset versioning.

- Added new dataset loading functions and deprecated the old ones.

- New vignette with guidelines for contribution and dataset formatting.

                 Changes in version 1.5.2 (2022-07-04)                  

- Updated wrapper functions

- Updated test functions

- Updated BiocViews and documentation

- Added LICENSE file

                 Changes in version 1.4.1 (2022-06-29)                  

- Added deprecation note for .onLoad functions

[NxtIRFdata](/packages/NxtIRFdata)
----------

                 Changes in version 1.3.2 (2022-10-25)                  

- Introducing SpliceWiz for Bioc 3.16 (which replaces NxtIRFcore)

                 Changes in version 1.3.1 (2022-10-07)                  

- More helpful message when ExperimentHub() fails to load or fails to
  retrieve
  Mappability files

[RforProteomics](/packages/RforProteomics)
--------------

                       Changes in version 1.35.1                        

- Remove failing code chunk.

[scpdata](/packages/scpdata)
-------

                       Changes in version 1.5.5
                       
- Added package stcker

                       Changes in version 1.5.4  
                       
- Fix: fixed character encoding errors

                       Changes in version 1.5.3   
                       
- tests: make sure each dataset is available
- docs: don't run examples. They take too long and crash the checks
                       
                       Changes in version 1.5.2      
                       
- Added CITATION
                       
                       Changes in version 1.5.1    
                       
- brunner2022: added dataset
- derks2022: added dataset
- add license statements
- leduc2022: added dataset
                       
[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.9.18                        

BUG FIXES

- Fixed a bug related to edgeR::filterByExpr() inside of
  registration_pseudobulk().

- Moved the min_ncells filtering step to registration_pseudobulk()
  rather than registration_wrapper() since you should drop low ncells
  before using edgeR::filterByExpr().

                       Changes in version 1.9.15                        

BUG FIXES

- Fixed some bugs in registration_stats_anova() in cases where we only
  had two different unique values to compute F-statistics with, when
  we need at least
  3.

- Made some parts of registration_stats_anova() and
  registration_stats_pairwise() more flexible.

- registration_model() now provides a more informative error message
  when you have an empty factor level, thus leading to a non-full rank
  model matrix.

                       Changes in version 1.9.12                        

NEW FEATURES

- Added functions for computing the modeling statistics used by the
  spatial registration process. See registration_wrapper() and related
  functions.

- Added a function for using the output of layer_stat_cor() and for
  labeling the clusters. This can help interpret the spatial
  registration results. See annotate_registered_clusters() for more
  details.

                       Changes in version 1.9.11                        

BUG FIXES

- Fixed bugs in gene_set_enrichment() for reverse = TRUE reported by
  @sparthib.

- Added a reverse option on the shiny app under the gene set
  enrichment tab, that we tested with the example spe data.

                       Changes in version 1.9.10                        

SIGNIFICANT USER-VISIBLE CHANGES

- Improved the automatic color palette selector when you switch
  discrete variables. It also now supports the ManualAnnotation
  option.

- Discrete variable (cluster) legend is no longer duplicated under the
  clusters interactive tab.

- You can now search the model test, which helps if you have lots of
  tests to choose from (this most likely occurs when you are looking
  at the pairwise results).

                        Changes in version 1.9.9                        

SIGNIFICANT USER-VISIBLE CHANGES

- Made the shiny application more memory efficient in different areas.

- Changed the default point_size from 1.25 to 2.

- Added the option to show or hide the spatial images on the grid
  panels in the shiny web application. Turn off by default since it is
  more efficient.

                        Changes in version 1.9.5                        

BUG FIXES

- Fix https://github.com/LieberInstitute/spatialLIBD/issues/41.
  Reported by @abspangler13. Now the gene selector changes
  automatically when you change the 'model results' (model type) or
  'model test' inputs. The gene selector is now only shown inside the
  'model boxplots' panel since it only affects that one.

                        Changes in version 1.9.4                        

BUG FIXES

- Fix https://github.com/LieberInstitute/spatialLIBD/issues/40.
  Reported by @Erik-D-Nelson.

                        Changes in version 1.9.3                        

BUG FIXES

- Added a more informative error message when 'stats' does not have
  ENSEMBL gene IDs as the rownames(). Reported by @abspangler13 and
  @sparthib at
  https://github.com/LieberInstitute/spatialLIBD/issues/33#issuecomment-1137544893

[STexampleData](/packages/STexampleData)
-------------

                 Changes in version 1.5.0 (2022-04-27)                  

- reformat datasets to SpatialExperiment version 1.5.3

[TargetSearchData](/packages/TargetSearchData)
----------------

                       Changes in version 1.36.0                        

NEW FEATURES

- Helper functions are added to the package. These functions are a
  simple
  interface to retrieve the data files provided by this package. The
  functions
  are mostly relevant in TargetSearch examples.


NEWS from existing Workflows
===================================


[GeoMxWorkflows](/packages/GeoMxWorkflows)
--------------

                        Changes in version 1.4.0                        

- update links

[SingscoreAMLMutations](/packages/SingscoreAMLMutations)
---------------------

                 Changes in version 1.12.0 (2022-04-05)                 

- Updated data download to download latest version of TCGA data

- Using the matched version (v36) of GENCODE for the data

Deprecated and Defunct Packages
===============================

Twenty eight software packages were removed from this release (after being deprecated
in Bioc 3.15): 
ABAEnrichment, Autotuner, BioPlex, CAnD, caOmicsV, clonotypeR, CountClust,
diffloop, GCSConnection, GCSFilesystem, GenoGAM, genphen, gprege, networkBMA,
Onassis, perturbatr, ppiStats, ProteomicsAnnotationHubData, PSICQUIC, PubScore,
Rgin, RmiR, RpsiXML, ScISI, SLGI, Sushi, tofsims, TSRchitect

Please note:  phemd and CHETAH, previously announced as deprecated in
3.15, have been updated and remain in Bioconductor.

Thirty software packages are deprecated in this release and will be removed in Bioc 3.17:
AffyCompatible, BAC, BitSeq, BrainSABER, bridge, cellTree, coexnet, conclus,
ctgGEM, CytoTree, DEComplexDisease, flowCL, flowUtils, gaia, gpart, inveRsion,
IsoGeneGUI, iteremoval, MACPET, PoTRA, rama, Rcade, RNASeqR, scAlign, scMAGeCK,
sojourner, TCGAbiolinksGUI, TDARACNE, TimeSeriesExperiment, TraRe, tspair, XCIR

Three experimental data packages were removed from this release (after being
deprecated in BioC 3.15):
DREAM4, MSstatsBioData, ppiData

Two experimental data packages are deprecated in this release and will be
removed in Bioc 3.17:
gatingMLData, RNASeqRData

One annotation packages was removed from this release (after being deprecated
in Bioc 3.15): 
MafH5.gnomAD.v3.1.1.GRCh38_3.13.1.tar.gz

No annotation packages were deprecated in this release and will be removed in
Bioc 3.17.

One workflow package was removed from this release (after being deprecated in
Bioc 3.15)
proteomics

No workflow packages were deprecated in this release and will be removed in
3.17.
