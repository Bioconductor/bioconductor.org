April 27, 2022

Bioconductors:

We are pleased to announce Bioconductor 3.15, consisting of
2140 software packages, 410 experiment data packages, 909 annotation
packages, 29 workflows and 8 books.

There are 78 new software packages, 5 new data experiment packages,
8 new annotation packages, no new workflows, no new books, and many updates and
improvements to existing packages.

Bioconductor 3.15 is compatible with R 4.2.0, and is supported on Linux,
64-bit Windows, and Intel 64-bit macOS 10.13 (High Sierra) or higher.
We do not currently support arm64 so arm64 Mac users who wish to install
Bioconductor Mac binary packages must install the Intel 64-bit build of
R available on CRAN. This release will include updated Bioconductor
[Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.15](#getting-started-with-bioconductor-315)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from new and existing software packages](#news-from-new-and-existing-software-packages)
* [NEWS from new and existing data experiment packages](#news-from-new-and-existing-data-experiment-packages)
* [NEWS from new and existing workflows](#news-from-new-and-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.15
======================================

To update to or install Bioconductor 3.15:

1. Install R 4.2.0. Bioconductor 3.15 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 78 new software packages in this release of Bioconductor.

- [APL](/packages/APL) APL is a package developed for computation of
  Association Plots (AP), a method for visualization and analysis of
  single cell transcriptomics data. The main focus of APL is the
  identification of genes characteristic for individual clusters of
  cells from input data. The package performs correspondence analysis
  (CA) and allows to identify cluster-specific genes using
  Association Plots. Additionally, APL computes the
  cluster-specificity scores for all genes which allows to rank the
  genes by their specificity for a selected cell cluster of interest.

- [ASURAT](/packages/ASURAT) ASURAT is a software for single-cell
  data analysis. Using ASURAT, one can simultaneously perform
  unsupervised clustering and biological interpretation in terms of
  cell type, disease, biological process, and signaling pathway
  activity. Inputting a single-cell RNA-seq data and knowledge-based
  databases, such as Cell Ontology, Gene Ontology, KEGG, etc., ASURAT
  transforms gene expression tables into original multivariate
  tables, termed sign-by-sample matrices (SSMs).

- [bandle](/packages/bandle) The Bandle package enables the analysis
  and visualisation of differential localisation experiments using
  mass-spectrometry data. Experimental method supported include
  dynamic LOPIT-DC, hyperLOPIT, Dynamic Organellar Maps, Dynamic PCP.
  It provides Bioconductor infrastructure to analyse these data.

- [beer](/packages/beer) BEER implements a Bayesian model for
  analyzing phage-immunoprecipitation sequencing (PhIP-seq) data.
  Given a PhIPData object, BEER returns posterior probabilities of
  enriched antibody responses, point estimates for the relative
  fold-change in comparison to negative control samples, and more.
  Additionally, BEER provides a convenient implementation for using
  edgeR to identify enriched antibody responses.

- [biodbExpasy](/packages/biodbExpasy) The biodbExpasy library
  provides access to Expasy ENZYME database, using biodb package
  framework.  It allows to retrieve entries by their accession
  number. Web services can be accessed for searching the database by
  name or comments.

- [biodbMirbase](/packages/biodbMirbase) The biodbMirbase library is
  an extension of the biodb framework package, that provides access
  to miRBase mature database. It allows to retrieve entries by their
  accession number, and run specific web services. Description: The
  biodbMirbase library provides access to the miRBase Database, using
  biodb package framework.

- [biodbNcbi](/packages/biodbNcbi) The biodbNcbi library provides
  access to the NCBI databases CCDS, Gene, Pubchem Comp and Pubchem
  Subst, using biodb package framework. It allows to retrieve entries
  by their accession number. Web services can be accessed for
  searching the database by name or mass.

- [biodbNci](/packages/biodbNci) The biodbNci library is an extension
  of the biodb framework package. It provides access to biodbNci, a
  library for connecting to the National Cancer Institute (USA)
  CACTUS Database. It allows to retrieve entries by their accession
  number, and run specific web services.

- [BOBaFIT](/packages/BOBaFIT) This package provides a method to
  refit and correct the diploid region in copy number profiles. It
  uses a clustering algorithm to identify pathology-specific normal
  (diploid) chromosomes and then use their copy number signal to
  refit the whole profile.  The package is composed by three
  functions: DRrefit (the main function), ComputeNormalChromosome and
  PlotCluster.

- [borealis](/packages/borealis) Borealis is an R library performing
  outlier analysis for count-based bisulfite sequencing data. It
  detectes outlier methylated CpG sites from bisulfite sequencing
  (BS-seq). The core of Borealis is modeling Beta-Binomial
  distributions. This can be useful for rare disease diagnoses.

- [CBEA](/packages/CBEA) This package implements CBEA, a method to
  perform set-based analysis for microbiome relative abundance data.
  This approach constructs a competitive balance between taxa within
  the set and remainder taxa per sample. More details can be found in
  the Nguyen et al. 2021+ manuscript. Additionally, this package adds
  support functions to help users perform taxa-set enrichment
  analyses using existing gene set analysis methods.  In the future
  we hope to also provide curated knowledge driven taxa sets.

- [cellxgenedp](/packages/cellxgenedp) The cellxgene data portal
  (https://cellxgene.cziscience.com/) provides a graphical user
  interface to collections of single-cell sequence data processed in
  standard ways to 'count matrix' summaries. The cellxgenedp package
  provides an alternative, R-based inteface, allowind data discovery,
  viewing, and downloading.

- [CNVMetrics](/packages/CNVMetrics) The CNVMetrics package
  calculates similarity metrics to facilitate copy number variant
  comparison among samples and/or methods. Similarity metrics can be
  employed to compare CNV profiles of genetically unrelated samples
  as well as those with a common genetic background. Some metrics are
  based on the shared amplified/deleted regions while other metrics
  rely on the level of amplification/deletion. The data type used as
  input is a plain text file containing the genomic position of the
  copy number variations, as well as the status and/or the log2 ratio
  values. Finally, a visualization tool is provided to explore
  resulting metrics.

- [cogeqc](/packages/cogeqc) cogeqc aims to facilitate systematic
  quality checks on standard comparative genomics analyses to help
  researchers detect issues and select the most suitable parameters
  for each data set. cogeqc can be used to asses: i. genome assembly
  quality with BUSCOs; ii. orthogroup inference using a protein
  domain-based approach and; iii. synteny detection using synteny
  network properties. There are also data visualization functions to
  explore QC summary statistics.

- [comapr](/packages/comapr) comapr detects crossover intervals for
  single gametes from their haplotype states sequences and stores the
  crossovers in GRanges object. The genetic distances can then be
  calculated via the mapping functions using estimated crossover
  rates for maker intervals. Visualisation functions for plotting
  interval-based genetic map or cumulative genetic distances are
  implemented, which help reveal the variation of crossovers
  landscapes across the genome and across individuals.

- [coMethDMR](/packages/coMethDMR) coMethDMR identifies genomic
  regions associated with continuous phenotypes by optimally
  leverages covariations among CpGs within predefined genomic
  regions. Instead of testing all CpGs within a genomic region,
  coMethDMR carries out an additional step that selects co-methylated
  sub-regions first without using any outcome information. Next,
  coMethDMR tests association between methylation within the
  sub-region and continuous phenotype using a random coefficient
  mixed effects model, which models both variations between CpG sites
  within the region and differential methylation simultaneously.

- [CompoundDb](/packages/CompoundDb) CompoundDb provides
  functionality to create and use (chemical) compound annotation
  databases from a variety of different sources such as LipidMaps,
  HMDB, ChEBI or MassBank. The database format allows to store in
  addition MS/MS spectra along with compound information. The package
  provides also a backend for Bioconductor's Spectra package and
  allows thus to match experimetal MS/MS spectra against MS/MS
  spectra in the database. Databases can be stored in SQLite format
  and are thus portable.

- [COTAN](/packages/COTAN) Statistical and computational method to
  analyze the co-expression of gene pairs at single cell level. It
  provides the foundation for single-cell gene interactome analysis.
  The basic idea is studying the zero UMI counts' distribution
  instead of focusing on positive counts; this is done with a
  generalized contingency tables framework. COTAN can effectively
  assess the correlated or anti-correlated expression of gene pairs.
  It provides a numerical index related to the correlation and an
  approximate p-value for the associated independence test. COTAN can
  also evaluate whether single genes are differentially expressed,
  scoring them with a newly defined global differentiation index.
  Moreover, this approach provides ways to plot and cluster genes
  according to their co-expression pattern with other genes,
  effectively helping the study of gene interactions and becoming a
  new tool to identify cell-identity marker genes.

- [crisprBase](/packages/crisprBase) Provides S4 classes for general
  nucleases, CRISPR nucleases and base editors. Several
  CRISPR-specific genome arithmetic functions are implemented to help
  extract genomic coordinates of spacer and protospacer sequences.
  Commonly-used CRISPR nuclease objects are provided that can be
  readily used in other packages. Both DNA- and RNA-targeting
  nucleases are supported.

- [crisprBowtie](/packages/crisprBowtie) Provides a user-friendly
  interface to map on-targets and off-targets of CRISPR gRNA spacer
  sequences using bowtie. The alignment is fast, and can be performed
  using either commonly-used or custom CRISPR nucleases. The
  alignment can work with any reference or custom genomes. Both DNA-
  and RNA-targeting nucleases are supported.

- [crisprBwa](/packages/crisprBwa) Provides a user-friendly interface
  to map on-targets and off-targets of CRISPR gRNA spacer sequences
  using bwa. The alignment is fast, and can be performed using either
  commonly-used or custom CRISPR nucleases. The alignment can work
  with any reference or custom genomes. Currently not supported on
  Windows machines.

- [crisprScore](/packages/crisprScore) Provides R wrappers of several
  on-target and off-target scoring methods for CRISPR guide RNAs
  (gRNAs). The following nucleases are supported: SpCas9, AsCas12a,
  enAsCas12a, and RfxCas13d (CasRx). The available on-target cutting
  efficiency scoring methods are RuleSet1, Azimuth, DeepHF, DeepCpf1,
  enPAM+GB, and CRISPRscan. Both the CFD and MIT scoring methods are
  available for off-target specificity prediction. The package also
  provides a Lindel-derived score to predict the probability of a
  gRNA to produce indels inducing a frameshift for the Cas9 nuclease.
  Note that DeepHF, DeepCpf1 and enPAM+GB are not available on
  Windows machines.

- [cytoMEM](/packages/cytoMEM) MEM, Marker Enrichment Modeling,
  automatically generates and displays quantitative labels for cell
  populations that have been identified from single-cell data. The
  input for MEM is a dataset that has pre-clustered or pre-gated
  populations with cells in rows and features in columns. Labels
  convey a list of measured features and the features' levels of
  relative enrichment on each population. MEM can be applied to a
  wide variety of data types and can compare between MEM labels from
  flow cytometry, mass cytometry, single cell RNA-seq, and spectral
  flow cytometry using RMSD.

- [DepInfeR](/packages/DepInfeR) DepInfeR integrates two
  experimentally accessible input data matrices: the drug sensitivity
  profiles of cancer cell lines or primary tumors ex-vivo (X), and
  the drug affinities of a set of proteins (Y), to infer a matrix of
  molecular protein dependencies of the cancers (ß). DepInfeR
  deconvolutes the protein inhibition effect on the viability
  phenotype by using regularized multivariate linear regression. It
  assigns a “dependence coefficient” to each protein and each sample,
  and therefore could be used to gain a causal and accurate
  understanding of functional consequences of genomic aberrations in
  a heterogeneous disease, as well as to guide the choice of
  pharmacological intervention for a specific cancer type, sub-type,
  or an individual patient. For more information, please read out
  preprint on bioRxiv: https://doi.org/10.1101/2022.01.11.475864.

- [DifferentialRegulation](/packages/DifferentialRegulation)
  DifferentialRegulation is a method for detecting differentially
  regulated genes between two groups of samples (e.g., healthy vs.
  disease, or treated vs. untreated samples), by targeting
  differences in the balance of spliced and unspliced mRNA
  abundances, obtained from single-cell RNA-sequencing (scRNA-seq)
  data. DifferentialRegulation accounts for the sample-to-sample
  variability, and embeds multiple samples in a Bayesian hierarchical
  model. In particular, when reads are compatible with multiple genes
  or multiple splicing versions of a gene (unspliced spliced or
  ambiguous), the method allocates these multi-mapping reads to the
  gene of origin and their splicing version. Parameters are inferred
  via Markov chain Monte Carlo (MCMC) techniques
  (Metropolis-within-Gibbs).

- [EpiCompare](/packages/EpiCompare) EpiCompare is used to compare
  and analyse epigenetic datasets for quality control and
  benchmarking purposes. The package outputs an HTML report
  consisting of three sections: (1. General metrics) Metrics on peaks
  (percentage of blacklisted and non-standard peaks, and peak widths)
  and fragments (duplication rate) of samples, (2. Peak overlap)
  Percentage and statistical significance of overlapping and
  non-overlapping peaks. Also includes upset plot and (3. Functional
  annotation) functional annotation (ChromHMM, ChIPseeker and
  enrichment analysis) of peaks. Also includes peak enrichment around
  TSS.

- [epimutacions](/packages/epimutacions) The package includes some
  statistical outlier detection methods for epimutations detection in
  DNA methylation data. The methods included in the package are
  MANOVA, Multivariate linear models, isolation forest, robust
  mahalanobis distance, quantile and beta. The methods compare a case
  sample with a suspected disease against a reference panel (composed
  of healthy individuals) to identify epimutations in the given case
  sample. It also contains functions to annotate and visualize the
  identified epimutations.

- [extraChIPs](/packages/extraChIPs) This package builds on existing
  tools and adds some simple but extremely useful capabilities for
  working with ChIP-Seq data. The focus is on detecting differential
  binding windows/regions. One set of functions focusses on
  set-operations retaining mcols for GRanges objects, whilst another
  group of functions are to aid visualisatino of results. Coercion to
  tibble objects is also included.

- [fastreeR](/packages/fastreeR) Calculate distances, build
  phylogenetic trees or perform hierarchical clustering between the
  samples of a VCF or FASTA file. Functions are implemented in Java
  and called via rJava. Parallel implementation that operates
  directly on the VCF or FASTA file for fast execution.

- [GBScleanR](/packages/GBScleanR) GBScleanR is a package for quality
  check, filtering, and error correction of genotype data derived
  from next generation sequcener (NGS) based genotyping platforms.
  GBScleanR takes Variant Call Format (VCF) file as input. The main
  function of this package is `estGeno()` which estimates the true
  genotypes of samples from given read counts for genotype markers
  using a hidden Markov model with incorporating uneven observation
  ratio of allelic reads. This implementation gives robust genotype
  estimation even in noisy genotype data usually observed in
  Genotyping-By-Sequnencing (GBS) and similar methods, e.g. RADseq.
  The current implementation accepts genotype data of a diploid
  population at any generation of multi-parental cross, e.g.
  biparental F2 from inbred parents, biparental F2 from outbred
  parents, and 8-way recombinant inbred lines (8-way RILs) which can
  be refered to as MAGIC population.

- [GenomicInteractionNodes](/packages/GenomicInteractionNodes) The
  GenomicInteractionNodes package can import interactions from bedpe
  file and define the interaction nodes, the genomic interaction
  sites with multiple interaction loops. The interaction nodes is a
  binding platform regulates one or multiple genes. The detected
  interaction nodes will be annotated for downstream validation.

- [GenProSeq](/packages/GenProSeq) Generative modeling for protein
  engineering is key to solving fundamental problems in synthetic
  biology, medicine, and material science. Machine learning has
  enabled us to generate useful protein sequences on a variety of
  scales. Generative models are machine learning methods which seek
  to model the distribution underlying the data, allowing for the
  generation of novel samples with similar properties to those on
  which the model was trained. Generative models of proteins can
  learn biologically meaningful representations helpful for a variety
  of downstream tasks. Furthermore, they can learn to generate
  protein sequences that have not been observed before and to assign
  higher probability to protein sequences that satisfy desired
  criteria. In this package, common deep generative models for
  protein sequences, such as variational autoencoder (VAE),
  generative adversarial networks (GAN), and autoregressive models
  are available. In the VAE and GAN, the Word2vec is used for
  embedding. The transformer encoder is applied to protein sequences
  for the autoregressive model.

- [ggmanh](/packages/ggmanh) Manhattan plot and QQ Plot are commonly
  used to visualize the end result of Genome Wide Association Study.
  The "ggmanh" package aims to keep the generation of these plots
  simple while maintaining customizability. Main functions include
  manhattan_plot, qqunif, and thinPoints.

- [GRaNIE](/packages/GRaNIE) Genetic variants associated with
  diseases often affect non-coding regions, thus likely having a
  regulatory role. To understand the effects of genetic variants in
  these regulatory regions, identifying genes that are modulated by
  specific regulatory elements (REs) is crucial. The effect of gene
  regulatory elements, such as enhancers, is often cell-type
  specific, likely because the combinations of transcription factors
  (TFs) that are regulating a given enhancer have celltype specific
  activity. This TF activity can be quantified with existing tools
  such as diffTF and captures differences in binding of a TF in open
  chromatin regions. Collectively, this forms a gene regulatory
  network (GRN) with cell-type and data-specific TF-RE and RE-gene
  links. Here, we reconstruct such a GRN using bulk RNAseq and open
  chromatin (e.g., using ATACseq or ChIPseq for open chromatin marks)
  and optionally TF activity data. Our network contains different
  types of links, connecting TFs to regulatory elements, the latter
  of which is connected to genes in the vicinity or within the same
  chromatin domain (TAD). We use a statistical framework to assign
  empirical FDRs and weights to all links using a permutation-based
  approach.

- [hermes](/packages/hermes) Provides classes and functions for
  quality control, filtering, normalization and differential
  expression analysis of pre-processed RNA-seq data. Data can be
  imported from `SummarizedExperiment` as well as `matrix` objects
  and can be annotated from BioMart. Filtering for genes without too
  low expression or containing required annotations, as well as
  filtering for samples with sufficient correlation to other samples
  or total number of reads is supported. The standard normalization
  methods including `cpm`, `rpkm` and `tpm` can be used, and `DESeq2`
  as well as `voom` differential expression analyses are available.

- [lineagespot](/packages/lineagespot) Lineagespot is a framework
  written in R, and aims to identify SARS-CoV-2 related mutations
  based on a single (or a list) of variant(s) file(s) (i.e., variant
  calling format). The method can facilitate the detection of
  SARS-CoV-2 lineages in wastewater samples using next generation
  sequencing, and attempts to infer the potential distribution of the
  SARS-CoV-2 lineages.

- [LinTInd](/packages/LinTInd) When we combine gene-editing
  technology and sequencing technology, we need to reconstruct a
  lineage tree from alleles generated and calculate the similarity
  between each pair of groups. FindIndel() and IndelForm() function
  will help you align each read to reference sequence and generate
  scar form strings respectively. IndelIdents() function will help
  you to define a scar form for each cell or read. IndelPlot()
  function will help you to visualize the distribution of deletion
  and insertion. TagProcess() function will help you to extract
  indels for each cell or read. TagDist() function will help you to
  calculate the similarity between each pair of groups across the
  indwells they contain. BuildTree() function will help you to
  reconstruct a tree. PlotTree() function will help you to visualize
  the tree.

- [Macarron](/packages/Macarron) Macarron is a workflow for the
  prioritization of potentially bioactive metabolites from
  metabolomics experiments. Prioritization integrates strengths of
  evidences of bioactivity such as covariation with a known
  metabolite, abundance relative to a known metabolite and
  association with an environmental or phenotypic indicator of
  bioactivity. Broadly, the workflow consists of stratified
  clustering of metabolic spectral features which co-vary in
  abundance in a condition, transfer of functional annotations,
  estimation of relative abundance and differential abundance
  analysis to identify associations between features and
  phenotype/condition.

- [MBECS](/packages/MBECS) The Microbiome Batch Effect Correction
  Suite (MBECS) provides a set of functions to evaluate and mitigate
  unwated noise due to processing in batches. To that end it
  incorporates a host of batch correcting algorithms (BECA) from
  various packages. In addition it offers a correction and reporting
  pipeline that provides a preliminary look at the characteristics of
  a data-set before and after correcting for batch effects.

- [mbOmic](/packages/mbOmic) The mbOmic package contains a set of
  analysis functions for microbiomics and metabolomics data, designed
  to analyze the inter-omic correlation between microbiology and
  metabolites. Integrative analysis of the microbiome and metabolome
  is the aim of mbOmic. Additionally, the identification of
  enterotype using the gut microbiota abundance is
  preliminaryimplemented.

- [MetaboAnnotation](/packages/MetaboAnnotation) High level functions
  to assist in annotation of (metabolomics) data sets. These include
  functions to perform simple tentative annotations based on mass
  matching but also functions to consider m/z and retention times for
  annotation of LC-MS features given that respective reference values
  are available. In addition, the function provides high-level
  functions to simplify matching of LC-MS/MS spectra against spectral
  libraries and objects and functionality to represent and manage
  such matched data.

- [MobilityTransformR](/packages/MobilityTransformR)
  MobilityTransformR collects a tool set for effective mobility scale
  transformation of CE-MS/MS data in order to increase
  reproducibility. It provides functionality to determine the
  migration times from mobility markers that have been added to the
  analysis and performs the transformation based on these markers.
  MobilityTransformR supports the conversion of numeric vectors,
  Spectra-objects, and MSnOnDiskExp.

- [Motif2Site](/packages/Motif2Site) Detect binding sites using
  motifs IUPAC sequence or bed coordinates and ChIP-seq experiments
  in bed or bam format. Combine/compare binding sites across
  experiments, tissues, or conditions. All normalization and
  differential steps are done using TMM-GLM method. Signal
  decomposition is done by setting motifs as the centers of the
  mixture of normal distribution curves.

- [MSA2dist](/packages/MSA2dist) MSA2dist calculates pairwise
  distances between all sequences of a DNAStringSet or a AAStringSet
  using a custom score matrix and conducts codon based analysis. It
  uses scoring matrices to be used in these pairwise distance
  calcualtions which can be adapted to any scoring for DNA or AA
  characters. E.g. by using literal distances MSA2dist calcualtes
  pairwise IUPAC distances.

- [MsBackendMsp](/packages/MsBackendMsp) Mass spectrometry (MS) data
  backend supporting import and handling of MS/MS spectra from NIST
  MSP Format (msp) files. Import of data from files with different
  MSP *flavours* is supported. Objects from this package add support
  for MSP files to Bioconductor's Spectra package. This package is
  thus not supposed to be used without the Spectra package that
  provides a complete infrastructure for MS data handling.

- [MuData](/packages/MuData) Save MultiAssayExperiments to h5mu files
  supported by muon and mudata. Muon is a Python framework for
  multimodal omics data analysis. It uses an HDF5-based format for
  data storage.

- [netZooR](/packages/netZooR) PANDA(Passing Attributes between
  Networks for Data Assimilation) is a message-passing model to
  reconstruction gene regulatory network. It integrates multiple
  sources of biological data, including protein-protein interaction
  data, gene expression data, and sequence motif information to
  reconstruct genome-wide, condition-specific regulatory
  networks.[(Glass et al. 2013)]. LIONESS(Linear Interpolation to
  Obtain Network Estimates for Single Samples) is a method to
  estimate sample-specific regulatory networks by applying linear
  interpolation to the predictions made by existing aggregate network
  inference approaches. CONDOR(COmplex Network Description Of
  Regulators)is a bipartite community structure analysis tool of
  biological networks, especially eQTL networks, including a method
  for scoring nodes based on their modularity contribution.[(Platig
  et al. 2016). ALPACA(ALtered Partitions Across Community
  Architectures) is a method for comparing two genome-scale networks
  derived from different phenotypic states to identify
  condition-specific modules.[(Padi and Quackenbush 2018)]. This
  package integrates pypanda--the Python implementation of PANDA and
  LIONESS(https://github.com/davidvi/pypanda),the R implementation of
  CONDOR(https://github.com/jplatig/condor) and the R implementation
  of ALPACA (https://github.com/meghapadi/ALPACA) into one workflow.
  Each tool can be call in this package by one function, and the
  relevant output could be accessible in current R session for
  downstream analysis.

- [nnSVG](/packages/nnSVG) Method for scalable identification of
  spatially variable genes (SVGs) in spatially-resolved
  transcriptomics data. The method is based on nearest-neighbor
  Gaussian processes and uses the BRISC algorithm for model fitting
  and parameter estimation. Allows identification and ranking of SVGs
  with flexible length scales across a tissue slide or within spatial
  domains defined by covariates. Scales linearly with the number of
  spatial locations and can be applied to datasets containing
  thousands or more spatial locations.

- [OGRE](/packages/OGRE) OGRE calculates overlap between user defined
  genomic region datasets. Any regions can be supplied i.e. genes,
  SNPs, or reads from sequencing experiments. Key numbers help
  analyse the extend of overlaps which can also be visualized at a
  genomic level.

- [omicsViewer](/packages/omicsViewer) omicsViewer visualizes
  ExpressionSet (or SummarizedExperiment) in an interactive way. The
  omicsViewer has a separate back- and front-end. In the back-end,
  users need to prepare an ExpressionSet that contains all the
  necessary information for the downstream data interpretation. Some
  extra requirements on the headers of phenotype data or feature data
  are imposed so that the provided information can be clearly
  recognized by the front-end, at the same time, keep a minimum
  modification on the existing ExpressionSet object. The pure
  dependency on R/Bioconductor guarantees maximum flexibility in the
  statistical analysis in the back-end. Once the ExpressionSet is
  prepared, it can be visualized using the front-end, implemented by
  shiny and plotly. Both features and samples could be selected from
  (data) tables or graphs (scatter plot/heatmap). Different types of
  analyses, such as enrichment analysis (using Bioconductor package
  fgsea or fisher's exact test) and STRING network analysis, will be
  performed on the fly and the results are visualized simultaneously.
  When a subset of samples and a phenotype variable is selected, a
  significance test on means (t-test or ranked based test; when
  phenotype variable is quantitative) or test of independence
  (chi-square or fisher’s exact test; when phenotype data is
  categorical) will be performed to test the association between the
  phenotype of interest with the selected samples. Additionally,
  other analyses can be easily added as extra shiny modules.
  Therefore, omicsViewer will greatly facilitate data exploration,
  many different hypotheses can be explored in a short time without
  the need for knowledge of R. In addition, the resulting data could
  be easily shared using a shiny server. Otherwise, a standalone
  version of omicsViewer together with designated omics data could be
  easily created by integrating it with portable R, which can be
  shared with collaborators or submitted as supplementary data
  together with a manuscript.

- [ompBAM](/packages/ompBAM) This packages provides C++ header files
  for developers wishing to create R packages that processes BAM
  files. ompBAM automates file access, memory management, and
  handling of multiple threads 'behind the scenes', so developers can
  focus on creating domain-specific functionality. The included
  vignette contains detailed documentation of this API, including
  quick-start instructions to create a new ompBAM-based package, and
  step-by-step explanation of the functionality behind the example
  packaged included within ompBAM.

- [PanomiR](/packages/PanomiR) PanomiR is a package to detect miRNAs
  that target groups of pathways from gene expression data. This
  package provides functionality for generating pathway activity
  profiles, determining differentially activated pathways between
  user-specified conditions, determining clusters of pathways via the
  PCxN package, and generating miRNAs targeting clusters of pathways.
  These function can be used separately or sequentially to analyze
  RNA-Seq data.

- [pareg](/packages/pareg) Compute pathway enrichment scores while
  accounting for term-term relations. This package uses a regularized
  multiple linear regression to regress differential expression
  p-values obtained from multi-condition experiments on a pathway
  membership matrix. By doing so, it is able to incorporate
  additional biological knowledge into the enrichment analysis and to
  estimate pathway enrichment scores more robustly.

- [protGear](/packages/protGear) A generic three-step pre-processing
  package for protein microarray data. This package contains
  different data pre-processing procedures to allow comparison of
  their performance.These steps are background correction, the
  coefficient of variation (CV) based filtering, batch correction and
  normalization.

- [PSMatch](/packages/PSMatch) The PSMatch package helps proteomics
  practitioners to load, handle and manage Peptide Spectrum Matches.
  It provides functions to model peptide-protein relations as
  adjacency matrices and connected components, visualise these as
  graphs and make informed decision about shared peptide filtering.
  The package also provides functions to calculate and visualise MS2
  fragment ions.

- [qmtools](/packages/qmtools) The qmtools (quantitative metabolomics
  tools) package provides basic tools for processing quantitative
  metabolomics data with the standard SummarizedExperiment class.
  This includes functions for imputation, normalization, feature
  filtering, feature clustering, dimension-reduction, and
  visualization to help users prepare data for statistical analysis.
  Several functions in this package could also be used in other types
  of omics data.

- [qsvaR](/packages/qsvaR) The qsvaR package contains functions for
  removing the effect of degration in rna-seq data from postmortem
  brain tissue. The package is equipped to help users generate
  principal components associated with degradation. The components
  can be used in differential expression analysis to remove the
  effects of degradation.

- [RAREsim](/packages/RAREsim) Haplotype simulations of rare variant
  genetic data that emulates real data can be performed with RAREsim.
  RAREsim uses the expected number of variants in MAC bins - either
  as provided by default parameters or estimated from target data -
  and an abundance of rare variants as simulated HAPGEN2 to
  probabilistically prune variants. RAREsim produces haplotypes that
  emulate real sequencing data with respect to the total number of
  variants, allele frequency spectrum, haplotype structure, and
  variant annotation.

- [Rbwa](/packages/Rbwa) Provides an R wrapper for BWA alignment
  algorithms. Both BWA-backtrack and BWA-MEM are available.
  Convenience function to build a BWA index from a reference genome
  is also provided. Currently not supported for Windows machines.

- [RCX](/packages/RCX) Create, handle, validate, visualize and
  convert networks in the Cytoscape exchange (CX) format to standard
  data types and objects. The package also provides conversion to and
  from objects of iGraph and graphNEL. The CX format is also used by
  the NDEx platform, a online commons for biological networks, and
  the network visualization software Cytocape.

- [rgoslin](/packages/rgoslin) The R implementation for the Grammar
  of Succint Lipid Nomenclature parses different short hand notation
  dialects for lipid names. It normalizes them to a standard name. It
  further provides calculated monoisotopic masses and sum formulas
  for each successfully parsed lipid name and supplements it with
  LIPID MAPS Category and Class information. Also, the structural
  level and further structural details about the head group, fatty
  acyls and functional groups are returned, where applicable.

- [rifi](/packages/rifi) 'rifi' analyses data from rifampicin time
  series created by microarray or RNAseq. 'rifi' is a transcriptome
  data analysis tool for the holistic identification of transcription
  and decay associated processes. The decay constants and the delay
  of the onset of decay is fitted for each probe/bin. Subsequently,
  probes/bins of equal properties are combined into segments by
  dynamic programming, independent of a existing genome annotation.
  This allows to detect transcript segments of different stability or
  transcriptional events within one annotated gene. In addition to
  the classic decay constant/half-life analysis, 'rifi' detects
  processing sites, transcription pausing sites, internal
  transcription start sites in operons, sites of partial
  transcription termination in operons, identifies areas of likely
  transcriptional interference by the collision mechanism and gives
  an estimate of the transcription velocity. All data are integrated
  to give an estimate of continous transcriptional units, i.e.
  operons. Comprehensive output tables and visualizations of the full
  genome result and the individual fits for all probes/bins are
  produced.

- [RolDE](/packages/RolDE) RolDE detects longitudinal differential
  expression between two conditions in noisy high-troughput data.
  Suitable even for data with a moderate amount of missing
  values.RolDE is a composite method, consisting of three independent
  modules with different approaches to detecting longitudinal
  differential expression. The combination of these diverse modules
  allows RolDE to robustly detect varying differences in longitudinal
  trends and expression levels in diverse data types and experimental
  settings.

- [rprimer](/packages/rprimer) Functions, workflow, and a Shiny
  application for visualizing sequence conservation and designing
  degenerate primers, probes, and (RT)-(q/d)PCR assays from a
  multiple DNA sequence alignment. The results can be presented in
  data frame format and visualized as dashboard-like plots. For more
  information, please see the package vignette.

- [sccomp](/packages/sccomp) A robust and outlier-aware method for
  testing differential tissue composition from single-cell data. This
  model can infer changes in tissue composition and heterogeneity,
  and can produce realistic data simulations based on any existing
  dataset. This model can also transfer knowledge from a large set of
  integrated datasets to increase accuracy further.

- [seqArchR](/packages/seqArchR) seqArchR enables unsupervised
  discovery of _de novo_ clusters with characteristic sequence
  architectures characterized by position-specific motifs or
  composition of stretches of nucleotides, e.g., CG-richness.
  seqArchR does _not_ require any specifications w.r.t. the number of
  clusters, the length of any individual motifs, or the distance
  between motifs if and when they occur in pairs/groups; it directly
  detects them from the data. seqArchR uses non-negative matrix
  factorization (NMF) as its backbone, and employs a chunking-based
  iterative procedure that enables processing of large sequence
  collections efficiently. Wrapper functions are provided for
  visualizing cluster architectures as sequence logos.

- [single](/packages/single) Accurate consensus sequence from
  nanopore reads of a DNA gene library. SINGLe corrects for
  systematic errors in nanopore sequencing reads of gene libraries
  and it retrieves true consensus sequences of variants identified by
  a barcode, needing only a few reads per variant. More information
  in preprint doi: https://doi.org/10.1101/2020.03.25.007146.

- [SPOTlight](/packages/SPOTlight) `SPOTlight`provides a method to
  deconvolute spatial transcriptomics spots using a seeded NMF
  approach along with visualization tools to assess the results.
  Spatially resolved gene expression profiles are key to understand
  tissue organization and function. However, novel spatial
  transcriptomics (ST) profiling techniques lack single-cell
  resolution and require a combination with single-cell RNA
  sequencing (scRNA-seq) information to deconvolute the spatially
  indexed datasets. Leveraging the strengths of both data types, we
  developed SPOTlight, a computational tool that enables the
  integration of ST with scRNA-seq data to infer the location of cell
  types and states within a complex tissue. SPOTlight is centered
  around a seeded non-negative matrix factorization (NMF) regression,
  initialized using cell-type marker genes and non-negative least
  squares (NNLS) to subsequently deconvolute ST capture locations
  (spots).

- [sSNAPPY](/packages/sSNAPPY) A single sample pathway pertrubation
  testing methods for RNA-seq data. The method propagate changes in
  gene expression down gene-set topologies to compute single-sample
  directional pathway perturbation scores that reflect potentail
  directions of changes.Perturbation scores can be used to test
  significance of pathway perturbation at both individual-sample and
  treatment levels.

- [standR](/packages/standR) standR is an user-friendly R package
  providing functions to assist conducting good-practice analysis of
  Nanostring's GeoMX DSP data. All functions in the package are built
  based on the SpatialExperiment object, allowing integration into
  various spatial transcriptomics-related packages from Bioconductor.
  standR allows data inspection, quality control, normalization,
  batch correction and evaluation with informative visualizations.

- [STdeconvolve](/packages/STdeconvolve) STdeconvolve as an
  unsupervised, reference-free approach to infer latent cell-type
  proportions and transcriptional profiles within multi-cellular
  spatially-resolved pixels from spatial transcriptomics (ST)
  datasets. STdeconvolve builds on latent Dirichlet allocation (LDA),
  a generative statistical model commonly used in natural language
  processing for discovering latent topics in collections of
  documents. In the context of natural language processing, given a
  count matrix of words in documents, LDA infers the distribution of
  words for each topic and the distribution of topics in each
  document. In the context of ST data, given a count matrix of gene
  expression in multi-cellular ST pixels, STdeconvolve applies LDA to
  infer the putative transcriptional profile for each cell-type and
  the proportional representation of each cell-type in each
  multi-cellular ST pixel.

- [TEKRABber](/packages/TEKRABber) TEKRABber is made to provide a
  user-friendly pipeline for comparing orthologs and transposable
  elements (TEs) between two species. It considers the orthology
  confidence between two species from BioMart to normalize expression
  counts and detect differentially expressed orthologs/TEs. Then it
  provides one to one correlation analysis for desired orthologs and
  TEs. There is also an app function to have a first insight on the
  result. Users can prepare orthologs/TEs RNA-seq expression data by
  their own preference to run TEKRABber following the data structure
  mentioned in the vignettes.

- [terraTCGAdata](/packages/terraTCGAdata) Leverage the existing open
  access TCGA data on Terra with well-established Bioconductor
  infrastructure. Make use of the Terra data model without learning
  its complexities. With a few functions, you can copy / download and
  generate a MultiAssayExperiment from the TCGA example workspaces
  provided by Terra.

- [tomoseqr](/packages/tomoseqr) `tomoseqr` is an R package for
  analyzing Tomo-seq data. Tomo-seq is a genome-wide RNA tomography
  method that combines combining high-throughput RNA sequencing with
  cryosectioning for spatially resolved transcriptomics. `tomoseqr`
  reconstructs 3D expression patterns from tomo-seq data and
  visualizes the reconstructed 3D expression patterns.

- [TREG](/packages/TREG) RNA abundance and cell size parameters could
  improve RNA-seq deconvolution algorithms to more accurately
  estimate cell type proportions given the different cell type
  transcription activity levels. A Total RNA Expression Gene (TREG)
  can facilitate estimating total RNA content using single molecule
  fluorescent in situ hybridization (smFISH). We developed a
  data-driven approach using a measure of expression invariance to
  find candidate TREGs in postmortem human brain single nucleus
  RNA-seq. This R package implements the method for identifying
  candidate TREGs from snRNA-seq data.

- [UCell](/packages/UCell) UCell is a package for evaluating gene
  signatures in single-cell datasets. UCell signature scores, based
  on the Mann-Whitney U statistic, are robust to dataset size and
  heterogeneity, and their calculation demands less computing time
  and memory than other available methods, enabling the processing of
  large datasets in a few minutes even on machines with limited
  computing power. UCell can be applied to any single-cell data
  matrix, and includes functions to directly interact with
  SingleCellExperiment and Seurat objects.

- [updateObject](/packages/updateObject) A set of tools built around
  updateObject() to work with old serialized S4 instances. The
  package is primarily useful to package maintainers who want to
  update the serialized S4 instances included in their package. This
  is still work-in-progress.

- [xcore](/packages/xcore) xcore is an R package for transcription
  factor activity modeling based on known molecular signatures and
  user's gene expression data. Accompanying xcoredata package
  provides a collection of molecular signatures, constructed from
  publicly available ChiP-seq experiments. xcore use ridge regression
  to model changes in expression as a linear combination of molecular
  signatures and find their unknown activities. Obtained, estimates
  can be further tested for significance to select molecular
  signatures with the highest predicted effect on the observed
  expression changes.


New Data Experiment Packages
=====================

There are 5 new data experiment packages in this release of Bioconductor.

- [crisprScoreData](/packages/crisprScoreData) Provides an interface
  to access pre-trained models for on-target and off-target gRNA
  activity prediction algorithms implemented in the crisprScore
  package. Pre-trained model data are stored in the ExperimentHub
  database. Users should consider using the crisprScore package
  directly to use and load the pre-trained models.

- [epimutacionsData](/packages/epimutacionsData) This package
  includes the data necessary to run functions and examples in
  epimutacions package. Collection of DNA methylation data. The
  package contains 2 datasets: (1) Control ( GEO: GSE104812), (GEO:
  GSE97362) case samples; and (2) reference panel (GEO: GSE127824).
  It also contains candidate regions to be epimutations in 450k
  methylation arrays.

- [healthyControlsPresenceChecker](/packages/healthyControlsPresenceChecker)
  A function that reads in the GEO accession code of a gene
  expression dataset, retrieves its data from GEO, and checks if data
  of healthy controls are present in the dataset. It returns true if
  healthy controls data are found, and false otherwise.  GEO: Gene
  Expression Omnibus. ID: identifier code. The GEO datasets are
  downloaded from the URL <https://ftp.ncbi.nlm.nih.gov/geo/series/>.

- [VectraPolarisData](/packages/VectraPolarisData) Provides two
  multiplex imaging datasets collected on Vectra instruments at the
  University of Colorado Anschutz Medical Campus. Data are provided
  as a Spatial Experiment objects. Data is provided in tabular form
  and has been segmented and phenotyped using Inform software. Raw
  .tiff files are not included.

- [xcoredata](/packages/xcoredata) Provides data to use with xcore
  package.


New Annotation Packages
=====================

There are 8 new annotation packages in this release of Bioconductor.

- [BSgenome.Cjacchus.UCSC.calJac4](/packages/BSgenome.Cjacchus.UCSC.calJac4)
  Full genome sequences for Callithrix jacchus (Marmoset) as provided
  by UCSC (calJac4, May 2020) and wrapped in a BSgenome object.

- [BSgenome.CneoformansVarGrubiiKN99.NCBI.ASM221672v1](/packages/BSgenome.CneoformansVarGrubiiKN99.NCBI.ASM221672v1)
  Full genome sequences for Cryptococcus neoformans var. grubii KN99
  (assembly ASM221672v1 assembly accession GCA_002216725.1).

- [BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0](/packages/BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
  The T2T-CHM13v2.0 assembly (accession GCA_009914755.4), as
  submitted to NCBI by the T2T Consortium, and wrapped in a
  BSgenome object. Companion paper: "The complete sequence of a
  human genome" by Nurk S, Koren S, Rhie A, Rautiainen M, et al.
  Science, 2022.

- [BSgenome.Mfascicularis.NCBI.6.0](/packages/BSgenome.Mfascicularis.NCBI.6.0)
  Full genome sequences for Macaca fascicularis (Crab-eating
  macaque) as provided by NCBI (assembly Macaca_fascicularis_6.0,
  assembly accession GCA_011100615.1) and stored in Biostrings objects.

- [JASPAR2022](/packages/JASPAR2022)
  JASPAR is an open-access database containing manually
  curated, non-redundant transcription factor (TF) binding
  profiles for TFs across six taxonomic groups. In this 9th
  release, we expanded the CORE collection with 341 new profiles
  (148 for plants, 101 for vertebrates, 85 for urochordates, and
  7 for insects), which corresponds to a 19% expansion over the
  previous release. To search thisdatabases, please use the
  package TFBSTools (>= 1.31.2).

- [MafH5.gnomAD.v3.1.2.GRCh38](/packages/MafH5.gnomAD.v3.1.2.GRCh38)
  Store minor allele frequency data from the Genome
  Aggregation Database (gnomAD version 3.1.2) for the human
  genome version GRCh38.

- [SNPlocs.Hsapiens.dbSNP155.GRCh38](/packages/SNPlocs.Hsapiens.dbSNP155.GRCh38)
  SNP locations and alleles for Homo sapiens extracted from NCBI dbSNP
  Build 155. The 948,979,291 SNPs in this package were extracted from
  the RefSNP JSON files for chromosomes 1-22, X, Y, and MT, located at
  https://ftp.ncbi.nih.gov/snp/latest_release/JSON/ (these files were
  created by NCBI on May 25, 2021). These SNPs can be "injected" in
  BSgenome.Hsapiens.NCBI.GRCh38 or BSgenome.Hsapiens.UCSC.hg38.

- [UCSCRepeatMasker](/packages/UCSCRepeatMasker)
  Store UCSC RepeatMasker AnnotationHub resource metadata.
  Provide provenance and citation information for UCSC
  RepeatMasker AnnotationHub resources. Illustrate in a vignette
  how to access those resources.

New Workflow Packages
=====================

There are no new workflow packages.


New Books
=====================

There are no new online books.


NEWS from new and existing Software Packages
===================================

[ADaCGH2](/packages/ADaCGH2)
-------

                 Changes in version 2.35.1 (2022-04-19)                 

- Added required code from ffbase (as ffbase no longer in CRAN)
  (ffbase authors added to authors' list). ffbase no longer
  a dependency.

[affxparser](/packages/affxparser)
----------

                 Changes in version 1.67.1 (2022-03-23)                 

SIGNIFICANT CHANGES

- This packages requires R (>= 4.0.0) when build on MS Windows. This is
  due to the added support for UCRT on MS Windows, which is required
  for
  the upcoming R 4.2.0.

SOFTWARE QUALITY

- Updates to build package from source on MS Windows with UCRT.  Thanks
  to Tomas Kalibera for the contribution.

- Now registering native routines - apparently never happened before.

                 Changes in version 1.67.0 (2021-10-27)                 

- The version number was bumped for the Bioconductor devel version,
  which is
  now BioC 3.15 for R-devel.

[airpart](/packages/airpart)
-------

                        Changes in version 1.3.4                        

- Paper published in Bioinformatics...

[alevinQC](/packages/alevinQC)
--------

                       Changes in version 1.11.1                        

- Added support for reading alevin-fry output

[AlpsNMR](/packages/AlpsNMR)
-------

                 Changes in version 3.5.1 (2022-04-07)                  

- plot_interactive now accepts an overwrite argument to avoid asking
the user interactively
- Improve nmr_detect_peaks_tune_snr to tune the SNR threshold with
the
right other parameters
- Documentation improvements
- Split Peak_detection page into smaller and more specific pages
- Let the user choose how code is parallellized, as suggested by
BiocParallel documentation.
- Replace furr/future parallellization loops with BiocParallel.
Provides a warning in case a future::plan() has been set.
- Demote Imports to Suggests: SummarizedExpriment, S4Vectors,
ggrepel,
GGally
- Remove dependencies: tidyselect, assertthat, plyr, furrr
- Add download_MTBLS242() function to help download the data for the
tutorial
- Skip mixOmics test if affected by
https://github.com/mixOmicsTeam/mixOmics/pull/199
- Fix auto setting of the baseline threshold for the peak detection

[AnnotationHub](/packages/AnnotationHub)
-------------

                        Changes in version 3.3.0                        

NEW FEATURES

- (3.3.6) Serialized S4 hub resources (AnnotationHub and ExperimentHub)
  are now passed thru updateObject() at load-time

USER-VISIBLE MODIFICATIONS

- (3.3.8) Add instructions for creating a hub shared across multiple
  users

- (3.3.2) Remove TESTING option as only needed to expose devel orgdb

- (3.3.2) Change filter for orgdbs. orgdbs at release time will be
  stamped
  with devel (to be release) and then manually have biocversion added
  for the
  upcoming new devel. Filtering then based on biocversion number. This
  will
  expose devel orgdbs as soon as generated

BUG CORRECTION

- (3.1.9) Fix broken test. Identical not appropriate. all present
  appropriate

[AnnotationHubData](/packages/AnnotationHubData)
-----------------

                       Changes in version 1.25.0                        

SIGNIFICANT UPDATES

- 1.25.7 Update recipes to upload to azure. NonStandardOrgDb release
  recipe
  updated

- 1.25.6 Update recipes to upload to azure. TwoBit ensembl and release
  recipes for standard TxDb and OrgDb updated

NEW FEATURES

- 1.25.5 Add helper function to upload to azure

MODIFICATIONS

- 1.25.2 Changed makeAnnotationHubMetadata to point to Azure instead of
  AWS

[AnVIL](/packages/AnVIL)
-----

                        Changes in version 1.8.0                        

NEW FEATURES

- (v 1.7.4) add avworkflow_configuration_*() functions for
manipulating workflow configurations, and a vignette describing use.

- (v 1.7.5) add avdata_import() to import 'REFERENCE DATA' and 'OTHER
DATA' tables.

- (v 1.7.9) export repository_stats() to summarize binary package
availability.

USER VISIBLE CHANGES

- (v 1.7.4) Deprecate avworkflow_configuration(),
avworkflow_import_configuration().

- (v 1.7.4) Update Dockstore md5sum.

- (v 1.7.5) avdata() is re-implemented to more faithfully report only
'REFERENCE DATA' and 'OTHER DATA' workspace attributes; previously,
other attributes such as the description and tags (from the
workspace landing page) were also reported.

BUG FIXES

- (v 1.7.4) avworkflow_files() and avworkflow_localize() do not fail
when the workflow has produced no files.

- (v 1.7.6) improve handling of authentication token for gcloud
utilities.

                       Changes in version 1.7.13                        

BUG FIXES

- Correct gcloud_project() when user environment variable set.
https://github.com/Bioconductor/AnVIL/pull/52

                        Changes in version 1.6.6                        

BUG FIXES

- Correct gsutil_pipe() argument mis-match, see
https://support.bioconductor.org/p/9141780/

[APAlyzer](/packages/APAlyzer)
--------

                 Changes in version 1.9.4 (2022-02-28)                  

- Fixed typos in vignettes

                 Changes in version 1.9.3 (2022-02-27)                  

- Fixed the issues caused by the lastest ensembl GRCm39 GTF.

- Fixed typos in vignettes

                 Changes in version 1.9.2 (2022-01-02)                  

- Added supporting of multi-condition design in `APAdiff`.

                 Changes in version 1.9.1 (2021-11-28)                  

- Added `MultiTest` parameter for APAdiff.

- Improved the performance and fixed issues in PAS2GEF.

- Update documents and authorships.

[aroma.light](/packages/aroma.light)
-----------

                 Changes in version 3.25.1 (2022-04-21)                 

BUG FIXES

- Fixed partial argument name used in iwpca().

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.19.3                        

- fix the email typo.

                       Changes in version 1.19.2                        

- Change the method to fix the issue that NA is generated for
conservation scores when call gscores.

                       Changes in version 1.19.1                        

- Fix the issue that NA is generated for conservation scores when
call
gscores.

[atena](/packages/atena)
-----

                 Changes in version 1.2.0 (2022-04-21)                  

USER VISIBLE CHANGES

- Higher accuracy in TE quantification for TEtranscripts and Telescope
  methods.

- Improved EM step running time.

- New AnnotationHub resource has been added: UCSCRepeatMasker.

- Implemented function to retrieve and parse TE annotations.

[AUCell](/packages/AUCell)
------

                        Changes in version 1.17                         

- New function: AUCell_run()

- Support for DelayedArray and Sparse matrices

[BANDITS](/packages/BANDITS)
-------

                       Changes in version 1.11.1                        

- "prior_precision" function made faster (prior computation based on a
  sub-set of the genes only)

[bandle](/packages/bandle)
------

                       Changes in version 1.0.0                       

version 0.99.8

- minor bioc changes

version 0.99.7

- speed up package compiling
- fixing compiling issues

version 0.99.6

- fix pointer to self in compiled code

version 0.99.5

- remove timeout issues

version 0.99.4

- fix bioconductor review suggestions

version 0.99.2

- bioconductor stable version initialized

version 0.0

- initial package commit

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.7.1 (2021-11-03)                  

- Add FixNu argument to fix scaling parameters to scran normalisation
  factors
  when WithSpikes=FALSE

- Exclude features with absolute estimates of fold change/difference
  less
  than the DE threshold when performing differential expression
  analysis
  with BASiCS_TestDE.

[BayesSpace](/packages/BayesSpace)
----------

                        Changes in version 1.5.1                        

Minor improvements and fixes

- Update readVisium to work properly when feature names include
whitespaces
- Fix rounding of w_i sum in MCMC
- Add details to vignette related to spatialEnhance

                        Changes in version 1.5.0                        

New Bioconductor devel (3.15)

- Version numbering change with Bioconductor version bump

[BEclear](/packages/BEclear)
-------

                 Changes in version 2.11.1 (2022-03-31)                 

- replaced dependency outliers with dixonTest

[beer](/packages/beer)
----

                 Changes in version 0.99.8 (2022-02-14)                 

- Fixed minor bugs with outdated arguments in the vignette.

                 Changes in version 0.99.7 (2022-02-14)                 

- Changed `bp.param` argument to be called `BPPARAM` to be consistent
  with BiocParallel arguments.

                 Changes in version 0.99.6 (2022-02-09)                 

- Corrected `future` reference in `README` to `BiocParallel`

- Changed ordering of arguments to put those without defaults first.

                 Changes in version 0.99.3 (2022-01-31)                 

- Documentation
  * Corrected typos, formatting errors, and missing links.
  * Modified `DESCRIPTION` to include `URL` and remove `Collates`
  * Removed directions for installing `PhIPData` in the `README` and
  Vignette.
  * Added package citations for `edgeR`.
  * Added package man page.
  * Additional details have been added to the documentation for
  `phipseq_model.bugs`.
  * Added citation.

- Code
  * Standardized function names to use `camelCase`.
  * Standardized function arguments to use `dot.case`.
  * Though `edgeROne()` and `brewOne()` are not meant to be directly
  interfaced by the user, the two function names have been changed to
  reflect that they are still exported.
  * Renamed `edgeR()` to `runEdgeR()` to be more descriptive.
  * Corrected `paste0()` that should have been `file.path()`.
  * Converted to `BiocParallel` from `future` for parallelization
  support.

                 Changes in version 0.99.0 (2022-01-13)                 

- Submitted to Bioconductor

[benchdamic](/packages/benchdamic)
----------

                 Changes in version 1.1.1 (2021-03-21)                  

- Corrected a bug in DA_ALDEx2 function

- Updated DA_ALDEx2 function manual

- Added more authors in the citation

                 Changes in version 1.1.0 (2021-10-26)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_14 branch

[BindingSiteFinder](/packages/BindingSiteFinder)
-----------------

                        Changes in version 1.1.2                        

- coverageOverRanges() can be allowed to produce NAs for uneven ranges
  in the
  output

                        Changes in version 1.1.1                        

- coverageOverRanges() matches the order of input ranges and output
  matrix for
  options merge_all_replicates and merge_replicates_per_condition

[bioCancer](/packages/bioCancer)
---------

                       Changes in version 1.23.03                       

- host cgdsr

                       Changes in version 1.23.02                       

- Update NEWS

                       Changes in version 1.23.01                       

- get cgdsr from github

[BiocCheck](/packages/BiocCheck)
---------

                        Changes in version 1.32                         

NEW FEATURES

- Add package metadata to main report for easier diagnostics

- `<pkgname>.BiocCheck` folder, created above the package folder,
  includes
  the full report and NAMESPACE suggestions, if available.

- Add check to find any stray `<pkgname>.BiocCheck` folders

- Update doc links and recommendations for additional information in
  report

- Update `BiocCheck` report to be more brief by only noting the
  conditions;
  details are included in the full report

BUG FIXES AND MINOR IMPROVEMENTS

- Initialize default verbose value (FALSE) for internal reference
  object

- Flag only hidden '.RData' files as bad files and allow 'myData.RData'
  (@hpages, #155)

- Improve internal handling of condition messages with unified
  mechanism

- Internal improvements to `BiocCheck` mechanism: export `.BiocCheck`
  object
  which contains all conditions, log list, and method for writing to
  JSON

- Update to changes in R 4.2 `--no-echo` flag

- Make use of `lib.loc` to helper functions that install and load the
  checked package

- (1.31.36) Reduce function length count slightly by removing empty
  lines.

- (1.31.35) Restricted files in `inst` will be flagged with a `WARNING`
  instead of an `ERROR`

- (1.31.32) Account for S3 print methods when checking for `cat` usage

- (1.31.31) Single package imports in the NAMESPACE were breaking the
  code
  to get all package imports.

- (1.31.29) Include other import fields from NAMESPACE file when
  checking
  consistency between imports in DESCRIPTION/NAMESPACE.

- (1.31.27) Update and clean up unit tests.

- (1.31.26) Improve load test for the package being checked.

- (1.31.25) Exclude GitHub URLs that end in HTML from external data
  check.

- (1.31.23) Internal updates to the `require` and `library` check.

- (1.31.22) Remove old code related to running `BiocCheck` on the
  command
  line and update `BiocCheck` documentation.

- (1.31.21) Remove redundant `=` from message to avoid `=` assignment.

- (1.31.20) Add line feed to "Checking function lengths..." message

- (1.31.18) Packages should not download files when loaded or attached.

- (1.31.17) Using '=' for assignment should be avoided and '<-' should
  be
  used instead for clarity and legibility.

- (1.31.16) Note the use of `cat` and `print` outside of show methods.

- (1.31.15) Check for pinned package versions in the `DESCRIPTION` file
  denoted by the use of `==`.

- (1.31.14) Enhancements to internal helper functions and
  `BiocCheckGitClone`

- (1.31.13) Revert move to new package checks. Update Bioc-devel
  mailing
  list check to fail early when not in BBS environment.

- (1.31.12) Move Bioc-devel mailing list and support site registration
  checks to new package checks.

- (1.31.10) Various internal improvements to `BiocCheck` and the
  identification of the package directory and name.

- (1.31.6) Use a more reliable approach to identify package name from
  the
  `DESCRIPTION` file.

- (1.31.5) Fixed bug in the case where the `VignetteBuilder` field in
  the
  a package's `DESCRIPTION` has more than one listed.

- (1.31.3) Add `BioCbooks` repository url to
  `checkIsPackageNameAlreadyInUse`, `VIEWS` file is pending.

- (1.31.2) Fix logical length > 1 error in `checkImportSuggestions`
  (@vjcitn, #141)

- (1.31.1) Simplify check for function lengths; remove excessive dots.

[BiocFileCache](/packages/BiocFileCache)
-------------

                         Changes in version 2.3                         

ENHANCEMENT

- (2.3.4) Add instructions for a shared cache across multiple users of
  a system

- (2.3.2) Add direct SQL calls for certain retrieval functions to speed
  up
  access time. This will speed up the bfcquery function as well as any
  function tha utilized the underlying .sql_get_field, .get_all_rids or
  .get_all_web_rids.

- (2.3.1) Add @LTLA solution for making bfcrpath thread safe

[BioCor](/packages/BioCor)
------

                       Changes in version 1.19.1                        

- Added alternative text to plots on vignettes.
- Reactivated again the comparison with GOSemSim

[BiocParallel](/packages/BiocParallel)
------------

                        Changes in version 1.30                         

USER VISIBLE CHANGES

- (v 1.29.1) Report first remote error in its entirety.
  https://github.com/Bioconductor/BiocParallel/issues/165

- (v 1.29.4) Add bpresult() (extract result vector from return value
  of tryCatch(bplapply(...))) and allow direct use of
  tryCatch(bplapply(...))
  return value as arugment to bplapply(BPREDO= ...). Closes #157

- (v 1.29.8) The default timeout for worker computation changes
  from 30 days to .Machine$integer.max (no timeout), allowing for
  performance improvements when not set.

- (v 1.29.11) The timeout for establishing a socket connection is
  set to getOption("timeout") (default 60 seconds).

- (v 1.29.15) Check for and report failed attempts to open SnowParam
  ports.

- (v 1.29.18) add bpfallback= option to control use of `lapply()`
  (fallback) when 0 or 1 workers are available.

- (v 1.29.19) add bpexportvariables= option to automatically
  export global variables, or variables found in packages on the
  search path, in user-provided `FUN=` functions.

BUG FIXES

- (v 1.29.2) Fix regression in use of debug() with SerialParam.
  https://github.com/Bioconductor/BiocParallel/issues/128

- (v 1.29.3) Fix regression in progress bar display with bplapply().
  https://github.com/Bioconductor/BiocParallel/issues/172

- (v 1.29.5) Fix default seed generation when user has non-default
  generator. https://github.com/Bioconductor/BiocParallel/pull/176

- (v 1.29.9) Fix validity when workers, specified as character(),
  are more numerous than (non-zero) tasks.
  https://github.com/Bioconductor/BiocParallel/pull/181

[BiocStyle](/packages/BiocStyle)
---------

                       Changes in version 2.24.0                        

BUG FIXES

- Fixed incompatibility with recent versions of the longtable latex
  package
  when producing a PDF vignette.  This bug manifested with the message
  'LaTeX Error: Missing \begin{document}.' when attempting to knit.
  (https://github.com/Bioconductor/BiocStyle/issues/89)

[biocViews](/packages/biocViews)
---------

                       Changes in version 1.63.0                        

BUG FIX

- (1.63.1) Fix bug with Authors@R parsing for making VIEWS

[biodb](/packages/biodb)
-----

                 Changes in version 1.3.3 (2022-04-02)                  

- Explain how to use custom CSV file in vignette.

                 Changes in version 1.3.2 (2022-03-11)                  

- Correct ext pkg upgrade: do not generate C++ example files again.

- Change default vignette name to package name.

- When upgrading a extension package, do not generate C++ files again.

- Accept an unknown total in Progress class.

- isSearchableByField() accepts field.type param now.

                 Changes in version 1.3.1 (2021-12-10)                  

- Remove custom cache folder setting in vignette.

[biodbChebi](/packages/biodbChebi)
----------

                 Changes in version 1.1.1 (2022-03-15)                  

- Change vignette name. Set author name.

- Remove code already in biodb package.

- Upgrade maintenance files.

[biodbExpasy](/packages/biodbExpasy)
-----------

                 Changes in version 0.99.0 (2022-03-15)                 

- Submitted to Bioconductor

[biodbHmdb](/packages/biodbHmdb)
---------

                 Changes in version 1.1.2 (2022-03-13)                  

- Change vignette name.

- Update HMDB extract zip file used for vignette and testing.

[biodbKegg](/packages/biodbKegg)
---------

                 Changes in version 1.1.1 (2022-03-17)                  

- Upgrade maintenance files.

- Correct instantiation example.

[biodbLipidmaps](/packages/biodbLipidmaps)
--------------

                 Changes in version 1.1.1 (2022-03-17)                  

- Update maintenance files.

- Add log output to ouput file during check.

[biodbMirbase](/packages/biodbMirbase)
------------

                 Changes in version 0.99.1 (2022-04-01)                 

- Remove getNbEntries() example in vignette. No such web service exists
  in
  miRBase.

                 Changes in version 0.99.0 (2022-03-22)                 

- Submitted to Bioconductor

[biodbNcbi](/packages/biodbNcbi)
---------

                 Changes in version 0.99.6 (2022-03-24)                 

- Corrected example inside vignette.

                 Changes in version 0.99.0 (2022-03-17)                 

- Submitted to Bioconductor

[biodbNci](/packages/biodbNci)
--------

                 Changes in version 0.99.0 (2022-03-22)                 

- Submitted to Bioconductor

[biodbUniprot](/packages/biodbUniprot)
------------

                 Changes in version 1.1.1 (2022-03-17)                  

- Update maintenance files.

- Correct instantiation in example.

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.52.0                        

BUG FIXES

- Stop reporting message about the use of https when using useEnsembl()
  with a 'version' argument.

- Use virtualSchemaName provided by a Mart, rather than simply
  "default".
  This caused issues with the Ensembl Plants Mart.

[BiRewire](/packages/BiRewire)
--------

                       Changes in version 3.27.5                        

- Update support tsne-> Rtsne

                       Changes in version 3.27.1                        

- Added a constraint flag for excluding positive and negative arcs
  between two nodes

[BridgeDbR](/packages/BridgeDbR)
---------

                        Changes in version 2.5.0                        

NEW FEATURES

- Updated to BridgeDb 3.0.13

BUG FIXES

- Small typo fix in a code example in the Vignette

[bugsigdbr](/packages/bugsigdbr)
---------

                        Changes in version 1.2.0                        

- importBugSigDB accepts Zenodo DOIs and Github hashes to obtain
defined data releases or devel to obtain the latest version (new
argument version).
- Ontology-based queries for experimental factors and body sites (new
functions getOntology and subsetByOntologies)
- Compilation of meta-signatures from individual signatures for one
body site and one condition at a time, weighted by sample size
(new function getMetaSignatures)

[CAGEr](/packages/CAGEr)
-----

                        Changes in version 2.2.0                        

BUG FIXES

- Restore the `correctSystematicG` option in `getCTSS()`.  See #61.

- Restore object class consistency in if / else statement in private
  function
  `bam2CTSS`.  Fixes #49.

- Restore proper CTSS conversion from BAM files (bug introduced in
  v1.34.0)
  while fixing issue #36.

- Ensure Tag Clusters have a Seqinfo. Fixes #63.

[CAMERA](/packages/CAMERA)
------

                       Changes in version 1.51.1                        

USER VISIBLE CHANGES

- Move to mzML for example data and vignette following the removal of
  mzData support in mzR

[canceR](/packages/canceR)
------

                       Changes in version 1.29.03                       

- host cgdsr functions

                       Changes in version 1.29.01                       

- get cgdsr from github

[cbaf](/packages/cbaf)
----

                 Changes in version 1.18.0 (2022-04-24)                 

New Features

- Packge now uses cBioPortalData package to communicate with cBioPorta,
  as cgdsr package is deprecated.

- Package now supports RNA-seq data with z-scores relative to normal
  samples.

- Terms updated.

- Minor improvemets

[CBEA](/packages/CBEA)
----

                       Changes in version 0.99.3                        

NEW FEATURES

- Created an output type object (CBEAout). This is an S3 type object
that is essentially a list that incorporates the final score matrix
as well as other diagnostic details.

SIGNIFICANT USER-VISIBLE CHANGES

- Due to the new feature above, now instead of getting a tibble,
users
would have to extract the scores out either using the provided
function or use a custom approach based on the list-of-list format
of the CBEAout objects.
- Implemented tidy and glance methods to deal with CBEAout objects

BUG FIXES

None

                       Changes in version 0.99.2                        

NEW FEATURES

- Added an option (parametric) to specify whether the null is
estimated parametrically or via pure permutation. To support this,
an option (n_perm) was also added to specify the number of
permutations. A warning will be added if parametric is FALSE but
n_boot is small (< 100)
- Added option (parallel_backend) to specify the parallel backend of
the loop using BiocParallel
- Argument control now allow for a special slot titled fix_comp that
can specify which component of the two-component mixture
distribution to fix during the adjustment process.

SIGNIFICANT USER-VISIBLE CHANGES

- Combined raw argument with output argument for the CBEA function to
specify returning raw CBEA scores (without any distribution fitting
and transformation).
- New and revamped vignettes
- Significant reduction in dependencies, including removing native
support for phyloseq

BUG FIXES

None

                       Changes in version 0.99.1                        

NEW FEATURES

None

SIGNIFICANT USER-VISIBLE CHANGES

None

BUG FIXES

- Removed .Rproj file to conform with Bioconductor error

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- Added a complete functionality to perform CBEA from scratch with
bundled data set

SIGNIFICANT USER-VISIBLE CHANGES

None

BUG FIXES

None

[cBioPortalData](/packages/cBioPortalData)
--------------

                        Changes in version 2.8.0                        

New features

- Auth token string or file can now be included in the cBioPortal
function.
- The check_build argument can be set to FALSE for alternative APIs,
e.g., KidsFirst, when using cBioPortalData
- queryGeneTable translates gene IDs ('hugoGeneSymbols' <>
'entrezGeneIds') via the API service
- getDataByGenes supersedes getDataByGenePanel
- getStudies() replaces data('studiesTable') to discover study IDs

Bug fixes and minor improvements

- Fixed issue where the by argument was not passed to getDataByGenes
in internal calls
- Add names to metadata elements that originate from GISTIC datasets.

[cbpManager](/packages/cbpManager)
----------

                        Changes in version 1.3.2                        

- Oncotree functionality is now avaliable also in the sample tab.
- Major modifications of the Mutation tab. Mutation data is now
editable.

[celda](/packages/celda)
-----

                 Changes in version 1.11.0 (2022-03-31)                 

- Improvments to decontX vignette
- Added ability to subsample to speed up perplexity calculations
- Added ability to use batch parameter with the raw matrix in decontX

[ceRNAnetsim](/packages/ceRNAnetsim)
-----------

                 Changes in version 1.7.3 (2022-15-03)                  

- New functions, find_affected_nodes and find_targeting nodes

- Fixed test files

- Fixed documentations

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.29.6                        

- update jianhong's email

                       Changes in version 3.29.5                        

- update the pipeline documentation to avoid the error when converting
  the
  GRanges object to data.frame.

                       Changes in version 3.29.3                        

- handle the error "non standard bed file."

                       Changes in version 3.29.2                        

- fix the error: trying to generate an object from a virtual class
  "DataFrame".

                       Changes in version 3.29.1                        

- add subGroupComparison parameter to getEnrichedGO and
  getEnrichedPATH.

[ChIPQC](/packages/ChIPQC)
------

                       Changes in version 1.31.2                        

- bugfixes

- handle larger plots

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.32.2                        

- update vignette

                       Changes in version 1.31.4                        

- readPeakFile now supports .broadPeak and .gappedPeak files
(2021-12-17, Fri, #173)

                       Changes in version 1.31.3                        

- bug fixed of determining promoter region in minus strand
(2021-12-16, Thu, #172)

                       Changes in version 1.31.1                        

- bug fixed to take strand information (2021-11-10, Wed, #167)

[chromstaR](/packages/chromstaR)
---------

                       Changes in version 1.21.1                        

BUG FIXES

- Fixed ENSEMBL host for rnorvegicus_gene_ensembl. Using default host
  now instead of may2012.archive.ensembl.org.

[CIMICE](/packages/CIMICE)
------

                        Changes in version 1.3.1                        

Overview:

- Removed dependency with relations

[ClassifyR](/packages/ClassifyR)
---------

                        Changes in version 3.0.0                        

- 
  Now supports survival models and their evaluation, in addition
  to existing classification functionality.

- 
  Cross-validation no longer requires specific annotations like
  data set name and classifier name. Now, the user can specify
  any characteristics they want and use these as variables to
  group by or change line appearances. Also, characteristics like
  feature selection name and classifier name are automatically
  filled in from an internal table.

- 
  Ease of use greatly inproved by crossValidate function which
  allows specification of classifiers by a single keyword.
  Previously, parameter objects such as SelectParams and
  TrainParams had to be explicitly specified, making it
  challenging for users not familar with S4 object-oriented
  programming.

- 
  Basic multi-omics data integration functionality available via
  crossValidate which allows combination of different tables.
  Pre-validation and PCA dimensionality techniques provide a fair
  way to compare high-dimensional omics data with low-dimensional
  clinical data. Also, it is possible to simply concatenate all
  data tables.

- 
  Model-agnostic variable importance calculated by training when
  leaving out one selected variable at a time. Turned off by
  default as it substantially increases run time. See
  doImportance parameter of ModellingParams for more details.

- 
  Parameters specifying the cross-validation procedure and data
  modelling formalised as CrossValParams and ModellingParams
  classes.

- 
  Feature selection can now be done either based a on
  resubstitution metric (i.e. train and test on the training
  data) or a cross-validation metric (i.e. split the training
  data into training and testing partitions to tune the selected
  features). all feature selection functions have been converted
  into feature ranking functions, because the selection procedure
  is a feature of cross-validation.

- 
  All function and class documentation coverted from manually
  written Rd files to Roxygen format.

- 
  Human Reference Interactome (binary experimental PPI) included
  in bundled data for pairs-based classification. See ?HuRI for
  more details.

- 
  Performance plots can now do either box plots or violin plots.
  Box plot remains the default style.

[cleanUpdTSeq](/packages/cleanUpdTSeq)
------------

                       Changes in version 1.33.2                        

- Merge the codes by Haibo

                       Changes in version 1.33.1                        

- fix a bug for one line GRanges reported by Haibo

[clusterProfiler](/packages/clusterProfiler)
---------------

                        Changes in version 4.3.4                        

- fix enrichGO , gseGO and groupGO when keyType = 'SYMBOL' &&
readable=TRUE(2022-4-9, Sat)

                        Changes in version 4.3.3                        

- parse GAF file to prepare GO annotation data (esp for proteomic
study) (2022-03-08, Tue, #397, #418, #421, #442)
- bug fixed in compareCluster() (2022-01-27, Thu, #424)

                        Changes in version 4.3.2                        

- bug fixed in extract_params() (2022-01-12, Wed, #392, @amcdavid)
- make simplify() works for gseGO() in compareCluster()
- support formula interface for GSEA methods in compareCluster()
(2022-01-04, Tue, @altairwei, #416)

                        Changes in version 4.3.1                        

- compareCluster() supports GSEA algorithm (2021-12-11, Sat)
- update error message of download.KEGG.Path() and
download.KEGG.Module()(2021-11-21, Sun)
- update simplify() function to support ont = ALL (2021-10-27, Wed)

[clustifyr](/packages/clustifyr)
---------

                 Changes in version 1.7.3 (2022-03-09)                  

- `vec_out` option for directly getting classification results as a
  vector, to be inserted into other metadata/workflow

- Maintainer change

[CNVMetrics](/packages/CNVMetrics)
----------

                        Changes in version 0.1.6                        

NEW FEATURES

- calculateLog2ratioMetric() method enables log2 ratio metric
calculation using similar workflow than
calculateOverlapRegionsMetric() method.

SIGNIFICANT USER-VISIBLE CHANGES

- plotMetric() replaces plotOverlapMetric() method. The new method
can
plot all metrics (state call metrics and log2 ratio metrics).
- Vignette section 'Workflow for metrics calculated using the level
of
amplification/deletion' is complete.
- New citing section in README and vignette refering to published
F1000Research poster
(http://www.doi.org/10.7490/f1000research.1118704.1).
- Instead of calculating distance, log2 ratio metrics are calculated
distance-based metrics (1/(1+distance)).

BUG FIXES

- None

                        Changes in version 0.1.4                        

NEW FEATURES

- None

SIGNIFICANT USER-VISIBLE CHANGES

- New website https://krasnitzlab.github.io/CNVMetrics/index.html
associated to package.
- Vignette section 'Workflow for metrics calculated using CNV status
calls' is complete.

BUG FIXES

- None

                        Changes in version 0.1.3                        

NEW FEATURES

- None

SIGNIFICANT USER-VISIBLE CHANGES

- plotOneOverlapMetric() method has a new argument silent=TRUE so
that
the plot is not drawn by default.

BUG FIXES

- plotOneOverlapMetric() method now uses sample distance for
clustering as default clustering method.

                        Changes in version 0.1.2                        

NEW FEATURES

- plotOneOverlapMetric() method enables plotting result of
overlapping
metric calculation.

SIGNIFICANT USER-VISIBLE CHANGES

- calculateOverlapRegionsMetric() method changed to
calculateOverlapMetric().

BUG FIXES

- None

                        Changes in version 0.1.1                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- calculateOverlapRegionsMetric() method enables calculation of
similarity metrics using overlapping amplified/deleted regions.

SIGNIFICANT USER-VISIBLE CHANGES

- None.

BUG FIXES

- None

[cola](/packages/cola)
----

                        Changes in version 2.1.1                        

- add hierarchical consensus partitioning

                        Changes in version 2.0.1                        

- `predict_classes()`: support svm/random forest methods for class
  label prediction.

[coMET](/packages/coMET)
-----

                 Changes in version 1.27.2 (2022-04-18)                 

- Update ENSEMBL functions for https

- Remove DNase_UCSC functions

- Add information about font.family for plotTrack and other Gviz
  functions

- Add showtext library in coMET.Rnw vignette

[coMethDMR](/packages/coMethDMR)
---------

                       Changes in version 0.99.9                        

We have been working on bug fixes and formatting/documentation
changes
as requested during the Bioconductor review. See these requests here:
https://github.com/Bioconductor/Contributions/issues/2064

Breaking Changes

- In the CpGsInfoAllRegions(), CpGsInfoOneRegion(), and
GetCpGsInRegion() functions, note that we have added a new argument
region_gr (in the second position) which allows for the input of a
GRanges object to these functions. This will cause errors in
existing code if that code uses positional argument matching. Any
existing code that matches arguments by name will be unaffected.

                       Changes in version 0.99.2                        

Major Changes

Bolstered documentation across all package help / manual files

Bug Fixes

                       Changes in version 0.99.1                        

Migrated all clean files to this repository. All development history
in
https://github.com/TransBioInfoLab/coMethDMR_old

                     Changes in version 0.0.0.9001                      

Major Changes

Included EPIC arrays annotation (#1) Added a check that: (#5) The
rownames of the beta matrix are probe IDs, and There are only numeric
columns in the beta matrix

Bug Fixes

[compcodeR](/packages/compcodeR)
---------

                       Changes in version 1.31.1                        

- Implement simulations according to the phylogenetic poisson log
  normal model

- Implement length-varying simulations

- New PhyloCompData object for phylogenetic informed simulation data

- Add possible length normalization to DESeq2 and limma analyses

- Add phylolm analyses

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

                       Changes in version 2.11.1                        

- add a global option `ht_opt$COLOR` to control colors for continuous
  color mapping.

- `annotation_label` can be an `expression` object.

- `recycle_gp()`: now consider when n = 0.

- `anno_block()`: add `align_to` argument.

- add `anno_text_box()` and `grid.text_box()`.

- add `show_name` argument in `anno_empty()`.

- the validation of annotations in `HeatmapAnnotation()` is simplified.

- add `anno_numeric()`.

- when `rect_gp = gpar(type = "none")`, `use_raster` is enforced to be
  FALSE.

- "global variables" outside `cell_fun`/`layer_fun` are aotumatially
  identified and saved locally.

[CompoundDb](/packages/CompoundDb)
----------

                        Changes in version 0.99                         

Changes in version 0.99.12

- Add mass2mz method for CompDb databases.

Changes in version 0.99.11

- Add peaksVariables method.

Changes in version 0.99.10

- Add parameter columns to peaksData.

Changes in version 0.99.9

- Add parameter dbFile to createCompDb and add an example on how to
create a CompDb database from custom input.

Changes in version 0.99.8

- Add citation.

Changes in version 0.99.7

- Add bug reports link to DESCRIPTION.

Changes in version 0.99.6

- MsBackendCompDb extends Spectra::MsBackendCached instead of
Spectra::MsBackendDataFrame.

Changes in version 0.99.5

- No updates, just version bump to cause a new build.

Changes in version 0.99.4

- Address more comments from @jianhong.

Changes in version 0.99.3

- Fix BiocCheck warnings.

Changes in version 0.99.2

- Fix BiocCheck warnings.

Changes in version 0.99.1

- Address comments/change requests from @jianhong.

Changes in version 0.99.0

- Preparing for Bioconductor submission.

                         Changes in version 0.9                         

Changes in version 0.9.4

- Add deleteIon and deleteSpectra allowing to delete ions or spectra.

Changes in version 0.9.3

- insertIons supports adding additional database columns.

Changes in version 0.9.2

- Add instertSpectra method to add MS/MS spectra from a Spectra
object
to the database.

Changes in version 0.9.1

- Add IonDb constructor methods.
- Expand documentation and examples.
- Add and fix unit tests.

Changes in version 0.9.0

- Add IonDb class as extension of CompDb (to allow adding ion
information to the database) and the functionalities to create such
object.
- Add insertIon to allow adding new ions to an IonDb object
- Add ionVariables, ions functions to access the ions data in the
database.
- Add filters: IonIdFilter, IonAdductFilter, IonMzFilter,
IonRtFilter.

                         Changes in version 0.8                         

Changes in version 0.8.1

- Import spectra type (MS level) and precursor type from MoNa.

Changes in version 0.8.0

- Rename database table name compound into ms_compound issue #74.

                         Changes in version 0.7                         

Changes in version 0.7.0

- Remove mass2mz and mz2mass function in favour of the functions
implemented in MetaboCoreUtils.

                         Changes in version 0.6                         

Changes in version 0.6.6

- Import compounds method from ProtGenerics.

Changes in version 0.6.5

- Add parameter onlyValid to compound_tbl_sdf to allow importing of
only valid elements issue #69.

Changes in version 0.6.4

- Add additional filters: MassFilter, FormulaFilter, InchiFilter and
InchikeyFilter.

Changes in version 0.6.3

- Add metadata, spectraVariables and compoundVariables functions.

Changes in version 0.6.2

- Support creation of databases without specifying the organism.
- Ensure database columns are mapped correctly to Spectra variable
names.

Changes in version 0.6.1

- Add SpectrumIdFilter to support filtering by spectrum IDs.

Changes in version 0.6.0

- Rename column names: compound_name -> name, mass -> exactmass,
inchi_key -> inchikey.

                         Changes in version 0.5                         

Changes in version 0.5.0

- Replace as.list with peaksData.
- Replace asDataFrame with spectraData.

                         Changes in version 0.4                         

Changes in version 0.4.3

- Updated to match new LIPID MAPS field names.

Changes in version 0.4.2

- Fix bug in as.list,MsBackendCompDb which returned a SimpleList
instead of a list.

Changes in version 0.4.0

- Rename method spectraData for MsBackendCompDb into asDataFrame
(adapting to the changes in Spectra).

                         Changes in version 0.3                         

Changes in version 0.3.2

- Import also smiles from SDF files.

Changes in version 0.3.1

- Move package Spectra from Depends to Imports

Changes in version 0.3.0

- Change from MSnbase to RforMassSpectrometry packages (Spectra and
MsCoreUtils).
- Store MS/MS spectra in two tables, msms_spectrum and
msms_spectrum_peak.

                         Changes in version 0.2                         

Changes in version 0.2.3

- Add instrument and precursor_mz spectra data columns (issue #32).

Changes in version 0.2.2

- Add adduct information from Jan Stanstrup's commonMZ package.
- Add matchWithPpm function to match numeric values allowing for a
small difference.
- Add adducts function to retrieve adduct definitions.
- Add mass2mz and mz2mass to convert between mass and m/z for
provided
adducts.
- Add annotateMz method to annotate m/z values.

Changes in version 0.2.1

- Change field collision_energy to character to support values from
MoNa (issue #31).
- Add functions import_mona_sdf and msms_spectra_mona functions to
enable import of spectrum data from MoNa SDF files (issue #30).
- Add support for MoNa SDF files (issue #30).

Changes in version 0.2.0

- Add hasMz,Spectrum and hasMz,Spectra methods to look for m/z values
within spectra (issue #28).
- Add MsmsMzRangeMinFilter and MsmsMzRangeMaxFilter (issue #29).
- Re-use Spectra object from MSnbase.
- Add supportedFilters,CompDb method.

                         Changes in version 0.1                         

Changes in version 0.1.1

- Add precursorMz, precursorCharge, precursorIntensity,
acquisitionNum, scanIndex, peaksCount, msLevel, tic, ionCount,
collisionEnergy, fromFile, polarity, smoothed, isEmpty, centroided
and isCentroided methods for Spectrum2List.

Changes in version 0.1.0

- Add expandMzIntensity function.
- Add spectra method to extract spectra from the CompDb database.
- Add functionality to store MS/MS spectra in a CompDb database (m/z
and intensity values stored as BLOB).
- Add functionality to load MS/MS spectra from HMDB xml files.

[conclus](/packages/conclus)
-------

                 Changes in version 1.3.4 (2022-03-28)                  

- Made the following significant changes
- Add parameter removeNoSymbol in the method normalizeCountMatrix
to filter genes doesn't have SYMBOL.
- clusterMarkers slot become topMarkers slot, accessors names
change in the same way.0

[CoreGx](/packages/CoreGx)
------

                        Changes in version 2.0.0                        

- The @cell slot has become the @sample slot. Associated generics and
accessor methods have been renamed, then aliased to their old names.
As such, old code should still work as expected, but will in fact be
calling different S4 methods.
- Added the @treatment slot to the CoreSet-class

[COTAN](/packages/COTAN)
-----

                       Changes in version 0.99.12                       

Release before the official release of Bioc 3.15. Main changes:

- The way in which the COEX matrix is estimated and stored is changed
to occupy less space and run faster.

                       Changes in version 0.99.10                       

Initial Bioconductor release

[crisprBase](/packages/crisprBase)
----------

                        Changes in version 1.0.0                        

- New package crisprBase, Base functions and classes for CRISPR
  gRNA design.

[crisprBowtie](/packages/crisprBowtie)
------------

                        Changes in version 1.0.0                        

- New package crisprBowtie, Bowtie-based alignment of CRISPR gRNA
  spacer sequences.

[crisprScore](/packages/crisprScore)
-----------

                        Changes in version 1.0.0                        

- New package crisprScore, for on-target and off-target scoring
  of guide RNAs (gRNAs).

[csdR](/packages/csdR)
----

                        Changes in version 1.1.3                        

Added citation to the csdR article which is now printed.

                        Changes in version 1.1.2                        

Made it explicit in the documentation that missing values are not
allows
and wrote a test for this case.

                        Changes in version 1.1.1                        

Made some minor modifications to the README and the vignette.

[customCMPdb](/packages/customCMPdb)
-----------

                 Changes in version 1.5.1 (2022-01-31)                  

- Add DrugAge build 4 database

[cytomapper](/packages/cytomapper)
----------

                 Changes in version 1.7.2 (2022-04-20)                  

- Fixes for next release

                 Changes in version 1.7.1 (2021-12-28)                  

- Added compImage function for channel spillover compensation

[cytoMEM](/packages/cytoMEM)
-------

                 Changes in version 0.99.2 (2022-04-11)                 

- Updated show methods for build_heatmaps function
- Updated assignment from = to <-

                 Changes in version 0.99.1 (2022-03-21)                 

- Updated class membership checks with is() and == / != to is(x,
'class') for S4 classes
- Version bump for Bioconductor

                 Changes in version 0.99.0 (2022-03-18)                 

- Submitted to Bioconductor

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

[dasper](/packages/dasper)
------

                        Changes in version 1.5.1                        

NEW FEATURES

- Fix bug related to coverage_norm(), ensuring that regions used to
normalise coverage have the same seqlevels as the inputted
junctions.

[dce](/packages/dce)
---

                 Changes in version 1.3.5 (2022-04-14)                  

- Update package data df_pathway_statistics

[dearseq](/packages/dearseq)
-------

                 Changes in version 1.7.2 (2022-04-13)                  

- reintroduced CompQuadForm::davies() as a longer alternative when
the
"saddlepoint" method in survey::pchisqsum() fails for computing
quadratic form asymptotic p-values

                 Changes in version 1.7.1 (2022-02-07)                  

- added a spaghetti plot functionality

[decoupleR](/packages/decoupleR)
---------

                         Changes in version 2.2                         

Changes

- Changed example mat and net to toy examples.

- Changed test data to toy data.

                         Changes in version 2.1                         

Changes

- likelihood param is deprecated, from now on, weights (positive or
negative) should go to the mor column of network. Methods will still
run if likelihood is specified, however they will be set to 1.

- Added minsize argument to all methods, set to 5 by default. Sources
containing less than this value of targets in the input mat will be
removed from the calculations.

- Changed default behavior of the decouple function. Now if no
methods
are specified in the statistics argument, the function will only run
the top performers in our benchmark (mlm, ulm and wsum). To run all
methods like before, set statistics to 'all'. Moreover, the argument
consensus_stats has been added to filter statistics for the
calculation of the consensus score. By default it only uses mlm, ulm
and norm_wsum, or if statistics=='all' all methods returned after
running decouple.

- viper method:

- Now properly handles weights in mor by normalizing them to -1
and +1.

- ulm/mlm/udt/mdt methods:

- Changed how they processed the input network. Before the model
matrix only contained the intersection of features between mat
and network's targets, now it incorporates all features coming
from mat ensuring a more robust prediction. Prediction values
may change slightly from older versions.
- Deprecated sparse argument.

- ora method:

- Now takes top 5% features as default input instead of 300 up and
bottom features.
- Added seed to randomly break ties

- consensus method:

- No longer based on RobustRankAggreg. Now the consensus score is
the mean of the activities obtained after a double tailed
z-score transformation.

- Discarded filter_regulons function.

- Moved major dependencies to Suggest to reduce the number of
dependencies needed.

- Updated README by adding:

- Kinase inference example
- Graphical abstract
- Manuscript and citation
- New vignette style

- Updated documentation for all methods.

New features

- Added wrappers to easily query Omnipath, one of the largest
data-bases collecting prior-knowledge resources. Added these
functions:

- show_resources: shows available resources inside Omnipath.
- get_resource: gets any resource from Omnipath.
- get_dorothea: gets the DoRothEA gene regulatory network for
transcription factor (TF) activity estimation. Note: this
version is slightly different from the one in the package
dorothea since it contains new edges and TFs and also weights
the interactions by confidence levels.
- get_progeny: gets the PROGENy model for pathway activity
estimation.

- Added show_methods function, it shows how many statistics are
currently available.

- Added check_corr function, it shows how correlated regulators in a
network are. It can be used to check for co-linearity for mlm and
mdt.

- Added new error for mlm when co-variables are co-linear (regulators
are too correlated to fit a model).

Bugfixes

- wmean and wsum now return the correct empirical p-values.

- ulm, mlm, mdt and udt now accept matrices with one column as input.

- Results from ulm and mlm now correctly return un-grouped.

- Methods correctly run when mat has no column names.

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

[DEGreport](/packages/DEGreport)
---------

                       Changes in version 1.31.6                        

- Migrate gather to pitvot_longer

- Fix warning in melt in degcovariates

- Fix significants function, issue with !!!sym

                       Changes in version 1.31.3                        

- Fix error in example when using `head` after pipe

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.22.0                        

DEPRECATED AND DEFUNCT

- The following stuff is now defunct after being deprecated in previous
  versions of the package:
  - blockGrid(): replaced with defaultAutoGrid()
  - rowGrid(): replaced with rowAutoGrid()
  - colGrid(): replaced with colAutoGrid()
  - multGrids(): replaced with defaultMultAutoGrids()
  - linearInd(): replaced with Mindex2Lindex()
  - viewportApply(): replaced with gridApply()
  - viewportReduce(): replaced with gridReduce()
  - getRealizationBackend(): replaced with getAutoRealizationBackend()
  - setRealizationBackend(): replaced with setAutoRealizationBackend()
  - RealizationSink(): replaced with AutoRealizationSink()

BUG FIXES

- Small tweak to updateObject() method for DelayedArray objects (see
  commit abcd154).

[DepecheR](/packages/DepecheR)
--------

                       Changes in version 1.11.4                        

- sPLSDA function testing excluded from BioConductor, as it works
  locally, but
  fails at BioConductor.

                 Changes in version 1.11.3 (2022-03-12)                 

- Bug fix in dColorPlot affecting character vectors and factors used as
  input.

                 Changes in version 1.11.2 (2022-03-11)                 

- Searching for unidentified bug in test for dSlsda function.

                 Changes in version 1.11.1 (2022-03-05)                 

- Bug fix in dColorPlot, solving problems with individual continuous
  vectors

- microClust, an internal function, has been updated to include the
  counting of
  donors in addition to median and mean calculations among the
  neighbors.

- The addition of nUniqueNeihgDons, in conjunction to the
  aforementioned change
  to the microClust function.

- Minor updates to dSplsda function testing, to hopefully remove the
  problem of
  intermittent problems in the testing on the Linux plaform.

[DepInfeR](/packages/DepInfeR)
--------

                 Changes in version 0.99.7 (2022-04-05)                 

- Modify package based on the comments from the second round revision

                 Changes in version 0.99.4 (2022-03-19)                 

- Modify package according to comments from Bioconductor team

                 Changes in version 0.99.0 (2022-01-13)                 

- Submitted to Bioconductor

[DEWSeq](/packages/DEWSeq)
------

                 Changes in version 1.9.2 (2022-04-19)                  

- add function to filter read counts based on max. position count

[DExMA](/packages/DExMA)
-----

                        Changes in version 1.3.1                        

- missGenesImput function added

- SampleKnn imputation method for unmeasured genes added

[DiffBind](/packages/DiffBind)
--------

                         Changes in version 3.6                         

- Change default spike-in normalization to NATIVE (TLE/TMM)

- Change spike-in normalization to use reference library sizes

- Fix bug involving Called columns in reports

- Maintain FRiP calculations

- Improve dba.plotProfile sample labels

- Fix issue with package:parallel

- Fix error/warning in dba.blacklist relating to non-matching
  chromosome names

- Update man page for dba.report to clarify how fold changes are
  reported when bNormalized=FALSE.

- Add some new test conditions to GenerateDataFiles.R

[distinct](/packages/distinct)
--------

                        Changes in version 1.7.2                        

- differential testing allowed, even when 1 gene only is provided

                        Changes in version 1.7.1                        

- when removing nuisance covariates, the previous linear model was
  replaced with a linear mixed effect model with: i) fixed effects for
  the nuisance covariates, and ii) random effects for the samples. This
  properly accounts for the correlation structure of cells belonging to
  the same sample.

[dittoSeq](/packages/dittoSeq)
--------

                         Changes in version 1.8                         

- Minor Feature Add: 'randomize' option for 'order' input of
'dittoDimPlot()' and 'dittoScatterPlot()'

[DMRcate](/packages/DMRcate)
-------

                        Changes in version 2.8.1                        

- CITATION updated.

[DOSE](/packages/DOSE)
----

                       Changes in version 3.21.2                        

- enable setReadable for compareCluster(GSEA algorithm)
result(2021-12-13, Mon)
- update the default order of GSEA result (2021-12-09, Thu)
- if p.adjust is identical, sorted by abs(NES)

                       Changes in version 3.21.1                        

- upate DisGeNET and NCG data (2021-11-14, Sun)
- DisGeNET v7: 21671 genes, 30170 diseases and 1134942
gene-disease associations
- 194515 variants, 14155 diseases and 369554 variant-disease
associations
- NCG v7: 3177 cancer genes, 130 diseases and 6095 gene-disease
associations

[edgeR](/packages/edgeR)
-----

                       Changes in version 3.38.0                        

- 
  New argument 'keep.EList' for voomLmFit() to store the
  normalized log2-CPM values and voom weights.

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.15.4                        

- update treeplot: support passing rel object to offset and
offset_tiplab (2022-04-24, Sun)

                       Changes in version 1.15.3                        

- export `drag_network' (2022-03-07, Mon)
- update cnetplot.enrichResult to be supported by
drag_network(2022-3-6, Sun)
- add function drag_network to drag the nodes of networks (2022-2-25,
Fri)
- fix a bug in goplot: goplot.gseaResult need setType slot instead of
ontology slot (2022-2-22, Tue)
- return gg object instead of print it in
dotplot.compareClusterResult() (2022-01-05, Wed, @altairwei, #160)

                       Changes in version 1.15.2                        

- add label_format_tiplab and label_format_cladelab parameters for
treeplot(2021-12-24, Fri)
- support treeplot of compareCluster(GSEA algorithm)
result(2021-12-13, Mon)
- support visualization of compareCluster(GSEA algorithm)
result(2021-12-11, Sat)
- support scientific notation for gseaplot2(2021-12-4, Sat)

                       Changes in version 1.15.1                        

- fixed R check by importing utils

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.19.10                       

- Fix issue in `getGeneRegionTrackForGviz` if no transcript was found
  in the
  provided genomic range (issue
  https://github.com/jorainer/ensembldb/issues/132).

                       Changes in version 2.19.9                        

- `seqlevelsStyle` allows to provide a custom mapping data frame.

                       Changes in version 2.19.8                        

- Require package Biostrings in the coordinate mapping vignette.

                       Changes in version 2.19.7                        

- `listColumns` does no longer report columns from the metadata
  database table
  (issue https://github.com/jorainer/ensembldb/issues/128).

                       Changes in version 2.19.6                        

- Fixes in `proteinToGenome`.

                       Changes in version 2.19.5                        

- Add support for additional `tx_is_canonical` column and add
  `TxIsCanonicalFilter` (issue
  https://github.com/jorainer/ensembldb/issues/123).

                       Changes in version 2.19.4                        

- Fix issue with quotes in transcript names when importing transcript
  tables.

                       Changes in version 2.19.2                        

- Restore backward compatibility (issue
  https://github.com/jorainer/ensembldb/issues/122).

                       Changes in version 2.19.1                        

- Add database column tx_name to store the external transcript names.

- Update TxNameFilter to allow filtering on this database column.

[ensemblVEP](/packages/ensemblVEP)
----------

                       Changes in version 1.37.0                        

- add support for Ensembl release 105

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.3.6 (2022-02-16)                  

- significant speed-up (1.3.5)

- methylation patterns

                 Changes in version 1.3.2 (2021-12-24)                  

- more efficient data handling (XPtr instead of Rcpp::wrap'ping)

[EpiCompare](/packages/EpiCompare)
---------

                     Changes in version 0.99.3

New Features

- New functions with examples/unit tests:
  - `import_narrowPeak`: Import narrowPeak files, with automated header annotation using metadata from ENCODE.\
  - `gather_files`: Automatically peak/picard/bed files and read them in as a list of `GRanges` objects.\
  - `write_example_peaks`: Write example peak data to disk.
- Update *.gitignore*
- Update *.Rbuildignore*

                     Changes in version 0.99.1

New features

- New parameter in EpiCompare:
  - `genome_build`: Specify the genome build, either "hg19" or "hg38". This parameter is also included in `plot_chromHMM`, `plot_ChIPseeker_annotation`, `tss_plot` and `plot_enrichment`.

                     Changes in version 0.99.0

New Features

- `EpiCompare` submitted to Bioconductor.

[epigraHMM](/packages/epigraHMM)
---------

                        Changes in version 1.3.2                        

- Exporting function estimateTransitionProb to estimate transition
probabilities from a sequence of states of a Markov chain

                        Changes in version 1.2.1                        

- Bug fix of output paths to handle paths with '.'

[epimutacions](/packages/epimutacions)
------------

                       Changes in version 0.99.33                       

- Bugs and notes (if possible) fixed

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[escape](/packages/escape)
------

                        Changes in version 1.4.1                        

- Version number and small edits for bioconductor compliance

- Removed singscore method

- Added UCell functions internally so they are compatible with
  Bioconductor

[EWCE](/packages/EWCE)
----

                        Changes in version 1.3.3                        

New features

- method argument from orthogene::create_background and
orthogene::convert_orthologs is now passed up as an argument to EWCE
functions to give users more control. "homologene" chosen as default
for all functions. "homologene" has fewer species than "orthogene"
but doesnt need to import data from the web. It also has more 1:1
mouse:human orthologs.
- Include notes on mismatches between GitHub documentation and
current
Bioc release version.
- Allow bin_specificity_into_quantiles to set specificity matrix name
produced.
- Merge GHA workflow yamls into one.

Bug fixes

- Add try({}) and error=TRUE to avoid "polygon edge not found" error
in vignettes.

                        Changes in version 1.3.1                        

New features

- Major changes: Pull Request from bschilder_dev branch.
- All functions can now use lists and CellTypeDatasets (CTD) from any
species and convert them to a common species (human by default) via
orthogene.
- Automated CTD standardisation via standardise_ctd.
- Can handle (sparse) matrices.
- Can create CTD from very large datasets using DelayedArray object
class.
- All functions automatically create appropriate gene backgrounds
given species.
- More modular, simplified vignettes.
- Additional gene pre-filtering options (DESeq2, MAST, variance
quantiles).
- New/improved plotting functions (e.g. plot_ctd).
- Added example bootstrapping enrichment results as extdata to speed
up examples (documented in data.R). Accessed via
EWCE::example_bootstrap_results().
- Replaced GHA workflow with check-bioc to automatically: run R-CMD
checks, run BiocCheck, and rebuild/deploy pkgdown site.
- Parallelised functions:
- drop_uninformative_genes
- generate_celltype-data
- bootstrap_enrichment_test
- Added tests (multiple functions tests per file to reduce number of
times ewceData files have to be downloaded):
- test-DelayedArray
- test-merge_sce
- test-get_celltype_table
- test-list_species
- test-run_DGE
- test-check_percent_hits
- Added function is_32bit() to all tests to ensure they don't get run
twice on Windows OS.
- Added GitHub Actions workflows:
- check-bioc-docker.yml: Runs CRAN/Bioc checks, rebuilds and
pushes pkgdown website, runs and uploads test coverage report,
- dockerhub.yml: Builds Bioconductor Docker container with EWCE
installed, runs CRAN checks and (if checks are successful)
pushes container to neurogenomicslab DockerHub.
- Removed docs folder, as the documentation website comes from the
gh-pages branch now, and is automatically built by GHA workflow
after each push to main branch.
- Added new exported function fix_celltype_names to help with
standardising celltype names in alignment with standardise_ctd.
- generate_bootstrap_plots_for_transcriptome: Now supports any
species
(not just mouse or human).
- Converts CTD and DGE table (tt) into output_species gene
symbols.
- Automatically generates appropriate gene background.
- Faster due to now having the option to only generate certain
plot types.
- Provide precomputed results from ewce_expression_data via new
example_transcriptome_results function.
- Reduced build runtime and oversized vignettes by not evaluating
certain code chunks.
- Prevent extended vignette from running entirely.
- @return documentation for internal functions.
- Added more installation checks to GHA.
- Fixed inconsistent naming of unit test files: test_ ==> test-
- Removed DGE args in drop_uninformative_genes for now until we run
benchmarking to see how each affects the EWCE results.
- Make bootstrap_plots function internal.
- Add report on how orthogene improve within- and across-species gene
mappings in extended vignette.
- Record extra info in standardise_ctd output:
- "species": both input_species and output_species
- "versions": of EWCE, orthogene, and homologene

[exomePeak2](/packages/exomePeak2)
----------

                 Changes in version 1.7.0 (2022-04-20)                  

- Improved accuracy by introducing new computational models of GC
  content bias correction and IP efficiency correction.

- Improved peak calling stability by applying Poisson test under
  default setting.

- Read count method is improved.

- Multi-step functions and quantification mode are deprecated.

- Unnecessary function parameters are removed.

[ExperimentHub](/packages/ExperimentHub)
-------------

                        Changes in version 2.3.0                        

USER VISIBLE MODIFICATIONS

- (2.3.5) Update documentaion for shared hub across multiple users

[ExperimentHubData](/packages/ExperimentHubData)
-----------------

                       Changes in version 1.21.0                        

MODIFICATIONS

- 1.21.2 Changed makeExperimentHubMetadata to point to Azure instead of
  AWS

[fastreeR](/packages/fastreeR)
--------

                        Changes in version 1.0.0                        

- New package fastreeR, Phylogenetic, Distance and Other
  Calculations on VCF and Fasta Files.

                 Changes in version 0.99.7 (2022-04-02)                 

- Updates and corrections after the 1st Bioconductor review.

                 Changes in version 0.99.6 (2022-03-26)                 

- Drop function dist2hist.

                 Changes in version 0.99.5 (2022-03-25)                 

- Update vignette's use of dist2hist.

                 Changes in version 0.99.4 (2022-03-25)                 

- Update java backend (createHistogram).

                 Changes in version 0.99.3 (2022-03-25)                 

- Update java backend.

                 Changes in version 0.99.2 (2022-03-25)                 

- Update README.md to inform about JDK>=8 requirement.

- Update java dependencies (jfreechart-1.5.3).

                 Changes in version 0.99.1 (2022-03-24)                 

- Update vignette to use BiocFileCache so that sample vcf and
  fasta downloads not get repeated needlessly.

                 Changes in version 0.99.0 (2022-03-21)                 

- Submitted to Bioconductor.

[fgsea](/packages/fgsea)
-----

                       Changes in version 1.21.1                        

- fix a reproducibility problem
  (https://github.com/ctlab/fgsea/issues/110)

[FindIT2](/packages/FindIT2)
-------

                 Changes in version 1.0.3 (2021-12-29)                  

- fix bug when chrosome length is too small (calcRP_coverage
function)

                 Changes in version 1.0.2 (2021-11-23)                  

- fix bug when type in colData is factor (integtate_replicates
function)

                 Changes in version 1.0.1 (2021-11-09)                  

- simplify loadPeakFile function

[fishpond](/packages/fishpond)
--------

                        Changes in version 2.2.0                        

- New vignette demonstrating allelic analysis at isoform, TSS, or
gene-level. See more details below.
- New import functions for equivalence class analysis of Salmon or
alevin data, written by Jeroen Gilis. See salmonEC() and alevinEC()
man pages.
- New plotAllelicGene() and plotAllelicHeatmap() functions for
plotting results from allelic expression analysis.
- New makeTx2Tss() helper function for allelic analysis.
- Now importFromAllelicCounts() can take a GRanges object as the
tx2gene argument, so that ranges will be distributed to the
rowRanges of the output SummarizedExperiment.
- Adding shiftX argument to plotInfReps() for numeric x variable, to
help with overplotting.
- Re-organized package for new pkgdown homepage:
https://mikelove.github.io/fishpond

[flowAI](/packages/flowAI)
------

                       Changes in version 1.25.3                        

- Added option to use IQR to make the flow rate check less sensitive.

[flowSpecs](/packages/flowSpecs)
---------

                 Changes in version 1.9.2 (2022-02-15)                  

- Adding a hexBin parameter to the oneVsAllPlot function

- Adding the dependency hexbin.

                 Changes in version 1.9.1 (2021-12-18)                  

- Introduction of a dirName parameter to the oneVsAllPlot function.

- Increased versatility for the madFilter function.

[fmrs](/packages/fmrs)
----

                        Changes in version 2.0.0                        

IMPROVEMENTS SINCE LAST RELEASE

- The package is rewritten using .Call function.
- The codes for Weibull distribution are improved.

BUG FIXES

- Several bugs are fixed which caused the results to be different for
the same analysis.

[fobitools](/packages/fobitools)
---------

                        Changes in version 1.3.2                        

- Fix bugs in vignettes

[FRASER](/packages/FRASER)
------

                        Changes in version 1.6.1                        

- Fixing quantile filtering defaults (#28)

- Require min expression in 5% instead of 95% of the samples

- Require min expression on both sides of the junction

- Align FRASER package with DROP pipeline (#24)

- Move temp directory from tempdir() to working directory getwd()

- Improve visualizations and Improve documentation

- Improve internal object validation

- Minor bugfixes

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.32.0                        

UTILITIES

- optimize using the utilities of the Matrix package for sparse
  matrices

[GeneNetworkBuilder](/packages/GeneNetworkBuilder)
------------------

                       Changes in version 1.37.2                        

- Add layout button for htmlwidgets

                       Changes in version 1.37.1                        

- Fix the version issue of STRINGdb

[GENESIS](/packages/GENESIS)
-------

                       Changes in version 2.25.5                        

- Add new test option "fastSMMAT" to assocTestAggregate.

                       Changes in version 2.25.4                        

- Fixed a bug in pcrelate when running with both multiple sample
  blocks and multiple cores

                       Changes in version 2.25.3                        

- Fixed a bug in assocTestSingle for computing saddle-point
  approximation (SPA) p-values for variants (i.e. when using test
  = Score.SPA) when the input null model is not a mixed models
  (i.e. when fitNullModel was run with cov.mat = NULL , or all
  variance components converged to 0) and the family = "binomial"

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 2.0.0                        

New features

- GeneTonic now offers the possibility to upload a GeneTonicList at
runtime. This makes it possible to use the app as a server-like
dashboard, which runs by default on no dataset provided, and
populates its components upon successfully providing the data as
expected
- The GeneSpector functionality in the Welcome panel provides a means
to explore any gene in the expression set, coloring and grouping by
any experimental covariate of interest
- It is possible to enter a set of genes and genesets in the
Bookmarks
panel, and these can be doubled checked against the available
features of the current GeneTonicList - this, combined to the upload
functionality, makes it possible to easily compare different gtl
objects
- The GeneTonic app has a button to export the currently provided
dataset - regardless of the input format - as a GeneTonicList. This
is especially useful if one is providing the individual components
(dds, res_de, res_enrich, annotation_obj) and would like to obtain
the correct serialized object
- gs_upset adds the possibility to represent the results of
enrichment
analyses as upset plots, with the option to decorate them with
DE-related information

Other notes

- The manuscript about GeneTonic is now published on BMC
Bioinformatics at
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04461-5
- the citation item has been updated accordingly
- The jittered position in the gene_plot is now completely
reproducible, by setting a seed internally in the jitter generating
function

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.32.0                        

DEPRECATED AND DEFUNCT

- releaseName() is now defunct on GenomeDescription objects.

- Remove fetchExtendedChromInfoFromUCSC() and available.species().
  Both were defunct in BioC 3.14.

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.32.0                        

- No significant changes in this version.

[GenomicDataCommons](/packages/GenomicDataCommons)
------------------

                       Changes in version 1.20.0                        

New features

- gdcdata has an ellipses argument to download data from the legacy
archive, e.g., legacy = TRUE (#84, @LiNk-NY)
- missing (is MISSING) and !missing (NOT MISSING) operations
implemented for filtering queries, see vignette (#96, @LiNk-NY)
- gdc-client version can be validated against last known good version
based on data release (#99, @LiNk-NY)

Bug fixes and minor improvements

- gdc_clinical uses readr::type_convert to handle columns with
inconsistent types from the API.
- update examples in documentation and vignette based on new data
release

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.48.0                        

NEW FEATURES

- Add proteinToGenome() generic and methods. Loosely modeled on
  ensembldb::proteinToGenome().

- Add extendExonsIntoIntrons().

SIGNIFICANT USER-VISIBLE CHANGES

- Use useEnsembl() instead of useMart() wherever possible.

- makeTxDbFromGFF() and makeTxDbFromGRanges() now recognize and
  import features of type "gene_segment", "pseudogenic_gene_segment",
  "antisense_RNA", "three_prime_overlapping_ncrna", or "vault_RNA",
  as transcripts.
  Note that, according to the Sequence Ontology, the "gene_segment"
  and "pseudogenic_gene_segment" terms are NOT offsprings of the
  "transcript" term via the **is_a** relationship but we still treat
  them as if they were because that's what Ensembl does in their GFF3
  files!

DEPRECATED AND DEFUNCT

- Add warning about FeatureDb-related functionalities no longer being
  actively maintained.

- disjointExons() is now defunct after being deprecated in BioC 3.13.

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.48.0                        

NEW FEATURES

- Add subtract() for subtracting a set of genomic ranges from a GRanges
  object. This is similar to bedtools subtract.

- Add 'na.rm' argument to makeGRangesFromDataFrame().

DEPRECATED AND DEFUNCT

- Remove the GenomicRangesList() constructor. This constructor got
  deprecated in BioC 3.10 and defunct in BioC 3.13.

BUG FIXES

- Make sure promoters() works on GPos objects.

[GenomicScores](/packages/GenomicScores)
-------------

                        Changes in version 2.8.0                        

USER VISIBLE CHANGES

- Added support to the latest version 3.1.2 of gnomAD MAF data, stored
  in the package MafH5.gnomAD.v1.1.2.GRCh38.

BUG FIXES

- Bugfix on a problem caused by unordered input genomic ranges of width
  larger than 1, whose result was return in the wrong order; see
  https://github.com/rcastelo/GenomicScores/issues/18

- Bugfix when accessing multiple populations of scores stored on an
  HDF5 backend.

[GenomicSuperSignature](/packages/GenomicSuperSignature)
---------------------

                        Changes in version 1.4.0                        

- [Bug Fix]
  - `findStudiesInCluster`: single-element clusters return the correct
  study now, instead of `null`

- [Major]
  - `findStudiesInCluster`: the output includes PC number of the
  participating
  study in the cluster and the variance explained by them. With
  `studyTitle=FALSE`,
  the output will be a data frame, not a character vector.
  - `subsetEnrichedPathways`: the new argument `include_nes` is added.
  If it
  is set to `TRUE`, the output will include NES from GSEA.
  - `getRAVInfo` and `getStudyInfo`: two new functions to extract basic
  metadata for RAVs and studies, respectively.
  - We characterize RAVs that are harder to interpret with the
  currently available
  information associated with them (more detail can be found
  bit.ly/RAVmodel_characterization).
  For the following functions, any output including those RAVs will
  have a
  defulat message: `meshTable`, `drawWordcloud`, `heatmapTable`,
  `validatedSignatures`.
  You can snooze this by setting `filterMessage = FALSE`
  - New required parameter `RAVmodel` for the following functions:
  `annotatePC`,`heatmapTable` and `validatedSignatures`.
  - New accessor `version` is available to check the version of
  RAVmodel
  - New function `availableRAVmodel` will output the different versions
  of
  RAVmodels available for downloading now.
  - `getModel` function now takes a new argument, `version`, to specify
  the
  version of RAVmodel to download.

- [Minor]
  - miniRAVmodel is updated

[GenProSeq](/packages/GenProSeq)
---------

                 Changes in version 0.99.0 (2021-09-20)                 

- submission to Bioconductor

[GeoDiff](/packages/GeoDiff)
-------

                        Changes in version 1.1.2                        

- transfer maintainership

                        Changes in version 1.1.1                        

- change default values in fitNBth and NBthmDE wrapper functions to
be
consistent

                        Changes in version 1.1.0                        

- accepted by Bioconductor

[GeomxTools](/packages/GeomxTools)
----------

                       Changes in version 2.99.2                        

- documentation updates

                       Changes in version 2.99.1                        

- documentation updates

                        Changes in version 2.1.9                        

- documentation updates

                        Changes in version 2.1.8                        

- subset objects in tests & examples to decrease time

                        Changes in version 2.1.7                        

- coercions to Seurat and SpatialExperiment are S3 coercions

                        Changes in version 2.1.6                        

- added updateObject capability

                        Changes in version 2.1.5                        

New features:

- Handle multiple PKC file versions for a single module
- Only common probes in all versions for each module will be used
- Default behavior most recent PKC version for resolving probe
assignments
- User can override default with a specified version file name
- Added SystematicName and GeneID from PKC to feature metadata

                        Changes in version 2.1.4                        

- Add code for grubbs test from deprecated outliers package

                        Changes in version 2.1.3                        

New features:

- Allow protein NGS experiment data reading
- New slot, analyte, added to refer to analyte type

                        Changes in version 2.1.2                        

Bug fixes:

- Bug fix in outlier testing

                        Changes in version 2.1.1                        

New features:

- Add coercion to Seurat and SpatialExperiment

                        Changes in version 2.1.0                        

- No changes from 1.1.4

[GEOquery](/packages/GEOquery)
--------

                       Changes in version 2.63.2                        

- Added a NEWS.md file to track changes to the package.

[ggmanh](/packages/ggmanh)
------

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[ggmsa](/packages/ggmsa)
-----

                        Changes in version 1.1.4                        

updated the way smooth is invoked on simplot(2022-01-03, Mon)

added smoothed curve on simplot.(2021-12-17, Fri)

                        Changes in version 1.1.3                        

fixed the typo in "posHighligthed", and changed it to snake_case
"position_highlight" from camelCase "posHighligthed" (2021-12-13,
Mon)

                        Changes in version 1.1.2                        

fixed the assignment error on line 155 'seqlogo.R'

                        Changes in version 1.1.1                        

fixed error: using || instead of | on 110 lines in geom_msa.R

[ggtree](/packages/ggtree)
------

                        Changes in version 3.3.3                        

- geom_striplab() that supports aes() mapping (2022-04-22, Fri, #493)
- to.bottom parameter introduced in geom_hilight() to allow the
highlight layer was added into the lowest layer stack (2022-04-22,
Fri, #492)

                        Changes in version 3.3.2                        

- mv identify() method to 'ggfun' (2022-04-01, Fri)
- update identify.gg() to support 'ggplot' object and +xlim()
- update man files (2022-03-23, Wed, #489)

                        Changes in version 3.3.1                        

- use graph layouts to visualize tree (2021-12-10, Fri, #460, #461)
- igraph layout
- graphlayouts:
https://cran.r-project.org/web/packages/graphlayouts/index.html
- scale_color_subtree and scale_colour_subtree to color subtree via
taxa group information (e.g., cutree, or kmeans) (2021-12-01, Wed)
- set na.value = 'white' in msaplot() (2021-10-29, Fri)

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                        Changes in version 1.5.4                        

- set orientation argument automatically, when it is not be provided
and the geom has the parameter. (2022-03-24, Thu)

                        Changes in version 1.5.3                        

- update the citation format. (2022-01-28, Fri)

                        Changes in version 1.5.2                        

- import geom_text , which was used in the axis of geom_fruit with
do.call. (2021-11-24, Wed)

                        Changes in version 1.5.1                        

- update the lower limit of reference range in normalization precess.
(2021-11-19, Fri)

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.21.1                        

- Avoid eval-parse in load_OrgDb() (2022-01-10, Mon)

[graphite](/packages/graphite)
--------

                 Changes in version 1.41.2 (2022-04-21)                 

- Updated all pathway data.

[GSVA](/packages/GSVA)
----

                        Changes in version 1.50                         

BUG FIXES

- Bugfix for https://github.com/rcastelo/GSVA/issues/54 to force
  filtering genes with constant expression behaving the same regardless
  of the delayed or non-delayed nature of the data container.

[Harman](/packages/Harman)
------

                       Changes in version 1.23.1                        

- Code syntax tweaks to remove warnings raised in later version
  of R.

- Increase the forceRand unit test from 10,000 to 30,000
  iterations.  This solves an occasional unit test failure due to
  sampling.

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.24.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Improve error reporting in internal helper .h5openlocalfile()

BUG FIXES

- Make sure updateObject() handles very old HDF5ArraySeed instances.

[hermes](/packages/hermes)
------

                       Changes in version 0.99.4                        

- Corrected authors.

                       Changes in version 0.99.2                        

- First public release of the hermes package.
- Submission to BioConductor.

Miscellaneous

- New utility function cut_quantile for cutting a numeric vector into
quantiles.
- New utility function cat_with_newline for concatenating and
printing
with newline.
- New check function check_proportion which checks for a single
proportion.
- Better legends on the genes barplot and the correlation heatmap.
- Improved vignette layout using the BioConductor style.

                        Changes in version 0.1.1                        

New Features

- New function draw_scatterplot to produce scatterplots of two genes
or gene signatures.
- New function draw_boxplot for boxplots of gene expression values.
- New function draw_barplot for barplots of dichotomized gene
expression counts into two or three percentile categories.
- New helper function wrap_in_mae that wraps a single
SummarizedExperiment object into an MAE object.
- New method rename that makes renaming columns of rowData and
colData
as well as assay names in existing SummarizedExperiment objects much
easier, as a step before converting to HermesData.
- New method lapply that allows user to apply a function on all
experiments in a MultiAssayExperiment.
- New method isEmpty that checks whether a SummarizedExperiment
object
is empty.
- New gene filtering option n_top in the calc_pca function, which
allows filtering genes with greatest variability across samples.
- New class GeneSpec for specification of genes or gene signatures,
see ?gene_spec for simple construction. Inclusion of gene signature
functions colPrinComp1 and colMeanZscores to supplement standard
column statistics functions.
- New helper function col_data_with_genes which extracts the sample
variables saved in colData together with selected gene information
as a combined data set.
- New helper function inner_join_cdisc which joins genetic with CDISC
data sets.

Bug Fixes

- normalize() now also works when the hermes package is not loaded,
i.e. you can use it with hermes::normalize().
- correlate() now also works when there are factor variables in the
sample variables of the HermesData object.
- add_quality_flags() does no longer return NA as the technical
failure flags for the samples if there is only a single gene
contained in the input, but instead a vector of FALSE to ensure
correct downstream functionality.

Miscellaneous

- Updated LICENCE and README with new package references.
- The multi_assay_experiment now contains HermesData experiments,
different patient IDs, one experiment with normalized assays, and
multiple samples per patient in one experiment.
- The main HermesData example is now saved in the package as
hermes_data, and the previous summarized_experiment is still
available. Note that patient IDs have been changed in the new
version to align with the multi_assay_experiment.
- Renaming of required rowData and colData columns to be more
consistent with standards and use lowercase snake-case names.
- Annotation querying and setting is now more flexible in that it
also
allows to query more annotations than the required ones.
- Instead of gene starts and ends, the total length of gene exons is
now used as the annotation column size. Corresponding queries from
BioMart are used to return this gene size.
- df_char_to_factor has been deprecated (and can still be used with a
warning) and replaced with df_cols_to_factor, which also converts
logical variables to factor variables.
- When providing SummarizedExperiment objects containing
DelayedMatrix
assays to the HermesData() constructor, these are silently converted
to matrix assays to ensure downstream functionality.

                        Changes in version 0.1.0                        

- First internal release of the hermes package, which contains
classes, methods and functions to import, quality-check, filter,
normalize, and analyze RNAseq counts data for differential
expression.
- hermes is a successor of the rnaseqTools R package. The core
functionality is built on the BioConductor ecosystem, especially the
SummarizedExperiment class. New users should first begin by reading
the "Introduction to hermes" vignette to become familiar with the
hermes concepts.

New Features

- Import RNAseq count data into the hermes ready format.
- Annotate gene information from the Ensembl database via biomaRt.
- Add quality control (QC) flags to genes and samples.
- Filter and subset the data set.
- Normalize the counts.
- Produce descriptive plots.
- Perform principal components analysis.
- Produce a templated QC Rmd report.
- Perform differential expression analysis.

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.32.0                        

- fix the issue on Win32 because of using deprecated
  tbb::task_scheduler_init

                       Changes in version 1.30.2                        

- require GCC >= v8.0 for compiling the AVX-512VPOPCNTDQ intrinsics

- fix `hlaGDS2Geno()` when loading a SeqArray GDS file

[HiCDCPlus](/packages/HiCDCPlus)
---------

                 Changes in version 1.3.1 (2022-01-23)                  

- Fixed bug that prevents using different window sizes on TopDom

[hpar](/packages/hpar)
----

                        Changes in version 1.37                         

Changes in version 1.37.1

- Update to HPA release 21.0 <2021-11-22 Mon>

[HPAStainR](/packages/HPAStainR)
---------

                 Changes in version 1.4.1 (2021-11-23)                  

- Quick fix for inclusion of

[HPiP](/packages/HPiP)
----

                 Changes in version 1.1.1 (2021-11-29)                  

- Updated the vignette examples
- Added clustering detection algorithms

[HubPub](/packages/HubPub)
------

                 Changes in version 1.3.0 (2021-11-24)                  

- Update CreateAHubPackageVignette to use Azure instead of AWS

- Update AzureStor for auth_header

[idr2d](/packages/idr2d)
-----

                        Changes in version 1.8.1                        

- added grid to Imports

- fixed citation

- updated DOI

[IgGeneUsage](/packages/IgGeneUsage)
-----------

                 Changes in version 1.9.3 (2022-03-31)                  

- input data checks extended

[illuminaio](/packages/illuminaio)
----------

                 Changes in version 3.24.0 (2021-05-19)                 

- The version number was bumped for the Bioconductor release version,
  which is
  now BioC 3.13 for R (>= 4.0.3).

                 Changes in version 3.22.0 (2020-10-27)                 

- The version number was bumped for the Bioconductor release version,
  which is
  now BioC 3.12 for R (>= 4.0.0).

                 Changes in version 0.37.0 (2021-10-27)                 

- The version number was bumped for the Bioconductor devel version,
  which is
  now BioC 3.15 for R-devel.

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.1.9 (2022-04-19)                  

- Bug fix: patchDetection now works on SpatialExperiment object

- Fix for release

                 Changes in version 1.1.8 (2022-04-03)                  

- Adjusted read_cpout function to newest pipeline

                 Changes in version 1.1.7 (2022-03-30)                  

- readSCEfromTXT: only read in last metal occurrence in .txt file name

                 Changes in version 1.1.6 (2022-01-12)                  

- Bug fix: correctly index graphs when reading in steinbock data

                 Changes in version 1.1.5 (2022-01-07)                  

- More specific access of metal tags

                 Changes in version 1.1.4 (2021-12-23)                  

- Added option to restrict the maximum distance between neighbors in
  delaunay triangulation

                 Changes in version 1.1.3 (2021-12-15)                  

- Bug fix: divide number of interactions by total number of cells for
  classic and patch interaction counting

                 Changes in version 1.1.2 (2021-11-28)                  

- Moved all helper functions to single script

                 Changes in version 1.1.1 (2021-11-16)                  

- Bug fix in testInteractions: consider now all possible combinations
  of cell types

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.27.6                        

- CHANGES
  * bugfix in default scales

                       Changes in version 1.27.5                        

- NEW FEATURES
  * Default-Scale for immunoMeta-objects

                       Changes in version 1.27.4                        

- CHANGES
  * problems with immunoMeta-objects with 1 cluster

                       Changes in version 1.27.3                        

- CHANGES
  * fixed problems with NA values in FCS data
  * default BD plot scales removed from meta.process

                       Changes in version 1.27.2                        

- CHANGES
  * class and parameter information added to methods weights and mu

                       Changes in version 1.27.1                        

- NEW FEATURES
  * added method cells for immunoClust-object
  * improvement in plot.immunoMeta
  * for.sample option in events.immunoClust

[infercnv](/packages/infercnv)
--------

                 Changes in version 1.11.2 (2022-02-03)                 

- Sort observation groups names when storing the list of indices so
  they are plotted in order, making it easier to change the sorting on
  the final figure.

- Fix plotting error on Windows by adding a check that bitmapType
  option is not null before checking what it is to prevent error in
  comparison.

- Fix for some group labels being cut off in width at the bottom of the
  figure

- Add conversion of expression data in sparse matrix format to dense
  matrix when splitting references in n groups as parallelDist does not
  handle them.

                 Changes in version 1.11.1 (2021-11-08)                 

- Link "plot_chr_scale" and "chr_lengths" options from plot_cnv() to
  run().

- Added parameter k_obs_group to plot_per_group so that it can be
  applied to each of the subsequent group plotting calls if desired.

- Added a check to plot_cnv when using the dynamic_resize option that
  the increased size is not above the max size allowed when using Cairo
  as the graphical back end on Linux. If it is, set the size to the
  highest allowed instead.

- Replaced all calls to dist() with parallelDist(num_threads).
  (suggestion from @WalterMuskovic)

- Updated the code to run the Leiden clustering to be able to handle
  bigger datasets (work around R not having long vectors implemented)
  and be more efficient.

- Add a warning to infercnv object creation if the number of cells is
  more than half of the setting for scientific notation in R, as it may
  cause an issue while using as.hclust(phylo_obj) in the Leiden
  subclustering step.

- Updated plot_cnv method and imports to handle sparse matrices.

- Fix main heatmap drawing which in some cases caused the plotting to
  be out of field in the output.

- Fix method that compares arguments with backups with changes in R
  4.1.1

                 Changes in version 1.10.1 (2021-11-08)                 

- Fix missing colnames in subcluster information when using the Leiden
  method (used downstream by add_to_seurat).

- Fix "plot_per_group" to handle infercnv objects with NULL clustering
  information (mainly to be able to plot using existing results but
  changing the annotations).

[InPAS](/packages/InPAS)
-----

                        Changes in version 2.3.1                        

- merge changes made by Haibo.

[InteractiveComplexHeatmap](/packages/InteractiveComplexHeatmap)
-------------------------

                        Changes in version 1.3.1                        

- depends on a more recent version of ComplexHeatmap.

- fixed the control icons when `compact = TRUE`.

- add two columns of "row_label" and "column_label" in the output of
  `selectPosition()`,
  `selectArea()`, `selectByLabels()`.

- `makeInteractiveComplexHeatmap()`: add two new arguments:
  `show_cell_fun`/`show_layer_fun`
  which controls whether show graphics made by cell_fun/layer_fun on
  the main heatmap.

- `default_click_action()`: numbers show three non-zero digits.

- `htShiny()`: add `app_options` argument.

[InterCellar](/packages/InterCellar)
-----------

                 Changes in version 2.2.0 (2022-04-06)                  

- Updated paper citation in README file

[IntEREst](/packages/IntEREst)
--------

                       Changes in version 1.20.0                        

NEW FEATURES

- limitRanges parameter in interest() and interest.sequential()
  functions allow for targetted analysis. It only loads reads
  that map to the dfined coordinates.

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.30.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Like the DataFrame class defined in the S4Vectors package, classes
  SimpleDataFrameList, CompressedDataFrameList,
  SimpleSplitDataFrameList,
  and CompressedSplitDataFrameList, are now virtual. This completes the
  replacement of DataFrame with DFrame announced in September 2019.
  See:
  https://www.bioconductor.org/help/course-materials/2019/BiocDevelForum/02-DataFrame.pdf

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.5.4 (2022-04-20)                  

MAJOR CHANGES

- ISAnalytics has now a new “dynamic vars system” to allow more
flexibility on user inputs, view the dedicated vignette with
vignette("workflow_start", package="ISAnalytics")
- All package functions were reviewed to work properly with this
system

NEW FEATURES

- gene_frequency_fisher() is a new function of the analysis family
that allows the computation of Fisher’s exact test p-values on gene
frequency - fisher_scatterplot() is the associated plotting function
- top_targeted_genes() is a new function of the analysis family that
produces the top n targeted genes based on the number of IS
- NGSdataExplorer() is a newly implemented Shiny interface that
allows
the exploration and plotting of data
- zipped examples were removed from the package to contain size. To
compensate, the new function generate_default_folder_structure()
generates the standard folder structure with package-included data
on-demand
- transform_columns() is a new utility function, also used internally
by other exported functions, that allows arbitrary transformations
on data frame columns

MINOR CHANGES

- remove_collisions() now has a dedicated parameter to specify how
independent samples are identified
- compute_near_integration_sites() now has a parameter called
additional_agg_lambda() to allow aggregation of additional columns
- CIS_grubbs() now signals if there are missing genes in the refgenes
table and eventually returns them as a df
- outlier_filter() is now able to take multiple tests in input and
combine them with a given logic. It now also produces an HTML
report.
- Several functions now use data.table under the hood
- Color of the strata containing IS below threshold can now be set in
integration_alluvial_plot()

BUG FIXES

- Fixed a minor bug in import_Vispa2_stats() - function failed when
passing report_path = NULL
- Fixed minor issue in circos_genomic_density() when trying to use a
pdf device

DEPRECATED FUNCTIONS

- unzip_file_system() was made defunct in favor of
generate_default_folder_structure()
- cumulative_count_union() was deprecated and its functionality was
moved to cumulative_is()

                 Changes in version 1.5.3 (2022-01-13)                  

MINOR CHANGES

- Added arguments fragmentEstimate_column and
fragmentEstimate_threshold in HSC_population_size_estimate().
Slightly revised filtering logic.
- Updated package logo and website

                 Changes in version 1.5.2 (2021-12-14)                  

NEW (MINOR)

- Added function to check for annotation problems in IS matrices

MINOR CHANGES

- Added argument max_workers in function remove_collisions()
- Updated default functions for aggregate_metadata()
- Added annotation issues section in import matrices report

FIXES

- Fixed minor issue in internals for file system alignment checks
- Fixed minor issue in internal call to import_Vispa2_stats() from
import_association_file()
- Added safe computation of sharing in remove_collisions(): if
process
fails function doesn’t stop

                 Changes in version 1.5.1 (2021-10-28)                  

FIXES

- Attempt to fix issues with parallel computation on Windows for some
plotting functions

[iSEEu](/packages/iSEEu)
-----

                        Changes in version 1.7.1                        

- Fix bug causing .setCachedCommonInfo() to cache NULL value for
valid.rowdata.names information of DynamicMarkerTable panel instead
of empty character vector.

[IsoformSwitchAnalyzeR](/packages/IsoformSwitchAnalyzeR)
---------------------

                Changes in version 1.17.04 (2022-01-06)                 

- Update type: minor.

- Fixed a bug in 1.17.03

                Changes in version 1.17.03 (2022-01-06)                 

- Update type: minor.

- Version bump due to correction in stable branch causing the 1.15 ->
  1.17 bump.

- Fixed date for the last update.

- Fixed a problem with the use of pairwiseAlignment in
  analyzeSwitchConsequences() that could cause jaccard similarities to
  be somewhat wrong.

- Fixed a problem with the switchPlot and transcriptPlot where the
  color of the transcript would be grey instead of red.

- Updates which prepare IsoformSwitchAnalyzeR for future updates

                Changes in version 1.17.02 (2021-10-01)                 

- Update type: minor.

- Various error message updates

- Fixed a problem where importGTF() could have seqLevel problems after
  removals.

- addORFfromGTF() was updated to give better error messages.

                Changes in version 1.17.01 (2021-09-01)                 

- Update type: minor.

- Version bump due to Bioconductor release.

- preFilter() now applies the gene expression cutoff to both conditions
  instead of the overall average.

- analyzePFAM() was updated to reflect recent updates to the tidyverse
  read_fwf function. It furhtermore now better distinguishes tap
  seperated and fixed with files.

[karyoploteR](/packages/karyoploteR)
-----------

                        Changes in version 1.21.0

BUG FIXES

- Fixed wrong positioning of chromosome names in some edge cases (github
  issue #114).

- kpAddLabels now work on zoomed in regions (github issue #112).

[kebabs](/packages/kebabs)
------

                       Changes in version 1.29.1                        

- fixed C++ code to work with newer version of Rcpp (enforced
  STRICT_R_HEADERS)

- updated URLs and DOIs (now requires R version >= 3.3.0)

                       Changes in version 1.29.0                        

- new branch for Bioconductor 3.15 devel

[LACE](/packages/LACE)
----

                         Changes in version 2.0                         

- Implementation of a more complete framework.

- Implementation of dynamic user interface.

[limma](/packages/limma)
-----

                       Changes in version 3.52.0                        

- 
  New function readSampleInfoFromGEO().

- 
  Allow the `trend` argument of eBayes() to specify a general
  variance covariate.

- 
  More detailed checking and error messages when the `object`
  input to lmFit() is a data.frame.

- 
  Improve checking for incorrectly specified `contrasts` argument
  in contrasts.fit().

[lineagespot](/packages/lineagespot)
-----------

                 Changes in version 0.99.0 (2021-09-21)                 

- Submitted to Bioconductor

- Added a NEWS file to track changes of the package

[LinTInd](/packages/LinTInd)
-------

                       Changes in version 0.99.1                        

Bugfixes

- Define the default parameters in IndelIdents()
- Cancels reading of unnecessary dependent packages

                       Changes in version 0.99.0                        

First built

- Submitted to Bioconductor

                     Changes in version 0.0.0.9000                      

NEWS.md setup

- added NEWS.md creation with newsmd

[maftools](/packages/maftools)
--------

                       Changes in version 2.12.00                       

(GitHub master branch/BC 3.15 RC)

NEW FUNCTIONS

- gtMarkers, prepAscat, prepAscat_t, segmentLogR provides interface
to
copy number analysis with ASCAT.
- plotMosdepth and plotMosdepth_t processes output generated by
Mosdepth and performs copy number analysis with DNAcopy
- cancerhotspotsAggr Aggregates cancerhotspots reports

ENHANCEMENTS

- Added support for plotting adjusted p-values in
somaticInteractions.
Issue: 813
- Added support for the protein structure BECN2 Issue: 696

[MAGeCKFlute](/packages/MAGeCKFlute)
-----------

                       Changes in version 1.99.0                        

- Avoid downloading reference files

- Perform enrichment analysis based on human pathways by default

[MassSpecWavelet](/packages/MassSpecWavelet)
---------------

                 Changes in version 1.61.3 (2022-04-05)                 

- Add signal to Suggests. Replace sav.gol() with
signal::sgolayfilt().
Closes #2
- Minor documentation improvements

                 Changes in version 1.61.2 (2022-04-04)                 

- Remove xcms and caTools. Those are great packages but we don't use
them directly in MassSpecWavelet.
- Replace Sweave with RMarkdown vignette and update styles with
BiocStyle

                 Changes in version 1.61.1 (2022-04-01)                 

- Change Maintainer to Sergio Oller
- Add NEWS.md file
- Fix warning (error on R>=4.2) when cwt() has a matrix in the
wavelet
argument
- Fix warning due to partial matching of arguments inside cwt()
- Move waveslim from Depends to Suggests
- Make all calls to recommended packages qualified (e.g. sd ->
stats::sd)
- Remove empty sections in man/ files.
- Register C routines

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.3.8 (2022-04-12)                  

- change unit tests after changing namespace in 1.3.7

                 Changes in version 1.3.7 (2022-04-12)                  

- change namespace in module_measuredValues_missingValues.R

                 Changes in version 1.3.6 (2022-04-08)                  

- set internally parameters (probs, batchColumn) in normalizeAssay and
  batchCorrectionAssay to default values when the parameters are NULL

                 Changes in version 1.3.5 (2022-04-04)                  

- fix error message in cvFeaturePlot after updating packages

                 Changes in version 1.3.4 (2022-02-10)                  

- fix bug in hoeffDValues for pivot_wider after updating tidyr
  (version 1.2.0)

                 Changes in version 1.3.3 (2022-01-25)                  

- harmonize clustering method in distShiny for columns/rows

                 Changes in version 1.3.2 (2021-12-09)                  

- add import of txt and xlsx files for function maxQuant

- change rounding in mosaic that the plot shows more detailed numbers

                 Changes in version 1.3.1 (2021-12-01)                  

- use make.names for character vectors and colnames(colData(se)) in the
  functions for dimension reduction plot, drift plot, ECDF plot, mosaic
  plot,
  features along variable histogram, UpSet plot

- add method "log" in function transformationAssay

- add function spectronaut to upload spectronaut files

[MBECS](/packages/MBECS)
-----

                       Changes in version 0.99.13                       

- Adjusted SVD code to more straightforward implementation.

- Fixed variable ordering for LM-variance calculations.

                       Changes in version 0.99.12                       

- Reduced size of data-sets in testing procedures to stop
  unit-test timeouts.

- Added examples to data objects and reporting functions.

                       Changes in version 0.99.11                       

- Cleaned repository.

- GitHub Readme includes Installation instructions and brief
  workflow tutorial.

                       Changes in version 0.99.10                       

- Added BiocViews: ReportWriting, Visualization, Normalization
  and QualityControl.

- Inlcuded 'importFrom' for all ggplot2 functionality.

- Required packages now live in Imports section.

- Switched to 'aes_string' to remove visible binding NOTEs in
  plot functions.

- Included local variables to shut up the remaining NOTEs.

- Moved 'match.arg' choices to function heads.

- Moved code for generation of mock-up dummy-data from data-raw
  into the data.R file.

- Added testing for dummy data.

- Fixed highlighting in vignette.

- Fixed gray-area issue in box-plots.

                       Changes in version 0.99.8                        

- Included NEWS file.

- Added Bioconductor installation instructions in vignette.

- Added package man-page.

- Set lazyData to false.

- Added URL and BugReports fields to Description.

- Remove bapred and permute packages from suggest because they
  are no longer required.

- Also removed pals package, but included reference because I use
  the tableau color-scheme.

- All packages that are required for execution are now 'Depends'
  instead of 'Suggests'.

- Revised vignette content, formatting and spelling for more
  convenient user experience.

- Fixed test for percentile normalization to actually test
  PN-function.

- Prelim report now checks if clr/tss transformed values are
  present prior to calculating them.

- Report-names/directories can now be changed from default by
  user.

- Reporting functions included in testing.

- Formatted code to adhere to 4-space indentation and 80
  characters width requirements. For the most part.

- Included code for generation of dummy-data.

[mbOmic](/packages/mbOmic)
-----

                        Changes in version 0.99.3

UPDATE MAIN CLASSES
- remove mbSet S4 class
- add Set virtual class
- add bSet and mSet extend to Set class

ADD NEW FUNCTIONS
- add function to identify enterotype

[memes](/packages/memes)
-----

                        Changes in version 1.2.5                        

- fixed error in ame_compare_heatmap_methods that triggers when
plotting without providing a group argument.

                        Changes in version 1.2.4                        

- updated NAMESPACE to fix R CMD CHECK note.

                        Changes in version 1.2.3                        

- fixed a bug in importTomTomXML where tomtom list column would
contain missing data if tomtom was run using multiple database
sources as input.

                        Changes in version 1.2.2                        

- fixed a bug in runStreme causing failures on data import for STREME
version >= 5.4.1

                        Changes in version 1.2.1                        

- fixed a bug in runStreme causing failures when using
BStringSetLists
as input

[meshr](/packages/meshr)
-----

                        Changes in version 2.0.2                        

- CITATION was updated

[metabCombiner](/packages/metabCombiner)
-------------

                        Changes in version 1.5.2                        

- Bug Fix in fit_model():
  - rtx & rty parameter issue resolved

                        Changes in version 1.5.1                        

- New function updateTables():
  - user changes to combinedTable report
  - inclusion of non-intersected features

- Changes in metabCombiner():
  - handling of missing features enabled (i.e. when update() is used)
  - new argument "impute"
  - rtOrder argument bug fix

- Changes in metabCombine():
  - new arguments "union" & "impute"

- Changes in fit_model/ fit_gam()/ fit_loess():
  - new arguments rtx & rty

- Changes in batchCombine():
  - new argument "union"
  - expected multiple (2 or more) datasets as input
  - end message added

- Changes in labelRows()/ reduceTable()/ labelRowParam()/
  reduceTableParam():
  - new argument 'useID'
  - 'resolveConflicts' and 'remove' set to TRUE in reduceTable()

- Changes to metabData objects:
  - "extra"" column count added to show() message

- Package functions updated to handle "group 0"

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                        Changes in version 0.99                         

Changes in 0.99.15

- Highlight query and target spectra in different colors for
validateMatchedSpectra.
- query and/or target of type SummarizedExperiment supported for
Matched objects.
- MatchedSummarizedExperiment class removed.
- query and/or target of type QFeatures supported for Matched
objects.
- Support SummarizedExperiment and QFeatures for both query and
target
parameters in matchValues.

Changes in 0.99.14

- Improve plotly-based mirror plots in validateMatchedSpectra.

Changes in 0.99.13

- Fix issue about matchedData not working for result objects of
matchValues, Mz2MassParam and matchValues, Mz2MassRtParam (issue
#69).

Changes in 0.99.12

- Update plotly-based mirror plots in validateMatchedSpectra.

Changes in 0.99.11

- Change matchMz into matchValues (issue #65).

Changes in 0.99.10

- Add validateMatchedSpectra for manual inspection and validation of
an MatchedSpectra object.

Changes in 0.99.9

- Add setBackend for MatchedSpectra objects.

Changes in 0.99.8

- Add matchMz, Mz2MassParam and matchMz, Mz2MassRtParam. (issue #56).

Changes in 0.99.7

- Add formula matching functions.

Changes in 0.99.5

- Add parameter ... to plotSpectraMirror.
- Definitions of "score", "score_rt" changed to be the difference
(with sign) between query and target m/z or retention time
respectively.
- "ppm_error" becomes error without sign.

Changes in 0.99.4

- Add matches m/z error (variable "ppm_error") to the Matched object
returned by matchMz.

Changes in 0.99.3

- Address Herve's comments.

                         Changes in version 0.2                         

Changes in 0.2.11

- Fix calculation of correct number of rows/columns of the plot in
plotSpectraMirror.

Changes in 0.2.10

- Add parameter toleranceRt to CompareSpectraParam to enable
retention
time-based pre-filtering (issue #35).

Changes in 0.2.9

- Add support for manually defined adducts to Mass2MzParam (issue
#41).

Changes in 0.2.8

- Add parameter THRESHFUN_REVERSE to MatchForwardReverseParam to
allow
filtering results on forward and reverse score (issue #37).

Changes in 0.2.7

- Performance improvement in matchSpectra if no precursor m/z filter
is used (issue #38).
- Report number of matching peaks in
matchSpectra,MatchForwardReverseParam (issue #36).

Changes in 0.2.6

- Fix bug in matchSpectra that was wrongly calculating the acceptable
m/z difference if tolerance was > 0 (issue #34). Fix proposed by
Hugo Varet (@hvaret).

Changes in 0.2.5

- Improve performance of matchMz.
- Rename queryColumn and targetColumn to queryColname and
targetColname.

Changes in 0.2.4

- Support data.frame, DataFrame and matrix in matchMz.
- Add addMatches and filterMatches functions.

Changes in 0.2.3

- Fixes in MatchedSpectra.

Changes in 0.2.2

- Add MatchedSummarizedExperiment.

Changes in 0.2.1

- Rename TargetMass2MzParam to Mass2MzParam.

Changes in 0.2.0

- Add support for matching m/z against m/z and m/z in addition to
retention times to matchMz.

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                         Changes in version 1.3                         

MetaboCoreUtils 1.3.8

- Support for heavy isotopes in
countElements/pasteElements/calculateMass (issue #53).

MetaboCoreUtils 1.3.7

- Fix bug in containsElements function (issue #51).

MetaboCoreUtils 1.3.6

- Add functions for Kendrick mass defects.

MetaboCoreUtils 1.3.5

- Add missing unit tests.

MetaboCoreUtils 1.3.4

- Add calculateMass function.

MetaboCoreUtils 1.3.3

- Vectorized versions for chemical mass functions.

MetaboCoreUtils 1.3.2

- Function for conversion of migration time to effective mobility in
CE-MS.

MetaboCoreUtils 1.3.1

- Add isotope detection functionality.

[MetNet](/packages/MetNet)
------

                 Changes in version 1.13.2 (2022-03-04)                 

- add function addSpectralSimilarity and allow to add a MS2 similarity
  matrix (contribution by Liesa Salzer)

- adjust the functions threshold and combine to be able to deal with
  MS2 similarities (contribution by Liesa Salzer)

- adjust the vignette to the changes imposed by the new function
  addSpectralSimilarity (contribution by Liesa Salzer)

                 Changes in version 1.13.1 (2022-02-11)                 

- update unit tests, e.g. remove ggm for as.data.frame, set R=1000 for
  bayes

[mia](/packages/mia)
---

                         Changes in version 1.3                         

- name change: testForExperimentCrossCorrelation to
  testExperimentCrossCorrelation

- getExperimentCrossCorrelation: Filtering disabled by default, option
  to suppress warnings

- bugfix: taxonomyTree gave error if taxa were agglomerated at highest
  level (taxa name mismatch)

- bugfix: subsampleCounts errors if no samples are found after
  subsampling

- added loadFromMetaphlan

- renamed calculateUniFrac to calculateUnifrac

- added na.rm option to getTopTaxa function

- bugfix: makeTreeSEFromPseq -- orientation of assay is taken into
  account

- bugfix: getExperimentCrossCorrelation's "matrix"" mode works with
  features named equally

- bugfix: getExperimentCrossCorrelation's calculates correlations
  correctly with features named equally

- getExperimentCrossCorrelation name changed to
  getExperimentCrossAssociation

- getExperimentCrossAssociation: user's own function supported, sort in
  mode == table enabled

- getExperimentCrossAssociation: added MARGIN & paired options,
  efficiency of algorithm improved

[miaViz](/packages/miaViz)
------

                         Changes in version 1.3                         

- Bugfix for plotting functions to support sparseMatrix and other assay
  types (2021-12-31)

[microbiome](/packages/microbiome)
----------

                 Changes in version 1.17.2 (2022-02-15)                 

- Bug fix error in transform when taxa_are_rows is FALSE

- Merge microbiomeutilites functionality
  * Convert phyloseq slots to tibbles
  * Combine_otu_tax joins otu_table and tax_table
  * Merged add_besthit a fine tuned version from format_to_besthit
  * Merged psmelt2 a fine tuned version from phy_to_ldf
  * Bug fix in transform method alr

- alr transformation added

                 Changes in version 1.17.1 (2022-01-11)                 

- Fixed bug in plot_core

- bfratio function removed

[microbiomeMarker](/packages/microbiomeMarker)
----------------

                        Changes in version 1.1.2                        

- Development version on Bioconductor
- Use 3rd version of testthat to fix test error (use expect_snapshot()
  rather than expect_known_ouput).
- Add two new arguments in plot_heatmap() scale_by_row and
  annotation_col to improve heatmap viaualization, #52.
- Set slot marker_table to NULL if no marker was identified.
- Add new import function import_picrust2() to import prediction
  functional table from PICRUSt2, and all DA functions support for
  PICRUSt2 output data.
- Keep color consistent between legend and plot in cladogram, #42.
- Add a new argument clade_label_font_size in plot_cladogram() to
  specify font size of clade label, #49.

                 Changes in version 1.1.1 (2020-03-07)                  

- Add a para only_marker in plot_cladogram to specify whether only
show the markers or all features in the cladogram.
- Fix a bug in run_test_multiple_groups(), error group names for
enrich groups (2021-10-12, #48).
- Fix a bug in plot_abundance(), error var name of effect size in
marker_table (2021-10-17, #47).

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.7.11                        

- add mp_plot_diff_cladogram to plot the result of mp_diff_analysis.
(2022-04-19)

                       Changes in version 1.7.10                        

- optimizing the mp_aggregate_clade and mp_balance_clade.
(2022-04-13)

                        Changes in version 1.7.9                        

- add mp_cal_pd_metric to calculate the related phylogenetic
diversity
metrics. (2022-04-02) including NRI, NTI, PD, PAE, HAED, EAED, IAC.
- add mp_balance_clade to calculate the balance score of internal
nodes according to their tip nodes abundances. (2022-03-22)
- add extract_binary_offspring to find the descendant
tip/internal/all
(with type parameter) nodes. (2022-03-17)
- add mp_aggregate_clade and mp_diff_clade to calculate and test the
abundance (differential signals) of internal node according to their
tip nodes abundance. (2022-03-16)
- add mp_select_as_tip and fix the bug of mp_diff_analysis with
specific tip.level (not OTU) argument. (2022-03-02, Mon)
- fix the replace_na bug of new tidyr. (2022-03-04, Fri)
- update mp_import_metaphlan to better parse the output of
MetaPhlAn2.
(2022-03-11, Fri)
- update the mp_cal_abundance to return the tbl_df contained numeric
type sample metadata. (2022-03-11, Fri)

                        Changes in version 1.7.8                        

- supporting multiple group names and supporting numeric type for
.group of mp_plot_alpha. (2022-02-15, Tue)
- supporting multiple group names for .group of mp_plot_abundance
when
plot.group=TRUE. (2022-02-14, Mon)
- fix the width of mp_plot_abundance with geom="flowbar".
(2022-02-01,
Tue)

- related issue

                        Changes in version 1.7.7                        

- fix the bug about the constant variables within groups in lda of
MASS (2022-01-27, Thu)

                        Changes in version 1.7.6                        

- fix the issue, that kingdom level of taxonomy information contains
k__ or K__, which is unknown annotation in kingdom. (2021-01-14,
Fri)
- add mp_extract_taxatree and mp_extract_otutree (alias of
mp_extract_tree). (2022-01-14, Fri)

                        Changes in version 1.7.5                        

- add the message for the not integers in mp_cal_alpha. (2021-12-31,
Fri)
- remove the features which variance of their abundance is zero
before
identify different taxa. (2021-12-28, Tue)
- add bar option in mp_plot_abundance, default is flowbar, the other
options are bar and heatmap
- use corrected relative eigenvalues when the eigenvalues has
negative
values. (2021-12-27, Mon)
- add new distmethod from hopach. (2021-12-27, Mon)
- update tax_table without required phyloseq. (2021-12-20, Mon)
- update mp_diff_analysis to support the factor type group (.group
specified). (2021-12-20, Mon)

                        Changes in version 1.7.4                        

- update taxatree<- and otutree<- which will extract the intersection
between the tip labels of input treedata and the rownames of MPSE.
(2021-12-14, Tue)
- add taxonomy<- for MPSE to assign the taxonomy information, which
will be converted to taxatree automatically. (2021-12-14, Tue)
- add tax_table for MPSE and return taxonomyTable defined in
phyloseq.
(2021-12-14, Tue)
- update mp_import_metaphlan to better parse the output of
MetaPhlAn2.
(2021-11-30, Tue)

                        Changes in version 1.7.3                        

- update 'mp_plot_abundance' (2021-11-24, Wed)
- support heatmap by setting geom.
- .group supports multiple characters and .sec.group will be
removed in the next version.
- update mp_plot_diff_res (2021-11-23, Tue)
- support otutree and taxatree class by setting tree.type.
- support multiple layout types of tree by setting layout.
- support adjusting the gap between panel and width of panel by
setting offset.abun, pwidth.abun, offset.effsize, pwidth.effsize
- support whether display the relative abundance of group instead
of sample by setting group.abun=TRUE or sample number > 50
- add mp_plot_diff_res to visualize the result of mp_diff_analysis.
(2021-11-22, Mon)

                        Changes in version 1.7.2                        

- speed up the mp_cal_abundance, mp_cal_venn and mp_cal_upset with
dtplyr. (2021-11-18, Thu)
- update the guide of x axis of ggside in mp_plot_ord. (2021-11-15,
Mon)
- update mp_plot_abundance to visualize the abundance of taxonomy
from
high (bottom) to low (top). (2021-11-15, Mon)
- support multiple annotation rows or cols of heatmap of mp_plot_dist
with .group=c(group1, group2), and add set_scale_theme to adjust the
scale or theme of subplot of heatmap. (2021-11-10, Wed)
- fix the issue when the taxonomy info is removed with select.
(2021-11-09, Tue)
- update print for MPSE class. (2021-11-09, Tue)
- update otutree<- for support phylo class. (2021-11-09, Tue)
- speed up the integration of mp_cal_dist result with action="add".
(2021-11-09, Tue)
- update as.MPSE for biom class to support parsing the metadata of
sample. (2021-11-09, Tue)

                        Changes in version 1.7.1                        

- fix the issue when using filter only return a assays contained one
feature (nrow=1). (2021-11-05, Fri)
- fix the error of rownames<- when rownames of MPSE is NULL.
(2021-11-04, Thu)
- update 'message' or 'stop error message' when the 'Abundance'
cannot
be rarefied in some functions, such as mp_cal_alpha, mp_cal_venn,
mp_cal_upset, mp_cal_abundance and mp_cal_NRT_NTI. (2021-10-29, Fri)
- introduce .sec.group argument to specify the second group name in
mp_plot_abundance, if it is provided, the nested facet will be
displayed. (2021-11-02, Tue)

[miloR](/packages/miloR)
-----

                 Changes in version 1.3.1 (2022-01-07)                  

- Fix bug in findNhoodGroupMarkers to merge on gene IDs explicitly
- Fix bug in makeNhoods to include index cell in nhoods() matrix
- Introduce graph-based neighbourhood definition - allows full
compatibility with graph-only batch correction and graphs
constructed by third-party tools
- Introduce graph-based spatial FDR correction to obviate the need
for
any distance calculations
- Add vignette to describe the use and application of contrasts in
testNhoods
- Patch to correct SpatialFDR with sparse nhoods where density is ~0

[mirTarRnaSeq](/packages/mirTarRnaSeq)
------------

                 Changes in version 1.3.2 (2022-04-06)                  

- suggest SPONGE rather than import

                 Changes in version 1.3.1 (2022-03-26)                  

- added support for SPONGE

[mistyR](/packages/mistyR)
------

                        Changes in version 1.4.0                        

- Release version for Bioconductor 3.15. See changes for 1.3.x.

                         Changes in version 1.3                         

- Switched from Louvain to Leiden algorithm for community detection
(requires igraph >= 1.2.7).
- The metamodel is now build by ridge regression. Intercept p-values
are not calculated, the values are set to NA for backwards
compatibility.
- Unique value error for cv folds is downgraded to warning.
- Prefix can be added to the column names generated by the fucntions
add_juxtaview and add_paraview. This allows modeling the maker
expression by its own juxtaview and paraview.
- Bug fixes.

                        Changes in version 1.2.1                        

- Fixed a separator issue in results aggregation and signature
generation that might cause issues with variable names containing
"_".

[MLP](/packages/MLP)
---

                       Changes in version 1.43.1                        

- getGeneSets: add 'Rhesus' to 'species' (available only for 'GOBP',
  'GOMF', 'GOCC' or 'KEGG')

[MobilityTransformR](/packages/MobilityTransformR)
------------------

                        Changes in version 0.99                         

- Function for conversion of migration time to effective mobility in
CE-MS.

- Submitted to Bioconductor

[monaLisa](/packages/monaLisa)
--------

                        Changes in version 1.1.2                        

- added citation to the Bioinformatics publication

                        Changes in version 1.1.1                        

- added link to pre-print manuscript on biorXiv to README.md
- added warning to bin(..., minAbsX = val) if adjusted zero-bin
breaks
deviate more than 20% from val
- added doPlot argument to plotMotifHeatmaps to select if heatmaps
should be plotted or just generated and returned
- added LICENSE.md file
- expanded monaLisa.Rmd vignette with illustration on how to do a
binary or single set motif enrichment analysis
- expanded on collineairty in regression in the
selecting_motifs_with_randLassoStabSel.Rmd vignette and the choice
of parameter values in stability selection.
- updated the results.binned_6mer_enrichment_LMRs.rds and
results.binned_motif_enrichment_LMRs.rds files stored in monaLisa
under the current version of the package.

[Motif2Site](/packages/Motif2Site)
----------

                       Changes in version 0.99.5                        

- Resolve problem with overwriting BiocStyle in vignette

                       Changes in version 0.99.4                        

- Remove examples from data.Rd

                       Changes in version 0.99.3                        

- Adding help page data.Rd

                       Changes in version 0.99.2                        

- Update readme by Bioconductor installation instructions using
BiocManager
- Update description by adding R>=4.1 as dependency
- Update description by adding BugReports link
- Update description by changing lazyDate to false
- Update Vignettes by fixing typos, abstract, and author
- Update Vignettes by using BiocStyle and keeping code lines under 80
character
- Improve coding style
- Fixing some of the build notes
- Adding Motif2Site.Rd with a complete example

                       Changes in version 0.99.1                        

- All HOX-mouse samples are replaced with FUR-Ecoli ones
- extdata folder size is 4 MB as E. coli data is much smaller
- Vignettes was updated by E. coli data
- Examples were updated by E. coli data
- README were updated by E. coli data
- The build time has been decreased substantially

                       Changes in version 0.99.0                        

- The first version of Motif2Site was submitted to Bioconductor

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.39.4                        

- Change the fontfamily from 'mono,Courier' to 'mono' in vignettes.

                       Changes in version 1.39.3                        

- add XMatrix format for importMatrix.

                       Changes in version 1.39.2                        

- Accept scales for y-axis for logo when ic.scale is FALSE.

                       Changes in version 1.39.1                        

- Accept user defined x-axis for logo.

[MouseFM](/packages/MouseFM)
-------

                 Changes in version 1.4.2 (2022-02-28)                  

- Citation adapted

                 Changes in version 1.4.1 (2021-11-13)                  

- Set biomaRt to url of archived version GRCm38 to be compatible with
  the package

- Citation added

[MQmetrics](/packages/MQmetrics)
---------

                        Changes in version 1.3.3                        

- Removed unused dependency.

                        Changes in version 1.3.2                        

- Reporthing the fixed Modifications if the mqpar.xml is present.

                        Changes in version 1.3.1                        

- Added median to the PlotAndromedaScore().
- Remove warnings from the vignettes.
- Now MQmetrics can read the mqpar.xml file. This one must be located
on the same directory as the combined folder.
- Reporting the number of threads used by MaxQuant if found the
mqpar.xml file.

[msa](/packages/msa)
---

                       Changes in version 1.27.2                        

- applied patch to allow msa to work with the new Windows UCRT
  toolchain

                       Changes in version 1.27.1                        

- workaround for problems running texi2dvi() on R 4.2.0; those occurred
  during package checks when running some examples and the vignette
  code

- updated URLs and DOIs (now requires R version >= 3.3.0)

- fixed msaConvert() function to now work well with newer versions of
  the
  'ape' package (now requires at least version 5.2)

                       Changes in version 1.27.0                        

- new branch for Bioconductor 3.15 devel

[MSA2dist](/packages/MSA2dist)
--------

                 Changes in version 0.99.3 (2022-03-18)                 

Major changes

Minor improvements and bug fixes

- Fix RcppThread::LdFlags warning

                 Changes in version 0.99.2 (2022-01-28)                 

Major changes

Minor improvements and bug fixes

- Added RcppThread::ProgressBar

                 Changes in version 0.99.1 (2022-01-27)                 

Major changes

- Changed version number into 0.99.2

Minor improvements and bug fixes

- Changed URL links in DESCRIPTION

                 Changes in version 0.99.0 (2021-12-22)                 

Major changes

- Changed version number into 0.99.1

- Changed name from distSTRING into MSA2dist

- Submitted to Bioconductor

Minor improvements and bug fixes

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                         Changes in version 1.3                         

Changes in 1.3.5

- Add parameter columns to peaksData.

Changes in 1.3.4

- Import and use spectraVariableMapping method from Spectra.

Changes in 1.3.3

- Map comment to spectra variable.

Changes in 1.3.2

- Use in addition unit tests from the Spectra package.
- Add filterPrecursorMzValues and filterPrecursorMzRange.

Changes in 1.3.1

- MsBackendMassbankSql extends Spectra::MsBackendCached to re-use the
general caching mechanism provided by that backend.

[MsBackendMgf](/packages/MsBackendMgf)
------------

                         Changes in version 1.3                         

Changes in 1.3.3

- Import coreSpectraVariables from Spectra.

Changes in 1.3.2

- Adapt the spectraVariableMapping to Spectra version >= 1.5.8.

Changes in 1.3.1

- Run additional unit test suits defined in Spectra.

[MsBackendMsp](/packages/MsBackendMsp)
------------

                        Changes in version 0.99                         

Changes in 0.99.4

- Import coreSpectraVariables from Spectra.

Changes in 0.99.2

- Small updates and changes.

Changes in 0.99.1

- Address review comments.

                        Changes in version 0.98                         

Changes in 0.98.2

- Vignette added.

Changes in 0.98.1

- Improve export method to support multi-value fields and
user-provided mappings.

Changes in 0.98.0

- Add additional spectra variable mappings.
- Call unit tests from the Spectra package.

[MsBackendRawFileReader](/packages/MsBackendRawFileReader)
----------------------

                         Changes in version 1.1                         

- Call the Spectra test suite.

- Add export snippet for MGF export in vignette file.

- Add .top_n peak list filter example in vignette file.

- Add test case for comparing .top_n(n=10) function with Proteome
  Discoverer Software v2.5 workflow using scan 9594.

- Change rawrr dependency to v1.3.5 to benefit from monoisotopic
  mZ values.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                         Changes in version 1.7                         

MsCoreUtils 1.7.5

- Function bin gains parameter returnMids to choose whether or not
bin
mid-points should be returned in the result list.

MsCoreUtils 1.7.4

- Fix ppm to always return a positive value (issue #94).

MsCoreUtils 1.7.3

- Add citation.

MsCoreUtils 1.7.2

- Use Matrix::colSums() by default to handle sparce 'Matix' and
'matrix' adjacency matrices.

MsCoreUtils 1.7.1

- New aggregate_by_matrix() that uses an adjacency matrix to
aggregate
quantitative features.
- Set colnames to the outputs of aggregate_by_matrix() and
aggregate_by_vector() to make sure that these are always set and not
reply on the underlying function.

MsCoreUtils 1.7.0

- New Bioc devel version

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.21                         

Changes in 2.21.7

- Fix mz calculation in calculateFragments for neutral losses with a
charge > 1 (see issue 573).

Changes in 2.21.6

- Allow different orientation of axis labels in image2 (see issue
571).

Changes in 2.21.5

- Adapt to changes in mzR version 2.29.3.

Changes in 2.21.4

- Add transformIntensity function for Chromatogram and MChromatograms
objects.

Changes in 2.21.3

- Fix bug (issue #561).

Changes in 2.21.2

- Fix bug in compareChromatograms that creates a non-symmetric
similarity matrix.

Changes in 2.21.1

- Change default for MSnbase fast load variable: set to TRUE also on
macOS.

[msPurity](/packages/msPurity)
--------

                       Changes in version 1.21.2                        

- Fixes due to mz and intensity being named in columns for mzR and XCMS

- Fixes for connection{base} file() opening changes (no longer accept
  w+a in function in R v4.2)

- Update tests for the above

- Update a reference of grpid in averageXFragSpectra functions more
  explicit (no change in functionality)

- Typo fix in vignette

                       Changes in version 1.21.1                        

- Bugfix for frag4feature for XCMS 3 compatability
  https://github.com/computational-metabolomics/msPurity/pull/93

- Remove imports that are no longer used

[MSstatsConvert](/packages/MSstatsConvert)
--------------

                        Changes in version 1.5.1                        

- Fixed bugs related to data types.
- Fixed bug that resulted in missing TechReplicate column in output.
- Added Philosopher converter.
- Added methods for MSstatsValidated objects.

[MSstatsTMT](/packages/MSstatsTMT)
----------

                 Changes in version 2.2.7 (2022-02-18)                  

- Minor change: extend PhilosophertoMSstatsTMTFormat function to
  have multiple types of input

                 Changes in version 2.2.6 (2022-02-14)                  

- Major change: add PhilosophertoMSstatsTMTFormat function as
  converter for outputs from Philosopher

                 Changes in version 2.2.5 (2021-10-25)                  

- Minor change: add different point shape to dataProcessPlotsTMT
  as indicator of imputed values

                 Changes in version 2.2.3 (2021-10-06)                  

- Minor change: fix the bug when df.prior is infinite

[MuData](/packages/MuData)
------

                       Changes in version 0.99.0                        

- submitted package to Bioconductor

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.22.0                        

Bug fixes and minor improvements

- Add data("miniACC") to examples after removing lazy loading.
- Class definition prototypes defined for cleaner extensibility
(@hpages, #306).
- Doc and internal improvments to MultiAssayExperimentToMAF
- synAssay and nonSynAssay now require exact assay names in
MultiAssayExperimentToMAF

[multiHiCcompare](/packages/multiHiCcompare)
---------------

                 Changes in version 1.13.0 (2022-04-21)                 

- Match NEWS and DESCRIPTION versioning

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.3.18                        

New features

- Can now handle general remote sumstats not just IEU GWAS
- More column header mappings

                       Changes in version 1.3.17                        

New features

- Clean up of column header mapping file, including FREQUENCY given
priority over MAF and addition of new CHR mappings.

                       Changes in version 1.3.15                        

Bug fixes

- Handle cases for multi-trait GWAS when P columns exists separate to
the trait specific P value so that when renaming occurs there isn't
two P columns. Inputted P column will be renamed to 'P_input'
- Issue where 'check allele flip' wasn't running when the sumstats
had
all SNP IDs missing and incorrect direction of A1/A2 and effect
columns has now been fixed.

                       Changes in version 1.3.14                        

New features

- liftover
- Now exported function.
- Added args for more user flexibility.
- Uses GenomeInfoDb::mapGenomeBuilds to standardise build names.
- Warns users when mapped builds do not match one of the
conversion options.
- Choice to output as data.table or GRanges.
- Added units tests for exported version.
- standardise_sumstats_column_headers_crossplatform
- Exported as standardise_header while keeping the original
function name as an internal function (they call the same code).
- Added unit tests for exported version.
- Added chunks to *Getting started` vignette
- liftover tutorial
- "Quick formatting" of headers and file formats.

Bug fixes

- check_pos_se: Remove extra message() call around string.
- check_signed_col: Remove extra message() call around string.
- write_sumstats
- Added extra round of sorting when tabix_index=TRUE because this
is required for tabix.

                       Changes in version 1.3.13                        

New Features

- Additional mappings for CHR
- Make A1, A2 upper-case

Bug fixes

- Bug fix for dealing with imputing SNP ID when there are indels

                       Changes in version 1.3.11                        

New Features

- MungeSumstats can now handle Indels better. It will:
- Not impute the RS ID of a SNP for an Indel
- Not remove the Indel based on the RS ID not being present in the
SNP ref dataset.
- Not remove the Indel if it has the same base-pair location as a
SNP in the sumstats.
- Can now handle vcfs with extensions .vcf.tsv, .vcf.tsv.gz and
.vcf.tsv.bgz

Bug fixes

- For non-bi-allelic SNP runs, no longer remove duplicated SNPs based
on their base-pair position or their RS ID.

                        Changes in version 1.3.9                        

New Features

- Exported functions. Added examples and unit tests:
- compute_nsize
- standardise_sumstats_column_headers_crossplatform
- formatted_example
- New arguments:
- standardise_sumstats_column_headers_crossplatform: Added arg
uppercase_unmapped to to allow users to specify whether they
want make the columns that could not be mapped to a standard
name uppercase (default=TRUE for backcompatibility). Added arg
return_list to specify whether to return a named list (default)
or just the data.table.
- formatted_example: Added args formatted to specify whether the
file should have its colnames standardised. Added args sorted to
specify whether the file should sort the data by coordinates.
Added arg return_list to specify whether to return a named list
(default) or just the data.table.
- Removed codecode.yml and _pkgdown.yml files (no longer necessary).
- Added Issues templates for Bugs and Feature requests.
- Added .datatable.aware=TRUE to .zzz as extra precaution.
- vcf2df: Documented arguments.
- Made v2 of hex sticker: inst/hex/hex.png

Bug fixes

- Regenerated the gh-pages branch after it accidentally got deleted.
- Remove temporary docs/ folder.
- Updated GitHub Actions.
- Updated Dockerfile so it doesn't run checks (this is now take care
of by the GHA workflow).
- Added Windows-specific folders to .Rbuildignore.
- Made to_GRanges.R and to_VRanges.R file names lowercase to be
congruent with function names.

                        Changes in version 1.3.7                        

Bug fixes

- Bug in checking for bad characters in RSID fixed

                        Changes in version 1.3.6                        

New Features

- Columns Beta and Standard Error can now be imputed. However note
that this imputation is an approximation so could have an effect on
downstream analysis. Use with caution.

                        Changes in version 1.3.5                        

Bug fixes

- Flipping of Odds Ratio corrected (1/OR rather than -1*OR)

                        Changes in version 1.3.4                        

Bug fixes

- Issue downloading chain file resolved

                        Changes in version 1.3.3                        

New Features

- More mappings added to default mapping file.

                        Changes in version 1.3.2                        

Bug fixes

- Previously rsids with characters added (e.g. rs1234567w) would
cause
an error when checking for the rsid on the reference genome. This
has been fixed and the correct rsid will now be imputed from the
reference genome for these cases.

                        Changes in version 1.3.1                        

New Features

- import_sumstats: Create individual folders for each GWAS dataset,
with a respective logs subfolder to avoid overwriting log files when
processing multiple GWAS.
- parse_logs: New function to convert logs from one or more munged
GWAS into a data.table.
- list_sumstats: New function to recursively search for local summary
stats files previously munged with MungeSumstats.
- Added new dataset inst/extdata/MungeSumstats_log_msg.txt to test
logs files.
- Added unit tests for list_sumstats and parse_logs.
- Added new Docker vignette.
- Updated GHA workflows using r_workflows.
- Remove docs/ folder as the website will now be pushed to the
gh-pages branch automatically by new GHA workflow.
- Made documentation in README more clear and concise.
- Added checks for p-values >1 or <0 via args convert_large_p and
convert_neg_p, respectively. These are both handled by the new
internal function check_range_p_val, which also reports the number
of SNPs found meeting these criteria to the console/logs.
- check_small_p_val records which SNPs were imputed in a more robust
way, by recording which SNPs met the criteria before making the
changes (as opposed to inferred this info from which columns are 0
after making the changes). This function now only handles
non-negative p-values, so that rows with negative p-values can be
recorded/reported separately in the check_range_p_val step.
- check_small_p_val now reports the number of SNPs <= 5e-324 to
console/logs.
- Unit tests have been added for both check_range_p_val and
check_small_p_val.
- parse_logs can now extract information reported by
check_range_p_val
and check_small_p_val.
- New internal function logs_example provides easy access to log file
stored in inst/extdata, and includes documentation on how it was
created.
- Both check_range_p_val and check_small_p_val now use #'
@inheritParams format_sumstats to improve consistency of
documentation.

Bug fixes

- Reduced vignette sizes.
- Removed usage of suppressWarnings where possible.
- Deleted old .Rproj file and hidden folder (contained large files).
- Configured .Rproj so it doesn't store large data files.
- Fix badger issues:
https://github.com/GuangchuangYu/badger/issues/34
- Prevent test-index_tabix.R from running due to errors (for now).

                        Changes in version 1.3.0                        

New Features

- Version bump to align with Bioconductor release 3.14.

[muscat](/packages/muscat)
------

                        Changes in version 1.9.3                        

- bug fix in pbDS(): drop samples w/o any detected features,
  otherwise edgeR::calcNormFactors() fails when lib.size 0

                        Changes in version 1.8.1                        

- bug fix in prepSim(): removal of genes with NA coefficients
  was previously not propagated to the dispersion estimates

- bug fix in test-resDR.R: set 'min_cells = 0' to assure that
  everything is being tested, otherwise unit tests could fail

[mzR](/packages/mzR)
---

                       Changes in version 2.29.4                        

- Re-apply fix for compile error on clang by Kurt Hornik, closes #263

- Remove text in DESCRIPTION hinting at the RAMP wrapper for mzData
  removed in 2.29.3

                       Changes in version 2.29.3                        

- Update to Proteowizard 3_0_21263

- Removed RAMP backend, dropping ability to read mzData

- header always returns a data.frame even for a single scan.

                       Changes in version 2.29.2                        

- Cleanup in build files

                       Changes in version 2.29.1                        

- Pwiz backend partially re-written to avoid segfault on macOS
  (https://github.com/sneumann/xcms/issues/422).

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 2.2.0                        

- Added heatmap argument to plot_gene(), plot_region() and
plot_granges(). This adds a read-heatmap to the plot.
- Added cluster_regions() function to perform k-means clustering on a
table of genomic regions based on methylation profile.
- Added median averaging method for trends in plot_gene(),
plot_region() and plot_granges(). This can be changed using the new
avg_method argument, default is mean.
- Added filter_methy() function to create a filtered methylation
file.
- Added region_methy_stats() to obtain average methylation fractions
of specific regions.
- Added methy_to_edger() direct conversion wrapper around
methy_to_bsseq() and bsseq_to_edger().
- Added palette argument to plot_gene(), plot_region() and
plot_granges() to allow custom colour palettes.
- Fixed bsseq_to_edger() failing when regions argument was used.
- Fixed heatmaps not staying in a single column when more than 2
groups were present.

[NanoStringNCTools](/packages/NanoStringNCTools)
-----------------

                 Changes in version 1.3.1 (2022-01-12)                  

- Enable compatibility with Gen 2.5 RCCs

                 Changes in version 1.3.0 (2021-10-26)                  

- Initial Bioconductor devel 3.14 version

[ndexr](/packages/ndexr)
-----

                       Changes in version 1.17.0                        

- **UPDATE: Using the RCX package for working with networks.**
  **Deprecated Functions:**

- *rcx_fromJSON:* `RCX::readJSON()`

- *rcx_toJSON:* `RCX::toCX()`

- *rcx_aspect_toJSON:* `rcx_aspect_toJSON`

- *rcx_new:* `RCX::createRCX()`

- *rcx_asNewNetwork:* `RCX::createRCX()`

- *rcx_updateMetaData:* `RCX::updateMetaData()`

- *print.RCX:* `RCX::print.RCX()`

- *rcx_toRCXgraph:* `RCX::toIgraph()`

- *rcxgraph_toRCX* `RCX::fromIgraph()`

[nnSVG](/packages/nnSVG)
-----

                 Changes in version 0.99.0 (2022-03-02)                 

- version for submission to Bioconductor

[NormalyzerDE](/packages/NormalyzerDE)
------------

                       Changes in version 1.13.2                        

- Fixed leastRepCount setting of zero for statistics report

[nullranges](/packages/nullranges)
----------

                        Changes in version 1.1.4                        

- Needed to drop features that have 0 width after trimming in
bootRanges.
- Made the validity test for bootRanges only look for iter.

                        Changes in version 1.1.1                        

- Change to factor-Rle output for bootRanges to simplify the
downstream plyranges.

[NxtIRFcore](/packages/NxtIRFcore)
----------

                 Changes in version 1.1.1 (2022-01-11)                  

- Bugfix for NxtSE constructor.

- Bugfix for Consistency filter: previously an "average" filter was
  used, such
  that upstream / downstream filter triggers counted for 0.5, whereas
  it added
  1.0 if both up/downstream filters were triggered. From 1.1.1 onwards,
  1.0 is
  added when either upstream or downstream consistency filter is
  triggered.

- Added two new annotation-based filters: Terminus and ExclusiveMXE.
  See
  ?NxtFilter for details

- Annotated retained introns `RI` are defined by any intron that is
  completely
  overlapped by a single exon of any transcript.
  They are calculated as binary events, i.e. as PSI between RI and
  specific
  spliced intron, and do not consider overlapping splice events
  (unlike `IR` events, which are calculated for all other constitutive
  introns)

                        Changes in version 1.1.0                        

- Initial devel release for Bioconductor 3.15

[OGRE](/packages/OGRE)
----

                       Changes in version 0.99.8                        

- added coverage plot

                       Changes in version 0.99.6                        

- added GUI

                       Changes in version 0.99.5                        

- added AnnotationHub support

                       Changes in version 0.99.4                        

- lazy loading to false

[orthogene](/packages/orthogene)
---------

                        Changes in version 1.1.5                        

BUG FIXES

- map_orthologs_babelgene
- Add "Bad credentials" check for piggyback.
- Add use_old as an optional arg so I can switch to more recent
versions of babelgene::orthologs_df if need be.
- Use updated built-in babelgene::orthologs_df by default.
- Throw error if trying to map between two non-human species.
- Filter support==NA mappings by default, not but support>=2 like
babelgene does by default (even when
babelgene::orthologs(min_support = 1)).
- See here for discussion of discrepancies with babelgene
maintainer: https://github.com/igordot/babelgene/issues/2

NEW FEATURES

- Removed aggregate_rows_delayedarray as it wasn't being used and was
far less efficient than the other methods anyway (which are also
compatible with DelayedArray matrices anyway). * New unit tests:
- load_data
- aggregate_mapped_genes(method='stat')
- sparsity

                        Changes in version 1.1.4                        

BUG FIXES

- Remove source_all as it included a library call.

                        Changes in version 1.1.3                        

NEW FEATURES

- Update GHA

BUG FIXES

- Fix failing benchmarking tests.

                        Changes in version 1.1.2                        

BUG FIXES

- convert_orthologs(method="babelgene") now gets gene mappings from
all_genes_babelgene instead babelgene::orthologs (which doesn't seem
to work very well, despite being dedicated for this purpose).
- map_species:
- Avoid running this function redundantly when nested in multiple
layers of other functions.
- common_species_names_dict now return "scientific_name" by
default, instead of "taxonomy_id"
- Match map_species method to whatever method is being used in the
function it's wrapped within, to avoid dropping species due to
naming differences.
- Add "id" column (e.g. "celegans") to all org databases to
enhance their searchability.
- Add map_species_check_args.
- Ensure proper method-specific output_format when passing species to
other functions.

NEW FEATURES

- plot_orthotree: Automated plotting of phylogenetic trees with 1:1
ortholog report annotations. Includes several subfunctions:
- prepare_tree (exported): Read, prune and standardise a
phylogenetic tree.
- gather_images (internal): More robust way to find and import
valid phylopic silhouettes. Will make PR requests to rphylopic
and ggimage/ggtree to include this functionality.
- Added unit tests for report_orthologs, especially when
method="babelgene".
- GitHub Actions:
- Merge both GHA workflows into one, as implemented in templateR.
- Added citation info to README.
- Save all_genes_babelgene ortholog data to orthogene-specific cache
instead of tempdir to avoid re-downloading every R session.

                        Changes in version 1.1.1                        

BUG FIXES

- Made GHA less dependent on hard-coded R/bioc versions.

                        Changes in version 1.1.0                        

NEW FEATURES

- Now on Bioconductor release 3.14.
- Docker containers automatically built and pushed to DockerHub via
GitHub Actions.
- Dockerfile provided to build and check any R package efficiently
with AnVil.
- CRAN checks and Bioc checks run via GitHub Actions.
- Added documentation on using Docker container to README.
- Documentation website now automatically built via GitHub Actions.
- Code coverage tests now automatically run and uploaded via GitHub
Actions.

[PanomiR](/packages/PanomiR)
-------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[pareg](/packages/pareg)
-----

                Changes in version 0.99.22 (2022-02-18)                 

- Submitted to Bioconductor

[pcaExplorer](/packages/pcaExplorer)
-----------

                       Changes in version 2.22.0                        

Other notes

- get_annotation_orgdb() gains an additional argument,
key_for_genenames, which defaults to "SYMBOL". This should not
change the behavior of the function, if not specified, but
accommodates for the use of annotation packages where the
information has been encoded differently (e.g. org.Sc.sgd.db where
the info is contained in the "ORF" column)

Bug fixes

- pcaplot correctly returns the values for the percent of explained
variance, which were correctly displayed on the plot but not stored
as they should in the attribute slot

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.9.2 (2022-04-18)                  

- GUI corrections and improvements (use of `bslib`)

                 Changes in version 1.9.1 (2021-12-18)                  

- BiocCheck format update

[pengls](/packages/pengls)
------

                        Changes in version 1.1.2                        

- Allow penalty.factor to be user supplied

                        Changes in version 1.1.1                        

- Also apply spatial correction to intercept

[phantasus](/packages/phantasus)
---------

                       Changes in version 1.15.7                        

- Using relative paths in counts metadata

                       Changes in version 1.15.5                        

- Fix windows-related bug in tests

                       Changes in version 1.15.3                        

- Major rework of RNA-seq counts support from external hdf5 files (like
  ARCHS4)

- Fix typing bug

                       Changes in version 1.15.1                        

- Option to omit gene version IDs (like in ENSEMBL) before conversion

- Annotations are trasnmitted with type information, which fixes many
  bugs

[PhIPData](/packages/PhIPData)
--------

                 Changes in version 1.3.3 (2022-02-04)                  

- Converted to using `BiocFileCache` for storing aliases and peptide
  libraries.

                 Changes in version 1.3.2 (2022-02-02)                  

- Corrected typo in sample name error.

- Changed `paste0()` calls to `file.path()` calls.

[PhosR](/packages/PhosR)
-----

                        Changes in version 1.5.1                        

- Development version

                        Changes in version 1.4.1                        

- Bug fix in PCA plot variance

[PhyloProfile](/packages/PhyloProfile)
------------

                        Changes in version 1.8.6                        

- option to identify sequence source

                        Changes in version 1.8.5                        

- added functions for exporting plot settings

                        Changes in version 1.8.4                        

- fixed error reading taxonNamesReduced.txt that contains "#"

- added functions for import and export taxonomy DB

                        Changes in version 1.8.3                        

- fixed bug of group comparison function

                        Changes in version 1.8.2                        

- fixed loading cluster from config file

                        Changes in version 1.8.1                        

- do not show pfam links for smart domains and vice versa

[Pigengene](/packages/Pigengene)
---------

                Changes in version 1.21.40 (2022-04-15)                 

Changes in existing functions

- Habil changed the default value from
  hu.mouse(host="useast.ensembl.org", ...) to
  hu.mouse(host="www.ensembl.org", ...) to prevent a check error
  on Bioconductor.

                Changes in version 1.21.36 (2021-11-16)                 

Changes in existing functions

- Habil added the doReturNetworks argument to
  one.step.pigengene().

                Changes in version 1.21.34 (2021-11-12)                 

Changes in existing functions

- Habil renamed identify.modules() to determine.modules().

                Changes in version 1.21.30 (2021-11-12)                 

New functions

- Neda exported identify.modules(), make.filter(), and
  apply.filter() functions.

[plotgardener](/packages/plotgardener)
------------

                       Changes in version 1.1.18                        

NEW FEATURES

- hicTriangles and hicRectangles can now be annotated with
annoDomains
or annoPixels if they are flipped.

                       Changes in version 1.1.17                        

NEW FEATURES

- plotIdeogram can now accept custom colors with a fill parameter.
Colors can be specified with a named or unnamed vector. To see which
stains are being assigned which colors, look inside the ideogram
object.

                       Changes in version 1.1.16                        

BUG FIXES

- ENTREZ IDs obtained from AnnotationDbi::select() are subset just
for
ENTREZID column when determining default gene priorities,
eliminating dplyr incompatible types error.
- All plus and minus strand gene name label parsing in plotGenes is
now carried out only if there is a non-zero number of that strand's
genes.

                       Changes in version 1.1.15                        

- Citation linked for plotgardener publication in Bioinformatics.

                       Changes in version 1.1.14                        

BUG FIXES

- plotSignal yrange parsing for negative scores now has fixed the
typo
on line 418 from "score2" to "score".

                       Changes in version 1.1.13                        

BUG FIXES

- plotSignal default yrange parsing now catches the invalid 0,0 range
and no longer throws a viewport related error.

                       Changes in version 1.1.12                        

BUG FIXES

- readHic and functions related to the reading of .hic files now
leaves the chromosome input formatted as is (e.g. "chr1" and "1").
Functions will throw an error if the input chromosome is not found
in the chromosomes listed in the .hic file.

                       Changes in version 1.1.11                        

BUG FIXES

- annoDomains coordinates fixed for plotHicRectangle.
- Clipping logic for plotPairsArches now clips arches both on left
and
right side of plot.
- Subsetting plotPairs logic fixed to match plotPairsArches.

NEW FEATURES

- clip.noAnchors parameter in plotPairsArches allows for inclusion or
clipping of arches that do not have anchors in the given genomic
region.
- plotPairsArches now allows for column name input to designate
archHeights.

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

                       Changes in version 1.5.16                        

- MSnbase::MSnSet class has been replaced by the SummarizedExperiment
class
- Color scale for all plots set to viridis (without yellow)
- All output tables provided as tibble istead of matrix or data.frame
- Allow users to select specific covariates and their position
(importance) in the model for PomaLimma, PomaUnivariate(method =
"anova"), and PomaOddsRatio
- Add SD to PomaUnivariate output tables
- Reduce dependencies
- Update vignettes
- Update documentation
- Some other major and minor improvements
- Some minor bugs and typos fixed
- New "loading plot" in PCA
- Compute p-values and FDR in PomaCorr

[preciseTAD](/packages/preciseTAD)
----------

                 Changes in version 1.5.0 (2022-04-21)                  

- Match NEWS and DESCRIPTION versioning

[progeny](/packages/progeny)
-------

               Changes in version 2020-10-14 (2020-10-14)               

- Model matrices are not accessed in the local and not in the global
enviroment

               Changes in version 2020-09-01 (2020-09-01)               

- Fixed issue with rownames when using Progeny with Permutations
function

               Changes in version 2020-06-09 (2020-06-09)               

- Website: Google Analytics

               Changes in version 2020-04-27 (2020-04-27)               

- PROGENy website development

Major update with the following main points:

- Added the mouse model matrix containing 14 pathways

- The human model matrix extended to 14 pathways

- Added the following functions: progenyPerm, progenyScatter,
progenySavePlots, getModel

- Added tests and test data

- Added the vignette for usage the PROGENy on single-cell RNA-seq
data

- Added functionality to work with Seurat objects

[pRolocGUI](/packages/pRolocGUI)
---------

                         Changes in version 2.5                         

CHANGES IN VERSION 2.5.4

- Fixed issue #108 by removing scrollX known issue in DT package

CHANGES IN VERSION 2.5.3

- Update README

CHANGES IN VERSION 2.5.2

CHANGES IN VERSION 2.5.1

- fix #109 bug in pRolocVis when 'markers' is missing from fvarLabels
- fontawesome deprecated switch to cogs for sidebar
- fix bug in colour menu #110 attributed from shinyWidgets changing
the class name from checkboxGroupButtons to checkbox-group-buttons

CHANGES IN VERSION 2.5.0

- New version for Bioc

[ProteoDisco](/packages/ProteoDisco)
-----------

                        Changes in version 1.1.3                        

- Added citation to manuscript.

[protGear](/packages/protGear)
--------

                Changes in version 0.99.546 (2022-03-25)                

- Updated the shiny app path
- Added Bioconductor installation in the app

                Changes in version 0.99.545 (2022-03-23)                

- Updated the NEWS file
- Removed other non Bioconductor files

                Changes in version 0.99.544 (2022-03-17)                

- Removed non Bioconductor files
- Added Suggests in Description
- Updated R requirement to 4.2
- Added the data documentation

                Changes in version 0.99.543 (2022-02-14)                

- Submitted to Bioconductor for review

                Changes in version 0.99.55 (2022-04-08)                 

- removed vignette.R and csv file

                 Changes in version 0.99.1 (2021-12-15)                 

- Submitted to Bioconductor

[ProtGenerics](/packages/ProtGenerics)
------------

                       Changes in version 1.27.1                        

- Add addProcessing generic <2022-01-04 Tue>

- Add new adjacencyMatrix generic <2021-12-11 Sat>

                       Changes in version 1.27.0                        

- New Bioc devel version

[psichomics](/packages/psichomics)
----------

                       Changes in version 1.20.2                        

- Bug fix: fix limma-trend approach when using newer versions of
limma
to calculate average gene expression

                       Changes in version 1.20.1                        

- Bug fix: allow to perform correlation analysis after being
performed
once

[PSMatch](/packages/PSMatch)
-------

                        Changes in version 0.99                         

Changes in 0.99.5

- Fix mz calculation in calculateFragments for neutral losses with a
charge > 1 (ported from MSnbase - see issue 573).

Changes in 0.99.4

- Set seed in the ConnectedComponents unit test to stop random errors
after clustering.

Changes in 0.99.3

- Fix bug in describePeptides() (close #11).

Changes in 0.99.2

- Describe the ConnectedComponents() return value.
- Add/update installation instructions.

Changes in 0.99.1

- Fix typo and improve documentation.

Changes in 0.99.0

- Prepare package for Bioconductor submission.

[PureCN](/packages/PureCN)
------

                        Changes in version 2.2.0                        

NEW FEATURES

- Added chunks parameter to Coverage.R and
  calculateBamCoverageByInterval
  to reduce memory usage (#218)

SIGNIFICANT USER-VISIBLE CHANGES

- When base quality scores are found in the VCF, they are now used to
  calculate the minimum number of supporting reads (instead of assuming
  a
  default BQ of 30). By default BQ is capped at 50 and variants below
  25
  are ignored. Set min.supporting.reads to 0 to turn this off (#206).

- More robust annotation of intervals with gene symbols

- Remove chromosomes not present in the centromeres GRanges object;
  useful
  to remove altcontigs somehow present (should not happen with
  intervals
  generated by IntervalFile.R)

BUGFIXES

- Fixed an issue with old R versions where factors were not converted
  to
  strings, resulting in numbers instead of gene symbols

- Fix for a crash when there are no off-target reads in off-target
  regions
  (#209).

- Fixed parsing of base quality scores in Mutect 2.2

- Fixed crash in GenomicsDB parsing when there were no variants in
  contig
  (#225)

[QDNAseq](/packages/QDNAseq)
-------

                 Changes in version 1.31.0 (2021-10-27)                 

RELEASE

- The version number was bumped for the Bioconductor release version,
  which is
  now BioC 3.15 for R-devel.

[QFeatures](/packages/QFeatures)
---------

                         Changes in version 1.5                         

QFeatures 1.5.2

- fix: implemented an updateObject() method for QFeatures objects.

QFeatures 1.5.1

- Document the use of peptide/protein adjacency matrices in
aggregateFeatures() and new adjacencyMatrix() accessor.

QFeatures 1.5.0

- New devel version (Bioc 3.15)

[qmtools](/packages/qmtools)
-------

                 Changes in version 0.99.3 (2022-04-12)                 

- The package has been accepted

- Updated the installation instruction in the README.md

                 Changes in version 0.99.2 (2022-04-11)                 

- Updated the package documents to respond to the 2nd round of
  Bioconductor
  review

                 Changes in version 0.99.1 (2022-03-17)                 

- Renamed the package as "qmtools"

- Made the significant changes overall according to the Bioconductor
  review

- Added the `removeFeatures` function to filter uninformative features
  from the
  data

- Added the `clusterFeatures` function to identify a group of features
  from the
  same originating compound

                 Changes in version 0.99.0 (2022-01-03)                 

- Submitted to Bioconductor (The package was previously named as
  "poplin")

[qpgraph](/packages/qpgraph)
-------

                        Changes in version 2.30                         

BUG FIXES

- Fixed NAMESPACE issues

- Fixed calls to some Fortran LAPACK functions to comply with Writing R
  Extensions §6.6.1

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

                       Changes in version 1.13.2                        

- Fixed guide = FALSE to guide = "none".

                       Changes in version 1.13.1                        

- Fixed bug in computeDiffStats.

[qsea](/packages/qsea)
----

                       Changes in version 1.21.2                        

- increased upper bound for enrichment model, allowing for a steeper
  increase in enrichment with CpG density. (Github PR \#9)

                       Changes in version 1.21.1                        

- Bugfixes:
  - makeTable with one window (github PR \#8)

[qsvaR](/packages/qsvaR)
-----

                       Changes in version 0.99.0                        

NEW FEATURES

- This is the initial version of the second iteration of the qSVA
framework, which was initially described by Jaffe et al, PNAS, 2017
https://doi.org/10.1073/pnas.1617384114.

[Qtlizer](/packages/Qtlizer)
-------

                 Changes in version 1.8.1 (2022-02-28)                  

- Citation adapted

[RAREsim](/packages/RAREsim)
-----

                  Changes in version 0.99.4

- Fixed format of afs_afr and nvariant_afr data

                  Changes in version 0.99.3

- Fixed "Installing the Package" in the vignette

                  Changes in version 0.99.2

- Reformatted vignette
- Addressed BiocCheck notes

                  Changes in version 0.99.0

- Submitted to Bioconductor

[rawrr](/packages/rawrr)
-----

                  Changes in version 1.3 (2022-03-19)                   

- Add barebone mode para in readSpectrum #43.

- Add rawrr namespace in help pages.

- Add 'Monoisotopic M/Z:' from TrailerExtraHeaderInformation as
  column to rawrr::readIndex function.

[Rbwa](/packages/Rbwa)
----

                        Changes in version 1.0.0                        

- New package Rbwa, R wrapper for BWA aligner

[RcisTarget](/packages/RcisTarget)
----------

                        Changes in version 1.15                         

- New function: showLogo() shows the motif enrichment table as HTML.

- Fix for maxRank checks: Now takes into account number of
  genes/regions in the database.

[RCM](/packages/RCM)
---

                       Changes in version 1.11.3                        

- For the unconstrained models: fit feature models one by one and
Gram-Schmidt orthogonalize and center afterwards, rather than using
Lagrange multipliers and huge Jacobian matrices. This will use less
memory and speed up computations, but may yield slightly different
solutions. Nothing changes for the constrained models.

                       Changes in version 1.11.2                        

- Explicitly import stats::model.matrix, and only load necessary VGAM
functions

                       Changes in version 1.11.0                        

- Added FAQ section in vignette with first frequent question on
number
of samples not shown.
- Fixed bugs for plots of data with missing values, and added
tests.

[Rcpi](/packages/Rcpi)
----

                 Changes in version 1.31.1 (2022-04-25)                 

Improvements

- Fixed a build error on macOS in the devel branch due to
dependencies
not available.

[RCy3](/packages/RCy3)
----

                       Changes in version 2.16.0                        

- Faster selectAll* functions

- Add a new vignette about cloud notebooks with RCy3

- Doc fixes:
  - Conflicting Brightness/Contrast documentation, #172
  - Conflicting Opacity documentation, #173
  - updateAnnotationText cleanup, #177

- New functions:
  - createView
  - selectAll

- Bug fixes:
  - addAnnotationShape customShape can only add rectangle, #160
  - setEdgeLineWidthMapping issue, #164
  - openAppStore function opens the Cytoscape App Store 404 web page,
  #169
  - groupAnnotation cleanup, #175

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.9.1 (2021-11-05)                  

- Adapted default FDR threshold in "plot_heatmap" to match other
  functions.

[recountmethylation](/packages/recountmethylation)
------------------

                        Changes in version 1.5.1                        

- Adds QC functions for BeadArray metrics and log M/U signals

- Adds data and accessor function for cross-reactive CpGs

- Adds vignette showing how to do power analysis with pwrEWAS

- Adds vignette showing how to infer genetic ancestry using
  GLINT/EPISTRUCTURE

- Adds vignette showing how to do nearest neighbors search using
  a search index

- Adds functions for feature hashing, search index construction,
  and KNN search

[RedeR](/packages/RedeR)
-----

                        Changes in version 2.0.0                        

- Major upgrade of the user interface.

[rgoslin](/packages/rgoslin)
-------

                        Changes in version 0.99                         

Please note that this Bioconductor version is based on Goslin
version 2.0.0. See the Goslin repository for more details.

Changes in 0.99.1

- The column names within the data frames returned from the parse*
methods now use column names with dots instead of spaces. This makes
it easier to use the column names unquoted within other R
expressions.
- All parse* methods now return data frames.
- The Messages column has been added to capture parser messages. If
parsing succeeds, this will contain NA and Normalized.Name will
contain the normalized lipid shorthand name.
- Parser implementations have been updated to reflect the latest
lipid
shorthand nomenclature changes. Please see the Goslin repository for
more details.
- Exceptions in the C++ part of the library are captured as warnings
in R. However, if you parse multiple lipid names, exceptions will
not stop the parsing process.

[rgsepd](/packages/rgsepd)
------

                        Changes in version 1.27                         

- library(hash) is going away in 2022, so we need to replace that
  functionality with native R 4.2 environments. Shouldn't impact users.

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.40.0                        

NEW FEATURES

- Added H5R functions for working with object and dataset region
  references.

- The HDF5 N-Bit filter has been enabled with via the function
  H5Pset_nbit().  This can be combined with H5Tset_precision() to
  compress integer and floating-point datasets.

CHANGES

- Argument 'cset' to h5createAttribute() and h5writeAttribute()
  have been deprecated.  The 'encoding' argument should be used
  going forward.  This ensures consistency with h5create() and
  h5write().

BUG FIXES

- The documentation for the 'encoding' argument to
  h5createDataset() and h5writeDataset() stated 'UTF-8' was a
  valid option, however this would produce an error. This has now
  been fixed. (Thanks to @ilia-kats for identifying this,
  https://github.com/grimbough/rhdf5/pull/101)

- Fixed bug where an uninitialized value was used in the C code
  underlying h5dump() potentially causing crashes.

- Addressed issue in h5dump() and h5ls() that falsely declared
  there were duplicated groups when used on a file with external
  links (Thanks to @acope3 for reporting this,
  https://github.com/grimbough/rhdf5/issues/107).

[Rhdf5lib](/packages/Rhdf5lib)
--------

                        Changes in version 1.18                         

New features

- Package now includes precompiled libraries for Windows built with the
  UCRT toolchain for R-4.2

- Swap bundled version of SZIP for LIBAEC.  This now reflects the
  official
  HDF5 group releases.

- The HDF5 configure option "--disable-sharedlib-rpath" is now exposed
  during
  package installation
  (thanks to Ben Fulton @benfulton,
  https://github.com/grimbough/Rhdf5lib/pull/39)

[Rhtslib](/packages/Rhtslib)
-------

                       Changes in version 1.28.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Acknowledge site-wide or user-specified Makevars file if present.

[RiboDiPA](/packages/RiboDiPA)
--------

                        Changes in version 1.3.1                        

- Functions to visualize tracks through genome browser igvr added

- Functions to export ribo-seq tracks to external genome browser added

[ribosomeProfilingQC](/packages/ribosomeProfilingQC)
-------------------

                        Changes in version 1.7.1                        

- export readLen for estimatePsite.

[RMassBank](/packages/RMassBank)
---------

                        Changes in version 3.5.1                        

- Switch to using mzML files in the vignette after mzR dropped support
  for mzData

[RnaSeqSampleSize](/packages/RnaSeqSampleSize)
----------------

                 Changes in version 2.5.1 (2021-12-02)                  

NEW FEATURES

- Support defining ratio of sample size between two groups by parameter
  k

BUG FIXES

- Change website URL in vignette

- Not-Run some example codes

[RnBeads](/packages/RnBeads)
-------

                       Changes in version 2.13.2                        

- Removed support for RefFreeEWAS, since the package is not supported
  anymore

[RolDE](/packages/RolDE)
-----

                       Changes in version 0.99.4                        

- Updated version dependency to R to 4.2.0.

                       Changes in version 0.99.3                        

- Added Bioconductor installation instructions in the Vignette.

- Fixed some bugs related to changing coding practices in v. 0.99.2

                       Changes in version 0.99.2                        

- Added a NEWS file.

- Added Bioconductor installation instructions in README.

- Removed separate licence file. Using GPL-3 licence.

- Added information for the included datasets.

- Added table of contents for the vignette.

- Updated the RolDE main functions documentation.

- Improved coding practices to match more Bioconductor style.

                       Changes in version 0.99.1                        

- Submitted to Bioconductor.

[rols](/packages/rols)
----

                        Changes in version 2.23                         

CHANGES IN VERSION 2.23.2

- Remove failing CVParams() examples

CHANGES IN VERSION 2.23.1

- Don't show empty termDesc(trm) in vignette.

CHANGES IN VERSION 2.23.0

- New devel version

[ropls](/packages/ropls)
-----

                       Changes in version 1.27.8                        

MINOR MODIFICATION

- minor vignette correction

                       Changes in version 1.27.6                        

MINOR MODIFICATION

- minor vignette correction

                       Changes in version 1.27.4                        

MINOR MODIFICATION

- minor vignette update

                       Changes in version 1.27.2                        

MINOR MODIFICATION

- minor documentation update

[rpx](/packages/rpx)
---

                        Changes in version 2.3.0                        

rpx 2.3.3

- Don't run caching example as it replies on PRIDE which is failing
too often for the tests to pass on all systems.

rpx 2.3.2

- Provide fix for PRIDE migration and annotation discrepancies (see
issue #17). This fix will however lead to re-downloading some cached
files due to different URLs. The fix might only be temporary, based
on if/how PRIDE will address their inconsistencies.

rpx 2.3.1

- Use BiocFileCache::bfcrpath to get correct rpath, irrespective of
absolute or relative.

rpx 2.3.0

- New bioconductor devel branch

[rrvgo](/packages/rrvgo)
-----

                        Changes in version 1.7.1                        

- Deprecate org.Pf.plasmo.db

[Rsubread](/packages/Rsubread)
--------

                       Changes in version 2.10.0                        

- 
  Added inbuilt RefSeq annotation for mm39 (mouse genome Build
  39).

- 
  Streamlined the mapping and counting processes in cellCounts.

- 
  Added support for processing dual-index 10x data in cellCounts.

[rWikiPathways](/packages/rWikiPathways)
-------------

                       Changes in version 1.16.0                        

- HTTPS updates

- Update vignette to include new clusterProfiler functions

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.34.0                        

NEW FEATURES

- Implement subassignment of TransposedDataFrame objects.

SIGNIFICANT USER-VISIBLE CHANGES

- DataFrame is now a virtual class! This completes the replacement of
  DataFrame with DFrame announced in September 2019. See:
  https://www.bioconductor.org/help/course-materials/2019/BiocDevelForum/02-DataFrame.pdf

BUG FIXES

- Avoid spurious warnings when DataFrame() is supplied a DataFrameList.

- Fix bug in combineRows(DataFrame(), DataFrame(ref=IRanges(1:2, 10))).

- Fix display of TransposedDataFrame objects with more than 11 rows.

- Fix isEmpty() on ordinary lists and derivatives.

- Fix handling of nested DataFrames in combineUniqueCols().

- Make sure internal helper lowestListElementClass() does not lose
  the "package" attribute of the returned class (fixes issue #103).

[satuRn](/packages/satuRn)
------

                        Changes in version 1.3.1                        

- Bug fix: allow for fit errors to be propagated as NA results
(github
issue 15 by @jgilis)
- Bug fix: handle experimental designs with empty factor levels
correctly (github issue 16 by @XueyiDong)
- Bug fix: identify transcripts that are the only expressed
transcript
of a gene and set NA results (github issue 17 by @jgilis)
- Bug fix: handle extreme z-scores in testDTU with diagplot2 option
- Enhancement: plotDTU now allows for sparseMatrix input

[SCArray](/packages/SCArray)
-------

                        Changes in version 1.4.0                        

- new functions `scMEX2GDS()` and `scHDF2GDS()`

                        Changes in version 1.2.1                        

- new Overview.Rmd in the vignettes

[scater](/packages/scater)
------

                       Changes in version 1.24.0                        

- Remove diffusion map functions that relied on destiny.

- Add point.padding,force args to plotReducedDim; passed to
  geom_text_repel.

- Add warning about unused use_dimred argument in runTSNE.

[scDblFinder](/packages/scDblFinder)
-----------

                 Changes in version 1.9.11 (2022-04-16)                 

- fixed larger kNN size

                        Changes in version 1.9.9                        

- improved amulet reimplementation

- added clamulet and scATAC vignette

                 Changes in version 1.9.1 (2021-11-02)                  

- added reimplementation of the amulet method for scATAC-seq

[scPCA](/packages/scPCA)
-----

                 Changes in version 1.9.1 (2022-01-19)                  

- Updating copyright years
- Updating citation information

[scTensor](/packages/scTensor)
--------

                        Changes in version 2.4.1                        

- A vignette modified.

[sechm](/packages/sechm)
-----

                 Changes in version 1.3.1 (2022-04-15)                  

- fixed a few bugs

- better support for default arguments

- removed deprecated arguments

- added meltSE

[seqArchR](/packages/seqArchR)
--------

                      Changes in version 0.99.339                       

- seqArchR available on Bioconductor

                       Changes in version 0.99.0                        

New features

- viz_seqs_acgt_mat() function can now add a legend via new arguments
add_legend and use_legend.
- Package name change from archR to seqArchR

Breaking changes

- (User-facing) Function name archR_set_config() changed to
set_config()
- (User-facing) Function name viz_seqs_acgt_mat_from_seqs() changed
to
viz_seqs_acgt_mat()
- (User-facing) Function runArchRUI() that launched a Shiny app is
moved to a different package coming up in the future

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.36.0                        

NEW FEATURES

- new functions `seqUnitCreate()`, `seqUnitSubset()` and
  `seqUnitMerge()`

- new functions `seqFilterPush()` and `seqFilterPop()`

- new functions `seqGet2bGeno()` and `seqGetAF_AC_Missing()`

- new function `seqGetData(, "$dosage_sp")` for a sparse matrix of
  dosages

- the first argument 'gdsfile' can be a file name in `seqAlleleFreq()`,
  `seqAlleleCount()`, `seqMissing()`

- new function `seqMulticoreSetup()` for setting a multicore cluster
  according to a numeric value assigned to the argument 'parallel'

UTILITIES

- allow opening a duplicated GDS file ('allow.duplicate=TRUE') when the
  input is a file name instead of a GDS object in `seqGDS2VCF()`,
  `seqGDS2SNP()`, `seqGDS2BED()`, `seqVCF2GDS()`, `seqSummary()`,
  `seqCheck()` and `seqMerge()`

- remove the deprecated '.progress' in `seqMissing()`,
  `seqAlleleCount()`
  and `seqAlleleFreq()`

- add `summary.SeqUnitListClass()`

- no genotype and phase data nodes from `seqSNP2GDS()` if SNP dosage
  GDS
  is the input

BUG FIXES

- `seqUnitApply()` works correctly with selected samples if 'parallel'
  is
  a non-fork cluster

- `seqVCF2GDS()` and `seqVCF_Header()` work correctly if the VCF header
  has
  white space

- `seqGDS2BED()` with selected samples for sex and phenotype
  information

- buf fix in `seqGDS2VCF()` if there is no integer genotype

[seqcombo](/packages/seqcombo)
--------

                       Changes in version 1.17.1                        

- update docs (2021-12-15, Wed)
- remove codes that were incorporated in ggmsa

[SEtools](/packages/SEtools)
-------

                 Changes in version 1.9.2 (2022-04-16)                  

- removed functions that have been moved to the sechm package

[signatureSearch](/packages/signatureSearch)
---------------

                 Changes in version 1.9.6 (2022-01-28)                  

- Addressed warning for non-ascii in cell_info2

                 Changes in version 1.9.5 (2022-01-26)                  

- Updated lincs_pert_info2 and cell_info2

                 Changes in version 1.9.3 (2021-12-17)                  

- Supported searching against the newest LINCS 2020 beta database in
devel version
- Modified gess_* functions to support adding customized compound
annotation table to the GESS result table.

                 Changes in version 1.9.2 (2021-12-06)                  

- Move eh to .onLoad function

[simplifyEnrichment](/packages/simplifyEnrichment)
------------------

                        Changes in version 1.5.2                        

- add `keyword_enrichment_from_GO()`

                        Changes in version 1.5.1                        

- word cloud supports perform enrichment on keywords

- value_fun in binary_cut() now takes 1-AUC as default

[singleCellTK](/packages/singleCellTK)
------------

                        Changes in version 2.5.2                        

- Added Seurat report functions
- Added TSCAN trajectory analysis functions
- Refactored EnrichR wrapper function (runEnrichR)
- Added new cut-offs for DE functions
- Other refactors and bug fixes

                 Changes in version 2.5.1 (2022-03-31)                  

- Added SoupX method for decontamination (runSoupX)
- Added useReducedDim parameter for DE analysis and Heatmap
- Added Differential Abundance section to the tutorials
- Fixed Mitochondrial gene list
- Other refactors and bug fixes

                 Changes in version 2.4.1 (2021-12-22)                  

- Added new function for DEG volcano plot (plotDEGVolcano)
- Added new function for plotting pathway scores (plotPathway)
- Added Pathway Analysis section to the tutorials
- Added seed parameter to several functions and UI for
reproducibility
- Updated R console and GUI tutorials to match each other
- Fixed console logging in the GUI

[sitePath](/packages/sitePath)
--------

                       Changes in version 1.11.2                        

- More function availability for objects

                       Changes in version 1.11.1                        

- Wrapper function for finding fixation and parallel sites

- Core number set to 1 will disable multiprocessing

[SNPRelate](/packages/SNPRelate)
---------

                       Changes in version 1.30.0                        

- return a object of S3 class "snpgdsGRMClass" in `snpgdsGRM()` instead
  of
  a list when `with.id=TRUE`

[sparrow](/packages/sparrow)
-------

                         Changes in version 1.2                         

Enhancements

- The default "zero-centering" logic is updated in mgheatmap2 when
col
isn't specified, but recenter is (backported to release 3.14)

Bug Fixes

- calculateIndividualLogFC is updated to handle situations when
$genes
data.frame has column names that collide with statistics generated
from differential expression, like pval, padg, and AveExpr. Thanks
to @sandersen12 for the bug report.

[SpatialDecon](/packages/SpatialDecon)
------------

                 Changes in version 1.5.0 (2021-10-27)                  

- No changes from 1.3.0

[SpatialExperiment](/packages/SpatialExperiment)
-----------------

                 Changes in version 1.5.3 (2022-02-28)                  

- rename SpatialImage class to VirtualSpatialImage

                 Changes in version 1.5.2 (2022-01-09)                  

- relocate and deprecate spatialData/Names

                 Changes in version 1.5.1 (2021-12-15)                  

- improved coercion methods from SingleCellExperiment to
  SpatialExperiment

- add new methods for image rotation/mirroring

- add path argument to imgSource()

- user flexibility whether to provide outs/ directory in
  read10xVisium()

- documentation updates in show methods

- additional documentation updates

- update title and description in DESCRIPTION

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 2.1.1 (2022-04-06)                  

- Developed the new functionality co-visualization of bulk and single
  cell data through auto-matching/coclustering, i.e. source bulk
  tissues are matched/assigned to single cells automatically through
  coclustering. This feature is implemented in both command line and
  Shiny app with testing data provided.

- Developed optimization functions for coclustering workflow with
  testing data provided.

- Co-visualization through manual matching was implemented in command
  line.

[spatzie](/packages/spatzie)
-------

                        Changes in version 1.0.1                        

- added reference to NAR paper

[Spectra](/packages/Spectra)
-------

                         Changes in version 1.5                         

Changes in 1.5.20

- Add parameters ppm and tolerance to PrecursorMzParam (for neutral
loss calculation) and add option filterPeaks = "removePrecursor".

Changes in 1.5.19

- Improved the bin method.

Changes in 1.5.18

- Set default for parameter columns in peaksData,Spectra and
peaksData,MsBackend to c("mz", "intensity").

Changes in 1.5.17

- Add peaksVariables method and add parameter columns (or ...) to
peaksData.
- Add columns parameter to the peaksData method of
MsBackendDataFrame,
MsBackendMzR and MsBackendHdf5peaks.

Changes in 1.5.16

- Fix issue in neutralLoss that would prevent calculation of neutral
loss spectra if

Changes in 1.5.15

- Fix typo in MZ delta plot title.

Changes in 1.5.14

- Add coreSpectraVariables function to export the core spectra
variables and their expected data types.

Changes in 1.5.13

- Fix figure sizes in vignette.

Changes in 1.5.12

- Add neutralLoss method and first algorithm to calculate neutral
loss
spectra.

Changes in 1.5.11

- Fix neutral loss example in the vignette.

Changes in 1.5.10

- Add citation.

Changes in 1.5.9

- Add examples for combineSpectra to the vignette.

Changes in 1.5.8

- Add spectraVariableMapping generic.

Changes in 1.5.7

- Add missing export of the filterPrecursorMz method.

Changes in 1.5.6

- Add filterPrecursorMzValue method which allows to filter using
multiple precursor m/z values (issue #230).
- Fix unit test suite.

Changes in 1.5.5

- Add a testing framework allowing to run standardized unit tests for
new MsBackend implementations (issue #186).

Changes in 1.5.4

- Add the MsBackendCached backend.

Changes in 1.5.3

- Only calculate number of peaks per spectra if the processing queue
of the Spectra is not empty. Otherwise call the backend's
implementation (issue MsBackendSql #31).

Changes in 1.5.2

- Small documentation update (related to MsCoreUtils issue #87).
- New countIdentifications() function.
- Add filterFourierTransformArtefacts function to remove fast fourier
artefact peaks seen on e.g. Orbitrap instruments (issue #223).

Changes in 1.5.1

- Don't read header information when importing peaks matrix on macOS.

[SpectralTAD](/packages/SpectralTAD)
-----------

                 Changes in version 1.11.0 (2022-04-21)                 

- Match NEWS and DESCRIPTION versioning

[splatter](/packages/splatter)
--------

                 Changes in version 1.20.0 (2022-04-27)                 

- 
  The splatPop simulation is now published
  doi.org/10.1186/s13059-021-02546-1!

- Improved initalisation of Params objects (from Wenjie Wang)

- 
  Improved fitting of dropout in splatEstimate()
  
  • Better initialisation of fitting as suggested by the
  InferCNV package
  
  • Additional fallback method

- Bug fixes for the splat simulation

- 
  Bug fixes for the the splatPop simulation (from Christina
  Azodi)

[SplicingFactory](/packages/SplicingFactory)
---------------

                        Changes in version 1.3.1                        

- Citation update
- Other minor corrections

[SPOTlight](/packages/SPOTlight)
---------

                       Changes in version 0.99.1                        

- text

                       Changes in version 0.99.0                        

- initial submission to Bioc devel v3.15

[sSNAPPY](/packages/sSNAPPY)
------

                       Changes in version 0.99.1

- Updated vignette to use pre-computed output
- Allow removal of isolated nodes in `plot_gs_network`

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[standR](/packages/standR)
------

                       Changes in version 0.99.0                        

- First release of the package.

[statTarget](/packages/statTarget)
----------

                       Changes in version 1.25.4                        

- Remove the RGtk2 in Description, which should be installed manually
  from the old source

                       Changes in version 1.24.1                        

- Remove the GUI

[STdeconvolve](/packages/STdeconvolve)
------

                        Changes in version 0.99.10

- R version for Bioconductor and `package` branch: (>= 4.1)
- R version in `devel` branch in GitHub: (>= 3.6)

                        Changes in version 0.99.0

- Submitted to Bioconductor

[struct](/packages/struct)
------

                        Changes in version 1.7.1                        

- fix as.code generic

- improve as.code for all objects

[structToolbox](/packages/structToolbox)
-------------

                        Changes in version 1.7.2                        

- internal updates to PLS

- some PLS charts have been renamed for consistency with other methods

                        Changes in version 1.6.1                        

- hotfix vector_norm now correctly normalises samples to length 1

- added vector_norm_tests

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.26.0                        

DEPRECATED AND DEFUNCT

- readKallisto() is now defunct after being deprecated in BioC 3.12.

[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.7.14                        

- Various improvements for GenRearrScen, improves consistency and
output formatting
- Major bugfix for ProtWeaver methods using dendrogram objects
- ProtWeaver now correctly guards against non-bifurcating dendrograms
in methods that expect it

                       Changes in version 1.7.13                        

- Introduces new ProtWeaver class to predict functional association
of
genes from COGs or gene trees. This implements many algorithms
commonly used in the literature, such as MirrorTree and Inverse
Potts Models.
- predict(ProtWeaverObject) returns a ProtWeb class with information
on predicted associations.
- Adds BlastSeqs to run BLAST queries on sequences stored as an
XStringSet or FASTA file.

                       Changes in version 1.7.12                        

- Updates to ExtractBy function. Methods and inputs simplified and
adjusted, and significant improvements to speed.

                       Changes in version 1.7.11                        

- Updated NucleotideOverlaps to now correctly registers hits in genes
with a large degree of overlap with the immediately preceding gene.
- Fixed aberrant behavior in BlockExpansion where contigs with zero
features could cause an error in expansion attempts.

                       Changes in version 1.7.10                        

- BlockReconciliation now allows for setting either block size or
mean
PID for reconciliation precedence.

                        Changes in version 1.7.9                        

- Added retention thresholds to BlockReconciliation.

                        Changes in version 1.7.8                        

- BlockExpansion cases corrected for zero added rows.

                        Changes in version 1.7.7                        

- Improvements to BlockExpansion and BlockReconciliation functions.

                        Changes in version 1.7.5                        

- Began integration of DECIPHER's ScoreAlignment function.

                        Changes in version 1.7.4                        

- Fixed a bug in PairSummaries function.

                        Changes in version 1.7.3                        

- Added BlockExpansion function.

                        Changes in version 1.7.2                        

- Adjustment in how PairSummaries handles default translation tables
and GFF derived gene calls.

                        Changes in version 1.7.1                        

- Large changes under the hood to PairSummaries.
- Failure to accurately assign neighbors in some cases should now be
fixed.
- Extraction of genomic features is now faster.
- OffSetsAllowed argument now defaults to FALSE. This argument may be
dropped in the future in favor of a more complex function
post-summary.
- Small edits to SequenceSimilarity

[systemPipeShiny](/packages/systemPipeShiny)
---------------

                       Changes in version 1.5.10                        

Major Change

- Redesign of the welcome page. Old content is moved to about. Now
the
welcome page is more clear.

- Adapt SPR to the 2.1.x version.

- Add more instructing images to the workflow module.
- A warning message is added if spsOption("demo", TRUE), to let
people know some workflow templates will fail if they use the
demo server to run jobs.

Minor Change

- Fix some text typo, links.
- Add more figures as instructions in different modules/tabs.
- Text/links fixed in workflow module.

[TADCompare](/packages/TADCompare)
----------

                 Changes in version 1.5.0 (2022-04-21)                  

- Match NEWS and DESCRIPTION versioning

[TargetSearch](/packages/TargetSearch)
------------

                       Changes in version 1.52.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Deprecate parameter 'libId' and replace it with 'libID'
  The parameter `libId` conflicts with the method libId, and in other
  functions we already use a variable `libID` that correspond with the
  library identifier. Therefore, we replace it for consistency.
  This change affects the functions `plotRIdev`, `plotSpectra` and

BUG FIXES

- plotRIdev: use near equality for comparison.
  The function compares the selective or top masses against the columns
  of the intensity matrix. However, their lengths are variable, so the
  for `if` condition threw a warning if the lengths didn't match.

- Refactor functions binsearch and find_peaks.
  Minor optimization to the function binsearch in which the starting
  RT scan is found directly. This will reduce a couple of CPU cycles.

- Extra CDF integrity checks. Check that the length of the variables
  are all greater than zero and replace logical OR operator for its
  longer form.

                       Changes in version 1.50.1                        

BUG FIXES

- Fixes the gcc warning produced by passing a pointer and with an
  incorrect
  size to the function `swapb`. This went undiscovered for years
  because
  the code is only executed in big-endian machines.

[TBSignatureProfiler](/packages/TBSignatureProfiler)
-------------------

                        Changes in version 1.7.1                        

Bug Fixes

- Fixed gene in RESPONSE5 (PNN to RP11-295G20.2) in TBsignatures and
TBcommon objects. (Stanley M. Kimbung)

Major Changes

- If any signatures in the object used with runTBsigProfiler() have
<2
genes present in the given sample, the signatures will not be
scored. This may affect existing scripts.
- Fixed a bug in the singscore algorithm called from
runTBsigProfiler() that would not allow for the scoring of
user-provided signatures.

Minor Changes

- Reorganized code for OriginalModel.R for clarity.
- Fixed the TB_hiv data to remove unnecessary factor level of Disease
metadata.
- Added 4 new signatures (Tabone_OD_11/TB12,
Tabone_RES_25/EarlyRESP-TB25, Tabone_RES_27/TREAT-TB27, Long_RES_10)
- Updated the website interface
- Changed HGNChelper installation to be checked during profiling if
update_genes = TRUE
- Reorganized code in mkAssay() for clarity. The output_name argument
is now appended to all output assays, whereas previously it was only
appended to the log of the input assay.
- Fixed the row numbers of existing sigAnnotData and
common_sigAnnotData objects, and added code to update them after new
signatures are added.

[TCGAutils](/packages/TCGAutils)
---------

                       Changes in version 1.16.0                        

New features

- The UUIDhistory function allows users to map old UUIDs to new UUIDs
according to the latest data release for UUIDs that were affected
and no longer query-able.
- The slides argument has been added to the filenameToBarcode
function
for translating slide file names into barcodes. Currently, the API
returns all barcodes of the associated case ID.
- Add sections in the vignette regarding GDC Data Updates and UUID
history lookup

Minor changes and bug fixes

- Update examples in package to new GDC Data Release, see vignette.
- Use AnnotationHub to download chain file in main vignette.
- Slide file names now resolve to a single TCGA barcode in
filenameToBarcode (Thanks @hermidalc)
- Improved error messages and documentation for
makeGRangesListFromExonFiles

[TFEA.ChIP](/packages/TFEA.ChIP)
---------

                       Changes in version 1.15.2                        

New ChIP-Gene databases available

- Using ReMap2022 collections for human and mouse.

- Adding cell specific regulatory regions predicted with
  ABC-Enhancer-Gene-Prediction (doi:10.1038/s41588-019-0538-0).

New Features

- New database format (older databases are still compatible). The
  format
  consist of a list containing two elements:
  - Gene Keys: vector of gene IDs
  - ChIP Targets: list of vectors, one per ChIP-seq experiment,
  containing the putative targets assigned. Each target is coded as
  its position in the vector 'Gene Keys'.

- Database generation has been streamlined by joining together the
  functions GR2tfbs_db() and makeTFBSmatrix() into one,
  makeChIPGeneDB().

New default database

- The TF-Gene database included with TFEA.ChIP was built using ReMap's
  ChIP-seq collection (v. 2022) and GeneHancer's Double Elite
  regulatory
  regions (v. 4.8). Because of memory limits, the internal database
  included
  in TFEA.ChIP can only store a fraction of the 8000+ ChIP-seq
  experiments in
  the colection. We selected the 926 ChIP-seq experiments done in
  ENCODE
  Project's Common Cell Types.
  To download the full database, as well as other ready-to-use
  databases
  generated for TFEA.ChIP, visit:
  https://github.com/LauraPS1/TFEA.ChIP_downloads

[tomoda](/packages/tomoda)
------

                        Changes in version 1.5.0                        

- Fixed bug in normalizeTomo()

[tomoseqr](/packages/tomoseqr)
--------

                 Changes in version 0.99.0 (2022-03-28)                 

- Submitted to Bioconductor

[topdownr](/packages/topdownr)
--------

                        Changes in version 1.17                         

Changes in version 1.17.3

- Depending on mzR 2.27.5.
- Fix unit test for .readSpectrum to adapt to new mzR 2.27.5.

Changes in version 1.17.2

- Fix roxygen2 warnings.

Changes in version 1.17.1

- New version for Bioc 3.15 (devel)
- Adapt to new DFrame.

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.31.4                        

- Fix the issue in windows 2022.

                       Changes in version 1.31.2                        

- Handle the issue with long tail.

                       Changes in version 1.31.1                        

- Move the heatmap legend to yaxis.

[transformGamPoi](/packages/transformGamPoi)
---------------

                         Changes in version 1.1                         

- Implement faster scaling with size factors for acosh and log-based
transformations
- Implement analytic Pearson residuals

[transite](/packages/transite)
--------

                       Changes in version 1.12.1                        

- Fix typos

[treeio](/packages/treeio)
------

                       Changes in version 1.19.2                        

- update offspring() to work as child(). Actually they are using the
same function with different default (child(type = "children") and
offspring(type="all")) (2022-03-16, Wed)
- update child() to support different types ("children", 'tips',
'internal', 'external', 'all') (2022-03-09, Wed, #75)
- write.beast allows treedata object only contains phylo slot, then
it
will equivalent to write.nexus (2022-02-23, Wed)

                       Changes in version 1.19.1                        

- bug fixed in groupClade.treedata to return a treedata object
instead
of phylo (2021-11-12, Fri)

[TREG](/packages/TREG)
----

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- Your main changes to a function foo() or parameter param.

BUG FIXES

- Your bug fixes. See more details at
http://bioconductor.org/developers/package-guidelines/#news.

[tricycle](/packages/tricycle)
--------

                        Changes in version 1.3.4                        

- Update tricycle citation.

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.14.0                        

- Allow GTF specification in linkedTxome to be a serialized
  GRanges file (a file path to a .rda or .RData file). This
  bypasses some apparent issue where makeTxDbFromGFF fails
  while makeTxDbFromGRanges works.

- Up to GENCODE 40 (H.s.), M29 (M.m), and Ensembl 106

[UCell](/packages/UCell)
-----

                        Changes in version 2.0.0                        

- Update code to pass all BioC checks.

- The function ScoreSignatures_UCell() and StoreRankings_UCell()
  accept directly sce objects.

- Takes custom BiocParallel::bpparam() object as input to specify
  parallelisation.

                        Changes in version 1.3.1                        

- Restructure code to conform to BioC standards.

- Switch from future to BiocParallel to parallelize jobs.

- Add support for SingleCellExperiment - new function
  ScoreSignatures_UCell_sce() interacts directly with sce
  objects.

- Signatures cannot be larger than maxRank parameter.

- Do not rank more genes (maxRank) than there are in the input
  matrix.

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.14.0                        

NEW FEATURES

- enrich_motifs(mode, pseudocount): Choose whether to count motif hits
  once per sequence, and whether to add a pseudocount for P-value
  calculation.

- New function, meme_alph(): Create MEME custom alphabet definition
  files.

- merge_similar(return.clusters): Return the clusters without merging.

- convert_motifs(): MotifDb-MotifList now available as an output
  format.

MINOR CHANGES

- enrich_motifs(): RC argument now defaults to TRUE, increased max
  significance values, no.overlaps now defaults to TRUE. Additional
  columns showing the motif consensus sequence and percent of sequences
  with hits are now included.

- scan_sequences(RC): Only print a warning if RC=TRUE for non-DNA/RNA
  motifs.

- Reduced the size of the message when a pseudocount is added to
  motifs.

                       Changes in version 1.12.4                        

BUG FIXES

- convert_motifs(): Properly handle TFBSTools class motifs with '*' as
  their strand. This was achieved by making the universalmotif object
  creator tolerant to using '*' as a user input. Thanks to David Oliver
  for the bug report (#22).

                       Changes in version 1.12.3                        

BUG FIXES

- scan_sequences(): Previously this function did not account for the
  fact
  that duplicate sequence names are allowed within XStringSet objects.
  To
  better keep track of which sequence hits are associated with, an
  additional sequence.i column has been added which keeps track of the
  sequence number. This change also fixes a knock-on issue with
  enrich_motifs(), where sequences with duplicate names did not
  contribute
  to the count of sequences containing hits. Thanks to Alexandre Blais
  for
  mentioning this issue.

                       Changes in version 1.12.2                        

BUG FIXES

- shuffle_sequences(..., method="markov"): Previously the returning
  sequences were longer by 1.

                       Changes in version 1.12.1                        

BUG FIXES

- DNA ambiguity letters can be used with create_motif() when
  alphabet="DNA" is specified. Previously ambiguity letters only worked
  when alphabet was not specified.

[updateObject](/packages/updateObject)
------------

                        Changes in version 1.0.0                        

- First version of the package that is ready for general use.

[variancePartition](/packages/variancePartition)
-----------------

                       Changes in version 1.25.13                       

- Fix compatibility issue with lme4 1.1.29
- reported
https://github.com/GabrielHoffman/variancePartition/issues/51

                       Changes in version 1.25.12                       

- in makeContrastsDream(), fix issue where terms with colon cause and
error

                       Changes in version 1.25.11                       

- fix bug in dream() for variables with NA values
- improve handling of invalid contrasts in makeContrastsDream()

                       Changes in version 1.25.9                        

- for getContrast() and makeContrastsDream() make sure formula
argument is a formula and not a string

                       Changes in version 1.25.8                        

- small bug fixes

                       Changes in version 1.25.7                        

- dream() now drops samples with missing data gracefully

                       Changes in version 1.25.6                        

- fix small plotting bug in plotStratify() and plotStratifyBy()

                       Changes in version 1.25.5                        

- add getTreat() to evaluate treat()/topTreat() seamlessly on results
of dream()

                       Changes in version 1.25.4                        

- in dream() set default ddf = "adaptive", which uses "KR" for less
than 12 samples
- all functions default tp BPPARAM=SerialParam()
- add eBayes() to vignette for dream()

                       Changes in version 1.25.3                        

- add genes argument to plotPercentBars()

                       Changes in version 1.25.2                        

- change plotPercentBars() to use generic S4

                       Changes in version 1.25.1                        

- update handling of weights in voomWithDreamWeights() and add
applyQualityWeights()

[VaSP](/packages/VaSP)
----

                 Changes in version 1.7.2 (2022-01-07)                  

- Removed the visualizing function splicePlot due to the dependent
  package Sushi deprecated.

[velociraptor](/packages/velociraptor)
------------

                        Changes in version 1.5.2                        

- Remove column names of reduced dimension representation before
velocity embedding.

                        Changes in version 1.5.1                        

- Add example for scvelo.params argument.

[vissE](/packages/vissE)
-----

                        Changes in version 1.4.0                        

- Added the adjusted rand index as a measure of gene-set overlap (now
the default measure of overlap).
- Added feature to only plot gene-sets that are marked in the
plotMsigNetwork function.
- Added function to identify and prioritise gene-set clusters
(findMsigClusters)

[weitrix](/packages/weitrix)
-------

                        Changes in version 1.7.1                        

- counts_shifts and counts_proportions now have a "typecast"
argument,
allowing use of memory-efficient matrix types.

[xcms](/packages/xcms)
----

                       Changes in version 3.17.6                        

- Rewrite code to subset features and chromatographic peaks. This
  results in a
  perfomance improvement for `filterFile` and similar functions.

- Add parameter `expandMz` to `featureChromatograms`
  (https://github.com/sneumann/xcms/issues/612).

                       Changes in version 3.17.5                        

- Change the way the m/z value for a chromatographic peak is determined
  by
  centWave: if a ROI contains more than one peak for one scan
  (spectrum) an
  intensity-weighted m/z is reported for that scan. The m/z of the
  chromatographic peak is then calculated based on these reported m/z
  values for
  each scan (spectrum). In the original version the mean m/z for a scan
  was
  reported instead. As a result, m/z values of chromatographic peaks
  are now
  slightly different but are expected to be more accurate. See
  https://github.com/sneumann/xcms/issues/590 for more details.

                       Changes in version 3.17.4                        

- Add `transformIntensity` method.

- Fix issue when calling `chromPeakSpectra` or `featureSpectra` on an
  object
  that contains also files with only MS1 spectra
  (https://github.com/sneumann/xcms/issues/603).

                       Changes in version 3.17.2                        

- Use mzML instead of mzData files in testing and vignettes,
  since mzR drop mzData reading and msdata package will drop mzData
  files as well

                       Changes in version 3.17.1                        

- Fix bug in feature grouping by EIC correlation that would return a
  non-symmetric similarity matrix.

- Fix error message from issue
  [584](https://github.com/sneumann/xcms/issues/584).

[zellkonverter](/packages/zellkonverter)
-------------

                        Changes in version 1.6.0                        

Major changes

- Added support for multiple *basilisk* environments with
  different *anndata* versions. Users can now specify the
  environment to use with options in readH5AD() and writeH5AD().
  To faciliate this some exported objects where converted to
  functions but this should only effect developers.

- Updated the default environment to use *anndata* v0.8.0. This
  is a major update and files written with v0.8.0 cannot be read
  by previous *anndata* versions. This was the motivation for
  supporting multiple environments and users can select the
  previous environment with *anndata* v0.7.6 if compatibility is
  required.

- Standardise naming in AnnData2SCE(). Column names of data
  frames and names of list items will now be modified to match R
  conventions (according to make.names()). When this happens a
  warning will be issued listing the modifications. This makes
  sure than everything in the created SingleCellExperiment is
  accessible.

Minor changes

- Allow data.frame's stored in varm to be converted in
  SCE2AnnData()

- Minor updates to the vignette and other documentation.

- Updates to tests to match the changes above.


NEWS from new and existing Data Experiment Packages
===================================

[crisprScoreData](/packages/crisprScoreData)
---------------

                       Changes in version 0.99.0                        

- Package submission

[curatedMetagenomicData](/packages/curatedMetagenomicData)
----------------------

                        Changes in version 3.4.0                        

- curatedMetagenomicData now contains 20,533 samples from 90 studies
- A total of 251 samples added since Bioconductor 3.14 (October 2021)
- Studies added since Bioconductor 3.14 (October 2021):
- FrankelAE_2017 (39 samples)
- LeeKA_2022 (165 samples)
- PetersBA_2019 (27 samples)
- WindTT_2020 (20 samples)
- Both "short" and "NCBI" row names were re-validated against NCBI
Taxonomy

[depmap](/packages/depmap)
------

                         Changes in version 1.9                         

Changes in version 1.9.2

- 22Q1 data added for crispr, copyNumber, TPM, mutationCalls and
metadata datasets. New datasets were added, including gene_summaries
and achilles describing Depmap Achilles screens and gene
essentiality probabilities, respectively. New loading functions were
created for these datasets. Newer versions for the other datasets
were not released.

Changes in version 1.9.1

- 21Q4 data added for crispr, copyNumber, TPM, mutationCalls and
metadata datasets. Newer versions for the other datasets were not
released.

Changes in version 1.9.0

- New devel version for Bioc 3.15

Changes in version 1.7.1

- 21Q3 data added for crispr, copyNumber, TPM, mutationCalls and
metadata datasets. Newer versions for the other datasets were not
released.
- CERES CRISPR data has been deprecated and has been replaced with
Chronos CRISPR dependency in 21Q3 and all future releases. For more
information, see:
https://cancerdatascience.org/blog/posts/ceres-chronos/

Changes in version 1.7.0

- New Bioc devel release.

[dorothea](/packages/dorothea)
--------

                 Changes in version 1.7.1 (2022-04-07)                  

- Removed outdated vignettes, now instead we point to decoupleR and
decoupler-py's most up-to-date vignettes.

[epimutacionsData](/packages/epimutacionsData)
----------------

                       Changes in version 0.99.7                        

- Bugs and notes (if possible) fixed

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[FlowSorted.Blood.EPIC](/packages/FlowSorted.Blood.EPIC)
---------------------

                 Changes in version 1.99.1 (2022-01-04)                 

- Made the following significant changes
  o The IDOL libraries are automatically selected
  o The slot counts is now only for cell counts, if available
  o Projections are based on cell proportions
  o estimateCellCounts2 supports additional packages and extended
  deconvolution
  o The reference library is called using libraryDataGet

[healthyControlsPresenceChecker](/packages/healthyControlsPresenceChecker)
------------------------------

                Changes in version 0.99.17 (2021-11-29)                 

- Fixed YAML

                Changes in version 0.99.15 (2021-11-26)                 

- Fixed additional minor issues

                Changes in version 0.99.13 (2021-11-24)                 

- Fixed additional minor issues

                Changes in version 0.99.12 (2021-11-22)                 

- Fixed minor issues

                Changes in version 0.99.11 (2021-11-19)                 

- Added unit tests

                Changes in version 0.99.10 (2021-11-18)                 

- Edited vignettes

                 Changes in version 0.99.9 (2021-11-18)                 

- Addressed points of the review and added a README

                 Changes in version 0.99.1 (2021-10-29)                 

- First release

                 Changes in version 0.1.6 (2021-11-02)                  

- Added print of the percentages of the elements of the healthy
  controls and of the other classes

[msdata](/packages/msdata)
------

                       Changes in version 0.35.3                        

- Adding CE-MS test data, thanks to Liesa Salzer

[pRolocdata](/packages/pRolocdata)
----------

                       Changes in version 1.33.1                        

- Add hyperLOPIT2015_se data

- Add mulvey2015_se and mulvey2015norm_se data

- Regenerated README to include new datasets

[scpdata](/packages/scpdata)
-------

                         Changes in version 1.3                         

scpdata 1.3.1

- williams2020: added 2 new datasets (need to be sent to EH)
- all: added supplementary metadata fields in metadata.csv.
- zhu2019EL: split protein data into protein_intensity and
protein_iBAQ
- dou2019: modified peptide data by excluding low-quality PSMs

scpdata 1.3.0

- New devel (Bioc 3.15)

[signatureSearchData](/packages/signatureSearchData)
-------------------

                 Changes in version 1.9.3 (2021-12-16)                  

- Add LINCS2 database with file name of lincs2020.h5

                 Changes in version 1.9.2 (2021-12-06)                  

- Add dest_path parameter to getCmapCEL function

[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.7.19                        

SIGNIFICANT USER-VISIBLE CHANGES

- Documentation of the layer-level data panel at run_app() has been
significantly increased. You can now also visualize more than 2
reduced dimensions computed on the pseudo-bulk level data
(layer-level for the Maynard et al, Nature Neurosci, 2021 data).
- Users can now control the font and point size on the reduced
dimension plots, as well as the overall font size on the model
boxplots.
- Image edit scenarios you might be interested in for having a
uniform
color background image are now documented; for example if you want a
white or black background, or actually any valid R color name or
color HEX value.

                       Changes in version 1.7.18                        

SIGNIFICANT USER-VISIBLE CHANGES

- run_app() now offers the option to chose any of the
paletteer::paletteer_d color palettes for discrete variables.
- Polychrome has been replaced as a dependency by paletteer. Note
that
Polychrome::palette36 is still the default.
- run_app() now looks for columns that end with '_colors' in their
name which can be used to pre-specify colors for any companion
variables. For example if you have spe$my_groups and
spe$my_groups_colors then the second one can specify the colors that
will be used for visualizing spe$my_groups. This makes specifying
default colors more flexible than before, and the user is still free
to change them if necessary.

                       Changes in version 1.7.17                        

BUG FIXES

- Fix bugs in layer_boxplot() where it was too specific to the
Maynard
et al 2021 data. We have made it more flexible now.
- Made the y-axis space more dynamic in gene_set_enrichment_plot()
and
layer_matrix_plot().

                       Changes in version 1.7.16                        

BUG FIXES

- Fixed a bug in sig_genes_extract() when there's only one set of
t-statistics or F statistics to extract.

                       Changes in version 1.7.12                        

SIGNIFICANT USER-VISIBLE CHANGES

- The visualization functions vis_*() of SpatialLIBD in this version
match the Bioconductor 3.15 version of SpatialExperiment. Note that
if you used SpatialExperiment::read10xVisium(), the names of the
spatial coordinates changed at
https://github.com/drighelli/SpatialExperiment/commit/6710fe8b0a7919191ecce989bb6831647385ef5f
and thus you might need to switch them back if you created your
SpatialExperiment object before this change. You can do so with
spatialCoordsNames(spe) <- rev(spatialCoordsNames(spe)).
read10xVisiumWrapper() uses SpatialExperiment::read10xVisium()
internally, so this change on SpatialExperiment would then also
affect you.

                       Changes in version 1.7.11                        

NEW FEATURES

- Now layer_stat_cor() has the top_n argument which can be used for
subsetting the marker genes prior to computing the correlation as
part of the spatial registration process.

                       Changes in version 1.7.10                        

NEW FEATURES

- Added the add_key() function to reduce code duplication and resolve
https://github.com/LieberInstitute/spatialLIBD/issues/31.

                        Changes in version 1.7.9                        

NEW FEATURES

- This version is now compatible with the bioc-devel version of
SpatialExperiment where spatialData() was deprecated. Details at
https://github.com/LieberInstitute/spatialLIBD/pull/29/files.

                        Changes in version 1.7.7                        

BUG FIXES

- Fixed a bug where the using the left-mouse click was not working
for
annotating individual spots under the "gene (interactive)" tab.

                        Changes in version 1.7.6                        

NEW FEATURES

- vis_gene_p(), vis_clus_p() and all related functions now have an
argument point_size which lets you control how big the points are
plotted. This can be useful for visualization purposes.
- The shiny app now has an input controlling the point size. If you
increase it to say 5, then if you zoom in the clusters (interactive)
panel, you can see larger spots when zooming in.
- These features are related to
https://github.com/LieberInstitute/spatialLIBD/issues/28 although
the spot diameter is still not the true spot diameter. However, now
you have more flexibility for visualizing the spots.

                        Changes in version 1.7.5                        

NEW FEATURES

- Expanded the Using spatialLIBD with 10x Genomics public datasets
vignette to show how you can deploy your web application. See
https://libd.shinyapps.io/spatialLIBD_Human_Lymph_Node_10x/ for the
live example.

                        Changes in version 1.7.4                        

BUG FIXES

- vis_gene() and vis_grid_gene() now support geneids that are found
in
the rownames(spe). This makese these functions more flexible.
- vis_grid_gene() and vis_grid_clus() now have the sample_order
argument which gives you more control in case you want to plot a
subset of samples. This should also reduced the memory required as
discovered at
https://github.com/LieberInstitute/spatialDLPFC/issues/45.

                        Changes in version 1.7.3                        

NEW FEATURES

- Added support for more than one background picture per sample. This
was done through the new argument image_id. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/25.
- Added options for side by side visualization of the background
image
and the clusters or gene expression values in the static versions.
Resolves https://github.com/LieberInstitute/spatialLIBD/issues/19.
- Allow changing the transparency level of the spots with the alpha
argument. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/20.
- Add support for image manipulation with the magick package. Adds
functions img_edit(), img_update() and img_update_all() as well as
new features on the web application. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/21.
- Added support for more control over the gene color scale and in the
web application also added support for reversing the order of the
scale. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/22 and
https://github.com/LieberInstitute/spatialLIBD/issues/23.
- Added export_cluster() and import_cluster() to help export/import
clustering results instead of having to save large spe objects when
exploring different clustering methods.
- Added locate_images() and add_images() for adding non-standard
images to a spe object.

                        Changes in version 1.7.2                        

BUG FIXES

- Fixed an issue introduced by newer versions of shiny. This version
of spatialLIBD works with shiny version 1.7.1, though it's likely
backwards compatible. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/24.
- Fix an issue where as.data.frame(colData(spe)) uses check.names =
TRUE by default and then changes the column names unintentionally.

                        Changes in version 1.7.1                        

NEW FEATURES

- Added read10xVisiumWrapper() and related functions that make it
easier to read in the SpaceRanger output files and launch a shiny
web application using run_app(). These new functions read in the
analysis output from SpaceRanger by 10x Genomics, in particular, the
clustering and dimension reduction (projection) results.

[STexampleData](/packages/STexampleData)
-------------

                 Changes in version 1.3.3 (2022-01-31)                  

- add new datasets ST_mouseOB, SlideSeqV2_mouseHPC

- reformat datasets to SpatialExperiment version 1.5.2

[tuberculosis](/packages/tuberculosis)
------------

                        Changes in version 1.2.0                        

- GEO Series GSE126614 is now available through tuberculosis
- GEO Series GSE152532 is now available through tuberculosis
- GEO Series GSE174552 is now available through tuberculosis
- GEO Series GSE183912 is now available through tuberculosis
- GEO Series GSE184241 is now available through tuberculosis
- GEO Series GSE190024 is now available through tuberculosis
- GEO Series GSE190850 is now available through tuberculosis


NEWS from new and existing Workflows
===================================

[GeoMxWorkflows](/packages/GeoMxWorkflows)
--------------

                        Changes in version 1.1.2                        

- Update install instructions, decrease R version, add links to other
packages, ensure compatibility with updated GeomxTools

                        Changes in version 1.1.1                        

- Bug Fix: Header in vignette, remove self contained exception


Deprecated and Defunct Packages
===============================

Twenty one software packages were removed from this release (after being deprecated
in Bioc 3.14): 
affyPara, ALPS, alsace, BrainStars, dualKS, ENCODExplorer,
ENVISIONQuery, FindMyFriends, GeneAnswers, gramm4R, KEGGprofile,
MSGFgui, MSGFplus, MSstatsTMTPTM, PanVizGenerator, predictionet, RGalaxy,
scClassifR, slinky, SRGnet, SwimR

Please note:  destiny and MouseFM, previously announced as deprecated in
3.14, fixed their packages and remained in Bioconductor.

Twenty nine software packages are deprecated in this release and will be removed in Bioc 3.16:
ABAEnrichment, Autotuner, CAnD, caOmicsV, CHETAH, clonotypeR, CountClust,
diffloop, GCSConnection, GCSFilesystem, GenoGAM, genphen, gprege, networkBMA,
Onassis, perturbatr, phemd, ppiStats, ProteomicsAnnotationHubData, PSICQUIC,
PubScore, Rgin, RmiR, RpsiXML, ScISI, SLGI, Sushi, tofsims, TSRchitect


Three experimental data packages were removed from this release (after being
deprecated in BioC 3.14):
ABAData, brainImageRdata, tcgaWGBSData.hg19

Please note: PCHiCdata and RITANdata previously announced as deprecated in
3.14, fixed their packages and remained in Bioconductor.

Three experimental data packages are deprecated in this release and will be
removed in Bioc 3.16:
DREAM4, MSstatsBioData, ppiData

One annotation packages was removed from this release (after being deprecated
in Bioc 3.14): 
org.Pf.plasmo.db

One annotation package was deprecated in this release and will be removed in
Bioc 3.16:
MafH5.gnomAD.v3.1.1.GRCh38_3.13.1.tar.gz

No workflow packages were removed from this release (after being deprecated in
Bioc 3.14)

One workflow package was deprecated in this release and will be removed in
3.16:
proteomics
