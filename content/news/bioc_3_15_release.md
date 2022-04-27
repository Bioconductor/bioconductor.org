April 27, 2022

Bioconductors:

We are pleased to announce Bioconductor 3.15, consisting of
2123 software packages, 409 experiment data packages, 909 annotation
packages, 29 workflows and 8 books.

There are 61 new software packages, 4 new data experiment packages,
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

There are 61 new software packages in this release of Bioconductor.

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

- [epimutacions](/packages/epimutacions) The package includes some
  statistical outlier detection methods for epimutations detection in
  DNA methylation data. The methods included in the package are
  MANOVA, Multivariate linear models, isolation forest, robust
  mahalanobis distance, quantile and beta. The methods compare a case
  sample with a suspected disease against a reference panel (composed
  of healthy individuals) to identify epimutations in the given case
  sample. It also contains functions to annotate and visualize the
  identified epimutations.

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

- [MBECS](/packages/MBECS) The Microbiome Batch Effect Correction
  Suite (MBECS) provides a set of functions to evaluate and mitigate
  unwated noise due to processing in batches. To that end it
  incorporates a host of batch correcting algorithms (BECA) from
  various packages. In addition it offers a correction and reporting
  pipeline that provides a preliminary look at the characteristics of
  a data-set before and after correcting for batch effects.

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

- [standR](/packages/standR) standR is an user-friendly R package
  providing functions to assist conducting good-practice analysis of
  Nanostring's GeoMX DSP data. All functions in the package are built
  based on the SpatialExperiment object, allowing integration into
  various spatial transcriptomics-related packages from Bioconductor.
  standR allows data inspection, quality control, normalization,
  batch correction and evaluation with informative visualizations.

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

There are 4 new data experiment packages in this release of Bioconductor.

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

TODO


NEWS from new and existing Data Experiment Packages
===================================

TODO


NEWS from new and existing Workflows
===================================

TODO


NEWS from new and existing books
===================================

No new NEWS to report


Deprecated and Defunct Packages
===============================

TODO

