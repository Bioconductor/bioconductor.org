October 27, 2021

Bioconductors:

We are pleased to announce Bioconductor 3.14, consisting of 
2083 software packages, 408 experiment data packages,  904 annotation
packages, 29 workflows and 8 books.

There are 89 new software packages, 13 new data experiment packages,
10 new annotation packages, 1 new workflow, no new books, and many updates and
improvements to existing packages; Bioconductor 3.14 is compatible with R 4.1.1,
and is supported on Linux, 32- and 64-bit Windows, and macOS 10.14.6 Mojave
or higher.  This release will include an updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3]
for details and downloads.

Bioconductor used Microsoft Azure VMs during our 3.14 release process for a
critical part of our branching process for software packages. These VMs are
available to Bioconductor through our collaboration with the Microsoft Genomics
team.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.14](#getting-started-with-bioconductor-314)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from new and existing software packages](#news-from-new-and-existing-software-packages)
* [NEWS from new and existing data experiment packages](#news-from-new-and-existing-data-experiment-packages)
* [NEWS from new and existing workflows](#news-from-new-and-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.14
======================================

To update to or install Bioconductor 3.14:

1. Install R 4.1.1. Bioconductor 3.14 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).

New Software Packages
=====================

There are 89 new software packages in this release of Bioconductor.

- [atena](/packages/atena) Quantify expression of transposable
  elements (TEs) from RNA-seq data through different methods,
  including ERVmap, TEtranscripts and Telescope. A common interface
  is provided to use each of these methods, which consists of
  building a parameter object, calling the quantification function
  with this object and getting a SummarizedExperiment object as
  output container of the quantified expression profiles. The
  implementation allows one to quantify TEs and gene transcripts in
  an integrated manner.

- [benchdamic](/packages/benchdamic) Starting from a microbiome
  dataset (16S or WMS with absolute count values) it is possible to
  perform several analysis to assess the performances of many
  differential abundance detection methods. A basic and standardized
  version of the main differential abundance analysis methods is
  supplied but the user can also add his method to the benchmark. The
  analyses focus on 4 main aspects: i) the goodness of fit of each
  method's distributional assumptions on the observed count data, ii)
  the ability to control the false discovery rate, iii) the within
  and between method concordances, iv) the truthfulness of the
  findings if any apriori knowledge is given. Several graphical
  functions are available for result visualization.

- [BindingSiteFinder](/packages/BindingSiteFinder) Precise knowledge
  on the binding sites of an RNA-binding protein (RBP) is key to
  understand (post-) transcriptional regulatory processes. Here we
  present a workflow that describes how exact binding sites can be
  defined from iCLIP data. The package provides functions for binding
  site definition and result visualization. For details please see
  the vignette.

- [biodbChebi](/packages/biodbChebi) The biodbChebi library provides
  access to the ChEBI Database, using biodb package framework. It
  allows to retrieve entries by their accession number. Web services
  can be accessed for searching the database by name, mass or other
  fields.

- [biodbHmdb](/packages/biodbHmdb) The biodbHmdb library is an
  extension of the biodb framework package that provides access to
  the HMDB Metabolites database. It allows to download the whole HMDB
  Metabolites database locally, access entries and search for entries
  by name or description. A future version of this package will also
  include a search by mass and mass spectra annotation.

- [biodbKegg](/packages/biodbKegg) The biodbKegg library is an
  extension of the biodb framework package that provides access to
  the KEGG databases Compound, Enzyme, Genes, Module, Orthology and
  Reaction. It allows to retrieve entries by their accession numbers.
  Web services like "find", "list" and "findExactMass" are also
  available. Some functions for navigating along the pathways have
  also been implemented.

- [biodbLipidmaps](/packages/biodbLipidmaps) The biodbLipidmaps
  library provides access to the Lipidmaps Structure Database, using
  biodb package framework. It allows to retrieve entries by their
  accession number, and run web the services lmsdSearch and
  lmsdRecord.

- [biodbUniprot](/packages/biodbUniprot) The biodbUniprot library is
  an extension of the biodb framework package. It provides access to
  the UniProt database. It allows to retrieve entries by their
  accession number, and run web service queries for searching for
  entries.

- [BioPlex](/packages/BioPlex) The BioPlex package implements access
  to the BioPlex protein-protein interaction networks and related
  resources from within R. Besides protein-protein interaction
  networks for HEK293 and HCT116 cells, this includes access to CORUM
  protein complex data, and transcriptome and proteome data for the
  two cell lines. Functionality focuses on importing the various data
  resources and storing them in dedicated Bioconductor data
  structures, as a foundation for integrative downstream analysis of
  the data.

- [bugsigdbr](/packages/bugsigdbr) The bugsigdbr package implements
  convenient access to bugsigdb.org from within R/Bioconductor. The
  goal of the package is to facilitate import of BugSigDB data into
  R/Bioconductor, provide utilities for extracting microbe
  signatures, and enable export of the extracted signatures to plain
  text files in standard file formats such as GMT.

- [BUSseq](/packages/BUSseq) BUSseq R package fits an interpretable
  Bayesian hierarchical model---the Batch Effects Correction with
  Unknown Subtypes for scRNA seq Data (BUSseq)---to correct batch
  effects in the presence of unknown cell types. BUSseq is able to
  simultaneously correct batch effects, clusters cell types, and
  takes care of the count data nature, the overdispersion, the
  dropout events, and the cell-specific sequencing depth of scRNA-seq
  data. After correcting the batch effects with BUSseq, the corrected
  value can be used for downstream analysis as if all cells were
  sequenced in a single batch. BUSseq can integrate read count
  matrices obtained from different scRNA-seq platforms and allow cell
  types to be measured in some but not all of the batches as long as
  the experimental design fulfills the conditions listed in our
  manuscript.

- [cageminer](/packages/cageminer) This package aims to integrate
  GWAS-derived SNPs and coexpression networks to mine candidate genes
  associated with a particular phenotype. For that, users must define
  a set of guide genes, which are known genes involved in the studied
  phenotype. Additionally, the mined candidates can be given a score
  that favor candidates that are hubs and/or transcription factors.
  The scores can then be used to rank and select the top n most
  promising genes for downstream experiments.

- [CellBarcode](/packages/CellBarcode) This package performs Cellular
  DNA Barcode (genetic lineage tracing) analysis. The package can
  handle all kinds of DNA barcodes, as long as the barcode within a
  single sequencing read and has a pattern which can be matched by a
  regular expression. This package can handle barcode with flexible
  length, with or without UMI (unique molecular identifier). This
  tool also can be used for pre-processing of some amplicon data such
  as CRISPR gRNA screening, immune repertoire sequencing and meta
  genome data.

- [Cepo](/packages/Cepo) Defining the identity of a cell is
  fundamental to understand the heterogeneity of cells to various
  environmental signals and perturbations. We present Cepo, a new
  method to explore cell identities from single-cell RNA-sequencing
  data using differential stability as a new metric to define cell
  identity genes. Cepo computes cell-type specific gene statistics
  pertaining to differential stable gene expression.

- [cfDNAPro](/packages/cfDNAPro) cfDNA fragment size metrics are
  important features for utilizing liquid biopsy in tumor early
  detection, diagnosis, therapy personlization and monitoring.
  Analyzing and visualizing insert size metrics could be time
  intensive. This package intends to simplify this exploration
  process, and it offers two sets of functions for data
  characterization and data visualization.

- [cliProfiler](/packages/cliProfiler) An easy and fast way to
  visualize and profile the high-throughput IP data. This package
  generates the meta gene profile and other profiles. These profiles
  could provide valuable information for understanding the IP
  experiment results.

- [Cogito](/packages/Cogito) Biological studies often consist of
  multiple conditions which are examined with different laboratory
  set ups like RNA-sequencing or ChIP-sequencing. To get an overview
  about the whole resulting data set, Cogito provides an automated,
  complete, reproducible and clear report about all samples and basic
  comparisons between all different samples. This report can be used
  as documentation about the data set or as starting point for
  further custom analysis.

- [csdR](/packages/csdR) This package contains functionality to run
  differential gene co-expression across two different conditions.
  The algorithm is inspired by Voigt et al. 2017 and finds Conserved,
  Specific and Differentiated genes (hence the name CSD). This
  package include efficient and variance calculation by bootstrapping
  and Welford's algorithm.

- [CyTOFpower](/packages/CyTOFpower) This package is a tool to
  predict the power of CyTOF experiments in the context of
  differential state analyses. The package provides a shiny app with
  two options to predict the power of an experiment: i. generation of
  in-sicilico CyTOF data, using users input ii. browsing in a grid of
  parameters for which the power was already precomputed.

- [cytoKernel](/packages/cytoKernel) cytoKernel implements a
  kernel-based score test to identify differentially expressed
  features in high-dimensional biological experiments. This approach
  can be applied across many different high-dimensional biological
  data including gene expression data and dimensionally reduced
  cytometry-based marker expression data. In this R package, we
  implement functions that compute the feature-wise p values and
  their corresponding adjusted p values. Additionally, it also
  computes the feature-wise shrunk effect sizes and their
  corresponding shrunken effect size. Further, it calculates the
  percent of differentially expressed features and plots
  user-friendly heatmap of the top differentially expressed features
  on the rows and samples on the columns.

- [deconvR](/packages/deconvR) This package provides a collection of
  functions designed for analyzing deconvolution of the bulk
  sample(s) using an atlas of reference omic signature profiles and a
  user-selected model. Users are given the option to create or extend
  a reference atlas and,also simulate the desired size of the bulk
  signature profile of the reference cell types.The package includes
  the cell-type-specific methylation atlas and, Illumina Epic B5
  probe ids that can be used in deconvolution. Additionally,we
  included BSmeth2Probe, to make mapping WGBS data to their probe IDs
  easier.

- [DelayedTensor](/packages/DelayedTensor) DelayedTensor operates
  Tensor arithmetic directly on DelayedArray object. DelayedTensor
  provides some generic function related to Tensor
  arithmetic/decompotision and dispatches it on the DelayedArray
  class. DelayedTensor also suppors Tensor contraction by einsum
  function, which is inspired by numpy einsum.

- [Dino](/packages/Dino) Dino normalizes single-cell, mRNA sequencing
  data to correct for technical variation, particularly sequencing
  depth, prior to downstream analysis. The approach produces a matrix
  of corrected expression for which the dependency between sequencing
  depth and the full distribution of normalized expression; many
  existing methods aim to remove only the dependency between
  sequencing depth and the mean of the normalized expression. This is
  particuarly useful in the context of highly sparse datasets such as
  those produced by 10X genomics and other uninque molecular
  identifier (UMI) based microfluidics protocols for which the
  depth-dependent proportion of zeros in the raw expression data can
  otherwise present a challenge.

- [dStruct](/packages/dStruct) dStruct identifies differentially
  reactive regions from RNA structurome profiling data. dStruct is
  compatible with a broad range of structurome profiling
  technologies, e.g., SHAPE-MaP, DMS-MaPseq, Structure-Seq,
  SHAPE-Seq, etc. See Choudhary et al, Genome Biology, 2019 for the
  underlying method.

- [easier](/packages/easier) This package provides a workflow for the
  use of EaSIeR tool, developed to assess patients' likelihood to
  respond to ICB therapies providing just the patients' RNA-seq data
  as input. We integrate RNA-seq data with different types of prior
  knowledge to extract quantitative descriptors of the tumor
  microenvironment from several points of view, including composition
  of the immune repertoire, and activity of intra- and extra-cellular
  communications. Then, we use multi-task machine learning trained in
  TCGA data to identify how these descriptors can simultaneously
  predict several state-of-the-art hallmarks of anti-cancer immune
  response. In this way we derive cancer-specific models and identify
  cancer-specific systems biomarkers of immune response. These
  biomarkers have been experimentally validated in the literature and
  the performance of EaSIeR predictions has been validated using
  independent datasets form four different cancer types with patients
  treated with anti-PD1 or anti-PDL1 therapy.

- [enhancerHomologSearch](/packages/enhancerHomologSearch) Get ENCODE
  data of enhancer region via H3K4me1 peaks and search homolog
  regions for given sequences. The candidates of enhancer homolog
  regions can be filtered by distance to target TSS. The top
  candidates from human and mouse will be aligned to each other and
  then exported as multiple alignments with given enhancer.

- [epistack](/packages/epistack) The epistack package main objective
  is the visualizations of stacks of genomic tracks (such as, but not
  restricted to, ChIP-seq, ATAC-seq, DNA methyation or genomic
  conservation data) centered at genomic regions of interest.

- [FindIT2](/packages/FindIT2) This package implements functions to
  find influential TF and target based on different input type. It
  have five module: Multi-peak multi-gene annotaion(mmPeakAnno
  module), Calculate regulation potential(calcRP module), Find
  influential Target based on ChIP-Seq and RNA-Seq data(Find
  influential Target module), Find influential TF based on different
  input(Find influential TF module), Calculate peak-gene or peak-peak
  correlation(peakGeneCor module). And there are also some other
  useful function like integrate different source information,
  calculate jaccard similarity for your TF.

- [FLAMES](/packages/FLAMES) Semi-supervised isoform detection and
  annotation from both bulk and single-cell long read RNA-seq data.
  Flames provides automated pipelines for analysing isoforms, as well
  as intermediate functions for manual execution.

- [genomicInstability](/packages/genomicInstability) This package
  contain functions to run genomic instability analysis (GIA) from
  scRNA-Seq data. GIA estimates the association between gene
  expression and genomic location of the coding genes. It uses the
  aREA algorithm to quantify the enrichment of sets of contiguous
  genes (loci-blocks) on the gene expression profiles and estimates
  the Genomic Instability Score (GIS) for each analyzed cell.

- [GeoDiff](/packages/GeoDiff) A series of statistical models using
  count generating distributions for background modelling, feature
  and sample QC, normalization and differential expression analysis
  on GeoMx RNA data. The application of these methods are
  demonstrated by example data analysis vignette.

- [GEOexplorer](/packages/GEOexplorer) GEOexplorer is a Shiny app
  that enables exploratory data analysis and differential gene
  expression of gene expression analysis on microarray gene
  expression datasets held on the GEO database. The outputs are
  interactive graphs that enable users to explore the results of the
  analysis. The development of GEOexplorer was made possible because
  of the excellent code provided by GEO2R (https:
  //www.ncbi.nlm.nih.gov/geo/geo2r/).

- [ggmsa](/packages/ggmsa) A visual exploration tool for multiple
  sequence alignment and associated data. Supports MSA of DNA, RNA,
  and protein sequences using 'ggplot2'. Multiple sequence alignment
  can easily be combined with other 'ggplot2' plots, such as
  phylogenetic tree Visualized by 'ggtree', boxplot, genome map and
  so on. More features: visualization of sequence logos, sequence
  bundles, RNA secondary structures and detection of sequence
  recombinations.

- [ggspavis](/packages/ggspavis) Visualization functions for
  spatially resolved transcriptomics datasets stored in
  SpatialExperiment format. Includes functions to create several
  types of plots for data from from spot-based (e.g. 10x Genomics
  Visium) and molecule-based (e.g. seqFISH) technological platforms.

- [HPiP](/packages/HPiP) HPiP (Host-Pathogen Interaction Prediction)
  uses an ensemble learning algorithm for prediction of host-pathogen
  protein-protein interactions (HP-PPIs) using structural and
  physicochemical descriptors computed from amino acid-composition of
  host and pathogen proteins.The proposed package can effectively
  address data shortages and data unavailability for HP-PPI network
  reconstructions. Moreover, establishing computational frameworks in
  that regard will reveal mechanistic insights into infectious
  diseases and suggest potential HP-PPI targets, thus narrowing down
  the range of possible candidates for subsequent wet-lab
  experimental validations.

- [imcRtools](/packages/imcRtools) This R package supports the
  handling and analysis of imaging mass cytometry and other highly
  multiplexed imaging data. The main functionality includes reading
  in single-cell data after image segmentation and measurement, data
  formatting to perform channel spillover correction and a number of
  spatial analysis approaches. First, cell-cell interactions are
  detected via spatial graph construction; these graphs can be
  visualized with cells representing nodes and interactions
  representing edges. Furthermore, per cell, its direct neighbours
  are summarized to allow spatial clustering. Per image/grouping
  level, interactions between types of cells are counted, averaged
  and compared against random permutations. In that way, types of
  cells that interact more (attraction) or less (avoidance)
  frequently than expected by chance are detected.

- [iPath](/packages/iPath) iPath is the Bioconductor package used for
  calculating personalized pathway score and test the association
  with survival outcomes. Abundant single-gene biomarkers have been
  identified and used in the clinics. However, hundreds of oncogenes
  or tumor-suppressor genes are involved during the process of
  tumorigenesis. We believe individual-level expression patterns of
  pre-defined pathways or gene sets are better biomarkers than single
  genes. In this study, we devised a computational method named iPath
  to identify prognostic biomarker pathways, one sample at a time. To
  test its utility, we conducted a pan-cancer analysis across 14
  cancer types from The Cancer Genome Atlas and demonstrated that
  iPath is capable of identifying highly predictive biomarkers for
  clinical outcomes, including overall survival, tumor subtypes, and
  tumor stage classifications. We found that pathway-based biomarkers
  are more robust and effective than single genes.

- [m6Aboost](/packages/m6Aboost) This package can help user to run
  the m6Aboost model on their own miCLIP2 data. The package includes
  functions to assign the read counts and get the features to run the
  m6Aboost model. The miCLIP2 data should be stored in a GRanges
  object. More details can be found in the vignette.

- [MAI](/packages/MAI) A two-step approach to imputing missing data
  in metabolomics. Step 1 uses a random forest classifier to classify
  missing values as either Missing Completely at Random/Missing At
  Random (MCAR/MAR) or Missing Not At Random (MNAR). MCAR/MAR are
  combined because it is often difficult to distinguish these two
  missing types in metabolomics data. Step 2 imputes the missing
  values based on the classified missing mechanisms, using the
  appropriate imputation algorithms. Imputation algorithms tested and
  available for MCAR/MAR include Bayesian Principal Component
  Analysis (BPCA), Multiple Imputation No-Skip K-Nearest Neighbors
  (Multi_nsKNN), and Random Forest. Imputation algorithms tested and
  available for MNAR include nsKNN and a single imputation approach
  for imputation of metabolites where left-censoring is present.

- [metapone](/packages/metapone) The package conducts pathway testing
  from untargetted metabolomics data. It requires the user to supply
  feature-level test results, from case-control testing, regression,
  or other suitable feature-level tests for the study design. Weights
  are given to metabolic features based on how many metabolites they
  could potentially match to. The package can combine positive and
  negative mode results in pathway tests.

- [methylclock](/packages/methylclock) This package allows to
  estimate chronological and gestational DNA methylation (DNAm) age
  as well as biological age using different methylation clocks.
  Chronological DNAm age (in years) : Horvath's clock, Hannum's
  clock, BNN, Horvath's skin+blood clock, PedBE clock and Wu's clock.
  Gestational DNAm age : Knight's clock, Bohlin's clock, Mayne's
  clock and Lee's clocks. Biological DNAm clocks : Levine's clock and
  Telomere Length's clock.

- [miaSim](/packages/miaSim) Microbiome time series simulation with
  generalized Lotka-Volterra model, Self-Organized Instability (SOI),
  and other models. Hubbell's Neutral model is used to determine the
  abundance matrix. The resulting abundance matrix is applied to
  SummarizedExperiment or TreeSummarizedExperiment objects.

- [microbiomeMarker](/packages/microbiomeMarker) To date, a number of
  methods have been developed for microbiome marker discovery based
  on metagenomic profiles, e.g. LEfSe. However, all of these methods
  have its own advantages and disadvantages, and none of them is
  considered standard or universal. Moreover, different programs or
  softwares may be development using different programming languages,
  even in different operating systems. Here, we have developed an
  all-in-one R package microbiomeMarker that integrates commonly used
  differential analysis methods as well as three machine
  learning-based approaches, including Logistic regression, Random
  forest, and Support vector machine, to facilitate the
  identification of microbiome markers.

- [MicrobiomeProfiler](/packages/MicrobiomeProfiler) This is an
  R/shiny package to perform functional enrichment analysis for
  microbiome data. This package was based on clusterProfiler.
  Moreover, MicrobiomeProfiler support KEGG enrichment analysis, COG
  enrichment analysis, Microbe-Disease association enrichment
  analysis, Metabo-Pathway analysis.

- [mitoClone2](/packages/mitoClone2) This package primarily
  identifies variants in mitochondrial genomes from BAM alignment
  files. It filters these variants to remove RNA editing events then
  estimates their evolutionary relationship (i.e. their phylogenetic
  tree) and groups single cells into clones. It also visualizes the
  mutations and providing additional genomic context.

- [monaLisa](/packages/monaLisa) Useful functions to work with
  sequence motifs in the analysis of genomics data. These include
  methods to annotate genomic regions or sequences with predicted
  motif hits and to identify motifs that drive observed changes in
  accessibility or expression. Functions to produce informative
  visualizations of the obtained results are also provided.

- [mosbi](/packages/mosbi) This package is a implementation of
  biclustering ensemble method MoSBi (Molecular signature
  Identification from Biclustering). MoSBi provides standardized
  interfaces for biclustering results and can combine their results
  with a multi-algorithm ensemble approach to compute robust ensemble
  biclusters on molecular omics data. This is done by computing
  similarity networks of biclusters and filtering for overlaps using
  a custom error model. After that, the louvain modularity it used to
  extract bicluster communities from the similarity network, which
  can then be converted to ensemble biclusters. Additionally, MoSBi
  includes several network visualization methods to give an intuitive
  and scalable overview of the results. MoSBi comes with several
  biclustering algorithms, but can be easily extended to new
  biclustering algorithms.

- [MsBackendRawFileReader](/packages/MsBackendRawFileReader)
  implements a MsBackend for the Spectra package using Thermo Fisher
  Scientific's NewRawFileReader .Net libraries. The package is
  generalizing the functionality introduced by the rawrr package
  (Kockmann T. et al. (2020) <doi:10.1101/2020.10.30.362533>) Methods
  defined in this package are supposed to extend the Spectra
  Bioconductor package.

- [MSstatsLiP](/packages/MSstatsLiP) Tools for LiP peptide and
  protein significance analysis. Provides functions for
  summarization, estimation of LiP peptide abundance, and detection
  of changes across conditions. Utilizes functionality across the
  MSstats family of packages.

- [NanoTube](/packages/NanoTube) NanoTube includes functions for the
  processing, quality control, analysis, and visualization of
  NanoString nCounter data. Analysis functions include differential
  analysis and gene set analysis methods, as well as postprocessing
  steps to help understand the results. Additional functions are
  included to enable interoperability with other Bioconductor
  NanoString data analysis packages.

- [netOmics](/packages/netOmics) netOmics is a multi-omics networks
  builder and explorer. It uses a combination of network inference
  algorithms and and knowledge-based graphs to build multi-layered
  networks. The package can be combined with timeOmics to incorporate
  time-course expression data and build sub-networks from multi-omics
  kinetic clusters. Finally, from the generated multi-omics networks,
  propagation analyses allow the identification of missing biological
  functions (1), multi-omics mechanisms (2) and molecules between
  kinetic clusters (3). This helps to resolve complex regulatory
  mechanisms.

- [NeuCA](/packages/NeuCA) NeuCA is is a neural-network based method
  for scRNA-seq data annotation. It can automatically adjust its
  classification strategy depending on cell type correlations, to
  accurately annotate cell. NeuCA can automatically utilize the
  structure information of the cell types through a hierarchical tree
  to improve the annotation accuracy. It is especially helpful when
  the data contain closely correlated cell types.

- [nullranges](/packages/nullranges) Modular package for generation
  of sets of ranges representing the null hypothesis. These can take
  the form of bootstrap samples of ranges (using the block bootstrap
  framework of Bickel et al 2010), or sets of control ranges that are
  matched across one or more covariates. nullranges is designed to be
  inter-operable with other packages for analysis of genomic overlap
  enrichment, including the plyranges Bioconductor package.

- [NxtIRFcore](/packages/NxtIRFcore) Interactively analyses Intron
  Retention and Alternative Splicing Events (ASE) in RNA-seq data.
  NxtIRF quantifies ASE events in BAM files aligned to the genome
  using a splice-aware aligner such as STAR. The core quantitation
  algorithm relies on the IRFinder/C++ engine ported via Rcpp for
  multi-platform compatibility. In addition, NxtIRF provides
  convenient pipelines for downstream analysis and publication-ready
  visualisation tools.

- [ODER](/packages/ODER) The aim of ODER is to identify previously
  unannotated expressed regions (ERs) using RNA-sequencing data. For
  this purpose, ODER defines and optimises the definition of ERs,
  then connected these ERs to genes using junction data. In this way,
  ODER improves gene annotation. Gene annotation is a staple input of
  many bioinformatic pipelines and a more complete gene annotation
  can enable more accurate interpretation of disease associated
  variants.

- [orthogene](/packages/orthogene) orthogene is an R package for easy
  mapping of orthologous genes across hundreds of species. It pulls
  up-to-date interspecies gene ortholog mappings across 700+
  organisms. It also provides various utility functions to map common
  objects (e.g. data.frames, gene expression matrices, lists) onto
  1:1 gene orthologs from any other species.

- [pairkat](/packages/pairkat) PaIRKAT is model framework for
  assessing statistical relationships between networks of metabolites
  (pathways) and an outcome of interest (phenotype). PaIRKAT queries
  the KEGG database to determine interactions between metabolites
  from which network connectivity is constructed. This model
  framework improves testing power on high dimensional data by
  including graph topography in the kernel machine regression
  setting. Studies on high dimensional data can struggle to include
  the complex relationships between variables. The semi-parametric
  kernel machine regression model is a powerful tool for capturing
  these types of relationships. They provide a framework for testing
  for relationships between outcomes of interest and high dimensional
  data such as metabolomic, genomic, or proteomic pathways. PaIRKAT
  uses known biological connections between high dimensional
  variables by representing them as edges of ‘graphs’ or ‘networks.’
  It is common for nodes (e.g. metabolites) to be disconnected from
  all others within the graph, which leads to meaningful decreases in
  testing power whether or not the graph information is included. We
  include a graph regularization or ‘smoothing’ approach for managing
  this issue.

- [pengls](/packages/pengls) Combine generalised least squares
  methodology from the nlme package for dealing with autocorrelation
  with penalised least squares methods from the glmnet package to
  deal with high dimensionality. This pengls packages glues them
  together through an iterative loop. The resulting method is
  applicable to high dimensional datasets that exhibit
  autocorrelation, such as spatial or temporal data.

- [plotgardener](/packages/plotgardener) Coordinate-based genomic
  visualization package for R. It grants users the ability to
  programmatically produce complex, multi-paneled figures. Tailored
  for genomics, plotgardener allows users to visualize large complex
  genomic datasets and provides exquisite control over how plots are
  placed and arranged on a page.

- [ProteoDisco](/packages/ProteoDisco) ProteoDisco is an R package to
  facilitate proteogenomics studies. It houses functions to create
  customized (mutant) protein databases based on user-submitted
  genomic variants, splice-junctions, fusion genes and manual
  transcript sequences. The flexible workflow can be adopted to suit
  a myriad of research and experimental settings.

- [rGenomeTracks](/packages/rGenomeTracks) rGenomeTracks package
  leverages the power of pyGenomeTracks software with the
  interactivity of R. pyGenomeTracks is a python software that offers
  robust method for visualizing epigenetic data files like
  narrowPeak, Hic matrix, TADs and arcs, however though, here is no
  way currently to use it within R interactive session. rGenomeTracks
  wrapped the whole functionality of pyGenomeTracks with additional
  utilites to make to more pleasant for R users.

- [RiboCrypt](/packages/RiboCrypt) R Package for interactive
  visualization and browsing NGS data. It contains a browser for both
  transcript and genomic coordinate view. In addition a QC and
  general metaplots are included, among others differential
  translation plots and gene expression plots. The package is still
  under development.

- [RLSeq](/packages/RLSeq) RLSeq is a toolkit for analyzing and
  evaluating R-loop mapping datasets. RLSeq serves two primary
  purposes:
  (1) to facilitate the evaluation of dataset quality
  (2) to enable R-loop analysis in the context of publicly-available
  data sets from RLBase. The package is intended to provide a simple
  pipeline, called with the `RLSeq()` function, which performs all
  main analyses. Individual functions are also accessible and provide
  custom analysis capabilities. Finally an HTML report is generated
  with `report()`.


- [rmspc](/packages/rmspc) The rmspc package runs MSPC (Multiple
  Sample Peak Calling) software using R. The analysis of ChIP-seq
  samples outputs a number of enriched regions (commonly known as
  "peaks"), each indicating a protein-DNA interaction or a specific
  chromatin modification. When replicate samples are analyzed,
  overlapping peaks are expected. This repeated evidence can
  therefore be used to locally lower the minimum significance
  required to accept a peak. MSPC uses combined evidence from
  replicated experiments to evaluate peak calling output, rescuing
  peaks, and reduce false positives. It takes any number of
  replicates as input and improves sensitivity and specificity of
  peak calling on each, and identifies consensus regions between the
  input samples.

- [scanMiR](/packages/scanMiR) A set of tools for working with miRNA
  affinity models (KdModels), efficiently scanning for miRNA binding
  sites, and predicting target repression. It supports scanning using
  miRNA seeds, full miRNA sequences (enabling 3' alignment) and
  KdModels, and includes the prediction of slicing and TDMD sites.
  Finally, it includes utility and plotting functions (e.g. for the
  visual representation of miRNA-target alignment).

- [scanMiRApp](/packages/scanMiRApp) A shiny interface to the scanMiR
  package. The application enables the scanning of transcripts and
  custom sequences for miRNA binding sites, the visualization of
  KdModels and binding results, as well as browsing predicted
  repression data. In addition contains the IndexedFst class for fast
  indexed reading of large GenomicRanges or data.frames, and some
  utilities for facilitating scans and identifying enriched
  miRNA-target pairs.

- [scAnnotatR](/packages/scAnnotatR) The package comprises a set of
  pretrained machine learning models to predict basic immune cell
  types. This enables all users to quickly get a first annotation of
  the cell types present in their dataset without requiring prior
  knowledge. scAnnotatR also allows users to train their own models
  to predict new cell types based on specific research needs.

- [scatterHatch](/packages/scatterHatch) The objective of this
  package is to efficiently create scatterplots where groups can be
  distinguished by color and texture. Visualizations in computational
  biology tend to have many groups making it difficult to distinguish
  between groups solely on color. Thus, this package is useful for
  increasing the accessibility of scatterplot visualizations to those
  with visual impairments such as color blindness.

- [scReClassify](/packages/scReClassify) A post hoc cell type
  classification tool to fine-tune cell type annotations generated by
  any cell type classification procedure with semi-supervised
  learning algorithm AdaSampling technique. The current version of
  scReClassify supports Support Vector Machine and Random Forest as a
  base classifier.

- [scShapes](/packages/scShapes) We present a novel statistical
  framework for identifying differential distributions in single-cell
  RNA-sequencing (scRNA-seq) data between treatment conditions by
  modeling gene expression read counts using generalized linear
  models (GLMs). We model each gene independently under each
  treatment condition using error distributions Poisson (P), Negative
  Binomial (NB), Zero-inflated Poisson (ZIP) and Zero-inflated
  Negative Binomial (ZINB) with log link function and model based
  normalization for differences in sequencing depth. Since all four
  distributions considered in our framework belong to the same family
  of distributions, we first perform a Kolmogorov-Smirnov (KS) test
  to select genes belonging to the family of ZINB distributions.
  Genes passing the KS test will be then modeled using GLMs. Model
  selection is done by calculating the Bayesian Information Criterion
  (BIC) and likelihood ratio test (LRT) statistic.

- [scTreeViz](/packages/scTreeViz) scTreeViz provides classes to
  support interactive data aggregation and visualization of single
  cell RNA-seq datasets with hierarchies for e.g. cell clusters at
  different resolutions. The `TreeIndex` class provides methods to
  manage hierarchy and split the tree at a given resolution or across
  resolutions. The `TreeViz` class extends `SummarizedExperiment` and
  can performs quick aggregations on the count matrix defined by
  clusters.

- [segmenter](/packages/segmenter) Chromatin segmentation analysis
  transforms ChIP-seq data into signals over the genome. The latter
  represents the observed states in a multivariate Markov model to
  predict the chromatin's underlying states. ChromHMM, written in
  Java, integrates histone modification datasets to learn the
  chromatin states de-novo. The goal of this package is to call
  chromHMM from within R, capture the output files in an S4 object
  and interface to other relevant Bioconductor analysis tools. In
  addition, segmenter provides functions to test, select and
  visualize the output of the segmentation.

- [sparrow](/packages/sparrow) Provides a unified interface to a
  variety of GSEA techniques from different bioconductor packages.
  Results are harmonized into a single object and can be interrogated
  uniformly for quick exploration and interpretation of results.
  Interactive exploration of GSEA results is enabled through a shiny
  app provided by a sparrow.shiny sibling package.

- [spatialDE](/packages/spatialDE) SpatialDE is a method to find
  spatially variable genes (SVG) from spatial transcriptomics data.
  This package provides wrappers to use the Python SpatialDE library
  in R, using reticulate and basilisk.

- [spatzie](/packages/spatzie) Identifies motifs that are
  significantly co-enriched from enhancer-promoter interaction data.
  While enhancer-promoter annotation is commonly used to define
  groups of interaction anchors, spatzie also supports co-enrichment
  analysis between preprocessed interaction anchors.  Supports BEDPE
  interaction data derived from genome-wide assays such as HiC,
  ChIA-PET, and HiChIP. Can also be used to look for differentially
  enriched motif pairs between two interaction experiments.

- [spiky](/packages/spiky) spiky implements methods and model
  generation for cfMeDIP (cell-free methylated DNA
  immunoprecipitation) with spike-in controls. CfMeDIP is an
  enrichment protocol which avoids destructive conversion of scarce
  template, making it ideal as a "liquid biopsy," but creating
  certain challenges in comparing results across specimens, subjects,
  and experiments. The use of synthetic spike-in standard oligos
  allows diagnostics performed with cfMeDIP to quantitatively compare
  samples across subjects, experiments, and time points in both
  relative and absolute terms.

- [surfaltr](/packages/surfaltr) Cell surface proteins form a major
  fraction of the druggable proteome and can be used for
  tissue-specific delivery of oligonucleotide/cell-based
  therapeutics. Alternatively spliced surface protein isoforms have
  been shown to differ in their subcellular localization and/or their
  transmembrane (TM) topology. Surface proteins are hydrophobic and
  remain difficult to study thereby necessitating the use of TM
  topology prediction methods such as TMHMM and Phobius. However,
  there exists a need for bioinformatic approaches to streamline
  batch processing of isoforms for comparing and visualizing
  topologies. To address this gap, we have developed an R package,
  surfaltr. It pairs inputted isoforms, either known alternatively
  spliced or novel, with their APPRIS annotated principal
  counterparts, predicts their TM topologies using TMHMM or Phobius,
  and generates a customizable graphical output. Further, surfaltr
  facilitates the prioritization of biologically diverse isoform
  pairs through the incorporation of three different ranking metrics
  and through protein alignment functions. Citations for programs
  mentioned here can be found in the vignette.

- [svaNUMT](/packages/svaNUMT) svaNUMT contains functions for
  detecting NUMT events from structural variant calls. It takes
  structural variant calls in GRanges of breakend notation and
  identifies NUMTs by nuclear-mitochondrial breakend junctions. The
  main function reports candidate NUMTs if there is a pair of valid
  insertion sites found on the nuclear genome within a certain
  distance threshold. The candidate NUMTs are reported by events.

- [svaRetro](/packages/svaRetro) svaRetro contains functions for
  detecting retrotransposed transcripts (RTs) from structural variant
  calls. It takes structural variant calls in GRanges of breakend
  notation and identifies RTs by exon-exon junctions and insertion
  sites. The candidate RTs are reported by events and annotated with
  information of the inserted transcripts.

- [synapsis](/packages/synapsis) Synapsis is a Bioconductor software
  package for automated (unbiased and reproducible) analysis of
  meiotic immunofluorescence datasets. The primary functions of the
  software can 
  i) identify cells in meiotic prophase that are
  labelled by a synaptonemal complex axis or central element protein,
  ii) isolate individual synaptonemal complexes and measure their
  physical length,
  iii) quantify foci and co-localise them with
  synaptonemal complexes, 
  iv) measure interference between
  synaptonemal complex-associated foci. The software has applications
  that extend to multiple species and to the analysis of other
  proteins that label meiotic prophase chromosomes. The software
  converts meiotic immunofluorescence images into R data frames that
  are compatible with machine learning methods. Given a set of
  microscopy images of meiotic spread slides, synapsis crops images
  around individual single cells, counts colocalising foci on strands
  on a per cell basis, and measures the distance between foci on any
  given strand.


- [tanggle](/packages/tanggle) Offers functions for plotting split
  (or implicit) networks (unrooted, undirected) and explicit networks
  (rooted, directed) with reticulations extending. 'ggtree' and using
  functions from 'ape' and 'phangorn'. It extends the 'ggtree'
  package [@Yu2017] to allow the visualization of phylogenetic
  networks using the 'ggplot2' syntax. It offers an alternative to
  the plot functions already available in 'ape' Paradis and Schliep
  (2019) <doi:10.1093/bioinformatics/bty633> and 'phangorn' Schliep
  (2011) <doi:10.1093/bioinformatics/btq706>.

- [TargetDecoy](/packages/TargetDecoy) A first step in the data
  analysis of Mass Spectrometry (MS) based proteomics data is to
  identify peptides and proteins. With this respect the huge number
  of experimental mass spectra typically have to be assigned to
  theoretical peptides derived from a sequence database. Search
  engines are used for this purpose. These tools compare each of the
  observed spectra to all candidate theoretical spectra derived from
  the sequence data base and calculate a score for each comparison.
  The observed spectrum is then assigned to the theoretical peptide
  with the best score, which is also referred to as the peptide to
  spectrum match (PSM). It is of course crucial for the downstream
  analysis to evaluate the quality of these matches. Therefore False
  Discovery Rate (FDR) control is used to return a reliable list
  PSMs. The FDR, however, requires a good characterisation of the
  score distribution of PSMs that are matched to the wrong peptide
  (bad target hits). In proteomics, the target decoy approach (TDA)
  is typically used for this purpose. The TDA method matches the
  spectra to a database of real (targets) and nonsense peptides
  (decoys). A popular approach to generate these decoys is to reverse
  the target database. Hence, all the PSMs that match to a decoy are
  known to be bad hits and the distribution of their scores are used
  to estimate the distribution of the bad scoring target PSMs. A
  crucial assumption of the TDA is that the decoy PSM hits have
  similar properties as bad target hits so that the decoy PSM scores
  are a good simulation of the target PSM scores. Users, however,
  typically do not evaluate these assumptions. To this end we
  developed TargetDecoy to generate diagnostic plots to evaluate the
  quality of the target decoy method.

- [transformGamPoi](/packages/transformGamPoi) Variance-stabilizing
  transformations help with the analysis of heteroskedastic data
  (i.e., data where the variance is not constant, like count data).
  This package provide two types of variance stabilizing
  transformations: (1) methods based on the delta method (e.g.,
  'acosh', 'log(x+1)'), (2) model residual based (Pearson and
  randomized quantile residuals).

- [traviz](/packages/traviz) traviz provides a suite of functions to
  plot trajectory related objects from Bioconductor packages. It
  allows plotting trajectories in reduced dimension, as well as
  averge gene expression smoothers as a function of pseudotime.
  Asides from general utility functions, traviz also allows plotting
  trajectories estimated by Slingshot, as well as smoothers estimated
  by tradeSeq. Furthermore, it allows for visualization of Slingshot
  trajectories using ggplot2.

- [TRESS](/packages/TRESS) This package is devoted to analyzing
  MeRIP-seq data. Current functionality is for detection of
  transcriptome-wide m6A methylation regions. The method is based on
  hierarchical negative binomial models.

- [tripr](/packages/tripr) TRIP is a software framework that provides
  analytics services on antigen receptor (B cell receptor
  immunoglobulin, BcR IG | T cell receptor, TR) gene sequence data.
  It is a web application written in R Shiny. It takes as input the
  output files of the IMGT/HighV-Quest tool. Users can select to
  analyze the data from each of the input samples separately, or the
  combined data files from all samples and visualize the results
  accordingly.

- [txcutr](/packages/txcutr) Various mRNA sequencing library
  preparation methods generate sequencing reads specifically from the
  transcript ends. Analyses that focus on quantification of isoform
  usage from such data can be aided by using truncated versions of
  transcriptome annotations, both at the alignment or
  pseudo-alignment stage, as well as in downstream analysis. This
  package implements some convenience methods for readily generating
  such truncated annotations and their corresponding sequences.

- [VAExprs](/packages/VAExprs) A fundamental problem in biomedical
  research is the low number of observations, mostly due to a lack of
  available biosamples, prohibitive costs, or ethical reasons. By
  augmenting a few real observations with artificially generated
  samples, their analysis could lead to more robust and higher
  reproducible. One possible solution to the problem is the use of
  generative models, which are statistical models of data that
  attempt to capture the entire probability distribution from the
  observations. Using the variational autoencoder (VAE), a well-known
  deep generative model, this package is aimed to generate samples
  with gene expression data, especially for single-cell RNA-seq data.
  Furthermore, the VAE can use conditioning to produce specific cell
  types or subpopulations. The conditional VAE (CVAE) allows us to
  create targeted samples rather than completely random ones.

- [veloviz](/packages/veloviz) VeloViz uses each cell’s current
  observed and predicted future transcriptional states inferred from
  RNA velocity analysis to build a nearest neighbor graph between
  cells in the population. Edges are then pruned based on a cosine
  correlation threshold and/or a distance threshold and the resulting
  graph is visualized using a force-directed graph layout algorithm.
  VeloViz can help ensure that relationships between cell states are
  reflected in the 2D embedding, allowing for more reliable
  representation of underlying cellular trajectories.

New Data Experiment Packages
=====================

There are 13 new data experiment packages in this release of Bioconductor.

- [curatedTBData](/packages/curatedTBData) The curatedTBData is an R
  package that provides standardized, curated tuberculosis(TB)
  transcriptomic studies. The initial release of the package contains
  49 studies. The curatedTBData package allows users to access
  tuberculosis trancriptomic efficiently and to make efficient
  comparison for different TB gene signatures across multiple
  datasets.

- [easierData](/packages/easierData) Access to internal data required
  for the functional performance of easier package and exemplary
  bladder cancer dataset with both processed RNA-seq data and
  information on response to ICB therapy generated by Mariathasan et
  al. "TGF-B attenuates tumour response to PD-L1 blockade by
  contributing to exclusion of T cells", published in Nature, 2018
  [doi:10.1038/nature25501](https://doi.org/10.1038/nature25501). The
  data is made available via
  [`IMvigor210CoreBiologies`](http://research-pub.gene.com/IMvigor210CoreBiologies/)
  package under the CC-BY license.

- [GSE103322](/packages/GSE103322) Single cell RNA-Seq data for 5902
  cells from 18 patients with oral cavity head and neck squamous cell
  carcinoma available as GEO accession [GSE103322]
  (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322).
  GSE103322 data have been parsed into a SincleCellExperiment object
  available in ExperimentHub.

- [GSE159526](/packages/GSE159526) 19 term and 9 first trimester
  placental chorionic villi and matched cell-sorted samples ran on
  Illumina HumanMethylationEPIC DNA methylation microarrays. This
  data was made available on GEO accession
  [GSE159526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159526).
  Both the raw and processed data has been made available on
  \code{ExperimentHub}. Raw unprocessed data formatted as an
  RGChannelSet object for integration and normalization using minfi
  and other existing Bioconductor packages. Processed normalized data
  is also available as a DNA methylation \code{matrix}, with a
  corresponding phenotype information as a \code{data.frame} object.

- [nullrangesData](/packages/nullrangesData) Provides datasets for
  the nullranges package vignette, in particular example datasets for
  DNase hypersensitivity sites (DHS), CTCF binding sites, and CTCF
  genomic interactions. These are used to demonstrate generation of
  null hypothesis feature sets, either through block bootstrapping or
  matching, in the nullranges vignette.  For more details, see the
  data object man pages, and the R scripts for object construction
  provided within the package.

- [NxtIRFdata](/packages/NxtIRFdata) NxtIRFdata is a companion
  package for NxtIRF, which is an IRFinder- based R package for
  Intron Retention and Alternative Splicing quantitation for RNA-seq
  BAM files. NxtIRFdata contains Mappability files required for the
  generation of human and mouse references. NxtIRFdata also contains
  a synthetic genome reference and example BAM files used to test
  NxtIRF. BAM files are based on 6 samples from the Leucegene dataset
  provided by NCBI Gene Expression Omnibus under accession number
  GSE67039.

- [plotgardenerData](/packages/plotgardenerData) This is a
  supplemental data package for the plotgardener package. Includes
  example datasets used in plotgardener vignettes and example raw
  data files. For details on how to use these datasets, see the
  plotgardener package vignettes.

- [RLHub](/packages/RLHub) | RLHub provides a convenient interface to
  the processed data provided within RLSuite, a tool-chain for
  analyzing R-loop-mapping data sets. The primary purpose of RLHub is
  to serve the processed data sets required by the RLSeq R package
  and the RLBase web service. Additionally, RLHub provides a
  stand-alone R interface to these data, benefiting users who are
  addressing questions related to R-loop regions (RL-Regions),
  R-loop-binding proteins (RLBPs), R-loop co-localizing factors, and
  the differences between R-loop-mapping methods. The full
  data-generating protocol is found here:
  https://github.com/Bishop-Laboratory/RLBase-data.

- [scanMiRData](/packages/scanMiRData) This package contains
  companion data to the scanMiR package. It contains `KdModel` (miRNA
  12-mer binding affinity models) collections corresponding to all
  human, mouse and rat mirbase miRNAs. See the scanMiR package for
  details.

- [scATAC.Explorer](/packages/scATAC.Explorer) This package provides
  a tool to search and download a collection of publicly available
  single cell ATAC-seq datasets and their metadata. scATAC-Explorer
  aims to act as a single point of entry for users looking to study
  single cell ATAC-seq data. Users can quickly search available
  datasets using the metadata table and download datasets of interest
  for immediate analysis within R.

- [spatialDmelxsim](/packages/spatialDmelxsim) Spatial allelic
  expression counts from Combs & Fraser (2018), compiled into a
  SummarizedExperiment object. This package contains data of allelic
  expression counts of spatial slices of a fly embryo, a Drosophila
  melanogaster x Drosophila simulans cross. See the CITATION file for
  the data source, and the associated script for how the object was
  constructed from publicly available data.

- [TabulaMurisSenisData](/packages/TabulaMurisSenisData) This package
  provides access to RNA-seq data generated by the Tabula Muris Senis
  project via the Bioconductor project. The data is made available
  without restrictions by the Chan Zuckerberg Biohub. It is provided
  here without further processing, collected in the form of
  SingleCellExperiment objects.

- [tuberculosis](/packages/tuberculosis) The tuberculosis
  R/Bioconductor package features tuberculosis gene expression data
  for machine learning. All human samples from GEO that did not come
  from cell lines, were not taken postmortem, and did not feature
  recombination have been included. The package has more than 10,000
  samples from both microarray and sequencing studies that have been
  processed from raw data through a hyper-standardized, reproducible
  pipeline.

New Annotation Packages
=====================

There are 10 new annotation packages in this release of Bioconductor.

- [chromhmmData](/packages/chromhmmData) Annotation files of the formatted
  genomic annotation for ChromHMM. Three types of text files are included the
  chromosome sizes, region coordinates and anchors specifying the transcription
  start and end sites. The package includes data for two versions of the genome
  of humans and mice.
  
- [CTCF](/packages/CTCF) Genomic coordinates of predicted CTCF binding sites
  with motif MA0139.1 (Jaspar), in BED format. With strand orientation
  (directionality of binding). Human (hg19, hg38) and mouse (mm9, mm10)
  genomes. The binding sites were detected using the FIMO tool of the MEME suite
  using default settings. Extra columns include motif name (MA0139.1), score,
  p-value, q-value, and the motif sequence.
  
- [excluderanges](/packages/excluderanges) Genomic coordinates of problematic
  genomic regions that should be avoided when working with genomic data. GRanges
  of exclusion regions (formerly known as blacklisted), centromeres, telomeres,
  known heterochromatin regions, etc. (UCSC 'gap' table data). Primarily for
  human and mouse genomes, hg19/hg38 and mm9/mm10 genome assemblies.
  
- [GeneSummary](/packages/GeneSummary) This package provides long description of
  genes collected from the RefSeq database. The text in "COMMENT" section
  started with "Summary" is extracted as the description of the gene. The long
  text descriptions can be used for analysis such as text mining.
  
- [ontoProcData](/packages/ontoProcData) This package manages rda files of
  multiple ontologies that are used in the ontoProc package. These ontologies
  were originally downloaded as owl or obo files and converted into Rda
  files. The files were downloaded at various times but most of them were
  downloaded on June 22 2021.
  
- [rGenomeTracksData](/packages/rGenomeTracksData) rGenomeTracksData is a
  collection of data from pyGenomeTracks project. The purpose of this data is
  testing and demonstration of rGenomeTracks. This package include 14 sample
  file from different genomic and epigenomic file format.
  
- [scAnnotatR.models](/packages/scAnnotatR.models) Pretrained models for
  scAnnotatR package. These models can be used to automatically classify several
  (immune) cell types in human scRNA-seq data.
  
- [synaptome.data](/packages/synaptome.data) The package provides access to the
  copy of the Synaptic proteome database. It was designed as an accompaniment
  for Synaptome.DB package. Database provides information for specific genes and
  allows building the protein-protein interaction graph for gene sets, synaptic
  compartments, and brain regions.
  
- [synaptome.db](/packages/synaptome.db) The package contains local copy of the
  Synaptic proteome database. On top of this it provide a set of utility R
  functions to query and analyse its content. It allows extraction of
  information for specific genes and building the protein-protein interaction
  graph for gene sets, synaptic compartments, and brain regions.
  
- [TxDb.Athaliana.BioMart.plantsmart51](/packages/TxDb.Athaliana.BioMart.plantsmart51)
  Exposes an annotation databases generated from BioMart by exposing these as
  TxDb objects. This package is for Arabidopsis thaliana (taxID: 3702). The
  BioMart plantsmart release number is 51.
  
New Workflow Packages
=====================

There is 1 new workflow package in this release of Bioconductor. 

- [GeoMxWorkflows](/packages/GeoMxWorkflows) Workflows for use with
  NanoString Technologies GeoMx Technology.  Package provides
  bioconductor focused workflows for leveraging existing packages
  (e.g. GeomxTools) to process, QC, and analyze the data.


New Books
=====================

There are no new online books.


NEWS from new and existing Software Packages
===================================


[a4Classif](/packages/a4Classif)
---------

                       Changes in version 1.41.1                        

- vignette: select subset ALL data without replacement

[affxparser](/packages/affxparser)
----------

                 Changes in version 1.65.3 (2021-09-22)                 

SOFTWARE QUALITY

- Making sure all pathnames are of length 100 or shorter.

                 Changes in version 1.65.2 (2021-09-22)                 

SOFTWARE QUALITY

- Now properly registering native routines.

                 Changes in version 1.65.1 (2021-09-09)                 

BUG FIXES

- The package did not install on macOS with the M1 chip with error:
  "use of
  undeclared identifier 'finite'; did you mean 'isfinite'?". This issue
  goes
  back to 2014, when macOS produced "warning: 'finite' is deprecated:
  first
  deprecated in OS X 10.9 [-Wdeprecated-declarations]. isOk =
  finite(x);".
  Patched by using isfinite() instead of finite().

[AlpsNMR](/packages/AlpsNMR)
-------

                 Changes in version 3.3.4 (2021-09-16)                  

- Fix issue with PCA plots not working as expected
- Ensure NMRExperiment names are not duplicated in a dataset (closes
#44)
- Fix issue with some title file formatting in Bruker samples (closes
#46)
- Export groups in to_ChemoSpec
- License since AlpsNMR was released has alwayd been MIT as stated in
the bioinformatics paper

[ANCOMBC](/packages/ANCOMBC)
-------

                 Changes in version 1.3.2 (2021-08-05)                  

- Add rmarkdown to the suggests field

                 Changes in version 1.3.1 (2021-08-05)                  

- Add warnings for small sample size.

[AnnotationHub](/packages/AnnotationHub)
-------------

                        Changes in version 3.1.0                        

MAJOR UPDATES

- (3.1.2) In accordance with the deprecated caching location, upgraded
  to
  error/defunct from warning/deprecated in preparaion for removal of
  dependency next release

USER-VISIBLE MODIFICATIONS

- (3.1.1) If there is a duplicate entry in the hub cache for a
  resource, hub
  code will no longer produce an ERROR. If the duplicate resource was
  not
  requested the duplicate is ignored. If the duplicate resource is
  requested,
  produce a warning for corrupt cache and continue with first found
  entry.

- (3.1.3) Fix typo in message display

- (3.1.5) Deprecate the `display,Hub-method`

BUG CORRECTION

- (3.1.7) Fix ERROR message for out-dated orgDbs

[AnVIL](/packages/AnVIL)
-----

                        Changes in version 1.6.0                        

NEW FEATURES

- (v. 1.5.5) add repository() to return the binary repository
location, if available.

- (v. 1.5.7) drs_stat() and drs_cp() support signed URLs

USER VISIBLE CHANGES

- (v. 1.5.2) drs_stat() uses multiple cores (on non-Windows) to
enhance performance

- (v. 1.5.6) install() delegates to BiocManager::install(), providing
more flexibility (e.g., installing from GitHub) and robustness.

- (v. 1.5.7) drs_stat() returns fields more selectively.


[AnVILPublish](/packages/AnVILPublish)
------------

                        Changes in version 1.4.0                        

New Features

- (v. 1.3.2) Support _bookdown.yml -- name and order vignettes

[APAlyzer](/packages/APAlyzer)
--------

                 Changes in version 1.7.3 (2021-08-01)                  

- Fixed the typo in vignettes.

                 Changes in version 1.7.2 (2021-07-31)                  

- Fixed the issues of using 3'most bam file (generate using
  ThreeMostPairBam) in PASEXP_IPA.

                 Changes in version 1.7.1 (2021-07-09)                  

- Fixed the missing SS column issues in PAS2GEF.

[aroma.light](/packages/aroma.light)
-----------

                 Changes in version 3.23.1 (2021-08-19)                 

DOCUMENTATION

- Update several citation URLs that were either broken or redirects
  elsewhere.

[artMS](/packages/artMS)
-----

                 Changes in version 1.10.2 (2021-07-13)                 

- Bug fix affecting {artmsAnalysisQuantifications()}

- Allow {artmsAnalysisQuantifications()} to process previous versions
  of artMS

- Change Extension of {artms_sessionInfo_quantification} file from
  {.txt} to {.log}

- New parameters available in {artmsQuantification()}, including:
  - Parameter {printTables}. Default {TRUE}, prints tables. FALSE
  otherwise.
  - Parameter {return_results_object}. If TRUE, it returns a list of
  data frames with MSstats results. Default is FALSE
  - If both {printTables} and {printPDF} are FALSE, then
  {return_results_object} becames TRUE and a list of data frames is
  returned

- Change default parameters of the configuration files: less verbose

                 Changes in version 1.10.1 (2021-06-30)                 

- Addressing major changes in MSstats
  - R version larger than 4.1 is now required
  - Fractions: the option "Fractions" is removed from the configuration
  file. If fractions are present, the user must include a "Fraction"
  column in the keys.txt file, which artMS will detect automatically.
  - {artmsQualityControlEvidenceBasic}: fraction parameter no longer
  required (automatically detected from the keys file)
  - {keys.txt}: use {Fraction} instead of {FractionKey}

- External packages used exclusively by the
  {artmsAnalysisQuantifications} function are not required. Those
  packages will have to be installed before running this function.

- {artmsAvgIntensityRT}: argument {species} is not longer required

- Example datasets:
  - {artms_data_ph_evidence}: the size has been significantly reduced
  (bioconductor requirement). Only two biological replicates and 1/20
  of the lines selected randomly. Including only 36 columns from the
  original evidence file
  - {artms_data_ph_msstats_results}: output from running
  {artmsQuantification} on the full version of the evidence file,
  including 4 biological replicates (instead of the reduced version
  available in the package)
  - {artms_data_ph_msstats_modelqc}: output from running
  {artmsQuantification} on the full version of the evidence file,
  including 4 biological replicates (instead of the reduced version
  available in the package)

[ASpli](/packages/ASpli)
-----

                        Changes in version 2.3.1                        

BUG FIXES

- Function vecMin inside .filterJunctionBySample threw an error if
  every element in counts vector was larger than the filter. Now works
  correctly using pmin instead.

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.17.1                        

- Fix the issue that NA is generated when no data available for
TSSEscore.

[atena](/packages/atena)
-----

                Changes in version 0.99.36 (2021-08-01)                 

USER VISIBLE CHANGES

- Submission of the first version to the Bioconductor project.

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.5.7 (2021-10-05)                  

- Fix tests for BiocParallel behaviour.

                 Changes in version 2.5.6 (2021-10-05)                  

- Revert 2.5.5 changes; add expectation to empty tests.

                 Changes in version 2.5.5 (2021-08-24)                  

- Bugfix tests with change in BiocParallel behaviour

                 Changes in version 2.5.4 (2021-08-24)                  

- Ensure ordering of genes in chains is consistent with input data
  when using divide and conquer

                 Changes in version 2.5.3 (2021-08-11)                  

- Add GeneExponent and CellExponent settings for divide and conquer.

                 Changes in version 2.5.2 (2021-08-11)                  

- Add better spacing around EFDR message.

- Fix ggplot2 `guide=FALSE` warnings

- Remove an errant browser() call

                 Changes in version 2.5.1 (2021-05-19)                  

- Add error condition for unnamed cells.

[BayesSpace](/packages/BayesSpace)
----------

                        Changes in version 1.3.1                        

Minor improvements and fixes

- Update documentation for getRDS(), mcmcChain(), and
spatialEnhance().

                        Changes in version 1.3.0                        

Minor improvements and fixes

- Added information to documentation.
- getRDS() updated with new URL.

[benchdamic](/packages/benchdamic)
----------

                 Changes in version 0.99.4 (2021-10-15)                 

- Changed stats::coef to stats4::coef

                 Changes in version 0.99.3 (2021-10-15)                 

- Added the plotIt = FALSE option for plotRMSE function

- Added some more explanations into the "Goodness of Fit" chapter

- fitNB, fitZINB, fitHURDLE, fitZIG, and fitDM functions now work with
  both phyloseq object or count matrices

                 Changes in version 0.99.2 (2021-10-11)                 

- Changed dependency from R 4.0.0 to 4.1.0

- Added example for iterative_ordering() function

- Removed the usage of @ to access object's slot in the vignette

- Removed a suppressMessages() and a suppressWarnings() from
  createTIEC()

- Added verbosity to createTIEC() function

- Added NEWS file

                 Changes in version 0.99.1 (2021-10-04)                 

- Removed unnecessary files

                 Changes in version 0.99.0 (2021-09-29)                 

- Bumped version for submission to Bioconductor

- Added README.md file

[BgeeDB](/packages/BgeeDB)
------

                       Changes in version 2.18.1                        

- Possibility to retrieve single cell full length RNA-Seq from
  Bgee 15.0 and after

- Possibility to filter data based on sex and strain for Bgee
  15.0 and after

- Possibility to filter on cellTypeId for single cell full length
  RNA-Seq for Bgee 15.0 and after

- Added pValue in the download files

[BindingSiteFinder](/packages/BindingSiteFinder)
-----------------

                       Changes in version 0.99.10                       

- coverageOverRanges() now supports mean and sum as combination method.
  Dpending
  on the returnOption, mean/ sum are computed over ranges or
  replicates.

                       Changes in version 0.99.9                        

- BSFDataSet() and BSFDataSetFromBigWig() now check the path to the
  bigwig
  files in the meta data for potential duplicates

- coverageOverRanges() now supports also ranges with different width,
  if
  returnOption = `merge_positions_keep_replicates`

- Fix bug in makeBindingSites(); The minWidth parameter is now
  implemented as
  true lower boundary (>= instead of >). The default has changed from 2
  to 3.

- Fix description in makeBindingSites(); The minCrosslinks parameter
  describes
  the number of positions covered by crosslink events, instead of the
  total
  number of crosslinks.

- Updated color scheme in rangeCoveragePlot(); and changed position of
  indicator
  box

- Updated visual of reproducibiliyCutoffPlot() function

                       Changes in version 0.99.8                        

- Updated coverageOverRange(), Function now does support different
  output
  formats, summarizing the coverage differently over range, replicates
  and
  condition

                       Changes in version 0.99.1                        

- Fix bugs for Bioconductor submission

                 Changes in version 0.99.0 (2021-05-15)                 

- Submitted to Bioconductor

[BiocCheck](/packages/BiocCheck)
---------

                        Changes in version 1.29                         

NEW FEATURES

- (1.29.10) Check for `Sys.setenv` and
  `suppressWarnings`/`suppressMessages`

- (1.29.8) Check for `sessionInfo` / `session_info` in vignette code.

- (1.29.5) Check for installation calls in vignette code.

- (1.29.1) Check for `install()` function calls in R code.

BUG FIXES

- (1.29.14) Various internal improvements to the codebase.

- (1.29.12) Checks on class membership code now include `is() ==`
  grammar.

- (1.29.6) Use appropriate input (`pkgdir`) to internal checking
  functions.

- (1.29.3) Add unit tests for legacy function searches.

- (1.29.2) rename internal function from checkIsPackageAlreadyInRepo to
  checkIsPackageNameAlreadyInUse

[BiocFileCache](/packages/BiocFileCache)
-------------

                         Changes in version 2.1                         

MAJOR UPDATES

- (2.1.1) Change caching location warning/deprecation to an ERROR in
  preparation for removal of dependency next release.

[BiocParallel](/packages/BiocParallel)
------------

                        Changes in version 1.28                         

USER VISIBLE CHANGES

- (v 1.27.3) Setting `progressbar = TRUE` for SnowParam() or
  MulticoreParam() changes the default value of `tasks` from 0 to
  `.Machine$integer.max`, so that progress on each element of `X` is
  reported.

- (v 1.27.3) `tasks` greater than `length(X)` are set to
  `length(X)`. Thus `.Machine$integer.max`, for instance, assures
  that each element of `X` is a separate task.

- (v 1.27.5) Use of random numbers is robust to the distribution
  of jobs across tasks for SerialParam(), SnowParam(), and
  MulticoreParam(), for both bplapply() and bpiterate(), using the
  RNGseed= argument to each *Param(). The change is NOT backward
  compatible -- users wishing to exactly reproduce earlier results
  should use a previous version of the package.

- (v 1.27.8) Standardize SerialParam() construct to enable setting
  additional fields. Standardize coercion of other BiocParallelParam
  types
  (e.g., SnowParam(), MulticoreParam()) to SerialParam() with
  as(., "SerialParam").

- (v. 1.27.9) By defualt, do _not_ only run garbage collection
  after every call to FUN(), except under MulticoreParam(). R's
  garbage collection algorithm only fails to do well when forked
  processes (i.e., MulticoreParam) assume that they are the only
  consumers of process memory.

- (v 1.27.11) Developer-oriented functions bploop.*() arguments
  changed.

- (v 1.27.12) Ignore set.seed() and never increment the global random
  number stream. This reverts a side-effect of behavior introduced in
  v.
  1.27.5 to behavior more consistent with version 1.26.

- (v 1.27.16) Better BPREDO support for previously started BPPARAM, and
  'transient' BPPARAM without RNGseed.

BUG FIXES

- (v 1.27.10) Typo in coercion to SerialParam when only a single worker
  specified. https://github.com/Bioconductor/BiocParallel/issues/151

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.12.0                        

NEW FEATURES

- `biocRevDepEmail` sends an email to several downstream maintainers to
  notify them of a deprecated package.

- `biocBuildEmail` allows deprecation notices using the template in the
  `inst` folder. Use `templatePath()` to see available templates by
  their
  location.

- `PackageStatus` indicates whether a package is slated for
  'Deprecation' by checking the `meat-index.dcf` file.

- `pkgDownloadStats` provides the download statistics table for a
  particular package.

BUG FIXES

- `biocDownloadStats` includes all types of Bioconductor packages

- `biocBuildReport` improved to work on old-rel, release, and devel
  Bioconductor versions

[BiocStyle](/packages/BiocStyle)
---------

                       Changes in version 2.22.0                        

USER-VISIBLE CHANGES

- Added section to HTML vignette regarding accessibility considerations
  when creating plots and figures.

BUG FIXES AND IMPROVEMENTS TO HTML STYLE

- Improved navigation for screen readers by adding role='link' and
  tabindex='0' to elements in the table of contents.  TOC navigation
  can also
  be followed by pressing "Enter" when selected in addition to clicking
  with
  the cursor.

- Addressed further styling issues for code blocks, introduced by
  changes
  in rmarkdown.  This time re-introducing padding around \<pre\> tags.

[biocthis](/packages/biocthis)
--------

                        Changes in version 1.3.8                        

SIGNIFICANT USER-VISIBLE CHANGES

- use_bioc_github_action() has been updated to match as much as
possible the changes in r-lib/actions up to the latest commit
https://github.com/r-lib/actions/commit/630f4c9d8b813f45d0327a2fc20eb264fd518450.

                        Changes in version 1.3.4                        

NEW FEATURES

- use_bioc_github_action() is now more robust in preventing tcltk
errors thanks to this pull request by Ben Laufer
https://github.com/lcolladotor/biocthis/pull/19.

                        Changes in version 1.3.2                        

NEW FEATURES

- use_bioc_github_action() now uses the AnVIL-powered package
binaries, which greatly speed up the dependency installation steps
in the docker (Linux) GitHub Actions builds. Details are available
in Nitesh Turaga's BioC2021 slides
https://github.com/nturaga/bioc2021-bioconductor-binaries.

[biocViews](/packages/biocViews)
---------

                       Changes in version 1.61.0                        

ENHANCEMENT

- (1.61.1) Added Spatial, SpatialData, SpatialWorkflow to distinguish
  from SigleCell

[biodb](/packages/biodb)
-----

                 Changes in version 1.1.16 (2021-10-19)                 

- Update documentation.

- Add ORCID for Alexis.

                 Changes in version 1.1.15 (2021-10-18)                 

- Correct getUrlContent() to handle binary files.

                 Changes in version 1.1.14 (2021-10-17)                 

- Factorize RCurl calls into global functions.

- Disable warnings in calls to readLines() when using base::url().

- Catch errors when trying to set locale.

- Correct some tests.

                 Changes in version 1.1.13 (2021-10-12)                 

- Make custom persistent the cache the default, following slowness with
  BiocFileCache in biodbHmdb.

- Remove useless bib refs in vignettes/references.bib.

- Add session info in vignettes.

- Add ORCID, URL and BugReports.

- Add an install section in main vignette.

                 Changes in version 1.1.12 (2021-10-09)                 

- Switch back to custom implementation of persistent cache, following
  errors with BiocFileCache on Windows and also slowness with HMDB.

                 Changes in version 1.1.11 (2021-10-07)                 

- Decompose test test.collapseRows() because of error on Bioconductor
  not reproduced on local computer.

                 Changes in version 1.1.10 (2021-09-30)                 

- Disable UniProt request test: it fails (result is NA) for reason
  unknown only on Bioconductor Linux server during "R CMD check". Works
  fine on local computer.

                 Changes in version 1.1.9 (2021-09-28)                  

- Correct handling of wrong URL with base::url().

                 Changes in version 1.1.8 (2021-09-28)                  

- Correct bug of UniProt request on Windows.

                 Changes in version 1.1.7 (2021-09-23)                  

- Ignore build folder when building package.

- Update documentation.

- Correct setting of R_ENVIRON_USER when building.

                 Changes in version 1.1.6 (2021-09-14)                  

- Update documentation.

                 Changes in version 1.1.5 (2021-09-13)                  

- Correct bug in return type of BiodbRequestScheduler::sendRequest().

- Correct encoding of test reference filenames.

                 Changes in version 1.1.4 (2021-09-12)                  

- Allow to set the test reference folder directly into
  runGenericTests(). This
  is now necessary for running generic tests in extension packages.

                 Changes in version 1.1.3 (2021-09-12)                  

- Set package name when calling runGenericTests() in order to find test
  ref
  files the correct way, by calling system.file().

                 Changes in version 1.1.2 (2021-09-09)                  

- Use BiocFileCache for the persistent cache system.

- Switch to R6.

- Define do...() private methods to be redefined in subclasses, instead
  of
  redefining public methods defined inside super class.

- Use now local entry files for testing parsing of entry fields.

                 Changes in version 1.1.1 (2021-06-10)                  

- Allow skipping of some fields when testing searchForEntries().

- Move test reference entries folder from tests/testthat/res to
  inst/testref.

- Move long tests folder from tests/long to longtests and enable
  Bioconductor
  long tests.

                 Changes in version 1.0.4 (2021-06-09)                  

- Bug fix: correct call to logger in BiodbPersitentCache class.

                 Changes in version 1.0.3 (2021-05-26)                  

- Bug fix: correct generic test of searchForEntries(), allowing testing
  with NA
  value.

                 Changes in version 1.0.2 (2021-05-23)                  

- Bug fix: correct return type of searchForEntries(), which now returns
  always
  a character vector and never NULL.

                 Changes in version 1.0.1 (2021-05-20)                  

- Bug fix: correct some calls to logging functions that raised a
  warning.

[biodbChebi](/packages/biodbChebi)
----------

                 Changes in version 0.99.6 (2021-10-12)                 

- Rename vignette.

- Add session info and install section in vignette.

- Use importFrom.

                 Changes in version 0.99.5 (2021-10-07)                 

- Write description on multiple lines.

- Put comment banners to indicate public and private sections in R6
  class.

                 Changes in version 0.99.4 (2021-09-28)                 

- Correct list index in vignette.

                 Changes in version 0.99.3 (2021-09-28)                 

- Complete all chapters in README.

- Remove 'foo' names in vignette.

                 Changes in version 0.99.2 (2021-09-23)                 

- Define R_BUILD_TAR in makefile to build on UNIX.

                 Changes in version 0.99.1 (2021-09-23)                 

- Upgrade maintenance files.

- Ignore build folder.

                 Changes in version 0.99.0 (2021-09-16)                 

- Submitted to Bioconductor

[biodbHmdb](/packages/biodbHmdb)
---------

                 Changes in version 0.99.4 (2021-10-18)                 

- When reading the zipped database, take the only XML file available,
  ignoring other files.

                 Changes in version 0.99.2 (2021-10-12)                 

- Add ORCID, URL and BugReports.

- Add session info and install sections in vignette.

- Corrected indentations in vignette.

                 Changes in version 0.99.1 (2021-09-29)                 

- Slight corrections for Bioconductor submission.

                 Changes in version 0.99.0 (2021-07-16)                 

- Submitted to Bioconductor

[biodbKegg](/packages/biodbKegg)
---------

                 Changes in version 0.99.4 (2021-10-13)                 

- Merge both vignettes into a single one.

- Remove messages from magick package.

                 Changes in version 0.99.3 (2021-10-12)                 

- Add install section in vignette.

- Use importFrom.

- Add ORCID, URL and BugReports.

                 Changes in version 0.99.2 (2021-10-12)                 

- Rename main vignette file.

- Add reference.

- Add session info in vignettes.

                 Changes in version 0.99.1 (2021-09-28)                 

- Made some corrections for submitting to Bioconductor.

                 Changes in version 0.99.0 (2021-09-16)                 

- Submitted to Bioconductor

[biodbLipidmaps](/packages/biodbLipidmaps)
--------------

                 Changes in version 0.99.3 (2021-10-18)                 

- Remove deprecated slow frequency of request send.

                 Changes in version 0.99.2 (2021-10-12)                 

- Add ORCID, URL and BugReports.

- Add citation into vignette.

- Add install and session info sections to vignette.

                 Changes in version 0.99.1 (2021-09-29)                 

- Complete README.

- Slight corrections for Bioconductor submission.

- Add generated documentation files.

                 Changes in version 0.99.0 (2021-09-16)                 

- Submitted to Bioconductor.

[biodbUniprot](/packages/biodbUniprot)
------------

                 Changes in version 0.99.4 (2021-10-12)                 

- Add session info and install section in vignette.

- Rename vignette.

- Use importFrom.

- Add ORCID, URL and BugReports.

                 Changes in version 0.99.3 (2021-10-11)                 

- Show progress when filtering results in geneSymbolToUniprotIds().

                 Changes in version 0.99.2 (2021-10-03)                 

- Use biodb 1.1.10.

                 Changes in version 0.99.1 (2021-09-23)                 

- Ignore build folder when building package.

- Add documentation.

                 Changes in version 0.99.0 (2021-09-16)                 

- Submitted to Bioconductor

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.50.0                        

MINOR CHANGES

- useMart() and listMarts() will warn users if using http to access
  Ensembl.  https will be enforced by Ensembl from late 2021.

BUG FIXES

- Address issue where checking the list of Ensembl Archives would stop
  all queries from working if the main www.ensembl.org site was
  unavailable.

- Fix bug introduced in getSequence() where asking for flanking
  sequences
  resulted in an invalid query.

- The argument 'host' is no longer ignored in useEnsembl() (Thanks to
  forum
  user "A" - https://support.bioconductor.org/p/9139019/)

[BioPlex](/packages/BioPlex)
-------

                        Changes in version 1.0.0                        

- Initial release of the BioPlex package

[biotmle](/packages/biotmle)
-------

                       Changes in version 1.17.0                        

- Removal of `future` and `doFuture` for simplification of
  parallelization. All
  control of parallel computation now done through `BiocParallel`.

[bugsigdbr](/packages/bugsigdbr)
---------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- Created package.

[BUSseq](/packages/BUSseq)
------

                Changes in version 0.99.18 (2020-06-01)                 

- Importing "assay<-" in the NAMESPACE

                Changes in version 0.99.17 (2020-06-01)                 

- Correcting the import of dependencies

                Changes in version 0.99.16 (2020-06-01)                 

- Avoiding the overlapping when loading dependencies

                Changes in version 0.99.15 (2020-06-01)                 

- Removing the redundant print and summary function

                Changes in version 0.99.14 (2020-06-01)                 

- Modifying the output of BUSseq_MCMC as a SingleCellExperiment object
  and allowing the input as that object as well

                Changes in version 0.99.13 (2020-05-15)                 

- Modifying the running example

                Changes in version 0.99.12 (2020-05-15)                 

- Revising the random number generator for MacOS

                Changes in version 0.99.11 (2020-05-15)                 

- Debugging for MacOS

                Changes in version 0.99.10 (2020-05-15)                 

- Continuing to debug for MacOS

                 Changes in version 0.99.9 (2020-05-15)                 

- Modifying the warnings when building on MacOS

                 Changes in version 0.99.8 (2020-05-15)                 

- Correcting the usage of OS macro

                 Changes in version 0.99.7 (2020-05-14)                 

- Updating th source code for MacOS

                 Changes in version 0.99.6 (2020-05-13)                 

- Downsizing the example dataset to avoid reaching the time limit of
  check

                 Changes in version 0.99.5 (2020-05-12)                 

- Adding the macro for Mac in the source code

- Shortening the vignettes

- Adding Watched Tags BUSseq to my profile

                 Changes in version 0.99.4 (2020-05-12)                 

- Correct R version dependency

- Solve the warnings by missing parenthesees

                 Changes in version 0.99.3 (2020-05-08)                 

- Set LazyData as TRUE for building vignette

                 Changes in version 0.99.2 (2020-05-08)                 

- Correctly pushing by Git

                 Changes in version 0.99.1 (2020-05-07)                 

- Revised the warnings by R CMD check, including adding parentheses in
  src and R dependency in DESCRIPTION

                 Changes in version 0.99.0 (2020-05-04)                 

- Submitted to Bioconductor

[cageminer](/packages/cageminer)
---------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[CAGEr](/packages/CAGEr)
-----

                        Changes in version 2.0.0                        

BACKWARDS-INCOMPATIBLE CHANGES

- The `CAGEset` class is removed.

- Accessors using plain `data.frame` formats are removed.

- Modifier functions return the object instead of silently modifying it
  in the global environment.  Thus, you can use R pipes (`|>`).

- Removed export function that unconditionally wrote files to the
  working directory, such as `exportCTSStoBedGraph`.  As a replacement
  a
  new converter is provided, `exportToBrowserTrack`, that produces
  `UCSCData` objects that the user can then wite to files using the
  _rtracklayer_ package.

- Removed the `extractExpressionClass` function, as the information is
  easily accessible from within the `CTSS` or `ConsensusClusters`
  objects.

NEW FEATURES

- `CTSSnormalizedTpmGR` and `CTSStagCountGR` now accept `"all"` as a
  special
  sample name to return them all in a `GRangesList` object.

- Plot interquantile width and expression clusters with ggplot2.

BUG FIXES

- Corrected the `getExpressionProfiles` function to accept CAGEexp
  objects.

- Updated the `exampleCAGEexp` object to contain expression classes.

- Restore paraclu support for CAGEexp objects.

- Corrected a bug in `aggregateTagClusters` that was causing
  mislabelling
  of tag clusters (PR#42).

- Prevent `plotReverseCumulatives` from crashing when values are not in
  range. (PR#43).

[Cardinal](/packages/Cardinal)
--------

                       Changes in version 2.11.3                        

BUG FIXES

- Fix strange behavior from random number generation in R >= 4.1.1

                       Changes in version 2.11.2                        

BUG FIXES

- Fix reference naming scheme for binning and alignment methods

                       Changes in version 2.11.1                        

BUG FIXES

- Use as(x, 'DFrame') instead of as(x, 'DataFrame')

- Fix logical length > 1 error in 'segmentationTest()'

[cbaf](/packages/cbaf)
----

                 Changes in version 1.16.0 (2021-10-14)                 

New Features

- Terms updated.

[cBioPortalData](/packages/cBioPortalData)
--------------

                        Changes in version 2.6.0                        

New features

- A study's build status can be obtained from getStudies(), which has
replaced data('studiesTable').
- Partial loading of data files supported. A warning is emitted when
a
data file is not able to be loaded in cBioDataPack.
- cBioPortalData checks the data(studiesTable) to verify that study
datasets are building, otherwise provide a message in interactive
sessions.

[celda](/packages/celda)
-----

                 Changes in version 1.9.3 (2021-10-04)                  

- Fixed bug in checking background matrix with decontX
- Switched to using Github Actions for Continuous Integration
- Fixed plotting bugs in celda results reports
- Speed up final step in decontX when creating final decontaminated
matrix

                 Changes in version 1.9.2 (2021-07-19)                  

- Added a NEWS.md file to track changes to the package.
- Added new tutorials and documentation generated with pkgdown.
- Removed warnings in plotRPC functions.
- Added use of "displayName" to several functions that show feature
names.
- Minor bug fix when the input matrix was sparse and contained
non-integer values.
- Several improvements to plotting functions.

[Cepo](/packages/Cepo)
----

                        Changes in version 0.0.1                        

- Added a NEWS.md file to track changes to the package.

[CeTF](/packages/CeTF)
----

                        Changes in version 1.5.4                        

- Update CITATION

                        Changes in version 1.4.3                        

- Fix netConditionsPlot function bug related to factors in data.frame

[cfDNAPro](/packages/cfDNAPro)
--------

                       Changes in version 0.99.1                        

- Remove bugs in plotSingleGroup.R
- Documentation improvements.

                       Changes in version 0.99.0                        

- Now cfDNAPro supports bam file as input for data characterisation.
- Coding style improvements.
- Documentation improvements.
- Submitted to Bioconductor.

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.27.7                        

- Fix the error "!anyNA(m32) is not TRUE" in seqlevelsStyle is not
  handled.

                       Changes in version 3.27.6                        

- Fix the error "pvalue" undefined columns selected in enrichmentPlot

                       Changes in version 3.27.5                        

- update documentation of pipeline.rmd

                       Changes in version 3.27.4                        

- use formatSeqnames function to handle the error from seqlevelsStyle:
  "cannot switch some of GRCm38's seqlevels from NCBI to UCSC style"

                       Changes in version 3.27.3                        

- use formatSeqnames function to handle the error from seqlevelsStyle:
  "!anyNA(m31) is not TRUE "

                       Changes in version 3.27.2                        

- add keepExonsInGenesOnly for genomicElementDistribution

- add upstream and downstream for assignChromosomeRegion function when
  define promoter and downstream regions.

                       Changes in version 3.27.1                        

- Add the possibility of find overlaps by percentage covered of
  interval
  for function findOverlapsOfPeaks

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.29.2                        

- extend functions for plotting peak profiles to support other types
of bioregions (2021-10-15, Fri, @MingLi-292, #156, #160, #162, #163)

                       Changes in version 1.29.1                        

- add example for seq2gene function (2021-05-21, Fri)

[CIMICE](/packages/CIMICE)
------

                        Changes in version 1.1.3                        

Overview:

- I/O Improvements

New functionalities:

- Fixed and improved input in CAPRIpop dataset format

                        Changes in version 1.1.1                        

Overview:

- I/O Improvements

New functionalities:

- Fixed bugs related with byrow options when reading CAPRI
  formatted datasets

- Improved VisNetwork output to include information about samples
  inside nodes

[ClassifyR](/packages/ClassifyR)
---------

                       Changes in version 2.14.0                        

- 
  Upsampling and downsampling to equalise class sizes added.

[cliProfiler](/packages/cliProfiler)
-----------

                 Changes in version 0.99.0 (2021-01-22)                 

- Submitted to Bioconductor

[cliqueMS](/packages/cliqueMS)
--------

                        Changes in version 1.7.0                        

- Bug fixes: Change in example of computeCliques function.

[clonotypeR](/packages/clonotypeR)
----------

                       Changes in version 1.32.0                        

DEPRECATION NOTICE

- Note that clonotypeR is depreacted.  Please adopting it or using a
  different package.

[clusterProfiler](/packages/clusterProfiler)
---------------

                        Changes in version 4.1.4                        

- import yulab.utils (2021-08-20, Fri)

                        Changes in version 4.1.3                        

- Remove Human Gut Microbiome dataset as the functionalities are
provided in https://github.com/YuLab-SMU/MicrobiomeProfiler
(2021-08-15, Sun)

                        Changes in version 4.1.2                        

- update citation and DESCRIPTION (2021-08-15, Sun)
- update kegg_species.rda and allow online download using KEGG api
(2021-08-14, Sat)

                        Changes in version 4.1.1                        

- add citation (new paper published on The Innovation) (2021-07-04,
Sun)

[clustifyr](/packages/clustifyr)
---------

                 Changes in version 1.5.2 (2021-10-04)                  

- `clustify_lists()` support for output of overlapping genes
  (`details_out = TRUE`)

- Added truncated mean and trimean modes to `average_clusters()`

                 Changes in version 1.5.1 (2021-08-04)                  

- `clustify_lists()` support for uneven number of markers

- Deprecated SeuratV2 support

[CNVfilteR](/packages/CNVfilteR)
---------

                        Changes in version 1.7.4                        

MINOR

- Minor fix: removed ` from Vignette

                        Changes in version 1.7.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- CNVfilteR specfically supports VCFs produced by Torrent Variant
  Caller

                        Changes in version 1.7.2                        

BUG FIXES

- Fixed: parsing some VCFs produced an unexpected error

MINOR

- Minor bug fixed: stop message was not completely shown on
  loadSNPsFromVCF()

                        Changes in version 1.7.1                        

MINOR

- CITATION udpated

[cola](/packages/cola)
----

                        Changes in version 1.9.1                        

- add `uniquely_high_in_one_group` method in `get_signatures()`.

- add `compare_partitions()`.

- parallel computing is implemented with foreach + doParallel

                        Changes in version 1.9.0                        

- use row/column* family functions in `adjust_matrix()` to reduce the
  memory
  usage as well as improve the speed.

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

                        Changes in version 2.9.4                        

- fixed a bug of missing right annotation legends for vertically
  concatenated heatmaps.

- `Legend()`: support `border` to be set to `asis`.

- Rasterization: the default maximal size for temporary image is set to
  30000 px (both for width and height).

- add a new argument `beside` in `anno_barplot()` to position bars
  beside each other.

- add `plot()` method for `Heatmap` and `HeatmapList` classes.

- add `anno_customize()`.

                        Changes in version 2.9.3                        

- `pheatmap()`/`heatmap()`/`heatmap.2()`: set default of run_draw to
  FALSE.

- throw error when the heatmaps (list) are already initialized by
  draw() when adding them.

- set `wrap = TRUE` in `grid.grabExpr()` when capturing the legend
  objects.

- `make_comb_mat()`: support `GRangesList` object as input.

- legends: fixed a bug of the grid heights were not correctedly
  calculated.

- discrete annotations: neighbour grids are merged into one single grid
  if they have the
  same values.

- `anno_barplot()`: allows to add numbers on top of bars.

- `UpSet()`: axis labels are automatically formated for genomic
  coordinates.

- `AnnotationFunction()`: add a new argument `cell_fun`.

- When the dendrogram height is zero, the corresponding viewport has
  scale (0, 1).

                        Changes in version 2.9.2                        

- fixed a bug of `bg_col` for transposed matrix in `UpSet()`.

- print warnings if names of annotations have different orders from the
  matrix row/column names.

                        Changes in version 2.9.1                        

- fixed a bug of editing gTree object where the list element "just" has
  been
  changed to "justification" in recent R versions.

[conclus](/packages/conclus)
-------

                Changes in version 1.1.002 (2021-08-31)                 

Made the following significant changes

  -  Added an internal function .retrieveClustersNumberK to suggest the clusters number to use in lusterCellsInternal().

  - Added suggestedClustersNumber slot in scRNAseq class to retrieve
this suggested clusters number.

  - Added accessors getSuggestedClustersNumber and
setSuggestedClustersNumber.

  - Updated all the objects used in examples and tests to have this
slot.

[CoreGx](/packages/CoreGx)
------

                        Changes in version 1.5.7                        

- Add TreatmentResponseExperiment class, a simple wrapper around
LongTable to make the class syntax more domain specific
- Add CoreSet2 structure to support creation of CoreSets with the
modified class structure introducted in BioC 3.13
- CoreSets can now be made with treatment combination experiments via
the TreatmentResponseExperiment class!

                        Changes in version 1.5.6                        

- Fix bug in LongTable -> data.table coerce method that was causing
rows of some assays to be dropped (closes issue #)

                        Changes in version 1.5.5                        

- Fix bug in .distancePointLine where function fails with no
intercept
specified (Issue #120)
- Added support for aggregating an assay inside of a LongTable class
object
- Some in-progress updates to the CoreSet constructor which will be
completed for the Fall release
- Fixed an error in treatmentNames example
- Fixed roxygen2 documentation warnings about S4 method documentation
- Overhauled LongTable coerce methods to use the LongTableDataMapper
class instead of the deprecated 'LongTable.config' attribute

                        Changes in version 1.5.4                        

- Fix bug in $<- and [[<- methods where value was returned instead of
updated object
- Fix bug in .sanitize input caused by length > 1 coercing to logical
vector

                        Changes in version 1.5.3                        

- Fix bug in connectivityScore caused by length > 1 coercing to
logical vector; this should fix errors in RadioGx and PharmacoGx
vignettes that were caused by failed R CMD build

                        Changes in version 1.5.2                        

- Add subsetBySample method for CoreSet object; this is the first
step
in modularizing the subset methods for reuse in dependent packages
- Added a CoreSet-utils documentation section to document subset,
intersect, combine and other set operations for a CoreSet object.

                        Changes in version 1.5.1                        

- Fixed some spelling errors and incorrect code chunk configurations
in the LongTable vignette
- Fix bug in .rebuildProfiles where the function fails if
replicate_id
is assigned as a rowID column in the LongTable in @sensitivity

                        Changes in version 1.5.0                        

- Bioconductor spring 2021 release
- Added the DataMapper abstract class
- Added the LongTableDataMapper concrete class
- Added the metaConstruct method, for making an S4 object from a
sub-class of DataMapper
- Updated LongTable vignette with documentation for the DataMapper
and
LongTableDataMapper
- Refactored various methods to work with a LongTable in
@sensititivty
- Refactored various methods to work with a MultiAssayExperiment in
@molecularProfiles

[csdR](/packages/csdR)
----

                       Changes in version 0.99.6                        

- Bugfix: Open MP was not working correctly because of missing
compiler flags. For this reason, the Makevars file has been created.
- Calculation of column ranks now uses matrixStats::colRanks instead
of an apply statement with base::rank.

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package. This is the
first public version of the package.

[cTRAP](/packages/cTRAP)
-----

                       Changes in version 1.12.0                        

Web server support (optimised to run in ShinyProxy)

- cTRAP(): new global interface with all cTRAP functionality in one
place
- Sessions can be created and loaded via a token or a RDS file
- Session data is automatically saved in a RDS file to a folder in
the working directory named based on the current session token
- Long-running tasks can be performed in the background using the
Celery task manager via Flower's REST API and their output is
automatically loaded in the corresponding session
- Loading icon in navigation menu when Shiny is busy
- Use the faster and efficient file format from R package qs instead
of RDS:
- Faster download and loading of pre-processed remote files
(compound molecular descriptors and gene expression and drug
sensitivity associations)

Bug fixes and minor improvements

- convertGeneIdentifiers() replaces convertENSEMBLtoGeneSymbols():
- Use AnnotationHub to convert to gene symbols (instead of biomaRt
that has been unstable)
- loadENCODEsamples():
- New argument to select folder where to download data
- analyseDrugSetEnrichment():
- Cross-match more compounds between datasets by discarding
non-alphanumeric characters and ignoring case
- Fix incorrect columns used for each dataset when merging
datasets
- Visual interface:
- Fix crash when plotting dataset comparison using values with too
many zeroes for density estimation
- Add progress bars for slower tasks
- Fix crash when using shiny 1.7.0 (avoid malformed, custom UI
elements)
- Drug set enrichment analysis interface:
- Show all drug sets available to (down)load
- Show loading indicator when loading different drug sets
- Hide "leading edge" column of the results by default


[cytoKernel](/packages/cytoKernel)
----------

                       Changes in version 0.99.0                        

- Submitted marr 0.99.0 to Bioconductor on October 2, 2020.
- Added a NEWS.md file to track changes to the package.

[cytomapper](/packages/cytomapper)
----------

                 Changes in version 1.5.4 (2021-09-17)                  

- It is not required anymore to specify exactly the right colours

                 Changes in version 1.5.3 (2021-09-16)                  

- Added option to read in .h5 files

                 Changes in version 1.5.2 (2021-09-15)                  

- Added description on how to handle images with couplet/patchwork

                 Changes in version 1.5.1 (2021-05-19)                  

- Bugfix: erroneous dimension setting when legend=NULL

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

[DaMiRseq](/packages/DaMiRseq)
--------

                        Changes in version 2.6.0                        

- Major: We fixed a bug in DaMiR.ModelSelect. Now optimal models
  are correctly selected;

- Major: Now users can plot specific graphs in DaMiR.Allplot and
  we added new plots;

- Minor: We modified the color scale in corrplot

[dasper](/packages/dasper)
------

                        Changes in version 1.3.2                        

NEW FEATURES

- Fix bugs within plot_sashimi() and enable the visualization of raw
junction counts.

[dearseq](/packages/dearseq)
-------

                 Changes in version 1.5.1 (2021-09-01)                  

- switching from CompQuadForm::davies() to the "saddlepoint" method
from survey::pchisqsum() for computing quadratic form asymptotic
p-values

[decoupleR](/packages/decoupleR)
---------

                        Changes in version 2.0.0                        

Changes

- Some method's names have been changed to make them easier to
identify:

- pscira now is called Weighted Sum (wsum).
- mean now is called Weighted Mean (wmean).
- scira now is called Univariate Linear Model (ulm).

- The column name for tf in the output tibbles has been changed to
source.

- Updated documentation for all methods.

- Updated vignette and README.

- decouple function now accepts order mismatch between the list of
methods and the list of methods's arguments.

- Moved benchmark branch to a separate repository as its own package:
https://github.com/saezlab/decoupleRBench

New features

- New methods added:

- Fast Gene Set Enrichment Analysis (fgsea).
- AUCell.
- Univariate Decision Tree (udt).
- Multivariate Decision Tree (mdt).
- Multivariate Linear Model (mlm).

- New decoupleR manuscript repository:
https://github.com/saezlab/decoupleR_manuscript

- New consensus score based on RobustRankAggreg::aggregateRanks()
added when running decouple with multiple methods.

- New statistic corr_wmean inside wmean.

- Methods based on permutations or statistical tests now return also
a
p-value for the obtained score (fgsea, mlm, ora, ulm, viper, wmean
and wsum).

- New error added when network edges are duplicated.

- New error added when the input matrix contains NAs or Infs.

                        Changes in version 1.1.0                        

New features

All new features allow for tidy selection. Making it easier to
evaluate
different types of data for the same method. For instance, you can
specify the columns to use as strings, integer position, symbol or
expression.

Methods

- New decouple() integrates the various member functions of the
decoupleR statistics for centralized evaluation.

- New family decoupleR statists for shared documentation is made up
of:

- New run_gsva() incorporate a convinient wrapper for
GSVA::gsva().
- New run_mean() calculates both the unnormalized regulatory
activity and the normalized (i.e. z-score) one based on an
empirical distribution.
- New run_ora() fisher exact test to calculate the regulatory
activity.
- New run_pscira() uses a logic equivalent to run_mean() with the
difference that it does not accept a column of likelihood.
- New run_scira() calculates the regulatory activity through the
coefficient $\beta_1$ of an adjusted linear model.
- New run_viper() incorporate a convinient wrapper for
viper::viper().

Converters

- New functions family convert_to_ variants that allows the
conversion
of data to a standard format.
- New convert_to_() return the entry without modification.
- New convert_to_gsva() return a list of regulons suitable for
GSVA::gsva().
- New convert_to_mean() return a tibble with four columns: tf,
target, mor and likelihood.
- New convert_to_ora() returns a named list of regulons; tf with
associated targets.
- New convert_to_pscira() returns a tibble with three columns: tf,
target and mor.
- New convert_to_scira() returns a tibble with three columns: tf,
target and mor.
- New convert_to_viper() return a list of regulons suitable for
viper::viper()

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

[DegNorm](/packages/DegNorm)
-------

                        Changes in version 1.3.1                        

- plot_coverage function has an option "samples" added to allow user to
  visualize a subset of samples.

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.20.0                        

BUG FIXES

- Fix long-standing bugs in dense2sparse():
  - mishandling of NAs/NaNs in input
  - 1D case didn't work

[DelayedTensor](/packages/DelayedTensor)
-------------

                        Changes in version 1.0.0                        

- Package released

[DepecheR](/packages/DepecheR)
--------

                 Changes in version 1.9.1 (2021-10-16)                  

- depecheCoFunction updated to remove traces of the old dAllocate
  functionality

[DExMA](/packages/DExMA)
-----

                        Changes in version 1.1.2                        

- Weights added to the Stouffer's method

[DIAlignR](/packages/DIAlignR)
--------

                        Changes in version 2.1.4                        

- Using non-linear fit to constrain similarity matrix.

- Recalculate intensity for given peaks.

- Direct alignment to root is added.

- Multiple strategies for getting distance matrix and
  tree-agglomeration is added.

- Remove peaks if not aligned.

                        Changes in version 2.1.0                        

- Added support for IPF based Post-translation Modification alignment.

- Added Minimum Spanning Tree based alignment for reference-free
  approach.

- Reintegrate area under the peak for low-confidence peaks.

- Supporting Savitzky-Golay smoothing while calculating area.

- Added input to select peptides which to align.

- Speed improvement in progressive alignment by storing new featues in
  lists instead of data frames.

[DiffBind](/packages/DiffBind)
--------

                        Changes in version 3.3.4                        

- Fix bFlip issues

- Fix bug where mode not passed to summarizeOverlaps

- Add $inter.feature config parameter

                        Changes in version 3.3.2                        

- Re-compute FDR when fold change other than 0 is specified

- Remove most gc() calls for performance

- Roll in bugfixes

[Dino](/packages/Dino)
----

                       Changes in version 0.99.5                        

- Data reformatting

                       Changes in version 0.99.4                        

- Revised formatting for Bioconductor submission.
- Bug fix

                       Changes in version 0.99.3                        

- Reduce size of tutorial data

                      Changes in version 0.99.1                        

- Updated formatting for Bioconductor submission.
- Reduced size of sample dataset
- Updated default number of cores to 2

                       Changes in version 0.99.0                        

- Formatted to submit to Bioconductor.
- Add SingleCellExperiment functionality.
- Update vignette with SingleCellExperiment example.

[distinct](/packages/distinct)
--------

                        Changes in version 1.4.1                        

- fixed bug in log2_FC.R;

- added rmarkdown as `Suggests` in DESCRIPTION.

[dittoSeq](/packages/dittoSeq)
--------

                         Changes in version 1.6                         

- Vignette Update: Added a 'Quick-Reference: Seurat<=>dittoSeq'
section.
- Build & Test Infrastructure Update: Removed Seurat dependency from
all build and test materials by removing Seurat code from the
vignette and making all unit-testing of Seurat interactions
conditional on both presence of Seurat and successful SCE to Seurat
cnversion.
- Bug Fixes: 

  1. Fixed dittoFreqPlot calculation machinery to properly
target all cell types but only necessary groupings for every sample.
Removed the 'retain.factor.levels' input because proper calculations
treat 'var'-data as a factor, and groupings data as non-factor. 
  2. Allowed dittoHeatmap() to properly 'drop_levels' of annotations by
ensuring 'annotation_colors' is not populated with colors for empty
levels which would be dropped. 3- Made 'do.label' machinery of
scatter plots robust to NAs.

[DOSE](/packages/DOSE)
----

                       Changes in version 3.19.4                        

- update clusterProfiler citation (2021-09-30, Thu)
- upate error message of enricher_internal (2021-9-3, Fri)

                       Changes in version 3.19.3                        

- upate DisGeNET and NCG data (2021-8-16, Mon)

                       Changes in version 3.19.2                        

- bug fixed, change 'is.na(path2name)' to 'all(is.na(path2name))'
(2021-06-21, Mon)

                       Changes in version 3.19.1                        

- add dr slot to compareClusterResult, enrichRestul and
gseaResult(2021-5-21, Fri)

[dStruct](/packages/dStruct)
-------

                 Changes in version 0.99.3 (2021-05-23)                 

- Ensured that all global variables are well-defined in the namespace.

                 Changes in version 0.99.2 (2021-05-22)                 

- Revised to address comments by the Bioconductor reviewer.

- dStruct now uses the IRanges object from Bioconductor.

- All function names follow camel case.

- More descriptions of functions.

- Added a test to check validity of code when running dStruct in the
  proximity_assisted mode.

                 Changes in version 0.99.1 (2021-04-11)                 

- Fixed errors and warnings from checks by bioc-issue-bot.

                 Changes in version 0.99.0 (2021-04-07)                 

- Submitted to Bioconductor

[easier](/packages/easier)
------

                        Changes in version 0.9.0                        

- Added a NEWS.md file to track changes to the package.

[eisaR](/packages/eisaR)
-----

                        Changes in version 1.5.2                        

- Add options to collapse introns by gene and restrict introns to
feature ranges in getFeatureRanges.

[EnhancedVolcano](/packages/EnhancedVolcano)
---------------

                        Changes in version 1.12                         

- added max.overlaps and min.segment.length to provide further control
  over
  connectors. max.overlaps replaces maxoverlapsConnectors, but both can
  still
  be used for legacy purposes

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.13.2                        

- mv ep_str_wrap to yulab.utils::str_wrap (2021-10-13, Wed)
- adjust the order of legends for dotplot, emapplot, cnetplot and
treeplot(2021-10-8, Fri)
- update treeplot: add "dotplot" and "heatmap" panels for
treeplot(2021-9-15, Wed)
- update dotplot: enable size parameter applicable to other columns
of
compareClusterResult(2021-9-17, Fri)
- enable label_format parameter for heatplot (2021-09-01, Wed)
- add get_ggrepel_segsize function to set segment.size value for
ggrepel(2021-08-29, Sun)
- update ep_str_wrap (2021-08-28, Sat)
- cnetplot now works with a named list (2021-08-23, Mon;
clusterProfiler#362)

                       Changes in version 1.13.1                        

- use aplot::plot_list instead of cowplot::plot_grid (2021-06-13, Sun
- add color_category and color_gene parameters for
cnetplot(2021-6-11,
Fri)
- Enables showCategory parameter to support character input in
dotplot.compareClusterResult(2021-6-10, Thu)

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.17.4                        

- Fix issue with extracting 5' or 3' UTRs for transcript without UTRs.

                       Changes in version 2.17.3                        

- Make parameter `port` optional in the script to create EnsDb
  databases.

                       Changes in version 2.17.2                        

- Disable ideogram plotting in vignettes.

                       Changes in version 2.17.1                        

- Fix error when importing uncompressed GTF files.

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.1.9 (2021-09-19)                  

- very fast end memory-efficient BAM loading using HTSlib
  - for now reads paired-end BAM only

- min.baseq to reduce the effect of sequencing errors

- very fast Fisher Exact from HTSlib

- old code removed

                 Changes in version 1.1.0 (2021-05-21)                  

- released at bioconductor

[epigraHMM](/packages/epigraHMM)
---------

                        Changes in version 1.0.8                        

- Minor bug fix in maxStepProb documentations

                        Changes in version 1.0.7                        

- Exporting maxStepProb, which compute the MLE of initial and
transition probabilities of a K-state HMM, as well as
simulateMarkovChain, which simulates a Markov chain of length 'n'
given a matrix of transition probabilities

                        Changes in version 1.0.6                        

- Minor bug fix in controlEM documentation

                        Changes in version 1.0.5                        

- Minor bug fix in callPatterns and info function (explict import of
S4Vectors::mcols and utils::tail).

- Exporting expStep function, which implements the E-step of EM
algorithm (forward-backward & Viterbi algorithm) for a K-state HMM.

                        Changes in version 1.0.4                        

- Adding function callPatterns to exp[ort] combinatorial patterns (or
posterior probabilities) associated with a given set of genomic
regions.

- Adding function info to print summary statistics from epigraHMM
output. This function will print the model's BIC, log-likelihood,
and combinatorial patterns associated with mixture model components.

- Adding new example dataset helas3 with ENCODE ChIP-seq data from
broad epigenomic marks H3K27me3, H3K36me3, and EZH2.

- Adding option to prune combinatorial patterns associated with rare
states. See vignette for details.

- In differential peak calling, epigraHMM now exports combinatorial
pattern table. See vignette for details.

- Improvement of the vignette to clarify epigraHMM's use of
blacklisted regions and gap tracks.

                        Changes in version 1.0.3                        

- Minor updates in the NEWS file as well as the README page.

                        Changes in version 1.0.2                        

- epigraHMM now exports a function called segmentGenome that segments
a given genome (e.g. 'mm10') into non-overlapping genomic windows
while considering gap tracks and blacklisted regions.

                        Changes in version 1.0.1                        

- Minor fix in the package DESCRIPTION file and version numbers


[escape](/packages/escape)
------

                        Changes in version 1.3.1                        

- Aligning versions to the current bioconductor release

- Added DietSeurat() call in vignette to prevent issues

[esetVis](/packages/esetVis)
-------

                       Changes in version 1.19.1                        

- default for typePlot: fix issue length > 1 in coercion to logical

- fix for ggplot2 >= 3.3.4: replace guides(fill = FALSE) by guides(fill
  = 'none')

- fix few notes check SummarizedExperiment + ggvis

[ExperimentHub](/packages/ExperimentHub)
-------------

                        Changes in version 2.1.0                        

MAJOR UPDATES

- (2.1.1) In accordance with the deprecated caching location, upgraded
  to
  error/defunct from warning/deprecated in preparaion for removal of
  dependency next release

[FamAgg](/packages/FamAgg)
------

                       Changes in version 1.21.2                        

- Allow writing of HaploPainter input files without a HaploPainter
  installation.

                       Changes in version 1.21.1                        

- Add information on kinship-based relatedness to kinship sum test
  results.

- Keep memory consumption constant. This allows arbitrary long
  simulation
  runs without running out of memory. Fixes issue #22. Due to the
  solution
  of that problem, histograms and densities reported by the simulation
  functions may slightly deviate from comparable former runs on the
  same
  data.

[famat](/packages/famat)
-----

                 Changes in version 1.3.1 (2021-10-13)                  

- Fix tests

                 Changes in version 1.3.0 (2020-11-27)                  

- Rshiny modifications for figure in paper

[fCCAC](/packages/fCCAC)
-----

                       Changes in version 1.19.7                        

- Adds html vignette.

[fgsea](/packages/fgsea)
-----

                       Changes in version 1.19.4                        

- plotGseaTable now accepts units vector for column widths

                       Changes in version 1.19.2                        

- Fixed fora() failing to run on a single pathway

- Fixed problems random gene set generation for large k (issue #94)

- Changed default eps to 1e-50

[FindIT2](/packages/FindIT2)
-------

                Changes in version 0.99.13 (2021-08-11)                 

- convert geom_density into geom_hist in plot_annoDistance function

                Changes in version 0.99.11 (2021-07-16)                 

- add findIT_enrichWilcox function
- delete findIT_enrichInShuffle function
- rename findIT_enrichInAll to findIT_enrichFisher

                Changes in version 0.99.10 (2021-07-15)                 

- move all shiny function in FindIT2 into InteractiveFindIT2

                 Changes in version 0.99.0 (2021-06-27)                 

- Submitted to Bioconductor

[fishpond](/packages/fishpond)
--------

                        Changes in version 2.0.0                        

- New loadFry() function, written by Dongze He with
  contributions from Steve Lianoglou and Wes Wilson.
  loadFry() helps users to import and process
  alevin-fry quantification results. Can process
  spliced, unspliced and ambiguous counts separately
  and flexibly. Has specific output formats designed
  for downstream use with scVelo or velocity analysis.
  See ?loadFry for more details.

- Adding correlation tests: Spearman or Pearson
  correlations of a numeric covariate with the
  log counts, or with the log fold changes across
  pairs. The Spearman correlation test with counts
  was already implemented in the original SAMseq
  method as response type = "Quantitative".
  For new functionality see 'cor' argument in the
  ?swish man page.

- Adding importAllelicCounts() to facilitate importing
  Salmon quantification data against a diploid
  transcriptome. Can import either as a 'wide'
  format or as 'assays'. Leverages tximeta().
  For gene-level summarization, importAllelicCounts()
  can create an appropriate tx2gene table
  with the necessary a1 and a2 suffices,
  and it will automatically set txOut=FALSE, see
  ?importAllelicCounts for more details.

- Added a 'q' argument to plotInfReps to change the
  intervals when making point and line plots.

- Switched the legend of plotInfReps so that
  reference levels will now be on the bottom,
  and non-reference (e.g. treatment) on top.

                       Changes in version 1.99.18                       

- Added helper functionality to importAllelicCounts,
  so it will create an appropriate tx2gene table
  with the necessary a1 and a2 suffices,
  and it will automatically set txOut=FALSE.

- Added a 'q' argument to plotInfReps to change the
  intervals when making point and line plots.

- Switched the legend of plotInfReps so that
  reference levels will now be on the bottom,
  and non-reference (e.g. treatment) on top.

- Added loadFry() to process alevin-fry
  quantification result. Can process spliced,
  unspliced and ambiguous counts separately
  and flexibly.

                       Changes in version 1.99.15                       

- Adding correlation tests: Spearman or Pearson
  correlations of a numeric covariate with the
  log counts, or with the log fold changes across
  pairs. The Spearman correlation test with counts
  was already implemented in the original SAMseq
  method as response type = "Quantitative".
  For new functionality see 'cor' argument in the
  ?swish man page.

- Adding importAllelicCounts() to facilitate importing
  Salmon quantification data against a diploid
  transcriptome. Can import either as a 'wide'
  format or as 'assays'. Leverages tximeta().

                        Changes in version 1.9.6                        

- Specifying ties.method in matrixStats::rowRanks.

                        Changes in version 1.9.1                        

- Added importAllelicCounts() with options for importing
  Salmon quantification on diploid transcriptomes.

[FLAMES](/packages/FLAMES)
------

                 Changes in version 0.99.0 (2021-04-15)                 

- Submitted to Bioconductor

[FlowSOM](/packages/FlowSOM)
-------

                       Changes in version 2.1.24                        

- Added UpdateMetaclusters function, removed RelabelMetaclusters,
  ReassignMetaclusters and Reordermetaclusters functions.

- Updated CheatSheet

- Added code from UpdateDerivedValues for metaclustersMFIs to
  UpdateFlowSOM

                       Changes in version 2.1.23                        

- Added checkNames = FALSE in MetaclusterMFIs

                       Changes in version 2.1.22                        

- Reordered code in UpdateDerivedValues, RelabelMetaclusters,
  ReorderMetaclusters and ReassignMetaclusters.

                       Changes in version 2.1.21                        

- Added ReorderMetaclusters, to reorder the metacluster levels.

                       Changes in version 2.1.20                        

- Updated PlotManualBars, Plot2DScatters and FlowSOMmary so that it
  works
  with relabeled metaclusters.

                       Changes in version 2.1.19                        

- AggregateFlowFrames accepts channels and markers

- AggregateFlowFrames now gives a warning when files do not contain the
  same
  number of channels

- AggregateFlowFrames now gives warnings when files do not contain the
  same
  markers

- Bugfix in AggregateFlowFrames now works when one channel is given

- Bugfix in PlotFileScatters now works when one channel is given

- Added silent parameter in PlotFileScatters to stop messages

- PlotFileScatters supports channels and markers now

- Add info to FlowSOM object: date when flowSOM object is made, FlowSOM
  verion
  and arguments given to FlowSOM call

- Fixed bug in PlotManualBars

- Added silent parameter in NewData. GetFeatures' silent parameter now
  also
  surpresses message from NewData (more concrete: ReadInput)

                       Changes in version 2.1.17                        

- Added ReassignMetaclusters, to rename or split metaclusters

- Fixed issue where a lot of warnings were printed in FlowSOMmary

- PlotFilescatters now makes filenames unique if they are not and the
  function now works with output of AggregateFlowFrames

                       Changes in version 2.1.16                        

- PlotManualBars allows input of NewData function

                       Changes in version 2.1.15                        

- Fixed warnings with ggtexttable in FlowSOMmary

                       Changes in version 2.1.13                        

- Added RelabelMetaclusters

- PlotFileScatters now has a parameter to change the y-axis label to
  markers
  and/or channels (yLabel)

- Now TRUE/FALSE vector is accepted as input in GetMarkers/GetChannels

                       Changes in version 2.1.11                        

- Added example to AddAnnotation

- Added example to NClusters, NMetaclusters

- Changed examples that used fsom to flowSOM.res

- Added textColor and textSize to AddLabels and PlotNumbers, PlotLabels

- PlotNumbers can plot clusters and metaclusters with parameter "level"

- In GetFeatures, the population parameter is changed to level

- Added GetCounts and GetPercentages to get counts or percentages
  respectively
  per cluster or metacluster

- FlowSOMmary doesn't crash anymore with a column with the same values
  in
  heatmap

- Included a print function for FlowSOM class

- Fixed bug in PlotManualBars

- PlotMarker also accept multiple markers now

                        Changes in version 2.1.8                        

- Solved issue when matrix with no column was given to the SOM function

                        Changes in version 2.1.5                        

- Scale parameter in FlowSOM function defaults to FALSE.

- FlowSOM wrapper function now returns the FlowSOM object instead of a
  list
  containing the FlowSOM object and a metaclustering

- The metaclustering is now found as an element in the flowSOM object.
  Also the
  number of metaclusters and the MFI values are stored and can be
  accessed by
  the NMetaclusters() and GetMetaclusterMFIs() functions.

- If you want to reuse FlowSOM objects generated by previous versions,
  you can use the UpdateFlowSOM function.

- FlowSOM now uses nClus = 10 as default instead of maxMeta = 10

- FlowSOM now makes use of ggplot2 for plotting. PlotFlowSOM provides
  the
  main structure, and has parameters to adapt nodeSize, view (grid, MST
  or some
  own layout matrix), ... PlotStars etc build on this by adding
  additional
  layers to the ggplot object. This also allows to easily incorporate
  multiple
  plots in all layout-tools such as ggarrange, cowplot, patchwork, ...

- GetChannels/GetMarkers can now also take a FlowSOM object as input
  instead of
  a flowFrame.
  New functions:

- To easily generate a clear summary of the model with multiple plots,
  you
  can now use the FlowSOMmary function, which creates a pdf file.

- GetFeatures allows to map new files (internally using the NewData
  function)
  and can return cluster counts, percentages and MFI values for each
  individual
  sample.

- PlotFileScatters can be useful to get an overview of potential batch
  effects
  before running the FlowSOM algorithm

[fobitools](/packages/fobitools)
---------

                        Changes in version 1.1.9                        

- Minor changes in package vignette titles.
- Minor changes in parse_fobi() function.
- Update all data files to FOBI 1.5 version.

fobitools 1.0.0

- Released to Bioconductor 3.13.

[gCrisprTools](/packages/gCrisprTools)
------------

                        Changes in version 2.0.0                        

MAJOR UPDATE Adding the Following Functionality:

- Format, Library, and Testing Improvements 
  
  o Enable processing of libraries with 1:many reagent:target assignments 
  
  o Standardization and clarification of Annotation objects and symbol/identifier
relationships 
  
  o Implementation of factored quantile normalization for timecourse screens 
  
  o Introduction of the simpleResult format and integration with associated functions 

  o Conditional testing framework for quantifying and visualizing signal agreement between contrasts

- Transition to gene set enrichment testing via Sparrow

  o Implement wrappers and provide recommendations fopr geneset enrichment testing in pooled screens 

  o Implementation of GREAT-style pathway mapping for libraries with heterogenous target:gene mappings
  
  o Summarization tools for comparing enrichment signals across
screens

- New Visualization and Interpretation Tools

  o Signal Summary Barchart (Single or Multiple Contrasts) 
  
  o Waterfall reagent/target/pathway visualization (Single Contrast) 
  
  o Contrast comparison plots:

    - Concordance at the Top (CAT)
    
    - Probability Space scatter plots
    
    - UpSet plots with conditional overlap framework

[GENESIS](/packages/GENESIS)
-------

                       Changes in version 2.23.9                        

- Subset covariance matrix to specified samples when sample.id
  argument is passed to fitNullModel when called with an
  AnnotatedDataFrame

                       Changes in version 2.23.8                        

- Added option for recessive and dominant coding to
  assocTestSingle.

                       Changes in version 2.23.7                        

- Implement MatrixGenotypeReader method for pcair by writing a
  temporary GDS file.

                       Changes in version 2.23.5                        

- assocTestSingle, assocTestAggregate, admixMap, and pcrelate use
  the BiocParallel package for parallel execution on blocks of
  variants.

                       Changes in version 2.23.4                        

- For assocTestAggregate, the total number of genotypes in a
  single iterator element (NxM where N=number of samples and
  M=number of variants) may be >2^31.

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 1.6.0                        

New features

- GeneTonic can now accept the input of clusterProfiler's gene set
enrichment analysis functions (gseGO and GSEA), as implemented in
the shake_gsenrichResult() function
- Below each plot and interactive widget, we provide a button that
opens up a modal window where the code required to reproduce that
output is shown as a snippet. These can be readily copied in
extended reports or used to document the exploratory process.

Other notes

- The manuscript about GeneTonic is now available on bioRxiv at
https://www.biorxiv.org/content/10.1101/2021.05.19.444862v1 - the
citation item has been updated accordingly
- GeneTonic's Shiny app now uses the latest version of bs4Dash, which
introduced some breaking changes. Most elements should be now
available as they were in the original implementation

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.30.0                        

NEW FEATURES

- Register NCBI assemblies:
  - mRatBN7.2
  - UMICH_Zoey_3.1
  - Callithrix_jacchus_cj1700_1.1
  - MU-UCD_Fhet_4.1 (GCA_011125445.2)

- Register UCSC genomes:
  - rn7
  - canFam5
  - calJac4

SIGNIFICANT USER-VISIBLE CHANGES

- UCSC hg38 genome is now based on GRCh38.p13 instead of GRCh38.p12

- UCSC mm10 genome is now based on GRCm38.p6 instead of GRCm38

- seqlevelsStyle() setter now issues a warning when some seqlevels
  cannot be switched.



[GenomicAlignments](/packages/GenomicAlignments)
---------------

                      Changes in version 1.30.0                        

- No changes in this version

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.46.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Small update to makeTxDbFromGFF()/makeTxDbFromGRanges():
  makeTxDbFromGFF() and makeTxDbFromGRanges() now recognize and import
  features of type protein_coding_gene (or of type any offspring of
  the protein_coding_gene term) as genes. This was achieved by adding
  protein_coding_gene + its offsprings to
  GenomicFeatures:::.GENE_TYPES.

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.46.0                        

- No changes in this version


[GenomicSuperSignature](/packages/GenomicSuperSignature)
---------------------

                        Changes in version 1.2.0                        

- [Bug Fix]
  - `n` argument of `annotatePC` was hard-coded. Now it can return
  different number of enriched pathways.
  - `abs` argument of `annotatePC` was fixed.
  - Fix wrongfully assigned variable within `plotAnnotatedPCA`
  function.

- [Major]
  - `drawWordcloud` has a new argument `droplist`.
  - Argument name for `plotAnnotatedPCA` is changed from `PCs` to
  `PCnum`.
  - New argument `studyTitle` for `findStudiesInCluster` function.

- [Minor]
  - Description of the package is updated.
  - If non-existing index is provided for any function, it will return
  with
  the error message.

[GeoDiff](/packages/GeoDiff)
-------

                       Changes in version 0.99.5                        

Revisions

- fix platform build/check errors

                       Changes in version 0.99.4                        

Revisions

- address reviewer's comments

                       Changes in version 0.99.3                        

Revisions

- remove seed

                       Changes in version 0.99.2                        

Revisions

- reduce the file size for NBthDEmod2.rda and kidney.rda
- update the unit tests
- change the maintainer

                       Changes in version 0.99.1                        

Revisions

- modify the vignette for Bioc submission

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package
- Package template creation

[GeomxTools](/packages/GeomxTools)
----------

- No changes from 1.1.4

[ggmsa](/packages/ggmsa)
-----

                     Changes in version 2021-09-22                      

- 0.99.1 update DESCRIPTION and NEWS files (2021-09-28, Tue)
- 0.99.2 add documentation for row data in extdata/inst and clean up
code (2021-09-29, Wed)
- 0.99.3 remove some vignettes from master (build on the gh-pages
branch) (2021-10-1, Fri)
- 0.99.4 remove 'stringr' package from 'Imports' (2021-10-11, Mon)
- 0.99.5 make the consensus_views compatible ggtreeExtra and add
package description. (2021-10-21, Thu)

                       Changes in version 0.0.10                        

- update default color schemes in lower part of the SeqDiff plot
(2021-08-20, Fri)

                        Changes in version 0.0.9                        

- import R4RNA to fix R check (2021-08-03, Tue)

                        Changes in version 0.0.8                        

- bugfix: fix variable names error in color_scheme. (2021-07-29, Thu)
- The migration of sequence recombination functionality from seqcombo
package. (2021-07-20, Tue)

                        Changes in version 0.0.7                        

- added gghelix() and geom_helix().(2021-04-1, Thu)
- added option to show the fill legend.(2021-03-23, Tue)
- added a error message to remind that "sequences must have unique
names".(2021-03-18, Thu)
- added ggSeqBundle() to plot Sequence Bundles for MSAs based ggolot2
(2021-03-18, Thu)

                        Changes in version 0.0.6                        

- supports linking ggtreeExtra. (2021-01-21, Thu)
- bugfix: reversed sequence in 'tree + geom_facet(font)' .
(2021-01-21, Thu)
- bugfix: partitioning error when the sequence starting point greater
than 1. (2021-01-21, Thu)
- bugfix: generates continuous x-axis labels for each panel.
(2021-01-21, Thu)
- supports customize colors custom_color. (2020-12-28, Mon)

                        Changes in version 0.0.5                        

- added a new view called by_conservation.(2020-12-22, Tue)
- added a new color scheme Hydrophobicity and a new parameter
border.(2020-12-21, Mon)
- rewrite the function facet_msa().(2020-12-03, Thu)
- Debug: tree + geom_facet(geom_msa()) does not work.(2020-12-03,
Thu)
- added a new function geom_msaBar().(2020-12-03, Thu)
- added a new parameter ignore_gaps used in consensus
views.(2020-10-09, Fri)
- debug in consensus views (2020-10-05, Mon)
- added consensus views (2020-9-30, Wed)
- added new colors LETTER and CN6 provided by ShixiangWang.issues#8

                        Changes in version 0.0.4                        

- fixed warning message in msa_data.R (2020-4-26, Sun)
- added ggplot_add methods for geom_*() (2020-4-24, Fri)
- added a parameter seq_name in ggmsa() (2020-4-23, Thu)
- added a new function facet_msa() --> break down the MSA (2020-4-17,
Fri)
- added a parameter posHighlighted in ggmsa() (2020-4-17, Fri)
- created a new layer geom_asterisk() to optimized geom_seed()
(2020-4-11, Sta)
- added new functions available_colors(), available_fonts() and
available_msa() (2020-3-30, Thu)
- added a new function geom_seed() --> highlight the seed region in
miRNA sequences (2020-3-27, Fri)
- added a new function ggmotif()--> plot sequence motifs
independently
(2020-3-23, Tue)
- added a Monospaced Font DroidSansMono (2020-3-23, Mon)

                        Changes in version 0.0.3                        

- release of v=0.0.3 (2020-03-16, Mon)
- added a new function geom_GC() --> plot GC content in MSA
(2020-02-28, Fri)
- added a new function geom_seqlogo() --> plot plot sequence motifs
in
MSA (2020-02-14, Fri)
- used a proportional scaling algorithm (2020-01-08, Wed)

                        Changes in version 0.0.2                        

- support plot sequence logo (2019-12-25, Wed)
- added three fonts：helvetical, times_new_roman, mono (2019-12-21,
Sta)
- ~~added three fonts：serif_font, Montserrat_font, roboto_font
(2019-12-17, Tue)~~
- added internal outline polygons (2019-12-15, Sun)
- bug fixed of tidy_msa
- import seqmagick for parsing fasta
- tidy_msa for converting msa file/object to tidy data frame
(2019-12-09, Mon)

                        Changes in version 0.0.1                        

- initial CRAN release (2019-10-17, Thu)
- removed from CRAN on 2021-08-17

[ggspavis](/packages/ggspavis)
--------

                 Changes in version 0.99.0 (2021-08-08)                 

- version for submission to Bioconductor

[ggtree](/packages/ggtree)
------

                        Changes in version 3.1.6                        

- geom_cladelab now supports extend parameter as in geom_cladelabel
(2021-10-14, Thu, @xiangpin, #446)
- geom_hilight supports fill linear gradient colour and round rect
background (2021-10-11, Mon; @xiangpin, #449, #444)
- work with negative edge lengths (hclust may generate negative tree
heights) (2021-09-29, Wed; @xiangpin, #441, #445)

                        Changes in version 3.1.5                        

- ggdensitree with align.tips=TRUE sets max x to 0 (2021-09-26, Sun;
@brj1, #437, #439)
- custom column headers for gheatmap (2021-09-15, Wed,
@matt-sd-watson, #434)
- bug fixed of nudge_x parameter in geom_segment2 (2021-09-03, Fri;
@xiangpin, #433)

                        Changes in version 3.1.4                        

- introduce align parameter in geom_hilight (2021-08-30, Mon;
@xiangpin, #431)
- the data parameter in geom_facet now accepts function as input
(2021-08-22, Sun; @xiangpin, #430)
- import ggfun and yulab.utils (2021-08-20, Fri)
- allow using options(layout.radial.linetype) to set linetype of
radial layout (either 'strainght' or 'curved') (2021-08-13, Fri;
@xiangpin, #427)

                        Changes in version 3.1.3                        

- data argument in geom_tiplab and position argument in geom_tree
(2021-08-10, Tue; #426, @xiangpin)
- geom_hilight and geom_cladelab supports function as input data
(2021-07-28, Wed; #421, @xiangpin)
- td_mutate for mutating tree data
- geom_tiplab supports fontface aesthetic (2021-07-06, Tue;
@xiangpin)

                        Changes in version 3.1.2                        

- calculate branch mid point for unrooted layout tree (2021-06-11,
Fri)
- branch.y and branch.x
- geom_range supports aes mapping (2021-06-04, Fri)

                        Changes in version 3.1.1                        

- bug fixed in geom_range (2021-06-01, Tue)
- https://github.com/YuLab-SMU/ggtree/pull/410
- now geom_nodelab has a node="internal" parameter. (2021-05-31, Mon)
- if node = "external", it equivalent to `geom_tiplab
- if node = "all", it equivalent to list(geom_tiplab(),
geom_nodelab())

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                        Changes in version 1.3.6                        

- fix the issue of gridlines when some data is removed. (2021-08-25,
Wed)

                        Changes in version 1.3.5                        

- update citation of ggtreeExtra (2021-08-25, Wed).

                        Changes in version 1.3.4                        

- fix the compute_aes to better compatible with ggplot2 (>=3.3.4)
(2021-08-09, Mon)
- The ggplot2 (>=3.3.4) introduced computed_mapping.

                        Changes in version 1.3.3                        

- update reference. (2021-06-08, Tue)
- fix vector logical check. (201-06-11, Fri)
- c(TRUE, TRUE) && c(TRUE, TRUE) is not allowed in devel
environment of bioconductor

                        Changes in version 1.3.1                        

- data argument of geom_fruit support function input. (2021-05-20,
Thu)
- the argument of axis.params and grid.params can be assigned by
intermediate variables. (2021-05-26, Wed)
- https://github.com/YuLab-SMU/ggtreeExtra/issues/9

                        Changes in version 1.3.0                        

- new version release, and bump new devel version (1.3.0).
(2021-05-20, Thu)

[glmGamPoi](/packages/glmGamPoi)
---------

                         Changes in version 1.5                         

- Choose a more reasonable scale for global overdispersion estimate

- Make code more robust accidental internal NA's

- Add fallback mechanism in case the Fisher scoring fails to converge.
  Instead of returing NA, try again using the BFGS algorithm.

- Better error message if the design contains NA's

[GOfuncR](/packages/GOfuncR)
-------

                       Changes in version 1.13.1                        

USER-LEVEL CHANGES

- update GO-graph (version 01-May-2021)

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.19.1                        

- TCSS method (@qibaiqi, #35; 2021-08-02, Mon)

- GOSemSim 2.18.0

- Bioconductor 3.13 release

[gscreend](/packages/gscreend)
--------

                        Changes in version 1.8.0                        

- Improved parallelization implementation for a faster analysis
  performance.

- resolved knitr error June 2021

[GSEABase](/packages/GSEABase)
--------

                        Changes in version 1.56                         

SIGNIFICANT USER-VISIBLE CHANGES

- goSlim() does not truncate Terms
  (https://github.com/Bioconductor/GSEABase/issues/5)

[GSVA](/packages/GSVA)
----

                        Changes in version 1.50                         

BUG FIXES

- Bugfix for https://github.com/rcastelo/GSVA/issues/54 to force
  filtering genes with constant expression behaving the same regardless
  of the delayed or non-delayed nature of the data container.

[HDF5Array](/packages/HDF5Array)
---------

                Changes in version 1.22.0
                
- No changes in this version

[HDTD](/packages/HDTD)
----

                 Changes in version 1.27.1 (2021-07-28)                 

- Updated CITATION FILE.

- Added rmarkdown to Suggests.

[HGC](/packages/HGC)
---

                 Changes in version 1.1.3 (2021-07-06)                  

- Modify the plotting functions.

                 Changes in version 1.1.2 (2021-06-14)                  

- Add APIs for Seurat pipeline and igraph pipeline.

                 Changes in version 1.1.1 (2021-05-27)                  

- Modify the package manual.

                 Changes in version 1.1.0 (2021-05-20)                  

- The development version 1.1.0 is available in Bioconductor.

[hiAnnotator](/packages/hiAnnotator)
-----------

                       Changes in version 1.27.1                        

- Suggest 'markdown' package in DESCRIPTION and utilize hg19 as default
  freeze where applicable.

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.30.0                        

- add the support of Intel AVX-512VPOPCNTDQ intrinsics (faster than
  AVX512BW)

[HiLDA](/packages/HiLDA)
-----

                 Changes in version 1.7.6 (2020-10-13)                  

- Trigger bioc build

                 Changes in version 1.7.5 (2020-10-08)                  

- Update R version dependency from 3.6 to 4.1

                 Changes in version 1.7.4 (2020-10-08)                  

- Fix a typo

- Add citation file

                 Changes in version 1.7.3 (2020-10-07)                  

- Prepare for release

                 Changes in version 1.7.2 (2020-10-07)                  

- Fix warnings in unit tests

                 Changes in version 1.7.1 (2020-10-07)                  

- Update the man files

[hiReadsProcessor](/packages/hiReadsProcessor)
----------------

                       Changes in version 1.29.1                        

- Suggest 'markdown' package in DESCRIPTION

[hpar](/packages/hpar)
----

                        Changes in version 1.35                         

Changes in version 1.35.1

- Suggest rmarkdown

Changes in version 1.35.0

- New Bioc devel version

[HPAStainR](/packages/HPAStainR)
---------

                        Changes in version 1.2.1                        

- Caught an error where not having `Not detected` column breaks the
  function

- Included 'rmarkdown' package in Suggests

[HPiP](/packages/HPiP)
----

                Changes in version 0.99.17 (2021-09-16)                 

- Updated the vignette examples
- Changed the get_positivePPI() function

                Changes in version 0.99.16 (2021-09-15)                 

- Updated the vignette Figures
- Changed the pred_ensembel() function

[hypeR](/packages/hypeR)
-----

                        Changes in version 1.9.1                        

- msigdb_download() filters by distinct gene symbols within genesets
(relevant to msigdb >=7.4.1)

                        Changes in version 1.9.0                        

- Version bump for bioconductor

[ideal](/packages/ideal)
-----

                       Changes in version 1.18.0                        

Bug fixes

- Fixed the unexpected behaviour when decorating the signature
heatmap
with colData elements

Other notes

- Some alignments of the UI elements in the ideal() app, to harmonize
the content of the page
- Updated the icons in the user interface to match recent changes in
the names from FontAwesome

[IgGeneUsage](/packages/IgGeneUsage)
-----------

                 Changes in version 1.7.24 (2021-10-14)                 

- minor model update: top level hyperpriors on alpha gamma par.
  means not needed anymore. Compare inst/extdata/zibb.stan vs.
  inst/extdata/zibb_original.stan

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 0.99.9 (2021-10-22)                 

- Removed certain unit tests for 64-bit windows

                 Changes in version 0.99.8 (2021-10-11)                 

- Added patch detection method

                 Changes in version 0.99.7 (2021-10-10)                 

- Added all unit tests

- Fixed read_steinbock x/y axis defaults

                 Changes in version 0.99.0 (2021-09-14)                 

- Bioconductor submission

                 Changes in version 0.3.12 (2021-09-01)                 

- clean up for Bioconductor submission

                 Changes in version 0.3.11 (2021-08-19)                 

- Added flip_x, flip_y argument for plotSpatial function

- readSCEfromTXT does not require spot names anymore

- knn graph construction can be pruned by distance

- added

                 Changes in version 0.3.10 (2021-08-18)                 

- Added countInteractions function

- Added testInteractions function

                 Changes in version 0.3.9 (2021-07-29)                  

- Added buildSpatialGraph function

- Added plotSpatial function

- Added aggregateNeighbors function

                 Changes in version 0.3.8 (2021-06-30)                  

- Adjusted default parameters for read_steinbock function

- Added updated test data

                 Changes in version 0.3.7 (2021-06-14)                  

- added read_cpout function, docs and tests

                 Changes in version 0.3.6 (2021-06-04)                  

- added read_steinbock function, docs and tests

                 Changes in version 0.3.5 (2021-05-20)                  

- unit tests and docs for filterPixels

                 Changes in version 0.3.4 (2021-05-18)                  

- added readImagefromTXT function and tests

                 Changes in version 0.3.3 (2021-05-17)                  

- unit tests for binAcrossPixels

                 Changes in version 0.3.2 (2021-05-16)                  

- adjusted plotSpotHeatmap function and unit test

                 Changes in version 0.3.1 (2021-05-15)                  

- readSCEfromTXT accepts list and path

                 Changes in version 0.3.0 (2021-05-07)                  

- Added helper functions for spillover correction

- Removed redundant functions

                 Changes in version 0.2.0 (2019-11-28)                  

- The functions calculateSummary, plotCellCounts, and plotDist have
  been added

                 Changes in version 0.1.0 (2019-09-17)                  

- initial commit

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.25.2                        

- CHANGES
  * bugfix immunoMeta contructor from single immunoClust-object

                       Changes in version 1.25.1                        

- CHANGES
  * bugfix im subset.immunoMeta-object

[immunotation](/packages/immunotation)
------------

                 Changes in version 1.1.1 (2021-08-09)                  

- Resolved error due to MHC reference list reading.

[infercnv](/packages/infercnv)
--------

                 Changes in version 1.8.1 (2020-08-16)                  

- Fix name generation error for step 15.

- Handle annotation names as strings even if they look like numbers.

[InteractiveComplexHeatmap](/packages/InteractiveComplexHeatmap)
-------------------------

                        Changes in version 1.1.4                        

- For the three div blocks of heatmap widgets, now `display:table-cell`
  is used so that
  the positions of divs won't change when changing the size of the
  browser window.

- Add a new vignette "Share interactive heatmaps to collaborators".

                        Changes in version 1.1.3                        

- fontawesome is directly from the fontawesome package.

- also inherit row_names_gp and column_names_gp from the complete
  heatmap

- content of js and css for specific heatmap is directly add to html
  instead of
  containing as files

                        Changes in version 1.1.2                        

- Add `save` argument in `htShiny()`.

                        Changes in version 1.1.1                        

- add new argument `sub_heatmap_cell_fun` and `sub_heatmap_layer_fun`
  to only set cell_fun
  or layer_fun for sub-heatmaps.

[InterCellar](/packages/InterCellar)
-----------

                 Changes in version 2.0.0 (2021-10-20)                  

- Added upload of up to 6 datasets that can be analyzed in parallel
- Added CellChat and ICELLNET as supported tools
- Added multiple conditions comparison in data-driven analysis
- Added significant functional terms in multiple conditions
- Added output folder and file download in automatically generated
folders
- Added weighting for interactions
- Added Network in gene-verse
- Updated About page

[iPath](/packages/iPath)
-----

                        Changes in version 0.99                         

NEW FEATURES

- Initial review.

                       Changes in version 0.99.0                        

- Revise required files and format the code style.

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.28.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Replace dim(), nrow(), and ncol() methods for DataFrameList objects
  with
  dims(), nrows(), and ncols() methods.

DEPRECATED AND DEFUNCT

- Deprecate dim(), nrow(), and ncol() methods for DataFrameList objects
  in favor of the new dims(), nrows(), and ncols() methods.

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.3.9 (2021-10-25)                  

FIXES

- Fixed issues with function that make use of BiocParallel that
sometimes failed on Windows platform

                 Changes in version 1.3.7 (2021-10-20)                  

NEW

- Added new feature iss_source()

FIXES

- Fixed minor issues in data files refGenes_mm9 and function
compute_near_integrations()

                 Changes in version 1.3.6 (2021-10-05)                  

NEW

- Added new feature purity_filter()

FIXES

- Fixed small issue in printing information in reports

                 Changes in version 1.3.5 (2021-09-21)                  

MAJOR CHANGES

- Reworked is_sharing() function, detailed usage in vignette
vignette("sharing_analyses", package = "ISAnalytics")

NEW

- New function cumulative_is()
- New function for plotting sharing as venn/euler diagrams
sharing_venn()

                        Changes in version 1.3.4                        

FIXES/MINOR UPDATES

- Fixed issue in tests that lead to broken build
- Slightly modified included data set for better examples

                 Changes in version 1.3.3 (2021-07-30)                  

MAJOR CHANGES

- Completely reworked interactive HTML report system, for details
take
a look at the new vignette vignette("report_system", package =
"ISAnalytics")
- Old ISAnalytics.widgets option has been replaced by
ISAnalytics.reports
- In remove_collisions(), removed arguments seq_count_col,
max_rows_reports and save_widget_path, added arguments quant_cols
and report_path (see documentation for details)

MINOR CHANGES

- import_single_Vispa2Matrix() now allows keeping additional
non-standard columns
- compute_near_integrations() is now faster on bigger data sets
- Changed default values for arguments columns and key in
compute_abundance()
- compute_near_integrations() now produces only re-calibration map in
*.tsv format
- CIS_grubbs() now supports calculations for each group specified in
argument by
- In sample_statistics() now there is the option to include the
calculation of distinct integration sites for each group (if
mandatory vars are present)

NEW FUNCTIONALITY

- Added new plotting function circos_genomic_density()

FIXES

- Fixed minor issue with NA values in alluvial plots

DEPRECATIONS

- import_parallel_Vispa2Matrices_interactive() and
import_parallel_Vispa2Matrices_auto() are officially deprecated in
favor of import_parallel_Vispa2Matrices()

OTHER

- The package has now a more complete and functional example data set
for executable examples
- Reworked documentation

                 Changes in version 1.3.2 (2021-06-28)                  

FIXES

- Corrected issues in man pages

                 Changes in version 1.3.1 (2021-06-24)                  

NEW FUNCTIONALITY

- is_sharing computes the sharing of IS between groups
- sharing_heatmap allows visualization of sharing data through
heatmaps
- integration_alluvial_plot allows visualization of integration sites
distribution in groups over time.
- top_abund_tableGrob can be used in combination with the previous
function or by itself to obtain a summary of top abundant
integrations as an R graphic (tableGrob) object that can be combined
with plots.

MINOR UPDATES

- Added more default stats functions to default_stats
- Added optional automatic conversion of time points in months and
years when importing association file
- Minor fixes in generate_Vispa2_launch_AF

[iSEE](/packages/iSEE)
----

                        Changes in version 2.5.3                        

- Replace icons with fontawesome 5 versions

                        Changes in version 2.5.2                        

- Bugfix for conversion of categorical columns with too many levels
to
numeric ones
- Bugfix for heatmap crashing if columns were ordered by a selection
that was not shown

                        Changes in version 2.5.1                        

- Bugfix for removed panels showing up among the selectable ones.

[iSEEu](/packages/iSEEu)
-----

                        Changes in version 1.5.2                        

- Fix bug causing AggregatedDotPlot to crash when a column selection
was transferred.
- Fix bug in retrieving a feature set

                        Changes in version 1.5.1                        

- Fix spelling typo in man page.

[IsoformSwitchAnalyzeR](/packages/IsoformSwitchAnalyzeR)
---------------------

                Changes in version 1.15.02 (2021-09-14)                 

- Update type: minor.

- Various error message updates

                Changes in version 1.15.01 (2021-09-01)                 

- Update type: minor.

- Version bump due to Bioconductor release.

- preFilter() now applies the gene expression cutoff to both conditions
  instead of the overall average.

- analyzePFAM() was updated to reflect recent updates to the tidyverse
  read_fwf function. It furhtermore now better distinguishes tap
  seperated and fixed with files.

[isomiRs](/packages/isomiRs)
-------

                       Changes in version 1.21.1                        

FIX

- rmarkdown dependcies

- fix bug that will create difference with older versions (#19)

[kebabs](/packages/kebabs)
------

                       Changes in version 1.27.2                        

- fix in evaluatePrediction() (avoid length>1 logical condition)

                       Changes in version 1.27.1                        

- fix to avoid warnings arising from compiled code

- fix in kebabsDemo() to avoid build warnings; executing this function
  now does
  not have any side effects anymore (previously, data objects were
  loaded into the
  global workspace)

                       Changes in version 1.27.0                        

- new branch for Bioconductor 3.14 devel

[limma](/packages/limma)
-----

                       Changes in version 3.50.0                        

- 
  Improve help pages for fry() and barcodeplot().

- 
  Revise Section 18.1.10 of User's Guide so that coloring of X
  and Y genes is consistent between plots.

- 
  Fix bug in the MDS class method of plotMDS() (introduced in the
  previous release) when new `dim.plot` values are set.

- 
  Fix bug in read.idat() when `annotation` contains NA values.

[LOBSTAHS](/packages/LOBSTAHS)
--------

                       Changes in version 1.19.2                        

- Eliminated some duplication resulting from relict "trimmed_databases"
  directory

                       Changes in version 1.18.1                        

- Updated authors and maintainers

- Formalized support for new lipid classes including bile salts, wax
  esters,
  quinones, etc.

[LoomExperiment](/packages/LoomExperiment)
--------------

                       Changes in version 1.10.1                        

NEW FEATURES

- (v 1.10.1) add '/attrs' to metadata, if present

BUG FIXES

- (v 1.10.1) don't drop all metadata when ReducedDims is not present

[LRBaseDbi](/packages/LRBaseDbi)
---------

                        Changes in version 2.4.0                        

- LRBase.XXX.eg.db-type annotation packages for organisms were
  deprecated.

- The distribution of each SQLite file has been changed to the
  AnnotationHub-style.

- By using LRBaseDbi, we can convert SQLite files acquired by the
  AnnotationHub's query function into LRBase objects, which can then be
  used for analysis using LRBase.XXX.eg.db as before.

- makeLRBasePackage function was deprecated.

- .loadLRBaseDbiPkg was deprecated.

[m6Aboost](/packages/m6Aboost)
--------

                 Changes in version 0.99.0 (2021-05-27)                 

- Submitted to Bioconductor

[Maaslin2](/packages/Maaslin2)
--------

                        Changes in version 1.7.1                        

- Update tutorial data files

- Update knitr dependency for bioconductor tests

[maftools](/packages/maftools)
--------

                        Changes in version 3.14                         

NEW FUNCTIONS

- sampleSwaps Given a list BAM files, the function genotypes known
SNPs and identifies potentially related samples.

BUG FIXES

- Return mutCountMatrix output as a matrix Issue: 769

ENHANCEMENTS

- Added showPct argument to oncoplot. Issue: Issue: 771
- Silently return sample order from oncoplot and coOncoplots Issue:
771

[marr](/packages/marr)
----

                        Changes in version 3.13                         

- added feature-label output on May 10, 2021

                        Changes in version 3.12                         

- added feature-label output on April 27, 2021

- added cutoff value on March 06, 2021

[MatrixGenerics](/packages/MatrixGenerics)
--------------

                        Changes in version 1.5.4                        

- Sync API with matrixStats v0.60.1.

                        Changes in version 1.5.2                        

- Sync API with matrixStats v0.60.0.

                        Changes in version 1.5.1                        

- Fix problem with function environment of fallback mechanism
  (<URL:
  https://github.com/Bioconductor/MatrixGenerics/issues/25> and
  <URL: https://github.com/Bioconductor/MatrixGenerics/pull/26>).
  Make sure that packages can use MatrixGenerics with the ::
  notation to call functions from sparseMatrixStats and
  DelayedMatrixStats.

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.1.2 (2021-09-06)                  

- take sample IDs for shinyQC from colnames(se)

- take feature IDs for shinyQC from rownames(se)

- fix error in report.Rmd (change input for
  create_boxplot to se)

                 Changes in version 1.1.1 (2021-08-27)                  

- fix bug in biocrates and maxQuant function

[matter](/packages/matter)
------

                 Changes in version 1.19.1 (2021-10-20)                 

BUG FIXES

- Fix error in 'cbind()' on two 'matter_vec' objects

[megadepth](/packages/megadepth)
---------

                        Changes in version 1.3.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- Make install_megadepth more transparent for users. Interactive
sessions will now prompt users to agree to the install Megadepth at
the printed locations.

[memes](/packages/memes)
-----

                        Changes in version 1.1.3                        

- fixed a bug in runTomTom where setting norc = TRUE failed on data
import

                        Changes in version 1.1.1                        

- runFimo now returns NULL and prints a message when text = FALSE and
FIMO detects no matches instead of throwing a cryptic error message

[MeSHDbi](/packages/MeSHDbi)
-------

                       Changes in version 1.30.0                        

- MeSH-related packages (MeSH.XXX.eg.db, MeSH.db, MeSH.AOR.db, and
  MeSH.PCR.db) were deprecated.

- The distribution of each SQLite file has been changed to the
  AnnotationHub-style.

- By using MeSHDbi, we can convert SQLite files acquired by the
  AnnotationHub's query function into MeSH objects, which can then be
  used for analysis using MeSH-related packages as before.

- makeGeneMeSHPackage was deprecated.

- .loadMeSHDbiPkg was deprecated.

[meshes](/packages/meshes)
------

                       Changes in version 1.19.3                        

- cache mesh db and table (2021-09-01, Wed)

                       Changes in version 1.19.2                        

- import yulab.utils (2021-8-19, Thu)

                       Changes in version 1.19.1                        

- remove MeSH.db package and use AnnotationHub to get MeSHDb
databases
(2021-8-13, Fri)

[meshr](/packages/meshr)
-----

                        Changes in version 2.0.0                        

- Specification changed as "AnnotationHub-style"

- Dependencies of MeSH.db, MeSH.AOR.db, MeSH.PCR.db, MeSH.Hsa.eg.db,
  MeSH.Aca.eg.db, MeSH.Bsu.168.eg.db, MeSH.Syn.eg.db were removed

- Vignette changed for the specification change

- All datasets removed

[metabCombiner](/packages/metabCombiner)
-------------

                        Changes in version 1.3.2                        

Changes to labelRows

- new resolveConflicts + rtOrder argument added for automated
  resolution of
  conflicting feature pair alignment fows in combinedTable, with
  embedded
  C function (resolveRows.c)

- new argument rtOrder paired with resolveConflicts determines if RT
  order
  consistency is expected when resolving alignment conflicts

- duplicate column names for {labels, subgroup, alt} removed and no
  longer
  allowed in resulting combinedTable

Changes to metabCombiner

- new check for metabCombiner object inputs: labelRows must be called
  before
  aligning a metabCombiner object with a new dataset

- resolveConflicts method applied to metabCombiner object to eliminate
  conflicting alignments/ attain 1-1 matches for all features

- new rtOrder argument controlling resolveConflicts similar to above

                        Changes in version 1.3.1                        

Changes to fit_gam()/ fit_loess

- messages changed from using cat() to using message()

Changes to metabCombiner() & combinedTable / featdata slots

- rowID column now the first column in both tables instead of having
  this
  information be the row names

Changes to calcScores() / evaluateParams() / labelRows

- updated to reflect above changes to keep rowID identical between
  combinedTable & featdata

Changes to write2file

- combinedTable & featdata merged by rowID column when metabCombiner
  objects
  used as input

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                         Changes in version 1.1                         

MetaboCoreUtils 1.1.1

- Add correctRindex function.
- Add isotopologue function to group isotopologues in MS spectra.

[metaseqR2](/packages/metaseqR2)
---------

                 Changes in version 1.5.1 (2021-06-14)                  

NEW FEATURES

- None.

BUG FIXES

- Fix in annotation backwards compatibility.

- Fix bigGenePred gene annotation track generation.

[methylclock](/packages/methylclock)
-----------

                       Changes in version 0.99.27                       

- New package in Bioconductor 3.14 release

[MetNet](/packages/MetNet)
------

                 Changes in version 1.11.5 (2021-10-12)                 

- add assays in structural based on columns in transformations that are
  defined by the var argument in structural, adjacency matrices of type
  will be stored in the AdjacencyMatrix object defined in the columns
  of
  transformation

- implement structural that it can also calculate mass differences of 0
  for
  undirected networks

                 Changes in version 1.11.4 (2021-09-09)                 

- shift calculation of as.data.frame(am) in mz_summary

- several fixes of typos in the comments and vignette

- fix rtCorrection function

                 Changes in version 1.11.3 (2021-08-30)                 

- change calculation of mass differences, use the differences between
  (M_2+ppm)-(M_1-ppm) and (M_2-ppm)-(M_1+ppm) instead of
  (M_1-ppm)-(M_1) and
  (M_2+ppm)-(M_1) for querying against the list of transformations

                 Changes in version 1.11.2 (2021-08-30)                 

- change error message in test_combine

[mia](/packages/mia)
---

                  Changes in version 1.1 (2021-06-04)                   

- split transformCounts into transformSamples and transformFeatures

- added log_modulo_skewness as a diversity index

- added functions for summarizing dominant taxa information

- added wrapper for adding dominant taxa information to colData

- added specialized subsetting function for subsetting by prevalence
  (subsetByPrevalentTaxa/subsetByRareTaxa)

- added mapTaxonomy

- added estimateDivergence

- bugfix: makePhyloseqFromTreeSummarizedExperiment checks now for
  rowTree be compatible

- bugfix: meltAssay supports Matrix types

- bugfix: meltAssay is able to include rowData also when there are
  duplicated rownames

- added subsampleCounts for Subsampling/Rarefying data

[miaSim](/packages/miaSim)
------

                 Changes in version 0.99.0 (2021-09-29)                 

- Three simulation models and three functions are added

- Submitted to Bioconductor

[miaViz](/packages/miaViz)
------

                         Changes in version 1.1                         

- Added plotAbundanceDensity (2021-06-23)

[microbiome](/packages/microbiome)
----------

                 Changes in version 2.1.2 (2020-07-01)                  

- Core heatmap labeling improved

- aggregate_top_taxa deprecated

- bimodality and potential_analysis functions fixed

                 Changes in version 2.1.1 (2020-04-06)                  

- Added overlap function

                 Changes in version 1.14.1 (2021-09-29)                 

- Removed categorical method from associate function

[microbiomeMarker](/packages/microbiomeMarker)
----------------

                 Changes in version 0.99.0 (2021-09-01)                 

- Submitted to Bioconductor

[MicrobiomeProfiler](/packages/MicrobiomeProfiler)
------------------

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                        Changes in version 1.5.9                        

- add include.rownames to control whether consider the OTU as
taxonomy
feature table in diff_analysis and get_alltaxadf or tip labels in
as.treedata. (2021-10-19, Tue)
- fix rename bug, rename the taxonomy names can work now.
(2021-10-12,
Tue)
- introduce trimSample in mp_rrarefy to check whether to remove the
samples that do not have enough abundance. (2021-10-11, Mon)
- update MPSE to allow assays supporting data.frame or DFrame class.
(2021-10-08, Fri)
- update mp_plot_ord to suppress the message of the third depend
package. (2021-10-08, Fri)

                        Changes in version 1.5.8                        

- fix the bug of AsIs list class in unnest for the tidyr (>= 1.1.4).
(2021-10-01, Fri)
- add mp_aggregate function. (2021-09-26, Sun)

                        Changes in version 1.5.7                        

- fix bug of mp_plot_upset. (2021-09-10, Fri)
- update the mp_plot_ord. (2021-09-13, Mon)

                        Changes in version 1.5.6                        

- convert the type of first element of assays to matrix to compatible
with DESeqDataSet of DESeq2, test_differential_abundance of
tidybulk. (2021-09-09, Thu)
- update show and print for format output of MPSE class. (2021-09-08,
Wed)
- update mp_cal_abundance use new tidytree. (2021-09-07, Tue)
- introduce include.lowest parameter in mp_filter_taxa. (2021-09-07,
Tue)

                        Changes in version 1.5.5                        

- update mp_plot_ord to display the bioplot for result of cca, rda
and
envfit. (2021-09-06, Mon)
- update the vignettes of MicrobiotaProcess. (2021-09-04, Sat)
- return updated MPSE object after the mp_diff_analysis is done with
action="add". (2021-08-31, Fri)
- then the taxtree and otutree with the result of different
analysis can be extracted with mp_extract_tree.
- fix issue print for one line of MPSE and update mp_plot_ord to
display the side boxplot. (2021-08-31, Tue)
- add mp_plot_ord for MPSE or tbl_mpse object after one of
mp_cal_pca,
mp_cal_pcoa, mp_cal_rda, mp_cal_nmds, mp_cal_rda, mp_cal_cca,
mp_cal_dca or mp_envfit has been run with action='add'. (2021-08-30,
Mon)
- add mp_plot_dist for MPSE or tbl_mpse object after mp_cal_dist is
performed with action="add". (2021-08-28, Sat)
- add mp_plot_abundance, mp_plot_alpha, mp_plot_rarecurve,
mp_plot_venn, mp_plot_upset for MPSE after the corresponding
mp_cal_abundance, mp_cal_alpha, mp_cal_rarecurve, mp_cal_venn,
mp_cal_upset are performed with action="add". (2021-08-27, Fri)
- fix the issue when the rowname or colnames of SummarizedExperiment
is NULL for as.MPSE. (2021-08-26, Thu)

                        Changes in version 1.5.4                        

- fix the rownames of assays and colnames of colData to identical for
SummarizedExperiment(1.23.3). (2021-08-26, Thu)
- add mp_extract_refseq for MPSE object. (2021-08-25, Wed)
- update as.MPSE for SummarizedExperiment object. (2021-08-24, Tue)
- add mp_filter_taxa to drop the taxa that low abundance and low
occurrences. (2021-08-24, Tue)
- add colData<- and left_join for MPSE. (2021-08-23, Mon)
- fix mutate for MPSE object.
- don't import the parse_taxonomy_greengenes and parse_taxonomy_qiime
from phyloseq. (2021-08-17, Tue)
- add as.MPSE for TreeSummarizedExperiment class. (2021-08-17, Tue)
- add mp_import_metaphlan to parsing the output of MetaPhlAn.
(2021-08-12, Thu)
- add treefile argument to import the tree of MetaPhlAn3
(mpa_v30_CHOCOPhlAn_201901_species_tree.nwk) (2021-08-13, Fri)
- update the print of MPSE object via pillar package. (2021-08-06,
Fri)
- update mp_extract_dist by introducing .group argument to return a
tbl_df for visualization. (2021-08-04, Wed)
- add taxatree, taxatree<-, otutree, otutree<-, refseq, refseq<- for
MPSE. (2021-08-04, Wed)
- add mp_extract_rarecurve to extract the result of mp_cal_rarecurve
from MPSE or tbl_mpse object. (2021-08-04, Wed)
- add mp_stat_taxa to count the number and total number taxa for each
sample at different taxonomy levels (Kingdom, Phylum, Class, Order,
Family, Genus, Species, OTU). (2021-08-03, Tue)

                        Changes in version 1.5.3                        

- rename mp_extract_abundance to mp_extract_assays from MPSE or
tbl_mpse. (2021-07-31, Sat)
- update the method to save the result of mp_cal_clust by introducing
action argument. (2021-07-29, Thu).
- update as.phyloseq for MPSE or tbl_mpse object. (2021-07-28, Wed)
- add mp_diff_analysis for MPSE or tbl_mpse object. (2021-07-27, Tue)
- add dr_extract for the visualization of the result of ordination.
(2021-07-26, Mon)
- comment out the function for phyloseq and add rd of the function
for
MPSE or tbl_mpse. (2021-07-24, Sat)
- update the function to parsing the result of rda, cca, envfit.
(2021-07-23, Fri)
- add tidydr to convert the result of reduce dimension to tbl_df
- such pca, pcoa, nmds, rda, cca. (2021-07-22, Thu)
- optimize the print for MPSE. (2021-07-22, Thu)

                        Changes in version 1.5.2                        

- add mp_mantel and mp_mrpp for MPSE or tbl_mpse object. (2021-07-19,
Mon)

- add mp_envfit and update mp_cal_dist to support the distance
calculation with continuous environment factors and rename
mp_cal_adonis to mp_adonis, mp_cal_anosim to mp_anosim. (2021-07-17,
Sat)

- add mp_cal_rda, mp_cal_cca, mp_cal_adonis and mp_cal_anosim for
MPSE
or tbl_mpse object. (2021-07-16, Fri)

- add mp_cal_dca, mp_cal_nmds and mp_extract_internal_attr.
(2021-07-15, Thu)

- add mp_cal_pca, mp_cal_pcoa and mp_extract_abundance. (2021-07-14,
Wed)

- add mp_cal_clust to perform the hierarchical cluster analysis of
samples and mp_extract_dist to extract the dist object from MPSE
object or tbl_mpse object. (2021-07-13, Thu)

- add mp_cal_dist to calculate the distance between samples with MPSE
or tbl_mpse object. (2021-07-12, Mon)

- add mp_extract_sample, mp_extract_taxonomy, mp_extract_feature to
extract the sample, taxonomy and feature (OTU) information and
return tbl_df or data.frame. (2021-07-09, Fri)

- add mp_cal_venn to build the input for venn plot (2021-07-09, Fri)

- mp_cal_rarecurve add action argument to control whether the result
will be added to MPSE and tbl_mpse or return directly. (2021-07-08,
Thu)

- add mp_cal_upset to get the input of ggupset. (2021-07-08, Thu)

- add mp_extract_tree to extract the otutree or taxatree from MPSE or
tbl_mpse object. (2021-07-07, Wed)

- add pull and slice to support the MPSE object. (2021-07-06, Tue)

- add mp_cal_rarecurve to calculate the rarecurve of each sample with
MPSE or tbl_mpse. (2021-07-06, Tue)

- add mp_cal_abundance to calculate the relative abundance of each
taxonomy class with MPSE or tbl_mpse. (2021-07-05, Mon)

- add mp_decostand provided several standardization methods for MPSE,
tbl_mpse and grouped_df_mpse. (2021-07-04, Sun)

- add mp_import_qiime to parse the output of qiime old version.
(2021-07-03, Sat)

- add taxatree slot to MPSE. (2021-06-30, Wed)

- add mp_cal_alpha function for MPSE or mpse object. (2021-07-01,
Thu)

- add rownames<- to support renaming the names of feature.
(2021-07-01, Thu)

- add mp_import_qiime2 and mp_import_dada2 to parse the output of
dada2 or qiime2 and return MPSE object. (2021-07-02, Fri)

- update print information for MPSE, tbl_mpse and grouped_df_mpse.
(2021-06-29, Tue)

- add [ to the accessors of MPSE. (2021-06-29, Tue)

- use MPSE object. (2021-06-28, Mon)

- add as.MPSE to convert phyloseq or tbl_mpse to MPSE class.
- Formatted output.

- tidy framework for MPSE object.
- as_tibble to convert MPSE and phyloseq to tbl_mpse. (2021-06-28,
Mon)
- filter to subset a data frame from tbl_mpse. (2021-06-28, Mon)
- group_by to do some data operations on groups for tbl_mpse.
(2021-06-28, Mon)
- arrange to order the rows of a data frame for tbl_mpse.
(2021-06-28, Mon)
- mutate to adds new variables and preserves existing ones for
tbl_mpse. (2021-06-28, Mon)
- select to select variables in tbl_mpse. (2021-06-28, Mon)
- distinct to select only unique/distinct rows in tbl_mpse.
(2021-06-28, Mon)
- rename to rename the variable names in tbl_mpse. (2021-06-28,
Mon)
- nest to create a list-column of tbl_mpse, it will convert
tbl_mpse to tbl_mpse_nest. (2021-06-28, Mon)
- unnest to convert the tbl_mpse_nest to tbl_mpse. (2021-06-28,
Mon)
- as.treedata to convert tbl_mpse to treedata, then we can explore
the data with treedata.
- add tiplevel argument to control whether use OTU as tip
label, default is OTU. (2021-06-28, Mon)
- left_join to mutate joins based the left tbl_mpse structure.
(2021-06-28, Mon)
- changed clustplotClass to treedata. (2021-06-28, Tue)
- add mp_rrarefy method to rarefy species richness. (2021-06-29, Tue)
- it supports MPSE, tbl_mpse, grouped_df_mpse object via wrapping
vegan::rrarefy.
- ~~update as.MPSE and as.treedata for grouped_df_mpse object.
(2021-06-29, Tue)~~ ~~- This feature is useful to explore the
microbiome data in taxa tree.~~ This feature has been replaced by
the taxatree slot

                        Changes in version 1.5.1                        

- add ellipse_linewd and ellipse_lty in ggordpoint to control the
width and line type of ellipse line. (2021-05-24, Mon)
- fixed the regular expression match for the internal function to
print the results of diff_analysis. (2021-06-06, Sun)
- add filter function to filter the result of diff_analysis.
(2021-06-07, Mon)
- more accessor function for result of diff_analysis. (2021-06-07,
Mon)
- head
- tail
- [
- [[]]
- $
- dim
- add get_NRI_NTI to calculate the NRI and NTI. (2021-06-08, Tue)

[midasHLA](/packages/midasHLA)
--------

                        Changes in version 1.1.0                        

- Bioconductor release!

                       Changes in version 1.0.13                        

- adds new HLA-KIR interaction A03_A11_KIR3DL2 defined as (A03 | A11)
& KIR3DL2.

                       Changes in version 1.0.12                        

- fixes bug causing MiDAS subsetting to break omnibus testing.

                       Changes in version 1.0.11                        

- runMiDAS inheritance_model argument is no longer by defaut
'additive'. Now it is required to specify desired inheritance model,
when appplicable.

                       Changes in version 1.0.10                        

- fix bug causing Bw6 groups to be counted twice in hla_NK_lingads
experiment.

                        Changes in version 1.0.9                        

- fix bug causing runMiDAS errors when statistical model evaluated
with a warrning.

                        Changes in version 1.0.8                        

- fixed bug causing filterByVariables and filterByFrequency to strip
omnibus groups from target experiment.

                        Changes in version 1.0.7                        

- fixed bug causing HWETest filtration to strip omnibus groups from
target experiment

                        Changes in version 1.0.6                        

- removed unused expression dictionaries

                        Changes in version 1.0.4                        

- In frequency calculations the "NA"s were counted as non-carriers,
this has been changed such that "NA" samples are now omitted.

                        Changes in version 1.0.3                        

- warnings and errors occuring upon model evaluation are now
summarized into more readable form

                        Changes in version 1.0.2                        

- fixed problem vignettes index entry values, preventing vignettes
from being build

                        Changes in version 1.0.1                        

- fixed bug in summariseAAPos, where argument specifying AA position
didn't consider AA position numbering starting from negative
positions
- frequencies in getFrequencies output are no longer formatted as
percentages
- kableResults scroll box height can now be adjusted
- omnibus result columns: dof, residue were renamed to df, residues
- missing Bw6 references were added to allele_HLA_Bw dictionary
- new inheritance model has been added the overdominance

[miloR](/packages/miloR)
-----

                 Changes in version 1.1.0 (2021-10-12)                  

- Fix bug in testNhoods to use user-specific reduced dimensions
- Vignettes now include set rownames() to avoid confusion
- Numerous doc-string typo fixes

[minfi](/packages/minfi)
-----

                        Changes in version 1.39                         

- v1.39.1 Initial support for the Allergy and Asthma array.

- v1.39.2 More support for the Allergy and Asthma array.

- v1.39.3 Bug fix that prevented R CMD build from working.

[miQC](/packages/miQC)
----

                 Changes in version 1.1.5 (2021-08-24)                  

- Updated citation

                 Changes in version 1.1.3 (2021-05-27)                  

- Added new option for model_type, one_dimensional

- Added new filtering parameters, keep_all_below_boundary and
  enforce_left_cutoff

- Added demonstrations of new models and parameters to vignette

[mirTarRnaSeq](/packages/mirTarRnaSeq)
------------

                 Changes in version 1.1.3 (2021-08-22)                  

- Changes made in the vignette

                 Changes in version 1.1.2 (2021-08-22)                  

- Link example data on zenodo in vignette

                 Changes in version 1.1.1 (2021-06-03)                  

- Addition to miRanda Files

[mistyR](/packages/mistyR)
------

                        Changes in version 1.2.0                        

- Release version for Bioconductor 3.14. See changes for 1.1.x.

                         Changes in version 1.1                         

- Added funtions for view manipulation, including view filtering and
marker selection.
- Added functions for performance, contribution and importance
signature extraction from results.
- Aggregation and signature generation is generalized for samples
with
non-identical targets by working on the intersection.
- Modeling of intraview can be bypassed.
- Added families of distances to calculate paraview.
- Paraview can exlude measurements within a used defined zone of
indifference around each spatial unit.
- Improved plotting control.
- Complete test coverage.

IMPORTANT: R2 is now reported in percentages for intra, multi and
gain.
Collecting results from running mistyR < 1.1.11 will lead to
miscalcuation of gain.R2. Update the performance.txt files by
multiplying the values in columns intra.R2 and multi.R2 by 100.

                        Changes in version 1.0.3                        

- Fixed display of messages and progress during view generation.
- Improved plotting control and display.
- Fixed handling of NaN in results.
- Vignette output switched from BiocStyle to rmarkdown for pdfs due
to
BiocStyle issue.

                        Changes in version 1.0.2                        

- Bugfix: models built with different parameters stored and retrieved
from the same cache file.
- Avoid calls to os-dependent file.info in tests.

                        Changes in version 1.0.1                        

- Bugfix: passing arguments to ranger.
- Warnings on clearing nonexistent cache folders and tests of
performance.
- Increased test coverage.

[mixOmics](/packages/mixOmics)
--------

                       Changes in version 6.18.0                        

new features / enhancements / changes

- new function plotMarkers to visualise the selected features in
block
analyses (see #134)
- auroc title now fixed (#135)
- cimDiablo takes trim argument to customise outlier filtering (#136)
- plotIndiv.pca default shape set to 16
- circosPlot & network now support blocks with similar feature names
- circosPlot now has methods for block.spls objects
- circosPlot now takes new formal and advanced args for
customisation.
See ?circosPlot.

bug fixes

- plotVar legend colour mismatch bug fixed
- plotDiablo error undefined variable (Y) fixed
- nipals initialisation bug with high-variance high-NA rate column
fixed
- cim bug with high NA rate data fixed using imputation by the column
mean

[monaLisa](/packages/monaLisa)
--------

                       Changes in version 0.99.5                        

- Updated R/monaLisa-package.R file

                       Changes in version 0.99.4                        

- Suppressed warnings from matchPWM (due to presence of Ns) in
regression vignette

                       Changes in version 0.99.3                        

- Updated README.md file

                       Changes in version 0.99.2                        

- Added fixes to the regression vignette
- Addressed failing test in calcBinnedKmerEnr

                       Changes in version 0.99.1                        

- Added examples where missing for exported functions
- Harmonized function naming (anno_seqlogo -> annoSeqlogo,
sample_random_regions -> sampleRandomRegions)
- Clarified details on Pearson residual calculation
- Adapted documentation for new version of BiocParallel
- Harmonized return values from plot functions
- Added legend position and size arguments to plotSelectionProb()

                       Changes in version 0.99.0                        

- Preparation for Bioconductor submission

                        Changes in version 0.2.0                        

- Added a NEWS.md file to track changes to the package.

[mosbi](/packages/mosbi)
-----

                        Changes in version 1.0.0                        

Major changes

- First public version of mosbi. Future changes will be reported
here.

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.37.5                        

- Change the style of motifPile.

- Fix the bug of the ylab grid.

                       Changes in version 1.37.4                        

- Fix the bug that the rownames were not checked for alignment.

                       Changes in version 1.37.3                        

- Fix the bug method argument of matAlign is ignored.

                       Changes in version 1.37.2                        

- Improve importMatrix to read the tags from file.

                       Changes in version 1.37.1                        

- handle the error "failed to load cairo DLL"

[MQmetrics](/packages/MQmetrics)
---------

                        Changes in version 1.1.8                        

- Fixed test that was resulting in error duet to version 1.1.6
updated
way to read the files.

                        Changes in version 1.1.7                        

- Removed parentheses from the news that was causing issues in
Bioconductor.

                        Changes in version 1.1.6                        

- Bug fixed: If a table was missing, the report was not generated.

                        Changes in version 1.1.5                        

- The plot generated by PlotPTM() will now indicate (in the legend
title) whether the Post-Translational modifications have been
aggregated or not as a result of the parameter aggregate_PTMs.

- The function PlotPeptidesIdentified() and PlotProteinsIdentified()
now return a plot containing Missing Values, Frequency of Identified
by Match Between Runs and Frequency of identified by MS/MS. With
this, the funciton PlotIdentificationType() becomes obsolete.

- The function PlotProteinCoverage() now can take as input multiple
UniprotID in a vector format.

- The function PlotPTMAcrossSamples() now can take as input multiple
PTM_of_Interest in a vector format.

- Change in the function make_MQCombined() to read faster the tables
and reducing the overall time required to generate a report.

                        Changes in version 1.1.4                        

- The function PlotProteinCoverage() now reports the coverage
individually in each plot rather than the total protein coverage.

- MQmetrics now is adapted to MaxQuant v.2.x, since the column names
are different than in MaxQuant v.1.x. MQmetrics will detect the
MaxQuant version used and read the columns accordingly.

- Enhanced error message in PlotiRT() and PlotiRTScore() when irt
peptides are note found. Enhanced error message for
PlotProteinCoverage() when the UniprotID is not found.

                        Changes in version 1.1.3                        

- Added new function PlotPTMAcrossSamples(), it takes as input one
PTM
of interest and shows its intensities across the samples. This
function is similar to PlotPTM() but in more detail.
- In the function PlotPTM() a parameter combine_same_residue_ptms has
been added. It combines multiple PTMs happening in the same residue
such as: Dimethyl (KR), Trimethyl (KR).

                        Changes in version 1.1.2                        

- Fixed units of time MaxQuantAnalysisInfo() when experiment lasting
longer than a day.
- Added new line to MaxQuantAnalyssInfo() showing when the experiment
ended.
- Improved aesthethics in the plots from PlotCombinedDynamicRange()
and PlotAllDynamicRange().

                        Changes in version 1.1.1                        

- Added pagination to PlotiRT() and PlotiRTScore().
- Updated vignette style to Bioconductor's.
- Improved aesthethics of PlotProteinOverlap() and PlotPCA().

[msa](/packages/msa)
---

                       Changes in version 1.25.3                        

- further changes to get rid of compiler warnings

                       Changes in version 1.25.2                        

- removed build/ directory from repo to avoid installation problems

                       Changes in version 1.25.1                        

- update of gc

- minor changes to get rid of compiler warnings

                       Changes in version 1.25.0                        

- new branch for Bioconductor 3.14 devel

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                         Changes in version 1.1                         

Changes in 1.1.4

- Fix wrong database column name for collision energy.

Changes in 1.1.3

- Change SQL queries to increase performance.

Changes in 1.1.2

- Fix bug in show,MsBackendMassbankSql.

Changes in 1.1.1

- Cache precursor m/z to allow faster queries for spectral matching.

[MsBackendMgf](/packages/MsBackendMgf)
------------

                         Changes in version 1.1                         

Changes in 1.1.3

- Fix issue with MGF files lacking peak data.

Changes in 1.1.2

- Export precursor charge in the expected format (issue #16).

Changes in 1.1.1

- Add an example to the vignette describing how to export only
selected spectra variables to the MGF file.

[MsBackendRawFileReader](/packages/MsBackendRawFileReader)
----------------------

                        Changes in version 1.0.0                        

- First bioconductor release of the MsBackendRawFileReader
  package.

[msImpute](/packages/msImpute)
--------

                        Changes in version 1.3.0                        

- Users can now specify the rank of the model to fit by msImpute

- Added mspip for identification transfer between runs using
  Maxquant results (Beta phase only)

- Added evidenceToMatrix which creates limma compatible objects
  from MaxQuant evidence table

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.19                         

Changes in 2.19.2

- Fix plot with type = "XIC" to create an empty plot for samples
without data points (issue #549).

Changes in 2.19.1

- Add compareChromatograms method.

[msPurity](/packages/msPurity)
--------

                       Changes in version 1.19.2                        

- XCMS 3 compatability update (M-R-JONES)
  https://github.com/computational-metabolomics/msPurity/pull/91

                       Changes in version 1.19.1                        

- Bug fix for flagRemove (full width was not calculated as expected)

[msqrob2](/packages/msqrob2)
-------

                        Changes in version 1.1.1                        

- Fix filtering steps in vignette: now using
QFeatures::filterFeatures()

[MSstatsConvert](/packages/MSstatsConvert)
--------------

                        Changes in version 1.2.1                        

- Fixed a bug related to SRM inputs.

[MSstatsTMT](/packages/MSstatsTMT)
----------

                 Changes in version 2.2.3 (2021-10-06)                  

- Minor change: fix the bug when df.prior is infinite

                 Changes in version 2.2.2 (2021-09-21)                  

- Major change: extend groupComparisonTMT() function to cover
  repeated measures design

- Allow flexible order of condition in dataProcessPlotsTMT.

- Fix bug in Condition label in dataProcessPlotsTMT.

- Improve MSstatsTestSingleProteinTMT() by directly reading
  lmerTest output. This may make statistics slightly different
  due to different numeric accuracy

- fix bug when condition name contains 'group'

- change the x-axis order in profile plot

                 Changes in version 2.0.1 (2021-06-14)                  

- update comments of PD converter function

- fix bug in proteinSummarization() function when MBimpute = F

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.20.0                        

Bug fixes and minor improvements

- Avoid dropping experiments with repeated calls to subsetByColData
and remove harmonization (@cvanderaa, #302)
- getWithColData suppresses messages from natural subsetting
operations by default with verbose = FALSE (@bhagwataditya, #301)
- getWithColData was using the old default (drop = TRUE) and causing
an error when the experiment is empty (@danielinteractive, #300).
- Calls to the internal .harmonize operation are reduced to increase
memory efficiency, when identical experiment colnames present
(@LTLA, #299).
- subsetByColData now errors on subscript vectors longer than the
nrow
of the colData (previously a warning).
- colData<- includes a check for identical row names. If so, direct
replacement of the colData occurs without harmonization.
- Added a warning when an empty sampleMap is provided in the
constructor function which may cause unexpected behavior.
Documentation is updated to include more details on the sampleMap
input.

[multiHiCcompare](/packages/multiHiCcompare)
---------------

                 Changes in version 2.0.0 (2020-11-23)                  

- Update selection of significant results in the 'topDirs' function.
  Major change,
  results will be more strict compared to pre-2.0.0 version

- Removed 'BLMA' and 'metap' dependency, added 'aggregation' dependency

- P-values are aggregated using 'max' by default. I.e., for a region
  differentially
  interacting with multiple regions, a maximum p-value will be
  selected. Fisher,
  Lancaster, and Sidak methods are also available

- Harmonize counting of significant regions using 'p.adj_cutoff' only
  ('alpha'
  cutoff removed)

- The 'manhattan' function is harmonized with 'topDirs'

- The 'perm_test' function is harmonized with 'topDirs'

- Update vignettes to match functions

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.1.27                        

Bug fixes

- validate_parameters can now handle ref_genome=NULL
- .tsv.gz no longer assigned suffix .tsv.
- Made code width <80 characters.
- Changed to_GRanges/to_GRanges functions to all-lowercase functions
(for consistency with other functions).
- Set nThread=1 in data.table test functions.

New Features

- Added tests for get_genome_builds
- Added early check for making sure the directory save_path is in was
actually created (as opposed to finding out at the very end of the
pipeline).
- Tabix-indexing now available for tabular output data.
- read_header and read_sumstats now both work with .bgz files.

                       Changes in version 1.1.26                        

New Features

- Extra mappings for FRQ column, see data("sumstatsColHeaders") for
details

                       Changes in version 1.1.23                        

New Features

- format_sumstats(FRQ_filter) added so SNPs can now be filtered by
allele frequency
- Mapping file now has mappings for allele frequency (AF) to FRQ
- VCF files with AF in INFO column e.g. 'AF=...' now converted to AF
column
- format_sumstats(frq_is_maf) check added to infer if FRQ column
values are minor/effect allele frequencies or not. frq_is_maf allows
users to rename the FRQ column as MAJOR_ALLELE_FRQ if some values
appear to be major allele frequencies

                       Changes in version 1.1.19                        

New Features

- get_genome_builds() can now be called to quickly get the genome
build without running the whole reformatting.
- format_sumstats(compute_n) now has more methods to compute the
effective sample size with "ldsc", "sum", "giant" or "metal".
- format_sumstats(convert_ref_genome) now implemented which can
perform liftover to GRCh38 from GRCh37 and vice-versa enabling
better cohesion between different study's summary statistics.

                       Changes in version 1.1.11                        

Bug fixes

- check_no_rs_snp can now handle extra information after an RS ID. So
if you have rs1234:A:G that will be separated into two columns.
- check_two_step_col and check_four_step_col, the two checks for when
multiple columns are in one, have been updated so if not all SNPs
have multiple columns or some have more than the expected number,
this can now be handled.
- Extra mappings for the FRQ column have been added to the mapping
file

New Features

- check_multi_rs_snp can now handle all punctuation with/without
spaces. So if a row contains rs1234,rs5678 or rs1234, rs5678 or any
other punctuation character other than , these can be handled.
- format_sumstats(path) can now be passed a dataframe/datatable of
the
summary statistics directly as well as a path to their saved
location.
- Input summary statistics with A0/A1 corresponding to ref/alt can
now
be handled by the mappign file as well as A1/A2 corresponding to
ref/alt.

                        Changes in version 1.1.2                        

New Features

- import_sumstats reads GWAS sum stats directly from Open GWAS. Now
parallelised and reports how long each dataset took to import/format
in total.
- find_sumstats searches Open GWAS for datasets.
- compute_z computes Z-score from P.
- compute_n computes N for all SNPs from user defined smaple size.
- format_sumstats(ldsc_format=TRUE) ensures sum stats can be fed
directly into LDSC without any additional munging.
- read_sumstats, write_sumstas, and download_vcf functions now
exported.
- format_sumstats(sort_coordinates=TRUE) sorts results by their
genomic coordinates.
- format_sumstats(return_data=TRUE) returns data directly to user.
Can
be returned in either data.table (default), GRanges or VRanges
format using format_sumstats(return_format="granges").
- format_sumstats(N_dropNA=TRUE) (default) drops rows where N is
missing.
- format_sumstats(snp_ids_are_rs_ids=TRUE) (default) Should the SNP
IDs inputted be inferred as RS IDs or some arbitrary ID.
- format_sumstats(write_vcf=TRUE) writes a tabix-indexed VCF file
instead of tabular format.
- format_sumstats(save_path=...) lets users decide where their
results
are saved and what they're named.
- When the save_path indicates it's in tempdir(), message warns users
that these files will be deleted when R session ends.
- Summary of data is given at the beginning and the end of
format_sumstats via report_summary().
- Readability of preview_sumstats() messages improved.
- New checks standard error (SE) must >0 and BETA (and other effect
columns) must not equal 0:
format_sumstats(pos_se=TRUE,effect_columns_nonzero=TRUE)
- Log directory containing all removed SNPs is now available and can
be changed to a different directory by setting:
format_sumstats(log_folder_ind=TRUE,log_folder=tempdir())
- All imputed data can now be identified with a column in the output
using: format_sumstats(imputation_ind=TRUE)
- Users can now input their own mapping file to be used for the
column
header mapping in place of data(sumstatsColHeaders). See
format_sumstats(mapping_file = mapping_file).

Bug fixes

- CHR column now standardised (X and Y caps, no "chr" prefix).
- Allele flipping done on a per-SNP basis (instead of whole-column).
- Allele flipping now includes FRQ column as well as effect columns.
- The effect allele is now interpreted as the A2 allele consistent
with IEU GWAS VCF approach. A1 will always be the reference allele.
- read_vcf upgraded to account for more VCF formats.
- check_n_num now accounts for situations where N is a character
vector and converts to numeric.

                        Changes in version 1.1.1                        

Bug fixes

- Preprint publication citation added.

[muscat](/packages/muscat)
------

                        Changes in version 1.7.2                        

- bug fix in prepSim(): removal of NA coefficients and
  subsetting of the input SCE was previously out of synch

[musicatk](/packages/musicatk)
--------

                 Changes in version 1.3.1 (2021-10-10)                  

- Updated version number to match Bioconductor

                 Changes in version 1.2.1 (2021-08-08)                  

- Fixed bug in fonts for some plots
- Added documentation site generated with pkgdown
- Added Shiny UI for interactive analysis of mutational signatures

[mzR](/packages/mzR)
---

                       Changes in version 2.27.1                        

- Add missing atomic_count_sync.hpp BH file (see PR #248 by vjcitn)

                       Changes in version 2.27.0                        

- New Bioc devel version

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 2.0.0                        

- Major changes to plot_agg_regions().
- Features of plot_agg_regions() and
plot_agg_regions_sample_grouped() merged into one interface.
- Regions now specified using single table.
- Changed plot_regions() default window proportion to 0.
- Changed default theme from theme_bw() to theme_tufte().
- Added Megalodon data import instructions to "Importing Data"
vignette.
- Added scico palette defaults for heatmaps. These are colourblind
friendly.
- Added check for 0 length queries which would cause program to hang
indefinitely.
- Added setters for NanoMethResult attributes methy, samples and
exons.
- Added MDS and PCA plots.
- Added vignette for using external annotation and dimensionality
reduction.
- Added binary thresholding for plot_gene(), plot_region() and
plot_agg_regions().
- Added regions argument to bsseq_to_edger() to calculate aggregate
counts over features rather than per site.

[NanoStringNCTools](/packages/NanoStringNCTools)
-----------------

                 Changes in version 1.1.2 (2020-09-28)                  

- Update license

                 Changes in version 1.1.1 (2020-08-13)                  

- Documentation updates

- Handle new parameters from ggiraph update

                 Changes in version 1.1.0 (2020-05-19)                  

- Initial Bioconductor devel 3.14 version

[NanoTube](/packages/NanoTube)
--------

                       Changes in version 0.99.0                        

- Pre-release version

[nearBynding](/packages/nearBynding)
-----------

                        Changes in version 1.3.3                        

Changes

- update description to accommodate changes in knitr/rmarkdown packages

- update dummy data files to new format so vignette can run

                        Changes in version 1.3.1                        

Changes

- include ability to iteratively acquire and visualize Bkg std error

- new arguments to choose whether to calculate Bkg std error

- accomodate StereoGene version 2.22, esp. seeding of analysis

- better run analysis on track files outside local directory

[netDx](/packages/netDx)
-----

                        Changes in version 1.5.3                        

- Moved RCy3, scater, clusterExperiment and netSmooth to "Suggests" to
  reduce dependency burden

- Sped up vignettes by limiting all to binary classification and
  limiting number of layers

- Removed TL;DR from vignettes as usefulness in question but
  maintainance high.
  Developers notes:

- Added Dockerfile and Github Actions for automated testing

- GHA auto-generates a Docker image with netDx which gets pushed to
  shraddhapai/netdx_devenv

                        Changes in version 1.5.2                        

- Added wrapper functions for ease-of-use. Includes:

- getResults() to plot results of running the predictor

- getPSN() for creating and visualizing integrated PSN

- confusionMatrix() to visualize confusion matrix

- tSNEPlotter() to visualize tSNE of integrated PSN (doesn't require
  Cytoscape)

- Added CITATION file with citations to netDx methods and software
  paper

                        Changes in version 1.5.1                        

- Adding support for Java 16.

- Disabling CNV-based vignette to allow other three vignettes to run
  without causing build timeout on devel system

[netresponse](/packages/netresponse)
-----------

                 Changes in version 1.53.2 (2021-07-28)                 

- rmarkdown added as a dependency to fix Bioc build

[ngsReports](/packages/ngsReports)
----------

                        Changes in version 1.9.3                        

- Bug fixes for importing macs2 logs

                        Changes in version 1.9.2                        

- Bug fixes for later versions of ggplot2

                        Changes in version 1.8.1                        

- Bug fix in .makeSidebar

[normr](/packages/normr)
-----

                       Changes in version 1.18.2                        

- Changed deprecated 'GenomeInfoDb::fetchExtendedChromInfoFromUCSC' to
  'GenomeInfoDb::getChromInfoFromUCSC' in R/methods.R

                       Changes in version 1.18.1                        

- Fixing R 4.1.0 _R_CHECK_LENGTH_1_LOGIC2 error in
  tests/testthat/utils.R:applyMap
  by using inherits() instead of class() to account for hadleyverse

[nullranges](/packages/nullranges)
----------

                        Changes in version 1.0.0                        

- nullranges is released on Bioconductor! the package offers the
creation of null genomic feature sets, either through sampling from
a pool in order to match covariates with a particular focal set, or
via block bootstrapping of features optionally with respect to a
genome segmentation. Critically, nullranges is designed as a modular
package, solely for the purpose of generating null feature sets, and
to be used in conjunction with another package for calculating
overlaps, such as GenomicRanges or plyranges. Let us know your
comments, suggestions or feedback on Bioconductor support site or
through GitHub Issues.

[NxtIRFcore](/packages/NxtIRFcore)
----------

                Changes in version 0.99.12 (2021-10-26)                 

- Fixed bug in MakeSE() and CoordToGR()

                Changes in version 0.99.10 (2021-10-20)                 

- Accounts for when NxtIRFdata cannot fetch data from ExperimentHub

                 Changes in version 0.99.9 (2021-10-20)                 

- Added GetCoverageBins()

- Add warning in IRFinder if coordinate sorted BAM file takes too long
  to run.

- Fixed missing coverage data at both ends of plot track.

                 Changes in version 0.99.8 (2021-10-13)                 

- Fixed memory leak when writing COV files

                 Changes in version 0.99.6 (2021-10-12)                 

- Added GetCoverageRegions() which calculates the mean coverage of each
  region
  in a given GRanges object

- Added BAM2COV() which calculates and creates COV files from BAM files

                 Changes in version 0.99.2 (2021-10-07)                 

- Fixed bug in Find_FASTQ

                 Changes in version 0.99.0 (2021-09-29)                 

- Bioconductor Release

[ODER](/packages/ODER)
----

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[OmnipathR](/packages/OmnipathR)
---------

                        Changes in version 3.2.0                        

- New resource: PrePPI
- Configuration option to completely disable logging to file
- New vignettes: Path reconstruction, and Bioconductor 2021 Workshop
- Many minor bugfixes and improvements

[OncoSimulR](/packages/OncoSimulR)
----------

                 Changes in version 3.1.2 (2021-10-06)                  

- Better test of poset transformation

                 Changes in version 3.1.1 (2021-10-06)                  

- XOR, AND, OR dependencies: plots of DAGs honor all possible values.

- Few miscell minor changes.

[openCyto](/packages/openCyto)
--------

                        Changes in version 3.11                         

Enhancements

- Allow control of mininum number of events for each step of
GatingTemplate #298

Bug Fixes

- Handle a few more edge cases in .improvedMinDensity
- Added a fix to the density estimate used by gate_tautstring

                        Changes in version 3.10                         

API Changes

Simple renaming

- gate_flowClust_1d -> gate_flowclust_1d
- gate_flowClust_2d -> gate_flowclust_2d
- tautStringGate -> gate_tautstring
- templateGen -> gh_generate_template
- add_pop_init -> gs_add_gating_method_init
- add_pop -> gs_add_gating_method
- remove_pop -> gs_remove_gating_method
- get.helperGates -> gs_get_helpergates
- toggle.helperGates -> gt_toggle_helpergates
- delete.helperGates -> gs_delete_helpergates
- getNodes -> gt_get_nodes
- getParent -> gt_get_parent
- getChildren -> gt_get_children
- getGate -> gt_get_gate
- gating -> gt_gating
- registerPlugins -> register_plugins
-
Classes and methods no longer exported

- registerGatingFunction
- gtMethod
- gtPopulation
- polyFunctions
- ppMethod

Bug Fixes

- Some minor fixes to gt_toggle_helpergates
- Fix "positive" argument to gate_mindensity to match doc

[oposSOM](/packages/oposSOM)
-------

                       Changes in version 2.11.1                        

- multiple bugfixes e.g. caused by tsne package

[ORFik](/packages/ORFik)
-----

                       Changes in version 1.13.7                        

SIGNIFICANT USER-VISIBLE CHANGES

- Massive improvement in speed of coveragePerTiling

- Improved p-shifting analysis (also added verbose output)

- Added possible optimization for annotation

- Rewritten vignettes

[orthogene](/packages/orthogene)
---------

                       Changes in version 0.99.9                        

New Features

- Replaced R-CMD GHA with bioc-check GHA.
- Added new badges.

Fixes

- Adjusted vignette yamls to make resulting htmls smaller.

                       Changes in version 0.99.8                        

New Features

- orthogene now supports DelayedArray objects as gene_df input.
- create_background now uses all_genes when all 3 species are the
same.

                       Changes in version 0.99.7                        

New Features

- Added new function create_background.
- Added new function infer_species.
- report_orthologs and convert_orthologs can now handle cases where
input_species is the same as output_species.
- Add internal function get_all_orgs to easily list all organisms
from
different packages.
- Added all_genes method "babelgene".

Fixes

- report_orthologs no longer throws error due to not finding
tar_genes.

                       Changes in version 0.99.6                        

Fixes

- Allow all messages to be suppressed in report_orthologs.

                       Changes in version 0.99.3                        

New Features

- License switched to GPL-3 (to be compliant with Bioc).
- New method "babelgene" added to convert_orthologs.

                       Changes in version 0.99.2                        

- License switched to GPL3 (>=3).

Fixes

- GenomeInfoDbData now required.

                        Changes in version 0.1.0                        

New Features

- orthogene released to Bioconductor.

[padma](/packages/padma)
-----

                  Changes in version 1.3 (2021-10-14)                   

- Minor bug change for lists of data.frame vs matrix for ExperimentList

[pairkat](/packages/pairkat)
-------

                 Changes in version 0.99.0 (2021-09-21)                 

- Submitted to Bioconductor

[pcaExplorer](/packages/pcaExplorer)
-----------

                       Changes in version 2.20.0                        

Other notes

- the tables in the PCA2GO tab panel can be compacted only if they
are
computed via the pca2go function (offline) - at runtime,
limmaquickpca2go is used and no compaction is required
- if an annotation is provided with a column gene_id, these values
are
actually overwriting the rownames (makes the object more robust with
respect to its provenance)

[PCAtools](/packages/PCAtools)
--------

                        Changes in version 2.6.0                        

- added max.overlaps and min.segment.length to provide further control
  over
  connectors. max.overlaps replaces maxoverlapsConnectors, but both can
  still
  be used for legacy purposes

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.7.1 (2021-06-21)                  

- Error due to change in behavior for default axis label in ggplot2
  histogram

- GGplot2 'guide' depreciation warning

[pengls](/packages/pengls)
------

                        Changes in version 0.1.0                        

[PharmacoGx](/packages/PharmacoGx)
----------

                        Changes in version 2.5.3                        

- Added PharmacoSet2 constructor to allow creation of PSets with
updated class definition introducted in BioC 3.13
- The sensitivity slot is now required to be a
TreatmentResponseExperiment
- The molecularProfiles slot is now required to be a
MultiAssayExperiment
- The original constructor and all accessors remain in the package
for
full backwards compatibility

                        Changes in version 2.5.2                        

- Fix: remove 'fdr' item from geneDrugSensitivity return vector

                        Changes in version 2.5.1                        

- Fix: reverted GDSCsmall.rda and CCLEsmall.rda to original data
format; they were accidentally pushed with MultiAssayExperiments in
@molecularProfiles

                        Changes in version 2.5.0                        

- Spring Bioconductor release!

[philr](/packages/philr)
-----

                       Changes in version 1.19.1                        

USER-VISIBLE CHANGES

- Added support for TreeSummarizedExperiment class

- Added support for phyloseq class

- Updated vignette

- Changed main argument name from df to x

INTERNAL

- Implemented philr as S3 method

[PhIPData](/packages/PhIPData)
--------

                 Changes in version 1.1.1 (2021-08-05)                  

- Added propReads() function to calculate the proportion of sample
reads pulled by each peptide.
- Added coercion function to convert a PhIPData object to a DataFrame
(and consequently also data.frames and tibbles).
- Changed the storage of alias, library, and beads-only indicators to
package environment variables rather than global variables.

[PhosR](/packages/PhosR)
-----

                        Changes in version 1.2.1                        

- Included Cell Reporta citation
- rainbow colour palette is used for plotQC function

[PhyloProfile](/packages/PhyloProfile)
------------

                        Changes in version 1.7.8                        

- only mainInput is required in config file

                        Changes in version 1.7.7                        

- use config file for input files, refspec selection and othes
  settings

                        Changes in version 1.6.6                        

- turn off auto sizing for large number of taxa or genes (>=
  10.000)

                        Changes in version 1.6.5                        

- fixed filter for "species relation"

- added midpoint colors

                        Changes in version 1.6.2                        

- increase default font size for profile and domain plot

- make links to DBs: ncbi taxonomy, ncbi protein, uniprot, pfam,
  smart

                        Changes in version 1.6.1                        

- modified taxonomy ranks (ropensci#875)

- improved rank indexing function (ropensci#874)

- improved x-axis label (#116)

[Pigengene](/packages/Pigengene)
---------

                Changes in version 1.19.50 (2021-10-14)                 

New functions

- Neda added the identify.modules(), make.filter(), and
  apply.filter() functions, but not exported them yet.

                Changes in version 1.19.30 (2021-09-08)                 

Bug Fixes

- The averaged.network() function does not have the nodes
  argument in bnlearn Version >=4.7, and thus this argument was
  removed.

                Changes in version 1.19.24 (2021-08-06)                 

New functions

- The get.enriched.pw() function is added.

Bug Fixes

- A bug fix in the gene.mapping() function that used to occur
  when we had multiple output databases.

                Changes in version 1.19.10 (2021-06-25)                 

Changes in existing functions

- The message.if() function can now write the message in a text
  file.

                 Changes in version 1.19.8 (2021-05-25)                 

Bug Fixes

- The C50 plot function seems to have different behavior when the
  number of Labels is 2. Habil reverse the color to fix the
  resulting bug.

[plgem](/packages/plgem)
-----

                       Changes in version 1.65.1                        

- new feature:
  --allow parameter `Prefix' in `run.plgem' to be passed down to
  `plgem.fit' when `writeFiles=TRUE'

- bug fix:
  --only check existence of `fittingEvalFileName' when both
  `fittingEval' and
  `plot.file' are TRUE in `plgem.fit'

- new feature:
  --added parameter `Prefix' to `run.plgem' to be passed down to
  `plgem.write.summary' when `writeFiles=TRUE'

- bug fix:
  --fixed error occurring when running `run.plgem' with `plotFile=TRUE'

[plotgardener](/packages/plotgardener)
------------

                       Changes in version 0.99.11                       

BUG FIXES

o `:` added back to readHic for `strawr` region parsing.
o `plotHicSquare` subsetting fixed for off diagonal regions.

NEW FEATURES

o All Hi-C functions now allow input of remote Hi-C files.

                       Changes in version 0.99.9                        

This package was previously called BentoBox.

                       Changes in version 0.99.0                        

NEW FEATURES

o Version bump to 0.99.0 for Bioconductor package submission.
o `bb_mapColors` function for users to map a vector to a palette
of colors.
o `linecolor` parameter in `bb_plotPairs`, `bb_plotPairsArches`,
and `bb_plotRanges` now accepts a single value, a vector of colors,
a `colorby` object, or the value "fill".

                      Changes in version 0.0.0.14                       

BUG
FIXES

o R version requirement changed to (R >= 4.1.0) for proper plot
placement.

NEW
FEATURES

o `colorby` object now has a `scalePerRegion` parameter to scale
numerical
`colorby` data to the range of data in a plotted genomic region.

                      Changes in version 0.0.0.13                       

SIGNIFICANT USER-VISIBLE CHANGES

o `bb_plotManhattan` `fill` paramete now accepts a single value,
a vector of colors, or a `colorby` object.

                      Changes in version 0.0.0.12                       

SIGNIFICANT USER-VISIBLE CHANGES

o `colorby` constructor now includes optional palette specification.
o `bb_plotPairs`, `bb_plotPairsArches`, and `bb_plotRanges` `fill`
parameter
now accepts a single value, a vector of colors, or a `colorby`
object.

BUG FIXES

o `GInteractions` assembly match checking moved before dataframe
conversion.

                      Changes in version 0.0.0.11                       

SIGNIFICANT USER-VISIBLE CHANGES

o Data moved to `plotgardenerData` package.
o Default genome assembly updated to "hg38".

BUG FIXES

o Streamlined parameter parsing and data reading logic.

                      Changes in version 0.0.0.10                       

NEW FEATURES

o Added unit tests with `testthat`.
o `bb_annoDomains` function addition.
o `bb_plotSignal` vertical orientation.

                       Changes in version 0.0.0.9                       

NEW FEATURES

o Added a `NEWS` file to track changes to the package.

BUG FIXES

o Updated viewport parsing for package `grid` version 4.1.0.

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.1.1 (2021-06-09)                  

- Implemented use of Reference well for control channel

[podkat](/packages/podkat)
------

                       Changes in version 1.25.1                        

- adjusted NAMESPACE to account for changes in BiocGenerics package

                       Changes in version 1.25.0                        

- new branch for Bioconductor 3.14 devel

[polyester](/packages/polyester)
---------

                       Changes in version 1.99.3                        

- NB function now exported

- note that version 1.99.3 on GitHub was version 1.1.0 on Bioconductor.

                       Changes in version 1.99.2                        

- bug fix in fragment generation (last 2 bases of transcript were never
  sequenced)

[preciseTAD](/packages/preciseTAD)
----------

                 Changes in version 1.3.3 (2021-09-28)                  

- The default for MinPts DBSCAN parameter has been changed to 100

                 Changes in version 1.3.2 (2021-09-21)                  

- Vignette edits

                 Changes in version 1.3.1 (2021-07-19)                  

- Change y of x.y.z version number to comply with the release

- Add MD as a maintainer

[proActiv](/packages/proActiv)
--------

                        Changes in version 1.3.3                        

- plotPromoters is deprecated

- Implemented three new functions
  - plotPCA: performs principal component analysis
  - plotHeatmap: heatmap visualization
  - integrateProactiv: integrate different proActiv runs

[procoil](/packages/procoil)
-------

                       Changes in version 2.21.1                        

- fix in plot() method

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

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.33                         

[ProteoDisco](/packages/ProteoDisco)
-----------

                         Changes in version 1.0                         

- Initial release of ProteoDisco (v1.0.0).

[ProtGenerics](/packages/ProtGenerics)
------------

                       Changes in version 1.25.1                        

- add bin, compareChromatograms and compareSpectra

                       Changes in version 1.25.0                        

- Bioc devel (3.14) version bump

[psichomics](/packages/psichomics)
----------

                       Changes in version 1.20.0                        

- Alternative splicing event annotations:
- Support new annotations from VAST-TOOLS for multiple species,
including mouse, zebrafish, fruit fly, chicken, frog, C. elegans
and A. thaliana
- Automatically create cache directory if downloading splicing
annotations for the first time
- Gene, transcript and protein annotation (visual interface):
- Automatically set species/genome based on selected annotation
- Improve species and genome selection

                       Changes in version 1.18.6                        

ShinyProxy support

- psichomics(): add argument shinyproxy; when set to TRUE, change set
of options to viably run in ShinyProxy
- Avoid automatically closing the app
- Replace custom file browsers with shiny's versions
- Fix issues with progress bar
- Include ShinyBS JavaScript library
- Upload files using default file browser input and allow to upload a
ZIP folder instead of a folder

Bug fixes

- Fix issues with Shiny 1.7.0:
- Change icon names as required by newest versions of font-awesome
- Avoid modifying Shiny tags using generic positions
- Load recount data (visual interface):
- Improve responsiveness of project selection
- Splicing annotation (visual interface):
- Load annotation only when opening the splicing quantification
tab
- Automatically select hg38 if using GTEx v8 data
- PCA plot (visual interface):
- Automatically plot PCA after calculating PCA scores
- Copy-edit text
- Fix specific errors related with PCA analysis of only one group
- Diagrams of alternative splicing events:
- Fix wrong coloring of reference exon used for AFE and A5SS
events
- Transcript plot:
- Orange region (the reference exon) is now on top of blue region

                       Changes in version 1.18.5                        

- Diagrams of alternative splicing events:
- Fix wrong coloring of reference exon used for AFE and A5SS
events
- Transcript plot:
- Avoid alternative regions from overlapping
- Fix loading twice when selecting a new event (visual interface)

                       Changes in version 1.18.4                        

- psichomics(): fix visual interface not launching
- getGtexReleases() not properly retrieving whether future GTEx
releases (9 and higher) are available
- Remove warning related with TCGA data when MD5 checks fail

                       Changes in version 1.18.3                        

- plotSplicingEvent(): avoid opening browser window in
non-interactive
contexts
- Fix Bioconductor build report's timeout when creating vignettes on
Windows

                       Changes in version 1.18.2                        

- Fix issues with unit tests in Bioconductor

                       Changes in version 1.18.1                        

- Fix issues with unit tests in Bioconductor
- Improve Docker images:
- Simplify Dockerfile
- Store Docker images in GitHub Container Registry
- Automatically build latest Docker image based on last release
version
- Improve README with install instructions for GitHub and Docker
- Minor copy-editing of user-provided data tutorial
- Fix minor spelling issues

[ptairMS](/packages/ptairMS)
-------

                 Changes in version 1.1.4 (2021-09-23)                  

- removal of the loaded prtset data to avoid any local path problems

                 Changes in version 1.1.1 (2021-09-13)                  

- Added graphical interface

[PureCN](/packages/PureCN)
------

                        Changes in version 2.0.0                        

NEW FEATURES

- Report median absolute pairwise difference (MAPD) of tumor vs normal
  log2
  ratios in runAbsoluteCN

- Improved mapping bias estimates: variants with insufficient
  information
  for position-specific fits (default 3-6 heterozygous variants)
  are clustered and assigned to the most similar fit

- Make Cosmic.CNT INFO field name customizable

SIGNIFICANT USER-VISIBLE CHANGES

- Cleanup of naming of command line arguments (will throw lots of
  deprecated
  warnings, but was long overdue)

- More robust alignment of on- and off-target tumor vs normal log2
  ratios.
  Ratios are shifted so that median difference of neighboring
  on/off-target
  pairs is 0. This should fix spurious segments consisting of only on-
  or
  off-target regions in high quality samples where those minor off-sets
  sometimes exceeded the noise.

- Added min.variants argument to runAbsoluteCN

- Added PureCN version to runAbsoluteCN results object (ret$version)

- Addressed observed over-segmentations in very clean data:
  - Do not attempt two-step segmentation in PSCBS when off-target noise
  is
  still very small (< 0.15, min.logr.sdev in runAbsoluteCN)
  - Increase automatically determined undo.SD in all segmentation
  functions
  when noise is very small (< min.logr.sdev)
  - min.logr.sdev is now accessible in PureCN.R via --min-logr-sdev

- Added pairwise sample distances to normalDB output object helpful for
  finding noisy samples or batches in normal databases

- Do not error out readCurationFile when CSV is missing and directory
  is not writable when re-generating it (#196)

- Add segmentation parameters as attributes to segmentation data.frame

- Added min.betafit.rho and max.betafit.rho to calculateMappingBias*

- Made --normal_panel in PureCN.R defunct

- Added GATK/Picard header with sequence lengths to interval file,
  added readIntervalFile function to parse it

BUGFIXES

- Fix for crash when --normal_panel in NormalDB.R contained no variants
  (#180).

- Fix for crash when rtracklayer failed to parse --infile in
  FilterCallableLoci.R (#182)

- More robust parsing of VCF with missing GT field (#184)

- Fix for bug and crash when mapping bias RDS file contains variants
  with
  multiple alt alleles (#184)

- Added missing dependency 'markdown'

- Fix for crash when only a small number of off-target intervals pass
  filters (#190)

- Fix for crash when PSCBS segmentation was selected without VCF file
  (#190)

- Fix for crash when Hclust segmentation was selected without
  segmentation
  file (#190)

- Fix for crashes when not many variant pass filters (#192, #195)

- Fix for crash when provided segmentation does not have chromosomes
  in common with VCF (#192) or does not provide all chromosomes present
  in
  the coverage file (#192)

[qcmetrics](/packages/qcmetrics)
---------

                        Changes in version 1.31                         

qcmetrics 1.31.1

- Suggest rmarkdown to fix build error

[QDNAseq](/packages/QDNAseq)
-------

                 Changes in version 1.29.6 (2021-10-23)                 

BUG FIXES

- segmentBins() would *report* the sample names as "NA" in output
  messages
  if the sample name contained hyphens, or other symbols automatically
  replaced by data.frame(..., check.names = TRUE).  This was a harmless
  bug.

                 Changes in version 1.29.5 (2021-10-20)                 

PERFORMANCE

- All internal row and column-based matrixStats calls now avoids
  overhead
  from handling row and column names.

MISCELLANEOUS

- Moved 'future' from Imports to Suggests.

                 Changes in version 1.29.4 (2021-10-16)                 

DOCUMENTATION

- Vignette now illustrate parallelization using the 'multisession'
  future
  strategy, instead of the deprecated 'multiprocess' strategy.

DEPRECATION AND DEFUNCT

- Argument 'seeds' of segmentBins() is defunct.  It has been deprecated
  and
  ignored since QDNAseq 1.21.3 (September 2019).

                 Changes in version 1.29.3 (2021-10-04)                 

NEW FEATURES

- Now argument 'logTransform' of exportBins() is ignored if 'type' =
  "calls".

- Now exportBins() returns the pathname to the files written.

SOFTWARE QUALITY

- Test code coverage was increased from 42% to 52%.

- Add package test for exportBins().

BUG FIXES

- exportBins(fit, format = "seg", ...) and format = "vcf" would merge
  segments with equal copy-number calls if they were interweaved with
  copy-neutral segments.

- exportBins(fit, format = "seg", ...) and format = "vcf" produced an
  obscure
  error with messages "Error in dimnames(x) <- dn : length of
  'dimnames' [2]
  not equal to array extent" for samples with no copy-number
  abberations.

- exportBins(fit, format = "seg", file = ...) and format = "vcf" did
  not
  respect argument 'file' but instead wrote files of its own names to
  the
  current working directory.

- exportBins() would corrupt option 'scipen'.  Now it is left
  unchanged.

KNOWN ISSUES

- callBins() produces warnings on "Recycling array of length 1 in
  vector-
  array arithmetic is deprecated. Use c() or as.vector() instead." in
  R (>= 3.4.0).  This is a problem in the package 'CGHcall' dependency
  and
  is something that needs to be fixed there.  For further details,
  please
  see https://github.com/tgac-vumc/CGHcall/issues/2.

                 Changes in version 1.29.2 (2021-09-22)                 

SOFTWARE QUALITY

- Test code coverage was increased from 32% to 39%.

- Added package tests for binReadCounts().

BUG FIXES

- binReadCounts() would fail when specifying argument 'chunkSize'.  The
  fix
  was to require 'future' package version 1.22.1 or newer.

                 Changes in version 1.29.1 (2021-08-26)                 

SOFTWARE QUALITY

- Add package test for binReadCounts().

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.3.0                        

QFeatures 1.3.6

- New feat3 example data to demonstrate and test more complex
AssayLinks structure.
- Improved the plot,QFeautres function to avoid cluttering of nodes.
- Adapted the visualization vignette using feat3.

QFeatures 1.3.5

- Add plot,QFeatures and visualisation vignette.

QFeatures 1.3.4

- Fixed bug that produced invalid AssayLinks when using filterNA.

QFeatures 1.3.3

- Improved validity checks on AssayLinks
- Fixed the subsetting of AssayLinks to ensure consistent data

QFeatures 1.3.2

- Add logo to package
- Fix class coercion error (see #b9ce7f1e9)

QFeatures 1.3.1

- Added rbindRowData: a function to select variables in the rowData
and bind it in a single DataFrame
- Added rowData<-: this new method replaces replaceRowDataCols to
offer a more standardize functionality.
- Added a new section in the QFeatures vignette to expand on how to
manipulate the metadata within a QFeatures object

QFeatures 1.3.0

- New devel version (Bioc 3.14)

[qsmooth](/packages/qsmooth)
-------

                 Changes in version 1.9.2 (2021-09-01)                  

- Added Hmisc to Imports

                 Changes in version 1.9.1 (2021-06-24)                  

- Added qsmoothGC function (Contributed from Koen Van den Berge)

[QuasR](/packages/QuasR)
-----

                       Changes in version 1.34.0                        

USER-VISIBLE CHANGES

- removed automatic downloading and installation of BSgenome references

- added option to parallelize qExportWig

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.18.0                        

New features

- Sparse matrices of dgCMatrix type can now be coerced to
RaggedExperiment when rownames are coercible to GRanges

Bug fixes and minor improvements

- 'counts' set as the default name for the values in mcols after
coercion from dgCMatrix

[ramr](/packages/ramr)
----

                 Changes in version 1.1.2 (2021-05-28)                  

- reviewers' suggestions were implemented, docs updated, typos fixed

- new methods for simulation of AMR and test data sets

- NULL as a default for data.samples (to use all)

- doRNG is used to ensure reproducibility during parallel computing

- a couple of new defaults for previously required parameters for easy
  usage

                 Changes in version 1.1.0 (2021-05-21)                  

- released at bioconductor


[rawrr](/packages/rawrr)
-----

                  Changes in version 1.1 (2021-05-31)                   

- Improved error handling of system2 call in .rawrrSystem2Source
  by logging stdout and stderr and make them available from the R
  console.

- Added helper function .checkReaderFunctions.

- Use pipe |> in vignette.

[RCy3](/packages/RCy3)
----

                       Changes in version 2.14.0                        

- Cleaned up dependencies, dramatically reducing RCy3 package
  installation time

- New functions:
  - add and update Annotations
  - uniqueList parameter added to edgeNameToEdgeSUID and
  nodeNameToNodeSUID, #139

- Bug fixes:
  - loadTableData now works with tibbles, #143
  - sandboxSendTo now works with cys and png, #138
  - .verifySupportedVersions fixed comparisons, #152

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.7.5 (2021-09-28)                  

- Fixed documentation of "plot_heatmap"

                 Changes in version 1.7.4 (2021-09-24)                  

- Added new plotting function "plot_heatmap"

                 Changes in version 1.7.3 (2021-09-14)                  

- Updated vignette to introduce the "open_reactome" command

                 Changes in version 1.7.2 (2021-09-09)                  

- Fixed timeout issue during large requests in
  perform_reactome_analysis

                 Changes in version 1.7.1 (2021-06-09)                  

- Fixed bug in scRNA-seq vignette.

[recount](/packages/recount)
-------

                       Changes in version 1.19.2                        

BUG FIXES

- Fix a bug in geo_info() for reading files on Windows where a
trailing \r was added to all variables.
- Avoid the implicit list embedding of S4 objects is deprecated
warning that was noted at
https://github.com/leekgroup/recount/runs/3286046827?check_suite_focus=true#step:20:1417.

[recount3](/packages/recount3)
--------

                        Changes in version 1.3.9                        

BUG FIXES

- Resolved https://github.com/LieberInstitute/recount3/issues/7.

                        Changes in version 1.3.7                        

NEW FEATURES

- Added the create_hub() function for creating UCSC track hub
configuration files for using the UCSC Genome Browser to explore the
recount3 BigWig base-pair coverage files.

                        Changes in version 1.3.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- rowRanges(rse_gene)$score is now rowRanges(rse_gene)$bp_length to
make it easier to use recount::getTPM() and recount::getRPKM() with
recount3 objects. Resolves
https://github.com/LieberInstitute/recount3/issues/4.

BUG FIXES

- Updated read_metadata() based on
https://github.com/LieberInstitute/recount3/issues/5. The empty
metadata will be dropped with a warning in situations like that.

[RegEnrich](/packages/RegEnrich)
---------

                        Changes in version 1.3.1                        

- Change the email from w.tao-2@umcturecht.nl to
weiyangtao1513@gmail.com

[rfaRm](/packages/rfaRm)
-----

                 Changes in version 1.5.3 (2021-08-04)                  

- Updated formatting of consensus secondary structure queries to match
  changes made by Rfam

[RGMQL](/packages/RGMQL)
-----

                       Changes in version 1.12.2                        

NEW FEATURES

- None

SIGNIFICANT USER-VISIBLE CHANGES

- None

DEPRECATED AND DEFUNCT

- None

BUG FIXES

- changed implementation show_all_metadata() for better preformance

                       Changes in version 1.12.1                        

NEW FEATURES

- Added new function show_all_metadata()

SIGNIFICANT USER-VISIBLE CHANGES

- removed is_GMQL from read_gmql function
  The entire dataset must have the right folder structure in order to
  works correctly <dataset_name> ---> <files>

- Swap order of arguments 'dir_out' and 'name' of the collect()
  function so now the latter comes before the former.

DEPRECATED AND DEFUNCT

- None

BUG FIXES

- fixed the remote processing

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.38.0                        

NEW FEATURES

- Added support for reading attributes where the datatype is
  either a 64-bit or unsigned 32-bit integer.

- Added many functions for working with file creation property
  lists.  (Thanks to @ilia-kats for the contribution,
  https://github.com/grimbough/rhdf5/pull/95)

- Added support for variable length and UTF-8 encoded string
  datasets.  (Thanks to Aaron Lun @LTLA for the contribution,
  https://github.com/grimbough/rhdf5/pull/88)

CHANGES

- Documentation switched to roxygen2

BUG FIXES

- h5createDataset() now prints a warning if a chunk dimension
  exceeds the maximum size of that dimension and automatically
  sets the corresponding chunk dimension to the maxiumum.
  (Thanks to Eric Kernfeld @ekernf01 for the report,
  https://github.com/grimbough/rhdf5/issues/97)

[Rhdf5lib](/packages/Rhdf5lib)
--------

                        Changes in version 1.16                         

Bug fixes

- AR and RANLIB programs used to compile R are now also used to compile
  the HDF5 library. This resolves issue when the default versions found
  on a
  system are incompatible with options used to build R.
  (thanks to @miesav, https://github.com/grimbough/Rhdf5lib/pull/41)

- Fixed issue in Windows installation introduced by upstream changes
  to libcurl distributed by rwinlibs.
  (https://github.com/grimbough/Rhdf5lib/pull/42)

[RLSeq](/packages/RLSeq)
-----

                       Changes in version 0.99.11                       

Package

- Base class
- RLRanges stores R-loop data and RLSeq results.
- Workflow -The wrapper function RLSeq() performs all analysis steps.
However, individual functions are also accessible.
- Model
- The model for predicting sample quality label is updated to the
latest version, incorporating 231 samples in training. This
model showed high accuracy (.9304) on a test set of 115 samples.

Status

Pre-release. This version is prior to the first official release
(v1.0.0), anticipated in Oct. 2021.

[rmspc](/packages/rmspc)
-----

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[RNAmodR](/packages/RNAmodR)
-------

                 Changes in version 1.7.1 (2021-07-27)                  

- internal bugfix

[rols](/packages/rols)
----

                        Changes in version 2.21                         

CHANGES IN VERSION 2.21.1

- Fix failing unit test

CHANGES IN VERSION 2.21.0

- New devel version for Bioc 3.14

[ROSeq](/packages/ROSeq)
-----

                       Changes in version 1.99.01                       

- no bug, but removed extra files from man folder

[RPA](/packages/RPA)
---

                 Changes in version 1.49.1 (2021-07-28)                 

- rmarkdown added to dependency to fix bioc builds

[rpx](/packages/rpx)
---

                         Changes in version 2.1                         

rpx 2.1.12

- Annotate additional experiments as returning errors (see issues for
full list) <2021-10-07 Thu>

rpx 2.1.11

- Annotate PXD012095 as returning an error (see #12) <2021-10-04 Mon>

rpx 2.1.10

- Caching PRIDE sitemap with all PXD projects.
- New pxinstruments(), pxptms() and pxprotocols() accessors.

rpx 2.1.9

- Internal function to tally local project vs full PX.

rpx 2.1.8

- New PXDataset2 class with richer interface and more stable data
downloading functions. PXDataset and PXDataset2 work transparently
and PXDataset2 is now default.

- Add deprecation notice in PXDataset() constructor.

rpx 2.1.7

- Improve object creating printout.
- Improve pxCachedProjects() output.

rpx 2.1.6

- Update installation instruction in README.md.
- Cache location is now also stored inside the PXDataset object.
- New pxCacheInfo() function and use it to show caching info in
show,PXDataset.

rpx 2.1.5

- Improve documentation.

- Check for cached PXDataset object validity.

rpx 2.1.4

- Fix bug in PXDdataset internal data storage.

- New pxCachedProjects() function that return the cached projects.

- Fixes and improvements in the documentation.

rpx 2.1.3

- PXDatasets are also cached upon creation and retrieved from cache
next time they are generated.

- cache is now returned by rpxCache().

rpx 2.1.2

- New feature: PXDataset objects now query and store data (ref, tax,
files, url) when generated instead of fetching these on the fly
every time. This new feature has been added to circumvent the issues
with data access (see #5).

rpx 2.1.1

- pxannouced() paused (see #7).

[rrvgo](/packages/rrvgo)
-----

                        Changes in version 1.5.4                        

- `scatterPlot()` doesn't warn anymore that we're using a deprecated
  parm to remove the guide

                        Changes in version 1.5.1                        

- `calculateSimMatrix()` now allows using arbitrary keys from Orgdb
  packages. Credit: illumination-k. Thanks!

[Rsamtools](/packages/Rsamtools)
---------

                        Changes in version 2.10                         

DEPRECATED AND DEFUNCT

- (v 2.9.1) Deprecate applyPileups() in favor of pileup().

[RTCGAToolbox](/packages/RTCGAToolbox)
------------

                       Changes in version 2.24.0                        

Bug fixes and minor improvements

- The deprecated functionality vignette moved to the
vignettes/analysis/ folder

[rWikiPathways](/packages/rWikiPathways)
-------------

                       Changes in version 1.14.0                        

- None

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.32.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Subsetting a DataFrame object by row names no longer uses partial
  matching.

[scanMiR](/packages/scanMiR)
-------

                Changes in version 0.99.25 (2021-06-25)                 

- fixed a recent bug preventing the recognition of many slicing sites

[scater](/packages/scater)
------

                       Changes in version 1.22.0                        

- Rename colour_columns_by in plotHeatmap to color_columns_by to
  match other arguments.

- Add color_rows_by and row_annotation_colors arguments to
  plotHeatmap, similar to analogous column arguments.

- Change text_by annotations in plotReducedDim to use
  geom_text_repel from ggrepel.

[scDblFinder](/packages/scDblFinder)
-----------

                 Changes in version 1.7.3 (2021-07-26)                  

- scDblFinder now includes both cluster-based and random modes for
  artificial doublet generation

- thresholding has been streamlined

- default parameters have been optimized using benchmark datasets

- added the `directDblClassification` method

[scone](/packages/scone)
-----

                 Changes in version 1.17.1 (2021-07-21)                 

- Improved PsiNorm vignette.

- Fix bug that prevented PSINORM_FN() to be exported.

[SCOPE](/packages/SCOPE)
-----

                        Changes in version 1.5.2                        

- remove loading warnings

                        Changes in version 1.5.1                        

- update vignette

[scp](/packages/scp)
---

                        Changes in version 1.3.3                        

- docs: included QFeatures plot in the vignette
- docs: created a vignette about advanced usage of scp

                        Changes in version 1.3.2                        

- feat: computeSCR now allows for user supplied function that will
summarize the values from multiple samples and multiple carrier.
- docs: used more standard variable names in scp vignette.
- docs: created a QFeatures recap vignette

                        Changes in version 1.3.1                        

- refactor: deprecated rowDataToDF. This function is now replaced by
QFeatures::rbindRowData.

                         Changes in version 1.3                         

                        Changes in version 1.3.0                        

- New devel (Bioc 3.14)

[scPCA](/packages/scPCA)
-----

                 Changes in version 1.7.3 (2021-10-07)                  

- Removing more tests attempting to verify that parallelized outputs
perfectly match their serial counterparts.

                 Changes in version 1.7.2 (2021-09-15)                  

- Removing tests checking that sequential and parallel calls to
scPCA() produce identical outputs when BiocParallel's SerialParam()
is used. This due to new handing of random number generation in
BiocParallel version 1.28.

[scReClassify](/packages/scReClassify)
------------

                       Changes in version 0.99.8                        

- Included citation for the package

                       Changes in version 0.99.7                        

- Removed global assignment in multiAdaSampling.

                       Changes in version 0.99.6                        

- Revised further to address comments/feedbacks

                       Changes in version 0.99.5                        

- Fixed NEWS to match the version bump

                       Changes in version 0.99.4                        

- Minor change to example in matPC function

                       Changes in version 0.99.3                        

- Correct for any spelling errors in all documentations
- Revising the package to address the comments from Bioconductor
review
- Cleaned the GSE87795 subset data to a SingleCellExperiment object

                       Changes in version 0.99.2                        

- version bump with submission to Bioconductor

                       Changes in version 0.99.1                        

- Cleaning repository to pass BiocCheck

                       Changes in version 0.99.0                        

- Initial submission for Bioconductor
- Code formats/documentations have been revised to meet Bioconductor
requirements

                        Changes in version 0.1.1                        

- Significant updates and addition to the documentations
- Major changes to code writing style to comply with BioC
- Updated example data

[scShapes](/packages/scShapes)
--------

                        Changes in version 0.1.0                        

- New package scShapes, for identifying and modeling distribution
shapes of single-cell RNA sequencing data.

[scTensor](/packages/scTensor)
--------

                        Changes in version 2.4.0                        

- Regularizer parameter (L2_A/L1_A) was added in cellCellDecomp()
  ("ntd", "ntd2").

- Multilinear CX Decompotision was added in cellCellDecomp() ("cx").

- convertNCBIGeneID is removed.

- The vignettes were modified.

- Support of LRBase.XXX.eg.db-type packages is completely deprecated

[segmenter](/packages/segmenter)
---------

                       Changes in version 0.99.00                       

- Added a NEWS.md file to track changes to the package.

[selectKSigs](/packages/selectKSigs)
-----------

                        Changes in version 1.5.1                        

- Prepare for release

                        Changes in version 1.4.1                        

- Prepare for release

[seqcombo](/packages/seqcombo)
--------

                       Changes in version 1.15.1                        

- import yulab.utils (2021-08-20, Fri)
- mv seqdiff and simplot to ggmsa package

[sigFeature](/packages/sigFeature)
----------

                       Changes in version 1.11.2                        

New functionality

- 
  Changes in version 1.11.2 (2021-09-14) Update the DESCRIPTION
  file, vignettes file and the NEWS file.

[signatureSearch](/packages/signatureSearch)
---------------

                 Changes in version 1.7.3 (2021-08-24)                  

- Improved get_targets function by supporting different output
format.

                 Changes in version 1.7.2 (2021-06-14)                  

- Supported automatic downloads of ExperimentHub cached files.

[signeR](/packages/signeR)
------

                       Changes in version 1.19.1                        

- migrate to PMCMRplus package

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.3.2 (2021-10-24)                  

- Added summary table into the cellQC report
- Improved formatting in QC report
- Added functions getDEGTopTable() & plotBatchCorrCompare()
- Other refactors and bug fixes

                 Changes in version 2.3.1 (2021-10-15)                  

- Several bug fixes

                 Changes in version 2.2.2 (2021-10-10)                  

- Several enhancements, refactors, and bug fixes to the UI
- Refactor documentation and pkgdown site
- Added tutorials for R console analysis
- Updates to the UMAP generation in the SCTK-QC pipeline
- Addition of VAM to Pathway prediction tab
- Bug fix to the mitochondrial gene set functions

[sitePath](/packages/sitePath)
--------

                        Changes in version 1.9.4                        

- Fix: special case when tree root has more than one lineage path

                        Changes in version 1.9.3                        

- Fix: invalid parallel mutations at divergent node.

- Update DESCRIPTION, README and vignettes.

                        Changes in version 1.9.2                        

- Allow partially plot 'lineagePath'.

- Improved 'plotMutSites' function for 'lineagePath'.

                        Changes in version 1.9.1                        

- Add 'useSites' argument to 'setSiteNumbering' function.

- First 'stable' path as default 'lineagePath'.

- Enable plot functions for 'parallelSites'.

[snifter](/packages/snifter)
-------

                 Changes in version 1.4.0 (2021-09-09)                  

- Add pca, partial_pca, pca_center, pca_scale, pca_dims arguments in
line with Rtsne::Rtsne

[sparrow](/packages/sparrow)
-------

                         Changes in version 1.0                         

Enhancements

- Released to Bioconductor
- Adds support for use of BiocSet as a means by which users can bring
their genesets to -- or take them from -- sparrow.
- Improvements to the corplot() functionality contributed by by
Arkadiusz Gladki (@gladki). Users can specify the size of the text
reported in the bottom half of the pair plot, and spurious/annoying
warnings that were produced after a a totally valid call are no
longer produced.

Breaking Changes from Pre-release

- First two parameters in ora() function have been swapped so that
the
first parameter (x) is the object (data.frame) to run an over
representation analysis against, and the second parameter is the
GeneSetDb.
- scoreSingleSamples no longer drops features in y that are not found
in the GeneSetDb used for scoring. This was changed so that gsva and
ssGSEA scores match the scores produced by a normal GSVA::gsva call.
You can set the drop.unconformed = TRUE to retain the older
behavior.

[SpatialDecon](/packages/SpatialDecon)
------------

                 Changes in version 1.3.0 (2021-09-28)                  

- added S4 wrappers for Seurat and GeoMxSet objects

- added custom profile matrix generation from single cell data

- added ~75 profile matrices avaliable to download for human and mouse

[SpatialExperiment](/packages/SpatialExperiment)
-----------------

                 Changes in version 1.3.2 (2021-07-27)                  

- spatialData moved from colData to int_colData

- restructuring of vignette and added imgData section

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 1.99.0 (2021-10-14)                 

- Implemented overlaying feature: overlay template images in raster
  format with SHMs, where charcoal and transparency options are
  provided.

- Implemented Spatial Single Cell functionality: co-visualize single
  cells and bulk tissues by placing single-cell embedding plots (PCA,
  tSNE, UMAP) and SHMs side by side, cell cluster assignments are
  defined in the app or provided by users, in dimensionality plots
  shapes are not restricted to 6, etc.

- Spatial enrichment was synced to R command line.

- Implemented Sigle and Multiple search mode for gene IDs.

[spatzie](/packages/spatzie)
-------

                       Changes in version 0.99.7                        

- added direct support for GenomicInteractions objects

- improved vignettes

- removed biomaRt requests

- added support for R 4.1

- added 'Transcription' Bioconductor Views annotation

                       Changes in version 0.99.0                        

- initial pre-release

[Spectra](/packages/Spectra)
-------

                         Changes in version 1.3                         

Changes in 1.3.11

- Fix error message in setBackend (issue #217).

Changes in 1.3.10

- Fix bug in plotSpectra and plotSpectraMirror that would cause an
error if the number of peaks in a spectrum was 1 and labels were
provided.

Changes in 1.3.9

- New features: joinSpectraData() now check for duplicated keys in x
(throws an error) and y (thows a warning).

Changes in 1.3.8

- New features: plotMzDelta() function to M/Z delta QC (ported from
MSnbase).

Changes in 1.3.7

- Add fix from MSnbase (issue #170) to Spectra: on macOS require
reading also the spectrum header before reading the peaks data.

Changes in 1.3.6

- Documentation updates for combineSpectra and combinePeaks.

Changes in 1.3.5

- filterMzValues supports also removing peaks matching specified m/z
values (issue #209).

Changes in 1.3.4

- Add list of additional R packages and repositories providing
MsBackend backends to the vignette.

Changes in 1.3.3

- Move generics for bin and compareSpectra to ProtGenerics.

Changes in 1.3.2

- Add parameter f to filterPrecursorScan to fix issue #194.

Changes in 1.3.1

- Add estimatePrecursorIntensity function (issue #202).

[spiky](/packages/spiky)
-----

                        Changes in version 0.1.0                        

- Initial addition of functions and tests

[splatter](/packages/splatter)
--------

                 Changes in version 1.18.0 (2021-10-27)                 

- 
  Updates to the splatPop simulation (from Christina Azodi)
  
  • Added functionality to simulate directly from empirical
  values
  
  • Added eqtl.coreg parameter to splatPop
  
  • Fixed a bug where too many cells were simulated in splatPop
  with multiple batches
  
  • Fixed duplicate cell names in splatPopSimulate

- 
  Improved checks for group.prob in SplatParams

- 
  Automatically rescale group.prob during setting if it doesn't
  sum to 1

[SplicingFactory](/packages/SplicingFactory)
---------------

                        Changes in version 1.0.3                        

- Bug fix for example data.
- Other minor changes.

                        Changes in version 1.0.2                        

- Bug fix for single sample analysis.

                        Changes in version 1.0.1                        

- Updates to URL, email and citation.

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

[struct](/packages/struct)
------

                        Changes in version 1.5.3                        

- fix variable_meta assignment

                        Changes in version 1.5.2                        

- use ontology slots instead of STATO

[structToolbox](/packages/structToolbox)
-------------

                        Changes in version 1.5.7                        

- improve NA handling in fold_change computations

                        Changes in version 1.5.5                        

- change t-test outputs to data.frame

- fix NA bug in feature_boxplot chart

- use new ontology system in place of stato

                        Changes in version 1.5.2                        

- add outputs to auto-generated documentation

- fix fold change threshold using median (#56)

- add fold change using means (#57)

- HSD param "unbalanced" is no longer ignored

                        Changes in version 1.5.1                        

- fixed broken paired tests

                        Changes in version 1.4.2                        

- fix fold change threshold using median (#56)

- HSD param "unbalanced" is no longer ignored

                        Changes in version 1.4.1                        

- fixed broken paired tests

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.24.0                        

NEW FEATURES

- Add 'checkDimnames' argument to SummarizedExperiment() constructor
  function

- Add showAsCell() method for SummarizedExperiment objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Check the assay dimnames at SummarizedExperiment construction time:
  The SummarizedExperiment() constructor function now raises an error
  if one of the supplied assays has rownames and/or colnames that don't
  match those of the SummarizedExperiment object to construct.

[surfaltr](/packages/surfaltr)
--------

                 Changes in version 0.99.6 (2021-10-25)                 

- Fixed text for consistency in package name

                 Changes in version 0.99.5 (2021-10-21)                 

- Fixed text for consistency in package name

                 Changes in version 0.99.4 (2021-10-12)                 

                 Changes in version 0.99.3 (2021-10-11)                 

- Fixed .gitignore file to fix build error

                 Changes in version 0.99.2 (2021-09-14)                 

- Fixed bug in script

                 Changes in version 0.99.0 (2021-08-24)                 

- Submitted to Bioconductor

[svaNUMT](/packages/svaNUMT)
-------

                 Changes in version 0.99.0 (2021-04-28)                 

- Submitted to BioConductor

[svaRetro](/packages/svaRetro)
--------

                 Changes in version 0.99.0 (2021-04-28)                 

- Submitted to BioConductor

[synapsis](/packages/synapsis)
--------

               Changes in version 2021-07-02 (2021-07-02)               

Changes:

- Added this news file

[synapter](/packages/synapter)
--------

                        Changes in version 2.17                         

Changes in version 2.17.3

- put commented code chunk back in

Changes in version 2.17.2

- fix error in vignette

Changes in version 2.17.1

- Update Laurent's email address

[SynExtend](/packages/SynExtend)
---------

                        Changes in version 1.5.4                        

- Added the function SubSetPairs that allows for easy trimming of
predicted pairs based on conflicting predictions and / or prediction
statistics.
- Added the function EstimageGenomeRearrangements that generates
rearrangement scenarios of large scale genomic events using the
double cut and join model.

                        Changes in version 1.5.3                        

- Added the function SequenceSimilarity and made improvements to
runtime in DisjointSet.

                        Changes in version 1.4.1                        

- Fixed a small bug in consensus scores in PairSummaries where
features facing on different strands had their score computed
incorrectly.

[systemPipeShiny](/packages/systemPipeShiny)
---------------

                       Changes in version 1.3.15                        

New Feature

- In workflow module, workflow designer (step 3), added two new
parameter when creating a new step, mandatory and place of
execution. These are new features added in systemPipeR 1.27.27.

Minor Change

- Bump version requirements of spsComps, systemPipeR, systemPipeRdata

- In global.R, now use spsOption with the .list argument to set up
options instead of the base options function.

- Replace includeMardown() by markdown(readLines()) so we don't need
additional {markdown} package as dependency.

Bug Fix

- Fix links and image urls that were not working or changed.

                       Changes in version 1.3.10                        

New Feature

- Add code display buttons to most plots that will show code to
reproduce the plot.

- Add two args buttonType and placeholder to dynamicFile, now users
can specify what bootstrap color the button is and use placeholder
to specify initial text on the upload bar.

- Enhanced the original shiny fileInput, now users can also
specify icon and button bootstrap colors for "server" mode in
dynamicFile.

Major Change

- Redesign of a few steps in Workflow module. The new version of
{systemPipeR} fundamentally changed how the workflow will be run. To
sync to this new version, WF module has to been redesigned. Major
change happens on workflow step selection. This requires users to
install systemPipeR > 1.27.10

- New methods to initiate the WF project

- New workflow plot

- New step selection mechanism

- New step editing functionalities

Minor Change

- For RNAseq module, the dendrogram plot library changed from
{ggtree}
package to {ape}. {ggtree} is not very compatible with Shiny under
current version. Plot cannot be created, always error, but no error
outside Shiny. An issue has submitted to Shiny on Github. We may
switch back to ggtree when this is fixed.
- Small UI optimization for RNAseq module.
- Fixed some typo in different tabs.

Bug Fix

- #85 fix dynamicFile icon not working

- Also add some icon validation code

- Fix the admin server tabs get loaded twice. Added a flag to prevent
this from happening.

                        Changes in version 1.3.0                        

- Update version number to 1.3.0 per Bioconductor regulation.

[tanggle](/packages/tanggle)
-------

                        Changes in version 0.99                         

NEW FEATURES

- initial release. The main functions are ggsplitnet and ggevonet to
  visualize split (or implicit) networks (unrooted, undirected) and
  explicit
  networks (rooted, directed) with reticulations.

[TAPseq](/packages/TAPseq)
------

                 Changes in version 1.5.1 (2021-08-27)                  

- Fixed a bug in createPrimerTrack()

- Added citation

[TargetSearch](/packages/TargetSearch)
------------

                       Changes in version 1.50.0                        

NEW FEATURES

- FindAllPeaks: Allow for asymmetric RT deviations. Formerly, the
  window
  search parameter was plus o minus a tolerance; now it can be
  different on either
  side of the expected RT.

- ncdf4_convert_from_path: New flag to convert CDF files recursively.

- checkRimLim: show multiple samples at the same time, as opposed to a
  single sample in previous versions.

BUG FIXES

- Make sure that the assertion that checks for NULL or NA is operating
  in a
  scalar. For vectors use another assertion.

- Code clean-up. Remove unneeded files.

[TBSignatureProfiler](/packages/TBSignatureProfiler)
-------------------

                        Changes in version 1.5.0                        

- 68 signatures currently available

Bug Fixes

- Fixed incorrect names of signatures in Tbcommon and
common_sigAnnotData objects

Major Changes

- Added Zimmer_RES_3, Gong_OD_4, Bloom_RES_268, and Bloom_RES_558
- Added Sivakumaran_11 and Mendelsoh_RISK_11 signatures
- Added Estevez_133, Estevez_259, LauxdaCosta_OD_3, and Maertzdorf_15
signatures to the package
- Added Chen_HIV_4, Gliddon_HIV_3, Gliddon_2_OD_4, Kulkarni_HIV_2,
and
Heycken_FAIL_22 signatures
- Added a COVIDsignatures object to the package that can be used to
profile COVID-19 gene transcript signatures, thanks to collaborator
Dylan Sheerin (WEHI)
- Added functions to evaluate some signatures using their original
models from Johnson lab member Xutao Wang

Minor Changes

- Added mention of COVIDsignatures object to main vignette and
website
- Included the OG models tutorial on website (Xutao Wang)
- Added addTBsignature() to more easily facilitate updating
signatures
in package
- Added pROC option to obtain confidence intervals on AUC values as
part of \code{tableAUC()}
- Added citation for newly published paper

[TCGAbiolinks](/packages/TCGAbiolinks)
------------

                       Changes in version 2.21.1                        

- Function GDCPrepare for TARGET-ALL-P3 fixed

- Function getMC3MAF fixed

[TCGAutils](/packages/TCGAutils)
---------

                       Changes in version 1.14.0                        

Minor changes and bug fixes

- UUIDtoBarcode with the from_type = "file_id" argument now returns
the IDs in the proper order when more than one UUID is input.
- Update makeGRangesListFromCopyNumber examples with new names from
API e.g., 'associated_entities.entity_submitter_id'

[TOAST](/packages/TOAST)
-----

                        Changes in version 1.6.1                        

- Correct the citation error.

[topdownr](/packages/topdownr)
--------

                        Changes in version 1.15                         

- New version for Bioc 3.15 (devel)

Changes in version 1.15.1

- Add rmarkdown to Suggests:; see
https://github.com/yihui/knitr/issues/1864 for details [2021-07-27].

[ToxicoGx](/packages/ToxicoGx)
--------

                        Changes in version 1.3.4                        

- Fix bug in tests when run on Windows due to uninherited namespace
imports for testthat::context and testthat::expect_equal inside a
bplapply call

                        Changes in version 1.3.3                        

- Debugging BioC build ERROR caused by updates to CoreGx

                        Changes in version 1.3.2                        

- Fix a bug in computeLimmaDiffExpr where subsetting a ToxicoSet
doesn't subset the protocolData of a SummarizedExperiment, causing
coercing to an ExpressionSet inside the function to fail
- For now just deleting protocolData from the metadata of the
SummarizedExperiment, but will eventually need to be fixed upstream
in ORCESTRA

                        Changes in version 1.3.1                        

- Molecular profile data is now subset in test to keep package size
down

                        Changes in version 1.3.0                        

- Spring 2021 Bioconductor release!

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.29.8                        

- Fix the a typo in read hic data.

                       Changes in version 1.29.7                        

- Fix the bug 'breaks' are not unique

                       Changes in version 1.29.6                        

- Add smooth curve to the tracks.

                       Changes in version 1.29.5                        

- Improve gene track plots.

                       Changes in version 1.29.4                        

- Fix the issue for auto-rescale lolliplot by emphasizing exon region
  when there are continues exons.

                       Changes in version 1.29.3                        

- Add the possibility to lolliplot to emphasize exon or intron region.

                       Changes in version 1.29.2                        

- Fix the yaxis when user supplied yaxis is greater than max scores.

                       Changes in version 1.29.1                        

- Update documentation lolliplot for rescale parameter.

[transformGamPoi](/packages/transformGamPoi)
---------------

                         Changes in version 0.1                         

- Add clipping functionality to residual_transform()
- Add check against residual_type argument in transformGamPoi()
- Fix bug in acosh_transform() related to sparse input and on_disk =
FALSE
- Change default of overdispersion_shrinkage to TRUE if
overdispersion
= TRUE for acosh_transform() and shifted_log_transform()

                        Changes in version 0.1.0                        

- Initial release of transformGamPoi on GitHub
https://github.com/const-ae/transformGamPoi

[transomics2cytoscape](/packages/transomics2cytoscape)
--------------------

                        Changes in version 1.2.2                        

BUG FIX

- Improved tab delimiter processing, argument processing for
  RCy3::getEdgeInfo()

                        Changes in version 1.2.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- There was a bug in version 1.2.0 with the release update of
  Bioconductor.
  That's fixed in 1.2.1.

[TreeAndLeaf](/packages/TreeAndLeaf)
-----------

                        Changes in version 1.5.1                        

- Object transformation functions condensed into the
  `treeAndLeaf` function, to be made automatically [2020-08-24].

- TreeAndLeaf function became more automated.

- The igraph object manipulation was upgraded, through the use of
  RedeR functions.

[treeio](/packages/treeio)
------

                       Changes in version 1.17.2                        

- allow additional parameter to pass to drop.tip methods (2021-06-23,
Wed, @xiangpin, #62)
- as.phylo and as.treedata for data.frame (2021-06-12, Sat)
- as.ultrametric method to force a tree to be ultrametric
(2021-06-09,
Wed)
- introduce force.ultrametric parameter in read.mcmctree

                       Changes in version 1.17.1                        

- read.mcmctree for PAML MCMCTree result (2021-06-04, Fri)

[TRESS](/packages/TRESS)
-----

                 Changes in version 0.1.0 (2021-07-03)                  

- Submitted to Bioconductor

[tripr](/packages/tripr)
-----

                 Changes in version 0.99.0 (2021-08-15)                 

- Submitted to Bioconductor
- Added a NEWS.md file to track changes to the package.

[TVTB](/packages/TVTB)
----

                 Changes in version 1.19.1 (2021-08-30)                 

Bug fix

- Substitute || by |.

[txcutr](/packages/txcutr)
------

                       Changes in version 0.99.0                        

NEW FEATURES

- Staged for Bioconductor submission.

SIGNIFICANT USER-VISIBLE CHANGES

- None.

BUG FIXES

- None.

                        Changes in version 0.3.2                        

NEW FEATURES

- Improved vignette

SIGNIFICANT USER-VISIBLE CHANGES

- None.

BUG FIXES

- None.

                        Changes in version 0.3.1                        

NEW FEATURES

- Compressed outputs.
- Tests for proper handling of transitive merging. Overlaps that
merge
A -> B and B -> C, but not A -> C, will output A -> C and B -> C.
That is, transitivity is applied and the final output will always
use the distal most transcript in a chain as the final output.

SIGNIFICANT USER-VISIBLE CHANGES

- All export*() methods now include automatic detection of .gz
filenames, which toggles the use of compressed (gzip) exports.

BUG FIXES

- None.

                        Changes in version 0.3.0                        

NEW FEATURES

- Merge table generation and exporting.

SIGNIFICANT USER-VISIBLE CHANGES

- Adds generateMergeTable() and exportMergeTable() for creating a
merge table for transcripts that are not separated by a thresholded
distance. Such files can be used by transcript quantification tools
to specify what transcripts should be merged.

BUG FIXES

- None.

                        Changes in version 0.2.2                        

NEW FEATURES

- None.

SIGNIFICANT USER-VISIBLE CHANGES

- None.

BUG FIXES

- The BPPARAM was not being passed through to internal bplapply
calls.

                        Changes in version 0.2.1                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- Provide more control over parallel execution.

SIGNIFICANT USER-VISIBLE CHANGES

- The truncateTxome() method includes an optional BPPARAM with which
users can pass a specific BiocParallelParam. If not provided, it
will respect the result of BiocParallel::bpparam(), which can be
globally set using BiocParallel::register().

BUG FIXES

- None.

                        Changes in version 0.2.0                        

NEW FEATURES

- Adds deduplication behavior. Note that deduplication does not
exclude transcripts from different genes.

SIGNIFICANT USER-VISIBLE CHANGES

- The truncateTxome() method now deduplicates transcripts spanning
identical ranges after truncation.

BUG FIXES

- None.

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.12.0                        

NEW FEATURES

- New function, sequence_complexity(): Using either the
  Wootton-Federhen,
  Trifonov, or DUST algorithms, calculate sequence complexity in
  sliding
  windows. A version for small arbitrary strings is also provided:
  calc_complexity().

- New function, mask_ranges(): Similarly to mask_seqs(), mask specific
  positions in a XStringSet object by replacing the letters with a
  specific
  filler character.

- New function, motif_range(): Get the min/max range of possible
  logodds
  scores for a motif.

- New function, calc_windows(): Utility function for calculating
  coordinates for sliding windows.

- New function, window_string(): Utility function for retrieving
  sliding
  windows in a string.

- New function, slide_fun(): Utility function which wraps
  window_string()
  and vapply() together.

- motif_pvalue(method): P-values and scores can now be calculated
  dynamically instead of exhaustively, substantially increasing both
  speed
  and accuracy for bigger jobs. The previous exhaustive method can
  still be
  used however, as the dynamic method does not allow non-finite values
  and
  thus must be pseudocount-adjusted.

- scan_sequences(calc.pvals, calc.qvals, motif_pvalue.method,
  calc.qvals.method): The calc.pvals argument defaults to TRUE.
  The P-value calculation method now defaults to dynamic P-values (the
  previous method was an exhaustive calculation), though this can be
  changed via motif_pvalue.method. Additionally, adjusted P-values can
  be
  calculated as either BH, FDR or a Bonferroni-adjusted P-value. More
  details can be found in the Sequence Searches vignette.

- write_homer(threshold, threshold.type): Finer control over the final
  motif logodds threshold included with the written motif is now
  available, using the style of argument parsing from scan_sequences().
  The previous logodds_threshold argument is now deprecated and set to
  NULL, but if set (e.g. an older script is being re-run) then the old
  behaviour of write_homer() will be used.

MINOR CHANGES

- New global option, options(pseudocount.warning): Disable the message
  printed when a motif is pseudocount-adjusted.

- Slight performance gains in get_bkg() window code.

- motif_pvalue(): Clarify that, indeed, background probabilities are
  taken
  into account when calculating P-values from score inputs. The
  background
  adjustment takes place during the initial conversion to PWM.

- motif_pvalue(): When bkg.probs are provided, use those when
  converting to
  a PWM.

- scan_sequences(): The default threshold is now 0.0001 (using
  threshold.type = "pvalue").

- The axis text in view_motifs() is now black instead of grey.

- create_motif(): When a named background vector is provided, it is
  sorted according to the alphabet characters.

- scan_sequences(): Check that the sequences aren't shorter than the
  motifs.

- print.universalmotif_df: Changed warning message when subsetting to
  an
  incomplete universalmotif_df object. Also added a way to turn off
  informative messages/warnings via the boolean
  universalmotif_df.warning
  global option.

- Miscellaneous changes and additions to the vignettes and various
  function manual pages.

                       Changes in version 1.10.2                        

BUG FIXES

- read_homer() now correct parses enrichment P-value and logodds score.

                       Changes in version 1.10.1                        

BUG FIXES

- Restore temporarily disabled ggtree(layout="daylight") example in the
  MotifComparisonAndPvalues.Rmd vignette, as tidytree is now patched.

- Fixed some awkwardness in view_motifs() panel spacing and title
  justification.

[VAExprs](/packages/VAExprs)
-------

                 Changes in version 0.99.0 (2021-06-07)                 

- submission to Bioconductor

[velociraptor](/packages/velociraptor)
------------

                        Changes in version 1.3.1                        

- Add typing_extensions to environment

[veloviz](/packages/veloviz)
-------

                     Changes in version 0.0.0.9000                      

- Added a NEWS.md file to track changes to the package.

[ViSEAGO](/packages/ViSEAGO)
-------

                         Changes in version 1.7                         

- Ensembl2GO() biomart update

[vissE](/packages/vissE)
-----

                        Changes in version 1.2.0                        

- Added weighting to text-mining analysis. Word clouds can now
incorporate statistics.
- Network plotting now performed using ggraph.
- Removal of excessive warnings produced when performing text-mining.
- Added PPI exploration of clusters

[weitrix](/packages/weitrix)
-------

                        Changes in version 1.5.1                        

- weitrix_sd_confects now has an option to drop the assumption of
normally distributed weighted residuals.

[wppi](/packages/wppi)
----

                 Changes in version 1.1.3 (2021-10-05)                  

- The workflow calculates Protein-Protein Interaction weights and
scores genes
- Database knowledge is automatically fetched from OmniPath, Gene
Ontology and Human Phenotype Ontology
- Submitted to Bioconductor

[xcms](/packages/xcms)
----

                       Changes in version 3.15.5                        

- Disable testing on windows i386, providing some speedup

- Disable parallel processing on Windows, causing an issue in testthat
  on BioC build check

                       Changes in version 3.15.4                        

- Fix in `plot` with `type = "XIC"` to plot an empty plot if no data is
  present.

- Skip re-indexing of peaks to features if not necessary. This results
  in
  performance improvements for MS1 only data.

                       Changes in version 3.15.3                        

- Add `manualFeatures` allowing to manually define and add features to
  an
  `XCMSnExp` object.

- Add `plotChromatogramsOverlay` function to support plotting of
  multiple EICs
  from the same sample into the same plot (eventually stacked).

- Add feature grouping by EIC similarity: `EicSimilarityParam`.

- Import `compareChromatograms` from `MSnbase`.

- Add feature grouping by similar retention time: `SimilarRtimeParams.

- Add feature grouping by similarity of feature abundances across
  samples:
  `AbundanceSimilarityParam`.

- Add feature grouping methodology based on `MsFeatures`.

                       Changes in version 3.15.2                        

- Fix LC-MS/MS vignette.

                       Changes in version 3.15.1                        

- Compatibility fix for nls() in R >= 4.1, contributed by Rick Helmus.

[zellkonverter](/packages/zellkonverter)
-------------

                        Changes in version 1.4.0                        

- Add arguments to control how slots are converted in
  AnnData2SCE() and SCE2AnnData(). Each slot can now be fully
  converted, skipped entirely or only selected items converted.

- Add support for converting the raw slot to an altExp in
  AnnData2SCE()

- Add recursive conversion of lists in AnnData2SCE()

- Add progress messages to various functions. These can be
  controlled by function arguments or a global variable.

- Add long tests for various public datasets. This should help to
  make the package more robust

- Fix bug in converting dgRMatrix sparse matrices

- Correctly handle DataFrame objects stored in adata.obsm

NEWS from new and existing Data Experiment Packages
===================================


[benchmarkfdrData2019](/packages/benchmarkfdrData2019)
--------------------

                 Changes in version 1.7.1 (2021-06-18)                  

- Added `rmarkdown` to Suggests in DESCRIPTION to resolve changes in
  `knitr`

[BioImageDbs](/packages/BioImageDbs)
-----------

                        Changes in version 1.0.2                        

NEW FEATURES

- Updated the vignettes.

                        Changes in version 1.0.1                        

NEW FEATURES

- Modified the file format from Rda to rds.

[bodymapRat](/packages/bodymapRat)
----------

                 Changes in version 1.9.1 (2021-06-18)                  

- Added `rmarkdown` to Suggests in DESCRIPTION to resolve changes in
  `knitr`

[curatedMetagenomicData](/packages/curatedMetagenomicData)
----------------------

                        Changes in version 3.2.0                        

- The curatedMetagenomicData() function now has a rownames argument:

- "long", the default character string derived from MetaPhlAn3

- "short", the NCBI Taxonomy species name from the CHOCOPhlAn
  database

- "short" row names are validated against NCBI Taxonomy with
  taxize

- "NCBI", the NCBI Taxonomy ID from the CHOCOPhlAn database

- "NCBI" row names are validated against NCBI Taxonomy with
  taxize

- rowData becomes NCBI Taxonomy ID numbers instead of taxa
  names

- The sparse matrix data structure was switched from dgTMatrix to
  dgCMatrix

- A few studies were reprocessed because of a minor error related to
  MetaPhlAn3

- Changes inside the package were made to address bugs discovered by
  users

- The combined_metadata object has been removed

[curatedTBData](/packages/curatedTBData)
-------------

                 Changes in version 0.99.4 (2021-10-18)                 

- Update vignette output

- Edit DESCRIPTION file

                 Changes in version 0.99.3 (2021-10-07)                 

- Remove analysis/graphical functions from the master branch

- Accepted by Bioconductor

                 Changes in version 0.99.2 (2021-09-24)                 

- Remove cached files from git history

- Modify packages based on reviewer's comment

                 Changes in version 0.99.0 (2021-09-16)                 

- The curatedTBData package collects 49 transcriptomic studies

- Package vignette is updated to Rmd syntax and uses BiocStyle

- All data is reprocessed in R (v4.1)

- Move all data to ExperimentHub

- Added a NEWS.md file to track changes to the package

- Submitted to Bioconductor

[depmap](/packages/depmap)
------

                Changes in version 1.7.1                  

- 21Q3 data added for crispr, copyNumber, TPM, mutationCalls and
        metadata datasets. Newer versions for the other datasets were
        not released.
        
- CERES CRISPR data has been deprecated and has been replaced with
        Chronos CRISPR dependency in 21Q3 and all future releases. For
        more information, see:
        https://cancerdatascience.org/blog/posts/ceres-chronos/


[dorothea](/packages/dorothea)
--------

                 Changes in version 1.5.2 (2021-10-13)                  

- Daniel Dimitrov is assigned as the new maintainer

                 Changes in version 1.4.2 (2021-10-08)                  

- Fixed lazy data warning

- Improved test coverage

                 Changes in version 1.4.1 (2021-05-25)                  

- Rebuild all regulons

- Fixed ambiguously mode of regulation in mouse regulons

[easierData](/packages/easierData)
----------

                       Changes in version 0.99.0                        

- All set for the Bioconductor submission!

                        Changes in version 0.9.0                        

- Getting ready for the submission to Bioconductor

- Added a NEWS.md file to track changes to the package.

[imcdatasets](/packages/imcdatasets)
-----------

                 Changes in version 1.1.1 (2021-05-25)                  

- Added documentation (dataset list in README, print out examples in
  help pages)

[LRcellTypeMarkers](/packages/LRcellTypeMarkers)
-----------------

                        Changes in version 3.13                         

- Add MSigDB datasets

- Add note in vignettes

[microbiomeDataSets](/packages/microbiomeDataSets)
------------------

                 Changes in version 1.1.5 (2021-09-08)                  

- GrieneisenTS data added

- HintikkaXO data added

- Minor fixes

[msigdb](/packages/msigdb)
------

                        Changes in version 1.2.0                        

- added MSigDB v7.4

- removed gene-sets from the "archived" category from all collections

- removed direct object referencing functions (e.g.
  msigdb.v7.2.hs.SYM()). Objects should be retrieved using the
  getMsigdb() function only or by querying the ExperimentHub

- added IMEx PPI data

[mtbls2](/packages/mtbls2)
------

                       Changes in version 1.23.4                        

- temporarily disable all vignettes until CAMERA issue is fixed

                       Changes in version 1.23.2                        

- switch to mzML files in the ISA-Tab metadata

- temporarily disable the vignette for the old xcms interface causing a
  build failure

                       Changes in version 1.23.1                        

- add mzML versions of mzData files converted by OpenMS FileConverter

[NxtIRFdata](/packages/NxtIRFdata)
----------

                 Changes in version 0.99.3 (2021-09-27)                 

- Now uses BiocFileCache to store a copy of ExperimentHub resources.
  Allows for faster subsequent recalls

                 Changes in version 0.99.0 (2021-09-17)                 

- Submitted to Bioconductor

[ptairData](/packages/ptairData)
---------

                 Changes in version 1.1.1 (2021-03-05)                  

- simulation files added

[RforProteomics](/packages/RforProteomics)
--------------

                       Changes in version 1.31.1                        

- Fix mztab file name (required for new rpx)

                       Changes in version 1.31.0                        

- New version for Bioc 3.14 (devel)

[RLHub](/packages/RLHub)
-----

                       Changes in version 0.99.6                        

- RLHub provides convenient access to the processed datasets
        within RLBase. A full step-by-step protocol for replicating RLHub is
        found here https://github.com/Bishop-Laboratory/RLBase-data


[RnBeads.mm10](/packages/RnBeads.mm10)
------------

                        Changes in version 2.1.1                        

- Temporarily rolled Ensemble version back to 75 to eliminate numerous
  non-coding transcripts

- Methylation array information based on the latest manifest file

                        Changes in version 2.1.0                        

- Release May 2021

- Added annotation for the Mouse Methylation Bead Chip (thanks to Maxi
  Schoenung for his great contribution!).

- Updated SNP information to GRCm38.p4, the last version available
  through NCBI.

[scATAC.Explorer](/packages/scATAC.Explorer)
---------------

                       Changes in version 0.99.3                        

- removed usage of paste() function in error handling functions

                       Changes in version 0.99.2                        

- Added NEWS.md file to track changes to the package.

- Formatted functions to shorten lines

- Removed usage of T/F in test cases

- Removed usage of 'paste' in condition signals

                       Changes in version 0.99.1                        

- added tests for queryATAC function

- added new error checking to queryATAC and fetchATAC functions

[scpdata](/packages/scpdata)
-------

                        Changes in version 1.1.1                        

- Added the schoof2021 dataset <2021-10-06>


[SingleCellMultiModal](/packages/SingleCellMultiModal)
--------------------

                        Changes in version 1.6.0                        

New features

- scMultiome version 1.0.1 provides the 10X format for RNAseq data.

Bug fixes and minor improvements

- Updates to seqFISH vignette and documentation.

- Updated to changes in SummarizedExperiment where assayDimnames are
  checked.

- scNMT defaults to version '1.0.0's QC filtered cells. For unfiltered
  cells see version section in ?scNMT.

[spatialDmelxsim](/packages/spatialDmelxsim)
---------------

                       Changes in version 0.99.9                        

- Submitting to Bioconductor as a ExperimentHub data package.

[TabulaMurisSenisData](/packages/TabulaMurisSenisData)
--------------------

                       Changes in version 0.99.1                        

- Add infoOnly argument to get details about download size

                       Changes in version 0.99.0                        

- Add documentation

- Prepare for Bioconductor submission

                        Changes in version 0.1.0                        

- Add a NEWS.md file to track changes to the package.

[TMExplorer](/packages/TMExplorer)
----------

                        Changes in version 1.2.1                        

- Added a summary column to the metadata table and SingleCellExperiment
  metadata

[tuberculosis](/packages/tuberculosis)
------------

                        Changes in version 1.0.0                        

- tuberculosis is now available in Bioconductor Release 3.14

NEWS from new and existing Workflows
===================================


[GeoMxWorkflows](/packages/GeoMxWorkflows)
--------------

                       Changes in version 0.99.4                        

- Updated workflow image

                       Changes in version 0.99.3                        

- Removed empty DCC file from example dataset
- Adjusted QC plotting histogram function to remove user-defined
limits

                       Changes in version 0.99.2                        

- Version update for bioconductor build review

                       Changes in version 0.99.1                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package
- Removed DESCRIPTION and Namespace conflicts
- Updated styling of vignette (GeomxTools_RNA-NGS_Analysis.Rmd) to fit with BioC

NEWS from new and existing books
===================================

No new NEWS to report

Deprecated and Defunct Packages
===============================


Forty seven software packages were removed from this release (after being deprecated
in Bioc 3.13): 
AffyExpress, affyQCReport, AnnotationFuncs, ArrayTools, bigmemoryExtras,
BiocCaseStudies, CancerMutationAnalysis, ChIPSeqSpike, CompGO, CoRegFlux,
CrossICC, cytofast, DBChIP, dexus, EasyqpcR, EDDA, eisa, ELBOW, ExpressionView,
FlowRepositoryR, genoset, HCABrowser, HCAExplorer, HCAMatrixBrowser, Imetagene,
mdgsa, metagenomeFeatures, methyAnalysis, MSEADbi, OutlierD,
pcot2, PCpheno, Polyfit, POST, RchyOptimyx, RDAVIDWebService, RNAither,
RNAprobR, rnaSeqMap, SAGx, samExploreR, seqplots, simulatorZ, SSPA, ToPASeq,
XBSeq, yaqcaffy

Please note:  CexoR and IntramiRExploreR, previously announced as deprecated in
3.13, fixed their packages and remained in Bioconductor.

Twenty three software are deprecated in this release and will be removed in Bioc 3.15:
affyPara, ALPS, alsace, BrainStars, destiny, dualKS, ENCODExplorer,
ENVISIONQuery, FindMyFriends, GeneAnswers, gramm4R, KEGGprofile, MouseFM,
MSGFgui, MSGFplus, MSstatsTMTPTM, PanVizGenerator, predictionet, RGalaxy,
scClassifR, slinky, SRGnet, SwimR

Eleven experimental data packages were removed this release (after being
deprecated in BioC 3.13):
ceu1kg, ceu1kgv, ceuhm3, cgdv17, dsQTL, facsDorit, gskb, hmyriB36, JctSeqData,
MAQCsubsetAFX, yri1kgv

Five experimental data packages are deprecated in this release and will be
removed in Bioc 3.15:
ABAData, brainImageRdata, PCHiCdata, RITANdata, tcgaWGBSData.hg19

Eighty seven annotation packages were removed from this release (after being deprecated
in Bioc 3.13): 
12 LRBase.XXX.eg.db packages (replaced with AHLRBaseDbs), 
MafDb.gnomAD.r3.0.GRCh38, MafH5.gnomAD.r3.0.GRCh38, 73 MeSH.XXX.eg.db packages 
(replaced with AHMeSHDbs)

One annotation package was deprecated in this release and will be removed in
Bioc 3.15:
org.Pf.plasmo.db

One workflow package was removed from this release (after being deprecated
in Bioc 3.13):
eQTL

No workflow packages were deprecated in this release.