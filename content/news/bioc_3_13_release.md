May 20, 2021

Bioconductors:

We are pleased to announce Bioconductor 3.13, consisting of 
1993 software packages, 395 experiment data packages, xxx annotation
packages, and 28 workflows.

There are 133 new software packages, 22 new data experiment packages,
7 new annotation packages, 1 new workflow, no new books, and many updates and
improvements to existing packages; Bioconductor 3.13 is compatible with R 4.1.0,
and is supported on Linux, 32- and 64-bit Windows, and macOS 10.14.6 Mojave
or higher.  This release will include an updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3]
for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.13](#getting-started-with-bioconductor-313)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from new and existing software packages](#news-from-new-and-existing-software-packages)
* [NEWS from new and existing data experiment packages](#news-from-new-and-existing-data-experiment-packages)
* [NEWS from new and existing workflows](#news-from-new-and-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.13
======================================

To update to or install Bioconductor 3.13:

1. Install R 4.1.0. Bioconductor 3.13 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).

New Software Packages
=====================

There are 133 new software packages in this release of Bioconductor.


- [airpart](/packages/airpart) Airpart identifies sets of genes
  displaying differential cell-type-specific allelic imbalance across
  cell types or states, utilizing single-cell allelic counts. It
  makes use of a generalized fused lasso with binomial observations
  of allelic counts to partition cell types by their allelic
  imbalance. Alternatively, a nonparametric method for partitioning
  cell types is offered. The package includes a number of
  visualizations and quality control functions for examining single
  cell allelic imbalance datasets.

- [autonomics](/packages/autonomics) This package offers a generic
  and intuitive solution for cross-platform omics data analysis. It
  has functions for import, preprocessing, exploration, contrast
  analysis and visualization of omics data. It follows a tidy,
  functional programming paradigm.

- [awst](/packages/awst) We propose an Asymmetric Within-Sample
  Transformation (AWST) to regularize RNA-seq read counts and reduce
  the effect of noise on the classification of samples. AWST
  comprises two main steps: standardization and smoothing. These
  steps transform gene expression data to reduce the noise of the
  lowly expressed features, which suffer from background effects and
  low signal-to-noise ratio, and the influence of the highly
  expressed features, which may be the result of amplification bias
  and other experimental artifacts.

- [barcodetrackR](/packages/barcodetrackR) barcodetrackR is an R
  package developed for the analysis and visualization of clonal
  tracking data. Data required is samples and tag abundances in
  matrix form. Usually from cellular barcoding experiments,
  integration site retrieval analyses, or similar technologies.

- [biodb](/packages/biodb) The biodb package provides access to
  standard remote chemical and biological databases (ChEBI, KEGG,
  HMDB, ...), as well as to in-house local database files (CSV,
  SQLite), with easy retrieval of entries, access to web services,
  search of compounds by mass and/or name, and mass spectra matching
  for LCMS and MSMS. Its architecture as a development framework
  facilitates the development of new database connectors for local
  projects or inside separate published packages.

- [BioNERO](/packages/BioNERO) BioNERO aims to integrate all aspects
  of biological network inference in a single package, including data
  preprocessing, exploratory analyses, network inference, and
  analyses for biological interpretations. BioNERO can be used to
  infer gene coexpression networks (GCNs) and gene regulatory
  networks (GRNs) from gene expression data. Additionally, it can be
  used to explore topological properties of protein-protein
  interaction (PPI) networks. GCN inference relies on the popular
  WGCNA algorithm. GRN inference is based on the "wisdom of the
  crowds" principle, which consists in inferring GRNs with multiple
  algorithms (here, CLR, GENIE3 and ARACNE) and calculating the
  average rank for each interaction pair. As all steps of network
  analyses are included in this package, BioNERO makes users avoid
  having to learn the syntaxes of several packages and how to
  communicate between them. Finally, users can also identify
  consensus modules across independent expression sets and calculate
  intra and interspecies module preservation statistics between
  different networks.

- [BloodGen3Module](/packages/BloodGen3Module) The BloodGen3Module
  package provides functions for R user performing module repertoire
  analyses and generating fingerprint representations. Functions can
  perform group comparison or individual sample analysis and
  visualization by fingerprint grid plot or fingerprint heatmap.
  Module repertoire analyses typically involve determining the
  percentage of the constitutive genes for each module that are
  significantly increased or decreased. As we describe in
  details;https://www.biorxiv.org/content/10.1101/525709v2 and
  https://pubmed.ncbi.nlm.nih.gov/33624743/, the results of module
  repertoire analyses can be represented in a fingerprint format,
  where red and blue spots indicate increases or decreases in module
  activity. These spots are subsequently represented either on a
  grid, with each position being assigned to a given module, or in a
  heatmap where the samples are arranged in columns and the modules
  in rows.

- [bnem](/packages/bnem) bnem combines the use of indirect
  measurements of Nested Effects Models (package mnem) with the
  Boolean networks of CellNOptR. Perturbation experiments of
  signalling nodes in cells are analysed for their effect on the
  global gene expression profile. Those profiles give evidence for
  the Boolean regulation of down-stream nodes in the network, e.g.,
  whether two parents activate their child independently (OR-gate) or
  jointly (AND-gate).

- [BumpyMatrix](/packages/BumpyMatrix) Implements the BumpyMatrix
  class and several subclasses for holding non-scalar objects in each
  entry of the matrix. This is akin to a ragged array but the
  raggedness is in the third dimension, much like a bumpy surface -
  hence the name. Of particular interest is the BumpyDataFrameMatrix,
  where each entry is a Bioconductor data frame. This allows us to
  naturally represent multivariate data in a format that is
  compatible with two-dimensional containers like the
  SummarizedExperiment and MultiAssayExperiment objects.

- [CAEN](/packages/CAEN) With the development of high-throughput
  techniques, more and more gene expression analysis tend to replace
  hybridization-based microarrays with the revolutionary
  technology.The novel method encodes the category again by employing
  the rank of samples for each gene in each class. We then consider
  the correlation coefficient of gene and class with rank of sample
  and new rank of category. The highest correlation coefficient genes
  are considered as the feature genes which are most effective to
  classify the samples.

- [cbpManager](/packages/cbpManager) This R package provides an R
  Shiny application that enables the user to generate, manage, and
  edit data and metadata files suitable for the import in cBioPortal
  for Cancer Genomics. Create cancer studies and edit its metadata.
  Upload mutation data of a patient that will be concatenated to the
  data_mutation_extended.txt file of the study. Create and edit
  clinical patient data, sample data, and timeline data. Create
  custom timeline tracks for patients.

- [CelliD](/packages/CelliD) CelliD is a clustering-free multivariate
  statistical method for the robust extraction of per-cell gene
  signatures from single-cell RNA-seq. CelliD allows unbiased cell
  identity recognition across different donors, tissues-of-origin,
  model organisms and single-cell omics protocols. The package can
  also be used to explore functional pathways enrichment in single
  cell data.

- [cellmigRation](/packages/cellmigRation) Import TIFF images of
  fluorescently labeled cells, and track cell movements over time.
  Parallelization is supported for image processing and for fast
  computation of cell trajectories. In-depth analysis of cell
  trajectories is enabled by 15 trajectory analysis functions.

- [censcyt](/packages/censcyt) Methods for differential abundance
  analysis in high-dimensional cytometry data when a covariate is
  subject to right censoring (e.g. survival time) based on multiple
  imputation and generalized linear mixed models.

- [CIMICE](/packages/CIMICE) CIMICE is a tool in the field of tumor
  phylogenetics and its goal is to build a Markov Chain (called
  Cancer Progression Markov Chain, CPMC) in order to model tumor
  subtypes evolution. The input of CIMICE is a Mutational Matrix, so
  a boolean matrix representing altered genes in a collection of
  samples. These samples are assumed to be obtained with single-cell
  DNA analysis techniques and the tool is specifically written to use
  the peculiarities of this data for the CMPC construction.

- [CNVgears](/packages/CNVgears) This package contains a set of
  functions to perform several type of processing and analysis on
  CNVs calling pipelines/algorithms results in an integrated manner
  and regardless of the raw data type (SNPs array or NGS). It
  provides functions to combine multiple CNV calling results into a
  single object, filter them, compute CNVRs (CNV Regions) and
  inheritance patterns, detect genic load, and more. The package is
  best suited for studies in human family-based cohorts.

- [CNViz](/packages/CNViz) CNViz takes probe, gene, and segment-level
  log2 copy number ratios and launches a Shiny app to visualize your
  sample's copy number profile. You can also integrate loss of
  heterozygosity (LOH) and single nucleotide variant (SNV) data.

- [ComPrAn](/packages/ComPrAn) This package is for analysis of SILAC
  labeled complexome profiling data. It uses peptide table in
  tab-delimited format as an input and produces ready-to-use tables
  and plots.

- [conclus](/packages/conclus) CONCLUS is a tool for robust
  clustering and positive marker features selection of single-cell
  RNA-seq (sc-RNA-seq) datasets. It takes advantage of a consensus
  clustering approach that greatly simplify sc-RNA-seq data analysis
  for the user. Of note, CONCLUS does not cover the preprocessing
  steps of sequencing files obtained following next-generation
  sequencing. CONCLUS is organized into the following steps:
  Generation of multiple t-SNE plots with a range of parameters
  including different selection of genes extracted from PCA. Use the
  Density-based spatial clustering of applications with noise
  (DBSCAN) algorithm for idenfication of clusters in each generated
  t-SNE plot. All DBSCAN results are combined into a cell similarity
  matrix. The cell similarity matrix is used to define "CONSENSUS"
  clusters conserved accross the previously defined clustering
  solutions. Identify marker genes for each concensus cluster.

- [condiments](/packages/condiments) This package encapsulate many
  functions to conduct a differential topology analysis. It focuses
  on analyzing an 'omic dataset with multiple conditions. While the
  package is mostly geared toward scRNASeq, it does not place any
  restriction on the actual input format.

- [CONSTANd](/packages/CONSTANd) Normalizes a data matrix `data` by
  raking (using the RAS method by Bacharach, see references) the
  Nrows by Ncols matrix such that the row means and column means
  equal 1. The result is a normalized data matrix `K=RAS`, a product
  of row mulipliers `R` and column multipliers `S` with the original
  matrix `A`. Missing information needs to be presented as `NA`
  values and not as zero values, because CONSTANd is able to ignore
  missing values when calculating the mean. Using CONSTANd
  normalization allows for the direct comparison of values between
  samples within the same and even across different
  CONSTANd-normalized data matrices.

- [cosmosR](/packages/cosmosR) COSMOS (Causal Oriented Search of
  Multi-Omic Space) is a method that integrates phosphoproteomics,
  transcriptomics, and metabolomics data sets based on prior
  knowledge of signaling, metabolic, and gene regulatory networks. It
  estimated the activities of transcrption factors and kinases and
  finds a network-level causal reasoning. Thereby, COSMOS provides
  mechanistic hypotheses for experimental observations across
  mulit-omics datasets.

- [CTDquerier](/packages/CTDquerier) Package to retrieve and
  visualize data from the Comparative Toxicogenomics Database
  (http://ctdbase.org/). The downloaded data is formated as
  DataFrames for further downstream analyses.

- [cyanoFilter](/packages/cyanoFilter) An approach to filter out
  and/or identify phytoplankton cells from all particles measured via
  flow cytometry pigment and cell complexity information. It does
  this using a sequence of one-dimensional gates on pre-defined
  channels measuring certain pigmentation and complexity. The package
  is especially tuned for cyanobacteria, but will work fine for
  phytoplankton communities where there is at least one cell
  characteristic that differentiates every phytoplankton in the
  community.

- [CytoGLMM](/packages/CytoGLMM) The CytoGLMM R package implements
  two multiple regression strategies: A bootstrapped generalized
  linear model (GLM) and a generalized linear mixed model (GLMM).
  Most current data analysis tools compare expressions across many
  computationally discovered cell types. CytoGLMM focuses on just one
  cell type. Our narrower field of application allows us to define a
  more specific statistical model with easier to control statistical
  guarantees. As a result, CytoGLMM finds differential proteins in
  flow and mass cytometry data while reducing biases arising from
  marker correlations and safeguarding against false discoveries
  induced by patient heterogeneity.

- [dce](/packages/dce) Compute differential causal effects (dce) on
  (biological) networks. Given observational samples from a control
  experiment and non-control (e.g., cancer) for two genes A and B, we
  can compute differential causal effects with a (generalized) linear
  regression. If the causal effect of gene A on gene B in the control
  samples is different from the causal effect in the non-control
  samples the dce will differ from zero. We regularize the dce
  computation by the inclusion of prior network information from
  pathway databases such as KEGG.

- [decoupleR](/packages/decoupleR) Transcriptome profiling followed
  by differential gene expression analysis often leads to lists of
  genes that are hard to analyze and interpret. Downstream analysis
  tools can be used to summarize deregulation events into a smaller
  set of biologically interpretable features.  In particular, methods
  that estimate the activity of transcription factors (TFs) from gene
  expression are commonly used. It has been shown that the
  transcriptional targets of a TF yield a much more robust estimation
  of the TF activity than observing the expression of the TF itself.
  Consequently, for the estimation of transcription factor
  activities, a network of transcriptional regulation is required in
  combination with a statistical algorithm that summarizes the
  expression of the target genes into a single activity score. Over
  the years, many different regulatory networks and statistical
  algorithms have been developed, mostly in a fixed combination of
  one network and one algorithm. To systematically evaluate both
  networks and algorithms, we developed decoupleR , an R package that
  allows users to apply efficiently any combination provided.

- [DeepPINCS](/packages/DeepPINCS) The identification of novel
  compound-protein interaction (CPI) is important in drug discovery.
  Revealing unknown compound-protein interactions is useful to design
  a new drug for a target protein by screening candidate compounds.
  The accurate CPI prediction assists in effective drug discovery
  process. To identify potential CPI effectively, prediction methods
  based on machine learning and deep learning have been developed.
  Data for sequences are provided as discrete symbolic data. In the
  data, compounds are represented as SMILES (simplified
  molecular-input line-entry system) strings and proteins are
  sequences in which the characters are amino acids. The outcome is
  defined as a variable that indicates how strong two molecules
  interact with each other or whether there is an interaction between
  them. In this package, a deep-learning based model that takes only
  sequence information of both compounds and proteins as input and
  the outcome as output is used to predict CPI. The model is
  implemented by using compound and protein encoders with useful
  features. The CPI model also supports other modeling tasks,
  including protein-protein interaction (PPI), chemical-chemical
  interaction (CCI), or single compounds and proteins. Although the
  model is designed for proteins, DNA and RNA can be used if they are
  represented as sequences.

- [DelayedRandomArray](/packages/DelayedRandomArray) Implements a
  DelayedArray of random values where the realization of the sampled
  values is delayed until they are needed. Reproducible sampling
  within any subarray is achieved by chunking where each chunk is
  initialized with a different random seed and stream. The usual
  distributions in the stats package are supported, along with
  scalar, vector and arrays for the parameters.

- [DExMA](/packages/DExMA) performing all the steps of gene
  expression meta-analysis without eliminating those genes that are
  presented in almost two data sets. It provides the necessary
  functions to be able to perform the different methods of gene
  expression meta-analysis. In addition, it contains functions to
  apply quality controls, download GEO data sets and show graphical
  representations of the results.

- [diffUTR](/packages/diffUTR) The diffUTR package provides a uniform
  interface and plotting functions for limma/edgeR/DEXSeq -powered
  differential bin/exon usage. It includes in addition an improved
  version of the limma::diffSplice method. Most importantly, diffUTR
  further extends the application of these frameworks to differential
  UTR usage analysis using poly-A site databases.

- [dir.expiry](/packages/dir.expiry) Implements an expiration system
  for access to versioned directories. Directories that have not been
  accessed by a registered function within a certain time frame are
  deleted. This aims to reduce disk usage by eliminating obsolete
  caches generated by old versions of packages.

- [drugTargetInteractions](/packages/drugTargetInteractions) Provides
  utilities for identifying drug-target interactions for sets of
  small molecule or gene/protein identifiers. The required
  drug-target interaction information is obained from a local SQLite
  instance of the ChEMBL database. ChEMBL has been chosen for this
  purpose, because it provides one of the most comprehensive and best
  annotatated knowledge resources for drug-target information
  available in the public domain.

- [epialleleR](/packages/epialleleR) Epialleles are specific DNA
  methylation patterns that are mitotically and/or meiotically
  inherited. This package calls hypermethylated epiallele frequencies
  at the level of genomic regions or individual cytosines in
  next-generation sequencing data using binary alignment map (BAM)
  files as an input. Other functionality includes computing the
  empirical cumulative distribution function for per-read beta
  values, and testing the significance of the association between
  epiallele methylation status and base frequencies at particular
  genomic positions (SNPs).

- [epidecodeR](/packages/epidecodeR) epidecodeR is a package capable
  of analysing impact of degree of DNA/RNA epigenetic chemical
  modifications on dysregulation of genes or proteins. This package
  integrates chemical modification data generated from a host of
  epigenomic or epitranscriptomic techniques such as ChIP-seq,
  ATAC-seq, m6A-seq, etc. and dysregulated gene lists in the form of
  differential gene expression, ribosome occupancy or differential
  protein translation and identify impact of dysregulation of genes
  caused due to varying degrees of chemical modifications associated
  with the genes. epidecodeR generates cumulative distribution
  function (CDF) plots showing shifts in trend of overall log2FC
  between genes divided into groups based on the degree of
  modification associated with the genes. The tool also tests for
  significance of difference in log2FC between groups of genes.

- [epigraHMM](/packages/epigraHMM) epigraHMM provides a set of tools
  for the analysis of epigenomic data based on hidden Markov Models.
  It contains two separate peak callers, one for consensus peaks from
  biological/technical replicates, and and one for differial peaks
  from multi-replicate multi-condition experiments. For the latter,
  window-specific posterior probabilities associated with read count
  enrichment for every possible combinatorial pattern are provided.

- [EWCE](/packages/EWCE) Used to determine which cell types are
  enriched within gene lists. The package provides tools for testing
  enrichments within simple gene lists (such as human disease
  associated genes) and those resulting from differential expression
  studies. The package does not depend upon any particular Single
  Cell Transcriptome dataset and user defined datasets can be loaded
  in and used in the analyses.

- [FEAST](/packages/FEAST) Cell clustering is one of the most
  important and commonly performed tasks in single-cell RNA
  sequencing (scRNA-seq) data analysis. An important step in cell
  clustering is to select a subset of genes (referred to as
  “features”), whose expression patterns will then be used for
  downstream clustering. A good set of features should include the
  ones that distinguish different cell types, and the quality of such
  set could have significant impact on the clustering accuracy. FEAST
  is an R library for selecting most representative features before
  performing the core of scRNA-seq clustering. It can be used as a
  plug-in for the etablished clustering algorithms such as SC3,
  TSCAN, SHARP, SIMLR, and Seurat. The core of FEAST algorithm
  includes three steps: 1. consensus clustering; 2. gene-level
  significance inference; 3. validation of an optimized feature set.

- [fedup](/packages/fedup) An R package that tests for enrichment and
  depletion of user-defined pathways using a Fisher's exact test. The
  method is designed for versatile pathway annotation formats (eg.
  gmt, txt, xlsx) to allow the user to run pathway analysis on custom
  annotations. This package is also integrated with Cytoscape to
  provide network-based pathway visualization that enhances the
  interpretability of the results.

- [fgga](/packages/fgga) Package that implements the FGGA algorithm.
  This package provides a hierarchical ensemble method based ob
  factor graphs for the consistent GO annotation of protein coding
  genes. FGGA embodies elements of predicate logic, communication
  theory, supervised learning and inference in graphical models.

- [flowGraph](/packages/flowGraph) Identifies maximal differential
  cell populations in flow cytometry data taking into account
  dependencies between cell populations; flowGraph calculates and
  plots SpecEnr abundance scores given cell population cell counts.

- [fobitools](/packages/fobitools) A set of tools for interacting
  with Food-Biomarker Ontology (FOBI). A collection of basic
  manipulation tools for biological significance analysis, graphs,
  and text mining strategies for annotating nutritional data.

- [GenomicDistributions](/packages/GenomicDistributions) If you have
  a set of genomic ranges, this package can help you with
  visualization and comparison. It produces several kinds of plots,
  for example: Chromosome distribution plots, which visualize how
  your regions are distributed over chromosomes; feature distance
  distribution plots, which visualizes how your regions are
  distributed relative to a feature of interest, like Transcription
  Start Sites (TSSs); genomic partition plots, which visualize how
  your regions overlap given genomic features such as promoters,
  introns, exons, or intergenic regions. It also makes it easy to
  compare one set of ranges to another.

- [GenomicSuperSignature](/packages/GenomicSuperSignature) This
  package contains the index, which is the Replicable and
  interpretable Axes of Variation (RAV) extracted from public RNA
  sequencing datasets by clustering and averaging top PCs. This
  index, named as RAVindex, is further annotated with MeSH terms and
  GSEA. Functions to connect PCs from new datasets to RAVs, extract
  interpretable annotations, and provide intuitive visualization, are
  implemented in this package. Overall, this package enables
  researchers to analyze new data in the context of existing
  databases with minimal computing resources.

- [GEOfastq](/packages/GEOfastq) GEOfastq is used to download fastq
  files from the European Nucleotide Archive (ENA) starting with an
  accession from the Gene Expression Omnibus (GEO). To do this,
  sample metadata is retrieved from GEO and the Sequence Read Archive
  (SRA). SRA run accessions are then used to construct FTP and aspera
  download links for fastq files generated by the ENA.

- [GeomxTools](/packages/GeomxTools) Tools for NanoString
  Technologies GeoMx Technology.  Package provides functions for
  reading in DCC and PKC files based on an ExpressionSet derived
  object.  Normalization and QC functions are also included.

- [geva](/packages/geva) Statistic methods to evaluate variations of
  differential expression (DE) between multiple biological
  conditions. It takes into account the fold-changes and p-values
  from previous differential expression (DE) results that use
  large-scale data (*e.g.*, microarray and RNA-seq) and evaluates
  which genes would react in response to the distinct experiments.
  This evaluation involves an unique pipeline of statistical methods,
  including weighted summarization, quantile detection, cluster
  analysis, and ANOVA tests, in order to classify a subset of
  relevant genes whose DE is similar or dependent to certain
  biological factors.

- [granulator](/packages/granulator) granulator is an R package for
  the cell type deconvolution of heterogeneous tissues based on bulk
  RNA-seq data or single cell RNA-seq expression profiles. The
  package provides a unified testing interface to rapidly run and
  benchmark multiple state-of-the-art deconvolution methods. Data for
  the deconvolution of peripheral blood mononuclear cells (PBMCs)
  into individual immune cell types is provided as well.

- [hca](/packages/hca) This package provides users with the ability
  to query the Human Cell Atlas data repository for single-cell
  experiment data. The `projects()`, `files()`, `samples()` and
  `bundles()` functions retrieve summary information on each of these
  indexes; corresponding `*_details()` are available for individual
  entries of each index. File-based resources can be downloaded using
  `files_download()`. Advanced use of the package allows the user to
  page through large result sets, and to flexibly query the
  'list-of-lists' structure representing query responses.

- [HGC](/packages/HGC) `HGC` (short for Hierarchical Graph-based
  Clustering) is a R package for conducting hierarchical clustering
  on large-scale single-cell RNA-seq (scRNA-seq) data. The key idea
  is to construct a dendrogram of cells on their shared nearest
  neighbor (SNN) graph. `HGC` provides functions for building cell
  graphs and for conducting hierarchical clustering on the graph.

- [HiCDCPlus](/packages/HiCDCPlus) Systematic 3D interaction calls
  and differential analysis for Hi-C and HiChIP. The HiC-DC+
  (Hi-C/HiChIP direct caller plus) package enables principled
  statistical analysis of Hi-C and HiChIP data sets – including
  calling significant interactions within a single experiment and
  performing differential analysis between conditions given replicate
  experiments – to facilitate global integrative studies. HiC-DC+
  estimates significant interactions in a Hi-C or HiChIP experiment
  directly from the raw contact matrix for each chromosome up to a
  specified genomic distance, binned by uniform genomic intervals or
  restriction enzyme fragments, by training a background model to
  account for random polymer ligation and systematic sources of read
  count variation.

- [HubPub](/packages/HubPub) HubPub provides users with functionality
  to help with the Bioconductor Hub structures. The package provides
  the ability to create a skeleton of a Hub style package that the
  user can then populate with the necessary information. There are
  also functions to help add resources to the Hub package metadata
  files as well as publish data to the Bioconductor S3 bucket.

- [immunotation](/packages/immunotation) MHC (major
  histocompatibility complex) molecules are cell surface complexes
  that present antigens to T cells.  The repertoire of antigens
  presented in a given genetic background largely depends on the
  sequence of the encoded MHC molecules, and thus, in humans, on the
  highly variable HLA (human leukocyte antigen) genes of the
  hyperpolymorphic HLA locus. More than 28,000 different HLA alleles
  have been reported, with significant differences in allele
  frequencies between human populations worldwide. Reproducible and
  consistent annotation of HLA alleles in large-scale bioinformatics
  workflows remains challenging, because the available reference
  databases and software tools often use different HLA naming
  schemes. The package immunotation provides tools for consistent
  annotation of HLA genes in typical immunoinformatics workflows such
  as for example the prediction of MHC-presented peptides in
  different human donors. Converter functions that provide mappings
  between different HLA naming schemes are based on the MHC
  restriction ontology (MRO). The package also provides automated
  access to HLA alleles frequencies in worldwide human reference
  populations stored in the Allele Frequency Net Database.

- [interacCircos](/packages/interacCircos) Implement in an efficient
  approach to display the genomic data, relationship, information in
  an interactive circular genome(Circos) plot. 'interacCircos' are
  inspired by 'circosJS', 'BioCircos.js' and 'NG-Circos' and we
  integrate the modules of 'circosJS', 'BioCircos.js' and 'NG-Circos'
  into this R package, based on 'htmlwidgets' framework.

- [InteractiveComplexHeatmap](/packages/InteractiveComplexHeatmap)
  This package can easily make heatmaps which are produced by the
  ComplexHeatmap package into interactive applications. It provides
  two types of interactivities: 1. on the interactive graphics
  device, and 2. on a Shiny app. It also provides functions for
  integrating the interactive heatmap widgets for more complex Shiny
  app development.

- [InterCellar](/packages/InterCellar) InterCellar is implemented as
  an R/Bioconductor Package containing a Shiny app that allows users
  to interactively analyze cell-cell communication from scRNA-seq
  data. Starting from precomputed ligand-receptor interactions,
  InterCellar provides filtering options, annotations and multiple
  visualizations to explore clusters, genes and functions. Finally,
  the user can define interaction-pairs modules and link them to
  significant functional terms from Pathways or Gene Ontology.

- [IRISFGM](/packages/IRISFGM) Single-cell RNA-Seq data is useful in
  discovering cell heterogeneity and signature genes in specific cell
  populations in cancer and other complex diseases. Specifically, the
  investigation of functional gene modules (FGM) can help to
  understand gene interactive networks and complex biological
  processes. QUBIC2 is recognized as one of the most efficient and
  effective tools for FGM identification from scRNA-Seq data.
  However, its availability is limited to a C implementation, and its
  applicative power is affected by only a few downstream analyses
  functionalities. We developed an R package named IRIS-FGM
  (integrative scRNA-Seq interpretation system for functional gene
  module analysis) to support the investigation of FGMs and cell
  clustering using scRNA-Seq data. Empowered by QUBIC2, IRIS-FGM can
  identify co-expressed and co-regulated FGMs, predict
  types/clusters, identify differentially expressed genes, and
  perform functional enrichment analysis. It is noteworthy that
  IRIS-FGM also applies Seurat objects that can be easily used in the
  Seurat vignettes.

- [KBoost](/packages/KBoost) Reconstructing gene regulatory networks
  and transcription factor activity is crucial to understand
  biological processes and holds potential for developing
  personalized treatment. Yet, it is still an open problem as
  state-of-art algorithm are often not able to handle large amounts
  of data. Furthermore, many of the present methods predict numerous
  false positives and are unable to integrate other sources of
  information such as previously known interactions. Here we
  introduce KBoost, an algorithm that uses kernel PCA regression,
  boosting and Bayesian model averaging for fast and accurate
  reconstruction of gene regulatory networks. KBoost can also use a
  prior network built on previously known transcription factor
  targets. We have benchmarked KBoost using three different datasets
  against other high performing algorithms. The results show that our
  method compares favourably to other methods across datasets.

- [lisaClust](/packages/lisaClust) lisaClust provides a series of
  functions to identify and visualise regions of tissue where spatial
  associations between cell-types is similar. This package can be
  used to provide a high-level summary of cell-type colocalization in
  multiplexed imaging data that has been segmented at a single-cell
  resolution.

- [LRcell](/packages/LRcell) The goal of LRcell is to identify
  specific sub-cell types that drives the changes observed in a bulk
  RNA-seq differential gene expression experiment. To achieve this,
  LRcell utilizes sets of cell marker genes acquired from single-cell
  RNA-sequencing (scRNA-seq) as indicators for various cell types in
  the tissue of interest. Next, for each cell type, using its marker
  genes as indicators, we apply Logistic Regression on the complete
  set of genes with differential expression p-values to calculate a
  cell-type significance p-value. Finally, these p-values are
  compared to predict which one(s) are likely to be responsible for
  the differential gene expression pattern observed in the bulk
  RNA-seq experiments. LRcell is inspired by the
  LRpath[@sartor2009lrpath] algorithm developed by Sartor et al.,
  originally designed for pathway/gene set enrichment analysis.
  LRcell contains three major components: LRcell analysis, plot
  generation and marker gene selection. All modules in this package
  are written in R. This package also provides marker genes in the
  Prefrontal Cortex (pFC) human brain region, human PBMC and nine
  mouse brain regions (Frontal Cortex, Cerebellum, Globus Pallidus,
  Hippocampus, Entopeduncular, Posterior Cortex, Striatum, Substantia
  Nigra and Thalamus).

- [MACSr](/packages/MACSr) The Model-based Analysis of ChIP-Seq
  (MACS) is a widely used toolkit for identifying transcript factor
  binding sites. This package is an R wrapper of the lastest MACS3.

- [MAGAR](/packages/MAGAR) "Methylation-Aware Genotype Association in
  R" (MAGAR) computes methQTL from DNA methylation and genotyping
  data from matched samples. MAGAR uses a linear modeling stragety to
  call CpGs/SNPs that are methQTLs. MAGAR accounts for the local
  correlation structure of CpGs.

- [MatrixQCvis](/packages/MatrixQCvis) Data quality assessment is an
  integral part of preparatory data analysis to ensure sound
  biological information retrieval. We present here the MatrixQCvis
  package, which provides shiny-based interactive visualization of
  data quality metrics at the per-sample and per-feature level. It is
  broadly applicable to quantitative omics data types that come in
  matrix-like format (features x samples). It enables the detection
  of low-quality samples, drifts, outliers and batch effects in data
  sets. Visualizations include amongst others bar- and violin plots
  of the (count/intensity) values, mean vs standard deviation plots,
  MA plots, empirical cumulative distribution function (ECDF) plots,
  visualizations of the distances between samples, and multiple types
  of dimension reduction plots. Furthermore, MatrixQCvis allows for
  differential expression analysis based on the limma (moderated
  t-tests) and proDA (Wald tests) packages. MatrixQCvis builds upon
  the popular Bioconductor SummarizedExperiment S4 class and enables
  thus the facile integration into existing workflows. The package is
  especially tailored towards metabolomics and proteomics mass
  spectrometry data, but also allows to assess the data quality of
  other data types that can be represented in a SummarizedExperiment
  object.

- [memes](/packages/memes) A seamless interface to the MEME Suite
  family of tools for motif analysis. 'memes' provides data aware
  utilities for using GRanges objects as entrypoints to motif
  analysis, data structures for examining & editing motif lists, and
  novel data visualizations. 'memes' functions and data structures
  are amenable to both base R and tidyverse workflows.

- [MetaboCoreUtils](/packages/MetaboCoreUtils) MetaboCoreUtils
  defines metabolomics-related core functionality provided as
  low-level functions to allow a data structure-independent usage
  across various R packages. This includes functions to calculate
  between ion (adduct) and compound mass-to-charge ratios and masses
  or functions to work with chemical formulas. The package provides
  also a set of adduct definitions and information on some
  commercially available internal standard mixes commonly used in MS
  experiments.

- [metapod](/packages/metapod) Implements a variety of methods for
  combining p-values in differential analyses of genome-scale
  datasets. Functions can combine p-values across different tests in
  the same analysis (e.g., genomic windows in ChIP-seq, exons in
  RNA-seq) or for corresponding tests across separate analyses (e.g.,
  replicated comparisons, effect of different treatment conditions).
  Support is provided for handling log-transformed input p-values,
  missing values and weighting where appropriate.

- [methylscaper](/packages/methylscaper) methylscaper is an R package
  for processing and visualizing data jointly profiling methylation
  and chromatin accessibility (MAPit, NOMe-seq, scNMT-seq, nanoNOMe,
  etc.). The package supports both single-cell and single-molecule
  data, and a common interface for jointly visualizing both data
  types through the generation of ordered representational
  methylation-state matrices. The Shiny app allows for an interactive
  seriation process of refinement and re-weighting that optimally
  orders the cells or DNA molecules to discover methylation patterns
  and nucleosome positioning.

- [mia](/packages/mia) mia implements tools for microbiome analysis
  based on the SummarizedExperiment, SingleCellExperiment and
  TreeSummarizedExperiment infrastructure. Data wrangling and
  analysis in the context of taxonomic data is the main scope.
  Additional functions for common task are implemented such as
  community indices calculation and summarization.

- [miaViz](/packages/miaViz) miaViz implements plotting function to
  work with TreeSummarizedExperiment and related objects in a context
  of microbiome analysis. Among others this includes plotting tree,
  graph and microbiome series data.

- [midasHLA](/packages/midasHLA) MiDAS is a R package for
  immunogenetics data transformation and statistical analysis. MiDAS
  accepts input data in the form of HLA alleles and KIR types, and
  can transform it into biologically meaningful variables, enabling
  HLA amino acid fine mapping, analyses of HLA evolutionary
  divergence, KIR gene presence, as well as validated HLA-KIR
  interactions. Further, it allows comprehensive statistical
  association analysis workflows with phenotypes of diverse
  measurement scales. MiDAS closes a gap between the inference of
  immunogenetic variation and its efficient utilization to make
  relevant discoveries related to T cell, Natural Killer cell, and
  disease biology.

- [miloR](/packages/miloR) This package performs single-cell
  differential abundance testing. Cell states are modelled as
  representative neighbourhoods on a nearest neighbour graph.
  Hypothesis testing is performed using a negative bionomial
  generalized linear model.

- [mina](/packages/mina) An increasing number of microbiome datasets
  have been generated and analyzed with the help of rapidly
  developing sequencing technologies. At present, analysis of
  taxonomic profiling data is mainly conducted using
  composition-based methods, which ignores interactions between
  community members. Besides this, a lack of efficient ways to
  compare microbial interaction networks limited the study of
  community dynamics. To better understand how community diversity is
  affected by complex interactions between its members, we developed
  a framework (Microbial community dIversity and Network Analysis,
  mina), a comprehensive framework for microbial community diversity
  analysis and network comparison. By defining and integrating
  network-derived community features, we greatly reduce
  noise-to-signal ratio for diversity analyses. A bootstrap and
  permutation-based method was implemented to assess community
  network dissimilarities and extract discriminative features in a
  statistically principled way.

- [miQC](/packages/miQC) Single-cell RNA-sequencing (scRNA-seq) has
  made it possible to profile gene expression in tissues at high
  resolution.  An important preprocessing step prior to performing
  downstream analyses is to identify and remove cells with poor or
  degraded sample quality using quality control (QC) metrics.  Two
  widely used QC metrics to identify a ‘low-quality’ cell are (i) if
  the cell includes a high proportion of reads that map to
  mitochondrial DNA encoded genes (mtDNA) and (ii) if a small number
  of genes are detected. miQC is data-driven QC metric that jointly
  models both the proportion of reads mapping to mtDNA and the number
  of detected genes with mixture models in a probabilistic framework
  to predict the low-quality cells in a given dataset.

- [mirTarRnaSeq](/packages/mirTarRnaSeq) mirTarRnaSeq R package can
  be used for interactive mRNA miRNA sequencing statistical analysis.
  This package utilizes expression or differential expression mRNA
  and miRNA sequencing results and performs interactive correlation
  and various GLMs (Regular GLM, Multivariate GLM, and Interaction
  GLMs ) analysis between mRNA and miRNA expriments. These
  experiments can be time point experiments, and or condition
  expriments.

- [mistyR](/packages/mistyR) mistyR is an impolementation of the
  Multiview Intercellular SpaTialmodeling framework (MISTy). MISTy is
  an explainable machine learning framework for knowledge extraction
  and analysis of single-cell, highly multiplexed, spatially resolved
  data. MISTy facilitates an in-depth understanding of marker
  interactions by profiling the intra- and intercellular
  relationships. MISTy is a flexible framework able to process a
  custom number of views. Each of these views can describe a
  different spatial context, i.e., define a relationship among the
  observed expressions of the markers, such as intracellular
  regulation or paracrine regulation, but also, the views can also
  capture cell-type specific relationships, capture relations between
  functional footprints or focus on relations between different
  anatomical regions. Each MISTy view is considered as a potential
  source of variability in the measured marker expressions. Each
  MISTy view is then analyzed for its contribution to the total
  expression of each marker and is explained in terms of the
  interactions with other measurements that led to the observed
  contribution.

- [moanin](/packages/moanin) Simple and efficient workflow for
  time-course gene expression data, built on publictly available
  open-source projects hosted on CRAN and bioconductor. moanin
  provides helper functions for all the steps required for analysing
  time-course data using functional data analysis: (1) functional
  modeling of the timecourse data; (2) differential expression
  analysis; (3) clustering; (4) downstream analysis.

- [ModCon](/packages/ModCon) Collection of functions to calculate a
  nucleotide sequence surrounding for splice donors sites to either
  activate or repress donor usage. The proposed alternative
  nucleotide sequence encodes the same amino acid and could be
  applied e.g. in reporter systems to silence or activate cryptic
  splice donor sites.

- [MQmetrics](/packages/MQmetrics) The package MQmetrics (MaxQuant
  metrics) provides a workflow to analyze the quality and
  reproducibility of your proteomics mass spectrometry analysis from
  MaxQuant.Input data are extracted from several MaxQuant output
  tables, and produces a pdf report. It includes several
  visualization tools to check numerous parameters regarding the
  quality of the runs.It also includes two functions to visualize the
  iRT peptides from Biognosysin case they were spiked in the samples.

- [MsBackendMassbank](/packages/MsBackendMassbank) Mass spectrometry
  (MS) data backend supporting import and export of MS/MS library
  spectra from MassBank record files. Different backends are
  available that allow handling of data in plain MassBank text file
  format or allow also to interact directly with MassBank SQL
  databases. Objects from this package are supposed to be used with
  the Spectra Bioconductor package. This package thus adds MassBank
  support to the Spectra package.

- [MsBackendMgf](/packages/MsBackendMgf) Mass spectrometry (MS) data
  backend supporting import and export of MS/MS spectra data from
  Mascot Generic Format (mgf) files. Objects defined in this package
  are supposed to be used with the Spectra Bioconductor package. This
  package thus adds mgf file support to the Spectra package.

- [MsFeatures](/packages/MsFeatures) The MsFeature package defines
  functionality for Mass Spectrometry features. This includes
  functions to group (LC-MS) features based on some of their
  properties, such as retention time (coeluting features), or
  correlation of signals across samples. This packge hence allows to
  group features, and its results can be used as an input for the
  `QFeatures` package which allows to aggregate abundance levels of
  features within each group. This package defines concepts and
  functions for base and common data types, implementations for more
  specific data types are expected to be implemented in the
  respective packages (such as e.g. `xcms`). All functionality of
  this package is implemented in a modular way which allows
  combination of different grouping approaches and enables its re-use
  in other R packages.

- [msqrob2](/packages/msqrob2) msqrob2 provides a robust linear mixed
  model framework for assessing differential abundance in MS-based
  Quantitative proteomics experiments. Our workflows can start from
  raw peptide intensities or summarised protein expression values.
  The model parameter estimates can be stabilized by ridge
  regression, empirical Bayes variance estimation and robust
  M-estimation. msqrob2's hurde workflow can handle missing data
  without having to rely on hard-to-verify imputation assumptions,
  and, outcompetes state-of-the-art methods with and without
  imputation for both high and low missingness. It builds on QFeature
  infrastructure for quantitative mass spectrometry data to store the
  model results together with the raw data and preprocessed data.

- [MSstatsLOBD](/packages/MSstatsLOBD) The MSstatsLOBD package allows
  calculation and visualization of limit of blac (LOB) and limit of
  detection (LOD). We define the LOB as the highest apparent
  concentration of a peptide expected when replicates of a blank
  sample containing no peptides are measured. The LOD is defined as
  the measured concentration value for which the probability of
  falsely claiming the absence of a peptide in the sample is 0.05,
  given a probability 0.05 of falsely claiming its presence. These
  functionalities were previously a part of the MSstats package. The
  methodology is described in Galitzine (2018)
  <doi:10.1074/mcp.RA117.000322>.

- [multiSight](/packages/multiSight) multiSight is an R package
  providing an user-friendly graphical interface to analyze your omic
  datasets in a multi-omics manner based on Stouffer's p-value
  pooling and multi-block statistical methods. For each omic dataset
  you furnish, multiSight provides classification models with feature
  selection you can use as biosignature: (i) To forecast phenotypes
  (e.g. to diagnostic tasks, histological subtyping), (ii) To design
  Pathways and gene ontology enrichments (Over Representation
  Analysis), (iii) To build Network inference linked to PubMed
  querying to make assumptions easier and data-driven.

- [mumosa](/packages/mumosa) Assorted utilities for multi-modal
  analyses of single-cell datasets. Includes functions to combine
  multiple modalities for downstream analysis, perform MNN-based
  batch correction across multiple modalities, and to compute
  correlations between assay values for different modalities.

- [MungeSumstats](/packages/MungeSumstats) The *MungeSumstats*
  package is designed to facilitate the standardisation of GWAS
  summary statistics. It reformats inputted summary statisitics to
  include SNP, CHR, BP and can look up these values if any are
  missing. It also removes duplicates across SNPs.

- [NanoStringNCTools](/packages/NanoStringNCTools) Tools for
  NanoString Technologies nCounter Technology.  Provides support for
  reading RCC files into an ExpressionSet derived object.  Also
  includes methods for QC and normalizaztion of NanoString data.

- [nempi](/packages/nempi) Takes as input an incomplete perturbation
  profile and differential gene expression in log odds and infers
  unobserved perturbations and augments observed ones. The inference
  is done by iteratively inferring a network from the perturbations
  and inferring perturbations from the network. The network inference
  is done by Nested Effects Models.

- [ORFhunteR](/packages/ORFhunteR) The ORFhunteR package is a R and
  C++ library for an automatic determination and annotation of open
  reading frames (ORF) in a large set of RNA molecules. It
  efficiently implements the machine learning model based on
  vectorization of nucleotide sequences and the random forest
  classification algorithm. The ORFhunteR package consists of a set
  of functions written in the R language in conjunction with C++. The
  efficiency of the package was confirmed by the examples of the
  analysis of RNA molecules from the NCBI RefSeq and Ensembl
  databases. The package can be used in basic and applied biomedical
  research related to the study of the transcriptome of normal as
  well as altered (for example, cancer) human cells.

- [PDATK](/packages/PDATK) Pancreatic ductal adenocarcinoma (PDA) has
  a relatively poor prognosis and is one of the most lethal cancers.
  Molecular classification of gene expression profiles holds the
  potential to identify meaningful subtypes which can inform
  therapeutic strategy in the clinical setting. The Pancreatic Cancer
  Adenocarcinoma Tool-Kit (PDATK) provides an S4 class-based
  interface for performing unsupervised subtype discovery,
  cross-cohort meta-clustering, gene-expression-based classification,
  and subsequent survival analysis to identify prognostically useful
  subtypes in pancreatic cancer and beyond. Two novel methods,
  Consensus Subtypes in Pancreatic Cancer (CSPC) and Pancreatic
  Cancer Overall Survival Predictor (PCOSP) are included for
  consensus-based meta-clustering and overall-survival prediction,
  respectively. Additionally, four published subtype classifiers and
  three published prognostic gene signatures are included to allow
  users to easily recreate published results, apply existing
  classifiers to new data, and benchmark the relative performance of
  new methods. The use of existing Bioconductor classes as input to
  all PDATK classes and methods enables integration with existing
  Bioconductor datasets, including the 21 pancreatic cancer patient
  cohorts available in the MetaGxPancreas data package. PDATK has
  been used to replicate results from Sandhu et al (2019)
  [https://doi.org/10.1200/cci.18.00102] and an additional paper is
  in the works using CSPC to validate subtypes from the included
  published classifiers, both of which use the data available in
  MetaGxPancreas. The inclusion of subtype centroids and prognostic
  gene signatures from these and other publications will enable
  researchers and clinicians to classify novel patient gene
  expression data, allowing the direct clinical application of the
  classifiers included in PDATK. Overall, PDATK provides a rich set
  of tools to identify and validate useful prognostic and molecular
  subtypes based on gene-expression data, benchmark new classifiers
  against existing ones, and apply discovered classifiers on novel
  patient data to inform clinical decision making.

- [PFP](/packages/PFP) An implementation of the pathway fingerprint
  framework that introduced in paper "Pathway Fingerprint: a novel
  pathway knowledge and topology based method for biomarker discovery
  and characterization". This method provides a systematic
  comparisons between a gene set (such as a list of differentially
  expressed genes) and well-studied "basic pathway networks" (KEGG
  pathways), measuring the importance of pathways and genes for the
  gene set. The package is helpful for researchers to find the
  biomarkers and its function.

- [PhenoGeneRanker](/packages/PhenoGeneRanker) This package is a
  gene/phenotype prioritization tool that utilizes multiplex
  heterogeneous gene phenotype network. PhenoGeneRanker allows
  multi-layer gene and phenotype networks. It also calculates
  empirical p-values of gene/phenotype ranking using random
  stratified sampling of genes/phenotypes based on their connectivity
  degree in the network.
  https://dl.acm.org/doi/10.1145/3307339.3342155.

- [PhIPData](/packages/PhIPData) PhIPData defines an S4 class for
  phage-immunoprecipitation sequencing (PhIP-seq) experiments.
  Buliding upon the RangedSummarizedExperiment class, PhIPData
  enables users to coordinate metadata with experimental data in
  analyses. Additionally, PhIPData provides specialized methods to
  subset and identify beads-only samples, subset objects using virus
  aliases, and use existing peptide libraries to populate object
  parameters.

- [planet](/packages/planet) This package contains R functions to
  infer additional biological variables to supplemental DNA
  methylation analysis of placental data. This includes inferring
  ethnicity/ancestry, gestational age, and cell composition from
  placental DNA methylation array (450k/850k) data. The package comes
  with an example processed placental dataset.

- [PoDCall](/packages/PoDCall) Reads files exported from 'QuantaSoft'
  containing amplitude values from a run of ddPCR (96 well plate) and
  robustly sets thresholds to determine positive droplets for each
  channel of each individual well. Concentration and normalized
  concentration in addition to other metrics is then calculated for
  each well. Results are returned as a table, optionally written to
  file, as well as optional plots (scatterplot and histogram) for
  both channels per well written to file. The package includes a
  shiny application which provides an interactive and user-friendly
  interface to the full functionality of PoDCall.

- [POWSC](/packages/POWSC) Determining the sample size for adequate
  power to detect statistical significance is a crucial step at the
  design stage for high-throughput experiments. Even though a number
  of methods and tools are available for sample size calculation for
  microarray and RNA-seq in the context of differential expression
  (DE), this topic in the field of single-cell RNA sequencing is
  understudied. Moreover, the unique data characteristics present in
  scRNA-seq such as sparsity and heterogeneity increase the
  challenge. We propose POWSC, a simulation-based method, to provide
  power evaluation and sample size recommendation for single-cell RNA
  sequencing DE analysis. POWSC consists of a data simulator that
  creates realistic expression data, and a power assessor that
  provides a comprehensive evaluation and visualization of the power
  and sample size relationship.

- [ppcseq](/packages/ppcseq) Relative transcript abundance has proven
  to be a valuable tool for understanding the function of genes in
  biological systems. For the differential analysis of transcript
  abundance using RNA sequencing data, the negative binomial model is
  by far the most frequently adopted. However, common methods that
  are based on a negative binomial model are not robust to extreme
  outliers, which we found to be abundant in public datasets. So far,
  no rigorous and probabilistic methods for detection of outliers
  have been developed for RNA sequencing data, leaving the
  identification mostly to visual inspection. Recent advances in
  Bayesian computation allow large-scale comparison of observed data
  against its theoretical distribution given in a statistical model.
  Here we propose ppcseq, a key quality-control tool for identifying
  transcripts that include outlier data points in differential
  expression analysis, which do not follow a negative binomial
  distribution. Applying ppcseq to analyse several publicly available
  datasets using popular tools, we show that from 3 to 10 percent of
  differentially abundant transcripts across algorithms and datasets
  had statistics inflated by the presence of outliers.

- [ptairMS](/packages/ptairMS) This package implements a suite of
  methods to preprocess data from PTR-TOF-MS instruments (HDF5
  format) and generates the 'sample by features' table of peak
  intensities in addition to the sample and feature metadata (as a
  single ExpressionSet object for subsequent statistical analysis).
  This package also permit usefull tools for cohorts management as
  analyzing data progressively, visualization tools and quality
  control. The steps include calibration, expiration detection, peak
  detection and quantification, feature alignment, missing value
  imputation and feature annotation. Applications to exhaled air and
  cell culture in headspace are described in the vignettes and
  examples. This package was used for data analysis of Gassin Delyle
  study on adults undergoing invasive mechanical ventilation in the
  intensive care unit due to severe COVID-19 or non-COVID-19 acute
  respiratory distress syndrome (ARDS), and permit to identfy four
  potentiel biomarquers of the infection.

- [quantiseqr](/packages/quantiseqr) This package provides a
  streamlined workflow for the quanTIseq method, developed to perform
  the quantification of the Tumor Immune contexture from RNA-seq
  data. The quantification is performed against the TIL10 signature
  (dissecting the contributions of ten immune cell types), carefully
  crafted from a collection of human RNA-seq samples. The TIL10
  signature has been extensively validated using simulated, flow
  cytometry, and immunohistochemistry data.

- [ramr](/packages/ramr) ramr is an R package for detection of
  low-frequency aberrant methylation events in large datasets
  obtained by methylation profiling using array or high-throughput
  bisulfite sequencing. In addition, package provides functions to
  visualize found aberrantly methylated regions (AMRs), and to
  generate sets of all possible regions to be used as reference sets
  for enrichment analysis.

- [rawrr](/packages/rawrr) This package wraps the functionality of
  the RawFileReader .NET assembly. Within the R environment, spectra
  and chromatograms are represented by S3 objects (Kockmann T. et al. (2020) 
  <doi:10.1101/2020.10.30.362533>). The package provides basic
  functions to download and install the required third-party
  libraries. The package is developed, tested, and used at the
  Functional Genomics Center Zurich, Switzerland <https://fgcz.ch>.

- [Rbec](/packages/Rbec) Rbec is a adapted version of DADA2 for
  analyzing amplicon sequencing data from synthetic communities
  (SynComs), where the reference sequences for each strain exists.
  Rbec can not only accurately profile the microbial compositions in
  SynComs, but also predict the contaminants in SynCom samples.

- [RCSL](/packages/RCSL) A novel clustering algorithm and toolkit
  RCSL (Rank Constrained Similarity Learning) to accurately identify
  various cell types using scRNA-seq data from a complex tissue. RCSL
  considers both lo-cal similarity and global similarity among the
  cells to discern the subtle differences among cells of the same
  type as well as larger differences among cells of different types.
  RCSL uses Spearman’s rank correlations of a cell’s expression
  vector with those of other cells to measure its global similar-ity,
  and adaptively learns neighbour representation of a cell as its
  local similarity. The overall similar-ity of a cell to other cells
  is a linear combination of its global similarity and local
  similarity.

- [ReactomeContentService4R](/packages/ReactomeContentService4R)
  Reactome is a free, open-source, open access, curated and
  peer-reviewed knowledgebase of bio-molecular pathways. This package
  is to interact with the Reactome Content Service API. Pre-built
  functions would allow users to retrieve data and images that
  consist of proteins, pathways, and other molecules related to a
  specific gene or entity in Reactome.

- [ReactomeGraph4R](/packages/ReactomeGraph4R) Pathways, reactions,
  and biological entities in Reactome knowledge are systematically
  represented as an ordered network. Instances are represented as
  nodes and relationships between instances as edges; they are all
  stored in the Reactome Graph Database. This package serves as an
  interface to query the interconnected data from a local Neo4j
  database, with the aim of minimizing the usage of Neo4j Cypher
  queries.

- [RiboDiPA](/packages/RiboDiPA) This package performs differential
  pattern analysis for Ribo-seq data. It identifies genes with
  significantly different patterns in the ribosome footprint between
  two conditions. RiboDiPA contains five major components including
  bam file processing, P-site mapping, data binning, differential
  pattern analysis and footprint visualization.

- [RLassoCox](/packages/RLassoCox) RLassoCox is a package that
  implements the RLasso-Cox model proposed by Wei Liu. The RLasso-Cox
  model integrates gene interaction information into the Lasso-Cox
  model for accurate survival prediction and survival biomarker
  discovery. It is based on the hypothesis that topologically
  important genes in the gene interaction network tend to have stable
  expression changes. The RLasso-Cox model uses random walk to
  evaluate the topological weight of genes, and then highlights
  topologically important genes to improve the generalization ability
  of the Lasso-Cox model. The RLasso-Cox model has the advantage of
  identifying small gene sets with high prognostic performance on
  independent datasets, which may play an important role in
  identifying robust survival biomarkers for various cancer types.

- [SANTA](/packages/SANTA) This package provides methods for
  measuring the strength of association between a network and a
  phenotype. It does this by measuring clustering of the phenotype
  across the network (Knet). Vertices can also be individually ranked
  by their strength of association with high-weight vertices (Knode).

- [satuRn](/packages/satuRn) satuRn provides a higly performant and
  scalable framework for performing differential transcript usage
  analyses. The package consists of three main functions. The first
  function, fitDTU, fits quasi-binomial generalized linear models
  that model transcript usage in different groups of interest. The
  second function, testDTU, tests for differential usage of
  transcripts between groups of interest. Finally, plotDTU visualizes
  the usage profiles of transcripts in groups of interest.

- [ScaledMatrix](/packages/ScaledMatrix) Provides delayed computation
  of a matrix of scaled and centered values. The result is equivalent
  to using the scale() function but avoids explicit realization of a
  dense matrix during block processing. This permits greater
  efficiency in common operations, most notably matrix
  multiplication.

- [SCArray](/packages/SCArray) Provides large-scale single-cell
  RNA-seq data manipulation using Genomic Data Structure (GDS) files.
  It combines dense and sparse matrices stored in GDS files and the
  Bioconductor infrastructure framework (SingleCellExperiment and
  DelayedArray) to provide out-of-memory data storage and large-scale
  manipulation using the R programming language.

- [scClassifR](/packages/scClassifR) The package comprises a set of
  pretrained machine learning models to predict basic immune cell
  types. This enables all users to quickly get a first annotation of
  the cell types present in their dataset without requiring prior
  knowledge. scClassifR also allows users to train their own models
  to predict new cell types based on specific research needs.

- [sechm](/packages/sechm) sechm provides a simple interface between
  SummarizedExperiment objects and the ComplexHeatmap package. It
  enables plotting annotated heatmaps from SE objects, with easy
  access to rowData and colData columns, and implements a number of
  features to make the generation of heatmaps easier and more
  flexible. These functionalities used to be part of the SEtools
  package.

- [shinyepico](/packages/shinyepico) ShinyÉPICo is a graphical
  pipeline to analyze Illumina DNA methylation arrays (450k or EPIC).
  It allows to calculate differentially methylated positions and
  differentially methylated regions in a user-friendly interface.
  Moreover, it includes several options to export the results and
  obtain files to perform downstream analysis.

- [SingleMoleculeFootprinting](/packages/SingleMoleculeFootprinting)
  SingleMoleculeFootprinting is an R package providing functions to
  analyze Single Molecule Footprinting (SMF) data. Following the
  workflow exemplified in its vignette, the user will be able to
  perform basic data analysis of SMF data with minimal coding effort.
  Starting from an aligned bam file, we show how to perform quality
  controls over sequencing libraries, extract methylation information
  at the single molecule level accounting for the two possible kind
  of SMF experiments (single enzyme or double enzyme), classify
  single molecules based on their patterns of molecular occupancy,
  plot SMF information at a given genomic location

- [sitadela](/packages/sitadela) Provides an interface to build a
  unified database of genomic annotations and their coordinates
  (gene, transcript and exon levels). It is aimed to be used when
  simple tab-delimited annotations (or simple GRanges objects) are
  required instead of the more complex annotation Bioconductor
  packages. Also useful when combinatorial annotation elements are
  reuired, such as RefSeq coordinates with Ensembl biotypes. Finally,
  it can download, construct and handle annotations with versioned
  genes and transcripts (where available, e.g. RefSeq and latest
  Ensembl). This is particularly useful in precision medicine
  applications where the latter must be reported.

- [SOMNiBUS](/packages/SOMNiBUS) This package aims to analyse
  count-based methylation data on predefined genomic regions, such as
  those obtained by targeted sequencing, and thus to identify
  differentially methylated regions (DMRs) that are associated with
  phenotypes or traits. The method is built a rich flexible model
  that allows for the effects, on the methylation levels, of multiple
  covariates to vary smoothly along genomic regions. At the same
  time, this method also allows for sequencing errors and can adjust
  for variability in cell type mixture.

- [SplicingFactory](/packages/SplicingFactory) The SplicingFactory R
  package uses transcript-level expression values to analyze splicing
  diversity based on various statistical measures, like Shannon
  entropy or the Gini index. These measures can quantify transcript
  isoform diversity within samples or between conditions.
  Additionally, the package analyzes the isoform diversity data,
  looking for significant changes between conditions.

- [Summix](/packages/Summix) This package contains the Summix method
  for estimating and adjusting for ancestry in genetic summary allele
  frequency data. The function summix estimates reference ancestry
  proportions using a mixture model. The adjAF function produces
  ancestry adjusted allele frequencies for an observed sample with
  ancestry proportions matching a target person or sample.

- [supersigs](/packages/supersigs) Generate SuperSigs (supervised
  mutational signatures) from single nucleotide variants in the
  cancer genome. Functions included in the package allow the user to
  learn supervised mutational signatures from their data and apply
  them to new data. The methodology is based on the one described in
  Afsari (2021, ELife).

- [systemPipeTools](/packages/systemPipeTools) systemPipeTools
  package extends the widely used systemPipeR (SPR) workflow
  environment with an enhanced toolkit for data visualization,
  including utilities to automate the data visualizaton for analysis
  of differentially expressed genes (DEGs). systemPipeTools provides
  data transformation and data exploration functions via
  scatterplots, hierarchical clustering heatMaps, principal component
  analysis, multidimensional scaling, generalized principal
  components, t-Distributed Stochastic Neighbor embedding (t-SNE),
  and MA and volcano plots. All these utilities can be integrated
  with the modular design of the systemPipeR environment that allows
  users to easily substitute any of these features and/or custom with
  alternatives.

- [tLOH](/packages/tLOH) tLOH, or transcriptomicsLOH, assesses
  evidence for loss of heterozygosity (LOH) in pre-processed spatial
  transcriptomics data. This tool requires spatial transcriptomics
  cluster and allele count information at likely heterozygous
  single-nucleotide polymorphism (SNP) positions in VCF format. Bayes
  factors are calculated at each SNP to determine likelihood of
  potential loss of heterozygosity event. Two plotting functions are
  included to visualize allele fraction and aggregated Bayes factor
  per chromosome. Data generated with the 10X Genomics Visium Spatial
  Gene Expression platform must be pre-processed to obtain an
  individual sample VCF with columns for each cluster. Required
  fields are allele depth (AD) with counts for reference/alternative
  alleles and read depth (DP).

- [TrajectoryGeometry](/packages/TrajectoryGeometry) Given a time
  series or pseudo-times series of gene expression data, we might
  wish to know: Do the changes in gene expression in these data
  exhibit directionality?  Are there turning points in this
  directionality.  Do different subsets of the data move in different
  directions?  This package uses spherical geometry to probe these
  sorts of questions.  In particular, if we are looking at (say) the
  first n dimensions of the PCA of gene expression, directionality
  can be detected as the clustering of points on the
  (n-1)-dimensional sphere.

- [TrajectoryUtils](/packages/TrajectoryUtils) Implements low-level
  utilities for single-cell trajectory analysis, primarily intended
  for re-use inside higher-level packages. Include a function to
  create a cluster-level minimum spanning tree and data structures to
  hold pseudotime inference results.

- [TraRe](/packages/TraRe) TraRe (Transcriptional Rewiring) is an R
  package which contains the necessary tools to carry out several
  functions. Identification of module-based gene regulatory networks
  (GRN); score-based classification of these modules via a rewiring
  test; visualization of rewired modules to analyze condition-based
  GRN deregulation and drop out genes recovering via cliques
  methodology. For each tool, an html report can be generated
  containing useful information about the generated GRN and
  statistical data about the performed tests. These tools have been
  developed considering sequenced data (RNA-Seq).

- [Travel](/packages/Travel) Creates a virtual pointer for R's ALTREP
  object which does not have the data allocates in memory. The
  pointer is made by the file mapping of a virtual file so it behaves
  exactly the same as a regular pointer. All the requests to access
  the pointer will be sent to the underlying file system and
  eventually handled by a customized data-reading function. The main
  purpose of the package is to reduce the memory consumption when
  using R's vector to represent a large data. The use cases of the
  package include on-disk data representation, compressed vector(e.g.
  RLE) and etc.

- [treekoR](/packages/treekoR) treekoR is a novel framework that aims
  to utilise the hierarchical nature of single cell cytometry data to
  find robust and interpretable associations between cell subsets and
  patient clinical end points. These associations are aimed to
  recapitulate the nested proportions prevalent in workflows
  inovlving manual gating, which are often overlooked in workflows
  using automatic clustering to identify cell populations. We
  developed treekoR to: Derive a hierarchical tree structure of cell
  clusters; measure the proportions to parent (proportions of cells
  each node in the tree relative to the number of cells belonging its
  parent node), in addition to the proportions to all (proportion of
  cells in each node relative to all cells); perform significance
  testing using the calculated proportions; and provide an
  interactive html visualisation to help highlight key results.

- [tricycle](/packages/tricycle) The package contains functions to
  infer and visualize cell cycle process using Single Cell RNASeq
  data. It exploits the idea of transfer learning, projecting new
  data to the previous learned biologically interpretable space. We
  provide a pre-learned cell cycle space, which could be used to
  infer cell cycle time of human and mouse single cell samples. In
  addition, we also offer functions to visualize cell cycle time on
  different embeddings and functions to build new reference.

- [ttgsea](/packages/ttgsea) Functional enrichment analysis methods
  such as gene set enrichment analysis (GSEA) have been widely used
  for analyzing gene expression data. GSEA is a powerful method to
  infer results of gene expression data at a level of gene sets by
  calculating enrichment scores for predefined sets of genes. GSEA
  depends on the availability and accuracy of gene sets. There are
  overlaps between terms of gene sets or categories because multiple
  terms may exist for a single biological process, and it can thus
  lead to redundancy within enriched terms. In other words, the sets
  of related terms are overlapping. Using deep learning, this pakage
  is aimed to predict enrichment scores for unique tokens or words
  from text in names of gene sets to resolve this overlapping set
  issue. Furthermore, we can coin a new term by combining tokens and
  find its enrichment score by predicting such a combined tokens.

- [VarCon](/packages/VarCon) VarCon is an R package which converts
  the positional information from the annotation of an single
  nucleotide variation (SNV) (either referring to the coding sequence
  or the reference genomic sequence). It retrieves the genomic
  reference sequence around the position of the single nucleotide
  variation. To asses, whether the SNV could potentially influence
  binding of splicing regulatory proteins VarCon calcualtes the
  HEXplorer score as an estimation. Besides, VarCon additionally
  reports splice site strengths of splice sites within the retrieved
  genomic sequence and any changes due to the SNV.

- [vissE](/packages/vissE) This package enables the interpretation
  and analysis of results from a gene set enrichment analysis using
  network-based and text-mining approaches. Most enrichment analyses
  result in large lists of significant gene sets that are difficult
  to interpret. Tools in this package help build a similarity-based
  network of significant gene sets from a gene set enrichment
  analysis that can then be investigated for their biological
  function using text-mining approaches.

- [wppi](/packages/wppi) Protein-protein interaction data is
  essential for omics data analysis and modeling. Database knowledge
  is general, not specific for cell type, physiological condition or
  any other context determining which connections are functional and
  contribute to the signaling. Functional annotations such as Gene
  Ontology and Human Phenotype Ontology might help to evaluate the
  relevance of interactions. This package predicts functional
  relevance of protein-protein interactions based on functional
  annotations such as Human Protein Ontology and Gene Ontology, and
  prioritizes genes based on network topology, functional scores and
  a path search algorithm.

- [XNAString](/packages/XNAString) The XNAString package allows for
  description of base sequences and associated chemical modifications
  in a single object. XNAString is able to capture single stranded,
  as well as double stranded molecules. Chemical modifications are
  represented as independent strings associated with different
  features of the molecules (base sequence, sugar sequence, backbone
  sequence, modifications) and can be read or written to a HELM
  notation. It also enables secondary structure prediction using
  RNAfold from ViennaRNA. XNAString is designed to be efficient
  representation of nucleic-acid based therapeutics, therefore it
  stores information about target sequences and provides interface
  for matching and alignment functions from Biostrings package.

New Data Experiment Packages
=====================

There are 22 new data experiment packages in this release of Bioconductor.

- [BeadSorted.Saliva.EPIC](/packages/BeadSorted.Saliva.EPIC) Raw data
  objects used to estimate saliva cell proportion estimates in
  ewastools. The FlowSorted.Saliva.EPIC object is constructed from
  saples assayed by Lauren Middleton et. al. (2021).

- [BioImageDbs](/packages/BioImageDbs) The package provides a
  bioimage dataset for the image analysis using machine learning and
  deep learning. The dataset includes microscopy imaging data with
  supervised labels. The data is provided as R list data that can be
  loaded to Keras/tensorflow in R.

- [DExMAdata](/packages/DExMAdata) Data objects needed to allSameID()
  function of DExMA package. There are also some objects that are
  necessary to be able to apply the examples of the DExMA package,
  which illustrate package functionality.

- [emtdata](/packages/emtdata) This package provides pre-processed
  RNA-seq data where the epithelial to mesenchymal transition was
  induced on cell lines. These data come from three publications
  Cursons et al. (2015), Cursons etl al. (2018) and Foroutan et al.
  (2017). In each of these publications, EMT was induces across
  multiple cell lines following treatment by TGFb among other
  stimulants. This data will be useful in determining the regulatory
  programs modified in order to achieve an EMT. Data were processed
  by the Davis laboratory in the Bioinformatics division at WEHI.

- [ewceData](/packages/ewceData) This package provides reference data
  required for ewce. Expression Weighted Celltype Enrichment (EWCE)
  is used to determine which cell types are enriched within gene
  lists. The package provides tools for testing enrichments within
  simple gene lists (such as human disease associated genes) and
  those resulting from differential expression studies. The package
  does not depend upon any particular Single Cell Transcriptome
  dataset and user defined datasets can be loaded in and used in the
  analyses.

- [GenomicDistributionsData](/packages/GenomicDistributionsData) This
  package provides ready to use reference data for
  GenomicDistributions package. Raw data was obtained from ensembldb
  and processed with helper functions. Data files are available for
  the following genome assemblies: hg19, hg38, mm9 and mm10.

- [GSE13015](/packages/GSE13015) Microarray expression matrix
  platform GPL6106 and clinical data for 67 septicemic patients and
  made them available as GEO accession
  [GSE13015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13015).
  GSE13015 data have been parsed into a SummarizedExperiment object
  available in ExperimentHub. This data data could be used as an
  example supporting BloodGen3Module R package.

- [imcdatasets](/packages/imcdatasets) The imcdatasets package
  provides access to publicly available IMC datasets. IMC is a
  technology that enables measurement of > 40 proteins from tissue
  sections. The generated images can be segmented to extract single
  cell data. Datasets typically consist of three elements: a
  SingleCellExperiment object containing single cell data, a
  CytoImageList object containing multichannel images and a
  CytoImageList object containing the cell masks that were used to
  extract the single cell data from the images.

- [LRcellTypeMarkers](/packages/LRcellTypeMarkers) This is an
  external ExperimentData package for LRcell. This data package
  contains the gene enrichment scores calculated from scRNA-seq
  dataset which indicates the gene enrichment of each cell type in
  certain brain region. LRcell package is used to identify specific
  sub-cell types that drives the changes observed in a bulk RNA-seq
  differential gene expression experiment. For more details, please
  visit: https://github.com/marvinquiet/LRcell.

- [MACSdata](/packages/MACSdata) Test datasets from the MACS3 test
  examples are use in the examples of the `MACSr` package. All 9
  datasets are uploaded to the `ExperimentHub`. The original data can
  be found at: https://github.com/macs3-project/MACS/.

- [methylclockData](/packages/methylclockData) Collection of 9
  datasets, andrews and bakulski cord blood, blood gse35069, blood
  gse35069 chen, blood gse35069 complete, combined cord blood, cord
  bloo d gse68456, gervin and lyle cord blood, guintivano dlpfc and
  saliva gse48472". Data downloaded from
  [meffil](https://github.com/perishky/meffil/). Data used to
  estimate cell counts using Extrinsic epigenetic age acceleration
  (EEAA) method Collection of 12 datasets to use with MethylClock
  package to estimate chronological and gestational DNA
  methylationwith estimators to use wit different methylation clocks

- [microbiomeDataSets](/packages/microbiomeDataSets)
  microbiomeDataSets is a collection of microbiome datasets loaded
  from Bioconductor'S ExperimentHub infrastructure. The datasets
  serve as reference for workflows and vignettes published adjacent
  to the microbiome analysis tools on Bioconductor. Additional
  datasets can be added overtime and additions from authors are
  welcome.

- [MouseThymusAgeing](/packages/MouseThymusAgeing) This package
  provides data access to counts matrices and meta-data for
  single-cell RNA sequencing data of thymic epithlial cells across
  mouse ageing using SMARTseq2 and 10X Genommics chemistries. Access
  is provided as a data package via ExperimentHub. It is designed to
  facilitate the re-use of data from Baran-Gale _et al._ in a
  consistent format that includes relevant and informative meta-data.

- [msigdb](/packages/msigdb) This package provides the Molecular
  Signatures Database (MSigDB) in a R accessible objects. Signatures
  are stored in GeneSet class objects form the GSEABase package and
  the entire database is stored in a GeneSetCollection object. These
  data are then hosted on the ExperimentHub. Data used in this
  package was obtained from the MSigDB of the Broad Institute.
  Metadata for each gene set is stored along with the gene set in the
  GeneSet class object.

- [ObMiTi](/packages/ObMiTi) The package provide RNA-seq count for 2
  strains of mus musclus; Wild type and Ob/Ob. Each strain was
  divided into two groups, and each group received either chow diet
  or high fat diet. RNA expression was measured after 12 weeks in 7
  tissues.

- [preciseTADhub](/packages/preciseTADhub) An experimentdata package
  to supplement the preciseTAD package containing pre-trained models
  and the variable importances of each genomic annotation used to
  build the model parsed into list objects and available in
  ExperimentHub. In total, preciseTADhub provides access to n=84
  random forest classification models optimized to predict
  TAD/chromatin loop boundary regions and stored as .RDS files. The
  value, n, comes from the fact that we considered l=2 cell lines
  {GM12878, K562}, g=2 ground truth boundaries {Arrowhead, Peakachu},
  and c=21 autosomal chromosomes {CHR1, CHR2, ..., CHR22} (omitting
  CHR9). Furthermore, each object is itself a two-item list
  containing: (1) the model object, and (2) the variable importances
  for CTCF, RAD21, SMC3, and ZNF143 used to predict boundary regions.
  Each model is trained via a "holdout" strategy, in which data from
  chromosomes {CHR1, CHR2, ..., CHRi-1, CHRi+1, ..., CHR22} were used
  to build the model and the ith chromosome was reserved for testing.
  See https://doi.org/10.1101/2020.09.03.282186 for more detail on
  the model building strategy.

- [ptairData](/packages/ptairData) The package ptairData contains two
  raw datasets from Proton-Transfer-Reaction Time-of-Flight mass
  spectrometer acquisitions (PTR-TOF-MS), in the HDF5 format. One
  from the exhaled air of two volunteer healthy individuals with
  three replicates, and one from the cell culture headspace from two
  mycobacteria species and one control (culture medium only) with two
  replicates. Those datasets are used in the examples and in the
  vignette of the ptairMS package (PTR-TOF-MS data pre-processing).
  There are also used to gererate the ptrSet in the ptairMS data :
  exhaledPtrset and mycobacteriaSet

- [scpdata](/packages/scpdata) The package disseminates mass
  spectrometry (MS)-based single-cell proteomics (SCP) datasets. The
  data were collected from published work and formatted using the
  `scp` data structure. The data sets contain quantitative
  information at spectrum, peptide and/or protein level for single
  cells or minute sample amounts.

- [SimBenchData](/packages/SimBenchData) The SimBenchData package
  contains a total of 35 single-cell RNA-seq datasets covering a wide
  range of data characteristics, including major sequencing
  protocols, multiple tissue types, and both human and mouse sources.

- [SingleMoleculeFootprintingData](/packages/SingleMoleculeFootprintingData)
  This Data package contains data objcets relevanat for the
  SingleMoleculeFootprinting package. More specifically, it contains
  one example of aligned sequencing data (.bam & .bai) necessary to
  run the SingleMoleculeFootprinting vignette. Additionally, we
  provide data that are essential for some functions to work
  correctly such as BaitCapture() and SampleCorrelation().

- [STexampleData](/packages/STexampleData) Collection of spatially
  resolved transcriptomics datasets in SpatialExperiment Bioconductor
  format, for use in examples, demonstrations, tutorials, and other
  purposes. The datasets have been sourced from various publicly
  available sources, and cover several technological platforms.

- [TENxVisiumData](/packages/TENxVisiumData) Collection of Visium
  spatial gene expression datasets by 10X Genomics, formatted into
  objects of class SpatialExperiment. Data cover various organisms
  and tissues, and include: single- and multi-section experiments, as
  well as single sections subjected to both whole transcriptome and
  targeted panel analysis. Datasets may be used for testing of and as
  examples in packages, for tutorials and workflow demonstrations, or
  similar purposes.


New Annotation Packages
=====================

There are 7 new annotation packages in this release of Bioconductor.

- [AHLRBaseDbs](/packages/AHLRBaseDbs) Supplies AnnotationHub with `LRbaseDb`
  Ligand-Receptor annotation databases for many species. All the SQLite files
  are generated by our Snakemake workflow
  [lrbase-workflow](https://github.com/rikenbit/lrbase-workflow). For the
  details, see the README.md of lrbase-workflow.

- [AHMeSHDbs](/packages/AHMeSHDbs) Supplies AnnotationHub with `MeSHDb` NIH MeSH
  annotation databases for many species. All the SQLite files and metadata.csv
  are generated by our Snakemake workflow
  [mesh-workflow](https://github.com/rikenbit/mesh-workflow).

- [AHPathbankDbs](/packages/AHPathbankDbs) The package provides a comprehensive
  mapping table of metabolites and proteins linked to PathBank pathways. The
  tables include HMDB, KEGG, ChEBI, CAS, Drugbank, Uniprot IDs. The tables are
  provided for each of the 10 species ("Homo sapiens", "Escherichia coli", "Mus
  musculus", "Arabidopsis thaliana", "Saccharomyces cerevisiae", "Bos taurus",
  "Caenorhabditis elegans", "Rattus norvegicus", "Drosophila melanogaster", and
  "Pseudomonas aeruginosa"). These table information can be used for Metabolite
  Set (and other) Enrichment Analysis.

- [AHPubMedDbs](/packages/AHPubMedDbs) Supplies AnnotationHub with some
  preprocessed sqlite, tibble, and data.table datasets of PubMed. All the
  datasets are generated by our Snakemake workflow
  [pubmed-workflow](https://github.com/rikenbit/pubmed-workflow). For the
  details, see the README.md of pubmed-workflow.

- [gwascatData](/packages/gwascatData) This package manages a text file in cloud
  with March 30 2021 snapshot of EBI/EMBL GWAS catalog.This simplifies access to
  a snapshot of EBI GWASCAT.  More current images can be obtained using the
  gwascat package.

- [MafH5.gnomAD.v3.1.1.GRCh38](/packages/MafH5.gnomAD.v3.1.1.GRCh38) Store minor
  allele frequency data from the Genome Aggregation Database (gnomAD version
  3.1.1) for the human genome version GRCh38.

- [Orthology.eg.db](/packages/Orthology.eg.db) Orthology mapping package, based
  on data from NCBI, using NCBI Gene IDs and Taxonomy IDs.

New Workflow Packages
=====================

There is 1 new workflow package in this release of Bioconductor. 

- [ExpHunterSuite](/packages/ExpHunterSuite) The ExpHunterSuite
  implements a comprehensive protocol for the analysis of
  transcriptional data using established *R* packages and combining
  their results. It covers all key steps in DEG detection, CEG
  detection and functional analysis for RNA-seq data. It has been
  implemented as an R package containing functions that can be run
  interactively. In addition, it also contains scripts that wrap the
  functions and can be run directly from the command line.


New Books
=====================

There are no new online books.


NEWS from new and existing Software Packages
===================================


[ACE](/packages/ACE)
---

                 Changes in version 1.9.3 (2021-01-28)                  

- Now defaulting to exclude sex chromosomes from model fitting

- Also included sgc argument in twosamplecompare

- Data frame output of ACEcall and twosamplecompare are now restricted
  to selected chromosomes

                 Changes in version 1.9.1 (2021-01-15)                  

- accommodating fitting of chromosomes with only a single germline copy
  (e.g. X and Y in males)

- added the option to specify which cellularities to include in
  squaremodel

- option to save readCounts-object in runACE

[airpart](/packages/airpart)
-------

                       Changes in version 0.0.99                        

- Submitting to Bioconductor...

[AlpsNMR](/packages/AlpsNMR)
-------

                        Changes in version 3.1.5                        

- Removed warning about future_options deprecation

                        Changes in version 3.1.4                        

- bug fix loading bruker files

                 Changes in version 3.1.3 (2020-11-19)                  

- Added instructions to follow a longer tutorial
- nmr_pca_outliers_plot modified to show names in all boundaries of
the plot

                 Changes in version 3.1.2 (2020-11-04)                  

- Bug fix related with Bioconductor Renviron variable
R_CHECK_LENGTH_1_CONDITION

                 Changes in version 3.1.1 (2020-10-30)                  

- Modified order of autor list

                 Changes in version 3.1.0 (2020-10-22)                  

- Package accepted in bioconductor

[ANCOMBC](/packages/ANCOMBC)
-------

                 Changes in version 1.1.5 (2021-03-09)                  

- Bugs fix: fix the bug of inflated p-values and inconsistent output
  formats.

                 Changes in version 1.1.4 (2021-02-28)                  

- Bug fix: fix the bug when metadata contains only a single variable
  and some
  samples were removed with the minimum library size cutoff.

                 Changes in version 1.1.3 (2021-02-19)                  

- Add a warning message for the case of the small number of taxa.

                 Changes in version 1.1.2 (2020-12-08)                  

- Bug fix: fix the bug that the sampling fraction estimate will return
  a
  single number instead of a vector.

                 Changes in version 1.1.1 (2020-11-20)                  

- Integrating with functions from the microbiome package.

[AnnotationDbi](/packages/AnnotationDbi)
-------------

                       Changes in version 1.54.0                        

NEW FEATURES

- There is a new replacement package for the Inparanoid orthology
  packages, called Orthology.eg.db

- This package uses NCBI orthology data to map NCBI Gene IDs between
  species using the usual select() interface

MODIFICATIONS

- UniGene data have been removed from OrgDb and ChipDb packages

- Gene type data have been added to OrgDb and ChipDb packages

                       Changes in version 1.53.0                        

NEW FEATURES

- organismKEGGFrame() provides a data.frame of species names and KEGG
  orgs

MODIFICATIONS

- Conversion from KEGG.db (which was deprecated) to KEGGREST
  - This involved changing how KEGG Ids checked when creating a
  KEGGFrame
  object.

- Removal of any use of Inparanoid data/structures since the
  hom.*.inp.db
  packages have been removed. The data is extremely outdated and there
  is no
  replacement for it.
  - Removed the InparanoidDb object and all of the methods

[AnnotationForge](/packages/AnnotationForge)
---------------

                       Changes in version 1.34.0                        

NEW FEATURES

- Removed UniGene from OrgDb and ChipDb packages

- Added Gene Type table to OrgDb and ChipDb packages

- Added functionality to build Orthology.eg.db package which maps NCBI
  Gene IDs between species

- RSQLite deprecated usage of dbGetQuery for database altering
  statements; updated to use dbExecute instead

[AnnotationHub](/packages/AnnotationHub)
-------------

                       Changes in version 2.99.0                        

MAJOR UPDATES

- (2.99.0) The default caching location has changed. Instead of
  rappdirs::user_cache_dir using tools::R_user_dir. To avoid
  conflicting
  caches, a user will have to manage an old cache location before
  proceeding. Information for handling an old cache location is
  provided in
  the vignette.

- (2.99.0) Another major change, a default caching location is
  automatically
  created in a non interactive session instead of using a temporary
  location. In an interactive session, a user is still prompted for
  permission.

                       Changes in version 2.23.0                        

USER-VISIBLE MODIFICATIONS

- (2.23.2) Create a new all encompassing vignette that references both
  ExpeirmentHub and AnnotationHub. Reference this one vignette in all
  four
  related packages instead of trying to maintain multiple vigenttes
  that were
  essentially the same. This also involves removing
  CreateAnAnnotationPackage

MODIFICATIONS

- (2.23.1) Fixed ERROR message to better indicate vignette
  troubleshooting
  document and fixed reference in Troubleshooting vignette. These
  ERRORs are
  triggered by both AnnotationHub and ExperimentHub so clarified the
  Troubleshooting document is in AnnotationHub.

[AnnotationHubData](/packages/AnnotationHubData)
-----------------

                       Changes in version 1.21.0                        

MODIFICATIONS

- 1.21.9 Add PNG as valid source type

- 1.21.4 Removed vignette for creating annotation hub package.
  Reference and
  refer to single vignette in AnnotationHub

- 1.21.3 Tags for database now combination of biocViews and meta$Tags.
  Also
  checks for valid AnnotationHub or AnnotationHubSoftware biocViews.

- 1.21.2 Add mtx.gz as valid source type

BUG CORRECTION

- 1.21.3 Fixed bug to run make*HubMetadata using "."

INTERNAL BUG CORRECTION

- 1.21.1 misplaced ! clause

REMOVED

- 1.21.5 Removed BioPax. url no longer valid. Resources were old and
  never
  used beyond first addition

[AnVIL](/packages/AnVIL)
-----

                        Changes in version 1.4.0                        

NEW FEATURES

- (v 1.3.1) support Rawls() service (more fine-grained implementation
/ extension of the 'Terra()' orchestration API).

- (v 1.3.2) introduce avworkspace_*() functions for viewing and
updating workflow configurations.

- (v 1.3.3) introduce avnotebooks_() functions for managing notebooks
on workspaces and runtimes.

- (v 1.3.11) introduce avtable_paged() for page-wise access to tables

- (v 1.3.14) introduce avworkspace_clone() for cloning existing
workspaces.

- (v 1.3.21) avworkspaces() returns a tibble of available workspaces.

- (v 1.3.24) gsutil_rsync() supports a regular expresion exclude = to
exclude files from synchronization.

- (v 1.3.24) avworkflow_localize() copies workflow control and / or
output files to the local disk.

USER VISIBLE CHANGES

- (v 1.3.1) service functions have signatures like fun(x, ...,
.__body__ = list(y)), where x is a argument for the 'URL' of the
RESTful interface, and y is an argument for the 'BODY' of POST and
similar requests. The ... provide backward compatibility, and is
used to populate elements of .__body__; the full interface is
required when URL and BODY have identically named arguments.

- (v 1.3.10, 1.3.11) return 'entity' column with name 'table_id',
rather than 'name'.

- (v 1.3.22) localize() / delocalize() warn when dry = TRUE, so that
lack of localization is more apparent.

- (v 1.3.24) gsutil_stat() returns a tibble summaring bucket status,
rather than character().

- (v 1.3.30) Add Referer: header to all Leonardo requests

BUG FIXES

- (v 1.3.6) when .__body__ consists of 1 argument, it is represented
as an unnamed set.

- (v 1.3.7) allow positional matching for .__body__ arguments

- (v. 1.2.1 / 1.3.31) drs_stat() returns a single record per URL when
multiple hashes available.

[AnVILPublish](/packages/AnVILPublish)
------------

                        Changes in version 1.2.0                        

New Features

- Include README.md on Workspace landing page (thanks Vince Carey).

Bug Fixes

- as_workspace() correctly passes an unboxed 'description' attribute
when setting the dashboard.

[APAlyzer](/packages/APAlyzer)
--------

                 Changes in version 1.5.5 (2021-01-31)                  

- Added ThreeMostPairBam to support paired-end bam.

                 Changes in version 1.5.4 (2021-01-10)                  

- Fixed the bug in PASEXP_3UTR.

                 Changes in version 1.5.3 (2021-01-05)                  

- Updated the link of PolyA_DB.

                 Changes in version 1.5.1 (2020-12-07)                  

- Updated Imports and authors.

[artMS](/packages/artMS)
-----

                 Changes in version 1.8.4 (2021-05-12)                  

- Several adjustments in artmsAnalysisQuantifications function:
  - Use NAs instead of 0 when reshaping the data for missing values
  (-log2fc-wide.txt and -log2fc-long.txt files)
  - plotClusteringAnalysis only available when 3 or more comparisons
  available in the configuration file

- R > v4.0.0 is now required

- Update documentation & vignette

- Code cleaning

                 Changes in version 1.8.3 (2021-04-05)                  

- artmsProtein2SiteConversion now supports uniprot id isoforms (thanks
  Emily King)

- Update vignette to make clearer how to provide several protein ids in
  "normalization_reference" (thanks Olga Schubert)

                 Changes in version 1.8.2 (2021-03-18)                  

- artmsProtein2SiteConversion: New PTM available as argument, a new
  user defined `PTM:XXX:yy`. Check documentation to find out more

- QC plots: by default, all qc plots are now output to a folder
  directory, by type (qc-basic, qc-extended, qc_summary)

- artMS working directory: artMS will create all the folders and
  subfolders relative to the working directory. No need to specify the
  full path to the working directory, but the user must set the working
  directory: setwd("/path/to/working/directory/")

- Configuration file data object updates:

- "output" the user can add a folder where would like to have the
  output results file. For example, "output:
  results_202003/example-results.txt" would create the "resutls_202003"
  folder (if it does not exist) with all the results files available
  there

- "LFC" (log2fc) updated to range -0.58 to 0.58, i.e., a fold change
  larger than a 1.5 (instead of 2 as before)

- Update and improve documentation and vignette

- Several bug fixes

                 Changes in version 1.8.1 (2020-10-27)                  

- Update "plotPCA" message

- Update documentation and vignettes

- Code cleaning

[ASpli](/packages/ASpli)
-----

                        Changes in version 2.1.3                        

NEW CITATION 

- ASpli was published on Bioinformatics

- https://doi.org/10.1093/bioinformatics/btab141

FEATURES

- binGenome now displays a progress bar when binning.

BUG FIXES

- gbCounts identifies correctly NA in junction's name

- minAnchor in jCounts was hardcoded so it had no impact on analysis.
  It is being passed correctly now.

                         Changes in version 2.1                         

NEW FUNCTIONS AND FUNCTIONALITIES

- New locus plot helper function .plotGenePattern.

- plots coverage, junctions and architecture

FEATURES

- Enhancement to the vignettes.

- Quick start section was modify in order to use gtf and bam files
  provided by ASpli package

- Add new section ASpli overview in vignette

BUG FIXES

- binGenome function assigns correctly bins located at the start and
  the end of each gene

- binGenome function calculates correctly gene range overlap

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.15.11                       

- Break the limitation of sequence length must have ends less than or
equal to 536870912.

                       Changes in version 1.15.10                       

- fix the issue that idxstatsBam return values with "*"

                       Changes in version 1.15.9                        

- Add rmarkdown as suggest package.

                       Changes in version 1.15.8                        

- update documentation for the case when no BSgenome object is
available.

                       Changes in version 1.15.7                        

- fix the NA values for TSSEscore when infinite value is in the data.

                       Changes in version 1.15.6                        

- fix the missing link of documentation for rtracklyaer:import.

                       Changes in version 1.15.5                        

- remove duplicates when shift reads.

                       Changes in version 1.15.4                        

- Fix the issue when empty object input into exportBamFile.

                       Changes in version 1.15.3                        

- Reuse header when exportBamFile in splitGAlignmentsByCut function.

                       Changes in version 1.15.2                        

- Fix the tag MC in exportBamFile function.

                       Changes in version 1.15.1                        

- write exportBamFile function to replace rtracklayer::export.bam.

[AUCell](/packages/AUCell)
------

                        Changes in version 1.13                         

- New function: plotEmb_rgb() and getSetNames()

[autonomics](/packages/autonomics)
----------

                 Changes in version 0.99.0 (2020-08-05)                 

- Submitted to Bioconductor

[awst](/packages/awst)
----

                       Changes in version 0.99.1                        

NEW FEATURES

- Added awst and gene_filter methods for SummarizedExperiment class.
- Added tests for the new methods.

SIGNIFICANT USER-VISIBLE CHANGES

- Matrix methods of awst and gene_filter now return a gene by sample
matrix.

BUG FIXES

- None.

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- Added a vignette to illustrate the main usage of the package.

SIGNIFICANT USER-VISIBLE CHANGES

- None.

BUG FIXES

- Fixed gene_filter to work with NA values.

[barcodetrackR](/packages/barcodetrackR)
-------------

                 Changes in version 0.99.0 (2021-02-22)                 

- Starting move to bioconductor
- Moving larger extdata to barcodetrackRData github

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.3.4 (2021-04-18)                  

- Add missing import from scran

                 Changes in version 2.3.3 (2021-04-14)                  

- Version bump to trigger new build

                 Changes in version 2.3.2 (2021-04-14)                  

- Bug fixes in handling of divide and conquer inference.

                 Changes in version 2.3.1 (2020-12-14)                  

- scaling of `mu.mu` in `.EmpiricalBayesMu` to match the scale given by
  spike-ins

- scaling of `mu0` in `.BASiCS_MCMC_Start` to match the scale given by
  spike-ins
  when using an EB prior for mu

- lower minimum tolerance `mintol_mu` (1e-5 instead of 1e-3) as a
  default value
  in `.BASiCS_MCMC_ExtraArgs`

[basilisk](/packages/basilisk)
--------

                        Changes in version 1.4.0                        

- Support installation from Python package directories on the
  file system.

- Clean Conda package directories during a system installation to
  reduce disk usage.

[basilisk.utils](/packages/basilisk.utils)
--------------

                        Changes in version 1.4.0                        

- Avoid caching the installer when performing a system
  installation in installConda(). Otherwise, cache in the
  external directory to avoid requiring/polluting BiocFileCache's
  cache.

- Officially give up on Windows 32-bit support in installConda().

- Migrated activateEnvironment() back here, from basilisk.

- Added cleanConda() utility to clean the Conda environment.

- Added setCondaPackageDir() to set the Conda package cache
  directory.

[batchelor](/packages/batchelor)
---------

                        Changes in version 1.8.0                        

- Migrate findMutualNN() to BiocNeighbors.

- Support d=NA in multiBatchPCA() for more convenient disabling
  of PCA in calling functions.

- Bugfix for d=NA with specified subset.row= in fastMNN().

- Added the applyMultiSCE() function to easily apply functions
  across main/alternative Experiments from multiple
  SingleCellExperiment inputs.

- Added the mnnDeltaVariance() function to compute diagnostics
  from the variances of the differences between MNN pairs.

- Added the quickCorrect() function to quickly perform
  intersection, normalization, feature selection and correction.

- Added some clustering-based diagnostics (clusterAbundanceVar(),
  clusterAbundanceTest() and compareMergedClusters()) from the
  OSCA book.

- File-backed matrices are now realized into memory prior to
  multiBatchPCA().

[BayesSpace](/packages/BayesSpace)
----------

                        Changes in version 1.1.3                        

Minor improvements and fixes

- Added information to documentation.
- getRDS() updated with new URL.

                        Changes in version 1.1.2                        

Minor improvements and fixes

- clusterPlot() accepts character vectors and factors as arguments to
label.

                        Changes in version 1.1.1                        

Minor improvements and fixes

- spatialPreprocess() uses exact rather than approximate PCA by
default.

                        Changes in version 1.1.0                        

New Bioconductor devel (3.13)

- Version numbering change with Bioconductor version bump

[beachmat](/packages/beachmat)
--------

                        Changes in version 2.8.0                        

- Improve the efficiency of sparse row subsetting in
  non-DelayedArray rowBlockApply().

- Avoid overhead of DelayedArray block processing when
  DelayedArray is pristine and the type is supported.

- Migrated whichNonZero() from scuttle.

- Added toCsparse() to make it easier to convert SparseArraySeeds
  to CsparseMatrixes.

- Added realizeFileBackedMatrix() to, well, realize a
  DelayedMatrix with file-backed components.

[BEclear](/packages/BEclear)
-------

                 Changes in version 2.7.1 (2020-11-03)                  

- fixed some warnings from R Check and BiocCheck

[BgeeCall](/packages/BgeeCall)
--------

                        Changes in version 1.7.3                        

- Add approaches based on pValue and qValue to generate
  present/absent calls

- Default approach to generate calls is pValue

- Add function ```merging_libraries``` allowing to merge calls
  per condition

[BiocCheck](/packages/BiocCheck)
---------

                        Changes in version 1.27                         

BUG FIX

- (1.27.17) Update support site watched tags. tags are case insensitive

- (1.27.15) Reporting checking of vignette despite package type
  (@lshep, #136)

- (1.27.9) Allow portability of child Rmd documents via parseFile

- (1.27.3) Correct check for if package already exists in CRAN/Bioc

- (1.27.3) Correct check for single colon use

- (1.27.2) Correct path to R license database file by calling
  R.home('share').

NEW FEATURES

- (1.27.16) Check vignettes for all package types (@lshep, #136)

- (1.27.12) Check for `LazyData: TRUE` in the DESCRIPTION (@lshep,
  #128)

- (1.27.11) R version dependency check in the 'DESCRIPTION' is now a
  'NOTE' (@lshep, #126)

- (1.27.10) Check for 'error' and other keywords in signaler functions,
  'message', 'warning', and 'stop' etc. (@hpages, #125)

- (1.27.8) Check for 'tests' entry in '.Rbuildignore'

- (1.27.7) Removed BiocCheck and BiocCheckGitClone installation
  scripts;
  recommended usage is `BiocCheck()`

- (1.27.6) Check that a user has the package name in watched tags of
  support site

- (1.27.5) Check for 'paste' / 'paste0' in signaler functions,
  'message',
  'warning', and 'stop' (@LiNk-NY, #64)

- (1.27.4) Check for downloads from external resources (github, gitlab,
  bitbucket, dropbox; @LiNk-NY, #75)

- (1.27.1) Check that licenses do not exclude classes of users,
  e.g., non-academic users.

[BiocFileCache](/packages/BiocFileCache)
-------------

                        Changes in version 1.99                         

MAJOR UPDATES

- (1.99.0) The default caching location has changed. Instead of
  rappdirs::user_cache_dir using tools::R_user_dir. To avoid
  conflicting
  caches, a user will have to manage an old cache location before
  proceeding. Information for handling an old cache location is
  provided in
  the vignette.

- (1.99.0) Another major change, a default caching location is
  automatically
  created in a non interactive session instead of using a temporary
  location. In an interactive session, a user is still prompted for
  permission.

- (1.99.0) An enviornment variable may be set system wide or user wide
  to
  control the default caching location: BFC_CACHE. Note: do not use R
  variables or command line export to set this variable. It must be set
  system
  wide or user wide for reproducibility in future R sessions or else it
  must
  be specified upon ever usage. It must be set before calling
  library(BiocFileCache) to take effect.

- (1.99.0) Fixes partial argument matching error in SQL function
  SQLExecute

- (1.15.1) Added file locking for thread-safe SQL operations. Thanks
  for the PR @LTLA

BUG FIX

- (1.99.6) cleanbfc() incorrect format string; see
  https://github.com/Bioconductor/BiocFileCache/issues/31

[BiocNeighbors](/packages/BiocNeighbors)
-------------

                       Changes in version 1.10.0                        

- Migrated findMutualNN() from batchelor.

- Vendored the RcppAnnoy headers for greater reproducibility.

- Added a distance="Cosine" option for all algorithms.


[BiocParallel](/packages/BiocParallel)
------------

                        Changes in version 1.26                         

USER VISIBLE CHANGES

- (v 1.25.2) bpvalidate() gains an argument to control warning /
  error / silent signaling, and returns a 'BPValidate' object.

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.10.0                        

NEW FEATURES

- `biocPkgRanges` allows for easy identification of package statuses
  from
  the build report for a specified range of packages (ordered
  alphabetically)

- `biocBuildEmail` provides core-team functionality for sending email
  notifications to package maintainers

SIGNIFICANT USER-VISIBLE CHANGES

- `biocBuildEmail` allows for saving a credentials file for email
  authentication via the `credFile` argument

- `setCache` uses `tools::R_user_dir("BiocPkgTools", "cache")` instead
  of
  `rappdirs::user_cache_dir`

BUG FIXES

- `biocBuildReport` accounts for some packages whose `DESCRIPTION` file
  is
  malformed

- `biocBuildReport` updated to changes in the build report format

[BiocStyle](/packages/BiocStyle)
---------

                       Changes in version 2.20.0                        

- Addressed styling issues for output code blocks, introduced by
  changes
  in rmarkdown version 2.7.0
  (https://github.com/Bioconductor/BiocStyle/issues/86)

[biocthis](/packages/biocthis)
--------

                       Changes in version 1.1.10                        

NEW FEATURES

- Added the function use_bioc_coc() as requested by Luke Zappia et
al.

                        Changes in version 1.1.9                        

NEW FEATURES

- Now use_bioc_github_action() has a docker argument which controls
whether to build a docker image at the end of the GHA workflow (only
on Linux) as requested by Kévin Rue-Albrecht.

                        Changes in version 1.1.7                        

BUG FIXES

- Switch to match usethis 2.0.1 which changed a lot of the internal
code in biocthis.

                        Changes in version 1.1.4                        

NEW FEATURES

- Switched from knitcitations to RefManageR given the discussion at
https://github.com/cboettig/knitcitations/issues/107.

                        Changes in version 1.1.3                        

BUG FIXES

- Updated to use usethis version 2.0.0 or newer. Check the following
for more information on usethis changes:
https://twitter.com/JennyBryan/status/1337858543404285952?s=20.

[biocViews](/packages/biocViews)
---------

                       Changes in version 1.59.0                        

ENHANCEMENT

- (1.57.3) Add biocViews term DifferentialDNA3DStructure

- (1.57.2) Add CRAN packages to reverse dependency list

- (1.57.1) Add biocViews term Chlamydomonas_reinhardtii

[biodb](/packages/biodb)
-----

                Changes in version 0.99.11 (2021-05-17)                 

- Change example inside BiodbConfig to avoid misinterpretation of
  `set('cache.directory', '~/my.biodb.cache')`, leading to believe that
  some files are written inside USER HOME folder when running the
  example.

                Changes in version 0.99.10 (2021-05-07)                 

- Correct documentation of C++ function.

                 Changes in version 0.99.9 (2021-05-07)                 

- Solving some NOTES from BiocCheck.

- Correct example in Progress class.

                 Changes in version 0.99.8 (2021-05-06)                 

- Correct template travis.yml for extensions: missing deps install, run
  all checks.

- Improve template Makefile for extensions.

- Add missing test-cpp.R template file for running C++ tests from
  testthat.

- Limit by default the entries to test to one entry inside generic
  tests.

- Improve vignettes.

- Rename default vignette into "biodb.Rmd".

                 Changes in version 0.99.7 (2021-04-27)                 

- Renamed Biodb class into BiodbMain in order to avoid "Rd warning:
  Previous alias or file overwritten by alias: biodb" on Windows
  platform.

- Implement newInst() global function for creating new BiodbMain
  instance.

                 Changes in version 0.99.6 (2021-04-27)                 

- Rebuilding doc.

                 Changes in version 0.99.5 (2021-04-27)                 

- Add missing parameters documentation for runGenericTests().

                 Changes in version 0.99.4 (2021-04-27)                 

- Move long tests to separate directory "long".

                 Changes in version 0.99.3 (2021-04-27)                 

- Added "biodb" as watched tag on my profile on support site
  https://support.bioconductor.org/.

                 Changes in version 0.99.2 (2021-04-27)                 

- Switch to MassBank extract for testing MassCsvFile and MassSqlite
  connectors.

                 Changes in version 0.99.1 (2021-04-26)                 

- Remove Git files refused by BiocCheckGitClone.

- Use CHECK_RENVIRON in local tests.

- Correct condition in BiodbEntryFields::getRealName() that did not
  pass check.

- Add all doc files man/*.Rd for BiocCheck run on
  http://bioconductor.org.

                 Changes in version 0.99.0 (2021-04-22)                 

- Submitted to Bioconductor

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.48.0                        

NEW FEATURES

- getSequence() now allows the cache to be turned off via the
  'useCache'
  argument.

- Automatic detection of SSL issues with Ensembl, and appropriate
  settings
  applied to httr functions used by biomaRt.

BUG FIXES

- Addressed issue with getSequence() and ID types that are not
  available on the 'sequences' page.  This could result in truncated
  sequences being returned from a query.

- getBM() would fail if it found a cache entry, but the file was
  corrupted.
  Invalid entries are now detected and deleted if encountered.

[BioNERO](/packages/BioNERO)
-------

                 Changes in version 0.99.0 (2021-03-05)                 

- Submitted package to Bioconductor.

[biosigner](/packages/biosigner)
---------

                       Changes in version 1.19.2                        

MINOR MODIFICATION

- getMset: minor correction in documentation

[BloodGen3Module](/packages/BloodGen3Module)
---------------

                       Changes in version 0.99.38                       

- submission to Bioconductor

[bluster](/packages/bluster)
-------

                        Changes in version 1.2.0                        

- 
  Previously zero-weight edges are now assigned a nominal
  positive weight in makeSNNGraph().

- 
  Added MbkmeansParam() to wrap mini-batch k-means from mbkmeans.

- 
  Added SOMParam() to wrap self-organizing map implementation
  from kohonen.

- 
  Added AffinityParam() to wrap the affinity propagation code
  from apcluster.

- 
  Added DbscanParam() to provide a custom DBSCAN implementation
  with automatic eps choice.

- 
  Added PamParam() to wrap the PAM implementation from cluster.

- 
  Added ClaraParam() to wrap the CLARA implementation from
  cluster.

- 
  Added AgnesParam() to wrap the agglomerative nesting method
  from cluster.

- 
  Added DianaParam() to wrap the divisive analysis method from
  cluster.

- 
  Added clusterSweep() to easily perform parameter sweeps via
  clusterRows().

- 
  Added linkClusters() to find relationships between clusters in
  different clusterings.

- 
  Added compareClusterings() to compute similarities between
  multiple clusterings.

- 
  Added nestedClusters() to quantify the degree of nesting across
  two clusterings.

- Moved objects into objects$kmeans for KmeansParam() when
  full=TRUE.

- Moved objects into objects$hclust for HclustParam() when
  full=TRUE.

- Added clusterRMSD() to compute the root-mean-squared-deviation
  for each cluster.

[BridgeDbR](/packages/BridgeDbR)
---------

                        Changes in version 2.1.2                        

BUG FIXES

- Fixed the link to the webpage to download mapping files

[BUSpaRse](/packages/BUSpaRse)
--------

                 Changes in version 1.5.2 (2020-12-04)                  

- Updated code that queries Ensembl for new version of biomaRt.
- Make sure there're no NAs in tr2g_intron when chrs_only = TRUE.
- Bypassed issue in BSgenome::getSeq in get_velocity_files when
genome
is DNAStringSet.


[CAGEr](/packages/CAGEr)
-----

                       Changes in version 1.34.0                        

BUG FIXES

- Reform the CTSS class.  New accessor: `CTSS()` (with no dot).

- Correct a class error when loading BAM files. (Closes #36).

- Use the BSgenome object from the main environment if available.

[CAMERA](/packages/CAMERA)
------

                       Changes in version 1.47.1                        

BUG FIXES

- Thanks to Lain (INRAe) and Team W4M for this PR, fixing argument
  passing
  to annotatedDiffreport()

[CARNIVAL](/packages/CARNIVAL)
--------

                        Changes in version 2.1.0                        

- Changed API: introduced separate functions for different flavours
of
CARNIVAL.
- A more convenient way to work with CARNIVAL parameters. Parameters
setup is possible through jsons and function calls.
- Better file naming to prevent concurrent file writing when running
several instances of CARNIVAL.
- An easy way to add more solvers.
- A possibility to tune CARNIVAL setups for cplex: manual addition of
parameters is possible through jsons.
- General improvements and refactoring of the code, removed multiple
duplications.
- Removed multiple experimental conditions in a matrix form (was not
used).
- Inputs are transformed to vectors (except prior knowledge network).
- MetaInfo is saved: runId, parsed data (internal data
representation), start time, CARNIVAL flavour.
- Reading from preparsed data and/or lp file is possible.

[cbaf](/packages/cbaf)
----

                 Changes in version 1.14.0 (2021-04-28)                 

New Features

- Heatmaps can now be stored as PDF files.

- Functions now show shorter and more clearer messages.

- Package no longer requires Java Runtime Environment for storing excel
  files
  and it is compatible with older 32 bit operating systems.

- log z-scores provided by cgdsr are used instead of z-scores.

[cBioPortalData](/packages/cBioPortalData)
--------------

                        Changes in version 2.4.0                        

New features

- Vignettes now include additional information (#38, @lwaldron)
- getDataByGenePanel deprecated for getDataByGenes which handles
input
of both gene panels and genes
- cBioPortalData now allows for gene inputs as either Entrez IDs or
Hugo symbols (#24, @jucor) and sampleIds input
- When gene inputs are provided, the by argument has to agree with
the
type of genes provided (either be entrezGeneId or hugoGeneSymbol).

Bug fixes and minor improvements

- Fixed an issue where the labels in the metadata from cBioDataPack
were missing ('LICENSE' and 'Fusion'; #37)
- loadStudy allows cleanup=TRUE for removing files after untar-ing
- Published article now available with citation("cBioPortalData")

[cbpManager](/packages/cbpManager)
----------

                        Changes in version 0.1.1                        

New features

- 'Validation' tab allows to validate created files.
- Improved usability by changing descriptions and adding interactive
tours

                        Changes in version 0.1.0                        

New features

- much of the functionality available, in a proof of concept format.

                        Changes in version 0.0.1                        

New features

- backbone of the project started!

[celda](/packages/celda)
-----

                 Changes in version 1.7.7 (2021-04-12)                  

- Added handling for sparse matrices

                 Changes in version 1.7.6 (2021-04-04)                  

- Added functions for creating HTML reports

- Fixed bug in decontX plotting

                 Changes in version 1.7.4 (2021-03-09)                  

- Enable input of raw/droplet matrix into decontX to estimate ambient
  RNA

[CelliD](/packages/CelliD)
------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[cellmigRation](/packages/cellmigRation)
-------------

                Changes in version 0.99.11 (2021-05-11)                 

- Addressed all points from the Bioconductor review
- Bug fixes and documentation update for Bioconductor release

                 Changes in version 0.99.0 (2020-09-02)                 

- Submitted to Bioconductor

[chipenrich](/packages/chipenrich)
----------

                       Changes in version 2.16.0                        

- Transition to Kai Wang as maintainer.

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.25.6                        

- Fix a bug in estLibSize introduced by last push.

                       Changes in version 3.25.5                        

- Add LazyDataCompression in description

                       Changes in version 3.25.4                        

- Add choice endMinusStart to annotatePeakInBatch.

                       Changes in version 3.25.3                        

- fix the missing link of documentation for rtracklyaer:import.

                       Changes in version 3.25.2                        

- update documentation.

- update findEnhancers to for known interaction data

                       Changes in version 3.25.1                        

- fix the bug for genomicElementDistribution when the peak length is
  zero.

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.27.4                        

- bug fixed in determine downstream gene (2021-04-27, Thu)
- https://github.com/YuLab-SMU/ChIPseeker/pull/148
- getBioRegion now supports '3UTR' and '5UTR' (2021-03-30, Tue)
- https://github.com/YuLab-SMU/ChIPseeker/pull/146

                       Changes in version 1.27.3                        

- add two parameter, cex and radius, to plotAnnoPie (2021-03-12, Fri)
- https://github.com/YuLab-SMU/ChIPseeker/pull/144

                       Changes in version 1.27.2                        

- bug fixed of getGenomicAnnotation (2021-03-03, Wed)
- https://github.com/YuLab-SMU/ChIPseeker/issues/142

                       Changes in version 1.27.1                        

- Add support for EnsDb annotation databases in annotatePeak.
- https://github.com/YuLab-SMU/ChIPseeker/pull/120

[ChromSCape](/packages/ChromSCape)
----------

                        Changes in version 1.1.3                        

Major Changes

- Support "multi-feature" analysis, e.g. parallel analysis of
multiple
features (bins, peaks or gene) on the same object.

- New "Coverage" tab & functions generate_coverage_tracks() and
plot_coverage_BigWig() to generate cluster coverage tracks and
interactively visualise loci/genes of interest in the
application.

- New inter- and intra-correlation violin plots to vizualise cell
correlation distribution between and within clusters.

- New normalization method : TF-IDF combined with systematic
removal of PC1 strongly correlated with library size.

- Simple 'Copy Number Alteration' approximation & visualization
using 'calculate_CNA' function for genetically re-arranged
samples, provided one or more control samples.

- New generate_analysis() & generate_report() functions to run a
full-on ChromSCape analysis and/or generate an HTML interactive
report of an existing analysis.

- Supports 'custom' differential analysis to find differential
loci between a subset of samples and/or clusters.

- New pathway overlay on UMAP to visualize cumulative pathways
signal directly on cells.

- Now supports 'Fragment Files' input (e.g. from 10X cell ranger
scATAC pipeline), using a wrapper around 'Signac' package
FeatureMatrix() function.

- New 'Contribution to PCA' plots showing most contributing
features and chromosome to PCA.

- Restructuration of the ChromSCape directory & faster
reading/saving of S4 objects using package 'qs'.

Minor Changes

- RAM optimisation & faster pearson cell-to-cell correlations with
'coop' package, and use of 'Rcpp' for as_dist() RAM-efficient
distance calculation.

- Faster correlation filtering using multi-parallel processing.

- plot_reduced_dim now supports gene input to color cells by gene
signal.

- All plots can now be saved in High Quality PDF files.

- Changed 'geneTSS' to 'genebody' with promoter extension to
better reflect the fact that mark spread in genebodies.

- Possibility to rename samples in the application.

- Downsampling of UMAPs & Heatmaps for fluider navigation.

- Changed 'total cell percent based' feature selection to manual
selection of top-covered features, as the previous was srongly
dependent on the experiment size.

- Faster sparse SVD calculation.

- Faster differential analysis using pairWise Wilcoxon rank test
from 'scran' package.

[CIMICE](/packages/CIMICE)
------

                       Changes in version 0.99.0                        

Overview:

- First commit.

New functionalities:

- Input dataset read an creation

- CIMICE analysis and CPMC inference

- Output data visualization

[circRNAprofiler](/packages/circRNAprofiler)
---------------

                        Changes in version 1.5.3                        

- Removed citr package from DESCRIPTION

[cleanUpdTSeq](/packages/cleanUpdTSeq)
------------

                       Changes in version 1.29.1                        

- rewrite the package by Haibo

[clusterProfiler](/packages/clusterProfiler)
---------------

                       Changes in version 3.99.1                        

- Add new data set, DE_GSE8057, which contains DE genes obtained from
GSE8057 (2020-03-08, Mon)

                       Changes in version 3.99.0                        

- Add KEGG enrichment analysis of Human Gut Microbiome data
(2021-02-20, Sat)

                       Changes in version 3.19.1                        

- setting default timeout to 300 for downloads (2021-02-05, Fri)
- fixed download method setting
- capable of setting KEGG download method via
options(clusterProfiler.download.method = METHOD) (2020-12-31, Thu)

[clustifyr](/packages/clustifyr)
---------

                 Changes in version 1.3.3 (2021-02-28)                  

- Launch shiny app with `run_clustifyr_app()`

- Plot and GO for most divergent ranks in correlation of query vs
  reference

                 Changes in version 1.3.2 (2021-02-25)                  

- `build_atlas()` for combining references

- More Q&A

                 Changes in version 1.3.1 (2020-12-26)                  

- Q&A section

- Now defaults to top 1000 variable genes in Seurat (including v4)

- Bug fixes

[CNVfilteR](/packages/CNVfilteR)
---------

                        Changes in version 1.5.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- loadCNVcalls() does not check cnvs.file names by default when loading
  (read.csv()) the cnvs.file

- loadCNVcalls() allows optional check.names.cnvs.file parameter

- Vignette updated

MINOR

- Added rmarkdown to Suggets in DESCRIPTION file

                        Changes in version 1.5.1                        

BUG FIXES

- Bug fixed: SNVs were not being correctly loaded after last
  Bioconductor update

SIGNIFICANT USER-VISIBLE CHANGES

- plotVariantsForCNV() allows two new parameters for customize legend
  visualization

- plotAllCNVs() allows 'genome' parameter to work with different
  genomes

MINOR

- Minor vignette fixes

- Other minor fixes

[CNViz](/packages/CNViz)
-----

                 Changes in version 0.99.3 (2021-04-14)                 

- Added a NEWS.md file to track changes to the package
- Submitted to Bioconductor

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

[coMET](/packages/coMET)
-----

                 Changes in version 1.23.1 (2021-05-16)                 

- Update datasets with Biomart
  For example:
  data(allIG)
  allIG@biomart@httr_config <- list()
  save(allIG,file="XXXX")

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

                       Changes in version 2.7.10                        

- `anno_simple()`: text symbols can have nchar > 1.

- `anno_text()`: add `show_name` argument.

                        Changes in version 2.7.9                        

- add `frequencyHeatmap()`.

- add `Heatmap3D()`.

                        Changes in version 2.7.8                        

- add `cluster_between_groups()`.

- add `graphics` argument in `anno_block()`.

                        Changes in version 2.7.7                        

- discrete numeric legend labels are in correct order now.

- parallel is implemented with foreach + doParallel

- expression is properly processed for discrete legends

- `adjust_dend_by_x()`: simplified the representation of units.

- number of split can be the same as number of matrix rows/columns.

                        Changes in version 2.7.6                        

- `Legend()`: add a new argument `grob`.

                        Changes in version 2.7.5                        

- `anno_block()`: add `labels_offset` and `labels_just`.

- `anno_lines()`: `show_points` can be a vector.

- `pheatmap()`: support `kmeans_k`.

                        Changes in version 2.7.4                        

- add `save_last` option in `ht_opt()`.

                        Changes in version 2.7.1                        

- `normalize_comb_mat()`: add `full_comb_sets` and `complement_set`
  arguments to control
  full sets of combination sets.

- adjust the space of column title according to ggplot2.

- `Legend()`: for title_position == "lefttop", the title position is
  adjusted.

- Legends are automatically adjusted according to the device size when
  resizing the device.

- `Legend()`: add `interval_dist` to control the distance of two
  neighbouring breaks.

- Fixed a bug that it crashes Rstudio

- `make_comb_mat()`: print warning messages when there are NA values in
  the matrix.

- temporary solution for woking under retina display with Rstudio

- add `bin_genome()` and `normalize_genomic_signals_to_bins()`

- print messages if directly sending `anno_*()` functions to
  `top_annotation` or similar arguments.

- `pheatmap()`: set heatmap name to " " so that there is no legend
  title by default.

- also translate `stats::heatmap()` and `gplots::heatmap.2()`.

- move all code for interactive heatmap to InteractiveComplexHeatmap
  package.

[ComPrAn](/packages/ComPrAn)
-------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[conclus](/packages/conclus)
-------

                Changes in version 0.99.343 (2021-04-09)                

- Submission to Bioconductor

- NAMESPACE

- Imported BiocFileCache package and R_user_dir() of tools package
- Exported new function conclusCacheClear()

- DESCRIPTION

- Removed LazyData: true.

- DataFormatting.R

- Added a caching system for retrieveFromGEO()
- Created conclusCacheClear() to delete the cache
- Updated documentation

- loadDataset.R

- Updated documentation
- Simplified nested "if" in loadDataset.R

- methods-normalization.R

- Modified the use of getBM to retrieve only genes of the count
matrix (instead of the all database)

- test_setters.R/test_getters.R

- Modified documentation for loadDataOrMatrix()
- Simplified nested "if" in loadDataset.R

- test_loadData.R

- Adapted the unit tests to the new format of coldata and rowdata

- test_scRNAseq-methods.R

- Changed the experiment name to "Light_Experience"
- Used tempdir() for output directory

- inst

- Added inst/script to generate data on inst/extdata
- New data generated

- vignette

- Used tempdir() for output directory
- Specified other parameters in the first example of runCONCLUS
for the very small dataset
- Replaced old paths by new ones

[condiments](/packages/condiments)
----------

                 Changes in version 0.99.0 (2021-03-01)                 

- Submitted to Bioconductor

[CONSTANd](/packages/CONSTANd)
--------

                       Changes in version 0.99.6                        

- Update vignette: add MA plotting and DEA sections, and example of bad
  MA plot

                       Changes in version 0.99.5                        

- depends: R4.1

- updated BugReports link

                       Changes in version 0.99.4                        

- negative input values are not allowed

- warn and automatically replace zero input values by NA

- bugfix in warning message formatting

                       Changes in version 0.99.3                        

- bugfix in float specification of warning message

                       Changes in version 0.99.2                        

- bugfix remove stray Rproj file

                       Changes in version 0.99.1                        

- bugfix in documentation examples

                       Changes in version 0.99.0                        

- initial release

[cosmosR](/packages/cosmosR)
-------

                 Changes in version 0.99.2                  

- Submitted to bioRxiv
- Release of github page
- Submitted to Bioconductor

[csaw](/packages/csaw)
----

                       Changes in version 1.26.0                        

- filterWindowsGlobal() finally behaves correctly for
  variable-width data=.

- normFactors() and normOffsets() accept DGEList inputs in their
  object= and se.out= arguments.

- mergeResults() and friends now default to taking tab= from the
  mcols() of the inputs.

[cTRAP](/packages/cTRAP)
-----

                       Changes in version 1.10.0                        

Improvements to graphical interface functions:

- New launchDrugSetEnrichmentAnalysis() function to analyse drug set
enrichment and visualize respective results
- launchCMapDataLoader():
- Now allows to load multiple CMap perturbation types
simultaneously
- Keep selected timepoint, dosage and cell line options when
selecting another perturbation type
- Add bubble plot of CMap perturbation types
- launchResultPlotter():
- Now allows to view tables below specific plots and
drag-and-select those plots to filter data in those same tables
- When plotting targeting drugs and similar perturbations, update
available columns and correctly use user-selected column to plot
- launchMetadataViewer() now correctly parses values from Input
attributes as numeric

Major changes

- prepareCMapPerturbations(): directly set perturbation type, cell
line, timepoint and dosage conditions as arguments
- rankSimilarPerturbations() and predictTargetingDrugs():
- Avoid redundant loading of data chunks, slightly decreasing run
time
- Lower memory footprint when using NCI60's gene expression and
drug sensitivity association (now available in HDF5 files) by
loading and processing data in chunks
- Faster GSEA-based score calculation (up to 4-7 times faster)
- New threads argument allows to set number of parallel threads
(not supported on Windows)
- New chunkGiB argument allows to set size of data chunks when
reading from supported HDF5 files (decreases peak RAM usage)
- New verbose argument allows to increase details printed in the
console
- prepareDrugSets(): allow greater control on the creation of bins
based on numeric columns, including the setting of maximum number of
bins per column and minimum bin size
- analyseDrugSetEnrichment() and plotDrugSetEnrichment(): allow to
select columns to use when comparing compound identifiers between
datasets

Bug fixes and minor changes

- filterCMapMetadata(): allow filtering CMap metadata based on
multiple perturbation types
- prepareDrugSets(): fix issues with 3D descriptors containing
missing
values
- plot():
- Fix wrong labels when plotting targetingDrugs objects
- Avoid printing "NA" in labels identifying metadata for
perturbations
- plotTargetingDrugsVSsimilarPerturbations():
- Fix highlighting of plot points depending whether drug activity
is directly proportional to drug sensitivity
- Include rug plot
- When subsetting a perturbationChanges or an
expressionDrugSensitivityAssociation object, passing only one
argument extracts its columns as in previous versions of cTRAP
(similarly to when subsetting a data.frame)
- analyseDrugSetEnrichment(): for the resulting table, the name of
the
first column was renamed from pathway to descriptor

[customCMPdb](/packages/customCMPdb)
-----------

                 Changes in version 1.1.0 (2020-10-02)                  

- Initial version


[CytoGLMM](/packages/CytoGLMM)
--------

                 Changes in version 0.99.0 (2021-02-19)                 

- Submitted to Bioconductor

[cytomapper](/packages/cytomapper)
----------

                 Changes in version 1.3.6 (2021-04-23)                  

- scaleImages accepts numeric vector value

                 Changes in version 1.3.5 (2021-04-07)                  

- Added measureObjects function

                 Changes in version 1.3.4 (2021-03-20)                  

- Support on disk representation of images

                 Changes in version 1.3.3 (2021-01-24)                  

- Added snapshot tests for shiny

- Support win32 again

                 Changes in version 1.3.2 (2021-01-12)                  

- Updated citation

                 Changes in version 1.3.1 (2020-12-01)                  

- Allow thick border contours

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

                       Changes in version 1.29.1                        

- add citation.

[DAMEfinder](/packages/DAMEfinder)
----------

                        Changes in version 1.3.1                        

- change cmdscale.out for eigen.vectors in methyl_MDS_plot

                        Changes in version 1.3.0                        

- update with new Bioc version

[dce](/packages/dce)
---

                 Changes in version 0.99.0 (2021-01-25)                 

- Submitted to Bioconductor

[decompTumor2Sig](/packages/decompTumor2Sig)
---------------

                 Changes in version 2.7.3 (2021-04-01)                  

- Fix: changed example for adjustSignaturesForRegionSet() to limit
  memory usage
  (previous example produced error during check on Windows for arch
  'i386').

- Fix: made sure that the extension of genomic regions by half the
  sequence
  pattern (needed for the adjustSignaturesForRegionSet() function) does
  not
  result in out-of-bounds regions.

                 Changes in version 2.7.2 (2021-03-26)                  

- Updated readAlexandrovSignatures() to read the COMSIC signature
  format v3.2
  (published in March 2021).

                 Changes in version 2.7.1 (2021-03-21)                  

- Updated readAlexandrovSignatures() to add the possibility to read
  COSMIC
  signatures of version 3.1 directly from an Excel file (as provided on
  the
  COSMIC website).

- Added possibility to adjust/normalize mutational signatures to
  specific
  subsets of the human genome (defined by means of GRanges objects).
  The
  adjustment/normalization is performed accoring to the nucleotide
  frequencies
  in the specified regions (with respect to nucleotide frequncies in
  the reference sequences, e.g., the whole genome, for which the
  signatures
  were derived in the first place).

[decoupleR](/packages/decoupleR)
---------

                       Changes in version 0.99.0                        

New features

- All new features allow for tidy selection. Making it easier to
evaluate different types of data for the same method. For instance, you can
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

[DeepPINCS](/packages/DeepPINCS)
---------

                 Changes in version 0.99.0 (2021-03-21)                 

- submission to Bioconductor

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

                       Changes in version 1.27.1                        

- Fix: Export n() from
  dplyrs.[#27](https://github.com/lpantano/DEGreport/issues/27#issuecomment-819201968)

[DelayedMatrixStats](/packages/DelayedMatrixStats)
------------------

                       Changes in version 1.14.0                        

- Fix for missing na.rm= argument in *AvgsPer*Set functions.

- DelayedMatrixStats no longer has a hard requirement on
  HDF5Array or BiocParallel.

- Correct handling of drop= by quantile functions (<URL:
  https://github.com/PeteHaitch/DelayedMatrixStats/pull/71>).

- Fix 2 issues with how the center argument is handled (<URL:
  https://github.com/PeteHaitch/DelayedMatrixStats/pull/69>).

[DelayedRandomArray](/packages/DelayedRandomArray)
------------------

                        Changes in version 1.0.0                        

- 
  New package DelayedRandomArray, for delayed generation of
  random numbers.

[densvis](/packages/densvis)
-------

                        Changes in version 1.2.0                        

- Use umap-learn python library instead of densmap-learn. Add umap
function for non density-preserving umap.
- Add normalize argument to dens-SNE.

[DepecheR](/packages/DepecheR)
--------

                 Changes in version 1.7.2 (2021-03-25)                  

- Major internal changes to the depeche function, with two user
  consequences:
  o The dualDepeche option is deprecated, as it made the function very
  heavy
  to maintain and was not flexible enough to be of great use.
  o The interface to dAllocate is much improved, allowing for smooth
  allocation of new data to an established model, which makes large
  dualDepeche runs, constructed outside of the function, possible in a
  more
  versatile way than previously.

                 Changes in version 1.7.1 (2021-02-02)                  

- Adding the option of not scaling the data within the depeche function

- Condensing the code for the depeche scaling procedure

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.31.16                       

- Turning off outlier replacement with glmGamPoi fitting.

                       Changes in version 1.31.15                       

- Added 'saveCols' in results() and lfcShrink() to pass
  metadata columsn to output.

                       Changes in version 1.31.13                       

- Allow additional arguments to be passed to data-accessing
  functions in integrateWithSingleCell().

                       Changes in version 1.31.2                        

- Fixed interface with glmGamPoi so that normalizationFactors
  can be used. Thanks to Michael Schubert for spotting this
  and to Constantin Ahlmann-Eltze for pointing out the fix.

[DEWSeq](/packages/DEWSeq)
------

                 Changes in version 1.5.2 (2021-05-05)                  

- Update vignette for dispersion estimation

                 Changes in version 1.4.2 (2020-09-30)                  

- Fix figure scaling issue in vignette

[DExMA](/packages/DExMA)
-----

                       Changes in version 0.99.0                        

- DExMA release.

[DIAlignR](/packages/DIAlignR)
--------

                       Changes in version 1.3.18                        

- Using data.table instead of data.frame for modify-in-place.

- Added support for Metabolomics DIA data.

- Hierarchical clustering based alignment.

- Create a child run (features + chromatograms) from two parents.

- Added support for sqMass files.

                        Changes in version 1.3.5                        

- Supporting transition level intensity for SAINTq

- Added support for pyopenms.

- Added support for sqMass files.

- Added parallelization using BiocParallel.

- Using context-specific qvalues to determine reference.

- Alignment is done over multipeptide instead of multiprecursor.

- Precursors of a peptides are forced to have same feature-RT.

- Savitzky Golay smoothing in C++.

- Fast and light global alignment functions.

[DiffBind](/packages/DiffBind)
--------

                         Changes in version 3.2                         

- New type of plot: dba.plotProfile()

- Can mix single-end and paired-end bam files

- Various bug fixes

[diffuStats](/packages/diffuStats)
----------

                       Changes in version 1.10.2                        

- Added helper functions to compute the exact moments, so that the
user can characterise the systematic biases in the diffusion scores

                       Changes in version 1.10.1                        

- Fixed issue in Rcpp code, due to deprecation in arma, info here
- Updated some warnings and notes from BiocCheck()
- Updated readme to link to bioconductor and the journal publications

[diffUTR](/packages/diffUTR)
-------

                Changes in version 0.99.27 (2021-04-06)                 

- major speed-up or gene-level calculations

- fixed bug using the wrong default coefficient with DEXSeq

                Changes in version 0.99.13 (2021-03-04)                 

- fixed misnamed variable bug when creating annotation from ensembldb

- formatting and renaming changes to conform with Bioc standards

                Changes in version 0.99.10 (2021-02-05)                 

- submitted to BioConductor

- improvement on limma::diffSplice

- differential 3' UTR usage

[dittoSeq](/packages/dittoSeq)
--------

                         Changes in version 1.4                         

- Added 1 new Visualization function: 'dittoFreqPlot()'.

- Added interaction with 'rowData' of SE and SCEs via a
  'swap.rownames' input, e.g. to simplify provision of 'var's via
  symbols vs IDs.

- Improved & expanded 'split.by' capabilities by: 1- adding them
  to 'dittoBarPlot()', 'dittoDotPlot()', and
  'dittoPlotVarsAcrossGroups()'; 2- adding 'split.adjust' input
  to all functions for passing adjudstments to underlying
  'facet_grid()' and 'facet_wrap()' calls; 3- adding
  'split.show.all.others' input to 'dittoDimPlot()' and
  'dittoScatterPlot()' to allow the full spectrum of points,
  rather than just points excluded with 'cells.use', to be shown
  as light gray in the background of all facets; 4- Bug fix:
  splitting now works with labeling of Dim/Scatter plots, with
  label position calculated per facet, and without affecting
  facet order.

- Improved 'dittoPlot()'-plotting engine (also effects
  'dittoPlotVarsAcrossGroups()', and 'dittoFreqPlot()') by: for
  y-axis plotting, 1- extended geom dodging to also work on
  jitters when 'color.by' is used to add subgroupings & 2- added
  a 'boxplot.lineweight' control option; for x-axis / ridge
  plotting, 1- added an alternative histogram-shaping option (Try
  'ridgeplot.shape = "hist"') & 2- improved use of white space
  via a new 'ridgeplot.ymax.expansion' input.

- Standardized output logic so that 'do.hover = TRUE' will lead
  to plotly conversion even when 'data.out = TRUE'.

- 'dittoHeatmap()': 'order.by' can also now accept multiple
  gene/metadata names to order by & bug fix: when given an
  integer vector, that vector will be used directly to set the
  order of heatmap columns.

- 'dittoBarPlot()': grouping & 'var' order control improved via
  addition of a 'retain.factor.levels' input.

[DOSE](/packages/DOSE)
----

                       Changes in version 3.17.1                        

- support setting seed for fgsea method if e.g. gseGO(seed = TRUE)
(2020-10-28, Wed)
- https://github.com/YuLab-SMU/DOSE/issues/45

[DropletUtils](/packages/DropletUtils)
------------

                       Changes in version 1.12.0                        

- Added BPPARAM= to read10xCounts() for parallelized reading of
  multiple samples.

- Gave all the *Ambience() functions better names, and
  soft-deprecated the current versions.

- Added ambientContribSparse() to estimate the ambient
  contribution under sparsity assumptions.

- Added cleanTagCounts() to remove undesirable barcodes from tag
  count matrices.

- Converted all matrix-accepting functions to S4 generics to
  support SummarizedExperiment inputs.

- emptyDrops() will now coerce all DelayedArray inputs into
  wrapped SparseArraySeeds.

- Setting test.ambient=TRUE in emptyDrops() will no longer alter
  the FDRs compared to test.ambient=FALSE. Added test.ambient=NA
  to retain back-compatible behavior.

- Bugfix for correct use of redefined lower when by.rank= is set
  in emptyDrops().

- Added a constant.ambient=TRUE option to hashedDrops() to better
  support experiments with very few HTOs.

[drugTargetInteractions](/packages/drugTargetInteractions)
----------------------

                 Changes in version 0.99.0 (2021-01-11)                 

- Submitted to Bioconductor

[Dune](/packages/Dune)
----

                 Changes in version 1.3.01 (2020-11-06)                 

- Dune now accepts multiple metrics

- Dune now uses the NMI by default

[dupRadar](/packages/dupRadar)
--------

                       Changes in version 1.21.2                        

- New pkgdown documentation

[easyRNASeq](/packages/easyRNASeq)
----------

                       Changes in version 2.27.1                        

- DESeq dependency removal

- Added extra warning about RPKM usage

- Removed Defunct functions

[edgeR](/packages/edgeR)
-----

                       Changes in version 3.34.0                        

- 
  New function featureCounts2DGEList() that converts results from
  Rsubread::featureCounts() to DGELists.

- 
  Remove the "ndim" argument of plotMDS.DGEList().

- 
  read10X() now counts the number of comment lines in mtx files
  and skips those lines when reading in the data.

- 
  Fix a bug in voomLmFit() whereby zeros were sometimes
  incorrectly identified due to floating point errors.

- 
  The "bcv" method of plotMDS.DGEList() is scheduled to be
  deprecated in a future release of edgeR.

[EnhancedVolcano](/packages/EnhancedVolcano)
---------------

                        Changes in version 1.10                         

- added functionality to parse expressions in labels via parseLabels
  (TRUE/FALSE)

- over-rides ggrepel's new default value for max.overlaps via
  introduction of
  maxoverlapsConnectors = 15

- user can now specify a direction for connectors via
  directionConnectors

- removed labhjust and labvjust

- added pCutoffCol (via Andrea Grioni)

[EnMCB](/packages/EnMCB)
-----

                        Changes in version 1.2.3                        

NEW FEATURES

- Update vignettes.

- This is a major release update for EnMCB package.

- We add new options for selecting the correlation methods.

- We add mboost algorithm in our ensemble prediction model.

SIGNIFICANT USER-VISIBLE CHANGES

- Delete unused data.

- Correct the parameters' names.

BUG FIXES

- Some fixes related to test functions.

[EnrichedHeatmap](/packages/EnrichedHeatmap)
---------------

                       Changes in version 1.21.1                        

- if the matrix is non-negative, after smoothing, negative values are
  reset to zero.

[EnrichmentBrowser](/packages/EnrichmentBrowser)
-----------------

                       Changes in version 2.22.0                        

- GO gene sets: option for hierarchical annotation (new argument
  `hierarchical` for function `getGenesets`)

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.11.3                        

- Reconstruct the emapplot function and replace emapplot_cluster by
emapplot(group_category = TRUE)
- fix bug in emapplot_cluster.enrichResult when the number of cluster
is 2 (2021-2-24, Wed)
- fix bug in treeplot: The legend is not the right size (2021-2-6,
Sat).
- fix dotplot for label_format parameter doesn't work(2021-2-3, Wed).
- fix bug in gseaplot2(2021-1-28, Thu)

                       Changes in version 1.11.2                        

- update document (2021-1-7, Thu)
- update dotplot: replace ggsymbol::geom_symbol with
ggstar::geom_star(2021-1-6, Wed)
- add parameter shadowtext for three functions: emapplot,
emapplot_cluster and cnetplot. (2021-1-5, Tue)
- update dotplot: supports the use of shapes and line colors to
distinguish groups (2021-1-3, Sun)
- add treeplot function (2020-12-29, Tue)
- rename function get_ww to get_similarity_matrix (2020-12-29, Tue)
- move the emapplot related functions to emapplot_utilities.R
- fix bug in emapplot and cnetplot when enrichment result is one line
(2020-12-26, Sat)
- fix pairwise_termsim for the bug of repeated filtering of
showCategory(2020-12-23, Wed)
- fix showCategory for cnetplot, emapplot, emapplot_cluster when
showCategory is a vector of term descriptions

                       Changes in version 1.11.1                        

- add orderBy and decreasing parameters for ridgeplot() (2020-11-19,
Thu)
- https://github.com/YuLab-SMU/enrichplot/pull/84/
- update emapplot_cluster() to label cluster in center by default and
use ggrepel if setting repel = TRUE (2020-11-08, Mon)
- https://github.com/YuLab-SMU/enrichplot/pull/81
- add a label_format parameter to support formatting label
(2020-10-28, Wed)
- if provided with a numeric value will simply string wrap by
default
- if provided with a function will instead set labels =
user_defined_function() within the scale function
- https://github.com/YuLab-SMU/enrichplot/pull/73

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.15.3                        

- Fix missing declaration of rmarkdown.

                       Changes in version 2.15.2                        

- Ensure remote gzipped files are handled properly by `ensDbFromGtf`.

                       Changes in version 2.15.1                        

- Add new field canonical_transcript to the gene table reporting the ID
  of the
  gene's canonical transcript.

[ensemblVEP](/packages/ensemblVEP)
----------

                       Changes in version 1.33.0                        

- add support for Ensembl release 102/103/104

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 0.99.0 (2021-04-09)                 

- R>=4.0 for submission

- removed unused dependencies

- correct work of generateVcfReport (although SNV only)

- unmatched reads are at the end of generateBed* output now

- compiles and works on Apple Silicon (native ARM64 R)

- fully documented methods

- fully covered with tests and examples

- comprehensive vignettes

                 Changes in version 0.4.0 (2021-03-08)                  

- going public

- CX report now includes only the context present in more than 50% of
  the reads

- generateVcfReport (capable of dealing with SNVs only for now)

- added documentation to some of the methods

- added several examples

- added sample data for amplicon and capture NGS

- added some tests based on sample data

- README.md

                 Changes in version 0.3.9 (2021-01-19)                  

- fast C++ CIGAR parser to lay queries in reference space

- new method to extract base frequences: generateBaseFreqReport

                 Changes in version 0.3.7 (2021-01-12)                  

- lots of refactoring again

- CX report sub now uses boost::container::flat_map (additional 2x
  speedup)

- removed dplyr as a dependence, whole package uses data.table now

                 Changes in version 0.3.5 (2021-01-09)                  

- lots of refactoring

- new method: preprocessBam() to save time on loading/preprocessing

- new C++ sub for CX report with std::map summary (5-10x speedup)

                 Changes in version 0.3.2 (2021-01-06)                  

- first attempt to stablilize API (generateCytosineReport and
  generateBedReport)

- temporary method for ECDF (generateBedEcdf)

- uploaded to GitHub

                 Changes in version 0.3.1 (2020-01-01)                  

- heavy refactoring, many internal methods added

- C++ functions for nearly all bottlenecks (pending fast: cigar,
  summary, genome loading)

                 Changes in version 0.2.1 (2020-12-21)                  

- made this second iteration of epialleleR a usable package

[epigraHMM](/packages/epigraHMM)
---------

                       Changes in version 0.99.0                        

- First release of epigraHMM.

- It is now possible to add normalizing offsets via addOffsets.

- epigraHMM now uses hdf5 files to store all intermediate data during
computation of the EM algorithm. Intermediate data include
window-based HMM and mixture model posterior probabiltiies, and
forward-backward probabilities. This change leads to a better memory
utilization upon convergece.


[escape](/packages/escape)
------

                        Changes in version 1.0.1                        

- Removed ggrepel, rlang, and factoextra dependencies.

- Updated Seurat package switch

- Switch the way counts are processed by first eliminating rows with 0
  expression in the sparse matrix before converting to a full matrix

[evaluomeR](/packages/evaluomeR)
---------

                        Changes in version 1.7.5                        

- Minor bug fix for Mac compilation.

                 Changes in version 1.7.4 (2020-06-01)                  

- Added 'scale' parameter to plotMetricsCluster method.

                 Changes in version 1.7.3 (2020-04-16)                  

- Clusterboot interfaces can be set through 'cbi' parameter to quality
  methods.
  It takes one the following values: "kmeans", "clara", "clara_pam",
  "hclust", "pamk", "pamk_pam".

                 Changes in version 1.7.2 (2020-04-14)                  

- Stability analysis and quality analysis will not stuck in bootstrap.

                 Changes in version 1.7.1 (2020-04-13)                  

- Clusterboot interfaces can be set through 'cbi' parameter to
  stability methods.
  It takes one the following values: "kmeans", "clara", "clara_pam",
  "hclust", "pamk", "pamk_pam".

[EWCE](/packages/EWCE)
----

                        Changes in version 1.0.0                        

New Features

- EWCE v1.0 on Bioconductor replaces the defunct EWCE v1.3.0
available
on Bioconductor v3.5.
- EWCE has been rendered scalable to the analysis of large datasets
- drop_uninformative_genes() has been expanded to allow the
utilisation of differential expression approaches
- EWCE can now handle SingleCellExperiment (SCE) objects or other
Ranged SummarizedExperiment (SE) data types and as input as well as
the original format, described as a single cell transcriptome (SCT)
object.

Deprecated & Defunct

- The following functions have been renamed to use underscore in
compliance with Bioconductor nomenclature:
check.ewce.genelist.inputs
cell.list.dist
bootstrap.enrichment.test
bin.specificity.into.quantiles
bin.columns.into.quantiles
add.res.to.merging.list
prepare.genesize.control.network
prep.dendro
get.celltype.table
calculate.specificity.for.level
calculate.meanexp.for.level
generate.celltype.data
generate.bootstrap.plots
generate.bootstrap.plots.for.transcriptome
fix.bad.mgi.symbols
fix.bad.hgnc.symbols
filter.genes.without.1to1.homolog
ewce.plot
cells.in.ctd
drop.uninformative.genes

[exomePeak2](/packages/exomePeak2)
----------

                 Changes in version 1.3.7 (2021-03-10)                  

- Parameter `parallel` is changed in the exomePeak2() and
  exomePeakCalling() functions; the parameter now enables user to
  configure specific number of cores used in the parallel computation
  (default = 1).

- To avoid the potential confusion for the downstream analysis, the
  default settings for the parameter `log2FC_cutoff` in functions
  exomePeak2() and exomePeakCalling() are changed from 1 to 0. The
  adjustment should have very little effect on the peak calling result.

- The naming of peaks in the output file is now sorted by their genomic
  order.

- A maximum for peak width is added now, which is by default
  100*fragment_length. Such a higher bound can significantly improve
  the results of DRACH motif finding for m6A-Seq.

                 Changes in version 1.3.5 (2021-02-07)                  

- Improved the grammar and details in the DESCRIPTION file and the
  vignettes.

- When performing the difference analysis using the function
  exomePeak2(), the sequencing depth of the interactive GLM will be
  estimated on the background features, which by default are the
  disjoint regions of the peaks detected on the exons. Tests on real
  data revealed that the background approach can make the differential
  methylation directions more in-line with the expectation of the
  perturbed protein regulator. Previously, the background sequencing
  depth estimation can only be realized in the multiple-step functions
  but not in exomePeak2().

                 Changes in version 1.3.4 (2021-02-05)                  

- The options `consistent_peak`, `consistent_log2FC_cutoff`,
  `consistent_fdr_cutoff`, `alpha`, and `p0` are deprecated from the
  functions exomePeak2() and exomePeakCalling(). The consistent_peak
  option was implemented to reproduce the consistent peak calling
  algorithm in the old package exomePeak, and its performance is
  significantly lower than the NB GLM derived methods according to our
  recent tests. Hence, the consistency related functionalities are
  removed in the later versions of exomePeak2.

                 Changes in version 1.3.3 (2021-02-03)                  

- Fix the bug of not merging the overlapping exons when the transcript
  annotation have no overlapping transcripts, this can happen when a
  very small annotation is provided.

[ExperimentHub](/packages/ExperimentHub)
-------------

                       Changes in version 1.99.0                        

MAJOR UPDATES

- (1.99.0) The default caching location has changed. Instead of
  rappdirs::user_cache_dir using tools::R_user_dir. To avoid
  conflicting
  caches, a user will have to manage an old cache location before
  proceeding. Information for handling an old cache location is
  provided in
  the vignette.

- (1.99.0) Another major change, a default caching location is
  automatically
  created in a non interactive session instead of using a temporary
  location. In an interactive session, a user is still prompted for
  permission.

                       Changes in version 1.17.0                        

- (1.17.1) Removed vignette for creating annotation hub package.
  Reference
  and refer to single vignette in AnnotationHub

[ExperimentHubData](/packages/ExperimentHubData)
-----------------

                       Changes in version 1.17.0                        

MODIFICATIONS

- 1.17.2 Removed vignette for creating annotation hub package.
  Reference
  and refer to single vignette in AnnotationHub

- 1.17.1 Tags for database now combination of biocViews and meta$Tags.
  Also
  checks for valid AnnotationHub or AnnotationHubSoftware biocViews.

BUG CORRECTION

- 1.17.1 Fixed bug to run make*HubMetadata using ".". Fixed in
  AnnotationHubData. bumped dependency

[ExperimentSubset](/packages/ExperimentSubset)
----------------

                 Changes in version 1.1.0 (2021-05-09)                  

- Added support for TreeSummarizedExperiment and SpatialExperiment
classes

[ExploreModelMatrix](/packages/ExploreModelMatrix)
------------------

                        Changes in version 1.3.2                        

- Enable MathJax in tour

[FamAgg](/packages/FamAgg)
------

                       Changes in version 1.19.1                        

- Add kinshipPairs function

[famat](/packages/famat)
-----

                 Changes in version 1.3.0 (2020-11-27)                  

- Rshiny modifications for figure in paper

[fcoex](/packages/fcoex)
-----

                 Changes in version 1.5.5 (2021-05-12)                  

- Changed fcoex object to store a dgCMatrix instead of a dataframe for
  expression and discretized expression.

- Removed options of 3+ class multiclass discretize (may lead to
  backwards incompatibility) as they would
  be incompatible with better memory handling of dgCMatrix system.

[FEAST](/packages/FEAST)
-----

                 Changes in version 0.99.0 (2021-04-01)                 

- Submitted to Bioconductor

[fedup](/packages/fedup)
-----

                 Changes in version 0.99.7 (2021-04-24)                 

- added info to explain example data in README and vignettes

- made minor changes to exported function example comments

- plotFemap + parameterized several hard-coded variables

                 Changes in version 0.99.6 (2021-02-24)                 

- Updated package data script paths

                 Changes in version 0.99.5 (2021-02-24)                 

- Trimmed external data size

                 Changes in version 0.99.4 (2021-03-24)                 

- updated package datasets (geneSingle, geneDouble, geneMulti)

- created 3 vignettes to describe package implementation using
  each dataset

- runFedup + runs analysis on an input list rather than single
  test vector + fold enrichment calculation evaluates 0 for 0/0
  instances + enriched pathways defined as fold enrichment ≥ 1
  (instead of > 1)

- writeFemap + writes EM-formatted tables for list of fedup
  results

- plotFemap + implements tryCatch() to return NULL if Cytoscape
  is not running locally

- prepInput + new function to prepare input gene list

                 Changes in version 0.99.3 (2021-02-24)                 

- Updating R version dependency to 4.1

                 Changes in version 0.99.2 (2021-02-24)                 

- Untracking system files

                 Changes in version 0.99.1 (2021-02-23)                 

- Version bump

                 Changes in version 0.99.0 (2021-02-17)                 

- Submitted to Bioconductor

[FGNet](/packages/FGNet)
-----

                        Changes in version 3.26                         

- Queries to DAVID are no longer supported.

[FilterFFPE](/packages/FilterFFPE)
----------

                 Changes in version 1.1.2 (2020-11-11)                  

- Fix bugs in keeping extra reads when filtering with minMapBase

- Allow skipping filtering with minMapBase

[fishpond](/packages/fishpond)
--------

                        Changes in version 1.8.0                        

- Added note in vignette about how to deal with estimated
  batch factors, e.g. from RUVSeq or SVA. Two strategies are
  outlined: either discretizing the estimate batch factors
  and performing stratified analysis, or regressing out the
  batch-associated variation using limma's removeBatchEffect.
  Demonstation code is included.

[flowAI](/packages/flowAI)
------

                       Changes in version 1.21.6                        

- the flow rate check is less sensitive now. Note that the default
  value of alphaFR is now 0.1 instead of 0.01

[flowCL](/packages/flowCL)
------

                       Changes in version 1.29.1                        

- Endpoint no longer exists. Please email Justin at
  justinmeskas@gmail.com to request a fix.

[flowGraph](/packages/flowGraph)
---------

                 Changes in version 0.99.0 (2020-09-30)                 

- Submitted to Bioconductor

[FlowSOM](/packages/FlowSOM)
-------

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

                       Changes in version 0.99.56                       

- New package vignette "Use case ST000291".
- New package vignette "Use case ST000629".
- Remove AppVeyor and Travis CI and move codecov to GitHub Actions.
- Updated vignette "Dietary text annotation".
- Improvements to the annotate_foods() function.

                       Changes in version 0.99.41                       

- Switch from the sigora CRAN package to the fgsea Bioconductor
package to perform enrichment analysis.
- Added new function named msea to perform GSEA using FOBI.

                       Changes in version 0.99.38                       

- Addressing Vince Carey (package reviewer from Bioconductor)
comments
and suggestions.

                       Changes in version 0.99.35                       

- Enabling the option to include a FOBI table instead of downloading
it from GitHub.
- Added FOBI table in package data.

                       Changes in version 0.99.29                       

- Fix Bioconductor Single Package Builder errors and warnings:

- fobitools.Rproj file removed
- Rd line widths set to less than 100 characters
- Examples added in ora.R
- R version dependency updated from 3.6.0 to 4.1.
- Use TRUE/FALSE instead of T/F in ora.R

                       Changes in version 0.99.24                       

- Added vignettes.
- pkgdown updated.
- Submitted to Bioconductor!

                       Changes in version 0.99.18                       

- Added a function called fobi_graph to generate FOBI graphs.

                       Changes in version 0.99.12                       

- Added a NEWS.md file to track changes to the package.

[FRASER](/packages/FRASER)
------

                        Changes in version 1.2.1                        

- Add merging of external counts

- Add publication

- Minor bugfixes

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.28.0                        

NEW FEATURES

- new function `exist.gdsn()`

- new function `is.sparse.gdsn()`

UTILITIES

- LZ4 updated to v1.9.3 from v1.9.2

- XZ is updated to v5.2.5 from v5.2.4

- `apply.gdsn()`: work around with factor variables if less-than-32-bit
  integers are stored

- a new component 'is.sparse' in `objdesp.gdsn()`

- `options(gds.verbose=TRUE)` to show additional information

                       Changes in version 1.26.1                        

UTILITIES

- comply with the R devel (> v4.0.3) to work with factor variables in
  `apply.gdsn()`

[GeneExpressionSignature](/packages/GeneExpressionSignature)
-----------------------

                       Changes in version 1.37.0                        

- Add a NEWS.md file to track changes to the package.
- Add inst/CITATION file to customise the citation.
- Add README.md, CODE_OF_CONDUCT.md, and create a biocthis-style
GitHub Actions workflow.
- Add a new function PGSEA from the PGSEA to remove the dependency on
PGSEA, since PGSEA was deprecated in Bioconductor version 3.12 and
removed from 3.13.
- Format code using styler and biocthis, redocument package using
Roxygen2.
- Rewrite vignette using knitr, rmarkdown, BiocStyle.
- Fix BiocCheck errors and warnings.

[GENESIS](/packages/GENESIS)
-------

                       Changes in version 2.21.5                        

- Added functions to compute variant-specific inflation factors.

                       Changes in version 2.21.4                        

- Added the option to perform a fast approximation to the score
  standard error in assocTestSingle. New function
  nullModelFastScore prepares a null model to be used with this
  option.

                       Changes in version 2.21.1                        

- Updated structure of fitNullModel objects. Null model objects
  with the previous structure will be automatically updated with
  a warning, but you may want to consider rerunning
  `fitNullModel` if you plan to use an older null model with the
  current version.

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 1.4.0                        

New features

- The main function GeneTonic() gains an extra parameter, gtl - this
can be used to provided a named list object where a single parameter
is passed (e.g. after loading in a single serialized object), while
the functionality stays unaltered. The same gtl parameter is also
exposed in other functions of the package - see the vignette for
some examples, or check the documentation of each specific function.
To create this object in a standardized manner, the function
GeneTonic_list() is now available.

- A new function to perform fuzzy clustering (following the
implementation of DAVID) is added - see gs_fuzzyclustering(). It
returns a table with additional information on the cluster of
genesets and the status of each set in the group.

- The ggs_backbone() function can extract the bipartite graph
backbone
from the Gene-Geneset graph, this can be further explored below the
main element in the Gene-Geneset panel. Once the backbone is
created, you are one step away from checking out the genes that act
as "hubs" in the Gene-Geneset graph, and possibly identify the nodes
playing an essential role based on their connectivity.

- A new function, signature_volcano(), adds a signature volcano plot
to the Gene-Geneset panel. This plot displays the genes of a chosen
geneset in color, while the remaining genes of the data are shown as
shaded dots in the background. The color and transparency of the
displayed genes can be chosen by the user, as well as the option to
display the gene names of all genes in the geneset.

- gs_summary_overview() can also generate bar plots instead of the
default segment-dot (lollipop) plots.

- A new function, summarize_ggs_hubgenes(), builds a DT datatable for
the Gene-Geneset panel. This table lists the individual genes of the
input data and their respective degree in the Gene-Geneset graph.
Furthermore, action buttons linking to the NCBI, GeneCards and GTEx
databases are included for each gene.

- gene_plot() gains the extra labels_display argument to control
whether the labels are at all shown; now the display of the labels
is also respecting the jitter of the points

Other notes

- gs_heatmap() has now the possibility to set the arguments to the
call to heatmap generating function, via ellipsis
- gs_heatmap() handles the colors in a consistent way over the
different executions, without relying on the random palettes
provided by the Heatmap's annotation functionality - could have been
misleading if encountering too similar hues are randomly picked
- the plots obtained via gs_mds() and gs_volcano() now always display
the line segments for the data points to be labeled (increasing the
readability - as "matching back the label to the drawed circle" -
thanks for the suggestion!)

[GENIE3](/packages/GENIE3)
------

                       Changes in version 1.13.3                        

- The regulators can now be provided as a list with different
  regulators for each gene/feature.

[GenomicOZone](/packages/GenomicOZone)
------------

                        Changes in version 1.5.1                        

- Fixed a bug caused by a deprecated function sjstats::eta_sq().
  + In 'MD_perform_zoning.R', function 'sjstats::eta_sq()' is
  deprecated. Use 'lsr::etaSquared()' instead.
  + In 'DESCRIPTION', removed 'sjstats' from 'Imports'; added 'lsr'
  into 'Imports'; removed 'GEOquery' from 'Suggests'.
  + No longer downloading the data from GEO in the vignettes. Added the
  data file in inst/extdata. The data is only used by the vignettes.
  + 'NAMESPACE' is re-auto-generated by Roxygen.

[GenomicScores](/packages/GenomicScores)
-------------

                        Changes in version 2.4.0                        

USER VISIBLE CHANGES

- The gscores() function now returns the SeqInfo from the input GScores
  object.

- Improvements on the shiny web app.

[GenomicSuperSignature](/packages/GenomicSuperSignature)
---------------------

                        Changes in version 1.0.0                        

- Initial release of the 'GenomicSuperSignature' package

[GEOfastq](/packages/GEOfastq)
--------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

                        Changes in version 0.6.5                        

- Resolved NOTES from BiocCheck::BiocCheck()
- added unit tests
- replaced sapply with vapply
- shortened lines
- added a NEWS.md file to track changes to the package.

[geva](/packages/geva)
----

                 Changes in version 0.99.0 (2021-02-01)                 

- Submitted to Bioconductor

[ggtree](/packages/ggtree)
------

                        Changes in version 2.5.3                        

- optimize text angle in geom_cladelab (2021-05-10, Mon)
- https://github.com/YuLab-SMU/ggtree/pull/396

                        Changes in version 2.5.2                        

- extend 'continuous' parameter to support 4 possible values, i.e.,
'none' to disable continuous transition, 'color' (or 'colour') to
enable continuous color transition, 'size' to enable continuous size
(branch thickness) transition and 'all' to enable continuous color
and size transition (2021-04-07, Wed)
- https://github.com/YuLab-SMU/ggtree/pull/385
- https://github.com/YuLab-SMU/ggtree/pull/387
- extendto argument for geom_hilight now compatible with
'inward_circular' and 'dendrogram' layouts (2021-02-25, Thu)
- https://github.com/YuLab-SMU/ggtree/pull/379

                        Changes in version 2.5.1                        

- update man file of geom_rootpoint (2021-01-08, Fri)
- label and offset.label introduced in geom_treescale layer
(2020-12-23, Wed)
- https://github.com/YuLab-SMU/ggtree/pull/360
- geom_rootedge supports reversed x (2020-12-17, Thu)
- https://github.com/YuLab-SMU/ggtree/pull/358
- geom_nodelab() now supports circular layout (2020-11-26, Thu)
- https://github.com/YuLab-SMU/ggtree/issues/352
- https://github.com/YuLab-SMU/ggtree/pull/353
- branch size can be grandualy changed (2020-10-29, Thu)
- https://github.com/YuLab-SMU/ggtree/pull/349

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                       Changes in version 1.1.12                        

- import ggtree to pass BiocCheck. (2021-05-14, Fri)

                       Changes in version 1.1.11                        

- fix a bug to solve the problem (variable of x has NA). (2021-05-13,
Thu)
- https://github.com/YuLab-SMU/ggtreeExtra/issues/8

                       Changes in version 1.1.10                        

- fix a bug for compute_aes ( This is to support mapping aesthetics
(x, not y in aes of geom_fruit) to functions of variables).
(2021-05-10, Mon)

                        Changes in version 1.1.9                        

- don't inherit aes (global aes from the ggtree). (2021-04-28, Wed)
- subset in mapping also supports the data from ggtree object.
(2021-04-29, Thu)
- support mapping aesthetics (x, not y in aes of geom_fruit) to
functions of variables. (2021-05-07, Fri)

                        Changes in version 1.1.8                        

- add new position functions: position_jitterx and
position_jitterdodgex. (2021-04-23, Fri)
- update vignettes. (2021-04-25, Sun)

                        Changes in version 1.1.7                        

- update man and vignettes. (2021-04-06, Tue)

                        Changes in version 1.1.6                        

- check whether the value of x is numeric to avoid warnings when x is
factor. (2021-02-24, Wed)
- remove axis of first geom_tile of vignettes, since the axis of this
layer is meaningless. (2021-02-24, wed)

                        Changes in version 1.1.5                        

- support title of panel. (2021-02-03, Wed)
- https://github.com/YuLab-SMU/ggtreeExtra/issues/7
- add citation info. (2021-02-04, Thu)
- don't use svg dev. (2021-02-04, Thu)

                        Changes in version 1.1.4                        

- add position_points_jitterx and position_raincloudx for geom of
ggridges. (2021-01-20, Wed)
- supports geom_msa of ggmsa. (2021-01-21, Thu)
- specific position method for specific geom method automatically.
(2021-01-27, Wed)

                        Changes in version 1.1.3                        

- support multiple density plot from geom of ggridges. (2020-12-31,
Thu) geom_density_ridges, geom_density_ridges2,
geom_density_ridges_gradient, geom_ridgeline,
geom_ridgeline_gradient.

                        Changes in version 1.1.2                        

- support subset in mapping, but the data should also be provided.
(2020-11-30, Mon)
- add default position methods for common geometric functions.
(2020-12-18, Fri)

[GlobalAncova](/packages/GlobalAncova)
------------

                        Changes in version 4.9.1                        

- removed function 'GAKEGG' for testing collection of pathways, due to
  deprecation of KEGG package

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.17.1                        

- bug fixed according to the update of GO.db (2020-10-29, Thu)
- https://github.com/YuLab-SMU/GOSemSim/issues/32

[granulator](/packages/granulator)
----------

                 Changes in version 0.99.0 (2021-03-30)                 

- Submitted to Bioconductor

[graphite](/packages/graphite)
--------

                 Changes in version 1.37.1 (2020-05-04)                 

- Removed Biocarta and NCI pathways.

- Added WikiPathways pathways.

- Updated all pathway data.

[GSEABase](/packages/GSEABase)
--------

                        Changes in version 1.54                         

SIGNIFICANT USER-VISIBLE CHANGES

- orgDb no longer provide Unigene information. Remove support.

[GSVA](/packages/GSVA)
----

                        Changes in version 1.40                         

USER VISIBLE CHANGES

- The vignette has been rewritten in R Markdown to produce an HTML
  vignette page and make it shorter and faster to produce.

- Development of a shiny app available through the function 'igsva()'.

BUG FIXES

- Replaced fastmatch::fmatch() by IRanges::match,CharacterList-method
  after disscussion at https://github.com/rcastelo/GSVA/issues/39 to
  avoid the row names of an input expression matrix being altered by
  fastmatch::fmatch() adding an attribute.

- Fixed wrong call to .mapGeneSetsToFeatures() when gene sets are given
  in a GeneSetCollection object.

[Heatplus](/packages/Heatplus)
--------

                       Changes in version 2.99.4                        

- Improved label- and axis handling for panels with continuous
  covariates

- Formally deprecated heatmap_2 and heatmap_plus

- Refactored to switch to devtools tool chain

[Herper](/packages/Herper)
------

                 Changes in version 1.1.1 (2021-02-19)                  

- Include option to build miniconda in import from yml, and update
  vignette

[HGC](/packages/HGC)
---

                 Changes in version 0.99.3 (2021-05-14)                 

- Remove build/ folder.

                 Changes in version 0.99.2 (2021-04-20)                 

- Remove some abundant files.

                 Changes in version 0.99.1 (2021-04-18)                 

- Remove some abundant files.

                 Changes in version 0.99.0 (2021-04-10)                 

- The first version 0.99.0 is submitted to Bioconductor

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.28.0                        

- `hlaPredict()` returns the dosage of HLA alleles when
  type="response+dosage", and
  `hlaPredict()` returns the best guess and dosages by default

- a new option "Pos+Allele" in `hlaPredict()`, `hlaGenoCombine()`,
  `hlaGenoSwitchStrand()`, `hlaSNPID()` and `hlaCheckSNPs()` for
  matching genotypes
  by positions, reference and alternative alleles; it is particularly
  useful when
  the training and test set are both matched to the same reference
  genome,
  e.g., 1000 Genomes Project

- `hlaGDS2Geno()` supports SeqArray GDS files

- a new option 'maf' in `hlaAttrBagging()` and
  `hlaParallelAttrBagging()`

- 'pos.start' and 'pos.end' are replaced by 'pos.mid' in
  `hlaFlankingSNP()` and
  `hlaGenoSubsetFlank()`

- new function `hlaAlleleToVCF()` for converting the imputed HLA
  classical alleles
  to a VCF file

                       Changes in version 1.26.1                        

- the hlaAttrBagging object can be removed in garbage collection
  without
  calling `hlaClose()`

- enable internal GPU API

- improved multithreaded performance compared with v1.26.0

[HiCDCPlus](/packages/HiCDCPlus)
---------

                Changes in version 0.99.14 (2021-04-13)                 

- Fixed bug that prevents passing negative values to .hic files

                Changes in version 0.99.13 (2021-03-23)                 

- Added RE-agnostic features into construct_features

                Changes in version 0.99.12 (2021-02-19)                 

- Submitted to Bioconductor

[HilbertCurve](/packages/HilbertCurve)
------------

                       Changes in version 1.21.1                        

- add "max_freq" method for mean_mode so that in each interval, value
  corresponded to
  the maximal frequency is selected.

[HPAanalyze](/packages/HPAanalyze)
----------

                         Changes in version 1.9                         

- Changes in version 1.9.5
  + hpaVisTissue now plots all tissues by default
  + hpaDownload shortcuts for downloadList arg can now be combined
  + Bug fixes and performance improvements

- Changes in version 1.9.4
  + hpaVisPatho can now facet by cancers or genes
  + Vignettes update
  + Bug fixes and performance improvements

- Changes in version 1.9.3
  + Update the default color to be more accessible (based on
  viridis::magma)
  + Updated hpaSubset and hpaListParam to work properly with new data

- Changes in version 1.9.2
  + hpaDownload can now download every downloadable datasets from HPA.
  + Update built-in data to HPA version 20.1.

- Changes in version 1.9.1
  + Update built-in data to HPA version 20.

- Changes in version 1.9.0
  + Starting devel for Bioconductor 3.13

[hpar](/packages/hpar)
----

                        Changes in version 1.33                         

Changes in version 1.33.2

- Update to HPA release 20.0 <2020-11-24 Thu>

Changes in version 1.33.1

- Update to HPA release 19.3 (2020.03.06) <2020-10-26 Mon>
- New release for Bioconductor devel 3.12 <2020-10-26 Mon>
- Updated R version (>= 3.5.0) in Depends field <2020-10-26 Mon>
- Added the Secretome data <2020-10-28 Wed>
- Updated the documentation and docs folder for pkgdown <2020-10-29
Thu>
- Automated allHparData() <2020-10-29 Thu>

Changes in version 1.33.0

- New Bioc devel version

[HPAStainR](/packages/HPAStainR)
---------

                        Changes in version 1.0.4                        

- The first update in response to F1000 comments

- Added dating arguments to HPA_data_downloader

- Now when you have save_file set to true, the file will be stamped
  with
  the date it was downloaded. This allows for reproducible access tot
  he
  file you created.

- To deal with dated files we now have the arguments
  `version_date_normal` and `version_date` cancer, which allows users
  to
  select which downloaded HPA data they want to use based on the date
  of
  downloading. The default to "last" which will search for the latest
  version of the files in the `save_location`

- The last addition is `force_download`, an argument that will overide
  the function's default usage of the local files in case you want to
  download a more recent version of the HPA data.

- Change to HPAStainR and by extension shiny_HPAStainR

- Added Fisher's Exact Test to HPAStainR's main function

- Due to the inconsistency that exists in using simulated p-values in
  `chisq.test` the default new test if Fisher's exact and an argument
  `test_type` has been added so users can pick between the two.

- Changed the output of the p-values to numeric instead of a character.

                 Changes in version 1.0.3 (2021-02-03)                  

- Changed read.table() in HPA_data_downloader.R to data.table's fread()

                 Changes in version 1.0.2 (2021-25-20)                  

- Changed section of code crashing due to dplyr update.

                 Changes in version 1.0.1 (2020-11-20)                  

- testthat() HPA_data_downloader.R failed due to spelling change

- HPA changed "unfavourable" to "unfavorable"

- The testthat() has been changed to reflect their change

- More updates soon after all reviews from F1000 are in

[HTSFilter](/packages/HTSFilter)
---------

                       Changes in version 1.31.1                        

- -- Remove all references and functionality related to the
  deprecated DESeq package.
  -- Vignette has been updated to Rmarkdown from Sweave.

[HubPub](/packages/HubPub)
------

                 Changes in version 0.99.0 (2021-04-23)                 

- Submitted to Bioconductor

[ideal](/packages/ideal)
-----

                       Changes in version 1.16.0                        

New features

- It is now possible to export the processed data (count data,
results, functional enrichment table) into a combined list object,
which can be seamlessly fed into GeneTonic
(http://bioconductor.org/packages/release/bioc/html/GeneTonic.html).
- Adjusted the behavior for the modeling with LRT, with extra
notifications to inform the user on how to make the most out of the
functionality.

Other notes

- The manuscript of ideal is now published in BMC Bioinformatics! The
citation file has been updated accordingly.
- A full round of styler has been applied to the codebase.

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.23.6                        

- CHANGES
  * code cleaning

                       Changes in version 1.23.3                        

- CHANGES
  * fixes a testthat misspelling

                       Changes in version 1.23.2                        

- CHANGES
  * tolerance sufficience in cell-subclustering not required for first
  model
  refinement test

                       Changes in version 1.23.1                        

- NEW FEATURES
  * added methods clusterDist, clusterProb, clusterCoeff for
  immunoMeta-object
  * starting with unit tests

[immunotation](/packages/immunotation)
------------

                 Changes in version 0.99.7 (2021-05-04)                 

- Increased number of times to try to connect to webresource

                 Changes in version 0.99.6 (2021-05-04)                 

- Adressed comments from Bioconductor review

- included a mro.obo.gz file which can be read without unzipping

                 Changes in version 0.99.1 (2021-01-01)                 

- Submitted to Bioconductor

[infercnv](/packages/infercnv)
--------

                 Changes in version 1.7.2 (2020-05-05)                  

- New dependencies : RANN, leiden, phyclust

- Added new partition method for subclustering that uses the Leiden
  algorithm based on a K-nn adjacency matrix.

- Changed the default subclustering method to leiden which is much
  faster than the random trees method.

- Split Bayesian filtering step in two steps, one that runs the model,
  and one that applies the filter threshold. This allows updating the
  threshold without having to rerun the whole model.

- Fix what groupings of references the subclustering is done on.

- Updated expectations of the internally stored clustering information
  when plotting references to allow for results obtained with version
  between ~1.3 and this one to be plotted.

- Fix add_to_seurat method to work when no seurat object is provided
  after the reordering fix.

- Fix the random trees subclustering applying a different method of
  centering to the data between the hclust stored in
  infercnv_obj@tumors_subclusters$hc and the splits in
  infercnv_obj@tumors_subclusters$subclusters.

- Make denoising step figure only be generated if plot_steps is true.
  (it is identical to final figure that is plotted based on a different
  option, making it redundant)

[Informeasure](/packages/Informeasure)
------------

                 Changes in version 1.1.1 (2020-10-30)                  

- make changes

[InPAS](/packages/InPAS)
-----

                       Changes in version 1.99.4                        

- update get_PAscore2.

                       Changes in version 1.99.3                        

- add rmarkdown as suggests.

                       Changes in version 1.99.2                        

- fix a bug if utr3 list is empty.

                       Changes in version 1.99.1                        

- add dontrun for getGCandMappability doc.

                       Changes in version 1.99.0                        

- merge Haibo's code.

[InteractiveComplexHeatmap](/packages/InteractiveComplexHeatmap)
-------------------------

                       Changes in version 0.99.10                       

- add `response` argument so that the server can only respond to one
  event from UI.

                       Changes in version 0.99.9                        

- output can be floating along with mouse positions.

                       Changes in version 0.99.8                        

- click and hover won't conflict with brush.

                       Changes in version 0.99.7                        

- In the sub-heatmap, it allows to remove rows and columns from the
  four sides.

                       Changes in version 0.99.0                        

- Submit to Bioconductor

[InterCellar](/packages/InterCellar)
-----------

                     Changes in version 0.0.0.9000                      

- Added a NEWS.md file to track changes to the package.

[IRISFGM](/packages/IRISFGM)
-------

                       Changes in version 0.99.8                        

- Submitted to Bioconductor

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.1.11 (2021-05-11)                 

NEW FUNCTIONALITY

- HSC_population_size_estimate and HSC_population_plot allow
estimates
on hematopoietic stem cell population size
- Importing of Vispa2 stats per pool now has a dedicated function,
import_Vispa2_stats
- outlier_filter and outliers_by_pool_fragments offer a mean to
filter
poorly represented samples based on custom outliers tests

VISIBLE USER CHANGES

- The argument import_stats of aggregate_metadata is officially
deprecated in favor of import_Vispa2_stats
- aggregate_metadata is now a lot more flexible on what operations
can
be performed on columns via the new argument aggregating_functions
- import_association_file allows directly for the import of Vispa2
stats and converts time points to months and years where not already
present
- File system alignment of import_association_file now produces 3
separate columns for paths
- separate_quant_matrices and comparison_matrix now do not require
mandatory columns other than the quantifications - this allows for
separation or joining also for aggregated matrices

FIXES

- Fixed a minor issue in CIS_volcano_plot that caused duplication of
some labels if highlighted genes were provided in input

                 Changes in version 1.1.10 (2021-04-08)                 

FIXES

- Fixed issue in compute_near_integrations: when provided
recalibration map export path as a folder now the function works
correctly and produces an automatically generated file name
- Fixed issue in aggregate_metadata: now paths to folder that
contains
Vispa2 stats is looked up correctly. Also, VISPA2 stats columns are
aggregated if found in the input data frame independently from the
parameter import_stats.

IMPROVEMENTS

- compute_abundance can now take as input aggregated matrices and has
additional parameters to offer more flexibility to the user. Major
updates and improvements also on documentation and reproducible
examples.
- Major improvements in function import_single_Vispa2Matrix: import
is
now preferentially carried out using data.table::fread greatly
speeding up the process - where not possible readr::read_delim is
used instead
- Major improvements in function import_association_file: greatly
improved parsing precision (each column has a dedicated type),
import report now signals parsing problems and their location and
signals also problems in parsing dates. Report also includes
potential problems in column names and signals missing data in
important columns. Added also the possibility to give various file
formats in input including *.xls(x) formats.
- Function top_integrations can now take additional parameters to
compute top n genes for each specified group
- Removed faceting parameters in CIS_volcano_plot due to poor
precision (easier to add faceting manually) and added parameters to
return the data frame that generated the plot as an additional
result. Also, it is now possible to specify a vector of gene names
to highlight even if they're not above the annotation threshold.

MINOR

- ISAnalytics website has improved graphic theme and has an
additional
button on the right that leads to the devel (or release) version of
the website
- Updated vignettes

FOR DEVS ONLY

- Complete rework of test suite to be compliant to testthat v.3

                 Changes in version 1.1.9 (2021-02-17)                  

FIXES

- Fixed minor issues in internal functions with absolute file paths &
corrected typos

                 Changes in version 1.1.8 (2020-02-15)                  

FIXES

- Fixed minor issues in internal functions to optimize file system
alignment

                 Changes in version 1.1.7 (2020-02-10)                  

FIXES

- Fixed minor issues in import_association_file when checking
parameters

                 Changes in version 1.1.6 (2020-02-06)                  

UPGRADES

- It is now possible to save html reports to file from
import_parallel_Vispa2Matrices_auto and
import_parallel_Vispa2Matrices_interactive, remove_collisions and
compute_near_integrations

FIXES

- Fixed sample_statistics: now functions that have data frame output
do not produce nested tables. Flat tables are ready to be saved to
file or can be nested.
- Simplified association file check logic in remove_collisions: now
function blocks only if the af doesn't contain the needed columns

                 Changes in version 1.1.5 (2020-02-03)                  

UPGRADES

- Upgraded import_association_file function: now file alignment is
not
mandatory anymore and it is possible to save the html report to file
- Updated vignettes and documentation

                 Changes in version 1.1.4 (2020-11-16)                  

UPGRADES

- Greatly improved reports for collision removal function
- General improvements for all widget reports

                 Changes in version 1.1.3 (2020-11-10)                  

FIXES

- Further fixes for printing reports when widgets not available
- Added progress bar to collision processing in remove_collisions
- Updated vignettes

NEW

- Added vignette "Using ISAnalytics without RStudio support"

                 Changes in version 1.1.2 (2020-11-05)                  

FIXES

- Fixed missing restarts for non-blocking widgets

                 Changes in version 1.1.1 (2020-11-04)                  

FIXES

- Functions that make use of widgets do not interrupt execution
anymore if errors are thrown while producing or printing the widgets
- Optimized widget printing for importing functions
- If widgets can't be printed and verbose option is active, reports
are now displayed on console instead (needed for usage in
environments that do not have access to a browser)
- Other minor fixes (typos)
- Bug fixes: fixed a few bugs in importing and recalibration
functions
- Minor fix in import_association_file file function: added multiple
strings to be translated as NA

IMPORTANT NOTES

- Vignette building might fail due to the fact that package
"knitcitations" is temporarily unavailable through CRAN
- ISAnalytics is finally in release on bioconductor!

[iSEE](/packages/iSEE)
----

                       Changes in version 2.3.14                        

- Allow modification of font sizes for row and column names in
ComplexHeatmapPlot.
- Bugfix for assignment of annotation colors in ComplexHeatmapPlot.

                       Changes in version 2.3.13                        

- Avoid partial name matching in .getCachedCommonInfo.
- Deprecated iSEEOptions in favor of panelDefaults (for
construction-time globals) and registerAppOptions (for runtime
globals).

                       Changes in version 2.3.12                        

- Added an .allowableColorByDataChoices generic for downstream panels
to control ColorBy*Data choices.

                       Changes in version 2.3.11                        

- Cleaned up tours for Tables and the ComplexHeatmapPlot.

                       Changes in version 2.3.10                        

- Document and export the .getDotPlotColorHelp utility.
- Bugfix for the RowDotPlot color tour.

                        Changes in version 2.3.9                        

- Add a distributed tour attached to each individual UI element.
- Bugfix for ordering of selected columns in ComplexHeatmapPlot.

                        Changes in version 2.3.8                        

- Use shiny::MockShinySession$new() to simulate Shiny session
objects.

                        Changes in version 2.3.7                        

- Bugfix for missing import of geom_density_2d

                        Changes in version 2.3.6                        

- Bugfix for graceful deprecation of old parameters in various
constructors.

                        Changes in version 2.3.5                        

- Added functionality to use multiple row/column selections as a
factor on the axes, for faceting or for coloring.
- Moved selection transparency setter into the "Visual parameters"
box.
- Deprecated SelectionEffect="Color" in favor of ColorBy="Column
selection" and ColorBy="Row selection".
- Deprecated SelectionColor as the coloring for selections is
determined using colDataColorMap() instead.
- Deprecated SelectionEffect="Restrict" in favor of
ColumnSelectionRestrict and RowSelectionRestrict.
- Deprecated ColumnSelectionType and ColumnSelectionSaved (ditto for
rows) as all active/saved selections are now transmitted.

                        Changes in version 2.3.4                        

- Fix wiring of button observer to open vignette.

                        Changes in version 2.3.3                        

- Edge-case bugfix for correct cleaning of zero-row/column
SummarizedExperiments.

                        Changes in version 2.3.2                        

- Added the cleanDataset() generic to ensure all names in the
SummarizedExperiment are present and unique.

                        Changes in version 2.3.1                        

- Bugfix to the heatmap color selection for near-zero length ranges.

[iSEEu](/packages/iSEEu)
-----

                        Changes in version 1.3.5                        

- Switch to registration for storing DE Panel options, via
registerPValuePatterns and related functions.

                        Changes in version 1.3.4                        

- Support in-memory feature set collections and their statistics via
registerFeatureSetCollections.

                        Changes in version 1.3.3                        

- Redistributed documentation from panel tours to UI-specific tours.

                        Changes in version 1.3.2                        

- Tour-related patch to fix the builds for the time being.

                        Changes in version 1.3.1                        

- Added the MarkdownBoard panel to show arbitrary Markdown-formatted
content.
- Eliminate duplicates in available fields, as these break
selectizes.

[IsoformSwitchAnalyzeR](/packages/IsoformSwitchAnalyzeR)
---------------------

                Changes in version 1.13.07 (2021-05-06)                 

- Update type: minor.

- importGTF() and importRdata() was updated to handle the rare cases of
  mixed stranded and unstranded isoforms (unstanded are now
  discareded).

- addORFfromGTF() was updated to better repport if no or only small
  number of ORFs were added.

- Various maintainance updates.

                Changes in version 1.13.06 (2021-04-09)                 

- Update type: minor.

- The runtimes repported by isoformSwitchTestDEXSeq() was updated to
  also consider the number of transcripts analysed.

- analyzeORF() was updated to enable analysis with
  analyzeNovelIsoformORF() when no overlaps were found.

- switchPlot was fixed so the alphas argument now work.

- Various updates of warning, descriptions and error messages.

- extractSequence() was updated to remove the terminal stop codon if it
  is included in the annoation.

- extractSequence() was updated to produce evenly sized files when
  alsoSplitFastaFile=TRUE.

- analysORF no longer allows identification of truncated ORFs.

                Changes in version 1.13.05 (2021-01-07)                 

- Update type: Major.

- analyseORF was updated with the orfMethod
  "longest.AnnotatedWhenPossible" a hybrid between "longes" and
  "longestAnnotated". See ?analyseORF for details.

- importGTF, importRdata and analyseORF was updated to also annoate the
  source of the ORF annoations. analyzeCPAT and analyzeCPC2 was updated
  to also changes these if removeNoncodinORFs = TRUE.

- To enable better ORF analysis the addORFfromGTF() and
  analyzeNovelIsoformORF() functions were added to
  IsoformSwitchAnalyzeR. These should be used instead of analyzeORF().
  These function also annotate the source of the ORF annoations. See
  vignette for description of why these are preferable.

- analyseORF() was updated with an additional method for ORF detection:
  "longest.AnnotatedWhenPossible".

- the getCDS() function and CDSSet class was removed for the user as
  addORFfromGTF() + analyzeNovelIsoformORF() provides a better way to
  analyse ORFs.

- Downstream functions relying on ORF data now checks that all isoforms
  have been assessed for ORFs. These are extractSequence(),
  analyzeSwitchConsequences(), switchPlotTranscript() and switchPlot().

- isoformSwitchAnalysisPart1() and isoformSwitchAnalysisPart2() was
  also updated to support the new ORF annotation scheme.

- The usage of isoformSwitchAnalysisPart1() and
  isoformSwitchAnalysisPart2() was made less complex by removing many
  arguments passed to sub-functions thereby relying more on default
  arguments.

- importRdata() was updated to import the "refrence gene_ids" instead
  of StringTie gene_ids (for all annotated genes).

- the StringTie annotation rescue in importRdata() was updated to use
  "refrence gene_ids" instead of "refrence gene_names" thereby fixing
  problems with closely spaced genes, that have the same gene name,
  which was merged by StringTie.

- importGTF() now also imports ref_gene_id from StringTie gtf to enable
  the above mentioned updates to importRdata(). If not pressent it will
  duplicate gene_name instead.

- extractGeneExpression() was updated to allow easy output of gene
  annoation.

- isoformToGeneExp() was updated to use rowsum() instead of a tidyverse
  implementation as it is much faster for large datasets.

- The result of importRdata()'s estimateDifferentialGeneRange option
  now repports the condition names in accordance with the rest of
  IsoformSwitchAnalyzeR.

- Removed mentions of StringTie2 as it has been merged into StringTie.

- Documentation and vignette was updated accordingly.

                Changes in version 1.13.04 (2020-12-10)                 

- Update type: minor.

- Vignette update.

                Changes in version 1.13.03 (2020-12-08)                 

- Update type: minor.

- Vignette update.

                Changes in version 1.13.02 (2020-12-07)                 

- Update type: minor.

- Description update.

- Update of vignette with regards to running on analysis Gallaxy.

                Changes in version 1.13.01 (2020-10-29)                 

- Update type: minor.

- Version bump due to Bioconductor release.

- Fixed an error in importRdata() that could cause trouble when fixing
  StringTie annotation. Thanks to @yaccos for identifying the problem.

- Fixed an edge-case senario where the estimation of DTU in
  importRdata() caused an error.

- analyseSignalP() was updated to handle cases where no signal peptides
  were found with a warning instead of an error.

[isomiRs](/packages/isomiRs)
-------

                       Changes in version 1.19.2                        

FIX

- Add affiliation of GE

- Remove funs from dplyr to avoid future errors

- Fix column_to_names error when there are rownames

[kebabs](/packages/kebabs)
------

                       Changes in version 1.26.0                        

- release as part of Bioconductor 3.13

                       Changes in version 1.25.2                        

- removed 'register' from macro in ksort.h in order to avoid warnings
  on Mac OS

                       Changes in version 1.25.1                        

- minor fix in MismatchC.cpp

                       Changes in version 1.25.0                        

- new branch for Bioconductor 3.13 devel

[KnowSeq](/packages/KnowSeq)
-------

                 Changes in version 1.4.5 (2021-04-15)                  

- SVA batch effect correction improved

- MAD outliers detection fixed

- KnowSeq report updated to include the changes of SVA and MAD
  modifications

- Cross-Validation DEGs Extraction implemented
  Further versions

- Incorporation of RUV to batch effect methods

[limma](/packages/limma)
-----

                       Changes in version 3.48.0                        

New functionality

- 
  Explicitly setting `weights=NULL` in a call to lmFit() no
  longer over-rides the `weights` value found in `object`.
  Default settings in lmFit() changed from `ndups=` and
  `spacing=1` to `ndups=NULL` and `spacing=NULL`, although this
  doesn't change function behavior from a user point of view.

- 
  A number of improvements to duplicateCorrelation() to make the
  results more robust and to make the interface consistent with
  lmFit().  duplicateCorrelation() sets `weights` same as
  lmFit(). Setting `weights=NULL` in the function call no longer
  overwrites weights found in `object`.  duplicateCorrelation()
  now checks whether the block factor is spanned by the design
  matrix. If so, it returns intrablock correlations of zero with
  a warning. Previously this usage error was not specifically
  trapped and could lead to correlations that were or NA or close
  to 1 depending on floating point errors.
  duplicateCorrelation() now issues a simplified message when
  design is not of full rank and uses message() to do so instead
  of cat().  In terms of output, duplicateCorrelation() now
  bounds the genewise correlations away from the upper and lower
  bounds by 0.01 so that the correlation matrix will always be
  positive-definite.  There is also a fix to the value returned
  by duplicateCorrelation() when no blocks or duplicates are
  present.

- 
  New argument `fc` for treat() so that the fold-change threshold
  can optionally be specified on the fold-change scale rather
  than as a log2-fold-change.

- 
  plotMDS() no longer calls cmdscale() but instead performs the
  necessary eigenvector computations directly. Proportion of
  variance explained by each dimension is now computed and is
  optionally added to the dimension labels.  The `ndim` argument
  is now removed. All eigenvectors are now stored so that
  plotMDS.MDS does not need to recompute them when different
  dimensions are plotted.

- 
  New arguments `path` and `bgxpath` for read.idat().
  read.idat() now checks for gzipped IDAT files and, if detected,
  gives an informative error message.  read.idat() now checks for
  existence of input files before calling illuminaio read
  functions.

Other code improvements

- 
  The average log-expression column written by write.fit() now
  has column heading "AveExpr" to match the output from
  topTable().  The column was previously called "A".

Documentation

- 
  Additional documentation for the `design` argument of `lmFit`
  using the term "samples" instead of "arrays" and mentioning
  that the design matrix defaults to `object$design` when that
  component is not NULL.

- 
  duplicateCorrelation() help page revised including new code
  example.

- 
  Help page for voom() now explains that the design matrix will
  be set from the `group` factor of the DGEList object if
  available.

- 
  coolmap() help page now clarifies which heatmap.2() arguments
  are reserved and which can be included in the coolmap call.

Bug fixes

- 
  Fix typo in voom() warning when negative counts are detected.

[LoomExperiment](/packages/LoomExperiment)
--------------

                       Changes in version 1.10.0                        

BUG FIXES

- (v 1.9.1) use BiocIO rather than rtracklayer for import(),
  export(), and LoomFile() definitions.

[LRcell](/packages/LRcell)
------

                       Changes in version 0.99.7                        

- Fix typo

                       Changes in version 0.99.6                        

- Edit Vignette

                       Changes in version 0.99.5                        

- Add PBMC dataset from ExperimentHub

                       Changes in version 0.99.4                        

- Debug test-LRcell.R and fix the problem due to the data uploaded
- Add PBMC datasets information in utils.R
- Edit description

                       Changes in version 0.99.3                        

- change both dependency back to R>=4.1
- R>=3.6 generates a warning

                       Changes in version 0.99.2                        

- change LRcellTypeMarkers dependency to R>=3.6
- change LRcell dependency to R>=3.6

                       Changes in version 0.99.1                        

- remove the LRcell.Rproj
- change the .gitignore file

                       Changes in version 0.99.0                        

- version 0.99.0 released
- Submitted to Bioconductor

[Maaslin2](/packages/Maaslin2)
--------

                        Changes in version 1.5.1                        

- Update log from base 10 to base 2.

- ZICP is now deprecated
  (https://cran.r-project.org/web/packages/cplm/NEWS)

- SLM is removed in favor of a future R2 functionality for all models

- Fitted values are returned along with residuals

- Extracted random effects are also returned

[maftools](/packages/maftools)
--------

                        Changes in version 2.8.0                        

NEW FUNCTIONS

- cancerhotspots Genotype known cancer hotspots from the tumor BAM
file
- bamreadcounts extract nucleotide counts for targeted variants from
the BAM file.
- maftools now natively loads TCGA cohorts. tcgaAvailable and
tcgaLoad
will display and load the desired cohorts.
- Added MAF constructor function
- Added maf2mae for converting MAF to MultiAssayExperiment class
objects Issue: 640 293 Discussion: 285
- Added plotProtein and mafbarplot

ENHANCEMENTS

- Added protein domains for the gene ALMS1. Issue: 705
- Added titv_col argumtn to oncoplot. Issue: 702
- Added protein domains for the gene FAM205A. Issue: 701
- oncoplot can now summarize variant_classifications similar to
cBioPortal style. Issue: 686
- Added pathway support for mafCompare() or clinicalEnrichment().
Issue: 681
- Added default title for side and topbar plots to oncoplot. Issue:
682
- Added annotationOrder argument to coOncoplot. Issue: 676
- Added plot argument to survGroup. Thank you OmarElAshkar PR: 674
- Added rmFlags argument to read.maf. Issue: 668
- Added path_order argument to oncoplot for custom ordering of
pathways on oncoplot.
- Added geneMar argument to coBarplot. Issue: 260

BUG FIXES

- coOncoplot not allowing more than one additional feature. Issue:
675

[marr](/packages/marr)
----

                        Changes in version 1.1.2                         

- added feature-label output on May 10, 2021

- added feature-label output on April 27, 2021

- added cutoff value on March 06, 2021

                        Changes in version 1.1.1                        

- hard-coding of alpha fixed on December 21, 2020.

[MatrixGenerics](/packages/MatrixGenerics)
--------------

                        Changes in version 1.2.1                        

- Sync API with matrixStats v0.58.0.

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                Changes in version 0.99.12 (2021-05-18)                 

- replace xlsx by openxlsx

                Changes in version 0.99.11 (2021-05-10)                 

- rename function normalize to normalizeAssay

- rename function transform to transformAssay

- rename function batch to batchCorrectionAssay

- rename function impute to imputeAssay

                Changes in version 0.99.10 (2021-05-06)                 

- bump version to trigger building

                 Changes in version 0.99.9 (2021-04-29)                 

- add hexbin in Suggests

- fix bug in MAplot that plot is displayed properly

                 Changes in version 0.99.8 (2021-04-28)                 

- set required version for S4Vectors to >= 0.29.15

                 Changes in version 0.99.7 (2021-04-28)                 

- add version number of dependencies in Description file

                 Changes in version 0.99.6 (2021-04-27)                 

- add MatrixQCvis to Watched Tags on the Bioconductor support site

                 Changes in version 0.99.5 (2021-04-27)                 

- reduce file size of vignette by using partial_bundle for driftPlot

                 Changes in version 0.99.4 (2021-04-26)                 

- reduce package dependencies
  - remove magick
  - use stats::cmdscale instead of ape::pcoa
  - remove MsCoreUtils
  - remove preprocessCore
  - remove Matrix

- add explained variance for PCoA

- add se argument in create_boxplot that allows for ordering the
  samples

- use ggplotly for driftPlot

- allow flexible addition of samples in MA-plot based on a supplied
  character
  vector of sample names

- return SummarizedExperiment when exiting the shiny application

- add function maxQuant that allows for creation of
  SummarizedExperiment
  objects from maxQuant output (.xlsx files)

                 Changes in version 0.99.3 (2021-03-18)                 

- reduce file size of vignette by using partial_bundle for plotly
  figures

                 Changes in version 0.99.2 (2021-03-18)                 

- reduce resolution of images in vignette to reduce file size

                 Changes in version 0.99.1 (2021-03-17)                 

- reduce file size of vignette

                 Changes in version 0.99.0 (2021-03-12)                 

- shinyQC including visualizations/functionality for
  - histogram of sample types,
  - information on number of missing/measured values
  - information on (intersecting, disjoint) sets for missing/measured
  values
  - barplot and violin plot for (count/intensity) values
  - visualization to detect drifts/trends in (count/intensity) values
  - coefficients of variation for samples,
  - mean-sd plots,
  - MA plots,
  - empirical cumulative distribution function,
  - visualizations of distances between samples,
  - intensities of features and coefficients of variation of features,
  - dimension reduction plots (PCA, PCoA, NMDS, tSNE, UMAP)
  - differential expression

- write functions for data manipulation and plots

- write tests for these functions

- create UI and server modules for shinyQC

- write tests for UI and server modules

- load different UI elements depending on the type of data (if the data
  contains missing values or is complete)

- load different UI if the SummarizedExperiment is loaded on start of
  shinyQC
  or not

[matter](/packages/matter)
------

                 Changes in version 1.17.1 (2020-11-27)                 

BUG FIXES

- Fix 'apply()' signatures for R 4.1

[megadepth](/packages/megadepth)
---------

                        Changes in version 1.1.5                        

NEW FEATURES

- Added process_junction_table (an R-implementation of
https://github.com/ChristopherWilks/megadepth#junctions) which
parses the output of read_junction_table into a STAR-compatible
format.

[memes](/packages/memes)
-----

                       Changes in version 0.99.11                       

- Add support for STREME with runStreme(). STREME will supercede
DREME
in a future MEME Suite release.

                       Changes in version 0.99.10                       

- Fixed a bug where paths weren't correctly expanded when used as
database entry under certain conditions

                       Changes in version 0.99.8                        

- Removed inline r call in integrative_analysis vignette to fix issue
on bioc build machine

                       Changes in version 0.99.7                        

- Version bump to force pkg rebuild

                       Changes in version 0.99.6                        

- added list S3 method for plot_sequence_heatmap so now named lists
of
sequences can be passed natively to this function.
- updated ChIP-seq vignette to demonstrate this

                       Changes in version 0.99.5                        

- added plot_sequence_heatmap for making heatmaps of sequence lists
- Added significantly more explanation to the ChIP-seq vignette
- renamed ame_plot_heatmap -> plot_ame_heatmap for consistency

                        Changes in version 0.1.2                        

- runFimo() skip_matched_sequence default is now FALSE. Set this to
TRUE if fimo takes a long time to run, then use add_sequence() to
add it back if needed.
- runTomTom() dist default is now ed (changed from pearson).

                        Changes in version 0.1.0                        

- Removed as_universalmotif_df(), as_universalmotif(), and
update_motifs().
- These functions are replaced by universalmotif::to_df(),
universalmotif::to_list(), and universalmotif::update_motifs()
- runDreme and runTomTom results are now returned in
universalmotif_df
format (behaves just like a data.frame)
- The motif column of universalmotif_df objects can no longer be
called directly in universalmotif operations like
view_motifs(df$motif). Use to_list() for this behavior instead.
- To support this change, the pvalue, evalue, and qvalue columns
are renamed pval, eval, and qval. The same is true for tomtom
output columns match_pvalue -> match_pval, best_match_pvalue ->
best_match_pval, etc.
- Updated example datasets to use unviversalmotif_df type
- ame_plot_heatmap ranking issue is resolved, plots now sort
correctly
- Added remove_duplicate_motifs and has_duplicate_motifs for
detecting
and removing duplicated matrices in a universalmotif object or
data.frame
- Overhauled the Tidying Motifs vignette for more extensive EDA and a
demo of deduplication
- Updated the flyFactorSurvey_cleaned.meme example database to
reflect new changes to the vignette

[metabCombiner](/packages/metabCombiner)
-------------

                        Changes in version 1.1.3                        

Package additions

- metabCombine(): main package workflow wrapper function

- Parameter list functions for loading defaults of main workflow
  methods

Changes to labelRows

- "conflict" argument replaced with "delta", with default value (0.2)

- default value of "minScore" argument increased to 0.5

Changes to calcScores

- default argument values for A (75), B (10), C (0.25)

Changes to fit_loess

- new argument "control" for controlling loess fit (list argument)
  - eliminated "iterLoess" parameter

Changes to selectAnchors

- default for "tolrtq" argument changed from 0.5 to 0.3

                 Changes in version 1.1.2 (2020-12-28)                  

Changes to fit_gam()/ fit_loess

- new outlier detection method based on boxplot / IQR added

- new argument: outlier, which accepts "MAD" or "boxplot" as a value

- altered argument names: "frac" to "prop", "ratio" to "coef"

- documentation and minor code changes to main and supporting functions

Changes to plot_fit

- new options for showing, hiding, or highlighting (with a legend)
  outliers

- new arguments: outlier, which accepts "show" / "s", "remove" / "r",
  or
  "highlight" / "h" as arguments ; ocol, outlier point color if outlier
  argument set to "highlight" / "h"
  - remove.outliers argument removed

Other changes

- new test case for fit_gam()

                 Changes in version 1.1.1 (2020-12-02)                  

Bug Fixes

- combinedTable check for missing group values (Issue #7)

- calcScores / evaluateParams groups bug (Issue #8)

- Warning for column names with bracket characters "{ ( [ ] ) }" (Issue
  #9)

- QCol bug (Issue #10)

New Functionality

- getExtra method added to metabCombiner objects (Issue #11)

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                        Changes in version 0.99                         

MetaboCoreUtils 0.99.1

- Add [M+H-2(H2O)]+ adduct definition.

MetaboCoreUtils 0.99.0

- Add package vignette and prepare for Bioconductor submission.

[metabolomicsWorkbenchR](/packages/metabolomicsWorkbenchR)
----------------------

                        Changes in version 1.1.1                        

- bug fixes

- added some helper functions

[metagene2](/packages/metagene2)
---------

                 Changes in version 1.7.2 (2020-12-17)                  

- Removed support for NCIS normalization, which relied on a now
  deprecated package.

[metapod](/packages/metapod)
-------

                        Changes in version 1.0.0                        

- 
  New package metapod, for meta-analyses of p-values from
  differential tests.

[metaseqR2](/packages/metaseqR2)
---------

                 Changes in version 1.3.13 (2021-02-20)                 

NEW FEATURES

- Ability to use a different sample list with different condition names
  but
  same sample set (or subset) with a previous gene model with same
  samples.

BUG FIXES

- None.


[methylGSA](/packages/methylGSA)
---------

                        Changes in version 1.9.1                        

- Enhancement in package vignette

[methylInheritance](/packages/methylInheritance)
-----------------

                       Changes in version 1.15.1                        

BUG FIXES AND IMPROVEMENTS

- Fix bug causing problem when using loadConvergenceData() with a
  different directory for the observation file and the permutation
  files.

[methylKit](/packages/methylKit)
---------

                       Changes in version 1.17.5                        

IMPROVEMENTS AND BUG FIXES

- fread.gzipped: skip header rows in tabix file to fix
  https://github.com/al2na/methylKit/issues/226

                       Changes in version 1.17.4                        

IMPROVEMENTS AND BUG FIXES

- .setMethylDBNames: correct possible methylDBclass from "methylDB"
  to "methylRawDB"

- fread.gzipped: disable skipping of decompression in fread.gzipped
  - can cause serious issues as investigated in
  https://groups.google.com/g/methylkit_discussion/c/UruFjvX89B4/m/vV2Qnd8NEAAJ
  and explained in https://github.com/al2na/methylKit/issues/222

                       Changes in version 1.17.3                        

IMPROVEMENTS AND BUG FIXES

- export methylRawListDB and methylRawList constructors

                       Changes in version 1.17.2                        

IMPROVEMENTS AND BUG FIXES

- update vignette:
  - add note about use of system.file function
  - reduce number of rows, when showing row ordering
  - add faq section about merging methylRaw into methylRawList update
  vignette

                       Changes in version 1.17.1                        

IMPROVEMENTS AND BUG FIXES

- readmethylDB: check if file exists before trying to read

- loading tabix files:
  - improve reading tabix files with header
  - update tests for loading tabix files

[MetNet](/packages/MetNet)
------

                 Changes in version 1.9.5 (2021-05-04)                  

- add functionality to adjust for multiple testing in correlation

- return symmetric matrices when returning ppm ranges in structural

                 Changes in version 1.9.4 (2021-04-30)                  

- add font in mz_vis to mono

- split the example on how to use filter in mz_summary from the
  visualisation

                 Changes in version 1.9.3 (2021-04-28)                  

- introduce AdjacencyMatrix S4 class, derived from
  SummarizedExperiment, to
  store the adjacency matrices. AdjacencyMatrix objects can be of type
  structural, statistical, and combine

- adjust the documentation and tests for AdjacencyMatrix objects

- add the functions mz_summary and mz_vis (contribution of Liesa
  Salzer)

                 Changes in version 1.9.2 (2021-03-24)                  

- add section on structual matrix generation for directed=FALSE

                 Changes in version 1.9.1 (2021-02-20)                  

- fix typos in the vignette

[mia](/packages/mia)
---

                 Changes in version 0.99.0 (2021-03-19)                 

- Submitted to Bioconductor

[miaViz](/packages/miaViz)
------

                 Changes in version 0.99.0 (2021-03-19)                 

- Submitted to Bioconductor

[microbiome](/packages/microbiome)
----------

                 Changes in version 2.1.2 (2020-07-01)                  

- Core heatmap labeling improved

- aggregate_top_taxa deprecated

- bimodality and potential_analysis functions fixed

                 Changes in version 2.1.1 (2020-04-06)                  

- Added overlap function

[microbiomeExplorer](/packages/microbiomeExplorer)
------------------

                        Changes in version 1.0.3                        

- added option to split a taxonomy column (as obtained via qiime)
  within the application

                        Changes in version 1.0.2                        

- fixed bug in correlation analysis on phenotypes with NaN values

- prevented application crash on loading ill-formatted feature data

                        Changes in version 1.0.1                        

- adjusted to be compatible with shinyjs 2.0.0

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.3.11                        

- fill ggclust bug to map color and shape. (2021-05-12, Wed)

                        Changes in version 1.3.9                        

- more layouts for ggdiffclade. (2021-04-16, Fri)

- remove unclassified, ambiguous taxonomy names. (2021-03-30, Tue)
- https://github.com/YuLab-SMU/MicrobiotaProcess/issues/23

                        Changes in version 1.3.8                        

- rename files of code. (2021-03-23, Tue)
- add aliases for ggbartaxa and ggdiffbartaxa. (2021-03-23, Tue)

                        Changes in version 1.3.7                        

- optimized import for installation. (2021-03-15, Mon)
- add sampleLevels in ggbartax to adjust the order of axis label.
(2021-03-12, Fri)
- update import_qiime2 to avoid error when only feature table is
provided. (2021-02-26, Fri)

                        Changes in version 1.3.6                        

- convert svg dev to pdf dev. (2021-02-04, Thu)

                        Changes in version 1.3.5                        

- fix an error for example of ggrarecurve. (2021-01-07, Thu)
factorNames="Group" to factorNames="group"

                        Changes in version 1.3.4                        

- fix a bug for numeric sample name. (2020-11-26, Thu)
- geom_tiplab also support circular layout, so remove geom_tiplab2.
(2020-11-26, Thu)

                        Changes in version 1.3.3                        

- add as.treedata for taxonomyTable class. (2020-11-23, Mon)

                        Changes in version 1.3.2                        

- ggrarecurve can be set color with variable of group for each
samples. (2020-11-11, Tue)
- using shadow=FALSE and providing factorNames
- https://github.com/YuLab-SMU/MicrobiotaProcess/issues/21
- add get_rarecurve to avoid repeated calculation when displaying
rare
curve. (2020-11-17, Tue)
- rareres <- get_rarecurve(obj, chunks=400) p <-
ggrarecurve(rareres)

                        Changes in version 1.3.1                        

- ggordpoint add showsample to show the labels of sample.
(2020-10-29,
Thu)
- the point of ggordpoint use the points of ggstar. (2020-10-30, Fri)
- to obtain the dynamic arguments of diff_analysis, the call was
changed to someparams. someparams contained the arguments used in
other functions. (2020-11-09, Mon)
- https://github.com/YuLab-SMU/MicrobiotaProcess/issues/20

[midasHLA](/packages/midasHLA)
--------


                       Changes in version 0.99.15                       

- fixes bug causing LRT to report wrong number of tested residues.

                       Changes in version 0.99.13                       

- changes version number to comply with Bioconductor package
submission guidelines.
- rename package to midasHLA; midas is already in use.


[miloR](/packages/miloR)
-----

                 Changes in version 0.99.1 (2021-03-13)                 

- Fix model normalisation bug - now using TMM normalisation by
default. Log(M_s) offset can be used by passing norm.method="logMS"
to testNhoods.

                 Changes in version 0.99.0 (2021-03-04)                 

- Submitted to Bioconductor

[miQC](/packages/miQC)
----

                 Changes in version 0.99.9 (2021-03-24)                 

- Added new function, plotMetrics

- mixtureModel now throws a warning if flexmix has not identified two
  distributions

- Vignette includes instructions for assessing whether to use miQC on a
  dataset

                 Changes in version 0.1.0 (2021-02-10)                  

- Submitted to Bioconductor

[miRSM](/packages/miRSM)
-----

                     Changes in version 1.9.4-1.9.5                     

- Add citation information <2021-04-15, Thus>

                        Changes in version 1.9.3                        

- Add sponge module (SM) method <2021-02-02, Thue>

                     Changes in version 1.9.1-1.9.2                     

- Update miRSM.R <2020-11-26, Thus>

[miRspongeR](/packages/miRspongeR)
----------

                    Changes in version 1.17.1-1.17.2                    

- Update src <2020-11-09, Mon>.

[mirTarRnaSeq](/packages/mirTarRnaSeq)
------------

                Changes in version 0.99.11 (2021-05-14)                 

- Changed News file

                       Changes in version 0.99.10                       

- shrunk repo

                       Changes in version 0.99.9                        

- updated with suggested notes

                       Changes in version 0.99.8                        

- fixed check error

                       Changes in version 0.99.7                        

- dealt with download issue

                       Changes in version 0.99.6                        

- dealt with vignette issue

                       Changes in version 0.99.5                        

- fixed documentation warnings for datasets and warnings for R CMD file
  size

                       Changes in version 0.99.4                        

- fixed documentation warnings for datasets

                       Changes in version 0.99.3                        

- fixed documentation warnings for code

                       Changes in version 0.99.2                        

- initiated


[mistyR](/packages/mistyR)
------

                       Changes in version 0.99.11                       

- Fixed a bug in Nystrom approximation for creating paraview.
- Added a suite of tests with high coverage.
- Cleaner cache control.

                       Changes in version 0.99.9                        

- Revisions requested by Bioconductor
- Removed magrittr from dependencies. Reexported pipe operator.
Removed from examples.
- Caching turned off by default.
- Added parameter for appending performance and coefficient files
in run_misty.
- Internal functions in views.R are now explicit.
- Removed alternatives and additional examples for function use
from documentation.
- Removed redundant messages, escalated to warnings where
required.
- run_misty cleans up empty cache directories.
- Vignette compatibility with the new release of
SpatialExperiment.
- Remove Seurat vignette from package due to missing hdf5r binary
for R 4.1 on CRAN for MacOS.
- Other changes
- README.md figures moved to the cloud.
- All passed paths and cache location are normalized.
- Changes in the mistyR vignette to reflect changes to insilico
evaluation from the paper.

                       Changes in version 0.99.0                        

- Version with vignettes ready to submit to Bioconductor.

                        Changes in version 0.1.0                        

- Initial beta release of mistyR (named as MISTy) with function
documentation.

[MLInterfaces](/packages/MLInterfaces)
------------

                       Changes in version 1.71.1                        

- rda is defunct and all references are removed

- bug in hclustWidget forbade use of more features than samples -- this
  is fixed

[MLP](/packages/MLP)
---

                       Changes in version 1.39.1                        

- use Imports rather than Depends

- use Authors@R

- KEGG.db (not available from BioC >= 3.13) -> KEGGREST

- getGeneSets: only return descriptions for selected pathways

- update example gene set

[ModCon](/packages/ModCon)
------

                 Changes in version 0.99.0 (2018-05-15)                 

- Submitted to Bioconductor

[MOFA2](/packages/MOFA2)
-----

                        Changes in version 1.1.9                        

- integrated MEFISTO into MOFA2

- Improve interoperability with Seurat and SingleCellExperiment

- Improve memory usage and training time by optionally replacing
  float64 arrays with float32 arrays

- mofapy2 has been updated to version 0.6.0 and now it has its own
  repository

                        Changes in version 1.1.4                        

- integrated MEFISTO into MOFA2

- Improving Python interface with basilisk

- Sample metadata can be incorporated to the MOFA object before and
  after training

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.35.2                        

- Add rmarkdown as suggest package.

[mpra](/packages/mpra)
----

                       Changes in version 1.12.1                        

- Fix argument checking in MPRASet construction to allow users to
  not have to specify barcode or eseq.

- Fix bug with ordering of eids in log ratio object with
  aggregate="none" option

[MQmetrics](/packages/MQmetrics)
---------

                       Changes in version 0.99.4                        

- The function generateReport() now takes the parameter
name_output_file to name the output pdf file.
- Added information regarding peptides, and peptides/protein ratio in
report tables.
- Added the function PlotProteinPeptideRatio() to visualize a
comparison between the proteins identified and the ratio
Peptide/Proteins among Experiments.

                       Changes in version 0.99.3                        

- Functions now use a general input resulting from make_MQCombined().
- Shortened lines in the code.
- Usage of seq_len rather than 1: ...
- Usage of ultiples of 4 spaces for line indents.
- Removed LazyData: TRUE.

                       Changes in version 0.99.2                        

- Added a NEWS.md file to track changes to the package.

[msa](/packages/msa)
---

                       Changes in version 1.24.0                        

- release as part of Bioconductor 3.13

                       Changes in version 1.23.1                        

- updated texshade.sty to newest version

                       Changes in version 1.23.0                        

- new branch for Bioconductor 3.13 devel

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                        Changes in version 0.99                         

Changes in 0.99.4

- Fix an issue in which selectSpectraVariables,MsBackendMassbankSql
would fail if subsetted to a single variable/column in @localData
issue #29

Changes in 0.99.3

- Drop names on the peaksData.

Changes in 0.99.2

- Minor updates and changes as requested during package review
process.

Changes in 0.99.1

- Directly call package internal functions using ::: if used within
bplapply to avoid errors on Windows.

Changes in 0.99.0

- Prepare for Bioconductor submission.

[MsBackendMgf](/packages/MsBackendMgf)
------------

                        Changes in version 0.99                         

Changes in 0.99.3

- Export readMgf function.

Changes in 0.99.2

- Update installation description in the vignette.

Changes in 0.99.1

- Directly call internal function from MsBackendMgf to avoid parallel
processing error (function not found) on Windows.

               
[MsCoreUtils](/packages/MsCoreUtils)
-----------

                         Changes in version 1.3                         

Changes in 1.3.3

- Add join_gnps and gnps to allow calculation of GNPS spectra
similarity scores.

Changes in 1.3.2

- Add rt2numeric(), rt2character() and formatRt().
- New impute_fun() function for user-provide imputation function.

Changes in 1.3.1

- Add Josep Badia Aparicio as a contributor

Changes in 1.3.0

- New Bioc devel versin

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.17                         

Changes in 2.17.7

- Bump to force propagation on Bioc

Changes in 2.17.6

- Import ProcessingStep from ProtGenerics.

Changes in 2.17.5

- Import quantify generic from ProtGenerics.
- Fix defaults in calculateFragments man page (see issue #537).

Changes in 2.17.4

- Re-implement impute,MSnSet-method with
MsCoreUtils::impute_matrix().
- Re-implement normalise,MSnSet-method with
MsCoreUtils::normalize_matrix().
- Illustrate ReporterIons class using TMT16.
- Re-implement formatRt() using MsCoreUtils::formatRt().
- Use MsCoreUtils::robustSummary() (deleted code in this package).
- Use MsCoreUtils::medianPolish().

Changes in 2.17.3

- Add normalize,Chromatogram method.
- Add filterIntensity,Chromatogram method.
- Add align,Chromatogram method.
- Use two digits when printing retention time seconds

Changes in 2.17.2

- Add TMT16 reporter ions (contribute by Miguel Cosenza-Contreras)

Changes in 2.17.1

- Use ProtGenerics::calculateFragments generic
- Document differences between spectrumId (spectrumID),
acquitisionNum
and scanIndex.

Changes in 2.17.0

- New Bioc devel version

[MSPrep](/packages/MSPrep)
------

                 Changes in version 1.1.1 (2020-03-27)                  

- Adds random forest imputation
- Fixes param validation bug

[msPurity](/packages/msPurity)
--------

                       Changes in version 1.17.2                        

- Update dev to match bug fixes in master

                       Changes in version 1.17.1                        

- Update dev to match bug fixes in master

                       Changes in version 1.16.2                        

- Author list updated

                       Changes in version 1.16.1                        

- Update of createDatabase to record all intra average spectra in
  database

- Add license and copyright info to code

- Add github workflow CI (and subsequent formatting updates to pass
  tests)

[msqrob2](/packages/msqrob2)
-------

                       Changes in version 0.99.6                        

- Update authors in Description file

                       Changes in version 0.99.5                        

- Fix standard errors on the model parameter estimates by msqrobLmer
when using doQR = TRUE

                       Changes in version 0.99.4                        

- Minor update vignette. Replace eval=FALSE in one R chunk so that
the
code is evaluated.

                       Changes in version 0.99.3                        

- Resolved notes on signalers
- Resolved examples with lines wider than 100 characters
- Avoiding sapply and 1:...

                       Changes in version 0.99.2                        

- Added a NEWS.md file to track changes to the package.
- Changed dependency from R (>= 4.0) to R (>= 4.1)
- Updated citation file

Changes in version 0.99.1

- Submission to Bioc

[MSstatsConvert](/packages/MSstatsConvert)
--------------

                        Changes in version 1.1.2                        

- Updated Spectronaut converter to allow annotation in input file.

                        Changes in version 1.1.1                        

- Minor bug fixes
- Updated MSstatsLogsSettings function to allow differents settings
for different packages

[MSstatsTMT](/packages/MSstatsTMT)
----------

                 Changes in version 2.0.0 (2021-05-14)                  

- Refactor the pacakge to make it modulized

                 Changes in version 1.8.2 (2020-12-17)                  

- Update progress bar

- Update groupComparisonTMT to avoid reusing the local function
  copied from lmer pacakge

                 Changes in version 1.8.1 (2020-12-10)                  

- Add citation of the MSstatsTMT paper

- Fix the bug in groupComparisonTMT() due to the update of
  dependent pacakge

- Fix the bug in MedianPolish summarization

- proteinSummarization(): replace the zero values with NA before
  and after peptide normalization

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.18.0                        

New features

- saveHDF5MultiAssayExperiment allows users to save data from most
classes (excluding RaggedExperiment) into a single H5 file (ctb @hpages)
- Support for maftools conversion has been added as
MultiAssayExperimentToMAF (ctb @PoisonAlien)
- renameColname and renamePrimary provide renaming facilities for
column names in experiments and rownames in the colData,
respectively
- Users can now rename some or all the column names in experiments
using colnames(mae) <- value
- When replacing colData or experiments (including [[<-), new
rownames
and colnames (respectively) are checked against existing values and
an error is given when none match
- Using List objects to replace the data in the ExperimentList is now
supported
- splitAssay allows users to separate / split columns across assays
- makeHitList is a facilitator function to create the logical lists
that are required as input to splitAssay
- drop argument when subsetting a MultiAssayExperiment is now FALSE
by
default

Bug fixes and minor improvements

- Updated the constructor function to auto-populate rownames in
colData when it is missing (@LTLA, #287)
- The metadata now includes names of dropped experiments
- Updated validity checks to support array-like classes
- Dropped experiments are tracked in the metadata

[multiGSEA](/packages/multiGSEA)
---------

                 Changes in version 1.1.7 (2021-05-10)                  

- Remove Biocarta database from vignette - it's no longer
  supported by graphite R package.

- List all available databases in error message.

                 Changes in version 1.1.5 (2021-04-19)                  

- Typos fixed.

- Add information and reference about metrics to calculate the
  feature ranks.

                 Changes in version 1.1.4 (2020-12-09)                  

- Add the correct citation of the BMC Bioinformatics article.

                 Changes in version 1.1.3 (2020-11-19)                  

- Prevent warning due to import of two different select()
  functions.

                 Changes in version 1.1.2 (2020-11-07)                  

- Speed-up the metabolite ID mapping between different ID formats

                 Changes in version 1.1.1 (2020-11-05)                  

- Bug fix. Prevent duplicated pathway title to cause an error.

- Include new function to enumerate duplicted pathway titles with
  trailing numbers.

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

[multiSight](/packages/multiSight)
----------

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[MungeSumstats](/packages/MungeSumstats)
-------------

                        Changes in version 1.0.0                        

New Features

- MungeSumstats released to Bioconductor

[muscat](/packages/muscat)
------

                        Changes in version 1.5.2                        

- added edgeR::calcNormFactors() step in prepSim()

- added argument 'dd' to simData() specifying
  whether or not to simulate 2 groups

- prepSim() and simData() now support simulation of "singular" design
  (no samples, no clusters), as well as only samples/clusters

- simData() defaults to simulating as many samples as available
  in order to avoid re-use (duplication) of reference samples

                        Changes in version 1.5.1                        

- significant speed-up of aggregateData() by replacing usage
  of rowX() over a list with scuttle::summarizeAssayByGroup()

- added options use "prop.detected" and "num.detected"
  as summary statistic (argument 'fun') in aggregateData()

- added parallelization support in aggregateData() and pbDS() through
  argument BBPARAM
  (passed to scater::sumCountsAcrossCells() and BiocParallel::bplapply,
  respectively)

- aggregateData() now stores the number of cells that went into
  aggregation under
  int_colData(.)$n_cells (vs. metadata(.)$n_cells) to support automated
  subsetting

- replaced argument n_threads with BPPARAM throughout all
  parallelizable functions (aggregateData(), pbDS(), mmDS())

- bug fix in prepSim(): the function previously failed when
  cluster/sample/group_id cell metadata columns were non-factors

- bug fix in resDS(): cpm = TRUE previously didn't handle
  missing cluster-sample combinations correctly

[musicatk](/packages/musicatk)
--------

                 Changes in version 1.1.2 (2021-02-03)                  

RELEASE

- Added new functionality for analysis including clustering sample
exposures using a variety of methods (e.g. k-means), differential
exposure using a variety of methods (e.g., Wilcoxon rank-sum,
Kruskal Wallis, GLM), and heatmap creation. Fixed a bug causing
incorrect counting of mutations during creation of indel tables.

[MutationalPatterns](/packages/MutationalPatterns)
------------------

                 Changes in version 3.1.4

* Plot lesion segregation now has background strips.
* Plot_rainfall can now remove more chromosome prefixes from the labels. 
chromosome(_), group(_) and chrom can now be removed in addition to chr. All are case insensitive.
* Improved the description in the vignette of choosing the cut-off for strict refitting.
* Fixed spelling mistakes.

                 Changes in version 3.1.2

* Plot_lesion_segregation has been improved. It can now plot multiple samples at the same time.
Users can also specify which chromosomes they want to plot. The plot now also contains colour and
the ratio of the mutations on the chromosomal strands is 
visualised by a horizontal line per chromosome.

                 Changes in version 3.1.2

* The plot_correlation_bootstrap output now shows less dominant colours for the
absent signatures.
* The rename_nmf_signatures function now adds the "-like" suffix to renamed signatures,
by default.
* The plot_river function now shows the number of mutations per sample.
* The read_vcfs_as_granges function now has a "predefined_dbs_mbs" argument, for when
dbs and mbs variants are already defined in the vcfs. This prevents merging of snvs.
* The get_indel_context function now looks for repeats both to the left and to the right of the
indel. Previously it only looked to the right of the indel.

 
 
[mzR](/packages/mzR)
---

                       Changes in version 2.25.5                        

- Fix compile error on clang-11 reported (and fixed!) by Kurt Hornik,
  closes #244

                       Changes in version 2.25.4                        

- Add dependency "rmarkdown" to "suggests:"

                       Changes in version 2.25.3                        

- Ensure `header` for CDF returns columns with correct data type.

                       Changes in version 2.25.2                        

- Fix issue #238: ensure `header` call returns the same columns for all
  backends.

                       Changes in version 2.25.1                        

- Bump version to trigger new build using latest Rcpp

                       Changes in version 2.25.0                        

- New Bioc devel version

[NADfinder](/packages/NADfinder)
---------

                       Changes in version 1.15.2                        

- add rmarkdown as suggest package.

                       Changes in version 1.15.1                        

- fix the missing link for rtracklayer::export.

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 1.1.4                        

- Added palette argument to aggregate plots
- Added exons_to_genes() function to convert exon annotation to gene
annotation
- Added plot_granges_heatmap() function to use GRanges for plotting
heatmaps

                        Changes in version 1.1.3                        

- Fixed group handling for list region input in plot_agg_regions()
- Fixed unused window size argument in plot_region_heatmap()
- Fixed error when reads overlap in name and position for internal
function StatLM()

                        Changes in version 1.1.2                        

- Changed example dataset exon annotations from all genomic exons to
just those contained in data.
- Fixed methylation heatmap to no longer be hard coded for Peg3.
- Added plot_region_heatmap() as analogue to plot_region().
- Fixed plot_agg_regions_sample_grouped() to use group column of
NanoMethViz::samples(x) rather than haplotype.
- Added unit tests.

                        Changes in version 1.1.1                        

- Added methylation heatmap via plot_gene_heatmap().
- Fixed gene_anno() in plot_gene() for argument so FALSE actually
turns off gene annotation.
- Added warning for cpp11 versions <0.2.5 which may cause memory
crashes when trying to import methylation data.
- Added cpp11 version dependency to address tidyverse/readr#1145.
- Added query methylation by gene using query_methy_gene().

[NanoStringNCTools](/packages/NanoStringNCTools)
-----------------

                 Changes in version 0.99.5 (2020-12-21)                 

- Fix the issues from the build report

                 Changes in version 0.99.4 (2020-12-18)                 

- Fix the issues from the build report

                 Changes in version 0.99.3 (2020-12-18)                 

- Fix the issues from the build report

                 Changes in version 0.99.2 (2020-12-18)                 

- Fix the issues from the build report

                 Changes in version 0.99.1 (2020-12-17)                 

- Fix the issues from the build report

                 Changes in version 0.99.0 (2020-12-15)                 

- Submitted to Bioconductor

[ncRNAtools](/packages/ncRNAtools)
----------

                 Changes in version 1.1.1 (2021-05-04)                  

- Fixed a bug causing endless queries in Ubuntu 20.04 when a firewall
  was present

[nearBynding](/packages/nearBynding)
-----------

                        Changes in version 1.1.4                        

Changes

- remove variable assignment bugs from symmetryContext

[netboost](/packages/netboost)
--------

                 Changes in version 1.99.4 (2021-04-15)                 

- Removal of MCUPGMA-dependencies for smaller networks.

                 Changes in version 1.99.0 (2021-03-08)                 

- Fully rank based extension
  (netboost(...,robust_PCs=TRUE,filter_method="spearman",method="spearman")).

- Full support of the non-parametric version.

[ngsReports](/packages/ngsReports)
----------

                        Changes in version 1.7.3                        

- Simplifed version of plotFastqcPCA. Now groups are an optional
  factor. No clustering is performed

                        Changes in version 1.7.2                        

- Deprecated runFastQC

                        Changes in version 1.7.1                        

- Added macs2 callpeak logs to importNgsLogs

                        Changes in version 1.6.1                        

- Added asPercent to plotAlignmentSummary

- Added the ability to assign new values via fqName<-

[NormalyzerDE](/packages/NormalyzerDE)
------------

                        Changes in version 1.9.1                        

- Version sync with Bioconductor

[NuPoP](/packages/NuPoP)
-----

                       Changes in version 1.99.0                        

- Included predNuPoP_chem function which predicts the nucleosome
  positioning based on profiles trained based on chemical maps

- Vignette file has been created with R markdown

- Added NEWS file to document version bump

[OmaDB](/packages/OmaDB)
-----

                 Changes in version 2.7.1 (2020-11-11)                  

- fix problem after OMA API version change

[OmnipathR](/packages/OmnipathR)
---------

                       Changes in version 2.99.19                       

- NicheNet pipeline works
- Fixed an error which resulted value 1 in the n_references columns
even for records without references

                       Changes in version 2.99.17                       

- Improved quality filtering of intercell networks

                       Changes in version 2.99.16                       

- Improved API for translate_ids
- Quality filtering of intercell networks

                       Changes in version 2.99.13                       

- Functions to access KEGG Pathways
- More robust access to UniProt (in case of network failures)

                       Changes in version 2.99.11                       

- Improved downloader backends

                       Changes in version 2.99.8                        

- OBO parser
- Gene Ontology access, functions to work with ontology trees
- Database manager

                       Changes in version 2.99.7                        

- New interactions query parameter: loops
- Fixed many caching bugs

                       Changes in version 2.99.6                        

- All downloaders attempt 3 times
- New resources: Human Phenotype Ontology and Gene Ontology
annotations

                 Changes in version 2.99.0 (2021-03-08)                 

- Many new resources besides OmniPath: BioPlex, ConsensusPathDB,
EVEX,
Guide to Pharmacology (IUPHAR/BPS), Harmonizome, HTRIdb, InWeb
InBioMap, Pathway Commons, Ramilowski et al. 2015, RegNetwork,
ReMap, TF census, TRRUST and Vinayagam et al. 2011
- Methods for converting network elements to Bio Model Analyzer (BMA)
motifs
- NicheNet workflow
- Some igraph related methods: finding all paths up to certain
length;
extracting the giant component of a graph
- Caching
- Logging
- Configuration handling

[OncoSimulR](/packages/OncoSimulR)
----------

                         Changes in version 3.0                         

- MAJOR change: frequency dependent fitness available.

- Removed v.1 functionality.

- Multiple initMutants.

- Added MAGELLAN's sources and functionality from MAGELLAN.

                Changes in version 2.99.93 (2021-04-30)                 

- Fixed bugs and improved testing of rfitness with
  three-element scale vector.

                Changes in version 2.99.92 (2021-04-27)                 

- rfitness: scale can take a three-element vector.

- Vignette: examples (not run) for deviations from SSWM.

                 Changes in version 2.99.9 (2021-04-22)                 

- Fixed date typo in one citation.

- We were inconsistent, allowing some examples of one gene.

- Readme: nem.

                 Changes in version 2.99.8 (2021-01-01)                 

- Removing unused code.

- Long tests: no longer using v.1.

                 Changes in version 2.99.7 (2020-12-30)                 

- Vignette: fixed two missing refs and add seed in two examples.

                 Changes in version 2.99.6 (2020-12-30)                 

- No longer v.1 functionality.

- Slightly faster vignette.

                 Changes in version 2.99.5 (2020-12-18)                 

- Random timeouts when building in tokay2 (Windows, BioC); fixing
  seed in vignette.

                 Changes in version 2.99.4 (2020-12-17)                 

- Clean up of C++ code.

                 Changes in version 2.99.3 (2020-12-13)                 

- Remove unnecessary (and cluttering) output and irrelevant
  warnings when running tests.

- Decrease execution time of longer running examples in man (Rd) files.

- Decrease time of vignette.

                 Changes in version 2.99.2 (2020-12-11)                 

- Failed on test on Mac.

                 Changes in version 2.99.1 (2020-12-10)                 

- Bump version number for BioC, so it will become version 3.0.0 in
  next release.

- Latest version of exprtk.

                Changes in version 2.21.995 (2020-12-09)                

- Vignette: rewrote most FDF examples using names (not numbers)
  for fitness specification.

                Changes in version 2.21.994 (2020-12-08)                

- Can start simulation from arbitrary configuration: multiple init
  mutants (and multispecies functionality).

- Freq-dep-fitness does not need to have a WT in fitness tables.

- Bumped version (to 2.21.xyz) for new BioC devel.

[onlineFDR](/packages/onlineFDR)
---------

                        Changes in version 2.0.0                        

MODIFICATIONS

- added Rcpp code for online testing algorithms

- added online batch algorithms of Zrnic et al. [2020]

- added Storey-BH algorithm

- added setBound function

- updated vignette and pkgdown site

- updated unit tests

- updated references

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

                        Changes in version 2.9.1                        

- version bump to match Bioconductor version

[pathview](/packages/pathview)
--------

                       Changes in version 1.31.3                        

- fix bug in mol.sum with single gene data reported by
  easygsea. Similarly add drop=F to a few other lines in mol.sum and
  node.map.

                       Changes in version 1.31.2                        

- solve the check error due to the change related to KEGGEdgeSubtype
  in KEGGgraph package (version 1.51.1).

                       Changes in version 1.31.1                        

- korg now include 6833 KEGG species or 1588 new species beyond 2017.

[PCAtools](/packages/PCAtools)
--------

                        Changes in version 2.4.0                        

- added DESeq2 section to vignette intro

- permit that users can now specifiy just a single PC for
  plotloadings()

- improved ellipse functionality with addition of parameters
  ellipseType,
  ellipseLevel, and ellipseSegments

- over-rides ggrepel's new default value for max.overlaps via
  introduction of
  maxoverlapsConnectors = 15

- user can now specify a direction for connectors via
  directionConnectors

- removed labhjust and labvjust

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.5.4 (2021-04-15)                  

- mzML parameter parsing optimized

                 Changes in version 1.5.3 (2021-03-30)                  

- Additional calculation of raw peak area without smoothing (peakArea
  Raw)

- Alignment of peak integration algo across functions, integrateFIR()
  is now resilient to missing scans

                 Changes in version 1.5.2 (2021-01-31)                  

- Change package alias due to overwritten .rd file warning on Windows
  (non case-sensitive name and path)

- Move to Github Actions continuous integration

                 Changes in version 1.5.1 (2021-01-19)                  

- Unittests updated to comply with new r-devel all.equal() environment
  checks behaviour

[PepsNMR](/packages/PepsNMR)
-------

                 Changes in version 1.9.1 (2021-01-14)                  

CI TOOLS

- Switched from travis to GitHub actions

BUG CORRECTION 

- Debugged Normalization with type.norm = "firstquartile"

[PFP](/packages/PFP)
---

                       Changes in version 0.99.12                       

- Remove bugs in *.R, date: 2020.12.19.

                       Changes in version 0.99.11                       

- We update the package after the first review.

                       Changes in version 0.99.7                        

- 1.  The error reporting problem of THRESH_slot =NULL in rank_PFP is
- solved in PFP-class.R.
- 2.  Implement grouping of update path information from web pages.
- Add the get_pathway_info function, the get_PFPRefnet function
changes the ID
- column for the pathway_info (with the species name added, such as
Hsa04010),
- and changes the trans_graph2pfprefnet.r operation for the
pathway_info
- in the load_kegg_graphs.R.
- 3.  Add pathway_info_hsa.rdata to data.

                       Changes in version 0.99.6                        

- We update the package after the first review.

                       Changes in version 0.99.5                        

- We update R version dependency from 4.0.0 to 4.1.

                       Changes in version 0.99.4                        

- We subscribe to the bioc-devel mailing list

                       Changes in version 0.99.3                        

- Remove bugs in *.R, date: 2020.11.19.

                       Changes in version 0.99.2                        

- Remove bugs in *.R, date: 2020.11.8.

                       Changes in version 0.99.1                        

- Remove bugs in *.R

                       Changes in version 0.99.0                        

- We build a package that implements the pathway fingerprint
framework.

[PhIPData](/packages/PhIPData)
--------

                 Changes in version 0.99.0 (2021-02-17)                 

- Submitted to Bioconductor

[PhosR](/packages/PhosR)
-----

                        Changes in version 1.1.9                        

- Fixes in documentation, examples and styles to meet BioC
requirements

                        Changes in version 1.1.8                        

- Fixes in documentation and examples to meet BioC requirements

                        Changes in version 1.1.7                        

- Fixed getSPS to expect residue information.

                        Changes in version 1.1.6                        

- Included warning signs to creating PhosphoExperiment object without
attributes.

                        Changes in version 1.1.5                        

- Fixed the warning message generated in plotQC for dendrograms.
(Warning message: Vectorized input to element_text() is not
officially supported.).

                        Changes in version 1.1.4                        

- Fixed the parameter of plotQC from cols to grps

                        Changes in version 1.1.3                        

- Fixed a bug in ptImpute function
- Edited medianScaling documentation to be more specific

[PhyloProfile](/packages/PhyloProfile)
------------

                       Changes in version 1.4.12                        

- domain plot works with group ID containing pipe

[Pigengene](/packages/Pigengene)
---------

                Changes in version 1.17.10 (2021-02-02)                 

Bug Fixes

- When building the Ubuntu, Neda got an error message on building
  the vignette "argument is of length zero". It was due to her
  old version of BiocStyle, 2.14.4. To solve this issue,
  Pigengene now depends BiocStyle Version >= 2.19.1. Even 2.18.1
  might be enough with R 4.0.3.

[planet](/packages/planet)
------

                       Changes in version 0.99.4                        

- Combined vignettes into a single one
- Updated references

                       Changes in version 0.99.3                        

- Update reference for cell composition
- Accepted into bioconductor, will be released in next cycle

                       Changes in version 0.99.2                        

- Vignettes updated with new function names
- Removed minfi from examples since installation fails on some R
versions
- Moved .orig files to inst/script
- Moved data-raw/ contents to inst/script
- Moved man/figures/logo.png to inst/figures/logo.png, updated README
- Lazydata set to false

                       Changes in version 0.99.0                        

- Improved vignette figure quality
- Using github-actions check from biocthis
- Removed viridis from suggests
- Rewrote many data man pages

                        Changes in version 0.3.0                        

- Final preparations for Bioconductor submission
- Changed license to => GPL-2
- More detailed attribution to glmnet
- Renamed all functions to camel case convention, including
functions,
data, tests, documentation, vignettes
- Added planet-deprecated.R
- pl_infer_ethnicity() and pl_infer_age() now deprecated, replaced
with predictEthnicity and predictAge, respectively.


[pmp](/packages/pmp)
---

                        Changes in version 1.3.4                        

- Added option to constrain spar values for QCRSC.

                        Changes in version 1.3.2                        

- Added inputs/outputs for pqn.

- Fixed glog plot bug.

                        Changes in version 1.2.1                        

- Account for missing injections in QCRSC.

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 0.99.6 (2021-04-06)                 

- Updated NEWS.Rd)

                 Changes in version 0.99.5 (2020-12-09)                 

- Corrected name of an author)

                 Changes in version 0.99.4 (2020-12-09)                 

- Changed R version dependency to >= 4.1 (R-devel)

                       Changes in version 0.99.3                        

- Minor code changes as suggested by Bioconductor reviewer

                 Changes in version 0.99.2 (2020-11-11)                 

- Changed arguments in examples to run faster

                 Changes in version 0.99.1 (2020-11-11)                 

- Removed .Rproj file from repository

                 Changes in version 0.99.0 (2020-11-10)                 

- Submitted to Bioconductor

[podkat](/packages/podkat)
------

                       Changes in version 1.24.0                        

- release as part of Bioconductor 3.13

                       Changes in version 1.23.2                        

- updates and documentation to vignette in order to adapt to newer
  version of
  Illumina's TruSeq DNA Exome library prep kit

- changed default col.names in readRegionsFromBedFile(); corresponding
  update of
  help page

                       Changes in version 1.23.1                        

- re-created genome data objects (old data objects had become
  incompatible
  with newer 'BSgenome' version)

                       Changes in version 1.23.0                        

- new branch for Bioconductor 3.13 devel

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

                       Changes in version 1.1.15                        

- Remove reshape2 and Biobase packages from Imports
- Implement viridis palette for PomaBoxplots, PomaDensity, and
PomaMultivariate
- Update mixOmics output names in PomaMultivariate
- New package description and biocViews
- Bug fixed in PomaMultivariate

                        Changes in version 1.1.8                        

- PomaNorm and PomaImpute warnings when methos parameter is missing
passed to messages
- Minor bugs fixed
- Minor changes in plots style
- plotly package used in PomaVolcano switched to Suggests
- Bioconductor badge table added to README

POMA 1.0.0

- Released to Bioconductor 3.12

[POWSC](/packages/POWSC)
-----

                        Changes in version 0.99                         

NEW FEATURES

- Initial review.

                       Changes in version 0.99.0                        

- Revise required files and format the code style.

[proActiv](/packages/proActiv)
--------

                       Changes in version 1.1.22                        

- Update package-down site

- Update NEWS for Bioconductor 3.13

                       Changes in version 1.1.18                        

- Gene expression data is now stored in the assays of the
  summarizedExperiment
  object returned by proActiv to facilitate easier filtering of the
  summarizedExperiment object. The metadata slot is now empty.

- Plotting promoter activity: Implementation of boxplotPromoters
  function to
  plot boxplots of absolute promoter activity, relative promoter
  activity, and
  gene expression.

- Identification of alternative promoters: Implementation of
  getAlternativePromoters, used to identify promoters that may exhibit
  alternative
  usage.

                       Changes in version 1.1.15                        

- Implement getAlternativePromoters for identifying alternative
  promoters

- Implement boxplotPromoters for visualizing promoter usage

                        Changes in version 1.1.6                        

- Enforce condition vector to following naming conventions


[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.31                         

Changes in version 1.31.3

- lopims() function moved to lgatto/lopims package on GitHub

Changes in version 1.31.2

- Update dunkley2006params, andy2011params and MartInterfaces objects

Changes in version 1.31.1

- Suggest magick (needed to build vignette)

Changes in version 1.31.0

- New devel version (Bioc 3.13)

[ProtGenerics](/packages/ProtGenerics)
------------

                       Changes in version 1.23.9                        

- Added new uniqueMsLevel generic <2021-04-08 Thu>

                       Changes in version 1.23.8                        

- Added new filterPrecursorCharge generic <2021-04-06 Tue>

                       Changes in version 1.23.7                        

- Add ProcessingStep object and related methods (moved from Spectra)

                       Changes in version 1.23.6                        

- Added quantify generics (moved from MSnbase) <2021-01-02 Sat>

                       Changes in version 1.23.5                        

- new alignRt generic.

                       Changes in version 1.23.4                        

- new virtual Param class.

                       Changes in version 1.23.3                        

- new filterIntensity generic.

                       Changes in version 1.23.2                        

- new compounds generic.

                       Changes in version 1.23.1                        

- new calculateFragments generic (moved from MSnbase)

[psichomics](/packages/psichomics)
----------

                       Changes in version 1.18.0                        

- discardLowCoveragePSIvalues(): improve performance (2x faster)
- When quantifying or loading PSI values, psichomics discards
splicing
events whose junctions:
- (1) are not present in junction quantification data or
- (2) have low numbers of reads across all samples. The number of
discarded events is now displayed during PSI quantification.
- When normalizing gene expression, support converting gene
identifiers to gene symbol names for any species with OrgDb data
- Allow to use plotSplicingEvent() directly with a PSI table to plot
diagrams of the alternative splicing events within
- Distribution plots (plotDistribution()):
- Plot sample distributions in density, violin or box plots
(argument type in command-line interface)
- Hide rug plot when showing 500 or more values (by default) to
avoid performance issues
- Add jitter to rug plot (helps to discern numerous points)
- Display interquantile range (IQR) per group in the tooltip
- Add subtitles (argument subtitle in command-line interface)
- Allow to invert axes (argument invertAxes in command-line
interface)

Bug fixes and minor changes

- loadLocalFiles(): print elapsed time after loading local files
- filterGroups() was modified to return a character vector whose
names
are original names (instead of groups) and include an attribute
Groups with the respective group of each value (together with their
colour, if available):
- Distribution plots in the graphical interface now show sample
name in the tooltip for each sample
- Importing VAST-TOOLS annotation:
- Fix parsing of VAST-TOOLS intron retention events as skipped
exon for certain annotations (including Hs2)
- Support diagrams for intron retention events with full exon
coordinates
- Graphical user interface improvements:
- Avoid "Matching subjects to their samples/Matching process
concluded" loop
- Timeout GTEx data retrieval after 3 seconds without server
response
- Avoid unnecessary messages when loading Firebrowse interface
- Fix intron retention diagrams in distribution plots not
displaying introns
- When searching for specific splicing events, fix results based
on the wrong genomic coordinates
- Improve tutorials

[ptairMS](/packages/ptairMS)
-------

                 Changes in version 0.99.5 (2021-04-08)                 

- folder src-x64 and src-i386 in gitignore

                 Changes in version 0.99.4 (2021-04-08)                 

- vignette: eval=FALSE in installation chuck

                 Changes in version 0.99.3 (2021-04-01)                 

- added BugReports to description files

- changed in vigette: remove echo = FASLSE and added installation
  section

- remove dontrun from examples

- avoided slot() and created generic method for slot access and setters

- modyfied the two datSet documentaion

- added documentation in inst/script

                 Changes in version 0.99.0 (2021-03-12)                 

- Submitted to Bioconductor

[PureCN](/packages/PureCN)
------

                       Changes in version 1.22.0                        

NEW FEATURES

- calculateNormalDatabase now suggests an off-target interval width
  that
  minimizes noise while keeping the resolution as high as possible

- Added support for GATK4 CollectAllelicCounts output as alternative
  to Mutect

- Added segmentationGATK4 to use GATK4's segmentation function
  ModelSegments

SIGNIFICANT USER-VISIBLE CHANGES

- Added min.total.counts filter to filterIntervals to remove
  intervals with low number of read counts in combined tumor and
  normal.
  Useful especially for off-target filtering in highly efficient assays
  where standard filters keep too many high variance regions.

- Changed default of min.mappability in preprocessIntervals for
  on-target
  intervals to 0.6 (from 0.5)

- Added min.mappability also to filterIntervals so that more
  conservative
  cutoffs can be tested after normalDB generation

- PSCBS: 1.20.0 two-step segmentation slightly tweaked in that only
  high quality on-target intervals (high mappability and low PoN noise)
  are used in the first segmentation

- Added --skipgcnorm flag to Coverage.R to skip GC-normalization

- Added AF.info.field option to calculateMappingBiasGatk4 for
  non-standard
  GenomicsDB imports

- If segmentation functions add breakpoints within baits, these
  breakpoints are now moved to the beginning or end of that bait to
  avoid
  that a single bait is assigned to two segments

- Dx.R now always generates a _signatures.csv file with --signatures,
  even
  if insufficient number of mutations

- Removed defunct calculateIntervalWeights function

BUGFIXES

- Fix for nonsensical error message when VCF does not contain germline
  variants (#166).

- Fix for various issues related to the seqlevelsStyle function (e.g.
  #171)

- Fix for crash in calculateMappingBiasGatk4 when not all samples had
  a single variant call on a particular chromosome (chrY)

- Fix related to annotating mapping bias with triallelic sites and
  GenomicsDB

- Fixed an issue in Mutect 1.1.7 data in which good SNPs were ignored
  (#174)

[qcmetrics](/packages/qcmetrics)
---------

                        Changes in version 1.29                         

qcmetrics 1.29.1

- Clean up (remove .gitignores for branch-specific .gitingore,
README.Rmd).
- Remove dependency on yaqcaffy (that is being deprecated).

qcmetrics 1.29.0

- The Bioc devel version

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.1.0                        

QFeatures 1.1.4

- Added replaceRowDataCols and removeRowDataCols, two functions to
streamline manipulation of rowData within a QFeature object.

QFeatures 1.1.3

- Added countUniqueFeatures, a function to count the number of unique
features per sample.

QFeatures 1.1.2

- Manually install preprocessCore (see
https://github.com/Bioconductor/bioconductor_docker/issues/22 for
details) to use quantile normalisation in vignette and tests.

- Update vignette to show normalize() and logTransform() directly on
a
QFeatures object and reference the QFeaturesWorkshop2020 workshop
and WSBIM2122 chap 8.

QFeatures 1.1.1

- Fix typo in vignette
- Improve first paragraph in intro vignette.

QFeatures 1.1.0

- New Bioconductor devel version

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

                        Changes in version 1.9.4                        

- Added documentation for new data set.

                        Changes in version 1.9.3                        

- New Feature: Added the function IRSnorm, which performs normalisation
  of
  intensities across multiple TMT runs using a common reference sample
  (Internal
  Reference Scale).

                        Changes in version 1.9.2                        

- Bug Fix: Update getContrastResults test reference

- Fix bug in is_validScalingFunction - switch from `are_equal` to
  `identical`

                        Changes in version 1.9.1                        

- Fix bug in getContrastRestults where file was written to `.txt`, file
  name
  should now be correctly formed from the contrast name

[quantiseqr](/packages/quantiseqr)
----------

                       Changes in version 0.99.0                        

- All set for Bioconductor submission!

                        Changes in version 0.0.1                        

- Added a NEWS.md file to track changes to the package.

[RadioGx](/packages/RadioGx)
-------

                        Changes in version 1.0.1                        

- Abstract further functionality from RadioGx into CoreGx dependency
package
- Harmonize function and API design to be in line with PharmacoGx R
package

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.16.0                        

New features

- sparseAssay and compactAssay now support sparseMatrix outputs from
the Matrix package

[ramr](/packages/ramr)
----

                 Changes in version 0.99.2 (2020-12-21)                 

- code cleanup for bioconductor

- initial submission to bioconductor


[randRotation](/packages/randRotation)
------------

                 Changes in version 1.3.3 (2020-11-18)                  

- Defunct function df_estimate()

                 Changes in version 1.2.1 (2020-11-17)                  

- pFdr adapted according to (Phipson and Smyth 2010)

[Rbowtie](/packages/Rbowtie)
-------

                       Changes in version 1.31.1                        

NEW FEATURES

- updated bowtie to version 1.3.0 (bugfixes, improved compiler support,
  removed support for colorspace alignments)

[RcisTarget](/packages/RcisTarget)
----------

                        Changes in version 1.11                         

- importRankings() updates for extended compatibility:
  * indexCol by default is now 'NULL', which will take the first
  column.
  If it is set to any other column, that one will be used instead.

[Rcwl](/packages/Rcwl)
----

                 Changes in version 1.7.12 (2021-02-24)                 

- Use Basilisk to manage python dependencies (cwltool).

- Rename of 'cwlParam' into 'cwlProcess', 'cwlStepParam' into
  'cwlWorkflow'.

- Add the support of wrapping R functions into Rcwl tools.

[RcwlPipelines](/packages/RcwlPipelines)
-------------

                 Changes in version 1.7.7 (2021-02-24)                  

- Moved all recipes to https://github.com/rworkflow/RcwlRecipes.

- Added new core functions: cwlUpdate, cwlSearch and cwlLoad.

- Core functions return a 'cwlHub' object.

[RCy3](/packages/RCy3)
----

                       Changes in version 2.12.0                        

- New support for cloud-hosted Jupyter notebooks!
  - Jyputer notebooks running RCy3 in the cloud can communicate with
  local Cytoscape instances
  - Includes a Sandbox mechanism to manage file transfers between cloud
  and local dir

- New functions:
  - addAnnotationText
  - getAnnotationList
  - deleteAnnotation
  - setters for Filter, Model Propagation and Catchup delays
  - many sandbox-related functions
  - getStyleMapping
  - getAllStyleMappings
  - paletteColorRandom
  - paletteColorBrewer* (33 in total!)
  - internal mapping value generators for color, opacity, dimension,
  line styles, arrows and shapes

- New parameters:
  - overwrite_file added to export functions
  - apply added to create_*_filter functions
  - ndex.url and ndex.version added to CyNDEx functions

- Consistent interchangeable handling of node|edge|network names and
  SUIDs

- createNetworkFromDataframes plays nice with tibbles

- Bug fixes:
  - .edgeNameToEdgeSUID revamped to better handle duplicate edge names

[ReactomeContentService4R](/packages/ReactomeContentService4R)
------------------------

                         Changes in version 1.0                         

- First release

[ReactomeGraph4R](/packages/ReactomeGraph4R)
---------------

                     Changes in version 0.0.0.9000                      

- Added a NEWS.md file to track changes to the package.

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.5.2 (2021-04-16)                  

- Added check to ensure that clustering was performed in Seurat objects
  prior to the pathway analysis.

                 Changes in version 1.5.1 (2021-04-01)                  

- Fixed bug in perform_reactome_analysis: Error messages were not
  displayed correctly.

[rebook](/packages/rebook)
------

                        Changes in version 1.2.0                        

- 
  Added the rmd2id() function to easily determine the ID for each
  chapter.

- 
  Link to the originating chapter for the set-up code in
  extractCached().

- 
  Expose collapseStart() and collapseEnd() for manual creation of
  collapsible chunks.

- 
  Added scrapeReferences() to scrape a bookdown book for
  references for external use.

- 
  Added configureBook() to configure a Bioconductor package as a
  book deployment.

- 
  Added link() to rapidly link to references in a configured
  Bioconductor book package.

- 
  Added extractFromPackage() to extract objects from Rmarkdown
  files in installed packages.

- 
  Added createRedirects() to redirect from old, deprecated pages
  to their new locations.

[recount3](/packages/recount3)
--------

                        Changes in version 1.1.7                        

BUG FIXES

- read_counts() now reads the gene/exon counts for every sample as
numeric instead of integer in order to support count values that
exceed the 32bit integer threshold (such as 2447935369). Previously,
read_counts() would report tiny fractions for such large numbers.
This bug was reported by Christopher Wilks.

                        Changes in version 1.1.6                        

BUG FIXES

- Now project_homes() reads in a text file from
recount3_url/organism/homes_index which enables support for custom
URLs such as http://snaptron.cs.jhu.edu/data/temp/recount3test.

                        Changes in version 1.1.4                        

BUG FIXES

- Fixed project_homes(), available_projects() and available_samples()
to support using non-standard recount3_urls where the user knows
that are the project_homes() for their organism of choice. This fix
enables users to create their own custom recount3-like webservers
and access their data using the functions in this package. This fix
introduces the argument available_homes to both available_projects()
and available_samples(). This bug was reported by Christopher Wilks.

                        Changes in version 1.1.3                        

NEW FEATURES

- Added expand_sra_attributes() that was contributed by Andrew E
Jaffe. This function expands the SRA attributes stored in a given
SRA study, which makes it easier to use that data. However, it makes
it harder to merge studies and thus should be used with caution.

                        Changes in version 1.1.1                        

BUG FIXES

- Fixed locate_url() for GTEx & TCGA BigWig files.

[recountmethylation](/packages/recountmethylation)
------------------

                        Changes in version 1.1.4                        

- Improves User's Guide, fixes typos, new citation, adds
  disclaimer text, updates servermatrix() chunk.

- Updates examples to further limit downloads. Function
  get_rmdl() uses download = FALSE,

- Updates Data Analyses vignette. Uses reduced dpi for images in
  Data Analyses vignette to limit package size. Uses updated
  metadata file name.

- Compresses new metadata v.0.0.2 files to limit package size.

- Uses uniform metadata file label for v.0.0.1 file.

                        Changes in version 1.1.3                        

- Added v.0.0.2 database compilation files to server
  (recount.bio/data) and revised recountmethylation functions for
  cross-platform support.  The new files reflect IDAT downloads
  completed in Nov 2020 from GEO/GDS, including the first
  compilations of EPIC/HM850K arrays.

- Added `platform` argument in relevant `getdb` functions.

- Added `which.platform` argument to `get_rmdl`

- Added new function `smfilt` to filter server data table for
  newest compilation files, accounting for platform in the
  filename.

- Cleaned up and shoretened the `servermatrix` function. This now
  handles RCurl call for "dn" (originally from `get_rmdl`) when
  handling condition `dn = NULL`.

- Updated the User's Guide to fix typos, reflect v.0.0.2 samples,
  and add a download troubleshoot section. Added numeric
  citations format, removed evaluation of validation section due
  to possible package build failure from bad internet connection.

- Updated ExperimentHub file metadata script and table to add new
  v.0.0.2 compilation files.

- Renamed sample metadata directory to "gsm_metadata" to avoid
  confusion with "metadata.csv" file table for hubs.


[regutools](/packages/regutools)
---------

                        Changes in version 1.3.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- The link provided inside connect_database() has been updated in
order to connect with the latest version of regulonDB v10.8.

                        Changes in version 1.3.1                        

BUG FIXES

- Fixed URL for the database. Resolves
https://github.com/ComunidadBioInfo/regutools/issues/33.

[RGMQL](/packages/RGMQL)
-----

                       Changes in version 1.11.1                        

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

                       Changes in version 2.36.0                        

NEW FEATURES

- Added additional hyberslab selection functions introduced in
  HDF5 1.10.7 (H5Scombine_hyperslab, H5Scombine_select,
  H5Sget_select_npoints).

- Support for read access to files in S3 buckets now includes
  Windows.

- Added function h5deleteAttribute().

BUG FIXES

- Addressed issue where messages printed when loading .Rprofile
  were breaking detection of the rhdf5filters package.
  (https://github.com/grimbough/rhdf5/issues/81)

[rhdf5filters](/packages/rhdf5filters)
------------

                        Changes in version 1.4.0                        

USER VISIBLE CHANGES

- Compression libraries updated:

- blosc: 1.16.3 🠪 1.20.1
- lz4: 1.8.3 🠪 1.9.2
- zstd: 1.3.8 🠪 1.4.5

- Added LZF filter

BUG FIXES

- Fixed some missing references to CFLAGS, CPPFLAGS & LDFLAGS in
package compilation
(https://github.com/grimbough/rhdf5filters/pull/4)
- Improved CPU detection on non-x86 architecture to allow compilation
to proceed with default instruction sets, rather than failing to
install. (https://github.com/grimbough/rhdf5filters/issues/3)
- Address issue in compilation where message printed while processing
.Rprofile where added to LDFLAGS
(https://github.com/grimbough/rhdf5filters/issues/11)(#11)

[Rhdf5lib](/packages/Rhdf5lib)
--------

                        Changes in version 1.14                         

New features

- Updated internal version of HDF5 to 1.10.7

Bug fixes

- CPPFLAGS used to build R are now used during HDF compilation.

- Package configure script will now stop if it encounters errors when
  compiling HDF5.  This should make diagnosing issues easier.

[RiboDiPA](/packages/RiboDiPA)
--------

                       Changes in version 0.99.3                        

- Renamed R functions by following bioconductor convention

                       Changes in version 0.99.2                        

- Initial version for review

[ribosomeProfilingQC](/packages/ribosomeProfilingQC)
-------------------

                        Changes in version 1.3.4                        

- add rmarkdown as suggest package.

                        Changes in version 1.3.3                        

- add github action.

                        Changes in version 1.3.2                        

- fix the issue if there is softclip in the mapping reads and the
reads length is smaller than shfit range.

                        Changes in version 1.3.1                        

- keep the raw counts for countReads.

[RLassoCox](/packages/RLassoCox)
---------

                 Changes in version 0.99.0 (2020-10-21)                 

- Submitted to Bioconductor

[rmelting](/packages/rmelting)
--------

                        Changes in version 1.7.1                        

- Fix minor error in links.

[RNAmodR](/packages/RNAmodR)
-------

                 Changes in version 1.5.4 (2021-01-23)                  

- add plot type "points" to plotCompare functions

                 Changes in version 1.5.3 (2021-01-12)                  

- bugfix for zero-length annotation

- bugfix for settings function

- bugfix for plotting functions

                 Changes in version 1.5.2 (2020-12-12)                  

- bugfix for names plotted by plotCompareByCoord

[RNAmodR.AlkAnilineSeq](/packages/RNAmodR.AlkAnilineSeq)
---------------------

                 Changes in version 1.5.1 (2021-01-12)                  

- bugfix for settings function

[RNAmodR.RiboMethSeq](/packages/RNAmodR.RiboMethSeq)
-------------------

                 Changes in version 1.5.1 (2021-01-12)                  

- bugfix for settings function

[RnBeads](/packages/RnBeads)
-------

                        Changes in version 2.9.3                        

- Changed some of the default options values

- Some bugfixes

                        Changes in version 2.9.2                        

- Implemented function rnb.execute.pOOBAH. Thanks to Nathan Steenbuck.

- Update of contact information in DESCRIPTION

                        Changes in version 2.9.1                        

- Improved function intensities.by.color. Thanks to Nathan Steenbuck.

[rols](/packages/rols)
----

                        Changes in version 2.19                         

CHANGES IN VERSION 2.19.4

- Replace failing unit tests after GO:0032801 changed.

CHANGES IN VERSION 2.19.3

- Fix error for empty query results
- New as.data.frame OlsSearch S3 method

CHANGES IN VERSION 2.19.2

- Fix failing unit test (number of obsolete terms in test_Terms.R)

CHANGES IN VERSION 2.19.1

- Fix pagination of ebi.ac.uk results contributed by Andrew Clugston
(AndrewC160) - see issue #27.
- For failing unit test (number of obsolete terms in test_Terms.R)

[ropls](/packages/ropls)
-----

                       Changes in version 1.23.12                       

MINOR MODIFICATION

- minor documentation update

                       Changes in version 1.23.8                        

MINOR MODIFICATION

- vignette: minor update

                       Changes in version 1.23.6                        

MINOR MODIFICATION

- 'view' method: minor update of documentation

                       Changes in version 1.23.4                        

BUG FIXED

- 'view' method: display of row names

                       Changes in version 1.23.2                        

BUG FIXED

- 'view' method: display of row names (reversed in the previous
  versions)

[ROSeq](/packages/ROSeq)
-----

                       Changes in version 1.99.01                       

- no bug, but removed extra files from man folder

[rpx](/packages/rpx)
---

                       Changes in version 1.99.0                        

rpx 1.99.8

- Remove rappdirs from Suggests.

rpx 1.99.7

- Add markdown to Suggests.

rpx 1.99.6

- Switching to tools::R_user_dir() from rappdirs::user_cache_dir() to
set package cache directory, following changes in BiocFileCache (see
for details).

rpx 1.99.5

- User can pass an own cache when downloading files <2021-04-03 Sat>

rpx 1.99.4

- Automatically detect/fix correct url <2021-03-22 Mon>

rpx 1.99.3

- Use MS:1002852 to get URL when PRIDE:0000411 returns empty string.

rpx 1.99.2

- Don't ask for cache location in non-interactive mode. This forces
the usage of the default cache location (rather than usage of a Rtmp
directory) and enables reusage of a persistent cache over different
build/check cycles.

rpx 1.99.1

- Using BiocFileCache to cache PXD files that are downloaded.

Changes in version 1.27.0

- New devel version

[rrvgo](/packages/rrvgo)
-----

                        Changes in version 1.3.1                        

- Bug fix for obsolete GO terms


[Rsamtools](/packages/Rsamtools)
---------

                         Changes in version 2.8                         

NEW FEATURES

- (v 2.7.2) idxstatsBam works on remote (e.g., http://) files and
  reports
  unmapped ('seqnames' equal to *) reads. See
  https://support.bioconductor.org/p/9136222.

[Rsubread](/packages/Rsubread)
--------

                        Changes in version 2.6.0                        

- 
  Improved the speed of cellCounts and also reduced its memory
  use.

- 
  Added a parameter 'umi.cutoff' to cellCounts to call all the
  cells that had a total UMI count greater than the specified
  threshold.

- 
  Added support for FASTQ-format read input in CellCounts.

[rSWeeP](/packages/rSWeeP)
------

                 Changes in version 1.3.2 (2020-11-24)                  

- baseMatrix error fix


[RTCGAToolbox](/packages/RTCGAToolbox)
------------

                       Changes in version 2.22.0                        

New features

- makeSummarizedExperimentFromGISTIC transferred from TCGAutils

Bug fixes and minor improvements

- biocExtract now merges datasets from the same platforms
- Re-worked and simplified getData method for FirehoseData and
FirehoseGISTIC
- Remove missing ranges when creating GRanges and GRangesList from
DataFrame

[RTN](/packages/RTN)
---

                       Changes in version 2.16.0                        

- Regular maintenance, stable release.

[RTNduals](/packages/RTNduals)
--------

                       Changes in version 1.16.0                        

- Regular maintenance, stable release.

[RTNsurvival](/packages/RTNsurvival)
-----------

                       Changes in version 1.16.0                        

- Regular maintenance, stable release.

[Rtpca](/packages/Rtpca)
-----

                 Changes in version 1.1.2 (2020-11-27)                  

- Significant improvement in performance of all major functions

[rWikiPathways](/packages/rWikiPathways)
-------------

                       Changes in version 1.12.0                        

- Removed getColoredPathway

- New function: writeGMT

- Updated URLs to BridgeDb datasources.tsv (again)

- Bug fix: downloadPathwaysArchive works with redirected urls

                       Changes in version 1.11.4                        

- Updated URLs to BridgeDb datasources.tsv

                       Changes in version 1.11.3                        

- New features
  - new readGMT and readGMTnames methods

[SAIGEgds](/packages/SAIGEgds)
--------

                        Changes in version 1.6.0                        

- `seqFitNullGLMM_SPA()` can use imputed dosages directly without
  converting the dosages to the best-guess genotypes

- new function `glmmHeritability()` for approximate heritability
  estimates

[sampleClassifier](/packages/sampleClassifier)
----------------

                       Changes in version 1.16.0                        

BUG FIXES

- fixed the error in classifyProfile.svm() and
  classifyProfile.rnaseq.svm() by
  changing the vector of Profile names in the call of svm() function
  from character to
  factor.

[SANTA](/packages/SANTA)
-----

                 Changes in version 0.99.0 (2021-01-06)                 

- Re-submitted to Bioconductor

[satuRn](/packages/satuRn)
------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- Your main changes to a function foo() or parameter param.

BUG FIXES

- Your bug fixes. See more details at
http://bioconductor.org/developers/package-guidelines/#news.

[SBGNview](/packages/SBGNview)
--------

                        Changes in version 1.5.1                        

- Massive updates and rewriting of SBGNview, by Kovidh Vegesna and
  Weijun Luo.

[scAlign](/packages/scAlign)
-------

                 Changes in version 1.5.60 (2021-05-18)                 

- Saving to tempdir instead of user home

[SCArray](/packages/SCArray)
-------

                        Changes in version 1.0.0                        

- initial version of SCArray

[scater](/packages/scater)
------

                       Changes in version 1.20.0                        

- runMDS can use user-supplied function for calculating the
  distance matrix. runMDS can optionally store the distance
  matrix computed. runMDS also stores the eig and GOF fields of
  the object returned.

- Made the handling of center, scale, color and limits similar in
  plotDots, plotHeatmap, and plotGroupedHeatmap

- Add use_fitsne argument to runTSNE allowing the use of fast
  interpolated t-SNE in place of Barnes-Hut t-SNE.

[scClassifR](/packages/scClassifR)
----------

                 Changes in version 0.99.5 (2020-12-30)                 

- Adapted code based on Bioconductor review

[scDblFinder](/packages/scDblFinder)
-----------

                 Changes in version 1.5.11 (2021-01-19)                 

- scDblFinder now provides doublet enrichment tests

- doublet generation and default parameters have been further optimized

[scone](/packages/scone)
-----

                 Changes in version 1.15.1 (2021-05-06)                 

- Add PsiNorm normalization method and its wrapper.

- Add vignette that describes how to use PsiNorm.

[scp](/packages/scp)
---

                         Changes in version 1.1                         

scp 1.1.6

- feature: readSCP now allows for a suffix argument to better
customize the sample names. <2021-03-17>

scp 1.1.5

- deprecation: thanks to the new normalization method in
medianCVperCell, 'computeMedianCV_SCoPE2' is now deprecated and
should no longer be used. <2021-02-19>
- feat: added a new normalization method to medianCVperCell. The
SCoPE2 normalization method can now reproduce the results from
SCoPE2. <2021-02-19>
- docs: improved vignette <2021-02-16>
- feat: added a rowDataName argument to computeSCR <2021-02-08>

scp 1.1.4

- fix: removed bug in vignette header <2021-02-06>
- data: update the example data with the latest release of SCoPE2
<2021-02-06>
- feat: added removeEmptyCol argument in readSCP to automatically
remove columns that contain only NA's <2021-02-06>

scp 1.1.3

- docs: improved the manual page for pep2qvalue and the corresponding
section in the vignette. <2021-01-23>
- refactor: reimplemented the computeFDR to catch up with the new
release of SCoPE2. computeFDR was renamed to pep2qvalue. This is
more in line with the theory. Also adapted the unit tests.
<2021-01-23>

scp 1.1.2

- docs: improved the description of the scp data structure in the
vignette <2021-01-05>
- refactor: renamed the groupCol to groupBy and pepCol to PEP in
computeFDR. <2020-12-08>
- refactor: renamed computeMedianCV to computeMedianCV_SCoPE2 and
deprecated the function. The function will be preserved for backward
compatibility with the replication of the SCoPE2 analysis (Specht et
al. 2020). Instead, a new function is implemented and called
medianCVperCell. See issue#7 for more information <2020-12-07>

scp 1.1.1

- Fix news file

scp 1.1.0

- New devel (Bioc 3.13)

[scPCA](/packages/scPCA)
-----

                 Changes in version 1.5.3 (2021-03-14)                  

- Updating plotting issue in vignette: comparison of cPCA and scPCA
loadings.
- Adding pkgdown site.
- Moving ScaledMatrix to "imports" section of DESCRIPTION.

                 Changes in version 1.5.2 (2020-12-21)                  

- Adding LTLA/ScaledMatrix to "Remotes" section of DESCRIPTION.

                 Changes in version 1.5.1 (2020-12-17)                  

- scPCA() and other internal functions may now take advantage of the
ScaledMatrix object class. This allows more computationally
efficient contrastive covariance matrix estimation when analyzing
large datasets.
- safeColScale() now used MatrixGenerics to handle feature
standardization.

[scran](/packages/scran)
-----

                       Changes in version 1.20.0                        

- All deprecated functions from the previous release are now
  defunct.

- Added a simplify= option to quickSubCluster() to get the
  cluster assignments directly.

- Deprecated combinePValues() as this is replaced by
  metapod::combineParallelPValues().

- getClusteredPCs() now uses bluster::clusterRows() by default.

- decideTestsPerLabel() now automatically detects pval.field= if
  not supplied.

- Added the clusterCells() wrapper around bluster functionality.

- Removed the option to pass a matrix in design= from
  pseudoBulkDGE().

- Migrated all normalization-related functions
  (computeSumFactors(), calculateSumFactors(), cleanSizeFactors()
  and computeSpikeFactors()) to a better home in scuttle.
  Soft-deprecated existing functions.

- Modified getTopHVGs() to accept a SingleCellExperiment and
  compute the DataFrame with modelGeneVar().

- Added fixedPCA() to compute a PCA with a fixed number of
  components, a la scater::runPCA() (but without requiring
  scater).

- Modified denoisePCA() so that it now complains if subset.row=
  is not provided.

- Modified all pairwise* functions so that the p-value from
  direction="any" is derived from the two p-values from the
  one-sided tests. This is necessary for correctness with all
  choices of lfc= and block=, at the cost of conservativeness
  when block=NULL and lfc is large.

[scRepertoire](/packages/scRepertoire)
------------

                        Changes in version 1.2.3                        

- Changed the access of the sample data to github.io repo:
  readRDS(url("https://ncborcherding.github.io/vignettes/scRepertoire_example.rds"))

                        Changes in version 1.2.2                        

- Removed Startrac-based functions in order to pass build on
  Bioconductor.
  DEPRECATED AND DEFUNCT

- Deprecate StartracDiversity()

                        Changes in version 1.2.0                        

SUBMITTED

- The first version of scRepertoire submitted to Bioconductor.

SIGNIFICANT USER-VISIBLE CHANGES

- Added support for SingleCellExperiement format.

DEPRECATED AND DEFUNCT

- Deprecate combineSeurat in favor or combineExpression().

- Deprecate seurat2List in favor of expression2List().

                        Changes in version 1.1.2                        

- Clonal Overlap Coefficient issue fixed, was comparing unique barcodes
  and not clonotypes

- Added function checkBlanks to remove list elements without
  clonotypes, this prevents errors for visualizations

- Re-added Startrac metrics by stripping down the package and adding it
  piecemeal

- Heavily modified dependencies to reduce total number

[scuttle](/packages/scuttle)
-------

                        Changes in version 1.2.0                        

- Migrated whichNonZero() to beachmat.

- Bugfixes for factor-based colData aggregation in
  aggregateAcrossCells(). Added proper support for Vectors.

- Bugfix for correct response to use.altexps= in
  perCellQCMetrics(), perFeatureQCMetrics().

- Added a normalize.all= option to normalizeCounts(). Removed
  unnecessary warning when down.target= is not specified. Exposed
  the default size.factors= in the SingleCellExperiment method.

- Modified the SingleCellExperiment method of logNormCounts() so
  that manually specified size factors do not apply to
  alternative Experiments. Only relevant if size.factors= and
  use.altexps= are specified.

- Deprecated use.altexps= in favor of applySCE() in
  logNormCounts() and aggregateAcrossCells().

- Renamed addPerCellQC() and addPerFeatureQC() to
  addPerCellQCMetrics() and addPerCellFeatureMetrics(), for
  consistency. Soft-deprecated the old functions.

- Moved most of quickPerCellQC() functionality into the new
  perCellQCFilters() function. Repurposed the former to directly
  return a filtered SummarizedExperiment object.

- Migrated scran's normalization-related functions into this
  package. Added pooledSizeFactors(), computePooledFactors(),
  cleanSizeFactors() and computeSpikeFactors().

- Added transform="asinh" to normalizeCounts() and
  logNormCounts() for inverse hyperbolic transformations of
  CITE-seq data.

- Modified isOutlier() to now return outlier.filter objects.
  These are simply logical vectors that preseve the "thresholds"
  attribute upon subsetting.

- Migrated correctGroupSummary() from scater, to compute
  corrected versions of group-level summary statistics.

[SDAMS](/packages/SDAMS)
-----

                       Changes in version 1.11.5                        

- return beta and gamma as matrix

                       Changes in version 1.11.4                        

- update vignette with SC RNA sequencing data

                       Changes in version 1.11.2                        

- re-upload files under Version 1.11.1

                       Changes in version 1.11.1                        

- update reference manual with single-cell RNA sequencing data

[sechm](/packages/sechm)
-----

                 Changes in version 1.0.0 (2021-04-19)                  

- sechm moved from the SEtools package

- generalized annotations

[seq2pathway](/packages/seq2pathway)
-----------

                       Changes in version 1.23.1                        

- Upgraded datasets to hg38 and mm10

- Upgraded python script from python2 to python3

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.32.0                        

NEW FEATURES

- new option 'ret.idx' in `seqSetFilter()` for unsorted sample and
  variant
  indices

- new option 'ret.idx' in `seqSetFilterAnnotID()` for unsorted variant
  index

- rewrite the function `seqSetFilterPos()`: new options 'ref' and
  'alt',
  'multi.pos=TRUE' by default

- new option 'packed.idx' in `seqAddValue()` for packing an indexing
  variable

- new option 'warn' in `seqSetFilter()` to enable or disable the
  warning

- new functions `seqNewVarData()` and `seqListVarData()` for
  variable-length data

UTILITIES

- allow no variant in `seqApply()` and `seqBlockApply()`

- the list object returned from `seqGetData()` always have names if
  there
  are more than one input variable names

BUG FIXES

- `seqGDS2VCF()` should output "." instead of NA in the FILTER column

- `seqGetData()` should support factor when '.padNA=TRUE' or
  '.tolist=TRUE'

- fix `seqGDS2VCF()` with factor variables

- `seqSummary(gds, "$filter")` should return a data frame with zero row
  if
  'annotation/filter' is not a factor

[SeqGate](/packages/SeqGate)
-------

                 Changes in version 1.1.1 (2021-01-22)                  

- Updated vignette and documentation to fix typos

[SeqGSEA](/packages/SeqGSEA)
-------

                 Changes in version 1.31.1 (2020-11-30)                 

- Depends on DESeq2 (previously DESeq)

[SEtools](/packages/SEtools)
-------

                 Changes in version 1.5.3 (2021-04-19)                  

- plotting functions moved to the `sechm' package

[SharedObject](/packages/SharedObject)
------------

                        Changes in version 1.5.0                        

- NEW FEATURES:
  + Support sharing character vectors
  + Support creating empty shared objects
  + new paramters `shareAttributes` and `minLength` in package options
  + new paramter `depth` in function is.shared
  + R level developer API supports vector input
  + add a safer memory check method before sharing an object on Linux
  + Chinese vignette

- CHANGES:
  + The package will not provide C APIs in this version
  + The C++ APIs are not compatible with the old version
  + is.shared function reports the sharing status of object
  attributes by appending "Shared" to the end of attribute names.

[shinyepico](/packages/shinyepico)
----------

                 Changes in version 0.99.9 (2021-01-07)                 

- Bug fixed: DMR plot group colors

                 Changes in version 0.99.0 (2020-11-13)                 

- Submitted to Bioconductor

[ShortRead](/packages/ShortRead)
---------

                        Changes in version 1.50                         

NEW FEATURES

- (v 1.49.1) `as(., "QualityScaledDNAStringSet")` propagates
  names. See https://github.com/Bioconductor/ShortRead/issues/3.

- (v 1.49.1) implement`as(., "DNAStringSet")`. See
  https://github.com/Bioconductor/ShortRead/issues/3.

BUG FIXES

- (v 1.49.3) report invalid FastqStreamer / FastqSampler on
  re-serialized objects. Fixes
  https://github.com/Bioconductor/ShortRead/issues/5.

[signatureSearch](/packages/signatureSearch)
---------------

                 Changes in version 1.5.3 (2021-04-12)                  

- Created cellNtestPlot function to visualize number of compounds
tested in cell types along with primary site information
- Added get_treat_info function to get the treatment information in
reference database including pert, cell, pert_type columns.
- Supported setReadable function and readable argument in TSEA
functions to convert Entrez id to gene Symbols in the itemID column
in the enrichment result table.

- Supported dtnetplot on Reactome pathway

                 Changes in version 1.5.2 (2021-02-22)                  

- Supported 3 enrichment methods in TSEA on Reactome pathway

[SimFFPE](/packages/SimFFPE)
-------

                 Changes in version 1.3.1 (2020-11-11)                  

- Speed up

- Modified more defaults based on publicaly available FFPE samples

- Added new parameter sameChrProp for a better simulation

- Modified the model of SRCR distance for adjcent ss-DNA combination:
  from normal distribution to log-normal distribution

- Renaming of parameters

- Fixed the bug of producing chimeric reads without SRCR

- Fixed the bug of producing read duplicates at low coverage

[simplifyEnrichment](/packages/simplifyEnrichment)
------------------

                        Changes in version 1.1.4                        

- add `export_to_shiny_app()`

- add `simplifyGOFromMultipleLists()`

                        Changes in version 1.1.2                        

- add `anno_word_cloud()` function

[SingleCellExperiment](/packages/SingleCellExperiment)
--------------------

                       Changes in version 1.14.0                        

- Added the unsplitAltExps() function to reverse the effect of
  splitting alternative experiments.

- Added a mainExpName() getter and setter to remember the name of
  the main experiment.

- splitAltExps() will now store the chosen ref as the
  mainExpName.

- swapAltExp() now discards the promoted experiment from the list
  of alternative experiments in the output. It will also exchange
  the colData between the swapped experiments when
  withColData=TRUE. These changes assist in achieving
  reversibility of the output.

- Added applySCE() to conveniently apply a function to the main
  and alternative Experiments.

- Added withDimnames= to reducedDim<-() and reducedDims<-().  If
  TRUE, these methods now emit warnings on observing incompatible
  row names in value.

- Respect any metadata passed in with value in reducedDims<-()
  and altExps<-().

- Added the reduced.dim.matrix class to preserve attributes
  inside the reducedDims during subsetting/combining.

- Setting withColData=TRUE in altExp() and altExps() will now
  prepend colData(x) to the output colData.

- Added withDimnames= to altExp<-() and altExps<-().  If TRUE,
  these methods now emit warnings on observing incompatible
  column names in value. Also added withColData=, which will now
  reverse the prepending in the getter if the left-most columns
  are the same as colData(x). (If not the same, a warning is
  emitted.)

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.1.2 (2021-05-13)                  

- Added diffAbundanceFET and plotClusterAbundance function

- Linked Shiny UI help buttons to new online help pages

- Expanded convertSCEtoSeurat() function to copy additional data

- Updated and merged pkgdown documentation

- Added HTML reports for Seurat curated workflow and marker finding

- Refactor of Normalization UI

- Added tagging system for matrix type

- Several bug fixes

- Added generic wrapper functions for normalization, dimensionality
  reduction and feature selection

                 Changes in version 2.0.1 (2021-01-07)                  

- Added cell type labeling functional, wrapping SingleR method

- Added cell type labeling UI under differential expression tab

- Added marker identification in Seurat workflow

[SingleR](/packages/SingleR)
-------

                        Changes in version 1.6.0                        

- Relaxed the requirements for consistent row names in
  combineRecomputedResults().

- Support sparse DelayedArray inputs in classifySingleR().

- Parallelize over labels instead of rows in
  aggregateReference(), with minor changes in the setting of the
  seed. Restrict the PCA to the top 1000 most highly variable
  genes, for speed.

[sitadela](/packages/sitadela)
--------

                 Changes in version 0.0.1 (2021-02-07)                  

NEW FEATURES

- First release

[sitePath](/packages/sitePath)
--------

                        Changes in version 1.7.8                        

- Change default 'minSNP' value for 'parallelSites' function.

                        Changes in version 1.7.7                        

- Fix: Add 'rmarkdown' in 'Suggests'.

- Only plot paths with duplication in number for 'sneakPeek' function.

                        Changes in version 1.7.6                        

- Bug fix: empty groups produced by 'groupTips' function.

- Create 'paraFixSites' and 'fixationIndel' functions.

                        Changes in version 1.7.5                        

- Treat 'phyMSAmatched' object as 'lineagePath' class for simplicity.

- Improved multiprocessing.

- Fix: repeated cluster name by '.assignClusterNames' internal
  function.

- Rename 'allSitesPos' to 'allSitesName'.

- Add 'plotMutSites' support for 'lineagePath' and 'fixationSites'
  objects.

                        Changes in version 1.7.4                        

- 'cl.cores' option for turning multiprocessing on and off.

                        Changes in version 1.7.3                        

- Fix: inability to get position of all the sites.

                        Changes in version 1.7.2                        

- Fix: missing export for 'plot.phyMSAmatched' function.

- Fix: 'addMSA' function unable to handle 'treedata' object.

- Multiprocess for 'addMSA' and 'sitesMinEntropy' function.

                        Changes in version 1.7.1                        

- Guess sequence type based on ATCG proportion for 'addMSA'.


[slingshot](/packages/slingshot)
---------

                        Changes in version 2.0.0                        

- Added a NEWS.md file to track changes to the package.

- Changed default output of most functions from SlingshotDataSet to
PseudotimeOrdering and added conversion functions between them and
SingleCellExperiment.

- getLineages now relies on createClusterMST

- Removed plotGenePseudotime

- added as.df option to slingCurves and slingMST, which provide
relevant information for plotting as data.frame objects. This should
help those plotting Slingshot results with ggplot, especially for
our traffic package.

- updated all documentation.

[Spaniel](/packages/Spaniel)
-------

                        Changes in version 1.5.0                        

New Features

- Import options for 10X Visium data
- Updated vignette for 10X Visium data

Bug Fixes

- Corrected dependency issue that was causing build to fail

[SpatialExperiment](/packages/SpatialExperiment)
-----------------

                Changes in version 1.1.700 (2021-05-05)                 

- updating roles in DESCRIPTION file

- adding GeneExpression biocVie

- switching LazyData to false

- fixing return value in show method documentation

                Changes in version 1.1.434 (2021-26-02)                 

- fixing subset method according to SummarizedExperiment generic
  definition

                Changes in version 1.1.432 (2021-15-02)                 

- fixing documentation on latex errors

            Changes in version 1.1.430-1.1.431 (2021-12-02)             

- adding BumpyMatrix as suggests

- update documentation

                Changes in version 1.1.429 (2021-12-02)                 

- adding cd_keep = TRUE binds all the colData to the spatialData

- fixing bug for cd_keep with multiple elements

                Changes in version 1.1.428 (2021-12-02)                 

- adding itemize to assays vignette item

- correcting typo cd_keep->cd_bind in spatialData documentation

                Changes in version 1.1.427 (2021-09-02)                 

- fixing tenx vignette

                Changes in version 1.1.426 (2021-09-02)                 

- fixing read10xVisium example with data parameter

                Changes in version 1.1.425 (2021-08-02)                 

- restoring data parameter in read10xVisium

- missing itemize in combine documentation

                Changes in version 1.1.424 (2021-07-02)                 

- cleaning documentation

- removing spatialImage-methods.R file

                Changes in version 1.1.423 (2021-02-02)                 

- fixing documentation issues on imgData

                Changes in version 1.1.422 (2021-29-01)                 

- removing ve data because of local image problem
  (using example(read10xVisium) instead)

                Changes in version 1.1.421 (2021-28-01)                 

- fixing ve data local image problem

                 Changes in version 1.1.42 (2021-21-01)                 

- implementing new SpatialExperiment class

- spatialData slot

- imgData and image handling methods (HLC)

                 Changes in version 1.1.6 (2021-02-04)                  

- removed additional slots in the SPE class definition
  (i.e. spatialData and spatialCoordsNames)

- spatialCoords stored in int_colData() = numeric matrix

- spatialDataNames stored in int_metadata()
  = character vector specifying a subset of colData()

- the SPE show method does not include spatialData in colData;
  instead, spatialData/CoordNames are printed separately

- the SPE constructor now allows specification of
  spatialData/Coords/-Names
  where -Names can be a subset of the supplied colData();
  spatialData/Coords are thus optional

- cbind() now allows duplicated sample_ids,
  which are made unique with a message

- consistent usage of "spe" for SpatialExperiment objects across
  all examples (previously, both ve and se were used as well)

- fixed cache/path-related error on windows in SpatialImage unit tests

- added unit-tests of SpatialExperiment class validity

- imgData field in int_metadata is now required
  to exist (but can be an empty DFrame)

- colData<- protects sample_id & spatialDataNames
  fields; spatialData<- protects colData

                 Changes in version 1.1.5 (2021-31-03)                  

- version bump to x.y.z format with .z increment

- general code-style revision to keep to Bioc guidelines including,
  e.g.

- usage of accessors (and not @)

- keeping to a 80-character limit

- spaces around logical operators (but not function arguments)

- in-line { for function definitions, if-else statements etc.

- re-ordering of roxygen2 documentation to be consistent across scripts

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 1.1.6 (2021-05-17)                  

- The spatial heatmaps are able to maitain outline (stroke) widths
  defined in aSVG. The stroke widths can be updated with
  "update_feature".

- Text in aSVGs can be independent from features, i.e. have separate
  colors, stroke widths, etc.

- Only "g", "path", "rect", "ellipse", "use", and "title" elements are
  allowed in aSVG files, and other elements will raise errors or
  warnings. The "use" element should not in "g", and "g" elements
  should not have the "transform" attribute with a "matrix" value.

- Extraction of shape coordinates includes two alternative methods. If
  one fails, the other will be used by default. So in most cases though
  some coordinates are missing due to irregular shapes, spatial
  heatmaps can still be created.

- Added spatiotemporal example of rice coleoptile to Shiny app and
  vignette.

- Downloaded gene expression data are cached.

- In "spatial_hm", the "tis.trans" argument was replaced by "ft.trans";
  width and height arguments were removed, since the spatial heatmap
  aspect ratio is set the same with original aSVG.

- Shiny app: the function "shiny_all" was renamed to "shiny_shm"; the
  app is able to take HDF5 database backend, which contains data and
  aSVG files, and the HDF5 database can also be uploaded on the user
  interface; re-organized user interface: Landing Page (includes links
  to different app instances), Spatial Heatmap (includes sub-tabs of
  Image, Interactive, Video, Matrix Heatmap, Interactive Network),
  Spatial Enrichment (identifies spatial feature-specific genes),
  About; new functionality introduced: auto-completion search box, URLs
  of specific app states can be bookmarked, full screen, scrolling
  height, tooltip, one-to-multiple re-matching of spatial features,
  fixed aspect ratio (SHMs are not squeezed in the case of multiple
  aSVGs), metadata column and link column in the data matrix; code was
  organized in modules, etc.

[Spectra](/packages/Spectra)
-------

                         Changes in version 1.1                         

Changes in 1.1.20

- Fix concatenating empty spectra (issue #200).

Changes in 1.1.19

- New filterPrecursorCharge() method.

Changes in 1.1.18

- Define plotSpectraMirror as a method.

Changes in 1.1.17

- Fix issue #187.
- Add function concatenateSpectra to allow concatenating Spectra
objects and list of Spectra objects.

Changes in 1.1.16

- Support arbitrary spectra variables to be passed to the functions
provided/added with addProcessing; issue #182.

Changes in 1.1.15

- Pass spectras' precursor m/z to the MAPFUN in compareSpectra; issue
#171.
- Add joinPeaksGnps to perform a peak matching between spectra
similar
to the one performed in GNPS (issue #171).

Changes in 1.1.14

- Support plotting of empty spectra (issue 175).

Changes in 1.1.13

- Move ProcessingStep to ProtGenerics.

Changes in 1.1.12

- Fix show method for Spectra to list only the 3 most recent
processing steps (issue 173).
- Add processingLog function to display the log messages of all
processing steps of a Spectra object.

Changes in 1.1.11

- Add support for ... to pickPeaks and smooth (issue 168).

Changes in 1.1.10

- Import filterIntensity from ProtGenerics.

Changes in 1.1.9

- Fix label in plotSpectra.

Changes in 1.1.8

- filterIntensity supports passing of additional parameters to the
used filter function (issue 164).

Changes in 1.1.7

- Fix bug in show,ProcessingStep (issue 162).

Changes in 1.1.6

- New joinSpectraData() function.

Changes in 1.1.5

- Add [[,Msbackend and [[<-,MsBackend methods (issue 149).
- Add [[,Spectra and [[<-,Spectra methods.

Changes in 1.1.4

- Fix issue with labelCol in plotSpectra (issue #157).

Changes in 1.1.3

- Implement a generic Spectra,ANY constructor replacing
Spectra,DataFrame and Spectra,character.

Changes in 1.1.2

- Fix problem in export to mzML files that failed for empty spectra
(issue #145)

Changes in 1.1.1

- Round retention time in figure titles.
- Document differences between spectrumId (spectrumID),
acquisitionNum
and scanIndex.

Changes in 1.1.0

- New Bioc devel version

[SpectralTAD](/packages/SpectralTAD)
-----------

                 Changes in version 1.7.1 (2020-11-22)                  

- Fix bedpe output for Juicebox

[splatter](/packages/splatter)
--------

                 Changes in version 1.16.0 (2020-05-20)                 

- 
  Substantial updates to the splatPop simulation (from Christina
  Azodi)
  
  • Added ability to simulate data with complex multiplexed
  sequencing designs
  
  • Added simulation of “conditional” effects, where a subset
  of DE and eQTL effects are applied to only a subset of
  individuals (e.g. disease vs. healthy samples)
  
  • Added the ability to simulate different numbers of cells
  for each sample, sampled from a gamma distribution.
  
  • Updates to the splatPop vignette describing these changes

- 
  Logical matrices should now be handled correctly when
  minimising output SingleCellExperiment objects

- 
  Other minor fixes

[SplicingFactory](/packages/SplicingFactory)
---------------

                       Changes in version 0.99.18                       

- R version update

                       Changes in version 0.99.17                       

- Minor corrections

                       Changes in version 0.99.16                       

- Minor corrections

                       Changes in version 0.99.15                       

- Adjust sample dataset creating script to the use of
SummarizedExperiment

                       Changes in version 0.99.14                       

- Documentation updates

                       Changes in version 0.99.13                       

- SummarizedExperiment input for calculate_difference
- Changed SE_assay parameter name to assayno
- Documentation and example updates

                       Changes in version 0.99.12                       

- Documentation and example updates

                       Changes in version 0.99.11                       

- SummarizedExperiment corrections
- Vignette updates

                       Changes in version 0.99.10                       

- SummarizedExperiment output for calculate_diversity
- Vignette updates
- Example dataset updates

                       Changes in version 0.99.9                        

- Documentation and code formatting updates.

                       Changes in version 0.99.8                        

- Examples for diversity calculation functions.

                       Changes in version 0.99.7                        

- More flexible single calculate_entropy function instead of separate
naive and Laplace entropy.
- P-value correction method can be set by user.
- Formatting corrections.

                       Changes in version 0.99.5                        

- Updated DESCRIPTION.
- Formatting corrections.
- Corrected bug in data description.
- Updated methods for obtaining class of input.
- Added verbose argument for functions.

                       Changes in version 0.99.4                        

- Code and documentation formatting corrections.
- Updates to example dataset.

                       Changes in version 0.99.3                        

- SummarizedExperiment input type updated in calculate_diversity, new
argument: SE_assay.

- Documentation, vignette updated.

                       Changes in version 0.99.2                        

- SummarizedExperiment input type instead of ExpressionSet.

                       Changes in version 0.99.1                        

- Correction: unnecessary file removed.

                       Changes in version 0.99.0                        

- Submitted to Bioconductor.

[statTarget](/packages/statTarget)
----------

                       Changes in version 1.21.1                        

- To update for coCV output

[structToolbox](/packages/structToolbox)
-------------

                        Changes in version 1.4.0                        

- added VIP summary chart

- add ellipse plotting options to pca_scores_plot

- mv_sample_filter can be used in train/predict mode

- Documentation updates

- Minor issue fixes

- Add median method for fold change

- Update PQN with new inputs/outputs due to changes in pmp package

- Update SBC with new inputs/outputs due to changes in pmp package

- Update MTBLS79 due to changes in SBC

[Summix](/packages/Summix)
------

                 Changes in version 0.99.9 (2021-02-19)                 

- Fixed NEWS text formatting

                 Changes in version 0.99.8 (2021-02-11)                 

- Fixed multiple issues in review for development

- No longer using for loops in function code

[supersigs](/packages/supersigs)
---------

                       Changes in version 0.99.0                        

- Adhere to Bioconductor submission guidelines
- Add simplify_signature function


[SWATH2stats](/packages/SWATH2stats)
-----------

                       Changes in version 1.21.2                        

BUG FIXES

- add rmarkdown to DESCRIPTION

                       Changes in version 1.21.1                        

UPDATE

- update axis.text and theme_bw

                       Changes in version 1.21.0                        

NEW FEATURES

- SWATH2stats in BioC 3.13 development release

[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.3.15                        

- Changes to concensus score in PairSummaries.

                       Changes in version 1.3.14                        

- Major changes to the PairSummaries function and minor changes to
NucleotideOverlaps, ExtractBy, and FindSets. Adjustments to the
model that PairSummaries calls on to predict PIDs.

                       Changes in version 1.3.13                        

- ExtractBy function has been added. Allows extraction of feature
sequences into XStingSets organized by the a PairSummaries object or
the single linkage clusters implied by pairings within the
PairSummaries objects.
- DisjointSet function added to extract single linkage clusters from
a
PairSummaries object.

                       Changes in version 1.3.12                        

- PairSummaries now computes 4-mer distance between predicted pairs.

                       Changes in version 1.3.11                        

- PairSummaries now returns a column titled Adjacent that provides
the
number of directly adjacent neighbor pairs to a predicted pair. Gap
filling code adjusted.
- The function FindSets has been added and performs single linkage
clustering on a pairs list as represented by vectors of integers
using the Union-Find algorithm. Long term this function will have a
larger wrapper function for user ease of access but will remain
exposed.

                       Changes in version 1.3.10                        

- NucleotideOverlap now passes it's GeneCalls object forward,
allowing
PairSummaries to forego inclusion of that object as an argument.

                        Changes in version 1.3.9                        

- Minor vignette and suggested package changes.

                        Changes in version 1.3.8                        

- PairSummaries now allows users to fill in specific matching gaps in
blocks of predicted pairs with the arguments AllowGaps and
OffSetsAllowed.

                        Changes in version 1.3.7                        

- Adjustments to progress bars in both PairSummaries and
NucleotideOverlap.
- PID prediction models in PairSummaries adjusted.

                        Changes in version 1.3.6                        

- Contig name matching has been implemented. Scheme expects users to
follow NCBI contig naming and gff formats, accepting contig names
from gffs directly, and removing the first whitespace and everything
thereafter from FASTA headers. Contig name matching can be disabled
if users wish, using the argument AcceptContigNames, but ensuring
that the correct contigs in GeneCalls objects are matched to the
appropriate contigs in Synteny objects are then the user's
responsibility.

                        Changes in version 1.3.5                        

- PairSummaries now translates sequences based on transl_table
attributes provided by gene calls
- PairSummaries now uses a generic model for predicting PID
- gffToDataFrame now parses out the transl_table attribute

                        Changes in version 1.3.2                        

- Minor changes to NucleotideOverlap
- Major changes to PairSummaries - can now take in objects of class
Genes build by the DECIPHER function FindGenes()

[systemPipeShiny](/packages/systemPipeShiny)
---------------

                       Changes in version 1.1.40                        

Major change

- Add is_demo option: only affect workflow module right now. Lock
users inside a temp folder when using the WF module and give users a
new temp folder every time they refresh. This will prevent directory
exist problem if many users are using a same deploy instance.

- Add welcome_guide option: whether to enable the welcome guide which
highlights the guide dropdown menu.

- Rewrite welcome tab with a gallery to show all SPS features.

- loadDF, dynamicFile and dynamicFileServer added back to this
mainframe work package instead of spsComps, because these
dependencies have already been using in SPS. Leave these functions
in spsComps will introduce extra dependencies, and these functions
are not too frequently used outside the framework.

Minor change

- Option warning_toast now also checks if you are on "local" mode.

- Deleted some unwanted entries in reference generating yaml file.

- Fix some typos.

- More informative error message when the config file cannot be found
for spsOptions

- Add some .onLoad methods so users can use the spsOption to get
default values on package load.

- Updated essquise functions

- Add more guides.

- Removed the scroll to top button by shinyDashboardPlus, we have our
own "go top" button.

- Add assertions to spsInit.

- Add some screenshots to readme.

Bug fix

- Fix a bug when that loads the server twice
- Fix some default option values
- Fix a bug on addResourcePath when the working directory and app
directory is not the same.
- Fix links not working caused by website change
- Fix code in spsInit overwrite all current SPS options.
- Fix errors on admin page when server stats cannot be found, better
text and warning messages
- Fix new version of essquise introduced errors
- Fix a warning in vroom due to the column type problem

                       Changes in version 1.1.35                        

Major change

- Login feature added:

- Users can choose whether to enable the login or not in
global.RSPS options.

- There are also the login loading screen feature which can be
turned on and off.

- There are 3 different login loading screens right now and users
can interact with them.

- Website updated. https://systempipe.org/sps

- Updates on the admin panel:

- App information: added live charts for CPU, temperature, and RAM

- User control: admins now can add/delete/change users directly
from this tab, instead of only from command line.

Minor change

- Addtarget="_blank" to all external links in the app, so when they
are clicked, it will open in a new tab.

Bug fix

- FIx bugs due to login page caused server not loading

- Add 1s delay in javascript after login to resize the page so the
dashboard can be displayed normal.

- Fix a table rendering bug in workflow cwl step markdown text.

                       Changes in version 1.1.30                        

Major change

- new spsAccount class. This class is associated with login
management
, which allows users to create/delete user/admin accounts, change
password, change roles.

- Deprecated the spsPlotContainer class since we rewrite the Canvas
feature and move to a separate package {drawer}.

- New spsCoreTabReplace, which allows users to overwrite the default
core tabs.

- A lot more SPS options.

- Users can now choose whether to load or not load certain tabs on
start, even for default core tabs. Combining the
spsCoreTabReplace function, now users can customize everything
of the original app.

- Users can change the app title, and logo image.

- Admin panel added to app. Users now can visit the admin panel by
adding "?user_definded_string" to the end of the url. Default is
"admin". Login with an admin account is required. Users can use the
spsAccount class to add/change an admin account before starting the
app.

- App information: a tab displays current SPS app server
information, like CPU, RAM, size, etc.

- User control: a tab to see account information of current SPS
app.

- Changed the way to install modules. Default modules, workflow,
RNAseq and quick ggplot dependency packages are not installed by
default, unless you use dependency = TRUE in installation command.
It means all these dependencies are moved from Imports to the
Suggests class. This helps to save quite some time on SPS package
installation. Users install these packages based on their needs.
When users loads these modules but depend packages are not fully
installed, app will not crash, instead, a warning message with
install instructions will be displayed on both console and app UI.

- Based on the module installation change, module loading methods are
also changed. Module server functions are only called if users set
the option to load them. In previous versions, the server functions
are still loaded, just hide the unloaded module UI. This saves a lot
of time on app starting time, roughly from > 10s to < 3s if none of
the default modules are loaded.

Bug fix

- update all links to our new website: https://systempipe.org/sps

- Fix some bugs in the guide system

                       Changes in version 1.1.20                        

Major change

- 3 default modules complete: workflow, RNAseq,quick ggplot. Details
of these modules updated in our website:
<https://systempipe.org/sps>{.uri}.

- Separation of SPS smaller functions into 3 different packages. We
hope these packages can help people in their own Shiny app, or other
R projects.

- {spsComps}: SPS components, all new Shiny custom components and
utility functions that are used in Shiny server side.

- {drawer}: the redesign of Canvas, purely front-end image editing
tool.

- {spsUtil}: SPS utilities, general useful utility functions that
can be used with/without Shiny.

- Redesigned the new tab feature. Now users use spsNewTab function to
create their new custom visualization tab. The old newSpsTab
function is deprecated. Easier syntax and templates are used. By
default it will use the "simple" template which wraps 90% of the
shiny code from users so they can focus on the plotting code. There
is also the "full" template which expose all the Shiny code to
users.

- New spsEzUI and spsEzServer functions are used in the "simple"
template to wrap complex Shiny code.

- New spsOptDefaults, which prints out all the default SPS options
and
current values of these options on console.

- New notification system. Developers can write some notifications
which stores in a remote location and when app starts, it will try
to download and parse this file to notifications messages to
broadcast to users. This way, developers can send messages to users
often without re-deploy the app. The notification will appear on the
top right corner.

- The interactive guide is back. After a few versions of tests, we
added the guide system back. This time, developers can customize
their own guides. A guide_content.R file is created when a SPS
project initialize. It is stored in R of folder relate to the
project root. The guide will also be displayed on the app top right
corner.

Minor change

- updated all unit test to testthat v3 format.

Bug fix

- fix bugs due to shiny updates to 1.6.0

- Fix all bugs caused by {shinydashboardPlus} v2.0 updates.

                       Changes in version 1.1.10                        

Changes made from 1.1.0 to 1.1.05

Workflow module R session

- Now workflow module R session uses a background child R process,
which runs independently to the parent R session which runs shiny.
- So the shiny will be not blocked while code is running in the
background (you can still click other buttons when the child session
is busy) -- synchronous and non-blocking. A child indicator is also
placed in the UI, updates every second.
- The UI design of R session is similar to Rstudio. Four panels,
source code, console, log (specific to SPR), and plots.
- Standard out/error and plots are captured in the workflow folder.
Users can download them in the bundle on step 5 Workflow Run.
- Plots will be displayed on the plots panel. Now supports plots that
opens R device (base and ggplot), html widget plots are not
supported as this moment.
- A new shiny disconnection popup for SPS. Besides the gray layer on
shiny disconnection, a panel will be displayed to users to indicate
the problem. Similar to what shows on a shiny server, but more
informative and also works locally.
- Results of this session can be downloaded by closing the session
and
go back to step 5 of workflow module and there is a button to
download all in a zipped file.

RNAseq module

- redesigned the UI and server logic. Plots for DEG analysis and
Canvas connections.
- {SummarizedExperiment} supports. Now it returns
SummarizedExperiment
objects to global environment once the normalization or DEG
calculation is done.

General UI

- Added a "Go Top" button on the right bottom corner, clicking on
this
button will automatically scroll to the top of the page. This button
only shows up when client has > 50px scroll height.

Workflow module CWL tab

- Now the CWL file and CWL input file can be be edited. The edits
will
be imported to CWL parser every one second. Now this is a very
useful place to test or write new CWL scripts.
- Now this tab has a dynamically rendered dropdown panel which allows
users to choose which column for the targets table to map to the
variables in CWL input file.

Workflow module fully functioning

- Now you can run a full example workflow in SPS by choosing the
"Example" option on workflow setup step.

- Other systemPipeR preconfiged workflows will cause problems because
formatting issues that will cause errors in systemPipeR::runWF
function, not introduced by SPS. Please wait the updates on
systemPipeR to fix this. You can still use SPS to prepare files for
all workflows. That means, step 1-4 will work, step 5 will give you
errors if you choose a workflow which is not "Example".

Rework on the workflow part

- All 3 tabs merged into the main tab

- changed config tab to CWL tab

- added support for the running wf as a sub tab

- Now the main tab has 5 subtabs, they are all connected.

- Better guidelines for users, step-like usage, can't reach other
steps if a previous step is not completed.

- Original snapshot management drop down page changed to running
workflow session. This session will lock users to a unique page,
they can't interactive other app parts on the page(working directory
changed), to prevent app crash due to wd change.

Other changes

- A new UI component spsTimeline : horizontal timeline UI unit, has
status, can be updated on server by updateSpsTimeline.

- A new UI bsHoverPopover: enhanced high level function of
bsPlus::HoverPopover, additional JS used to make the popover work on
buttons as well.

- Fixed some link problems in renderDesc. Better links in renderDesc,
enlarged and spacing animation for links.

- Change on about page

- The news is now rendered on about tab in the app

- reduced developer content on about page.

- changed developer emails to github links.

Change on visualization

- RNAseq part is now only in one tab as big module: users upload the
targets file and a raw count table, and make different plots in
subtabs.

- This introduced a lot of dependencies, will decide later if we
keep as it is or separate it to spsBio.

[TargetSearch](/packages/TargetSearch)
------------

                       Changes in version 1.48.0                        

NEW FEATURES

- Function `ri_data_extract` allows for a time range for each searched
  m/z,
  instead of a single range for all masses.

- Man-pages typos and clarifications. No more user-significant changes.

BUG FIXES

- Add extra assertions on `ncdf4_convert`.

- The dependency package `ncdf4` should be on `Imports` rather than
  `Depends`
  on the DESCRIPTION file.


[TCC](/packages/TCC)
---

                       Changes in version 1.31.1                        

- changed default DE estimation method from exactTest to GLM-based test
  (glmQLFit, glmQLFTest) when using edgeR.

- removed DE estimation method for no-replicates dataset.

[TCGAutils](/packages/TCGAutils)
---------

                       Changes in version 1.12.0                        

New features

- makeSummarizedExperimentFromGISTIC has been moved to RTCGAToolbox.
- splitAssays now deprecated for TCGAsplitAssays to avoid conflict
with MultiAssayExperiment::splitAssays

Minor changes and bug fixes

- Properly identifies genome annotation (hg*) in oncoPrintTCGA
- qreduceTCGA now works with updates to seqlevelsStyle where genome
annotation include patch versions when available

[ternarynet](/packages/ternarynet)
----------

                       Changes in version 1.35.1                        

- added replica exchange MCMC parallel fitting algorithm

[TFutils](/packages/TFutils)
-------

                       Changes in version 1.11.1                        

- retrieve_lambert_main has revised URL for Cell supplemental xlsx
file

[TOAST](/packages/TOAST)
-----

                        Changes in version 1.5.1                        

- Correct F-test for multiple level testing.


[topdownr](/packages/topdownr)
--------

                        Changes in version 1.13                         

- New version for Bioc 3.13 (devel)

Changes in version 1.13.1

- as(..., "NCBSet") now treats neutral losses and modifications as
bonds as well.
- readTopDownFiles gains a new argument customModifications to allow
user-defined modifications. Suggestion and first implementation by
Maša Babović masab@bmb.sdu.dk [2021-03-15].

[ToxicoGx](/packages/ToxicoGx)
--------

                        Changes in version 1.1.0                        

- Continue to abstract functionality into CoreGx
- Add additional plotting functions such as grouped boxplots
- Extend coverage of unit tests to >90%
- Implement a faster version of drugPertubationSignature
- Add additional plotting functions
- Include scripts for differential expression analysis and GSEA of
toxico-genomic data (limma)

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.27.15                       

- Update documentation of geneModelFromTxdb

                       Changes in version 1.27.14                       

- Add rmarkdown into Suggests

- update importScSeqScore.

                       Changes in version 1.27.13                       

- Fix the size by number when read from file.

                       Changes in version 1.27.12                       

- add decontructor to hic.cpp.

                       Changes in version 1.27.10                       

- split the vignette to multiple files.

                       Changes in version 1.27.9                        

- plot back to back Interaction Data Track

                       Changes in version 1.27.8                        

- fix the typo in hic.cpp

                       Changes in version 1.27.7                        

- figure out the error in dyn.load trackViewer.so

                       Changes in version 1.27.6                        

- add support for .hic and .cool for importGInteraction

                       Changes in version 1.27.5                        

- add importGInteraction

                       Changes in version 1.27.4                        

- change the re-sample method for viewTracks

                       Changes in version 1.27.3                        

- Update documentation.

- add label_on_feature for lolliplot

                       Changes in version 1.27.2                        

- Fix the bug for pie plot of dandelion.plot when introduce
  label.parameters.

                       Changes in version 1.27.1                        

- Fix the bug that if all scores are greater than 10 and all scores are
  integer.

[tradeSeq](/packages/tradeSeq)
--------

                 Changes in version 1.5.02 (2021-01-21)                 

- Major update on associationTest, where the contrasts no longer rely
  on the knots but rather rely on a new nPoints argument, that
  specifies the number of points to use per lineage in the contrast.
  The associationTest also has a new argument contrastType that allows
  to use three different contrast types to do the test. See the docs on
  associationTest for more details.

[TrajectoryGeometry](/packages/TrajectoryGeometry)
------------------

                       Changes in version 0.99.12                       

Major changes

- Removed top level doc/ folder

                       Changes in version 0.99.11                       

Major changes

- depends R changed to >= 4.1

                       Changes in version 0.99.10                       

Major changes

- NEWS.md file added
- Formatting of SingleCellTrajectoryAnalysis.Rmd improved

[transomics2cytoscape](/packages/transomics2cytoscape)
--------------------

                        Changes in version 1.1.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Separation of only one function into two

- Introducing new input data formats to the above functions

[TraRe](/packages/TraRe)
-----

                       Changes in version 0.99.13                       

- No changes

                       Changes in version 0.99.12                       

- minor bugs fixed

                       Changes in version 0.99.11                       

- T/F changed by TRUE/FALSE

                       Changes in version 0.99.10                       

- More evaluated chunks in the vignettes

                       Changes in version 0.99.9                        

- xlsx files changed to html, library removed

- generate_cliques structure has changed, new plot function.

                       Changes in version 0.99.8                        

- Logfile added to runrewiring method.

- Message structure of preparerewiring method has been cleaned.

                       Changes in version 0.99.7                        

- Rewiring method outputs every supermodule

                       Changes in version 0.99.6                        

- vignette files are locally available (instead of github)

- SummarizedExperiment objects are allowed to use

- code style and structure has been improved

                       Changes in version 0.99.5                        

- openxlsx to xlsx library changed

                       Changes in version 0.99.4                        

- class() functions removed.

                       Changes in version 0.99.3                        

- excel_generation function added

- vignette file updated

                       Changes in version 0.99.2                        

- rewiring method part parallelized

                       Changes in version 0.99.1                        

- vignette file added

- minor bugs fixed

                       Changes in version 0.99.0                        

- automatic modules selection in rewiring test.

- rewiring html report file added.

- unit tests added.

[Travel](/packages/Travel)
------

                 Changes in version 0.99.0 (2020-12-30)                 

- Submitted to Bioconductor

[treeio](/packages/treeio)
------

                       Changes in version 1.15.6                        

- optimized read.nhx for large tree file (2021-03-12, Fri)

- https://github.com/YuLab-SMU/treeio/pull/51

                       Changes in version 1.15.5                        

- read.beast.newick and write.beast.newick for importing and
exporting
newick text with metadata in BEAST style (2021-03-11, Thu)
- https://github.com/YuLab-SMU/treeio/pull/50

                       Changes in version 1.15.4                        

- support parsing tree qza file from qiime2 (2020-03-01, Mon)
- https://github.com/YuLab-SMU/treeio/pull/46/files

                       Changes in version 1.15.3                        

- support parsing phyloxml (2021-02-04, Thu)
- https://github.com/YuLab-SMU/treeio/pull/44

                       Changes in version 1.15.2                        

- bug fixed of parsing nhx, now compatible with missing nhx tag
(2020-11-19, Thu)
- https://github.com/YuLab-SMU/treeio/pull/40

                       Changes in version 1.15.1                        

- remove magrittr::%<>% as it throw error of 'Error: could not find
function "%>%<-"' (2020-11-19, Thu)

[treekoR](/packages/treekoR)
-------

                       Changes in version 0.99.3                        

- First version of package submitted to BioConductor

[tricycle](/packages/tricycle)
--------

                        Changes in version 0.99                         

- Initial release.

- Added preprint to CITATION, vignette.

[TSCAN](/packages/TSCAN)
-----

                       Changes in version 1.30.0                        

- Migrated createClusterMST() to the TrajectoryUtils package.

- Modified orderCells() to return a more informative
  PseudotimeOrdering object.

- Handle pseudotime matrices in testPseudotime() by testing each
  path separately. Support inclusion of custom row.data= in each
  output DataFrame.

[ttgsea](/packages/ttgsea)
------

                 Changes in version 0.99.0 (2020-09-30)                 

- submission to Bioconductor

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.10.0                        

- Added more tximeta() messaging about specifying the
  'source' in linkedTxome. Essentially, this triggers
  GTF processing behavior that users may want to avoid,
  and so specifying a string other than "Ensembl" may be
  preferred. Also added note to vignette.

- Added note to vignette about alevin import with tximeta
  where the 'tgMap' step requires gene IDs and not gene
  symbols.

- Added hashes for:
  GENCODE 38 (H.s.), M27 (M.m), and Ensembl 104;
  GENCODE 37 (H.s.), M26 (M.m), and Ensembl 103;
  GENCODE 36 and Ensembl 102.

- Fixed a bug around AnnotationHub pulldown when using RefSeq
  as the source.
  * Fixed a bug where multiple parsed Ensembl GTF TxDb would be
  added to the BiocFileCache with the same rname.

                       Changes in version 1.9.11                        

- Added hashes for GENCODE 38 (H.s.), M27 (M.m), and
  Ensembl 104.

                        Changes in version 1.9.6                        

- Using tools::R_user_dir instead of rappdirs, in line with
  changes in BiocFileCache..

                        Changes in version 1.9.5                        

- Added hashes for GENCODE 37 (H.s.), M26 (M.m), and
  Ensembl 103.

                        Changes in version 1.9.2                        

- Added hashes for GENCODE 36 and Ensembl 102.

[tximport](/packages/tximport)
--------

                       Changes in version 1.19.4                        

- 'ignoreAfterBar' and txOut=TRUE will now strip the characters
  after '|' on the rownames of the output matrices.

[UMI4Cats](/packages/UMI4Cats)
--------

                        Changes in version 1.1.9                        

- Allow selection of number of reads to load from FastQ file in prep
  and split
  functions. Default: 1M reads (1e9).

- Minor fixes BiocCheck.

                        Changes in version 1.1.8                        

- Use `query_regions` to select the restriction fragments to use for
  differential
  testing in `diffWaldUMI4C`.

                        Changes in version 1.1.7                        

- Fixed bug with limits of the log2 OR values when plotting
  differential
  windows.

                        Changes in version 1.1.6                        

- Fixed bug in creation of gene annotation when an exon belongs to more
  than one transcript.

- Fixed bug when only one sample is provided (no name in assay
  columns).

                        Changes in version 1.1.5                        

- Fixed bug in selection of reference UMI4C sample when more than 1
  sample
  has the same number of total UMIs. Now it will select the first one.
  The selection of sample to use as reference can be overriden by the
  `ref_umi4c` argument.

                        Changes in version 1.1.4                        

- Fixed bug in domainogram plotting where white color was not aligned
  with 0 log2 FC.

                        Changes in version 1.1.3                        

- Fixed bug when providing a reference sample to use for normalizing
  UMI
  counts.

                        Changes in version 1.1.2                        

- Fixed bug when `cut_pos`!=0 that generated a gap in the digested
  genome
  object (see issue #8)

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.10.0                        

NEW FEATURES

- A new data structure, universalmotif_df, has been made available.
  This
  allows for motifs to be manipulated as one would a data.frame object.
  The to_df() function is used to generate this stucture from lists of
  motifs. The update_motifs() function is used to apply changes to the
  actual motifs, and to_list() returns the actual motifs. Note that
  this is
  only meant as an option for more conveniently manipulating motif
  slots of
  multiple motifs simultaneously before returning them to a list; the
  universalmotif_df structure cannot be used in the various
  universalmotif
  functions. Additionally, requires_update() can be used to ascertain
  whether motifs are out of date in a universalmotif_df object. Many
  thanks
  to @snystrom for discussions and significant contributions.

- view_motifs(): the universalmotif package now relies entirely on its
  own code to generate the polygon data used by ggplot2 to plot motifs,
  meaning the ggseqlogo import has been dropped. A number of new
  options
  are now available, including plotting multifreq logos and finer
  control
  over letter spacing. An effort has been made to ensure that the
  default
  behaviour of the function be unchanged from previous versions. This
  change should also allow for easier fixing of bugs and flexibility
  for
  future additions or changes.

- New function, merge_similar(): identify and merge similar motifs in a
  list of motifs. Essentially, a wrapper around compare_motifs(),
  hclust(),
  cutree(), and merge_motifs().

- New function, view_logo(): plot logos with matrix input instead of
  motif
  object input.  Arbitrary column heights and multi-character letters
  are
  allowed.

- New function, average_ic(): calculate the average information content
  for a list of motifs.

- trim_motifs(..., trim.from): trim from both directions or just one.

- shuffle_sequences(..., window, window.size, window.overlap): shuffle
  sequences iteratively in windows of specified size.

- scan_sequences(..., return.granges): optionally return a GRanges
  object.

- scan_sequences(..., no.overlaps, no.overlaps.by.strand,
  no.overlaps.strat): remove overlapping hits after scanning,
  preventing
  overlapping hits by the same motifs from being returned. This can
  optionally be done per strand. Either the first hit or the highest
  scoring hit can be preserved per set of overlapping hits. These new
  arguments can also be used in enrich_motifs().

- scan_sequences(..., respect.strand): whether to scan the sequence
  strands
  according to the motif strand slot. Only applicable for DNA/RNA
  motifs.
  This option is also available in enrich_motifs().

MINOR CHANGES

- Some additions and clean-up to documentation and vignettes.

- Support for MotIV-pwm2 formatted motifs has been dropped, as the
  package
  is no longer a part of the current Bioconductor version.

- read_matrix()/write_matrix(): the sep argument can now be NULL (no
  seperators.)

- The Rdpack dependency has been dropped.

- merge_motifs(): single-motif input now simply returns the motif
  instead
  of throwing an error.

- view_motifs(..., dedup.names): now TRUE by default. Furthermore, the
  make.unique() function is now used to deduplicate names.

- compare_motifs(..., method): the default comparison method has been
  changed back to PCC.

BUG FIXES

- Using create_motif() with a single character no longer throws an
  error.

- Generating random motifs with filled multifreq slots now works.


[VarCon](/packages/VarCon)
------

                 Changes in version 0.99.0 (2018-05-15)                 

- Submitted to Bioconductor

[VaSP](/packages/VaSP)
----

                 Changes in version 1.2.5 (2021-02-14)                  

- Changed package name to VaSP from vasp.

                 Changes in version 1.2.1 (2021-01-16)                  

- Added the citation and improved some codes.

[velociraptor](/packages/velociraptor)
------------

                        Changes in version 1.1.6                        

- Move sanity check vignette to inst/.

                        Changes in version 1.1.5                        

- Add Michael Stadler to package authors.

                        Changes in version 1.1.4                        

- Fix typo in documentation.

                        Changes in version 1.1.3                        

- Add vignette subdirectory with sanity checks.

                        Changes in version 1.1.2                        

- Add functions plotVelocity and plotVelocityStream.

                        Changes in version 1.1.1                        

- Refresh cached environments.

                        Changes in version 1.1.0                        

- Bioconductor release 1.1.0.

[vissE](/packages/vissE)
-----

                       Changes in version 0.99.0                        

- submitted to bioconductor

[wppi](/packages/wppi)
----

                 Changes in version 0.99.8 (2021-05-07)                 

- The workflow calculates Protein-Protein Interaction weights and
scores genes
- Database knowledge is automatically fetched from OmniPath, Gene
Ontology and Human Phenotype Ontology
- Submitted to Bioconductor

[xcms](/packages/xcms)
----

                       Changes in version 3.13.8                        

- Fix plotQC() for XCMSnExp objects

                       Changes in version 3.13.7                        

- Add `featureArea` function to extract the m/z-rt region for features.

- Fix `featureSpectra` function.

- Re-add the LC-MS/MS vignette.

- Feature: plotQC() supports XCMSnExp objects now

                       Changes in version 3.13.6                        

- Fix issue #545: skip second centWave run with CentWavePredIsoParam in
  regions
  of interest with undefined peak boundaries/scan ranges.

- Temporarily remove the LC-MS/MS vignette (until MsBackendMgf is added
  to
  Bioconductor).

                       Changes in version 3.13.5                        

- Add `filterChromPeaks` method to filter chromatographic peaks in a
  `XChromatogram` or `XChromatograms` object.

- Add `filterChromPeaks` method for `XCMSnExp` (issue #541).

- Support return of `Spectra` objects by `chromPeakSpectra`,
  `featureSpectra`
  and `reconstructChromPeakSpectra`.

- Support extraction of MS1 spectra with `chromPeakSpectra`.

- Support extraction of the spectrum with the largest total signal or
  largest
  base peak signal in `chromPeakSpectra`.

- Add support for extraction of spectra for selected/individual
  peaks/features
  using the `peaks` and `features` parameter in `chromPeakSpectra` and
  `featureSpectra`, respectively.

                       Changes in version 3.13.4                        

- Import `Param` object from `ProtGenerics`.

- Import `filterIntensity`, `normalize` and `alignRt` for
  `Chromatogram` and
  `MChromatograms` from `MSnbase`.

                       Changes in version 3.13.3                        

- `align,Chromatogram` gains new method `"none"` which will only keep
  values
  with identical retention times. For `method = "matchRtime"` the (much
  faster)
  matching function `closest` from the `MsCoreUtils` package is used.

- Method `correlate,Chromatogram` gains parameter `useIntensitiesAbove`
  to
  perform the correlation only with values larger than this threshold
  (avoiding thus high correlation because of many 0-values).

- Add method `filterIntensity,Chromatogram` that allows to filter a
  chromatogram
  object keeping only data points with an intensity above a user
  provided
  threshold.

                       Changes in version 3.13.2                        

- Add new function `manualChromPeaks` allowing to manually add and
  integrate
  chromatographic peaks.

                       Changes in version 3.13.1                        

- Support subsetting of `XChromatograms` with `drop = FALSE`.

[YAPSA](/packages/YAPSA)
-----

                       Changes in version 1.17.2                        

- We provide a new function LCD_extractCohort_callPerPID() which also
  belongs to
  the LCD family and which performs the detection of signatures at
  cohort-wide
  level, but re-runs the actual computation of the exposures per-PID
  with only
  the signatures identified in the cohort-wide calling. The ovall
  wrapper
  function LCD_complex_cutoff_combined() now also calls the new
  function and
  stores the result in the returned list with item name
  extractCohort_callPerPID

                       Changes in version 1.17.1                        

- Introduction of an input parameter minimumNumberOfAlterations for the
  functions LCD_complex_cutoff_perPID(), LCD_complex_cutoff_consensus()
  and
  LCD_complex_cutoff_combined(). If a sample has less mutations than
  this
  cutoff, a warning is issued. By default, this values is set to 25 and
  may
  be a good choice for analysis of SNV mutational signatures. For
  analysis of
  Indel mutational signatures, a better choice is 20.

[zellkonverter](/packages/zellkonverter)
-------------

                        Changes in version 1.2.0                        

- Update *anndata* and other Python dependencies, now using
  *anndata* v0.7.6

- Improved conversion checks for all slots in AnnData2SCE()

- Enable return conversion of the varm slot in AnnData2SCE()

- Avoid converting obsp and varp to dense matrices in
  AnnData2SCE()

- AnnData2SCE() should now always return dgCMatrix matrices when
  assays are sparse

- More consistent conversion of metadata to uns in SCE2AnnData()

- Handle conversion of list columns in colData and rowData in
  SCE2AnnData()

- Better support for converting *anndata* SparseDataset arrays

- Improved support for conversion of HDF5 backed AnnData objects

- Better support for writing DelayedArray assays in writeH5AD()

- Store X_name in AnnData2SCE() for use by SCE2AnnData() and add
  an X_name argument to AnnData2SCE() and readH5AD()

- Add a compression argument to writeH5AD()

- Export zellkonverterAnnDataEnv for use by other packages


NEWS from new and existing Data Experiment Packages
===================================


[BioImageDbs](/packages/BioImageDbs)
-----------

                       Changes in version 0.99.3                        

NEW FEATURES

- Package released

[chipenrich.data](/packages/chipenrich.data)
---------------

                       Changes in version 2.16.0                        

- Transition to Kai Wang as maintainer.

[curatedMetagenomicData](/packages/curatedMetagenomicData)
----------------------

                        Changes in version 3.0.0                        

- curatedMetagenomicData now contains 20,283 samples from 86 studies

- A total of 10,084 samples added since Bioconductor 3.10 (October 2019)

-  Studies added since Bioconductor 3.10 (October 2019):

    - AsnicarF_2021 (1098 samples)
    - BrooksB_2017 (408 samples)
    - ChuDM_2017 (86 samples)
    - DeFilippisF_2019 (97 samples)
    - GhensiP_2019 (113 samples)
    - GuptaA_2019 (60 samples)
    - HallAB_2017 (259 samples)
    - HMP_2019_ibdmdb (1628 samples)
    - HMP_2019_t2d (296 samples)
    - IjazUZ_2017 (94 samples)
    - KaurK_2020 (31 samples)
    - KeohaneDM_2020 (117 samples)
    - LassalleF_2017 (23 samples)
    - LifeLinesDeep_2016 (1135 samples)
    - LokmerA_2019 (57 samples)
    - MehtaRS_2018 (928 samples)
    - NagySzakalD_2017 (100 samples)
    - PasolliE_2019 (112 samples)
    - RosaBA_2018 (24 samples)
    - RubelMA_2020 (175 samples)
    - SankaranarayananK_2015 (37 samples)
    - ShaoY_2019 (1644 samples)
    - ThomasAM_2019_c (80 samples)
    - VilaAV_2018 (355 samples)
    - WampachL_2018 (63 samples)
    - WirbelJ_2018 (125 samples)
    - YachidaS_2019 (616 samples)
    - YassourM_2016 (36 samples)
    - YassourM_2018 (271 samples)
    - ZhuF_2020 (171 samples)
    
- All raw data has been reprocessed with MetaPhlAn3 & HUMAnN3
- The `curatedMetagenomicData()` method has been refactored for efficiency

    - It now returns SummarizedExperiment/TreeSummarizedExperiment objects
    - Sample metadata always stays up to date and is updated weekly
    - It is now the primary (and only) means to access data
    
- The `mergeData()` method has been refactored for accuracy and efficiency
- The `returnSamples()` method has been added for returns across studies
- The `sampleMetadata` object replaces the `combined_metadata` object
- The `combined_metadata` object will be removed in the next release
- A number of methods have moved directly to defunct status:

    - `cmdValidVersions()`
    - `getMetaphlanTree()`
    - `ExpressionSet2MRexperiment()`
    - `ExpressionSet2phyloseq()`
    
- All named accessors (e.g. `HMP_2012.pathcoverage.stool()`) have become defunct

    - These were very hard to maintain and document; the package is now simpler
    - The `curatedMetagenomicData()` method replaces all named accessors

[curatedTCGAData](/packages/curatedTCGAData)
---------------

                       Changes in version 1.14.0                        

New features

- The version argument now allows users to select either 1.1.38 or
  2.0.1.

- Version 2.0.1 includes RNASeq2Gene data as RSEM TPM gene expression
  values (#38, @mherberg).

- Genomic information updated for RaggedExperiment type data objects
  where '37' is now 'GRCh37' (#40, @vjcitn).

- Datasets (e.g., OV, GBM) that contain multiple assays that could be
  merged are now provided as merged assays (#27, @lwaldron).

- The vignette now includes sections on how to use the
  TCGAprimaryTumors and getWithColData functions.

Bug fixes and minor improvements

- mRNAArray assays now return matrix type data instead of DataFrame
  (#31, @lgeistlinger, @vjcitn).

[depmap](/packages/depmap)
------

                      Changes in version 1.5.1
                      
-  20Q4 data added for `crispr`, `copyNumber`, `TPM`, `mutationCalls` and `metadata` datasets. Newer versions for the other datasets were not released.


           
[DExMAdata](/packages/DExMAdata)
---------

                       Changes in version 0.99.0                        

- DExMAdata release.

[dorothea](/packages/dorothea)
--------

                 Changes in version 1.3.2 (2021-03-09)                  

- Added pancancer regulons for application in cancer.

                 Changes in version 1.3.1 (2021-02-08)                  

- Fixed bug in Seurat's related unit tests due to Seurats package
  update to version 4.0. s@assays$dorothea@misc is now list(), before
  it was NULL.

[emtdata](/packages/emtdata)
-------

                       Changes in version 0.99.0                        

- submitted to bioconductor

[GSE13015](/packages/GSE13015)
--------

                       Changes in version 0.99.11                       

- submission to Bioconductor

[HCAData](/packages/HCAData)
-------

                        Changes in version 1.8.0                        

Other notes

- The functions to load the data gain a as.sparse parameter, to
  control whether the underlying HDF5 dataset should be treated as
  sparse or not.

[imcdatasets](/packages/imcdatasets)
-----------

                 Changes in version 0.99.9 (2021-04-20)                 

- Added data from Zanotelli et al. Mol Syst Biol 16:e9798(2020)

                 Changes in version 0.99.8 (2021-04-15)                 

- Allow on disk storage of images and masks

                 Changes in version 0.99.7 (2021-03-24)                 

- Improved documentation

- pkgdown website

                 Changes in version 0.99.6 (2021-03-23)                 

- Added data from Jackson, Fischer et al. Nature 578,615–620(2020)

                 Changes in version 0.99.0 (2020-11-12)                 

- Extended vignette

- Added function documentation

- Added dataset documentation

- Formatted the package for Bioconductor submission

                 Changes in version 0.1.0 (2020-11-02)                  

- Initial commit

- Creation of the imcdatasets package

- Added the damond-pancreas-2019 dataset

[LRcellTypeMarkers](/packages/LRcellTypeMarkers)
-----------------

                       Changes in version 0.99.3                        

- Add human PBMC data

- Add LRcell related information in vignettes

                       Changes in version 0.99.2                        

- change dependency back to R >=4.1

- R>=3.6 triggers warning

                       Changes in version 0.99.1                        

- change dependency to R >=3.6

                       Changes in version 0.99.0                        

- version 0.99.0 released

- Submitted to Bioconductor

[microbiomeDataSets](/packages/microbiomeDataSets)
------------------

                 Changes in version 0.99.0 (2021-01-21)                 

- Submitted to Bioconductor

[MouseThymusAgeing](/packages/MouseThymusAgeing)
-----------------

                 Changes in version 0.99.5 (2021-04-21)                 

- Submitted to Bioconductor

[msigdb](/packages/msigdb)
------

                       Changes in version 0.99.0                        

- submitted to bioconductor

[ptairData](/packages/ptairData)
---------

                 Changes in version 0.99.8 (2021-04-08)                 

- ptairData Watched Tags added to Bioconductor Support Site User
  Profile

                 Changes in version 0.99.7 (2021-04-01)                 

- extended description

- added BugReports to description files

- completed NEWS file

- added inst/script/script.R

- added section installation and sessionInfo to the vignette

                 Changes in version 0.99.0 (2021-03-05)                 

- Submitted to Bioconductor

[RforProteomics](/packages/RforProteomics)
--------------

                       Changes in version 1.29.2                        

- Suggest rpx version 1.99.2 or later (to make use of caching and
  avoid repeated downloads).

- Remove the shinyMA function.

- Delete the code from the rTANDEM section, and only mention the
  package.

- Remove synapter(data) suggestion.

                       Changes in version 1.29.1                        

- Specify MSnID::peptides() (see also
  https://github.com/vladpetyuk/MSnID/issues/12).

[SBGNview.data](/packages/SBGNview.data)
-------------

                        Changes in version 1.5.1                        

- Major updates of SBGNview.data, by Kovidh Vegesna and Weijun Luo.

[scpdata](/packages/scpdata)
-------

                        Changes in version 0.99.3   
                        
- Updated documentation <2020-05-03>
- Added `liang2020_hela` datasets <2020-04-23>

                        Changes in version 0.99.2                        
                        
- Removed remaining tilde (U+223C) in man pages <2020-01-09>

                        Changes in version 0.99.1

- Removed tilde (U+223C) in man pages <2020-01-09>
                        
                        Changes in version 0.99.0         
                        
- Bioconductor submission <2020-01-06>

[scRNAseq](/packages/scRNAseq)
--------

                        Changes in version 2.6.0                        

- Added the Bacher T cell dataset.

- Added the Bhaduri organoid dataset.

- Added the Darmanis brain dataset.

- Added the Ernst spermatogenesis dataset.

- Added the Fletcher olfactory dataset.

- Added the Giladi HSC dataset.

- Added the He organ atlas dataset.

- Added the Jessa brain dataset.

- Added the Nowakowski cortex dataset.

- Added the Pollen glia dataset.

- Added the Zeisel nervous system dataset.

- Added the Zhao immune liver dataset.

- Added the Zhong prefrontal cortex dataset.

- Added the Bunis HSPC dataset (Dan Bunis).

[SimBenchData](/packages/SimBenchData)
------------

                 Changes in version 0.99.1 (2021-03-05)                 

- Submitted to Bioconductor

[SingleCellMultiModal](/packages/SingleCellMultiModal)
--------------------

                        Changes in version 1.4.0                        

New features

- SingleCellMultiModal function allows the combination of multiple
  multi-modal technologies.

- GTseq data from Macaulay et al. (2015) now available (@lgeistlinger)

- SCoPE2 data from Specht et al. now available thanks to @cvanderaa
  (#26)

- scMultiome provides PBMC from 10X Genomics thanks to @rargelaguet

Bug fixes and minor improvements

- Metadata information (function call and call to technology map)
  included in SingleCellMultiModal

- scNMT includes the original call in the MultiAssayExperiment
  metadata

- Improved and edited Contributing Guidelines for clarity

[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.3.19                        

SIGNIFICANT USER-VISIBLE CHANGES

- spatialLIBD has been updated to work with SpatialExperiment version
  1.1.701 which will be released as part of Bioconductor 3.13. This
  changes internal code of spatialLIBD which will work with any
  objects created with SpatialExperiment version 1.1.700.

                       Changes in version 1.3.16                        

SIGNIFICANT USER-VISIBLE CHANGES

- The citation information has changed now that spatialLIBD has a
  bioRxiv pre-print at
  https://www.biorxiv.org/content/10.1101/2021.04.29.440149v1.

                       Changes in version 1.3.15                        

SIGNIFICANT USER-VISIBLE CHANGES

- We now use plotly::toWebGL() to make the web application more
  responsive.

                       Changes in version 1.3.14                        

SIGNIFICANT USER-VISIBLE CHANGES

- The documentation and help messages shown in the web application
  have been revamped and improved.

                       Changes in version 1.3.12                        

NEW FEATURES

- We added a new vignette that shows how you can use spatialLIBD with
  any 10x Genomics Visium dataset processed with spaceranger. The
  vignette uses the publicly available human lymph node example from
  the 10x Genomics website.

                        Changes in version 1.3.3                        

NEW FEATURES

- Overall the package has been updated to use SpatialExperiment
  version 1.1.427 available on Bioconductor 3.13 (bioc-devel). Several
  functions were re-named such as sce_image_gene_p() now has a shorter
  name vis_gene_p(). This update also changes these visualization
  functions to ONLY support SpatialExperiment objects instead of the
  original modified SingleCellExperiment objects.

- Updated citation information to reflect that
  https://doi.org/10.1038/s41593-020-00787-0 is now public. Also added
  a link on the README to
  https://doi.org/10.6084/m9.figshare.13623902.v1 for the manuscript
  high resolution images.

[STexampleData](/packages/STexampleData)
-------------

                 Changes in version 0.99.0 (2021-03-28)                 

- Initial submission to Bioconductor

[TENxVisiumData](/packages/TENxVisiumData)
--------------

                       Changes in version 0.99.0                        

- initial package submission of 13 Visium spatial
  gene expression datasets by 10X Genomics

[TMExplorer](/packages/TMExplorer)
----------

                        Changes in version 1.1.2                        

- Data is now hosted on FigShare.

- Added new datasets: GSE150430, GSE154778, GSE125969, GSE134520,
  GSE123366

                        Changes in version 1.1.1                        

- Now uses BiocFileCache to download data.



NEWS from new and existing Workflows
===================================

[ExpHunterSuite](/packages/ExpHunterSuite)
--------------

                       Changes in version 0.99.11                       

- Improved documentation and code style <2021-02-21, Sun>


Deprecated and Defunct Packages
===============================


Sixty Five software packages were removed from this release (after being deprecated
in Bioc 3.12): 
adaptest, ArrayTV, BioSeqClass, CHARGE, chimera, CNVtools,
CorMut, DESeq, explorase, flowFit, flowSpy, flowType, focalCall, FourCSeq,
FunciSNP, GeneticsDesign, GenRank, GGBase, GGtools, GOFunction, gQTLBase,
gQTLstats, hicrep, ImpulseDE, ImpulseDE2, joda, JunctionSeq, LINC, Logolas,
mcaGUI, metaArray, metaseqR, methVisual, methyvim, Mirsynergy, MmPalateMiRNA,
MOFA, MotIV, NarrowPeaks, netbenchmark, netReg, OGSA, OmicsMarkeR, pathprint,
PathwaySplice, PGA, PGSEA, plrs, prada, Prize, Rariant, reb, Roleswitch,
rTANDEM, sapFinder, scsR, shinyTANDEM, sigaR, signet, simpleaffy,
spotSegmentation, Starr, SVAPLSseq, TxRegInfra, xps


Forty nine software are deprecated in this release and will be removed in Bioc 3.14:
AffyExpress, affyQCReport, AnnotationFuncs, ArrayTools, bigmemoryExtras,
BiocCaseStudies, CancerMutationAnalysis, CexoR, ChIPSeqSpike, CompGO, CoRegFlux,
CrossICC, cytofast, DBChIP, dexus, EasyqpcR, EDDA, eisa, ELBOW, ExpressionView,
FlowRepositoryR, genoset, HCABrowser, HCAExplorer, HCAMatrixBrowser, Imetagene,
IntramiRExploreR, mdgsa, metagenomeFeatures, methyAnalysis, MSEADbi, OutlierD,
pcot2, PCpheno, Polyfit, POST, RchyOptimyx, RDAVIDWebService, RNAither,
RNAprobR, rnaSeqMap, SAGx, samExploreR, seqplots, simulatorZ, SSPA, ToPASeq,
XBSeq, yaqcaffy


Fourteen experimental data packages were removed this release (after being
deprecated in BioC 3.12):
flowFitExampleData, FunciSNP.data, geuvPack, geuvStore2, GGdata, methyvimData,
mitoODEdata, Mulder2012, pathprintGEOData, pcaGoPromoter.Hs.hg19,
pcaGoPromoter.Mm.mm9, pcaGoPromoter.Rn.rn4, waveTilingData, yriMulti

Eleven experimental data packages are deprecated in this release and will be
removed in Bioc 3.14:
ceu1kg, ceu1kgv, ceuhm3, cgdv17, dsQTL, facsDorit, gskb, hmyriB36, JctSeqData,
MAQCsubsetAFX, yri1kgv


Fourteen annotation packages were removed from this release (after being deprecated
in Bioc 3.12): 
hom.At.inp.db, hom.Ce.inp.db, hom.Dm.inp.db, hom.Dr.inp.db,
hom.Hs.inp.db, hom.Mm.inp.db, hom.Rn.inp.db, hom.Sc.inp.db, KEGG.db,
MeSH.Eco.55989.eg.db, MeSH.Eco.ED1a.eg.db, MeSH.Eco.IAI39.eg.db,
MeSH.Eco.UMN026.eg.db, MeSH.Eqc.eg.db

Eighty seven annotation packages are deprecated in this release and will be 
removed in Bioc 3.14:
12 LRBase.XXX.eg.db packages (replaced with AHLRBaseDbs), 
MafDb.gnomAD.r3.0.GRCh38, MafH5.gnomAD.r3.0.GRCh38, 73 MeSH.XXX.eg.db packages 
(replaced with AHMeSHDbs)

No workflow packages were removed in this release.

One workflow package is deprecated in this release to be removed in 3.14:
eQTL
