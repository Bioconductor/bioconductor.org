April 16, 2025

**Bioconductor:**

We are pleased to announce Bioconductor 3.21, consisting of
2341 software packages, 432 experiment data packages, 928 annotation
packages, 30 workflows and 5 books.

There are 72 new software packages, 3 new data experiment packages,
no new annotation packages, no new workflows, no new books, and many updates and
improvements to existing packages.

Bioconductor 3.21 is compatible with R 4.5, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher, macOS arm64 and Linux
arm64. This release will also include updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.21](#getting-started-with-bioconductor-321)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.21
======================================

To update to or install Bioconductor 3.21

1. Install R 4.5. Bioconductor 3.21 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 72 new software packages in this release of Bioconductor.

- [alabaster.sfe](/packages/alabaster.sfe) Builds upon the existing
  ArtifactDB project, expending alabaster.spatial for language
  agnostic on disk serialization of SpatialFeatureExperiment.

- [barbieQ](/packages/barbieQ) The barbieQ package provides a series
  of robust statistical tools for analysing barcode count data
  generated from cell clonal tracking (i.e., lineage tracing)
  experiments. In these experiments, an initial cell and its
  offspring collectively form a clone (i.e., lineage). A unique
  barcode sequence, incorporated into the DNA of the inital cell, is
  inherited within the clone. This one-to-one mapping of barcodes to
  clones enables clonal tracking of their behaviors. By counting
  barcodes, researchers can quantify the population abundance of
  individual clones under specific experimental perturbations.
  barbieQ supports barcode count data preprocessing, statistical
  testing, and visualization.

- [BatchSVG](/packages/BatchSVG) `BatchSVG` is a feature-based
  Quality Control (QC) to identify SVGs on spatial transcriptomics
  data with specific types of batch effect. Regarding to the spatial
  transcriptomics data experiments, the batch can be defined as
  "sample", "sex", and etc.The `BatchSVG` method is based on binomial
  deviance model (Townes et al, 2019) and applies cutoffs based on
  the number of standard deviation (nSD) of relative change in
  deviance and rank difference as the data-driven thresholding
  approach to detect the batch-biased outliers.

- [beachmat.tiledb](/packages/beachmat.tiledb) Extends beachmat to
  initialize tatami matrices from TileDB-backed arrays. This allows
  C++ code in downstream packages to directly call the TileDB C/C++
  library to access array data, without the need for block processing
  via DelayedArray. Developers only need to import this package to
  automatically extend the capabilities of beachmat::initializeCpp to
  TileDBArray instances.

- [bedbaser](/packages/bedbaser) A client for BEDbase. bedbaser
  provides access to the API at api.bedbase.org. It also includes
  convenience functions to import BED files into GRanges objects and
  BEDsets into GRangesLists.

- [biocmake](/packages/biocmake) Manages the installation of CMake
  for building Bioconductor packages. This avoids the need for
  end-users to manually install CMake on their system. No action is
  performed if a suitable version of CMake is already available.

- [BreastSubtypeR](/packages/BreastSubtypeR) BreastSubtypeR is an R
  package that provides a collection of methods for intrinsic
  molecular subtyping of breast cancer. It includes subtyping methods
  for nearest centroid-based subtyping (NC-based) and single sample
  predictor (SSP-based), along with tools for integrating clinical
  data and visualizing results.

- [BulkSignalR](/packages/BulkSignalR) Inference of ligand-receptor
  (LR) interactions from bulk expression (transcriptomics/proteomics)
  data, or spatial transcriptomics. BulkSignalR bases its inferences
  on the LRdb database included in our other package,
  SingleCellSignalR available from Bioconductor. It relies on a
  statistical model that is specific to bulk data sets. Different
  visualization and data summary functions are proposed to help
  navigating prediction results.

- [CARDspa](/packages/CARDspa) CARD is a reference-based
  deconvolution method that estimates cell type composition in
  spatial transcriptomics based on cell type specific expression
  information obtained from a reference scRNA-seq data. A key feature
  of CARD is its ability to accommodate spatial correlation in the
  cell type composition across tissue locations, enabling accurate
  and spatially informed cell type deconvolution as well as refined
  spatial map construction. CARD relies on an efficient optimization
  algorithm for constrained maximum likelihood estimation and is
  scalable to spatial transcriptomics with tens of thousands of
  spatial locations and tens of thousands of genes.

- [CCAFE](/packages/CCAFE) Functions to reconstruct case and control
  AFs from summary statistics. One function uses OR, NCase, NControl,
  and SE(log(OR)). The second function uses OR, NCase, NControl, and
  AF for the whole sample.

- [chevreulPlot](/packages/chevreulPlot) Tools for plotting
  SingleCellExperiment objects in the chevreulPlot package. Includes
  functions for analysis and visualization of single-cell data.
  Supported by NIH grants R01CA137124 and R01EY026661 to David
  Cobrinik.

- [chevreulProcess](/packages/chevreulProcess) Tools analyzing
  SingleCellExperiment objects as projects. for input into the
  Chevreul app downstream. Includes functions for analysis of single
  cell RNA sequencing data. Supported by NIH grants R01CA137124 and
  R01EY026661 to David Cobrinik.

- [chevreulShiny](/packages/chevreulShiny) Tools for managing
  SingleCellExperiment objects as projects. Includes functions for
  analysis and visualization of single-cell data. Also included is a
  shiny app for visualization of pre-processed scRNA data. Supported
  by NIH grants R01CA137124 and R01EY026661 to David Cobrinik.

- [clustSIGNAL](/packages/clustSIGNAL) clustSIGNAL: clustering of
  Spatially Informed Gene expression with Neighbourhood Adapted
  Learning. A tool for adaptively smoothing and clustering gene
  expression data. clustSIGNAL uses entropy to measure heterogeneity
  of cell neighbourhoods and performs a weighted, adaptive smoothing,
  where homogeneous neighbourhoods are smoothed more and
  heterogeneous neighbourhoods are smoothed less. This not only
  overcomes data sparsity but also incorporates spatial context into
  the gene expression data. The resulting smoothed gene expression
  data is used for clustering and could be used for other downstream
  analyses.

- [CPSM](/packages/CPSM) The CPSM package provides a comprehensive
  computational pipeline for predicting the survival probability of
  cancer patients. It offers a series of steps including data
  processing, splitting data into training and test subsets, and
  normalization of data. The package enables the selection of
  significant features based on univariate survival analysis and
  generates a LASSO prognostic index score. It supports the
  development of predictive models for survival probability using
  various features and provides visualization tools to draw survival
  curves based on predicted survival probabilities. Additionally, SPM
  includes functionalities for generating bar plots that depict the
  predicted mean and median survival times of patients, making it a
  versatile tool for survival analysis in cancer research.

- [crumblr](/packages/crumblr) Crumblr enables analysis of count
  ratio data using precision weighted linear (mixed) models.  It uses
  an asymptotic normal approximation of the variance following the
  centered log ration transform (CLR) that is widely used in
  compositional data analysis.  Crumblr provides a fast, flexible
  alternative to GLMs and GLMM's while retaining high power and
  controlling the false positive rate.

- [crupR](/packages/crupR) An R package that offers a workflow to
  predict condition-specific enhancers from ChIP-seq data. The
  prediction of regulatory units is done in four main steps: Step 1 -
  the normalization of the ChIP-seq counts. Step 2 - the prediction
  of active enhancers binwise on the whole genome. Step 3 - the
  condition-specific clustering of the putative active enhancers.
  Step 4 - the detection of possible target genes of the
  condition-specific clusters using RNA-seq counts.

- [dandelionR](/packages/dandelionR) dandelionR is an R package for
  performing single-cell immune repertoire trajectory analysis, based
  on the original python implementation. It provides the necessary
  functions to interface with scRepertoire and a custom
  implementation of an absorbing Markov chain for pseudotime
  inference, inspired by the Palantir Python package.

- [DeconvoBuddies](/packages/DeconvoBuddies) Funtions helpful for
  LIBD deconvolution project. Includes tools for marker finding with
  mean ratio, expression plotting, and plotting deconvolution
  results. Working to include DLPFC datasets.

- [DNAcycP2](/packages/DNAcycP2) This package performs prediction of
  intrinsic cyclizability of of every 50-bp subsequence in a DNA
  sequence. The input could be a file either in FASTA or text format.
  The output will be the C-score, the estimated intrinsic
  cyclizability score for each 50 bp sequences in each entry of the
  sequence set.

- [ELViS](/packages/ELViS) Base-resolution copy number analysis of
  viral genome. Utilizes base-resolution read depth data over viral
  genome to find copy number segments with two-dimensional
  segmentation approach. Provides publish-ready figures, including
  histograms of read depths, coverage line plots over viral genome
  annotated with copy number change events and viral genes, and
  heatmaps showing multiple types of data with integrative clustering
  of samples.

- [G4SNVHunter](/packages/G4SNVHunter) G-quadruplexes (G4s) are
  unique nucleic acid secondary structures predominantly found in
  guanine-rich regions and have been shown to be involved in various
  biological regulatory processes. G4SNVHunter is an R package
  designed to rapidly identify genomic sequences with G4-forming
  potential and accurately screen user-provided single nucleotide
  variants (also applicable to single nucleotide polymorphisms) that
  may destabilize these structures. This enables users to screen key
  variants for further experimental study, investigating how these
  variants may influence biological functions, such as gene
  regulation, by altering G4 formation.

- [geyser](/packages/geyser) Lightweight Expression displaYer
  (plotter / viewer) of SummarizedExperiment object in R. This
  package provides a quick and easy Shiny-based GUI to empower a user
  to use a SummarizedExperiment object to view (gene) expression
  grouped from the sample metadata columns (in the `colData` slot).
  Feature expression can either be viewed with a box plot or a
  heatmap.

- [h5mread](/packages/h5mread) The main function in the h5mread
  package is h5mread(), which allows reading arbitrary data from an
  HDF5 dataset into R, similarly to what the h5read() function from
  the rhdf5 package does. In the case of h5mread(), the
  implementation has been optimized to make it as fast and
  memory-efficient as possible.

- [HiCParser](/packages/HiCParser) This package is a parser to import
  HiC data into R. It accepts several type of data: tabular files,
  Cooler `.cool` or `.mcool` files, Juicer `.hic` files or HiC-Pro
  `.matrix` and `.bed` files. The HiC data can be several files, for
  several replicates and conditions. The data is formated in an
  InteractionSet object.

- [imageTCGA](/packages/imageTCGA) A Shiny application to explore the
  TCGA Diagnostic Image Database.

- [islify](/packages/islify) This software is meant to be used for
  classification of images of cell-based assays for neuronal surface
  autoantibody detection or similar techniques. It takes imaging
  files as input and creates a composite score from these, that for
  example can be used to classify samples as negative or positive for
  a certain antibody-specificity. The reason for its name is that I
  during its creation have thought about the individual picture as an
  archielago where we with different filters control the water level
  as well as ground characteristica, thereby finding islands of
  interest.

- [jazzPanda](/packages/jazzPanda) This package contains the function
  to find marker genes for image-based spatial transcriptomics data.
  There are functions to create spatial vectors from the cell and
  transcript coordiantes, which are passed as inputs to find marker
  genes. Marker genes are detected for every cluster by two
  approaches. The first approach is by permtuation testing, which is
  implmented in parallel for finding marker genes for one sample
  study. The other approach is to build a linear model for every
  gene. This approach can account for multiple samples and backgound
  noise.

- [Lheuristic](/packages/Lheuristic) The Lheuristic package
  identifies scatterpots that follow and L-shaped, negative
  distribution. It can be used to identify genes regulated by
  methylation by integration of an expression and a methylation
  array. The package uses two different methods to detect expression
  and methyaltion L- shapped scatterplots. The parameters can be
  changed to detect other scatterplot patterns.

- [limpa](/packages/limpa) Quantification and differential analysis
  of mass-spectrometry proteomics data, with probabilistic recovery
  of information from missing values. Estimates the detection
  probability curve (DPC), which relates the probability of
  successful detection to the underlying expression level of each
  peptide, and uses it to incorporate peptide missing values into
  protein quantification and into subsequent differential expression
  analyses. The package produces objects suitable for downstream
  analysis in limma. The package accepts peptide-level data with
  missing values and produces complete protein quantifications
  without missing values. The uncertainty introduced by missing value
  imputation is propagated through to the limma analyses using
  variance modeling and precision weights. The package name "limpa"
  is an acronym for "Linear Models for Proteomics Data".

- [LimROTS](/packages/LimROTS) Differential expression analysis is a
  prevalent method utilised in the examination of diverse biological
  data. The reproducibility-optimized test statistic (ROTS) modifies
  a t-statistic based on the data's intrinsic characteristics and
  ranks features according to their statistical significance for
  differential expression between two or more groups (f-statistic).
  Focussing on proteomics and metabolomics, the current ROTS
  implementation cannot account for technical or biological
  covariates such as MS batches or gender differences among the
  samples. Consequently, we developed LimROTS, which employs a
  reproducibility-optimized test statistic utilising the limma
  methodology to simulate complex experimental designs. LimROTS is a
  hybrid method integrating empirical bayes and
  reproducibility-optimized statistics for robust analysis of
  proteomics and metabolomics data.

- [maaslin3](/packages/maaslin3) MaAsLin 3 refines and extends
  generalized multivariate linear models for meta-omicron association
  discovery. It finds abundance and prevalence associations between
  microbiome meta-omics features and complex metadata in
  population-scale epidemiological studies. The software includes
  multiple analysis methods (including support for multiple
  covariates, repeated measures, and ordered predictors), filtering,
  normalization, and transform options to customize analysis for your
  specific study.

- [MetaboDynamics](/packages/MetaboDynamics) MetaboDynamics is an
  R-package that provides a framework of probabilistic models to
  analyze longitudinal metabolomics data. It enables robust
  estimation of mean concentrations despite varying spread between
  timepoints and reports differences between timepoints as well as
  metabolite specific dynamics profiles that can be used for
  identifying "dynamics clusters" of metabolites of similar dynamics.
  Provides probabilistic over-representation analysis of KEGG
  functional modules and pathways as well as comparison between
  clusters of different experimental conditions.

- [miaDash](/packages/miaDash) miaDash provides a Graphical User
  Interface for the exploration of microbiome data. This way, no
  knowledge of programming is required to perform analyses. Datasets
  can be imported, manipulated, analysed and visualised with a
  user-friendly interface.

- [mist](/packages/mist) mist (Methylation Inference for Single-cell
  along Trajectory) is a hierarchical Bayesian framework for modeling
  DNA methylation trajectories and performing differential
  methylation (DM) analysis in single-cell DNA methylation (scDNAm)
  data. It estimates developmental-stage-specific variations,
  identifies genomic features with drastic changes along pseudotime,
  and, for two phenotypic groups, detects features with distinct
  temporal methylation patterns. mist uses Gibbs sampling to estimate
  parameters for temporal changes and stage-specific variations.

- [mitology](/packages/mitology) mitology allows to study the
  mitochondrial activity throught high-throughput RNA-seq data. It is
  based on a collection of genes whose proteins localize in to the
  mitochondria. From these, mitology provides a reorganization of the
  pathways related to mitochondria activity from Reactome and Gene
  Ontology. Further a ready-to-use implementation of MitoCarta3.0
  pathways is included.

- [MotifPeeker](/packages/MotifPeeker) MotifPeeker is used to compare
  and analyse datasets from epigenomic profiling methods with motif
  enrichment as the key benchmark.  The package outputs an HTML
  report consisting of three sections: (1. General Metrics) Overview
  of peaks-related general metrics for the datasets (FRiP scores,
  peak widths and motif-summit distances).  (2. Known Motif
  Enrichment Analysis) Statistics for the frequency of user-provided
  motifs enriched in the datasets.  (3. De-Novo Motif Enrichment
  Analysis) Statistics for the frequency of de-novo discovered motifs
  enriched in the datasets and compared with known motifs.

- [mspms](/packages/mspms) This package provides functions for the
  analysis of data generated by the multiplex substrate profiling by
  mass spectrometry for proteases (MSP-MS) method. Data exported from
  upstream proteomics software is accepted as input and subsequently
  processed for analysis. Tools for statistical analysis,
  visualization, and interpretation of the data are provided.

- [MSstatsBioNet](/packages/MSstatsBioNet) A set of tools for network
  analysis using mass spectrometry-based proteomics data and network
  databases. The package takes as input the output of MSstats
  differential abundance analysis and provides functions to perform
  enrichment analysis and visualization in the context of prior
  knowledge from past literature. Notably, this package integrates
  with INDRA, which is a database of biological networks extracted
  from the literature using text mining techniques.

- [OSTA.data](/packages/OSTA.data) 'OSTA.data' is a companion package
  for the "Orchestrating Spatial Transcriptomics Analysis" (OSTA)
  with Bioconductor online book. Throughout OSTA, we rely on a set of
  publicly available datasets that cover different sequencing- and
  imaging-based platforms, such as Visium, Visium HD, Xenium (10x
  Genomics) and CosMx (NanoString). In addition, we rely on scRNA-seq
  (Chromium) data for tasks, e.g., spot deconvolution and label
  transfer (i.e., supervised clustering). These data been deposited
  in an Open Storage Framework (OSF) repository, and can be queried
  and downloaded using functions from the 'osfr' package. For
  convenience, we have implemented 'OSTA.data' to query and retrieve
  data from our OSF node, and cache retrieved Zip archives using
  'BiocFileCache'.

- [pathMED](/packages/pathMED) PathMED is a collection of tools to
  facilitate precision medicine studies with omics data (e.g.
  transcriptomics). Among its funcionalities, genesets scores for
  individual samples may be calculated with several methods. These
  scores may be used to train machine learning models and to predict
  clinical features on new data. For this, several machine learning
  methods are evaluated in order to select the best method based on
  internal validation and to tune the hyperparameters. Performance
  metrics and a ready-to-use model to predict the outcomes for new
  patients are returned.

- [PICB](/packages/PICB) piRNAs (short for PIWI-interacting RNAs) and
  their PIWI protein partners play a key role in fertility and
  maintaining genome integrity by restricting mobile genetic elements
  (transposons) in germ cells. piRNAs originate from genomic regions
  known as piRNA clusters. The piRNA Cluster Builder (PICB) is a
  versatile toolkit designed to identify genomic regions with a high
  density of piRNAs. It constructs piRNA clusters through a stepwise
  integration of unique and multimapping piRNAs and offers
  wide-ranging parameter settings, supported by an optimization
  function that allows users to test different parameter combinations
  to tailor the analysis to their specific piRNA system. The output
  includes extensive metadata columns, enabling researchers to rank
  clusters and extract cluster characteristics.

- [poem](/packages/poem) This package provides a comprehensive set of
  external and internal evaluation metrics. It includes metrics for
  assessing partitions or fuzzy partitions derived from clustering
  results, as well as for evaluating subpopulation identification
  results within embeddings or graph representations. Additionally,
  it provides metrics for comparing spatial domain detection results
  against ground truth labels, and tools for visualizing spatial
  errors.

- [Polytect](/packages/Polytect) Polytect is an advanced
  computational tool designed for the analysis of multi-color digital
  PCR data. It provides automatic clustering and labeling of
  partitions into distinct groups based on clusters first identified
  by the flowPeaks algorithm. Polytect is particularly useful for
  researchers in molecular biology and bioinformatics, enabling them
  to gain deeper insights into their experimental results through
  precise partition classification and data visualization.

- [QRscore](/packages/QRscore) In genomics, differential analysis
  enables the discovery of groups of genes implicating important
  biological processes such as cell differentiation and aging.
  Non-parametric tests of differential gene expression usually detect
  shifts in centrality (such as mean or median), and therefore suffer
  from diminished power against alternative hypotheses characterized
  by shifts in spread (such as variance). This package provides a
  flexible family of non-parametric two-sample tests and K-sample
  tests, which is based on theoretical work around non-parametric
  tests, spacing statistics and local asymptotic normality
  (Erdmann-Pham et al., 2022+ [arXiv:2008.06664v2]; Erdmann-Pham,
  2023+ [arXiv:2209.14235v2]).

- [RbowtieCuda](/packages/RbowtieCuda) This package provides an R
  wrapper for the popular Bowtie2 sequencing read aligner, optimized
  to run on NVIDIA graphics cards. It includes wrapper functions that
  enable both genome indexing and alignment to the generated indexes,
  ensuring high performance and ease of use within the R environment.

- [ReducedExperiment](/packages/ReducedExperiment) Provides
  SummarizedExperiment-like containers for storing and manipulating
  dimensionally-reduced assay data. The ReducedExperiment classes
  allow users to simultaneously manipulate their original dataset and
  their decomposed data, in addition to other method-specific outputs
  like feature loadings. Implements utilities and specialised classes
  for the application of stabilised independent component analysis
  (sICA) and weighted gene correlation network analysis (WGCNA).

- [RFLOMICS](/packages/RFLOMICS) R-package with shiny interface,
  provides a framework for the analysis of transcriptomics,
  proteomics and/or metabolomics data. The interface offers a guided
  experience for the user, from the definition of the experimental
  design to the integration of several omics table together. A report
  can be generated with all settings and analysis results.

- [Rigraphlib](/packages/Rigraphlib) Vendors the igraph C source code
  and builds it into a static library. Other Bioconductor packages
  can link to libigraph.a in their own C/C++ code. This is intended
  for packages wrapping C/C++ libraries that depend on the igraph C
  library and cannot be easily adapted to use the igraph R package.

- [rigvf](/packages/rigvf) The IGVF Catalog provides data on the
  impact of genomic variants on function. The `rigvf` package
  provides an interface to the IGVF Catalog, allowing easy
  integration with Bioconductor resources.

- [RUCova](/packages/RUCova) Mass cytometry enables the simultaneous
  measurement of dozens of protein markers at the single-cell level,
  producing high dimensional datasets that provide deep insights into
  cellular heterogeneity and function. However, these datasets often
  contain unwanted covariance introduced by technical variations,
  such as differences in cell size, staining efficiency, and
  instrument-specific artifacts, which can obscure biological signals
  and complicate downstream analysis. This package addresses this
  challenge by implementing a robust framework of linear models
  designed to identify and remove these sources of unwanted
  covariance. By systematically modeling and correcting for technical
  noise, the package enhances the quality and interpretability of
  mass cytometry data, enabling researchers to focus on biologically
  relevant signals.

- [scHiCcompare](/packages/scHiCcompare) This package provides
  functions for differential chromatin interaction analysis between
  two single-cell Hi-C data groups. It includes tools for imputation,
  normalization, and differential analysis of chromatin interactions.
  The package implements pooling techniques for imputation and offers
  methods to normalize and test for differential interactions across
  single-cell Hi-C datasets.

- [scQTLtools](/packages/scQTLtools) This package specializes in
  analyzing and visualizing eQTL at the single-cell level. It can
  read gene expression matrices or Seurat data, or
  SingleCellExperiment object along with genotype data. It offers a
  function for cis-eQTL analysis to detect eQTL within a given range,
  and another function to fit models with three methods. Using this
  package, users can also generate single-cell level visualization
  result.

- [SEraster](/packages/SEraster) SEraster is a rasterization
  preprocessing framework that aggregates cellular information into
  spatial pixels to reduce resource requirements for spatial omics
  data analysis. SEraster reduces the number of spatial points in
  spatial omics datasets for downstream analysis through a process of
  rasterization where single cells’ gene expression or cell-type
  labels are aggregated into equally sized pixels based on a
  user-defined resolution. SEraster is built on an R/Bioconductor S4
  class called SpatialExperiment. SEraster can be incorporated with
  other packages to conduct downstream analyses for spatial omics
  datasets, such as detecting spatially variable genes.

- [shinyDSP](/packages/shinyDSP) This package is a Shiny app for
  interactively analyzing and visualizing Nanostring GeoMX Whole
  Transcriptome Atlas data. Users have the option of exploring a
  sample data to explore this app's functionality. Regions of
  interest (ROIs) can be filtered based on any user-provided
  metadata. Upon taking two or more groups of interest, all pairwise
  and ANOVA-like testing are automatically performed. Available
  ouputs include PCA, Volcano plots, tables and heatmaps. Aesthetics
  of each output are highly customizable.

- [Site2Target](/packages/Site2Target) Statistics implemented for
  both peak-wise and gene-wise associations. In peak-wise
  associations, the p-value of the target genes of a given set of
  peaks are calculated. Negative binomial or Poisson distributions
  can be used for modeling the unweighted peaks targets and
  log-nromal can be used to model the weighted peaks. In gene-wise
  associations a table consisting of a set of genes, mapped to
  specific peaks, is generated using the given rules.

- [smoppix](/packages/smoppix) Test for univariate and bivariate
  spatial patterns in spatial omics data with single-molecule
  resolution. The tests implemented allow for analysis of nested
  designs and are automatically calibrated to different biological
  specimens. Tests for aggregation, colocalization, gradients and
  vicinity to cell edge or centroid are provided.

- [sosta](/packages/sosta) sosta (Spatial Omics STructure Analysis)
  is a package for analyzing spatial omics data to explore tissue
  organization at the anatomical structure level. It reconstructs
  anatomically relevant structures based on molecular features or
  cell types. It further calculates a range of metrics at the
  structure level to quantitatively describe tissue architecture. The
  package is designed to integrate with other packages for the
  analysis of spatial omics data.

- [spacexr](/packages/spacexr) Spatial-eXpression-R (spacexr) is a
  package for analyzing cell types in spatial transcriptomics data.
  This implementation is a fork of the spacexr GitHub repo
  (https://github.com/dmcable/spacexr), adapted to work with
  Bioconductor objects. The original package implements two
  statistical methods: RCTD for learning cell types and CSIDE for
  inferring cell type-specific differential expression. Currently,
  this fork only implements RCTD, which learns cell type profiles
  from annotated RNA sequencing (RNA-seq) reference data and uses
  these profiles to identify cell types in spatial transcriptomic
  pixels while accounting for platform-specific effects. Future
  releases will include an implementation of CSIDE.

- [SpatialExperimentIO](/packages/SpatialExperimentIO) Read in
  imaging-based spatial transcriptomics technology data. Current
  available modules are for Xenium by 10X Genomics, CosMx by
  Nanostring, MERSCOPE by Vizgen, or STARmapPLUS from Broad
  Institute. You can choose to read the data in as a
  SpatialExperiment or a SingleCellExperiment object.

- [spatialFDA](/packages/spatialFDA) spatialFDA is a package to
  calculate spatial statistics metrics. The package takes a
  SpatialExperiment object and calculates spatial statistics metrics
  using the package spatstat. Then it compares the resulting
  functions across samples/conditions using functional additive
  models as implemented in the package refund. Furthermore, it
  provides exploratory visualisations using functional principal
  component analysis, as well implemented in refund.

- [SplineDV](/packages/SplineDV) A spline based scRNA-seq method for
  identifying differentially variable (DV) genes across two
  experimental conditions. Spline-DV constructs a 3D spline from 3
  key gene statistics: mean expression, coefficient of variance, and
  dropout rate. This is done for both conditions. The 3D spline
  provides the “expected” behavior of genes in each condition. The
  distance of the observed mean, CV and dropout rate of each gene
  from the expected 3D spline is used to measure variability. As the
  final step, the spline-DV method compares the variabilities of each
  condition to identify differentially variable (DV) genes.

- [SVP](/packages/SVP) SVP uses the distance between cells and cells,
  features and features, cells and features in the space of MCA to
  build nearest neighbor graph, then uses random walk with restart
  algorithm to calculate the activity score of gene sets (such as
  cell marker genes, kegg pathway, go ontology, gene modules,
  transcription factor or miRNA target sets, reactome pathway, ...),
  which is then further weighted using the hypergeometric test
  results from the original expression matrix. To detect the
  spatially or single cell variable gene sets or (other features) and
  the spatial colocalization between the features accurately, SVP
  provides some global and local spatial autocorrelation method to
  identify the spatial variable features. SVP is developed based on
  SingleCellExperiment class, which can be interoperable with the
  existing computing ecosystem.

- [TaxSEA](/packages/TaxSEA) TaxSEA is an R package for Taxon Set
  Enrichment Analysis, which utilises a Kolmogorov-Smirnov test
  analyses to investigate differential abundance analysis output for
  whether there are alternations in a-priori defined sets of taxa
  from five previously published databases (BugSigDB, MiMeDB,
  GutMGene, mBodyMap and GMRepoV2). TaxSEA takes as input a list of
  taxonomic identifiers (e.g. species names, NCBI IDs etc.) and a
  rank (E.g. fold change, correlation coefficient). TaxSEA be applied
  to any microbiota taxonomic profiling technology (array-based, 16S
  rRNA gene sequencing, shotgun metagenomics & metatranscriptomics
  etc.) and enables researchers to rapidly contextualize their
  findings within the broader literature to accelerate interpretation
  of results.

- [TENET](/packages/TENET) TENET identifies key transcription factors
  (TFs) and regulatory elements (REs) linked to a specific cell type
  by finding significantly correlated differences in gene expression
  and RE methylation between case and control input datasets, and
  identifying the top genes by number of significant RE DNA
  methylation site links. It also includes many additional tools to
  aid in visualization and analysis of the results, including plots
  displaying and comparing methylation and expression data and RE DNA
  methylation site link counts, survival analysis, TF motif searching
  in the vicinity of linked RE DNA methylation sites, custom TAD and
  peak overlap analysis, and UCSC Genome Browser track file
  generation. A utility function is also provided to download
  methylation, expression, and patient survival data from The Cancer
  Genome Atlas (TCGA) for use in TENET or other analyses.

- [terapadog](/packages/terapadog) This package performs a Gene Set
  Analysis with the approach adopted by PADOG on the genes that are
  reported as translationally regulated (ie. exhibit a significant
  change in TE) by the DeltaTE package. It can be used on its own to
  see the impact of translation regulation on gene sets, but it is
  also integrated as an additional analysis method within
  ReactomeGSA, where results are further contextualised in terms of
  pathways and directionality of the change.

- [TrIdent](/packages/TrIdent) The `TrIdent` R package automates the
  analysis of transductomics data by detecting, classifying, and
  characterizing read coverage patterns associated with potential
  transduction events. Transductomics is a DNA sequencing-based
  method for the detection and characterization of transduction
  events in pure cultures and complex communities. Transductomics
  relies on mapping sequencing reads from a viral-like particle
  (VLP)-fraction of a sample to contigs assembled from the metagenome
  (whole-community) of the same sample. Reads from bacterial DNA
  carried by VLPs will map back to the bacterial contigs of origin
  creating read coverage patterns indicative of ongoing transduction.

- [visiumStitched](/packages/visiumStitched) This package provides
  helper functions for working with multiple Visium capture areas
  that overlap each other. This package was developed along with the
  companion example use case data available from
  https://github.com/LieberInstitute/visiumStitched_brain.
  visiumStitched prepares SpaceRanger (10x Genomics) output files so
  you can stitch the images from groups of capture areas together
  with Fiji. Then visiumStitched builds a SpatialExperiment object
  with the stitched data and makes an artificial hexogonal grid
  enabling the seamless use of spatial clustering methods that rely
  on such grid to identify neighboring spots, such as PRECAST and
  BayesSpace. The SpatialExperiment objects created by visiumStitched
  are compatible with spatialLIBD, which can be used to build
  interactive websites for stitched SpatialExperiment objects.
  visiumStitched also enables casting SpatialExperiment objects as
  Seurat objects.

- [vmrseq](/packages/vmrseq) High-throughput single-cell measurements
  of DNA methylation allows studying inter-cellular epigenetic
  heterogeneity, but this task faces the challenges of sparsity and
  noise. We present vmrseq, a statistical method that overcomes these
  challenges and identifies variably methylated regions accurately
  and robustly.

- [XAItest](/packages/XAItest) XAItest is an R Package that
  identifies features using eXplainable AI (XAI) methods such as SHAP
  or LIME. This package allows users to compare these methods with
  traditional statistical tests like t-tests, empirical Bayes, and
  Fisher's test. Additionally, it includes a system that enables the
  comparison of feature importance with p-values by incorporating
  calibrated simulated data.

- [xCell2](/packages/xCell2) xCell2 provides methods for cell type
  enrichment analysis using cell type signatures. It includes three
  main functions - 1. xCell2Train for training custom references
  objects from bulk or single-cell RNA-seq datasets. 2.
  xCell2Analysis for conducting the cell type enrichment analysis
  using the custom reference. 3. xCell2GetLineage for identifying
  dependencies between different cell types using ontology.

- [XeniumIO](/packages/XeniumIO) The package allows users to readily
  import spatial data obtained from the 10X Xenium Analyzer pipeline.
  Supported formats include 'parquet', 'h5', and 'mtx' files. The
  package mainly represents data as SpatialExperiment objects.

New Data Experiment Packages
=====================

There are 3 new data experiment packages in this release of Bioconductor.

- [humanHippocampus2024](/packages/humanHippocampus2024) This is an
  ExperimentHub Data package that helps to access the
  spatially-resolved transcriptomics and single-nucleus RNA
  sequencing data. The datasets are generated from adjacent tissue
  sections of the anterior human hippocampus across ten adult
  neurotypical donors. The datasets are based on
  [spatial_hpc](https://github.com/LieberInstitute/spatial_hpc)
  project by Lieber Institute for Brain Development (LIBD)
  researchers and collaborators.

- [muSpaData](/packages/muSpaData) Data package containing a
  multi-sample multi-group spatial dataset in SpatialExperiment
  Bioconductor object format.

- [TENET.ExperimentHub](/packages/TENET.ExperimentHub) ExperimentHub
  package containing datasets for use in the TENET package's vignette
  and function examples. These include a variety of different objects
  to illustrate different datasets used in TENET functions. Where
  applicable, all datasets are aligned to the hg38 human genome.

New Annotation Packages
=====================

There are no new annotation packages in this release of Bioconductor.

New Workflow Packages
=====================

There are no new workflow packages in this release of Bioconductor.

New Online Books
=====================

There are no new books in this release of Bioconductor.

NEWS from existing Software Packages
===================================


[ADaCGH2](/packages/ADaCGH2)
-------

                 Changes in version 2.47.2 (2025-04-08)                 

- Updated NEWS: Removed dependency on GLAD:
  GLAD was failing in 2025-02 on all platforms on BioC 3.21;
  it is now an Enhances.

                 Changes in version 2.47.1 (2025-02-19)                 

- Removed dependency on GLAD. GLAD fails in all platforms on BioC 3.21.

[affxparser](/packages/affxparser)
----------

                 Changes in version 1.79.1 (2025-01-06)                 

Bug Fixes

- The package failed to compile in R-devel (to become R 4.5.0). This
i
because the C-level API of R is being tidied up and now requiring
the R_ and Rf_ prefixes to be used. For example, R_Calloc() should
be used instead of Calloc(), and Rf_allocVector() should be used
instead of allocVector().

                 Changes in version 1.79.0 (2024-10-30)                 

Notes

- The version number was bumped for the Bioconductor develop version,
which is now Bioconductor 3.21 for R (>= 4.5).

[alabaster.base](/packages/alabaster.base)
--------------

                        Changes in version 1.8.0                        

- Added createDedupSession(), addObjectToDedupSession() and
  checkObjectInDedupSession() utilities, that allow saveObject()
  method developers to check for duplicated objects and avoid
  re-saving them. Most of this functionality was migrated and
  generalized from the alabaster.matrix package.

- Added the cloneDirectory() function to easily clone the
  contents of an existing directory, either by copying or linking
  to the files. This should be used for creating lightweight
  on-disk representations of duplicated objects.

- Added the absolutizePath() function to obtain an absolute file
  path. This is occasionally necessary to ensure that references
  to deduplicated resources are robust to changes to the working
  directory.

- Coerce numeric_version instances to strings or character
  vectors before saving them in saveObject(). This includ such
  objects encountered inside lists or DataFrames.

- Added getSaveEnvironment() and related functions to keep track
  of the R environment used to save each object. This facilitates
  debugging and enables corrective actions to be taken by loading
  functions when encountering files created by buggy versions of
  any package in the alabaster stack.

- Added preliminary support for custom variable length string
  (VLS) arrays, which compress the heap for more efficient
  storage than HDF5's VLS implementation. This is enabled in the
  saveObject() methods for standalone character vectors as well
  as character vectors in DataFrames and lists via the relevant
  character.vls= option. Setting character.vls=NULL will
  automatically switch between fixed-length string datasets and
  VLS arrays based on the estimated storage. (Defaults to FALSE
  for back-compatibility.)

[alabaster.matrix](/packages/alabaster.matrix)
----------------

                        Changes in version 1.8.0                        

- Deprecated the WrapperArraySeed in favor of direct inheritance
  from DelayedUnaryIsoOp.

- Support deduplication of arrays and sparse matrices via the new
  array.dedup.session= argument in their corresponding
  saveObject() methods. This also allows deduplication of
  external DelayedArray seeds if they have already been processed
  outside a DelayedArray in the same session.

- Added a DelayedArray.force.external= argument to the
  saveObject() method for DelayedArrays. This is a more
  convenient way of forcing seeds to be saved as external arrays
  for deduplication.

- Support the new custom variable length string format for
  storing character arrays in HDF5. This is controlled by the
  array.character.vls= option in the array's saveObject() method.


[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.20.0                        

USER VISIBLE CHANGES

- (v 1.19.5) Added host slot to Service class slot to show any
subdomains in the API host URL.

- (v 1.19.1) Deprecated av*, gcloud, etc. functions are now defunct;
see *-defunct documentation pages.

BUG FIXES AND MINOR IMPROVEMENTS

- (v 1.19.9) Trigger updates based on changes to Dockstore API

- (v 1.19.7) Updated NEWS.md and GitHub Actions to automate version
updates in Dockstore.

- (v 1.19.6) Remove test meant for rapiclient client package.

- (v 1.19.4) Update Dockstore API version

- (v 1.19.3) Remove examples and tests for defunct functions.

- (v 1.19.2) Use gcloud_exists from the AnVILGCP package.

[AnVILGCP](/packages/AnVILGCP)
--------

                        Changes in version 1.2.0                        

New features

- Remove defunct drs_stat, drs_access_url, and drs_cp functions
- Defunct bucket_location and bucket arguments

Bug fixes and minor improvements

- Filter out "prerequisites" element from response (#103)
- Clean up unused authenticate_* code
- Remove httr2 from Suggests

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.31.1                        

- Update email address.

[ATACseqTFEA](/packages/ATACseqTFEA)
-----------

                        Changes in version 1.9.1                        

- update email address.


[bandle](/packages/bandle)
------

                        Changes in version 1.11                         

Changes in version bandle 1.11.1

- Use C++14

Changes in version bandle 1.11.0

- push to Bioc devel


[basilisk.utils](/packages/basilisk.utils)
--------------

                       Changes in version 1.20.0                        

- Detect the BASILISK_EXTERNAL_FALLBACK_R environment to use an
  external fallback R environment.

- Updated to the latest version of the Miniforge installer
  (24.11.3-0).

- Updated the version of reticulate in the fallback R environment
  to 1.40.0.

- Expose the default version of the Miniforge installation via
  the new defaultMiniforgeVersion constant.

[BatchQC](/packages/BatchQC)
-------

                        Changes in version 2.2.5                        

Minor Changes

- Added sva batch correction method (for unknown variation
correction)

                        Changes in version 2.2.4                        

Major Changes

- Updated variation ratio analysis to be log transformed
- Added umap exploratory option

                        Changes in version 2.2.3                        

Major Changes

- Added limma as a batch correction
- Added umap plot option
- Added ellipses to pca plot

Minor Changes

- Set p-val plot x scale to always be 0 to 1
- Removed Intercept from p-val violin plots

                        Changes in version 2.2.2                        

Minor Changes

- Updated SE object upload to allow assays of any name (no longer
require one assay to be called "counts")

                        Changes in version 2.2.1                        

Minor Changes

- Added pval information to the DESeq2 binomial evaluation

[beachmat](/packages/beachmat)
--------

                       Changes in version 2.24.0                        

- Automatically cast integer/logical NAs to their double
  equivalents in the matrices wrapped by initializeCpp(). This
  can be disabled by setting .check.na=FALSE.

- Added a initializeCpp() method for the ConstantArraySeed class.

- Using a TileDBArraySeed in initializeCpp() now automatically
  dispatches to beachmat.tiledb.

[benchdamic](/packages/benchdamic)
----------

                 Changes in version 1.13.1 (2024-11-20)                 

- Bug-fix: changed an internal parameter of DA_ANCOM() to make it works
  after the new release

                 Changes in version 1.13.0 (2024-10-29)                 

- Bump x.y.z version to odd y following creation of RELEASE_3_20 branch

                 Changes in version 1.12.1 (2024-11-20)                 

- Porting the changes of devel version 1.13.1 to release

[BiocCheck](/packages/BiocCheck)
---------

                       Changes in version 1.44.0                        

NEW FEATURES

- Add `lifeCycle`, `deprecate_warn`, and `deprecate_stop` to the list
  of
  functions for the `.Deprecated` / `.Defunct` check.

- Support checking `\()` syntax as anonymous R functions.

- Exclude comment-only lines in function lengths check (@hpages, #220)

- Add `checkFndPerson()` to validate Authors@R / Authors fields in the
  `DESCRIPTION` file. A corresponding `NOTE` is shown in the vignette
  (#215)

- Refactor function length check to count only coding lines and not
  comment
  lines (@hpages, #220)

BUG FIXES AND MINOR IMPROVEMENTS

- Check file sizes only for source directories and in
  `BiocCheckGitClone`
  (@hpages, #219)

- Use R version from `BiocManager:::.version_map` for
  `checkRVersionDependency` (@helenalc, #223)

- Use `checkInstalled()` to verify packages listed in `Suggests:`

- Suppress warning in `.hasPkg()` when the package is not installed

- Use `INFO` messages instead of bullet points in output

[BiocFileCache](/packages/BiocFileCache)
-------------

                        Changes in version 2.15                         

BUG FIX

- (2.15.1) Fix issue trying to create multiple caches in single R
  session

[BiocGenerics](/packages/BiocGenerics)
------------

                       Changes in version 0.54.0                        

NEW FEATURES

- Add longForm() S4 generic.
  See https://github.com/waldronlab/MultiAssayExperiment/issues/333 for
  some context.

- Add paste2() S4 generic, with methods defined for ordinary vectors
  and arrays. Also add add_prefix() and add_suffix(), two simple
  wrappers
  around paste2() provided for convenience and code readability.

- Define setequal() S4 generic with generics::setequal as default
  method.

SIGNIFICANT USER-VISIBLE CHANGES

- Add CRAN package generics to Depends field. The default methods for
  S4 generics union(), intersect(), and setdiff() now are
  generics::union(), generics::intersect(), and generics::setdiff(),
  respectively. See '?BiocGenerics::setops' for more information.

[BiocHubsShiny](/packages/BiocHubsShiny)
-------------

                        Changes in version 1.8.0                        

Significant user-visible changes

- Save the table of result to a text file instead of Rds.
- Use a helper function to link images to the vignette from either
online or local sources
- Add a README.md file

Bug fixes and minor improvements

- Added a clipboard icon to the shinyAce editor to allow copying the
code to the clipboard.
- Include package anchors in documentation links
- Import write.table from the utils package

[BiocParallel](/packages/BiocParallel)
------------

                        Changes in version 1.42                         

BUG FIXES

- (1.41.3) Patch a compiler error under Apple clang 17.0.0; see
  <https://github.com/Bioconductor/BiocParallel/issues/271>

- (1.41.2) Set options (via `do.call()`) in package namespace,
  rather than base environment / search path. Thanks Charlotte
  Soneson!
  <https://github.com/Bioconductor/BiocParallel/issues/266>

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.26.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Defunct `pkgType = "all"` argument in `biocDownloadStats` and
  `biocRevDepEmail`

- Added `R (>= 4.1.0)` requirement and replaced usage of `%>%`
  (magrittr)
  with native pipe `|>` throughout the package

- Removed dependency on the `magrittr` package

BUG FIXES AND MINOR IMPROVEMENTS

- Filter URLs with helper function that checks availability in
  `biocBuildReportDB`

- Always use `.CRAN_WEB_URL` in `CRAN_pkg_rds_url()` helper

- Updated `getPkgYearsInBioc()` and documentation

- Reformatted output from `anacondaDownloadStats()`

[biocthis](/packages/biocthis)
--------

                       Changes in version 1.17.4                        

NEW FEATURES

- use_bioc_github_action() now creates a
.github/workflows/bioc-check.yml file that has the option
'bioc_version' set to "bioc-release" by default. It can also take
the values "bioc-devel" or a specific Bioconductor version number in
the X.YY format such as "3.20". This updated GitHub Actions workflow
aims to be as static as possible across Bioconductor release
versions. Meaning that in the future you won't need to use
use_bioc_github_action() again to update the GHA workflow after
every Bioconductor release. This was motivated by the "dynamic
matrix generation" documentation available at
https://runs-on.com/github-actions/the-matrix-strategy/#dynamic-matrix-generation.
For a detailed explanation of the changes in this version, check the
LIBD RStats club presentation "[2025-03-21] biocthis v1.17.4 updated
GitHub Actions workflow". The video is available at
https://www.youtube.com/watch?v=bzzPBt3Mz0A and the notes at
https://docs.google.com/document/d/1z8xkC_3kAsGlpF_UyM9YfV-08o1lTer7xoBz2lQ1IjY/edit?usp=sharing.

                       Changes in version 1.17.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- use_bioc_github_action() now uses simpler code for caching R
packages across GitHub Action runs. This code won't need to be
updated across versions, meaning that it'll be easier to update the
.github/workflows/check-bioc.yml file across Bioconductor versions.
Also, @gaborcsardi's recent commit to r-lib/actions as noted at
https://github.com/r-lib/actions/issues/912#issuecomment-2667950006
gave me a hint on how to simplify code across macOS/winOS and linux
for the caching of R packages. Also, note that thanks to
https://github.com/r-lib/remotes/commit/0e4e23051041d9f1b15a5ab796defec31af6190d
we will soon be able to re-enable automatic installation of linux
system dependencies thanks to remotes::system_requirements("ubuntu",
"24.04") being supported in the near future. Finally, we no longer
need to specify RSPM although there are detailed instructions at
https://packagemanager.posit.co/client/#/repos/bioconductor/setup on
how to do so.

                       Changes in version 1.17.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- use_bioc_description() now has report_bioc set to FALSE by default
to match the Bioconductor guidelines listed at
https://contributions.bioconductor.org/description.html#description-bugreport.
This was brought up in a December 2024 package review at
https://github.com/Bioconductor/Contributions/issues/3503#issuecomment-2551233199.
- Similarly, use_bioc_pkg_templates() now adds the dev/ directory to
the main .gitignore file. This was also brought up in the same
package review.
- use_bioc_vignette() no longer creates a vignettes/.gitignore file.
The rendered .R and .html files are instead ignored on the package
main .gitignore file. This complies with the request from that same
package review while also keeping in line with the behavior from
usethis::use_vignette()
https://github.com/r-lib/usethis/blob/a653d6e05f9172772cea1055f8415cda2f26de69/R/vignette.R#L11-L12.
- use_bioc_vignette()'s template no longer tracks how much time was
used to render the vignette, nor shows the code used for obtaining
the .R file with knit(tangle = TRUE), nor the date the vignette was
generated. This was also brought up in the same package review from
December 2024, as well as in a second one from October 2024 at
https://github.com/Bioconductor/Contributions/issues/3501#issuecomment-2408081535.

                       Changes in version 1.17.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- use_bioc_github_action() now uses the actions/cache@v4 instead of
v3.

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.64.0                        

BUG FIXES

- Version 1.10 of {httr2} changed how URLs are parsed, and this broke
  some biomaRt functionality.  This has been patched.
  (Backported to biomaRt 2.62.1)

[BioNAR](/packages/BioNAR)
------

                        Changes in version 1.9.2                        

- Basic clustering is component wise now. That means that before
clustering algorithm application graph split into weak components,
to which clustering algorithms applied independently. If algorithm
fails on some component, whole component considered as a cluster.

[Biostrings](/packages/Biostrings)
----------

                       Changes in version 2.76.0                        

NEW FEATURES

- Add functions update_DNA_palette(), update_RNA_palette(),
  update_AA_palette(), and update_B_palette(), to let users
  set their own color palettes for XString/XStringSet objects.
  Author: Aidan Lakshman

SIGNIFICANT USER-VISIBLE CHANGES

- Bring back the BLOSUM* and PAM* scoring matrices from pwalign but now
  they display a warning that their new home is the pwalign package and
  that they will soon be removed from the Biostrings package.

DEPRECATED AND DEFUNCT

- Deprecate XStringPartialMatches objects. They've been a stale
  work-in-progress since 2008 so time to get rid of them!

- Formally deprecate all the functions that went to pwalign in BioC
  3.19.

[BreastSubtypeR](/packages/BreastSubtypeR)
--------------

                        Changes in version 1.0.0                        

Major changes

- Initial release of the BreastSubtypeR package.
- Implements NC-based and SSP-based molecular subtyping methods for
breast cancer.
- Provides an R Shiny app for users unfamiliar with R.

Features

- Comprehensive Intrinsic Subtyping for Breast Cancer: Integrates
multiple published intrinsic subtyping methods, including NC-based
approaches like the original PAM50 (Parker et al., J Clin
Oncol, 2009) and SSP-based methods like AIMS (Paquet et al., J Natl
Cancer Inst, 2015).
- Multi-Method Subtyping Functionality: Simultaneously predicts
breast
cancer intrinsic subtypes using a variety of validated methods for
comparative analysis.
- AUTO Mode Feature: Evaluates the distribution of ER and HER2 status
in the test cohort to automatically select subtyping methods that
align with the cohort's characteristics, ensuring compatibility with
method-specific assumptions for greater accuracy and reliability.
- Optimized Gene Mapping: Optimizes gene mapping using Entrez IDs to
maximize the inclusion of genes across subtyping methods.
- Streamlined Input and Output: Provides standardized input/output
formats to ensure smooth integration with other gene expression
analysis tools.
- User-Friendly Shiny App Interface: Features a web-based GUI that
runs entirely locally, ensuring data privacy with no online sharing,
ideal for users who prefer a visual interface over R scripting.

[BridgeDbR](/packages/BridgeDbR)
---------

                       Changes in version 2.17.1                        

- Updated to BridgeDb 3.0.28
- Added CiTO annotation to the vignette (not visible in the HTML)
- Added session info to the vignette

[BSgenome](/packages/BSgenome)
--------

                       Changes in version 1.76.0                        

DEPRECATED AND DEFUNCT

- Formally deprecate BSgenome::forgeMasksFile, BSgenome::forgeSeqFiles,
  BSgenome::forgeSeqlengthsRdaFile, BSgenome::forgeSeqlengthsRdsFile,
  BSgenome::forgeBSgenomeDataPkg, and
  BSgenome::forgeMaskedBSgenomeDataPkg.
  All these functions are now defined in the BSgenomeForge package.

[BSgenomeForge](/packages/BSgenomeForge)
-------------

                        Changes in version 1.8.0                        

- No significant changes in this version.

[bugsigdbr](/packages/bugsigdbr)
---------

                       Changes in version 1.13.6                        

- Add domain support for signatures

                       Changes in version 1.13.2                        

- Update functions to accommodate PMID as Study ID


[Cardinal](/packages/Cardinal)
--------

                 Changes in version 3.8.3 (2024-12-18)                  

BUG FIXES

- Fix error in 'selectROI()' for whole-experiment images

                 Changes in version 3.8.2 (2024-12-11)                  

SIGNIFICANT USER-VISIBLE CHANGES

- Update 'peakAlign()' with 'binfun' parameter

- Previous (3.6) behavior corresponds to 'binfun="median"'

BUG FIXES

- Fix error in 'peakAlign()' for high resolution peaks

- Set minimum relative alignment resolution of 0.5 ppm

- Set minimum absolute alignment resolution of 0.0001 mz

- Fix 'peakProcess()' not processing centroided 'MSImagingArrays'

                 Changes in version 3.8.1 (2024-12-10)                  

BUG FIXES

- Fix error in 'simulateSpectra()' when 'mz' is unsorted

- Add 'resolution' and 'fmax' to 'simulateImage()'

[CatsCradle](/packages/CatsCradle)
----------

                        Changes in version 1.1.2                        

- Function getSubsetComponents added. This is designed to dectect the
components of a gene subset in the case where median complement
distance detects clustering.
- Function edgeLengthsAndCellTypePairs rewritten for faster runtime.

                        Changes in version 1.1.1                        

- Default behaviour of nbhdsAsEdgesToNbhdsAsList changed to exclude
self. This fixes a bug leading to inflated values of MoransI.

[cBioPortalData](/packages/cBioPortalData)
--------------

                       Changes in version 2.20.0                        

New features

- Add copyNumberData helper function to cBioPortalData
- Add studyId argument to the fetchData function
- Use molecularProfiles to determine mutation or copy number
molecular
profile type; previously this was deduced from the data name
- Add example for filtering discrete copy number data

Bug fixes and minor improvements

- Caching calls removed from cBioPortalData. API requests will always
be performed. cBioDataPack will continue to cache data files.
- Filter out NULL parameters in internal .invoke_fun() API calls
- Fix edge case in clinical data parsing when only a single file is
present
- Update long unit tests expectations to reflect current study
coverage (80%)
- Various documentation and example improvements

[celda](/packages/celda)
-----

                 Changes in version 1.22.1 (2024-10-29)                 

- Fixed issue with enrichR not being loaded

[cfTools](/packages/cfTools)
-------

                        Changes in version 1.8.0                        

NEW FEATURES

- Add the `PlotFractionPie()` function, which generates a pie chart for
  a
  vector of cfDNA fractions.

MINOR IMPROVEMENTS

- Updated output column names for clarity and adjusted file format for
  consistency across functions.

- Suppressed non-critical TensorFlow warning messages.

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.41.1                        

- Update Jianhong's email address.

[ClassifyR](/packages/ClassifyR)
---------

                       Changes in version 3.12.0                        

- 
  Metric calculation restored to permutation level by default,
  modifiable by parameter.

- 
  crissCrossValidate has TOP mode.

- 
  easyHard convenience function for calculating sample-specific
  accuracy and regressing covariates to it.

- 
  Parameter tuning grid no longer executed by default.

- 
  crossValidate's formal parameters reduced.

[cleanUpdTSeq](/packages/cleanUpdTSeq)
------------

                       Changes in version 1.45.1                        

- Update email address.

[clusterProfiler](/packages/clusterProfiler)
---------------

                       Changes in version 4.15.2                        

- more general regular pattern to remove species suffix in KEGG
pathway name (2025-02-27, Thu)
- remove input duplicated genes in groupGO() (2024-11-29, Fri, #741)

                       Changes in version 4.15.1                        

- simplify() keeps the most informative term if there exist multiple
terms that meets the criteria (2024-11-29, Fri, #744)
- add 'RichFactor', 'FoldEnrichment' and 'zScore' in enrichDAVID()
result (2024-11-12, Tue)
- update DAVID Web Service URL to make enrichDAVID() work properly
(2024-11-09, Sat)
- https://davidbioinformatics.nih.gov/content.jsp?file=WS.html
- add new citation (2024-11-07, Thu)
- https://doi.org/10.1016/j.xinn.2024.100722

[ClustIRR](/packages/ClustIRR)
--------

                       Changes in version 1.5.41                        

- Function for community inspection/decoding

                       Changes in version 1.5.36                        

- Two function for visualization of results: get_honeycombs and
  get_beta_violins

[cola](/packages/cola)
----

                       Changes in version 2.13.1                        

- register functions are put inside `.onAttach()`

[CompoundDb](/packages/CompoundDb)
----------

                        Changes in version 1.11                         

Changes in version 1.11.2

- Import extractByIndex() from ProtGenerics.

Changes in version 1.11.1

- Complete unit test coverage.

[COTAN](/packages/COTAN)
-----

                        Changes in version 2.7.5                        

Solved issue with use of parallelDist::parDist() by allowing a
fall-back
to stats::dist(). This was needed to address failures running
parDist()
on Linux-aarch64 machines

                        Changes in version 2.7.4                        

Solved bug causing errors while using torch with a CPU device

Ensure the drop out cluster from cellsUnifromClustersing() [-1] keeps
its name if it has not been merged at the end of the function
mergeUniformCellsClusters()

Stopped using broken BioConductor PCAtools::pca: using
BioSingular::runPCA() instead

Added new utility function asClusterization() that takes any input
representing a clusterization (factor, vector or data.frame) and
makes
it into a factor. This function is now used in all functions taking
in a
clusterization to standardize the given input

Added initialIteration as input parameter to clusterization functions
so
to avoid override of partial data when the functions are being
restarted

                        Changes in version 2.7.3                        

Added new function clusterGeneContingencyTables(): it returns
observed
and expected contingency tables for a given gene and group of cells
(cluster)

                        Changes in version 2.7.2                        

Added possibility to specify genes' selection method used in the
cellsUniformClustering() function

Improved "AdvancedGDIUniformityCheck" by adding a new check about 99%
quantile

Solved issue in clustersMarkersHeatmapPlot(): now passing in a genes
not
in the data-set will result in a corresponding gray column instead of
an
error

Added possibility to specify clean() thresholds in the functions
proceedToCoex() and automaticCOTANObjectCreation()

                        Changes in version 2.7.1                        

Improved function clustersMarkersHeatmapPlot(): now its shows a
column
for each marker gene and the shown content is more expressive

Marked the function clustersDeltaExpression() as legacy: it has been
replaced with the function DEAOnClusters() in the package

Fixed minor bug in class AdvancedGDIUniformityCheck regarding third
check: was testing third highest GDI value instead of second

[CRISPRball](/packages/CRISPRball)
----------

                         Changes in version 1.4                         

- Bug Fixes:
1.  plot_rank(): Fixed silent bug with h.id suffix being the same as
plot_volcano().

[CytoMDS](/packages/CytoMDS)
-------

                         Changes in version 1.3                         

CytoMDS 1.3.7

- changed package citation (and in README.md) to peer reviewed
journal
article

CytoMDS 1.3.6

- ggplotMarginalDensities(): corrected a bug when channels were
specified as marker names.

CytoMDS 1.3.5

- removed contraint (max = 3) on nb of plotly tool_tips variables
(pDataForAdditionalLabelling)

CytoMDS 1.3.4

- corrected plotly tool_tips (pDataForAdditionalLabelling)

CytoMDS 1.3.3

- DistSum class to store distance matrices computed as the sum of
marker contributions
- ggplotDistFeatureImportance() now can be used to create a stacked
bar ggplot object, displaying feature importance in a distance
matrix (extracted from a DistSum object)

CytoMDS 1.3.2

- added pointSize argument to ggplotSampleMDS()

CytoMDS 1.3.1

- corrected bug with nPoints() method of MDS class, due to slightly
modified behaviour of dist S3 class in R4.5.

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

                       Changes in version 1.45.1                        

- update email address.

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

[DegCre](/packages/DegCre)
------

                        Changes in version 1.3.0                        

NEW FEATURES

- Again updated collapseDegCreToGene to run much faster. No longer
accepts a useParallel argument. Does not use BiocParallel.

BUG FIXES

- Fixed bug in convDegCreResListToCreGeneScoreGR that called a
deprecated function.

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.34.0                        

NEW FEATURES

- Implement delayed paste2() for DelayedArray objects. This makes
  BiocGenerics::add_prefix() and BiocGenerics::add_suffix() work
  out-of-the box on DelayedArray objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Binary operations between a DelayedArray object and an array now
  raise an error. This is temporary and only to prevent users from
  trying to use a tile object (new in S4Arrays 1.7.3) on the left or
  right side of the operation, and get the wrong result.

- Revisit implementation of the path() getter (and setter) for
  DelayedOp.
  The new implementation is recursive so will acknowledge path()
  methods
  defined for custom DelayedOp derivatives like the path() getter
  defined
  for ReloadedArraySeed objects in the alabaster.matrix package.

DEPRECATED AND DEFUNCT

- SparseArraySeed objects, OLD_extract_sparse_array(), and
  read_sparse_block() are now defunct after being deprecated in BioC
  3.20.

- Get rid of previously defunct stuff. This is:
  - blockGrid(): replaced with defaultAutoGrid()
  - rowGrid(): replaced with rowAutoGrid()
  - colGrid(): replaced with colAutoGrid()
  - multGrids(): replaced with defaultMultAutoGrids()
  - viewportApply(): replaced with gridApply()
  - viewportReduce(): replaced with gridReduce()
  - getRealizationBackend(): replaced with getAutoRealizationBackend()
  - setRealizationBackend(): replaced with setAutoRealizationBackend()
  All these functions went defunct 2.5 years ago in BioC 3.15.

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.47.4                        

- Adding `pasilla` extdata to DESeq2 package for build stability.

[DESpace](/packages/DESpace)
-------

                       Changes in version 1.99.1                        

- add detection of differential spatial patterns (DSP) across
  conditions

- rename functions for consistency

- enable FeaturePlot for non-Visium/ST data

[dominoSignal](/packages/dominoSignal)
------------

                        Changes in version 1.0.5                        

- circos_ligand_receptor will not fail in cases where the rl_map of a
domino objects includes ligands for a receptor where the ligands are
not present in the expression matrix. A message is returned if
ligands from the rl_map had to be excluded from the plot

                        Changes in version 1.0.4                        

- refactorization of circos_ligand_receptor function to remove use of
grepl-based regular expressions to reformat cell and molecule names
from a domino object into a data frame for plotting a receptor
circos plot. Component functions for creating the data frame of
ligand expression centered on a single receptor and to render this
data frame as a circos plot are now included in the package as
"obtain_circos_expression" and "render_circos_ligand_receptor",
respectively.

                        Changes in version 1.0.3                        

- Added functions for calculation of mean gene expression among
components of a complex using purrr functions

                        Changes in version 1.0.2                        

Pkgdown Site Customization Scripts

- Restored _pkgdown.yml and index scripts that specify site building
parameters for https://FertigLab.github.io/dominoSignal to improve
formatting.

                        Changes in version 1.0.1                        

Receptor Complex Bugfix

- Resolved issue where if create_domino was run with complexes=TRUE
and no complexes were found to have active signaling, the full
signaling matrix would be replaced by a NULL value.

GitHub Actions

- Restoration of .github/workflows scripts for automatic build
checks.

[dreamlet](/packages/dreamlet)
--------

                        Changes in version 1.4.1                        

- March 25, 2025
- improve checking for sample_id and cluster_id in
aggregateToPseudoBulk()

[DropletUtils](/packages/DropletUtils)
------------

                       Changes in version 1.28.0                        

- Use a more stable algorithm for identifying the knee point in
  barcodeRanks(). The new algorithm is based on maximizing the
  distance from a line between the plateau and the inflection
  point. Previously, we tried to minimize the signed curvature
  but this was susceptible to many local minima due to the
  instability of the empirical second derivative, even after
  smoothing.

- Set alpha=Inf as the default for testEmptyDrops(). This is
  motivated by the realization that an underestimated alpha can
  still yield anticonservative p-values and is not universally
  safer than alpha=Inf. Defaulting alpha=Inf is preferable as it
  is at least correct in the expected case of multinomial
  sampling.

[DuplexDiscovereR](/packages/DuplexDiscovereR)
----------------

                        Changes in version 1.1.2                        

- Date: 2025-03-12
- Removed BiocCheck files

                        Changes in version 1.1.1                        

- Date: 2024-12-30
- Added added raw abundance values and duplex counts for 2x2
contingency table to the p-value calculation output

[EBSeq](/packages/EBSeq)
-----

                        Changes in version 2.5.2                        

- Using c++14 for the updated [BH]
  (https://cran.r-project.org/web/packages/BH/index.html) dependency,
  and restrict BH version <= 1.87.0-1 to prevent higher version C++
  requirement in the future

[edgeR](/packages/edgeR)
-----

                 Changes in version 4.6.0 (2025-04-16)                  

- 
  edgeR now depends on limma >= 3.63.6.

- 
  New diffSplice() method for DGEGLM objects. The output object
  has the same format as that from the diffSplice() method in the
  limma package for MArrayLM objects. The new method is intended
  to replace diffSpliceDGE(), although the latter is kept for the
  current release for backward compatibility. The new
  diffSplice() method is slightly slower than diffSpliceDGE(),
  but has improved power and receiver operating curve while still
  controlling the FDR correctly.

- 
  Corrections to diffSpliceDGE() so that it provides more
  rigorous control of the FDR. The numerator degrees of freedom
  used for the gene-level tests is increased by 1.

- 
  getOffset() now gives an error if the data object doesn't
  contain library sizes.

- 
  The entries in the package DESCRIPTION file have been revised
  and reformated.

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.27.5                        

- able to scale pie size for 'compareClusterResult' (2025-03-11, Tue,
#308, #311)

                       Changes in version 1.27.4                        

- adjust pie size and category label position in cnetplot()
(2025-01-08, Wed, #306)
- clean up code (2024-12-20, Fri)

                       Changes in version 1.27.3                        

- scale pies and add pie legend in emapplot() (2024-12-12, Thu, #304)
- a safe way to extract gene sets in ridgeplot() (2024-12-12, Thu,
#303)

                       Changes in version 1.27.2                        

- emapplot() now allows passing color to a specific color, e.g.,
color
= "black" (2024-11-29, Fri, #300)
- bug fixed in emapplot()
- size_category now works for pie node (2024-11-29, Fri, #301)
- legend of term nodes will be retained when group = TRUE
(2024-11-29, Fri, #300)
- supports passing ID to 'showCategory' in ridgeplot() (2024-11-06,
Wed, #295)
- enhancement of cnetplot() (2024-11-06, Wed)
- 'node_label' can be a vector of selected items/genes to specify
the items to be displayed (#293)
- 'node_label' can be 'exclusive' to label genes that are uniquely
belongs to categories (#253)
- 'node_label' can be 'share' to label genes that are share
between categories (#253)
- 'node_label' can be, e.g. '> 1' or '< 1', to label genes that
have log2FC values larger or smaller than the threshold (#253)
- supports using ggtangle::geom_cnet_label() to label items/genes
in independent layer (#194, #266, #267)
- fixed ridgeplot() when selecting a specific gene set and plotting
non-core genes (2024-11-06, Wed, #298)

                       Changes in version 1.27.1                        

- add 'ID' parameter in goplot() (2024-10-30, Wed)
-
https://github.com/YuLab-SMU/enrichplot/issues/292#issuecomment-2445788948

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.15.3 (2025-03-05)                 

- test coverage for namespaces

- new ref

                 Changes in version 1.15.1 (2024-11-16)                 

- can override alignment endness when loading BAM

[EpiCompare](/packages/EpiCompare)
----------

                       Changes in version 1.11.3                        

- Update minimum version for dependency IRanges. See issue #155.

                       Changes in version 1.11.1                        

Bug Fixes

- prepare_output_build
- Fix support for mouse genome builds (mm9 and mm10).

Documentation

- Epicompare.Rmd
- Mention support for mouse genome builds.

[EpiDISH](/packages/EpiDISH)
-------

                       Changes in version 2.23.1                        

- Added unified lifecourse 19 immune cell-type reference matrix
  for 450k and EPIC arrays.

[escape](/packages/escape)
------

                        Changes in version 2.2.4                        

UNDERLYING CHANGES

- moved dependency from msigdbr to msigdb
- getGeneSets() now locally caches the gene sets to improve speed of
repeated use
- getGeneSets() now only supports Mouse or Human

                        Changes in version 2.2.3                        

UNDERLYING CHANGES

- fixed handling of groups parameter and data splitting in
escape.matrix()
- improved efficiency of internal .split_data.matrix()

                        Changes in version 2.2.2                        

UNDERLYING CHANGES

- fix performNormalization() conditional statements
- fix performNormalization() rescaling for per gene set calculations

                        Changes in version 2.2.1                        

[EWCE](/packages/EWCE)
----

                       Changes in version 1.15.1                        

Bug fixes

- Remove print statements, clean up CTD creation steps.

[extraChIPs](/packages/extraChIPs)
----------

                       Changes in version 1.11.2                        

- Changed handling of arguments in plotting functions, setting NULL
as
the primary default value

[fenr](/packages/fenr)
----

                        Changes in version 1.4.2                        

- Added evidence code column to GO-term mapping table. It can be used
to filter mapping based on their quality. See
https://geneontology.org/docs/guide-go-evidence-codes for
explanation.

                        Changes in version 1.4.1                        

- Attempted to fix a bizarre error message on Bioconductor's test
machines with older version of MacOS. Windows and Linux are not
affected; my laptop running Sequoia 5.2 does not show show errors. I
suspect a memory leak in older systems. The error vector memory
limit of 64.0 Gb reached, see mem.maxVSize() happened in the
function parse_kegg_genes(), a flat-file parser for KEGG. It
occurred around the call tidyr::separate(), which I replaced with an
alternative approach. Will see if the error is fixed.

[fgsea](/packages/fgsea)
-----

                       Changes in version 1.33.3                        

- update code for using `msigdbr` in geseca vignette

[flowcatchR](/packages/flowcatchR)
----------

                       Changes in version 1.42.0                        

Other notes

- Some changes in the source of the documentation, providing anchors
to all function calls. This avoids the new note in R CMD check in
the new major release

[fmrs](/packages/fmrs)
----

                        Changes in version 2.0.1                        

IMPROVEMENTS SINCE LAST RELEASE

- Non-mixture of regression models are now added to the package.

BUG FIXES

- Several bugs are fixed.

                        Changes in version 2.0.0                        

IMPROVEMENTS SINCE LAST RELEASE

- The package is rewritten using .Call function.
- The codes for Weibull distribution are improved.

BUG FIXES

- Several bugs are fixed which caused the results to be different for
the same analysis.

[fobitools](/packages/fobitools)
---------

                       Changes in version 1.15.1                        

- Fix bugs in vignettes

[gDNAx](/packages/gDNAx)
-----

                 Changes in version 1.6.0 (2025-04-14)                  

USER VISIBLE CHANGES

- Handling of metadata from non-model organisms has been improved.
  Thanks to discussions on the support site
  (https://support.bioconductor.org/p/9160897/) and GitHub
  (https://github.com/rcastelo/gDNAx/issues/5).

- Added more documentation to the output of the 'getDx()' method.

[gDRcore](/packages/gDRcore)
-------

               Changes in version 2025-03-31 (2025-03-31)               

- extend the logic of functions for getting annotation from dt to
fill
missing cols with unknown

               Changes in version 2025-03-10 (2025-03-10)               

- move map_conc_to_standardized_conc to gDRutils package

               Changes in version 2025-03-05 (2025-03-05)               

- update assay names for combo data in get_assays_per_pipeline_step
function

               Changes in version 2025-01-13 (2025-01-13)               

- switch to gDRutils::remove_drug_batch function

               Changes in version 2024-12-18 (2024-12-18)               

- fix melt error after changed intersect behaviour

               Changes in version 2024-11-15 (2024-11-15)               

- fix melt in annotation function to fix Bioc error

[gDRimport](/packages/gDRimport)
---------

               Changes in version 2025-04-10 (2025-04-10)               

- replace with_mock with with_mocked_bindings

               Changes in version 2025-03-26 (2025-03-26)               

- unify support for using OncotreeLineage for PRISM level 5 and
level 6

               Changes in version 2025-03-12 (2025-03-12)               

- add support for using OncotreeLineage as a tissue annotation if
available in PRISM data

               Changes in version 2025-02-18 (2025-02-18)               

- update extracting Tissue information for cell lines in PRISM data

               Changes in version 2025-02-07 (2025-02-07)               

- extract Tissue information for cell lines in PRISM data

               Changes in version 2025-02-04 (2025-02-04)               

- update error messages, to make them less confusing

               Changes in version 2025-01-11 (2025-01-11)               

- take into account run in public PRISM data

               Changes in version 2024-12-09 (2024-12-09)               

- support Duration in old private PRISM data format

               Changes in version 2024-12-04 (2024-12-04)               

- support old private PRISM data format

               Changes in version 2024-11-06 (2024-11-06)               

- fix invalid moa for PRISM data

               Changes in version 2024-11-05 (2024-11-05)               

- synchronize Bioconductor and GitHub versioning

[gDRstyle](/packages/gDRstyle)
--------

               Changes in version 2025-02-11 (2025-02-11)               

- make different default lintr config for different versions of lintr

               Changes in version 2025-02-08 (2025-02-08)               

- switch off return_linter in default lintr config

               Changes in version 2025-02-04 (2025-02-04)               

- move BiocStyle to Depends

               Changes in version 2025-01-20 (2025-01-20)               

- update notes.json

               Changes in version 2024-11-05 (2024-11-05)               

- synchronize Bioconductor and GitHub versioning

[gDRutils](/packages/gDRutils)
--------

               Changes in version 2025-03-26 (2025-03-26)               

- fix default parameter in get_settings_from_json

               Changes in version 2025-03-19 (2025-03-19)               

- add get_gDR_session_info function

               Changes in version 2025-03-18 (2025-03-18)               

- move split_big_table_for_xlsx from gDRsearch2 package

               Changes in version 2025-03-12 (2025-03-12)               

- improve logic in get_assay_req_uniq_cols

               Changes in version 2025-03-07 (2025-03-07)               

- add support for combination experiment in cap_assay_infinities
- move map_conc_to_standardized_conc from gDRcore package

               Changes in version 2025-02-21 (2025-02-21)               

- refactor average_biological_replicates

               Changes in version 2025-02-14 (2025-02-14)               

- update default value of capping_fold param in cap_assay_infinities

               Changes in version 2025-02-05 (2025-02-05)               

- add support for dropping masked values in the assay data

               Changes in version 2025-02-03 (2025-02-03)               

- keep 'replicate' column as additional perturbation in
get_additional_variables

               Changes in version 2025-01-30 (2025-01-30)               

- add support for unifying metadata in convert_se_assay_to_dt
function

               Changes in version 2025-01-24 (2025-01-24)               

- add cap_assay_infinities

               Changes in version 2025-01-13 (2025-01-13)               

- refactor remove_drug_batch

               Changes in version 2024-12-10 (2024-12-10)               

- make split_SE_components working correctly for sa assay data,
modified with avearge_biological_duplicates

               Changes in version 2024-12-09 (2024-12-09)               

- minor improvement in the logic of average_biological_replicates
(new
blacklisted column)

               Changes in version 2024-12-02 (2024-12-02)               

- refactor set_unique_* functions

               Changes in version 2024-11-05 (2024-11-05)               

- add get_env_var helper

- synchronize Bioconductor and GitHub versioning

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.44.0                        

UTILITIES

- fix the C++ error with Apple clang version 17: no
  'std::basic_string<unsigned short>' & 'std::basic_string<unsigned
  int>'

                       Changes in version 1.42.2                        

UTILITIES

- tweak the C Macro to support Linux MUSL

- fix according to the C compiler C23

[gemma.R](/packages/gemma.R)
-------

                        Changes in version 4.0.0                        

- get_annotation_children and get_annotation_parents are added.
get_child_terms is deprecated in favor of get_annotation_children

[GeneNetworkBuilder](/packages/GeneNetworkBuilder)
------------------

                       Changes in version 1.49.2                        

- Add warning message if there is no stringApp installed.

                       Changes in version 1.49.1                        

- Update email address.

[GENESIS](/packages/GENESIS)
-------

                       Changes in version 2.37.1                        

- Add option "imputed" to admixMap to allow testing imputed
  dosages.

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 3.2.0                        

Bug fixes

- The tooltip functionality in visualizing the ggs_graph() output now
handles correctly the text for the geneset description. Thanks to
@thomas-keller for spotting this and for the fix!

Other notes

- gene_plot() defaults now to NULL in the intgroup parameter, which
translates into using the first colData item

[GenomAutomorphism](/packages/GenomAutomorphism)
-----------------

                        Changes in version 1.8.1                        

- Add new a new function: 'automorphism_prob', which applies a
Dirichlet-Multinomial Modelling (in a Bayesian framework) to compute
the posterior probability of each type of mutational event.

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.44.0                        

NEW FEATURES

- Register the following NCBI assemblies:
  - GAculeatus_UGA_version5 (Stickleback)
  - GadMor_May2010 and gadMor3.0 (Atlantic cod)
  - Orenil1.1 and O_niloticus_UMD_NMBU (Nile tilapia)
  - CriGri_1.0, C_griseus_v1.0, and CriGri-PICRH-1.0 (Chinese hamster)

- Register the following UCSC genomes:
  - gadMor1 (Atlantic cod), based on NCBI assembly GadMor_May2010
  - oreNil2 (Nile tilapia), based on NCBI assembly Orenil1.1
  - criGri1 (Chinese hamster), based on NCBI assembly C_griseus_v1.0

SIGNIFICANT USER-VISIBLE CHANGES

- Use slightly better heuristic in getChromInfoFromNCBI() for guessing
  the
  circular sequences of non-registered assemblies.

BUG FIXES

- Fix getChromInfoFromNCBI("panpan1.1"). See commit 3be5818.

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.44.0                        

- No changes in this version.

[GenomicDataCommons](/packages/GenomicDataCommons)
------------------

                       Changes in version 1.32.0                        

Bug fixes and minor improvements

- Minor updates to unit tests and GitHub Actions

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.60.0                        

DEPRECATED AND DEFUNCT

- Formally deprecate all the functions that went to txdbmaker in BioC
  3.19.

[GenomicInteractionNodes](/packages/GenomicInteractionNodes)
-----------------------

                       Changes in version 1.11.3                        

BUG FIXES

- Fix multiple typos.

                       Changes in version 1.11.2                        

BUG FIXES

- Fix the error "invalid 'scipen'".

                       Changes in version 1.11.1                        

BUG FIXES

- Update email address.

[GenomicPlot](/packages/GenomicPlot)
-----------

                        Changes in version 1.5.3                        

- Modified handle_bed to allow skipping of first n rows

                        Changes in version 1.5.2                        

- Modified draw_matrix_heatmap to handle momogeneous matrix

                        Changes in version 1.5.1                        

- Modified draw_matrix_heatmap to allow it to handle extreme cases

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.60.0                        

- No changes in this version.

[geomeTriD](/packages/geomeTriD)
---------

                       Changes in version 1.1.11                        

NEW FEATURES

- Add createTADGeometries function.

                       Changes in version 1.1.10                        

NEW FEATURES

- Add pointCluster.

                        Changes in version 1.1.9                        

BUG FIXES

- Fix multiple typos.

                        Changes in version 1.1.8                        

BUG FIXES

- Fix the issues "Uncaught SyntaxError: Expected ',' or ']' after array
  element in JSON" if data.frame is too big.

                        Changes in version 1.1.7                        

BUG FIXES

- Fix the issues of measurment marker do not rotate and unexpected
  resize window when export.

                        Changes in version 1.1.6                        

BUG FIXES

- Export pdf with visible layers only.

                        Changes in version 1.1.5                        

NEW FEATURES

- Add pdf exporter.

                        Changes in version 1.1.4                        

BUG FIXES

- Fix the input check issue in kabsch.

                        Changes in version 1.1.3                        

BUG FIXES

- Update email address.

                        Changes in version 1.1.2                        

NEW FEATURES

- Add resolution for view3dStructure.

                        Changes in version 1.1.1                        

NEW FEATURES

- Add autocomplete input for search key words.

[GEOquery](/packages/GEOquery)
--------

                 Changes in version 2.99.0 (2024-10-01)                 

New Features

- RNAseq data support for GEOquery. Now you can use RNASeq
quantification data prepared by NCBI.
- Basic search in GEO database. Now you can search for datasets in
GEO
database using GEOquery.
- browseGEO() function to open a web browser with a GEO accession.

Breaking changes

- getGEO() now returns a list of SummarizedExperiment objects. This
is
a breaking change from previous versions of GEOquery. If you are
using GEOquery in a script, you will need to update your code to
reflect this change.

Bug Fixes or Improvements

Not an exhaustive list, but some highlights:

- Using httr2 instead of curl for better control over HTTP requests.
- Removed dead gunzip code.

[getDEE2](/packages/getDEE2)
-------

                       Changes in version 1.17.3                        

- Better handling of exact query matches to bundle names.

                       Changes in version 1.17.1                        

- Added support for O sativa (rice) and Z mays (corn) in the DEE2
  database.
  Thanks to Dr Wen-Dar Lin, Bioinformatician at IPMB, Academia Sinica,
  Taiwan, for contributing
  these data.

[ggseqalign](/packages/ggseqalign)
----------

                        Changes in version 1.1.2                        

- Updated installation instructions in Vignette to reflect
Bioconductor release with correct version number

                        Changes in version 1.1.0                        

- Automatic bump from Bioconductor

                        Changes in version 1.0.2                        

- Updated installation instructions in Vignette to reflect
Bioconductor release

                        Changes in version 1.0.1                        

- Updated installation instructions in README to reflect Bioconductor
release

[ginmappeR](/packages/ginmappeR)
---------

                        Changes in version 1.2.3                        

ENHANCEMENT

- .testEquals function revised.

                        Changes in version 1.2.2                        

ENHANCEMENT

- UniProt similar genes tests revised for macOS platforms.

                        Changes in version 1.2.1                        

ENHANCEMENT

- KEGG to NCBI and UniProt similar genes tests revised.

[gINTomics](/packages/gINTomics)
---------

                        Changes in version 1.4.0                        

- New class comparison analysis. Now the package run a single
  analysis for all classes together and then highlights the
  results for differentially expressed genes for all the
  available comparisons.

- New Differentially Methylated genes analysis available.

- New Differential Copy Number Variations analysis available.

- Added DEGs section in the shiny app.

- Added expression insights sections in the shiny app.

- General improvements.

[glmGamPoi](/packages/glmGamPoi)
---------

                  Changes in version 1.19 (2024-11-26)                  

- glmGamPoi is using the new tatami library instead of the old beachmat
  interface.
  These libraries helped to work with different kind of matrices
  (dense, HDF5-backed).
  One important change is that glmGamPoi now requires C++17 (instead of
  C++14).
  (PR#66,thanks Aaron Lun)

- The input data can be sparse! Building on top of tatami allows me to
  support arbitrary
  matrix types such as sparse matrices and DelayedArrays.

- Fix bug regarding matrix columns in `pseudobulk`
  (https://github.com/const-ae/lemur/issues/20).

- Add work-around for
  https://github.com/Bioconductor/DelayedArray/issues/123

[GNOSIS](/packages/GNOSIS)
------

                 Changes in version 1.99.0 (2023-09-04)                 

- Made the following significant changes
  o added functionality to select and upload cBioPortal study
  o deprecated ability to save R script with executed code

- Submitted to Bioconductor

[goSorensen](/packages/goSorensen)
----------

                       Changes in version 1.10.0                        

In this updated version, unless the user specifies otherwise in the
onlyEnriched argument of the enrichedIn function, the enrichment
matrix
of GO terms includes only those terms that are enriched in at least
one
of the lists being compared, excluding all terms that are not
enriched
in any of the lists. This optimization significantly reduces
computation
times for enrichment analysis compared to previous versions, such as
when generating enrichment contingency tables using buildEnrichTable.

Additionally, with the latest versions of Bioconductor and its
updated
packages, we have revised the naming conventions and outputs for
results
obtained from the primary goSorensen functions. These changes are
summarized in two new bullet points, offering clearer guidance and
improved usability for package users compared to earlier versions.

[graphite](/packages/graphite)
--------

                 Changes in version 1.53.1 (2024-04-14)                 

- Updated all pathway data.

[GSVA](/packages/GSVA)
----

                        Changes in version 2.20                         

USER VISIBLE CHANGES

- When the input expression data is provided as either a
  'SummarizedExperiment', 'SingleCellExperiment' or a
  'SpatialExperiment' object, the default assay is now either one
  called 'logcounts' or the first one, otherwise.

- Added a method called 'spatCor()' to efficiently compute spatial
  autocorrelation using Moran's I statistic for an input
  'SpatialExperiment' object, using an inverse squared distance weight
  matrix as default, or an inverse distance weight matrix as an
  alternative. It also tests for spatial autocorrelation assuming
  normality. This functionality is added to help detecting
  autocorrelated gene sets using GSVA enrichment scores.

- After one cycle of deprecation and two cycles of being defunct, the
  classical non-object-oriented way of calling the function 'gsva()'
  has been completely removed. Now 'gsva()' can only be used by
  constructing first a parameter object.

- Fixed URLs in vignette and DOI refs in documentation.

- The vignette has been updated with a new section illustrating how to
  import gene sets from GMT files with the function 'readGMT()' and the
  section on "Differential expression at pathway level" has been
  updated using the bulk RNA-seq data from Costa et. al (2021).

BUG FIXES

- Bugfix in gsvaEnrichment() when input data contains missing values.

- Bugfix in checking and detecting missing values in input expression
  data.

- Bugfix in the PLAGE method when using the argument BPPARAM with a
  value other than Serial

Param(). It resolves issue #220

- Bugfix when tau is not double for the GSVA method.

- Bugfix in the C code to avoid potential memory leaks.

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.36.0                        

NEW FEATURES

- The package now has a vignette that discusses/compares performance of
  various HDF5-based representations in the context of single cell
  analysis.

SIGNIFICANT USER-VISIBLE CHANGES

- h5mread() and other low-level HDF5 manipulation functions have moved
  from the HDF5Array package to the new h5mread package and the former
  now depends on the latter.

DEPRECATED AND DEFUNCT

- All the OLD_extract_sparse_array() and read_sparse_block() methods
  are
  gone (the corresponding generics went defunct in DelayedArray
  0.33.1).

BUG FIXES

- Make sure writeTENxMatrix() can handle a matrix-like object with more
  than 2^31 - 1 nonzero values.

[HoloFoodR](/packages/HoloFoodR)
---------

                        Changes in version 1.1.1                        

Date: 2025-04-03

- Improved MetaboLights support

                        Changes in version 1.1.0                        

Date: 2024-10-24

- Add addMGnify function

[HuBMAPR](/packages/HuBMAPR)
-------

                        Changes in version 1.0.8                        

- Update the vignette plots by adding more disguishable colors

                        Changes in version 1.0.7                        

- Error check

                        Changes in version 1.0.6                        

- Updated NEW.md

                        Changes in version 1.0.5                        

- Updated dataset_metadata() function to accommodate with API changes
- Reformatted vignette with with more technical details based on
comments
- Added summary lists in README

                        Changes in version 1.0.4                        

- Updated dataset_metadata() function to accommodate with API changes
- Reformatted vignette with with more technical details based on
comments
- Added summary lists in README

- Updated functions to accommodate with API changes

[HubPub](/packages/HubPub)
------

                 Changes in version 1.15.3 (2024-01-21)                 

- Update CreateAHubPackageVignette to include R package code for data
  upload

                 Changes in version 1.15.1 (2024-01-13)                 

- Update CreateAHubPackageVignette to use S3 instructions for data
  upload

[ideal](/packages/ideal)
-----

                        Changes in version 2.2.0                        

Other notes

- Some changes in the source of the documentation, providing anchors
to all function calls. This avoids the new note in R CMD check in
the new major release

[immApex](/packages/immApex)
-------

                        Changes in version 1.0.5                        

UNDERLYING CHANGES

- getIMGT() checks for availability of IMGT website
- Expanded Unit Tests
- Unit Tests and Vignette now evaluate for proper python installation
overall

                        Changes in version 1.0.4                        

UNDERLYING CHANGES

- Optional testthat variationalSequences() evaulate presence of Keras

                        Changes in version 1.0.3                        

UNDERLYING CHANGES

- Drop evaluation of variationalSequences() example

                        Changes in version 1.0.2                        

UNDERLYING CHANGES

- Vignette includes eval of keras installation for certain chunks

                        Changes in version 1.0.1                        

UNDERLYING CHANGES

- Fix issue with optimizer call in variationalSequences()
- Move package to keras3
- Version Bump to be consistent with Bioconductor release

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.39.4                        

- -CHANGES
  * sanitize parameter desciption in immunoMeta-object

                      Changes in version 1.39.2-3                       

- CHANGES
  * fixes BioC 3.21 compiler errors

                       Changes in version 1.39.1                        

- CHANGES
  * in subset.immunoClust keep NAs in label of subset-object
  * respect already transformed FCS-data in Default_Scales

[InPAS](/packages/InPAS)
-----

                       Changes in version 2.15.2                        

- Update email address.

                       Changes in version 2.15.1                        

- Add citation.

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.42.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- New implementations of the viewMins(), viewMaxs(), viewSums(),
  viewMeans(), viewWhichMins(), viewWhichMaxs(), viewRangeMins(),
  and viewRangeMaxs() methods for RleViews objects.
  This is a complete rewrite (including of the C code behind these
  methods)
  that addresses many long-standing (and unnoticed-for-years) bugs in
  the
  previous implementation. This new implementation can also be
  significantly
  faster than the previous one e.g. 5x faster or more for RleViews
  objects
  with millions of views.

BUG FIXES

- Fix multiple bugs in viewMins(), viewMaxs(), viewSums(), viewMeans(),
  viewWhichMins(), viewWhichMaxs(), viewRangeMins(), and
  viewRangeMaxs()
  methods for RleViews objects. See above for the details.

- Fix bug in disjointBins() method for IntegerRanges derivatives when
  the input object has names and is not sorted.
  See https://github.com/Bioconductor/IRanges/issues/54 for the
  details.

[iSEE](/packages/iSEE)
----

                       Changes in version 2.19.4                        

- Added a parameter to control the title displayed in the browser
tab,
currently defaulting to "iSEE". Thanks to @RiboRings for the
implementation!

                       Changes in version 2.19.3                        

- Added an option to have no legend shown in the dotplot panel class

                       Changes in version 2.19.1                        

- Changed the processing of the title for the app - if provided
explicitly, now it can also correctly handle HTML content tags and
more. Thanks to @RiboRings for spotting this, resolving #681

[iSEEtree](/packages/iSEEtree)
--------

                         Changes in version 1.1                         

- Updated .multiSelectionResponsive

- Conformed to iSEE v2.19.2

- Created TreePlot family

- Added typical tree operations

- Created GraphPlot family

- Added ScreePlot panel

- Added PrevalencePlot panel

[Lheuristic](/packages/Lheuristic)
----------

                        Changes in version 1.0.0                        

- Package added to Bioconductor

- Bioc-submission branch merged with master

[limma](/packages/limma)
-----

                 Changes in version 3.64.0 (2025-04-16)                 

- 
  In topSplice(), the default for `test` is now "F" instead of
  "simes".

- 
  diffSplice() is now an S3 generic function with a method for
  MArrayLM objects.

- 
  New argument `legacy` for diffSplice().

- 
  diffSplice() now adjusts transcript-level p-values using the
  full Simes method, instead of a modified Simes method leaving
  out the largest p-value for each gene that has been used until
  now. This will cause the Simes p-values returned by
  diffSplice() to be slightly more conservative than before.

- 
  diffSplice() now returns Bonferroni adjusted p-values as an
  alternative to Simes method.

- 
  The `gene.genes` data.frame returned by diffSplice() now has
  row names set to the gene IDs. Previously the last exon ID for
  each gene was used.

- 
  diffSplice() now allows for exons with df.residual=0 and
  sigma=NA.

- 
  Bug fix to diffSplice() when fitFDistUnequalDF1() was called
  with `robust=TRUE`, i.e., whenever `robust=TRUE` and the number
  of transcripts differs between genes.

- 
  Bug fix to genewise variance pooling in diffSplice() when the
  residual df are unequal. Previously the exon-wise variances
  were not weighted by df when pooling.

- 
  New argument `directional` for camera() and cameraPR().

- 
  vooma() now allows missing values in both `y` and `predictor`,
  provided there are no entirely missing rows and provided all
  non-NA y values have a non-NA predictor.

- 
  voomaLmFit() now allows missing values in `y` and `predictor`,
  provided there are no entirely missing rows and provided all
  non-NA y values have a non-NA predictor.

- 
  voomaLmFit() now allows `predictor` to be either a matrix or a
  column vector with genewise values.

- 
  Edit voomaLmFit() help page. Add Li (2024) PhD Thesis as
  reference.

[maaslin3](/packages/maaslin3)
--------

                        Changes in version 3.0.1                        

- Changed data augmentation for logistic models

- Replaced iterative renormalization with median comparisons

- Replaced group, OMP, and GOMP models with contrast tests and
  ANOVA-style comparisons

                        Changes in version 3.0.0                        

- Updated dependencies

- Modularized previous code

- Added prevalence/logistic models

- Added data augmentation for logistic models

- Added iterative renormalization for compositionally

- Added spike-in references for compositionally

- Added group, OMP, and GOMP models

[MassSpecWavelet](/packages/MassSpecWavelet)
---------------

                 Changes in version 1.73.1 (2024-12-18)                 

- Fix off-by-one error. Thanks to Aixiang Jiang, Ivan Krylov, Lluís
Revilla, Duncan Murdoch and Vincent Carey for reporting, suggesting
a fix and reaching out.

[matter](/packages/matter)
------

                        Changes in version 2.9.1                        

BUG FIXES

- Change C-level 'Free()' calls to 'R_Free()'

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                        Changes in version 1.11                         

Changes in 1.11.1

- Fix reporting of matched peaks with MatchForwardReverseParam and
gnps-matching/similarity calculation (issue #119).

[methyLImp2](/packages/methyLImp2)
----------

                 Changes in version 1.3.1 (2025-02-04)                  

- Fixing wrong counter in the loop.

[methylKit](/packages/methylKit)
---------

                       Changes in version 1.33.3                        

IMPROVEMENTS AND BUG FIXES

- remove import of "key<-" function deprecated since data.table 1.17.0

- update makeMethTabix call to remove temporary file after processing

- refactor destrand function to handle non-base resolutions and
  streamline file cleanup

                       Changes in version 1.33.2                        

IMPROVEMENTS AND BUG FIXES

- processBismarkAln: fix cigar buffer storage issue which led to faulty
  methylation calls (#334)

- remove zlibbioc dependency, zlibbioc is now deprecated in Bioc 3.21
  and will be removed in Bioc 3.22

                       Changes in version 1.33.1                        

IMPROVEMENTS AND BUG FIXES

- bedgraph: fix bedgraph error due to NA handling in meth.diff export
  (#331)

[mia](/packages/mia)
---

                        Changes in version 1.15                         

- subsetBy*: added update.tree argument

- agglomerateBy*: Add na.rm option for excluding NA counts

- Added getPERMANOVA and addPERMANOVA functions to calculate PERMANOVA

- update.tree = TRUE by default

- Add functions for classifying taxa based on their average abundance

- Add addPrevalence function that adds results of getPrevalence to
  rowData

- Added support for dimred to getCrossAssociation

- Add wrapper for PhILR transformation

- Support rarefaction when applying unifrac

- Added getMDS and addMDS: wrappers for scater::calculateMDS and
  scater::runMDS

- Fix tree merging and agglomeration

- getReducedDimAttribute: function for fetching elements from
  attributes of reducedDim

- Improve decontam functions

- Support precalculated dissimilarity matrix in dbRDA

- Calculate standard, binary Jaccard index by default


[miaViz](/packages/miaViz)
------

                        Changes in version 1.15                         

- plotAbundance: Improved visualization of sample metadata

- plotScree: Method for creating scree plots

- plotRDA: Now it is possible to add only specified vectors

- plotMediation: Method to visualise results of mediation analysis

- plotRDA: Fix error that occurred when variable names were similar

- Added plotHistogram

- Added plotBarplot

[MIRit](/packages/MIRit)
-----

                        Changes in version 1.3.1                        

Now it is possible to use a local copy of miRTarBase with the
getTargets() function.


[mitch](/packages/mitch)
-----

                       Changes in version 1.19.4                        

- New networkplot visualisation has been added to help interpret
  enrichment results.

                       Changes in version 1.19.3                        

- Fixed a breaking bug in mitch_plots() for 2D analysis.

[MLP](/packages/MLP)
---

                       Changes in version 1.55.1                        

- addGeneSetDescription: fix issue KEGG

- doc: remove description-like items in itemize (R-devel 11/10/2023)

[monaLisa](/packages/monaLisa)
--------

                       Changes in version 1.13.3                        

- switched to ggplot2 with all plotting functions (except
plotMotifHeatmaps that is still based on ComplexHeatmaps)
- new arguments introduced in plotStabilityPaths() to allow for
labels
to be shown at the end of the stability paths if desired
- changed some argument names in plotStabilityPaths() and
plotStabilityPaths(). See the help pages of these function for more
details
- updated the stability selection vignette to showcase labeling
predictors
- exposed glmnet arguments in randLassoStabSel()

                       Changes in version 1.13.2                        

- minor fix to ensure that heatmap bin legends are ordered in the
same
way as the heatmap columns

[Moonlight2R](/packages/Moonlight2R)
-----------

                        Changes in version 1.5.5                        

Summary

- fixed mistake in the documentation of the TFinfluence function

                        Changes in version 1.5.4                        

Summary

- added missing documentation file

                        Changes in version 1.5.3                        

Summary

- Filter only missense mutations of MAF in TFinfluence.
- Switched to RaSp/FoldX consensus as a default of MAVISp in
TFinfluence.
- Corrected the code of 'mutation_available' column of TFinfluence
output to correctly check the presence of the mutation in MAVISp
data.

                        Changes in version 1.5.2                        

Summary

- Added additional columns in TFinfluence output that describe the
availabilty of the TF and its mutation in MAVISp

                        Changes in version 1.5.1                        

Summary

- Added TFinfluence function, that implements a secondary layer that
checks the influence of destabilizing mutations on transcription
factors as mechanistic explanation for expression changes

                        Changes in version 1.5.0                        

Summary

- version bump due to release of Bioconductor 3.20

[mosdef](/packages/mosdef)
------

                        Changes in version 1.4.0                        

- gene_plot() defaults now to NULL in the intgroup parameter, which
translates into using the first colData item

[Motif2Site](/packages/Motif2Site)
----------

                       Changes in version 1.11.1                        

- Address the Error in Bioconductor build
- Move the source github from Bioinference to fls-bioinformatics-core


[motifStack](/packages/motifStack)
----------

                       Changes in version 1.51.1                        

- Update email address.

[motifTestR](/packages/motifTestR)
----------

                        Changes in version 1.3.5                        

- Introduced prior counts for enrichment testing

                        Changes in version 1.3.3                        

- Added simSeq for generating random sequences

[MPAC](/packages/MPAC)
----

                        Changes in version 1.1.6                        

- updated pseudocode in vignettes
- added new functions
- pltConMtf(): to plot consensus pathway submodules
- pltMtfPrtIPL(): to plot a heatmap of IPLs on proteins from
consensus pathway submodules
- pltSttKM(): to plot Kaplan-Meier curve for patient samples
stratified by given protein(s)' pathway states

                        Changes in version 1.1.5                        

- use 'Multi-omic Pathway Analysis of Cells' for MPAC

                        Changes in version 1.1.4                        

- updated text in DESCRIPTION

                        Changes in version 1.1.3                        

- Vignettes
- added examples on pltOvrHm()
- added pseudocode for MPAC
- changed MPAC's full name

                        Changes in version 1.1.2                        

- pltOvrHm(): use lowercase 'p' instead of 'P' in heatmap colorbar
legend
- MPAC's full name is changed to 'Multi-omic Pathway Analysis of
Cell'

                        Changes in version 1.1.1                        

- added new functions
- pltOvrHm(): to plot a heatmap of over-represented gene sets for
clustered samples
- ppRunPrd(): to prepare required files to run PARADIGM
- updated existing functions
- conMtf(): use decompose() instead of decompose.graph() for
igraph 2.0
- ppPermInp(): increased default permutations from 3 to 100
- ppRealInp() and ppRnaInp(): added an option rna_n_sd to set
standard deviation range to define RNA state
- pltNeiStt(): fixed bugs and adjusted output figure height
- colPermIPL() & fltByPerm(): added a threads option to
parallelize
- runPrd() & runPermPrd(): when no sampleids is provided, all the
samples in input real_se and perml will be used.
- subNtw(): entities with NA values in all samples in the input
fltmat will be ignored.
- updated tests
- test-conMtf() & test-subNtw(): use as_edgelist() instead of
get.edgelist() for igraph 2.0
- test-conMtf() & test-ovrGMT(): add upgrade_graph() on extdata

[msa](/packages/msa)
---

                       Changes in version 1.39.5                        

- workaround to avoid segfault caused by ClustalOmega

                       Changes in version 1.39.4                        

- fix of gc 8.2.8 source for Windows platform

                       Changes in version 1.39.3                        

- upgraded gc library to version 8.2.8

- changed compiler standard to gnu++11 for ClustalOmega compilation on
  Linux/Mac

- fixed a bug regarding passing file names as inputs

                       Changes in version 1.39.2                        

- further fix to cope with the upcoming change
  to C23 as default in R 4.5.0

                       Changes in version 1.39.1                        

- added USE_C17 to DESCRIPTION as SystemRequirements to cope with the
  upcoming change
  to C23 as default in R 4.5.0

[MSA2dist](/packages/MSA2dist)
--------

                 Changes in version 1.11.2 (2025-03-28)                 

- update DRESCRIPTION

                 Changes in version 1.11.1 (2024-11-14)                 

BUG FIXES

- fixed pal2nal
- fixed aastring2aabin to return alignment

NEW FEATURES

- added globalDeletionAA

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                        Changes in version 1.15                         

Changes in 1.15.3

- Fix import from MassBank txt files:
- force ordered m/z values
- skip import of deprecated records
- introduce sanity checks on spectra variables.

Changes in 1.15.2

- Import extractByIndex() from ProtGenerics.

Changes in 1.15.1

- Complete unit test coverage.

[MsBackendMetaboLights](/packages/MsBackendMetaboLights)
---------------------

                         Changes in version 1.1                         

Changes in 1.1.4

- Remove debug message.

Changes in 1.1.3

- Add mtbls_delete_cache() to delete locally cached files for a
specified MetaboLights ID.
- Change unit tests to only remove selected content instead of wiping
the full cache.

Changes in 1.1.2

- Fetch and cache each data file individually.
- Retry retrieval of data from MetaboLights up to 3 times before
throwing an error message, with an increasing time period between
attempts.

[MsBackendMgf](/packages/MsBackendMgf)
------------

                        Changes in version 1.15                         

Changes in 1.15.3

- Add support for peak annotations in MGF files.

Changes in 1.15.2

- Fix missing links in the documentation.

Changes in 1.15.1

- Complete unit test coverage.

[MsBackendMsp](/packages/MsBackendMsp)
------------

                        Changes in version 1.11                         

Changes in 1.11.1

- Complete unit test coverage.

[MsBackendSql](/packages/MsBackendSql)
------------

                         Changes in version 1.7                         

Changes in 1.7.4

- Fix handling of parameter msLevel. in filterRt().

Changes in 1.7.3

- Add new peaks data storage mode blob2 which stores the full peaks
matrix as a single entity to the database table.
- Change default peaks data storage mode to blob2.
- Add parameter peaksStorageMode to database creation function
allowing to select the new peaks data storage mode.
- Add mz() and intensity() methods.
- Small performance improvements for peaksData() using fmatch() from
the fastmatch package and avoiding to re-order the results if not
needed.
- Performance improvement for blob and blob2 storage modes by using
xda = FALSE in serialize().
- Complete unit test coverage.

Changes in 1.7.2

- Import extractByIndex from ProtGenerics.
- Add missing links in documentation.

Changes in 1.7.1

- Complete unit tests to cover all code lines.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                        Changes in version 1.19                         

MsCoreUtils 1.19.2

- Fix reduce() to merge also ranges with the same end.

MsCoreUtils 1.19.1

- Add new reduce() function to reduce overlapping numeric ranges to
disjoint ranges.

[MsDataHub](/packages/MsDataHub)
---------

                         Changes in version 1.7                         

MsDataHub 1.7.4

- Fix documentation.

MsDataHub 1.7.3

- Add cRAP contaminant databases.

MsDataHub 1.7.2

- Fix space in latest resource URLs.

MsDataHub 1.7.1

- DIA-NN report files from Ai et al (2025) for adult and iPSC-derived
cardiomyocytes single-cell proteomics data.

[MsExperiment](/packages/MsExperiment)
------------

                         Changes in version 1.9                         

MsExperiment 1.9.1

- Small fixes and update to MsBackendSql version >= 1.7.3.

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.33                         

MSnbase 2.33.5

- Remove methods-fragments.R (is now part of PSMatch).

MSnbase 2.33.4

- Fix bug in MzTab() (see issue #608).

MSnbase 2.33.3

- Fix failing unit test.

MSnbase 2.33.2

- Fix the pheplus1 code chunk in MSnbase-demo vignette to handle the
new Rdispo::getIsotope() list return value.

MSnbase 2.33.1

- Add functionality to convert a Spectra object to a MSpectra.
- Suggest pRolocdata (>= 1.43.2.1) (that has some extdata, needed to
other packages' vignettes).

MSnbase 2.33.0

- New Bioc devel version


[msqrob2](/packages/msqrob2)
-------

                       Changes in version 1.15.1                        

- fix: fixed fnames in cptac vignette and moved away from msdata to
rely on MsDataHub instead

                        Changes in version 1.15                         

                       Changes in version 1.15.0                        

- New Bioconductor devel release 3.21

                       Changes in version 1.14.1                        

- Added more flexibility in the msqrobLM function for models with
missing variables
- Fixed full rank check for models without ridge in certain cases

[MSstatsTMT](/packages/MSstatsTMT)
----------

                 Changes in version 2.14.1 (2024-10-30)                 

- Added interactive plotting functions using Plotly.

- Integrated support for ProteinProspector data conversion via
  MSstatsConvert.

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.34.0                        

New features

- Replace longFormat generic and methods with longForm (@lgatto,
#333)
- Use h5mread::h5writeDimnames instead of HDF5Array::h5writeDimnames
- Add note to warning message in saveHDF5MultiAssayExperiment to use
loadHDF5MultiAssayExperiment to reliably load data

Bug fixes and minor improvements

- Add package anchors to links in documentation
- Optimize subsetByColData (@leopoldguyot, #334)
- Use BiocBaseUtils::checkInstalled to check for Suggests packages
- Additional checks for j and ... in double bracket [[
MultiAssayExperiment method

[multistateQTL](/packages/multistateQTL)
-------------

                 Changes in version 1.99.5 (2025-03-17)                 

- Added support for additional arguments in plotting functions to be
passed to ComplexHeatmap.

                 Changes in version 1.99.0 (2025-03-07)                 

- Patches to make multistateQTL compatible with major change to
internal representations in QTLExperiment.

                 Changes in version 1.3.2 (2025-01-28)                  

- Fix up a test with a rounding error
- Remove bug in getComplete where setting verbose = TRUE would lead
to
an error
- Update replaceNAs test

                 Changes in version 1.3.1 (2025-01-14)                  

- Fix up a test that was failing due to machine tolerance

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.15.3                        

Bug fix

*Updated retrieval of IEU OpenGWAS to new approach requiring login.
Also
updated to use the IEU OpwnGWAS R package as a dependency.

                       Changes in version 1.15.2                        

New features

*FAQ Website updated.

                       Changes in version 1.15.1                        

New features

*FAQ Website added.

[mzR](/packages/mzR)
---

                       Changes in version 2.41.4                        

- Use "rtime" and "intensity" as column names for the data.frame
  returned
  by `chromatogram()`. Fixes #304.

                       Changes in version 2.41.3                        

- Fixed off-by-one indexing issue with index access for
  chromatogramHeader vs
  chromatogram. Fixes #302. Thanks Nils Hoffmann !

                       Changes in version 2.41.2                        

- Report also electron beam energy (MS:1003410) for EAD data.

                       Changes in version 2.41.1                        

- Fix compilation error with stricter compiler checks

[NADfinder](/packages/NADfinder)
---------

                       Changes in version 1.31.1                        

- Update email address.

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 3.4.0                        

- Fixed modbam_to_tabix() to not delete existing directories when
provided as the output file. Previously it would delete the
directory and create a file with the same name, now it will fail if
the output is a directory.

[NanoTube](/packages/NanoTube)
--------

                       Changes in version 1.13.1                        

- processNanostringData now reads RCC files recursively from the
parent directory, so that it can handle RCC files within subfolders.

[nullranges](/packages/nullranges)
----------

                       Changes in version 1.13.1                        

- Removing unevaluated code in vignette

[OmnipathR](/packages/OmnipathR)
---------

                        Changes in version 4.0.0                        

Features

- ID translation ambiguity analysis
- KEGG REST API client
- New metabolomics resources: HMDB, RaMP, Chalmers Sysbio GEMs,
STITCH,
- Metabolite-protein interactions from MetalinksDB
- Metabolite ID translation using RaMP and HMDB data

Technical

- Rewritten and improved parameter processing for OmniPath queries
- Rewritten downloaders based on curl and httr2
- Fine control over curl handlers
- Detailed log messages about HTTP requests
- Diagnostic log messages: session info, libraries, versions,
platform
info, curl options, HTTP timings, headers
- Option to include curl debug log in the log file
- Robust TCP keep-alive parameters that hopefully fix
rest.uniprot.org
dowloads

[OncoSimulR](/packages/OncoSimulR)
----------

                 Changes in version 4.9.1 (2025-03-04)                  

- Fix three Notes in Rmd files.

[ontoProc](/packages/ontoProc)
--------

                        Changes in version 2.1.7                        

- update the module dependencies in basilisk.R

- add quickOnto() for convenient retrieval of OIRDS serialized
  ontologyIndex instances in cache

[OutSplice](/packages/OutSplice)
---------

                 Changes in version 1.7.1 (2025-02-03)                  

- Added warning message for only having 1 Tumor Sample

- Updated vignette file

[pairedGSEA](/packages/pairedGSEA)
----------

                        Changes in version 1.7.0                        

- Ensure compatibility with updated msigdbr and new R version
- Fixed rare bug in aggregate_pvalue

[pathview](/packages/pathview)
--------

                       Changes in version 1.47.1                        

- fix "R CMD check" error caused by some problem in org.EcK12.eg.db

- remove org.EcK12.eg.db from Suggests for BioC build check

- korg now include 10718 KEGG species or 951 new species beyond March
  2024.

[pcaExplorer](/packages/pcaExplorer)
-----------

                        Changes in version 3.2.0                        

Other notes

- Some changes in the source of the documentation, providing anchors
to all function calls. This avoids the new note in R CMD check in
the new major release

[Pedixplorer](/packages/Pedixplorer)
-----------

                        Changes in version 1.3.4                        

- Add example of interactivness in vignette
- Fix label adjusting position in plot functions
- Fix arrow size in ggplot

                        Changes in version 1.3.1                        

- Add support for .ped, .tsv files in data import
- is_informative independent from useful_inds
- Use directly columns names from fill instead of the mods columns
- Move to R 4.4 and Bioc 3.20
- Fix unittests and update snapshots
- Change normalisation process to directly use id, dadid, momid,
famid, sex no more need for indId, fatherId, ...
- affection is now used as default affection modality columns that
will be used to generate affected
- status is replaced by deceased
- steril is replaced by fertility and corresponding symbols is added
for infertile and infertile_choice_na
- terminated sex code is replace by miscarriage new slot
- miscarriage, evaluated, consultand, proband, carrier, asymptomatic
and adopted are now recognize and use for plotting
- Argument order of Ped() as changed when using vectors. This choice
has been made for a better consistency across the package. Please
check that your argument are properly named (i.e. sex has been moved
after famid and avail after deceased).
- Shiny application is updated and improved (aesthetics, errors,
warnings, functionnalities).
- Add dateofbirth and dateofdeath to the Ped object
- Changee from round to signif for the precision argument
- Improve stability of test by adding and controlling the options()
and par() arguments in the unittests.
- Carrier symbols is proportional to the mean of the box size

[pgxRpi](/packages/pgxRpi)
------

                 Changes in version 1.3.3 (2025-02-12)                  

- Updated the "sample_count" extraction method from
"services/collations" to beacon count response, expanding counts to
include all available entities (analyses, biosamples, individuals)
and changing the type from "sample_count" to "counts" in pgxLoader.
- Enabled parallel queries across multiple resource domains (in
pgxmetaLoader).
- Optimized code for Beacon response mapping and updated YAML mapping
rules.

- Updated pgxfreqplot to set ylim as max(lowfreq, highfreq), ensuring
proper scaling when high-level frequency exceeds low-level
frequency.
- Replaced id parameter with the standard query path for Beacon
queries.
- Modified variant data retrieval in seg format to fetch directly
from
the service API.

                 Changes in version 1.3.1 (2025-01-14)                  

- Removed the pgxFilter function and incorporated its functionality
into pgxLoader with the "filtering_terms" type, following the Beacon
v2 response mapping.
- Updated vignette filenames for improved clarity.
- Made parameter checks more flexible for Beacon queries.

[PhyloProfile](/packages/PhyloProfile)
------------

                        Changes in version 2.0.0                        

- runtime significantly improved #105

                       Changes in version 1.20.4                        

- added t-SNE dimension reduction

- highlight genes and taxa for sub-plot

- use fastcluster for clustering and sorting the input taxa

                       Changes in version 1.20.3                        

- fixed bug with NCBI taxonomy link at the first time running

- improve UMAP 3D (hover info, dot zooming, same colors with 2D
  plot)

                       Changes in version 1.20.2                        

- hover option for UMAP

- UMAP 3D using plotly

- fixed bug with Point's info

                       Changes in version 1.20.1                        

- auto identify refspec

- fixed subsetting data using a list of gene IDs

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.15.1 (2024-12-19)                 

- Modified results table

[POMA](/packages/POMA)
----

                       Changes in version 1.17.6                        

- Allow grouping variable selection in PomaImpute
- New documentation

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.47                         

Changes in version 1.47.5

- Remove long title in mcmc_plot_probs

Changes in version 1.47.4

- Addition of UMAP as a dimensionality reduction method in plot2D
- In this version there have been several new arguments to the plot2D
function to improve visualisation. Please see below.
- A new argument bg has been added to plot2D to allow the use of
filled point characters see issue #153
- plot2D also has a new argument called palette to facilitate the
generation of nice colour schemes to work with filled point
characters
- The functions setStockbg, getStockbg, setUnknownbg and getUnknownbg
have been added for colour management
- The behaviour of the function addLegend when called with argument
where = "other". The call to dev.new() has been removed so that
instead of opening a new device and plotting the legend it instead
plots the legend on an empty plot in the current window.
- The default of grid in plot2D has been changed to FALSE.
- All plotting arguments now have an argument unknown for the
specification of how unlabelled points are defined. Default is still
"unknown".

Changes in version 1.47.3

- Deprecate setAnnotationParam, getAnnotationParams,
showGOEvidenceCodes, getGOEvidenceCodes, addGoAnnotations,
orderGoAnnotations and associated helper functions
- Remove vignette v04 for adding GO annotations
- Update getting started vignette v01 with clustDist (previously in
v04)
- Update the transfer learning vignettes to reflect the deprecation
of
GO functions

Changes in version 1.47.2

- Bump version to propage to landing page

Changes in version 1.47.1

- Lisa Breckels is now maintainer.

Changes in version 1.47.0

- New devel version

[ProtGenerics](/packages/ProtGenerics)
------------

                       Changes in version 1.39.2                        

- Move `applyProcessing()`, `processingChunksize()`,
  `processingChunkSize()<-`
  and `processingChunkFactor()` generics from the Spectra package.

                       Changes in version 1.39.1                        

- Add `extractByIndex()` generic.

[PSMatch](/packages/PSMatch)
-------

                        Changes in version 1.11                         

PSMatch 1.11.5

- Correct plotSpectraPTM annotations.

PSMatch 1.11.4

- Added plotSpectraPTM: a plotting function to visualise
post-translational modifications.

PSMatch 1.11.3

- Deprecated addFragments. The use of labelFragments is endorsed
instead. See PR #20.

PSMatch 1.11.2

- Replace calculateFragments with calculateFragments2. See PR #19.

PSMatch 1.11.1

- New calculateFragments2 function includes fixed and variable
modifications to fragments ions. See PR #16.

PSMatch 1.11.0

- New devel version

[PureCN](/packages/PureCN)
------

                       Changes in version 2.14.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Do not filter allosome coverage for male samples.

BUGFIXES

- Fixed large runtime of readAllelicCounts. Thanks luyh-xp (#378).

- --bai was ignored in Coverage.R except for lists (#272). Thanks
  @nspies-carisls and @lbeltrame.

[pwalign](/packages/pwalign)
-------

                        Changes in version 1.4.0                        

- No significant changes in this version.

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.17                         

QFeatures 1.17.5

- Remove superfluous message when filtering with keep = TRUE (see
issue 231).

QFeatures 1.17.4

- Optimisation of aggregateFeatures in the case of multiple assays
aggregation.

QFeatures 1.17.3

- Fix bug in nNA() on empty assays (see #174).

QFeatures 1.17.2

- Fix fnames argument in readQFeatures(). The fnames argument is not
passed to readSummarizedExperiment() anymore, but used at in
readQFeatures(). Rownames are now set after splitting. Given that
rownames must be unique and that this was enforced with
make.unique(), the previous behaviour misnamed features that should
get the same name.

QFeatures 1.17.1

- Fix: solved readQFeatures() bug and adapted unit tests.
- Docs: created new vignette about reading data with QFeatures

QFeatures 1.17.0

- New Bioconductor 3.21 (devel) release


[QTLExperiment](/packages/QTLExperiment)
-------------

                 Changes in version 1.99.1 (2024-03-05)                 

- Major update to the underlying representation of the slots. o Slots
now use elementMetadata and colData. o No longer use internal col
and row Data. o Almost all functions have been updated and testing
is working but there may still be some edge cases that I have
missed. o Better consistency with other child classes of
SingleCellExperiment.

                 Changes in version 1.5.0 (2024-10-30)                  

- Added a test for object validity error messages

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.32.0                        

Bug fixes and minor improvements

- Added package anchors to links in documentation

[RAIDS](/packages/RAIDS)
-----

                        Changes in version 1.5.1                        

SIGNIFICANT USER-VISIBLE CHANGES

o Intagration of Rsamtools

[ramr](/packages/ramr)
----

                 Changes in version 1.15.2 (2025-03-06)                 

- complete rewrite of getAMR and simulateData for speed and robustness

                 Changes in version 1.15.1 (2024-12-17)                 

- zero-and-one inflated beta distribution by gamlss.dist::BEINF

- cleaner plots

[Rarr](/packages/Rarr)
----

                         Changes in version 1.7                         

- Added path() method for ZarrArray class that returns the location
of
the zarr array root.
- Removed used of non-API call SETLENGTH in C code.
- Small changes to compilation of internal blosc libraries to cope
with the C23 compiler becoming the default in R-4.5.0

[rawDiag](/packages/rawDiag)
-------

                         Changes in version 1.3                         

- adapt to rawrr changes

[rawrr](/packages/rawrr)
-----

                        Changes in version 1.15                         

- Add Leo as author.

- Replace `mono` with `dotnet`; Introduced a statically linked
  runtime `rawrr` assembly for Windows, Linux, and macOS on x64
  architectures, streamlining installation and usage #75.

- Consistently normalize paths #76.

[RCy3](/packages/RCy3)
----

                       Changes in version 2.28.0                        

- No changes since version 2.26.0.

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.21.1 (2024-11-07)                 

- Added new functions to perform quantitative comparative pathway
  analyses for scRNA-seq data.

[RegionalST](/packages/RegionalST)
----------

                        Changes in version 1.5.6                        

- Fix bugs when plotting figures with BayesSpace

[RESOLVE](/packages/RESOLVE)
-------

                       Changes in version 1.10.0                        

- Major code update, implementing signatures-based clustering and
  associations to signatures exposures.

- Package released in Bioconductor 3.21.

[rGREAT](/packages/rGREAT)
------

                        Changes in version 2.9.1                        

- change NCBI api to v2.

- use `exp(phyper(..., log.p = TRUE))` to get rid of zero p-values.
  Thank @andvon
  for the code contribution.

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.52.0                        

CHANGES

- Removed used of the non-API `STRING_PTR` function when
  permuting string datasets.

Bug fixes

- Fixed overflow when trying to write a large number of very long
  strings.  Backported to 2.50.1.  (Thanks to Aaron Lun @LTLA for
  reporting this https://github.com/grimbough/rhdf5/issues/151).

- Addressed potential issued caused by undefined compiler
  behaviour when arrays of length zero were being created in C
  code.  Backported to 2.50.2 (Thanks to Prof B. Ripley for
  reporing this.)

[rhdf5filters](/packages/rhdf5filters)
------------

                       Changes in version 1.20.0                        

- Small changes to compilation of internal blosc and vbz libraries to
cope with the C23 compiler becoming the default in R-4.5.0

BUG FIXES

- Ensure PKG_CPPFLAGS are passed to all compilation step. This was an
issue in the case where R was installed via conda and package
binaries were not available. Backported to version 1.18.1.

[Rhtslib](/packages/Rhtslib)
-------

                        Changes in version 3.4.0                        

- No significant changes in this version.

[RiboDiPA](/packages/RiboDiPA)
--------

                       Changes in version 1.15.1                        

- Rcpp function centerCPP.cpp was modified so it can take reads aligned
  to more than 10 consecutive regions in the genome. Ribo-seq read
  should typically align to no more than three consecutive blocks in
  the genome.

[ribosomeProfilingQC](/packages/ribosomeProfilingQC)
-------------------

                       Changes in version 1.19.2                        

- Export the reads distribution with considering the precedence.

                       Changes in version 1.19.1                        

- Update email address.

[ROTS](/packages/ROTS)
----

                        Changes in version 2.0.0                        

- Added support for linear and mixed-effects models

[rpx](/packages/rpx)
---

                        Changes in version 2.15                         

rpx 2.15.1

- Fix error due to remote side change (another one).

rpx 2.15.0

- New devel version (Bioc 3.21)

[Rsubread](/packages/Rsubread)
--------

                       Changes in version 2.22.0                        

- Improved identification and reporting of junction counting
  functionality in featureCounts.

[rWikiPathways](/packages/rWikiPathways)
-------------

                       Changes in version 1.28.0                        

- None

[S4Arrays](/packages/S4Arrays)
--------

                        Changes in version 1.8.0                        

NEW FEATURES

- arep_times() and arep_each(): multidimensional versions of
  base::rep( , times=) and base::rep( , each=).

- as_tile(): a convenient way to control the direction of recycling in
  the context of arithmetic and other binary operations between an
  array
  and a vector, or between two arrays of distinct dimensions.
  Still work-in-progress!

- kronecker() methods for Array objects that work out-of-the-box on
  Array
  derivatives that support [ and *.

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.46.0                        

NEW FEATURES

- Add 'margin' argument to wmsg() internal helper.

SIGNIFICANT USER-VISIBLE CHANGES

- Expose ranges_mapper() and positions_mapper() thru the S4Vectors C
  interface.

[SCArray](/packages/SCArray)
-------

                       Changes in version 1.16.0                        

- update according to the SparseArray package

[SCArray.sat](/packages/SCArray.sat)
-----------

                        Changes in version 1.8.0                        

- update according to DelayedArray & SparseArray

[scater](/packages/scater)
------

                       Changes in version 1.36.0                        

- Add min.value,max.value arguments to plotReducedDim to enable
  truncating colour scales using a numeric value or a quantile
  (eg "q10").

[scDotPlot](/packages/scDotPlot)
---------

                        Changes in version 1.1.1                        

- Enabled clustering for a group with annotations for a single
feature


[scoup](/packages/scoup)
-----

                        Changes in version 1.1.1                        

- Improved flexibility by allowing specification of mutation rate
  parameters.

- Adjusted the `fixed` input in `codonCoeffs` function.

- Reversed the mutation matrix flexibility due to expensive
  computational cost.

- Updated help files to aptly emphasise episodic in-place of
  "deterministic".

- Revised the `codonFreq` function so that it uses `softmax` identity.

[scp](/packages/scp)
---

                        Changes in version 1.17                         

scp 1.17.1

- Deprecate aggregateFeaturesOverAssays, use
QFeatures::aggregateFeatures instead.

scp 1.17.0

- New Bioconductor 3.21 devel

[scrapper](/packages/scrapper)
--------

                         Changes in version 1.2                         

- Added the aggregateAcrossGenes() function, to compute an
  aggregate expression value for gene sets.

- Added compute.cohens.d=, compute.delta.mean= and
  compute.delta.detected= options to scoreMarkers().

- Support top=Inf in chooseHighlyVariableGenes(). Also added the
  bound= argument to set a hard upper/lower bound.

- Bugfix for correct filtering with block= in the various
  filter*QcMetrics() functions when not all blocking levels are
  present.

- Bugfix to clusterKmeans() to respect the user-supplied seed= in
  relevant initialization methods.

- Added a return.graph= option to return the SNN graph from
  runAllNeighborSteps().

- Added a testEnrichment() function for quick and dirty gene set
  enrichment testing.

- Modified runPca() so that it caps number= to the maximum number
  of available PCs.

- Added an analyze() function that provides a one-click approach
  for analyzing single-cell data.

- Added a reportGroupMarkerStatistics() function to combine all
  marker statistics for a single group into one data frame.

- Switch clusterGraph() to use C++ wrappers around the igraph
  community detection algorithms via Rigraphlib. This replaces
  the dependency on the igraph R package.

- Added a delayed= option to avoid wrapping the matrix in a
  DelayedArray in normalizeCounts().

[scRepertoire](/packages/scRepertoire)
------------

                        Changes in version 2.3.2                        

UNDERLYING CHANGES

- Fixed issue with denominator in getCirclize()
- Fixing chain issue with clonalCompare() - expanded assertthat
statement

                        Changes in version 2.2.1                        

- Rebasing for the purposes of Bioconductor version 2.2.0

NEW FEATURES

- Added support for BCRs for loading ParseBio sequences.
- Added quietBCRgenes() and quietTCRgenes() from Ibex and Trex and
quietVDJgenes() as a convenience that runs both. The functions
filter out known TCR and/or BCR gene signatures.

UNDERLYING CHANGES

- Added Seurat to the Suggests field in the DESCRIPTION file.

[scRNAseqApp](/packages/scRNAseqApp)
-----------

                       Changes in version 1.7.14                        

- Add function extractFragmentNameMapList.

                       Changes in version 1.7.13                        

- Add search for token.

                       Changes in version 1.7.12                        

- Fix the token issue.

                       Changes in version 1.7.11                        

- Fix the issue if the data list is not updated for the main page.

                       Changes in version 1.7.10                        

- Fix multiple typos.

                        Changes in version 1.7.9                        

- Update the footnote.

                        Changes in version 1.7.8                        

- Remove the import of X11Fonts.

                        Changes in version 1.7.7                        

- Update email address.

                        Changes in version 1.7.6                        

- Add fragmentNameMap to createDataset.

- Add fragmentNameMap to createDataset.

                        Changes in version 1.7.5                        

- Optimize read and save ATAC peak matrix.

                        Changes in version 1.7.4                        

- Add module ExplorerAcc.

                        Changes in version 1.7.3                        

- check parameters before creating dataset.

                        Changes in version 1.7.2                        

- add binSize parameter for atac bigwig files.

                        Changes in version 1.7.1                        

- fix the bug if no signals for atac coverage and invalid filenames.

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.48.0                        

NEW FEATURES

- `seqAddValue()`: use bit1 for a logical vector; new argument
  'use_float32' for storing double

- new argument 'start' in `seqResetVariantID()`

- new argument 'digest' in `seqRecompress()` to add MD5 hash codes

- `seqGetData(, "$chromosome")` returns chromosome codes in an object
  of
  'S4Vectors::Rle'

- `seqGetData(, .tolist=NA)` returns an extended list defined in
  IRanges
  (e.g., IntegerList) when it is applicable

- `seqListVarData(, useList=TRUE)` returns an extended list defined in
  IRanges

UTILITIES

- Tweak display in `seqResetVariantID()`

- use `crayon::silver()` instead of `crayon::blurred()` in the display
  since
  RStudio blurs the screen output

BUG FIXES

- `seqBlockApply()` should recover the filter when the user-specified
  function fails

                       Changes in version 1.46.2                        

UTILITIES

- `seqVCF_Header()` allows multiple cores to calculate the total number
  of
  variants when 'getnum=TRUE' (the Rsamtools package should be
  installed);
  `seqVCF2GDS()` is faster when obtaining the number of variants for
  splitting files.

- new 'variant_count' in `seqVCF2GDS()` to specify the number of
  variants
  in the VCF file when it is known or an approximation is known; it is
  only
  applicable when multiple cores are used. If 'variant_count' is
  specified,
  counting the number of variants will be skipped.


[SingleCellAlleleExperiment](/packages/SingleCellAlleleExperiment)
--------------------------

                        Changes in version 1.3.2                        

Fixing error in vignette regarding the knee plot in 3.21 devel / R4.5
Removing the fitted line from knee plot

                        Changes in version 1.3.1                        

Fixing error in vignette regarding the knee plot in 3.21 devel

                        Changes in version 1.3.0                        

No changes, devel version bump due to new bioconductor 3.21 devel
branch

[SingleCellExperiment](/packages/SingleCellExperiment)
--------------------

                       Changes in version 1.30.0                        

- Preserve NAs in the input matrices to rowPairs() and colPairs()
  converting them to SelfHits.

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.16.1 (2024-02-12)                 

- Fixed decontX/soupX webapp error
- Fixed error in reportCellQC
- Fixed additional bugs with Seurat V5 integration

[SingleR](/packages/SingleR)
-------

                       Changes in version 2.10.0                        

- Switch to scrapper for DE detection when de.method="t" or
  de.method="wilcox" in trainSingleR(). This should give similar
  results to but is faster than the previous scran functions.

- Switch to scrapper for variance calculation, PCA and clustering
  in aggregateReference(). This should be faster than the
  previous BiocSingular and stats::kmeans functions, and avoids
  the need to consider the seed.

- Added a configureMarkerHeatmap() function to perform all the
  calculations used by plotMarkerHeatmap(). This allows users to
  re-use the calculations with a custom visualization for the
  expression values.

- Automatically remove duplicated gene names in trainSingleR().
  This avoids matching to the wrong gene after identifying
  markers from the reference dataset.

- Report scores as a DataFrame of nested Dataframes in
  combineRecomputedResults(). Each inner DataFrame corresponds to
  a reference and contains the identity of the best label and the
  recomputed score in that reference. This is simpler and more
  efficient than the previous "expanded with NA" format.

- Report the deltas (i.e., difference between the best and
  next-best scores) in combineRecomputedResults().

- Separate the missingness check arguments in SingleR() with the
  new check.missing.test= and check.missing.ref= options. The
  former is disabled by default, to avoid an unnecessary
  missingness check in the vast majority of test cases.

- Added a fine.tune.combined= option to classifySingleR(), to
  disable fine-tuning only when combining multiple references.
  This allows users to recover pre-2.8.0 behavior and trades
  accuracy for speed.

- Removed the deprecated combineCommonResults() function.


[smartid](/packages/smartid)
-------

                        Changes in version 1.3.2                        

- Update batch param in top_markers function.

                        Changes in version 1.3.1                        

- Update top_markers function to allow batch correction for glm
method.

[sparrow](/packages/sparrow)
-------

                        Changes in version 1.14                         

Enhancements

- Updated to support refactored msigdbr package (version >= 10.0.0),
which now only includse the human Hallmark genesets internally. The
msigdbr package relies on an external (non CRAN) msigdbdf package,
for the rest of the (now more recent) definitions of the MSigDB gene
set collections.

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

[SparseArray](/packages/SparseArray)
-----------

                        Changes in version 1.8.0                        

NEW FEATURES

- Make rowsum() and colsum() work on any dsparseMatrix derivative and
  not
  just on a dgCMatrix object. For example now they also accept a
  dgRMatrix
  or dgTMatrix object.

- Make crossprod(), tcrossprod(), and %*% work between a SparseMatrix
  object and an ordinary vector.

- Arith operations between a SparseArray (or NaArray) object 'x' and an
  atomic vector 'y' are no longer restricted to the latter being a
  single
  value: now the length of the latter can also be 'dim(x)[[1]]' or a
  divisor of it.

BUG FIXES

- Make sure that SparseArray(x) and SVT_SparseArray(x) always propagate
  the
  dimnames.

- Fix bug in Math and Math2 operations on SVT_SparseArray and NaArray
  objects.

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 2.13.1 (2024-11-23)                 

- Added project wibsite in vignettes.

[Spectra](/packages/Spectra)
-------

                        Changes in version 1.17                         

Change in 1.17.10

- Accept labels argument as a list instead of a character in the
plotting functions.

Change in 1.17.9

- Allow parameter msLevel. = integer() for filterRt() to filter
spectra of all MS levels. This was msLevel. = uniqueMsLevels(),
which, depending on the backend, can be computationally intense. Add
related unit tests to the unit test suite.

Change in 1.17.8

- Add parameter return.type to peaksData().

Change in 1.17.7

- Add the spectraVariableMapping<- generic method.

Change in 1.17.6

- Add new fillCoreSpectraVariables() function that allows to add
eventually missing core spectra variables (with the correct data
type) to a data frame.

Change in 1.17.5

- Move generics processingChunkSize(), processingChunkFactor() and
applyproceesing() to ProtGenerics. Required ProtGenerics version
1.39.2 or higher. These were moved to be able to implement them in
the Chromatograms package.

Change in 1.17.4

- Import extractByIndex() from ProtGenerics.

Change in 1.17.3

- Fix cbind2() unit test for backends that fails if the number of
spectra in the tested backend is (by chance) equal to 4.

Change in 1.17.2

- Add cbind2() method to easily add multiple spectraVariables and
their content to the spectraData of a Spectra object. See also issue
#342

Changes in 1.17.1

- Refactor containsMz() to support chunk-wise processing.

[SPIAT](/packages/SPIAT)
-----

                        Changes in version 1.9.0                        

Development version on Bioconductor 3.21.

[SpotSweeper](/packages/SpotSweeper)
-----------

                        Changes in version 1.3.2                        

Minor Changes

- Broadened Compatibility: Updated all functions to use inherits(spe,
"SpatialExperiment") instead of checking class(spe) directly. This
change ensures that derived classes (e.g., SpatialFeatureExperiment)
are also supported, improving flexibility and ease of use.

New Features

- Added the 'flagVisiumOutliers()' function to identify and flag
systematic outlier spots in Visium datasets. This feature enhances
data quality by allowing users to efficiently detect and exclude
problematic spots from downstream analyses.

                        Changes in version 1.3.1                        

Major Changes

- Function Renaming: The function plotQC has been renamed to
plotQCmetrics to better reflect its purpose. The new function
plotQCmetrics should be used moving forward. This change improves
clarity in the package’s API by specifying that this function is
designed for plotting QC metrics.

New Features and Enhancements

- shape argument: Added a shape argument to findArtifacts, allowing
users to specify the neighborhood shape as either "hexagonal" or
"square" for local variance calculations. This enhancement provides
flexibility for different spatial arrangements in spatial
transcriptomics data.

- Updated n_order parameter: Renamed the n_rings parameter to n_order
in the findArtifacts function to better describe its purpose of
specifying the N-order neighbors for local mitochondrial variance
calculations.

- Parallelization: Added a workers argument for parallel processing
using BiocParallel in both localOutlier and localVariance functions.
This allows for faster computation, particularly on larger datasets.

Deprecations

- plotQC Function Deprecated: The plotQC function is now deprecated.
While it remains available for backward compatibility, users are
encouraged to transition to plotQCmetrics. Calling plotQC will
display a warning, reminding users of the deprecation.

This change is backward compatible; existing code using plotQC will
still work but will show a warning. We recommend updating your code
to
use plotQCmetrics to avoid any issues in future versions where plotQC
may be removed.

[sRACIPE](/packages/sRACIPE)
-------

                       Changes in version 1.99.0                        

- New interaction types added

- Convergence testing for deterministic simulations added

- Ornstein-Uhlenbeck noise option added

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

[structToolbox](/packages/structToolbox)
-------------

                       Changes in version 1.19.2                        

- fix pca_biplot

                       Changes in version 1.19.1                        

- fix ANOVA output with interaction terms

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.38.0                        

- No changes in this version.


[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.19.11                       

- Fixes bug where ExoLabel wouldn't handle multiple files correctly
- ExoLabel now reads in networks roughly twice as fast
- ExoLabel now correctly skips negative weighted edges
- Various improvements to ExoLabel output formatting
- Enables experimental tuning of distance scaling in ExoLabel
attenuation

                       Changes in version 1.19.10                       

- Fixes critical bug in ExoLabel that previously caused large
networks
to always report a single large cluster
- ExoLabel no longer allows negative weights to become positive via
attenuation (positive weights can still become negative)

                       Changes in version 1.19.9                        

- Bugfixes to ExoLabel
- ExoLabel now defaults to use_fast_sort=TRUE
- Various QoL changes for ExoLabel(..., return_table=TRUE)
- ExoLabel now correctly reports timing

                       Changes in version 1.19.8                        

- ExoLabel has a new parameter to tune the performance of hop-length
attenuation.
- Documentation and formatting updates.

                       Changes in version 1.19.7                        

- ExoLabel now uses hop-length attenuation to mitigate formation of
massive communities.
- ExoLabel no longer supports inflation, since attentuation does a
better job handling this without introducing additional parameters.
- ExoLabel will now print a lot less when running with verbose=TRUE
in
non-interactive mode -- just as informative, but less junk caused by
lots of unrendered carriage returns.

                       Changes in version 1.19.6                        

- fixes multiple bugs in EstimRearrScen

                       Changes in version 1.19.5                        

- fixes bug in GeneVector.EvoWeaver that could affect DNA-based
analyses

                       Changes in version 1.19.4                        

- fixes bug that prevented building on Windows
- adds multiple clustering support for ExoLabel

                       Changes in version 1.19.3                        

- ExoLabel will no longer crash when given a network lacking a
trailing newline
- Various internal improvements and code refinements

                       Changes in version 1.19.2                        

- Lots of bug fixes for ExoLabel
- ExoLabel now reports disk consumption during execution

                       Changes in version 1.19.1                        

- ExoLabel is even faster due to in-place external sort for faster
file I/O
- Other quality of life improvements to ExoLabel

                       Changes in version 1.19.0                        

- First development version of Bioconductor 3.21

[systemPipeShiny](/packages/systemPipeShiny)
---------------

                       Changes in version 1.18.0                        

Major Change

- Update to systemPipeR 2.12 and systemPipeRdata 2.10
- Add support for Single Cell RNAseq, BLAST, and Cheminformatics
workflows which are new in systemPipeR.

Minor Change

- use esquisse_container instead of esquisseContainer due to the new
version of {esquisse}.

Bug Fix

- Fix shiny 1.10 update that caused some buttons not working.
- Fix welcome page logo animation not working, now remove the
animation.
- Fix checkModulePkgs_internal had wrong logic when all packages are
installed. Now return the correct value.

[tadar](/packages/tadar)
-----

                        Changes in version 1.6.0                        

- Deprecated dar_val arg of assignFeatureDar()

[TargetSearch](/packages/TargetSearch)
------------

                       Changes in version 2.10.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The method `libId` is deprecated. Use the fully equivalent function
  `makeIndex`. This method was mostly used internally, so it should not
  cause disruption.

NEW FEATURES

- New method `libUID` which sets or gets the metabolite library unique
  identifiers. Previously, there was no method to change them after
  import.

BUG FIXES

- Use `SET_INTEGER_ELT` to set an integer in the C code. The wrong
  macro was used for unknown reasons. This generated an error on newer
  R.

INTERNAL

- Replace the memory allocation functions due to STRICT_R_HEADERS being
  enabled (Calloc, Free => R_Calloc, R_Free)

[TBSignatureProfiler](/packages/TBSignatureProfiler)
-------------------

                       Changes in version 1.19.0                        

Minor Changes

- Updated gene name update mechanism in runTBsigProfiler to sum genes
with identical names

Bug Fixes

- Fixed issue with GSVA needing a non-data.frame object in
runTBsigProfiler


[TENxIO](/packages/TENxIO)
------

                       Changes in version 1.10.0                        

Bug fixes and minor improvements

- Removed projection from the show method
- Fixed a missing anchor in @inheritParams documentation
- The TENxH5 constructor properly validates remote files with
.validateRanges
- imported slotNames from the methods package

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

[topconfects](/packages/topconfects)
-----------

                       Changes in version 1.23.2                        

- Signed color palette in plot_confects_me2().

                       Changes in version 1.23.1                        

- Add plot_confects_me2().

[topdownr](/packages/topdownr)
--------

                       Changes in version 1.29.1                        

- Adapt .calculateFragments to changes in PSMatch::calculateFragments
introduced in PSMatch:PR19.
- Replace partial matching arguments with complete argument names.

                        Changes in version 1.29                         

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.43.7                        

- Fix a potential bug in reading sparse cool file.

                       Changes in version 1.43.6                        

- Update email address.

                       Changes in version 1.43.5                        

- Add line and histogram to tracktype for trackstlye.

                       Changes in version 1.43.4                        

- Remove the clip for annotation track.

                       Changes in version 1.43.3                        

- Add documentation to remove y label.

                       Changes in version 1.43.2                        

- Make it possible to plot annotation for peak regions.

                       Changes in version 1.43.1                        

- Make it possible to plot link for back2back interaction data.

[transmogR](/packages/transmogR)
---------

                        Changes in version 1.3.8                        

- Added shiftByVar() to produce shifted co-ordinates which match
those
after incorporation of variants

                        Changes in version 1.3.1                        

- Changed default behaviour of digestSalmon() to exclude 'TPM' and
'effectiveLength' assays, which are now optional via the
'extra_assays' argument


[tximeta](/packages/tximeta)
-------

                       Changes in version 1.25.1                        

- Added skipRanges to summarizeToGene which allows summarization of
assay data when skipMeta was used and therefore ranges should not be
used / output in summarization.

[UCSC.utils](/packages/UCSC.utils)
----------

                        Changes in version 1.4.0                        

- No significant changes in this version.

[UniProt.ws](/packages/UniProt.ws)
----------

                       Changes in version 2.48.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Switched to using the `httr2` package for making HTTP requests

- Reorganized and renamed documentation for the `UniProt.ws` class and
  its
  methods

- Converted unit testing framework from `testthat` to `tinytest`

- Added the `README.md` file from the vignette

- Improved handling of API responses: empty results with code 200 are
  now caught and handled properly

- The `collapse` argument in `queryUniProt` now shows both `OR` and
  `AND`
  options

BUG FIXES AND MINOR IMPROVEMENTS

- Added an example to `queryUniProt()`

- Fixed formatting issues in unit tests and corrected indentation in
  code

- Adjusted internal logic to support queries returning larger result
  sets

- Re-added `Sys.sleep()` call to respect rate limits

- Set `n` argument in internal use of `head()`

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.24.2                        

BUG FIXES

- Incorrect version number in last NEWS update.

                       Changes in version 1.24.1                        

BUG FIXES

- Fixed another missed Rcpp bounds check warning.

[updateObject](/packages/updateObject)
------------

                       Changes in version 1.12.0                        

- No changes in this version.


[ViSEAGO](/packages/ViSEAGO)
-------

                        Changes in version 1.21                         

- GOterms_heatmap static image update

- GOclusters_heatmap static image update

- show_heatmap static image update

[VisiumIO](/packages/VisiumIO)
--------

                        Changes in version 1.4.0                        

New features

- Added support for CytAssist images (@ZheFrench, #8)
- Enabled use of alternative readers via VisiumIO.csvreader option
- Included metadata in SpatialExperiment outputs

Bug fixes

- Resolved issue where *_feature_bc_matrix folder checks were
incorrectly triggered for non-mtx formats (@estellad, #4)
- Fixed missing package anchor in @param documentation entries

[xcms](/packages/xcms)
----

                         Changes in version 4.5                         

Changes in version 4.5.4

- Replace usage of deprecated (and removed) class NAnnotatedDataFrame
with AnnotatedDataFrame.
- Fix a bug in manualChromPeaks() that caused an error when only a
single chrom peak was added.

Changes in version 4.5.3

- Address issue #765: peak detection on chromatographic data: report
a
chromatogram's "mz", "mzmin" and "mzmax" as the mean m/z and lower
and upper m/z in the chromPeaks() matrix.
- Fix calculation of the correlation coefficient for peak shape
similarity with an idealized bell shape (beta) during gap filling
for centWave-based chromatographic peak detection with parameter
verboseBetaColumns = TRUE.
- Add chromPeakSummary generic (issue #705).
- Add chromPeakSummary() method to calculate the beta quality
metrics.
- Add c() method to combine multiple XcmsExperiment objects into one.
- Add a method to coerce from XCMSnExp to XcmsExperiment objects.

Changes in version 4.5.2

- Small update to featureSpectra() and chromPeakSpectra() to allow
addition of chromPeaks() and featuresDefinitions() columns to be
added to the Spectra output.
- Tidied the xcms vignette, to order the filtering of features and
remove the outdated normalisation paragraph.In depth discussion on
this subject can be found on metabonaut.

Changes in version 4.5.1

- Fix compile errors with R-4.5

[Xeva](/packages/Xeva)
----

                       Changes in version 1.23.1                        

- The datasets previously available have been deprecated.
  `downloadXevaSet` will no longer work.

[zellkonverter](/packages/zellkonverter)
-------------

                       Changes in version 1.18.0                        

New features

- 
  Add minimal support for SpatialExperiment objects to
  writeH5AD() and SCE2AnnData(). This stores the spatial
  coordinates in a obsm item named "spatial" as expected by the
  *squidpy* Python package. (PR from @mcmero)

Major changes

- 
  Add environment for *anndata* v0.11.4. This is now the default
  environment for the Python reader/writer.

- 
  Modify SCE2AnnData() to covert sparse matrices to dgRMatrix
  when they are. This mostly applies to assays and should be more
  compatible with what is expected by Python packages.

Minor changes

- 
  Add a testload argument to basiliskRun() calls which may help
  with problems creating Python environments

- 
  Updates to documentation and tests

Bug fixes

- 
  Improve handling of missing row or column names in
  SCE2AnnData()

[zenith](/packages/zenith)
------

                        Changes in version 1.8.1                        

- March 19, 2025
- fix get_MSigDB() to handle change to msigdbr
- now use msigdbdf


NEWS from existing Data Experiment Packages
===================================

[MetaScope](/packages/MetaScope)
---------

                        Changes in version 3.21                         

Major changes

- All MetaBLAST functionality and references removed.

[pRolocdata](/packages/pRolocdata)
----------

                       Changes in version 1.45.1                        

- Lisa Breckels is now maintainer.

[scRNAseq](/packages/scRNAseq)
--------

                       Changes in version 2.22.0                        

- Added the listPaths() function to list available paths to
  subdatasets. This is now referenced by an error message in
  fetchDataset() when a path= is required but not supplied.


[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.19.12                       

NEW FEATURES

- @lahuuki added guide_point_size as a argument to vis_clus() and
  vis_grid_clus(), which allows controlling the size of the points in
  the legends for the discrete variable plots. See
  https://github.com/LieberInstitute/spatialLIBD/pull/104 for more
  details.

                       Changes in version 1.19.11                       

BUG FIXES

- @Nick-Eagles fixed a bug for specifying colors in vis_clus() and
  related functions when NA values are present in the cluster
  variable. For more details check
  https://github.com/LieberInstitute/spatialLIBD/pull/102.

                       Changes in version 1.19.10                       

NEW FEATURES

- @lahuuki added the gene_name argument to sig_genes_extract() and
  sig_genes_extract_all() to make these functions more flexible. See
  https://github.com/LieberInstitute/spatialLIBD/pull/101 for details.

                       Changes in version 1.19.9                        

BUG FIXES

- @lahuuki resolved the bug we discovered on a dataset being analyzed
  by @Nick-Eagles as we documented at
  https://github.com/LieberInstitute/spatialLIBD/issues/98. Basically,
  an internal function used by layer_stat_cor_plot() was having issues
  with cluster names that were just numbers such as "1", "2", ...,
  "13" where annotate_registered_clusters() could result in a "9/13"
  annotation but then that would match with "1", "3", "9", and "13"
  instead of just "9" and "13", thus introducing incorrect "X"
  annotations on layer_stat_cor_plot(). @lahuuki resolved this at
  https://github.com/LieberInstitute/spatialLIBD/pull/99.

                       Changes in version 1.19.8                        

NEW FEATURES

- vis_gene() now has the cap_percentile argument as implemented by
  @Nick-Eagles. This allows you to cap the expression values at a
  certain percentiles, which can be useful to exclude super high
  outlier values which are common with Visium HD data. See
  https://github.com/LieberInstitute/spatialLIBD/pull/97 for details.

BUG FIXES

- Fixed colors on layer_stat_cor_plot() thanks to @lahuuki. Details at
  https://github.com/LieberInstitute/spatialLIBD/pull/94.

                       Changes in version 1.19.7                        

BUG FIXES

- Removed code that made the app very slow to load using HDF5Array
  objects. See
  https://github.com/LieberInstitute/spatialLIBD/commit/b80d92c3271a6ad92859f79a3bc343f77bad9bf2
  for details.

- Fixed a bug on the Z-score calculation in multi_gene_z_score(). See
  https://github.com/LieberInstitute/spatialLIBD/commit/2d17ea2c3d1b73b38fbd50503765052c8487b9b1
  for details. Implemented by @Nick-Eagles.

- Fix a bug in the documentation of run_app() for the example stitched
  data.

- No longer point to Twitter, instead point to Bluesky. This is for
  the package main README file as well as the default app
  documentation files.

                       Changes in version 1.19.5                        

BUG FIXES

- Fixed internal errors in add_qc_metrics() on the
  scuttle::isOutlier() function calls.

                       Changes in version 1.19.4                        

NEW FEATURES

- @lahuuki fully re-implemented gene_set_enrichment_plot() using
  ComplexHeatmap::Heatmap(). This new version has several new
  arguments that allow adding more annotation to the resulting
  heatmap. See https://github.com/LieberInstitute/spatialLIBD/pull/93
  for more details. This also means that layer_matrix_plot() has been
  removed from the package since it previously served as a helper
  function for gene_set_enrichment_plot().

                       Changes in version 1.19.3                        

BUG FIXES

- Resolved https://github.com/LieberInstitute/spatialLIBD/issues/90
  which made add_key() too strict and would create issues with
  export_cluster(). Reported by @lahuuki and @manishabarse.

                       Changes in version 1.19.2                        

BUG FIXES

- Merged https://github.com/LieberInstitute/spatialLIBD/pull/92 by
  @lahuuki. This fixes
  https://github.com/LieberInstitute/spatialLIBD/issues/72 and
  https://github.com/LieberInstitute/spatialLIBD/issues/48 by making
  registration_pseudobulk() more robust. The original issues were
  reported by @boyiguo1 and @berniejmulvey.

                       Changes in version 1.19.1                        

NEW FEATURES

- Merged https://github.com/LieberInstitute/spatialLIBD/pull/91 by
  @lahuuki. This pull request fully re-implemented
  layer_stat_cor_plot() with a version that uses
  ComplexHeatmap::Heatmap() internally. It also adds support for
  incorporating the automatic annotation results from
  annotate_registered_clusters(). NOTE that the max argument was
  renamed to color_max, as well as min to color_min. Also, the default
  for min used to be -max and now for color_min the default is the
  min() correlation observed. The default for max was 0.81 and the
  default for color_max() is the max() observed correlation.

- run_app() was also updated to match the updated in
  layer_stat_cor_plot() and now has 2 new inputs for controlling the
  annotation process with annotate_registered_clusters(). It also
  allows downloading a CSV file with the annotation results.

[STexampleData](/packages/STexampleData)
-------------

                 Changes in version 1.15.2 (2025-03-10)                 

- bug fixes: re-built objects Visium_humanDLPFC, Visium_mouseCoronal,
  Janesick_breastCancer_Chromium, Janesick_breastCancer_Visium

Deprecated and Defunct Packages
===============================

**SOFTWARE:**

Twenty software packages were removed from this release (after being deprecated
in Bioc 3.20):

- ATACCoGAPS, BiocOncoTK, biodbExpasy, biodbKegg, brainflowprobes, BRGenomics, CellaRepertorium, microbiomeMarker, MQmetrics, nanotatoR, netOmics, Pi, polyester, psygenet2r, RandomWalkRestartMH, rDGIdb, Risa, RNAinteract, single, SummarizedBenchmark
  
- Please note:  HTqPCR, previously announced as deprecated in 3.20, has been updated and remains in Bioconductor. 

Thirty nine software packages are deprecated in this release and will be removed in Bioc 3.22:

- AneuFinder, BEARscc, CBEA, chromstaR, coMET, crossmeta, dce, DeProViR, DIAlignR, Director, erma, GeneGeneInteR, genoCN, gespeR, girafe, GraphPAC, HTSeqGenie, iPAC, MAGeCKFlute, netDx, NeuCA, PanViz, pareg, paxtoolsr, PICS, PING, QuartPAC, RBioinf, ReactomeContentService4R, RGMQL, Rtreemix, SpacePAC, staRank, STdeconvolve, supraHex, synapter, trigger, TypeInfo, zlibbioc

**EXPERIMENT DATA:** 

Two experimental data packages were removed from this release (after being
deprecated in BioC 3.20):

-  DmelSGI, RNAinteractMAPK

Three experimental data packages are deprecated in this release and will be
removed in Bioc 3.22:

- benchmarkfdrData2019, parathyroidSE, synapterdata

**ANNOTATION DATA:** 

No annotation packages were removed from this release.

Three annotation packages are deprecated in this release and will be removed in
3.22:

- mirbase.db, targetscan.Hs.eg.db, targetscan.Mm.eg.db


**WORKFLOWS:** 

No workflow packages were removed from this release.

One workflow package was deprecated in this release:

- BiocMetaWorkflow

**BOOKS:**

No books were removed from this release.

No books were deprecated in this release.
