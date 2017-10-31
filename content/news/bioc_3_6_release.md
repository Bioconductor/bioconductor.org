October 31, 2017

Bioconductors:

We are pleased to announce Bioconductor 3.6, consisting of 1473
software packages, 326 experiment data packages, and 911 annotation
packages.

There are 100 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.6 is compatible with R 3.4.2,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release will include an updated Bioconductor [Amazon Machine Image][1]
and [Docker containers][2].

Visit [https://bioconductor.org][3]
for details and downloads.

[1]: https://bioconductor.org/help/bioconductor-cloud-ami/
[2]: https://bioconductor.org/help/docker/
[3]: https://bioconductor.org

Contents
--------

* [Getting Started with Bioconductor 3.6](#getting-started-with-bioconductor-35)
* [New Software Packages](#new-software-packages)
* [NEWS from new and existing packages](#news-from-new-and-existing-packages)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.6
======================================

To update to or install Bioconductor 3.6:

1. Install R >=3.4.2.  Bioconductor 3.6 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 100 new software packages in this release of Bioconductor.

- [amplican](https://bioconductor.org/packages/amplican) `amplican`
  performs alignment of the amplicon reads, normalizes gathered data,
  calculates multiple statistics (e.g. cut rates, frameshifts) and
  presents results in form of aggregated reports. Data and statistics
  can be broken down by experiments, barcodes, user defined groups,
  guides and amplicons allowing for quick identification of potential
  problems.

- [ANF](https://bioconductor.org/packages/ANF) This package is used
  for complex patient clustering by integrating multi-omic data
  through affinity network fusion.

- [anota2seq](https://bioconductor.org/packages/anota2seq) anota2seq
  provides analysis of translational efficiency and differential
  expression analysis for polysome-profiling and ribosome-profiling
  studies (two or more sample classes) quantified by RNA sequencing
  or DNA-microarray. Polysome-profiling and ribosome-profiling
  typically generate data for two RNA sources; translated mRNA and
  total mRNA. Analysis of differential expression is used to estimate
  changes within each RNA source (i.e. translated mRNA or total
  mRNA). Analysis of translational efficiency aims to identify
  changes in translation efficiency leading to altered protein levels
  that are independent of total mRNA levels (i.e. changes in
  translated mRNA that are independent of levels of total mRNA) or
  buffering, a mechanism regulating translational efficiency so that
  protein levels remain constant despite fluctuating total mRNA
  levels (i.e. changes in total mRNA that are independent of levels
  of translated mRNA). anota2seq applies analysis of partial variance
  and the random variance model to fulfill these tasks.

- [apeglm](https://bioconductor.org/packages/apeglm) apeglm provides
  Bayesian shrinkage estimators for effect sizes for a variety of GLM
  models, using approximation of the posterior for individual
  coefficients.

- [AUCell](https://bioconductor.org/packages/AUCell) AUCell allows to
  identify cells with active gene sets (e.g. signatures, gene
  modules...) in single-cell RNA-seq data. AUCell uses the "Area
  Under the Curve" (AUC) to calculate whether a critical subset of
  the input gene set is enriched within the expressed genes for each
  cell. The distribution of AUC scores across all the cells allows
  exploring the relative expression of the signature. Since the
  scoring method is ranking-based, AUCell is independent of the gene
  expression units and the normalization procedure. In addition,
  since the cells are evaluated individually, it can easily be
  applied to bigger datasets, subsetting the expression matrix if
  needed.

- [BASiCS](https://bioconductor.org/packages/BASiCS) Single-cell mRNA
  sequencing can uncover novel cell-to-cell heterogeneity in gene
  expression levels in seemingly homogeneous populations of cells.
  However, these experiments are prone to high levels of technical
  noise, creating new challenges for identifying genes that show
  genuine heterogeneous expression within the population of cells
  under study. BASiCS (Bayesian Analysis of Single-Cell Sequencing
  data) is an integrated Bayesian hierarchical model to perform
  statistical analyses of single-cell RNA sequencing datasets in the
  context of supervised experiments (where the groups of cells of
  interest are known a priori, e.g. experimental conditions or cell
  types). BASiCS performs built-in data normalisation (global
  scaling) and technical noise quantification (based on spike-in
  genes). BASiCS provides an intuitive detection criterion for highly
  (or lowly) variable genes within a single group of cells.
  Additionally, BASiCS can compare gene expression patterns between
  two or more pre-specified groups of cells. Unlike traditional
  differential expression tools, BASiCS quantifies changes in
  expression that lie beyond comparisons of means, also allowing the
  study of changes in cell-to-cell heterogeneity. The latter are
  quantified via a biological over-dispersion parameter that measures
  residual over-dispersion (with respect to Poisson sampling) after
  normalisation and technical noise removal.

- [beachmat](https://bioconductor.org/packages/beachmat) Provides a
  consistent C++ class interface for a variety of commonly used
  matrix types, including sparse and HDF5-backed matrices.

- [BiocSklearn](https://bioconductor.org/packages/BiocSklearn) This
  package provides interfaces to selected sklearn elements, and
  demonstrates fault tolerant use of python modules requiring
  extensive iteration.

- [bnbc](https://bioconductor.org/packages/bnbc) Tools to normalize
  (several) Hi-C data from replicates.

- [cbaf](https://bioconductor.org/packages/cbaf) This package
  contains functions that allow analysing and comparing various gene
  groups from different cancers/cancer subgroups easily. So far, it
  is compatible with RNA-seq, microRNA-seq, microarray and
  methylation datasets that are stored on cbioportal.org.

- [CEMiTool](https://bioconductor.org/packages/CEMiTool) The CEMiTool
  package unifies the discovery and the analysis of coexpression gene
  modules in a fully automatic manner, while providing a
  user-friendly html report with high quality graphs. Our tool
  evaluates if modules contain genes that are over-represented by
  specific pathways or that are altered in a specific sample group.
  Additionally, CEMiTool is able to integrate transcriptomic data
  with interactome information, identifying the potential hubs on
  each network.

- [ChIPanalyser](https://bioconductor.org/packages/ChIPanalyser)
  Based on a statistical thermodynamic framework, ChIPanalyser tries
  to produce ChIP-seq like profile. The model relies on four
  consideration: TF binding sites can be scored using a Position
  weight Matrix, DNA accessibility plays a role in Transcription
  Factor binding, binding profiles are dependant on the number of
  transcription factors bound to DNA and finally binding energy
  (another way of describing PWM's) or binding specificity should be
  modulated (hence the introduction of a binding specificity
  modulator). The end result of ChIPanalyser is to produce profiles
  simulating real ChIP-seq profile and provide accuracy measurements
  of these predicted profiles after being compared to real ChIP-seq
  data. The ultimate goal is to produce ChIP-seq like profiles
  predicting ChIP-seq like profile to circumvent the need to produce
  costly ChIP-seq experiments.

- [chromswitch](https://bioconductor.org/packages/chromswitch)
  Chromswitch implements a flexible method to detect chromatin state
  switches between samples in two biological conditions in a specific
  genomic region of interest given peaks or chromatin state calls
  from ChIP-seq data.

- [chromVAR](https://bioconductor.org/packages/chromVAR) Determine
  variation in chromatin accessibility across sets of annotations or
  peaks. Designed primarily for single-cell or sparse chromatin
  accessibility data, e.g. from scATAC-seq or sparse bulk ATAC or
  DNAse-seq experiments.

- [ClusterJudge](https://bioconductor.org/packages/ClusterJudge)
  ClusterJudge implements the functions, examples and other software
  published as an algorithm by Gibbons, FD and Roth FP. The article
  is called "Judging the Quality of Gene Expression-Based Clustering
  Methods Using Gene Annotation" and it appeared in Genome Research,
  vol. 12, pp1574-1581 (2002). See package?ClusterJudge for an
  overview.

- [coexnet](https://bioconductor.org/packages/coexnet) Extracts the
  gene expression matrix from GEO DataSets (.CEL files) as a
  AffyBatch object. Additionally, can make the normalization process
  using two different methods (vsn and rma). The summarization (pass
  from multi-probe to one gene) uses two different criteria (Maximum
  value and Median of the samples expression data) and the process of
  gene differentially expressed analisys using two methods (sam and
  acde). The construction of the co-expression network can be
  conduced using two different methods, Pearson Correlation
  Coefficient (PCC) or Mutual Information (MI) and choosing a
  threshold value using a graph theory approach.

- [consensusOV](https://bioconductor.org/packages/consensusOV) This
  package implements four major subtype classifiers for high-grade
  serous (HGS) ovarian cancer as described by Helland et al. (PLoS
  One, 2011), Bentink et al. (PLoS One, 2012), Verhaak et al. (J Clin
  Invest, 2013), and Konecny et al. (J Natl Cancer Inst, 2014). In
  addition, the package implements a consensus classifier, which
  consolidates and improves on the robustness of the proposed subtype
  classifiers, thereby providing reliable stratification of patients
  with HGS ovarian tumors of clearly defined subtype.

- [cytolib](https://bioconductor.org/packages/cytolib) This package
  provides the core data structure and API to represent and interact
  with the gated cytometry data.

- [DASC](https://bioconductor.org/packages/DASC) The package is used
  for identifying batches in high-dimensional dataset.

-
  [DelayedMatrixStats](https://bioconductor.org/packages/DelayedMatrixStats)
  A port of the 'matrixStats' API for use with DelayedMatrix objects
  from the 'DelayedArray' package. High-performing functions
  operating on rows and columns of DelayedMatrix objects, e.g. col /
  rowMedians(), col / rowRanks(), and col / rowSds(). Functions
  optimized per data type and for subsetted calculations such that
  both memory usage and processing time is minimized.

- [DEP](https://bioconductor.org/packages/DEP) This package provides
  an integrated analysis workflow for robust and reproducible
  analysis of mass spectrometry proteomics data for differential
  protein expression or differential enrichment. It requires tabular
  input (e.g. txt files) as generated by quantitative analysis
  softwares of raw mass spectrometry data, such as MaxQuant or
  IsobarQuant. Functions are provided for data preparation,
  filtering, variance normalization and imputation of missing values,
  as well as statistical testing of differentially enriched /
  expressed proteins. It also includes tools to check intermediate
  steps in the workflow, such as normalization and missing values
  imputation. Finally, visualization tools are provided to explore
  the results, including heatmap, volcano plot and barplot
  representations. For scientists with limited experience in R, the
  package also contains wrapper functions that entail the complete
  analysis workflow and generate a report. Even easier to use are the
  interactive Shiny apps that are provided by the package.

- [diffuStats](https://bioconductor.org/packages/diffuStats) Label
  propagation approaches are a widely used procedure in computational
  biology for giving context to molecular entities using network
  data. Node labels, which can derive from gene expression,
  genome-wide association studies, protein domains or metabolomics
  profiling, are propagated to their neighbours in the network,
  effectively smoothing the scores through prior annotated knowledge
  and prioritising novel candidates. The R package diffuStats
  contains a collection of diffusion kernels and scoring approaches
  that facilitates their computation and benchmarking.

- [DMCHMM](https://bioconductor.org/packages/DMCHMM) DMCHMM is a
  novel profiling tool for identifying differentially methylated CpG
  sites using Hidden Markov Model in bisulfite sequencing data.

- [Doscheda](https://bioconductor.org/packages/Doscheda) Doscheda
  focuses on quantitative chemoproteomics used to determine protein
  interaction profiles of small molecules from whole cell or tissue
  lysates using Mass Spectrometry data. The package provides a shiny
  application to run the pipeline, several visualisations and a
  downloadable report of an experiment.

- [EpiDISH](https://bioconductor.org/packages/EpiDISH) EpiDISH is a R
  package to infer the proportions of a priori known cell subtypes
  present in a sample representing a mixture of such cell-types.
  Inference proceeds via one of 3 methods (Robust Partial
  Correlations-RPC, Cibersort (CBS), Constrained Projection (CP)), as
  determined by user.

- [epivizrChart](https://bioconductor.org/packages/epivizrChart) This
  package provides an API for interactive visualization of genomic
  data using epiviz web components. Objects in R/BioConductor can be
  used to generate interactive R markdown/notebook documents or can
  be visualized in the R Studio's default viewer.

- [esATAC](https://bioconductor.org/packages/esATAC) This package
  provides a framework and complete preset pipeline for the
  quantification and analysis of ATAC-seq Reads. It covers a complete
  workflow starting from raw sequence reads, over creation of
  alignments and quality control report, to the quantification of
  genomic regions of interest. The package is managed by dataflow
  graph, and users can also build their own pipeline easily and
  flexibly.

- [GA4GHshiny](https://bioconductor.org/packages/GA4GHshiny)
  GA4GHshiny package provides an easy way to interact with data
  servers based on Global Alliance for Genomics and Health (GA4GH)
  genomics API through a Shiny application. It also integrates with
  Beacon Network.

- [GENIE3](https://bioconductor.org/packages/GENIE3) This package
  implements the GENIE3 algorithm for inferring gene regulatory
  networks from expression data.

- [HiCcompare](https://bioconductor.org/packages/HiCcompare)
  HiCcompare provides functions for joint normalization and
  difference detection in multiple Hi-C datasets. HiCcompare operates
  on processed Hi-C data in the form of chromosome-specific chromatin
  interaction matrices. It accepts three-column tab-separated text
  files storing chromatin interaction matrices in a sparse matrix
  format which are available from several sources. HiCcompare is
  designed to give the user the ability to perform a comparative
  analysis on the 3-Dimensional structure of the genomes of cells in
  different biological states.`HiCcompare` differs from other
  packages that attempt to compare Hi-C data in that it works on
  processed data in chromatin interaction matrix format instead of
  pre-processed sequencing data. In addition, `HiCcompare` provides a
  non-parametric method for the joint normalization and removal of
  biases between two Hi-C datasets for the purpose of comparative
  analysis. `HiCcompare` also provides a simple yet robust
  permutation method for detecting differences between Hi-C datasets.

- [InterMineR](https://bioconductor.org/packages/InterMineR)
  Databases based on the InterMine platform such as FlyMine, modMine
  (modENCODE), RatMine, YeastMine, HumanMine and TargetMine are
  integrated databases of genomic, expression and protein data for
  various organisms. Integrating data makes it possible to run
  sophisticated data mining queries that span domains of biological
  knowledge. This R package provides interfaces with these databases
  through webservices. It makes most from the correspondence of the
  data frame object in R and the table object in databases, while
  hiding the details of data exchange through XML or JSON.

- [IntramiRExploreR](https://bioconductor.org/packages/IntramiRExploreR)
  Intra-miR-ExploreR, an integrative miRNA target prediction
  bioinformatics tool, identifies targets combining expression and
  biophysical interactions of a given microRNA (miR). Using the tool,
  we have identified targets for 92 intragenic miRs in D.
  melanogaster, using available microarray expression data, from
  Affymetrix 1 and Affymetrix2 microarray array platforms, providing
  a global perspective of intragenic miR targets in Drosophila.
  Predicted targets are grouped according to biological functions
  using the DAVID Gene Ontology tool and are ranked based on a
  biologically relevant scoring system, enabling the user to identify
  functionally relevant targets for a given miR.

- [IrisSpatialFeatures](https://bioconductor.org/packages/IrisSpatialFeatures)
  IrisSpatialFeatures reads the output of the PerkinElmer inForm
  software and calculates a variety of spatial statistics. In
  addition to simple counts, it can derive average nearest neighbors
  for each cell-type and interaction summary profiles for each
  celltype. These statistics are derived across images, both overall
  and regions of interest as defined by user defined masks.

- [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR)
  IsoformSwitchAnalyzeR enabels identification and analysis of
  isoform switches with predicted functional consequences (such as
  gain/loss of protein domains etc) from quantification by Kallisto,
  Salmon, Cufflinks/Cuffdiff, RSEM etc.

- [iterClust](https://bioconductor.org/packages/iterClust) A
  framework for performing clustering analysis iteratively.

- [ivygapSE](https://bioconductor.org/packages/ivygapSE) Define a
  SummarizedExperiment and exploratory app for Ivy-GAP glioblastoma
  image, expression, and clinical data.

- [JASPAR2018](https://bioconductor.org/packages/JASPAR2018) Data
  package for JASPAR 2018. To search this databases, please use the
  package TFBSTools (>= 1.15.6).

- [M3C](https://bioconductor.org/packages/M3C) A central task in
  genomic data analyses for stratified medicine is class discovery
  which is accomplished through clustering. However, an unresolved
  problem with current clustering algorithms is they do not test the
  null hypothesis and derive p values. To solve this, we developed a
  novel hypothesis testing framework that uses consensus clustering
  called Monte Carlo Consensus Clustering (M3C). M3C use a multi-core
  enabled Monte Carlo simulation to generate a distribution of
  stability scores for each number of clusters using null datasets
  with the same gene-gene correlation structure as the real one.
  These distributions are used to derive p values and a beta
  distribution is fitted to the data to cheaply estimate p values
  beyond the limits of the simulation. M3C improves accuracy, allows
  rejection of the null hypothesis, removes systematic bias, and uses
  p values to make class number decisions. We believe M3C deals with
  a major pitfall in current automated class discovery tools.

- [MetaCyto](https://bioconductor.org/packages/MetaCyto) This package
  provides functions for preprocessing, automated gating and
  meta-analysis of cytometry data. It also provides functions that
  facilitate the collection of cytometry data from the ImmPort
  database.

- [methimpute](https://bioconductor.org/packages/methimpute) This
  package implements functions for calling methylated and
  unmethylated regions and estimate variability among a population of
  samples.

- [methInheritSim](https://bioconductor.org/packages/methInheritSim)
  Simulate a multigeneration methylation case versus control
  experiment with inheritance relation using a real control dataset.

- [methyvim](https://bioconductor.org/packages/methyvim) This package
  provides facilities for differential methylation analysis based on
  variable importance measures (VIMs), a class of statistical target
  parameters that arise in causal inference. The estimation and
  inference procedures provided are nonparametric, relying on
  ensemble machine learning to flexibly assess functional
  relationship among covariates and the outcome of interest. These
  tools can be applied to differential methylation at the level of
  CpG sites, with valid inference after multiple hypothesis testing.

- [mfa](https://bioconductor.org/packages/mfa) MFA models genomic
  bifurcations using a Bayesian hierarchical mixture of factor
  analysers.

- [microbiome](https://bioconductor.org/packages/microbiome)
  Utilities for microbiome analysis.

- [MIRA](https://bioconductor.org/packages/MIRA) MIRA measures the
  degree of "dip" in methylation level surrounding a regulatory site
  of interest, such as a transcription factor binding site, for
  instances of that type of site across the genome which can then be
  used to infer regulatory activity.

- [miRBaseConverter](https://bioconductor.org/packages/miRBaseConverter)
  A comprehensive tool for converting and retrieving the miRNA Name,
  Accession, Sequence, Version, History and Family information in
  different miRBase versions. It can process a huge number of miRNAs
  in a short time without other depends.

- [miRmine](https://bioconductor.org/packages/miRmine) miRmine
  database is a collection of expression profiles from different
  publicly available miRNA-seq datasets, Panwar et al (2017) miRmine:
  A Database of Human miRNA Expression, prepared with this data
  package as RangedSummarizedExperiment.

- [miRsponge](https://bioconductor.org/packages/miRsponge) This
  package provides several functions to study miRNA sponge (also
  called ceRNA or miRNA decoy), including popular methods for
  identifying miRNA sponge interactions, and the integrative method
  to integrate miRNA sponge interactions from different methods, as
  well as the functions to validate miRNA sponge interactions, and
  infer miRNA sponge modules, conduct enrichment analysis of modules,
  and conduct survival analysis of modules.

- [motifmatchr](https://bioconductor.org/packages/motifmatchr)
  Quickly find motif matches for many motifs and many sequences.
  Wraps C++ code from the MOODS motif calling library, which was
  developed by Pasi Rastas, Janne Korhonen, and Petri Martinmäki.

- [mpra](https://bioconductor.org/packages/mpra) Tools for data
  management, count preprocessing, and differential analysis in
  massively parallel report assays (MPRA).

- [MSstatsQC](https://bioconductor.org/packages/MSstatsQC) MSstatsQC
  is an R package which provides longitudinal system suitability
  monitoring tools for proteomic experiments.

- [multiMiR](https://bioconductor.org/packages/multiMiR) A collection
  of microRNAs/targets from external resources, including validated
  microRNA-target databases (miRecords, miRTarBase and TarBase),
  predicted microRNA-target databases (DIANA-microT, ElMMo,
  MicroCosm, miRanda, miRDB, PicTar, PITA and TargetScan) and
  microRNA-disease/drug databases (miR2Disease, Pharmaco-miR VerSe
  and PhenomiR).

- [ndexr](https://bioconductor.org/packages/ndexr) This package
  offers an interface to NDEx servers, e.g. the public server at
  http://ndexbio.org/. It can retrieve and save networks via the API.
  Networks are offered as RCX object and as igraph representation.

- [Onassis](https://bioconductor.org/packages/Onassis) A package that
  allows the annotation of text with ontology terms (mainly from OBO
  ontologies) and the computation of semantic similarity measures
  based on the structure of the ontology between different annotated
  samples.

- [oncomix](https://bioconductor.org/packages/oncomix) This package
  helps identify mRNAs that are overexpressed in subsets of tumors
  relative to normal tissue. Ideal inputs would be paired
  tumor-normal data from the same tissue from many patients (>15
  pairs). This unsupervised approach relies on the observation that
  oncogenes are characteristically overexpressed in only a subset of
  tumors in the population, and may help identify oncogene candidates
  purely based on differences in mRNA expression between previously
  unknown subtypes.

- [oneSENSE](https://bioconductor.org/packages/oneSENSE) A graphical
  user interface that facilitates the dimensional reduction method
  based on the t-distributed stochastic neighbor embedding (t-SNE)
  algorithm, for categorical analysis of mass cytometry data. With
  One-SENSE, measured parameters are grouped into predefined
  categories, and cells are projected onto a space composed of one
  dimension for each category. Each dimension is informative and can
  be annotated through the use of heatplots aligned in parallel to
  each axis, allowing for simultaneous visualization of two
  catergories across a two-dimensional plot. The cellular occupancy
  of the resulting plots alllows for direct assessment of the
  relationships between the categories.

- [ontoProc](https://bioconductor.org/packages/ontoProc) Support
  harvesting of diverse bioinformatic ontologies, making particular
  use of the ontologyIndex package on CRAN. We provide snapshots of
  key ontologies for terms about cells, cell lines, chemical
  compounds, and anatomy, to help analyze genome-scale experiments,
  particularly cell x compound screens.  Another purpose is to
  strengthen development of compelling use cases for richer
  interfaces to emerging ontologies.

- [openPrimeR](https://bioconductor.org/packages/openPrimeR) An
  implementation of methods for designing, evaluating, and comparing
  primer sets for multiplex PCR. Primers are designed by solving a
  set cover problem such that the number of covered template
  sequences is maximized with the smallest possible set of primers.
  To guarantee that high-quality primers are generated, only primers
  fulfilling constraints on their physicochemical properties are
  selected. A Shiny app providing a user interface for the
  functionalities of this package is provided by the 'openPrimeRui'
  package.

- [openPrimeRui](https://bioconductor.org/packages/openPrimeRui) A
  Shiny application providing methods for designing, evaluating, and
  comparing primer sets for multiplex polymerase chain reaction.
  Primers are designed by solving a set cover problem such that the
  number of covered template sequences is maximized with the smallest
  possible set of primers. To guarantee that high-quality primers are
  generated, only primers fulfilling constraints on their
  physicochemical properties are selected.

- [OPWeight](https://bioconductor.org/packages/OPWeight) This package
  perform weighted-pvalue based multiple hypothesis test and provides
  corresponding information such as ranking probability, weight,
  significant tests, etc . To conduct this testing procedure, the
  testing method apply a probabilistic relationship between the test
  rank and the corresponding test effect size.

- [panelcn.mops](https://bioconductor.org/packages/panelcn.mops) CNV
  detection tool for targeted NGS panel data. Extension of the
  cn.mops package.

- [PathwaySplice](https://bioconductor.org/packages/PathwaySplice)
  Pathway analysis of alternative splicing would be biased without
  accounting for the different number of exons associated with each
  gene, because genes with higher number of exons are more likely to
  be included in the 'significant' gene list in alternative splicing.
  PathwaySplice is an R package that: (1) performs pathway analysis
  that explicitly adjusts for the number of exons associated with
  each gene (2) visualizes selection bias due to different number of
  exons for each gene (3) formally tests for presence of bias using
  logistic regression (4) supports gene sets based on the Gene
  Ontology terms, as well as more broadly defined gene sets (e.g.
  MSigDB) or user defined gene sets (5) identifies the significant
  genes driving pathway significance (6) organizes significant
  pathways with an enrichment map, where pathways with large number
  of overlapping genes are grouped together in a network graph

- [phenopath](https://bioconductor.org/packages/phenopath) PhenoPath
  infers genomic trajectories (pseudotimes) in the presence of
  heterogeneous genetic and environmental backgrounds and tests for
  interactions between them.

- [progeny](https://bioconductor.org/packages/progeny) This package
  provides a function to infer pathway activity from gene expression
  using PROGENy. It contains the linear model we inferred in the
  publication "Perturbation-response genes reveal signaling
  footprints in cancer gene expression".

- [PROPS](https://bioconductor.org/packages/PROPS) This package
  calculates probabilistic pathway scores using gene expression data.
  Gene expression values are aggregated into pathway-based scores
  using Bayesian network representations of biological pathways.

- [Rbowtie2](https://bioconductor.org/packages/Rbowtie2) This package
  provides an R wrapper of the popular bowtie2 sequencing reads
  aligner and AdapterRemoval, a convenient tool for rapid adapter
  trimming, identification, and read merging.

- [restfulSE](https://bioconductor.org/packages/restfulSE) This
  package provides functions and classes to interface with remote
  data stores by operating on SummarizedExperiment-like objects.

- [rexposome](https://bioconductor.org/packages/rexposome) Package
  that allows to explore the exposome and to perform association
  analyses between exposures and health outcomes.

- [rhdf5client](https://bioconductor.org/packages/rhdf5client)
  Provides functionality for reading data from h5serv server from
  within R.

- [Rhdf5lib](https://bioconductor.org/packages/Rhdf5lib) Provides C
  and C++ hdf5 libraries.

- [RProtoBufLib](https://bioconductor.org/packages/RProtoBufLib) This
  package provides the headers and static library of Protocol buffers
  2.6.0 for other R packages to compile and link against.

- [RTNsurvival](https://bioconductor.org/packages/RTNsurvival)
  RTNsurvival is a tool for integrating regulons generated by the RTN
  package with survival information. For a given regulon, the
  2-tailed GSEA approach computes a differential Enrichment Score
  (dES) for each individual sample, and the dES distribution of all
  samples is then used to assess the survival statistics for the
  cohort. There are two main survival analysis workflows: a Cox
  Proportional Hazards approach used to model regulons as predictors
  of survival time, and a Kaplan-Meier analysis assessing the
  stratification of a cohort based on the regulon activity. All plots
  can be fine-tuned to the user's specifications.

- [runibic](https://bioconductor.org/packages/runibic) This package
  implements UbiBic algorithm in R. This biclustering algorithm for
  analysis of gene expression data was introduced by Zhenjia Wang et
  al. in 2016. It is currently considered the most promising
  biclustering method for identification of meaningful structures in
  complex and noisy data.

- [RVS](https://bioconductor.org/packages/RVS) Rare Variant Sharing
  (RVS) implements tests of association and linkage between rare
  genetic variant genotypes and a dichotomous phenotype, e.g. a
  disease status, in family samples. The tests are based on
  probabilities of rare variant sharing by relatives under the null
  hypothesis of absence of linkage and association between the rare
  variants and the phenotype and apply to single variants or multiple
  variants in a region (e.g. gene-based test).

- [Scale4C](https://bioconductor.org/packages/Scale4C) Scale4C is an
  R/Bioconductor package for scale-space transformation and
  visualization of 4C-seq data. The scale-space transformation is a
  multi-scale visualization technique to transform a 2D signal (e.g.
  4C-seq reads on a genomic interval of choice) into a tesselation in
  the scale space (2D, genomic position x scale factor) by applying
  different smoothing kernels (Gauss, with increasing sigma). This
  transformation allows for explorative analysis and comparisons of
  the data's structure with other samples.

- [scfind](https://bioconductor.org/packages/scfind) Recently a very
  large collection of single-cell RNA-seq (scRNA-seq) datasets has
  been generated and publicly released. For the collection to be
  useful, the information must be organized in a way that supports
  queries that are relevant to researchers. `scfind` builds an index
  from scRNA-seq datasets which organizes the information in a
  suitable and compact manner so that the datasets can be very
  efficiently searched for either cells or cell types in which a
  given list of genes is expressed.

- [scmap](https://bioconductor.org/packages/scmap) Single-cell
  RNA-seq (scRNA-seq) is widely used to investigate the composition
  of complex tissues since the technology allows researchers to
  define cell-types using unsupervised clustering of the
  transcriptome. However, due to differences in experimental methods
  and computational analyses, it is often challenging to directly
  compare the cells identified in two different experiments. scmap is
  a method for projecting cells from a scRNA-seq experiment on to the
  cell-types identified in a different experiment.

- [SCnorm](https://bioconductor.org/packages/SCnorm) This package
  implements SCnorm — a method to normalize single-cell RNA-seq data.

- [scoreInvHap](https://bioconductor.org/packages/scoreInvHap)
  scoreInvHap can get the samples' inversion status of known
  inversions. scoreInvHap uses SNP data as input and requires the
  following information about the inversion: genotype frequencies in
  the different haplotypes, R2 between the region SNPs and inversion
  status and heterozygote genotypes in the reference. The package
  include this data for two well known inversions (8p23 and 17q21.31)
  and for two additional regions.

- [scPipe](https://bioconductor.org/packages/scPipe) to process
  single cell RNA-seq data from fastq to gene counting matrix. it can
  process data generated by CEL-seq, MARS-seq, Drop-seq and
  SMART-seq.

- [seqCAT](https://bioconductor.org/packages/seqCAT) The seqCAT
  package uses variant calling data (in the form of VCF files) from
  high throughput sequencing technologies to authenticate and
  validate the source, function and characteristics of biological
  samples used in scientific endeavours.

- [seqcombo](https://bioconductor.org/packages/seqcombo) Provides
  useful functions for visualizing sequence recombination and virus
  reassortment events.

- [SeqSQC](https://bioconductor.org/packages/SeqSQC) The SeqSQC is
  designed to identify problematic samples in NGS data, including
  samples with gender mismatch, contamination, cryptic relatedness,
  and population outlier.

- [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment)
  Defines a S4 class for storing data from single-cell experiments.
  This includes specialized methods to store and retrieve spike-in
  information, dimensionality reduction coordinates and size factors
  for each cell, along with the usual metadata for genes and
  libraries.

- [slalom](https://bioconductor.org/packages/slalom) slalom is a
  scalable modelling framework for single-cell RNA-seq data that uses
  gene set annotations to dissect single-cell transcriptome
  heterogeneity, thereby allowing to identify biological drivers of
  cell-to-cell variability and model confounding factors.

- [SPONGE](https://bioconductor.org/packages/SPONGE) This package
  provides methods to efficiently detect competitive endogeneous RNA
  interactions between two genes. Such interactions are mediated by
  one or several miRNAs such that both gene and miRNA expression data
  for a larger number of samples is needed as input.

- [stageR](https://bioconductor.org/packages/stageR) The stageR
  package allows automated stage-wise analysis of high-throughput
  gene expression data. The method is currently published as a
  preprint at http://biorxiv.org/content/early/2017/02/16/109082

- [tenXplore](https://bioconductor.org/packages/tenXplore) Perform
  ontological exploration of scRNA-seq of 1.3 million mouse neurons
  from 10x genomics.

- [TFARM](https://bioconductor.org/packages/TFARM) It searches for
  relevant associations of transcription factors with a transcription
  factor target, in specific genomic regions. It also allows to
  evaluate the Importance Index distribution of transcription factors
  (and combinations of transcription factors) in association rules.

- [TFHAZ](https://bioconductor.org/packages/TFHAZ) It finds
  trascription factor (TF) high accumulation DNA zones, i.e., regions
  along the genome where there is a high presence of different
  transcription factors. Starting from a dataset containing the
  genomic positions of TF binding regions, for each base of the
  selected chromosome the accumulation of TFs is computed. Three
  different types of accumulation (TF, region and base accumulation)
  are available, together with the possibility of considering, in the
  single base accumulation computing, the TFs present not only in
  that single base, but also in its neighborhood, within a window of
  a given width. Then, high accumulation DNA zones are defined as
  regions of the genome formed by contiguous bases with accumulation
  equal or higher than a given threshold. In addition, some functions
  are provided in order to analyze, visualize and compare results
  obtained with different input parameters.

- [TMixClust](https://bioconductor.org/packages/TMixClust)
  Implementation of a clustering method for time series gene
  expression data based on mixed-effects models with Gaussian
  variables and non-parametric cubic splines estimation. The method
  can robustly account for the high levels of noise present in
  typical gene expression time series datasets.

- [TnT](https://bioconductor.org/packages/TnT) A R interface to the
  TnT javascript library (https://github.com/ tntvis) to provide
  interactive and flexible visualization of track-based genomic data.

- [topdownr](https://bioconductor.org/packages/topdownr) The topdownr
  package allows automatic and systemic investigation of fragment
  conditions. It creates Thermo Orbitrap Fusion Lumos method files to
  test hundreds of fragmentation conditions. Additionally it provides
  functions to analyse and process the generated MS data and
  determine the best conditions to maximise overall fragment
  coverage.

- [transcriptogramer](https://bioconductor.org/packages/transcriptogramer)
  R package for transcriptional analysis based on transcriptograms, a
  method to analyze transcriptomes that projects expression values on
  a set of ordered proteins, arranged such that the probability that
  gene products participate in the same metabolic pathway
  exponentially decreases with the increase of the distance between
  two proteins of the ordering. Transcriptograms are, hence, genome
  wide gene expression profiles that provide a global view for the
  cellular metabolism, while indicating gene sets whose expression
  are altered.

- [trena](https://bioconductor.org/packages/trena) Methods for
  reconstructing transcriptional regulatory networks, especially in
  species for which genome-wide TF binding site information is
  available.

- [vulcan](https://bioconductor.org/packages/vulcan) Vulcan (VirtUaL
  ChIP-Seq Analysis through Networks) is a package that interrogates
  gene regulatory networks to infer cofactors significantly enriched
  in a differential binding signature coming from ChIP-Seq data. In
  order to do so, our package combines strategies from different
  BioConductor packages: DESeq for data normalization, ChIPpeakAnno
  and DiffBind for annotation and definition of ChIP-Seq genomic
  peaks, csaw to define optimal peak width and viper for applying a
  regulatory network over a differential binding signature.

- [zFPKM](https://bioconductor.org/packages/zFPKM) Perform the zFPKM
  transform on RNA-seq FPKM data. This algorithm is based on the
  publication by Hart et al., 2013 (Pubmed ID 24215113). Reference
  recommends using zFPKM > -3 to select expressed genes. Validated
  with encode open/closed chromosome data. Works well for gene level
  data using FPKM or TPM. Does not appear to calibrate well for
  transcript level data.

- [zinbwave](https://bioconductor.org/packages/zinbwave) Implements a
  general and flexible zero-inflated negative binomial model that can
  be used to provide a low-dimensional representations of single-cell
  RNA-seq data. The model accounts for zero inflation (dropouts),
  over-dispersion, and the count nature of the data. The model also
  accounts for the difference in library sizes and optionally for
  batch effects and/or other covariates, avoiding the need for
  pre-normalize the data.

NEWS from new and existing packages
===================================

[ABAEnrichment](https://bioconductor.org/packages/ABAEnrichment)
-------------

Changes in version 1.8.0:

NEW FEATURES

- two new statistical tests are available (in addition to
  hypergeometric and Wilcoxon rank-sum test): 1) binomial test when
  genes are associated with two counts (e.g. number of amino-acid
  changes in two species since a common ancestor); 2) 2x2 contingency
  table test when genes are associated with with four counts (e.g.
  number of non-synonymous or synonymous variants that are fixed
  between or variable within species.)

- new 'silent' option represses all output to the screen except for
  warnings and error messages

USER-LEVEL CHANGES

- default 'genes' input is now a dataframe with gene-identifiers in
  column 1 and gene-associated variables, like scores or counts in the
  other columns (the old vector-input will still work for
  hypergeometric and Wilcxon rank-sum test)

- sort aba_enrich(...)[[1]] output also on structure-ID (after sorting
  on age category and FWER)

- sort genes in 'get_annotated_genes' output (mixedsort)

- rename aba_enrich() column 'times_FWER_under_0.05' to 'n_significant'

IMPROVEMENTS

- performance improvements for aba_enrich and get_annotated_genes

- modified generation of temporary files allows parallel processing
  with mclapply

[ABSSeq](https://bioconductor.org/packages/ABSSeq)
------

Version: 2017-07-18
Date: 2017-07-18
Category: Updating! Adapt aFold model for complex design with linear
        model (R funtion: ABSSeqlm
Text:

Version: 2017-06-02
Date: 2017-06-02
Text: 1) Updating! Introduce a new approach for normalization as qtotal
        (R funtion: qtotalNormalized) 2) aFold now moderates fold
        changes upon overall dispersion and gene-specific dispersion

Version: 2017-02-14
Date: 2017-02-14
Category: Updating! Adding new method: aFold to call DE via fold-change
Text:

Version: 2016-02-29
Date: 2016-02-29
Category: Updating! Enable the paired comparison by seting up 'paired'
        parameter in 'ABSDataSet' object
Text:

Version: 2015-08-24
Date: 2015-08-24
Category: This is the second version of ABSSeq
Text: 1) Use NB distribution insead of GP to model the counts
        difference between conditions 2) Introduce an efficient outlier
        detection method: moderated MAD 3) Introduce penalty for
        dispersion estimation 4) Provide moderation of log2 fold-change
        at expression level and gene-specific dispersion

[affycoretools](https://bioconductor.org/packages/affycoretools)
-------------

Version: 1.50.0
Category: Significant user-visible changes
Text:

Version: 1.50.0
Category: Converted runRomer and outputRomer to S4, in the process
        re-factoring
Text:

Version: 1.50.0
Category: they now accept either an ExpressionSet and MArrayLM object
        (for
Text:

Version: 1.50.0
Category: Affymetrix microarray analyses), a DGEList and DGEGLM object
        (for
Text:

Version: 1.50.0
Category: RNA-Seq analysis using edgeR) or an EList and MArrayLM object
        (for
Text:

Version: 1.50.0
Category: either Agilent one-color microarray analysis or limma-voom
        RNA-Seq
Text:

Version: 1.50.0
Category: analysis
Text:

Version: 1.50.0
Category: affycoretools v1.43.1
Text:

Version: 1.50.0
Category: Changes
Text:

Version: 1.50.0
Category: Added new functionality to automatically annotate
        ExpressionSets
Text:

Version: 1.50.0
Category: using either ChipDb (e.g., hugene10sttranscriptcluster.db)
        packages
Text:

Version: 1.50.0
Category: pdInfo (e.g., pd.hugene.1.0.v1) packages, or user-supplied
Text:

Version: 1.50.0
Category: data.frames.
Text:

[ALDEx2](https://bioconductor.org/packages/ALDEx2)
------

Changes in version 1.9.2:

NEW FEATURES (JRW, GBG

- added ability to auto-choose low variance high relative abundance
  features as the basis

- added new class and generic definitions to get the features used as
  basis

- the getDenom function can be exported

- updated documentation for above

- version bump

Changes in version 1.8.1:

NEW FEATURES (TPQ

- rennovated aldex function

- new 'test = iterative' uses results of one t-test as 'denom' input to
  a second

- large improvements to function documentation

- rennovated aldex.ttest function

- "progress bar" tracks progress across Monte-Carlo instances

- made aldex.ttest function faster (~300% for 10,000 features)

- now uses Wilcox signed rank test when 'paired.test = TRUE'

- added aldex.clr method for signature 'matrix'

[amplican](https://bioconductor.org/packages/amplican)
--------

Changes in version 1.0.0:

NEW FEATURES

- ampliCan first release

[Anaquin](https://bioconductor.org/packages/Anaquin)
-------

Changes in version 2.0:

- Introduced Major changes to all R functions

- Prepared upcoming sequins for germline and cancer mutations

- Prepared upcoming sequins for metagenomic sequins

[AneuFinder](https://bioconductor.org/packages/AneuFinder)
----------

Version: 1.5.1
Category: SIGNIFICANT USER-LEVEL CHANGES
Text:

Version: 1.5.1
Category: BUG FIXES
Text: Can use BSgenome*NCBI* for GC correction now (worked only for
        BSgenome*UCSC* before).

Version: 1.5.1
Category: BUG FIXES
Text: Bugfix in bam2GRanges(..., pairedEndReads=TRUE) if both
        alignments are not on the same chromosome.

[AnnotationFilter](https://bioconductor.org/packages/AnnotationFilter)
----------------

Changes in version 1.1.2:

NEW FEATURES

- supportFilters returns a data.frame with filter class name and field.

[AnnotationHub](https://bioconductor.org/packages/AnnotationHub)
-------------

Changes in version 2.10.0:

NEW FEATURES

- AnnotationHub will now work offline utilizing argument 'localHub';
  will also use this option automatically if no internet connection is
  detected.

- Added new GDSResource class

- Added documentation for creating an AnnotationHub package

MODIFICATIONS

- Modified tags vector when passed to display to improve speed of
  display querying

- Moved readMetadataFromCsv from AnnotationHubData.

- Removed listResources and loadResource from AnnotationHub; not
  implemented and only valid in ExperimentHub

BUG FIXES

- Expose snapshot less than or equal to release date

- Force rebuild of index if index file corrput or out of date

[AnnotationHubData](https://bioconductor.org/packages/AnnotationHubData)
-----------------

Changes in version 1.8.0:

NEW FEATURES

- Instead of using dropbox or ftp to deliver contributed resources to
  Bioconductor Core, temporary access to Annotation-Contributor user on
  S3 is utilized.

MODIFICATIONS

- Modified readMetadataFromCsv; make RDataPath mandatory entry and if
  location_prefix is Bioconductor S3 bucket the Rdatapath must start
  with the package name

BUG FIXES

- Add garbage collection to fix twobit memory allocation error

- Fix files not deleting do to special characters in file names

- Import dbGetQuery from DBI

- Remove hard coded biocVersion in unit tests

[anota2seq](https://bioconductor.org/packages/anota2seq)
---------

Changes in version 0.99.2:

- renaming to adhere to Bioconductor naming style for classes,
  functions and data objects

- correct mistake in man page of anota2seq_data

- modify default values (no NULL default value for mandatory arguments)

- fix a bug with argument fileName in the plot functions

Changes in version 0.99.1:

- the latex2 function of BiocStyle is now deprecated in the devel
  version of the package. Replace it by latex()

- remove inst/doc and build/ folders

- update man pages

- change default value of assayNum in anota2seqDataSetFromSE to 1

- formatting of anota2seq show method (and minor text edits in the
  description of the object)

- fixed issues related to sortBy parameter in anota2seqSelSigGenes

- minor text edits in messages given to user

- fixed bug regarding threshold lines in anota2seqPlotFC

- fixed bug in anota2seqSelSigGenes when useRVM = FALSE and only one
  gene is to be selected

Changes in version 0.99.0:

- First submitted version

[apeglm](https://bioconductor.org/packages/apeglm)
------

Changes in version 1.0.0:

- apeglm is released on Bioconductor! Methods can be accessed from
  within DESeq2 by calling 'lfcShrink' with type="apeglm".

[arrayQualityMetrics](https://bioconductor.org/packages/arrayQualityMetrics)
-------------------

Changes in version 3.33.2:

USER VISIBLE CHANGES

- Removed reliance on SVGAnnotation as the BioCextra respository is
  being withdrawn.  Currently SVG plots are now not interactive.

[ASpli](https://bioconductor.org/packages/ASpli)
-----

Changes in version 1.3.12:

FEATURES

- Events named previously as 'as' are now names 'Undefined AS' for the
  sake of clarity.

- Enhancement to the vignettes.

BUG FIXES

- Added verbose and filterWithContrasted arguments to DUreport method.

Changes in version 1.3.11:

FEATURES

- Change the assigning of some 'ES*' events to 'as' events.

- Added verbose option of DUreportBinSplice method

- junctionDUReport methods changed to junctionDUreport to be consistent
  with other methods names.

BUG FIXES

- Many corrections in the vignette and man pages.

- ZNF410_E013.gr.pdf file in vignettes/images was changed to
  ZNF410_E013_gr.pdf to be correctly used in the latex code.

Changes in version 1.3.10:

DOCUMENTATION

- Vignette is updated to reflect changes in objects, methods and
  functions. Minor revisions are still required.

BUG FIXES

- PlotGenomicRegions now can make plots from unmerged bams correctly.

Changes in version 1.3.9:

BUG FIXES

- plotBins function now shows x-tick labels correctly.

Changes in version 1.3.8:

FEATURES

- PlotGenomicRegions do not requires an ASpliCounts object

BUG FIXES

- AsDiscover and DUreport can handle conditions with one samples.
  Results should be taken with care, but ASpli doesn't break.

Changes in version 1.3.7:

FEATURES

- Added functions to filter elements in ASpliCounts objects by read
  counts and read density.

Changes in version 1.3.6:

FEATURES

- geneMode is the default mode of offset estimation.

Changes in version 1.2.1:

BUG FIXES

- AsDiscover() no longer contains mixed row data for some IR and IR*
  bins in irPIR table.

- DUReport() calculates correctly junction ratio for junctions sharing
  their end in junction.counts table.

- DUReport() no longer fails when there are no junctions that share
  their start or end.

- reacCounts() requires an explicit targets argument, instead of taking
  it from the global environment (in which may not exist with that
  name).

[attract](https://bioconductor.org/packages/attract)
-------

Version: 1.29.1
Date: 2017-05-31
Category: Fix bug when using microarray data and the reactive gene sets
Text:

[AUCell](https://bioconductor.org/packages/AUCell)
------

Changes in version 0.99.6:

- Interface changes to conform to Bioconductor requirements. (Function
  names and parameters, usage of "GeneSet" from GSEABase...)

[BASiCS](https://bioconductor.org/packages/BASiCS)
------

Changes in version 1.1.0 (2017-08-16):

- Minor changes required to pass BiocCheck (used formatR and
  BiocChecks)

Changes in version 1.0.13 (2017-08-15):

- Minor changes to documentation

- Minor changes to vignette

Changes in version 1.0.12 (2017-08-14):

- 'SummarizedExperiment' class replaced by 'SingleCellExperiment'

- Updated vignette

Changes in version 1.0.11 (2017-08-06):

- BASiCS_VarianceDecomp : small bugs fixed (unused Data object)

- colnames fixed on cell-specific slots ChainSC and ChainRNA

- Rcpp::checkUserInterrupt() activated in C++ code for interruptions

- 'grDevices' added as dependency

- 'WithSpikes' argument removed from 'BASiCS_LoadChain' function

Changes in version 1.0.10 (2017-08-06):

- BASiCS_TestDE : changes in output messages

- Gene name added in plot method for BASiCS_Chain objects

Changes in version 1.0.9 (2017-08-04):

- newBASiCS_Data : added check to ensure SpikeInfo is a data.frame

- BASiCS_TestDE : default value for EpsilonM & EpsilonD has changed

- Vignette + examples updated accordingly

Changes in version 1.0.8 (2017-08-03):

- 'Data' arg removed from args in: 'HiddenHeaderDetectHVG_LVG' +
  'HiddenVarDecomp' + 'BASiCS_DetectHVG' + 'BASiCS_DetectLVG' +
  'BASiCS_VarThresholdSearchHVG' + 'BASiCS_VarThresholdSearchLVG'

- 'object' arg renamed as 'Chain' in 'BASiCS_DetectHVG' +
  'BASiCS_DetectLVG' + 'BASiCS_VarThresholdSearchHVG' +
  'BASiCS_VarThresholdSearchLVG'

- Examples in BASiCS_MCMC updated to use built-in chain

- Removal of 'smoothPlot' usage in plots for 'BASiCS_Summary' class
  (avoids issues with log)

- Volcano plots added in BASiCS_TestDE function

- Resolved issues with log scale in 'smoothScatter'

- Minor bugs fixed in BASiCS_TestDE

- 'EviThreshold' arg replaced by 'ProbThreshold' in 'BASiCS_HVG' +
  'BASiCS_LVG'

- Vignette updated accordingly

Changes in version 1.0.7 (2017-08-03):

- Multiple changes in BASiCS_TestDE: * Added examples of table viewing
  in doc * Added chain extracts to run examples * Results for
  differential dispersion to exclude differential mean genes * EFDR /
  EFNR control plot added * MA plots added

- Main hidden functions moved into individual files

- 'graphics' + 'KernSmooth' added as dependence to use 'smoothScatter'

- 'smooth' argument changed to 'SmoothPlot' in plot functions/methods

- Added reference :: to basic functions (e.g. boxplot, var, acf, etc)

- devtools::check() OK - 1 warning : BiocStyle in vignette

Changes in version 1.0.6 (2017-08-02):

- 'matrixStats' added within 'Imports' in DESCRIPTION (to use
  'colMedians' function)

- Minor changes function-by-function (revision of documentation + minor
  bugs) * BASiCS_TestDE : Improved offset plots

- Major changes function-by-function * BASiCS_TestDE : Ref/Test
  notation changed to Group1/2 + output table names changed (to
  Mean/Disp) + offset corrected objects as output + separated output
  tables Mean/Disp + post prob threshold moved to hidden function
  (avoids duplicated code)

Changes in version 1.0.5 (2017-08-01):

- Major changes function-by-function * BASiCS_TestDE : plotting
  functionality for offset added

- ACF plot added to BASiCS_Chain plots

Changes in version 1.0.4 (2017-08-01):

- Minor changes function-by-function (revision of documentation + minor
  bugs) * BASiCS_VarThresholdSearchHVG / LVG : EFDR criteria updated as
  in BASiCS_DetectHVG / LVG * BASiCS_TestDE : added examples

- Major changes function-by-function * BASiCS_D_TestDE : no longer in
  use message added

- Examples in BASiCS_MCMC updated accordingly

- Vignette updated accordingly

- Updated DESCRIPTION file (authors and description)

- Updated welcome message (with link to new wiki)

Changes in version 1.0.3 (2017-08-01):

- Minor changes function-by-function (revision of documentation + minor
  bugs) * BASiCS_DetectHVG / LVG: added vertical/horizontal lines in
  EFDR/EFND plot, subtasks moved to hidden functions (avoids repeated
  code)

- Major changes function-by-function * BASiCS_DetectLVG : evidence
  threhold to be defined by EFDR only * BASiCS_DetectHVG / LVG: fix on
  default prob threshold = 0.5 when optimal is < 0.5

- Examples in BASiCS_MCMC updated accordingly

- Vignette updated accordingly

Changes in version 1.0.2 (2017-08-01):

- Minor changes function-by-function (revision of documentation + minor
  bugs) * BASiCS_MCMC : 'MCMC_Output' replaced by 'Chain' in examples

- Major changes function-by-function * BASiCS_MCMC : chains stored a
  single .Rds files, with BASiCS_Chain format (rather than separate
  .txt files) * BASiCS_LoadChain : updated accordingly *
  BASiCS_DetectHVG : evidence threhold to be defined by EFDR only

- Vignette updated accordingly

Changes in version 1.0.1 (2017-07-31):

- Minor changes function-by-function (revision of documentation + minor
  bugs) * newBASiCS_Data : added ref to GB paper + default par value
  doc + Nils as author * BASiCS_Filter : added ref to GB paper +
  default par value doc * makeExampleBASiCS_Data : added ref to GB
  paper + default par value doc + Nils as author * BASiCS_MCMC : added
  ref to GB paper + Nils as author + fix print of system.time + updated
  examples

- Major changes function by function * makeExampleBASiCS_Data : `Case1`
  param replaced by `Example` (vignette + unit test updated
  accordingly) * BASiCS_MCMC : default value for PriorDelta set to
  'log-normal' * HiddenBASiCS_MCMC_Start : new starting values for
  delta (as in main repo)

Changes in version 1.0.0 (2017-07-29):

- Initial big merge with @nilseling changes (in preparation for
  Bioconductor submission)

- For reference, here is a summary of the changes implemented by
  @nilseling: * Creation of individual .R for each R function (but all
  'Hidden' functions are in one file) * BASiCS_Data class replaced by
  SummarizedExperiment class (Bioconductor) * Contact email updated to
  'cvallejos@@turing.ac.uk' * Methods for BASiCS_Data and BASiCS_D_Data
  objects removed (these classes no longer exist) * Adjust
  newBASiCS_Data function to create a SummarizedExperiment object *
  'cat' calls replaced by 'message' calls * Clean up of BASiCS_D_TestDE
  function (removed unused args) * Clean up of BASiCS_DetectHVG_LVG
  function (using SummarizedExperiment class) * Clean up of BASiCS_MCMC
  function (using SummarizedExperiment class) * Clean up of BASiCS_Sim
  function (using SummarizedExperiment class; extra arg 'muSpikes') *
  Clean up of BASiCS_DetectHVG_LVG function (using SummarizedExperiment
  class) * Clean up of BASiCS_VarianceDecomp function (using
  SummarizedExperiment class) * Clean up of
  BASiCS_VarThresholdSearchHVG function (using SummarizedExperiment
  class) * Clean up of BASiCS_DenoisedCounts function (using
  SummarizedExperiment class) * Clean up of BASiCS_DenoisedRates
  function (using SummarizedExperiment class) * Clean up of
  HiddenBASiCS_MCMC_Start function (using SummarizedExperiment class) *
  Clean up of makeExampleBASiCS_Data function (using
  SummarizedExperiment class) * Removed methods associated to
  BASiCS_D_MCMC class (no longer in use) * NAMESPACE + Rd files updated
  accordingly * DESCRIPTION file updated accordingly * Vignette file
  updated accordingly

Changes in version 0.99.14 (2017-10-19):

- Added checks to ensure input `SingleCellExperiment` object has all
  info

Changes in version 0.99.13 (2017-10-17):

- @neling email added to DESCRIPTION + author contributions

Changes in version 0.99.12 (2017-10-17):

- Update of NEWS with changes by @neling

- Small typo fix in class validity from @neling last commit

- Changes in `BASiCS_Chain` and `BASiCS_Summary` classess to contain a
  single `list` slot where all parameters are stored. This allows
  greater flexibility to incorporate ongoing developments of BASiCS but
  does not affect user interface.

- `BASiCS_DenoisedCounts`, `BASiCS_DenoisedRates`, `BASiCS_DetectHVG`,
  `BASiCS_DetectLVG`, `BASiCS_MCMC`, `BASiCS_TestDE`,
  `HiddenVarDecomp`, `newBASiCS_Chain` functions and all methods
  updated accordingly

- `ChainSC` and `ChainRNA` example objects updated accordingly

- Unit tests updated accordingly

- Wiki page updated accordingly

- `updateObject` generic method imported from `BiocGenerics`

- New generic methods definition moved to AllGenerics.R

- Change in class definition to store current packageVersion

- `BASiCS_LoadChain` function updated to deal with old class definition
  and .txt storage.

- Shorter welcome message.

- Sanity checks in `BASiCS_TestDE` moved to separate file.

- Updated code indentation in `BASiCS_TestDE`

- Warning in `BASiCS_DetectHVG` and `BASiCS_DetectLVG` when user
  provides posterior probability thresholds as input

- Unused argument removed from `HiddenEFDR` and `HiddenEFNR` functions

- Minor changes to pass R CMD check

- Check library passes R CMD BiocCheck

- Tidy up NAMESPACE

- Tidy up NEWS

- Merge Devel_class into master branch

Changes in version 0.99.10 (2017-10-17):

- Updated README with BioC installation instructions and @nturaga
  acknowledgment

Changes in version 0.99.9 (2017-10-16):

- Updated URL and BugReports in DESCRIPTION

Changes in version 0.99.8 (2017-10-10):

- version bump to trigger a new Bioconductor build

Changes in version 0.99.7:

- `inst` folder has been deleted

- Removal of `...` arg in `segments` - `plot` method for
  `BASiCS_Summary`

Changes in version 0.99.5 (2017-09-26):

- Change of default `PriorDelta` value in `BASiCS_MCMC`

- Unit tests modified accordingly

Changes in version 0.99.4 (2017-09-25):

- Format in `BASiCS_MCMC`

- Minor bug fixed in `BASiCS_DetectLVG` plots

- Added examples in `BASiCS_DenoisedCounts` + `BASiCS_DenoisedRates`

- Removal of `Example` argument in `makeExampleBASiCS_D_Data`

- Format in vignette

- R CMD CHECK passed

- R CMD BiocCheck passed

Changes in version 0.99.3 (2017-09-21):

- Extra unit test (parameter estimates for a given seed)

- Format + vectorisation in `HiddenThresholdSearchTestDE` function

- `HiddenProbDE` removed as no longer needed

- Additional unit test for differential test added

- Format for `BASiCS_TestDE` function

- Format for `HiddenBASiCS_MCMC_Start` function

- Format for `HiddenChecksBASiCS` function

- Format for `BASiCS_D_TestDE`

- Format for all S4 methods

- 2^ added to FC in `BASiCS_TestDE` TableDisp

- Format for data documentation

- Format for classes definition

- Format + vectorization for `BASiCS_Filter`

- Format for `BASiCS_LoadChain`

- Removal of second example in `makeExampleBASiCS_Data` (no longer
  required due to DataSC, DataRNA)

- Format for `makeExampleBASiCS_Data`

- Format for `BASiCS_Sim`

- Format + matrixStats usage in `BASiCS_VarianceDecomp`

- Unit test for HVG/LVG detection added

- Minor bug in `BASiCS_TestDE` resolved (related to colMeans2 use)

- Format + matrixStats usage in `HiddenVarDecomp`

- `message` replaced by `cat` in Methods.R

- Removal of no-longer-needed Hidden functions

- Format in `HiddenHeaderDetectHVG_LVG`

- Format in `HiddenPlot1DetectHVG_LVG` + `HiddenPlot2DetectHVG_LVG`

- Format + matrixStats usage in `BASiCS_DetectHVG` + `BASiCS_DetectLVG`

- Format in `BASiCS_VarThresholdSearchHVG`

- Format in `HiddenThresholdSearchDetectHVG_LVG`

- Calculations in `BASiCS_DenoisedRates` moved to C++ to increase speed

- Format in `BASiCS_DenoisedCounts`

- Edits in `HiddenBASiCS_MCMCcpp` and associated C++ functions *
  Dimensions for mu storage upto q0 only * Creation of generic storage
  values * Formating of long lines * Removal of no-longer-required
  functions and debug lines

- Edits in `HiddenBASiCS_MCMCcpp` and associated C++ functions *
  Optimization of MH updates

- Edits in `HiddenBASiCS_MCMCcpp` and associated C++ functions * Batch
  & no-batch versions merged in a single function * Avoids double
  transponse for chain storage

Changes in version 0.99.2:

- 'newBASiCS_Data': format + checks moved to `HiddenChecksBASiCS_Data`
  function

- `HiddenChecksBASiCS_Data`: use of colMeans2, rowMeans2

- Format in welcome message

- 'newBASiCS_Chain': format

Changes in version 0.99.0 (2017-08-28):

- Minor changes to complete Bioconductor submission checklist

- Version re-numbered according to Bioconductor guidelines

Changes in version 0.7.30 (2017-07-26):

- Typo fixed in 'BASiCS:::HiddenBASiCS_MCMC_Start' function (thanks to
  @bdulken for pointing out this issue)

Changes in version 0.7.29 (2017-07-24):

- 'BASiCS:::HiddenBASiCS_MCMC_Start' function modified regarding the
  definition of starting values for mu (to avoid adding +1 when it's
  not necessary)

- 'BASiCS:::HiddenBASiCS_MCMC_Start' modified to have better starting
  values for delta (this is based on eq 2 from Vallejos et al, 2016)

Changes in version 0.7.28 (2017-07-21):

- Fixed dependency to 'data.table' to prevent errors in
  'BASiCS_LoadChain' function(thanks to Yongchao Ge for pointing out
  this issue). This change affects the DESCRIPTION file as well as the
  'BASiCS_LoadChain' function.

Changes in version 0.7.27 (2017-05-09):

- Bug fixed in 'BASiCS_Filter' function to stop the code crashing when
  there are genes with zero counts across all cells.

Changes in version 0.7.26 (2017-02-09):

- Update of 'BASiCS_Filter' function to include a 'BatchInfo' argument
  (thanks to Dmitriy Zhukov for suggesting this)

- Vignette file updated accordingly.

- Update of 'newBASiCS_Data' function to prevent issues related to
  unused levels when the 'BatchInfo' argument is set as a 'factor'
  (thanks to Kevin Rue-Albrecht for pointing out this issue)

Changes in version 0.7.25 (2017-01-11):

- Fixed typo in vignette file (thanks to Dmitriy Zhukov for pointing
  this out)

Changes in version 0.7.24 (2016-11-23):

- Median changed to mean when calculating identifiability constrain for
  the no spike case

Changes in version 0.7.23 (2016-11-21):

- Cleaning-up of identifiability constrain options for the no spike
  case

Changes in version 0.7.22 (2016-11-16):

- Fixed bug for 'ConstrainType' = 5 (affects no spike case only)

Changes in version 0.7.21 (2016-11-15):

- 'ConstrainType' = 4, 5 added (exclude genes with capture in less than
  'ConstrainLimit'100% of the cells)

Changes in version 0.7.20 (2016-11-14):

- 'newBASiCS_Data' function modified to allow no spikes case

Changes in version 0.7.19 (2016-11-14):

- RefGene information to be stored as .txt file

Changes in version 0.7.18 (2016-11-08):

- Function 'BASiCS_LoadChain' added to simplify the loading of
  pre-computed chains

- Imports from 'testthat' library added to NAMESPACE file

- Small changes in 'BASiCS_D' classes to allow offset

Changes in version 0.7.17 (2016-11-07):

- Error fixed in the definition of RefGenes

Changes in version 0.7.16 (2016-11-07):

- ConstrainType = 3 (stochastic ref) re-introduced for no-spikes case

Changes in version 0.7.15 (2016-11-07):

- ConstrainType = 1 (untrimmed) re-introduced for no-spike case
  (illustration purposes only)

- RefGene removed from AR report

- Added console report of RefGene choice

Changes in version 0.7.14 (2016-11-03):

- Debug message removed from 'muUpdateNoSpikes' function.

Changes in version 0.7.13 (2016-11-02):

- Typo on the arguments of 'HiddenBASiCS_MCMCcppNoSpikes' function has
  been fixed ('Constrain' instead of 'ConstrainType')

- Validity check '(y(i) > 1e-3)' moved later on in
  'deltaUpdateNoSpikes' to avoid breaking the code when it is not
  really necessary!

Changes in version 0.7.12 (2016-11-01):

- Change in defaul value of 'PriorDelta' (v 0.7.9) has been reverted
  for reproducibility issues. Message added to the start of the
  function so that users are aware of making the change.

- PriorDelta = 'gamma' is now allowed for no spikes case

Changes in version 0.7.11 (2016-10-31):

- Comments added to clarify no spikes case in 'BASiCS_MCMC' function.

Changes in version 0.7.10 (2016-10-31):

- Minor changes in 'BASiCS_MCMC' to have more meaningful variable names

- Minor changes in 'BASiCS_MCMC' to make no-spikes case functional

Changes in version 0.7.9 (2016-10-31):

- C++ code cleaned regarding parameter updates (no spikes case)

- Updated email in description file and welcome message

- 'BASiCS_MCMC' modified so that "log-normal" is the default value for
  PriorDelta

Changes in version 0.7.8 (2016-08-15):

- ConstrainType = 4 implemented (no spikes case only)

Changes in version 0.7.7 (2016-08-12):

- ConstrainType = 3 implemented (no spikes case only)

Changes in version 0.7.6 (2016-08-11):

- As v 0.7.5 but with extra argument for constrain limit (no spikes
  case only)

Changes in version 0.7.5 (2016-08-10):

- Modified constrain to include 'expressed' genes only (no spikes case
  only)

Changes in version 0.7.4 (2016-08-09):

- Modified constrain to include 'expressed' genes only (no spikes case
  only)

Changes in version 0.7.3 (2016-08-08):

- Reference set to be the gene closest to the average (no spikes case
  only)

Changes in version 0.7.2 (2016-08-05):

- Same as v 0.7.2 but more efficient

Changes in version 0.7.1 (2016-08-05):

- Random reference in sequential updates for mu (no spikes case only)

Changes in version 0.7.0 (2016-08-02):

- Sequential updates for mu in MCMC (no spikes case only)

Changes in version 0.6.9 (2016-08-01):

- Dirichlet-based proposals for mu in MCMC (no spikes case only)

Changes in version 0.6.8 (2016-07-28):

- Sequential updates for mu in MCMC (no spikes case only)

Changes in version 0.6.7 (2016-07-28):

- Minor changes in the validity checks of 'BASiCS_D_Data' objects

Changes in version 0.6.6 (2016-07-27):

- First working version of no spikes case (constrained mu's)

Changes in version 0.6.5:

- Modified identifiability constrains for no spikes case.

Changes in version 0.6.4 (2016-06-01):

- Safety checks added to 'phiUpdateNoSpikes'

Changes in version 0.6.3 (2016-05-31):

- Small update in wellcome message.

Changes in version 0.6.2 (2016-05-31):

- Template for 'phiUpdateNoSpikes' has been added.

- Minor change in 'phiUpdate' for all 'HiddenBASiCS_MCMCcpp' functions
  (to avoid re-writing phi0)

- Faster adaptation steps for 'HiddenBASiCS_MCMCcppBatch'
  'HiddenBASiCS_MCMCcppNoSpikes'

- 'HiddenBASiCS_MCMCcppNoSpikes' completed for testing stage (not ready
  for users)

Changes in version 0.6.1 (2016-05-25):

- Argument 'WithSpikes' added to 'makeExampleBASiCS_Data', generating
  an example data with no spikes

- Modifications to validity checks of 'BASiCS_Data' class to allow no
  spikes case (SpikeInput = 1)

- 'show' methods for 'BASiCS_Data' class modified to deal with no spike
  case

- Checks for 'BASiCS_Data' object modified to require: either spikes or
  batches

- 'BASiCS_MCMC' modified for no spikes case (if/else statement only,
  not yet functional)

- 'HiddenBASiCS_MCMC_Start' modified for no spikes case

- Function 'HiddenBASiCS_MCMCcppNoSpikes' added. Update for 'phi' is
  pending.

- Small vignette file for no spikes case added

Changes in version 0.5.11 (2016-05-24):

- 'importMethodsFrom(scran, computeSumFactors)' added to NAMESPACE

Changes in version 0.5.10 (2016-05-20):

- Updated bug in 'HiddenBASiCS_MCMC_Start' function (thanks to Joanna
  Dreux for noticing this issue)

Changes in version 0.5.9 (2016-05-19):

- Dependency to 'scran' added

- 'HiddenBASiCS_MCMC_Start' modified to use 'scran' estimates as
  starting values. This allows faster convergence.

- Updated 'README.md' file to build and access vignette during
  installation

Changes in version 0.5.8 (2016-05-18):

- After profiling of c++ functions, minor edits to improve memory
  usage.  These do not affect user interface

- Small change in funtion 'makeExampleBASiCS_Data' to avoid 'NOTICE'
  message

Changes in version 0.5.7 (2016-04-27):

- Welcome message added

Changes in version 0.5.6 (2016-04-19):

- Fixed message in `BASiCS_D_Test` function. This does not affect
  functionality. Only the message that is printed in the console.

Changes in version 0.5.5 (2016-04-18):

- Small fix when checking the validity of a `BASiCS_D_Data` object.
  Thanks to Nils Eling.

- Fixed issue with documentation files.

Changes in version 0.5.4 (2016-04-15):

- Minor fixes to the documentation files

- Argument `GenesSelect` was added to function `BASiCS_D_TestDE` to
  allow user-defined lists of genes to be included in the comparison
  between cell types.

Changes in version 0.5.3 (2016-03-09):

- Slots `BatchInfoTest` and `BatchInfoRef` added to `BASiCS_D_Data`
  class.

- `newBASiCS_D_Data` function modified to incorporate `BatchInfoTest`
  and `BatchInfoRef` slots.

- `CombineBASiCS_Data` function modified to incorporate `BatchInfoTest`
  and `BatchInfoRef` slots.

- `show` method for `BASiCS_D_Data` class modified to incorporate
  `BatchInfoTest` and `BatchInfoRef` slots.

- `GeneNames` argument missing in `makeExampleBASiCS_D_Data` has been
  added.

- Slots `thetaTest` and `thetaRef` of `BASiCS_D_Chain` objects modified
  to accept matrices (necessary to allow multiple batches).

Changes in version 0.5.2 (2016-03-07):

- In function `BASiCS_D_TestDE`. Probability threshold set to 0.90 when
  EFDR fails to calibrate (simulations under the null)

- Fixed version number in description file.

Changes in version 0.5.1 (2016-03-07):

- In function `BASiCS_D_TestDE`. Probability threshold set to 0.95 when
  EFDR fails to calibrate (simulations under the null)

Changes in version 0.5.0 (2016-03-03):

- Extended package description to include comparisons between
  pre-specified populations of cells

- `BASiCS_DV_Data`, `BASiCS_DV_Chain` and `BASiCS_DV_Summary` classes
  replaced by `BASiCS_D_Data`, `BASiCS_D_Chain` and `BASiCS_D_Summary`
  classes, respectively.

- Slot 'offset' added to `BASiCS_D_Chain` and `BASiCS_D_Summary`
  classes.

- Slot 'probHPD' added to `BASiCS_D_Summary` class.

- Creation of `CombineBASiCS_Data` function to combine 2 independent
  `BASiCS_Data` objects into one `BASiCS_D_Data` object.

- Offset value added to the input of `show` method for `BASiCS_D_Chain`
  and `BASiCS_D_Summary` classes.

- Creation of `CombineBASiCS_Chain` function to combine 2 independent
  `BASiCS_Chain` objects into one `BASiCS_D_Chain` object.

- 'GeneNames' slot added to `BASiCS_Data` class

Changes in version 0.4.1 (2016-02-05):

- Fixed bug in vignette (usage of `newBASiCS_Data` function). Same bug
  has been fixed in documentation of class `BASiCS_Data`.

- Optional argument 'Start' added to `BASiCS_MCMC` function to allow
  running chains with a variety of user-defined starting point. In
  general, we do not advise to use this argument as the default option
  has been tuned to facilitate convergence.

- Updated documentation of 'BASiCS_DetectHVG' and
  'BASiCS_VarianceDecomp' function (removing 'GeneNames' argument)

- Updated documentation for 'BASiCS_DV_TestDE' function.

- Missing @description sections added to several .Rd files.

Changes in version 0.4.0 (2015-12-21):

- Preliminary functions for differential expression analyses (mean and
  over-dispersion) have been added.

Changes in version 0.3.5 (2015-11-27):

- Optional argument `PriorDelta` added to `BASiCS_MCMC` function, to
  allow optional use of `gamma` or `log-normal` priors for delta.

- In `BASiCS_MCMC` function, change of default value of `s2.mu` to 0.5.
  Allows better shrinkage in situations where a gene has zero counts
  across all cells.

- In `BASiCS_Data` class. Count matrices where a gene has zero counts
  across all cells are now allowed. However, a warning is returned in
  such case.  Only use such matrices when running differential
  expression results.

Changes in version 0.3.4 (2015-11-11):

- Fix of example code of 'newBASiCS_Data' function. Thanks to Simon
  Andrews for pointing out this issue.

Changes in version 0.3.3 (2015-11-05):

- Same method as before but improved performance of C++ code

Changes in version 0.3.2 (2015-10-23):

- Minor changes to the output of `BASiCS_MCMC` function (colnames of
  the elements related to the parameter $\theta$)

- Addition of extra optional parameter `ls.phi0` to `BASiCS_MCMC`
  function.  This is helpful on situations where the default value led
  to slow mixing to the chains related to the normalising constants
  $\phi_j$'s.

Changes in version 0.3.1 (2015-09-17):

- New slot GeneNames in 'BASiCS_Data' class

- Changes in the constructor `newBASiCS_Data` to allow easier
  construction of `BASiCS_Data` objects

- Minor changes to some functions' output to use the new GeneNames slot
  of 'BASiCS_Data' class

Changes in version 0.3.0 (2015-07-31):

- New slot BacthInfo in 'BASiCS_Data' class

- Batch-speficic technical variability parameters are allowed

- 'BASiCS_VarianceDecomp' modified to accommodate batch membership
  (including graphical output)

Changes in version 0.2.1 (2015-07-23):

- Argument 'GeneNames' has been added to functions
  'BASiCS_VarianceDecomp', 'BASiCS_DetectHVG', 'BASiCS_DetectLVG' so
  that users can specify gene labels or names that will be used for
  these functions's output.

Changes in version 0.2.0 (2015-06-24):

- New parametrization for mRNA content size factors $\phi_j$ that
  improves mixing of the MCMC algorithm (using adaptive Dirichlet
  proposals)

- Updated vignette

Changes in version 0.1.6 (2015-06-01):

- Replacement of parameter-specific 'displayChain's functions by a
  generic method displayChainBASiCS

- Replacement of parameter-specific 'displaySummary's functions by a
  generic method displaySummaryBASiCS

- In 'BASiCS_DetectHVG' and 'BASiCS_DetectLVG': constrain added to that
  prob >= 0.5 is required to detect HVG and LVG

Changes in version 0.1.5 (2015-05-27):

- Small typo in 'BASiCS_VarThresholdSearchHVG' and
  'BASiCS_VarThresholdSearchLVG' has been fixed.

Changes in version 0.1.4 (2015-05-21):

- 'BASiCS_Denoise' function added. It produces a table of normalised
  and denoised expression counts (after removing the effect of
  technical variation).

Changes in version 0.1.3 (2015-05-18):

- 'plotBASiCSDispVsExp': 'log="xy"' argument added

- 'BASiCS_Data': internal checks added (pre-filter of samples and
  transcripts)

- 'BASiCS_MCMC' c++ code: Rcpp::checkUserInterrupt() added for fast
  user-requested interruption of MCMC sampler

- S4 method 'plot' for 'BASiCS_Chain' signature: argument 'Column'
  replaced by 'Gene' and 'Cell'

- Replacement of 'plotBASiCSNormPhi' and 'plotBASiCSNormS' by S4 method
  'plot' ('BASiCS_Summary' signature)

- S4 method 'plot' for 'BASiCS_Summary' signature: does not include
  plots for all model parameters

- S4 method 'plot' for 'BASiCS_Summary' signature: does not include
  scatterplots of gene-specific and cell-specific model parameters
  (replacing 'plotBASiCSExpVsDisp')

- 'BASiCS_VarianceDecomp' has been made visible to users (former
  'HiddenVarDecomp' function) and output does now allow gene
  annotation.

- 'BASiCS_Filter' function added.

- Vignette has been updated according to the changes above.
  Instructions for installation were added.

- 'BASiCS_Filter' function added.

- Dependency to R (>= 3.1.0) has been added.

Changes in version 0.1.2 (2015-04-23):

- In 'BASiCS_Sim': change 'sum(phi)==n' by 'all.equal(sum(phi),n)'

Changes in version 0.1.1 (2015-04-21):

- Improved documentation

Changes in version 0.1 (2015-04-20):

- Original release

[beachmat](https://bioconductor.org/packages/beachmat)
--------

Changes in version 1.0.0:

- New package beachmat, which provides a C++ API for handling various R
  matrix types.

[BgeeDB](https://bioconductor.org/packages/BgeeDB)
------

Changes in version 2.3.2:

- Compatibility of the TopAnat part of the package with the new Bgee
  release 14 (www.bgee.org/bgee14/).

- Better management of potential issues with the ontology structure

Changes in version 2.3.1:

- Compatibility with the new Bgee release 14 (www.bgee.org/bgee14/)
  only for the processed data download part.

- The TopAnat part still temporarily uses data from Bgee release 13.2,
  until the webservice is ready.

- Retro compatibility with Bgee release 13.2 (www.bgee.org/bgee13/)

[Biobase](https://bioconductor.org/packages/Biobase)
-------

Changes in version 2.37:

BUG FIXES

- '$' completion on eSet and ExpressionSet works in RStudio.

[bioCancer](https://bioconductor.org/packages/bioCancer)
---------

Version: 1.6.00
Category: metamorphosis: bioCancer is a radiant.data extension
Text:

Version: 1.6.00
Category: reduce size of package by half 14 -> 7 mb
Text:

Version: 1.5.12
Category: improve Reactome_ui functions
Text:

Version: 1.5.11
Category: add switch button to ui_circomics
Text:

Version: 1.5.11
Category: improve circomics functions
Text:

Version: 1.5.09
Category: delete commented files and figures
Text:

Version: 1.5.09
Category: cleanup ui_radiant, /Rbis, /quant, /bioCancer
Text:

Version: 1.5.08
Category: switchButton
Text:

Version: 1.5.07
Category: data.row.names(row.names, rowsi, i)
Text:

Version: 1.5.07
Category: some row.names duplicated: 11,12,13,14 --> row.names NOT used
Text:

Version: 1.5.06
Category: dplyr::mutate_each() is deprecated
Text:

Version: 1.5.06
Category: dply::summarise_each() is deprecated
Text:

Version: 1.5.06
Category: replace BiocStyle by prettydoc
Text:

Version: 1.5.05
Category: Warning in formals(fun) : argument is not a function
Text:

Version: 1.5.05
Category: Warning in body(fun) : argument is not a function
Text:

Version: 1.5.04
Category: Fix setting plot size
Text:

Version: 1.5.03
Category: Change the color rang of Circular layout plot as in standards
Text:

Version: 1.5.02
Category: update Correlation Methods
Text:

Version: 1.5.01
Category: replace .libPath() by path.package('bioCancer') in portal.R
        file
Text:

[BiocCheck](https://bioconductor.org/packages/BiocCheck)
---------

Changes in version 1.14:

NEW FEATURES

- NOTE when maintainer subscription to bioc-devel mailing list cannot
  be checked (checking requires mailing list admin password).

BUG FIXES

- Use shell quotes to allow spaces in package paths

[BiocFileCache](https://bioconductor.org/packages/BiocFileCache)
-------------

Changes in version 1.1:

NEW FEATURES

- (v. 1.1.7) bfcmeta() and friends allows arbitrary metadata on records
  and ability to query bfcmeta data.

USER-VISIBLE CHANGES

- (v. 1.1.16) if httr::HEAD fails don't try httr::GET, just report not
  available

- (v. 1.1.7) bfcquery() syntax changed to grep(), rather than SQL
  'LIKE'.

- (v. 1.1.7) bfcquery() supports user-specified 'fields='.

- (v. 1.1.5) queryCount renamed to bfccount

- (v. 1.1.2) Add user specified extension to file in cache

- (v. 1.1.2) Use dbExecute/dbSendStatement rather than pasting query

- Files in cache retain basename and extension of original file

[BiocInstaller](https://bioconductor.org/packages/BiocInstaller)
-------------

Changes in version 1.28.0:

NEW FEATURES

- biocLite() supports full URLs, e.g., to archived Bioconductor
  packages.

- Support MRAN (Microsoft R) archives.

[BiocParallel](https://bioconductor.org/packages/BiocParallel)
------------

Changes in version 1.12:

BUG FIXES

- (v. 1.11.1) Change registered backend initialization to first
  invocation, rather than on load.

- (v 1.11.8) Ensure registry is initiailized before (public) use. Issue
  #65

NEW FEATURES

- (v. 1.11.2) bpiterate() gains a progress counter.

- (v. 1.11.5) ipclock(), etc: inter-process locks and counters

[BiocStyle](https://bioconductor.org/packages/BiocStyle)
---------

Changes in version 2.6.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Rename 'latex2', 'pdf_document2' and 'html_document2' to 'latex',
  'pdf_document' and 'html_document', respectively

- Rename argument to 'pdf_document' from 'use.unsrturl' to
  'use_unsrturl' for consistency with rmarkdown naming scheme

- Split figure/table captions on ". " rather than "." in order to
  minimize chances of hitting decimal point in numbers, etc.

BUG FIXES AND IMPROVEMENTS TO PDF STYLE

- Make knitr option 'fig.pos' functional in 'pdf_document'

- Improve formatting of TOC chapter entries
  (https://github.com/Bioconductor/BiocStyle/issues/35)

- Adjust width of non-floating figures according to the
  'fig.wide'/'fig.small' setting (reported by Aaron Lun)

- Robustify \texttt, e.g.  by enabling backticks and fixing the
  substitution of control spaces

- Update LaTeX code highlighting macros for pandoc output (reported by
  Vince Carey)

- Silence spurious warning in 'latex()' emitted when the source file is
  missing a final EOL

- Fix broken code in knitr hook adjusting page width

- Include package version field in R Markdown documents even if they
  are missing an abstract
  (https://github.com/Bioconductor/BiocStyle/issues/29)

- Correct way of enclosing non-floating figures in LaTeX 'adjustwidth'
  environment (https://github.com/Bioconductor/BiocStyle/issues/28)

BUG FIXES AND IMPROVEMENTS TO HTML STYLE

- Correct displaying of superscript text

- Allow modifying table layout in 'knitr::kable'
  (https://github.com/Bioconductor/BiocStyle/issues/30)

- Correct padding of nested lists, and indentation of document date and
  multi-line headings

- Fix a bug in the function matching author affiliations
  (https://support.bioconductor.org/p/98155/)

- Proper horizontal alignment of ggplotly graphics
  (https://support.bioconductor.org/p/97609/)

- Fix compatibility with 'runtime: shiny' document option

[biomaRt](https://bioconductor.org/packages/biomaRt)
-------

Changes in version 2.34.0:

NEW FEATURES

- Added the listEnsemblArchives() function.  This returns a table of
  the available Ensembl archives, and replaces the archive = TRUE
  argument to several functions, which was no longer working.

BUG FIXES

- The Ensembl BioMart server doesn't always respond well if queries
  with more than 500 filter values are submitted.  If a query that
  exceed this is detect biomaRt will now submit the query in batches
  and concatonate the result when completed.

MINOR CHANGES

- You can now provide a host with 'http://' at the start, or a trailing
  '/' (typically copy/pasted from a browser) and useMarts() will cope.

[BiRewire](https://bioconductor.org/packages/BiRewire)
--------

Changes in version 3.9.1:

- Minorbug fixed (order not preserved using graph.edgelist and infinite
  loop in matrix with just one entry)

[bnbc](https://bioconductor.org/packages/bnbc)
----

Changes in version 0.99:

- Initial release to Bioconductor.

[branchpointer](https://bioconductor.org/packages/branchpointer)
-------------

Version: 1.3.0
Category: Improved speed for functions
Text:

Version: 1.3.0
Category: Functions combined and renamed
Text:

Version: 1.3.0
Category: Updated model to single gbm
Text:

Version: 1.3.0
Category: Bug fixes
Text:

Version: 1.3.1
Category: Version Bump
Text:

[BrowserViz](https://bioconductor.org/packages/BrowserViz)
----------

Changes in version 1.9.8:

BUG FIXES

- Significant (3x) speedup.  A 5000-node, 6000-edge graph transmits to
  Cytoscape from R in about 20 seconds.

[bsseq](https://bioconductor.org/packages/bsseq)
-----

Changes in version 1.13:

- 1.13.6: Fix performance regression in BSmooth(). Thanks to Shan
  Andrews for the report (<URL:
  https://github.com/kasperdanielhansen/bsseq/pull/57>).

- 1.13.5: Fix major bug in combine() and combineList().  *This bug led
  to bad BSseq objects with incorrect methylation estimates due to
  incorrect 'M', 'Cov', 'coef', and 'se.coef' assays.* To be safe,
  BSseq objects created with versions 1.13.0 to 1.13.4 should be
  re-created using a newer version. More specifically, any BSseq
  objects created with combine() or combineList() should be re-created.
  Also, BSseq objects created using read.bismark() or read.bsmooth()
  with multiple 'files' should be re-created.  Thanks to Alejandro
  Reyes (@areyesq89) for the report (<URL:
  https://github.com/kasperdanielhansen/bsseq/pull/54>).

- 1.13.4: Fix performance regression in getMeth() and getCoverage()
  when 'regions' were supplied. Thanks to Alejandro Reyes (@areyesq89)
  for the report (<URL: https://support.bioconductor.org/p/97611/>).

- Moved vignettes to Rmd.

[CAGEr](https://bioconductor.org/packages/CAGEr)
-----

Changes in version 1.18.1:

- getCTSS() can now load data from large BAM files.  Before this
  update, there were "negative length vectors" errors when filtering
  out low quality reads from datasets with multiple dozens of millions
  of sequences.

[CAMERA](https://bioconductor.org/packages/CAMERA)
------

Changes in version 1.33.3:

BUG FIXES

- Fix getPeaklist not pulling rownames correctly, closing #22, thanks
  to Jan Stanstrup

Changes in version 1.33.2:

BUG FIXES

- Fix getPeaklist, updated getReducedPeaklist(), Thanks to Kristian
  Peters

Changes in version 1.33.1:

NEW FEATURES

- allow matchedfilter alternative intval in getPeaklist (J. Stanstrup)

- pass intval from Groupcorr to calcCaS (J. Stanstrup)

- getReducedPeaklist returns just one intensity per pspec (K. Peters)

[canceR](https://bioconductor.org/packages/canceR)
------

Version: 3.1
Text:

Version: 3.2
Text: dimensions levels will be plot. 1- dialogMetOption(): add
        "Circos" argument to make the difference between
        getMetDataMultipleGene() and getListMetData() 2- getGeneList():
        add rm("GeneListMSigDB"", envir="myGlobalEnv")

[Cardinal](https://bioconductor.org/packages/Cardinal)
--------

Changes in version 1.9.2 (2017-10-25):

BUG FIXES

- Corrected author name in all documentation

Changes in version 1.9.1 (2017-10-23):

BUG FIXES

- Fixed bug in package dependency 'matter' affecting the size of
  datasets that can be processed with 'batchProcess'

- Fixed bug where bin sizes for units='ppm' were twice as wide as they
  should be in 'readImzML' and 'reduceDimension'

[categoryCompare](https://bioconductor.org/packages/categoryCompare)
---------------

Version: 1.21.4
Text:

[cbaf](https://bioconductor.org/packages/cbaf)
----

Changes in version 1.0.0 (2017-09-10):

New Functions

- New availableData() functon to scan all the cancer studies to examine
  presence of RNA-seq, microRNA-seq, microarray(mRNA),
  microarray(miRNA) and methylation data.

- New obtainOneStudy() function to obtain and store the supported data
  for at least one group of genes across multiple subgroups of a cancer
  study. In addion, it can check whether or not all genes are included
  in different subgroups of a cancer study and, if not, looks for the
  alternative gene names.

- New obtainMultipleStudies() function to obtain and store the
  supported data for at least one group of genes across multiple cancer
  studies. It can check whether or not all genes are included in each
  cancer study and, if not, it looks for the alternative gene names.

- New automatedStatistics() function to calculate the statistics of the
  data obtained by obtainOneStudy() or obtainMultipleStudies()
  functions. Based on user's preference, these statistics can include
  frequency percentage, frequency ratio, mean value and median value of
  samples greater than specific value. Furthermore, it can look for the
  genes that comprise the highest values in each cancer and list the
  top 5 genes for frequency percentage, mean value and median value.

- New heatmapOutput() function to prepare heatmap for frequency
  percentage, mean value and median value data provided by
  automatedStatistics() function. Heatmaps for every gene group are
  stored in separate folder.

- New xlsxOutput() function to export the output of
  automatedStatistics() and the gene validation result of one of the
  obtainOneStudy() or obtainMultipleStudies() functions as an excel
  file. For every gene group, an excel file will be generated and
  stored in the same folder as heatmaps.

- New cleanDatabase() function to remove the created databases in the
  cbaf package directory. This helps users to obtain the fresh data
  from cbioportal.org.

- New processOneStudy() function to combine obtainOneStudy(),
  automatedStatistics(), heatmapOutput() and xlsxOutput() functions for
  the ease of use.

- New processMultipleStudies() function to combine
  obtainMultipleStudies(), automatedStatistics(), heatmapOutput() and
  xlsxOutput() functions for the ease of use.

[ChAMP](https://bioconductor.org/packages/ChAMP)
-----

Changes in version 2.8.8:

- Added "force" parameter in champ.load(), which can be used for
  "minfi" loading method, if your data comes from different arrays,
  force parameter would allow minfi's read.meth.exp function to extract
  their commmon probes and continue analysis.

Changes in version 2.8.7:

- Make champ.import() more robust for modifed csv file.

Changes in version 2.8.5:

- champ.DMP() works on numeric variable now

- champ.DMP() do pairewise comparison between each two categorical
  phenotypes.

- goseq was replaced by gometh.

- DMP.GUI() is modified heavily.

- SVD plot added legend.

- "minfi" loading method fixed.

- vignette of ChAMP and github Demo pages updated.

- Added some figures from GSE40279 in vignette.

Changes in version 2.8.3:

- Updated zzz.R file, which means the loading messages would be
  different.

- Fixed a warning in champ.load.Rd file.

- Fixed bug in champ.filter(), if filerDetP is false, update pd part
  would faile because of lacking of RemainSample variable.

Changes in version 2.8.2:

- Updated EPIC annotation to B4 version. The B4 version is downloaded
  from illumina website.

- Added one new parameter "method" in champ.load() function, which
  allows user to choose which method they want to use to read data.
  ChAMP or Minfi.

- champ.filter() has been totally recoded, now user can do any
  filtering on any data set they want. Merely champ.filter() is focused
  to take champ.import() result as input and generate filtered beta
  value for future analysis.

- Provide Whole New function champ.import() to read IDAT file to R,
  which is similar to minfi's read.meth.exp() function.

- Added more strict checking in champ.runCombat(), now
  champ.runCombat() would check if your variable and batches conflict
  with each other.

- Removed some useless code in champ.DMR() to make it faster.

Changes in version 2.8.1:

- Added impute option for champ.load().

- Add ProbeCutoff and SampleCutoff parameters in champ.load().

- Added Demo on github: In respond to our reviewer's question and to
  make users have better understanding on our package, we processed
  ChAMP fully on some data sets and saved all messages shown during
  processing. We upload these information to
  [github](https://github.com/JoshuaTian/ChAMPDemos).

[ChemmineR](https://bioconductor.org/packages/ChemmineR)
---------

Changes in version 2.30.0:

NEW FEATURES

- Plots can be generated using OpenBabel library (requires ChemmineOB
  package and OpenBabel)

- Functions to query data from PubChem directly.

[chipenrich](https://bioconductor.org/packages/chipenrich)
----------

Changes in version 2.2.0:

NEW FEATURES

- polyenrich now supports weighting peaks by signal value

- A hybrid test, hybridenrich() is available for those unsure of which
  test, between chipenrich() and polyenrich() to use.

- A function to join two different results files, hybrid.join(), and it
  will give an adjusted set of p-values and FDR-adjusted p-values using
  the two.

- A new approximation method using the Score test is available for
  quick results for chipenrich and polyenrich. Only recommended for
  significantly enriched results, and not depleted results. ~30x
  faster.

IMPROVEMENTS

- Several updates to the vignette

[ChIPpeakAnno](https://bioconductor.org/packages/ChIPpeakAnno)
------------

Changes in version 3.11.8:

- make featureAlignedExtendSignal accept GAlignments list object

Changes in version 3.11.6:

- fix the bugs Codoc mismatches from documentation object
  'peakPermTest'

Changes in version 3.11.4:

- fix the bugs when annoMcols is numeric values.

Changes in version 3.11.3:

- update the function featureAlignedSignal

Changes in version 3.11.2:

- add new function binOverGene, binOverRegions, plotBinOverRegions.

Changes in version 3.11.1:

- make FAQs availble

- Fix the bug in annoPeaks for warning of out-of-bound ranges

[ChIPQC](https://bioconductor.org/packages/ChIPQC)
------

Changes in version 1.13.2:

- Fixed bug when specifying consensus peakset as GRanges

[ChIPseeker](https://bioconductor.org/packages/ChIPseeker)
----------

Changes in version 1.13.1:

- fixed issue of naming intronList <2017-07-06, Thu> +
  https://github.com/GuangchuangYu/ChIPseeker/issues/57#issuecomment-313342399

[chromstaR](https://bioconductor.org/packages/chromstaR)
---------

Changes in version 1.3.1:

NEW FEATURES

- New parameter 'stepsize' allows sliding bins. This improves
  localization of peaks.

SIGNIFICANT USER-LEVEL CHANGES

- New default value for Chromstar(..., stepsize = 1/2 * binsize).

[ClassifyR](https://bioconductor.org/packages/ClassifyR)
---------

Changes in version 1.12.0:

- Alterations to make plots compatible with ggplot versions 2.2 and
  greater.

- calcPerformance can calculate some performance metrics for
  classification tasks based on datasets with more than two classes.

- Sample-wise metrics, like sample-specific error rate and
  sample-specific accuracy are calculated by calcPerformance and added
  to the ClassifyResult object, rather than by samplesMetricMap and
  being inaccessible to the end-user.

[cleaver](https://bioconductor.org/packages/cleaver)
-------

Changes in version 1.15.1 (2017-06-08):

- Don't consider the last aminoacid as cleavage site.
  `cleavageSites("ACK", custom="K")` results in `integer()` instead of
  `3` now. Thanks Apurva Hegde <ahegde@tgen.org> for reporting this
  issue.

- Remove superfluous "missedCleavages" argument in `cleavageSites`.

[clustComp](https://bioconductor.org/packages/clustComp)
---------

Changes in version 1.5.2:

- A new feature is added to flatVShier(); if the expanded version of
  the dendrogram is plotted, a coloured bar on the right hand side
  displays how genes are distributed across the flat clusters.

- New arguments ‘bar1.col’ and ‘bar2.col’ in flatVShier() allow tuning
  the coloured bars when the expanded dendrogram is plotted.

- Argument ‘weights’ in flatVSflat() is replaced by ‘flat1’ and ‘flat2’
  for analogy with flatVShier().

- New argument ‘greedy’ in flatVSflat() allows displaying the
  one-to-one mapping between superclusters as in flatVShier() after
  finding the best layout for the bi-graph.

- New argument ‘greedy.colours’ allows setting up the visualisation of
  the mapping between superclusters, as in flatVShier().

- The visualisation of superclusters in SCmapping() includes new labels
  with their sizes.

- Internal function drawTreeGraph() includes new arguments ‘flat.obj’,
  ‘bar1.col’ and ‘bar2.col’ to display the second coloured bar and
  control the colours.

[clusterExperiment](https://bioconductor.org/packages/clusterExperiment)
-----------------

Changes in version 1.3.7 (2017-10-24):

Changes

- Added function `tableClusters` for tabling clusters by name.

- Added largeDataset option to subsampleClustering

- allow `mergeInfo` to be called in `plotDendrogram` to use stored
  merge info in plotting.

Bugs

- Fixed bug in .makeMergeDendrogram in getFeatures

Changes in version 1.3.6 (2017-10-18):

Changes

- MAJOR CHANGE TO DEFINITION OF CLASS: Added slots to
  `ClusterExperiment` object so that the object keeps the information
  about the merging, and added corresponding helper functions (see
  documentation). This means previous versions made of a
  `ClusterExperiment` object will no longer be valid objects! Old,
  saved objects from previous versions must be manually adapted to have
  these slots (with `NULL` or `NA` values as appropriate).

- Add function `plotClustersWorkflow`, see documentation

- Add function `plotDimReduce` for plotting low-dim pca with points
  labeled by cluster.

- Add function `plotContrastHeatmap` to plot the significant genes that
  are result of getBestFeatures

- `getBestFeatures`: can now be run on result of mergeClusters without
  having to call makeDendrogram again for the merged clusters.

- Change argument to `plotClusters` and `plotBarplot` from `clusters`
  to `object`

- `makeDendrogram`: default in is now `dimReduce="mad"` to avoid
  accidentally calling it with `dimReduce="none"` (previous default),
  which can kill interactive R sessions.

- Changes to `plotHeatmap`: - default in `plotHeatmap` to
  `clusterSamplesData="dendrogramValue"`. - Changed handling of
  `clusterSamplesData` argument in `plotHeatmap` so that if argument
  equals `dendrogramValue` or `primaryCluster` and gives bad results /
  errors, will give warning and change argument appropriately.  - Added
  arguments to `plotHeatmap` allowing more user control regarding
  making blank lines for when have gene groupings.

- `mergeClusters` now aligns the colors from `mergeClusters` and
  `combineMany` internally.

- `mergeClusters` now suppresses warnings created by the estimation of
  the percentage non-null unless `showWarnings=TRUE`.

- `plotDendrogram`: - has new argument `clusterLabelAngle` allowing
  user to change the angle of the clusterLabel printed on top (when
  plot is of type "colorblock") - argument `labelType` has been changed
  to `plotType`

Bugs

- Arguments passed to `mergeClusters`, ClusterExperiment version, via
  the `...` command will now go first to the `mergeClusters` matrix
  version and then onto the plotting, as stated in documentation
  (plotting arguments were being ignored on the clusterExperiment
  version of the function).

Changes in version 1.3.4 (2017-09-28):

Changes

- Add argument `clusterLabels` to `plotClusters` to allow user to set
  clusterLabels without changing the object.

- Added data object of run of RSEC. Changed `lazyLoad` to `false`
  because loading this object on installation was creating errors.

- Updated vignette to more completely cover RSEC.

- Change default `clusterFunction` argument to `RSEC` to be
  "hierarchical01" rather than all 01 methods (previous default). This
  reduces load of making simple call.

- Updated validity checks to reduce memory load, and dropped
  `validObject` calls.

- Added function `getClusterManyParams` to parse the parameter values
  in the clusterLabels from clusterMany

- Removed old dependency on diagram (from previous vignette, no longer
  needed)

Bugs

- Fixed error in subsetting of clusterExperiment object when dendrogram
  is attached (previously didn't reset dendro_outbranch to NA)

Changes in version 1.3.3 (2017-09-07):

Changes

- Bug fix in clusterContrasts -- missing `match.arg` option for
  `outputType` argument

Changes in version 1.3.2 (2017-07-05):

Changes

- Default for `top.can` in seqCluster are changed to be `top.can=5`.

- makeDendrogram now has the default argument
  `ignoreUnassignedVar=TRUE` like in RSEC

- add ClusterFunction class and update all functions to work on this.
  All built in cluster functions are now given ClusterFunction Objects,
  so all built in clustering functions can now work for either
  `subsampleClustering` or `mainClustering`. This will also make it
  easier for a user to define their own ClusterFunction object and have
  it be used by functions like `clusterSingle`. This is a major change
  in how some of the underlying functions work, but should not impact
  common functions like `clusterMany` and `RSEC`. Some of the more
  notable changes in the arguments for programmers are: - `clusterD`
  and `clusterDArgs` have been changed to `mainClustering` and
  `mainClusterArgs`. This change was made to make these arguments more
  clear as to their role in the clustering workflow (and because the
  clusterD refered to clustering a dissimilarity but it has clustered
  either x or D for many versions now. ) - `seqCluster` and
  `clusterSingle` no longer take the argument `clusterFunction`.
  `clusterFunction` must be given via `mainClusterArgs` and
  `subsampleArgs` to be passed to `mainClustering` or
  `subsamplingCluster`, respectively. Now only the upper-level function
  `clusterMany` takes `clusterFunction` directly as an argument. -
  `mainClustering` (previously `clusterD`) and `subsampleClustering` no
  longer take `k` nor `alpha` as a direct argument. These arguments,
  like all arguments used directly by the cluster function, now need to
  be passed to the clustering function in a list via `clusterArgs`.  -
  The list of available built-in clustering routines provided by the
  package can now be accessed via `listBuiltInFunctions()`. The
  functions that are used for these functions are now available to the
  user via their ClusterFunction object that the user can access with
  the function `getBuiltInFunction`. (see ?listBuiltInFunctions)

- `hiearchical01` clustering now has a different default, namely to
  apply `as.dist` to the input `diss` in order to get a `dist` object,
  rather than `dist(1-diss)` which was previously the default for
  historical reasons. This is controlled by argument `whichHierDist`,
  and can be set to the previous version by passing
  `whichHierDist="dist"` to the `clusterArgs` argument in either
  `subsampleArgs` or `mainClusterArgs`, depending on where
  `hierarchical01` is being used.

- Spectral clustering is now available (`"spectral"`) via the `specc`
  function of `kernlab`.

- `clusterSingle` now only returns the dissimilarity matrix in the
  `coClustering` slot if `subsample=TRUE` in the call. Furthermore, for
  the resulting dissimilarity to replace an existing `coClustering`
  slot value, the user must request it by setting
  `replaceCoClustering=TRUE` in the call to `clusterSingle`.

- Removed default value for argument `proportion` in `combineMany`. The
  previous default was `proportion=1` but didn't match most common use
  cases and was accidentally called by upper functions like RSEC.

- If the `clusterFunction` argument is not given to `subsampleArgs` by
  the user explicitly, and the `clusterFunction` of `mainClusterArgs`
  is appropriate, it will be used for `subsampleClustering`; if the
  `clusterFunction` in `mainClusterArgs` is not appropriate (e.g.
  `subsampleClustering` needs a type `K` because `sequential=TRUE`),
  then the default for the `subsampleClustering` will be `'pam'`. This
  changes the previous behavior of `subsampleClustering` where the
  default was 'pam' in all cases where not explicitly given by the
  user. This change should have no impact on RSEC: since the
  `clusterFunction` for the `mainClustering` terms is a `'01'` type in
  RSEC and the `subsampleClustering` has to be type `'K'` when
  `sequential=TRUE`, it will revert to the default `"pam"` as before.

Bugs

- Fixed error so where if `clusterSingle` was called on existing
  clusterExperiment object it would overwrite the information of the
  existing `clusterExperiment` object.

- Fixed `RSEC` so now if rerun on existing `clusterExperiment` object,
  it will grab defaults from the matrix version (previously defaults
  were those of the underlying function, which were not always the
  same, e.g. `combineProportion` default was previously 1)

- Fixed `clusterMany` so now it explicitly sets `dimReduce="none"` in
  call to `clusterSingle`. Before, might have been calling all of the
  `dimReduce` defaults (i.e. all of them!).

- Fixed so gives error if whichClusters in `plotBarplot` doesn't match
  anything.

Changes in version 1.3.1 (2017-06-14):

Changes

- change how `plotHeatmap` handles visualizeData argument, so not
  required to have same number of genes as original, only same number
  of samples.

- Now if color of vectors given in `clusterLegend` does not have names,
  `plotHeatmap` will give them names matching the variable so that they
  will be used by `aheatmap` (previously would have left all colors
  white because do not have matching names).

- Large changes to how dendrograms are plotted by `plotDendrogram` and
  `mergeClusters`. This includes the ability to see the before and
  after clusterings along side the mergeClusters result, as well as a
  new slot added to the clusterExperiment class (`dendro_outbranch`).
  The names of several arguments to `mergeClusters` and
  `plotDendrogram` were changed for clarity: - `leaves` is now
  `leafType` in `plotDendrogram`.  - `plotType` is now `plotInfo` in
  `mergeClusters` - `doPlot` is now `plot` in `mergeClusters` -
  `leafType` is now an option for `mergeClusters` as well. - Now when
  `plotInfo` (previously `plotType`) is set to `none`, the plot is
  still drawn, but just no information about the merging is added to
  the plot. To not plot the dendrogram at all, set `plot=FALSE`. - The
  option `labelType` in either `plotDendrogram` or `mergeClusters`
  controls whether names (`name`) or rectangular color blocks
  corresponding to the cluster (`colorblock`) are put at the tips of
  the dendrogram to label the clusters/samples.

- added `dendroClusterIndex` that behaves similarly to
  `primaryClusterIndex`

- added ability to give `dendro` as charater option to `whichClusters`
  argument

- added `transformation<-` to be able to assign manually the
  transformation slot

- Move MAST into 'suggests' pacakge so that not need R 3.4 to run the
  package.

- Change calculation of PCA dimensionality reduction to use `svds` from
  `RSpectra` package to improve speed

Bugs

- Fixed bug in RSEC where `combineProportion` argument was being
  ignored (set to 1)

- Fixed bug in definition of `transform` so that extends existing
  generic rather than masking it.

Changes in version 1.3.0 (2017-05-24):

Changes

- `plotHeatmap` accepts `data.frame` or `ExpressionSet` objects for the
  data argument (calls `data.matrix` or `exprs` on object and sends to
  matrix version)

- Added `plotBarplot` to plot a barplot for 1 cluster or comparison of
  2 clusters along with tests.

- Added `whichClusters` argument to `clusterMatrix` to return only
  clusters corresponding to certain clusters. Mainly relevant for using
  arguments like `workflow` that are used by other commands (otherwise
  could just index the complete matrix manually...)

Bug fixes

- `plotHeatmap` now goes through the `clusterLegend` input and removes
  levels that do not exist in the sampleData; this was causing
  incorrect coloring when the `clusterLegend` had more (or less) levels
  that it assigned color to than the `sampleData` did (e.g. if
  `sampleData` was a subset of larger dataset upon which the original
  colors were assigned.) NOTE: that this now has the effect of NOT
  plotting all values in the clusterLegend if they are not represented
  in the data, thus changing the previous behavior of `plotHeatmap`
  legend.

- fixed bug in how `plotHeatmap` checked that the dimensions of
  user-supplied dendrogram match that of data (matrix version).

- fixed `convertClusterLegend` so when `output` is `matrixNames` or
  `matrixColors`, the resulting matrix has the `colnames` equal to
  cluster labels, like `clusterMatrix`.

- internal function .convertToNum now preserves names of input vector.

- fixed bug in plotting with merge clusters; previously if
  plotType="all", might not have been correctly plotted with the right
  internal node of the dendrogram.

[ClusterJudge](https://bioconductor.org/packages/ClusterJudge)
------------

Changes in version 0.99:

USER VISIBLE CHANGES

- added an Introductory vignette

[clusterProfiler](https://bioconductor.org/packages/clusterProfiler)
---------------

Changes in version 3.5.6:

- fixed R check <2017-09-28, Thu>

Changes in version 3.5.5:

- ko2name <2017-08-14, Mon>

- bioc git transition <2017-07-18, Tue>

Changes in version 3.5.4:

- change keytype parameter in enrichGO to keyType for consistent
  <2017-06-28, Wed> +
  https://github.com/GuangchuangYu/clusterProfiler/issues/91#issuecomment-311657933

- bug fixed of simplify for compareCluster result <2017-06-19, Mon> +
  https://github.com/GuangchuangYu/clusterProfiler/issues/67#issuecomment-307815402

Changes in version 3.5.3:

- fixed https://github.com/GuangchuangYu/clusterProfiler/issues/88
  <2017-05-24, Wed>

- accept background gene list via `universe` parameter in `enrichDAVID`
  <2017-05-16, Tue> +
  https://github.com/GuangchuangYu/clusterProfiler/issues/87

Changes in version 3.5.2:

- support keyType in enrichMKEGG <2017-05-10, Wed>

- bug fixed of maintaining input list order in plotting compareCluster
  result <2017-05-02, Tue> +
  https://support.bioconductor.org/p/95400/#95426

Changes in version 3.5.1:

- bug fixed of converting KEGG PATH ID to NAME using KEGG.db
  <2017-04-28, Fri> +
  https://github.com/GuangchuangYu/clusterProfiler/issues/85

[ClusterSignificance](https://bioconductor.org/packages/ClusterSignificance)
-------------------

Changes in version 1.5.3:

- Fixed bug in pcp 2D projection plot where colors did not match
  classification plot.

Changes in version 1.5.2:

- Fixed typo in vignette.

[CNEr](https://bioconductor.org/packages/CNEr)
----

Changes in version 3.5:

NEW FEATURES

- Add function orgKEGGIds2EntrezIDs to fetch the mapping between KEGG
  IDs and Entrez IDs

- Add function makeAxtTracks

- Add function addAncestorGO

Changes in version 3.4:

NEW FEATURES

- Updated CNE class for storing all the information about running the
  pipeline.

- Add read.rmMask.GRanges to read RepeatMasker .out file.

- Add read.rmskFasta to read soft repeat masked fasta file.

- Add the distribution plot of axt alignment matches.

- Add the distribution plot of CNE length.

- Add the syntenicDotplot for axt alignment and GRangePairs.

- When readAxt and readBed, seqinfo is kept when available.

- New fixCoordinates function makes the coordinates of Axt alignments
  always relative to positive strands.

- Add parallel subAxt for Axt alignment.

- Add the plot of genomic distribution of CNE.

- Add the function to make bed and bigwig files of CNEs.

BUG FIXES

- Instead of an error, an empty GRangePairs is returned when no CNEs
  identified.

Changes in version 3.3:

BUG FIXES

- Fix a bug caused by "format" in the blat step.

NEW FEATURES

- Add the pairwise whole genome alignment pipeline

- Add a new class "GRangePairs"

- "Axt" class is now based on "GRangePairs" class.

- readAncora for reading Ancora format CNE files.

[coexnet](https://bioconductor.org/packages/coexnet)
-------

Version: 0.1.0
Category: NEW FEATURES
Text: 18/04/17 --> changing the CCP function.  23/06/17 --> Released in
        Biocondonductor-devel 12/10/17 --> geneSymbol error test fixed.

[COMPASS](https://bioconductor.org/packages/COMPASS)
-------

Changes in version 1.15.1:

- Exported SimpleCOMPASS interface for use with count matrices

- Previous versions - added translate_marker_names API for use with
  SimpleCOMPASS

[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap)
--------------

Changes in version 1.15.1:

- random colors are generated by new `rand_color()` function in
  circlize package.

- add `density_param` in `densityHeatmap()` function

- annotations with duplicated names have no legends any more

- re-implement `grid.xaxis()` to draw axis labels rotated 90 degrees

- grids in discrete legend are arranged by rows if ncol > 1

- raster image is generated in an independent R session

- empty string in annotation or heatmap is mapped to NA

- annotation and heatmap legends can be merged into one column.

- change the default value of `row_names_max_width` and
  `column_names_max_height`

- default legend style for continuous values is changed to "continuous"

- add `grid.dendrogram2()` which draws dendrograms with uneven position
  for leaves

- move **dendextend** to Suggests field because it depends/imports
  rlang indirectly which has a `print.frame()` function and it will
  affect to print a `frame` object returned by `frameGrob()`.

- `decorate_*()` functions return to the viewport where they are
  called.

- fixed a bug that annotation names are drawn for all row slices.

- construct a valid path under Windows

[cosmiq](https://bioconductor.org/packages/cosmiq)
------

Changes in version 1.11.3 (2017-06-19):

- added vignettes/poster.Rmd for 13th Annual Conference of the
  Metabolomics Society <URL: metabolomics2017.org> (in portrait).

[CountClust](https://bioconductor.org/packages/CountClust)
----------

Changes in version 1.5.1:

- Release We have added a FitGoMpool() function that automatically
  performs GoM model with multiple starting points and outputs the one
  run with the most optimal BIC. Besides, removed
  switch_axis_position() as a dependency from cowplot as the function
  has been deprecated. r

[CrispRVariants](https://bioconductor.org/packages/CrispRVariants)
--------------

Changes in version 1.5.9:

- More comprehensive input checking in readsToTarget, removed redundant
  checking from CrisprSet initializer.

- Fix bug where sequences falsely called no variant if target is on the
  negative strand but positive strand reference given.

- Changed default SNV calling to 6 bases downstream instead of 5 to
  cover PAM

- Added tests for mismatched reference and target

Changes in version 1.5.8:

- Adds an option to filter variants by name when counting or plotting

Changes in version 1.5.7:

- Fixes major bug preventing filtering in plotFreqHeatmap

- Fixes bug in mergeChimeras if no chimeras mergeable

Changes in version 1.5.6:

- Autogenerate bam index for readsToTarget option chimeras = "ignore"

Changes in version 1.5.3:

- plotAlignments now accepts the same filtering arguments as
  variantCounts

Changes in version 1.5.1:

- New argument alleles in plotAlignments and plotFreqHeatmap for
  selecting which alleles to display or specifiying a plotting order

- Removed unnecessary fields from CrisprRun class

- New accessor function alns to get alignments from a CrisprSet

- Improvements to plotAlignments to avoid unnecessary symbols in legend

- Fix to header of plotFreqHeatmap when using type = "proportions" and
  providing sample order

- consensusSeqs now returns cigars as metadata by default

- Started indenting with four spaces at the start

[csaw](https://bioconductor.org/packages/csaw)
----

Changes in version 1.11.4:

- Removed support for paramList objects.

- Added option for normOffsets() to return SummarizedExperiment objects
  containing normalization data.

- Deprecated normalize() to avoid S4 method clashes.

- Moved scaling prior to control-based filtering into a new function,
  scaleControlFilter(), for greater modularity.

- Updated user's guide.

[cydar](https://bioconductor.org/packages/cydar)
-----

Changes in version 1.1.4:

- Added labelSpheres() function for labelling unannotated hyperspheres.

- Exported multiIntHist() for plotting multiple intensity histograms.

- Slight fix to spatialFDR(), which now computes the correct n-th
  nearest neighbour.

[cytofkit](https://bioconductor.org/packages/cytofkit)
--------

Changes in version 1.9.5 (2017-10-23):

MODIFICATIONS

- Standardise plot scales across samples in shiny App

BUG FIXES

- Fixed an issue with Phenograph failing if sampling with replacement
  was done. Now cytofkit tests if any FCS files have less events than
  specified in fixedNum argument.

Changes in version 1.9.4 (2017-09-27):

MODIFICATIONS

- Added function, cytof_clustermtrx(), to obtain marker expression
  values for cells in a given cluster.

- Included cytof_clustermtrx() into main cytofkit function.

- Cluster ID integers saved to fcs files as additional channels.

BUG FIXES

- Quick fix on FCS saving.

Changes in version 1.9.3 (2017-09-27):

MODIFICATIONS

- Added arguments to cytofkit() and cytof_cluster() to allow user to
  set a seed, for reproducible results.

BUG FIXES

- Users can now save .fcs from shiny output after renaming samples.

Changes in version 1.9.2 (2017-09-11):

MODIFICATIONS

- Added save data options for shiny app. Now you can select any
  combination of .fcs, .csv, and .rdata outputs.

- cytofkitShinyAPP() uses code from cytofkit_shinyAPP.R for
  visualisation on local machines, while ui.R and server.R can be used
  separately to host the app on servers.

BUG FIXES

- Fixed some issues with viewing and downloading the marker heatmap in
  the shiny app

- Fixed an inconsistency where expression data used for clustering was
  not selective for markers used for dimension reduction

Changes in version 1.8.3 (2017-09-06):

BUG FIXES

- cytofkitShinyAPP file size cap consistently set to 100mb

Changes in version 1.8.2 (2017-07-24):

NEW FEATURES

- cytofkitShinyAPP function now takes RData as argument to skip
  reuploading RData

- While choosing selected markers for dimension reduction and
  clustering, all markers can be visualised in the shiny app

- Added a "reset" button to ShinyApp to clear the session and start
  over

- Added Server file select button to allow browsing server files

MODIFICATIONS

- Combined both cytofkitShinyAPP functions into one

- ShinyApp uses actionButton for download instead of a downloadHandler

- ShinyApp displays what .RData is loaded into its reactive data

- ShinyApp lists markers in alphabetical order

- Marker selection done at dimensionality reduction stage rather than
  at raw data transformation, to allow all expression data to be
  visualised at later stages

- Updated examples and vignettes to account for updates done

- Updated maintainer email address

Changes in version 1.8.1 (2017-04-26):

NEW FEATURES

- fixed documentation warning for function cytofkitShinyAPP2

- updated my maintainer email address

[DAPAR](https://bioconductor.org/packages/DAPAR)
-----

Changes in version 1.9.15:

BUG FIXES

- When the aggregation step has been performed, the interface switches
  to the first tab of the 'Descriptive Statistics' in order to view
  informations aout the new dataset (the protein one).

- A new package (readxl) is used to read xls or xlxs files. In certain
  circumstances, the functions of the previsous package openxlsx is not
  able to decode properly Excel files.

- When converting a new (text or Excel) file in Prostar : the missing
  values were not registered as expected. Especially, they did not
  appear in blue in the table above the volcanoplot. Bug fixed

NEW FEATURES

- Implementation of a parallel version of the function which saves the
  (new) protein dataset after the aggregation step.

- Enhancement of the string-based filtering UI

- The automatic generation of an analysis report has been integrated in
  the Dataset Manager (menu 'Export'). It allows the user to download
  plots and parameters used in Prostar ont their dataset.

- Implementation of a parallel version of the function which saves the
  (new) protein dataset after the aggregation step.

- Added a Gene Ontology (GO) analysis module in Data Processing. This
  module allows to perform GO grouping and GO Enrichment

- Several plots are now based on the package highcharter which is a
  wrapper to the highcharts graphical library. It provides
  interactivity with the user.

[DASC](https://bioconductor.org/packages/DASC)
----

Changes in version 0.1.0:

- Requirements: Please install NMF, cvxclustr, Biobase R package before
  the installation of DASC package.

- Workflow: The function convexBatch() performs all steps for the
  detction of the batches in the dataset.

[DChIPRep](https://bioconductor.org/packages/DChIPRep)
--------

Changes in version 1.7.2:

- added betaPrior = TRUE to nbinomWaldTest in order to keep the
  shrinkage of the log2 fold changes. (The shrinkage is now disabled by
  default in DESeq2)

- Changed R dependency to R 3.4

- Added WholeGenome Bioc View

Changes in version 1.7.1:

- extended the importing section in the vignette a little bit

- fixed the robust mean function so that now plotting also works
  without replicates

- new function robust_mean that is used in plotting

- when importing matrices, it is now checked that their column names
  correspond to the sample IDs given

[debrowser](https://bioconductor.org/packages/debrowser)
---------

Changes in version 1.5.4:

- heatmap redblue fix

Changes in version 1.5.3:

- Figure caption fix

- Menu fix

Changes in version 1.5.2:

- img path build fix

[deepSNV](https://bioconductor.org/packages/deepSNV)
-------

Changes in version 1.99.3 (2013-07-25):

Updates

- A few changes to shearwater vignette

- Renamed arguments pi.gene and pi.backgr in makePrior()

Bugfixes

- Fixed bug in bf2Vcf() when no variant is called

Changes in version 1.99.2 (2013-07-11):

Updates

- Updated CITATION

- Added verbose option to bam2R to suppress output

- Changed mode() to "integer" for value of loadAllData()

Bugfixes

- Fixed bug when only one variant is called in bf2Vcf()

Changes in version 1.99.1 (2013-06-25):

Updates

- Using knitr for prettier vignettes

- Including shearwater vignette

Bugfixes

- fixed issues with deletions in bf2Vcf()

- makePrior() adds background on all sites

Changes in version 1.99.0 (2013-04-30):

Updates

- New shearwater algorithm

- Including VCF output through summary(deepSNV, value="VCF")

[DEGreport](https://bioconductor.org/packages/DEGreport)
---------

Version: 1.13.12
Text: 2017-10-11 Lorena Pantano <lorena.pantano@gmail.com> Feature:
        Return scatter plots between PCs and metadata in degCovariates.
        Feature: Use ConsensusClusterPlus to cluster genes with
        degPatterns.

Version: 1.13.11
Text: 2017-09-22 Lorena Pantano <lorena.pantano@gmail.com> Feature:
        significant works with DESeqResults class. Fix: log2 in degPlot
        wasn't active. Feature: Allow plot samples together or not in
        degCheckFactor. Feature: Migrate vignette to new BiocStyle.
        Fix: Automatic QC report. Reduce final report with most
        important figures.

Version: 1.13.10
Text: 2017-09-07 Lorena Pantano <lorena.pantano@gmail.com> Fix:
        Complete vignette with new functions.  Feature: Add DEGSet
        construct to accept other sources.  Feature: Adapt degQC to
        accept DEGSet object. Feature: Allow multiple group for degMB
        and degVB.  Feature: Add optional log2 for gene plotting.
        Feature: Plot correlation of shrunken vs unshruken log2fc.
        Feature: Allow to plot original MA plot.  Feature: Adapt
        summary of DESeq2Results to data.frame and compatible with
        markdown output, and multiple alpha values.  Fix: links in man
        pages.  Feature: Pass options to Heatmap in
        degCorCov.@vbarrera.  Feature: Add parameter to select top rows
        from DEGset.  Fix: Change method names to short words.
        Feature: degVolcano accepts DESeq2Results class.  Fix: degPCA
        print the correct PC number on x/ylabels.  Fix: Move NEWS to
        parent folder. Feature: Add method to get significant genes
        from DEGSet class. Feature: Add plotMA method to show shrunken
        effect. @vbarrera Fix: Move to testthat for examples.  Feature:
        Adding main class and methods to handle DEG output.  Fix: axis
        in degPCA now show the values.  Feature: Handle multiple
        contrasts/coefficient for DESeq2 results.

Version: 1.13.7
Text: 2017-08-08 Lorena Pantano <lorena.pantano@gmail.com> Feature: new
        function to analyze the correlation among covaritaes in
        metatdata Deprecation: all functions related to foldchange
        accuraty are removed.  Using lfcShrink much better now Style:
        Add more unit tests Feature: Accept SE like objects to degPlot
        and use better gene names if rowData has it Feature: Use
        plot_grid for degPattern and save plot Feature: Use text or
        point in degPCA Feature: Fix labels of degPlot Feature: Accept
        matrix in degPlot Fix: correctly handling rowData in SE objects
        for degPlot Fix: plot only legend if group > 1 Feature: More
        output for degPattern Style: change to lower-cases inside
        degCovariates function

Version: 1.13.6
Text: 2017-05-30 Lorena Pantano <lorena.pantano@gmail.com> Feature: Add
        degPCA plot from Rory Kirchner

Version: 1.13.5
Text: 2017-05-30 Lorena Pantano <lorena.pantano@gmail.com> Feature:
        Accept matrix for degWidePlot

Version: 1.13.4
Text: 2017-05-22 Lorena Pantano <lorena.pantano@gmail.com> Feature: Add
        degCovariates to calculate correlations between PCs from count
        data and covariates from metadata

Version: 1.13.3
Text: 2017-05-19 Lorena Pantano <lorena.pantano@gmail.com> Feature: Add
        degMDS for PCA like clustering figures Feature: Add labels
        parameters to degPlot

Version: 1.13.2
Text: 2017-05-05 Lorena Pantano <lorena.pantano@gmail.com> Feature: add
        degFilter to filter genes by group

Version: 1.13.1
Text: 2017-04-27 Lorena Pantano <lorena.pantano@gmail.com> Feature:
        change arrange for plotQC plots with cowplot Feature: Use
        theme_minimal inside degResults Feature: Change title for some
        sections in degResults

[DEP](https://bioconductor.org/packages/DEP)
---

Changes in version 1.0.0:

- First release of the DEP package for differential enrichment analysis
  of proteomics data.

[derfinder](https://bioconductor.org/packages/derfinder)
---------

Changes in version 1.11.8:

BUG FIXES

- Improved the documentation regarding an error when
  coverageInfo$position is NULL when running analyzeChr() [and
  indirectly running preprocessCoverage()]. See
  https://support.bioconductor.org/p/99400/ for details.

Changes in version 1.11.7:

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette now uses the new BiocStyle::html_document that was recently
  released.

Changes in version 1.11.4:

BUG FIXES

- regionMatrix() will now pass the hidden arguments 'species' and
  'currentStyle' to getRegionCoverage() so they can be used by
  extendedMapSeqlevels(). Related to
  https://support.bioconductor.org/p/95721/.

Changes in version 1.11.2:

BUG FIXES

- Improved documentation of extendedMapSeqlevels(). Related to
  https://support.bioconductor.org/p/95521/.

- Improved filterData() based on
  https://github.com/lcolladotor/derfinder/issues/38

[derfinderHelper](https://bioconductor.org/packages/derfinderHelper)
---------------

Changes in version 1.11.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette now uses the new BiocStyle::html_document that was recently
  released.

[derfinderPlot](https://bioconductor.org/packages/derfinderPlot)
-------------

Changes in version 1.11.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette now uses the new BiocStyle::html_document that was recently
  released.

[DESeq2](https://bioconductor.org/packages/DESeq2)
------

Changes in version 1.18.0:

- lfcShrink() offers alternative estimators type="apeglm" and
  type="ashr", making use of shrinkage estimators in the 'apeglm' and
  'ashr' packages, respectively. See ?lfcShrink for more details and
  appropriate references. The integration of these alternative
  shrinkage estimators is still in development. Additionally, the
  DESeqResults object gains priorInfo(res), which passes along details
  of the fitted prior on LFC.

- Factor levels using characters other than letters, numbers, '_' and
  '.' will print a message (not a warning or error) that it is
  recommended to restrict to these "safe characters". This follows a
  suggestion from the Bioconductor support site to avoid user errors.

[DEsubs](https://bioconductor.org/packages/DEsubs)
------

Version: 1.3.4
Text: the package and the output in the landing page.

Version: 1.3.3
Text: Placed external figure importing within chunks.

Version: 1.3.2
Text:

Version: 1.3.1
Text:

[DiffBind](https://bioconductor.org/packages/DiffBind)
--------

Changes in version 2.6.0:

- Feature changes
  
  * Change sortFun parameter in dba.plotHeatmap to default sd, with
  FALSE as option for no sorting
  
  * Change sortFun parameter in dba.plotHeatmap to default sd, with
  FALSE as option for no sorting
  
  * Sort peaks when adding directly via dba.peakset
  
  * don't add _ if no initString in dba.report
  
  * Internal feature: alternate peak counts: pv.resetCounts

- Bug Fixes
  
  * Fix bug in dba.plotHeatmap is all values in a row are zero
  
  * Change authors in vignette to conform to new standard
  
  * Fix single peak boundary conditions
  
  * Subset config$fragmentSize when masking
  
  * Update example peaks to match data objects
  
  * Fix bug in dba.plotHeatmap if all values in a row are zero
  
  * Bugfix when returning report as GRanges with only one site
  
  * Bugfix when plotting venns of consensus peaksets
  
  * Fix buffer overrun causing segfault on MacOS

[diffHic](https://bioconductor.org/packages/diffHic)
-------

Changes in version 1.9.9:

- Added extractPatch() function to count bin pairs in a specified area
  of the interaction space.

- Modified connectCounts() to eliminate warnings upon stranded entries,
  unknown chromosomes. All entries of input regions are now retained,
  though not necessarily in the input order. Also switched original
  metadata to NA when second.regions is an integer.

- Modified preparePairs() to be more generous when considering
  inward-facing reads if they overlap past each other.

- Fixed bug in savePairs() involving failure to swap other information
  when enforcing index ordering.

- Added mergeCMs() function to allow entry into the pipeline from
  ContactMatrix objects.

- Moved pre-processing scripts out of the package to the repository for
  the user's guide.

- Updated presplit_map.py to use new samtools API for sorting.

- Updated user's guide.

[diffuStats](https://bioconductor.org/packages/diffuStats)
----------

Changes in version 1.0.0:

- Five diffusion kernels available, they can be computed from an
  'igraph' object.

- Diffusion implementations divided between 'diffuse_raw' for
  deterministic scores and 'diffuse_mc' for permutation analysis, which
  is parallelised. In total, seven diffusion scores are accessible
  through the 'diffuse' function.

- Performance evaluation wrapped in the 'perf' function.

- Helper functions in helpers.R (to plot diffusion scores, to check if
  a kernel matrix is actually a kernel, to extract largest CC from a
  graph)

[Director](https://bioconductor.org/packages/Director)
--------

Changes in version 1.3.1:

- Improved notice of discovered feedback loops.

- Added example input with feedback loop to manual.

[DMCHMM](https://bioconductor.org/packages/DMCHMM)
------

Changes in version 0.99.10:

CHANGES SINCE LAST TIME

- Several issues are fixed.

- Some modification are done to expedise the process.

- Instead of the parallel package, the BioParallel package is used.

Changes in version 0.99.9:

CHANGES SINCE LAST TIME

- Some parameters in examples are changed to make the running time
  faster.

- The BSData-method is changed to cBSData-method.

- The BSDMCs-method is changed to cBSDMCs-method. BUG FIXES

- None reported.

Changes in version 0.99.0:

IMPROVEMENTS SINCE LAST TIME

- We have refined some of the codes to speed up the running the
  package.

FUTURE DEVELOPMENT

- Some functions will be rewitten in C.

BUG FIXES

- None reported so far.

[DOSE](https://bioconductor.org/packages/DOSE)
----

Changes in version 3.3.2:

- new project site using blogdown <2017-09-28, Thu>

- ridgeplot for gseaResult <2017-08-17, Thu>

Changes in version 3.3.1:

- throw warning instead of error when fail to `setReadable`.
  <2017-06-28, Wed> +
  https://github.com/GuangchuangYu/clusterProfiler/issues/91

- better msg when using wrong ID types in GSEA <2017-06-01, Thu> +
  https://support.bioconductor.org/p/96512/#96516

[EBImage](https://bioconductor.org/packages/EBImage)
-------

Changes in version 4.20.0:

NEW FEATURES

- 'abind()' method for combining Image arrays

- 'floodFill()' has been vectorized over its arguments 'pt' and 'col'
  allowing to specify multiple start points and different fill colors

SIGNIFICANT USER-VISIBLE CHANGES

- 'display()' "browser" mode has been updated to use the htmlwidgets
  infrastructure

- 'filter2()' does not require filter dimensions to be odd numbers when
  filter size equals image size

BUG FIXES

- fixed issues with 'normalize()'
  (https://github.com/aoles/EBImage/issues/25)

- various improvements to the 'clahe()' function

[edgeR](https://bioconductor.org/packages/edgeR)
-----

Changes in version 3.20.0:

- DGEList() sets genes and counts to have same row.names.

- topTags() preserves row.names.

- estimateDisp() uses 'y$design' if it exists.

- estimateDisp() doesn't use average log-CPM in the prior.df
  calculation if 'trend.method' is 'none'.

- estimateDisp() doesn't return trended.dispersion if 'trend.method' is
  'none'.

- Design matrix defaults to 'y$design' before 'y$samples$group' in all
  the gene set testing functions.

- New arg 'group' for mglmOneWay(). Results in slight speed improvement
  for glmFit().

- 'design' arg for predFC() is now compulsory.

- Switched 'coef.start' back to a vector in mglmOneGroup().

- New functions cpmByGroup() and rpkmByGroup().

- Renamed arg 'x' to 'y' in cpm() and rpkm().

- Restored null dispersion check in glmFit().

- Removed 'offset' arg from glmQLFit() to be consistent with glmFit().

- Exported CompressedMatrix subset operator.

- Refactored C++ code with greater C++11 support to use Rcpp.

- Streamlined input dimension checks in C++ code.

- Supported zero-row input to addPriorCounts() C++ code.

- Added cbind and rbind S3 methods for DGEList objects.

- Added 'Dims' as part of the compressedMatrix class.

- Added common methods for the compressedMatrix class.

- Register S3 methods for compressedMatrix.

- Added a case study of differential methylation analysis to the user's
  guide.

[EGSEA](https://bioconductor.org/packages/EGSEA)
-----

Changes in version 1.5.6 (2017-09-11):

- modified: GSVA invokation due to changes from GSVA developers on the
  return value and arguments

Changes in version 1.5.5 (2017-08-24):

- fixed: bug in egsea.ora where interpret files not generated
  correctly.

- fixed: summary plot based on ranking when using egsea.ora

- Renamed: entrezIDs in egsea.ora into geneIDs

- fixed: bug in egsea.ora when 'title' includes white spaces or special
  characters

- fixed: several minor bugs

Changes in version 1.5.4 (2017-08-10):

- changed: entrezIDs into geneIDs in buildCustomIdx and buildGMTIdx.

Changes in version 1.5.3 (2017-07-20):

- Replaced: underscore characters and dots with dashes in vignette,
  gsaplots.R and htmlUtils.R

Changes in version 1.5.2 (2017-07-18):

- changed: method name from getlogFCFromLMFit to runStandardLimmaDEA

- fixed: a minor bugin buildKEGGIdx

- Added: votep to combining p-values

- Changed: egsea.dir into report.dir

- Replaced: print statements with message, warning, stop where
  appropriate.

- Added: interactive stats table and interactive summary plots.

- Added: parameter 'interactive' to egsea main functions,
  generateReport and plotSummary

- Modified: vignette so that images that are not generated are replaced
  with a warning message

Changes in version 1.5.1 (2017-04-11):

- Merged: the documentation of relevenat functions into one help page

- Added: a function buildGMTIdx to build gene set collection index from
  a GMT file

- Added: a new function named egsea.ma to work with Microarray dataset

- Modified: egsea functions to generate all the pairwise comparisons
  (contrasts) when the contrast parameter is NULL. This is mainly done
  based on the primary factor of the design matrix that is defined by
  'group' parameter or column in voom.results$targets$group.

[ENmix](https://bioconductor.org/packages/ENmix)
-----

Changes in version 1.13.9:

- getBeta: bugfix for handling large dataset

Changes in version 1.13.7:

- updated User's guide

Changes in version 1.13.6:

- corrected User's guide

Changes in version 1.13.4:

- bug fix

Changes in version 1.13.3:

- add pipeline function mpreprocess

- bug fix

Changes in version 1.12.1:

- bug fix

- bug fix

[EnrichedHeatmap](https://bioconductor.org/packages/EnrichedHeatmap)
---------------

Changes in version 1.7.1:

- show_row_names is put in the function argument list

- we don't use a position with 0 (either +1 or -1)

- add `dist_by_closeness()`

- add `extract_anno_enriched()`

- add ticks on axes

- add a new vignette (compare row ordering methods)

[EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser)
-----------------

Changes in version 2.8.0:

- Major migration from ExpressionSet to SummarizedExperiment

[ensembldb](https://bioconductor.org/packages/ensembldb)
---------

Changes in version 2.1.12:

BUG FIXES

- Use new defaults from the IRanges package for arguments maxgap = -1L,
  minoverlap = 0L in transcriptsByOverlaps and exonsByOverlaps methods.

- Remove RSQLite warnings (issue #54).

Changes in version 2.1.11:

BUG FIXES

- ensDbFromGtf failed to parse header for GTF files with more than one
  white space.

Changes in version 2.1.10:

USER VISIBLE CHANGES

- supportedFilters returns a data frame with the filter class name and
  corresponding field (column) name.

Changes in version 2.1.9:

NEW FEATURES

- Support for global filters in an EnsDb object.

- Add filter function.

Changes in version 2.1.8:

NEW FEATURES

- New annotations available in EnsDb objects: gene.description and
  tx.tx_support_level.

- New TxSupportLevelFilter object.

[EpiDISH](https://bioconductor.org/packages/EpiDISH)
-------

Changes in version 0.99.0:

- Initial submission to Bioconductor.

- Added NEWS file.

- Added constraint mode parameter in CP mode.

[epiNEM](https://bioconductor.org/packages/epiNEM)
------

Version: 2017.01
Category: Github made public and submitted to Bioconductor
Text:

[epivizrChart](https://bioconductor.org/packages/epivizrChart)
------------

Changes in version 0.99.0:

- Bioc Devel Release

[epivizrData](https://bioconductor.org/packages/epivizrData)
-----------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrData>

[epivizrServer](https://bioconductor.org/packages/epivizrServer)
-------------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrServer>

[epivizrStandalone](https://bioconductor.org/packages/epivizrStandalone)
-----------------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrStandalone>

[esATAC](https://bioconductor.org/packages/esATAC)
------

Version: 1.0.0
Category: INITIAL RELEASE
Text:

Version: 1.0.0
Category: Preset pipelines are available for case study and
        case-control study
Text:

Version: 1.0.0
Category: All basic elements in preset pipeline are available for
        building customized pipeline.
Text:

[ExperimentHub](https://bioconductor.org/packages/ExperimentHub)
-------------

Changes in version 1.4.0:

NEW FEATURES

- ExperimentHub will now work offline utilizing argument 'localHub';
  will also use this option automatically if no internet connection is
  detected.

- Add new vignette for creating ExperimentHub packages

MODIFICATIONS

- Update AnnotationHub dependency; new resource class RDSResource

- move listResources and loadResources from AnnotationHub to here

BUG FIXES

- Fix typo in createHubAccessors with hard coded value

[ExperimentHubData](https://bioconductor.org/packages/ExperimentHubData)
-----------------

Changes in version 1.4.0:

MODIFICATIONS

- Moved vignette to ExperimentHub as more practical use there; created
  new vignette

- Allow specification of metadata file in inst/extdata to be used

BUG FIXES

- updated addResources

[FamAgg](https://bioconductor.org/packages/FamAgg)
------

Changes in version 1.5.3:

SIGNIFICANT USER-VISIBLE CHANGES

- plotPed with haplopaint plotting supports device = "txt" that just
  writes the data.frame for haplopaint and returns the file name, does
  however not call Haplopaint.

Changes in version 1.5.2:

NEW FEATURES

- New binomial test implemented (binomialTest and FABinTestResults
  object).

[GA4GHclient](https://bioconductor.org/packages/GA4GHclient)
-----------

Changes in version 1.2:

BUG FIXES

- Fix bug due to changes in dplyr package

[gage](https://bioconductor.org/packages/gage)
----

Changes in version 2.27.3:

- fixed bug in kegg.gsets function, which cause error when species="ko"
  (KEGG Orthology).

Changes in version 2.27.1:

- major expansion in korg, which now include both KEGG and NCBI
  taxonomy IDs, two more gene ID types, i.e. NCBI protein and uniprot
  IDs. In addition, Entrez or NCBI Gene IDs are discontinued for most
  prokaryotes.

- korg now include 4800 KEGG species, in the meantime, an updated
  version of korg is now checked out from Pathview Web server each R
  session when kegg.gsets function is called the first time. version
  2.20.1

- updated RNA-seq workflow vignette, especially step 1
  summarizeOverlaps, and correct tophat web link. version 2.19.1

- updated gage main vignette and RNA-seq workflow vignette. Add
  reminder on species and gene ID data consistence check to the former,
  and correct Cufflinks output format and web link. version 2.17.2

- updated khier to included newly added reference pathways. kegg.gsets
  can work with 397 pathways now. version 2.15.5

- updated korg to included over 80 newly added species, such as sheep,
  apple, mandarin orange etc. Pathview can work with 3050 species now.
  version 2.14.3

- revised the internal function gs.heatmap as to allow margin sizes
  ajustment for gene set heatmaps directly or through sigGeneSet.
  version 2.14.1

- revised "RNA-Seq Data Pathway and Gene-set Analysis Workflows" to
  reflect summarizeOverlaps() migration to GenomicAlignments package
  for Bioc 2.14. version 2.13.5

- add function go.gsets, which generates up-to-date GO gene sets for 19
  common species annotated in Bioconductor and more others by the
  users. version 2.13.3 (2.12.3)

- updated korg to included over 600 newly added species. kegg.gsets can
  work with 2970 species now.

- fixed typos in joint workflows with Cufflinks (page 13):
  range(exp.fc) to range(cuff.fc) version 2.12.2

- fixed typos in joint workflows with DESeq2 (page 11): cnts.kegg.p
  should be fc.kegg.p

[gdsfmt](https://bioconductor.org/packages/gdsfmt)
------

Changes in version 1.14.0:

- tweak error message

[GeneNetworkBuilder](https://bioconductor.org/packages/GeneNetworkBuilder)
------------------

Changes in version 1.19.1:

- Fix the problem when STRINGdb not work.

[GENESIS](https://bioconductor.org/packages/GENESIS)
-------

Changes in version 2.7.4:

- In fitNullMM, for Binomial and Poisson GLM families, the variance of
  the fixed effects will not be multiplied by the dispersion parameter
  RSS.

Changes in version 2.7.3:

- Change defaults in assocTestSeq*: Davies method for SKAT p-values,
  flat weights.

[geneXtendeR](https://bioconductor.org/packages/geneXtendeR)
-----------

Changes in version 1.3.5:

- Remove deprecated README.

Changes in version 1.3.4:

- Update citation.

Changes in version 1.3.3:

- New peaksMerge() function.

Changes in version 1.3.2:

- New allPeakLengths() function.

Changes in version 1.3.1:

- New cumulative line plot (cumlinePlot()) function.

- Output to distinct() function is a lot more visually appealing now
  (no more tab delimiter signs), yay!

- New differential gene ontology (diffGO()) function.

- New gene-GO network (makeNetwork()) function.

- New GO word cloud (makeWordCloud()) function.

- New barplot function for plotting word frequencies found within GO
  terms (plotWordFreq()).

- New significant/total peaks (hotspotPlot()) function.

- Add new meanPeakLength(), meanPeakLengthPlot(), and
  peakLengthBoxplot() functions.

[GENIE3](https://bioconductor.org/packages/GENIE3)
------

Version: 0.99.4
Text: Removed argument "seed", use set.seed instead Minor improvemtents

Version: 0.99.2
Category: INTERFACE CHANGES
Text: Functions and parameters renamed to CamelCase rather than
        containing dots (e.g. get.link.list into getLinkList)

Version: 0.99.2
Category: NEW FEATURES
Text: GENIE3 now accepts ExpressionSet, SummarizedExperiment and SCESet
        as input

[genomation](https://bioconductor.org/packages/genomation)
----------

Version: 1.9.7
Category: IMPROVEMENTS AND BUG FIXES
Text: bug fix relating to calculation of the last bin in scoreMatrixBin
        function for bin.op = "max"
        (https://github.com/BIMSBbioinfo/genomation/issues/166)

Version: 1.9.7
Category: Implemented by Bozena Mika-Gospodorz
Text:

Version: 1.9.6
Category: IMPROVEMENTS AND BUG FIXES
Text: bug fix relating to calculation of the last bin in scoreMatrixBin
        function
        (https://github.com/BIMSBbioinfo/genomation/issues/166)

Version: 1.9.6
Category: Implemented by Bozena Mika-Gospodorz
Text:

Version: 1.9.5
Category: NEW FUNCTIONS AND FEATURES
Text: enrichmentMatrix() function computes an enrichment of IP sample
        over IgG or input DNA control sample (issue:
        https://github.com/BIMSBbioinfo/genomation/issues/138)

Version: 1.9.5
Category: IMPROVEMENTS AND BUG FIXES
Text: removed a requirement of having a 'chr' string in chromosome
        names in checkBedValidity() function (issue:
        https://github.com/BIMSBbioinfo/genomation/issues/160)

Version: 1.9.5
Category: Implemented by Bozena Mika-Gospodorz
Text:

Version: 1.9.4
Category: NEW FUNCTIONS AND FEATURES
Text: C++ functions that create a matrix storing data with desirable
        number of bins for each window: - listSliceMean() - calls the
        binMean() function, - listSliceMedian() - calls the
        binMedian(), - listSliceMax() - calls the binMax(), -
        listSliceMin() - calls the binMin(), - listSliceSum() - calls
        the binSum().

Version: 1.9.4
Category: NEW FUNCTIONS AND FEATURES
Text: C++ functions that compute the value for each bin: - binMean() -
        computes a mean value, - binMedian() - computes a median value,
        - binMax() - computes a maximum values, - binMin() - computes a
        minumum values, - binSum() - computes a sum value.

Version: 1.9.4
Category: NEW FUNCTIONS AND FEATURES
Text: C++ function Median_c() - computes a median value from a vector.

Version: 1.9.4
Category: Implemented by Bozena Mika-Gospodorz
Text:

Version: 1.9.3
Category: IMPROVEMENTS AND BUG FIXES
Text: tests for different parameter combinations for ScoreMatrixBin

Version: 1.9.2
Category: NEW FUNCTIONS AND FEATURES
Text: c() function to append a ScoreMatrix as well as a ScoreMatrixList
        to an existing ScoreMatrixList (issue:
        https://github.com/BIMSBbioinfo/genomation/issues/151).
        Implemented by Bozena Mika-Gospodorz.

Version: 1.9.2
Category: IMPROVEMENTS AND BUG FIXES
Text: added a slot "names" to a ScoreMatrixList class

Version: 1.9.1
Category: IMPROVEMENTS AND BUG FIXES
Text: The knitrBootstrap dependecy is removed and following issues are
        fixed - https://github.com/BIMSBbioinfo/genomation/issues/135 -
        https://github.com/BIMSBbioinfo/genomation/issues/146 -
        https://github.com/BIMSBbioinfo/genomation/issues/147

[GenomeInfoDb](https://bioconductor.org/packages/GenomeInfoDb)
------------

Version: 1.14.0
Category: NEW FEATURES
Text:

Version: 1.14.0
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text:

Version: 1.14.0
Category: DEPRECATED AND DEFUNCT
Text: Remove 'force' argument from seqinfo() and seqlevels() setters
        (the argument got deprecated in BioC 3.5 in favor of new and
        more flexible 'pruning.mode' argument).

Version: 1.14.0
Category: BUG FIXES
Text: Add missing Y/chrY entry in seqlevel style db for Drosophila
        melanogaster and Rattus norvegicus.

[GenomicAlignments](https://bioconductor.org/packages/GenomicAlignments)
-----------------

Changes in version 1.14.0:

SIGNIFICANT USER-LEVEL CHANGES

- makeGAlignmentPairs() no more drops pairs with discordant seqnames.

- Change 'maxgap' and 'minoverlap' argument defaults in methods of the
  findOverlaps() so they adhere to the new argument defaults of the
  generic defined in IRanges 2.12.0. See NEWS file in the IRanges
  package for more information about this change.

DEPRECATED AND DEFUNCT

- Remove 'force' argument from seqinfo() and seqlevels() setters (the
  argument got deprecated in BioC 3.5 in favor of new and more flexible
  'pruning.mode' argument).

BUG FIXES

- Fix bug in pairing code of readGAlignmentPairs() when one mate in a
  pair is lost because of user-supplied filtering (e.g. mapqFilter=10).

[GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures)
---------------

Changes in version 1.30:

NEW FEATURES

- Add makeTxDbFromEnsembl() for creating a TxDb object by querying
  directly an Ensembl MySQL server. This seems to be faster and more
  reliable than makeTxDbFromBiomart().

- Improve makeTxDbFromBiomart() support for EnsemblGenomes marts
  fungal_mart, metazoa_mart, plants_mart, and protist_mart.

- makeTxDbFromGFF() and makeTxDbFromGRanges() now import the CDS phase.
  This required a change in the schema of the underlying SQLite db of
  TxDb objects. This is still a work-in-progress e.g. cdsBy(txdb,
  by="tx") still needs to be modified to return the phase info.

SIGNIFICANT USER-VISIBLE CHANGES

- The *ByOverlaps() functions now use the same 'maxgap' and
  'minoverlap' defaults as subsetByOverlaps().

DEPRECATED AND DEFUNCT

- Remove 'force' argument from seqinfo() and seqlevels() setters (the
  argument got deprecated in BioC 3.5 in favor of new and more flexible
  'pruning.mode' argument).

BUG FIXES

- exonicParts() and intronicParts() are now documented.

- Address a couple of issues pointed out by Matt Chambers in internal
  helpers get_organism_from_Ensembl_Mart_dataset() and
  .extractEnsemblReleaseFromDbVersion() used by makeTxDbFromBiomart().

- Fix internal utility .Ensembl_getMySQLCoreDir(). Was failing for some
  of the 69 datasets from the Ensembl mart, causing
  makeTxDbFromBiomart() to fail loopking up the organism scientific
  name and the chromosome lengths. Thanks to Matt Chambers for
  reporting this.

- Some tweaks and fixes needed to support RSQLite 2.0.

[GenomicRanges](https://bioconductor.org/packages/GenomicRanges)
-------------

Changes in version 1.30.0:

NEW FEATURES

- Support GPos-based GRangesList objects.

- Add 'na.rm' argument to binnedAverage().

SIGNIFICANT USER-LEVEL CHANGES

- Change 'maxgap' and 'minoverlap' defaults for findOverlaps() and
  family (i.e. countOverlaps(), overlapsAny(), and subsetByOverlaps()).
  This change addresses 2 long-standing issues: (1) by default
  zero-width ranges are not excluded anymore, and (2) control of
  zero-width ranges and adjacent ranges is finally decoupled (only
  partially though). New default for 'minoverlap' is 0 instead of 1.
  New default for 'maxgap' is -1 instead of 0. See ?findOverlaps for
  more information about 'maxgap' and the meaning of -1. For example,
  if 'type' is "any", you need to set 'maxgap' to 0 if you want
  adjacent ranges to be considered as overlapping.

- GPos now extends GRanges but with a ranges slot that must be an IPos
  object. Update "old" GPos objects with updateObject().

- Move pos() generic to IRanges package.

- Move rglist() generic to IRanges package.

- Rename GenomicRangesORmissing and GenomicRangesORGRangesList classes
  -> GenomicRanges_OR_missing and GenomicRanges_OR_GRangesList,
  respectively.

- Remove "seqinfo" method for RangesList objects.

- Remove "stack" method for GenomicRangesList objects.

DEPRECATED AND DEFUNCT

- Remove 'force' argument from seqinfo() and seqlevels() setters (the
  argument got deprecated in BioC 3.5 in favor of new and more flexible
  'pruning.mode' argument).

BUG FIXES

- nearest() and distanceToNearest() now call findOverlaps() internally
  with maxgap=0 and minoverlap=0. This fixes incorrect results obtained
  in some situations e.g. in the situation reported here:
  https://support.bioconductor.org/p/99369/ (zero-width ranges) but
  also in this situation: nearest(GRanges("chr1", IRanges(5, 10)),
  GRanges("chr1", IRanges(1, 4:5)), select="all") where the 2 ranges in
  the subject are *both* nearest to the 5-10 range.

- '$' completion on GenomicRanges works in RStudio.

- Minor tweaks to conversion from character to GRanges and reverse
  conversion.

[GenomicScores](https://bioconductor.org/packages/GenomicScores)
-------------

Changes in version 1.2.0:

USER VISIBLE CHANGES

- Added methods 'name()' and 'type()' for GScores objects.

- Enabled the retrieval of multiple score values per genomic position
  (e.g., as in CADD or M-CAP scores).

- Added method 'citation()' to fetch citation information for genomic
  scores.

- Added function 'makeGScoresPackage()' to create an annotation package
  from an AnnotationHub genomic scores resource.

- Added 'qfun()' and 'dqfun()' methods to fetch the quantization and
  dequantization functions from used to store and retrieved genomic
  scores.

- Added 'quantized' argument to the 'scores()' method to obtain
  quantized values if the user wants to dequantize the values him or
  herself.

- Fallback to local AnnotationHub when there is no internet connection
  to fetch genomic scores through AnnotationHub resources.

- Added 'MafDb' class, derived from 'GScores' to store and access minor
  allele frequency values. This was originally defined in the
  'VariantFiltering' package.

- The vignette has been updated to illustrate the use of some of the
  previous changes.

[GEOquery](https://bioconductor.org/packages/GEOquery)
--------

Version: 2.45.2
Text: Improvements: * GPL parsing 4-5x faster * GSM parsing 3x faster *
        GSEMatrix parsing is much smarter with respect to sample
        characteristics. In short, for GSEs where sample
        characteristics are actually used, the pData should have nice,
        neat column headers with the phenodata keys and values in the
        columns, including correct handling of missing values, etc.

Version: 2.45.1
Text: Bug fixes * getDirectoryListing fixed to deal with changes to
        NCBI server listing formats

Version: 2.45
Text: Improvements: * GDS parsing is 2-3x faster * GSEMatrix parsing is
        2-3x faster

[ggcyto](https://bioconductor.org/packages/ggcyto)
------

Version: 1.5.5
Category: add data replacement feature through %+% operator
Text:

[ggtree](https://bioconductor.org/packages/ggtree)
------

Changes in version 1.9.4:

- geom_hilight_encircle and geom_cladelabel2 <2017-09-12, Tue> +
  https://github.com/GuangchuangYu/ggtree/pull/152

- set_hilight_legend <2017-08-30, Wed>

- geom_motif for aligned motif <2017-08-22, Tue> +
  https://github.com/GuangchuangYu/ggtree/issues/148

- fixed `R CMD build` error: cannot stat 'ggtree/site_src/themes': No
  such file or directory <2017-08-22, Tue>

Changes in version 1.9.3:

- update to using !! in tidyr::gather for compatible with tidyr 0.7.0
  <2017-08-03, Thu>

- now geom_text2, geom_label2, geom_point2 and geom_segment2 work with
  ggplot2 <2017-08-01, Tue>

- update fortify.jplace to support number of placement (nplace)
  <2017-07-27, Thu>

Changes in version 1.9.2:

- add bg_line and height parameter in msaplot <2017-07-26, Wed> + use
  can set bg_line = FALSE and height = 1 to plot more beautiful
  alignment

- extend parameter in geom_cladebar <2017-07-26, Wed> +
  https://github.com/GuangchuangYu/ggtree/issues/142#issuecomment-317817995

- scaleClade works after calling viewClade <2017-07-20, Thu> +
  https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!topic/bioc-ggtree/QVSryszPaFY

- gheatmap support handling collapsed tree <2017-06-29, Thu> +
  https://github.com/GuangchuangYu/ggtree/issues/137

Changes in version 1.9.1:

- now mapping parameter will passed to segment layer in
  geom_tiplab(align=T) <2017-06-19, Mon>

- geom_cladelabel support `angle="auto"` for circular layout tree
  <2017-05-05, Fri>

[Glimma](https://bioconductor.org/packages/Glimma)
------

Changes in version 1.6.0:

- Added table to MDS plot.

- Changed encoding of javascript data to be more compact

- Fixed handling of top and gene.selection parameters in glMDSPlot

[GoogleGenomics](https://bioconductor.org/packages/GoogleGenomics)
--------------

Changes in version 2.0.0:

gRPC support and more authentication options

- Can use gRPC to access the entire method set in v1 API.

- Support for application default credentials.

- Support for GCE service account credentials.

[GOSemSim](https://bioconductor.org/packages/GOSemSim)
--------

Changes in version 2.3.1:

- new project site using blogdown <2017-09-28, Thu>

- speed up by pre-calculating GO similarities <2017-05-22, Mon> +
  contributed by Lluís Revilla Sancho +
  https://github.com/GuangchuangYu/GOSemSim/pull/13

[graphite](https://bioconductor.org/packages/graphite)
--------

Changes in version 1.23.7 (2017-10-19):

- Provide the URL of each pathway in its original database.

- New vignette describing pathway analysis of metabolic activities.

Changes in version 1.23.6 (2017-10-18):

- Removed DEGraph support.

Changes in version 1.23.5 (2017-10-13):

- runClipper gains an option to set the seed for random number
  generation.

Changes in version 1.23.3 (2017-10-12):

- Metabolites in pathways.

- Multiple views of each pathway: proteins only, metabolites only or
  both.

- New databases: SMPDB and PharmGKB.

- Removed support for RCytoscape (by default, use RCy3).

- Conversion of metabolite identifiers.

- Faster conversion of idenfiers on multicore systems.

[GraphPAC](https://bioconductor.org/packages/GraphPAC)
--------

Changes in version 1.19.1:

- Url corrections for new pdb.org website file structure.

[GRmetrics](https://bioconductor.org/packages/GRmetrics)
---------

Changes in version 1.3.3:

- Renamed column "GR" to "GRvalue" for GR value table to match other
  code.

- Fixed erroneous error message in division rate GR value calculation

Changes in version 1.3.2:

- Renamed some curve parameters for the relative cell count curve

- changed "IC" to "rel_cell".

- Changed GRdrawDRC function parameter option for dose response curves
  based on relative cell counts from 'metric = "IC"' to 'metric =
  "rel_cell"'.

- Changed GRinf (and Einf) value to the minimum of the mean GR values
  (relative cell count) at the two highest concentrations tested for
  the case of horizontal line fits (same as GRmax).

- Changed GRmax (and Emax) so that in the case of averaging over
  multiple conditions, the value taken is the minimum of the mean of GR
  values (relative cell count) at the two highest concentrations
  instead of the absolute minimum GR value at these concentrations.

- "duration" column in input changed to "treatment_duration" to match
  the nomenclature in other code.

- Added calculation of GR values from division rates for Case "C"

Changes in version 1.3.1:

- Added calculation of GR values from division rates for Case "A"

- now accepts columns "duration" and "division_time" instead of
  "cell_count__time0"

- Added a few parameters to calculation

- "control_cell_doublings" is the number of cell doublings that occur
  in the control population over the assay, calculated either from
  initial and final cell counts or given division rate and time of
  assay.

- "concentration_points" is the number of different concentrations used
  in an experiment

[GSEABase](https://bioconductor.org/packages/GSEABase)
--------

Changes in version 1.39:

BUG FIXES

- goSlim() did not correctly count duplicate identifiers.
  (https://support.bioconductor.org/p/100403/)

[GSVA](https://bioconductor.org/packages/GSVA)
----

Changes in version 1.26:

USER VISIBLE CHANGES

- Updated implementation of the option 'abs.ranking=TRUE' to use the
  original Kuiper statistic.

- Arguments 'rnaseq' and 'kernel' have been deprecated and replaced by
  a new argument 'kcdf'.

- Arguments 'no.bootstraps' and 'bootstrap.percent' have been
  deprecated.

- The return value with the default argument 'method="gsva"' has been
  simplified and it is not a list object anymore. Now the 'gsva()'
  function return always a matrix or an 'ExpressionSet' object, when
  the input expression data is also an 'ExpressionSet' object.

- The 'gsva()' function can now be used through a shiny app that runs
  through the function 'igsva()'.

[gtrellis](https://bioconductor.org/packages/gtrellis)
--------

Changes in version 1.9.1:

- Functions in GenomicRanges are imported

[GWASTools](https://bioconductor.org/packages/GWASTools)
---------

Changes in version 1.23.2:

- Add a function to coerce a GenotypeData object to a VCF object for
  use with VariantAnnotation.

Changes in version 1.23.1:

- Move ncdf4 from Imports to Suggests. Users who wish to use NetCDF
  files instead of GDS will have to install the ncdf4 package
  separately. This change eliminates the requirement to install the
  NetCDF library on Linux machines for users who plan to use GDS only.

[Heatplus](https://bioconductor.org/packages/Heatplus)
--------

Changes in version 2.23.1:

- Resolved issue with inconsistent plotting on the same device

[hiAnnotator](https://bioconductor.org/packages/hiAnnotator)
-----------

Changes in version 1.11.1:

- GenomicRanges package update adjustments & improved tests

[HIBAG](https://bioconductor.org/packages/HIBAG)
-----

Changes in version 1.14.0:

- modify the kernel to support the GPU extension

- add matching proportion to measure the similarity of SNP haplotypes
  between training and test datasets

- new function `hlaReportPlot()`

- the argument 'cl' in `predict.hlaAttrBagClass()`, `hlaPredict()` and
  `hlaParallelAttrBagging()` allows a numeric value for the number of
  cores

[HiTC](https://bioconductor.org/packages/HiTC)
----

Changes in version 1.21.0:

BUG FIXES

- Fix bug in getExpectedCountsMean for non-symmetrical data

- Deal with NA in getPearson function

- Fix bug in normLGF leading to non symmetrical matrices

[hpar](https://bioconductor.org/packages/hpar)
----

Changes in version 1.19.1:

- Re-generate data sets from scripts/getHpaData.R, as discrepancies
  with the data downloaded form the hpa site were documented by Martin
  Bush <2017-07-19 Wed>

Changes in version 1.19.0:

- new Bioconductor devel

[iCOBRA](https://bioconductor.org/packages/iCOBRA)
------

Changes in version 1.5.5:

- The default selection of the measure to use for ROC and FP curves
  have been changed to improve consistency. Now, the preferred order is
  pval, padj, score. The behaviour of previous versions can be obtained
  by setting prefer_pval = FALSE in calculate_performance(). Note that
  the type of measure that is used to calculate a certain performance
  value can be retrieved from the respective slots of the
  COBRAPerformance and COBRAPlot objects.

- Improved robustness in selection of measure to use for FDR/TPR and
  FDR/NBR curves, especially for methods where both p-values and scores
  are available.

- Additional validity checks for pval and padj slots

Changes in version 1.5.1:

- Added support for false sign rate calculations

[immunoClust](https://bioconductor.org/packages/immunoClust)
-----------

Version: 1.9.2
Text: CHANGES * Optional label parameter in meta-clustering to continue
        the meta-clustering with an initial cluster to meta-cluster
        assignment

Version: 1.9.1
Text: CHANGES * Minor improvements and additional option in plot and
        splom methods * Introduce of ALPHA option also for
        normalization precedures

[intansv](https://bioconductor.org/packages/intansv)
-------

Changes in version 1.17.1:

Notes

- Fix a small bug in function methodsMerge.

Changes in version 1.17.0:

Notes

- Update to the latest version of different tools.

[InteractionSet](https://bioconductor.org/packages/InteractionSet)
--------------

Changes in version 1.5.7:

- Removed c() method for InteractionSet, rbind method for
  GInteractions.

- Generalized ContactMatrix to allow any type of matrix-like object.

- Supported inflate() for GInteractions without specifying fill.
  Changed default fill for InteractionSet to missing.

- Separated subsetting and combining documentation into two different
  pages.

- Added convenience wrappers for resize(), narrow() and shift() on the
  GenomicRanges slot in all objects.

- Modified anchors() to return a list rather than GRangesList.

- Modified width() for GInteractions to return a list() rather than a
  DataFrame(). Also changed names of list element.

- Added the anchorIds() function for rapid extraction of anchor IDs.

- Removed the requirement for identical regions in pcompare() and
  match().

[IntEREst](https://bioconductor.org/packages/IntEREst)
--------

Changes in version 1.2.0:

NEW FEATURES

- DEXSeqIntEREst() runs DEXSeq differential exon usage or intron
  retention test on the SummarizedExperiment objects (resulted from
  interest() or interest.sequential() analysis).

- unionRefTr() performs union on the genes in reference data frame with
  overlapping introns/exons. The resulted data frame features from each
  repeating exon or intron a single copy only.

- deseqInterest() runs differential intron retention test adapted from
  the DESeq2 package.

- interestResultIntEx() builds a SummarizedExperiment object from an
  intron retention and an exon-exon junction object (resulted from
  interest(), interest.sequential() and/or InterestResult() functions).
  The average of the junction levels of the flanking exons are added to
  the SummerizedExperiment with the intron retention values.

SIGNIFICANT USER-VISIBLE CHANGES

- u12Index() supports intronExon parameter that provides the
  possibility to extract either the rows that represent the U12-type
  introns or the exons flanking to the U12-type introns from a
  SummerizedExperiment object.

BUG FIXES

- Correcting interestAnalyse.R and interestAnalyse.sequential.R so that
  interest() and interest.sequential() support bam files with single
  (unpaired) reads.

- Correct bpparam setting in interest().

[iPAC](https://bioconductor.org/packages/iPAC)
----

Changes in version 1.20.1:

- Url corrections for new pdb.org website file structure.

[IPO](https://bioconductor.org/packages/IPO)
---

Changes in version 1.7.5:

- added usage of clustering method FORK on unix-systems (thanks to
  Pablo Moreno)

- fixed bug in parallelization to prevent conflicts with package 'snow'

Changes in version 1.7.4:

- preceded parallel-functions with 'parallel::' to use right package

- fixed bug in function writeRScript using 'loess' retention time cor.

- decreased runtime for R CMD check IPO

Changes in version 1.7.3:

- added runnable examples

- decreased size of pictures in vignettes/rsmDirectory

- decreased runtime for unit-tests

- replaces expand.grid with expand.grid.subset (in utils.R)

Changes in version 1.7.2:

- bugfix: try to prevent error in calcPPS possibly caused by NAs

- replaced cat() and print() calls with message()

Changes in version 1.7.1:

- checking correlation of peak-shape with sinus curve (-pi/2 to
  pi*1.5), normal distribution or checkBorderIntensity

- findIsotopes.IPO renamed parameter checkBorderIntensity to
  checkPeakShape

- performance improvement calcPPS for checkPeakShape=FALSE

- calculating xcmsSet-object and respective PPS for each DoE. (PPS is
  not estimated from rsm anymore)

- additionally forwarding nSlaves for xcmsSet-function (also see
  getDefaultXcmsSetStartingParams())

Changes in version 1.7.0:

- added support for XCMS-method retcor.loess

- updated help files

- changed return value of getRGTVValues

- adapted unit tests

- parameter scanrange for XCMS-methods findPeaks can be set but not
  optimized

Changes in version 1.6.2:

- Updated the function getNormalizedResponse() to prevent NAs

Changes in version 1.6.1:

- Added installation script and installation description in vignette

Changes in version 1.6:

- Added support of CAMERA isotope identification (findIsotopes.CAMERA)

- selectivity of findIsotopes.IPO may be increased if
  checkBorderIntensity is set to TRUE: 'maxo' of each peak of an
  isotopologue must be three times higher than the intensities at
  'rtmin' and 'rtmax'

- simplified return value of calcPPS() to vector with meaningful names

- changed getDefaultXcmsSetStartingParams() for min_peakdwidth = c(12,
  28) and for ppm to c(17, 32)

Changes in version 1.5.7.1:

- using predict() to identify best levels and expand.grid to generate
  testdata for model

Changes in version 1.5.7:

- supporting single parameter optimization. Only basic version with
  redundant levels in consecutive DoEs

- removed integer-rounding in maximum focusing for all findPeaks
  parameters except prefilter(I) and steps

- Update documentation to match code

- Remove use of getwd() preventing absolute subdir paths

Changes in version 1.5.6:

- updated vignette

- generally using Central-Composite design instead of Box-Behnken
  design

- fixed bug when defining subdir=NULL in functions optimizeXcmsSet and
  optimizeRetGroup

- modified unit tests to handle versions > 1.5.6

- added function writeRScript to NAMESPACE export

- updated man for optimizeXcmsSet and optimizeRetGroup

Changes in version 1.5.5:

- increased recall of reliable peaks in calcPPS. Changed exponent for
  PPS calculation from 1.5 to 2

Changes in version 1.5.4.8:

- remove file lookup in optimizeXcmsSet, and leave that to xcmsSet()

Changes in version 1.5.4.7:

- explicitely use serial evaluation if nSlaves=1, to help debugging

Changes in version 1.5.4.6:

- remove dependency on Rmpi

Changes in version 1.5.4.5:

- added examples from msdata

- fix Depends, imports and library() and require() calls

Changes in version 1.5.4.4:

USER VISIBLE CHANGES

- packaged script

- changed method name attachparams to attachList

- changed method name calculateRGTV to getRGTVValues

- changed method name getDefaultStartingXcmsParams to
  getDefaultXcmsSetStartingParams

- changed method name typeCastFactor to typeCastParams

- changed method name writeRSkript to writeRScript

- changed the parameter name n_slaves to nSlaves

- resultIncreased: if last optimization score was 0, no isotopes have
  been found hence the dataset is not optimizable with IPO.

- added man files for attachList, calcPPS, combineParams, createModel,
  decode, decodeAll, encode, getBbdParameter, getCcdParameter,
  getDefaultRetCorCenterSample, getDefaultRetGroupStartingParams,
  getDefaultXcmsSetStartingParams, getNormalizedResponse,
  getRGTVValues, IPO-package, optimizeRetGroup, optimizeXcmsSet,
  startSlaves, toMatrix, typeCastParams

- removed xcmsSetsettingsAsString.R

- getResponses: now able to handle NULL value for slices parameter

- calcPPS: peaks with NA values are removed before isotopes
  identification

- if subdir is NULL, no rsm's are saved

- optimizeXcmsSet: lowere minimum value for min_peakwith from 5 to 3
  #IPO_V1.5.4.3: * LIP calculation in calcPPS fixed #IPO_V1.5.4.2: *
  added initial parameter check # * renamed all factor-variables to
  params # * in optimizeXcmsSet: - also look for mzML-files # - check
  if files were found # * bug in optimization for matchedFilter fixed;
  sigma and mzdiff have to be # definded later (combineFactors()) when
  sigma and step as well as steps are already known #IPO_V1.5.4.1:
  changes in calcPPS: # rt_window <- rt * 0.005 # rt_lower <-
  part_peaks[,"rt"] - rt_window # rt_upper <- part_peaks[,"rt"] +
  rt_window #IPO_V1.5.4: if bad_group == 0; bad_group = 1 && good_group
  += 1 #IPO_V1.5.3: no parameter for isotope detection.  #
  c13_peak[,"mz"] has to be within (mzmin + isotope_mass) and (mzmax +
  isotope_mass) # c13_peak[,"rt"] has to be within (rtmin +
  isotope_mass) and (rtmax + isotope_mass) #IPO_V1.5.: in
  RCSandGSIncreased: also used good_groups ^ 2 #IPO_V1.4.: vectorized
  isotope identification; # no intensity window, between intensity of
  max carbon and 1 #IPO_V1.3.: good_groups ^ 2 to increase recall

Changes in version 1.3.3:

User visible

- new `plot`-parameter for `optimizeXcmsSet` and `optimizeRetGroup` to
  control plotting (#51)

Changes in version 1.3.2:

- bug fix #50: correct peaks-matrix, to handle xcms bug
  (sneumann/xcms#220) for older xcms-versions

- test order of parameters to optimize

Changes in version 1.3.1:

- bug fix regarding conflict of BPPARAM and nSlaves arguments (thx to
  @lauzikaite)

[IRanges](https://bioconductor.org/packages/IRanges)
-------

Changes in version 2.12.0:

NEW FEATURES

- Add IPos objects for storing a set of integer positions where most of
  the positions are typically (but not necessarily) adjacent.

- Add coercion of a character vector or factor representing ranges
  (e.g.  "22-155") to an IRanges object, as well as "as.character" and
  "as.factor" methods for Ranges objects.

- Introduce overlapsRanges() as a replacement for "ranges" methods for
  Hits and HitsList objects, and deprecate the latter.

- Add "is.unsorted" method for Ranges objects.

- Add "ranges" method for Ranges objects (downgrade the object to an
  IRanges instance and drop its metadata columns).

- Add 'use.names' and 'use.mcols' args to ranges() generic.

SIGNIFICANT USER-VISIBLE CHANGES

- Change 'maxgap' and 'minoverlap' defaults for findOverlaps() and
  family (i.e. countOverlaps(), overlapsAny(), and subsetByOverlaps()).
  This change addresses 2 long-standing issues: (1) by default
  zero-width ranges are not excluded anymore, and (2) control of
  zero-width ranges and adjacent ranges is finally decoupled (only
  partially though).  New default for 'minoverlap' is 0 instead of 1.
  New default for 'maxgap' is -1 instead of 0. See ?findOverlaps for
  more information about 'maxgap' and the meaning of -1. For example,
  if 'type' is "any", you need to set 'maxgap' to 0 if you want
  adjacent ranges to be considered as overlapping.  Note that
  poverlaps() still uses the old 'maxgap' and 'minoverlap' defaults.

- subsetByOverlaps() first 2 arguments are now named 'x' and 'ranges'
  (instead of 'query' and 'subject') for consistency with the
  transcriptsByOverlaps(), exonsByOverlaps(), and cdsByOverlaps()
  functions from the GenomicFeatures package and with the
  snpsByOverlaps() function from the BSgenome package.

- Replace ifelse() generic and methods with ifelse2() (eager
  semantics).

- Coercion from Ranges to IRanges now propagates the metadata columns.

- Move rglist() generic from GenomicRanges to IRanges package.

- The "union", "intersect", and "setdiff" methods for Ranges objects
  don't act like endomorphisms anymore: now they always return an
  IRanges *instance* whatever Ranges derivatives are passed to them
  (e.g. NCList or NormalIRanges).

DEPRECATED AND DEFUNCT

- Deprecate "ranges" methods for Hits and HitsList objects (replaced
  with overlapsRanges()).

- Deprecate the "overlapsAny", "subsetByOverlaps", "coverage" and
  "range" methods for RangedData objects.

- Deprecate the universe() getter and setter as well as the 'universe'
  argument of the RangesList(), IRangesList(), RleViewsList(), and
  RangedData() constructor functions.

- Default "togroup" method is now defunct (was deprecated in BioC 3.3).

- Remove grouplength() (was deprecated in BioC 3.3 and replaced with
  grouplengths, then defunct in BioC 3.4).

BUG FIXES

- nearest() and distanceToNearest() now call findOverlaps() internally
  with maxgap=0 and minoverlap=0. This fixes incorrect results obtained
  in some situations e.g. in the situation reported here:
  https://support.bioconductor.org/p/99369/ (zero-width ranges) but
  also in this situation: nearest(IRanges(5, 10), IRanges(1, 4:5),
  select="all") where the 2 ranges in the subject are *both* nearest to
  the 5-10 range.

- Fix restrict() and reverse() on IRanges objects with metadata
  columns.

- Fix table() on Ranges objects.

- Various other minor fixes.

[IrisSpatialFeatures](https://bioconductor.org/packages/IrisSpatialFeatures)
-------------------

Changes in version 0.99.3 (2017-08-10):

- First version

[IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR)
---------------------

Version: 2017-10-25
Text: Fixed a small mistake in the documentation causing build warnings

Version: 2017-10-22
Text: isoformSwitchTestDRIMSeq() was updated to per default use
        dmFilter()

Version: 2017-10-22
Text: Small updates to documentation better explaining the
        functionalities from udate 0.99.12

Version: 2017-10-19
Text: Version bump for Bioconductor to keep up

Version: 2017-10-12
Text: importIsoformExpression() have been completely redesigned to
        utilize the tximport package as well as implementing the option
        for inter-library normalization of abundance (TxPM) values.

Version: 2017-10-12
Text: The vignette got a thorough workover - huge shoutout to Maria
        Dalby for the help!

Version: 2017-10-12
Text: isoformSwitchTestDRIMSeq() was extended to also include the
        dmFilter() functionality as part of the workflow.

Version: 2017-10-12
Text: The internal process calculating gene expression from isoform
        expression was cast as its own function: isoformToGeneExp().

Version: 2017-10-12
Text: Fixed an error that could cause problems when importing CDSs from
        a GTF file

Version: 2017-10-12
Text: Updated descriptions and other minor style issues.

Version: 2017-06-01
Text: Fixes some issue raised in the Bioconductor review To adhere to
        Bioconductor conventions the subset() method was removed and
        replaced by the subsetSwitchAnalyzeRlist() function.

Version: 2017-06-01
Text: The importIsoformExpression() function was updated to support
        import of Transcript Per Million (TxPM) as the relative
        abundance measure (Instead of TPM and RPKM/FPKM, which are
        discontinued) when importing data from Kallisto, Salmon and
        RSEM.

Version: 2017-06-01
Text: The isoformSwitchTestDRIMSeq() function was updated to make one
        linear model (one dmFit) instead of one model per pairwise
        comparison.

Version: 2017-06-01
Text: Small update to the switchPlot() functions to make it robust to
        NA annotation in non-essential data.

Version: 2017-06-01
Text: Added citation information since the article describing the R
        package was published: Vitting-Seerup et al. The Landscape of
        Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).

Version: 2017-05-24
Text: Fixes some issue raised in the Bioconductor review

Version: 2017-05-24
Text: Fixes a but introdued during the recent update in how pfam
        results were integrated.

Version: 2017-05-24
Text: Updates of the vignette for inproved readability.

Version: 2017-05-19
Text: MAJOR update which: 1) Introduces the iso_ref and gene_ref
        handles to all entires in the switchAnalyzeRlist which allows
        for easy integration of data across the different enteries.  2)
        Now offers full integration with the DRIMSeq tool which
        utilizises advanced linear models to identify significant
        changes in isoform usage at isoform level enabling robust
        analysis of more complex designs including batch effects. The
        integraiton is availabe via the isoformSwitchTestDRIMSeq()
        function.  3) Updates IsoformSwitchAnalyzeR to handle EBI's new
        server for running Pfam. 4) To enable the integration with
        DRIMSeq switchAnalyzeRlist object have been extended with: a)
        Isoform replciate count matrix. b) A design matrix.  5) The
        preFilter function have been updated with new functionalities
        and default cutoffs that are more suitable for use with
        DRIMSeq. See function documentation for details.  6) Implements
        suggested updates from Bioconductor reviewer This update is so
        large backward compatability is unfortunatly not feasiblie so
        all existing switchAnalyzeRlists will have to be remade. The
        extention of the switchAnalyzeRlist have also made a few
        changes in how to import data nessesary. Specifically: - The
        importRdata() function now take a replicate count matrix as
        it's main input and the replicate FPKM matrix is optional.  -
        The importBallgownData() function and it's accompanying
        "exampleRdata.RData" have been decapitated since it does not
        contain count information.  - The importIsoformExpression()
        function have been introduced to help with importing data from
        Kallisto, Salmon and RSEM. This function generates a isoform
        count matrix from the parent directory of the
        Kallisto/Salmon/RSEM analysis - which can easily be used with
        the importRdata() function to generate a switchAnalyzeRlist.
        Lastly the vignette have naturally been updated and improved
        accordingly.

Version: 2017-04
Text: Small incremental updates to ensure IsoformSwitchAnalyzeR lives
        up to all Bioconductor standards mostly consering how
        namespaces are organised and imported.

Version: 2017-04-18
Text: The following functionalities were added: - Enable filtering for
        significant switches in the preFilter() function.  - The
        extractGenomeWideAnalysis() function was extended with the
        "annotationToAnalyze" parameter enabling specification of which
        annotation types to analyze.  - The analyzeSwitchConsequences()
        function was extended to enable analysis of truncated protein
        (by supplying 'domain_length' to the 'consequencesToAnalyze'
        argument). - The analyzeSwitchConsequences() function was
        extended so the 'ntCutoff' also applies to TSS and TTS
        analysis. The following bugs were corrected: - A bug where
        importCufflinksCummeRbund() imported all genomic features of
        isoforms, including CDS etc, resulting in duplicated regions
        which caused problems for the intron retention analysis. This
        is only a problem for Cufflinks/Cuffdiff analysis where the
        refrence transcriptome contaied non-exon annotation (as defined
        in the type columns (column 3)) of the gtf file.  - A bug in
        the analyzePFAM() function that sometimes prevented
        IsoformSwitchAnalyzeR in correctly format the result file
        whereby the function could not run.  - The multi-threading
        option was removed since it was not supported by windows
        computers. We plan to bring it back in a later update.  - The
        option of manually supplying the start and stop codon sequences
        that the annotateORF() function should scan for in transcripts.
        Furthermorethe vignette was extended for enhanced usability.

[isomiRs](https://bioconductor.org/packages/isomiRs)
-------

Changes in version 1.5.5:

FIX

- Migrate vignette to new BiocStyle

- Remove unused function join_all

- Use parameter not integer number

- Using testthat for unit test

Changes in version 1.5.4:

FEATURES

- Better documentation for isoCorrection function. Add proper authors
  and citation.

Changes in version 1.5.3:

FEATURES

- Better colors for polar plot of isomiRs

FIX

- Fix notes during R CHECK for variables inside dplyr/ggplot functions

- Use roxygen2 for NAMESPACE

Changes in version 1.5.2:

FIX

- Remove TMB dependency

Changes in version 1.5.1:

FEATURES

- Add new polar figure to plot all isomiRs at the same time

- Add NLQO distribution to correct expression knowing sequencing bias
  [Argyropoulos et al, 2017]

- Improve data documentation

[IWTomics](https://bioconductor.org/packages/IWTomics)
--------

Version: 1.1.1
Text:

Version: 1.1.2
Text:

Version: 1.1.3
Text:

[JASPAR2018](https://bioconductor.org/packages/JASPAR2018)
----------

Changes in version 3.6:

NEW FEATURES

- JASPAR2018 data package

[JunctionSeq](https://bioconductor.org/packages/JunctionSeq)
-----------

Changes in version 1.7.5:

Structural revamp

- Changes to DESeq2 have made it necessary to copy C++ code from
  DESeq2. As a result: JunctionSeq once again requires compilation.

- Minor bugfixes.

- Added more useful error messages when encountering incompatible GFF /
  count files.

[KEGGprofile](https://bioconductor.org/packages/KEGGprofile)
-----------

Changes in version 1.19.3:

- Fix bugs in download_KEGGfile function, which was caused by the
  updates of KEGG web site.

Changes in version 1.19.1:

- Fix a bug in download_KEGGfile function, which was caused by the
  updates of KEGG web site.

[limma](https://bioconductor.org/packages/limma)
-----

Changes in version 3.34.0:

- read.idat() can now handle idat files in SNP format.

- New argument 'sort.by' for topSplice().

- diffSplice() now treats any NA geneid as a separate gene with one
  exon. Previously all NA geneids were treated as the same gene.

- Bug fixes for plotSplice() which, in some circumstances, was not
  highlight significant exons.

- Default value for 'style' in volcanoplot() changed to "p-value"
  instead of "p". This doesn't change the behavior of the function.

- volcanoplot() didn't work with MArrayLM objects produced by treat().
  Now fixed.

- The gene.info.file argument of alias2SymbolUsingNCBI() will now
  optionally accept a data.frame instead of a file name.

- alias2SymbolUsingNCBI() now always produces a data.frame.  Previously
  a single-column data.frame was returned as a vector.

- Bug fixes to camera() and cameraPR() when 'index' contains character
  row names.

- New argument 'restrict.universe' for the default kegga() method.  The
  default is now not to the restrict the universe to genes with KEGG
  annotation.

- goana.default() code is rewritten to more closely match
  kegga.default. Slight speed improvement. Prior probabilities are now
  computed using the restricted universe with GO annotation instead of
  the whole universe.

- Bug fix for kegga() when trend=TRUE or a prior.prob or covariate is
  specified. Previously the Down p-values were substantially too small
  and the Up p-values were too large.

- plotWithHighlights() checks whether 'status' is a factor and converts
  to character.

- Bug fixes for beadCountWeights(). Default is now set correct for
  'design' and the function now works correctly when 'y' is an EList
  object contain bead standard errors but not standard deviations.

[LymphoSeq](https://bioconductor.org/packages/LymphoSeq)
---------

Changes in version 1.4.1:

- Fixed bug in topSeqPlot

[maftools](https://bioconductor.org/packages/maftools)
--------

Changes in version 1.4.00:

NEW FUNCTIONS

- plotApobecDiff - plots differences between APOBEC enriched and non
  APOBEC enriched samples

- gisticOncoplot, gisticBubblePlot and gisticChromPlot to visualize
  GISTIC results

- somaticInteractions - to identify mutually exclusive/co-occuring gene
  sets

- genotypeMatrix - function to create genotype matrix

- mutCountMatrix - generate count matrix

SIGNIFICANT USER-LEVEL IMPROVEMENT

- Changes to MAF object: It now includes clinical data slot similar to
  PhenoData of expressionset objects.

- Changes to MAF object: Silent variants will be stored seperately in
  MAF object and won't be mixed with non-syn variants.

- Changes to MAF object: Oncomatrix is built on the fly whenever
  required, its no longer stored in MAF object.

- You can specify a manual list of variant classifications to be
  considered as non-synonymous via argument vc_nonSyn in read.maf

- Dropped mutExclusive function - Use somaticInteractions instead.

- Many sorting options and plotting improvements to oncplots.

- One can include q values from mutsig (or any similar program) as a
  side barplot in oncoplot

- rainfallPlot can detect hyper mutated genomic segments via
  ChangePoint detection method

- plotSignatures includes cosine similarity score and aetiology of
  detected signature

- readGistic has argument cnLevel to choose deep or shallow CN variants

- inferHeterogeneity includes MATH score in the plot

- Tumor_Sample_Barcodes remains as is; earlier '-' were converted to
  '.' in sample names

NON SIGNIFICANT CHANGES

- mafCompare output includes adjusted p-values

- trinucleotideMatrix output includes adjusted p-values for APOBEC
  enrichment

- added plot arguments to control title size and point size in
  lollipopPlot

BUG FIXES

- Major bug fix in signature analysis

- minor bug fixes in oncostrip

[matter](https://bioconductor.org/packages/matter)
------

Changes in version 1.3.8 (2017-10-26):

NEW FEATURES

- Added 'matter_fc' class for on-disk factors

SIGNIFICANT USER-VISIBLE CHANGES

- Added character encodings to 'matter_str' on-disk strings

BUG FIXES

- Fixed bug in 'matter_df' when subsetting by columns

- Fixed bug when coercing 'data.frame's to 'matter_df'

Changes in version 1.3.7 (2017-10-25):

BUG FIXES

- Fixed bug where default 'matter_df' chunksize was tiny

- Changed 'datamode' for 'matter_df' from 'list' to 'virtual'

Changes in version 1.3.6 (2017-10-20):

SIGNIFICANT USER-VISIBLE CHANGES

- Lists (ragged arrays) may now be heterogenous (elements must still be
  vectors)

- More memory-efficient initialization of data on disk when
  length(data) is 1

- Files are now created by default in constructors if they do not
  already exist

- When coercing with 'as.matter', names and/or dimnames are now
  retained

Changes in version 1.3.5 (2017-10-18):

SIGNIFICANT USER-VISIBLE CHANGES

- Improved length/dim defaults when data is given to constructors

BUG FIXES

- Trying to write to a read-only dataset now fails gracefully

- Fixed integer overflow when subsetting via subsetMatterMatrix

Changes in version 1.3.4:

BUG FIXES

- Fixed bug where NAs did not propogate correctly for integers

- It is now an error when 'matter()' cannot infer data structure

Changes in version 1.3.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Added new vignettes with benchmarks and examples on real data

BUG FIXES

- Fixed bug in delayed arithmetic operations with vectors

Changes in version 1.3.2:

NEW FEATURES

- Added class 'matter_list' for homogenous lists (ragged arrays)

- Added class 'matter_str' for on-disk strings

- Added class 'matter_df' for on-disk data frames

- Added 'as.matter' for function coercing to matter objects

Changes in version 1.3.1:

NEW FEATURES

- Updated base matter classes (matter_vec, matter_mat, matter_arr) to
  support 'Arith' and 'Compare' group generics

[MEAL](https://bioconductor.org/packages/MEAL)
----

Version: 1.8.0
Category: MAJOR CHANGES
Text:

Version: 1.8.0
Category: o Substitute AnalysisResults and AnalysisRegionResults by
        ResultSet
Text:

Version: 1.8.0
Category: o Substitute MethylationSet by GenomicRatioSet
Text:

Version: 1.8.0
Category: USER-VISIBLE CHANGES
Text:

Version: 1.8.0
Category: o Rename analysis functions
Text:

Version: 1.8.0
Category: o Create wrappers for each DMR methd
Text:

Version: 1.8.0
Category: NEW FEATURES
Text:

Version: 1.8.0
Category: o Add differences of variances analysis
Text:

Version: 1.8.0
Category: o Add new plot to simultaneously show all results of the same
        region
Text:

[meshes](https://bioconductor.org/packages/meshes)
------

Changes in version 1.3.1:

- new project site <2017-09-30, Fri>

[metabomxtr](https://bioconductor.org/packages/metabomxtr)
----------

Changes in version 1.11.1:

- *corrected bugs that occur when mixture models include at most 1
  predictor

[MetaboSignal](https://bioconductor.org/packages/MetaboSignal)
------------

Changes in version 1.7.12:

- MetaboSignal integrates KEGG metabolic and signaling networks with
  regulatory interactions reported in OmniPath and TRRUST.

- New functions: "MS2_ppiNetwork( )" and "MS2_mergeNetworks( )".

Changes in version 1.7.8:

- MetaboSignal includes a new function: "MS_getPathIds()" which allows
  retrieving the identifiers (IDs) of all metabolic and signaling KEGG
  pathways of a given organism.

- "MS_keggNetwork()" can now transform KEGG IDs into Entrez gene IDs
  (use expand_genes = TRUE, convert_entrez = TRUE).

- "MS_filterNetwork()" has been renamed as "MS_topologyFilter".

Changes in version 1.7.2:

- Major changes in function names:

- "MetaboSignal_matrix()" renamed as "MS_keggNetwork"

- "MetaboSignal_distances()" renamed as "MS_distances"

- "MetaboSignal_NetworkCytoscape()" renamed as
  "MS_shortestPathsNetwork"

- "MS_ToCytoscape()” renamed as "MS_exportCytoscape"

- "MS_ChangeNames()" renamed as "MS_changeNames"

- "MS_FilterNetwork()" renamed as "MS_filterNetwork"

- "MS_FindKEGG()" renamed as "MS_keggFinder"

- "MS_GetKEGG_GeneID" renamed as "MS_convertGene"

- "MS_NodeBW()" renamed as "MS_nodeBW"

- "MS_ReplaceNode()" renamed as "MS_replaceNode"

- "MS_RemoveNode()" renamed as "MS_removeNode"

- "MS_FindMappedNodes()" renamed as "MS_findMappedNodes"

- New functions:

- "MS_tissueFilter()": allows filtering a network based on tissue
  expression data.

- Functions removed:

- "MS_interactionType()": interactions can be now retrieved directly
  from MS_keggNetwork.

- Functionality changes:

- "MS_keggNetwork()": output networks are now formatted as 3-column
  matrices, where the third column indicates interaction type. Unlike,
  MetaboSignal_matrix, “MS_keggNetwork()" does not perform network
  filtering by tissue expression data.  Tissue filtering is now done
  with "MS_tissueFilter()". Also, the compound nodes of the networks
  are not only metabolites, but also drugs and glycans.

- "MS_distances()", "MS_shortestPathNetworks()": use network nodes
  (i.e. KEGG IDs) as source nodes, not entrez IDs or gene symbols.
  However, the function "MS_convertGene()" allows to easily convert
  entrez IDs or gene symbols into KEGG IDs.

Changes in version 1.7.1:

- Significant improvements in the tissue-filtering option from
  "MetaboSignal_matrix()".

- Cytoscape .txt files are now exported with a header.

[metaseqR](https://bioconductor.org/packages/metaseqR)
--------

Changes in version 1.17.4 (2017-10-06):

NEW FEATURES

- New option: utr.flank for counting reads in 3' UTR flank regions when
  analyzing Quant-Seq data.

- New gene filter: genes where x samples present less than y counts are
  excluded from statistical testing. x is determined as a fraction of
  available samples.

BUG FIXES

- Fixed bug in problematic backwards compatibility of strandedness in
  targets files which caused strandedness to be ignored (thans to
  Martin Reczko, BSRC 'Alexander Fleming').

[metavizr](https://bioconductor.org/packages/metavizr)
--------

Changes in version 1.3.20:

- Standalone mode introduced, a version of epiviz with reduced
  capabilities is now included as part of epivizr. The epiviz web app
  is run locally using 'httpuv's http server

- Add and remove seqinfo (e.g., chromosome info) to any epiviz session

Changes in version 1.3.11:

- Add NEWS file

- Update documentation on 'slideshow' function

Changes in version 1.3.10:

- Changed default on 'slideshow' to show all ranges

Changes in version 1.3.9:

- Added 'heatmapChart' convenience function

Changes in version 1.3.8:

- Fixed bug in 'startEpiviz' not sending 'seqName' parameter correctly

Changes in version 1.3.7:

- Fixed bug in 'EpivizBpData' that sent 'metadata' info in wrong format

Changes in version 1.3.6:

- Changes slots using lists in 'EpivizDeviceMgr' to environments to
  avoid crashing RStudio due to inspection of manager objects

Changes in version 1.3.5:

- Fails gracefully on daemonization request on Windows

- Deprecates the 'proxy' argument to 'startEpiviz'

Changes in version 1.3.4:

- Upgrading to Epiviz v2 webapp

[MetCirc](https://bioconductor.org/packages/MetCirc)
-------

Changes in version 1.5.0:

- no changes, check if package passes R CMD build and R CMD check
  without any error messages and vignette can be run without any errors
  [2017-10-20 Fri]

[methimpute](https://bioconductor.org/packages/methimpute)
----------

Version: 1.0.0
Category: INITIAL RELEASE
Text:

[methylInheritance](https://bioconductor.org/packages/methylInheritance)
-----------------

Changes in version 1.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Functions runObservation() and runPermutation() don't return the
  result anymore. The results are saved in RDS files. The results can
  be loaded through the loadAllRDSResults() function.

- New plotConvergence() function enable visualization of convergence

- New "saveInfoByGeneration" parameter which generates RDS files
  containing information about each generation for each permutation

BUG FIXES AND IMPROVEMENTS

- Major changes in parallel processing to limit memory consumption

[methylKit](https://bioconductor.org/packages/methylKit)
---------

Changes in version 1.3.8:

IMPROVEMENTS AND BUG FIXES

- changes to methCall(): more verbose output, reduce -Wunsigned
  warnings, check for fragmented chromosome order, write conversion
  stats to extra file

- add Tests to check for chromosome disorder

Changes in version 1.3.7:

IMPROVEMENTS AND BUG FIXES

- if methylDiff object contains NA in a qvalue/meth.diff column the
  getMethylDiff function returns an error. This is now fixed:
  https://github.com/al2na/methylKit/issues/83

- Changed how the methylBase version of regionCounts handle missing
  values.  Before, a missing value made the whole region/tile missing.
  Now, a missing value is treated as a site with zero coverage.
  Contributed by Karl Nordström:
  https://github.com/al2na/methylKit/pull/77

Changes in version 1.3.6:

- changes to methSeg() function: update function description by beeing
  more explicit about the sorting, sort(x,ignore.strand=TRUE) is now
  called once if needed, update tests

Changes in version 1.3.5:

IMPROVEMENTS AND BUG FIXES

- fix problem with methSeg() error occurring when GRanges is not sorted
  by position

Changes in version 1.3.4:

IMPROVEMENTS AND BUG FIXES

- fix problem with assocComp() on methylBaseDB

Changes in version 1.3.3:

IMPROVEMENTS AND BUG FIXES

- changes to calculateDiffMeth() function: add a message explaining the
  calculation procedure based on either two or more treatment groups,
  fixed bugs that prohibited the use of treatment groups other than 0/1
  in two group case, allow for more than two groups, update tests

- changes in methSeg() function: update function documentation to
  inform that provided Granges has to be sorted, add check if granges
  is sorted and contains at least one meta.col

Changes in version 1.3.2:

IMPROVEMENTS AND BUG FIXES

- quick fix for long filename issue: in unite() and
  calculateMethylDiff() instead of concatenated sample_ids, filename
  will be "methylBase" concatenated with either 13-char random string
  or provided suffix

- updated vignette: added two short FAQs

Changes in version 1.3.1:

IMPROVEMENTS AND BUG FIXES

- fix conversion rate error in methCall()

[microbiome](https://bioconductor.org/packages/microbiome)
----------

Changes in version 0.99.8 (2017-08-21):

- Added sample names in divergence output

- coreset option added in divergence

- plot_taxa_prevalence argument "detection" added

- transform argument removed from plot_composition for simplicity

- data/DynamicsIBD data set and associated documentation removed to
  shrink the package

- inst/extdata/qiita data sets and associated documentation removed to
  shrink the package

- Baxter data removed to save space

- Utilities that do not belong to the package moved to maintenance
  branch

Changes in version 0.99.6 (2017-07-21):

- Fixed bug in transform option "clr"

- "lineplot" option added in plot_composition

Changes in version 0.99.5 (2017-07-07):

- fixed vignette headers

- heat function now takes also labels as input for order.rows and
  order.cols

- added age_group and bmi_group for easy access to standard groupings
  of these variables

- package homepage host location changed to master:docs/

Changes in version 0.99.3 (2017-06-05):

- Reduced the number of dependencies

- Removed marginal functionality

Changes in version 0.99.1 (2017-05-15):

- Major rewrite of the package

- Switched to support the phyloseq class

- HTML vignette added

- A separate online tutorial added

- Removed less essential functionality and dependencies

Changes in version 0.99.0 (2014-09-14):

- First draft of the BioC release version

Changes in version 0.12.14 (2012-05-10):

- added pitchipdb in PITChip array name list in profiling.R

Changes in version 0.12.12 (2012-05-02):

- added MITChip level 0

- added separate get.phylogeny.MITChip function to include level 0 for
  MITChip

Changes in version 0.12.10 (2012-04-17):

- automatized the levels and methods output within run.profiling.script

- now providing also non-logarithmized data matrices in the
  run.profiling.script output

- added in run.profiling.script save.image(paste(params$wdir,
  "/tmp.RData", sep = "")) to enable the use of data from the current
  workspace during a session

- multcomp replaced by parallel in dependencies

- Removed "Use existing data from the workspace" (useDB) option.
  Problems when called from within a function. More clear and modular
  when the run.profiling.script is used for preprocessing, and any new
  plots are performed as separate steps based on the profiling output
  files if needed.

- Moved run.profiling.script plots to their own functions and added
  usage examples to vignette.

Changes in version 0.12.08 (2012-04-05):

- detailed protocol added to the vignette

- hierarchical clustering method moved from complete to ward in
  profiling script as it was in some of the previous versions

- 'level 0' added in MITChip profiling

- removed featureprofile from output files to save disk space,
  typically not needed

- added remove.nonspecific.oligos option to run.profiling.script and
  FetchHITChipAtlas functions (through get.phylogeny)

- summarize.probesets modified to include levels both in the form of
  "level.1" and "level 1"

Changes in version 0.12.06 (2012-04-04):

- fixed PITChip-specific issues

- removed parallel from dependencies to allow compatibility with
  R-2.12.2

Changes in version 0.12.05 (2012-03-29):

- removed background correction from profiling script
  (run.profiling.script)

- kept standard normalizations: none, minmax, quantile - others removed

- now using 16S phylogeny in profiling script. Other option removed

- polishing of the functions

- merged funcions from the phyloarray package

- removed minmax, quantnorm from profiling.R functions

Changes in version 0.12.04 (2012-03-28):

- removed the additional background correction step which resulted in
  considerable information loss based on external validations

Changes in version 0.12.01 (2012-03-23):

- PlotMatrix updates

Changes in version 0.3.63 (2012-12-24):

- separated MySQL/preprocessing functions into a distinct package,
  HITChipDB

Changes in version 0.3.43 (2012-09-28):

- added relative.abundance function

- polished estimate.diversity function

- added functions: richness, evenness, diversity

Changes in version 0.3.15 (2012-08-31):

- FetchHITChipAtlas and run.profiling.script validated against each
  other

- default imputation removed from run.profiling.script; this had
  considerable effect on normalized data

Changes in version 0.3.13 (2012-08-30):

- standardized the pipeline; output validated against previous atlas
  version

Changes in version 0.3.05 (2012-08-22):

- added MDS.classical and MDS.isometric to project.data

Changes in version 0.3.02 (2012-06-20):

- calculate.stability function added

- added L2->species option to levelmap function

- estimate.diversity function updated. Detection thresholds applied
  with richness and evenness

- estimate.min.threshold mode detection improved

- cross.correlate function added

Changes in version 0.3.01 (2012-06-19):

- draft version for GitHub

Changes in version 0.2.04 (2012-06-15):

- RMySQL removed from explicit dependencies

Changes in version 0.2.03 (2012-06-06):

- hitchip.phylodistance added

- cross-hyb control functions added

Changes in version 0.2.02 (2012-06-03):

- annotation functions added

Changes in version 0.2.01 (2012-05-31):

- removed absolute scale matrices from run.profiling.script output. The
  idea is that just one log10 file for each phylogenetic level is
  outputted, and the user can then use other reading/analysis routines
  to read the data, convert to original non-log domain if necessary,
  etc. This is to provide a unique output format to avoid mistakes.

- save.data and tree.display removed from run.profiling.script
  arguments

- combined get.phylogeny and get.phylogeny.MITChip

- moved plotting functions from run.profiling.script into a separate
  function profiling.hclust;

- removed the option to read profiling parameters from R file in
  ReadParameters function

- output file name phylogenyinfo.tab changed into
  phylogenyinfo-nspecies.tab

- made run.profiling.script modular with respect to preprocessing,
  plotting, data storage, and logging

- added instructions in the vignette on reading the final data

Changes in version 0.0.12 (2012-03-15):

- multicore removed from dependencies

- vignette updated

- R-2.12.1 required

Changes in version 0.0.10 (2012-03-06):

- standard profiling script included

- added selected.samples in fetch.sample.info

Changes in version 0.0.08 (2012-02-25):

- Atlas preprocessing finalized. For usage examples, see vignette.

Changes in version 0.0.05 (2012-02-22):

- atlas functions added and polished

Changes in version 0.0.02 (2012-02-06):

- profiling_010.r script functions and dependencies added

Changes in version 0.0.01 (2012-01-12):

- First version based on profiling script version 010

[minfi](https://bioconductor.org/packages/minfi)
-----

Changes in version 1.23:

- Fixed as() (coercion) from RGChannelSetExtended to RGChannelSet, to
  support the argument extended=TRUE in read.matharray(). The core
  issue is the new("RChannelSetExtended") is an invalid object because
  it does not have correct elements of the assay slot.  Instead of
  addressing this, I used a check for ncol=0, nrow=0 in the coercion
  function which asssumes the presence of correctly named assays.
  Original issue report by Stewart Morris <swmorris@exseed.ed.ac.uk>.

- Improved (made prettier) the printing of messages in read.metharray()
  and friends.

- Changed seqlevels(..., force = TRUE) to seqlevel(..., pruning.mode =
  "coarse").

[miRmine](https://bioconductor.org/packages/miRmine)
-------

Version: 0.99.8
Text:

Version: 0.99.1
Text:

Version: 0.99.0
Text:

[miRsponge](https://bioconductor.org/packages/miRsponge)
---------

Version: 1.1.1
Category: Update miRsponge.R <2017-10-01, Sun
Text:

Version: 1.1.0
Category: Update netModule function <2017-09-07, Thur
Text:

Version: 0.99.17
Category: Update DESCRIPTION <2017-08-23, Wed
Text:

Version: 0.99.16
Category: Update moduleSurvival function <2017-08-08, Tues
Text:

Version: 0.99.15
Category: Update test_miRsponge.R <2017-07-22, Sat
Text:

Version: 0.99.14
Category: Update test_miRsponge.R <2017-07-22, Sat
Text:

Version: 0.99.13
Category: Improve the code of miRsponge.R <2017-07-22, Sat
Text:

Version: 0.99.12
Category: Update related references in man folder <2017-07-20, Thur
Text:

Version: 0.99.11
Category: Improve the code of miRsponge.R <2017-07-20, Thur
Text:

Version: 0.99.10
Category: Update README.md and miRsponge.Rmd <2017-07-18, Tues
Text:

Version: 0.99.9
Category: Update moduleSurvival function <2017-07-14, Fri
Text:

Version: 0.99.8
Category: Update cernia.Rd and hermes.Rd <2017-07-13, Thur
Text:

Version: 0.99.7
Category: Update dtHybrid function <2017-07-13, Thur
Text:

Version: 0.99.6
Category: VERSION BUMP REQUIRED <2017-07-13, Thur
Text:

Version: 0.99.5
Category: Update DESCRIPTION file <2017-07-13, Thur
Text:

Version: 0.99.4
Category: Rename functions in R folder <2017-07-13, Thur
Text:

Version: 0.99.3
Category: Update LICENSE file <2017-07-13, Mon
Text:

Version: 0.99.2
Category: Improve the code format of miRsponge.R <2017-07-13, Mon>.
Text:

Version: 0.99.1
Category: Improve the code format of miRsponge.R <2017-07-13, Mon>.
Text:

Version: 0.99.0
Category: This is the first version of miRsponge package. If any bugs,
        please let me know. Contact Email: zhangjunpeng_411@yahoo.com
Text:

[monocle](https://bioconductor.org/packages/monocle)
-------

Changes in version 2.5.0:

- Various plot utilities for visualizing complex developmental
  trajectory (plot_complex_cell_trajectory), multi-way kinetic curve
  and multi-way heatmap (plot_multiple_branches_pseudotime,
  plot_multiple_branches_heatmap).

- ClusterCells now in principle supports clustering for 100 k+ cells
  which depends on a procedure of kNN based density peak clustering
  algorithm. densityClust package version 0.3 are required for users to
  use this functionality.

- New functions supporting importing scater or Seurat cds into monocle
  or vice versa (importCDS, exportCDS).

[motifcounter](https://bioconductor.org/packages/motifcounter)
------------

Changes in version 1.1.9:

- Removed compiler warnings

- Changes to the vignette

Changes in version 1.1.8:

- Internal refactoring

- Changed citation information

[MotifDb](https://bioconductor.org/packages/MotifDb)
-------

Changes in version 1.20:

NEW FEATURES

- Now 13 sources, 52 organisms and 8369 motifs

- New sources: HOCOMOCOv10, HOMER, jaspar2016, SwissRegulon

- New methods: motifToGene, geneToMotif, associateTranscriptionFactors

SIGNIFICANT USER-VISIBLE CHANGES

- The new methods mentioned above provide both strict and permissive
  mapping from motifs to their cognate transcription factors, thus
  supporting reproducible regulatory network construction.

[motifStack](https://bioconductor.org/packages/motifStack)
----------

Changes in version 1.21.7:

- add newline at the end of RUNX1.pfm

Changes in version 1.21.6:

- update the documentations

Changes in version 1.21.5:

- add class psam

- add function plotAffinityLogo

Changes in version 1.21.4:

- add function importMatrix

Changes in version 1.21.3:

- add one more option for motifPiles and motifCircos

Changes in version 1.21.2:

- fix a bug in plotMotifLogo for pcm

Changes in version 1.21.1:

- update the vignettes.

[mpra](https://bioconductor.org/packages/mpra)
----

Changes in version 0.99:

- Initial release to Bioconductor.

[MSnbase](https://bioconductor.org/packages/MSnbase)
-------

Changes in version 2.3.14:

- Use `normalizePath` to force absolute file paths in `readMSData`.

Changes in version 2.3.13:

- Add write support for MSnExp and OnDiskMSnExp objects allowing to
  save the MS data to mzML or mzXML files. <2017-09-15 Fri>

Changes in version 2.3.12:

- Keep `protocolData` in isobaric quantification; fixes #265

Changes in version 2.3.11:

- Amend `addIdentificationData` when sourceInfo reports multiple files
  and when scores are missing from the identification results (closes
  #261).

- Don't overwrite `processingData` slot when creating an `MSnSet`
  object (closes #264).

Changes in version 2.3.10:

- New `isCentroidedFromFile` function <2017-08-11 Fri>

- Add msLevel slot to Chromatogram object <2017-08-16 Wed>

- Add msLevel argument to chromatogram,MSnExp method <2017-08-16 Wed>

- `calculateFragments` now just calculate fragments for all `n - 1L`
  bonds (before it incorrectly adds an additional bond; fixes #248)
  <2017-08-20 Sun>

- Add `isEmpty` methods for `Chromatogram` and `Chromatograms` objects
  <2017-09-05 Tue>

- plot,Chromatogram[s] creates an empty plot and returns a warning if
  the Chromatogram[s] object is empty (issue #249) <2017-09-05 Tue>

Changes in version 2.3.9:

- Using new mzR::openIdfile backend to add identifcation data to raw
  and quantitative data (see issue #232) <2017-07-28 Fri>

- New utils functions: factorsAsStrings, makeCamelCase and
  reduce,data.frame <2017-07-29 Sat>

- Coerce mzRident to data.frames <2017-07-29 Sat>

- Add phenoData slot to Chromatograms class<2017-08-02 Wed>

- new readMzIdData function to read mzId files as data.frames (uses the
  new coerce,mzRident,data.frame method) <2017-08-03 Thu>

- new filterIdentificationDataFrame function to filter PSM data.frames
  as produced by readMzIdData. Also used in the addIdentificationData
  methods. <2017-08-03 Thu>

- readMSData has a new mode argument to set onDisk or inMemory
  (default) mode <2017-08-10 Thu>

Changes in version 2.3.8:

- New infrastructure for chromatogram data <2017-06-24 Sat>

- Change naming scheme for spectra: FFILEID.SSPECTRUMID, e.g.
  F01.S0001. Before it has been XSPECTRUMID.FILEID. The new naming
  scheme changes the order of the spectra. See #255 (and PR #256) for
  details <2017-06-25 Sun>.

Changes in version 2.3.7:

- export filterEmptySpectra

Changes in version 2.3.6:

- Brutally remove xic and chromatogram functions/methods, to be
  replaced by the Chromatogram[s] infrastructure <2017-06-15 Thu>

Changes in version 2.3.5:

- Fix superscript syntax in demo vignette <2017-06-14 Wed>

Changes in version 2.3.4:

- Use the injection time from mzR (see PR #109) which in now in seconds
  (was in milliseconds) <2017-06-13 Tue>

Changes in version 2.3.3:

- Rewrite `getColsFromPattern` and `getRowsFromPattern` and add unit
  tests <2017-05-11 Thu>.

- Add `.filterNA` and rewrite `filterNA` for `matrix` and `MSnSet`
  <2017-05-11 Thu>.

- Convert main MSnbase-demo vignette to Rmd/html <2017-05-27 Sat>

- `naplot` gains a `reorderRows` and `reorderColumns` argument
  <2017-06-05 Mon>.

Changes in version 2.3.2:

- Rewrite `utils.clean`. It now keeps just the zeros in the direct
  neighbourhood (see #210) <2017-05-04 Thu>.

Changes in version 2.3.1:

- Introduce "NTR" method for `combineFeatures` <2017-04-26 Wed>.

- Rewrite `nQuants` and `featureCV` to avoid returning of rows for
  empty factors; see PR #208 for details <2017-04-28 Fri>.

- Add $,pSet method to easily access columns in the phenoData (see
  #203)

Changes in version 2.3.0:

- Version bump for new Bioc devel.

[MSnID](https://bioconductor.org/packages/MSnID)
-----

Changes in version 1.11.1:

- Switch from mclapply to parLapply in filter optimizaition. This fixes
  a bug "Assertion failure at kmp_runtime.cpp(6480): __kmp_thread_pool
  == __null." Also allows to run parallel on Window OS

[msPurity](https://bioconductor.org/packages/msPurity)
--------

Changes in version 1.3.9:

- Added very basic SIMS stitch compatibility

- pcalc can handle NAs

- Update of purityX to handle obiwarp RT correction (requires recording
  the RT RAW at an earlier step)

- bug fix for when library spectra is bigger than target spectra
  (thanks Martin)

Changes in version 1.3.1:

- Add spectral matching functionality for LC-MS/MS

[MSstats](https://bioconductor.org/packages/MSstats)
-------

Version: 3.9.7
Text: NEW FEATURES - cluster (default=1) is no longer available for
        groupComparison function, due to memory issue.

Version: 3.8.6
Text: NEW FEATURES - can cluster (default=1) for dataProcess and
        groupComparison function (Thanks John!!)

Version: 3.8.5
Text: BUG FIXES - PDtoMSstatsFormat : three options are added for
        outputs from different versions of PD. (Thanks to Felipe!)  -
        which.quantification - which.proteinid - which.sequence

Version: 3.8.4
Date: 2017-08-28
Text: BUG FIXES - SkylinetoMSstatsFormat : DDA case lost ‘StandardType’
        column after summing peaks. Fixed. (Thanks, Nick)

Version: 3.8.3
Date: 2017-07-13
Text: NEW FEATURES - SpectronauttoMSstatsFormat : if PG.Qvalue is
        available, filter out the data with greater than 0.01.  -
        dataProcessPlots, groupComparisonPlots, modelBasedPlots : with
        address=FALSE option, one plot a time can be drawn in the panel
        and won't be saved as in pdf.  BUG FIXES - designSampleSize :
        fix the calculation of variance (Thanks, Tsung-Heng) -
        SkylinetoMSstatsFormat : when Condition and BioReplicate
        columns are NA, there was issue for merge with annotation.  -
        SkylinetoMSstatsFormat : fix the bug to recognize the protein
        with one peptide only for the option:
        'removeProtein_with1Peptide = TRUE' - dataProcess : when
        cutoff.lower is negative, with maxQuantileforCensored option +
        censoredInt='0', zero log2 endogenous intensity should be
        censored.  - ProgenesistoMSstatsFormat : handle inputs with
        some limited columns. such as no Spectral.counts columns.

Version: 3.8.2
Date: 2017-04-21
Text: NEW FEATURES - required ‘Fraction’ information in annotation for
        pre-processing - dataProcess function is updated for merge
        fractions BUG FIXES - warning message during dataProcessPlots
        for profile plot is not shown anymore.

[multiClust](https://bioconductor.org/packages/multiClust)
----------

Changes in version 6-19-16:

- -Package version pushed to 1.0.3 -Updated package title and abstract
  to reflect that of publication in Cancer Informatics DOI:
  10.4137/CIN.S38000.

Changes in version 3-2-16:

- -Package version pushed to 0.99.6 -Added Biobase, GEOquery, and
  preprocessCore to package suggests -Minor revisions in the package
  vignette

Changes in version 2-13-16:

- -Package version changed to 0.99.5 -Added option to specify FDR
  cutoff when using the Adaptive GMM method in the number_probes
  function. -Updated documentation in the vignette

Changes in version 2-11-16:

- -Packaged changed to version 0.99.4 -Revised the code in the vignette
  -Added explanation of using RNA-seq data with package in vignette
  -Revised code documentation for probe_ranking function

[MWASTools](https://bioconductor.org/packages/MWASTools)
---------

Changes in version 1.1.1:

- MWASTools includes new functions: "MWAS_heatmap()",
  "MWAS_KEGG_pathways()", "MWAS_KEGG_network()", and
  "MWAS_KEGG_shortestpaths()".

[mzR](https://bioconductor.org/packages/mzR)
---

Changes in version 2.11.11:

- Fix problem in writeMSData: ensure precursor data is saved even if
  precursor scan is not available (see MSnbase issue #245).

Changes in version 2.11.10:

- Fix problem that can cause a SEGFAULT in writeMSData/copyWriteMSData
  when MS data with spectra linking to missing precursor scans is saved
  (issue #129).

Changes in version 2.11.9:

- Update peaks man page with details about spectrumId, acquisitionNum
  and seqNum

Changes in version 2.11.8:

- Add contributions guide with code of conduct.

- Update installation instructions for Mac.

- Report the spectrum ID in the header data.frame (column spectrumId).

- Fix in copyWriteMSData and writeMSData ensuring that MSn data is
  correctly

- Import pwiz r11174 fix for mzML without <componentList> (see #113).

Changes in version 2.11.7:

- Nothing yet.

- Import fix by Brian Pratt (pwiz r11174) for mzML without
  <componentList> Another way to fix #113

- Removing mz5 support from manual page, as currently unsupported.

Changes in version 2.11.6:

- runInfo returns the run start time stamp from files providing this
  information (mzML files).

Changes in version 2.11.5:

- writeMSData and copyWriteMSData functions enabling to export MS data
  to mzML or mzXML files.

Changes in version 2.11.4:

- Use full TMT file pattern to select a single file

Changes in version 2.11.3:

- Read ion injection time from mzML files and add it to the data.frame
  returned by the header function.

Changes in version 2.11.2:

- New getScanHeaderInfo and getAllScanHeaderInfo implementations for
  the pwiz backend. Fixes issue #106 and issue #216 in MSnbase.

Changes in version 2.11.1:

- Change default I/O backend from Ramp to pwiz.

Changes in version 2.11.0:

- Bioc devel 3.6

[NanoStringQCPro](https://bioconductor.org/packages/NanoStringQCPro)
---------------

Changes in version 1.9.1 (2017-10-11):

- Fixed vignette build issues following recent overhaul of BiocStyle
  (cf.
  https://github.com/Bioconductor/BiocStyle/blob/master/announcement.md
  (91e6565 on Aug 9)).

[ndexr](https://bioconductor.org/packages/ndexr)
-----

Changes in version 0.99.10:

- **Breaking changes of the function names!**

- Changed the naming of the functions: dots are replaced by
  underscores, so that the functions may not be confused with S3
  methods.  The naming convention itself doesn't change!  E.g. the new
  name of the function "ndex.network.update.aspect" is now
  "ndex_network_update_aspect"

- Moved packages httr, jsonlite, plyr and tidyr from "Depends" field to
  "Imports" (Import into NAMESPACE not necessary, because external
  package function names are always explicitly qualified)

- Fixed a bug in ngraph_fromRCX, which prevented the resulting ngraph
  object from node attributes to be set correctly. This also lead the
  following warnings: "Warning in vattrs[[name]][index] <- value :
  number of items to replace is not a multiple of replacement length"

- Changed functions from passing a quoted string (e.g. host =
  "ndexConf$connection$host") as default parameter to using the actual
  object (e.g. host = ndexConf$connection$host)

- Made some minor changes to NDExConnection object; removed an
  undesired warning.

- Implement a print method for the classes "NDExConnection" and "RCX"
  to provide the user with a useful summary of this complicated object

- changed to the usage of message() rather than cat() for verbose
  outputs, where it wasn't already done

- Specified a single 'Maintainer' representing the primary contact for
  maintenance issues related to this package.

- Removed the .gitignore file from this directory

- Changed installation instructions in the vignette to Bioconductor

Changes in version 0.99.0:

- Pushing version number to 0.99.0 according to Bioconductor checklist

- Initial submission to Bioconductor

[normr](https://bioconductor.org/packages/normr)
-----

Version: 1.3.1
Category: FEATURES
Text: `summary()` output more concise

Version: 1.3.1
Category: FEATURES
Text: `getEnrichment()` takes `F` for specifying a desired foreground
        component for standardization of the enrichment

Version: 1.3.1
Category: FEATURES
Text: `iterations` argument allows running multiple fits with different
        starting values

Version: 1.3.1
Category: FEATURES
Text: T Filter threshold can now be specified with `minP` argument for
        less stringent filtering

Version: 1.3.1
Category: FEATURES
Text: Added a normR overview scheme to the vignette

Version: 1.3.1
Category: FEATURES
Text: Qvalue computation improvement by specifying range for prop. of
        H_0

Version: 1.3.1
Category: BUGFIXES
Text: Fixed the erroneous tiling of supplied GR objects in
        char,char,GRanges

Version: 1.3.1
Category: BUGFIXES
Text: `diffR()` is more robust by doing a label-switched fit also

Version: 1.3.0
Text:

Version: 3.6
Text:

Version: 3.5
Text:

Version: 3.4
Text:

[Onassis](https://bioconductor.org/packages/Onassis)
-------

Changes in version 0.99.8:

- SRAdbHandler

- GEOHandler

- EntityFinder-methods

- similarity-methods

- AllGenerics

- AllClasses

- AllFunctions

- CMdict-methods

- Init_methods

- onLoad

- SupportingFunctions

- CMoptions-methods

[oncomix](https://bioconductor.org/packages/oncomix)
-------

Changes in version 0.99.0:

- First version (0.99.0) of oncomix is available.

- Paper to follow soon!

[OncoSimulR](https://bioconductor.org/packages/OncoSimulR)
----------

Changes in version 2.7.2 (2017-09-27):

- genot_to_adj_mat in C++.

- fast_peaks (for no backmutation cases).

- Better explanation and testing of peaks and valleys.

- Clarified simOGraph transitive reduction.

- Better handling of ti corner cases.

- Magellan reading fuctions adapted to output of newer (as of 2017-07)
  version of Magellan.

- sorting gene names in allGenotypes_to_matrix.

- sampledGenotypes: genotype names with sorted gene names.

[oneSENSE](https://bioconductor.org/packages/oneSENSE)
--------

Changes in version 0.99.1 (2017-07-04):

MODIFICATION

- Formatting changes

- Added Vignette

- Removed assignments to global variables

Changes in version 0.99.0:

MODIFICATION

- renamed package name from onesense to oneSENSE

[ontoProc](https://bioconductor.org/packages/ontoProc)
--------

Version: 1.0.1
Text:

[oposSOM](https://bioconductor.org/packages/oposSOM)
-------

Version: 1.15.1
Text:

[OPWeight](https://bioconductor.org/packages/OPWeight)
--------

Changes in version 1.0.0:

- First stable release with Bioconductor (2017-06)

[pathview](https://bioconductor.org/packages/pathview)
--------

Changes in version 1.17.7:

- updated combineKEGGnodes.R, i.e. changed the stop error to a warning
  message when a group node include different node types, i.e. both
  gene and compound. This problematic KEGG node definition does exist
  in pathway 04136 Autophagy - other, and caused error when calling
  pathview with kegg.native = F.

- fixed bug in mol.sum introduced at 1.15.1, i.e. indexing using
  which(eff.idx). This problem only affect direct call of mol.sum with
  sum.method other than "sum" and "mean".

Changes in version 1.17.1:

- major expansion in korg, which now include both KEGG and NCBI
  taxonomy IDs, two more gene ID types, i.e. NCBI protein and uniprot
  IDs. In addition, Entrez or NCBI Gene IDs are discontinued for most
  prokaryotes.

- korg now include 4800 KEGG species, in the meantime, an updated
  version of korg is now checked out from Pathview Web server each time
  pathview package is loaded.

[PathwaySplice](https://bioconductor.org/packages/PathwaySplice)
-------------

Changes in version 0.99.0:

- Initial version(2016-12-09)

[philr](https://bioconductor.org/packages/philr)
-----

Changes in version 1.3.1:

USER-VISIBLE CHANGES

- Squashed Bug in philr and philrInv handling of vector input

[piano](https://bioconductor.org/packages/piano)
-----

Changes in version 1.18.0:

DOCUMENTATION

- Switched to accommodate documentation in r-files using roxygen

BUG FIXES

- Fix print addInfo for GSC bug

Changes in version 1.16.4:

BUG FIXES

- Fix R CMD check NOTE on calling require(snowfall)

Changes in version 1.16.3:

BUG FIXES

- Fix BiocGenerics::sd and stats::sd collision when loading

Changes in version 1.16.2:

BUG FIXES

- Fix error: "Error in if (maxP > -minP) { : missing value where
  TRUE/FALSE needed" that appeared occasionally during runGSA for
  method gsea. This happened when the running sum was only positive or
  only negative, leading to that the if-statement broke due to missing
  values.

Changes in version 1.16.1:

BUG FIXES

- Fix biomaRt Rchunk in vignette to avoid build error

[Pigengene](https://bioconductor.org/packages/Pigengene)
---------

Changes in version 1.3.8 (2017-09-13):

Changes in existing functions

- Order of conditions in pheatmap.type can now be determined by the
  user.

Changes in version 1.3.6 (2017-08-20):

Changes in existing functions

- A bug in gene.mapping () function fixed to better map probe IDs.

Changes in version 1.3.4 (2017-07-31):

Changes in existing functions

- The doTranspose argument added to the heatmap.type() function.

[polyester](https://bioconductor.org/packages/polyester)
---------

Version: 1.99.3
Text: NB function now exported

Version: 1.99.3
Text: note that version 1.99.3 on GitHub was version 1.1.0 on
        Bioconductor.

Version: 1.99.2
Text: bug fix in fragment generation (last 2 bases of transcript were
        never sequenced)

Version: 1.99.1
Text:

[pRoloc](https://bioconductor.org/packages/pRoloc)
------

Changes in version 1.17.5:

- Filtering for unique features when running plot2D with t-SNE method
  <2017-10-15 Sun>

Changes in version 1.17.4:

- Added new (private) dimred function that computes dimensionality
  reduction <2017-06-05 Mon>

- Add F1000research workflow to citations <2017-06-22 Thu>

- Classification functions now return the classification score matrix
  for all classes as a single column in fData, rather that each class
  as its own fData column. <2017-09-01 Fri>

Changes in version 1.17.3:

- Convert vignettes to Rmd with html output <2017-05-25 Thu>

- Import, rather than suggest Rtsne <2017-05-25 Thu>

Changes in version 1.17.2:

- phenoDisco speed improvements and added support for t-SNE <2017-05-19
  Fri>

Changes in version 1.17.1:

- Support Rtsne's new pca_center and pca_scale arguments <2017-05-02
  Tue>

Changes in version 1.17.0:

- Version bump for Bioc devel 3.6

[pRolocGUI](https://bioconductor.org/packages/pRolocGUI)
---------

Changes in version 1.11.2:

- Fix links in vignette to point to new html pRoloc vignettes
  <2017-05-30 Tue>

Changes in version 1.11.1:

- Avoid computing dimensionality reduction at every reactive rendering,
  assuring that other, slower methods, in particular t-SNE, can be used
  <2017-05-20 Sat>

Changes in version 1.11.0:

- New version for Bioc devel 3.6

[Prostar](https://bioconductor.org/packages/Prostar)
-------

Changes in version 1.9.15:

BUG FIXES

- When the aggregation step has been performed, the interface switches
  to the first tab of the 'Descriptive Statistics' in order to view
  informations aout the new dataset (the protein one).

- Implementation of a parallel version of the function which saves the
  (new) protein dataset after the aggregation step.

- Disable the extra row appearing in the metadata table when
  convertinga text file to a MSnSet file.

- Disable the extra row appearing in the metadata table when
  convertinga text file to a MSnSet file.

- A new package (readxl) is used to read xls or xlxs files. In certain
  circumstances, the functions of the previsous package openxlsx is not
  able to decode properly Excel files.

- When converting a new (text or Excel) file in Prostar : the missing
  values were not registered as expected. Especially, they did not
  appear in blue in the table above the volcanoplot. Bug fixed

- A bug occured when the user load successively several datasets in
  Prostar. The previous ones were note correctly erased and this has
  lead to side effects. This bug is now fixed

NEW FEATURES

- Enhancement of the string-based filtering UI

- The automatic generation of an analysis report has been integrated in
  the Dataset Manager (menu 'Export'). It allows the user to download
  plots and parameters used in Prostar ont their dataset.

- Added a Gene Ontology (GO) analysis module in Data Processing. This
  module allows to perform GO grouping and GO Enrichment.

- Several plots are now based on the package highcharter which is a
  wrapper to the highcharts graphical library. It provides
  interactivity with the user.

[ProtGenerics](https://bioconductor.org/packages/ProtGenerics)
------------

Changes in version 1.9.1:

- add writeMSData <2017-09-20 Wed>

[PSEA](https://bioconductor.org/packages/PSEA)
----

Changes in version 1.11.1 (2017-06-09):

- A new vignette shows how we applied PSEA to deconvolute expression
  from RNA mixtures (originally presented in the Nature Methods paper).

[PureCN](https://bioconductor.org/packages/PureCN)
------

Changes in version 1.8.0:

NEW FEATURES

- Support for off-target reads in copy number normalization and
  segmentation

- Added mutation burden calculation

- More robust mapping bias estimation

- Added support for CNVkit coverage files (*.cnn, *.cnr)

- IntervalFile.R can annotate targets with gene symbols and
  automatically convert chromosome naming styles

- Better artifact filtering by using normalDB more efficiently

- Support for mappability scores

- Coverage calculation can now include duplicates

- calculateBamCoverageByInterval now provides fragment counts and
  duplication rates

- findBestNormal pooling now fragment count based, not coverage based

- Experimental support for GATK4

- predictSomatic now reports posterior probabilites of minor segment
  copy numbers, flags segments if copy numbers are unreliable

- Targets can be annotated with multiple gene symbols (comma separated)

- Code cleanups (switch to GRanges where possible, switch to optparse
  in command line tools)

API CHANGES

- Due to novel optimizations of provided bait intervals, we highly
  recommend to regenerate the interval files and normal databases and
  recalculate all coverages from BAM files

- New functions: annotateTargets, callMutationBurden

- Defunct functions: createSNPBlacklist, getDiploid, autoCurateResults,
  readCoverageGatk

- min.normals defaults to 2 (changed from 4) in setMappingBiasVcf

- normalDB.min.coverage defaults to 0.25 (changed from 0.2) in
  filterTargets

- log.ratio.calibration defaults to 0.1 (from 0.25) in runAbsoluteCN;
  now relative to purity, not log-ratio noise

- Removed gc.data from filterTargets since gc_bias is now added to
  tumor coverage

- dropped purecn.output from correctCoverageBias (no two-pass anymore)

- Coverage.R argument --gatkcoverage renamed to --coverage

- Dropped GC-normalization functionality in NormalDB, since this is now
  conveniently done in Coverage.R

- Renamed PureCN.R --outdir argument to --out. Can now specify a file
  prefix as in GATK. Filenames are thus not forced to sample id
  anymore.  If --out is a directory, it will behave like before and
  will use out/sampleid_suffix as filename.

[qcmetrics](https://bioconductor.org/packages/qcmetrics)
---------

Changes in version 1.15.2:

- Merge PR by aoles to re-use BiocStyle (see commit
  381d01d1c4dcced040f77a447962cd7865cb7417 for details) <2017-10-27
  Fri>

Changes in version 1.15.1:

- Don't user BiocStyle for vignette to allow to knit documents from
  within the vignette. This was needed to fix the error at build time.
  <2017-10-25 Wed>

Changes in version 1.15.0:

- Bioconductor devel

[qpgraph](https://bioconductor.org/packages/qpgraph)
-------

Changes in version 2.12:

USER VISIBLE CHANGES

- Updated citation data to include the published work in Castelo and
  Roverato (2017) on path weights.

[qsea](https://bioconductor.org/packages/qsea)
----

Changes in version 1.3.2:

New feature

- select active BSgenome masks for pattern density estimation

Bugfix

- adaptions to changed behavior of GenomicRanges

Changes in version 1.3.1:

- transition to Github

[QuartPAC](https://bioconductor.org/packages/QuartPAC)
--------

Changes in version 1.8.1:

- Bug fix to support new pdb file locations.

[Rbowtie2](https://bioconductor.org/packages/Rbowtie2)
--------

Changes in version 1.0.0:

INITIAL RELEASE

- bowtie2 version 2.3.2

- AdapterRemoval version 2.2.1a

[RCyjs](https://bioconductor.org/packages/RCyjs)
-----

Changes in version 1.9.8:

BUG FIXES

- Significant (3x) speedup.  A 5000-node, 6000-edge graph transmits to
  Cytoscape from R in about 20 seconds.

[ReactomePA](https://bioconductor.org/packages/ReactomePA)
----------

Changes in version 1.21.3:

- new project site using blogdown <2017-09-30, Fri>

Changes in version 1.21.2:

- correct typo in vignette <2017-05-18, Thu>

Changes in version 1.21.1:

- update vignette according to the change of reactome.db (pathway ID
  was changed) <2017-04-28, Fri>

[recount](https://bioconductor.org/packages/recount)
-------

Changes in version 1.3.13:

BUG FIXES

- Changed reproduce_ranges() since disjoint exons are more useful than
  reduced exons for downstream analyses.

Changes in version 1.3.12:

NEW FEATURES

- Added the function read_counts().

Changes in version 1.3.9:

SIGNIFICANT USER-VISIBLE CHANGES

- Added citations for
  https://www.biorxiv.org/content/early/2017/06/03/145656 and
  https://f1000research.com/articles/6-1558/v1 as well as mentions to
  them in the vignette.

Changes in version 1.3.7:

SIGNIFICANT USER-VISIBLE CHANGES

- add_predictions() was bumped to version 0.0.05

Changes in version 1.3.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette now uses the new BiocStyle::html_document that was recently
  released.

Changes in version 1.3.2:

NEW FEATURES

- coverage_matrix() now has two new arguments: scale and round. Use
  scale = FALSE to get raw coverage counts, which you can then scale
  with scale_counts(). scale is set to TRUE by default, so the counts
  are scaled to a library size of 40 million reads. round is set to
  FALSE by default, but can be set to TRUE if you want to get integer
  counts, just as in the default of scale_counts().

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Changed the default version argument of add_predictions() to
  'latest'. Internally, that's still 0.0.03.

[RedeR](https://bioconductor.org/packages/RedeR)
-----

Changes in version 1.26.0:

- Regular maintenance, stable release.

[regioneR](https://bioconductor.org/packages/regioneR)
--------

Changes in version 1.9.2:

NEW FEATURES

- Simplified the interface of toGRanges for simpler use when manually
  creating GRanges. Now toGRanges("chr1", 10, 20) is valid.

BUG FIXES

- Multiple minor bug fixes

[regionReport](https://bioconductor.org/packages/regionReport)
------------

Changes in version 1.11.6:

SIGNIFICANT USER-VISIBLE CHANGES

- Changed the default style to BiocStyle::html_document to mirror
  recent changes in BiocStyle.

Changes in version 1.11.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette now uses the new BiocStyle::html_document that was recently
  released.

Changes in version 1.11.2:

BUG FIXES

- Fixed the citations.

- Fixed DESeq2Report() for limma-results so that it will properly cite
  limma.

Changes in version 1.11.1:

SIGNIFICANT USER-VISIBLE CHANGES

- DESeq2Report() can now be used with other software if their results
  are made to look like DESeq2 results. For example, with limma-voom
  results.

BUG FIXES

- Made the DESeq2Report() more robust in case rlog() fails initially.

[rGREAT](https://bioconductor.org/packages/rGREAT)
------

Changes in version 1.9.1:

- change url according to changes on GREAT website

[Rhdf5lib](https://bioconductor.org/packages/Rhdf5lib)
--------

Changes in version 1.0:

New features

- Updated internal version of HDF5 to 1.8.19

Bug fixes

- Switched Windows compilation script to use CMake rather than
  configure. This seems to have solved problems in opening multiple
  files using different access modes.

Changes in version 0.99:

- Package submitted

[RnaSeqSampleSize](https://bioconductor.org/packages/RnaSeqSampleSize)
----------------

Changes in version 1.9.1 (2017-08-03):

BUG FIXES

- Fix a bug about estimating genes in pathway power in
  est_power_distribution function, which was caused by update of KEGG
  web site.

[RnBeads](https://bioconductor.org/packages/RnBeads)
-------

Changes in version 1.9.4:

- You can now retrieve custom region sets from the RnBeads resource
  using the rnb.load.annotation.from.db function

- improved import from GEO

- bigff disk dump is now the default option for dealing with disk-based
  big matrices

- Several bugfixes and performance improvements in Greedycut,
  normalization, Bisulfite data loading and others

- Added differential variability analysis in addition to differential
  methylation

Changes in version 1.9.3:

- LOLA support for differentially methylated regions

- various LOLA utitility functions

- What used to be the "differential.enrichment" option is now called
  "differential.enrichment.go" to reflect the distinction between GO
  and LOLA enrichment analyses

- Missing methylation data can now be imputed (using mean, random or
  nearest-neighbor imputation)

- Gender prediction now supports sequencing and EPIC data

- Added a function for combining two array-based datasets from
  (potentially) different platforms.

- Several bugfixes and performance improvements

[rols](https://bioconductor.org/packages/rols)
----

Changes in version 2.5.6:

- bump version

Changes in version 2.5.5:

- Fix failing test do to different order <2017-10-23 Mon>

Changes in version 2.5.4:

- Fix bug in pagesize limit when reporting children, ancestors, ...
  (see issue #25) <2017-10-17 Tue>

- New term(s) to data.frame coersion methods <2017-10-14 Sat>

Changes in version 2.5.3:

- Add import Ontology in NAMESPACE

Changes in version 2.5.2:

- Update to latest BiocStyle <2017-09-01 Fri>

- Quick fix for issue #24 (reported upsteams) <2017-09-02 Sat>

- Use Ontology generic from BiocGenerics <2017-09-03 Sun>

Changes in version 2.5.1:

- Fix unit test <2017-06-20 Tue>

Changes in version 2.5.0:

- Bioconductor devel 3.6

[ROTS](https://bioconductor.org/packages/ROTS)
----

Changes in version 1.6.0:

- Added support for testing multiple groups

- Bug fixes

[rpx](https://bioconductor.org/packages/rpx)
---

Changes in version 1.13.4:

- Restore reference unit test

Changes in version 1.13.3:

- No commit

Changes in version 1.13.2:

- Fix errors due to changes on PX side. PXD000001 reference remains
  untested - see https://twitter.com/lgatt0/status/885091284142239744
  for details. <2017-07-12 Wed>

Changes in version 1.13.1:

- Using xml2 <2017-06-13 Tue>

- Use functions instead of methods <2017-06-13 Tue>

Changes in version 1.13.0:

- Bioc devel 3.6

Changes in version 1.12.1:

- Using xml2 (backported from devel) <2017-06-13 Tue>

- Use functions instead of methods (backported from devel) <2017-06-13
  Tue>

[Rqc](https://bioconductor.org/packages/Rqc)
---

Changes in version 1.12:

BUG FIXES

- Fix bug when generating report on Windows system

- Using default graphic device for report generation

[rSFFreader](https://bioconductor.org/packages/rSFFreader)
----------

Version: 0.99.0
Category: Seattle time
Text:

Version: 0.99.0
Category: Marc will contact you off-tracker for the follow-up
Text:

Version: 0.99.0
Category: MS: If the version is set to 0.5, the package still goes to
        the release branch, or stays in devel
Text:

Version: 0.99.0
Category: If it goes to release, I'm cool with that, it is still under
        construction, but functional
Text:

[Rsubread](https://bioconductor.org/packages/Rsubread)
--------

Version: 1.28.0
Category: NEW FEATURES
Text:

Version: 1.28.0
Category: o New functions: promoterRegions() and txUnique
Text:

Version: 1.28.0
Category: o New parameter in featureCounts() for counting long reads -
        'isLongRead
Text:

Version: 1.28.0
Category: o New parameter in featureCounts() for counting reads by read
        groups - 'byReadGroup
Text:

Version: 1.28.0
Category: o New parameter in featureCounts() for requiring minimum
        fraction of overlapping bases in each feature -
        'fracOverlapFeature
Text:

Version: 1.28.0
Category: o Assignment results for each read can be added to the
        provided SAM/BAM files by featureCounts ('reportReads' option
Text:

Version: 1.28.0
Category: o Align() and subjunc() produce a mapping summary including
        percentages of uniquely mapped reads, multi-mapping reads and
        unmapped reads
Text:

[RTCA](https://bioconductor.org/packages/RTCA)
----

Changes in version 2009-07-13:

- combineRTCA(list): Additional column is renamed into Plate. The vlues
  is evaluated from list item names. When the list has no name, an
  integer index beginning from 1 is used. Special attentions to list
  partially with names is noted in the documentation.

- parseRTCA(file, dec=".",phenoData, skipWell,...): Example is added in
  the documentation how to import pre-configured phenoData. Details
  section in the documentation is re-written to describe the process of
  parsing.

- RTCA-class: Experiment ID added to RTCA class

- Makefile: add Makefile to simplify common tasks like check and
  install

- plotGridEffect: takes 'column' instead of 'col' as mode parameter,
  and renders the mode as the title of the legend. Documentation
  updated.

- plotRTCA: is removed from the package and is substituted by the plot
  function.

[RTN](https://bioconductor.org/packages/RTN)
---

Changes in version 1.16.0:

- Improved annotation slots for TNA/TNI classes.

[RTNduals](https://bioconductor.org/packages/RTNduals)
--------

Changes in version 1.2.0:

- Improved workflow integration with derivative packages (RTN /
  RTNsurvival).

[RTNsurvival](https://bioconductor.org/packages/RTNsurvival)
-----------

Changes in version 1.0.0:

- 1st Bioconductor release of RTNsurvival [2017-03-01].

[runibic](https://bioconductor.org/packages/runibic)
-------

Version: 0.99.14
Category: Updating tutorial, adding citation
Text:

Version: 0.99.13
Category: Documentation to the package was updated.
Text:

Version: 0.99.9
Category: Including required biocViews. Refactoring R code
Text:

Version: 0.99.8
Category: Code optimizations
Text:

Version: 0.99.6
Category: Parallel version using OpenMP
Text:

Version: 0.99.5
Category: Fixing discretize() function bug
Text:

Version: 0.99.2
Category: Improved performance of functions in package
Text:

Version: 0.99.2
Category: Removing some errors and warnings
Text:

Version: 0.99.1
Category: Compatibility with SummarizedExperiment
Text:

Version: 0.99.1
Category: Providing documentation and examples
Text:

Version: 0.99.0
Category: Development of the runibic package prototype completed
Text:

Version: 0.99.0
Category: First upload to Bioconductor
Text:

Version: 0.90.1
Category: Working version of UniBic algorithm
Text:

[S4Vectors](https://bioconductor.org/packages/S4Vectors)
---------

Changes in version 0.16.0:

NEW FEATURES

- Introduce FilterResults as generic parent of FilterMatrix.

- Optimized subsetting of an Rle object by an integer vector. Speed up
  is about 3x or more for big objects with respect to BioC 3.5.

SIGNIFICANT USER-VISIBLE CHANGES

- coerce,list,DataFrame generates "valid" names when list has none.
  This ends up introducing an inconsistency between DataFrame and
  data.frame but it is arguably a good one. We shouldn't rely on
  DataFrame() to generate variable names from scratch anyway.

BUG FIXES

- Fix showAsCell() on data-frame-like and array-like objects with a
  single column, and on SplitDataFrameList objects.

- Calling DataFrame() with explict 'row.names=NULL' should block
  rownames inference.

- cbind.DataFrame() ensures every argument is a DataFrame, not just
  first.

- rbind_mcols() now is robust to missing 'x'.

- Fix extractROWS() for arrays when subscript is a RangeNSBS.

- Temporary workaround to make the "union" method for Hits objects work
  even in the presence of another "union" generic in the cache (which
  is the case e.g. if the user loads the lubridate package).

- A couple of (long-time due) tweaks and fixes to "unlist" method for
  List objects so that it behaves consistently with "unlist" method for
  CompressedList objects.

- Modify Mini radix C code to accomodate a bug in Apple LLVM version
  6.1.0 optimizer. [commit 241150d2b043e8fcf6721005422891baff018586]

- Fix match,Pairs,Pairs() [commit
  a08c12bf4c31b7304d25122c411d882ec52b360c]

- Various other minor fixes.

[SC3](https://bioconductor.org/packages/SC3)
---

Changes in version 1.5.5 (2017-10-19):

- `scater` now works based on `SingleCellExperiment` class

- `sc3` slot has moved to the `metadata` slot of `SingleCellExperiment`
  class

- `scater` class is completely deprecated from `SC3`

- All `scater` functionality can work on `SingleCellExperiment` class

- Clustering logic hasn't been changed at all

[Scale4C](https://bioconductor.org/packages/Scale4C)
-------

Version: 1.0.0
Category: First release of Scale4C
Text:

[scater](https://bioconductor.org/packages/scater)
------

Changes in version 1.5.11:

- Complete refactoring of the package to use the SingleCellExperiment
  class

[scfind](https://bioconductor.org/packages/scfind)
------

Changes in version 0.99.0:

- The first version of scmap is submitted to Bioconductor.

[scmap](https://bioconductor.org/packages/scmap)
-----

Changes in version 0.99.2 (2017-07-05):

- Fixed issues for the first round of reviews at Bioconductor.

Changes in version 0.99.0 (2017-06-28):

- The first version of scmap is submitted to Bioconductor.

[scone](https://bioconductor.org/packages/scone)
-----

Changes in version 1.1.3 (2018-10-19):

- Modified COR scores to reflect R^2 for regression on UV/WV/QC rather
  than max correlation.

- Bug fix for shiny report.

- Updated vignette.

Changes in version 1.1.2:

- Moved many Imports: packages to Suggests: to reduce impact of loading
  scone

Changes in version 1.1.1:

- Added fast ziber imputation method

[scPipe](https://bioconductor.org/packages/scPipe)
------

Changes in version 0.99.20 (2017-09-22):

- scPipe now supports SingleCellExperiment class and use it as the base
  class

- add two functions `plot_demultiplex` and `plot_UMI_dup`

- scPipe support the offical bam tags for cell barcode and UMI

Changes in version 0.99.0 (2017-07-28):

- Package prepared for Bioconductor submission.

[scran](https://bioconductor.org/packages/scran)
-----

Changes in version 1.5.13:

- Supported parallelization in buildSNNGraph(), overlapExprs() with
  BPPARAM options.

- Forced zero-derived residuals to a constant value in
  correlatePairs(), overlapExprs().

- Allowed findMarkers() to return IUT p-values, to identify uniquely
  expressed genes in each cluster. Added options to specify the
  direction of the log-fold changes, to focus on upregulated genes in
  each cluster.

- Fixed bug in correlatePairs() when per.gene=TRUE and no spike-ins are
  available.  Added block.size argument to control caching.

- Switched all C++ code to use the beachmat API. Modified several
  functions to accept ANY matrix-like object, rather than only base
  matrix objects.

- quickCluster() with method="igraph" will now merge based on
  modularity to satisfy min.size requirements. Added max.size option to
  restrict the size of the output clusters.

- Updated the trendVar() interface with parametric, method arguments.
  Deprecated the trend="semiloess" option in favour of parametric=TRUE
  and method="loess". Modified the NLS equation to guarantee
  non-negative coefficients of the parametric trend.  Slightly modified
  the estimation of NLS starting parameters. Second d.f. of the fitted
  F-distribution is now reported as df2 in the output.

- Modified decomposeVar() to automatically use the second d.f. when
  test="f".

- Added option in denoisePCA() to return the number of components or
  the low-rank approximation. The proportion of variance explained is
  also stored as an attribute in all return values.

- Fixed a variety of bugs in mnnCorrect().

[SeqArray](https://bioconductor.org/packages/SeqArray)
--------

Changes in version 1.18.0:

NEW FEATURES

- progress information: showing overall running time when completed

- new variable names "$ref" and "$alt" can be used in `seqGetData()`
  and `seqBlockApply()`

- new argument '.progress' in `seqDigest()`

- new argument 'ref.allele' in `seqAlleleCount()`

- new variable name "$chrom_pos_allele" can be used in `seqGetData()`
  and `seqBlockApply()`

UTILITIES

- move VariantAnnotation to the suggest field from the import field

- remove an unused argument '.list_dup' in `seqBlockApply()`

- slightly improve the computational efficiency of `seqAlleleFreq()`
  and `seqAlleleCount()` when 'ref.allele=0'

- `seqGetData(f, "$chrom_pos")` outputs characters with the format
  'chromosome:position' instead of 'chromosome_position'

BUG FIXES

- fix the unexpected behaviors in `seqSetFilter(, action="push")` and
  `seqSetFilter(, action="push+intersect")`

- fix a bug in `seqGetData(f, "$dosage")` when the number of unique
  alleles at a site greater than 3
  (https://github.com/zhengxwen/SeqArray/issues/21)

- fix a bug in `seqSNP2GDS()` for inverted genotypes during importing
  data from SNP GDS files
  (https://github.com/zhengxwen/SeqArray/issues/22)

- fix an issue of no phase data in `seqExport()`

[seqCAT](https://bioconductor.org/packages/seqCAT)
------

Changes in version 1.0.0:

FEATURES

- Create single nucleotide variant (SNV) profiles from RNA/DNA-seq
  samples

- Characterise the biological equivalency and difference between
  samples

- Evaluate putative impacts of SNVs differing between samples

- Investigate and validate known variants and specific genomic regions

- Authenticate cell lines with a known SNV profile or the COSMIC
  database

[seqcombo](https://bioconductor.org/packages/seqcombo)
--------

Changes in version 0.99.11:

- geom_genotype <2017-08-29, Tue>

Changes in version 0.99.10:

- geom_hybrid <2017-08-17, Thu>

Changes in version 0.99.9:

- hybrid_plot <2017-06-30, Fri>

Changes in version 0.0.3:

- more parameters for plot, by, xlab, color, fill etc.

Changes in version 0.0.2:

- add vignette

Changes in version 0.0.1:

- initial version with plot method for nucleotide differences between
  twwo aligned sequences <2016-11-16, Wed>

[SeqSQC](https://bioconductor.org/packages/SeqSQC)
------

Version: 0.99.10
Text: -- updated the man for `sampleQC` and `LoadVfile`.

Version: 0.99.9
Text: Updated the `output` argument in function of `sampleQC` and
        `LoadVfile`, by using the default value of `sampleqc` in
        working directory. In vignette and the example code of these 2
        functions, temporary directory is used for the directory to
        save the QC results.

Version: 0.99.8
Text: vignette, using a tempdir() for vignette output.

Version: 0.99.8
Text: Vignette: Additional explanation added for the option to do a
        wrapper all-inclusive QC or specific QC steps.

Version: 0.99.8
Text: documentation for argument "output" in "LoadVfile" and
        "sampleQC", to indicate that the `dirname(output)` would be
        used as the directory to save the other QC results in `.txt`
        and plots in `.pdf`. The default value is added as
        `file.path(tempdir(), "sampleqc")`.

Version: 0.99.8
Text: sampleQC.R: Arguments of `plotting=TRUE, results=TRUE` added and
        documented, giving the option for users to whether or not
        output the QC results and plots. The problem list will always
        in output for the simplest sample QC result and summary.

Version: 0.99.8
Text: the package documentation updated, by adding a `@seealso` for
        links of main functions of SeqSQC.

Version: 0.99.7
Text: Added checking for input in subsetGDS.

Version: 0.99.7
Text: Removed commented code for subsetGDS.

Version: 0.99.7
Text: SeqSQCclass renamed as SeqSQC.

Version: 0.99.7
Text: mergeGDS.R. using lapply for some repeated functions like
        "read.gdsn(index.gdsn())". In SexCheck, Inbreeding, IBDCheck,
        PCACheck.

Version: 0.99.7
Text: sampleQC: 1. check the class of vfile for SeqSQC file or
        vcf/plink file, instead of using vfile / sfile.

Version: 0.99.7
Text: plotting function writing separately (not export), wrapper with
        plotQC().

Version: 0.99.7
Text: Removed IBDRemoveAll.R. Will check IBD for all sample pairs
        (including benchmark samples) by adding argument "all" in
        "IBDRemove.R".

Version: 0.99.7
Text: debugged warning messages from "LoadVfile" and using "rbokeh".

Version: 0.99.6
Text: modified SeqSQCclass to SeqSQC (not in vignette..).

Version: 0.99.5
Text: Installed BiocGenerics locally.

Version: 0.99.4
Category: Load benchmark data from ExperimentHub
Text:

Version: 0.99.4
Category: o added in the R script, will download only once with first
        time running of the LoadVfile() or sampleQC
Text:

Version: 0.99.4
Category: Bioconductor submission check:
Text:

Version: 0.99.4
Category: o Added unit test
Text:

Version: 0.99.4
Category: o Added a NEWS file to keep track of changes
Text:

Version: 0.99.4
Category: o Added zzz.R to fix the no visible binding for global
        functions or variables.
Text:

Version: 0.99.4
Category: o Added the "example_sub.vcf" for 1000 lines of variants to
        run as example in the package vignette
Text:

Version: 0.99.4
Category: o Added accessor methods for SeqSQCclass data structure to
        get the slots of "gdsfile" and "QCresult
Text:

Version: 0.99.4
Category: Vignettes:
Text:

Version: 0.99.4
Category: o Added bioconductor installation and library load section in
        the vignette.
Text:

Version: 0.99.4
Category: o Added runnable example vcf file added in
        "inst/extdata/example_sub.vcf", with 1000 lines of variants.
Text:

Version: 0.99.4
Category: MAN:
Text:

Version: 0.99.4
Category: o added package documentation for dataset, class, methods and
        constructor
Text:

[SeqVarTools](https://bioconductor.org/packages/SeqVarTools)
-----------

Changes in version 1.15.3:

- alleleFrequency method accounts for sex when computing frequency for
  X and Y chromosomes.

Changes in version 1.15.2:

- Added iterator classes: SeqVarBlockIterator, SeqVarRangeIterator,
  SeqVarWindowIterator, SeqVarListIterator.

- Creating a SeqVarData object with missing sample or variant
  annotation will store 0-column data frames in sampleData or
  variantData, instead of duplicating sample.id and variant.id.

- Added methods to return variant data in expanded form, with one row
  per alternate allele.

Changes in version 1.15.1:

- Following SeqArray, remove dependency on VariantAnnotation

- Add generic isSNV (replacing previous import of this generic from
  VariantAnnotation)

[ShortRead](https://bioconductor.org/packages/ShortRead)
---------

Changes in version 1.35:

SIGNIFICANT USER-VISIBLE CHANGES

- Reads up to 2M bases can be parsed

[SIMLR](https://bioconductor.org/packages/SIMLR)
-----

Changes in version 1.3.1 (2017-09-15):

- Implemented estimation of number of clusters from data

Changes in version 1.2.1 (2017-05-29):

- Fix Windows memory usage

[SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment)
--------------------

Changes in version 0.99.4:

- New package SingleCellExperiment, for representation of single-cell
  genomics data.

[slalom](https://bioconductor.org/packages/slalom)
------

Changes in version 0.99.0:

- R/C++ package implementing f-scLVM model submitted to Bioconductor

[SNPRelate](https://bioconductor.org/packages/SNPRelate)
---------

Changes in version 1.12.0:

- new arguments 'with.sample.id' and 'with.snp.id' in
  `snpgdsSNPRateFreq()`

Changes in version 1.10.2:

- progress information: showing overall running time when completed

- An unexpected exception in a thread could result in deadlock: the
  current implementation shows error information and exits the R
  session

[splatter](https://bioconductor.org/packages/splatter)
--------

Changes in version 1.1.8 (2017-10-13):

- Now published in Genome Biology!

- Converted to the SingleCellExperiment object

- Added new simulations: BASiCS, mfa, PhenoPath, ZINB-WaVE

- Added batch effects to the Splat simulation. This required a change
  to the SplatParams object.

- Improved scDD estimation

- Added and improved comparison functions

- Improved default Splat parameters and estimation

- Improvements to the Lun2Params object

- Added addGeneLength function

- Updated simulation references

- Various other minor updates and bug fixes

[SPONGE](https://bioconductor.org/packages/SPONGE)
------

Changes in version 0.99.1:

- Support for ExpressionSet

- Extensive unit testing

- Simplified dependencies

- Code improvements

- Methods return data frames by default

- Default log level is now ERROR

[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
--------------------

Version: 1.8.0
Category: NEW FEATURES
Text: Add 'chunk_dim' and 'level' arguments to
        saveHDF5SummarizedExperiment().

Version: 1.8.0
Category: NEW FEATURES
Text: Add coercion from ExpressionSet to SummarizedExperiment.

Version: 1.8.0
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text:

Version: 1.8.0
Category: DEPRECATED AND DEFUNCT
Text: Remove 'force' argument from seqinfo() and seqlevels() setters
        (the argument got deprecated in BioC 3.5 in favor of new and
        more flexible 'pruning.mode' argument).

Version: 1.8.0
Category: BUG FIXES
Text: Coercion from SummarizedExperiment to RangedSummarizedExperiment
        was losing the metadata columns. Fixed now.

Version: 1.8.0
Category: BUG FIXES
Text: Fix cbind() and rbind() of SummarizedExperiment objects when some
        of the assays are DataFrame or data.frame objects.

Version: 1.8.0
Category: BUG FIXES
Text: '$' completion on SummarizedExperiment works in RStudio and on
        RangedSummarizedExperiment.

[SWATH2stats](https://bioconductor.org/packages/SWATH2stats)
-----------

Changes in version 1.7.7:

BUG FIXES

- update vignette to prevent misunderstanding what data to load

Changes in version 1.7.6:

BUG FIXES

- manual convert4PECA

Changes in version 1.7.5:

BUG FIXES

- library not recognized in vignette

Changes in version 1.7.4:

NEW FEATURES

- add function convert4PECA

Changes in version 1.7.3:

NEW FEATURES

- add option to keep all columns to disaggregate function

Changes in version 1.7.2:

NEW FEATURES

- add option to no not remove decoys to filter_proteotypic_peptides(),
  filter_on_max_peptides(), and filter_on_min_peptides()

Changes in version 1.7.1:

NEW FEATURES

- SWATH2stats in BioC 3.6 development release

Changes in version 1.6.1:

NEW FEATURES

- SWATH2stats in BioC 3.5 release

[TarSeqQC](https://bioconductor.org/packages/TarSeqQC)
--------

Changes in version 1.7.1:

- CITATION file was added to include the new publication describing
  TarSeqQC application on real datasets

[TFARM](https://bioconductor.org/packages/TFARM)
-----

Changes in version 0.99.8:

MAJOR CORRECTIONS

- two datasets were renamed (from I to TF_Imp and from p to p_TFs)

- some loops were vectorized in the vignettes

[TFBSTools](https://bioconductor.org/packages/TFBSTools)
---------

Changes in version 3.5:

NEW FEATURES

- Add function to parse MEME output.

- Add parallel computing of searchSeq.

Changes in version 3.4:

BUG FIXES

- Fix a bug in PFMSimilarity.

- Fix an error when there are multiple classes for motif matrx.

NEW FEATURES

- readJASPARMatrix function to add the JASPAR format file.

Changes in version 3.3:

BUG FIXES

- Adapt the runMEME to work with meme 4.10.x version.

- Fix the scientific notation in run_MEME

- Better error handling of MEME wrappe

[tofsims](https://bioconductor.org/packages/tofsims)
-------

Version: 099.1
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text: changed function behvaiour in the whole package from call-by-ref
        to call-by value. Adjusted accordingly all examples and the
        vignette.

Version: 099.1
Category: INTERNALS
Text: depends now on ProtGenerics from which it uses 'mz'

Version: 099.1
Category: INTERNALS
Text: exchanged various print() with message()

Version: 099.1
Category: BUGFIXES
Text:

[topdownr](https://bioconductor.org/packages/topdownr)
--------

Changes in version 0.99.0 (2017-10-05):

- First public release.

[TPP](https://bioconductor.org/packages/TPP)
---

Changes in version 3.5.0:

- New Bioconductor release candidate! Chenges in version 3.5.xx
  (May-xx-2017):

- Importing the whole tidyverse instead of ggplot2, dplyr etc.
  individually

- Shift tidyverse and Biobase to "depends" so that they are
  automatically availabe for downstream analyses.

[trackViewer](https://bioconductor.org/packages/trackViewer)
-----------

Changes in version 1.13.11:

- fix the bug maxgap could not be less than 0 in findOverlaps

Changes in version 1.13.10:

- fix the bug documentation of ideogramPlot is different from the
  codes.

Changes in version 1.13.9:

- fix the bug when using findOverlaps, when 'type' is "any", at least
  one of 'maxgap' and 'minoverlap' must be set to its default value

Changes in version 1.13.8:

- automatic set the font size in ideogramPlot

Changes in version 1.13.7:

- fix the bug in ideogramPlot when data range is identical for heatmap
  layout.

Changes in version 1.13.6:

- add draggable tracks in browseTracks function

Changes in version 1.13.5:

- add arrow tools in browseTracks function

Changes in version 1.13.4:

- add label tools in browseTracks function

Changes in version 1.13.3:

- fix bugs in browseTracks function

Changes in version 1.13.2:

- add more columns to GRoperator function

- add browseTracks function

[transcriptogramer](https://bioconductor.org/packages/transcriptogramer)
-----------------

Changes in version 1.0.0:

- Release of the transcriptogramer package at Bioconductor
  [YYYY-MM-DD].

[treeio](https://bioconductor.org/packages/treeio)
------

Changes in version 1.1.2:

- new project site using blogdown <2017-09-28, Thu>

Changes in version 1.1.1:

- parse mlc file without dNdS <2017-08-31, Thu> +
  https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!topic/bioc-ggtree/hTRj-uldgAg

- better implementation of merge_tree <2017-08-31, Thu>

[TRONCO](https://bioconductor.org/packages/TRONCO)
------

Changes in version 2.8.1:

- Minor fix on documentation

[TVTB](https://bioconductor.org/packages/TVTB)
----

Changes in version 1.3.2 (2017-09-02):

Bug fix

- Updated vignette output format for compatibility with BiocStyle (>=
  2.5.19); one vignette turned to PDF to respect the package size limit
  in the new context.

Changes in version 1.3.1 (2017-05-30):

Bug fix

- Fixed VCF import from multiple files in _Shiny_ application.

[VariantFiltering](https://bioconductor.org/packages/VariantFiltering)
----------------

Changes in version 1.14:

USER VISIBLE CHANGES

- The 'MafDb' class and methods have been moved to the 'GenomicScores'
  package.

- 'VariantFilteringResults' objects also store now a
  'GenomeDescription' object that can be fetched using the method
  'referenceGenome()'.

- Changes in the internal storage of 'VariantFilteringResults' objects
  that allow one to serialize this type of object with 'saveRDS()'.

[wiggleplotr](https://bioconductor.org/packages/wiggleplotr)
-----------

Changes in version 1.1.2:

- New makeManhattanPlot() function makes it easier visualise
  association p-values together with read coverage and transcript
  structure plots.

[xcms](https://bioconductor.org/packages/xcms)
----

Changes in version 2.99.10:

BUG FIXES

- Fix #230: Failing vignettes on Windows.

Changes in version 2.99.9:

USER VISIBLE CHANGES

- Chromatographic peak detection uses adjusted retention times on an
  aligned XCMSnExp object (issue #213, #208).

- New parameter msLevel for processHistory,XCMSnExp.

- New parameter keepAdjustedRtime for filterMsLevel,XCMSnExp,
  dropChromPeaks, XCMSnExp and dropFeatureDefinitions,XCMSnExp.

- Add parameter msLevel to chromatogram,XCMSnExp method (issue #205).

- Obiwarp alignment is now performed on one MS level and adjustment is
  applied to all MS levels (issue #214).

- Add function plotMsData to plot intensity against retention time and
  m/z against retention time for a MS slice in one sample.

- Add argument msLevel = 1L to extractMsData method (issue #223).

- New applyAdjustedRtime function to consolidate the alignment results,
  i.e. replace the raw retention times in the XCMSnExp with the
  adjusted retention times.

- [,XCMSnExp method gains argument keepAdjustedRtime to allow keeping
  adjusted retention times in the sub-setting.

- Implement spectrapply,XCMSnExp to ensure returned results use
  adjusted retention times (if present).

- [[,XCMSnExp method returns a Spectrum object with adjusted retention
  time, if the XCMSnExp contains adjusted retention times.

- Argument 'sampleGroups' is mandatory for 'PeakDensityParam' (issue
  #228).

BUG FIXES

- Fix #191: Excessive memory use in fillPeaks.

- Fix #220: peaks matrix is missing column "sample" if no peaks were
  found in the first sample.

- Fix #222: findChromPeaks does not return an XCMSnExp object filtered
  to a single MS level despite peak detection is performed on a single
  level.

- Fix problem in plotMsData causing wrong colors to be used to label
  the data points.

Changes in version 2.99.8:

BUG FIXES

- Replace xcmsMSn Rnw with Rmd vignette to fix Windows build errors.

Changes in version 2.99.7:

BUG FIXES

- Fix #201: Warnings: 'readMSData2' is deprecated, thanks to L. Gatto.

- Merge with BioC git after transition

Changes in version 2.99.6:

NEW FEATURES

- calibrate,XCMSnExp method that allows to calibrate chromatographic
  peaks.

USER VISIBLE CHANGES

- Export phenoDataFromPaths function (issue $195).

- Add arguments mz and rt to featureDefinitions method allowing to
  extract features within the specified ranges.

- Increase n for the density function call in group density-based
  correspondence by 2.

- Replace xcmsDirect.Rnw with rmarkdown-based vignette using the new
  user interface.

BUG FIXES

- issue #196: removed the unnecessary requirement for same-dimension
  profile matrices in adjustRtime,XCMSnExp,ObiwarpParam.

- issue #194: fixes in retcor.obiwarp: 1) subset raw data if scanrange
  != NULL. 2) if the mz range of the two files to be aligned differ,
  expand them correctly. Depending on the profStep and the mz
  values/ranges the matrices were not expanded correctly.

- Potential problems in the plotChromPeakDensity function.

Changes in version 2.99.5:

USER VISIBLE CHANGES

- Re-enable sleep parameter in findPeaks.centWave and
  findPeaks.matchedFilter.

Changes in version 2.99.4:

NEW FEATURES

- Add plotChromPeaks function to plot the definition (rt and mz range)
  of detected chromatographic peaks of one file into the mz-rt plane.

- Add plotChromPeakImage function to plot the number of detected peaks
  along the retention time axis per file as an image plot.

USER VISIBLE CHANGES

- Move Chromatogram class and functionality to the MSnbase package

- Add argument msLevel to the findChromPeaks method to allow
  (chromatographic) peak detection also on MS level > 1.

BUG FIXES

- Polarity information was not read from mzXML files (issue #192).

Changes in version 2.99.3:

BUG FIXES

- issue #188: determine file type from file content if file ending not
  known.

Changes in version 2.99.2:

BUG FIXES

- issue #181: problem when isCentroided,Spectrum method returns NA
  because of too few peaks in a spectrum. Fixed by checking in such
  cases all spectra in the file.

- issue #184: add parameter sleep to do_groupChromPeaks_density
  function to be backwards compatible with the old group.density code.

Changes in version 2.99.1:

NEW FEATURES

- extractMsData to extract raw MS data as a data.frame (issue #120).

BUG FIXES

- issue #175: an error is now thrown if no peak group was identified
  for peak group retention time correction.

- issue #178: scanrange was collapsed when the adjusted range was
  reported (pull request by Jan Stanstrup).

- issue #180: error when both parameters method and smooth are provided
  in the retcor method.

Changes in version 2.99.0:

NEW FEATURES

- plotChromatogram and highlightChromPeaks functions.

- plotChromPeakDensity function.

- clean method for Chromatogram classes.

USER VISIBLE CHANGES

- Change default for ppm parameter in chromPeaks method to 0.

- extractChromatograms supports extraction of multiple rt and mz
  ranges.

- New parameter missing for extractChromatograms allowing to specify
  the intensity value to be used for rts for which no signal is
  available within the mz range.

- extractChromatograms returns Chromatograms of length equal to the
  number of scans within the specified rt range, even if no signals are
  measured (intensity values are NA).

Changes in version 1.53.1:

BUG FIXES

- Increase parameter n for the density call in the peak density
  correspondence method. This enables to separate neighboring peaks
  using small n (issue #161). Thanks to Jan Stanstrup.

[xps](https://bioconductor.org/packages/xps)
---

Changes in version 3.4:

VERSION xps-1.37.2

- configure.in file - unix line endings

VERSION xps-1.37.1

- update INSTALL and README file

[yamss](https://bioconductor.org/packages/yamss)
-----

Changes in version 1.3:

- Updating citation information.

- Changes related to moving Bioconductor from SVN to GIT.

[zFPKM](https://bioconductor.org/packages/zFPKM)
-----

Changes in version 0.99.18:

- initial release

- Methods adapted from primary manuscript

[zinbwave](https://bioconductor.org/packages/zinbwave)
--------

Changes in version 0.99.10 (2017-10-23):

- Added AIC and BIC to decide number of factors

- Added function to compute observational weights for DE

Changes in version 0.99.7 (2017-07-17):

- `zinbwave()` now returns a `SingleCellExperiment` object.

Changes in version 0.99.6 (2017-07-05):

- Fixed bug in zinb.loglik.matrix to avoid Inf values

Changes in version 0.99.5 (2017-07-03):

- Added function computeDevianceResiduals() to compute residuals

- Added function imputeZeros() to use the model to impute technical
  zeros

- Added function zinbwave() to perform dimensionality reduction

- Changed vignette to illustrate new functions

Changes in version 0.99.4 (2017-06-08):

- Switch from parallel to BiocParallel

- Improvements to code efficiency, e.g., avoid copying ZinbModel
  objects

Changes in version 0.99.3 (2017-05-31):

- More informative Description: field

- Improved documentation

- Added getAlpha, getBeta, and getGamma accessor functions

- Improved show() method

- Vectorized code in zinbSim()

- Fixed bug that introduced an error when initializing an object with
  empty X or V

- Added tests on numerical correctness

Changes in version 0.99.2 (2017-05-12):

- New formula interface for SummarizedExperiment

- Add t-SNE example to vignette

Changes in version 0.99.0 (2017-05-09):

- Bumped version for submission to Bioconductor

- Minor changes to compile vignette

- Required R version 3.4

Changes in version 0.1.4 (2017-04-10):

- Improved documentation

- zinbSim now produces matrix of J x n dimension

- Removed unnecessary dependency on clusterExperiment

- Additional tests

- Better default for epsilon

Changes in version 0.1.2 (2017-03-29):

- New vignette

- Change name to zinbwave

Changes in version 0.1.0 (2016-09-23):

- Introducing S4 class zinbModel

- Major restructuring of the initialization and optimization (see
  vignette)

- Method zinbFit to fit a model

- Many other methods, including zinbInitialize, zinbOptimize, zinbSim


Deprecated and Defunct Packages
===============================

Nine software packages were removed from this release (after being
deprecated in BioC 3.5): pdmclass, coRNAi, GENE.E, mmnet, AtlasRDF,
CopyNumber450k, saps, MeSHSim, GEOsearch.

Nine software packages (BioMedR, ddgraph, EWCE, HCsnip, stepwiseCM,
domainsignatures, iontree, oneChannelGUI, RCytoscape ) are deprecated
in this release and will be removed in BioC 3.7.

One experimental data package (CopyNumber450kData) was removed from this
release (after being deprecated in BioC 3.5).

No experimental data packages are deprecated in this release.

