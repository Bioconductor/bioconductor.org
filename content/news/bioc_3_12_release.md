October 28, 2020

Bioconductors:

We are pleased to announce Bioconductor 3.12, consisting of 
1939 software packages, 397 experiment data packages, 967 annotation
packages, and 28 workflows.

There are 124 new software packages, 9 new data experiment packages,
2 new annotation packages, 1 new workflow, and many updates and improvements
to existing packages; Bioconductor 3.12 is compatible with R 4.0.0,
and is supported on Linux, 32- and 64-bit Windows, and macOS 10.14.6 Mojave
or higher.  This release will include an updated Bioconductor [Amazon Machine
Image][1] and [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3]
for details and downloads.

[1]: /help/bioconductor-cloud-ami/
[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.11](#getting-started-with-bioconductor-311)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [NEWS from new and existing software packages](#news-from-new-and-existing-software-packages)
* [NEWS from new and existing data experiment packages](#news-from-new-and-existing-data-experiment-packages)
* [NEWS from new and existing workflows](#news-from-new-and-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.12
======================================

To update to or install Bioconductor 3.12:

1. Install R 4.0.0. Bioconductor 3.12 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).

New Software Packages
=====================

There are 124 new software packages in this release of Bioconductor.

- [ADImpute](/packages/ADImpute) Single-cell RNA sequencing
  (scRNA-seq) methods are typically unable to quantify the expression
  levels of all genes in a cell, creating a need for the
  computational prediction of missing values (‘dropout imputation’).
  Most existing dropout imputation methods are limited in the sense
  that they exclusively use the scRNA-seq dataset at hand and do not
  exploit external gene-gene relationship information. Here we
  propose two novel methods: a gene regulatory network-based approach
  using gene-gene relationships learnt from external data and a
  baseline approach corresponding to a sample-wide average. ADImpute
  can implement these novel methods and also combine them with
  existing imputation methods (currently supported: DrImpute, SAVER).
  ADImpute can learn the best performing method per gene and combine
  the results from different methods into an ensemble.

- [aggregateBioVar](/packages/aggregateBioVar) For single cell
  RNA-seq data collected from more than one subject (e.g. biological
  sample or technical replicates), this package contains tools to
  summarize single cell gene expression profiles at the level of
  subject. A SingleCellExperiment object is taken as input and
  converted to a list of SummarizedExperiment objects, where each
  list element corresponds to an assigned cell type. The
  SummarizedExperiment objects contain aggregate gene-by-subject
  count matrices and inter-subject column metadata for individual
  subjects that can be processed using downstream bulk RNA-seq tools.

- [AlpsNMR](/packages/AlpsNMR) Reads Bruker NMR data directories both
  zipped and unzipped. It provides automated and efficient signal
  processing for untargeted NMR metabolomics. It is able to
  interpolate the samples, detect outliers, exclude regions,
  normalize, detect peaks, align the spectra, integrate peaks, manage
  metadata and visualize the spectra. After spectra proccessing, it
  can apply multivariate analysis on extracted data. Efficient
  plotting with 1-D data is also available. Basic reading of 1D
  ACD/Labs exported JDX samples is also available.

- [ANCOMBC](/packages/ANCOMBC) ANCOMBC is a package for normalizing
  the microbial absolute abundance data due to unequal sampling
  fractions across samples, and identifying taxa (e.g. phyla,
  families, genera, species, etc.) that are differentially abundant
  with respect to the covariate of interest (e.g. study groups)
  between two or more groups of multiple samples.

- [AnVILBilling](/packages/AnVILBilling) AnVILBilling helps monitor
  AnVIL-related costs in R, using queries to a BigQuery table to
  which costs are exported daily.  Functions are defined to help
  categorize tasks and associated expenditures, and to visualize and
  explore expense profiles over time. This package will be expanded
  to help users estimate costs for specific task sets.

- [AnVILPublish](/packages/AnVILPublish) Use this package to create
  or update AnVIL workspaces from resources such as R / Bioconductor
  packages. The metadata about the package (e.g., select information
  from the package DESCRIPTION file and from vignette YAML headings)
  are used to populate the 'DASHBOARD'. Vignettes are translated to
  python notebooks ready for evaluation in AnVIL.

- [bambu](/packages/bambu) bambu is a R package for multi-sample
  transcript discovery and quantification using long read RNA-Seq
  data. You can use bambu after read alignment to obtain expression
  estimates for known and novel transcripts and genes. The output
  from bambu can directly be used for visualisation and downstream
  analysis such as differential gene expression or transcript usage.

- [BayesSpace](/packages/BayesSpace) Tools for clustering and
  enhancing the resolution of spatial gene expression experiments.
  BayesSpace clusters a low-dimensional representation of the gene
  expression matrix, incorporating a spatial prior to encourage
  neighboring spots to cluster together. The method can enhance the
  resolution of the low-dimensional representation into "sub-spots",
  for which features such as gene expression or cell type composition
  can be imputed.

- [BiocIO](/packages/BiocIO) Implements `import()` and `export()`
  standard generics for importing and exporting biological data
  formats. `import()` supports whole-file as well as chunk-wise
  iterative import. The `import()` interface optionally provides a
  standard mechanism for 'lazy' access via `filter()` (on row or
  element-like components of the file resource), `select()` (on
  column-like components of the file resource) and `collect()`. The
  `import()` interface optionally provides transparent access to
  remote (e.g. via https) as well as local access. Developers can
  register a file extension, e.g., `.loom` for dispatch from
  character-based URIs to specific `import()` / `export()` methods
  based on classes representing file types, e.g., `LoomFile()`.

- [biocthis](/packages/biocthis) This package expands the usethis
  package with the goal of helping automate the process of creating R
  packages for Bioconductor or making them Bioconductor-friendly.

- [bluster](/packages/bluster) Wraps common clustering algorithms in
  an easily extended S4 framework. Backends are implemented for
  hierarchical, k-means and graph-based clustering. Several utilities
  are also provided to compare and evaluate clustering results.

- [BrainSABER](/packages/BrainSABER) The Allen Institute for Brain
  Science provides an RNA sequencing (RNA-Seq) data resource for
  studying transcriptional mechanisms involved in human brain
  development known as BrainSpan. BrainSABER is an R package that
  facilitates comparison of user data with the various developmental
  stages and brain structures found in the BrainSpan atlas by
  generating dynamic similarity heatmaps for the two data sets. It
  also provides a self-validating container for user data.

- [CellaRepertorium](/packages/CellaRepertorium) Methods to cluster
  and analyze high-throughput single cell immune cell repertoires,
  especially from the 10X Genomics VDJ solution. Contains an R
  interface to CD-HIT (Li and Godzik 2006). Methods to visualize and
  analyze paired heavy-light chain data. Tests for specific
  expansion, as well as omnibus oligoclonality under hypergeometric
  models.

- [cfDNAPro](/packages/cfDNAPro) cfDNA fragment size metrics are
  important features for utilizing liquid biopsy in tumor early
  detection, diagnosis, therapy personlization and monitoring.
  Analyzing and visualizing insert size metrics could be time
  intensive. This package intends to simplify this exploration
  process, and it offers two sets of functions for data
  characterization and data visualization.

- [ChromSCape](/packages/ChromSCape) ChromSCape - Chromatin landscape
  profiling for Single Cells - is a ready-to-launch user-friendly
  Shiny Application for the analysis of single-cell epigenomics
  datasets (scChIP-seq, scATAC-seq, scCUT&Tag, ...) from aligned data
  to differential analysis & gene set enrichment analysis. It is
  highly interactive, enables users to save their analysis and covers
  a wide range of analytical steps: QC, preprocessing, filtering,
  batch correction, dimensionality reduction, vizualisation,
  clustering, differential analysis and gene set analysis.

- [corral](/packages/corral) Correspondence analysis (CA) is a matrix
  factorization method, and is similar to principal components
  analysis (PCA). Whereas PCA is designed for application to
  continuous, approximately normally distributed data, CA is
  appropriate for non-negative, count-based data that are in the same
  additive scale. The corral package implements CA for dimensionality
  reduction of a single matrix of single-cell data, as well as a
  multi-table adaptation of CA that leverages data-optimized scaling
  to align data generated from different sequencing platforms by
  projecting into a shared latent space. corral utilizes sparse
  matrices and a fast implementation of SVD, and can be called
  directly on Bioconductor objects (e.g., SingleCellExperiment) for
  easy pipeline integration. The package also includes the option to
  apply CA-style processing to continuous data (e.g., proteomic TOF
  intensities) with the Hellinger distance adaptation of CA.

- [customCMPdb](/packages/customCMPdb) This package serves as a query
  interface for important community collections of small molecules,
  while also allowing users to include custom compound collections.

- [CytoTree](/packages/CytoTree) A trajectory inference toolkit for
  flow and mass cytometry data. CytoTree is a valuable tool to build
  a tree-shaped trajectory using flow and mass cytometry data. The
  application of CytoTree ranges from clustering and dimensionality
  reduction to trajectory reconstruction and pseudotime estimation.
  It offers complete analyzing workflow for flow and mass cytometry
  data.

- [dasper](/packages/dasper) The aim of dasper is to detect aberrant
  splicing events from RNA-seq data. dasper will use as input both
  junction and coverage data from RNA-seq to calculate the deviation
  of each splicing event in a patient from a set of user-defined
  controls. dasper uses an unsupervised outlier detection algorithm
  to score each splicing event in the patient with an outlier score
  representing the degree to which that splicing event looks
  abnormal.

- [DegNorm](/packages/DegNorm) This package performs degradation
  normalization in bulk RNA-seq data to improve differential
  expression analysis accuracy.

- [densvis](/packages/densvis) Implements the density-preserving
  modification to t-SNE and UMAP described by Narayan et al. (2020)
  <doi:10.1101/2020.05.12.077776>. The non-linear dimensionality
  reduction techniques t-SNE and UMAP enable users to summarise
  complex high-dimensional sequencing data such as single cell RNAseq
  using lower dimensional representations. These lower dimensional
  representations enable the visualisation of discrete
  transcriptional states, as well as continuous trajectory (for
  example, in early development). However, these methods focus on the
  local neighbourhood structure of the data. In some cases, this
  results in misleading visualisations, where the density of cells in
  the low-dimensional embedding does not represent the
  transcriptional heterogeneity of data in the original
  high-dimensional space. den-SNE and densMAP aim to enable more
  accurate visual interpretation of high-dimensional datasets by
  producing lower-dimensional embeddings that accurately represent
  the heterogeneity of the original high-dimensional space, enabling
  the identification of homogeneous and heterogeneous cell states.
  This accuracy is accomplished by including in the optimisation
  process a term which considers the local density of points in the
  original high-dimensional space. This can help to create
  visualisations that are more representative of heterogeneity in the
  original high-dimensional space.

- [escape](/packages/escape) A bridging R package to facilitate gene
  set enrichment analysis (GSEA) in the context of single-cell RNA
  sequencing. Using raw count information, Seurat objects, or
  SingleCellExperiment format, users can perform and visualize GSEA
  across individual cells.

- [ExperimentSubset](/packages/ExperimentSubset) Experiment objects
  such as the SummarizedExperiment or SingleCellExperiment are data
  containers for one or more matrix-like assays along with the
  associated row and column data. Often only a subset of the original
  data is needed for down-stream analysis. For example, filtering out
  poor quality samples will require excluding some columns before
  analysis. The ExperimentSubset object is a container to efficiently
  manage different subsets of the same data without having to make
  separate objects for each new subset.

- [famat](/packages/famat) Famat is made to collect data about lists
  of genes and metabolites provided by user, and to visualize it
  through a Shiny app. Information collected is: - Pathways
  containing some of the user's genes and metabolites (obtained using
  a pathway enrichment analysis). - Direct interactions between
  user's elements inside pathways. - Information about elements
  (their identifiers and descriptions). - Go terms enrichment
  analysis performed on user's genes. The Shiny app is composed of: -
  information about genes, metabolites, and direct interactions
  between them inside pathways. - an heatmap showing which elements
  from the list are in pathways (pathways are structured in
  hierarchies). - hierarchies of enriched go terms using Molecular
  Function and Biological Process.

- [FilterFFPE](/packages/FilterFFPE) This package finds and filters
  artificial chimeric reads specifically generated in next-generation
  sequencing (NGS) process of formalin-fixed paraffin-embedded (FFPE)
  tissues. These artificial chimeric reads can lead to a large number
  of false positive structural variation (SV) calls. The required
  input is an indexed BAM file of a FFPE sample.

- [flowCut](/packages/flowCut) Common techinical complications such
  as clogging can result in spurious events and fluorescence
  intensity shifting, flowCut is designed to detect and remove
  technical artifacts from your data by removing segments that show
  statistical differences from other segments.

- [fmrs](/packages/fmrs) Provides parameter estimation as well as
  variable selection in Finite Mixture of Accelerated Failure Time
  Regression and Finite Mixture of Regression Models. Furthermore,
  this package provides Ridge Regression and Elastic Net.

- [FScanR](/packages/FScanR) 'FScanR' identifies Programmed Ribosomal
  Frameshifting (PRF) events from BLASTX homolog sequence alignment
  between targeted genomic/cDNA/mRNA sequences against the peptide
  library of the same species or a close relative. The output by
  BLASTX or diamond BLASTX will be used as input of 'FScanR' and
  should be in a tabular format with 14 columns. For BLASTX, the
  output parameter should be: -outfmt '6 qseqid sseqid pident length
  mismatch gapopen qstart qend sstart send evalue bitscore qframe
  sframe'. For diamond BLASTX, the output parameter should be:
  -outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend
  sstart send evalue bitscore qframe qframe.

- [GCSFilesystem](/packages/GCSFilesystem) Mounting a Google Cloud
  bucket to a local directory. The files in the bucket can be viewed
  and read as if they are locally stored. For using the package, you
  need to install GCSDokan on Windows or gcsfuse on Linux and MacOs.

- [getDEE2](/packages/getDEE2) Digital Expression Explorer 2 (or DEE2
  for short) is a repository of processed RNA-seq data in the form of
  counts. It was designed so that researchers could undertake
  re-analysis and meta-analysis of published RNA-seq studies quickly
  and easily. As of April 2020, over 1 million SRA datasets have been
  processed. This package provides an R interface to access these
  expression data. More information about the DEE2 project can be
  found at the project homepage (http://dee2.io) and main publication
  (https://doi.org/10.1093/gigascience/giz022).

- [ggtreeExtra](/packages/ggtreeExtra) 'ggtreeExtra' extends the
  method for mapping and visualizing associated data on phylogenetic
  tree using 'ggtree'. These associated data can be mapped to
  circular layout, fan layout, or other layout tree built by 'ggtree'
  with the grammar of 'ggplot2'.

- [GSEAmining](/packages/GSEAmining) Gene Set Enrichment Analysis is
  a very powerful and interesting computational method that allows an
  easy correlation between differential expressed genes and
  biological processes. Unfortunately, although it was designed to
  help researchers to interpret gene expression data it can generate
  huge amounts of results whose biological meaning can be difficult
  to interpret. Many available tools rely on the hierarchically
  structured Gene Ontology (GO) classification to reduce reundandcy
  in the results. However, due to the popularity of GSEA many more
  gene set collections, such as those in the Molecular Signatures
  Database are emerging. Since these collections are not organized as
  those in GO, their usage for GSEA do not always give a
  straightforward answer or, in other words, getting all the
  meaninful information can be challenging with the currently
  available tools. For these reasons, GSEAmining was born to be an
  easy tool to create reproducible reports to help researchers make
  biological sense of GSEA outputs. Given the results of GSEA,
  GSEAmining clusters the different gene sets collections based on
  the presence of the same genes in the leadind edge (core) subset.
  Leading edge subsets are those genes that contribute most to the
  enrichment score of each collection of genes or gene sets. For this
  reason, gene sets that participate in similar biological processes
  should share genes in common and in turn cluster together. After
  that, GSEAmining is able to identify and represent for each
  cluster: - The most enriched terms in the names of gene sets (as
  wordclouds) - The most enriched genes in the leading edge subsets
  (as bar plots). In each case, positive and negative enrichments are
  shown in different colors so it is easy to distinguish biological
  processes or genes that may be of interest in that particular
  study.

- [GSgalgoR](/packages/GSgalgoR) A multi-objective optimization
  algorithm for disease sub-type discovery based on a non-dominated
  sorting genetic algorithm. The 'Galgo' framework combines the
  advantages of clustering algorithms for grouping heterogeneous
  'omics' data and the searching properties of genetic algorithms for
  feature selection. The algorithm search for the optimal number of
  clusters determination considering the features that maximize the
  survival difference between sub-types while keeping cluster
  consistency high.

- [GWAS.BAYES](/packages/GWAS.BAYES) This package is built to perform
  GWAS analysis for selfing species. The research related to this
  package was supported in part by National Science Foundation Award
  1853549.

- [GWENA](/packages/GWENA) The development of high-throughput
  sequencing led to increased use of co-expression analysis to go
  beyong single feature (i.e. gene) focus. We propose GWENA (Gene
  Whole co-Expression Network Analysis) , a tool designed to perform
  gene co-expression network analysis and explore the results in a
  single pipeline. It includes functional enrichment of modules of
  co-expressed genes, phenotypcal association, topological analysis
  and comparison of networks configuration between conditions.

- [HCAMatrixBrowser](/packages/HCAMatrixBrowser) The HCAMatrixBrowser
  queries the HCA matrix endpoint to download expression data and
  returns a standard Bioconductor object. It uses the LoomExperiment
  package to serve matrix data that is downloaded as HDF5 loom
  format.

- [Herper](/packages/Herper) Many tools for data analysis are not
  available in R, but are present in public repositories like conda.
  The Herper package provides a comprehensive set of functions to
  interact with the conda package managament system. With Herper
  users can install, manage and run conda packages from the comfort
  of their R session. Herper also provides an ad-hoc approach to
  handling external system requirements for R packages. For people
  developing packages with python conda dependencies we recommend
  using basilisk
  (https://bioconductor.org/packages/release/bioc/html/basilisk.html)
  to internally support these system requirments pre-hoc.

- [HPAStainR](/packages/HPAStainR) This package is built around the
  HPAStainR function. The purpose of the HPAStainR function is to
  query the visual staining data in the Human Protein Atlas to return
  a table of staining ranked cell types. The function also has
  multiple arguements to personalize to output as well to include
  cancer data, csv readable names, modify the confidence levels of
  the results and more. The other functions exist exlcusively to
  easily acquire the data required to run HPAStainR.

- [hummingbird](/packages/hummingbird) A package for detecting
  differential methylation. It exploits a Bayesian hidden Markov
  model that incorporates location dependence among genomic loci,
  unlike most existing methods that assume independence among
  observations. Bayesian priors are applied to permit information
  sharing across an entire chromosome for improved power of
  detection. The direct output of our software package is the best
  sequence of methylation states, eliminating the use of a
  subjective, and most of the time an arbitrary, threshold of p-value
  for determining significance. At last, our methodology does not
  require replication in either or both of the two comparison groups.

- [idpr](/packages/idpr) ‘idpr’ aims to integrate tools for the
  computational analysis of intrinsically disordered proteins (IDPs)
  within R. This package is used to identify known characteristics of
  IDPs for a sequence of interest with easily reported and dynamic
  results. Additionally, this package includes tools for IDP-based
  sequence analysis to be used in conjunction with other R packages.

- [ILoReg](/packages/ILoReg) ILoReg is a tool for identification of
  cell populations from scRNA-seq data. In particular, ILoReg is
  useful for finding cell populations with subtle transcriptomic
  differences. The method utilizes a self-supervised learning method,
  called Iteratitive Clustering Projection (ICP), to find cluster
  probabilities, which are used in noise reduction prior to PCA and
  the subsequent hierarchical clustering and t-SNE steps.
  Additionally, functions for differential expression analysis to
  find gene markers for the populations and gene expression
  visualization are provided.

- [infinityFlow](/packages/infinityFlow) Pipeline to analyze and
  merge data files produced by BioLegend's LEGENDScreen or BD Human
  Cell Surface Marker Screening Panel (BD Lyoplates).

- [Informeasure](/packages/Informeasure) This package compiles most
  information measures currently available: mutual information,
  conditional mutual information, interaction information, partial
  information decomposition and part mutual information. Using gene
  expression profile data, all these estimators can be employed to
  quantify nonlinear dependence between variables in biological
  regulatory network inference. The first estimator is used to infer
  bivariate network while the last four estimators are dedicated to
  analyze trivariate networks.

- [ISAnalytics](/packages/ISAnalytics) In gene therapy, stem cells
  are modified using viral vectors to deliver the therapeutic
  transgene and replace functional properties since the genetic
  modification is stable and inherited in all cell progeny. The
  retrieval and mapping of the sequences flanking the virus-host DNA
  junctions allows the identification of insertion sites (IS),
  essential for monitoring the evolution of genetically modified
  cells in vivo. A comprehensive toolkit for the analysis of IS is
  required to foster clonal trackign studies and supporting the
  assessment of safety and long term efficacy in vivo. This package
  is aimed at (1) supporting automation of IS workflow, (2)
  performing base and advance analysis for IS tracking (clonal
  abundance, clonal expansions and statistics for insertional
  mutagenesis, etc.), (3) providing basic biology insights of
  transduced stem cells in vivo.

- [lefser](/packages/lefser) lefser is an implementation in R of the
  popular "LDA Effect Size (LEfSe)" method for microbiome biomarker
  discovery. It uses the Kruskal-Wallis test, Wilcoxon-Rank Sum test,
  and Linear Discriminant Analysis to find biomarkers of groups and
  sub-groups.

- [marr](/packages/marr) marr (Maximum Rank Reproducibility) is a
  nonparametric approach that detects reproducible signals using a
  maximal rank statistic for high-dimensional biological data. In
  this R package, we implement functions that measures the
  reproducibility of features per sample pair and sample pairs per
  feature in high-dimensional biological replicate experiments. The
  user-friendly plot functions in this package also plot histograms
  of the reproducibility of features per sample pair and sample pairs
  per feature. Furthermore, our approach also allows the users to
  select optimal filtering threshold values for the identification of
  reproducible features and sample pairs based on output
  visualization checks (histograms).

- [megadepth](/packages/megadepth) This package provides an R
  interface to Megadepth by Christopher Wilks available at
  (<https://github.com/ChristopherWilks/megadepth>). It is
  particularly useful for computing the coverage of a set of genomic
  regions across bigWig or BAM files. With this package, you can
  build base-pair coverage matrices for regions or annotations of
  your choice from BigWig files. Megadepth was used to create the raw
  files provided by <https://bioconductor.org/packages/recount3>.

- [MesKit](/packages/MesKit) MesKit provides commonly used analysis
  and visualization modules based on mutational data generated by
  multi-region sequencing (MRS). This package allows to decipher ITH,
  infer metastatic routes as well as uncover the underlying process
  of mutagenesis. Shiny application was also developed for a need of
  GUI-based analysis.  As a handy tool, MesKit can facilitate the
  understanding of cancer cell evolution and its relevance to cancer
  therapeutics.

- [metabCombiner](/packages/metabCombiner) This package aligns
  LC-HRMS metabolomics datasets acquired from biologically similar
  specimens analyzed under similar, but not necessarily identical,
  conditions. Two peak-picked and aligned metabolomics feature tables
  (consisting of m/z, rt, and per-sample abundance measurements, plus
  optional identifiers & adduct annotations) are accepted as input.
  The package outputs a combined table of feature pair alignments,
  organized into groups of similar m/z, and ranked by a similarity
  score. Input tables are assumed to be acquired using similar (but
  not necessarily identical) analytical methods.

- [metabolomicsWorkbenchR](/packages/metabolomicsWorkbenchR) This
  package provides functions for interfacing with the Metabolomics
  Workbench RESTful API. Study, compound, protein and gene
  information can be searched for using the API. Methods to obtain
  study data in common Bioconductor formats such as
  SummarizedExperiment and MultiAssayExperiment are also included.

- [MethReg](/packages/MethReg) Epigenome-wide association studies
  (EWAS) detects a large number of DNA methylation differences, often
  hundreds of differentially methylated regions and thousands of
  CpGs, that are significantly associated with a disease, many are
  located in non-coding regions. Therefore, there is a critical need
  to better understand the functional impact of these CpG
  methylations and to further prioritize the significant changes.
  MethReg is an R package for integrative modeling of DNA
  methylation, target gene expression and transcription factor
  binding sites data, to systematically identify and rank functional
  CpG methylations. MethReg evaluates, prioritizes and annotates CpG
  sites with high regulatory potential using matched methylation and
  gene expression data, along with external TF-target interaction
  databases based on manually curation, ChIP-seq experiments or gene
  regulatory network analysis.

- [microbiomeExplorer](/packages/microbiomeExplorer) The
  MicrobiomeExplorer R package is designed to facilitate the analysis
  and visualization of marker-gene survey feature data. It allows a
  user to perform and visualize typical microbiome analytical
  workflows either through the command line or an interactive Shiny
  application included with the package. In addition to applying
  common analytical workflows the application enables automated
  analysis report generation.

- [MOFA2](/packages/MOFA2) The MOFA2 package contains a collection of
  tools for running and analysing MOFA models.

- [MOGAMUN](/packages/MOGAMUN) MOGAMUN is a multi-objective genetic
  algorithm that identifies active modules in a multiplex biological
  network. This allows analyzing different biological networks at the
  same time. MOGAMUN is based on NSGA-II (Non-Dominated Sorting
  Genetic Algorithm, version II), which we adapted to work on
  networks.

- [MouseFM](/packages/MouseFM) This package provides methods for
  genetic finemapping in inbred mice by taking advantage of their
  very high homozygosity rate (>95%).

- [MSEADbi](/packages/MSEADbi) Interface to construct annotation
  package for MSEA (MSEA.XXX.pb.db). The program design is same as
  Bioconductor LRBaseDbi or MeSHDbi pacakge, and the usage is also
  the same as these packages.

- [msImpute](/packages/msImpute) MsImpute is a package for imputation
  of peptide intensity in proteomics experiments. It additionally
  contains tools for MAR/MNAR diagnosis and assessment of distortions
  to the probability distribution of the data post imputation.
  Currently, msImpute completes missing values by low-rank
  approximation of the underlying data matrix.

- [MSPrep](/packages/MSPrep) Package performs summarization of
  replicates, filtering by frequency, several different options for
  imputing missing data, and a variety of options for transforming,
  batch correcting, and normalizing data.

- [MSstatsConvert](/packages/MSstatsConvert) MSstatsConvert provides
  tools for importing reports of Mass Spectrometry data processing
  tools into R format suitable for statistical analysis using the
  MSstats and MSstatsTMT packages.

- [MSstatsPTM](/packages/MSstatsPTM) MSstatsPTM provides general
  statistical methods for quantitative characterization of
  post-translational modifications (PTMs). Typically, the analysis
  involves the quantification of PTM sites (i.e., modified residues)
  and their corresponding proteins, as well as the integration of the
  quantification results. MSstatsPTM provides functions for
  summarization, estimation of PTM site abundance, and detection of
  changes in PTMs across experimental conditions.

- [MSstatsTMTPTM](/packages/MSstatsTMTPTM) Tools for Post
  Translational Modification (PTM) and protein significance analysis
  in shotgun mass spectrometry-based proteomic experiments with
  tandem mass tag (TMT) labeling. The functions in this package
  should be used after PTM/protein summarization. They can be used to
  both plot the summarized results and model the summarized datasets.

- [MultiBaC](/packages/MultiBaC) MultiBaC is a strategy to correct
  batch effects from multiomic datasets distributed across different
  labs or data acquisition events. MultiBaC is the first Batch effect
  correction algorithm that dealing with batch effect correction in
  multiomics datasets. MultiBaC is able to remove batch effects
  across different omics generated within separate batches provided
  that at least one common omic data type is included in all the
  batches considered.

- [multicrispr](/packages/multicrispr) This package is for designing
  Crispr/Cas9 and Prime Editing experiments. It contains functions to
  (1) define and transform genomic targets, (2) find spacers (4)
  count offtarget (mis)matches, and (5) compute Doench2016/2014
  targeting efficiency. Care has been taken for multicrispr to scale
  well towards large target sets, enabling the design of large
  Crispr/Cas9 libraries.

- [multiGSEA](/packages/multiGSEA) Extracted features from pathways
  derived from 8 different databases (KEGG, Reactome, Biocarta, etc.)
  can be used on transcriptomic, proteomic, and/or metabolomic level
  to calculate a combined GSEA-based enrichment score.

- [musicatk](/packages/musicatk) Mutational signatures are
  carcinogenic exposures or aberrant cellular processes that can
  cause alterations to the genome. We created musicatk (MUtational
  SIgnature Comprehensive Analysis ToolKit) to address shortcomings
  in versatility and ease of use in other pre-existing computational
  tools. Although many different types of mutational data have been
  generated, current software packages do not have a flexible
  framework to allow users to mix and match different types of
  mutations in the mutational signature inference process. Musicatk
  enables users to count and combine multiple mutation types,
  including SBS, DBS, and indels. Musicatk calculates replication
  strand, transcription strand and combinations of these features
  along with discovery from unique and proprietary genomic feature
  associated with any mutation type. Musicatk also implements several
  methods for discovery of new signatures as well as methods to infer
  exposure given an existing set of signatures. Musicatk provides
  functions for visualization and downstream exploratory analysis
  including the ability to compare signatures between cohorts and
  find matching signatures in COSMIC V2 or COSMIC V3.

- [NanoMethViz](/packages/NanoMethViz) NanoMethViz is a toolkit for
  visualising methylation data from Oxford Nanopore sequencing. It
  can be used to explore methylation patterns from reads derived from
  Oxford Nanopore direct DNA sequencing with methylation called by
  callers including nanopolish, f5c and megalodon. The plots in this
  package allow the visualisation of methylation profiles aggregated
  over experimental groups and across classes of genomic features.

- [ncRNAtools](/packages/ncRNAtools) ncRNAtools provides a set of
  basic tools for handling and analyzing non-coding RNAs. These
  include tools to access the RNAcentral database and to predict and
  visualize the secondary structure of non-coding RNAs. The package
  also provides tools to read, write and interconvert the file
  formats most commonly used for representing such secondary
  structures.

- [nearBynding](/packages/nearBynding) Provides a pipeline to discern
  RNA structure at and proximal to the site of protein binding within
  regions of the transcriptome defined by the user. CLIP
  protein-binding data can be input as either aligned BAM or
  peak-called bedGraph files. RNA structure can either be predicted
  internally from sequence or users have the option to input their
  own RNA structure data. RNA structure binding profiles can be
  visually and quantitatively compared across multiple formats.

- [Nebulosa](/packages/Nebulosa) This package provides a enhanced
  visualization of single-cell data based on gene-weighted density
  estimation. Nebulosa recovers the signal from dropped-out features
  and allows the inspection of the joint expression from multiple
  features (e.g. genes). Seurat and SingleCellExperiment objects can
  be used within Nebulosa.

- [NewWave](/packages/NewWave) A model designed for dimensionality
  reduction and batch effect removal for scRNA-seq data. It is
  designed to be massively parallelizable using shared objects that
  prevent memory duplication, and it can be used with different
  mini-batch approaches in order to reduce time consumption. It
  assumes a negative binomial distribution for the data with a
  dispersion parameter that can be both commonwise across gene both
  genewise.

- [Omixer](/packages/Omixer) Omixer - an R package for multivariate
  and reproducible randomization with lab-friendly sample layouts.
  Omixer ensures optimal sample distribution across batches with
  well-documented methods, and can output intuitive sample sheets for
  the wet lab if needed.

- [padma](/packages/padma) Use multiple factor analysis to calculate
  individualized pathway-centric scores of deviation with respect to
  the sampled population based on multi-omic assays (e.g., RNA-seq,
  copy number alterations, methylation, etc). Graphical and numerical
  outputs are provided to identify highly aberrant individuals for a
  particular pathway of interest, as well as the gene and omics
  drivers of aberrant multi-omic profiles.

- [pageRank](/packages/pageRank) Implemented temporal PageRank
  analysis as defined by Rozenshtein and Gionis. Implemented
  multiplex PageRank as defined by Halu et al. Applied temporal and
  multiplex PageRank in gene regulatory network analysis.

- [PeacoQC](/packages/PeacoQC) This is a package that includes
  pre-processing and quality control functions that can remove margin
  events, compensate and transform the data and that will use
  PeacoQCSignalStability for quality control. This last function will
  first detect peaks in each channel of the flowframe. It will remove
  anomalies based on the IsolationTree function and the MAD outlier
  detection method. This package can be used for both flow- and mass
  cytometry data.

- [periodicDNA](/packages/periodicDNA) This R package helps the user
  identify k-mers (e.g. di- or tri-nucleotides) present periodically
  in a set of genomic loci (typically regulatory elements). The
  functions of this package provide a straightforward approach to
  find periodic occurrences of k-mers in DNA sequences, such as
  regulatory elements. It is not aimed at identifying motifs
  separated by a conserved distance; for this type of analysis,
  please visit MEME website.

- [PhosR](/packages/PhosR) PhosR is a package for the comprenhensive
  analysis of phosphoproteomic data. There are two major components
  to PhosR: processing and downstream analysis. PhosR consists of
  various processing tools for phosphoproteomics data including
  filtering, imputation, normalisation, and functional analysis for
  inferring active kinases and signalling pathways.

- [pipeComp](/packages/pipeComp) A simple framework to facilitate the
  comparison of pipelines involving various steps and parameters. The
  `pipelineDefinition` class represents pipelines as, minimally, a
  set of functions consecutively executed on the output of the
  previous one, and optionally accompanied by step-wise evaluation
  and aggregation functions. Given such an object, a set of
  alternative parameters/methods, and benchmark datasets, the
  `runPipeline` function then proceeds through all combinations
  arguments, avoiding recomputing the same step twice and compiling
  evaluations on the fly to avoid storing potentially large
  intermediate data.

- [POMA](/packages/POMA) POMA introduces a structured, reproducible
  and easy-to-use workflow for the visualization, pre-processing,
  exploratory and statistical analysis of mass spectrometry data. The
  main aim of POMA is to enable a flexible data cleaning and
  statistical analysis processes in one comprehensible and
  user-friendly R package. This package also has a Shiny app version
  that implements all POMA functions. See
  https://github.com/pcastellanoescuder/POMAShiny.

- [preciseTAD](/packages/preciseTAD) preciseTAD provides functions to
  predict the location of boundaries of topologically associated
  domains (TADs) and chromatin loops at base-level resolution. As an
  input, it takes BED-formatted genomic coordinates of domain
  boundaries detected from low-resolution Hi-C data, and coordinates
  of high-resolution genomic annotations from ENCODE or other
  consortia. preciseTAD employs several feature engineering
  strategies and resampling techniques to address class imbalance,
  and trains an optimized random forest model for predicting
  low-resolution domain boundaries. Translated on a base-level,
  preciseTAD predicts the probability for each base to be a boundary.
  Density-based clustering and scalable partitioning techniques are
  used to detect precise boundary regions and summit points. Compared
  with low-resolution boundaries, preciseTAD boundaries are highly
  enriched for CTCF, RAD21, SMC3, and ZNF143 signal and more
  conserved across cell lines. The pre-trained model can accurately
  predict boundaries in another cell line using CTCF, RAD21, SMC3,
  and ZNF143 annotation data for this cell line.

- [proActiv](/packages/proActiv) Most human genes have multiple
  promoters that control the expression of different isoforms. The
  use of these alternative promoters enables the regulation of
  isoform expression pre-transcriptionally. Alternative promoters
  have been found to be important in a wide number of cell types and
  diseases. proActiv is an R package that enables the analysis of
  promoters from RNA-seq data. proActiv uses aligned reads as input,
  and generates counts and normalized promoter activity estimates for
  each annotated promoter. In particular, proActiv accepts junction
  files from TopHat2 or STAR or BAM files as inputs. These estimates
  can then be used to identify which promoter is active, which
  promoter is inactive, and which promoters change their activity
  across conditions.

- [QFeatures](/packages/QFeatures) The QFeatures infrastructure
  enables the management and processing of quantitative features for
  high-throughput mass spectrometry assays. It provides a familiar
  Bioconductor user experience to manages quantitative data across
  different assay levels (such as peptide spectrum matches, peptides
  and proteins) in a coherent and tractable format.

- [RadioGx](/packages/RadioGx) Computational tool box for
  radio-genomic analysis which integrates radio-response data,
  radio-biological modelling and comprehensive cell line annotations
  for hundreds of cancer cell lines. The 'RadioSet' class enables
  creation and manipulation of standardized datasets including
  information about cancer cells lines, radio-response assays and
  dose-response indicators. Included methods allow fitting and
  plotting dose-response data using established radio-biological
  models along with quality control to validate results. Additional
  functions related to fitting and plotting dose response curves,
  quantifying statistical correlation and calculating area under the
  curve (AUC) or survival fraction (SF) are included. For more
  details please see the included documentation, references, as well
  as: Manem, V. et al (2018) <doi:10.1101/449793>.

- [rebook](/packages/rebook) Provides utilities to re-use content
  across chapters of a Bioconductor book. This is mostly based on
  functionality developed while writing the OSCA book, but
  generalized for potential use in other large books with heavy
  compute. Also contains some functions to assist book deployment.

- [recount3](/packages/recount3) Explore and download data from the
  recount project available at
  https://jhubiostatistics.shinyapps.io/recount3/. Using the recount3
  package you can download RangedSummarizedExperiment objects at the
  gene, exon or exon-exon junctions level, the raw counts, the
  phenotype metadata used, the urls to the sample coverage bigWig
  files. The RangedSummarizedExperiment objects can be used by
  different packages for performing differential expression analysis.
  Using data from the recount3 project you can perform
  annotation-agnostic differential expression analyses as described
  at http://doi.org/TODO.

- [recountmethylation](/packages/recountmethylation) Access
  cross-study compilations of DNA methylation array databases.
  Database files can be downloaded and accessed using provided
  functions. Background about database file types (HDF5 and
  HDF5-SummarizedExperiment), SummarizedExperiment classes, and
  examples for data handling, validation, and analyses, can be found
  in the package vignettes. Note the disclaimer on package load, and
  consult the main manuscript for further info.

- [RegEnrich](/packages/RegEnrich) This package is a pipeline to
  identify the key gene regulators in a biological process, for
  example in cell differentiation and in cell development after
  stimulation. There are four major steps in this pipeline: (1)
  differential expression analysis; (2) regulator-target network
  inference; (3) enrichment analysis; and (4) regulators scoring and
  ranking.

- [ResidualMatrix](/packages/ResidualMatrix) Provides delayed
  computation of a matrix of residuals after fitting a linear model
  to each column of an input matrix. Also supports partial
  computation of residuals where selected factors are to be preserved
  in the output matrix. Implements a number of efficient methods for
  operating on the delayed matrix of residuals, most notably matrix
  multiplication and calculation of row/column sums or means.

- [Rfastp](/packages/Rfastp) Rfastp is an R wrapper of fastp
  developed in c++. fastp performs quality control for fastq files.
  including low quality bases trimming, polyX trimming, adapter
  auto-detection and trimming, paired-end reads merging, UMI
  sequence/id handling. Rfastp can concatenate multiple files into
  one file (like shell command cat) and accept multiple files as
  input.

- [RIPAT](/packages/RIPAT) RIPAT is developed as an R package for
  retroviral integration sites annotation and distribution analysis.
  RIPAT needs local alignment results from BLAST and BLAT. Specific
  input format is depicted in RIPAT manual. RIPAT provides RV
  integration pattern analysis result as forms of R objects, excel
  file with multiple sheets and plots.

- [rnaEditr](/packages/rnaEditr) RNAeditr analyzes site-specific RNA
  editing events, as well as hyper-editing regions. The editing
  frequencies can be tested against binary, continuous or survival
  outcomes. Multiple covariate variables as well as interaction
  effects can also be incorporated in the statistical models.

- [RnaSeqSampleSize](/packages/RnaSeqSampleSize) RnaSeqSampleSize
  package provides a sample size calculation method based on negative
  binomial model and the exact test for assessing differential
  expression analysis of RNA-seq data. It controls FDR for multiple
  testing and utilizes the average read count and dispersion
  distributions from real data to estimate a more reliable sample
  size. It is also equipped with several unique features, including
  estimation for interested genes or pathway, power curve
  visualization, and parameter optimization.

- [rsemmed](/packages/rsemmed) A programmatic interface to the
  Semantic MEDLINE database. It provides functions for searching the
  database for concepts and finding paths between concepts. Path
  searching can also be tailored to user specifications, such as
  placing restrictions on concept types and the type of link between
  concepts. It also provides functions for summarizing and
  visualizing those paths.

- [Rtpca](/packages/Rtpca) R package for performing thermal proximity
  co-aggregation analysis with thermal proteome profiling datasets to
  analyse protein complex assembly and (differential) protein-protein
  interactions across conditions.

- [sangeranalyseR](/packages/sangeranalyseR) This package builds on
  sangerseqR to allow users to create contigs from collections of
  Sanger sequencing reads. It provides a wide range of options for a
  number of commonly-performed actions including read trimming,
  detecting secondary peaks, and detecting indels using a reference
  sequence. All parameters can be adjusted interactively either in R
  or in the associated Shiny applications. There is extensive online
  documentation, and the package can outputs detailed HTML reports,
  including chromatograms.

- [SCATE](/packages/SCATE) SCATE is a software tool for extracting
  and enhancing the sparse and discrete Single-cell ATAC-seq Signal.
  Single-cell sequencing assay for transposase-accessible chromatin
  (scATAC-seq) is the state-of-the-art technology for analyzing
  genome-wide regulatory landscapes in single cells. Single-cell
  ATAC-seq data are sparse and noisy, and analyzing such data is
  challenging. Existing computational methods cannot accurately
  reconstruct activities of individual cis-regulatory elements (CREs)
  in individual cells or rare cell subpopulations. SCATE was
  developed to adaptively integrate information from co-activated
  CREs, similar cells, and publicly available regulome data and
  substantially increase the accuracy for estimating activities of
  individual CREs. We demonstrate that SCATE can be used to better
  reconstruct the regulatory landscape of a heterogeneous sample.

- [scCB2](/packages/scCB2) scCB2 is an R package implementing CB2 for
  distinguishing real cells from empty droplets in droplet-based
  single cell RNA-seq experiments (especially for 10x Chromium). It
  is based on clustering similar barcodes and calculating Monte-Carlo
  p-value for each cluster to test against background distribution.
  This cluster-level test outperforms single-barcode-level tests in
  dealing with low count barcodes and homogeneous sequencing library,
  while keeping FDR well controlled.

- [scDataviz](/packages/scDataviz) In the single cell World, which
  includes flow cytometry, mass cytometry, single-cell RNA-seq
  (scRNA-seq), and others, there is a need to improve data
  visualisation and to bring analysis capabilities to researchers
  even from non-technical backgrounds. scDataviz attempts to fit into
  this space, while also catering for advanced users. Additonally,
  due to the way that scDataviz is designed, which is based on
  SingleCellExperiment, it has a 'plug and play' feel, and
  immediately lends itself as flexibile and compatibile with studies
  that go beyond scDataviz. Finally, the graphics in scDataviz are
  generated via the ggplot engine, which means that users can 'add
  on' features to these with ease.

- [SCFA](/packages/SCFA) Subtyping via Consensus Factor Analysis
  (SCFA) can efficiently remove noisy signals from consistent
  molecular patterns in multi-omics data. SCFA first uses an
  autoencoder to select only important features and then repeatedly
  performs factor analysis to represent the data with different
  numbers of factors. Using these representations, it can reliably
  identify cancer subtypes and accurately predict risk scores of
  patients.

- [scp](/packages/scp) Utility functions for manipulating,
  processing, and analyzing mass spectrometry-based single-cell
  proteomics (SCP) data. The package is an extension to the
  'QFeatures' package designed for SCP applications.

- [scRepertoire](/packages/scRepertoire) scRepertoire was built to
  process data derived from the 10x Genomics Chromium Immune
  Profiling for both T-cell receptor (TCR) and immunoglobulin (Ig)
  enrichment workflows and subsequently interacts with the popular
  Seurat and SingleCellExperiment R packages.

- [scuttle](/packages/scuttle) Provides basic utility functions for
  performing single-cell analyses, focusing on simple normalization,
  quality control and data transformations. Also provides some helper
  functions to assist development of other packages.

- [SeqGate](/packages/SeqGate) Filtering of lowly expressed features
  (e.g. genes) is a common step before performing statistical
  analysis, but an arbitrary threshold is generally chosen. SeqGate
  implements a method that rationalize this step by the analysis of
  the distibution of counts in replicate samples. The gate is the
  threshold above which sequenced features can be considered as
  confidently quantified.

- [simplifyEnrichment](/packages/simplifyEnrichment) A new method
  (binary cut) is proposed to effectively cluster GO terms into
  groups from the semantic similarity matrix. Summaries of GO terms
  in each cluster are visualized by word clouds.

- [snifter](/packages/snifter) Provides an R wrapper for the
  implementation of FI-tSNE from the python package openTNSE. See
  Poličar et al. (2019) <doi:10.1101/731877> and the algorithm
  described by Linderman et al. (2018)
  <doi:10.1038/s41592-018-0308-4>.

- [SpatialDecon](/packages/SpatialDecon) Using spatial or bulk gene
  expression data, estimates abundance of mixed cell types within
  each observation. Based on "Advances in mixed cell deconvolution
  enable quantification of cell types in spatially-resolved gene
  expression data", Danaher (2020). Designed for use with the
  NanoString GeoMx platform, but applicable to any gene expression
  data.

- [SpatialExperiment](/packages/SpatialExperiment) Defines S4 classes
  for storing data for spatial experiments. Main examples are
  reported by using seqFISH and 10x-Visium Spatial Gene Expression
  data. This includes specialized methods for storing, retrieving
  spatial coordinates, 10x dedicated parameters and their handling.

- [spatialHeatmap](/packages/spatialHeatmap) The spatialHeatmap
  package provides functionalities for visualizing cell-, tissue- and
  organ-specific data of biological assays by coloring the
  corresponding spatial features defined in anatomical images
  according to a numeric color key.

- [Spectra](/packages/Spectra) The Spectra package defines an
  efficient infrastructure for storing and handling mass spectrometry
  spectra and functionality to subset, process, visualize and compare
  spectra data. It provides different implementations (backends) to
  store mass spectrometry data. These comprise backends tuned for
  fast data access and processing and backends for very large data
  sets ensuring a small memory footprint.

- [SPsimSeq](/packages/SPsimSeq) SPsimSeq uses a specially designed
  exponential family for density estimation to constructs the
  distribution of gene expression levels from a given real RNA
  sequencing data (single-cell or bulk), and subsequently simulates a
  new dataset from the estimated marginal distributions using
  Gaussian-copulas to retain the dependence between genes. It allows
  simulation of multiple groups and batches with any required sample
  size and library size.

- [systemPipeShiny](/packages/systemPipeShiny) systemPipeShiny (SPS)
  extends the widely used systemPipeR (SPR) workflow environment with
  a versatile graphical user interface provided by a Shiny App. This
  allows non-R users, such as experimentalists, to run many
  systemPipeR’s workflow designs, control, and visualization
  functionalities interactively without requiring knowledge of R.
  Most importantly, SPS has been designed as a general purpose
  framework for interacting with other R packages in an intuitive
  manner. Like most Shiny Apps, SPS can be used on both local
  computers as well as centralized server-based deployments that can
  be accessed remotely as a public web service for using SPR’s
  functionalities with community and/or private data. The framework
  can integrate many core packages from the R/Bioconductor ecosystem.
  Examples of SPS’ current functionalities include: (a) interactive
  creation of experimental designs and metadata using an easy to use
  tabular editor or file uploader; (b) visualization of workflow
  topologies combined with auto-generation of R Markdown preview for
  interactively designed workflows; (d) access to a wide range of
  data processing routines; (e) and an extendable set of
  visualization functionalities. Complex visual results can be
  managed on a 'Canvas Workbench’ allowing users to organize and to
  compare plots in an efficient manner combined with a session
  snapshot feature to continue work at a later time. The present
  suite of pre-configured visualization examples. The modular design
  of SPR makes it easy to design custom functions without any
  knowledge of Shiny, as well as extending the environment in the
  future with contributions from the community.

- [TADCompare](/packages/TADCompare) TADCompare is an R package
  designed to identify and characterize differential Topologically
  Associated Domains (TADs) between multiple Hi-C contact matrices.
  It contains functions for finding differential TADs between two
  datasets, finding differential TADs over time and identifying
  consensus TADs across multiple matrices. It takes all of the main
  types of HiC input and returns simple, comprehensive, easy to
  analyze results.

- [tidySingleCellExperiment](/packages/tidySingleCellExperiment)
  tidySingleCellExperiment is an adapter that abstracts the
  'SingleCellExperiment' container in the form of tibble and allows
  the data manipulation, plotting and nesting using 'tidyverse'

- [tidySummarizedExperiment](/packages/tidySummarizedExperiment)
  tidySummarizedExperiment is an adapter that abstracts the
  'SingleCellExperiment' container in the form of tibble and allows
  the data manipulation, plotting and nesting using 'tidyverse'

- [TileDBArray](/packages/TileDBArray) Implements a DelayedArray
  backend for reading and writing dense or sparse arrays in the
  TileDB format. The resulting TileDBArrays are compatible with all
  Bioconductor pipelines that can accept DelayedArray instances.

- [TimiRGeN](/packages/TimiRGeN) TimiRGeN (Time Incorporated miR-mRNA
  Generation of Networks) is a novel R package which functionally
  analyses and integrates time course miRNA-mRNA differential
  expression data. This tool can generate small networks within R or
  export results into cytoscape or pathvisio for more detailed
  network construction and hypothesis generation. This tool is
  created for researchers that wish to dive deep into time series
  multi-omic datasets. TimiRGeN goes further than many other tools in
  terms of data reduction. Here, potentially hundreds of thousands of
  potential miRNA-mRNA interactions can be whittled down into a
  handful of high confidence miRNA-mRNA interactions effecting a
  signalling pathway, across a time course.

- [tomoda](/packages/tomoda) This package provides many easy-to-use
  methods to analyze and visualize tomo-seq data. The tomo-seq
  technique is based on cryosectioning of tissue and performing
  RNA-seq on consecutive sections. (Reference: Kruse F, Junker JP,
  van Oudenaarden A, Bakkers J. Tomo-seq: A method to obtain
  genome-wide expression data with spatial resolution. Methods Cell
  Biol. 2016;135:299-307. doi:10.1016/bs.mcb.2016.01.006) The main
  purpose of the package is to find zones with similar
  transcriptional profiles and spatially expressed genes in a
  tomo-seq sample. Several visulization functions are available to
  create easy-to-modify plots.

- [ToxicoGx](/packages/ToxicoGx) Contains a set of functions to
  perform large-scale analysis of toxicogenomic data, providing a
  standardized data structure to hold information relevant to
  annotation, visualization and statistical analysis of toxicogenomic
  data.

- [transomics2cytoscape](/packages/transomics2cytoscape)
  transomics2cytoscape generates a file for 3D transomics
  visualization by providing input that specifies the IDs of multiple
  KEGG pathway layers, their corresponding Z-axis heights, and an
  input that represents the edges between the pathway layers. The
  edges are used, for example, to describe the relationships between
  kinase on a pathway and enzyme on another pathway. This package
  automates creation of a transomics network as shown in the figure
  in Yugi.2014 (https://doi.org/10.1016/j.celrep.2014.07.021) using
  Cytoscape automation (https://doi.org/10.1186/s13059-019-1758-4).

- [UMI4Cats](/packages/UMI4Cats) UMI-4C is a technique that allows
  characterization of 3D chromatin interactions with a bait of
  interest, taking advantage of a sonication step to produce unique
  molecular identifiers (UMIs) that help remove duplication bias,
  thus allowing a better differential comparsion of chromatin
  interactions between conditions. This package allows processing of
  UMI-4C data, starting from FastQ files provided by the sequencing
  facility. It provides two statistical methods for detecting
  differential contacts and includes a visualization function to plot
  integrated information from a UMI-4C assay.

- [uncoverappLib](/packages/uncoverappLib) a Shiny application
  containing a suite of graphical and statistical tools to support
  clinical assessment of low coverage regions.It displays three web
  pages each providing a different analysis module: Coverage
  analysis, calculate AF by allele frequency app and binomial
  distribution.

- [velociraptor](/packages/velociraptor) This package provides
  Bioconductor-friendly wrappers for RNA velocity calculations in
  single-cell RNA-seq data. We use the basilisk package to manage
  Conda environments, and the zellkonverter package to convert data
  structures between SingleCellExperiment (R) and AnnData (Python).
  The information produced by the velocity methods is stored in the
  various components of the SingleCellExperiment class.

- [VERSO](/packages/VERSO) Mutations that rapidly accumulate in viral
  genomes during a pandemic can be used to track the evolution of the
  virus and, accordingly, unravel the viral infection network. To
  this extent, sequencing samples of the virus can be employed to
  estimate models from genomic epidemiology and may serve, for
  instance, to estimate the proportion of undetected infected people
  by uncovering cryptic transmissions, as well as to predict likely
  trends in the number of infected, hospitalized, dead and recovered
  people. VERSO is an algorithmic framework that processes variants
  profiles from viral samples to produce phylogenetic models of viral
  evolution. The approach solves a Boolean Matrix Factorization
  problem with phylogenetic constraints, by maximizing a
  log-likelihood function. VERSO includes two separate and subsequent
  steps; in this package we provide an R implementation of VERSO STEP
  1.

- [VplotR](/packages/VplotR) The pattern of digestion and protection
  from DNA nucleases such as DNAse I, micrococcal nuclease, and Tn5
  transposase can be used to infer the location of associated
  proteins. This package contains useful functions to analyze
  patterns of paired-end sequencing fragment density. VplotR
  facilitates the generation of V-plots and footprint profiles over
  single or aggregated genomic loci of interest.

- [wpm](/packages/wpm) This is a shiny application for creating
  well-plate plans. It uses a backtracking-inspired algorithm to
  place samples on plates based on specific neighborhood constraints.

- [zellkonverter](/packages/zellkonverter) Provides methods to
  convert between Python AnnData objects and SingleCellExperiment
  objects. These are primarily intended for use by downstream
  Bioconductor packages that wrap Python methods for single-cell data
  analysis. It also includes functions to read and write H5AD files
  used for saving AnnData objects to disk.

New Data Experiment Packages
=====================

There are 9 new data experiment packages in this release of Bioconductor.

- [celldex](/packages/celldex) Provides a collection of reference
  expression datasets with curated cell type labels, for use in
  procedures like automated annotation of single-cell data or
  deconvolution of bulk RNA-seq.

- [clustifyrdatahub](/packages/clustifyrdatahub) References made from
  external single-cell mRNA sequencing data sets, stored as average
  gene expression matrices. For use with clustifyr
  <https://bioconductor.org/packages/clustifyr> to assign cell type
  identities.

- [DropletTestFiles](/packages/DropletTestFiles) Assorted files
  generated from droplet-based single-cell protocols, to be used for
  testing functions in DropletUtils. Primarily intended for storing
  files that directly come out of processing pipelines like 10X
  Genomics' CellRanger software, prior to the formation of a
  SingleCellExperiment object. Unlike other packages, this is not
  designed to provide objects that are immediately ready for
  analysis.

- [FieldEffectCrc](/packages/FieldEffectCrc) Processed RNA-seq data
  for 1,139 human primary colorectal tissue samples across three
  phenotypes, including tumor, normal adjacent-to-tumor, and healthy,
  available as Synapse ID syn22237139 on synapse.org. Data have been
  parsed into SummarizedExperiment objects available via
  ExperimentHub to facilitate reproducibility and extension of
  results from Dampier et al. (PMCID: PMC7386360, PMID: 32764205).

- [MethylSeqData](/packages/MethylSeqData) Base-level (i.e.
  cytosine-level) counts for a collection of public bisulfite-seq
  datasets (e.g., WGBS and RRBS), provided as SummarizedExperiment
  objects with sample- and base-level metadata.

- [NanoporeRNASeq](/packages/NanoporeRNASeq) The NanoporeRNASeq
  package contains long read RNA-Seq data generated using Oxford
  Nanopore Sequencing. The data consists of 6 samples from two human
  cell lines (K562 and MCF7) that were generated by the SG-NEx
  project. Each of these cell lines has three replicates, with 1
  direct RNA sequencing data and 2 cDNA sequencing data. Reads are
  aligned to chromosome 22 (Grch38) and stored as bam files. The
  original data is from the SG-NEx project.

- [SCATEData](/packages/SCATEData) SCATEData is an ExperimentHub
  package for SCATE which is a software tool for extracting and
  enhancing the sparse and discrete Single-cell ATAC-seq Signal.

- [timecoursedata](/packages/timecoursedata) This data package
  contains timecourse gene expression data sets. The first dataset,
  from Shoemaker et al, consists of microarray samples from lung
  tissue of mice exposed to different influenzy strains from 14
  timepoints. The two other datasets are leaf and root samples from
  sorghum crops exposed to pre- and post-flowering drought stress and
  a control condition, sampled across the plants lifetime.

- [TMExplorer](/packages/TMExplorer) This package provides a tool to
  search and download a collection of tumour microenvironment
  single-cell RNA sequencing datasets and their metadata. TMExplorer
  aims to act as a single point of entry for users looking to study
  the tumour microenvironment at the single cell level. Users can
  quickly search available datasets using the metadata table and then
  download the ones they are interested in for analysis.

New Annotation Packages
=====================

There are two new annotation packages in this release of Bioconductor.

- [metaboliteIDmapping](/packages/metaboliteIDmapping) The package provides a
    comprehensive mapping table of nine different Metabolite ID formats and
    their common name. The data has been collected and merged from four publicly
    available source, including HMDB, Comptox Dashboard, ChEBI, and the graphite
    Bioconductor R package.
    
- [geneplast.data](/packages/geneplast.data) The package geneplast.data provides
    an interface for obtaining input data used in the analyses pipelines from
    geneplast package. Objects containing species, phylogenetic trees, and
    orthology information of eukaryotes from different orthologous databases are
    provided

New Workflow Packages
=====================

There is 1 new workflow package in this release of Bioconductor. 

- [BP4RNAseq](/packages/BP4RNAseq) An automated pipe for reproducible
  RNA-seq analysis with the minimal efforts from researchers. The
  package can process bulk RNA-seq data and single-cell RNA-seq data.
  You can only provide the taxa name and the accession id of RNA-seq
  data deposited in the National Center for Biotechnology Information
  (NCBI). After a cup of tea or longer, you will get formated gene
  expression data as gene count and transcript count based on both
  alignment-based and alignment-free workflows.


NEWS from new and existing Software Packages
===================================


[a4](https://bioconductor.org/packages/a4)
--

                       Changes in version 1.37.2                        

- add imports utils packageDescription

[a4Base](https://bioconductor.org/packages/a4Base)
------

                       Changes in version 1.37.3                        

- fix specification of link for fct with different name help page

                       Changes in version 1.37.2                        

- use Imports rather than Depends + use roxygen2 for documentation

- use Authors@R

- replace toptable (deprecated in limma 3.36.0) by limma:::.topTableT

                       Changes in version 1.37.1                        

- spectralMap: fix legend win-32

[a4Classif](https://bioconductor.org/packages/a4Classif)
---------

                       Changes in version 1.37.1                        

- use Imports rather than Depends + use roxygen2 for documentation

- use Authors@R

- add vignette, examples

[a4Core](https://bioconductor.org/packages/a4Core)
------

                       Changes in version 1.37.1                        

- use Imports rather than Depends + use roxygen2 for documentation

- use Authors@R

- add vignette

[a4Preproc](https://bioconductor.org/packages/a4Preproc)
---------

                       Changes in version 1.37.1                        

- use Imports rather than Depends + use roxygen2 for documentation

- use Authors@R

- add vignette

[a4Reporting](https://bioconductor.org/packages/a4Reporting)
-----------

                       Changes in version 1.37.1                        

- use Imports rather than Depends + use roxygen2 for documentation

- use Authors@R

- add vignette

[ADImpute](https://bioconductor.org/packages/ADImpute)
--------

                        Changes in version 2.1.0                        

- Added a NEWS.md file to track changes to the package.

[AffiXcan](https://bioconductor.org/packages/AffiXcan)
--------

                        Changes in version 1.7.1                        

BUG FIXED

- Null eigenvectors for TBA with null variance or NA values in the
  training
  set are returned correctly

[aggregateBioVar](https://bioconductor.org/packages/aggregateBioVar)
---------------

                 Changes in version 0.99.0 (2020-09-03)                 

- Submitted to Bioconductor

[alevinQC](https://bioconductor.org/packages/alevinQC)
--------

                        Changes in version 1.5.2                        

- Added selected summary distributions to the reports, and a function
to plot histograms of arbitrary numeric columns of the cbTable

[AlphaBeta](https://bioconductor.org/packages/AlphaBeta)
---------

                 Changes in version 1.2.3 (2020-06-06)                  

- Fixed bugs.

- Compatible with new version of dplyr,data.tables

                 Changes in version 1.2.1 (2020-04-20)                  

- Fixed bugs.

- Made the following significant changes
  + Using ape package to make tree.
  + Added some plots.

- Added SOMA functions.

- Plotting estimates

[AlpsNMR](https://bioconductor.org/packages/AlpsNMR)
-------

                 Changes in version 2.99.5 (2020-10-14)                 

- Bug in bp_kfold_VIP_analysis solved
- Several packages moved from import to depends
- Reexport of some functions removed
- to_rDolphin_blood code reorganized
- Typos removed from tutorial
- norm_pqn_diagnostic$norm_factor used in tutorial instead of plot it

                 Changes in version 2.99.4 (2020-09-28)                 

- Warning in plot_interactive function added
- Suppressed other warnings of plot_interactive function

                 Changes in version 2.99.3 (2020-09-21)                 

- sapply calls changed for vapply
- Bioconductor installation instructions included
- MIT license removed
- LazyData: TRUE removed
- Excessive print statements removed from vignettes
- sessionInfo() added to end of vignettes
- Created inst/script directoy to describe inst/extdata source and
creation #TODO falta rellenar el archivo
- Commented out code removed

                 Changes in version 2.99.2 (2020-08-26)                 

- AlpsNMR.Rproj removed from git repository
- Reduced demo dataset to avoid package size > 5 MB
- Modified introduction to alpsnmr vignette and some tests to work
with reduced demo dataset

                 Changes in version 2.99.1 (2020-08-25)                 

- AlpsNMR.Rproj added to gitignore
- Modified examples to avoid create files in main package folder

                 Changes in version 2.99.0 (2020-08-24)                 

- Added bootstrap and permutation method and some plots related to it
- Minor modifications for bioconductor submision

                Changes in version 2.5.9002 (2020-05-25)                

- Changes to pass BiocCheck
- Added permutation test and permutation test plot to
nmr_data_analysis

                Changes in version 2.4.9002 (2020-05-13)                

- Changes to pass checks for R4

                     Changes in version 2.3.3.9002                      

- NIHS_specific removed
- Tests coverage up to 30%
- Update of save_profiling_plots
- Add tutorial
- Remotes installation
- nmr_diagnose is deprecated. Since nmr_diagnose was only used for
getting extra normalization information, it was been replaced with
nmr_normalize_extra_info that offers a less confusing name.

                     Changes in version 2.3.3.9001                      

- Add nmr_identify_regions_cell function
- Add documentation of HMDB_cell
- Vignettes updated
- New functions to apply multilevel statistics
- Update of README file

                        Changes in version 2.3.3                        

- Change of nmr_identify_regions_blood function
- Add nmr_identify_regions_urine function
- Add documentation of HMDB_urine
- Add computes_peak_width_ppmfunction for
nmr_integrate_peak_positions
- New get_integration_with_metadata
- Vignettes updated
- New functions to apply machine learning to proccessed datasets

                        Changes in version 2.3.2                        

- Inclusion of baseline removal using assymetric least squares
- Change the baselineThresh to NULL so it is autodetected
- Vignettes updated including baseline removal
- Bug correction in nmr_baseline_threshold
- Elimination of package vignettes (there is an error to be solved
there)
- New nmr_identify_regions function
- Add documentation of HMDB_blood
- New files_to_rDolphin function

                     Changes in version 2.3.1.9000                      

- Rename package from NIHSnmr to AlpsNMR

                        Changes in version 2.3.1                        

- Change SNR.Th value from 3 to 4 in pipeline_example.R
- Update installation instructions
- Last version form Sergio (changes not significant since 2.3.0)

                        Changes in version 2.3.0                        

- Improve installation instructions under R<3.5
- nmr_peak_detection_tune_snr function added.
- Minor bug fixes

                        Changes in version 2.2.0                        

- Improve installation instructions
- Clarify Add metadata vignette
- Add normalization diagnostics
- Add some data analysis helpers
- Enable parallellization for sample loading, peak detection and data
analysis helpers
- Do not set negative area values to zero, to avoid biasing variances
- Add metadata from a single tidy excel function
- Add nmr_diagnose to get and set diagnostic information
- Add nmr_diagnose support to nmr_normalize
- Minor bug fixes

                        Changes in version 2.1.0                        

- Documentation improvements
- nmr_dataset_peak_table object for peak detection results

                        Changes in version 2.0.0                        

- Too many changes to be listed here. Check the vignette for a
summary
of all the features. Use browseVignettes("NIHSnmr").

                        Changes in version 1.2.0                        

Breaking changes

- Rename injection_id to NMRExperiment.

- nmr_dataset_load and nmr_dataset_save now use readRDS and saveRDS
instead of load and save. This is the right approach to serialize
single R objects. If you need a script to convert previously saved
datasets (created using nmr_dataset_save) please use
NIHSnmr:::nmr_dataset_load_old_and_save("your_old_file.RData",
"your_old_file.RDS") to convert the files. Sorry for the
inconvenience, but the sooner we fix this the better.

- filter to select a subset of samples from an nmr_dataset object has
been adapted to dplyr >= 0.7.4. Unless you used the .dots argument
in your calls there is no need to change anything. This means we now
use a tidy evaluation syntax for filter.

- nmr_get_metadata() returns always a data frame / tibble, even when
only a single column is requested. It also always includes the
"NMRExperiment" column.

- nmr_dataset object has two tables metadata and metadata_ext. The
metadata_ext table includes all the metadata we add with
nmr_add_metadata while metadata has the internal metadata
(acquisition parameters, etc). Please use
nmr_get_metadata(nmr_dataset) instead of nmr_dataset$metadata.

Other changes

- Remove workaround to dplyr issue:
https://github.com/tidyverse/dplyr/issues/2203 (Sergio Oller
reported and fixed the issue, dplyr-0.7.0 is fixed)

- The Bruker title file has quite a free format definition. A title
file can contain lines like "Field value" or "Field value ;" or
simply "value". The heuristics to parse the title file have been
improved.

- Depend on tidyr 0.8.1. tidyr 0.8.0 had a bug that we reported (and
for which we also provided a fix):
https://github.com/tidyverse/tidyr/pull/419

- nmr_get_metadata gives a warning if the user asks for metadata
columns that are missing.

- New nmr_integrate_regions function.

- nmr_normalize accepts pqn normalization.

[ANCOMBC](https://bioconductor.org/packages/ANCOMBC)
-------

                 Changes in version 0.99.5 (2020-10-19)                 

- Update README

                 Changes in version 0.99.4 (2020-09-26)                 

- Update README

                 Changes in version 0.99.3 (2020-09-13)                 

- Using the GlobalPatterns dataset in the examples

                 Changes in version 0.99.2 (2020-09-05)                 

- Integrating with phyloseq-class experiment-level object

                 Changes in version 0.99.1 (2020-08-23)                 

- Addressed reviewer's comments

                 Changes in version 0.99.0 (2020-08-10)                 

- Submitted to Bioconductor

[AnnotationHub](https://bioconductor.org/packages/AnnotationHub)
-------------

                       Changes in version 2.21.0                        

BUG FIX

- (2.21.5) Fix documentation for setting AnnotationHubOptions

- (2.21.3) Fix printing of proxy when present

USER-VISIBLE MODIFICATIONS

- (2.21.6) Make internet connection test less stringent

- (2.21.4) Add link to github for reporting issues

- (2.21.2) Update to reference hubs@bioconductor.org for help

- (2.21.1) Update .tidyGRanges to account for incorrect or missing
  genomes

[AnnotationHubData](https://bioconductor.org/packages/AnnotationHubData)
-----------------

                       Changes in version 1.19.0                        

INTERNAL BUG CORRECTION

- 1.19.2 Update Metadata from Ensembl function to use
  GenomeInfoDb:::fetch_species_index_from_Ensembl_FTP instead of
  parsing the
  file path

- 1.19.1 misplaced ! clause

[annotatr](https://bioconductor.org/packages/annotatr)
--------

                       Changes in version 1.16.0                        

USER-FACING CHANGES

- Export expand_annotations(), tidy_annotations(), and
  subset_order_tbl().

BUGFIXES

- Fix incorrect shortcut search for HMMs.

[AnVIL](https://bioconductor.org/packages/AnVIL)
-----

                        Changes in version 1.2.0                        

NEW FEATURES

- (v 1.1.3) introduce .deprecated flag in operations() / tags();
don't
include deprecated APIs by default; warn on use of deprecated APIs.

- (v 1.1.4) add repositories() to return binary (if available),
Bioconductor, and CRAN repository paths.

- (v 1.1.6) provide md5sum as check on service version.

- (v 1.1.9) add avfiles_*() for managing workspace bucket files.

- (v 1.1.15) add avtable_import_set() to create subsets of tables,
following the Terra data model.

- (v 1.1.16) add avruntimes(), avworkspace_jobs() to query for
runtimes and jobs associated with the active billing account.

- (v 1.1.17) add avdisks() to query for persistent disks associate
with the active billing account.

- (v 1.1.21) add avworkflow_*() for interacting with workflow jobs
and
outputs.

[AnVILBilling](https://bioconductor.org/packages/AnVILBilling)
------------

                       Changes in version 0.0.12                        

- pass CMD check and BiocCheck

                       Changes in version 0.0.11                        

- removes browse_reck parameters, encapsulates authentication better

                       Changes in version 0.0.10                        

- New parameter to browse_reck() that defaults to NOT running bq_auth
explicitly. Now only one authentication per session is performed if
do_auth is FALSE

                        Changes in version 0.0.9                        

- Plots in plot tab are now plotly for pointwise segmented display,
and simple cumulative display

                        Changes in version 0.0.7                        

- No warnings on R CMD check or BiocCheck; one ERROR concerning
support site registration
- Includes browse_reck(), an app that tabulates expenses over time
- Updates to vignette, including demonstrative static graphics
- It seems that the page_size parameter does not properly propagate
to
BigQueryConnection; documenting this will take time; for now avoid
querying for long intervals of time

[AnVILPublish](https://bioconductor.org/packages/AnVILPublish)
------------

                       Changes in version 0.0.10                        

- Add 'best practices' and rationale for Rmarkdown-to-jupyter
notebook
conversion.

                        Changes in version 0.0.9                        

- Create a notebook '00-<<workspace name>>' to install package / book
dependencies specified in the original source.
- Don't link to vignettes from the DASHBOARD, since the namespace
changes in cloned workspaces.
- as_workspace(..., create = FALSE, update=FALSE) now evaluates code,
silently.

                        Changes in version 0.0.8                        

- Support collections of Rmd files that are not packages, e.g.,
bookdown sites.
- Add R / Bioconductor version to dashboard

                        Changes in version 0.0.7                        

- Revise Rmd-to-ipynb work flow

- Don't evaluate code chunks (avoids including output in notebook,
and side-effects because rmarkdown::render does not start a
separate process)
- Insert metadata to use the R kernel. jupytext can do this more
elegantly, but does from .md renders code chunks and
pre-formatted rather than evaluation cells, and from .Rmd does
not process markdown well enough, e.g., not suppporting
[foo][]-style links when the definition is elsewhere in the
document.

                        Changes in version 0.0.6                        

- Added a NEWS.md file to track changes to the package.

- Extensive interface renaming

- as_workspace() (formerly package_source_as_workspace())
- as_notebook() (formerly vignettes_to_notebooks())
- add_access() (formerly bioconductor_user_access())

[APAlyzer](https://bioconductor.org/packages/APAlyzer)
--------

                 Changes in version 1.3.3 (2020-07-20)                  

- Fixed hg38_REF.RData.

                 Changes in version 1.3.2 (2020-07-20)                  

- Updated the CITATION.

                 Changes in version 1.3.1 (2020-07-19)                  

- Fixed the interval issues in hg38_REF.RData and mm10_REF.RData.

- Updated the CITATION.

[artMS](https://bioconductor.org/packages/artMS)
-----

                 Changes in version 1.6.7 (2020-10-22)                  

- The following plots are now deprecated

- .clustering.log2fcSign.all-zoom.pdf

- .clustering.log2fc.all-zoom.pdf

- Fix bugs

                 Changes in version 1.6.6 (2020-10-21)                  

- Convert new MaxQuant format of PTMs to the old format

- MSstats messages are not displayed by default when using
  artmsQuantification.
  The user can enable MSstats messages by selecting "display_msstats =
  TRUE"

- Prevent artmsWriteConfigYamlFile() from overwriting an existing
  configuration
  file unless the user allows it (overwrite = TRUE)

- printPDF now available in all functions printing plots to pdf, which
  means that
  notebooks can be used and print all plots. Default is still printPDF
  = TRUE

- Fix bugs

                 Changes in version 1.6.5 (2020-05-20)                  

- Fix ggplot warnings (caused by NA values)

- Fix artmsAnalysisQuantification reproducibility plots

- Improves artmsQualityControlEvidenceBasic() correlation matrix
  clustered plot

- Fix pca01.pdf plot

- New pca04.pdf plot (dot plot)

- artmsAnalysisQuantifications check point: check if sufficient data is
  available

                 Changes in version 1.6.4 (2020-05-12)                  

- Fix Quality Control functions to handle a small number of runs (less
  than 5)

- New argument "printPDF" for artmsQualityControlSummaryExtended, to
  select
  whether to print plots to PDFs (default = TRUE)

- Vignette: example plots added

                 Changes in version 1.6.3 (2020-05-06)                  

- Bug Fixes affecting artmsAnalysisQuantifications()

                 Changes in version 1.6.2 (2020-05-05)                  

- Fix NEWS formatting

- Update vignette with AC options

                 Changes in version 1.6.1 (2020-04-29)                  

- Fix NEWS formatting

- Update vignette with AC options

[ASpli](https://bioconductor.org/packages/ASpli)
-----

                       Changes in version 1.99.3                        

NEW FUNCTIONS AND FUNCTIONALITIES

- New parameters for gbCounts and jCounts functions:

- libType: to specify single end ("SE") or paired end ("PE") library

- strandMode:
  - 0: strand of the pair is always *.
  - 1: strand of the pair is strand of its first alignment
  - 2: strand of the pair is strand of its last alignment.

                       Changes in version 1.99.1                        

NEW FUNCTIONS AND FUNCTIONALITIES

- For counting: gbCounts and jCounts. We add library type and
  strand-specificity parameters.

- For DU estimation: DUreport.norm and DUreport.offset

- For prepare DU reports: filterSignals, gbDUreport, integrateSignals,
  jDUreport, splicingReport and junctionsPJU

- For printing reports: exportIntegratedSignals, writeJDU and
  writeSplicingReport

FUNCTIONS DEPRECATED

- loadBAM: will be replaced by gbCounts. Bams files are proccesed one
  by one, according target object. It improves running time and memory
  usage

- readCounts: it is replaced by gbCounts

- AsDiscover: it is replaced by jCounts

- DUreport: it is replaced by gbDUreport and jDUreport

- mergeBinDUAS: has no direct replacement.

- plotGenomicRegions: has not direct replacement. There are new
  functions for exportiing results into HTML pages.

[AssessORF](https://bioconductor.org/packages/AssessORF)
---------

                 Changes in version 1.7.1 (2020-10-23)                  

- Genome viewer plots will now show the forward strand above the
  reverse strand, and reading frames will now be laid out in 1 through
  6 order moving down the genome viewer (before it was reversed)

                 Changes in version 1.7.0 (2020-10-21)                  

- Can now toggle interactivity of the genome viewer when plotting
  Assessment objects

- Can now zoom into an initial range with the genome viewer when
  plotting Assessment objects

- Documented the genome viewer in the plot.Assessment man page

- Added a genome viewer example to the vignette

- Fixed a bug that was causing an error when users attempted to zoom
  out ten-fold with the genome viewer plot
  # AssessORF 1.3

[ATACseqQC](https://bioconductor.org/packages/ATACseqQC)
---------

                       Changes in version 1.13.9                        

- fix the issue that plotCorrelation heatmap is scaled by row.

                       Changes in version 1.13.8                        

- throw an error if not enought nucleosome free read nor
mononucleosome reads for training.

                       Changes in version 1.13.7                        

- fix a bug introduced by matchPWM by paste ^ and $ into exclude
sequence name.

                       Changes in version 1.13.6                        

- update documentation of plotFootprints.

                       Changes in version 1.13.4                        

- fix a formular for TSSE score.

                       Changes in version 1.13.3                        

- fix a bug the after shift, the index is not changed.

                       Changes in version 1.13.2                        

- change the normalization method by library size for
factorFootprints
for user defined group samples.

                       Changes in version 1.13.1                        

- Add documentation to decribe the format of bindingSites for
factorFootprints.

[Autotuner](https://bioconductor.org/packages/Autotuner)
---------

                 Changes in version 1.3.0 (2020-09-24)                  

- Updated the required version of dependencies

- Added flag for plot peaks courtesy of John Bouranis

[BANDITS](https://bioconductor.org/packages/BANDITS)
-------

                        Changes in version 1.4.1                        

- parallel::makeCluster() updated to "sequential" strategy; default
  strategy gets stuck with RStudio on macOS with R 4.x.

[BASiCS](https://bioconductor.org/packages/BASiCS)
------

                 Changes in version 2.1.16 (2020-10-22)                 

- Updates `.BASiCS_MCMC_Start` to use different starting values when an
  empirical Bayes prior is assigned to mu

- Adds `log_scale` parameter to `.EmpiricalBayesMu`

- Fixes small bug in `.BASiCS_MCMC_Start`

- Updates unit tests linked to the change in `.BASiCS_MCMC_Start`

- Adds fixed seeds to `test_divide_and_conquer.R`

                 Changes in version 2.1.14 (2020-10-08)                 

- Add support for divide and conquer inference with
  `BASiCS_DivideAndConquer`
  function. Add support for this function to `BASiCS_MCMC`

- Add `BASiCS_MockSCE` function to create a SingleCellExperiment
  object with arbitrary dimensions.

- Fix Makevars warning relating to PKG_CXXFLAGS

                 Changes in version 2.1.13 (2020-10-05)                 

- Add PriorParam argument to `BASiCS_MCMC`

                 Changes in version 2.1.12 (2020-10-01)                 

- Remove unused DelayedArray call from tests

                 Changes in version 2.1.11 (2020-09-30)                 

- Updated documentation for `Threads` in `BASiCS_MCMC` to indicate
  default value.

                 Changes in version 2.1.10 (2020-09-30)                 

- Documents `Threads` argument.

                 Changes in version 2.1.9 (2020-09-28)                  

- Add `Threads` argument to `BASiCS_MCMC`, allowing parallelisation of
  MCMC
  updates across cells or genes using openMP.

- Reduce BatchInfo requirement for no spikes sampler to a warning, from
  an
  error.

                 Changes in version 2.1.8 (2020-09-28)                  

- Minor edit in the vignette to recommend the use of a seed for
  reproducible
  results when using `BASiCS_MCMC`.

                 Changes in version 2.1.7 (2020-09-23)                  

- Bugfix in format method for ResultVG class.

- For `BASiCS_PlotDE`, if only one plot is produced, the plot is
  returned as
  a bare ggplot object rather than using `cowplot::plot_grid`

- Added `EpsilonThreshold` argument for `BASiCS_DetectHVG` and
  `BASiCS_DetectLVG`.
  This uses a threshold on epsilon values to identify LVGs and HVGs.
  Also adds a `MinESS` argument to `BASiCS_DetectHVG` and
  `BASiCS_DetectLVG`.

- Deprecated `newBASiCS_Data` and the `format` methods for Results
  classes.
  For the latter, use `as.data.frame` instead.

- Improved documentation for BASiCS_ResultVG class

- `BASiCS_MCMC` now computes and stores ESS as an attribute of
  parameter
  matrices when sampling has completed.

- Make use of Parameter argument consistent across several methods,
  where Param
  and Which were used before.

- Extended BASiCS_DiagPlot to enable other diagnostic measures

- Moved to store molecule counts in `rowData(altExp(sce))` rather than
  `metadata(sce)`

- Remove dependency to KernSmooth

- Set some values to NA if they relate to genes captured in < 2 cells

- Add option to exclude spike-ins in `BASiCS_DenoisedCounts`.
  Also use a slightly more principled calculation internally.

- See issues #182, #39, #91, #169, #181, #178, #202, #190, #173

                 Changes in version 2.1.6 (2020-08-05)                  

- Prevents usage of `BASiCS_VarThresholdSearchVG` when using residual
  overdispersion parameters

                 Changes in version 2.1.5 (2020-08-05)                  

- Minor change on the definition of `OrderVariable` within
  `BASiCS_DetectVG`

                 Changes in version 2.1.4 (2020-08-05)                  

- Extended unit tests to verify the validity of input parameters for
  `BASiCS_DetectLVG` and `BASiCS_DetectHVG`.

- Updated documentation for `BASiCS_DetectVG` to account for new
  default values
  and the type of output (now object with class `BASiCS_ResultVG`).

- Unit tests for LVG updated as new default value for `ProbThreshold`
  is used
  whenever the EFDR calibration fails (2/3 as opposed to 0.5).

- Minor changes in `BASiCS_VarThresholdSearch` to avoid errors in the
  presence
  or absence of epsilon. This function is likely to be deprecated.

- Updates default probability threshold in `.ThresholdSearch` when EFDR
  calibration fails. This is now set to `ProbThreshold` (not 0.9)

                 Changes in version 2.1.3 (2020-07-29)                  

- New default values for `PercentileThreshold` (= NULL), `VarThreshold`
  (= NULL)
  and `ProbThreshold` (= 2/3) in `BASiCS_DetectVG`). The latter is
  compatible with
  the default values in `BASiCS_TestDE`.

- Changes default behaviour of `BASiCS_DetectLVG`, `BASiCS_DetectHVG`
  in
  accordance with the changes introduced for `BASiCS_DetectVG`.

                 Changes in version 2.1.2 (2020-07-29)                  

- Moves auxiliary functions related to HVG/LVG analysis to
  `utils_VG.R`. This
  includes `.VGPlot`, `.VGGridPlot`, `.HeaderDetectHVG_LVG` and `.VG`.

- Updated warning messages in `.HeaderDetectHVG_LVG`

- Updated rules for `EFDR` and `ProbThreshold` parameters in
  `BASiCS_DetectVG`.
  If `EFDR = NULL`, EFDR calibration is not used and posterior
  probability
  threshold is set to be equal to `ProbThreshold`. This behaviour is
  consistent
  with the rules used in `BASiCS_TestDE`.

- Creates `.CheckProbEFDR` to perform validity checks for
  `ProbThreshold` and
  `EFDR` in all tests (HVG/LVG and DE).

- Cleans code within `.ThresholdSearch` and `BASiCS_DetectVG`

- Removes `Threshold` from the output provided by `BASiCS_DetectVG`

- Removes duplicated plotting code in `BASiCS_DetectVG`. Now calls
  `BASiCS_PlotVG`.

- Changes `.VGGridPlot` and `.VGPlot` to use `ggplot2:theme_classic()`

                 Changes in version 2.1.1 (2020-07-28)                  

- Bugfix in BASiCS_DetectHVG and BASiCS_DetectLVG.
  Epsilon/Delta,Sigma values were erroneously not returned.

- More detailed colouring of differential expression plots.

- Bugfix in labels for DE plots.

[basilisk](https://bioconductor.org/packages/basilisk)
--------

                        Changes in version 1.2.0                        

- Added support for different Conda channels in the
  BasiliskEnvironment() constructor.

- Added locking to setupBasiliskEnv() for safe parallel
  construction of environments.

- Ensure that environments are always activated before use in
  useBasiliskEnv().

[basilisk.utils](https://bioconductor.org/packages/basilisk.utils)
--------------

                        Changes in version 1.2.0                        

- Migrated most environment-related functions to basilisk.

- Added locking to installConda() for safe parallel lazy Conda
  installations.

- Switched to the latest Miniconda3 installer.

[batchelor](https://bioconductor.org/packages/batchelor)
---------

                        Changes in version 1.6.0                        

- Allow regressBatches() to operate without batch= when design=
  is provided. Added d= and related options to conveniently
  perform a PCA on the ResidualMatrix.

- Added correct.all= option to all correction functions for
  consistency.

- Added a deferred=TRUE default to multiBatchPCA and its callers,
  to encourage use of deferred matrix multiplication for speed.

- Switched default PCA algorithm in multiBatchPCA to IrlbaParam.

- Added add.single= mode for endomorphic addition of correction
  results in correctExperiments().

[BayesSpace](https://bioconductor.org/packages/BayesSpace)
----------

                       Changes in version 0.99.8                        

Minor improvements and fixes

- spatialCluster() and spatialEnhance() now use a faster
implementation of the multivariate normal density that reduces
runtime by approximately 40%.

                       Changes in version 0.99.7                        

Minor improvements and fixes

- Documentation examples now use fewer iterations in order to reduce
the runtime of R CMD check.
- In qTune(), the min_rep and max_rep parameters have been replaced
with burn.in and nrep, respectively, to be consistent with
spatialCluster().

                       Changes in version 0.99.6                        

New features

- getRDS() gains a cache parameter. When TRUE, the RDS is cached
locally using BiocFileCache.

Minor improvements and fixes

- Addressed reviewer concerns
(https://github.com/Bioconductor/Contributions/issues/1624)
- Updated stop/warning/message statements to remove redundancies
and unnecessary use of paste().
- Removed inline conditional statements.
- Cache downloaded RDS in getRDS() (see above).
- spatialCluster() and spatialEnhance() handle the edge case where
only one iteration is kept after excluding burn-in.
- The coda::mcmc object returned by mcmcChain() now specifies the
thinning interval used in enhanced objects.
- spatialCluster() and spatialEnhance() now include platform-specific
defaults for the gamma parameter.
- Minor internal refactoring.

                       Changes in version 0.99.5                        

Minor improvements and fixes

- In spatialCluster() and spatialEnhance(), setting burn.in equal to
nrep now raises an error.

                       Changes in version 0.99.4                        

New features

- enhanceFeatures() now takes an nrounds parameter that corresponds
to
the same parameter in xgboost. If nrounds is set to 0, we
automatically tune the parameter using a train/test split for
improved feature prediction.
- spatialCluster() and spatialEnhance() both gain a burn.in parameter
specifying the number of MCMC iterations to exclude when aggregating
cluster labels and enhanced PCs.
- In clusterPlot(), label now accepts factors and vectors of strings,
in addition to numeric vectors or a column name in colData.
- Additional vignettes provided for reproducing the analyses of the
melanoma, dorsolateral prefrontal cortex, and squamous cell
carcinoma datasets presented in the bioRxiv manuscript.

Minor improvements and fixes

- The internal layout of subspots is now correctly oriented
(accounting for vertical flip of spot coordinates) when using
spatial plot functions on enhanced Visium data.
- In spatialEnhance(), PCs are now averaged over the MCMC iterations
(excluding the burn-in period).
- In enhanceFeatures(), negative expression is now clipped to 0.
- spatialPreprocess() now adds a boolean is.HVG column to rowData.
- In featurePlot(), additional arguments to geom_polygon() are
correctly passed through.

                       Changes in version 0.99.3                        

Minor improvements and fixes

- Updated README.md to include system requirements, additional
installation details, and link to vignette with demonstration of
package functions, per journal guidelines.

                       Changes in version 0.99.2                        

Minor improvements and fixes

- spatialEnhance() incorrectly added row offset to spot column
coordinate when generating subspot colData, and vice versa. This
resulted in subspots being reflected over y=x in spatial plots, and
has been fixed.
- Figures in the demonstration vignette have been updated with this
fix.

                       Changes in version 0.99.1                        

Minor improvements and fixes

- Removed Maintainer field from DESCRIPTION to adhere to Bioconductor
guidelines.

                       Changes in version 0.99.0                        

New features

- Initial Bioconductor submission

[BgeeCall](https://bioconductor.org/packages/BgeeCall)
--------

                        Changes in version 1.5.3                        

- Add function ```generate_slurm_indexes``` generating kallisto
  indexes on a cluster with slurm queing system

- Add function ```generate_slurm_calls``` generating expression
  calls on a cluster with slurm queing system

- Update documentation

- Add ?BgeeCall documentation with link to main functions

[BgeeDB](https://bioconductor.org/packages/BgeeDB)
------

                       Changes in version 2.14.1                        

- Solve download timeout with Windows

[bigPint](https://bioconductor.org/packages/bigPint)
-------

                 Changes in version 1.5.1 (2020-08-23)                  

plotClusters() function now has option for superimposing pairs or
whole
data using the showPairs parameter.

[BiocCheck](https://bioconductor.org/packages/BiocCheck)
---------

                        Changes in version 1.25                         

DEPRECATION ANNOUNCEMENT

- (1.25.11) R CMD BiocCheck and R CMD BiocCheckGitClone are deprecated.
  The
  recommended way to run the functions is within R.

NEW FEATURES

- (1.25.4) Check for single colon typos when using qualified imports
  pkg::foo()

- (1.25.1) Check for warning/notes on too-brief a Description: field
  (@federicomarini, #65)

- (1.25.8) Validate ORCID IDs (if any) in the DESCRIPTION file
  (@LiNk-NY, #97)

- (1.25.10) Add check for properly formatted CITATION file

- (1.25.12) Add NOTE to change dontrun to donttest

BUG FIXES

- (1.25.14) The ORCID ID check now accepts IDs with a X at the end.

- (1.25.9) All packages including infrastructure require a vignette

- Usage of donttest and dontrun in manual pages tagged with the keyword
  'internal' will no longer trigger a NOTE (@grimbough, #59)

- (1.25.7) Adding the sessionInfo at the end of the vignette (@llrs)

USER SIGNIFICANT CHANGES

- (1.25.3) Require Aurhors@R format over Author/Maintainer fields in
  the
  DESCRIPTION file. This has been upgraded to an ERROR.

- (1.25.2) Resolve
  <https://github.com/Bioconductor/BiocCheck/issues/57>: Suggest
  styler over formatR for automatic code re-formatting.

- (1.25.5) Add warning to new package versions with non-zero x version
  (@mtmorgan, #101)

[BiocParallel](https://bioconductor.org/packages/BiocParallel)
------------

                        Changes in version 1.24                         

BUG FIXES

- (v.1.23.1) bpvalidate() detects variables defined in parent
  environments; warns on use of global variables.

- (v.1.23.2) bplapply() runs gc() after each evaluation of `FUN()`, so
  that workers do not accumulate excessive memory allocations (memory
  on a per-process basis is not excessive, but cluster-wise could be).
  See
  https://github.com/Bioconductor/BiocParallel/pull/124:x

[BiocSet](https://bioconductor.org/packages/BiocSet)
-------

                         Changes in version 1.3                         

BUG FIX

- (1.3.5) Improve `GeneSetCollection` functionality.

- (1.3.3) Fix failing tests for `kegg_set()`.

- (1.3.1) Fix failing tests due to new version of GO.db.

NEW FEATURES

- (1.3.6) Made the switch from import/export functionality in
  rtracklayer
  to BiocIO.

- (1.3.4) Extending our export functionality to allow for OBO files.
  Also
  create functions to disply relationships of an OBOSet object.

- (1.3.2) Extending our import functionality to allow for OBO files.

[BiocSingular](https://bioconductor.org/packages/BiocSingular)
------------

                        Changes in version 1.6.0                        

- Migrated the ResidualMatrix class to its own ResidualMatrix
  package.

[biocthis](https://bioconductor.org/packages/biocthis)
--------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- Added bioc_style() which provides a partial Bioconductor coding
style compatible with styler.
- Added usethis-style functions use_bioc_citation(),
use_bioc_description(), use_bioc_github_action(),
use_bioc_issue_template(), use_bioc_news_md(),
use_bioc_pkg_templates(), use_bioc_readme_rmd(), use_bioc_support()
and use_bioc_vignette(). These functions provide
Bioconductor-friendly alternatives to several functions in the
usethis package.
- use_bioc_github_action() allows you to use a Bioconductor-friendly
GitHub Actions workflow for checking your Bioconductor package (or
one that depends on Bioconductor packages). Check the vignettes for
details on its features as well as the developer notes.

[biocViews](https://bioconductor.org/packages/biocViews)
---------

                       Changes in version 1.57.0                        

NEW FEATURES

- (1.57.5) New views AnnotationHubSoftware and ExperimentHubSoftware

BUG FIX

- (1.57.4) In NEWS generation, fix formatting.

- (1.57.1) In VIEWS generation fix but with vignetteTitles. When
  combinding
  different format types could potentially remove vignette titles that
  ended
  with "RNA,". Do strict start of and end of string check for
  formatting.

[biomaRt](https://bioconductor.org/packages/biomaRt)
-------

                       Changes in version 2.46.0                        

BUG FIXES

- getLDS() now detects if trying to use datasets from different Marts
  and reports this to the user.

[BioMM](https://bioconductor.org/packages/BioMM)
-----

                        Changes in version 1.5.7                        

- Improved description in the tutorial (BioMM() function returns the
  prediction scores instead of metrics)

- Updated BioMMstage2pred() and BioMM() functions to get prediction
  scores instead of performance metrics

- Updated getMetrics() to cope with the 'probability' predicted score
  as well.

                        Changes in version 1.5.4                        

- Improved description in the tutorial

                        Changes in version 1.5.2                        

- updated R functions for reporting results and plotting.

- added another omics data (gene expression) and KEGG pathway.

- updated tutorial by demonstrating BioMM with another omics data (gene
  expression) and KEGG pathway annotation.

[biosigner](https://bioconductor.org/packages/biosigner)
---------

                       Changes in version 1.17.2                        

MINOR MODIFICATION

- getMset: minor correction in documentation

[BrainSABER](https://bioconductor.org/packages/BrainSABER)
----------

                 Changes in version 0.99.0 (2020-02-14)                 

- *submitted to Bioconductor

[breakpointR](https://bioconductor.org/packages/breakpointR)
-----------

                        Changes in version 1.7.2                        

- Added new function to allow for removal of spikes in read coverage
  present in composite files.

                        Changes in version 1.7.1                        

- Added new deltaW calculator that allows to define multiplication of
  the original window size to
  use for deltaW calculations. One can invoke this by setting
  'multi.sizes' parameter.

[BRGenomics](https://bioconductor.org/packages/BRGenomics)
----------

                        Changes in version 1.1.3                        

- Increased safeguards against automatic vector coercion of matrices
and dataframes
- Increased tolerance of the precision workaround in
getStrandedCoverage() (workaround releveant to GenomicRanges issue
#39; will be obviated in next major Bioconductor release)
- Various minor efficiency & formatting changes

                        Changes in version 1.1.1                        

- Merge changes from 0.99.33 into release branch
- Bug workaround for getStrandedCoverage() and dependent methods for
getting coverage of normalized data (apparent bug in
IRanges::coverage() for weighting by normalized values)
- Bug fix in getCountsByPositions() for getting counts over an
unstranded region with expand_ranges=TRUE

[BridgeDbR](https://bioconductor.org/packages/BridgeDbR)
---------

                       Changes in version 1.99.1                        

BUG FIXES

- no more downloading of databases in the examples

                       Changes in version 1.99.0                        

BUG FIXES

- Added missing NEWS section for 1.23.1

- BridgeDb 3.0.0-SNAPSHOT (2020-10-18)

                       Changes in version 1.23.1                        

SIGNIFICANT USER-LEVEL CHANGES

- no longer supports 32bit on Windows

[CAGEr](https://bioconductor.org/packages/CAGEr)
-----

                       Changes in version 1.31.4                        

BUG FIXES

- Update end coordinates before start coordinates in the function
  `.aggregateTagClustersGR()`.  This should stop triggering
  "'width(x)' cannot contain negative integers" errors.

                       Changes in version 1.31.3                        

BUG FIXES

- Correct `.make.consensus.clusters` internal function.  This should
  stop
  triggering "Consensus clusters must not overlap with each other"
  errors.

                       Changes in version 1.31.2                        

BUG FIXES

- Allow empty `CTSS.chr` objects.

- Correct `plotInterquantileWidth()` to really use consensus clusters
  when passed the argument `clusters = "consensusClusters".

- Fix failures on `CAGEexp` objects containing only one sample.

[CAMERA](https://bioconductor.org/packages/CAMERA)
------

                       Changes in version 1.45.2                        

BUG FIXES

- Thanks to Sebastien Hutinet @shutinet for fixing the adduct mass of
  CF3COOH.
  Includes a corresponding change in unit tests. Closes #61

[Cardinal](https://bioconductor.org/packages/Cardinal)
--------

                 Changes in version 2.7.2 (2020-10-21)                  

SIGNIFICANT USER-VISIBLE CHANGES

- For 'mzAlign()', the 'ref' parameter now expects a vector
  of reference m/z-values rather than a complete spectrum

                        Changes in version 2.7.1                        

BUG FIXES

- Fixed issue where 'spatialDGMM()' would sometimes
  fail for features with singular segmentations

- Suppressed warnings on 'Mclust()' initialization
  to 'spatialDGMM()' caused by R 4.0 changes

- Fixed pixel/feature mapping in 'spatialDGMM()' metadata

[cbaf](https://bioconductor.org/packages/cbaf)
----

                 Changes in version 1.12.0 (2020-04-08)                 

New Features

- Terms are updated. Package can recongize even more cancer studies!

- Improvements for methylation analysis.

- If desired genes are entered as a vector, they are converted to a
  list without returning an error.

[cBioPortalData](https://bioconductor.org/packages/cBioPortalData)
--------------

                        Changes in version 2.2.0                        

New features

- studiesTable includes additional columns pack_build and api_build
to
indicate to the user which datasets have been successfully build as
MultiAssayExperiment objects. Users will be notified when a dataset
reported as not building is input to the cBioDataPack function.
- Add sampleIds argument to getDataByGenePanel as part of cache
re-work
- Allow more flexibility in the hostname when accessing the API with
cBioPortal (@inodb, #16)
- cBioDataPack downloads from a more robust repository (AWS S3;
@inodb, #22)
- removePackCache and removeDataCache now remove data from the user's
cache based on inputs to respective functions (cBioDataPack and
cBioPortalData)

Bug fixes and minor improvements

- Attempt to merge additional clinical data files from tarballs in
cBioDataPack.
- Switch to using read.delim instead of read_tsv internally to avoid
assigning NA to chromosome column
- Use 'PATIENT_ID' when available to determine if experiment data is
provided in the tarball files.
- Add tests using testthat
- Update and include percentages of studies successfully imported
using cBioDataPack and cBioPortalData in the documentation
- Fix read-in when identifiers are numeric instead of character
(@jucor, #27)
- Include pagination parameters in geneTable function (@xinwei-sher,
#29)

[CellaRepertorium](https://bioconductor.org/packages/CellaRepertorium)
----------------

                 Changes in version 0.99.0 (2020-09-28)                 

Submitted to bioconductor.

[CellTrails](https://bioconductor.org/packages/CellTrails)
----------

                        Changes in version 1.7.1                        

- Compatible with R 4.0.0

- Bugfixes:
  - Removed isSpike
  - Resolved error with extracting latent space information at
  'showTrajInfo'
  - Resolved inverted color palette for 'plotManifold' and 'plotMaps'
  - Resolved machine precision issue for 'findStates'

[CeTF](https://bioconductor.org/packages/CeTF)
----

                        Changes in version 1.0.2                        

- Fix netConditionsPlot function bugs

[cfDNAPro](https://bioconductor.org/packages/cfDNAPro)
--------

                       Changes in version 0.99.1                        

- Remove bugs in plotSingleGroup.R
- Documentation improvements.

                       Changes in version 0.99.0                        

- Now cfDNAPro supports bam file as input for data characterisation.
- Coding style improvements.
- Documentation improvements.
- Submitted to Bioconductor.

[ChAMP](https://bioconductor.org/packages/ChAMP)
-----

                       Changes in version 2.18.2                        

- Added parallel running for ebGSEA function. And allow it to
  return enriched gene list.

[Chicago](https://bioconductor.org/packages/Chicago)
-------

                 Changes in version 1.17.1 (2020-09-08)                 

- Fixed a small bug that affected distance function estimation in some
  cases, particularly when used with four-cutter enzymes.

[ChIPpeakAnno](https://bioconductor.org/packages/ChIPpeakAnno)
------------

                       Changes in version 3.23.12                       

- fix the bug for genomicElementDistribution when the order of
  SortedByQueryHits is incorrect.

                       Changes in version 3.23.11                       

- use seqlevelsStyle to reformat the seqlevels for annotation.

                       Changes in version 3.23.10                       

- update documentation for genomicElementDistribution

                       Changes in version 3.23.9                        

- add new function enrichmentPlot, genomicElementDistribution,
  genomicElementUpSetR, and methagenePlot to improve visualization.

                       Changes in version 3.23.8                        

- update documentation for findOverlapsOfPeaks.

                       Changes in version 3.23.7                        

- update documentation for findOverlapsOfPeaks.

                       Changes in version 3.23.6                        

- change parameter from otherCount to otherCounts to makeVennDiagram
  function.

                       Changes in version 3.23.5                        

- add plot parameter to makeVennDiagram function.

                       Changes in version 3.23.4                        

- update README file.

                       Changes in version 3.23.3                        

- use roxygen2 to generate the help file.

- move multiple package from Imports to Suggests.

                       Changes in version 3.23.2                        

- fix the issue for new paste output.

                       Changes in version 3.23.1                        

- remove dependence of Rfast

[ChIPQC](https://bioconductor.org/packages/ChIPQC)
------

                       Changes in version 1.25.1                        

- Update to maintain compatibility with DiffBind 3.0

[ChromSCape](https://bioconductor.org/packages/ChromSCape)
----------

                 Changes in version 0.99.0 (2020-10-23)                 

- Submitted to Bioconductor

[chromstaR](https://bioconductor.org/packages/chromstaR)
---------

                       Changes in version 1.14.2                        

BUG FIXES

- Corrected format check for experiment.table: Spaces are not excepted
  any longer, because they lead to downstream errors.

                       Changes in version 1.14.1                        

BUG FIXES

- Compatibility update: Replaced class() checks with is().

[CHRONOS](https://bioconductor.org/packages/CHRONOS)
-------

                       Changes in version 1.17.2                        

- Introduced the use of rJava in getLinearSubpath().

- Introduced futures in class LinearPaths.

[circRNAprofiler](https://bioconductor.org/packages/circRNAprofiler)
---------------

                        Changes in version 1.3.1                        

Added files to fix bugs

[cleanUpdTSeq](https://bioconductor.org/packages/cleanUpdTSeq)
------------

                       Changes in version 1.27.1                        

- update the classifier dataset

[cleaver](https://bioconductor.org/packages/cleaver)
-------

                 Changes in version 1.27.1 (2020-06-05)                 

- Fix cleavage rule for glutamyl endopeptidase in documentation
  Thanks to Mariia Chernigovskaya
  <mariia.chernigovskaya@medisin.uio.no> for
  reporting this error.

[clusterProfiler](https://bioconductor.org/packages/clusterProfiler)
---------------

                       Changes in version 3.17.5                        

- update [[.compareClusterResult (2020-10-14, Wed)

                       Changes in version 3.17.3                        

- internal suports of enrichment analyses using WikiPathways
(2020-09-09, Wed)
- enrichWP for ORA analysis
- gseWP for GSEA analysis
- get_wp_organisms for listing supported organisms
- read.gmt.wp for parsing gmt file downloaded from wikiPathways

                       Changes in version 3.17.2                        

- use libcurl if capable (2020-09-08, Tue)
- https://github.com/YuLab-SMU/clusterProfiler/pull/290

                       Changes in version 3.17.1                        

- bug fixed of extract_params (2020-08-18, Tue)
- https://github.com/YuLab-SMU/clusterProfiler/issues/282

[clustifyr](https://bioconductor.org/packages/clustifyr)
---------

                 Changes in version 1.1.2 (2020-09-21)                  

- USCS cell browser reference building

- Tutorial update

- Bug fixes

                 Changes in version 1.1.0 (2020-05-21)                  

- Bioc release

- Bug fixes

[CNVPanelizer](https://bioconductor.org/packages/CNVPanelizer)
------------

                       Changes in version 1.21.2                        

- Moving BiocGenerics at DESCRIPTION from Suggests to Imports

[cola](https://bioconductor.org/packages/cola)
----

                        Changes in version 1.5.5                        

- add `enforce` argument in `get_signatures()`.

- `subset` can be a vector of indices in
  `consensus_partition_by_down_sampling()`.

                        Changes in version 1.5.3                        

- `predict_classes()` is speeded up 2x.

- `get_signatures()`: add `top_signatures` argument to control the
  number of
  top signatures.

                        Changes in version 1.5.2                        

- add a `DownSamplingConsensusPartition` class and corresponding
  methods.

- add back `HierarchicalPartition` class and corresponding methods.

                        Changes in version 1.5.1                        

- add a `predict_classes()` function.

- add the cola analysis for Golub dataset as a data object in the
  package.

- automatically install the "suggested" packages.

[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap)
--------------

                        Changes in version 2.5.6                        

- `ht_shiny()`: add argument `app`.

- `grid.dendrogram()`: change the recursive implementation with
  iterations.

- change default raster device to `CairoPNG`.

- `Heatmap()`: If the discrete `col` covers more than the levels in the
  matrix,
  the full color set is still saved, which means, in
  `heatmap_legend_param` you
  can set `at` that are not all in the matrix but are in the `col`.

- padding of the whole plot and spaces of column titles are adjusted to
  fit ggplot2

- add `row_gap` and `column_gap` in `Legend()`.

- `oncoPrint()`: now draw legends the same as `alter_fun`.

- add a new function `attach_annotation()`.

- legends for row annotations can be grouped with column annotation
  legends.

- annotation name allows rotations.

                        Changes in version 2.5.5                        

- still draw the legend when all values are NA in an annotation.

- add `show_fraction` argument in `anno_oncoprint_barplot()` function
  to show the fractions
  of mutations instead of the counts.

- `pheatmap()`: improve the setting of `color` and `breaks`.

- `ht_opt$TITLE_PADDING` can be set with a unit of length two.

- `HeatmapAnnotation()`: remove colors that are not in the annotations.

- `pheatmap()`: fixed a bug when length(breaks) = length(color) + 1

- `pheatmap()`: legend breaks are centered to zero if the matrix is
  scaled.

- `pheatmap()`: color mapping is symmetric to zero when scale is set.

- support ragg package to write temporary png files

- `densityHeatmap()`: column dendrogram is reordered by column means
  for ks method.

                        Changes in version 2.5.4                        

- fixed a bug where slice clusters were wrongly reordered.

- `Heatmap()`: add `border_gp` argument.

- Legends are nicely placed.

- `anno_block()`: allows to set height and width.

- support better rasterization.

- support setting graphics on dendrogram nodes.

- Add a new vignette "interactive heatmap"

- `Legends()`: fixed a bug of mixtype "legend" to "Legend".

- now assign correct envir to `decorate_dend()`.

- `pheatmap()`: check `NA` in the matrix.

- `grid.dendrogram()`: consider branches with height zero.

- checking the dimension of the matrix and the nobs of annotations when
  adding them.

                        Changes in version 2.5.3                        

- add `selectArea()`/`selectPosition()` which allows interactively
  select a region from
  the heatmaps.

- export the heatmap as a shiny app!!!

- `col` in `Heatmap()` accpets a `ColorMapping` object.

- `default_col()`: print a message if there are outliers in the matrix.

- `discrete_legend_body()`: adjust ncol and nrow if there are empty
  rows and columns in the layout.

- `anno_image()`: fixed a bug that images are not reordered.

- `anno_mark()`: now expression is correctly supported.

- `anno_zoom()`: order of index in `panel_fun` is adjusted to the order
  in the heatmap

- `list_to_matrix()`: convert elements to characters.

- print messages for `anno_mark()`, `anno_zoom()`, `draw_legend()` (if
  legends are wrapped)
  if working under RStudio.

                        Changes in version 2.5.2                        

- translate pheatmap to Heatmap

- `upset_top_annotation()` and `upset_right_annotation()`: the names of
  the annotations
  are changed to `intersection_size`, `set_size` and `union_size`.

- `list_components()`: adds `pattern` argument.

                        Changes in version 2.5.1                        

- A temporary solution of the sum of two complicated units (in temp.R).

[CoreGx](https://bioconductor.org/packages/CoreGx)
------

                        Changes in version 1.1.5                        

- Implemented a new class, the LongTable, to store the results of a
treatment response experiment. This class provides a flexible and
fast data storage object which can be subclassed for use in other R
packages.
- Added vignette documenting LongTable accessors and usage of the new
object.

                        Changes in version 1.0.2                        

- Bug fix: suppress warnings thrown by piano::runGSA inside the
connectivitScore function

                        Changes in version 1.0.1                        

- Updated the CoreGx vignette to include more information on
extending
the CoreSet class for use in other treatment-response experiments.

[CRISPRseek](https://bioconductor.org/packages/CRISPRseek)
----------

                       Changes in version 1.29.2                        

- added parameters such as predIndelFreq to allow the prediction of
  indels and their frequecies for Cas9 targeted sites

                       Changes in version 1.29.1                        

- added parameter calculategRNAefficacyForOfftargets, default to TRUE.

[csaw](https://bioconductor.org/packages/csaw)
----

                       Changes in version 1.24.0                        

- Accept a list of BamFiles in bam.files= for all functions.

[cTRAP](https://bioconductor.org/packages/cTRAP)
-----

                         Changes in version 1.8                         

Interactive functions for loading data and analysing results

- New Shiny-based graphical interface functions:
- launchDiffExprLoader(): load differential expression data
- launchCMapDataLoader(): load CMap data
- launchResultPlotter(): view and plot data results
- launchMetadataViewer(): check metadata of a given object

Major changes

- downloadENCODEknockdownMetadata(): metadata is automatically saved
to a file in order to avoid downloading metadata every time this
function is run
- plotTargetingDrugsVSsimilarPerturbations():
- automatically look for matching compounds in multiple columns of
both datasets
- allow to manually select columns on which to merge datasets
- prepareDrugSets(): drug sets based on numeric molecular descriptors
are now prepared using evenly-distributed intervals
- Simplify tutorial

                        Changes in version 1.6.1                        

- listExpressionDrugSensitivityAssociation() lists available gene
expression and drug sensitivity associations
- First argument of rankSimilarPerturbations() and
predictTargetingDrugs() changed name from diffExprGenes to input and
now accepts:
- Named numeric vector containing differential gene expression
values with gene symbols as names, as before;
- Character vector containing a custom gene set to test for
enrichment (only to use with GSEA).
- In rankSimilarPerturbations() and predictTargetingDrugs(), when
performing gsea method, allow to set different gene set size for top
up- and down-regulated genes with geneSize argument:
- e.g. geneSize=c(100, 200) creates gene sets from the top 100 up-
and top 200 down-regulated genes
- using geneSize=c(150, 150) or geneSize=150 is equivalent
- Plotting:
- plot() now supports plotting predictTargetingDrugs() results for
a given drug, e.g. plot(targetingDrugs, "1425")
- plot() nows allows to set plot title with argument title
- plot() now plots results based on available methods instead of
trying to plot based on results from spearman method only
- GSEA plots now support two or less gene hits
- GSEA plots now support plotting of multiple perturbations
- GESA plots now show the first and last values of ranked genes
- plotDrugSetEnrichment() now returns a list whose names are drug
set names
- as.table() improvements:
- Return cell identifiers and gene information (if available and
as needed)
- Support predictTargetingDrugs() results
- Return results ordered as found on input

Bug fixes and minor changes

- downloadENCODEknockdownMetadata() now correctly retrieves metadata
following a change in the metadata content from ENCODE
- Fix bugs when rendering GSEA plots due to deprecated functions in
ggplot2
- Improve tutorial
- Copy-edit CMap-related console messages
- Copy-edit function documentation

[customCMPdb](https://bioconductor.org/packages/customCMPdb)
-----------

                 Changes in version 1.1.0 (2020-10-02)                  

- Initial version

                 Changes in version 0.99.5 (2020-10-09)                 

- Addressed hard links in functions like "processDrugage"

                 Changes in version 0.99.0 (2020-09-30)                 

- Consolidated two packages into one named as customCMPdb
- Supported customizing compound annotation database
- Add customized annotation table
- Delete
- List and set default
- Updated query interface
- Supported selecting annotation resources
- Supported querying added customized annotations

[cytomapper](https://bioconductor.org/packages/cytomapper)
----------

                 Changes in version 1.1.6 (2020-10-23)                  

- Dropped 32-bit Windows support

                 Changes in version 1.1.5 (2020-10-11)                  

- Prepared for Bioc 3.12 release

- Started unit testing the shiny app

                 Changes in version 1.1.4 (2020-09-13)                  

- Allow channel-specific inputRange inputs for normalisation

                 Changes in version 1.1.3 (2020-09-12)                  

- Extended vignette

- Changed package title

                 Changes in version 1.1.2 (2020-07-17)                  

- Added shiny app to package

                        Changes in version 1.1.1                        

- Images are no longer re-normalized after channel merging (> 3
  channels)

- Instead images are clipped at 1 leading to brighter colours

- The same happens for colour merging when colouring masks by feature
  expression

[CytoML](https://bioconductor.org/packages/CytoML)
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

[CytoTree](https://bioconductor.org/packages/CytoTree)
--------

                 Changes in version 0.99.6 (2020-06-19)                 

- Remove INFO/WARNING/ERROR tags

- Add vignettes

                 Changes in version 0.99.5 (2020-06-18)                 

- Re-build this package

                 Changes in version 0.99.4 (2020-06-05)                 

- Change examples and provide use case

                 Changes in version 0.99.3 (2020-06-02)                 

- Update for the comments from the reviewer

                 Changes in version 0.99.2 (2020-05-24)                 

- Update R version to 4.0

                 Changes in version 0.99.1 (2020-05-10)                 

- Fixed some warnings in BiocCheck

                 Changes in version 0.99.0 (2020-05-10)                 

- First commit

[dagLogo](https://bioconductor.org/packages/dagLogo)
-------

                       Changes in version 1.27.9                        

- update the hyperlink of p.adjust.methods in documentation.

                       Changes in version 1.27.8                        

- add adjust p-value for testDAU function.

                       Changes in version 1.27.7                        

- optimize the label position of markers for logo.

                       Changes in version 1.27.6                        

- add markers for logo.

                       Changes in version 1.27.5                        

- allow multiple species for prepareProteome.

                       Changes in version 1.27.3                        

- fix the typo in dagLogo documentation.

                       Changes in version 1.27.2                        

- add availableSchemes function.

- fix the x-axis.

                       Changes in version 1.27.1                        

- adjust the Depends, Imports and Suggests packages

[DAMEfinder](https://bioconductor.org/packages/DAMEfinder)
----------

                        Changes in version 1.1.3                        

- Add CITATION

- Remove "bad chromosome" filtering

- Argument `build` in extract_bams() removed

                        Changes in version 1.1.2                        

- Remove vcfR dependency and add VariantAnnotation

                        Changes in version 1.1.1                        

- Fix typo in split_bams that excluded some chromosomes

[dasper](https://bioconductor.org/packages/dasper)
------

                       Changes in version 0.99.2                        

NEW FEATURES

- Merge documentation into one man page for junction, coverage and
outlier processing functions to reduce runtime of roxygen examples.

                       Changes in version 0.99.1                        

NEW FEATURES

- Change outlier_detect() to using basilisk for interfacing into
python replacing reticulate.

                       Changes in version 0.99.0                        

NEW FEATURES

- Converted dasper into a Bioconductor-friendly format using
biocthis.
- Added junction_load(), which loads raw junction data from
RNA-sequencing into an RangedSummarizedExperiment object. Includes
an option to allow download of user-specified control junctions.
- Added junction_annot(), which uses information from reference
annotation and the strand of a junction to classify junctions as
"annotated", "novel_acceptor", "novel_donor", "novel_exon_skip",
"novel_combo", "ambig_gene" and "unannotated".
- Added junction_filter(), which filters junctions by their count,
width, annotation or if they overlap a set of user-defined regions.
- Added junction_norm(), which normalises raw junction counts (into a
proportion-spliced-in) by dividing the counts of each junction by
the total number of counts in it's associated cluster.
- Added junction_process(), a wrapper function for all "junction_"
prefixed functions except junction_load().
- Added junction_score(), which scores patient junctions based on the
extent their counts deviate from a control count distribution.
- Added coverage_norm(), which will load and normalise coverage for
exonic/intronic regions corresponding to each junction.
- Added coverage_score(), which scores coverage associated with each
junction based on it's deviation from control coverage
distributions.
- Added coverage_process(), a wrapper function for all "coverage_"
prefixed functions.
- Added outlier_detect(), which uses the junction scores and coverage
scores as input into an unsupervised outlier detection algorithm to
find the most outlier-looking junctions in each sample.
- Added outlier_aggregate(), which aggregates the junction-level
outlier data to a cluster-level.
- Added outlier_process(), a wrapper function for all "outlier_"
prefixed functions.
- Added plot_sashimi(), which enables the visualisation of junction
data across genes/transcripts or regions of interest.

[decompTumor2Sig](https://bioconductor.org/packages/decompTumor2Sig)
---------------

                 Changes in version 2.4.1 (2020-07-27)                  

- Removed dependency of package CRAN vcfR (archived on 2020-07-05),
  using
  functions of Bioconductor package VariantAnnotation instead

- Improved the mutation filtering so that multiallelic SNVs aren't
  excluded
  when loading tumor genomes from a VCF file

- Updated citation and affiliation information

- Added consistency check for reference genome and genome annotation

- Improved error messages

[deepSNV](https://bioconductor.org/packages/deepSNV)
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

[DegNorm](https://bioconductor.org/packages/DegNorm)
-------

                       Changes in version 0.99.13                       

- bug fix and version bump .

                       Changes in version 0.99.12                       

- version bump.

                       Changes in version 0.99.11                       

- fix a bug in plot_coverage function and make changes in color scheme.

                       Changes in version 0.99.10                       

- Made additional required changes from review

                       Changes in version 0.99.9                        

- Made most required changes from review

                       Changes in version 0.99.8                        

SUBMIT TO BIOCONDUCTOR  FOR REVIEW

- Corrected warnings and errors

[DEGreport](https://bioconductor.org/packages/DEGreport)
---------

                       Changes in version 1.25.1                        

- Fix: load dplyr::n() to avoid error

[DelayedArray](https://bioconductor.org/packages/DelayedArray)
------------

                       Changes in version 0.16.0                        

NEW FEATURES

- Added 'as.sparse' argument to read_block() (see ?read_block) and to
  AutoRealizationSink() (see ?AutoRealizationSink).

- SparseArraySeed objects now can hold dimnames. As a consequence
  read_block() now also propagates the dimnames to sparse blocks,
  not just to dense blocks.

- Matrix multiplication is now sparse-aware via sparseMatrices.

- Added is_sparse<- generic (with methods for HDF5Array/HDF5ArraySeed
  objects only, see ?HDF5Array in the HDF5Array package).

- Added viewportApply() and viewportReduce() to the blockApply()
  family.

- Added set_grid_context() for testing/debugging callback functions
  passed
  to blockApply() and family.

SIGNIFICANT USER-VISIBLE CHANGES

- Renamed first write_block() argument 'x' -> 'sink'

- Renamed:
  RealizationSink() -> AutoRealizationSink()
  get/setRealizationBackend() -> get/setAutoRealizationBackend()
  blockGrid() -> defaultAutoGrid()
  row/colGrid() -> row/colAutoGrid()

- Improved support of sparse data:
  - Slightly more efficient coercion from SparseArraySeed to
  dgCMatrix/lgCMatrix (small speedup and memory footprint reduction).
  This provides a minor speedup to the sparse aware block-processed
  row/col summarization methods for DelayedMatrix objects when the
  object is sparse. (These methods are: row/colSums(), row/colMeans(),
  row/colMins(), row/colMaxs(), and row/colRanges(). The methods
  defined
  in DelayedMatrixStats are not sparse aware yet so are not affected.)
  - Made the following block-processed operations on DelayedArray
  objects
  sparse aware: anyNA(), which(), max(), min(), range(), sum(), prod(),
  any(), all(), and mean(). With a typical 50%-60% speedup when the
  DelayedArray object is sparse.
  - Implemented a bunch of methods to operate natively on
  SparseArraySeed
  objects. Their main purpose is to support the above i.e. to support
  block processed methods for DelayedArray objects like sum(), mean(),
  which(), etc... when the object is sparse. Note that more are needed
  to also support the sparse aware block-processed row/col
  summarization
  methods for DelayedMatrix objects so we can finally ditch the costly
  coercion from SparseArraySeed to dgCMatrix/lgCMatrix that they
  currently
  rely on.

- The utility functions for retrieving grid context for the current
  block/viewport should now be called with no argument (previously
  one needed to pass the current block to them). These functions are
  effectiveGrid(), currentBlockId(), and currentViewport().

- DelayedArray now depends on the MatrixGenerics package.

BUG FIXES

- Various fixes and improvements to block processing of sparse logical
  DelayedMatrix objects (e.g. DelayedMatrix object with a lgCMatrix
  seed from thr Matrix package).

- Fix extract_sparse_array() inefficiency on dgCMatrix and lgCMatrix
  objects.

- Switch matrix multiplication to bplapply2() from bpiterate() to fix
  error handling.

[DelayedMatrixStats](https://bioconductor.org/packages/DelayedMatrixStats)
------------------

                       Changes in version 1.12.0                        

- Dispatch to sparseMatrixStats for sparse seeds that do not have
  their own methods (<URL:
  https://github.com/PeteHaitch/DelayedMatrixStats/pull/65>).

- Fix center= handling for all affected functions (<URL:
  https://github.com/PeteHaitch/DelayedMatrixStats/pull/65>).

- DelayedMatrixStats now imports the generics from
  MatrixGenerics. Thanks to Aaron Lun resolving this (<URL:
  https://github.com/PeteHaitch/DelayedMatrixStats/pull/62>).

[DeMixT](https://bioconductor.org/packages/DeMixT)
------

                        Changes in version 1.6.0                        

- 
  MPI for parallel computing is avaliable under R 4.0.0 for linux
  and Mac OS platforms.

- 
  Gene expression data of normal tissues (Lung, Prostate and
  Thyroid) from the GTEx study are included.

- 
  Rename DeMixT_S1 function to DeMixT_DE.

[densvis](https://bioconductor.org/packages/densvis)
-------

                 Changes in version 0.99.0 (2020-09-24)                 

- First version of the package.

[DepecheR](https://bioconductor.org/packages/DepecheR)
--------

                 Changes in version 1.5.4 (2020-09-17)                  

- Correction of the information about the content of the output from
  depeche.

                 Changes in version 1.5.3 (2020-07-02)                  

- Correcting the p-adjustments, so that it in fact uses
  Benjamini-Hochberg, and

- not the more conservative Hochberg, as default.

                 Changes in version 1.5.2 (2020-06-05)                  

- Introducing samplingSubset in the depeche function

- Bug fixes in dOptPenalty

- Small text updates to main vinjette and examples, without code
  implications.

                 Changes in version 1.5.1 (2020-05-18)                  

- Introduction of neighSmooth - a generalization of groupStatPlot.

- Bug fix and simplification of dOptPenalty termination criteria.

- Re-tidying of the code base.

[DESeq2](https://bioconductor.org/packages/DESeq2)
------

                       Changes in version 1.30.0                        

- Major overhaul of dispersion estimation and GLM estimation
  functions from Constantin Ahlmann-Eltze, which will allow use of
  the glmGamPoi package from within DESeq2, in particular relevant
  for single-cell datasets.  DESeq() can be directed to use
  glmGamPoi for dispersion and GLM fitting by specifying
  fitType="glmGamPoi". The glmGamPoi estimation is much faster
  than original DESeq2 estimation for single-cell datasets,
  e.g. for 30,000 cells, calling glmGamPoi was 13x faster than
  original DESeq2. In addition, the dispersion estimation is more
  accurate for genes with many small counts, as found in
  single-cell datasets.
  See glmGamPoi manuscript for details on methods,
  doi: 10.1101/2020.08.13.249623.

- Added integrateWithSingleCell(), written by Kwame Forbes,
  which directs user to a menu of single-cell datasets
  available on Bioconductor and downloads/loads the one
  chosen by the user for further analysis visualization.
  (Interactive only)

[DEWSeq](https://bioconductor.org/packages/DEWSeq)
------

                 Changes in version 1.3.0 (2020-09-30)                  

- modified authors and description sections in DESCRIPTION

[DIAlignR](https://bioconductor.org/packages/DIAlignR)
--------

                        Changes in version 1.1.5                        

- ananlyteFDR to limit features for multipeptide.

- Removed decoy from features.

- BugFix alignment if all intensities are zero.

- Fixed printed stats.

[DiffBind](https://bioconductor.org/packages/DiffBind)
--------

                         Changes in version 3.0                         

- 
  This is a major release of DiffBind and the new version will be
  3.0.
  
  The main upgrade involves how the modelling and normalization
  are done.  DiffBind now supports models and contrasts of
  arbitrary complexity using either/both DESeq2 or/and edgeR, as
  well as a myriad of normalization options.
  
  NB:
  
  The previous methods for modelling are maintained for backward
  compatibility, however they are not the default.  To repeat
  earlier analyses, dba.contrast() must be called explicitly with
  design=FALSE. See ?DiffBind3 for more information.
  
  The default mode for dba.count() is now to center around
  summits (resulting in 401bp intervals).  To to avoid
  recentering around summits as was the previous default, set
  summits=FALSE (or to a more appropriate value).
  
  Normalization options have been moved from dba.analyze() to the
  new interface function dba.normalize(). Any non-default
  normalization options must be specified using dba.normalize().

- Summary of Changes:
  
  • dba.analyze():
  
  • Automatic mode can start at any point, including from a
  sample sheet, and continue default analysis
  
  • Remove normalization options bSubControl,
  bFullLibrarySize, filter, and filterFun from
  dba.analyze(), now set in dba.normalize().
  
  • Add support to analyze using full model design formula.
  
  • Update DESeq2 analysis method.
  
  • Update edgeR analysis method.
  
  • Moved edgeR bTagwise parameter to $config option
  
  • Remove support for DESeq analysis method.
  
  • Add ability to retrieve DESeq2 or edgeR object
  
  • dba.contrast():
  
  • Add design parameter to set design formula
  
  • Add contrast parameter to specify contrast using
  variety of methods
  
  • Add reorderMeta parameter to set factor value order
  
  • Add bGetCoefficients parameter to get design
  coefficient names to use in contrasts
  
  • NEW: dba.normalize():
  
  • Support TMM, RLE, and Library size noramlization for
  both DESeq2 and edgeR
  
  • Support background bin normalization using csaw
  
  • Support offsets and loess fit normalization
  
  • Support spike-in normalization with combined or
  separate reference genomes
  
  • Support parallel factor normalization
  
  • bSubControl default depend on presence of Greylist
  
  • dba.count():
  
  • Change default to summits=250; to avoid recentering
  around summits,must set summits=FALSE
  
  • Default for bUseSummarizeOverlaps in dba.count is now
  TRUE
  
  • Automatically detect single-end/paired-end in dba.count
  
  • Automatically index unindexed bams in dba.count and
  dba.normalize
  
  • move bSubControl parameter
  
  • Default score is now new score DBA_SCORE_NORMALIZED
  
  • Add minCount parameter to dba.count(), default now 0
  instead of 1
  
  • Filtering peak by read count thresholds only available
  in dba.count()
  
  • Fix bug in dba.count() with user-supplied peakset and
  summits=n
  
  • NEW: dba.blacklist():
  
  • Apply ENCODE blacklist
  
  • Automatically detect reference genome for blacklist
  
  • Apply Greylists
  
  • Generate Greylists from controls using GreyListChIP
  package
  
  • Plotting changes:
  
  • Add loess fit line to dba.plotMA()
  
  • Add ability in dba.plotMA() to plot aribitrary samples
  (without contrast).
  
  • Add mask parameter to dba.plotBox()
  
  • Support negative scores, eg Fold changes in
  report-based objects, to enable fold-change heatmaps.
  
  • Removed bCorPlot as a parameter to dba(), dba.count(),
  and dba.analyze(). Use config.
  
  • dba.show() / print changes:
  
  • Updated dba.show() and print() to deal with designs and
  different contrast types
  
  • Add ability to retrieve design formula in dba.show()
  
  • Removed bUsePval parameter in dba.show()
  
  • Added constant variable DBA_READS to access library
  sizes
  
  • Vignette and help pages:
  
  • Replace multi-factor analysis section
  
  • Add extensive normalization section
  
  • Add blacklist/greylist section.
  
  • Add pike-in and parallel normalization examples
  
  • Add DiffBind3 help page and vignette section with
  information on backward compatibility.
  
  • Update technical details sections
  
  • General updates to all sections
  
  • Add GenerateDataFiles.R to package
  
  • Various bugfixes and cosmetic changes.

[Director](https://bioconductor.org/packages/Director)
--------

                 Changes in version 1.15.1 (2020-10-04)                 

- Updated DESCRIPTION and NEWS formatting.

[DiscoRhythm](https://bioconductor.org/packages/DiscoRhythm)
-----------

                 Changes in version 1.5.6 (2020-07-22)                  

- Visual change to screen plot (No longer apply alpha to "after").

- Bug fixes.

[distinct](https://bioconductor.org/packages/distinct)
--------

                        Changes in version 1.0.3                        

- FC and log2-FC calculation added (via 'log2_FC' function);

- 'top_results' function edited to sort results by p-value and log2-FC.

                        Changes in version 1.0.1                        

- substantial speed-up in `distinct_test` (~7 times faster);

- parallel computing introduced;

- modelling of covariates introduced, via a design matrix;

- two new functions to visualize results: `plot_cdfs` and
  `plot_densities`.

[dittoSeq](https://bioconductor.org/packages/dittoSeq)
--------

                         Changes in version 1.2                         

- Added 3 New Visualization Functions, `dittoDotPlot()`,
  `dittoDimHex()` & `dittoScatterHex()`.

- Expanded SummarizedExperiment compatibility across the entire
  toolset.

- Added ComplexHeatmap integration to `dittoHeatmap()`,
  controlled by a new input, `complex`.

- Added Rasterization for improved image editor compatibility of
  complex plots. (See the dedicated section in the vignette for
  details.)

- Added `labels.split.by` input & `do.contour`, `contour.color`,
  and `contour.linetype` inputs to scatter/dim-plots.

- Added `order` input to scatter/dim-plots for control of
  plotting order.

- Added `metas` input for displaying such data with
  `dittoHeatmap()`.

- Added `adjustment` input to `meta()`, which works exactly as in
  `gene()` (but this is not yet implemented within data grab of
  visualiation functions).

- Added `adj.fxn` input to `meta()` aand `gene()` for added
  control of how data might be adjusted (but this is not yet
  implemented within data grab of visualiation functions).

- Replaced (deprecated) `highlight.genes` input with
  `highlight.features` in `dittoHeatmap()`.

- Replaced (deprecated) `OUT.List` input with `list.out` for all
  `multi_*` plotters.

[DMCFB](https://bioconductor.org/packages/DMCFB)
-----

                        Changes in version 1.3.1                        

New Features

- Parallel reading of data is added to readBismark-method.

Removed

- export.bed is removed due to warning. Will be modified and added in
the next release.

[DMCHMM](https://bioconductor.org/packages/DMCHMM)
------

                       Changes in version 1.11.1                        

CHANGES

- Parallel reading is added to readBismark-method.

BUG FIX

- Some bugs are fixed in methHMMCMC-method to avoid creating
infinity.

[DOSE](https://bioconductor.org/packages/DOSE)
----

                       Changes in version 3.15.4                        

- update setReadable and geneInCategory methods for
compareClusterResult object (2020-10-12, Mon)

                       Changes in version 3.15.3                        

- allow passing additional parameters to fgsea (2020-10-09, Fri)
- https://github.com/YuLab-SMU/DOSE/pull/40
- add termsim and method slots to compareClusterResult, enrichRestul
and gseaResult
- https://github.com/YuLab-SMU/DOSE/pull/39

                       Changes in version 3.15.2                        

- update NCG and DGN data (2020-10-09, Thu)

[drawProteins](https://bioconductor.org/packages/drawProteins)
------------

                        Changes in version 1.9.0                        

- Default addition of chain where it is absent.

[DropletUtils](https://bioconductor.org/packages/DropletUtils)
------------

                       Changes in version 1.10.0                        

- Migrated downsampleMatrix() to scuttle with a re-export.

- Added features= to downsampleReads() for per-feature-set
  downsampling.

- Added matrix support for y= and ambient= in maximumAmbience().

- Added controlAmbience() for easy estimation of ambient
  contamination with control features.

- Added removeAmbience() function to remove the ambient solution
  from a count matrix, mostly for aesthetics.

- Report library index and feature type in output of
  read10xMolInfo().

- Support subsetting by library index/type in functions that use
  the molecule information file, such as swappedDrops() and
  chimericDrops().

- Added by.rank= option to estimateAmbience() and emptyDrops(),
  for estimation of the ambient profile by excluding barcodes
  with the largest totals.

- Added exclude.from= option to barcodeRanks(), to avoid problems
  with instability at low ranks for knee/inflection calculations
  (contributed by Stefano Mangiola).

- Minor bugfix in barcodeRanks() calculation of the knee point.
  Note that this affects the default choice of retain= in
  emptyDrops().

- Split off HTO ambience inferences into a separate
  inferAmbience() function.

- Added support for combinatorial barcodes in hashedDrops().

[easyRNASeq](https://bioconductor.org/packages/easyRNASeq)
----------

                       Changes in version 2.25.1                        

- Ported changes from 2.24.1

                       Changes in version 2.24.1                        

- Final documentation fix and removal of defunct RangedData unit tests.

[EDASeq](https://bioconductor.org/packages/EDASeq)
------

                       Changes in version 2.23.1                        

- Removed coercion methods to `CountDataSet` following `DESeq`
  deprecation.

- Updated vignette to show how to use `EDASeq` with `DESeq2`.

[edgeR](https://bioconductor.org/packages/edgeR)
-----

                       Changes in version 3.32.0                        

- 
  cpm.default() and rpkm.default() now accept offset.

- 
  scaleOffset() now accepts CompressedMatrix offset and accounts
  for norm.factors.

- 
  Revise the lowess trend fitting in voomLmFit() to downweight
  genes with exact zeros and hence fewer df to estimate the
  variance.

- 
  Add as.data.frame method for DGEList class.

- 
  Change default choice for refColumn in calcNormFactors() with
  method="TMMwsp". The new method chooses the column with the
  largest sum of sqrt-counts.

- 
  processAmplicons() can now accommodate data from newer screens
  that use a staggered primer design.

- 
  Fixed a bug that diffSpliceDGE() accept more than one coef. It
  now gives a warning if more than one coef or contrast is
  supplied. It only uses the first.

                       Changes in version 3.30.2                        

- 
  New function voomLmFit() that combines the limma voom-lmFit
  pipeline with loss of residual df due to zero counts as for
  glmQLFit(). The new function is more robust to zero counts than
  running voom() and lmFit() separately. The new function allows
  sample quality weights and intra-block correlations to be
  estimated it incorporates the functionality of
  duplicateCorrelation() and voomWithQualityWeights() as well.

- 
  New function SE2DGEList() to convert a SummarizedExperiment
  object into a DGEList object.

- 
  S3 methods for SummarizedExperiment objects are added to the
  following functions: aveLogCPM(), calcNormFactors(), cpm(),
  cpmByGroup(), estimateDisp(), filterByExpr(), glmFit(),
  glmQLFit(), plotMD(), plotMDS(), predFC(), rowsum(), rpkm(),
  rpkmByGroup() and sumTechReps().

- 
  New cpm and rpkm methods for DGEGLM and DGELRT objects.

- 
  New function effectiveLibSizes() to extract normalized library
  sizes from an edgeR data object or fitted model object.

- 
  Add as.data.frame methods for DGEExact and DGELRT objects and
  remove the 'optional' argument from as.data.frame.TopTags().

- 
  readBismark2DGE() now forces 'files' to be character vector.

- 
  Add warning messages when filterByExpr() is used without
  specifying group or design.

- 
  Add warning message when calcNormFactors() is applied to
  DGEList object containing an offset matrix.

- 
  Rewrite User's Guide Section 3.5 on Multilevel Experiments so
  that the code is valid regardless of the number of subjects in
  each disease group.

[EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano)
---------------

                         Changes in version 1.8                         

- added functionality to encircle variables of interest

- added option to remove arrowheads on connectors

- added option to rasterise images via ggrastr::geom_point_rast
  (Benjamin Ostendorf)

- changed axis.text.y = element_text(..., vjust = 1.0) to 0.5 (Benjamin
  Ostendorf)

[ENmix](https://bioconductor.org/packages/ENmix)
-----

                       Changes in version 1.25.9                        

- bugfix, updated user manual

[EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser)
-----------------

                       Changes in version 2.20.0                        

- New function `import` to import results from differential expression
  analysis with limma, edgeR, and DESeq2

- New function `showAvailableSpecies` to list supported species
  for a gene set database of choice (GO, KEGG, MSigDB, Enrichr, ...)

- New function `showAvailableCollections` to list provided
  gene set collections for a supported species of a gene set database
  of
  choice (GO, KEGG, MSigDB, Enrichr, ...)

- Gene sets: obtaining and caching of gene sets for different gene ID
  types (new argument `gene.id.types` for function `getGenesets`)

- GO gene sets: filter by GO evidence codes (new argument `evid` for
  function `getGenesets`)

- Including NEAT among nbea methods

[enrichplot](https://bioconductor.org/packages/enrichplot)
----------

                        Changes in version 1.9.5                        

- fix wordcloud_i (2020-10-15, Thu)
- Remove similarity calculation from emapplot

                        Changes in version 1.9.4                        

- implement pairwise_termsim to calculate similarity of enriched
terms
(2020-10-09, Fri)
- https://github.com/YuLab-SMU/enrichplot/pull/67
- change parameters to be more consistent
- https://github.com/YuLab-SMU/enrichplot/pull/62

                        Changes in version 1.9.3                        

- add node_label_size parameter to adjust the size of node label in
emapplot function (2020-09-18, Fri)

                        Changes in version 1.9.2                        

- add function emapplot_cluster (2020-09-01, Tue)

[ensemblVEP](https://bioconductor.org/packages/ensemblVEP)
----------

                       Changes in version 1.32.0                        

- add support for Ensembl release 100/101

[epivizrServer](https://bioconductor.org/packages/epivizrServer)
-------------

                       Changes in version 999.999                       

- This NEWS file is only a placeholder. The version 999.999 does
  not really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrServer>

[escape](https://bioconductor.org/packages/escape)
------

                       Changes in version 0.99.9                        

- Changing Seurat dependency, updated vignette

                       Changes in version 0.99.8                        

- Edited getSignificance ANOVA model call

                       Changes in version 0.99.7                        

- Edited getSignificance fit call to match documentation

                       Changes in version 0.99.6                        

- Edited match.args() in getSignificance

                       Changes in version 0.99.5                        

- Edited match.args() in getSignificance

                       Changes in version 0.99.4                        

- Added match.args() to getSignificance

- Changed stop() to message()

- Modified getSignficance to allow for ANOVA and T.test

                       Changes in version 0.99.3                        

- Updated link in description of getGeneSets.

                       Changes in version 0.99.2                        

- *Fixed a parenthesis, yeah a parenthesis. (In enrichIt() call I
  edited for 99.1)

                       Changes in version 0.99.1                        

- Removed parallel call in gsva() and added biocparallel

- Changed cores = 4 to cores = 2 in the vignette

                       Changes in version 0.99.0                        

- Preparing for bioconductor submission

[ExperimentHub](https://bioconductor.org/packages/ExperimentHub)
-------------

                       Changes in version 1.15.0                        

BUG FIXES

- (1.15.2) Proxy message formatting

USER-VISIBLE CHANGES

- (1.15.4) less stringent internet check

- (1.15.3) Add link for github issue reporting

- (1.15.1) Add hubs@bioconductor.org email for help

[ExperimentSubset](https://bioconductor.org/packages/ExperimentSubset)
----------------

                 Changes in version 0.99.0 (2020-10-02)                 

PRE-RELEASE

- Pre-release version of the ExperimentSubset package to be submitted
to Bioconductor

[ExploreModelMatrix](https://bioconductor.org/packages/ExploreModelMatrix)
------------------

                        Changes in version 1.1.4                        

- Updated interface to make better use of screen space
- Expanded explanation of statistical concepts

[FamAgg](https://bioconductor.org/packages/FamAgg)
------

                       Changes in version 1.17.1                        

- Fix issue #20 (binary trait drops names upon conversion).

[famat](https://bioconductor.org/packages/famat)
-----

                 Changes in version 0.99.9 (2020-10-26)                 

- Updated examples data

                 Changes in version 0.99.8 (2020-10-24)                 

- Fixed unit tests

                 Changes in version 0.99.7 (2020-10-19)                 

- Removed UniprotR package

                 Changes in version 0.99.6 (2020-10-19)                 

- Updated Biocmanager, devel version

                 Changes in version 0.99.5 (2020-10-19)                 

- Fixed unit tests

                 Changes in version 0.99.4 (2020-10-16)                 

- Removed biomart package

                 Changes in version 0.99.3 (2020-10-16)                 

- Updated NEWS file

                 Changes in version 0.99.2 (2020-10-16)                 

- Changed loops by apply functions
- Updated datasets documentation
- Updated DESCRIPTION file

                 Changes in version 0.99.1 (2020-10-01)                 

- Fixed testthat issues

                 Changes in version 0.99.0 (2020-09-24)                 

- Updated version for Bioconductor submission

               Changes in version 0.0.0.9000 (2020-09-01)               

- creation of the package

[FGNet](https://bioconductor.org/packages/FGNet)
-----

                        Changes in version 3.23                         

- Kegg-related functionalities have been removed.

[fgsea](https://bioconductor.org/packages/fgsea)
-----

                       Changes in version 1.15.2                        

- Faster perturbate thanks to Nikolay Budin

- Cleaner P-value and error estimations

[FilterFFPE](https://bioconductor.org/packages/FilterFFPE)
----------

                 Changes in version 0.99.3 (2020-10-08)                 

- Fix error in the algorithm under multiple threading

- Speed up the filtering process using minMapBase

                 Changes in version 0.99.0 (2020-08-20)                 

- Submitted to Bioconductor

[fishpond](https://bioconductor.org/packages/fishpond)
--------

                        Changes in version 1.6.0                        

- Added makeInfReps() to create pseudo-inferential replicates
  via negative binomial simulation from mean and variance
  matrices. Note: the mean and the variance provide the
  _inferential_ distribution per element of the count matrix.
  See preprint for details, doi: 10.1101/2020.07.06.189639.

- Added splitSwish() and addStatsFromCSV(), which can be used
  to distribute running of Swish across a number of jobs
  managed by `Snakemake`. See vignette for description of
  a suggested workflow. For a large single-cell dataset
  with mean and variance summaries of inferential uncertainty,
  splitSwish() avoids generating the inferential replicate
  counts until the data has been split into smaller pieces and
  sent to different jobs, then only the necessary summary
  statistics are gathered and q-values computed by
  addStatsFromCSV().

- plotInfReps() gains many new features to facilitate plotting of
  inferential count distributions for single cells, as quantified
  with alevin and imported with tximport. E.g. allow for numeric
  `x` argument plus grouping with `cov` for showing
  counts over pseudotime across groups of cells. Also added
  `applySF` argument which can be used to divide out a
  size factor, and the `reorder` argument which will re-order
  the samples/cells within groups by the count. plotInfReps()
  will draw boxplots with progressively thinner visual features
  as the number of cells grows to make the plots still legible.

                        Changes in version 1.5.2                        

- First version of makeInfReps(), to create pseudo-infReps
  via negative binomial simulation from set of mean and
  variance matrices in the assays of the SummarizedExperiment.

[flowCut](https://bioconductor.org/packages/flowCut)
-------

                        Changes in version 1.0.0                        

NEW FEATURES

- (v. 0.99.0) This is the submitted version of the package.
- (v. 0.99.6) Added many crashes to avoid crashes.
- (v. 0.99.7) Added vignette files.
- (v. 0.99.8) Added "Was the file run twice" in data table. Edited
data file to be consistent with fixes from other versions. Manual
pages updated.
- (v. 0.99.11) Added UseCairo parameter since docker requires no
Cairo
when creating the pngs.
- (v. 0.99.12) Added UnifTimeCheck parameter to allow the user to
adjust the value.
- (v. 0.99.13) Changed default UseCairo to F. Added UseCairo and
UnifTimeCheck parameter to second run code. Changed email.
- (v. 0.99.14) Added AlwaysClean parameter to allow for not skipping
if nice. Also added to the man pages.
- (v. 0.99.16) Added IgnoreMonotonic paramater to ignore the
monotonically increasing in time flagging test.
- (v. 0.99.17) Changed license.
- (v. 0.99.18) Added monotonic time fix option.

BUG FIXES

- (v 0.99.1) Small fixes from the bioconductor bot review process.
- (v 0.99.2) Review fixes.
- (v 0.99.9) Fixed a bug that occured if the first channel being
cleaned was the first one. This only happens if there is not FSC or
SSC channels.
- (v 0.99.10) Fixed bug where time tests fails on second run of
flowCut
- (v 0.99.15) Small bug fix with monotonically increasing channels.
- (v 0.99.19) Fixed formatting issues in vignette.
- (v 0.99.20) Fixed formatting issues in vignette.
- (v 0.99.21) Fixed spacing issues.
- (v 0.99.22) changed 1: to seq_len
- (v 0.99.23) changed codes formatting
- (v 0.99.24) changed codes formatting
- (v 0.99.25) More formatting and code optimization edits, changed
NEWS to NEWS.md
- (v 0.99.26) resolved installation issues
- (v 0.99.27) Small bug fix with monotonically increasing channels.

[flowFP](https://bioconductor.org/packages/flowFP)
------

                       Changes in version 1.47.0                        

- The update to R4.0.2 caused warnings like:
  Warning messages:
  1: In .Call("bin_level", fcs@exprs, model@.tmp_tags,
  model@split_axis[[level]], :
  converting NULL pointer to R NULL.
  The c code was updated to avoid returning a null pointer.

- Updated citation syntax was incorporated in the CITATION file.

[flowSpecs](https://bioconductor.org/packages/flowSpecs)
---------

                 Changes in version 1.3.3 (2020-09-03)                  

- Changes to two tests due to minor errors

- Excluding the recommendation to use flowVS, due to its deprecation.

                 Changes in version 1.3.2 (2020-05-15)                  

- Addition of peakNorm and associated test.

                 Changes in version 1.3.1 (2020-05-15)                  

- Bug fix in arcTrans.

- Correcting the way specUnmix exchanges exprs objects.

[fmrs](https://bioconductor.org/packages/fmrs)
----

                       Changes in version 0.99.3                        

IMPROVEMENTS SINCE LAST RELEASE

- Compatibility with bioconductor is added.

                       Changes in version 0.99.0                        

- Package moved to bioconductor.

BUG FIXES

- Several bugs are fixed.

[FRASER](https://bioconductor.org/packages/FRASER)
------

                        Changes in version 1.1.6                        

- Use proper S3/S4 methods to share functions between packages

- Minor API changes due to S3/S4 changes (e.g fds -> object)

- Switch from psiSite to theta

- Improved documentation

- Minor bugfixes

                        Changes in version 1.1.3                        

- Update and adjust injectOutlier and hyperParameter functions

- Option to compute z scores in logit space or not

- Add cap value [0.01,0.99] to logit function

- Use pairedEnd counting with Rsubread

- Correct assayName pajd -> padj

- Minor bugfixes

                        Changes in version 1.1.2                        

- Option to consider only the standard chromosomes in the counting

- Option to include additional columns from mcols(fds) in the result
  table

- Annotation of junctions with corresponding gene names/ids now
  produces
  an additional column in mcols(fds) that contains further gene
  names/ids
  if the junction overlaps with multiple genes

- Minor bugfixes

                        Changes in version 1.1.1                        

- Bugfix correcting the strand specific counting for paired-end reads

[FScanR](https://bioconductor.org/packages/FScanR)
------

                 Changes in version 0.99.7 (2020-09-08)                 

- ReSubmitted to Bioconductor

                 Changes in version 0.99.6 (2020-08-27)                 

- ReSubmitted to Bioconductor

                 Changes in version 0.99.5 (2020-08-27)                 

- ReSubmitted to Bioconductor

                 Changes in version 0.99.4 (2020-08-27)                 

- ReSubmitted to Bioconductor
  o Remove duplicates of detected PRF events from output

                 Changes in version 0.99.3 (2020-08-26)                 

- Submitted to Bioconductor

[gage](https://bioconductor.org/packages/gage)
----

                       Changes in version 2.39.3                        

- fixed warning on "library(GO.db)" in go.gsets.R. Now "import" from
  GO.db, instead of "suggest" it.

                       Changes in version 2.39.2                        

- fixed error caused by class(exprs) == "data.frame" in
  gagePrep.R. class(exprs) now returns a vector of length 2, which
  caused the error.

[GAPGOM](https://bioconductor.org/packages/GAPGOM)
------

                 Changes in version 1.5.1 (2020-07-01)                  

Changed

- Added cre to new maintainer in description file.

[GCSConnection](https://bioconductor.org/packages/GCSConnection)
-------------

                        Changes in version 1.1.6                        

- add gcs_rm function

                        Changes in version 1.1.5                        

- Add gcloud into the default authentication process

                        Changes in version 1.1.4                        

- gcs_is_requester_pays supports uri

- Better print format in gcs_get_cloud_auth

                        Changes in version 1.1.3                        

- FileClass object can show the file URL now

- no warning will be given if `gcs_dir` find a non-standard file path

- Fix some word issues: all xxx_url functions are renamed to xxx_uri

                        Changes in version 1.1.2                        

- Support Requester Pays

- Support `~` symbol to go to the bucket root

- Support conversion from Folder/File class to character

[gdsfmt](https://bioconductor.org/packages/gdsfmt)
------

                       Changes in version 1.24.1                        

UTILITIES

- 'show' option in `print.gds.class()` for array preview

[GENESIS](https://bioconductor.org/packages/GENESIS)
-------

                       Changes in version 2.19.7                        

- Change default value of small.samp.correct in pcrelate to TRUE.

- Add options to remove NxN matrices from a null model (function
  nullModelSmall and fitNullModel argument return.small).

- Add check for collinearity in covariates.

                       Changes in version 2.19.6                        

- Add test options "BinomiRare" and "CMP" to assocTestSingle and
  assocTestAggregate.

                       Changes in version 2.19.5                        

- Add function jointScoreTest to perform a joint score test of a
  set of variants using a null model and a matrix of genotype
  dosages.

                       Changes in version 2.19.4                        

- Add function effectAllele to return the effect allele for
  association tests.

                       Changes in version 2.19.1                        

- Force design matrices to be non-sparse.

[GeneTonic](https://bioconductor.org/packages/GeneTonic)
---------

                        Changes in version 1.2.0                        

New features

- The geneset distillery is officially open! GeneTonic offers
functionality to aggregate together gene sets into overarching
biological themes, based on a network-based refinement of the
enrichment map. Corresponding graphical functionalities are also
extended to accommodate meta-genesets. An efficient implementation
for the Markov clustering on graph objects is also provided

- GeneTonic can now receive the input of many other tools for
functional enrichment analysis - this includes the output (text
export) of DAVID (shake_davidResult), enrichr (from website and via
the package, with shake_enrichrResult), fgsea (shake_fgseaResult),
and g:Profiler (with shake_gprofilerResult, which can handle the
textual output from the website, as well the one from the call to
the gost in gprofiler2)

- An export button to a SummarizedExperiment object for iSEE and its
underlying machinery has been added. If the visualization options in
GeneTonic are not exactly what you would expect, you might find an
excellent venue in the iSEE framework

Other notes

- Added an additional mechanism for safe fails when not finding the
GO
Term and searching for the definition - this could happen e.g. when
the term becomes outdated and is removed from the GO.db package, or
also mistyped if entered by hand at some point.
- gs_heatmap has a new parameter, plot_title, to override the title
to
be displayed and set it to any custom string
- It is now possible to save a snapshot of the graphs created with
visNetwork
- The Gene Box now also contains links to the GTEx portal for the
selected feature
- export_to_sif enables to export a graph object to a text file,
encoded with the SIF format

[GenomeInfoDb](https://bioconductor.org/packages/GenomeInfoDb)
------------

                       Changes in version 1.26.0                        

NEW FEATURES

- The seqlevelsStyle() getter and setter now support style "RefSeq"
  when
  the underlying genome is known.

- Register a bunch of new NCBI assemblies and UCSC genomes. Use
  registered_NCBI_assemblies() and registered_UCSC_genomes() to get the
  lists of supported NCBI assemblies and UCSC genomes.

SIGNIFICANT USER-VISIBLE CHANGES

- The seqlevelsStyle() getter and setter do a better job when the
  underlying genome is known. The new behaviors address two
  long-standing
  shortcomings of the old behaviors:
  - In general, the seqlevelsStyle() setter didn't know how to rename
  the scaffolds in an object. Now it does.
  - Also, for some assemblies that use unconventional chromosome naming
  conventions (e.g. Macaca_fascicularis_5.0), the seqlevelsStyle()
  getter was not able to detect the naming style and the
  seqlevelsStyle()
  setter was not able to rename the chromosomes. Now they both do the
  right thing.
  These improvements address these shortcomings but only in the
  situation
  where the underlying genome is known e.g. when 'unique(genome(x))' is
  "macFas5" or "Macaca_fascicularis_5.0". When the underlying genome is
  not known, nothing has changed.

DEPRECATED AND DEFUNCT

- Deprecate releaseName() method for GenomeDescription objects.

BUG FIXES

- Small fix to getChromInfoFromNCBI().

[genomeIntervals](https://bioconductor.org/packages/genomeIntervals)
---------------

                       Changes in version 1.45.2                        

- Ported changes from 1.44.2

                       Changes in version 1.45.1                        

- Ported changes from 1.44.1

                       Changes in version 1.44.2                        

- Updated the maintainer email address

- Cleared the imports

                       Changes in version 1.44.1                        

- Updated the vignette (data.frame() default is not a factor anymore
  for character vectors in R4)

- Removed the defunct RangedData usage

[GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures)
---------------

                       Changes in version 1.42.0                        

NEW FEATURES

- Implement a restricted seqinfo() setter for TxDb objects that
  supports
  altering only the seqlevels and/or genome of the object, but not its
  seqlengths or circularity flags. This is all we need to make the
  improved
  seqlevelsStyle() setter work on TxDb objects (see below).

SIGNIFICANT USER-VISIBLE CHANGES

- The seqlevelsStyle() getter and setter do a better job when the
  underlying genome is known. See NEWS file in the GenomeInfoDb package
  for more information.

[GenomicOZone](https://bioconductor.org/packages/GenomicOZone)
------------

                        Changes in version 1.3.1                        

- Fixed a bug of always using Enssembl US server to annotate the genes.

- Updated the description in the DESCRIPTION file.

- Updated the citations.

[GenomicRanges](https://bioconductor.org/packages/GenomicRanges)
-------------

                       Changes in version 1.42.0                        

NEW FEATURES

- Add nearestKNeighbors() method for GenomicRanges derivatives.

- coverage() now supports 'method="naive"'. This is in addition to the
  already supported methods "sort" and "hash". This new method is a
  slower
  version of the "hash" method that has the advantage of avoiding
  floating
  point artefacts in the no-coverage regions of the numeric-Rle object
  returned by coverage() when the weights are supplied as a numeric
  vector
  of type 'double'. See "FLOATING POINT ARITHMETIC CAN BRING A
  SURPRISE"
  example in '?coverage' in the IRanges package.

[getDEE2](https://bioconductor.org/packages/getDEE2)
-------

                       Changes in version 0.99.16                       

- Functions `getDee2Metadata` and `queryDee2` are now called
  `getDee2Metadata` and `queryDEE2` respectively to be consistent
  with the other functions.

- Fixed a bug with some samples that have a # in the name.
  Thanks to @uilnauyis for the suggestion.

- New function getDEE2_bundle to fetch entire project data from
  http://dee2.io/huge/

                       Changes in version 0.99.9                        

- New function se() which constructs SummarizedExperiment object

                        Changes in version 0.0.2                        

- Some slight change to the vignette

[ggtree](https://bioconductor.org/packages/ggtree)
------

                        Changes in version 2.3.7                        

- add label_pad() function to add padding characters to taxa labels
(2020-10-09,
Fri)
-
https://groups.google.com/g/bioc-ggtree/c/INJ0Nfkq3b0/m/lXefnfV5AQAJ
- add family parameter to geom_tiplab()

                        Changes in version 2.3.6                        

- new layouts, roundrect and ellipse
- https://github.com/YuLab-SMU/ggtree/pull/344
- https://github.com/YuLab-SMU/ggtree/pull/346
- fortify() method for treedataList object (2020-09-20, Sun)
- vexpand() and ggexpand() to expand plot limit by ratio of plot
range
(2020-09-18, Fri)
- geom_cladelab(), an updated version of geom_cladelabel that
supports
aes mapping (2020-09-17, Thu)
- https://github.com/YuLab-SMU/ggtree/pull/342

                        Changes in version 2.3.5                        

- td_unnest() which return a function to flatten ggtree plot data
(2020-09-14, Mon)
- https://yulab-smu.top/treedata-book/chapter12.html#td_unnest
- update geom_hilight to support geom_hilight(data = mydata, node =
selected_node) (2020-09-03, Thu)
- Defunct geom_nodelab2() (2020-09-02, Wed)
- geom_tiplab() and geom_nodelab() support geom = "shadowtext"
- td_filter() which return a function to subset ggtree plot data in
geom layers (2020-08-29, Sat)
- https://yulab-smu.top/treedata-book/chapter12.html#td_filter
- update man files of geom_rootedge and geom_point2
- update geom_hilight to support geom_hilight(data = tbl_tree, node =
selected_node). (2020-09-03, Thu)

                        Changes in version 2.3.4                        

- zoomClade and geom_zoom_clade to zoom in selected clade
(2020-08-04,
Tue)
- these two functions are wrapper function of ggforce::facet_zoom
- update facet_labeller according to the change of ggplot2
(2020-07-28, Tue)
- defunct set_hilight_legend as now geom_hilight supports aesthetic
mapping and can generate legend automotically
- remove annotation_image, phylopic and subview as they were defunct
for quite a long time. Users should refer to the ggimage package if
they want to annotate tree with image or subplots.
- as_ylab parameter added in geom_tiplab(), which supports displaying
tip labels as y-axis label and only works for rectangular and
dendrogram layouts
- hexpand to expand x limits by ratio of x range and supports both
direction (1 for rhs and -1 for lhs) (2020-07-27, Mon)

                        Changes in version 2.3.3                        

- add type parameter in geom_hilight, default is auto, optional rect
to rectangular layer, encircle to encircle layer and comment
original geom_hilight, and support subset in aesthetic. (2020-07-23,
Thu)
- update geom_hilight to support aesthetic mapping (2020-07-22, Wed)
- update geom_taxalink to support aesthetic mapping (2020-07-20, Mon)
- layout_inward_circular for layout_circular() + scale_x_reverse()
(2020-07-16, Thu)

                        Changes in version 2.3.2                        

- update geom_taxalink to support circular layout tree (2020-07-13,
Mon).

                        Changes in version 2.3.1                        

- fortify method for pvclust object (2020-06-21, Mon)
- add dot parameters for color or size of geom_hilight and more
detail
messages of warnings for extendto. (2020-06-16, Tue)
- modified the angle of clade labels. Added horizontal parameter to
control whether set clade labels to horizontal. When the parameter
was set to FALSE, it will be useful for the layouts in coord_polar,
such as circular, fan, radial. To better view the clade labels,
their angles has been adjusted. (2020-06-15, Mon)
- bug fixed in getYcoord_scale_category (2020-05-13, Wed)

[ggtreeExtra](https://bioconductor.org/packages/ggtreeExtra)
-----------

                       Changes in version 0.99.0                        

the 0.99.0 or 0.99.x version mean I am submitting it to Bioconductor. (20200710, Fri)

- 0.99.1 change svg of dev to png, set the dpi to 300. (20200714,
Tue)
- 0.99.2 geom_axis_text support the single column axis. (20200717,
Fri)
- 0.99.3 support inward_circular tree and geom_axis_text was build by
position_identityx. (20200724, Fri)
- 0.99.4 fix the constant of barplot bug. (2020-07-25, Sat)
- 0.99.5 better support the layer when x axis is reverse.
(2020-07-29,
Wed)
- 0.99.6 support adjusting the angle of geom_text. (2020-07-31, Fri)
- 0.99.7 revise stylistic comment of R code and examples.
(2020-08-04,
Tue)
- 0.99.8 change the formatting of code chunk. (2020-08-04, Tue)

0.99.9

- modified the color aesthetics for geom_boxplot and geom_violin.
(2020-08-14, Fri)
- keep all column of data of tree when it merge with external data.
(2020-08-14, Fri)

0.99.10

- add geom_ringline to create the grid line of external ring layers.
(2020-08-21, Fri)
- the addbrink and other argument control the line of margin has been
removed. (2020-08-21, Fri)
- add pseudo axis line in geom_axis_text. (2020-08-21, Fri)

0.99.11

- geom_axis_text and geom_ringline are removed. (2020-08-24, Mon)
- user can use axis.params=list(add.axis=TRUE) and
grid.params=list(add.grid=TRUE) to add the axis and grid lines of
external layers, respectively. (2020-08-24, Mon)

0.99.12

- add upper and lower grid lines of y. (2020-08-26, Wed)

0.99.13

- update the method of ggplot_add.layer_fruits, when the offset
between different fruit_plot is different, will use each offset in
each fruit_plot. (2020-08-31, Mon)

0.99.14

- update normxy to support dendrogram layout. (2020-09-02, Wed)

0.99.15

- add add.another.axis in geom_fruit to add another axis.
(2020-09-04,
Fri)
- update normxy to fix the bug when negative values are present.
(2020-09-04, Fri)
- update ggplot_add method of geom_fruit to support the orientation
which x axis of the external is from bottom to top, when layout is
dendrogram. (2020-09-04, Fri)

0.99.16

- change add.axis in axis.params from TRUE or FALSE to x or y or xy.
(2020-09-05, Sat)
- remove add.grid in grid.params and default of grid.params is NULL.
(2020-09-07, Mon)
- update method of axis tick and remove nbreaks, the breaks will be
calculate by pretty. (2020-09-07, Mon)
- use substitute to allow list of axis.params or grid.params has
empty
argument, eg grid.params=list(color="black",). (2020-09-07, Mon)

0.99.17

- change add.axis to axis in axis.params. (2020-09-08, Tue)

0.99.18

- add nbreak in axis.params of geom_fruit, it will be sent to n of
pretty to generate the desired number of intervals of axis.
(2020-09-08, Tue)

0.99.19

- modified the namespace, remove geom_vline, add geom_segment.
(2020-09-11, Fri)

                        Changes in version 0.0.1                        

- add vignettes. (20200707, Tue)

                       Changes in version 0.0.0.9                       

- support add the axis text of extra layer of geom_fruit using
geom_axis_text. (20200630, Tue)
- optimize normalization of extra layer data. (20200701, Wed)
- fixed bug: changed layer to layers. (20200703, Fri)

                       Changes in version 0.0.0.8                       

- support data of NULL in geom_fruit, and user can add data by %<+%
of
ggtree. (20200628, Sun)

                       Changes in version 0.0.0.7                       

- add geom_fruit_list to support add the same position for multi
layers. (20200622, Mon)

                       Changes in version 0.0.0.6                       

- automatically detect the 'position'. (20200612)

                       Changes in version 0.0.0.5                       

- change pratio parameter to pwidth and remove tippoint parameter.
- change geom_add function to geom_totree function. (20200610)

                       Changes in version 0.0.0.4                       

- support tip point for geom_star or geom_point
- support geom_boxplot and geom_violin (20200606)

                       Changes in version 0.0.0.3                       

- add marginal line, and when x is character, the new x normalized
should be started with zero. (20200601)

                       Changes in version 0.0.0.2                       

- The distance between the panel and tree can be adjusted using the
"offset". The value of associate panel were normalized in the range
of x of tree. The width can be adjusted using the "pratio". The
"offset" and "pratio" are the ratio related to tree. (20200529)

                       Changes in version 0.0.0.1                       

- first version to github (20200528)

[Glimma](https://bioconductor.org/packages/Glimma)
------

                        Changes in version 2.0.0                        

- Brand new backend and API with major changes

- Plots can now be embedded in html reports

- Plots can now be saved

- Data from gene annotation and counts can now be saved

- Added many plot customisation options in MDS plot

[glmGamPoi](https://bioconductor.org/packages/glmGamPoi)
---------

                         Changes in version 1.1                         

- Remove dual likelihood functions for overdispersion estimation.
  Instead merge functionality into conventional_***. This should
  cause no user facing changes, however should make it easier to
  maintain the package

- Make conventional_score_function_fast() more robust to extreme
  inputs. Avoid numerically imprecise subtractions and employ
  bounds based on series expansions for very small input

- If dispersion estimate quits because there is no maximum or
  all y are 0, return iterations = 0

- Add limits (1e-16 / 1e16) for nlminb estimates of the
  dispersion. This protects against errors due to NA's in
  the conventional_likelihood_fast

- Automatically set 'size_factors = FALSE' for input with
  0 or 1 row. This will change the estimated beta, but not the
  mu's

- Rename gampoi_overdispersion_mle() -> overdispersion_mle()

- Store data in the object returned by glm_gp()

- Remove Y from the interface of residuals.glmGamPoi, because
  I can just get it directly from fit$data

- Add function test_de() that does a quasi-likelihood ratio
  test to detect differentially expressed genes

- Add functionality to make a pseudobulk test directly
  from test_de() by aggregating the data around one column

- In group-wise beta estimation, fall back to optimize()
  if the Newton method fails

- Change the default size factor estimation method from
  "poscounts" to "normed_sum" and provide an easy way to
  call scran::calculateSumFactors()

- New "global" mode for dispersion estimation

[GNET2](https://bioconductor.org/packages/GNET2)
-----

                        Changes in version 1.5.5                        

- Add function for similarity scores between given labels for samples
  and the clustering in each predicted module.

                        Changes in version 1.5.2                        

- extract_edges() now returns an adjacency matrix for interations and
  their scores.

[GOSemSim](https://bioconductor.org/packages/GOSemSim)
--------

                       Changes in version 2.15.2                        

- new site, https://yulab-smu.top/biomedical-knowledge-mining-book/
for documentation (2020-09-04, Fri)
- update vignette
- update data/gotbl

                       Changes in version 2.15.1                        

- bug fixed of IC method when input IDs contain invalid terms.
(2020-07-25, Sat)

[GPA](https://bioconductor.org/packages/GPA)
---

                        Changes in version 1.1.0                        

- Made the following significant changes

- Incorporate ShinyGPA developed by Emma Kortemeier as a part of the R
  package.

- Add fitAll() to fit GPA models for possible pairs of phenotypes.

- Add shinyGPA() to open the shiny app for the interactive pleiotropy
  visualization.

[graphite](https://bioconductor.org/packages/graphite)
--------

                 Changes in version 1.35.2 (2020-10-12)                 

- Updated all pathway data.

[GreyListChIP](https://bioconductor.org/packages/GreyListChIP)
------------

                       Changes in version 1.22.0                        

- Add black lists for C. Elegans ce11, mouse dm6.

- Update existing black lists to version 2, or version 3 for
  human GRCh37, GRCh38.

[GSEAmining](https://bioconductor.org/packages/GSEAmining)
----------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[GSgalgoR](https://bioconductor.org/packages/GSgalgoR)
--------

                 Changes in version 0.99.3 (2020-10-23)                 

- fix several issues acording to bioc revision

- add accessor functions to galgo.Obj object

                 Changes in version 0.99.0 (2020-09-01)                 

- Remove gpuR support

- fix problem with small population

- Bioconductor Submission

                 Changes in version 0.5.0 (2020-08-02)                  

- Rename to GSGalgoR

                 Changes in version 0.4.2 (2020-07-21)                  

- Pre-release before Bioconductor submission

[Gviz](https://bioconductor.org/packages/Gviz)
----

                       Changes in version 1.33.1                        

NEW FEATURES

- score in sashimi plots can be transformed using function specified
  either in
  sashimiTransformation or transfromation (latter applied also to
  coverage)

BUG FIXES

- none

[gwascat](https://bioconductor.org/packages/gwascat)
-------

                       Changes in version 2.21.7                        

USER VISIBLE CHANGES

- Oct 19 2020 -- peculiar content in CHR_ID and CHR_POS cause
  truncation.  Improved read_tsv call
  by setting col_types

                        Changes in version 2.21                         

USER VISIBLE CHANGES

- April 30 2020 -- use BiocFileCache to manage retrieval from EBI

- May 2 2020 -- ebicat_2020_04_30 is a sample of 50000 records from a
  full retrieval

- May 2 2020 -- many data() elements moved to inst/legacy, LazyData
  turned off

[GWASTools](https://bioconductor.org/packages/GWASTools)
---------

                       Changes in version 1.35.2                        

- Replace read.table and write.table with much faster
  data.table::fread and data.table::fwrite in functions that
  convert to and from text. The exception is
  createAffyIntensityFile and
  checkIntensityFile(affy.inten=TRUE), which still use read.table
  to preserve the behavior of removing lines commented with the
  "#" character.

[HCAMatrixBrowser](https://bioconductor.org/packages/HCAMatrixBrowser)
----------------

                        Changes in version 1.0.0                        

New features

- HCAMatrixBrowser finally on Bioconductor!
- HCAMatrixBrowser uses OpenAPI Specification version 2 and
rapiclient
to provide R API representations.
- loadHCAMatrix provides users with matrix data given a set of
'bundle_fqids'
- Filtering on the main API object is supported see HCAMatrix for
details
- Representations in all formats is supported this includes (.csv,
.mtx, and .loom)
- Caching implemented using BiocFileCache
- MTX format support provided by Martin @mtmorgan
- Support for LOOM files provided by LoomExperiment
- CSV format files given as a tibble list

Bug fixes and minor improvements

- Updated vignettes to include API changes
- Allow for singleton bundle_fqid queries for v0 endpoint
(@dvantwisk,
#2)

[HDF5Array](https://bioconductor.org/packages/HDF5Array)
---------

                       Changes in version 1.18.0                        

NEW FEATURES

- Add 'as.sparse' argument to h5mread(), HDF5Array(), HDF5ArraySeed(),
  writeHDF5Array(), saveHDF5SummarizedExperiment(), and
  HDF5RealizationSink().
  Even though it won't change how the data is stored in the HDF5 file
  (data will still be stored the usual dense way), the 'as.sparse'
  argument allows the user to control whether the HDF5 dataset should
  be considered sparse (and treated as such) or not. More precisely,
  when HDF5Array() is called with 'as.sparse=TRUE', the returned object
  will be considered sparse i.e. blocks in the object will be loaded as
  sparse objects during block processing. This should lead to less
  memory usage and hopefully overall better performance.

- Add is_sparse() setter for HDF5Array and HDF5ArraySeed objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Change default value of 'verbose' argument from FALSE to NA for
  writeHDF5Array(), saveHDF5SummarizedExperiment(), and
  writeTENxMatrix().

BUG FIXES

- Fix handling of logical NAs in h5mread().

- Fix bug in saveHDF5SummarizedExperiment() when 'chunkdim' is
  specified.

[Herper](https://bioconductor.org/packages/Herper)
------

                 Changes in version 0.99.2 (2020-10-19)                 

- Updates for Bioconductor review

                 Changes in version 0.99.0 (2020-09-17)                 

- Submitted to Bioconductor

[HIBAG](https://bioconductor.org/packages/HIBAG)
-----

                       Changes in version 1.26.0                        

- users can interrupt the model building in an interactive R session

- remove `hlaErrMsg()` since it is never used

- a new option 'nthread' in `hlaAttrBagging()` as a complement to
  `hlaParallelAttrBagging()`

- kernel version 1.5: generates the same training model as v1.4,
  but 2-6x faster, by taking advantage of Intel AVX, AVX2 and AVX512
  intrinsics

- new function `hlaSetKernelTarget()` to automatically select the CPU
  target the algorithm is optimized for

[HiCcompare](https://bioconductor.org/packages/HiCcompare)
----------

                 Changes in version 1.11.1 (2020-08-14)                 

- Add CITATION, also to README.md file

- Update NEWS

- Update DESCRIPTION
  o Add Mikhail Dozmorov as Maintainer
  o Add URL and BugReports fields

- Revert breaking changes
  o KRnormalization.R, line 75, Z = rk/v; rho_km1 = t(rk) %*% Z;
  o hic_compare.R, line 99, A.min <- ceiling(mean(A_q10))

                 Changes in version 1.11.0 (2020-02-11)                 

- New method to read in .cool files, cooler2bedpe()

- New method for comparison, hic_compare()

[HPAStainR](https://bioconductor.org/packages/HPAStainR)
---------

                 Changes in version 0.99.2 (2020-09-02)                 

- R code

- Changed high stringency to include Enhanced and Supported removing
  Approved

                 Changes in version 0.99.1 (2020-08-06)                 

- Responded to comments from Bioconductor reviewer

- DESCRIPTION

- Moved `tibble` and `shiny` to imports

- Vignettes

- Made installation section

- Added table of contents

- Updated vignette text to agree with code

- Included sessionInfo()

- R code

- Updated `HPA_data_downloader` to fit bioconductor standards
  for accesing a website.

- Reduced number of lines >80

- Reduced number of indents not multiples of 4 spaces

                 Changes in version 0.99.0 (2020-07-15)                 

- Package prepared for submission to Bioconductor

- Updated to coding practices

- No Biocheck Errors

- Updated imports and examples

[hypeR](https://bioconductor.org/packages/hypeR)
-----

                       Changes in version 1.05.03                       

- rctbl_build() now wraps hyp objects into unlabled multihyp objects
- rctbl_build() nested tables shows the number of enriched genesets

                       Changes in version 1.05.02                       

- Correct file extensions (.rmd/.html) from output of hyp_to_rmd()
- Relative paths are now supported by hyp_to_rmd()

                       Changes in version 1.05.01                       

- Fixed bug for wrong column names from enrichr_available()
- Added the first shiny module for geneset selection

                       Changes in version 1.05.00                       

- Version bump for bioconductor
- Fixed hyp_dots(merge=TRUE) bug where some genesets were not showing
- Added support for fetching non-human Enrichr libraries (e.g. Yeast,
Fly, Worm, Fish)
- Better reporting through rctbl_hyp() and rctbl_mhyp()

[ideal](https://bioconductor.org/packages/ideal)
-----

                       Changes in version 1.14.0                        

Other notes

- Replaced dependency from d3heatmap with the functionality of
heatmaply

[idr2d](https://bioconductor.org/packages/idr2d)
-----

                        Changes in version 1.2.2                        

- skipped Python-dependent tests on i386 systems (64-bit version of
  Python is required)

[IgGeneUsage](https://bioconductor.org/packages/IgGeneUsage)
-----------

                 Changes in version 1.3.3 (2020-10-15)                  

- More flexible model for microbial sample analysis provided in
  inst folder, however, no interface is provided to it.

[illuminaio](https://bioconductor.org/packages/illuminaio)
----------

                 Changes in version 0.31.1 (2020-07-15)                 

UPDATES

- Error messages produced by readIDAT() when failing decrypt no longer
  appends
  an extra newline at the end.

- readIDAT() uses explicit stringsAsFactors=FALSE internally.

- readIDAT() no longer keeps two file connections open at the same.

DOCUMENTATION

- The BibTeX URL reported by citation(package="illuminaio") was broken.

BUG FIXES

- readBGX() would leave an open connection if there was an file reading
  error.

                 Changes in version 0.31.0 (2020-04-27)                 

NOTES

- The version number was bumped for the Bioconductor develop version,
  which is
  now BioC 3.12 for R (>= 4.0.0).

[ILoReg](https://bioconductor.org/packages/ILoReg)
------

                 Changes in version 0.99.0 (2020-06-25)                 

- The first version 0.99.0 is submitted to Bioconductor

[immunoClust](https://bioconductor.org/packages/immunoClust)
-----------

                       Changes in version 1.21.2                        

- CHANGES
  * plot of no clusters
  * soves the deprecated warning from flowCore
  * code cleaning

                       Changes in version 1.21.1                        

- CHANGES
  * authors email contact

[infercnv](https://bioconductor.org/packages/infercnv)
--------

                 Changes in version 1.5.3 (2020-10-23)                  

- Add check that detectCores() doesn't return NA before comparing. In
  case it does, just use the value provided as an option directly.

                 Changes in version 1.5.2 (2020-10-23)                  

- Added reordering of cells in metadata exported to a seurat object so
  that it always matches, in case the cells are not sorted in the same
  order in the data provided to infercnv and in seurat.

                 Changes in version 1.5.1 (2020-10-06)                  

- Fix to reload in cases where comparing NULL/NA.

- Fixed issue in denoising when trying to reload results from step 18
  or 19.

- Fixed MCMC Diagnostic plots by adding diagnostic generation.

- Update included data objects to contain additional option slot, and
  prevent common name collisions.

- Fully rename data objects and name of the vars they provide.

- Fix reference plotting not having access to the actual subclustering
  information but that of the previous provided data object (that was
  renamed to avoid name collision by mistake like this one).

- Added checks in add_to_seurat methods that there are gains/losses
  found when taking the top hits.

- Added check that output to write in add_to_seurat top regions is not
  null and output an empty file without erroring if it is.

- Change to remove genes in "chr_exclude" from counts before doing the
  read level filtering of cells.

- Added minimum read count requirement per cell of 1 after removal of
  "chr_exclude" genes so that there are no divisions by 0 when
  normalizing.

- Changed default min_max_counts_per_cell to select cells with at least
  100 counts by default.

- Fix to plotting for HMM coloring of heatmap when the full range of
  values are not present.

- Fix to plotting when no reference groups are used to not produce
  warnings.

[Informeasure](https://bioconductor.org/packages/Informeasure)
------------

                 Changes in version 0.99.8 (2020-10-06)                 

- make changes

                 Changes in version 0.99.7 (2020-10-06)                 

- make changes

                 Changes in version 0.99.6 (2020-09-21)                 

- make changes

                 Changes in version 0.99.5 (2020-09-20)                 

- Fix the the issue reported in the build report where used "1:lenght""
  instead of "seq_len"

- Use "<-" for assignment rather than "="

- Fix indentation

                 Changes in version 0.99.4 (2020-08-14)                 

- Fix the warning of unnecessary package calls

                 Changes in version 0.99.3 (2020-08-14)                 

- Fix the warning of unnecessary package calls

                 Changes in version 0.99.2 (2020-08-14)                 

- Fix the warning of unnecessary package calls

                 Changes in version 0.99.1 (2020-08-13)                 

- Made the following significant changes

- Added Informeasure.Rproj into the .gitignore file

- Replaced Author/Maintainer with Authors@R in the DESCRIPTION file

- Added a NEWS file

- Added a .R file in tests/ directory

- Update R version dependency from 3.5.0 to 4.0

- Use TRUE/FALSE instead of T/F in PMI.plugin()

                 Changes in version 0.99.0 (2020-08-09)                 

- Submitted to Bioconductor

[intansv](https://bioconductor.org/packages/intansv)
-------

                       Changes in version 1.29.0                        

Notes

- Fix build error on Windows and Mac.

[IRanges](https://bioconductor.org/packages/IRanges)
-------

                       Changes in version 2.24.0                        

NEW FEATURES

- coverage() now supports 'method="naive"'. This is in addition to the
  already supported methods "sort" and "hash". This new method is a
  slower
  version of the "hash" method that has the advantage of avoiding
  floating
  point artefacts in the no-coverage regions of the numeric-Rle object
  returned by coverage() when the weights are supplied as a numeric
  vector
  of type 'double'. See "FLOATING POINT ARITHMETIC CAN BRING A
  SURPRISE"
  example in '?coverage'.

DEPRECATED AND DEFUNCT

- Removed RangedData class and anything related to RangedData objects.

BUG FIXES

- Fix bug in list element recycling.

[ISAnalytics](https://bioconductor.org/packages/ISAnalytics)
-----------

                Changes in version 0.99.15 (2020-10-22)                 

- Minor fix in import_association_file file function: added multiple
strings to be translated as NA

                Changes in version 0.99.14 (2020-10-21)                 

- Minor fixes in tests

                Changes in version 0.99.13 (2020-10-19)                 

NEW FEATURES

- Added analysis functions CIS_grubbs and cumulative_count_union
- Added plotting functions CIS_volcano_plot

                Changes in version 0.99.12 (2020-10-04)                 

NEW FEATURES

- Added analysis function sample_statistics

SIGNIFICANT USER-VISIBLE CHANGES

- aggregate_values_by_key has a simplified interface and supports
multi-quantification matrices

MINOR CHANGES

- Updated vignettes
- import_parallel_Vispa2Matrices_interactive and
import_parallel_Vispa2Matrices_auto now have an option to return a
multi-quantification matrix directly after import instead of a list

                Changes in version 0.99.11 (2020-09-21)                 

NEW FEATURES

- Added analysis functions threshold_filter, top_integrations
- Added support for multi-quantification matrices in
compute_abundance

MINOR FIXES

- Fixed bug in comparison_matrix that ignored custom column names
- Fixed issues in some documentation pages

                Changes in version 0.99.10 (2020-09-14)                 

ISanalytics is officially on bioconductor!

NEW FEATURES

- Added analysis functions comparison_matrix and
separate_quant_matrices
- Added utility function as_sparse_matrix
- Added package logo

SIGNIFICANT USER-VISIBLE CHANGES

- Changed algorithm for compute_near_integrations
- Added support for multi-quantification matrices to
remove_collisions
- Added usage of lifecycle badges in documentation: users can now see
if a feature is experimental/maturing/stable etc

MINOR FIXES

- Added fix for import_single_Vispa2Matrix to remove non significant
0
values

                 Changes in version 0.99.9 (2020-09-01)                 

NEW FEATURES

- Added functionality: aggregate functions
- Added vignette on aggregate functions
- Added recalibration functions
- Added first analysis function (compute_abundance)

SIGNIFICANT USER-VISIBLE CHANGES

- Dropped structure ISADataFrame: now the package only uses standard
tibbles
- Modified package documentation

                 Changes in version 0.99.8 (2020-08-12)                 

- Submitted to Bioconductor

[iSEE](https://bioconductor.org/packages/iSEE)
----

                       Changes in version 2.1.27                        

- Minor edits to the API documentation.

                       Changes in version 2.1.26                        

- Disable the import button on the dimnames modal when transmitter
has
not made a selection.
- Separate the maximum number of factor levels for colors from other
applications.

                       Changes in version 2.1.25                        

- Support a named vector in the SearchColumns field of the Table
subclasses.

                       Changes in version 2.1.24                        

- Enable custom saving of the application state via the new
saveState=
argument.
- Switch colormap getters to use an internal cache to avoid conflicts
with user entries.

                       Changes in version 2.1.23                        

- ExperimentColorMap inherits from Annotated.

                       Changes in version 2.1.22                        

- Split and rename scripts for test setup.

                       Changes in version 2.1.21                        

- Minor fix to unit test.

                       Changes in version 2.1.20                        

- Export even more internal utilities for re-use in downstream
packages.
- Minor fixes to the ComplexHeatmapPlot documentation.

                       Changes in version 2.1.19                        

- Export more internal utilities for re-use in downstream packages.
- Added a generic to define the selection effect UI.
- Fixes to the ComplexHeatmapPlot observers, most obviously for the
assay choice.
- Fixes to the underlying reactive framework to avoid bugs due to
unresponsiveness.

                       Changes in version 2.1.18                        

- Bugfix to properly support dynamic classes on landing page.

                       Changes in version 2.1.17                        

- .refineParameters() for FeatureAssayPlot, SampleAssayPlot protects
the x-axis choice against absent metadata.
- CSS classes for each panel are now defined at app run-time, to make
it easier to write landing pages without specifying initial= in
iSEE().

                       Changes in version 2.1.16                        

- Bugfix for initialization of the ColorByFeatureDynamicSource UI
element in ColumnDotPlots.

                       Changes in version 2.1.15                        

- Allow labelling of the medoid for each level of a discrete variable
in a scatter-type DotPlot.
- Turned on validity checks during [[<- assignment into Panel
classes.

                       Changes in version 2.1.14                        

- Right-aligned the help icon for individual panel tour.
- Added tooltip for mouseovers on DotPlot panels.
- Allowed custom annotation about selected table row to to be
displayed in Table panels.

                       Changes in version 2.1.13                        

- Refactored the heatmap feature selection modal to be reusable for
selecting rows or columns in other contexts.

                       Changes in version 2.1.12                        

- Added panel-specific tours via the .definePanelTour() generic.
- Generalized the HiddenColumns mechanism to all Table subclasses.

                       Changes in version 2.1.11                        

- Added HiddenColumns slot to hide columns in ColumnDataTables and
RowDataTables.

                       Changes in version 2.1.10                        

- Avoided transmitting multiple selections for Tables when the number
of columns change.
- Hid irrelevant UI elements for Table multiple selections.
- Generalized information about the number of selected rows/columns
in
each panel.
- Streamlined the .defineOutput signature.

                        Changes in version 2.1.9                        

- Improved documentation for generics and classes.
- Defined .exportOutput method for the ComplexHeatmapPlot class.
- Added missing documentation for slots and methods in the
ComplexHeatmapPlot class.
- Enforced sensible defaults for the dynamic selection setting.
- Extract assays with dimnames for correct indexing.

                        Changes in version 2.1.8                        

- Fixed control of legend point size for continuous covariates.
- Extended control of legend point size for violin plots and Hinton
plots.

                        Changes in version 2.1.7                        

- Added control of legend point size under the "Text" category of teh
"Visual parameters" box.

                        Changes in version 2.1.6                        

- Fixed bug for sizeBy observers.

                        Changes in version 2.1.5                        

- Fixed initialization of panel size to current value when the
"Organize panels" window is closed and re-opened.
- Fixed removal of last panel from the interface.

                        Changes in version 2.1.4                        

- Added progress bar when exporting panel outputs.
- Fixed missing section in createLandingPage() man page.

                        Changes in version 2.1.3                        

- Fixed handling of logical > 1 when processing the CustomRowsText
slot in the ComplexHeatmapPlot constructor.

                        Changes in version 2.1.2                        

- Added a new vignette to describe panel links.
- Fixed documentation for *DynamicSource slots.
- Fixed reception of single selection from plot at initialization.
- Removed deprecated functionality.

                        Changes in version 2.1.1                        

- Added vignette documenting the use of out-of-memory matrices for
big
data.
- Added TENxPBMCData to Suggests:.

[iSEEu](https://bioconductor.org/packages/iSEEu)
-----

                        Changes in version 1.1.9                        

- Ensure that global parameters only affect panels during
construction.

                        Changes in version 1.1.8                        

- Added the AggregatedDotPlot panel to show marker-based dot plots.

                        Changes in version 1.1.7                        

- Improved safety and correctness of the calculation of the number of
DEGs.

                        Changes in version 1.1.6                        

- Version bump to trigger reinstallation with new iSEE class
definitions.

                        Changes in version 1.1.5                        

- Added panel-specific tours for all panel classes via the
.definePanelTour() generic.

                        Changes in version 1.1.4                        

- Generalized the DE-related globals to work as patterns rather
directly specifying the acceptable names.
- Align DynamicMarkerTable's treatment of getTableExtraFields() with
the globals strategy.

                        Changes in version 1.1.3                        

- Overhauled handling of global parameters for greater consistency.
- Added the LogFCLogFCPlot to plot two DE comparisons against each
other.
- Switched to KEGGREST to get the names of pathways.

                        Changes in version 1.1.2                        

- Replaced GeneSetTable with the more general FeatureSetTable.
Improved handling of arbitrary feature sets.
- Renamed DifferentialStatisticsTable to the more appropriate
DynamicMarkerTable. Support inclusion of extra fields from the
rowData.
- Global parameters now only affect construction of MAPlots and
VolcanoPlots.

                        Changes in version 1.1.1                        

- Improved documentation of the ReducedDimensionHexPlot methods.

[IsoCorrectoR](https://bioconductor.org/packages/IsoCorrectoR)
------------

                        Changes in version 1.6.2                        

- as response to a post from Am Zimmerman on the bioconductor support
  site (https://support.bioconductor.org/p/132345/), log-space
  factorials were used to calculate internal correction probabilities.

[IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR)
---------------------

                Changes in version 1.11.11 (2020-10-13)                 

- Update type: minor.

- Fixed the mistake in importIsoformExpression() introduced in last
  updated

                Changes in version 1.11.10 (2020-10-13)                 

- Update type: minor.

- Update of importIsoformExpression() fix the import of
  countsFromAbundance

                 Changes in version 1.11.9 (2020-10-12)                 

- Update type: Minor.

- More updates regarding namespace and dependencies

                 Changes in version 1.11.8 (2020-10-09)                 

- Update type: Minor.

- Description update regarding to namespace

                 Changes in version 1.11.7 (2020-09-30)                 

- Update type: Minor.

- Various documentation updates

- The plotting order of the sub-plots of switchPlot() was changed to
  avoid problems when having long isoform names.

- analyzeSignalP() was updated to be more robust at handling SignalP5
  data where very few predictions were done.

                 Changes in version 1.11.6 (2020-09-17)                 

- Update type: Minor.

- Updated namespace.

                 Changes in version 1.11.5 (2020-09-14)                 

- Update type: Minor.

- Updated example code in importIsoformExpression()

- Updated namespace

                 Changes in version 1.11.4 (2020-09-10)                 

- Update type: Medium.

- importRdata() was updated to give examples of sequence names when no
  overlap between fasta file and expression data was found.

- importRdata() was updated to try and rescue missing gene_name
  annoations (must likely due to novel transcripts) and split merged
  genes (a problem often occuring when doing transcript assembly with
  tools such as Cufflinks/StringTie).

- importRdata() and importCufflinksFiles() was udated with an option to
  print a guesstimate on the number of genes with differential isoform
  usage.

- isoformToGeneExp() and importGTF() was updated to look for the
  annotation problems fixed by importRdata() and give warnings if
  pressent.

- The example data (from individual files) was updated to include CDS.

- extractGeneExpression(), a function that extracts gene level
  counts/expression from a switchAnalyzeRlist was introduced.

- prepareSalmonFileDataFrame() and importSalmonData() was introduced.
  Jointly these functions enable import of Salmon data via tximeta
  thereby omitting the manual integration of annotation data
  (gtf/fasta)

- isoformSwitchAnalysisPart2() was updated to only do enrichment
  analysis if enough events were found.

- preFilter() was updated to apply the gene expression to both
  conditions instead of the average across all samples thereby better
  filtering out untrustworthy genes.

- extractSplicingGenomeWide() and extractConsequenceGenomeWide() was
  updated to handle missing values when calculating summary statistics.

- extractSplicingEnrichment(), extractSplicingEnrichmentComparison(),
  extractConsequenceEnrichment(),
  extractConsequenceEnrichmentComparison() was updated to use
  binom.test() instead of prop.test() as this test more apporpriate
  when analyzing smaller number of events.

- extractSplicingEnrichment() and extractSplicingEnrichmentComparison()
  was updated to have more easily interpretable descriptions.

- extractSwitchOverlap() was updated to also plot overlap in isoform
  switches and now allows for control of which venn diagrams to make.

- switchPlot() was updated to only consider the dIFcutoff when
  classifying the "increased/decrease/unchanged usage".

- switchPlotGeneExp(), switchPlotIsoExp(), switchPlotIsoUsage() and
  switchPlotTranscript() was updated to enable return of the ggplot2
  object (instead of printing it)

- all extractConsequence*() and extractSplicing*() functions was
  updated to enable return of the ggplot2 object (instead of printing
  it)

- switchPlotTopSwitches() was updated to also have the onlySigIsoforms
  argument.

- Various documentation updates.

                 Changes in version 1.11.3 (2020-05-20)                 

- Update type: Minor.

- importRdata was updated to handle GTF files with a lot of additional
  information.

- analyzeSignalP and analyzePFAM was updated to fix a problem with
  handling multiple files.

- analyzePFAM was updated to attempt to handle both fixed width files
  and broken fixed width files.

- isoformSwitchTestDEXSeq() and isoformSwitchTestDRIMSeq() was updated
  with the "keepIsoformInAllConditions" argument which allows data for
  an isoform to be kept in all comparisons even if it is only deemed
  significant in one comparison. TRUE by default.

- importCufflinksData() was updated to handle when
  pathToSplicingAnalysis was not used.

                 Changes in version 1.11.2 (2020-05-12)                 

- Update type: Minor.

- A bug was corrected in extractSequence() which caused an error:
  "object 'filterShortAALength' not found".

- In extractSequence() the minimul length kept when using the
  "removeShortAAseq" argument was raised to 11 amino acids to match the
  pfam website.

- An error in importRdata was fixed to re-enable removal of isoforms
  only found in annotation.

                 Changes in version 1.11.1 (2020-05-05)                 

- (Version bump due to Bioconductor release).

- Update type: Minor.

- Update of code comments in importRdata()

- Update of printed message in importIsoformExpression()

- Update of analyzePFAM() to enable more robust import of fwf files
  with and without headers included. It also handles the mistake in
  pfam files with regards to the fixed width of files when the
  "coiled-coil" type are included.

[isomiRs](https://bioconductor.org/packages/isomiRs)
-------

                       Changes in version 1.17.2                        

FIX

- Fix n() error from dplyr. Now importing in NAMESPACE.

                       Changes in version 1.17.1                        

FIX

- Fix error in tibble when selecting a column. Now it needs unlist to
  coarce to a vector.

[KEGGprofile](https://bioconductor.org/packages/KEGGprofile)
-----------

                       Changes in version 1.31.1                        

- Fix bugs in genes_kept to deal with NA in data

[LACE](https://bioconductor.org/packages/LACE)
----

                        Changes in version 1.0.1                        

- From Travis-CI to Github Actions.

[lefser](https://bioconductor.org/packages/lefser)
------

                        Changes in version 1.0.0                        

- lefser is the R implementation to the LEfSE method for microbiome
marker discovery.
- It uses the Kruskal-Wallis test, Wilcoxon-Rank Sum test, and Linear
Discriminant Analysis to find biomarkers in groups and sub-group
blocks.
- See the method paper here: https://doi.org/10.1186/gb-2011-12-6-r60

[limma](https://bioconductor.org/packages/limma)
-----

                       Changes in version 3.46.0                        

New functionality

- 
  Add new function chooseLowessSpan().

- 
  fitFDist() with a non-NULL `covariate` now fits a smoother
  trend (with fewer spline knots) than before unless there are at
  least 30 observations. This also affects squeezeVar() and
  eBayes() with `trend=TRUE`.

- 
  Add new argument `output.style` to weightedLowess().

- 
  write.fit() now outputs a blank column heading for the row
  names so that the number of column names in the delimited
  output file will agree with the number of columns (similar to
  write.csv in the base package).

Code improvements

- 
  Subsetting and contrasts.fit() now work on MArrayLM objects
  even if the cov.coefficients component is absent.

- 
  volcanoplot() now prompts user to run eBayes() if the object
  doesn't contain p-values.

- 
  Add checks and more informative error messages to getEAWP() and
  lmFit() when the data object is missing, NULL or has zero rows.

- 
  Subsetting for RGList, MAList, EListRaw, EList, MArrayLM or
  TestResults objects now accepts arguments other than `i` or `j`
  but ignores them without an error. This is relevant if a user
  adds a `drop` argument by analogy with matrix subsetting.
  Previously that produced an error.

- 
  Subsetting of TestResults now requires two arguments, same as
  all of the other data object classes listed above. Previously
  single index subsetting (for example object[i]) was allowed.

- 
  Improved treatment of zero weights by weighted.median().

- 
  More consistent use of rep.int(), rep_len() and rep()
  throughout the package.

- 
  Replace NA with NA_real_ where appropriate.

- 
  Add a message to topTableF() to warn that the function is
  obsolete and will be removed in a future version of limma.
  Remove usages of topTableF() from the User's Guide.

- 
  Remove obsolete function toptable(), which has been officially
  deprecated in favor of topTable() since 1 Feb 2018.

Documentation

- 
  Revise and expand the weightedLowess help page.

- 
  Add more explanation about "heirarchical" method in the details
  section of the decideTests() help page.

- 
  Add Note to voom help page to emphasise that voom is not an
  appropriate tool to analyse expresssion measures that have been
  normalized for library size.

- 
  Revise the help page for plotMDS() to clarify that the function
  can produce PCA plots as well as PCoA plots.

- 
  Add example to EList-class help page.

- 
  Edits to the eBayes help page.

- 
  Edit concept entries for help files (for reference manual
  index).

Bug fixes

- 
  Revise write.fit() so that column names of the output file are
  correct when `digits=NULL` and the `fit` object contains only
  one coefficient column. Previously the same contrast name was
  repeated as the column name for the coefficients, t-statistics
  and p.values, making them hard to distinguish.

- 
  Bug fix to arrayWeights() when `y` contains NAs, the design
  matrix has several columns and some but not all genes have no
  residual df.

[lipidr](https://bioconductor.org/packages/lipidr)
------

                         Changes in version 2.4                         

- Added imputation function for untargeted lipidomics datasets

- Added function to boxplot enriched lipidsets related to chain length
  and unsaturation

[LRBaseDbi](https://bioconductor.org/packages/LRBaseDbi)
---------

                        Changes in version 2.0.0                        

- pkgtitle and pkgdescription options were added in makeLRBasePackages

[Maaslin2](https://bioconductor.org/packages/Maaslin2)
--------

                        Changes in version 1.3.2                        

- Resolve issue with x labels missing from boxplots for metadata
  variables without levels and increase max jpgs written.

- Update the check for variables with more than two levels that will
  require a reference provided by the user for the model and the
  boxplots to be more strict (ignore UNK in the level count, don't
  check if it is a factor, and check to see if all of the levels are
  numeric to ignore continuous variables).

- Fix ZINB error.

                        Changes in version 1.3.1                        

- Add random effects capability to all the other non-LM models.

- Add variance filtering option with default set to 0.

- Normalization is now performed after filtering and N.not.zero is
  calculated on the untransformed space. Also includes minor edits to
  synchronize with the manuscript.

- Add reference option required for fixed effects variables with more
  than two levels. Reference used in model and for boxplots.

- Update heatmap to include all categorical levels.

[maftools](https://bioconductor.org/packages/maftools)
--------

                        Changes in version 3.12                         

BUG REPORTS

- Incorrect deduplication cases while validating MAFs. Issue: 623
- Fix repeated domain labels in lollipopPlot2. Issue: 614
- Fix Copy number labels and coloring in coBarPlot. Issue: 609
- Fix coloring issues for annotations in oncoplot in R < 4.0. Issue:
599

ENHANCEMENTS

- Updated TCGA TMB. All variants from MC3 results are restricted and
harmonized to Agilent 35.8 MB capture kit. See ?tcgaCompare for
details. Issue: 612
- Added SIMC1 protein to domain database. Issue: 616
- Added sampleOrder1 and sampleOrder2 arguments to coOncoplot. Issue:
592
- Added gene_mar outer_mar argument to coOncoplot. Issue: 260
- Added sample size warning messages to mafCompare. Issue: 602
- Added cohortFontSize and axisFontSize to tcgaCompare

NEW FUNCTIONS

- Added filterMaf function.

DEPRECATED

- signatureEnrichment will be deprecated in future. Currently kept
for
legacy purpose with a warning message. Issue: 607 615

[MAGeCKFlute](https://bioconductor.org/packages/MAGeCKFlute)
-----------

                        Changes in version 1.9.1                        

- Turn off the Depmap analysis in vignettes because it's time consuming

- Use count in EnrichedView when NES is not available

                        Changes in version 1.8.1                        

- Update the standard gene names

- Debug the color issue in the ScatterView

- Remove the usage of scale_color_npg and data.table::melt from the
  DensitiView

- Remove nPerm parameter from the enrich.GSE

- Download GO-gene relationship from NCBI ftp folder

- Move multiple imports to suggests

- Replace view_allpath with pathview.top in FluteMLE to allow users
  to set the number of pathways for pathview visualization.

[marr](https://bioconductor.org/packages/marr)
----

                       Changes in version 0.99.4                        

- code indentation updated.
- Renamed the class and constructor from marr to Marr on
October 23, 2020.

                       Changes in version 0.99.3                        

- Bioconductor single package building error

                       Changes in version 0.99.2                        

- Expanded Description and changed the binwidth parameter in
histograms on October 16, 2020.

                       Changes in version 0.99.1                        

- Added Authors@R to the description file.

                       Changes in version 0.99.0                        

- Submitted marr 0.99.0 to Bioconductor on October 2, 2020.

[MatrixGenerics](https://bioconductor.org/packages/MatrixGenerics)
--------------

                        Changes in version 1.2.0                        

- Add drop and type to generic signature of [row|col]Quantiles
  (<URL:
  https://github.com/Bioconductor/MatrixGenerics/pull/14>).

- Sync API with matrixStats v0.57.0 (<URL:
  https://github.com/Bioconductor/MatrixGenerics/issues/17>).

- Add default methods with user-friendly fallback mechanism
  (<URL:
  https://github.com/Bioconductor/MatrixGenerics/pull/16>).
  Suggested packages are now loaded the first time a
  MatrixGenerics' generic is called (e.g. the first time
  MatrixGenerics::colVars() is called). With this new approach,
  if the user passes a _dgCMatrix_ object and if
  sparseMatrixStats is already loaded, will 'just work' and the
  fallback mechanism won't try to load anything.

- Dispatch on methods for matrix objects when table objects are
  supplied (<URL:
  https://github.com/Bioconductor/MatrixGenerics/pull/15>)

[matter](https://bioconductor.org/packages/matter)
------

                       Changes in version 1.15.1                        

BUG FIXES

- Remove '...' from 'which()'

[megadepth](https://bioconductor.org/packages/megadepth)
---------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- Documentation website is now available at
http://LieberInstitute.github.io/megadepth/. It gets updated with
every commit on the master branch (bioc-devel) using GitHub Actions
and pkgdown.

[MesKit](https://bioconductor.org/packages/MesKit)
------

                       Changes in version 0.99.11                       

Changes

- Preparations for Bioconductor release 4.0.

- Bugfixes:
  + Fix 'plotCNA' to not show 'chrSilent'.

- Other changes:
  + Update exmaples.
  + Removed useless datasets.

                 Changes in version 0.99.0 (2019-06-01)                 

Changes

- Package created.

[metabCombiner](https://bioconductor.org/packages/metabCombiner)
-------------

                 Changes in version 0.99.0 (2020-09-11)                 

- Submitted to Bioconductor

[metabolomicsWorkbenchR](https://bioconductor.org/packages/metabolomicsWorkbenchR)
----------------------

                       Changes in version 0.99.0                        

- Bioconductor submission

[metaMS](https://bioconductor.org/packages/metaMS)
------

                       Changes in version 1.25.1                        

- Unit Test bug fix

                       Changes in version 1.25.0                        

- NEW function SearchNIST added (Only Windows OS)

[MetaNeighbor](https://bioconductor.org/packages/MetaNeighbor)
------------

                        Changes in version 1.9.1                        

- Add several parameters to control fast version of MetaNeighbor and
  MetaNeighborUS

- Add function to pre-train MetaNeighbor models.

- Add preprocessing function to merge SingleCellExperiment objects.

- Add utility functions related to meta-clusters, cluster graphs and
  visualization.

[metaseqR2](https://bioconductor.org/packages/metaseqR2)
---------

                 Changes in version 1.1.21 (2020-08-05)                 

NEW FEATURES

- Creation of UCSC track and assembly hubs through the provision of a
  FASTA
  genome and/or assembly file for a non-supported or custom organism.

BUG FIXES

- Fixed a problem with rownames when importing an existing counts
  table.

                 Changes in version 1.1.19 (2020-07-29)                 

NEW FEATURES

- None

BUG FIXES

- Fixed several problems with respect to custom GTF annotation import,
  especially when the GTF is compliant with minor GTF standards.

- Fixed a few warnings resulting in errors in R>4.0.

                 Changes in version 1.1.18 (2020-07-03)                 

NEW FEATURES

- Simple quantification, QC and reporting is now allowed without having
  to define a contrast or statistical tests to be performed.

- In case of only a contrast is defined but no statistical tests
  required
  only fold changes are calculated.

BUG FIXES

- Fixed undefined class/function leftovers from NOISeq integration.

                 Changes in version 1.1.17 (2020-06-29)                 

NEW FEATURES

- Integrated functionalities from deprecated packages DESeq and
  NOISeq in order to maintain PANDORA

BUG FIXES

- Fixed report crash with the new readTargets wrapper function.

                 Changes in version 1.1.16 (2020-06-22)                 

NEW FEATURES

- Added support for reading targets in YAML

BUG FIXES

- Fixed "rnacomp" plot trying to import to JSON db while not run.

- Fixed report crash when "exportWhere" was NA.

- Other minor bug fixes

                 Changes in version 1.1.15 (2020-06-12)                 

NEW FEATURES

- Allow simple differential exon analysis

BUG FIXES

- Fixed annotation breaks from UCSC (rn, danRer).

- Fixed crash when using excludeList (due to attributes being gone
  after
  list subsetting).

                 Changes in version 1.1.13 (2020-05-25)                 

NEW FEATURES

- Addition of the harmonic mean p-value combination method
  (Wilson, 2019 - https://doi.org/10.1073/pnas.1814092116)

- Addition of a new dataset (hg19pvalues) containing a set of
  gene p-values from Giakountis et al.,
  https://doi.org/10.1016/j.celrep.2016.05.038 to be used in the
  main vignette as well as a playground for p-value combination

BUG FIXES

- Major bug fix! The order of p-values returned by statDss was
  inconsistent with the rest of the stat* functions.

- Other minor bug fixes.

[MetCirc](https://bioconductor.org/packages/MetCirc)
-------

                 Changes in version 1.19.0 (2020-07-02)                 

- adapt to changes in MSnbase (renaming of Spectra to MSpectra)

[MethReg](https://bioconductor.org/packages/MethReg)
-------

                        Changes in version 0.1.0                        

- Added a NEWS.md file to track changes to the package.

[methrix](https://bioconductor.org/packages/methrix)
-------

                       Changes in version 1.4.00                        

- read_bedgraphs() supports bedgraphs files from "Bismark_cov",
  "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap"

- write_bedgraphs() supports output to multiBed and metilene file
  formats

- write_bigwigs() exports methrix object as bigWigs

                       Changes in version 1.2.10                        

- BSgenome ref build parsing error. Issue: #24

- Show which CpGs are missing when loading bedGraphs. Issue: #22

- Faster HDF5 processing. PR: #21

- include SeqStyle option in Manual write_bedgraphs() function. PE: #20

- added the argument for SeqlevelsStyle and trackline in bedgraph. PR:
  #19

[methylKit](https://bioconductor.org/packages/methylKit)
---------

                       Changes in version 1.15.3                        

NEW FUNCTIONS AND FEATURES

- add verbosity argument to methCall to allow for
  complete silencing of output, this will be set by
  processBismarkAln which gained 'verbose' argument

IMPROVEMENTS AND BUG FIXES

- update R requirement to (>= 3.5.0)

- processBismarkAln:
  - add verbose argument
  - update tests
  - update docs

- joinsegmentneighbours:
  - Catch case when only one groups is there

- regionCounts:
  - update docs: mention that the output will be sorted by default
  - Fixes #195

- vignette:
  - add faq paragraph about regionCount to the vignette explaining
  why we order input region
  - update paragraph about methyldackel when to use
  it instead of processBismarkAln
  - Fixes #162

- add extensive test for bedgraph export

- update and export tests for getMethylationStats

- add extensive tests for getCoverageStats

- add extensive tests for normalizCoverage

- update checks for dataSim:
  - if sample ids given check if length(sample.ids) == # replicates

- update tests for getMethylDiff

- unite destranding:
  - make sure high numbers are not used in scientific notation
  - fixes #129
  - switch to data.table based destranding as it is faster and
  more memory efficient for larger number of regions

- methseg:
  - catch error when joined segments produce only single range
  - Fixes #195

- calculateDiffmeth:
  - check if output is equal for memory and tabix file
  - Fixes #189

                       Changes in version 1.15.2                        

IMPROVEMENTS AND BUG FIXES

- check.dbdir: catch if path equals root ("/") dir

- methRead:
  - add missing context argument for call
  to .procBismarkCytosineReport
  - fixes https://github.com/al2na/methylKit/issues/198

[microbiome](https://bioconductor.org/packages/microbiome)
----------

                 Changes in version 2.1.2 (2020-07-01)                  

- Core heatmap labeling improved

- aggregate_top_taxa deprecated

- bimodality and potential_analysis functions fixed

                 Changes in version 2.1.1 (2020-04-06)                  

- Added overlap function

[microbiomeExplorer](https://bioconductor.org/packages/microbiomeExplorer)
------------------

                        Changes in version 1.0.0                        

- Release

[MicrobiotaProcess](https://bioconductor.org/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.1.13                        

- removed retrieve_seq and mapply_retrieve_seq function, since these
need internet. Which might cause time out when check. (2020-10-16,
Fri)

                       Changes in version 1.1.12                        

- modified a bug in diff_analysis.phyloseq: change tax_table(ps) to
ps@tax_table to avoid generate error when tax_table is NULL.
(2020-10-15, Thu)

                       Changes in version 1.1.11                        

- update the examples of drop_taxa. (2020-10-14, Wed)

                       Changes in version 1.1.10                        

- update ggdiffclade to support data.frame input when reduce is TRUE.
(2020-08-28, Fri)

                        Changes in version 1.1.9                        

- update ggordpoint to fit the usage when user want to set mapping by
manually. (2020-08-25, Tue)

                        Changes in version 1.1.8                        

- get_taxadf, get_alltaxadf and diff_analysis has supported function
datasets or other type datasets. (2020-08-17, Mon)

                        Changes in version 1.1.7                        

- bugfix: cladetext argument has been omitted in ggdiffclade, now it
has been fixed. (2020-08-14, Fri)
- deprecated argument: the size argument controlled the width of line
of tree has been deprecated. The linewd replace it (2020-08-14,
Fri).

                        Changes in version 1.1.6                        

- removeUnkown argument has been replaced with removeUnknown in
ggdiffbox, ggeffectsize, ggdifftaxbar and ggdiffclade. (2020-08-12,
Wed)
- class argument has been replaced with classgroup in diff_analysis.
(2020-08-12, Wed)
- add inward_circular layout in ggdiffclade. (2020-08-12, Wed)

                        Changes in version 1.1.5                        

- ggdifftaxbar supports png, tiff format. (2020-08-10, Mon)
- add stop information to state the class argument in diff_analysis.
(2020-08-10, Mon)

                        Changes in version 1.1.4                        

- add tax_table information to result of get_taxadf. (2020-08-07,
Fri)

                        Changes in version 1.1.3                        

- change according to dplyr (v=1.0.0) (2020-08-05, Wed)
- remove rename_ and group_by_
- modified the angle to 90 in ggdiffclade when layout is slanted or
rectangular (2020-08-05, Wed)

                        Changes in version 1.1.2                        

- fix a bug. When the first rank taxa level (Kingdom) is NA.

[miRSM](https://bioconductor.org/packages/miRSM)
-----

                     Changes in version 1.7.2-1.7.3                     

- Update miRSM.R <2020-09-18, Fri>

                        Changes in version 1.7.1                        

- Update module_Coexpress function <2020-08-21, Fri>

[mitch](https://bioconductor.org/packages/mitch)
-----

                       Changes in version 1.0.12                        

- A fix to problems in >3D plots caused by contrast names

[mixOmics](https://bioconductor.org/packages/mixOmics)
--------

                       Changes in version 6.14.0                        

new features / enhancements / changes

- circosPlot: The radial location of feature names can now be
cutomised using var.adj
- added plot and print methods for nipals ouput (#87)
- all Discriminat Analyses now run solely on mode=regression (#79)
- cim argument change: threshold replaced by cutoff
- nipals and pca with missing values allow skipping reconstitution of
the input matrix
- tune.block.splsda now allows random number seed also for parallel
processing (#72)
- New biplot methods for the pca family (#90)

bug fixes

- plotIndiv: Legend bug which misspecified the groups resolved
- plotIndiv: Legends now ordered as inputted, and not alphabetically
- plot method issue for spca resolved
- plotLoadings.spca bug with var.names now fixed (#81)
- ipca deprecation warning fixed

[MOFA2](https://bioconductor.org/packages/MOFA2)
-----

                       Changes in version 0.99.0                        

- submitted package to Bioconductor

[MOGAMUN](https://bioconductor.org/packages/MOGAMUN)
-------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor.

[motifStack](https://bioconductor.org/packages/motifStack)
----------

                       Changes in version 1.33.14                       

- remove reverse complement alignment for RNA motifs.

                       Changes in version 1.33.13                       

- narrow the predefined font.

                       Changes in version 1.33.12                       

- fix the issue for pcm2pfm if the column counts are not identical.

                       Changes in version 1.33.11                       

- update the markers to allow label + line/rect.

                       Changes in version 1.33.10                       

- add alignment function for AA motifs.

                       Changes in version 1.33.9                        

- fix the issue for AA motif of clusterMotifs

                       Changes in version 1.33.8                        

- fix the issue of symbol position for predefined font

                       Changes in version 1.33.7                        

- update check ghostscript

                       Changes in version 1.33.6                        

- fix the bug in pcm2pssm

                       Changes in version 1.33.5                        

- fix the bug in pcm2pssm

                       Changes in version 1.33.4                        

- add pssm class

                       Changes in version 1.33.3                        

- replace MotIV by matalign.

- add cutoffPval for motifSignature function.

                       Changes in version 1.33.2                        

- add parameter environment for coloredSymbols

                       Changes in version 1.33.1                        

- change grImport, grImport2 and MotIV to Suggests.

- use roxygen2 to generate help files.

[MouseFM](https://bioconductor.org/packages/MouseFM)
-------

                 Changes in version 0.99.9 (2020-08-23)                 

- Check package bugfix in getURL function

                 Changes in version 0.99.8 (2020-08-13)                 

- Getter/Setter for URL of backend server

                 Changes in version 0.99.7 (2020-06-08)                 

- Updated database

- New authors

                 Changes in version 0.99.1 (2020-05-22)                 

- Updated .gitignore

                 Changes in version 0.99.0 (2020-05-21)                 

- Submission to Bioconductor

                 Changes in version 0.1.0 (2020-05-04)                  

- Creation

[msa](https://bioconductor.org/packages/msa)
---

                       Changes in version 1.21.1                        

- changed msaClustalW() examples to run smoothly on Windows with R
  4.0.x

- added warning to msaClustalW() help page regarding cluster="upgma" on
  Windows

                       Changes in version 1.21.0                        

- new branch for Bioconductor 3.12 devel

[MsCoreUtils](https://bioconductor.org/packages/MsCoreUtils)
-----------

                         Changes in version 1.1                         

Changes in 1.1.7

- Rewrite c("left", "right", "inner", "outer") join in C <2020-10-06
Tue>.

Changes in 1.1.6

- Rewrite closest in C <2020-09-24 Thu>.
- Fix #65 and #66.

Changes in 1.1.5

- Add ... to functions to join and compare peaks; see also #131.

Changes in 1.1.4

- Change references to Feature to QFeatures <2020-07-14 Tue>
- Ensure closest accept just argument tolerance of length 1 or
length(x); see also #61, PR #62 <2020-08-07 Thu>.
- The tolerance argument in closest should now be of length 1 or of
length(x) (was length(table) before) <2020-08-20 Thu>.

Changes in 1.1.3

- For an empty table closest and common return a vector of length x
with NA or FALSE, respectively (instead of 1 and TRUE). Fixes #55
<2020-06-18 Thu>.
- closest and common ignore NA in table <2020-06-19 Fri>.
- Fix rbindFill for single data.frame or DataFrame as input
<2020-06-23 Tue>.

Changes in 1.1.2

- New colCounts() aggregation function <2020-05-27 Wed>.

Changes in 1.1.1

- Add some popular distance/similarity metrices: ndotproduct
neuclidean navdist nspectraangle; see also PR #33.

- Add deprecation note to dotproduct <2020-05-22 Fri>.

Changes in 1.1.0

- Bioconductor devel version (Bioc 3.12)

[MSEADbi](https://bioconductor.org/packages/MSEADbi)
-------

                       Changes in version 0.99.0                        

- Package released

[msImpute](https://bioconductor.org/packages/msImpute)
--------

                       Changes in version 0.99.26                       

- update doc for msImpute

                       Changes in version 0.99.25                       

- fix typo in msImpute man page

                       Changes in version 0.99.24                       

- Bug fix in the internal function l2bary

                       Changes in version 0.99.23                       

- selectFeatures and msImpute now use information theoretic
  approaches to find informative features for MAR/MNAR diagnosis
  and estimation of optimal rank, respectively.

- lambda in msImpute is now estimated from the data, using the
  bayesian interpretation of this shrinkage operator.

- msImpute can be run in three modes: "v1" is the original
  implementation of softImpute-als algorithm, "v2" is the
  enhanced low-rank estimation implemented in this version
  update, "v2-mnar" is adaptation of low-rank models for MNAR
  data. More details about methods in documentation.

                       Changes in version 0.99.22                       

- Submitted to Bioconductor

[MSnbase](https://bioconductor.org/packages/MSnbase)
-------

                        Changes in version 2.15                         

Changes in 2.15.7

- Fixed some consistency in the processing reporting (contributed by
@stanstrup).

Changes in 2.15.6

- Add extractSpectraData function to extract all data from an MSnExp
or Spectra object as a DataFrame to be used with the Spectra
package.

Changes in 2.15.5

- combineSpectraMovingWindow gains parameter ppm to allow definition
of m/z relative differences.

Changes in 2.15.4

- Add splitByFile,OnDiskMSnExp method to provide a more efficient
splitting.

Changes in 2.15.3

- Rename Spectra class to MSpectra.
- Rename Chromatograms class to MChromatograms.

Changes in 2.15.2

- Add MSnbase2 citation entry <2020-05-15 Fri>

Changes in 2.15.1

- Fix: don't read mzTab twice <2020-05-05 Tue>

Changes in 2.15.0

- Bioconductor devel version (Bioc 3.12)

[MSPrep](https://bioconductor.org/packages/MSPrep)
------

                 Changes in version 0.99.0 (2020-10-02)                 

- Submitted to Bioconductor

[msPurity](https://bioconductor.org/packages/msPurity)
--------

                       Changes in version 1.15.1                        

- Update dev to match bug fixes in master

                       Changes in version 1.14.1                        

- Update of rdpc algorithm (see
  https://github.com/computational-metabolomics/msPurity/issues/78)

- Update of align algorithm (see
  https://github.com/computational-metabolomics/msPurity/issues/79)

- Fix for spectralMatching of type 'scan' previously incorrectly
  outputing no matches

[MSstatsConvert](https://bioconductor.org/packages/MSstatsConvert)
--------------

                       Changes in version 0.99.0                        

- Added logging utitilies: MSstatsLogsSettings,
MSstatsSaveSesionInfo.
- Added functions for importing and cleaning data: MSstatsImport,
MSstatsClean.
- Added functions for preprocessing: MSstatsMakeAnnotation,
MSstatsPreprocess, MSstatsBalancedDesign.
- Added vignette.

                        Changes in version 0.0.1                        

- Added a NEWS.md file to track changes to the package.

[MSstatsPTM](https://bioconductor.org/packages/MSstatsPTM)
----------

                 Changes in version 0.99.0 (2020-09-28)                 

- Submitted to Bioconductor

[MSstatsTMT](https://bioconductor.org/packages/MSstatsTMT)
----------

                 Changes in version 1.6.6 (2020-10-13)                  

- Fix the bug in converters due to fractions with same mean, sum
  and max values

- Fix the bug in converters due to summaryforMultipleRows

- Fix the bug in OpemMS converter due to duplicated rows

                 Changes in version 1.6.3 (2020-06-05)                  

- Allow NA in the annotation file

                 Changes in version 1.6.2 (2020-06-02)                  

- Fix the bug in proteinSummarization() function and make sure
  the input to dataProcess is data.frame

                 Changes in version 1.6.1 (2020-05-10)                  

- Update groupComparisonTMT() to make predictions for every
  protein

[MSstatsTMTPTM](https://bioconductor.org/packages/MSstatsTMTPTM)
-------------

                 Changes in version 0.99.0 (2018-09-21)                 

- Submitted to Bioconductor

[MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.16.0                        

New features

- Coercion methods from list/List to MultiAssayExperiment method now
available.

Bug fixes and minor improvements

- Provide more details in documentation for mergeReplicates
- Improved documentation for accessor function return values, helper
function examples (@llrs, #281)
- Fixed bug when using longFormat with character assay matrices
(@jonocarroll, #282)

[multicrispr](https://bioconductor.org/packages/multicrispr)
-----------

                 Changes in version 0.99.0 (2020-04-24)                 

- Submitted to Bioconductor

[multiGSEA](https://bioconductor.org/packages/multiGSEA)
---------

                       Changes in version 0.99.0                        

- Package added to Bioconductor

[multiHiCcompare](https://bioconductor.org/packages/multiHiCcompare)
---------------

                 Changes in version 1.7.1 (2020-08-28)                  

- Track NEWS

- Bump version

- Update R version dependency from 3.5.0 to 4.0

- Add exported objects to NAMESPACE that caused BiocCheck error

- Add CITATION, also to README.md file

- Update DESCRIPTION
  o Add Mikhail Dozmorov as Maintainer
  o Add URL and BugReports fields

[musicatk](https://bioconductor.org/packages/musicatk)
--------

                Changes in version 0.99.10 (2020-10-02)                 

PRE-RELEASE

- Pre-release version of the musicatk package to be submitted to
Bioconductor

[NanoMethViz](https://bioconductor.org/packages/NanoMethViz)
-----------

                        Changes in version 1.0.0                        

- Initial Bioconductor release

[NBAMSeq](https://bioconductor.org/packages/NBAMSeq)
-------

                 Changes in version 1.5.1 (2020-05-04)                  

- Updated citation

[ncRNAtools](https://bioconductor.org/packages/ncRNAtools)
----------

                 Changes in version 0.99.9 (2020-10-01)                 

- Reverted changes introduced in 0.99.6

- Set SSL security level to 1 in Linux environment to bypass a bug in
  OpenSSL 1.1.1

                 Changes in version 0.99.6 (2020-09-29)                 

- Changed all URLs to http to prevent issues when accessing websites
  via HTTPS due to a bug in OpenSSL 1.1.1

                 Changes in version 0.99.4 (2020-09-07)                 

- Introduced changes suggested during Bioconductor review

                 Changes in version 0.99.0 (2020-07-24)                 

- Submitted to Bioconductor

[Nebulosa](https://bioconductor.org/packages/Nebulosa)
--------

                       Changes in version 0.99.94                       

- Color scale based now on scale_color_viridis_c()

                       Changes in version 0.99.93                       

- Stable version for bioconductor release

                       Changes in version 0.99.92                       

- Sets up BiocCheck requirements

[netDx](https://bioconductor.org/packages/netDx)
-----

                        Changes in version 1.1.4                        

Changes

- New functionality to smooth mutations over interaction, starting from
  sparse

somatic mutations

- BiocFileCache usage update

[NewWave](https://bioconductor.org/packages/NewWave)
-------

                 Changes in version 0.99.5 (2020-09-08)                 

- Modification to optimd function to prevent decreasing of the
  likelihood function when appling mini-batch with commonwise
  dispersion

                 Changes in version 0.99.0 (2020-09-08)                 

- Submit to bioconductor

                 Changes in version 0.2.0 (2020-09-07)                  

- Matrix and Delayed array framework complete running

[ngsReports](https://bioconductor.org/packages/ngsReports)
----------

                        Changes in version 1.5.5                        

- Added importSJ for importing SJ.out.tab files from the aligner STAR

- Added plotType = "residuals" for plotSeqContent()

                        Changes in version 1.5.4                        

- Allowed for ignoring basename() calls in all importNgsLogs functions

- Added cumulative GC plot to plotGcContent() as plotType = "cdf"

- Changed plotting option name from cumulative to cdf for
  plotSeqLengthDistn()

- Announced deprecation of runFastQC()

                        Changes in version 1.5.3                        

- Fixed bug in plotOverrep

[Omixer](https://bioconductor.org/packages/Omixer)
------

                        Changes in version 0.99                         

- Initial release to Bioconductor.

- Added NEWS file.

- Changed image paths in vignette.

[OmnipathR](https://bioconductor.org/packages/OmnipathR)
---------

                 Changes in version 1.99.0 (2020-10-03)                 

- New tests
- Many little bugfixes
- Updated package metadata
- Preparation for Bioconductor 3.12

                 Changes in version 1.3.4 (2020-08-04)                  

- License and password handling
- Can be directed to different server by options

                 Changes in version 1.3.1 (2020-06-05)                  

- Addition of the package website.
- Modification of the pdf vignette for an html one.
- Major refactoring of functions.
- Addition of the intercellular network function.
- Modification for the intercellular categories.

[OncoSimulR](https://bioconductor.org/packages/OncoSimulR)
----------

                 Changes in version 2.19.2 (2020-06-03)                 

- Updated ctb: added Magellan and exprTk and clarified contributions.

                 Changes in version 2.19.1 (2020-05-05)                 

- Occasional failures of @test.sample-prob.R#42
  in Windows 386

                 Changes in version 2.19.0 (2020-05-05)                 

- Bumped version to match current Biocdevel.

[onlineFDR](https://bioconductor.org/packages/onlineFDR)
---------

                        Changes in version 1.8.0                        

MODIFICATIONS

- harmonised inputs to algorithms

- added unit tests

- updated authors

[openCyto](https://bioconductor.org/packages/openCyto)
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

[OrganismDbi](https://bioconductor.org/packages/OrganismDbi)
-----------

                       Changes in version 1.32.0                        

BUG FIXES

- (v. 1.31.1) Load OrganismDb objects even if not on search
  path. See https://support.bioconductor.org/p/134141/

[OUTRIDER](https://bioconductor.org/packages/OUTRIDER)
--------

                        Changes in version 1.7.1                        

- Move to S3/S4 methods to be compatible with FRASER

- Due to the S3/S4 changes minor changes in the argument names happend
  mainly ods -> object or x

- Bugfixes: #29

[padma](https://bioconductor.org/packages/padma)
-----

                  Changes in version 0.99 (2018-05-15)                  

- Substantial changes to prepare for Bioconductor submission
  o added S4 classes / methods / generics
  o integration with MultiAssayExperiment

                  Changes in version 0.1 (2020-01-09)                   

- Initial development version

[pageRank](https://bioconductor.org/packages/pageRank)
--------

                 Changes in version 0.99.0 (2020-09-21)                 

- Submitted to Bioconductor

[PAST](https://bioconductor.org/packages/PAST)
----

                 Changes in version 1.4.1 (2020-04-30)                  

- bug fixes and new minimum R version

[pathview](https://bioconductor.org/packages/pathview)
--------

                       Changes in version 1.29.1                        

- developmental version now available on GitHub at:
  https://github.com/datapplab/pathview

- fixed warning on "requireNamespace(pkg.name)" in geneannot.map.R
  and sim.mol.data.R. Now "import" from org.Hs.eg.db, instead of
  "suggest" it. Also import more functions and class from
  AnnotationDbi.

- prevented potential error caused by class(exprs) == "data.frame"
  or "class()==" if conditions in multiple functions. class(exprs)
  now returns a vector of length 2, which caused the error.

[pcaExplorer](https://bioconductor.org/packages/pcaExplorer)
-----------

                       Changes in version 2.16.0                        

Other notes

- Replaced dependency from d3heatmap with the functionality of
heatmaply

[PCAtools](https://bioconductor.org/packages/PCAtools)
--------

                        Changes in version 2.2.0                        

- added support to overlay variable loadings arrows on biplot

- better positioning of PC labels in triangular pairs plot

- added functionality to encircle groups of samples

- added functionality to draw stat ellipses with confidence intervals

- users can now add legend titles and these will be added by default

[peakPantheR](https://bioconductor.org/packages/peakPantheR)
-----------

                 Changes in version 1.3.4 (2020-10-23)                  

- Add missing function export for GUI

                 Changes in version 1.3.3 (2020-10-11)                  

- Shiny Graphical User Interface

                 Changes in version 1.3.2 (2020-09-29)                  

- Unittest numerical precision correction for Ubuntu 20.04.1 LTS

                 Changes in version 1.3.1 (2020-09-21)                  

- Retention time correction procedure based on reference compounds

- Exponentially Modified Gaussian (EMG) peak fitting

[PepsNMR](https://bioconductor.org/packages/PepsNMR)
-------

                 Changes in version 1.7.1 (2020-07-21)                  

BUG CORRECTION 

- Debugged ReadFid when the type of raw data is double (i.e. if
  DTYPA==2)

[periodicDNA](https://bioconductor.org/packages/periodicDNA)
-----------

                        Changes in version 0.3.2                        

- IMPORTANT:

- Implemented data-raw for reproducibility

- MINOR:

- Changed xlim of norm. distr. plot in plotPeriodicityResults()

                 Changes in version 0.3.1 (2020-05-05)                  

- IMPORTANT:

- rollmean(k=3) is now applied before normalisation as well, on
the raw distribution vector
- plotPeriodicityResults() output returns one single plot (with
cowplot)
- getPeriodicityTrack() now returns the Rle
- Improved plotting functions -now show shuffled for
plotPeriodicityResults()
- Added ggplot2 theming

- MINOR:

- Changed many variable names (all to snake_case)
- sampleGRanges is now full-fledged function (GRanges,
DNAStringSet, character and BSgenome methods)
- sampleGenome is an alias for sampleGRanges.character
- Added sacCer3 to getPeriodicity BSgenomes
- Added DNAString method for getPeriodicity
- Added a vignette describing the internal steps
- Clarified user-level functions in README
- Added ce11_TSSs data
- Renamed generateperiodicitytrack as getPeriodicityTrack
- Renamed variables in getFPI and getPeriodicity
- Created a utility char2BSgenome()

                 Changes in version 0.3.0 (2020-05-03)                  

- Added tests
- Added getFPI function
- cleaned-up functions
- cleaned-up function dependencies
- Added toy data
- Added vignette

                 Changes in version 0.2.1 (2020-03-04)                  

- Added Travis build check
- Simplified README.md

[PGA](https://bioconductor.org/packages/PGA)
---

                       Changes in version 1.19.1                        

- set2key, set2keyv, and key2 of data.table were deprecated and had
  been removed.

[phantasus](https://bioconductor.org/packages/phantasus)
---------

                        Changes in version 1.9.5                        

- Support for preloaded sessions

                        Changes in version 1.9.3                        

- Changed Dockerfile to nginx+rApache (default ports are changed)

[PharmacoGx](https://bioconductor.org/packages/PharmacoGx)
----------

                       Changes in version 2.1.12                        

- Added experimental support for a new class, the LongTable, for
storing the sensitivity data in a PharmacoSet.
- Because this is not well tested, we have left not updated the PSets
available via the downloadPSets function.
- Instead we have provided a convenient function,
convertSensitivitySlotToLongTable, which takes in a PharmacoSet
object, converts the data in the @sensitivty slot to a LongTable and
returns an updated PharmacoSet object.
- The LongTable class will be used in the future to allow
PharmacoSets
to store treatment response experiments with multiple drugs or
cell-lines, greatly expanding the variety of data which can be
stored in a PharmacoSet object.
- For more details on the LongTable class, please see the vignette in
the CoreGx package.

[PhosR](https://bioconductor.org/packages/PhosR)
-----

                        Changes in version 0.1.5                        

- Significant updates and addition to the documentation
- Major changes to code writing style to comply with BioC
- Updated TC example data

[PhyloProfile](https://bioconductor.org/packages/PhyloProfile)
------------

                       Changes in version 1.3.11                        

- Add filter number of co-orthologs to parseInfoProfile()

                       Changes in version 1.3.10                        

- Calculate percentage of present taxa after filtering of var1
  and var2

- Fixed filter when working with high taxonomy ranks

                        Changes in version 1.2.8                        

- Added new NCBI taxonomy ranks (e.g. biotype, isolate,
  pathogroup, ...)

- Added function to reset taxonomy data

                        Changes in version 1.2.7                        

- Solved problem with new NCBI taxonomy rank "clade" by replace
  them with "norank"

                        Changes in version 1.2.6                        

- Fixed bug customized profile of subset of taxa not clickable

                        Changes in version 1.2.4                        

- Fixed bug checking invalid taxon IDs

                        Changes in version 1.2.1                        

- Fixed bug in rankIndexing and processNcbiTaxonomy

- Improved check for invalid input taxon IDs

[Pigengene](https://bioconductor.org/packages/Pigengene)
---------

                Changes in version 1.15.20 (2020-09-28)                 

Changes in existing functions

- When automatically computing the threshold in the consensus
  function, we now do not allow it to be more than 1 to be
  compatible with bnlearn 4.6.1.

                Changes in version 1.15.16 (2020-09-22)                 

Changes in existing functions

- The modules can now be determined in the module.heatmap
  function using the new mes argument.

                Changes in version 1.15.14 (2020-09-08)                 

Bug Fixes

- combine.networks now really removes the big TOM file when
  doRemoveTOM=TRUE.

                Changes in version 1.15.12 (2020-06-22)                 

Changes in existing functions

- The DiseaseColors argument in the plot.pigengene function can
  now be set automatically.

Bug Fixes

- pheatmap.type now works fine even when the number of samples in
  a condition is only 1.

[pipeComp](https://bioconductor.org/packages/pipeComp)
--------

                Changes in version 0.99.43 (2020-07-20)                 

- there is now the possibility to continue runs despite errors, and
  enlarged merging capacities

                Changes in version 0.99.27 (2020-04-29)                 

- the plotting functions for the scRNAseq clustering pipeline
  (`scrna_evalPlot_DR` and `scrna_evalPlot_clust`) have been replaced
  by more flexible, pipeline-generic functions (`evalHeatmap`) and a
  silhouette-specific plotting function (`scrna_evalPlot_silh`). The
  general heatmap coloring scheme has also been changed to make
  meaningful changes clearer.

                Changes in version 0.99.23 (2020-04-21)                 

- Added SVA-DEA pipeline and vignette

- Standardized results heatmaps for any PipelineDefinition

                 Changes in version 0.99.7 (2020-04-01)                 

- Submitted to Bioconductor

                 Changes in version 0.99.3 (2020-02-07)                 

- made important changes to the scRNAseq pipeline, greatly simplifying
  the format of the output. As a consequence, results produced with
  older version of the package are not anymore compatible with the
  current version's aggregation and plotting functions.

[plyranges](https://bioconductor.org/packages/plyranges)
---------

                        Changes in version 1.9.3                        

- minor spelling and layout fixes to vignette,
- @PeteHaitch corrected table layout in vignette

                        Changes in version 1.9.2                        

- minor documentation fixes

                        Changes in version 1.9.1                        

- Spencer Nystrom has made several significant contributions to the
join_nearest family of functions:

1.  join_nearest_(x, y, ..., distance = TRUE) family of functions
now takes a new argument, distance, which allows the user to add
a column for the distance of the nearest range y to that in x.
2.  add_nearest_distance_(x, y, ...) family of functions, which will
add a new metadata column to the x ranges object which contains
the distance to its nearest neighbor in y. If there are no
nearest neighbors, the new column will be given a missing value.

[polyester](https://bioconductor.org/packages/polyester)
---------

                       Changes in version 1.99.3                        

- NB function now exported

- note that version 1.99.3 on GitHub was version 1.1.0 on Bioconductor.

                       Changes in version 1.99.2                        

- bug fix in fragment generation (last 2 bases of transcript were never
  sequenced)

[POMA](https://bioconductor.org/packages/POMA)
----

                       Changes in version 0.99.45                       

- PomaOutliers bug fixed
- New references in PomaRankProd help
- PomaLimma and PomaUnivariate bug related with one covariate
analysis
fixed
- The elbow method to calculate the optimum number of clusters has
been added in PomaClust function

                       Changes in version 0.99.37                       

- POMA has been accepted to Bioconductor!
- PomaRankProd minor bug fixed
- Authors updated
- Bioconductor logo

                       Changes in version 0.99.33                       

- All Bioconductor issues in the review process have been addressed
- POMA EDA vignette added

                       Changes in version 0.99.16                       

- BiocCeck requirements fixed
- Examples added in functions
- pkgdown files removed from master branch

                       Changes in version 0.99.0                        

- POMA is now submitted to Bioconductor!
- All features implemented work as expected
- All tests finished
- Achieved desired coverage (>95%)
- Vignettes and documentation are ready for the first release

[preciseTAD](https://bioconductor.org/packages/preciseTAD)
----------

                 Changes in version 0.99.6 (2020-08-30)                 

- Major simplification of the preciseTAD function, many arguments are
  kept default

- Major vignette rewrite

- Add pkgdown web site

- Update DESCRIPTION
  o Package description

- Start tracking NEWS

                 Changes in version 0.99.4 (2020-07-14)                 

- Preparing for release

                 Changes in version 0.99.0 (2020-06-16)                 

- Initial version

[proActiv](https://bioconductor.org/packages/proActiv)
--------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[proDA](https://bioconductor.org/packages/proDA)
-----

                         Changes in version 1.3                         

- Escape result coefficient names if result_names() is called and they
  contain potentially
  problematic characters (see issue #4)

[progeny](https://bioconductor.org/packages/progeny)
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

[pRolocGUI](https://bioconductor.org/packages/pRolocGUI)
---------

                        Changes in version 1.99                         

CHANGES IN VERSION 1.99.2

- fixed warning in profiles plot
- fixed colour matching in factet plot

CHANGES IN VERSION 1.99.1

- Apps updated to new shinydashboardPlus style
- Options to save plots as high res images
- User set colourpicker for creating figures
- Support for data export to csv
- Support batch import searching
- PCA app has been replaced with explore app
- New aggregation app
- Classify app removed

[psichomics](https://bioconductor.org/packages/psichomics)
----------

                       Changes in version 1.14.4                        

- Improve file browser dialog:
- Now correctly allows to select a single file in macOS
- Now correctly allows to select a directory in Windows
- Correctly set default location when opening a file in Linux
- Internal improvements
- Data loading:
- Improve performance when loading compressed data (GZ and BZ2
files)
- Increase optional verbosity while loading files
- Speed up loading of VAST-TOOLS' inclusion levels
- Fix loading of VAST-TOOLS' inclusion levels sometimes returning
an error regarding duplicated row names
- Fix issues with missing values as alternative splicing event
complexity in VAST-TOOLS' inclusion levels
- Copy-edit tutorials on loading user-provided data and using the
command-line version

                       Changes in version 1.14.3                        

- Fix unit tests due to changes in Ensembl API response

                       Changes in version 1.14.2                        

Support for loading more data formats

- getSplicingEventData(): get a table with all information for a
given
dataset of alternative splicing quantification data
- VAST-TOOLS inclusion levels and gene expression tables:
- Import VAST-TOOLS' output files (inclusion levels, cRPKMs and
gene read count tables); in the visual interface, use the
User-provided data loading panel and load a specific alternative
splicing quantification or load the folder containing such file;
in the command-line interface, use loadLocalFiles()
- To parse splicing event from VAST-TOOLS in the command-line
interface, use for instance parseSplicingEvent("HsaEX0007927",
data=VASTTOOLSpsi)
- Support for importing inclusion tables with arbitrary alternative
splicing event identifiers (information for these events will not be
available, such as event type and cognate gene)
- SRA metadata:
- Automatically load metadata from the SRA Run Selector (usually
comes in files named SraRunTable.txt) avoiding the need to use
prepareSRAmetadata() first
- GTEx data loading (loadGtexData()):
- Support for loading GTEx data from V8, V7, V6 or V4 releases
- Organise GTEx data into folders named based on release version
- Copy-edit tutorial on loading user-provided data

New features

- plotSplicingEvent():
- Alternative splicing diagrams now render automatically instead
of showing SVG code (if printing more than one event, a table is
displayed with event identifiers and respective diagrams)
- Plot intron retention events (e.g. from VAST-TOOLS)
- plotLibrarySize(): plot library size from gene expression data
- Alternative splicing quantification filtering (visual interface):
- New panel to allow filtering already loaded alternative splicing
quantification based on event types, samples, PSI statistics
(such as median, variance and range) and cognate genes (some
filtering steps may not be available when using user-provided
tables with inclusion levels)
- Filter data based on individual PSI values (also available when
quantifying alternative splicing)
- Filter VAST-TOOLS events based on its quality scores for read
coverage
- Toggle specific statistics based on alternative splicing
quantification values when filtering
- Preview the original and filtered events in an easy-to-use plot
- Gene expression filtering:
- Toggle specific statistics based on gene expression values when
filtering (visual interface)
- Arguments used in filterGeneExpr() and normaliseGeneExpression()
are now returned as attributes of those functions
- Differential analysis:
- Tooltips of volcano plots now include diagram of alternative
splicing events
- Gene, transcript and protein annotation (visual interface):
- Suggest genes based on loaded data
- Change species and assembly when fetching information
- Improve interface of PubMed article query, including persistent
user-inputted tags even when changing selected gene
- Improve relevance of results from PubMed articles
- Much faster calculation of row/column-wise means, medians,
variances
and ranges (helpful for plotting statistics of large datasets)
- Include Dockerfile to create Docker images based on code revisions

Bug fixes and minor changes

- Progress bar is now animated as in previous versions
- Data loading (visual interface):
- Allow to discard selected files in file browser input elements
- When creating groups by sample index/identifiers, suggest sample
names not only obtained from sample information, but also from
alternative splicing quantification and gene expression datasets
- "Browse..." button now opens file browser to select folder where
data is stored (as expected) in GTEx and SRA panels
- Only pre-create groups of genes (based on literature-based gene
lists) if at least one of its genes exists in any of the loaded
datasets
- Correctly parse gene symbols containing underscores
- Fix gene expression summary plots not showing in specific
situations
- Fix library size plot not working properly and causing rendering
issues
- Fix settings used to quantify alternative splicing not showing
up
- Show errors raised while reading a file (e.g. if file is too big
for available memory)
- Show alert if no GTEx data options are selected
- Show filename of the file used to load gene expression and
alternative splicing data from GTEx and SRA
- Show helpful context messages in panel interfaces
- Improve visual interface and minor copy-editing
- Local data loading:
- Support loading data from GTEx V8 or previous releases
- Fix bad formatting of help tooltips when using shiny 1.4.0 or
newer (visual interface)
- Gene expression filtering:
- Hide message regarding the usage of no design matrix
- Alternative splicing annotation:
- Try to load cached alternative splicing annotation if a timeout
occurs
- Include gene symbols in custom annotations if available
- Account for possible filename changes when parsing annotations
from VAST-TOOLS, rMATS, SUPPA and MISO
- Alternative splicing event selection (visual interface):
- Include loading indicator while searching for events
- Decrease number of operations performed after selecting an event
- Fix crash when changing to an alternative splicing
quantification dataset without the selected splicing event
- Allow to search using an event identifier directly
- Show event identifier instead of prettier identifier to avoid
confusion
- Data grouping (visual interface):
- Simplify group selection interface
- Fix suggested attributes and index/identifiers in group creation
not being cleared when changing to datasets where such data is
unavailable (thus showing the attributes/index/identifiers of
the previous dataset)
- Show an alert when there is no data to create groups
- Dimensionality reduction:
- performPCA() and performICA() now directly raise errors (instead
of simply capturing and returning them)
- PCA: loading plots now show parsed information of events
(cognate gene, event type and genomic position) if available
- PCA: the correct splicing event is now selected when clicking on
any loadings in the loadings plot (visual interface)
- Differential analyses (visual interface):
- Use group colours in density plots of differential expression
table
- Fix occasional crash when performing differential analyses with
different number of groups during the same session
- Tooltips of volcano plots are now properly positioned in
high-resolution screens
- The correct splicing event is now selected when clicking the
density plots or survival curves in the table
- When creating a group from differential splicing results,
correctly set cognate genes of alternative splicing events
- diffAnalyses(): deprecated psi argument was now removed
- Distribution plots (plotDistribution()):
- After hiding all plot series, hide Y axis (rug plots of the
different groups have arbitrary Y values to easily distinguish
them)
- Rug plots of gene expression density plots are now placed near
the X axis as expected
- Correlation analyses (visual interface):
- Warn when selecting genes that are not available in the selected
gene expression dataset (instead of crashing the app)
- print() extended to better display information on gene list
objects;
e.g. print(getGeneList())
- Fix issues when installing the package:
- Fix error when unit testing in R 4.0 or higher (strings in data
frames are not converted to factors by default)
- Fix R CMD check warning of Unicode symbol translation in Windows
- Fix comparing signed and unsigned integers in Rcpp functions

                       Changes in version 1.14.1                        

Fix unit tests for R 4.0

[PureCN](https://bioconductor.org/packages/PureCN)
------

                       Changes in version 1.20.0                        

NEW FEATURES

- Support for GATK4 GenomicsDB import for mapping bias calculation

- Added --additionaltumors to PureCN.R to provide coverage files
  from additional biopsies from the same patient when available

- PSCBS segmentation now identifies on-target breakpoints first when
  off-target is noisy, thus boosting sensitivity in on-target regions

- Beta-binomial model in runAbsoluteCN now uses the fits in mapping
  bias
  database. We plan to set this as default in upcoming versions and
  appreciate feedback.

SIGNIFICANT USER-VISIBLE CHANGES

- We now check if POP_AF or POPAF is -log10 scaled as new Mutect2
  versions
  do.

- Added support for GERMQ info field containing Phred-scaled germline
  probabilities.

- Detect Mutect2 VCF more reliably

- Updated Mutect2 failure flags: "strand_bias", "slippage",
  "weak_evidence",
  "orientation", "haplotype"

- Removed defunct normal.panel.vcf.file from setMappingBiasVcf

- Removed defunct interval.weight.file from segmentationPSCBS,
  segmentationCBS and processMultipleSamples

- Made calculateIntervalWeights defunct

- Changed default of min.normals in calculateMappingBiasVcf/Gatk4 to 1
  from 2

- Changed default of --signature_databases to
  "signatures.exome.cosmic.v3.may2019" (v3 instead of v2)

- Now warn if recommended -funsegmentation is not used

- Added parallel option for callAmplificationsInLowPurity

- callMutationBurden now uses all non-filtered targets as callable
  region
  when callable is not provided

- plotAbs in chromosome mode now displays wider range of log2 ratios
  (makes it possible to examine outliers)

- Moved vcf.field.prefix from predictSomatic to runAbsoluteCN since it
  now
  adds more fields like prior somatic and mapping bias to the VCF

- Changed default of runAbsoluteCN min.ploidy to 1.4

BUGFIXES

- Fix for crash with CNVkit input when log-ratio contained highly
  negative
  outliers

- Fixed a bug in preprocessIntervals/IntervalFile.R when input
  contained
  overlapping and stranded intervals

- Fix for crash when GC-correction is attempted on empty coverage (for
  example
  off-target region without any off-target reads)

- Fix for crash when VCF FA field contained missing values

- Fix for a bug in callAmplificationsInLowPurity that can cause a wrong
  chromosome percentile

[QFeatures](https://bioconductor.org/packages/QFeatures)
---------

                        Changes in version 0.99                         

QFeatures 0.99.4

- Fix: improved nNA with new implementation and additional unit tests
<2020-10-23 Fri>

QFeatures 0.99.3

- New feature: the longFormat function returns a long DataFrame with
quantitative data along with metadata (see #116) <2020-10-8 Thu>
- New feature: the rowData method returns a list containing the
rowData for all assays (see #86) <2020-09-16 Wed>
- Keep colnames when reading a single column assay (see #108)
<2020-09-09 Wed>

QFeatures 0.99.2

- Added infIsNA.
- Add Christophe Vanderaa as an author.

QFeatures 0.99.1

- Address comments Bioconductor review (see submission issue for
details.

QFeatures 0.99.0

- Bioconductor submission

[qPLEXanalyzer](https://bioconductor.org/packages/qPLEXanalyzer)
-------------

                        Changes in version 1.7.3                        

- Fix limits of scale in corrPlot to c(0, 1)

- Fix bug in regress intensities where control columns were not being
  removed

                        Changes in version 1.7.1                        

- Change Vignette to HTML

- Add textsize argument to corrPlot

- Allow missing values in scaling normalisations

[QuasR](https://bioconductor.org/packages/QuasR)
-----

                       Changes in version 1.30.0                        

NEW FEATURES

- added binSize argument for qProfile

[RadioGx](https://bioconductor.org/packages/RadioGx)
-------

                        Changes in version 1.0.1                        

- Abstract further functionality from RadioGx into CoreGx dependency
package
- Harmonize function and API design to be in line with PharmacoGx R
package

                        Changes in version 1.0.0                        

- RadioGx archived on CRAN for submission to Bioconductor
- Molecular profile data has been migrated from ExpressionSet to
SummarizedExperiment to allow incorporation of multiple assays as
well as to be in line with Bioconductor best practices
- Wrote a vignette to highlight simple usage of the package

[RaggedExperiment](https://bioconductor.org/packages/RaggedExperiment)
----------------

                       Changes in version 1.14.0                        

Bug fixes and minor improvements

- Update package to changes in MatrixGenerics (new location of the
rowRanges generic)

[randRotation](https://bioconductor.org/packages/randRotation)
------------

                 Changes in version 1.1.1 (2020-09-23)                  

- Function contrastModel() added to support contrast rotation

- Vignette extended

- Deprecate function df_estimate()

[rBiopaxParser](https://bioconductor.org/packages/rBiopaxParser)
-------------

                       Changes in version 2.29.1                        

- Removed functionality for transitive reduction, which was imported
  from nem package (dropped out of Bioconductor).

[RCy3](https://bioconductor.org/packages/RCy3)
----

                       Changes in version 2.10.0                        

- New functions:
  - getCurrentStyle, #15

- Added Sys.sleep to buggy CyREST steps in...
  - createNetworkFromDataFrames, #98
  - importNetworkFromFile
  - importNetworkFromNDEx
  - exportNetworkToNDEx
  - updateNetworkInNDEx
  - importFilters
  - create***Filter
  - applyFilter
  - updateStyleMapping
  - setVisualPropertyDefault

- Refactored getNetworkViewSuid

- Handled special 404 cases in .cyError

- Bug fix #94: added base.url param

- Overhauled error handling and messaging

[ReactomeGSA](https://bioconductor.org/packages/ReactomeGSA)
-----------

                 Changes in version 1.3.6 (2020-09-01)                  

- Improved analyse_sc_clusters for SingleCellExperiment objects to
  define cell grouping using factors and a single character string
  specifying the metada field.

                 Changes in version 1.3.4 (2020-09-01)                  

- Fixed bug in analyse_sc_clusters when processing SingleCellExperiment
  objects.

                 Changes in version 1.3.3 (2020-05-26)                  

- Updated documentation

                 Changes in version 1.3.2 (2020-05-26)                  

- Added workaround for SEC_E_ILLEGAL_MESSAGE error on older Windows
  systems.

                 Changes in version 1.3.1 (2020-04-29)                  

- Added option to hide non-significant pathway in plot_correlations

[reconsi](https://bioconductor.org/packages/reconsi)
-------

                        Changes in version 1.1.2                        

- Revert to nonparametric null estimation
- Random tiebreaking is by default off for reconsi(), but on for
testDAA()

                        Changes in version 1.1.1                        

- Renamed "permZvals" argument to "resamZvals", which is more
accurate
- Allow pi0 to be provided rather than estimated
- Actualized return arguments

[recount3](https://bioconductor.org/packages/recount3)
--------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- Documentation website is now available at
http://LieberInstitute.github.io/recount3/. It gets updated with
every commit on the master branch (bioc-devel) using GitHub Actions
and pkgdown.

[recountmethylation](https://bioconductor.org/packages/recountmethylation)
------------------

                       Changes in version 0.99.0                        

- Added getdb functions for database file download and load

- Added servermatrix function to get latest database file
  metadata

- Added User's Guide and Data Analyses vignettes

- Added metadata and data_analyses.RData to /inst/extdata

                        Changes in version 0.01                         

- Added key query/accessor functions.

- Added package vignette.

[recoup](https://bioconductor.org/packages/recoup)
------

                 Changes in version 1.17.1 (2020-05-03)                 

NEW FEATURES

- New annotation based on an SQLite database. Rebuild must be performed
  otherwise the application will break. Download of a ready made
  annotation
  is also available (see vignette).

- More supported genomes

- Quick, high-level profile based on RPM instead of coverage

- Speed and code improvements

BUG FIXES

- Fixed bug resulting in running the same operation multiple times when
  chunking was requested.

[RegEnrich](https://bioconductor.org/packages/RegEnrich)
---------

                       Changes in version 0.99.19                       

- Fix bugs in .regSEA function.

                       Changes in version 0.99.18                       

- Update R version dependency to 4.0.0.

                       Changes in version 0.99.17                       

- Fix bugs of regenrich_diffExpr function.
- Add an example for %>%.
- Use Authors@R [cre] designation.

                       Changes in version 0.99.16                       

- Import magrittr

                       Changes in version 0.99.15                       

- Fix a bug in regenrich_rankScore function.
- Reexport pipe %>%.

                       Changes in version 0.99.14                       

- Add \donttest tags.
- Optimize COEN function.
- Replace bplapply by lapply in pickSoftThreshold2 function.

                       Changes in version 0.99.13                       

- Fix bugs of plotSoftPower on Windows.
- Hide documentations of COEN and GRN functions.
- Remove doParallel package from imports.

                       Changes in version 0.99.12                       

- Remove .Rpoj file.

                       Changes in version 0.99.11                       

- Remove \dontrun tags from all examples.
- Remove lazyData: true from DESCRIPTION.
- Formatting vignettes using BiocStyle.
- Replace NEWS file by NEWS.md file to track changes to the package.
- Some packages in Imports are moved to Depends in DESCRIPTION file.
- DeaSet class inherits SummarizedExperiment class.
- RegenrichSet class inherits DeaSet class.
- TopNetwork class inherits BiocSet class.
- topResult and allResult slots of Enrich object are tibble rather
than data frame.
- Score class inherits tibble.
- The show methods for DeaSet, TopNetwork, Enrich, Score, and
RegenrichSet object have been optimized.

                       Changes in version 0.99.1                        

- Submission to Bioconductor.

[regionReport](https://bioconductor.org/packages/regionReport)
------------

                       Changes in version 1.23.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- No longer suggest DESeq since this package has been deprecated.

[regutools](https://bioconductor.org/packages/regutools)
---------

                        Changes in version 1.1.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- The package CITATION file has been updated now that the online
paper
on regutools is available at
https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btaa575/5861528
with DOI 10.1093/bioinformatics/btaa575 and PMID: 32573705.

DOCUMENTATION UPDATES

- roxygen version update was included and documentation was build
under the new roxygen version.

                        Changes in version 1.1.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- The package CITATION file has been updated now that the pre-print
on
regutools is available at
https://www.biorxiv.org/content/10.1101/2020.04.29.068551v1 with DOI
https://doi.org/10.1101/2020.04.29.068551.

[ResidualMatrix](https://bioconductor.org/packages/ResidualMatrix)
--------------

                        Changes in version 1.0.0                        

- New package ResidualMatrix, containing the ResidualMatrix class
  migrated from BiocSingular.

[rfaRm](https://bioconductor.org/packages/rfaRm)
-----

                 Changes in version 1.1.10 (2020-10-13)                 

- Added missing item (alignmentEndPositionQuery) to the output list in
  the documentation of rfamSequenceSearch

                 Changes in version 1.1.9 (2020-10-13)                  

- Fixed a bug that caused sequence searches with sequences longer than
  10000 bases to crash

                 Changes in version 1.1.6 (2020-10-01)                  

- Set SSL security level to 1 in Linux environment to bypass a bug in
  OpenSSL 1.1.1

                 Changes in version 1.1.5 (2020-09-30)                  

- Reverted changes introduced in 1.1.4

                 Changes in version 1.1.4 (2020-09-29)                  

- Changed all URLs to http to prevent issues when accessing websites
  via HTTPS due to a bug in OpenSSL 1.1.1

                 Changes in version 1.1.3 (2020-07-18)                  

- Fixed a bug that caused sequence searches to crash when applying clan
  competition filter and no associated clan was found for some Rfam
  families

                 Changes in version 1.1.2 (2020-07-15)                  

- Fixed a bug that caused sequence searches to crash when applying clan
  competition filter and a match was found in the - strand

                 Changes in version 1.1.1 (2020-07-08)                  

- Added the possibility to filter hits of a sequence search by clan
  competition

- Sequences larger than 10000 nucleotides can now be provided as input
  for sequence-based searches of the Rfam database

[rGREAT](https://bioconductor.org/packages/rGREAT)
------

                       Changes in version 1.21.1                        

- parseRegionGeneAssociationFile(): parsing also considers positions
  with no sign.

[rgsepd](https://bioconductor.org/packages/rgsepd)
------

                        Changes in version 1.21                         

- changing some less than signs to less than or equal to, to ensure we
  get results.

- planning a refactor to handle other gene ID systems and example
  datasets. not done

[rhdf5](https://bioconductor.org/packages/rhdf5)
-----

                       Changes in version 2.34.0                        

NEW FEATURES

- Added support for read access to files in Amazon S3 buckets
  (currently only available on non-Windows platforms).

- Included read and write support for dynamic compression filters
  distributed in rhdf5filters.

CHANGES

- All datasets written with h5write() now have the attribute
  rhdf5-NA.OK added to them.  This is used to indicate that rhdf5
  was used to create the file and that the user does not need to
  be informed that specific values will be mapped to NA in R.

BUG FIXES

- Fix bug in H5Dget_storage_size() where the wrong C function was
  called.

- NA values in logical datatypes are now preserved when written
  and read back into R
  (https://github.com/grimbough/rhdf5/issues/58).

- Fixed error when trying to write a vector containing only empty
  strings (https://github.com/grimbough/rhdf5/issues/60).

- h5ls() and h5dump() no longer crash when given a file
  containing recursive or duplicated groups
  (https://github.com/grimbough/rhdf5/issues/48).

- Reading compound datasets with at least one 8-bit integer field
  now works (https://github.com/grimbough/rhdf5/issues/71).

- Fixed problem when writing a data.frame containing a column of
  raw values. These columns were ommitted when creating a
  compound dataset.

- Patch early UNPROTECT() when reading a Enum type that could
  cause a segmentation fault
  (https://github.com/grimbough/rhdf5/issues/73)

[rhdf5filters](https://bioconductor.org/packages/rhdf5filters)
------------

                        Changes in version 1.2.0                        

BUG FIXES

- Now passed R CMD config LDFLAGS to the compilation of the BLOSC
filter.

[ribosomeProfilingQC](https://bioconductor.org/packages/ribosomeProfilingQC)
-------------------

                        Changes in version 1.1.6                        

- fix the issue in shiftReads for reads width and readsEndPlot for
seqlevels.

                        Changes in version 1.1.5                        

- Update documentation about the R version requirement.

                        Changes in version 1.1.4                        

- Update authorship.

                        Changes in version 1.1.3                        

- Fix the error "CIGAR is empty after qnarrowing".

                        Changes in version 1.1.2                        

- update vignette to fix the typos.
- update plotTranscript to give more informative warning message.

                        Changes in version 1.1.1                        

- update vignette by adding background colors for better
understanding
the pipeline.

[RNAAgeCalc](https://bioconductor.org/packages/RNAAgeCalc)
----------

                 Changes in version 1.1.2 (2020-09-20)                  

- Added citation and README

                 Changes in version 1.1.1 (2020-06-24)                  

- Replaced "DESeq" with "DESeq2" throughout the package

[rnaEditr](https://bioconductor.org/packages/rnaEditr)
--------

                       Changes in version 0.99.1                        

2020-10-01

Added rnaEditr.Rproj into .gitignore to fix the building error.

                       Changes in version 0.99.0                        

2020-09-22

We are planning to submit to Bioconductor by this Friday, and here is
the link to the issue on gitlab where we resolve BiocCheck() ERRORs,
WARNINGs, and NOTEs:
https://gitlab.com/Jennyzly2016/rnaEditr/-/issues/64

[RNAmodR](https://bioconductor.org/packages/RNAmodR)
-------

                 Changes in version 1.3.5 (2020-08-30)                  

- bugfix for GRangesList as annotation data. Only non-overlapping
  elements are allowed

                 Changes in version 1.3.1 (2020-04-29)                  

- added stats function to access details about number of reads used for
  analysis

[RnaSeqSampleSize](https://bioconductor.org/packages/RnaSeqSampleSize)
----------------

                 Changes in version 1.99.2 (2020-10-22)                 

BUG FIXES

- Change codes to make the package pass all Rcmd check with R 4.0.3

                 Changes in version 1.19.2 (2020-10-02)                 

BUG FIXES

- Change codes to make the package pass all Rcmd check with R 4.0

                 Changes in version 1.19.1 (2020-09-30)                 

BUG FIXES

- Fix a bug in est_power_distribution, which make results incorrect
  when power is very large.

                 Changes in version 1.15.1 (2018-11-12)                 

BUG FIXES

- Fix a bug in power estiamtion (cumsumBorder.cpp), which make results
  incorrect when k not equal to one.

                 Changes in version 1.9.1 (2017-08-03)                  

BUG FIXES

- Fix a bug about estimating genes in pathway power in
  est_power_distribution function, which was caused by update of KEGG
  web site.

                 Changes in version 0.99.8 (2015-04-12)                 

BUG FIXES

- Fix a bug in example of optimize_parameter function, which was caused
  by update of heatmap3 package. Thanks Dan.

                 Changes in version 0.99.7 (2015-04-06)                 

BUG FIXES

- Fix a bug in example of convertId function, which was caused by
  update of BioMart database. Thanks Dan.

                 Changes in version 0.99.5 (2015-04-06)                 

SIGNIFICANT USER-VISIBLE CHANGES

- Minor changes for genes with less than 1 read count;

                       Changes in version 0.99.4                        

NEW FEATURES

- The user friendly web interface was improved

BUG FIXES

- Bug fixed: Max sample size for distribution based sample size
  estimation;

- Bug fixed: Keep consistent for all function names;

SIGNIFICANT USER-VISIBLE CHANGES

- Other improvement based on the comments from Bioconductor reviewer;

                 Changes in version 0.99.3 (2014-11-23)                 

SIGNIFICANT USER-VISIBLE CHANGES

- New parameter: countFilterInRawDistribution,
  selectedGeneFilterByCount;

- Other improvement based on the comments from Bioconductor reviewer;

                       Changes in version 0.99.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette was switched to a BiocStyle;

- Other improvement based on the comments from Bioconductor reviewer;

                 Changes in version 0.99.1 (2014-10-19)                 

SIGNIFICANT USER-VISIBLE CHANGES

- Parameter was used to determine the power of genes below minAveCount;

- Some recommendations from BiocCheck were improved;

                 Changes in version 0.99.0 (2014-10-16)                 

SIGNIFICANT USER-VISIBLE CHANGES

- Submit to Bioconductor.

                 Changes in version 0.5.0 (2014-10-15)                  

NEW FEATURES

- The datasets from TCGA were moved to RnaSeqSampleSizeData package;

- The examples were improved;

- The vignette was generated;

                 Changes in version 0.4.0 (2014-10-10)                  

NEW FEATURES

- Rcpp package was used to make some functions faster;

- The package structure was improved;

- The web interface was improved;

                 Changes in version 0.3.0 (2014-08-01)                  

NEW FEATURES

- Beta approximate method was improved and used in most of the cases;

- Read count and dispersion distribution estimation function was
  improved;

- Estimation power or sample size by read count and dispersion
  distribution was improved.

                 Changes in version 0.2.0 (2014-06-01)                  

NEW FEATURES

- The function of estimation power or sample size by read count and
  dispersion distribution was imported and improved.

                 Changes in version 0.1.7 (2014-03-30)                  

SIGNIFICANT USER-VISIBLE CHANGES

- An approximate estimation method was used to decrease the running
  time when average read count (lambda0) was larger than 20;

NEW FEATURES

- A plotting power curve function plot_power_curve was provided;

                 Changes in version 0.1.6 (2014-03-19)                  

SIGNIFICANT USER-VISIBLE CHANGES

- The code was improved to fit the web interface;

NEW FEATURES

- The function est_power was provided for power estimation;

BUG FIXES

- A bug was fixed;

                 Changes in version 0.1.5 (2014-03-16)                  

SIGNIFICANT USER-VISIBLE CHANGES

- The C code for R function dnbinom was improved to decrease the
  running time greatly;

- The powers for different N in estimating sample size were returned to
  prepare power curve;

                 Changes in version 0.1.4 (2014-03-01)                  

NEW FEATURES

- A user friendly web interface was provided at
  http://cqs.mc.vanderbilt.edu/shiny/RnaSeqSampleSize/;

                 Changes in version 0.1.3 (2014-02-17)                  

SIGNIFICANT USER-VISIBLE CHANGES

- The algorithm was improved to decrease the memory usage;

                 Changes in version 0.1.2 (2014-02-04)                  

SIGNIFICANT USER-VISIBLE CHANGES

- A more detail parameter to sample size table was provided;

- The algorithm was improved to decrease the running time and memory
  usage;

                 Changes in version 0.1.1 (2014-02-04)                  

SIGNIFICANT USER-VISIBLE CHANGES

- The algorithm was improved to decrease the running time greatly;

                 Changes in version 0.1.0 (2014-01-24)                  

SIGNIFICANT USER-VISIBLE CHANGES

- A more detail parameter to sample size table was provided to make
  sample size estimation faster;

- An estimated time to perform sample size estimation will be
  displayed;

BUG FIXES

- A bug fixed when sample size N equals to 1;

                 Changes in version 0.0.1 (2014-01-10)                  

SIGNIFICANT USER-VISIBLE CHANGES

- First version.

[RnBeads](https://bioconductor.org/packages/RnBeads)
-------

                        Changes in version 2.7.1                        

- Added conversion from RnBeadRawSet to RGChannelSet

[rols](https://bioconductor.org/packages/rols)
----

                        Changes in version 2.17                         

CHANGES IN VERSION 2.17.4

- Fix failing unit test <2020-08-26 Wed>

CHANGES IN VERSION 2.17.3

- Fix failing unit test (Matched error messsage changed) <2020-07-01
Wed>

CHANGES IN VERSION 2.17.2

- Fix failing unit test (Terms changed) <2020-06-09 Tue>

CHANGES IN VERSION 2.17.1

- Fix failing unit test (GO description changed) <2020-05-01 Fri>

CHANGES IN VERSION 2.17.0

- Bioconductor devel version (Bioc 3.12)

[ROSeq](https://bioconductor.org/packages/ROSeq)
-----

                       Changes in version 1.99.01                       

- no bug, but removed extra files from man folder

                       Changes in version 1.1.12                        

- updated vignettes and readme file and example in main file
- Publication (Preprint available at BioRxiv): is added

[rrvgo](https://bioconductor.org/packages/rrvgo)
-----

                        Changes in version 1.1.3                        

- Warn user when no scores are provided and thus falling back to term's
  size as the score to rank terms

                        Changes in version 1.1.2                        

- Fix bug that would break `reduceSimMatrix()` when no scores were
  provided

[rScudo](https://bioconductor.org/packages/rScudo)
------

                        Changes in version 1.5.1                        

USER VISIBLE CHANGES

- Updated citation

- An error is returned if there are duplicated sample names in the
  expression data

[rsemmed](https://bioconductor.org/packages/rsemmed)
-------

                       Changes in version 0.99.0                        

- Initial release to Bioconductor

[Rsubread](https://bioconductor.org/packages/Rsubread)
--------

                        Changes in version 2.4.0                        

- 
  The 'isPairedEnd' parameter in featureCounts() is now used to
  check if the type of input reads (paired-end or single-end) is
  correctly specified.

- 
  A new parameter 'countReadPairs' was added to featureCounts to
  specify if read pairs should be counted for paired-end read
  data.

- 
  Changes to the input parameters and output of cellCounts()
  function. CellCounts will generate a Sample Sheet for samples
  included in the scRNA-seq data, based on the sample index set
  name provided by the user. Structure of the List object
  returned by cellCounts() is also simplified. cellCounts() now
  also outputs a BAM file and a gzipped FASTQ file including raw
  reads for each sample.

[RTCA](https://bioconductor.org/packages/RTCA)
----

                     Changes in version 2009-07-13                      

- combineRTCA(list): Additional column is renamed into Plate. The vlues
  is evaluated from list item names. When the list has no name, an
  integer index beginning from 1 is used. Special attentions to list
  partially with names is noted in the documentation.

- parseRTCA(file, dec=".",phenoData, skipWell,...): Example is
  added in the documentation how to import pre-configured
  phenoData. Details section in the documentation is re-written to
  describe the process of parsing.

- RTCA-class: Experiment ID added to RTCA class

- Makefile: add Makefile to simplify common tasks like check and
  install

- plotGridEffect: takes 'column' instead of 'col' as mode
  parameter, and renders the mode as the title of the
  legend. Documentation updated.

- plotRTCA: is removed from the package and is substituted by the
  plot function.

[RTCGAToolbox](https://bioconductor.org/packages/RTCGAToolbox)
------------

                       Changes in version 2.20.0                        

New features

- Added the RNASeq2Gene slot to the FirehoseData class. This data
type
mainly obtains RNASeq v2 raw_counts from the pipeline
(scaled_estimates also available; @mherberg #39)
- Added an accmini example dataset as obtained from getFirehoseData
- getLinks function shows the user some file provenance based on data
requested
- Newly deprecated functions: getDiffExpressedGenes,
getCNGECorrelation, getSurvival and getReport. It is no longer
possible for the maintainer to update these functions in a way that
would benefit users. A transfer of responsibility would be required,
i.e. to another package.
- Vignettes are updated to reflect changes in the codebase.

Bug fixes and minor improvements

- Improvements to internal functions for converting tabular data to
Bioconductor classes
- Missing (NA) seqnames are removed when converting to
RaggedExperiment
- Remove static file dependencies from GitHub and use text inside
function (@DavisWeaver, #34)
- Added default values to helper for making SummarizedExperiment
datasets
- Coerce sample names to character when in (the rare) case they're
numeric
- Added an ellipsis argument to biocExtract for specifying the
names.field in tabular data that will correspond to the row names of
a SummarizedExperiment

[Rtpca](https://bioconductor.org/packages/Rtpca)
-----

                 Changes in version 0.99.2 (2020-06-05)                 

- Submitted to Bioconductor

[rWikiPathways](https://bioconductor.org/packages/rWikiPathways)
-------------

                       Changes in version 1.10.0                        

- Defuncted getColoredPathway

- Doc fix: examples fixed and sample gmt included

                        Changes in version 1.8.5                        

- Doc fix: added reference to clusterProfiler

                        Changes in version 1.8.4                        

- Doc fix: fixed Pathway Analysis vignette

                        Changes in version 1.8.3                        

- Bug fix: GMT parsing #16

                        Changes in version 1.8.2                        

- Doc fix: vignette with BridgeDbR

- Bug fix: fixed rjson replacement

[S4Vectors](https://bioconductor.org/packages/S4Vectors)
---------

                       Changes in version 0.28.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Replaced DataTable class with RectangularData class.

- Replaced DataTable_OR_NULL with DataFrame_OR_NULL class.

- Add parallel_slot_names() generic and methods for Vector derivatives.
  This replaces vertical_slot_names().  The concept of "vertical" and
  "horizontal" slots is now a RectangularData concept only i.e. only
  RectangularData derivatives should define vertical_slot_names() and
  horizontal_slot_names() methods. For RectangularData derivatives that
  are also Vector derivatives, one of the two methods should typically
  be defined as a synonym of parallel_slot_names().  For example
  horizontal_slot_names() now returns parallel_slot_names() on a
  DataFrame derivative and vertical_slot_names() will return
  parallel_slot_names() on a SummarizedExperiment derivative.

- makeClassinfoRowForCompactPrinting() is now exported.

- showAsCell() now trims strings that are > 22 characters.

- Small tweak to show() method for Rle objects. Now it uses
  showAsCell()
  instead of as.character() for more compact display of the run values
  of the Rle object, and for consistency with other show() methods
  (e.g.
  with method for DataFrame objects).

DEPRECATED AND DEFUNCT

- parallelSlotNames() is now defunct (after being deprecated in BioC
  3.11).

BUG FIXES

- Fix coercion from SimpleList to DataFrame.

- Fix bug in showAsCell() on ordinary data frames.

- Make sure showAsCell() works on a list of non-subsettable objects.

[SAIGEgds](https://bioconductor.org/packages/SAIGEgds)
--------

                        Changes in version 1.3.3                        

- set-based tests: burden, ACAT-V

                        Changes in version 1.2.2                        

- update the citation

- work around gcc-10

[sangeranalyseR](https://bioconductor.org/packages/sangeranalyseR)
--------------

                       Changes in version 0.99.1                        

- Base class: SangerReads is designed to store each forward/revers
  reads.

                        Changes in version 0.1.0                        

- Project starts.

[SCATE](https://bioconductor.org/packages/SCATE)
-----

                        Changes in version 1.0.0                        

- submission to Bioconductor

[scCB2](https://bioconductor.org/packages/scCB2)
-----

                       Changes in version 0.99.31                       

- Added extra filtering for barnyard data (a mixture of multiple
species).

                       Changes in version 0.99.30                       

- Added thresholding of number of top genes to use in testing steps.
This avoids high number of false positives in ultra-high dimensional
datasets, e.g. 10x barnyard data.

                       Changes in version 0.99.29                       

- Updated functions of calculating Pearson correlation in sparse
matrix. New functions are more accurate.

                       Changes in version 0.99.25                       

- Changed R dependency back to 4.0.0.

                       Changes in version 0.99.24                       

- Removed direct download files in vignettes to speed up building.
- Changed R dependencies so that package can be built in older R
versions.
- Changed vignette title to match the citation.

                       Changes in version 0.99.23                       

- Minor edit of README.md

                       Changes in version 0.99.22                       

- Received comments from Bioconductor reviewer.
- Moved BiocStyle from Imports: to Suggests:.
- Modified vignettes. Deleted GitHub installation instruction.
- Changed default value of Ncores to be 2.
- Changed output of CB2FindCell to be an SummarizedExperiment object.
Changed GetCellMat accordingly. Changed examples and testthat
accordingly.

                       Changes in version 0.99.21                       

- Updated testthat scripts.
- Minor bug fix in Read10xRaw.R. Added as.numeric when building
sparse
matrix.

                       Changes in version 0.99.20                       

- Updated citation.

                       Changes in version 0.99.19                       

- Minor edits on package vignettes.

                       Changes in version 0.99.18                       

- Minor bug fix.
- Minor edits on package vignettes.

                       Changes in version 0.99.17                       

- Initial version preparing to submit to Bioconductor.

                       Changes in version 0.99.16                       

- Change function name Read10X to be Read10xRaw, Read10Xh5 to be
Read10xRawH5.
- Change parameter name PrintProg to be verbose.
- Minor edits on parameter description of function GetCellMat.
- Added a quick function QuickCB2 by combining all necessary
functions
into one. Input: directory of raw data. Output: filtered matrix, or
a Seurat object containing filtered matrix.

                       Changes in version 0.99.15                       

- When dividing barcodes into groups with similar barcode counts, the
last group will be combined with second last group if the number of
barcodes in the last group is less than half of that in the second
last group.

                       Changes in version 0.99.14                       

- Change parameter names to match those in the paper: retain ->
upper.
background -> lower.

                       Changes in version 0.99.13                       

- Rounded baseline threshold.
- Minor changes on vignette.
- Bug fix when retain threshold is larger than any barcode.

                       Changes in version 0.99.12                       

- Minor changes on vignette.
- Baseline (2 * background) clustering threshold will be printed out.

                       Changes in version 0.99.11                       

- Customized knee point calculation function.
- CB2FindCell will print out retain threshold in standard output.
- Minor bug fixes when there is no cluster.

                       Changes in version 0.99.0                        

- Created initial github page.
- Created package.
- Created vignettes.
- Checked package... passed R CMD check and BiocCheck.
- Ready to submit to Bioconductor.

[scDataviz](https://bioconductor.org/packages/scDataviz)
---------

                        Changes in version 1.0.0                        

- Submitted to Bioconductor

[scDblFinder](https://bioconductor.org/packages/scDblFinder)
-----------

                 Changes in version 1.3.25 (2020-10-26)                 

- scDblFinder has important improvements on speed, robustness and
  accuracy

- in additional to doublet calls, scDblFinder reports the putative
  origin (combination of clusters) of doublets

                 Changes in version 1.3.19 (2020-08-06)                 

- scDblFinder now hosts the doublet detection methods formerly part of
  `scran`

[scDD](https://bioconductor.org/packages/scDD)
----

                 Changes in version 1.13.1 (2020-08-17)                 

- When testing for difference in proportion of zeroes, if the Bayesian
  logistic
  regression fails to converge for a given gene (for example, if there
  is perfect
  separation between conditions), then NA will be returned (instead of
  throwing
  an error.

[SCFA](https://bioconductor.org/packages/SCFA)
----

                 Changes in version 0.99.0 (2020-08-13)                 

- Submitted to Bioconductor

[scp](https://bioconductor.org/packages/scp)
---

                        Changes in version 0.99                         

scp 0.99.4

- Update installation instructions <2020-10-14 Wed>

scp 0.99.3

- fix: solved 'invalid subsetting' issue <2020-10-14 Wed>
- Adapted the vignette to remove warnings and fix missing PCA plot.
<2020-10-14 Wed>
- README.md: extended the installation guide, providing both a stable
and a devel installation. <2020-10-13 Tue>
- Removed the LazyLoad from the DESCRIPTION file and adapted the data
loading (eg data(scp1) to data("scp1")) <2020-10-13 Tue>
- Documentation: added data collection description for the 3 example
datasets <2020-10-13 Tue>

scp 0.99.2

- fix: computeFDR can handle missing values (see issue #12)
<2020-10-02 Fri>

scp 0.99.1

- Maintainer subscribed to bioc-devel mailing list
- Removed infIsNA, the implementation was moved to the QFeatures
packages

scp 0.99.0

- Bioconductor submission

[scPCA](https://bioconductor.org/packages/scPCA)
-----

                 Changes in version 1.3.10 (2020-10-16)                 

- Implementing suggested improvements from Aaron Lun

                 Changes in version 1.3.9 (2020-10-12)                  

- scPCA now accepts DelayedMatrix objects as target and background
datasets.

                 Changes in version 1.3.8 (2020-09-01)                  

- Minor bug fixes

                 Changes in version 1.3.6 (2020-08-30)                  

- Fixed issue where n_centers was required when only one penalty and
contrast term were provided
- Users can now pass factors and character vectors to the clusters
argument.

                 Changes in version 1.3.5 (2020-08-18)                  

- Fixed citations in docs
- Provided more detailed warning when RSpectra::eigs_sym fails to
converge
- Included argumetns in scPCA to control RSpectra::eigs_sym
convergence: error tolerance and max number of iterations

                 Changes in version 1.3.4 (2020-08-12)                  

- Replaced calls to base::eigen by RSpectra::eigs_sym to speed up
eigendecompositions of contrastive covariance matrices. cPCA is now
performed much more quickly when only whishing to compute a handful
of leading contrastive principal components.
- Replaced calls to stats::cov by coop::covar to speed up computation
of large sample covariance matrices.
- In future updates, we'd like to explore using the DelayedArray
framework to support the analysis of larger datasets.

                 Changes in version 1.3.3 (2020-08-08)                  

- The n_centers argument no longer matters when When the contrasts
argument is of length 1 and the penalty term is set to 0.
- Users can now pass in their own cluster labels

                 Changes in version 1.3.2 (2020-08-05)                  

- Updated scPCA function documentation
- Corrected pelling mistakes

[scran](https://bioconductor.org/packages/scran)
-----

                       Changes in version 1.18.0                        

- Deprecated coassignProbs() as this is replaced by
  bluster::pairwiseRand()

- Deprecated boostrapCluster() as this is replaced by
  bluster::bootstrapStability().

- Deprecated gene.names= in the various pairwise* functions as
  being out of scope.

- Added the testLinearModel() function to obtain inferences from
  a linear model.

- Modified pseudoBulkDGE() to use formulas/functions in the
  design= argument. Allow contrast= to be a character vector to
  be run through makeContrasts().

- Added the pseudoBulkSpecific() function to test for
  semi-label-specific DEGs in pseudo-bulk analyses.

- Added the summaryMarkerStats() function to compute some basic
  summary statistics for marker filtering.

- Modified row.data= in findMarkers() to support list inputs.
  Added a add.summary= option to easily include summary
  information.

- Modified combineVar() and combineCV2() to support list inputs.

- Deprecated doubletCells() as this is replaced by
  scDblFinder::computeDoubletDensity().

- Deprecated doubletCluster() as this is replaced by
  scDblFinder::findDoubletClusters().

- Deprecated doubletRecovery() as this is replaced by
  scDblFinder::recoverDoublets().

- Added sparse-optimized variance calculations to modelGeneVar(),
  modelGeneCV2() and related functions, which may result in
  slight changes to the results due to numeric precision.

- Exported combineBlocks() to assist combining of block-wise
  statistics in other packages.

- Added lowess= and density.weights= options to fitTrendVar() to
  rescue overfitted curves.

- Raised an error in denoisePCA() upon mismatches in the matrix
  and technical statistics.

[scRepertoire](https://bioconductor.org/packages/scRepertoire)
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

                       Changes in version 0.99.18                       

- Updated author information in the vignette

                       Changes in version 0.99.17                       

- Updated NEWS formatting

- Edited DESCRIPTION to SingleCellExperiment R package

- Updated information in the vignette

                       Changes in version 0.99.16                       

- Added getCirclize()

                       Changes in version 0.99.15                       

- Modified numerator for index function

                       Changes in version 0.99.14                       

- Removed bracket from indexing function

                       Changes in version 0.99.13                       

- Added exportTable to remaining viz functions

- Modified morisita index to correct error

                       Changes in version 0.99.12                       

- Reducing the size of the screp_example to fulfill < 5 mB requirement.
  Randomly samples 100 cells and removed RNA counts from Seurat object

                       Changes in version 0.99.11                       

- Updated compareClonotype to allow for clonotype comparisons

                       Changes in version 0.99.10                       

- Bioconductor did not detect the version update.

                       Changes in version 0.99.9                        

- Bioconductor had no love - changed the Seurat package to imports
  instead of required, see if that will address the compiling issue
  that results in a killed: 9 error.

                       Changes in version 0.99.8                        

- Passed checks on system, let's see how much bioconductor hates it

                       Changes in version 0.99.7                        

- But really this time, changed the colData import

                       Changes in version 0.99.6                        

- Changed colData import

                       Changes in version 0.99.5                        

- Added screp_example data to package

- Added visVgene function for visualizing the distribution of V genes
  in TCR

- Added support for monocle to combineExpression function

- Updated documentation for combineTCR() and combineBCR()

- Updated documentation to utilize SingleCellExperiment formats

- Updated Vignette to utilize SingleCellExperiment formats

- Added Author information to vignette

- Add intro and conclusion to vignette

- Removed html knitted vignette

- Removed descriptive code snippets

                       Changes in version 0.99.4                        

- Modified expression2List() to allow for variables across meta data

                       Changes in version 0.99.1                        

- Changed R (>= 3.6) to R (>= 4.0)

                       Changes in version 0.99.0                        

- Changed DESCRIPTION version to 0.99.0

- Removed file seurat_example.rda, accidentally committed

- Deleted git attributes

- reduced Seurat object size for alluvialClonotype in vignette

- Changed the alluvialClonotype assessment to account for only 1
  condition

[scTensor](https://bioconductor.org/packages/scTensor)
--------

                        Changes in version 2.0.0                        

- Extended to use the version 2.0.0 of LRBase.XXX.eg.db-type packages

- Omitted typing Enter-key many times to perform
  example('cellCellReport')

- lr.evidence option was added in cellCellRanks() and cellCellDecomp()
  to select ligand-receptor databases to construct CCI-tensor (cf.
  Evidence code: https://github.com/rikenbit/lrbase-workflow)

- The L-R evidence was embeded in the HTML report.

- The bug related to the hyper-link was fixed in .hyperLinks.

- Auto library installation/loading for MeSH.XXX.eg.db. Note that this
  function is based on the NCBI Taxonomy ID embedded in METATDATA table
  in the sqlite3 file of LRBase.XXX.eg.db.

- The order of the parameters of cellCellSetting were changed as
  cellCellSetting(sce, lrbase, label, lr.evidence="all", color=NULL),
  and color parameter was changed as a optional parameter. If it is not
  specified, the color is automatically selected.

- The way to specify the celltype label was changed in cellCellSetting.

- The vignettes were modified.

- The convertNCBIGeneID is deprecated. The same functionality can be
  available as scTGIF::convertRowID instead.

                       Changes in version 1.4.1-3                       

- A bug was fixed in .cellCellDecomp.Halpern()

[scTGIF](https://bioconductor.org/packages/scTGIF)
------

                        Changes in version 1.2.1                        

- A bug is fixed in convertRowID()

[scuttle](https://bioconductor.org/packages/scuttle)
-------

                        Changes in version 1.0.0                        

- Split off scuttle from scater by migrating all
  non-visualization code from the latter.

- Began transition to dot-separate argument names from original
  snake case format.

- Added a geometricSizeFactors() function, deprecated
  geometric=TRUE in librarySizeFactors().

- Single-object downsampling in downsampleBatches() now behaves
  more consistently with multi-object downsampling.

[SDAMS](https://bioconductor.org/packages/SDAMS)
-----

                        Changes in version 1.9.2                        

- update description for single-cell data in reference manual.

                        Changes in version 1.9.1                        

- add example for single-cell data.

[semisup](https://bioconductor.org/packages/semisup)
-------

                 Changes in version 1.14.0 (2020-05-05)                 

- Updated documentation

[SeqArray](https://bioconductor.org/packages/SeqArray)
--------

                       Changes in version 1.29.2                        

UTILITIES

- show a warning when an unsorted index is used in `seqSetFilter()`

- show a message if `seqVCF_Header()` fails

- a new option 'chr_prefix' in `seqGDS2VCF()`

BUG FIXES

- `seqVCF_Header()` fixes 'contig' in the header of VCF if there are
  different fields

                       Changes in version 1.28.1                        

BUG FIXES

- `seqRecompress(, verbose=FALSE)` works correctly

- `seqSetFilter(, action="push+set")` should not reset the filter
  before setting a new filter

[SeqGate](https://bioconductor.org/packages/SeqGate)
-------

                 Changes in version 0.99.3 (2020-09-18)                 

- Updated vignette and documentation regarding SummarizedExperiment

- Updated NAMESPACE regarding SummarizedExperiment

                 Changes in version 0.99.2 (2020-09-08)                 

- Changed code to take SummarizedExperiment objects as input and output

- Changed tests for this new code

                 Changes in version 0.99.1 (2020-09-01)                 

- Added BugReports field to the DESCRIPTION file

- Fixed typo in checkForReplicates function name

- Replaced cat by message, warning or stop

- Suppressed unappropriate keywords (e.g., ERROR, WARNING) in messages

- Replaced external datasets by generated datasets in tests

- Suppressed rm(list=list()) in tests

                 Changes in version 0.99.0 (2020-07-09)                 

- Submitted to Bioconductor

[seqLogo](https://bioconductor.org/packages/seqLogo)
-------

                       Changes in version 1.55.2                        

NEW FEATURES

- Added accessor function to access pwm, ic and consensus

- Added possibility to specify colors (by setting fill)

- Added support for RNA logos

BUG FIXES

- Replaced `class` statements with `is(..., "class")`

[SeqVarTools](https://bioconductor.org/packages/SeqVarTools)
-----------

                       Changes in version 1.27.1                        

- Add 'parallel' argument to methods.

[sevenbridges](https://bioconductor.org/packages/sevenbridges)
------------

                       Changes in version 1.19.2                        

Improvements

- Transfer package maintainership.

                       Changes in version 1.19.1                        

Bug Fixes

- Fixed an API issue on invalid JSON when running tasks.
- Fixed a subtle column class mismatch issue when binding data frames
since R 4.0.0.

Improvements

- Use data.table::rbindlist() when possible to increase the data
frame
binding performance.

[SharedObject](https://bioconductor.org/packages/SharedObject)
------------

                       Changes in version 1.3.18                        

- Make the package functions process-safe

                        Changes in version 1.3.8                        

- Support complex

- Support pairlist

- Use memcpy to copy the data when possible

                        Changes in version 1.3.7                        

- change the parameters of is.shared function

- `mustWork = FALSE` by default, no error will be given
  when sharing a nonsharable object

- Support environment object

                        Changes in version 1.3.6                        

- Support S4 object

- new dispatching method(signiture "ANY") for share, unshare and
  is.shared functions

[ShortRead](https://bioconductor.org/packages/ShortRead)
---------

                        Changes in version 1.48                         

NEW FEATURES

- (v 1.47.1) `countFastq()` counts the number of records,
  nucleotides, and quality scores in one or several fastq files.

[signatureSearch](https://bioconductor.org/packages/signatureSearch)
---------------

                 Changes in version 1.2.4 (2020-08-14)                  

- Supported defining gene set database from score matrix by setting
higher, lower, as well as padj cutoffs for gCMAP and Fisher GESS
methods

                 Changes in version 1.2.2 (2020-07-11)                  

- Supported converting gmt file to HDF5 file (01 matrix) as gene set
reference database for gCMAP and Fisher GESS methods

[simplifyEnrichment](https://bioconductor.org/packages/simplifyEnrichment)
------------------

                       Changes in version 0.99.5                        

- in the clustering, considered when the size of cluster is only 1.

- add `partition_by_kmeanspp()`.

- `col` can be set as a vector of colors.

                       Changes in version 0.99.4                        

- support more ontologies.

- support to calculate similarity matrix by gene overlap.

- node partition fun changed to pam

[SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment)
--------------------

                       Changes in version 1.12.0                        

- Added the rowSubset() function as a standard location for a row
  subset.

- Added colPairs() and rowPairs() to store pairwise information
  (e.g., for graphs).

- Added method specifications for S4Vectors compatibility.

[singleCellTK](https://bioconductor.org/packages/singleCellTK)
------------

                 Changes in version 2.0.0 (2020-10-16)                  

- Added quality control (empty droplet detection, doublet detection,
  etc) functionality

- Ability to import data from varying preprocessing tools

- Ability to export SingleCellExperiment object as varying file types
  (flat file, Python anndata)

- Added functions for visualization of data

- New CellViewer functionality in UI

- Improvements to differential expression, now includes DESeq2, limma,
  ANOVA

- Incorporates Seurat workflow

[SingleR](https://bioconductor.org/packages/SingleR)
-------

                        Changes in version 1.4.0                        

- Migrated all of the dataset getter functions to the celldex
  package.

- Streamlined the vignette to point to the book at <URL:
  http://bioconductor.org/books/devel/SingleRBook/>.

- Added a restrict= argument to trainSingleR() and SingleR() to
  easily restrict to a subset of features.

- Deprecated the method= argument in SingleR().

- Protect against accidental data.frames in ref= or test= in all
  functions.

[sitePath](https://bioconductor.org/packages/sitePath)
--------

                       Changes in version 1.5.25                        

- Bug fix: use 'geom_point2' instead of 'geom_tippoint' to avoid error.

                       Changes in version 1.5.24                        

- Add sequence type option for for DNA and amino acid.

- Deprecate 'multiFixationSites' function.

- Use 'y' argument as mutation label option in 'plot.sitePath'
  function.

- Finer lineage resolving method used in 'lineagePath' function.

                       Changes in version 1.5.23                        

- Create 'groupTips' functions to replace 'as.list' functions for
  'fixationSites' and 'fixationPath' object.

- Create 'sitesMinEntropy' function to output raw result of entropy
  minimization.

- Create 'parallelSites' function and other functions for its return
  object
  such as 'plotSingleSite' and 'as.data.frame'

                       Changes in version 1.5.22                        

- Fix missing newline when printing 'phyMSAmatched' object.

- Create 'as.list.fixationSites' for retrieving grouped tips.

- Remove 'tipname' option in 'as.data.frame.fixationSites'.

                       Changes in version 1.5.21                        

- Fix wrong group name in some corner cases.

- Use 'ggtree' for 'plotSingleSite'.

                       Changes in version 1.5.20                        

- Further fix the merging issue in 'fixationSites'.

                       Changes in version 1.5.19                        

- Speed up 'SNPsites'.

                       Changes in version 1.5.18                        

- Bifurcation check for the phylogenetic tree and force bifurcation.

- Fix path merging issue in 'fixationSites'.

                       Changes in version 1.5.17                        

- Import 'aes' and 'theme' from 'ggplot2'.

                       Changes in version 1.5.16                        

- Allow turning off mutation label for 'plot.fixationSites' while
  legend of cluster name becomes compulsory.

- Import 'scale_color_manual' from 'ggplot2'.

- Update vignette.

                       Changes in version 1.5.15                        

- Bug fix: NA in cluster name.

                       Changes in version 1.5.14                        

- Add 'as.treedata' function for 'fixationSites'.

                       Changes in version 1.5.13                        

- Hierarchical naming of the clusters.

                       Changes in version 1.5.12                        

- Establish 'phyMSAmatched' S3 class for better encapsulation.

                       Changes in version 1.5.11                        

- Deprecate 'multiFixationSites' function.

                        Changes in version 1.5.8                        

- Wrap mutations text in 'plot.fixationSites'.

- Remove 'color' argument for 'plot.fixationSites' as the number of
  groups
  is usually unknown.

                        Changes in version 1.5.7                        

- Left padding with 0 for the cluster name.

- Add mutation label when plot 'fixationSites'.

- Add 'as.data.frame' function for 'fixationSites'.

                        Changes in version 1.5.6                        

- Add 'sitewiseClusters' function and plot function for its
  visualization.

                        Changes in version 1.5.5                        

- Add 'plotMutSites' function to visualize mutations of each tree tip.

                        Changes in version 1.5.4                        

- Use 'ggtree' for 'plot.lineagePath'.

- More informative plot for 'sneakPeek' and add 'lineagePath' function
  for its return.

                        Changes in version 1.5.3                        

- Add 'as.phylo.fixationSites' function that represent site fixations
  as simplified phylgenetic tree.

                        Changes in version 1.5.2                        

- Add 'minEffectiveSize' in 'plot.fixationSites' for filtering
  small sized tip clusters.

                        Changes in version 1.5.1                        

- Add 'plot.fixationSites' function.

                        Changes in version 1.4.1                        

- Fix: broken link in the DESCRIPTION.

[snifter](https://bioconductor.org/packages/snifter)
-------

                 Changes in version 0.99.0 (2020-07-17)                 

- First version of the package.

[SNPRelate](https://bioconductor.org/packages/SNPRelate)
---------

                       Changes in version 1.24.0                        

- definition of IBS in the `snpgdsIBS` help file

[sparseMatrixStats](https://bioconductor.org/packages/sparseMatrixStats)
-----------------

                 Changes in version 1.0.5 (2020-05-17)                  

- Fix links in documentation to get rid of WARNINGS

                 Changes in version 1.0.4 (2020-05-17)                  

- Fix bugs in colTabulates

- Update documentation to avoid warnings in build on Windows

                 Changes in version 1.0.3 (2020-05-17)                  

- Fix bug in colAnys and colAlls if value = TRUE

                 Changes in version 1.0.2 (2020-05-10)                  

- Fix bug in colMaxs, colMins related to missing values

                 Changes in version 1.0.1 (2020-05-08)                  

- Fix bugs in colMaxs, colMins, colAnys

- Fix bug in colLogSumExps

[SparseSignatures](https://bioconductor.org/packages/SparseSignatures)
----------------

                        Changes in version 2.0.0                        

- Migration from Travis-CI to Github Actions.

- Major refactoring.

[SpatialDecon](https://bioconductor.org/packages/SpatialDecon)
------------

                 Changes in version 0.99.0 (2020-10-02)                 

- Submitted to Bioconductor

[spatialHeatmap](https://bioconductor.org/packages/spatialHeatmap)
--------------

                 Changes in version 0.99.0 (2020-09-21)                 

- Submitted the first version to Bioconductor

[Spectra](https://bioconductor.org/packages/Spectra)
-------

                        Changes in version 0.99                         

Changes in 0.99.11

- Re-add mz and intensity as core spectra variables.

Changes in 0.99.10

- Fix in spectraData<-,Spectra to avoid removing m/z and intensity
values (issue #146).
- Add default implementations of filter functions for MsBackend.

Changes in 0.99.9

- Fix in Spectra,character constructor to ensure the backend is
changed even if source inherits from backend (issue #143).

Changes in 0.99.8

- combineSpectra applies data processing steps in the processing
queue
prior to combination (issue #140).

Changes in 0.99.7

- Fix problem in dropNaSpectraVariables that would also drop m/z and
intensity values for most backends (issue #138.

Changes in 0.99.6

- Support intensity in filterIntensity method to be a function to
enable peak intensity-based filtering of spectra (issue #126).

Changes in 0.99.5

- Add filterMzRange and filterMzValues to filter spectra based on an
m/z range or a list of target m/z values, respectively.

Changes in 0.99.4

- Add export,MsBackendMzR to export spectra data to mzML or mzXML
file(s).
- Add an export,MsBackend method to allow backends to take care of
data export.
- Refactor export,Spectra to use the MsBackend class to export the
data.
- Change parameter source in Spectra,character to MsBackendMzR and
set
parameter backend = source. Thus by default, the import backend will
also be used to store the data.

Changes in 0.99.3

- Replace lapply,Spectra with spectrapply,Spectra.

Changes in 0.99.2

- Replace asDataFrame,MsBackend with spectraData,MsBackend.
- Replace asDataFrame<-,MsBackend with spectraData<-,MsBackend.
- Replace as.list,MsBackend with peaksData,MsBackend.
- Replace replaceList<-,MsBackend with peaksData<-,MsBackend.
- Replace as.list,Spectra with peaksData,Spectra and add methods to
coerce a Spectra to a list or SimpleList.

Changes in 0.99.0

- Add reset method.
- Add processing by chunk to compareSpectra.

[splatter](https://bioconductor.org/packages/splatter)
--------

                 Changes in version 1.14.0 (2020-10-28)                 

- Add the splatPop simulation. This is a extension to the splat
  simulation contributed by Christina Azodi and Davis McCarthy
  that adds population effects. It allows you to specify
  relatedness between individuals and generate cell-type specific
  eQTL effects.

- Add a batch.rmEffect parameter to the Splat simulation. This
  allows generation of a paired simulation without any batch
  effects.

- Add a new minimiseSCE function which can be used to remove
  unneeded information from simulation output (or any
  SingleCellExperiment)

- All simulations now return sparse assay matrices by default
  when they would be smaller than the equivalent dense matrix.
  This is controlled by a new sparsify argument.

- Users will now be automatically prompted to install packages if
  they try to use a simulation for which the suggested
  dependencies are not available

[splineTimeR](https://bioconductor.org/packages/splineTimeR)
-----------

                 Changes in version 1.17.1 (2020-07-18)                 

- Fixed bug. Limit for counter of inner loop in function splinePlot
  corrected

- Signigicant changes
  o Log output to console in function splinePlot
  o package NEWS file added

[SPsimSeq](https://bioconductor.org/packages/SPsimSeq)
--------

                       Changes in version 0.99.0                        

- First private GitHub document
- Major revision including gene/featurewise dependence
- Minor revision on documentation and adding new features
- Optimizing outputs, further enhancement, fixing bugs, and
addressing
issues
- Added unit tests

[statTarget](https://bioconductor.org/packages/statTarget)
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

                       Changes in version 1.19.1                        

- To fixed bugs for coCV function

[struct](https://bioconductor.org/packages/struct)
------

                        Changes in version 1.1.2                        

- improved 'show' output for objects

- allow ANY for entities

                        Changes in version 1.1.1                        

- added citations slot to struct classes

- added corresponding citations method

- added method to get/set seq_in slot

- as.SummarizedExepriment now works correctly

- using seq_in now works for sequences with more than 2 steps

[structToolbox](https://bioconductor.org/packages/structToolbox)
-------------

                        Changes in version 1.1.2                        

- Documentation updates

- Bug fixes

                        Changes in version 1.0.1                        

- Fix HCA bug

[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.20.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- SummarizedExperiment now depends on the MatrixGenerics package.

- DelayedArray was moved from Depends to Imports.

DEPRECATED AND DEFUNCT

- Deprecated readKallisto().

BUG FIXES

- Avoid triggering copies of the assays in assays() getter.

- Fix long-standing bug in dim() method for Assays objects.

- Fix assays(x) <- SimpleList(). Before that fix this operation was
  turning
  SummarizedExperiment object (or derivative) 'x' into an invalid
  object.

[SWATH2stats](https://bioconductor.org/packages/SWATH2stats)
-----------

                       Changes in version 1.19.3                        

UPDATE

- Add link to website

                       Changes in version 1.19.2                        

BUG FIXES

- remove grid from Description

                       Changes in version 1.19.1                        

NEW FEATURES

- validate_columns

update

- update documentation to roxygen2 (many thanks to Ashton Trey Belew
  (abelew) for starting this and doing most of the work) (IS#1, PR#2,
  PR#3)

- Some parameters were renamed to omit having a dot. E.g. rm.decoy is
  now called rm_decoy. If you want the previous names use SWATH2stats
  up to version 1.19.0

                       Changes in version 1.19.0                        

NEW FEATURES

- SWATH2stats in BioC 3.12 development release

[TADCompare](https://bioconductor.org/packages/TADCompare)
----------

                Changes in version 0.99.35 (2020-05-15)                 

- Proper formatting of the NEWS file

                Changes in version 0.99.33 (2020-04-27)                 

- Submitted to Bioconductor

[TargetSearch](https://bioconductor.org/packages/TargetSearch)
------------

                       Changes in version 1.46.0                        

NEW FEATURES

- New function `ri_data_extract` to extract peaks from RI files. It
  works
  similar to `FindAllPeaks` but uses different (simpler) input
  parameters,
  comparable to ncdf4_data_extract.

- New function `ri_plot_peak` to plot peaks from RI files, built upon
  `ri_data_extract`. It can be used as an alternative to `plotPeakRI`
  as it
  has a simple interface.

- New function `ncdf4_plot_peak`. An alternative to function
  `plotPeakSimple`
  with a simple interface to plot peaks from NetCDF format 4. This
  function
  supersedes `plotPeakSimple`.

BUG FIXES

- Remove unneeded ICO file.

- Man pages improvements. Mostly grammar and spelling changes.

SIGNIFICANT USER-VISIBLE CHANGES

- The function `peakPlotSimple` is considered deprecated, and its
  use should be avoided. Use the function `ncdf4_plot_peak` instead.

- The parameter `column` in many columns is now `NULL` by default. To
  change the column names use the global option `TS_RI_columns` instead

[TCC](https://bioconductor.org/packages/TCC)
---

                       Changes in version 1.29.1                        

- delete 'deseq' option, since DESeq is deprecated and will be removed
  from Bioconductor

[TCGAutils](https://bioconductor.org/packages/TCGAutils)
---------

                       Changes in version 1.10.0                        

New features

- correctBuild attempts to provide the official name of a particular
human genome build to agree with changes in GenomeInfoDb
- isCorrect checks that the build name matches the official name

Minor changes and bug fixes

- Documentation improvements to simplifyTCGA
- Improvements to findGRangesCols to locate ranged columns in a
DataFrame
- Fixed a bug in UUIDtoBarcode where only the first record was
returned (#26, @DarioS)
- Fixed a bug in filenameToBarcode when multiple inputs were used
(#22, @DarioS)

[tofsims](https://bioconductor.org/packages/tofsims)
-------

                        Changes in version 099.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- changed function behvaiour in the whole
  package from call-by-ref to call-by value. Adjusted
  accordingly all examples and the vignette.

INTERNALS

- depends now on ProtGenerics from which it uses 'mz'

- exchanged various print() with message()

[tomoda](https://bioconductor.org/packages/tomoda)
------

                 Changes in version 0.99.0 (2020-09-15)                 

- Submitted to Bioconductor

[topconfects](https://bioconductor.org/packages/topconfects)
-----------

                        Changes in version 1.5.3                        

- Small tweaks to better display weitrix_confect outputs.

- Less digits when print()ed.

                        Changes in version 1.5.2                        

- nest_confects now copes with NA p-values.

                        Changes in version 1.5.1                        

- Fix limits calculation for confects_plot.

[ToxicoGx](https://bioconductor.org/packages/ToxicoGx)
--------

                        Changes in version 1.1.0                        

- Continue to abstract functionality into CoreGx
- Add additional plotting functions such as grouped boxplots
- Extend coverage of unit tests to >90%
- Implement a faster version of drugPertubationSignature
- Add additional plotting functions
- Include scripts for differential expression analysis and GSEA of
toxico-genomic data (limma)

                        Changes in version 1.0.0                        

- Package archived on CRAN
- Package submitted to Biocondcutor
- Modified package to depend on updated CoreGx
- All molecularProfiles are now SummarizedExperiment instead of
ExpressionSet
- Abstracted some additional functions to CoreGx

                        Changes in version 0.1.2                        

- Updated downloadTSet function to use published Zenodo DOIs to
retrieve data
- Modified rankGeneDrugsPerturbation to fix a bad unit conversion
which would return concentrations in the wrong unit

                        Changes in version 0.1.1                        

- Bug Fix: Regenerated TGGATESsmall (sample dataset) to fix make a
result in the vignette consistent with previous releases.

                        Changes in version 0.1.0                        

- Rewrote plots using ggplot2 to improve aesthetics
- Also can now extend plotting functions using standard ggplot2
syntax
- Improved package documentation

                        Changes in version 0.0.1                        

- Minimal package submitted

[TPP](https://bioconductor.org/packages/TPP)
---

                       Changes in version 3.17.6                        

- fix warnings due to unused argument of select(!!!syms(...)) statement
  during histogram generation for
  reference data

                       Changes in version 3.17.5                        

- fix bugs and warnings in executable examples

                       Changes in version 3.17.4                        

- fix errors in build report on Bioconductor

- fix warnings due to conversion of non-numeric values during 2D-TPP
  import

- fix warnings in ggplot command

- upgrade deprecated dplyr functions

- make syntax of testthat checks more consistent across files

                      Changes in version 3.17.2-3                       

- Fixed bugs and warnings after update to dplyr v1.0.0.

                       Changes in version 3.17.1                        

- Fixed bug in tpp2dCreateTPPTRreferenece upon user request (#10)

- removed adding spline fit column from 2DTPP output table

- fixed bug in tpp2dCreateTPPTRReference and made example in vignette
  work

- removed leftover parameters in tpp2dCreateTPref function

- final fix of tpp2dCreateTPPTRreference and call in vignette now works

                       Changes in version 3.17.0                        

- New Bioconductor release candidate

[TPP2D](https://bioconductor.org/packages/TPP2D)
-----

                        Changes in version 1.5.5                        

- new visualization options

- volcanoplot: `plot2dTppVolcano`

- heatmap of fold changes: `plot2dTppFcHeatmap`

                 Changes in version 1.4.1 (2020-06-20)                  

- bug fix in moderated F statistic computation - old version is
  unneccessarily stringent

[trackViewer](https://bioconductor.org/packages/trackViewer)
-----------

                       Changes in version 1.25.4                        

- remove http_status from documentation.

- add function to split the lollipop plot into multiLayers.

                       Changes in version 1.25.3                        

- Fix the issue in AddArrowMark that grid changed the unit id.

                       Changes in version 1.25.2                        

- Fix the issue if there is interrupt of the internet connection to
  generate vignette.

                       Changes in version 1.25.1                        

- Provide more clearly warning or error message if the input is not
  sorted or contain NA values.

[tradeSeq](https://bioconductor.org/packages/tradeSeq)
--------

                 Changes in version 1.3.19 (2020-10-16)                 

- Merged `conditions` branch into `master`. Now the `master` branch can
  therefore handle multiple conditions. There is a new `conditions`
  argument in `fitGAM` to handle that, and a condition-specific
  smoother will be fitted for each lineage. There is also a new test,
  `conditionTest`, which tests for DE between conditions within a
  lineage. The way this is done under the hood is exactly like the
  `patternTest`.

[transite](https://bioconductor.org/packages/transite)
--------

                        Changes in version 1.6.3                        

- added limits argument for spectrum plot x-values to harmonize color
  scale
  of multiple spectrum plots

- replaced Rcpp::RcppArmadillo::sample() with now fixed Rcpp::sample()

                        Changes in version 1.6.2                        

- updated paper references

- fixed minor formatting issues

- incremented Roxygen2 version

[transomics2cytoscape](https://bioconductor.org/packages/transomics2cytoscape)
--------------------

                        Changes in version 1.0.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Package introduced.

NEW FEATURES

- Package introduced.

[treeio](https://bioconductor.org/packages/treeio)
------

                       Changes in version 1.13.1                        

- as_tibble for pvclust (2020-06-22, Mon)
- as.phylo and as.treedata methods for pvclust object (2020-06-21,
Sun)

[TRONCO](https://bioconductor.org/packages/TRONCO)
------

                       Changes in version 2.21.1                        

- From Travis-CI to Github Actions

[TSCAN](https://bioconductor.org/packages/TSCAN)
-----

                       Changes in version 1.28.0                        

- Added createClusterMST() to create a cluster-based MST from a
  variety of inputs, migrated from the scran package.

- Added reportEdges() to report edge coordinates for plotting.

- Added mapCellsToEdges() to map cells to the closest edge on the
  MST.

- Added orderCells() to compute a pseudotemporal ordering from
  mapped cells.

- Added quickPseudotime() to wrap MST construction and ordering
  into a single call.

- Added testPseudotime() to test for DE genes along one or more
  paths through a MST.

- Added the rowmean() utility to compute column means for row
  groupings.

- Added perCellEntropy() to compute per-cell entropies across
  various matrix types.

[tximeta](https://bioconductor.org/packages/tximeta)
-------

                        Changes in version 1.8.0                        

- Added 'fromDb' argument to addIds() to allow IDs to be
  added from the associated TxDb/EnsDb instead of the org
  package (which is used by default).
  Feature suggestion from Kristoffer Vitting-Seerup.

- Added function retrieveCDNA() that will download or load
  a cached version of the transcript sequences used for
  quantification. Note that the returned sequences are not
  ordered or matched to the rows of the SummarizedExperiment
  object. Feature suggestion from Kristoffer Vitting-Seerup.

- Added function addCDS() that will add CDS ranges for coding
  transcripts (and fills in original ranges for non-coding),
  as well as a 'coding' column as a logical indicator.
  Feature suggestion from Kristoffer Vitting-Seerup.

- Added option that environmental variable TXIMETA_HUB_CACHE
  can be used to set tximeta's cache location, to avoid
  prompting the user on the first run of tximeta().

- tximeta() will now pull GENCODE TxDb from AnnotationHub
  when it is listed there (only Homo sapiens are at this
  point in time). Thanks to Leonardo Collado-Torres for
  the suggestion!

- Now summarizeToGene() will add a column tx_ids, which is a
  CharacterList of the transcript IDs.

                       Changes in version 1.7.13                        

- Added 'fromDb' argument to addIds() to allow IDs to be
  added from the associated TxDb/EnsDb instead of the org
  package (which is used by default).
  Feature suggestion from Kristoffer Vitting-Seerup.

                       Changes in version 1.7.12                        

- Added function retrieveCDNA() that will download or load
  a cached version of the transcript sequences used for
  quantification. Note that the returned sequences are not
  ordered or matched to the rows of the SummarizedExperiment
  object. Feature suggestion from Kristoffer Vitting-Seerup.

                       Changes in version 1.7.11                        

- Added function addCDS() that will add CDS ranges for coding
  transcripts (and fills in original ranges for non-coding),
  as well as a 'coding' column as a logical indicator.
  Feature suggestion from Kristoffer Vitting-Seerup.

                       Changes in version 1.7.10                        

- Added option that environmental variable TXIMETA_HUB_CACHE
  can be used to set tximeta's cache location, to avoid
  prompting the user on the first run of tximeta().

                        Changes in version 1.7.9                        

- tximeta() will now pull GENCODE TxDb from AnnotationHub
  when it is listed there (only Homo sapiens are at this
  point in time). Thanks to Leonardo Collado-Torres for
  the suggestion!

                        Changes in version 1.7.6                        

- Now summarizeToGene() will add a column tx_ids, which is a
  CharacterList of the transcript IDs.

                        Changes in version 1.7.1                        

- Updated to Ensembl release 100, GENCODE 34/M25.

[tximport](https://bioconductor.org/packages/tximport)
--------

                       Changes in version 1.18.0                        

- Code cleanup for deprecated functions. tximport reads in
  only alevin version 0.14.0 or greater. For older data, use
  previous versions of tximport.

[UMI4Cats](https://bioconductor.org/packages/UMI4Cats)
--------

                       Changes in version 0.99.21                       

- Improvements in vignette and documentation.

                       Changes in version 0.99.20                       

- Avoid duplication of fragment end calculation in plotDifferential for
  DESeq2 results.

                       Changes in version 0.99.19                       

- Minor documentation changes to pass Bioconductor checks.

                       Changes in version 0.99.18                       

- Fix DESeq2 example.

- Reduce size of installed example tsv.gz count files.

                       Changes in version 0.99.17                       

- Re-run example datasets and use again links from figshare.

                       Changes in version 0.99.16                       

IMPORTANT

- This version is not compatible with the .tsv.gz files created by
  previous versions. You will need to run `contactsUMI4C()` again to
  generate
  updated .tsv.gz files.

UPDATES

- Improved grouping arguments for UMI-4C objects: now creates a new
  UMI-4C
  object that can be accessed using `groupsUMI4C(umi4c)$condition`.
  This allows
  retaining replicate information in the main UMI4C object while
  allowing
  plotting grouped trends stored in `groupsUMI4C()`.

- Added new statistical test using DESeq2:
  `differentialNbinomWaldTestUMI4C()`.

                       Changes in version 0.99.15                       

- Fixed duplicated read number in read id (.singlePrepUMI4C) (see issue
  #5).

- Changed example download urls to gattaca server.

                       Changes in version 0.99.14                       

- Fixed bug where adaptive smoothed trend was normalized twice (see
  issue #4).

                       Changes in version 0.99.11                       

- Uploaded example datasets urls in `downloadUMI4CexampleData()` to a
  more
  stable and permanent location (figshare.com).

                       Changes in version 0.99.10                       

- Avoid running long and redundant examples, already tested in the
  vignette to
  avoid TIMEOUT build error.

                       Changes in version 0.99.9                        

- Add data object `ex_ciita_umi4c` to use in examples and reduce check
  running
  times.

                       Changes in version 0.99.8                        

- Update package vignette to clarify the origin of the different sample
  files
  used to exemplify a workflow using the UMI4Cats package.

                       Changes in version 0.99.7                        

- Added unit tests using `testthat`.

- Use BiocFileCache to download sample files.

- Use tempdir() for demo purposes both in vignette and examples.

- Added inst\scripts to describe how the sample data was generated.

- Other minor changes to comply with Bioconductor review
  (see
  https://github.com/Pasquali-lab/UMI4Cats/issues/2#issue-637249954)

                       Changes in version 0.99.6                        

- Delete downloaded and intermediate folders when building vignette.

- Added `UMI4Cats_index` to .Rbuildignore to prevent ERRORs and
  WARNINGs in BioCCheck.

                       Changes in version 0.99.5                        

- Increased speed of `getViewpointCoordinates()` by allowing
  pre-selection
  of viewpoint chromosome using `sel_seqname` argument.

- Added reduced fastq files in extdata and allow downloading of reduced
  bowtie index to increase vignette building speed.

                       Changes in version 0.99.4                        

- Added `.Rproj` files to .gitignore

                       Changes in version 0.99.3                        

- Changed example in vignette and manuals to *CIITA*.

- Added viewpoint name in `plotTrend()`.

- Improved multi-panel plotting of `plotUMI4C()`.

                       Changes in version 0.99.2                        

- Allow `ref_umi4c` to be used as reference for plotting colors,
  domainogram
  and differential analysis (not only for normalization).

- Fixed error when using `sampleID` as `grouping` variable in
  `makeUMI4C()`.

- Fixed bug in `results()` when `fomat=data.frame` and `ordered=TRUE`.

- Improved visualization of differential regions reconverting `Inf` and
  `-Inf`
  to maximum and minimum (respectively) odd's ratio values.

- Add more functionality details in the `Analyzing UMI-4C data with
  UMI4Cats`
  vignette.

                       Changes in version 0.99.1                        

- Fixed error in function `createGeneAnnotation` and `plotGenes` that
  occurs
  when there are no genes in the region or a gene has multiple
  identifiers.

- Fixed duplicated generics definition for `SummarizedExperiment`
  objects to
  avoid error when reloading the package.

- Fixed error when `bait_exclusion` is set to 0.

- Added possibility to specify the sample to use as reference for
  normalization
  (`ref_umi4c` argument in `makeUMI4C`).

- Now the `grouping` variable in `makeUMI4C()` is used more upstream in
  the
  analysis. For using different grouping variables, user must create
  different
  `UMI4C` objects.

- Fixed bug where sometimes bait coordinates in the output tsv file are
  `NA`.

- `statsUMI4C` now also outputs a stats summary table in
  `wk_dir/logs/stats_summary.txt`.

- Improve function documentation.

- Improve pkgdown UMI4Cats site.

- Rewrite and improve the `Analyzing UMI-4C data with UMI4Cats`
  vignette.

                       Changes in version 0.99.0                        

- First public release of UMI4Cats.

- Added a `NEWS.md` file to track changes to the package.

[uncoverappLib](https://bioconductor.org/packages/uncoverappLib)
-------------

                     Changes in version 0.0.0.9000                      

setup

- added NEWS.md creation

[universalmotif](https://bioconductor.org/packages/universalmotif)
--------------

                        Changes in version 1.8.0                        

NEW FEATURES

- scan_sequences()/enrich_motifs() can now be used to scan/enrich for
  gapped motifs. A new section has been added to the
  SequenceSearches.Rmd
  vignette.

- scan_sequences(..., use.gaps), enrich_motifs(..., use.gaps): ignore
  motif
  gap information.

- read_meme(), write_meme(): now fully support custom alphabets.

- prob_match(), prob_match_bkg(): calculate the probability of a motif
  match based on background frequencies of the motif object or provided
  values, respectively.

- enrich_motifs(), get_matches(), get_scores(), motif_pvalue(),
  motif_score(), scan_sequences(), score_match(): new allow.nonfinite
  parameter, allowing for the functions to work even if non-finite
  values
  are present in the motif PWM.

- read_matrix(..., comment): allows for comments to be ignored in motif
  files.

- write_matrix(..., digits): control the number of digits to use for
  writing motif positions.

- New mask_seqs() utility function: inject hard masks into sequences.

- scan_sequences(..., warn.NA), enrich_motifs(..., warn.NA): new option
  which can disable warnings from non-standard letters being detected
  in the input sequences.

- get_bkg(..., window, window.size, window.overlap): new options for
  calculating sequence background in windows.

- get_bkg(..., merge.res): new option to return background information
  for
  individual sequences.

- scan_sequences(..., calc.pvals): new option to calculate P-values for
  sequence hits. This is merely automating using the results from
  scan_sequences() to calculate P-values manually with motif_pvalue().

- view_motifs(..., show.positions, show.positions.once, show.names):
  new
  options for customizing the look of plotted motifs.

MINOR CHANGES

- read_matrix(..., positions): added partial argument matching.

- create_sequences(), shuffle_sequences(), motif_pvalue(): the c++
  random
  engine has been changed from std::default_random_engine to
  std::mt19937.
  This should allow for the same rng.seed value to result in the same
  output regardless of OS.

- score_match() has been vectorized (alongside new prob_match()
  function).

- The ape and ggtree packages are now no longer imported and must be
  installed seperately in order to use motif_tree().

- The processx package is no longer imported and must be installed
  seperately in order to use run_meme().

- The pseudocount slot is now shown when universalmotif class objects
  are
  printed.

- get_bkg(): the list.out and as.prob options have been disabled. To
  simplify the function, the only possible output (exception: if
  to.meme is
  not NULL) is a DataFrame showing both counts and probabilities.

- Changed the default look of motifs plotted by view_motifs().

- General documentation cleanup.

BUG FIXES

- Changing motif backgrounds with `[<-` will now make sure to set
  correct
  vector names.

- get_bkg() will now correctly ignore non-standard letters and letters
  missing from the provided alphabet during counting.

                        Changes in version 1.6.4                        

BUG FIXES

- cbind(): do not ignore the pseudocount slot.

- Fixed typo in IntroductionToSequenceMotifs.Rmd.

- Fixed U() function in IntroductionToSequenceMotifs.Rmd, no longer
  returns
  NA values if 0s are present.

- read_cisbp(): no parsing errors for motifs with missing/partial
  header
  info.

                        Changes in version 1.6.3                        

BUG FIXES

- scan_sequences(): commented out WIP code for scanning gapped motifs.

                        Changes in version 1.6.2                        

BUG FIXES

- motif_tree(): 'daylight' layout is no longer disabled.

                        Changes in version 1.6.1                        

BUG FIXES

- summarise_motif(): properly retrieves altname slot. Contribution from
  Spencer Nystrom (https://github.com/bjmt/universalmotif/pull/9).

- read_meme(): for LIKE type alphabets, make sure PROTEIN-LIKE is
  understood as being AA.

[variancePartition](https://bioconductor.org/packages/variancePartition)
-----------------

                       Changes in version 1.19.20                       

- fix bug discovered when the number of features is less than the
  number of chunks in iterBatch()

                       Changes in version 1.19.19                       

- simple bug fixes to pass R CMD check

                       Changes in version 1.19.18                       

- simplify calcVarPart for lm and lmer.  Add compatibility for glm

- Simplify checkModelStatus.merMod to allow formula (A|B) where A is
  continuous

- remove unused "adjust" arguments for clarity

                       Changes in version 1.19.17                       

- add get_prediction() for results of lm()

- improve documentation of get_prediction()

                       Changes in version 1.19.16                       

- in canCorPairs() change statistic used to summarize CCA to Cramer's
  V.  The difference is very subtle, but is now based on first
  principles.

- in dream, check that data is a data.frame

- dream() defaults to computeResiduals=TRUE for compatability with
  zenith

                       Changes in version 1.19.14                       

- fix issue with residuals() where examples fail

                       Changes in version 1.19.13                       

- fix issues with residuals()
  - https://github.com/GabrielHoffman/variancePartition/issues/18

- fix issue exporting eBayes, topTable, etc

                       Changes in version 1.19.12                       

- Improve documentation for contrasts in dream.Rmd

- check that contrasts sum to zero in plotContrasts.

                       Changes in version 1.19.11                       

- in voomWithDreamWeights() fix issue with not defining design
  - https://github.com/GabrielHoffman/variancePartition/issues/17

                       Changes in version 1.19.10                       

- in voomWithDreamWeights() fix issue with returning design matrix

- better error if counts can't be converted to matrix
  - https://github.com/GabrielHoffman/variancePartition/issues/15

                       Changes in version 1.19.7                        

- Round numbers in plotContrasts()

- fix issues with strings are passed to formula arguments

                       Changes in version 1.19.6                        

- New gives meaning full error message for dream(), etc when variable
  is not found in data.

                       Changes in version 1.19.5                        

- Better error catching when running fitVarPartModel() with fxn that
  fails

- add get_prediction() function
  - the following code now can be run in parallel
  fitList = fitVarPartModel( Y, ~ (1|Batch), data, fxn = function(fit){
  B = variancePartition::get_prediction(fit, ~(1|Batch))
  fit@resp$y - B
  }, BPPARAM=SnowParam(3))

                       Changes in version 1.19.4                        

- Update vignette #3, and update documentation of REML argument

                       Changes in version 1.19.3                        

- add new FAQ.Rmd

                       Changes in version 1.19.2                        

- canCorPairs() now returns NA correlation when two variables have
  no overlapping observed values

- plotCorrMatrix() now handles NA correlation values

                       Changes in version 1.19.1                        

- Bump to next Bioconductor version

                       Changes in version 1.18.3                        

- Improve documentation

- move location of eBayesFMT code

                       Changes in version 1.18.2                        

- Clean up some code and add documentation

- document ebayesFMT

                       Changes in version 1.18.1                        

- Clean up some code and add documentation

- compute effective degrees of freedom for each model

[VariantAnnotation](https://bioconductor.org/packages/VariantAnnotation)
-----------------

                       Changes in version 1.36.0                        

NEW FEATURES

- ref<-, alt<-, qual<- and filt<- allow replacement value length
  recycling

[velociraptor](https://bioconductor.org/packages/velociraptor)
------------

                       Changes in version 0.99.9                        

- Converted various functions to S4 generics for easier use with
SingleCellExperiment objects.

                       Changes in version 0.99.8                        

- Trigger new build to repeat ExperimentHub download.

                       Changes in version 0.99.7                        

- Delete empty line to force cache update. See
https://github.com/rubocop-hq/rubocop/pull/4342#issuecomment-305449759.

                       Changes in version 0.99.6                        

- Set autoscale=FALSE in the call to scvelo function
velocity_embedding to avoid issue related to Qt and plotting.

                       Changes in version 0.99.5                        

- Trigger new build to check if Windows issue resolved itself.

                       Changes in version 0.99.4                        

- Trigger new build to check whether TIMEOUT issue on Windows is
reproducible.

                       Changes in version 0.99.3                        

- Explicitly declare all Conda dependencies for scvelo.

                       Changes in version 0.99.2                        

- Add hexsticker.

                       Changes in version 0.99.1                        

- Remove .Rproj file from git repository.

                       Changes in version 0.99.0                        

- First submission to Bioconductor.

[VERSO](https://bioconductor.org/packages/VERSO)
-----

                        Changes in version 1.0.0                        

- Package released in October 2020.

[ViSEAGO](https://bioconductor.org/packages/ViSEAGO)
-------

                         Changes in version 1.3                         

- print graph with orca

- Ensembl2GO() biomart update

- show_heatmap() upgrade

- upset print update

- merge_enrich_terms upgrade pvalue cutoff

- merge_enrich_terms globale upgrade

- GOterms_heatmap remove row side colors text and correct showIC column

- vignettes update

- annotate() update for uniprot

- fgsea support

[VplotR](https://bioconductor.org/packages/VplotR)
------

                 Changes in version 0.5.0 (2020-05-25)                  

- IMPORTANT:

- Created plotProfile function

- MINOR:

- Changed color scale to scico package, roma scale

                 Changes in version 0.4.2 (2020-05-10)                  

- IMPORTANT:

- Improved nucleosomeEnrichment and background computation
- Treat strands separately in computeVmat
- Added plot theming

- MINOR:

- Clarified normalization approaches
- Improved documentation
- Make computeVmat a bit faster by removing some intermediate
steps

[wavClusteR](https://bioconductor.org/packages/wavClusteR)
----------

                       Changes in version 2.23.0                        

- Deprecated CWT-based cluster boundary identification

[weitrix](https://bioconductor.org/packages/weitrix)
-------

                       Changes in version 1.1.10                        

- Fix bug with weitrix_confects due to [[ ]] <- NULL deleting
elements
from a list instead of storing NULL in the list.

                        Changes in version 1.1.9                        

- Peaks were sometimes in reverse order in APA example, data file
updated.
- Try to get rid of an odd new build error about stack usage by using
serial processing in vignettes.

                        Changes in version 1.1.8                        

- Use geom_bin2d in weitrix_calplot scatterplots.

                        Changes in version 1.1.7                        

- Add mu_min, mu_max arguments to weitrix_calibrate_all.

                        Changes in version 1.1.6                        

- Use glm2, which is less prone to optimization failure.
- Auto-disable parallel processing if X11 device is open.
- Add SLAM-Seq vignette.

                        Changes in version 1.1.5                        

- well_knotted_spline for natural splines with good choice of knots.

                        Changes in version 1.1.4                        

- weitrix_confects for differential testing.
- weitrix_rms_confects now called weitrix_sd_confects.
- weitrix_calplot now shows mean trend and mean +/- standard
deviation
trend.

                        Changes in version 1.1.3                        

- weitrix_calibrate_all now includes a simple scaling factor to
account for residuals from a fitted model being smaller than
residuals from the true model.
- weitrix_calplot blue guidelines are similarly adjusted.
- Add weitrix_rms_confects to find rows with confidently excessive
variation.

                        Changes in version 1.1.2                        

- Vignettes use weitrix_calibrate_all, demonstrate weitrix_calplot.
- weitrix_calplot now uses sqrt(weight)*residual on y axis.

                        Changes in version 1.1.1                        

- Switch calibration from linear model on log dispersions to using a
gamma GLM.
- Add weitrix_calibrate_all function for very flexible calibration.
- Add weitrix_calplot to examine quality of weights.

[wpm](https://bioconductor.org/packages/wpm)
---

                       Changes in version 0.99.13                       

- Added a zzz.R file for the welcome message when loading the package
in an R session
- Corrections in package vignette and Help tab of the shiny app
- Revised all documentation of functions

                       Changes in version 0.99.12                       

- fixed the lack of display of images in the help tab
- fixed inconsistencies regarding the setting of quotes when
importing
data
- added an infobox for plate dimensions compatibility with the number
of samples

                       Changes in version 0.99.11                       

- fixed issue #10
- fixed issue #22
- changed the display of the datatable in the Results tab
- updated Vignette and Help page about toy dataset (CSV section)

                       Changes in version 0.99.10                       

- corrected unit tests for CSV import
- corrected README R commands for convertCSV section
- corrected Help tab for the upload file section
- corrected the Parameters section in the Vignette

                 Changes in version 0.99.9 (2020-06-18)                 

- Resolved crashes regarding bad settings when importing CSV
- Parameter sections are now collapsible
- created module for number of iterations
- defined the launch for the browser on port 8000
- changed the order of CSV import section in the parameters Panel
- some upload parameters are now available only if the CSV file is
correctly imported
- added preview of CSV file when importing in the UI

                 Changes in version 0.99.8 (2020-06-11)                 

- Modified structure of the Home tab in the shiny application
- Added new module for the help tab.
- Created the help.md file for the Help tab of the app.
- Added CSS for tables, and images.
- Changed the "blank" word to "buffer" throughout the package.
- Changed the "not random" term to "fixed" throughout the package.
- Updated images and text in the tutorial vignette.

                 Changes in version 0.99.7 (2020-06-05)                 

- Added the convertSE function to manage SummarizedExperiments
objects
in command line version
- Added the grouping factor option for data import (managed for both
shiny app and command line). Now the user can specify the column
name corresponding to the grouping factor wpm has to use for
backtracking.
- Modified app_ui structure regarding the project title input
- Updated the tutorial vignette's content
- Updated README file
- Revised unit tests for import functions
- The WPM package version output in Home panel is now obtained with
packageVersion()
- Created a new module for the integration of markdown files in the
shiny application.

                 Changes in version 0.99.6 (2020-06-02)                 

- Added function checkWpmInputs to control the correct use of the
wrapperWpm function (command line use).
- Updated the README file and the tutorial vignette explaining how to
use WPM using command line.

                 Changes in version 0.99.5 (2020-06-02)                 

- Added functions for importing CSV files and ExpressionSet / MSnSet
objects
- Added the wrapper function allowing to use wpm in command line
- Added support for eSet and MSnSet objects as input of the wrapper
function
- Added unit tests for import functions

                 Changes in version 0.99.4 (2020-05-28)                 

- Fix the TIMEOUT error during R CMD CHECK

                 Changes in version 0.99.3 (2020-05-28)                 

- All the R code is now in the /R directory.
- Modified R code structure by creating additional modules:
- mod_home, fusion between the Home panel and the Help Panel
- mod_data_export, a specific module for data export
- mod_plate_dimensions
- mod_special_wells, mainly to avoid code redundancy and gain
readability
- Added unit test for the convertVector2Df function
- Added Rd examples for convertVector2Df and drawMap functions

                       Changes in version 0.99.2                        

- Added a NEWS.md file to track changes to the package
- Modified code and files structure of the package
- Added Roxygen comments to functions

[xcms](https://bioconductor.org/packages/xcms)
----

                       Changes in version 3.11.8                        

- Disable parallel processing in vignettes.

                       Changes in version 3.11.7                        

- More efficient splitting data per file especially for larger data
  sets.

- Disable parallel processing in examples.

                       Changes in version 3.11.6                        

- Add `FilterIntensityParam` to filter chromatographic peaks on
  intensity
  (issue #502).

- Add `estimatePrecursorIntensity` function to determine the precursor
  intensity
  for MS2 spectra from the neighboring MS1 spectra.

                       Changes in version 3.11.4                        

- Change from `Spectra` and `Chromatograms` to `MSpectra` and
  `MChromatograms`
  from MSnbase version >= 2.15.3.

                       Changes in version 3.11.3                        

- `reconstructChromPeakSpectra`: report also polarity and
  `precusorIntensity`.

- `reconstructChromPeakSpectra`: ensure a retention time is reported
  for
  reconstructed MS2 spectra (issue #485).

- Change default for `expandRt` to `0` in
  `reconstructChromPeakSpectra`.

- Fix error in `refineChromPeaks,MergeNeighboringPeaksParam` if no
  peaks found
  to be merged.

                       Changes in version 3.11.2                        

- Add `fillChromPeaks,ChromPeakAreaParam` to base the area from which
  missing
  peak data should be filled-in on the actually detected
  chromatographic peaks
  of a feature.

- Potential fix for issue #481: function should no longer throw an
  error because
  retention times are of length 0.

- More efficient splitting of processing which should increase the
  speed of
  the findChromPeaks, refineChromPeaks, reconstructChromPeakSpectra and
  chromPeakSpectra calls.

                       Changes in version 3.11.1                        

- Fix issue #471: conversion from `XCMSnExp` to `xcmsSet` looses
  phenodata
  (thanks to Andris Jankevics for reporting and providing a solution).

- Add `normalize` method for `Chromatogram` and `Chromatograms`
  objects.

- `featureChromatograms` gets new parameter `n` and `value` to extract
  EICs
  only from the top n samples with highest intensities.

- `filterFile` gets new parameter `keepFeatures` to support retaining
  correspondence results even if a data set is filtered by file.

- Export the virtual `Param` class.

- Add filterColumnsIntensityAbove method for Chromatograms object that
  allows
  to select columns (samples) of an Chromatograms object for which
  intensities
  of its chromatographic data are higher than a threshold.

- Add removeIntensity method for Chromatogram, Chromatograms,
  XChromatogram
  and XChromatograms objects allowing to *remove* intensities based on
  different
  criteria.

- Add correlate method for Chromatograms allowing to correlate multiple
  chromatograms with each other.

[zellkonverter](https://bioconductor.org/packages/zellkonverter)
-------------

                        Changes in version 1.0.0                        

- Accepted into Bioconductor for Release 3.12

- 
  zellkonverter provides methods to convert between Python
  AnnData objects and SingleCellExperiment objects. These are
  primarily intended for use by downstream Bioconductor packages
  that wrap Python methods for single-cell data analysis. It also
  includes functions to read and write H5AD files used for saving
  AnnData objects to disk.

[zinbwave](https://bioconductor.org/packages/zinbwave)
--------

                 Changes in version 1.11.6 (2020-07-18)                 

- Fixed a bug in the initialization of beta_j.

- Fixed a bug in zinbsurf.

- Changed zinbwave default to `K=2`.

- Fix bug in the initialization of W.


NEWS from new and existing Data Experiment Packages
===================================


[chipseqDBData](https://bioconductor.org/packages/chipseqDBData)
-------------

                        Changes in version 1.6.0                        

- Return BamFile objects from all getter functions, to properly
  capture index files with different names from the BAM files.

[clustifyrdatahub](https://bioconductor.org/packages/clustifyrdatahub)
----------------

                 Changes in version 0.99.4 (2020-09-02)                 

- Add mouse atlas



                 Changes in version 0.99.0 (2020-07-02)                 

- Cleanup for bioc

[curatedTCGAData](https://bioconductor.org/packages/curatedTCGAData)
---------------

                       Changes in version 1.12.0                        

Bug fixes and minor improvements

- Output dataset options as table when dry.run is enabled in the main
  function.

- Check for RaggedExperiment dependency when loading data that uses
  the data representation (@vjcitn, #39)

[depmap](https://bioconductor.org/packages/depmap)
------

                         Changes in version 1.3                         

Changes in version 1.3.2

- 20Q3 data added for crispr, copyNumber, TPM, mutationCalls and
metadata datasets. Newer versions for the other datasets were not
released.

Changes in version 1.3.1

- 20Q2 data added
- Included new proteomic dataset
- expression variable in TPM dataset changed to rna_expression



[dorothea](https://bioconductor.org/packages/dorothea)
--------

                 Changes in version 1.1.2 (2020-10-08)                  

- Changed TF census from TFclass to the more recent version from
  Lambert et al.. Information of mode of regulation for each TF
  (activator, supressor, dual) is still taken from Garcia-Alonso et
  al..

- Updated deprecated gene symbols to their latest version with the
  limma package (version 3.44.3).

- Shifted viper package from suggest to depends in the DESCRIPTION
  file.

- Added a further argument specifially for run_viper().Seurat to
  select a specific assay name to extract the normalized gene
  expression values from.

                 Changes in version 1.1.1 (2020-09-02)                  

- Export df2regulon function

- Improved documentation (added gh page URL to DESCRIPTION)

                 Changes in version 1.0.1 (2020-08-13)                  

- Improved package documentation

- Updated link to 10x genomics data set in single-cell vignette

- Fixed tests related to Seurat and SCE class

[dsQTL](https://bioconductor.org/packages/dsQTL)
-----

                        Changes in version 2.17                         

USER VISIBLE CHANGES

- ch2locs (retrievable via dsQTL::getSNPlocs) has been changed at about
  1850 locations
  where rs numbers had been associated with hg19 addresses; the dsQTL
  regions are hg18
  as are all the chr2... SNP addresses.  Previously the discoverable rs
  numbers used
  in the Chicago distribution from
  http://eqtl.uchicago.edu/dsQTL_data/GENOTYPES/ had
  be mapped via SNPlocs...20111119, but now they come directly from the
  Chicago text file.

[FieldEffectCrc](https://bioconductor.org/packages/FieldEffectCrc)
--------------

                 Changes in version 0.99.2 (2020-09-20)                 

- (0.99.2) Added PMCID and PMID in reference to original study

- (0.99.2) Updated citation to include full reference

- (0.99.2) Minor improvements in readability of R code

                 Changes in version 0.99.1 (2020-07-24)                 

- (0.99.1) Simplified vignette Installation code

- (0.99.1) Updated vignette typical-filter code chunk to evaluate as R
  code

                 Changes in version 0.99.0 (2020-07-14)                 

- (0.99.0) Passed R CMD check and R CMD BiocCheck without errors or
  warnings

- (0.99.0) Submitted to Bioconductor

                 Changes in version 0.1.0 (2020-07-03)                  

- (0.1.0) Begin building FieldEffectCrc package

[MethylSeqData](https://bioconductor.org/packages/MethylSeqData)
-------------

                        Changes in version 1.0.0                        

- New package MethylSeqData, providing DNA methylation sequencing
  datasets.

[pRolocdata](https://bioconductor.org/packages/pRolocdata)
----------

                       Changes in version 1.27.1                        

- add data from Kozik et al. (2020)

[RegParallel](https://bioconductor.org/packages/RegParallel)
-----------

                        Changes in version 1.8.0                        

- added funtionality for include survey weights via survey::svyglm

- fixed bug whereby, in ceratin situations, an incorrect number of
  blocks would be calculated

[RforProteomics](https://bioconductor.org/packages/RforProteomics)
--------------

                       Changes in version 1.27.1                        

- Remove deprecated IPPD dependency

[SCATEData](https://bioconductor.org/packages/SCATEData)
---------

                        Changes in version 1.0.0                        

- submission to Bioconductor

[scRNAseq](https://bioconductor.org/packages/scRNAseq)
--------

                        Changes in version 2.4.0                        

- Added the Zilionis lung dataset (Jens Preussner).

- Added the Hermann spermatogenesis dataset (Charlotte Soneson).

- Added the Mair and Kotliarov PBMC datasets (Stephany Orjuela).

- Added the Stoeckius cell hashing dataset.

- Added the Wu kidney snRNA-seq dataset.

- Added the Hu cortex snRNA-seq dataset.

- Added spike-in concentrations to the altExp rowData for various
  datasets (Alan O'Callaghan).

[signatureSearchData](https://bioconductor.org/packages/signatureSearchData)
-------------------

                 Changes in version 1.2.1 (2020-08-14)                  

- Edited vignette containing 'DEGs and Cutoffs Definition' subsections
  to document how to defining query DEGs and gene set reference
  database by setting LFC scores and FDR cutoffs.

[SingleCellMultiModal](https://bioconductor.org/packages/SingleCellMultiModal)
--------------------

                        Changes in version 1.2.0                        

New features

- CITEseq function, vignette, and 'cord_blood' data available
  (@drighelli, #18)

- Include seqFISH function, vignette, and 'mouse_visual_cortex' data
  (v1 and v2 from @drighelli, #14)

- New 'mouse_gastrulation' dataset released (version "2.0.0").

- Use version argument to indicate the mouse_gastrulation data version

- The data includes all cells not only the ones that passed the QC of
  all three 'omics (thanks @rargelaguet, @ajabadi).

Bug fixes and minor improvements

- Caching mechanism uses tools::R_user_dir and not rappdirs.

- Improved display of available data using ExperimentHub metadata.

- Improved documentation explaining versioning differences.

- Contribution guidelines available at
  https://github.com/waldronlab/SingleCellMultiModal/wiki/Contributing-Guidelines

- Default version argument in scNMT function now set to "2.0.0"
  (version "1.0.0" still available)

[spatialLIBD](https://bioconductor.org/packages/spatialLIBD)
-----------

                        Changes in version 1.1.5                        

NEW FEATURES

- fetch_data() takes the data from sce object and creates a
  VisiumExperiment object containing these data. VisiumExperiment
  object can be obtained with fetch_data("ve").

                        Changes in version 1.1.4                        

NEW FEATURES

- fetch_data() now uses BiocFileCache() when downloading the data from
  Dropbox.

[timecoursedata](https://bioconductor.org/packages/timecoursedata)
--------------

                 Changes in version 0.99.0 (2020-05-15)                 

- Submitted to Bioconductor

[TMExplorer](https://bioconductor.org/packages/TMExplorer)
----------

                       Changes in version 0.99.6                        

- Renamed Vignette to TMExplorer ("tutorial" was too generic)

- More man page updates

                       Changes in version 0.99.5                        

- Updated man pages and cleaned up some code.

                 Changes in version 0.99.3 (2020-07-03)                 

- Fixed a bug when downloading data on Windows

                 Changes in version 0.99.2 (2020-07-03)                 

- Fixed missing cell type labels in some datasets

                 Changes in version 0.99.1 (2020-06-25)                 

- Now uses SingleCellExperiment objects

                 Changes in version 0.99.0 (2020-06-22)                 

- Submitted to Bioconductor


NEWS from new and existing Workflows
===================================


[fluentGenomics](https://bioconductor.org/packages/fluentGenomics)
--------------

                        Changes in version 1.0.1                        

- fix downloading of ATAC seq files for windows


Deprecated and Defunct Packages
===============================

Fifty Six software packages were removed from this release (after being deprecated
in Bioc 3.11): affypdnn, AnalysisPageServer, anamiR, BayesPeak, bgafun, biosvd, birta, CALIB,
CAMTHC, cellGrowth, chroGPS, cobindR, CTDquerier, CVE, DChIPRep, DEDS,
DupChecker, FEM, gCMAP, gCMAPWeb, geecc, Genominator, IdMappingAnalysis,
IdMappingRetrieval, IPPD, kimod, LMGene, lol, LVSmiRNA, M3D, manta,
MaxContrastProjection, MCRestimate, MergeMaid, mitoODE, MoPS, motifRG, MTseeker,
nem, PAPi, pcaGoPromoter, pint, plw, PowerExplorer, proteoQC, QUALIFIER, readat,
RefNet, RIPSeeker, SANTA, scfind, splicegear, sRAP, triform, Vega, waveTiling

Sixty eight software are deprecated in this release and will be removed in Bioc 3.13:
adaptest, ArrayTV, BioSeqClass, CGEN, CHARGE, chimera, CNVtools, CorMut, DESeq,
explorase, flowFit, flowSpy, flowType, flowVS, focalCall, FourCSeq, FunciSNP,
GeneticsDesign, GenRank, GGBase, GGtools, GOFunction, gQTLBase, gQTLstats,
hicrep, ImpulseDE, ImpulseDE2, joda, JunctionSeq, LINC, Logolas, mcaGUI,
metaArray, metaseqR, methVisual, methyvim, Mirsynergy, MmPalateMiRNA, MOFA,
MotIV, NarrowPeaks, netbenchmark, netReg, OGSA, OmicsMarkeR, pathprint,
PathwaySplice, PGA, PGSEA, plrs, prada, Prize, Rariant, reb, Roleswitch,
rTANDEM, sampleClassifier, sapFinder, scsR, shinyTANDEM, sigaR, signet,
simpleaffy, spotSegmentation, Starr, SVAPLSseq, TxRegInfra, xps

Two experimental data packages were removed in this release (after being
deprecated in BioC 3.11): MTseekerData, RIPSeekerData

Sixteen experimental data packages are deprecated in this release and will be
removed in Bioc 3.13: flowFitExampleData, FunciSNP.data, geuvPack, geuvStore2,
GGdata, methyvimData, mitoODEdata, Mulder2012, pathprintGEOData,
pcaGoPromoter.Hs.hg19, pcaGoPromoter.Mm.mm9, pcaGoPromoter.Rn.rn4,
RNAinteractMAPK, sampleClassifierData, waveTilingData, yriMulti

Nine annotation packages were removed this release: 
hom.At.inp.db, hom.Ce.inp.db, hom.Dm.inp.db, hom.Dr.inp.db,
hom.Hs.inp.db, hom.Mm.inp.db, hom.Rn.inp.db, hom.Sc.inp.db, KEGG.db. 

No annotation packages are deprecated in this release and will be removed in
Bioc 3.13.

No workflow package were removed in this release (after being deprecated in BioC
3.11).

No workflow packages were deprecated in this release.