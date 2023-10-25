October 25, 2023

Bioconductors:

We are pleased to announce Bioconductor 3.18, consisting of
2266 software packages, 429 experiment data packages, 918 annotation
packages, 30 workflows and 4 books.

There are 69 new software packages, 10 new data experiment packages,
8 new annotation packages, no new workflows, 1 new books, and many updates and
improvements to existing packages.

Bioconductor 3.18 is compatible with R 4.3, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher and
macOS arm64. This release will also include updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.18](#getting-started-with-bioconductor-318)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [NEWS from existing workflows](#news-from-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.18
======================================

To update to or install Bioconductor 3.18:

1. Install R 4.3. Bioconductor 3.18 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 69 new software packages in this release of Bioconductor.

- [adverSCarial](/packages/adverSCarial) adverSCarial is an R Package
  designed for generating and analyzing the vulnerability of
  scRNA-seq classifiers to adversarial attacks. The package is
  versatile and provides a format for integrating any type of
  classifier. It offers functions for studying and generating two
  types of attacks, single gene attack and max change attack. The
  single gene attack involves making a small modification to the
  input to alter the classification. The max change attack involves
  making a large modification to the input without changing its
  classification. The package provides a comprehensive solution for
  evaluating the robustness of scRNA-seq classifiers against
  adversarial attacks.

- [alabaster.files](/packages/alabaster.files) Save common
  bioinformatics file formats within the alabaster framework. This
  includes BAM, BED, VCF, bigWig, bigBed, FASTQ, FASTA and so on. We
  save and load additional metadata for each file, and we support
  linkage between each file and its corresponding index.

- [beachmat.hdf5](/packages/beachmat.hdf5) Extends beachmat to
  support initialization of tatami matrices from HDF5-backed arrays.
  This allows C++ code in downstream packages to directly call the
  HDF5 C/C++ library to access array data, without the need for block
  processing via DelayedArray. Some utilities are also provided for
  direct creation of an in-memory tatami matrix from a HDF5 file.

- [BioCartaImage](/packages/BioCartaImage) The core functionality of
  the package is to provide coordinates of genes on the BioCarta
  pathway images and to provide methods to add self-defined graphics
  to the genes of interest.

- [BiocBook](/packages/BiocBook) A BiocBook can be created by authors
  (e.g. R developers, but also scientists, teachers, communicators,
  ...) who wish to 1) write (compile a body of biological and/or
  bioinformatics knowledge), 2) containerize (provide Docker images
  to reproduce the examples illustrated in the compendium), 3)
  publish (deploy an online book to disseminate the compendium), and
  4) version (automatically generate specific online book versions
  and Docker images for specific Bioconductor releases).

- [CaDrA](/packages/CaDrA) Performs both stepwise and backward
  heuristic search for candidate (epi)genetic drivers based on a
  binary multi-omics dataset. CaDrA's main objective is to identify
  features which, together, are significantly skewed or enriched
  pertaining to a given vector of continuous scores (e.g.
  sample-specific scores representing a phenotypic readout of
  interest, such as protein expression, pathway activity, etc.),
  based on the union occurence (i.e. logical OR) of the events.

- [CardinalIO](/packages/CardinalIO) Fast and efficient reading and
  writing of mass spectrometry imaging data files. Supports imzML and
  Analyze 7.5 formats. Provides ontologies for mass spectrometry
  imaging.

- [CCPlotR](/packages/CCPlotR) CCPlotR is an R package for
  visualising results from tools that predict cell-cell interactions
  from single-cell RNA-seq data. These plots are generic and can be
  used to visualise results from multiple tools such as Liana,
  CellPhoneDB, NATMI etc.

- [CDI](/packages/CDI) Single-cell RNA-sequencing (scRNA-seq) is
  widely used to explore cellular variation. The analysis of
  scRNA-seq data often starts from clustering cells into
  subpopulations. This initial step has a high impact on downstream
  analyses, and hence it is important to be accurate. However, there
  have not been unsupervised metric designed for scRNA-seq to
  evaluate clustering performance. Hence, we propose clustering
  deviation index (CDI), an unsupervised metric based on the modeling
  of scRNA-seq UMI counts to evaluate clustering of cells.

- [cfdnakit](/packages/cfdnakit) This package provides basic
  functions for analyzing shallow whole-genome sequencing (~0.3X or
  more) of cell-free DNA (cfDNA). The package basically extracts the
  length of cfDNA fragments and aids the vistualization of
  fragment-length information. The package also extract
  fragment-length information per non-overlapping fixed-sized bins
  and used it for calculating ctDNA estimation score (CES).

- [ClustIRR](/packages/ClustIRR) ClustIRR is a quantitative method
  for clustering of immune receptor repertoires (IRRs). The algorithm
  identifies groups of T or B cell receptors (TCRs or BCRs) with
  similar specificity by comparing their sequences. ClustIRR uses
  graphs to visualize the specificity structures of IRRs.

- [compSPOT](/packages/compSPOT) Clonal cell groups share common
  mutations within cancer, precancer, and even clinically normal
  appearing tissues. The frequency and location of these mutations
  may predict prognosis and cancer risk. It has also been well
  established that certain genomic regions have increased sensitivity
  to acquiring mutations. Mutation-sensitive genomic regions may
  therefore serve as markers for predicting cancer risk. This package
  contains multiple functions to establish significantly mutated
  hotspots, compare hotspot mutation burden between samples, and
  perform exploratory data analysis of the correlation between
  hotspot mutation burden and personal risk factors for cancer, such
  as age, gender, and history of carcinogen exposure. This package
  allows users to identify robust genomic markers to help establish
  cancer risk.

- [CuratedAtlasQueryR](/packages/CuratedAtlasQueryR) Provides access
  to a copy of the Human Cell Atlas, but with harmonised metadata.
  This allows for uniform querying across numerous datasets within
  the Atlas using common fields such as cell type, tissue type, and
  patient ethnicity. Usage involves first querying the metadata table
  for cells of interest, and then downloading the corresponding cells
  into a SingleCellExperiment object.

- [CytoPipelineGUI](/packages/CytoPipelineGUI) This package is the
  companion of the CytoPipeline package. It provides GUI's (shiny
  apps) for the visualization of flow cytometry data analysis
  pipelines that are run with CytoPipeline. Two shiny applications
  are provided, i.e. an interactive flow frame assessment and
  comparison tool and an interactive scale transformations
  visualization and adjustment tool.

- [DCATS](/packages/DCATS) Methods to detect the differential
  composition abundances between conditions in singel-cell RNA-seq
  experiments, with or without replicates. It aims to correct bias
  introduced by missclaisification and enable controlling of
  confounding covariates. To avoid the influence of proportion change
  from big cell types, DCATS can use either total cell number or
  specific reference group as normalization term.

- [decontX](/packages/decontX) This package contains implementation
  of DecontX (Yang et al. 2020), a decontamination algorithm for
  single-cell RNA-seq, and DecontPro (Yin et al. 2023), a
  decontamination algorithm for single cell protein expression data.
  DecontX is a novel Bayesian method to computationally estimate and
  remove RNA contamination in individual cells without empty droplet
  information. DecontPro is a Bayesian method that estimates the
  level of contamination from ambient and background sources in
  CITE-seq ADT dataset and decontaminate the dataset.

- [demuxSNP](/packages/demuxSNP) This package assists in
  demultiplexing scRNAseq data using both cell hashing and SNPs data.
  The SNP profile of each group os learned using high confidence
  assignments from the cell hashing data. Cells which cannot be
  assigned with high confidence from the cell hashing data are
  assigned to their most similar group based on their SNPs. We also
  provide some helper function to optimise SNP selection, create
  training data and merge SNP data into the SingleCellExperiment
  framework.

- [dreamlet](/packages/dreamlet) Recent advances in single
  cell/nucleus transcriptomic technology has enabled collection of
  cohort-scale datasets to study cell type specific gene expression
  differences associated disease state, stimulus, and genetic
  regulation. The scale of these data, complex study designs, and low
  read count per cell mean that characterizing cell type specific
  molecular mechanisms requires a user-frieldly, purpose-build
  analytical framework. We have developed the dreamlet package that
  applies a pseudobulk approach and fits a regression model for each
  gene and cell cluster to test differential expression across
  individuals associated with a trait of interest. Use of
  precision-weighted linear mixed models enables accounting for
  repeated measures study designs, high dimensional batch effects,
  and varying sequencing depth or observed cells per biosample.

- [easylift](/packages/easylift) The easylift package provides a
  convenient tool for genomic liftover operations between different
  genome assemblies. It seamlessly works with Bioconductor's GRanges
  objects and chain files from the UCSC Genome Browser, allowing for
  straightforward handling of genomic ranges across various genome
  versions. One noteworthy feature of easylift is its integration
  with the BiocFileCache package. This integration automates the
  management and caching of chain files necessary for liftover
  operations. Users no longer need to manually specify chain file
  paths in their function calls, reducing the complexity of the
  liftover process.

- [enrichViewNet](/packages/enrichViewNet) This package enables the
  visualization of functional enrichment results as network graphs.
  First the package enables the visualization of enrichment results,
  in a format corresponding to the one generated by gprofiler2, as a
  customizable Cytoscape network. In those networks, both gene
  datasets (GO terms/pathways/protein complexes) and genes associated
  to the datasets are represented as nodes. While the edges connect
  each gene to its dataset(s). The package also provides the option
  to create enrichment maps from functional enrichment results.
  Enrichment maps enable the visualization of enriched terms into a
  network with edges connecting overlapping genes.

- [fenr](/packages/fenr) Perform fast functional enrichment on
  feature lists (like genes or proteins) using the hypergeometric
  distribution. Tailored for speed, this package is ideal for
  interactive platforms such as Shiny. It supports the retrieval of
  functional data from sources like GO, KEGG, Reactome, and
  WikiPathways. By downloading and preparing data first, it allows
  for rapid successive tests on various feature selections without
  the need for repetitive, time-consuming preparatory steps typical
  of other packages.

- [gatom](/packages/gatom) This package implements a metabolic
  network analysis pipeline to identify an active metabolic module
  based on high throughput data. The pipeline takes as input
  transcriptional and/or metabolic data and finds a metabolic
  subnetwork (module) most regulated between the two conditions of
  interest. The package further provides functions for module
  post-processing, annotation and visualization.

- [gDNAx](/packages/gDNAx) Provides diagnostics for assessing genomic
  DNA contamination in RNA-seq data, as well as plots representing
  these diagnostics. Moreover, the package can be used to get an
  insight into the strand library protocol used and, in case of
  strand-specific libraries, the strandedness of the data.
  Furthermore, it provides functionality to filter out reads of
  potential gDNA origin.

- [gDR](/packages/gDR) Package is a part of the gDR suite. It
  reexports functions from other packages in the gDR suite that
  contain critical processing functions and utilities. The vignette
  walks through the full processing pipeline for drug response
  analyses that the gDR suite offers.

- [gDRcore](/packages/gDRcore) This package contains core functions
  to process and analyze drug response data. The package provides
  tools for normalizing, averaging, and calculation of gDR metrics
  data. All core functions are wrapped into the pipeline function
  allowing analyzing the data in a straightforward way.

- [gDRimport](/packages/gDRimport) The package is a part of the gDR
  suite. It helps to prepare raw drug response data for downstream
  processing. It mainly contains helper functions for
  importing/loading/validating dose-response data provided in
  different file formats.

- [gDRstyle](/packages/gDRstyle) Package fills a helper package role
  for whole gDR suite. It helps to support good development practices
  by keeping style requirements and style tests for other packages.
  It also contains build helpers to make all package requirements
  met.

- [gDRutils](/packages/gDRutils) This package contains utility
  functions used throughout the gDR platform to fit data, manipulate
  data, and convert and validate data structures. This package also
  has the necessary default constants for gDR platform. Many of the
  functions are utilized by the gDRcore package.

- [GenomicPlot](/packages/GenomicPlot) Visualization of next
  generation sequencing (NGS) data is essential for interpreting
  high-throughput genomics experiment results. 'GenomicPlot'
  facilitates plotting of NGS data in various formats (bam, bed, wig
  and bigwig); both coverage and enrichment over input can be
  computed and displayed with respect to genomic features (such as
  UTR, CDS, enhancer), and user defined genomic loci or regions.
  Statistical tests on signal intensity within user defined regions
  of interest can be performed and represented as boxplots or bar
  graphs. Parallel processing is used to speed up computation on
  multicore platforms. In addition to genomic plots which is suitable
  for displaying of coverage of genomic DNA (such as ChIPseq data),
  metagenomic (without introns) plots can also be made for RNAseq or
  CLIPseq data as well.

- [gg4way](/packages/gg4way) 4way plots enable a comparison of the
  logFC values from two contrasts of differential gene expression.
  The gg4way package creates 4way plots using the ggplot2 framework
  and supports popular Bioconductor objects. The package also
  provides information about the correlation between contrasts and
  significant genes of interest.

- [ggkegg](/packages/ggkegg) This package aims to import, parse, and
  analyze KEGG data such as KEGG PATHWAY and KEGG MODULE. The package
  supports visualizing KEGG information using ggplot2 and ggraph
  through using the grammar of graphics. The package enables the
  direct visualization of the results from various omics analysis
  packages.

- [ggsc](/packages/ggsc) Useful functions to visualize single cell
  and spatial data. It supports both 'SingleCellExperiment' and
  'Seurat' objects. It also supports visualizing the data using
  grammar of graphics implemented in 'ggplot2'.

- [GloScope](/packages/GloScope) This package aims at representing
  and summarizing the entire single-cell profile of a sample. It
  allows researchers to perform important bioinformatic analyses at
  the sample-level such as visualization and quality control. The
  main functions Estimate sample distribution and calculate
  statistical divergence among samples, and visualize the distance
  matrix through MDS plots.

- [GNOSIS](/packages/GNOSIS) GNOSIS incorporates a range of R
  packages enabling users to efficiently explore and visualise
  clinical and genomic data obtained from cBioPortal. GNOSIS uses an
  intuitive GUI and multiple tab panels supporting a range of
  functionalities. These include data upload and initial exploration,
  data recoding and subsetting, multiple visualisations, survival
  analysis, statistical analysis and mutation analysis, in addition
  to facilitating reproducible research.

- [HarmonizR](/packages/HarmonizR) An implementation, which takes
  input data and makes it available for proper batch effect removal
  by ComBat or Limma. The implementation appropriately handles
  missing values by dissecting the input matrix into smaller matrices
  with sufficient data to feed the ComBat or limma algorithm. The
  adjusted data is returned to the user as a rebuild matrix. The
  implementation is meant to make as much data available as possible
  with minimal data loss.

- [HERON](/packages/HERON) HERON is a software package for analyzing
  peptide binding array data. In addition to identifying significant
  binding probes, HERON also provides functions for finding epitopes
  (string of consecutive peptides within a protein). HERON also
  calculates significance on the probe, epitope, and protein level by
  employing meta p-value methods.  HERON is designed for obtaining
  calls on the sample level and calculates fractions of hits for
  different conditions.

- [hicVennDiagram](/packages/hicVennDiagram) A package to generate
  high-resolution Venn and Upset plots for genomic interaction data
  from HiC, ChIA-PET, HiChIP, PLAC-Seq, Hi-TrAC, HiCAR and etc. The
  package generates plots specifically crafted to eliminate the
  deceptive visual representation caused by the counts method.

- [hoodscanR](/packages/hoodscanR) hoodscanR is an user-friendly R
  package providing functions to assist cellular neighborhood
  analysis of any spatial transcriptomics data with single-cell
  resolution. All functions in the package are built based on the
  SpatialExperiment object, allowing integration into various spatial
  transcriptomics-related packages from Bioconductor. The package can
  result in cell-level neighborhood annotation output, along with
  funtions to perform neighborhood colocalization analysis and
  neighborhood-based cell clustering.

- [iNETgrate](/packages/iNETgrate) The iNETgrate package provides
  functions to build a correlation network in which nodes are genes.
  DNA methylation and gene expression data are integrated to define
  the connections between genes. This network is used to identify
  modules (clusters) of genes. The biological information in each of
  the resulting modules is represented by an eigengene. These
  biological signatures can be used as features e.g., for
  classification of patients into risk categories. The resulting
  biological signatures are very robust and give a holistic view of
  the underlying molecular changes.

- [iSEEde](/packages/iSEEde) This package contains diverse
  functionality to extend the usage of the iSEE package, including
  additional classes for the panels or modes facilitating the
  analysis of differential expression results. This package does not
  perform differential expression. Instead, it provides methods to
  embed precomputed differential expression results in a
  SummarizedExperiment object, in a manner that is compatible with
  interactive visualisation in iSEE applications.

- [iSEEindex](/packages/iSEEindex) This package provides an interface
  to any collection of data sets within a single iSEE
  web-application. The main functionality of this package is to
  define a custom landing page allowing app maintainers to list a
  custom collection of data sets that users can selected from and
  directly load objects into an iSEE web-application.

- [iSEEpathways](/packages/iSEEpathways) This package contains
  diverse functionality to extend the usage of the iSEE package,
  including additional classes for the panels or modes facilitating
  the analysis of pathway analysis results. This package does not
  perform pathway analysis. Instead, it provides methods to embed
  precomputed pathway analysis results in a SummarizedExperiment
  object, in a manner that is compatible with interactive
  visualisation in iSEE applications.

- [IsoBayes](/packages/IsoBayes) IsoBayes is a Bayesian method to
  perform inference on single protein isoforms. Our approach infers
  the presence/absence of protein isoforms, and also estimates their
  abundance; additionally, it provides a measure of the uncertainty
  of these estimates, via: i) the posterior probability that a
  protein isoform is present in the sample; ii) a posterior credible
  interval of its abundance. IsoBayes inputs liquid cromatography
  mass spectrometry (MS) data, and can work with both PSM counts, and
  intensities. When available, trascript isoform abundances (i.e.,
  TPMs) are also incorporated: TPMs are used to formulate an
  informative prior for the respective protein isoform relative
  abundance. We further identify isoforms where the relative
  abundance of proteins and transcripts significantly differ. We use
  a two-layer latent variable approach to model two sources of
  uncertainty typical of MS data: i) peptides may be erroneously
  detected (even when absent); ii) many peptides are compatible with
  multiple protein isoforms. In the first layer, we sample the
  presence/absence of each peptide based on its estimated probability
  of being mistakenly detected, also known as PEP (i.e., posterior
  error probability). In the second layer, for peptides that were
  estimated as being present, we allocate their abundance across the
  protein isoforms they map to. These two steps allow us to recover
  the presence and abundance of each protein isoform.

- [lemur](/packages/lemur) Fit a latent embedding multivariate
  regression (LEMUR) model to multi-condition single-cell data. The
  model provides a parametric description of single-cell data
  measured with complex experimental designs. The parametric model is
  used to (1) align conditions, (2) predict log fold changes between
  conditions for all cells, and (3) identify cell neighborhoods with
  consistent log fold changes. For those neighborhoods, a
  pseudobulked differential expression test is conducted to assess
  which genes are significantly changed.

- [MICSQTL](/packages/MICSQTL) Our pipeline, MICSQTL, utilizes
  scRNA-seq reference and bulk transcriptomes to estimate cellular
  composition in the matched bulk proteomes. The expression of genes
  and proteins at either bulk level or cell type level can be
  integrated by Angle-based Joint and Individual Variation Explained
  (AJIVE) framework. Meanwhile, MICSQTL can perform cell-type-specic
  quantitative trait loci (QTL) mapping to proteins or transcripts
  based on the input of bulk expression data and the estimated
  cellular composition per molecule type, without the need for single
  cell sequencing. We use matched transcriptome-proteome from human
  brain frontal cortex tissue samples to demonstrate the input and
  output of our tool.

- [Moonlight2R](/packages/Moonlight2R) The understanding of cancer
  mechanism requires the identification of genes playing a role in
  the development of the pathology and the characterization of their
  role (notably oncogenes and tumor suppressors). We present an
  updated version of the R/bioconductor package called MoonlightR,
  namely Moonlight2R, which returns a list of candidate driver genes
  for specific cancer types on the basis of omics data integration.
  The Moonlight framework contains a primary layer where gene
  expression data and information about biological processes are
  integrated to predict genes called oncogenic mediators, divided
  into putative tumor suppressors and putative oncogenes. This is
  done through functional enrichment analyses, gene regulatory
  networks and upstream regulator analyses to score the importance of
  well-known biological processes with respect to the studied cancer
  type. By evaluating the effect of the oncogenic mediators on
  biological processes or through random forests, the primary layer
  predicts two putative roles for the oncogenic mediators: i) tumor
  suppressor genes (TSGs) and ii) oncogenes (OCGs). As gene
  expression data alone is not enough to explain the deregulation of
  the genes, a second layer of evidence is needed. We have automated
  the integration of a secondary mutational layer through new
  functionalities in Moonlight2R. These functionalities analyze
  mutations in the cancer cohort and classifies these into driver and
  passenger mutations using the driver mutation prediction tool,
  CScape-somatic. Those oncogenic mediators with at least one driver
  mutation are retained as the driver genes. As a consequence, this
  methodology does not only identify genes playing a dual role (e.g.
  TSG in one cancer type and OCG in another) but also helps in
  elucidating the biological processes underlying their specific
  roles. In particular, Moonlight2R can be used to discover OCGs and
  TSGs in the same cancer type. This may for instance help in
  answering the question whether some genes change role between early
  stages (I, II) and late stages (III, IV). In the future, this
  analysis could be useful to determine the causes of different
  resistances to chemotherapeutic treatments.

- [MSstatsBig](/packages/MSstatsBig) MSstats package provide tools
  for preprocessing, summarization and differential analysis of mass
  spectrometry (MS) proteomics data. Recently, some MS protocols
  enable acquisition of data sets that result in larger than memory
  quantitative data. MSstats functions are not able to process such
  data. MSstatsBig package provides additional converter functions
  that enable processing larger than memory data sets.

- [MultiRNAflow](/packages/MultiRNAflow) Our R package MultiRNAflow
  provides an easy to use unified framework allowing to automatically
  make both unsupervised and supervised (DE) analysis for datasets
  with an arbitrary number of biological conditions and time points.
  In particular, our code makes a deep downstream analysis of DE
  information, e.g. identifying temporal patterns across biological
  conditions and DE genes which are specific to a biological
  condition for each time.

- [multiWGCNA](/packages/multiWGCNA) An R package for deeping mining
  gene co-expression networks in multi-trait expression data.
  Provides functions for analyzing, comparing, and visualizing WGCNA
  networks across conditions. multiWGCNA was designed to handle the
  common case where there are multiple biologically meaningful sample
  traits, such as disease vs wildtype across development or
  anatomical region.

- [nipalsMCIA](/packages/nipalsMCIA) Computes Multiple Co-Inertia
  Analysis (MCIA), a dimensionality reduction (jDR) algorithm, for a
  multi-block dataset using a modification to the Nonlinear Iterative
  Partial Least Squares method (NIPALS) proposed in (Hanafi et. al,
  2010). Allows multiple options for row- and table-level
  preprocessing, and speeds up computation of variance explained.
  Vignettes detail application to bulk- and single cell- multi-omics
  studies.

- [orthos](/packages/orthos) orthos decomposes RNA-seq contrasts,
  for example obtained from a gene knock-out or compound treatment
  experiment, into unspecific and experiment-specific components.
  Original and decomposed contrasts can be efficiently queried
  against a large database of contrasts (derived from ARCHS4,
  https://maayanlab.cloud/archs4/) to identify similar experiments.
  orthos furthermore provides plotting functions to visualize the
  results of such a search for similar contrasts.

- [partCNV](/packages/partCNV) This package uses a statistical
  framework for rapid and accurate detection of aneuploid cells with
  local copy number deletion or amplification. Our method uses an EM
  algorithm with mixtures of Poisson distributions while
  incorporating cytogenetics information (e.g., regional deletion or
  amplification) to guide the classification (partCNV). When
  applicable, we further improve the accuracy by integrating a Hidden
  Markov Model for feature selection (partCNVH).

- [phantasusLite](/packages/phantasusLite) PhantasusLite – a
  lightweight package with helper functions of general interest
  extracted from phantasus package. In parituclar it simplifies
  working with public RNA-seq datasets from GEO by providing access
  to the remote HSDS repository with the precomputed gene counts from
  ARCHS4 and DEE2 projects.

- [plasmut](/packages/plasmut) A Bayesian method for quantifying the
  liklihood that a given plasma mutation arises from clonal
  hematopoesis or the underlying tumor. It requires sequencing data
  of the mutation in plasma and white blood cells with the number of
  distinct and mutant reads in both tissues. We implement a Monte
  Carlo importance sampling method to assess the likelihood that a
  mutation arises from the tumor relative to non-tumor origin.

- [plyinteractions](/packages/plyinteractions) Operate on
  GInteractions objects as tabular data using dplyr-like verbs.
  The functions and methods in plyinteractions provide a
  grammatical approach to manipulate GInteractions, to facilitate
  their integration in genomic analysis workflows.

- [QTLExperiment](/packages/QTLExperiment) QLTExperiment defines an
  S4 class for storing and manipulating summary statistics from QTL
  mapping experiments in one or more states. It is based on the
  'SummarizedExperiment' class and contains functions for creating,
  merging, and subsetting objects. 'QTLExperiment' also stores
  experiment metadata and has checks in place to ensure that
  transformations apply correctly.

- [raer](/packages/raer) Toolkit for identification and statistical
  testing of RNA editing signals from within R. Provides support for
  identifying sites from bulk-RNA and single cell RNA-seq datasets,
  and general methods for extraction of allelic read counts from
  alignment files. Facilitates annotation and exploratory analysis of
  editing signals using Bioconductor packages and resources.

- [RAIDS](/packages/RAIDS) This package implements specialized
  algorithms that enable genetic ancestry inference from various
  cancer sequences sources (RNA, Exome and Whole-Genome sequences).
  This package also implements a simulation algorithm that generates
  synthetic cancer-derived data. This code and analysis pipeline was
  designed and developed for the following publication: Belleau, P et
  al. Genetic Ancestry Inference from Cancer-Derived Molecular Data
  across Genomic and Transcriptomic Platforms. Cancer Res 1 January
  2023; 83 (1): 49–58.

- [regionalpcs](/packages/regionalpcs) Functions to summarize DNA
  methylation data using regional principal components. Regional
  principal components are computed using principal components
  analysis within genomic regions to summarize the variability in
  methylation levels across CpGs. The number of principal components
  is chosen using either the Marcenko-Pasteur or Gavish-Donoho method
  to identify relevant signal in the data.

- [RegionalST](/packages/RegionalST) This package analyze spatial
  transcriptomics data through cross-regional analysis. It selects
  regions of interest (ROIs) and identifys cross-regional cell
  type-specific differential signals. The ROIs can be selected using
  automatic algorithm or through manual selection. It facilitates
  manual selection of ROIs using a shiny application.

- [RNAseqCovarImpute](/packages/RNAseqCovarImpute) The
  RNAseqCovarImpute package implements multiple imputation of missing
  covariates and differential gene expression analysis by: 1)
  Randomly binning genes into smaller groups, 2) Creating M imputed
  datasets separately within each bin, where the imputation predictor
  matrix includes all covariates and the log counts per million (CPM)
  for the genes within each bin, 3) Estimating gene expression
  changes using voom followed by lmFit functions, separately on each
  M imputed dataset within each gene bin, 4) Un-binning the gene sets
  and stacking the M sets of model results before applying the
  squeezeVar function to apply a variance shrinking Bayesian
  procedure to each M set of model results, 5) Pooling the results
  with Rubins’ rules to produce combined coefficients, standard
  errors, and P-values, and 6) Adjusting P-values for multiplicity to
  account for false discovery rate (FDR).

- [roastgsa](/packages/roastgsa) This package implements a variety of
  functions useful for gene set analysis using rotations to
  approximate the null distribution. It contributes with the
  implementation of seven test statistic scores that can be used with
  different goals and interpretations. Several functions are
  available to complement the statistical results with graphical
  representations.

- [Rvisdiff](/packages/Rvisdiff) Creates a muti-graph web page which
  allows the interactive exploration of differential expression
  results. The graphical web interface presents results as a table
  which is integrated with five interactive graphs: MA-plot, volcano
  plot, box plot, lines plot and cluster heatmap. Graphical aspect
  and information represented in the graphs can be customized by
  means of user controls. Final graphics can be exported as PNG
  format.

- [SARC](/packages/SARC) Imports a cov/coverage file (normalised read
  coverages from BAM files) and a cnv file (list of CNVs - similiar
  to a BED file) from WES/ WGS CNV (copy number variation) detection
  pipelines and utilises several metrics to weigh the likelihood of a
  sample containing a detected CNV being a true CNV or a false
  positive. Highly useful for diagnostic testing to filter out false
  positives to provide clinicians with fewer variants to interpret.
  SARC uniquely only used cov and csv (similiar to BED file) files
  which are the common CNV pipeline calling filetypes, and can be
  used as to supplement the Interactive Genome Browser (IGV) to
  generate many figures automatedly, which can be especially helpful
  in large cohorts with 100s-1000s of patients.

- [scDesign3](/packages/scDesign3) We present a statistical
  simulator, scDesign3, to generate realistic single-cell and spatial
  omics data, including various cell states, experimental designs,
  and feature modalities, by learning interpretable parameters from
  real data. Using a unified probabilistic model for single-cell and
  spatial omics data, scDesign3 infers biologically meaningful
  parameters; assesses the goodness-of-fit of inferred cell clusters,
  trajectories, and spatial locations; and generates in silico
  negative and positive controls for benchmarking computational
  tools.

- [scider](/packages/scider) scider is an user-friendly R package
  providing functions to model the global density of cells in a slide
  of spatial transcriptomics data. All functions in the package are
  built based on the SpatialExperiment object, allowing integration
  into various spatial transcriptomics-related packages from
  Bioconductor. After modelling density, the package allows for
  serveral downstream analysis, including colocalization analysis,
  boundary detection analysis and differential density analysis.

- [simona](/packages/simona) The package implements a rich set of
  methods for semantic similarity analysis on bio-ontologies. They
  include methods for information contents, similarities between two
  terms as well as similarities between two groups of terms. It also
  implements visualizations on DAGs.

- [tadar](/packages/tadar) This package provides functions to
  standardise the analysis of Differential Allelic Representation
  (DAR). DAR compromises the integrity of Differential Expression
  analysis results as it can bias expression, influencing the
  classification of genes (or transcripts) as being differentially
  expressed. DAR analysis results in an easy-to-interpret value
  between 0 and 1 for each genetic feature of interest, where 0
  represents identical allelic representation and 1 represents
  complete diversity. This metric can be used to identify features
  prone to false-positive calls in Differential Expression analysis,
  and can be leveraged with statistical methods to alleviate the
  impact of such artefacts on RNA-seq data.

- [TSAR](/packages/TSAR) This package automates analysis workflow for
  Thermal Shift Analysis (TSAS) data. Processing, analyzing, and
  visualizing data through both shiny applications and command lines.
  Package aims to simplify data analysis and offer front to end
  workflow, from raw data to multiple trial analysis.


New Data Experiment Packages
=====================

There are 10 new data experiment packages in this release of Bioconductor.

- [cfToolsData](/packages/cfToolsData) The cfToolsData package
  supplies the data for the cfTools package. It contains two
  pre-trained deep neural network (DNN) models for the cfSort
  function. Additionally, it includes the shape parameters of beta
  distribution characterizing methylation markers associated with
  four tumor types for the CancerDetector function, as well as the
  parameters characterizing methylation markers specific to 29
  primary human tissue types for the cfDeconvolve function.

- [gDRtestData](/packages/gDRtestData) R package with internal
  dose-response test data. Package provides functions to generate
  input testing data that can be used as the input for gDR pipeline.
  It also contains RDS files with MAE data processed by gDR.

- [HCATonsilData](/packages/HCATonsilData) This package provides
  access to the scRNA-seq, scATAC-seq, multiome, CITE-seq and spatial
  transcriptomics (Visium) data generated by the tonsil cell atlas in
  the context of the Human Cell Atlas (HCA). The data is provided via
  the Bioconductor project in the form of SingleCellExperiments.
  Additionally, information on the whole compendium of identified
  cell types is provided in form of a glossary.

- [HiBED](/packages/HiBED) Hierarchical deconvolution for extensive
  cell type resolution in the human brain using DNA methylation. The
  HiBED deconvolution estimates proportions up to 7 cell types
  (GABAergic neurons, glutamatergic neurons, astrocytes, microglial
  cells, oligodendrocytes, endothelial cells, and stromal cells) in
  bulk brain tissues.

- [multiWGCNAdata](/packages/multiWGCNAdata) Stores expression
  profiling data from experiments compatible with the multiWGCNA R
  package. This includes human postmortem microarray data from
  patients and controls (GSE28521), astrocyte Ribotag RNA-seq data
  from EAE and wildtype mice (GSE100329), and mouse RNA-seq data from
  tau pathology (rTg4510) and wildtype control mice (GSE125957).
  These data can be accessed using the ExperimentHub workflow (see
  multiWGCNA vignettes).

- [orthosData](/packages/orthosData) orthosData is the companion
  ExperimentData package to the orthos R package for mechanistic
  studies using differential gene expression experiments. It provides
  functions for retrieval from ExperimentHub and local caching of the
  models and datasets used internally in orthos.

- [raerdata](/packages/raerdata) raerdata is an ExperimentHub package
  that provides a collection of files useful for demostrating
  functionality in the raer package. Datasets include 10x genomics
  scRNA-seq, bulk RNA-seq, and paired whole-genome and RNA-seq data.
  Additionally databases of human and mouse RNA editing sites are
  provided.

- [smokingMouse](/packages/smokingMouse) This is an ExperimentHub
  package that provides access to the data at the gene, exon,
  transcript and junction level used in the analyses of the
  smokingMouse project. See
  https://github.com/LieberInstitute/smokingMouse_Indirects. This
  datasets contain the expression counts of genes, transcripts, exons
  and exon-exon junctions across 208 mice samples from pup and adult
  brains and adult blood. They also contain relevant information of
  these samples and features, such as conditions, QC metrics and if
  they were used after filtering steps and also if the features were
  differently expressed in the different experiments.

- [SpatialDatasets](/packages/SpatialDatasets) This is a collection
  of publically available spatial omics datasets. Where possible we
  have curated these datasets as either SpatialExperiments,
  MoleculeExperiments or CytoImageLists and included annotations of
  the sample characteristics.

- [TumourMethData](/packages/TumourMethData) TumourMethData collects
  tumour methylation data from a variety of different tumour types
  (and also matching normal samples where available) and produced
  with different technologies (e.g. WGBS, RRBS and methylation
  arrays) and provides them as RangedSummarizedExperiments. This
  facilitates easy extraction of methylation data for regions of
  interest across different tumour types and studies.

New Annotation Packages
=====================

There are 8 new annotation packages.

- [AlphaMissense.v2023.hg19](/packages/AlphaMissense.v2023.hg19)
  Store Google DeepMind AlphaMissense v2023 hg19 pathogenicity scores
  AnnotationHub Resource Metadata. Provide provenance and citation
  information for Google DeepMind AlphaMissense v2023 hg19
  pathogenicity score AnnotationHub resources. Illustrate in a
  vignette how to access those resources.

- [AlphaMissense.v2023.hg38](/packages/AlphaMissense.v2023.hg38)
  Store Google DeepMind AlphaMissense v2023 hg38 pathogenicity scores
  AnnotationHub Resource Metadata. Provide provenance and citation
  information for Google DeepMind AlphaMissense v2023 hg38
  pathogenicity score AnnotationHub resources. Illustrate in a
  vignette how to access those resources.

- [cadd.v1.6.hg19](/packages/cadd.v1.6.hg19) Store University of
  Washington CADD v1.6 hg19 pathogenicity scores AnnotationHub
  Resource Metadata. Provide provenance and citation information for
  University of Washington CADD v1.6 hg19 pathogenicity score
  AnnotationHub resources. Illustrate in a vignette how to access
  those resources.

- [cadd.v1.6.hg38](/packages/cadd.v1.6.hg38) Store University of
  Washington CADD v1.6 hg38 pathogenicity scores AnnotationHub
  Resource Metadata. Provide provenance and citation information for
  University of Washington CADD v1.6 hg38 pathogenicity score
  AnnotationHub resources. Illustrate in a vignette how to access
  those resources.

- [HPO.db](/packages/HPO.db) Human Phenotype Ontology (HPO) was
  developed to create a consistent description of gene products with
  disease perspectives, and is essential for supporting functional
  genomics in disease context. Accurate disease descriptions can
  discover new relationships between genes and disease, and new
  functions for previous uncharacteried genes and alleles.We have
  developed the [DOSE](https://bioconductor.org/packages/DOSE/)
  package for semantic similarity analysis and disease enrichment
  analysis, and DOSE import an Bioconductor package DO.db to get
  the relationship(such as parent and child) between MPO terms. But
  DO.db hasn't been updated for years, and a lot of semantic
  information is
  [missing](https://github.com/YuLab-SMU/DOSE/issues/57). So we
  developed the new package HPO.db for Human Human Phenotype
  Ontology annotation.

- [JASPAR2024](/packages/JASPAR2024) JASPAR
  (https://testjaspar.uio.no/) is a widely-used open-access database
  presenting manually curated high-quality and non-redundant
  DNA-binding profiles for transcription factors (TFs) across taxa.
  In this 10th release and 20th-anniversary update, the CORE
  collection has expanded with 329 new profiles. We updated three
  existing profiles and provided orthogonal support for 72 profiles
  from the previous release UNVALIDATED collection. Altogether, the
  JASPAR 2024 update provides a 20 percent increase in CORE profiles
  from the previous release. A trimming algorithm enhanced profiles
  by removing low information content flanking base pairs, which were
  likely uninformative (within the capacity of the PFM models) for
  TFBS predictions and modelling TF-DNA interactions. This release
  includes enhanced metadata, featuring a refined classification for
  plant TFs structural DNA-binding domains. The new JASPAR
  collections prompt updates to the genomic tracks of predicted
  TF-binding sites in 8 organisms, with human and mouse tracks
  available as native tracks in the UCSC Genome browser. All data are
  available through the JASPAR web interface and programmatically
  through its API and the updated Bioconductor and pyJASPAR packages.
  Finally, a new TFBS extraction tool enables users to retrieve
  predicted JASPAR TFBSs intersecting their genomic regions of
  interest.

- [MPO.db](/packages/MPO.db) We have developed the human disease
  ontology R package HDO.db, which provides the semantic relationship
  between human diseases. Relying on the DOSE and GOSemSim packages
  we developed, we can carry out disease enrichment and semantic
  similarity analyses. Many biological studies are achieved through
  mouse models, and a large number of data indicate the association
  between genotypes and phenotypes or diseases.  The study of model
  organisms can be transformed into useful knowledge about normal
  human biology and disease to facilitate treatment and early
  screening for diseases. Organism-specific genotype-phenotypic
  associations can be applied to cross-species phenotypic studies to
  clarify previously unknown phenotypic connections in other species.
  Using the same principle to diseases can identify genetic
  associations and even help to identify disease associations that
  are not obvious. Therefore, as a supplement to HDO.db and DOSE, we
  developed mouse phenotypic ontology R package MPO.db.

- [SomaScan.db](/packages/SomaScan.db) An R package providing
  extended biological annotations for the SomaScan Assay, a
  proteomics platform developed by SomaLogic Operating Co., Inc. The
  annotations in this package were assembled using data from public
  repositories. For more information about the SomaScan assay and its
  data, please reference the 'SomaLogic/SomaLogic-Data' GitHub
  repository.

New Workflow Packages
=====================

There are no new workflow packages in this release of Bioconductor.

New Online Books
=====================

There is one new online book.

- [BiocBookDemo](https://bioconductor.org/books/3.18/BiocBookDemo/) This package has been
  created using the BiocBook package. It serves as a demo of a
  BiocBook online book. Read BiocBook package documentation to
  know more about BiocBooks.

NEWS from existing Software Packages
===================================


[adverSCarial](/packages/adverSCarial)
------------

                Changes in version 0.99.54 (2023-10-22)                 

- Ready for production

                 Changes in version 0.99.1 (2023-04-05)                 

- Modify package to remove warnings and notes from BiocCheck

                 Changes in version 0.99.0 (2023-02-05)                 

- Submitted to Bioconductor

[affxparser](/packages/affxparser)
----------

                 Changes in version 1.73.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor develop version,
which is now Bioconductor 3.18 for R (>= 4.4.0).

[ALDEx2](/packages/ALDEx2)
------

                        Changes in version 1.33                         

NEW FEATURES (MPN, GBG

- aldex.clr: now takes a gamma parameter to incorporate scale modelling

- aldex.makeScaleMatrix: new method to make an explicit scale model

- all p-values calculated are now posterior p-values with consistent
  sign

[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.14.0                        

NEW FEATURES

- (v 1.13.1) Add paged support for large tables in avtable_import()
and avtable_import_set().

- (v 1.13.2) Only show avtable_paged() and avtable_import*() progress
bar in interactive() sessions

- (v 1.13.4) Report messages when avtable_import_status() contains
one. https://github.com/Bioconductor/AnVIL/issues/79

- (v 1.13.3) Use 'op' when .avworkflow_response() calls
avstop_for_status(). https://github.com/Bioconductor/AnVIL/issues/80

- (v 1.13.7) Check requester pays for destination URIs when using
gsutil_cp (@smgogarten, #82)

USER VISIBLE CHANGES

- (v 1.13.8) Update documentation on updating workflow
configurations.
(@amstilp, #84)

[AnVILPublish](/packages/AnVILPublish)
------------

                       Changes in version 1.12.0                        

Bug Fixes

- (v. 1.11.3) Link Rmd (and / or ipynb) vignettes to dashboard

[apeglm](/packages/apeglm)
------

                       Changes in version 1.23.1                        

- Sometimes prefit.beta was NaN, now properly handled in apeglm().
  Resolves some "Error in optimHess" issue propagated to DESeq2.

[aroma.light](/packages/aroma.light)
-----------

                 Changes in version 3.31.1 (2023-06-30)                 

Documentation

- Update redirecting and broken URLs.

- Fix R CMD check notes on "Escaped LaTeX specials: &".

                 Changes in version 3.31.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor devel version,
which is now Bioconductor 3.18 for R (>= 4.4.0).

[ASpli](/packages/ASpli)
-----

                       Changes in version 2.11.1                        

BUG FIXES

- Solve bug in readCounts. Undefined seqnames produced invalid
  junctions.orders

- Change package mantainer

                       Changes in version 2.10.1                        

BUG FIXES

- Solves jCounts() Error in av[at] <- a[at] : NAs are not allowed in
  subscripted assignments bug with large junction possitions.

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.25.2                        

- Add reNormalizeByDistalSig parameter to plotFootprints function.

                       Changes in version 1.25.1                        

- Fix the error 'number of columns of matrices must match' for
TSSEscore.

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.13.3 (2023-06-16)                 

- Reverts `LazyData` to `false` based on BiocCheck feedback

                 Changes in version 2.13.2 (2023-06-16)                 

- Changes `LazyData` to `true` to avoid errors

                 Changes in version 2.13.1 (2023-06-16)                 

- Fixes error in tests associated to the behaviour of
  `matrixStats::colMedians`

[benchdamic](/packages/benchdamic)
----------

                 Changes in version 1.7.5 (2023-10-13)                  

- Updated DA_Seurat() for the new version

- Waiting for bug-fix: DA_ANCOM() (random effects not working)

- Waiting for bug-fix: DA_mixMC() (multilevel analysis not working)

- Bug-fix: correct dataset name in unit tests

- Bug-fix: in DA_ALDEx2() replaced unlist() with as.vector()

- Minor documentation updates related to dependencies

                 Changes in version 1.7.4 (2023-07-22)                  

- Reducing the number of comparisons in vignette to reduce build time

                 Changes in version 1.7.3 (2023-07-11)                  

- Bug-fix: Minor changes in vignette

                 Changes in version 1.7.2 (2023-07-09)                  

- Bug-fix for DA_Maaslin2() function

                 Changes in version 1.7.1 (2023-07-07)                  

- New methods: linDA, Maaslin2, ZicoSeq

- New method: mixMC (temporarily unavailable)

- Add 'alpha' parameter in DA_ANCOM()

                 Changes in version 1.7.0 (2023-04-25)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_18 branch

                 Changes in version 1.6.4 (2023-07-22)                  

- Porting the changes of devel version 1.7.4 to release

                 Changes in version 1.6.3 (2023-07-11)                  

- Porting the changes of devel version 1.7.3 to release

                 Changes in version 1.6.2 (2023-07-09)                  

- Porting the changes of devel version 1.7.2 to release

                 Changes in version 1.6.1 (2023-07-08)                  

- Porting the changes of devel version 1.7.1 to release

[BindingSiteFinder](/packages/BindingSiteFinder)
-----------------

                       Changes in version 1.7.12                        

- Updated the vignette to feature differential binding analysis

                       Changes in version 1.7.11                        

- assignToGenes() no requires binding sites to fully overlapp the gene
  range to be
  assigned to that gene

- added plotBsMA() and plotBsVolcano() to visualize differential
  binding output

                       Changes in version 1.7.10                        

- added calculateBsBackground() and calculateBsFoldChange() core
  function for
  differntial binding analysis

                        Changes in version 1.7.9                        

- fixed a bug in combineBSF() where meta data was not correctly merged

                        Changes in version 1.7.8                        

- added combineBSF() function to combine two or more object of type
  BSFDataSet

                        Changes in version 1.7.7                        

- Fixed a bug in exportToBED() which caused the export to fail if the
  last function
  exectuted was calculateSignalToFlankScore()

- Updated rangeCoveragePlot() to work with clipCoverage() function

- Added the clipCoverage() function as improoved function calculate
  coverage

                        Changes in version 1.7.6                        

- Update assignToGenes() to not require a gene annotation when options
  'remove'
  or 'keep' is selected

                        Changes in version 1.7.5                        

- Update vignette to include all new options

- Exchanged default test object to fit new class definition

- Let estimateBsWidth() fail more gracefull when no maximum can be
  found

- Fix name space dependencies

                        Changes in version 1.7.4                        

- Added region length based normalization to
  transcriptRegionSpectrumPlot()

- Changed minWidth default from (3 -> 2) in makeBindingSites()

- Allow object subsetByChr() to handle multiple chromosomes

- Added calculateSignalToFlankScore() function

- Added binding site definedness plot

- Added a 'local' version to estimateBsWidth() if no maximum can be
  found on
  global level

- Added a 'sensitivity' mode to estimateBsWidth()

                        Changes in version 1.7.3                        

- Started major rework

- Restructuring of class definition

- Added BSFind() as core function

- Added workflow functions pureClipGlobalFilter(), estimateBsWidth(),
  pureClipGeneWiseFilter(), assignToGenes(),
  assignToTranscriptRegions()

- Reworked makeBindingSites(), reproducibilityFilter(),
  annotateWithScore()

- Added plotting functions processingStepsFlowChart(),
  pureClipGlobalFilterPlot(),
  estimateBsWidthPlot(), duplicatedSitesPlot(),
  mergeCrosslinkDiagnosticsPlot(),
  makeBsSummaryPlot(), reproducibilityFilterPlot(),
  reproducibilitySamplesPlot(),
  reproducibilityScatterPlot(), geneOverlapsPlot(),
  targetGeneSpectrumPlot(),
  transcriptRegionOverlapsPlot(), transcriptRegionSpectrumPlot(),
  bindingSiteDefinednessPlot()

                        Changes in version 1.7.2                        

- Added further input checks to reproducibilityFilter() function

                        Changes in version 1.7.1                        

- Fix Namesspace issues

                        Changes in version 1.6.1                        

- Fix bugs in colorPalette option

- Added custom coloring

[bioCancer](/packages/bioCancer)
---------

                       Changes in version 1.29.05                       

- Comment pickGO during running examples: It works manually.

                       Changes in version 1.29.01                       

- Correct items format in NEWS file

- Import needed packages

- Add @method section in documentation of cgdsr methods

- Update links and reference in the vignette

[BiocBaseUtils](/packages/BiocBaseUtils)
-------------

                        Changes in version 1.4.0                        

New features

- Added isScalarLogical for completeness; identical to isTRUEorFALSE.

[BiocCheck](/packages/BiocCheck)
---------

                       Changes in version 1.38.0                        

BUG FIXES AND MINOR IMPROVEMENTS

- Exclude data docs with `\\format` tags in addition to 'package' docs
  when
  checking for `\\value` / `@return` in documentation.

- Resolve unknown macro warnings when using Rdpack (@LiNK-NY, #196)

- Improve the `read.dcf` operation that looks for any deprecated
  packages
  in both Bioconductor release and devel versions.

- Remove overwrite prompt when updating cached resources in deprecated
  packages check.

[BiocFileCache](/packages/BiocFileCache)
-------------

                         Changes in version 2.9                         

USER VISIBLE CHANGE

- (2.9.1) Add documentation for operating behind a proxy

[BiocHubsShiny](/packages/BiocHubsShiny)
-------------

                        Changes in version 1.2.0                        

Bug fixes and minor improvements

- Updated the NEWS.md file formatting

[BioCor](/packages/BioCor)
------

                        Changes in version 1.26                         

- Added support for plots of similarities

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.20.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- `biocRevDepEmail` accepts a vector of deprecated packages for whose
  reverse dependencies are to be notified.

- `biocDownloadStats` now filters out packages that are not part of the
  `pkgType` option and strictly returns stats for Bioconductor packages

- When package stats are not available for a particular year,
  `biocDownloadStats` will warn about the year there is no data for.

BUG FIXES

- Warn when a download stats URL is not available and filter out in
  `pkgDownloadStats`.

- `pkgDownloadRank` calculates the denominator for ranks using the
  `VIEWS`
  file; matches Bioc badges more closely.

[BiocStyle](/packages/BiocStyle)
---------

                       Changes in version 2.30.0                        

BUG FIXES

- Addressed issue where code chunks in HTML output did not have the
  correct background colour

- Fixed incompatibility with ragged2e LaTeX package distributed in
  TexLive 2023
  (https://github.com/Bioconductor/BiocStyle/issues/105)

[biocthis](/packages/biocthis)
--------

                       Changes in version 1.11.3                        

BUG FIXES

- Fixed internal code on use_bioc_news_md(), use_bioc_readme_rmd(),
and use_bioc_vignette() after usethis changed some of their
un-exported functions that biocthis relies on.

                       Changes in version 1.11.2                        

BUG FIXES

- Ignore remotes::system_requirements("ubuntu", "20.04") for now
since
that leads to a JSON error. See
https://github.com/lcolladotor/biocthis/issues/41 but also
https://github.com/LieberInstitute/spatialLIBD/commit/edc8b72505af097895dcbf35887df28da8122e3c.

                       Changes in version 1.11.1                        

BUG FIXES

- Resolved https://github.com/lcolladotor/biocthis/issues/40 reported
by @lmweber with help from @bschilder and @LiNk-NY noted at
https://github.com/neurogenomics/rworkflows/issues/58. Basically,
there's no longer a need to explicitly list the AnVIL repositories
to benefit from those gains.

[biodbHmdb](/packages/biodbHmdb)
---------

                 Changes in version 1.7.1 (2023-09-20)                  

- Fix truncation of HMDB Metabolites XML file.

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.58.0                        

USER VISIBLE CHANGES

- getSequence() will now provide a more informative error message if
  requesting a flanking sequence and not provided with an upstream or
  downstream
  range.

- Remove references to the uswest mirror, which has now been retired
  (https://www.ensembl.info/2023/01/13/retirement-of-ensembl-us-west-aws-mirror/)

[BioNAR](/packages/BioNAR)
------

                         Changes in version 1.3                         

- Take into account edge weights in clustering algorithms and
centrality measure calculations.
- Add calculation of the DYNAMO perturbation pattern from signed
weight directed networks proposed in Santolini,M. and Barabasi,A.-L.
(2018) PNAS 169, 201720589
- Decoupled from synaptome.db and synaptome.data packages. All code,
related to graph building from the synaptome.db data is moved to
synaptome.db package.

                        Changes in version 1.2.1                        

- Allow analysis of directed graphs and add four new centrality
measures specific for the directed graphs.
- Modify annotation functions to allow annotate nodes not only by its
name propety but by any other arbitrary property if its value
uniquely identify the node.
- Added new framework for the visualisation and analysis of
enrichment
results.

[biosigner](/packages/biosigner)
---------

                       Changes in version 1.29.2                        

MINOR MODIFICATION

- minor update of the 'show' method

[Biostrings](/packages/Biostrings)
----------

                       Changes in version 2.70.0                        

NEW FEATURES

- Character set of AAString/AAStringSet/AAStringSetList objects is now
  enforced (a long-due feature). Thanks to Aidan Lakshman
  <ahl27@pitt.edu>
  for implementing this.

[BridgeDbR](/packages/BridgeDbR)
---------

                       Changes in version 2.11.2                        

BUG FIXES

- Fixed build/test process by including a long-running example into
the dontrun{} environment

                       Changes in version 2.11.1                        

NEW FEATURES

- Migrated NEWS to NEWS.md
- Applied BioC code styling
- Added templates for issues and feature requests
- Added a method exists(mapper, xref) that checks if the xref is
found
in the mapping file

SIGNIFICANT USER-VISIBLE CHANGES

- loadDatabase() now reports an error if the given location does not
exist

BUG FIXES

- Fixed backwards compatibility of map()
- Clarified in the vignette the output getMatchingSources()
- Fixed returning an empty data frame when no mappings are found for
map()

[BSgenomeForge](/packages/BSgenomeForge)
-------------

                        Changes in version 1.2.0                        

NEW FEATURES

- Two improvements to forgeBSgenomeDataPkgFromNCBI() when the NCBI
  assembly is registered in the GenomeInfoDb package:
  - The assembly name can be specified via the 'assembly_accession'
  argument instead of its GenBank or RefSeq accession.
  - The 'organism' argument no longer needs to be specified.

[CaDrA](/packages/CaDrA)
-----

                       Changes in version 0.99.1                        

- Published in https://github.com/montilab/CaDrA
- CaDrA-Shiny repo for web interface:
https://github.com/montilab/CaDrA-shiny
- Web Interface: https://cadra.bu.edu/

[CAGEr](/packages/CAGEr)
-----

                        Changes in version 2.8.0                        

BUG FIXES

- Correct quantile positions, which were shifted by one base. This
bug
may have been introduced in version 1.22 or later.
- Ensure cluster objects are properly sorted. Fixes #79, introduced
in
version 2.6.0 and causing crashes or incorrect quantile
calculations.
- Apply fix for #77 (aggregateTagClusters losing TCs), which slipped
out of 2.6.0 because of Git branch mixup.
- Fix consensus cluster coordinates, where the maxDist padding was
erroneously remaining in some parts of the computation.
- Corrected on-the-fly cumulative sum computation for consensus
clusters when sample = NULL. The bug was causing incorrectly short
quantile ranges.
- Force the cluster names to stay sorted, to avoid a bug
desynchronising quantile information and genome coordinates.
- Fix accidental deletion of the record of the tag clustering method
(#96).

NEW FEATURES

- Allow URLs to files in getCTSS() (Fixes #50).
- Accelerated the computation of quantile position by ~20 times.
- New resetCAGEexp() function.
- New flagByUpstreamSequences() function.
- The annotateCTSS and annotateConsensusClusters function gain a
upstream and a downstream parameter to change the width of promoter
regions.

[Cardinal](/packages/Cardinal)
--------

                 Changes in version 3.3.5 (2023-10-19)                  

BUG FIXES

- Allow 'guess.max=Inf' in 'readImzML()'

                 Changes in version 3.3.4 (2023-10-19)                  

SIGNIFICANT USER-VISIBLE CHANGES

- Bin peaks when importing centroid spectra with 'readImzML()'

- Replaced all functions from deprecated 'sp' package

                        Changes in version 3.3.3                        

BUG FIXES

- Fixed I/O bugs introduced by matter v2.3.13 changes

                        Changes in version 3.3.2                        

BUG FIXES

- Fixed I/O bugs introduced by matter v2.3.11 changes

- Other I/O bugs fixed by matter v2.3.13 changes

                        Changes in version 3.3.1                        

BUG FIXES

- Merged 3.2.1 fixes to resolve R CMD check warnings

                        Changes in version 3.2.1                        

BUG FIXES

- Cleaned up escaped LaTeX specials in documentation

- Fixed 'sprintf()' => 'snprintf()' warning in C code

[CardinalIO](/packages/CardinalIO)
----------

                       Changes in version 0.99.2                        

NEW FEATURES

- Added support for additional binary data arrays

SIGNIFICANT USER-VISIBLE CHANGES

- Casefold checksums and UUIDs before comparison

                       Changes in version 0.99.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Added package-level 'CardinalIO' help page

                       Changes in version 0.99.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Update 'Analyze75' and 'ImzML' classes for Bioconductor

- Prepare for Bioconductor submission


[cBioPortalData](/packages/cBioPortalData)
--------------

                       Changes in version 2.14.0                        

Bug fixes and minor improvements

- Remove error when studyId build success is unknown (@vlaufer, #69)
- Parse tokens with different formats, e.g. genie token (@ZWael, #70)
- Clean up duplicate rows created in the sampleMap generated from the
data
- Fixed warning trigger when empty molecularProfileId datasets are
found

[CCPlotR](/packages/CCPlotR)
-------

                 Changes in version 0.99.3 (2023-10-06)                 

- Made changes to address Bioconductor reviews

- added unit tests

- added validity checks and helpful error messages for function
  arguments

- set Lazydata to false

                 Changes in version 0.99.2 (2023-09-19)                 

- Added extra optional paramaters to cc_circos and cc_network functions

                 Changes in version 0.99.1 (2023-08-18)                 

- Submitted to Bioconductor

[CellBarcode](/packages/CellBarcode)
-----------

                        Changes in version 1.7.0                        

New features

- bc_extract_sc_sam: function to extract single-cell barcodes from
SAM
file
- bc_extract_sc_fastq: function to extract single-cell barcodes from
fastq file
- bc_plot_count: function to plot the count feature of the barcodes
such as barcode reads versus total reads, reads per UMI
distribution.

Incompatible changes

bc_extract now only accepts a vector or list, and returns BarcodeObj,
instead of returning data.frame when the input is a single sample.

[cellxgenedp](/packages/cellxgenedp)
-----------

                         Changes in version 1.6                         

- (v 1.5.2) use CELLxGENE 'Discover' API, changing column names of
some return values. See 'API changes' of the 'Discover and download
datasets...' vignette.

[cfdnakit](/packages/cfdnakit)
--------

                 Changes in version 0.99.4 (2023-05-28)                 

- fix required issues and suggesstion except code formatting including
  lien space and line length.

                 Changes in version 0.99.0 (2022-11-09)                 

- Submitted to Bioconductor

[cfTools](/packages/cfTools)
-------

                        Changes in version 1.2.0                        

NEW FEATURES

- Add the `cfSort()` function, which is the first supervised
  tissue deconvolution approach with deep learning models.

- Publish the `cfToolsData` package to supply models and
  marker files for `cfTools`.

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.35.3                        

- update documentation for the annotatePeakInBatch function

                       Changes in version 3.35.1                        

- fix the seqlevelsStyle check for custom species.

[circRNAprofiler](/packages/circRNAprofiler)
---------------

                       Changes in version 1.15.2                        

- Possibility to use a maximum of 3 unsupported/‘other’ tools.

[ClassifyR](/packages/ClassifyR)
---------

                        Changes in version 3.6.0                        

- 
  Balancing of non-censored event times across folds.

- 
  Multiview methods with combinations of views have improved
  plots that switch to UpSet axis.

- 
  Multiview methods have feature extractors to retain chosen
  features for analysis.

[clusterProfiler](/packages/clusterProfiler)
---------------

                        Changes in version 4.9.5                        

- fixed R check (2023-10-18, Wed)

                        Changes in version 4.9.4                        

- use check_installed() to check package dependency (2023-09-08, Fri,
#621)
- use yread() in WikiPathway utilities (2023-09-07, Thu)

                        Changes in version 4.9.3                        

- enrichKEGG() and gseKEGG() now supports organism = 'cpd' to accept
KEGG Compound ID (2023-08-31, Thu)
- gson_cpd() and gson_ko()
- use yulab.utils::yread() to parse file (2023-08-15, Tue)
- supports Pathways Common (2023-08-02, Wed, #613)

                        Changes in version 4.9.2                        

- append_kegg_category() function to add KEGG pathway category
information to KEGG enrichment result and now it is the default
behavior of enrichKEGG() and gseKEGG() (2023-07-12, Wed)
- parse KEGG Pathway Category information (2023-07-11, Tue)
- mv parse_gff() to GOSemSim::read.gaf() and re-export (2023-07-10,
Mon)
- mv buildGOmap() to `GOSemSim::buildGOmap() and re-export

                        Changes in version 4.9.1                        

- getPPI() to query PPI network from 'stringdb' (2023-05-15, Mon)
- getTaxID() and getTaxInfo() functions to query taxonomy information
(2023-05-14, Sun)

[ClustIRR](/packages/ClustIRR)
--------

                     Changes in version 0.99.28-29                      

- Accepted in Bioconductor

                   Changes in version 0.99.20-0.99.26                   

- Functions added for: graph building, joining, plotting

- S4 object accessors added

                       Changes in version 0.99.16                       

- Changes after first round of review done

- Now we have two versions (old version=1 removed, unused)

- S4 object output of clust_irr

                       Changes in version 0.99.4                        

- Build errors/warnings fixed

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[CNVfilteR](/packages/CNVfilteR)
---------

                       Changes in version 1.15.2                        

MINOR

- Fixed bug: CopyNumberPlots requires additional parameter

[CNVMetrics](/packages/CNVMetrics)
----------

                        Changes in version 1.5.1                        

NEW FEATURES

- processSim() method generates simulated samples with copy number
profiles derived from a specific sample.

SIGNIFICANT USER-VISIBLE CHANGES

- None

BUG FIXES

- None

[compcodeR](/packages/compcodeR)
---------

                       Changes in version 1.37.3                        

- Removed support for baySeq since it has been deprecated

- Added runComparisonShiny function

[compSPOT](/packages/compSPOT)
--------

                       Changes in version 0.99.10                       

- Edits made based on reviewer comments

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[COTAN](/packages/COTAN)
-----

                        Changes in version 2.1.8                        

Made passing clusterizations to COTAN functions more easy: now all
functions that take a COTAN object and a clusterization as input
parameters can also take a clusterization name

Added time-stamps to log entries when written on a log file

Fixed bug in the clustersMarkersHeatmapPlot function when given a
clusterization not matching the latest added to the COTAN object

Fixed issue with the highest possible resolution in
seuratClustering()
function, needed when large datasets must be split in many clusters

                        Changes in version 2.1.7                        

Added new flag to the function cleanPlots() to suppress evaluation of
the PCA on the normalized data. In particular, this allows to reduce
significantly time spent within the function checkClusterUniformity()

Added initialResolution parameter to cellsUniformClustering(): it
allows
users to specify the initial resolution used in the calls to
Seurat::FindClusters() method. It now uses the same default as Seurat

Added new method estimateNuLinearByCluster() that calculates nu
ensuring
that its average is 1.0 in each given cluster

                        Changes in version 2.1.6                        

Added function reorderClusterization(): it reorders the given
clusterization so that near clusters have also near labels

The functions cellsUniformClustering() and
mergeUniformCellsClusters()
now return the result of this new function

Separated p-value calculations from DEAOnClusters() into the new
function pValueFromDEA(). Those data.frames are no longer part of the
list returned by the functions DEAOnClusters() and
mergeUniformCellsClusters()

Added function getClusters() to retrieve the wanted clusterization
from
the cells' meta-dataset

Added function calculateGenesCE(): it returns the cross-entropy
between
the expected absence of a gene reading against the observed state

Fixed minor issue with logThis() to file: it was always appending a
new
line even when appendLF was set to FALSE

Now checkClusterUniformity() returns more GDI stats like the
percentage
of genes above threshold or the last percentile of the GDI values

Revamped mergeUniformCellsClusters() to select in order all the the
most
likely candidates pairs of clusters to merge. Provided new user
parameter to balance the merging of most possible candidates versus
the
time spent doing so

Improved dropGenesCells() method: it now retains all meta-data
information that is not related to the results of the other methods

Added zoomed UDE plot to cleanPlots() return. It suggests a possible
cut
level for low UDE cells

                        Changes in version 2.1.5                        

Improved mergeUniformCellsClusters(): now it attempts to merge more
clusters pairs

Now errors in the seuratClustering() function are interpreted as
remaining cells not-clustered "-1". This applies mostly to cases when
Seurat finds only singlets

Added flag calcCoex to proceedToCoex() and
automaticCOTANObjectCreation() functions to allow user not to spend
time
calculating the genes' COEX when not needed

Solved potential issue in the clustersMarkersHeatmapPlot() regarding
clusters' labels

Added new internal function niceFactorLevels() that ensures all the
factors' levels will have labels with the same length, via padding
the
integers values with '0' and string values with '_'

Relaxed tolerance on tests comparing against saved data

                        Changes in version 2.1.4                        

Speed-up by use of parallelDist::parDist() to calculate distances
instead of stats::dist()

Fixed regression tests failing on non-Linux architectures

                        Changes in version 2.1.3                        

Completed function clustersMarkersHeatmapPlot()

Added new utility function normalizeNameAndLabels()

Added mergeClusters() and multiMergeClusters() functions

Added support to conditions in cells' meta-data

Now clusterizations are stored as factors

Fixed COTAN::validity method in AllClasses.R

                        Changes in version 2.1.2                        

Fixed bug in proceedToCoex() in cases when saveObj == TRUE

                        Changes in version 2.1.1                        

Updated README.md and NEWS.md

Renamed methods dealing with housekeeping genes and fully-expressed
cells to use the more proper names fully-expressed genes and
fully-expressing cells

Added possibility to users to set the cutoff and thresholds used by
the
clean and related methods

                        Changes in version 2.1.0                        

First release in Bioconductor 3.18

[CTdata](/packages/CTdata)
------

                         Changes in version 1.1                         

CTdata 1.1.5

- Assay from CT_methylation_in_tissues() converted from a tibble to a
matrix.
- Suggest and load SingleCellExperiment in vignette.
- Update vignette figure.

CTdata 1.1.4

- New data: scRNAseq_HPA and testis_sce.
- Updated data: CT_methylation_in_tissues,
CT_mean_methylation_in_tissues, TCGA_CT_methylation, CT_genes and
CCLE_correlation_matrix.

CTdata 1.1.3

- Load SummarizedExperiment in vignette.

CTdata 1.1.2

- Suggest SummarizedExperiment.

CTdata 1.1.1

- Only display a small subset of CCLE_correlation_matrix() in the
vignette.

CTdata 1.1.0

- New Bioconductor devel release.

[cytomapper](/packages/cytomapper)
----------

                 Changes in version 1.13.1 (2023-09-26)                 

- loadImages: fixed pattern argument to only test on the actual files
  and not the path.

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

[CytoPipeline](/packages/CytoPipeline)
------------

                         Changes in version 1.1                         

CytoPipeline 1.1.5

- areSignalCols(), are FluoCols() can now accept a flowSet as input,
on top of a flowFrame
- applyScaleTransforms() processing step has been improved (can take
flowSet as input, checks channel concordance between transList and
data object)

CytoPipeline 1.1.4

- updated fcs files, that are at the source of OMIP021Samples dataset

CytoPipeline 1.1.3

- in subSample(), renamed parameter 'nSamples' into 'nEvents', and
added possibility for passing unused parameters, in order to support
the use of the function as a processing step. Also amended the
function as to keep the original order of the events (keep
chronology). Finally, adds a 'keepOriginalCellIDs' parameter
(default=TRUE).
- simplified the arguments of execute() related to the storage of the
results after last pre-processing step.

CytoPipeline 1.1.2

- storage of phenoData into cache upon execution of CytoPipeline
object (and back into CytoPipeline object when re-built from cache)
- changed the default behaviour of estimateScaleTransforms() so that
the default method for scatter channels is now "none" instead of
"linearQuantile"
- changed default behaviour of ggplotEvents() and
ggplotFilterEvents(), when logicle scale is used but no logicle
parameters provided, these are now estimated using
flowCore::estimateLogicle(), instead of explicit default values

CytoPipeline 1.1.1

- tiny modifications to support upgrade to Bioc 3.18

CytoPipeline 0.99

CytoPipeline 0.99.6

- corrected the OMIP021Samples fcs data in order to keep the original
file name
- bug correction: error message on execution with no sample file
- added phenoData slot in CytoPipeline object
- updated readSamples() to allow passing a pData parameters
- updated compensateFromMatrix() to allow passing a mapping based on
a
pData variable
- updated readSamples() to allow selecting a random number of samples
and removed selectSamples()
- vignette with demo and links to videos

CytoPipeline 0.99.5

- reactivated unit tests for ggplot2 objects
- added man page for CytoPipeline package
- a few modifs in the vignette related to Bioc review process
- replaced withr::local_tempdir() by base::tempdir()
- removed extraneous whitespaces in CytoPipeline show() method
- removed LazyData: true in DESCRIPTION file
- replaced paste0(path, "/", filename) by file.path(path, filename)
- updated License field in DESCRIPTION file

CytoPipeline 0.99.4

- improved CytoPipeline constructors (experimentName and sampleFiles
are now parameters of all constructor version)
- centralized the production of standard outputs during pipeline
execution, set all tuning parameters in execute() instead of slots
in CytoPipeline object.

CytoPipeline 0.99.3

- some minor changes for BiocCheck()

CytoPipeline 0.99.2

- removed dependencies to a number of packages, moved corresponding
implementations of CytoProcessingSteps (wrappers) into
CytoPipelineUtils package

CytoPipeline 0.99.1

- Maintenance due to Bioc version change (3.17)
- removed use of openCyto::gate_tail() (disappeared w/o deprecation),
replaced by flowDensity::deGate()
- implemented export of pre-processed file (writeFlowFrame as a
CytoProcessingStep implementation)
- extended readSampleFiles : mapping between channels and markers
- selectRandomSamples (new CytoProcessing step implementation)

CytoPipeline 0.99.0

- Prior to Bioconductor submission

[CytoPipelineGUI](/packages/CytoPipelineGUI)
---------------

                        Changes in version 0.99                         

CytoPipelineGUI 0.99.2

- added ShinyApps as additional BiocView

CytoPipelineGUI 0.99.1

- Modified man pages and vignettes to address comments raised during
Bioconductor submission

CytoPipelineGUI 0.99.0

- First version

[cytoviewer](/packages/cytoviewer)
----------

                 Changes in version 1.1.3 (2023-10-19)                  

- enlarged default display area and changed default color settings

                 Changes in version 1.1.2 (2023-10-09)                  

- vignette updates

                 Changes in version 1.0.1 (2023-07-27)                  

- added validation messages and fixes for metadata overlay

- added more tests and validation checks

[DCATS](/packages/DCATS)
-----

                 Changes in version 0.99.7 (2023-05-25)                 

- Added CITATION.cff file

                 Changes in version 0.99.6 (2023-05-25)                 

- Updated license information and changed the version of R
requirement
to R (>= 4.1.0)

                 Changes in version 0.99.5 (2023-05-12)                 

- Updated version to trigger a new biuld in Bioconductor

                 Changes in version 0.99.4 (2023-05-11)                 

- Made changes recommended in Bioconductor review process

                 Changes in version 0.99.3 (2023-04-21)                 

- Added NEWS file
- Made changes recommended in Bioconductor review process

- Made changes recommended in Bioconductor review process

[dearseq](/packages/dearseq)
-------

                 Changes in version 1.13.3 (2023-06-16)                 

- bug fix on Linux: default number of cores for parallel computations
is now detect_cores(logical=FALSE) - 1

                 Changes in version 1.12.1 (2023-05-28)                 

- weights_var2test_condi is enforced to FALSE for the permutation
test

[decontX](/packages/decontX)
-------

                 Changes in version 0.99.5 (2023-10-19)                 

- First submission to Bioconductor after review.

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

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.28.0                        

NEW FEATURES

- Add coercion from DelayedArray to SparseArray.

- Add efficient rowVars/colVars methods for DelayedMatrix objects.
  These methods, like all other row/col summarization methods
  implemented
  in the DelayedArray package, use block processing and can handle
  blocks
  of **arbitrary** geometry, that is, they can handle a grid of class
  ArbitraryArrayGrid (the most general type of grid).

- Add 'useNames' arg to row/colMins, row/colMaxs, row/colRanges, and
  row/colVars methods for DelayedMatrix objects.

- Add 'current_viewport' argument to set_grid_context().

SIGNIFICANT USER-VISIBLE CHANGES

- DelayedArray now depends on S4Arrays and SparseArray.

- Some improvements to the rowMeans/colMeans methods for DelayedMatrix
  objects.

[DelayedMatrixStats](/packages/DelayedMatrixStats)
------------------

                        Changes in version 1.23                         

- colAnyMissings() and rowAnyMissings() are deprecated in favour
  of the long-recommended colAnyNAs() and rowAnyNAs(). These
  deprecated functions will be made defunct in the next release
  cycle.

[demuxSNP](/packages/demuxSNP)
--------

                        Changes in version 0.99                         

Submitted to Bioconductor

Features

Functions to create a supervised demultiplexing classification model
using cell hashing and SNPs to aid in:

- SNP selection
- creating training data using demuxmix
- adding SNPs data to SingleCellExperiment objects
- cell assignment

[DepecheR](/packages/DepecheR)
--------

                 Changes in version 1.17.4 (2023-08-14)                 

- Adding collapse package to the description file.

- Correction of the generateSparseData function.

- Set.seed removed from two places and now non-present in the package.

                 Changes in version 1.17.3 (2023-07-26)                 

- Correction to groupProbPlot, as the last version threw out the most
  significant
  cells.

                 Changes in version 1.17.2 (2023-07-14)                 

- Correction to scale in groupProbPlot.

                 Changes in version 1.17.1 (2023-07-13)                 

- A mode option is added to the neighSmooth function.

- The scale of the groupProbPlot output is corrected, to -100, over 50
  to 100.

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.42.0                        

- collapseReplicates() now noisier (warning) about other assays.

- plotPCA() prints out the `ngenes` setting.

- Added pcsToUse for plotPCA (idea from Vince Carey)

- Added test that SE exist for lfcShrink (in case of
  glmGamPoi fitType).

                       Changes in version 1.41.13                       

- collapseReplicates() now noisier (warning) about other assays.

- plotPCA() prints out the `ngenes` setting.

                       Changes in version 1.41.8                        

- Use of MatrixGenerics

                       Changes in version 1.41.4                        

- Fixed dispersionFunction argument issue.

                       Changes in version 1.41.2                        

- Added pcsToUse for plotPCA (idea from Vince Carey)

                       Changes in version 1.41.1                        

- Added test that SE exist for lfcShrink (in case of
  glmGamPoi fitType).

[DESpace](/packages/DESpace)
-------

                        Changes in version 1.0.2                        

- sample-specific covariates allowed

                        Changes in version 1.0.1                        

- bug fixed in vignettes (giving error: "there is no package called
  'reticulate'")

[DifferentialRegulation](/packages/DifferentialRegulation)
----------------------

                       Changes in version 1.98.0                        

- bulk RNA-seq data allowed: the package now works with both bulk and
  single-cell RNA-seq data.

[dittoSeq](/packages/dittoSeq)
--------

                        Changes in version 1.14                         

- Feature Extensions:
1.  'dittoDotPlot()' & 'dittoPlotVarsAcrossGroups()': Improved
'group.by' ordering control via retention of factor levels and
addition of a new 'groupings.drop.unused' input to control
retention of empty levels.
2.  'dittoHeatmap()': Targeting Seurat clusters with the "ident"
shortcut now works for the 'annot.by' input of 'dittoHeatmap()'.
- Bug Fixes:
1.  'dittoHeatmap()': Fixed a feature exclusion check in
'dittoHeatmap()' meant to remove features without any non-zero
values. Previously, it removed all features with a mean of zero,
which blocked plotting from pre-scaled data.
2.  'dittoDimPlot()' & 'getReductions()': Eliminated cases where
'getReductions()' did not return NULL for 'object's containing
zero dimensionality reductions. This fix also improves
associated error messaging of 'dittoDimPlot()' where such cases
were missed.

[DMRcate](/packages/DMRcate)
-------

                       Changes in version 2.15.1                        

- DMR.plot() updated: ellipsis removed; biomaRt gene tracks now used;
  collapseTranscripts="meta"; exonAnnotation="symbol"; overlapping
  regions plotted as optional extra

- goregion() ontology changed to KEGG and for hypomethylated DMRs only

[DNAcopy](/packages/DNAcopy)
-------

                       Changes in version 1.75.5                        

- Fixed the bug in passing weights for weighted segmentation

                       Changes in version 1.75.4                        

- Updated the reference and source for Coriell data

                       Changes in version 1.75.3                        

- Added init.c to register native (Fortran) routines and to disable
symbol search
- Gzipped cytoBand.tab and converted default.DNAcopy.bdry to rda file
in data directory

                       Changes in version 1.75.2                        

- Added a NEWS.md file to track changes to the package.
- changed all dfloat in fortran to dble (Ripley email for
CRAN/clinfun)

[DOSE](/packages/DOSE)
----

                       Changes in version 3.27.3                        

- update TERM2NAME() to return term if corresponding name not found.
(2023-10-09, Mon)

                       Changes in version 3.27.2                        

- use 'MPO.db' and 'HPO.db' to support phenotype ontology for mouse
and human (2023-06-30, Fri)

                       Changes in version 3.27.1                        

- options(enrichment_force_universe = TRUE) will force enrichment
analysis to intersect the universe with gene sets (2023-05-03, Wed)
- use inherits to judge the class of objects (2022-11-20, Sun)
- test whether slot in GSON object is NULL (e.g., GSON@keytype) when
assigning it to enrichment result (2022-11-07, Mon)

[dreamlet](/packages/dreamlet)
--------

                       Changes in version 0.99.28                       

- Sept 5, 2026
- Update error handling for processAssays() and fitVarPart()

                       Changes in version 0.99.26                       

- August 18, 2026
- Update error handling and documentation

                       Changes in version 0.99.23                       

- August 8, 2023
- run styler::style_pkg()

                       Changes in version 0.99.22                       

- August 8, 2023
- dreamletCompareClusters() now allows cell-level covariates in
response to https://github.com/GabrielHoffman/dreamlet/issues/11
- Fix code for Bioconductor submission

                       Changes in version 0.99.21                       

- July 17, 2023
- Improve functionality and documentation of dreamlet::residuals()

                       Changes in version 0.99.20                       

- June 29, 2023
- in processAssays() use voomWithDreamWeights(..., span="auto") to
estimate the lowess tuning parameter
- rare error in merge_metadata() when a cell type is not observed
for all donors.

                       Changes in version 0.99.19                       

- June 28, 2023
- in dreamlet() fix issue when contrasts are specified and formula
includes variable from metadata()

                       Changes in version 0.99.18                       

- June 27, 2023
- add assays argument to buildClusterTreeFromPB()

                       Changes in version 0.99.14                       

- June 16, 2023
- bug fix in processAssays() when assays is dropped

                       Changes in version 0.99.13                       

- May 31, 2023
- improved error reporting
- Compatibility with variancePartition v2.0.5 (renamed 1.31.1)

                       Changes in version 0.99.12                       

- May 24, 2023
- required zellkonverter (>= 1.10.1) to avoid issues with previous
version
- issue solved by
https://github.com/theislab/zellkonverter/blob/b56718d113327020c024e188d9ac67ea57eaf35d/R/AnnData2SCE.R#L351

                       Changes in version 0.99.11                       

- May 12, 2023
- Compatibility with variancePartition v2.0.1

                       Changes in version 0.99.10                       

- April 24, 2023
- Fix issue in raised in Bioconductor submission:
https://github.com/Bioconductor/Contributions/issues/2955#issuecomment-1498070980
- Compatibility with variancePartition v2.0.0

                       Changes in version 0.99.6                        

- March 29, 2023
- fix topTable() for dreamletResult in the case where one or more
cells didn't estimate the coefficient of interest

                       Changes in version 0.99.3                        

- March 23, 2023
- reduce compiler time
- add computeNormCounts() and computeLogCPM()

                       Changes in version 0.99.1                        

- March 20, 2023
- fix Biocondcutor submission based on
https://github.com/Bioconductor/Contributions/issues/2955#issuecomment-1476037237

                       Changes in version 0.99.0                        

- March 15, 2023
- Biocondcutor submission


[dupRadar](/packages/dupRadar)
--------

                       Changes in version 1.30.3                        

- Fix warning associated with documentation

- Ignore installation of license file

                       Changes in version 1.30.2                        

- Correctly added roxygen2 imports to remove warnings on package checks

[easylift](/packages/easylift)
--------

                        Changes in version 1.0.0                        

- Bioconductor release 3.18

                       Changes in version 0.99.9                        

- Bug fixes and further improvements to the code.

                       Changes in version 0.99.1                        

- Updated the code and documentation to reflect the changes suggested
by the Bioconductor team.

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[EBSeq](/packages/EBSeq)
-----

                        Changes in version 2.0.0                        

- A new implementation of the core functions EBTest and EBMultiTest,
  which scales up the computation for big number of conditions to
  compare.

[edgeR](/packages/edgeR)
-----

                        Changes in version 4.0.0                        

- 
  New statistical methods implemented in glmQLFit() to ensure
  accurate estimation of the quasi-dispersion for data with small
  counts. The new method computes adjusted residual deviances
  with adjusted degrees of freedom to improve the chisquare
  approximation to the residual deviance. The new methodology
  includes the new argument 'top.proportion' for glmQLFit() to
  specify the proportion of highly expressed genes used to
  estimate the common NB dispersion used in the new method. The
  output DGEGLM object contains new components `leverage`,
  `unit.deviance.adj`, `unit.df.adj`, `deviance.adj`,
  `df.residual.adj` and `working.dispersion`. The new method can
  be turned on `legacy=FALSE`. By default, glmQLFit() will give
  the same results as in previous releases of edgeR.

- 
  New argument 'covariate.trend' for glmQLFit() to allow a
  user-specified covariate for the trended prior used to estimate
  the quasi-dispersions.

- 
  The gene set testing functions roast(), mroast(), fry(),
  camera() and romer() now have S3 methods for DGEGLM objects.

- 
  The edgeR Introductory vignette is converted from Sweave and
  pdf to Rmd and html.

- 
  Revised help pages for filterByExp() and catchSalmon().

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.21.3                        

- set_enrichplot_color(), a helper function to set colors
(2023-09-13,
Wed)
- change default color: from c("red", "blue") to c("#e06663",
"#327eba")
- use check_installed() to check package dependency (2023-09-08, Fri,
#254)

                       Changes in version 1.21.2                        

- introduce 'facet' parameter in dotplot() method for
compareClusterResult. If facet = "intersect", the dots will be
separated by enriched pathway intersection among clusters. It can
set to other variable that can be used for splitting the figure
(e.g., "category" for KEGG results) (2023-08-21, Mon)

                       Changes in version 1.21.1                        

- fixed cnetplot.compareClusterResult() for only contains one cluster
(2023-05-24, Wed, #243)

[enrichViewNet](/packages/enrichViewNet)
-------------

                       Changes in version 0.99.2                        

NEW FEATURES

- The man pages are respecting 80 character width.

                       Changes in version 0.99.1                        

NEW FEATURES

- The 'Installation' and 'Introduction' sections of the vignette have
been updated.

                       Changes in version 0.99.0                        

NEW FEATURES

- The new 'createEnrichMap()' function enables the creation of an
enrichment map from enrichment results.


[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.25.1                        

- Skip reading gtf file from Ensembl ftp server in unit test.

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.9.8 (2023-09-29)                  

- creates sample BAMs

- linearized MHL

                 Changes in version 1.9.4 (2023-07-03)                  

- methylation calls for bwa-meth, etc

                 Changes in version 1.9.2 (2023-06-21)                  

- both paired-end and single-end alignments

[epigraHMM](/packages/epigraHMM)
---------

                        Changes in version 1.9.1                        

- Fix normalization of log-probabilities in cpp code to avoid
underflow

[epistack](/packages/epistack)
--------

                 Changes in version 1.7.1 (2022-07-21)                  

- New parameters in plotEpistack(): rel_widths and rel_heights, to
adjust the relative widths and heights of the panels.
- In plotEpistack(), the boxMetric panel is now a bit higher.
- plotMetric() now have a ylab parameter, exposed in plotEpistack()

- left-most panel in plotEpistack() is a bit wider
- fixed a code typo in vignette

[escheR](/packages/escheR)
------

                        Changes in version 1.2.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add generic functions to support SingleCellExperiment object by
providing values to the argument dimred, and data.frame object
- Modify the generic function make_escheR.SpatialExperiment to
support
SpatialExperiment that is beyond Visium
- Update and add new vignette to show how to work with
SingleCellExperiment to visualize dimred and color palette for
bi-variate display.

                        Changes in version 1.1.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- The dependency on the package spatialLIBD are removed
- Revise the README file to add explanations for Gestalt Principles
to
help users grasp the big picture idea without going to read the
manuscript.

[esetVis](/packages/esetVis)
-------

                       Changes in version 1.27.1                        

- rbokeh plot - fixes error:
  - when no aesthetic variable is specified
  - dataPlotWithAnnotationWthtNA not found

- add GO.db to the 'Suggests'

[EWCE](/packages/EWCE)
----

                        Changes in version 1.9.3                        

Bug fixes

- generate_bootstrap_plots
- Use stored gene_data object whenever possible.
- Only show filtered celltypes.

                        Changes in version 1.9.2                        

Bug fixes

- generate_bootstrap_plots
- Missing available parameters for check_ewce_genelist_inputs()
call.

                        Changes in version 1.9.1                        

Bug fixes

- drop_uninformative_genes
- Hash out DGE options. Somehow these got re-exposed to users in
Bioc>=3.16.

[ExperimentHub](/packages/ExperimentHub)
-------------

                        Changes in version 2.9.0                        

DOCUMENTATION

- (2.9.1) Document access behind a proxy

[extraChIPs](/packages/extraChIPs)
----------

                       Changes in version 1.5.14                        

- Added p_mu0 to output of fitAssayDiff()
- Added respectLevels and filtering to plotProfileHeatmap()

                       Changes in version 1.5.13                        

- Added handling of unquoted column names to most plotting functions
- Added passing of specific columns to dualFilter
- Added drop to addDiffStatus

                       Changes in version 1.5.12                        

- Added plotGrlCol()

                       Changes in version 1.5.11                        

- Added defineSeqinfo()

                       Changes in version 1.5.10                        

- Added bed format to importPeaks()

                        Changes in version 1.5.8                        

- Added control of side-axis label position for plotProfileHeatmap()
- Added option to return merged key-value ranges for mergeByHMP()

                        Changes in version 1.5.7                        

- Matched DiffBind and csaw settings for fitAssayDiff()
- Added min_win to all merging functions
- Added n_max to getProfileData()

                        Changes in version 1.5.6                        

- Added mapGrlCols()

                        Changes in version 1.5.5                        

- Added plotPairwise() and addDiffStatus()

                        Changes in version 1.4.2                        

- Added Fixed-width vignette & edited sliding window
- Added se and peaks as example data for man pages and vignettes
- Included defineRegions()
- Enabled plotHFCG() without Ideogram tracks
- Enabled use of offsets for normalisation in fitAssayDiff()

[fastreeR](/packages/fastreeR)
--------

                 Changes in version 1.5.2 (2023-08-24)                  

- Update java backend to BioInfoJavaUtils-1.4.0 (haploid GT in
  vcf).

                 Changes in version 1.5.1 (2023-04-30)                  

- Update java backend to BioInfoJavaUtils-1.3.1 (fasta header
  name).

                 Changes in version 1.5.0 (2023-04-27)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_17
  branch.

[fenr](/packages/fenr)
----

                       Changes in version 0.99.7                        

- Minor changes to prepare for Bioconductor release
- Reverting temporarily to readr version 1 to circumvent a vroom
1.6.4
bug

                       Changes in version 0.99.6                        

In response to reviewer's comments

- The wording in the vignette was adjusted to more clearly convey the
purpose of the package to users
- Rewritten the description in DESCRIPTION file to clearly convey the
purpose of the package to users

BioPlanet seems defunct

- Removed BioPlanet for good, as their webpage is continuously down
and the maintainer is not responding

Minor adjustments to speed up building and testing

- Removed KEGG from interactive example to speed up vignette building
(GO and Reactome are sufficient for a simple example)
- Replaced yeast with simpler organisms in Wiki and KEGG tests to
speed up testing
- Replaced yeast with simpler organisms in Wiki and KEGG examples to
speed up checking

                       Changes in version 0.99.5                        

- Major overhaul following comments from Bioconductor's reviewer.

                       Changes in version 0.99.4                        

- Taking BiocCheck new warnings into account: adding @return to data
roxygens.

                       Changes in version 0.99.3                        

- Continuing issues with access to BioPlanet. fetch_bp example is now
marked donotrun and testing fetch_bp is removed to ensure smooth
build and check even when BioPlanet server is down.

                 Changes in version 0.99.2 (2023-05-24)                 

- BioPlanet's tripod.nih.gov SSL certificate seems to be fixed, so
reversing to the original read_csv code.

                 Changes in version 0.99.1 (2023-04-25)                 

- BioPlanet database vanished from internet and there is no sign of
it
coming back. Removing all BioPlanet-related code and replacing
BioPlanet with GO in the vignette and examples (this, alas, makes it
longer to check).
- OK, it is back, but I keep GO examples and vignettes.
- Minor improvements to documentation.

                 Changes in version 0.99.0 (2023-04-20)                 

- Pre-release Bioconductor version.


[fishpond](/packages/fishpond)
--------

                        Changes in version 2.8.0                        

- Corrected a bug in the two-group interaction (without pairing)
functionality, when the groups were imbalanced, as identified by
Samuel Chen. Fixes GitHub issue #35.

                        Changes in version 2.7.1                        

- Corrected a bug in the two-group interaction (without pairing)
functionality, when the groups were imbalanced, as identified by
Samuel Chen. Fixes GitHub issue #35.

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

[gcatest](/packages/gcatest)
-------

                 Changes in version 2.1.9 (2023-07-26)                  

- Commented out one more strict test (NA deviances) that fail too
often on bioconductor.

                 Changes in version 2.1.8 (2023-07-18)                  

- Commented out two more strict tests (for non-negative deviances)
that fail too often on bioconductor.

                 Changes in version 2.1.7 (2023-06-20)                  

- Commented out various excessive tests against glm, which differ
more
often than expected due to poor or lack of convergence.
- Removed unused LaTeX package dependencies from vignette to prevent
errors restricted to specific testing platforms.
- Fixed ..density.. deprecation warning in vignette plot.

                 Changes in version 2.1.6 (2023-05-25)                  

- Version bump for bioconductor devel.

                 Changes in version 2.0.6 (2023-05-25)                  

- README.md upgraded links from http to https
- Minor doc reformatting automatically performed by roxygen2.

                 Changes in version 2.0.5 (2021-06-18)                  

- Lots of minor changes for Bioconductor update.
- DESCRIPTION:
- Updated to Authors@R.
- Lengthened "Description" paragraph.
- Increased R dependency from 3.2 to 4.0.
- Reformatted this NEWS.md slightly to improve its automatic
parsing.
- Added examples for function delta_deviance_lf.
- Updated vignette to reflect that lfa::read.bed has been
deprecated in favor of genio::read_plink and BEDMatrix objects.
- Updated README.md, including corrections to examples.
- Updated citations:
- README.md: only had GCATest paper link, now has full
citation and also full LFA citation.
- Vignette: used to point to LFA arXiv preprint, now points to
published paper.
- inst/CITATION: didn't exist! Now includes both LFA and
GCATest papers.
- Added LICENSE.md.
- Internal changes:
- All unexported functions are now prefixed with a period.
- Replaced 1:x with seq_len(x) several functions.
- Reformatted all code with package reformatR and otherwise
match Bioconductor guidelines.

               Changes in version 2.0.4.9000 (2021-05-13)               

- Function delta_deviance_lf debugged case where either LF0 or LF1 is
a column matrix. Previously these 1-column matrices were getting
dropped to a vector incorrectly, which resulted in the mysterious
error message "Error: argument is of length zero". This 1-column
case is not typically observed in gcatest, but is common in the
reverse-dependent jackstraw package.

               Changes in version 2.0.3.9000 (2021-05-11)               

- Added function delta_deviance_lf, which calculates the delta
deviance from two logistic models and the genotype matrix data. This
function is a more general version of gcat.stat (which uses the new
function internally), to essentially consider models that differ by
more than one degree of freedom. It was written in particular for an
external application in mind, namely the jackstraw package.
- Internal function assoc_snp was renamed to delta_deviance_snp_lf
and
its last argument changed to match that of delta_deviance_lf
(alternative logistic factors instead of trait).

               Changes in version 2.0.2.9000 (2021-03-01)               

- Added internal tests for deviance calculations against stats::glm.
- Deviance code (internal delta_deviance_snp) now returns NA instead
of stopping when an "impossible" case is encountered (when the
genotype x is non-zero but the fitted probabilities under either
null or alternative model are zero, or the alternative allele dosage
(x-2) has the same problem). These cases are clearly model fitting
failures, and can arise for common ill-defined problems,
particularly under binary adjustment variables passed to gcat
together with rare variants; these individual cases are not handled
any better by stats::glm, so it seemed most sensible to return NA at
such loci and not stop.

               Changes in version 2.0.1.9000 (2021-02-16)               

- Documentation updates:
- Fixed links to functions, in many cases these were broken
because of incompatible mixed Rd and markdown syntax (now
markdown is used more fully).

               Changes in version 2.0.0.9000 (2020-11-13)               

Major overhaul from last version (1.3.2, last updated 2016-10-06).
Visible differences are support for BEDMatrix and fewer cases in
which
association p-values are NA. Internally there was major code
restructuring, and added unit tests for all functions.

- User-facing changes: Functions gcat/gcatest/gcat.stat

- added support for BEDMatrix objects for the genotype matrix X.
- This consumes lower memory when the number of loci m is very
large, so it enables analysis of larger datasets.
- Fixed some cases where the test statistic (the delta deviance)
and ultimately the p-values were NA or NaN and are no longer
missing.
- One common case is when fitted probabilities were zero or
one, which used to lead to NaN deviances when their correct
contribution was instead zero (because the limit of p*log(p)
as p goes to zero is zero, not 0 * (-Inf) = NaN).
- Other NA and NaN cases are avoided in the lfa function
af_snp (fixed in lfa 2.0.0.9000, 2020-09-18) used to
estimate the individual-specific allele frequencies used
here to compute the delta deviance. However, in rare cases
the logistic regression in af_snp fails to converge or there
are other problems, resulting in NA values propagated to
GCATest's test statistic and p-values.
- Otherwise, the new delta deviance code (function
delta_deviance_snp) is more numerically-stable than before.

- Internal changes

- Separated R functions into one source file each.
- Added more input checks to all functions.
- Added .gitignore files from another project.
- Added unit tests for all functions using testthat.
- Removed internal assoc C code
- Previously only used for genotype data without missingness
(so practically not on real datasets)
- Was entirely redundant with lfa::af_snp, which is now called
in all cases instead.
- Had bugs concerning handling of p == 0 or 1 cases that are
better handled in assoc_snp R code
- Minor scattered changes solely to pass latest R CMD check
requirements.

[gDNAx](/packages/gDNAx)
-----

                 Changes in version 0.99.6 (2023-05-26)                 

USER VISIBLE CHANGES

- Submission of the first version to the Bioconductor project.

[gDR](/packages/gDR)
---

                 Changes in version 0.99.9 (2023-10-18)                 

- adjust NEWS to Bioc format

                 Changes in version 0.99.8 (2023-08-17)                 

- update README

                 Changes in version 0.99.7 (2023-08-17)                 

- clean-up
- update documentation
- simplify Docker-based installation

                 Changes in version 0.99.6 (2023-05-22)                 

- format the vignette with BiocStyle

                 Changes in version 0.99.5 (2023-05-11)                 

- fix related with data.table

                 Changes in version 0.99.4 (2023-04-20)                 

- switch to OSI license

                 Changes in version 0.99.3 (2023-04-13)                 

- update documentation (Bioc-compatibility)

                 Changes in version 0.99.2 (2023-04-07)                 

- update maintainer

                 Changes in version 0.99.1 (2023-04-04)                 

- update requirements (gDRcore)

                 Changes in version 0.99.0 (2023-03-23)                 

- switch to Bioc compatible versioning


[gDRcore](/packages/gDRcore)
-------

                Changes in version 0.99.43 (2023-10-17)                 

- adjust NEWS to Bioc format

                Changes in version 0.99.42 (2023-10-05)                 

- bump version of gDRtestData
- fix bug with merging controls in triple combo with additional
perturbations

                Changes in version 0.99.41 (2023-09-25)                 

- add support for adding custom annotations inside input files
- improve the performance

                Changes in version 0.99.40 (2023-09-25)                 

- fix bug with subsetting wrong combo matrix value

                Changes in version 0.99.39 (2023-09-19)                 

- extend the logic for matching missing controls

                Changes in version 0.99.38 (2023-09-12)                 

- set Drug3 as an official tertiary drug in the experiment

                Changes in version 0.99.37 (2023-09-04)                 

- fill NA by average values when there is no match with plate

                Changes in version 0.99.36 (2023-09-01)                 

- fill NA during aggregation of ref and trt data with mean

                Changes in version 0.99.35 (2023-08-17)                 

- fix issue with missing subsetting Day0 data

                Changes in version 0.99.34 (2023-08-16)                 

- update logic for supporting manifest and template files sharing the
same column

                Changes in version 0.99.33 (2023-08-10)                 

- update annotation column names for cell line annotation as per
changes in the gDRutils

                Changes in version 0.99.32 (2023-07-25)                 

- extended logic for supporting cols with dash, e.g. additional
perturbations with "-"

                Changes in version 0.99.31 (2023-07-19)                 

- update the logic for handling warnings in the pipeline

                Changes in version 0.99.30 (2023-07-13)                 

- fix issue with wrong merging of data.tables without nested
confounders

                Changes in version 0.99.29 (2023-07-07)                 

- add information about source type for cases without metric data
- refactor the logic for splitting raw data from metadata (get rid of
iterative approach)

                Changes in version 0.99.28 (2023-07-05)                 

- update logic for parallel computing

                Changes in version 0.99.27 (2023-06-29)                 

- optimize unit tests

                Changes in version 0.99.26 (2023-06-29)                 

- remove backward compatibility for old data model

                Changes in version 0.99.25 (2023-06-27)                 

- fix bug with missing rownames in normalized assay

                Changes in version 0.99.24 (2023-06-19)                 

- update logic for merging data.table objects

                Changes in version 0.99.23 (2023-06-13)                 

- replace order with data.table::setorder

                Changes in version 0.99.22 (2023-06-09)                 

- switch from merge to [[ for data.table objects

                Changes in version 0.99.21 (2023-06-06)                 

- switch from aggregate to data.table

                Changes in version 0.99.20 (2023-06-07)                 

- switch from zoo::rollmean to data.table::frollmean

                Changes in version 0.99.19 (2023-06-06)                 

- replaced reshape2 functions by functions from data.table

                Changes in version 0.99.18 (2023-05-31)                 

- fix managing of mixed types of raw data

                Changes in version 0.99.17 (2023-05-29)                 

- fix bug with subsetting data for calculating isobologram

                Changes in version 0.99.16 (2023-05-22)                 

- format the vignette with BiocStyle

                Changes in version 0.99.15 (2023-05-16)                 

- fix related with data.table

                Changes in version 0.99.14 (2023-05-15)                 

- rename excess to x to unify colnames in assay data

                Changes in version 0.99.13 (2023-05-10)                 

- refactor normalization_types in combo-specific assays

                Changes in version 0.99.12 (2023-05-09)                 

- utilize gDRutils::apply_bumpy_function in fit_SE

                Changes in version 0.99.11 (2023-05-05)                 

- fix bug with swapping untreated/vehicle values

                Changes in version 0.99.10 (2023-05-04)                 

- fix bug with data.table

                 Changes in version 0.99.9 (2023-04-21)                 

- utilize gDRutils::apply_bumpy_function in average_SE

                 Changes in version 0.99.8 (2023-04-20)                 

- switch to OSI license

                 Changes in version 0.99.7 (2023-04-19)                 

- fix bug with replacing vehicle to untreated values

                 Changes in version 0.99.6 (2023-04-19)                 

- moved wrapper fuctions from gDRtestData

                 Changes in version 0.99.5 (2023-04-18)                 

- update dependencies
- add fix for bioc-devel (correct sorting in merge test)

                 Changes in version 0.99.4 (2023-04-17)                 

- fix namespacing issue in examples
- add R 4.2 as dependency
- fix examples for normalize_SE

                 Changes in version 0.99.3 (2023-04-12)                 

- add logic for retrieving raw data from assay data

                 Changes in version 0.99.2 (2023-04-07)                 

- update maintainer

                 Changes in version 0.99.1 (2023-04-04)                 

- bugfix for the logic in 'cleanup_metadata'

                 Changes in version 0.99.0 (2023-03-24)                 

- make the package Bioc-compatible


[gDRimport](/packages/gDRimport)
---------

                Changes in version 0.99.25 (2023-10-17)                 

- adjust NEWS to Bioc format

                Changes in version 0.99.24 (2023-10-02)                 

- add functions & unit tests for converting gDR MAE to PSet

                Changes in version 0.99.23 (2023-09-22)                 

- correct plate size calculation

                Changes in version 0.99.22 (2023-09-20)                 

- set barcode as character in the manifest file

                Changes in version 0.99.21 (2023-09-14)                 

- disable support for 'xls' file format due to crashes

                Changes in version 0.99.20 (2023-08-25)                 

- refactor subsetting of data.table using colname

                Changes in version 0.99.19 (2023-07-19)                 

- update warning messages

                Changes in version 0.99.18 (2023-07-02)                 

- add BiocStyle

                Changes in version 0.99.17 (2023-06-27)                 

- add exception entry for invalid average dose-response data

                Changes in version 0.99.16 (2023-06-23)                 

- increase compression level (Tavor_2020.qs; < 5MB limit)

                Changes in version 0.99.15 (2023-06-22)                 

- replaced rds with qs

                Changes in version 0.99.14 (2023-06-13)                 

- switch from merge to [[

                Changes in version 0.99.13 (2023-06-05)                 

- replaced reshape2 functions by functions from data.table

                Changes in version 0.99.12 (2023-05-24)                 

- format the vignette with BiocStyle

                Changes in version 0.99.11 (2023-05-16)                 

- data.frame => data.table switch (next round of changes)

                Changes in version 0.99.10 (2023-05-04)                 

- switch from tibble to data.table in excel files

                 Changes in version 0.99.9 (2023-04-25)                 

- refactor tibble, data.frame --> data.table

                 Changes in version 0.99.8 (2023-04-20)                 

- switch to OSI license

                 Changes in version 0.99.7 (2023-04-20)                 

- clean-up vignette

                 Changes in version 0.99.6 (2023-04-19)                 

- add object S4 gdr_test_data

                 Changes in version 0.99.5 (2023-04-19)                 

- mocked PSets tests

                 Changes in version 0.99.4 (2023-04-17)                 

- add R 4.2 as dependency
- bugfix for Pset-related tests and examples (reset identifiers)

                 Changes in version 0.99.3 (2023-04-13)                 

- add minor improvements (BiocCheck compatibility)

                 Changes in version 0.99.2 (2023-04-11)                 

- add support for PharmacoGx

                 Changes in version 0.99.1 (2023-04-07)                 

- update maintainer

                 Changes in version 0.99.0 (2023-03-24)                 

- make the package Bioc-compatible


[gDRstyle](/packages/gDRstyle)
--------

                Changes in version 0.99.22 (2023-10-17)                 

- adjust NEWS to Bioc format

                Changes in version 0.99.21 (2023-10-02)                 

- add options to skip tests/lintering in checkPackage

                Changes in version 0.99.20 (2023-08-08)                 

- add deploy trigger to workflow template

                Changes in version 0.99.19 (2023-06-15)                 

- fix pattern for finding *.R files
- lintr R files from 'inst/shiny' (if present)

                Changes in version 0.99.18 (2023-06-09)                 

- add reshape2 to lintr config

                Changes in version 0.99.17 (2023-05-10)                 

- add check for data.frame-related functions
- update package versioning rules

                Changes in version 0.99.16 (2023-05-04)                 

- ignore note for exported functions without examples
- handle properly BiocCheck notes with mulitple lines (notes to be
ignored)

                Changes in version 0.99.15 (2023-05-02)                 

- ignore note for 50 lines per function in biocCheck

                Changes in version 0.99.14 (2023-04-27)                 

- removed CRAN check from biocCheck

                Changes in version 0.99.13 (2023-04-21)                 

- add check for BiocCheck's notes

                Changes in version 0.99.12 (2023-04-20)                 

- switch to OSI license

                Changes in version 0.99.11 (2023-04-17)                 

- avoid dependencies upgrade
- add examples check

                Changes in version 0.99.10 (2023-04-17)                 

- update style guide (package doc)

                 Changes in version 0.99.9 (2023-04-17)                 

- add R 4.2 as a dependency

                 Changes in version 0.99.8 (2023-04-13)                 

- fix format in NEWS.md

                 Changes in version 0.99.7 (2023-04-07)                 

- update maintainer

                 Changes in version 0.99.6 (2023-04-07)                 

- update the license

                 Changes in version 0.99.5 (2023-04-06)                 

- update maintainer

                 Changes in version 0.99.4 (2023-04-05)                 

- remove unstable test

                 Changes in version 0.99.3 (2023-04-05)                 

- update examples

                 Changes in version 0.99.2 (2023-04-04)                 

- update examples
- switch to lintr::linters_with_defaults
- add 'test_mode' parameter in installAllDeps

                 Changes in version 0.99.1 (2023-04-04)                 

- change location of NEW.md file

                 Changes in version 0.99.0 (2023-03-24)                 

- downgrade version to make it Bioconductor compatible
- fix unit tests


[gDRutils](/packages/gDRutils)
--------

                Changes in version 0.99.34 (2023-10-18)                 

- adjust NEWS to Bioc format

                Changes in version 0.99.33 (2023-10-09)                 

- add support for flattening averaged assays

                Changes in version 0.99.32 (2023-09-22)                 

- fix bug in the case of conc=0 for evaluating efficacy

                Changes in version 0.99.31 (2023-09-19)                 

- add wide_structure param to convert_mae_assay_to_dt

                Changes in version 0.99.30 (2023-09-08)                 

- updated experimentalist description in schema

                Changes in version 0.99.29 (2023-09-05)                 

- add Replicate as a new identifier

                Changes in version 0.99.28 (2023-09-05)                 

- improve the logic of standardize_MAE to keep SE-specific metadata
and be able to revert standardization

                Changes in version 0.99.27 (2023-08-01)                 

- keep unchanged names in DataFrame

                Changes in version 0.99.26 (2023-08-01)                 

- tidy code

                Changes in version 0.99.25 (2023-06-27)                 

- add assert for missing rownames

                Changes in version 0.99.24 (2023-06-22)                 

- replaced RDS with qs

                Changes in version 0.99.23 (2023-06-20)                 

- fix check in R 4.3

                Changes in version 0.99.22 (2023-06-12)                 

- switch from merge to [[

                Changes in version 0.99.21 (2023-06-12)                 

- replace order with data.table::setorder
- add support for custom identifiers in merge_SE

                Changes in version 0.99.20 (2023-06-07)                 

- switch from aggregate to data.table

                Changes in version 0.99.19 (2023-06-06)                 

- replaced reshape2 functions by functions from data.table

                Changes in version 0.99.18 (2023-05-22)                 

- format the vignette with BiocStyle

                Changes in version 0.99.17 (2023-05-22)                 

- fix related with data.table
- remove .get_treated_conditions and .get_untreated_conditions

                Changes in version 0.99.16 (2023-05-18)                 

- add support for merging combination-data assays

                Changes in version 0.99.15 (2023-05-12)                 

- update after unifying normalization types

                Changes in version 0.99.14 (2023-05-12)                 

- fix lintr

                Changes in version 0.99.13 (2023-05-09)                 

- removed cotreatment entry from EXPERIMENT_GROUPS

                Changes in version 0.99.12 (2023-05-09)                 

- fix bug in convert_mae_assay_to_dt

                Changes in version 0.99.11 (2023-05-08)                 

- refactor code with single ampersand in if statements

                Changes in version 0.99.10 (2023-04-28)                 

- change order of untreated tags

                 Changes in version 0.99.9 (2023-04-24)                 

- changed data.frame to data.table

                 Changes in version 0.99.8 (2023-04-20)                 

- switch to OSI license

                 Changes in version 0.99.7 (2023-04-20)                 

- clean-up vignette

                 Changes in version 0.99.6 (2023-04-18)                 

- extend the logic of apply_bumpy_function

                 Changes in version 0.99.5 (2023-04-17)                 

- add R 4.2 as a dependency

                 Changes in version 0.99.4 (2023-04-14)                 

- fix examples

                 Changes in version 0.99.3 (2023-04-13)                 

- make linter happy

                 Changes in version 0.99.2 (2023-04-12)                 

- add licence

                 Changes in version 0.99.1 (2023-04-07)                 

- update maintainer

                 Changes in version 0.99.0 (2023-03-28)                 

- downgrade version to make it Bioconductor compatible


[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.36.1                        

UTILITIES

- `gdsfmt:::.reopen()` allows forking

[gemma.R](/packages/gemma.R)
-------

                        Changes in version 2.0.0                        

- Breaking change to get_dataset_differential_expression_analyses
function in order to return annotations for contrasts with multiple
characteristics.

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 2.6.0                        

Bug fixes

- describe_gtl() correctly extracts the number of up and down
regulated genes from the DE results
- Fortified the behavior of gs_scores() to handle cases where only
one
gene would be included in the signature to plot

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.38.0                        

NEW FEATURES

- Register the following NCBI assemblies:
  - bStrHab1.2.pri (Kakapo)
  - a few Salmo salar (Atlantic salmon) assemblies
  - a few African elephant (Loxodonta africana) assemblies

BUG FIXES

- Switch from HTTP to HTTPS for requests to *.ucsc.edu

- Remove library() calls from inst/registered/UCSC_genomes/*.R files

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.38.0                        

- No changes in this version.

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.54.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Fix flaw in heuristic for inferring exon ranks in
  makeTxDbFromGRanges():
  makeTxDbFromGRanges() can guess the exon ranks of a given transcript
  either
  (a) based on their position on the chromosome
  (b) or by looking at the suffixes of the exon ids (e.g. .1, .2,
  etc...)
  if any
  Previously (i.e. in GenomicFeatures < 1.53.1) it was trying (b)
  first,
  and would fall back on (a) if (b) failed.
  Starting with GenomicFeatures 1.53.1, it does the opposite: it tries
  (a)
  first, and uses (b) only for exons for which (a) failed.
  This change addresses the problem that the suffixes of the exon ids
  cannot
  be trusted to infer the ranks of exons located on the minus strand.
  See
  https://github.com/Bioconductor/GenomicFeatures/issues/59 for an
  example.

- Try to more clearly distinguish between CDS and CDS parts in
  documentation.

BUG FIXES

- Switch from HTTP to HTTPS for requests to *.ucsc.edu

[GenomicPlot](/packages/GenomicPlot)
-----------

                       Changes in version 0.99.15                       

NEW FEATURES

None

SIGNIFICANT CHANGES

- Provide link to all external and internal data in manual
- Add examples to all plot functions

BUG FIXES

None

                       Changes in version 0.99.14                       

NEW FEATURES

None

SIGNIFICANT CHANGES

- Make txdb available in /inst/data to avoid repetitive generation of
it in examples
- Store gf5_meta and gf5_genomic in /data so that they can be
loaded by calling data()
- Create /R/data.R for data documentations

BUG FIXES

None

                       Changes in version 0.99.13                       

NEW FEATURES

None

SIGNIFICANT CHANGES

- Add 'title' to draw_combo_plot function arguments.

BUG FIXES

- Fixed a bug in plot_peak_annotation, such that txdb$user_genome is
treated as a vector of strings rather than a single string

                       Changes in version 0.99.12                       

NEW FEATURES

- Add PCA plot in plot_bam_correlation.
- The function plot_5parts_metagene can generate profile plots for
both 5 parts (Promoter, 5'UTR, CDS, 3'UTR, TTS) and 3 parts
(Promoter, Gene, TTS) metagene.

SIGNIFICANT CHANGES

- The function plot_3parts_metagene is removed.
- Add setImportParams function to provide default import
parameters.
- Add saveRds option to import parameters to control saving of
imported data, the default is FALSE.
- Add data type and missing file checking for function arguments

BUG FIXES

None

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.54.0                        

NEW FEATURES

- Add 'ignore.strand' arg to gaps() method for GenomicRanges objects.

BUG FIXES

- Small tweaks to GRanges() constructor and to coercion from GPos
  to GRanges.

[GenomicScores](/packages/GenomicScores)
-------------

                       Changes in version 2.14.0                        

USER VISIBLE CHANGES

- Added four new pathogenicity score sets and corresponding metadata
  for hg19 and hg38: AlphaMissense.v2023.hg19,
  AlphaMissense.v2023.hg38, cadd.v1.6.hg19 and cadd.v1.6.hg38.

- Updated vignette.

[gg4way](/packages/gg4way)
------

                       Changes in version 0.99.2                        

- Documentation updates

                       Changes in version 0.99.1                        

- Revised based on:
https://github.com/Bioconductor/Contributions/issues/3129

                       Changes in version 0.99.0                        

- Pre-release version of the package.

[ggkegg](/packages/ggkegg)
------

                 Changes in version 0.99.3 (2023-08-25)                 

- Added new files

                 Changes in version 0.99.2 (2023-08-25)                 

- Remove the unnecessary file to pass R CMD CHECK

                 Changes in version 0.99.1 (2023-08-25)                 

- Revising the codes based on the Bioconductor review

                 Changes in version 0.99.0 (2023-06-27)                 

- Submitted to Bioconductor

[ggsc](/packages/ggsc)
----

                       Changes in version 0.99.11                       

- support density visualization for single and spatial transcriptomic
data (2023-10-18, Wed)

                       Changes in version 0.99.10                       

- on Bioconductor (2023-10-16, Mon)
- add \value session in the reexports.Rd (2023-10-15, Sun)
- add package level man page and update vignette (2023-10-14, Sat)
- add examples in Rd to satisfy BiocCheck (2023-09-18, Mon, #7)
- sc_dim_count() function to generate a barplot from a dimension
reduction plot (sc_dim() plot) to visualize the number of cells for
each clusters (2023-09-13, Wed)
- add 'biocViews' in DESCRIPTION required by Bioconductor

                       Changes in version 0.99.0                        

- compatible with 'SingleCellExperiment' (2023-09-05, Tue, #5)
- using S4 OOP to reorganize the functions (2023-09-05, Tue, #4)
- rename the package to 'ggsc' as there is a package called 'scplot'
in CRAN
- add H&E image to sc_spatial() (#3)


[ggtree](/packages/ggtree)
------

                        Changes in version 3.9.1                        

- use rlang::check_installed() to check whether the suggested pkg is
installed (2023-08-11, Fri, #580)
- allows using linewidth parameter (synonyms for size) in ggtree()
(2023-07-15, Sat, #574)
- bug fixed in setting branch.length = "none" to plot 'hclust' object
(2023-07-15, Sat, #574)

[ggtreeDendro](/packages/ggtreeDendro)
------------

                        Changes in version 1.3.2                        

- update according to the change of tidytree (2023-08-18, Fri)

                        Changes in version 1.3.1                        

- autoplot method for 'dendro' object (ggdendro::dendro_data()
output)
(2023-03-02, Thu)
- see also https://github.com/YuLab-SMU/treeio/pull/95

[glmGamPoi](/packages/glmGamPoi)
---------

                  Changes in version 1.13 (2023-07-03)                  

- Implement a likelihood ratio test based on the Chi-squared
  distribution, if
  `test_de` is called after setting `overdispersion_shrinkage = FALSE`.
  Note that
  this test is less reliable than than the quasi-likelihood F test that
  is run
  for `overdispersion_shrinkage = TRUE`.

[GloScope](/packages/GloScope)
--------

                 Changes in version 0.99.5 (2023-10-23)                 

- Fix typos in README and CITATION files

                 Changes in version 0.99.4 (2023-10-03)                 

- Improve example data and test cases

- Update documentation

- Fix Monte Carlo sampling edge case bug

                 Changes in version 0.99.3 (2023-09-17)                 

- Address feedback from Bioconductor review

- Catch edge case where more GMM components specified by user than
  cells

                 Changes in version 0.99.2 (2023-08-09)                 

- Fix issues in Bioconductor build report

                 Changes in version 0.99.1 (2023-07-18)                 

- Update ROxygen examples and test cases to use SingleCellExperiment
  data

                 Changes in version 0.99.0 (2023-06-01)                 

- This is the pre-submission version.

[GNOSIS](/packages/GNOSIS)
------

                 Changes in version 1.99.0 (2023-09-04)                 

- Made the following significant changes
  o added functionality to select and upload cBioPortal study
  o deprecated ability to save R script with executed code

- Submitted to Bioconductor

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.27.3                        

- use check_installed() to check package dependency (2023-09-12, Tue,
#43)

                       Changes in version 2.27.2                        

- read.blast2go() to parse 'blast2go' result (2023-07-10, Mon)
- move buildGOmap() and read.gaf() from 'clusterProfiler'
(2023-07-10,
Mon)

                       Changes in version 2.27.1                        

- semantic similarity measurement support for MPO (2023-04-06, Thu)
- TCSS semantic similarity measurement support for DO and MPO
(2023-04-06, Thu)

[goSorensen](/packages/goSorensen)
----------

                        Changes in version 1.3.0                        

We add the following functions

- allBuildEnrichTable: Given k lists of genes, it generates the
k(k–1)/2 contingency tables of joint enrichment for all possible
pairs of lists.
- allEquivTestSorensen: Accepts the objects created with
allBuildEnrichTable as its first argument and quickly obtains the
k(k–1)/2 equivalence tests.
- sorenThreshold: Implements an algorithm allowing the computation of
the “equivalence threshold” dissimilarities matrix for all the
k(k–1)/2 tests. According to what is indicated in the arguments,
the results of this function are stored in a matrix for a specific
ontology and level or a list with a matrix for more than one
ontology and/or level.
- hclustThreshold: Generates an object of class "hclust". For a
specific ontology and level, plots a dendrogram where all the
k(k–1)/2 comparisons are joined at the height of their respective
"equivalence threshold" dissimilarity.
- allHclustThreshold: Performs the same calculations as
hclustThreshold but for the specified GO ontologies and levels (all
three ontologies PB, MF, and CC and levels from 2 to 10 by default)

In addition:

- We improve the vignette "An introduction to the goSorensen package"
by implementing the new functions mentioned above and updating the
results of the examples with the results of the latest version of
Bioconductor.

[GRaNIE](/packages/GRaNIE)
------

              Changes in version 1.5.2-1.5.3 (2023-08-20)               

New features

- say hello to a new function filterConnectionsForPlotting() that can
be used to include or exclude particular connections from the stored
eGRN for visualization purposes only (!). Note that this filter only
applies to visualization and enables a flexible system to visually
explore particular features of the stored eGRN. THis is particularly
handy when the eGRN is large. For more details, see the help pages
of the new function.
- similarly, the function visualizeGRN() now by default only
visualizes connections that are marked as such (the result from
filterConnectionsForPlotting()) - that is, it excludes connections
that the user beforehand excluded from plotting. This allows to
specifically plot only part of the eGRN network and explore specific
T&F regulons, for example, a feature that before was not so easy to
do.
- It is now possible to integrate SNP data into GRaNIE via the new
function addSNPData(). For more information, see the Package
vignette.

                 Changes in version 1.5.1 (2023-06-19)                  

- version jump due to new Bioconductor development cycle

New features and stability improvements

- we replaced biomaRt for the full genome annotation retrieval in
addData with a different approach that is more reliable, as we had
more and more issues with biomaRt in the recent past. While using
the old biomaRt approach is still an option, the default is now to
use the AnnotationHub package from Bioconductor. This makes GRaNIE
overall more stable and less reliant on biomaRt due to the strict
timeouts and query size restrictions.

[graphite](/packages/graphite)
--------

                 Changes in version 1.47.1 (2023-10-11)                 

- Removed Clipper support.

- Updated all pathway data.

[GreyListChIP](/packages/GreyListChIP)
------------

                       Changes in version 1.32.1                        

- Change maintainer email from Rory to Matt

[GSVA](/packages/GSVA)
----

                        Changes in version 1.50                         

USER VISIBLE CHANGES

- The API has changed. The main function remains under the same name
  'gsva()', but the way in which is called is different. From this
  release, it has three parameters only: the first is a parameter
  object whose class depends on the method to be used, the second is a
  flag to set verbosity and the third controls the parallelization of
  the calculations. The old way of using 'gsva()' has been deprecated,
  which means that during this release, the user may still use the old
  API, but will get a deprecation warning message. In the next release,
  the old way of using 'gsva()' will become defunct and prompt an
  error. Please check the help page of 'gsva()' for details.

BUG FIXES

- Bugfix for https://github.com/rcastelo/GSVA/issues/88 to correctly
  deal with a GeneSetCollection object as input gene sets, when the
  input expression data is a SingleCellExperiment object.

- Bugfix for https://github.com/rcastelo/GSVA/issues/90 to enable
  working with long vectors in the calls to C code by the GSVA
  algorithm.

[HarmonizR](/packages/HarmonizR)
---------

                       Changes in version 0.99.1                        

NEW FEATURES

- HarmonizR now accepts S4 summarized experiment input and will return
  the batch effect adjusted result as such

SIGNIFICANT OTHER CHANGES

- HarmonizR now imports SummarizedExperiment

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.30.0                        

NEW FEATURES

- Add 'dim' and 'sparse.layout' args to H5SparseMatrixSeed().

SIGNIFICANT USER-VISIBLE CHANGES

- HDF5Array now imports S4Arrays.

[HERON](/packages/HERON)
-----

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.38.0                        

- fix a compiler warning of "unused-but-set-variable" on Apple ARM
  chips

                       Changes in version 1.36.3                        

- the output of `hlaPredict(, type="response+prob")` includes dosages

- new arguments 'ret.dosage', 'ret.postprob', 'max.resolution' and
  'rm.suffix' in `hlaPredMerge()`

- new arguments 'allele.list' and 'prob.cutoff' in `hlaAlleleToVCF()`
  for more possible outputs

- `hlaAlleleToVCF()` accepts a list of 'hlaAlleleClass' as the first
  argument: output multiple 'hlaAlleleClass' objects to a single VCF
  file

[HiCcompare](/packages/HiCcompare)
----------

                 Changes in version 1.23.1 (2023-06-03)                 

- Fix cooler2sparse function when a cool file has only a single
  chromosome.
  Thanks to @junjunlab, @agata-sm,
  https://github.com/dozmorovlab/HiCcompare/issues/28

- Update cool file handling in the vignette.

[HiCDOC](/packages/HiCDOC)
------

                 Changes in version 1.3.2 (2023-09-22)                  

- Move license to a LICENSE file, as asked by open access repo

                 Changes in version 1.3.1 (2023-09-06)                  

- Add a cyclicLoessSpan parameter to speed technical bias normalization
  to
  speed the technical bias normalization if set (in
  `normalizeTechnicalBiases()`)

- Remove 2 dependencies (ggpubr and ggExtra)

[hoodscanR](/packages/hoodscanR)
---------

                       Changes in version 0.99.5                        

Remove codes that are commented out.

                       Changes in version 0.99.4                        

user perspective/experience

- documentation update about output values
- output variable/names parameter added to many functions
- validity checks for spatialCoords are added to related functions
- for loops are gone.

                       Changes in version 0.99.3                        

code

- sapply change to vapply.
- validity checks on inputs are now added to all functions using
either is(spe, "SpatialExperiment") or is(m, "matrix") etc.

vignette

- Explanation of plots and algorithm inside different functions are
now added to the vignette briefly.
- BiocStyle are now used.
- BiocStyle to hyperlink packages are now used.
- Fixed typo

documentation

- data description of the test data has been added.
- cross-link key function have added to the package help page.
- the equation and theory of the calcMetrics are now added in the
help
page.

other

- Bioconductor installation instructions is now added to the readme
file.
- code coverage is now almost 85%.
- The LICENSE placeholders is now edited.

                       Changes in version 0.99.0                        

First submission to Bioconductor.

[HPAanalyze](/packages/HPAanalyze)
----------

                        Changes in version 1.19                         

- Changes in version 1.19.0
  + Starting devel for Bioconductor 3.18

- Changes in version 1.19.1
  + Bug fix: Correct url for rna_tissue_consensus

[hypeR](/packages/hypeR)
-----

                        Changes in version 2.0.1                        

- Dots plots are now explicitly sized with size_by=c("genesets",
"significance", "none")

                        Changes in version 2.0.0                        

- Version bump for bioconductor

[IgGeneUsage](/packages/IgGeneUsage)
-----------

                        Changes in version 1.15                         

- model for gene usage (GU) analysis

- model for multi-condition analysis, ANOVA like

[illuminaio](/packages/illuminaio)
----------

                 Changes in version 0.43.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor devel version,
which is now Bioconductor 3.18 for R (>= 4.4.0).

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.7.8 (2023-10-19)                  

- spatialCoords are not initialised with rownames anymore

                 Changes in version 1.7.7 (2023-10-02)                  

- Bug fix: in rare cases the testInteractions p-values were not
  correctly computed due to machine precision issues

                 Changes in version 1.7.6 (2023-09-26)                  

- Updated example data in "inst" to newest version

                 Changes in version 1.7.5 (2023-09-14)                  

- Added more internal validity checks for 'read_steinbock'

- Changed the way messages in vroom are silenced

- Stop allowing Object numbers to be different between intensities and
  regionprops

- Stop allowing missing files in regionprops and neighbors

                 Changes in version 1.7.4 (2023-08-09)                  

- Bug fix: correctly setting the "aspect.ratio" argument in plotSpatial
  to fix the physical units of the x- and y-axis

                 Changes in version 1.7.3 (2023-06-20)                  

- More info on reproducibility using testInteractions

                 Changes in version 1.7.1 (2023-06-07)                  

- More info on reproducibility using detectCommunity

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.33.6                        

- CHANGES
  * introducing method option for SON/ormalization
  * a bit clarification in E/M-calls interface

                       Changes in version 1.33.2                        

- CHANGES
  * SON/ormalization back to Version 1.31.9

                       Changes in version 1.33.1                        

- CHANGES
  * normalization variants for combine clustering

[iNETgrate](/packages/iNETgrate)
---------

                Changes in version 0.99.120 (2023-06-29)                

General

- After the iNETgrate package was added to Bioconductor, Habil
  removed it from the Bitbucket and Github repositories.

                Changes in version 0.99.102 (2023-06-12)                

Changes in existing functions

- In prepareSurvival and cleanAllData function, the default
  values changes as follows: riskLow="Low" and riskHigh="High",
  respectively.

                Changes in version 0.99.100 (2023-06-09)                

Changes in existing functions

- The Labels argument is removed from analyzeSurvival.  and the
  riskCol argument is added.

                Changes in version 0.99.96 (2023-06-08)                 

Changes in existing functions

- The abnormalityCol argument is removed from accelFailAnalysis.

- The otherLabel argument is removed from prepareSurvival and
  cleanAllData functions.

                Changes in version 0.99.94 (2023-06-08)                 

Changes in existing functions

- The computeEigengenes, computeEigenloci, and sample2atient
  functions were renamed to computEigengenes, computEigenloci,
  and sample2pat, respectively.

- The makeNetwork function now saves the result for each mu in a
  separate folder.

                Changes in version 0.99.74 (2023-05-24)                 

General

- Data were compressed using "xz" to save 20% space.

Changes in existing functions

- The downloaData function is not exported anymore because it is
  specific to the legacy TCGA data, which is no longer available.

                Changes in version 0.99.70 (2023-05-23)                 

Changes in existing functions

- Habil renamed downloadData to downloaData and updated the
  data.category argument of GDCquery.

- Habil removed the legacy option from the donwlowData function
  because it is no longer supported by TCGAbiolinks.

                Changes in version 0.99.64 (2023-05-11)                 

General

- Habil commented the installation command line in the vignette.

                Changes in version 0.99.50 (2023-04-14)                 

General

- Habil changed function names to camelCase i.e, from
  function.name to functionName.

                Changes in version 0.99.46 (2023-04-11)                 

General

- Habil changed getwd() to tempdr() in the docs.

- Habil changed most of occurrences of cat() to paste().

- Habil changed most of occurrences of print() to message() or
  Pigengene::message.if().

- Habil removed paste() inside all stop() to warning().

- Habil replaced dontrun with donttest.

                Changes in version 0.99.30 (2023-04-06)                 

General

- Tested on R 4.3.

- TCGAbiolinks (<= 2.24.3) because of issues with newer versions.

                 Changes in version 0.5.42 (2021-10-31)                 

General

- Being prepared for submission to Bioconductor.

- Created.

[IntEREst](/packages/IntEREst)
--------

                       Changes in version 1.26.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- IntSpan method is changed to count reads that span the whole
  intron precisely. Previous these were counted separayted for 5'
  end and 3' end of the introns.

- repeatsTableToFilter is ignored for IntSpana nd ExSkip as the
  skipped regions of the reads (denoted by N or sometimes D in
  the CIGAR of bam/sam file ) are taken into count for these
  benchmarks (not the mapped regions which are dentoed with Ms in
  the cigar).

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.36.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add link to revElements() in man page for reverse().

BUG FIXES

- Fix is.unsorted() methods for Compressed[Integer|Numeric]List
  objects (they were never working since their introduction years
  ago).

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.11.2 (2023-07-26)                 

FIXES

- Fixed issues in html report for outlier filtering - reported
incorrect numbers due to missing conversion in percentage
- Fixed warnings for bslib::nav deprecation
- Fixed minor issue in default_af_transform(), transformation failed
if NAs were present in the columns

                 Changes in version 1.11.1 (2023-05-09)                 

FIXES

- Fixed broken tests with new updates in underlying packages

[iSEE](/packages/iSEE)
----

                       Changes in version 2.13.5                        

- Add generic .isBrushable to support panels that are not DotPlot
extensions.

                       Changes in version 2.13.4                        

- Fix COLORMAP bug introduced in 2.13.3.

                       Changes in version 2.13.3                        

- Let app maintainer define colormap in the landing page.

                       Changes in version 2.13.2                        

- Fix bug introduced in 2.11.2 (DataBoxOpen would apply also to
Visual
parameters box)

                       Changes in version 2.13.1                        

- Define missing methods for generics for custom tables (issues #608
and #612)

[iSEEde](/packages/iSEEde)
------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[iSEEhub](/packages/iSEEhub)
-------

                        Changes in version 1.3.2                        

BUG FIXES

- Cleaner version of bug fix introduced in version 1.3.1.

                        Changes in version 1.3.1                        

BUG FIXES

- Fix cleaning of missing rowData and colData.

[iSEEindex](/packages/iSEEindex)
---------

                       Changes in version 0.99.13                       

- Added possibility to inject custom header and footer in landing
page.

                       Changes in version 0.99.12                       

- Added control over inclusion and location of default initial
configuration amongst available choices.

                       Changes in version 0.99.11                       

- Added second example data set to replace copy of first one.

                       Changes in version 0.99.10                       

- Added support for custom tours through initial state configuration
scripts.

                       Changes in version 0.99.9                        

- Disable GitHub CommonMark extensions and demonstrate
target="_blank"
for issue #41.

                       Changes in version 0.99.8                        

- Fixed duplication of UI output displaying information about initial
configuration.

                       Changes in version 0.99.7                        

- Fix package man page.

                       Changes in version 0.99.6                        

- Added man page for package.

                       Changes in version 0.99.5                        

- Changed options(width=120) in vignettes when displaying session
info.

                       Changes in version 0.99.4                        

- Removed export of internal function in NAMESPACE.

                       Changes in version 0.99.3                        

- Deduplicated choice of initial configuration.

                       Changes in version 0.99.2                        

- Added \value section to fix BiocCheck::BiocCheck() WARNING.

                       Changes in version 0.99.1                        

- Trigger a new build on the Bioconductor Single Package Builder.

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[iSEEpathways](/packages/iSEEpathways)
------------

                       Changes in version 0.99.2                        

MINOR UPDATES

- Added missing depedency to DESCRIPTION file (Bioconductor WARNING).
- Moved set.seed() out of functions (Bioconductor WARNING).
- Added \value section to man pages (Bioconductor WARNING).

                       Changes in version 0.99.1                        

MINOR UPDATES

- Added BiocViews (Bioconductor NOTE).

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[IsoformSwitchAnalyzeR](/packages/IsoformSwitchAnalyzeR)
---------------------

                       Changes in version 2.01.14                       

- Update type: Minor.

- Update for Bioconductor

                       Changes in version 2.01.13                       

- Update type: Minor.

- Fixed a problem with analysis of 5_utr_seq_similarity in
  analyzeSwitchConsequences()

- importRdata() was updated to handle sva analysis better

- importRdata() was updated by removing the addIFmatrix argument as the
  IF matrix is now alwasy needed

- importRdata() had it's detectAndCorrectUnwantedEffects argument
  updated to

- isoformSwitchTestDEXSeq was updated to not batch correct IF values as
  this is already done by importRdata

- Various documentation updates

                       Changes in version 2.01.12                       

- Update type: Minor.

- Update of switchPlot() to turn off topology plotting

- Update of importRdata() to better handle datasets with no replicates

                       Changes in version 2.01.11                       

- Update type: Minor.

- importRdata() was updated to fix problem with fasta import.

                       Changes in version 2.01.10                       

- Update type: Minor.

- Updated satuRn version requirement

- Updated importRdata() to allow skipping sva analysis incoperation.

- Updated importRdata() documentation accordingly.

- Updated importRdata() documentation to better describe the
  switchAnalyzeRlist created.

- Updated isoformSwitchTestSatuRn() to be more robust to various id
  types.

                       Changes in version 2.01.09                       

- Update type: Minor.

- Updated importRdata() to also handle when there are to few samples to
  run SVA.

                       Changes in version 2.01.08                       

- Update type: Minor.

- Updated importRdata() to use more stringent filtering (inspired by
  edgeR::filterByExpr()) before running SVA. Output in final
  switchAnalyzeRlist is not affected (aka that have not been filtered).

- Updated importRdata() to also handle when to many SVAs are found.

- Updated importRdata() to also handle when there are to few samples to
  run SVA.

                       Changes in version 2.01.07                       

- Update type: Minor.

- Fixed an edgecase bug in importRdata()

                       Changes in version 2.01.06                       

- Update type: Minor.

- Fixed an bug in isoformSwitchAnalysisPart2() that could result in
  problem when running without toplogy analysis.

- Introduced a better error message in analyzeORF().

                       Changes in version 2.01.05                       

- Update type: Minor.

- Updated switchPlotTranscript() to make a message instaed of an error
  when plotTopology=TRUE but isoform topology had not beed added.

- More detailed descriptions of analyzeDeepTMHMM() and
  analyzeDeepLoc2() added to the vignette.

                       Changes in version 2.01.04                       

- Update type: Minor.

- Fix to handle duplicated levels

                       Changes in version 2.01.03                       

- Update type: Minor.

- Fixes to accomodate dplyr updates

                       Changes in version 2.01.02                       

- Update type: Minor.

- Fixed a problem with batch correction in importRdata()

                        Changes in version 2.1.2                        

- Update type: Minor.

- Documentation update for Bioconductor

                       Changes in version 2.01.01                       

- Update type: Major.

- createSwitchAnalyzeRlist() was removed. All users should instead use
  importRdata().

- importRdata() now automatically detects un-annoated covariates in
  data via the sva package.

- importRdata() now automatically corrects abundance and isoform
  fractions for unwanted covariates (both used supplied and those found
  via sva).

- Accordingly all batch correction functionallity in the
  isoformSwitchTestDEXSeq() function was removed.

- isoformSwitchTestSatuRn() was introduced. This test uses satuRn for
  switch identification which works extremely well for larger sample
  sizes. Huge thanks to Jeroen Gilis making this functionality and the
  pull request!

- Accordingly the suboptimal isoformSwitchTestDRIMSeq function have
  been removed. All documentation was updated accordingly.

- IsoformSwitchAnalyzeR now depends on the R package pfamAnalyzeR for
  analyzing pfam domain isotypes.

- analyzeSignalP() was updated to support import of results predicted
  with SignalP6.

- analyzeDeepTMHMM() was introduced to add topological predictions to
  the switchAnalyzeRList.

- analyzeDeepLoc2() was introduced to add predictions of sub-cellular
  localization to the switchAnalyzeRList.

- analyzeIUPred2A() was tested against with result files from IUPred3
  and seem to work.

- analyzeSwitchConsequences() was updated to predict a number of new
  consequences based on the new annoation described above.

- analyzeSwitchConsequences()'s AaFracCutoff default was updated from
  0.5 to 0.8 resulting in more lenient differenceses being identified.

- extractSubCellShifts() was introduced to enable a deeper analysis of
  changes in sub-cellular localization due to isoform switches.

- Vignette was updated to recomend IsoQuant instead of TALON for long
  read data.

- analyzePFAM() was updated to import envelope (instead of alignment)
  coordinates as currently recomended. In practice this is a minor
  change for most domains.

- Example data was updated to reflect new annoation and consequences
  that can be predicted

- Various code corrections and improvements

- Various documentation improvements

                        Changes in version 2.1.1                        

- Update type: Minor.

- Update for Bioconductor

[KEGGREST](/packages/KEGGREST)
--------

                       Changes in version 1.42.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- `keggCompounds` lists compound IDs for a given pathway
  (@KristinaRiemer,
  #6).

BUG FIXES

- Update URL path in `.get.kegg.url` from `tmp` to `kegg` subfolder.

[lefser](/packages/lefser)
------

                       Changes in version 1.12.0                        

Significant user-visible changes

- The checkAbundances argument in lefser() checks that data are as
relative abundances and warns if otherwise (@LiNk-NY @sdgamboa, #28)
- relativeAb helper function available to convert data (@LiNk-NY)
- Deprecate the expr argument and use relab (short for relative
abundances)
- Add group labels to lefserPlot (@LiNk-NY #25, @asyakhl #31)
- 'Interoperating with phyloseq' section added to the vignette (#16)

[lfa](/packages/lfa)
---

                 Changes in version 2.1.11 (2023-06-20)                 

- Commented out excessive test for internal function .lreg against
glm, which differ more often than expected due to poor or lack of
convergence.
- Removed unused LaTeX package dependencies from vignette to prevent
errors restricted to specific testing platforms.

                 Changes in version 2.1.10 (2023-05-25)                 

- Version bump for bioconductor devel.

[limma](/packages/limma)
-----

                       Changes in version 3.58.0                        

- 
  New argument `covariate` for vooma(), which allows the variance
  trend to depend on extra covariates in addition to average
  log-expression.

- 
  New argument 'group' to simplify the typical calling sequence
  for removeBatchEffect().  The `group` argument provides an
  alternative and possibly simpler way to specify the
  experimental conditions to be preserved.

- 
  treat() now works when there are NA coefficients. Previously it
  gave an error if any estimated coefficient was NA.

- 
  Fix behaviour of topTable() with `adjust.method=NULL` so that
  it gives "BH" adjustment as documented. Previously the argument
  was passed to p.adjust(), which gave "holm" adjustment.

- 
  Convert Introduction vignette from Sweave/pdf from to
  Rmarkdown/html.  The Introduction vignette is also somewhat
  expanded.

- 
  cameraPR() has been addeed to the `10GeneSetTests.Rd` help
  topic.

- 
  The pdf limma User's Guide is now formally listed as a
  vignette.

- 
  Non-user-visible changes: limma now imports statmod rather than
  suggesting it.  BiocStyle has been added to the Suggests field
  of DESCRIPTION.

[maftools](/packages/maftools)
--------

                       Changes in version 2.18.0                        

NEW FUNCTIONS

- pathways and plotPathwaysfor summarizing & visualizing pathways
Issue: 956
- coGisticChromPlot for plotting two GISTIC objects side-by-side. PR
by biosunsci 954
- readGistic can take gistic output directory as an input. PR by
biosunsci 954

BUG FIXES

- Bug fixes while processing custom pathways
- Bug fix in oncoplot for drawing borders. 958
- Bug fix in plotSignatures for hardcoded axis limits. 949
- Bug fix in mafSurvival legend when samples argument is give. 937
- Bug fix in subsetMaf while handling only CNV events. 908
- Error handling when no deep/shallow CNV events found. 899
- Bug fix in oncoplot for duplicated values in gene list. 889

ENHANCEMENTS

- Added argument collapsePathway to oncoplot. Issue: 956
- Improved annovarToMaf with better handling of indels and
Variant_Type. Issue: 940
- Include absolute contribution of each signature in
extractSignatures
output. Issue: 939
- Added tsbToPIDs for custom names in oncoplot. Issue: Issue: 936
- Added DSEL protein to the database. Issue: 933
- Added MUC3A protein to the database. Issue: 932
- Added showOnlyPathway argument to oncoplot
- Added pathdb argument to PlotOncogenicPathways. Issue: 923
- Emit warnings when fishers test can not be performed during
somaticInteractions. Issue: 921
- Added leftMar and topMar arguments to somaticInteractions. Issue:
913
- Added toptBarLims argument to oncoplot. Issue: 910
- Added data argument to lollipopPlot function. Issue: 894
- Added sortByM1 and sortByM2 argument to coOncoplot. Issue: 888
- Added arguments leftBarVline, leftBarVlineCol, rightBarVline,
rightBarVlineCol topBarHline topBarHlineCol to oncoplot. Issue: 874
- Added revPal argument to somaticInteractions. Issue: 859
- Fix legend and color codes for numeric annotations in oncoplot.
Issue: 363

[magpie](/packages/magpie)
------

                        Changes in version 1.1.3                        

- Added a NEWS.md file to track changes to the package.

[MatrixGenerics](/packages/MatrixGenerics)
--------------

                       Changes in version 1.13.1                        

- Apply matrixStats 1.0.0 breaking change for useNames's default:
  All generic functions and methods defined in MatrixGenerics now
  set the useNames argument to TRUE instead of NA by default. See
  https://github.com/Bioconductor/MatrixGenerics/issues/31 for a
  discussion of this change.

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.9.1 (2023-07-19)                  

- fix bug in scree plot (subset the dimensionReduction object)

[matter](/packages/matter)
------

                       Changes in version 2.3.18                        

BUG FIXES

- Fixed bug in 'kdtree()' stack size allocation

- Cleaned up a few C++ compiler warnings

                       Changes in version 2.3.17                        

SIGNIFICANT USER-VISIBLE CHANGES

- Parameter 'filename' in 'struct()' is now 'path'

                       Changes in version 2.3.16                        

SIGNIFICANT USER-VISIBLE CHANGES

- Added 'checksum' method for character vectors (i.e., files)

                       Changes in version 2.3.15                        

BUG FIXES

- Fixed erroneous warning in 'matter_list' or 'matter_arr'
  constructors with existing files but no 'data' argument

                       Changes in version 2.3.14                        

SIGNIFICANT USER-VISIBLE CHANGES

- Added 'append' argument to matter constructors

- Simplified 'sparse_arr' implementation by removing
  shared index representation (which was not really used)

- When writing to an output file, 'chunkApply()' and friends
  now write entire chunks instead of each element

- Added coercions from 'matter_list' to 'matter_arr', etc.

BUG FIXES

- No longer truncate existing files when creating a matter object
  if nonexistent files are also included in the path

                       Changes in version 2.3.13                        

BUG FIXES

- Fixed translation from aliased C types to R types

                       Changes in version 2.3.12                        

NEW FEATURES

- Added 'image' method for formulas

                       Changes in version 2.3.11                        

SIGNIFICANT USER-VISIBLE CHANGES

- Changed specification of 'atoms' data types

- Use 'int16', int32', 'uint32', etc. for integer types

- Use 'float32' and 'float64' for floating point typess

- Added aliases so 'short', 'int', 'double', etc., still work

                       Changes in version 2.3.10                        

NEW FEATURES

- Added 'rowMaj' S4 generic and methods

SIGNIFICANT USER-VISIBLE CHANGES

- Remove defunct S4 generics

BUG FIXES

- Fixed bug in 'approx1' for "max" and "min" interpolation

                        Changes in version 2.3.9                        

NEW FEATURES

- Added 'vizi_pixels' for image plotting

- Added 'vizi_voxels' for 3D image plotting

SIGNIFICANT USER-VISIBLE CHANGES

- In estnoise_xxx() functions, renamed argument 'width' to 'n'

                        Changes in version 2.3.8                        

NEW FEATURES

- Added PLS (partial least squares)

- Added 'pls_nipals' (NIPALS)

- Added 'pls_simpls' (SIMPLS)

- Added 'pls_kernel' (kernel algorithms #1 and #2)

- Added OPLS (orthogonal PLS)

- Added 'opls_nipals' (NIPALS)

- Added FastMap projection with 'fastmap()'

- Added 'rowDists()' and 'colDists()'

- Added 'convolve_at()'

SIGNIFICANT USER-VISIBLE CHANGES

- Added regularization to NMF functions for stability

                        Changes in version 2.3.7                        

NEW FEATURES

- Added NMF (nonnegative matrix factorization)

- Added 'nnmf_mult()' (multiplicative updates)

- Added 'nnmf_als()' (alternating least squares)

- Exported 'prcomp_lanczos()'

SIGNIFICANT USER-VISIBLE CHANGES

- Accessing matter objects now checks for user interrupts

                        Changes in version 2.3.6                        

NEW FEATURES

- Added 'rowdist()' and 'coldist()'

- Added 'rowdist_at()' and 'coldist_at()'

- Added parallelization support for matrix multiplication

                        Changes in version 2.3.5                        

NEW FEATURES

- Added 'enhance_hist()' for histogram equalization

- Added 'enhance_adapt()' for CLAHE

BUG FIXES

- Fixed handling of missing values in smoothing filters

                        Changes in version 2.3.4                        

NEW FEATURES

- Added CWT based peak detection with
  'findpeaks_cwt()', 'findridges()' and 'cwt()'

- Added 'estnoise_quant()'

- Added 'filt1_conv()'

- Added 'filt2_conv()'

                        Changes in version 2.3.3                        

NEW FEATURES

- Added 'trans2d()' for affine transformations

- Added 'warp2_trans()' for transformation-based
  image registration using 'optim()'

                        Changes in version 2.3.2                        

NEW FEATURES

- More new signal processing features!

- Added 'approx2()' for 2D signal resampling
  and interpolation of scattered data

- Added 2D filtering including: 'filt2_ma()', 'filt2_gauss()',
  'filt2_bi()', 'filt2_adapt()', 'filt2_guide'

- Added nonlinear diffusion 'filt1_diff()' and 'filt2_diff()'

- Added the traditional Savitzky-Golay filter 'filt1_sg()'

- Added 'kdsearch()' and 'kdtree()' for K-dimensional searches

SIGNIFICANT USER-VISIBLE CHANGES

- Updated package DESCRIPTION

                        Changes in version 2.3.1                        

NEW FEATURES

- Lots of new signal processing features!

- Added 'approx1()' for signal resampling and interpolation

- Added 1D filtering including: 'filt1_ma()', 'filt1_gauss()',
  'filt1_bi()', 'filt1_adapt()', 'filt1_guide()', 'filt1_pag()'

- Added 1D warping and alignment including: 'warp1_loc()',
  'warp1_dtc()', 'warp1_cow()'

- Added continuum estimation including: 'estbase_loc()',
  'estbase_hull()', 'estbase_snip()', 'estbase_med()'

- Added local noise estimation including: 'estnoise_sd()',
  'estnoise_mad()', 'estnoise_diff()',
  'estnoise_filt()'

- Added peak processing including: 'findpeaks()',
  'peakwidths()', 'peakareas()', 'binpeaks()', 'mergepeaks()'

- Added 'downsample()' for signal and time series visualization

- Added 'vizi_plot' visualization methods

BUG FIXES

- Fix missing <cstdint> include for gcc 13.1.1 compatibility

[mbQTL](/packages/mbQTL)
-----

                 Changes in version 1.1.4 (2023-08-16)                  

- changes made on vignette

                 Changes in version 1.1.3 (2023-08-16)                  

- changes made on vignette

                 Changes in version 1.1.2 (2023-08-15)                  

- changes made on vignette

                 Changes in version 1.1.1 (2023-08-15)                  

- changes made on vignette and fixed bug

[metabCombiner](/packages/metabCombiner)
-------------

                       Changes in version 1.11.1                        

- m/z group size limit eliminated (previous: 10000 maximum)

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                         Changes in version 1.5                         

Changes in 1.5.9

- Addition of global function createStandardMixes

Changes in 1.5.8

- Fix .randomize_grouping to prevent collapsing of matrix when input
in a single column

Changes in 1.5.7

- Add function .group_standards_iteration and .randomize_grouping to
allow iteration through matrix of standards and group them if they
are dissimilar enough.

Changes in 1.5.6

- Fix issue in the vignette. Thanks @RemyDeB for the fix.

Changes in 1.5.5

- Update objects to the new definitions in Spectra version 1.11.10.

Changes in 1.5.4

- Add functions targetIndex and queryIndex to extract the indices of
the matched pairs query-target.
- Add examples and a section to the vignette explaining their use.

Changes in 1.5.3

- Add support to matchValues for matching between data.frame and
Spectra objects.

Changes in 1.5.2

- Fix vignette, examples and unit tests using QFeatures.
- Import query from AnnotationHub.

Changes in 1.5.1

- Add possibility to select the spectra variable for retention time
matching in matchSpectra (issue #98).

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                         Changes in version 1.9                         

MetaboCoreUtils 1.9.4

- Add function mclosest (issue #20).

MetaboCoreUtils 1.9.3

- isotopologues checks if provided m/z values are increasingly
ordered.

MetaboCoreUtils 1.9.2

- countElements returns NA for invalid elements instead of silently
dropping them ( PR #65).

MetaboCoreUtils 1.9.1

- countElements, subtractElements and addElements returns NA if an
input arguments is NA (issue #61, PR #62).

[metabolomicsWorkbenchR](/packages/metabolomicsWorkbenchR)
----------------------

                       Changes in version 1.11.1                        

- further vignette cache improvements

- use httptest to cache responses when building vignette

[metaseqR2](/packages/metaseqR2)
---------

                 Changes in version 1.13.2 (2023-10-04)                 

NEW FEATURES

- None.

BUG FIXES

- Removed the deprecated baySeq from PANDORA

[methylKit](/packages/methylKit)
---------

Version: 1.27.1
Category: IMPROVEMENTS AND BUG FIXES
Text: fix R cmd check warnings: - fix "Undocumented arguments in
        documentation object" by documenting ‘...’, ‘treatment’
        arguments for methylRawList and methylRawListDB constructors -
        methCall: replace variable length char with vector of chars to
        fix the warning: "variable length arrays are a C99 feature
        [-Wvla-extension]" - methCall: simplify the parsing of cigar
        string and replace calls to deprecated std::sprintf with
        std::snprintf

Version: 1.27.1
Category: IMPROVEMENTS AND BUG FIXES
Text: add test for processBismarkAln to check that reading sam and bam

Version: 1.27.1
Category: returns the same object for the same content
Text:

[mia](/packages/mia)
---

                         Changes in version 1.9                         

- loadFromMetaphlan: Bugfix, not all files include ID column.

- cluster: added wrapper for bluster's clusterRows function

- Added loadFromHumann

- calculateDMM: deprecated/updated outdated functions

- Added Tengeler2020 dataset

- *RDA & *CCA: calculate also statistical significance

- altExp support for meltAssay

- Deprecate mergeRows, mergeCols, agglomerateByRank,
  agglomerateByPrevalence

- Removed getAbundanceSample and getAbundanceFeature

- Updated test to avoid warning from deprecated functions

- Export mergeFeaturesByRank

- *RDA & *CCA: scores parameter for specifying output

- Improve mergeSEs and loadFromBiom

- Faith index: bugfix

[miaViz](/packages/miaViz)
------

                         Changes in version 1.9                         

- Updated plotDMN to work with newest mia version

- Added plotCCA and plotRDA functions

[MicrobiomeProfiler](/packages/MicrobiomeProfiler)
------------------

                        Changes in version 1.7.1                        

- Change maintainer from Meijun to Guangchuang

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.13.3                        

- fix the issue when assays is dgCMatrix other sparse matrix class.
(2023-09-12, Tue)
- using internal functions to convert dist object to tbl_df or tbl_df
to dist in mp_cal_dist and fix a bug of mp_extract_abundance when
rowData contains list metadata information. (2023-08-21, Mon)
- update the doc of rmun argument in mp_plot_abundance to avoid
misunderstanding. (2023-08-18, Fri)
https://github.com/YuLab-SMU/MicrobiotaProcess/issues/99
- update mp_import_humann_regroup() to keep the abundance of
contributed taxa in each sample with keep.contribute.abundance=TRUE.
(2023-08-15, Tue)
- use rlang::check_installed() to check if a suggested package is
installed, which will offer to install the package before
continuing. (2023-08-02, Wed)
- introduce order.by.feature argument in mp_plot_abundance to adjust
the sample order. (2023-07-24, Mon)
- update the format of citation and suppress the message introduced
by
tidytree. (2023-07-14, Fri)

                       Changes in version 1.13.2                        

- update mp_plot_ord to display the result of mp_adonis with
show.adonis = TRUE. (2023-06-21, Wed)
- using theme_blinds of ggfun. (2023-06-20, Tue)
- add mp_import_humann_regroup function to parsing the output of
humann_regroup_table. (2023-05-15, Mon)
- add fortify method for MPSE object. (2023-05-18, Thu)
- fix a bug of mp_plot_diff_res when ggnewscale updated to 0.4.9.
(2023-05-30, Tue)

                       Changes in version 1.13.1                        

- fix a bug for the abundance calculation with force = TRUE and
relative = FALSE. (2023-04-28, Fri)

[MICSQTL](/packages/MICSQTL)
-------

                       Changes in version 0.99.17                       

- Fix file name in vignette.

MICSQTL 0.99.16

- Enhance PGD to handle a wider range of omics data inputs.

                       Changes in version 0.99.15                       

- Update deconv and ajive by integrating a novel feature selection
approach and improved deconvolution method.

                       Changes in version 0.99.14                       

- Update deconv by allowing proteins selected by loadings.

                       Changes in version 0.99.13                       

- Update ajive_decomp by adding option to output loadings.

                       Changes in version 0.99.12                       

- Update package loading in vignette.

                       Changes in version 0.99.11                       

- Add @return to data documentation.

                       Changes in version 0.99.10                       

- Merge cns_plot to ajive.
- Remove extra space.
- Downsize data in example code.

                       Changes in version 0.99.9                        

- Update email address.

                       Changes in version 0.99.8                        

- Remove unused lines and update importFrom.
- Add unit tests.
- Update abstract and intro.

                       Changes in version 0.99.7                        

- Add iteration option in deconv function.

                       Changes in version 0.99.6                        

- Remove unwanted lines.
- Update documentation.

                       Changes in version 0.99.5                        

- Update pkg loading.

                       Changes in version 0.99.4                        

- Reformat function and documentation.

                       Changes in version 0.99.3                        

- Resolve errors in different OS systems.
- Fix unrecognized unicode.

                       Changes in version 0.99.2                        

- Add a NEWS file to track changes to the package.

[mistyR](/packages/mistyR)
------

                         Changes in version 1.9                         

                        Changes in version 1.8.1                        

- Requires dplyr >=1.1.0 due to change of sorting order in arrange()

[MoleculeExperiment](/packages/MoleculeExperiment)
------------------

                 Changes in version 1.1.4 (2023-10-06)                  

- show method has been improved, and package is less verbose.

                 Changes in version 1.1.3 (2023-08-22)                  

- feature_name is not called feature_id for consistency with
segment_id.

                 Changes in version 1.1.2 (2023-07-26)                  

- IMPORTANT: readCosmx and readMerscope can now also handle reading
in
and standardising boundary information!!

[monaLisa](/packages/monaLisa)
--------

                        Changes in version 1.7.1                        

- allow modification of heatmap graphical parameters by forwarding
...
argument in plotMotifHeatmaps to all calls to
ComplexHeatmap::Heatmap

[Moonlight2R](/packages/Moonlight2R)
-----------

                       Changes in version 0.99.14                       

Summary

- made test for GRN and arm64 less stringent to account for
architectural differences

                       Changes in version 0.99.13                       

Summary

- refactored examples to make them faster

                       Changes in version 0.99.12                       

Summary

- added data loadings in functions

                       Changes in version 0.99.11                       

Summary

- added missing data loadings in vignette
- changed handling of null cases in PRAToTibble

                       Changes in version 0.99.10                       

Summary

- updated code style in all functions

                       Changes in version 0.99.9                        

Summary

- fixed hardcoded plot title in plotFEA
- updated reference in docs to published paper
- fixed excessively long example line in GLS function
- moved globalVariables calls after function definitions
- turned LazyData to false and issued connected fixes

                       Changes in version 0.99.8                        

Summary

- replaced number sequences generation (e.g. using seq(n) instead of
1:n) in all functions

                       Changes in version 0.99.7                        

Summary

- added bindings for global variables in majority of functions

                       Changes in version 0.99.6                        

Summary

- switched \dontrun to \donttest in some examples

- fixed vignette to have fewer eval=FALSE chunks

                       Changes in version 0.99.5                        

Summary

- added several checks for correctness of main function arguments

                       Changes in version 0.99.4                        

Summary

- changed instances of sapply to vapply

                       Changes in version 0.99.3                        

Summary

- added tests with testthat

- updated following example data: dataFEA, dataGRN, dataURA, dataPRA
and cscape_somatic_output

- added following example data: dataURA_plot, dataGRN_no_noise

                       Changes in version 0.99.2                        

Summary

- added GLS (Gene Literature Search) function

- fixed library problems with vignettes

                       Changes in version 0.99.1                        

Summary

- removed package documentation from Rdata.R

                       Changes in version 0.99.0                        

Summary

- first release of Moonlight2R.

New features (added or significantly changed respect to MoonlightR)

- DMA Driver mutation analysis

- plotDMA Creates one or more heatmap of the output from DMA

- plotMoonlight Creates heatmap of Moonlight Gene Z-scores for
selected genes

- EncodePromoters Experimentially verified promoter sites

- LOC_protein Level of consequence protein

- LOC_translation Level of consequence translation

- LOC_transcription Level of consequence transcription

- NCG Network of Cancer Genes 7.0

- dataPRA output from PRA function

- dataMAF Mutation data from TCGA-LUAD

- dataDMA Output from DMA function

- cscape_somatic_output Cscape-somatic annotations of TCGA-LUAD

- DEG_Mutations_Annotations Differentially expressed genes's
Mutations

- Oncogenic_mediators_mutation_summary Oncogenic Mediators Mutation
Summary

- moonlight Function to run moonlight pipeline

[MoonlightR](/packages/MoonlightR)
----------

                       Changes in version 1.26.1                        

BUG FIX

- changed arguments for downloading datasets following updates to
  TCGAbiolinks

- added suggested package to fix install crash
  FIRST VERSION - FEATURES
  * dataFilt Gene Expression (Rnaseqv2) data from TCGA LUAD
  * dataGRN GRN gene regulatory network output
  * dataURA Output example from function Upstram Regulator Analysis
  * DEGsmatrix DEG Differentially expressed genes
  * DiseaseList Information on 101 biological processes
  * DPA DPA
  * EAGenes Information about genes
  * FEA FEA
  * GDCprojects Information on GDC projects
  * geneInfo Information about genes for normalization
  * GEO_TCGAtab Information on GEO data (and overlap with TCGA)
  * getDataGEO getDataGEO
  * getDataTCGA getDataTCGA
  * GRN Generate network
  * GSEA GSEA
  * knownDriverGenes Information on known cancer driver gene from
  COSMIC
  * listMoonlight Output list from Moonlight
  * LPA LPA
  * moonlight moonlight pipeline
  * MoonlightR MoonlightR
  * PEA PEA
  * plotCircos plotCircos
  * plotFEA plotFEA
  * plotNetworkHive plotNetworkHive: Hive network plot
  * plotURA plotURA: Upstream regulatory analysis heatmap plot
  * PRA Pattern Recognition Analysis (PRA)
  * tabGrowBlock Information growing/blocking characteristics for 101
  selected biological processes
  * URA URA Upstream Regulator Analysis

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.45.1                        

- Fix the issue of 'This picture was not generated by Cairo graphics'
  for importSVG.

[msa](/packages/msa)
---

                       Changes in version 1.33.2                        

- update of Makevars: added -lpthread to PKG_LIBS in order to make sure
  that
  package also builds correctly on Bioconda

                       Changes in version 1.33.1                        

- update of msaConsensusSequence() and msaConsensusSequence() methods
  to
  account for recent change in function Biostrings::consensusMatrix()

                       Changes in version 1.33.0                        

- new branch for Bioconductor 3.18 devel

[MSA2dist](/packages/MSA2dist)
--------

                 Changes in version 1.5.2 (2023-05-24)                  

Major changes

Minor improvements and bug fixes

- added aa2selfscore

                 Changes in version 1.5.1 (2023-05-22)                  

Major changes

Minor improvements and bug fixes

- additional option to use asymmetric score matrix (symmetric=FALSE)
- changed rcpp_distSTRING.cpp
- changed rcpp_pairwiseDeletionAA.cpp
- changed rcpp_pairwiseDeletionDNA.cpp
- changed aastring2dist.R
- changed dnastring2dist.R

[MsBackendSql](/packages/MsBackendSql)
------------

                         Changes in version 1.1                         

Changes in 1.1.5

- Add dbconn methods for MsBackendSql and MsBackendOfflineSql.

Changes in 1.1.4

- Improve performance of createMsBackendSqlDatabase by using also
parallel processing for the peaksData call.

Changes in 1.1.3

- Add support for setBackend to MsBackendOfflineSql.

Changes in 1.1.2

- Mention in documentation that MsBackendSql can not be saved to
disk.
- Expand vignette adding related documentation.

Changes in 1.1.1

- Fix for filterRt avoiding to filter if range is infinite.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                        Changes in version 1.13                         

MsCoreUtils 1.13.1

- Add functions entropy and nentropy.

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.27                         

MSnbase 2.27.1

- Fix declarations of centroided/smoothed setters for OnDiskMSnExp
objects (from Hervé Pagès via Github).

MSnbase 2.27.0

- New devel

[msPurity](/packages/msPurity)
--------

                       Changes in version 1.27.1                        

- createMSP fix - now uses the median precursor MZ and precursor RT in
  the MSP file

[MsQuality](/packages/MsQuality)
---------

                         Changes in version 1.1                         

Changes in version 1.1.3

- update tests after update of OBO file / rmzqc

Changes in version 1.1.2

- Fix implementation of mzAquisitionRange
- add numberEmptyScans in qualityMetrics function
- add unit tests for export in rmzqc format
- add interpretation aid of metrics in vignette

Changes in version 1.1.1

- move msdata from Suggests to Imports
- rename function rtDuration to chromatographyDuration
- rename function rtOverTicQuantiles to ticQuartersRtFraction
- create function numberEmptyScans
- add attributes (MS QC terms) to the output of the Spectra metrics
functions if the output matches the described term
- add rmzqc to IMPORTS
- add functionality to export quality metrics as in rmzqc format
- adjust documentation to newest version of PSI MS CV obo file
- add to the vignette information on how the metrics are calculated
- add argument filterEmptySpectra to remove entries of length 0 or
that have intensity 0, implement the argument in the functions
calculateMetricsFromOneSampleSpectra, calculateMetricsFromSpectra,
calculateMetricsFromMsExperiment, and calculateMetrics

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.28.0                        

New features

- Dropped experiments are no longer kept in the metadata slot. They
can be seen with drops() (@LTLA, #323).

Bug fixes and minor improvements

- Checking colnames in sampleMap vs ExperimentList is more robust by
only comparing unique and sorted values in each.

[MultiRNAflow](/packages/MultiRNAflow)
------------

                       Changes in version 0.99.9                        

Suggestions and remarks from a bioconductor team member

- modification of outputs in order to replace list ouputs by SE class
object
- Vignette (Running_analysis_with_MultiRNAflow.Rmd)
- modification because ouputs modification
- More unit tests (now 197 unit tests in 31 files), coverage = 60/100

                       Changes in version 0.99.8                        

Suggestions and remarks from a bioconductor team member

- Correction of DESCRIPTION file
- add a new function to improve inputs in others function
- modification of R functions depending on previous inputs
- Vignette (Running_analysis_with_MultiRNAflow.Rmd)
- code-style to highlight (function, variable, package names)
- keeping includegraphics only for the introduction part
- More unit tests (now 45 unit tests in 25 files)
- add the R file MultiRNAflow-package.R for man folder

                       Changes in version 0.99.7                        

Suggestions and remarks from a bioconductor team member

- Correction of DESCRIPTION file
- Vignette (Running_analysis_with_MultiRNAflow.Rmd)
- code-style to highlight (function, variable, package names)
- More details about the vignette
- eval=FALSE replaced by eval=TRUE
- More unit tests

                       Changes in version 0.99.6                        

Bug fixes

- Necessary correction by the automated single package builder of
bioconductor.org.

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[multiWGCNA](/packages/multiWGCNA)
----------

                       Changes in version 0.99.4                        

- Added BiocStyle to suggests

                       Changes in version 0.99.3                        

- Added suggested importFroms in R check
- Added Coexpression Line graph function documentation and to
vignette

                       Changes in version 0.99.2                        

- All ERRORS and WARNINGS from R CMD CHECK and BiocCheck have been
addressed

                       Changes in version 0.99.1                        

- Added a NEWS.md file to track changes to the package.

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.9.19                        

New features

- infer_eff_direction parameter added so user can decide whether to
run the check

Bug fix

- Typo in unit test for infer effect direction.
- IEU GWAS unit tests updated to account for server outages.

                       Changes in version 1.9.18                        

Bug fix

- Fixed column header mappings
- Made all uncorrected header names uppercase and removed
duplicates
- "TOTALSAMPLESIZE" now maps to "N" instead of "NSTUDY"
- "MAJORALLELE", "MAJOR_ALLELE", "MAJOR-ALLELE", and "MAJOR
ALLELE" now map to "A1" instead of "A2"
- Removed the mappings for "OR-A1", "OR.A1", "OR_A1", and "BETA1"
because MSS assumes that A2 is the effect allele
- Removed mappings for "A1FREQ", "A1FRQ", "AF1",
"FREQ.A1.1000G.EUR", "FREQ.A1.ESP.EUR",
"FREQ.ALLELE1.HAPMAPCEU", "FREQ1", "FREQ1.HAPMAP", and "FRQ_A1"
because MSS defines "FRQ" to be the allele frequency of A2
- Removed mappings for "CHR36", "BASE_GRCH36", "POSITION36",
"POSGRCH36", "BASEGRCH36", "POS36", "POS GRCH36", "POS.GRCH36",
"POS-GRCH36", and "POS_GRCH36" because MSS does not support the
GRCh36 genome build
- Removed the ambiguous mapping "NMISS" -> "N" because "NMISS" can
refer to the number of samples with missing data
- Removed the ambiguous mapping "WEIGHT" -> "N" because "WEIGHT"
can refer to coefficient weights
- Fixed inference of allele where ambiguous (A1, A2) naming used (see
infer_effect_column.R for code) but in short:
- Three checks now made to infer which allele the effect/frequency
information relates to. See infer_effect_column.R for further
details.
- See get_eff_frq_allele_combns.R for how effect/frequency columns
that infer the allele are captured in the mapping file

New features

- New column header mappings:
- "VARIANT_ID" and "RSIDS" --> "SNP"
- "P_BOLT_LMM" --> "P"
- "NCASES" --> "N_CAS"
- "N_EFFECTIVE", "N_INFORMATIVE", and "TOTAL_N" --> "N"
- "HET_P" --> "HETPVAL"
- "HET_ISQ" --> "HETISQT"
- "ALL_AF" --> "FRQ"
- "DIRECT" --> "DIRECTION"
- "ALT_EFFSIZE" --> "BETA"
- "INFORMATIVE_ALT_AC" --> "AC"

                       Changes in version 1.9.17                        

Bug fix

- Cases checking ref genome where there are no indels would sometimes
cause an error when joining. This resolved this issue.

                       Changes in version 1.9.16                        

New features

- flip_frq_as_biallelic parameter added enabling frequencies of
non-bi-allelic SNPs to be flipped as if they were bi-allelic (1 -
frequency) i.e. ignoring the frequencies of other alternative
alleles (assuming these will be negligible). Note this will not be
done as default as it is not fully correct but may be useful for
some users.

                       Changes in version 1.9.15                        

Bug fix

- Fix for imputation column when imputing RS ID from CHR:BP. Avoids
crash and ensures correct identification of imputed SNPs.
- Avoid running compute_nsize function when no imputation is wanted
by
user - also avoids message output in this situation.

                       Changes in version 1.9.14                        

Bug fix

- Fix reporting of genome-wide sign variants before formatting.

                       Changes in version 1.9.13                        

Bug fix

- In check_bp_range ensure that the BP column is numeric.

                       Changes in version 1.9.12                        

Bug fix

- In check_no_rs_snp the order of operations had to be reversed to
ensure all values were present before sorting column headers when
imputation_ind=TRUE and imputing rsIDs.

                       Changes in version 1.9.11                        

New features

- The rmv_chrPrefix parameter in format_sumstats() has been replaced
with the new chr_style parameter, which allows users to specify
their desired chromosome name style. The supported chromosome styles
are "NCBI", "UCSC", "dbSNP", and "Ensembl" with "Ensembl" being the
default.
- check_chr() now automatically removes all SNPs with nonstandard CHR
entries (anything other than 1-22, X, Y, and MT in the Ensembl
naming style).

                       Changes in version 1.9.10                        

Bug fix

- Better method to detect vcf files - looks for vcf in extension not
in name.

                        Changes in version 1.9.9                        

Bug fix

- Check ref genome change - if not match found for either genome
build, an error will now be thrown.
- Checks has been added so that if chrom col has chr as a prefix,
this
will be removed before testing genome build.

                        Changes in version 1.9.8                        

Bug fix

- Bug fix when using imputation_ind with NA in chr column.

                        Changes in version 1.9.7                        

New features

- ignore_multi_trait parameter added which will ignore any
multi-trait
p-values if set to TRUE. By default it is false to maintain the
current default running conditions for MSS.

[muscat](/packages/muscat)
------

                       Changes in version 1.15.1                        

- bug fix in 'pbDS': too stringent filtering causing no genes in any
  clusters
  to be tested previously resulted in a 'subscript out of bounds'
  error;
  execution is stopped and an informative error thrown instead.

- bug fix in 'mmDS': 'dream' (new version?) wouldn't recognize
  model variables provided as data; fixed via adding 'as.formula()'.

- "analysis" vignette: replaced suspended 'dplyr' function
  'top_n' with 'slice_min' when filtering for top DS hits;
  fixed some typos; updated preprint to journal reference.

                       Changes in version 1.15.0                        

- Bioconductor release v3.17

[mzR](/packages/mzR)
---

                       Changes in version 2.35.1                        

- fix compilation on Fedora 38 / R-4.3.0. Thanks to Christian Iseli!
  Closes #282

[NanoTube](/packages/NanoTube)
--------

                        Changes in version 1.7.2                        

- runLimmaAnalysis now allows optional arguments, which are passed to
limma::lmFit.

                        Changes in version 1.7.1                        

- The codeclass.retain option now allows runLimmaAnalysis() to be run
using a CodeClass/CodeClasses specified by the user, instead of
automatically removing non-endogenous genes. See
help(runLimmaAnalysis) for details.
- processNanostringData() can now handle a vector of .rcc files, in
addition to the previous options for loading NanoString data.

[ngsReports](/packages/ngsReports)
----------

                        Changes in version 2.3.5                        

- added support for rnaseqc metrics files

[nipalsMCIA](/packages/nipalsMCIA)
----------

                 Changes in version 0.99.7 (2023-10-07)                 

Minor Changes

- Updated readme to reflect MAE changes.
- Updated citation.
- Fixed bug in documentation for nipals_multiblock.

                 Changes in version 0.99.6 (2023-09-09)                 

Major changes

- Switched primary input to nipals_multiblock() to a
MultiAssayExperiment object.
- nipals_multiblock() now outputs an object of the NipalsResult
class.
- Converted all downstream analysis functions to work with the
NipalsResult class.

Minor improvements and bug fixes

- Added simple_mae() function to convert a list of dataframes to a
MultiAssayExperiment object.
- Fixed missing \value fields in man page for NipalsResult class.

                 Changes in version 0.99.5 (2023-06-22)                 

Major changes

- Bumping version number to trigger re-build following bug in
BiocCheck.

                 Changes in version 0.99.4 (2023-06-01)                 

Major changes

- Fixed data corruption issues from v0.99.3

                 Changes in version 0.99.3 (2023-06-01)                 

Major changes

- Changed the eigenvalue calculation in nipals_iter() to compute the
variance of the global score at each deflation step. Prior versions
used an SVD method to compute the singular values of the deflated
data matrix directly.
- Fixed bug in projection_plot() where there was a mismatch between
color labels and plotting order.
- Added parameter to nipals_multiblock that specifies whether the
samples are in the rows or the columns

Minor improvements and bug fixes

- Changed vignette styling from rmdformats to BiocStyle and added
installation sections to all of the vignettes.
- Removed empty helper.R
- Added significantly more unit testing.
- Fixed bug in projection_plot() when metadata was provided but no
color_col was selected.
- Renamed the associated output of col_preproc_method in
nipals_multiblock. The metadata field is also now available in the
output independent of whether metadata is provided in the input.
- Added checks for consistency in sample names across data blocks and
metadata.

                 Changes in version 0.99.2 (2023-03-25)                 

Major changes

- Shrank the vignettes sizes (especially Vignette 2).

Minor improvements and bug fixes

- Restructured Vignette 2 to be more streamlined and have more
explanations.
- Add an additional single cell data file to the repository using
piggyback.

                 Changes in version 0.99.1 (2023-02-26)                 

Major changes

- Included support for MultiAssayExperiment in nipals_multiblock.
- Improved access to existing data objects.

Minor improvements and bug fixes

- Made get_colors() more flexible for different color palette
options.

                 Changes in version 0.99.0 (2022-10-21)                 

- Added single cell data to the repository using piggyback.

[nnSVG](/packages/nnSVG)
-----

                 Changes in version 1.5.3 (2023-05-30)                  

- use model formula without intercept for non-spatial model - to enable
  weighted model

[NormalyzerDE](/packages/NormalyzerDE)
------------

                       Changes in version 1.19.7                        

- Remove RCmdrMisc as dependency, resolving a fatal crash
  caused by updated signature of RCmdrMisc CV calculation function

- Updating code to remove deprecation warnings for ggplot2 code

[oncoscanR](/packages/oncoscanR)
---------

                        Changes in version 1.4.0                        

Versions on bioconductor and Github have been merged and, from now
on,
will have consistent version numbers.

[OncoSimulR](/packages/OncoSimulR)
----------

                 Changes in version 4.3.3 (2023-07-17)                  

- One test failing on Linux aarch64: fixed
  (by commenting out, which is the right thing to do).

                 Changes in version 4.3.2 (2023-06-03)                  

- Allow passing a seed to MAGELLAN's random fitness landscapes
  in rfitness function.

- Fixed bug in rfitness when using a 3-element vector for scaling,
  and one or more fitness values had the same value as WT.

[Organism.dplyr](/packages/Organism.dplyr)
--------------

                       Changes in version 1.30.0                        

BUG FIXES

- src_ucsc() failed to correctly handle new 'hs1' resources (T2T
  genomes) for 'human'

[orthogene](/packages/orthogene)
---------

                        Changes in version 1.7.2                        

Bug fixes

-

                        Changes in version 1.7.1                        

New features

- remove_all_nas
- Can handle multiple cols.
- plot_orthotree
- New arg clades_rotate
- New func: rotate_clades

Bug fixes

- add_columns
- Handle both vectors and columns.
- sort_rows_func
- Handle both vectors and columns.
- filter_gene_df
- Avoid coercing single-col dataframe into vector.
- Flagged in #34
- Fix test-report_orthologs
- Recognize either Gene.symbol or input_gene cols depending on
when ortholog conversion was done.
- Fix test-convert_orthologs
- Line 99 test had wrong number of cols.
- map_genes_planosphere
- Add backup download strategy.

                        Changes in version 1.7.0                        

New features

- Bump version.

[orthos](/packages/orthos)
------

                       Changes in version 0.99.4                        

- Temporary inactivate failing Windows tests

                       Changes in version 0.99.3                        

- Pinning more packages

                       Changes in version 0.99.2                        

- Pinning more packages
- updated stick with Bioconductor url

                       Changes in version 0.99.1                        

- make all function examples runnable
- added scheduled cron runs on R-CMD-check.yaml
- added inst/scripts to describe extdata generation

                       Changes in version 0.99.0                        

- prepare package for submission to Bioconductor.

                        Changes in version 0.1.0                        

- Added a NEWS.md file to track changes to the package.

[OutSplice](/packages/OutSplice)
---------

                 Changes in version 1.0.1 (2023-06-21)                  

- Fixed issue where nonsignificant events would be listed in
  ASE.types when no under-expressed outliers were found

[pairedGSEA](/packages/pairedGSEA)
----------

                        Changes in version 1.1.2                        

- Fix bug where DEXSeq fails if there are NAs in the metadata

                        Changes in version 1.1.1                        

- Improved how surrogate variables are transferred in the design
formulas

                        Changes in version 1.1.0                        

- paired_ora now adjusts gene background for analysis to reduce bias
- paired_ora now also does an ora combining the genes from the two
analyses
- Added plotting mode paired = TRUE to plot_ora for the new
paired_ora
analysis type
- Added baseMean to gene-level aggregation output

[partCNV](/packages/partCNV)
-------

                       Changes in version 0.99.2                        

Changes

- Added unit tests.

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.15.3 (2023-10-03)                 

- correct `plotAnnotationDiagnosticMultiplot` axes limit after change
  in ggplot2 behaviour following a rotation

                 Changes in version 1.15.2 (2023-10-02)                 

- ROI, uROI and FIR input checks for NA in rtMin, rtMax, mzMin and
  mzMax

                 Changes in version 1.15.1 (2023-08-06)                 

- corrections in vignettes

[phantasus](/packages/phantasus)
---------

                       Changes in version 1.21.3                        

- Advanced normalizations: TMM & voom

- Complex differential expression design for limma & DESeq2

- Rework of count meta files (breaks backward compatibility)

[phenomis](/packages/phenomis)
--------

                        Changes in version 1.3.8                        

- writing (SummarizedExperiment & MAE): metadata.l argument now
  included

                        Changes in version 1.3.6                        

- writing (SummarizedExperiment): now saves the metadata as an
  additional _metadata.rds file

                        Changes in version 1.3.4                        

- hypotesting: minor documentation update

                        Changes in version 1.3.2                        

- UTF-8 file encoding for writing (utils::write.table)

[PhyloProfile](/packages/PhyloProfile)
------------

                       Changes in version 1.16.0                        

- extended options in config file

- option to specify IP and port

- able to generate domain plot from main/customized profile tab

- highlight multiple genes/taxa; division lines for taxon group
  (#110)

- simplify gene IDs in detailed plot; show taxonomy hierarchy
  (#110)

- option to filter features in architecture plot (#110)

- option to show all input taxa (#110)

- added more info to domain plot (e-value, bit-score, pHMM)

- option to specify host and port for runPhyloProfile() fn

- fixed bug when lowest rank is subspecies

                       Changes in version 1.14.5                        

- improved highlighted duplicated ortho IDs

                       Changes in version 1.14.4                        

- fixed bug ordering genes

- auto identify lowest rank

                       Changes in version 1.14.2                        

- option to use user defined taxonomy DB (#124)

                       Changes in version 1.14.1                        

- option for highlight duplicated ortholog IDs

- option for cluster profiles based on ortholog IDs

- speedup by pre-calculating taxonomy tree (#123)

- disable rank selection after plotting

- turn on profile clustering by default

- changed "ftp://" to "https://" for taxdmp.zip

- fixed #126

[Pigengene](/packages/Pigengene)
---------

                Changes in version 1.27.16 (2023-06-16)                 

Changes in existing functions

- Scaling is done in the compute.pigengene() function to avoid
  the following error that Sogand detailes on her lano: Error in
  La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'

                 Changes in version 1.27.2 (2023-05-12)                 

Changes in existing functions

- The compute.pigengene() function does not warn about almost
  constant genes.

[planet](/packages/planet)
------

                        Changes in version 1.9.2                        

                        Changes in version 1.9.1                        

[plasmut](/packages/plasmut)
-------

                       Changes in version 0.99.7                        

SIGNIFICANT USER-VISIBLE CHANGES

- Initial plasmut package release:

- Adds importance sampling method for Bayesian inference through
  importance_sampler()

[plotgardener](/packages/plotgardener)
------------

                        Changes in version 1.7.7                        

BUG FIXES

- params parsing removes non-standard chromosomes when identifying
gene transcripts

                        Changes in version 1.7.6                        

BUG FIXES

- colorby logic parsing is patched.

                        Changes in version 1.7.5                        

NEW FEATURES

- plotPairs, plotPairsArches, and plotRanges allow NA fill.

                        Changes in version 1.7.4                        

BUG FIXES

- Removed extra page creation with annoYaxis axisLine = FALSE.

                        Changes in version 1.7.3                        

BUG FIXES

- Removed blank page creation with pdf() calls for all major plotting
functions.

                        Changes in version 1.7.2                        

BUG FIXES

- Fixed plotMultiSignal width and height parsing bug.

                        Changes in version 1.7.1                        

NEW FEATURES

- plotManhattan y-scales can be reversed for Miami plot-style
layouts.

                        Changes in version 1.7.0                        

Version bump for Bioconductor 3.17 release.

[plyinteractions](/packages/plyinteractions)
---------------

                       Changes in version 0.99.6                        

- Added validity functions for new classes
- Added tests (~90% coverage at this point)
- Improved doc:
- Fixed lines with > 80 characters
- Moved generics and classes to separate R files

                       Changes in version 0.99.5                        

- Initiated NEWS.md
- Internal refactoring:
- Removed paste calls in condition signals
- Moved classes and generic definitions in AllClasses.R and
AllGenerics.R
- Fixed R minimum working version (4.3.)
- LazyData set to FALSE following Bioc recommendations
- Added package-level documentation
- Improved vignettes

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.9.3 (2023-05-15)                  

- Added PoDCall publication citation

                 Changes in version 1.9.2 (2023-05-10)                  

- Bugfixes

                 Changes in version 1.9.1 (2023-05-03)                  

- Made PoDCall compatible with new BioRad software QX Manager

[polyester](/packages/polyester)
---------

                       Changes in version 1.99.3                        

- NB function now exported

- note that version 1.99.3 on GitHub was version 1.1.0 on Bioconductor.

                       Changes in version 1.99.2                        

- bug fix in fragment generation (last 2 bases of transcript were never
  sequenced)

[pRolocGUI](/packages/pRolocGUI)
---------

                        Changes in version 2.11                         

CHANGES IN VERSION 2.11.0

- New version for Bioc devel

CHANGES IN VERSION 2.11.1

- Fix bug in DT table when fData column is a matrix see issue #117

[ProtGenerics](/packages/ProtGenerics)
------------

                       Changes in version 1.33.1                        

- rename uniqueMsLevel() to uniqueMsLevels()

[PureCN](/packages/PureCN)
------

                        Changes in version 2.8.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Make processMultipleSamples temporarily defunct because the
  copynumber package was removed from Bioconductor

- Make it possible to specify saveRDS version to make output files
  readable by old R versions prior to 3.6.0 (#255)

BUGFIXES

- Fixed an issue with callLOH and --model-homozygous (#254)

- Fixed crash when VCF contained NAs in base quality scores (#249)

- Fixed wrong check for outdir write permissions in Coverage.R and
  NormalDB.R (#258)

- Fixed inverted return codes in couple of scripts (#284)

- Fixed an issue in callAlterations where the id argument was largely
  ignored (#292)

- Fixed broken support for GenomicsDB-R from their developer branch
  (#296)

- Fixed crash in plotAbs (#260)

- Fixed an issue with gene-level calls when annotation contained
  non-official symbols found on multiple chromosomes (#298)

- Fixed a wrongly formatted error message when no germline database
  information was found (#302)

- Fixed a crash when DB field in VCF only contains NAs(#301)

- Added --min-base-quality argument for PureCN.R (#320)

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.11                         

QFeatures 1.11.2

- Update message to fix test upon recent changes in MAE.

QFeatures 1.11.1

- Update nNA() and filterNA() and man pages to clarify percentages
and
proportions (see #189).

QFeatures 1.11.0

- New Bioconductor 3.18 (stable) release

[qsvaR](/packages/qsvaR)
-----

                        Changes in version 1.5.3                        

BUG FIXES

- Fixed the error messages displayed by k_qsvs() to handle different
types of situations. We implemented this update with @HediaTnani,
@reneegf, and @lahuuki.

                        Changes in version 1.5.2                        

BUG FIXES

- Fixed a bug in qSVA() which was not passing sig_transcripts to
getDegTx(). Related to
https://github.com/LieberInstitute/qsvaR/issues/29.
- Fixed the documentation to highlight when users should use
set.seed() to ensure the reproducibility of their results. Related
to https://github.com/LieberInstitute/qsvaR/issues/28.

                        Changes in version 1.5.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Hedia Tnani is now the maintainer of qsvaR.

[raer](/packages/raer)
----

                       Changes in version 0.99.12                       

- Changes made to prepare for bioc submission

                       Changes in version 0.99.11                       

- Replaced base R fisher test with c-wrapper to call htslib fisher
test, which speeds up execution many fold.

                       Changes in version 0.99.10                       

- The options to write to tabix indexed output files have been
removed
from pileup_sites() as they have limited utility and introduce
unwanted code complexity.

                       Changes in version 0.99.9                        

- The genomic-unstranded option for the library-type argument in
FilterParam() has been renamed to unstranded, and the unstranded
option has been removed.

                       Changes in version 0.99.8                        

- Function arguments involving a fasta file have been renamed to all
be fasta

                       Changes in version 0.99.7                        

- added a single cell specific AEI calculation (calc_scAEI())

                       Changes in version 0.99.6                        

- added method to count base consensus base when counting UMIs with
pileup_cells() using the sum of base qualities to select consensus.

                       Changes in version 0.99.5                        

- pileup_cells() now allows for multiple alleles to be queried at a
site.

- Fixed an indexing bug in pileup_cells() that misassigned sites to
counts.

                       Changes in version 0.99.4                        

- annot_snps will now compare the site allele to the SNP allele and
return a column snp_matches_site indicating if the site matches the
snp.

- added new function, find_scde_sites() to identify differentially
editing sites in single cell data using fishers exact tests.

- pileup_cells now respects the min_depth and min_variant_reads
FilterParameters.

                       Changes in version 0.99.3                        

- support BamFile and BamFileList inputs to pileup_sites() and
pileup_cells(), which provides an option to provide custom BAI index
file names.

                       Changes in version 0.99.2                        

- rename prep_for_de() and perform_de() to make_de_object() and
find_de_sites().

                       Changes in version 0.99.1                        

- default values for edit_from and edit_to for calc_edit_frequency()
have been changed to A and G respectively.

- renamed type argument in perform_de to test and removed type
argument in prep_for_de

                       Changes in version 0.99.0                        

- added support for processing multiple BAM files with calc_AEI().

- Dropped minimally used bad_reads and reads parameters from
pileup_sites()

- Added utility to screen scRNA-seq bam files for regions with
oligo-dT mispriming (find_mispriming_sites()).

- add option to query ref and alt SNP alleles

- added tests for SummarizedExperiment filtering approaches

- added a strand bias stat sor using approach from GATK
(StrandOddsRatio), and a confidence score calc_confidence() from
SAILOR pipeline.

- 'N' bases in read or reference are ignored

- Removed outdated or unused functionality:

- bed indexing (indexBed and related C code)
- bam tag indexing (build_tag_index, show_tag_index, get_tag_bam,
)
- bam tag index based single cell approach (sc_editing)
- bam tag indexing C code from bri (src/bri/*)
- sparse matrix merging for merge_pileups().
- unneeded utilities (filter_by_coverage)
- Remaining (and mostly unused) Rcpp code
- Removed fastmap, Rcpp, zlibbioc, RColorBrewer, and BiocGenerics
dependencies
- Removed system requirements for C libraries used by bri

- The bed indexing used in pileup_sites() has been replaced with the
region indexing approach from pileup_cells().

- pileup_sites() now requires a GRanges object rather than a bed
file.
The bedfile parameter has been removed and replaced with a sites
parameter.

- Renamed Ref and Var output columns to REF and ALT and nVar was
renamed to nAlt. This provides consistency with VCF format and
consistency across pileup_cells() and pileup_sites() function calls

- pileup_cells() gained functionality to process multiple smart-seq2
style bam files.

- Changed filterParam argument in pileup_sites and pileup_cells to
param for simplicity.

- Added FilterParam to exclude multi-allelic sites
report_multiallelic, or exclude reporting a variant in the Var assay
based on allelic frequency (min_allelic_freq).

- The bam_flags parameter used in pileup_sites and pileup_cells has
been moved into the FilterParam class.

- The bedindex parameter for pileup_sites has been removed. This
option is not needed at the user level and is planned to be replaced
by the regional indexing used in pileup_cells().

- Added FilterParam option to trim reads based on fractional distance
from 5' (ftrim_5p) or 3' end (ftrim_3p).

- Incorporated RBPZ and VDB statistics from bcftools, now returned as
rowData columns when calling pileup_sites.

- A RangedSummarizedExperiment object is now directly returned from
pileup_sites. Using merge_pileups is no longer necessary and is not
an exported function.

- Renamed get_pileup to pileup_sites and create_se to merge_pileups

- Rename remove_clustered_variants, remove_multiallelic, and
remove_splice_variants to filter_* for consistency.

- Rewrote and renamed the single cell editing function sc_editing to
pileup_cells(). pileup_cells() does not require sorting and index by
cell barcode, uses a new format to specify sites to query and
requires providing the reference and alternate alleles of interest,
writes to disk in a sparse matrix compatible format to reduce memory
usage, and should have more performance as there is no need to query
a fasta index.

- Implemented method to collapse reads with duplicate UMIs.

- Added option to filter sites in pileup based on number of reads
containing a variant (#54)

- Added a NEWS.md file to track changes to the package.

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.26.0                        

New features

- Added a contributed vignette for ASCAT workflows (@Lydia-King,
#28).

[RAIDS](/packages/RAIDS)
-----

                       Changes in version 0.99.15                       

SIGNIFICANT USER-VISIBLE CHANGES

- Updating installation section in vignette.

                       Changes in version 0.99.14                       

SIGNIFICANT USER-VISIBLE CHANGES

- Adding missing author David Tuveson.
- Updating BiocViews terms.

                       Changes in version 0.99.13                       

SIGNIFICANT USER-VISIBLE CHANGES

- Update in Reference GDS vignette.

                       Changes in version 0.99.12                       

SIGNIFICANT USER-VISIBLE CHANGES

- Seven new loadable objects are available in the package.
- The new readSNVVCF() function enable the use of VCF SNP files as input for the runExomeAncestry() and runRNAAncestry() functions.

                       Changes in version 0.99.11                       

SIGNIFICANT USER-VISIBLE CHANGES

- Update main vignette.

                       Changes in version 0.99.10                       

SIGNIFICANT USER-VISIBLE CHANGES

- Update main vignette.

                       Changes in version 0.99.9                        

SIGNIFICANT USER-VISIBLE CHANGES

- Better documentation for the runRNAAncestry() function.

                       Changes in version 0.99.8                        

SIGNIFICANT USER-VISIBLE CHANGES

- Some examples have been updated in the documentation.

                       Changes in version 0.99.7                        

SIGNIFICANT USER-VISIBLE CHANGES

- The new runRNAAncestry() function executes most steps leading to the ancestry inference call on a specific RNA profile.
- A vignette describing the content of the Reference GDS files has been created.
- More parameter names have been changed to follow the camelCase style.

                       Changes in version 0.99.6                        

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette now referes to the generic formatted reference GDS
rather than 1KG GDS file to showcase that the software is not
dependant of the 1KG GDS file. Any refence dataset can be used as
long as the dataset is formatted into a GDS file.

                       Changes in version 0.99.5                        

SIGNIFICANT USER-VISIBLE CHANGES

- The function documentation has been improved.
- New vignette has been created. The vignette covers the steps done
by the runExomeAncestry() function.

                       Changes in version 0.99.4                        

SIGNIFICANT USER-VISIBLE CHANGES

- More parameter names have been changed to follow the camelCase
style.
- The function documentation has been improved.
- A wrapper function runExomeAncestry() is now available.

                       Changes in version 0.99.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- More parameter names have been changed to follow the camelCase
style.

BUG FIXES

- The warning related to the package man page has been removed.

                       Changes in version 0.99.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Most parameter names have been changed to follow the camelCase
style.

                       Changes in version 0.99.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- The new runExomeAncestry() function encapsulates multiple ancestry
inference steps in one command.

BUG FIXES

- Ensure GDS file is closed before using stop() in the
addPhase1KG2SampleGDSFromFile() function.

                       Changes in version 0.99.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The number of visible functions has been limited to simplify usage.


[Rarr](/packages/Rarr)
----

                         Changes in version 1.1                         

- Fixed bug when reading an array if the fill value in .zarray was
null.
- Addressed bug in makevars where Rarr.so could be compiled before
libblosc.a was ready. Also backported to Rarr 1.0.2. (Thanks to
Michael Sumner for reporting this issue:
https://github.com/grimbough/Rarr/issues/5)
- Corrected issue where fixed length string datatypes would be
written
with null terminators, resulting in strings that were one byte
longer than the dtype value written in the .zarray metadata. Also
backported to Rarr 1.0.3.
- Added support for reading and writing the fixed length Unicode
datatype, and for reading variable length UTF-8 datatype.

[rawrr](/packages/rawrr)
-----

                 Changes in version 1.9.2 (2023-10-24)                  

- Download RawFileReader DLL’s from
  <https://github.com/thermofisherlsms/RawFileReader> #66.

[Rcpi](/packages/Rcpi)
----

                 Changes in version 1.37.1 (2023-06-27)                 

Improvements

- Migrate rcdk to a runtime soft dependency using rlang, to avoid
build time issues.
- Migrate unit tests from using RUnit to testthat.
- Use GitHub Actions workflows for R CMD check and building the
pkgdown website.
- Fix URLs that are broken or moved in the documentation.

[recount3](/packages/recount3)
--------

                       Changes in version 1.11.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Windows users can now use http://duffel.rail.bio/recount3 again!
The
switch to httr::http_error() resolved the duffel access problem for
Windows users.

                       Changes in version 1.11.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Switched from using RCurl::url.exists() to !httr::http_error().

[RedeR](/packages/RedeR)
-----

                        Changes in version 2.6.0                        

- Major code refactoring.

- Improved readability.

- Reduced code complexities.

- Improved compatibility with standard igraph attributes.

[regionalpcs](/packages/regionalpcs)
-----------

                 Changes in version 0.99.0 (2023-08-25)                 

NEW FEATURES

- Added compute_regional_pcs() function to compute principal
components for given gene regions.

- Included gene region annotations dataset gene_annots.

IMPROVEMENTS

- Improved performance of compute_regional_pcs() by using a more
efficient algorithm.

- Updated vignette with more in-depth explanations and examples.

BUG FIXES

- Fixed bug in compute_regional_pcs() that produced incorrect results
for certain edge cases.

- Corrected typos in man pages.

DEPRECATED AND DEFUNCT

- Deprecated old_function_name() in favor of new_function_name().

DOCUMENTATION

- Added a detailed vignette for working with compute_regional_pcs().

- Improved man pages with runnable examples for all exported
functions.


[RegionalST](/packages/RegionalST)
----------

                       Changes in version 0.99.8                        

- Add unit tests

- Add input checking functions

                       Changes in version 0.99.2                        

- 
  New package RegionalST, for regional analysis of spatial
  transcriptomics data.

[rgoslin](/packages/rgoslin)
-------

                        Changes in version 1.5.0                        

Please note that this Bioconductor version is based on Goslin
version 2.0.0. See the Goslin repository and Goslin C++ repository
for
more details.

BioConductor 3.18 - Changes in 1.5.0

Improvements

- Better handling of mediators
- Translating gangliosides into new nomenclature structure
- Updated HMDB grammar for parsing mediators in lipid names: e.g.,
PA(P-16:0/LTE4)
- Parsing of adducts with heavy labeled isotopes now possible

Bug Fixes

- Minor bug fixes
- Removed wrong Lyso classification for SPB, SPBP, LHexCer, LSM

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.46.0                        

CHANGES

- Added support for reading nullable booleans and integers from
  the AnnData specification.  h5read() will detect these
  automatically an attempt to cooerce them to the appropriate R
  data format.

[rhdf5filters](/packages/rhdf5filters)
------------

                       Changes in version 1.14.0                        

CHANGES

- The package will now test for system libraries for several
compression tools and will use those libraries rather than compiling
from source if they are found.

[Rhdf5lib](/packages/Rhdf5lib)
--------

                        Changes in version 1.24                         

Bug fixes

- Fixed compilation issue when there was a space in the library path.
  (thanks to Pius Martinn @albert180 for reporting this,
  https://github.com/grimbough/rhdf5/issues/128)

- The configure script now detect the AR and RANLIB associated with R.
  This fixes an issue installing the package from source when using R
  provided by conda.
  (thanks to Jim Jeffers @jimjeffers for reporting this,
  https://github.com/grimbough/Rhdf5lib/issues/52)

[ribosomeProfilingQC](/packages/ribosomeProfilingQC)
-------------------

                       Changes in version 1.13.2                        

- fix the bug the estimatePsite did not pass ignore.seqlevelsStyle to
assignReadingFrame function.

                       Changes in version 1.13.1                        

- add ignore.seqlevelsStyle parameter to provide the flexibility of
seqlevels.

[RNAseqCovarImpute](/packages/RNAseqCovarImpute)
-----------------

                Changes in version 0.99.12 (2023-10-13)                 

- limmavoom_imputed_data_list now returns results for all variables in
  model formula

                Changes in version 0.99.11 (2023-10-06)                 

- License to GPL-3

                Changes in version 0.99.10 (2023-10-06)                 

- No longer exports internal functions

- License to GPL-3 + file LICENSE

- Updated vignette introduction

                 Changes in version 0.99.8 (2023-09-18)                 

- Replaced RNAseqCovarImpute_data.RData with two individual rda files
  for example_DGE and example_data

- Added documentation for the package (?RNAseqCovarImpute)

- Added unit test structure and some basic unit tests

                 Changes in version 0.99.7 (2023-09-16)                 

- Parallel now implemented with BiocParallel

- Added validity tests

                 Changes in version 0.99.0 (2023-04-23)                 

- Submitted to Bioconductor

[RnBeads](/packages/RnBeads)
-------

                       Changes in version 2.19.0                        

- Fixed issue with logical having a length greater than 1

[rols](/packages/rols)
----

                        Changes in version 2.29                         

CHANGES IN VERSION 2.29.1

- Fix unit tests.

CHANGES IN VERSION 2.29.0

- New devel version (Bioc 3.18)

[rpx](/packages/rpx)
---

                         Changes in version 2.9                         

rpx 2.9.1

- Remove the generated subdir in the ftp_url when creating the
PXDataset object <2023-09-26 Tue> (see issue #25).

[Rsubread](/packages/Rsubread)
--------

                       Changes in version 2.16.0                        

- 
  Add support for ARM64 platforms.

- 
  Add support for CBCL format in cellCounts.

[RTCGAToolbox](/packages/RTCGAToolbox)
------------

                       Changes in version 2.32.0                        

New features

- The functions getCNGECorrelation, getDiffExpressedGenes, and
getSurvival have been removed from the package.

[Rvisdiff](/packages/Rvisdiff)
--------

                 Changes in version 0.99.4 (2023-09-27)                 

- DESeq2 and limma packages changed from 'Imports' to 'Suggests'

                 Changes in version 0.99.3 (2023-09-26)                 

- More checks on the input data

- More arguments to provide the names of columns

- Paths constructed using file.path() to be OS-agnostic

- Added package help-(man-)page (to be viewed via ?Rvisdiff)
  summarizing the
  packages functionality

- In the vignette, added BiocStyle functions to hyperlink external
  packages

- In the vignette, added code-style in order to distinguish/highlight
  anything
  R-related (package/function/argument/variable name), including
  Rvisdiff

- implemented more comprehensive unit testing

                 Changes in version 0.99.2 (2023-09-21)                 

- Maintainer subscribed to the Bioc-Devel mailing list

                 Changes in version 0.99.1 (2023-09-21)                 

- R version changed

- Added readme

                 Changes in version 0.99.0 (2023-08-04)                 

- Submitted to Bioconductor

[rWikiPathways](/packages/rWikiPathways)
-------------

                       Changes in version 1.22.0                        

- Reimplementation of every function, replacing web service with static
  JSON

- getRecentChanges now returns data.frame of pathways per last-edit
  date

- findPathwaysByText has new param, "field" to optionally specify which
  fields to search

- getPathwayHistory now opens pathway commit history in browser

- getPathwayInfo now returns all pathways if param is left NULL

- New functions:
  - listCommunities
  - getPathwaysByCommunity
  - getPathwayXXXsByCommunity
  - getCounts
  - findPathwaysByOrcid
  - getCurationStatus

- Deprecated functions:
  - getCurationTags
  - getCurationTagNames
  - getXXXByCurationTag
  - wikipathwaysAPI
  - wikipathwaysGET

[S4Arrays](/packages/S4Arrays)
--------

                        Changes in version 1.2.0                        

NEW FEATURES

- Add abind() generic + default method.

- Add drop() method and dim() setter for Array objects.

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.40.0                        

NEW FEATURES

- Subscript now can be any 1D array-like object when subsetting a
  Vector derivative.

SIGNIFICANT USER-VISIBLE CHANGES

- Drop empty "DataFrame and DataFrameList objects" section from
  S4VectorsOverview.Rnw vignette.

BUG FIXES

- Fix incorrect C-level if statement in Rle_utils.c

[SAIGEgds](/packages/SAIGEgds)
--------

                        Changes in version 2.2.0                        

- fix "Matrix-deprecated" when calling `as(<dsCMatrix>, "dgCMatrix")`


[SCArray](/packages/SCArray)
-------

                       Changes in version 1.10.0                        

- update for DelayedArray (>= v0.27.2)

                        Changes in version 1.8.4                        

- `rowMeans()`, `rowSums()`, `colMeans()`, `colSums()` with row or
  column
  names

- override S4 functions rowDiffs(), colDiffs(), rowSdDiffs(),
  colSdDiffs(),
  rowVarDiffs(), colVarDiffs(), rowLogSumExps(), colLogSumExps()

                        Changes in version 1.8.3                        

- progress bar in `scHDF2GDS()`

- update `%*%` with SC_GDSMatrix

- update `rbind()` and `cbind()` with SC_GDSMatrix

- `runPCA()` on SC_GDSMatrix

[SCArray.sat](/packages/SCArray.sat)
-----------

                        Changes in version 1.0.3                        

- update `FoldChange.SCArrayAssay()` and `FindMarkers.SCArrayAssay()`

- fix `subset.SCArrayAssay()`

- update `RunICA.SCArrayAssay()`, `RunSPCA.SCArrayAssay()`,
  `RunLDA.SCArrayAssay()`, `RenameCells.SCArrayAssay()`,
  `merge.SCArrayAssay()`

[scBubbletree](/packages/scBubbletree)
------------

                        Changes in version 1.3.1                        

- vignette checks, minor modification

[scDesign3](/packages/scDesign3)
---------

                 Changes in version 1.0.1 (2023-01-20)                  

- Submitted to Bioconductor

                 Changes in version 0.99.5 (2023-07-15)                 

- Update the MVN sampling

[scider](/packages/scider)
------

                       Changes in version 0.99.0                        

First release.

[scp](/packages/scp)
---

                        Changes in version 1.11                         

scp 1.11.3

- feat: added readSCPfromDIANN() that creates a QFeatures object from
DIANN output tables.

scp 1.11.2

- Nothing yet.

scp 1.11.2

- feat: added reportMissingValues(), jaccardIndex(),
cumulativeSensitivityCurve() and predictSensitivity() to facilitate
reporting missing values. The vignette is also adapted with the new
functionality.
- docs: created vignette about reporting missing values in SCP
- fix failing unit test.

scp 1.11.1

- Updated citation

scp 1.11.0

- New Bioconductor 3.18 (devel) release

[scPCA](/packages/scPCA)
-----

                 Changes in version 1.15.1 (2023-06-19)                 

- Address useNames issue in colSds() that caused tests to throw
warnings.

[screenCounter](/packages/screenCounter)
-------------

                        Changes in version 1.2.0                        

- Added support for counting dual barcodes via
  countDualBarcodes().

- Refactored countSingleBarcodes() to support arbitrary numbers
  of substitutions, insertions and deletions.

- Simplified countComboBarcodes() at the expense of dropping
  support for edits inside variable regions.

[scRNAseqApp](/packages/scRNAseqApp)
-----------

                       Changes in version 1.1.10                        

- Fix the download button for the explorer module.

                        Changes in version 1.1.9                        

- Fix the repeat retreive for reference by createAppConfig.

                        Changes in version 1.1.8                        

- Change the filepath for credential file.

                        Changes in version 1.1.7                        

- Fix the bug in scInit for datafolder.

                        Changes in version 1.1.6                        

- Fix the download button which download the fixed values.

- Add multiple subset buttons.

                        Changes in version 1.1.5                        

- List the available datasets by users privilege.

                        Changes in version 1.1.4                        

- Fix the bug if the data is removed where user is explorering.

                        Changes in version 1.1.3                        

- Fix the bug introduced by changing actionBution to checkboxInput for
  modBubbleHeatmap.

                        Changes in version 1.1.2                        

- Add multiple layer of violin plot

- Add sample order controls for stats plot

                        Changes in version 1.1.1                        

- Fix the bug that `coverage` function is not imported.

[SDAMS](/packages/SDAMS)
-----

                       Changes in version 1.21.1                        

- update maintainer email address

[sechm](/packages/sechm)
-----

                 Changes in version 1.9.4 (2023-07-14)                  

- row ordering now changed in the heatmap rather than matrix

- fixed a few bugs, added a few extras

- enabled list of features as input

[seqArchRplus](/packages/seqArchRplus)
------------

                        Changes in version 1.1.3                        

New

- (User-facing) Enables choosing the out device type for motif
heatmaps

                        Changes in version 1.1.2                        

New

- (User-facing) Ability to generate HTML reports to view large
combined panels for multiple processed samples as scrollable
carousels

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.42.0                        

UTILITIES

- new option 'write.rsid' in `seqGDS2BED()`

                       Changes in version 1.40.1                        

BUG FIXES

- `seqAddValue(gdsfile, varnm="position")` works correctly

[shinyMethyl](/packages/shinyMethyl)
-----------

                       Changes in version 1.37.2                        

- Removing failing test

[signeR](/packages/signeR)
------

                        Changes in version 2.3.3                        

- Minor documentation fixes

                        Changes in version 2.3.2                        

- New function: genCountMatrixFromMAF, to generate a count matrix from
  a MAF file

- support for reading MAF files on signeRFlow

- support for utilizing any BSgenome on signeRFlow

                        Changes in version 2.3.1                        

- Moved the source repository to TojalLab

- Fixed issue when running with fixed signatures matrix

- Reduced the maxeval time in rare cases when there is no convergence

[signifinder](/packages/signifinder)
-----------

                        Changes in version 1.4.0                        

- Add glioCellStateSign function to compute the glioblastoma
  cellular states defined by Neftel C. et al. Cell (2019).

- Add whichAssay argument to signature functions to allow the
  user to specify which assay to use for the signature
  computation.

- Users can now plot also other signatures not computed with
  signifinder when using the heatmapSignPlot, correlationSignPlot
  and ridgelineSignPlot functions to compare them with the
  signatures computed with signifinder.

                        Changes in version 1.2.1                        

- Add evaluationSignPlot function to show some technical
  information of the signatures computed.

- Add nametype argument to geneHeatmapSignPlot function to allow
  more gene name ID in data.

- The vignette now contains an example with a single-cell dataset
  and an example with a spatial transcriptomics dataset.

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.10.1 (2023-07-26)                 

- Added function for bubble plot
- In SCTK-QC pipeline, added support for batch processing multiple
inputs
- In SCTK-QC pipeline, added support for importing and exporting
AnnData objects
- In SCTK-QC pipeline, fixed a bug causing YAML output files to be
empty
- Update the SCTK-QC tutorial
- Fixed bug in combineSCE causing it to create multiple copies of row
or column data

[SparseArray](/packages/SparseArray)
-----------

                        Changes in version 1.2.0                        

NEW FEATURES

- Add aperm() method for SVT_SparseArray objects.

- Add abind() method for SparseArray objects.

- Add dim() setter for SVT_SparseArray objects.

- Introduce nzwhich() generic and method for SparseArray derivatives.
  Also provide a default method for ordinary arrays and other
  array-like
  objects.

- Implement 'Logic' ops on SVT_SparseArray objects.

- Implement 'Math'/'Math2' ops on SVT_SparseArray objects of type
  "double".

- 'Compare' ops now support SVT_SparseArray objects of type() "raw"
  or "complex".

- All matrixStats methods (except row/colMedians()) now work on
  multidimensional SVT_SparseArray objects and support the 'dims'
  argument, like the row/colSums() and row/colMeans() functions in base
  R.

- Add row/colAnys() + row/colAlls() + row/colAnyNAs() + row/colProds()
  methods for SVT_SparseArray objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Rename nzvals() slot getter (for COO_SparseArray objects) ->
  nzdata().
  Also reintroduce nzvals() as a fast way to get 'x[nzwhich(x)]' on a
  sparse array-like object 'x'.

- Re-implement all matrixStats methods (except row/colMedians()) for
  SVT_SparseArray objects in C.

[sparseMatrixStats](/packages/sparseMatrixStats)
-----------------

                        Changes in version 1.13                         

- Make sparseMatrixStats compatible with matrixStats release v1.0.0.
  In particular change 'useNames'to 'TRUE' by default.

- Add fast path for 'rowSums2(x, cols = logical_vector)'

- Add useNames parameter to all functions

- fix incomplete method signature of rowQuantiles

[spaSim](/packages/spaSim)
------

                        Changes in version 1.3.1                        

- Added newly published paper.

                        Changes in version 1.3.0                        

Development version in Bioconductor 3.18

[SpatialExperiment](/packages/SpatialExperiment)
-----------------

                 Changes in version 1.11.2 (2023-09-01)                 

- move DropletUtils package to Suggests

                 Changes in version 1.11.1 (2023-08-21)                 

- add methods to rotate/mirror spatial coordinates and objects

[SpatialFeatureExperiment](/packages/SpatialFeatureExperiment)
------------------------

                        Changes in version 1.2.1                        

- Fixed bug in .check_features and .symbol2id where "symbol" column
is
hard coded

[Spectra](/packages/Spectra)
-------

                        Changes in version 1.11                         

Changes in 1.11.11

- Fix issue with filterFourierTransformArtefacts function (see issue
#302). Thanks Adriano Rutz for reporting.

Changes in 1.11.10

- peaksData,MsBackendMemory returns a data.frame if additional peak
variables (in addition to "mz" and "intensity") are requested. For
columns = c("mz", "intensity") (the default) a list of matrix is
returned.
- peaksData,Spectra returns either a matrix or data.frame and ensures
the peak data is correctly subset based on the lazy evaluation
processing queue.
- $,Spectra to access peak variables ensures the lazy evaluation
queue
is applied prior to extracting the values.
- applyProcessing correctly subsets and processes all peak variables
depending on the processing queue.
- spectraData<-,Spectra throws an error if processing queue is not
empty and values for peaks variables should be replaced.
- $<-,Spectra throws an error if processing queue is not empty and a
peaks variable is going to be replaced.
- Add full support for additional peaks variables to
MsBackendDataFrame.

Changes in 1.11.9

- Add filterPrecursorPeaks to allow filtering peaks within each
spectrum with m/z values relative to the precursor m/z of the
spectrum.

Changes in 1.11.8

- Add an example to the vignette describing how spectral similarity
scores from the msentropy package can be used with compareSpectra.

Changes in 1.11.7

- Fix in compareSpectra to also pass parameters ppm and tolerance to
the peak similarity calculation functions FUN: this allows to use
custom similarity function with integrated mapping of peaks.
- Add joinPeaksNone to skip the peak matching in compareSpectra if
the
similarity scoring function performs its own peak matching.
- Only use parallel processing in setBackend,Spectra if both backends
support it.

Changes in 1.11.6

- Add filterPrecursorMaxIntensity function.
- Add filterPrecursorIsotopes function.

Changes in 1.11.5

- Add scalePeaks function (see issue #291).

Changes in 1.11.4

- Import uniqueMsLevels from ProtGenerics.

Changes in 1.11.3

- Rename combinePeaks for lists of peak matrices into
combinePeaksData.
- Add combinePeaks generics.
- Add combinePeaks,Spectra to combine peaks within each spectrum in a
Spectra.

Changes in 1.11.2

- Add deisotopeSpectra and reduceSpectra functions.

Changes in 1.11.1

- Add example for filtering precursor m/z peaks from fragment spectra
to the vignette.

[SpectralTAD](/packages/SpectralTAD)
-----------

                 Changes in version 1.16.1 (2023-07-04)                 

- Adjust installation url

[SPIAT](/packages/SPIAT)
-----

                        Changes in version 1.3.5                        

SIGNIFICANT USER-VISIBLE CHANGES

- Added a new argument margin_dist for define_structure(). Specifying
the margin width with microns instead of layers of cells.

                        Changes in version 1.3.4                        

BUG FIXES

- Added legends to plot_cell_categories() plots when layered = TRUE.
- Added a distance parameter to fix the dimension error in
calculate_spatial_autocorrelation().

                        Changes in version 1.3.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- Fixed thresholding bug in predict_phenotypes().

NOTES

- Added Shiny App (reading data) link.

                        Changes in version 1.3.2                        

NOTES

- Fixed typo in citation and SPIAT overview diagram.
- Moved the following packages from Imports to Suggests: alphahull,
plotly.

                        Changes in version 1.3.1                        

- Added citation to the newly published paper.

                        Changes in version 1.3.0                        

Development version on Bioconductor 3.18.

[splatter](/packages/splatter)
--------

                 Changes in version 1.26.0 (2023-10-25)                 

- 
  Fixed a bug in splatSimPathDE() where DE factors were not
  adjusted based on the path origin (path.from parameter). This
  affected paths where the path origin was not the simulation
  origin (i.e.  path.from != 0), particularly when the path DE
  was minimal.  With this fix paths should no longer drift
  towards the origin.

[SpliceWiz](/packages/SpliceWiz)
---------

                 Changes in version 1.3.2 (2023-06-05)                  

- Bugfix: error when running featureCounts wrapper due to non-numeric
  assignment of single / paired end reads

                 Changes in version 1.3.1 (2023-05-20)                  

- Bugfix: fixed - static plot coverages did not show when
  reverseGenomeCoords = TRUE

- Bugfix: error when running featureCounts with overwrite = TRUE

[sSNAPPY](/packages/sSNAPPY)
-------

                 Changes in version 1.4.4 (2023-10-09)                  

- Remove databases that are no longer supported: pathbank, panther,
pharmgkb, smpdb

                 Changes in version 1.4.3 (2023-07-31)                  

- Updated the perturbation scoring step to account for the
orientiaton
of topology matrices for KEGG pathways (row - Downstream genes;
column - Upstream genes)
- Add prefix parameter to the weight_ss_fc function to allow
user-specified prefix that are other than "ENTREZID:"

                 Changes in version 1.4.2 (2023-07-12)                  

- Updated the permtuation strategy so that the
generate_permuted_scores function now construct all possible
permuted pairs by default
- Updated the plot_community function so KEGG pathways that are not
assigned to categories will be ignored in community labeling

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

                       Changes in version 1.13.5                        

- improve consistency between methods in fold change computations

- better handling of NA in fold change computations

                       Changes in version 1.13.4                        

- fix PLSDA predicted group assignment

- add option to PLSDA to use probability for yhat for predictions

- update tests

                       Changes in version 1.13.3                        

- d-ratio equations changed to match Broadhurst et al (2018)

- add tests for d-ratio

- fix PCA eigenvalues calculation

                       Changes in version 1.13.2                        

- d-ratio filter now correctly removes features

- update descriptions/documentation for missing value filters

- add tests for d-ratio

                       Changes in version 1.13.1                        

- hotfix vector_norm now correctly normalises samples to length 1

- re-sync with Bioc devel

- fix broken documentation

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.32.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette gains a new section about interactive visualization with
  iSEE.

[SVMDO](/packages/SVMDO)
-----

                 Changes in version 1.1.1 (2023-06-28)                  

- Removal of unnecessary variables in differential expression and
classification analyses


[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.13.8                        

- Updates to all EvoWeaver documentation files
- Fixed small bug in PhyloDistance causing Method='JRF' to return
similarity rather than the distance
- Fixed small bug in TreeDistance.EvoWeaver resulting in an
inconsistent calculation of score when using TreeMethods='JRF'

                       Changes in version 1.13.7                        

- Small fixes

                       Changes in version 1.13.6                        

- ProtWeaver and ProtWeb have been renamed to EvoWeaver and EvoWeb,
respectively
- New sequence level method for EvoWeaver
- Various small internal updates to EvoWeaver

                       Changes in version 1.13.5                        

- Minor changes to SelectByK and vignette

                       Changes in version 1.13.4                        

- New predictor PAPV.ProtWeaver to calculate p-values for
presence/absence profiles.
- ContextTree now uses MirrorTree with species tree correction and
p/a
overlap correction
- Updates to documentation

                       Changes in version 1.13.3                        

- predict.ProtWeaver now supports multiple algorithms at once (ex.
predict(ew, Method=c("Jaccard", "Hamming")))
- Documentation for ProtWeaver and associated methods has been
updated
to match recent updates.

                       Changes in version 1.13.2                        

- FastQFromSRR function added as a convenience wrapper for the
SRAtoolkit function fastq-dump.

                       Changes in version 1.13.1                        

- SuperTree now works directly with dist objects, providing better
performance and scaling
- Updates to simMat objects
- No longer throw a warning when initialized in RStudio
- Formatting is cleaner and supports larger object names
- Updates to NVDC.ProtWeaver
- Now supports amino acid sequences using the DNAseqs=FALSE
argument
- Now calculates a p-value-weighted score
- Adds MakeBlastDb function to create a BLAST database from R, plus
associated documentation updates
- Smaller fixes to some ProtWeaver methods
- predict.ProtWeaver no longer returns using invisible (this was
annoying and unneccessary)
- APC correction for MutualInformation.ProtWeaver removed to allow
for parallelization
- MirrorTree.ProtWeaver now works correctly with
MTCorrection="speciestree"
- CorrGL.ProtWeaver now uses Fisher's Exact Test for p-values
rather than the R value of spearman correlation
- Many internal performance improvements
- ProtWeaver almost entirely uses dist objects rather than matrix,
saving significantly on memory
- faster Cophenetic function implemented internally
- Copied internal .Call('cophenetic') from DECIPHER to SynExtend
to avoid potential namespace issues
- Small fixes to remove some notes from BiocCheck::BiocCheck()
- Variety of small updates to pass BiocCheck

[syntenet](/packages/syntenet)
--------

                        Changes in version 1.3.5                        

NEW FEATURES

- Added function run_last() to run alternative BLAST search. Vignette
was updated accordingly.

                        Changes in version 1.3.4                        

BUG FIXES

                        Changes in version 1.3.3                        

                        Changes in version 1.3.2                        

BUG FIXES

- Strand information is now preserved in the output of
process_input()
(for users who want to plot synteny).

                        Changes in version 1.3.1                        

NEW FEATURES

- Added function collapse_protein_ids() to replace protein IDs in
sequence names (equivalent to FASTA headers) with gene IDs. If there
are multiple protein for the same gene, onlt the longest is kept.
Vignette was updated accordingly.

[tadar](/packages/tadar)
-----

                       Changes in version 0.99.0                        

- Submitted to Bioconductor.

[TargetSearch](/packages/TargetSearch)
------------

                        Changes in version 2.4.0                        

NEW FEATURES

- `plotPeakRI`: new option to plot by RT. This, however, needs the
  parameter `dev` because the relationship between RI and RT is
  unknown.

- `FindPeaksAll`: new option to search by RT. Basically this change
  implements the functionality of `plotPeakRI`.

- C code: refactor or most part of the code that deals with file
  parsing,
  in particular the text parser. The parser is now based on regular
  expressions, which should be more robust in detecting errors.

- Tests: though is not relevant to the user, TargetSearch include more
  unit testing for several internal functions.

BUG FIXES

- Add assertion that the option `TS_RI_columns` has three elements and
  that it is a character or integer vector.

[TBSignatureProfiler](/packages/TBSignatureProfiler)
-------------------

                       Changes in version 1.14.0                        

Bug Fixes

- Fixed bug for signatureBoxplot by removing quotations around
variables. This issue was introduced with the newest ggplot2
version. Thanks to Arthur VanValkenburg for identifying the issue
and solution.
- Removed Zimmer_RES_3 which was an identical signature as the
previously published Sweeney_OD_3.
- Fixed error in example for compare_algs caused by the signature
being used.

Minor Changes

- Added Vargas_18 and Vargas_45 signatures
(doi: 10.1371/journal.pcbi.1010770)

[TCGAbiolinks](/packages/TCGAbiolinks)
------------

                       Changes in version 2.29.1                        

- Removing support to legacy archive since it will be shutdown by GDC
  soon.

- When saving files we will not include folders prefix
  legacy/harmonized anymore

[TCGAutils](/packages/TCGAutils)
---------

                       Changes in version 1.22.0                        

Bug fixes and minor improvements

- UUIDtoBarcode returns barcodes consistent with Genomic Data Commons
API update

[TENxIO](/packages/TENxIO)
------

                        Changes in version 1.4.0                        

Bug fixes and minor improvements

- Skip unit tests when remote H5 access is not configured.

[terraTCGAdata](/packages/terraTCGAdata)
-------------

                        Changes in version 1.6.0                        

New features

- findTCGAworkspaces has been removed from the package.

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

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.37.15                       

- Fix reading valid pairs typo.

                       Changes in version 1.37.14                       

- update the lollipopPlot documentation to fix the issue
  unable to find an inherited method for function 'content' for
  signature 'response'.

                       Changes in version 1.37.13                       

- fix the feature legend space.

                       Changes in version 1.37.12                       

- add ARA function.

- fix the missleading color legend of interaction data.

- fix the bug in addInteractionAnnotation.

                       Changes in version 1.37.11                       

- fix the minimal font of optimized style.

                       Changes in version 1.37.10                       

- add function GIoperator.

                       Changes in version 1.37.9                        

- add function addInteractionAnnotation for Interaction data.

                       Changes in version 1.37.8                        

- change the filter condiction for tads annotation of Interaction data.

                       Changes in version 1.37.7                        

- add tads annotation for Interaction data.

                       Changes in version 1.37.6                        

- add lollipop_style_switch_limit for lollipop plot.

                       Changes in version 1.37.5                        

- add border_color for Interaction data.

                       Changes in version 1.37.4                        

- update the documentation for lollipopPlot for changes of the snp
  label.

                       Changes in version 1.37.3                        

- use strawr to replace the local C++ script.

                       Changes in version 1.37.2                        

- remove the CXX_STD = CXX11 from straw Makevars.

                       Changes in version 1.37.1                        

- add stop message when there is negative values for grid.pie plot.

[transomics2cytoscape](/packages/transomics2cytoscape)
--------------------

                       Changes in version 1.11.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Added pbapply progress bar at the "getEdgeInfo" points.

DATA or DOCUMENT CHANGES

- Wrote the time it takes to complete the process in my environment for
  each use case.

[treeio](/packages/treeio)
------

                       Changes in version 1.25.4                        

- reexport as.phylo.hclust_node() to fix as.phylo.pvclust() issue
(2023-8-25, Fri, #110)

                       Changes in version 1.25.3                        

- add find.hclust.igraph() method to hierarchical clustering graph
nodes (2023-08-11, Fri, #105)
- update spt() and as.phylo.igraph() to consider edge attributes
(#105)
- move tree operation methods to the 'tidytree' package so that this
package is focus on input, output and object conversion
- fixed issue in parse BEAST file that contains negative branch
length
(2023-08-03, Thu, #106)

                       Changes in version 1.25.2                        

- spt method to find shortest path tree (2023-07-14, Fri, #102)
- update old-style 'CITATION' from citEntry() to bibentry()
(2023-07-14, Fri, #102)
- bug fixed in as.treedata() for tbl_df object (2023-07-14, Fri,
#101)
- keep.tip() method to remove all tips excepts the selected tips
(2023-07-13, Thu, #100)
- better support of converting 'igraph' object to 'phylo' object
(2023-07-12, Wed, #99)

                       Changes in version 1.25.1                        

- bug fixed in read.nhx() when metadata contains both character and
numeric (e.g., AAA111) (2023-05-31, Wed, #97)

[TSAR](/packages/TSAR)
----

                       Changes in version 0.99.0                        

- First Draft for BioConductor review, completed on August 15, 2023

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.20.0                        

- Add argument to summarizeToGene(): assignRanges that takes either
"range" (default) or "abundant", and determines the ranges that are
attached to the SE (rowRanges). Note that this new argument does not
affect the data aggregation at all (counts and abundance are
summarized to gene by tximport). The default behavior of
summarizeToGene() returns ranges that correspond to the range of the
isoforms of the gene, that is the leftmost basepair to the rightmost
basepair of any isoform. The non-default "abundant" instead returns
the range of the most abundant isoform in the data, averaging over
samples. Information about the choice of range is included in mcols
- Added support for piscem-infer: type="piscem" also auto-detected
from file ending.

                       Changes in version 1.19.8                        

- Added support for piscem-infer: type="piscem" also auto-detected
from file ending.

                       Changes in version 1.19.6                        

- Fixed genome build for mouse M26 and higher to GRCm39, thanks to
Charlotte Soneson.

[tximport](/packages/tximport)
--------

                       Changes in version 1.30.0                        

- Support for piscem-infer: use `type="piscem"`.

[UniProt.ws](/packages/UniProt.ws)
----------

                       Changes in version 2.42.0                        

NEW FEATURES

- `pageSize` and `n` arguments added to `queryUniProt` to expose
  underlying API request defaults. It is recommended to set the
  pageSize
  to a large value e.g. 500 for large queries.

BUG FIXES AND MINOR IMPROVEMENTS

- Support for directly mapping to 'Ensembl' IDs with `select`.

- Fixed issue with pagination with large queries (over 25 results) in
  `queryUniProt` (@jdreyf, #23)

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.18.1                        

BUG FIXES

- Fixed compilation flags causing errors on linux.

[updateObject](/packages/updateObject)
------------

                        Changes in version 1.6.0                        

BUG FIXES

- Minor tweak to low-level utilities update_rds_file()
  and update_rda_file().

[variancePartition](/packages/variancePartition)
-----------------

                        Changes in version 2.0.5                        

- May 31, 2023
- fix convergence issues
- fix initialization of lmer() fit
- use 1 OMP thread internally, then restore to original value

                        Changes in version 2.0.4                        

- May 30, 2023
- When running dream(), ensure model convergence using second fitting
with Nelder_Mead to avoid edge cases where the approximate hessian
from lmerTest::as_lmerModLT() has a negative eigenvalue
- fix issue in get_prediction() returning NA values when variables
modeled as categorical and levels are omitted
- fix issue in voomWithDreamWeights() when some genes don't converge
- retry lmer() model fit with another optimizer after it fails
convergence test.

                        Changes in version 2.0.3                        

- May 13, 2023
- fix vcov()

                        Changes in version 2.0.2                        

- May 17, 2023
- add matrix argument to mvTest()

                        Changes in version 2.0.1                        

- May 12, 2023
- mvTest() now shrinks covariance using the Schafer-Strimmer method
- vcovSqrt() returns the matrix whose cross product gives the vcov()
result from fits with dream()

                        Changes in version 2.0.0                        

- April 20, 2023
- Major code refactoring to:
- improve code reuse
- simplify debugging and maintaining code
- simplify addition of new features
- improve error handling
- some linear mixed model analyses are 50% faster
- enable additional features for dreamlet package that depends
heavily on variancePartition.

                       Changes in version 1.31.22                       

- Oct 19, 2023
- fix handling of variables with missing data
- return fit$genes properly

                       Changes in version 1.31.21                       

- Oct 16, 2023
- handle weights properly when the linear mixed model fails for some
genes
- lmFit() and
- in iterRows() set scale = FALSE as default
- in voomWithDreamWeights(), scale in input weights and weights in
side fitVarPartModel()
- get weights estimated most similar to voomLmFit()
- in dream() use rescaleWeights = FALSE to get sigma estimates
compatable with lmFit()

                       Changes in version 1.31.20                       

- Sept 26, 2023
- allow weights to be a matrix in voomWithDreamWeights()

                       Changes in version 1.31.19                       

- Sept 22, 2023
- add rescaleWeightsAfter argument to voomWithDreamWeights()

                       Changes in version 1.31.18                       

- Sept 5, 2023
- improved error handling for fitVarPartModel(),
fitExtractVarPartModel(), and voomWithDreamWeights()

                       Changes in version 1.31.16                       

- August 18, 2023
- in dream(), if "Kenward-Roger" is specified but gives covariance
matrix that has poor condition number or is not positive definite,
then fall back to "Satterthwaite" for hypothesis testing in linear
mixed models
- Update documentation, and reformat code

                       Changes in version 1.31.15                       

- August 10, 2023
- fit = dream() now returns fit$loglik (the log-likelihood for each
gene), and fit$edf (the effective degreees of freedom for each gene)

                       Changes in version 1.31.13                       

- August 7, 2023
- fix bug in calcVarPart() where weights was ignored in some cases
- add additional tests to check this

                       Changes in version 1.31.12                       

- July 3, 2023
- makeContrastsDream() converts NA contrasts to NULL

                       Changes in version 1.31.11                       

- setting voomWithDreamWeights(..., span="auto") now estimates tuning
parameter from data using fANCOVA::loess.as()

                       Changes in version 1.31.10                       

- filterInputData() now ensures EList contains a matrix

                       Changes in version 1.31.9                        

- Fix issue in mvTest() when specifying features with strings

                       Changes in version 1.31.8                        

- Fix error message when linear mixed model fails

                       Changes in version 1.31.7                        

- update mvTest() to run in parallel

                       Changes in version 1.31.6                        

- update mvTest() to include Hotelling T2 test and LS.empirical()

                       Changes in version 1.31.1                        

- June 4, 2023
- Rename for Bioconductor compatability

[Voyager](/packages/Voyager)
-------

                        Changes in version 1.3.1                        

- Removed functions and arguments deprecated in 1.2.0

                        Changes in version 1.2.6                        

- Fixed bug in plotColGraph when one out of multiple samples is
plotted.
- Allow 16 bit images in spatial plotting functions.
- Removed adespatial from Suggests as it's only used as a reference
in
unit tests and it got removed from CRAN.

                        Changes in version 1.2.5                        

- Use imgRaster getter rather than the S4 no-no of @image to get
images to plot, as the latter will no longer work as of SFE 1.2.3
that wraps SpatRaster images when saving RDS. Reading RDS won't
unwrap so images need to be unwrapped when they're needed.

                        Changes in version 1.2.4                        

- Remove useNames = NA warning when calling MULTISPATI; the warning
comes from generic of colVars.
- Use algebraic eigenvalues for MULTISPATI when either nfposi or
nfnega is 0
- Added bins_contour argument to moranPlot to change the number of
bins in cell density contours

                        Changes in version 1.2.3                        

- Fix bug when plotting a feature with illegal name alongside another
feature with legal name
- Make sure runBivariate and calculateBivariate use gene symbols in
results even if Ensembl IDs are specified when swap_rownames is set
- Change secondary sequential palette in the light theme to YlOrRd so
it's more distinguishable from the Blues primary palette at low
values

                        Changes in version 1.2.2                        

- Some minor bugs: runBivariate gets correct feature names when only
feature1 is specified and swap_rownames is used to show gene symbol
- Correct output for cross variogram maps for only one pair of genes
- Added default_attr to localmoran_bv's SFEMethod
- Don't plot attribute when localResult is a vector and there's no
default attr
- When plotting multiple features, the panels follow the same order
the features are specified
- Allow illegal characters in names of colData and reducedDims in
plots
- Plot only one component in spatialReducedDim with the components
argument
- Deprecate plotColDataBin2D and plotRowDataBin2D

[weitrix](/packages/weitrix)
-------

                       Changes in version 1.13.1                        

- Use read.csv rather than read_csv in vignettes, as read_csv was
causing a hard-to-reproduce error when building vignettes.

[xcms](/packages/xcms)
----

                       Changes in version 3.99.6                        

- Add method to coerce a `XcmsExperiment` to a `xcmsSet` (issue #696).

- Support providing only `mz` or `rt` also for
  `chromatogram,MsExperiment`.

                       Changes in version 3.99.5                        

- Only `mz` or `rt` need to be provided for `chromatogram`.

                       Changes in version 3.99.4                        

- Add `chromPeakChromatograms` function to extract (EIC) chromatograms
  for
  chromatographic peaks.

                       Changes in version 3.99.3                        

- Small fixes in the *direct injection* vignette.

- Add parameter `isolationWindowTargetMz` to the `chromatogram`
  function for
  `MsExperiment` and `XcmsExperiment` to ensure MS2 chromatographic
  data is
  extracted from the MS2 spectra containing fragments of the compound
  of
  interest.

                       Changes in version 3.99.2                        

- Add the `xmse` data set representing an `XcmsExperiment` object.

- Update the *compounding* vignette to use the new objects.

- Add `loadXcmsData` to load test data objects (and fix/update paths).

- Add `groupFeatures` methods for `XcmsExperiment`.

- Fix issue in `featureArea` for `XcmsExperiment`.

- Update main vignette to use and describe the new data objects.

- Add `findChromPeaksIsolationWindow` method for `MsExperiment` and
  `XcmsExperiment`.

- Make `reconstructChromPeakSpectra` a method.

- Add `reconstructChromPeakSpectra` implementation for
  `XcmsExperiment`.

- Add `filterIsolationWindow` for `MsExperiment` and `XcmsExperiment`
  to filter
  spectra (and eventually chromatographic peaks) based on the isolation
  window.

- Update the LC-MS/MS vignette adding also an example how to deisotope
  SWATH
  MS2 spectra.

                       Changes in version 3.99.1                        

- `featureSummary` and `overlappingFeatures` gain support for
  `XcmsExperiment`.

- Fix in `featureChromatograms` to ensure a valid object is returned.

                       Changes in version 3.99.0                        

- Add `XcmsExperiment` and support for `MsExperiment`/`Spectra`: add
  all
  functionality for a full xcms processing on a `MsExperiment` object.

- Fix issue in `refineChromPeaks` with `MergeNeighboringPeaksParam`
  where a
  wrong apex position was considered in the evaluation whether
  candidate peaks
  should be merged (would only happen for merging of > 2 candidate
  peaks).

- Re-write the `reconstructChromPeakSpectra` for DIA data analysis to
  fix an
  issue with chromatographic peaks in overlapping SWATH isolation
  windows and
  generally to improve performance.

[YAPSA](/packages/YAPSA)
-----

                       Changes in version 1.26.7                        

- Some URLs are reformatted in the documentations and vignettes.

- Data is only loaded into the function environment with `data(...,
  envir = environment())`.

- Namespace imports are adjusted to fixed the function name conflict.

[zellkonverter](/packages/zellkonverter)
-------------

                       Changes in version 1.12.0                        

Major changes

- 
  Add environments for *anndata* v0.9.2 and v0.10.2. Version
  0.10.20 is now the default envrionment for the Python
  reader/writer.

Minor changes

- 
  Changes for compatibility with *rhdf5* v2.45.1 including enum
  types that simplifies reading of nullable types in the native R
  reader

- 
  Dimensions are now passed correctly when converting the raw
  slot

- 
  Backed sparse matrices are now converted in AnnData2SCE()

[zenith](/packages/zenith)
------

                        Changes in version 1.3.1                        

- user can specify organism

NEWS from existing Data Experiment Packages
===================================


[cfToolsData](/packages/cfToolsData)
-----------

                        Changes in version 1.0.0                        

- New Package Release

- 1st version of the package

[curatedTCGAData](/packages/curatedTCGAData)
---------------

                       Changes in version 1.24.0                        

Bug fixes and minor improvements

- Create an on-the-fly sampleMap for RNASeq2GeneNorm* data version
  2.1.1. Data source has munged colnames and sample maps were not
  updated in the latest upload (#59, @LiNk-NY)

[gDRtestData](/packages/gDRtestData)
-----------

                Changes in version 0.99.21 (2023-09-18)                 

- adjust NEWS to Bioc format

                Changes in version 0.99.20 (2023-09-04)                 

- update testdata

                Changes in version 0.99.19 (2023-06-23)                 

- BiocStyle added to dependency

                Changes in version 0.99.18 (2023-06-22)                 

- replaced rds files with qs

                Changes in version 0.99.17 (2023-06-15)                 

- switch from merge to [[

                Changes in version 0.99.16 (2023-05-29)                 

- update datasets

                Changes in version 0.99.15 (2023-05-24)                 

- format the vignette with BiocStyle

                Changes in version 0.99.14 (2023-05-15)                 

- fix related with data.table

                Changes in version 0.99.13 (2023-05-15)                 

- update testdata as per changes in excess assays

                Changes in version 0.99.12 (2023-05-11)                 

- update testdata as per new data model

                Changes in version 0.99.11 (2023-04-25)                 

- changed data.frame to data.table

                Changes in version 0.99.10 (2023-04-24)                 

- removing redundant files

                 Changes in version 0.99.9 (2023-04-20)                 

- fix warning in Bioc check

                 Changes in version 0.99.8 (2023-04-20)                 

- switch to OSI license

                 Changes in version 0.99.7 (2023-04-19)                 

- update testdata

                 Changes in version 0.99.6 (2023-04-18)                 

- moved wrappers to gDRcore

                 Changes in version 0.99.5 (2023-04-17)                 

- update packages version

                 Changes in version 0.99.4 (2023-04-17)                 

- add R 4.2 as a dependency

                 Changes in version 0.99.3 (2023-04-14)                 

- update testdata

                 Changes in version 0.99.2 (2023-04-13)                 

- improve documentation (Bioc compatibility)

                 Changes in version 0.99.1 (2023-04-07)                 

- update maintainer

                 Changes in version 0.99.0 (2023-03-31)                 

- preparing package for Bioc submission

- fix examples

[HCATonsilData](/packages/HCATonsilData)
-------------

                       Changes in version 0.99.0                        

New features

- Dataset access provided for final (v2) data submission

- 7 additional tonsils from young and old adults added during
  revision.

- Improved package vignette to document the new changes.

- Users can now download Visium data into SpatialExperiment object.

- Vignette explains how to download ATAC and Multiome datasets from
  Zenodo as Seurat objects.

- Glossary available for final set of cell types, with links to key
  references.

Other notes

- HCATonsilData is now submitted to Bioconductor!

- Users can still access data from version 1 (preprint)

[msigdb](/packages/msigdb)
------

                       Changes in version 1.10.0                        

- added MSigDB v7.5.1, v2022.1, and v2023.1

- added functions to retrieve pre-computed IDFs

[orthosData](/packages/orthosData)
----------

                       Changes in version 0.99.4                        

- Addition of main package documentaion

                       Changes in version 0.99.3                        

- ASmall additions to manual and vignette

                       Changes in version 0.99.2                        

- Small additions to manual and vignette

                       Changes in version 0.99.1                        

- Update of Zenodo links

                        Changes in version 0.1.0                        

- Initial version

[raerdata](/packages/raerdata)
--------

                 Changes in version 0.99.0 (2023-05-18)                 

- init

[RforProteomics](/packages/RforProteomics)
--------------

                       Changes in version 1.39.2                        

- Clearn up .gitignore and .Rbuildignore (in the hope to fix the
  'figure not found error')

                       Changes in version 1.39.1                        

- Change paths to figures.

                       Changes in version 1.39.0                        

- New Bioc devel version.

[scMultiome](/packages/scMultiome)
----------

                        Changes in version 1.1.1                        

- renamed assay from "logcounts" to "normalizedCounts" in reprogramseq
  and "counts" to "normalizedCounts" in the hematopoiesis, prostate
  and colon datasets

- specify version number of ExperimentHub and remove version number of
  R

[SingleCellMultiModal](/packages/SingleCellMultiModal)
--------------------

                       Changes in version 1.14.0                        

New features

- The ontomap function provides a reference table of ontology IDs and
  cell names by data type available in the package.

- scRNAseq colData added to cord_blood and peripheral_blood datasets
  provided by the CITEseq function. (@drighelli)

Bug fixes and minor improvements

- When using HDF5 as format input in scMultiome, the filtering of file
  paths obtained from ExperimentHub has been fixed.

- Using BiocBaseUtils internally to handle assertions and checks.

[smokingMouse](/packages/smokingMouse)
------------

                       Changes in version 0.99.0                        

NEW FEATURES

- Provides access to smokingMouse project objects. Check here for the
  code and data generation.


[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.13.4                        

NEW FEATURES

- Added fetch_data("spatialDLPFC_Visium_example_subset") which is a
  subset of 3 samples with only the lowres images that can be used for
  example / tutorial purposes.

                       Changes in version 1.13.2                        

NEW FEATURES

- Louise A. Huuki-Myers @lahuuki added a vignette explaining the
  spatial registration process and all related functions. See
  https://github.com/LieberInstitute/spatialLIBD/pull/46 for the full
  pull request.

[TumourMethData](/packages/TumourMethData)
--------------

                       Changes in version 0.99.0                        

- Submission to Bioconductor

NEWS from existing Workflows
===================================

[seqpac](/packages/seqpac)
------

                        Changes in version 1.1.1                        

- First hard release August 2021

                        Changes in version 1.0.4                        

- Major updates to accomodate package specific tests, and
  devtools/BiocCheck

- S4 compatability

- merge_lanes can now merge flowcell lanes

- make_conv can now generate conversion tables between for example
  UCSC, NCBI and Ensembl chromosome names

                        Changes in version 1.0.3                        

- Streamlined PAC generation and annotation

- Vignette update

                        Changes in version 1.0.2                        

- The fundation of functions for sequence-based counting and annotation
  is set.

                        Changes in version 1.0.1                        

- First github version in 2020

- Working version for constructing PAC objects (S3)



Deprecated and Defunct Packages
===============================

Thirty three software packages were removed from this release (after being deprecated
in Bioc 3.17):
alpine, ArrayExpressHTS, ASpediaFI, BiocDockerManager, ChIC, chromswitch,
copynumber, CopywriteR, dasper, epihet, GAPGOM, GeneAccord, genotypeeval,
maanova, metavizr, MethCP, MIGSA, MIMOSA, NanoStringQCPro, NBSplice, netboxr,
NxtIRFcore, ODER, pkgDepTools, PrecisionTrialDrawer, proBatch, proFIA,
pulsedSilac, savR, sigPathway, STAN, TarSeqQC, tscR

Please note:  gcatest and lfa, previously announced as deprecated in 3.17, has been
updated and remain in Bioconductor.

Forty nine software packages are deprecated in this release and will be removed in Bioc 3.19:
baySeq, BGmix, bigPint, biodbMirbase, BioMM, biomvRCNS, Clonality, CSSP, deco,
DeepBlueR, DMRforPairs, exomeCopy, fcoex, gaggle, GCSscore, genbankr, GISPA,
GOsummaries, GRridge, HPAStainR, imageHTS, LineagePulse, logitT, LowMACA,
LPEadj, macat, mAPKL, mbOmic, MEIGOR, Metab, MSstatsSampleSize, multiSight,
netbiov, OmicsLonDA, PFP, plethy, pwrEWAS, qrqc, Ringo, RNAdecay, SCATE, SEPIRA,
seqbias, seqCNA, SISPA, snapCGH, sscore, Travel, trena


Two experimental data packages were removed from this release (after being
deprecated in BioC 3.17):
alpineData, plasFIA

Ten experimental data packages are deprecated in this release and will be
removed in Bioc 3.19:
ccTutorial, ChIC.data, DLBCL, mAPKLData, MAQCsubsetILM, MIGSAdata, pwrEWAS.data,
SCATEData, seqCNA.annot, stjudem

One annotation packages was removed from this release (after being deprecated
in Bioc 3.17).
MafH5.gnomAD.v3.1.1.GRCh38

No annotation packages were deprecated in this release and will be removed in
Bioc 3.19.

No workflow packages were removed from this release (after being deprecated in
Bioc 3.17).

No workflow packages were deprecated in this release and will be removed in
3.19.

No books were removed from this release (after being deprecated in
Bioc 3.17).

No books were deprecated in this release and will be removed in
3.19.
