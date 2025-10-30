October 30, 2025

**Bioconductor:**

We are pleased to announce Bioconductor 3.22, consisting of
2361 software packages, 435 experiment data packages, 926 annotation
packages, 29 workflows and 6 books.

There are 59 new software packages, 6 new data experiment packages,
2 new annotation packages, no new workflows, 1 new book, and many updates and
improvements to existing packages.

Bioconductor 3.22 is compatible with R 4.5, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher, macOS arm64 and Linux
arm64. This release will also include updated Bioconductor [Docker containers][2].

Note: Currently Bioconductor does not have active daily windows builders. We
will be testing the windows binaries generated through the r-universe
system. These updated windows binaries will be available shortly after this
release. We appreciate your patience as we make them available.

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.22](#getting-started-with-bioconductor-322)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.22
======================================

To update to or install Bioconductor 3.22

1. Install R 4.5. Bioconductor 3.22 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 59 new software packages in this release of Bioconductor.

- [anansi](/packages/anansi) Studies including both microbiome and
  metabolomics data are becoming more common. Often, it would be
  helpful to integrate both datasets in order to see if they
  corroborate each others patterns. All vs all association is
  imprecise and likely to yield spurious associations. This package
  takes a knowledge-based approach to constrain association search
  space, only considering metabolite-function pairs that have been
  recorded in a pathway database. This package also provides a
  framework to assess differential association.

- [anglemania](/packages/anglemania) anglemania extracts genes from
  multi-batch scRNA-seq experiments for downstream dataset
  integration. It shows improvement over the conventional usage of
  highly-variable genes for many integration tasks. We leverage
  gene-gene correlations that are stable across batches to identify
  biologically informative genes which are less affected by batch
  effects. Currently, its main use is for single-cell RNA-seq dataset
  integration, but it can be applied for other multi-batch downstream
  analyses such as NMF.

- [anndataR](/packages/anndataR) Bring the power and flexibility of
  AnnData to the R ecosystem, allowing you to effortlessly manipulate
  and analyse your single-cell data. This package lets you work with
  backed h5ad and zarr files, directly access various slots (e.g. X,
  obs, var), or convert the data into SingleCellExperiment and Seurat
  objects.

- [AWAggregator](/packages/AWAggregator) This package implements an
  attribute-weighted aggregation algorithm which leverages
  peptide-spectrum match (PSM) attributes to provide a more accurate
  estimate of protein abundance compared to conventional aggregation
  methods. This algorithm employs pre-trained random forest models to
  predict the quantitative inaccuracy of PSMs based on their
  attributes. PSMs are then aggregated to the protein level using a
  weighted average, taking the predicted inaccuracy into account.
  Additionally, the package allows users to construct their own
  training sets that are more relevant to their specific experimental
  conditions if desired.

- [batchCorr](/packages/batchCorr) From the perspective of
  metabolites as the continuation of the central dogma of biology,
  metabolomics provides the closest link to many phenotypes of
  interest. This makes metabolomics research promising in teasing
  apart the complexities of living systems. However, due to
  experimental reasons, the data includes non-biological variation
  which limits quality and reproducibility, especially if the data is
  obtained from several batches. The batchCorr package reduces
  unwanted variation by way of between-batch alignment, within-batch
  drift correction and between-batch normalization using
  batch-specific quality control samples and long-term reference QC
  samples. Please see the associated article for more thorough
  descriptions of algorithms.

- [Bioc.gff](/packages/Bioc.gff) Parse GFF and GTF files using C++
  classes. The package also provides utilities to read and write GFF3
  files. The GFF (General Feature Format) format is a tab-delimited
  file format for describing genes and other features of DNA, RNA,
  and protein sequences. GFF files are often used to describe the
  features of genomes.

- [blase](/packages/blase) BLASE is a method for finding where bulk
  RNA-seq data lies on a single-cell pseudotime trajectory. It uses a
  fast and understandable approach based on Spearman correlation,
  with bootstrapping to provide confidence. BLASE can be used to
  "date" bulk RNA-seq data, annotate cell types in scRNA-seq, and
  help correct for developmental phenotype differences in bulk
  RNA-seq experiments.

- [CalibraCurve](/packages/CalibraCurve) CalibraCurve is a
  computational tool designed to generate calibration curves for
  targeted mass spectrometry-based quantitative data. It is
  applicable to various omics disciplines, including proteomics,
  lipidomics, and metabolomics. The package also offers
  functionalities for data and calibration curve visualization and
  concentration prediction from new datasets based on the established
  curves.

- [CBN2Path](/packages/CBN2Path) CBN2Path package provides a unifying
  interface to facilitate CBN-based quantification, analysis and
  visualization of cancer progression pathways.

- [cellmig](/packages/cellmig) High-throughput cell imaging
  facilitates the analysis of cell migration across many wells
  treated under different biological conditions. These workflows
  generate considerable technical noise and biological variability,
  and therefore technical and biological replicates are necessary,
  leading to large, hierarchically structured datasets, i.e., cells
  are nested within technical replicates that are nested within
  biological replicates. Current statistical analyses of such data
  usually ignore the hierarchical structure of the data and fail to
  explicitly quantify uncertainty arising from technical or
  biological variability. To address this gap, we present cellmig, an
  R package implementing Bayesian hierarchical models for migration
  analysis. cellmig quantifies condition- specific velocity changes
  (e.g., drug effects) while modeling nested data structures and
  technical artifacts. It further enables synthetic data generation
  for experimental design optimization.

- [Chromatograms](/packages/Chromatograms) The Chromatograms packages
  defines an efficient infrastructure for storing and handling of
  chromatographic mass spectrometry data. It provides different
  implementations of *backends* to store and represent the data. Such
  backends can be optimized for small memory footprint or fast data
  access/processing. A lazy evaluation queue and chunk-wise
  processing capabilities ensure efficient analysis of also very
  large data sets.

- [cigarillo](/packages/cigarillo) CIGAR stands for Concise
  Idiosyncratic Gapped Alignment Report. CIGAR strings are found in
  the BAM files produced by most aligners and in the AIRR-formatted
  output produced by IgBLAST. The cigarillo package provides
  functions to parse and inspect CIGAR strings, trim them, turn them
  into ranges of positions relative to the "query space" or
  "reference space", and project positions or sequences from one
  space to the other. Note that these operations are low-level
  operations that the user rarely needs to perform directly. More
  typically, they are performed behind the scene by higher-level
  functionality implemented in other packages like Bioconductor
  packages GenomicAlignments and igblastr.

- [Coralysis](/packages/Coralysis) Coralysis is an R package
  featuring a multi-level integration algorithm for sensitive
  integration, reference-mapping, and cell-state identification in
  single-cell data. The multi-level integration algorithm is inspired
  by the process of assembling a puzzle - where one begins by
  grouping pieces based on low-to high-level features, such as color
  and shading, before looking into shape and patterns. This approach
  progressively blends the batch effects and separates cell types
  across multiple rounds of divisive clustering.

- [CSOA](/packages/CSOA) Cell Set Overlap Analysis (CSOA) is a tool
  for calculating per-cell gene signature scores in an scRNA-seq
  dataset. CSOA constructs a set for each gene in the signature,
  consisting of the cells that highly express the gene. Next, all
  overlaps of pairs of cell sets are computed, ranked, filtered and
  scored. The CSOA per-cell score is calculated by summing up all
  products of the overlap scores and the min-max-normalized
  expression of the two involved genes. CSOA can run on a Seurat
  object, a SingleCellExperiment object, a matrix and a dgCMatrix.

- [DeeDeeExperiment](/packages/DeeDeeExperiment) DeeDeeExperiment is
  an S4 class extending the SingleCellExperiment class, designed to
  integrate and manage omics analysis results. It introduces two
  dedicated slots to store Differential Expression Analysis (DEA)
  results and Functional Enrichment Analysis (FEA) results, providing
  a structured approach for downstream analysis.

- [dmGsea](/packages/dmGsea) The R package dmGsea provides efficient
  gene set enrichment analysis specifically for DNA methylation data.
  It addresses key biases, including probe dependency and varying
  probe numbers per gene. The package supports Illumina 450K, EPIC,
  and mouse methylation arrays. Users can also apply it to other
  omics data by supplying custom probe-to-gene mapping annotations.
  dmGsea is flexible, fast, and well-suited for large-scale
  epigenomic studies.

- [DNEA](/packages/DNEA) The DNEA R package is the latest
  implementation of the Differential Network Enrichment Analysis
  algorithm and is the successor to the Filigree Java-application
  described in Iyer et al. (2020). The package is designed to take as
  input an m x n expression matrix for some -omics modality (ie.
  metabolomics, lipidomics, proteomics, etc.) and jointly estimate
  the biological network associations of each condition using the
  DNEA algorithm described in Ma et al. (2019). This approach
  provides a framework for data-driven enrichment analysis across two
  experimental conditions that utilizes the underlying correlation
  structure of the data to determine feature-feature interactions.

- [DOtools](/packages/DOtools) This package provides functions for
  creating various visualizations, convenient wrappers, and
  quality-of-life utilities for single cell experiment objects. It
  offers a streamlined approach to visualize results and integrates
  different tools for easy use.

- [DspikeIn](/packages/DspikeIn) Provides a reproducible and modular
  workflow for absolute microbial quantification using spike-in
  controls. Supports both single spike-in taxa and synthetic
  microbial communities with user-defined spike-in volumes and genome
  copy numbers. Compatible with 'phyloseq' and
  'TreeSummarizedExperiment' (TSE) data structures. The package
  implements methods for spike-in validation, preprocessing, scaling
  factor estimation, absolute abundance conversion, bias correction,
  and normalization. Facilitates downstream statistical analyses with
  'DESeq2', 'edgeR', and other Bioconductor-compatible methods.
  Visualization tools are provided via 'ggplot2', 'ggtree', and
  related packages. Includes detailed vignettes, case studies, and
  function-level documentation to guide users through experimental
  design, quantification, and interpretation.

- [FinfoMDS](/packages/FinfoMDS) F-informed MDS is a new
  multidimensional scaling-based ordination method that configures
  data distribution based on the F-statistic (i.e., the ratio of
  dispersion between groups with shared or differing labels).

- [GCPtools](/packages/GCPtools) Lower-level functionality to
  interface with Google Cloud Platform tools. 'gcloud' and 'gsutil'
  are both supported. The functionality provided centers around
  utilities for the AnVIL platform.

- [goatea](/packages/goatea) Geneset Ordinal Association Test
  Enrichment Analysis (GOATEA) provides a 'Shiny' interface with
  interactive visualizations and utility functions for performing and
  exploring automated gene set enrichment analysis using the 'GOAT'
  package. 'GOATEA' is designed to support large-scale and
  user-friendly enrichment workflows across multiple gene lists and
  comparisons, with flexible plotting and output options.
  Visualizations pre-enrichment include interactive 'Volcano' and
  'UpSet' (overlap) plots. Visualizations post-enrichment include
  interactive geneset dotplot, geneset treeplot, gene-effectsize
  heatmap, gene-geneset heatmap and 'STRING' database of
  protein-protein-interactions network graph. 'GOAT' reference: Frank
  Koopmans (2024) <doi:10.1038/s42003-024-06454-5>.

- [gVenn](/packages/gVenn) Tools to compute and visualize overlaps
  between gene sets or genomic regions. Venn diagrams with
  proportional areas are provided, while UpSet plots are recommended
  for larger numbers of sets. The package supports GRanges and
  GRangesList inputs, and integrates with analysis workflows for
  ChIP-seq, ATAC-seq, and other genomic interval data. It generates
  clean, interpretable, and publication-ready figures.

- [HiCaptuRe](/packages/HiCaptuRe) Capture Hi-C is a set of
  techniques that enable the detection of genomic interactions
  involving regions of interest, known as baits. By focusing on
  selected loci, these approaches reduce sequencing costs while
  maintaining high resolution at the level of restriction fragments.
  HiCaptuRe provides tools to import, annotate, manipulate, and
  export Capture Hi-C data. The package accounts for the specific
  structure of bait–otherEnd interactions, facilitates integration
  with other omics datasets, and enables comparison across samples
  and conditions.

- [HiCPotts](/packages/HiCPotts) The HiCPotts package provides a
  comprehensive Bayesian framework for analyzing Hi-C interaction
  data, integrating both spatial and genomic biases within a
  probabilistic modeling framework. At its core, HiCPotts leverages
  the Potts model (Wu, 1982)—a well-established graphical model—to
  capture and quantify spatial dependencies across interaction loci
  arranged on a genomic lattice. By treating each interaction as a
  spatially correlated random variable, the Potts model enables
  robust segmentation of the genomic landscape into meaningful
  components, such as noise, true signals, and false signals. To
  model the influence of various genomic biases, HiCPotts employs a
  regression-based approach incorporating multiple covariates:
  Genomic distance (D): The distance between interacting loci,
  recognized as a fundamental driver of contact frequency. GC-content
  (GC): The local GC composition around the interacting loci, which
  can influence chromatin structure and interaction patterns.
  Transposable elements (TEs): The presence and abundance of
  repetitive elements that may shape contact probability through
  chromatin organization. Accessibility score (Acc): A measure of
  chromatin openness, informing how accessible certain genomic
  regions are to interaction. By embedding these covariates into a
  hierarchical mixture model, HiCPotts characterizes each
  interaction’s probability of belonging to one of several latent
  components. The model parameters, including regression
  coefficients, zero-inflation parameters (for ZIP/ZINB
  distributions), and dispersion terms (for NB/ZINB distributions),
  are inferred via a MCMC sampler. This algorithm draws samples from
  the joint posterior distribution, allowing for flexible posterior
  inference on model parameters and hidden states. From these
  posterior samples, HiCPotts computes posterior means of regression
  parameters and other quantities of interest. These posterior
  estimates are then used to calculate the posterior probabilities
  that assign each interaction to a specific component. The resulting
  classification sheds light on the underlying structure:
  distinguishing genuine high-confidence interactions (signal) from
  background noise and potential false signals, while simultaneously
  quantifying the impact of genomic biases on observed interaction
  frequencies. In summary, HiCPotts seamlessly integrates spatial
  modeling, bias correction, and probabilistic classification into a
  unified Bayesian inference framework. It provides rich posterior
  summaries and interpretable, model-based assignments of interaction
  states, enabling researchers to better understand the interplay
  between genomic organization, biases, and spatial correlation in
  Hi-C data.

- [HVP](/packages/HVP) HVP is a quantitative batch effect metric that
  estimates the proportion of variance associated with batch effects
  in a data set.

- [Ibex](/packages/Ibex) Implementation of the Ibex algorithm for
  single-cell embedding based on BCR sequences. The package includes
  a standalone function to encode BCR sequence information by amino
  acid properties or sequence order using tensorflow-based
  autoencoder. In addition, the package interacts with
  SingleCellExperiment or Seurat data objects.

- [igblastr](/packages/igblastr) The igblastr package provides
  functions to conveniently install and use a local IgBLAST
  installation from within R. IgBLAST is described at
  <https://pubmed.ncbi.nlm.nih.gov/23671333/>. Online IgBLAST:
  <https://www.ncbi.nlm.nih.gov/igblast/>.

- [iModMix](/packages/iModMix) The iModMix network-based method
  offers an integrated framework for analyzing multi-omics data,
  including metabolomics, proteomics, and transcriptomics data,
  enabling the exploration of intricate molecular associations within
  heterogeneous biological systems.

- [iscream](/packages/iscream) BED files store ranged genomic data
  that can be queried even when the files are compressed. iscream can
  query data from BED files and return them in muliple formats:
  parsed records or their summary statistics as data frames or
  GenomicRanges objects, and matrices as matrix, GenomicRanges, or
  SummarizedExperiment objects. iscream also provides specialized
  support for importing methylation data.

- [linkSet](/packages/linkSet) Provides a comprehensive framework for
  representing, analyzing, and visualizing genomic interactions,
  particularly focusing on gene-enhancer relationships. The package
  extends the GenomicRanges infrastructure to handle paired genomic
  regions with specialized methods for chromatin interaction data
  from Hi-C, Promoter Capture Hi-C (PCHi-C), and single-cell ATAC-seq
  experiments. Key features include conversion from common
  interaction formats, annotation of promoters and enhancers,
  distance-based analyses, interaction strength metrics, statistical
  modeling using CHiCANE methodology, and tailored visualization
  tools. The package aims to standardize the representation of
  genomic interaction data while providing domain-specific functions
  not available in general genomic interaction packages.

- [LipidTrend](/packages/LipidTrend) "LipidTrend" is an R package
  that implements a permutation-based statistical test to identify
  significant differences in lipidomic features between groups. The
  test incorporates Gaussian kernel smoothing of region statistics to
  improve stability and accuracy, particularly when dealing with
  small sample sizes. This package also includes two plotting
  functions for visualizing significant tendencies in 1D and 2D
  feature data, respectively.

- [looking4clusters](/packages/looking4clusters) Enables the
  interactive visualization of dimensional reduction, clustering, and
  cell properties for scRNA-Seq results. It generates an interactive
  HTML page using either a numeric matrix, SummarizedExperiment,
  SingleCellExperiment or Seurat objects as input. The input data can
  be projected into two-dimensional representations by applying
  dimensionality reduction methods such as PCA, MDS, t-SNE, UMAP, and
  NMF. Displaying multiple dimensionality reduction results within
  the same interface, with interconnected graphs, provides different
  perspectives that facilitate accurate cell classification. The
  package also integrates unsupervised clustering techniques, whose
  results that can be viewed interactively in the graphical
  interface. In addition to visualization, this interface allows
  manual selection of groups, labeling of cell entities based on
  processed meta-information, generation of new graphs displaying
  gene expression values for each cell, sample identification, and
  visual comparison of samples and clusters.

- [markeR](/packages/markeR) markeR is an R package that provides a
  modular and extensible framework for the systematic evaluation of
  gene sets as phenotypic markers using transcriptomic data. The
  package is designed to support both quantitative analyses and
  visual exploration of gene set behaviour across experimental and
  clinical phenotypes. It implements multiple methods, including
  score-based and enrichment approaches, and also allows the
  exploration of expression behaviour of individual genes. In
  addition, users can assess the similarity of their own gene sets
  against established collections (e.g., those from MSigDB),
  facilitating biological interpretation.

- [MetaDICT](/packages/MetaDICT) MetaDICT is a method for the
  integration of microbiome data. This method is designed to remove
  batch effects and preserve biological variation while integrating
  heterogeneous datasets. MetaDICT can better avoid overcorrection
  when unobserved confounding variables are present.

- [miaTime](/packages/miaTime) The `miaTime` package provides tools
  for microbiome time series analysis based on
  (Tree)SummarizedExperiment infrastructure.

- [MSstatsResponse](/packages/MSstatsResponse) Tools for detecting
  drug-protein interactions and estimating IC50 values from
  chemoproteomics data. Implements semi-parametric isotonic
  regression, bootstrapping, and curve fitting to evaluate compound
  effects on protein abundance.

- [mutscan](/packages/mutscan) Provides functionality for processing
  and statistical analysis of multiplexed assays of variant effect
  (MAVE) and similar data. The package contains functions covering
  the full workflow from raw FASTQ files to publication-ready
  visualizations. A broad range of library designs can be processed
  with a single, unified interface.

- [notame](/packages/notame) Provides functionality for untargeted
  LC-MS metabolomics research as specified in the associated protocol
  article in the 'Metabolomics Data Processing and Data
  Analysis—Current Best Practices' special issue of the Metabolites
  journal (2020). This includes tabular data preprocessing and
  quality control, uni- and multivariate analysis as well as quality
  control visualizations, feature-wise visualizations and results
  visualizations. Raw data preprocessing and functionality related to
  biological context, such as pathway analysis, is not included.

- [notameStats](/packages/notameStats) Provides univariate and
  multivariate statistics for feature prioritization in untargeted
  LC-MS metabolomics research.

- [notameViz](/packages/notameViz) Provides visualization
  functionality for untargeted LC-MS metabolomics research. Includes
  quality control visualizations, feature-wise visualizations and
  results visualizations.

- [omicsGMF](/packages/omicsGMF) omicsGMF is a Bioconductor package
  that uses the sgdGMF-framework of the \code{sgdGMF} package for
  highly performant and fast matrix factorization that can be used
  for dimensionality reduction, visualization and imputation of omics
  data. It considers data from the general exponential family as
  input, and therefore suits the use of both RNA-seq (Poisson or
  Negative Binomial data) and proteomics data (Gaussian data). It
  does not require prior transformation of counts to the log-scale,
  because it rather optimizes the deviances from the data family
  specified. Also, it allows to correct for known sample-level and
  feature-level covariates, therefore enabling visualization and
  dimensionality reduction upon batch correction. Last but not least,
  it deals with missing values, and allows to impute these after
  matrix factorization, useful for proteomics data. This Bioconductor
  package allows input of SummarizedExperiment, SingleCellExperiment,
  and QFeature classes.

- [peakCombiner](/packages/peakCombiner) peakCombiner, a fully R
  based, user-friendly, transparent, and customizable tool that
  allows even novice R users to create a high-quality consensus peak
  list. The modularity of its functions allows an easy way to
  optimize input and output data. A broad range of accepted input
  data formats can be used to create a consensus peak set that can be
  exported to a file or used as the starting point for most
  downstream peak analyses.

- [PMScanR](/packages/PMScanR) Provides tools for large-scale protein
  motif analysis and visualization in R. PMScanR facilitates the
  identification of motifs using external tools like PROSITE's
  ps_scan (handling necessary file downloads and execution) and
  enables downstream analysis of results. Key features include
  parsing scan outputs, converting formats (e.g., to GFF-like
  structures), generating motif occurrence matrices, and creating
  informative visualizations such as heatmaps, sequence logos (via
  seqLogo/ggseqlogo). The package also offers an optional Shiny-based
  graphical user interface for interactive analysis, aiming to
  streamline the process of exploring motif patterns across multiple
  protein sequences.

- [SanityR](/packages/SanityR) a Bayesian normalization procedure
  derived from first principles. Sanity estimates expression values
  and associated error bars directly from raw unique molecular
  identifier (UMI) counts without any tunable parameters.

- [scafari](/packages/scafari) Scafari is a Shiny application
  designed for the analysis of single-cell DNA sequencing (scDNA-seq)
  data provided in .h5 file format. The analysis process is
  structured into the four key steps "Sequencing", "Panel",
  "Variants", and "Explore Variants". It supports various analyses
  and visualizations.

- [scGraphVerse](/packages/scGraphVerse) A package for inferring,
  comparing, and visualizing gene networks from single-cell RNA
  sequencing data. It integrates multiple methods (GENIE3, GRNBoost2,
  ZILGM, PCzinb, and JRF) for robust network inference, supports
  consensus building across methods or datasets, and provides tools
  for evaluating regulatory structure and community similarity.
  GRNBoost2 requires Python package 'arboreto' which can be installed
  using init_py(install_missing = TRUE). This package includes
  adapted functions from ZILGM (Park et al., 2021), JRF (Petralia et
  al., 2015), and learn2count (Nguyen et al. 2023) packages with
  proper attribution under GPL-2 license.

- [scLANE](/packages/scLANE) Our scLANE model uses truncated power
  basis spline models to build flexible, interpretable models of
  single cell gene expression over pseudotime or latent time. The
  modeling architectures currently supported are Negative-binomial
  GLMs, GEEs, & GLMMs. Downstream analysis functionalities include
  model comparison, dynamic gene clustering, smoothed counts
  generation, gene set enrichment testing, & visualization.

- [Seqinfo](/packages/Seqinfo) The Seqinfo class stores the names,
  lengths, circularity flags, and genomes for a particular collection
  of sequences. These sequences are typically the chromosomes and/or
  scaffolds of a specific genome assembly of a given organism.
  Seqinfo objects are rarely used as standalone objects. Instead,
  they are used as part of higher-level objects to represent their
  seqinfo() component. Examples of such higher-level objects are
  GRanges, RangedSummarizedExperiment, VCF, GAlignments, etc...
  defined in other Bioconductor infrastructure packages.

- [SETA](/packages/SETA) Tools for compositional and other
  sample-level ecological analyses and visualizations tailored for
  single-cell RNA-seq data. SETA includes functions for taxonomizing
  celltypes, normalizing data, performing statistical tests, and
  visualizing results. Several tutorials are included to guide users
  and introduce them to key concepts. SETA is meant to teach users
  about statistical concepts underlying ecological analysis methods
  so they can apply them to their own single-cell data.

- [shinybiocloader](/packages/shinybiocloader) Add a Bioconductor
  themed CSS loader to your shiny app. It is based on the
  shinycustomloader R package. Use a spinning Bioconductor note
  loader to enhance your shiny app loading screen. This package is
  intended for developer use.

- [SmartPhos](/packages/SmartPhos) To facilitate and streamline
  phosphoproteomics data analysis, we developed SmartPhos, an R
  package for the pre-processing, quality control, and exploratory
  analysis of phosphoproteomics data generated by MaxQuant and
  Spectronaut. The package can be used either through the R command
  line or through an interactive ShinyApp called SmartPhos Explorer.
  The package contains methods such as normalization and
  normalization correction, transformation, imputation, batch effect
  correction, PCA, heatmap, differential expression, time-series
  clustering, gene set enrichment analysis, and kinase activity
  inference.

- [SpaceTrooper](/packages/SpaceTrooper) SpaceTrooper performs
  Quality Control analysis using data driven GLM models of
  Image-Based spatial data, providing exploration plots, QC metrics
  computation, outlier detection. It implements a GLM strategy for
  the detection of low quality cells in imaging-based spatial data
  (Transcriptomics and Proteomics). It additionally implements
  several plots for the visualization of imaging based polygons
  through the ggplot2 package.

- [spARI](/packages/spARI) The R package used in the manuscript
  "Spatially Aware Adjusted Rand Index for Evaluating Spatial
  Transcritpomics Clustering".

- [SpectriPy](/packages/SpectriPy) The SpectriPy package allows
  integration of Python-based MS analysis code with the Spectra
  package. Spectra objects can be converted into Python MS data
  structures. In addition, SpectriPy integrates and wraps the
  similarity scoring and processing/filtering functions from the
  Python matchms package into R.

- [SPICEY](/packages/SPICEY) SPICEY (SPecificity Index for Coding and
  Epigenetic activitY) is an R package designed to quantify cell-type
  specificity in single-cell transcriptomic and epigenomic data,
  particularly scRNA-seq and scATAC-seq. It introduces two
  complementary indices: the Gene Expression Tissue Specificity Index
  (GETSI) and the Regulatory Element Tissue Specificity Index
  (RETSI), both based on entropy to provide continuous, interpretable
  measures of specificity. By integrating gene expression and
  chromatin accessibility, SPICEY enables standardized analysis of
  cell-type-specific regulatory programs across diverse tissues and
  conditions.

- [STADyUM](/packages/STADyUM) STADyUM is a package with
  functionality for analyzing nascent RNA read counts to infer
  transcription rates. This includes utilities for processing
  experimental nascent RNA read counts as well as for simulating
  PRO-seq data. Rates such as initiation, pause release and landing
  pad occupancy are estimated from either synthetic or experimental
  data. There are also options for varying pause sites and including
  steric hindrance of initiation in the model.

- [stPipe](/packages/stPipe) This package serves as an upstream
  pipeline for pre-processing sequencing-based spatial
  transcriptomics data. Functions includes FASTQ trimming, BAM file
  reformatting, index building, spatial barcode detection,
  demultiplexing, gene count matrix generation with UMI
  deduplication, QC, and revelant visualization. Config is an
  essential input for most of the functions which aims to improve
  reproducibility.

- [SuperCellCyto](/packages/SuperCellCyto) SuperCellCyto provides the
  ability to summarise cytometry data into supercells by merging
  together cells that are similar in their marker expressions using
  the SuperCell package.


New Data Experiment Packages
=====================

There are 6 new data experiment packages in this release of Bioconductor.

- [AWAggregatorData](/packages/AWAggregatorData) The AWAggregatorData
  package contains the data associated with the AWAggregator R
  package. It includes two pre-trained random forest models, one
  incorporating the average coefficient of variation as a feature,
  and the other one not including it. It also contains the PSMs in
  Benchmark Set 1~3 derived from the psm.tsv output files generated
  by FragPipe, which are used to train the random forest models.

- [CENTREprecomputed](/packages/CENTREprecomputed) Interface and
  documentation for the Experiment Hub records needed by the CENTRE
  Bioconductor software package. The Experiment Hub records contains
  the precomputed fisher combined p-values, CRUP correlations.
  Additionally, the records hold ChIP-seq and RNA-seq data used for
  the example of the software package.

- [ChIPDBData](/packages/ChIPDBData) Provides curated gene target
  databases derived from ChIP-seq datasets, formatted as ChIPDB
  objects for use with TFEA.ChIP.

- [DoReMiTra](/packages/DoReMiTra) DoReMiTra is an R data package
  providing access to curated transcriptomic datasets related to
  blood radiation, with a focus on neutron, x-ray, and gamma ray
  studies. It is designed to facilitate radiation biology research
  and support data exploration and reproducibility in radiation
  transcriptomics. All datasets are provided as SummarizedExperiment
  objects, allowing seamless integration with the Bioconductor
  ecosystem.

- [iModMixData](/packages/iModMixData) Provides example datasets for
  the iModMix package, including gene, protein, and metabolite
  partial correlation matrices derived from ccRCC4 and
  FloresData_K_TK studies. The data are preprocessed and ready to use
  for testing, demonstrating iModMix workflows, and exploring
  correlation networks.

- [nmrdata](/packages/nmrdata) Provides example one-dimensional
  proton NMR spectra of murine urine samples collected before and
  after bariatric or sham surgery (Roux-en-Y gastric bypass). The
  data are adapted from Jia V Li et al. (2011), "Metabolic surgery
  profoundly influences gut microbial-host metabolic cross-talk",
  Gut, 60(9), 1214–1223. <doi:10.1136/gut.2010.234708>. This package
  serves as example data for metabolomics analysis and teaching
  purposes.


New Annotation Packages
=====================

There are 2 new annotation packages in this release of Bioconductor.

- [CENTREannotation](/packages/CENTREannotation) This is an
  AnnotationHub package for the CENTRE Bioconductor software package.
  It contains the GENCODE version 40 annotation and ENCODE Registry
  of candidate cis-regulatory elements (cCREs) version 3. All for
  Human hg38 genome.

- [org.Hbacteriophora.eg.db](/packages/org.Hbacteriophora.eg.db)
  Provides genome-wide annotation for Heterorhabditis bacteriophora,
  primarily based on mapping using custom gene identifiers. This
  OrgDb annotation package is intended for use with
  AnnotationDbi-based tools and supports querying of gene identifiers
  and related metadata.

New Workflow Packages
=====================

There are no new workflow packages in this release of Bioconductor.

New Online Books
=====================

There is one new books in this release of Bioconductor.

- [OSTA](/packages/OSTA) This package contains source files for the
  "Orchestrating Spatial Transcriptomics Analysis with Bioconductor"
  online book. This book provides interactive examples and discussion
  on key principles of computational analysis workflows for spatial
  transcriptomics data using Bioconductor in R. The book contains
  chapters describing individual analysis steps as well as extended
  workflows, each with examples including R code and datasets.


NEWS from existing Software Packages
===================================


[alabaster.base](/packages/alabaster.base)
--------------

                       Changes in version 1.10.0                        

- Bugfix to cloneDirectory() to work with relative/recursive
  symlinks in the directory to be cloned.

[alevinQC](/packages/alevinQC)
--------

                       Changes in version 1.25.1                        

- Minor fixes to tests to be more robust to possible future changes
in
ggplot2

[AlphaMissenseR](/packages/AlphaMissenseR)
--------------

                        Changes in version 1.6.0                        

- (v. 1.5.4) Avoid gghalves in 'benchmarking' vignette for
compatilibity with updated ggplot2.
- (v. 1.5.4) ensure af_prediction() returns a single entry exactly
matching accession number.

[AlpsNMR](/packages/AlpsNMR)
-------

                 Changes in version 4.11.1 (2025-09-24)                 

- Compatibility with ggplot2-4.0.

[amplican](/packages/amplican)
--------

                       Changes in version 1.31.2                        

- major changes to amplican way of finding primers, it should be way
  faster and more robust

- some small fixes to make findLQR more robust

[anansi](/packages/anansi)
------

                       Changes in version 0.99.0                        

Submission to Bioconductor

[anglemania](/packages/anglemania)
----------

                       Changes in version 0.99.5                        

- Initial submission to Bioconductor.

[anndataR](/packages/anndataR)
--------

- Initial submission to Bioconductor

[AnnotationHub](/packages/AnnotationHub)
-------------

                       Changes in version 3.99.0                        

MAJOR UPDATES

- 3.99.0 Update downloading method from httr to httr2

BUG CORRECTION

- 3.99.1 Fix bug preventing status message with removeResource

- 3.99.2 Add ability to pass a config argument to download resources

- 3.99.3 Add ability to pass a progress argument to download resources

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.33.1                        

- add option 'ATACseqQC.bigFile'

[AWAggregator](/packages/AWAggregator)
------------

                       Changes in version 0.99.1                        

- Initial Bioconductor submission


[BANDITS](/packages/BANDITS)
-------

                       Changes in version 1.24.1                        

- 'create_data' updated to prevent users from providing an unnamed
  'eff_len' object;

[barbieQ](/packages/barbieQ)
-------

                        Changes in version 1.1.2                        

- Comment out code that runs bartools functions in vignette.

                        Changes in version 1.1.1                        

- Fixed issues in the vignette.

- Fixed bug in the "testBarcodeSignif.R".

- Added more examples to "testBarcodeSignif.R".

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.21.1 (2025-09-17)                 

- Updates divide and conquer approach to use an effect size as a metric
  to
  decide whether the partions are balanced (instead of a p-value)

- Minor typos and technical fixes

[basilisk](/packages/basilisk)
--------

                       Changes in version 1.22.0                        

- Replace conda with PyEnv + virtualenv for Python management and
  environment creation. This aims to avoid GLIBCXX version errors
  from mismatches in system libraries between conda and R.
  Package developers should now directly use PyPI package
  versions in their packages= specification instead of the conda
  package versions. Administrators also have more power and
  flexibility to override the Python instance used by
  setupBasiliskEnv.

[batchCorr](/packages/batchCorr)
---------

                       Changes in version 0.99.5                        

- Submitted to Bioconductor

[BatchQC](/packages/BatchQC)
-------

                       Changes in version 2.5.14                        

Minor Changes

- Added voom as a normalization method

Bug Fixes

- Corrected limma correction to allow 1 covariate

                       Changes in version 2.5.13                        

Major Changes

- Updated Vignette to reflect updates to BatchQC (distribution
checks,
kBET, UMAP)
- Updated vignette example data to be TB data

Minor Changes

- Updated UMAP example to include covar option
- Changed location of kBET tab

                       Changes in version 2.5.12                        

Minor Changes

- Fixed sva batch correction function
- Added covar to UMAP

                       Changes in version 2.5.11                        

Major Changes

- Allowed AIC function to take multiple covariates as input

                       Changes in version 2.5.10                        

Major Changes

- Added a new kBET tab for measuring batch effect presence

                        Changes in version 2.5.9                        

Major Changes

- Added AIC under batch correction tab

                        Changes in version 2.5.8                        

Major Changes

- Added svaseq as a batch correction method

                        Changes in version 2.5.7                        

Minor Changes

- Added a data onject with subject IDs and BMI info to match the
original study with the curatedTBData info and provide missing BMI
info

                        Changes in version 2.5.6                        

Major Changes

- Added Kruskal-Wallis test as a differential expression analysis
method

                        Changes in version 2.5.5                        

Major Changes

- Added ANOVA as a differential expression analysis method

                        Changes in version 2.5.4                        

Major Changes

- Redesigned tab flow; upload tab reduced, batch correction tab
created, and lambda stat and neg binomial distribution check moved
to appropriate tabs

Minor Changes

- Updated lambda statistic to provide a recommendation and to check
for balanced experimental design
- Abbreviated Normalization usage description

Bug Fix

- Updated "counts" language to "assay" for uploaded data

                        Changes in version 2.5.3                        

Major Changes

- Added lambda statistic

Bug Fix

- Corrected log to be log(x+1) and changed CPM(x+1) to CPM(x)

                        Changes in version 2.5.2                        

Major Changes

- Added edgeR as a normalization method
- Added edgeR as a differential expression analysis method

                        Changes in version 2.5.1                        

Major Changes

- Added TB data example

Minor Changes

- Changed sapply to vapply in nb check

[beachmat](/packages/beachmat)
--------

                       Changes in version 2.26.0                        

- Added the getExecutor() and Rtatami::set_executor() functions.
  This ensures that unknown matrices are parallelized safely in
  downstream packages.  Also updated instructions for developers
  to set the parallel executor in their .onLoad().

[bedbaser](/packages/bedbaser)
--------

                        Changes in version 1.1.1                        

- Update bedbase to 0.10.6
- Update examples
- Remove bigbed

[bettr](/packages/bettr)
-----

                        Changes in version 1.5.3                        

- Expand vignette with more details about the app

                        Changes in version 1.5.2                        

- Add possibility to specify the default weight for metrics

                        Changes in version 1.5.1                        

- Minor fixes to tests to be more robust to possible future changes
in
ggplot2

[Bioc.gff](/packages/Bioc.gff)
--------

                       Changes in version 0.99.5                        

- Initial Bioconductor submission.

[BiocBaseUtils](/packages/BiocBaseUtils)
-------------

                       Changes in version 1.12.0                        

Bug fixes and minor improvements

- The error message from checkInstalled now reports required packages
in addition to the installation text.

[BiocCheck](/packages/BiocCheck)
---------

                       Changes in version 1.46.0                        

NEW FEATURES

- Checks support `qmd` vignettes and require that a
  `SystemRequirements`
  field be present and list the `quarto` CLI in the DESCRIPTION file.

- All chunks in vignettes should have labels (@vjcitn, #222)

- Flag export all regex patterns in `NAMESPACE` (@lshep, #230)

BUG FIXES AND MINOR IMPROVEMENTS

- `Rd` pages that do not have `\usage` sections are excluded from
  `\value`
  checks (@hpages, #225)

- `Rd` pages that have `\keyword{internal}` are excluded from `\value`
  checks

- Avoid overwritting `data` document types for `\value` checks

- Vignettes with no code chunks were incorrectly flagged in chunk label
  checks

- Chunks in Rnw vignettes with no labels are now correctly identified

- `sessionInfo` checks in vignettes occur only when code chunks are
  present

- BiocCheck() options are now ordered alphabetically in documentation
  to
  make it easier to find specific options (@Bisaloo, #229).

[BiocFileCache](/packages/BiocFileCache)
-------------

                       Changes in version 2.99.0                        

MAJOR UPDATES

- Update downloading method from httr to httr2

BUG FIX

- (2.99.1) Soft failure for if HEAD does not work. Will still download
  and
  add to cache without caching information. Additionally fix bug so
  config can
  be passed down to HEAD call from bfcneedsupdate, bfcadd, and
  bfcrpath.

- (2.99.4) Fix bug to pass proxy down to HEAD call

NEW FEATURE

- (2.99.2) Add progress argument to bfcadd/bfcdownload to control if a
  progress bar is displayed during an interactive session

[BiocGenerics](/packages/BiocGenerics)
------------

                       Changes in version 0.56.0                        

BUG FIXES

- Fix updateObject() infinite recursion on an environment that contains
  itself.

- Small tweak to updateObject(<data.frame>).

[BiocNeighbors](/packages/BiocNeighbors)
-------------

                        Changes in version 2.4.0                        

- Updated the C++ interfaces to match updates to the underlying
  knncolle libraries.

[BioCor](/packages/BioCor)
------

                        Changes in version 1.34                         

- Remove example on vignettes with deprecated package
targetscan.Hs.eg.db
- Fix documentation issue
- Fix deprecation message related to testthat 3rd edition.

[BiocParallel](/packages/BiocParallel)
------------

                        Changes in version 1.44                         

TRANSFERRING MAINTAINERSHIP

- (1.43.4) Transfer maintainership to Jiefei Wang

NEW FEATURES

- (1.43.2) Lazy-load S4 classes in workers; see
  <https://github.com/Bioconductor/BiocParallel/issues/262>
  <https://github.com/Bioconductor/BiocParallel/pull/270>

- (1.43.1) Support NULL return value from bplapply; see
  <https://github.com/Bioconductor/BiocParallel/issues/267>
  <https://github.com/Bioconductor/BiocParallel/pull/269>

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.28.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- `biocVIEWSdb` parses the package metadata VIEWS file

- `pkgBiocRevDeps` now has `recursive` argument

- `repositoryStats` supports in-house repositories with `local=TRUE`

- The dependency report template now uses `pkgndep` to display a
  dependency
  heatmap

BUG FIXES AND MINOR IMPROVEMENTS

- Enhance `summary` method for `biocrevdeps` class output from
  `pkgBiocRevDeps`

[biocViews](/packages/biocViews)
---------

                       Changes in version 1.77.2                        

ENHANCEMENT

- (1.77.2) Add build_db_from_source, build_meta_aliases_db,
  build_meta_rdxrefs_db

- (1.77.3) Add biocViews term GPU

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.66.0                        

BREAKING CHANGES

- the `https` argument in `useEnsembl()` has been removed, as announced
  in
  version 2.50.0. Connections to Ensembl mirrors now always use
  `https://`, as
  enforced by a redirect from the server.

USER VISIBLE CHANGES

- This package documentation is now available as a pkgdown website at
  https://huber-group-embl.github.io/biomaRt/.

- The `exportFASTA()` function is now roughly 5x faster.

BUG FIXES

- `getGene()` no longer errors with "This function only works when used
  with
  the ensembl BioMart.", that was always returned, even when using the
  ensembl
  BioMart. This bug was introduced in version 2.27.1 (Bioconductor
  3.3).
  Thanks to [Tobias Hoch](https://github.com/toobiwankenobi) for the
  report.

- This package now explicitly depends on R >= 4.1.0. This was
  effectively
  already the case since version 2.58, in line with Bioconductor
  policy, but
  not explicitly documented in `DESCRIPTION`.

- `getBM()` no longer silently convert alleles "T" and "F" to logical
  `TRUE`
  and `FALSE`. This rarely happened when a single allele,
  coincidentally named
  "T" or "F", was returned in a column. Thanks to
  [Sina Rüeger](https://github.com/sinarueeger) for the report.

- `getBM()` now internally handles redirects and automatically
  resubmits the
  query to the new location. This fixes issues when using versioned
  forms such
  as https://eg59-plants.ensembl.org/, which redirects to
  <https://may2024-plants.ensembl.org/index.html>. Thanks to Edoardo
  Bertolini
  for the report.

INTERNAL CHANGES

- Coding style throughout the package has been harmonized using the air
  tool.
  Contributors using RStudio, Positron or VS Code should have their
  code styled
  automatically on save.

- A duplicated definition of the `.getEnsemblSSL()` internal function
  has been
  removed.

- Static analysis via the lintr package is now performed on each push
  and PR.
  It should mostly be invisible to users but might result in slightly
  increased
  performance in some cases.

- This package now uses roxygen2 to generate documentation and
  `NAMESPACE`.

- A `R CMD check` `NOTE` about a missing import has been resolved.

- This package continuous integration setup now errors on `R CMD check`
  `WARNING`s. This is checked on each push and PR.

- The digest and rappdirs dependencies have been removed in favour of
  functions provided by the base R tools package.

- Examples have been reviewed and fixed where necessary.

- testthat edition 3 (instead of previously edition 2) is now used for
  unit
  tests.

[biosigner](/packages/biosigner)
---------

                       Changes in version 1.37.4                        

MINOR MODIFICATION

- minor documentation update

                       Changes in version 1.37.2                        

MINOR MODIFICATION

- SpikePos dataset from the BioMark package now included in biosigner

[Biostrings](/packages/Biostrings)
----------

                       Changes in version 2.78.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The package now depends on Seqinfo instead of GenomeInfoDb.

DEPRECATED AND DEFUNCT

- All the functions that went to the pwalign package in BioC 3.19 and
  got formally deprecated in Biostrings in BioC 3.21 are now defunct.

- After being deprecated in BioC 3.19, XStringPartialMatches objects
  and matchprobes() are now defunct.

- Remove needwunsQS() (was defunct in BioC 3.19).

BUG FIXES

- Fix: minor bug with input verification in `write.phylip`; clean up
  formatting (#124). By Aidan Lakshman.

[blase](/packages/blase)
-----

                        Changes in version 1.0.0                        

- Initial release.

[BreastSubtypeR](/packages/BreastSubtypeR)
--------------

                        Changes in version 1.1.8                        

Highlights

- Paper published in NAR Genomics and Bioinformatics (2025), Editor’s
Choice (DOI: 10.1093/nargab/lqaf131).
- Support for raw RNA-seq counts (requires gene lengths).
- iBreastSubtypeR refresh: cleaner UX, smarter AUTO guidance,
consistent exports.
- Refined, data-driven thresholds for ER/HER2 skew detection in AUTO.
- Broader input-validation across the subtyping pipeline.

Bug Fixes

- AUTO: fixed ER/HER2 strata messaging for ER+/HER2− and ER−/HER2−
cohorts. The spurious message “ssBC.v2 for samples: ERnegHER2neg” no
longer appears in ER+/HER2− cohorts.
- AUTO: in pure HER2+ cohorts, ssBC.v2 is now preferred as intended;
ssBC is not listed or subset-indexed there.
- AUTO: ssBC / ssBC.v2 sample-subsetting and messages now occur only
when that method is actually selected.
- PAM50 variants no longer error with calibration = "None" or
"External".
- ssBC variants handle datasets with <50 PAM50 genes more robustly.
- Fixed data.frame issue (“check.names matched by multiple
arguments”)
via safe builders.
- Removed duplicate “Subtype == BS_*class” columns in downloads.
- Safer ROR merges on PatientID with clearer notifications.
- Eliminated jsonlite named-vector warning in plotting.

Changes

- Documentation: clarified guidance for subtype-specific cohorts:
- ER/HER2-defined cohorts (ER+/HER2−, ER−/HER2−, ER+/HER2+,
ER−/HER2+): NC-based → ssBC.v2 only; plus SSP-based (AIMS,
sspbc).
- ER-only (ER+ or ER−) and TNBC: NC-based → ssBC and/or ssBC.v2
(subject to minimum sizes); plus SSP-based.
- AUTO clarifications (simulation-based defaults):
- ER balance gate: lower_ratio = 0.39, upper_ratio = 0.69.
- Minimum sizes: ER+ total = 15, ER− total = 18, TN total = 18.
- Subgroup gates (rounded): ER+ subgroups (HER2+ / HER2−) use
round(15 / 2); ER− subgroups (HER2+ / HER2−) use round(18 / 2).
- Notes: these thresholds gate method eligibility; they do not
force a consensus call. Details in the vignette (“AUTO logic”).

Shiny (iBreastSubtypeR)

- New hero card + method chips (NC, SSP, ROR, AUTO) and “Why AUTO?”
explainer modal.
- AIMS: 4-class toggle disabled (AIMS is 5-class only).
- AUTO/AIMS/SSPBC: Full-metrics export selector hidden.
- Clear per-method help + PAM50 calibration notes.
- Small UX polish: file dialogs + Mapping preserve scroll position.

Exports

- Two modes: Calls only and Full metrics (incl. ROR when available).
- Standard column names: Call_5class / Call_4class (internal BS /
BS.Subtype mapped).
- AUTO: Calls = selected-k table (+ entropy); Full = selected-k
merged
with other-k (suffix _4class / _5class).
- Filenames encode method, class, and mode (e.g.,
results-PAM50-5class-full.txt).

Behavior & Validation

- ROR only for NC methods when TSIZE/NODE exist and are numeric (auto
coercion + warnings).
- Stricter checks for ssBC subgroups (ER/HER2/TN presence and
coding).
- AUTO preflight warns on missing or mis-coded ER/HER2.

Core Features (unchanged)

- Unified NC (PAM50, cIHC, PCAPAM50, ssBC) and SSP (AIMS, SSPBC)
subtyping under one API.
- BS_Multi with cohort-aware AUTO selection; Mapping for
platform-agnostic preprocessing.
- Interactive Shiny app with Bioconductor-friendly exports.

Upgrade Notes

- Raw RNA-seq counts are supported from v1.1.3 onward (requires gene
lengths).
- If you previously parsed BS / BS.Subtype, switch to Call_5class /
Call_4class.
- Package API unchanged.

[BridgeDbR](/packages/BridgeDbR)
---------

                       Changes in version 2.19.3                        

BUG FIXES

- Sorts the vignette to list the tutorial first

                       Changes in version 2.19.2                        

BUG FIXES

- Better package install instructions in the new secondary vignette

                       Changes in version 2.19.1                        

NEW FEATURES

- Added a compact resource identifier example to the vignette
- Added a second vignette explaining how sec2pri identifier mappings
databases can be used

BUG FIXES

- Handle Bioregistry.io ChEBI compact identifiers better

[BSgenome](/packages/BSgenome)
--------

                       Changes in version 1.78.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The package now depends on the new Seqinfo package instead of
  GenomeInfoDb.

[BSgenomeForge](/packages/BSgenomeForge)
-------------

                       Changes in version 1.10.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add Seqinfo to Depends.

BUG FIXES

- Make sure internal helper create_2bit_BSgenome_datapkg() can
  move the 2bit file across partitions without getting hit by
  an "Invalid cross-device link" error.

[CAGEr](/packages/CAGEr)
-----

                       Changes in version 2.16.0                        

NEW FEATURES

- Pass normalised expression values to quickEnhancers(). Closes #132.
- New correlationMatrix() function to compute correlation matrix
without plotting it. Closes #129.
- Save power law alpha value in the reference.slope element of the
meta slot of the ggplot2 object produced by
plotReverseCumulatives().
- Support for importing CTSS data from bigWig files.

BUG FIXES

- Stop passing StitchedGPos to quickEnhancesr(). Closes #133.
- Migrate from GenomeInfoDb to Seqinfo where needed.

[CalibraCurve](/packages/CalibraCurve)
------------

                       Changes in version 0.99.0                        

Initial submission to BioConductor.

[CAMERA](/packages/CAMERA)
------

                       Changes in version 1.65.1                        

BUG FIXES

- Fix a few more logic expressions causing "the condition has length >
  1" type issues

[Cardinal](/packages/Cardinal)
--------

                 Changes in version 3.11.1 (2025-10-14)                 

NEW FEATURES

- Added 'contrastTest()' method for post-hoc contrasts on 'MeansTest'
  objects

- New 'ContrastTest' class for storing contrast analysis results

- Added support for 'use_lmer=TRUE' in 'meansTest()' to use
  lmerTest::lmer

- Integration with 'emmeans' package for estimated marginal means and
  contrasts

- Support for various contrast methods: pairwise, trt.vs.ctrl, custom
  contrasts

- Added 'emm_adjust' parameter for multiple comparison adjustments

SIGNIFICANT USER-VISIBLE CHANGES

- When 'use_lmer=TRUE' in 'meansTest()', likelihood ratio tests are not
  performed; use 'contrastTest()' for post-hoc comparisons instead

- Suggests 'lmerTest' and 'emmeans' packages for enhanced functionality

[CatsCradle](/packages/CatsCradle)
----------

                        Changes in version 1.3.3                        

- Example data ligandReceptorResults.rda saved in compressed format.

                        Changes in version 1.3.2                        

- Function getLigandReceptorPairsInPanel rewritten for faster
runtime.
- Analytical option to calculate pvalues added to
computeNeighbourEnrichment. This has a faster runtime.
- Function performLigandReceptorAnalysis rewritten for faster runtime
by making use of sparse matrices.
- Analytical option to calculate pvalues added to
performLigandReceptorAnalysis. This has a faster runtime.
- Function plotLRDotplot added to visualise ligand-receptor
interactions returned by performLigandReceptorAnalysis.

                        Changes in version 1.3.1                        

- smallXenium@images$fov@misc = list() field added to smallXenium
example data for compatability with Seurat updates.

[cBioPortalData](/packages/cBioPortalData)
--------------

                       Changes in version 2.22.0                        

New features

- Defunct removeDataCache for cBioPortalData as caching mechanism is
no longer used. Run cBioCache() to get the cache location for
deletion. Note that cBioDataPack will continue to cache data files
at that location.

[cellmig](/packages/cellmig)
-------

                       Changes in version 0.99.14                       

- Initial submission to Bioconductor

[cellxgenedp](/packages/cellxgenedp)
-----------

                        Changes in version 1.13                         

BUG FIXES

- (v. 1.12.1 / 1.13.1) Do not attempt to 'unlist' any facet column in
datasets or collections

[ceRNAnetsim](/packages/ceRNAnetsim)
-----------

                 Changes in version 1.21.0 (2025-08-10)                 

- ggplot2 uses a new class of the objects, do some tests were deprecated.  
- The citation errors were fixed.

[chevreulPlot](/packages/chevreulPlot)
------------

- Initial Bioconductor submission

[chevreulProcess](/packages/chevreulProcess)
---------------

- Initial Bioconductor submission

[chevreulShiny](/packages/chevreulShiny)
-------------

- Initial Bioconductor submission

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.43.2                        

- Fix an issue for the corrdinates changes of GRCh37 annotation.

                       Changes in version 3.43.1                        

- Fix an issue in genomicElementDistribution when convert
  seqlevelsStyle.

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.45.2                        

- new cache mechanism from 'yulab.utils' (2025-10-15, Wed)

[Chromatograms](/packages/Chromatograms)
-------------

- Initial Bioconductor submission

[CircSeqAlignTk](/packages/CircSeqAlignTk)
--------------

                 Changes in version 1.11.2 (2025-05-09)                 

SIGNIFICANT USER-VISIBLE CHANGES

- Remove `get_slot_contents` function. Use `slot` function instead.

- Change default aligner to HISAT2 from Bowtie2.

- Change repository.

- Update documentation.

[ClassifyR](/packages/ClassifyR)
---------

                       Changes in version 3.14.0                        

- 
  Precision pathways gains an all-combinations-of-all-subsets
  mode.

[cleaver](/packages/cleaver)
-------

                 Changes in version 1.47.1 (2025-10-14)                 

- Move package to codeberg.org.

[ClustIRR](/packages/ClustIRR)
--------

                        Changes in version 1.7.7                        

- Bugfix, model improvement

                        Changes in version 1.7.3                        

- Function for quantifying community purity wrt node features

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

                       Changes in version 2.25.2                        

- `make_layout,HeatmapList-class`: `row_order`/`column_order` is reset
  if `cluster_rows`/`cluster_columns`
  is set to TRUE in `draw()`.

                       Changes in version 2.25.1                        

- `adjust_heatmap_list()`: fixed a bug where slice gaps were repeatedly
  used to calculate the heatmap width/height. #1240

- `names<-,HeatmapAnnotation()`: `label` slot is correctedly assigned
  with the new name.

[consensusSeekeR](/packages/consensusSeekeR)
---------------

                       Changes in version 1.37.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- The package documentation has been updated.

[Coralysis](/packages/Coralysis)
---------

                        Changes in version 0.99.0                        

First release:

- Functions:
- AggregateDataByBatch()
- RunParallelDivisiveICP()
- ReferenceMapping()
- GetCellClusterProbability()
- SummariseCellClusterProbability()
- CellClusterProbabilityDistribution()
- BinCellClusterProbability()
- TabulateCellBinsByGroup()
- CellBinsFeatureCorrelation()
- FindClusterMarkers()
- FindAllClusterMarkers()
- MajorityVotingFeatures()
- RunPCA()
- RunTSNE()
- RunUMAP()
- PCAElbowPlot()
- PlotDimRed()
- PlotExpression()
- PlotClusterTree()
- VlnPlot()
- HeatmapFeatures()
- PrepareData()
- GetFeatureCoefficients()
- Vignettes:
- Integration
- Reference-mapping
- Cell States

[COTAN](/packages/COTAN)
-----

                        Changes in version 2.9.7                        

Fixed underflow error appearing when the dispersion is very small

Improved the test.dataset to be somewhat more realistic

                        Changes in version 2.9.6                        

Upgraded dispersion solver to use Newton-Raphson method Minor
variations
of the model parameters are expected

Improved the torch coex code so to use less memory and do less
calculations

                        Changes in version 2.9.5                        

Reduced significantly the memory foot-print while use multi-process
solvers Also improved speed of GDI and data-reduction calculations,
by
going multi-process there too

Improved function cellsUniformClustering() mechanism in case of
larger
datasets when resolution ceiling was causing unnecessary slowdown and
in
cases early termination too

                        Changes in version 2.9.4                        

Solved minor issue with maxIterations argument in the function
cellsUniformClustering(): it was doing two extra loops

Fixed bug with function isCoexAvailable() when used with COTAN
objects
created with older versions

                        Changes in version 2.9.3                        

Fixed issue with function getSelectedGenes() in cases when asked for
HGDI selection and dataset GDI is low overall: it was returning too
few
and, in some cases, even no genes

Added possibility to choose which data matrix to use for
clusterizations
and UMAP plots. Also now it is possible to use a projection onto the
main COEX matrix eigen-vectors as alternative to the PCA as
dimensionality reduction tool

                        Changes in version 2.9.2                        

Solved issues with breaking changes from zeallot package

Expanded calculateLikelihoodOfObserved() to support more formulas
related to likelihood

Unified code that produces UMAP plots: now they all run UMAPPlot()
that
in turn calls the Seurat::RunUMAP() function

Fixed issue with use of parallelDist::parDist() on Linux aarch64
machines: it crashed irrecoverably due to an invalid alignment error.
Now code calls the normal dist() function on such machines

                        Changes in version 2.9.1                        

Solved issue with reordering of merged clusterization: flag
keepMinusOne
was off and label "-1" was moved around

Added checker class type to the data.frame obtained after conversion.
It
will make loading checker results from file easier

Improved README.md with stable links to devel and release versions in
BioConductor

[countsimQC](/packages/countsimQC)
----------

                       Changes in version 1.27.1                        

- Minor maintenance updates

[CPSM](/packages/CPSM)
----

                        Changes in version 1.1.4                        

Updates

- Updated Lasso_PI_scores_f.R, MTLR_pred_model_f.R,
Univariate_sig_features_f.R, km_overlay_plot_f.R,
mean_median_surv_barplot_f.R, surv_curve_plots_f.R,
predict_survival_risk_group_f.R to provide customization for user in
cross-validation , Font-size and figure styling etc.
- Updated MTLR_pred_model_f.R to compute IBS (integrated brier score)
- Updated vignettes/CPSM.Rmd to provide improve description and
workflow
- Updated tests/testthat files to incoprate updated changes
- Updated DESCRIPTION file:
- Bumped version from 1.1.3 to 1.1.4.

                        Changes in version 1.1.3                        

Updates

- Added CC0 license for processed data in the /data folder to comply
with data-sharing requirements.
- Updated km_overlay_plot_f.R to improve figure
- Updated DESCRIPTION file:
- Added contributors (ctb) to Authors@R.
- Bumped version from 1.1.2 to 1.1.3.

                        Changes in version 1.1.2                        

##Improvement

- Updated MTLR_pred_model_f.R, Lasso_PI_scores_f.R,
mean_median_surv_barplot_f.R,
Nomogram_generate_f.R,surv_curve_plots_f.R script to remove
undesired package
- Updated km_overlay_plot_f.R to improve figure
- Updated DESCRIPTION file include required biocview terms

Version update

- Bumped version from 1.1.1 to 1.1.2.

                        Changes in version 1.1.1                        

##Improvement

- Updated MTLR_pred_model_f function to make it simple

Version update

- Bumped version from 1.1.0 to 1.1.1.

                        Changes in version 1.1.0                        

New Features

- predict_survival_risk_group_f():
- Added a new function to build Random Forest-based survival risk
group classifiers using selected molecular features.
- Automatically selects the best-performing model based on
prediction accuracy across various ntree values (10 to 1000).
- Outputs include predicted risk groups, prediction probabilities,
misclassification summary, best model, and tree-wise performance
matrix.
- km_overlay_plot_f():
- Introduced a new function to generate Kaplan-Meier survival
plots overlayed with individual test sample survival curves.
- Useful for visualizing predicted test sample risk in the context
of population-level survival groups.
- Includes annotated output with sample ID, predicted group, and
probability.

Updated data files, added new files

Improvement

- Updated MTLR_pred_model_f function to compute MAE and improve the
structure of function.

Documentation

- Updated Description file
- Updated vignette file CPSM.Rmd:
- Included usage details and examples for new functions.
- Enhanced step-by-step explanation of prediction and
visualization workflows.

Version update

- Bumped version from 1.0.0 to 1.1.0

[crisprBase](/packages/crisprBase)
----------

                       Changes in version 1.13.2                        

- Added PAM weights to the SPG nuclease.

                       Changes in version 1.13.1                        

- Added the option of not scaling editing weights in BaseEditor
  function.

[crumblr](/packages/crumblr)
-------

                        Changes in version 1.0.7                        

- Feb 1, 2022
- fix vst()
- add plots plotScatterDensity() and meanSdPlot() for VST
- add vignette vst.Rmd

                        Changes in version 1.0.6                        

- Jan 31, 2022
- add logo, example data

                        Changes in version 1.0.5                        

- Dec 27, 2021
- add dmn.mle()
- use generic crumblr() function

                        Changes in version 1.0.4                        

- Dec 23, 2021
- Update vignette
- export clr()

                        Changes in version 1.0.3                        

- Dec 20, 2021
- fix vst bug

                        Changes in version 1.0.2                        

- Dec 20, 2021
- add vst

                        Changes in version 1.0.1                        

- Oct 29, 2021
- handle data.frame

[CSOA](/packages/CSOA)
----
                 Changes in version 0.99.0 (2025-07-12)                 

- Submitted to Bioconductor

[CTexploreR](/packages/CTexploreR)
----------

                         Changes in version 1.5                         

CTexploreR 1.5.2

- Remove documentation for internal functions.

CTexploreR 1.5.1

- Move some tests to longtests.

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

                         Changes in version 1.9                         

CytoPipeline 1.9.4

- improved aggregateAndSample() added a setup parameter, such that
when it is equal to "forceNEvent", the nb of events selected per
flow frame is not always balanced, to always obtain nTotalEvents
(except when the total nb of available events, taking all
flowFrames, is less than nTotalEvents).

CytoPipeline 1.9.3

- added sampleDisplayNames()

CytoPipeline 1.9.2

- sample files can now have duplicate base names (provided full paths
are different)
- pData<- is now more liberal.

1.  It can accept new pData containing more rows than existing sample
names (the corresponding subset of pData is taken).
2.  It can accept pData with row names pointing to either sample file
full paths or base file names
3.  It can accept pData with no row names provided the number of rows
correspond to the number of sample files. Row names are then set by
default to sample file base names (if unique), or sample file full
paths.

CytoPipeline 1.9.1

- upgraded to GHA cache v4

[CytoPipelineGUI](/packages/CytoPipelineGUI)
---------------

                         Changes in version 1.7                         

CytoPipelineGUI 1.7.2

- using sample display names instead of base names or full paths, as
provided from CytoPipeline 1.9.3

CytoPipelineGUI 1.7.1

- migrated to GHA cache v4

[dar](/packages/dar)
---

                        Changes in version 1.5.3                        

Bug Fixes

- https://github.com/MicrobialGenomics-IrsicaixaOrg/dar/issues/117

[DeconvoBuddies](/packages/DeconvoBuddies)
--------------

                        Changes in version 1.1.5                        

BUG FIXES

Correct (est_prop/est_prop_test) data usage after changes in v1.1.4

                        Changes in version 1.1.4                        

BUG FIXES

Package BisqueRNA is no longer a suggested dependency (no longer
available on CRAN).

Vignette "Deconvolution Benchmark in Human DLPFC" now loads
pre-computed
est_prop data instead of running Bisque deconvolution.

                        Changes in version 1.1.3                        

NEW FEATURES

- plot_gene_express() now has the option to use "free_y" axis.

                        Changes in version 1.1.2                        

NEW FEATURES

- plot_gene_express() now has the option to change plot_type to
'violin' or 'boxplot'

                        Changes in version 1.1.1                        

BUG FIXES

- get_mean_ratio() was initially coercing sparse matrices into in
memory matrices. We resolved this by using
MatrixGenerics::rowMeans() and MatrixGenerics::rowMedians(). This
issue was reported by @cyntsc.

[DeeDeeExperiment](/packages/DeeDeeExperiment)
----------------

- Initial Bioconductor Submission

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

                       Changes in version 1.45.1                        

- Fixed labs() function calls in multiple R files:
  - methods.R (2 occurrences in degMean and degVar functions)
  - clustering.R (2 occurrences)
  - volcanoPlot.R (1 occurrence)

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.36.0                        

DEPRECATED AND DEFUNCT

- Get rid of SparseArraySeed objects, OLD_extract_sparse_array(), and
  read_sparse_block(). These were deprecated in BioC 3.20 and defunct
  in BioC 3.21.

[DEP](/packages/DEP)
---

                       Changes in version 1.31.1                        

- 
  Fixed broken documentation links to Msnbase impute methods.

- 
  Fixed crash at the launch of the shiny apps.

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.50.0                        

- New method UPSHOT for lfcThreshold/altHypothesis.

- Fixes for lesser-used functionality, e.g. useT, weights
  with numeric contrast and the SE computation, summary
  printing.
  CHANGES IN VERSION 1.49.6

- Fixes for combined use of useT and greaterAbs in results(),
  thanks to a PR from Raphael (@rrross)

                       Changes in version 1.49.4                        

- The wrong weights matrix was being used when recomputing the
  SE within results() for the numeric-style contrast. Fixed.

- glmGamPoi and parallel error instead of warning.

                       Changes in version 1.49.1                        

- Resolving an issue where lfcThreshold was changed to 0 in
  the summary() print-out after lfcShrink() even if s-values
  were not computed.

[DFplyr](/packages/DFplyr)
------

                        Changes in version 1.4.0                        

NEW FEATURES

- [ and rbind methods now respect groups
- Added left_join(), right_join(), inner_join(), and full_join()

SIGNIFICANT USER-VISIBLE CHANGES

- rename2 deprecated; rename is now a DataFrame method on S4Vectors
generic

[dittoSeq](/packages/dittoSeq)
--------

                       Changes in version 1.21.1                        

- Minor changes to ensure compatibility with ggplot2 version 4.0.0
changes. Notably, the 'MASS' package is now added as a suggested
companion package, and error messages suggesting its installation
were added for features which rely on functionality from this
package. This change is concurrent with ggplot2 itself downgrading
the tool from an imported requirement to a suggested addition.
- Minor changes to ensure compatibility with ComplexHeatmap in how
the
'drop_levels' input of 'dittoHeatmap()' is handled. See function
documentation for details.

[dmGsea](/packages/dmGsea)
------

                        Changes in version 0.99                         

- Initial release to Bioconductor.

[DMRcaller](/packages/DMRcaller)
---------

Changes in Version 1.42.1:

- Added functions to read Oxford Nanopore methylation data from Dorado bam files
- Added functions to detect Partially Methylated Domains PMDs
- Added functions to detect Variably methylated domains (VMDs) using ONT data
- Added functions to detect Variably Methylated Regions (VMRs) using ONT data
- Added functions to detect co-methylation using ONT data
- Updated parallel computing to use BiocParallel

[DOtools](/packages/DOtools)
-------

- Initial submission to Bioconductor

[DropletUtils](/packages/DropletUtils)
------------

                       Changes in version 1.30.0                        

- Further improvements to the knee point detection in
  barcodeRanks(). This now uses a curve-tracing algorithm where
  the knee is defined inside the window with the highest
  curvature. We also refine inflection point to report the
  earlier location if multiple inflection points are present on
  the curve.

- Improve the efficiency of Matrix Market file parsing in
  read10xCounts(). We support multi-threaded parsing and can
  directly read into dgCMatrices or SVT_SparseMatrices without
  intermediate. This is based on the eminem C++ library from the
  assorthead package.

[DspikeIn](/packages/DspikeIn)
--------

- Initial submission to Bioconductor

[DuplexDiscovereR](/packages/DuplexDiscovereR)
----------------

                        Changes in version 1.3.4                        

- Date: 2025-09-29
- Added experimental functionality for prioritisation of ambiguous
annotation to annotateGi()

                        Changes in version 1.3.1                        

- Date: 2025-07-04
- Added trimming of the PE alignments relative to the ligation point.
- For the PE reads processed by STAR, the type of the junction is
detected.
- To narrow down the hybrid prediction, user can trim reads -
N_trim nt will be left adjacent to the ligation point.
- Added missing arguments to the function runDuplexDiscoverer which
runs the whole pipeline. Keys for annotation features, trimming and
minimum chimeric length cutoffs have defaults, also can be set by
user to be passed to the functions called from inside the
runDuplexDiscoverer.

[dupRadar](/packages/dupRadar)
--------

                       Changes in version 1.39.1                        

- Feature: add tmpDir parameter to analyzeDuprates (credit: Adam
  Talbot)
  Resolves performance issues on FUSE-based filesystems by allowing
  temporary files to be written to fast local storage instead of CWD

[edgeR](/packages/edgeR)
-----

                 Changes in version 4.8.0 (2025-10-01)                  

- 
  The `adaptive.span` argument of voomLmFit() is now set to
  `TRUE` by default, meaning that the span is chosen by
  chooseLowessSpan().

- 
  The default value of `span` for estimateDisp() and WLEB() is
  now chosen by chooseLowessSpan().

- 
  New function catchOarfish().

- 
  New argument `path` for catchRSEM().

- 
  glmTreat() no longer writes a message about switching to
  glmLRT() when `lfc=0`.

- 
  Revise help page for diffSplice.DGEGLM().

- 
  Revise help page for filterByExpr().

[eisaR](/packages/eisaR)
-----

                       Changes in version 1.21.1                        

- various code improvements detected by lintr

[ELViS](/packages/ELViS)
-----

                       Changes in version 1.1.12                        

- Minor bug fix in test script

                       Changes in version 1.1.11                        

- Added error handling code regarding conda env creation error

                       Changes in version 1.1.10                        

- Reflected changes made in 1.1.9 to test scripts
- If conda env creation fails in samtools_reticulate mode, it
defaults
to Rsamtools for base-level read depth calculation.

                        Changes in version 1.1.9                        

- Using reticulate to install and use conda environment

                        Changes in version 1.1.8                        

- Fixed some incompatibility issue with basilisk package. basilisk
does not support conda package installation at least for now. So
partly reverted to reticulate package for samtools installation.

                     Changes in version 1.1.1-1.1.7                     

- Build error fix

                        Changes in version 1.1.0                        

- Released in Bioconductor

[EnrichmentBrowser](/packages/EnrichmentBrowser)
-----------------

                       Changes in version 2.40.0                        

- getGenesets now supports MSigDB's new designated mouse collections
  named
  'MH' (hallmark), 'M1'(genomic position), 'M2' (curated databases),
  and so on.

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.29.4                        

- remove deprecated aes_string/aes_ (2025-10-23, Thu, #332)

                       Changes in version 1.29.3                        

- bug fixed of cnetplot for CompareClusterResult (2025-09-13, Sat,
#329)
- color gene according to the gene cluster info
- bug fixed in pie scale label (2025-07-14, Mon, #328)

                       Changes in version 1.29.2                        

- update treeplot with two parameters, leave_fontsize and
clade_fontsize (2025-07-12, Sat, #324, #325)
- remove the fontsize parameter as it only works for
clade_fontsize
- 'log10' transform for pvalue color scale by default (2025-07-12,
Sat, #316)
- introduce new parameters in gseaplot2() (2025-07-12, Sat)
- pvalue_table_columns
- pvalue_table_rownames
- https://github.com/YuLab-SMU/clusterProfiler/issues/774

                       Changes in version 1.29.1                        

- throw error in goplot() if ontology is not one of the 'MF', 'CC' or
'BF' (2025-04-28, Mon, clusterProfiler#768)

[enrichViewNet](/packages/enrichViewNet)
-------------

                        Changes in version 1.7.3                        

NEW FEATURES

- The new functions createBasicEmapAsIgraph(),
  createEnrichMapMultiBasicAsIgraph(), generates enrichment maps in an
  igraph format.

                        Changes in version 1.7.2                        

BUG FIXES AND IMPROVEMENTS

- All functions using ggplot2 are using updated ggplot2 functions.

                        Changes in version 1.7.1                        

BUG FIXES

- The package uses a fixed version of gprofiler database in the unit
  tests.

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.33.2                        

- Update description to create `EnsDb` databases using the docker
  image.

[EpiCompare](/packages/EpiCompare)
----------

                       Changes in version 1.13.2                        

Bug Fixes

- [BREAKING CHANGE] Temporarily disable upset_plot option as one of
the dependencies, ComplexUpset, is not yet updated to support
ggplot 4.4.0.
- Fix deprecated ggplot arguments.

Miscellaneous

- add_download_button:
- No longer changes to TRUE if run_all is TRUE.

[EpipwR](/packages/EpipwR)
------

                        Changes in version 1.3.2                        

Major updates

- Users now have the ability to provide custom methylation data
- Users now have the ability to specify a non-normal distribution for
phenotype data
- A new function is available for pwrEWAS users to "translate" a
pwrEWAS distribution of effect sizes into EpipwR parameters. Please
note that all pwrEWAS effect size distributions are truncated normal
distributions centered at 0

[epiregulon](/packages/epiregulon)
----------

                        Changes in version 2.0.0                        

- New function in the workflow, optimizeMatacellNumber, which
estimates the optimal value of the cellNum parameter passed to
calculateP2G.
- Using scrapper::clusterKmeans to implement deterministic algorithm
of determination k-mean cluster. It is used to create metacells in
the calculateP2G function.
- Added to the output of calculateP2G function empirical p-values and
FDR. The null distribution is calculated based on random links of
peaks to the genes from other chromosomes.
- Default ChIP-seq output from tfBinding is changed to version 2
which
imposes more stringent cutoffs. We have increased the number of
unique reads to 20M (previously 10M in version 1) and the number of
peaks passing FDR < 1e-5 to 1000 peaks (previously 100 peaks in
version 1). For ENCODE data, we now remove any samples with
Audit.NOT_COMPLIANT or Audit.ERROR. Flags.
- The argument method to calculateActivity has been deprecated.
- Bug corrected in addLogFC. p-values were mistakenly transformed
from
log10 when it should be transformed from natural log when logFC_ref
or/and logFC_condition was specified

                        Changes in version 1.5.1                        

- checking for the duplicated gene names in the input gene expression
SingleCellExperiment
- validation of the regulon object passed to pruneRegulon,
addWeights,
and calculateActivity

[epistasisGA](/packages/epistasisGA)
-----------

                       Changes in version 1.11.2                        

- Add reference to maternal/fetal GADGETS paper.
- Update E-GADGETS and Vignette

[escape](/packages/escape)
------

                 Changes in version 2.5.5 (2025-06-11)                  

Bug fix & enhanced functionality

- Enable color.by for both metadata columns and features (other gene
sets)
- Introduce summarise.by argument for geyserEnrichment()
- Enable scaling if color.by is another gene.set. Enable scaling for
dgCMatrix

                 Changes in version 2.5.4 (2025-06-05)                  

Bug fixes

- Fix conversion wide-to-long format for heatmapEnrichment()
- Fix issue with t() call on sparse matrices.

                 Changes in version 2.5.3 (2025-05-19)                  

Highlights

- Streamlined code-base – major internal refactor for clarity, speed
and a ~20 % smaller dependency tree.
- Consistent, flexible visualisation API across all plotting helpers.
- Robust unit-test suite (>250 expectations) now ships with the
package.

New & enhanced functionality

| Area | Function(s) | What changed |
|------|-------------|--------------|
| **Visualisation** | `ridgeEnrichment()` | *True gradient* coloring mode for numeric `color.by`; optional per-cell rugs; quantile median line; fixed grey-fill bug |
| | `densityEnrichment()` | accepts new `rug.height`; ~4× faster ranking routine using `MatrixGenerics::rowMeans2`; cleaner two-panel layout via **patchwork** |
| | `gseaEnrichment()` | new `rug.height`; clearer legend showing ES/NES/ *p*; internal vectorised ES calculation |
| | `splitEnrichment()` | rewritten: split violins when `split.by` has 2 levels, dodged violins otherwise; inline boxplots; auto Z-scaling; palette helper |
| | `scatterEnrichment()` | density-aware points (via **ggpointdensity**), hex-bin alternative, optional Pearson/Spearman overlay, continuous or discrete color mapping |
| **Dimensionality reduction** | `performPCA()` / `pcaEnrichment()` | uses `irlba::prcomp_irlba()` automatically for large matrices; stores eigen-values/contribution in `misc`; `add.percent.contribution` now always respected |
| **Scoring backend** | `escape.matrix()` / `.compute_enrichment()` | lazy loading of heavy back-ends (*GSVA*, *UCell*, *AUCell*); unified `.build_gsva_param()`; drops empty gene-sets up-front |
| **Normalization** | `performNormalization()` | chunk-wise expressed-gene scaling (memory-friendly); accepts external `scale.factor`; optional signed log-transform; returns object with assay `<assay>_normalized` |
| **Gene-set retrieval** | `getGeneSets()` | downloads now cached under `tools::R_user_dir("escape", "cache")`; graceful KEGG append; clearer error for non-human/mouse requests |

Performance & dependency reductions

- Replaced plyr, stringr, rlang usage with base-R helpers; these
packages are now Suggests only.
- Common color and label utilities (.colorizer(), .colorby(),
.orderFunction()) removed redundant tidyverse imports.
- Internal matrices split/chunked with new .split_* helpers to cap
memory during parallel scoring/normalization.

Bug fixes

- Gradient mode in ridgeEnrichment() no longer produces grey fills
when the chosen gene-set is mapped to color.by.
- pcaEnrichment() axis labels correctly include variance contribution
when display.factors = FALSE.
- .grabDimRed() handles both Seurat v5 and <v5 slot structures; fixes
missing eigen-values for SCE objects.
- escape.matrix() respects min.size = NULL (no filtering) and handles
zero-overlap gene-sets gracefully.
- Global variable declarations consolidated – eliminates R CMD check
NOTES regarding na.omit, value, etc.

Documentation

- DESCRIPTION rewritten – heavy packages moved to Suggests; added
explicit Config/reticulate for BiocParallel.
- escape.gene.sets data object now fully documented with source,
usage, and reference.

                 Changes in version 2.4.1 (2025-03-05)                  

- Version bump to align with Bioconductor release cycle.
- escape.matrix() now silently removes gene-sets with zero detected
features.

[EWCE](/packages/EWCE)
----

                       Changes in version 1.17.1                        

- Update maintainer (Alan -> Hiru)
- Fix deprecated ggplot arguments

[ExperimentHub](/packages/ExperimentHub)
-------------

                       Changes in version 2.99.0                        

MAJOR UPDATES

- Update downloading method from httr to httr2

BUG FIXES

- (2.99.1) Add ability to pass config to download resources

- (2.99.2) Add ability to pass progress to download resources

- (2.99.4) Fix call to BiocManager::install. use `update=FALSE` instead
  of
  `suppressUpdates=TRUE`

- (2.99.6) Update installed.packages to requireNamespace

[extraChIPs](/packages/extraChIPs)
----------

                       Changes in version 1.12.1                        

Bug Fixes

- Patched for ggplot2 v4.0.0
- Switched to SimpleUpset for Upset plots

[faers](/packages/faers)
-----

                        Changes in version 1.5.1                        

- Fixed error in handle_setopt(h, ...) caused by unsupported option
multi_timeout. Requires curl >= 6.0.0.

[fenr](/packages/fenr)
----

                        Changes in version 1.8.2                        

- Stopped using cache in tests, as it might lead to failures in case
of a corrupted cache.


[fgsea](/packages/fgsea)
-----

                       Changes in version 1.35.7                        

- fix handling of zero stat values with gseaParam=0 in plotEnrichment

                       Changes in version 1.35.3                        

- plotGesecaTable shows a legend

- plotCoregulationProfileImage for in-situ profiling

- change default color scheme for GESECA plots

- fix GESECA plots for Seurat 5

                       Changes in version 1.35.1                        

- more reproducibility between platforms (fix #170), but the exact
  results differ from previous versions.

[FinfoMDS](/packages/FinfoMDS)
--------

                       Changes in version 0.99.2                        

- Initial Bioconductor submission.

[FLAMES](/packages/FLAMES)
------

                 Changes in version 2.3.3 (2025-06-18)                  

- Gzipped FASTQ outputs to reduce file size

- Support adding minimap2 indexed genome (as `genome_mmi`)

- Added crew controllers support, steps can now be submitted to slurm
  etc.

                 Changes in version 2.3.1 (2025-05-22)                  

- Refactored pipelines to S4 classes

- Included minimap2, oarfish binaries

- Added a permutation based DTU method in sc_DTU_analysis

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

[FRASER](/packages/FRASER)
------

                        Changes in version 2.4.6                        

- Fixing a bug in the annotation, where we miss the annotation for
  junction,
  where we do not have any splice site available at the time of
  annotation.
  Thanks to @AtaJadidAhari

                        Changes in version 2.4.5                        

- We fixed a bug affecting unstranded paired-end data. We were counting
  read
  fragments instead of read pairs. The counted split reads with the
  fixed
  method are, on average, lower by around 20%. If you have such data,
  please
  rerun your entire analysis from the start (make sure you have
  recount=TRUE).

                        Changes in version 2.4.4                        

- Bugfix: reading fragments instead of read pairs for unstranded
  paired-end data.
  If you have such kind of data, please rerun the whole analysis from
  the start (make sure you have recount=TRUE).

                        Changes in version 2.4.3                        

- Fix bug in cache

                        Changes in version 2.4.2                        

- Add Optimal Hard Thresholding to determine optimal encoding dimension
  #69

- Add singular value plot for optimal encoding dimension

- Add updateSeqlevelsStyle to change seqlevels without internet
  connection

[G4SNVHunter](/packages/G4SNVHunter)
-----------

                        Changes in version 1.1.4                        

- Attempting to resolve build issues introduced by the uniquely
genius
Taishan compilation environment

                        Changes in version 1.1.3                        

- Remove .o files

                        Changes in version 1.1.2                        

- Tweaked the vignette

                        Changes in version 1.1.1                        

- Deprecated functions (still supported but scheduled for removal in
a
future release):

- checkSNV()
- SNVImpactG4() | replaced by G4VarImpact()
- filterSNVImpact() | replaced by filterVarImpact()
- plotSNVImpact() | replaced by plotVarImpact()
- plotImpactSeq() | replaced by plotImpactedG4()

- Updated functions:

- G4HunterDetect() | now includes global metadata in the returned
object

- Removed functions:

- validateG4HunterParams() | replaced by
validateG4HunterDetectInputs()

- New functions:

- exportG4()
- exportMutG4()
- filterVarImpact()
- G4VarImpact()
- loadVariant()
- plotVarImpact()
- plotImpactedG4()

                        Changes in version 1.1.0                        

- Bioconductor auto-bumped devel version to 1.1.0 without code
changes.

[GCPtools](/packages/GCPtools)
--------

                        Changes in version 1.0.0                        

- Initial Bioconductor release!

[gDR](/packages/gDR)
---

               Changes in version 2025-10-30 (2025-10-30)               

- synchronize Bioconductor and GitHub versioning

               Changes in version 2025-06-11 (2025-06-11)               

- update installation instructions

 
[gDRcore](/packages/gDRcore)
-------

               Changes in version 2025-10-30 (2025-10-30)               

- synchronize Bioconductor and GitHub versioning

               Changes in version 2025-08-12 (2025-08-12)               

- fix usage of ifelse

               Changes in version 2025-07-28 (2025-07-28)               

- replace row/col_fittings in source for matrix Metrics assay

               Changes in version 2025-07-15 (2025-07-15)               

- swap drugs to unify combination order

               Changes in version 2025-07-14 (2025-07-14)               

- update headers in fit_SE

               Changes in version 2025-06-10 (2025-06-10)               

- switch from lapply into gDRutils::loop in create_SE

               Changes in version 2025-05-22 (2025-05-22)               

- add functions for reannotating SummarizedExperiment and
MultiAssayExperiment objects
- update vignettes

               Changes in version 2025-05-06 (2025-05-06)               

- remove support for masked values

[gDRimport](/packages/gDRimport)
---------

               Changes in version 2025-10-30 (2025-10-30)               

- synchronize Bioconductor and GitHub versioning

               Changes in version 2025-08-12 (2025-08-12)               

- fix usage of ifelse

               Changes in version 2025-07-17 (2025-07-17)               

- remove run suffix in the PRISM data

               Changes in version 2025-07-10 (2025-07-10)               

- add support for new format of private PRISM data

               Changes in version 2025-05-15 (2025-05-15)               

- add support for depmap public PRISM data

[gDRstyle](/packages/gDRstyle)
--------

               Changes in version 2025-10-30 (2025-10-30)               

- synchronize Bioconductor and GitHub versioning

               Changes in version 2025-08-12 (2025-08-12)               

- fix usage of ifelse

               Changes in version 2025-08-04 (2025-08-04)               

- added support for skipping vignette building in run_tests.sh

               Changes in version 2025-05-26 (2025-05-26)               

- add check for pkgdown in the checkPackage function

[gDRutils](/packages/gDRutils)
--------

               Changes in version 2025-10-30 (2025-10-30)               

- synchronize Bioconductor and GitHub versioning

               Changes in version 2025-09-30 (2025-09-30)               

- add predict_smooth_from_combo function

               Changes in version 2025-08-12 (2025-08-12)               

- fix usage of ifelse

               Changes in version 2025-08-05 (2025-08-05)               

- fix convert_se_assay_to_dt to display appropriate columns
- add flag for capping values in convert_se_assay_to_custom_dt

               Changes in version 2025-08-01 (2025-08-01)               

- improve logic in split_SE_components

               Changes in version 2025-07-28 (2025-07-28)               

- replace row/col_fittings in source for matrix Metrics assay

               Changes in version 2025-07-14 (2025-07-14)               

- add "maxlog10Concentration" and "N_conc" to headers list

               Changes in version 2025-07-01 (2025-07-01)               

- refactor average_biological_replicates_dt to properly average data

               Changes in version 2025-06-24 (2025-06-24)               

- add support in merge_MAE for merging mixed experiments

               Changes in version 2025-06-09 (2025-06-09)               

- update loop to support batch approach

               Changes in version 2025-06-04 (2025-06-04)               

- update dependencies

               Changes in version 2025-05-26 (2025-05-26)               

- update documentation

               Changes in version 2025-05-21 (2025-05-21)               

- add support for combo data in merge_SE
- add merge_MAE fun

               Changes in version 2025-05-13 (2025-05-13)               

- add append arg to .set_SE_metadata

               Changes in version 2025-05-12 (2025-05-12)               

- improve logic in cap_assay_infinities

               Changes in version 2025-05-07 (2025-05-07)               

- update cap_assay_infinities
- add additional grouping columns in cap_assay_infinities

               Changes in version 2025-04-28 (2025-04-28)               

- add result entries in headers_list for combination data

               Changes in version 2025-04-17 (2025-04-17)               

- refactor get_gDR_session_info to display used packages

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.46.0                        

UTILITIES

- fix a strange bug in `append.gdsn(, val=character())` in MacOS when
  compiled by Apple clang version 17


[gemma.R](/packages/gemma.R)
-------

                        Changes in version 4.0.0                        

- get_annotation_children and get_annotation_parents are added.
get_child_terms is deprecated in favor of get_annotation_children

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 3.4.0                        

Other notes

- The enrichment_map() function now also stores as attribute the
original values of the property used to encode the color for the
nodes. This change has been made to simplify the transfer of
information to any gg-based visual representation where the
value-to-color encoding is best left handled in an automated manner
(i.e. without enforcing a color-centric attribute, as it is needed
within visNetwork)
- The message content returned by the shaker_ functions is not
correctly matching the tools they were referring to, thanks
@NajlaAbassi for the fix!
- The ggs_backbone() function has been adapted to reflect the latest
changes in version >= 3.0.0 of the backbone package. Thanks to
@zpneal for the heads up and to @edo98811 for spotting this out
(fixes https://github.com/federicomarini/GeneTonic/issues/67)

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.46.0                        

NEW FEATURES

- Started to include a small "chrominfo db" in the package with
  chromosome
  info for the most commonly used NCBI and UCSC genome assemblies e.g.
  GRCh38.p14, GRCh38.p13, GRCm39, hg38, mm39, etc...
  This makes calls like getChromInfoFromNCBI("GRCh38.p14") or
  getChromInfoFromUCSC("hg38") work offline, plus now they are fast
  and reliable.

SIGNIFICANT USER-VISIBLE CHANGES

- The Seqinfo class was moved from the GenomeInfoDb package to the
  new Seqinfo package.
  IMPORTANT WARNING: This means that:
  - Seqinfo-holding S4 objects that were serialized before this change
  (i.e. prior to Bioconductor 3.22) are no longer valid objects in
  Bioconductor >= 3.22. However they can be fixed with updateObject().
  - This also means that Seqinfo-holding S4 objects serialized with
  Bioconductor >= 3.22 are no longer valid objects in Bioconductor
  < 3.22. These objects can **not** be fixed with updateObject().

- The GenomeDescription class was also moved from the GenomeInfoDb
  package to the new Seqinfo package.

- The following functions have moved to the new Seqinfo package:
  - constructor functions Seqinfo() and GenomeDescription()
  - seqinfo(), seqnames(), seqlevels(), seqlevels0()
  - seqlengths()
  - isCircular()
  - genome()
  - orderSeqlevels(), rankSeqlevels()
  - restoreSeqlevels()
  - sortSeqlevels(), seqlevelsInUse()
  - commonName(), provider(), providerVersion()
  - releaseDate(), bsgenomeName()

- The GenomeInfoDb package now depends on new Seqinfo package.

- GenomeInfoDbData was moved from Imports to Suggests (the package is
  only used by the loadTaxonomyDb() function).

BUG FIXES

- Fix mapping between UCSC genome dm3 and NCBI genome assembly for
  Drosophila melanogaster Release 5. See commit 3775ca4.

- Fix getChromInfoFromNCBI("Pan_troglodytes-2.1.4"). See commit
  4f3c737.

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.46.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The low-level CIGAR utilities and sequenceLayer() function are now in
  their own package, the cigarillo package, but with different names.
  Calling a low-level CIGAR utility function (or sequenceLayer())
  defined
  in GenomicAlignments still works but issues a deprecation warning
  that
  indicates the name of the new function to use from cigarillo.
  For example:
  > cigarOpTable(c("18M4S", "8M2X5M2I10M"))
  M I D N S H P = X
  [1,] 1 0 0 0 1 0 0 0 0
  [2,] 3 1 0 0 0 0 0 0 1
  Warning message:
  In call_new_fun_in_cigarillo("cigarOpTable", "tabulate_cigar_ops", :
  cigarOpTable() is formally deprecated in GenomicAlignments >= 1.45.5
  and replaced with the tabulate_cigar_ops() function from the new
  cigarillo package

- The package now depends on Seqinfo instead of GenomeInfoDb.

DEPRECATED AND DEFUNCT

- All the low-level CIGAR utilities and the sequenceLayer() function
  have moved to the new cigarillo package and are now deprecated in
  GenomicAlignments.

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.62.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Drop mirbase.db from Suggests (package got removed from BioC 3.22
  after being deprecated in BioC 3.21). As a consequence, function
  microRNAs() is now defunct.

- The package now depends on the new Seqinfo package instead of
  GenomeInfoDb.

DEPRECATED AND DEFUNCT

- microRNAs() is defunct.

- All the functions that went to txdbmaker in BioC 3.19 and got
  formally
  deprecated in GenomicFeatures in BioC 3.21 are now defunct.

[GenomicPlot](/packages/GenomicPlot)
-----------

                        Changes in version 1.7.2                        

- Add Seqinfo to Imports
- Move GenomeInfoDb from Imports to Suggests

                        Changes in version 1.7.1                        

- Updated test function
- Bump up version number to trigger Bioconductor build

                        Changes in version 1.7.0                        

- Import seqlevelsStyle from GenomeInfoDb
- Bump up version number

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.62.0                        

NEW FEATURES

- Add makeGPosFromDataFrame().

SIGNIFICANT USER-VISIBLE CHANGES

- The package depends on Seqinfo instead of GenomeInfoDb.

- The package no longer depends on XVector.

[geomeTriD](/packages/geomeTriD)
---------

                       Changes in version 1.13.21                       

IMPROVEMENT

- improve the reset button for tjviewer

                       Changes in version 1.13.20                       

IMPROVEMENT

- update the documentations in DESCRIPTION.

                       Changes in version 1.13.19                       

IMPROVEMENT

- update the documentations.

- add a reset button for tjviewer

                       Changes in version 1.13.17                       

BUG FIXES

- update the author list.

                       Changes in version 1.13.15                       

NEW FEATURES

- view3dStructure for multiple sequences.

                       Changes in version 1.13.14                       

NEW FEATURES

- Add scale bar in threejsViewer.

                       Changes in version 1.13.13                       

BUG FIXES

- Fix the issue for fill_NA if there are multiple NA values in dist
  matrix.

                       Changes in version 1.13.12                       

NEW FEATURES

- Add function to call TADs and compartments for spatial distance.

                       Changes in version 1.13.11                       

NEW FEATURES

- Add DSDC distance method for cellDistance.

                       Changes in version 1.13.10                       

IMPROVEMENT

- Fix multiple issues for PDF exporter.

                       Changes in version 1.13.9                        

NEW FEATURES

- Add NMI and ARI distance method for cellDistance.

                       Changes in version 1.13.8                        

IMPROVEMENT

- Add lwd controlor for GenomicInteraction links.

                       Changes in version 1.13.7                        

IMPROVEMENT

- Change the default value of normalized for TSS measurement to false.

                       Changes in version 1.13.6                        

NEW FEATURES

- Export the function cellDistance.

                       Changes in version 1.13.5                        

IMPROVEMENT

- Multiprocess for cellClusters.

                       Changes in version 1.13.4                        

NEW FEATURES

- Add function cellClusters.

                       Changes in version 1.13.2                        

BUG FIXES

- Move all geos into a sub properties in tjviewer.

                       Changes in version 1.13.1                        

NEW FEATURES

- Add hide/show button for each element in tjviewer.

[GeomxTools](/packages/GeomxTools)
----------

                       Changes in version 3.12.1                        

- Added pathway information to feature metadata
- Allowed more input types for pheno data (data.frames, lab
worksheet)
during readNanoStringGeoMxSet
- Bug fix to qcProteinSignal figure axis limit when one or more
protein have 0 expression
- Updated tests for Mac and ggplot2 update

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

[gg4way](/packages/gg4way)
------

                        Changes in version 1.6.1                        

- Updates for new major version of ggplot2 and other minor
improvements

[ggmsa](/packages/ggmsa)
-----

                       Changes in version 1.15.1                        

- replace ggalt::geom_xspline() with ggfun::geom_xspline()
(2017-07-12, Sat)

[ggsc](/packages/ggsc)
----

                        Changes in version 1.7.1                        

- compatible with new ggplot2 > 3.5.2 (2025-07-02, Wed)

[ggspavis](/packages/ggspavis)
--------

                 Changes in version 1.15.8 (2025-10-05)                 

- add Visium HD and imaging-based visualization to vignette; add
  sticker

                 Changes in version 1.15.6 (2025-05-11)                 

- renamed functions plotSpatial() and plotSpatialQC() to plotCoords()
  and
  plotObsQC() to avoid function name clash with imcRtools package

                 Changes in version 1.15.1 (2025-05-07)                 

- support for Visium HD data; renamed functions plotSpots() /
  plotSpotQC() to
  plotSpatial() / plotSpatialQC() (contributions by Estella Dong)

[ggtree](/packages/ggtree)
------

                       Changes in version 3.99.2                        

- update gheatmap, geom_range and msaplot to compatible with ggplot2
v=4.0.0 (2025-10-16, Thu, #668, #672, #674)
- interactive ggtree (2025-09-16, #662)
- incorporated iggtree, https://github.com/YuLab-SMU/iggtree
- Added support for XStringSet in msaplot (2027-07-13, Sun, #631)

                       Changes in version 3.17.1                        

- consistent with ggplot2 v=4.0.0 (2025-07-12, Sat, #657)

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                       Changes in version 1.19.1                        

- supporting the ggplot2_v4. (2025-08-08, Fri)

[GloScope](/packages/GloScope)
--------

                 Changes in version 1.99.1 (2025-09-27)                 

- Changed default number of GMM components to 5,10,15,20

                 Changes in version 1.99.0 (2025-09-20)                 

- Changed name of argument to `gloscope` and `gloscope_proportion` from
  `dist_mat` to `dist_metric` to be more consistent in argument names

- Changed name of function `gloscope_proportion` to `gloscopeProp` for
  consistent naming conventions

- `gloscopeProp` now will calculated the total variation distance
  between cell-type proportions (`dist_metric="TV"`)

- `gloscopeProp` now will add the value `ep` to all counts, not just
  those that are zero.

- Changed `plotMDS` to be more versatile. ARGUMENTS HAVE CHANGED.
  Instead of `group_id`, user can now give any combination of
  `color_by` and `shape_by` to indicate variables to use for the
  ggplot. The function has also removed unnecessary ggplot commands,
  leaving those to the user to chose.

- Added function `plotHeatmap` as a wrapper to `pheatmap` function for
  plotting a heatmap of the GloScope divergence matrix

- Added functions `getMetrics`, `bootstrap_gloscope`, `bootCI` to
  perform hypothesis testing and confidence intervals for testing
  differences between groups.

- Added function `plotCI` to plot output of `bootCI`

                 Changes in version 1.7.2 (2025-09-17)                  

- Added functionality for custom GMM and kNN parameters in fitting

- Added aditional test cases and update documentation

                 Changes in version 1.7.1 (2025-07-23)                  

- Changed default number of GMM components to 6,8,...,24

[GNOSIS](/packages/GNOSIS)
------

                 Changes in version 1.99.0 (2023-09-04)                 

- Made the following significant changes
  o added functionality to select and upload cBioPortal study
  o deprecated ability to save R script with executed code

- Submitted to Bioconductor

[GOexpress](/packages/GOexpress)
---------

                       Changes in version 1.43.1                        

BUG FIXES

- Remove now unused argument 'curl' from all calls to biomaRt::getBM().

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.35.2                        

- bug fixed that introduced in previous commit (2025-10-03, Fri)
- https://github.com/YuLab-SMU/DOSE/issues/86#issuecomment-3361787836

                       Changes in version 2.35.1                        

- alternative url for load_onto() (2025-08-12, Tue)

[graphite](/packages/graphite)
--------

                 Changes in version 1.55.2 (2025-10-22)                 

- Reinstated clipper support.

                 Changes in version 1.55.1 (2025-10-22)                 

- Automatic cleanup of stale pathway data.

- Updated all pathway data.

[GSVA](/packages/GSVA)
----

                        Changes in version 2.40                         

USER VISIBLE CHANGES

- Input needs not row/feature/gene names in the expression data
  anymore, and gene sets can be directly integer values indexing rows
  in the expression data. This is mostly useful for quick testing
  purposes.

- gsvaEnrichment(plot="ggplot") outputs now ggplot2 4.0 compliant
  plots.

- Input (sparse) expression data for the GSVA method can be now also
  provided as a SVT_SparseMatrix object, defined in the SparseArray
  package, either directly or within a SingleCellExperiment object.
  This allows for analysing single-cell data sets with more than 2^31
  nonzero values.

- Input expression data for the GSVA method can be now also provided as
  a DelayedArray matrix storing the values using an HDF5 on-disk
  backend, either directly or within a SingleCellExperiment object.

- Added a new vignette illustrating the use of the GSVA method with
  single-cell RNA-seq data.

[gVenn](/packages/gVenn)
-----

- Initial Submission to Bioconductor

[h5mread](/packages/h5mread)
-------

                        Changes in version 1.2.0                        

- No changes in this version.

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.38.0                        

- No changes in this version.

[hermes](/packages/hermes)
------

                       Changes in version 1.13.1                        

Bug Fixes

- Fixed issue with draw_boxplot for new ggplot2 version.

[HiCPotts](/packages/HiCPotts)
--------

- Initial Submission to Bioconductor

[HoloFoodR](/packages/HoloFoodR)
---------

                        Changes in version 1.3.1                        

- getMetaboLights: Improve data import.

[HPAanalyze](/packages/HPAanalyze)
----------

                        Changes in version 1.27                         

- Changes in version 1.27.1
  + Rework hpaDownload for more resilience against changing APIs.
  + Update function documentations
  + Update built-in data and internal gene name lookup table to v24.0
  + Small bug fixes.

- Changes in version 1.27.2
  + Update vignettes
  + Optimize package size
  + Better code comments on download and vis functions

[HuBMAPR](/packages/HuBMAPR)
-------

                        Changes in version 1.3.2                        

- Error corrections

                        Changes in version 1.3.1                        

- Update the functions to display the information properly

[Ibex](/packages/Ibex)
----

- Initial Submission to Bioconductor

[iCOBRA](/packages/iCOBRA)
------

                       Changes in version 1.37.1                        

- Minor fixes to tests to be more robust to possible future changes in
  ggplot2

- Minor rewriting to fix some check notes

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.15.5 (2025-10-20)                 

- fixed machine imprecision in tests

- fixed test_plotSpatialContext

                 Changes in version 1.15.4 (2025-10-08)                 

- added readImagefromTIFF function and tests

                 Changes in version 1.15.3 (2025-09-26)                 

- Added general zscore output to 'testInteractions()'

- Changed name "histocat" method to "conditional" method in
  countInteractions and testInteractions

                 Changes in version 1.15.2 (2025-06-18)                 

- Added "interaction" method to countInteractions and testInteractions

- Added plotInteractions function

                       Changes in version 1.15.1                        

- Replaced minDistToCells function with distToCells function

[immApex](/packages/immApex)
-------

                        Changes in version 1.3.7                        

UNDERLYING CHANGES

- Ensure buildNetwork() returns a symmetric distance matrix and drops
nonstructural zeros

                        Changes in version 1.3.6                        

UNDERLYING CHANGES

- Rebuilt propertyEncoder() and onehotEncoder() to use C++ backend
with encodeSequences.cpp
- Rebuilt buildNetwork() via C++ integration (fastEditEdges.cpp)
- Motif quantification in calculateMotif() use C++
(calculateMotif.cpp)
- Reduce overall dependencies on external packages
- Improved speed of generateSequences(), mutateSequences(),
adjacencyMatrix(), tokenizeSequences(), geometricEncoder
- Added group.by argument to getIR()
- Defunct variationalSequences()
- Removed keras3/tensorflow suggests/dependencies
- Removed baslisk environment, python no longer required
- Updated C++ backend to C++17
- Robust checks before running getIMGT() for website functioning
- Removing special unicode symbols

NEW FEATURES

- calculateEntropy() added to calculate the positional entropy along
a
biological sequence
- Added basic diversity metrics (exported) to support
calculateEntropy()
- calculateFrequency() added to calculate the positional frequency
along a biological sequence
- calculateMotif() added to get motif quantification of sequences
- calculateGeneUsage() added for single/paired gene enumeration
- Added scaleMatrix() for comprehensive scale/transformation
functions
- Added summaryMatrix() for fast summarization of matrix values

BUGS

- Fixed NT position issue with inferCDR()
- Fixed generateSequence() range issue for single-length sequences
- Fixed binary frequency issue in positionalEncoder()
- Fixed buildNetwork() for relative thresholding returns relative
value

                        Changes in version 1.2.2                        

- Updated unit tests and vignette check
- Converted package to basilisk from reticulate

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.41.1                        

- -CHANGES
  * experimental trail for sub-clustering using kmeans

[iModMix](/packages/iModMix)
-------

                        Changes in version 1.0.0                        

- **Initial Release**
  - Introduction of `iModMix`, a novel network-based approach that
  integrates multi-omics data, including metabolomics, proteomics, and
  transcriptomics data.
  - Ability to incorporate both identified and unidentified
  metabolites.
  - Functions for loading and preprocessing data, performing partial
  correlation analysis, hierarchical clustering, and module eigengene
  calculation.
  - Integration of multiple omics datasets to reveal significant
  associations with phenotypes of interest.
  - Visualization tools for dendrograms, module assignments, and
  network plots.
  - Enrichment analysis using gene symbols with support for various
  databases.
  Key Features
  - `load_data()`: Function to load and preprocess expression data.
  - `partial_cors()`: Function to perform partial correlation analysis
  using Graphical Lasso.
  - `hierarchical_cluster()`: Function to perform hierarchical
  clustering using TOM or partial correlation matrix.
  - `cluster_assignments()`: Function to assign clusters to features
  based on hierarchical clustering.
  - `Eigengenes()`: Function to calculate module eigengenes.
  - `perform_classification()`: Function to relate modules to
  phenotypes using statistical analysis.
  - `Modules_correlation()`: Function to integrate multiple omics
  datasets and calculate correlations.
  - `Assigment_genes_enrichr()`: Function to perform enrichment
  analysis using gene symbols.
  - Created documentation and examples for the functions.
  - Added vignette to demonstrate the usage of the package.

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.44.0                        

NEW FEATURES

- Add makeIRangesFromDataFrame(). Also add:
  - coercion from data.frame or DataFrame to IRanges,
  - support for IRanges(<data.frame>) and IRanges(<DataFrame>).

- Define Views() methods for integer and numeric vectors (they require
  the XVector package).

- Start experimenting with multithreading in findOverlaps().

SIGNIFICANT USER-VISIBLE CHANGES

- Significantly reduce findOverlaps() memory footprint. See commit
  247ed4e.

- Some clarification in nearest()'s documentation.

BUG FIXES

- Small tweak to viewApply() method for Views objects.

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.19.3 (2025-07-21)                 

UPDATE

- colour update for .alluvial_plot

                 Changes in version 1.19.2 (2025-07-17)                 

UPDATE

- vegan package switched from Suggests to Imports
- added refGene_hg38 for CIS_grubbs test

                 Changes in version 1.19.1 (2025-07-09)                 

UPDATE

- Changed again color scale for alluvial_plot

                 Changes in version 1.19.0 (2025-07-09)                 

GENERAL UPDATE

- Changed color scale for alluvial_plot and made some general updates

[iscream](/packages/iscream)
-------

- Initial Bioconductor submission

[iSEE](/packages/iSEE)
----

                       Changes in version 2.21.1                        

- Fixed dimensions of exported figures.

[islify](/packages/islify)
------

                        Changes in version 1.1.1                        

- Introduction of the problematicUnevenDistribution flag in the core
  islify function.

[IsoBayes](/packages/IsoBayes)
--------

                        Changes in version 1.6.1                        

- we added a normalization of intensity values (they add to 10,000)

[isomiRs](/packages/isomiRs)
-------

                       Changes in version 1.37.1                        

- Remove functions that were using targetScan since it is not available
  anymore.

[jazzPanda](/packages/jazzPanda)
---------

                        Changes in version 1.1.2                        

Updated on Oct 14, 2025

- Minor bug fix for Bioconductor devel.

                        Changes in version 1.1.1                        

Updated on June 24, 2025

Changes to the package include:

- function get_vectors():

- Parameters w_x and w_y have been removed. Spatial windows are
now automatically determined per sample using a 5% buffer around
the coordinate range.
- Improved handling of genes that are only detected in a subset of
samples. If a gene is not present in a sample, a zero-filled
vector is used for that sample to maintain consistent vector
lengths across samples.

- function create_genesets():

- Parameters w_x and w_y have been removed. Spatial windows are
now computed per sample using a 5% buffer based on input
coordinates.
- Improved handling of genes that are only detected in a subset of
samples. If a gene is not present in a sample, a zero-filled
vector is used for that sample to maintain consistent vector
lengths across samples.

- function compute_permp():

- Parameters w_x and w_y are no longer required. Spatial
boundaries are now automatically derived per sample using a 5%
buffer.

                        Changes in version 1.0.1                        

Updated on Sep 1, 2025

Changes to the package include:

- function get_vectors():

- Buffer used for auto-detecting spatial windows has been changed
from 5% to a fixed value of 1e-6.
- Added a new parameter to allow exporting the calculated bounding
box for each sample.

- function compute_permp():

- Buffer for automatic spatial boundaries is now set to 1e-6.
- Fixed spatial bin indexing issues to ensure correct alignment
when calcualting spatial vectors with count matrix.

- Backported updates from version 1.1.1 on Bioconductor devel

[kebabs](/packages/kebabs)
------

                       Changes in version 1.43.1                        

- updated DESCRIPTION file to new format (Authors@R)

                       Changes in version 1.43.0                        

- new branch for Bioconductor 3.22 devel

[KEGGREST](/packages/KEGGREST)
--------

                       Changes in version 1.48.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Defunct `mark.pathway.by.objects` function (no longer supported).

[lefser](/packages/lefser)
------

                       Changes in version 1.20.0                        

Significant user-visible changes

- Removed blockCol and groupCol arguments in lefser function.

[limma](/packages/limma)
-----

                 Changes in version 3.66.0 (2025-10-28)                 

- 
  New argument `upshot` for treat().

- 
  New argument `span` for fitFDistUnequalDF1(), squeezeVar() and
  eBayes().

- 
  The default for `adaptive.span` in voom() and
  voomWithQualityWeights() is now `TRUE`. All the voom functions
  now use chooseLowessSpan() by default to choose the `span`.

- 
  New argument `xlab` for plotSplice() to allow labeling in terms
  of transcripts or other splicing events instead of exons.

- 
  The help pages for diffSplice(), topSplice() and plotSplice()
  have been updated to emphasize differential transcript usage
  instead of differential exon usage or differential splicing.
  "differential usage" concept tags have been added to the
  appropriate help pages, creating an entry in the reference
  manual index.

- 
  makeContrasts() now uses the contrast names if `contrasts` is a
  named vector.

- 
  Fix typo in error message from makeContrasts() when the
  variable names are not syntactically valid.

- 
  Bug fix to contrastAsCoef() when the `contrast` matrix is not
  of full column rank.

- 
  The columns of the `covariates` matrix input to
  removeBatchEffect() are now mean-corrected.

- 
  removeBatchEffect() now gives a message if `design` is not
  specified.

- 
  Bug fix for topTable() when `confint=TRUE` and either `p.value`
  is less than 1 or `lfc` is greater than zero.

- 
  Revise help page for lmFit(), removing out-of-date mention of
  default correlation.

- 
  Minor edits to poolVar() and update of help page.

- 
  Complete remaining dates for release versions in NEWS.Rd.

[limpa](/packages/limpa)
-----

                 Changes in version 1.2.0 (2025-10-28)                  

- 
  New function readFragPipe() to read FragPipe output.

- 
  New function readMaxQuant() to read MaxQuant output.

- 
  readDIANN() now allows user-specified input column names and
  "parquet" is now the default format.

- 
  Bug fixes to readSpectronaut() when `sep` is other than a tab
  or when there is only one peptide annotation column.

- 
  Revisions to all the read functions. readDIANN(),
  readFragPipe(), readMaxQuant(), and readSpectronaut() all use
  `fread()` with `data.table=FALSE` now to read ascii delimited
  files.

- 
  New argument `logit.p` for dztbinom(), which allows the
  function to avoid subtractive cancellation for success
  probabilities close to 1. dpc() now uses the new argument when
  calling dztbinom().

- 
  New function dpcON() to estimate the DPC assuming the
  oberved-normal model. Should give same output as dpc() but is
  somewhat faster. dpcON() calls a new function,
  pztbinomSameSizeLogitPBothTails(), which computes left and
  right tail p-values for zero-truncated binomial quantiles.

- 
  Revisions to dpcCN() to improve convergence. To improve speed,
  the function now randomly selects a subset of rows of data to
  process, with the number of rows given by the new argument
  `subset`. The default number of outer iterations is decreased
  to 2.

- 
  dpcImpute() now sets row.names. Sections on dpcImpute() and
  imputeByExpTilt() have been added to the package vignette.

- 
  dpcQuant() now calls dpcImpute() automatically when all
  proteins have just one peptide.

- 
  New function plotPeptides() to plot the peptide-level
  log-intensities for one protein.

- 
  Fix simProteinDataSet() processing of `prop.missing` argument.

- 
  Cite the limpa bioRxiv preprint (Li et al, bioRxiv 2025) in the
  vignette and help pages.

- 
  Expand Quantification section in the vignette to give more math
  details.

[LipidTrend](/packages/LipidTrend)
----------

                 Changes in version 0.99.0 (2025-02-25)                 

- Submitted to Bioconductor

[LRBaseDbi](/packages/LRBaseDbi)
---------

                       Changes in version 2.18.1                        

- A vignette modified.

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


[Macarron](/packages/Macarron)
--------

                       Changes in version 1.13.2                        

- Changed NAs to zeros in calES.R and calQval.R

- Features common across AVA, Qval, and ES are prioritized

                       Changes in version 1.13.1                        

- Kept NAs as NAs in calES.R and calQval.R

[markeR](/packages/markeR)
------

- Initial submission to Bioconductor 

[MassSpecWavelet](/packages/MassSpecWavelet)
---------------

                 Changes in version 1.75.1 (2025-05-03)                 

- identifyMajorPeaks():

- Properly handle SNR.Th=0 and excludeBoundariesSize=0 cases to
bypass those filtering criteria.
- Refactor and include additional peak details such as
peakRidgeLengthScale, that can be used to tune ridgeLength;
peakNoise, useful to understand the role of w.g. the
winSize.noise argument; and selInd, useful to understand which
filtering rule is impacting a potential peak from not being
identified.

- Fix NOTE on escaping special LaTeX characters

[mastR](/packages/mastR)
-----

                        Changes in version 1.9.1                        

- Add Bioinformatics citation for mastR package, add sticker and
remove BisqueRNA dependency.

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.17.1 (2025-06-26)                 

- fix error with respect to setExperimentHubOption in vignette

[matter](/packages/matter)
------

                       Changes in version 2.11.2                        

NEW FEATURES

- Add 'as_plotly()' for coercing vizi plots

                       Changes in version 2.11.1                        

BUG FIXES

- Fix S3/S4 method issues with 'bplapply()'

- Improve to 'sgmix()' numeric stability

[memes](/packages/memes)
-----

                       Changes in version 1.17.1                        

- Fixes an error when reading XML files with jagged evalue metadata
entries

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                        Changes in version 1.13                         

Changes in 1.13.3

- Expand documentation for CompareSpectraParam to also mention use of
msentropy or GNPS-like similarity calculation.

Changes in 1.13.2

- For plotSpectraMirror,MatchedSpectra method, we now pass the
available ppm and tolerance value that were used in the creation of
the object.

Changes in 1.13.1

- Add parameter matchedPeaksCount to CompareSpectraParam that enables
reporting of the number of matched peaks for matchSpectra() with
CompareSpectraParam. The result are reported in a column
"matched_peaks_count".

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                        Changes in version 1.17                         

Changes in 1.17.1

- Fix adductFormula() for negative adducts. and added respective unit
test.

[MetaboDynamics](/packages/MetaboDynamics)
--------------

                       Changes in version 1.99.12                       

- contains a vignette describing the package workflow with a data
frame input
- altered plot_cluster and plot_ORA visualization for patchwork
option
enabling easier visualization
- all functions require columns names "metabolite","condition","time"
for consistency throughout the package
- get_ORA_annotations() function (retrieving KEGG IDs hierarchies)
added back to package
- no more samples in estimates dynamics as probability of clustering
solution is implemented now as bootstrapping from model posterior in
function cluster_dynamics
- estimates_dynamics() now returns a list with: estimated metabolite
abundance (mu), estimated standard deviation of metabolite abundance
(sigma), estimated pooled standard deviation per metabolite and dose
(lambda), differences in metabolite abundances between time points,
euclidean distance between metabolite dynamics of different
conditions
- plot_estimates() additionally visualizes the euclidean distance
between conditions of metabolite specific dynamics
- cluster_dynamics() provides clustering solution of mean estimates
of
mu as well as bootstraps clustering solutions
- plot_cluster provides bubbletree, cluster identity, dynamics plots
as well as patchwork plot (combining bubbletree, cluster identity
and dynamics plots)
- plot_ORA() has now option to be added to bubbletree obtained by
plot_cluster
- ORA_hypergeometric() requires now a data frame annotationg
metabolites to KEGG IDs (column names "metabolite" and "KEGG")
- internal adding of ANOVA model with euclidean distance estimation
between doses and ANOVA model that integrates cell count uncertainty
- internal adding of simulated cell counts to
data("longitudinalMetabolomics")
- differences between time points are now ordered and return is a
list
of plots

                        Changes in version 1.0.2                        

Minor bug fix in vignette that caused errors in package checks

                        Changes in version 1.0.1                        

bug fix: fit_dynamics_model() can now correctly handle provided
variable
names

[methylInheritance](/packages/methylInheritance)
-----------------

                       Changes in version 1.33.1                        

BUG FIXES AND IMPROVEMENTS

- Ensure functions using ggplot2 are using update functions

[MetNet](/packages/MetNet)
------

                 Changes in version 1.27.4 (2025-10-24)                 

- adjust unit tests for aracne

                 Changes in version 1.27.3 (2025-10-23)                 

- add 1 as diagonal values in adj in vignette

                 Changes in version 1.27.2 (2025-10-10)                 

- extend unit test to larger data table in test_methods.R

                 Changes in version 1.27.1 (2025-10-02)                 

- adjust unit tests after ggplot2 changes

[mia](/packages/mia)
---

                        Changes in version 1.17                         

- Added add.dimred parameter to meltSE() (2025-04-29)

- Added transformation method that replaces values <= threshold with NA
  (2025-05-09)

- importMetaPhlAn: Preserve features without taxonomy

- Added agglomerateByModule()

- Added inverse rank normalization (2025-09-10)

[miaDash](/packages/miaDash)
-------

                        Changes in version 1.1.5                        

- Added Cluster tab

                        Changes in version 1.1.4                        

- Switch from csv to tsv support

- Improve interoperability with MGnify datasets

                        Changes in version 1.1.3                        

- Introduced tab for Quality Control

                        Changes in version 1.1.2                        

- Added importers for HUMAnN, QIIME2 and Mothur

- Improved biom importer

[miaTime](/packages/miaTime)
-------
- Initial submission to Bioconductor

[miaViz](/packages/miaViz)
------

                        Changes in version 1.17                         

- Improve *Histogram and *Barplot by adding options for filling and
  facetting (2025-04-24)

- Added plotBoxplot function (2025-04-27)

- plotSeries: Improve by adding functionality for sample groups
  (2025-04-28)

[microbiome](/packages/microbiome)
----------

                 Changes in version 1.32.2 (2025-07-17)                 

- Fixed NEWS file

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.20.1                        

- rm gghalves since it is not compatible with ggplot2 4.0.0.
(2025-09-17, Wed)

[MIRit](/packages/MIRit)
-----

                        Changes in version 1.5.1                        

The link for downloading miRTarBase v10 has been updated.

[miRSM](/packages/miRSM)
-----

                        Changes in version 2.5.4                        

- Update s4vd and subspace <2025-06-22, Sun>

                     Changes in version 2.5.1-2.5.3                     

- Update linkcomm <2025-06-15, Sun>

[miRspongeR](/packages/miRspongeR)
----------

                       Changes in version 2.13.5                        

- Update Groundtruth <2025-10-16, Thur>.

                    Changes in version 2.13.1-2.13.4                    

- Update linkcomm method <2025-06-15, Sun>.

[mitology](/packages/mitology)
--------

                        Changes in version 1.2.0                        

- New enrichMito function to perform a mitochondrial enrichment
  analysis of a gene list.

- New gseaMito function to perform a mitochondrial GSEA of a gene
  list.

- New mitoTreePoint function draw a circular dotplot on the
  mitochondrial gene set tree.

[mobileRNA](/packages/mobileRNA)
---------

                        Changes in version 1.5.2                        

- Updated several functions including RNAfeatures.

[monaLisa](/packages/monaLisa)
--------

                       Changes in version 1.15.3                        

- minor updates to adjust to breaking changes in ggplot2 version 4

                       Changes in version 1.15.2                        

- minor updates to adjust to functionality migrated from GenomeInfoDb
to new Seqinfo package

[Moonlight2R](/packages/Moonlight2R)
-----------

                        Changes in version 1.7.4                        

Summary

- Updated GLS() and GMA() examples to reduce the check timings

                        Changes in version 1.7.3                        

Summary

- fixed error in GSEA() module due to issues in DOSE package devel
version

                        Changes in version 1.7.2                        

Summary

- fixed error in plotDMA and other plotting functions due to updates
in tidyHeatmap package

                        Changes in version 1.7.1                        

Summary

- fixed error in GLS function due to updates in easyPubMed package

                        Changes in version 1.7.0                        

- version bump due to release of Bioconductor 3.21

[mosdef](/packages/mosdef)
------

                        Changes in version 1.6.0                        

- Slight edits to have better behaved de_volcano() and go_volcano()
functions, with filters/looks depending on the adjusted p-values,
and a better management of the color values to be used when plotting
the different subsets

[MotifPeeker](/packages/MotifPeeker)
-----------

                        Changes in version 1.1.2                        

Miscellaneous

- Only consider unique peaks in peak count.
- Only return unique peaks in segregate_seqs().

Bug Fixes

- Sanitise peak input names before running FIMO.

                        Changes in version 1.1.1                        

Bug Fixes

- Fix importing MACS peak files with no peak names.
- Update tests for messager().

Miscellaneous

- Add citation to pre-print.
- Switch to the official repo for rworkflows.

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.53.2                        

- Speedup motifHclust

                       Changes in version 1.53.1                        

- Support import Hocomoco format.

[motifTestR](/packages/motifTestR)
----------

                        Changes in version 1.5.6                        

- Added GLM-Poisson model for enrichment testing

                        Changes in version 1.5.1                        

- Added simMultiMotifs

[MPAC](/packages/MPAC)
----

                        Changes in version 1.3.1                        

- updated citation to point to the journal Bioinformatics

[msa](/packages/msa)
---

                       Changes in version 1.41.3                        

- fix in ClustalOmega to better determine whether OpenMP is available

- updated DESCRIPTION file to new format (Authors@R)

                       Changes in version 1.41.2                        

- fix of memory leaks in Muscle (thanks for GitHub user Rong-Zh)

                       Changes in version 1.41.1                        

- fix of compilation issues caused by ClustalOmega on arm64 systems

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                        Changes in version 1.17                         

Changes in 1.17.1

- Allow columns/fields "comment" and "data_processing_comment" to
have
multiple values.

[MsBackendMetaboLights](/packages/MsBackendMetaboLights)
---------------------

                         Changes in version 1.3                         

Changes in 1.3.2

- Support backendMerge() (issue #10).

Changes in 1.3.1

- Fix download/sync of files with BiocFileCache.

[MsBackendSql](/packages/MsBackendSql)
------------

                         Changes in version 1.9                         

Changes in 1.9.4

- Add additional unit tests to check results are similar to the
reference implementation.

Changes in 1.9.3

- Small performance improvement in longForm() implementation: avoid
using findMatches() to order the result if not needed (i.e., if
duplicated spectra are requested).

Changes in 1.9.2

- Add a longForm() method for MsBackendSql to use an SQL query to
extract the data in long from the database.

Changes in 1.9.1

- Add support for blob2 storage mode for duckdb databases.

[MsExperiment](/packages/MsExperiment)
------------

                        Changes in version 1.11                         

Changed in 1.11.1

- Fix readMsExperiment(): if sampleData has row names, they have to
match the names of the data files.

[MsFeatures](/packages/MsFeatures)
----------

                        Changes in version 1.17                         

MsFeatures 1.17.1

- Add new groupSimilarityMatrixTree() function.

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.35                         

2.35.2

- Second fix error from adding Selenocysteine and Pyrrolysine to
PSMatch AA table.
- Import PSmatch >= 1.13.3.

2.35.1

- Fix error from adding Selenocysteine and Pyrrolysine to PSMatch AA
table.

2.35.0

- New devel version

[msqrob2](/packages/msqrob2)
-------

                        Changes in version 1.17                         

msqrob 1.17.1

- Fixed NEWS file
- Fixed vignette upon renaming of longFormat to longForm

msqrob 1.17.0

- New Bioconductor devel release 3.22

[MsQuality](/packages/MsQuality)
---------

                         Changes in version 1.9                         

Changes in version 1.9.3

- change unit tests for filterEmptySpectra after update

Changes in version 1.9.2

- resolve NAMESPACE issue of unexported objects

Changes in version 1.9.1

- fix bug due to update in rmzqc package, use $new for constructing
R6
objects

[MSstatsResponse](/packages/MSstatsResponse)
---------------

                 Changes in version 1.0.0 (2025-09-22)                  

- Submitted to Bioconductor

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.36.0                        

New features

- Defunct longFormat generic and methods replaced with longForm
- Added hasRowData method check for assays
- subsetByRowData and intersectByRowData allow filtering of rows
based
on values in assay rowData (#228)

Bug fixes and minor improvements

- Properly document i subsetByRow and [ bracket methods.
- Match subsetByRow method signatures to generic.

[MultiRNAflow](/packages/MultiRNAflow)
------------

                        Changes in version 1.7.2                        

- Vignette (MultiRNAflow_vignette-knitr.Rnw)
- We use now use.unsrturl=TRUE

                        Changes in version 1.7.1                        

- R functions of the packages
- myYlabelNorm() (included in the R function
DATAplotExpressionGenes()): there were an error when the RPKM
normalization were used.
- GSEAQuickAnalysis():
- replace size aesthetic by linewidth because : Using size
aesthetic for lines was deprecated in ggplot2 3.4.0.
- correction errors: print(gManhattan) to
print(gproList$gManhattan)
- MFUZZclustersNumber(): slight modification for correcting one
warning "Removed 1 row containing missing values or values
outside the scale range (geom_point())."
- DEplotBarplot(): in the section "Examples", we replace
c("Spe.Pos", "Spe.Neg", "Other") by c("UpRegulated",
"DownRegulated", "Other").
- DEplotBarplotFacetGrid(): in the section "Examples", we add
Melt.Dat.2 <- stats::aggregate(Nb.Spe.DE~., data=Melt.Dat.2,
sum) so Melt.Dat.2 is more representative of the data used by
our function. We also modify the function to prevent recent
warnings from ggplot2 (.data[["<col_names>"]]).
- Vignette (MultiRNAflow_vignette-knitr.Rnw)
- We add our publication from Bioinformatics
- we correct new errors produced by latex output
- problem with the input max.level of the function str()
- we correct the following new error 'Illegal parameter number
in definition of \NewValue'
- we use '\printbibliography' and '\addbibresource' in the
preamble for the 'cleaner' biblatex interface.

[multistateQTL](/packages/multistateQTL)
-------------

                 Changes in version 2.1.2 (2025-10-07)                  

- Updated callSignificance so that it works with one state.

                 Changes in version 2.1.1 (2025-07-14)                  

- Moved data.table to imports instead of depends.
- Updated testing to avoid error with ggplot2 version 4.0.0.

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.17.5                        

Bug fix

- Before inferring effect column ensure A1 and A2 are uppercase.

                       Changes in version 1.17.3                        

Bug fix

- Graciously handle unsinking messages when run fails.

                       Changes in version 1.17.2                        

Bug fix

- Bug in long tests utilising ieugwasr - swap to new token system

                       Changes in version 1.17.1                        

New features

- MungeSumstats can now handle local versions of dbSNP in tarball
format (using the dbSNP_tarball parameter in format_sumstats). This
was enabled to help with dbSNP versions >=156, after the decision to
no longer provide dbSNP releases as bioconductor packages. dbSNP 156
tarball is available here:http://149.165.171.124/SNPlocs/.

[muscat](/packages/muscat)
------

                       Changes in version 1.23.2                        

- fixed unit test with edge-case failure

- added missing '@importFrom grDevices hcl.colors'

- fixed 'BiocParallel::SerialParam progressbar' arg

                       Changes in version 1.23.1                        

- added bulk-based hypothesis weighing

                       Changes in version 1.23.0                        

- Bioconductor 3.21 release

[musicatk](/packages/musicatk)
--------

                 Changes in version 2.2.1 (2025-07-15)                  

- Removed conclust dependency
- Fixed plotting issues with indels
- Cosmetic changes to k value and plotting code

[mutscan](/packages/mutscan)
-------

- Initial Bioconductor submission

[mzR](/packages/mzR)
---

                       Changes in version 2.43.3                        

- Make electronBeamEnergy an optional header field.

                       Changes in version 2.43.2                        

- More fixes around NaN

                       Changes in version 2.43.1                        

- Update CV, adds multiple terms

[NanoStringNCTools](/packages/NanoStringNCTools)
-----------------

                 Changes in version 1.16.2 (2025-09-26)                 

- Bug fixes for dependency package updates

[notame](/packages/notame)
------

                         Changes in version 1.0                         

- split notameViz and notameStats from notame

- removed MetaboSet-support

- included SummarizedExperiment-support

- parallelization with BiocParallel

- included website with rworkflows (reference index and vignette also
  cover notameViz and notameStats)

- removed perform_pairwise_t_test() and perform_paired_t_test(),
  replaced by perform_t_test()

- removed perform_pairwise_non_parametric() and
  perform_wilcoxon_signed_rank(), replaced by perform_non_parametric()

- removed align_batches(), replaced by batchCorr::alignBatches()

- removed normalize_batches(), replaced by
  batchCorr::normalizeBatches()

- other removed functions: install_dependencies(), merge_exprs(),
  plot_sample_suprahex(), prop_found(), prop_na(), finite helpers such
  as finite_mean() and simple transformations such as exponential()

- removed plotting in correct_drift(), replaced by save_dc_plots()

- removed plotting in cluster_features(), replaced by
  visualize_clusters()

- de-exported low-level functionality for drift correction and feature
  clustering

- renamed read_from_excel() to import_from_excel()

- renamed clean_stats_results() to summarize_results()

- renamed visualizations() to save_QC_plots()

- renamed merge_objects() to merge_notame_sets()

- renamed example_set to toy_notame_set and removed other example data

- reworked vignettes

- reworked plotting functionality in mixOmics_pls() and other mixOmics
  wrappers

- updated flag_contaminants to use whole number instead of ratio in
  threshold and include more metrics to flag by

- updated muvr_analysis() to use MUVR2 instead of MUVR package

                       Changes in version 0.99.0                        

- submitted to Bioconductor

[notameStats](/packages/notameStats)
-----------

- Initial Bioconductor submission

[notameViz](/packages/notameViz)
---------

- Initial Bioconductor submission

[omicsGMF](/packages/omicsGMF)
--------

- Initial Bioconductor submission

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

                 Changes in version 4.11.1 (2025-06-09)                 

- Example using WT in interventions.

[orthogene](/packages/orthogene)
---------

                       Changes in version 1.15.02                       

Bug fixes

- Fix #46

New features

- Disable blank Issues
- Update inst/orthogene manual

                       Changes in version 1.15.01                       

Bug fixes

- Update ggplot2 facet calls to avoid deprecation errors.
- Fix all roxygen note warnings.

[orthos](/packages/orthos)
------

                        Changes in version 1.7.3                        

- Update environment definitions

[OSTA.data](/packages/OSTA.data)
---------

                 Changes in version 1.1.2 (2025-08-28)                  

- add sticker

                 Changes in version 1.1.1 (2025-07-21)                  

- change Oliveira from biorxiv to Nature Genetics

[OUTRIDER](/packages/OUTRIDER)
--------

                       Changes in version 1.26.2                        

- Add Optimal Hard Thresholding to determine optimal encoding dimension
  #69

- Add singular value plot for optimal encoding dimension

- Add highlight.label.size to method-plot.R

[peakCombiner](/packages/peakCombiner)
------------

- Initial Bioconductor submission

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.23.3 (2025-09-27)                 

- update unittest for `ggplot2` 4.0.0

- update vignette figures location

[Pedixplorer](/packages/Pedixplorer)
-----------

                        Changes in version 1.5.6                        

- Fix legend plotting when ggplot_gen is TRUE
- Reduce kindepth() warning and error messages to show necessary
information.
- Fix lwd usage instead of cex in ped_to_legdf()

                        Changes in version 1.5.5                        

- Add label_cex and label_size to control id, date and label text
position and size
- Fix csv reading with incoherent number of columns in rows
- Fix check box align_parents in ped_shiny()
- Fix auto_hint() error when couple present twice
- Fix stress computation when individual present more than twice
- Add ind_max_warning and ind_max_error for a better control of
pedigree size to plot in shiny application

                        Changes in version 1.5.4                        

- Add best_hint computation to the shiny app and split it in
different
subfunction
- Remove little used library from Imports
- Change shinyWidgets::pickerInput to shiny::selectInput
- Remove 0 from default na_strings
- Update citation to use published article in Bioinformatics
- Add error for _ present in any id columns
- Fix shiny app notification title and slider intermediary value
selection
- Fix full scale data reaffection on same column
- Improve position and size of proband label and assymptomatic symbol
for logo

                        Changes in version 1.5.3                        

- Add short explanation of hints() usage in Pedigree alignment
details
vignette.

                        Changes in version 1.5.2                        

- Fix resizing plot in shiny module by adding dependency
- Update website and readme

                        Changes in version 1.5.1                        

- Separate plot management in a shiny module
- Add plot_resize shiny module
- Add help messages using shinyhelper
- Add plink_to_pedigree function to convert plink files to Pedigree
object
- Add 0 as one of the default missing identifier value
- Set ggplot generation as independent without plotting on current
device
- All draw_* functions return a layer instead of a ggplot object
- Add back message for individuals not plotted
- Add autocompletion of missing twins relationship with
complete_twins()

[pengls](/packages/pengls)
------

                       Changes in version 1.15.2                        

- Include citation to publication through CITATION file

                       Changes in version 1.15.1                        

- Return named coefficients

[pgxRpi](/packages/pgxRpi)
------

                 Changes in version 1.5.3 (2025-05-15)                  

- Added filter_pattern parameter to pgxLoader to enable keyword-based
search of available filter terms.
- Improved warning messages to include HTTP status codes for failed
requests, and clarified domain information in queries involving
multiple domains.
- Simplified domain field description by retaining only the domain
address.

[PharmacoGx](/packages/PharmacoGx)
----------

                       Changes in version 3.13.1                        

- Updated only the sessionInfo fields example data sets due to errors

[PhyloProfile](/packages/PhyloProfile)
------------

                        Changes in version 2.0.5                        

- auto adjust plot hight for taxonomic labels

- consider 90 percent of input taxa for core gene ident. by
  default

                        Changes in version 2.0.4                        

- show co-orthologs in overview plot

                        Changes in version 2.0.3                        

- option to download 3D DIM plot as HTML

                        Changes in version 2.0.2                        

- data selection in 3D DIM plot

[PMScanR](/packages/PMScanR)
-------

- Initial Bioconductor submission

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.17.1 (2025-07-08)                 

- Added alternative way to set threshold in instances when
  threshold would be returned as NA

[podkat](/packages/podkat)
------

                       Changes in version 1.41.1                        

- updated DESCRIPTION file to new format (Authors@R)

                       Changes in version 1.41.0                        

- new branch for Bioconductor 3.22 devel

[poem](/packages/poem)
----

                 Changes in version 1.1.3 (2025-10-21)                  

- Fixed bug of metrics argument not being passed with graph input

                        Changes in version 1.1.2                        

- lowMemory argument added to getFuzzyPartitionMetrics()

- returnElementPairAccuracy argument removed from
  getFuzzyPartitionMetrics(),
  getFuzzyPartitionGlobalMetrics() and getFuzzyPartitionClassMetrics()

- SpatialAccuracy metric renamed to nsAccuracy for consistency

- Bug fixed in getEmbeddingGlobalMetrics() and
  getFuzzyPartitionMetrics()

[procoil](/packages/procoil)
-------

                       Changes in version 2.37.1                        

- updated DESCRIPTION file to new format (Authors@R)

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.49                         

Changes in version 1.49.3

- Bump version for Bioc devel

Changes in version 1.49.2

- Change link to codeberg/lpsvm-tl-code
- An optimised version of tagm_mcmc was implemented. It should be
faster and more memory efficient

Changes in version 1.49.1

- Add citation for new f1000 subcellular workflow
- Make defunct setAnnotationParam, getAnnotationParams,
showGOEvidenceCodes, getGOEvidenceCodes, addGoAnnotations,
orderGoAnnotations and associated helper functions

Changes in version 1.49.0

- New Bioc release

[PRONE](/packages/PRONE)
-----

                        Changes in version 1.3.2                        

- Now major changes, only applied same changes for Bioconductor
release 3.22 version than for version 3.21 -> PRONE 1.2.2


[PSMatch](/packages/PSMatch)
-------

                        Changes in version 1.13                         

PSMatch 1.13.2

- Add CITATION to the pre-print
(https://doi.org/10.31219/osf.io/62v9p_v2)
- Add possibility to split protein groups in PSM data with
makeAdjacencyMatrix (see issue #30)
- Removed showDetails() from setMethod("show", "PSM") as discussed in
issue #30
- Add USI parameter to plotSpectraPTM.
- Correct filterPsmFdr output message
- Update Fragments vignette.
- Improve labelFragments()runtime (see issue #25).
- Add Selenocysteine and Pyrrolysine to getAminoAcids().

PSMatch 1.13.1

- Refactor and improve plotSpectraPTM().

PSMatch 1.13.0

- New Bioconductor devel.

[PureCN](/packages/PureCN)
------

                       Changes in version 2.16.0                        

BUGFIXES

- Fixed a crash when GATK log-ratios and a normal coverage file were
  provided (#386)

- Fixed a few code warnings (e.g. #387)

[pwalign](/packages/pwalign)
-------

                        Changes in version 1.6.0                        

- No changes in this version.

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.19                         

QFeatures 1.19.4

- readQFeatures() doesn't ignore colData$quantCols (PR #234)
- Create unique precursor ids and use them in joinAssays() (PR #223)
- joinAssays() allows joining based on fcol (PR #247)

QFeatures 1.19.3

- Fixed joinAssays with fcol so that assayLinks are correctly created
(PR 247).

QFeatures 1.19.2

- Aggregate features optimisation.

QFeatures 1.19.1

- Fix duplicated fnames bug (see issue #237).
- Use longForm() methods (replacing longFormat(), now defunct).
- Add longForm,SummarizedExperiment, supporting rowvars and colvars,
as longForm,QFeatures does. Without this method, uses longForm,ANY
that doesn't.
- Defined a QFeatures object type (bulk or single-cell) in the
instance's metadata slot (developer features).

QFeatures 1.19.0

- New devel version

[qpgraph](/packages/qpgraph)
-------

                        Changes in version 2.44                         

BUG FIXES

- Switch some imports from GenomeInfoDb to the new Seqinfo package.


[QTLExperiment](/packages/QTLExperiment)
-------------

                 Changes in version 2.1.3 (2025-10-27)                  

- Changed the data loading function sumstats2qtle(). o Added a flag
to
skip checking for multiple association tests for combinations of
state, feature and variant. o Replaced a left_join with a match
which should be faster.

                 Changes in version 2.1.2 (2025-09-30)                  

- Added support for QTLExperiment objects with only one state, i.e.
one column. o sumstats2qtle() can read in data when input has only
one path. o mockQTLE() allows nStates = 1.

                 Changes in version 2.1.1 (2025-07-14)                  

- subset function is now a generic.

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.34.0                        

Bug fixes and minor improvements

- Added vignette chunk labels to both vignettes for easier debugging.

[RAIDS](/packages/RAIDS)
-----

                        Changes in version 1.7.3                        

BUG FIXES

o Change default parameters for snpgdsIBDKING() function

[Rarr](/packages/Rarr)
----

                         Changes in version 1.9                         

New features

- New functions to work with Zarr attributes have been added:
- read_zarr_attributes() reads Zarr v2 and v3 attributes
- write_zarr_attributes() only supports writing Zarr v2 attributes
for now.
- This package now has a pkgdown website, available at
https://huber-group-embl.github.io/Rarr/.
- Zarr v3 arrays are now supported for reading metadata via
zarr_overview().

Breaking changes

- zarr_overview(as_data_frame = TRUE) now returns information more in
line with the output of zarr_overview(as_data_frame = FALSE). In
particular:
- a new endianness column has been added to indicate the byte
order of the array data.
- the nchunks column is now a list column specifying the number of
chunks in each dimension, rather than a single integer giving
the total number of chunks.

Minor improvements

- An explicit error message is now given when attempting to read a
Zarr array version 3. This version will be supported in a future
release of Rarr.

Bug fixes

- .url_parse_other() now accounts for port numbers in host name and
colons in S3 buckets.
- writeZarrArray() now allows writing character arrays, and no longer
errors complaining about null 'nchar' argument value. Default of
'nchar' is now NULL.
- writeZarrArray() no longer silently and incorrectly fills the last
rows/columns when dim is not divisible by chunk_dim.
- The object name is no longer repeated (e.g., name.zarrname.zarr)
when writing a Zarr array to a file in the current working
directory.
- Invalid URLs for examples with S3 storage in read_zarr_array() and
zarr_overview() have been updated.
- read_zarr_array() no longer errors on arrays with numeric values
other than float, int, uint and complex.
- zarr_overview() now returns an explicit error message when the
.zarray file is absent

Internal changes

- Coding style throughout the package has been harmonized using the
air tool. Contributors using RStudio, Positron or VS Code should
have their code styled automatically on save.
- Continuous integration checks have been made stricter by setting
biocCheck() error level to "error" rather than "never", and R CMD
check error level to "warning" rather than "error".
- Static analysis via the lintr package is now performed on each push
and PR. It should mostly be invisible to users but might result in
slightly increased performance in some cases.
- The superseded httr dependency has been replaced with the lighter
curl package, thus reducing the total number of dependencies for the
package from 42 to 40.
- The unused stringr dependency has been removed, reducing the total
number of dependencies for the package from 40 to 38.
- A minor PROTECT()/UNPROTECT() imbalance in the C code, exposed by
rchk, has been fixed. It is not likely to cause problems in
real-world situations but it could theoretically lead to crashes in
some cases.
- Argument path in internal function read_array_metadata() has been
renamed to zarr_path for consistency with other internal functions
- Some internal functions have been renamed with a leading dot, in
line with the officially recommended style for Bioconductor
packages.
- This package now uses testthat instead of tinytest as a testing
framework. This comes with more utilities to handle snapshot tests
and mocked tests.
- Function calls are now counted in tests to ensure we don't
repeatedly perform a task (in particular, an expensive I/O task)
more often than necessary.

[RCM](/packages/RCM)
---

                       Changes in version 1.25.2                        

- Fixed issues with citation file

                       Changes in version 1.25.1                        

- Correct method dispatching for buildConfMat
- Addressing some BiocCheck notes

[RCy3](/packages/RCy3)
----

                       Changes in version 2.30.0                        

- fix clearNodePropertyBypass, #229

[recount3](/packages/recount3)
--------

                       Changes in version 1.19.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- We have verified that the data hosted by IDIES at JHU is accessible
via either https://idies.jhu.edu/recount3/data or
https://data.idies.jhu.edu/recount3/data/. See
https://github.com/LieberInstitute/recount3/issues/58 for the
details.

                       Changes in version 1.19.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- recount3 no longer suggests using interactiveDisplayBase::display()
as that package will be deprecated on bioc-devel. While the
interactive sub-setting of a table was a nice little feature, it was
not an important feature. Related to
https://github.com/LieberInstitute/recount3/issues/57.

                       Changes in version 1.19.1                        

BUG FIXES

- Fixed a bug in file_retrieve() as noted at
https://github.com/LieberInstitute/recount3/issues/56. The bug and
solution to it was reported by @LiNk-NY.

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.54.0                        

CHANGES

- Disabled dataset timestamps, which were applied for compound
  datasets but not other types.  Back ported to 2.52.1.  (Thanks
  to Luke Zappia for reporting this
  https://github.com/Huber-group-EMBL/rhdf5/issues/160)

- Exposed HDF5 API functions H5Oget_info(), H5Fget_intent()
  H5Dget_num_chunks().

- Added function h5readTimestamps to extract object timestamp
  metadata.

- The previously deprecated `cset` argument in
  h5createAttribute() has been removed entirely. Please use the
  `encoding` argument instead.

[rhinotypeR](/packages/rhinotypeR)
----------

               Changes in version 2025-10-21 (2025-10-21)               

Fixed

- Example and test errors related to internal helpers
(compareAndColorSequences(), compute_cache()) now handled safely
using \dontrun{} and updated documentation.
- Removed outdated License stub issue and corrected DESCRIPTION
metadata formatting.
- Corrected internal tests referencing deprecated
testthat::local_null_device().

Improved

- CountSNPs() re-implementation: directly computes pairwise SNP
counts
without relying on model-based rounding (e.g., from MSA2dist). No
rounding artifacts

               Changes in version 2025-10-01 (2025-10-01)               

Added

- New function: alignToRefs(), which aligns user sequences to the
packaged rhinovirus prototype references using msa::msa() (ClustalW,
ClustalOmega, or MUSCLE), with an option to trim alignments to the
non-gap span of a chosen reference.

Improved

- SNPeek() and plotAA(): now support zooming to specific genomic
regions and highlighting individual sequences of interest.
- assignTypes(): now reports the genetic distance to the nearest
prototype even when the query is "unassigned".
- Distance calculation: switched to using MSA2dist for pairwise
distance computation. This provides broader model support — access
to all substitution models in ape::dist.dna, in addition to the
"IUPAC" model.
- Data access: prototype and example datasets are now exported as
packaged objects (rhinovirusVP4, rhinovirusPrototypesVP4) instead of
being accessed via system.file(), simplifying workflows.
- Improved documentation:
- Improved explanations and runnable examples added where
feasible.
- Visualization functions return results invisibly to avoid
console clutter, while still allowing advanced users to capture
outputs.
- Code style and formatting improved (e.g.,consistent internal
helpers).

[Rhtslib](/packages/Rhtslib)
-------

                        Changes in version 3.6.0                        

- No significant changes in this version.

[RNAdecay](/packages/RNAdecay)
--------

                       Changes in version 1.28.1                        

- maintenance - bugs introduced from ggplot2 4.0.0 update fixed.

[RnBeads](/packages/RnBeads)
-------

                       Changes in version 2.27.1                        

- Fixed issue with nv probes. Now they are filtered during the
  preprocessing.

[rrvgo](/packages/rrvgo)
-----

                       Changes in version 1.21.3                        

- Add new parameter to keep in the output of `calculateSimMatrix()`
  terms without Information Content (IC) (reported by GitHub user
  savagedude3. Thank you!)

- Fixed warning calling `calculateSimMatrix()` with default value for
  the `semdata` parameter (accurately detailed by GitHub user lahuuki.
  Thank you!)

                       Changes in version 1.21.1                        

- Add warning when using children=TRUE in `reduceSimMatrix()` and no
  keytype named GOALL in the Orgdb object

[rWikiPathways](/packages/rWikiPathways)
-------------

                       Changes in version 1.29.1                        

- README: cite latest WikiPathways paper

- README: updated a few links (including a relative link)

- README: removed a link to an unmaintained file

- removed .travis.yml

[S4Arrays](/packages/S4Arrays)
--------

                       Changes in version 1.10.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Drop dependency on the crayon package.

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.48.0                        

NEW FEATURES

- Add C callable new_SortedByQueryHits_from_IntAEAE.

SIGNIFICANT USER-VISIBLE CHANGES

- Replace C callable new_DataFrame with new_DFrame.

- Convert Rnw slides (rendered as PDF) to Rmd slides (rendered as
  HTML) (#128). By Marcel Ramos Pérez.

[SAIGEgds](/packages/SAIGEgds)
--------

                        Changes in version 2.9.2                        

- improve the work load balancing in `seqAssocGLMM_SPA()` according to
  SeqArray_1.49.6

- new option 'verbose.pval' to improve the display in
  `seqAssocGLMM_SPA()`

[SanityR](/packages/SanityR)
-------

- Initial Bioconductor submission.

[scafari](/packages/scafari)
-------

- Initial Bioc submission

[scDblFinder](/packages/scDblFinder)
-----------

                 Changes in version 1.23.1 (2025-07-17)                 

- fixed compatibility issues with xgboost version>=3

[scDiagnostics](/packages/scDiagnostics)
-------------

                        Changes in version 1.4.0                        

- Add functionality to process SCE objects for PCA computation via
the
new processPCA function.
- Add functionality to downsample SCE objects for diagnostic
functions.
- Add new diagnostic functions (calculateTopLoadingGeneShifts,
compareMarkers and calculateMMDPValue).
- Add graph integration diagnostic function in replacement of nearest
neighbor diagnostic, calculateGraphIntegration.
- Improve regressPC function and plot method, which can now also
regress against cell types and batches.
- Improve normalization for plotMarkerExpression diagnostic function.
- Improve user control and options for plot methods.
- Update vignettes to reflect all new changes.

[scDotPlot](/packages/scDotPlot)
---------

                        Changes in version 1.2.1                        

- Updates for new major version of ggplot2

[scGraphVerse](/packages/scGraphVerse)
------------

o Initial Bioconductor Submission

o Project Funding

- Supported by the National Centre for HPC, Big Data and Quantum
Computing under the European Union – Next Generation EU – CN00000013
(CUP: B93C22000620006).

[scider](/packages/scider)
------

                 Changes in version 1.8.0 (2025-10-30)                  

- Update help pages
- New function findNbrsGrid() to construct a neighbour list from grid
coordinates
- New function findNbrsSNN() to construct a SNN neighbour list from
assay
- New function findNbrsSpatial() to construct a distance-based
neighbour list from cell coordinates
- New function getNiche() to build a niche assay based on the profile
of neighbouring cells
- New function getClusters() to cluster cells in spe using graph
methods
- New functions globalMoran() and localMoran() to calculate Moran's I
statistics
- New function normalizeAssay() to perform spatial data normalization
- New function plotDR() for plotting reduced dimensions
- New functions plotImage() and plotLISA() for visualization
- New function realignVisium() to scale and straighten out Visium
coordinates
- New functions runPCA() and runUMAP() for performing dimensionality
reduction

[scLANE](/packages/scLANE)
------

- Initial Bioconductor Submission

[scp](/packages/scp)
---

                        Changes in version 1.19                         

scp 1.19.5

- Update Seqinfo-holding serialized object(s) to reflect new location
of Seqinfo class definition (the class was recently moved from
GenomeinfoDb to the new Seqinfo package). Hervé Pagès

scp 1.19.4

- Update scplainer reference.

scp 1.19.3

- docs: explicitly mention scplainer in the man pages

scp 1.19.2

- feat: improved error message when pattern is not found in
computeSCR
- feat: use typeMetadata (based on QFeatues 1.19.1)

scp 1.19.1

- fix: moved from the deprecated longFormat() to longForm()
- doc: fixed readSingleCellExperiment example (see #94)
- doc: minor documentation recompilation

scp 1.19.0

- New Bioconductor 3.22 release

[scQTLtools](/packages/scQTLtools)
----------

                 Changes in version 1.1.4 (2025-05-29)                  

- Updated vignettes and README.

                 Changes in version 1.1.3 (2025-05-26)                  

- Updated dataset names.
- Updated code comments and documentation.

[scrapper](/packages/scrapper)
--------

                        Changes in version 1.4.0                        

- Multiple changes to scoreMarkers():
  
  • Setting all.pairwise= to an integer will now return the top
  markers by effect size for each pairwise comparison. This
  is more efficient than creating the 3D array of effects and
  then ranking the genes in each comparison.
  
  • The min-rank summary is now limited to the top
  min.rank.limit= genes, for efficiency. All other genes will
  have their reported min-rank set to a maximum placeholder.
  
  • Added compute.group.mean= and compute.group.detected=
  options to disable reporting of the group-wise means.

- Multiple changes to testEnrichment():
  
  • Setting universe=NULL will automatically define the
  universe as the union of all genes.
  
  • Return value is now a data frame with the number of
  overlapping genes and the size of each gene set.

- Updated correctMnn() to the new algorithm used by the
  mnncorrect C++ library. This aims to improve correction
  accuracy and reduce sensitivity to the choice of reference.

- Bugfix to runAllNeighborSteps() to return an empty list when no
  step is requested.

- Set bound=0 by default in chooseHighlyVariableGenes(), for more
  convenient use with residuals.

- Modified scoreGeneSet() to return weights as a data frame
  containing the row indices of the genes in the set.

- Support character and logical vectors in sets= for
  aggregateAcrossGenes().

- Automatically remove empty clusters from the output of
  clusterKmeans().

- Supported blocking and related arguments for weighted averages
  of distances in scaleByNeighbors().

[scRepertoire](/packages/scRepertoire)
------------

                        Changes in version 2.5.7                        

BUG FIXES

- Fixed chain handling for BCR genes in percentGeneUsage() and
propagates to wrappers: percentGenes, percentVJ and vizGenes

                        Changes in version 2.5.6                        

BUG FIXES

- Fixed passing group.by for combineBCR()

                        Changes in version 2.5.5                        

- Update unit tests for ggplot2 v4

                        Changes in version 2.5.3                        

BUG FIXES

- Fixed clonalProportion calculation to use grouping properly during
combineExpression()
- Fixed immunarch support for exportClones() - TRA/Light chain column
handling

                        Changes in version 2.5.2                        

- Added support for mouse genes in quietBCRgenes() and
quietTCRgenes()

                        Changes in version 2.5.1                        

- Introduced pairwise calculations to StartracDiversity()
- Internal function conversion for clonalSizeDistribution() - remove
cubature, truncdist and VGAM from dependencies.
- Increased speed of clonalSizeDistribution()

                        Changes in version 2.5.0                        

UNDERLYING CHANGES

- Update/Improve code for loadContigs()
- Consolidated support for discrete AIRR formats under the umbrella
of
AIRR
- Added "tcrpheno" and immunarch to exportClones()
- Converted exportClones() to base R to reduce dependencies
- Added dandelionR vignette to pkgdown site
- Added tcrpheno vignette to pkgdown site
- percentAA() refactored to minimize dependencies and use immApex
calculateFrequency()
- positionalEntropy() refactored to minimize dependencies and use
immApex calculateEntropy()
- clonalDiversity() refactored for performance - it now calculates a
single diversity metric at a time and includes new estimators like
"gini", "d50", and supports hill numbers).
- percentKmer() refactored to use immApex calculateMotif for both aa
and nt sequences. No longer calculates all possible motifs, but only
motifs present.
- clonalCluster() now allows for dual-chain clustering, V/J
filtering,
normalized or straight edit distance calculations, and return of
clusters, igraph objects or adjacency matrix
- combineBCR() offers single/dual chain clustering, aa or nt
sequences, adaptive filtering of V and J genes and normalized or
straight edit distance calculations
- percentGeneUsage() now is the underlying function for
percentGenes(), percentVJ(), and vizGenes() and allows for percent,
proportion and raw count quantification.
- Added common theme (internal .themeRepertoire()) to all plots and
allow users to pass arguments to it

BUG FIXES

- clonalCompare() issue with plotting a 0 row data frame now errors
with message
- clonalScatter() group.by/axes call now works for non-single-cell
objects
- Fixed issue with NULL and "none" group.by in combineExpression()
- Allowing multi groupings via x.axis and group.by in
clonalDiversity()

[scRNAseqApp](/packages/scRNAseqApp)
-----------

                       Changes in version 1.9.26                        

- Fix the bug by sanitize the file names when load ATAC track.

                       Changes in version 1.9.25                        

- add sum,max to plot type of Sunburst plot panel.

                       Changes in version 1.9.24                        

- Add zoom in for plot.

                       Changes in version 1.9.23                        

- Fix a bug for file accessibility check.

                       Changes in version 1.9.22                        

- Add more color themes

                       Changes in version 1.9.21                        

- support for dark mode

                       Changes in version 1.9.20                        

- Add support for config.dcf file.

                       Changes in version 1.9.19                        

- Add close and delete button for issues.

                       Changes in version 1.9.18                        

- Update css for issues.

                       Changes in version 1.9.17                        

- Add css for issues.

                       Changes in version 1.9.16                        

- Add markdown support for issues.

                       Changes in version 1.9.15                        

- Add edit and reply button for issues.

                       Changes in version 1.9.14                        

- Fix the test when there is old files in the tmpfolder.

                       Changes in version 1.9.13                        

- Handle the file accessibility issue.

                       Changes in version 1.9.12                        

- Fix the bug if the script.js does not update.

                       Changes in version 1.9.11                        

- Add issues tab.

                       Changes in version 1.9.10                        

- Update the homepage to include datasets information.

                        Changes in version 1.9.9                        

- Update the unit test.

                        Changes in version 1.9.8                        

- Filter cells by a range of value.

                        Changes in version 1.9.7                        

- Allow single gene heatmap plot.

                        Changes in version 1.9.6                        

- Add progress bar for atac signal plot.

                        Changes in version 1.9.3                        

- Adjust the atac peak links plot space.

                        Changes in version 1.9.2                        

- Add edge links.

                        Changes in version 1.9.1                        

- Fix the issue for ATAC track plot if there is no 'type' in peak
  annotation
  and no 'score' in links.

[scTensor](/packages/scTensor)
--------

                       Changes in version 2.18.1                        

- A vignette modified.

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.50.0                        

NEW FEATURES

- new '.status_file' & '.proc_time' in `seqParallel()` for controlling
  the display of the status of child processes

- enhance workload balancing within `seqParallel()`

- new 'balancing=NA' in `seqAlleleFreq()`, `seqAlleleCount()`,
  `seqGetAF_AC_Missing()`, `seqSetFilterCond()`, `seqUnitFilterCond()`
  (using workload balancing by default)

- new '.balancing' in `seqApply()` and `seqBlockApply()`

- new 'alt' & 'ns' in `seqGetAF_AC_Missing()`

- new 'dosage' in `seqExampleFileName()`

UTILITIES

- `seqDigest()` allows a GDS file name in the first argument

- a new option 'parallel' in `seqDigest()`

- `seqVCF2GDS()` and `seqBED2GDS()`: save the gds node 'chromosome'
  before 'position'

- `seqAddValue(, varnm="annotation/filter", desp=...)` adds
  'Description'
  via the argument 'desp'

- tweak the display of `seqGDS2BED()`

- `seqBlockApply()` and `seqApply()` allows a GDS file name in the
  first
  argument

- `seqGDS2VCF()` generates an indexing file (.tbi) according to the
  vcf.gz output file, requiring the Rsamtools package

- the first argument of `seqOptimize()` can be a SeqVarGDSClass object

BUG FIXES

- Previous `seqGetAF_AC_Missing(, minor=FALSE)` returns the AF for the
  alt. alleles and the AC for the ref. alleles, when a dosage GDS file
  is
  used; a new 'alt' is used to specify the ref. or alt. alleles.

- Fix '##INFO=<>' in the output of `seqGDS2VCF()` when the INFO
  variable is
  added by users

[Seqinfo](/packages/Seqinfo)
-------

                        Changes in version 1.0.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The Seqinfo class was moved from the GenomeInfoDb package to the
  new Seqinfo package.
  IMPORTANT WARNING: This means that:
  - Seqinfo-holding S4 objects that were serialized before this change
  (i.e. prior to Bioconductor 3.22) are no longer valid objects in
  Bioconductor >= 3.22. However they can be fixed with updateObject().
  - This also means that Seqinfo-holding S4 objects serialized with
  Bioconductor >= 3.22 are no longer valid objects in Bioconductor
  < 3.22. These objects can **not** be fixed with updateObject().

[SETA](/packages/SETA)
----

- Initial Bioconductor Release

                 Changes in version 0.99.3 (2025-09-15)                 

Bioconductor Preparation

- Updated package for Bioconductor review
- Added package-level documentation (?SETA)
- Bioconductor installation instructions
- Added proper error handling for optional dependencies
- Created NEWS file for change tracking

Documentation Improvements

- Added detailed motivation and comparison sections to introductory
vignette
- Vignette descriptions with step-by-step explanations
- Added backticks for code formatting throughout vignettes
- Added BiocStyle dependency for proper vignette formatting

Package Structure

- Created R/SETA-package.R with roxygen2 documentation
- Added proper error handling for SeuratObject and
SingleCellExperiment dependencies
- Updated DESCRIPTION with BiocStyle dependency
- package-level help with overview

                 Changes in version 0.99.1 (2025-08-26)                 

New Features

- Initial submission of SETA package for compositional analysis of
single-cell RNA-seq data
- Added compositional transform functions (CLR, ALR, ILR, balance)
- Implemented latent space analysis methods (PCA, PCoA, NMDS)
- Created multi-resolution analysis capabilities with hierarchical
taxonomies
- Added distance calculation functions for compositional data
- Included vignettes with educational content

Functions Added

- setaCounts(): Extract cell-type count matrices from single-cell
objects
- setaTransform(): Apply various compositional transforms
- setaCLR(): Centered log-ratio transformation
- setaALR(): Additive log-ratio transformation
- setaILR(): Isometric log-ratio transformation
- setaBalance(): User-defined balance transformations
- setaPercent(): Percentage transformation
- setaLogCPM(): Log counts-per-million transformation
- setaLatent(): Dimensionality reduction methods
- setaDistances(): Calculate compositional distances
- setaTaxonomyDF(): Create hierarchical taxonomies
- taxonomy_to_tbl_graph(): Convert taxonomies to graph objects
- setaMetadata(): Extract sample-level metadata
- resolveGroup(): Resolve group specifications for balance transforms

Documentation

- Added package documentation
- Created three detailed vignettes:
- Introductory vignette with basic workflow
- Comparing samples vignette with statistical analysis
- Reference frames vignette with multi-resolution analysis
- Added function-level documentation with examples
- Included mathematical foundations and references

Data

- Added mock data functions for testing and examples:
- mockSeurat(): Mock Seurat object
- mockSCE(): Mock SingleCellExperiment object
- mockCount(): Mock count matrix
- mockLong(): Mock long-form data frame

Dependencies

- Core dependencies: dplyr, MASS, Matrix, stats, tidygraph, rlang
- Suggested dependencies: BiocStyle, caret, corrplot, ggplot2,
ggraph,
knitr, methods, patchwork, reshape2, rmarkdown, SeuratObject,
Seurat, SummarizedExperiment, SingleCellExperiment,
TabulaMurisSenisData, tidyr, tidytext, testthat

Bioconductor Integration

- Prepared package for Bioconductor submission
- Added BiocStyle dependency for vignette formatting
- Updated README with Bioconductor installation instructions
- Created package-level documentation for ?SETA
- Added NEWS file for change tracking

Testing

- Test suite covering all major functions
- Tests for error handling and edge cases
- Validation of mathematical correctness for transforms
- Mock data testing for integration with single-cell objects

[shinybiocloader](/packages/shinybiocloader)
---------------

                        Changes in version 1.0.0                        

- Initial Bioconductor submission

[shinyDSP](/packages/shinyDSP)
--------

                        Changes in version 1.1.2                        

Bug fixes for peer review

- in the setup nav panel, count and annotation tables actually show
the first 100 rows.
- added custom BH-adjusted P value input.
- removed 'turbo' as an option for heatmap colormap because it's not
a
part of hcl.colors(). Instead made all 115 palettes in hcl.colors()
available for selection.
- updated the secondary vignette to show RUV4 normalization,
differential gene expression analysis, plotting Volcano and
heatmaps, and code chunk execution times.

[signeR](/packages/signeR)
------

                       Changes in version 2.11.2                        

- Move VariantAnnotation from Depends to Imports

                       Changes in version 2.11.1                        

- Remove survivalAnalysis dependency

- Fixes for some ggplot2 warnings

- Update C++ req to C++14

[simPIC](/packages/simPIC)
------

                 Changes in version 1.5.3 (2025-04-28)                  

- Extended simPICcountclass object to include: nGroups, batch,
differential accessibility, and bcv parameters.
- Updated estimateBCV function in simPICestimate.
- Updates in simPICsimulate:
- simPICsimBatchEffects
- simPICsimBatchCellMeans
- simPICsimulatemultiDA
- simPICsimulateBCVmeans
- simPICsimulateTrueCountsGroups
- Updated vignette to demonstrate new functionality--- multiple
cell-types
- Removed redundancy in citation in DESCRIPTION

                 Changes in version 1.5.2 (2025-04-21)                  

- Added package logo to vignette.

                 Changes in version 1.5.1 (2025-04-21)                  

- Fixing typo in vignette, changed lognormal to lognormal-gamma

[SingleCellAlleleExperiment](/packages/SingleCellAlleleExperiment)
--------------------------

                        Changes in version 1.5.1                        

Adressing 3.22 devel error for failing test after DropletUtils update
Minor changes regarding Links set in documentation

[SingleCellSignalR](/packages/SingleCellSignalR)
-----------------

                  Changes in version 2.0 (2024-09-17)                   

- Major Code redesign : Major update based on the use of S4
  classes from BulkSignalR.

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.18.1 (2025-07-01)                 

- Updated enrichR examples

[sketchR](/packages/sketchR)
-------

                        Changes in version 1.5.2                        

- Add vignette with workflow examples

[SmartPhos](/packages/SmartPhos)
---------

 - Initial Bioconductor Submission

                 Changes in version 0.99.0 (2025-02-24)                 

- SmartPhos package with shinyApp SmartPhos Explorer

- create MultiAssayExperiment object from Spectronaut & MaxQuant files

- different visualization and analysis available.

- possible analysis:

- data transformation and normalization

- data imputation

- batch effect correction

- principal component analyis

- differential expression analyis

- time-series clustering

- geneset enrichment analyis

- phopho-signature enrichment

- kinase activity inference

- visualizations:

- intensity plot

- pca plot

- heatmaps

- volacano plot for differential expression analysis

- p value histogram

- time series plot

- geneset enrichment plots

- kinase activity inference

[smoothclust](/packages/smoothclust)
-----------

                 Changes in version 1.5.7 (2025-09-15)                  

- optimizations to speed up runtime

[smoppix](/packages/smoppix)
-------

                        Changes in version 1.1.9                        

- Added citation file to refer to the preprint directly
- Make sure vignettes are build upon installation and mentioned in
readme

                        Changes in version 1.1.8                        

- Export some convenience functions for other packages

                        Changes in version 1.1.7                        

- Remove memory intensive looping over features in data preparation
for mixed-effect modelling
- Depend on Rfast package for fast rowsorting
- Include tikzpicture as schematic overview in readme

                        Changes in version 1.1.6                        

- Removing functionality related to Moran’s I
- Speed up and reduce memory requirements of linear (mixed) model
fitting by preparing design matrices only once
- Correct bug in weighting for cell-wise nn distributions

                        Changes in version 1.1.5                        

- Export loadBalanceBplappy for use in other packages
- Mention Bioconductor installer in README

                        Changes in version 1.1.4                        

- Optimize addWeightFunction

                        Changes in version 1.1.3                        

- Some speedups in the buildDataFrame function
- With reference to bioRxiv preprint

                        Changes in version 1.1.2                        

- Calculate nnPair PIs only between desired features

                        Changes in version 1.1.1                        

- Allow to pass own colours to plotExplore

[SNPRelate](/packages/SNPRelate)
---------

                       Changes in version 1.44.0                        

UTILITIES

- tweak the display (by showing # of variants) when calculating the
  allele
  counts and frequencies of a SeqArray GDS file

- tweak the display of date

- import RhpcBLASctl to control the number of threads in BLAS and
  LAPACK
  in `snpgdsPCA()`

- set the threshold 'missing.rate=0.01' by default in
  `snpgdsIBDMoM()`, `snpgdsIBDMLE()`, `snpgdsIBDKING()`,
  `snpgdsDiss()`,
  `snpgdsGRM()`, `snpgdsFst()`, `snpgdsIndivBeta()`, `snpgdsIBS()`,
  `snpgdsIBSNum()`, `snpgdsLDpruning()`, `snpgdsPCA()`,
  `snpgdsEIGMIX()`

[sosta](/packages/sosta)
-----

                        Changes in version 1.1.4                        

- the argument imageCol is now not longer needed in reconstruction
functions
- assingCellsToStructures for large datasets and 1 sample

                        Changes in version 1.1.3                        

- added sticker and preprint reference

                        Changes in version 1.1.2                        

- switch from GenomeInfoDb to new Seqinfo package

                        Changes in version 1.1.1                        

- added new function minCellTypeStructDist
- added possibility to reconstruct everything around an object
(complement = TRUE)
- relaxed assertion of having at least two cells of each cell type in
reconstruction to having at least one of one cell type

                        Changes in version 1.0.1                        

- fix problems in when assigning cell to structures

[SpaceMarkers](/packages/SpaceMarkers)
------------

                        Changes in version 2.0.0                        

- Added support for directed cell-cell interaction (see
calculate_gene_scores_directed, calculate_influence)
- Functions for supporting ligand-receptor interactions based on
directed cell-cell interactions
- Added ligand–receptor and gene-set scoring utilities:
calculate_lr_scores, calculate_gene_set_score, and
calculate_gene_set_specificity.
- Added support for handling Visium HD
- plotting functions using circlize for showing ligand-receiver
interactions.
- API and naming standardization (backwards-incompatible) Major
public
functions renamed to snake_case for consistency (examples:
calcInfluence → calculate_influence, getInteractingGenes →
get_interacting_genes, getPairwiseInteractingGenes →
get_pairwise_interacting_genes, getSpatialFeatures →
get_spatial_features, getSpatialParameters → get_spatial_parameters)

 
[SpaceTrooper](/packages/SpaceTrooper)
------------

- Initial Bioconductor Submission

[spARI](/packages/spARI)
-----

- Submitted to Bioconductor

[SparseArray](/packages/SparseArray)
-----------

                       Changes in version 1.10.0                        

NEW FEATURES

- Add the two following capabilities to subassignment of
  SVT_SparseArray
  objects:
  - The right value (a.k.a. replacement value) now can be an ordinary
  array.
  - Using a long linear subscript (i.e. long integer or numeric vector)
  is now supported. Note that this is also enabled for NaArray objects.

[spatialFDA](/packages/spatialFDA)
----------

                       Changes in version 1.1.21                        

- allowing for the explicit passing of correction in spatialInference
and crossSpatialInference to fix issues with calculating Lcross.
- passing messages with verbose logical for better usage
- bug fix in tests with removed subsetby argument

                       Changes in version 1.1.20                        

- only returning the transformed curves in spatialInference plus
removing the assays and the rowData for faster computations.

                       Changes in version 1.1.19                        

- fixing plotFbPlot ylims across conditions to have the same range-

                       Changes in version 1.1.18                        

- removed the argument subsetby from spatialInference and
crossSpatialInference as this will always be the same as image_id in
this function
- indicated in the documentation that ... in spatialInference and
crossSpatialInference is passed to both spatstat.explore functions
and to refund::pffr.

                       Changes in version 1.1.17                        

- added an assertion that the maximum radius considered must always
be
smaller than the maximum window length of each image.
- added a test for this case
- add a heuristic to check for a convenient and good rMax to consider
- added heuristic to the vignettes

                       Changes in version 1.1.16                        

- adjusted that all failed calcMetricPerFov output NA in the
dimension
of rSeq, clearly stating which values in metricRes where NA.
- return the raw non-filtered metricResRaw in spatialInference as
well
as the filtered metricRes for comparison.
- implemented option to add different QC metrics of the functional
model and plot them with different shapes in plotCrossHeatmap
- adapted the overview vignette to plot the raw spatial statistics
curves
- accelerated all examples and tests with functional GAM via
algorithm
= "bam"

                       Changes in version 1.1.15                        

- added a test to check the edf values and order from the RSE
calculation. Made the RSE code more robust.
- added progressbar in crossSpatialInference.
- adapted the overview vignette to reflect the recent changes.

                       Changes in version 1.1.14                        

- change the sum of mean residuals per condition to be the residual
standard error. The main reason is to have a more comparable
estimator of model fit by dividing by the degrees of freedom per
condition. The degrees of freedom are the no. of data points per
condition (no. of curves * data points per curve)

- the sum of the estimated degrees of freedom of the parameters per
condition as provided by the summary.pffr() output.

                       Changes in version 1.1.13                        

- provide the condition-wise sum of mean residuals across the curves.
This is a measure to judge the quality of the fitted functional GAM
by condition.
- allow for an automated threshold selection delta in
spatialInference
based on the mean over all images of the minimum nearest neighbour
distance between points. Added test for this use case

                       Changes in version 1.1.12                        

- changed the thresholding in calcMetricPerFov from having at least 1
point per (sub)point pattern to at least two in order to avoid
having images with way to uncertain curves. In order to accommodate
for this change some tests and examples had to be adapted
- spatialInference now prints the adjusted R squared to quickly
assess
overall model fit. As this is just one summary, still rigorous
assessment via qq plots and autocorrelation should be performed.

                       Changes in version 1.1.11                        

- Change in the weighting procedure in spatialInference and no
calculation of weights in calcMetricPerFov. It is now possible to
define custom weights, by equal weights, by the total number points,
and for multitype processes by the smaller point pattern (min) and
the larger point pattern (max)
- Option to extract the mean functional coefficient of the functional
GAM over the domain $r$ as a measure of the effect size in
extractCrossInferenceData
- plotCrossHeatmap creates now a bubble plot with the effect size
being the colour and the size of the bubble being the p-value

                       Changes in version 1.1.10                        

- p-value adjustment option in plotCrossHeatmap and filtering of
coefficients for the heatmap. Fix some small global variable
defintion errors

                        Changes in version 1.1.9                        

- enabling Fisher's variance-stabilising transformation in
spatialInference as it is recommended for $G$ functions.

                        Changes in version 1.1.8                        

- adding a constant to the log in plotCrossHeatmap to avoid NULL
values being plotted

                        Changes in version 1.1.7                        

- Fixing error in plotCrossHeatmap when the model is NULL.

                        Changes in version 1.1.6                        

- Error handling in spatialInference when one condition has no images
with curves.
- Wrote as well a new test for this case

                        Changes in version 1.1.5                        

- New overview vignette, showing the usability of spatialInference
- Plotting of a heatmap from extracted from crossSpatialInference
- Adding of the package sticker to the README

                        Changes in version 1.1.4                        

- functionalGAM accepts flexible naming of the intercept column.

                        Changes in version 1.1.3                        

- changed the default of the GAM estimation to be a Gaussian with log
link for positivty of the response

                        Changes in version 1.1.2                        

- wrote a new convenience function crossSpatialInference which makes
the estimation across a range of cell types easier.

                        Changes in version 1.1.1                        

- wrote a new convenience function spatialInference which makes the
FDA estimation of the spatial statistics functions simpler
- in order to do this various changes to prepData were introduced
- furthermore, small syntactic changes for the retrieval of the
functional intercept were added to functionalGam and plotMdl
- The vignette was simplified to reflect the changes in prepData

[Spectra](/packages/Spectra)
-------

                        Changes in version 1.19                         

Change 1.19.11

- Fix in longForm() method for MsBackendMemory if only peaks
variables
are extracted.

Change 1.19.10

- Add default implementation of methods for the MsBackend class
(issue
#338). This reduced code redundancy in Spectra and simplifies
implementation of new MsBackend classes considerably.

Change 1.19.9

- Add implementations of spectraVariableMapping() and
spectraVariableMapping<- for MsBackend and Spectra objects.

Change 1.19.8

- Add default implementations for MsBackend for replacement methods
centroided<-, collisionEnergy<-, dataOrigin<-,
isolationWindowLowerMz<-, isolationWindowTargetMz<-,
isolationWindowUpperMz<-, msLevel<-, polarity<-, precursorMz<-,
rtime<- and smoothed<-. This reduces the amount of methods required
to be implemented for new MsBackend classes.

Change 1.19.7

- Add longForm,MsBackend and longForm,Spectra methods (issue #360).

Change 1.19.6

- Add a precursorPurity() function.

Change 1.19.5

- Addition of fragmentGroupIndex() function. This function generates
a
integer index grouping MSn (MS level > 1) spectra with their
corresponding (unique) MS1 spectra based on acquisition order.
- Fix on cbind2() so that the output is a Spectra object and not the
backend.

Change 1.19.4

- Fix passing of BPPARAM parameter in peaksData() method for Spectra.
Now is being passed down to .peaksapply().

Change in 1.19.3

- Fix in estimatePrecursorMz() that would ignore parameter BPPARAM.

Change in 1.19.2

- Fix export of data using MsBackendMzR: save also the precursor scan
number for MSn data (spectra variable "precScanNum").

Change in 1.19.1

- Add precScanNum() method for MsBackendCached.

[SpectraQL](/packages/SpectraQL)
---------

                         Changes in version 1.3                         

Changes in 1.3.1

- Add citation.

[SpectriPy](/packages/SpectriPy)
---------

- Initial Submission to Bioconductor 

[SPICEY](/packages/SPICEY)
------

- Initial Submission to Bioconductor

[STADyUM](/packages/STADyUM)
-------

Initial submission to BioConductor.

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

[stPipe](/packages/stPipe)
------

 - Bioconductor submission.

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.40.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Depends on the new Seqinfo package instead of GenomeInfoDb.

[SuperCellCyto](/packages/SuperCellCyto)
-------------

- Initial submission to Bioconductor


[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.21.2                        

- Deprecated consensus clustering removed from ExoLabel
- ExoLabel now supports argument header=c(FALSE, TRUE, integer(1L))
to
skip the first 0,1,n lines of each file (respectively)
- ExoLabel now supports gzip-compressed files as input
- Some examples in man pages have been made faster

                       Changes in version 1.21.1                        

- Internal code reorganization for ExoLabel
- Further arguments to SuperTree now correctly passthrough to
Treeline
- Some deprecated arguments to EvoWeaver have been removed

                       Changes in version 1.21.0                        

- First development version of Bioconductor 3.22

[TargetSearch](/packages/TargetSearch)
------------

                       Changes in version 2.12.0                        

NEW FEATURES

- Return intensities and masses as attributes in the output of
  `RIcorrect` and `RImatrix`.

BUG FIXES

- `RIcorrect` cannot open file under certain conditions. This occurs
  when
  the input files are NC4 and the `writeCDF4path` option is set to
  write
  to a different path. Same if they are TargetSearch files, but lack
  the
  correct extension.

INTERNAL

- Add new tests for functions `RIcorrect` and `RImatrix`.

[TaxSEA](/packages/TaxSEA)
------

                 Changes in version 1.1.8 (2025-08-12)                  

- Added output for all taxon sets
- Fixed spelling in GutMGene_producers_of_Phenyalanine

                 Changes in version 1.1.7 (2025-08-11)                  

- Added Gut-Brain modules database from Valles-Colomer et al. Nature
Microbiology. 2019.

                 Changes in version 1.1.6 (2025-08-08)                  

- Bug fixes

                 Changes in version 1.1.5 (2025-08-06)                  

- Fixed bug in BacDive Database

                 Changes in version 1.1.4 (2025-08-05)                  

- Added BacDive Database
- Added Blossum and VANISH taxon sets in disease associations
category

[TCGAutils](/packages/TCGAutils)
---------

                       Changes in version 1.30.0                        

Significant User-visible changes

- Deprecated mirbase.db package affects mirToRanges function.

Bug fixes and minor improvements

- Use BiocBaseUtils::checkInstalled to check for suggested packages.

[TENET](/packages/TENET)
-----

                        Changes in version 1.2.0                        

- Add options to generate violin plots instead of boxplots

- Add support for Jenks natural breaks and custom grouping to
step7TopGenesSurvival

- Add support for custom regions and MotifDb filtering options to
step7LinkedDNAMethylationSitesMotifSearching

- step7TopGenesUserPeakOverlap now accepts peak names, generates
results for each individual peak rather than each peak file, and
accepts multiple input objects or directories

- External genomic region files can now be provided as a list of
multiple files and/or directories containing them, not just a
directory

- Fix TENETSavedSizePlot bug that caused warnings and failure to
respect an explicitly specified size

- Fix TENETCacheAllData to fail if any dataset does not exist on the
hub

- Significant documentation improvements

- Do not suggest installing the development version of Bioconductor
to
use the stable version of the package, since it is no longer
necessary

- Update installation instructions for Ubuntu 24.04 and add a missing
required package

- Add a dependency on ggplot2 4.0 because the example
MultiAssayExperiment object now requires it

- Add CITATION file

                        Changes in version 1.1.1                        

- Merge survival bug fix from Bioconductor 3.21 release branch; see
below

                        Changes in version 1.1.0                        

- Release generated automatically by Bioconductor - no changes

[TFEA.ChIP](/packages/TFEA.ChIP)
---------

                       Changes in version 1.29.4                        

New Features

- Added a meta-analysis function to combine results across datasets.

- Introduced a new wrapper function that performs the full enrichment
  workflow in one step, reducing manual intervention.

Enhancements

- Refactored and optimized several core functions for improved
  performance.

- Improved plotting functions for clearer and more customizable
  visualizations.

Data Updates

- Updated the ChIPDB object using newly available regulatory
  region-to-gene annotations from the ENCODE rE2G method and CREdb,
  enhancing the biological relevance and accuracy of TF-gene
  associations.

- Integrated new datasets from ReMap.

New default database

- The TF-gene database bundled with TFEA.ChIP was built using the ReMap
  2022 ChIP-seq collection and regulatory links from the ENCODE rE2G
  method. Due to memory constraints, the internal database includes
  only a subset of the 8,000+ available ChIP-seq experiments.
  Specifically, we selected 926 experiments performed in the ENCODE
  Project’s Common Cell Types, prioritizing broadly representative and
  high-quality datasets.

To download the full database, as well as other ready-to-use databases generated for TFEA.ChIP, visit

- <https://github.com/yberda/ChIPDBData>

- <https://github.com/LauraPS1/TFEA.ChIP_downloads>

[topdownr](/packages/topdownr)
--------

                        Changes in version 1.31                         

Changes in version 1.31.1

- Move package to codeberg.org.

[topGO](/packages/topGO)
-----

                       Changes in version 2.62.0                        

- The vignette of topGO is now sporting its new RMarkdown dress,
which
ideally will make it easier to navigate within the modern browsers
- The documentation of the package is now leveraging the roxygen
framework, which should simplify the process of editing it in-sync
with the code
- On the sidelines of the whole development process, Federico has
been
added as a co-maintainer of the package, joining Adrian (the
original creator of the package) to possibly propel the creation of
additional new features/improvement of the performance

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.45.2                        

- Use yaxis to control draw or not for the interaction legend.

- Fix a bug in reading mcool file.

                       Changes in version 1.45.1                        

- Fix a bug in reading mcool file.

[transmogR](/packages/transmogR)
---------

                        Changes in version 1.5.3                        

New Features

- Exported the previously internal function overdispFromBoots() for
standalone importing of bootstraps

                        Changes in version 1.5.2                        

Bug Fixes

- Switched from ComplexUpset to SimpleUpset as dependency

[treeclimbR](/packages/treeclimbR)
----------

                        Changes in version 1.5.1                        

- Add GenomeInfoDb to Suggests


[TVTB](/packages/TVTB)
----

                 Changes in version 1.35.2 (2025-08-05)                 

Unit test fix

- Fix unit test to match new S7 structure in ggplot2.

[txcutr](/packages/txcutr)
------

                       Changes in version 1.15.2                        

NEW FEATURES

- None.

SIGNIFICANT USER-VISIBLE CHANGES

- Removed deprecation notice.
- Updated author email.

BUG FIXES

- Added GenomeInfoDbData to Suggests (Issue #23).

[txdbmaker](/packages/txdbmaker)
---------

                        Changes in version 1.6.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add Seqinfo to Depends and move GenomeInfoDb to Imports.

- Drop mirbase.db from Suggests (package got removed from BioC 3.22
  after being deprecated in BioC 3.21). As a consequence, utility
  function supportedMiRBaseBuildValues() and argument 'miRBaseBuild'
  are now defunct (the latter is an argument of various functions in
  the package).

DEPRECATED AND DEFUNCT

- Utility function supportedMiRBaseBuildValues() and argument
  'miRBaseBuild' are defunct (the latter is an argument of various
  functions in the package).

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.28.0                        

- Addition of makeLinkedTxpData() a lightweight version of
makeLinkedTxome().
- makeLinkedTxome() will now accept a digest and indexName, as an
alternative to pointing to a indexDir.
- Addition of three new functions. importData(), inspectDigests(),
updateMetadata() allowing import of oarfish (v0.9.0) quantification
files when annotated and novel reference transcript sets have been
mixed.
- TODO: As of tximeta v1.28.0, some more work is needed to resolve
how
the add*() and retrieve*() functions will work, as well as
summarizeToGene() for data imported with importData().

                       Changes in version 1.27.11                       

- Addition of makeLinkedTxpData() a lightweight version of
makeLinkedTxome().

                       Changes in version 1.27.6                        

- makeLinkedTxome() will now accept a digest and indexName, as an
alternative to pointing to a indexDir.

                       Changes in version 1.27.4                        

- Addition of three new functions. importData(), inspectDigests(),
updateMetadata() allowing import of oarfish (v0.9.0) quantification
files when annotated and novel reference transcript sets have been
mixed. TODO: As of tximeta v1.28.0, some more work to do is the
development of linkedTxpData and resolving how the add*() and
retrieve*() functions will work, as well as summarizeToGene().
- Major code reorganiziation, splitting out alevin-processing code
and
moving a number of metadata sub-routines to metadata_helpers.R.

                       Changes in version 1.27.2                        

- GENCODE 49 (H.s.), M38 (M.m), and Ensembl 115 (Sep 2025)

                       Changes in version 1.27.1                        

- GENCODE 48 (H.s.), M37 (M.m), and Ensembl 114 (May 2025)

[UCell](/packages/UCell)
-----

                       Changes in version 2.13.3                        

- New parameter missing_genes to handle missing genes from
  signatures (either 'impute' or 'skip')

- New demo to discuss important parameters of the algorithm

- Reformatted scoring function, using gene indices instead of
  string matches

[UCSC.utils](/packages/UCSC.utils)
----------

                        Changes in version 1.6.0                        

- No significant changes in this version.

[UniProt.ws](/packages/UniProt.ws)
----------

                       Changes in version 2.50.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- `query` input in `queryUniProt` can now be a named list; clarified in
  documentation

BUG FIXES AND MINOR IMPROVEMENTS

- `collapse` in `queryUniProt` does not need spaces padding, argument
  is not
  case-sensitive

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.28.0                        

BUG FIXES

- Quiet a couple of warnings triggered during testing the
  universalmotif
  initializer function. These warnings should never trigger in
  real-world
  scenarios.


[updateObject](/packages/updateObject)
------------

                       Changes in version 1.14.0                        

- No changes in this version.

[UPDhmm](/packages/UPDhmm)
------

                        Changes in version 1.5.2                        

- Fix vignette for LatEX error.

                        Changes in version 1.5.1                        

- Update vignette

                        Changes in version 1.5.0                        

- Added parallelization option BPPARAM in calculateEvents() and some
new columns in the output.
- Added new vignette: "Preprocessing Guide".
- Update User Guide vignette
- Add 3 new functions markRecurentRegions(),
identifyRecurrentRegions() and collapseEvents()

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

                       Changes in version 1.39.2                        

- update checks to handle voom(...,adaptive.span)

                       Changes in version 1.39.1                        

- version bump

[velociraptor](/packages/velociraptor)
------------

                       Changes in version 1.19.1                        

- Switch import of OS-detecting function to basilisk.
- Update scvelo to 0.3.3 (PyPI) for all operating systems.

[visiumStitched](/packages/visiumStitched)
--------------

                        Changes in version 1.1.1                        

NEW FEATURES

- build_SpatialExperiment() (and add_array_coords()) takes a new
parameter algorithm whose default is now "LSAP". The older
implementation is available via the "Euclidean" algorithm. The now
LSAP approach is a significant improvement in how new array
coordinates are assigned, preventing all duplicate and most empty
mappings within a single capture area, leading to improvements in
downstream applications like clustering.

[xcms](/packages/xcms)
----

                         Changes in version 4.7                         

Changes in version 4.7.3

- Fix conversion from XcmsExperiment objects to XCMSnExp (issue
#807).
- Add parameter columns to chromPeakData() to allow extraction of
selected columns.
- Add support for parameter skipFilled to featureSpectra() for
XcmsExperiment and XcmsExperimentHdf5 to avoid extraction of spectra
for gap-filled chromatographic peaks.
- Fix: add missing import of the spectraData() method from Spectra.
- Add parameter minMzWidthPpm to featureArea() and
ChromPeakAreaParam() allowing to define a minima guaranteed m/z
width of the feature areas respectively regions to integrate data
from for gap-filling (issue #756).

Changes in version 4.7.2

- Change default for [,XcmsExperiment and [,XcmsExperimentHdf5
keepAdjustedRtime from FALSE to TRUE.
- Fix issue in [,XcmsExperiment: the subset operation did not revert
the retention times of chromatographic peaks with parameter
keepAdjustedRtime = FALSE (issue #801).
- Fix issue in chromPeakSpectra for XcmsExperiment with duplicated
chrom peak IDs (issue #796).
- Fix issue in chromPeakSpectra() for XcmsExperiment that added
column
(sample) names to the resulting XChromatograms object with newer
versions of the MsExperiment package (>= 1.11.1); issue #803.

Changes in version 4.7.1

- Change the naming convention of chromatographic peaks to include
also the MS level after "CP": "CP1001" instead of "CP001" for a
chromatographic peak in MS level 1 or "CP2001" instead of "CP001"
for a chromatographic peak in MS level 2.
- Massively reduce main memory demand of xcms with the new on-disk
storage mode of the preprocessing results: the new
XcmsExperimentHdf5 result object stores all preprocessing results in
an HDF5 file.
- Optimization and performance improvements for extraction of
chromatographic data. This includes using MsCoreUtils::reduce().
- Restructure and clean-up of documentation.
- Don't export unnecessary get/set methods for Param classes.
- Throw an error if parameters mz or rt for chromatogram() contain
missing values.

[zellkonverter](/packages/zellkonverter)
-------------

                       Changes in version 1.20.0                        

New features

- 
  Add an environment for *anndata* v0.12.3. This is now the
  default.

- 
  Add a message recommending the anndataR package

Major changes

- 
  Update Python environments for compatibility with the new
  basilisk system

Minor changes

- 
  Adapt to changes in DelayedArray v0.35.2


NEWS from existing Data Experiment Packages
===================================


[ChIPDBData](/packages/ChIPDBData)
----------

- first version of the package


[curatedMetagenomicData](/packages/curatedMetagenomicData)
----------------------

                       Changes in version 3.16.1                        

- Updates to rownames = "short" due to downstream changes in mia
  (@TuomasBorman, #309)

- Use convertToPhyloseq in vignette instead of
  makePhyloseqFromTreeSummarizedExperiment (@cmirzayi, #310)

- Use rep_named instead of deprecated list_along (@LiNk-NY)

[curatedTCGAData](/packages/curatedTCGAData)
---------------

                       Changes in version 1.32.0                        

Bug fixes and minor improvements

- Added chunk labels to vignette

- Improvements to GitHub Actions workflows

[DoReMiTra](/packages/DoReMiTra)
---------

- Initial Bioconductor submission


[iModMixData](/packages/iModMixData)
-----------

                       Changes in version 0.99.0                        

Initial Release

[mCSEAdata](/packages/mCSEAdata)
---------

                       Changes in version 1.29.1                        

- Added support for EPICv2 platform.

[MetaScope](/packages/MetaScope)
---------

                        Changes in version 3.23                         

Bug fixes

- fixed the metascope_id() function to export the location of the id
  file

- Made significant improvements in memory usage in Metascope_id

                        Changes in version 3.22                         

Bug Fixes

- Fixed filter compression to take a pre-formatted character string in
  csv format

- Fixed minor issues in convert_animalcules() when non-linked taxonomy
  IDs are present

- Updated convert_animalcules to correct a bug in the final
  consolidation of reads, where more read counts were being added than
  in the original sample

[nmrdata](/packages/nmrdata)
-------

- Initial submission to Bioconductor.

- Provides curated 1D ^1H NMR spectra of rat (murine) urine samples
  from a bariatric surgery study.

- Includes two datasets:

- Processed dataset (bariatric) accessible via loadBariatric().

- Raw Bruker NMR experiment folders accessible via getRawExpDir().

[RforProteomics](/packages/RforProteomics)
--------------

                       Changes in version 1.47.2                        

- Use latest MSnbase

                       Changes in version 1.47.1                        

- Remove DEP suggests.

- Update vignettes: remove some references to old/outdated packages,
  add note about RforMS.


[SingleCellMultiModal](/packages/SingleCellMultiModal)
--------------------

                       Changes in version 1.22.0                        

Bug fixes and minor improvements

- Updated vignette authorship information

- Added sticker to README.md

- Enhancements to GitHub Actions workflows

[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.21.6                        

NEW FEATURES

- registration_pseudobulk() will now create the rowData()$gene_search
  if rowData()$gene_id and rowData()$gene_name are present. This makes
  the output be more in sync with requirements for run_app(). See
  https://github.com/LieberInstitute/spatialLIBD/pull/115 for details.
  Implemented by @lahuuki.

                       Changes in version 1.21.4                        

BUG FIXES

- @manishabarse fixed a bug in registration_wrapper() which led to
  NULL results for anova. See
  https://github.com/LieberInstitute/spatialLIBD/pull/110 for details.

                       Changes in version 1.21.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- @lahuuki edited registration_pseudobulk() to make it more compatible
  with scran::pseudobulkDGE().See
  https://github.com/LieberInstitute/spatialLIBD/pull/108 for details.

                       Changes in version 1.21.2                        

NEW FEATURES

- @lahuuki added the vis_image() for plotting just the histology
  image. See https://github.com/LieberInstitute/spatialLIBD/pull/107
  for details.

                       Changes in version 1.21.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- The documentation of registration_pseudobulk() has been expanded to
  further explain the logcounts() assay. That is, to highlight that
  this assay contains log2 CPM values computed with edgeR::cpm() and
  not log2 library-size normalized counts (as computed with
  scuttle::logNormCounts()). See
  https://support.bioconductor.org/p/9161754 and
  https://github.com/LieberInstitute/spatialLIBD/issues/106 by
  @kinnaryshah for more details.

[STexampleData](/packages/STexampleData)
-------------

                 Changes in version 1.17.1 (2025-10-05)                 

- add sticker

[TENET.ExperimentHub](/packages/TENET.ExperimentHub)
-------------------

                        Changes in version 1.1.2                        

- Correct Location_Prefix as per Bioconductor staff email

                        Changes in version 1.1.1                        

- Regenerate the example MultiAssayExperiment object with ggplot2 4.0,
  because it cannot load objects saved by previous ggplot2 versions.
  Since the reverse is also true, add a dependency on ggplot2 4.0.

- Correct the name of a dataset loaded from TENET in make-data.R,
  which was changed late in the development process but not changed in
  this script

- Update the link to the example TAD dataset to use the Internet
  Archive since the original URL is no longer accessible

- Do not suggest installing the development version of Bioconductor to
  use the stable version of the package, since it is no longer
  necessary

- Update installation instructions for Ubuntu 24.04 and add a missing
  required package

- Add CITATION file

- Add TENET to Suggests field in DESCRIPTION

                        Changes in version 1.1.0                        

- Release generated automatically by Bioconductor - no changes


Deprecated and Defunct Packages
===============================

**SOFTWARE:**

Thirty nine software packages were removed from this release (after being deprecated
in Bioc 3.21):

- AneuFinder, BEARscc, CBEA, chromstaR, coMET, crossmeta, dce, DeProViR, DIAlignR, Director, erma, GeneGeneInteR, genoCN, gespeR, girafe, GraphPAC, HTSeqGenie, iPAC, MAGeCKFlute, netDx, NeuCA, PanViz, pareg, paxtoolsr, PICS, PING, QuartPAC, RBioinf, ReactomeContentService4R, RGMQL, Rtreemix, SpacePAC, staRank, STdeconvolve, supraHex, synapter, trigger, TypeInfo, zlibbioc

Forty software packages are deprecated in this release and will be removed in Bioc 3.23:

- bgx, BiGGR, biodbHmdb, biodbNcbi, biodbNci, biodbUniprot, ccmap, CellScore, CINdex, cisPath, DEP, gpuMagic, Harshlight, hiAnnotator, hiReadsProcessor, hypeR, interactiveDisplay, interactiveDisplayBase, lapmix, LinTInd, lute, MADSEQ, netZooR, oppti, PhenStat, qckitfastq, ReactomeGraph4R, Repitools, Rfastp, rGADEM, rRDP, SARC, seqArchR, seqArchRplus, seqTools, Streamer, TitanCNA, TransView, traviz, XNAString


**EXPERIMENT DATA:** 

Three experimental data packages were removed from this release (after being
deprecated in BioC 3.21):

- benchmarkfdrData2019, parathyroidSE, synapterdata

Five experimental data packages are deprecated in this release and will be
removed in Bioc 3.23:

- AneuFinderData, chromstaRData, curatedCRCData, RGMQLlib, rRDPData

**ANNOTATION DATA:** 

Three annotation packages were removed from this release (after being deprecated
in Bioc 3.22).

- mirbase.db, targetscan.Hs.eg.db, targetscan.Mm.eg.db

No annotation packages are deprecated in this release.


**WORKFLOWS:** 

One workflow packages were removed from this release (after being deprecated in
Bioc 3.22).

- BiocMetaWorkflow

One workflow package is deprecated in this release and will be removed in Bioc 3.23:

- spicyWorkflow

**BOOKS:**

No books were removed from this release.

No books were deprecated in this release.
