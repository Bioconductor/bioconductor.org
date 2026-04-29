April 29, 2026

**Bioconductor:**

We are pleased to announce Bioconductor 3.23, consisting of
2418 software packages, 437 experiment data packages, 928 annotation
packages, 28 workflows and 8 books.

There are 94 new software packages, 5 new data experiment packages,
no new annotation packages, no new workflows, 2 new book, and many updates and
improvements to existing packages.

Bioconductor 3.23 is compatible with R 4.6, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher, macOS arm64 and Linux
arm64. This release will also include updated Bioconductor [Docker containers][2].

Note: Currently Bioconductor does not have active daily windows builders. We
will be testing the windows and mac binaries generated through the r-universe
system. These updated windows and mac binaries will be available shortly after this
release. We appreciate your patience as we make them available.

Thank you to everyone for your contribution to Bioconductor.

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.23](#getting-started-with-bioconductor-323)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.23
======================================

To update to or install Bioconductor 3.23:

1. Install R 4.6. Bioconductor 3.23 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 94 new software packages in this release of Bioconductor.

- [Aerith](/packages/Aerith) Visualisation of peptide isotopic peaks
  and SIP peptide spectra match (PSM). Filtration of high quality
  PSM. Accurate isotopic abundance calculation of peptide and
  metabolites. Visualisation of SIP proteomics results.

- [annoLinker](/packages/annoLinker) Fast annotation of genomic peaks
  using DNA interaction data by constructing interaction networks
  with igraph, where peaks overlapping any node in a connected
  subgraph are annotated with all genes in that subgraph. The
  annotation evidence could be visualized as either a network graph
  or a genomic track integrated with gene annotation information.

- [asuri](/packages/asuri) The ASURI (Analysis of SUrvival and
  patients RIsk prediction based on gene signatures) package
  discovers marker genes that are related to risk prediction
  capabilities and to a clinical variable of interest. It uses two
  main steps, including subsampling glmnet and unicox. The package
  implements robust functions to discover survival markers related to
  a clinical phenotype and to predict a risk score, allowing to study
  the patient's risk based on the gene signatures. Several plots are
  provided to visualise the relevance of the genes, the risk score,
  and patient stratification, as well as a robust version of the
  Kaplan-Meier curves.

- [atacInferCnv](/packages/atacInferCnv) The package prepares input
  scATAC-seq data and adapts for copy number variance profiling with
  InferCNV package usage. It has also various paramters to control
  the analysis (e.g. external normal reference usage, meta-cells, bin
  size, etc) and custom plot visualizations.

- [BatChef](/packages/BatChef) This package implements a variety of
  methods for batch correction in single-cell RNA sequencing
  (scRNA-seq) data. It incorporates quantitative metrics (e.g.
  Wasserstein distance, Adjusted Rand Index) to evaluate their
  performance. Furthermore, the package assists users in identifying
  and applying the optimal method for specific datasets.

- [Battlefield](/packages/Battlefield) Battlefield is a Swiss-army
  toolkit originally developed to define and extract spatial spots
  from specific tissue regions—such as front regions, niche borders,
  invasive margins, and cluster interfaces—using spatial
  transcriptomics data or clustered tissue maps. It has since been
  extended to support trajectory selection and layer inspection, and
  now provides a collection of low-level utilities for spatial
  transcriptomics analysis. These utilities are primarily intended to
  be reused within higher-level analytical packages. It is designed
  to work with sequencing-based platforms such as Visium at several
  resolutions and Visium HD(binned).

- [betterChromVAR](/packages/betterChromVAR) A much faster analytical
  implementation of chromVAR, with additional features, used to infer
  TF activity from (bulk or single-cell) ATAC-seq data and motif
  annotations (or binding probabilities). The package also includes
  the CVnorm normalization method based on the chromVAR logic.

- [BiocAzul](/packages/BiocAzul) Represents the OpenAPI v2 Azul API
  as an R object for performing requests. The infrastructure uses the
  AnVIL and rapiclient packages. Users can connect to either the
  AnVIL or Human Cell Atlas Data Explorers.

- [BiocBuildReporter](/packages/BiocBuildReporter) This package reads
  remote parquet files that have processed Bioconductor build report
  logs. Users may query the tables directly for specific information
  or use pre-defined helper functions for common queries. The logs
  processed are from https://bioconductor.org/checkResults/. In the
  future we will extend this package out to include processing of
  r-universe logs.

- [BiocMaintainerApp](/packages/BiocMaintainerApp) This package
  allows interactive viewing of package maintainer information. The
  Bioconductor Package Maintainer Application sends yearly
  verification emails to accept Bioconductor policies; this
  application also depicts maintainer status on opting in and if the
  email is deemed valid.

- [BiocPkgDash](/packages/BiocPkgDash) This package provides an
  interactive Shiny dashboard for Bioconductor package maintainers.
  It visualizes various package statuses, metadata, and development
  metrics, offering insights into package health and activity. This
  tool aims to support maintainers of multiple packages by filtering
  packages via maintainer email.

- [carnation](/packages/carnation) Highly interactive & modular shiny
  app to explore three facets of RNA-Seq analysis: differential
  expression (DE), functional enrichment and pattern analysis.
  Several visualizations are implemented to provide a wide-ranging
  view of data sets. For DE analysis, we provide PCA plot, MA plot,
  Upset plot & heatmaps, in addition to a highly customizable gene
  plot.  Seven different visualizations are available for functional
  enrichment analysis, and we also support gene pattern analysis.
  Genes of interest can be tracked across all modules using the gene
  scratchpad. In addition, carnation provides an integrated platform
  to manage multiple projects and user access that can be run on a
  central server to share with collaborators.

- [CellMentor](/packages/CellMentor) Implements supervised cell
  type-aware non-negative matrix factorization (NMF) for dimensional
  reduction in single-cell RNA sequencing analysis. The package
  provides methods for incorporating cell type information into the
  dimensionality reduction process, enabling improved visualization
  and downstream analysis of single-cell data while preserving
  biological structure. CellMentor employs a unique loss function
  that simultaneously minimizes variation within known cell
  populations while maximizing distinctions between different cell
  types, enabling effective transfer of learned patterns from labeled
  reference datasets to new unlabeled data.

- [ClonalSim](/packages/ClonalSim) ClonalSim generates realistic
  mutational profiles of tumor samples with hierarchical clonal
  structure. It simulates founder, shared, and private mutations with
  biologically realistic noise models including intra-tumor
  heterogeneity (Beta distribution) and technical sequencing noise
  (negative binomial depth variation, binomial read sampling, base
  errors). The package is designed for benchmarking variant callers,
  testing clonal deconvolution algorithms, and teaching tumor
  heterogeneity concepts.

- [ClusterGVis](/packages/ClusterGVis) Provides a streamlined
  workflow for clustering and visualizing gene expression patterns,
  particularly from time-series RNA-Seq and single-cell experiments.
  The package is designed to integrate seamlessly within the
  Bioconductor ecosystem by operating directly on standard data
  classes such as `SummarizedExperiment` and `SingleCellExperiment`.
  It implements common clustering algorithms (e.g., k-means, fuzzy
  c-means) and generates a suite of publication-ready visualizations
  to explore co-expressed gene modules. Functions are also included
  to facilitate the visualization of clustering results derived from
  other popular tools.

- [CompensAID](/packages/CompensAID) The CompensAID is an automated
  quality control tool, which determines for each marker combination
  in the FCS file, whether there a potential presence of reference
  errors. Such reference errors, which represent themselves in the
  form of skewed populations, are detected by integrating the
  Secondary Stain Index (SSI) score. Marker combinations with an SSI
  < 1 are flagged by CompensAID.

- [CrcBiomeScreen](/packages/CrcBiomeScreen) A developed and
  benchmarked reproducible machine learning framework for
  microbiome-based colorectal cancer (CRC) screening. By
  systematically evaluating normalization strategies, taxonomic
  resolutions, and class imbalance handling. This R package allows
  users to apply the full pipeline or selectively run specific
  components depending on their analytical needs. It establishes a
  scalable foundation for developing interpretable microbiome-based
  screening tools to support early CRC detection. This approach could
  be easily implemented in a national screening programme, to improve
  early detection rates for this disease.

- [damidBind](/packages/damidBind) The damidBind package provides a
  straightforward formal analysis pipeline to analyse and explore
  differential DamID binding, gene transcription or chromatin
  accessibility between two conditions. The package imports processed
  data from DamID-seq experiments, either as external raw files in
  the form of binding bedGraphs and GFF/BED peak calls, or as
  internal lists of GRanges objects. After optionally normalising
  data, combining peaks across replicates and determining
  per-replicate peak occupancy, the package links bound loci to
  nearby genes.  For RNA Polymerase DamID data, the package
  calculates occupancy over genes, and optionally calcualates the FDR
  of significantly-enriched gene occupancy. damidBind then uses
  either limma (for conventional log2 ratio DamID binding data) or
  NOIseq (for counts-based CATaDa chromatin accessibility data) to
  identify differentially-enriched regions, or differentially
  epxressed genes, between two conditions. The package provides a
  number of visualisation tools (volcano plots, Gene Ontology
  enrichment plots via ClusterProfiler and proportional Venn diagrams
  via BioVenn for downstream data exploration and analysis. An
  powerful, interactive IGV genome browser interface (powered by
  Shiny and igvShiny) allows users to rapidly and intuitively assess
  significant differentially-bound regions in their genomic context.

- [decemedip](/packages/decemedip) The R package decemedip is a novel
  computational paradigm developed for inferring the relative
  abundances of cell types and tissues measure by methylated DNA
  immunoprecipitation sequencing (MeDIP-Seq). This paradigm allows
  using reference data from other technologies such as microarray or
  WGBS.

- [DenoIST](/packages/DenoIST) DenoIST identifies and removes
  contamination in Image-based Spatial Transcriptomics data, using a
  transposed poisson mixture model with local neighbourhood offsets
  to infer genes that are likely to be due to neighbourhood
  contamination rather than endogenous expression.

- [dominatR](/packages/dominatR) dominatR is an R package for
  quantifying and visualizing feature dominance in datasets.
  dominatR applies concepts drawn from physics such as center of mass
  and shannon's entropy to effectively visualize features (e.g.
  genes) that are present within a specific context or condition. The
  package integrates, dataframes, matrices and SummerizedExperiment
  objects and is able to perform common genomic normalization
  methods. The key aspect is the generation of plots that serve to
  highlight context-relevant feature dominance.

- [DOTSeq](/packages/DOTSeq) Differential open reading frame (ORF)
  translation analysis framework for ribosome profiling (Ribo-seq)
  with matched RNA-seq. Implements (i) Differential ORF Usage (DOU),
  a beta-binomial generalized linear model that models the expected
  proportion of Ribo-seq versus RNA-seq reads mapping to each ORF
  within a gene, and (ii) ORF-level Differential Translation
  Efficiency (DTE), a negative binomial GLM that capture changes in
  translation efficiency of individual ORFs across experimental
  conditions. Supports ORF-level read summarization for bulk and
  single-cell Ribo-seq.

- [drugfindR](/packages/drugfindR) This package provides a convenient
  way to access the LINCS Signatures available in the iLINCS
  database. These signatures include Consensus Gene Knockdown
  Signatures, Gene Overexpression signatures and Chemical Perturbagen
  Signatures. It also provides a way to enter your own transcriptomic
  signatures and identify concordant and discordant signatures in the
  LINCS database.

- [epiRomics](/packages/epiRomics) Integrates various levels of
  epigenomic information, including ChIP-seq, histone modification,
  ATAC-seq, and RNA-seq data. Regulatory network analysis uses
  combinatory approaches to infer regions of significance, such as
  enhancers. Downstream analysis identifies co-occurrence of
  epigenomic data at regions of interest. Visualization functions
  display multi-track genomic views with signal overlays. Please
  contact <ammawla@ucdavis.edu> for suggestions, feedback, or bug
  reporting.

- [epiSeeker](/packages/epiSeeker) This package implements functions
  to analyze multi-omics epigenetic data. Data of fragment type and
  base type are supported by epiSeeker. It provides functions to
  retrieve the nearest genes around the peak, annotate genomic region
  of the peak, statistical methods to estimate the significance of
  overlap among peak data sets, and motif analysis. It incorporates
  the GEO database for users to compare their own dataset with those
  deposited in the database. The comparison can be used to infer
  cooperative regulation and thus can be used to generate hypotheses.
  Several visualization functions are implemented to summarize the
  coverage of the peak experiment, average profile and heatmap of
  peaks binding to TSS regions, genomic annotation, distance to TSS,
  overlap of peaks or genes, and the single-base resolution
  epigenetic data by considering the strand, motif, and additional
  information.

- [ExpoRiskR](/packages/ExpoRiskR) ExpoRiskR provides tools for
  exposure-aware multi-omics risk modeling in translational and
  environmental health studies. The package aligns sample identifiers
  across exposure and multi-omics blocks, performs lightweight
  preprocessing, and fits exposure-adjusted association models to
  build interpretable microbe–metabolite networks. It also computes
  simple exposure perturbation summaries and generates
  publication-ready visualizations. Workflows support both
  matrix-based inputs and SummarizedExperiment objects.

- [fastRanges](/packages/fastRanges) High-performance interval
  overlap and join operations for 'IRanges' and 'GenomicRanges'. The
  package provides deterministic multithreaded overlap computation,
  reusable subject indexes for repeated queries, and join helpers
  that keep range metadata in a consistent output grammar.

- [fourSynergy](/packages/fourSynergy) fourSynergy is an ensemble
  algorithm leveraging synergies among the existing 4C-seq algorithms
  r3C-seq, peakC, r.4cker and fourSig. It uses a weighted voting
  approach to perform improved interaction calling. fourSynergy
  supports also differential interaction calling.

- [fRagmentomics](/packages/fRagmentomics) A user-friendly R package
  that enables the characterization of each cfDNA fragment
  overlapping one or multiple mutations of interest, starting from a
  sequencing file containing aligned reads (BAM file). fRagmentomics
  supports multiple mutation input formats (e.g., VCF, TSV, or string
  "chr:pos:ref:alt" representation), accommodates one-based and
  zero-based genomic conventions, handles mutation representation
  ambiguities, and accepts any reference file and species in FASTA
  format. For each cfDNA fragment, fRagmentomics outputs its size,
  its 3' and 5' sequences, and its mutational status. Optionally,
  when users set apply_bcftools_norm = TRUE, fRagmentomics invokes
  the external command-line tool bcftools norm to left-align and
  normalize variants. If bcftools is not found on the system PATH
  while this option is enabled, the function errors. The package does
  not install external software; see the INSTALL file for per-OS
  instructions.

- [fraq](/packages/fraq) High-throughput extensible toolkit for
  processing FASTQ data. The goal of this package is to empower users
  to quickly build out small programmatic 'kernels' to define any
  FASTQ processing task they may need. Builds on Intel TBB’s flow
  graph to orchestrate concurrent I/O and data processing; throughput
  can be as fast as compression and disk speed allows. The package
  also ships with a suite of predefined kernels for common FASTQ
  tasks.

- [GenomicCoordinates](/packages/GenomicCoordinates) Extends string
  parsing capabilities for genomic coordinates, supporting various
  formats including comma-separated numbers, space-delimited
  coordinates, and automatic detection of GRanges, GPos, and
  GInteractions objects.

- [glycoTraitR](/packages/glycoTraitR) GlycoTraitR is an R package
  for analyzing glycoproteomics data, particularly
  glycopeptide-spectrum matches (GPSMs). It supports results
  generated by the pGlyco3 and Glyco-Decipher search engines. The
  package parses glycan structures, computes monosaccharide
  compositions and structural traits, and performs differential
  analysis of glycan heterogeneity. It constructs trait-by-PSM
  matrices stored in a SummarizedExperiment object, supports
  user-defined structural motifs, and provides visualization
  utilities for interpreting glycan trait changes.

- [GOaGO](/packages/GOaGO) GO-a-GO annotates Gene Ontology terms that
  are enriched in a given set of gene pairs. The enrichment is
  calculated from a permutation test for overrepresentation of gene
  pairs that are associated with a shared term. Such gene pairs are
  counted for the original set of gene pairs and compared against
  randomized sets in which the structure of the pairs is preserved,
  but the gene identities (including the associated terms) are
  permuted.

- [GOfan](/packages/GOfan) GOfan provides an intuitive and compact
  visualization of Gene Ontology (GO) enrichment results using a
  sunburst layout inspired by SynGO, preserving hierarchical
  relationships among GO terms and allowing color-based encoding of
  information such as p-values or gene counts. By converting complex
  GO DAGs into clean, circular representations, it allows researchers
  to quickly grasp the hierarchical structure and biological
  significance of enriched terms. The interactive and customizable
  visualizations facilitate exploration of key GO categories,
  enhancing interpretation and presentation of enrichment analyses.

- [GraphExperiment](/packages/GraphExperiment) GraphExperiment
  provides users and developers with an S4 class that extends
  `SingleCellExperiment` by offering infrastructure to store and
  retrieve networks (`igraph` objects) representing how features are
  associated with each other. The class was designed to store
  networks inferred from high-dimensional quantitative data,
  including gene coexpression networks (GCNs), gene regulatory
  networks (GRNs), and co-abundance networks (from proteomics and
  metabolomics), as well as networks inferred from other types of
  data (e.g., protein-protein interactions).

- [GSABenchmark](/packages/GSABenchmark) GSABenchmark is a package
  designed for benchmarking scRNA-seq gene set analysis (scGSA)
  methods. It provides both traditional and novel benchmark metrics,
  as well as visualization tools. Currently, GSABenchmark supports 17
  scGSA methods.

- [hammers](/packages/hammers) hammers is a utilities suite for
  scRNA-seq data analysis compatible with both Seurat and
  SingleCellExperiment. It provides simple tools to address tasks
  such as retrieving aggregate gene statistics, finding and removing
  rare genes, performing representation analysis, computing the
  center of mass for the expression of a gene of interest in
  low-dimensional space, and calculating silhouette and
  cluster-normalized silhouette.

- [HiSpaR](/packages/HiSpaR) Provides R bindings for HiSpa, a
  hierarchical Bayesian model for inferring three-dimensional
  chromatin structures from Hi-C contact matrices using Markov Chain
  Monte Carlo (MCMC) sampling. The package implements a cluster-based
  hierarchical approach that efficiently handles large-scale Hi-C
  datasets. It uses Rcpp and RcppArmadillo for efficient C++
  integration with the original HiSpa C++ implementation, enabling
  fast computation of chromatin structure inference through parallel
  MCMC sampling.

- [HistoImagePlot](/packages/HistoImagePlot) Create side-by-side
  visualizations of tissue thumbnail image and HoverNet cell
  segmentation with colored cell type labels. Functionality
  automatically retrieves the thumbnail image associated with a
  HoverNet JSON file and overlays the segmentation data. This package
  is intended for researchers working with histopathological images,
  facilitating exploratory analysis, and integrates with the
  imageFeatureTCGA Bioconductor package.

- [ImageArray](/packages/ImageArray) ImageArray provides a framework
  for on-disk and in-memory image arrays, specifically for pyramidal
  images stored in HDF5, Zarr and life sciences image file formats
  (OME Bio-Formats).

- [imageFeatureTCGA](/packages/imageFeatureTCGA) The package imports
  data from HoverNet, and ProvGigaPath pipelines. Pipeline output
  data are hosted in a self-owned online repository. Package
  functionality conveniently incorporates pipeline data into existing
  MultiAssayExperiment instances from curatedTCGAData.

- [imageTCGAutils](/packages/imageTCGAutils) Utility functions for
  working with CONCH data, listing remote files. One function assigns
  HoverNet nuclei to ProvGigaPath tiles with a scale factor to align
  coordinates. Provides internal utility functions for
  'imageFeatureTCGA' and most functions are not meant for end users.

- [immLynx](/packages/immLynx) A comprehensive toolkit that bridges
  popular Python-based immune repertoire analysis tools and Hugging
  Face protein language models into the R environment. Provides
  unified interfaces for TCR distance calculations (tcrdist3),
  sequence generation probability (OLGA), selection inference
  (soNNia), clustering (clusTCR), protein embeddings (ESM-2),
  metaclone discovery (metaclonotypist). Fully compatible with the
  scRepertoire and immApex ecosystem for single-cell immune
  repertoire analysis.

- [immReferent](/packages/immReferent) Provides a consistent
  interface for downloading, storing, and accessing immune receptor
  (TCR/BCR) and HLA sequences from IMGT, IPD-IMGT/HLA, and OGRDB
  (AIRR-C). Supports export to popular analysis tools including
  MiXCR, TRUST4, Cell Ranger, and IgBLAST. This package serves as a
  core dependency for immunogenomics packages, ensuring reliable and
  high-quality sequence access with local caching for
  reproducibility.

- [jvecfor](/packages/jvecfor) Drop-in replacement for
  BiocNeighbors::findKNN using the jvecfor Java library, which builds
  on the jvector library to leverage the Java Vector API for portable
  SIMD acceleration across AVX2, AVX-512, and ARM NEON hardware.
  jvecfor/jvector implements HNSW-DiskANN approximate search and
  VP-tree exact search. The package achieves approximately 2x speedup
  over Annoy-based search at n >= 50K cells while returning output
  structurally identical to BiocNeighbors, making it suitable for
  seamless integration into existing Bioconductor single-cell
  workflows. Convenience wrappers delegate shared nearest-neighbor
  (SNN) and k-nearest-neighbor (KNN) graph construction to the
  bluster package.

- [LACHESIS](/packages/LACHESIS) This package provides modalities to
  analyze tumor evolution from whole genome sequencing data. In
  particular, it provides estimates of mutation densities at genomic
  segments and uses these to time the origin of the tumor.

- [lcmsPlot](/packages/lcmsPlot) lcmsPlot is an R package designed
  for visualising Liquid Chromatography-Mass Spectrometry (LC-MS)
  data with publication-ready high-quality plots. The package enables
  users to generate and customise chromatograms, mass traces,
  spectra, and more with fine-tuned aesthetics and annotation
  options.

- [leapR](/packages/leapR) leapR is a package that identifies
  pathways that are enriched across diverse 'omics experiments. It
  leverages any tabular expression data (proteomics, transcriptomics)
  using the `SummarizedExperiment` object. It works with any pathway
  in the .gct file format.

- [lncRna](/packages/lncRna) Provides a complete workflow for the
  identification, analysis, and functional annotation of long
  non-coding RNAs (lncRNAs) from RNA-Seq data. The package includes
  functions for filtering transcripts from GTF files, evaluating the
  performance of multiple coding potential prediction tools (e.g.,
  CPC2, PLEK, CPAT), and summarizing their agreement. It enables
  systematic performance analysis of individual tools, "at least N"
  tool consensus, and all possible tool combinations. Functional
  analysis is supported through the identification of potential cis-
  and trans-acting interactions with protein-coding genes, followed
  by enrichment analysis. Results can be visualized using a variety
  of plots, including radar plots, clock plots, and interactive
  Sankey diagrams.

- [LRDE](/packages/LRDE) Provides hurdle negative binomial models for
  differential expression analysis with long-read RNA-Seq data.

- [MDSvis](/packages/MDSvis) This package implements visulization of
  Multi Dimensional Scaling (MDS) results.

- [MeLSI](/packages/MeLSI) MeLSI (Metric Learning for Statistical
  Inference) is a novel machine learning method for microbiome data
  analysis that learns optimal distance metrics to improve
  statistical power in detecting group differences. Unlike
  traditional distance metrics (Bray-Curtis, Euclidean, Jaccard),
  MeLSI adapts to the specific characteristics of your dataset to
  maximize separation between groups. The method uses an ensemble of
  weak learners to identify which microbial features drive group
  differences, providing both improved statistical power and
  biological interpretability through feature importance weights.

- [MetaboAnnotatoR](/packages/MetaboAnnotatoR) Performs feature
  annotations on LC-MS All-ion fragmentation datasets using fragment
  ion libraries.

- [metabom8](/packages/metabom8) Tools for 1D NMR metabolomics
  workflows, including import and preprocessing of Bruker
  experiments, multivariate modeling (PCA, PLS, OPLS) and model
  analytics and validation (y-permutations, cv-anova).
  Performance-critical routines are implemented in C++ and use the
  Armadillo and Eigen linear algebra libraries to improve runtime.

- [MetaProViz](/packages/MetaProViz) MetaProViz can analyse standard
  metabolomics and exometabolomics data (CoRe). It performs
  pre-processing including feature filtering, missing value
  imputation, normalisation and outlier detection. It performs
  functional analysis including differential metabolite analysis
  (DMA), clustering based on regulatory rules (MCA) and contains
  different visualisation methods to extract biological interpretable
  graphs and saves them in a publication ready format.

- [MutSeqR](/packages/MutSeqR) Standard methods for analysis of
  mutation data following error- corrected sequencing (ECS) for the
  purpose of mutagencity assessment. Functions include importing the
  mutation lists provided by a variant caller, and a set of
  analytical tools for statistical testing and visualization of
  mutation data; comparison to COSMIC and/or germline signatures;
  etc.

- [OAtools](/packages/OAtools) Provides a suite of R functions to
  analyze gene expression experiments on the OpenArray real-time PCR
  platform. OAtools fits logistic regressions to fluorescence curves
  to distinguish between real amplification and false positives.
  OAtools supports data import, analysis, and visualization through
  plots and a dynamic HTML report.

- [parati](/packages/parati) Infers maternal and paternal transmitted
  and non-transmitted alleles from phased trio genotype data. The
  package supports SNP-level analyses of genetic nurture and
  transgenerational effects. It interoperates with Bioconductor VCF
  infrastructure through support for VariantAnnotation::VCF objects
  and returns R objects for downstream analysis.

- [plaid](/packages/plaid) PLAID (Pathway Level Average Intensity
  Detection) is an ultra-fast method to compute single-sample
  enrichment scores for gene expression or proteomics data. For each
  sample, plaid computes the gene set score as the average intensity
  of the genes/proteins in the gene set. The output is a gene set
  score matrix suitable for further analyses.

- [PlinkMatrix](/packages/PlinkMatrix) This package provides a
  DelayedArray interface for plink bed files. There is support for
  interfacing to plink genotype data via RangedSummarizedExperiment.
  Example data from the GEUVADIS project (internationalgenome.org)
  are used for demonstration.

- [posDemux](/packages/posDemux) Demultiplexing and filtering
  utilities intended for reads with combinatorial barcodes (i.e.
  PETRI-seq and SPLiT-seq). The demultiplexer algorithm uses the
  position of the segments to extract and compare the barcodes with
  the reference (whitelist). A Shiny application is provided to
  interactively select cutoffs for which barcode combinations to
  keep.

- [postNet](/packages/postNet) A tool that enables in silico
  identification, integration, and modeling of mRNA features that
  influence post-transcriptional regulation of gene expression at a
  transcriptome-wide scale.

- [proBatch](/packages/proBatch) These tools facilitate batch effects
  analysis and correction in high-throughput experiments. It was
  developed primarily for mass-spectrometry proteomics (DIA/SWATH),
  but could also be applicable to most omic data with minor
  adaptations. The package contains functions for diagnostics
  (proteome/genome-wide and feature-level), correction (normalization
  and batch effects correction) and quality control. Non-linear
  fitting based approaches were also included to deal with complex,
  mass spectrometry-specific signal drifts.

- [PTMods](/packages/PTMods) An interface to the community supported
  database for amino acid/protein modifications using mass
  spectrometry.

- [queeems](/packages/queeems) Biological inferences obtained from
  molecular data are only as good as the extent of evolutionary
  signatures retained in the genetic data. Techniques available to
  quantify these signatures are largely targeted towards phylogeny
  reconstruction and they often rely on adhoc hypothesis tests of
  significance. I present a Bayesian function that assesses whether a
  set of genetic sequences are saturated. That is, it is useful for
  determining whether the evolutionary information in the sequences
  has eroded with time. Site specific Bayes factors are generated
  with respect to codon bases to allow for straightforward
  applications in extensive computational biology inquiries,
  including natural selection analyses.

- [RankMap](/packages/RankMap) RankMap is a fast and scalable tool
  for reference-based cell type annotation of single-cell and spatial
  transcriptomics data. It uses ranked gene expression and
  multinomial regression to achieve robust predictions, even with
  partial gene coverage. Compatible with Seurat,
  SingleCellExperiment, and SpatialExperiment objects, RankMap offers
  flexible preprocessing and significantly faster runtime than tools
  like SingleR, Azimuth, and RCTD.

- [RBedMethyl](/packages/RBedMethyl) Bioconductor-native
  infrastructure for handling large nanoporetech modkit bedMethyl
  pileup files from ONT data using HDF5Array and DelayedArray.

- [Rega](/packages/Rega) The European Genome-phenome Archive (EGA)
  provides long-term storage and controlled sharing of personally
  identifiable genetic data. The Rega package offers a streamlined
  and extensible R interface to the EGA API, facilitating the
  programmatic upload of metadata. GEO-like Excel submission template
  is provided as a default method of organizing submission metadata.

- [RFGeneRank](/packages/RFGeneRank) Tools to harmonize bulk RNA-seq
  matrices, optionally apply batch correction, and train
  cross-validated classification models using ranger, glmnet, or
  xgboost. Supports leakage-safe feature selection, permutation
  importance, SHAP-based interpretability, and calibration methods
  (Platt or isotonic). Provides stability metrics across folds,
  embeddings (PCA/UMAP), ROC visualization, SHAP dependence plots,
  and tidy ranked-gene tables for downstream analysis.

- [RNAshapeQC](/packages/RNAshapeQC) RNAshapeQC provides
  coverage-shape-based quality control (QC) metrics for mRNA-seq and
  total RNA-seq data. It supports per-gene pileup construction from
  BAM files as well as toy datasets for quick-start examples. The
  package implements protocol-specific metrics, including decay rate
  (DR), degradation score (DS), mean coverage depth (MCD), window
  coefficient of variation (wCV), area under the curve (AUC), and
  shape-based sample-level indices. RNAshapeQC also includes
  HPC-friendly functions for per-gene batch processing and
  cross-study pileup generation. This package enables interpretable,
  protocol-specific QC assessments for diverse RNA-seq workflows.

- [scConform](/packages/scConform) Builds prediction interval for
  cell type annotation using conformal inference and conformal risk
  control. It provides two main methods. The first one gives
  prediction intervals with coverage guarantees based on standard
  conformal inference. The second one instead gives hierarchical
  prediction intervals that are consistent with the cell ontology.

- [scECODA](/packages/scECODA) The scECODA R package provides a
  complete workflow for the analysis and visualization of
  compositional data, primarily focusing on cell type proportions
  derived from single-cell data. It implements specialized methods,
  such as the Centered Log-Ratio (CLR) transformation, to properly
  analyze proportional data while avoiding the biases introduced by
  the compositional constraint. The package encapsulates data
  management, transformation, and analysis into a single
  SummarizedExperiment object, offering downstream tools for
  dimensionality reduction via PCA, calculating critical metrics like
  the Adjusted Rand Index (ARI) and Modularity to quantify sample
  grouping quality, and generating high-quality visualizations like
  heatmaps and scatter plots.

- [scLang](/packages/scLang) scLang is a suite for package
  development for scRNA-seq analysis. It offers functions that can
  operate on both Seurat and SingleCellExperiment objects. These
  functions are primarily aimed to help developers build tools
  compatible with both types of input.

- [scPassport](/packages/scPassport) Stamps Seurat,
  SingleCellExperiment, and SummarizedExperiment objects with a
  persistent metadata passport. For Seurat objects the passport is
  stored in the misc slot; for SingleCellExperiment and
  SummarizedExperiment objects it is stored in the metadata slot.
  Tracks animal info, experiment details, lineage (parent/child
  relationships), RDS registry numbers, processing logs, and custom
  fields. Includes an interactive Shiny gadget to fill and update the
  passport, and a read mode to print the full passport to console.
  The passport persists inside the RDS file with no external files
  needed.

- [scToppR](/packages/scToppR) scToppR provides an easy-to-use API
  wrapper for the ToppGene web platform, used for gene ontology and
  functional enrichment research. The package also integrates
  visualization tools, making it a convenient tool directly
  connecting ToppGene to code-based workflows in R. The tool can also
  easily save results into different formats.

- [scTypeEval](/packages/scTypeEval) scTypeEval provides tools to
  evaluate and validate cell type classifications in single-cell
  transcriptomics when ground truth labels are limited or
  unavailable. Results are organized in an S4 object that integrates
  processed data, dimensional reductions, dissimilarity assays, and
  consistency metrics computed across samples. The workflow includes
  preprocessing and feature selection, principal component analysis,
  computation of dissimilarity matrices, internal validation metrics
  (for example, silhouette-based summaries), and visualization
  utilities to inspect heatmaps and PCA plots. Functions support
  common single-cell containers and enable comparison of clustering
  and labeling strategies across datasets.

- [SEMPLR](/packages/SEMPLR) SEMPLR computes transcription factor
  binding affinity scores for genomic positions and genetic variants.
  Scores are computed from SNP Effect Matrices (SEMs) produced by
  SEMpl. 223 pre-computed SEMs are included with the package or
  custom sets can be provided. Enrichment can be tested among sets of
  genomic positions to determine if transcription factor binding
  events occur more often than expected. Comparing binding affinity
  scores between alleles can reveal differences in transcription
  factor binding with genetic variation. This package also includes
  several visualization functions to view scores both on the motif
  and variant/position level.

- [Seqtometry](/packages/Seqtometry) This package provides functions
  used in Seqtometry (Kousnetsov et al. 2024), a method for analyzing
  single cell (scRNA-seq or scATAC-seq) data via signature (gene set)
  enrichment scores. The Seqtometry scores may be useful for
  annotating or characterizing cells, either in a flow cytometry like
  workflow (where scores are standalone features used for progressive
  partitoning as described in the Seqtometry publication) or in a
  cluster-based workflow (as features of clusters). The exported
  impute function (a port of Python's MAGIC-impute, van Dijk et al.
  2018), may also be useful for single cell analysis on its own.

- [sfi](/packages/sfi) Data analysis for Single File Injections(SFIs)
  mode LC-MS analysis. In SFIs mode, pooled samples are initially
  injected to serve as reference peaks for subsequent analyses.
  Repeated injections of individual samples are then performed at
  fixed time intervals using isocratic elution. This package provides
  the functions to analyze data from SFIs mode including peak picking
  and peak reassignment.

- [singIST](/packages/singIST) Provides with toolkits to implement a
  full singIST analysis with pseudobulked Seurat objects of disease
  models and human data.

- [SMTrackR](/packages/SMTrackR) The package uses exogenous enzyme
  imprinted information to map protein-DNA binding on individual
  sequenced DNA molecules. For example, GpC methyltransferase, CpG
  methyltransferase, and Adenine methyltransferases. Public datasets
  from such assays are compiled into tracks, and hosted at public
  servers like Galaxy for their seamless access by this package.

- [SpatialArtifacts](/packages/SpatialArtifacts) SpatialArtifacts
  provides a data-driven two-step workflow to identify, classify, and
  handle spatial artifacts in spatial transcriptomics data. The
  package combines median absolute deviation (MAD)-based outlier
  detection with morphological image processing (fill, outline, and
  star patterns) to detect edge and interior artifacts. It supports
  multiple platforms including 10x Genomics Visium (standard and HD),
  allowing for consistent quality control across different spatial
  resolutions.

- [SpiecEasi](/packages/SpiecEasi) Estimate networks from the
  precision matrix of compositional microbial abundance data.

- [SpliceImpactR](/packages/SpliceImpactR) Works by taking in
  processed data from the HIT Index and/or rMATS and identifying how
  differentially used alternative RNA processing events lead to
  changes in protein function through various means. Primarily this
  is done through protein similarity, functional protein domain
  analysis, and domain-domain interaction changes. Notably, we both
  identify alterantive RNA processing event 'swaps' across condition
  and are able to perform holistic analyses regarding the impact of
  different RNA processing events.

- [splicelogic](/packages/splicelogic) Translate differential
  transcript usage results into discrete splice events.

- [SpNeigh](/packages/SpNeigh) SpNeigh provides methods for
  neighborhood-aware analysis of spatial transcriptomics data. It
  supports boundary detection, spatial weighting (centroid- and
  boundary-based), spatially informed differential expression using
  spline-based models, and spatial enrichment analysis via the
  Spatial Enrichment Index (SEI). Designed for compatibility with
  Seurat objects, SpatialExperiment objects and spatial data frames,
  SpNeigh enables interpretable, publication-ready analysis of
  spatial gene expression patterns.

- [staRgate](/packages/staRgate) An R-based automated gating pipeline
  for flow cytometry data designed to mimic the manual gating
  strategy of defining flow biomarker positive populations relative
  to a unimodal background population to include cells with varying
  intensities of marker expression. The pipeline’s main feature is a
  flexible density-based gating strategy capable of capturing varying
  scenarios based on marker expression patterns to analyze a
  29-marker flow panel that characterizes T-cell lineage,
  differentiation, and functional states.

- [StatescopeR](/packages/StatescopeR) StatescopeR is an R wrapper
  around Statescope, a computational framework designed to discover
  cell states from cell type-specific gene expression profiles
  inferred from bulk RNA profiles.

- [tidyexposomics](/packages/tidyexposomics) The tidyexposomics
  package is designed to facilitate the integration of exposure and
  omics data to identify exposure-omics associations. We structure
  our commands to fit into the tidyverse framework, where commands
  are designed to be simplified and intuitive. Here we provide
  functionality to perform quality control, sample and exposure
  association analysis, differential abundance analysis, multi-omics
  integration, and functional enrichment analysis.

- [tidyprint](/packages/tidyprint) Provides customized print methods
  for 'SummarizedExperiment' objects to enhance readability and
  usability within a tidy workflow. It offers consistent,
  tidyverse-aligned console displays, including alternative tibble
  abstractions for large genomic data to improve discoverability and
  interpretation. The package also includes unified, contextual
  messaging utilities intended for the 'tidyomics' ecosystem.

- [toppgene](/packages/toppgene) The ToppGene Suite is a one-stop
  portal for gene list enrichment analysis and candidate gene
  prioritization based on functional annotations and protein
  interactions network.  Although the ToppCluster web application
  provides convenient graphical access to the ToppGene Suite, the
  OpenAPI 3.0 compliant interface of ToppGene is better suited for
  automation and reproducibility.  This package includes Bioconductor
  class interfaces and biological examples.

- [VISTA](/packages/VISTA) The VISTA (Visualization and Integrated
  System for Transcriptomic Analysis) platform streamlines
  differential expression workflows by wrapping DESeq2 and edgeR into
  a SummarizedExperiment-based container with consistent metadata.
  The package includes visualization utilities, MSigDB enrichment
  helpers, and optional deconvolution support to simplify interactive
  exploration of RNA-seq experiments.

- [wavFeatExt](/packages/wavFeatExt) Provides tools for simulating
  copy-number alteration (CNA) profiles, applying a non-decimated
  Haar wavelet transform to genomic signals, and extracting
  wavelet-derived features for use in supervised learning. Multiple
  machine learning methods including lasso and elastic-net
  regularisation, random forest, partial least squares, neural
  networks and k-nearest neighbours are implemented to train
  predictive models from genomic feature vectors. The workflow
  enables end-to-end analysis from CNA simulation to feature
  extraction and classification.

- [ZarrArray](/packages/ZarrArray) The ZarrArray package leverages
  the Rarr package to bring Zarr datasets in R as DelayedArray
  objects. The main class in the package is the ZarrArray class. A
  ZarrArray object is an array-like object that represents a Zarr
  dataset in R. ZarrArray objects are DelayedArray derivatives and
  therefore support all operations (delayed or block-processed)
  supported by DelayedArray objects.


New Data Experiment Packages
=====================

There are 5 new data experiment packages in this release of Bioconductor.

- [DMRsegaldata](/packages/DMRsegaldata) Data package providing
  example DNA methylation files used in the DMRsegal vignette and
  examples. Includes a sorted beta matrix as a tab-delimited,
  bgzip-compressed file and a matching phenotype table. The data
  contains 10 healthy and 10 cancer samples, and preprocessing has
  already been performed on the beta values.

- [dominatRData](/packages/dominatRData) dominatRData is a data
  package useful for showcasing dominatR examples. dominatR is an R
  package for quantifying and visualizing feature dominance in
  datasets. dominatR makes use of entropy-based triangular
  projections and compositional comparison metrics.

- [EMTscoreData](/packages/EMTscoreData) Provides 12 single-cell
  RNA-seq datasets profiling epithelial–mesenchymal transition (EMT)
  in human cancer cell lines (MCF7, OVCA420, DU145, and A549) under
  TGF-beta stimulation, kinase inhibition, and time-course
  conditions, as reported by Cook DP and Vanderhyden BC (2020). The
  datasets are distributed via ExperimentHub as SingleCellExperiment
  objects.

- [HumanRetinaLRSData](/packages/HumanRetinaLRSData) Dataset package
  containing gene and isoform count matrices, and sample metadata for
  long-read direct cDNA sequencing of human retinal organoids, 2D
  retinal ganglion cell (RGC) cultures, and flowthrough fractions
  from H9 and EP1 iPSC cell lines. Data were generated using Oxford
  Nanopore Technology (ONT) direct cDNA sequencing and mapped to the
  GRCh38 reference genome (GENCODE v46 annotation). The package
  provides accessor functions returning SummarizedExperiment objects
  for gene-level counts, isoform-level counts, and a matrix of
  allele-specific expression (ASE) gene counts. Data files are stored
  in flat CSV format in an Open Science Framework (OSF) repository
  and cached locally via BiocFileCache.

- [MutSeqRData](/packages/MutSeqRData) Experimental data for use with
  the MutSeqR vignette and examples. This dataset is taken from
  LeBlanc et al., 2022. 24 MutaMouse animals were exposed to one of
  three doses of benzo&#91;a&#93;pyrene or a vehicle control for 28 days by
  oral gavage. 28 days after the end of the exposure, bone marrow of
  the femurs was harvested from euthanized animals. DNA extraction
  was conducted via DNeasy Blood and Tissue kit. DNA samples were
  sequenced using TwinStrand's Duplex Sequencing on the Mouse
  Mutagenesis Panel at > 10,000 depth. The Mouse Mutagenesis Panel
  comprises 20 2.4kb genomic targets with one located on each mouse
  autosome (two on chromosome 1). Pre-processing of sequence reads
  was redone since publication using an updated version of
  TwinStrand's Mutagenesis App (v. 3.20.1) which produced tabular
  mutation data files for each sample. Data contained herein are only
  those required for running MutSeqR examples and vignette.

New Annotation Packages
=====================

There are no new annotation packages in this release of Bioconductor.

New Workflow Packages
=====================

There are no new workflow packages in this release of Bioconductor.

New Online Books
=====================

There are two new online books in this release of Bioconductor.

- [OMA](/packages/OMA) This is a reference cookbook for **Microbiome
  Data Science** with R and Bioconductor.

- [scrapbook](/packages/scrapbook) Demonstrates some basic analysis
  of single-cell RNA-seq data with scrapper. This includes quality
  control, normalization, feature selection, various forms of
  dimensionality reduction, clustering into subpopulations, and
  detection of marker genes. It is intended for users who already
  have some familiarity with R and want to get hands-on with some
  basic single-cell analyses.

NEWS from existing Software Packages
===================================

[ADaCGH2](/packages/ADaCGH2)
-------

                 Changes in version 2.51.3 (2026-04-16)                 

- Fixed minor notes and warnings in r-universe.

                 Changes in version 2.51.2 (2026-04-14)                 

- Fixed "Remove Maintainer field. Use Authors@R &#91;cre&#93; designation"
  error from r-universe.

                 Changes in version 2.51.1 (2026-01-31)                 

- Fixed (most) NOTES from R CMD check --as-cran

[alabaster.base](/packages/alabaster.base)
--------------

                       Changes in version 1.12.0                        

- Added a cloneFile() utility to clone a single file, by copying,
  hard-linking or (relative) symlinking.

[alevinQC](/packages/alevinQC)
--------

                       Changes in version 1.27.1                        

- Remove explicit C++11 system dependency

[anndataR](/packages/anndataR)
--------

                        Changes in version 1.2.0                        

New features

- 
  Add initial Zarr support for reading and writing .zarr stores

- 
  Allow manually setting chunk size for HDF5 writes

- 
  Add auto-chunking for HDF5 writes to improve performance

- 
  Optimise sparse matrix reading performance by avoiding
  Matrix::sparseMatrix and constructing objects manually

Bug fixes

- 
  Handle unnamed SingleCellExperiment assays in as_AnnData() by
  automatically assigning names with a warning

- 
  Improve warnings when items fail to convert between objects and
  add new warnings and checks

Documentation

- 
  Updates to vignettes

- 
  Fix broken links

- 
  Update author details

Testing

- 
  Add a benchmarking system to monitor function performance

- 
  Improvements to roundtrip tests

- 
  Improvements to HDF5 tests

- 
  Improvements to GitHub Actions CI

[annoLinker](/packages/annoLinker)
----------

                       Changes in version 0.99.7                        

- Initial Acceptance to Bioconductor

[AnnotationHubData](/packages/AnnotationHubData)
-----------------

                       Changes in version 1.41.0                        

NEW FEATURES

- 1.41.1 Added Parquet as acceptable source type

[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.24.0                        

NEW FEATURES

- (v 1.23.13) Add package options to bypass hardcoded URLs.

- (v 1.23.12) Use keyring to store auth.json file.

USER VISIBLE CHANGES

- (v 1.23.10) Use sha256 instead of md5sum for service validation.

- (v 1.23.4) Use native pipe operator |> instead of magrittr pipe.

BUG FIXES AND MINOR IMPROVEMENTS

- (v 1.23.11) Update Dockstore API version to 1.19.2.

- (v 1.23.10) Fix typo in Service docs and use rmarkdown backticks
for
code in documentation.

- (v 1.23.8) Use qualified gcloud_access_token to avoid NAMESPACE
conflicts.

- (v 1.23.7) Update APIs and track openapi.yaml in automatic commits.

- (v 1.23.6) Update Dockstore API version to 1.18.2.

- (v 1.23.5) Update R version dependency to R-devel.

- (v 1.23.4) Add code chunk labels to vignettes.

- (v 1.23.3) Add URL and BugReports fields to documentation and
DESCRIPTION.

[AnVILAz](/packages/AnVILAz)
-------

                        Changes in version 1.6.0                        

Bug fixes and minor improvements

- Sanitize error messages in .az_do() to prevent sensitive
information
in URLs.
- Added unit tests for utility functions.

[AnVILBase](/packages/AnVILBase)
---------

                        Changes in version 1.6.0                        

- Use GCPtools::gcloud_exists() for checking gcloud availability.

[AnVILGCP](/packages/AnVILGCP)
--------

                        Changes in version 1.6.0                        

Deprecated and defunct

- Move gcloud functions to DEFUNCT status as they have been fully
migrated to GCPtools

Bug fixes and minor improvements

- Use GCPtools::gcloud_exists for internal checks
- Use with_mocked_bindings in tests

[AnVILPublish](/packages/AnVILPublish)
------------

                       Changes in version 1.22.0                        

New Features

- (v. 1.21.1) Export create_workspace() for creating a new AnVIL
workspace without necessarily populating it from an R package.

[APAlyzer](/packages/APAlyzer)
--------

                 Changes in version 1.25.2 (2026-04-05)                 

- Fixed compatibility with R-devel and Bioconductor 3.22/3.23.

- Fixed test failure caused by ggplot2 4.0+ S7 object type change.

- Replaced deprecated element_line(size=) with
  element_line(linewidth=).

- Replaced library(HybridMTest) calls with requireNamespace(); moved
  HybridMTest to Suggests.

- Migrated makeTxDbFromGRanges import from GenomicFeatures to
  txdbmaker.

- Added txdbmaker, rlang, and S4Vectors to Imports.

- Fixed undefined global variables (Print, queryHits, download.file,
  etc.).

- Added .Rbuildignore to exclude .travis.yml from build.

- Updated License field in DESCRIPTION.

[asuri](/packages/asuri)
-----

                        Changes in version 1.0.0                        

- Launch version

[atacInferCnv](/packages/atacInferCnv)
------------

                        Changes in version 1.0.0                        

- Initial Acceptance to Bioconductor

[atena](/packages/atena)
-----

                 Changes in version 1.18.0 (2026-04-27)                 

BUG FIXES

- Fix in OneCodeToFindThemAll() when the 'dictionary' argument is not
  NULL.

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.23.1 (2026-04-09)                 

- Bug fix to replace deprecated `arma::is_finite` by `std::isfinite`

[BatChef](/packages/BatChef)
-------

                       Changes in version 0.99.12                       

- Initial Acceptance to Bioconductor

[BatchQC](/packages/BatchQC)
-------

                        Changes in version 2.7.2                        

Minor Changes

- Changed maintainer from Jessica Anderson to Yaoan Leng

Minor Changes

- Added lambda example to Intro vignette

[Battlefield](/packages/Battlefield)
-----------

                       Changes in version 0.99.03                       

- Submitted to Bioconductor

[beachmat](/packages/beachmat)
--------

                       Changes in version 2.28.0                        

- Added the .unknown.action= option to initializeCpp(), to
  determine whether to print a message, warning or error when
  using the unknown matrix fallback.

- Consolidated row/column variants of various utilities into a
  single function, i.e., tatami.sums() for row/column sums,
  tatami.nan.counts() for NaN counts.

- Added tatami.medians() to compute the row- and column-wise
  medians.

- Added tatami.sums.by.group() to compute the row- and
  column-wise sums for groups of columns/rows, respectively.

[BEclear](/packages/BEclear)
-------

                 Changes in version 2.27.1 (2026-01-02)                 

- migration from dependency futile.logger to logger

[bedbaser](/packages/bedbaser)
--------

                        Changes in version 1.3.5                        

- Incorporated test_request parameter
- Improved how bedbaser constructs BED file url

                        Changes in version 1.3.1                        

- Update bedbase to 0.12.1
- Replace download.file with curl::curl_download
- Fix test for bb_bed_text_search

[beer](/packages/beer)
----

                 Changes in version 1.15.1 (2025-11-04)                 

- Updated NEWS versioning to reflect Bioconductor updates.

- Updated mislabeled progressr::progress dependency.

[betterChromVAR](/packages/betterChromVAR)
--------------

                        Changes in version 0.0.1                        

- Initial Acceptance to Bioconductor


[BiocAzul](/packages/BiocAzul)
--------

                        Changes in version 1.0.0                        

New features

- Initial Bioconductor submission!

[BiocBuildReporter](/packages/BiocBuildReporter)
-----------------

                       Changes in version 0.99.0                        

SIGNIFICANT NEW FEATURES

- First Bioconductor release.

[BiocCheck](/packages/BiocCheck)
---------

                       Changes in version 1.48.0                        

BUG FIXES AND MINOR IMPROVEMENTS

- Package tarball size check updated to have a max of 10 MB (@LiNk-NY,
  #234).

[BiocMaintainerApp](/packages/BiocMaintainerApp)
-----------------

                       Changes in version 0.99.0                        

SIGNIFICANT NEW FEATURES

- First Bioconductor release.

[BiocNeighbors](/packages/BiocNeighbors)
-------------

                        Changes in version 2.6.0                        

- Enable manual vectorization for Annoy and HNSW, if supported by
  the compiler. This improves speed at the cost of exact
  reproducibility across different toolchains/architectures. We
  consider this trade-off to be acceptable given that many other
  floating-point calculations already have non-identical
  precision across machines.

- Accept non-ordinary matrices in X= and query= arguments for all
  functions. This requires beachmat to be installed.

- Added .check.nonfinite= option to skip checks for non-finite
  values in all functions. This is useful when the user can
  guarantee that missing or infinite values are not present.

- Moved the BNPARAM= argument closer to the front of the argument
  list of buildIndex(). This is more intuitive when examining the
  generic arguments.

- Added saveIndex() and loadIndex() functions to save and load
  indices from disk.

[BiocPkgDash](/packages/BiocPkgDash)
-----------

                       Changes in version 0.99.33                       

- Initial Acceptance to  Bioconductor

[BiocSingular](/packages/BiocSingular)
------------

                       Changes in version 1.28.0                        

- Bugfix for center=TRUE and scale=TRUE when there are fewer than
  two rows.

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.68.0                        

BREAKING CHANGES

- The `port` argument, if unspecified and if no scheme is provided in
  the `host`
  argument, now defaults to `443` (the standard port for `https`),
  instead of
  `80` (the standard port for `http`). This fixes issue #155 reported
  by
  Lori Shepherd.

- For security reasons, biomaRt no longer tries to downgrade SSL
  settings in
  case of HTTPS connection issues. It is expected that most clients and
  servers are
  now set up to deal properly with HTTPS. In case of persistent issues,
  users are
  still able to set the SSL settings manually via `setEnsemblSSL()`.

MINOR IMPROVEMENTS

- All requests sent by biomaRt now include a custom user agent linking
  to the
  source GitHub repository, and indicating the package version.

- The defunct `archive` argument in `listMarts()` and `useMart()` has
  been
  removed. It has been deprecated since Bioconductor 3.7 and defunct
  since
  Bioconductor 3.8 (2018).

INTERNAL CHANGES

- Tests and some internals are now more robust to Ensembl Biomart being
  not
  responsive or unstable.

[biomformat](/packages/biomformat)
----------

                       Changes in version 1.39.17                       

BUG FIXES

- Added explicit S4 dispatch method for biom_data() with signature
  c("biom", "character", "missing"). Previously, calling biom_data()
  with only a character rows= argument (e.g. biom_data(x,
  rows="GG_OTU_3"))
  matched the catch-all ("biom", "character", "ANY") method, which then
  forwarded the still-missing columns argument to the next dispatch
  layer,
  causing the error: "argument 'columns' is missing, with no default".
  The new method defaults columns to 1:ncol(x) and dispatches cleanly.
  Fixes the vignette build ERROR on Bioconductor (chunk
  named_subsetting,
  biomformat.Rmd lines 402-413).
  R CMD check: 0 ERRORs, 0 WARNINGs, 4 pre-existing NOTEs.

                       Changes in version 1.39.16                       

BUG FIXES

- Updated Joseph N. Paulson's email address in Authors@R to
  joseph.paulson@yale.edu.
  R CMD check: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs.

                       Changes in version 1.39.15                       

BUG FIXES

- Replaced deprecated Author:/Maintainer: DESCRIPTION fields with the
  modern Authors@R: format using person() with role=c("aut","cre") for
  the maintainer. This was flagged as an ERROR by BiocCheck >= 1.46 and
  was the root cause of the "invalid email" complaint from Bioconductor
  —
  the field itself was structurally invalid, not just the address
  string.
  Updated maintainer email to joey711@gmail.com (reachable address).
  Added .BiocCheck to .Rbuildignore.
  R CMD check: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs.
  BiocCheck: 1 ERROR (Watched Tags on support site — manual web
  action),
  3 WARNINGs (version parity, no BioC deps, missing \value),
  14 NOTEs (style/cosmetic).

                       Changes in version 1.39.14                       

DOCUMENTATION / VIGNETTE

- Added two new vignette sections to vignettes/biomformat.Rmd:
  * "Constructing a BIOM from R data": end-to-end make_biom() workflow
  showing how to build a biom object from a count matrix and a
  data.frame taxonomy table with list-valued columns (the dada2
  pattern),
  write it with write_biom(), and read it back. Includes a callout
  directing large-dataset users to write_hdf5_biom() (addressing Issue
  #8), with a guarded HDF5 code chunk. Directly addresses Issues #4,
  #6, #9 for users who land on the vignette.
  * "Subsetting biom_data() by name": demonstrates the character-vector
  row/column subsetting interface of biom_data() that was previously
  undocumented in the vignette.
  R CMD check: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs.

                       Changes in version 1.39.13                       

DOCUMENTATION

- write_biom(): added @details section documenting the 2^31-1 byte
  size limitation that arises because jsonlite serialises the entire
  BIOM object to a single R character string. Users with large datasets
  (thousands of samples or features) are now explicitly directed to
  write_hdf5_biom() which has no such size constraint.
  Closes GitHub Issue #8.
  R CMD check: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs.

                       Changes in version 1.39.12                       

NEW FEATURES / TESTS

- Issue #7 regression test: added inst/extdata/zero_col_hdf5.biom, a
  minimal 3x3 HDF5 BIOM fixture where the middle sample ("ZeroSamp")
  has
  all-zero counts. The fixture is generated by the companion script
  inst/extdata/create_zero_col_hdf5.R using rhdf5 directly.
  Added two tests in test-hdf5-write.R:
  * "zero-column HDF5 fixture: biom_data() returns correct 3x3 matrix"
  —
  verifies that ZeroSamp is all zeros and the other two columns have
  the expected non-zero values. Validates that generate_matrix()
  (rewritten in v1.39.8 to use sparseMatrix directly from CCS triplets)
  correctly handles the case where indptr&#91;j&#93; == indptr&#91;j+1&#93;.
  * "zero-column HDF5 fixture: make_biom() + write_hdf5_biom()
  round-trip"
  — verifies that write_hdf5_biom() -> read_biom() preserves the
  all-zero column.
  Both tests are guarded with skip_if_not_installed("rhdf5").
  R CMD check: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs.

                       Changes in version 1.39.11                       

BUG FIXES

- make_biom(): fix Issue #4 — NULL id was serialised as {} (empty JSON
  object)
  by jsonlite. make_biom() now substitutes "No Table ID" when id is
  NULL,
  matching the BIOM spec and ensuring write_biom() -> read_biom() is a
  lossless round-trip. validObject() now succeeds on the round-tripped
  object.

- make_biom(): fix Issue #6 — when observation_metadata (or
  sample_metadata)
  is a data.frame with list columns (e.g. a "taxonomy" column holding
  character vectors of rank assignments, as produced by dada2), the
  metadata
  was serialised as a bare JSON array (&#91;&#91;...&#93;&#93;) instead of a named JSON
  object
  ({"taxonomy":&#91;...&#93;}). The BIOM spec and downstream tools (phyloseq
  import_biom(), the Python biom library) require the named-object
  form.
  Root cause: as.matrix() on a list-column data.frame produces a
  list-matrix;
  as.list(row) then collapses field names. Fix: detect list-column
  data.frames
  and build per-row metadata directly as named lists, bypassing
  as.matrix().
  Closes GitHub Issue #6. Also resolves user-support Issue #9 (dada2 ->
  biom
  -> write_biom workflow now works correctly).
  R CMD check: 0 ERRORs, 0 NOTEs, 2 pre-existing vignette WARNINGs.

                       Changes in version 1.39.10                       

NEW FEATURES

- Updated vignette with four new sections:
  * HDF5 (BIOM v2) read and write via write_hdf5_biom() / read_biom(),
  including JSON-to-HDF5 conversion.
  * Tidy long-format output via as.data.frame() and as_tibble.biom(),
  with purrr-style summarisation examples (Shannon diversity,
  per-sample total counts) and base-R fallbacks.
  * SummarizedExperiment interoperability via
  biom_to_SummarizedExperiment() and as(x, "SummarizedExperiment"),
  showing assay(), colData(), and rowData() access.
  * Session info section.

- All new vignette sections are guarded with requireNamespace() so
  the vignette builds cleanly without optional dependencies
  (rhdf5, tibble, purrr/dplyr, SummarizedExperiment/S4Vectors).

                       Changes in version 1.39.9                        

NEW FEATURES

- write_hdf5_biom(x, biom_file): new exported function that serialises
  a
  biom object to the BIOM v2 HDF5 format.  Writes both the sample-major
  and observation-major compressed-sparse representations required by
  the
  spec, plus all sample and observation metadata.  Requires rhdf5
  (Bioconductor); a clear error is raised if it is absent.
  read_biom() -> write_hdf5_biom() -> read_biom() is a lossless
  round-trip
  for both count data and metadata.

BUG FIXES

- Moved rhdf5 from Imports to Suggests.  The package loads and all JSON
  BIOM functionality works without rhdf5 installed; HDF5 read/write
  simply
  stops with an informative message when rhdf5 is absent.

                       Changes in version 1.39.8                        

BUG FIXES / PERFORMANCE

- generate_matrix(): rewrote HDF5/BIOM-v2 matrix reconstruction to
  build
  a sparse Matrix directly from the CCS (indptr/indices/data) triplets
  stored in the HDF5 file, instead of first constructing a dense
  base::matrix via sapply() and then converting.  For large datasets
  this
  avoids an O(n_obs * n_samples) allocation.  The return value (list of
  named vectors, one per observation) is unchanged so all downstream
  code
  is unaffected.  Also handles the edge case of an all-zero matrix
  (length(data) == 0) explicitly.

                       Changes in version 1.39.7                        

NEW FEATURES

- as.data.frame.biom(): new S3 method that converts a biom object to a
  long-format (tidy) data.frame with one row per (feature, sample)
  pair.
  Columns: feature_id, sample_id, count, plus any sample and
  observation
  metadata columns appended via left-join. Pure base R, no tidyverse
  dependency.

- as_tibble.biom(): thin wrapper around as.data.frame.biom() that
  returns
  a tibble. Requires the 'tibble' package (Suggests only). Call via
  tibble::as_tibble(x) or as_tibble.biom(x) directly.

DEPENDENCY CHANGES

- Added tibble to Suggests (optional; only needed for
  as_tibble.biom()).

                       Changes in version 1.39.6                        

NEW FEATURES

- biom_to_SummarizedExperiment(): new exported function that converts a
  biom object into a SummarizedExperiment, placing the count/value
  matrix
  in assay("counts"), sample metadata in colData(), and feature
  metadata
  in rowData(). Both colData and rowData are S4Vectors::DataFrame
  objects.
  When a biom object carries no metadata (the accessor returns NULL),
  an
  empty DataFrame with correct row/col names is used, ensuring the SE
  is
  always valid. No hard dependency is introduced: SummarizedExperiment
  and
  S4Vectors are listed in Suggests only.

- as(x, "SummarizedExperiment"): S4 coercion method registered at load
  time when SummarizedExperiment is available, delegating to
  biom_to_SummarizedExperiment().

DEPENDENCY CHANGES

- Added SummarizedExperiment and S4Vectors to Suggests.

TESTS

- New tests/testthat/test-SE.R with 6 tests covering: return class,
  assay content, colData content, rowData content, S4 coercion, NULL
  metadata, and SE dimension/dimname correctness.
  All tests are skipped gracefully when SummarizedExperiment is not
  installed.

                       Changes in version 1.39.5                        

DEPENDENCY CHANGES

- Removed plyr (>= 1.8) from Imports entirely. plyr is unmaintained and
  every use in this codebase now has a direct base-R equivalent.
  This eliminates an unmaintained dependency and reduces install
  footprint.

- Bumped R dependency from >= 3.2 to >= 4.1, ensuring modern base-R
  idioms
  (including the native pipe |>) are available.

- Replaced import(Matrix) (whole-namespace) with selective
  importFrom(Matrix, Matrix, sparseMatrix, drop0) in NAMESPACE.
  Replaced import(methods) with selective importFrom(methods, ...).
  Follows CRAN/Bioconductor best practices; prevents namespace
  pollution.

USER-VISIBLE CHANGES

- biom_data(), sample_metadata(), observation_metadata(): the parallel=
  argument is now a no-op with a deprecation warning when passed as
  TRUE.
  The plyr-backed parallel execution it previously enabled no longer
  exists.
  Existing code that passes parallel=FALSE (the default) is unaffected.

INTERNAL CHANGES

- make_biom(): replaced plyr::alply() with lapply(seq_len(nrow(...)))
  for building per-row named metadata lists.

- biom_data() dense path: replaced plyr::laply() with
  do.call(rbind, lapply(...)).

- biom_data() sparse numeric path: replaced plyr::ldply(x$data) with
  do.call(rbind, lapply(x$data, as.data.frame)).

- biom_data() sparse unicode path: replaced plyr::ldply(x$data,
  function...)
  with do.call(rbind, lapply(x$data, function(e) ...)).

- extract_metadata(): replaced plyr::llply()/plyr::ldply() with
  lapply() / do.call(rbind, lapply(...)).

                       Changes in version 1.39.4                        

BUG FIXES

- biom_data(): fixed data-corruption bug on both dense and sparse BIOM
  paths
  where subsetting to a single row or single column silently collapsed
  the
  result into a dimensionless named vector, discarding dim(),
  rownames(), and
  colnames(). This caused downstream tools (notably
  phyloseq::import_biom())
  to fail silently or produce incorrect OTU tables.
  Fix: on the sparse path, added drop = FALSE to the matrix subsetting
  call
  (m&#91;rows, columns, drop = FALSE&#93;). On the dense path, the laply()
  result is
  now immediately reshaped with matrix(m, nrow, ncol) before Matrix()
  coercion,
  ensuring a 2-D object is always returned regardless of dimension
  lengths.
  Closes GitHub PR #12 (https://github.com/joey711/biomformat/pull/12):
  "Fix biom_data() when dealing with 1-taxon and 1-sample BIOM data"
  Supersedes GitHub PR #11
  (https://github.com/joey711/biomformat/pull/11):
  "Fix unidentical biom output by make_biom()"

- biom_data(): simplified the post-subsetting naming block. Both paths
  now
  always produce a 2-D object, so rownames() and colnames() are applied
  unconditionally on all code paths (the previous is.null(dim(m))
  branch
  is no longer needed).

TESTS

- Added regression tests in tests/testthat/test-IO.R covering all new
  behaviour introduced across v1.39.2-v1.39.4:
  * sparse single-row subset retains dim() and dimnames (PR #12)
  * sparse single-col subset retains dim() and dimnames (PR #12)
  * sparse single-cell (1x1) subset retains dim()
  * full sparse matrix unaffected by drop = FALSE fix
  * read_biom() routes HDF5 fixture cleanly without jsonlite warning
  (Issue #14, PR #16)
  * read_biom() correctly classifies all JSON and HDF5 extdata fixtures

                       Changes in version 1.39.3                        

BUG FIXES

- read_biom(): replaced the fragile JSON-first / HDF5-fallback try()
  chain
  with a deterministic magic-bytes router. The function now reads the
  first
  4 bytes of the file; if the HDF5 signature (\x89 H D F) is detected
  it
  routes exclusively to read_hdf5_biom() and never invokes jsonlite.
  JSON
  files route exclusively to jsonlite. This eliminates the confusing
  "lexical error: invalid char in json text ... <89>HDF" warning that
  users
  encountered when HDF5 files were accidentally passed through the JSON
  parser first.
  Closes GitHub Issue #14
  (https://github.com/joey711/biomformat/issues/14):
  "Unable to read HDF5 biom file"
  Supersedes GitHub PR #16
  (https://github.com/joey711/biomformat/pull/16):
  "Improve handling of HDF5 BIOM files"
  Partially addresses GitHub Issue #5
  (https://github.com/joey711/biomformat/issues/5)
  and Issue #3 (https://github.com/joey711/biomformat/issues/3):
  fatal C-level aborts when reading large or malformed HDF5 files are
  now
  caught by the new tryCatch() wrapper in read_hdf5_biom() and
  re-emitted
  as informative R-level warnings instead of crashing the session.

- read_hdf5_biom(): added requireNamespace("rhdf5", quietly = TRUE)
  guard.
  If the rhdf5 package is not installed, or if the underlying HDF5
  system
  libraries are absent (common on stripped BBS nodes or end-user
  machines
  without libhdf5), the function now emits a clear R-level warning
  identifying the missing dependency and returns NULL invisibly,
  instead of
  producing a fatal C-level abort.

- read_hdf5_biom(): wrapped the h5read() call in tryCatch() so that any
  C-level or system-library error is caught and re-emitted as an
  informative
  R warning, keeping the R session alive.

DEPENDENCY CHANGES

- Bumped Matrix dependency from (>= 1.2) to (>= 1.7-0) in DESCRIPTION.
  This pins the package against the post-SuiteSparse ABI break and
  prevents
  runtime crashes caused by binary incompatibilities in upstream
  sparse-matrix
  libraries on BBS nodes running R 4.4+.

                       Changes in version 1.39.2                        

BUG FIXES

- Fixed fatal test ERROR under testthat >= 3.0.0: replaced all
  deprecated
  expect_that(x, is_true()), expect_that(x, is_identical_to(y)), and
  expect_that(x, is_a("cls")) calls in tests/testthat/test-IO.R with
  their
  modern equivalents (expect_true(), expect_identical(), expect_is(),
  expect_true(is(x, "cls"))). The removed helpers caused the entire
  test
  suite to ERROR on current Bioconductor BBS nodes, which was the
  primary
  trigger for the deprecation warning.
  Addresses GitHub Issue #17
  (https://github.com/joey711/biomformat/issues/17):
  "Bioconductor failure and risk of deprecation"

- Added missing importFrom(stats, setNames) and importFrom(utils,
  packageVersion) directives to NAMESPACE, resolving "no visible global
  function definition" NOTEs from R CMD check.

[blase](/packages/blase)
-----

                        Changes in version 1.2.0                        

Major changes

- "confident_mapping" has now been renamed as "strong_mapping" to
indicate that it is not a robust statistical confidence. The
confident_mapping() getter function has been renamed to
strong_mapping()
- plot_mapping_result_heatmap() annotate_confidence parameter has
been
renamed to annotate_strong.

Minor improvements and bug fixes

- Now possible to use BLASE to create pseudotime bins from multiple
trajectories at once.
- Speeds up by cell pseudotime bin assignment.
- Possible to apply different metrics for judging the best mapping.
Spearman is still recommended.
- It is now possible to rename a bulk sample for a MappingResult
object.
- MappingResult heatmaps no longer convert the Y-axis titles (bulk
sample name) to character, allowing numeric style scales.

[bnbc](/packages/bnbc)
----

                        Changes in version 1.33                         

- Replaced HiCBricks with direct rhdf5 calls in getChrCGFromCools
  to fix a segfault in HDF5 cleanup that was crashing the
  vignette build. HiCBricks is no longer a dependency.

[BreastSubtypeR](/packages/BreastSubtypeR)
--------------

                        Changes in version 1.3.2                        

Highlights (from v1.1.3 onward)

- Paper published in NAR Genomics and Bioinformatics (2025), Editor’s
Choice (DOI: 10.1093/nargab/lqaf131).
- Support for raw RNA-seq counts (requires gene lengths).
- iBreastSubtypeR refresh: cleaner UX, smarter AUTO guidance,
consistent exports.

Enhancements

- ssBC/ssBC.v2: singleton subgroup robustness: Subgroups with n=1 no
longer error:
- Keeps matrix shape (drop=FALSE) and hardens dimnames/types.
- Primary path: original sspPredict().
Fallback: nearest-centroid (Spearman) when needed.
- If there are 0 common PAM50 genes, returns NA labels with shaped
distances/dist.RORSubtype to avoid downstream errors.
- ROR computation guarded for incomplete inputs.
- SSPBC output now "full": BS_sspbc() and Shiny "sspbc" runs return a
full metrics table (not calls-only).
- Exports map core label columns to the standard names
(Call_5class / Call_4class when applicable).
- Shiny: "Load example data…" button
- One-click load of a small demo dataset from inst/RshinyTest/ to
explore the UI without uploads.\
- Shows a notification on success; users can immediately run
Preprocess & map and analyses.
- AUTO preflight UI (Shiny): Now detects cohort kind (TN, ER/HER2,
ER-only, HER2-only) and shows compact stats:
- ER/HER2 subgroups: ER+/HER2-, ER-/HER2-, ER+/HER2+, ER-/HER2+\
- TN cohorts: TN vs nonTN\
- Readiness uses the same minimums used by AUTO (sourced
programmatically; no duplicated thresholds).
- Shorter notifications.
-Routine toasts (e.g., "Step 1 complete. Proceed to Step 2.") now
auto-dismiss sooner to reduce UI clutter.
- Phenodata normalization (Mapping): Accepts flexible ER/HER2/TN
encodings and normalizes to canonical forms (ER+/ER-, HER2+/HER2-,
TN/nonTN). Ambiguous HER2="2+" remains as-is and raises a warning.

Bug fixes

- TN cohorts + ssBC: BS_Multi() now respects TN cohorts when methods
are specified manually; ssBC/ssBC.v2 switch to s = "TN" / "TN.v2"
when a TN column indicates a TN cohort. Falls back to s = "ER" /
"ER.v2" otherwise.
- AUTO: Fixed a crash in BS_Multi(methods="AUTO") when ER and/or HER2
contained missing values (NA).
- AUTO internals: fixed variable name typo (samples_ERHER2.icd).
- Mapping(): Robust ENTREZID coercion (from as.character() to
as.integer() with suppressed warnings).
- cIHC.itr: outList$distances now returned as numeric matrix.

Shiny

- Surface method warnings as toasts:
- Runs are wrapped in a warning handler; package warnings (e.g.,
ssBC.v2 singleton fallbacks) appear as yellow notifications.\
- Warnings include subgroup, n, and example sample IDs for quick
triage.
- Shiny preflight reset:
- Fixed a stale cohort summary after switching data sources
(manual uploads ↔ example). The preflight panel now revalidates
once inputs change.

Developer notes

- Added lightweight internal logger ._msg() and replaced scattered
message() calls in AUTO to standardize package output without
affecting CRAN/Bioc checks.

Documentation

- README/vignette: brief note on the example-data button and expected
file locations.
- Mapping(): Column metadata clarified. Added explicit requirements
for receptor fields used by AUTO and ER/HER2/TN-dependent methods
(ssBC, cIHC/cIHC.itr, PCAPAM50) and for ROR covariates (TSIZE, NODE
as numeric 0/1). Documented preferred coding and automatic
normalization behavior.

Compatibility Notes

- SSPBC "full" output keeps previous columns for calls; additional
metrics may appear.

Upgrade Notes

- Raw RNA-seq counts are supported from v1.1.3 onward (requires gene
lengths).
- If you previously parsed BS / BS.Subtype, switch to Call_5class /
Call_4class.
- Package API unchanged.

[bsseq](/packages/bsseq)
-----

                        Changes in version 1.47                         

- Making `read.modkit()` and `read.modbam2bed()` defunct.

[CalibraCurve](/packages/CalibraCurve)
------------

                        Changes in version 1.1.4                        

- readDataTable() can now import a file containing data for multiple
substances

                        Changes in version 1.1.3                        

- small fix in documentation

                        Changes in version 1.1.2                        

- Allow calculation of calibration curves with only one replicate per
concentration level.
- Fix bug when input data is not ordered by concentration column
- Allow legend and add title to response factor plots
- Add new Logo to readme
- Fix bug when unequal numbers of replicates per concentration level
are present

                        Changes in version 1.1.1                        

Bugfixes

[Cardinal](/packages/Cardinal)
--------

                 Changes in version 3.13.4 (2026-04-24)                 

Major changes

- Converted NEWS to NEWS.md

Major changes

- Removed deprecated functions (color palettes, etc.)
- Removed defunct functions (color palettes, etc.)
- Deprecated 'fetch()' and 'flash()'

                 Changes in version 3.13.2 (2026-04-22)                 

Bug fixes

- Fix 'mass.range' bug in 'convertMSImagingExperiment2Arrays()'

                 Changes in version 3.13.1 (2026-03-18)                 

Major changes

- Remove dependency on EBImage

                 Changes in version 3.12.1 (2026-02-12)                 

Bug fixes

- Add na.rm parameter to 'meansTest()'
- Fix spectrumRepresentation logic in 'writeImzML()'
- Rename 'bi' to 'bilateral' in smoothing functions

[carnation](/packages/carnation)
---------

                       Changes in version 0.99.10                       

- carnation submitted to bioconductor

[CDI](/packages/CDI)
---

                        Changes in version 1.9.0                        

- Package accepted by Bioconductor (2023-06-30)

[CellMentor](/packages/CellMentor)
----------

                 Changes in version 0.99.1 (2024-10-27)                 

- Initial Acceptance to Bioconductor

[cellmig](/packages/cellmig)
-------

                        Changes in version 1.1.8                        

- Unit tests fixed for Bioconductor release

                        Changes in version 1.1.7                        

- Stable version used for publication


[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.45.3                        

- Fix the issue for biomaRt in vignettes.

                       Changes in version 3.45.2                        

- Remove microRNAs from annoGR

                       Changes in version 3.45.1                        

- Remove microRNAs from assignChromosomeRegion.

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.47.1                        

- fixed issue in 'test-txdb.R' as 'TxDb.Hsapiens.UCSC.hg19.knownGene'
changes its transcript ID from UCSC (e.g., uc002qsd.4) to Ensembl
(e.g., ENST00000487630.1_3) (2025-11-04, Tue)

[Chromatograms](/packages/Chromatograms)
-------------

                         Changes in version 1.1                         

Changes in 1.1.8

- Improve performance of matchRtime().

- Fix compareChromatograms(): ... arguments (e.g. tolerance) are now
routed to MAPFUN or FUN based on their formal parameters, preventing
errors when FUN = cor received unknown arguments.

Changes in 1.1.7

- Improve performance of .prepare_spectra_input(): spectra are now
pre-filtered to the non-overlapping union of EIC retention time
ranges using MsCoreUtils::reduce() and Spectra::filterRanges(), and
peak data is loaded in a single peaksData() call instead of separate
mz() and intensity() calls. This reduces I/O and memory usage,
especially for file-backed backends.

Changes in 1.1.6

- Add compareChromatograms() and matchRtime() for pairwise similarity
of chromatographic intensity profiles.

Changes in 1.1.5

- Add peakBoundary() method for Chromatograms objects. Determines the
retention time boundaries of the tallest peak in each chromatogram
using MsCoreUtils::valleys() to locate flanking valleys, with a
threshold-based fallback. Returns a matrix with left_boundary and
right_boundary columns.

Changes in 1.1.4

- setBackend() now clears the processing queue after switching
backend, preventing queued processing steps from being applied twice
(once during the data transfer and again on subsequent peaksData()
calls).

- Fix setBackend() parallel branch (used for ChromBackendMzR) to
correctly apply queued processing steps to each chunk before
transferring data to the new backend.

- Major performance improvement in .process_peaks_data() for
ChromBackendSpectra. Key optimizations include: pre-extracting mz()
and intensity() as plain R lists (avoiding slow SimpleNumericList
indexing), global retention time pre-filtering, a fast path for
TIC/BPC cases, and using findInterval() with cumsum() for m/z range
lookups. Combined, setBackend() showed 9x speed up for 1000
chromatograms.

- Add optimized intensity() and rtime() accessors for
ChromBackendMemory using direct [[ extraction instead of the slower
&#91;, col, drop&#93; path through peaksData().

- Further accessor optimizations: intensity(), rtime(), and lengths()
on Chromatograms now bypass peaksData() dispatch when the processing
queue is empty. peaksData() on ChromBackendMemory uses a fast [[
path for single-column requests. Direct lengths() methods added for
all backends using nrow() instead of going through intensity().

- Replace do.call(rbind, ...) with data.table::rbindlist() in
chromExtract() for ChromBackendMemory, ChromBackendMzR, and
ChromBackendSpectra for faster row-binding of many data.frames.

- Replace replicate(n, .EMPTY_PEAKS_DATA, simplify = FALSE) with
rep(list(.EMPTY_PEAKS_DATA), n) across backends to avoid repeated
expression evaluation overhead.

Changes in 1.1.3

- Add filterEmptyChromatograms() function to remove empty
chromatograms (i.e., chromatograms without peaks) from a
Chromatograms or ChromBackend object.

- Add concatenateChromatograms() function and c() method to combine
multiple Chromatograms objects into a single object. Also add
split() method to split a Chromatograms object based on a grouping
factor.

- Add extrapolate parameter to imputePeaksData() (default FALSE).
When
TRUE, leading/trailing NA values outside the range of observed data
are extrapolated. When FALSE (default), only interpolation is
performed and edge NA values remain as NA.

Changes in 1.1.2

- Fix peaksData() for ChromBackendSpectra to return data in the
correct row order when multiple chromatograms share the same
chromSpectraIndex. This bug caused setBackend() to produce
mismatched chromData and peaksData when converting from
ChromBackendSpectra to ChromBackendMemory with objects containing
multiple EICs.

Changes in 1.1.1

- Aligned the package with the Bioconductor 3.22 release.
- Expanded the vignette to cover ChromBackendSpectra usage,
chromatogram extraction with chromExtract(), and imputation
workflows via imputePeaksData().
- Added spectraSortIndex() for ChromBackendSpectra to compute the
desired retention-time order on demand, avoiding the need to keep
on-disk Spectra objects sorted in memory.

[ClonalSim](/packages/ClonalSim)
---------

                       Changes in version 0.99.6                        

Initial Acceptance to Bioconductor

[ClusterGVis](/packages/ClusterGVis)
-----------

                 Changes in version 0.99.0 (2025-08-01)                 

Initial Acceptance to Bioconductor

[clusterProfiler](/packages/clusterProfiler)
---------------

                       Changes in version 4.19.8                        

- fix: map non-ENTREZID universe in enrichGO (2026-04-22, Wed)

                       Changes in version 4.19.7                        

- interpret(), interpret_agent(), and interpret_hierarchical() now
use
aisdk's global default model when model = NULL, so users can switch
the package-wide default with aisdk::set_model() while still
overriding per call with an explicit model argument (2026-03-31,
Tue)

                       Changes in version 4.19.6                        

- update ko2name() to robustly parse KO names via KEGG REST, support
vector input with deduplication, and return NA when NAME is missing
(2026-02-25, Wed)
- bug fixed for plot.interpret (2026-02-05, Thu)

                       Changes in version 4.19.5                        

- interpret() prompt optimized with 'Comparative Analysis' and 'Rule
of Exclusion' to better distinguish cell types with shared functions
(e.g. NK vs CD8+ T cells) using specific marker genes (2026-01-22,
Thu)
- fixed a bug in interpret() where empty interpretation results were
returned due to incorrect list structure handling in
process_enrichment_input (2026-01-22, Thu)
- interpret() now considers specific marker genes in cell type
annotation to avoid key markers being overshadowed by general
pathways (2026-01-21, Wed)
- optimize enrichGO() to avoid memory boom when keyType is not
ENTREZID (2026-01-21, Wed, #805)
- gson_GO_local() to support local GO annotation by adding ancestral
terms (2026-01-21, Wed)
- interpret() implements a gene-based fallback mode for clusters with
no enriched pathways, ensuring comprehensive analysis (2026-01-20,
Tue)
- plot() method for interpretation object to visualize the
LLM-inferred regulatory network using ggtangle (2026-01-20, Tue)
- interpret_agent() supports multi-agent system (Deep Mode) for
interpretation (2026-01-20, Tue)
- Agent Cleaner: Filters noise and selects relevant pathways
- Agent Detective: Identifies key regulators and functional
modules using PPI/TF data
- Agent Synthesizer: Synthesizes findings into a coherent
narrative
- interpret() supports 'Knowledge-Guided Interpretation' (2026-01-20,
Tue)
- add_ppi parameter to integrate PPI network and identify hub
genes
- gene_fold_change parameter to incorporate expression levels
- Mixed-source enrichment analysis support (e.g. Pathways + TFs)
for causal integration
- LLM-guided network refinement to output core regulatory networks
- interpret_hierarchical() for hierarchical interpretation (e.g.
Major
-> Minor clusters) (2026-01-20, Tue)
- interpret() now supports prior parameter for reference-guided
interpretation (e.g. from SingleR/scGPT) (2026-01-20, Tue)

                       Changes in version 4.19.4                        

- interpret() now supports task parameter to specify the task:
'interpretation', 'annotation' and 'phenotyping' (2025-01-18, Sat)
- interpret() supports enrichResult, gseaResult, compareClusterResult
and list of enrichment results (2025-01-18, Sat)

                       Changes in version 4.19.3                        

- instead of packing KEGG cache data in the package, we now download
it from https://yulab-smu.top/clusterProfiler (2025-12-15, Mon)
- add github action to automatically update KEGG cache data
(2025-12-09, Tue)
- use 'enrichit' as engine for enrichment analysis (2025-12-07, Sun)

                       Changes in version 4.19.2                        

- use 'quarto' as vignette engine (2025-11-20, Thu)
- update KEGG cache data (2025-11-20, Thu, #792)
- number of KEGG pathway with category information: 582
- Number of species: 11344

                       Changes in version 4.19.1                        

- bug fixed in enrichPC (2025-11-01, Sat, #789)

[ClustIRR](/packages/ClustIRR)
--------

                       Changes in version 1.9.25                        

- move from CRAN's blaster to Bioconductors rBLAST

                       Changes in version 1.9.17                        

- quantiative comparison community abundance with normalized rank
  abundance
  distributions (NRADs) implemented in get_nrad()

                       Changes in version 1.9.10                        

- bigfixes in joint graph building
  * when clustering both CDR3a/CDR3b, CDR3b edges between-repertoires
  got deleted
  * one between-repertoire edge between two clones with double edges
  was deleted

                        Changes in version 1.9.9                        

- get_blosum() bugfix -> identical CDR3 did not form edges in certain
  cases

- added test-get_blosum.R

                        Changes in version 1.9.2                        

- K-nearest-neighbor approximation implemented

- get_cosine_similarity to compare community abundance overlap

[CNVMetrics](/packages/CNVMetrics)
----------

                       Changes in version 1.15.2                        

BUG FIXES

- The parameter names of the generic functions have been updated to
fit requirement.

                       Changes in version 1.15.1                        

NEW FEATURES

- The documentation has been updated.

[CompensAID](/packages/CompensAID)
----------

                       Changes in version 0.99.7                        

- Initial Acceptance to Bioconductor

[CompoundDb](/packages/CompoundDb)
----------

                        Changes in version 1.15                         

Changes in version 1.15.4

- Support filling missing columns/spectra variables in HMDB import.

Changes in version 1.15.3

- Use data.table::rbindlist() to combine MS spectra for HMDB import
improving its performance.

Changes in version 1.15.2

- Change internal mapping of the precursorIntensity spectra variable
to a database table name from "precursor_intensity" to
"precursorIntensity".

Changes in version 1.15.1

- Add addJoinDefinition() function to define relationships between
core database tables and newly added tables.

[COTAN](/packages/COTAN)
-----

                       Changes in version 2.11.4                        

Added new vignette for Differential Expression Analysis

Fixed minor bug in clustersMarkersHeatmapPlot() about genes
duplication
in markers' list

Refactored vignette into 3 separate files for improved theme focus

Allowed the function calculatePValue() to run on multiple cores

Added parameter minimumUTClusterSize to cellsUniformClustering()
function

                       Changes in version 2.11.3                        

Solved various linting issues and BiocCheck warning

Now the function canUseTorch() is exported

                       Changes in version 2.11.2                        

Improved conversions to/from SingleCellExperiment objects

Solved inconsistencies in cell/genes names between raw data and
meta-data

Improvements for plot functions:

- fixed label placements genes' plot in cleanPlots()
- added genesPercentagePlot() as generalization of
mitochondrialPercentagePlot()

                       Changes in version 2.11.1                        

Fixed issue with isCoexAvailable() for cases when the global
meta-data
have multiple columns

Fixed issue with dispersion solver

Dropped dependency from gghalves package as it is deprecated: using
ggdist instead

                       Changes in version 2.11.0                        

Moved to new version

[CPSM](/packages/CPSM)
----

                        Changes in version 1.3.1                        

Updates

- Enhanced Lasso_PI_scores_f.R to compute the feature Stability Index
and updated Univariate_sig_features_f.R to include the Proportional
Hazards (PH) violation test for selected features.
- Updated MTLR_pred_model_f.R and predict_survival_risk_group_f.R to
generate results for each cross-validation fold.
- Improved train_test_normalization_f.R by adding an additional
filter
to remove features with zero variance in at least 80% of training
samples.
- Revised vignettes/CPSM.Rmd with improved descriptions and a clearer
workflow.
- Updated tests/testthat files to reflect recent changes.
- Updated the DESCRIPTION file:
- Version bumped from 1.3.0 to 1.3.1.

                        Changes in version 1.3.0                        

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
- Bumped version from 1.2.0 to 1.3.9.

[CrcBiomeScreen](/packages/CrcBiomeScreen)
--------------

                 Changes in version 0.99.0 (2025-10-21)                 

Overview

CrcBiomeScreen is an R package designed to streamline
microbiome-based
colorectal cancer (CRC) screening workflows.
It provides standardized functions for preprocessing, taxonomic data
handling, machine learning model training,
and cross-cohort validation — supporting reproducible and
interpretable
microbiome analysis for biomarker discovery.

This version marks the first public release of the package, submitted
to
Bioconductor.

New Features

Data Processing

- SplitTaxas() now automatically detects taxonomy separators
(supports |, ., _, and ;) to handle common formats from MetaPhlAn,
QIIME, etc.
Adds:

- OriginalTaxa column to retain the raw taxonomy string.
- Special handling for uncultured and unclassified taxa
(e.g. converts f__Rikenellaceae|g__unclassified →
Rikenellaceae_unclassified).

- KeepTaxonomicLevel() filters data at a user-defined rank (e.g.,
Genus, Family, Species)
and automatically collapses lower-level abundances.
Handles multiple nested unclassified levels (e.g.,
D_2__Clostridia.D_3__uncultured.D_4__uncultured).

Machine Learning

- TrainModel() serves as a unified interface for training Random
Forest (RF) and XGBoost classifiers.

- Integrates with internal preprocessing and class-weight
balancing.
- Uses withr::with_seed() for local reproducibility instead of
set.seed().
- Adds support for both training and external validation
workflows.

- EvaluateRF() and EvaluateXGBoost() now return standardized
performance metrics
(AUC, accuracy, recall, F1) and store model outputs in the main
object structure.

- ValidateModelOnData() supports model evaluation across independent
datasets.

Data Object Structure

- Introduced the CrcBiomeScreenObject class to store:
- Absolute and relative abundance tables
- Taxonomic annotations
- Model training results
- Validation data and predictions
- Visual outputs (e.g., ROC, variable importance)

This ensures data provenance and reproducibility across the full
workflow.

Utility and Visualization

- Added built-in plotting functions for:
- Model ROC curves (PlotAUC)

Documentation and Testing

- Added comprehensive vignette:

- Demonstrates preprocessing, genus-level filtering, model
training, and validation.
- Includes examples for multiple taxonomy formats and classifier
comparisons.

- Implemented unit tests under tests/testthat/ for key components:

- Taxa splitting
- Level filtering
- Model reproducibility

Technical and Compliance Updates

- Adopted MIT license.
- Added BugReports and URL fields in DESCRIPTION.
- Removed unnecessary system files (.Rproj, .DS_Store).
- Replaced all instances of set.seed() with withr::with_seed() for
Bioconductor compliance.
- Reduced hard-coded dependencies and improved optional imports.

Future Plans

- Add support for additional classifiers.
- Integrate feature interpretation.
- Provide reproducible benchmarking across public CRC datasets.

Version 0.99.13

- Fixed vignette and example issues
- Temporarily disabled XGBoost execution due to compatibility issues
with caret
- Improved documentation and S4 accessor usage

Version 0.99.14(2026-04-10)

- runnable toy workflow added to vignette
- examples updated for TrainModels, ValidateModelOnData, RunScreening
- documentation for included datasets completed

Version 0.99.15(2026-04-14)

- Standardized the format of the vignette and documentation

Maintainer: Chengxin Li (University of Leeds)
Date: 2025-10-21
License: MIT
Repository: https://github.com/omicsForestry/CrcBiomeScreen

[crisprBwa](/packages/crisprBwa)
---------

                       Changes in version 1.15.1                        

- Removed bbsoptions file

[crisprDesign](/packages/crisprDesign)
------------

                       Changes in version 1.13.10                       

- Added optional argument substitutions to addEditedAlleles"

                       Changes in version 1.13.9                        

- Fixed condaEnv argument for addOnTargetScores

                       Changes in version 1.13.8                        

- SNPs can now be annotated from a GRanges object.

                       Changes in version 1.13.6                        

- Fixed splicing junction aa annotation when low scores

                       Changes in version 1.13.4                        

- Fixed the accessor editedAlleles.

                       Changes in version 1.13.1                        

- Added more annotation columns to the addEditedAlleles function.

[crisprScore](/packages/crisprScore)
-----------

                       Changes in version 1.15.2                        

- Conda environment now have to be built outside of crisprScore.

- The tracer argument for RuleSet3 is now fixed; the previous
  version of crisprScore was now taking into account the tracer,
  and the Chen2013 tracer was always used.

                       Changes in version 1.15.1                        

- The scores developed in Python 2 are no longer available:
  DeepSpCas9, DeepCpf1, Azimuth, and the crisprai scores.

[crisprShiny](/packages/crisprShiny)
-----------

                        Changes in version 1.7.2                        

- Removed deprecated Azimuth scores.

[CSOA](/packages/CSOA)
----

                 Changes in version 1.1.1 (2025-11-11)                  

- Changed the implementation of multiple correction testing methods and
  the
  default to Benjamini-Hochberg.

- Added an error check for gene signatures containing repeated genes

[CytoMDS](/packages/CytoMDS)
-------

                         Changes in version 1.7                         

CytoMDS 1.7.2

- added ggplotVolcano()

CytoMDS 1.7.1

- fixed unti tests with new version of ggplot2

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

                        Changes in version 1.11                         

CytoPipeline 1.11.1

- Bug correction for scale transfo pipelines with sample specific
data: implemented a work-around in readSampleFiles(): forced flowSet
sample names to full file paths as this was not done automatically
by flowCore::read.flowSet().
- estimateScaleTransforms() now has ... in parameter list, allowing
to
pass additional parameters to flowCore::estimateLogicle().

[CytoPipelineGUI](/packages/CytoPipelineGUI)
---------------

                         Changes in version 1.9                         

CytoPipelineGUI 1.9.1.

- bug correction: since version 1.7.2, it was no more possible to
display flowFrames from the scale transformation pipeline part =>
this has been corrected
- bug correction: logicle transform now uses all (t,m,w,a)
parameters.
Previously, t (maximum scale) was always set to a default value
hence providing wrong visualization => this has been corrected.

[damidBind](/packages/damidBind)
---------

                       Changes in version 0.99.14                       

- Initial Acceptance toBioconductor.


[decemedip](/packages/decemedip)
---------

                       Changes in version 0.99.8                        

- Initial Acceptance to Bioconductor

[DeconvoBuddies](/packages/DeconvoBuddies)
--------------

                        Changes in version 1.3.3                        

NEW FEATURE

- plot_gene_express() now has option to add label points with
parameter label_points.

                        Changes in version 1.3.2                        

NEW FEATURES

- findMarkers_1vAll() now only returns standardized log fold-change
values, which cuts default run time roughly in half. A new parameter
raw_logFC has been added; when TRUE, it yields the old behavior of
returning both versions of the log fold-change.
- findMarkers_1vAll() can now be parallelized with near-linear
speedup
via a new BPPARAM parameter.

                        Changes in version 1.3.1                        

NEW FEATURES

- get_mean_ratio() has a new BPPARAM parameter allowing
parallelization. Run time without parallelization was also
marginally improved.

[DeeDeeExperiment](/packages/DeeDeeExperiment)
----------------

                        Changes in version 1.2.0                        

- The dea and fea slots now use S4Vectors::SimpleList instead of base
list

- The DEA original objects are no longer stored in the dea slot. They
are now stored under metadata(dde)$singlecontrast or
metadata(dde)$multicontrast, and each contrast entry in the dea slot
contains a pointer to its original object

- Added helper functions to simplify inserting DE results from
muscat::pbDS() and MArrayLM object with multiple contrast into a dde
object

- Added a helper function to export results stored in dea and fea
slots as excel files

- Improved warning messages

- A new vignette showcasing how to apply DeeDeeExperiment to a
single-cell dataset

- Updating citation since the DeeDeeExperiment paper is out!

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

[DeMixT](/packages/DeMixT)
------

                        Changes in version 2.0.0                        

- 
  Added DeMixNB: a negative binomial deconvolution method for
  sparse count data.

[DenoIST](/packages/DenoIST)
-------

                 Changes in version 0.99.4 (2026-03-17)                 

- Initial Acceptance to Bioconductor

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.51.3                        

- Numerous partial matching fixes from Hugo Gruson.

[DESpace](/packages/DESpace)
-------

                        Changes in version 2.3.2                        

- add quasi-likelihood F tests

[DFplyr](/packages/DFplyr)
------

                        Changes in version 1.6.0                        

NEW FEATURES

- Fixed column-wise slicing
- Aligned to S4Vectors updates

SIGNIFICANT USER-VISIBLE CHANGES

- Grouped DataFrame now appears as a new class GroupedDataFrame with
its own dispatched methods

[DMRcaller](/packages/DMRcaller)
---------

                       Changes in version 1.42.1                        

- Added functions to read Oxford Nanopore methylation data from Dorado
  bam files

- Added functions to detect Partially Methylated Domains PMDs

- Added functions to detect Variably methylated domains (VMDs) using
  ONT data

- Added functions to detect Variably Methylated Regions (VMRs) using
  ONT data

- Added functions to detect co-methylation using ONT data

- Updated parallel computing to use BiocParallel

[dominatR](/packages/dominatR)
--------

                       Changes in version 0.99.4                        

- Initial Acceptance to Bioconductor

[dominoSignal](/packages/dominoSignal)
------------

                        Changes in version 1.6.0                        

New Features

- Added print() and show() methods for linkage_summary objects to
provide concise summary output.

Bug Fixes

- Fixed create_rl_map_cellphonedb() handling of partner B complex
mappings and gene assignment.
- Fixed gene_network() to avoid repeated prefixing of outgoing
cluster
names and to correctly subset outgoing signaling matrices.
- Fixed gene_network() to only include receptor and TF nodes that are
associated with a ligand if OutgoingSignalingClust is used.
- Fixed gene_network() to accumulate ligand expression across
clusters
for ligand node scaling.
- Fixed signaling_network() to assign undefined (NA) vertex sizes to
0
when scaling by signaling.
- Fixed dom_linkages() with by_cluster = TRUE and link_type =
"tf-receptor" to return clust_tf_rec.
- Fixed dom_signaling(cluster = ...) to return the selected cluster
matrix via list indexing.

Documentation

- Added figure alt text to images in vignettes for accessibility.
- Updated pkgdown and vignette links to use working URLs.
- Updated README/index documentation links and citation text to
current release metadata.


[DOSE](/packages/DOSE)
----

                        Changes in version 4.5.1                        

- mv get_organism to 'GOSemSim' (2025-12-07, Sun)
- update enrichment functions to use 'enrichit' (2025-12-05, Fri)
- mv enricher_internal and GSEA_internal to the 'enrichit' package
(2025-12-05, Fri)
- renamed to ora_gson and gsea_gson with updates
- mv helper functions to 'enrichit'

[DOtools](/packages/DOtools)
-------

                        Changes in version 1.1.8                        

All notable changes to DOtools will be documented in this file.

&#91;1.1.8&#93; - 2026-04-13

Added

- Added citation from Bioinformatics Advances since it got published.

Changed

- DO.UMAP: Is now able to plot a DensityPlot for provided genes.

&#91;1.1.2&#93; - 2026-01-13

Added

- DO.HeatmapFC: Heat map function for foldchange visualisation,
includes statistics about significance.
- DO.Barplot: New barplot function with more testing methods.

Changed

- Major adjustments to DO.VlnPlot and DO.Barplot. These functions now
take multiple testing methods and correction methods. Will be
applied also to DO.BoxPlot.
- Major bug fixes for comparability with Seurat's newest version.

Removed

- DO.BarplotWilcox: Replaced by DO.Barplot

[DOTSeq](/packages/DOTSeq)
------

                       Changes in version 0.99.8                        

- Accepted into Bioconductor 

[drawProteins](/packages/drawProteins)
------------

                        Changes in version 4.0.0                        

- Changing tests based on an updated version of ggplot2 to 3.5.0 which
  changes
  the length of the ggplot object.
  ARGUMENT

- Replace size with linewidth in function draw_chain() following a
  depreciation in ggplot.

[dreamlet](/packages/dreamlet)
--------

                        Changes in version 1.9.2                        

- fixed failed check

[DropletUtils](/packages/DropletUtils)
------------

                       Changes in version 1.32.0                        

- Updated downsampleReads() to use the new downsampling algorithm
  from scuttle. This is more efficient but slightly changes the
  results compared to the previous version.

[drugfindR](/packages/drugfindR)
---------

- Initial Acceptance to Bioconductor

[edgeR](/packages/edgeR)
-----

                 Changes in version 4.10.0 (2026-04-25)                 

- 
  catchOarfish() now uses nanoparquet::read_parquet() instead of
  arrow::read_parquet().

- 
  New function DGEListFromTximport() to assemble a DGEList object
  from the list object created by tximport::tximport.

- 
  New function DGEListFromTximeta() to assemble a DGEList object
  from the list object created by tximeta::tximeta.

- 
  glmQLFit() with `legacy=FALSE` now caps any preset dispersions
  to be <= 4. Previously, large preset dispersions could lead to
  a segmentation error.

- 
  glmQLFit() with `legacy=FALSE` no longer automatically
  estimates the NB dispersion from trended values found in the
  DGEList object, if present. The results will now be the same
  whether the DGEList contains dispersion estimates or not.

- 
  New function sampleWeights() to estimate sample weights from a
  matrix of adjusted unit deviances.

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.31.5                        

- cnetplot.compareClusterResult() now supports categorySizeBy for
category pie sizing and aligns docs with ggtangle::cnetplot()
semantics (2026-04-22, Wed)
- ridgeplot now supports stat parameter (default is 'density_ridges'
and can be changed to 'binline') (2026-04-01, Wed, #343)
- manhattan plot for enriched result (2026-03-26, Thu)
- update roxygen document to use markdown syntax (2026-03-02, Mon)
- bug fixed in xy-lab format in ssplot() (2026-03-02, Mon)
- bug fixed in formula supports in dotplot() (2026-02-26, Thu)

                       Changes in version 1.31.4                        

- fix cnetplot() S3 generic/method consistency warnings (2026-01-14,
Wed)
- fix treeplot() column selection bug when color variable equals size
variable (2026-01-14, Wed)
- fix fortify.compareClusterResult() warnings about missing imports
and global variables (2026-01-14, Wed)
- remove plyr and use dplyr in method-fortify.R (2026-01-14, Wed)
- fixed treeplot() issue where pairwise_termsim() with method="JC"
produced unnamed similarity matrix, causing "undefined column
selected" error (2025-01-14)
- fixed fortify.compareClusterResult() warning "NAs introduced by
coercion" when Cluster names are not numeric (2025-01-14)
- bug fixed in barplot() as fortify() generic in ggplot2 checks for
unused arguments in ... (2026-01-14, Wed)
- remove categorySize parameter in cnetplot() (2026-01-14, Wed)
- bug fixed in goplot() as GOSemSim uses cache (2026-01-13, Tue)
- also fix gotbl object not found issue (2026-01-13, Tue)
- re-export geneID, geneInCategory and gseaScores from 'enrichit'
(2026-01-12, Mon)
- update documentation: fix typos, grammar errors and use modern
markdown syntax (2026-01-12, Mon)
- bug fixed in update_n() if showCategory is a vector of term names
(2026-01-08, Thu)
- avoid the "condition has length > 1" error in outer() by using
Vectorize() (2026-01-08, Thu)

                       Changes in version 1.31.3                        

- use 'enrichit' package (2025-12-07, Sun)
- optimize source code (2025-12-02, Tue)
- error handling functions imported from 'yulab.utils' (2025-12-01,
Mon)

                       Changes in version 1.31.2                        

- add 'fc_threshold' parameter to cnetplot (2025-11-30, Sun, #338)
- requires 'ggtangle' v>= 0.0.9
- update all line width aes mapping from 'size' to 'linewidth'
(2025-11-30, Sun)
- add 'node_label_size' parameter for emapplot (2025-11-30, Sun)
- remove emapplot parameters, 'group', 'group_style' and
'label_group_style' (#339)
- add 'showTop' parameter to limit number of genes shown in
heatplot()
and distinguish tip point size variable for treeplot() through
internal parameter size_var (2025-11-23, Sat, #335)

                       Changes in version 1.31.1                        

- import ggfun::%<+% (2025-11-18, Tue)
- update ssplot(), treeplot() and get_wordcloud() (2025-11-15, Sat)
- change set_enrichplot_color(transform = 'identity') as default
behavior (2025-11-11, Tue)
- now it only sets the color scale without changing the transform
method
- explicitly set transform = 'log10' in dotplot
- use 'quarto' as vignette engine (2025-11-11, Tue)
- use set_enrichplot_color(transform = 'identity') in heatplot
(2025-11-11, Tue)
- use set_enrichplot_color(transform = 'identity') in cnetplot
(2025-11-05, Wed)

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.19.3 (2026-01-25)                 

- checks conventional base probabilities in long-read data

- accuracy comparison vs modkit

                 Changes in version 1.19.2 (2026-01-19)                 

- limit/clip reads to targets during BAM loading

- filtering by number of context sites

- new defaults for filtering

                 Changes in version 1.19.1 (2026-01-02)                 

- API change for thresholding in all 'generate*Report' functions

- filtered out reads are dropped and not counted as hypomethylated
  anymore

[epigraHMM](/packages/epigraHMM)
---------

                       Changes in version 1.19.1                        

- epigraHMM examples and the vignette no longer use the chromstaRData
package due to its deprecation

[epiregulon](/packages/epiregulon)
----------

                        Changes in version 2.0.1                        

Fixed bugs that caused addWeights to fail when aggregateCells = TRUE.
When aggregateCells is set to TRUE, cells are aggregated by
calculating
the mean of the original cell features rather than their sums. Fixed
a
warning about the class of the cellNum argument passed to
calculateP2G,
which was raised incorrectly.

[epiRomics](/packages/epiRomics)
---------

                       Changes in version 0.99.5                        

Changes

- Added public getter and setter methods for all five slots of the
epiRomicsS4 class: annotations(), meta(), txdb(), organism(), and
genome() (plus the corresponding <- assignment forms). Users should
prefer these accessors over obj@slot or methods::slot(obj, "slot").
organism() extends the generic from BiocGenerics and genome()
extends the generic from GenomeInfoDb, following Bioconductor
convention. Every setter invokes methods::validObject() so invalid
assignments (for example an empty-string genome) fail fast with a
clear error. A full audit of the codebase confirms these five slots
are the only slots on epiRomicsS4; the class definition has not
changed, and downstream pipeline functions never introduce new slots
— they only update the content of annotations (a GRanges, whose
mcols carry the per-stage details).
- Added ?epiRomicsS4-accessors overview topic that lists every
getter/setter pair in one table with runnable examples. The
?epiRomicsS4-class page cross-references each individual accessor
via @seealso and documents each slot with an explicit "Access via
..." note pointing at the matching getter.
- Refactored both vignettes to demonstrate the new accessor API.
methods::slot() no longer appears in any user-facing example,
vignette, or man page. getting-started-with-epiRomics.Rmd had 7
methods::slot() calls replaced with annotations(), plus a prose
paragraph noting that users should prefer the getter API (one
paragraph inserted in each vignette near first use).
articles/full-analysis-with-epiRomics.Rmd had 22 methods::slot()
calls replaced.
- Every roxygen @examples block that previously used
methods::slot(obj, "...") now uses the matching getter (touched
files: R/synthetic-data.R, R/enhanceosomes.R, R/enhancers.R,
R/regions_of_interest.R).
- Swapped the four runtime methods::slot(database, "...") calls
inside
make_example_enhanceosome() in R/synthetic-data.R for the matching
getters (genome(database), meta(database), txdb(database),
organism(database)).
- Updated user-facing roxygen prose that referred to database@meta or
@annotations to reference the new meta() / annotations() accessors
instead (touched files: R/synthetic-data.R, R/chromatin_states.R,
R/enhancers.R, R/benchmark_enhancer_predictor.R,
R/regions_of_interest.R).
- Added tests/testthat/test-accessors.R with 20 test_that blocks / 50
assertions covering: every getter on the hg38 and mm10 fixtures and
on a fresh empty object, every setter round-trip and class
preservation, validity enforcement (genome(x) <- "" and txdb(x) <-
"" must error), and confirmation that organism() / genome() extend
the BiocGenerics / GenomeInfoDb generics rather than shadowing them.
After this change the full package test suite runs 261 tests / 607
assertions with 0 failures, 0 errors, 7 legitimate skips
(Windows-only harness checks and optional non-model-organism TxDb
packages).
- Added a one-time session tip to the rename-warning helper
(.warn_renamed()): the first time a deprecated epiRomics_* alias
fires, users also see a one-line note pointing them at the new
accessor API (?epiRomicsS4-accessors).
- Fixed the test-deprecated-aliases.R harness, which was still
written
against the pre-0.99.4 .Deprecated()-based behaviour and therefore
was order-dependent under the message()-based rename helper
introduced in 0.99.4. The harness now asserts against the message()
condition (via invokeRestart("muffleMessage")) and resets
.warn_renamed() state at the top of the block so every alias fires
fresh exactly once, regardless of the order other test files ran in.
No package code was changed by this test fix.
- Quality gates for 0.99.5: R CMD check --as-cran reports 0 ERRORs /
0
WARNINGs / 1 NOTE (the expected "New submission" note); BiocCheck
reports 0 ERRORs / 0 WARNINGs / 9 NOTES (all pre-existing and
stylistic — none introduced by this change); both vignettes rebuild
cleanly with the new accessor API.

                       Changes in version 0.99.4                        

Changes

- Replaced .Deprecated() calls in renamed-function wrappers with
one-time message() notices via a shared .warn_renamed() helper.
Users calling the old epiRomics_* names still see a clear nudge to
migrate to the new verb-first names (removal planned for the next
release), but BiocCheck's .Deprecated / .Defunct usage warning no
longer fires.

                       Changes in version 0.99.3                        

Highlights

- On-demand example data: cache_data() downloads example datasets
(~1.3 GB) via BiocFileCache, keeping the package tarball under 5 MB.
has_cache() and get_cache_path() check and locate cached data.

- Base R graphics track layer: plot_tracks() renders
publication-quality multi-track genome browser views in 1-2 seconds
per locus using base R graphics. ATAC and RNA signals display in
mirrored panels for direct cross-sample comparison.

- Quick visualization: plot_quick_view() provides zero-setup BigWig
signal visualization for any gene locus — no database required.

- Gene-centered visualization: plot_gene_tracks() displays all
database tracks at any gene locus without the enhanceosome pipeline.

- Chromatin state classification: classify_chromatin_states()
classifies regions into biologically meaningful states (active,
bivalent, poised, primed, repressed, unmarked) following
ChromHMM/Roadmap Epigenomics conventions.

- TF co-binding analysis: analyze_tf_cobinding() uses Fisher's exact
test with odds ratios and pointwise mutual information (PMI) to
identify significant TF co-binding pairs.

API

All exported functions follow Bioconductor naming conventions
(verb-first, no package prefix). Legacy names remain exported as
.Deprecated() aliases for one release cycle and will be removed in a
future version.

Function names:

- epiRomics_build_dB() → build_database()
- epiRomics_quick_view() → plot_quick_view()
- epiRomics_track_layer() → plot_tracks()
- epiRomics_track_layer_fast() → plot_tracks_fast()
- epiRomics_track_layer_gene() → plot_gene_tracks()
- epiRomics_putative_enhancers() → find_putative_enhancers()
- epiRomics_enhancers_co_marks() → find_enhancers_by_comarks()
- epiRomics_enhanceosome() → find_enhanceosomes()
- epiRomics_enhancers_filter() → filter_enhancers()
- epiRomics_filter_accessible() → filter_accessible_regions()
- epiRomics_chromatin_states() → classify_chromatin_states()
- epiRomics_chromatin_states_categories() →
chromatin_state_categories()
- epiRomics_annotate_putative() → annotate_enhancers()
- epiRomics_tf_cobinding() → analyze_tf_cobinding()
- epiRomics_tf_overlap() → analyze_tf_overlap()
- epiRomics_enhancer_predictor_to_ref() →
benchmark_enhancer_predictor()
- epiRomics_regions_of_interest() → get_regions_of_interest()
- epiRomics_cache_data() → cache_data()
- epiRomics_cache_path() → get_cache_path()
- epiRomics_has_cache() → has_cache()

Parameter names (no deprecation shim — update calls directly):

- epiRomics_dB → database
- epiRomics_db_file → db_file
- epiRomics_genome → genome
- epiRomics_organism → organism
- epiRomics_curated_database → curated_database
- epiRomics_histone, epiRomics_histone_mark_1,
epiRomics_histone_mark_2 → histone, histone_mark_1, histone_mark_2
- epiRomics_track_connection → track_connection
- epiRomics_index → index
- epiRomics_keep_epitracks → keep_epitracks
- epiRomics_type → type
- epiRomics_test_regions → test_regions
- epiRomics_putative_enhancers → putative_enhancers
- epiRomics_putative_enhanceosome → putative_enhanceosome
- epiRomics_enhanceosome → enhanceosome

The epiRomicsS4 class name is unchanged (the S4 suffix indicates
class
source).

- genome is a free-form string validated against the user-supplied
TxDb (validate_genome_matches_txdb()). No hardcoded whitelist; any
organism with a TxDb.* and org.*.db package is supported.

Dependencies

- TxDb.Hsapiens.UCSC.hg38.knownGene,
TxDb.Mmusculus.UCSC.mm10.knownGene, org.Hs.eg.db, org.Mm.eg.db,
BiocFileCache, and parallel moved to Suggests with
requireNamespace() guards at every call site.
- Full @importFrom audit across all R files; NAMESPACE regenerated.

Examples and vignettes

- Ships a ~0.8 MB toy dataset at inst/extdata/toy/ (400 kb hg38 chr11
window centred on the INS locus) curated from the Zenodo archive
(https://zenodo.org/records/19189987). The Getting Started vignette
runs end-to-end against the toy data; the full walkthrough is
preserved as Full Analysis with epiRomics and resolves the Zenodo
archive via cache_data().
- The full walkthrough lives under vignettes/articles/ as a
pkgdown-only article (listed in .Rbuildignore). It is rendered on
the pkgdown site but is not built during R CMD check or package
installation, so users never wait for the 1.3 GB Zenodo download at
install time. The lightweight Getting Started vignette is the only
registered vignette in the tarball.
- Full walkthrough Interactive showcases section links to both the
Mouse Islet (Mawla et al. 2023) and Human Islet (Mawla &
Huising 2021) Shiny browsers.
- Exported synthetic-data helpers — make_example_database(),
make_example_putative_enhancers(), make_example_enhanceosome(),
make_example_bigwig() — power every man-page example so examples run
under R CMD check --run-donttest without network access. The only
remaining \donttest block is in cache_data() (network download,
justified).
- Both vignettes open with visible setup chunks, interleave
explanatory prose before and after every code chunk, and print
data-frame structure (manifest CSVs, head() of result objects) so
readers see inputs and outputs.
- Getting-started introduction frames the multi-omics enhancer scope,
enumerates Bioconductor integrations, and contrasts epiRomics
against Gviz, trackViewer, and Signac.
- vignettes/references.bib holds package citations; vignettes use
&#91;@mawla2023; @mawla2021&#93; inline references.

Documentation and installation

- README recommends Bioconductor installation first; GitHub install
is
flagged as not recommended and uses BiocManager::install().
- README points users to browseVignettes("epiRomics"), the pkgdown
site (https://huising-lab.github.io/epiRomics/), and the
Bioconductor landing page.
- pkgdown site restored: _pkgdown.yml and a dedicated
.github/workflows/pkgdown.yaml rebuild the site on every push to
main and on every release.
- Top-level doc/ directory removed from the source tree.
- Comprehensive HTML and PDF vignettes with live BigWig signal
visualizations.
- .onAttach() startup message with version, cache status, and
citation
info.

Speed

- Parallel overlap counting via parallel::mclapply() when available.
- Pre-split annotations by type for O(1) lookup in
find_enhanceosomes().
- Vectorized chromatin state classification and signal aggregation.
- BigWig RDS caching for memory-efficient large-scale analyses.

CI/CD

- Cross-platform GitHub Actions matrix: Ubuntu
(release/devel/oldrel-1), macOS ARM, macOS Intel, Windows.
- BiocCheck validation on every push/PR.
- Test coverage with full integration tests via cached example data.

Housekeeping

- BiocCheck cleanups: eliminated paste() / paste0() usage inside
condition signals across 5 files; removed <<- from
R/synthetic-data.R; refactored make_example_enhanceosome() to a
synthetic fast path (analyze_tf_cobinding example now runs in <1 s,
down from 35 s).
- Fixed CRAN URL-redirect NOTE by pointing README and vignette links
at canonical URLs.
- Aligned README R-version minimum to DESCRIPTION Depends (R ≥
4.5.0).
- Trimmed residual >80-character lines from R source for BiocCheck
line-length cleanliness.
- Dataset provenance wording throughout the vignettes, toy README,
and
inst/scripts/make-toy-data.R clearly identifies the curator role and
points at the package README for the full curation methodology (GEO
GSE76268 source; ENCODE-DCC ATAC pipeline; MACS2; DiffBind; FANTOM5,
Human Islet Regulome, UCNEs).

[epiSeeker](/packages/epiSeeker)
---------

                       Changes in version 0.99.12                       

- Initial Acceptance to Bioconductor

[escape](/packages/escape)
------

                        Changes in version 2.7.3                        

BUG FIXES

- Robust Seurat v5/v3 assay handling: Fixed compatibility with latest
SeuratObject where the slot argument to GetAssayData() is now
defunct. All internal accessors (.cntEval(), .pull.Enrich()) now use
the layer API for SeuratObject >= 5.0.0 and fall back to slot only
for older versions.
- Updated .adding.Enrich() to use the SeuratObject package version
(not sc@version) when selecting CreateAssay5Object vs
CreateAssayObject, preventing mismatches on objects created across
Seurat versions.
- Fixed performPCA() error ("Enrichment matrix must be numeric") when
enrichment data is returned as a sparse matrix from the Seurat v5
layer API.

ENHANCEMENTS

- CI/CD: Added devel branch to R-CMD-check and test-coverage workflow
triggers
- Automated releases: Added GitHub Actions workflow to create
releases
from version tags (push v* tags to trigger)

                        Changes in version 2.6.2                        

NEW FEATURES

- Added .themeEscape() internal theme function for consistent
visualization styling across all plotting functions

ENHANCEMENTS

- Seurat v5 compatibility: Updated .cntEval() to detect SeuratObject
version and use layer argument instead of deprecated slot argument
for SeuratObject >= 5.0.0
- Consistent theming: Applied unified theme styling across all
visualization functions (ridgeEnrichment(), splitEnrichment(),
geyserEnrichment(), heatmapEnrichment(), scatterEnrichment(),
pcaEnrichment(), densityEnrichment(), gseaEnrichment(),
enrichItPlot())
- Improved densityEnrichment(): Added plot title showing gene set
name, alphanumeric sorting of group labels, and improved rug segment
styling

DOCUMENTATION

- Reformatted roxygen2 documentation across all exported functions
for
consistency
- Standardized use of \code{}, \itemize{}, \enumerate{}, \strong{},
and \emph{} tags
- Replaced Unicode characters with ASCII equivalents for better
portability

                 Changes in version 2.6.1 (2025-10-31)                  

Bioconductor Release 3.22

BUG FIXES

- Fixed densityEnrichment() interaction with GSVA package through
compute.gene.cdf — corrected boolean argument for internal CDF
computation

[EventPointer](/packages/EventPointer)
------------

                         Changes in version 4.0                         

- Update version of EventPoint to work with BAM files

- Use of bootstrap statistics in novel events

- Update of vignette

[EWCE](/packages/EWCE)
----

                       Changes in version 1.19.1                        

Bug fixes

- bootstrap_enrichment_test / check_species
- Validate species arguments before background generation.
- Fail fast when fewer than the minimum number of input hits are
supplied, so invalid inputs return errors instead of upstream
warnings.

[ExpoRiskR](/packages/ExpoRiskR)
---------

                       Changes in version 0.99.3                        

- Initial Acceptance to Bioconductor

[fastreeR](/packages/fastreeR)
--------

                        Changes in version 2.7.1                        

- Windowed / streaming VCF distance and tree output. Emit one
distance
matrix (or Newick tree) per genomic window of N base pairs or per N
consecutive variants for VCF2DIST and VCF2TREE.

                 Changes in version 2.5.0 (2026-02-02)                  

- Added embedding-based distance calculation for VCF files
(backend 2.5.0).

                 Changes in version 2.3.0 (2025-12-01)                  

- Implement reading compressed VCF input (gz, bzip2, xz).

                 Changes in version 2.2.0 (2025-11-01)                  

- Major update: enhanced streaming, faster bootstrapping, and backend
v2.2.0 integration.

[fenr](/packages/fenr)
----

                        Changes in version 1.9.2                        

- Corrected response to the error message from the responsive
WikiPathways server, so it does not through error when on_error =
"warn".

                        Changes in version 1.8.1                        

- Top gene ontology URL does not seem to be accessible directly
anymore; changed the URL in fetch_go_genes_go() accordingly.

[fgsea](/packages/fgsea)
-----

                       Changes in version 1.37.1                        

- Hash-based resolution of tied gene sets in the multilevel algorithm
  (by Nikita Golikov), properly fixes #151

- P-avlues for ES close to 1 are more accurate now

- Internally, the gene weights are now converted to integers (with
  scaling), which might introduce minor result discrepancies in rare
  cases

[FlowSOM](/packages/FlowSOM)
-------

                       Changes in version 2.19.4                        

- Removed #include &#60;PrtUtil.h&#62; from the C script, as this is not
  supported anymore, and I don't think we actually still used it.

                       Changes in version 2.19.3                        

- Attempt to fix git issue with double file name

                       Changes in version 2.19.2                        

- In NewData, only give warning messages when paramter in fsom object
  is TRUE

[flowSpecs](/packages/flowSpecs)
---------

                 Changes in version 1.25.2 (2026-03-17)                 

- +Bug fixes for version created yesterday

                 Changes in version 1.25.1 (2026-03-16)                 

- +Updates for S8 data, but with more general applicability:
  +specMatCalc now allows for a custom set of exclusion parameters
  +specMatCalc now allows for exclusion of the scatter gating step

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

[fourSynergy](/packages/fourSynergy)
-----------

                       Changes in version 0.99.12                       

- Initial Acceptance to Bioconductor

[fRagmentomics](/packages/fRagmentomics)
-------------

                Changes in version 0.99.12 (2026-03-31)                 

Pre-release version on Bioconductor devel.


[fraq](/packages/fraq)
----

                 Changes in version 0.99.2 (2026-02-12)                 

- Initial Acceptance to Bioconductor

[FRASER](/packages/FRASER)
------

                        Changes in version 2.6.1                        

- Internal: use moved makeTxDbFromGFF function from txdbmaker instead
  from GenomicFeatures.

[gcatest](/packages/gcatest)
-------

                 Changes in version 2.11.1 (2026-01-28)                 

- Updated vignette from old Sweave to modern R markdown.

[gDNAx](/packages/gDNAx)
-----

                 Changes in version 1.10.0 (2026-04-27)                 

USER VISIBLE CHANGES

- Changed default useRMSK=TRUE to useRMSK=FALSE in gDNAdx().

BUG FIXES

- Fix in the calculations of the strandedness() function, and its
  corresponding unit tests.

[gDR](/packages/gDR)
---

               Changes in version 2025-10-30 (2025-10-30)               

- synchronize Bioconductor and GitHub versioning

               Changes in version 2025-06-11 (2025-06-11)               

- update installation instructions


[gDRcore](/packages/gDRcore)
-------

               Changes in version 2026-04-27 (2026-04-27)               

- fix C++ compilation error on R-devel (STRING_PTR → STRING_PTR_RO)

               Changes in version 2026-04-13 (2026-04-13)               

- migrate from qs to qs2 package (qs::qread → qs2::qs_read, qs::qsave
→ qs2::qs_save)
- update checkpoint and generated data file extensions from .qs to
.qs2

               Changes in version 2026-02-18 (2026-02-18)               

- fix handling of numeric identifiers in annotation functions

               Changes in version 2026-02-16 (2026-02-16)               

- fix stack imbalance warnings during byte-compilation

               Changes in version 2026-02-02 (2026-02-02)               

- remove duplicated code from map_references in map_untreated

               Changes in version 2025-11-27 (2025-11-27)               

- add support for Incucyte data

               Changes in version 2025-11-17 (2025-11-17)               

- remove empty clids from day0 data during merging raw data


[gDRimport](/packages/gDRimport)
---------

               Changes in version 2026-04-23 (2026-04-23)               

- add support for generation of day0 template for tdd files

               Changes in version 2026-04-13 (2026-04-13)               

- migrate from qs to qs2 package (qs::qread → qs2::qs_read, qs::qsave
→ qs2::qs_save)
- update test data and reference file extensions from .qs to .qs2

               Changes in version 2026-04-07 (2026-04-07)               

- add support for excel format of new envision data

               Changes in version 2026-02-25 (2026-02-25)               

- add support for multi-plate new envision format

               Changes in version 2026-02-16 (2026-02-16)               

- fix stack imbalance warnings during byte-compilation

               Changes in version 2025-12-27 (2025-12-27)               

- extend support for csv EnVision format

               Changes in version 2025-12-17 (2025-12-17)               

- add support broad_id in PRISM data
- add suport for combo data in private prism data

               Changes in version 2025-12-10 (2025-12-10)               

- make barcode detection more robust (Incucyte)

               Changes in version 2025-11-26 (2025-11-26)               

- add support for Incucyte data

               Changes in version 2025-11-17 (2025-11-17)               

- add support for new EnVision format



[gDRutils](/packages/gDRutils)
--------

               Changes in version 2026-04-18 (2026-04-18)               

- migrate from qs to qs2 package (qs::qread → qs2::qs_read, qs::qsave
→ qs2::qs_save)
- update batch file extension from .qs to .qs2

               Changes in version 2026-04-14 (2026-04-14)               

- add support for metadata in merge_MAE

               Changes in version 2026-03-23 (2026-03-23)               

- standardize_MAE standardizes also internal identifiers

               Changes in version 2026-02-12 (2026-02-12)               

- remove_drug_batch supports atomic vectors as input

               Changes in version 2025-12-02 (2025-12-02)               

- convert_se_assay_to_dt supports merging additional variables

               Changes in version 2025-11-27 (2025-11-27)               

- add support for the time-course experiment

               Changes in version 2025-11-03 (2025-11-03)               

- update merge_SE function to merge drugs with different batches
together


[gemma.R](/packages/gemma.R)
-------

                        Changes in version 4.0.0                        

- get_annotation_children and get_annotation_parents are added.
get_child_terms is deprecated in favor of get_annotation_children

[GenomicPlot](/packages/GenomicPlot)
-----------

                        Changes in version 1.9.2                        

- Add chromInfo parameter support for non model organisms (custom
organisms)
- Bump up version number to trigger Bioconductor build

                        Changes in version 1.9.1                        

- Fix bugs caused by plyranges update
- Bump up version number to trigger Bioconductor build

[geomeTriD](/packages/geomeTriD)
---------

                        Changes in version 1.5.1                        

IMPROVEMENT

- add citation.

[GEOquery](/packages/GEOquery)
--------

                       Changes in version 2.99.1                        

Bug Fixes

- Fixed error when parsing GSE matrix files with malformed or empty
lines between sample metadata (e.g., GSE425). Sample lines are now
extracted directly using pattern matching to avoid issues with
irregular file formatting.

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

[ggtree](/packages/ggtree)
------

                        Changes in version 4.1.2                        

- make gheatmap() robust to missing tree labels in data row names by
warning and auto-filling missing rows with NA instead of stopping
(2026-03-23, Mon)
- improve gheatmap() input validation and custom column label
remapping (2026-03-09, Mon)
- add native paired-tree/tanglegram support via ggdoubletree(),
fortify.tanglegram(), fortify.cophylo(), and geom_tanglelink(),
including support for raw trees, ggtree objects, deterministic
crossing minimization, and optimize_side = "right"|"left"|"both"
regression tests (2026-03-09, Mon)
- add a deterministic tree_and_leaf layout inspired by the
TreeAndLeaf
package, with native ggtree() integration and regression tests
(2026-03-09, Mon)
- add a height parameter to collapse() supporting NULL,
ggplot2::rel(), and absolute numeric heights for collapsed
triangles, while restoring the original layout correctly on expand()
(2026-03-09, Mon, #409)
- optimize the daylight layout by caching subtree topology and
removing repeated tidyverse reshaping while preserving the existing
geometry; add regression tests for cached and uncached paths
(2026-03-09, Mon)
- fix custom layout handling in ggtree() when layout functions return
an xy component (2026-03-09, Mon)
- make geom_hilight() and geom_cladelabel() more compatible with
recent ggplot2 changes, and strengthen regression tests (2026-03-09,
Mon)
- modernize roxygen documentation using markdown syntax (2026-01-12,
Mon)
- add 'data' argument to geom_rootedge() (2026-01-10, Sat, #696,
#697)
- more options in ggtree_set_interactive() (2025-12-19, Fri)
- update default_aes to work with from_theme
- geom_highlight and geom_taxalink (2025-12-18, Thu, #694)
- fix issue of geom_taxalink by adding 'outward' argument
(2025-12-18,
Thu, #692)
- print() method for 'ggtree' object (2025-12-17, Wed)
- ggtree_set_interactive() and ggtree_unset_interactive() to set and
unset interactive mode (2025-12-17, Wed)
- if interactive mode is set and there are no interactive
attributes, it will throw a warning message, which can be
disabled if ggtree_set_interactive(check=FALSE)

                        Changes in version 4.1.1                        

- bug fixed in geom_striplab (2025-10-30, Thu, #677)

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                       Changes in version 1.21.0                        

- Bioconductor 3.22 released, and Bioconductor 3.23 (devel) bump.
(2025-10-31, Fri)

[GloScope](/packages/GloScope)
--------

                 Changes in version 2.1.1 (2025-09-27)                  

- Fix so handles if no column names which was creating error with
  latest version of `mclust`

[glycoTraitR](/packages/glycoTraitR)
-----------

                 Changes in version 0.99.0 (2025-12-04)                 

Initial Bioconductor submission

[GNOSIS](/packages/GNOSIS)
------

                 Changes in version 1.99.0 (2023-09-04)                 

- Made the following significant changes
  o added functionality to select and upload cBioPortal study
  o deprecated ability to save R script with executed code

- Submitted to Bioconductor

[GOfan](/packages/GOfan)
-----

                       Changes in version 0.99.2                        

- Initial Acceptance to Bioconductor

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.37.2                        

- refactor .initial() to use a controlled environment for robust data
loading (2026-01-15, Thu)
- enhance cache reading robustness in computeIC() and
wangMethod_internal() (2026-01-15, Thu)
- update load_onto() documentation (2026-01-12, Mon)
- enhance godata() input validation (2026-01-12, Mon)
- refactor internal ontology type checking (2026-01-12, Mon)
- fix cache integration in .initial and refactor consumers to use
get_cache_element (2026-01-12, Mon)
- fix load_onto() using cache incorrectly (2026-01-12, Mon)
- use cache mechanisms defined in 'yulab.utils' for caching
(2026-01-12, Mon)
- update URLs, fix typos and refine roxygen documentation
(2026-01-12,
Mon)

                       Changes in version 2.37.1                        

- use cache mechanisms defined in 'yulab.utils' for caching
(2025-12-08, Mon)
- update roxygen to use 'markdown' syntax and unify duplicated docs
into templates (in 'man-roxygen' folder) (2025-12-07, Sun)
- move get_organism() from 'DOSE' (2025-12-07, Sun)
- move load_OrgDb() to 'yulab.utils' (2025-12-05, Fri)

[GraphExperiment](/packages/GraphExperiment)
---------------

                       Changes in version 0.99.0                        

NEW FEATURES

- Initial Acceptance to Bioconductor

[graphite](/packages/graphite)
--------

                 Changes in version 1.57.1 (2026-04-26)                 

- Updated all pathway data.

[GSABenchmark](/packages/GSABenchmark)
------------

                 Changes in version 0.99.0 (2025-11-05)                 

- Submitted to Bioconductor

[GSVA](/packages/GSVA)
----

                         Changes in version 2.6                         

USER VISIBLE CHANGES

- Performance improvement in the calculation of GSVA scores from ranks.

- Replaced trend=gssizes by trend=sqrt(gssizes) in the vignette that
  illustrates the limma-trend pipeline with GSVA scores.

- Added support for ssGSEA, PLAGE and Z-score to work with large
  datasets stored using a DelayedArray backend, such as HDF5.

- Moving from scuttle/scran used in the single-cell vignette to
  scrapper/bluster to fix deprecation warnings.

- *Param() functions report more informative errors when the input
  expression data or gene sets do not correspond to one of the expected
  object classes.

- Added new functions saveHDF5GSVAranks() and loadHDF5GSVAranks() to
  first save on disk the output of the gsvaRanks() function, and
  second, load that output from disk to callthe gsvaScores() function
  with it.

- Improved unit test coverage to > 85%.

BUG FIXES

- Fix on the handing of missing data with the GSVA method. Bug reported
  in issue #257.

- Fix on the propagation of NA values with the GSVA method.

- Fix in detecting rows with constant (nonzero) values.

- Fix in kcdf="Gaussian" when input is SparseArray.

- Fix in fetching rows from SVT_SparseArray matrices.

- Fix on handling NaN values, which were causing a core dump, they are
  now treated as NA values. Bug reported in issue #258.

- Removed non-API calls from C code.

- Code reorganization and multiple improvements and vulnerability
  fixes.

[gVenn](/packages/gVenn)
-----

                        Changes in version 1.1.1                        

New features

- Add bg parameter to saveViz() for controlling plot background
color,
including transparent backgrounds. Users can now save plots with bg
= "transparent" for use in presentations or publications requiring
transparent backgrounds.
- Add hex sticker logo created using the hexSticker R package.
- Update graphical abstract highlighting gVenn's overlap
visualization
and extraction capabilities

Documentation

- Add example to vignette demonstrating transparent background export
using bg = "transparent" parameter in saveViz()

[hammers](/packages/hammers)
-------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[Harman](/packages/Harman)
------

                       Changes in version 1.39.0                        

- Fix Linux build failure under C++20 by correcting the
  CFactorials constructor declaration in src/CMapSelectKFromN.h.

- Updated Makevars to stop producing deprecation warnings.

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.48.0                        

- remove "System requirements specified C++11" in DESCRIPTION

- New optimized C codes for AArch64 (ARM 64-bit) architecture

- improve Makevars.win: replace `-ffixed-xmm16..31` with
  `-fno-asynchronous-unwind-tables` for AVX-512 compilation on Windows,
  allowing the GCC compiler to use all 32 XMM/ZMM registers for better
  performance


[hipathia](/packages/hipathia)
--------

                 Changes in version 3.11.2 (2026-03-25)                 

- Fist version using Zenodo as repository for necessary files

[HiSpaR](/packages/HiSpaR)
------

                 Changes in version 0.99.0 (2026-01-16)                 

- Initial Bioconductor submission

- Implements hierarchical Bayesian inference of 3D chromatin structures

- Supports MCMC sampling with cluster-based hierarchical approach

- Includes example Hi-C contact matrix dataset (su1_contact_mat)

- Provides efficient C++ implementation via Rcpp and RcppArmadillo

- OpenMP support for parallel computation

[HistoImagePlot](/packages/HistoImagePlot)
--------------

                        Changes in version 1.0.0                        

New features

- Initial Acceptance to Bioconductor 

[HPAanalyze](/packages/HPAanalyze)
----------

                        Changes in version 1.29                         

- Changes in version 1.29.1
  + Update built-in data and internal gene name lookup table to v25.0

[Ibex](/packages/Ibex)
----

                        Changes in version 1.1.1                        

- Aligned version with Bioconductor release
- Switched default branch to devel for Bioconductor compatibility
- Updated CI workflows to target devel branch
- Converted NEWS to NEWS.md format
- Added automated GitHub Release workflow via tags
- Ibex_matrix() now accepts character vectors of amino acid sequences
directly
- Removed rlang from Imports, added lifecycle
- As per basilisk documentation:
- Add .BBSoptions with UnsupportedPlatforms: win32
- Add configure and configure.win scripts
- Add Docker infrastructure with Dockerfile and
.devcontainer/devcontainer.json
- Improved testthat compatibility across platforms
- Improve adherence to verbosity arguments

[imageFeatureTCGA](/packages/imageFeatureTCGA)
----------------

                        Changes in version 1.0.0                        

New features

- Added a NEWS.md file to track changes to the package.

[imageTCGAutils](/packages/imageTCGAutils)
--------------

                        Changes in version 1.0.0                        

New features

- Initial Bioconductor release!

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.17.2 (2026-04-24)                 

- Replaced deprecated aggregateAcrossCells function from scuttle
  package with scrapper implementation

                 Changes in version 1.17.1 (2026-03-30)                 

- update function page/vignette

[immApex](/packages/immApex)
-------

                        Changes in version 1.5.4                        

BUG FIXES

- Fixed vignette build failure when IMGT network data is unavailable
(e.g., on Bioconductor build servers)
- Fixed operator precedence bug in inferCDR() reference validation
check

                        Changes in version 1.5.3                        

UNDERLYING CHANGES

- getIMGT() now uses immReferent as a backend for downloading and
caching IMGT reference sequences
- Removed frame, max.retries, and verbose parameters from getIMGT()
- Added refresh parameter to getIMGT() for controlling cache behavior
- Removed hash, httr, and rvest package dependencies (now handled by
immReferent)
- Replaced .parseSpecies() with .mapSpecies() for immReferent species
name compatibility

                        Changes in version 1.5.2                        

- improved error message for the propertyEncoder() function when the
suggested Peptides package is not installed.

                        Changes in version 1.5.0                        

- Bioconductor devel branch update to match versions


[immLynx](/packages/immLynx)
-------

                       Changes in version 0.99.4                        

- Initial Bioconductor submission

[immReferent](/packages/immReferent)
-----------

                       Changes in version 0.99.7                        

New Features

- Initial Acceptance to Bioconductor

[InPAS](/packages/InPAS)
-----

                       Changes in version 2.18.1                        

- Remove plyranges::filter.

[IntEREst](/packages/IntEREst)
--------

                       Changes in version 1.36.0                        

DEPRECATED AND DEFUNCT

- miRBaseBuild parameter is removed from referencePrepare
  function as it was defunct in txdbmaker package.

[iscream](/packages/iscream)
-------

                        Changes in version 1.1.8                        

INTERNAL

- Remove an unused debug string for printing the locus ID from
src/query_all.cpp

DOCUMENTATION

- Use blockquote the "Getting Started" vignette to call out the note
about chromosome formatting in input regions - chr1 vs 1.

                        Changes in version 1.1.7                        

BUG FIXES

- Fix htslib_version() to return the version as a string instead of
simply printing it to allow programmatic checks for libdeflate
presence.

- iscream now informs users whether they have libdeflate-enabled
htslib on package load.

                        Changes in version 1.1.6                        

BUG FIXES

- Fix a check for missing input files by correcting checking the
filepath vector length

INTERNAL

- Use inherits(x, "data.frame") instead of "data.frame" %in%
class(x)"

                        Changes in version 1.1.5                        

BREAKING CHANGES to summarize_regions() and summarize_meth_regions():

- The summary output now has position as chromosome, start, and end
columns instead of just a feature column with the position string.
This allows direct conversion to GRanges objects without having to
use rownames. It is still possible to use chromosome names as the
input query regions ("chr1"), the start and end in those cases will
be NA.

- The feature column is now populated only if the input regions are
named, if a vector, or the feature_col argument is set to a column
for data.frame/GRanges inputs.

- The count column from summarize_meth_regions() is now named count
instead of cpg_count.

- These functions now return a data.table instead of a data.frame to
allow for faster in-place modifications and consistency with tabix()
output. Use data.table::setDF() to convert to a data.frame in-place.

BUG FIXES

- Rownames are no longer set for the summarize_regions() or
summarize_meth_regions() output using the input region strings as
these were repeating and not unique. It was only possible because
they were set from the C++ scope instead of R.

ENHANCEMENTS

- summarize_regions() and summarize_meth_regions() now support the
following new functions, inspired by bedtools map:

- Population standard deviation (pstddev)
- First element (first)
- Last element (last)
- Mode (mode)
- Anti-mode (antimode)
- Absolute min (absmin)
- Absolute max (absmax)
- Count of unique values (count_distinct)

- get_granges_string() can now extract names from its mcols using a
feature_col argument as in get_df_string() - names were previously
pulled only from names(). This means GRanges inputs to
summarize_regions() can have the feature_col in its mcols rather
than just as its names.

DOCUMENTATION

- Update Zenodo URLs in vignettes to the updated record:
https://zenodo.org/records/18089082. The genes.bed file used in
vignette("performance") now contains the gene names that the
vignette references.

INTERNAL

- Refactored summarize_regions to collect similar validation and
regions parsing functions

- tabix() writes the input regions of interest to disk only once
instead of doing it for every file

                        Changes in version 1.1.4                        

BUG FIX

- Reduced the minimum R version from 4.5 to 4.4

                        Changes in version 1.1.3                        

DOCUMENTATION

- Document how to make use of strand information with tabix()
- Improve wording of tabix()'s' aligner and col.names documentation

                        Changes in version 1.1.2                        

- Fix typo in error message from thread count checks

                        Changes in version 1.1.1                        

BUG FIXES

- Fix configure to correctly check htslib version when rhtslib is not
used and not issue spurious warnings

[iSEE](/packages/iSEE)
----

                       Changes in version 2.23.1                        

- Use linewidth to set width of lines in plots, require ggplot2 >=
3.4
- Allow sampling resolution and guide sizes to be non-integer values
- Remove Windows line endings when validating inputs in the ace
editor
- Add line to load ggplot2 in the exported code
- Fix bug in the reporting of collapsible box status. Addresses #701

[IsoformSwitchAnalyzeR](/packages/IsoformSwitchAnalyzeR)
---------------------

                       Changes in version 2.11.1                        

- Update type: Minor.

- Updated deprecated functions and minor modifications

- Updated version to sync with Bioconductor

- Added Elena as contributed to package

[isomiRs](/packages/isomiRs)
-------

                       Changes in version 1.39.1                        

FIX

- Fix ggplot2 compatibility issues with labs() function - removed
  deprecated list() syntax

- Fix ggplot2 deprecation warning by replacing aes_string() with modern
  aes() and tidy evaluation

[jvecfor](/packages/jvecfor)
-------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor.
- Initial implementation of fastFindKNN(), drop-in replacement for
BiocNeighbors::findKNN using the jvecfor Java backend (HNSW-DiskANN
and VP-tree).
- Added fastMakeSNNGraph() and fastMakeKNNGraph() convenience
wrappers
delegating graph construction to the bluster package.
- Supports euclidean, cosine, and dot_product distance metrics.
- Multi-threaded index build and search via Java ForkJoinPool.
- Optional Product Quantization (PQ) for approximately 4-8x search
speedup.

[KEGGREST](/packages/KEGGREST)
--------

                       Changes in version 1.52.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The `mark.pathway.by.objects` function has been removed.

[LACHESIS](/packages/LACHESIS)
--------

                       Changes in version 0.99.5                        

NEW FEATURES

o Initial Acceptance to Bioconductor

[lcmsPlot](/packages/lcmsPlot)
--------

                       Changes in version 0.99.20                       

- Initial Acceptance to Bioconductor 

[leapR](/packages/leapR)
-----

                       Changes in version 0.99.6                        

- Built-in Plotting Function Created.

[lefser](/packages/lefser)
------

                       Changes in version 1.22.0                        

Minor improvements and bug fixes

- Fixed issue identifying classes in the ldaFunction effect size
calculation (@nnmbr, #88)
- Fix default colors arguments in lefserPlot(), lefserPlotFeat(), and
lefserPlotClad()
- Added label.font.color and ... (ellipses) parameters to
lefserPlot()
function for increased flexibility in customizing feature labels.
- Added label.font.face parameter to lefserPlot() function to allow
customization of font face for feature labels (e.g., "plain",
"italic", "bold", "bold.italic")

[lfa](/packages/lfa)
---

                 Changes in version 2.11.3 (2026-01-29)                 

- Function sHWE for BEDMatrix inputs now deletes temporary files
after
they are no longer needed

                 Changes in version 2.11.2 (2026-01-28)                 

- Updated vignette from old Sweave to modern R markdown.

                 Changes in version 2.11.1 (2026-01-28)                 

- Function sHWE greatly reduced memory handling for BEDMatrix inputs
by writing random genotypes to temporary plink1 BED files.
- Packages BEDMatrix and genio are now in "imports" (used to be
"suggests").
- Added option m_chunk to write temporary random genotypes in
large chunks for favorable I/O.
- Functions switched from deprecated to defunct: read.bed,
read.tped.recode, model.gof, center
- Package utils removed from "imports"

[limma](/packages/limma)
-----

                 Changes in version 3.68.0 (2026-04-28)                 

- 
  New arguments `offset` and `offset.prior` for voom(). voom()
  now supports offsets in DGEList objects.

- 
  New argument `treat.lfc` for topSplice() to allow TREAT
  fold-change thresholding with diffSplice() output.

- 
  In normalizeCyclicLoess(), default for `adaptive.span` becomes
  `TRUE`, and the `span` is chosen using the default settings of
  chooseLoessSpan().

- 
  The function topTableF(), which has been deprecated for many
  years, is now formally removed.

[limpa](/packages/limpa)
-----

                 Changes in version 1.4.0 (2026-04-29)                  

- 
  A limpa documentation page has been created at
  https://github.com/SmythLab/limpa/ containing a number of case
  study analyses of real datasets. This URL has been addded to
  the package DESCRIPTION file and the package vignette now has a
  "Further documentation" section that links out to this site.

- 
  Change of terminology from "peptides" to "precursors" as the
  default unit of input data throughout the package. More
  discussion of precursors vs peptides and other levels of data
  in the vignette.

- 
  Further clarification in the vignette that limpa can process
  any level of data for which intensities are available,
  including fragments, precursors, peptides, proteins, or PTMs.
  limpa can summarize precursors to either peptides or proteins
  and can conduct differential analyses at any level. The
  vignette now cites the Corso et al (2026) preprint as an IP-MS
  example. limpa can also explore isoform-specific expression by
  conducting differential usage analyses of peptides or PTMs
  within proteins, and there is a case study of such an analysis
  at https://github.com/SmythLab/limpa/.

- 
  As part of the above change, the output column
  `genes$NPeptides` from dpcQuant() nas been renamed to
  `genes$NPrec`, and the progress message from dpcQuant() and
  dpcQuantByRows() now says "Precursors" instead of "Peptides"
  for the number of rows. (The column name and message are
  generic: the column name will still be `NPrec` and the message
  will still say "Precursors" even if fragments are being
  processed instead of precursors.) Percentage progress has also
  been added to the progress message.

- 
  A new function EListFromLongFormatFile() has been added to read
  intensities from any long format file, and readDIANN() and
  readSpectronaut() now call this function instead of handling
  their own read code. The input file can be either an
  uncompressed delimited text file or a Parquet file, and the
  file type is detected automatically from the file name. The new
  function uses the nanoparquet package to read Parquet files
  whereas readDIANN() previously used arrow. The new function
  adds a number of new arguments and changes the names of some
  existing arguments. It allows filtering by any number of
  columns and allows reading of observation-level covariates. All
  three functions now allow `feature.column` to be a vector of
  arbitrary length, which allows columns to be combined to
  created a unique precursor or fragment ID. Amongst other
  things, this functionality facilitates reading of fragment
  level data and reading of MSstats format files. All three
  functions have a new argument `verbose` to optionally output
  progress information and the number of intensity values that
  have been filtered.

- 
  readDIAN() has a new argument list to match
  EListFromLongFormatFile(). The old argument `precursor.column`
  is renamed to `feature.column`, `qty.column` to
  `intensity.column`, and `extra.columns` to
  `annotation.columns`. Default setting for `q.columns` is
  updated.

- 
  readSpectronaut() has a new argument list to match
  EListFromLongFormatFile(). The old argument `precursor.column`
  is renamed to `feature.column`, `qty.column` to
  `intensity.column`, and `extra.columns` to
  `annotation.columns`. New arguments `filter.columns` and
  `censor.value` give more flexibility over filtering.

- 
  New function readSpectronautRunInfo(), which is called by
  readSpectronaut() to include run (sample) information in the
  `targets` component of the output.

- 
  New function plotAveVsMis() to plot average observed
  log-intensity per row vs the number of missing values.

- 
  dpc() has been replaced by a new function that combines the
  functionalities of dpcON() and dpcCN(), with the CN method as
  the default. The original dpc() function has been renamed to
  dpcLegacy().

- 
  dpcCN() now uses a systematic rather than a random sample of
  rows. Amongst other things, this ensures that the result does
  not change if the rows of `y` are reordered. The default value
  for `subset` is increased to 2000 and the default value for
  `verbose` is now `FALSE`.

- 
  The upper limit of the y-axis in the plot from plotDPC() is now
  extended slightly to allow for jittering of the detection
  proportions. Previously some of the jittered proportions above
  1 were not shown on the plot.

- 
  New argument `show.start` for plotDPC() to control whether the
  starting curve from dpc() is shown on the plot as well as the
  final plot. Default is now `FALSE`, whereas previously the
  starting curve was always shown.

- 
  dpcImpute() renamed to dpcQuantByRow() (although old function
  name will continue to work for backward compatibility).

- 
  When operating on an EList, dpcQuant() now removes the protein
  ID annotation column from the output object so that the IDs
  appear only as row.names. dpcQuant() also stores the DPC and
  the hyperparameters in the output EList object.

- 
  New function filterByDetection(), to choose proteins detected
  in at least a specified number of samples.

- 
  Axis labels produced by plotPeptides() and plotProtein() are
  now user-settable. plotPeptides() now plots sample names, in
  the same style as plotProtein(), and the y-axis is now labelled
  "Log-intensities" instead of "Log-expression".

- 
  New argument `las` for plotProtein() and plotPeptides() to
  allow vertical orientation of sample names on x-axis.

- 
  readFragPipe() and readMaxQuant() are now mentioned in the
  vignette.

- 
  Acknowledgement section added to vignette.

- 
  A number of programming efficiencies have been introduced to
  various functions.

[LimROTS](/packages/LimROTS)
-------

                       Changes in version 1.3.24                        

- Added differential expression analysis for repeated measures using
the dream linear mixed model.

                       Changes in version 1.3.22                        

- Incorporating repeated measurement analysis.

                       Changes in version 1.3.17                        

- Add survival analysis functionality with support for Cox
proportional hazards models and competing risks regression.

                       Changes in version 1.3.13                        

- Add copyright statements to R project files.

                       Changes in version 1.3.10                        

- Change the license from Artistic-2.0 to GPL-2.0-or-later.

                        Changes in version 1.3.6                        

- Adding Citation and Disclaimer Information.

[lncRna](/packages/lncRna)
------

                 Changes in version 0.99.3 (2026-03-11)                 

- Initial Acceptance to Bioconductor

[LRDE](/packages/LRDE)
----

                       Changes in version 0.99.7                        

- Initial Bioconductor submission of LRDE.

- Implements a hurdle negative binomial (Hurdle-NB) model for
  differential expression analysis of long-read RNA-seq data.

- Includes:
  - Data preparation (prepareDGE)
  - Size factor estimation (sizeFactorsEst)
  - Tag-wise dispersion estimation (tagwiseEst)
  - Differential expression testing via:
  * Likelihood Ratio Test (hurdle_LRT)
  * Wald Test (hurdle_Wald_Test)

- Supports matrix, data.frame, and SummarizedExperiment inputs.

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


[MAPFX](/packages/MAPFX)
-----

                 Changes in version 1.7.1 (2025-12-20)                  

- Corrected the version number in the inst/NEWS file (0.99.9 to 1.7.1).

[markeR](/packages/markeR)
------

                        Changes in version 1.1.2                        

Minor Changes

- Moved Python bridge scripts from inst/python/ to a top-level
python/
directory, as these are supplementary scripts not part of the R
package itself.
- Added requirements.txt to the python/ directory listing all needed
Python dependencies (rpy2, pandas, numpy, and optionally ipython and
jupyter) for easier environment setup.
- Removed redundant code snippets from the Python bridge scripts.

                        Changes in version 1.1.1                        

- Added p.adjust.method parameter across all functions performing or
depending on multiple testing correction, allowing users to specify
any correction method supported by stats::p.adjust(), beyond the
default Benjamini-Hochberg FDR.
- Added Python bridge scripts in inst/python/ for users who wish to
call markeR from a Python environment via rpy2. Includes a tutorial
workflow script and a generic command-line wrapper capable of
invoking any exported markeR function. See inst/python/README.md for
installation and usage.

[mastR](/packages/mastR)
-----

                       Changes in version 1.11.2                        

- Fix the problems caused by SeuratObject v5.

                       Changes in version 1.11.1                        

- Fix the problems caused by Seurat v5 and added a new vignette.

[matter](/packages/matter)
------

                       Changes in version 2.13.5                        

BUG FIXES

- Remove 'rescale_ref()' outside domain limits error

                       Changes in version 2.13.4                        

SIGNIFICANT USER-VISIBLE CHANGES

- Removes character vector support from 'bsearch()'

BUG FIXES

- Remove 'STRING_PTR' usage from C/C++ code

                       Changes in version 2.13.3                        

BUG FIXES

- Backport 'DATAPTR()' support for older R versions

                       Changes in version 2.13.2                        

BUG FIXES

- Remove 'DATAPTR()' usage from C/C++ code

                       Changes in version 2.13.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add argument for passing list of vectors to 'chunk_mapply()'

- Add argument for passing list of vectors to 'chunked_list()'

[MDSvis](/packages/MDSvis)
------

                       Changes in version 0.99.0                        

- Initial Acceptance to Bioconductor

[MeLSI](/packages/MeLSI)
-----
 
                       Changes in version 0.99.0                        

Initial Bioconductor Submission

- Initial submission of MeLSI (Metric Learning for Statistical
Inference) to Bioconductor
- Novel machine learning method for microbiome beta diversity
analysis
- Learns optimal distance metrics to improve statistical power in
detecting group differences
- Comprehensive validation against standard methods (Bray-Curtis,
Euclidean, Jaccard)
- Robust ensemble learning approach with conservative pre-filtering
- Validated on real microbiome datasets with proper Type I error
control
- Provides feature importance weights for biological interpretability
- Includes helper functions for CLR transformation and visualization
(VIP plots, PCoA)
- Full integration with Bioconductor ecosystem (phyloseq, microbiome
packages)

[memes](/packages/memes)
-----

                       Changes in version 1.19.1                        

- Removes ggseqlogo dependency as it is being archived from CRAN,
standardizes on universalmotif::view_logo for any internals that did
not already use it (thanks @bjmt for the implementation hint)
- @karawoo fixed a CI bug

[meshes](/packages/meshes)
------

                       Changes in version 1.37.1                        

- use 'enrichit' as engine for enrichment analysis (2025-12-07, Sun)

[metabinR](/packages/metabinR)
--------

                        Changes in version 2.0.0                        

Breaking changes

- The three binning entry points — abundance_based_binning(),
composition_based_binning(), hierarchical_binning() — now return a
MetabinResult S4 object instead of a plain data.frame. The v1.x
tabular layout is available via as.data.frame(result).
- The Java backend no longer exposes a command-line interface; the
shaded JAR dropped the commons-cli dependency. R is the only
supported entry point.
- Minimum Java runtime is now 17 (previously 8).

New features

- New MetabinResult class with accessors assignments(), nClusters(),
parameters(), algorithm(), plus show() and as.data.frame() methods.
- Binning functions accept, in addition to file paths,
Biostrings::DNAStringSet, Biostrings::QualityScaledDNAStringSet, and
ShortRead::ShortReadQ objects. Non-file inputs are staged to a
tempfile before the Java backend runs.
- numOfThreads defaults to BiocParallel::bpworkers().
- metabinR_jvm_options() and the metabinR.jvm.flags option let users
customise JVM flags at load time without losing the default G1GC +
string-deduplication tuning.
- Parameter validation is centralised via checkmate; errors surface
through cli::cli_abort() with classed conditions (metabinR_error_*).

Internal

- Java sources build via Maven (java/metabinR/pom.xml,
tools/build-jar.sh); the shaded JAR is reproducible
(project.build.outputTimestamp pinned).
- Typed JNI entry points on MTxAB, MTxCB, MTxABxCB return
tab-separated assignments directly to R, removing the previous
round-trip through on-disk CLI output.

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                        Changes in version 1.15                         

Changes in 1.15.3

- Use Spectra::rbindlistWithRownames() to expand on the changes
introduced in version 1.15.2 supporting also data.frames with row
names.

Changes in 1.15.2

- Use data.table::rbindlist() for merging lists of data.frame's to
improve performance.

Changes in 1.15.1

- Remove dependency on the msdata package: get test data from
MsDataHub.

[MetaboAnnotatoR](/packages/MetaboAnnotatoR)
---------------

                       Changes in version 0.99.14                       

- Initial Acceptance to Bioconductor

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                        Changes in version 1.19                         

Changes in 1.19.3

- Add functionality to translate between naming schemes used by
different software.

Changes in 1.19.2

- Add beta_values() function for chromatographic peak quality
assessment.

Changes in 1.19.1

- Update unit tests to replace deprecated functions.

[MetaboDynamics](/packages/MetaboDynamics)
--------------

                       Changes in version 2.1.104                       

- bug fix in estimates_dynamics

                       Changes in version 2.1.103                       

- bug fix in estimates_dynamics

                       Changes in version 2.1.101                       

- bug fix in fit_dynamics_model()
- new message in diagnostics_dynamics

                       Changes in version 2.1.100                       

- bug fix in plot_PPC for raw_plus_counts_model
- estimates_dynamics does not return 'sigma' or 'lambda' anymore

                       Changes in version 2.1.99                        

- added model_option in fit_dynamics_model(): either
"sd_per_time_point" (only option before) or "sd_per_condition" (less
sigmas, more robust with few replicates)

                        Changes in version 2.1.2                        

- bug fix in plot_cluster: correct sorting of time points

                        Changes in version 2.1.1                        

- cluster_dynamics bug fix: can handle all metabolite names

[metabom8](/packages/metabom8)
--------

                       Changes in version 0.99.10                       

- Initial Acceptance to Bioconductor

[MetaProViz](/packages/MetaProViz)
----------

                Changes in version 3.99.40 (2026-02-18)                 

- New cluster_pk() function for clustering prior knowledge terms
- New viz_graph() function for network graph visualization
- SummarizedExperiment support extended to viz_superplot()
- Bug fixes in processing(), viz_volcano(), dma(), and ORA functions
- NA handling improvements for correlation analysis
- License changed from GPLv3 to BSD 3-clause
- Code quality improvements for Bioconductor compliance
- Added ggraph dependency
- Updated documentation and argument descriptions

                Changes in version 3.99.10 (2025-10-28)                 

- SummarizedExperiment support as default input and output
- First candidate for Bioconductor submission

                 Changes in version 3.99.1 (2025-10-08)                 

- Fixed Bioconductor timeout issues
- Code refactoring to meet Bioconductor quality standards

                 Changes in version 3.99.0 (2025-09-30)                 

- General refactoring to meet Bioconductor checks
- RCMD check improvements
- Updated GitHub actions workflows


[metaseqR2](/packages/metaseqR2)
---------

                 Changes in version 1.23.1 (2026-03-04)                 

NEW FEATURES

- None.

BUG FIXES

- Moved several js libraries locally because of policy changes in
  highcharts

[methyLImp2](/packages/methyLImp2)
----------

                 Changes in version 1.7.1 (2025-11-12)                  

- Removing documentation of internal functions from the manual.

[MetMashR](/packages/MetMashR)
--------

                        Changes in version 1.5.1                        

- Rebuild documentation
- move cowplot to Suggests
- add ggplot2 to namespace imports

[mia](/packages/mia)
---

                        Changes in version 1.19                         

- Fixed a bug related to agglomerateByPrevalence and referenceSeq
  agglomeration (2025-12-08)

- Use ecodive in UniFrac (1.19.2, 2025-12-09)

- Added binning transformation (1.19.3, 2026-02-23)

- Added support for reducedDim in clustering (1.19.5, 2026-03-06)

[miaViz](/packages/miaViz)
------

                        Changes in version 1.19                         

- plotRDA: Now plotting works with interaction term (2025-11-09)

- plotBoxplot: Added option to add p-values (2026-01-07)

[microbiome](/packages/microbiome)
----------

                 Changes in version 1.32.2 (2025-07-17)                 

- Fixed NEWS file

[MicrobiomeProfiler](/packages/MicrobiomeProfiler)
------------------

                       Changes in version 1.17.1                        

- use 'enrichit' as engine for enrichment analysis (2025-12-07, Sun)

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.23.1                        

- fixed the issue of mp_import_qiime2 due to the new biomformat.
(2026-04-02, Thu)

                       Changes in version 1.23.0                        

- Bioconductor 3.22 released, and Bioconductor 3.23 (devel) bump.
(2025-10-31, Fri)

[MIRit](/packages/MIRit)
-----

                        Changes in version 1.7.4                        

This new release of MIRit fixes the rendering of the vignette and
updates the documentation website. GitHub worflows have also been
updated.

                        Changes in version 1.7.3                        

This version fixes some minor bugs in the integration analysis, and
adds
references to the new manuscript in the documentation, README,
CITATION,
and vignette. A new setTargets() function has also been included.

[miRSM](/packages/miRSM)
-----

                        Changes in version 2.7.1                        

- Update miRSM.R <2026-01-15, Thus>

[miRspongeR](/packages/miRspongeR)
----------

                       Changes in version 2.15.2                        

- Update moduleDEA.Rd <2026-01-11, Sun>.

                       Changes in version 2.15.1                        

- Update miRspongeR.Rmd <2025-12-23, Tues>.

[mobileRNA](/packages/mobileRNA)
---------

                        Changes in version 1.7.0                        

- Added Namespace dependency to DESCRIPTION Imports/Depends entries:
‘Seqinfo’

[monaLisa](/packages/monaLisa)
--------

                       Changes in version 1.17.1                        

- clarifications of supported families for the stability selection

[Moonlight2R](/packages/Moonlight2R)
-----------

                        Changes in version 1.9.4                        

Summary

- Made p-value tests for FEA(method='fgsea') more qualitative to
address low reproducibility of gfsea()

                        Changes in version 1.9.3                        

Summary

- Fixed tests for FEA with "fgsea" as method and added new parameter
('seed') to the FEA() function for setting random seed and control
stochastic perumation-based fgsea scores and p-values.

                        Changes in version 1.9.2                        

Summary

- Added new dataset for proteomics data and implemented fgsea method
in the FEA function

                        Changes in version 1.9.1                        

Summary

- Removed hotfix in GSEA() module

                        Changes in version 1.9.0                        

Summary

- Updated version for BioC release

[MotifPeeker](/packages/MotifPeeker)
-----------

                        Changes in version 1.3.3                        

Miscellaneous

- Add unsupported platform to DESCRIPTION.

                        Changes in version 1.3.2                        

Bug Fixes

- Bound maximum bootstrapping population size.

Misclaneous

- Update CITATION: MotifPeeker is now published on Bioinformatics
Advances!

                        Changes in version 1.3.1                        

Section Rework

- Motif-summit distances
- Remove distances by peak count.
- Reword descriptions
- &#91;NEW&#93; Add bootstrapping to visualise the distribution of
motif-summit distances.

New Features

- Add support for directly importing "mm10" and "mm39" mouse genome
builds.
- read_motif_file can now accept an universalmotif object, returns
the
same object.
- Input datasets section now reports datasets in a tabular form.

Documentation

- Improve instructions for importing genome builds.
- Improve instructions for importing motifs.

[motifTestR](/packages/motifTestR)
----------

                        Changes in version 1.7.2                        

- Changed default for min_score to 50% instead of 80%

[MouseFM](/packages/MouseFM)
-------

                 Changes in version 1.21.1 (2026-03-22)                 

- Improved handling of Ensembl BioMart host and version mismatches

- Added more informative error messages for BioMart connection and
  genome version issues

- Prevented online example execution during package checks

[msa](/packages/msa)
---

                       Changes in version 1.43.1                        

- fix of bug that deleted sequence names if the sequences were read
  from
  FASTA file (all three alignment algorithms)

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                        Changes in version 1.19                         

Changes in 1.19.3

- Use rbindlistWithRownames() for faster concatenation of results.

Changes in 1.19.2

- Refactor metaDataBlocks() adding parameters for the individual
metadata blocks simplifying configuration of data import.
- Fix importing issues of individual optional blocks. Import was
validated using MassBank releases 2025.05.1 and 2025.10. Several
fields/spectra variables now return a list to support records with
multiple values per field. These are: "name", "chrom_solvent",
"comment", "data_processing_comment", "sample",
"data_processing_reanalyze", "data_processing_whole".
- Map field MS$DATA_PROCESSING: FIND_PEAK to a spectra variable with
name "data_processing_find_peak" (was "data_processing_find").

Changes in 1.19.1

- Avoid error when "N/A" is reported as retention time.
- Performance improvement: use data.table::rbindlist() for merging
and
importing records from files in MassBank format.

[MsBackendMetaboLights](/packages/MsBackendMetaboLights)
---------------------

                        Changes in version 1.5.3                        

- Use MsCoreUtils::retry() instead of the internal implementation.
- Improve fail-safe mechanism when caching MS data from MetaboLights.

                        Changes in version 1.5.2                        

- Add helper functions mtbls_assay_data(), mtbls_sample_data() and
mtbls_metadata().

                        Changes in version 1.5.1                        

- Increase number of download retries and increase waiting time in
between.
- Export the retry() function.

                         Changes in version 1.5                         

[MsBackendMgf](/packages/MsBackendMgf)
------------

                        Changes in version 1.19                         

Changes in 1.19.1

- Use data.table::rbindlist() instead of MsCoreUtils::rbindFill() to
combine individual spectra's data.frames into a single one. This can
have performance improvements, in particular for large MGF files.

[MsBackendMsp](/packages/MsBackendMsp)
------------

                        Changes in version 1.15                         

Changes in 1.15.1

- Use data.table::rbindlist() to combine individual spectra into one
resulting data.frame/DataFrame in readMsp(). This improves the
performance of readMsp() as well as backendInitialize() in
particular for large MGF files.

Changes in 1.11.1

- Complete unit test coverage.

[MsBackendSql](/packages/MsBackendSql)
------------

                        Changes in version 1.11                         

Changes in 1.11.3

- Remove dependency from the msdata package: load test data from the
MsDataHub package.

Changes in 1.11.2

- Small update to the internal function to extract the spectra data:
avoid fmatch() call if not needed.
- Pass sorted spectra IDs to the SQL query, which can result in a
tiny
performance gain.

Changes in 1.11.1

- For storage mode of peaks data in long form (i.e., one row per
peak), create a incremental (unique) peak_id_ database column if the
database requires that (e.g. for DuckDb).
- longForm() on a database requiring peak_id_ (e.g. DuckDb) order
results based on the peak_id_ column.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                        Changes in version 1.23                         

MsCoreUtils 1.23.10

- Fix partial argument match warnings.
- Fix bug in impute_mixed() and allow for two MARGINs for mixed
imputation.
- Remove the ... in impute_mixed()`.

MsCoreUtils 1.23.9

- Handle all type of error in retry(), not only simpleError.

MsCoreUtils 1.23.8

- Add parameters warningsAsErrors and verbose to retry().

MsCoreUtils 1.23.7

- Add retry() function.

MsCoreUtils 1.23.6

- Fix unsupported Unicode characters in documentation.
- Fix join_gnps() to support also type being different from "outer".

MsCoreUtils 1.23.5

- rbindFill(): improve performance when joining only matrices.

MsCoreUtils 1.23.4

- gnps() and join_gnps() use C implementations for modified cosine
similarity calculation. The original R implementations are available
as gnps_r() and join_gnps_r(). See issue #131 for discussion and
performance comparison.
- FastCosine spectral similarity calculation implementation:
gnps_chain_dp().

MsCoreUtils 1.23.3

- Add parameter matchedPeaksCount to spectra similarity/distance
functions to report, in addition to the similarity, also the number
of matched peaks.

MsCoreUtils 1.23.2

- Fix Found non-API call to R: ‘SETLENGTH’ in C_reduce, C_join_inner,
and C_join_outer by using lengthgets() instead (see issue 136).
- Use Rf_ prefix for R API functions and set R_NO_REMAP to avoid use
of accidentially using non-Rf_ prefixed functions in future.

MsCoreUtils 1.23.1

- Fix RF imputation, that now needs dimnames.

MsCoreUtils 1.23.0

- New devel version

[MsDataHub](/packages/MsDataHub)
---------

                       Changes in version 1.11.5                        

- Fix metadata

                       Changes in version 1.11.4                        

- Add DDA and DIA proteomics data to illustrate
QFeatures::readQFeatures().
- Various fixes.

                       Changes in version 1.11.3                        

- Add Boekweg et al. (2022) SCP and bulk mzML files, and Sage
identification TSV files (see ?Boekweg2022)
- Added Guillaume Deflandre as contributor.

                       Changes in version 1.11.2                        

- Add MS3TMT data files.

                       Changes in version 1.11.1                        

- Add CE-MS test data files (see PR #11).
- Add description of the data sets to the package vignette.

                        Changes in version 1.11                         

[MsExperiment](/packages/MsExperiment)
------------

                        Changes in version 1.13                         

Changes in 1.13.1

- Load test data from MsDataHub and remove dependency on the msdata
package

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.37                         

2.37.3

- Bump version to make sure the latest updates get published on Bioc.

2.37.2

- Fixing vignette (see #613)
- use MsDataHub for TMT and PestMix1_DDA data.

2.37.1

- Fixing tests

2.37.1

- Remove rols dependency.
- Remove pryr dependency.

[msPurity](/packages/msPurity)
--------

                       Changes in version 1.37.3                        

- Update within the tests only (not within the R code base): library
  database path for lvl spectral matching testing (required after
  removal of the local SQLite database from msPurityData).

                       Changes in version 1.37.2                        

- Bug fix for subsetting the grouped data frame in the
  `average_xcms_grouped_msms_indiv` function

                       Changes in version 1.37.1                        

- Removed library_spectra.db direct reference from msPurityData package
  (too large for Bioconductor)

- spectralMatching() now downloads default library database from Zenodo
  (https://zenodo.org/records/18700802) and uses BiocFileCache with MD5
  verification

- spectral_matching() (deprecated function) now requires explicit
  library_db_pth parameter; default library no longer available

- Added BiocFileCache to Imports for automatic caching of remote
  databases

[msqrob2](/packages/msqrob2)
-------

                        Changes in version 1.19                         

msqrob 1.19.1

- Add functions to calculate log-normalisation factors that can be
used with sweep function
- Add function msqrobCollect to collect results tables
- Add function plotVolcano to make volcanoplots for hypothesisTest
tables
- Add function createPairwiseContrasts to generate contrasts for all
pairwise comparisons
- Illustrate these new functionalities in the vignettes
- New Vignettes for DIA (DIA-NN and Spectronaut)

msqrob2 1.19.2

- Update cptac vignette (fread with argument check.names = TRUE)

msqrob2 1.19.3

- Update vignettes (arguments fread check.names = TRUE and integer64 = "double") to avoid issues with readQFeatures
- Update README to point to msqrob2book

[MsQuality](/packages/MsQuality)
---------

                        Changes in version 1.11                         

1.11.3 (2026-03-10)

- Implementation of per-chromatogram metrics for Chromatograms
objects: maxIntensity, intensityMean, intensitySd,
intensityQuartiles, intensityRange, peakCount, baselineIntensity,
signalToNoiseRatio, xicFwhm, peakBoundary, peakWidth,
gaussianSimilarity, peakProminence
- Shared functions (chromatographyDuration, rtAcquisitionRange,
rtIqr,
areaUnderTic, areaUnderTicRtQuantiles, ticQuantileRtFraction,
msSignal10xChange, medianTicRtIqr, numberEmptyScans,
calculateMetrics) now support both Spectra and Chromatograms
- msSignal10xChange gains minIntensity parameter for noise
robustness.
- Renamed argument filterEmptySpectra to filterEmptyObject
- Fixed calculateMetrics() data.frame output
- Fixed shinyMsQuality() documentation example

Changes in version 1.11.2 (2026-03-02)

- replace msdata by MsDataHub package (issue #19)
- adjust all examples, unit tests and vignette to the MsDataHub
package

Changes in version 1.11.1 (2026-02-17)

- added conversion to uri from local file to the export function in
transformIntoMzQC (contribution by Helge Hecht, PR #20)

Changes in version 1.11.0 (2025-10-29)

- bump version to 1.11.0

[MSstatsConvert](/packages/MSstatsConvert)
--------------

                       Changes in version 1.22.0                        

Protein Turnover Support

- DIA-NN: DIANNtoMSstatsFormat now accepts a labeledAminoAcids
parameter to enable protein turnover workflows. Provide the
single-letter amino acid codes that carry the SILAC label (e.g.
c("K") or c("K", "R")). The converter automatically distinguishes
heavy and light peptide forms — either from a Channel column in
DIA-NN 2.x exports, or by parsing SILAC modification tags (e.g.
(SILAC-K-H)) in DIA-NN 1.x exports — and populates the
IsotopeLabelType column required by MSstats downstream.

- Spectronaut: SpectronauttoMSstatsFormat now accepts a heavyLabels
parameter to classify peptides in protein turnover experiments.
Supply the label names as they appear in square brackets in the
peptide sequence (e.g. c("Lys6") or c("Lys6", "Arg10")). Peptides
are automatically assigned as heavy (H), light (L), or unlabeled.
Any novel label name reported by Spectronaut (e.g. "Leu6", "Phe10")
is supported.

- Spectronaut: A new peptideSequenceColumn parameter lets you specify
which column in your Spectronaut export contains the peptide
sequence. This is especially useful for protein turnover reports,
which may use a different column layout than a standard Spectronaut
output.

- Spectronaut: Protein turnover reports that lack fragment-level
columns (e.g. FFrgLossType, FFrgIon) are now handled automatically —
previously these required manual pre-processing before import.

- Fraction selection: The algorithm that selects the best fraction
for
each peptide feature now correctly accounts for heavy and light
isotope channels in protein turnover data, ensuring that
light-channel coverage is prioritized when choosing fractions. When
multiple fractions have identical coverage, the fraction with the
highest mean intensity is chosen.

DIA-NN Quality Control: Anomaly Detection

- DIANNtoMSstatsFormat now supports automated data quality scoring
via
an isolation-forest anomaly detection model (calculateAnomalyScores
= TRUE). When enabled, each precursor-level feature receives an
anomaly score that can be used downstream in MSstats to down-weight
measurements that look like outliers. Users supply the quality
metric columns to use as model features (anomalyModelFeatures) and,
optionally, a run order table (runOrder) to engineer temporal trend
features that detect gradual instrument drift across an experiment.

Spectronaut: MS1 Quantification

- SpectronauttoMSstatsFormat now accepts "FG.MS1Quantity" (and other
raw Spectronaut column names) as the intensity argument, in addition
to the existing "PeakArea" and "NormalizedPeakArea" options.

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.38.0                        

Bug fixes and minor improvements

- Fixed bug in getWithColData when colData from i object did not
align
with sampleMap(mae) colname order (@aboyoun / James Bonaffini, #346)

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.19.4                        

Bug fix

- convert p-values from 0 to 1e-300 when imputing z-score

                       Changes in version 1.19.3                        

Bug fix

- Add more documentation for chr checks

                       Changes in version 1.19.2                        

Bug fix

- liftover() couldn't previously handle CHR outside of 1-22, added
functionality following the same logic as format_sumstats()

                       Changes in version 1.19.1                        

Bug fix

- Update approach to download from ieugwasr.

[muscat](/packages/muscat)
------

                       Changes in version 1.25.4                        

- update Gilis et al. citation for DD

- remove 'cowplot' in favor of 'patchwork'

- move 'IHW', 'DESeq2' and 'sctransform' to 'Suggests:'

- omit 'data.table' and 'purrr' everywhere and remove from 'Imports:'

- #102 fix: use correct CPM values in 'resDS' when 'method=limma-voom'

                       Changes in version 1.25.3                        

- up R version dependency to ≥ 4.6

- rename 'edgeR::calcNormFactors' to 'normLibSizes'

                       Changes in version 1.25.2                        

- bug fix in '.mm_dream' (underlying 'mmDS' methods 'poisson' and
  'nbinom'):
  transform size factors from linear- to log-scale before being used as
  offset

                       Changes in version 1.25.1                        

- fix #146: support for sparse matrices by 'pbDS'

- updated Germain et al. citation for BBHW

[MutSeqR](/packages/MutSeqR)
-------

                 Changes in version 0.99.4 (2025-12-09)                 

Preparing for Bioconductor release. We address comments from
Bioconductor reviewers. Major changes include parameter validation
and
vectoriation of functions. filter_mut() "snv_in_germ_mnv" parameter
was
fixed such that it now only filters out overlapping snvs if their
variation matches that of the germline mnv (previously it was blanket
removing all snvs that overlap with germline mnvs). plot_lollipop()
can
now colour the plot by different subtype resolutions (previously just
base_6). Fixed plot_spectra() axes labels after clustering
(previously
not showing labels).

                 Changes in version 0.99.3 (2025-10-10)                 

Fix dependencies for Bioconductor release.

                 Changes in version 0.99.2 (2025-09-15)                 

Preparing for Bioconductor release. This change removes bmd_toxicR
from
the package. ToxicR dependency is not supported by Bioconductor. See
ToxicR_archive branch for bmd_toxicr function. Modifies
signature_fitting() to use SigProfilerMatrixGenerator python
dependency
rather than R dependency.

                 Changes in version 0.99.1 (2025-08-13)                 

Preparing for Bioconductor release. This change adds MutSeqRData to
suggests, and alters how the examples are run (now depends on the
ExperimentHub accessions).

                 Changes in version 0.99.0 (2025-06-19)                 

Initial public version.

Major changes

- Added filter_mut() to workflow: germline identification via
vaf_cutoff, region filtering, and depth correction now occur here
instead of the import functions.
- calculate_mut_freq() is renamed to calculate_mf().
- calculate_mf() no longer requires depth; users may:
1.  calculate depth from mutation data,
2.  supply a separate depth table, or
3.  omit depth entirely (only mutation counts returned).
- correct_depth option moved to calculate_mf().
- plot_spectra(), plot_trinucleotide(), and spectra_comparison() now
use mf_data instead of raw mutations.
- Output options added: VCF, FASTA, SigProfiler-compatible format,
Excel workbook.
- Example dataset (~44MB) added.

New features

- render_report() added for standardized summary reporting.

Other

- Removed custom_regions parameter; replaced by generalized regions
argument.
- Public release 🎉
- See the vignette for details.

[mzR](/packages/mzR)
---

                       Changes in version 2.45.1                        

- Load test files from the MsDataHub package and remove dependency on
  the
  msdata package.

- Catch error in spectrum ID parsing.

- Remove support for C++11 and C++14 from R-devel (see issue #315)

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 3.8.0                        

- Addes sorting to the heatmap so reads appear from most to least
methylated. This is done by calculating the average methylation
across the region for each read and sorting by this value.

[netboost](/packages/netboost)
--------

                 Changes in version 2.19.3 (2026-02-12)                 

- Consensus network calling.

                 Changes in version 2.19.2 (2026-02-11)                 

- Speed up of filter calculation.

[NormalyzerDE](/packages/NormalyzerDE)
------------

                       Changes in version 1.29.2                        

- Update maintainer list, adding Måns Zamore as an author

[notame](/packages/notame)
------

                        Changes in version 1.1.0                        

New features

- All feature data columns that differ between batches are now retained
  when merging batches: the most common flag value across batches is
  included in the "Flag" column in the merged object

- Added function `read_from_msdial()` which reads MS-DIAL output files
  (tab-separated .txt files) into a SummarizedExperiment object

- Improved type conversion when importing data with notame functions

Bug fixes

- Fixed an issue where metadata was not properly retained when joining
  new data with join_rowData() and join_colData()

[notameViz](/packages/notameViz)
---------

                        Changes in version 1.1.0                        

New features

- plot_injection_lm() now utilizes limma for p-value calculation,
  resulting in improved performance (instant plotting)

- Merging QC plots in save_QC_plots() does not require external
  software in Windows/Linux

Bug fixes

- save_QC_plots() now correctly uses "QC" as default color option

- Fixed an issue that prevented using custom color scale in
  plot_effect_heatmap()

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

[omXplore](/packages/omXplore)
--------

                        Changes in version 1.5.1                        

- 1.5.1 - Replaced highcharter by plotly

- This is the first release of omXplore.

[OncoSimulR](/packages/OncoSimulR)
----------

                 Changes in version 4.13.4 (2026-04-16)                 

- Fixed miscell minor notes and warnings from
  r-universe.

                 Changes in version 4.13.3 (2026-04-14)                 

- Upgraded exprtk to version from 2025-01-01,
  commit cebb369.

- And proper use of all.equal with check.attributes = FALSE!

                 Changes in version 4.13.2 (2026-04-14)                 

- Fixed the "More than one VignetteEngine specified"
  error in r-universe.

- Added unexported complete_fitness_landscape function.

- Added a few magellan epistasis tests.

- Removed an equality comparison of floats.

                 Changes in version 4.13.1 (2026-01-31)                 

- Fixed NOTES in "R CMD check --as-cran"

[orthogene](/packages/orthogene)
---------

                       Changes in version 1.17.4                        

Bug fixes

- run_benchmark_once(): replace message(e) with
message(conditionMessage(e)) inside the tryCatch error handlers.
Passing a condition object to message() re-signals it, and under
testthat's calling handlers an error condition signaled this way
escapes the tryCatch and is reported as a test failure — even though
tryCatch did catch it. This caused test-run_benchmark to fail on
GitHub Actions whenever the g:Profiler API returned malformed JSON,
even though run_benchmark_once() was supposed to record NA and move
on.

                       Changes in version 1.17.3                        

Bug fixes

- Make examples of convert_orthologs(), map_orthologs(), and
all_genes() use method = "homologene" so they don't depend on the
g:Profiler API. The remote g:Profiler service has been
intermittently returning JSON with literal NaN tokens, which
jsonlite cannot parse, causing the Bioconductor R CMD check for
examples to fail.
- Tests that incidentally exercised the default g:Profiler path
(test-aggregate_mapped_genes, test-convert_orthologs,
test-infer_species) now request method = "homologene" explicitly.
- Tests that intentionally cover the g:Profiler backend
(test-all_genes_gprofiler, the method = "gprofiler" block in
test-convert_orthologs) now skip gracefully with a clear message
when the API is unavailable, rather than halting the whole suite.
- test-run_benchmark: aggregate mean() checks now use na.rm = TRUE so
a transient g:Profiler failure doesn't NA-poison the benchmark
assertions.
- test-infer_species: tolerate ties in top_match so that when two
species are equally well-matched (e.g. human and monkey via
babelgene), the test passes as long as the expected species is among
them.

                       Changes in version 1.17.2                        

New features

- all_genes: enable method/species-specific caching.
- run_benchmark_once: enable caching
- cache_dir(): new function for consistency
- map_orthologs_gprofiler: BIG upgrade! Can now rapidly query large
gene lists via batching and parallelzation without server-side
timeout errors.
- Improve test coverage.

[OUTRIDER](/packages/OUTRIDER)
--------

                       Changes in version 1.28.1                        

- Adapt to new testthat behaviour (#87)

[parati](/packages/parati)
------

                 Changes in version 0.99.8 (2026-04-12)                 

- Initial Acceptance to Bioconductor

[pcaExplorer](/packages/pcaExplorer)
-----------

                        Changes in version 3.6.0                        

Other notes

- pcaplot() labels now all samples with a segment connector in grey
to
increase its visibility on a bw/classic theme

[Pedixplorer](/packages/Pedixplorer)
-----------

                        Changes in version 1.7.1                        

- Fix slices in circfun()
- Add contributing guidelines
- Replace dplyr pipe %>% by |>

[phantasus](/packages/phantasus)
---------

                       Changes in version 1.31.1                        

- dedicated docker server script to reduce memory usage

- fixes in annotation parsing

[philr](/packages/philr)
-----

                       Changes in version 1.37.1                        

USER-VISIBLE CHANGES

- Updating documentation (and vignette) to make input format for philr
  function clearer.

[PIUMA](/packages/PIUMA)
-----

                        Changes in version 1.8.0                        

- PIUMA now identifies clusters and can be easily integrated in the
  Seurat
  workflow

[plaid](/packages/plaid)
-----

                       Changes in version 0.99.0                        

INITIAL RELEASE

- plaid: ultra-fast single-sample enrichment scoring using rank-based
  statistics

- replaid methods: fast implementations of popular enrichment methods
  (ssGSEA, GSVA, singscore, UCell, AUCell, scSE)

- dualGSEA: statistical testing for differential enrichment with
  dualGSEA

- Full Bioconductor integration with SummarizedExperiment,
  SingleCellExperiment, and BiocSet support

- GMT utilities for reading, writing, and converting gene set files

- Performance optimizations including sparse matrix support and
  parallel
  processing

[PlinkMatrix](/packages/PlinkMatrix)
-----------

                       Changes in version 0.99.5                        

- Initial version

[PLSDAbatch](/packages/PLSDAbatch)
----------

                       Changes in version 1.99.1                        

- Date: 2026-01-06
- Text: Fix the figure display issue
- Details: Fix the issue where Figures 2 and 3 are not shown in the
vignette.

                       Changes in version 1.99.0                        

- Date: 2025-12-26
- Text: improve usability
- Details:
- Added a mode argument to PLSDA_batch().
- Added a criterion argument to linear_regres() to select P-values
from the optimal model based on the specified criterion.
- Added a return.model arugument to linear_regres() to reduce
memory usage when set to FALSE.
- Extended Scatter_Density() to support any multivariate method
that returns component scores, including PCA and PLS, with
corresponding arguments updated.
- Added lighten() and darken() functions for enhanced color
generation.
- Refined multiple functions to improve usability.
- Updated the vignette accordingly.

[plyranges](/packages/plyranges)
---------

                       Changes in version 1.32.0                        

- added new functionality for slice_*, as in dplyr, including head,
tail, max, min and sample. These include the dplyr arguments such as
n, prop, and weight_by (for sampling).
- added join_mcols_left(x,y) for GRanges x and tabular y
- added distance=TRUE option for join_overlap_*

                       Changes in version 1.31.5                        

- added distance=TRUE option for join_overlap_*
- added join_mcols_left(x,y) for GRanges x and tabular y

[pmp](/packages/pmp)
---

                       Changes in version 1.23.1                        

- fix broken test due to changes in RF imputation

[posDemux](/packages/posDemux)
--------

                       Changes in version 0.99.9                        

- Initial Acceptance to Bioconductor 

[postNet](/packages/postNet)
-------

                       Changes in version 0.99.9                        

- Initial submission to Bioconductor

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.51                         

Changes in version 1.51.1

- Bump version for Bioc devel

Changes in version 1.51.0

- New version for devel

[pRolocGUI](/packages/pRolocGUI)
---------

                       Changes in version 2.21.2                        

- New version for Bioc release

                       Changes in version 2.21.1                        

- Fix bug in explore app

                        Changes in version 2.21                         

                       Changes in version 2.21.0                        

- New version for Bioc devel 3.23

[psichomics](/packages/psichomics)
----------

                       Changes in version 1.36.1                        

- GTEx data loading (loadGtexData()):
- Support for loading GTEx V10 data
- Documentation:
- Remove confusing instructions on launching psichomics from
README
- Improve Docker instructions on README
- Update external function documentation
- Update license copyright years
- Update email contact
- Bug fixes:
- Fix interface not loading for collapse panels
- Fix eBayes error on graphical interface
- Fix unit tests when testing shiny buttons

[PSMatch](/packages/PSMatch)
-------

                        Changes in version 1.15                         

PSMatch 1.15.3

- Adjusted documentation on addCarbamidomethyl = TRUE in
calculateFragments() set by default.
- Corrected plotSpectraPTM() relative PTMods dependencies,
identifications are highlighted in bold in the USI. Added parameters
to call addFixed and addVariable within plotSpectraPTM().
- Add PTMods dependency and thus positional modifications in
calculateFragments() (see issue 38)

PSMatch 1.15.2

- Use MsDataHub instead of msdata.

PSMatch 1.15.1

- Fix minor check notes (Rd link and IRanges import)

PSMatch 1.15.0

- New devel version

[PTMods](/packages/PTMods)
------

                        Changes in version 0.99                         

- Initial Acceptance to Bioconductor

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.21                         

QFeatures 1.21.4

- Add new readQFeautres vignette.

QFeatures 1.21.3

- Add readQFeaturesFromDIANN(multiplexing = 'dimethyl') example.
- New replaceColnames() function (see Issue 258)

QFeatures 1.21.2

- Add support for dimethyl multiplexing to readQFeaturesFromDIANN()
(see PR #249), contributed by Karolína Kryštofová.
- Aggregate redundant messages in aggregateFeatures() (PR #255).
- Add a progress bar to aggregateFeatures() and readQFeatures() (PR
#255).

QFeatures 1.21.1

- Link to RforMS contribution guide.

QFeatures 1.21.0

- New devel version

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

                       Changes in version 1.29.3                        

- Update maVolPlot man page

                       Changes in version 1.29.2                        

- Update to getContrasts man page example to account for changes in
  limma::topTable

                       Changes in version 1.29.1                        

- Update to getContrasts to account for changes in limma::topTable

- Update to maVolPlot to account fo changes in limma::topTable

[QTLExperiment](/packages/QTLExperiment)
-------------

                 Changes in version 2.3.1 (2026-01-06)                  

- Updated documentation for sumstats2qtle() to remove warning.

[QUBIC](/packages/QUBIC)
-----

                       Changes in version 1.39.0                        

- Runtime dependency on biclust is removed

[queeems](/packages/queeems)
-------

                       Changes in version 0.99.5                        

- Initial Acceptance to Bioconductor 

[RankMap](/packages/RankMap)
-------

                       Changes in version 0.99.1                        

Initial submission to Bioconductor

New Features

- Fast, robust, and scalable reference-based cell type annotation
using multinomial regression on ranked expression matrices.
- Supports both single-cell and spatial transcriptomics data.
- Compatible with Seurat, SingleCellExperiment, and SpatialExperiment
objects.
- Core function RankMap() provides a streamlined pipeline for
preprocessing, model training, and prediction.
- Customizable preprocessing: top-K gene masking, optional binning,
expression weighting, and scaling.
- Additional functions:
- computeRankedMatrix() – generate ranked matrices
- trainRankModel() – train multinomial GLM
- predictRankModel() – apply trained model to query data
- evaluatePredictionPerformance() – assess accuracy
- Optimized for large datasets with significantly faster runtime than
SingleR, Azimuth, and RCTD.

[Rarr](/packages/Rarr)
----

                        Changes in version 1.99                         

Breaking changes

- The DelayedArray backend (writeZarrArray() and ZarrArray()
functions) has been migrated to a separate, dedicated package. This
reduces the number of dependencies from 37 to 24. This also greatly
improves performance in for the standard case (when the DelayedArray
backend is not used).
- write_zarr_array() now writes Zarr v3 by default. Writing Zarr v2
is
still possible by explicitly setting the argument zarr_version = 2.

New features

- Zarr v3 arrays with data types and codecs that already existed in
v2
can now be read via read_zarr_array(), and written via
write_zarr_array().
- Zarr v3 consolidated metadata is now returned by zarr_overview(),
the same way it was already previously done for v2 consolidated
metadata.
- More data types are available when writing Zarr arrays:
- boolean / logical
- int8
- int16
- int64 (up to values that can be represented as R integers)
- uint8
- uint16
- uint32 (up to values that can be represented as R integers)
- uint64 (up to values that can be represented as R integers)
- float32 / single
- Scalar arrays (i.e., arrays with zero dimensions) can now be read.
Thanks to Artür Manukyan for the bug report.
- Zarr attributes can now be read by passing an s3 URL directly as
the
first argument of read_zarr_attributes(). This makes
read_zarr_attributes() consistent with read_zarr_array() and
zarr_overview().
- "Simple" structured data types (i.e., only one level of nesting and
no arrays) can now be read from Zarr v2 arrays.
- simplifyVector = FALSE is added to fromJSON in
read_zarr_attributes(), thus attributes of both local and s3 zarr
stores are read identically.
- The dimension_names optional field is support in both v2 (not
strictly part of the spec) and v3. It is mapped to
names(dimnames(.)) in R.
- NA_real_ is now an allowed fill value in write_zarr_array() when
writing numeric arrays, following a request from Hervé Pagès.
- Fill values stored as their byte representation are now understood
when reading Zarr arrays.
- write_zarr_array() now supports writing NA_character_, which means
it is possible to preserve NAs when roundtriping an R character
array, based on a request from Hervé Pagès.

Minor improvements

- There is now a dedicated vignette describing the supported Zarr
features in Rarr, available at
https://huber-group-embl.github.io/Rarr/articles/features.html. This
makes it more easily discoverable on the Bioconductor landing page.
- Rarr initializes empty/missing chunks only once per read operation,
which significantly improves performance when reading arrays with
many missing chunks.
- Reading fixed-length string and unicode arrays is now ~20% faster.
- The shape and chunks fields in v2 metadata are now always encoded
as
JSON arrays, even when they contain a single element. This makes
Rarr more compatible with other Zarr implementations. Thanks to
Artür Manukyan for the bug report and pull request.
- Empty zarr arrays (i.e., arrays with shape and chunks equal zero)
can now be written.
- Compression for writing Zarr arrays now default to zstd rather than
zlib. zstd achieves similar or better compression levels while being
much faster at compressing (= writing Zarr arrays) and decompressing
(= reading Zarr arrays). This matches the default used by Zarr
Python implementation.
- write_zarr_array() now fails early with an explicit error message
when x is not an array.

Bug fixes

- Rarr is now fully compatible with big endian platforms.
- ZSTD decompression now also works in case where we cannot guess a
priori the buffer size from the data type, such as when using
variable length strings. Thanks to Artür Manukyan for the bug report
and test data.
- zarr_overview() no longer fails on consolidated metadata containing
uncompressed arrays. This was introduced in
https://github.com/Huber-group-EMBL/Rarr/pull/45. Thanks to Sharla
Gelfand for reporting the issue and providing test data.
- the fill_value is now correctly interpreted when reading Zarr v2
string or unicode arrays. This is visible for example when trying to
read missing chunks from such arrays. Thanks to Artür Manukyan for
the bug report.

Internal changes

- Some internal changes are preparing the transition to support Zarr
v3:
- "C" and "F" fill orders are now handled via a codec mechanism,
which also supports a wider range of transpose operations.
- The endian configuration is now handled via a codec.
- A GitHub Actions workflow has been added to occasionally test this
package on a big endian platform.
- Bundled libraries have been updated:
- blosc 1.20.1 -> 1.21.6
- snappy 1.1.1 -> 1.2.2
- zstd 1.5.5 -> 1.5.7
- lz4 1.9.2 -> 1.10.0
- Resizable vector in C code for compression now uses the official
exported R C API, instead of internal R functions.
- The const qualifier is used where appropriate in the C code.

[RBedMethyl](/packages/RBedMethyl)
----------

                       Changes in version 0.99.0                        

- Initial Bioconductor submission with core functionality:
- Disk-backed import of modkit bedMethyl files via readBedMethyl()
with HDF5Array storage.
- Field discovery with bedMethylFields().
- Per-site methylation fraction via beta().
- Assay-based filtering with subsetBy() and filterByCoverage().
- Genomic subsetting by interval and GRanges with subsetByRegion()
and [.
- Region-level summaries with summarizeByRegion().
- Coercion to RangedSummarizedExperiment and optional BSseq.
- Vignette and runnable examples.

[Rbwa](/packages/Rbwa)
----

                       Changes in version 1.15.3                        

- Added Config field in DESCRIPTION file for newest BioC
  deployment

                       Changes in version 1.15.1                        

- Changed DFLAGS to LDFLAGS in Makefile.

[RCX](/packages/RCX)
---

                 Changes in version 1.14.2 (2026-03-19)                 

- Fix: Previously, toIgraph required the RCX object to contain an
edges aspect, making it impossible to convert node-only networks to
igraph. The function now handles edgeless networks.

                 Changes in version 1.14.1 (2026-01-20)                 

- Fix: cyGroups follows a different naming structure in the JSON
export. Updated the export function accordingly. Also the CyGroup
ids have to be present in RCX object, a check for those references
has been added.

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.25.1 (2026-01-28)                 

- Adapted code due to changes in Seurat V5

[ReactomePA](/packages/ReactomePA)
----------

                       Changes in version 1.55.1                        

- use 'enrichit' as engine for enrichment analysis (2025-12-07, Sun)

[recount3](/packages/recount3)
--------

                       Changes in version 1.21.1                        

BUG FIXES

- Switched from pryr to lobstr. This was done by @gpertea.

[RedisParam](/packages/RedisParam)
----------

                       Changes in version 1.14.0                        

- (1.13.1) Use 'logger' rather than 'futile.logger'.

[Rega](/packages/Rega)
----

                       Changes in version 0.99.3                        

- Initial public release.
- Provides excel template for filling in submission data
- Implements the core workflow:
- default_parser() – Parses submission data from excel template
- create_client() - Creates API client based on specification
- new_submission() – Submits parsed data to EGA through created
API client
- Provides methods for data validation default_validator()
- Includes vignette for basic usage

[RFGeneRank](/packages/RFGeneRank)
----------

                       Changes in version 0.99.4                        

- Initial Acceptance to Bioconductor

[rGREAT](/packages/rGREAT)
------

                       Changes in version 2.13.1                        

- fix docs

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.56.0                        

Breaking changes

- The default value for the `bit64conversion` argument in
  `H5Aread()` and `H5Dread()` has been set explicitly to `"int"`.
  Relying on a missing value or `"default"` results in a warning
  and will be disallowed in the next release cycle.  Users should
  either specify `bit64conversion="int"` explicitly, or omit this
  argument entirely if they wish to use the provided default.

- The `native` argument to `H5Pcreate()` now triggers an explicit
  deprecation warning and will be removed in the next release
  cycle. This argument was already documented as defunct (without
  warning) in Bioconductor 3.10.

- This package now uses HDF5 version 1.14.6 via Rhdf5lib 1.33.3.

Bug fixes

- The value ordering (row-major vs column-major) is not correct
  when reading unsigned integers with native=TRUE. Thanks to
  Jared Lumpe for the report and detailed reproducible example:
  https://github.com/Huber-group-EMBL/rhdf5/issues/157.

- Minor memory protection issues identified by rchk has been
  resolved.  This might prevent some rare crashes.

- Subsetting and assignment via `[` and `[`<-` operators using a
  variable to define the indices now works as expected. This was
  due to an issue in the timing of evaluation of the index
  arguments. Thanks to Michael Schubert and Malcolm Perry for the
  report and reproducible example:
  https://github.com/Huber-group-EMBL/rhdf5/issues/69.

- Subsetting with `drop=FALSE` no longer errors. Thanks to
  Michael Schubert for the report:
  https://github.com/Huber-group-EMBL/rhdf5/issues/68.

Minor changes

- Documentation for the `name` argument as well as the error
  message when trying to work on a non-existing object now
  include a hint to use `h5ls()` to list the objects in a file.
  Inspired from a suggestion from Frederik Ziebell.

- A warning about void pointer arithmetic in H5Aread_helper_ENUM
  (C code) has been resolved. This increases portability across
  different compilers.

- The H5type argument in h5createAttribute() and
  h5createDataset() now behaves more consistently across both
  functions. Both functions now accept either a HDF5 data type id
  as a string, or the output of H5Tcopy(). Thanks to Luke Zappia
  for the report:
  <https://github.com/Huber-group-EMBL/rhdf5/issues/162>.

- h5version() now also reports the version of the Rhdf5lib
  package.

Internal changes

- Code coverage reports via codecov have been restored.

- The pkgdown documentation website for this package have been
  restored.

- Static analysis via the lintr package is now continuously
  changed on every changes via GitHub Actions.

- Number of skipped tests (when testing, e.g., features only
  available in certain situations) is now explicitly reported.

- GitHub Actions now run on ARM64 versions of all operating
  systems.

- rchk is now part of our continuous integration setup to detect
  potential memory issues in the C code.

[rhdf5filters](/packages/rhdf5filters)
------------

                       Changes in version 1.24.0                        

BUG FIXES

- This package now compiles on -std=c23. This was an issue in the
bundled blosc library. Thanks to Alexander Bontempo for reporting
the issue.
- A R CMD check NOTE about bashims has been resolved.

MINOR CHANGES

- This package no longer bundles the bzip2 compression library, and
uses the system version instead. bzip2 is required by R itself so it
is expected to be available on all systems.

[Rhdf5lib](/packages/Rhdf5lib)
--------

                         Changes in version 2.0                         

Breaking changes

- Build system has been changed from Make to CMake, via the biocmake
  Bioconductor package. Thanks to Aaron Lun (@LTLA) for the initial
  work
  on the migration.

- The HDF5 version bundled by Rhdf5lib has been upgraded from 1.10.7
  to 1.14.6.

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

[RNAshapeQC](/packages/RNAshapeQC)
----------

                       Changes in version 0.99.10                       

- Initial Acceptance to Bioconductor

[rols](/packages/rols)
----

                        Changes in version 3.7.1                        

- Fix bug in rols::Ontologies() (see #48)

                         Changes in version 3.7                         

[SAIGEgds](/packages/SAIGEgds)
--------

                       Changes in version 2.12.0                        

- fix a C++ error that 'span' is ambiguous


[scAnnotatR](/packages/scAnnotatR)
----------

                 Changes in version 1.17.1 (2026-01-15)                 

- Adapted to address changes in Seurat object structure.

[SCArray](/packages/SCArray)
-------

                       Changes in version 1.20.0                        

- fix the package anchors for all Rd \link{}

[SCArray.sat](/packages/SCArray.sat)
-----------

                       Changes in version 1.12.0                        

- update according to Seurat v5

[scConform](/packages/scConform)
---------

                       Changes in version 0.99.1                        

- New package scConform, for uncertainty quantification for cell type
annotation using conformal inference

[scDiagnostics](/packages/scDiagnostics)
-------------

                        Changes in version 1.6.0                        

- Renamed gene shift function for consistency (previously
calculateTopLoadingGeneShifts())
- Added gene specification parameter to calculateGeneShifts()
- Improved calculateGeneShifts() function, plot method, and color
scheme

[scECODA](/packages/scECODA)
-------

                       Changes in version 0.99.7                        

- Adopted SummarizedExperiment data structure

[scider](/packages/scider)
------

                 Changes in version 1.10.0 (2026-04-25)                 

- Update help pages
- New function getHVG() to find top HVGs
- New function plotLISAscatter() to generate a LISA scatter plot
- New functions readVisium(), readVisiumHD(), readXenium(),
readProseg() for reading in data of different types
- Use Matrix package for matrix operations
- Update gridDensity() to handle different types of input
- Speed up hexDensity computation

[scifer](/packages/scifer)
------

                       Changes in version 1.12.1                        

- fix minor bug that crashed quality_report when some sequences had
problematic quality and had NA or Inf values
- this affects summarise_quality and summarise_abi_file

[scLANE](/packages/scLANE)
------

                        Changes in version 1.1.2                        

- Trying to fix Seurat slot vs. layer deprecation error once again.

                        Changes in version 1.1.1                        

- Fixed a deprecation error in geneProgramScoring().

[scLang](/packages/scLang)
------

                 Changes in version 0.99.0 (2025-12-27)                 

- Submitted to Bioconductor

[scoup](/packages/scoup)
-----

                        Changes in version 1.5.1                        

- Extended the vignettes/exData folder to incorporate reviewers
  suggestions.

- Updated the paper, inst and vignettes folders in response to editor's
  queries.

- Fine-tuned paper and necessary files as part of final JOSS
  requirements.

- Updated citation.

[scp](/packages/scp)
---

                        Changes in version 1.21                         

scp 1.21.1

- Update Chris' email address.

scp 1.21.0

- New Bioconductor devel

[scPassport](/packages/scPassport)
----------

                       Changes in version 0.99.0                        

New Features

- Initial Bioconductor submission version.
- seuratPassport(): interactive Shiny gadget to stamp a Seurat object
with a persistent metadata passport stored in @misc$passport.
- read_passport(): prints the full passport and processing log to
console.
- log_step(): appends a processing step entry to
@misc$processing_log.
- Core passport and log logic implemented in C++ via Rcpp for
performance.
- Lineage tracking: parent/child relationships propagated
automatically when parent argument is supplied to seuratPassport().
- Custom fields: any extra key-value metadata can be added through
the
Shiny popup and persists inside the .rds file.

[scran](/packages/scran)
-----

                        Changes in version 1.40                         

- Bugfix for the AUC calculations in pairwiseWilcox() and
  scoreMarkers(), to avoid integer overflow with large numbers of
  ties.

- Deprecate many historical functions in favor of their
  replacements in the scrapper package.

[scrapper](/packages/scrapper)
--------

                        Changes in version 1.6.0                        

- Multiple changes to clusterGraph():
  
  • Updates to use the latest version of the igraph library.
  
  • method="leiden" now supports the ER objective function in
  leiden.objective=.
  
  • Non-zero status are now directly reported as errors, status
  is removed from the output.

- Multiple additions to scoreMarkers():
  
  • Support use of quantiles to average statistics across
  blocks, via the new block.average.policy= and
  block.quantile= arguments.
  
  • Allow calculation of arbitrary quantiles as a summary
  statistic for the effect sizes of each group, via the new
  compute.summary.quantiles= argument.
  
  • Output now contains the number of rows, row names and group
  IDs for easier interpretation by downstream functions.
  
  • Added top.index.only= argument to report only the indices
  of the top genes when all.pairwise= is an integer.

- Added arguments to summarizeEffects() to disable calculation of
  certain summary statistics, e.g., by setting
  compute.summary.median=FALSE. Also support calculation of
  arbitrary quantiles as described for scoreMarkers().

- suggest*QcThresholds() functions now return a block.ids field
  that preserves the type of block=. This enables easier matching
  in the corresponding filter*QcMetrics() function.

- Multiple changes to modelGeneVariances():
  
  • Allow the use of quantiles to average statistics across
  blocks via the new block.average.policy= and
  block.quantile= arguments.
  
  • Return a block.ids field that preserves the type of block=
  for easier matching.

- Multiple changes to runPca():
  
  • Added a subset= argument that performs a PCA on a subset of
  features but still reports the full rotation matrix for all
  features.
  
  • Return a block.ids field that preserves the type of block=
  for easier matching.

- clusterKmeans() now emits a warning upon convergence failure.
  This can be disabled by setting warn=FALSE.

- Simplified the automatic naming performed by combineFactors()
  when factors is unnamed.

- Added SummarizedExperiment-compatible wrappers for all
  functions. These can be identified by their *.se suffix, e.g.,
  aggregateAcrossCells.se(). Some functions execute multiple
  steps for convenience, e.g., quickRnaQc.se() runs all of the
  RNA-related QC functions.

- Added analyze.se() to replace the analyze() function. The
  former is more convenient it stores most outputs in a
  SingleCellExperiment and formats the marker tables correctly.
  The latter is now deprecated.

- compute*QcMetrics() functions now return a DataFrame for more
  convenient inspection and manipulation.

- All functions that previously returned a base data.frame will
  now return a DataFrame from S4Vectors. This is a little
  prettier and we already load S4Vectors anyway.

- Removed checks on num.neighbors= when pre-computed neighbor
  search results are supplied in runUmap(), buildSnnGraph() and
  subsampleByNeighbors(). Any differences between num.neighbors=
  and the number of neighbors in the pre-computed results are now
  ignored.

- Implemented the LogNormalizedMatrix class for delayed
  log-transformation. Modified normalizeCounts() to return
  instances of this class when log=TRUE and delayed=TRUE. This
  allows downstream C++ code to use a more efficient
  normalization wrapper via initializeCpp(). Note that there are
  some very minor differences in the log-normalized values when
  computed in R via extract_array() and those in C++ via tatami
  extraction.

- Added the countGroupsByBlock() function to count the number of
  cells within each combination of group and block. This is
  primarily intended for diagnosing batch effects.

- Added options to control which statistics are computed in
  aggregateAcrossCells(). Also added bindings to compute the
  per-group medians via the compute.median= option.

[scRepertoire](/packages/scRepertoire)
------------

                        Changes in version 2.7.3                        

BUG FIXES

- Fixed combineBCR() assigning the same group.by value to all cells
instead of per-barcode values.
- Fixed clonalCluster() failing when group.by produces a single group
by correcting the condition for adding the group column to the edge
list.
- Fixed test for CTstrict clustering pattern to match the new default
chain = "IGH" behavior in combineBCR().

                        Changes in version 2.7.2                        

NEW FEATURES

- New clonalBin() function to bin clones by frequency or proportion
without requiring a single-cell object. Adds clonalFrequency,
clonalProportion, and cloneSize columns to the output of
combineTCR(), combineBCR(), or combineExpression(). Supports custom
bin thresholds, optional grouping by metadata variable, and chain
filtering.
- New vizCirclize() function for quick chord diagram visualization
without manual circlize code. Supports directional arrows, custom
colors, and sector annotations.
- getCirclize() major enhancements:
- Multi-level hierarchical grouping via group.by accepting a
vector of columns
- New method parameter: "unique", "abundance", "jaccard",
"overlap"
- New symmetric parameter for directional flow analysis
- New include.metadata parameter returning sector statistics
- New filtering options: min.shared, top.links, filter.sectors
- Built-in color palette generation with palette parameter
- alluvialClones() major enhancements:
- New top.clones, min.freq, highlight.clones, and highlight.color
parameters
- Visual customization: stratum.width, flow.alpha, show.labels,
label.size
- New order.strata parameter for controlling level ordering within
each stratum
- Enhanced export.table output now includes freq, prop, and rank
columns
- combineBCR() defaults the clustering call to "IGH" instead of
"both"

API CHANGES

- Soft-deprecated camelCase arguments across all exported functions
in
favor of dot.notation (e.g., cloneCall to clone.call, exportTable to
export.table, cloneSize to clone.size, filterNA to filter.na,
addLabel to add.label, clonalSplit to clonal.split). All deprecated
arguments will continue to work with a deprecation warning until
version 3.0.0.

BUG FIXES

- Fixed combineExpression() failing with "undefined columns selected"
when input data already contained clonalFrequency/clonalProportion
columns from a prior clonalBin() call.

                        Changes in version 2.6.2                        

UNDERLYING CHANGES

- Expanded functionality in combineBCR() and clonalCluster():
- New metrics beyond normalized Levenshtein edit distances
- Support for raw and normalized-based calculations
- Support for distance matrices to allow for alignment
- Added support for declaring chains = "IGL", "IGK", or "Light" to
get
all light chains in downstream quantification.

BUG FIXES

- Fixed handling of multiple chains in combineBCR(), specifically in
formatting CTstrict.

                        Changes in version 2.6.1                        

Update to match Bioconductor Release 3.22 on 2025/10/30.

BUG FIXES

- Fixed order.by issue in positionalProperty().
- Fixed individual chain call for combineExpression().
- Fixed issue with removing kmer with ";" in percentKmer().

[scRNAseqApp](/packages/scRNAseqApp)
-----------

                       Changes in version 1.11.26                       

- Add datasets available date.

                       Changes in version 1.11.25                       

- Add search link for species.

                       Changes in version 1.11.24                       

- Add 'download' button for webstats.

                       Changes in version 1.11.23                       

- Add 'add' and 'delete' method for lasso.

                       Changes in version 1.11.22                       

- Add desity plot for sunburst and deconvolution.

                       Changes in version 1.11.21                       

- Fix the issue for idconv.

                       Changes in version 1.11.20                       

- Add cell deconvolution plot.

                       Changes in version 1.11.19                       

- Fix a bug in molecule plot for discrete_scale.

                       Changes in version 1.11.17                       

- Using coord_cartesian but not xlim/ylim.

                       Changes in version 1.11.16                       

- Add molecule module.

                       Changes in version 1.11.15                       

- Fix the bugs in database links.

                       Changes in version 1.11.14                       

- Add pan controler and resizable divider.

                       Changes in version 1.11.13                       

- fix a bug for double click events in the user mode.

                       Changes in version 1.11.12                       

- save cell segmentation for Seurat object.

                       Changes in version 1.11.11                       

- add help function addCellBorders, addCellLinks and
  addBackgroundImage.

                       Changes in version 1.11.10                       

- add background plot.

                       Changes in version 1.11.9                        

- add cell segmentaion plot.

                       Changes in version 1.11.8                        

- remove the force gene symbol table update.

                       Changes in version 1.11.7                        

- fix the bug when subset the cells by percentage.

                       Changes in version 1.11.6                        

- Add xlim and ylim for cellinfo and gene exprs plot give the
  possibility
  of cross dataset comparison.

                       Changes in version 1.11.5                        

- replace 'size' by 'linewidth' for element_line in theme.

                       Changes in version 1.11.4                        

- Keep the height controler for interactive plot in explorer.

- Add perentage filter for cell info plot.

                       Changes in version 1.11.3                        

- Add apply lasso selection for explorer.

                       Changes in version 1.11.2                        

- Add group B for reduction in explorer.

                       Changes in version 1.11.1                        

- Add download function for lasso.

[scToppR](/packages/scToppR)
-------

                       Changes in version 0.99.5                        

- Initial Acceptance to Bioconductor

[scTypeEval](/packages/scTypeEval)
----------


                Changes in version 0.99.21 (2024-01-27)                 

New Features

- Initial Bioconductor submission with ground-truth-agnostic cell
type
evaluation framework
- Support for multiple input formats: Seurat, SingleCellExperiment,
and count matrices
- Multiple dissimilarity methods:
- Pseudobulk-based distances (Euclidean, Cosine, Pearson)
- Wasserstein distance on single-cell distributions
- Reciprocal classification approaches
- Comprehensive internal validation metrics:
- Silhouette scores
- Neighborhood Purity
- Ward's clustering consistency
- Orbital medoid distances
- Average similarity measures
- 2-label silhouette analysis
- Gene selection methods:
- Highly variable genes (HVG) detection
- Cell-type-specific marker identification
- Custom gene list support
- Dimensional reduction support (PCA, pre-computed embeddings)
- Cross-sample and cross-study benchmarking capabilities
- Customizable visualization tools (heatmaps, MDS, PCA plots)
- Hierarchical clustering analysis

Improvements

- Comprehensive documentation with vignettes
- Examples for all 19 exported functions
- Support for optional dependencies (transformGamPoi, glmGamPoi)
- Efficient euclidean distance computation with custom C
implementation
- Parallel processing support via BiocParallel

Bug Fixes

- Proper handling of sparse matrix formats
- Robust batch effect handling
- Consistent results across different label granularities

Documentation

- Main vignette: comprehensive tutorial with real-world examples
- Quick start guide: minimal workflow for rapid evaluation
- Full API documentation for all exported functions

[scuttle](/packages/scuttle)
-------

                        Changes in version 1.22                         

- Improved the efficiency of the downsampling algorithm in
  downsampleMatrix(), which will change the results slightly. In
  addition, new options are added to control the output type
  (integer/double) and class (dense matrix or SVT_SparseMatrix).

- Removed "median" from the default statistics= in
  summarizeAssayByGroup(). This enables default use of the fast
  code path for improved performance.

- Deprecate all functions that have more efficient counterparts
  in scrapper.

[SEMPLR](/packages/SEMPLR)
------

                       Changes in version 0.99.0                        

- Ready for Bioconductor

- Your bug fixes. See more details at
http://bioconductor.org/developers/package-guidelines/#news.

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.52.0                        

NEW FEATURES

- new 'header' in `seqGDS2VCF()` to specify whether exporting the
  metadata header to a VCF file or not

- set 'nosample=TRUE' in `seqGDS2VCF()` to exclude the sample-level
  data

- new 'verbose.progress' in `seqGDS2VCF()` to control the progress
  display

- `seqGet2bGeno()` allows a GDS file name in the first argument

- new 'parallel' in `seqGet2bGeno()` for parallel loading

- new 'SeqArray:::process_block_index' and
  'SeqArray:::process_block_count'
  used in children processes when balancing work

- new 'ploidy' in `seqVCF2GDS()` and `seqBCF2GDS()`

UTILITIES

- Reduce memory usage in `seqParallel()` by avoiding the transfer of
  unused
  data during work balancing in parallel

- minor fix in `seqAsVCF()` when loading VariantAnnotation

- fix the package anchors for all Rd \link{}

- Updated `seqSetFilterPos()` is significantly faster than the old
  version

BUG FIXES

- fix `seqApply(..., margin="by.sample")` dimension mismatch error when
  multiallelic variants are present

                       Changes in version 1.50.1                        

UTILITIES

- Reduce memory usage in `seqMerge()` when there are many INFO
  variables

[Seqtometry](/packages/Seqtometry)
----------

- Initial Acceptance to Bioconductor

[sfi](/packages/sfi)
---

                       Changes in version 0.99.5                        

- Initial Acceptance to Bioconductor

[SingleCellExperiment](/packages/SingleCellExperiment)
--------------------

                       Changes in version 1.34.0                        

- counts(), logcounts() and other named getters/setters will now
  raise a warning if i= is specified, rather than silently
  treating it as the next argument to assay() (usually
  withDimnames=).

[singleCellTK](/packages/singleCellTK)
------------

                       Changes in version 2.21.1                        

- Updated depreciated parameters in Seurat function calls

[SingleR](/packages/SingleR)
-------

                       Changes in version 2.14.0                        

- Added hint.sce= to trainSingleR() to nudge users towards
  setting a more appropriate de.method= for single-cell
  references. Similarly, suggest the use of aggr.ref=TRUE for
  large single-cell references.

- Added a warning when names of results= and trained= are
  different in combineRecomputedResults().

- Updated to the latest version of singlepp, which provides
  direct support for sparse test/reference matrices. This also
  eliminates the need for BiocNeighbors so all BNPARAM= arguments
  are now soft-deprecated.

- Added score.args= to configureMarkerHeatmap() to pass along to
  scoreMarkers(), mainly for block= and threshold= options.

- Added center= and average= options to plotMarkerHeatmap(), for
  centering of genes and averaging of labels, respectively. Also
  pass along score.args=.

- Support block= in getClassicMarkers() to define a blocking
  factor within a single reference matrix. This serves as an
  alternative to supplying multiple reference matrices.

[smartid](/packages/smartid)
-------

                        Changes in version 1.7.3                        

- Major memory and performance optimization for cal_score() and
top_markers(). On a 20,000 gene x 100,000 cell sparse input peak
memory drops from roughly 100 GB to a few GB, and top_markers() with
the default gaussian() family runs in seconds instead of hours.
- Breaking change: cal_score() no longer stores the intermediate tf,
idf and iae matrices in metadata() by default. Callers that relied
on metadata(se)$tf / $idf / $iae (introduced in v1.1.1) must now
pass return.intermediate = TRUE. When the flag is TRUE, the stored
idf / iae for labelled methods (prob, rf) are now compact G x K
matrices (columns = unique labels); expand with md$idf&#91;,
as.character(label)&#93; to recover the legacy per-cell form.
- Internal refactor: labelled idf_prob, idf_rf, iae_prob and iae_rf
helpers now return a compact G x K matrix. cal_score() composes the
final score through per-group column-block multiplication, avoiding
the materialisation of full G x N intermediates.
- cal_score() no longer forces dense conversion of dgCMatrix inputs;
the score assay stays sparse throughout the pipeline when the input
is sparse.
- tf(), idf_hdb(), iae_hdb() and all IAE helpers now preserve
dgCMatrix sparsity by routing column scaling through
Matrix::Diagonal and replacing the densifying x&#91;x < 0&#93; <- 0 pattern
with pmax0_offset().
- top_markers_glm() has a vectorised closed-form least-squares fast
path for the default gaussian() + identity link; non-gaussian
families or rank-deficient designs automatically fall back to the
legacy per-gene glm() loop with no behaviour change.
- top_markers_abs() aggregates directly on the scored matrix via
sparseMatrixStats::rowMeans2 / rowMedians / rowMads, removing the
intermediate wide data.frame that previously reached tens of GB.
- scale_mgm() caches per-group column indices and collapses the
two-step (expr - mgm) / (sds + 1e-8) into a single broadcast.
- The multi = TRUE branch of the labelled IDF/IAE helpers switched
from an O(G * K^2) apply() to an O(G * K) top-1 + top-2 trick via
the new rowwise_notin_max() helper.
- New inst/bench/benchmark_smartid.R micro-benchmark script; new
tests/testthat/test-numerical-equivalence.R pins cal_score() and
top_markers() outputs to a frozen pre-refactor snapshot at 1e-10
tolerance.
- No new dependencies; the refactor relies entirely on Matrix,
sparseMatrixStats and base R.

                        Changes in version 1.7.2                        

- Fix dplyr defunct.

                        Changes in version 1.7.1                        

- Add bioRxiv citation.

[smoppix](/packages/smoppix)
-------

                        Changes in version 1.2.3                        

- Switch to mgcv::scasm for monotonically decreasing p-splines

                        Changes in version 1.2.2                        

- Adapt code for negative hypergeometric from the orphaned extraDistr
package and take it in-house

                        Changes in version 1.2.1                        

- Anticipate coming changes to lme4

[SMTrackR](/packages/SMTrackR)
--------

                       Changes in version 0.99.8                        

- Initial Acceptance to Bioconductor

[SNPRelate](/packages/SNPRelate)
---------

                       Changes in version 1.46.0                        

UTILITIES

- fix the package anchors for all Rd \link{}

[sosta](/packages/sosta)
-----

                        Changes in version 1.3.4                        

- removed non working GitHub actions

                        Changes in version 1.3.3                        

- clarified difference between shapeMetric and totalShapeMetrics.
- updated FOV border characterization in vignette.

                        Changes in version 1.3.2                        

- changed minBoundaryDistances to allow for distance to FOV border
characterization
- added "Structure boundary vs FOV boundary" section to vignette.

                        Changes in version 1.3.1                        

- added argument to simulate noise in createPointPatternTissue
- improved intensity value threshold estimation for reconstruction in
.intensityThreshold

[SpaceTrooper](/packages/SpaceTrooper)
------------

                        Changes in version 1.1.7                        

- fixing author name typo

                        Changes in version 1.1.6                        

- adding new vignette about SpaceTrooper utilities

                        Changes in version 1.1.5                        

- minor enhancement of cosmx input reading, polygonsCol added
- minor documentation clarification on QCScore computation
- cleaned LICENSE stub and creating LICENSE.md file

                        Changes in version 1.1.4                        

- updating documentation along multiple functions
- adding graphical abstract
- updating README
- updating Vignettes and providing better data examples

                        Changes in version 1.1.3                        

- fixing verbose message in QC functions.

                        Changes in version 1.1.2                        

- adding volume instead of area for MERFISH.
- refactoring log2CountArea to log2DensitySignal
- added size and alpha arguments to plotCellsFovs
- added scaleBar argument to plotCellsFovs, plotCentroids,
plotPolygons and plotZoomFovsMap

                        Changes in version 1.0.1                        

- fixing minor bugs in AspectRatio internal computation for
technology
missing it.
- fixing merfish reading where colData where not properly sorted and
sync with assay cells.

[SpatialArtifacts](/packages/SpatialArtifacts)
----------------

                       Changes in version 0.99.10                       

- Initial Acceptance to Bioconductor

[SpatialDecon](/packages/SpatialDecon)
------------

                 Changes in version 1.20.1 (2025-12-22)                 

- Bug fix in vignette Seurat creation

[spatialFDA](/packages/spatialFDA)
----------

                        Changes in version 1.3.3                        

- rewrote summary.functionalGam and added a new coef.functionalGam to
have no breaking changes in spatialFDA due to changes in
refund::pffr naming.
- this makes spatialFDA again compatible with the OSTA chapter.

                        Changes in version 1.3.2                        

- covariance matrices are now corrected for cluster-correlations
along
the domain $r$ with sandwich-estimators by default - thank you
@fabianscheipl for implementing this in refund::pffr
- models can now also be fit with gamm4
- Experimental: implemented an AR(1) penalty along the domain $r$ in
bam fits. This is to some extent already adressed by the
cluster-robust sandwich estimators so might be removed in the
future.

                        Changes in version 1.3.1                        

- option to square-root transform the weights in the functional Gam
- changes to cutoff calculations at the saturating end of $G$
functions
- fix as.formula calls by passing an empty.env call
- cleaned up spatialInference from exploratory code.


[Spectra](/packages/Spectra)
-------

                        Changes in version 1.21                         

Change 1.21.7

- Add parameter direction to shiftPeaks() allowing to define whether
peaks should be shifted left or right (default).

Change 1.21.6

- Add shiftPeaks() function.
- Add parameter onlyCore = FALSE to dropNaSpectraVariables().

Change 1.21.5

- Fix potential issue by introducing concatenation of data.frames
with
data.table::rbindlist() in version 1.21.3: new function
rbindlistWithRownames() added which uses rbindlist() but preserves
also the row names of the data.frames.
- Fix unit tests for joinPeaksGnps() to follow recent updates in
MsCoreUtils.

Change 1.21.4

- Refactor compareSpectra() to support spectra similarity functions
returning the similarity score and the number of peak pairs on which
the score was calculated. This allows to correctly report the number
of peak pairs used by e.g. the modified cosine (GNPS) similarity
score. See also issue #350.

Change 1.21.3

- Use data.table::rbindlist() for merging of MsBackendMemory
instances
and for backendInitialize() for MsBackendMzR. This adds data.table
as a dependency but improves performance for the above mentioned
functionality.
- Fix bug in fragmentGroupIndex() that would return the indices in a
wrong order.

Change 1.21.2

- Replace msdata with MsDataHub package.

Change 1.21.1

- Fix precursorPurity() for empty spectra.

[SpectraQL](/packages/SpectraQL)
---------

                         Changes in version 1.5                         

Changes in 1.5.1

- Load test and example files from MsDataHub dropping the dependency
from msdata.

[SpectriPy](/packages/SpectriPy)
---------

                         Changes in version 1.1                         

Changes in 1.1.8

- Implement pyspec_copy_on_replace() to enable copying/cloning the MS
data in Python for any data replacement operation performed in R
(issue #91).

Changes in 1.1.7

- Improved $<- operation: only replace the values for the selected
spectra variable.
- Improved peaksData<-, mz<- and intensity<- implementations: only
replace the respective data, but not the full spectra data
(including metadata).
- Improved lengths() implementation: retrieves the number of peaks
directly in Python.
- Documentation updates.

Changes in 1.1.6

- Automatic renaming and remapping of upper case or camelCase spectra
variables to snake_case metadata fields in backendInitialize(),
spectraData<- and $<-.
- Full support of the Spectra test suite.
- Add spectraNames<- method for MsBackendPy and support
getting/setting spectrum names. They are stored in a metadata field
(spectra variable) spectrum_name.

Changes in 1.1.5

- Fix $<-, spectraData()<- and related replacement methods for
Spectra/MsBackendPy of length 1.

Changes in 1.1.4

- Implement a setBackend() method for MsBackendPy with a parameter
applyProcessing to allow applying the processing queue and storing
the modified peaks data to Python.

Changes in 1.1.3

- Introduce new ModifiedCosineHungarian() and ModifiedCosineGreedy()
parameters to match changes in matchms version 0.32.

Changes in 1.1.2

- Use matchms 0.31.

[SpiecEasi](/packages/SpiecEasi)
---------

                       Changes in version 1.99.5                        

BUG FIXES

- Fixed vignette rebuild error with pulsar >= 0.3.13: resolve the
estimation function before passing to pulsar() - this avoids lazy
match.fun evaluation that fails on PSOCK cluster workers.

INFRASTRUCTURE

- Removed obsolete CXX_STD = CXX11 from src/Makevars (R >= 4.6
defaults to C++20)
- Fixed .Rbuildignore to exclude leftover vignettes/figure directory
from knitr

                       Changes in version 1.99.4                        

BUG FIXES

- Fixed vignette build failure caused by makeCluster(4, type =
"SOCK")
requiring unavailable snow package
- Replaced broken snow cluster example with reference to batch.pulsar
with conffile='snow'
- Set batch mode vignette chunks to eval=FALSE to avoid build
failures
on systems without batchtools infrastructure

                       Changes in version 1.99.3                        

BUG FIXES

- Fixed BugReports URL in DESCRIPTION (removed double slash)
- Replaced sapply with vapply in R/fitdistr.R for improved type
stability
- Fixed incomplete final line in .Rbuildignore file

IMPROVEMENTS

- Added BiocManager installation instructions to main vignette
- Added runnable examples to man pages for getOptInd, robustPCA,
stars.roc, and triu functions
- Cleaned up promotional language in NEWS entries
- Removed unused inst/extdata files to reduce package size

INFRASTRUCTURE

- Added .registration = TRUE to useDynLib in NAMESPACE
- Updated package version to 0.99.2 for Bioconductor submission
- Verified vignette runtime compliance with Bioconductor guidelines

                       Changes in version 1.99.2                        

NEW FEATURES

- Added comprehensive Bioconductor package infrastructure
- Added automated GitHub Actions workflow for Bioconductor continuous
integration
- Added test coverage reporting and automated testing
- Added Bioconductor-style vignettes and documentation
- Added biocViews categorization: Software, Microbiome,
Metagenomics, GraphAndNetwork, NetworkInference
- Added comprehensive support documentation and issue templates
- Implemented proper code of conduct and contribution guidelines

INFRASTRUCTURE IMPROVEMENTS

- Added Bioconductor-specific GitHub Actions workflow
(check-bioc.yml)
- Implemented automated R CMD check with Bioconductor standards
- Added BiocCheck integration for package validation
- Enhanced test coverage with automated reporting
- Added lifecycle badge indicating stable package status

BUG FIXES

- Resolved compliance issues identified during BiocCheck
- Enhanced error handling and package robustness

[splatter](/packages/splatter)
--------

                 Changes in version 1.36.0 (2026-04-29)                 

- Replace deprecated scuttle functions with scrapper

- Deprecate mfa simulation functions as the mfa package is
  deprecated

- Adjust help messages in warnings

- Adjust failing lun2Estimate tests

[SpliceImpactR](/packages/SpliceImpactR)
-------------

                       Changes in version 0.99.4                        

- Initial Acceptance to Bioconductor

[splicelogic](/packages/splicelogic)
-----------

                       Changes in version 0.99.0                        

Submit to Bioconductor 3.23

[SpNeigh](/packages/SpNeigh)
-------


                       Changes in version 0.99.0                        

Initial Bioconductor submission

New Features

- Introduced the RunSpatialDE() function for spatial differential
expression using spatial weights.
- Added ComputeSpatialEnrichmentIndex() for enrichment
quantification.
- Implemented GetBoundary(), GetRingRegion(), and RemoveOutliers()
for
spatial region detection.
- Added PlotSpatialExpression() and other visualization tools.

[StatescopeR](/packages/StatescopeR)
-----------

                      Changes in version 26-7-2025                      

- Initial Acceptance to Bioconductor 

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

                       Changes in version 1.23.2                        

- remove stato classes; these were deprecated since 1.5.2 and are no
  longer used.

- ontology is still available via ontology slots.

- rols package is deprecated; using internal functions instead.

- replaced ontology cache with online/offline mode

- added ontology output for model sequences

[structToolbox](/packages/structToolbox)
-------------

                       Changes in version 1.23.2                        

- remove deprecated stato classes; use ontology slots instead

                       Changes in version 1.23.1                        

- general maintenance

- fix PLSDA probability calculation

[SVP](/packages/SVP)
---

                        Changes in version 1.3.1                        

- Bioconductor 3.22 release, and Bioconductor 3.23 (devel) bump.
(2025-10-31, Fri)

[TargetSearch](/packages/TargetSearch)
------------

                       Changes in version 2.14.0                        

DOCUMENTATION

- Fix typo in `ri_plot_peak` man page.

[TaxSEA](/packages/TaxSEA)
------

                        Changes in version 1.3.3                        

- Add ssTaxSEA() function for single sample enrichment testing and
documentation

                        Changes in version 1.3.2                        

- Add taxon_rank_sets() helper with documentation and tests
- Add functionality to interact with tse objects

[tidyexposomics](/packages/tidyexposomics)
--------------

                       Changes in version 0.99.16                       

Initial Bioconductor submission.

- Implements a full exposome-omics pipeline including:
- Quality control
- Association analysis
- Differential abundance analysis
- Multi-omics integration analysis
- Functional enrichment analysis
- Vizualization

[TileDBArray](/packages/TileDBArray)
-----------

                       Changes in version 1.22.0                        

- Workaround for TileDBRealizationSink to properly write
  higher-dimensional dense arrays.

[topdownr](/packages/topdownr)
--------

                        Changes in version 1.33                         

Changes in version 1.33.1

- Adapt to changes in ggplot2.
- Adapt .calculateFragments to changes in PSMatch::calculateFragments
introduced in PSMatch:PR41.

[toppgene](/packages/toppgene)
--------

                 Changes in version 0.99.2 (2026-02-16)                 

- Initial Acceptance to Bioconductor

[TPP](/packages/TPP)
---

                       Changes in version 3.39.1                        

- Avoid deprecated dplyr::tbl_df() calls in
  biobroom::tidy.ExpressionSet() by changing output to data.frame for
  compatibility with dplyr ≥ 1.0.0.

- Deprecated parallelization and nCores parameter in fitting and
  plotting functions to improve maintainability and robustness to
  dependency changes.

- Fixed ggplot2 deprecation warning to maintain compatibility with
  recent ggplot2 versions

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.47.1                        

- Fix the mis-match y-axis for dandelion.plot when max score is smaller
  than 1.

[transcriptR](/packages/transcriptR)
-----------

                       Changes in version 1.39.4                        

- Fixed a failing unit test caused by an update to the
  TxDb.Hsapiens.UCSC.hg19.knownGene library. Replaced the external
  dependency with the internal annot dataset to ensure stability.

- Moved e1071 from Imports to Suggests

- Updated the roxygen2 version and rebuilt all manual pages.

- Fixed incorrect roxygen documentation for internal datasets.


[TSAR](/packages/TSAR)
----

                        Changes in version 1.9.3                        

- Introduction of cubic spline with beta-knots for smoother curve
fitting, completed on December, 2025.

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.30.0                        

- Updates to vignette to mention edgeR's new function
DGEListFromTximeta.
- Fix the 'call dbDisconnect() when finished' message by explicitly
doing it on exit. Thanks to Gordon Smyth for noting this.

                       Changes in version 1.29.5                        

- Working on updateMetadata, gene_id had issues as CharacterList
type.
Now downgrading to character.

                       Changes in version 1.29.4                        

- Fixing sapply bug in inspectDigests and updateMetadata.

                       Changes in version 1.29.3                        

- Fixing github actions workflow re quarto issue.

[tximport](/packages/tximport)
--------

                       Changes in version 1.40.0                        

- Updates from Gordon Smyth and Pedro Baldoni
  to the vignette using new functions in edgeR.


[UPDhmm](/packages/UPDhmm)
------

                        Changes in version 1.7.1                        

NEW FEATURES

- Added a long-read processing example in a new vignette.

IMPROVEMENTS

- Refined recurrent annotation in markRecurrentRegions(). Recurrent
status is now determined using a minimum overlap fraction between
events and recurrent regions, preventing marginal overlaps from
being annotated as recurrent.

- More precise recurrent region definition in
identifyRecurrentRegions(). Replaced direct interval merging with an
overlap-aware, chromosome-wise agglomerative hierarchical clustering
approach. Recurrent regions are now defined based on overlap
similarity between events, with configurable clustering stringency.

BUG FIXES

- Corrected collapsed event size calculation in collapseEvents(). The
total size of collapsed events is now computed as the true genomic
span of the region rather than the sum of individual event sizes.

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


[velociraptor](/packages/velociraptor)
------------

                       Changes in version 1.21.3                        

- Update dependencies for MacOSXArm (again).

                       Changes in version 1.21.2                        

- Update dependencies for MacOSXArm.

                       Changes in version 1.21.1                        

- Pin version of all dependencies.

[VisiumIO](/packages/VisiumIO)
--------

                        Changes in version 1.8.0                        

Bug fixes and minor improvements

- Fixed issue when using sample_id parameter in TENxVisiumHD
(@michaelplynch, #17)
- Allow spacerangerOut inputs in TENxVisium to have an outs folder
and
optionally check for the spatial subfolder within it (@estellad,
#18).

[VISTA](/packages/VISTA)
-----

- Initial Acceptance to Bioconductor

[vsn](/packages/vsn)
---

                       Changes in version 3.79.1                        

- The deprecated `ggplot2::aes_string()` used in `meanSdPlot()` has
  been replaced with `ggplot2::aes()`.

- The error message when `reference@mu` doesn't have the expected
  length in `vsnMatrix()` no longer says
  "too few arguments" and now reports the actual and expected lengths.

[wavClusteR](/packages/wavClusteR)
----------

                       Changes in version 2.45.1                        

- Fixed imports following migration of `supportedUCSCtables()` and
  `makeTxDbFromUCSC()` to {txdbmaker}

[wavFeatExt](/packages/wavFeatExt)
----------

                       Changes in version 0.99.14                       

- Initial Acceptance to Bioconductor

[XAItest](/packages/XAItest)
-------

                        Changes in version 1.3.1                        

- Add package citation metadata for the Bioconductor landing page.
- Document exported helper functions and sync the XAI.test() manual.
- Import stats::runif in genScenario() to address package check
notes.

- Remove obsolete vignette; bump to 1.1.1

[xcms](/packages/xcms)
----

                         Changes in version 4.9                         

Changes in version 4.9.4

- Improve performance for some functions by using
Spectra::rbindlistWithRownames() to merge data.frames.
- Replace the msdata package for test data files with MsDatahub.

Changes in version 4.9.3

- Fix title of plotPrecursorIons() plot.
- Add parameter rtCenterFun to PeakDensityParam() to support
specifying the function to calculate the reported retention time of
a feature ("rtmed").

Changes in version 4.9.2

- Fix issue with refineChromPeaks() and MergeNeighboringPeaksParam
where in certain cases completely overlapping peaks were not merged
(issue #825).

Changes in version 4.9.1

- Add the citation and reference for the new 2025 xcms paper

[zellkonverter](/packages/zellkonverter)
-------------

                       Changes in version 1.22.0                        

Minor changes

- 
  Minor fixes to tests and documentation


NEWS from existing Data Experiment Packages
===================================


[curatedCRCData](/packages/curatedCRCData)
--------------

                       Changes in version 2.43.1                        

- Updated package maintainer to Levi Waldron
  (lwaldron.research@gmail.com).

- Refactored DESCRIPTION to use standard Authors@R fields, including
  maintainer ORCID and funding info (NCI 5U24CA180996).

- Fixed invalid biocViews categories to strictly use acceptable
  ExperimentData terms.

- Replaced the unnecessary affy dependency with Biobase in the package
  Suggests and the vignette.

- Added a new CITATION file pointing to the 2018 Genome Biology
  meta-analysis paper (PMID: 30253799).

- Initialized Git LFS to properly track .rda and .RData files.

- Removed deprecated zzz.R warning.

- Converted the vignette to RMarkdown (.Rmd).

[dominatRData](/packages/dominatRData)
------------

                       Changes in version 0.99.1                        

- Modification of the Description

- Modification of the vignettes!

                       Changes in version 0.99.0                        

- Initial Bioconductor submission.

[DoReMiTra](/packages/DoReMiTra)
---------

                        Changes in version 1.2.0                        

- adding a new function to list metadata values

- filtering by author name

- adding Tissue field in the metadata

[EMTscoreData](/packages/EMTscoreData)
------------

                       Changes in version 0.99.0                        

- submitted to bioconductor

[gDRtestData](/packages/gDRtestData)
-----------

              Changes in version 22023-05-15 (2023-05-15)               

- fix related with data.table

               Changes in version 2026-04-13 (2026-04-13)               

- migrate test data files from .qs to .qs2 format

               Changes in version 2025-10-30 (2025-10-30)               

- synchronize Bioconductor and GitHub versioning

               Changes in version 2025-08-12 (2025-08-12)               

- fix usage of ifelse

               Changes in version 2025-07-29 (2025-07-29)               

- update test data

               Changes in version 2025-04-16 (2025-04-16)               


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


[msdata](/packages/msdata)
------

                       Changes in version 0.51.2                        

- remove
  TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML.gz,
  available in MsDataHub.

- remove
  TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz.

- remove PestMix1_DDA.mzML, available in MsDataHub.

- remove PestMix1_SWATH.mzML, available in MsDataHub.

                       Changes in version 0.51.1                        

- Add deprecation startup message

[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.23.2                        

NEW FEATURES

- gene_set_enrichment() now returns a $GeneList as implemented by
  @lahuuki. See details at
  https://github.com/LieberInstitute/spatialLIBD/pull/119.

                       Changes in version 1.23.1                        

NEW FEATURES

- fetch_data() now provides access to all LFF_spatial_ERC project
  files associated with the preprint at
  https://doi.org/10.1101/2025.11.20.689483. @lahuuki added this
  information as part of
  https://github.com/LieberInstitute/spatialLIBD/pull/118.


Deprecated and Defunct Packages
===============================

**SOFTWARE:**

Thirty seven software packages were removed from this release (after being deprecated
in Bioc 3.22):

- bgx, BiGGR, biodbHmdb, biodbNcbi, biodbNci, biodbUniprot, ccmap, CellScore, CINdex, cisPath, DEP, gpuMagic, Harshlight, hiAnnotator, hiReadsProcessor, interactiveDisplay, interactiveDisplayBase, lapmix, LinTInd, lute, MADSEQ, oppti, PhenStat, qckitfastq, ReactomeGraph4R, Repitools, rGADEM, rRDP, SARC, seqArchR, seqArchRplus, seqTools, Streamer, TitanCNA, TransView, traviz, XNAString

Previously announced deprecated packages were fixed and restored to Bioconductor 3.23: 

- hypeR, netZooR, Rfastp

Fifty two software packages are deprecated in this release and will be removed in Bioc 3.24:

- APAlyzer, ballgown, bamsignals, barcodetrackR, basecallQC, biobroom, biodbChebi, BiRewire, BPRMeth, BubbleTree, ccrepe, CelliD, ChIPQC, cummeRbund, CuratedAtlasQueryR, debCAM, DeconRNASeq, geneXtendeR, GEOexplorer, GNET2, granulator, hca, hmdbQuery, IMAS, IONiseR, Melissa, MetaNeighbor, MethReg, mfa, microSTASIS, MineICA, motifcounter, MSPrep, nearBynding, netprioR, normr, Organism.dplyr, partCNV, RcisTarget, receptLoss, RgnTX, RiboProfiling, rols, RTCGA, scviR, SGCP, SigFuge, soGGi, spatzie, SQLDataFrame, supersigs, tLOH

**EXPERIMENT DATA:** 

Four experimental data packages were removed from this release (after being
deprecated in BioC 3.22):

- AneuFinderData, chromstaRData, curatedCRCData, RGMQLlib, rRDPData

Previously announced deprecated package was fixed and restored to Bioconductor 3.23: 

- curatedCRCData


Two experimental data packages are deprecated in this release and will be
removed in Bioc 3.24:

- curatedBreastData, Fletcher2013b

**ANNOTATION DATA:** 

No annotation annotation packages were removed from this release (after being deprecated
in Bioc 3.22).

Four annotation packages are deprecated in this release and will be removed in
Bioc 3.24:

- hpAnnot, humanCHRLOC, mouseCHRLOC, ratCHRLOC 

We will also be removing seven out-dated TxDb packages from Bioc 3.24:

- TxDb.Mmusculus.UCSC.mm10.ensGene, TxDb.Athaliana.BioMart.plantsmart22, TxDb.Athaliana.BioMart.plantsmart25, TxDb.Cfamiliaris.UCSC.canFam3.refGene, TxDb.Rnorvegicus.UCSC.rn4.ensGene, TxDb.Hsapiens.BioMart.igis, TxDb.Rnorvegicus.BioMart.igis


**WORKFLOWS:** 

One workflow package was removed from this release (after being deprecated in
Bioc 3.22).

- spicyWorkflow

No workflow packages have been deprecated in this release:

**BOOKS:**

No books were removed from this release.

No books were deprecated in this release.
