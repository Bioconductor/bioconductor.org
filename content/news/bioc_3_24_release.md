November 25, 2026

**Bioconductor:**

We are pleased to announce Bioconductor 3.24, consisting of
2521 software packages, 442 experiment data packages, 929 annotation
packages, 29 workflows and 8 books.

There are 149 new software packages, 8 new data experiment packages,
6 new annotation packages, 1 new workflow, no new books, and many updates and
improvements to existing packages.

Bioconductor 3.24 is compatible with R 4.6, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher, macOS arm64 and Linux
arm64. This release will also include updated Bioconductor [Docker containers][2].

Note: Bioconductor will be testing the propagation of source and binary
generation through the [r-universe system](https://bioc-release.r-universe.dev/builds). 

Thank you to everyone for your contribution to Bioconductor.

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.24](#getting-started-with-bioconductor-324)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.24
======================================

To update to or install Bioconductor 3.24:

1. Install R 4.6. Bioconductor 3.24 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 149 new software packages in this release of Bioconductor.

- [alabaster.treese](/packages/alabaster.treese) Save
  TreeSummarizedExperiment into file artifacts, and load them back
  into memory. This is a more portable alternative to serialization
  of such objects into RDS files. Each artifact is associated with
  metadata for further interpretation; downstream applications can
  enrich this metadata with context-specific properties.

- [AnnotationGx](/packages/AnnotationGx) A package for building,
  updating and querying an annotation database for pharmaco-genomic
  data. Provides API interfaces to various bioinformatic and
  chemoinformatic databases, as well as a set of functions for
  curating your data using the annotations retrieved from those APIs.

- [AnVILVRS](/packages/AnVILVRS) Process Variant Call Format (VCF)
  files and perform lookup operations on Genomic Variation
  Representation Service (GA4GH VRS) identifiers. The GA4GH VRS
  identifiers provide a standardized way to represent genomic
  variations, making it easier to exchange and share genomic
  information.

- [apbCatalogR](/packages/apbCatalogR) Looks up the fields of AnnData
  Proteomics Bridge (APB) results by meaning, such as a PEP or an
  unnormalized intensity, without naming vendor columns. Reimplements
  the apb-catalog API on the bundle apb-catalog releases: reviewed
  sets, apb2's rule levels and expected answers. Reads APB results,
  as h5ad, Parquet or DuckDB, into anndataR objects.

- [augere.core](/packages/augere.core) Provides utility functions for
  automatic generation of reports in the augere framework. This
  mostly involves parsing and dynamically editing Rmarkdown files,
  controlling inputs into the augere pipelines. Each pipeline
  generates and executes a self-contained Rmarkdown file for a common
  bioinformatics analysis; see downstream packages like augere.de and
  augere.screen for examples of some analysis pipelines.

- [augere.de](/packages/augere.de) Implements pipelines for
  generating differential expression analysis reports in the augere
  framework. This includes analyses with edgeR or voom-limma, with a
  variety of options for contrasts, blocking and covariates. Each
  pipeline function generates a self-contained Rmarkdown report with
  all of the steps required to reproduce the DE analysis.

- [augere.gsea](/packages/augere.gsea) Implements pipelines for
  generating gene set enrichment analysis reports in the augere
  framework. This includes various competitive and self-contained
  gene sets from a variety of Bioconductor packages. Each pipeline
  function generates a self-contained Rmarkdown report with all of
  the steps required to reproduce the gene set enrichment analysis.

- [augere.screen](/packages/augere.screen) Implements pipelines for
  generating functional screen analysis reports in the augere
  framework. This uses voom to test for differential abundance of
  barcodes with consolidation into gene-level statistics. Each
  pipeline function generates a self-contained Rmarkdown report with
  all of the steps required to reproduce the analysis.

- [augere.solo](/packages/augere.solo) Implements pipelines for
  generating single-cell analysis reports in the augere framework.
  This uses scrapper to execute routine steps such as quality
  control, normalization, feature selection, clustering and marker
  detection. We also implement a pipeline for automatic cell type
  annotation against a labelled reference with SingleR. Each pipeline
  function generates a self-contained Rmarkdown report with all of
  the steps required to reproduce its analysis.

- [BamScale](/packages/BamScale) Multithreaded sequential BAM
  processing built on top of the ompBAM C++ engine. BamScale provides
  user-friendly BAM read and scan interfaces designed for
  compatibility with existing Bioconductor workflows.

- [barmixR](/packages/barmixR) Implements the Bayesian modeling
  framework underlying the barmixR (BARcode MIXture analysis)
  platform for high-throughput quantitative analysis of
  genotype-specific treatment responses in pooled cancer cell
  populations. The package integrates barcode sequencing count data
  with volumetric measurements such as tumor volume (in vivo) or
  cellular confluency (in vitro) using hierarchical probabilistic
  models. Barcode counts are modeled with a Dirichlet–multinomial
  distribution to account for compositional sequencing data, while
  volumetric measurements are modeled using log-normal (tumor volume)
  or beta (confluency) likelihoods. Posterior inference is performed
  using Hamiltonian Monte Carlo through 'rstan'. The resulting
  posterior distributions enable estimation of clone-specific
  quantitative treatment resistance (QTR) together with uncertainty
  propagation from both sequencing and volumetric data. Additional
  functions provide posterior predictive checks, estimation of
  resistance ratios, treatment ranking, and visualization of
  resistance landscapes using violin plots and bubble heatmaps. The
  methods are designed for multiplexed lineage-tracing experiments in
  cancer research and were developed to analyze treatment resistance
  in gastrointestinal stromal tumors (GIST), but are broadly
  applicable to barcoding-based studies of treatment response and
  clonal dynamics across diverse cancer types.

- [BatchVaria](/packages/BatchVaria) Evaluates batch correction in
  omics data by quantifying how a correction redistributes variance,
  not only whether it removes it. Records correction and
  variance-decomposition provenance in the metadata of an ordinary
  SummarizedExperiment, and harmonises both the batch correction
  methods and the variance-profiling engines behind documented
  contracts, so that either set can be extended without modifying the
  package. Variance fractions are compositional, so removing variance
  associated with one term inflates the fraction attributed to every
  other. Each decomposition is therefore recorded alongside the total
  variance it is relative to and the share an unbalanced design
  cannot attribute to any single term, which is what distinguishes a
  correction that removed unwanted variation from one that merely
  removed variation.

- [BenchHub](/packages/BenchHub) The trio is the combination of a
  data set, a metric and supporting evidence which provides some best
  case scenario, if not the ground truth itself. BenchHub has data
  downloaders for FigShare, G.E.O., and ExperimentHub. Caching is
  used to avoid lengthy downloads after the first time a data set is
  accessed. The user may also specify their own data set and
  supporting evidence. The Benchmark Insights module provides
  functionality for comparing and contrasting the performance of
  alternative algorithms.

- [bHIVE](/packages/bHIVE) bHIVE implements an artificial immune
  network (AI-Net) for clustering and classification of
  high-dimensional data such as mass cytometry and single-cell
  profiles. A population of prototypes, called antibodies, evolves by
  clonal selection, somatic hypermutation and network suppression,
  and the learned prototypes summarize cell populations in the
  original marker units. The same fitted model clusters a dataset and
  labels new cells. bHIVE accepts SummarizedExperiment and
  SingleCellExperiment objects and integrates with bluster through
  bHIVEParam. A C++ backend computes affinities, mutations and
  suppression, and optional R6 modules add mechanisms such as
  germinal center selection, microenvironment-dependent mutation and
  isotype switching. The package also provides multilayer clustering,
  grid search and caret model wrappers.

- [BiocDuckDB](/packages/BiocDuckDB) Integration package providing
  high-level Parquet I/O functions and optimized methods for
  single-cell analysis workflows using DuckDB-backed data structures.
  Includes readParquet and writeParquet functions for seamless
  serialization of SummarizedExperiment, SingleCellExperiment,
  MultiAssayExperiment, and MultiAssaySpatialExperiment objects, plus
  SQL-optimized implementations of scran and scuttle methods for
  variance modeling, marker detection, QC metrics, and normalization.
  This package brings together DuckDBDataFrame, DuckDBArray,
  DuckDBGRanges, and DuckDBSpatial (optional) for complete
  Bioconductor integration.

- [BreedingExperiment](/packages/BreedingExperiment) Provides the
  BreedingExperiment class, a container that keeps marker genotypes,
  phenotypes, pedigree and marker coordinates together in one object,
  extending RangedSummarizedExperiment so that markers are
  addressable as genomic ranges. Relationship matrices central to
  quantitative genetics are computed as methods on the container: the
  pedigree numerator relationship matrix of Henderson (1976)
  <doi:10.2307/2529339>, the genomic relationship matrix of VanRaden
  (2008) <doi:10.3168/jds.2007-0980>, the dominance relationship
  matrix of Vitezica, Varona and Legarra (2013)
  <doi:10.1534/genetics.113.155176>, and the single-step matrix
  combining pedigree and genomic information of Legarra, Aguilar and
  Misztal (2009) <doi:10.3168/jds.2009-2061> and Christensen and Lund
  (2010) <doi:10.1186/1297-9686-42-2>, which allows genotyped and
  ungenotyped individuals to be analysed together as in Aguilar and
  others (2010) <doi:10.3168/jds.2009-2730>. Marker quality control,
  imputation, inbreeding coefficients and allele frequency summaries
  are included.

- [CanvasXpressBio](/packages/CanvasXpressBio) Bridges the
  CanvasXpress interactive visualization library to core Bioconductor
  data structures. Provides a generic to render a
  SummarizedExperiment, and any class that extends it such as
  SingleCellExperiment or DESeqDataSet, as an interactive
  CanvasXpress heatmap that carries its row and column annotations,
  together with a volcano helper for differential-expression result
  tables. Charts are reproducible CanvasXpress widgets usable in R, R
  Markdown and Shiny.

- [causalASCEND](/packages/causalASCEND) Causal ASCEND (Ancestral
  Scalable Causal discovEry via iNherited Descent) infers ancestral
  (cause-before-effect) relations among downstream molecular
  variables, such as gene expression, when a layer of upstream
  variables, such as genotypes or methylation sites, is known to
  precede them. Each conditional independence test conditions only on
  a small, iteratively refined set of nearest ancestors instead of
  the whole upstream layer, so the number of tests grows polynomially
  with the number of variables. The package accepts matrices,
  SummarizedExperiment objects or a MultiAssayExperiment, and returns
  the ancestral relations with the rule and evidence behind each one.

- [cellchatr](/packages/cellchatr) Fast identification of
  over-expressed genes and ligand-receptor interactions in
  single-cell RNA-seq data, the pre-processing steps of cell-cell
  communication analysis tools such as CellChat. Works directly on
  SummarizedExperiment and SingleCellExperiment objects and accepts
  any ligand-receptor table (for example the CellChatDB interaction
  table). The per-gene, per-cell-group Wilcoxon rank-sum tests are
  implemented in Rust via extendr and run in parallel with Rayon.

- [cellNexus](/packages/cellNexus) Provides access to a copy of the
  Human Cell Atlas, but with harmonised metadata. This allows for
  uniform querying across numerous datasets within the Atlas using
  common fields such as cell type, tissue type, and patient
  ethnicity. Usage involves first querying the metadata table for
  cells of interest, and then downloading the corresponding cells
  into one of several supported formats, including
  SingleCellExperiment, Seurat, or pseudobulk SummarizedExperiment
  objects. Note: The Human Cell Atlas data accessed through this
  package is subject to its own licensing terms (typically Creative
  Commons Attribution or similar open data licenses as specified by
  the Human Cell Atlas Data Use Agreement), which differ from the
  package license (GPL-2). See the package documentation for details.

- [CGRphylo2](/packages/CGRphylo2) An alignment-free phylogenetic
  analysis method for viral genomes using Chaos Game Representation
  (CGR), a technique based on statistical physics concepts. Viruses
  exhibit high mutation rates, facilitating rapid evolution and
  emergence of new species, subspecies, strains, and recombinant
  forms. Accurate classification is crucial for understanding viral
  evolution and therapeutic development. Traditional phylogenetic
  methods require sequence alignment, which is computationally
  intensive. CGRphylo2 addresses this by implementing CGR-based
  whole-genome comparison that is fast, accurate, and computationally
  efficient. The package successfully classifies closely related
  viral lineages (demonstrated on SARS-CoV-2 lineages A and B),
  identifies recombinants (such as the XBB variant), and
  distinguishes multiple strains simultaneously. It processes
  sequences 5-13.7x faster than alignment-based methods
  (Clustal-Omega) with linear computational scaling. As a k-mer based
  approach, it enables simultaneous comparison of numerous
  closely-related sequences of different lengths. The package creates
  frequency matrices for distance calculations and phylogenetic tree
  construction, with outputs compatible with standard formats (MEGA,
  PHYLIP, Newick). Methods are based on Thind and Sinha (2023)
  <doi:10.2174/0113892029264990231013112156>.

- [circaTime](/packages/circaTime) Provides a single
  SummarizedExperiment-based interface to published circadian phase
  (time-of-sample-collection) inference methods for bulk and
  donor-level pseudobulk transcriptomes, along with a shared
  cross-validation benchmarking harness, cross-dataset transfer
  testing, and a calibrated uncertainty layer: bootstrap and Bayesian
  confidence intervals validated for real coverage, a fixed or
  dynamically-selected clock-gene quality flag, and consensus phase
  estimation across every bundled engine. Does not implement a new
  phase-inference algorithm; unifies existing published methods that
  are otherwise scattered across incompatible standalone
  implementations.

- [CLAMP](/packages/CLAMP) CLAMP performs prior-informed latent
  variable decomposition of high-dimensional transcriptomic data. It
  integrates curated gene sets to learn biologically interpretable
  latent variables, supports file-backed matrices for large datasets,
  and provides tools for preprocessing, normalization, projection,
  and evaluation of latent structures. CLAMP is designed to scale to
  tens of thousands of samples, making it suitable for large public
  resources such as recount3 and ARCHS4. It enables researchers to
  uncover biologically meaningful patterns that connect genes,
  pathways, and complex traits in transcriptomics studies.

- [CMEnt](/packages/CMEnt) CMEnt implements a correlation-based
  method for identifying Differentially Methylated Regions (DMRs)
  from genomic seeds, commonly Differentially Methylated Positions
  (DMPs). The package expands regions around significant seeds
  considering both statistical significance and biological relevance
  of methylation changes. It supports array-based (450K, EPIC,
  EPICv2) and NGS methylation data through tabix-indexed files,
  provides DMR interaction analysis via motif similarity, and
  includes comprehensive visualization tools including circos plots
  and beta value heatmaps.

- [CONCERTDR](/packages/CONCERTDR) Provides a complete pipeline for
  drug repurposing analysis using the Connectivity Map (CMap) L1000
  database. Functionality includes database subsetting, reference
  profile extraction from GCTX files, signature matching via
  Kolmogorov-Smirnov, GSEA-based, cosine similarity, and other
  methods, result annotation against CMap compound metadata, z-score
  matrix extraction, and publication-quality clustermap visualization
  via ComplexHeatmap. Connectivity scoring methods (KS, GSEA, XCos,
  XSum, Zhang) are implemented internally based on the algorithms
  described in Lamb et al. (2006), Subramanian et al. (2005), Cheng
  et al. (2014), and Zhang & Gant (2008), removing the need for
  external GitHub-only dependencies. CamSum, a correlation-adjusted
  sum score with an analytic p-value, is also provided.

- [consensusRegions](/packages/consensusRegions) Merges the evidence
  carried by replicated ChIP-seq, ATAC-seq or CUT&RUN experiments
  into a single set of consensus regions. When several replicates
  agree on a position their p-values are combined, so that a peak too
  weak to survive on its own can still be retained. Fisher, Stouffer,
  Lancaster and rank-product combinations are available; the weighted
  ones accept per-replicate weights taken from FRiP scores, library
  sizes, or the peak sets themselves when no alignment file is at
  hand. Input can be anything from a full narrowPeak to a bare BED3:
  scores are rank-transformed when p-values are missing, and a
  weighted k-of-n presence rule takes over when even the score is
  absent. The threshold on the combined significance can be
  calibrated against a permutation null rather than fixed a priori.
  Quality control plots and export to SummarizedExperiment are
  included. The approach follows MSPC (Jalili et al., 2015). It
  differs from consensus sets built on occupancy, such as the
  minOverlap rule of DiffBind, in that a region is kept on the
  combined significance of the peaks that overlap it rather than on
  how many replicates called it, and from IDR, which ranks a pair of
  replicates by reproducibility, in that any number of replicates is
  reduced to one set of regions.

- [CorNetto](/packages/CorNetto) Builds knowledge-guided multi-omic
  correlation networks from normalized transcriptomic, proteomic, and
  metabolomic abundance data. Group-specific correlation networks are
  compared with the Fisher z-difference test, either across all
  within-omic pairs in one or every assay without generating
  cross-omic pairs, or across candidate edges supplied by a
  prior-knowledge network, and the resulting differential network is
  summarized into node-level rewiring scores that identify features
  whose interaction patterns change between groups. Rewiring scores
  can be compared with group-label permutation reference
  distributions, and networks can be restricted to pathway-focused
  neighbourhoods, converted to 'igraph' objects, or exported as
  Cytoscape-ready node and edge tables. The package is built on
  Bioconductor containers so that multi-omic assays and sample
  metadata are managed consistently throughout the workflow.

- [crossarrayEWAS](/packages/crossarrayEWAS) Estimation and
  region-calling routines for epigenome-wide association studies that
  combine Illumina HumanMethylation450 and MethylationEPIC samples
  and follow subjects over more than one visit. Probe-level effects
  are fitted within subject, so subject identity, array type and chip
  are absorbed rather than modelled; adjacent probes are grouped by
  within-subject co-methylation; regions are obtained by
  total-variation denoising with a distance-decaying penalty, which
  places boundaries where the effect track jumps instead of
  thresholding a smooth curve; and long-range blocks are called from
  a three-state hidden semi-Markov model whose transition
  probabilities decay with genomic distance. Selection stability and
  permutation nulls resample subjects, not samples, so a subject's
  visit series is never split.

- [ctdR](/packages/ctdR) Identifies chemicals significantly
  associated with a set of genes using data from the Comparative
  Toxicogenomics Database (CTD, <https://ctdbase.org>). Supports four
  enrichment methods through a unified interface: Over-Representation
  Analysis (ORA) via a hypergeometric test on stats::phyper(), Gene
  Set Enrichment Analysis (GSEA) via 'fgsea', the competitive
  gene-set test CAMERA via 'limma', and per-sample Gene Set Variation
  Analysis (GSVA) via 'GSVA'. Users download the CTD chemical-gene
  interactions file once, import it with import_CTD(), and then run
  enrichment_CTD() on either a gene list (ORA, GSEA) or an expression
  matrix (CAMERA, GSVA). This package does not bundle or redistribute
  CTD data: its example interaction file is synthetic, and the tables
  of the replication vignette hold counts and p-values only. Users
  must download the required data files directly from
  <https://ctdbase.org> and comply with the CTD data licensing terms
  (see <https://ctdbase.org/about/legal.jsp>).

- [CyFj11](/packages/CyFj11) Provides functions to parse FlowJo v11
  workspace files (.flowjo) and extract gating information,
  populations, and associated FCS files. Includes utilities for
  exporting data to FlowJo v10 workspace format. The package enables
  bidirectional data flow between FlowJo and R/Bioconductor,
  facilitating reproducible analysis workflows.

- [CySA](/packages/CySA) CySA provides an interactive Shiny-based
  tool for analyzing flow and mass cytometry data. It supports
  self-organizing map (SOM) based clustering, dimension reduction
  (t-SNE, UMAP, PCA), and statistical comparison of cell populations
  across sample groups. Users can select cell clusters interactively
  via 2D scatter plots, dendrograms, or SOM node grids, and explore
  marker expression patterns across the selected populations.

- [DaparToolshed](/packages/DaparToolshed) The package DaparToolshed
  is a Bioconductor distributed R package which provides all the
  necessary functions to analyze quantitative data from label-free
  proteomics experiments. It is an update of our previous package
  DAPAR and contains more functions to analyze the data and uses
  MultAssayExperiment and SummarizedExperiment data structures.
  Contrarily to most other similar R packages, it is endowed with
  rich and user-friendly graphical interfaces, so that no programming
  skill is required (see `Prostar` package).

- [DAssemble](/packages/DAssemble) The DAssemble package implements
  an ensemble framework for differential-abundance and
  differential-expression analysis across bulk RNA-seq, single-cell
  RNA-seq, and microbiome studies. It wraps a collection of popular
  DA/DE methods as core methods and combines them via Cauchy
  Combination Tests (CCT), optionally augmented by simple enhancers
  such as Wilcoxon, Kolmogorov–Smirnov, and presence–absence logistic
  regression.

- [daVis](/packages/daVis) This package contains utility functions to
  visualize the output from differential expression analysis. The
  input can be either a model or a list of top tables or a
  combination of the two. The model can be output from limma, edgeR
  or DESeq2.

- [dbSequence](/packages/dbSequence) Provides a DuckDB-backed
  infrastructure for working with large-scale genomic data that
  exceeds available memory. Supports lazy evaluation of genomic
  ranges and efficient overlap queries. Integrates with Bioconductor
  classes (GRanges, BiocIO) and provides plyranges-compatible API.

- [DEP](/packages/DEP) This package provides an integrated analysis
  workflow for robust and reproducible analysis of mass spectrometry
  proteomics data for differential protein expression or differential
  enrichment. It requires tabular input (e.g. txt files) as generated
  by quantitative analysis softwares of raw mass spectrometry data,
  such as MaxQuant or IsobarQuant. Functions are provided for data
  preparation, filtering, variance normalization and imputation of
  missing values, as well as statistical testing of differentially
  enriched / expressed proteins. It also includes tools to check
  intermediate steps in the workflow, such as normalization and
  missing values imputation. Finally, visualization tools are
  provided to explore the results, including heatmap, volcano plot
  and barplot representations. For scientists with limited experience
  in R, the package also contains wrapper functions that entail the
  complete analysis workflow and generate a report. Even easier to
  use are the interactive Shiny apps that are provided by the
  package.

- [DiPPER](/packages/DiPPER) Implements a Bayesian hierarchical model
  for differential prevalence analysis, tailored for microbiome
  studies. The method accommodates both cross-sectional and
  longitudinal designs, providing finite effect estimates and
  uncertainty intervals that are naturally adjusted for multiple
  testing via hierarchical shrinkage. Posterior inference is
  performed via Markov chain Monte Carlo (MCMC) using 'CmdStan' as
  the backend.

- [dnaEPICO](/packages/dnaEPICO) A modular and reproducible workflow
  for preprocessing and analysing Illumina DNA methylation array data
  from the EPICv2, EPIC, and 450K platforms. It supports quality
  control, probe filtering, cell-type deconvolution, phenotype
  preparation, CpG-wise generalised linear models with optional
  term-level omnibus F tests, longitudinal mixed-effects models using
  'lmerTest'/'lme4' or 'nlme', optional 'nlme' residual correlation
  structures, and automated reports. The workflow uses established
  Bioconductor tools, including 'minfi', 'ENmix', and 'wateRmelon',
  and runs locally or in high-performance computing environments.

- [DuckDBArray](/packages/DuckDBArray) Provides DuckDB-backed
  implementations of DelayedArray and DelayedMatrix classes for
  efficient, out-of-memory operations on large array and matrix
  datasets. The package includes DuckDBArray, DuckDBArraySeed, and
  DuckDBMatrix classes with optimized matrixStats methods. Supports
  sparse arrays, grid partitioning, and parallel processing for
  scalable genomic data analysis.

- [DuckDBDataFrame](/packages/DuckDBDataFrame) Provides DuckDB-backed
  implementations of core tabular data structures including
  DuckDBTable, DuckDBDataFrame, and DuckDBColumn. These classes
  extend Bioconductor's S4Vectors framework to enable efficient,
  out-of-memory operations on large tabular datasets using DuckDB's
  analytical database engine. The package includes SQL function
  discovery, connection management, and support for list columns and
  embeddings stored in ARRAY[] columns.

- [DuckDBGRanges](/packages/DuckDBGRanges) Provides DuckDB-backed
  implementations of GenomicRanges data structures for efficient,
  out-of-memory operations on large genomic interval datasets. The
  package includes DuckDBGRanges and DuckDBGRangesList classes with
  optimized range operations including restrict, flank, shift,
  findOverlaps, coverage, and distance calculations. Ideal for
  variant analysis, scATAC-seq peak analysis, and other genomic
  workflows requiring memory-efficient range operations on large
  datasets.

- [DuckDBSpatial](/packages/DuckDBSpatial) Provides DuckDB-backed
  implementations of spatial data structures for efficient
  out-of-memory operations on large spatial datasets. The package
  includes 68 sf-compatible spatial methods (st_intersects,
  st_contains, st_within, etc.) and GeoParquet 1.0 I/O support via
  writeGeoParquet(). Ideal for spatial analysis, genomic interval
  analysis, and integration with BiocDuckDB for
  MultiAssaySpatialExperiment workflows requiring memory-efficient
  operations on large datasets.

- [ECLIPSE](/packages/ECLIPSE) This package is used for the analysis
  of single-cell T cell receptor sequencing data (scTCR-seq). ECLIPSE
  maps chain pairings across a sample to predict the true clonotype
  of each cell, unlike other methods which only consider the chains
  that are sequenced in each cell. Seurat/SingleCellExperiment
  objects and Cell Ranger vdj/multi contig files are required as
  inputs.

- [EMMA](/packages/EMMA) EMMA is a package that provides a
  provenance-aware execution framework for functional enrichment
  analyses. It allows users to run enrichment analyses using existing
  tools and workflows while explicitly capturing analytical
  parameters during runtime, and returns the native enrichment
  results output together with structured metadata.

- [EMTscore](/packages/EMTscore) Epithelial-mesenchymal transition
  (EMT) is an important form of cellular plasticity that is fully or
  partially activated in several biological scenarios including
  development and disease progression. EMT involves altered
  expression of hundreds of protein-coding and non-protein-coding
  genes. Recent studies showed the prevalence of partial EMT in
  multiple processes such as various cancers and organ fibrosis,
  which necessitates rigorous quantification of the degree of EMT.
  While traditional gene set scoring methods such as gene set
  variation analysis have been used to generate EMT scores from omics
  data, multiple EMT scoring algorithms and EMT gene sets have been
  used by different groups without standardization. Furthermore,
  comparisons of EMT scores computed from different methods and/or
  different EMT gene sets are generally difficult due to both the
  context dependent nature of EMT and the lack of tools that
  comprehensively integrate varying components for EMT scoring. To
  address this problem, EMTscore enables users to select scoring
  methods from a list of previously used algorithms and EMT gene sets
  from a list of gene sets produced from different experiments.
  Several visualization methods are provided for making
  publication-quality plots of EMT scores from omics data. The
  package also implements a principal-component-analysis based method
  for scoring divergent EMT processes from a single dataset. Overall,
  EMTscore provides an integrated solution for assessing the degree
  and complexity of EMT from omics data, and paves the way for
  standardizing the comparison of EMT programs across multiple
  contexts.

- [enrichmet](/packages/enrichmet) Provides functions to perform
  metabolite-based pathway enrichment analysis using Fisher's exact
  test and MetSEA. Computes relative betweenness centrality (RBC) of
  input metabolites in a pathway-metabolite interaction network. The
  package generates visualizations including pathway dot plots,
  MetSEA plots, metabolite centrality plots, and network plots. It
  supports KEGG-based pathway inputs and generates exportable
  results.

- [exploreSE](/packages/exploreSE) Interactive exploration of a
  Summarized experiment and associated DE analyses. A shiny app to
  inspect and interact with pre-calculated differential expression
  and functional enrichment results, to facilitate comparisons
  between multiple models, and to include known associations, such as
  genetic risk factors. This is meant to be a tool for exploration
  and discussion, rather than analysis.

- [Fancy](/packages/Fancy) Infers microbial co-abundance networks by
  combining mutual information (via MRNET) and distance correlation
  into a single hybrid edge score. Bootstrap resampling provides
  frequency-based confidence, and a weighted scoring function ranks
  edges for downstream module analysis. Designed for
  metagenome-assembled genome (MAG) count tables.

- [fishash](/packages/fishash) Assigns guide RNAs or other genetic
  perturbations to cells in single-cell sequencing experiments using
  a one-sided Fisher's exact test. Implements an iterative refitting
  procedure to mitigate Simpson's paradox, supports multiple false
  discovery rate correction methods (Benjamini-Hochberg,
  Benjamini-Yekutueli, and Guo & Sarkar 2020), and provides
  simulation utilities for benchmarking demultiplexing methods.
  Results are returned as SummarizedExperiment objects.

- [gDRplots](/packages/gDRplots) gDRplots provides static
  visualization functions and dedicated helpers for the gDR drug
  response platform. Designed for high-throughput screening (HTS),
  pharmacogenomics, and chemical genomics research, it enables
  comprehensive analysis of single-agent, co-dilution, and drug
  combination cell viability assays. The package supports quality
  control (plate layouts, control mapping), data normalization
  inspection, dose-response curves, distribution boxplots, drug
  synergy heatmaps with overlaid isobolograms, and combination index
  plots. Additionally, it integrates with PRISM and DepMap datasets
  to generate volcano plots and correlation scatter plots for
  molecular feature associations.

- [GDSAnnotator](/packages/GDSAnnotator) Fast and memory-efficient
  manipulation of whole-genome variant annotation data stored in
  Genomic Data Structure (GDS) files. The annotated GDS (aGDS) format
  introduced by FAVOR stores the FAVOR annotations of the variants
  together with the genotypes in one SeqArray GDS file. This package
  generalizes the format to annotation GDS files of any source, with
  or without genotypes. Annotation resources and annotator outputs,
  including FAVOR, gnomAD, Ensembl VEP and SnpEff, are converted into
  compact SeqArray-compatible GDS files with random access. Variants
  can be given as chr-pos-ref-alt strings, data frames, GRanges
  objects, VCF files or objects, or SeqArray GDS files. They are
  annotated by querying the annotation files, and the annotations can
  be written back into a genotype GDS file. Summary statistics of the
  variants and their annotations are computed in a single pass over a
  file and rendered as a VEP/SnpEff-like report. The package is
  intended for large-scale whole-genome sequencing studies.

- [gdscloud](/packages/gdscloud) Provides read-only access to GDS
  (Genomic Data Structure) files stored on cloud storage services
  including Amazon S3, Google Cloud Storage (GCS), and Azure Blob
  Storage, as well as any HTTP/HTTPS URL. It extends the 'gdsfmt'
  package so that cloud URLs (http://, https://, s3://, gs://, az://)
  can be opened transparently, without downloading the whole file
  first. Only the blocks that are actually read are fetched, using
  HTTP Range requests via libcurl together with an in-memory
  least-recently-used block cache, so that random access to a remote
  GDS file behaves like access to a local one. Credentials are
  resolved from the usual environment variables of each service, or
  set per session and per URL prefix, and they can be exported to the
  workers of a parallel cluster.

- [geneClusterPattern](/packages/geneClusterPattern) An R package for
  constructing and comparing ortholog gene patterns around a focal
  gene across species. It enables the exploration of conserved gene
  neighborhoods (synteny) by organizing orthologous genes within a
  defined genomic window and aligning their relative positions across
  multiple genomes. It can retrieve homologous gene directly from
  Ensembl or it can incorporate user-provided ortholog group
  assignments.

- [gesel](/packages/gesel) Query the gesel indices to identify
  interesting gene sets. Users can test for enrichment in their own
  list of genes and/or search by text in the set names or
  descriptions. The indices can be hosted remotely or on a local
  filesystem.

- [GExPipe](/packages/GExPipe) Shiny application (GExPipe) for
  high-throughput genomic analysis of bulk RNA-seq and microarray
  data (e.g. from GEO). Four analysis types: RNA-seq only, microarray
  only, Merged (Both) (per-dataset normalize then one joint batch and
  limma DE), and Parallel DE then merge (separate pipelines through
  DE, then RNA-seq intersect microarray). Each RNA-seq and microarray
  box accepts one or more GSE IDs in the same run. Integrates with
  Bioconductor (GEOquery, Biobase, limma, DESeq2, edgeR,
  clusterProfiler) for download, QC, normalization, batch correction,
  differential expression, WGCNA, pathway enrichment, PPI, and
  machine learning. Uses common data structures (ExpressionSet,
  DGEList) for interoperability. For a full dependency tree
  (including STRINGdb + PPI helpers), use
  BiocManager::install("GExPipe", dependencies = TRUE). STRING data
  are downloaded on first PPI use (internet required); they cannot be
  bundled in the package. Microarray CEL normalization uses affy
  and/or oligo when supplementary CEL files are available.

- [GPlinksR](/packages/GPlinksR) GPlinksR constructs gene-peak
  regulatory networks for ATAC-RNA integration by combining
  enhancer-based, promoter-based, and proximity (closest-gene)
  mappings. The package accepts direct peak and gene vectors as well
  as container-based inputs through a wrapper for common Bioconductor
  object classes. Enhancer-gene links are obtained from the PEREGRINE
  enhancer-gene datasets provided by AnnoQ, while promoter and gene
  coordinates are retrieved from EnsDb.Hsapiens.v86.

- [grayleafspotr](/packages/grayleafspotr) Quantitative phenotyping
  of gray leaf spot (Magnaporthe oryzae) fungal colonies grown on
  petri dishes. The package segments colonies from time-lapse plate
  photographs using a bundled SmallUNet deep-learning model, extracts
  morphometric and texture features (area, eccentricity, crack
  coverage, radial profile), and provides tidy result objects
  together with template ggplot2 visualisations. Python dependencies
  are managed automatically through basilisk; no manual environment
  setup is required.

- [growkar](/packages/growkar) growkar provides Bioconductor-native
  infrastructure for high-throughput microbial growth phenotyping
  from plate-based assays such as optical density (OD) measurements.
  The package is built around SummarizedExperiment as its canonical
  data model, with OD measurements in assays, timepoint annotations
  in rowData, sample annotations in colData, and derived phenotypes
  in metadata. growkar supports data import, quality control,
  normalization, empirical growth-rate estimation, exponential-phase
  detection, and parametric modeling of growth dynamics. These
  workflows enable extraction of biologically meaningful phenotypes
  including lag time, growth rate, doubling time, and carrying
  capacity, and allow growth-derived phenotypes to be integrated with
  other Bioconductor analyses in microbial systems biology and
  multi-omics studies. The data-structure and analysis layers depend
  only on core Bioconductor infrastructure; graphing is optional and
  provided by 'plot_*' functions that are enabled when the suggested
  graphics package is installed.

- [GSEAlens](/packages/GSEAlens) GSEAlens provides an interactive
  exploration layer on top of standard Bioconductor RNA-seq
  workflows. It consumes fitted model objects from limma (MArrayLM)
  or DESeq2 (DESeqDataSet) as input and accepts expression matrices
  and sample metadata as SummarizedExperiment objects, ensuring
  interoperability with core Bioconductor data containers. For core
  computation, GSEAlens wraps clusterProfiler::GSEA() as its
  statistical engine (thereby inheriting the fgsea fast GSEA
  methodology) and draws on MSigDB gene set collections via the
  msigdbr package from CRAN; multi-contrast parallel computation is
  handled by future (future::multisession). Visualization output
  relies on Bioconductor graphics packages including enrichplot,
  ComplexHeatmap, and circlize, producing figures suitable for
  publication pipelines. The package also includes a built-in Shiny
  application for interactive exploration of enrichment results after
  DEG analysis, with the ability to export self-contained
  reproducible R scripts.

- [HTAN](/packages/HTAN) A client for the open-access data of the
  Human Tumor Atlas Network (HTAN), Phase 1. File, biospecimen,
  clinical and provenance metadata are queried locally from a
  versioned parquet snapshot of the ISB-CGC HTAN BigQuery release
  tables, through DuckDB and dbplyr, with no Google account.
  Open-access files are downloaded from the CRDC General Commons
  (anonymous GA4GH DRS) or Synapse (personal access token), and
  common HTAN file formats are read into Bioconductor classes such as
  SingleCellExperiment, SpatialExperiment and MultiAssayExperiment.

- [iconic](/packages/iconic) Provides a model selection workflow for
  causal inference with genetic instruments and negative controls in
  observational omics data. The package fits eight estimators of the
  natural direct and indirect effects (NDE/NIE), diagnoses which are
  valid for the user's data via instrument-strength checks and
  negative-control validity screens (A1, A2, A2'), stress-tests them
  against confounding and pleiotropy violations, and recommends the
  estimator most likely to be unbiased. A hybrid generative texture
  model (a torch GAN for sample-level structure and a Gaussian copula
  for the mediator panel) lets the sensitivity analysis mirror the
  marginal and joint structure of the user's cohort. The modular API
  comprises iconic_data() to standardize data, iconic_diagnose() to
  assess estimator eligibility, iconic_estimate() to fit all eligible
  estimators, iconic_sensitivity() to sweep instrument exogeneity,
  iconic_recommend() to rank estimators by identification strength
  and robustness, and iconic_prospect() to plan future studies.
  infer_confounding() estimates held-fixed confounding parameters
  from the data. Multi-omic cohorts are supplied as
  MultiAssayExperiment objects (omics panels as experiments,
  sample-level exposure, instruments and covariates as colData), and
  analysis data can be exported back to MultiAssayExperiment. Helper
  functions construct exposure instruments (polygenic scores from
  GWAS summary statistics or PGS Catalog scoring files), mediator
  instruments (GenomicRanges-based cis-eQTL scans and elastic-net
  predicted-expression composites from RangedSummarizedExperiment or
  matrix inputs), and negative-control panels (principal components
  of residualized omics matrices). Time-to-event outcomes are
  supported on the Cox log-hazard-ratio and
  restricted-mean-survival-time scales. Binary outcomes are supported
  on the log-odds-ratio (logistic two-stage predictor-substitution
  outcome stage) and risk-difference (linear probability model)
  scales. The torch package is required for the generative texture
  model (sensitivity and prospective analysis); all other
  functionality works without it. The methods are described in
  Bresnahan et al. (2026) <doi:10.64898/2026.08.26.26361466>.

- [immGLIPH](/packages/immGLIPH) An R implementation of the GLIPH and
  GLIPH2 algorithms for clustering T cell receptors (TCRs) predicted
  to bind the same HLA-restricted peptide antigen. Identifies
  specificity groups based on local (motif-based) and global
  (sequence-based) CDR3 similarities. Integrates with the
  scRepertoire ecosystem via immApex for single-cell immune
  repertoire analysis. Users should cite the original GLIPH algorithm
  papers: Glanville et al. (2017) <doi:10.1038/nature22976> and Huang
  et al. (2020) <doi:10.1038/s41587-020-0505-4>.

- [inferRecom](/packages/inferRecom) This package implements a
  nuclear family-based algorithm to detect homologous recombination.
  Using genotypes in PLINK binary format, recombination events are
  inferred from changes in the inheritance state of informative SNPs
  between siblings. Closely spaced crossovers that are unlikely to
  represent true recombination can be filtered out using genetic maps
  or a SNP-count filter. Additional functions are supplied to detect
  runs of homozygosity and to phase haplotypes at informative SNPs
  and intervals that are inferred identical-by-descent. Recombination
  intervals and runs of homozygosity are returned as `GRanges`
  objects, and phased genotypes are returned as a `VariantAnnotation`
  `VCF` object, for use with other Bioconductor packages and genome
  browsers.

- [IntegratedLearner](/packages/IntegratedLearner) Provides a unified
  interface for multi-omics prediction using early and late fusion
  for continuous, binary, multiclass, and survival outcomes. It
  supports both MultiAssayExperiment and PCL-style inputs, performs
  input validation and feature/sample harmonization across layers,
  and provides model fitting, prediction, plotting, and
  variable-importance utilities.

- [IntrinsicScan](/packages/IntrinsicScan) Scans continuous molecular
  features for intrinsic binary partitions without using phenotype
  labels to determine the split. For each feature, observations are
  divided into two one-dimensional states by selecting the split that
  minimizes the total within-state L1 loss, defined here as the sum
  of absolute distances between each observation and the medoid of
  its assigned state. Thus, unlike a squared-error objective,
  deviations are measured according to their absolute magnitude
  rather than being squared. Binary phenotype labels are used only
  after the partition is fixed to orient the two states, quantify
  phenotype concordance using orientation-free Youden J, perform
  exact conditional inference, and optionally assess resampling
  stability or fixed-partition permutation calibration. The method is
  intended for preprocessed continuous omics measurements and does
  not imply that a fitted partition represents true biological
  bimodality or discrete biological states.

- [iSEEid](/packages/iSEEid) This package provides a dedicated panel
  to assist the sample identification within a single iSEE
  web-application. The main functionality of this package can be
  adopted to efficiently select samples (e.g. cells) that need to be
  identified and annotated in the analysis workflow. With iSEEid, the
  iSEE framework fully qualifies to be an ideal companion for
  assigning appropriate labels to large omics datasets.

- [IsoNetAnalyzeR](/packages/IsoNetAnalyzeR) Compares predicted
  protein-interaction neighbourhoods of transcript isoforms. The user
  supplies isoform seeds labelled "up" or "down", for example the
  isoforms that gain or lose usage in a differential transcript usage
  analysis. Each seed set is expanded in a bundled isoform-resolved
  predicted interaction network, or in a user-supplied igraph
  network, and the resulting neighbourhoods are compared per gene and
  summarised across genes. The package also provides
  connected-component analysis, network diagnostics, plots, and
  over-representation and gene set enrichment analysis of
  neighbourhood statistics. Protein-identical transcripts are
  resolved to shared representatives, and nonsense-mediated and
  non-stop decay transcripts, annotated from Ensembl through
  AnnotationHub, are excluded by default. Results are S4Vectors
  DataFrames, and gene sets can be GSEABase collections such as those
  of the msigdb package. Statistical testing of isoform usage remains
  upstream of the package.

- [karioCaS](/packages/karioCaS) Provides a comprehensive framework
  for analyzing Kraken2 metagenomic classification using multiple
  confidence scores. karioCaS implements a robust comparative
  approach to evaluate taxonomic stability across different
  stringency levels, facilitating the identification of reliable
  microbiota components and the discovery of potentially novel or
  unrepresented species. A key innovation is its domain-specific
  analysis, separately processing Bacteria, Archaea, Eukaryota, and
  Viruses to ensure balanced visibility of all biological domains.
  Includes high-quality visualization tools for retention rates,
  shared taxa (Upset plots), and stability heatmaps.

- [KEGGemUP](/packages/KEGGemUP) KEGGemUP retrieves and renders
  interactively KEGG pathway graphs. The retrieval makes full use of
  the caching functionality to avoid unnecessary download processes.
  The interactivity of the graph is warranted via the visNetwork
  interface to the vis.js library, fully supporting bindings to be
  used in Shiny for further operations.

- [lantern](/packages/lantern) Ancestry-specific rare variant burden
  testing for admixed populations. Splits genotype dosages by local
  ancestry (as inferred by 'RFMix') into per-ancestry dosages, using
  either phased haplotypes or unphased genotypes, where ambiguous
  mixed-ancestry heterozygotes are split with shrinkage toward 1/2,
  and writes them to 'SeqArray' GDS files. Ancestry-stratified
  gene-level tests are then run with the 'SMMAT' framework of Chen et
  al. (2019) <doi:10.1016/j.ajhg.2018.12.012> as implemented in
  'GMMAT', and the per-ancestry p-values are combined with the Cauchy
  combination test of Liu and Xie (2020)
  <doi:10.1080/01621459.2018.1554485>.

- [LIPIDIFy](/packages/LIPIDIFy) Provides a comprehensive toolkit for
  end-to-end lipidomics data analysis, including missing value
  imputation, batch effect correction, normalization, differential
  abundance analysis using limma and edgeR, gene set enrichment
  analysis, and extensive visualization capabilities. Lipid names are
  automatically classified by class, subclass, and fatty-acid
  saturation. Features both an interactive Shiny interface for bench
  biologists and fully scriptable R functions for bioinformaticians.
  Supports flexible custom lipid classification schemes and
  user-defined enrichment sets.

- [MagellanNTK](/packages/MagellanNTK) MagellanNTK (Magellan
  Navigation ToolKit) is a collection of Shiny modules which provide
  the basis for a generic workflow manager for data processing tools
  with several steps. It works in combination with any package which
  implements data processing pipelines compliant with MagellanNTK
  such as Prostar2.

- [mesa](/packages/mesa) A package for the analysis of methylation
  enrichment sequencing data (e.g. MBD-seq or MEDIP-seq). This allows
  for the window-based evaluation of methylation levels (using the
  'qsea' package), including functions for determining differentially
  methylated regions. Many functions are provided for tidyverse style
  modification of the qseaSet objects originally defined in 'qsea'.

- [MetaPathNet](/packages/MetaPathNet) Provides tools to construct
  KEGG-based metabolic and signaling networks as edge lists for
  single-organism or cross-species analyses. The package supports
  identifier mapping, shortest-path and topology analyses, community
  detection, permutation testing, pathway over-representation
  analysis, and node annotation for host-microbiome studies. It also
  provides network visualisation in R and Cytoscape and supports
  extension of KEGG-based networks through additional reaction
  resources and user-defined reactions.

- [methylTFR](/packages/methylTFR) methylTFR is an extensive R
  package designed specifically for analyzing DNA methylation
  signatures within transcription factor binding sites (TFBS) across
  diverse genomic annotations. It enhances methylation analysis by
  calculating a bias-corrected 'deviation' in methylation levels,
  leveraging GC bias distribution and frequency data.  This package
  offers a variety of functions that facilitate the motif activity
  exhibiting diverse methylation patterns across distinct samples or
  cells.

- [mgcca](/packages/mgcca) Implements Generalized Canonical
  Correlation Analysis (GCCA) with missing individuals (whole missing
  rows; van de Velden and Bijmolt, 2006) for multi-block omics
  integration. The package operates directly on HDF5 matrices via the
  BigDataStatMeth C++ API to scale computations; R functions are thin
  wrappers around the C++ core.

- [MiCARA](/packages/MiCARA) Provides a comprehensive framework for
  multi-cohort microbiome differential abundance and residual
  correlation analysis. It integrates confounder-adjusted linear
  models with centered log-ratio transformed taxonomic and functional
  pathway profiles.  The package features robust preprocessing,
  automated covariate handling, multiple-testing correction, and
  quality control diagnostics for multi-omics data.  Interactive
  visualization tools are included to facilitate downstream
  interpretation of cross-domain interaction networks.

- [MRCCC](/packages/MRCCC) Implements MR-CCC, a Bayesian Mendelian
  randomization framework for inferring causal ligand-receptor
  signalling between cell types from population-scale single-cell
  expression and genotype data. For each ligand-receptor-pathway
  triplet, cis-eQTL genotypes instrument ligand expression in the
  sender cell type and receptor expression in the receiver cell type.
  A working model with a ligand-by-receptor interaction is then
  fitted to downstream pathway activity by a blocked Gibbs sampler. A
  spike-and-slab prior on the causal block yields a posterior
  inclusion probability that serves as the communication score. The
  retained draws support credible intervals, sign probabilities, the
  receptor level at which the ligand effect reverses sign, and
  Bayesian false discovery rate control across many triplets. Inputs
  are taken from SummarizedExperiment objects of donor-level
  expression and genotype dosages, and cis-eQTL instruments can be
  selected within the package.

- [MRManalyzeR](/packages/MRManalyzeR) Reproducible post-acquisition
  processing of targeted LC-MS/MS lipidomics and metabolomics results
  exported from Waters TargetLynx or from Skyline, or supplied as a
  plain sample-by-analyte matrix. Reads the quantitative result
  tables, integrates them with sample and feature metadata, and
  applies configurable signal filtering (SNR, LOD or LOQ), blank
  filtering, batch correction, normalisation, internal-standard and
  volume adjustment of vendor-reported concentrations, and
  missing-value imputation. Each step is an exported function acting
  on a struct DatasetExperiment; a complete workflow can additionally
  be driven from a single YAML configuration via runMRManalyzeR(),
  which also renders self-contained HTML data-quality and statistical
  reports. Chromatographic peak detection, peak integration and
  calibration-curve fitting from raw mass-spectrometry data are
  outside its scope.

- [MsBackendMassIVE](/packages/MsBackendMassIVE) MassIVE is one of
  the main public repositories for storage of metabolomics
  experiments. The MsBackendMassIVE package provides functionality to
  retrieve and represent mass spectrometry (MS) data from MassIVE.
  Data files are downloaded and cached locally avoiding repetitive
  downloads. MS data from metabolomics experiments can thus be
  directly and seamlessly integrated into R-based analysis workflows
  with the Spectra and MsBackendMassIVE package.

- [MsBackendMetabolomicsWorkbench](/packages/MsBackendMetabolomicsWorkbench)
  Metabolomics Workbench is one of the main public repositories for
  storage of metabolomics experiments. The
  MsBackendMetabolomicsWorkbench package provides functionality to
  retrieve and represent mass spectrometry (MS) data from
  Metabolomics Workbench. Data files are downloaded and cached
  locally avoiding repetitive downloads. MS data from metabolomics
  experiments can thus be directly and seamlessly integrated into
  R-based analysis workflows with the Spectra and
  MsBackendMetabolomicsWorkbench package.

- [MsBackendParquet](/packages/MsBackendParquet) Mass spectrometry
  (MS) data backend for the Spectra Bioconductor package that stores
  and retrieves spectra and peak data from on-disk Apache Parquet
  datasets. Spectra metadata together with m/z and intensity values
  are stored in a columnar Parquet format that supports compression,
  dictionary encoding and Hive-style partitioning. Datasets are
  written with the Apache Arrow R bindings and queried with DuckDB.
  Spectra variables are held in a bounded in-memory index so that
  filtering and metadata access are fast while peak data remains on
  disk.

- [MsExperimentStash](/packages/MsExperimentStash) The serialization
  mechanism of R allows to save and load R data objects in a binary
  format, that can however not be read by other programming languages
  or software. The MsExperimentStash package implements the
  serialization methods from the MsStash package for MsExperiment
  objects. Data can be stashed in a variety of file formats and
  types, including HDF5 and JSON-based formats defined by the
  Bioconductor *alabaster* package. The package provides also
  functionality to read a full MS experiment, including MS data and
  experimental metadata information, from the MetaboLights
  repository.

- [MsStash](/packages/MsStash) The serialization mechanism of R
  allows to save and load R data objects in a binary format, that can
  however not be read by other programming languages or software. The
  MsStash package defines classes and method to enable serializing
  and restoring, or importing, mass spectrometry data objects to and
  from language agnostic file formats. A variety of different file
  types, including HDF5 and JSON-based formats defined by the
  Bioconductor *alabaster* package are supported. The file type can
  be defined and configured through a second argument `param` of the
  save/read methods. Additional packages extending MsStash provide
  the implementations to serialize/restore MS data objects.

- [MSTree](/packages/MSTree) This package is used to generate a graph
  object from the output of chewBBACA pipeline
  (https://chewbbaca.readthedocs.io/en/latest/). Then, the generated
  graph object can be used to make a minimum spanning tree (MST). The
  minimum spanning tree can be customized using all the available
  arguments. This package consists of two functions: one to build the
  graph and another one for plotting.

- [MultiAssaySpatialExperiment](/packages/MultiAssaySpatialExperiment)
  Extends MultiAssayExperiment with spatial context elements (images,
  labels, points, shapes) for integrative analysis of spatially
  resolved multi-omics data. A spatial map links the columns of each
  assay to the instances of the spatial elements. Supports coercion
  to and from SpatialExperiment, SpatialFeatureExperiment and
  spatialdataR's SpatialData, and provides readers for Visium, Visium
  HD, Xenium, CosMx and MERSCOPE data.

- [MultiOmicsBridge](/packages/MultiOmicsBridge) MultiOmicsBridge
  provides an end-to-end, reproducible computational framework for
  integrative analysis of paired host transcriptomics (bulk RNA-seq)
  and gut microbiome (16S rRNA or shotgun metagenomics) data. The
  package addresses the lack of a unified Bioconductor workflow for
  this pairing by implementing five modules: (1) data harmonization
  and normalization with CLR transformation for microbiome
  compositional data and TMM/voom for RNA-seq; (2) joint
  dimensionality reduction via sparse multi-block PLS-DA (DIABLO);
  (3) multi-omics biomarker discovery through cross-omics correlation
  networks and sparse feature loadings; (4) integrated diagnostic
  classification comparing host-only, microbiome-only, and joint
  Random Forest models with stratified cross-validation; and (5)
  publication-quality visualization of integration results, biomarker
  networks, classifier comparisons, and feature flow diagrams. All
  functions operate natively on SummarizedExperiment and
  MultiAssayExperiment objects and return a structured MOBResult S4
  object. The package is validated on inflammatory bowel disease
  multi-omics data and designed with complex disease contexts
  (tuberculosis, HIV, EED) in mind.

- [MultipleAlignment](/packages/MultipleAlignment) The package
  implements a set of S4 classes (DNAMultipleAlignment,
  RNAMultipleAlignment, AAMultipleAlignment) for representing
  Multiple Sequence Alignments (MSA). The classes allow users to
  represent groups of aligned DNA, RNA or amino acid sequences as a
  single object. The package also provides functions to read/write
  such object from/to traditional MSA file formats including
  Stockholm and Clustal.

- [multipletR](/packages/multipletR) Detects human-mouse multiplets
  in patient-derived xenograft (PDX) single-cell RNA-seq data using
  an adaptive threshold method that does not assume a fixed species
  proportion. Takes a 10x CellRanger GEM classification file and
  returns the data with added multiplet classifications, with
  optional diagnostic plots and helpers to annotate a Seurat or
  SingleCellExperiment object with multiplet classifications or
  remove multiplets from it. The adaptive approach starts from a
  conservative central region of balanced droplets and expands its
  thresholds until the selected cells no longer resemble true
  multiplets, which handles the imbalanced species mixtures common in
  real PDX samples.

- [multipointR](/packages/multipointR) `multipointR` is a package to
  compare the distribution of cells in an image or cross images with
  point process models. On a single image level point process models
  (`ppm`) model the spatial distribution of a cell type point pattern
  as a function of spatial covariates while accounting for natural
  spacing of cells. The main model class considered in `multipointR`
  are inhomgoeneous Gibb's point process models. Across multiple
  images, users can either compare multiple univariate `ppm` models
  in a for loop or fit one joint model across all images with `mppm`.
  `multipointR` provides an interface between `SpatialExperiment` and
  `SpatialFeatureExperiment` objects and let's users flexibly define
  their own `ppm`/`mppm` models with R's formula interface.

- [NADIA](/packages/NADIA) Provides an end-to-end workflow for
  differential abundance analysis of protein quantification data,
  with particular attention to missing values in data-independent
  acquisition (DIA) proteomics. Quantification reports exported from
  Spectronaut, DIA-NN and Proteome Discoverer for TMT and label-free
  experiments (DDA/DIA) are converted into a common data structure,
  normalised, optionally batch-corrected, imputed and analysed using
  limma or limpa. The package includes individual and two-stage
  imputation strategies that treat values assumed to be missing at
  random and missing not at random separately. Data-driven functions
  evaluate and rank normalisation and imputation workflows using
  simulated missingness or spike-in datasets with known ground truth,
  allowing method selection to be justified for each dataset. Results
  can be explored through interactive Highcharts visualisations,
  static ggplot2 and ComplexHeatmap figures, and sortable tables.

- [normScore](/packages/normScore) Provides tools to evaluate and
  rank normalization methods for omics datasets using a composite
  score derived from multiple performance metrics. The package is
  designed to support systematic benchmarking and comparison of
  normalization strategies across datasets and experimental settings.
  It also includes utilities for summarizing results and visualizing
  normalization performance.

- [OmniCorr](/packages/OmniCorr) Provides tools for correlation-based
  integration and visualization of multi-omics datasets, including
  transcriptomics, metatranscriptomics, metagenomics, metaproteomics,
  and associated sample metadata. The package computes pairwise
  correlations between independently processed omics layers and
  provides flexible tools for organizing and visualizing cross-omics
  associations using aligned heatmaps and optional network-style
  visualizations. OmniCorr is designed to complement existing
  preprocessing and network inference workflows rather than replace
  them.

- [ontoProc2](/packages/ontoProc2) This package provides ontology
  facilities based on INCAtools Semantic SQL.  Tooling is provided to
  retrieve and cache SQLite databases representing curated
  ontologies. ontologyIndex ontology_index instances can be produced.
  Anthropic's Claude was used in the development of S7 classes and
  methods and in the production of documentation.

- [orthoRibbon](/packages/orthoRibbon) OrthoRibbon is a visualization
  toolkit for displaying orthologous gene relationships across
  species using ribbon-based diagrams. It enables intuitive
  comparison of conserved genes and genomic correspondence between
  chromosomes or genomic regions. The package supports flexible input
  formats and produces publication-quality figures for comparative
  genomics and evolutionary analysis.

- [PACE](/packages/PACE) PACE (Proximity-Associated Changes in
  Expression) is a hierarchical empirical Bayes framework that
  quantifies how a cell's gene expression changes with proximity to
  specific neighbouring cell types in imaging-based spatial
  transcriptomics such as Xenium and CosMx. It fits negative binomial
  mixed models with partial pooling across cell types, corrects for
  transcript contamination between adjacent cells with a per-cell
  ambient term, and decomposes expression variance into cell-type
  identity, spatial cell state, contamination, and residual
  components. Input and output follow SpatialExperiment conventions.

- [panoramic](/packages/panoramic) Quantifies and compares cell-type
  spatial colocalization across samples, patients, and conditions in
  spatial omics studies. Its primary statistic, edge-corrected local
  competition enrichment, measures excess or deficit local mixing for
  ordered cell-type pairs on a percentage-point scale. The package
  estimates within-sample uncertainty with spatial bootstrap
  procedures and propagates it through multilevel random-effects
  meta-analysis for pooled inference and group-level contrasts;
  classical spatial summary statistics are also supported.

- [panpathAudit](/packages/panpathAudit) Audits source FASTA
  sequences against named paths and coordinate-aware walks in GFA
  files. It reports exact recovery, divergence coordinates, digests,
  and optional base-level statistics. Graph topology is not
  validated.

- [path2path](/packages/path2path) Builds a directional
  pathway-by-pathway interaction map from genome-scale Perturb-seq
  differential-expression statistics by two stacked runs of preranked
  gene set enrichment analysis: once per perturbation over genes,
  then once per readout pathway over perturbations. Entry (j,k) of
  the resulting matrix measures whether perturbations targeting
  pathway j's members produce a coordinated transcriptional response
  in pathway k. The package includes target exclusion, effective-size
  filters, overlap masking, Stage 2 p-values with
  false-discovery-rate control, and a whole-pipeline permutation
  calibration with generalized-Pareto tail estimation.

- [pepVet](/packages/pepVet) Simulates proteolytic digestion, scores
  the resulting peptides for LC-MS/MS suitability, compares candidate
  enzymes, and reports digest quality at the protein level. Supports
  40 cleaver-compatible enzyme rules, workflow presets, peptide mass
  and pI calculations, sequence-local cleavage-efficiency
  annotations, and proteome-aware uniqueness scoring. Evaluates
  multi-FASTA files in batches with per-protein triage and
  proteome-level summaries. Exports peptide lists for Skyline and
  generic downstream tools and prints styled console reports.

- [pgen2gds](/packages/pgen2gds) Provides functions to convert files
  from the PLINK2 pgen format to SeqArray GDS, and from SeqArray GDS
  to the PLINK2 pgen format.

- [phylobar](/packages/phylobar) The phylobar package supports
  interactive visualization of microbiome data by allowing a stacked
  barplot to be constructed using a guiding taxonomic or phylogenetic
  hierarchy. The package provides a strategy for collapsing and
  expanding the hierarchy to different levels of resolution and then
  for interactively "painting" the stacked barplot by placing the
  mouse over different subtrees. This makes it possible to
  interactively test different color palettes at different resolution
  and identify taxonomic groups with interesting variation before
  settling on a final stacked barplot. One advantage of the approach
  is that multiple levels of taxonomic resolution can be compared at
  once within the same view.

- [PinPath](/packages/PinPath) PinPath enables flexible visualization
  of (omics) data onto pathways diagrams, allowing users to pinpoint
  where the relevant changes occur. It supports pathway diagrams from
  WikiPathways and KEGG, as well as custom GPML and KGML files. Data
  can be displayed on both native pathway layouts and network
  representations.

- [polyICSFlow](/packages/polyICSFlow) polyICSFlow systematically
  identifies all cytokine combinations detected in an intracellular
  cytokine staining (ICS) flow cytometry assay and quantifies
  polyfunctional antigen-specific responses to single or multiple
  antigens. The package also offers versatile plotting options for
  efficient data exploration and visualization of cytokine
  co-expression profiles.

- [PostChicago](/packages/PostChicago) PostChicago is a toolbox for
  visualizing and assessing the output from the CHiCAGO pipeline. The
  plots created by PostChicago show reads or CHiCAGO scores over
  different regions. PostChicago can integrate different experiments
  with other types of datasets and compare separate conditions.
  Detailed description can be found in Krall et al. (2026)
  <doi:10.64898/2026.10.08.757510>.

- [prolfqua](/packages/prolfqua) Covers the steps of a mass
  spectrometry based label-free quantification (LFQ) analysis:
  quality control, normalization, aggregation of peptides to
  proteins, modelling, hypothesis testing and sample size estimation.
  A single function, build_contrast_analysis(), fits a model per
  protein and computes contrasts with one of many methods: linear,
  robust and rank-based regression, 'limma' and 'limma'-voom, 'DEqMS'
  variance moderation, the detection probability model of 'limpa',
  Firth logistic and binomial models of detection, mixed models and
  ROPECA on peptide-level data, a hurdle model combining intensity
  and detection, and variants that impute groups without observations
  at the limit of detection. All methods return their results through
  the same interface and columns, so they can be compared and
  combined. Models simple designs with one explanatory variable as
  well as experiments with several factors. Developed while analysing
  dozens of proteomics datasets at the Functional Genomics Center
  Zurich.

- [ProteinBatcher](/packages/ProteinBatcher) ProteinBatcher provides
  utilities to read, validate, and harmonize quantitative proteomics
  result tables generated by DIA-NN from data-independent acquisition
  mass spectrometry experiments. The package enables batch processing
  of multiple DIA-NN output files, consistent handling of protein
  identifiers and associated metadata, and transformation of tabular
  data into analysis-ready formats. It facilitates downstream
  differential abundance analysis using linear modeling frameworks
  based on limma, with explicit support for multivariable
  experimental designs including batch effects and other confounding
  covariates, as well as correlation-aware designs via
  duplicateCorrelation.

- [PscanR](/packages/PscanR) R implementation of the 'Pscan'
  algorithm for motif enrichment analysis in regulatory DNA
  sequences, based on Zambelli et al. (2009)
  <doi:10.1093/nar/gkp464>. 'PscanR' scans a set of DNA sequences
  (e.g., promoters) from co-regulated or co-expressed genes, using
  transcription factor binding motifs from databases such as JASPAR.
  It identifies significantly enriched motifs, providing thus hints
  on which transcription factors could be common regulators of the
  genes studied, together with the location of their candidate
  binding sites in the sequences.

- [punKEGGer](/packages/punKEGGer) Provides tools to parse KEGG
  pathway KGML files into tidy graphs, expand metanodes (e.g.,
  multi-gene entries), annotate them with custom dictionaries, and
  produce consistent visualizations using ggraph. Enables integration
  of KEGG pathways with tidygraph-based workflows.

- [QFeaturesGUI](/packages/QFeaturesGUI) QFeaturesGUI is a suite of
  shiny apps that serve as graphical interfaces for the QFeatures
  package. The package provides three apps: import, process, and
  visualise.

- [quantMSImageR](/packages/quantMSImageR) Implements tools for
  processing and quantifying targeted DESI-MRM mass spectrometry
  imaging (MSI) datasets, extending the Cardinal package. Includes
  signal-to-noise filtering against background pixels,
  tissue/background separation, per-feature ion images, quantile
  heatmaps, quantification against on-slide or on-tissue calibration
  standards, and batch export of per-feature text images for external
  viewers.

- [RAMEN](/packages/RAMEN) RAMEN is an R package which goal is to
  identify Variable Methylated Loci (VML) in microarray DNA
  methylation data. Then, using genome and exposome data, RAMEN
  identifies which individual genetic (G) or environmental (E) or
  joint (G+E or GxE) model better explains the VML variability. This
  package provides estimates of genetic and environmental
  contribution for each VML, allowing researchers to gain insights
  into the complex interplay between genetics, environment and DNA
  methylation.

- [RBPEqBind](/packages/RBPEqBind) Simulates competitive binding of N
  RBPs to RNA sequences using an equilibrium binding model. Supports
  single sequences, concentration grids, and transcriptome-wide
  analysis. Designed for Bioconductor.

- [RBPSpecificity](/packages/RBPSpecificity) Provides tools to
  analyze RNA Binding Protein (RBP) binding specificities from
  high-throughput sequencing data. Functions include calculating
  K-mer enrichment, Inherent Specificity (IS), and Mutational
  Sensitivity (VS), along with visualization of IS and VS. For
  detailed methodology and applications, please refer to our
  manuscript (see URL field or vignette).

- [reglScatterplotR](/packages/reglScatterplotR) Lightweight
  'htmlwidgets' interface to the JavaScript 'regl-scatterplot'
  library, rendering millions of two-dimensional points in the
  browser via WebGL. Designed for exploratory visualisation of
  single-cell, spatial transcriptomics and other high-dimensional
  biological data. Features include synchronised pan/zoom across
  multiple plots, categorical and continuous colour mapping, lasso
  selection, range and category filtering, PNG/SVG/PDF export, and
  Shiny integration. The widget works in the 'RStudio' Viewer,
  standalone HTML files, Shiny applications and Jupyter notebooks
  running the 'IRkernel'.

- [RhyME](/packages/RhyME) Fits beta cosinor models to methylation
  measurements using penalized quasi-likelihood (PQL), an optional
  subject random intercept, adjustment covariates, and an optional
  empirical Bayes (EB) prior on log precision. Tests rhythmicity with
  a null-model score test and estimates rhythm parameters from a
  separate full-model fit.

- [RmzTabM](/packages/RmzTabM) The RmzTabM package provides functions
  to read and write files in the HUPO-PSI mzTab-M file format for
  reporting of results from metabolomics and lipidomics experiments.
  This package represents the reference implementation for mzTab-M
  import/export for R. Core, low level functions allow to format and
  parse the individual sections of the mzTab-M format. These
  functions are designed to be re-used by other packages implementing
  import/export of results in mzTab-M.

- [rnaSentry](/packages/rnaSentry) Orchestrates the statistical
  pipeline underlying prognostic biomarker discovery from bulk
  RNA-seq SummarizedExperiment objects: sample-identity and
  quality-control auditing, confounder and design-formula detection
  prior to differential expression, signature construction with
  explicit cross-validation stability reporting, survival modeling
  (Kaplan-Meier, log-rank, Cox proportional hazards), and external
  validation against public cohorts. Every automated statistical
  decision is reported with its supporting statistic and a
  plain-language rationale, and validation stages are gated rather
  than silently degraded when guardrail criteria are not met. The
  package records all guardrail decisions, flag ledgers, and
  validation results in a structured audit trail that can be rendered
  into a standalone HTML report.

- [RPBC](/packages/RPBC) Implements Region-Preserving Batch
  Correction (RPBC), an empirical-Bayes location-scale framework for
  paired region-resolved single-cell and spatial-omics data. RPBC
  preserves a patient-balanced population-average regional contrast
  while harmonizing patient-wide location shifts, patient-specific
  regional deviations, and patient-region scale heterogeneity. The
  package also provides Adaptive RPBC for group-size-dependent
  stabilization of patient-region scale variation and
  covariate-adjusted RPBC for preserving prespecified patient-level
  covariate-associated variation during harmonization. RPBC is
  designed to interoperate with Bioconductor data structures for
  single-cell and spatial experiments.

- [S4Cartographer](/packages/S4Cartographer) S4 Object-Oriented
  Programming is fundamental to the Bioconductor ecosystem, enabling
  modularity and code reuse across hundreds of packages. However, the
  multi-layered inheritance structures within and between packages
  can be difficult to navigate. S4Cartographer facilitates the
  exploration of these relationships by visualizing S4 class
  hierarchies as directed graphs. By representing classes as nodes
  and inheritance as edges, it provides a clear overview of complex
  dependency networks.

- [sbivar](/packages/sbivar) Sbivar implements a suite of tests for
  Spatial BIVARiate association across omics modalities, with
  possibly disjoint coordinate sets. Implemented tests are
  generalized additive models (GAMs), modified t-test, bivariate
  Moran's I and Gaussian processes (GPs). Both single images and
  replicated experiments can be analysed.

- [scBatchQC](/packages/scBatchQC) scBatchQC provides a hierarchical
  empirical Bayes framework for quality control in multi-sample,
  multi-batch single-cell RNA-seq experiments. Unlike per-sample QC
  tools, scBatchQC jointly models QC metric distributions (library
  size, gene count, mitochondrial fraction) and doublet rates across
  batches, enabling calibrated cell-level QC calls that account for
  batch structure. The package operates natively on
  SingleCellExperiment objects and returns augmented colData with
  per-cell QC flags and batch-adjusted doublet scores.

- [scCertify](/packages/scCertify) scCertify provides an explainable
  framework for assessing the confidence of cell-type annotations in
  single-cell RNA sequencing data. Confidence scores are computed
  using marker enrichment, neighborhood agreement, entropy-based
  uncertainty estimation, and ontology-aware label matching to
  support transparent interpretation of annotation quality.

- [scCompoundDE](/packages/scCompoundDE) scCompoundDE decomposes
  pseudo-bulk differential expression (DE) signals into two
  orthogonal components: transcriptional changes (cell-intrinsic
  expression shifts) and compositional changes (shifts in the
  relative abundance of cell subtypes). Standard pseudo-bulk DE tools
  confound these two sources of signal, producing spurious DE calls
  when subtype proportions differ between conditions. scCompoundDE
  fits per-subtype limma-voom models, estimates subtype proportion
  shifts, and uses a z-score-normalized decomposition to assign each
  gene a TC_ratio score — the fraction of its DE signal attributable
  to transcription versus composition. Genes are then classified as
  transcriptional (real biology), compositional (artifact), or mixed
  (requires caution). All functions operate natively on
  SingleCellExperiment objects and return a CDEResult S4 object that
  extends the standard DE output with full decomposition statistics.

- [scFastDE](/packages/scFastDE) scFastDE provides fast,
  donor-weighted pseudo-bulk differential expression analysis for
  multi-donor single-cell RNA-seq experiments. Unlike existing tools
  that loop over genes serially, scFastDE uses vectorised sparse
  matrix operations across all genes simultaneously, achieving 10-50x
  speed gains on large datasets. Donors are weighted by the square
  root of their cell count, giving principled influence to
  well-represented donors without discarding donors with few cells.
  Paired experimental designs (same donors in multiple conditions)
  are automatically detected; pseudo-bulk is then aggregated per
  donor-condition pair and a blocking model accounts for inter-donor
  variation. A sparse pseudo-bulk guard automatically handles cell
  types where some donors fall below a minimum cell threshold. All
  functions operate natively on SingleCellExperiment objects.

- [scOverlay](/packages/scOverlay) scOverlay provides
  publication-grade plotting helpers to display a foreground feature
  (gene expression or metadata) on top of a contextual background
  (clusters, cell types, samples, or a solid colour) in the same
  reduced dimension panel, each with independent colour scales (via
  ggnewscale). It integrates with SingleCellExperiment /
  SummarizedExperiment, supports flexible subsetting of foreground
  cells, multiple assays (e.g. counts, logcounts), and grid
  arrangements by user-defined grouping columns.

- [scProkaR](/packages/scProkaR) scProkaR provides a
  SingleCellExperiment-centred workflow for microbial single-cell RNA
  sequencing data. It standardises bacterial input objects from count
  matrices, 10x Genomics directories or Seurat objects, and computes
  prokaryote-specific quality-control metrics such as ribosomal RNA
  fraction. It performs batch integration with several backends and
  scores the resulting embeddings with scIB-inspired benchmarking
  metrics. It also supports pseudobulk differential expression with
  edgeR and Time Aware Trajectory Analysis (TATA), a time-informed
  trajectory abstraction for perturbation time courses.

- [selexprepR](/packages/selexprepR) Provides reproducible
  preprocessing for high-throughput SELEX sequencing data. Supports
  primer inference, extraction of random regions, sparse sequence
  counting across rounds, quality control, and provenance manifests.
  Public-study discovery and checksum-aware ENA FASTQ download
  planning are included for reproducible data acquisition.

- [seqRdna](/packages/seqRdna) Functions for working with rDNA data
  generated from sequencing experiments. The package includes
  functions to calculate the relative rDNA copy number from coverage
  data, as well as functions to calculate the absolute rDNA copy
  number using methylation data.

- [sigvar](/packages/sigvar) This package allows users to import
  mutational signature attribution (a.k.a. exposure) matrices and
  compute, visualize, and test their variabilities within and across
  samples.

- [SimiCviz](/packages/SimiCviz) Visualization and export utilities
  for SimiC and SimiCPipeline outputs. The package focuses on
  importing SimiC-style results (e.g., weights, AUC metrics) from
  pickle/CSV files, processing and generating publication-ready plots
  (networks, heatmaps, distributions) and tables saved into a
  reproducible, ordered directory hierarchy.

- [spammR](/packages/spammR) SpammR provides a suite of tools and a
  data structure to analyze diverse multiomics measurements in a
  spatial context. It enables the analysis and visualization of omics
  measurements paired with a specific image, each with spatial
  coordinates.

- [SpatialData.data](/packages/SpatialData.data) A companion package
  to access scverse-curated example SpatialData stores. The datasets
  span spatial omics technologies such as Visium, Visium HD, Xenium,
  MERFISH, MIBI-TOF, and imaging mass cytometry, and are distributed
  as zipped Zarr stores through the Bioconductor bucket on the NSF
  Open Storage Network (OSN) or the scverse spatialdata sandbox.
  Available datasets and their metadata can be listed, and each
  dataset is downloaded once, cached locally with 'BiocFileCache',
  and read into R as a SpatialData object via 'spatialdataR' package.

- [SpatialData.plot](/packages/SpatialData.plot) Visualization suit
  for 'spatialdataR'. Current functionality includes handling of
  multiscale 'images', visualizing 'labels', 'points', and 'shapes'.
  For the latter, POINT, POLYGON, and MULTIPOLYGON geometries are
  supported. For all but 'images', there is support to pass
  information from annotation 'tables' as aesthetics (e.g., 'colData'
  or 'assay' data).

- [spatialdataR](/packages/spatialdataR) R interface to
  Python/scverse's 'spatialdata' framework for unified spatial omics
  data handling. Adheres to OME-NGFF standards, providing lazy,
  on-disk representations for multiscale images and labels
  (ZarrArray), as well as points and shapes (DuckDB-backed tables).
  Includes handling of coordinate transformation systems and spatial
  utilities for cropping, masking, and querying. Integrates tabular
  annotations as 'SingleCellExperiment's (via 'anndataR'), enabling
  interoperable spatial omics workflows across R and Python.

- [SpectraStash](/packages/SpectraStash) The serialization mechanism
  of R allows to save and load R data objects in a binary format,
  that can however not be read by other programming languages or
  software. The SpectraStash package implements the serialization
  methods from the MsStash package for Spectra objects and its
  MsBackend classes. A variety of different file formats and types,
  including HDF5 and JSON-based formats defined by the Bioconductor
  *alabaster* package are supported. The file type can be defined and
  configured through a second argument `param` of the save/read
  methods.

- [strainspy](/packages/strainspy) Fits statistical models to detect
  microbial strains associated with a phenotype. For small sample
  sizes, models can be fitted using empirical Bayes maximum a
  posteriori (MAP) priors to improve detection power. Supports
  various experimental design scenarios, including repeated measures.

- [STRprofilerR](/packages/STRprofilerR) Compares short tandem repeat
  (STR) profiles to authenticate biomedical models such as cell lines
  and xenografts against the primary tissue samples they were derived
  from. Provides the Tanabe (Sorenson-Dice), Masters (versus query),
  and Masters (versus reference) similarity scores, heuristics for
  detecting sample mixing, flexible ingest of wide and long profile
  formats, comparison against user-supplied reference databases, and
  querying of the Cellosaurus knowledge base through the CLASTR REST
  API. A command line application is included for use in automated
  workflows, and a Shiny application for interactive queries.

- [SummarizedHeatmap](/packages/SummarizedHeatmap) Represents a
  heatmap as a SummarizedExperiment-derived object, holding a single
  assay together with its row and column annotations and any
  hierarchical clustering. Assay values, annotations, display order,
  and dendrograms stay synchronized as the object is subsetted or
  reclustered, so a heatmap's components cannot drift out of step
  with each other. Each visual component (the main heatmap,
  annotation strips, dendrograms) is returned as an ordinary ggplot2
  plot, and these can be composed, laid out, and customized using
  patchwork.

- [SwarnSeq](/packages/SwarnSeq) This R package performs differential
  expression and differential zero inflation analysis of single-cell
  RNA-seq (scRNA-seq) UMI counts data through adjusting cell capture
  efficiency.

- [thyroidBRS](/packages/thyroidBRS) Computes the BRAF-RAS Score
  (BRS), a 71-gene expression-based signature that quantifies whether
  a papillary thyroid carcinoma (PTC) sample resembles a
  BRAF-V600E-mutant or a RAS-mutant transcriptional profile,
  following Agrawal et al. (2014) <doi:10.1016/j.cell.2014.09.050>,
  Extended Experimental Procedures, section 14.1. It is a
  re-implementation from that description rather than the original
  code, which was never released: it reproduces the published
  classification (99.7 percent against driver-mutation status under
  cross-validation) but not the exact published score values. The
  original study only scored the subset of samples for which both
  exome and RNA sequencing data were available, since exome data was
  needed to define the two reference groups (but not to score
  additional samples). This package lets you fit the reference
  centroids from any cohort with known BRAF-mutant/RAS-mutant labels
  and then score new samples — including ones with expression data
  only — using just their gene expression. Expression can be supplied
  as a matrix, a SummarizedExperiment or an ExpressionSet.

- [TiDEomics](/packages/TiDEomics) TiDEomics provides an integrated
  workflow for multi-group time-course omics data analysis. It uses
  model-based differential expression, variance decomposition, and
  co-expression module analysis with WGCNA to characterise temporal,
  group, and group-specific temporal effects. Multi-dataset and
  multi-omics comparisons and consensus module identification can be
  performed using shared features. The package includes quality
  control, data processing, functional enrichment, and integrated
  visualisation. It supports data with missing values, and operates
  on SummarizedExperiment objects for interoperability with the
  Bioconductor ecosystem.

- [TraianProt](/packages/TraianProt) A proteomics data analysis
  platform that enables the analysis of both label-free and labeled
  data from Data-Dependent or Data-Independent Acquisition mass
  spectrometry mode, supporting MaxQuant, MSFragger, DIA-NN,
  ProteoScape, and Proteome Discoverer output formats. TraianProt
  provides a comprehensive suite of stepwise downstream analysis
  modules, which includes data filtering, normalization procedures,
  and missing value imputation strategies. The platform also
  incorporates robust statistical frameworks for differential
  expression testing, with peptide-spectrum match level correction,
  thereby enhancing the reliability of biological interpretations.

- [TSENAT](/packages/TSENAT) Quantifies and models isoform-usage
  complexity in RNA-seq data using Tsallis entropy, a scale-dependent
  diversity measure. By tuning the entropic index parameter (q),
  TSENAT examines transcriptome heterogeneity at different scales:
  rare variants (low q) or dominant isoforms (high q). It enables
  computing Tsallis entropy and Tsallis divergence from
  transcript-level estimates, comparing measures between conditions,
  testing for differences, and visualizing scale-dependent complexity
  via q-curves.

- [TSSr](/packages/TSSr) TSSr package provides a comprehensive
  workflow on TSS data starts from identification of accurate TSS
  locations, clustering TSSs within small genomic regions
  corresponding to core promoters, and transcriptional activity
  quantifications, as well as specialized downstream analyses
  including core promoter shape, cluster annotation, gene
  differential expression, core promoter shift. TSSr can take
  multiple formats of files as input, such as Binary Sequence
  Alignment Map (BAM) files (single-ended or paired-ended), Browser
  Extension Data (bed) files, BigWig files, ctss files or tss tables.
  TSSr also generates various types of TSS or core promoter track
  files which can be visualized in the UCSC Genome Browser or
  Integrative Genomics Viewer (IGV). TSSr also exports downstream
  analyses result tables and plots.  Multiple cores are supported on
  Linux or Mac platforms.

- [tTEscanR](/packages/tTEscanR) Translation elongation is dependent
  on codon-anticodon interactions, with suitable nucleotide pairing
  being essential for efficient translation. To quantify this
  relationship, we previously developed a computational pipeline
  (GitHub - wgao688/sc_tRNA_mRNA) that uses mRNA codon usage relative
  to tRNA anticodon availability as a proxy for theoretical
  translation efficiency (tTE). Here, we introduce tTEscanR, a
  powerful and user-friendly R-based package that extends this
  approach to quantify translation efficiency from both bulk and
  single-cell sequencing data. tTEscanR is a versatile tool for
  exploring translation efficiency in diverse cellular processes,
  disease mechanisms, and therapeutic development. It also features
  an advanced visualization module to generate high-quality plots,
  enhancing result interpretation and communication.

- [WithinHostExperiment](/packages/WithinHostExperiment) Provides a
  SummarizedExperiment-based container and workflow functions for
  analysing intra-host single nucleotide variants (iSNV) from
  pathogen deep sequencing data. Supports replicate-aware quality
  control with full provenance tracking, multi-caller import
  harmonisation, within-host diversity quantification, and
  standardised export to transmission bottleneck estimation tools.
  Designed as infrastructure that bridges upstream variant callers
  (iVar, LoFreq, Freebayes) with downstream analysis packages through
  a common data object with uncertainty-preserving QC.

- [woven](/packages/woven) Supervised multi-omics integration for
  block-missing ("ragged") data, where some subjects lack one or more
  whole omics modalities (assay blocks) rather than scattered
  individual values. WOVEN fits a shared latent space with SUMCOR
  (sum of correlations) multi-set canonical correlation analysis, a
  criterion that maximizes the summed pairwise correlation between
  modalities' latent scores, using only the anchor subjects observed
  in every modality; class labels enter through a label-augmented
  cross-covariance, and per-modality graph Laplacians regularize the
  fit. Every other subject is then scored from the modalities it has
  by a Nystrom extension, a way to place new points in an already
  fitted embedding (for WOVEN's linear embedding, multiplying the
  observed features by the learned projection matrices), so no
  missing block is imputed. The effective sample size (ESS)
  retention, the fraction of subjects a method can score, is
  therefore the fraction with at least one observed modality, rather
  than the complete-case fraction. woven_tune() chooses the
  regularization settings by cross-validated balanced error rate, and
  the fitted model is summarized by plots and a table of the top
  features per modality.

- [wSIR](/packages/wSIR) Weighted Sliced Inverse Regression (wSIR) is
  a supervised dimension reduction algorithm for spatial
  transcriptomics gene expression data. For a provided gene
  expression matrix and dataframe of each cell's spatial coordinates,
  wSIR creates a low-dimensional representation of the gene
  expression data that preserves the ability to predict spatial
  coordinates that was present in the gene expression data.
  Furthermore, wSIR provides interpretable loadings which allow for
  projection of new single-cell gene expression data into a
  low-dimensional space which preserves the spatial information
  present in the gene expression data.

New Data Experiment Packages
=====================

There are 8 new data experiment packages in this release of Bioconductor.

- [AnnotatedBCGEData](/packages/AnnotatedBCGEData) AnnotatedBCGEData
  is an R Bioconductor package that provides access to more than 100
  gene-expression datasets generated using microarray and RNA-
  sequencing technologies. We collected these datsets from public
  repositories, including Gene Expression Omnibus, ArrayExpress, and
  The Cancer Genome Atlas. Where feasible, we obtained and
  reprocessed the raw data using computational pipelines and
  annotations that were standardized and more updated than those used
  to process the data originally. Additionally, we annotated the
  sample metadata using oncology terms from the National Cancer
  Institute Thesaurus. This was necessary because different data
  submitters used different terminology to describe breast-cancer
  patients and samples. Mapping the metadata to ontology terms makes
  it possible for users to search the data using these standardized
  terms.

- [CLAMPData](/packages/CLAMPData) Curated CLAMP resources (latent
  variables, loadings, priors, and example matrices) for use with the
  CLAMP package, distributed via ExperimentHub.

- [DaparToolshedData](/packages/DaparToolshedData) Mass-spectrometry
  based UPS proteomics data sets from Ramus C, Hovasse A, Marcellin
  M, Hesse AM, Mouton-Barbosa E, Bouyssie D, Vaca S, Carapito C,
  Chaoui K, Bruley C, Garin J, Cianferani S, Ferro M, Dorssaeler AV,
  Burlet-Schiltz O, Schaeffer C, Coute Y, Gonzalez de Peredo A.
  Spiked proteomic standard dataset for testing label-free
  quantitative software and statistical methods. Data Brief. 2015 Dec
  17;6:286-94 and Giai Gianetto, Q., Combes, F., Ramus, C., Bruley,
  C., Coute, Y., Burger, T. (2016). Calibration plot for proteomics:
  A graphical tool to visually check the assumptions underlying FDR
  control in quantitative experiments. Proteomics, 16(1), 29-32.

- [GSE142512](/packages/GSE142512) Provides ExperimentHub access to
  processed DNA methylation resources derived from GEO accession
  GSE142512. The resources contain Illumina HumanMethylation450 and
  Infinium MethylationEPIC data represented as
  RangedSummarizedExperiment objects with beta, M, and CN assays,
  hg19 genomic ranges, and sample phenotype metadata.

- [GSE280465](/packages/GSE280465) Provides ExperimentHub access to a
  processed DNA methylation resource derived from GEO accession
  GSE280465. The resource contains adult Illumina Infinium
  MethylationEPIC version 2 data from 163 samples across four tissue
  types, represented as a RangedSummarizedExperiment with beta, M,
  and CN assays, hg38 genomic ranges, and sample phenotype metadata.

- [HuMMANet](/packages/HuMMANet) Provides access to curated paired
  human microbiome-metabolome study bundles available through
  Bioconductor ExperimentHub. The package organizes matched sample
  metadata, species abundance profiles, original and harmonized
  metabolite profiles, and study-specific annotation resources
  including metabolite identifier mappings, microbial producer
  annotations, disease associations, pathway annotations, and drug
  similarity matches. The resource is designed to support consistent
  and reproducible cross-study microbiome-metabolome analyses.

- [PscanRBackgrounds](/packages/PscanRBackgrounds) Provides
  ExperimentHub metadata and documentation for the
  promoter-background score distributions used by PscanR
  transcription factor binding motif enrichment analyses. Each of the
  105 backgrounds is a separate ExperimentHub resource, covering the
  JASPAR 2020, 2022 and 2024 CORE collections, seven genome
  assemblies of five organisms, and five promoter windows. The files
  are hosted on Zenodo, one record per JASPAR release, and archived
  together in a version-specific Zenodo record.

- [SpaMTPData](/packages/SpaMTPData) Provides lightweight,
  version-aware access to spatial metabolomics and multi-omics
  experiment data used in SpaMTP workflows. Resources cover mouse
  brain, mouse urinary bladder, simulated single-cell multi-omics,
  annotation-refinement, and human brain examples. Large files are
  retrieved through ExperimentHub or verified source-URL caches
  rather than stored in the package repository. A small deterministic
  synthetic SpatialExperiment with a paired transcriptome supports
  offline examples and interoperability tests. Native
  SpatialExperiment datasets are available in a separate resource
  release; historical archival objects retain their original
  representations. Conversion and analysis are handled by the SpaMTP
  software package. Native experimental files share a Zenodo version
  family with SpaMTPdb without coupling the packages at runtime.


New Annotation Packages
=====================

There are 6 new annotation packages in this release of Bioconductor.

- [geneslator](/packages/geneslator) Geneslator is a comprehensive R
  package that performs gene identifier conversion and ortholog
  mapping. The tool integrates multiple cross-organism databases
  (NCBI, Ensembl, UniProt, GO, KEGG, Reactome, Wikipathways) and
  organism-specific resources within a single, coherent framework.
  Geneslator supports many organisms, including human, mouse, rat,
  yeast, worm, fly, zebrafish and arabidopsis.

- [JASPAR](/packages/JASPAR) Retrieves the JASPAR databases of
  transcription factor binding profiles through AnnotationHub. JASPAR
  is an open-access database of curated, non-redundant DNA binding
  profiles of transcription factors, stored as position frequency
  matrices. Each JASPAR release from 2022 onwards is an SQLite file
  registered in AnnotationHub. JASPAR() downloads the release you ask
  for once, keeps it in the AnnotationHub cache and returns its path,
  which can be searched with TFBSTools or queried with RSQLite.
  getAvailableJASPARVersions() lists the releases registered in the
  hub.

- [methylTFRAnnotationHg38](/packages/methylTFRAnnotationHg38)
  Precomputed transcription factor binding sites, motif GC frequency
  tables and a genome-wide GC distribution for use with the methylTFR
  package. Binding sites are genome-wide motif matches extended to
  the footprint window methylTFR reads across; the GC tables record
  how each motif's sites distribute over genome-wide GC quintiles,
  which is what allows a deviation score to be corrected for sequence
  composition. The data are hosted on AnnotationHub and retrieved on
  first use; this package provides the accessors and the
  documentation.

- [methylTFRAnnotationMm10](/packages/methylTFRAnnotationMm10)
  Precomputed transcription factor binding sites, motif GC frequency
  tables and a genome-wide GC distribution for use with the methylTFR
  package. Binding sites are genome-wide motif matches extended to
  the footprint window methylTFR reads across; the GC tables record
  how each motif's sites distribute over genome-wide GC quintiles,
  which is what allows a deviation score to be corrected for sequence
  composition. The data are hosted on AnnotationHub and retrieved on
  first use; this package provides the accessors and the
  documentation.

- [TxDb.Athaliana.Ensembl.TAIR10](/packages/TxDb.Athaliana.Ensembl.TAIR10)
  TxDb.Athaliana.Ensembl.TAIR10 provides transcript, exon, CDS, and
  gene annotations for Arabidopsis thaliana based on the TAIR10
  reference genome and Ensembl Plants release 63. The TxDb object is
  built from the official Ensembl Plants GFF3 annotation and
  distributed through AnnotationHub. Standard GenomicFeatures
  accessors such as genes(), transcripts(), exonsBy(), and cdsBy()
  can be used to query and extract gene models.

- [TxDb.Osativa.MSU.RGAP7](/packages/TxDb.Osativa.MSU.RGAP7) Provides
  transcript, exon, and CDS annotations for Oryza sativa (rice) as a
  Bioconductor TxDb object, generated from the Rice Genome Annotation
  Project (RGAP) Release 7 / MSU7 GFF3 file (genome assembly
  IRGSP-1.0). The TxDb object is distributed through AnnotationHub
  and can be retrieved directly with the AnnotationHub package.
  Standard GenomicFeatures accessors such as genes(), transcripts(),
  exonsBy(), and cdsBy() can be used to obtain gene models keyed by
  MSU locus identifiers. The package also provides msu_to_rap() for
  translating MSU locus identifiers into RAP-DB identifiers.

New Workflow Packages
=====================

There is 1 new workflow packages in this release of Bioconductor.

- [savingBiocObjects](/packages/savingBiocObjects) Documentation and
  guidance on strategies for saving and sharing Bioconductor objects,
  covering RDS serialization, HDF5-backed storage, the alabaster
  ecosystem for cross-language interoperability, exporting genomic
  ranges to standard formats, and preserving object metadata.

New Online Books
=====================

There are no new online books in this release of Bioconductor.

NEWS from existing Software Packages
===================================


[a4Base](/packages/a4Base)
------

                       Changes in version 1.61.1                        

- replace limma topTableF by topTable

[AHMassBank](/packages/AHMassBank)
----------

                        Changes in version 1.13                         

Changes in AHMassBank 1.13.2

- Add MassBank release 2026.11.

Changes in AHMassBank 1.13.1

- Add MassBank releases 2025.05.1, 2025.10 and 2026.03 with their
data
files being retrieved from GitHub.

[alabaster.base](/packages/alabaster.base)
--------------

                       Changes in version 1.14.0                        

- Bugfix to avoid requiring write access to files when reading
  VLS structures in various read*() functions.

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(), respectively. The former have
  been soft-deprecated and scrubbed from documentation for quite
  some time, so it is unlikely that anyone is using them.

- Automatically dispatch to alabaster.sfe and alabaster.treese
  for SpatialFeatureExperiment and TreeSummarizedExperiment
  objects, respectively.

[alabaster.bumpy](/packages/alabaster.bumpy)
---------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(). This includes related functions
  like useBumpyHDF5().

[alabaster.files](/packages/alabaster.files)
---------------

                       Changes in version 1.12.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(). This includes the related
  Wrapper classes.

[alabaster.mae](/packages/alabaster.mae)
-------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject().

[alabaster.matrix](/packages/alabaster.matrix)
----------------

                       Changes in version 1.14.0                        

- Add an explicit saveObject() method for HDF5Arrays, to avoid
  creating an unnecessary no-op file via the DelayedArray method.

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(). This includes related functions
  like writeSparseMatrix(), .createRawArraySeed(), etc.

- Support scrapper's LogNormalizedMatrixSeed in
  storeDelayedObject(). Also renamed the ResidualMatrix type hint
  to avoid potential confusion about the identity of the
  suggested class.

- Added a ToBeRealizedArray class for finer control of delayed
  operation preservation in saveObject() with
  Delayed.preserve.ops=TRUE.

[alabaster.ranges](/packages/alabaster.ranges)
----------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(), respectively. The former have
  been soft-deprecated and scrubbed from documentation for quite
  some time, so it is unlikely that anyone is using them.

[alabaster.sce](/packages/alabaster.sce)
-------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(), respectively. The former have
  been soft-deprecated and scrubbed from documentation for quite
  some time, so it is unlikely that anyone is using them.

[alabaster.se](/packages/alabaster.se)
------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(), respectively. The former have
  been soft-deprecated and scrubbed from documentation for quite
  some time, so it is unlikely that anyone is using them.

[alabaster.spatial](/packages/alabaster.spatial)
-----------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(), respectively. The former have
  been soft-deprecated and scrubbed from documentation for quite
  some time, so it is unlikely that anyone is using them. We also
  deprecated related functions like stageSpatialImage() and
  loadSpatialImage().

[alabaster.string](/packages/alabaster.string)
----------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(), respectively. The former have
  been soft-deprecated and scrubbed from documentation for quite
  some time, so it is unlikely that anyone is using them.

[alabaster.treese](/packages/alabaster.treese)
----------------

                 Changes in version 0.99.0 (2026-10-29)                 

- Initial version submitted to Bioconductor.

[alabaster.vcf](/packages/alabaster.vcf)
-------------

                       Changes in version 1.14.0                        

- Officially deprecate stageObject() and loadObject() in favor of
  saveObject() and readObject(), respectively. The former have
  been soft-deprecated and scrubbed from documentation for quite
  some time, so it is unlikely that anyone is using them. We also
  deprecate the separate load/stage methods for VCFHeader
  objects.

[amplican](/packages/amplican)
--------

                       Changes in version 1.35.6                        

- Handle missing forward-primer offsets when extracting reverse-only
  reads.

- Preserve valid empty results when assembling fully filtered
  alignments.

- Extend upstream truth tests and mutation checks for assignment,
  orientation,
  read multiplicities, batching, sampling, resume, and parallel
  execution.

- Clarify that primer matching uses the best overlap alignment.

- Count paired control evidence once per read/event during
  normalization.
  The control-frequency threshold remains strictly greater than
  min_freq.

- Correct relative cut-site overlap when a nonzero cut buffer is used.

- Recognize zero-start relative events in primer-overlap filtering.

- Resolve disagreements between two clipped reads without strand bias.

- Aggregate off-target features and safety-cap counts once per read.
  Cached filtered/normalized event tables and summaries require
  regeneration.

- Extend core decision tests and the developer mutation audit.

                       Changes in version 1.35.5                        

- Correct mismatch reference and replacement labels around insertions,
  including label reassignment between reads. Existing event tables
  require
  regeneration to receive this correction.

[ANCOMBC](/packages/ANCOMBC)
-------

                  Changes in version 2.15 (2026-09-17)                  

- Reduce the running time of ancombc, ancombc2, secom_linear, and
  secom_dist without changing the results

- Add the 'conservative' argument to ancombc2, which specifies whether
  the pseudo-counts of the sensitivity analysis are added before or
  after the estimation of sampling fractions

- Fix the error in ancombc2 when the trend test is requested with the
  sensitivity analysis but without the global test

[anndataR](/packages/anndataR)
--------

                        Changes in version 1.4.0                        

New features

- Add support for writing Zarr v3 stores, using (i) zarr_format
  argument of write_zarr()/as_ZarrAnnData() or (ii) the
  anndataR.zarr_format option. New stores are written as Zarr v3
  by default.

- Add support for reading backed objects using DelayedArray
  matrices, including a backed argument for HDF5AnnData,
  ZarrAnnData and InMemoryAnnData, and conversion of backed
  AnnData objects to SingleCellExperiment/Seurat.

- Add support for nullable-string-array elements in H5AD and
  Zarr.

- Add a group argument to read_h5ad(), write_h5ad(),
  as_HDF5AnnData() and HDF5AnnData, to read and write an AnnData
  stored in a group inside an HDF5 file, e.g. a modality inside a
  .h5mu file.

- Add read_elem() and write_elem() to read and write a single
  element of an H5AD file or Zarr store, e.g. read_elem(path,
  "obsm/X_pca"), without reading the whole object. With backed =
  TRUE, read_elem() returns arrays and matrices as DelayedArray
  objects.

Minor changes

- Support the anndata.AnnData class name used by Python anndata >= 0.13

- Handle the None layer key used by Python anndata >= 0.13 which
  stores X as a layer

- Read the sharded Zarr v3 stores that Python anndata >= 0.13
  writes by default. This requires Rarr 2.1.45 or later.

- Check that the indices of sparse matrices fit their shape when
  reading H5AD files and Zarr stores, with an error that suggests
  the indices may not be 0-based.

- Keep NA in double and integer matrices written to H5AD instead
  of replacing it with NaN, so that NA and NaN can still be told
  apart when reading the file back in R. Python anndata reads
  both as NaN.

- Use zstd compression by default in write_zarr()

- Error when writing a logical matrix containing NA to H5AD or
  Zarr, instead of writing a nullable-boolean that neither R nor
  Python anndata can read back. Convert it to numeric first.

- Require rhdf5 2.57.13 or later, which fixes a spurious
  NA_character_ warning for numeric datasets containing NA, and
  Rarr 2.1.59 or later, which fixes reading and writing
  zstd-compressed Zarr arrays on R 4.6.2.

Bug fixes

- Fix create_zarr() so that Zarr stores can be created at
  relative paths and at paths containing regex metacharacters.

- Use the format of an existing Zarr store when writing to it and
  truncate the store when mode = "w", so that a store can no
  longer end up with a mix of Zarr v2 and v3 nodes.

[See the full NEWS](/packages/release/bioc/news/anndataR/NEWS) or run `news(package = "anndataR")`.


[AnnotationDbi](/packages/AnnotationDbi)
-------------

                       Changes in version 1.76.0                        

NEW FEATURES

- The `organismKEGGFrame` helper function is deprecated for
  `genomeKEGGFrame` due to upstream API changes (@LiNk-NY)

[AnnotationHub](/packages/AnnotationHub)
-------------

                        Changes in version 4.3.0                        

NEW FEATURES

- (4.3.3) Add dispatchclass MspSpectra

[annotatr](/packages/annotatr)
--------

                       Changes in version 1.40.0                        

NEW FEATURES

- Add builtin_annotations_table(), which shows which groups of
  annotations
  are available for which genome builds, one row per genome.

- Update the vignette: a table of available annotations, the new
  annotations (MANE, Ensembl canonical, cCREs, any TxDb), current data
  sources and installation, and session information.

- Add ENCODE candidate cis-regulatory element (cCRE) annotations from
  the
  SCREEN registry (version 4) for hg38 and mm10:
  '[genome]_ccre_[class]' for
  the 8 classes (PLS, pELS, dELS, CA-H3K4me3, CA-CTCF, CA-TF, CA, TF),
  and the
  '[genome]_ccres' shortcut. The id is the cCRE accession.

- Add Ensembl canonical transcript annotations, '[genome]_canonical_*'
  and
  the '[genome]_basiccanonical' shortcut, for hg38, mm39, rn7,
  danRer11, dm6,
  and oviariramb2: one transcript per gene, of every biotype, including
  intergenic. They are built from the Ensembl 113 EnsDbs in
  AnnotationHub,
  with sequences renamed to UCSC names. Genes on alternate haplotypes
  and fix
  patches are left out.

- All annotations have entrez_id and ensembl_id columns, next to
  gene_id
  and symbol. For gene annotations of every genome, symbol, entrez_id,
  and
  ensembl_id (without version) are filled in from the org.*.eg.db
  package,
  the EnsDb, or the MANE summary, whatever kind of gene ID the source
  uses.
  When a gene ID maps to more than one, the first is used. gene_id is
  unchanged. summarize_genes() includes the new columns.

- Add build_txdb_annotations(), which builds the gene annotation types
  (promoters, exons, introns, intergenic, etc.) from any TxDb or EnsDb,
  e.g.
  one made from a GENCODE or RefSeq GTF with
  txdbmaker::makeTxDbFromGFF().
  Annotations are named [genome]_[group]_[type], e.g.
  'mm39_gencode_promoters', and work with summarize_genes() and the
  plotting
  functions. Gene symbols come from an EnsDb or an optional OrgDb.

[See the full NEWS](/packages/release/bioc/news/annotatr/NEWS) or run `news(package = "annotatr")`.


[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.26.0                        

USER VISIBLE CHANGES

- (v 1.25.3) Add api_reference_sha256 argument to Service() for
validating the API reference with a SHA-256 digest;
api_reference_md5sum validates with an MD5 digest again. For
backward compatibility, a 64-character api_reference_md5sum is
treated as a SHA-256 digest.
- Service() only retrieves the API reference at api_reference_url
when
api_reference_md5sum, api_reference_sha256, or api_reference_version
is provided, avoiding an unnecessary download when there is nothing
to validate. An unreachable api_reference_url is then no longer an
error.

[AnVILGCP](/packages/AnVILGCP)
--------

                        Changes in version 1.8.0                        

New features

- avworkflow_jobs() includes methodConfigurationNamespace and
methodConfigurationName columns identifying the workflow submitted;
the new columns follow the existing ones (#4, @ amstilp)
- avtable_import() supports list columns, e.g., tibble(id = 1,
cohorts
= list(c("A", "B"))); elements are uploaded as JSON arrays, which
AnVIL imports as list attributes (#5, @ smgogarten)
- New avtable_nest() reshapes list attributes returned by avtable()
(<name>.items and <name>.itemsType columns) or avtable_paged() into
list columns named as in the workspace DATA tab; avtable() output is
unchanged

Bug fixes and minor improvements

- Remove internal gcloud SDK helpers and the
print.gcloud_sdk_result()
method duplicated from GCPtools, avoiding the 'Registered S3 method
overwritten' message when both packages are loaded.
- Workspace write operations (avtable_import(), avtable_import_set(),
avtable_delete(), avtable_delete_values(), avdata_import(),
avworkflow_run(), avworkflow_stop(), and
avworkflow_configuration_set()) check the user's workspace access
level and report an informative error instead of e.g., a generic
HTTP 500 error when write (or compute) permission is lacking (#10).
- Use all_of() in avtable_import_set() to avoid a tidyselect
deprecation warning.
- avtable_import_status() warnings explain that Terra reports
re-importing values identical to the existing table as an "Error"
('All entities failed to update. There were 0 errors in total.'; #8,
@ smgogarten)
- avtable_import_status() no longer re-queries jobs with status
"Error"

[AnVILPublish](/packages/AnVILPublish)
------------

                       Changes in version 1.24.0                        

New Features

- (v. 1.23.1) Add dry.run argument to add_access(), as_notebook(),
and
as_workspace() to allow local processing and previewing without
making remote changes in AnVIL.

- (v. 1.23.3) Add AnVILBase to Suggests.

Significant User-Visible Changes

- (v. 1.23.1) Remove legacy support for python notedown conversion. A
system installation of Quarto is now required to convert .Rmd to
.ipynb notebooks, and a clear error is thrown if it is not found.

- (v. 1.23.1) Add quarto-cli to SystemRequirements in DESCRIPTION.

- (v. 1.23.1) Add URL and BugReports fields to DESCRIPTION.

- (v. 1.23.1) Add funding (fnd) role and grant number to package
authors.

Bug Fixes

- (v. 1.23.1) Handle lowercase .rmd extensions in as_notebook().

[AnVILVRS](/packages/AnVILVRS)
--------

                        Changes in version 1.0.0                        

- Initial release of AnVILVRS package for translating variant
identifiers and retrieving allele frequency data using GA4GH VRS
standards.

[apbCatalogR](/packages/apbCatalogR)
-----------

- Initial Bioconductor Submission

[arrayQualityMetrics](/packages/arrayQualityMetrics)
-------------------

                       Changes in version 3.69.2                        

MINOR CHANGES

- RColorBrewer dependency has been replaced by base grDevices
  package.

BUG FIXES

- cleanPhenoData function no longer errors when passed a
  data.frame. Thanks to GitHub users `@dah77` and `@w2felix` for
  reporting this issue.

[artMS](/packages/artMS)
-----

                 Changes in version 1.31.1 (2026-11-16)                 

- Fix the unit test that made R CMD check fail on Bioconductor: MSstats
  now
  returns an extra {LABEL} column in the {ModelQC} output of
  {artmsQuantification()}.
  The test now checks that the columns used by artMS are present
  instead of an
  exact number of columns.

- Enrichment analysis now uses the {gprofiler2} package instead of
  {gProfileR}.
  {gProfileR} is deprecated and its server no longer accepts requests,
  so the
  enrichment in {artmsEnrichProfiler()}, {artmsEnrichLog2fc()} and
  {artmsAnalysisQuantifications()} was failing. The results keep the
  same column
  names as before, with two differences:
  - the {domain} column now uses the g:Profiler source names ({GO:BP},
    {KEGG}, {REAC}, ...)
  - the gProfileR-only columns {subgraph.number} and {relative.depth}
    are no longer included
    {OMIM} is no longer available in g:Profiler and is ignored.
    Install {gprofiler2} to use these functions:
    {install.packages("gprofiler2")}

- Fix typo in the enrichment error message

[ASpli](/packages/ASpli)
-----

                       Changes in version 2.23.4                        

OTHER CHANGES

- gbCounts: removed the explicit calls to gc() inside the per-sample
  loops. They accounted for ~90% of the running time (gbCounts on the
  example data went from ~40 s to ~3 s).

                       Changes in version 2.23.3                        

OTHER CHANGES

- NAMESPACE: utils::data is no longer imported (conflicted with
  BiocGenerics::data and gave a warning when loading the package);
  limma is imported again, as required by its presence in Depends.

- Rd cross-references to edgeR now include the package anchor.

                       Changes in version 2.23.2                        

BUG FIXES

- Solve errors "subscript out of bounds" in gbCounts, jDUreport and the
  splicing reports with limma >= 3.99, where limma::strsplit2 now
  splits by fixed strings instead of regular expressions. ASpli uses an
  internal helper (.strsplit2) with the previous behaviour.

OTHER CHANGES

- Depends on R (>= 3.5.0), required by the serialized example objects
  in inst/extdata.

                       Changes in version 2.23.1                        

FUNCTIONS DEFUNCT

- The functions of the ASpli 1.x workflow, deprecated since the
  introduction of the new pipeline, were removed. See ?"ASpli-defunct"
  for a guide to migrate old scripts. Older Bioconductor releases still
  contain them.

- DUreport, DUreportBinSplice: use gbDUreport (DUreport.norm and
  DUreport.offset remain available).

- junctionDUreport: use jDUreport.

- mergeBinDUAS and filterDU: use splicingReport, integrateSignals and
  filterSignals.

- plotGenomicRegions: use exportSplicingReports and
  exportIntegratedSignals.

- readCounts, loadBAM: use gbCounts. They are no longer exported
  (internal helpers of gbCounts).

- AsDiscover: use jCounts. It is no longer exported (internal helper of
  jCounts).

- aspliBamsExample was removed together with loadBAM.

OTHER CHANGES

- NAMESPACE now uses explicit exports instead of exportPattern.

- Added examples and \value sections to man pages, and the man page
  ASpli-defunct.

- The example functions (aspliCountsExample, aspliASexample,
  aspliDUexample1/2, aspliJunctionDUexample) were broken because they
  read files missing in inst/extdata. They now read the objects stored
  in inst/extdata/aspli_example_objects.rds (result of the workflow on
  the example data). New: aspliSplicingReportExample and
  aspliIntegratedSignalsExample.

- Runnable examples for DUreport.norm, DUreport.offset, filterSignals,
  plotBins, rds, subset and the accessors, using the example objects.

[See the full NEWS](/packages/release/bioc/news/ASpli/NEWS) or run `news(package = "ASpli")`.


[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.37.3                        

- Fix the issues in BiocCheck. To decrease the package size, the
datas
in extdata is random subsampled.

                       Changes in version 1.37.2                        

- import cigarillo package.

                       Changes in version 1.37.1                        

- remove the imported package preseqR and adjust the authorship.

[ATACseqTFEA](/packages/ATACseqTFEA)
-----------

                       Changes in version 1.15.2                        

- update the documentation.

                       Changes in version 1.15.1                        

- import cigarillo package.

[augere.core](/packages/augere.core)
-----------

                        Changes in version 1.0.0                        

- First release of the augere.core package.

[augere.de](/packages/augere.de)
---------

                        Changes in version 1.0.0                        

- First release of the augere.de package.

[augere.gsea](/packages/augere.gsea)
-----------

                        Changes in version 1.0.0                        

- First release of the augere.gsea package.

[augere.screen](/packages/augere.screen)
-------------

                        Changes in version 1.0.0                        

- First release of the augere.screen package.

[BamScale](/packages/BamScale)
--------

- Initial Bioconductor Submission

[BANDITS](/packages/BANDITS)
-------

                       Changes in version 1.28.1                        

- 'DRIMSeq' dependency removed;

[bandle](/packages/bandle)
------

                       Changes in version 1.17.1                        

- minor bug fixes

[Banksy](/packages/Banksy)
------

                        Changes in version 1.9.5                        

- Point the large-dataset documentation link at the SeuratWrappers
upstream repository rather than a development fork

                        Changes in version 1.9.4                        

- The lazy path in runBanksyPCA is now reproducible by default. With
no seed, both pca_backend options start the iterative solver from a
fixed vector, so repeated runs give identical embeddings and the
caller's RNG stream is left untouched. Supplying a seed draws a
random starting vector instead, reproducibly, for callers who want
to vary it. Previously the C++ backend called set.seed(42)
internally, overriding the user's seed, and the R backend was
unseeded

                        Changes in version 1.9.2                        

- Lazy PCA mode in runBanksyPCA (lazy=TRUE, now the default) computes
PCA via an implicit linear operator without materializing the full
BANKSY matrix, scaling to millions of cells with low memory usage
- C++ backend (pca_backend="cpp", default) implements irlba with
OpenMP-parallelized sparse matrix-vector products and column-wise H0
scaling that computes neighborhood statistics without materializing
the gene-by-cell product matrix
- The lazy path supports multi-sample analysis with per-group kNN,
within-group scaling, optional parallel kNN via mclapply, and
multiple lambda values in a single call with kNN computed once and
reused.
- The lazy path reads on-disk matrices (e.g. BPCells) without
coercing
them to memory, allowing datasets beyond R's 2^31 sparse non-zero
limit
- OpenMP regions in the C++ backend are gated on a minimum work size,
so small datasets are not slowed by thread creation when the thread
count exceeds the available cores. Override with BANKSY_OMP_MIN_WORK
- runBanksyPCA with the lazy path has been benchmarked to 25 million
cells. Tiled Xenium Prime 5K, lambda=0.2, k_geom=30, npcs=50,
per-sample kNN and scaling, on 8 CPUs:

| Cells | Features | Sparsity | Samples | Runtime | Peak memory |
|---:|---:|---:|---:|---:|---:| | 1,221,372 | 5,101 | 96.18% | 3 |
1.6
min | 22.1 GB | | 2,035,620 | 5,101 | 96.18% | 5 | 2.5 min | 38.9 GB
|
| 5,292,612 | 5,101 | 96.18% | 13 | 6.6 min | 83.9 GB | | 10,178,100
| 5,101 | 96.18% | 25 | 12.5 min | 163.0 GB | | 25,241,688 | 5,101
| 96.18% | 62 | 34.6 min | 162.4 GB |

Runtime scales approximately linearly in cell count (1.3 * n^0.97
minutes, n in millions). Memory scales linearly up to 10M; the 25M
run
uses on-disk BPCells storage, which caps resident memory at the cost
of
runtime.

[barmixR](/packages/barmixR)
-------

- Initial Bioconductor submission.

[BatChef](/packages/BatChef)
-------

                        Changes in version 1.1.5                        

In the new Bioconductor 3.24 release, we introduce support for two
additional Python integration methods: Scanorama and BBKNN. We have
also
updated the AnnData object structure.

Additional improvements:

- Fixed an R CMD check NOTE related to undefined global variables.
- Updated the prediction output.

                        Changes in version 1.0.2                        

- Update SVM for prediction

                        Changes in version 1.0.1                        

- Fix bugs and add visualization of prediction on vignette

[BatchQC](/packages/BatchQC)
-------

                        Changes in version 2.9.3                        

Minor Changes

- Fixed a Unicode character in run_kBET documentation that broke
checks

                        Changes in version 2.9.2                        

Minor Changes

- Updated some descriptions to support manuscript and provide info.
- Added additional details to method functionality in the vignette.

                        Changes in version 2.9.1                        

Minor Changes

- Updated discrepancy in roxygen comment for lambda stat
- Added headers for all vignette code chunks
- Corrected spacing and line indent issues for a number of man files

[BatchVaria](/packages/BatchVaria)
----------

- Initial Bioconductor submission

[beachmat](/packages/beachmat)
--------

                       Changes in version 2.30.0                        

- Added a num.groups= argument to tatami.sums.by.group(), so that
  we don't have to guess the number of groups.

- Improved efficiency of tatami.multiply() with new C++
  algorithms. This comes at the cost of not correctly handling
  non-finite values in multiplications involving sparse matrices.

- Added a tatami.variances() function to easily compute means and
  variances from an initialized pointer.

[BEclear](/packages/BEclear)
-------

                 Changes in version 2.29.2 (2026-09-14)                 

- Bug fixes
  - Batch effects are detected in both directions again, as 1.12.x did.
  - correctBatchEffect no longer fails with when glue isn't installed

- Tests
  - Test for batches shifted downward in calcSummary
  - Test for logging when logger formats with sprintf

                 Changes in version 2.29.1 (2026-06-07)                 

- Tests
  - Reintroduced integration tests for correctBatchEffect

[bedbaser](/packages/bedbaser)
--------

                        Changes in version 1.6.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Update bedbase to 0.15.0 (#32, @jwokaty)

NEW FEATURES

- Add a host argument to BEDbase() for using another API server, such
as a development or test server (#5, @nsheff)
- Add a show() method for BEDbase objects (#31, @lwaldron)
- Add after to bb_list_beds() for paging by the id of the last record
(#34, @jwokaty)

BUG FIXES AND MINOR IMPROVEMENTS

- Fix bb_list_beds() failing on later pages, and report bedhost
errors
clearly instead of failing or requiring xml2 (#34, @jwokaty)

[beer](/packages/beer)
----

                 Changes in version 1.15.1 (2025-11-04)                 

- Updated NEWS versioning to reflect Bioconductor updates.

- Updated mislabeled progressr::progress dependency.

[BenchHub](/packages/BenchHub)
--------

                        Changes in version 1.0.0                        

- Initial release.

[betterChromVAR](/packages/betterChromVAR)
--------------

                       Changes in version 1.99.1                        

New features

- The background expectation of each motif is now rescaled so that
its
global expectation equals that of the motif's peaks (neglible effect
in single-cell data, better agreement with chromVAR in bulk data).
- Motif interaction analysis, including functions moved from
weightedMotifAccess (and improved), as well as peak-level functions
and a vignette.
- fitDeviations() (moved from weightedMotifAccess)
- peakLevelMotifTests() and stepwisePeakLevelMotifTests() as
peak-logFC-based alternative to the analysis of deviations.
- Motif interaction analysis, including functions moved from
weightedMotifAccess (and improved), as well as peak-level functions
and a vignette.

                       Changes in version 1.1.14                        

New features

- ridge argument of betterChromVAR() and computeDeviationsAnalytic():
optional decorrelation of the motifs that share peaks (ridge
deconvolution based on the overlap between motifs after background
subtraction).
- getGradedMatches(): graded (instead of binary) motif matches,
giving
partial weights to weak matches according to their score.

Changes

- Removed the bias shrinkage options -- see downstream smoothing
instead.
- The bin-to-bin background probabilities of getBackgroundBins() are
now balanced (tol argument; tol=NULL restores the earlier behavior).
This slightly changes the results, removing a bias in mean z.
- Minor gains in efficiency.

                       Changes in version 1.1.13                        

- added pseudoCount and loo to make CVnorm more robust.
- added adaptive kernel width

                       Changes in version 1.1.12                        

New features

- smoothDeviationsKNN(): adaptive EB smoothing of single-cell
deviations across nearest-neighbor cells.

Bug fixes

- Fixed wrong denominator in computeDeviationsWeighted
- Fixed mistake in background variance (affecting the scale of the
z-scores): it is now the variance of the counts of each bin's
background peaks (as drawn by sampleBackgroundPeaks()), computed
before any shrinkage

Changes

- The z-scores (and things based on it, e.g. motif variability) now
include a finite-population correction (fpc=TRUE) of the null
variance.
- functions now also return an effective N per motif that is used in
normalizeDevsForSize for more accurate adjustment.

[See the full NEWS](/packages/release/bioc/news/betterChromVAR/NEWS) or run `news(package = "betterChromVAR")`.


[BgeeDB](/packages/BgeeDB)
------

                       Changes in version 2.39.1                        

- topAnat can use the expression calls of Bgee 16. These calls
  are made in anatomical structures possibly restricted to a cell
  type

- New argument callFileType of loadTopAnatData: the expression
  calls propagated by Bgee, default from Bgee 16, or the observed
  ones, that topGO propagates

- From Bgee 16, the classic Fisher test, without decorrelation,
  runs on the propagated calls: BgeeDB computes it itself,
  without topGO, that would propagate them again, with the same
  counts and p-values as the classic algorithm of topGO, in
  seconds rather than minutes for the species with many
  conditions. termStat, genesInTerm, sigGenes, makeTable and
  printTopAnatGraph work on its object. The decorrelation
  algorithms run on the observed calls. Any other combination
  throws an error

- From Bgee 16, topAnat_galaxy runs the classic algorithm on the
  propagated calls

- From Bgee 16, new argument background of topAnat: the gene list
  is restricted by default to the protein coding genes,
  downloaded from Bgee; "polyA", "observed" or a vector of gene
  IDs choose another background. The previous releases are not
  concerned

- New argument foregroundGenes of topAnat, to provide the genes
  of interest alone instead of a gene list: the gene list is
  built from the background

- From Bgee 16, a Bgee object created without any data type uses
  the data types of the release

- New function printTopAnatGraph, drawing in a PDF file the
  conditions significant in a topAnat analysis and their
  ancestors, with their names

- makeTable and printTopAnatGraph leave out by default the
  conditions redundant with a more precise one: those with the
  numbers of annotated and significant genes and the p-value of
  one of their sub-conditions. New argument keepRedundant to keep
  them, makeTable then telling them in a column redundant

- loadTopAnatData reports the answer of the server when it
  refuses a query

- loadTopAnatData reads the expression calls faster and with less
  memory: with fread, and by checking that a download is complete
  from its last bytes rather than by reading it whole. From Bgee
  16, the Bgee API sends the expression calls compressed, and
  they are kept compressed in the download directory, in a file
  ending with .tsv.gz

[bHIVE](/packages/bHIVE)
-----

- Initial Bioconductor Submission

[BindingSiteFinder](/packages/BindingSiteFinder)
-----------------

                       Changes in version 2.11.2                        

Bug fixes

- Condition names containing '_' caused errors in filterBsBackground()
  and
  reproducibilityFilter(). Conditions are now matched by sample id
  instead of being parsed from sample names.

- reproducibilityFilter() applied wrong cutoff values due to rounding
  errors

- Fixed duplicated transcriptRegion meta data column in
  assignToTranscriptRegions()

- Fixed missing colors in targetGeneSpectrumPlot() for more than 7 gene
  types

- Fixed typo in the title of duplicatedSitesPlot() (#20)

Changes

- Condition names are now validated when creating a BSFDataSet or
  calling
  setMeta(). Allowed are letters, numbers and one type of separator
  ('-', '_'
  or '.') that has to be the same for all conditions. Names with other
  characters (e.g. spaces) or mixed separators now give an error

- Replaced deprecated ggplot2 'size' with 'linewidth' for lines in all
  plots

Documentation

- Clarified in the vignette that overlaps = "flag" is only available in
  assignToTranscriptRegions(), not in assignToGenes()

- Use txdbmaker::makeTxDbFromGFF() explicitly in the example data
  script

- Removed links to the non-existent help topic 'overlaps'

                       Changes in version 2.11.1                        

- Fixed a bug in assignToGenes that caused duplicated meta data columns

- Updated documentation to roxygen 8 and re-documented

[Bioc.gff](/packages/Bioc.gff)
--------

                        Changes in version 1.4.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The ##species metadata added by import() now has a Rank element
next
to Taxonomy ID and Organism. A ##species URL that gives an organism
name is treated as a search term, and every matching taxon is
returned.

BUG FIXES AND MINOR IMPROVEMENTS

- Organism names are found again for the ##species NCBI Taxonomy
Browser URL. The Taxonomy Browser now sits behind a bot check, so
the lookup uses the NCBI E-utilities API instead of scraping the
page.
- export() to GFF and GTF files no longer fails with could not find
function "score" when BiocGenerics is not attached; score and strand
are now imported (#1, @jdrnevich).
- import() no longer fails with "replacement 'metadata' value must be
a list" when the NCBI lookup for ##species fails. It now warns,
gives the reason, and returns empty metadata.
- Replace rvest with jsonlite in Suggests.
- Use roxygen2 8.1.1 and regenerate NAMESPACE.
- Add descriptions to @describeIn tags and use @aliases tags in the
GFFFile-class documentation.
- Add tests for NCBI name searches and Rank, GTF export without
BiocGenerics attached, and GFF import when the NCBI lookup fails.

[BiocBook](/packages/BiocBook)
--------

                       Changes in version 1.11.7                        

New features

- enrich_llms_txt() names the Markdown copy of each page after the
page, <page>.md rather than quarto's <page>.llms.md, so that it
lives at the page's URL with .md in place of .html (e.g.
pages/Chapter-2.md for pages/Chapter-2.html). The links of llms.txt,
and those between the copies, follow.
- Every page of a book links to its Markdown copy: in its <head>
(<link rel="alternate" type="text/markdown">), for AI assistants,
and under its title ("View as Markdown"), for readers.
?BiocBook-llms says how to hide the latter.
- New books are created from template 1.1.3, whose biocbook GitHub
workflow writes a .nojekyll file at the root of the book's GitHub
Pages site. Otherwise, GitHub Pages runs Jekyll, which renders every
.md file into an HTML page of the same name, colliding with the
book's own pages. Books created from an earlier template need the
same change: in .github/workflows/biocbook.yml, add touch
llms-index/.nojekyll right after mkdir -p llms-index.

                       Changes in version 1.11.6                        

Bug fixes

- New books are created from template 1.1.2, whose biocbook GitHub
workflow no longer uses catchpoint/workflow-telemetry-action. Its
repository no longer exists, so every run failed before building
anything. Books created from an earlier template need the same
change: remove the "Collect Workflow Telemetry" step from
.github/workflows/biocbook.yml.

                       Changes in version 1.11.5                        

New features

- setup_python() installs the pip: section of inst/requirements.yml,
with the environment's own pip, once its conda packages are in
place. Wheels from PyPI carry their compiled dependencies under
names of their own, so the python packages a book imports can no
longer pick up a library the R session has already loaded. On the
Bioconductor builders, conda builds of h5py (and so cooler) failed
to import with OPENSSL_3.2.0 not found, and matplotlib >= 3.11 with
undefined symbol: hb_ft_font_get_ft_face. ?BiocBook-python explains
which packages to declare under pip:, and which to keep under conda,
such as command-line tools.
- setup_python() puts the command-line tools of the book's
environment
(e.g. samtools) on the PATH from the first chunk on, rather than
once python has started.
- llms.txt lists the packages of the pip: section after the conda
ones.
- New books are created from template 1.1.1, whose
inst/requirements.yml declares the python packages a book imports
under pip:, and keeps conda for python itself and for command-line
tools, rather than pinning matplotlib below 3.11.

[See the full NEWS](/packages/release/bioc/news/BiocBook/NEWS) or run `news(package = "BiocBook")`.


[BiocCheck](/packages/BiocCheck)
---------

                       Changes in version 1.50.0                        

NEW FEATURES

- Identify whether package source directories are live git clones or
  untarred build directories, and adjust check behavior accordingly

- Update vignette engine checks to recommend `VignetteBuilder` instead
  of
  `VignetteEngine` in the DESCRIPTION file (@LiNk-NY, #246).

- Add new check `checkValidDESCfields` to validate and flag unknown or
  non-standard fields in the `DESCRIPTION` file.

- Add new check `checkBiocType` to validate the `BiocType` field in the
  `DESCRIPTION` file against the allowed values: `Software`,
  `ExperimentData`, `Annotation`, `Workflow`, `Book`, and `Docker`.
  Invalid values no longer override the package type inferred from
  `biocViews`, e.g., when checking the package size. Common
  misspellings,
  such as `AnnotationData`, and incorrect capitalization prompt a
  suggestion of the valid value.

- Add check for S4 classes to verify they provide a non-derived default
  `show()` method.

- Write a machine-readable `00BiocCheck.json` report next to
  `00BiocCheck.log` in the `<PackageName>.BiocCheck` folder. It
  includes the
  session metadata, a `summary` count of the errors, warnings, and
  notes, the
  overall `status`, one structured `entries` record per condition
  raised
  (with the originating check function and file locations where
  reported),
  and the plain text report. See the vignette for the schema.

- Extend `checkRemotesUsage` to flag `basilisk` packages that disable
  version checks with `setBasiliskCheckVersions()` or that install
  Python
  packages from VCS / URL specifications, e.g., `git+https://...`
  (@Artur-man, #263).

BUG FIXES AND MINOR IMPROVEMENTS

- `checkS4ClassShowMethod` no longer requires a `show()` method for S4
  classes that extend an S3 or basic class, e.g., `contains = "array"`,
  since the default `show()` method prints the underlying data
  (@hpages, #268).

- `check1toN` no longer flags sequences with a numeric literal upper
  bound, e.g., `1:10`; only `1:n`-style sequences, which can be
  unexpectedly decreasing when `n` is 0, are reported (@Bisaloo, #256).

- Remove the Support Site registration and 'Watched Tags' checks as
  Bioconductor transitions away from the Support Site. The
  `no-check-bioc-help` option now only disables the Bioc-devel mailing
  list
  subscription check (@lshep, #267).

[See the full NEWS](/packages/release/bioc/news/BiocCheck/NEWS) or run `news(package = "BiocCheck")`.


[BiocDuckDB](/packages/BiocDuckDB)
----------

- Initial Bioconductor Submission

[BiocGenerics](/packages/BiocGenerics)
------------

                       Changes in version 0.60.0                        

NEW FEATURES

- Add the following S4 generics:
  - data() -- commit 93a0f55
  - detail() -- commit 0c73da8
  - group() -- commit 9f6939f
  - rotate() -- commit fab8ae0
  - scale() -- commit b221557
  - sequence() -- commit 54579af
  - transform() -- commit 54579af
  - version() -- commit 9263b0b

SIGNIFICANT USER-VISIBLE CHANGES

- Change interface of table() generic to table(x, ...).
  See commit 235a41a

BUG FIXES

- paste2() now propagates NAs, as it should.

[BiocNeighbors](/packages/BiocNeighbors)
-------------

                        Changes in version 2.8.0                        

- Added findResultsToMatrix(), findResultsToHits(),
  queryResultsToMatrix() and queryResultsToHits(), which convert
  the results of findKNN() to more amenable formats.

[BioCor](/packages/BioCor)
------

                        Changes in version 1.38                         

- Update according to latest BiocChecks recommendations.
- Fix a bug on combineScores: if a vector was all NA the calculation
was performed.
- Removed all sapply() calls.
- Used shortcut logical operators (&&, &#124;&#124;) whenever possible.

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.32.0                        

DEPRECATED AND DEFUNCT

- Removed defunct `pkgType = "all"` argument in `biocDownloadStats` and
  defunct `pkg` argument in `biocRevDepEmail`

SIGNIFICANT USER-VISIBLE CHANGES

- Updated R dependency to >= 4.6.0

BUG FIXES AND MINOR IMPROVEMENTS

- Cached Bioconductor files (VIEWS, download stats, etc.) are no longer
  treated as stale immediately after download. Staleness is now judged
  from the remote size, ETag and Last-Modified instead of
  `BiocFileCache::bfcneedsupdate()`, whose day-granularity `expires`
  check marked every file served with a short `max-age` as expired and
  forced a re-download on every call (#92).

- Tests that need the network are now all gated, and setting
  `BIOCPKGTOOLS_NO_NETWORK=true` skips them, so the suite runs offline
  in seconds (#90).

- `pkgBiocRevDeps(recursive = TRUE)` no longer returns forward
  (non-reverse) dependencies. Recursion now follows only strong
  dependencies, avoiding the `Suggests` / `Enhances` explosion that
  pulled
  in nearly the entire repository (#81).

- Switch `generateBiocPkgDOI` to use `httr2` instead of `httr`

- Switch `repositoryStats` to use `httr2` instead of `httr`

- Support `match.arg` for `event` argument in `generateBiocPkgDOI`

- Removed the dependency on `httr` by migrating `.url_exists` helper to
  `httr2`

- Removed duplicate internal function `generateBiocPkgDOI_REST`

- Raised unit test coverage from 18% across the package (#85). Fixed
  `.get_all_biocpkgs()` calling `BiocManager::version` unevaluated
  (always errored) and `checkDeps()` dropping to an atomic vector on
  subset, which broke `problemPage(dependsOn = ...)` entirely.

[biocthis](/packages/biocthis)
--------

                       Changes in version 1.23.1                        

BUG FIXES

- Fixed a bug about the regex related to r-universe. Thanks to
@jan-spr https://github.com/lcolladotor/biocthis/pull/61. Reviewed
by @nick-eagles.

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.70.0                        

NEW FEATURES

- More functions have gained a `useCache` argument, which defaults to
  `TRUE`,
  to reduce the load on the upstream server.

BREAKING CHANGES

- The deprecated Ensembl mirrors ("asia" and "useast") have been
  removed.

- This package now defaults to https://jun2026.archive.ensembl.org
  where it
  was previously defaulting to https://www.ensembl.org/. The June 2026
  (Ensembl
  116) release is the last release to include the Biomart utility.

- Modification of SSL settings via `setEnsemblSSL()` is now deprecated
  and
  will be removed in the next biomaRt release. Please get in touch if
  your
  workflow depends on this feature.

MINOR IMPROVEMENTS

- Unnecessary call to `listEnsemblArchives()` in `getBM(useCache =
  FALSE)`
  has been removed.

- A complete static dataset of previous Ensembl releases is now shipped
  with
  this package. This solves an issue where the upstream Ensembl website
  stopped
  listing all releases on their page. The dataset is directly
  accessible to
  users via `data("ensembl_releases")`. It is also possible to get the
  list of
  Ensembl archives without an internet connection via
  `listEnsemblArchives(fetch = FALSE)`.

INTERNAL CHANGES

- Vignettes now run again on continuous integration to ensure they stay
  up
  to date.

- Unused internal functions have been removed.

- Code duplication across functions have been reduced.

[biotmle](/packages/biotmle)
-------

                       Changes in version 1.37.1                        

- Removal of the `superheat` dependency, which has been archived on
  CRAN. The
  heatmap drawn by `heatmap_ic()` is now built directly with `ggplot2`.

- `heatmap_ic()` gains a `row_dendrogram` argument, which restores the
  row
  dendrogram previously available through `superheat`. Biomarkers are
  ordered
  by hierarchical clustering whether or not the dendrogram is drawn, so
  this
  argument affects only the display. It defaults to `FALSE`, matching
  the
  default of the `superheat` argument it replaces. Drawing the
  dendrogram
  introduces dependencies on `ggdendro` and `patchwork`.

- `heatmap_ic()` now returns a `ggplot` object, which can be stored and
  further
  modified in place, rather than drawing to the active graphics device
  as a
  side effect.

- `heatmap_ic()` gains a `scale` argument controlling whether each
  biomarker's
  contributions are standardized across observations before plotting;
  this
  replaces the argument of the same name previously forwarded to
  `superheat`.

- Arguments passed to `heatmap_ic()` via `...` are now forwarded to
  `ggplot2::geom_tile()` rather than to `superheat::superheat()`. Calls
  that
  relied on `superheat`-specific arguments (e.g., `left.label`,
  `clustering.method`, `row.dendrogram`) must be updated.

[BreastSubtypeR](/packages/BreastSubtypeR)
--------------

                       Changes in version 1.5.14                        

Performance

- sspPredict() avoids repeated ranking by computing Spearman
correlations in one matrix call per eligible panel. Panels with
non-finite values, constant columns or unsupported classes retain
the existing scalar calculation and its warning and error behavior.

                       Changes in version 1.5.13                        

API compatibility

- BS_ssBC() uses PatientID as a unique sample lookup key, not a
repeated biological-patient identifier. For ER, ER.v2, TN and TN.v2,
exactly one original column and one exact converted PatientID column
must provide the same complete, nonempty, unique labels as the assay
names in order. This intentionally tightens compatibility even for
excluded/empty rows and contradictions whose numerical results
agree. The converted-key check prevents metadata flattening or
partial matching from changing the key actually used. Supported
storage must satisfy both boundaries; the existing empty four-class
restrictions still apply. Guard refusals use
bsr_patientid_mismatch_error. Manual BS_Multi() propagates that
error when dispatch reaches ssBC; earlier preconditions such as a
missing PatientID column can still give an unclassed error. AUTO
routing and its early order check are unchanged, but selected
ssBC/ssBC.v2 wrappers can newly refuse duplicate original columns or
ambiguous converted keys that pass that precheck. Existing TSIZE
metadata preparation errors in .nc_pheno() occur first; this is not
a priority promise for missing clinical columns, invalid s, NODE
errors or missing receptors. No identifiers are realigned or
overwritten. Reconcile source records; internally consistent keys
cannot detect an upstream sample mix-up. Other standalone NC
wrappers are unchanged.

- Four-class (Subtype = TRUE) ssBC with s = "ER", "ER.v2" or "TN" no
longer stops with "undefined columns selected" when no sample
belongs to any selected receptor subgroup and the input is otherwise
well formed. Standalone BS_ssBC(), manual BS_Multi() and the AUTO
mixed-cohort route with no recognized ER/HER2 pair return the
existing five-class empty result plus missing BS.Subtype,
Call.Subtype and predictions.FourSubtype. This needs at least one
represented PAM50 gene, unique plain-vector receptor columns and a
single character or factor PatientID column whose labels are
identical to the assay column names in assay order. Other inputs
remain errors, with identity diagnostics as described above.
Receptor labels are matched literally and are not normalized;
alias-coded AUTO inputs still fail, now with the same error as the
corresponding five-class call. TN.v2, all five-class results and the
other methods' results, summaries and entropy are unchanged.

[See the full NEWS](/packages/release/bioc/news/BreastSubtypeR/NEWS) or run `news(package = "BreastSubtypeR")`.


[BrowserViz](/packages/BrowserViz)
----------

                       Changes in version 2.35.1                        

- Fixed unmatched braces in send() @param msg documentation (Rd NOTE).

- Added GitHub Actions CI workflow (R-CMD-check-bioc).

- Added automated push to Bioconductor devel on merge to master.

                       Changes in version 2.35.0                        

- Version bump due to Bioconductor 3.23 devel synchronization.

[BSgenomeForge](/packages/BSgenomeForge)
-------------

                       Changes in version 1.14.0                        

- No significant changes in this version.

[bugsigdbr](/packages/bugsigdbr)
---------

                       Changes in version 1.20.0                        

- importBugSigDB() defaults to BugSigDB release v1.4.0
(10.5281/zenodo.24379934): 11,939 signatures from 1,847 studies, up
from 7,425 signatures from 1,246 studies in v1.3.1.
- New Study columns: Curation method, Curation protocol, submittedBy,
submissionDate, Review state and Quality control tag.
- Releases exclude retracted studies, studies whose review found them
incorrect, studies with neither a PMID nor a DOI, and signatures
with no taxa resolved to an NCBI Taxonomy ID. version = "devel"
still includes them; retracted and incorrect studies are marked in
Quality control tag and Review state.
- Releases now include experiments without signatures from reviewed
studies, which record a negative taxonomic result and alpha
diversity. Their rows in importBugSigDB() have NA signature columns;
getSignatures() skips them.
- DOIs are normalized to the 10.xxxx/... form. Studies whose DOI was
written as a doi.org URL were previously left out of releases; they
are now included.
- Release v1.4.0 uses NCBI taxonomy from ncbi-taxon-db 2026.11.1
(NCBI
taxdump of 1 November 2026).

[CAGEr](/packages/CAGEr)
-----

                       Changes in version 2.19.2                        

BUG FIXES

- Fix plotReverseCumulatives(group = "sampleLabels"), which failed
with 'by' must specify a uniquely valid column since
S4Vectors 0.51.8 stopped de-duplicating column names in
as.data.frame() of a DataFrame: the sampleLabels column was selected
twice.
- Remove the utils::data import, which clashed with the new data()
generic of BiocGenerics and caused a replacing previous import
warning at load time. importPublicData() now calls utils::data()
explicitly.

                       Changes in version 2.19.1                        

BACKWARDS-INCOMPATIBLE CHANGES

- In aggregateTagClusters(), tpmThreshold now only selects the tag
clusters that define the boundaries of the consensus clusters. Two
erroneously applied filters are removed: the one that kept only the
CTSSs with a normalised signal above tpmThreshold in at least one
sample, and the one that dropped the CTSSs whose normalised signal
summed over all samples was not above tpmThreshold.
- The excludeSignalBelowThreshold argument of aggregateTagClusters()
is deprecated. It was erroneously discarding the individual CTSSs
below tpmThreshold in every sample.
- CustomConsensusClusters() no longer drops the CTSSs whose
normalised
signal summed over all samples is 1 or less. Its threshold and
nrPassThreshold arguments still apply.

BUG FIXES

- Fix vignette build with recent S4Vectors, whose table() method for
Rle objects does not accept the useNA argument anymore.

OTHER CHANGES

- Accelerated paraclu() ~10×, and distclu() and
aggregateTagClusters()
2–4×, with unchanged results. Clusters that do not overlap any CTSS
now stop the computation with an explicit error message.
- Accelerated exportToTrack() 10–30× when qLow and qUp are given,
with
unchanged results.

[CalibraCurve](/packages/CalibraCurve)
------------

                        Changes in version 1.3.1                        

- add warning if minReplicates is set to 1
- allows turning off log10 transformation of x-axis for plots
- allow usage of column names (as alternative to column numbers) when
reading data
- option "multiplot" and "all_in_one" now also available for response
factor plot
- add legend for response factor plot
- fix bug when multiple input files are used and at least one result
cannot be computed
- add option for "math" appearance of x- and y-axis labels

[canceR](/packages/canceR)
------

                       Changes in version 1.47.2                        

- add param for dialogSpecificMut

                       Changes in version 1.47.1                        

- update DESCRIPTION file and add pandoc dependency

- remove @usage tag of dialogSpecificMut as cause Warning

[CanvasXpressBio](/packages/CanvasXpressBio)
---------------

- Initial Bioconductor Submission

[Cardinal](/packages/Cardinal)
--------

                 Changes in version 3.15.1 (2026-11-18)                 

New features

- Added 'findIsotopesByCorr()' for grouping isotope peaks by m/z
spacing and image correlation, with a target-decoy score cutoff

- Added 'monoisotopes()' to filter findIsotopesByCorr() results to
one
feature per isotope group

[carnation](/packages/carnation)
---------

                        Changes in version 1.1.5                        

Updates for Bioconductor 3.24 release include the following:

- a new volcano plot module

- various plot updates:

- aspect ratio controls for MA, PCA and scatter plots.

- MA plot has additional alpha control and 'Refresh' button for
reduced reacitivity.

- Updated landing page

- Analysis table shows all available datasets
- More info (file size, mtime), improved selection behavior

                        Changes in version 1.1.2                        

- New volcano plot added

- various plot updates

- MA plot
- added aspect ratio and alpha controls
- new 'refresh' button for reduced reactivity
- Gene plot
- add option to show unnormalized counts. This is helpful when
the data in the dds element is pre-normalized, e.g. in
proteomics
- scatter plot
- raw p-values added as a plotting option
- plot/table selection panels now open upon selection
- aspect ratio control for PCA & scatter plots

- analysis table

- updated to show all available datasets
- obtain file info (size, mtime) and show in the table
- improved selection behavior

- misc fixes

- correctly reset gene scratchpad when loading data so previously
added genes don't show up again.
- the gene plot function now sets the 'text' aesthetic to factor()
to stop warnings about statistical transformation
- fromList.with.names now uses base R column indexing instead of
dplyr::relocate.

[CatsCradle](/packages/CatsCradle)
----------

                        Changes in version 1.7.1                        

- Bugfix to allow use of user-provided LRN.
- Remove check that species is human or mouse if LRN is provided by
user.
- Heirachical ordering of LRInteraction heatmap.

[causalASCEND](/packages/causalASCEND)
------------
- Initial Bioconductor Submission

[cBioPortalData](/packages/cBioPortalData)
--------------

                       Changes in version 2.26.0                        

New features

- Add structuralVariantData to download structural variant (e.g.,
gene
fusion) data. The endpoint is listed only in the 'internal' group of
the cBioPortal OpenAPI v3 specification and is therefore called
directly.
- Add genericAssayData and genericAssayMeta to download generic assay
data (e.g., arm-level copy number, PanCancer Atlas methylation,
phosphoprotein) and the corresponding entity identifiers (stableId)
- fetchData routes 'STRUCTURAL_VARIANT' and 'GENERIC_ASSAY' molecular
profiles to the new functions based on molecularAlterationType
- Add genericAssayStableIds argument to fetchData

Bug fixes and minor improvements

- fetchData uses all molecular profiles in the study (excluding
'GENERIC_ASSAY') when molecularProfileIds is not specified (NULL or
NA). Previously, no profiles were downloaded when using the default
molecularProfileIds in cBioPortalData and getDataByGenes.
- getDataByGenes returns data without an entrezGeneId column (i.e.,
structural variant and generic assay data) as is rather than
attempting to join gene metadata
- cBioPortalData excludes 'STRUCTURAL_VARIANT' and 'GENERIC_ASSAY'
profiles, which are not indexed by gene, and warns when they are
explicitly requested
- Use the SHA-256 digest of the API specification in cBioPortal() to
avoid the 'service version differs from validated version' warning;
AnVIL now validates the specification with SHA-256 rather than MD5
- Remove the version constraint on AnVIL in DESCRIPTION
- Add unit tests for structural variant, generic assay, and fetchData
profile routing

[celda](/packages/celda)
-----

                 Changes in version 1.29.3 (2026-10-09)                 

- The decontX vignette is now a short pointer to the decontX
package's
vignette, which has the full, maintained DecontX tutorial; celda no
longer keeps its own copy. The duplicate pkgdown article
(decontX_pbmc4k) was removed, and the site's DecontX menu entry now
links to the decontX vignette on Bioconductor.

                 Changes in version 1.29.2 (2026-10-09)                 

- Fixed the decontX vignette, which failed to build on Bioconductor
devel with "'groupClusters' not found in 'z': 6". It hard-coded
cluster numbers from decontX's heuristic clustering, which differ
between decontX versions and platforms. The vignette now maps each
cell type to the cluster where its marker genes are detected most
often, as the decontX package's own vignette does.

                 Changes in version 1.29.1 (2026-09-25)                 

- Replaced deprecated ggplot2::aes_string() (and size for line geoms)
with tidy-eval aes() / linewidth in the dimension-reduction, violin,
and perplexity plotting functions, silencing ggplot2 deprecation
warnings (#418).
- Migrated the DecontX implementation to the standalone decontX
Bioconductor package. celda's decontX(), decontXcounts(),
plotDecontXContamination(), plotDecontXMarkerExpression(),
plotDecontXMarkerPercentage(), and simulateContamination() are now
thin, backward-compatible wrappers around the decontX package (added
to Suggests).
- Removed celda's DecontX C++ backend (src/DecontX.cpp) and dropped
the DecontX-only dependencies MCMCprecision, scater, scran, and
DelayedArray from Imports (scater retained in Suggests for the
vignettes).
- Bioconductor release readiness: replaced T/F logical literals with
TRUE/FALSE in findMarkersTree.R (BiocCheck WARNING); added a URL
field and bumped the R dependency to R (>= 4.6.0) in DESCRIPTION;
set LazyData: false; added labels to unlabeled vignette knitr
chunks; and added the required url field to _pkgdown.yml.
- Declared Seurat in Suggests and namespaced the Seurat::Idents() and
Seurat::RunUMAP() calls in findMarkersTree()'s Seurat-object input
path, clearing the remaining R CMD check "no visible global
function" NOTE.

[cellchatr](/packages/cellchatr)
---------
- Initial Bioconductor Submission

[CellMentor](/packages/CellMentor)
----------

                 Changes in version 1.1.2 (2026-05-22)                  

Bug fixes

- calculate_performance() now passes as.factor(true_labels) and
as.factor(clusters) to aricode::NMI(). As of aricode 1.1.0
(currently shipping in Bioconductor 3.24), NMI() routes mixed-type
inputs through sort_pairs(), which calls as.integer() on character
labels — turning cell-type names such as "alpha" into NA and
aborting with "NA are not supported.". Coercing both sides to factor
keeps the computation well-defined and restores compatibility with
the Bioconductor build farm. This was the underlying cause of the
vignette failure that 1.1.1 surfaced but did not fix.

                 Changes in version 1.1.1 (2026-05-19)                  

Bug fixes

- cluster_nmf() now falls back to k-means clustering when the NMF
embedding has too few features or cells for SNN (variable-feature
selection on low-dim embeddings was failing under Seurat 5.x). The
number of k-means / hclust centers is now clamped to a feasible
range, so grid search no longer fails silently on small reference
splits.
- CellMentor() now stops with a clear message if every parameter
configuration in the grid search fails, instead of cascading to a
cryptic switch() error during final-model training.
- evaluate_params() warnings now include the underlying error
message,
so build logs surface the real cause of failures.
- initialize_wh() validates that method is a single character string
and reports a helpful error if not.

Documentation

- Updated citation to the published Nature Communications paper
(Hevdeli et al., 2026, doi:10.1038/s41467-025-67088-7).
- Added inst/CITATION so citation("CellMentor") returns a structured
bibliographic entry.

[cellNexus](/packages/cellNexus)
---------

- Initial Bioconductor Submission

[CGRphylo2](/packages/CGRphylo2)
---------

- Initial Bioconductor Submission 

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.47.4                        

- Fix the bug for as.data.frame(GRangesObj).

                       Changes in version 3.47.3                        

- Update the vignettes to stop trying biomaRt since it is not stable.

                       Changes in version 3.47.2                        

- Fix a bug in featureAlignedExtendSignal when input is
  GAlignmentPairs.

                       Changes in version 3.47.1                        

- Add warning message when the input is not unique for
  genomicElementDistribution.

[ChIPQC](/packages/ChIPQC)
------

                       Changes in version 1.49.3                        

- Add package anchors to \link{...} cross-references that target
  topics in other packages (GenomicRanges, DiffBind,
  BiocParallel), eliminating the ``Found Rd file(s) with Rd
  \link{} targets missing package anchors'' R CMD check NOTE.  9
  man/*.Rd files updated. No code change.

                       Changes in version 1.49.1                        

- Formally deprecate ChIPQCreport() per Bioconductor deprecation
  guidelines: the function now emits a .Deprecated() warning and
  returns invisibly. The Nozzle.R1 CRAN package it depended on is
  no longer available, so HTML report generation is no longer
  possible.  Per Bioconductor policy the function will be made
  defunct in the next release cycle. The original source is
  preserved in inst/legacy/report.r for reference. See
  ?'ChIPQC-deprecated' for replacements (QC accessor methods;
  deeptools for visual reports).

- Remove Nozzle.R1 from Imports: (no longer reachable from any
  executed code path).

- Add zzz.R to the Collate: field (R 4.6 strict Collate
  enforcement).

- Reduce R CMD check NOTE noise: import is from methods,
  geom_point from ggplot2, and multicoreWorkers/snowWorkers from
  BiocParallel; register the data-frame column names used in
  ggplot2::aes() calls (Sample, FRIP, FRIBL, etc.) via
  utils::globalVariables(). No behaviour change.

- Fix \itemize \to \describe in man/ChIPQC.Rd and
  man/ChIPQCsample-class.Rd: items used \item{label}{description}
  (a \describe semantic) inside an \itemize wrapper, producing
  roughly 40 ``Lost braces'' NOTEs. The rendered output also
  improves (proper definition list rather than malformed
  bullets).

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.49.5                        

- test-txdb.R no longer pins the tx_id -> gene mapping of
TxDb.Hsapiens.UCSC.hg19.knownGene. That mapping is annotation data
rebuilt for each Bioconductor release, and tx 70455 changed from
ENST00000487630.1_3 to ENST00000959756.2_2/3554 in BioC 3.24,
failing R CMD check on both Linux builders. Expected values are now
derived from the installed TxDb at test time. (2026-11-12, Thu)

                       Changes in version 1.49.4                        

- The vignette no longer evaluates the KEGG compareCluster() example
at build time. It calls enrichKEGG(), which fetches rest.kegg.jp
live; a transient outage of that endpoint aborted R CMD build on the
Bioconductor Linux builder (nebbiolo1, Bioc 3.23) while building the
vignettes. The example is kept as eval = FALSE with a precomputed
figure (figures/compareCluster_KEGG.png), consistent with the other
cached figures in the vignette. (2026-11-03, Tue)
- annotatePeak() keeps annotation, geneChr/geneStart/geneEnd, geneId
and transcriptId consistent with each other at transcript level
(issue #252). The transcript id of an exon/intron/UTR hit was taken
from names(genomicRegion)[subjectIndex], indexing the unlisted
ranges with the names of the GRangesList; this mostly returned NA
and occasionally a wrong transcript, so a peak could be reported
with the metadata of an unrelated transcript (even on another
chromosome). The ids are now expanded before indexing, and the
Promoter branch writes the feature id reset into anno instead of a
stale local variable, so promoter peaks keep the nearest-TSS
transcript. The same alignment is now applied at level = "gene": the
gene of the overlapping exon/intron/UTR is resolved with TXID2EG()
and used instead of the nearest gene, so
geneId/geneStart/geneEnd/distanceToTSS follow the annotation there
as well. (2026-10-01, Thu)
- seq2gene() no longer fails with $ operator is invalid for atomic
vectors when none of the queried regions overlaps an exon
(getGenomicAnnotation.internal() returns NA then); host-gene
extraction is skipped and the nearest/flanking genes are still
reported (issue #248). (2026-10-01, Thu)
- getTagMatrix()/plotPeakProf2() no longer fail with Error in
cursor:(cursor + seq - 1) : result would be too long a vector when
type = "body" is used with a flank extension shorter than 1kb (issue
#250). The share of bins of a flank is derived from its actual
length (500bp now contributes 5% of the bins instead of being
rounded down to zero), a non-empty flank always gets at least one
column, and the binning loops were replaced by edge-based averaging
of equal-width bins. This changes numbers slightly: the old code
divided the last bin of a section by count-1 instead of count
(about 10% too high for a 10 bp bin), and it put all leftover
positions of a section into its last bin instead of spreading them,
which can move that bin considerably (up to ~40% for a 7 bp
remainder in a 300 bp bin). The x-axis breaks of plotPeakProf2() are
derived from the same column layout. (2026-10-01, Thu)

[See the full NEWS](/packages/release/bioc/news/ChIPseeker/NEWS) or run `news(package = "ChIPseeker")`.


[Chromatograms](/packages/Chromatograms)
-------------

                         Changes in version 1.3                         

Change in 1.3.5

- Add funtions to generate ggplot figures: ggplotChromatograms() and
ggplotChromatogramsOverlay().

Changes in 1.3.4

- Add package logo.

- Add citation to the Analytical Chemistry publication.

Changes in 1.3.3

- Major ChromBackendSpectra performance improvements: chromExtract(),
Chromatograms(spectra), peaksData() and [ no longer re-validate the
wrapped Spectra on every call (which re-stated every backing file),
so they no longer scale with the number of files.

- Improve performance of peaksData() for ChromBackendSpectra with
overlapping chromatogram windows (e.g. from chromExtract()): each
spectrum is aggregated once and shared across the windows it falls
in, instead of once per window. Results are unchanged; the speed-up
grows with the number of overlapping windows.

- Improve performance of peakBoundary(): valleys flanking the apex
are
located by scanning outwards from it rather than scanning the whole
chromatogram. Results are unchanged.

- Order dataOrigin by first appearance when computing the spectra
sort
index, consistent with backendParallelFactor().

Changes in 1.3.2

- Change plotChromatograms() and plotChromatogramsOverlay() to
methods.
- Import compareChromatograms() from ProtGenerics.

Changes in 1.3.1

- Addition of logo and align to BioC 3.23 release.

- Fix ChromBackendSpectra spectraSortIndex test to shuffle spectra
first, avoiding a spurious failure when input data is already
sorted.

[CLAMP](/packages/CLAMP)
-----

- Initial Bioconductor Submission

[clusterExperiment](/packages/clusterExperiment)
-----------------

                       Changes in version 2.33.1                        

- `getBestFeatures` with `DEMethod="limma-voom"` (and functions that
  call it, such as `mergeClusters` and `RSEC` with count data) now uses
  `limma::voomLmFit` instead of `voom` followed by `lmFit`, following
  the limma authors' recommendation of `voomLmFit` over `voom`.
  `voomLmFit` reduces the residual degrees of freedom for genes that
  have zero counts in all of the samples in a cluster, so results for
  such genes can be noticeably different (typically less significant)
  than in previous versions. Because the shared variance estimates are
  affected, the results of other genes in the same analysis can also
  change. If no gene has zero counts in all of the samples of a
  cluster, results are the same as in previous versions. Requires limma
  (>= 3.99.3).

[clusterProfiler](/packages/clusterProfiler)
---------------

                       Changes in version 4.21.6                        

mechanism evidence: correctness and limits

- Pathway and mechanism reports show ledger-derived support,
conflict,
mixed, undetermined and noncomputable evidence before legacy grades,
with attempted/computable denominators and incomplete-coverage
flags; A/B/C rules and discovery scores are unchanged (2026-11-11,
Wed)
- Added an offline paired replay of frozen perturbation responses
across historical/fixed code and 8/4 versus 5/2 caps: it
distinguishes corrected coverage diagnoses from cap sensitivity
without claiming improved accuracy (2026-11-11, Wed)
- adjudicate() retains explicit mixed measurements and does not
select
a programme half when collection directions disagree; incomplete
half-statistics return unknown rather than crashing (2026-11-11,
Wed)
- interpret() preserves a supplied method-normalized effect even when
raw perturbation responses are mixed, without inventing an effect
when it is missing (2026-11-11, Wed)
- adjudicate() reads every measured perturbation row instead of the
first: row order no longer changes the effect; raw response
directions and method-normalized effects are summarized separately,
and opposing normalized effects yield no single effect (2026-11-11,
Wed)
- pathway_evidence() reports whether dataset resolution returned
data,
found nothing usable, or failed; a network outage or a machine
without the suggested gpsadb previously printed "No usable
perturbation data found in GPSAdb" and reason no-dataset, which is a
false statement about the database. Failures are now lookup-failed
with the error text in the dataset_lookup attribute, and a partial
failure is disclosed even when another candidate resolves
(2026-11-11, Wed)
- pathway_evidence() distinguishes missing fold changes from members
that are absent from them: with direction supplied by
pathway_direction or NES and no signed list, the row now says
no-fold-changes instead of claiming the members were not found
(2026-11-11, Wed)
- Subsetting pathway_evidence() by columns keeps the retained
evidence
ledger; a pure column reorder used to strip the attributes while
keeping the class, so mechanism_report() saw no candidate/dataset
evidence (2026-11-11, Wed)
- pathway_evidence() falls back to the other removal spellings GPSAdb
records (CRISPRi, knockout, siRNA, shRNA) when the legacy literal
"knockdown" query returns nothing; STAT1, recorded as CRISPRi plus
knockout only, had been reported as having no data. The search stops
at the first spelling that yields datasets, so methods are not
united; overexpression stays excluded (2026-11-11, Wed)

[See the full NEWS](/packages/release/bioc/news/clusterProfiler/NEWS) or run `news(package = "clusterProfiler")`.


[ClustIRR](/packages/ClustIRR)
--------

                       Changes in version 1.11.6                        

- Paper describing ClustIRR accepted in eLife.

[clustSIGNAL](/packages/clustSIGNAL)
-----------

                        Changes in version 1.5.1                        

- Minor updates to vignette content.
- Added option for using externally generated low embeddings of
smoothed data in the final clustering.
- Reduced the minimum number of centers used in bluster k-means
from 5000 to 3000 in the initial and final clustering steps.
- Improved overall code for better scaling on large datasets.

[CMEnt](/packages/CMEnt)
-----

- Initial Bioconductor Submission

[CONCERTDR](/packages/CONCERTDR)
---------
- Initial Bioconductor Submission

[CorNetto](/packages/CorNetto)
--------

- Initial Bioconductor Submission

[cosmosR](/packages/cosmosR)
-------

                       Changes in version 1.21.2                        

- Remove the legacy DoRothEA installation dependency. Classic
CARNIVAL
preprocessing now requires an explicit TF regulon; use a prepared
CollecTRI regulon for new work.
- Retire load_tf_regulon_dorothea() with a clear migration error.

                       Changes in version 1.21.1                        

- Declare Pandoc (>= 2.8) as a system requirement for rendering
vignettes.

[COTAN](/packages/COTAN)
-----

                       Changes in version 2.13.1                        

Updated package to support R version 4.6 and Roxygen 8.0

                       Changes in version 2.13.0                        

Moved to new version

[crisprDesign](/packages/crisprDesign)
------------

                       Changes in version 1.15.10                       

- Changed default of standard_chr_only (now TRUE).

                       Changes in version 1.15.9                        

- Implemented latest CFD score (2026 CRISPick score).

                       Changes in version 1.15.8                        

- Remove harmless error message when crisprBwa not installed.

                       Changes in version 1.15.7                        

- Another big time optimization for addEditedAlleles.

                       Changes in version 1.15.6                        

- Optimized addEditedAlleles for speed.

                       Changes in version 1.15.5                        

- Added addRepCoordinates function.

                       Changes in version 1.15.4                        

- Added parallelization to addEditedAlleles.

                       Changes in version 1.15.3                        

- Fixed splicing junctions calls for addEditedAlleles.

                       Changes in version 1.15.1                        

- Added addSpecificityProxyScore functionality.

[crisprScore](/packages/crisprScore)
-----------

                       Changes in version 1.17.1                        

- Implemented 2026 version of the CFD score.

[CRISPRseek](/packages/CRISPRseek)
----------

                       Changes in version 1.53.2                        

- updated the stored paired and unpaired OfftargetAnalysis.xlsx files
  for the Bioconductor 3.24 hg19 TxDb and org.Hs.eg.db.

                       Changes in version 1.53.1                        

- replaced mltools::one_hot in deepCpf1() with an internal A/C/G/T
  encoder and removed the mltools, dplyr, and data.table dependencies.

[crossarrayEWAS](/packages/crossarrayEWAS)
--------------

- Initial Bioconductor Submission.

[crumblr](/packages/crumblr)
-------

                        Changes in version 1.5.1                        

- July 6, 2026
- remove dependency on HMP package from CRAN

[ctdR](/packages/ctdR)
----

- Initial Bioconductor Submission

[CyFj11](/packages/CyFj11)
------

- Initial Bioconductor Submission

[CySA](/packages/CySA)
----
- Initial Bioconductor Submission

[CytoMDS](/packages/CytoMDS)
-------

                         Changes in version 1.9                         

CytoMDS 1.9.3

- refactored computeMetricMDS to make the automatic selection of
number of dimensions more efficient, while also handling non
increasing pseudo rsquares, and implementing verbose messages

CytoMDS 1.9.2

(no change - required update of unit tests)

CytoMDS 1.9.1

(no change)

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

[See the full NEWS](/packages/release/bioc/news/CytoML/NEWS) or run `news(package = "CytoML")`.


[CytoPipeline](/packages/CytoPipeline)
------------

                        Changes in version 1.13                         

CytoPipeline 1.13.1

(no change)

[CytoPipelineGUI](/packages/CytoPipelineGUI)
---------------

                        Changes in version 1.11                         

CytoPipelineGUI 1.11.1

(no change)

[cytoviewer](/packages/cytoviewer)
----------

                 Changes in version 1.13.3 (2026-08-31)                 

- improve zoom-preservation

                 Changes in version 1.13.2 (2026-08-31)                 

- fix cell-level/downloads after points integration

                 Changes in version 1.13.1 (2026-06-09)                 

- added point-level visualization

[dagLogo](/packages/dagLogo)
-------

                       Changes in version 1.51.1                        

- update DESCRIPTION.

- remove the set.seed.

[damidBind](/packages/damidBind)
---------

                        Changes in version 1.1.1                        

- NEW: Added additional diagnostic plots (mostly used to generate
supplementary figures in the damidBind manuscript):
plot_catada_mean_variance() and test_weighting_vs_bias_artifact()
- NEW: load_data_peaks() and load_data_genes() now support multiple
normalisation options: quantile normalisation ("quantile"), cyclic
LOESS norm ("loess"), and RPM (reads per million) norm ("rpm". The
default is "none" (no normalisation).
- NEW: a pre_scale option has been added to both load_data_peaks()
and
load_data_genes(), which will perform zero-preserving (uncentred)
scaling (via base R scale(center=FALSE, scale=TRUE)) prior to
normalisation. The default is FALSE.
- NEW: all input metadata and genome annotation metadata are saved in
the DamIDResults S4 object. Basic metadata are displayed with the
default show() method, and all are accessible via the metadata()
accessor.
- NEW: plot_volcano() now takes a new labels option, values of "all",
"highlight", or a character vector containing one or more group
names from the highlight list. Only this group will be labelled on
the plot. The original label configuration options still remain.
- NEW: plot_volcano() supports highlight groups of either loci or
gene
names, including different categories in the same plot. The default
"auto" behaviour matches a highlight group against whichever
definition has the most matches, but use the new highlight_by option
to enforce matching by either "gene_name" or "id" (i.e. locus).
- CHANGED: deprecated the quantile_norm option for load_data_peaks()
and load_data_genes().
- Added a citation to the current bioRxiv preprint on damidBind

[DaparToolshed](/packages/DaparToolshed)
-------------

- Initial bioconductor submission

[dar](/packages/dar)
---

                        Changes in version 2.0.0                        

Changes since dar 1.8.0.

Breaking changes

- R 4.6 or later is required. Differential-abundance recipes
configure
their targets and taxonomic rank with add_model(). The old selectors
and recipe slots are removed. Rebuild saved recipes from their
inputs.
- Use step_aldex2() for ALDEx2 and recipe_qc() for quality control.
step_aldex(), phy_qc() and step_metagenomeseq() are removed. ALDEx2
effects now use the fitted GLM coefficient.
- Recipes use a versioned JSON exchange format. Importing bakes into
an unprepared recipe requires prepare = TRUE.
- otu_table() returns counts with taxa in rows. Sample and taxonomy
accessors return the full tables. ComplexHeatmap, gplots, heatmaply
and UpSetR are optional.

Models and results

- Shared models support covariate adjustment, longitudinal contrasts
and compatible mixed designs. Added LinDA, ALDEx3 and radEmu, with
per-stage engine_args for advanced options.
- tidy_results(), tidy_abundance(), get_execution() and effect_plt()
expose per-method results, aggregated abundances, execution details
and effect estimates.
- Abundance plots support shared or union discoveries and explicit
taxa. Automatic selection requires a contrast when ambiguous.
Consensus outputs retain the contrast and direction of each
discovery.

Fixes and maintenance

- Random-intercept models accept singleton participants alongside
participants with repeated observations. Checks use observed groups
and still reject designs with one group or no repeated observations
(#189).
- Recipes validate steps before execution, run preprocessing in
declared order, reserve result column names and assign unique
automatic IDs (#187).
- DESeq2 uses native poscounts normalization and applies the
significance threshold and result options before shrinkage (#188).
- Corrected contrast handling in MaAsLin3 and LEfSe, robust inference
in ANCOM-BC2 and corncob, and depth-dependent zero handling in
LinDA.
- Five vignettes document the current workflow. Docker Hub
publication
runs separately from GitHub package checks (#190).

[DAssemble](/packages/DAssemble)
---------

- Initial Bioconductor Submission

[daVis](/packages/daVis)
-----
- Initial Bioconductor Submission

[dbSequence](/packages/dbSequence)
----------

- Initial Bioconductor Submission 

[debrowser](/packages/debrowser)
---------

                       Changes in version 1.41.2                        

Bug fixes

- Dropped the sodium dependency, which broke R CMD check on the
Bioconductor Linux builder (nebbiolo2 has no libsodium, so the check
failed with "Package suggested but not available: 'sodium'" before
running any other stage). Auth crypto now uses scrypt and openssl —
both already Imports: of shinymanager, and both present on every
Bioconductor builder.
- Password hashing moved from libsodium argon2id to
scrypt::hashPassword(), the same primitive shinymanager itself
stores.
- Per-user AI key encryption moved from libsodium secretbox to
AES-256-GCM with an explicit encrypt-then-MAC HMAC-SHA256 tag.
R's openssl package neither emits nor verifies a GCM tag, so
authentication is applied explicitly; encryption and MAC subkeys
are domain-separated.
- Existing password hashes and stored AI keys are not readable
across this change — see the migration note below.

Migration

- Password hashes and encrypted per-user AI keys written by 1.41.1 or
earlier use libsodium formats and cannot be read by 1.41.2. Any user
database created during D2 development must be re-created: users
need to reset their passwords and re-enter their AI provider keys.

                       Changes in version 1.41.1                        

Documentation overhaul

- Rewrote the user-guide vignette (vignettes/DEBrowser.Rmd) around
the
modernized six-step Data Prep wizard and the current top-level tabs
(Main Plots, QC Plots, Concordance, Enrichment, Tables). Added
coverage for the features the old guide missed — batch-correction
QC, cross-contrast concordance, GO/KEGG + GSEA enrichment,
reproducibility exports (R / R Markdown / Jupyter / methods
paragraph), bookmarking & sharing, theming and keyboard shortcuts,
and the optional AI interpretation panel.
- Kept the package small. The full illustrated user guide — with all
screenshots and worked examples — lives outside the package at
https://debrowser.readthedocs.io (the UMMS-Biocore/debrowser-docs
repo, served by Read the Docs). The vignette is a concise,
image-free reference that links there, and the README references the
same externally hosted screenshots, so no images are bundled in the
package tarball.
- Rewrote README.md to mirror the modern walkthrough at a shorter
altitude and point to the online guide.
- Removed the unused fontawesome and RColorBrewer entries from
Imports
(icons come via shiny::icon(); no brewer.pal/RColorBrewer:: usage),
clearing the R CMD check "declared Imports should be used" note.

[See the full NEWS](/packages/release/bioc/news/debrowser/NEWS) or run `news(package = "debrowser")`.


[decemedip](/packages/decemedip)
---------

                        Changes in version 1.1.3                        

- Fix deprecated Stan array syntax (real[], int y[N]) in
inst/stan/*.stan for compatibility with stanc3 (BioC 3.23 build).

[decompTumor2Sig](/packages/decompTumor2Sig)
---------------

                 Changes in version 2.99.0 (2026-07-07)                 

- IMPORTANT: major version due to a change of the following default
  parameters:
  1) The default input signatures for readAlexandrovSignatures() have
  been
  removed. The previous default input (URL for COSMIC singatures v2) is
  not
  available anymore (fix was needed). Moreover, version 2 is outdated
  and
  it is preferable to not have any default input rather than another
  version
  that will sooner or later be outdated, too.
  2) The default reference genome and gene annotations, used for
  several
  functions, have been updated from human genome hg19 to hg38.

- Additional features: the functions readAlexandrovSignatures(),
  readGenomesFromVCF() and readGenomesFromMPF() can now also take
  already
  loaded objects as input (VariantAnnotation::VCF or matrix
  representations of
  signatures or MPF data, respectively), rather than only file names or
  URLs
  from which to read the input.

- The provided example VCF and MPF file have been reduced from six
  example
  genomes to three example genomes to reduce the package size and the
  package compilation time. For the latter purpose, also a workflow
  example
  in the vignette has been reworked to speed up its computation.

[decontX](/packages/decontX)
-------

                       Changes in version 1.11.34                       

- The decontPro vignette is rewritten. It uses the cord blood data's
own cell-type labels (no Seurat clustering), explains every output,
shows how to tune the priors with before/after density plots, has a
template for estimating the ambient profile from empty droplets with
background, and covers run time, reproducibility, a FAQ and how to
cite DecontPro. Stan's progress no longer fills the page (#83).
- decontPro() help describes every element it returns and says that
only the ADT totals of background are used, so it must hold empty
droplets only. Its example is smaller and seeded (#83).
- dplyr is no longer suggested (#83).

                       Changes in version 1.11.33                       

- decontPro() arguments are now camelCase, like decontX's: cellType,
deltaSd, backgroundSd, and background for the empty-droplet counts.
The old names cell_type, delta_sd, background_sd and ambient_counts
keep working as permanent aliases; giving both names of an argument
is an error. Positional calls are unchanged. Two kinds of call that
used to work now stop with an error: an abbreviated old name (such
as ambient =, or background = for background_sd, since background
now names the empty-droplet counts), and a positional value after an
alias (such as cell_type = k, 1e-2) (#81).
- New decontPro() arguments iter (default 50000) and seed (default
12345) control Stan's fit; the defaults are the values used before,
so
results don't change. verbose = FALSE hides decontPro's messages
and Stan's console output, which is shown only if the fit fails
(#81).

- decontPro() deletes the temporary CSV that rstan writes, about 100
bytes per ADT x droplet, as soon as the fit ends; before, it stayed
until R exited (#46, #81).
- decontPro() help: deltaSd and backgroundSd are standard deviations
of the priors, not variances, and background is described apart from
the model's background component (#81).

                       Changes in version 1.11.32                       

- The decontX vignette shows the code that maps clusters to cell
types, reuses the clusters for its second run, and splits "Running
decontX" into runs without and with the background matrix, showing
how to import the raw matrix and what decontX needs from it (#43).
New sections cover multiple samples (#45) and a FAQ on data with one
cell type (#43) and doublet detection (#44). It also says how to
cite DecontX.
- The package help page, the README and DESCRIPTION cite both papers.
- Documentation fixes: contamination is a proportion between 0 and 1,
not a percentage; simulateContamination() gives its real delta
default (c(1, 10)) and lists everything it returns; decontXcounts()
has an example (#83).

[See the full NEWS](/packages/release/bioc/news/decontX/NEWS) or run `news(package = "decontX")`.


[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.54.0                        

- Fix scale of the fallback starting values for optim in
  fitNbinomGLMsOptim(), see 1.53.3 below (GitHub issue #134).

- Documentation in ?lfcShrink that type="apeglm" and
  type="normal" fit on the original counts, not the
  outlier-replaced counts, with instructions for using the
  replaced counts if desired.

- New vignette FAQ on random effects: DESeq2 fits fixed effects
  only; pointers to limma-voom with duplicateCorrelation and
  variancePartition's dream for mixed models.

- Added CONTRIBUTING.md.

                       Changes in version 1.53.3                        

- Fix scale of the fallback starting values for optim in
  fitNbinomGLMsOptim() (R/fitNbinomGLMs.R, line 356): the
  initial beta_mat is on the natural log scale and is now
  converted to log2 to match the optim objective. Only
  affects rows where IRLS gave NA or |beta| >= 30, with
  minimal impact on estimates. Thanks to Hien Nguyen Minh
  (@btvn-hien-nm) for reporting, addresses GitHub issue #134.

[DESpace](/packages/DESpace)
-------

                        Changes in version 2.5.1                        

- update citation: DESpace2 published in Bioinformatics
  (doi:10.1093/bioinformatics/btag450)

[DEXSeq](/packages/DEXSeq)
------

                       Changes in version 1.59.1                        

- Fixed a bug where quotes in data files could lead to an error in
  `DEXSeqDataSetFromHTSeq()`. Thanks to multiple reported in issue #15.

[DiffBind](/packages/DiffBind)
--------

                       Changes in version 3.23.6                        

- Restore dba.plotProfile(), which was disabled when the
  profileplyr package (and its soGGi dependency) became
  uninstallable. Profiles are now computed directly from the BAM
  files with EnrichedHeatmap::normalizeToMatrix() and rendered
  with EnrichedHeatmap / ComplexHeatmap, removing the dependency
  on profileplyr and soGGi.  All existing parameters and the ways
  of specifying samples and sites are supported, and the returned
  object is now a RangedSummarizedExperiment rather than a
  profileplyr object.

- Feature annotation of profiled sites is not provided by the new
  backend; the annotate argument is accepted but ignored.

- The dashed vertical line that EnrichedHeatmap draws at the
  profiled position is suppressed, matching the previous plots.
  It can be turned on by passing pos_line=TRUE.

- Scale the sample-name column titles down when the figure is too
  narrow to show them side by side, instead of letting them
  overlap.  Wider figures are unaffected, and an explicit
  column_title_gp still takes precedence.

- Choose a bitmap device type that works on the current system
  when rasterizing the heatmap bodies. ComplexHeatmap requests
  type="cairo", which fails on installations without a usable
  cairo (notably macOS without X11), previously making
  dba.plotProfile() stop with an opaque error. If no type works,
  the heatmaps are drawn un-rasterized instead.

                       Changes in version 3.23.5                        

- Fix bug where a user-supplied minCount in dba.count() was
  silently reset to 0 when summit re-centering was used (the
  default): the recursive re-centering count did not forward
  minCount. Now forwarded.

- Fix a latent typo in pv.vectors() (pv$minCount <- pv$minCount
  read from the just-rebuilt object instead of the saved value),
  which dropped minCount whenever an object was reconstructed. As
  pv.vectors() underlies many operations, minCount is now
  preserved consistently.  Default behavior (minCount = 0) is
  unchanged.

                       Changes in version 3.23.4                        

- Add package anchors to \link{...} cross-references that target
  topics in other packages (GenomicRanges, SummarizedExperiment,
  S4Vectors, DESeq2, GenomicAlignments, GreyListChIP, amap,
  gplots, lattice, systemPipeR), eliminating the ``Found Rd
  file(s) with Rd \link{} targets missing package anchors'' R CMD
  check NOTE. 14 man/*.Rd files updated. No code change.

[See the full NEWS](/packages/release/bioc/news/DiffBind/NEWS) or run `news(package = "DiffBind")`.


[DiPPER](/packages/DiPPER)
------

- Initial Bioconductor Submission

[distinct](/packages/distinct)
--------

                       Changes in version 1.24.2                        

- bug fixed in distinct_test when using 1 core without covariates

                       Changes in version 1.24.1                        

- T_obs added in results

[dnaEPICO](/packages/dnaEPICO)
--------

- Initial Bioconductor Submission

[dominoSignal](/packages/dominoSignal)
------------

                        Changes in version 2.0.0                        

dominoSignal v2 released as part of Bioconductor release 3.24 to
reflect
major (and breaking) changes.

                        Changes in version 1.7.1                        

BREAKING CHANGES

Changes are noted as breaking when a call on the same inputs gives a
different result or a call that used to work now errors. A number of
bug
fixes are included here (prefixed with [BF]), as the fixes alter the
behavior of functions.

Renamed or removed functions, arguments, and data

- Removed the use_clusters argument from create_domino(), as TRUE is
the only option. The counts, z_scores, and clusters arguments are
now required and no longer default to NULL.
- Renamed the min_tf_pval argument of build_domino() to max_tf_pval,
as it as a maximum p-value threshold. Objects built with previous
versions (which store min_tf_pval in build_vars) still pass validity
checks.
- Renamed the first argument of count_linkage() and
test_differential_linkages() from linkage_summary to link_summary to
avoid conflict with linkage_summary() class.
- Removed the test_name argument from test_differential_linkages(),
as
Fisher's exact test is the only option.
- Renamed the return argument of dom_network_items() to which_return
to avoid conflict with return() function.
- obtain_circos_expression() and render_circos_ligand_receptor() are
no longer exported. Use circos_ligand_receptor().
- Renamed elements of the PBMC example dataset: RNA_count_tiny is now
count_tiny and RNA_zscore_tiny is now zscore_tiny.
- Removed mock_linkage_summary(). An example linkage_summary() object
is now provided in the LinkageSummary dataset.

Changed defaults

- Changed the default of remove_rec_dropout in create_domino() to
FALSE to reflect current recommended usage.
- Changed the default of remove_rec_dropout in cor_scatter() to NULL,
which uses the value stored in the domino() object by
create_domino(). Objects created before v1.7.1 do not store this
value and now default to FALSE (previously TRUE).
- The plot_differential_linkages() colors now fun from the gradient
minimum color (red, now customizable) to the maximum color (grey,
now customizable) regardless of stat_ranking. Previously,
stat_ranking = "descending" colored the maximum value red.

Stricter input checks and new errors

- domino() objects now have validity checks: cell-level slots must
share the same cells in the same order, and stored parameters,
clust_de, and cor must be consistent within the object. Objects
created with earlier versions may fail validObject(), in which case
they should be recreated.

[See the full NEWS](/packages/release/bioc/news/dominoSignal/NEWS) or run `news(package = "dominoSignal")`.


[doppelgangR](/packages/doppelgangR)
-----------

                       Changes in version 1.42.0                        

Significant user-visible changes

- Parallel computation now uses the CRAN future framework instead of
BiocParallel. The BPPARAM argument of doppelgangR() is removed;
choose a backend with future::plan(). Supplying BPPARAM gives a
deprecation warning (#29, @lwaldron).

New features

- doppelgangR() accepts SummarizedExperiment objects, and lists of
them, converting them to ExpressionSet internally (#30, @lwaldron).

Bug fixes and minor improvements

- doppelgangR() is faster and uses less memory on large datasets
(#30,
@lwaldron).
- Fix a recycling warning when smoking guns are enabled (#29,
@lwaldron).
- Fix examples, build the vignette when suggested packages are
missing, and expand test coverage (#26, #27, #28, @lwaldron).

[DOSE](/packages/DOSE)
----

                        Changes in version 4.7.4                        

- gseDO(), gseNCG() and gseDisease() now forward the new reportNA
argument from enrichit::gsea_gson() through .... Setting
pvalueCutoff = 1, reportNA = TRUE returns the complete set of tested
pathways, including those whose p-values were set to NA by the
multilevel GSEA engine under unbalanced gene-level statistics
(previously such rows were silently dropped, even with pvalueCutoff
= 1). Requires enrichit >= 0.2.5.9000. (2026-10-04, Sun, #88)

                        Changes in version 4.7.3                        

- gseDisease() now exposes an eps argument and forwards it to
enrichit::gsea_gson(), matching the other GSEA wrappers in the
family (2026-09-22, Tue)
- restore explicit seed support on gseDO(), gseNCG() and gseDisease()
(2026-08-15, Fri)
- seed is forwarded to enrichit::gsea_gson(): set it to a number
(or TRUE for a fixed default seed) to get identical GSEA results
across runs, FALSE (default) draws a fresh seed on each run;
set.seed() before the call still works
- pvalueCutoff is now applied to both pvalue and p.adjust inside
enrichit::gsea_gson() (requires enrichit >= 0.2.2), matching the
historical DOSE double-filtering behavior
- add canonical disease interpretation workflow (2026-07-27, Mon)
- add interpretDisease() and doseInterpretResult as the canonical
seam for disease interpretation outputs, with stable result,
evidence, query, sources, parameters, and explanation components
- support human gene, ranked gene, and gene-set-list inputs,
including multi-query outputs with explicit query_id tracking
- add offline template explanations via explainDisease(..., method
= "template")
- add real-data-backed mouse-model interpretation support from
MGI_DiseaseMouseModel.rpt, HOM_MouseHumanSequence.rpt, and
MGI_GenePheno.rpt, keeping genotype/allele-defined mouse models
as public targets and phenotype profiles as evidence
- enable cross-species interpretation with target = "mouse_model"
and target = "both" when the required contract files are
supplied, preserving mixed disease and mouse_model target types
in one canonical object
- add optional LLM explanation adapter support with evidence-only
payloads, citation validation, unsupported-claim rejection, and
audit metadata hashes

                        Changes in version 4.7.2                        

- store remote NCG data as GSON (2026-07-15, Wed)
- generate NCG.gson.gz directly from the Network of Cancer Genes
download
- record the NCG homepage version and data access date in the GSON
metadata
- publish NCG.gson.gz to gh-pages and load it with
gson::read.gson()
- enforce ontology, species, and Entrez ID consistency (2026-07-15,
Wed)

[See the full NEWS](/packages/release/bioc/news/DOSE/NEWS) or run `news(package = "DOSE")`.


[drawProteins](/packages/drawProteins)
------------

                        Changes in version 4.0.0                        

- Changing tests based on an updated version of ggplot2 to 3.5.0 which
  changes
  the length of the ggplot object.
  ARGUMENT

- Replace size with linewidth in function draw_chain() following a
  depreciation in ggplot.

[DropletUtils](/packages/DropletUtils)
------------

                       Changes in version 1.34.0                        

- Bugfix to read10xCounts() when mtx.class = "SVT_SparseMatrix",
  to generate a valid SVT_SparseMatrix from a Matrix Market file
  containing empty columns.

[drugTargetInteractions](/packages/drugTargetInteractions)
----------------------

                Changes in version 1.21.15 (2026-10-13)                 

- buildGtoPdbDb() now builds from the public GtoPdb bulk download files
  (interactions.tsv and GtP_to_HGNC_mapping.tsv) instead of the
  /services/interactions REST endpoint, which requires a registered API
  key as of 2026 and returned HTTP 401 on any machine without one. No
  key is needed for the default build.

- downloadGtoPdb() and buildGtoPdbDb() gained useServices and apiKey
  arguments for building from GtoPdb's authenticated REST endpoint
  instead; apiKey defaults to the GTOPDB_API_KEY environment variable.
  A
  missing or rejected key is now reported as an authentication failure
  rather than an opaque bfcadd() error.

- The gtp_interactions table and gtoPdbTargetAnnot() now report
  pubmedId,
  the PubMed ID of each interaction's reference, in place of refIds,
  GtoPdb's internal reference IDs. The REST route cannot supply
  pubmedId.

- buildGtoPdbDb() gained a paths argument for building from local
  copies
  of the two bulk files. A small extract of both now ships in
  inst/extdata, so the GtoPdb examples and vignette section run without
  network access.

- A cached SQLite is reused only when its tables and columns are the
  ones the installed version of the package writes. The upstream
  release
  in the cache name is not a sufficient key on its own: BiocFileCache
  persists between runs on the Bioconductor builders, so a database
  built by an earlier package version survives an upgrade there and was
  returned under the same name, after which the query functions failed
  on columns they expect and cannot find. Renaming refIds to pubmedId
  is
  exactly that case. A database with any other layout is now rebuilt.
  Applies to the TTD, Repurposing Hub, GtoPdb and UniChem builds.

- buildGtoPdbDb() always returns an unnamed character(1), as
  documented;
  the cached-reuse and fresh-build paths disagreed before.

- A field named through 'fields' that a source documents but a given
  response does not carry is now returned as an NA column instead of
  raising an error. Previously a transient outage at ChEMBL, PubChem,
  DGIdb or Open Targets turned any fields = c(...) request into a hard
  error, which could fail a vignette build. A name no source documents
  is
  still rejected.

- Repurposing Hub downloads pass the supplemental CA bundle as a
  per-request option rather than setting CURL_CA_BUNDLE in the session.

- getSymEnsUp() loads its EnsDb package namespace instead of attaching
  it
  with require(), and reports a missing package.

- Added runnable examples throughout; most man pages now demonstrate
  their function without network access.

[DspikeIn](/packages/DspikeIn)
--------

                        Changes in version 1.3.2                        

Bug fixes

- validate_spikein_clade(): fixed the could not find function
"DNAStringSet" error in R CMD check (BioC 3.24). Alignment now uses
DECIPHER::AlignSeqs() instead of msa::msa(), whose S4 constructor
calls DNAStringSet() unqualified and fails when Biostrings is not
attached.
- validate_spikein_clade(): terminal branch lengths are now matched
to
the correct tips (previously taken in edge order, not tip order).
- validate_spikein_clade(): the reported clade bootstrap is now the
support of the sample clade's MRCA node rather than the maximum
support in the tree; support is computed on bipartitions because NJ
trees are unrooted.
- validate_spikein_clade(): edge-length annotations are now drawn
with
ape::edgelabels() (they were previously drawn on nodes).

Improvements

- validate_spikein_clade(): input validation, unique tip labels, tree
rooted on the first reference, plots recorded off-screen (no stray
devices or Rplots.pdf), new node_support and branch_table outputs.
- Removed the msa dependency (one fewer package in Imports).

[DuckDBArray](/packages/DuckDBArray)
-----------

- Initial Bioconductor Submission

[DuckDBDataFrame](/packages/DuckDBDataFrame)
---------------

- Initial Bioconductor Submission 

[DuckDBGRanges](/packages/DuckDBGRanges)
-------------

- Initail Bioconductor Submission

[DuckDBSpatial](/packages/DuckDBSpatial)
-------------

- Initial Bioconductor Submission

[easyreporting](/packages/easyreporting)
-------------

                 Changes in version 1.25.2 (2026-11-20)                 

- import the show generic for namespace loading with only stated
  dependencies

- handle missing, NA and empty optional author metadata

- test author parsing independently of the optional distill dependency

                 Changes in version 1.25.1 (2026-11-19)                 

- fix vignette metadata and add a concise S4 show method

- check the optional distill dependency without attaching it

- make documentation examples runnable using temporary report files

- replace bibliography example placeholders and guard HTML rendering on
  Pandoc

[EBImage](/packages/EBImage)
-------

                       Changes in version 4.56.0                        

MINOR CHANGES

- 'rotate()' is now defined as a method for both 'Image' and
  'array' classes, instead of a regular function. This enables
  the definition of a 'rotate()' method in other packages while
  avoiding name clashes.

- the 'image()' generic is no longer defined in EBImage as it is
  already defined in the 'graphics' package.

[ECLIPSE](/packages/ECLIPSE)
-------

- Initial Bioconductor submission.

[edgeR](/packages/edgeR)
-----

                 Changes in version 5.0.0 (2026-11-25)                  

- The genewise C kernels (GLM fitting, deviances, CPM,
  quasi-likelihood adjustment, exact test, diffSplice, etc.) are
  now parallelized with OpenMP. The number of threads is
  controlled by a new `nthreads` argument on the affected
  functions, or globally via `options(edgeR.nthreads=)`, and
  defaults to 1 (serial) so results are unchanged by default.

- The C backend for diffSplice() now exploits the bordered
  block-diagonal (arrowhead) structure of the splice design
  matrix, fitting the genewise null models in time that increases
  linearly, rather than cubicly, with the number of exons per
  gene. Results are unchanged. diffSplice() also gains `maxit`
  and `tol` arguments controlling the genewise GLM fits.

- voomLmFit() has been removed from edgeR and is now provided by
  the limma package.

- New functions binQLFit() and binQLFTest() to fit quasi-binomial
  generalized linear models, to estimate bias-adjusted deviances,
  and to test hypotheses. These functions are closely analogous
  to glmQLFit() and glmQLFTest() but for paired counts and
  differential proportion analyses instead of regular counts and
  differential abundance analyses.

- New functions binFit(), mBinOneWay(), and mBinIWLS() to fit
  binomial generalized models.

- New class DGEBin to store quasi-binomial generalized linear
  model fits.

- New function readBismark2PC() to read Bismark output into an
  edgeR PCList object.

- New class PCList to store paired counts, for example methylated
  and unmethylated reads, for the same genomic regions.

- New functions catchKallistoGene(), catchSalmonGene(), and
  catchOarfishGene() that produce gene-level quantifications from
  kallisto, Salmon or Oarfish output, collapsing the tx
  quantifications to gene-level counts and RTA dispersions. They
  automatically extract gene IDs and annotation from Gencode
  transcript names if they are present, or alternatively accept
  tx to gene mappings via a data.frame.

- New function splitGencodeTxNames() to unpack gene and
  transcript annotation found in Gencode transcript IDs.

- catchSalmon(), catchKallisto(), catchOarfish() and catchRSEM()
  will now detect and unpack Gencode transcript IDs if they are
  present. This provides gene IDs and other annotation that is
  stored in the output object.

- New argument `parent.dir` for catchSalmon(), catchSalmonGene(),
  catchKallisto(), catchOarfish() and catchRSEM() to more easily
  read all the samples available in a specified output directory.
  New argument `sample.dirs` for catchSalmon(),
  catchSalmonGene(), and catchKallisto(), that replaces `paths`
  and allows individual samples to be specified.

[See the full NEWS](/packages/release/bioc/news/edgeR/NEWS) or run `news(package = "edgeR")`.


[eisaR](/packages/eisaR)
-----

                       Changes in version 1.25.1                        

- Rewrite plotEISA to use ggplot2 (and return a ggplot object), set
the aspect ratio of the plot to 1 and add an identity line.
Additional arguments provided via ..., previously passed to plot(),
are ignored.

[ELMER](/packages/ELMER)
-----

                       Changes in version 2.37.3                        

Changes:

- New maintainer: Christopher Smolen. Tiago Chedraoui Silva remains
an
author.
- ELMER is no longer deprecated.

                       Changes in version 2.37.2                        

Changes:

- Fixed writing output files with current readr: the write functions
now use file instead of the removed path argument.
- Probe names are kept when converting probe GRanges to data frames.
- keepStandardChromosomes() is now imported from GenomeInfoDb.

[EMMA](/packages/EMMA)
----

- Initial Bioconductor Submission

[EMTscore](/packages/EMTscore)
--------
- Initial Bioconductor Submission

[enrichmet](/packages/enrichmet)
---------
- Initial Bioconductor Submission

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.99.8                        

- emapplot() / ssplot() no longer require the suggested ggforce
package, and the grouping controls are consolidated on node_label:
- when ggforce is not installed the group outline falls back to
ggplot2::stat_ellipse(). Grouping is on by default in ssplot(),
so the previous unconditional check_installed("ggforce") made
the default call abort on a machine without a merely suggested
package.
- node_label = "category_grouped" replaces group = TRUE (outline
the groups, label the categories); group = FALSE is node_label =
"none". group still works but now warns as deprecated.
- an unknown node_label is an error instead of silently drawing
nothing, and a non-default top_k / target_density that the
selected edge_filter does not use warns instead of being dropped
silently.
- group labels are now always added by emapplot() itself, so "all"
labels the groups exactly as "group" does, and the node_label
values together with their outlines and labels are tabulated in
?emapplot and ?ssplot. (2026-10-09, Fri)
- Align import_fgsea() test suite with the enrichit GSEA leading-edge
contract written for YuLab-SMU/DOSE#46:
- The rank-vector assertions in tests/testthat/test-converters.R
now check, in order: expect_type(result$rank, "integer"),
expect_false(anyNA(result$rank)), expect_true(all(rank >= 0L)).
- The final check carries an in-line comment documenting the
canonical rank == 0L no-signal sentinel — "no usable weighted
in-set signal for this pathway" (N_R == 0, empty gene-set
overlap, or a non-finite running score; matches enrichit >= 0.2.5.9008).
- Together the three assertions guard against silent type drift
(numeric 0 or NA_real_ leaking in for the documented integer
sentinel) and against accidental re-introduction of an
NA_integer_ sentinel that would collapse all(c(NA, 4) >= 0L) to
NA instead of a clean boolean. (2026-10-05, Mon, aligns with
YuLab-SMU/enrichit >= 0.2.6)

                       Changes in version 1.99.7                        

- add opt-in threshold, top_k and adaptive edge filters to emapplot()
and ssplot(); top_k keeps strongest per-term neighbors, while
adaptive targets a configurable graph density without changing the
similarity matrix or dimension-reduction input (2026-09-27, Sun)
- expose min_edge and size_edge directly in ssplot(), document
classical MDS as one available reduction choice, and add an optional
enrichplot_edge_diagnostic attribute plus emapplot_edge_density()
helper for retained-pair density diagnostics (2026-09-27, Sun)
- harden similarity-space selection for term IDs, preserve labels
through one- and two-term reductions, and validate edge/coordinate
controls before plotting (2026-09-27, Sun)

[See the full NEWS](/packages/release/bioc/news/enrichplot/NEWS) or run `news(package = "enrichplot")`.


[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.37.3                        

- Fix security issues reported by the Bioconductor Security Audit:

- `.valueForEnsDb()`: replace single quotes to avoid SQL injection.

- `ensDbQuery()` for `SeqNameFilter`: remove all single quotes from
  sequence
  names to avoid SQL injection; chromosome names are not expected to
  contain
  single quotes.

- `buildWhereForGRanges()`: remove all single quotes form sequence
  names to
  avoid SQL injection.

- `ensDbFromGtf()`: use a parameterized query for `gtfFilename` to
  avoid SQL
  injection.

- `.buildQuery()`: check parameter `order.type` for correct values to
  avoid
  potential SQL injection.

- `.getMetaDataValue()`: use a parameterized query to avoid potential
  SQL
  injection.

                       Changes in version 2.37.2                        

- Fix first security issues reported by the Bioconductor Security
  Audit:

- `fetchTablesFromEnsembl()`: fix potential command injection issue.

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.21.3 (2026-10-31)                 

- generateAmpliconReport now matches both ends of amplicons

                 Changes in version 1.21.2 (2026-09-05)                 

- improvments for plotPatterns

- acknowledgements in README

[epigraHMM](/packages/epigraHMM)
---------

                       Changes in version 1.21.1                        

- epigraHMM now uses csaw's regionCounts in the
epigraHMMDataSetFromBam function to compute read counts

[epimutacions](/packages/epimutacions)
------------

                       Changes in version 1.17.2                        

- Reduced the dependency footprint: removed the unused ensembldb
dependency and moved the optional annotation/visualization packages
(Homo.sapiens, Gviz, rtracklayer, AnnotationHub, ExperimentHub, the
UCSC TxDb packages, the Illumina manifest/annotation packages,
reshape2, purrr, ggrepel and gridExtra) from Imports to Suggests, as
they are already used conditionally via requireNamespace().
- Fixed the CITATION file so it can be read without a
declared-encoding warning (moved to bibentry() with ASCII-escaped
author names).
- Minor BiocCheck NOTE cleanups: labelled all vignette chunks and
replaced sapply()/1:n with vapply()/seq_len() in
plot_epimutations().

                       Changes in version 1.17.1                        

- Fixed plot_epimutations() failing with a non-numeric value column.
- Made the vignette annotation robust to Ensembl/biomaRt network
outages.
- Fixed the license DCF stub and an Rd \itemize NOTE.

[epiSeeker](/packages/epiSeeker)
---------

                        Changes in version 1.1.5                        

- test-txdb.R no longer pins the tx_id -> gene mapping of
TxDb.Hsapiens.UCSC.hg38.knownGene. That mapping is annotation data
rebuilt for each Bioconductor release, and tx 70455 changed from
ENST00000479182.5/91653 to ENST00000798887.1/102724542 in BioC 3.24,
failing R CMD check on both Linux builders. Expected values are now
derived from the installed TxDb at test time. (2026-11-12, Thu)

                        Changes in version 1.1.4                        

- annotateSeq() keeps annotation, geneChr/geneStart/geneEnd, geneId
and transcriptId consistent with each other at transcript level
(issue #252 of ChIPseeker). The transcript id of an exon/intron/UTR
hit was taken from names(genomicRegion)[subjectIndex], indexing the
unlisted ranges with the names of the GRangesList; this mostly
returned NA and occasionally a wrong transcript, so a peak could be
reported with the metadata of an unrelated transcript (even on
another chromosome). The ids are now expanded before indexing, and
the same alignment is applied at level = "gene", where the gene of
the overlapping exon/intron/UTR is resolved with TXID2EG() and used
instead of the nearest gene. (2026-10-01, Thu)
- seq2gene() no longer fails with $ operator is invalid for atomic
vectors when none of the queried regions overlaps an exon
(getGenomicAnnotation.internal() returns NA then); host-gene
extraction is skipped and the nearest/flanking genes are still
reported (issue #248 of ChIPseeker). (2026-10-01, Thu)
- downloadGEObedFiles()/downloadGSMbedFiles() rewrite the ftp:// urls
of gsminfo$supplementary_file to https:// before downloading and
report the underlying download error when a file cannot be fetched
(issue #254 of ChIPseeker). (2026-10-01, Thu)
- Documented the distance definition of flank_gene_distances reported
by annotateSeq(..., addFlankGeneInfo = TRUE): a distance of 0 means
that the peak overlaps the feature range - at level = "transcript"
the feature is the whole transcript, which is why most entries can
be 0 - while non-overlapping peaks get the signed distance to the
feature TSS (issue #235 of ChIPseeker). (2026-10-01, Thu)
- annotateSeq() no longer fails with Error: invalid subscript when
every peak is dropped, which typically happens because none of the
seqlevels of the peaks matches TxDb (chr1 vs NC_000001.11):
.get_distance_to_gene_end() no longer hands the SortedByQueryHits
object returned by follow() for an empty query to features[], and
annotateSeq() reports the seqlevels mismatch together with how to
align the chromosome names (seqlevelsStyle()). (2026-10-01, Thu)

[See the full NEWS](/packages/release/bioc/news/epiSeeker/NEWS) or run `news(package = "epiSeeker")`.


[escape](/packages/escape)
------

                        Changes in version 2.9.1                        

NEW FEATURES

- plaid backend: escape.matrix() and runEscape() gain backend =
c("native", "plaid"). Setting backend = "plaid" routes ssGSEA, GSVA,
UCell, and AUCell to the replaid.* family for large speed and memory
gains. plaid is a Suggests, so the native path is unaffected if it
is not installed. These are approximations, and the gap is wider
than the plaid documentation suggests — measured pooled Pearson
agreement with the native scores was roughly 0.85 (ssGSEA), 0.83
(UCell, AUCell) and 0.73 (GSVA) on a 2000-gene simulation, and lower
on narrower matrices. See ?escape.matrix and the vignette before
substituting one backend for the other.
- Three new methods: method = "PLAID", "singscore", and "scSE". These
have no native implementation and select the plaid backend
automatically.
- backend.args: a named list carrying method-specific tuning (alpha,
tau, rowtf, aucMaxRank, rmax, nsmooth, stats, chunk, removeLog2,
scoreMean). Names are validated against the target function, so
typos error instead of being silently swallowed. Note
backend.args$normalize is plaid's median normalization of the scores
and is unrelated to escape's normalize.
- input.assay: choose which expression matrix to score. "auto"
(default) keeps raw counts for the native backend and picks
log-normalized values for plaid. "counts" and "logcounts" map to the
correct layer for both Seurat and SummarizedExperiment-derived
objects.
- SpatialExperiment support is now explicit and tested, including
spatialCoords(), imgData(), and sample_id preservation through
runEscape().
- Enrichment scores carry an escape.backend attribute recording the
engine, method, and plaid version; runEscape() stores the same
record in metadata() (SCE) or Misc() (Seurat).

BUG FIXES

- escape.matrix(normalize = TRUE) failed on every
SummarizedExperiment-derived input with unable to find an inherited
method for function 'assay' for signature 'x = "NULL"'.
.pull.Enrich() used altExp(sc)[[name]], which indexes the colData of
the first altExp rather than selecting the named one, and returned
NULL. Now uses altExp(sc, name). This also fixes
performNormalization() and performPCA() on SingleCellExperiment and
SpatialExperiment objects. Reported by Cathal King, who also
provided the Xenium test data.
- Underscored gene-set names were silently dropped during
normalization. performNormalization() rewrote every _ to - to match
Seurat's feature-name coercion, which meant HALLMARK_*, GO_*, and
REACTOME_* sets matched nothing for matrix, SingleCellExperiment, or
SpatialExperiment input. Matching is now literal first with the
Seurat-mangled form as a fallback, and per-set scale factors are
aligned to the enrichment columns explicitly so a partial match can
never divide the wrong column.

[See the full NEWS](/packages/release/bioc/news/escape/NEWS) or run `news(package = "escape")`.


[EventPointer](/packages/EventPointer)
------------

                         Changes in version 4.0                         

- Update version of EventPoint to work with BAM files

- Use of bootstrap statistics in novel events

- Update of vignette

[EWCE](/packages/EWCE)
----

                       Changes in version 1.21.1                        

Bug fixes

- Import cbind/rbind from BiocGenerics instead of DelayedArray, which
no longer re-exports them in Bioc 3.24 (fixes install failure,
#106).

Miscellaneous

- Fix "Lost braces in \itemize" Rd NOTE in compute_gene_scores and
drop_uninformative_genes docs.
- Rebuild README.
- README falls back to badges without the Bioc release version when
rworkflows::use_badges() fails to scrape it (#108).

[ExperimentSubset](/packages/ExperimentSubset)
----------------

                Changes in version 1.23.50 (2026-11-10)                 

- The README and the vignette now install the development version
with
BiocManager::install("campbio/ExperimentSubset@devel"). Before,
BiocManager::install("campbio/ExperimentSubset") installed the
default branch, which holds the release version. The README also
links to the website and the workflow article, and shows
Bioconductor build, years, downloads and pkgdown badges next to
R-CMD-check. No change to the package

                Changes in version 1.23.49 (2026-11-09)                 

- The workflow article on the website now uses the 10x 5k PBMC
CITE-seq data set (RNA and 32 antibodies) instead of PBMC 4k, and
keeps Cell Ranger's raw matrix whole: all 6,794,880 barcodes, with
the antibodies in altExp("ADT"). Next to the RNA workflow, the same
object now holds an antibody workflow: decontPro and CLR on the
high-quality cells, stored in the subset's ADT while its counts stay
linked to the main object; the RNA cell-type labels checked against
CD3, CD4 and CD8a protein; and decontPro and CLR again on the T
cells alone, stored with unlink = c("ADT/decontPro", "ADT/clr"). RNA
decontamination now uses the decontX package rather than celda, and
the article keeps cells with up to 15% mitochondrial counts rather
than 10%, since this data set has more. No change to the package
(#93)

                Changes in version 1.23.48 (2026-11-08)                 

- The package code, tests and vignette now pass make lint with no
findings: spacing, braces, quotes, assignment, indentation and line
length follow the package's lint settings, which now also accept the
names of the exported constructors, the class generators and the
use.names and deparse.level arguments that methods take from their
generics. subsetSummary() loops with seq_along() and createSubset()
drops an internal variable that always equalled parentAssay. No
change in behavior (#37)

[See the full NEWS](/packages/release/bioc/news/ExperimentSubset/NEWS) or run `news(package = "ExperimentSubset")`.


[exploreSE](/packages/exploreSE)
---------
- Initial Bioconductor Submission

[ExpressionAtlas](/packages/ExpressionAtlas)
---------------

                        Changes in version 2.5.0                        

- Enable support for Single-cell Expression Atlas.

- Improve data searches.

- Add visualisation functions.

- Add two new vignettes.

[extraChIPs](/packages/extraChIPs)
----------

                       Changes in version 1.17.2                        

Bug Fixes

- Improved memory consumption for dualFilter()
- Deprecated keep.totals & set default bin size to be 10kb

                       Changes in version 1.17.1                        

Bug Fixes

- Better handled label background colours for plotSplitDonut
- Updated fitAssayDiff() for compatability with edgeR v4.0.0

[Fancy](/packages/Fancy)
-----

- Initial Bioconductor Submission

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

[fenr](/packages/fenr)
----

                       Changes in version 1.10.3                        

- A deprecation note was added about the dataset argument in
fetch_go() function due to Biomart being retired.

                       Changes in version 1.10.2                        

- A permanent fix to the species naming issue in Gene Ontology. Now
fenr operates with the new names (e.g. HUMAN-uniprot) but also
accepts the legacy names (such as goa_human) and translates them
using a build-in table. Users of legacy GO names should not notice
any difference. See Section "A note on Gene Ontology species naming"
in the package vignette for more details.
- A fix to fetch_kegg_species() due to a phased out API endpoint.

                       Changes in version 1.10.1                        

- Replaced fetch_go_species() - species designation discovery - by
scraping geneontology.org website with a fixed saved table. This is
because gene ontology completely changed the species website and the
naming convention. This is a temporary fix.

[fishash](/packages/fishash)
-------

- Initial Bioconductor submission.

[fishpond](/packages/fishpond)
--------

                       Changes in version 2.19.5                        

- loadFry() reads the alevin-fry count matrix with a new C reader. It
needs much less memory and time, and the output does not change.
Duplicate (row, column) entries, MatrixMarket files that are not
"coordinate", and NA values now give an error.

                       Changes in version 2.19.4                        

- loadFry() now reads the 10x Flex (multiplex) output of simpleaf
multiplex-quant and alevin-fry quant-probe. For a multiplex run, the
colData gets sample_name and cell_barcode, the two parts of the row
label <sample_name>_<cell_barcode>, from featureDump.txt. Without
featureDump.txt, a multiplex run gives a warning and gets no new
columns. Other outputs do not change.
- loadFry() has two new arguments. featureLevel = "probe" reads the
cell by probe matrix of --probe-mtx. probeRef gives the 10x
probe-set CSV of the run (by default probe_set.csv in fryDir); the
rowData is the CSV as it is, with the rows of the matrix probes.

                       Changes in version 2.19.1                        

- The swish vignette now loads the pre-built macro_txp_se and
macro_tx2gene objects from the macrophage package (>= 1.29.1), as
macrophage no longer ships the Salmon inferential replicates.
- The alevinEC() test is now skipped, as the alevin bfh.txt files
were
removed from tximportData (>= 1.41.1) to reduce package size.
alevinEC() currently has no test coverage.

[flowGate](/packages/flowGate)
--------

                 Changes in version 1.13.1 (2026-07-02)                 

- Implemented fixes due to the ggplot2 version 4 changes, with the
update ggcyto syntax changes
- Where appropiate, switched from import to importFrom tags
- Implemented gs_apply_gate_check to allow adjustment existing gates
- Additional unit tests to test suite

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

[fraq](/packages/fraq)
----

                        Changes in version 1.1.3                        

- Fix fraq_concat() discarding output written to a .fifo path when
the
reader is not keeping up. Buffered bytes are now drained before the
writer is closed, matching what fraq_run() already did.
- Fix a file descriptor leak when closing a .fifo writer failed.
- Detect forked R processes with RcppParallel::isProcessForkedChild()
instead of probing parallel::isChild(), so the serial fallback
triggers regardless of which package performed the fork. This raises
the minimum RcppParallel version to 6.1.1.
- PlainWriter and FraqfFileWriter now tolerate being closed twice,
like the gzip and zstd writers already did.

                        Changes in version 1.1.2                        

- Fix .zst outputs gaining a spurious trailing empty frame:
ZstdWriter
now flushes with ZSTD_flushStream and only ends the frame in
close().
- Close .fraq/.mem writers explicitly on the success path, so write
errors at close time are reported instead of silently discarded.

                        Changes in version 1.1.1                        

- Move ShortRead from Imports to Suggests and check for it at runtime
in ShortRead bridge helpers.
- Add TBB flow graph compatibility helpers so the package builds
against both classic TBB and oneTBB (fixes a build failure caused by
an RcppParallel/TBB update that removed tbb::flow::source_node and
made limiter_node::decrement private).
- Add fraq_run_r() for block-oriented R kernels that run on the main
R
thread while file I/O and compression can continue on nthreads TBB
worker threads.
- Add serial fallback paths for C++ and R kernels when nthreads = 1,
and force the serial path inside forked R processes to avoid using
TBB after fork().
- Update bundled zstd to 1.5.7.

[FRASER](/packages/FRASER)
------

                        Changes in version 2.10                         

- Import `cbind` from `BiocGenerics` instead of `DelayedArray` as it is
  not
  reexported anymore in BioC 3.24

[GCPtools](/packages/GCPtools)
--------

                        Changes in version 1.4.0                        

- Update examples and documentation to use a truly public Landsat
dataset (gs://gcp-public-data-landsat/...) instead of the restricted
genomics dataset.
- Update gsutil_pipe() examples to read structured text data with
readLines() rather than read.csv().
- Fix argument passing and shell quoting in gsutil_rsync() for the
exclude pattern.
- Fix shell quoting in gsutil_pipe().

[GDSAnnotator](/packages/GDSAnnotator)
------------

- Initial Bioconductor Submission

[gdscloud](/packages/gdscloud)
--------

- Initial Bioconductor Submission

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.50.0                        

NEW FEATURES

- new storage types "rle.uint8", "rle.int32" (or "rle.int"),
  "rle.int64",
  "rle.float32" (or "rle.float"), "rle.float64" (or "rle.real" and
  "rle.double") and "rle.string": run-length encoded arrays.  The
  values
  are stored as runs of one repeated value, runs of an arithmetic
  progression (a first value and a step, so that `seq(0, 1, by=0.01)`
  takes 18 bytes) and blocks of up to 128 values kept as they are, with
  a
  table of checkpoints for random access

- add.gdsn() and append.gdsn() accept an object of class "rle" returned
  by base::rle(), for a node of any storage type.  The elements are
  written through a buffer of a fixed size, so the memory needed does
  not
  depend on the number of elements

- new argument 'allow.block' in `append.gdsn()`: if TRUE (the default),
  a GDS node given as 'val' can be appended by copying its encoded data
  as a block when it has the same storage type as the target node;
  'allow.block=FALSE' decodes every element and encodes it again, so
  that the appended data are compressed according to the target node
  (e.g., with its own compression level)

SIGNIFICANT USER-VISIBLE CHANGES

- add.gdsn(node, name, val) without `storage` stores an object of class
  "rle" in a run-length encoded array ("rle.int32", "rle.float64",
  "rle.string" or "rle.uint8"), and read.gdsn() returns an ordinary
  vector.  It used to create a folder of `lengths` and `values`, which
  `storage="list"` still does

UTILITIES

- update XZ to v5.8.4 from v5.2.9 (liblzma is now under the BSD Zero
  Clause License); GDS files stay compatible in both directions

- update zlib to v1.3.1 from v1.2.13: the bundled zlib sources are the
  unmodified upstream release, so the local changes for C23 compilers
  are no longer needed

- `append.gdsn()` copies variable-length strings ("string", "cstring"
  and
  their 16-/32-bit variants) as a raw block when the source is a GDS
  node
  of the same storage type, instead of decoding and re-encoding one
  string at a time

- `append.gdsn()` copies a sparse array ("sp.int8", "sp.real64" and the
  other sparse types) as a raw block of records when the source is a
  GDS
  node of the same storage type, instead of decoding every element and
  encoding it again.  The cost now follows the size of the encoded data
  instead of the length of the array: appending a node of a hundred
  million elements holding a hundred thousand non-zero values is about
  87 times faster

[See the full NEWS](/packages/release/bioc/news/gdsfmt/NEWS) or run `news(package = "gdsfmt")`.


[gemma.R](/packages/gemma.R)
-------

                        Changes in version 4.0.0                        

- get_annotation_children and get_annotation_parents are added.
get_child_terms is deprecated in favor of get_annotation_children

[geneClusterPattern](/packages/geneClusterPattern)
------------------

- Initial Bioconductor Submission

[GeneNetworkBuilder](/packages/GeneNetworkBuilder)
------------------

                       Changes in version 1.55.1                        

- Fix the warning in BiocCheck.

[GenomAutomorphism](/packages/GenomAutomorphism)
-----------------

                       Changes in version 1.15.2                        

- Fix Bioconductor devel build failure from Biostrings >= 2.81.3
(Bioconductor 3.24). Adapt to Biostrings >= 2.81.3
(Bioconductor 3.24), where DNAMultipleAlignment,
AAMultipleAlignment, and readDNAMultipleAlignment moved to the new
MultipleAlignment package. Import from MultipleAlignment, declare
the S4 class via importClassesFrom(), add explicit
DNAMultipleAlignment methods, and regenerate bundled MSA datasets.
This fixes lazy-loading and example failures on the Bioconductor
build system. (Not applicable to Bioconductor 3.23 /
Biostrings 2.80.1.)

                       Changes in version 1.15.1                        

- Update math documentation notation and BiocCheck vignette fixes

                       Changes in version 1.15.0                        

- Fix a bug in the internal function 'weighted_manhattan' (used by
'codon_dist', 'codon_dist_matrix', and 'aminoacid_dist') where, for
\eqn{group = "Z4"}, the second coordinate weight ('w[2]') was
silently omitted and 'w[3]' was applied twice instead, contradicting
the documented weighted Manhattan distance formula. This changes the
numeric distances returned for \emph{group = "Z4"} (the default).
The bundled 'cdm_z64' dataset was regenerated accordingly.

- Documentation updated

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.50.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- getChromInfoFromEnsembl() queries the old Ensembl FTP.
  Following the total revamp of the Ensembl website and FTP servers
  that
  occurred in Summer 2026, getChromInfoFromEnsembl() now queries the
  old
  (legacy) Ensembl FTP sites that got moved to the EMBL-EBI FTP server.
  This is temporary until getChromInfoFromEnsembl() gets refactored to
  query the new live Ensembl website.

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.50.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Move SummarizedExperiment from Depends to Suggests.
  This reduces the number of hard dependencies (direct and indirect)
  from
  29 to 18, which in turn significantly speeds up
  library(GenomicAlignments)
  (speedup is approx. 2.2x).

DEPRECATED AND DEFUNCT

- Defunct CIGAR utilities and sequenceLayer(). These were deprecated
  in GenomicAlignments 1.45.5 (BioC 3.22).

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.66.0                        

- No significant changes in this version.

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.66.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Some changes to what as.data.frame.GenomicRanges() and
  as.data.frame.GPos() return.

[geomeTriD](/packages/geomeTriD)
---------

                        Changes in version 1.7.2                        

IMPROVEMENT

- Fix the missing value sections in man pages.

                        Changes in version 1.7.1                        

IMPROVEMENT

- update the curver connection calculation methods.

[GEOquery](/packages/GEOquery)
--------

                       Changes in version 2.82.0                        

Summary of changes since 2.80.0 (Bioconductor 3.23). Entries for
each 2.81.x development version below have the details.

Breaking changes

- getGEO() returns SummarizedExperiment objects by default for GSE
Series Matrix records (previously ExpressionSet). Pass returnType =
"ExpressionSet" for the old behavior, or convert existing results
with as_SummarizedExperiment(). SOFT-format results (GSE, GSM, GPL,
GDS objects) are unchanged (#168).

New features

- Single-cell data. getGEOSingleCell() and readGEOSingleCell() load
GEO single-cell supplementary files (10x Matrix Market, 10x HDF5,
AnnData .h5ad, and .rds Seurat/SingleCellExperiment objects) into
SingleCellExperiment, or Seurat with as = "Seurat". They work from a
GSE or a GSM, find files at the series or sample level, combine
samples with by = "sample", "platform" or "all", and attach each
sample's GEO metadata to colData. geoSingleCellManifest() and
geoSingleCellUnits() show what a study contains before downloading
(#158, #190, #195-#197, #210).
- Raw sequencing data. geoToSRA() maps a GSE or GSM to its SRA study,
experiment, run and sample accessions (#219).
- Private records. getGEO() and getGEOfile() take a reviewer access
token for embargoed GEO records (#154).
- Persistent cache. options(GEOquery.cache = TRUE) caches downloads
across sessions with BiocFileCache; geoCache() and clearGEOCache()
manage it. Off by default (#171).
- Faster search. searchGEO() uses NCBI's batched JSON esummary and
returns stable, structured columns (#223).
- Sturdier downloads. Files stream to disk, retry transient HTTP and
connection failures, and respect options(GEOquery.download.timeout)
(default 300 s) (#147, #173).
- Typed errors. Failures raise geoquery_error subclasses
(geoquery_download_error, geoquery_private_accession,
geoquery_parse_error, geoquery_bad_accession,
geoquery_ncbi_challenge), so they can be handled with tryCatch()
(#170, #230).
- getGEO() gains encoding for non-UTF-8 records (#148), and
getGEOSuppFiles() gains quiet (#68).

Deprecated

- GDS2MA() is deprecated in favor of GDS2eSet(). limma moves to
Suggests (#217).

Bug fixes

- getRNASeqData() works again for human series: NCBI now blocks the
shared gene annotation download with a reCAPTCHA page, so the
annotation is read from the GEO FTP mirror. Mouse series still get a
clear geoquery_ncbi_challenge error (#230, #235).
- Series Matrix records with duplicate feature IDs (e.g. GSE136400)
no
longer fail (#98).

[See the full NEWS](/packages/release/bioc/news/GEOquery/NEWS) or run `news(package = "GEOquery")`.


[gesel](/packages/gesel)
-----

- Initial Bioconductor Submission

[GExPipe](/packages/GExPipe)
-------

- Initial Bioconductor Submission

[geyser](/packages/geyser)
------

                        Changes in version 1.5.1                        

Document pandoc dependency. Increment version to match (?)
bioconductor.

[ggmsa](/packages/ggmsa)
-----

                       Changes in version 1.19.2                        

- AlleleMono glyphs updated to Allele Mono 1.100: thinner horizontal
strokes, overshoot on round letters, tapered diagonal vertices and
optical centring in the cell. The character set and the
glyph_polygon format are unchanged (2026-10-30, Fri)
- seqlogo(), geom_seqlogo(), geom_seed() and geom_GC() no longer
error
when ggmsa is loaded but not attached, e.g. when called as
ggmsa::seqlogo() from another package. Options set before
library(ggmsa) are now kept. issues#78 (2026-10-24, Sat)
- added a monospaced font AlleleMono (font = "AlleleMono"), designed
for reading DNA, RNA and protein sequences: easily confused letters
such as C/G, H/M/N/W, D/Q and O/0 are drawn apart. The glyph
polygons are generated from the Allele Mono font (SIL Open Font
License 1.1, https://github.com/aistanbulresearch/allele-mono)
(2026-10-21, Wed)
- Fixed custom_color: the color column was fetched with partial $
name
matching, so a data frame carrying an extra column (e.g. column)
silently dropped the color column; candidate names are now matched
explicitly and an unsupported column raises a clear error.
consensus_views = TRUE also forgot to pass custom_color to the color
scheme, so custom colors were silently ignored. issues#21 issues#54
issues#69 (2026-09-21, Mon)
- Fixed the missing fill legend: geom_msa() mapped fill = ~I(color),
which prevented scale_fill_manual from training any level, so
show.legend = TRUE produced nothing. It now maps fill = ~color; tile
colors are unchanged. issues#68 (2026-09-21, Mon)
- position_highlight no longer crops the alignment data. Highlighted
columns are drawn as a separate geom_tile layer instead of being
squeezed into the main layers, so facet_msa(), geom_seqlogo(),
geom_seed(), geom_GC() and geom_msaBar() added afterwards still see
the full alignment. issues#75 issues#76 (2026-09-21, Mon)
- facet_msa() now iterates over all layers and facets each one with a
shared origin, rather than assuming the alignment always lives in
layers[[1]] and layers[[2]]. issues#75 issues#76 (2026-09-21, Mon)
- Background tiles and highlight boxes are deduplicated per alignment
cell, so each cell is drawn once instead of once per glyph vertex.
A 3600-cell alignment drops from 392277 to 7260 <path> objects in
SVG output, which is what made Illustrator crawl; rendered ink is
unchanged. issues#74 (2026-09-21, Mon)
- geom_msaBar() is composable again: it used to hand back an <aplot>
object, which ggplot2 4.x has no + method for, so every later +
either errored (facet_msa()) or silently returned NULL (theme*(),
labs(), scales, further layers). It now returns a patchwork, and
geom_msaBar() and facet_msa() work in either order — the
conservation bar is split with the same field, so its columns stay
aligned with the tiles. issues#49 issues#72 (2026-09-21, Mon)

[See the full NEWS](/packages/release/bioc/news/ggmsa/NEWS) or run `news(package = "ggmsa")`.


[ggsc](/packages/ggsc)
----

                       Changes in version 1.11.4                        

- add sc_dim_geom_ring() and sc_dim_geom_ring_label() to draw
annotation rings around a dimension reduction plot without rescaling
the embedding: the cells are laid out in sectors (by default the
clusters used to colour the points) and each ring shows a metadata
variable as coloured arcs, so the cluster, the treatment and any
other track can be read around the same circle. Each metadata ring
gets its own legend (legend, legend.title); a ring that repeats the
variable already colouring the points stays silent by default so the
point legend is not duplicated, and multiple legend-bearing rings
use ggnewscale (2026-11-08, Sun)
- add sc_dim_fraction(), the bar-plot companion of the rings for the
aplot-combined style: it shows the same cell-type composition beside
the embedding (e.g. aplot::plot_list()) and reuses the embedding's
colour scale so a cell type has one colour across both views
(2026-11-08, Sun)

                       Changes in version 1.11.3                        

- sc_dim_geom_ellipse() and sc_dim_geom_label() now report that the
group column could not be determined when the plot has no
non-numeric column to group by (e.g. sc_dim() without colLabels());
the check indexed the plot data with the missing column name and
failed with attempt to select less than one element in get1index
before the message could be displayed (2026-10-10, Sat)
- sc_dim_geom_ellipse() only passes level to layer functions that
accept it, so the documented geom = ggforce::geom_mark_hull no
longer warns about an unknown parameter (2026-10-09, Fri)
- sc_dim_geom_label() places labels at the true group centre; the
reference point was averaged from a sampled ellipse outline, which
was offset from the centre and drifted with the ellipse level
(2026-10-09, Fri)
- sc_dim_geom_label(data = ) now reports that data must not be
supplied, instead of an unrelated message about the label mapping
(2026-10-09, Fri)

                       Changes in version 1.11.2                        

- compatible with new S4Vectors (2026-09-10, Thu)

                       Changes in version 1.11.1                        

- compatible with Seurat V5 (2026-05-17, Sun)

[ggtree](/packages/ggtree)
------

                        Changes in version 4.3.2                        

- fortify.tanglegram() sets tbl_tree once in the class vector instead
of twice, which broke any exact class() comparison (2026-10-27)
- tree coordinates are computed by the new treelayout package, which
ggtree now depends on; the layout code that used to live here is in
treelayout (2026-10-25)
- speed up the equal_angle layout and fortify() on large trees:
layoutEqualAngle() no longer calls offspring() once per node to
count subtree tips, and accumulates the geometry in numeric vectors
instead of assigning into a tibble row by row. The output is
unchanged; a 2000-tip tree drops from 20.6 s to 0.2 s.
getNodeEuclDistances() is vectorised (5.5 s to instant at 4000 tips)
(2026-10-24, Sat)
- regression tests comparing the fast tip count against offspring()
on
binary, polytomous and reordered trees (2026-10-24, Sat)

                        Changes in version 4.3.1                        

- new layout = "tidy" for ggtree() and geom_tree(), implementing the
non-layered tidy tree layout of Ploeg (2014, doi:10.1002/spe.2213),
which arranges branches to minimise vertical displacement
(2026-09-14, Mon, #710)
- this layout may place two tips on the same y; the tree drawing
stays valid, but panel-aligned layers (gheatmap(), msaplot())
now warn in that case
- gheatmap() and msaplot() now warn when the tips of the tree are not
on a unique, evenly spaced y, instead of silently drawing
overlapping rows (2026-09-14, Mon, #710)
- aes(subset = ) keeps the environment of the original aes() call
when
ggtree merges it with its own condition (& isTip, & !isTip, &
node==parent), so helper variables defined in the calling frame now
resolve instead of failing with object 'x' not found (2026-09-14,
Mon, #705)
- regression tests for aes(subset = ) in
geom_tiplab()/geom_tippoint(), including string and factor
comparison and the treedataList + facet_wrap() workflow (2026-09-14,
Mon, #705, #699)
- regression tests pinning the ggplot2 4.0 compatibility of
gheatmap()
(no data.frame leaking into the plot mapping, stacked heatmaps via
ggnewscale) and of geom_tiplab(align = TRUE) (2026-09-14, Mon, #686,
#700, #707)
- ggtree(tr, size = ) / geom_tree(size = ) no longer emit the
ggplot2 3.4.0 "size aesthetic for lines" deprecation warning; the
parameter is translated to linewidth inside stat_tree(), so the
documented idiom keeps working (2026-09-14, Mon, #684)

[See the full NEWS](/packages/release/bioc/news/ggtree/NEWS) or run `news(package = "ggtree")`.


[ggtreeExtra](/packages/ggtreeExtra)
-----------

                       Changes in version 1.23.1                        

- introduce text.height in axis.params parameter in geom_fruit.
(2026-06-04, Thu)

                       Changes in version 1.23.0                        

- Bioconductor 3.23 released, and Bioconductor 3.24 (devel) bump.
(2025-04-29, Wed)

[glmGamPoi](/packages/glmGamPoi)
---------

                  Changes in version 1.25 (2025-05-26)                  

- Change license from GPL-3 to MIT.

[GOfuncR](/packages/GOfuncR)
-------

                       Changes in version 1.33.5                        

OTHER

- GOfuncR is no longer deprecated

                       Changes in version 1.33.3                        

OTHER

- New maintainer: Christopher Smolen

                       Changes in version 1.33.2                        

BUG FIXES

- Fix installation with R >= 4.6: remove the unused internal function
  unlock_environment(), which used the non-API entry points ENVFLAGS
  and SET_ENVFLAGS (see https://github.com/sgrote/GOfuncR/issues/12)

- Add package anchor to the vioplot link in ?plot_anno_scores

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.39.5                        

- fix combineScores() handing back the maximum instead of the
best-match average when one of the two term sets holds a single
term. The early return for a one-row or one-column score matrix
applied max() to every method but avg(), so mgoSim(), geneSim(),
clusterSim() and mclusterSim() scored a single-term set 1 against
any set whose annotation contained that term, however many terms the
other side carried. The branch is the common route rather than an
edge case: drop = FALSE in subsetCombine() (#48) sends every
single-term gene through it, and 25% of the genes annotated in MF in
'org.Hs.eg.db' have exactly one term. Such a matrix now follows the
same definition as a larger one -- the row and column maxima
average, (max + sum) / (n + 1) for a single row or column -- which
leaves max and rcmax where they were, since both do equal the
maximum there. The result is identical for all five methods
on 12,500 two-dimensional score matrices; in a 24-gene mgeneSim()
run 134 pairs that read 1.0 took their proper value and no pair
between multi-term genes moved. The combine documentation now states
what BMA averages (#52, 2026-11-21, Sat)

                       Changes in version 2.39.4                        

- let godata() take a 'GSON' object of GO annotations as annoDb,
closing the ## to be supported branch of #44. The branch was already
present but could never fire -- it tested inherits(annoDb, 'gson')
while the S4 class is GSON -- so passing e.g. the output of
clusterProfiler::gson_GO() fell through to check_goAnno() and failed
with "annoDb as a data.frame should contains 'GO' and 'ONTOLOGY'
columns". GSON carries no ontology information, so ONTOLOGY is
recovered from 'GO.db' term by term; against
AnnotationDbi::toTable(GOTERM) on an 80,000-row sample of
gson_GO(org.Hs.eg.db, ont = "ALL") the resulting BP/CC/MF split
(51058/16942/12000) agrees with GO.db for every term. GO IDs held
without the GO: prefix are recognised, terms absent from 'GO.db' are
dropped with a warning, and a gene set containing no GO IDs at all
is now rejected instead of quietly yielding an empty dataset. 'gson'
goes to Suggests rather than Imports, because reading the slots
autoloads the namespace and the route should cost nothing to users
who do not take it (2026-10-08, Thu)

                       Changes in version 2.39.3                        

- take the set of terms to score in computeIC() from the installed
'GO.db' instead of from the gotbl table shipped in data/. That table
was last rebuilt in 2021 and 524 BP terms added to GO.db since then
are missing from it, so goSim() and geneSim() returned NA for them
even though their descendants carry annotations and their IC is
computable -- 167 of them are annotated in org.Hs.eg.db, plus 45 for
CC and 367 for MF (#33). The missing terms also biased every other
IC slightly, since they count towards the denominator but can never
contribute to an ancestor's descendant sum; about 0.3% of terms were
overestimated by more than 0.5 (2026-09-22, Tue)

[See the full NEWS](/packages/release/bioc/news/GOSemSim/NEWS) or run `news(package = "GOSemSim")`.


[GPlinksR](/packages/GPlinksR)
--------

- Initial Bioconductor Submission

[graphite](/packages/graphite)
--------

                 Changes in version 1.59.1 (2026-11-22)                 

- Updated all pathway data.

[grayleafspotr](/packages/grayleafspotr)
-------------

- Initial Bioconductor Submission

[growkar](/packages/growkar)
-------
- Initial Bioconductor Submission

[GSEAlens](/packages/GSEAlens)
--------

- Initial Bioconductor Submission

[GSVA](/packages/GSVA)
----

                         Changes in version 3.0                         

USER VISIBLE CHANGES

- Enrichment scores may differ from those of earlier versions of GSVA
  in the following cases: the GSVA method with sparse data, because of
  the new handling of zeros, and with missing values, because random
  walks now run only over the non-missing genes of each column; the
  PLAGE method, whose sign is now set as described below; and scores
  calculated on GPU hardware, which uses single precision.

- The GSVA algorithm has been updated to handle zeros in sparse data in
  such a way that enrichment scores are stable, and the scoring step is
  now between one and two orders of magnitude faster than in previous
  versions of the package.

- There is a new row-normalization method in the GSVA algorithm, the
  centered log ratio (CLR) transformation, which is computed in blocks
  of columns and scales better with the size of the input data than the
  original ECDF method. The row-normalization method can now be
  controlled through the argument rowNorm=c("auto", "ecdf", "clr",
  "none") of the parameter constructor function 'gsvaParam()', where
  the default "auto" selects CLR for sparse input processed with
  sparse=TRUE whose nonzero values are positive, and ECDF otherwise,
  while "none" skips row normalization.

- Functions 'gsvaRanks()' and 'gsvaScores()' have been deprecated and
  replaced by 'gsvaRowNorm()', 'gsvaColRanks()' and 'gsvaColScores()',
  which perform separately the three steps of the GSVA algorithm (row
  normalization, column rank transformation and gene set scoring), and
  take two new parameters 'first' and 'last' that restrict calculations
  to a range of rows or columns.

- A new "average" method has been added, through the parameter
  constructor function 'avgParam()', which currently calculates the
  arithmetic mean of expression values per gene set and column,
  handling missing values as in GSVA or ssGSEA.

- Z-scores and PLAGE scores are now calculated by column to help keep
  the memory consumption within the budget given in the 'maxmem'
  parameter. The sign of PLAGE scores may change with respect to
  earlier versions of GSVA, since instead of being the arbitrary
  default of the SVD, it is now set so that scores correlate positively
  with the mean expression of the gene set. PLAGE uses a fast power
  iteration for gene sets with more than 1,000 genes, and warns about
  gene sets too large for the memory budget.

- Several improvements in ssGSEA make it faster and consume less
  memory, including the matching of gene sets to ranked genes and
  computing R^alpha by blocks of columns.

- ssGSEA, Z-score and PLAGE were ignoring the 'ondisk' parameter and
  now take it into account.

- Gene identifiers in gene sets are mapped all at once, instead of one
  gene set at a time, which makes it substantially faster with large
  collections of gene sets whose gene identifiers differ from the type
  used in the input expression data.

- Rows whose values are all missing are now discarded during row
  filtering.

[See the full NEWS](/packages/release/bioc/news/GSVA/NEWS) or run `news(package = "GSVA")`.


[gVenn](/packages/gVenn)
-----

                       Changes in version 1.99.0                        

New features

- Add a mode argument to computeOverlaps() controlling how genomic
intervals are made non-redundant before they are classified. The
default, mode = "reduce", keeps the previous "reduce-then-classify"
behavior. The new mode = "disjoin" collapses each set on its own,
then partitions the union into non-overlapping segments with
GenomicRanges::disjoin(), so that every segment is covered by
exactly one combination of sets. Because the resulting intervals are
merged in one mode and disjoint in the other, the reduced_regions
element of GenomicOverlapResult is renamed to regions.
- plotVenn() now attaches eulerr's goodness-of-fit diagnostics
(stress, diagError, regionError, undrawnRegions) to the returned
plot as a "fit_diagnostics" attribute, and the fit itself as
"euler_fit". A message reports diagError and whether it falls above
or below the 1e-6 threshold of Micallef and Rodgers (2014), who
introduced the measure, and names the populated regions the diagram
gives no area to at all, which undrawnRegions records. Set the new
verbose = FALSE to silence it.
- Add plotVennError(), which redraws the diagram shaded by the signed
error of each region, so the regions a diagram misrepresents can be
read off the picture. A thin wrapper around eulerr::error_plot(). It
takes the plot returned by plotVenn() as input.
- Add an ignore.strand argument to
computeOverlaps()/computeGenomicOverlaps(), passed through to
GenomicRanges::reduce(), GenomicRanges::disjoin(), and
IRanges::overlapsAny(). Defaults to FALSE (previous, strand-aware
behavior is unchanged); set to TRUE to disregard strand when merging
or partitioning regions and when determining overlaps.
- computeOverlaps() now warns when two or more input genomic sets
share no chromosome name at all, a common symptom of mismatched
chromosome naming conventions (e.g. "chr1" vs "1") or of comparing
different genome assemblies, cases where overlaps would otherwise be
silently and permanently empty. Genome-assembly conflicts on a
shared chromosome name already error via
GenomicRanges::GRangesList(), unchanged.

Minor updates

- The right-hand annotation of plotUpSet() is now labeled according
to
the type of the input: "Region size" for a GenomicOverlapResult and
"Set size" for a SetOverlapResult.
- computeOverlaps() labels overlap categories faster at large numbers
of regions (~10x at 10^5-10^6 regions), by vectorizing the internal
defineCategories() helper instead of looping row by row.

[See the full NEWS](/packages/release/bioc/news/gVenn/NEWS) or run `news(package = "gVenn")`.


[h5mread](/packages/h5mread)
-------

                        Changes in version 1.6.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Small change in get_h5mread_returned_type() / h5mread().
  get_h5mread_returned_type() now returns "double" instead of "integer"
  if the H5 type of the dataset to read is a 32-bit **unsigned**
  integer
  type (e.g. H5T_NATIVE_UINT32).
  This means that h5mread() now returns an array of type() "double"
  instead of "integer" when reading such dataset. This avoids
  truncation
  of values that are >= 2^31.
  Note that the user can always force the returned type() to "integer"
  by
  calling h5mread() with 'as.integer=TRUE' if they know that the
  dataset
  does not contain values >= 2^31 or if they don't care that big values
  get truncated.

BUG FIXES

- Fix h5mread() on H5T_VARIABLE dataset with NULL pointers.
  Starting with rhdf5 2.57.10, NULL pointers in an H5T_VARIABLE dataset
  are interpreted as NAs. So h5mread() now:
  - works on an H5T_VARIABLE dataset with NULL pointers (used to
    crash!);
  - interprets NULL pointers as NAs like rhdf5::h5read() does.

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.42.0                        

NEW FEATURES

- Add group() methods for HDF5ArraySeed and H5SparseMatrixSeed objects.

- As per rhdf5's recommendation (starting with version 2.57.10),
  writeHDF5Array() and writeTENxMatrix() now use a variable-length
  string dataset instead of a fixed-size string dataset when the
  array- or matrix-like object to write is of type() character.

- Add support for H5SparseMatrix objects on S3 buckets.
  Constructors H5SparseMatrix(), TENxMatrix(), H5ADMatrix(), (and
  corresponding seed constructors H5SparseMatrixSeed(),
  TENxMatrixSeed(),
  and H5ADMatrixSeed()) now accept an H5File object in addition to
  the path to a local HDF5 file. This allows access to an HDF5 file
  stored in an Amazon S3 bucket. See '?H5SparseMatrix', '?TENxMatrix',
  and '?H5ADMatrix' for more information.

- Add 'type' argument to H5SparseMatrixSeed() constructor.
  The new optional 'type' argument lets the user specify the desired
  type()
  of the object to construct. Its default value is NA, in which case
  the
  type() of the returned object will be inferred from the H5 datatype
  of
  the HDF5 dataset.
  All the other matrix-oriented constructors, that is,
  H5SparseMatrix(),
  TENxMatrixSeed(), TENxMatrix(), H5ADMatrixSeed(), and H5ADMatrix(),
  also
  get the new 'type' argument.

- Add 'H5type' argument to writeTENxMatrix().

- Internal helper .load_h5ad_dimnames() can now find rownames and
  colnames
  that are tucked deep in the h5ad file i.e. in '/var/_index/values'
  (for
  the rownames) and in '/obs/_index/values' (for the colnames).

SIGNIFICANT USER-VISIBLE CHANGES

- Significant speedup of the type() getter for HDF5ArraySeed and
  H5SparseMatrixSeed objects. Thanks to Robert Castelo for reporting
  the inefficiency for H5SparseMatrixSeed objects and for proposing
  a fix. See #71.

- Add 'type' argument to getHDF5DumpChunkLength() and
  getHDF5DumpChunkDim().

- writeHDF5Array(x) now propagates chunkdim(x) if it's not NULL.

DEPRECATED AND DEFUNCT

- Deprecate setHDF5DumpChunkLength() in favor of
  setHDF5DumpChunkSize().

BUG FIXES

- Fix bug in coercion from HDF5RealizationSink to HDF5ArraySeed where
  it
  would not propagate type() under certain circumstances.

[See the full NEWS](/packages/release/bioc/news/HDF5Array/NEWS) or run `news(package = "HDF5Array")`.


[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.48.2                        

- new argument 'nthread' in `hlaPredict()` and `predict()` for
  specifying
  the number of threads, as a more intuitive alias of the argument 'cl'

- new argument 'position' in `hlaGDS2Geno()` for importing the SNPs at
  the
  specified base positions only

- fix `hlaGDS2Geno()` with a SeqArray GDS file: "0" is used for the A
  allele at a multi-allelic site (e.g., "0/G"), instead of a list of
  the
  alternative alleles (e.g., "A,C/G"); it works when there is no
  "annotation/id" node in the GDS file

- `hlaGDS2Geno()` uses the variant ID (SeqArray) or the SNP ID
  (SNPRelate)
  when the RS id is not available (i.e., NA, "" or "."), instead of an
  invalid SNP ID

- require SeqArray (>= 1.20.0)

- fix NaN propagation in `hlaPredict()` when a classifier produces zero
  posterior probabilities (e.g., due to large Hamming distances from
  allele/strand mismatch), which caused all predictions to fail

- fix segfault from calling R API (Rprintf, R_CheckUserInterrupt) on
  TBB
  worker threads during parallel prediction

                       Changes in version 1.48.1                        

- fix when "RCPP_PARALLEL_USE_TBB=0"

[HiCcompare](/packages/HiCcompare)
----------

                 Changes in version 1.35.2 (2026-09-30)                 

- Build `GInteractions` input tables from InteractionSet accessors
  instead of
  coercing S4 objects through `data.table`, restoring compatibility
  with
  Bioconductor 3.24.

                 Changes in version 1.35.1 (2026-09-16)                 

- Fix conversion of `GInteractions` inputs, resolving "no method for
  coercing
  this S4 class to a vector" errors.

[HiCDOC](/packages/HiCDOC)
------

                 Changes in version 1.15.4 (2026-09-18)                 

- No longer convert directly GInteractions into data.frame/data.table

                 Changes in version 1.15.1 (2026-05-28)                 

- Imports Data via HiCParser package

- Default filter for small chromosomes is now set to 0

[HiCPotts](/packages/HiCPotts)
--------

                 Changes in version 1.3.1 (2026-09-06)                  

- Parameter documentation and diagnostic messages describe individual
  criteria and thresholds. The fitted model and diagnostic calculations
  are unchanged.

- The function reference now runs a small worked analysis during
  vignette
  builds and includes session information. Argument descriptions are
  formatted for readability in both source and rendered documentation.

- No model changes: the three-component likelihood, Potts formulation,
  posterior target and official classification rule are unchanged.

- Faster computation: cached native likelihood evaluations, direct
  checkerboard allocation updates, shared ABC calibration and optional
  parallel execution of independent chains.

- Regression proposals now retain the tuned QR-whitened random-walk
  kernel
  through burn-in and report ordinary beta acceptance separately from
  label
  swaps, preventing a rejected replacement proposal from freezing a
  chain.

- Robust empirical-Bayes fitting now uses a short multi-chain pilot to
  build
  one shared frozen prior for all production chains, so split-R-hat
  compares
  chains with the same posterior target.

- diagnose_hicpotts_fit() no longer returns an overall verdict. The
  allocation_reliable, parameter_reliable, gamma_reliable, reliable and
  reliability_by_domain fields are removed: one flag collapsed twelve
  independent criteria into a single word, so a marginal result on any
  one
  gate read as a blanket failure. reliability_flags, holding one row
  per
  gate with its threshold, is the reported result.

- summarise_hicpotts_parameters() now judges the estimates it is about
  to
  report instead of requiring every gate to pass. require_reliable =
  TRUE
  rejects a summary only when a precondition leaves the estimates
  ill-defined -- chains on different posterior targets, an unmoved
  regression block, unaligned components, mixed coefficient scales or
  ill-conditioned covariates -- or when a reported parameter falls
  short of
  the requested effective sample size or exceeds the R-hat tolerance.
  Gamma near a boundary, low occupancy and fewer chains than
  recommended no
  longer block the table; they remain visible in reliability_flags. The
  returned table gains a resolved column marking which rows met the
  criteria.

[See the full NEWS](/packages/release/bioc/news/HiCPotts/NEWS) or run `news(package = "HiCPotts")`.


[hoodscanR](/packages/hoodscanR)
---------

                       Changes in version 1.11.1                        

Add the published hoodscanR paper as the preferred package citation.

[HPAanalyze](/packages/HPAanalyze)
----------

                        Changes in version 1.31                         

- Changes in version 1.31.2
  - Update built-in data and internal gene name lookup table to HPA
    version 25.1.

- Changes in version 1.31.1
  - Add hpaSearch() to query the HPA search API.
  - Add hpaXmlParse() to generically parse an entire imported xml
    document
    (based on the proteinatlas.xsd schema) into a flat, joinable list of
    tibbles, and a new vignette documenting it.
  - Add a new vignette documenting hpaSearch(), and cover hpaSearch()
    and
    hpaXmlParse() in the quick-start and in-depth vignettes.
  - Add testthat unit tests and clean up internal code.
  - Various BiocCheck and continuous integration fixes.

[HTAN](/packages/HTAN)
----

- Initial Bioconductor Submission

[iCOBRA](/packages/iCOBRA)
------

                       Changes in version 1.41.1                        

- Replace UpSetR dependency with SimpleUpset

[iconic](/packages/iconic)
------

- Initial Bioconductor Submission

[igblastr](/packages/igblastr)
--------

                        Changes in version 1.4.0                        

NEW FEATURES

- Add get_fwrcdr_widths_for_imgt_organism() to obtain the "real"
  FWR/CDR widths for a given IMGT organism, and use it in
  install_IMGT_germline_db() to compute the internal data that
  gets included in the germline db.
  See '?get_fwrcdr_widths_for_imgt_organism' for more information.

[See the full NEWS](/packages/release/bioc/news/igblastr/NEWS) or run `news(package = "igblastr")`.


[igvR](/packages/igvR)
----

                       Changes in version 1.33.3                        

- Fixed checks examples and styling

                       Changes in version 1.33.2                        

- Fixed BiocCheck ERRORs/WARNINGs: replaced Author/Maintainer fields
  with Authors@R [cre],
  added \value to enableMotifLogoPopups man page, fixed class() ==
  check to use is().

- Updated genomeSpec URL from deprecated S3 bucket to
  https://igv.org/genomes/genomes.json.

- Added MotifDb:: and seqLogo:: namespace qualifiers to resolve 'no
  visible binding' NOTEs.

- Replaced full MIT LICENSE text with DCF stub (YEAR/COPYRIGHT HOLDER).

- Added GitHub Actions CI workflow (R-CMD-check-bioc, Bioc 3.22, R
  4.5).

- Added Podman-based local CI (Containerfile, run_local_ci.R, make
  targets).

- Resolved all lintr violations in R/ and vignettes/ (1752 infix_spaces
  and line_length issues).

                       Changes in version 1.33.1                        

- Added CramTrack.R to Collate field in DESCRIPTION to fix build error.

                       Changes in version 1.33.0                        

- Version bump due to Bioconductor 3.23 devel synchronization.


[ImageArray](/packages/ImageArray)
----------

                        Changes in version 1.1.7                        

New features

- Compatibility with https://github.com/Huber-group-EMBL/romeo
- ImageArray class is now exported.
- axes slot is dropped from ImageArray class and now parsed from
names(scales(object)[[1]]) (which can still be parsed from
axes(object) as usual).
- scales of pyramid axes are now inverted, e.g. 1, 0.5, 0.125 are
now 1, 2, 4.

                        Changes in version 1.1.6                        

New features

- ImageArray function now accepts Bioformats compatible image pyramid
formats beyond ome.tiff and qptiff. See
https://bio-formats.readthedocs.io/en/latest/supported-formats.html
for more information.

                        Changes in version 1.1.5                        

New features

- Built-in create_zarr and create_zarr_group functions are now
replaced by the Rarr::write_zarr_group as of Rarr v2.1.27.

Bug fixes

- Fixing internal read_image method for reading magick bitmap image
class.

                        Changes in version 1.1.4                        

New features

- createImageArray function is now deprecated and replaced by
ImageArray constructor.
- ImageArray function now accepts arguments axes and scales that are
default axes ordering and list of scaling vectors for each pyramid
layers.

                        Changes in version 1.1.3                        

New features

- Novel delayed (or lazy) transformations for image pyramids are
introduced. Available transformations are: Translation, scaling,
rotation, and affine (#45).
- New extent function returns the physical space that the image
occupies.

                        Changes in version 1.1.2                        

New features

- Images with any combination of XYZCT dimensions are can now be
read.
See OME webpage for more information:
https://docs.openmicroscopy.org/ome-model/6.2.2/ome-tiff/specification.html#dimensionorder

[imageFeatureTCGA](/packages/imageFeatureTCGA)
----------------

                        Changes in version 1.2.0                        

New features

- Updated catalog to version 1.1.0 with additional datasets.

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.19.2 (2026-10-05)                 

- fixed malformed Rd markup in testInteractions, readSCEfromTIFF and
  aggregateNeighbors

- moved grid from Suggests to Imports

- require scrapper >= 1.6.0 for aggregateAcrossCells.se

- colData column names from steinbock regionprops now retain their
  original names (e.g. 'centroid-0' instead of 'centroid.0')

- updated as_directed, as_undirected and layout_with_sugiyama functions

                 Changes in version 1.19.1 (2026-05-08)                 

- fixed bug in plotSpotHeatmap function that had also affected
  RELEASE_3_23

[immApex](/packages/immApex)
-------

                        Changes in version 1.7.2                        

BUG FIXES

- Fixed sequenceEncoder() and sequenceDecoder() property mode
silently
mismatching amino-acid identity to property values for
crucianiProperties, MSWHIM, and ProtFP. These property sets' column
order in Peptides::AAdata is not already canonical (MSWHIM/ProtFP
are stored alphabetically; crucianiProperties has an isolated E/Q
swap), and .aa.property.matrix() was returning that raw storage
order instead of reordering to match the canonical amino-acid order
its callers assume positionally. For example, encoding "R" under
MSWHIM silently returned Cysteine's values. .aa.property.matrix()
now always returns canonical-AA-ordered columns, matching
calculateProperty()'s existing contract and the documented
sequenceEncoder()/sequenceDecoder() interface. Verified against
Peptides' own independent public scoring functions (mswhimScores(),
crucianiProperties(), protFP()).
- Fixed sequenceEncoder() and sequenceDecoder() pairing residues with
the wrong property values when property.set was combined with a
non-default sequence.dictionary. Both aligned the property matrix to
the canonical 20 amino acids regardless of the dictionary actually
in use, and neither checked the two agreed. .aa.property.matrix()
now takes the sequence.dictionary and aligns columns to it by name,
which also supports restricted alphabets. Requesting a residue a
scale has no values for is now an error naming that residue, rather
than a silently wrong matrix (this also covers "pK", which exists in
Peptides::AAdata but only defines 9 residues).

DOCUMENTATION

- Corrected sequenceDecoder()'s property.matrix documentation, which
described the matrix as 20 x P. It is P x 20, the transpose of
sequenceEncoder()'s property.matrix argument.

UNDERLYING CHANGES

- amino.acids was defined twice, in R/utils.R and again in
R/calculateProperty.R. Consolidated to a single definition in
R/aaa-constants.R, which collates before both.
- Fixed test-buildNetwork-star.R comparing igraph's arbitrary
component ID labels rather than the component partition itself. The
IDs depend on edge row order, which is not stable between expand =
"clique" and expand = "star", so the test failed intermittently
depending on merge order.

                        Changes in version 1.7.1                        

PERFORMANCE

buildNetwork() now scales to large, clonally expanded repertoires.
The
C++ engine was rewritten around candidate generation instead of
all-pairs comparison, with no new dependencies and identical edge
output
for the default settings:

[See the full NEWS](/packages/release/bioc/news/immApex/NEWS) or run `news(package = "immApex")`.


[immGLIPH](/packages/immGLIPH)
--------

- Initial Bioconductor Submission

[immLynx](/packages/immLynx)
-------

                        Changes in version 1.1.5                        

- Fixed runScXpand() mapping tests that failed on the Bioconductor
builder. They assumed HLA-DRA maps to several Ensembl IDs, which
held in org.Hs.eg.db 3.21 but not in later releases. The tests now
pull a multi-mapping symbol from the installed annotation database
and skip if none exists. Package code is unchanged.

                        Changes in version 1.1.4                        

- Pinned the Linux CI runner to ubuntu-24.04. The ubuntu-latest label
rolled over to 26.04, where the first Python call into immLynxEnv
segfaults inside py_call_impl. The environment pins python 3.9,
numpy 1.23.4 and torch 2.1.2 because tcrdist3, olga and clusTCR
require them, and that stack does not survive the newer toolchain.
The identical package set still passes on macOS and Windows.
Bioconductor's x86_64 Linux builder runs Ubuntu 24.04 LTS, so the
pin also puts CI back on the platform the package is released on.

                        Changes in version 1.1.3                        

- Added CITATION and information for article related to immLynx, see
inst/CITATION

                        Changes in version 1.1.2                        

- Added runScXpand(), which predicts T-cell clonal expansion from
gene
expression alone using scXpand's pretrained pan-cancer models.
Unlike the other wrapped tools it needs no receptor sequences, so it
works on datasets with no paired TCR sequencing. Inference only;
training and hyperparameter optimization stay in Python.

- Added listScXpandModels() to enumerate the available pretrained
models without building the Python environment.

- Added a dedicated scXpandEnv basilisk environment (Python 3.11, CPU
PyTorch, scxpand 0.4.6). It cannot share the existing environments
because scxpand needs Python 3.11 and torch 2.5. The first call
builds several gigabytes and downloads the selected model from
figshare.

- runScXpand() resolves gene identifiers to the Ensembl IDs scXpand's
models are indexed by, optionally mapping symbols through
org.Hs.eg.db. Ambiguous symbols (HLA-DRA alone maps to eight Ensembl
IDs) and collapsed duplicates are counted and reported rather than
resolved silently, because scXpand zero-fills genes it cannot find.

- When scRepertoire clone calls are present, runScXpand() derives
clone_id_size, median_clone_size and expansion per sample using
scXpand's 1.5x-median rule, so a gene-expression-only prediction can
be scored against the observed repertoire. Clone sizes are tabulated
fresh rather than read from clonalFrequency, which
combineExpression() computes under whatever grouping was in effect.

- Pretrained models are cached under tools::R_user_dir("immLynx",
"cache"). scXpand's own default would write a .scxpand_cache
directory into the current working directory, so runScXpand()
downloads the model as an explicit step with an explicit cache
location.

[See the full NEWS](/packages/release/bioc/news/immLynx/NEWS) or run `news(package = "immLynx")`.


[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.45.1                        

- CHANGES
  - added ICL_bias accessor for immunoClust- and immunoMeta-objects

[IntegratedLearner](/packages/IntegratedLearner)
-----------------

- Initial Bioconductor Submission

[InteractiveComplexHeatmap](/packages/InteractiveComplexHeatmap)
-------------------------

                       Changes in version 1.21.1                        

- `add_js_css_dep()`: directly paste JS code by `HTML()` instead of
  specifying `htmlDependency()`.

[IntrinsicScan](/packages/IntrinsicScan)
-------------

- Initial Bioconductor Submission

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.48.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Minor tweaks to as.data.frame.IPosRanges() and as.data.frame.IPos()

- Revisit show() method for MaskCollection objects.

BUG FIXES

- Fix handling of names and mcols in ranges() method for Views objects.

- Various small tweaks and fixes by Michael Chirico
  <chiricom@google.com>.
  See #65, #66, #61, and #67.

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.23.2 (2026-09-24)                 

UPDATE

- After fold-abundance filtering, remove_collisions() removes an
integration site from each independent sample represented by only
one row and retains it in every sample represented by at least two
rows. The same integration site may therefore remain in multiple
patients.
- Sequencing date is no longer used to resolve collisions. date_col
is
kept as an ignored argument for backward compatibility, and
ReplicateNumber is no longer required.

                 Changes in version 1.23.1 (2026-07-16)                 

UPDATE

- remove_collisions() now resolves collisions by first comparing
relative integration-site abundance across independent samples.
Abundance is measured with the sequence count column selected
through quant_cols["seqCount"], summed within each independent
sample.
- Observations at least fold_threshold times less abundant than the
maximum abundance observed for the same integration site are removed
before any temporal rule is applied. The default threshold is 10 and
the exact threshold is inclusive (max(seqCount) / seqCount >=
fold_threshold).
- The temporal rule is now applied only to collision observations not
resolved by the fold-abundance comparison.

[iscream](/packages/iscream)
-------

                        Changes in version 1.3.2                        

- Update package logo

                        Changes in version 1.3.1                        

DOCUMENTATION

- Note that there may be duplicate records returned when input
regions
overlap

[iSEEid](/packages/iSEEid)
------

- Initial Bioconductor Submission

[IsoNetAnalyzeR](/packages/IsoNetAnalyzeR)
--------------

- Initial Bioconductor Submission

[jvecfor](/packages/jvecfor)
-------

                        Changes in version 4.0.1                        

- Release of Bioconductor Version 3.24
- Update jfector backend to v4.0.1
- Examples, vignette and tests now skip when Java < 20 is on PATH,
instead of failing (fixes build ERROR on hosts with an older JVM)
- Look for Java in JAVA_HOME as well as on PATH, and use whichever is
Java >= 20


[karioCaS](/packages/karioCaS)
--------

- Initial Bioconductor Submission

[kebabs](/packages/kebabs)
------

                       Changes in version 1.47.1                        

- fixed bug in model selection that led to an error with new version of
  LiblineaR

- various improvements of package vignette

[KEGGemUP](/packages/KEGGemUP)
--------

- Initial Bioconductor Submission

[KEGGREST](/packages/KEGGREST)
--------

                       Changes in version 1.53.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- The "organism" KEGG endpoint has stopped working due to upstream
  changes,
  the "genome" endpoint can be used instead with `keggList("genome")`
  (@assaron, #30)

BUG FIXES AND MINOR IMPROVEMENTS

- Updated `color.pathway.by.objects` to align with the new KEGG Map
  Coloring
  web application specifications, using space-separated multi-query
  format.

- Added unit tests and examples for querying the KEGG Orthology
  database
  e.g., `keggLink("ko", "sce")` (@assaron, #38)

[knowYourCG](/packages/knowYourCG)
----------

                        Changes in version 1.9.7                        

- library(knowYourCG) prints a one-line description and how to cite the
  paper (suppressPackageStartupMessages() silences it).

- Array and Sequencing vignettes point to the kycg command-line tool
  and
  the web app (try it online, no install).

                        Changes in version 1.9.6                        

- BEHAVIOUR CHANGE: testEnrichmentSEA reports nQ as the size of the
  query
  and nD as the size of the database. For a numeric query the two were
  the
  other way round (nQ was the set, nD the CpGs with a value).

- testEnrichment warns when no query probe is in the universe, instead
  of
  returning all-zero overlaps silently.

- Databases dropped while loading are named in a warning: a requested
  group
  with no version for the query's platform, an unknown group, an empty
  set.

- Vignettes: no sesameDataCache("HM450.address"); sesameData 1.31.2
  lists
  that record with no ExperimentHub id yet, which stopped the devel
  build.

                        Changes in version 1.9.5                        

- testEnrichment: no "Unknown or uninitialised column" warnings on
  sequencing results; the estimate's 1024 / -1022 sentinels are
  documented.

- Vignettes: an install chunk for the Suggests the examples load; the
  Manhattan and feature-aggregation chunks cache their own data; MM285
  examples use the MM285 TFBS database; the HM450 Zenodo example uses
  an
  HM450 query; the sequencing quick start runs on the shipped single
  cell;
  the bedToCg example is self-contained; permutation chunks are seeded.

- Manual: loadDBs, testEnrichmentSpearman and df_master are internal;
  no
  developer notes or local paths in examples.

                        Changes in version 1.9.4                        

- bedToCg: the bedtools/yame pipeline reached sh -c unquoted, so sh ran
  a
  bare `yame` and the output file was empty. It is now quoted.

- testEnrichment2: the odds ratio is computed in double precision.
  Integer
  counts from YAME overflowed on genome-scale tables and gave NA
  estimates.

- testEnrichment2: an unmeetable min_overlap returns the empty result
  without a second warning.

- KYCG_plotEnrichAll: geom_segment(linewidth=), not the deprecated
  size=.

[See the full NEWS](/packages/release/bioc/news/knowYourCG/NEWS) or run `news(package = "knowYourCG")`.


[LACHESIS](/packages/LACHESIS)
--------

                        Changes in version 1.1.1                        

NEW FEATURES

- LACHESIS now allows de novo estimation of mutation rates from
  cohort data
- Inferred time axis updated
- Copy number segments with negative values are now filtered

[lantern](/packages/lantern)
-------

- Initial Bioconductor Submission

[lcmsPlot](/packages/lcmsPlot)
--------

                       Changes in version 1.1.29                        

- The Compound Discoverer node's compound table has a logical
has_name
column, TRUE for compounds Compound Discoverer has named. Unnamed
compounds are labelled with their formula or ID, so a Compounds
Query of has_name (or has_name & compound_rank <= 20) plots only the
named ones. Released as node 1.1.0, which needs lcmsPlot 1.1.29 or
later.

                       Changes in version 1.1.28                        

- Compound Discoverer node released with an installer: each version
of
the Compounds Chromatograms node (its own version, starting at 1.0.0
in VERSION) is published as a GitHub release,
CompoundsChromatograms.zip, with a Windows installer (install.bat).
The installer finds R 4.4.0 or later, installs lcmsPlot at the
version the node was released with (from its GitHub tag, as
Bioconductor's release may be older) and rawrr, copies the node into
each Compound Discoverer version found and writes the R and script
paths into node.json. Package versions are tagged vX.Y.Z, and node
releases cd-node-vX.Y.Z, by the release workflow.
- The node checks, when run from a release, that the installed
lcmsPlot is at least the version it was released with, and stops
with a clear message if not. It now ships the icons its node.json
refers to, and its default ExecutablePath points to an R that
lcmsPlot supports (4.4).

                       Changes in version 1.1.27                        

- Peak integration windows in the Compound Discoverer node (#41): the
Compounds Chromatograms scripting node now requests the Features per
File and Chromatogram Peaks tables, and
CompoundDiscovererNodeSource() takes each compound's peak window in
each file from them: the left and right retention times of the
reference chromatogram peak of the compound's reference ion, the
window Compound Discoverer integrated the area over, as the
.cdResult reader does. Previously the node exported no peak bounds,
so the window was the apex ± half the FWHM (about a third of the
integrated window) or the apex alone. The reference peak is the one
flagged Is Ref Peak, else isotope 0, else the one whose apex is
closest to the ion's retention time. Exports without these tables
keep the previous behaviour, and per-file feature exports without
Left/Right RT now also fall back to the FWHM. Workflows built with
the previous node.json need the node re-registered (or its Requested
Tables and Columns updated) to benefit.

                       Changes in version 1.1.26                        

- Compound Discoverer 3.5 result files: lp_compound_discoverer()
reads
.cdResult files written by Compound Discoverer 3.5, which name the
table linking compounds to their per-file instances
UnknownCompoundInstanceItemsConsolidatedUnknownCompoundItems instead
of ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems.
Previously such files failed with "no such table".

[See the full NEWS](/packages/release/bioc/news/lcmsPlot/NEWS) or run `news(package = "lcmsPlot")`.


[LEA](/packages/LEA)
---

                       Changes in version 3.25.2                        

- Correct genomic inflation factors in lfmm2.test() for multivariate
  GLM
  association tests: compute one factor per environmental covariate
  across loci,
  rather than one per locus.

- Calibrated p-values may change substantially; rerun affected
  analyses.

- Raw GLM statistics, uncalibrated p-values, single-covariate behavior,
  linear tests and full tests are unchanged.

[lefser](/packages/lefser)
------

                       Changes in version 1.26.0                        

Minor improvements and bug fixes

- Updated the vignette to use mia::convertFromPhyloseq() in place of
the deprecated mia::makeTreeSummarizedExperimentFromPhyloseq()
(@0xMuluh, #105)

                       Changes in version 1.24.0                        

Minor improvements and bug fixes

- Added CITATION.cff file

[lemur](/packages/lemur)
-----

                       Changes in version 1.11.2                        

- Remove the dependency on the harmony package. The maximum diversity
clustering step that align_harmony relies on is now implemented
natively in pure R (see
init_max_diversity_clustering/run_max_diversity_clustering),
avoiding breakage from upstream harmony changes. align_harmony's
behavior and arguments are unchanged.

[levi](/packages/levi)
----

                        Changes in version 2.0.0                        

Statistical corrections

- leviSingleCellInteractionTFCE() now tests per-donor interaction
  contrasts, (stim - ctrl) in a cell type minus (stim - ctrl) in the
  reference type, with a limma-trend one-sample moderated t and a
  within-donor sign-flip null. The previous refit of
  ~ donor + type * condition on swapped labels was valid only under the
  global null: with a shared condition effect and no interaction it
  rejected in 85% of null data sets (now 3%, attainable level 3.1%).
  Cell types are aligned on their shared genes before the contrasts are
  formed; stacking matrices filtered per type had mixed genes across
  rows or aborted.

- leviBulkGraphInference() refits the whole DESeq2 pipeline in every
  permutation by default (refit_dispersions = TRUE). Keeping the
  label-dependent size factors and dispersions fixed gave a
  family-wise error rate about three times the attainable level; that
  mode now warns that it has no level guarantee. Cook-outlier
  replacement is switched off so that observed and permuted fits use
  the
  same pipeline.

- leviGraphMoran() computes local P-values by conditional permutation
  (each node keeps its own score; its neighbours are drawn from the
  other scores), as in spdep and GeoDa, and centres the two-sided
  P-values on the null expectations. Total permutation called an
  isolated outlier with iid neighbours significant in about half of
  null
  data sets.

- Documentation of the node-label and rewiring tests (levi(),
  leviGraphGetisOrd(), leviGraphSpectrum(),
  leviGraphWeightedTopology(),
  leviGraphRewiringInference(), leviGraphTFCERewiring()) states that
  their nulls treat scores as exchangeable over nodes and reject under
  co-expression or degree association without any condition effect.

Behaviour changes (these alter landscapes and regional P-values

- Nodes without a measurement have weight zero in the landscape, and so
  do the edge midpoints that touch them. They used to enter the average
  at the neutral value with full weight, which pulled the surface to
  "no
  change" around every unmeasured gene. Unmeasured parts of the network
  now drop out of the average and out of the silhouette; the weights
  are
  the same in every permutation, and the result stores them in
  $metadata$support_weights. The effect grows with the share of missing
  nodes; a complete expression table gives the same result as before.

- Reciprocal edges (A-B and B-A) and self-loops count once in the
  landscape, as they already did in the graph tests. A network file
  that
  lists both directions no longer puts two support points on each edge
  midpoint.

[See the full NEWS](/packages/release/bioc/news/levi/NEWS) or run `news(package = "levi")`.


[limma](/packages/limma)
-----

                 Changes in version 4.0.0 (2026-10-30)                  

- voomLmFit() transferred to limma from the edgeR package.
  Previously, voomLmFit() needed to call edgeR::glmFit() in order
  to identify exact zero fitted values, but voomLmFit() now calls
  a new limma C routine 'poisfit' for the same purpose, removing
  the dependence on edgeR.

- For the purpose of estimating sample weights, voomLmFit() calls
  arrayWeights() with 'method="reml"' if there is no loss of df
  due to exact zero fitted values and 'method="genebygene"'
  otherwise.

- The voomLmFit help page has been expanded, while the voom and
  voomWithQualityWeights help pages now have preambles that point
  users to voomLmFit.

- voomLmFit() now gives priority to 'offset' over 'offset.prior'
  or 'lib.size' and no longer forces the rows of 'offset.prior'
  to add to zero.

- New argument 'contrasts' for voomLmFit(), lmFit(), lm.series()
  and gls.series(), which converts coefficient estimates and
  standard deviations into contrast estimates and standard
  deviations. This argument provides the same functionality as
  calling 'contrasts.fit' on the output object, but with exact
  rather than approximate standard deviations.

- New argument 'nthreads' for voomLmFit(), lmFit(), lm.series(),
  gls.series(), duplicateCorrelation(), and arrayWeights() to
  enable OpenMP parallelization of genewise fits.

- C code backends for lm.series() and gls.series().

- C code backend for duplicateCorrelation().
  duplicateCorrelation() previously called statmod::mixedModel2()
  to estimate genewise correlations, but this functionality is
  now absorbed into the C code. The new code is able to evaluate
  the genewise SVDs needed by mixedModel2() more efficiently than
  was possible for the univariate R function, resulting in
  substantial speed and memory improvements.

- C code backend for .arrayWeightsPrWtsREML(), which is called by
  arrayWeights() when 'method="reml"' and prior weights are set.

- Minor changes to weightedLowess() C code to fit in with the new
  C code elsewhere.

- New function normalizeBetweenSamples(). Similar to
  normalizeBetweenArrays() but works on ELists.

- camera() and cameraPR() now support gene sets as expression
  signatures, whereby a logFC or similar expected directional
  change is provided for each gene in the gene set.

- avearrays() methods for EList and MAList objects now sum the
  weights instead of averaging them. Other matrices are still
  averaged.

- Fix filtering error in Agilent rat case study in User's Guide.

[See the full NEWS](/packages/release/bioc/news/limma/NEWS) or run `news(package = "limma")`.


[LIPIDIFy](/packages/LIPIDIFy)
--------

- Initial Bioconductor Submission 
 
[LipidTrend](/packages/LipidTrend)
----------

                 Changes in version 1.3.1 (2026-08-24)                  

- Breaking change: removed own_contri from analyzeLipidRegion(); self-
  contribution is now fixed at 0.5.

- Breaking change: analyzeLipidRegion() now requires exactly 2 groups.

- Added regionalTestFC(): region-level paired t-test and fold-change
  for
  the Increase/Decrease regions.

- plotRegion1D()/plotRegion2D() now use regionalTestFC() to decide
  coloring
  when both directions have >5 significant features; new
  regional_p_cutoff
  argument (default 0.05) controls this.

                        Changes in version 1.3.0                        

- Removed the dependency on the MKmisc package.

- Updated RoxygenNote from 7.3.2 to 7.3.3.

[lisaClust](/packages/lisaClust)
---------

                       Changes in version 1.21.13                       

- lisaClust() put regions on the right cells only when cell IDs were
  unique across images; with cell IDs
  that repeat in each image (1, 2, ... per image), most cells were left
  without a region. Regions now go
  back to the cells by position, and results with unique cell IDs are
  unchanged.

- lisa() returns its rows in the order of `cells` rather than grouped
  by image, so its curves can be
  clustered and put back on a SingleCellExperiment or SpatialExperiment
  directly. Input that is already
  ordered by image gives the same matrix as before. Without a cellID
  column, the rows are named cell_1,
  cell_2, ... in the order of `cells`.

- lisa() and lisaClust() stop with a message naming the column when
  image IDs, cell types or coordinates
  are missing (a missing cell type crashed R; a missing image ID
  dropped the cell). Image IDs that are
  factors with unused levels work. lisaClust() checks the length of
  `cellWeights`.

- hatchingPlot() treats a numeric region column as categories instead
  of failing with "Continuous value
  supplied to a discrete scale".

- Deprecated arguments of lisa() are reported as lisa()'s rather than
  lisaClust()'s.

                       Changes in version 1.21.12                       

- Matrix is listed in Suggests (it is used only when `cellWeights` is a
  sparse Matrix), and inhomLocalK()
  documents its internal `membership` argument, clearing the two R CMD
  check warnings.

                       Changes in version 1.21.11                       

- regionBoxPlot() takes `subject`: the shares are then those of each
  patient's cells, all images pooled, one
  point per patient.

- regionMap() turns the region labels upright when there are more than
  six regions, so they do not overlap.

- The vignette is simpler: the example tumours are named rather than
  computed, the group medians come from
  spicyR's colTest() (which now reports them), patients are given as
  `subject`, and the region labels need no
  theme() adjustment. Needs spicyR 1.99.38.

                       Changes in version 1.21.10                       

- Tests pass with a single-threaded BLAS (R's reference BLAS, as on the
  Bioconductor build machines).

[See the full NEWS](/packages/release/bioc/news/lisaClust/NEWS) or run `news(package = "lisaClust")`.


[LRDE](/packages/LRDE)
----

                        Changes in version 1.1.2                        

- Improved numerical stability of Hurdle-NB parameter estimation.

- Added safeguards against non-finite likelihood values during
  BFGS optimization.

- Improved handling of extreme negative binomial parameters and
  zero-truncation probabilities.

- Added protection against zero or non-finite initial parameter
  values.

- Improved initialization of group-specific mean parameters in
  tag-wise dispersion estimation.

- Added probability clipping to prevent infinite logit
  transformations.

- Improved robustness of tag-wise dispersion estimation for sparse
  genes and genes with all-zero counts.

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

[MACSr](/packages/MACSr)
-----

                 Changes in version 1.21.4 (2026-10-21)                 

- Upgrade to MACS 3.1.0, which requires Python 3.12.

- hmmratac: add remove_dup and em_jump; deprecate keepduplicates.

- hmmratac: decoding_steps now defaults to 5000, as in MACS 3.1.0.

- Add threads argument (default 1) to callpeak, predictd, pileup,
  filterdup, randsample, refinepeak, hmmratac and callvar.

- callvar: deprecate np; use threads instead.

- pileup: fix the default format = "AUTO", which MACS3 pileup rejected.

- Add CITATION for the MACS3 paper.

                 Changes in version 1.21.3 (2026-10-07)                 

- Use bundled test data in examples and tests instead of ExperimentHub.

                 Changes in version 1.21.2 (2026-10-02)                 

- Use bundled test data in the vignette instead of ExperimentHub.

[MagellanNTK](/packages/MagellanNTK)
-----------

- Initial Bioconductor Submission

[mariner](/packages/mariner)
-------

                       Changes in version 1.13.1                        

New features:

Support for reading .cool and .mcool (cooler) files, contributed by
Sarah Parker.

- pullHicPixels() and pullHicMatrices() now accept .cool and .mcool
files in addition to .hic. File type is detected automatically, so
no change to existing code is required. Functions built on top of
these, such as pileupPixels(), pileupDomains(), pileupBoundaries(),
and calcLoopEnrichment(), therefore accept cooler files as well.

- coolStraw() extracts a sparse contact matrix from a cooler file,
mirroring the interface of strawr::straw().

- readCoolChroms(), readCoolBpResolutions(), and readCoolNormTypes()
report the chromosomes, resolutions, and normalizations available in
a cooler file.

Limitations:

- Only matrix = "observed" is supported for cooler files. Expected
and
observed/expected matrices remain .hic only.

- A single call cannot mix .hic and cooler files.

[markeR](/packages/markeR)
------

                        Changes in version 1.3.1                        

Minor Changes

- Updated the citation to the markeR paper published in NAR Genomics
and Bioinformatics.
- Updated the contact email address.
- Fixed runGSEA function.

[MassSpecWavelet](/packages/MassSpecWavelet)
---------------

                 Changes in version 1.79.2 (2026-07-05)                 

- tuneInPeakInfo(): Fix several crashes:

- A peak whose local refinement window contains no usable ridges
(e.g. a quiet/flat neighborhood, or a window near the edge of
the spectrum with too few scales to search) no longer crashes
with Error in strsplit(ridgeName, "_") : non-character argument;
the peak is now reported as unprocessed instead.
- Calling tuneInPeakInfo() with peakIndex/peakScale directly
(instead of majorPeakInfo) no longer crashes with Error in
names<-.

- Deprecated the internal, unexported functions getRidgeValue(),
i2u(), u2i(), mzInd2vRange(), mzV2indRange() and smoothDWT(). These
functions appear to be unused, both internally in MassSpecWavelet
and by known downstream packages, and are scheduled for removal in a
future release. If you rely on any of them (e.g. via :::), please
open an issue at https://github.com/zeehio/MassSpecWavelet/issues.

- Migrated the unit test suite from RUnit to testthat (3rd edition).
Tests now live in tests/testthat/ instead of inst/tests/, following
the standard testthat package layout. Significant increased in test
coverage.

                 Changes in version 1.79.1 (2026-07-05)                 

- identifyMajorPeaks(): Fix Error in strsplit(ridgeName, "_") :
non-character argument crash when no ridges are found (e.g. for a
flat/constant input signal). An empty (but well-formed) result is
now returned instead.

- getRidge(): Fix ridge name/lookup mismatches on high resolution
spectra (m/z indices at or above 1e5). Ridge indices were internally
promoted from integer to double partway through the tracing loop,
and R's default number-to-string coercion uses scientific notation
for large "round" doubles (e.g. 200000 -> "2e+05") but never for
integers. The same peak could therefore end up keyed under two
different strings, causing the internal ridge lookup to silently
miss and ridges to be dropped or duplicated. Ridge names are now
always built with a representation- independent formatting helper.
Thanks to drobertsicl for reporting the issue and suggesting the
fix. Closes #8

[mastR](/packages/mastR)
-----

                       Changes in version 1.13.1                        

- Update vignette name to avoid repeated names.

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.21.1 (2026-08-04)                 

- fix error by adding namespaces to shiny::testServer

[matter](/packages/matter)
------

                 Changes in version 2.15.1 (2026-11-16)                 

BUG FIXES

- Fix grouped var/sd in rowStats() when chunks are missing groups

[MDSvis](/packages/MDSvis)
------

                         Changes in version 1.1                         

MDSvis 1.1.1

- max file upload size is now a parameter in mdsvis_app()

[MeLSI](/packages/MeLSI)
-----

                        Changes in version 1.1.7                        

- Performance: the weak-learner metric optimizer
(optimize_weak_learner_robust, the hot path of ensemble metric
learning) is now a C++ kernel (melsi_opt_weak_learner in
src/melsi_fstat.cpp, via Rcpp). Profiling showed this
gradient-descent loop accounted for ~83% of total runtime; the C++
port runs it ~4.6x faster, for roughly ~2x end-to-end speedup,
composing with the 1.1.6 F-statistic kernel and BPPARAM parallelism.
- Reproducibility: results are bit-faithful to the R implementation.
The C++ loop draws its within/between-class pairs with R's own index
sampler (R_unif_index), so it consumes the global RNG stream
identically to the previous sample()-based loop. End-to-end
F-statistics reproduce 1.1.6 to floating-point precision (~1e-15
relative) with no p-value changes.

                        Changes in version 1.1.6                        

- Performance: the PERMANOVA pseudo-F statistic is now computed by a
fused C++ kernel (src/melsi_fstat.cpp, via Rcpp) that accumulates
squared distances in a single pass over sample pairs, instead of
materialising and squaring the full n x n distance matrix in R. This
speeds up the F-statistic by ~4x and the overall analysis by
roughly 1.3-2.2x (most for two-group / pairwise analyses), and
composes multiplicatively with BPPARAM parallelism.
- Numerical note: results are statistically equivalent but not
bit-identical to 1.1.5. Computing squared distances directly (no
sqrt round-trip) and summing in a different order changes
F-statistics at the ~1e-12 relative level. Across extensive testing
no permutation p-value changed; in rare boundary cases a p-value
could differ by a single permutation step (1/(n_perms+1)). Set a
larger n_perms if exact reproducibility of borderline p-values
matters.

                        Changes in version 1.1.5                        

- Dependencies: vegan moved from Imports to Suggests. The package
computes the PERMANOVA F-statistic directly and no longer calls
vegan::adonis2() at runtime, so vegan (and its transitive
dependencies) is no longer installed for users; it is needed only to
build the comparison section of the vignette.
- Cleanup: removed unused stats::aov, stats::as.dist, and
stats::t.test imports and an unreachable t-test pre-filtering branch
in the metric learner. Pre-filtering is handled by
apply_conservative_prefiltering() as before.
- UX: permutation progress is now shown with a single progress bar
instead of one message per permutation.
- UX: melsi() gains an optional seed argument for reproducible runs
and now returns an object of class "melsi" with a print() method
summarising the F-statistic, p-value, and top features. The returned
object is still a list, so existing $ access is unchanged.
- All statistical results (F-statistics, p-values, feature weights)
are numerically identical to 1.1.4; these changes do not alter the
method.

[See the full NEWS](/packages/release/bioc/news/MeLSI/NEWS) or run `news(package = "MeLSI")`.


[memes](/packages/memes)
-----

                       Changes in version 1.21.1                        

- Fixes a bug converting string sets to data.frame in newer versions
of bioc

[meshes](/packages/meshes)
------

                       Changes in version 1.39.2                        

- gseMeSH() now exposes an eps argument and forwards it to
enrichit::gsea_gson(), matching the other GSEA wrappers in the
family (2026-09-22, Tue)
- expose an explicit seed argument on gseMeSH() and forward it to
enrichit::gsea_gson() (2026-08-15, Fri)
- set it to a number (or TRUE for a fixed default seed) to get
identical GSEA results across runs, FALSE (default) draws a
fresh seed on each run; set.seed() before the call still works
- pvalueCutoff is applied to both pvalue and p.adjust inside
enrichit::gsea_gson() (requires enrichit >= 0.2.2), matching the
historical clusterProfiler/DOSE double-filtering behavior

                       Changes in version 1.39.1                        

- use Quarto for the package vignette and keep examples in the online
book (2026-06-30, Tue)
- update package description to better summarize supported data and
analyses (2026-06-30, Tue)
- improve performance of gene-to-MeSH lookup and MeSH
ancestor/offspring traversal (2026-06-30, Tue)
- validate MeSHDb resources before querying to report clearer errors
for incompatible AnnotationHub results (2026-06-30, Tue)

[metabinR](/packages/metabinR)
--------

                        Changes in version 2.2.0                        

New features

- bin_summary() reports the read count, proportion, and
assigned-distance summaries for each observed bin in a
MetabinResult.
- ambiguous_reads() identifies reads whose two closest finite bin
distances differ by no more than a chosen margin. For hierarchical
results, distances outside the read's parent abundance bin are
ignored.
- evaluate_bins() joins bin assignments to known read origins by ID
and reports a confusion table, per-bin purity, per-origin best-bin
recovery, and the adjusted Rand index. All measures use read counts.

[MetaboAnnotatoR](/packages/MetaboAnnotatoR)
---------------

                        Changes in version 1.1.4                        

ADDITIONS

- New function to re-annotate AIF and ISF pseudo-MS/MS spectra.

                        Changes in version 1.1.3                        

ADDITIONS

- New function to annotate spectra loaded as MSnExp objects.

                        Changes in version 1.1.1                        

ADDITIONS

- LC-MS in-source fragmentation annotation functionality added.

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                        Changes in version 1.21                         

Changes in 1.21.2

- translate() gains parameter nomatch to configure behavior when
names
can not be translated.

Changes in 1.21.1

- Add adductCharge() and standardizeSingleCharge() functions.

[MetaboDynamics](/packages/MetaboDynamics)
--------------

                       Changes in version 2.3.15                        

- in "fit_dynamics_model" argument "raw_plus_counts" changed to
"raw_plus_cell_counts"

                        Changes in version 2.3.9                        

- in "plot_estimates" argument "dynamics" changed to "MDVs"
- in "plot_ORA" OvE instead of p(OvE) are displayed, reducing
distortion due to log-transformation

                        Changes in version 2.3.8                        

Vignette polishing

                        Changes in version 2.3.7                        

- bug fix in prior settings including corrected explanation in
vignette

                        Changes in version 2.3.1                        

new functionality: priors in model "raw_plus_counts" are adjustable
by
user, compatible with older versions as default priors are previous
fixed priors

- if run on Summarized Experiment object, priors are stored in
metadata(data)[["priors"]]
- new vignette "MetaboDynamics_Guide_to_prior_setting"

[MetaPathNet](/packages/MetaPathNet)
-----------

- Initial BIoconductor Submission

[methylclock](/packages/methylclock)
-----------

                        Changes in version 2.0.0                        

(In Bioconductor devel as the 1.99.x series, following the convention
for a major rewrite; it becomes 2.0.0 at the next release.)

A ground-up rewrite. The package now computes 44 epigenetic clocks
behind a declarative registry, runs out-of-core on HDF5 through
BigDataStatMeth, and ships a toolbox to compare, interpret and
quality-control the results.

Architecture

- Declarative clock registry (clock_register(), clock_list(),
clock_info(), and the shipped clock_catalog table): one entry per
clock — target, engine, what it estimates, native platforms,
training tissue, licence and citation. Adding a clock is adding a
row, not editing the engine.
- Manifest-backed resource resolver (mcd_resource(), mcd_backends()):
coefficients are resolved through a configurable backend chain — a
local mirror, ExperimentHub, or a direct Zenodo download — and
normalized to one canonical shape on entry (normalize_coef()).
- Three predictor engines in this release: linear (including models
with squared-beta terms), counter (mitotic clocks, including the
parametric epiTOC2), and neural networks (AltumAge with weights from
HDF5; a compiled Bayesian network).
- A unified methylclock() entry point; DNAmAge() and DNAmGA() remain
as compatible wrappers.

Clocks

- 44 clocks: chronological (including the cross-platform Garma 2024
model, the only one here native to EPICv2), biological (PhenoAge),
pace of ageing (DunedinPACE), telomere length, gestational (cord
blood and placenta), neonatal (the four NEOage clocks), mitotic
counters, causal clocks, and nine trait/exposure EpiScores (opt-in
via target = "trait").
- Every clock declares the tissue it was trained on, the quantity it
estimates with its units (the four NEOage clocks report
post-menstrual and postnatal age in weeks; Levine, phenotypic age;
TL, telomere length in kilobases), and its native arrays. The arrays
are where a clock was built and validated, not where it can be
computed — CpG coverage decides that. Restrictively licensed clocks
are implemented and clearly marked.

Scale

- Input can be an in-memory matrix or an HDF5-backed matrix;
whole-array clocks (BLUP, EN) stream from disk without materializing
the array, and the heavy algebra runs through BigDataStatMeth.
mc_to_hdf5() converts a matrix or a text file once, to work on disk
from then on.
- Results follow the same policy: HDF5 as the durable store with an
in-memory cache, and persist() to keep a session result.

[See the full NEWS](/packages/release/bioc/news/methylclock/NEWS) or run `news(package = "methylclock")`.


[methylSig](/packages/methylSig)
---------

                       Changes in version 1.26.0                        

NEW FEATURES

- diff_methylsig() has local_disp and local_meth parameters, which
  choose
  whether local information (local_window_size > 0) is used for the
  dispersion estimate and degrees of freedom, the group methylation
  estimates and likelihood ratio, or both (issue #47). Both are TRUE by
  default, which gives the same results as before. These were the
  local.disp and local.meth options of methylSigCalc() before 0.5.0.

- diff_dss_test() has a covariate_percentiles parameter, which sets the
  percentiles used to group samples by a numerical
  methylation_group_column
  (issue #44). The default, c(25, 75), is the same grouping as before.

- diff_dss_test() checks the contrast before testing. A contrast of the
  wrong size gives an error with the number of values it needs and the
  columns of diff_fit$X they correspond to, in order, e.g.
  "1: (Intercept)" and "2: Typenormal". A message says what the
  contrast
  tests, e.g. "Testing Typenormal = 0", so a contrast that is the right
  size but tests the wrong coefficient is easier to catch.

- The contrast can be a vector as well as a matrix, and named values
  (or
  row names) are matched to the columns of diff_fit$X by name.
  Contrasts
  with NA values, a column of zeros, or linearly dependent columns give
  an error.

SIGNIFICANT USER-VISIBLE CHANGES

- diff_dss_test() with a numeric methylation_group_column reports
  methylation rates the other way around, so they agree with the
  covariate. The case group (meth_case) is now the samples at or above
  the
  upper covariate percentile, and the control group (meth_control) the
  samples at or below the lower one. meth_diff is high minus low, and
  direction is Hyper when methylation increases with the covariate,
  which
  is the sign of stat for a contrast on the covariate. Before, the low
  group was the case group, so a locus gaining methylation with the
  covariate was labeled Hypo with a negative meth_diff. Results with
  character or factor columns, and the test statistics and p-values,
  are
  unchanged.

- The example data sets (bsseq_stranded, bsseq_destranded,
  bsseq_multichrom, promoters_gr) are no longer lazy-loaded, as
  Bioconductor recommends. Load them with data(), e.g.
  data(promoters_gr, package = 'methylSig'), as the examples and
  vignettes
  already do.

[See the full NEWS](/packages/release/bioc/news/methylSig/NEWS) or run `news(package = "methylSig")`.


[methylumi](/packages/methylumi)
---------

                       Changes in version 2.59.3                        

Bug fixes

- methylumIDAT() reads IDAT files from per-slide subdirectories,
which
is the layout Illumina's own software produces
(<idatPath>/<Slide>/<Slide>_<Array>_Grn.idat). Both the
file-presence check and the code that actually opened the files
assumed a flat directory, so a nested run reported every sample as
missing. Flat layouts are unaffected. If the same basename appears
under more than one subdirectory the read now errors and lists the
candidates rather than picking one. getBarcodes() had the same blind
spot and is fixed too. (#19)

- diagnostics() works on objects with no color-channel annotation,
such as the bundled GoldenGate mldat. Three if (annotation(x) ==
"...") tests errored outright on a zero-length annotation, and past
those, the plotting loop split probes on a COLOR_CHANNEL column that
such objects do not have, so every panel would have been drawn from
zero probes. The long-reported call to the non-existent
plot.density() was real, on the same path, and is also fixed. (#23)

- methylumi.bgcorr() explains itself when asked for noob correction
on
an object that has no out-of-band intensities. Objects built from
GEO series matrices or GenomeStudio output via methylumiR() carry
only methylated/unmethylated/p-value data; OOB probes come only from
reading raw IDATs. This used to fail deep inside method dispatch
with "unable to find an inherited method for 'intensities.OOB'".
(#24)

Deprecated and defunct

- The gamma-family background corrections — method = "goob", "gamma"
and "mode" — have been removed. They called gamma.mle(),
gamma.mode() and gamma.integral() from rGammaGamma, a GitHub-only
package that is not available from CRAN or Bioconductor and was
never declared as a dependency, so these methods have errored for
years. Requesting one now produces an error saying so. method =
"noob" (the default) is unaffected. (#18)

Internal changes

- qc.probe.plot(), plotNAs() and plotProbeNAs() no longer use
ggplot2's deprecated qplot(). Plot output is unchanged except in
plotNAs(), where the segment layer now maps linewidth rather than
the also-deprecated size; scale_linewidth uses a linear palette
where scale_size used an area one, so line widths differ slightly.
(#40)

- A pkgdown site is published to
https://seandavi.github.io/methylumi/
from a GitHub Actions workflow. (#43)

[See the full NEWS](/packages/release/bioc/news/methylumi/NEWS) or run `news(package = "methylumi")`.


[MetMashR](/packages/MetMashR)
--------

                        Changes in version 1.7.3                        

- add zenodo_file to import files from Zenodo records
- vignettes now build without network access: online resources are
read from stored copies in inst/extdata/vignette_snapshots and all
lookups use cache_mode = "offline";
inst/scripts/make_vignette_snapshots.R rebuilds the stored copies
and refreshes the lookup caches
- revise all vignettes following review: new titles, clearer
explanations, worked examples, and counts computed in the vignettes
- pubchem_property_lookup replaces the deprecated PubChem property
names CanonicalSMILES and IsomericSMILES with ConnectivitySMILES and
SMILES, with a warning
- kegg_lookup caches ids with no match as NA so they are not queried
again
- rest_api lookups no longer rewrite the cache when nothing was
queried
- document model_train and model_predict with model_apply
- fix prioritise_columns assigning the last source tag to rows where
all columns are NA; the source is now NA
- fix ls_source truncating LipidName at /, + or - inside brackets
- fix mwb_refmet_database failing when the cache has duplicate
entries
for the RefMet download

                        Changes in version 1.7.2                        

- add pubchem_id_exchange for batch identifier conversion
- add classyfire_batch_lookup for batch ClassyFire queries
- add cts_lite_lookup for Chemical Translation Service queries
- add mwb_refmet_lookup and mwb_study_source for Metabolomics
Workbench
- add lipidmaps_database to import the full LIPID MAPS database
- add columns input to chebi_lookup
- add cache_mode input to rest_api lookups and kegg_lookup
- add unique input to select_max and select_min
- add case study vignette reimplementing the Metabolites Merging
Strategy
- vectorise mz_match and rt_match interval-overlap matching
- mspurity_source is now an lcms_table
- fix vertical_join matching_columns rename direction
- fix silent type coercion in database_lookup
- fix crash when importing an empty ls_source
- fix outside-label positioning in bar and pie charts
- fix BiocFileCache_database when the cache has duplicate entries
- use environment variable for GitHub API requests to avoid rate
limits

                        Changes in version 1.7.1                        

- add ChEBI lookup class

[MGnifyR](/packages/MGnifyR)
-------

                       Changes in version 1.99.0                        

Date: 2026-11-07

- Supports MGnify API v2

- Substantially faster: what took minutes takes seconds, and tens of
  thousands of analyses can be fetched in minutes

- Focuses on abundance tables: fetchStudies and fetchAnalyses return
  taxonomic and functional profiles as TreeSummarizedExperiment or
  MultiAssayExperiment

- Additional features of earlier versions, e.g., fetching files, are
  deprecated

                        Changes in version 1.9.1                        

Date: 2026-10-19

- getResult returns NULL with a warning instead of failing when no data
  could be fetched

[mia](/packages/mia)
---

                        Changes in version 1.21                         

- Fixed explained percentages on PCA axes for jointRPCA and hardened
  rank-deficient SVD handling in joint-RPCA (1.21.1, 2026-05-12)

- Generalise agglomerateByModule to non-binary numerical modules
  (1.21.2, 2026-05-31)

- Modify joint-RPCA to match with Gemelli implementation (1.21.3,
  2026-06-16)

- Implement Striped Unifrac algorithm (1.21.4, 2026-07-29)

- Added conversion functions between MAE and TreeSE (1.21.8,
  2026-09-22)

- Added applyByModule (1.21.11, 2026-10-02)

- Fixed getAlpha/addAlpha returning the first rarefaction round instead
  of
  the mean over niter rounds (1.21.12, 2026-10-19)

- getAlpha/addAlpha rarefy once per iteration and calculate all indices
  from
  the same rarefied data (1.21.13, 2026-10-13)

- Add Allen and Rao alpha diversity indices (1.21.14, 2026-10-22)

- Added getQCMetrics and addQCMetrics (1.21.15, 2026-10-26)

- Added getDysbiosis and addDysbiosis (1.21.16, 2026-10-26)

- Fixed Camargo evenness, which was computed on the wrong entries when
  samples contained zeroes and depended on feature order (1.21.17,
  2026-10-27)

- Fixed UniFrac truncating sample totals to integers, which gave NaN or
  wrong
  distances on relative abundances and percentages (1.21.18,
  2026-10-27)

- Added method = "sd" and "cv" to getTop (1.21.19, 2026-10-27)

- Formula handling in RDA, CCA and PERMANOVA no longer uses
  as.character() (1.21.20, 2026-10-27)

- transformAssay and rarefyAssay warn instead of failing when
  overwriting an assay (1.21.21, 2026-10-27)

- getMDS and addMDS use all features by default (1.21.22, 2026-10-27)

- Added getPCA and addPCA (1.21.23, 2026-10-27)

- Renamed getQCMetrics and addQCMetrics to getQC and addQC (1.21.24,
  2026-10-28)

- Added exporters (1.21.25, 2026-10-28)

- Fixed applyByModule failing on parallel workers in separate sessions,
  e.g. on Windows (1.21.26, 2026-10-29)

- transformAssay: method = "log" now computes natural log consistently
  with "log2" and "log10"; the previous vegan::decostand log
  transformation is available as method = "log_zero_replace" (1.21.27,
  2026-10-30)

- Renamed agglomerateByModule to agglomerateBySet (1.21.28, 2026-11-05)

- Renamed applyByModule to applyBySet (1.21.28, 2026-11-05)

- getPERMANOVA: binary=TRUE when Jaccard calculated (1.21.28,
  2026-11-05)

[See the full NEWS](/packages/release/bioc/news/mia/NEWS) or run `news(package = "mia")`.


[miaTime](/packages/miaTime)
-------

                 Changes in version 1.3.7 (2026-11-20)                  

- getBimodality: the coefficient uses excess kurtosis, as documented in
  Shade et al. 2014

                 Changes in version 1.3.6 (2026-11-10)                  

- HMP_2019_ibdmdb: removed three subjects without a follow-up sample

- SilvermanAGutData: added the sampling time Time.hr

- temporalMicrobiome20: abundances above 2.1e9 were missing; counts are
  now doubles

- hitchip1006: removed 11 duplicated samples

                 Changes in version 1.3.5 (2026-11-06)                  

- getStepwiseDivergence: Option to select reference samples based on
  time lag

                 Changes in version 1.3.4 (2026-11-06)                  

- Added minimum gap to getSamplingSummary

                 Changes in version 1.3.3 (2026-10-23)                  

- Added getSamplingSummary

                 Changes in version 1.3.2 (2026-10-22)                  

- Added getSurvival and addSurvival for Cox models

[miaViz](/packages/miaViz)
------

                        Changes in version 1.21                         

NEW FEATURES

- plotHeatmap: heatmap of features by samples straight from a TreeSE,
  as a
  ggplot2 object, with row and column facets from rowData/colData and
  an
  optional row tree

- plotDissimilarity: heatmap of a dissimilarity matrix between samples

- getHeatmapOrder and sort.method (in plotHeatmap and
  plotDissimilarity):
  sort the rows and columns of heatmaps

- plotOrdination for ordinations stored in reducedDims

- Joint-RPCA visualisations: plotJointRPCA, plotJointRPCACovariance,
  plotJointRPCAExplainedVariance and plotJointRPCANetwork

- plotRowTree, plotColTree: heatmap.by shows per-taxon values as a
  heatmap
  or as rings next to the tree tips

- plotSeries: each feature is drawn as its own line, as in matplot(),
  also
  without colour.by

- plotBoxplot: outlier.size, outlier.colour, outlier.alpha and
  outlier.shape
  style the boxplot outliers (shown with add.points = FALSE). The help
  page
  now lists all customisation arguments, grouped by plot layer

- New example datasets in multimedia_data_demo: se_relative,
  tse_relative and
  tse_pathway

SIGNIFICANT USER-VISIBLE CHANGES

- plotRDA, plotCCA: the vectors are scaled by a common factor so that
  the
  longest one reaches 75% of the extent of the samples (as in vegan
  biplots), instead of crowding at the origin. Their relative lengths
  and
  directions are unchanged; vec.scale = FALSE draws them at the
  original
  length, and a number sets the multiplier

- Colour-blind-friendly colours by default: discrete variables are
  coloured
  with the Okabe-Ito palette for up to 9 groups (black is used only for
  the
  9th), with the Paired palette for 10 to 12 groups and with viridis
  for
  more, instead of the Tableau palette

- plotOrdination and the functions built on it (e.g. plotJointRPCA)
  draw
  points with point.alpha = 0.7 instead of 0.4, so that light colours
  stay
  visible

[See the full NEWS](/packages/release/bioc/news/miaViz/NEWS) or run `news(package = "miaViz")`.


[MiCARA](/packages/MiCARA)
------

- Initial Bioconductor submission.

[microbiome](/packages/microbiome)
----------

                 Changes in version 1.35.2 (2026-10-25)                 

- Fixed Camargo evenness (#231). Zeroes were not removed, so the index
  was computed on the wrong entries and depended on the order of taxa.
  All Camargo values change; samples with no observed taxa now return
  NA instead of an error

- Documented that the zeroes argument of evenness() has no effect

                 Changes in version 1.35.1 (2026-08-10)                 

- Added support for TreeSummarizedExperiment and other
  SummarizedExperiment-derived objects
  - abundances(), meta() and taxa() now accept these objects alongside
    phyloseq objects, which extends support to most of the package
  - abundances() gains an assay.type argument, defaulting to the counts
    assay when present and otherwise the first assay
  - transform() stores its result as a new named assay and leaves the
    original assays untouched, following the mia transformAssay
    convention. It gains a name argument to control that assay name, and
    the result is read back with abundances(x, assay.type = name) rather
    than abundances(x)
  - Functions that manipulate the taxonomy table (aggregate_taxa,
    plot_composition, map_levels, psmelt2, the tibble utilities and the
    read_* family) remain phyloseq only

- We recommend switching from phyloseq to the TreeSummarizedExperiment
  based methods described in the Orchestrating Microbiome Analysis
  (OMA)
  online book, https://bioconductor.org/books/release/OMA/, which is
  where method development now takes place

- Passing a phyloseq object now emits a one-time deprecation notice per
  session, recommending the (Tree)SummarizedExperiment framework and
  the
  mia package. Silence it with
  options(microbiome.phyloseq_deprecation = FALSE)

- Fixed low_abundance returning an unnamed vector for phyloseq input

- Fixed collapse_replicates failing when given a single replicate field

[MicrobiomeProfiler](/packages/MicrobiomeProfiler)
------------------

                       Changes in version 1.19.2                        

- move the external-data builders for Disbiome and eggNOG from
data-raw/ into R/ so that they ship with the package and their tests
run under R CMD check instead of being skipped (2026-09-17, Thu)
- data-raw/build_disbiome.R and data-raw/build_eggnog.R are now
thin loaders, so the scheduled external-data update workflow is
unaffected
- stop advertising placeholder OG0001-style identifiers for eggNOG
(2026-09-15, Tue)
- the Example button now builds its input from the published
eggNOG artifact instead of hard-coded fixture IDs, and uses a
complete KEGG pathway that fits the analysis window, so the
generated example really produces a non-empty enrichment result
- the GSEA example also carries background identifiers from
outside the pathway, without which the permutation p-values
cannot be computed
- when the eggNOG artifact cannot be loaded the app reports the
problem instead of filling the input with identifiers that
cannot be analysed
- guard against an NA condition while the universe selector has
not rendered yet
- expose an explicit seed argument on all GSEA entry points —
gseKO(),
gseModule(), gseCOG(), gseMDA(), gseMBKEGG(), gseSMPDB(), gseHMDB()
and gseEggNOG() — forwarded to enrichit::gsea_gson() (2026-08-15,
Fri)
- set it to a number (or TRUE for a fixed default seed) to get
identical GSEA results across runs, FALSE (default) draws a
fresh seed on each run; set.seed() before the call still works
- pvalueCutoff is applied to both pvalue and p.adjust inside
enrichit::gsea_gson() (requires enrichit >= 0.2.2), matching the
historical clusterProfiler/DOSE double-filtering behavior

                       Changes in version 1.19.1                        

- setup GitHub Actions workflow for automated internal data updates
(KEGG, COG, SMPDB, HMDB, Disbiome) (2026-05-21, Thu)
- support shinyserver (2026-05-21, Thu, #3)
- import clusterProfiler to prevent R check error (2026-05-21, Thu)

[miRSM](/packages/miRSM)
-----

                     Changes in version 2.9.2-2.9.3                     

- Update Description <2026-09-03, Thus>

                        Changes in version 2.9.1                        

- Remove enirchDGN <2026-07-23, Thus>

[miRspongeR](/packages/miRspongeR)
----------

                    Changes in version 2.17.3-2.17.4                    

- Update DESCRIPTION <2026-09-03, Thur>.

                    Changes in version 2.17.1-2.17.2                    

- Remove the function of enrichDGN <2026-07-23, Thur>.

[mitology](/packages/mitology)
--------

                        Changes in version 1.6.0                        

- MitoGenesDB table updated with localization info.

- vignette updated to provide more info to use the package.

[mixOmics](/packages/mixOmics)
--------

                       Changes in version 6.38.0                        

Bioconductor 3.24, November 2026.

new features / enhancements / changes

- mixOmics is now licensed under AGPL-3 instead of GPL (>= 2) (#421)
- circosPlot is drawn by new code (#419)
- predict on multi-block models is 3 to 67 times faster, and perf and
tune.block.plsda about 2 times faster, with identical results (#369,
#383)
- perf and perf.assess on DIABLO models accept folds as a list (#400)
- mint.block.pls, mint.block.spls, mint.block.plsda and
mint.block.splsda return model objects with the same structure as
the block.* models, so predict works on them (#386)
- ?plot.perf.pls now has examples (#414)
- citation("mixOmics") gives the full author names and the DOI (#416)
- the GitHub repository moved from mixOmicsTeam/mixOmics to
mixOmics-org/mixOmics
- expanded test coverage for model utilities and internal helpers
(#381)
- the official Bioconductor package maintainer is now correctly set
as
Professor Kim-Anh Lê Cao, founding author and project lead of
mixOmics
- new title and description, reflecting the package as it is today
- mixOmics logo updated in the vignette
- ORCID iDs added for authors and contributors where known

bug fixes

- circosPlot no longer fails when a block, including the Y block, has
a single selected variable (#402, #324)
- circosPlot showed block names containing "chr" truncated and failed
with line = TRUE (#419)
- circosPlot with several blocks in blocks.link kept only about half
of the links (#419)
- plotVar no longer places points outside the correlation circle when
the plotted components are not orthogonal (#405)
- keepX and keepY of block.spls and mint.block.* models now give the
values used (#386)
- plotLoadings handles a vector xlim for contribution plots in the
graphics style (#377)
- the vignette link in the startup message returned 404; it now
points
to https://guides.mixomics.org/mixomics-guide-vignette/, and the
website and contact links use https
- ?tune.rcc gave the wrong default for BPPARAM; it is SerialParam()
- error messages now link to the mixOmics-org issue page instead of
the old mixOmicsTeam address
- help pages now link to https://mixomics.org, and references and
data
sources whose pages have moved or gone link to the paper's DOI

[See the full NEWS](/packages/release/bioc/news/mixOmics/NEWS) or run `news(package = "mixOmics")`.


[MLP](/packages/MLP)
---

                       Changes in version 1.61.1                        

- KEGG gene set description: KEGG list API update: organism -> genome

[mobileRNA](/packages/mobileRNA)
---------

                        Changes in version 1.9.3                        

- changed vignette rendering settings.

                        Changes in version 1.9.2                        

- updated mRNA mapping function; silent errors due to bugs and
updated
in dependencies.
- Updated RNAdifferentialAnalysis() fails at return: object 'groups'
not found.

                        Changes in version 1.9.1                        

- added baymobil integration

[MOFA2](/packages/MOFA2)
-----

                       Changes in version 1.24.0                        

New features

- create_mofa_from_SingleCellExperiment() can now load altExps as
additional views (alt_experiments).
- create_mofa_from_MultiAssayExperiment() can now select experiments,
assays and altExps (experiments, assays, alt_experiments).
- New mofa2() wrapper runs create_mofa() and prepare_mofa() in one
call, taking only the options that differ from the defaults.
- New plot_factor_mean_cor() heatmap to identify factors related to
intercept or library-size effects.

Significant user-visible changes

- create_mofa() and the SingleCellExperiment and MultiAssayExperiment
constructors have new arguments before groups; pass groups by name.
- create_mofa_from_SingleCellExperiment() argument assay is renamed
to
assays.
- Documentation and vignettes now recommend reticulate::py_require()
to set up mofapy2.

Bug fixes and minor improvements

- Updated to mofapy2 0.7.5, which fixes MEFISTO for some covariate
grids and adds pandas 3 compatibility.
- The warning about intercept factors now names the affected views
and
points to plot_factor_mean_cor().
- run_mofa() no longer warns when mofapy2 is a newer patch release
than the pinned version.

[Moonlight2R](/packages/Moonlight2R)
-----------

                       Changes in version 1.11.4                        

Summary

- Fix GMA() to take probe IDs from the names of the probe annotation
GRanges object, since as.data.frame() in recent GenomicRanges no
longer stores them as row names.

                       Changes in version 1.11.3                        

Summary

- Adapted plotFEA() function to parse dataFEA objects obtained with
both "ora" and "fgsea" methods

                       Changes in version 1.11.2                        

Summary

- Adapted getDataGEO() for the new GEOquery return type
"SummarizedExperiment"

                       Changes in version 1.11.1                        

Summary

- Fix FEA Moonlight Z-score to consider experimental logFC in the
score calculation and removal of fixed 'Exp.Log.Ratio' column from
DiseaseList collection

[MotifPeeker](/packages/MotifPeeker)
-----------

                        Changes in version 1.5.1                        

Miscellaneous

- Use ghcr.io/neurogenomics/bioc-meme:latest container for Ubuntu CI.

[MRManalyzeR](/packages/MRManalyzeR)
-----------
- Initial Bioconductor Submission

[msa](/packages/msa)
---

                       Changes in version 1.45.6                        

- correction of Imports in DESCRIPTION

- updated texshade.sty to newest version (1.29)

                       Changes in version 1.45.5                        

- updated vignette (some classes/methods are now imported from the
  MultipleAlignment package instead of the Biostrings package)

                       Changes in version 1.45.4                        

- fix of ClustalOmega source file squid/hsregex.c to maintain
  compatibility with
  newer versions of clang (thanks for GitHub user jeroen)

                       Changes in version 1.45.3                        

- adapted Windows makefile for ClustalOmega to ensure better
  compatibility with
  newer versions of clang

                       Changes in version 1.45.2                        

- fixes of some links in help pages (thanks to Hervé Pages)

                       Changes in version 1.45.1                        

- adapted imports to new MultipleAlignment package (Bioconductor 3.24)

- minor fix of link in help page

[MSA2dist](/packages/MSA2dist)
--------

                 Changes in version 1.17.9 (2026-11-18)                 

- update collapseXStringSet

                 Changes in version 1.17.8 (2026-06-26)                 

- shorten package title
- update vignette abstract

                 Changes in version 1.17.7 (2026-06-25)                 

BUG FIXES

- fixed Rd files

                 Changes in version 1.17.6 (2026-06-25)                 

- update DESCRIPTION

                 Changes in version 1.17.5 (2026-06-24)                 

BUG FIXES

- fixed vignette

                 Changes in version 1.17.4 (2026-06-24)                 

- update DESCRIPTION
- update dnastring2dist to handle unique sequence patterns

NEW FEATURES

- added collapseChar
- added collapseXStringSet
- added normpop
- added dnastring2weightedPi
- added rcpp_weightedPi
- added rcpp_weightedPi_pop
- added dnastring2dxy
- added rcpp_dxy_fst_pop
- added dist2dxy to calculate dxy per pop from dnastring2dist result

                 Changes in version 1.17.3 (2026-06-11)                 

- update DESCRIPTION

                 Changes in version 1.17.2 (2026-06-09)                 

- update DESCRIPTION

                 Changes in version 1.17.1 (2026-06-05)                 

- update DRESCRIPTION
- added MultipleAlignment package as dependency

[MsBackendMassIVE](/packages/MsBackendMassIVE)
----------------

- Initial Bioconductor Submission

[MsBackendMetaboLights](/packages/MsBackendMetaboLights)
---------------------

                       Changes in version 1.7.11                        

- Fix examples.

                       Changes in version 1.7.10                        

- Fix syncing data files replacing the metadata of all previously
cached files, which could then no longer be found offline (with
mtbls_cached_data_files() or offline = TRUE). The metadata is now
merged with the existing content (issue #29).
- Data files cached with previous versions that are no longer found
offline need to be synced once online again (e.g. with
mtbls_sync_data_files() using the same fileName) to become available
offline. Files are not downloaded again.

                        Changes in version 1.7.9                        

- Fix assay IDs being computed from the selected assays only
(assayName), causing files with the same name in different assays to
overwrite each other when assays were synced individually (issue
#27).

                        Changes in version 1.7.8                        

- Add sleep/delay calls to vignette to reduce failures on
Bioconductor
build systems.

                        Changes in version 1.7.7                        

- Internal refactoring to reduce need to load BiocFileCache
repeatedly.

                        Changes in version 1.7.6                        

- Add 5 second pause to examples.

                        Changes in version 1.7.5                        

- Fix installation instructions.
- Add 4 second pause to several unit tests trying to avoid
MetaboLights connection errors.

                        Changes in version 1.7.4                        

- Support/handle duplicated file names across multiple assays.

                        Changes in version 1.7.3                        

- Support parameter fileName in backendInitialize() to load only
selected file(s).

                        Changes in version 1.7.2                        

- Append the MetaboLights data set ID to locally stored/cached file
names to avoid overwriting files with same names in different data
sets.
- Enforce ordering of returned file names in alphabetic order.

                        Changes in version 1.7.1                        

- Fix issue white spaces in ftp link (issue [#16]).

                         Changes in version 1.7                         

[MsBackendMetabolomicsWorkbench](/packages/MsBackendMetabolomicsWorkbench)
------------------------------

- Initial Bioconductor Submission

[MsBackendMgf](/packages/MsBackendMgf)
------------

                        Changes in version 1.21                         

Changes in 1.21.1

- Fix TITLE of exported MGF file (issue #38).

[MsBackendMsp](/packages/MsBackendMsp)
------------

                        Changes in version 1.17                         

Changes in 1.17.2

- Small performance improvement through stringi.

Changes in 1.17.1

- Refactor the peak extraction to improve overall performance.

[MsBackendParquet](/packages/MsBackendParquet)
----------------

- Initial Bioconductor Submission

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                        Changes in version 1.25                         

MsCoreUtils 1.25.4

- Add group_mz_int() function to group peaks based on m/z similarity
considering also their intensity.

MsCoreUtils 1.25.3

- getImputeMargin() works without argument, and return the default
margins for all MsCoreUtils::impute_* function.
- Add link to QFeatures's Imputation vignette in the impute_matrix()
man page.

MsCoreUtils 1.25.2

- impute_mixed(): new split argument that controls of imputed values
are computed from the whole matrix of MAR and MNAR sub-matrices (see
issue 139).

MsCoreUtils 1.25.1

- impute_mixed(): new marArgs and mnarArgs arguments that allow to
prive different sets of arguments for the missing-at-random mar and
missing-not-at-random mnar imputation function.

MsCoreUtils 1.25.0

- New devel version.

[MsDataHub](/packages/MsDataHub)
---------

                        Changes in version 1.13                         

MsDataHub 1.13.4

- Suppress messages for the faahKO() function.
- Use GPL-3 and CC BY-SA licenses.

MsDataHub 1.13.3

- Add all files from the faahKO R package.

MsDataHub 1.13.2

- Move the Pest Mix and sciex example files to Zenodo.

MsDataHub 1.13.1

- Move Ai2025, OR11 and D19 files to Zenodo
(10.5281/zenodo.21900721).

MsDataHub 1.13.0

- New devel version.

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.39                         

MSnbase 2.39.6

- Use MsDataHub, completely drop msdata.

MSnbase 2.39.5

- Use selective import from ProtGenerics and BiocGenerics.

MSnbase 2.39.4

- Move more tests to longtests.

MSnbase 2.39.3

- Use long tests.

MSnbase 2.39.2

- Defunct deprecated functions.

MSnbase 2.39.1

- Use inherits() to check class of object.

MSnbase 2.39.0

- New devel version

[msPurity](/packages/msPurity)
--------

                       Changes in version 1.99.8                        

- Thomas N. Lawson, the original author of msPurity, is listed first
in the package authors.

                       Changes in version 1.99.7                        

- Daily package checks use compact SQLite examples and a bundled
spectral library, reducing check time and source package size.
Bioconductor's weekly long tests cover full-library matching and
compare results against the historical fixtures in
msPurityData 1.41.2.

                       Changes in version 1.99.6                        

- The new msPurity-msn-standards vignette follows L-tryptophan and
L-kynurenine from MS1 feature detection through MS2–MS4 purity
assessment, spectral trees and matching. It uses the suggested
package lcmsPlot to display chromatograms, detected peaks and
precursor purity.
- The general and LC-MS/MS matching vignettes now use Parquet studies
and xcms 4. The matching workflow compares the default score and
BLINK against a bundled MoNA library and reads the saved match
evidence from Parquet. Vignette HTML is generated from the R
Markdown sources during package builds.
- writeGraphML() returns the normalized output path on Windows after
the file has been written.
- The README and package description document the LC-MSn workflows
and
cite the BLINK method. citation("msPurity") includes both the
original package paper and the deep metabolome annotation workflow
paper.

                       Changes in version 1.99.5                        

- Build scan and precursor–fragment trees from measured spectra and
follow a branch with subsetTree(). Peak trees reuse the parent peaks
recorded by purityA() for MS3/MS4. For raw spectra, peak matching
uses 5 ppm by default. fragmentLinks() exposes measured masses and
their signed difference for fragment and candidate neutral-loss
searches. Plot or export a selected branch with plot(), as_igraph()
or writeGraphML() (#155).

                       Changes in version 1.99.4                        

- BlinkMatchingParam() selects BLINK fragment scoring for MS2–MS4
searches. Hits report approximate cosine and estimated matching-ion
count. Parquet evidence records both scores and the settings used
(#161).

- For Spectra queries, spectralMatching() reports an outcome for
every
input query. Retained hits identify shared scans and source files,
with the original measurements contributing to each stored consensus
(#160).

- spectralMatching() accepts experimental Spectra queries for
in-memory MS2–MS4 searches against experimental spectra or a
reference library. Hits carry the input positions, compatibility
reason, scores and rank (#159).

[See the full NEWS](/packages/release/bioc/news/msPurity/NEWS) or run `news(package = "msPurity")`.


[msqrob2](/packages/msqrob2)
-------

                        Changes in version 1.21                         

Changes since Bioconductor 3.23 (msqrob2 1.20).

New features

- Empirical Bayes variance trend: argument trend in msqrob(),
msqrobLm(), msqrobLmer(), msqrobAggregate() and msqrobHurdle()
(intensity component). trend = FALSE (default) or "none" keeps a
global prior variance. "mean" conditions the prior on mean
log-intensity (as limma::eBayes(trend = TRUE)). "count" conditions
it on log2 of the number of features (e.g. peptides) summarised into
each protein: the .n column of QFeatures::aggregateFeatures() for
protein-level data, the number of features per group for
msqrobAggregate() and msqrobLmer(featureGroups = ), or counts in
msqrobLm()/msqrobLmer(). This is similar in spirit to DEqMS. The
binomial models (msqrobGlm(), msqrobQB(), the count component of
msqrobHurdle()) have no trend.
- Designs that are rank deficient for a feature, e.g. when a
condition
or the reference level is not observed for a protein, are now fitted
(msqrobLm(), msqrobLmer(), also with robust = TRUE) instead of
returning fitError. Contrasts that are not estimable for that
feature are NA (using the estimability package).
- msqrobLmer() stores the estimated variance component ratios of each
fit in getModel(model)$theta (lme4's parameterisation and names).
- msqrobLmer() also accepts a data.frame as data.

Changes in results

- Robust fits without random effects or ridge (msqrob() and the
intensity component of msqrobHurdle() with such a formula, and
msqrobLm()) use the variance of the Huber M-estimator (as
MASS::summary.rlm(method = "XtX")) instead of the weighted least
squares variance. The latter treats the Huber weights as known and
underestimates the variance, which made robust tests liberal. Log2
fold changes are unchanged; standard errors are typically a few
percent larger. Robust mixed and ridge models (msqrobLmer(), and
thus msqrobAggregate()) keep the previous variance.
- Empirical Bayes variance moderation always estimates the prior with
limma's fitFDistUnequalDF1() (squeezeVar(legacy = FALSE)), on the
fitted features only. Previously limma switched to its moment
estimator fitFDist() as soon as one feature failed to fit, so the
method depended on fit failures and filtering. Posterior variances,
and hence p-values, change slightly. Requires limma >= 3.62.0.
- msqrobLm(): maxitRob now defaults to 1, as in msqrob() and the
other
functions (was 5). Results of msqrob() are not affected.

Bug fixes

- msqrobAggregate() (and msqrobLmer() with rowdata): when the formula
uses exactly one sample variable and one feature variable (e.g. ~
condition + (1 | Sequence)), the sample variable keeps its type (in
the previous release a factor became integer codes and a numeric
covariate was treated as a factor).

[See the full NEWS](/packages/release/bioc/news/msqrob2/NEWS) or run `news(package = "msqrob2")`.


[MsStash](/packages/MsStash)
-------

- Initial Bioconductor Submission

[MSstatsConvert](/packages/MSstatsConvert)
--------------

                       Changes in version 1.24.0                        

- Spectronaut: SpectronauttoMSstatsFormat now filters on q-values by
default (filter_with_Qvalue = TRUE), matching the other converters.
Pipelines that rely on the default and do not set the value
explicitly may produce different results.

[MSstatsResponse](/packages/MSstatsResponse)
---------------

                 Changes in version 1.4.0 (2026-11-25)                  

- Added a protein turnover workflow: calculateQCScore,
  calculateConfidence and classifyTurnoverProteins.

- Added the simulated_turnover_data example dataset and a protein
  turnover vignette.

- Added tracer enrichment correction to calculateTurnoverRatios.

- Updated calculatePeptideWeights to score coverage relative to
  each protein's peptides.

- Observations with weight 0 no longer affect isotonic fits.

- Fixed isotonic fits that failed when pooled observations all
  had weight 0.

[MSTree](/packages/MSTree)
------

- Initial Bioconductor Submission
 
[MultiAssaySpatialExperiment](/packages/MultiAssaySpatialExperiment)
---------------------------

- Initial Bioconductor Submission

[multiHiCcompare](/packages/multiHiCcompare)
---------------

                 Changes in version 1.31.1 (2022-07-09)                 

- Remove aggretation package dependence

- Add metap implementation of fisher p-value aggregation

- Add custom implementation of Sidak p-value aggregation

- Remove support for Lancaster p-value aggregation

[MultiOmicsBridge](/packages/MultiOmicsBridge)
----------------
- Initial Bioconductor Submission

[MultipleAlignment](/packages/MultipleAlignment)
-----------------

- Initial Bioconductor Submission

[multipletR](/packages/multipletR)
----------

- Initial Bioconductor Submission

[multipointR](/packages/multipointR)
-----------

- Initial Bioconductor Submission

[MultiRNAflow](/packages/MultiRNAflow)
------------

                       Changes in version 1.11.1                        

- R functions
- Updated 13 R functions by adding the required indentation after
#' in multi-line @importFrom directives to avoid the warning:
@importFrom must use a hanging indent to span multiple lines.
This warning was triggered when running devtools::document().
The 13 updated R functions are : DATAnormalization.R,
DATAplotBoxplotSamples.R, DATAplotExpression1Gene.R,
DEanalysisTime.R, DEplotAlluvial.R, DEplotBarplot.R,
DEplotBarplotFacetGrid.R, DEplotBarplotTime.R,
DEplotVolcanoMA.R, GSEAQuickAnalysis.R, HCPCanalysis.R,
MFUZZanalysis.R, MFUZZclustersNumber.R.
- tests/testthat
- Updated test-DEanalysisGroup.R to use withCallingHandlers() to
suppress deprecation warnings generated by the UpSetR R package.
These warnings are related to deprecated ggplot2 functions used
internally by UpSetR and are not caused by the MultiRNAflow
code.

[multistateQTL](/packages/multistateQTL)
-------------

                 Changes in version 2.5.1 (2026-10-08)                  

- Fixed error in the test for replaceNAs.
- Updated roxygen2 version

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.21.4                        

Bug fix

- P-values below ~4.9e-324, the smallest number R can store, are set
to 0. This now gives a message rather than happening silently
(#218):
- read_log_pval() counts the -log10 p-values above ~323.3 that
become P=0 and notes that the exact values are kept in the
-log10 column.
- check_small_p_val() counts text p-values (e.g. "1e-400", which
fread keeps as character) that become 0 whatever convert_small_p
is. It previously reported "0 p-values are <=5e-324" for these
and only converted P to numeric in the next check.
- read_log_pval() no longer errors when the -log10 p-value column
contains an NA ("missing value where TRUE/FALSE needed").
- The convert_small_p documentation now explains this limit.

                       Changes in version 1.21.3                        

Bug fix

- get_chain_file() (used by liftover()) now falls back to the bundled
Ensembl chain file whenever the download fails. Previously it only
did so if the failed download left an empty file behind, so a
failure that left no file stopped liftover() with a "No such file"
error (seen building the vignette on the Bioconductor nebbiolo2
builder). A failed UCSC download now gives a clear error.
- Chain files are now downloaded over https rather than ftp.

                       Changes in version 1.21.2                        

Bug fix

- pkgdown GitHub Action: skip the unused
SNPlocs.Hsapiens.dbSNP155.GRCh38 Suggests, free runner disk space
and disable the package cache, as installing every reference genome
Suggests ran the runner out of disk (#170).

                       Changes in version 1.21.1                        

Bug fix

- format_sumstats() with the default ref_genome=NULL (genome build
inference) errored in validate_parameters() since dbSNP tarball
support was added (#49, #114).
- check_chr() converted MT/chrMT to MTT, removing mitochondrial SNPs
when they were kept with rmv_chr (#116, #130, #158).
- check_on_ref_genome(): SNPs not on the reference genome were not
corrected or removed when indels were present, SNPs sharing an ID
with an indel were silently dropped and the indel removal could
exceed data.table's 2^31 join limit (#164).
- check_bp_range() now checks BP for chr-prefixed, 23/24 and lower
case CHR names (#144) and no longer silently drops rows with an
unknown CHR (#150).
- check_zscore() no longer errors when some p-values are NA (#212).
- infer_effect_column() no longer errors when swapping mixed case
allele column headers (#208).
- liftover() no longer deletes a gzipped local_chain file (#183) and
now applies rmv_chr (#211). Its documented rmv_chr default is
corrected.

[See the full NEWS](/packages/release/bioc/news/MungeSumstats/NEWS) or run `news(package = "MungeSumstats")`.


[muscat](/packages/muscat)
------

                       Changes in version 1.27.4                        

- clear 'scrapper'-related warnings, e.g., 'scater::runPCA >
  runPca.se',
  'scuttle::computeLibraryFactors > normalizeRnaCounts.se', etc.

                       Changes in version 1.27.3                        

- enable contrast-based 'lfcShrink' for 'pbDS' with DESeq2 (#151)

                       Changes in version 1.27.2                        

- added logFC shrinkage option and downstream arguments to 'pbDS'

- changed default edgeR normalization from TMM to TMMwsp

[musicatk](/packages/musicatk)
--------

                 Changes in version 2.7.19 (2026-11-15)                 

- R/benchmarking.R now passes lintr. This changes only formatting,
not
results, messages or the exported functions and their arguments
(#257).

                 Changes in version 2.7.18 (2026-11-13)                 

- benchmark() splits a composite signature (one matching several true
signatures) into pieces that keep its noise. The composite is fitted
to its true signatures with non-negative least squares, and each
mutation type's probability is divided among them in proportion to
the fit, so the "<true signature>_like" pieces add up to the
composite (a true signature with no weight in the fit gets no
piece). Each piece gets the composite's exposure times its share of
the composite, so Final RE equals the reconstruction error before
the split, also for composites that don't sum to 1. Before, a
composite was replaced by copies of the true signatures, with cosine
similarity 1, so a prediction with a composite scored better than
the correct one (#191). Composites of three or more signatures
(#189) and DBS78 composites (#190) work; before, they stopped with
an error. Summaries of predictions with composites change.
- The benchmark() summary's Num Found, Num Missing, Mean/Min/Max CS
and Sigs Found all describe the Final stage, with each true
signature counted once at its best cosine similarity, so Found +
Missing is the number of true signatures. Before, Found counted rows
of the Final comparison (a signature matched twice counted twice)
and Missing came from the Initial one (#198). sig_view_summary()
also uses each run's best cosine similarity, and lists the true
signatures in their own order.
- benchmark() merges each predicted signature into at most one group
of duplicates, that of the true signature it matches best. Before, a
signature matching two duplicated true signatures was merged into
both, counting its exposure twice (#188). After composites are
split, a piece and another predicted signature that match the same
true signature best are merged too, so the Final stage has one
signature per true signature. The summary counts each predicted
signature once: a duplicate matches a true signature that another
prediction also matches, a composite matches several, and Dup/Comp
is both; before, a signature in two duplicate groups was counted in
both. Predicted signatures whose names contain "like" are no longer
left out of duplicate merging.

[See the full NEWS](/packages/release/bioc/news/musicatk/NEWS) or run `news(package = "musicatk")`.


[MutSeqR](/packages/MutSeqR)
-------

                        Changes in version 1.1.2                        

- spectra_comparison() selects the existing G2 test with method =
"g2"
(the default) or the replicate-aware test with method = "rpdmm". The
summary report now uses the latter with sample-level counts and a
ridge penalty chosen by select_optimal_lambda() using leave-one-out
cross-validation. The replicate-aware, ridge-penalized
Dirichlet-multinomial model accounts for replicate heterogeneity and
overdispersion in high-count ECS data. The L2 penalty stabilizes
category probability estimates, and hybrid and penalized
likelihood-ratio tests are calibrated by parametric bootstrap. Both
methods share input validation, contrast handling, and count-matrix
preparation while retaining pooled G2 counts and per-sample RPDMM
counts.

- filter_mut() matches custom_filter_val as a literal substring (not
as a regular expression) in rows of custom_filter_col, including a
FILTER label among semicolon-separated labels, and warns when a
value matches no rows. Empty or missing custom_filter_val values now
produce an error. A custom_filter_col given without
custom_filter_val now stops, where it silently filtered on the
column's first value.

- spectra_comparison_rpdmm() is now internal. Use
spectra_comparison(method = "rpdmm") instead of calling it directly.

                        Changes in version 1.1.1                        

- import_vcf_data() reads FORMAT SITE_DP, the read depth at the
record's position that vartovcf writes, as the total_depth. Records
whose total_depth is missing take it from AD or DP.

- import_vcf_data() takes alt_depth from the alternate allele depths
of AD when a VCF has no VD, instead of setting every record's
alt_depth to 1. Records without an alternate allele get 0.

- import_vcf_data() no longer counts a missing AD value as 0 reads. A
record whose AD lacks the reference depth, as VarDict writes for MNV
and complex calls, takes its total_depth from DP, so its vaf is no
longer inflated toward 1.

- import_vcf_data() reads an INFO CONTEXT_3mer, as vartovcf writes
with --context-3mer, as the context, so no BS genome is needed.

- import_vcf_data() keeps the VCF FILTER value in a new filter
column,
which filter_mut() can act on with custom_filter_col = "filter".

- import_vcf_data() imports large VCFs several times faster: each
file's records are bound as data frames rather than VCF objects,
FORMAT fields with several values are spread without a per-record
loop, and each distinct REF/ALT pair and context is classified once.
Three 48,000-record gVCFs import in 4 seconds instead of 23, with
identical results.

[See the full NEWS](/packages/release/bioc/news/MutSeqR/NEWS) or run `news(package = "MutSeqR")`.


[mzR](/packages/mzR)
---

                       Changes in version 2.47.2                        

- Update controlled vocabularies to PSI-MS 4.2.14, UO
  releases/2026-10-12 and
  Unimod 17:02:2026. cv.hpp/cv.cpp are regenerated with cvgen from the
  last
  C++ ProteoWizard commit (a9ed4ce960) before its move to C#, see
  updateCvParams.sh.

[NADfinder](/packages/NADfinder)
---------

                       Changes in version 1.37.1                        

- Fix the issues in BiocCheck.

[NADIA](/packages/NADIA)
-----

- Initial Bioconductor Submission

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 3.9.4                        

- Fixed get_cgi_t2t() returning exon annotation instead of CpG
islands.
- Fixed CpG island annotations from get_cgi_hg19(), get_cgi_hg38(),
get_cgi_t2t(), get_cgi_mm10() and get_cgi_grcm39() using 0-based
starts. Starts are now 1-based like the rest of the package, so
every island starts 1 bp later than before.
- Fixed query_methy(simplify = FALSE) collapsing duplicated regions
into a single entry, and double counting their data for modBAM
input. Each input region now gets its own output entry.
- Fixed query_methy() on modBAM input returning data in the wrong
order when regions on the same chromosome were interleaved with
regions on other chromosomes.
- Fixed query_methy() on modBAM input failing when a sequence is
missing from the headers of only some BAM files. It now warns and
returns no data for the affected samples.
- Fixed query_methy() on modBAM input dropping regions with no reads.
They now return an empty table with the standard columns.
- Changed query_methy() on modBAM input to warn about sequences
missing from all BAM files and error if no queried sequence is
present, unless force = TRUE. This matches the behaviour for tabix
input.
- Changed query_methy() to require site_filter to be a single number
of at least 1.
- Fixed create_tabix_file() merging all modification types from
modkit
input into one indistinguishable set of calls. Modkit input is now
filtered to a single modification code, set with the new mod_code
argument (default 5mC, "m"). It errors if no calls match, or if
mod_code is given for non-modkit input.
- Fixed ModBamFiles() rejecting BAM files indexed with .csi instead
of
.bai.
- Fixed create_tabix_file() continuing with an unsorted file when the
system sort command failed. It now errors instead.
- Fixed create_tabix_file() default sample names keeping .tsv from
compressed inputs (e.g. sample1.tsv.gz became sample1.tsv).

                        Changes in version 3.9.3                        

- Fixed methy_to_bsseq() silently attributing methylation data to the
wrong samples when the sample annotation order differs from the
order samples appear in the methylation file.
- Fixed modBAM MM tags declaring multiple modification codes (e.g.
C+mh from dorado's 5mCG_5hmCG model) being silently dropped during
parsing.

                        Changes in version 3.9.2                        

- Fixed ModKit positions being off by 1.
- Fixed modbam_to_tabix() leaving the tabix index (.tbi) in the
temporary directory, making the converted file unqueryable.

[See the full NEWS](/packages/release/bioc/news/NanoMethViz/NEWS) or run `news(package = "NanoMethViz")`.


[ncRNAtools](/packages/ncRNAtools)
----------

                 Changes in version 1.23.3 (2026-07-22)                 

- Updated the secondary structure prediction functions to use new web
  services, following the discontinuation of the ncrna.org/rtools
  server. centroidFold now uses the ViennaRNA RNAfold web server,
  IPknot the IPknot web server at the Sato lab, centroidHomFold the
  Rfam batch search, and predictAlternativeSecondaryStructures the
  MC-Fold web server

[ndexr](/packages/ndexr)
-----

                       Changes in version 1.35.1                        

- **Maintenance: fixes for R CMD check and BiocCheck**

- added `inst/CITATION` with the Bioinformatics publication (Auer et
  al., 2018, doi:10.1093/bioinformatics/btx683)

- replaced `httr:::default_ua()` with an internal helper
  `ndex_default_ua()` using `utils::getFromNamespace()`

- imports: added `utils` and `RCX` to NAMESPACE, removed unused `plyr`
  and `tidyr` from Imports

- documentation: fixed broken `\link{ndex_connect}` references,
  replaced deprecated `@docType package` by `"_PACKAGE"`

- LICENSE file converted to the DCF stub required for `BSD_3_clause +
  file LICENSE`, full license text moved to `LICENSE.md`

- moved `ndex_api_config.yml` from `R/` to `data-raw/` and updated
  references in code documentation and vignette

- extended `.Rbuildignore` (`data-raw`, `doc`, `man-roxygen*`,
  `LICENSE.md`, Rd2pdf and BiocCheck output)

[netboost](/packages/netboost)
--------

                 Changes in version 2.21.5 (2026-06-29)                 

- nb_consensus()/nb_consensus_transfer() reworked. The module_retention
  criterion ("all"/"any"/"median"), the new min_pc1_cor threshold and
  robust_PCs are now
  included.

- Added network_type="signed" (keep only positively co-directional
  pairs) and
  consensus_method="median" for both adjacency and TOM aggregation.

- Filter-mask cache keys now include a data fingerprint, so caches are
  no
  longer reused across different datasets that share the same
  dimensions.

                 Changes in version 2.21.2 (2026-06-17)                 

- The bundled MC-UPGMA clustering suite (external C++/Perl/GNU-make
  tool) has
  been replaced by a portable in-package C++ reimplementation. netboost
  no
  longer requires GNU make, Bash, Perl or Gzip and now builds on
  Windows in
  addition to unix/macOS. Sparse-clustering output is byte-identical to
  previous versions for the supported single-round regime (edge count
  up to the
  former 1e7 heap size); very large inputs that previously triggered
  the
  external tool's multi-round mode may yield identically-structured
  trees with
  differently-numbered cluster IDs.

- SIMD (AVX) acceleration of the boosting step is now selected at
  runtime via
  the executing CPU's capabilities instead of being baked in at build
  time. The
  same binary runs correctly across heterogeneous CPU clusters (AVX and
  non-AVX nodes) and on non-x86 architectures; CPU-feature
  autodetection
  (autoconf) was removed.

- Removed the exported helper mcupgma_exec(), which only operated on
  the
  now-removed external clustering binaries.

[NetSAM](/packages/NetSAM)
------

                       Changes in version 1.53.1                        

BUG FIXES

- Wrapped the Ensembl BioMart queries in the mapToSymbol examples in
  try() so that a temporary outage of the Ensembl service no longer
  aborts 'R CMD check'.

- Replaced the removed biocLite() calls in the vignette installation
  instructions with BiocManager::install().

OTHER CHANGES

- Added biocViews terms to DESCRIPTION so that the package is
  categorized on the Bioconductor package index.

[normScore](/packages/normScore)
---------

- Initial Bioconductor Submission

[OmaDB](/packages/OmaDB)
-----

                 Changes in version 2.29.2 (2026-11-18)                 

- BUG FIXES
  - CITATION file uses bibentry() and includes the DOI of the
    publication

                 Changes in version 2.29.1 (2026-11-13)                 

- NEW FEATURES
  - caching of all API responses in memory, and optionally on disk with
    setOmaCache(persistent = TRUE). The persistent cache is cleared
    automatically when a new OMA release is available. clearOmaCache()
    empties the caches.
  - new functions getHOGs() (all HOGs / ancestral genes at a taxonomic
    level), getSynteny() (ancestral synteny), getSharedAncestry() and
    getOrthoXML(), providing similar functionality as the pyomadb package
  - getProtein() gains the attributes 'hog_derived_orthologs', 'xref'
    and
    'isoforms', and the argument rel_type to filter orthologs
  - getHOG() gains the argument attribute to retrieve the members, the
    Gene Ontology annotations and HOGs with similar phylogenetic profiles
  - getGenome() gains the attribute 'genes'
  - getTaxonomy() gains the argument collapse
  - as.data.frame() and print() methods for omadb_obj objects

- BUG FIXES
  - resolveURL() always returns a data.frame for endpoints that return
    a
    list of entries, also if the list contains a single entry (GitHub
    issue #4)
  - getProtein(attribute = 'homeologs') queried a non-existing endpoint
  - getTopGO() with format 'GO2geneID' and formatTopGO() for lists of
    proteins failed
  - getTree() keeps the full names of species containing special
    characters

- CHANGES IN RETURN TYPES AND BEHAVIOUR
  - all functions return results of consistent types: a single entity
    is
    an omadb_obj, a list of entities a data.frame - independent of the
    number of entries
  - getHOG() always returns a data.frame (one row per matching HOG)
  - loci are GRanges objects with the chromosome as seqnames and the
    OMA
    id as names (previously the OMA id was used as seqnames)
  - errors of the OMA Browser API (e.g. unknown ids) raise R errors
    with
    the explanation of the server instead of returning NULL
  - accessing an attribute of an omadb_obj that is loaded from the OMA
    Browser no longer modifies the object in the global environment;
    the data is cached instead
  - getObjectAttributes() returns a named character vector
  - getTopGO() requires topGO to be attached with library(topGO)
  - the deprecated functions getXref(), getAnnotation(),
    getGenomeAlignment() and getData() are now defunct

[See the full NEWS](/packages/release/bioc/news/OmaDB/NEWS) or run `news(package = "OmaDB")`.


[OmnipathR](/packages/OmnipathR)
---------

                        Changes in version 4.1.7                        

Bug fixes

- consensuspathdb_download and consensuspathdb_raw_table read the
current ConsensusPathDB file (version 35), which has 9 columns.
Before, they read it with 4 columns: the columns were shifted, and
the filters of min_score and complex_max_size did not work. The
UniProt ACs and Gene Symbols come from the ConsensusPathDB file,
there is no UniProt ID mapping query any more. The uniprot_a and
uniprot_b columns contain UniProt ACs instead of entry names
- consensuspathdb_download: in_complex is TRUE for records with more
than two participants, and self pairs come only from records of one
participant. Before, in_complex was TRUE for all records
- Cache: a cache file that can not be read, or that contains NULL, is
removed when loaded, and the data is downloaded or built again.
Before, loading such a file failed, e.g. "attempt to set an
attribute on NULL"
- Databases (get_db) are not saved to the cache if they are empty. On
build servers, the slow loaders return NULL while R checks
OmnipathR, and this NULL was saved and used by the later R processes
of the check
- Databases are not saved to the cache if UniProt or Ensembl failed
on
a build server, and the loader got an empty table
- uniprot_id_mapping_table: if UniProt rejects the job, the error
gives the message of UniProt. The default of
omnipathr.uniprot_idmapping_timeout is 120 seconds instead of 30

Technical

- R-universe builds (the environment variable MY_UNIVERSE is set)
count as build servers, and skip the tests of the live web service,
like the Bioconductor build machines
- The debug information at errors shows the R-universe URL

                        Changes in version 4.1.6                        

Bug fixes

- consensuspathdb_download: the UniProt ID mapping query used the ID
type labels of the retired UniProt API ("ACC", "GENENAME"), and
UniProt rejected it with HTTP 400
- uniprot_id_mapping_table: if the UniProt ID mapping job does not
finish within omnipathr.uniprot_idmapping_timeout, it stops with an
error that names the job and its status. Before, it went on, and
failed with "object 'result' not found"

[See the full NEWS](/packages/release/bioc/news/OmnipathR/NEWS) or run `news(package = "OmnipathR")`.


[omXplore](/packages/omXplore)
--------

                        Changes in version 1.7.2                        

- 1.5.1 - Replaced highcharter by plotly

- This is the first release of omXplore.

[onlineFDR](/packages/onlineFDR)
---------

                       Changes in version 2.21.3                        

MODIFICATIONS

- the help pages and vignettes now state the dependence conditions
  under
  which each procedure controls the FDR or FWER. In particular, LORD,
  SAFFRON, Alpha_investing and ADDIS do not in general control the FDR
  for
  positively dependent p-values (examples are given in
  vignette("theory")), and the original LOND (original = TRUE) also
  controls the FDR under positive dependence (PRDS)

- the help pages of ADDIS_spending and ADDIS_exhaustive state that FWER
  control is not guaranteed for positively dependent p-values, and
  those
  of Alpha_spending and bonfInfinite that they allow arbitrary
  dependence

- LONDstar and LORDstar with version = 'async' are documented to
  control
  the FDR when the null p-values are independent (Zrnic et al., 2021,
  Theorem 3); the 'dep' and 'batch' versions of the star functions are
  documented as guaranteed to control only the mFDR

- the help page of SAFFRONstar now states that version = 'async'
  follows
  Algorithm 3 of Zrnic et al. (2021) with asynchronous conflict sets,
  counting candidates among tests that started after each rejection
  became known (their Appendix Algorithm 7 counts candidates by
  decision
  time instead, which can give an undefined index); a test checks this
  against a reference implementation

                       Changes in version 2.21.2                        

MODIFICATIONS

- ADDIS_exhaustive now checks that lambda >= alpha * tau, as required
  by
  Fischer et al. (2024); this also ensures no p-value above lambda is
  rejected

- ADDIS_exhaustive now rejects alpha = 1, lambda = tau and gamma
  summing
  to more than 1, which previously gave invalid results without an
  error

- ADDIS now rejects lambda = tau, as required by Tian and Ramdas (2019)

- ADDIS with async = TRUE now follows Algorithm 3 of Tian and Ramdas
  (2019)
  exactly. Rejection times and candidate counts were previously based
  on the
  order in which tests started rather than on their decision times,
  which
  gave different (sometimes higher) test levels when tests overlapped

- LOND now checks that a user-supplied betai is at least as long as the
  number of p-values

[See the full NEWS](/packages/release/bioc/news/onlineFDR/NEWS) or run `news(package = "onlineFDR")`.


[ontoProc2](/packages/ontoProc2)
---------

- Initial Bioconductor Submission

[orthogene](/packages/orthogene)
---------

                       Changes in version 1.19.3                        

Test coverage

- Lifted overall test coverage from ~84% to ~92% by adding ~370 new
assertions across 7 new test files and 4 expanded ones. New tests
cover:
- check_species_babelgene, check_agg_args, check_keep_popular,
check_sparseMatrix, check_gene_map, check_gene_df_type,
check_gene_output, check_bool_args, check_agg_opts (helper
validation functions);
- map_orthologs_gprofiler non-chunked, chunked-serial, and
chunked-parallel paths plus null/empty-input and
same-source/target early returns and chunk_size validation;
- prepare_tree local-newick else-branch and age_max calibration
paths;
- aggregate_mapped_genes transpose=TRUE, within-species,
sort_rows, as_integers, and the early-exit "rename only" branch;
- infer_species test_species=NULL defaults and
method-as-test_species expansion plus make_plot=FALSE;
- create_background use_intersect=FALSE (union), user-supplied bg,
and within-species fast path;
- run_benchmark_once cache-hit, run_convert_orthologs=FALSE, and
within-species (is_human) branches;
- aggregate_rows_monocle3 mean and count aggregation;
- rotate_clades, is_gha, load_rda, cache_dir,
gconvert_target_opts, non121_strategy_opts, is_human,
is_matrix/is_sparse_matrix, is_converted, add_rowcol_names,
messager, format_species, methods_opts,
gene_input_opts/gene_output_opts, invert_dictionary, and
get_orgdb_gprofiler (use-local + use-remote) paths;
- Error/skip paths in filter_gene_df, aggregate_rows,
map_orthologs, all_genes, and many2many_rows
(aggregate_orthologs=FALSE).

                       Changes in version 1.19.2                        

Bug fixes

- test-map_genes: replace message(e) with
message(conditionMessage(e))
in the planarian (Schmidtea mediterranea) tryCatch error handler.
Same pattern as the 1.17.4 fix to run_benchmark_once(): passing a
condition object to message() re-signals it, and under testthat's
calling handlers an error condition signaled this way escapes the
tryCatch and is reported as a test failure. This caused
test-map_genes to fail on the GitHub Actions Windows runner whenever
the planosphere mapping file failed to download. Same fix has been
backported to RELEASE_3_23 as 1.18.1.

[See the full NEWS](/packages/release/bioc/news/orthogene/NEWS) or run `news(package = "orthogene")`.


[PACE](/packages/PACE)
----

- Initial Bioconductor submission.

[pairedGSEA](/packages/pairedGSEA)
----------

                       Changes in version 1.13.1                        

Bug fixes

- paired_ora(): the combined (paired) analysis now uses genes tested
for both expression and splicing as its universe; previously genes
without a splicing result were included. Missing paired results are
now filled in the paired columns instead of overwriting the splicing
relative risk and enrichment score. Affects paired ORA results
since 1.1.0 (Bioconductor 3.18); expression-only and splicing-only
ORA are unchanged.
- paired_diff(): test = "Wald" no longer fails; the reduced model is
only used for the likelihood ratio test.
- The reduced model for the likelihood ratio test now removes the
group term by name instead of the first term of the design.
- paired_diff(): surrogate variable estimation is now seeded
(sva_seed), so results are reproducible between runs.

Changes

- paired_diff(): the number of automatically estimated surrogate
variables is capped so that the model keeps enough residual degrees
of freedom (sva_min_residual_df). An optional upper limit can be set
with sva_max_n; when set, the number is also capped at the square
root of the number of samples. An exact request (sva_n) stops with
an informative error when it is unsafe instead of being reduced
silently.
- paired_diff(): SVA diagnostics (requested, estimated and retained
surrogate variables, residual degrees of freedom) are stored in
metadata(result)$pairedGSEA_sva.

[panoramic](/packages/panoramic)
---------

- Initial Bioconductor Submission

[panpathAudit](/packages/panpathAudit)
------------

- Initial Bioconductor Submission

[path2path](/packages/path2path)
---------

- Initial Bioconductor Submission

[pathlinkR](/packages/pathlinkR)
---------

                        Changes in version 1.99                         

- Added support for mouse data for all functions through argument
species. Impacted functions include:
- eruption()
- pathnetFoundation() and pathnetCreate()
- pathwayEnrichment()
- plotFoldChange()
- ppiBuildNetwork() and ppiEnrichNetwork()
- Data objects for mouse including:
- exampleDESeqResultsMM
- innateDbPPIMM
- keggDatabaseMM
- mappingFileMM
- pathwayCategoriesMM
- reactomeDatabaseMM
- sigoraDatabaseMM
- sigoraExamplesMM
- Human data objects are still available with the "HS" suffix, e.g.
exampleDESeqResultsHS is the same as previous exampleDESeqResults

                        Changes in version 1.9.2                        

- Added check for no genes passing filter in pathwayEnrichment

[PCAtools](/packages/PCAtools)
--------

                       Changes in version 2.26.0                        

- revert prior change to encircle functionality that was added due to
  `ggalt` being removed from CRAN.
  The prior functionality was implemented directly, so encircle
  behavior now work as it did in v2.22.0 and prior.

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.27.5 (2026-11-02)                 

- Adapt to new MSnBase / mzR behaviour

- BiocCheck fixes, updated documentation and dependencies versions

[Pedixplorer](/packages/Pedixplorer)
-----------

                        Changes in version 1.9.3                        

- Fix vignette failing in R.universe
- Fix linting and reduce complexity
- Fix duplicated imports

                        Changes in version 1.9.2                        

- Update website releases
- Fix News and Readme

                        Changes in version 1.9.1                        

- Update licence to GPLv3
- Fix contributing guidelines
- Update to Roxygen2 8.0.0

[pepVet](/packages/pepVet)
------

- Initial Bioconductor Submission

[pgen2gds](/packages/pgen2gds)
--------

- Initial Bioconductor Submission

[phenomis](/packages/phenomis)
--------

                       Changes in version 1.15.10                       

- minor updates

                       Changes in version 1.15.8                        

- reducing: update of the mzdiff_db.csv file

                       Changes in version 1.15.6                        

- minor corrections

                       Changes in version 1.15.4                        

- hypotesting: computation of 'fold change' now adapted to the
  transformation of the data (as a result, a previous call to
  'transforming' is now required)

- plot_hypotesting: new external method to plot volcano and
  boxplots/dotplots resulting from a statistical test

- methods for ExpressionSet and MultiDataSet objects now deprecated

                       Changes in version 1.15.2                        

- annotating method deprecated (call to ChEBI API not always
  responding)

[phylobar](/packages/phylobar)
--------

- Initial Bioconductor submission.

[PinPath](/packages/PinPath)
-------

- Initial Bioconductor Submission

[planttfhunter](/packages/planttfhunter)
-------------

                       Changes in version 1.13.1                        

This release adds two new classification schemes (TAPscan and
Plant-TFClass) and redesigns the API around a simpler workflow:
annotate_domains() -> classify_tfs() -> get_tf_counts().

NEW FEATURES

- TFs can now be classified using the classification scheme of
TAPscan
v4 (Petroll et al., 2025), in addition to PlantTFDB's scheme. The
TAPscan classification algorithm (tapscan_classify.pl v4.76, with
rules v82 and coverage values v11) was ported to R and reproduces
the results of the original Perl implementation. TAPscan
classifications include a variable TAP_class indicating whether each
family contains transcription factors (TF), transcriptional
regulators (TR), or putative TAPs (PT).

- TFs are also classified into superclasses, classes, and families
from the structure-based classification of Plant-TFClass
(Blanc-Mathieu et al., 2024). Plant-TFClass classifications are
inferred from TAPscan classifications using a correspondence table
between the two systems. When a TAPscan subfamily corresponds to
more than one Plant-TFClass family (TAPscan's AP2, ABI3/VP1, and
MADS), PlantTFDB subfamilies are used to resolve the ambiguity; when
that is not possible, all possible classifications are reported,
separated by ';'.

- PlantTFDB classifications now include both families (e.g., AP2/ERF)
and subfamilies (e.g., ERF), as in PlantTFDB's classification
scheme.

- New function annotate_domains() to annotate sequences with the
domains used by PlantTFDB and TAPscan. It has a threads argument to
run HMMER with multiple threads.

- New data sets tapscan_scheme (TAPscan's classification rules),
planttfclass_scheme (correspondence between TAPscan and
Plant-TFClass), and tap_class (classes, descriptions, and references
for TAPscan families).

BREAKING CHANGES

- annotate_pfam() was replaced by annotate_domains(), which returns a
list of two data frames (PlantTFDB and TAPscan) with variables Gene,
Domain, qlen, c_evalue, hmm_from, and hmm_to.

- classify_tfs() now takes the output of annotate_domains() as input
(data frames returned by annotate_pfam() are no longer accepted),
and returns a data frame with one row per gene and variables Gene,
PlantTFDB_family, PlantTFDB_subfamily, TAPscan_family,
TAPscan_subfamily, TAP_class, PlantTFClass_superclass,
PlantTFClass_class, and PlantTFClass_family. Note that the
classifications previously reported in the variable Family (e.g.,
ERF) are now in PlantTFDB_subfamily, and PlantTFDB_family contains
PlantTFDB families (e.g., AP2/ERF). Genes assigned to more than one
PlantTFDB family have families (and subfamilies) separated by ';'.

[See the full NEWS](/packages/release/bioc/news/planttfhunter/NEWS) or run `news(package = "planttfhunter")`.


[PLSDAbatch](/packages/PLSDAbatch)
----------

                        Changes in version 2.1.1                        

- Date: 2026-06-26
- Text: Fix percentile_norm()
- Details: Fix name matching issue causing generated data to contain
NAs

[pmp](/packages/pmp)
---

                       Changes in version 1.25.1                        

- fix user-supplied reference for PQN

- reset random seed state after knn

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.21.1 (2026-07-09)                 

- Added choice of droplet volume for calculating concentrations

- Added option to download DDES formatted results from shiny app

[poem](/packages/poem)
----

                        Changes in version 1.5.1                        

- SpatialExperiment, SummarizedExperiment, elsa, sp and spdep moved
  from
  Imports to Suggests, so library(poem) loads much faster (PR #4 by
  btraven00). Install them to use ELSA() or to pass a SpatialExperiment
  to getSpatialExternalMetrics()/getSpatialInternalMetrics().

- getSpatialInternalMetrics(level="class") now only computes the
  requested metrics.

[polyICSFlow](/packages/polyICSFlow)
-----------

- Initial Bioconductor submission.

[PostChicago](/packages/PostChicago)
-----------

- Initial Bioconductor submission

[powerTCR](/packages/powerTCR)
--------

                       Changes in version 1.33.1                        

- Replaced dependency on `evmix` with `mev` for estimation of
  generalized pareto distribution.

[proBatch](/packages/proBatch)
--------

                        Changes in version 2.2.0                        

NEW FEATURES

- Added design and metadata diagnostics, grouped missing-value
  filtering
  and plots, and t-SNE and UMAP plots for matrices and
  ProBatchFeatures.

- Added unified long/wide centering and batch correction, plus
  provider-aware registered transformations for ProBatchFeatures.

SIGNIFICANT USER-VISIBLE CHANGES

- Batch correction now errors on missing values by default; use "keep",
  "drop_features", or "fill" explicitly. Invalid or conflicting
  identifiers, assay selections, lineage, and result names now also
  error.

DEPRECATED

- Deprecated center_feature_batch_medians_df(),
  center_feature_batch_means_df(), center_feature_batch_medians_dm(),
  center_feature_batch_means_dm(), correct_with_ComBat_df(),
  correct_with_ComBat_dm(), correct_batch_effects_df(),
  correct_batch_effects_dm(), and correct_with_removeBatchEffect_dm();
  use
  center_feature_batch(), correct_with_ComBat(),
  correct_batch_effects(),
  or correct_with_removeBatchEffect().

BUG FIXES

- Fixed edge cases in batch correction, sample matching, CV and
  nonlinear
  fitting, date parsing, and diagnostic plots.

[prolfqua](/packages/prolfqua)
--------

- Initial Bioconductor Submission


[pRolocGUI](/packages/pRolocGUI)
---------

                        Changes in version 2.23                         

CHANGES IN VERSION 2.23.2

- Bump for Bioc

CHANGES IN VERSION 2.23.1

- Fix remap() so all datasets are projected onto the reference PCs;
proteins are now drawn in the correct place in the aggregation app
- Fix "Download plot" in the explore app when proteins are selected
- Fix fcol = NULL (and missing fcol) in the explore app
- Compare app: "Reset colours" now resets the colour pickers, the
class profiles download now contains both datasets, and long class
names containing special characters are truncated correctly

CHANGES IN VERSION 2.23.0

- New version for devel

[ProteinBatcher](/packages/ProteinBatcher)
--------------
- Initial Bioconductor submission.

[PscanR](/packages/PscanR)
------

- Initial Bioconductor Submission

[psichomics](/packages/psichomics)
----------

                       Changes in version 1.38.1                        

- Bug fix: handle empty inputs in performPCA() to prevent crashes

[PSMatch](/packages/PSMatch)
-------

                        Changes in version 1.17                         

PSMatch 1.17.5

- Added findBetterMatch().
- Import utils::data, suggest IRanges, address other notes.

PSMatch 1.17.4

- Fix CC BY-SA link.

PSMatch 1.17.3

- Change license to GPL-3.

PSMatch 1.17.2

- Added validatePSM() and its corresponding vignette and tests.
- Added data bundled with the package: spBoekweg and psmBoekweg
- Added paper reference in CITATION

PSMatch 1.17.1

- Added parameter z and allCharges in plotSpectraPTM() and
labelFragments() and set it to 1:precursorCharge(x) by default when
it used to only display charge 1 before. See this PR.

PSMatch 1.17.0

- New Bioconductor devel.

[PTMods](/packages/PTMods)
------

                         Changes in version 1.1                         

PTMods 1.1.2

- addVariableModifications() now counts how many modifications are
already present and links that to maxMods.
- Add CC BY-SA license in vignette and update link in README.

PTMods 1.1.1

- Add license note in README.
- Added pos parameter in addFixedModifications() to allow positional
modifications as they can be defined in some search engine outputs.

PTMods 1.1.0

- New Bioconductor devel.

[punKEGGer](/packages/punKEGGer)
---------

- Initial Bioconductor Submission

[pwalign](/packages/pwalign)
-------

                       Changes in version 1.10.0                        

- No significant changes in this version.

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.23                         

QFeatures 1.23.3

- Use GPL-3 and CC BY-SA license.
- Fix colnames replacement inconsistencies (#266).

QFeatures 1.23.2

- add aggregateSamples function

QFeatures 1.23.1

- New imputation vignette.
- Replace msdata by MsDataHub.

QFeatures 1.23.0

- New devel version.

[QFeaturesGUI](/packages/QFeaturesGUI)
------------

- Initial Bioconductor submission.

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

                       Changes in version 1.31.3                        

- Tidy dependencies: MSnbase is now only listed in Depends, and the
  unused
  Suggests statmod, testthat and vdiffr have been removed

- Declare the column names used in dplyr, tidyr and ggplot2 code as
  global
  variables, so R CMD check no longer reports them as undefined

                       Changes in version 1.31.2                        

- Fix computeDiffStats contrasts: group names containing "_vs_" are no
  longer
  rewritten, and contrasts such as "(A - B) - (C - D)" now work

- computeDiffStats: numeric batchEffect columns are now treated as
  factors

- Fix coefVar cumulative fraction plot to show the proportion of
  features
  with a CV less than or equal to x

- Replace deprecated tidyverse and ggplot2 code; now requires R >=
  4.1.0,
  ggplot2 >= 3.5.0 and purrr >= 1.0.0

- Use the native R pipe (|>) instead of the magrittr pipe (%>%);
  magrittr is
  no longer imported

- Correct and extend the vignette and help pages

                       Changes in version 1.31.1                        

- Fix computeDiffStats: pass the TechRep block to lmFit so the
  duplicateCorrelation
  consensus correlation is used (previously it was silently ignored)

- Document TechRep usage in computeDiffStats and convertToMSnset

[QRscore](/packages/QRscore)
-------

                        Changes in version 1.5.2                        

- Weights for zero-inflated models now follow the exact score of the
ZINB likelihood.

                        Changes in version 1.5.1                        

- Size factors now stay aligned with their samples in two-group
tests.
- Improved numerical stability of NB/ZINB weight estimation,
including
glm.nb dispersion fits and large dispersion values.
- Corrected the variance used in the mean test, improving power for
low-count genes.
- Updated the citation to the published article (Cell Reports
Methods, 2025).
- Restored the package vignette, which was no longer being built.

[quantiseqr](/packages/quantiseqr)
----------

                       Changes in version 1.21.3                        

- Fixed a division-by-zero error when rescaling immune cell fractions
after mRNA-content correction. When all estimated immune cell
fractions were zero, this step produced NaN values. The fix now
preserves zero immune cell estimates and assigns the remaining
fraction to "Other".

[quantMSImageR](/packages/quantMSImageR)
-------------

- Initial Bioconductor Submission
 
[raer](/packages/raer)
----

                       Changes in version 1.11.3                        

- Replaced scran::summaryMarkerStats() in find_scde_sites() to avoid
deprecation warnings from scuttle functions on Bioconductor devel.

[RAIDS](/packages/RAIDS)
-----

                       Changes in version 1.11.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Integration of ADMIXTURE
- New RAIDSparam object that contains all the parameters needed to
  run the pipeline
- Update vignette to showcase how to use the RAIDSparam object

[RAMEN](/packages/RAMEN)
-----

                        Changes in version 3.0.0                        

This version reflects RAMEN's incorporation to Bioconductor. Apart
from
the removal of the seed argument described below, the analysis
functionsreturn the same results as in 2.1.2, and the error and
warning
messages are unchanged.

- The seed argument was removed from selectVariables() and
nullDistGE(). Call set.seed() before either function instead. As a
consequence, nullDistGE() no longer hands the same seed to
selectVariables() in every permutation, so the cross-validation
folds now differ from one permutation to the next instead of being
held fixed. Null distributions computed with this version therefore
differ numerically from those computed with earlier ones.
- big fix: selectVariables() no longer fails when
environmental_matrix
holds a single variable and no covariates are given.
- findVML(), selectVariables() and lmGE() were split into internal
helper functions, so that each step of their work can be read on its
own. The three are now roughly half their former length. This is a
refactoring only: all three return exactly what they returned
before.
- Every package listed in Imports now has a matching importFrom()
entry in NAMESPACE. These are collected in R/RAMEN-package.R rather
than spread across the individual function files.
- The vignette calls functions by name instead of prefixing them with
package::, and every one of its chunks now has a label.
- The vignette Introduction describes how RAMEN fits into the
Bioconductor ecosystem and how it relates to other Bioconductor
packages.
- Package metadata was updated: biocViews terms were added, the R
dependency was raised to 4.6.0, and the Title and Description fields
were revised.

Claude Opus 5 was used to assist the function refactoring, using
specific prompts for the blocks that were to be moved to internal
helpers. All modifications were manually reviewed, corrected where
necessary, and fully validated by the package developer.


[RankMap](/packages/RankMap)
-------

                        Changes in version 1.1.2                        

- Update vignette

                        Changes in version 1.1.1                        

- Add topCellTypeGenes() function

[Rarr](/packages/Rarr)
----

                         Changes in version 2.1                         

Breaking changes

- (Minor:) s3_client = argument default value in read_zarr_array(),
zarr_overview(), etc. is now set to NULL instead of missing. In
practice, we expect this grant to be invisible to most users but it
makes it easier to pass down missing values in Rarr reverse
dependencies.
- The name and configuration options for the fixed-length-ascii (|S
in
Zarr v2) and fixed-length-ucs4 (&#60;U or &#62;U in Zarr v2) data types have
been updated to null_terminated_bytes and fixed_length_utf32
respectively to match their newly specified format in Zarr v3.
- Structured data types (record arrays) now always return lists as
the
internal elements, instead of vectors as previously. This allows
structured data types to contain different data types in a single
element.
- Unless data_type is specified explicitly, integers are now written
using the smallest possible bitsize based on the array x range in
write_zarr_array().
- Attributes can no longer be written before an array or group has
been initialized at the target location. This is visible in two
situations:
- write_zarr_attributes(zarr_version = 3L) fails if ran before
write_zarr_array()/ create_empty_zarr_array() or
write_zarr_group().
- create_empty_zarr_array() (and thus write_zarr_array()) fail if
attributes already exist at the target location.
- ZSTD compression and decompression now relies on base R
memCompress() and memDecompress() functions. The bundled ZSTD
library and the custom ZSTD compression/decompression of R objects
code will be removed in a future version of this package. Please get
in touch if your R installation does not come with ZSTD support and
you rely on our custom implementation.
- User-specified compression levels are now ignored when writing.
This
change is necessary to eventually enable writing on S3, and writing
sharded arrays. It allows results in better speed performance when
writing. A patch will be provided to base R to enable custom
user-specified compression levels again in the future.
- Functions related to the DelayedArray backend have been dropped
entirely after being deprecated in the previous release. To use the
DelayedArray backend, or lazy reading or Zarr datasets, look at the
ZarrArray package.
- Empty or missing attributes are now returned as NULL instead of an
empty list by read_zarr_attributes().
- Character arrays are now written as VLen-UTF8 by default if
write_zarr_array(data_type = ) argument is left empty or "character"
and no NAs are present.

[See the full NEWS](/packages/release/bioc/news/Rarr/NEWS) or run `news(package = "Rarr")`.


[rawrr](/packages/rawrr)
-----

                       Changes in version 1.21.1                        

- Provide native rawrr executables for linux-arm64, osx-arm64,
  and win-arm64 in addition to the x64 runtimes; the runtime is
  derived from the architecture of the running R session, so
  Apple Silicon no longer requires Rosetta #96 #97.

- Use Thermo Fisher RawFileReader 8.0.42 (.NET 8) #97.

- buildRawrrExe() gains the arguments runtime, allowing
  cross-builds for all six runtimes, and sourceUrl; it mirrors
  inst/rawrrassembly/build.sh: the RawFileReader NuGet packages
  are downloaded once into a local feed and passed via dotnet
  publish --source, leaving the user's global NuGet configuration
  untouched #97.

- Remove the internal helpers .addPackages() and
  .addNupkgSource().

- .isAssemblyWorking() reports the exit status and hints (e.g.,
  macOS quarantine attribute) instead of failing with 'argument
  is of length zero' if the executable cannot be started #96.

- installRawrrExe() downloads from the 1.21.1 binary tree.

- Add tests comparing reader output of sample.raw against
  reference values (scan header, trailer, file header) and tests
  for the runtime detection.

- Update INSTALL.

[RBedMethyl](/packages/RBedMethyl)
----------

                        Changes in version 1.1.3                        

- Updated auto realization backend to HDF5
- Resolved issue with dependency on DelayedArray, which is no longer
needed.
- Updated cbind to reflect the correct S4Array dependency

[rBLAST](/packages/rBLAST)
------

                        Changes in version 1.9.4                        

Changes

- Improved error handling in predict.BLAST(). Query FASTA files are now
  checked after writing, and write failures report the affected file
  and the
  underlying error.

- BLAST exit status is now checked before parsing results. Nonzero exit
  status and missing output files produce clear diagnostic messages
  with
  the query and output locations.

- predict.BLAST() now uses data.table::fread() for faster result
  parsing when
  the suggested data.table package is installed, with an automatic
  fallback
  to utils::read.table().

- Fixed blast_db_get() to honor the baseURL argument.

[RBPEqBind](/packages/RBPEqBind)
---------

- Initial Bioconductor Submission

[RBPSpecificity](/packages/RBPSpecificity)
--------------
- Initial Bioconductor Submission

[RCM](/packages/RCM)
---

                       Changes in version 1.29.0                        

- Update R version requirement
- Restyle to align with Bioconductor requirements
- Added r-universe testing

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.27.2 (2026-06-11)                 

- Replaced remaining calls to scater::aggregateAcrossCells

                 Changes in version 1.27.1 (2026-06-05)                 

- Replaced deprecated function scater::aggregateAcrossCells

[ReactomePA](/packages/ReactomePA)
----------

                       Changes in version 1.99.2                        

- bayes_reactome() default false_negative follows the enrichit change
(0.9, was 0.1): the old rate assumed hits carry 90% of every active
program's genes, which pushed every posterior to ~0 on realistic
Reactome data -- hierarchical pathway sets with 20-30 hits among
hundreds of genes. With the realistic rate the wrapper returns a
useful selection layer (see enrichit 0.2.7 NEWS for the measured
size of the difference). Requires enrichit >= 0.2.7 (2026-10-26,
Mon)

- gsePathway() now exposes eps and forwards extra arguments (...) to
enrichit::gsea_gson(), matching the other GSEA wrappers in the
family (2026-09-22, Tue, #47)

- gsePathway() and nsePathway() now expose an explicit seed parameter
for reproducible GSEA results, aligned with enrichit::gsea_gson():
set it to a number (or TRUE for a fixed default seed) to get
identical results across runs; FALSE (default) draws a fresh seed on
each run (2026-08-15, Fri)

- prepareReactomeNetwork() is ~40x faster: the per-pathway conversion
no longer round-trips through
graphite::convertIdentifiers()/pathwayGraph() (whose per-call
mapIds() overhead was hash/SQLite-bound and capped parallel scaling
at ~4-5 effective cores). It now builds a single UNIPROT->Entrez map
up front and converts the protein edges with vectorized lookups,
keeping mclapply() parallelization over pathways. The resulting
network is identical to the previous output (same nodes and edge
values; only node ordering may differ). The default cores is now 1
(serial, ~40 s for human; options(mc.cores=...) or an explicit
cores= still enables parallel) instead of detectCores() (2026-08-16,
Sat)

                       Changes in version 1.99.1                        

- nsePathway(): network-based set enrichment analysis (NSEA) on the
Reactome reaction graph using Random Walk with Restart (2026-07-08,
Wed)
- bayes_reactome(): Bayesian term selection for compressing redundant
Reactome pathway output via enrichit::bayes_enrich() (2026-07-08,
Wed)
- prepareReactomeNetwork(): utility to build and cache the Reactome
reaction graph for reuse in nsePathway() (2026-07-08, Wed)
- enrichPathway() and gsePathway() now accept a weight parameter for
weighted ORA/GSEA (2026-07-08, Wed)
- migrate caching from custom environment (reactome_env) to
yulab.utils::with_cache(); delete get_Reactome_DATA(), getALLEG()
(2026-07-08, Wed)
- gson_Reactome() now caches the final GSON object directly, with
data
cleaning inlined (2026-07-08, Wed)

[RegEnrich](/packages/RegEnrich)
---------

                       Changes in version 1.23.3                        

- Fix a bug caused by downloading data.

                       Changes in version 1.23.2                        

- Fix a bug caused by upgrade of GEOquery package (getGEO function)

                       Changes in version 1.23.1                        

- Fix a bug caused by upgrade of DOSE package (enricher_internal
function)

[regionalpcs](/packages/regionalpcs)
-----------

                       Changes in version 1.11.1                        

- Reduced repeated mapping scans and matrix conversions in
compute_regional_pcs().
- Added input checks and consistent handling of region-map columns in
summarize_region() and compute_regional_pcs().
- Moved RMTstat to Imports and dplyr and tibble to Suggests.
- Updated the vignette and function documentation, including input
scales, variance fractions, and component-selection assumptions.
- Added regression tests for both selection methods and regional
output order.

[reglScatterplotR](/packages/reglScatterplotR)
----------------

- Initial Bioconductor Submission

[rfaRm](/packages/rfaRm)
-----

                 Changes in version 1.25.3 (2026-11-17)                 

- Wait for completed sequence search results even when the Rfam API
  returns HTTP 200 early, and parse the checked response

- Fix SVG and sequence region retrieval by omitting unsupported Accept
  headers

- Fix PDB mapping JSON retrieval and select fields by name

- Updated maintainer email addresses

                 Changes in version 1.25.2 (2026-07-22)                 

- Imported readRNAMultipleAlignment from the new MultipleAlignment
  package

                 Changes in version 1.25.1 (2026-06-10)                 

- Updated sequence search to use the new Rfam batch search API

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.58.0                        

New features

- `H5Fopen()` now supports opening remote files on S3. Thanks to
  a PR contribution from Michael Sumner
  (https://github.com/mdsumner).

- New h5getAllChunkInfo() function iterates over all stored
  chunks of a chunked HDF5 dataset, returning a list with
  elements offset (chunk origins as an integer matrix),
  filter_mask, addr (file byte offsets), and size (compressed
  chunk sizes). Requires chunked storage; contiguous datasets
  raise an error Thanks to PR #217 from Michael Sumner
  (https://github.com/mdsumner).

- Variable length string datasets now support writing and reading
  NA_character_ values. These values will be read as the empty
  string "" in Python.

- It is no longer necessary to pass a file name to h5set_extent()
  when providing a dataset handle.

Deprecations & breaking changes

- The `H5Pget_version()` function was removed entirely after
  being deprecated in version 2.56.0.

- It is no longer possible to pass a missing value or "default"
  to the `bit64conversion` argument in `H5Aread()` and
  `H5Dread()`.

- The `native` argument to `H5Pcreate()` has been removed
  entirely after being deprecated in version 2.56.0.

- The `h5tid` argument to `H5Pset_lzf()` is now deprecated and
  will be removed in the next release.

- Writing `NA_character_` values to fixed-length string datasets
  is now deprecated and will be removed in a future release.

- Variable length strings is now the default behavior for
  `h5write()`, `h5writeDataset()` and `h5writeAttribute()`. This
  makes them consistent with the behaviour of `h5create()`,
  `h5createDataset()` and `h5createAttribute()`.

- Logical attributes are now always written as a TRUE/FALSE/NA
  enum even if no NA are present in the dataset, rather than
  TRUE/FALSE in some cases and TRUE/FALSE/NA in others. As a
  consequence, the checkForNA argument to h5writeAttribute() is
  now deprecated and will be removed in a future release.

- The callGeneric argument to h5read() is now deprecated and will
  be removed in the next version. Please get in touch if you rely
  on this feature. As far as we know, this feature has never been
  used.

Bug fixes

- H5T_enum datatypes (factors in R) are now read correctly when
  the drop argument in H5Dread() is set to TRUE.

- native = TRUE now works correctly on large arrays with more
  than 2^31-1 elements.

- BLOSC and LZF compression no longer fail when applied to
  variable length types. Thanks to Luke Zappia for the report and
  reproducible example.

[See the full NEWS](/packages/release/bioc/news/rhdf5/NEWS) or run `news(package = "rhdf5")`.


[rhdf5filters](/packages/rhdf5filters)
------------

                       Changes in version 1.26.0                        

BREAKING CHANGES

- H5Z_filter_blosc() and H5Z_filter_lzf() no longer expect
outbuf_size
in cd_values. This has no visible consequence for H5Z_filter_lzf()
but requires shifting the parameters that previously followed it in
H5Z_filter_blosc().

BUG FIXES

- This package now compiles on Windows ARM. Thanks to Jeroen Ooms for
reporting the issue and providing a patch.
- This package now compiles on WASM, making it available in webR.


[Rhtslib](/packages/Rhtslib)
-------

                       Changes in version 3.10.0                        

- No changes in this version.

[RhyME](/packages/RhyME)
-----

- Initial Bioconductor Submission

[RiboDiPA](/packages/RiboDiPA)
--------

                       Changes in version 1.21.1                        

- Fixed vignette build error with recent data.table: Bioconductor
  objects (GAlignments, GRanges, GRangesList) are now converted with
  as.data.frame() before as.data.table()

[ribosomeProfilingQC](/packages/ribosomeProfilingQC)
-------------------

                       Changes in version 1.25.1                        

- import cigarillo package

[rigvf](/packages/rigvf)
-----

                        Changes in version 1.5.3                        

- Dealt with changes in the Catalog API causing vignette issues

                        Changes in version 1.5.2                        

- Added portal_files() to query files from the IGVF Data Portal by
content type

[RJMCMCNucleosomes](/packages/RJMCMCNucleosomes)
-----------------

                       Changes in version 1.37.1                        

BUG FIXES

- plotNucleosomes() now plots the nucleosome positions when a
GRangesList with a single element is given.
- plotNucleosomes() legend is placed above the plot and each entry is
sized to its own label so long names are no longer truncated.
- splitAsList() is now imported from S4Vectors.

CODE IMPROVEMENTS

- Native routine registration is generated by
Rcpp::compileAttributes(); the hand-written init.c file has been
removed.
- Use of is() instead of class() == for class checks.
- Compliance with BiocCheck recommendations.
- Continuous integration and test coverage with GitHub Actions.

DOCUMENTATION

- New Installation section in the vignette.
- Corrections and clarifications in the vignette and the function
documentation.
- ORCID identifiers added for the authors, and URL and BugReports
fields updated.

[RmzTabM](/packages/RmzTabM)
-------

- Initial Bioconductor Submission.

[rnaSentry](/packages/rnaSentry)
---------

- Initial Bioconductor Submission

[RNAseqCovarImpute](/packages/RNAseqCovarImpute)
-----------------

                 Changes in version 1.11.1 (2026-11-13)                 

- Updated combine_rubins regression test to accommodate changes in
  predictive mean matching in mice 3.20.0

[rols](/packages/rols)
----

                         Changes in version 3.9                         

rols 3.9.3

- Point URL and BugReports at the current repository and
documentation
site. The previous URL no longer resolved.
- Stop excluding NEWS.md from the built package, so that package news
is picked up for the Bioconductor release announcements.

rols 3.9.2

- Remove the startup deprecation warning, left over from before the
package was transferred.
- olsTerm() now looks terms up by their identifier instead of
reconstructing their IRI from the ontology's fileLocation. The
latter only handled ontologies served from the EBI and OBO PURL
servers, and failed with "Unknown fileLocation" for 75 of the 283
ontologies in OLS4.
- Terms are no longer rejected as invalid when OLS4 omits an empty
field (such as ontology_prefix) from its response; the corresponding
slot is set to NA instead.
- Drop the 1000 row cap in OlsSearch(). OLS3 enforced it server-side,
OLS4 does not, and it was applied to the rows slot but not to the
query URL, leaving the two out of sync.
- Fix the OlsSearch() example, which failed R CMD check when
rbind()-ing pages of search results that had different columns
(contributed by Levi Waldron, lwaldron, in #52).
- Document that OLS4 omits empty fields from search responses, so the
columns of a coerced data.frame vary between queries.
- Use https rather than http for the search and properties endpoints.

rols 3.9.1

- Transfer package and maintenance to Sehyun Oh (shbrief). Thank you
for taking the package over!

[RPBC](/packages/RPBC)
----

- Initial Bioconductor Submission
 
[rprimer](/packages/rprimer)
-------

                       Changes in version 1.17.2                        

- Regenerated documentation with roxygen2 8.0.0.
- Updated imports following the split of MultipleAlignment from
Biostrings.

Changes in version 0.99.0 (2021-10-08)

- Package was submitted to Bioconductor Changes in version 0.99.5
(2021-11-12)
- Package was accepted to Bioconductor

[rpx](/packages/rpx)
---

                        Changes in version 2.21                         

rpx 2.21.2

- Fix unit test to address remote changes.

rpx 2.21.1

- Fix unit test to address remote changes.
- Update licence to GPL-3

rpx 2.21.0

- New devel version

[RVS](/packages/RVS)
---

                       Changes in version 1.35.0                        

- New option useUnaffected to the RVsharing function, a logical value
  indicating whether to include unaffected relatives (pedigree members
  with affected=0) in the computation of probabilities of rare variant
  sharing by carriers. The condition that the variant be seen in at
  least one affected relative remains.

- Correction of subject identifiers in ex.ped.mat data object.

- Corrections of wrong argument names in help files.

[S4Arrays](/packages/S4Arrays)
--------

                       Changes in version 1.14.0                        

NEW FEATURES

- Default abind() method now handles heterogeneity
  The default abind() method now accepts an heterogeneous mix of
  objects
  like any combination of ordinary arrays, Matrix/sparseMatrix
  derivatives,
  SparseArray derivatives, and DelayedArray derivatives.
  As a consequence rbind() and cbind() now also handles heterogeneity.

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.52.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Refactor as.data.frame.DataFrame(). See #138.

- Small tweaks to as.data.frame() methods for Vector and Rle objects.

DEPRECATED AND DEFUNCT

- S4Vectors:::anyMissing() is now defunct (was deprecated)

BUG FIXES

- Several bug fixes by various contributors. See #129, #136, #140,
  #144, and #145.

[SAIGEgds](/packages/SAIGEgds)
--------

                       Changes in version 2.14.0                        

- SPA root-finding hardening in `SPATest.cpp`: the saddlepoint equation
  solvers `Saddle_Prob()` and `Saddle_Prob_Fast()` now use a bracketed
  safeguarded Newton (rtsafe) fallback for `getroot_K1()` /
  `getroot_K1_fast()`. Since K1_adj is strictly increasing (dK1/dt = K2 > 0)
  the root is unique and bracketable, so the bisection safeguard
  guarantees
  convergence. The fallback is hybrid: the existing fast Newton runs
  first
  and rtsafe is invoked only when it reports non-convergence, so there
  is no
  added cost on the common path. This removes the rare non-convergence
  that
  silently reverts the variant to the anticonservative normal
  approximation

- `seqFitNullGLMM_SPA()` and `seqAssocGLMM_SPA()` support time-to-event
  (survival) outcomes via `trait.type="survival"` and a new
  `event.time`
  argument, fitting a Cox proportional-hazards frailty model through
  the
  Breslow Cox-Poisson equivalence and applying a Poisson saddlepoint
  approximation in the single-variant score test (the GATE method, Dey
  et
  al. 2022); the Cox risk-set information correction (W - U) is used in
  the
  variance-ratio estimation via a Woodbury identity

- a new vignette "Set-based Rare-Variant Association Analysis in
  SAIGEgds", covering Burden, SKAT, ACAT-V and ACAT-O tests with full
  and
  sparse GRMs

- a new vignette "Time-to-event (survival) analysis in SAIGEgds"

[sangeranalyseR](/packages/sangeranalyseR)
--------------

                       Changes in version 1.23.0                        

Performance

- BiocParallel-backed per-SangerRead construction loop. Replaces all
parallel::mclapply call sites with BiocParallel::bplapply. The
constructor now accepts a BPPARAM = bpparam() argument and a
.resolveBPPARAM(processorsNum, BPPARAM) helper maps the legacy
integer processorsNum argument onto the right backend (SerialParam
for 1 worker, MulticoreParam on Unix and SnowParam on Windows for
≥ 2 workers). Cross-platform parallel is now first-class.

- Lazy 3-frame amino-acid translation (new lazyAA = TRUE default).
When no refAminoAcidSeq is supplied, Biostrings::translate is no
longer called eagerly during SangerRead construction. The slots
primaryAASeqS{1,2,3} start as empty AAStrings; the new exported
accessors primaryAASeqS1(), primaryAASeqS2(), primaryAASeqS3()
compute on demand. Eliminates ~35 % of construction wall time on
protein-coding ABIF reads.

- Rcpp port of the peak-detection inner loop. New src/peakvalues.cpp
exports peakvalues_cpp and peakvalues_batch_cpp. The latter
processes all peak windows for one channel in a single .Call,
eliminating ~22,400 per-window R-to-C++ marshalling round-trips per
SangerAlignment build. Saves ~320 ms per build on the bundled
fixture (≈ 1.28× end-to-end).

Cumulative end-to-end timings on the bundled
Allolobophora_chlorotica/ACHLO fixture (8 reads, 4 contigs, mean
of 5 repetitions, single thread)

Full methodology and raw artifacts under
plans/05_e2e_validation_report.md, plans/06_scaling_summary.md, and
plans/07_rcpp_optimization_log.md.

- Removed redundant mclapply parallelism around nPairwiseDiffs and
oneAmbiguousColumn — the per-element work was sub-millisecond and
fork-setup overhead dominated. They are now serial lapply calls and
~5 % faster on small alignments.

New features

- chromatogram_plotly(obj, max_points = 8000, showtrim = FALSE,
colors
= "default") — interactive Plotly htmlwidget rendering of Sanger
chromatograms via scattergl (WebGL). Uniform-stride downsampling
caps points-per-channel; the returned widget carries a
downsample_info attribute reporting the original / rendered counts.
- globalTrimApp(SA) — Shiny gadget that exposes M1 / M2 trim sliders
across an entire SangerAlignment. Each "Apply" click calls
updateQualityParam(SA, ...) which cascades to every child read; live
previews of consensus length, contig count, and per-contig stats
update reactively. Returns the re-trimmed SangerAlignment on "Done".
- primaryAASeqS1(sr), primaryAASeqS2(sr), primaryAASeqS3(sr) — lazy
AA
accessors that return the cached slot when populated, otherwise
compute on demand via the existing calculateAASeq helper.

[See the full NEWS](/packages/release/bioc/news/sangeranalyseR/NEWS) or run `news(package = "sangeranalyseR")`.


[sbivar](/packages/sbivar)
------

- Initial Bioconductor Submission

[scanMiR](/packages/scanMiR)
-------

                       Changes in version 1.19.1                        

- fixed "object 'n' not found" bug

- fixed gapped vs m8 confusion in some corner cases

[SCArray](/packages/SCArray)
-------

                       Changes in version 1.21.3                        

- bug fix: `rowProds()` and `colProds()` returned NaN instead of 0 when
  the
  running product overflowed before a zero (Inf * 0); the product is
  now
  accumulated as mantissa and exponent, so it can not overflow, and the
  result no longer depends on the precision of 'long double'

- the unit tests compare products with matrixStats'
  `method="expSumLog"`,
  since the "direct" reference (base::prod) is platform-dependent on
  overflow; this fixes the test failure on x86_64 and aarch64 Linux

                       Changes in version 1.21.2                        

- register the native routines (src/init.c) and stop searching symbols
  dynamically; .Call() now checks the number of arguments

- bug fix: `colWeightedMeans()` ignored the weights on a non-transposed
  SC_GDSMatrix (it called the unweighted C kernel instead of
  c_colWMeans)

- bug fix: `rowProds()` and `colProds()` failed on a GDS-backed matrix
  because the parent method was called with the unsupported `method=`
  argument; both now use the native kernels in either orientation

- bug fix: `rowWeightedVars()`, `colWeightedVars()`, `rowWeightedSds()`
  and `colWeightedSds()` returned an incorrect weighted variance (the
  sum-of-squares term used the weight squared and a non-standard
  normalization); now matches matrixStats: sum_w(x-mu)^2 / (sum_w - 1)

- bug fix: `colWeightedMeans(., na.rm=TRUE)` on a sparse matrix used a
  wrong denominator (the full weight sum was not initialized)

- bug fix: memory protection leak in the C kernel `c_colProds`

- `rowLogSumExps()` and `colLogSumExps()` are now numerically stable
  (max-subtraction) and no longer overflow for large values

- `useNames=TRUE` is now the default for the row*/col* statistics,
  matching
  DelayedMatrixStats; names are kept unless `useNames=FALSE`.
  `rowRanges()`
  and `colRanges()` gain a `useNames` argument and set the row names

- `rowCollapse()` and `colCollapse()` no longer keep state in C
  file-scope
  globals (which were not reentrant and shared a column cursor between
  the
  two); the row map is held in an external pointer and the column
  offset is
  passed per block. This also removes a dangling-pointer hazard in the
  old
  colCollapse code

- new unit tests checking every overridden row*/col* statistic against
  matrixStats for both matrix orientations, with and without NAs, plus
  tests for log-sum-exp overflow, name handling, and
  rowCollapse/colCollapse

[SCArray.sat](/packages/SCArray.sat)
-----------

                       Changes in version 1.14.0                        

- fix `RunPCA()` on Seurat Assay5 objects with DelayedMatrix layers,
  which failed with Seurat (>= 5.6.0); new `RunPCA.Assay5()`,
  `RunICA.Assay5()` and `RunSLSI.Assay5()`

- fix `scNewSeuratGDS()` when the GDS file has assays other than
  "counts" and "logcounts"

[scBatchQC](/packages/scBatchQC)
---------

- Initial Bioconductor Submission

[scCertify](/packages/scCertify)
---------

- Initial Bioconductor Submission

[scCompoundDE](/packages/scCompoundDE)
------------

- Initial Bioconductor Submission

[scDblFinder](/packages/scDblFinder)
-----------

                Changes in version 1.27.13 (2026-11-12)                 

- fixed very slow processing due to scrapper forcing delayed output

                 Changes in version 1.27.9 (2026-10-05)                 

- fixed bug that led to wrong labeling of doublet origins

- changed the default dbr.sd, and set the importance of the expected
  doublet rate to decay across iterations

- added classifier for doublet origins, and a specialized vignette

- fixed multi-sample doublet stats

- made fallback on failed training not silent

- fixed issue with xgboost >=3.4

                 Changes in version 1.27.3 (2026-09-29)                 

- added option to manually set threshold

- added kNN-based score smoothing function

- added option to set xgboost nthreads

                 Changes in version 1.26.5 (2026-06-23)                 

- Changed default kNN method to Kmknn, with argument to set

                 Changes in version 1.26.4 (2026-06-11)                 

- Fixed deprecation warnings -> moved to scrapper

- scDblFinder results when using propMarkers>0 will not be identical.

                 Changes in version 1.26.1 (2026-05-26)                 

- fixed xgboost default arg change that was causing using all threads

[scDiagnostics](/packages/scDiagnostics)
-------------

                        Changes in version 1.8.0                        

- Added calculateReconstructionError() to detect out-of-distribution
anomalies using cell-type-specific PCA reconstruction errors.
- Replaced the 4 general-purpose vignettes with 4 vignettes built
around the package's manuscript: an introduction with a simulated
ground-truth example, a benchmark of
detectAnomaly()/calculateReconstructionError() against known ground
truth on the Zeisel brain dataset, and two case studies (COVID-19
PBMC scRNA-seq and MERFISH spatial colitis data) walking through the
project/detect/characterize diagnostic workflow.
- Added new example datasets (zeisel_reference_data,
zeisel_query_data, zeisel_benchmark_results, covid_reference_data,
covid_query_data, merfish_reference_data, merfish_query_data)
supporting the new vignettes.
- Added plot.calculateReconstructionErrorObject() featuring robust
visualization options (violin, boxplot, ridge, and ComplexHeatmap).
- Upgraded detectAnomaly() to resolve the curse of dimensionality by
allowing Isolation Forests to run on the union of query and
reference Highly Variable Genes (via n_hvgs) when pc_subset = NULL.
- Improved anomaly detection by switching default thresholding to a
dynamic, data-driven Median Absolute Deviation method
(threshold_method = "MAD", mad_multiplier = 2) across relevant
functions.

[scECODA](/packages/scECODA)
-------

                        Changes in version 1.1.6                        

- plot_heatmap: Handled error case if user provided no label_col
- Improved robustness: & -> && in plot_pca and plot_barplot
- Updated roxygen documentation format to newest version

                        Changes in version 1.0.1                        

- plot_boxplot: Fixed critical bug which broke the function.
- plot_boxplot: added parameter plot_signif to prevent clutter if
many
cell types and/or groups are plotted

[scFastDE](/packages/scFastDE)
--------

- Initial Bioconductor Submission 
 
[scider](/packages/scider)
------

                 Changes in version 1.12.0 (2026-11-20)                 

- Now requires R >= 4.4.0, for the base R operator %&#124;&#124;%.
- RBioFormats is now suggested rather than imported, so installing
scider no longer requires Java. readXenium() reads the morphology
image only when RBioFormats is installed.
- Fixed a bug in allocateCells() for contours of multiple cell types.
- Fixed a bug in mergeROI() where calling it a second time set the
ROIs merged by the earlier call to NA.
- Fixed a bug in plotLISA() where specifying 'overlay' always gave an
error.
- getClusters(), getSubClusters(), localMoran(), globalMoran() and
getNiche() now default to the most recent neighbour list of the type
they need. Previously they used whichever list was stored last.
- getClusters() and getSubClusters() no longer reset the user's
random
number stream. Their 'seed' argument now defaults to NULL, which
uses the current stream. withr is a new dependency.
- plotLISAscatter() now shows the global Moran's I in its title, as
returned by globalMoran().
- localMoran() now also returns the sum of the weights (W).

[scifer](/packages/scifer)
------

                       Changes in version 1.15.1                        

BUG FIXES

- Fix igblast() crashing with the cryptic error "argument of length
0"
when the igblast output contained no AIRR header line. The wrapper
now returns NULL gracefully and reports the captured igblast/conda
output so build-time failures are diagnosable. This fixes the
failing test-igblast.R check on the Bioconductor build machines.
- Make the "returns a data.frame object" test skip (instead of
erroring) when the igblast/conda environment is unavailable in the
test environment.
- Guard a division-by-zero in igblastwrap.py when zero sequences are
processed.

[scOverlay](/packages/scOverlay)
---------

- Initial Bioconductor Submission 
 
[scp](/packages/scp)
---

                        Changes in version 1.23                         

scp 1.23.4

- Fix CC BY-SA link.

scp 1.23.3

- Change maintainer.
- Update vdiffr images.
- Change license to GPL-3.

scp 1.23.2

- aggregateFeaturesOverAssays() is now defunct.

scp 1.23.1

- Update citations.

scp 1.23.0

- New Bioconductor devel

[scProkaR](/packages/scProkaR)
--------

- Initial Bioconductor Submission 
 
[scran](/packages/scran)
-----

                        Changes in version 1.42                         

- Finally removed the old convertTo() function, given that both
  edgeR and DESeq2 can operate directly on a
  SummarizedExperiment.

[scrapper](/packages/scrapper)
--------

                        Changes in version 1.8.0                        

- Respect default options from the underlying C++ libraries in
  most functions. Specifically, almost all R defaults are now set
  to NULL, and we only replace the defaults in the C++ library if
  a non-NULL argument is passed in R. This eliminates the need
  for manual synchronization between R and C++ when updating
  libraries with new defaults. We also add *Defaults() utility
  functions to report the C++ defaults for user inspection. This
  has some additional specific consequences for individual
  functions:
  
  - Updated default max.bias= in choosePseudoCount() to reduce
    log-transformation bias.
  
  - Updated default extra.work= in runPca() and scoreGeneSet()
    for better performance.
  
  - Updated default compute.median= in aggregateAcrossCells()
    to avoid unnecessary work.
  
  - Updated default use.min.width= in fitVarianceTrend() and
    modelGeneVariances() to improve trend fits for typical
    scRNA-seq datasets.

- Deprecated the num.neighbors= argument in runTsne(). Now, the
  number of neighbors is fully determined by perplexity=. For
  precomputed KNN results, an error is thrown if there are any
  mismatches between the number of supplied neighbors and
  tsnePerplexityToNeighbors(perplexity).

- Optionally add column names to the dimensions for runTsne.se(),
  runPca.se() and runUmap.se(). This is occasionally helpful for
  downstream plotting frameworks.

- Added fit.trend= option to skip trend fitting in
  modelGeneVariances().

- Check for IRLBA convergence in scoreGeneSet() and raise a
  warning if the algorithm did not converge.

- Bugfix to preserve assay names so that they match up to
  character subsets= in various compute*QcMetrics() functions.

- Added subsampleByPartition() to subsample cells within each
  partition.

- Improved support for experiments with spike-in transcripts:
  
  - Added centerSpikeInFactors() for convenient centering of
    size factors for endogenous genes and spike-in transcripts.
  
  - Added normalizeRnaCountsWithSpikeIns.se() for convenient
    normalization in the presence of spike-in transcripts.
  
  - Added chooseRnaHvgsWithSpikeIns.se() for convenient
    variance modelling against a spike-in trend.

- Renamed components.from.residuals= to center.scores.by.block=
  in runPca(), for consistency with the C++ code.

- Prevent an R session crash when an NA is supplied to block= and
  group= in functions like suggestRnaQcThresholds() and
  scoreMarkers(). Now, an error is thrown if a NA value is
  detected.

[See the full NEWS](/packages/release/bioc/news/scrapper/NEWS) or run `news(package = "scrapper")`.


[scRepertoire](/packages/scRepertoire)
------------

                        Changes in version 2.9.5                        

BUG FIXES

- clonalSizeDistribution() no longer aborts on a sample whose clones
are nearly all the same size. The discrete GPD likelihood raised an
error when the optimizer drove the scale out of range, which
discarded that threshold entirely; once every threshold in threshold
had failed, which.min() returned integer(0) and the next line
reported attempt to select less than one element in get1index. The
likelihood now returns an infinite penalty the optimizer can step
away from, so these samples fit normally. Affected combineTCR()
output on the bundled contig_list at four and five samples.
- The fallback that re-initializes the GPD after a failed fit ran
inside an error handler, so anything it raised escaped the
surrounding tryCatch() rather than marking that threshold
unfittable. A threshold with a single clone above it leaves var()
undefined and tripped this, which is why threshold = 1:5 failed
across the full eight-sample example.
- If no threshold in threshold can be fit, clonalSizeDistribution()
now says so and names the sample instead of surfacing a subscript
error.

UNDERLYING CHANGES

- Internal function conversion for clonalSizeDistribution() - removed
evmix from dependencies ahead of its CRAN deprecation. The
generalized Pareto density, distribution, quantile and random
generation functions, along with the threshold fit used to
initialize the spliced gamma-GPD model, are now internal. The
replacements reproduce the previous parameter estimates and distance
matrices exactly.

                        Changes in version 2.9.4                        

BUG FIXES

- clonalCluster(chain = "both") now resolves the two chains from the
data (Heavy/Light for BCR, TRA/TRB for TCR) instead of always asking
for TRA/TRB. The chain name drives how the V(D)J block is parsed, so
BCR data was read with the wrong layout: the heavy slot took the D
gene as its J call and the light slot took the C gene. With use.J =
TRUE (the combineBCR() default) the J filter was therefore matching
on the wrong gene. Affects combineBCR(chain = "both") as well.
- clonalCluster(group.by = ...) writes cluster IDs into the .Cluster
column. The column was renamed by position, and the grouping
variable occupies that position once group.by is set, so the output
held group labels and the real IDs were left in a trailing cluster
column.
- combineBCR() reads the cluster column from clonalCluster() by name.
It previously took the last column, which picked up an unrelated
column in the case where no sequence pair passes the threshold and
clonalCluster() returns the input untouched.

[See the full NEWS](/packages/release/bioc/news/scRepertoire/NEWS) or run `news(package = "scRepertoire")`.


[scRNAseqApp](/packages/scRNAseqApp)
-----------

                       Changes in version 1.13.37                       

- Add warning message when user using default password.

                       Changes in version 1.13.36                       

- Fix the warning about usage of set.seed issue.

                       Changes in version 1.13.35                       

- Make the download filename more safety.

                       Changes in version 1.13.34                       

- Add plots for multiple grouped paired genes.

                       Changes in version 1.13.33                       

- Count the vistors by removing the bots.

                       Changes in version 1.13.32                       

- Update the gene table only when file changed.

                       Changes in version 1.13.31                       

- Remove the duplicated notification for search.

                       Changes in version 1.13.30                       

- Using promise process to update the gene table.

                       Changes in version 1.13.29                       

- Add phyper pval for co-expression.

                       Changes in version 1.13.28                       

- Change the file format of molecule.

                       Changes in version 1.13.27                       

- Get IP from reqest string but not session.

                       Changes in version 1.13.26                       

- Fix a bug for splash screen.

                       Changes in version 1.13.24                       

- Allow user to close the splash screen.

- Add policy text.

                       Changes in version 1.13.23                       

- Update splash screen style.

                       Changes in version 1.13.22                       

- Add splash screen.

                       Changes in version 1.13.21                       

- Add gene score vs gene expression modules.

                       Changes in version 1.13.20                       

- Add stream plot for burstplot.

                       Changes in version 1.13.19                       

- Add stream plot for co-expression.

                       Changes in version 1.13.18                       

- Fix the bug in download button after remove the reactive values of
  ranges.

                       Changes in version 1.13.17                       

- Replace the binnedAverage by binAverage.

                       Changes in version 1.13.16                       

- Update the bigwig export methods.

                       Changes in version 1.13.15                       

- Fix plot for legend of co-expression.

                       Changes in version 1.13.14                       

- Fix the coordinates issue for init.

                       Changes in version 1.13.13                       

- Fix the coordinates issue for init.

[See the full NEWS](/packages/release/bioc/news/scRNAseqApp/NEWS) or run `news(package = "scRNAseqApp")`.


[scTensor](/packages/scTensor)
--------

                       Changes in version 2.23.3                        

- Fixed empty \item label in newCCSParams.Rd (R CMD check WARNING).

- Added manual page for BuildSignedCCI (R CMD check WARNING).

- Removed unused symTensor from Imports (R CMD check NOTE).

- Fixed 'no visible binding' NOTE for m and v in
  .simulateDropoutCounts.

                       Changes in version 2.23.2                        

- DGN enrichment (enrichDGN) is no longer available because enrichDGN()
  was removed from DOSE (>= 4.7.1) due to a DisGeNET licensing change.
  The dgnenrich parameter default is changed from TRUE to FALSE.
  Setting dgnenrich=TRUE now issues a warning and returns empty
  results.
  The DGN element is still present in the output for backward
  compatibility.

[scTypeEval](/packages/scTypeEval)
----------

                 Changes in version 1.1.2 (2026-07-16)                  

New Features

- Added aggregation, center, and scale parameters to run_pca()
function, allowing users to tune PCA computation for single-cell
and/or pseudobulk data.

Improvements

- Added comprehensive test coverage for the aggregation parameter in
run_pca(), including selective aggregation and error handling.
- Enhanced README documentation with link to external resources page
(tutorials, case studies, and metric explanations).

Bug Fixes

- Fixed incorrect name assignment in run_pca() when using selective
aggregation parameter—now correctly assigns names based on computed
aggregations rather than all data slots.

[seahtrue](/packages/seahtrue)
--------

                        Changes in version 1.7.2                        

- Fixed the BugReports URL in DESCRIPTION, which pointed to a
non-existent GitHub Pages path; it now points to the GitHub issue
tracker.

                        Changes in version 1.7.1                        

- Added input validation with clear error messages for
calculate_space().
- Documented the canonical injection-name contract; calculate_space()
now warns when the plot-required columns can't be produced.
- calculate_space() arguments OCR_var/ECAR_var renamed to
ocr_var/ecar_var for naming consistency (no change in behavior or
defaults).
- calculate_space() gains an optional atp_factors argument to supply
custom ATP conversion factors.
- Space/trajectory functions moved to their own source file.
- Added test coverage for revive_xfplate(), calculate_space(), and
the
space plots.
- Package hygiene: converted source files to ASCII-only, cleared R
CMD
check NOTEs, and switched internal messaging to cli.
- Regenerated the revive_output_donor_A example dataset with current
code, removing a stale dependency on the validate/settings packages
that broke R CMD check in environments without them.
- Switched the vignette output from BiocStyle::html_document to
rmarkdown::html_document and dropped BiocStyle from Suggests, so the
webR build (which cannot resolve Bioconductor-only packages) builds
without it.
- Rewrote the internal O2/pH tick-range quality check as a single
vectorized whole-plate pass instead of a per-well loop; identical
output, but fast enough to run under webR/wasm in the interactive
book.

[SEMPLR](/packages/SEMPLR)
------

                        Changes in version 1.1.2                        

- Adjusted enrichSEMs() function to only consider the highest scoring
binding orientation for each sequence/SEM combination in enrichment
calculations
- Documentation updates to scoreVariants(), enrichSEMs(), and the
vignette.

                        Changes in version 1.1.1                        

- Fixed bug in plotSEMMotifs() that was mislabeling plotted motifs
with scores for both reverse and forward orientations

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.54.0                        

NEW FEATURES

- new argument '.package' in `seqParallel()`, `seqApply()`,
  `seqBlockApply()` and `seqUnitApply()` for the packages loaded before
  calling the user-defined function, so that the functions in these
  packages can be used directly: they are loaded in each child process
  of a cluster or a `BiocParallelParam` object, or in the current R
  process before forking or when only one process is used

- `seqParallel()`, `seqUnitApply()` and `seqParallelSetup()` accept a
  `BiocParallelParam` object

- AARCH64 NEON SIMD optimization for vectorized internal functions,
  providing hardware-accelerated performance on Apple Silicon and other
  ARM64 platforms

UTILITIES

- `seqParallel()` with a `BiocParallelParam` object sends the jobs to
  the
  workers one by one and combines each result once it is received, as
  with a cluster object; the workers not started yet are started once
  for the whole call, and `.initialize` and `.finalize` are called once
  in each worker

- `seqGetData(, "position")`, "$chrom_pos", "$chrom_pos2" and
  "$chrom_pos_allele", and `seqApply()` with these variables no longer
  load the positions of all variants into memory: the positions are
  read
  from the GDS file via a cache of up to 131072 consecutive variants
  per file, which is efficient when the variants are accessed in order
  (e.g., `seqBlockApply()` and `seqUnitApply()`)

- `seqSetFilterChrom(, from.bp=, to.bp=)`, and so `seqSetFilter()` with
  GRanges or IRanges and `seqSetFilterPos()`, no longer load the
  positions of all variants into memory: only the minimum and maximum
  positions of every 131072 variants on a chromosome are kept, the
  chunks of variants out of the requested ranges are skipped without
  reading, and binary search is used in a chunk if its positions are
  sorted

- `seqFilterPush()` and `seqSetFilter(, action="push")`, also used
  internally by `seqBlockApply()`, `seqUnitApply()`, etc., store the
  pushed sample and variant filters as bit vectors instead of byte
  vectors, so each pushed filter uses 1/8 of the memory (1 bit per
  sample or variant) until it is popped back

- to reduce memory usage, `seqOpen()` no longer allocates the sample
  and variant filters (1 byte per sample or variant), since all samples
  and variants are selected; a filter is allocated only when it is
  needed (e.g., to set a filter or read data), and resetting the filter
  (e.g., `seqResetFilter()`) does not release the allocated memory,
  which is released when the file is closed

- the internal function `.buffer_position()` is removed, since the
  positions of all variants are no longer loaded into memory

- `seqVCF2GDS()` is faster on VCF files with many INFO fields

- `seqVCF2GDS()` parses diploid genotypes with single-digit alleles
  (e.g., "0|1" and "0/1") faster

[See the full NEWS](/packages/release/bioc/news/SeqArray/NEWS) or run `news(package = "SeqArray")`.


[Seqinfo](/packages/Seqinfo)
-------

                        Changes in version 1.4.0                        

- No significant changes in this version.

[seqRdna](/packages/seqRdna)
-------

- Initial Bioconductor Submission

[sesame](/packages/sesame)
------

                       Changes in version 1.31.12                       

- The Basic Usage vignette's note on the C implementation names it
  SeSAMe2 and links its current home (zhou-lab.github.io/sesame), gives
  the working conda line (sesame and yame from zhou-lab), and states
  the
  agreement with this package as measured on the shared annotation
  (QCDPB betas to about 1e-15 with sesameData >= 1.31.1).

                       Changes in version 1.31.11                       

- The startup message no longer carries the knowYourCG pointer added in
  1.31.10; the README and the Basic Usage vignette point to knowYourCG
  for enrichment of CpG sets instead.

                       Changes in version 1.31.10                       

- The startup message points at the knowYourCG package for enrichment
  of CpG sets.

                       Changes in version 1.31.9                        

- Documentation fixes from a review of the vignettes. The inference
  vignette shows its library(sesame) and caching setup, gives raw
  download URLs for the two clock files, states the mouse age the code
  returns (1.18 months) and re-identifies from ~/Downloads in both
  calls.
  The Basic Usage and QC vignettes define idat_dir as the built-in HM27
  IDATs. The non-human vignette drops an undefined tmp and uses prep
  code
  SHCDPB, as the Basic Usage table recommends for Mammal40. Ethnicity
  is
  no longer listed as an inference in DESCRIPTION or the vignette.

- sesameQCtoDF() works on HM27: the median beta statistics are unnamed,
  so an array without ch or rs probes no longer gives "row names
  contain
  missing values".

- matchDesign() (prep code M) stops with an install line when the
  suggested package RPMM is missing.

- inferSpecies() warns when fewer than 30% of the probes are detected,
  since the species call is then unreliable.

                       Changes in version 1.31.8                        

- inst/CITATION writes the one non-ASCII author name (Szab\u00f3) as an
  escape, so utils::readCitationFile() reads the file without a
  declared
  encoding. BiocCheck warned "Unable to read CITATION file".

[See the full NEWS](/packages/release/bioc/news/sesame/NEWS) or run `news(package = "sesame")`.


[sfi](/packages/sfi)
---

                        Changes in version 1.1.0                        

- Replaced the mzR dependency with RaMS for reading mzML files in
getmzml(). RaMS is a pure-R reader with no system-library
dependencies (no netCDF/HDF5), which makes installation lighter and
more portable. Output is unchanged: a data.frame with mz, intensity,
and rt (seconds) columns.

[SharedObject](/packages/SharedObject)
------------

                       Changes in version 1.27.1                        

- NEW FEATURES:
  - Use new R's C API
  - Fix a bug in the subset method for shared vectors that may cause a
    crash when the indices are out of range

[shinybiocloader](/packages/shinybiocloader)
---------------

                        Changes in version 1.4.0                        

BUG FIXES AND MINOR IMPROVEMENTS

- Enlarge the Bioconductor loader graphic (biocspin) by 25% and
center
it horizontally within the loader container
- Update installation command in README

[signeR](/packages/signeR)
------

                       Changes in version 2.15.2                        

- fix seqinfo subset error

- Add seqlength to example vcf to avoid seqinfo errors

[sigvar](/packages/sigvar)
------

- Initial Bioconductor Submission
 
[SimiCviz](/packages/SimiCviz)
--------

- Initial BIoconductor Submission

[simpleSeg](/packages/simpleSeg)
---------

                 Changes in version 1.15.2 (2026-10-14)                 

- Tiled segmentation of large images supports pca, cellBody =
  "discModel" or
  cell body markers, and transform "tissueMask". The PCA and the cell
  body
  models are fitted on pixels sampled across the whole image (up to
  100,000
  and about 1,000,000), chosen by their position so that the sample
  does not
  depend on the tiles or the number of cores.

- Image files are opened once per process, which makes reading tiles
  much
  faster.

- Channel names are read from the band names of multi-band TIFF files.

- Pixels equal to a file's "no data" value (such as 255 in 8-bit TIFFs
  written by terra) keep that value instead of becoming missing.

- Cell body models no longer fail when the cell body signal is
  constant.

                 Changes in version 1.15.1 (2026-10-14)                 

- cellBody = "discModel" and cell body markers work again for several
  images
  of the same size, and cell body markers work for a single image.

- An image where no nuclei pass the threshold gives a warning and an
  empty
  mask instead of an lm.fit error (#20).

- simpleSeg() checks that nucleus and tissue markers exist, and checks
  watershed, cellBody and cores; nucleus can be given as channel
  positions.

- normalizeCells() stops when imageID is missing and a per-image method
  is
  requested (#16, #22), and works with a single marker.

- Faster segmentation with lower memory use: each image is processed
  once in
  one task, on plain arrays, with fitted values from one QR fit.

- Parallel processing starts the BiocParallel workers once per call,
  runs
  BLAS with one thread per worker, and gives every image its own random
  number stream, so masks do not depend on the number of cores.

- simpleSeg() reads images from files: multi-channel files (OME-TIFF,
  QPTIFF, multi-page TIFF), folders of single-channel files (one image
  split
  by marker) and folders of such folders.

- New argument output: "masks", "cells" (a SpatialExperiment with the
  mean
  intensity of every channel, area and centroid of each cell, as from
  cytomapper::measureObjects()) or "none"; maskDir saves the masks as
  TIFF
  files. Paths return cells by default, images in R return masks.

[See the full NEWS](/packages/release/bioc/news/simpleSeg/NEWS) or run `news(package = "simpleSeg")`.


[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.23.6 (2026-10-08)                 

- Help pages link to other packages' topics with package anchors (for
example \link[SummarizedExperiment]{colData}), as R CMD check and
Bioconductor now require. No change to functionality.

                 Changes in version 2.23.5 (2026-10-08)                 

- runFastMNN() hides deprecation warnings raised inside batchelor's
fastMNN() and reducedMNN(), as runMNNCorrect() already did for
mnnCorrect(). batchelor still calls deprecated scuttle functions.
Other warnings are unchanged.

                 Changes in version 2.23.4 (2026-10-07)                 

- Fixed feature selection on the Celda tab of the Shiny app. Both the
Scran_modelGeneVar and SeuratFindHVG options failed: the first
called a function that does not exist, the second was never run
because its option value did not match, and both requested the top
genes with an argument that getTopHVG() does not accept. With
SeuratFindHVG, only the vst method works with Seurat 5. The selected
genes are now stored in the featureSubset alternative experiment
used by Celda; the main dataset keeps all genes, as with the None
option, instead of being reduced to them.
- Fixed starting the Shiny app with data, singleCellTK(inSCE = sce),
which failed for any object without alternative experiments.
- The SoupX help page in the app describes the new clustering method.

                 Changes in version 2.23.3 (2026-10-07)                 

- runModelGeneVar() now models gene variances with
scrapper::modelGeneVariances() instead of the deprecated
scran::modelGeneVar(). Means and total variances are unchanged, but
the biological component comes from a different trend fit, so the
ranking of highly variable genes changes (about 90-97% of the
top 500 genes are the same on the test datasets). Results that
depend on modelGeneVar HVGs, such as PCA, UMAP, and clusters, can
change.
- runSoupX() clusters cells (when no clusters are supplied) with the
steps of scran::quickCluster(), computed with scrapper and bluster.
The automatic clusters, and therefore SoupX's contamination
estimates, can change slightly.

                 Changes in version 2.23.2 (2026-10-07)                 

- Replaced functions deprecated in scuttle 1.22 and scran 1.40.
Results are unchanged (identical, or within floating-point error):
- runClusterSummaryMetrics(), plotSCEHeatmap(), and the TSCAN
functions aggregate with scrapper or base R instead of
scuttle::aggregateAcrossCells().
- scaterlogNormCounts() normalizes with
scrapper::normalizeRnaCounts.se().
- runPerCellQC() and sampleSummaryStats() compute the per-cell QC
metrics with scrapper::computeRnaQcMetrics() and singleCellTK
code (for the top-N, gene-set, and altExp columns) instead of
scater::addPerCellQC(). This is about 4 times faster on large
data. runPerCellQC() uses BPPARAM to set the number of threads.
The metadata packageVersion recorded by runPerCellQC() is now
singleCellTK's.

[See the full NEWS](/packages/release/bioc/news/singleCellTK/NEWS) or run `news(package = "singleCellTK")`.


[SingleR](/packages/SingleR)
-------

                       Changes in version 2.16.0                        

- Added top= option to combineRecomputedResults() to restrict the
  number of per-label markers used for integrated classification.
  This defaults to 10 to reduce the feature space for faster
  compute, but can be set to Inf to recover the previous results.

- Bugfix to SingleR() to respect non-default assay.type.ref=.

[smartid](/packages/smartid)
-------

                        Changes in version 1.9.2                        

- All the main wrapper functions keep the same usage.

Performance

- top_markers() no longer scales the full feature-by-cell matrix. Row
centring densifies a dgCMatrix by construction, and that dense
intermediate was the single largest allocation in the pipeline. The
scaling is now deferred to the reduced G x K / K x G statistic:
scale_mgm() applies a per-row affine map with a strictly positive
scale, and every statistic downstream of it is affine-equivariant,
so scaling the reduced result is the same operation on a matrix
thousands of times smaller. On a 5000 x 3000 input at 3% non-zero
this is 9-18x faster and uses 10-17x less memory for sparse input,
and roughly 1.2x faster with 1.4x less memory for dense input.
- mad and the 1-vs-max beta contrast are location-invariant, so they
take the scale factor only; mean, median and the raw label
coefficients take both terms.
- scale_mgm() itself is unchanged, including its return class. The
new
row_scaling_params() holds the centre and inverse SD both it and the
deferred path use, so the two cannot drift apart.
- Scores are unchanged. median and mad are bitwise identical; the
mean
and glm paths agree to within 1e-15, the round-off expected from
reassociating the arithmetic, and the selected marker sets were
identical across every scaling, aggregation, batch and softmax
combination tested.
- The deferral applies only where fit_label_betas() takes the
closed-form solve. A non-identity link is not affine in the
response, and a rank-deficient design (a batch or donor covariate
nested inside the label) is rejected by the closed form, so both
cases fall back to the previous materialised route with no change in
behaviour.

                        Changes in version 1.9.1                        

Bug fixes

- Completed the sparse-matrix support started in 1.7.1. The IAE
branch
had been converted but the IDF branch had not, so cal_score() still
aborted on a dgCMatrix for four of the eight idf_iae_methods() with
'x' must be an array of at least two dimensions. idf(), idf_m() and
idf_sd() now use sparseMatrixStats::rowSums2(), and idf_m() masks
with (expr > thres) * n_sub instead of ifelse(), which silently
drops dim on a sparse logical matrix. Affected methods: standard, m,
sd and hdb.
- gs_score_init() used base::colMeans() and therefore rejected all
sparse input, which also broke gs_score() and ova_score_boxplot().
It now uses sparseMatrixStats::colMeans2().

[See the full NEWS](/packages/release/bioc/news/smartid/NEWS) or run `news(package = "smartid")`.


[smoothclust](/packages/smoothclust)
-----------

                 Changes in version 1.9.5 (2026-07-26)                  

- simplify boundary density calculations

                 Changes in version 1.9.4 (2026-07-20)                  

- add library size normalization before spatial smoothing by default

- add log2 transformation after spatial smoothing by default

- accept optional precomputed size factors

                 Changes in version 1.9.2 (2026-05-18)                  

- rename metric to boundary density

- include normalization for number of nearest neighbors

- include optional composition adjustments

- return local and sample-level boundary density metric

[smoppix](/packages/smoppix)
-------

                        Changes in version 1.5.4                        

- Cleaning code with the help of the goodpractice suggestions

                        Changes in version 1.5.3                        

- Update citation file to main publication

                        Changes in version 1.5.2                        

- When building the vignette, use serial processing on windows

                        Changes in version 1.5.1                        

- Switch to SnowParam for multithreading on windows


[sosta](/packages/sosta)
-----

                        Changes in version 1.5.2                        

- update citation

[SpaceMarkers](/packages/SpaceMarkers)
------------

                        Changes in version 2.3.4                        

- Documentation-only release. The 2.3.3 NEWS entries below were
written after 2.3.3 had already been pushed to Bioconductor devel,
so this version exists to carry them upstream. No code changes
relative to 2.3.3.

                        Changes in version 2.3.3                        

Breaking changes

- Directed overlap scores now use target / source column names in
place of pattern / influence (relAbundance is unchanged). target is
the pattern whose own GMM hotspot is being examined; source
(prefixed near.) is the pattern whose influence zone overlaps into
that hotspot. Code that reads the old pattern / influence columns
from directed results must be updated.

Bug fixes

- SpaceMarkers() on a hand-built SpaceMarkersExperiment no longer
errors with "spatial_params(sme) is NULL and no visiumDir is stored"
— an SME assembled without load10X() now falls back to computing
default spatial parameters instead of failing.
- add_features(): fixed a positional misalignment when the SME and
the
features share the same barcodes in a different order. The NA mask
and the subsequent colData assignment were applied against the wrong
spot, so values were attached to the wrong barcodes and the wrong
spot was dropped. Features are now reordered to the SME's column
order first.
- load10X() / add_features(): closed an NA-padding gap left by the
earlier spatial-pattern NA handling.
- Directed scores no longer mix numeric and character gene labels.
When one direction of a pattern pair did not have two groups, the
all-NA fallback table in .calc_IM_scores() was built without row
names, so its gene column became row numbers ("1", "2", ...) while
the other direction used real gene names. The fallback tables now
carry gene names.
- calculate_fc_score(): reworked binning logic so low_bins and
high_bins can no longer come back empty.
- The parallel backend is no longer registered globally from inside a
closure, leaving a user's own registered BiocParallel backend
intact.

Internal

- getSpatialParameters and the image-creation path refactored;
AnnData
loader improvements; missing documentation added for several
functions.
- Tests skip gracefully when optional AnnData reader dependencies are
absent.
- No changes to the exported API in this release (NAMESPACE is
unchanged since 2.3.1).

[See the full NEWS](/packages/release/bioc/news/SpaceMarkers/NEWS) or run `news(package = "SpaceMarkers")`.


[SpaceTrooper](/packages/SpaceTrooper)
------------

                        Changes in version 1.3.3                        

- minor changes in vignettes to get lighter images

                        Changes in version 1.3.2                        

- Updaging vignette with details on xenium_small dataset

                        Changes in version 1.3.1                        

- Added the canonical computeQScore(), computeQScoreFlags(),
computeOutliersQScore(), and applyQScoreModel() APIs.
- Retained computeQCScore(), computeQCScoreFlags(),
computeOutliersQCScore(), and applyQCScoreModel() as deprecated
compatibility APIs. They preserve historical arguments, QC_score,
low_qcscore, low_threshold_qcscore, and QCScore_model* metadata.
- Restored historical plotting arguments (size, alpha, alphaNumbers,
csize, calpha, mapNumbersCol, and mapAlphaNumbers) as aliases for
the canonical plotting arguments. Canonical arguments take
precedence when both forms are supplied.
- Restored getModelFormula(formulaVars, verbose=FALSE) compatibility
while accepting metricList as a named replacement.
- Restored the public computeLambda(trainDF, modelFormula) signature,
complete-case filtering, and support for both qcscore_train and
QScore_train.
- User-supplied modelFormula values are now fitted without adding,
removing, or rebuilding terms. Supported subsets, additive formulas,
and selected interactions are preserved.
- Supported base QC metrics selected by a custom modelFormula now
determine the outliers, good and bad cells, and QScore_train used
for model training; omitted metrics no longer contribute to
training-label construction.
- Formula predictors are limited to log2SignalDensity, Area_um,
log2AspectRatio, and log2Ctrl_total_ratio. Unsupported predictors
and transformations now produce informative errors, and CosMx-only
border terms are checked against dataset technology.
- Model transfer now has a canonical public API, preserves training
matrix column order, supports validated custom formulas, and retains
the historical deprecated interface.
- Fixed the SpaceTrooper utilities vignette name and added the
bioRxiv
citation.
- Follow-up: intercept handling remains unchanged in this release and
should be reviewed consistently across training, lambda selection,
prediction, and model transfer.

[spammR](/packages/spammR)
------

- Initial Bioconductor submission.

[SpaNorm](/packages/SpaNorm)
-------

                        Changes in version 1.8.0                        

New Features

- SpaNorm() gains a BPPARAM argument to parallelise normalisation
across workers via BiocParallel, accelerating large datasets. It
defaults to BiocParallel::SerialParam() (no parallelisation), and
results are identical regardless of the backend used.

- SpaNorm() now normalises DelayedArray-backed count assays (e.g.
disk-backed via HDF5Array) block-wise, so out-of-core datasets are
processed without ever loading the full matrix into memory. Results
match the in-memory path.

- Exported fitNB(), which fits a per-gene negative binomial GLM over
an arbitrary design matrix using SpaNorm's IRLS engine (with
optional ridge regularisation and adjustable outlier winsorisation).
This exposes the model-fitting machinery for reuse independently of
SpaNorm's spatial model.

- fitNB() and calculateMu() gain an offset argument: a genes x cells
matrix added to the linear predictor with its coefficient fixed
at 1, so log(mu) = gmean + tcrossprod(alpha, W) + offset. Use it for
an effect that is already known rather than adding a column to the
design, whose coefficient would float and could absorb signal
correlated with that effect. The offset is subset alongside the
counts, applied to the dispersion estimation as well as the mean,
and sliced in step with the counts when the fit is carved into
gene-blocks. It defaults to NULL, which is a strict no-op, so all
existing behaviour is unchanged. The GPU path is implemented
backend-agnostically but has not yet been executed on an
accelerator.

- fitNB() also gains a psi argument for supplying per-gene NB
dispersions instead of estimating them. Supplied dispersions are
used as-is (no re-estimation, no winsorisation) and the outer
dispersion loop is bypassed: the coefficients come from a single
IRLS fit at the given psi, so edgeR::estimateDisp is never called
and the returned sampling factor carries no "dispersion" level. The
dispersions need not come from an identical design -- values
estimated on a design nested in (or equal to) the fitting design are
appropriate, e.g. pooled across a coarser model, which errs
conservative. psi composes with offset, and psi = NULL (the default)
is a strict no-op.

[See the full NEWS](/packages/release/bioc/news/SpaNorm/NEWS) or run `news(package = "SpaNorm")`.


[SparseArray](/packages/SparseArray)
-----------

                       Changes in version 1.14.0                        

NEW FEATURES

- rbind() and cbind() work on a mix of SparseMatrix derivatives
  and ordinary matrices.

BUG FIXES

- Fix the following stack imbalance warning:

  ```
  nms <- asNamespace("SparseArray")
  # Warning: stack imbalance in '<-', 2 then 4
  ```

- Small tweak to type() setter.
  Changing the type of an SVT_SparseArray object can replace non-one
  values
  with one values. For example this can happen when changing the type
  from
  double to integer or from complex to double. When this happens, some
  SVT
  leaves in the object might need to be turned into lacunar leaves.
  This change makes sure that changing the type of an SVT_SparseArray
  object
  switches the leaves that need it to lacunar representation.

[SpatialData.plot](/packages/SpatialData.plot)
----------------

- Initial Bioconductor Submission

[spatialdataR](/packages/spatialdataR)
------------

- Initial Bioconductor Submission

[spatialFDA](/packages/spatialFDA)
----------

                        Changes in version 1.5.2                        

- added option to apply the inverse link on the intercept function to
have it on the original scale for the pltMdl function. To recover
the old behaviour, set linkinverse to FALSE in pltMdl

                        Changes in version 1.5.1                        

- update version to match bioconductor git

[SpatialFeatureExperiment](/packages/SpatialFeatureExperiment)
------------------------

                       Changes in version 1.99.0                        

- Made affine transform specifications user-facing with the
AffineTransform S3 class
- Added option to use DuckDB for geometries in read functions

                       Changes in version 1.16.0                        

- Deprecated some functionalities in formatTxSpots(): split spots by
z-planes (and related arguments), split spots by cell compartments
in CosMX, and option to format the transcript spots for
colGeometires
- Deprecated geometryFeatureData except for annotGeometry because I
never used it for colGeometries
- Allow out of memory geometries and geometric operations with
DuckDBSpatial and DuckDBDataFrame
- Fixed bug in coercion method from SPE when the geometries aren't
flipped for SpatRaster image
- Breaking change: withDimNames in colGeometry getters and setters
now
defaults to FALSE instead of TRUE because DuckDBDataFrame does not
allow row names and must have withDimNames = FALSE

[spatialHeatmap](/packages/spatialHeatmap)
--------------

UPDATES

- Removed the rols package and developed custom internal functions for
  retrieving tissue ontology.
  
- Removed the scran
  and scuttle packages, and related functions are replaced with
  counterparts in bluster and scrapper respectively.

[Spectra](/packages/Spectra)
-------

                        Changes in version 1.23                         

Changes in 1.23.7

- Add filterNoisePeaks() function to support removing noise peaks
following the method from Dalla Valle et al. Anal. Chem. (2025).

Changes in 1.23.6

- Update citation information and add links to additional external
backends.

Changes in 1.23.5

- Fix reference unit test in the unit test suite for MsBackend.

Changes in 1.23.4

- Small fixes in the documentation.

Changes in 1.23.3

- Improve documentation for the use of SerialParam() vs SnowParam()
on
Windows.

Changes in 1.23.2

- Export the MsBackendHdf5Peaks class.

Changes in 1.23.1

- Add support for dataStorageBasePath() and dataStorageBasePath<- to
MsBackendHdf5Peaks.

[SpectralTAD](/packages/SpectralTAD)
-----------

                 Changes in version 1.29.2 (2026-07-19)                 

- Remove unavailable PRIMME dependency

[SpectraQL](/packages/SpectraQL)
---------

                         Changes in version 1.7                         

Changes in 1.7.1

- Fix unit test.

[SpectraStash](/packages/SpectraStash)
------------

- Initial Bioconductor Submission

[spicyR](/packages/spicyR)
------

                        Changes in version 2.1.1                        

- spicy(), cell method with soft cellWeights: the mean and variance of
  a label with itself (and of a group with
  its subtype) under random labelling now come by default from a closed
  form that corrects every pair of cells
  exactly for sampling without replacement (exact for up to two cells,
  second order in 1 / cells). Its cost no
  longer grows with cells x neighbourhood size^2, and each image's
  neighbour structure is built once per fit:
  all 784 weighted pairs of a 151-biopsy heart dataset take 14 s
  instead of 220 s. selfVariance = "exact" keeps
  the exact moments.

 
[SpiecEasi](/packages/SpiecEasi)
---------

                        Changes in version 2.1.1                        

DOCUMENTATION

- Set vignette YAML version fields to match the package version.

INFRASTRUCTURE

- Added venv to .Rbuildignore so local Python virtualenvs are not
included in package tarballs.

[splatter](/packages/splatter)
--------

                 Changes in version 1.38.0 (2026-11-25)                 

- Remove deprecated mfa simulation

[SpliceImpactR](/packages/SpliceImpactR)
-------------

                        Changes in version 1.1.5                        

Bug fixes

- The annotation cache-path test uses forward slashes on Windows,
matching the cache implementation and fixing R CMD check failures.

                        Changes in version 1.1.4                        

Bug fixes

- BioMart queries for release go to https://e<release>.ensembl.org,
which Ensembl forwards to that release's archive, instead of to the
archive host itself. With biomaRt before 2.70, new installations
could no longer connect to archive hosts ("Unable to contact any
Ensembl mirror"): biomaRt first downloads Ensembl's list of
archives, which is no longer available. Release 116 now works with
any biomaRt version. The data come from the same archives, so
results and cached features are unchanged.
- With biomaRt before 2.70, connection errors for an archive host
passed as ensembl_host explain this and suggest the e address.

                        Changes in version 1.1.3                        

Breaking changes and changed results

- get_protein_features() queries Ensembl 111 by default (was 109),
the
Ensembl release of the default GENCODE v45 annotation from
get_annotation(). Features from full runs change, and cached
features are fetched again once.
- ensembl_mirror is deprecated and ignored: Ensembl retired its
BioMart mirrors.

Bug fixes

- BioMart queries go to the Ensembl archive of release by its
explicit
host (releases 105 to 116). Ensembl 116 (June 2026) is the last
release with BioMart, and www.ensembl.org no longer serves it. The
new ensembl_host argument selects another host. Release 116 needs
biomaRt 2.70 or later (earlier versions send its queries to
www.ensembl.org), and connection errors name the host.

                        Changes in version 1.1.2                        

Breaking changes and changed results

- get_splicing_impact() selects transcript pairs with the ORF-aware
matcher (matching = "orf") by default; matching = "legacy" restores
the previous matcher. The ORF matcher checks each form's exact event
structure, never pairs a transcript with itself, labels pairs
without a structural match as approximate, and returns unresolved
comparisons in matching. Selected transcripts and pair counts
change, and matching is slower.

New features

- get_ranked_pairs(return_class = "S4") stores the selected pairs,
form rows and diagnostics in a SpliceImpactResult, where
get_splicing_impact() stores them, so the stepwise S4 workflow uses
the ORF matcher. The default still returns the list. get_pairs() S4
output sets metadata$matching = "legacy" and drops earlier ORF
diagnostics.

[See the full NEWS](/packages/release/bioc/news/SpliceImpactR/NEWS) or run `news(package = "SpliceImpactR")`.


[splicelogic](/packages/splicelogic)
-----------

                        Changes in version 1.2.0                        

- New event types: alternative transcription start and end sites
(find_atss(), find_ates()), spliced introns (find_si()), and
consecutive skipped and consecutive included exons (find_cse(),
find_cie()). The last two cover runs of two or more adjacent exons
skipped or included together, which find_se() / find_ie() do not
report: those report an exon only when it is skipped or included on
its own.
- New long-form aliases for every finder that spell out the event
type, e.g. find_skipped_exons() for find_se() or
find_alternative_transcription_start_sites() for find_atss().
- New input option: prepare_exons_by_partition() builds the input
from
two user-defined sets of transcripts (up and down). Unlike
prepare_exons(), it does not need DTU results.
- get_seq() extracts exon sequences against a genome.
- Plotting functions for event results: plot_transcripts(),
plot_transcripts_by_gene(), plot_transcripts_by_pair(),
plot_transcript_lfc() and event_palette().
- find_se(), find_ie(), find_mxe() and find_all_events() pass extra
arguments on to GenomicRanges::findOverlaps().
- Results carry event_type, event_tx_id, event_estimate and
event_<col> for any additional_columns.
- The default type is now "over"; "boundary" matches only the splice
site facing the candidate; and "in" is renamed "identical".

                        Changes in version 1.1.6                        

- Added find_cse() and find_cie() (long forms
find_consecutive_skipped_exons() and
find_consecutive_included_exons()) for runs of two or more adjacent
exons skipped or included together. Both are included in
find_all_events(), and event_palette() has a colour for each.
- Added mock data generators generate_ie(), generate_cse(),
generate_cie() and generate_si(), so every finder has its own
generator.
- Every generator now makes at most one event per gene: n_events is
the number of genes that get an event, and it warns and caps when
fewer genes have a candidate.

                        Changes in version 1.1.5                        

- Added plotting functions for event results: plot_transcripts(),
plot_transcripts_by_gene(), plot_transcripts_by_pair(),
plot_transcript_lfc() and event_palette() (require ggplot2).
- Added find_si() for spliced introns, the reverse of find_ri().
- The default type of find_se(), find_ie(), find_mxe() and
find_all_events() is now "over" (was "boundary").
- type = "boundary" now matches only the splice site of each flanking
exon that faces the candidate exon, on either strand.
- The "in" match type is renamed "identical".
- find_se(), find_ie(), find_mxe() and find_all_events() pass extra
arguments on to GenomicRanges::findOverlaps(), e.g. minoverlap = 5L.
- plyranges moved to Depends, so it is attached with splicelogic.

[See the full NEWS](/packages/release/bioc/news/splicelogic/NEWS) or run `news(package = "splicelogic")`.


[SpNeigh](/packages/SpNeigh)
-------

                        Changes in version 1.1.2                        

- Correct a typo in vignette

                        Changes in version 1.1.1                        

- Update README
- Update vignette

[SpotSweeper](/packages/SpotSweeper)
-----------

                        Changes in version 1.9.3                        

Bug fixes

- Fixed Windows parallel nearest-neighbor search in localOutliers()
and localVariance() (#27). The existing workers argument
automatically selects socket workers on Windows and multicore
workers on macOS/Linux; one worker continues to use serial
execution. No new arguments are needed.

                        Changes in version 1.9.2                        

Bug fixes

- Corrected the robust local z-score calculation in localOutliers().
Scores now use the focal spot value and the unscaled median absolute
deviation of neighboring spots, matching the published method.
Previous versions applied the normal-consistency scaling twice and
could associate a neighbor's score with the focal spot.
- Removed the spatialEco dependency, which is no longer needed for
local z-score calculation.

Seurat compatibility

- Reworked Seurat support around small, tested object adapters that
preserve spot names and ordering and use public SeuratObject
accessors.
- Added integration tests for localOutliers(), spatial-coordinate
extraction, metadata replacement, subsetting, and
flagVisiumOutliers() with real Seurat objects.
- Added explicit guidance when findArtifacts() is used with data
lacking mitochondrial signal, including Visium FFPE probe-panel
data.

[standR](/packages/standR)
------

                       Changes in version 1.17.6                        

- Update prepareSpatialDecon() to convert the selected log2 assay to
linear expression by default, preserving its normalization or batch
correction. Add isLog for linear inputs and recomputeCPM to
recompute linear CPM from raw counts with matching negative-probe
scaling. Raw metadata probes now require recomputation or raw-count
input rather than an independently normalized background.

                       Changes in version 1.17.4                        

- Fix installed-package checks for the DCC importer tests and qualify
stats::setNames() to avoid R CMD check notes.

                       Changes in version 1.17.3                        

- Add a standR-native readGeoMxFromDcc() importer for DCC/PKC
workflows and remove GeomxTools from Suggests.

                       Changes in version 1.17.2                        

- Add GeomxTools to Suggests so examples using
readGeoMxFromNanoStringGeoMxSet() can run in coverage workflows.

                       Changes in version 1.17.1                        

- Add readGeoMxFromNanoStringGeoMxSet() to convert
GeomxTools::readNanoStringGeoMxSet() / NanoStringGeoMxSet objects
from DCC/PKC workflows into SpatialExperiment objects for standR.

- Improve prepareSpatialDecon() so negative probes removed by
readGeoMx(rmNegProbe = TRUE) can be read from
metadata(spe)$NegProbes.

[statTarget](/packages/statTarget)
----------

                         Changes in version 2.0                         

NEW FEATURES

- New GUI
  - Mouse Hover for help information
  - .log file

- New Signal correction
  - Combat for QC-free Signal correction
  - QC-RFSC methods for metabolomics and proteomics data

- New feature slection
  - Random Forest and the Permutation based variable importance
    measures
  - new MDSplot for Random Forest
  - P-value based importance plot

- New data preprocessing
  - PQN/SUM/none normalization
  - center/none Scaling method

[strainspy](/packages/strainspy)
---------

- Initial Bioconductor Submission

[STRINGdb](/packages/STRINGdb)
--------

                       Changes in version 2.25.7                        

Bug fixes

- get_enrichment(), get_annotations() and get_geneset_description()
keep descriptions that contain #; they were cut there.
- get_homologs_besthits() uses every species of target_species_id; it
used only the first one. Species given as NA are left out.
- get_homology_graph(), and so remove_homologous_interactions(),
filter on the bitscore, as documented, instead of on the length of
the second protein. Their results change, also with the default
threshold.
- get_clusters() stops with the offered algorithms (fastgreedy,
walktrap and edge.betweenness) for any other value, such as
"spinglass", which its help page listed.
- With protocol = "http", the API methods work: all requests to
STRING's API use https, because STRING's redirect from http drops
the data of the request. The protocol applies only to the download
files.

Documentation

- ?get_png, ?plot_network, ?get_interactions, ?mp and the other
method
pages now open by name (except load: ?load stays base R's page).
- Corrected help pages: get_paralogs() returns the similarity among
the given proteins, not all paralogs; get_annotations(),
add_diff_exp_color() and get_ppi_enrichment() describe what they
really return; get_term_proteins() is marked as deprecated; the
class page describes every setting and every working method.
- Vignette: the clustering example uses the most significant genes
(it
used only up-regulated ones) and shows four clusters; the paralog
example compares the members of the p53 family; get_annotations()
leaves out only PubMed (KEGG terms are included); the figures fit
the page.

Planned for the next Bioconductor release (STRINGdb 3.0.0)

- inputGenes and preferredNames stay comma-separated text. The
planned
change to list columns is withdrawn, because write.csv() and similar
functions cannot write list columns.
- The caller_identity argument will stop with an error, and will be
removed in the release after.
- STRINGdb$new() without species will warn, and will stop in the
release after: give the species, e.g. species = 9606.
- When input_directory is not given, the downloaded files will be
kept
between R sessions, in the folder of tools::R_user_dir("STRINGdb",
"cache"). input_directory = "" keeps meaning a temporary folder.
- Without version, the newest STRING version known to the package
will
be used, without asking STRING; version = "latest" will ask STRING
for its current version.

[See the full NEWS](/packages/release/bioc/news/STRINGdb/NEWS) or run `news(package = "STRINGdb")`.


[STRprofilerR](/packages/STRprofilerR)
------------

- Initial Bioconductor Submission
 
[structToolbox](/packages/structToolbox)
-------------

                       Changes in version 1.25.2                        

- fix corr_coef p-values and coefficients for pearson

- fix nroot_transform ignoring the root argument

- fix knn_impute swapping row/column limits when imputing by samples

- fix ttest unpaired estimates labelled the wrong way round

- fix grid_search_1d optimum selection

- fix HSDEM keeping p-values for only one term

- fix bootstrap run() failing on a missing parameter

- fix linear_model model_predict erroring

- fix tSNE predicted output and init parameter

- fix log_transform, nroot_transform and mean_of_medians passing empty
  data when trained

- add drop_empty option to filter_na_count

- blank_filter warns when blank_label or qc_label match no samples

- fix equal_split only choosing the first rows of larger groups

- fix pairs_filter level check always passing

- fix PLSDA keeping training yhat and scores after predict

- fix SVM class_weights only working for a single weight

- fix empty test set when split training fraction rounds to 0

- fix mv_feature_filter_hist plotting wrong columns for within methods

- add tests for the above

                       Changes in version 1.25.1                        

- fix pqn_norm passing qc_frac as sample_frac

- fix pqn_norm ignoring ref_method

- fix pqn_norm predict not using the reference computed during training

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.44.0                        

NEW FEATURES

- Preserve names of combined metadata in combineRows/Cols (#98).
  Contributed by Aaron Lun <infinite.monkeys.with.keyboards@gmail.com>

[SummarizedHeatmap](/packages/SummarizedHeatmap)
-----------------

- Initial Bioconductor Submission

[SVP](/packages/SVP)
---

                        Changes in version 1.5.2                        

- add citation information. (2026-11-17, Tue)

                        Changes in version 1.5.1                        

- add p-value multiple correction for runGLOBALBV, runLISA and
runLOCALBV. (2026-08-18, Tue)

                        Changes in version 1.5.0                        

- Bioconductor 3.23 release, and Bioconductor 3.24 (devel) bump.
(2025-04-29, Wed)

[SwarnSeq](/packages/SwarnSeq)
--------

- Initial Bioconductor submission.


[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.25.2                        

- PhyloDistance now throws a warning if two trees with no labels in
common are compared.
- PhyloDistance now uses the Hungarian algorithm for Method="CI" so
that comparing a tree to itself correctly returns 0 distance.
- PhyloDistance now has a Exact argument, defaulting to TRUE to use
the Hungarian algorithm for Method="CI". Exact=FALSE will use the
previous greedy matching algorithm for backwards compatibility.

                       Changes in version 1.25.1                        

- Large rewrite of gff import, SummarizePairs and other functions,
added EvaluatePairs function and some other ancillary functions.

[SynMut](/packages/SynMut)
------

                 Changes in version 1.29.1 (2026-09-05)                 

- Bundled the component flowchart with the package and repaired its
  vignette image.

- Made the self-contained vignette independent of the remote flowchart
  URL.

- Updated README tutorial and development-install links and added
  maintainer checks.

[systemPipeShiny](/packages/systemPipeShiny)
---------------

                       Changes in version 1.23.1                        

Bug Fix

- Fix BiocCheck ERROR "Avoid references to external hosting
platforms": the About tab now reads the installed NEWS.md, and
module tabs use the bundled systemPipeR logo.
- Fix BiocCheck WARNINGs: the vignette no longer sets a global
eval=FALSE, and the R6 class man pages now have a Value section.
- Regenerate documentation with roxygen2 8.1.1.

[TADCompare](/packages/TADCompare)
----------

                 Changes in version 1.23.5 (2026-08-20)                 

- Remove residual PRIMME references

                 Changes in version 1.23.2 (2026-07-15)                 

- Fix z_thresh argument to default to 3 and used in the functions

                 Changes in version 1.23.1 (2026-07-15)                 

- Soften PRIMME dependency

[TaxSEA](/packages/TaxSEA)
------

                        Changes in version 1.5.8                        

New features

- ssTaxSEA() gains a reference that rare taxa cannot move. reference
"top_k" centers each sample on its top_k (default 20) most abundant
detected taxa instead of on every taxon in the table. With the
standard CLR, losing many rare taxa (after antibiotics or
chemotherapy, for example) lowers the reference and raises the score
of every surviving taxon, so hundreds of sets can increase together
with no change among the survivors. The top-k reference is unchanged
by losing rare taxa or by adding or removing rare and empty rows,
and rare taxa are still scored. A sample with fewer than top_k
detected taxa uses all of them, so zeros never enter the reference.
- ssTaxSEA() gains scale = "proportion", which divides each sample by
its total and applies pseudocount as a proportion (default 1e-6).
Counts, proportions and percentages such as MetaPhlAn output all
give the same result, and the pseudocount is the same proportion in
every sample whatever its depth.
- The default is unchanged (scale = "counts", reference = "all",
pseudocount 0.5), so existing results are reproduced exactly.

Bug fixes

- ssTaxSEA() now refuses percentages and other relative abundances on
the count scale, as it already did for proportions. Previously a
MetaPhlAn table in percent passed the check (its columns sum to 100,
not 1), and the pseudocount of 0.5, half a percent, swamped every
taxon rarer than that. Use scale = "proportion" for such tables.

                        Changes in version 1.5.7                        

Documentation

- Corrected the described output of TaxSEA() in the README and the
main vignette: it returns All_databases plus one data frame per
group of sources, each with 6 columns, not three data frames with 5.
- The README said BugSigDB was one of six databases used by default;
it has been optional since 1.5.4. It now also lists the
literature-collated sets (Valles-Colomer gut-brain modules, BloSSUM,
VANISH and others).
- Added ssTaxSEA() to the README, which did not mention single-sample
scoring at all, and corrected taxon_set_ranks() to taxon_rank_sets()
and data(TaxSEA_DB) to data(TaxSEA_db).
- The README now points at the guides on the package website, and the
website leads with them.

                        Changes in version 1.5.6                        

- Rewrote the package description: it now leads with what TaxSEA is
for and describes both approaches, TaxSEA() on ranked taxa and
ssTaxSEA() on a single sample. Also fixed two long-standing typos in
it ("alternations", and a missing "can be applied").

[See the full NEWS](/packages/release/bioc/news/TaxSEA/NEWS) or run `news(package = "TaxSEA")`.


[TDbasedUFE](/packages/TDbasedUFE)
----------

                       Changes in version 1.13.1                        

- Bundle the original 10,000-gene, six-sample count matrix.
- Use bundled data in QuickStart and tests.
- Remove unused readr, tximport and tximportData imports.

[thyroidBRS](/packages/thyroidBRS)
----------

- Initial Bioconductor Submission

[TiDEomics](/packages/TiDEomics)
---------

- Initial Bioconductor Submission


[tidyexposomics](/packages/tidyexposomics)
--------------

                        Changes in version 1.1.3                        

- Updated association functions to take formulas.
- Added functionality to increase range of plots.
- Increased code effeciency of enrichment analysis code.
- Updated sample/exposure cluster assignments to update the column
data and codebook, respectively.
- Updated multi-omics integration function to add factors to column
data.
- Added download functionality for single-cell example.

                        Changes in version 1.1.2                        

- Updated exposome score vignette to pull in suggested package.
- Updated error in MCIA factor naming conventions.
- Updated GO link to a live link to pull GO annotation data.
- Removing suggested package, mirt, from supplemental vignette.

                        Changes in version 1.1.1                        

- Added transform_sce_to_exposomicset() function to allow
SingleCellExperiment input.

[topdownr](/packages/topdownr)
--------

                        Changes in version 1.35                         

Changes in version 1.35.1

- Fix .groupByLabels for column-less input, e.g. Rle, that isn't
named
by as.data.frame anymore, see S4Vectors:79e0c50.
- Use the .data pronoun in ggplot2 calls to avoid the "no visible
binding for global variable" NOTE.
- Import IRanges::subject, it was resolved via the search path
before.
- Adapt roxygen comments to roxygen2 8.1.0 (hanging indent for
multi-line @importFrom, one line per @aliases, unambiguous links)
and re-document.

[toppgene](/packages/toppgene)
--------

                 Changes in version 1.1.2 (2026-07-31)                  

Bug fixes and improvements

- Fix dispatch of subset() for the CategoriesDataFrame object to the
appropriate superclass method. This addresses the failing unit test
in tests4classes.R` named "CategoriesDataFrame checks <- validity"
with the error message of:

Error in match.call(sys.function(which), sys.call(which),
expand.dots = FALSE, envir = sys.frame(parents[which])): invalid
'definition' argument`

- GitHub Actions CI fix of BiocCheck failing on
libgit2::git_repository_open_ext(). Add global config safe.directory
of the existing cloned package to be checked, identical to the
command run in "Checkout this git repository" CI step.

                 Changes in version 1.1.1 (2026-04-29)                  

Bug fixes and improvements

- Fix GitHub organization typo in URLs of the DESCRIPTION file.

- GitHub Actions CI fixes:

- Add back the "remotes" package that was dropped from the
bioconductor/bioconductor_docker:devel container dependencies.

- Putatively fix the test coverage sometimes hanging on
covr::package_coverage() by disabling parallel unit tests. The
rationale for this change was finding similar reports on covr
GitHub issues and the Posit community forum. The covr package
itself uses codecov for it's GitHub Actions and, therefore, does
not serve as useful reference CI code. Unable to reproduce this
fault with a local rootless container.

                 Changes in version 1.1.0 (2026-04-28)                  

Significant user-visible changes

- Increment development version following release of
Bioconductor 3.23.

                 Changes in version 1.0.2 (2026-08-03)                  

Bug fixes and improvements

- Backport fixes from 1.1.2 for CategoriesDataFrame dispatch of
subset() method (24c42a9) and for git repo ownership for BiocCheck()
during GitHub Actions (9659a0c).

                 Changes in version 1.0.1 (2026-04-30)                  

Bug fixes and improvements

- Backport fix from 1.1.1 of GitHub organization typo in URLs of the
DESCRIPTION file.

[TPP](/packages/TPP)
---

                       Changes in version 3.41.3                        

- Fixed R CMD check NOTES caused by undefined tidyverse/ggplot2
  symbols,
  and startup messages emitted from .onLoad, malformed roxygen/Rd
  markup.

- Fixed NEWS version headers so package NEWS can be parsed by R CMD
  check.

                       Changes in version 3.41.2                        

- Fixed package dependency checks by replacing methods::is() with
  inherits() and removing the unused mefa import.

                       Changes in version 3.41.1                        

- Removed the deprecated biobroom dependency and replaced
  biobroom::tidy.ExpressionSet() with an internal ExpressionSet
  tidying helper.

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.49.4                        

- Fix multiple warning in BiocCheck.

- Decrease the documentation size.

                       Changes in version 1.49.3                        

- Add break layout for lollipop plot.

                       Changes in version 1.49.2                        

- update documentation for track class.

                       Changes in version 1.49.1                        

- add help function getGeneModel.

[TraianProt](/packages/TraianProt)
----------

- Initial Bioconductor Submission

[treeio](/packages/treeio)
------

                       Changes in version 1.37.1                        

- read.beast() now reports a file that does not contain a tree
instead
of handing a text that is not Newick over to read.tree(), which
parsed it into a meaningless tree and could crash the R session,
e.g. the tabular output of MEGA (use read.mega_tabular() for it)
(2026-09-20, Sun)
- read.beast() is faster on a large file: the annotation of a node is
parsed by a function that no longer calls use_perl() and %>% for
every value and the tree it writes is not parsed a second time, a
file of 500 trees took 33s and now takes 25s (2026-09-20, Sun)
- read.phyloxml() is much faster and it is no longer limited to a
tree
nested less than 256 levels deep, the default limit of libxml2; the
clade tree is now walked with an explicit stack instead of a
recursion that bound all the nodes visited so far at every node, a
tree of 512 tips took 7s and now takes 0.4s (2026-09-20, Sun)
- edgeNum2nodeNum() is exported, it maps the edge_num of an
EPA/pplacer placement (a post-order traversal number) to the node
number of the reference tree, #31 (2026-09-20, Sun, #31)
- read.beast() now supports the UTREE keyword of an unrooted tree and
write.beast() no longer annotates a node without data (it wrote NULL
and could stop with object 'nl' not found), which is what happened
to the LSD2 timetree of IQ-TREE (2026-09-20, Sun, #111)
- read.mcmctree() now stores the 95% credibility interval of the node
age in a reltime_0.95_CI column (it used to be a column named 0.95)
so that it can be plotted with geom_range(range='reltime_0.95_CI',
center='reltime'), and the interval is kept as numbers (2026-09-20,
Sun, #13)
- the substitution is now reported with the node it cannot find a
sequence for instead of failing with seqA should have equal length
to seqB, and that message names the two sequences and their lengths
(2026-09-20, Sun, #91)
- read.iqtree() now says when the node labels hold a single support
value that cannot be split into SH-aLRT and UFBoot (e.g. the
standard bootstrap alone), and it no longer fails on a tree without
any branch support (2026-09-20, Sun, #114)
- read.iqtree() now reports that it cannot find a Newick tree in the
input instead of parsing the IQ-TREE report file (*.iqtree) into a
meaningless tree or crashing the R session (2026-09-20, Sun, #98)
- as.phylo() now uses the branch.length column by default so that the
branch lengths are not lost when a tree is re-rooted (2026-09-20,
Sun, #134)
- as.phylo() now uses the label column by default instead of the node
numbers, so that the tip and node labels are not lost when a tree is
re-rooted (2026-09-20, Sun, #120)
- read.beast() and read.mega() now support a TRANSLATE table with
non-consecutive keys (e.g. MEGA output); the tips are
numbered 1:Ntip and the node data is mapped accordingly (2026-09-20,
Sun, #132)

[See the full NEWS](/packages/release/bioc/news/treeio/NEWS) or run `news(package = "treeio")`.


[TSENAT](/packages/TSENAT)
------

- Initial Bioconductor Submission

[TSSr](/packages/TSSr)
----

- Initial Bioconductor Submission

[tTEscanR](/packages/tTEscanR)
--------

- Initial Bioconductor Submission

[txdbmaker](/packages/txdbmaker)
---------

                       Changes in version 1.10.0                        

NEW FEATURES

- Add 'kingdom' argument to makeTxDbFromEnsembl().

- makeTxDbFromUCSC("hg19", tablename="lincRNAsTranscripts") now is
  supported.

BUG FIXES

- Small fix to supportedUCSCtables().

- Get rid of 'makeTxDbFromUCSC()' spurious warning "genome version
  information is not available for this TxDb object".

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.32.0                        

- importData() now supports type="salmon" for mixed reference
transcript workflows, where salmon quantifies against an index that
combines annotated (e.g. GENCODE, Ensembl) and novel (e.g. de novo
assembled) transcripts. Because salmon does not record per-sub-index
digests itself, a companion Snakemake rule
(thelovelab/salmon-with-mixed-ref) computes them at index time, and
the new mixedDigest argument points importData() to the resulting
JSON file. This brings salmon users the same provenance tracking
previously available only for oarfish: the annotated portion of a
combined index is still recognized, so inspectDigests() and
updateMetadata() can attach transcript and gene metadata
automatically. See the "How it works" section of the README and the
"Mixed reference transcripts" section (subsection "salmon mixed
reference") of the vignette.
- New argument gencode_gtf_prefix in tximeta() to specify a
non-default GENCODE GTF, e.g. "primary_assembly" or "basic".
- Support for GENCODE 50 (H.s.), M39 (M.m.), and Ensembl 116.
- Bug fixes in updateMetadata() (prefer argument now respected,
missing character columns filled with NA) and in checkViaAHub().
- Tests and vignette updated for the reduced tximportData (>=
1.41.1).
- Added contributing guidelines (CONTRIBUTING.md) and updated
diagrams.

                       Changes in version 1.31.9                        

- Fix bug where updateMetadata() ignored its prefer argument, always
using the default order of registries.
- Fix bug in updateMetadata() where character columns added to
rowData
were filled with the string "NA" instead of NA for transcripts
without metadata.
- Tests and vignette updated for the reduced tximportData (>=
1.41.1),
which no longer contains the GENCODE v48 GTF.

                       Changes in version 1.31.6                        

- Fix bug in internal function checkViaAHub()

                       Changes in version 1.31.4                        

- GENCODE 50 (H.s.), M39 (M.m), and Ensembl 116 (June 2026). (except
for Drosophila which has a GTF issue)

                       Changes in version 1.31.3                        

- importData() now supports type="salmon" for mixed reference
transcript workflows. A new mixedDigest argument accepts the
filename (relative to each sample's quantification directory) of a
JSON file containing per-sub-index sha256 digests, produced by a
companion Snakemake rule using compute_fasta_digest. The digest
structure is compatible with the existing inspectDigests() and
updateMetadata() functions.

[See the full NEWS](/packages/release/bioc/news/tximeta/NEWS) or run `news(package = "tximeta")`.


[tximport](/packages/tximport)
--------

                       Changes in version 1.41.1                        

- Vignette, examples and tests updated for the smaller
  tximportData (>= 1.41.1), which now has 2 GEUVADIS samples
  and no longer includes inferential replicates (salmon Gibbs,
  kallisto bootstraps) or kallisto abundance.h5 files. These
  are archived at https://doi.org/10.5281/zenodo.22982575.
  Tests for inferential replicate and h5 import now write
  synthetic replicates to a temporary directory.

[UCell](/packages/UCell)
-----

                       Changes in version 2.17.1                        

- Accept any matrix-like inputs, includes DelayedArrays (credit
  to B. Demaille, see PR #55:
  <https://github.com/carmonalab/UCell/pull/55>)

- Add tests for main functions

- Default BPparam to SerialParam() for single-threaded runs

[UCSC.utils](/packages/UCSC.utils)
----------

                       Changes in version 1.10.0                        

- No significant changes in this version.

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.32.0                        

NEW FEATURES

- Added CWM (Contribution Weight Matrix) motifs for signed,
  unnormalised scores from TF-MoDISco / DeepLIFT-style attribution
  methods. Create with type = "CWM" and plot with use.type = "CWM".
  convert_type(cwm, "PPM") normalises absolute contributions per
  column, enabling downstream analyses; conversion from probabilistic
  types to CWM is not supported. MEME readers/writers accept
  CWM = TRUE, and generic matrix readers/writers accept type = "CWM".
  CWM input is opt-in, not auto-detected. See ?create_motif and
  ?convert_type.

- Added trim_cwm(), a CWM-aware trimmer that drops edge
  columns by absolute column sum. Defaults to the
  TF-MoDISco-lite fraction-of-peak rule (`trim.threshold =
  0.3`); an optional `abs.threshold` argument switches to an
  absolute-value cutoff. Walks inward from the chosen edges
  only (`trim.from = c("both", "left", "right")`),
  mirroring `trim_motifs()`'s edge-only semantics. Returns a
  trimmed CWM with the type tag preserved.

- Added compare_motifs_lite(), a DNA/RNA motif comparison function
  aligned with yamtk cmp (https://github.com/bjmt/yamtk). Uses
  per-column Pearson correlation and a TOMTOM-style empirical or
  parametric null, with Bonferroni and Benjamini-Hochberg adjustment.
  Returns a score/p-value/q-value matrix or a table of significant
  pairs. Score-only comparisons were ~40x faster than compare_motifs()
  on a 200x200 motif set at 4 threads. See ?compare_motifs_lite.

- Added scan_sequences_lite(), a minimalist DNA/RNA scanner using
  the same C++ core as scan_sequences(). It omits multifreq scoring,
  gapped motifs, q-values, exhaustive p-values, non-pvalue thresholds,
  respect.strand, and allow.nonfinite. On a 50-motif HOCOMOCOv11
  benchmark it was ~2-3x faster than scan_sequences(calc.pvals = TRUE)
  at 4 threads on equivalent work. See the Sequence Searches vignette
  for benchmarks, including comparisons with yamtk scan.

- Added view_motifs_lite() and motif_tree_lite(), DNA/RNA plotting
  alternatives using compare_motifs_lite(). Logos align to the
  highest-IC motif; trees use Pearson distances, (1 - score) / 2.
  Both reuse the existing plotting machinery. For amino-acid or
  custom-alphabet motifs, use view_motifs() and motif_tree().

- Added geom_motif(), annotate_motif(), geom_logo(), and
  annotate_logo(): ggplot2 layers that embed sequence logos into
  arbitrary plots. The annotate_*() helpers use fixed coordinates;
  geom_*() layers support mapped positions and different motifs per
  group. The *_motif variants accept universalmotif objects, while
  *_logo accepts row-named numeric matrices. Letter colours do not
  affect the plot's fill scale. See ?geom_motif.

- motif_tree_lite() now draws a sequence logo at each tip by default
  (tiplogo = TRUE), with aligned motif positions. The default layout
  changes from "circular" to "rectangular". Tip logos support only
  "rectangular" and "slanted" layouts and are unavailable for dist
  input. See ?motif_tree_lite for alignment and sizing options.

[See the full NEWS](/packages/release/bioc/news/universalmotif/NEWS) or run `news(package = "universalmotif")`.


[updateObject](/packages/updateObject)
------------

                       Changes in version 1.18.0                        

- No changes in this version.

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

[VISTA](/packages/VISTA)
-----

                        Changes in version 1.1.7                        

Internal

- The test that PCA and MDS plots do not depend on the order genes
are
supplied in no longer fails when an axis comes back negated.
Shuffling the input genes can make prcomp() return the mirror image
of the same embedding on some LAPACK/BLAS builds (it did on Linux
CI), which the test read as an order dependence. Axes are now
oriented with the same canonicalisation the 1.1.4 snapshots use, and
a new test checks that a sign flip passes while a moved point still
fails. No user-visible behaviour changes; only the test suite is
affected.

                        Changes in version 1.1.6                        

New features

- edgeR's quasi-likelihood F-test. create_vista(),
run_edger_analysis(), benchmark_vista_equivalence() and
validate_vista_deep() gain edger_test = c("LRT", "QL"). "QL" runs
glmQLFit()/glmQLFTest() (edgeR v4 pipeline), which the edgeR authors
recommend for bulk RNA-seq; the likelihood-ratio test treats the
dispersion as known and is liberal with few replicates. The test
used is recorded in cutoffs(v)$edger_test, and the report config
accepts edger_test:; for an edgeR object made before 1.1.6, which
does not record the test, the report states the likelihood-ratio
test and marks it as inferred. The equivalence check runs its own QL
reference. The vignettes now set edger_test.

Default changing in VISTA 1.4.0

- edger_test will default to "QL". It stays "LRT" in this release so
that no results change without warning. A call that runs edgeR
(method = "edger" or "both") without setting edger_test, or with
edger_test = NULL, prints a message, classed vista_default_change
and throttled like the deprecation warnings, saying so. Set
edger_test = "LRT" to keep current results after the change, or
edger_test = "QL" to use the new default now; either silences the
message. DESeq2 and limma-voom analyses are not affected. See
?VISTA-deprecated.

Bug fixes

- tximport abundance or lengths could be analysed as counts.
read_vista_counts(format = "tximport", counts_from = "abundance")
imported TPM, and counts_from = "length" transcript lengths, as if
they were counts, and edgeR and limma tested them without complaint.
Both now raise an error.

[See the full NEWS](/packages/release/bioc/news/VISTA/NEWS) or run `news(package = "VISTA")`.


[vsn](/packages/vsn)
---

                       Changes in version 3.81.1                        

- This package now includes a README and pkgdown website.

[WithinHostExperiment](/packages/WithinHostExperiment)
--------------------

- Initial Bioconductor Submission 

[woven](/packages/woven)
-----

- Initial Bioconductor Submission

[wSIR](/packages/wSIR)
----

- Initial Bioconductor Submission
  
[xcms](/packages/xcms)
----

                        Changes in version 4.11                         

Changes in version 4.11.5

- findChromPeaks() with verboseBetaColumns = TRUE calculates the peak
quality metrics considering also the retention times. Results will
thus be identical to those calculated with chromPeakSummary(). See
issue #849
- loadXcmsData() updates also $dataOrigin to ensure functions such as
fileNames() return the correct path/file names (issue #851).

Changes in version 4.11.4

- Add new vignette Plotting xcms results using lcmsPlot describing
how
the lcmsPlot package can be used to create ggplot2-based
visualizations of xcms preprocessing results (issue #844).

Changes in version 4.11.3

- Workaround for errors reading subsets of data from HDF5 files with
XcmsExperimentHdf5: original issue reported in rhdf5:
https://github.com/Huber-group-EMBL/rhdf5/issues/240 . Code can be
reverted once the problem is fixed upstream in rhdf5.

Changes in version 4.11.2

- Add support to return chromatographic data as a Chromatograms
object
with the chromatogram() method from MsExperiment, XcmsExperiment and
XcmsExperimentHdf5 objects.
- Add support for Chromatograms to the chromPeakChromatograms()
method
for XcmsExperiment and XcmsExperimentHdf5 objects.
- Add support for Chromatograms to the featureChromatograms() method
for XcmsExperiment and XcmsExperimentHdf5 objects.
- Add new methods featureChromPeaks() and featurePeakidx() to access
the mapping between features and chromatographic peaks.
- loadXcmsData() updates parameter class definitions in the object's
process history.
- Fix chromPeakData<- for XcmsExperiment not supporting a DataFrame
as
input (issue #839).

Changes in version 4.11.1

- refineChromPeaks,MergeNeighboringPeakParam: drop row names of
intermediate objects to avoid potentially misleading warning
message.

[ZarrArray](/packages/ZarrArray)
---------

                        Changes in version 1.2.0                        

NEW FEATURES

- Implement ZarrSparseMatrix objects. Analog to H5SparseMatrix and
  TENxMatrix objects from the HDF5Array package but for Zarr-based
  sparse matrices instead of HDF5-based.
  Co-authored-by: Artur-man <artur-man@hotmail.com>

- Add writeZarrSparseMatrix() for realizing a matrix-like object
  (typically sparse) as a ZarrSparseMatrix object.

- Add a new set of get/set_Zarr_realization_*() functions for
  controlling Zarr realization settings. They replace the
  old get/set_writeZarrArray_*() functions.

- Add support for bfloat and enable more Rarr test cases (#13).
  Authored by Hugo Gruson <git@hugogruson.fr>

SIGNIFICANT USER-VISIBLE CHANGES

- The default value for writeZarrArray()'s automatic arguments now
  is "auto" instead of NULL.

- writeZarrArray()'s argument 'zarr_version' now is an automatic
  argument.

DEPRECATED AND DEFUNCT

- Deprecate all the get/set_writeZarrArray_*() functions in favor of
  the new get/set_Zarr_realization_*() functions.

[zellkonverter](/packages/zellkonverter)
-------------

                       Changes in version 1.24.0                        

New features

- Add an environment for *anndata* v0.13.4. This is now the
  default.

NEWS from existing Data Experiment Packages
===================================


[AnnotatedBCGEData](/packages/AnnotatedBCGEData)
-----------------

- Initial Bioconductor submission.

[aracne.networks](/packages/aracne.networks)
---------------

                       Changes in version 1.39.1                        

- The networks are no longer bundled in the package. They are hosted
  on Zenodo (https://doi.org/10.5281/zenodo.22918956) and downloaded
  on demand, keeping the package size well below the 100 MB limit of
  the Bioconductor build system. The network files are identical to
  the ones distributed in previous versions.

- New function getRegulon() downloads a network, caches it locally
  with BiocFileCache and returns the regulon object. Replace
  data(regulonblca) with regulonblca <- getRegulon("blca").

- New function listRegulons() lists the available networks.

- The documentation of the individual networks is merged into
  ?regulons.

- The vignette is converted to R Markdown.

- Added unit tests.

[bugphyzz](/packages/bugphyzz)
--------

                        Changes in version 1.8.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- importBugphyzz() defaults to data release v1.2.0
  (10.5281/zenodo.24166279): 477,741 annotations for 93,900 taxa, 36
  attributes and 19 sources.

- The data are licensed CC BY-SA 4.0 from release v1.2.0, because they
  include annotations derived from Wikipedia; earlier releases are CC
  BY 4.0. The package stays Artistic-2.0.

- Imputed annotations are kept or dropped per attribute, by the
  lower 95% bound of a distance-matched Cohen's kappa, or R2 for
  numeric attributes (#292).

- Discrete signature names start with bugphyzz:, not bugphyz: (#256).

- excludeRarely and v are deprecated; use asrScoreThreshold,
  kappaThreshold and r2Threshold (#271).

- Levi Waldron is the maintainer (#257).

NEW FEATURES

- physiologies() also reads the scripted sources in bugphyzzExports
  sources/: the 2026 Bergey's Manual extraction, Madin et al. 2020,
  BacDive, The Microbe Directory and NCBI RefSeq assemblies (#275).

- New attribute gc content and new shape value vibrio.

- New source Wikipedia_2026 (confidence medium, evidence nas): traits
  stated in English Wikipedia introductions, with the new attributes
  halophily class, temperature class and ph class; skin is also a
  habitat (#296).

- attribute_synonyms.tsv maps source terms to the controlled
  vocabulary (#264).

- New arguments kappaThreshold, r2Threshold and asrScoreThreshold; a
  threshold below its default also loads the imputations the export
  left out (#292).

- Validation scores per attribute and per state ship as
  validation_attributes.tsv and validation_states.tsv.

BUG FIXES AND MINOR IMPROVEMENTS

- Cell length and width no longer collapse to the lower bound of their
  range (#275).

- Downloads are cached per version.

- Typos and missing terms in attributes.tsv are fixed (#281).

- Tests that download data skip when offline (#266).

[CLAMPData](/packages/CLAMPData)
---------

- Initial Bioconductor submission

[curatedBladderData](/packages/curatedBladderData)
------------------

                       Changes in version 1.50.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The 12 datasets are no longer stored inside the package. They are
  hosted on Zenodo and downloaded individually on first use, then
  cached locally with BiocFileCache. This shrinks the installed
  package from about 75 MB to a few MB.

- New exported function curatedBladderData(): call it with no
  arguments to list the available datasets, with one dataset name to
  get an ExpressionSet, or with several names to get a named list.
  test = TRUE loads small offline subsets bundled with the package,
  used by examples, tests and the vignette.

- Legacy data(GSE89_eset) access still works. It downloads through the
  same cache, deferred until the object is first used, but is
  deprecated and will be removed in a future release.

- inst/extdata/createEsetList.R now loads datasets through the getter.
  A new test.mode option in the patientselection config selects the
  offline subsets.

- The cache is package-specific, at
  tools::R_user_dir("curatedBladderData", "cache"). Downloads are
  verified against md5 checksums recorded in
  inst/extdata/zenodo-manifest.csv.

BUG FIXES

- biocViews listed OvarianCancerData; removed, and GEO added.

- The package-level help page described ovarian cancer, listed the
  curatedOvarianData author list, and named a different maintainer
  than DESCRIPTION. Corrected to match DESCRIPTION.

- Replaced the unused affy dependency with Biobase, which is where
  ExpressionSet, exprs() and pData() actually come from.

- createEsetList.R no longer errors when the filters exclude every
  dataset (which() on an empty list).

INTERNAL

- Added testthat tests and inst/scripts/ documenting how the Zenodo
  files, the offline fixtures and the data stubs are generated.

[curatedBreastData](/packages/curatedBreastData)
-----------------

                       Changes in version 2.41.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Package data has been successfully offloaded to Zenodo
  (https://zenodo.org/records/21415886) and is no longer directly
  embedded in the package. The data is now dynamically retrieved and
  cached via BiocFileCache.

NEW FEATURES

- Added new retrieval functions getCuratedBreastDataExprSetList() and
  getClinicalData() to access the remote dataset from Zenodo.

- Embedded a miniature test dataset archive in inst/extdata for rapid
  local testing and offline execution without remote dependencies.

- Added curatedBreastData_pData (pre-concatenated phenotype metadata
  from all datasets) to the package reference index, manual vignette,
  and documentation.

- Configured a new pkgdown documentation site to host the package
  reference manual and manual vignette.

BUG FIXES

- Fixed operator precedence bugs and replaced deprecated 1:... loop
  sequence generators with type-safe seq_along() and seq_len().

- Converted sapply calls to type-safe vapply() in collapseDupProbes().

- Updated DESCRIPTION metadata to include the latest R dependency
  (>= 4.6.0), missing biocViews, URL, and BugReports fields.

INTERNAL CHANGES

- Migrated all package documentation to roxygen2 and removed legacy
  .Rd manual pages.

- Converted the package vignette from Sweave (.Rnw) to R Markdown
  (.Rmd) for improved readability and compilation performance.

- Added README.md and _pkgdown.yml mirroring the design of
  MultiAssayExperiment.

- Configured GitHub Actions CI workflows for automated R CMD check and
  Bioconductor compliance validation (BiocCheck).

[curatedCRCData](/packages/curatedCRCData)
--------------

                       Changes in version 2.46.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The 34 datasets are no longer stored inside the package: they are
  hosted on Zenodo and downloaded individually on first use, then
  cached locally with BiocFileCache. This shrinks the installed
  package from ~400 MB to a few MB.

- New exported function curatedCRCData(): call it with no arguments to
  list available datasets, with one dataset name to get an
  ExpressionSet, or with several names to get a named list. test =
  TRUE loads small offline subsets bundled with the package (used by
  examples, tests, and the vignette).

- Legacy data(GSE39582_eset) access still works (it downloads through
  the same cache, deferred until the object is first used) but is
  deprecated and will be removed in a future release.

- inst/extdata/createEsetList.R now loads datasets through the getter;
  a new test.mode option in the patientselection config files selects
  the offline subsets.

- The cache is package-specific (tools::R_user_dir("curatedCRCData",
  "cache")); delete that directory to reclaim disk space. Downloads
  are verified against md5 checksums recorded in
  inst/extdata/zenodo-manifest.csv.

INTERNAL

- inst/scripts/data-integrity-report.Rmd now verifies size_bytes as
  well as md5. R/getData.R compares size_bytes against the cached file
  on every cache hit, so a wrong size in the manifest would have
  silently evicted and re-downloaded a valid file on every call.
  Caught by Copilot review on lima1/curatedBladderData#3.

- Replaced the never-executed RUnit scaffolding with testthat tests.

- CI restored to R CMD check + BiocCheck + pkgdown via the waldronlab
  reusable workflow (previously reduced to pkgdown-only because
  checking 400 MB of data exhausted runner memory).

- Removed external_data_store.txt and Namespace: auto; explicit
  NAMESPACE exporting only curatedCRCData().

DOCUMENTATION AND METADATA

- Repaired UTF-8 transcoding damage ("??" artifacts) in 32 dataset man
  pages: author names, journal references, and statistical symbols
  restored (ASCII equivalents).

- biocViews extended to ColonCancerData, MicroarrayData, RNASeqData,
  and GEO for better discoverability.

- Removed the empty FULLVcuratedCRCData_counts.csv placeholder.

[curatedMetagenomicData](/packages/curatedMetagenomicData)
----------------------

                       Changes in version 3.22.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- mergeData() and returnSamples() return a sparse assay (dgCMatrix)
  for gene_families, as curatedMetagenomicData() already does for a
  single study, instead of a dense matrix. For 48 samples from two
  studies the assay is 0.95 GB rather than 2.5 GB. Code that needs a
  dense matrix can call as.matrix() on the assay (#304, @amsegura;
  #318, @adamsorbie)

- rownames = "NCBI" gives bare NCBI Taxonomy IDs as row names, without
  the species: and genus: prefixes. Rows that share an ID are summed,
  with a message naming them. For "short" and "NCBI", clades that NCBI
  resolves no further than genus are dropped with a warning; rownames
  = "long" keeps every row (#316, @LiNk-NY)

BUG FIXES AND MINOR IMPROVEMENTS

- returnSamples() no longer mislabels columns. When a requested sample
  had no assay match and sampleMetadata was not in study order, the
  returned columns carried the names and colData of other samples
  (#354, @lwaldron)

- returnSamples() no longer drops a sample whose sample_id also occurs
  in another requested study (#354, @lwaldron)

- mergeData() returns a SummarizedExperiment input unchanged, instead
  of failing with subscript is out of bounds when given the result of
  returnSamples() (#321, @BinhongLiu)

- returnSamples() keeps only the requested samples of each study
  before merging, and mergeData() no longer builds intermediate data
  frames. Peak memory for four gene_families samples from two studies
  fell from 10.3 GB to 3.6 GB (#304, @amsegura)

- mergeData() and returnSamples() stop with a clear error when the
  merged gene_families assay would hold more than the 2^31 - 1
  non-zero values that a sparse matrix can store, instead of failing
  with p[length(p)] cannot exceed 2^31-1. The 5,474 samples from 45
  studies requested in #304 are over this limit (#304, @amsegura)

- rownames = "NCBI" keeps its row tree with mia 1.21.30 and later.
  Rows that share an NCBI Taxonomy ID are summed in this package
  instead of with mia::agglomerateByVariable(), which now removes the
  tree (#356, @lwaldron)

- The curatedMetagenomicData() help page and the vignette explain how
  the "long", "short" and "NCBI" row names differ (#342, @lwaldron)

- The vignette builds when packages in Suggests are not installed, and
  OmicsMLRepoR is removed from Suggests while it cannot be installed
  on Bioconductor 3.24 (#339, #340, @lwaldron)

- methods and Matrix are declared in Imports (#353, #355, @lwaldron)

- Unit tests updated for changes in mia (#314, @TuomasBorman)

- Added CITATION.cff (@LiNk-NY)

- Continuous integration uses the Bioconductor reusable GitHub
  workflow, with a coverage report (#338, #344, @lwaldron)

[curatedOvarianData](/packages/curatedOvarianData)
------------------

                       Changes in version 1.52.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The 30 datasets are no longer stored inside the package: they are
  hosted on Zenodo and downloaded individually on first use, then
  cached locally with BiocFileCache. This shrinks the installed
  package from ~400 MB to a few MB.

- New exported function curatedOvarianData(): call it with no
  arguments to list available datasets, with one dataset name to get
  an ExpressionSet, or with several names to get a named list. test =
  TRUE loads small offline subsets bundled with the package (used by
  examples, tests, and the vignette).

- Legacy data(TCGA_eset) access still works (it downloads through the
  same cache, deferred until the object is first used) but is
  deprecated and will be removed in a future release.

- inst/extdata/createEsetList.R now loads datasets through the getter;
  a new test.mode option in the patientselection config files selects
  the offline subsets.

- The cache is package-specific
  (tools::R_user_dir("curatedOvarianData", "cache")); delete that
  directory to reclaim disk space. Downloads are verified against md5
  checksums recorded in inst/extdata/zenodo-manifest.csv.

INTERNAL

- inst/scripts/data-integrity-report.Rmd now verifies size_bytes as
  well as md5. R/getData.R compares size_bytes against the cached file
  on every cache hit, so a wrong size in the manifest would have
  silently evicted and re-downloaded a valid file on every call.
  Caught by Copilot review on lima1/curatedBladderData#3.

- Replaced the never-executed RUnit scaffolding with testthat tests.

- Explicit NAMESPACE exporting only curatedOvarianData().

- Replaced the dead DFCI URL with the GitHub repository URL and added
  BugReports; removed the committed build/ PDF.

DOCUMENTATION AND METADATA

- Repaired UTF-8 transcoding damage ("??" artifacts) in 7 dataset man
  pages (author names, copyright marks, and statistical symbols).

- biocViews extended to ExpressionData, OvarianCancerData, miRNAData,
  ArrayExpress, and GEO for better discoverability.

[DaparToolshedData](/packages/DaparToolshedData)
-----------------

- Initial Bioconductor Submission

[davidTiling](/packages/davidTiling)
-----------

                       Changes in version 1.53.0                        

- The pre-built ExpressionSet object, previously available via
  data(davidTiling) has been removed, alongside the bundled CEL files
  in inst/celfiles. To get the data from now on, please use the new
  getDavidTilingData() function, which will download the raw CEL files
  from ArrayExpress (or use a local cache if they were downloaded
  previously), and convert them to the ExpressionSet object via
  tilingArray::readCel2eSet().

[FlowSorted.Blood.EPIC](/packages/FlowSorted.Blood.EPIC)
---------------------

                       Changes in version 2.17.2                        

- Added regression tests protecting adult-blood IDOL estimates.

- Hardened ExperimentHub reference retrieval and validation.

- Added direct custom RGChannelSet and RGChannelSetExtended support
  while
  retaining character-name compatibility.

- Avoided attaching reference-data packages during estimation.

- Moved minfi and ExperimentHub from Depends to Imports.

- Removed obsolete memory.limit usage and the unreachable nlme pathway.

- Declared MatrixGenerics::rowRanges explicitly.

- Removed unused probe statistics from the IDOL workflow.

- Improved vignette chunk labeling and execution.

                       Changes in version 2.17.1                        

- FlowSorted.CordBlood.450k was removed from Suggests because it was
  deprecated by Bioconductor.

- compositeCellType = "CordBlood" now gives an informative error
  directing
  users to FlowSorted.CordBloodCombined.450k.

- The maintainer ORCID and package metadata were updated.

- Added the missing colData<- namespace import.

[gDRtestData](/packages/gDRtestData)
-----------

               Changes in version 2026-11-17 (2026-11-17)               

- add the sotorasib-resistant NCI-H23 and NCI-H358 cell lines to the
  cell line annotations

               Changes in version 2026-10-24 (2026-10-24)               

- update the reference datasets to the fit types and p-values of
  gDRutils 1.11.22

               Changes in version 2026-10-20 (2026-10-20)               

- update the combination reference datasets to the current gDRcore,
  which restores CIScore and the single-agent arm fit

               Changes in version 2026-10-08 (2026-10-08)               

- update authors data

               Changes in version 2026-10-01 (2026-10-01)               

- add drug and cell line annotations for the gCSI Sotorasib workshop
  dataset

- use the Gnumber as the drug name for undisclosed compounds

               Changes in version 2026-09-17 (2026-09-17)               

- update reference datasets so that constant fits carry the corrected
  xc50 sign

               Changes in version 2026-07-07 (2026-07-07)               

- update documentation for PRISM dataset

               Changes in version 2026-06-22 (2026-06-22)               

- add PRISM public dataset

               Changes in version 2026-06-09 (2026-06-09)               

- add workshop drug and cell line annotations

               Changes in version 2026-06-02 (2026-06-02)               

- add PRISM public dataset

               Changes in version 2026-06-01 (2026-06-01)               


[geneLenDataBase](/packages/geneLenDataBase)
---------------

                       Changes in version 1.50.0                        

- The provided data objects have been resaved with a more efficient
  compression ("xz"), enabling a 10-15% reduction of the size taken up
  by the whole package - this change should not present any visible
  change to the end-user, as the objects are identical once loaded

[GSE142512](/packages/GSE142512)
---------

- Initial Bioconductor Submission

[GSE280465](/packages/GSE280465)
---------

- Initial Bioconductor Submission

[hapmapsnp6](/packages/hapmapsnp6)
----------

                       Changes in version 1.55.1                        

- The data is no longer distributed inside the package. It is deposited
  on Zenodo (doi:10.5281/zenodo.23990454) and served through
  ExperimentHub, which brings the source tarball down from 264 MB to a
  few kilobytes.

- USER VISIBLE CHANGE: data(crlmmResult) no longer works, because the
  package ships no data directory. Use crlmmResult() instead, which
  returns the same SnpSet. The first call downloads the object and
  stores it in the local ExperimentHub cache; later calls read from
  that cache.

  ```
  ## before
  data(crlmmResult)
  ## now
  crlmmResult <- crlmmResult()
  ```

- USER VISIBLE CHANGE: the three CEL files are now obtained through
  celFiles(), which returns their paths on disk. Code that used
  system.file("celFiles", package = "hapmapsnp6") should use
  celFiles(asDir = TRUE), which returns a directory holding the files.

  ```
  ## before
  path <- system.file("celFiles", package = "hapmapsnp6")
  cels <- list.celfiles(path, full.names = TRUE)
  ## now
  cels <- celFiles()
  ```

- crlmmResult(metadata = TRUE) and celFiles(metadata = TRUE) return the
  resource metadata without downloading anything.

- Both accessors pass further arguments on to ExperimentHub(), so that
  an alternate hub location can be used.

- The maintainer address was updated, and authorship now uses the
  Authors@R field.

[HuMMANet](/packages/HuMMANet)
--------

- Initial Bioconductor Submission

[macrophage](/packages/macrophage)
----------
                       Changes in version 1.29.1                        

- Reduced package size. The source tarball is now under 100 MB 
  The following files were removed. All of them are archived, with the full
  Salmon output for all 24 samples, on Zenodo:  The `with-inf-reps` branch
  on GitHub is a frozen snapshot of the package that still includes them.

  - The Gibbs inferential replicates
    (`inst/extdata/quants/*/aux_info/bootstrap/`, ~393 MB).
    `tximeta()` / `tximport()` on the files in `inst/extdata/quants` no
    longer imports inferential replicates.
  - `aux_info/ambig_info.tsv.gz` for each sample.
  - The full Gencode v29 GTF, `inst/extdata/gencode.v29.annotation.gtf.gz`.

- new files and datatsets

  - `inst/extdata/gencode.v29.annotation.slim.gtf.gz`: a slim Gencode v29
    GTF with only gene and transcript lines and 6 attributes (5 MB instead
    of 38 MB). It has no exon lines, so `tximeta::addExons()` will not
    work with it. It has a different filename from the full GTF on purpose:
    tximeta caches TxDbs by GTF filename, so this avoids clashing with a
    TxDb built from the full Gencode v29 GTF. Code that called
    `makeLinkedTxome(gtf=...)` with the full GTF path should switch to the
    slim filename. Built by `inst/scripts/makeSlimGTF.sh`.
  - `data(macro_txp_se)`: a transcript-level `SummarizedExperiment` with
    20 inferential replicates, for all transcripts of 2000 random expressed
    genes plus GBP1-GBP7. The `counts` and `infRep` assays are rounded to
    1 decimal. This is now used by the fishpond `swish` vignette.
  - `data(macro_tx2gene)`: a transcript-to-gene table for all Gencode v29
    transcripts, for use with
    `tximeta::summarizeToGene(se, skipRanges=TRUE, tx2gene=macro_tx2gene)`.
  - Both datasets are built by `inst/scripts/makeSumExp.R`. The script
    needs a macrophage install that still has the inferential replicates
    (<= 1.29.0, or the `with-inf-reps` branch).


[marinerData](/packages/marinerData)
-----------

                       Changes in version 1.13.1                        

NEW FEATURES

- Added new (m)cool files:

- LEUK_HEK_PJA27_inter_30.mcool

- LEUK_HEK_PJA30_inter_30.mcool

- LEUK_HEK_PJA27_inter_30_2500Kb.cool

- LEUK_HEK_PJA30_inter_30_2500Kb.cool


[msdata](/packages/msdata)
------

                       Changes in version 0.53.1                        

- removed sciex 20171016_POOL_POS_3_105-134.mzML and
  20171016_POOL_POS_1_105-134.mzML, now available in MsDataHub 1.13.2.

- removed MS3TMT11.mzML and MS3TMT10_01022016_32917.33481.mzML.gz,
  available in MsDataHub.

[oct4](/packages/oct4)
----

                       Changes in version 1.29.1                        

- Removed the Gibbs inferential replicates (aux_info/bootstrap/) from
  all 12 samples to reduce the package size to under 100 MB.
  num_bootstraps in each aux_info/meta_info.json was set to 0 so that
  tximport/tximeta do not look for the missing files.

[PscanRBackgrounds](/packages/PscanRBackgrounds)
-----------------
     
- Initial Bioconductor submission

[qPLEXdata](/packages/qPLEXdata)
---------

                       Changes in version 1.31.2                        

- Converted the vignette from Sweave/PDF (qPLEXdata.Rnw) to R
  Markdown/HTML
  (qPLEXdata.Rmd) using BiocStyle, removing the LaTeX build dependency.

                       Changes in version 1.31.1                        

- Fixed vignette build error: pass the contrast name (names(contrasts))
  rather than the named contrast vector to getContrastResults(),
  matching
  the current qPLEXanalyzer API.

[RforProteomics](/packages/RforProteomics)
--------------

                       Changes in version 1.51.1                        

- Remove defunct readMzTab reading. Replace by
  readSummarizedExperiment() and convert to MSnSet.

- Update startup message and point to RforMS.

[sesameData](/packages/sesameData)
----------

                       Changes in version 1.31.3                        

- The eight objects added in 1.31.1 now come from ExperimentHub, which
  assigned their ids on 2026-10-20: HM450/EPIC/EPICv2/MSA.address
  (EH10648-EH10651) and KYCG.<platform>.Mask.20261008
  (EH10652-EH10655).
  sesameDataGet() no longer falls back to Zenodo for them, and
  sesameDataCache() caches them with the rest.

                       Changes in version 1.31.2                        

- Documentation fixes from a review of the vignette. The setup chunk
  also caches genomeInfo.hg38, Mammal40.address and probeIDSignature,
  which the probe and annotation examples need. A sentence about a
  decoy = TRUE option in a
  function that does not exist is removed.

- metadata.csv describes EPIC.5.SigDF.normal as its five normal SigDFs
  (PrEC.Rep1, PrEC.Rep2, NAF.A, NAF.B, NAF.C), not an LNCaP
  segmentation. The ExperimentHub record is unchanged.

- The Zenodo fallback in sesameDataGet() downloads the file first, with
  a 600 s timeout. On a slow link (about 80 kB/s from HPC compute
  nodes)
  R's default 60 s cut the download short, so titles without an EHID
  (e.g. HM450.address, MSA.address) failed with "cannot read from
  connection".

                       Changes in version 1.31.1                        

- Added KYCG.{HM450,EPIC,EPICv2,MSA}.Mask.20261008, built from the
  InfiniumAnnotation v8.2 masks, so sesame's quality mask reads the
  same
  M_* tracks as SeSAMe2 on all four human arrays (sesame >= 1.31.6,
  zhou-lab/sesame#1).

- Rebuilt HM450.address, EPIC.address, EPICv2.address and MSA.address
  from the InfiniumAnnotation v8.3 orderings, so sesame reads the same
  probe design table as SeSAMe2. With these and the masks above, R's
  default QCDPB betas equal SeSAMe2's to ~1e-15 on all four arrays.

[SingleCellMultiModal](/packages/SingleCellMultiModal)
--------------------

                       Changes in version 1.26.0                        

Bug fixes and minor improvements

- Replaced local pull request workflow with reusable GitHub Actions
  wrapper

- Updated author, contributor, and funding metadata in DESCRIPTION and
  CITATION.cff

- Trimmed unnecessary software references from CITATION.cff

[SpaMTPData](/packages/SpaMTPData)
----------

- Initial Bioconductor Submission

[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.25.5                        

BUG FIXES

- An incorrect link for an object pulled by fetch_data() was
  corrected.

                       Changes in version 1.25.4                        

BUG FIXES

- An incorrect link for an object pulled by fetch_data() was
  corrected.

                       Changes in version 1.25.3                        

NEW FEATURES

- spatialLIBD now uses scrapper internally for pseudobulking the data.
  This was implemented by @manishabarse at as part of
  https://github.com/LieberInstitute/spatialLIBD/pull/123.

- fetch_data() was updated to release the habenulaAtlas project data.
  Implemented by @Nick-Eagles.

                       Changes in version 1.25.2                        

BUG FIXES

- Fixed app_server() to use the default_cluster information with a
  _colors suffix on the sce_layer object for choosing colors for the
  spatial registration and gene set enrichment tabs on the layer-level
  component of run_app(). This was done live during a LIBD RStats club
  session with input from @manishabarse, @lahuuki, and @Nick-Eagles.

                       Changes in version 1.25.1                        

NEW FEATURES

- fetch_data() now has example Xenium data you can access with type =
  "spe_xenium_example".

- vis_gene(), vis_gene_grid(), vis_clus(), vis_clus_grid() all now
  support Xenium data. This is partially done through the new internal
  functions vis_gene_c() and vis_clus_c().

- run_app() now supports Xenium data. So shiny-powered apps can now be
  built for Xenium (10x Genomics) data.

- These features were added by @lahuuki and @lcolladotor as part of
  https://github.com/LieberInstitute/spatialLIBD/pull/121.

[systemPipeRdata](/packages/systemPipeRdata)
---------------

                       Changes in version 2.17.1                        

- Large example data resources, including FASTQ files, annotation
  resources, and precomputed BAM files, are no longer distributed
  directly with the systemPipeRdata package. These resources are now
  hosted on Zenodo and downloaded and cached locally using
  BiocFileCache.

- Added getSPRdata() for retrieving the external example data
  resources. FASTQ and annotation resources are retrieved by default,
  while precomputed BAM files can optionally be retrieved with
  getSPRdata(bam = TRUE).

- Updated genWorkenvir() and genWorkdata() to retrieve the required
  example data through getSPRdata(). Their existing workflow setup
  behavior is otherwise unchanged.
  NEW FEATURES

- Command-line interface uses now Common Workflow Language (CWL)

- Introduction of new SYSargsList control class for managing workflows

- Seamless processing of R code sections as workflow steps

- Added technical reports to scientific reporting infrastructure

- Auto-generation of CWL files from command-line strings

- Plotting of workflow topology graphs

[tximportData](/packages/tximportData)
------------

REDUCED PACKAGE SIZE

- The package's example data has been cut from ~593 MB on disk (~388 MB
  compressed) to ~171 MB on disk (~96 MB compressed). The built tarball is
  now ~96 MiB.

  The full data as it was before this change (commit `02887f7`) is kept on
  the `with-inf-reps` branch and tag at
  https://github.com/mikelove/tximportData and archived on Zenodo:
  https://doi.org/10.5281/zenodo.22982575

  - Removed

    - `salmon_gibbs/`: Salmon 0.8.1 output with Gibbs samples (Ensembl v87,
      6 samples)
    - `kallisto_boot/`: kallisto 0.43.0 output with bootstraps (Ensembl v87,
      6 samples)
    - `tx2gene.ensembl.v87.csv`
    - `gencode/gencode.v48.annotation.gtf.gz`. tximeta can now link
      transcripts to genes with a DataFrame, so it no longer needs a GTF.
    - `alevin/{mouse1_unst_50,mouse1_LPS2_50}/alevin/bfh.txt` and
      `raw_cb_frequency.txt`
    - `kallisto/*/abundance.h5` (no bootstraps). `abundance.tsv.gz` is kept.
    - GEUVADIS samples `ERR188297`, `ERR188329`, `ERR188288` and `ERR188356`
      from `kallisto/`, `rsem/`, `salmon/` and `sailfish/`. `samples.txt` and
      `samples_extended.txt` now have 2 rows: `ERR188088` and `ERR188021`.

  - Kept

    - `salmon_ec/*/aux_info/eq_classes.txt` (used by fishpond)
    - `alevin/mouse1_unst_50_boot/` (used by eds)
    - `cufflinks/`, still with all 6 samples

  - Changes needed in downstream code

    - Code that reads `samples.txt` now gets 2 samples instead of 6. For example,
      replace `paste0("sample", 1:6)` with `paste0("sample", seq_along(files))`,
      and replace `rep(c("A","B"), each=3)` with `factor(c("A","B"))`.
    - For example data with inferential replicates (Gibbs samples or
      bootstraps), use the Zenodo record above or another data package such as
      `macrophage`.

Deprecated and Defunct Packages
===============================

**SOFTWARE:**

Forty six software packages were removed from this release (after being deprecated
in Bioc 3.23):

- APAlyzer, barcodetrackR, basecallQC, biobroom, biodbChebi, BiRewire, BPRMeth, BubbleTree, ccrepe, CelliD, cummeRbund, CuratedAtlasQueryR, debCAM, DeconRNASeq, geneXtendeR, GEOexplorer, GNET2, granulator, hca, hmdbQuery, IMAS, IONiseR, Melissa, MetaNeighbor, MethReg, mfa, microSTASIS, MineICA, motifcounter, MSPrep, nearBynding, netprioR, Organism.dplyr, partCNV, receptLoss, RgnTX, RiboProfiling, RTCGA, scviR, SGCP, SigFuge, soGGi, spatzie, SQLDataFrame, supersigs, tLOH

Previously announced deprecated packages were fixed and restored to Bioconductor 3.24: 

- ballgown, bamsignals, ChIPQC, hypeR, normr, RcisTarget, rols, shiny.gosling

Sixteen software packages are deprecated in this release and will be removed in Bioc 3.25:

- ChIPXpress, CSAR, GSCA, iModMix, InterCellar, MMDiff2, MoonlightR, netZooR, pathRender, RMassBank, scAnnotatR, sevenbridges, SwathXtend, systemPipeTools, tidytof, TSAR

**EXPERIMENT DATA:** 

Two experimental data packages were removed from this release (after being
deprecated in BioC 3.23):

- Fletcher2013a, Fletcher2013b

Previously announced deprecated package was fixed and restored to Bioconductor 3.24: 

- curatedBreastData


Thirty experimental data packages are deprecated in this release and will be
removed in Bioc 3.25:

- Affyhgu133aExpr, Affyhgu133Plus2Expr, Affymoe4302Expr, ccdata, ChIPXpressData, ConnectivityMap, CytoMethIC, FlowSorted.CordBlood.450k, furrowSeg, GeuvadisTranscriptExpr, iModMixData, ListerEtAlBSseq, mammaPrintData, MMDiffBamSubset, PasillaTranscriptExpr, pd.atdschip.tiling, RcisTarget.hg19.motifDBs.cisbpOnly.500bp, RMassBankData, RTCGA.clinical, RTCGA.CNV, RTCGA.methylation, RTCGA.miRNASeq, RTCGA.mRNA, RTCGA.mutations, RTCGA.PANCAN12, RTCGA.rnaseq, RTCGA.RPPA, Single.mTEC.Transcriptomes, SVM2CRMdata, tofsimsData

**ANNOTATION DATA:** 

Four annotation annotation packages were removed from this release (after being deprecated
in Bioc 3.23).

- hpAnnot, humanCHRLOC, mouseCHRLOC, ratCHRLOC 

The following out-dated annotation packages were deprecated in 3.24 and will be removed in
Bioc 3.25:

- TxDb.Mmusculus.UCSC.mm10.ensGene, TxDb.Athaliana.BioMart.plantsmart22, TxDb.Athaliana.BioMart.plantsmart25, TxDb.Cfamiliaris.UCSC.canFam3.refGene, TxDb.Cfamiliaris.UCSC.canFam4.refGene, TxDb.Ggallus.UCSC.galGal4.refGene, TxDb.Mmulatta.UCSC.rheMac3.refGene, TxDb.Ptroglodytes.UCSC.panTro4.refGene, TxDb.Rnorvegicus.UCSC.rn4.ensGene, TxDb.Hsapiens.BioMart.igis, TxDb.Rnorvegicus.BioMart.igis, TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts

**WORKFLOWS:** 

No workflow packages were removed from this release.

No workflow packages have been deprecated in this release.

**BOOKS:**

No books were removed from this release.

No books were deprecated in this release.
