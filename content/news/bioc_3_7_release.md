May 1, 2018

Bioconductors:

We are pleased to announce Bioconductor 3.7, consisting of 1560
software packages, 342 experiment data packages, 919 annotation
packages, and 21 workflows.

There are 98 new software packages, 16 new data experiment packages, 
2 new workflows, and many updates and improvements
to existing packages; Bioconductor 3.7 is compatible with R 3.5.0,
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

Getting Started with Bioconductor 3.7
======================================

To update to or install Bioconductor 3.7:

1. Install R >=3.5.0.  Bioconductor 3.7 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 98 new software packages in this release of Bioconductor.

- [adaptest](https://bioconductor.org/packages/adaptest)
  Data-adaptive test statistics represent a general methodology for
  performing multiple hypothesis testing on effects sizes while
  maintaining honest statistical inference when operating in
  high-dimensional settings. The utilities provided here
  extend the use of this general methodology to many common data
  analytic challenges that arise in modern computational and genomic
  biology.

- [ASICS](https://bioconductor.org/packages/ASICS) With a set of pure
  metabolite reference spectra, ASICS quantifies concentration of
  metabolites in a complex spectrum. The identification of
  metabolites is performed by fitting a mixture model to the spectra
  of the library with a sparse penalty. The method and its
  statistical properties are described in Tardivel et al. (2017)
  <doi:10.1007/s11306-017-1244-5>.

- [bcSeq](https://bioconductor.org/packages/bcSeq) This Rcpp-based
  package implements a highly efficient data structure and algorithm
  for performing alignment of short reads from CRISPR or shRNA
  screens to reference barcode library. Sequencing error are
  considered and matching qualities are evaluated based on Phred
  scores. A Bayes' classifier is employed to predict the originating
  barcode of a read. The package supports provision of user-defined
  probability models for evaluating matching qualities. The package
  also supports multi-threading.

- [BEARscc](https://bioconductor.org/packages/BEARscc) BEARscc is a
  noise estimation and injection tool that is designed to assess
  putative single-cell RNA-seq clusters in the context of
  experimental noise estimated by ERCC spike-in controls.

- [BiFET](https://bioconductor.org/packages/BiFET) BiFET identifies
  TFs whose footprints are over-represented in target regions
  compared to background regions after correcting for the bias
  arising from the imbalance in read counts and GC contents between
  the target and background regions. For a given TF k, BiFET tests
  the null hypothesis that the target regions have the same
  probability of having footprints for the TF k as the background
  regions while correcting for the read count and GC content bias.
  For this, we use the number of target regions with footprints for
  TF k, t_k as a test statistic and calculate the p-value as the
  probability of observing t_k or more target regions with footprints
  under the null hypothesis.

- [BiocOncoTK](https://bioconductor.org/packages/BiocOncoTK) Provide
  a central interface to various tools for genome-scale analysis of
  cancer studies.

- [BioNetStat](https://bioconductor.org/packages/BioNetStat) A
  package to perform differential network analysis, differential node
  analysis (differential coexpression analysis), network and
  metabolic pathways view.

- [CAGEfightR](https://bioconductor.org/packages/CAGEfightR) CAGE is
  a widely used high throughput assay for measuring transcription
  start site (TSS) activity. CAGEfightR is an R/Bioconductor package
  for performing a wide range of common CAGE data analysis tasks.
  Core functionality includes: import of CAGE TSSs (CTSSs), tag (or
  unidirectional) clustering for TSS identification, bidirectional
  clustering for enhancer identification, annotation with transcript
  and gene models, calculation of TSS shapes and quantification of
  CAGE expression as expression matrices.

- [ccfindR](https://bioconductor.org/packages/ccfindR) A collection
  of tools for cancer single cell RNA-seq analysis. Cell clustering
  and feature gene selection analysis employ maximum likelihood and
  Bayesian non-negative matrix factorization algorithm. Input data
  set consists of RNA count matrix, gene, and cell bar code
  annotations.  Analysis outputs are factor matrices for multiple
  ranks, quality measures (maximum likelihood) or evidence (Bayesian)
  with respect to rank.  The package includes utilities for
  downstream analyses, including meta-gene identification,
  visualization, and construction of rank-based trees for cell
  clusters.

- [CellScore](https://bioconductor.org/packages/CellScore) The
  CellScore package contains functions to evaluate the cell identity
  of a test sample, given a cell transition defined with a starting
  (donor) cell type and a desired target cell type. The evaluation is
  based upon a scoring system, which uses a set of standard samples
  of known cell types, as the reference set. The functions have been
  carried out on a large set of microarray data from one platform
  (Affymetrix Human Genome U133 Plus 2.0). In principle, the method
  could be applied to any expression dataset, provided that there are
  a sufficient number of standard samples and that the data are
  normalized.

- [CHARGE](https://bioconductor.org/packages/CHARGE) Identifies
  genomic duplications or deletions from gene expression data.

- [ChIC](https://bioconductor.org/packages/ChIC) Quality control
  pipeline for ChIP-seq data using a comprehensive set of quality
  control metrics, including previously proposed metrics as well as
  novel ones, based on local characteristics of the enrichment
  profile. The framework allows assessing quality of samples with
  sharp or broad enrichment profiles, whereas previously proposed
  metrics were not taking this into account. CHIC provides a
  reference compendium of quality control metrics and trained machine
  learning models for scoring samples.

- [ChIPSeqSpike](https://bioconductor.org/packages/ChIPSeqSpike)
  Chromatin Immuno-Precipitation followed by Sequencing (ChIP-Seq) is
  used to determine the binding sites of any protein of interest,
  such as transcription factors or histones with or without a
  specific modification, at a genome scale. The many steps of the
  protocol can introduce biases that make ChIP-Seq more qualitative
  than quantitative. For instance, it was shown that global histone
  modification differences are not caught by traditional downstream
  data normalization techniques. A case study reported no differences
  in histone H3 lysine-27 trimethyl (H3K27me3) upon Ezh2 inhibitor
  treatment. To tackle this problem, external spike-in control were
  used to keep track of technical biases between conditions.
  Exogenous DNA from a different non-closely related species was
  inserted during the protocol to infer scaling factors that enabled
  an accurate normalization, thus revealing the inhibitor effect.
  ChIPSeqSpike offers tools for ChIP-Seq spike-in normalization.
  Ready to use scaled bigwig files and scaling factors values are
  obtained as output. ChIPSeqSpike also provides tools for ChIP-Seq
  spike-in assessment and analysis through a versatile collection of
  graphical functions.

- [CTDquerier](https://bioconductor.org/packages/CTDquerier) Package
  to retrieve and visualize data from the Comparative Toxicogenomics
  Database (http://ctdbase.org/). The downloaded data is formated as
  DataFrames for further downstream analyses.

- [CytoDx](https://bioconductor.org/packages/CytoDx) This package
  provides functions that predict clinical outcomes using single cell
  data (such as flow cytometry data, RNA single cell sequencing data)
  without the requirement of cell gating or clustering.

- [ddPCRclust](https://bioconductor.org/packages/ddPCRclust) The
  ddPCRclust algorithm can automatically quantify the CPDs of
  non-orthogonal ddPCR reactions with up to four targets. In order to
  determine the correct droplet count for each target, it is crucial
  to both identify all clusters and label them correctly based on
  their position. For more information on what data can be analyzed
  and how a template needs to be formatted, please check the
  vignette.

- [DEComplexDisease](https://bioconductor.org/packages/DEComplexDisease)
  It is designed to find the differential expressed genes (DEGs) for
  complex disease, which is characterized by the heterogeneous
  genomic expression profiles. Different from the established DEG
  analysis tools, it does not assume the patients of complex diseases
  to share the common DEGs. By applying a bi-clustering algorithm,
  DECD finds the DEGs shared by as many patients. In this way, DECD
  describes the DEGs of complex disease in a novel syntax, e.g. a
  gene list composed of 200 genes are differentially expressed in 30%
  percent of studied complex disease. Applying the DECD analysis
  results, users are possible to find the patients affected by the
  same mechanism based on the shared signatures.

- [decontam](https://bioconductor.org/packages/decontam) Simple
  statistical identification of contaminating sequence features in
  marker-gene or metagenomics data. Works on any kind of feature
  derived from environmental sequencing data (e.g. ASVs, OTUs,
  taxonomic groups, MAGs,...). Requires DNA quantitation data or
  sequenced negative control samples.

- [DEScan2](https://bioconductor.org/packages/DEScan2) Integrated
  peak and differential caller, specifically designed for broad
  epigenomic signals.

- [DEsingle](https://bioconductor.org/packages/DEsingle) DEsingle is
  an R package for differential expression (DE) analysis of
  single-cell RNA-seq (scRNA-seq) data. It defines and detects 3
  types of differentially expressed genes between two groups of
  single cells, with regard to different expression status (DEs),
  differential expression abundance (DEa), and general differential
  expression (DEg). DEsingle employs Zero-Inflated Negative Binomial
  model to estimate the proportion of real and dropout zeros and to
  define and detect the 3 types of DE genes. Results showed that
  DEsingle outperforms existing methods for scRNA-seq DE analysis,
  and can reveal different types of DE genes that are enriched in
  different biological functions.

- [diffcoexp](https://bioconductor.org/packages/diffcoexp) A tool for
  the identification of differentially coexpressed links (DCLs) and
  differentially coexpressed genes (DCGs). DCLs are gene pairs with
  significantly different correlation coefficients under two
  conditions. DCGs are genes with significantly more DCLs than by
  chance.

- [diffcyt](https://bioconductor.org/packages/diffcyt) Statistical
  methods for differential discovery in high-dimensional cytometry
  (including flow cytometry, mass cytometry or CyTOF, and
  oligonucleotide-tagged cytometry) using high-resolution clustering
  and moderated tests.

- [dmrseq](https://bioconductor.org/packages/dmrseq) This package
  implements an approach for scanning the genome to detect and
  perform accurate inference on differentially methylated regions
  from Whole Genome Bisulfite Sequencing data. The method is based on
  comparing detected regions to a pooled null distribution, that can
  be implemented even when as few as two samples per population are
  available. Region-level statistics are obtained by fitting a
  generalized least squares (GLS) regression model with a nested
  autoregressive correlated error structure for the effect of
  interest on transformed methylation proportions.

- [DominoEffect](https://bioconductor.org/packages/DominoEffect) The
  functions support identification and annotation of hotspot residues
  in proteins. These are individual amino acids that accumulate
  mutations at a much higher rate than their surrounding regions.

- [drawProteins](https://bioconductor.org/packages/drawProteins) This
  package draws protein schematics from Uniprot API output. From the
  JSON returned by the GET command, it creates a dataframe from the
  Uniprot Features API. This dataframe can then be used by geoms
  based on ggplot2 and base R to draw protein schematics.

- [DropletUtils](https://bioconductor.org/packages/DropletUtils)
  Provides a number of utility functions for handling single-cell
  (RNA-seq) data from droplet technologies such as 10X Genomics. This
  includes data loading, identification of cells from empty droplets,
  removal of barcode-swapped pseudo-cells, and downsampling of the
  count matrix.

- [enrichplot](https://bioconductor.org/packages/enrichplot) The
  'enrichplot' package implements several visualization methods for
  interpreting functional enrichment results obtained from ORA or
  GSEA analysis. All the visualization methods are developed based on
  'ggplot2' graphics.

- [FELLA](https://bioconductor.org/packages/FELLA) Enrichment of
  metabolomics data using KEGG entries. Given a set of affected
  compounds, FELLA suggests affected reactions, enzymes, modules and
  pathways using label propagation in a knowledge model network. The
  resulting subnetwork can be visualised and exported.

- [GARS](https://bioconductor.org/packages/GARS) Feature selection
  aims to identify and remove redundant, irrelevant and noisy
  variables from high-dimensional datasets. Selecting informative
  features affects the subsequent classification and regression
  analyses by improving their overall performances. Several methods
  have been proposed to perform feature selection: most of them
  relies on univariate statistics, correlation, entropy measurements
  or the usage of backward/forward regressions. Herein, we propose an
  efficient, robust and fast method that adopts stochastic
  optimization approaches for high-dimensional. GARS is an innovative
  implementation of a genetic algorithm that selects robust features
  in high-dimensional and challenging datasets.

- [GateFinder](https://bioconductor.org/packages/GateFinder) Given a
  vector of cluster memberships for a cell population, identifies a
  sequence of gates (polygon filters on 2D scatter plots) for
  isolation of that cell type.

- [GDCRNATools](https://bioconductor.org/packages/GDCRNATools) This
  is an easy-to-use package for downloading, organizing, and
  integrative analyzing RNA expression data in GDC with an emphasis
  on deciphering the lncRNA-mRNA related ceRNA regulatory network in
  cancer. Three databases of lncRNA-miRNA interactions including
  spongeScan, starBase, and miRcode, as well as three databases of
  mRNA-miRNA interactions including miRTarBase, starBase, and miRcode
  are incorporated into the package for ceRNAs network construction.
  limma, edgeR, and DESeq2 can be used to identify differentially
  expressed genes/miRNAs. Functional enrichment analyses including
  GO, KEGG, and DO can be performed based on the clusterProfiler and
  DO packages. Both univariate CoxPH and KM survival analyses of
  multiple genes can be implemented in the package. Besides some
  routine visualization functions such as volcano plot, bar plot, and
  KM plot, a few simply shiny apps are developed to facilitate
  visualization of results on a local webpage.

- [GDSArray](https://bioconductor.org/packages/GDSArray) GDS files
  are widely used to represent genotyping or sequence data. The
  GDSArray package implements the `GDSArray` class to represent nodes
  in GDS files in a matrix-like representation that allows easy
  manipulation (e.g., subsetting, mathematical transformation) in
  _R_. The data remains on disk until needed, so that very large
  files can be processed.

- [GeneStructureTools](https://bioconductor.org/packages/GeneStructureTools)
  GeneStructureTools can be used to create in silico alternative
  splicing events, and analyse potential effects this has on
  functional gene products.

- [gep2pep](https://bioconductor.org/packages/gep2pep) Pathway
  Expression Profiles (PEPs) are based on the expression of pathways
  (defined as sets of genes) as opposed to individual genes. This
  package converts gene expression profiles to PEPs and performs
  enrichment analysis of both pathways and experimental conditions,
  such as "drug set enrichment analysis" and "gene2drug" drug
  discovery analysis respectively.

- [GOfuncR](https://bioconductor.org/packages/GOfuncR) GOfuncR
  performs a gene ontology enrichment analysis based on the ontology
  enrichment software FUNC. GO-annotations are obtained from
  OrganismDb or OrgDb packages ('Homo.sapiens' by default); the
  GO-graph is included in the package and updated regularly
  (10-Apr-2018). GOfuncR provides the standard candidate vs.
  background enrichment analysis using the hypergeometric test, as
  well as three additional tests: (i) the Wilcoxon rank-sum test that
  is used when genes are ranked, (ii) a binomial test that is used
  when genes are associated with two counts and (iii) a Chi-square or
  Fisher's exact test that is used in cases when genes are associated
  with four counts. To correct for multiple testing and
  interdependency of the tests, family-wise error rates are computed
  based on random permutations of the gene-associated variables.
  GOfuncR also provides tools for exploring the ontology graph and
  the annotations, and options to take gene-length or spatial
  clustering of genes into account. From version 0.99.14 on it is
  also possible to provide custom annotations and ontologies.

- [GSEABenchmarkeR](https://bioconductor.org/packages/GSEABenchmarkeR)
  The GSEABenchmarkeR package implements an extendable framework for
  reproducible evaluation of set- and network-based methods for
  enrichment analysis of gene expression data. This includes support
  for the efficient execution of these methods on comprehensive real
  data compendia (microarray and RNA-seq) using parallel computation
  on standard workstations and institutional computer grids. Methods
  can then be assessed with respect to runtime, statistical
  significance, and relevance of the results for the phenotypes
  investigated.

- [gsean](https://bioconductor.org/packages/gsean) Biological
  molecules in a living organism seldom work individually. They
  usually interact each other in a cooperative way. Biological
  process is too complicated to understand without considering such
  interactions. Thus, network-based procedures can be seen as
  powerful methods for studying complex process. However, many
  methods are devised for analyzing individual genes. It is said that
  techniques based on biological networks such as gene co-expression
  are more precise ways to represent information than those using
  lists of genes only. This package is aimed to integrate the gene
  expression and biological network. A biological network is
  constructed from gene expression data and it is used for Gene Set
  Enrichment Analysis.

- [hipathia](https://bioconductor.org/packages/hipathia) Hipathia is
  a method for the computation of signal transduction along signaling
  pathways from transcriptomic data. The method is based on an
  iterative algorithm which is able to compute the signal intensity
  passing through the nodes of a network by taking into account the
  level of expression of each gene and the intensity of the signal
  arriving to it. It also provides a new approach to functional
  analysis allowing to compute the signal arriving to the functions
  annotated to each pathway.

- [hmdbQuery](https://bioconductor.org/packages/hmdbQuery) Define
  utilities for exploration of human metabolome database, including
  functions to retrieve specific metabolite entries and data
  snapshots with pairwise associations
  (metabolite-gene,-protein,-disease).

- [iCNV](https://bioconductor.org/packages/iCNV) Integrative copy
  number variation (CNV) detection from multiple platform and
  experimental design.

- [igvR](https://bioconductor.org/packages/igvR) Access to igv.js,
  the Integrative Genomics Viewer running in a web browser.

- [IMMAN](https://bioconductor.org/packages/IMMAN) Reconstructing
  Interlog Protein Network (IPN) integrated from several Protein
  protein Interaction Networks (PPINs). Using this package,
  overlaying different PPINs to mine conserved common networks
  between diverse species will be applicable.

- [InTAD](https://bioconductor.org/packages/InTAD) The package is
  focused on the detection of correlation between expressed genes and
  selected epigenomic signals i.e. enhancers obtained from ChIP-seq
  data within topologically associated domains (TADs). Various
  parameters can be controlled to investigate the influence of
  external factors and visualization plots are available for each
  analysis step.

- [iSEE](https://bioconductor.org/packages/iSEE) Provides functions
  for creating an interactive Shiny-based graphical user interface
  for exploring data stored in SummarizedExperiment objects,
  including row- and column-level metadata. Particular attention is
  given to single-cell data in a SingleCellExperiment object with
  visualization of dimensionality reduction results.

- [iteremoval](https://bioconductor.org/packages/iteremoval) The
  package provides a flexible algorithm to screen features of two
  distinct groups in consideration of overfitting and overall
  performance. It was originally tailored for methylation locus
  screening of NGS data, and it can also be used as a generic method
  for feature selection. Each step of the algorithm provides a
  default method for simple implemention, and the method can be
  replaced by a user defined function.

- [kissDE](https://bioconductor.org/packages/kissDE) Retrieves
  condition-specific variants in RNA-seq data (SNVs,
  alternative-splicings, indels). It has been developed as a
  post-treatment of 'KisSplice' but can also be used with user's own
  data.

- [LineagePulse](https://bioconductor.org/packages/LineagePulse)
  LineagePulse is a differential expression and expression model
  fitting package tailored to single-cell RNA-seq data (scRNA-seq).
  LineagePulse accounts for batch effects, drop-out and variable
  sequencing depth. One can use LineagePulse to perform longitudinal
  differential expression analysis across pseudotime as a continuous
  coordinate or between discrete groups of cells (e.g. pre-defined
  clusters or experimental conditions). Expression model fits can be
  directly extracted from LineagePulse.

- [loci2path](https://bioconductor.org/packages/loci2path) loci2path
  performs statistics-rigorous enrichment analysis of eQTLs in
  genomic regions of interest. Using eQTL collections provided by the
  Genotype-Tissue Expression (GTEx) project and pathway collections
  from MSigDB.

- [MACPET](https://bioconductor.org/packages/MACPET) The MACPET
  package can be used for binding site analysis for ChIA-PET data.
  MACPET reads ChIA-PET data in BAM or SAM format and separates the
  data into Self-ligated, Intra- and Inter-chromosomal PETs.
  Furthermore, MACPET breaks the genome into regions and applies 2D
  mixture models for identifying candidate peaks/binding sites using
  skewed generalized students-t distributions (SGT). It then uses a
  local poisson model for finding significant binding sites. MACPET
  is mainly written in C++, and it supports the BiocParallel package.

- [MAGeCKFlute](https://bioconductor.org/packages/MAGeCKFlute)
  MAGeCKFlute is designed to surporting downstream analysis,
  utilizing the gene summary data provided through MAGeCK or
  MAGeCK-VISPR. Quality control, normalization, and screen hit
  identification for CRISPR screen data are performed in pipeline.
  Identified hits within the pipeline are categorized based on
  experimental design, and are subsequently interpreted by functional
  enrichment analysis.

- [martini](https://bioconductor.org/packages/martini) martini deals
  with the low power inherent to GWAS studies by using prior
  knowledge represented as a network. SNPs are the vertices of the
  network, and the edges represent biological relationships between
  them (genomic adjacency, belonging to the same gene, physical
  interaction between protein products). The network is scanned using
  SConES, which looks for groups of SNPs maximally associated with
  the phenotype, that form a close subnetwork.

- [mCSEA](https://bioconductor.org/packages/mCSEA) Identification of
  diferentially methylated regions (DMRs) in predefined regions
  (promoters, CpG islands...) from the human genome using Illumina's
  450K or EPIC microarray data. Provides methods to rank CpG probes
  based on linear models and includes plotting functions.

- [mdp](https://bioconductor.org/packages/mdp) The Molecular Degree
  of Perturbation webtool quantifies the heterogeneity of samples. It
  takes a data.frame of omic data that contains at least two classes
  (control and test) and assigns a score to all samples based on how
  perturbed they are compared to the controls. It is based on the
  Molecular Distance to Health (Pankla et al. 2009), and expands on
  this algorithm by adding the options to calculate the z-score using
  the modified z-score (using median absolute deviation), change the
  z-score zeroing threshold, and look at genes that are most
  perturbed in the test versus control classes.

- [MDTS](https://bioconductor.org/packages/MDTS) A package for the
  detection of de novo copy number deletions in targeted sequencing
  of trios with high sensitivity and positive predictive value.

- [MetaNeighbor](https://bioconductor.org/packages/MetaNeighbor)
  MetaNeighbor allows users to quantify cell type replicability
  across datasets using neighbor voting.

- [missRows](https://bioconductor.org/packages/missRows) The missRows
  package implements the MI-MFA method to deal with missing
  individuals ('biological units') in multi-omics data integration.
  The MI-MFA method generates multiple imputed datasets from a
  Multiple Factor Analysis model, then the yield results are combined
  in a single consensus solution. The package provides functions for
  estimating coordinates of individuals and variables, imputing
  missing individuals, and various diagnostic plots to inspect the
  pattern of missingness and visualize the uncertainty due to missing
  values.

- [MSstatsQCgui](https://bioconductor.org/packages/MSstatsQCgui)
  MSstatsQCgui is a Shiny app which provides longitudinal system
  suitability monitoring and quality control tools for proteomic
  experiments.

- [netSmooth](https://bioconductor.org/packages/netSmooth) netSmooth
  is an R package for network smoothing of single cell RNA sequencing
  data. Using bio networks such as protein-protein interactions as
  priors for gene co-expression, netsmooth improves cell type
  identification from noisy, sparse scRNAseq data.

- [OmaDB](https://bioconductor.org/packages/OmaDB) A package for the
  orthology prediction data download from OMA database.

- [omicplotR](https://bioconductor.org/packages/omicplotR) A Shiny
  app for visual exploration of omic datasets as compositions, and
  differential abundance analysis using ALDEx2. Useful for exploring
  RNA-seq, meta-RNA-seq, 16s rRNA gene sequencing with visualizations
  such as principal component analysis biplots (coloured using
  metadata for visualizing each variable), dendrograms and stacked
  bar plots, and effect plots (ALDEx2). Input is a table of counts
  and metadata file (if metadata exists), with options to filter data
  by count or by metadata to remove low counts, or to visualize
  select samples according to selected metadata.

- [omicsPrint](https://bioconductor.org/packages/omicsPrint)
  omicsPrint provides functionality for cross omic genetic
  fingerprinting, for example, to verify sample relationships between
  multiple omics data types, i.e. genomic, transcriptomic and
  epigenetic (DNA methylation).

- [ORFik](https://bioconductor.org/packages/ORFik) Tools for
  manipulation of RiboSeq, RNASeq and CageSeq data. ORFik is
  extremely fast through use of C, data.table and GenomicRanges.
  Package allows to reassign starts of the transcripts with the use
  of CageSeq data, automatic shifting of RiboSeq reads, finding of
  Open Reading Frames for the whole genomes and many more.

- [perturbatr](https://bioconductor.org/packages/perturbatr)
  perturbatr does stage-wise analysis of large-scale genetic
  perturbation screens for integrated data sets consisting of
  multiple screens. For multiple integrated perturbation screens a
  hierarchical model that considers the variance between different
  biological conditions is fitted. The resulting list of gene effects
  is then further extended using a network propagation algorithm to
  correct for false negatives.

- [phantasus](https://bioconductor.org/packages/phantasus) Phantasus
  is a web-application for visual and interactive gene expression
  analysis. Phantasus is based on Morpheus – a web-based software for
  heatmap visualisation and analysis, which was integrated with an R
  environment via OpenCPU API. Aside from basic visualization and
  filtering methods, R-based methods such as k-means clustering,
  principal component analysis or differential expression analysis
  with limma package are supported.

- [plyranges](https://bioconductor.org/packages/plyranges) A
  dplyr-like interface for interacting with the common Bioconductor
  classes Ranges and GenomicRanges. By providing a grammatical and
  consistent way of manipulating these classes their accessiblity for
  new Bioconductor users is hopefully increased.

- [pogos](https://bioconductor.org/packages/pogos) Provide simple
  utilities for querying bhklab PharmacoDB, modeling API outputs, and
  integrating to cell and compound ontologies.

- [PowerExplorer](https://bioconductor.org/packages/PowerExplorer)
  Estimate and predict power among groups and multiple sample sizes
  with simulated data, the simulations are operated based on
  distribution parameters estimated from the provided input dataset.

- [powerTCR](https://bioconductor.org/packages/powerTCR) This package
  provides a model for the clone size distribution of the TCR
  repertoire. Further, it permits comparative analysis of TCR
  repertoire libraries based on theoretical model fits.

- [RandomWalkRestartMH](https://bioconductor.org/packages/RandomWalkRestartMH)
  This package performs Random Walk with Restart on multiplex and
  heterogeneous networks. It is described in the following article:
  "Random Walk With Restart On Multiplex And Heterogeneous Biological
  Networks". https://www.biorxiv.org/content/early/2017/08/30/134734


- [RcisTarget](https://bioconductor.org/packages/RcisTarget)
  RcisTarget identifies transcription factor binding motifs (TFBS)
  over-represented on a gene list. In a first step, RcisTarget
  selects DNA motifs that are significantly over-represented in the
  surroundings of the transcription start site (TSS) of the genes in
  the gene-set. This is achieved by using a database that contains
  genome-wide cross-species rankings for each motif. The motifs that
  are then annotated to TFs and those that have a high Normalized
  Enrichment Score (NES) are retained. Finally, for each motif and
  gene-set, RcisTarget predicts the candidate target genes (i.e.
  genes in the gene-set that are ranked above the leading edge).

- [Rgin](https://bioconductor.org/packages/Rgin) C++ implementation
  of SConES.

- [RGMQL](https://bioconductor.org/packages/RGMQL) This package
  brings the GenoMetric Query Language (GMQL) functionalities into
  the R environment. GMQL is a high-level, declarative language to
  manage heterogeneous genomic datasets for biomedical purposes,
  using simple queries to process genomic regions and their metadata
  and properties. GMQL adopts algorithms efficiently designed for big
  data using cloud-computing technologies (like Apache Hadoop and
  Spark) allowing GMQL to run on modern infrastructures, in order to
  achieve scalability and high performance. It allows to create,
  manipulate and extract genomic data from different data sources
  both locally and remotely. Our RGMQL functions allow complex
  queries and processing leveraging on the R idiomatic paradigm. The
  RGMQL package also provides a rich set of ancillary classes that
  allow sophisticated input/output management and sorting, such as:
  ASC, DESC, BAG, MIN, MAX, SUM, AVG, MEDIAN, STD, Q1, Q2, Q3 (and
  many others). Note that many RGMQL functions are not directly
  executed in R environment, but are deferred until real execution is
  issued.

- [RNAdecay](https://bioconductor.org/packages/RNAdecay) RNA
  degradation is monitored through measurement of RNA abundance after
  inhibiting RNA synthesis. This package has functions and example
  scripts to facilitate (1) data normalization, (2) data modeling
  using constant decay rate or time-dependent decay rate models, (3)
  the evaluation of treatment or genotype effects, and (4) plotting
  of the the data and models. Data Normalization: functions and
  scripts make easy the normalization to the initial (T0) RNA
  abundance, as well as a method to correct for artificial inflation
  of Reads per Million (RPM) abundance in global assesements as the
  total size of the RNA pool deacreases. Modeling: Normalized data is
  then modeled using maximum likelihood to fit parameters. For making
  treatment or genotype comparisons (up to four), the modeling step
  models all possible treatement effects on each gene by repeating
  the modeling with constraints on the model parameters (i.e., the
  decay rate of treatments A and B are modeled once with them being
  equal and again allowing them to both vary independently). Model
  Selection: The AICc value is calculated for each model, and the
  model with the lowest AICc is chosen. Modeling results of selected
  models are then compiled into a single data frame. Graphical
  Plotting: a function is provided to easily visualize the data and
  the selected model using ggplot2 package functions.

- [RSeqAn](https://bioconductor.org/packages/RSeqAn) Headers from the
  SeqAn C++ library for easy of usage in R.

- [rWikiPathways](https://bioconductor.org/packages/rWikiPathways)
  Use this package to interface with the WikiPathways API.

- [scFeatureFilter](https://bioconductor.org/packages/scFeatureFilter)
  An R implementation of the correlation-based method developed in
  the Joshi laboratory to analyse and filter processed single-cell
  RNAseq data. It returns a filtered version of the data containing
  only genes expression values unaffected by systematic noise.

- [scmeth](https://bioconductor.org/packages/scmeth) Functions to
  analyze methylation data can be found here. Some functions are
  relevant for single cell methylation data but most other functions
  can be used for any methylation data. Highlight of this workflow is
  the comprehensive quality control report.

- [Sconify](https://bioconductor.org/packages/Sconify) This package
  does k-nearest neighbor based statistics and visualizations with
  flow and mass cytometery data. This gives tSNE maps"fold change"
  functionality and provides a data quality metric by assessing
  manifold overlap between fcs files expected to be the same. Other
  applications using this package include imputation, marker
  redundancy, and testing the relative information loss of lower
  dimension embeddings compared to the original manifold.

- [SDAMS](https://bioconductor.org/packages/SDAMS) This Package
  utilizes a Semi-parametric Differential Abundance analysis (SDA)
  method for metabolomics and proteomics data from mass spectrometry.
  SDA is able to robustly handle non-normally distributed data and
  provides a clear quantification of the effect size.

- [SEPIRA](https://bioconductor.org/packages/SEPIRA) SEPIRA (Systems
  EPigenomics Inference of Regulatory Activity) is an algorithm that
  infers sample-specific transcription factor activity from the
  genome-wide expression or DNA methylation profile of the sample.

- [seqsetvis](https://bioconductor.org/packages/seqsetvis) seqsetvis
  enables the visualization and analysis of multiple genomic
  datasets. Although seqsetvis was designed for the comparison of
  mulitple ChIP-seq datasets, this package is domain-agnostic and
  allows the processing of multiple genomic coordinate files
  (bed-like files) and signal files (bigwig files or bam pileups).

- [sevenC](https://bioconductor.org/packages/sevenC) Chromatin
  looping is an essential feature of eukaryotic genomes and can bring
  regulatory sequences, such as enhancers or transcription factor
  binding sites, in the close physical proximity of regulated target
  genes. Here, we provide sevenC, an R package that uses protein
  binding signals from ChIP-seq and sequence motif information to
  predict chromatin looping events. Cross-linking of proteins that
  bind close to loop anchors result in ChIP-seq signals at both
  anchor loci. These signals are used at CTCF motif pairs together
  with their distance and orientation to each other to predict
  whether they interact or not. The resulting chromatin loops might
  be used to associate enhancers or transcription factor binding
  sites (e.g., ChIP-seq peaks) to regulated target genes.

- [SIAMCAT](https://bioconductor.org/packages/SIAMCAT) Pipeline for
  Statistical Inference of Associations between Microbial Communities
  And host phenoTypes (SIAMCAT). A primary goal of analyzing
  microbiome data is to determine changes in community composition
  that are associated with environmental factors. In particular,
  linking human microbiome composition to host phenotypes such as
  diseases has become an area of intense research. For this, robust
  statistical modeling and biomarker extraction toolkits are
  crucially needed. SIAMCAT provides a full pipeline supporting data
  preprocessing, statistical association testing, statistical
  modeling (LASSO logistic regression) including tools for evaluation
  and interpretation of these models (such as cross validation,
  parameter selection, ROC analysis and diagnostic model plots).

- [signet](https://bioconductor.org/packages/signet) An R package to
  detect selection in biological pathways. Using gene selection
  scores and biological pathways data, one can search for
  high-scoring subnetworks of genes within pathways and test their
  significance.

- [singleCellTK](https://bioconductor.org/packages/singleCellTK) Run
  common single cell analysis directly through your browser including
  differential expression, downsampling analysis, and clustering.

- [singscore](https://bioconductor.org/packages/singscore) A simple
  single-sample gene signature scoring method that uses rank-based
  statistics to analyze the sample's gene expression profile. It
  scores the expression activities of gene sets at a single-sample
  level.

- [SparseSignatures](https://bioconductor.org/packages/SparseSignatures)
  Point mutations occurring in a genome can be divided into 96
  categories based on the base being mutated, the base it is mutated
  into and its two flanking bases. Therefore, for any patient, it is
  possible to represent all the point mutations occurring in that
  patient’s tumor as a vector of length 96, where each element
  represents the count of mutations for a given category in the
  patient. A mutational signature represents the pattern of mutations
  produced by a mutagen or mutagenic process inside the cell. Each
  signature can also be represented by a vector of length 96, where
  each element represents the probability that this particular
  mutagenic process generates a mutation of the 96 above mentioned
  categories. In this R package, we provide a set of functions to
  extract and visualize the mutational signatures that best explain
  the mutation counts of a large number of patients.

- [srnadiff](https://bioconductor.org/packages/srnadiff) Differential
  expression of small RNA-seq when reference annotation is not given.

- [SummarizedBenchmark](https://bioconductor.org/packages/SummarizedBenchmark)
  This package defines the BenchDesign and SummarizedBenchmark
  classes for building, executing, and evaluating benchmark
  experiments of computational methods. The SummarizedBenchmark class
  extends the RangedSummarizedExperiment object, and is designed to
  provide infrastructure to store and compare the results of applying
  different methods to a shared data set. This class provides an
  integrated interface to store metadata such as method parameters
  and software versions as well as ground truths (when these are
  available) and evaluation metrics.

- [TCGAutils](https://bioconductor.org/packages/TCGAutils) A suite of
  helper functions for checking and manipulating TCGA data including
  data obtained from the curatedTCGAData experiment package. These
  functions aim to simplify and make working with TCGA data more
  manageable.

- [TFEA.ChIP](https://bioconductor.org/packages/TFEA.ChIP) Package to
  analize transcription factor enrichment in a gene set using data
  from ChIP-Seq experiments.

- [TFutils](https://bioconductor.org/packages/TFutils) Package to
  work with TF data.

- [TissueEnrich](https://bioconductor.org/packages/TissueEnrich) The
  TissueEnrich package is used to calculate enrichment of
  tissue-specific genes in a set of input genes. For example, the
  user can input the most highly expressed genes from RNA-Seq data,
  or gene co-expression modules to determine which tissue-specific
  genes are enriched in those datasets. Tissue-specific genes were
  defined by processing RNA-Seq data from the Human Protein Atlas
  (HPA) (Uhlén et al. 2015), GTEx (Ardlie et al. 2015), and mouse
  ENCODE (Shen et al. 2012) using the algorithm from the HPA (Uhlén
  et al. 2015).The hypergeometric test is being used to determine if
  the tissue-specific genes are enriched among the input genes. Along
  with tissue-specific gene enrichment, the TissueEnrich package can
  also be used to define tissue-specific genes from expression
  datasets provided by the user, which can then be used to calculate
  tissue-specific gene enrichments.

- [Trendy](https://bioconductor.org/packages/Trendy) Trendy
  implements segmented (or breakpoint) regression models to estimate
  breakpoints which represent changes in expression for each
  feature/gene in high throughput data with ordered conditions.

- [tRNAscanImport](https://bioconductor.org/packages/tRNAscanImport)
  The package imports the result of tRNAscan-SE as a GRanges object.

- [TTMap](https://bioconductor.org/packages/TTMap) TTMap is a
  clustering method that groups together samples with the same
  deviation in comparison to a control group. It is specially useful
  when the data is small. It is parameter free.

- [TxRegInfra](https://bioconductor.org/packages/TxRegInfra) This
  package provides interfaces to genomic metadata employed in
  regulatory network creation, with a focus on noSQL solutions.
  Currently quantitative representations of eQTLs, DnaseI
  hypersensitivity sites and digital genomic footprints are assembled
  using an out-of-memory extension of the RaggedExperiment API.

- [vidger](https://bioconductor.org/packages/vidger) The aim of
  vidger is to rapidly generate information-rich visualizations for
  the interpretation of differential gene expression results from
  three widely-used tools: Cuffdiff, DESeq2, and edgeR.


New Data Experiment Packages
=====================

There are 16 new data experiment packages in this release of Bioconductor.

- [ASICSdata](https://bioconductor.org/packages/ASICSdata) 1D NMR
  example spectra and additional data for use with the ASICS package.
  Raw 1D Bruker spectral data files were found in the MetaboLights
  database (https://www.ebi.ac.uk/metabolights/, study MTBLS1).

- [BloodCancerMultiOmics2017](https://bioconductor.org/packages/BloodCancerMultiOmics2017)
  The package contains data of the Primary Blood Cancer Encyclopedia
  (PACE) project together with a complete executable transcript of
  the statistical analysis and reproduces figures presented in the
  paper "Drug-perturbation-based stratification of blood cancer" by
  Dietrich S, Oles M, Lu J et al., J. Clin. Invest. (2018)
  128(1):427-445. doi:10.1172/JCI93801.

- [ChIC.data](https://bioconductor.org/packages/ChIC.data) This
  package contains annotation and metagene profile data for the ChIC
  package.

- [CLLmethylation](https://bioconductor.org/packages/CLLmethylation)
  The package includes DNA methylation data for the primary Chronic
  Lymphocytic Leukemia samples included in the Primary Blood Cancer
  Encyclopedia (PACE) project. Raw data from the 450k DNA methylation
  arrays is stored in the European Genome-Phenome Archive (EGA) under
  accession number EGAS0000100174. For more information concerning
  the project please refer to the paper "Drug-perturbation-based
  stratification of blood cancer" by Dietrich S, Oles M, Lu J et al.,
  J. Clin. Invest. (2018) and R/Bioconductor package
  BloodCancerMultiOmics2017.

- [HDCytoData](https://bioconductor.org/packages/HDCytoData) Data
  package containing a set of high-dimensional cytometry data sets
  saved in SummarizedExperiment and flowSet Bioconductor object
  formats, including row and column meta-data describing samples,
  cell populations (clusters), and protein markers.

- [hgu133plus2CellScore](https://bioconductor.org/packages/hgu133plus2CellScore)
  The CellScore Standard Dataset contains expression data from a wide
  variety of human cells and tissues, which should be used as
  standard cell types in the calculation of the CellScore. All data
  was curated from public databases such as Gene Expression Omnibus
  (https://www.ncbi.nlm.nih.gov/geo/) or ArrayExpress
  (https://www.ebi.ac.uk/arrayexpress/). This standard dataset only
  contains data from the Affymetrix GeneChip Human Genome U133 Plus
  2.0 microarrays. Samples were manually annotated using the database
  information or consulting the publications in which the datasets
  originated. The sample annotations are stored in the phenoData slot
  of the expressionSet object. Raw data (CEL files) were processed
  with the affy package to generate present/absent calls (mas5calls)
  and background-subtracted values, which were then normalized by the
  R-package yugene to yield the final expression values for the
  standard expression matrix. The annotation table for the microarray
  was retrieved from the BioC annotation package hgu133plus2. All
  data are stored in an expressionSet object.

- [HMP16SData](https://bioconductor.org/packages/HMP16SData)
  HMP16SData is a Bioconductor ExperimentData package of the Human
  Microbiome Project (HMP) 16S rRNA sequencing data for variable
  regions 1–3 and 3–5. Raw data files are provided in the package as
  downloaded from the HMP Data Analysis and Coordination Center.
  Processed data is provided as SummarizedExperiment class objects
  via ExperimentHub.

- [mCSEAdata](https://bioconductor.org/packages/mCSEAdata) Data
  objects necessary to some mCSEA package functions. There are also
  example data objects to illustrate mCSEA package functionality.

- [MetaGxBreast](https://bioconductor.org/packages/MetaGxBreast) A
  collection of Breast Cancer Transcriptomic Datasets that are part
  of the MetaGxData package compendium.

- [MetaGxOvarian](https://bioconductor.org/packages/MetaGxOvarian) A
  collection of Ovarian Cancer Transcriptomic Datasets that are part
  of the MetaGxData package compendium.

- [MetaGxPancreas](https://bioconductor.org/packages/MetaGxPancreas)
  A collection of pancreatic Cancer transcriptomic datasets that are
  part of the MetaGxData package compendium.

- [RcisTarget.hg19.motifDBs.cisbpOnly.500bp](https://bioconductor.org/packages/RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
  RcisTarget databases: Gene-based motif rankings and annotation to
  transcription factors. This package contains a subset of 4.6k
  motifs (cisbp motifs), scored only within 500bp upstream and the
  TSS. See RcisTarget tutorial to download the full databases,
  containing 20k motifs and search space up to 10kbp around the TSS.

- [RGMQLlib](https://bioconductor.org/packages/RGMQLlib) A package
  that contains scala libraries to call GMQL from R used by RGMQL
  package. It contains a scalable data management engine written in
  Scala programming language.

- [TCGAbiolinksGUI.data](https://bioconductor.org/packages/TCGAbiolinksGUI.data)
  Supporting data for the TCGAbiolinksGUI package. It includes the
  following objects: glioma.gcimp.model, glioma.idhwt.model
  glioma.idhmut.model,glioma.idh.mode, probes2rm,
  maf.tumor,GDCdisease.

- [TENxBrainData](https://bioconductor.org/packages/TENxBrainData)
  Single-cell RNA-seq data for 1.3 million brain cells from E18 mice,
  generated by 10X Genomics.

- [tissueTreg](https://bioconductor.org/packages/tissueTreg) The
  package provides ready to use epigenomes (obtained from TWGBS) and
  transcriptomes (RNA-seq) from various tissues as obtained in the
  study (Delacher and Imbusch 2017, PMID: 28783152). Regulatory T
  cells (Treg cells) perform two distinct functions: they maintain
  self-tolerance, and they support organ homeostasis by
  differentiating into specialized tissue Treg cells. The underlying
  dataset characterises the epigenetic and transcriptomic
  modifications for specialized tissue Treg cells.


New Workflows
=====================

- [BiocMetaWorkflow](https://bioconductor.org/packages/BiocMetaWorkflow)
   Bioconductor Workflow describing how to use BiocWorkflowTools to work with a
   single R Markdown document to submit to both Bioconductor and F1000Research.

- [simpleSingleCell](https://bioconductor.org/packages/simpleSingleCell)
  This workflow implements a low-level analysis pipeline for scRNA-seq data
  using scran, scater and other Bioconductor packages. It describes how to
  perform quality control on the libraries, normalization of cell-specific
  biases, basic data exploration and cell cycle phase identification. Procedures
  to detect highly variable genes, significantly correlated genes and
  subpopulation-specific marker genes are also shown. These analyses are
  demonstrated on a range of publicly available scRNA-seq data sets.


NEWS from new and existing Software Packages
===================================

[ABSSeq](https://bioconductor.org/packages/ABSSeq)
------

Version: 2017-12-18
Date: 2017-12-18
Category: Updating! Refine qtotal (R funtion: qtotalNormalized
Text:

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

[ADaCGH2](https://bioconductor.org/packages/ADaCGH2)
-------

Changes in version 2.19.2 (2018-04-19):

- More on warning on Windows in linking from Rd files.

Changes in version 2.19.1 (2018-04-17):

- Prevent warning on Windows in linking from Rd files.

[adaptest](https://bioconductor.org/packages/adaptest)
--------

Changes in version 1.0.0:

- The first release of this package was made as part of Bioconductor
  3.7, in April 2018. The adaptest R package provides routines for the
  method first described in the the technical manuscript [1] and the
  software paper [2]: 1. Weixin Cai, Nima S. Hejazi, Alan E. Hubbard.
  Data-adaptive statistics for multiple hypothesis testing in
  high-dimensional settings. 2. Weixin Cai, Alan E. Hubbard, Nima S.
  Hejazi. adaptest: Data-Adaptive Statistics for High-Dimensional
  Testing in R.

[affycoretools](https://bioconductor.org/packages/affycoretools)
-------------

Version: 1.52.0
Category: Changed plotting for makeImages to use ggplot2 rather than
        lattice

Version: 1.52.0
Category: Added a new function, vennInLine, intended to allow for
        placement of

Version: 1.52.0
Category: Venn diagrams in Rmarkdown-derived HTML documents,
        particularly those

Version: 1.52.0
Category: that use the BiocStyle package for formatting. The Venn
        diagrams

Version: 1.52.0
Category: contain clickable links that will open HTML tables containing
        results

Version: 1.52.0
Category: for the genes that are found within that cell of the diagram

[ALDEx2](https://bioconductor.org/packages/ALDEx2)
------

Changes in version 1.10.1:

- fixed a bug whereby the order of samples affected the centring of the
  data

[amplican](https://bioconductor.org/packages/amplican)
--------

Changes in version 1.1.2:

- ampliCan supports now other aligners with SAM format eg. bwa mem
  (cigarsToEvents) and pair format eg. EMBOSS needleall (pairToEvents)

- HDR quantification

- more then 10x speed up of extracting events from alignments

[AneuFinder](https://bioconductor.org/packages/AneuFinder)
----------

Changes in version 1.7.2:

NEW FEATURES

- New method 'edivisive' available.

- Breakpoint detection available for both Aneufinder(...,
  strandseq=TRUE) and Aneufinder(..., strandseq=FALSE)

- A stepsize for a sliding window can be selected in addition to the
  binsize for method "HMM". This improves resolution of detected
  breakpoints.

- Breakpoints for Aneufinder(..., strandseq=TRUE) are reported with
  confidence intervals in folder BROWSERFILES.

- Breakpoint detection for Aneufinder(..., strandseq=TRUE) has an
  additional breakpoint refinement step which improves localization of
  breakpoints.

SIGNIFICANT USER-LEVEL CHANGES

- Reorganized output folder structure and added README.txt

- Renamed parameter "plot.SCE" to "plot.breakpoints".

- GC correction is now done with a loess-fit by default instead of the
  quadratic fit. This should improve accuracy.

- Default epsilon is now 'eps=0.01' (was 'eps=0.1' before) for method
  "HMM".

[AnnotationFilter](https://bioconductor.org/packages/AnnotationFilter)
----------------

Changes in version 1.3.1:

NEW FEATURES

- Add DoubleFilter

[AnnotationHub](https://bioconductor.org/packages/AnnotationHub)
-------------

Changes in version 2.12.0:

BUG FIXES

- Prompt for permission when downloading many (more than
  AnnotationHubOption("MAX_DOWNLOADS")) resources.

MODIFICATIONS

- Moved readMetadataFromCsv back to AnnotationHubData.

- Use AnnotationHubData::makeAnnotationHubMetadata to validate
  metadata.csv

[AnnotationHubData](https://bioconductor.org/packages/AnnotationHubData)
-----------------

Changes in version 1.10.0:

MODIFICATIONS

- Moved readMetadataFromCsv back to AnnotationHubData.

- Use AnnotationHubData::makeAnnotationHubMetadata to validate
  metadata.csv

- readMetadataFromCsv is now internal function

[annotatr](https://bioconductor.org/packages/annotatr)
--------

Changes in version 1.6.0:

NEW FEATURES

- Add support for chicken (galGal5).

USER-FACING CHANGES

- Add the ability to facet over two variables in plot_numerical().

- Add the ability to keep duplicate regions in summarize_categorical()
  and plot_categorical(). This is accomplished with the 'by' parameter
  in the former and by the 'x' and 'fill' parameters in the latter, and
  passing their contents into the '.dots' parameter of
  dplyr::distinct_().

- Make TxDb and OrgDb packages Suggests instead of Imports. NOTE: This
  saves space, but also requires downloading the appropriate packages
  as needed.

- Add list_env() function to the annotatr_cache environment to see what
  custom annotations have been read in and added to the cache.

BUGFIXES

- Replace dplyr::summarize_each_() with dplyr::summarize_at() in line
  with deprecation in the dplyr package.

- Prefix builtin_ functions with annotatr:: so that packages that
  Import annotatr don't encounter errors.

[anota2seq](https://bioconductor.org/packages/anota2seq)
---------

Changes in version 1.0.1:

- We changes semantics buffering down into buffered (mRNA up) and
  buffering up into buffered (mRNA down) - so the interpretation is
  more intuitive.

- Add customization options to anota2seqPlotFC and
  anotat2seqPlotPvalues functions

- Fixed a bug where an unitialized object "regModes" was being accessed
  causing an uninformative error message.

[aroma.light](https://bioconductor.org/packages/aroma.light)
-----------

Changes in version 3.9.1 (2017-12-19):

NEW FEATURES

- robustSmoothSpline() now supports using Tukey's biweight (in addition
  to already exising L1) estimators.  See argument 'method'.  Thanks to
  Aaron Lun at the Cancer Research UK Cambridge Institute for adding
  this feature.

Changes in version 3.9.0 (2017-10-30):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC 3.7 for R (>= 3.5.0).

[ASICS](https://bioconductor.org/packages/ASICS)
-----

Version: 0.99.0
Category: First version submitted to Bioconductor
Text:

Version: 0.99.0
Category: New features
Text: preprocessing functions: baseline correction, alignment,

Version: 0.99.0
Category: normalisation
Text: library of pure spectra management functions

Version: 0.99.0
Category: normalisation
Text: binning algorithm

Version: 0.99.0
Category: normalisation
Text: diagnosis tool to assess the quality of the quantification

Version: 0.99.0
Category: normalisation
Text: post-quantification statistical analysis functions: PCA, OPLS-DA

Version: 0.99.0
Category: and Kruskal-Wallis test
Text:

Version: 0.99.0
Category: Improvements
Text: user's guide

Version: 0.99.0
Category: Improvements
Text: spectra quantification

Version: 0.99.0
Category: Improvements
Text: parallel with BiocParallel package

Version: 0.99.0
Category: PREVIOUS VERSIONS ON CRAN
Text:

Version: 0.99.0
Category: New features
Text: data importation from Bruker files

Version: 0.99.0
Category: New features
Text: spectra quantification

Version: 0.99.0
Category: New features
Text: user's guide

[ASpli](https://bioconductor.org/packages/ASpli)
-----

Changes in version 1.5.1:

BUG FIXES

- readCounts supports identical gene names in different chromosomes.

- Fixed default value for minAnchor argument of readCounts code, man
  page and vignette.

[attract](https://bioconductor.org/packages/attract)
-------

Version: 1.31.1
Date: 2018-03-28
Category: Fix bug when custom gene set and expression set both have
        same type of gene IDs

[AUCell](https://bioconductor.org/packages/AUCell)
------

Changes in version 1.2.0:

- Added shiny app to modify thresholds

- Added support for sparse matrices (dgCMatrix)

- New function: AUCell_plotTSNE()

- AUC values are normalized to max 1 by default.

[BASiCS](https://bioconductor.org/packages/BASiCS)
------

Changes in version 1.1.63 (2018-04-27):

- Extra warning message added when `scran` normalisation fails to
  produce positive size factors

Changes in version 1.1.62 (2018-04-21):

- `Regression` is now a compulsory parameter in `BASiCS_MCMC`, hidden
  functions, documentation and unit tests updated accordingly.

- Minor updates in vignette to clarify default parameter in
  `newBASiCS_Data`

- `BASiCS_Filter` adapted to no-spikes case and minor bug resolved

- Added unit test for `BASiCS_Filter`

- Minor bug fixed for no-spikes case in `BASiCS_DenoisedCounts`

Changes in version 1.1.61:

- `Regression` is now a compulsory parameter in `BASiCS_MCMC`, hidden
  functions and SOME unit tests updated accordingly.

Changes in version 1.1.60 (2018-04-19):

- Updated DESCRIPTION file to introduce residual over-dispersion
  testing

Changes in version 1.1.59 (2018-04-19):

- Documentation of `RefFreq` in `BASiCS_Chain` class.

Changes in version 1.1.58 (2018-04-19):

- Updates in vignette to fully incorporate the regression and non-spike
  BASiCS models.

Changes in version 1.1.57 (2018-04-16):

- Minor changes in NEWS file

Changes in version 1.1.56 (2018-04-16):

- Minor changes in `BASiCS_MCMCcppNoSpikes.cpp` to avoid warnings in
  BioC

Changes in version 1.1.55 (2018-04-16):

- Speed-up of `BASiCS_TestDE` by avoiding unnecessary calculations

Changes in version 1.1.54 (2018-04-15):

- CITATION file updated to include latest BASiCS bioRxiv preprint

Changes in version 1.1.53 (2018-04-13):

- Fixed bug in plotting for TestDE function

- Added extra code in vignette to construct BASiCS_Data with no spikes

- Added plots to show changes in residual variability during testing

Changes in version 1.1.52 (2018-04-09):

- Added extra sentence to vignette describing the BASiCS_Data
  generation when SpikeInfo is not available

Changes in version 1.1.51 (2018-03-26):

- Adjusted BASiCS_TestDE function to output data.frames containing same
  genes.

Changes in version 1.1.50 (2018-03-26):

- Removed OffSet correction for phi in the BASiCS_TestDE function

Changes in version 1.1.49 (2018-03-22):

- Fixed link to `BASiCS_Chain` documentation in `BASiCS_showFit` doc

Changes in version 1.1.48 (2018-03-22):

- Excluded epsilon object and RefFreq object from validity test in
  Chain Class

Changes in version 1.1.47 (2018-03-21):

- Reduced number of genes in data examples

- Unit test updated accordingly

- Nils added as maintainer in DESCRIPTION

Changes in version 1.1.46 (2018-03-21):

- Split C++ code into individual functions

- One file per main function and individual files for utility function

- Created header files for main functions and utility functions

Changes in version 1.1.45 (2018-03-21):

- Messages within C++ code moved to functions

Changes in version 1.1.44 (2018-03-21):

- `HiddenVarDecomp` adapted to no-spikes case

- Unit test for `BASiCS_VarDecomp` expanded accordingly

- Specification of `colnames` and `rownames` added to `newBASiCS_Data`

- `colnames` and `rownames` are now specified as part of the output of
  `BASiCS_DenoisedCounts` and `BASiCS_DenoisedRates`

Changes in version 1.1.43 (2018-03-20):

- Batch information moved to `colData` instead of `metadata` when
  building a `SingleCellExperiment` object via the `newBASiCS_Data`
  function.  This facilitates data subsetting.

- `BASiCS_MCMC` and `HiddenBASiCS_MCMC_InputCheck` updated accordingly

- `SummarizedExperiment::colData` and `S4Vectors::DataFrame` added as
  imports in NAMESPACE

Changes in version 1.1.42 (2018-03-20):

- `OrderVariable` default value set as "GeneIndex" in `BASiCS_TestDE`

- Fixed usage of `OrderVariable` in `BASiCS_TestDE` regression case

- Associated unit test updated accordingly

- Missing denoised rates/counts check added to unit tests for no-spikes
  case

- CITATION file has been added

Changes in version 1.1.41 (2018-03-20):

- `BASiCS_MCMC`: Default value for `WithSpikes` set as FALSE when there
  are no spike-ins in available in the data

Changes in version 1.1.40 (2018-03-20):

- Fixes small bug in `BASiCS_DenoisedRates`, no-spikes case

- Adds denoised rates/counts check to unit tests for no-spikes case

Changes in version 1.1.39 (2018-03-20):

- Denoised counts/rates added to unit test for regression case

- `BASiCS_DenoisedCounts` adapted for no-spikes case

- `BASiCS_DenoisedRates` adapted for no-spikes case

Changes in version 1.1.38 (2018-03-20):

- Validity checks updated for `BASiCS_Chain` and `BASiCS_Summary`
  classes. This accounts for different parameter configurations

- `Summary` method to exclude `RefFreq` for the no-spikes case

- Unit tests fo no-spikes case updated accordingly

Changes in version 1.1.37 (2018-03-20):

- `phi` removed from the output of `BASiCS_MCMC` for no-spikes case

- Unit tests updated accordingly

- Checks for parameter names added to unit tests

Changes in version 1.1.36 (2018-03-18):

- Fixed bug in `newBASiCS_Data` to define a default value for `Tech`
  when spike-ins are not available. Thanks to Muad Abd El Hay (@Cumol)
  for pointing out this issue.

Changes in version 1.1.35 (2018-03-18):

- New unit tests (more modular): example data, starting values MCMC,
  MCMC sampler for fixed starting values (spikes/no-spikes;
  regression/no-regression)

Changes in version 1.1.34 (2018-01-31):

- `HiddenBASiCS_MCMC_Start` edited to account for new default parameter
  values in `scran::computeSumFactors` call

- Unit tests updated accordingly

- Imports/Suggests packages are now in alphabetic order

- Minor style changes suggested by BiocCheck

Changes in version 1.1.33 (2018-01-30):

- In `HiddenBASiCS_MCMC_Start`, `positive = TRUE` added to
  `scran::computeSumFactors` when estimated size factor contain invalid
  values. Thanks to Mike Morgan for suggesting this solution.

Changes in version 1.1.32 (2018-01-29):

- Fix to ensure that `BASiCS_MCMC` function adds cell labels for $\phi$

- `makeExampleBASiCS_Data` now produces valid `colnames` (cell labels)

- Fixed colnames of $\phi$ parameter in all data examples

- Extended validity test for `BASiCS_Chain` class

- `rownames` and `colnames` methods created for `BASiCS_Chain` class

- Added example and unit test for `BASiCS_VarianceDecomp` function

Changes in version 1.1.31 (2018-01-22):

- Updated documentation for `BASiCS_showFit` method

Changes in version 1.1.30 (2018-01-22):

- Commit to triger a new built

Changes in version 1.1.29 (2017-12-21):

- `showFit` generic renamed as `BASiCS_showFit` (to avoids potential
  conflicts)

Changes in version 1.1.28 (2017-12-20):

- Vignette header as in BiocStyle

- In `BASiCS_TestDE` function: `ProbThresholdE` argument renamed as
  `ProbThresholdR`

- In `BASiCS_TestDE` function: `ProbThresholdE` argument renamed as
  `ProbThresholdR`

- Updated documentation for example `BASiCS_Chain` objects

- Minor additional changes to the documentation

Changes in version 1.1.27 (2017-12-20):

- Minor typo resolved in `displayChainBASiCS` method

- Reduced size of example objects

- Unit tests updated accordingly

- In `BASiCS_TestDE` function: `PsiE` argument renamed as `EpsilonR`

Changes in version 1.1.26 (2017-12-20):

- Vignette bibliography moved to .bib file

- Added diagram to summarise different implementations

- Reduced size of example objects

- Unit tests updated accordingly

- `subset` method created for `BASiCS_Chain` objects

- Removal of parameter `lambda` from `BASiCS_Chain` objects (for
  storage)

- R/Methods.R + R/Classes.R updated accordingly

- `newBASiCS_Chain` updated accordingly.

- Code format for `showFit`

Changes in version 1.1.24 (2017-12-18):

- Reduced no of iterations in unit test for quicker testing

- library `hexbin` added to imports (required for `showFit`)

- Vignette edits (simplification of Quick Start)

Changes in version 1.1.23 (2017-12-18):

- Minor improvements to documentation

Changes in version 1.1.22 (2017-12-15):

- Documentation on ChainSCReg and ChainRNAReg datasets added

- Fixed showFit method

- Updated dpi option in vignette to reduce .html size

- Default value for `ConstrainProp` set to 0.20 (no-spikes)

- Unit tests updated accordingly (no-spikes)

Changes in version 1.1.21 (2017-12-13):

- `ConstrainProp` added as an optional argument for `BASiCS_MCMC`
  (no-spikes)

Changes in version 1.1.20 (2017-12-11):

- Updated dpi option in vignette

Changes in version 1.1.19 (2017-12-09):

- extended vignette to describe regression and non-spike case

- changed markdown to rmarkdown engine to compile vignette to include
  table of contents

Changes in version 1.1.18 (2017-12-08):

- Minor typo in `HiddenBASiCS_MCMC_NoSpikesParam` resolved

Changes in version 1.1.17 (2017-12-07):

- Stop scaling of size factors for no-spikes case
  (`HiddenBASiCS_MCMC_Start`)

- Unit testing updated accordingly

- Extra unit test to ensure that no-spikes results match for objects
  with/without spikes when setting `WithSpikes = FALSE`

- Modified constrain for no-spikes cases (to deal with non-zero genes)

- Default value of `ConstrainType` set to 1 in `BASiCS_MCMC` (no-spikes
  only)

Changes in version 1.1.16 (2017-11-29):

- Change in `BASiCS_MCMC` no-spikes so that only genes with zero total
  counts are excluded from the identifiability constrain

Changes in version 1.1.15 (2017-11-29):

- Typo fixed in `AtLeast2Cells` (`BASiCS_MCMC` function) definition for
  no-spikes + regression case when original data has spikes

Changes in version 1.1.14 (2017-11-29):

- Minor bug resolved for no-spikes case when original data has spikes

- New unit test to assess that case

Changes in version 1.1.13 (2017-11-28):

- Fixed storage of `RefGene` in no-spikes when `StochasticRef = FALSE`

Changes in version 1.1.12 (2017-11-28):

- Additional changes to improve usage of no-spikes case

Changes in version 1.1.11 (2017-11-24):

- Minor typo in the call to `HiddenBASiCS_MCMCcppNoSpikes` has been
  resolved

Changes in version 1.1.10 (2017-11-27):

- `HiddenBASiCS_MCMC_Start` updated to remove regression-related
  hyper-params

- `HiddenBASiCS_MCMC_ExtraArgs` updated to include regression-related
  hyper-par

- `BASiCS_MCMC` modified accordingly

- Storage of adaptive variances fixed for the no-spike case

Changes in version 1.1.9 (2017-11-24):

- `means` only updated when needed (regression case); unit test updated
  accordingly

- `deltaUpdateNoSpikes` + `deltaUpdateRegNoSpikes` removed as no longer
  required

- Simplification of terms in `deltaUpdateReg`

- Small changes to unify code format in C++ file

- Notation change: `phi` replaced by `s` in no-spikes implementation

- Unit test changed accordingly (no-spike)

- Minor style changes in `BASiCS_MCMC`

- `BASiCS_MCMC` function broken into smaller functions (easier to
  read). This creates the following hidden functions:
  `HiddenBASiCS_MCMC_InputCheck`, `HiddenBASiCS_MCMC_ExtraArgs`,
  `HiddenBASiCS_MCMC_NoSpikesParams`, `HiddenBASiCS_MCMC_OutputStore`
  and `HiddenBASiCS_MCMC_RefFreqStore`

- Merge between no-spikes and regression case completed (code to be
  tested)

- Minor change in storage of reference frequency (no spikes only)

Changes in version 1.1.8 (2017-11-21):

- `muUpdateRegNoSpikes` and `deltaUpdateRegNoSpikes` added to C++ code

- Clean-up of C++ code (repeated debug checks made into function)

- `lambdaUpdateReg` created as a separate function

- Vectorization of calculations related to regression implementation

- Global quantities (e.g. inv(V0)) taken outside the look (regression
  case)

- `sigma2UpdateReg` created as a separate function

- `betaUpdateReg` created as a separate function

- Fixed typo in sigma2 updates (regression case); unit test updated
  accordingly

Changes in version 1.1.7 (2017-11-20):

- Nils added as author in cpp file

- Updated unit test for no-spikes case (changes due to change in stoch
  ref)

- Change in `BASiCS_MCMC` so that only biological counts go to C++

- Added factorizations in `muUpdateReg` and `deltaUpdateReg` (C++ code)

Changes in version 1.1.6 (2017-11-15):

- Added prototype for regression + no-spikes sampler

Changes in version 1.1.5 (2017-11-15):

- Removal of collapsed sampler prototype (from C++ and R) - moved to
  `CollapsedSampler` branch

- Number of genes used for stochastic reference increased to 10% of
  genes (no spikes only)

Changes in version 1.1.4 (2017-11-12):

- Only include genes that are expressed in at least 2 cells per
  condition for differential residual dispersion testing (regression
  case only)

Changes in version 1.1.3 (2017-11-08):

- Collapsed sampler implemented (to be tested)

- Unit tests reverted to original case

- Version number bump

- Store `TableRef` on the `StoreDir` directory (no-spike case only)

Changes in version 1.1.2 (2017-11-08):

- Skeleton for sampler with integrated out s

Changes in version 1.1.1 (2017-11-08):

- 'ggplot2' added to imports

- Default value for SpikeInfo changed to NULL in `newBASiCS_Data`
  function

- Cleaning of `HiddenBASiCS_MCMCcppNoSpikes`

- Cleaning of full conditionals for no-spikes implementation

- Fixed bug in stochastic reference implementation

- Added unit test for no-spikes case (estimation)

- Minor tyle changes to merge regression case

- 'Eta' parameter removed from `Summary` method (regression case only)

- Corrected HPD interval calculation for epsilon and sigma2 (regression
  case)

- Added mark for 'ExcludedFromRegression' genes in `BASiCS_TestDE`

Changes in version 1.1.0 (2017-08-16):

- Minor changes required to pass BiocCheck (used formatR and
  BiocChecks)

- Version number bump after bioconductor 3.6 release

[bcSeq](https://bioconductor.org/packages/bcSeq)
-----

Version: 1.0.1
Category: fixed a bug for windows version that caused by pthread
Text:

Version: 1.0.2
Category: added helper functions to preprocess user files
Text: "trimRead" to trim adpators. "uniqueBar" to extract unique
        barcode sequence from the library file

Version: 1.0.3
Category: change "alignment" to "mapping" in the title
Text:

Version: 1.0.4
Category: change "for" to "in" in the title
Text:

[beachmat](https://bioconductor.org/packages/beachmat)
--------

Changes in version 1.1.13:

- Changed environment variable to BEACHMAT_RPATH for consistency with
  other packages.

- Added native support for transposition and subsetting in
  DelayedMatrix objects.

- Added support for chunk-by-chunk realization of otherwise unsupported
  matrices, including DelayedMatrix objects with other delayed
  operations.

- Added the get_const_col_indexed() method for input matrices,
  especially fast for sparse representations.

- Added the set_col_indexed() and set_row_indexed() methods for output
  matrices.

- Updated vignettes and expanded the test suite.

[bioCancer](https://bioconductor.org/packages/bioCancer)
---------

Version: 1.7.05
Text: omit importFrom clusterProfiler plot

Version: 1.7.04
Text: export Diagrammer and viNetwork to report

Version: 1.7.04
Text: reset image size of circomics to 1024px

Version: 1.7.03
Text: Save network widgets as HTML and png

Version: 1.7.03
Text: use swithc buttom for Networking Tab.

Version: 1.7.03
Text: setwd(~)/ for windows system by setwd(Sys.getenv("R_USER"))

Version: 1.7.02
Text: resolve conflicts renderMetabologram and renderCoffeewheel.
        redefine initCoffeewhell in /htmlwidges

Version: 1.7.02
Text: report circomics widget to markdown document

Version: 1.7.02
Text: Save static wheel as html and png (low resolution)

Version: 1.7.02
Text: need phantomJS to catupe html widget output as png file

Version: 1.7.01
Text: add helps to ? menu

Version: 1.7.01
Text: change stop message

Version: 1.7.01
Text: addResourcePath for figures

Version: 1.7.00
Category: metamorphosis: bioCancer is a radiant.data extension
Text:

Version: 1.7.00
Category: reduce size of package by half 14 -> 7 mb
Text:

[BiocCheck](https://bioconductor.org/packages/BiocCheck)
---------

Changes in version 1.16:

BUG FIXES

- handle interactive BiocCheck() arguments correctly

[BiocFileCache](https://bioconductor.org/packages/BiocFileCache)
-------------

Changes in version 1.3:

NEW FEATURES

- (v. 1.3.35) Save a post download processed file to cache.

- (v. 1.3.28) Add ask = TRUE argument to BiocFileCache().

- (v. 1.3.24) Add function makeBiocFileCacheFromDataFrame to convert
  data.frame to BiocFileCache

- (v. 1.3.19) etag now checked in addition to last modified time to
  determine if local version of file is current

- (v. 1.3.16) importbfc to load output of exportbfc

- (v. 1.3.15) exportbfc allows users to create exportable archive of
  bfc related files

- (v. 1.3.11) bfcisrelative checks for rtype='local' in addition to
  relative

- (v. 1.3.11) Helper function to check rtype for local but relative,
  and web and update if necessary

- (v. 1.3.10) Add function to check portability of BiocFileCache:
  bfcisrelative

- (v. 1.3.10) Add function to convert rpaths for portability:
  bfcrelative

- (v. 1.3.9) Optionally download web resource when adding to cache
  bfcadd(download=TRUE)

SCHEMA CHANGE

- (v. 1.3.36) expires added

- (v. 1.3.19) etag added

- (v. 1.3.9) Last modified time default is NA

USER-VISIBLE CHANGES

- (v. 1.3.38) ... argument exposed and pased to GET, this includes
  bfcadd which original ... was pasesed to file.copy. This use of ...
  is used in bfc functions bfcadd(), bfcupdate() and bfcdownload()
  which could potential download the file.

- (v. 1.3.30) bfcnew(), bfcadd(), accept vector arguments; performance
  improvements and

- (v. 1.3.25) bfcinfo and bfcquery will show full rpaths even if stored
  as relative

- (v. 1.3.22) prompt user only once when using default cache

- (v. 1.3.17) prompt user when overwriting exisiting file

- (v. 1.3.16) add importbfc to extract output of exportbfc and load bfc
  object

- (v. 1.3.16) added exportbfc to export bfc or subset of bfc

- (v. 1.3.16) bfcisportblae and bfcportbale removed for exportbfc

- (v. 1.3.14) upon bfc creation, option to update to most current
  schema

- (v. 1.3.13) bfcportable operates over all offending ids instead of
  asking individually

- (v. 1.3.12) bfcisrelative/bfcrelative changed to
  bfcisportbale/bfcportable

- (v. 1.3.11) web resource rpaths are stored as relative links

- (v. 1.3.9) If web resource was not downloaded, bfcneedsupdate is TRUE

- (v. 1.3.9) Local and non downloaded web will have last modified time
  as NA

- (v. 1.3.6) Update how default cache location is determined

- (v. 1.3.1) Expose GET::config argument to web resource functions

BUG CORRECTION

- (v. 1.3.40) correct trycatch with cache_info. cache_info bug
  workaround orginally output NA even if present, now manually grab
  values if present

- (v. 1.3.37) correct which functions update access_time

- (v. 1.3.26) patch for cache_info after returning etag and
  last_modified

- (v. 1.3.4) patch for cache_info bug

[BiocInstaller](https://bioconductor.org/packages/BiocInstaller)
-------------

Changes in version 1.30.0:

NEW FEATURES

- biocLite() supports github repositories using the remotes package,
  rather than devtools. This change should be transparent to end users.
  (From Peter Hickey
  https://github.com/Bioconductor/BiocInstaller/issues/4)

[BiocParallel](https://bioconductor.org/packages/BiocParallel)
------------

Changes in version 1.14:

BUG FIXES

- (v 1.13.1) bpiterate,serial-method does not unlist the result of FUN
  before passing to REDUCE.

[BiocStyle](https://bioconductor.org/packages/BiocStyle)
---------

Changes in version 2.8.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Add 'relative_path' argument to 'pdf_document'

- Add argument 'titlecaps' to 'html_document'

- Defunct deprecated functions

BUG FIXES

- Fix a bug in generating navigation links which resulted in an empty
  document when there were no subsections to merge

[biomaRt](https://bioconductor.org/packages/biomaRt)
-------

Changes in version 2.36.0:

BUG FIXES

- Patched problem returning the list of available datasets, if the
  description of one or more datasets included an apostrophe
  (introduced with new primate species in Ensembl).

- Caught scenario where ensemblRedirect=FALSE was still being ignored.

- Changed query submission when redirection is detected to cope with
  apparently new behaviour of the Ensembl mirrors.

MINOR CHANGES

- Increase query timeout limit to 5 minutes.

[biosigner](https://bioconductor.org/packages/biosigner)
---------

Changes in version 1.7.4:

INTERNAL MODIFICATIONS

- unit testing on Windows 32 bits disabled

Changes in version 1.7.2:

INTERNAL MODIFICATIONS

- modification in unit testing

[biotmle](https://bioconductor.org/packages/biotmle)
-------

Changes in version 1.4.0:

- An updated release of this package for Bioconductor 3.7, released
  April 2018.

- This release primarily implements minor changes, including the use of
  colors in the plots produced by the visualization methods.

[bnbc](https://bioconductor.org/packages/bnbc)
----

Changes in version 1.1.2:

- Removed unneeded so files.

[BPRMeth](https://bioconductor.org/packages/BPRMeth)
-------

Changes in version 1.6.0:

- Add additional observation models: "binomial", "bernoulli", "beta"
  and "gaussian".

- Add functionlity for Bayesian inference using mean-field variational
  inference for both inferring and clustering profiles.

- Add functionality for Fourier basis functions.

- Update plotting functionality to use ggplot2.

[BrowserViz](https://bioconductor.org/packages/BrowserViz)
----------

Changes in version 2.1.0:

NEW FEATURES

- Javascript code uses npm and webpack

- R code uses latest httpuv from RStudio, with async websocket windows
  support

[CAGEfightR](https://bioconductor.org/packages/CAGEfightR)
----------

Changes in version 0.99:

- Submitted CAGEfightR to Bioconductor

[CAGEr](https://bioconductor.org/packages/CAGEr)
-----

Changes in version 1.21.6:

- mergeCAGEsets accepts CAGEexp objects as input.

- CTSSes are now represented by GPos objects, which are more compact in
  memory and on display, while still inheriting from GRanges.

Changes in version 1.21.5:

- Simplified CAGEexp constructor, similar to the one for CAGEset
  objects.

- New removeStrandInvaders() function to count and remove
  strand-invasion artefacts (see Tang et al., 2013,
  doi:10.1093/nar/gks112).

- Tag clusters are sorted before being stored in the CAGEr objects.

- plotAnnot: - new "facet" option for running ggplot2::facet_wrap() on
  the plots. - better documentation for the "scopes" of the plot. -
  removed "customScope" argument: pass a function directly to "scope"
  argument instead.

- plotCorrelation2: - new methods for SummarizedExperiment, DataFrame,
  data.frame and matrix. - Corrected the correlation matrix output
  (previously the diagonal values were 0 and the other values were
  swapped).

- plotReverseCumulatives: - fitInRange accepts a value of "NULL" to
  turn off power law fitting. - New "legend" argument to remove legend
  when set to "FALSE". - Axis range and labels can be modified with
  xlab/ylab and xlim/ylim. - Ladders steps start on the values instead
  of being centered on.

- setColors: allow lowercase in color names.

Changes in version 1.21.4:

- Corrected a bug that was crashing CAGEset objects when loading more
  than one BAM file.

- Load BAM as gapped alignment with readGAlignments() instead of
  scanBam().  Without this correction, TSS position on minus strand is
  incorrect in case of indels in the read.

- Partial support for loading BAM data in CAGEexp objects.
  (correctSystematicG is not yet implemented)

- Added multicore processing in hanabi function.

Changes in version 1.21.3:

- Enforce syntactically valid sample names in CAGEexp class.

Changes in version 1.21.2:

- Use a plain DataFrame as rowData in seqNameTotals slots.

Changes in version 1.21.1:

BACKWARDS-INCOMPATIBLE CHANGES

- The plotting functions send their output to the graphical device
  instead of writing it to a file.  This makes their use more
  consistent with most plotting functions in R.

NEW FEATURES

- New "CAGEexp" class extending the MultiAssayExperiment class.  It
  stores expression data more efficiently than "CAGEset", and uses core
  Bioconductor typse natively.  For backwards compatibility it also
  support many of the original generic functions for "CAGEset" objects.

- New "CTSS" and "TagClusters" classes wrapping GRanges objects, for
  more type safety.

- New functions for quality controls such as plotAnnot() or
  hanabiPlot(). See the CAGEexp vignette for details.

- Data export as DESeqDataSet object for DESeq2 with the new
  "consensusClustersDESeq2" and "GeneExpDESeq2" functions.

- New "bedctss" format to load the FANTOM5 and FANTOM6 CAGE data.

- New "CAGEscanMolecule" format to load CAGEscan 3.0 data.

- Multicore parallelisation with BiocParallel instead of parallel.

- New function sampleList() to help looping on samples with lapply().

- New plotCorrelation2() function, faster than plotCorrelation()
  because it is plain black and white.

- Multicore loading of CTSS data in CAGEexp object.

OTHER CHANGES

- Example data "exampleCAGEexp", "exampleCAGEexp" and
  "exampleZv9_annot" are now is lazy-loaded.

- Passes R CMD check without errors or notes.

- NULL can be passed as genome name, to circumvent the requirement for
  a BSgenome object when actually not needing one.

- In CAGEexp objects, expression quantile positions are given relative
  to the cluster start site.

- For performance reasons, the positions of a quantile Q is now
  calculated as the position of the first base where cumulative
  expression is higher or equal to Q% of the total expression of a
  cluster.

DOCUMENTATION UPDATES

- Roxygen is used to generate the manual pages.

- A new vignette describes the "CAGEexp" class.

[canceR](https://bioconductor.org/packages/canceR)
------

Version: 1.11.01
Text: rm gamma = 1.5 argument for grDevices::rainbow function

Version: 3.1
Text:

Version: 3.2
Text: dimensions levels will be plot. 1- dialogMetOption(): add
        "Circos" argument to make the difference between
        getMetDataMultipleGene() and getListMetData() 2- getGeneList():
        add rm("GeneListMSigDB"", envir="myGlobalEnv")

[Cardinal](https://bioconductor.org/packages/Cardinal)
--------

Changes in version 1.11.2:

NEW FEATURES

- Added 'writeMSIData', 'writeImzML', and 'writeAnalyze' methods for
  writing MSI data to supported file formats

- Added support for on-disk 'processed' imzML (via argument
  'attach.only' in 'readImzML' method)

SIGNIFICANT USER-VISIBLE CHANGES

- Package 'matter' is used for all file I/O now

- Switched from using 'Hashmat' to using 'sparse_mat' class from
  'matter' for 'processed' imzML data

BUG FIXES

- Changed compiler settings for parsing XML so reading large imzML
  files should use much less memory now (but may take slightly longer
  for smaller files)

Changes in version 1.11.1 (2017-10-25):

SIGNIFICANT USER-VISIBLE CHANGES

- Use 'drop=NULL' from now on instead of 'drop=NA' to do endomorphic
  subsetting of 'SimageData' objects

BUG FIXES

- Use 'keys' and 'keys<-' generic from 'matter'

[cbaf](https://bioconductor.org/packages/cbaf)
----

Changes in version 1.2.0 (2018-04-26):

New Features

- heatmapOutput() can now determine the heatmap margines and column and
  row name sizes automatically.

- New image formats TIFF, JPG and BMP in addition to the previous PNG
  file format for heatmapOutput(). They can be chosen from
  processOneStudy() and processMultipleStudies or directly from
  heatmapOutput() function.

- heatmapOutput() now uses two methods for ranking the genes prior to
  generating heatmap(s). One of them is suited for finding genes that
  have unique high values in one or few cancer studies whereas the
  other method aids in detemining genes that possess high values in
  multiple / many cancers.

- If function argumnets are entered wrongly, more meaningful errors
  will appear.

- All functions are improved.

[ccfindR](https://bioconductor.org/packages/ccfindR)
-------

Changes in version 0.99.4 (2018-04-06):

Changes

- Split plot functionality from filter_genes() into plot_genes().

- Removed file name arguments from read_10x() and write_10x().

[CellNOptR](https://bioconductor.org/packages/CellNOptR)
---------

Changes in version 1.25.0 (2018-01-10):

- BUG FIX * fix gray boxes in plotOptimResultsPan * in
  plotOptimResultsPan when errors were greater than 99.9% color were
  white (instead of red) * c simulator bug fix: at time 0, node that
  are inhibited and measured were reset to 0 but inhibitors are off.
  Fix is simple: do not reset inhibitors where time is zero * Fixing
  issue with time 0 not being simulated properly (see
  https://github.com/cellnopt/CellNOptR/issues/6) This fixes regression
  bug following fix made in release 1.11.3 * node name in the sif file
  containing the word AND (e.g. ligand) will not result in an AND-node
  anymore * Fixing issues with matrix subsetting, when the subsetting
  converts the matrix to a vector * plotOptimResultsPan: in the
  computation of root-mean-square error(for coloring the background),
  the NA data is not counted in the number of data points

- CHANGES * Bioconductor's version of the package got merged with the
  Github's version leading to minor changes * readSIF reads only the
  unique interactions/lines from the SIF file * plotOptimResultsPan and
  plotCNOlist plots intermediate cue values (0,1) for CNORode add-on

- NEW FEATURES * add readErrors and writeErrors functions in CNOlist:
  read/write measurement error/variance from/to MIDAS files * toSBML()
  function writes the model to SBMLqual format * makeCNOlist() can
  import MIDAS with multiple, continuous cue levels * add new function
  called crossInhibitedData

[ChemmineR](https://bioconductor.org/packages/ChemmineR)
---------

Changes in version 2.32.0:

NEW FEATURES

- SDFDataTable function, allows viewing compound image and data in a
  web browser.

- read.SDFset can now skip over compounds with syntax errors in sdf
  files.

- The functions 'getIds', 'searchString', and 'searchSim' have been
  modified to query pubchem directly rather than going through the
  intermediate web service 'ChemmineTools'.

[chimeraviz](https://bioconductor.org/packages/chimeraviz)
----------

Changes in version 1.5.6:

CHANGES

- The interface for working with chimeraviz has changed to comply with
  https://github.com/jimhester/lintr. chimeraviz is the same as before,
  but methods have new names.

Changes in version 1.5.1:

NEW FEATURES

- Added support for the mm10 genome.

[ChIPanalyser](https://bioconductor.org/packages/ChIPanalyser)
------------

Changes in version 3.6:

paramters

- Data pre-processing Function

- rtracklayer import support

[chipenrich](https://bioconductor.org/packages/chipenrich)
----------

Changes in version 2.4.0:

NEW FEATURES

- A new function, peaks2genes(), to run the analysis up to, but not
  including, the enrichment testing. Useful for checking QC plots,
  check qualities of peak-to-gene assignments, and easier custom tests.

SIGNIFICANT USER-LEVEL CHANGES

- The hybridenrich() method now returns the same format as chipenrich()
  and polyenrich()

IMPROVEMENTS

- Vignette now describes all available gene sets

BUG FIXES

- Fixed multiAssign weighting method to use the correct weights.

[ChIPexoQual](https://bioconductor.org/packages/ChIPexoQual)
-----------

Version: 1.3.1
Text:

Version: 1.3.2
Text:

[ChIPpeakAnno](https://bioconductor.org/packages/ChIPpeakAnno)
------------

Changes in version 3.13.19:

- Fix the bug in assignChromosomeRegion for the percentage calculation
  when nucleotideLevel equal to true.

Changes in version 3.13.18:

- Fix the bug in assignChromosomeRegion for resetting seqlevels.

Changes in version 3.13.17:

- Fix the bug in annoScore by setting the ignore.strand for punion.

Changes in version 3.13.16:

- Fix the bug in assignChromosomeRegion for the integer overflow.

Changes in version 3.13.14:

- Fix the warnings under windows.

Changes in version 3.13.12:

- Change "class" function to "is" function.

Changes in version 3.13.11:

- Fix the help file for assignChromosomeRegion.

Changes in version 3.13.10:

- Improve binOverRegions.

Changes in version 3.13.9:

- Improve addMetadata.

Changes in version 3.13.8:

- Fix the warnings under windows.

Changes in version 3.13.6:

- Improve the functions featureAlignedSignal and
  featureAlignedExtendSignal

Changes in version 3.13.5:

- add minFragmentSize parameter to estFragmentLength

Changes in version 3.13.4:

- suppress the message from featureAlignedExtendSignal

- modified help files to make the samples run on local test.

Changes in version 3.13.3:

- fix the bug in binOverRegions and binOverGene

Changes in version 3.13.2:

- featureAlignedSignal generate outputs with rownames

Changes in version 3.13.1:

- change the test value from GNAI3;MIR197 to GNAI3 to avoid the error
  because of update of annotation database.

[ChIPseeker](https://bioconductor.org/packages/ChIPseeker)
----------

Changes in version 1.15.2:

- bug fixed for 'overlap = "all"' to consider strand information
  <2017-12-12, Tue>

Changes in version 1.15.1:

- define downstream distance via options(ChIPseeker.downstreamDistance
  = 3000) + https://support.bioconductor.org/p/103135/

[chromstaR](https://bioconductor.org/packages/chromstaR)
---------

Changes in version 1.5.1:

SIGNIFICANT USER-LEVEL CHANGES

- New column 'maxPostInPeak' containing the maximum posterior within
  each peak.

- Score in exported BED files is calculated as
  -10*log10(maxPostInPeak).

- 'changeFDR()' was renamed to 'changeMaxPostCutoff()'.

[ClassifyR](https://bioconductor.org/packages/ClassifyR)
---------

Changes in version 2.0.0:

- Broad support for DataFrame and MultiAssayExperiment data sets by
  feature selection and classification functions.

- The majority of processing is now done in the DataFrame method for
  functions that implement methods for multiple kinds of inputs.

- Elastic net GLM classifier and multinomial logistic regression
  classifier wrapper functions.

- Plotting functions have a new default style using a white background
  with black axes.

- Vignette simplified and uses a new mass cytometry dataset with
  clearer differences between classes to demonstrate classification and
  its performance evaluation.

[clusterExperiment](https://bioconductor.org/packages/clusterExperiment)
-----------------

Changes in version 1.99.4 (2018-04-19):

Changes: 

- Added support for hdf5 files stored in `assay` slot via the
  `HDF5Array` package

- Removed most defaults from RSEC arguments -- pull them from
  underlying functions' defaults.

Changes in version 1.99.3 (2018-04-17):

Changes

- Re-implemented subsampleClustering() and combineMany() to use C++.

- Method "adjP" in `mergeClusters` now allows for further requirement
  that gene have a minimal log-fold change ('logFCcutoff').

Bugs

- Fix bug in `setBreaks` (isPositive and isNegative variables)

- use `stringr::str_sort` to make sort of character values locale
  independent

Changes in version 1.99.2 (2018-03-22):

Bugs

- Fix defaultNDims so that returns minimum of 50 and the minimum
  dimension of data.

- Fix RSEC so still returns results if hit error after clusterMany.

Changes in version 1.99.0 (2018-02-15):

Changes

- MAJOR CHANGE TO DEFINITION OF CLASS: This version consists of a major
  update of how dimensionality reduction and filtering is done. The
  class has been updated to extend the new `SingleCellExperiment`
  class, which save the dimensionality reductions. Furthermore,
  calculating of per-gene statistics, which are usually used for
  filtering, are stored in `colData` of the object and can be easily
  accessed and used for repeated filtering without recalculating. This
  has created a massive change under-the-hood in functions that allow
  dimensionality reduction and filtering. Changes to function names are
  the following: - `transform` is now `transformData` - New functions
  `makeReducedDims` and `makeFilterStats` will calculate (and thus
  store) dimensionality reductions and statistics for filtering the
  data. - New function `filterData` will return the filtered data as a
  matrix - New functions `listBuiltInReducedDims` and
  `listBuiltInFilterStats` give the list of currently available
  functions for dimensionality reduction and filtering statistics,
  respectively. - Filtering on arbitrary statistics and user-defined
  dimensionality reduction can used in `clusterMany` and related
  functions, so long as they are saved in the appropriate slots of the
  object.

- Changed the following functions/arguments to be consistent with
  SingleCellExperiment naming conventions and improve distinction
  between terminology of cluster and clustering. - Capitalized
  constructor functions. Now: `ClusterFunction()` and
  `ClusterExperiment()` - `nPCADims` now changed to `nReducedDims` in
  clusterMany-related functions - `nVarDims` now changed to
  `nFilterDims` in clusterMany-related functions - `dimReduce` argument
  now changed to `reduceMethod` across functions - `ndims` to `nDims`
  in `clusterSingle` and `makeDendrogram` to keep consistency. -
  `plotDimReduce` to `plotReducedDims` - Changed `nClusters` to
  `nClusterings` to better indicate purpose of function. `nClusters`
  now gives the number of clusters per clustering. - `addClusters` to
  `addClusterings` and `removeClusters` to `removeClusterings`. New
  function `removeClusters` allows the user to actually ``remove" a
  cluster or clusters from a clustering by assigning samples in those
  clusters to `-1` value. - `clusterInfo()` to `clusteringInfo()`

- In addition these structural changes, the following enhancements are
  also included in this release - New function `plotClusterLegend` that
  will plot a legend for a clustering. - Color definition changes:
  `showBigPalette` has been replaced with `showPalette` and now can
  show any palette of colors. Adjusted color definitions of `seqPal2`
  and `seqPal4` to be completely symmetric around center. The colors in
  `bigPalette` have been changed and shuffled to reduce similar colors
  and `massivePalette` has been created by adding all of the non-grey
  colors (in random order) from `colors()` so that `plotClusters` will
  not run out of colors. - `getClusterManyParams`: now uses saved
  `clusterInfo` rather than more fragile `clusterLabels` to get
  parameters. The resulting output is formatted somewhat differently. -
  `ClusterExperiment`: removed `transformation` as a required argument.
  Now sets with default of `function(x){x}`. Allows argument
  `clusterLegend` to define the clusterLegend slot in the constructor.
  - `plotClustersWorkflow`: Argument `existingColors` in now takes
  arguments `ignore`,`all`,`highlightOnly` similar to `plotClusters` -
  `plotDendrogram`: Argument `nodeColors` now available.  Changed
  defaults so default is to do colorblock of samples. -
  `plotContrastHeatmap`: Argument `contrastColors` now available to
  assign colors to the contrasts. Genes are now ordered by fold-change
  within each contrast. - `plotClusters`: argument `existingColors` now
  allows for the option `firstOnly` - `makeDendrogram`: now allows
  option 'coCluster' to the argument `reduceMethod` indicating use of
  the coClustering matrix to build the dendrogram. `makeDendrogram` now
  also has a method for building a dendrogram from an arbitrary
  distance function - `clusterMatrix`: now returns cluster matrix with
  rownames corresponding to sample names. - `convertClusterLegend`: now
  takes argument `whichClusters`

Bugs

- converted automatic assignment of colors in `clusterLegend` to be
  based on `massivePalette` so won't run out on toy examples.

- fixed minor bugs in `plotHeatmap` so that will - handle factor with
  only one value in annotation - will plot annotation labels when there
  is `NA` in the annotation - no longer calls internal function
  `NMF:::vplayout` in making those labels, more robust

- fixed bug in how `plotClustersWorkflow` handled existing colors.

- Fixed so `diss` now passed to subsampling in calls to
  clusterSingle/clusterMany

- Fixed so `plotClusters` now will not give incomprehensible error if
  given duplicates of a color

- Fixed `plotDendrogram` so will not create blank plot.

[ClusterSignificance](https://bioconductor.org/packages/ClusterSignificance)
-------------------

Changes in version 1.7.3:

- Adding warn argument to pcp. Mostly for internal use. Fixes the
  warning with non default df argument repeating each perumtation.

- Minor plotting improvments including providing a fix for the plot
  margin error.

Changes in version 1.7.2:

- Tidying and fixing bugs in README. README now from README.Rmd to
  avoid similar problems in the future.

- Adding codecov badge to Readme

Changes in version 1.7.1:

- Updating author info and package URL.

[CNEr](https://bioconductor.org/packages/CNEr)
----

Changes in version 3.5:

NEW FEATURES

- Add function orgKEGGIds2EntrezIDs to fetch the mapping between KEGG
  IDs and Entrez IDs

- Add function makeAxtTracks

- Add function addAncestorGO

[CNVPanelizer](https://bioconductor.org/packages/CNVPanelizer)
------------

Changes in version 1.11.0:

- Add Shiny GUI

[coMET](https://bioconductor.org/packages/coMET)
-----

Changes in version 1.11.5 (2018-04-11):

- Update vignettes

Changes in version 1.11.4 (2018-03-28):

- Update examples in manual because change of functions

- Put some functions obselete because change of database (such as
  COSMIC, regulatory segment ENSEMBL, ISCA)

- remove R package, ggplot2, ggbio and trackviewer

Changes in version 1.11.3 (2018-02-25):

- Update vignette because of the update of BiocStyle

- Update two data (chipTFtrack.rda and genesGencodetrack.rda) because
  error formating

Changes in version 1.11.1 (2017-07-27):

- Update vignette because of the update of BiocStyle

- Add examples data that was not uploaded correctly the first time

[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap)
--------------

Changes in version 1.17.1:

- `Legend()`: add `by_row` argument to control the arrangement of
  legends if they are put in more than one columns

- `Legend()`: use `textGrob()` if the point symbol is text

- `grid.dendrogram()`: fix a bug that the dendrogram is wrong when
  row/column names have duplicated names.

- `anno_boxplot()`: axis rescaled when outline = FALSE

- `oncoPrint()`: rows are first ordered by total number of mutations
  and then ordered by number of samples that have mutations

- correctly reorder rows

- add `row_gap` argument for list of heatmaps

- `oncoPrint()`: add `j` and `i` as optional argument for `alter_fun`

[CountClust](https://bioconductor.org/packages/CountClust)
----------

Changes in version 1.6.1:

- Release Deprecated the cg_topics() and FitGoMpool() functions.
  Modified the FitGoM() function so that it has the flexibility to run
  multiple runs using the num_trials argument and also returns the BIC
  for each model fit. Modified the compGoM() function so that it can
  take either a topic model object, as well as a list of topic model
  objects. Also deprecated the base R graphics for the Structureplot
  and modified the StructureGGplot() function to take single label
  samples.

Changes in version 1.5.1:

- Release We have added a FitGoMpool() function that automatically
  performs GoM model with multiple starting points and outputs the one
  run with the most optimal BIC. Besides, removed
  switch_axis_position() as a dependency from cowplot as the function
  has been deprecated.

[CrispRVariants](https://bioconductor.org/packages/CrispRVariants)
--------------

Changes in version 1.7.12:

- Change CrisprSe$insertion_sites "idxs" column to read indices for
  reads with multiple insertions

Changes in version 1.7.10:

- CrisprSet method setCigarLabels now allows renaming labels

Changes in version 1.7.9:

- Updates mergeCrisprSets in accordance with new order of counting
  operations

Changes in version 1.7.7:

- Adds min and max for guide bounding box in plot

- Fixes bug in plotFreqHeatmap caused when using "group" with a single
  row count matrix

Changes in version 1.7.5:

- Update tests after changes to allele counting

- Added create.plot argument for plotFreqHeatmap with signature
  CrisprSet.

- Adds an option "style" to plotAlignments for colouring only mismatch
  nucleotides

- Changes to narrowAlignments for PacBio cigar format

- Bug fix in collapsePairs.  Only occurred when running outside of
  readsToTarget.

- Adds "alleles" accessor for relating variant labels to the truncated
  cigar strings

- Return unmergeable alignments instead of raising an error

- Minor code changes to make it easier to run a non-standard counting
  pipeline

- Code from initialisers split into separate files for easier
  readability

Changes in version 1.7.4:

- Added create.plot argument for plotFreqHeatmap with signature
  CrisprSet.

Changes in version 1.7.2:

- Allowing plotting arbitrarily many subsets of aligned regions

Changes in version 1.7.1:

- Reorganising plotAlignments code for allowing plotting subsets of the
  aligned regions.

- Minor change to transcript plot plotVariants to make background white
  not transparent.

[csaw](https://bioconductor.org/packages/csaw)
----

Changes in version 1.13.1:

- Fully removed support for paramList objects.

- Removed support for normalize(), modified default option for se.out=
  in normOffsert().

[cydar](https://bioconductor.org/packages/cydar)
-----

Changes in version 1.3.1:

- Bug fix to interpreSpheres() when making additional plots.

- Switched to custom colour calculation in plotCellIntensity().

[cytofkit](https://bioconductor.org/packages/cytofkit)
--------

Changes in version 1.11.3 (2018-01-19):

MODIFICATIONS

- Added select all buttons to other marker related plots

- Changed scaling and centering options to drop-down menus

- Updated maintainer emails

- Added transparency option to expression level plot points so points
  dont obscure each other too much. Point size can be adjusted to
  improve visibility as well

BUG FIXES

- Fixed handling of global scaling so that expression values don't fall
  outside the limits for the colour palette in the shiny color plots

- Improved the outlier removal code in color plots to be more robust
  for various marker expression patterns

- Shifted a line that compares fixedNum argument against the number of
  events so that it doesn't break when the 'min' mergeMethod sets
  fixedNum to NULL

Changes in version 1.11.2 (2017-12-15):

MODIFICATIONS TO EXPRESSION LEVEL SCATTERPLOT

- Marker selection changed to selectize style, choices "All Markers"
  and "All Markers (scaled)" removed.

- Added checkbox option to scale the legend and dot colours
  locally/globally

- Added checkbox option to scale and center expression values

- Added actionButton to select all markers (will update selectize
  choices to select all)

- Added actionButton to update plot after changing marker selection
  (otherwise will update with each marker added/removed)

BUG FIXES

- Changed limits of scatterplot scale from min-max to min-98th
  percentile to account for rare extreme outliers.

Changes in version 1.11.1 (2017-11-06):

MODIFICATIONS

- Added a select/deselect all checkbox to the sample selection panel
  for ease of selection.

- Added FlowSOM option and options to specify k for FlowSOM and
  Rphenograph in cytofkit_GUI.

- Added a popup dialog to cytofkit_GUI that warns if more than 10,000
  cells are being run with DensVM or isomap

BUG FIXES

- Switched the treetype used for nn2 from bd to kd, as bd was causing
  problems in Phenograph clustering.

[dada2](https://bioconductor.org/packages/dada2)
-----

Changes in version 1.7.8:

SIGNIFICANT USER-VISIBLE CHANGES

- nbases has replaced the nreads parameter in the learnErrors function.
  As suggested by the name, this controls the amount of data the
  machine learning uses by the total number of bases rather than the
  read count, which is more appropriate given the range of read-lengths
  in target applications.

- OMEGA_C has been set to 1e-40 by default. This means that
  error-correction is no longer performed on all reads, but instead
  just those a post-hoc probability less than OMEGA_C=1e-40. In
  practice this has a very small impact on final abundances.

BUG FIXES

- The memory usage and speed of assignSpecies on large datasets has
  significantly improved.

Changes in version 1.7.7:

NEW FEATURES

- Error-correction can now be modified, and turned off, by the OMEGA_C
  parameter, which controls the threshold at which reads that are
  inferred to contain errors are corrected (or not) to the sequence
  from which they are inferred to originate.

SIGNIFICANT USER-VISIBLE CHANGES

- The DADA2 options enabling the quick gapless alignment check, and
  extremely conservative greediness in the partioning method were
  turned on by default (GAPLESS=TRUE, GREEDY=TRUE). Some speedup in the
  core denoising algorithm.

- plotQualityProfile now includes a cumulative description of read
  length variation.

Changes in version 1.7.6:

NEW FEATURES

- A new and extremely conservative form of greediness in the core
  denoising algorithm was added, and can be turned on by setting the
  DADA2 option GREEDY=TRUE. This provides some speedup in the core
  denoising algorithm.

Changes in version 1.7.5:

NEW FEATURES

- The dada(...) function now accepts a list of "priors", i.e. sequences
  for which there is prior evidence they might be real. Input sequences
  that match one of the priors are evaluated against a relaxed
  threshold of statistical evidence (OMEGA_P instead of OMEGA_A), and
  can be detected even as singletons.

- The dada(...) function can perform "pseudo-pooling" with dada(...,
  pool="pseudo"). In pseudo-pooling, the input samples are denoised
  independently, then a set of sequences that appear in at least
  MIN_PREVALENCE samples are used as priors for a second and final
  round of sample inference. MIN_PREVALENCE=2 by default.

Changes in version 1.7.4:

NEW FEATURES

- A new fast screen for optimal gapless alignments in the core
  denoising algorithm was added, and can be turned on by setting the
  DADA2 option GAPLESS=TRUE. This provides some speedup in the core
  denoising algorithm.

Changes in version 1.7.3:

BUG FIXES

- Fixed an overflow bug on sequences 260nts or longer in the SSE=2
  code.

Changes in version 1.7.2:

SIGNIFICANT USER-VISIBLE CHANGES

- The DADA2 option enabling explicit 8-bit SSE vectorization in the C
  code was turned on by default (SSE=2). Some speedup in the core
  denoising algorithm.

Changes in version 1.7.1:

NEW FEATURES

- seqComplexity calculates the complexity of input sequences, and can
  be used to identify and filter out low-complexity sequences.

SIGNIFICANT USER-VISIBLE CHANGES

- The default minOverlap parameter of mergePairs was reduced from 20 to
  12, and the alignment parameters used during merging were altered to
  more strongly penalize mismatches and gaps, which improves merging
  performance in repetitive sequences.

[dagLogo](https://bioconductor.org/packages/dagLogo)
-------

Changes in version 1.17.2:

- Fix the bug for formatSequence.R when seq is not character.

Changes in version 1.17.1:

- add .globals environment.

[DaMiRseq](https://bioconductor.org/packages/DaMiRseq)
--------

Changes in version 1.3.7:

- DaMiRseq performs multi-class classification anlysis.

- The Stacking meta-learner can be composed by the user, setting the
  new parameter 'cl_type' of the DaMiR.EnsembleLearning() function. Any
  combination of the 8 classifiers is now allowed.

- If the dataset is imbalanced, a 'Down-Sampling' strategy is
  automatically applied.

- The DaMiR.FSelect() function has the new argument, called 'nPlsIter',
  which allows the user to have a more robust features set. In fact,
  several feature sets are generated by the bve_pls() fuction (embedded
  in DaMiR.FSelect()), setting 'nPLSIter' parameter greater than 1.
  Finally, an intersection among all the feature sets is performed to
  return those features which constantly occur in all runs. However, by
  default, 'nPlsIter = 1'.

- DaMiR.Allplot() accepts also 'matrix' objects as well as NA values
  (which are not plotted).

- The DaMiR.normalization() function estimates the dispersion, through
  the parameter 'nFitType'; as in DESeq2 package, the argument can be
  'parametric' (default), 'local' and 'mean'.

- In the DaMiR.normalization() function, the gene filtering is desabled
  if 'minCount = 0'.

- In the DaMiR.EnsembleLearning() function, the method for implementing
  the Logistic Regression has been changed to allow multi-class
  comparisons; instead of the native 'lm' function, 'bayesglm' method
  implemented in the caret 'train' function, properly set, is now used.

- The new parameter 'second.var' of the DaMiR.SV() function, allows the
  user to take into account a secondary variable of interest (factorial
  or numerical) that the user does not wish to correct for, during the
  sv identification.

[DAPAR](https://bioconductor.org/packages/DAPAR)
-----

Changes in version 1.11.15:

- See changes in version 1.11.13 of the package Prostar

[ddPCRclust](https://bioconductor.org/packages/ddPCRclust)
----------

Version: 0.99.0
Text:

[debrowser](https://bioconductor.org/packages/debrowser)
---------

Changes in version 1.6.6:

- HTML Header fixed

- A more generic file import tool that senses string columns developed

[decontam](https://bioconductor.org/packages/decontam)
--------

Changes in version 0.20.0:

NEW FEATURES

- The decontam R package is released! Read more in our preprint:
  https://www.biorxiv.org/content/early/2017/11/17/221499

[DeepBlueR](https://bioconductor.org/packages/DeepBlueR)
---------

Version: 1.4.1
Text:

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

[DEFormats](https://bioconductor.org/packages/DEFormats)
---------

Changes in version 1.8.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Preserve rownames in gene annotation when converting between
  'DGEList' and 'DESeqDataSet'

[DEGreport](https://bioconductor.org/packages/DEGreport)
---------

Version: 1.15.4
Text: 2017-03-27 Lorena Pantano <lorena.pantano@gmail.com> Fix: Fix
        typo in variable inside degClean Fix: Remove all columsn with
        NA values in degClean Feature: Plot only when degPatterns has
        only one gene. Thanks Amir Jassim.  Feature: Add geom_cor to
        plot correlation values to a ggplot2 plot.  Feature: Add
        eachStep option to degPattern to apply groupDifference to each
        time point and not only to the maximum and minimum values.
        Feature: Add covariates dendograme to degCovariates.  Fix:
        Wrong matrix in degPattern. Thanks Amir Jassim. Feature: Add
        option to filter genes in degPattern. Thanks Amir Jassim.
        Feature: Return raw and summarise table in degPattern Feature:
        Migrate to rmarkdown for vignette Feature: Return prcomp output
        when using degPCA Fix: Typo in degPattern function, and set up
        to FALSE the use of consensusCluster.  Fix: degPlot to be able
        to work with one gene.  Feature: Add the option to look for
        specific patterns, or genes as reference.  Feature: Return
        scaled values if scale==TRUE in degPattern.  Feature: Add
        values used in plots for degPattern function. Thanks to Amir
        Jassim.  Feature: Get significants for a list of DEGSet objects
        binding the tables first, calculating a new FDR, and aplying
        the filter as last step.
        https://support.bioconductor.org/p/104059/#104072

Version: 1.15.2
Text: 2017-01-08 Lorena Pantano <lorena.pantano@gmail.com> Feature: Add
        support to list for significant and recover full table.
        Feature: Add support to different shrinkage estimator.  Fix:
        Volcano plot was plotting wrong the shadows in the y-axis.
        Fix: Use correct option in DESeq2::results to count UP/DOWN
        genes.  Feature: Allow to ask for up/down genes. Thanks to
        Radhika Khetani.

Version: 1.15.1
Text: 2017-11-13 Lorena Pantano <lorena.pantano@gmail.com> Fix: Add
        checking point in degPCA Feature: Add function to plot basic
        expression signatures.

[DEP](https://bioconductor.org/packages/DEP)
---

Changes in version 1.1.5:

- Enabled differential testing with missing values.

- Added data.frame output functionality from plotting functions.

- Enabled multi-SE plotting in plot_normalization() and
  plot_imputation().

- Deprecated se2msn() and msn2se() functions. Use as('MSnSet') or
  as('SummarizedExperiment') from MSnbase package instead.

- Added a vignette on missing value handling in DEP.

- Added citation: Zhang, Smits, van Tilburg et al. Nature Protocols
  2018.

[derfinder](https://bioconductor.org/packages/derfinder)
---------

Changes in version 1.13.8:

BUG FIXES

- Fixed a unit test that was breaking version 1.13.7.

Changes in version 1.13.1:

NEW FEATURES

- Added an extra example to regionMatrix() in response to
  https://support.bioconductor.org/p/103591

[derfinderPlot](https://bioconductor.org/packages/derfinderPlot)
-------------

Changes in version 1.13.5:

BUG FIXES

- Fixed an issue in plotCluster() on how it was loading the
  hg19IdeogramCyto object from the biovizBase package.

Changes in version 1.13.4:

BUG FIXES

- Fixed an issue with a call to GenomicRanges::gaps() that affected how
  the introns were plotting in plotRegionCoverage() when the underlying
  data has a specifying start and end of the chromosome (that is, a
  seqinfo() with seqlengths specified). Thanks to Emily Burke for
  reporting this issue https://github.com/emilyburke.

[DESeq2](https://bioconductor.org/packages/DESeq2)
------

Changes in version 1.20.0:

- Added 'lfcThreshold' argument to lfcShrink() for use with
  type="normal" and type="apeglm". For the latter, lfcShrink() will
  compute FSOS s-values, for bounding when the LFC will be "false sign
  or small", where small is defined by lfcThreshold.

- Switching to a ~10x faster apeglm implementation for use in the
  lfcShrink() function.

- Beginning the deprecation of exploratory analysis of designs without
  replicates. Analysis of designs without replicates will be removed in
  the Oct 2018 release: DESeq2 v1.22.0, after which DESeq2 will give an
  error.

- Elevate 'minmu' to DESeq() as this proves useful for single cell
  applications and certain zero-inflated data.

- Elevate 'useT' to DESeq(), which will use (n - p) for the degrees of
  freedom of the t distribution, and if weights are provided, it will
  use the sum of weights as 'n'.

[DEsingle](https://bioconductor.org/packages/DEsingle)
--------

Changes in version 0.99.12:

- Documentation improvements.

Changes in version 0.99.9:

- Add Parallelization.

Changes in version 0.99.0:

- Package released.

[DiffBind](https://bioconductor.org/packages/DiffBind)
--------

Changes in version 2.8.0:

- Features * dba.report: add precision option to dba.report *
  dba.plot*: make plots use same precision as reports when thresholding
  * dba: add dir option to dba

- Documentation updates * Vignette: change default to
  bFullLibrarySize=TRUE in description of DESeq2 analysis * Vignette:
  update vignette to not change dir * dba.report: clean up description
  of bCalled in man page * dba.report: modify example inman page to be
  clearer * dba: update man page to not change dir * dba.save: dontrun
  example code for dba.save writing into LIB * dba: dontrun example
  code for dba setting wd to LIB

- Bug fixes * dba.report: Sort report by FDR instead of p_value *
  dba.peakset: fix bug when adding consensus peaks with chromosomes not
  present in some peaksets * dba.count: fix bug when recentering
  passed-in peaks * dba.report: fix bug when using filtering from
  dba.analyze * dba.count: fix bug when passing in peaks using factors
  * dba.count: fix bug caused by not registering one of the C routines
  correctly.

[diffcoexp](https://bioconductor.org/packages/diffcoexp)
---------

Changes in version 0.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- example data exprs.1 and exprs.2 are represented as matrices

BUG FIXES

- this package imports rather than depends on the following packages:
  stats, DiffCorr, psych, igraph, BiocGenerics

- messages are generated using message() instead of print() function.

- use is() function to test inheritance relationships between an object
  and a class.

- format NEWS file so that utils::news() parses the file.

Changes in version 0.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- diffcoexp(), coexpr(), and comparecor() accept SummarizedExperiment
  objects.

Changes in version 0.99.0:

SIGNIFICANT USER-VISIBLE CHANGES

- this package was given version number 0.99.0 and submitted to
  Bioconductor.

[diffcyt](https://bioconductor.org/packages/diffcyt)
-------

Version: 1.0.0
Text:

[diffHic](https://bioconductor.org/packages/diffHic)
-------

Changes in version 1.11.8:

- Extended prunePairs() to acknowledge restrict, discard and cap in
  param= argument.

- Extended getPairs() to acknowledge restrict, discard and cap in
  param= argument.

- Added restrict.regions= option to connectCounts(), squareCounts().

- Removed unnecessary normalize() export.

- Upgraded presplit_map.py, iter_map.py to run on Python 3 and to use
  Bio.SeqIO.parse().

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

[DMRcaller](https://bioconductor.org/packages/DMRcaller)
---------

Version: 1.11.1
Text: Added functions detect DMRs from biological replicates (beta
        regression) Added functions to compute/plot spatial correlation
        of methylation level Added function to extract GC content
        (useful for NOME-seq)

[dmrseq](https://bioconductor.org/packages/dmrseq)
------

Changes in version 0.99.11 (2018-04-05):

- dmrseq now provides support for detecting large-scale methylation
  blocks. To use this feature, specify `block=TRUE`, and increase the
  smoothing span parameters `minInSpan`, `bpSpan`, and `maxGapSmooth`.
  More details are provided in the documentation and vignette.

Changes in version 0.99.8 (2018-03-21):

- dmrseq no longer requires balanced, two-group comparisons. To run
  using a continuous or categorial covariate with more than two groups,
  simply pass in the name of a column in `pData` that contains this
  covariate. A continuous covariate is assmued if the data type in the
  `testCovariate` slot is continuous, with the exception of if there
  are only two unique values (then a two group comparison is carried
  out).

Changes in version 0.99.6 (2018-03-02):

- dmrseq is now available on Bioconductor devel (3.7) !

[DominoEffect](https://bioconductor.org/packages/DominoEffect)
------------

Changes in version 0.99.5:

NEW FEATURES

- add code to fail gracefully if uniprot is down

Changes in version 0.99.4:

NEW FEATURES

- add code to fail gracefully if ensembl is down

Changes in version 0.99.3:

NEW FEATURES

- add Citation of publication

Changes in version 0.99.2:

NEW FEATURES

- address comments from Bioconductor review

Changes in version 0.99.1:

NEW FEATURES

- new functions import_vcf.R and import_txdb.R to use standard
  Bioconductor classes

- minor other updates

Changes in version 0.99.0:

NEW FEATURES

- First Version of DominoEffect submitted to Bioconductor

[DOSE](https://bioconductor.org/packages/DOSE)
----

Changes in version 3.5.2:

- bug fixed of gseaScores <2018-04-18, Wed> +
  https://github.com/GuangchuangYu/DOSE/issues/23

- mv web site to https://guangchuangyu.github.io/software/DOSE

[drawProteins](https://bioconductor.org/packages/drawProteins)
------------

Version: 0.99.10
Category: ARGUEMENT
Text: show.legend added as argument to allow the option to show or not
        show the legend as per ggplot2.

Version: 0.99.9
Category: FEATURES
Text: New function called draw_recept_dom(). This function allows the
        drawing of the TOPO_DOM and TRANSMEM types of receptors. Data
        from TNFR1 and CD40 are included to demonstate the function.

Version: 0.99.8
Category: FEATURES
Text: New function called extact_transcripts. This function will ammend
        the data frame to allow each chain from the same UniProt
        accession number to to drawn separately. A vignette entitled
        drawProteins_extract_transcripts has been written to
        demonstrate.

Version: 0.99.8
Category: FEATURES
Text: LazyData is now false and NAMESPACE updated as per Bioconductor
        review.

Version: 0.98.3
Category: FEATURES
Text: New function called draw_canvas. This function was previously
        within draw_chains but has now been pulled out to allow the
        generation of a canvas separately from the chains. It did
        require quite a rewrite but I think it will make things more
        useful For example, it will allow the plotting of domains
        without chains which has the potential to be very useful.

Version: 0.98.2
Category: FEATURES
Text: Rename functions from geom to draw. E.g geom_chains is now
        geom_draw. This is because they weren't really geoms and using
        the word draw seem more helpful and a better reflection of the
        function.

Version: 0.98.2
Category: LAUNCH VERSION 0.98.1
Text:

Version: 0.98.2
Category: FEATURES
Text: Drawing protein schematics from Uniprot database with Accession
        numbers

[DropletUtils](https://bioconductor.org/packages/DropletUtils)
------------

Changes in version 0.99.0:

- New package DropletUtils, for handling droplet-based single-cell RNA
  sequencing data.

[easyRNASeq](https://bioconductor.org/packages/easyRNASeq)
----------

Changes in version 2.15.5:

- Ensured that counting happens in a strand specific way when stranded
  data is provided

- The easyRNASeq and all related methods are now defunct.

- Added a BamParam strandProtocol argument value to count reads on the
  reverse strand

- Removed calls to RangedData constructor defunct parameters in unit
  tests

- Removed dependencies to RnaSeqTutorial in unit tests

- Removed the easyRNASeq vignette. Replaced it with a knitr vignette -
  to be completed.

Changes in version 2.15.4:

- Removed dependencies to RnaSeqTutorial

- Removed calls to RangedData constructor defunct parameters

[EBImage](https://bioconductor.org/packages/EBImage)
-------

Changes in version 4.22.0:

BUG FIXES

- fixed compilation errors on Solaris

[edgeR](https://bioconductor.org/packages/edgeR)
-----

Changes in version 3.22.0:

- New function read10X() to read 10X Genomics files.

- New function nearestTSS() to find the nearest transcriptional start
  site (TSS) for given genomic loci.

- New function nearestReftoX() to find the element of a reference table
  that is closest to each element of an incoming vector.

- New function modelMatrixMeth() to construct design matrices for
  analysis of methylation data.

- New function filterByExpr() to filter low expression genes or
  features.

- New rowsum method for DGEList objects.

- nbinomUnitDeviance() now respects vectors.

- DGEList() takes 'group' from 'samples' only if samples has a column
  called group.

- decideTestsDGE() now includes a 'label' attribute, which allows more
  information row.names for the summary results table from
  decideTestsDGE() or decideTests().

- Design now defaults to y$design for all the gene set tests.

- More intuitive error messages from glmFit() when the arguments are
  not conformal.

- Update User's Guide to cite the Chen et al (2017) methylation
  workflow.

- Change glmTreat() default to lfc=log2(1.2).

- Fix incorrect implementation of weights in adjustedProfileLik().

- Bug fix to glmLRT() when there is just one gene but multiple
  contrasts.

- Bug fix to cpmByGroup().

[EGSEA](https://bioconductor.org/packages/EGSEA)
-----

Changes in version 1.6.1 (2017-12-03):

- added: exception handling code to avoid generateReport failure when
  generating GO graphs under some circumstances.

- fixed: a minor issue in writing CSV files.

[EnrichedHeatmap](https://bioconductor.org/packages/EnrichedHeatmap)
---------------

Changes in version 1.9.3:

- support visualize category signals

- add `discretize()` to transform continuous matrices to discrete
  matrices

- add one more vignette showing the usage of categorical signals

Changes in version 1.9.2:

- improvement on vignettes

- enriched_scores() directly applied to the normalized matrix

- row dendrogram is reordered by enriched scores if cluster_row is set
  to TRUE

[EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser)
-----------------

Changes in version 2.10.0:

- Adding scripts to inst/scripts to invoke the EnrichmentBrowser from
  the command line (for non-R users)

- GRN compilation: supporting additional pathway databases (via
  graphite)

- Caching for download of GO and KEGG gene sets (via BiocFileCache)

- Default output destination changed to
  rappdirs::user_data_dir("EnrichmentBrowser")

- Function names: deprecation of x.x notation - read.eset -> readSE -
  probe.2.gene.eset -> probe2gene - de.ana -> deAna -
  compile.grn.from.kegg -> compileGRN - ggea.graph -> ggeaGraph -
  make.example.data -> makeExampleData

[ensembldb](https://bioconductor.org/packages/ensembldb)
---------

Changes in version 2.3.14:

- New ProteinDomainSourceFilter and ProteinDomainIdFilter.

Changes in version 2.3.11:

- Fix for issue #75.

Changes in version 2.3.9:

- Add transcriptToCds and cdsToTranscript functions (issue #73).

Changes in version 2.3.8:

- Fix problem creating EnsDb from GTF file lacking exon IDs (issue #72,
  https://support.bioconductor.org/p/105536/).

Changes in version 2.3.7:

- Switch to EnsDb.Hsapiens.v86 in examples, unit tests and vignettes.

Changes in version 2.3.6:

BUG FIXES

- Fix issue #69 failing to map genomic coordinates to proteins if
  different genomic coordinates are mapped to the same
  transcripts/proteins.

Changes in version 2.3.5:

BUG FIXES

- Fix problem mapping genomic coordinates on named GRanges.

Changes in version 2.3.2:

NEW FEATURES

- Functionality to between genomic and transcript-relative or
  protein-relative coordinates: proteinToTranscript,
  genomeToTranscript, transcriptToGenome, transcriptToProtein,
  genomeToProtein.

Changes in version 2.3.1:

NEW FEATURES

- Functionality to map within-protein coordinates to genomic
  coordinates: proteinToGenome function.

[ensemblVEP](https://bioconductor.org/packages/ensemblVEP)
----------

Changes in version 1.22.0:

NEW FEATURES

- add support for Ensembl release 92

- add support for Ensembl release 91

MODIFICATIONS

- The default 'host' is not specified, defaulting to vep default of
  'ensembldb.ensembl.org'. The default previously was
  'useastdb.ensembl.org' Users in the US may find connection and
  transfer speeds quicker using the East coast mirror,
  'useastdb.ensembl.org'. It was updated because the mirror only
  supports current and current minus one vep version, and to bring the
  default in line with the default of vep.

[epiNEM](https://bioconductor.org/packages/epiNEM)
------

Version: 2017.01
Category: Github made public and submitted to Bioconductor
Text:

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

[erccdashboard](https://bioconductor.org/packages/erccdashboard)
-------------

Changes in version 1.13.1:

SIGNIFICANT USER-VISIBLE CHANGES

- replaced QuasiSeq with edgeR for differential expression testing, due
  to deprecation of QuasiSeq

- users will observe significant speed increases for DE testing

- users will observe changes in results from example data,
  Quasi-Likelihood methods are similar for edgeR and QuasiSeq, but not
  exactly identical, so P-value distributions are different, LODR
  estimates have changed.

- dispersion plot (dispPlot) is now generated via edgeR instead of
  QuasiSeq, using a similar Quasi-Likelihood dispersion shrinkage
  method.

BUG FIXES AND MINOR IMPROVEMENTS

- changed behavior so that dispPlot is not deleted if DE testing is
  rerun, instead increment file name to prevent overwriting.

[esATAC](https://bioconductor.org/packages/esATAC)
------

Version: 1.2.0
Category: Publication
Text:

Version: 1.2.0
Category: Wei, Zheng, et al. "esATAC: an easy-to-use systematic
        pipeline for ATAC-seq data analysis." Bioinformatics (2018
Text:

Version: 1.2.0
Category: https://doi.org/10.1093/bioinformatics/bty141
Text:

Version: 1.2.0
Category: Fix some known bugs
Text:

Version: 1.2.0
Category: Speed up motif scanning
Text:

Version: 2018.04.18
Category: Using the new method to scan motif in the genome, remove
        preset motif data. The new method
Text:

Version: 2018.04.18
Category: is from package "motifmatchr", it cost less than 2mins for a
        sample(test sample: SRR891268;
Text:

Version: 2018.04.18
Category: of course, depend on your computer performance
Text:

[FamAgg](https://bioconductor.org/packages/FamAgg)
------

Changes in version 1.7.4:

- Fix trait-related problem in binomial test.

Changes in version 1.7.3:

- Deprecating the probability test because of problems with the gap
  package in MS Windows.

[fCCAC](https://bioconductor.org/packages/fCCAC)
-----

Changes in version 1.5.1:

- Updated contact email.

[FELLA](https://bioconductor.org/packages/FELLA)
-----

Version: 0.99.6
Text: Small corrections in vignette

Version: 0.99.6
Text: Removed `rm` calls

Version: 0.99.5
Text: Small correction in vignette

Version: 0.99.5
Text: Version bump to see if bioc build is bugged

Version: 0.99.4
Text: Moved contents of `NEWS.md` to NEWS

Version: 0.99.4
Text: Deleted most of `data-raw` folder, moved the rest to
        `inst/script`

Version: 0.99.4
Text: Removed redundant `Authors` field in DESCRIPTION

Version: 0.99.4
Text: Removed `class(data) != "FELLA.DATA"` by using built-in
        `is.FELLA.DATA`

Version: 0.99.3
Text: Version bump (biomaRt down?)

Version: 0.99.2
Text: Fixed more doc links

Version: 0.99.1
Text: Fixed doc links (hopefully)

Version: 0.99.1
Text: Updated funding in vignette

Version: 0.99.1
Text: Travis will only test the devel branch

Version: 0.99.0
Text: Submission to Bioconductor

[fgsea](https://bioconductor.org/packages/fgsea)
-----

Changes in version 1.5.2:

- fgsea throws warning for rank ties

- fgsea throws warning for duplicate gene names

- Leading edge is now ordered by decreasing of absolute statistic value

Changes in version 1.5.1:

- Reproducibility fixes

- Added collapsePathway function to intelligently collapse overlapping
  pathways

- Added fgseaLabel function for label-permuting GSEA

[flowcatchR](https://bioconductor.org/packages/flowcatchR)
----------

Changes in version 1.14.0:

NEW FEATURES

- the engine for plotting the trajectories is now based on plotly,
  abandoning the older rgl-based solution

OTHER NOTES

- vignette: now also based on Rmarkdown, leaving the older Rnw-pdf be

- shortened example lines for better documentation

- better definition of import statements for better behaviour with
  other packages

- updated travis.yml for current setup

[flowPloidy](https://bioconductor.org/packages/flowPloidy)
----------

Changes in version 1.5.6 (2018-03-08):

User Visible Changes

- BUG FIX: the G2 peak of the B sample was not getting incorporated
  into model construction, which caused model fitting to fail on
  samples with histograms skewed towards the left.

- Updated DebrisModel documentation.

Internal Changes

- Fixed broken test.

Changes in version 1.5.3 (2018-01-17):

User Visible Changes

- The browser interface will no longer allow users to enter arbitrary
  text to select the number of samples. Only valid values, i.e., 1, 2,
  or 3, will be offered as choices. The old version would crash with a
  bad value for sample number.

- Improved file importing, so that low-quality samples for which peaks
  cannot be detected can still be imported.

- Added a new argument to FlowHist, so users can set the threshold
  below which the data is ignored when screening out debris. The
  default remains the same, 40, but users with very clean histograms,
  with peaks far to the left, can now lower this value if needed.

- Model fitting is now limited to the range of the data. That is, empty
  bins at the right (upper) end of the histogram will not be included
  when fitting the NLS. This addresses problems that generated inflated
  RCS values.

- The limits for the linearity parameter have been extended from
  1.9-2.1, to 1.5-2.5. This will improve/enable model fitting on
  samples where the linearity (the ratio of the G2/G1 peaks) was
  outside the range 1.9-2.1.

Changes in version 1.5.1 (2017-12-04):

Minor bug squashed

- Changes made upstream to `car::deltaMethod()` introduced a bug. This
  will be resolved as of `car` version 2.1-7. Until that version of
  `car` makes it's way into CRAN, calling `deltaMethod()` on an `nls`
  object will require the argument `vcov.` be explicitly set to `vcov`,
  a function in the stats package. This has been done in `flowPloidy`,
  and should be invisible to users. In addition, the bug is not present
  in the current R release, so it should not need to be backported to
  previous versions of flowPloidy.

[gage](https://bioconductor.org/packages/gage)
----

Changes in version 2.28.2:

- korg now include 5245 KEGG species. In addition an updated version of
  korg is now checked out from Pathview Web server by loading .rda file
  (inst4ead of read .tsv file) when needed.

[GARS](https://bioconductor.org/packages/GARS)
----

Changes in version 0.99.0:

- GARS package has been released!


[GDSArray](https://bioconductor.org/packages/GDSArray)
--------

Changes in version 0.99.11:

- example(GDSArray) added.

Changes in version 0.99.10:

- renamed "gdsNodes" method as "gdsnodes".

Changes in version 0.99.7:

- renamed "GDSlight" as "GDSFile".

- added accessors for "GDSFile".

- $ method check for valid gds node.

Changes in version 0.99.6:

- GDSArray seed setter and getter.

Changes in version 0.99.5:

- Implemented "GDSlight" class and $ method.

Changes in version 0.99.1:

- Added example for GDSArray-methods.Rd.

- Added BiocViews.

[gdsfmt](https://bioconductor.org/packages/gdsfmt)
------

Changes in version 1.16.0:

UTILITIES

- a new storage name 'single' in `add.gdsn()` for single-precision
  floating numbers

- improve the efficiency of bit2-unpacking when there are lots of zero

- `system.gds()` outputs 'POPCNT' flag if available

- enable the compression modes "LZMA.ultra", "LZMA.ultra_max",
  "LZMA_RA.ultra" and "LZMA_RA.ultra_max"

- show more compression information in `system.gds()`

BUG FIXES

- avoid the integer overflow when the compression rate is too small
  using LZMA_RA (e.g., <0.01%)

Changes in version 1.14.0-1.14.1:

UTILITIES

- tweak error messages in `apply.gdsn()`

- `cleanup.gds()` allows a file name with a prefix '~' which will be
  automatically replaced by the home directory

[GeneNetworkBuilder](https://bioconductor.org/packages/GeneNetworkBuilder)
------------------

Changes in version 1.21.2:

- Fix the problem for missing links in documentation.

Changes in version 1.21.1:

- Fix the problem for missing links in documentation.

[GENESIS](https://bioconductor.org/packages/GENESIS)
-------

Changes in version 2.9.3:

- New methods assocTestSingle and assocTestAggregate are refactors of
  assocTestMM and assocTestSeq/assocTestSeqWindow, respectively.
  assocTestSeq and assocTestSeqWindow are deprecated. assocTestMM is
  still used for GenotypeData objects, but will be deprecated in a
  future release. fitNullModel is a refactor of fitNullMM/fitNullReg
  and should be used with the new association test methods.

[GeneStructureTools](https://bioconductor.org/packages/GeneStructureTools)
------------------

Version: 0.99.0
Category: Pre-BioC submission
Text:

Version: 0.99.1
Category: Fixes for BioC submission
Text: notNMD is not required for vignettes

Version: 0.99.1
Category: Fixes for BioC submission
Text: whippet files are read in as a single object

Version: 0.99.1
Category: Fixes for BioC submission
Text: Fixed bug where leafcutter set creation snowballed

Version: 0.99.2
Text: Fix for fread()'ing in gzip files on windows

Version: 0.99.4
Text: Exon skipping and intron retention now works with manual
        coordinates again

Version: 0.99.4
Text: Documentation updates

[genomation](https://bioconductor.org/packages/genomation)
----------

Version: 1.11.3
Category: IMPROVEMENTS AND BUG FIXES
Text: fixed an error of data.table() function that occurs when
        GrangesList object contains unnamed windows

Version: 1.11.3
Category: in ScoreMatrixBin() function
Text:

Version: 1.11.2
Category: IMPROVEMENTS AND BUG FIXES
Text: improved the C++ function Median_c() to handle NAs

Version: 1.11.2
Category: IMPROVEMENTS AND BUG FIXES
Text: the following C++ functions return NAs instead of zeros if the
        length of the vector is smaller than the number of bins:
        binMean(), binMedian(), binMax(), binMin(), binSum()

Version: 1.11.2
Category: NEW FUNCTIONS AND FEATURES
Text:

Version: 1.11.2
Category: C++ functions that compute a desired value from a vector and
        handle NAs
Text: Mean_c() - computes a mean value,

Version: 1.11.2
Category: C++ functions that compute a desired value from a vector and
        handle NAs
Text: Max_c() - computes a maximum values,

Version: 1.11.2
Category: C++ functions that compute a desired value from a vector and
        handle NAs
Text: Min_c() - computes a minumum values,

Version: 1.11.2
Category: C++ functions that compute a desired value from a vector and
        handle NAs
Text: Sum_c() - computes a sum value.

Version: 1.11.1
Category: IMPROVEMENTS AND BUG FIXES
Text: bug fix relating to ScorematrixBin that returns all 1's when
        is.noCovNA=T
        (https://github.com/BIMSBbioinfo/genomation/issues/168)

Version: 1.11.1
Category: IMPROVEMENTS AND BUG FIXES
Text: xcoords argument for heatMatrix and multiHeatMatrix now can take
        character vectors.

Version: 1.11.1
Category: The character vectors will label the x-axis of heatmaps.
        Examples: xcoords=c("-2kb","0","2kb
Text:

[GenomicAlignments](https://bioconductor.org/packages/GenomicAlignments)
-----------------

Changes in version 1.16.0:

NEW FEATURES

- Add coercion from list to GAlignmentsList.

SIGNIFICANT USER-VISIBLE CHANGES

- Improve performance of [[<- on GAlignmentsList objects. This is a
  100x speedup or more on a big GAlignmentsList object.

BUG FIXES

- Remove spurious warning in summarizeOverlaps().

[GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures)
---------------

Changes in version 1.32.0:

NEW FEATURES

- The first argument of mapToTranscripts() and pmapToTranscripts() now
  can be a GPos object and a GPos object is returned in that case.

- Add 'use.names' argument to "transcripts", "exons", "cds, and
  "promoters" methods for TxDb objects.

- makeTxDbFromUCSC() now uses direct SQL queries (to the UCSC MySQL
  server at genome-mysql.soe.ucsc.edu) instead of
  rtracklayer::getTable() to fetch data from the Genome Browser. This
  avoids the issue reported here
  https://github.com/lawremi/rtracklayer/issues/5 . Another benefit is
  that direct SQL queries are much faster than rtracklayer::getTable().

SIGNIFICANT USER-VISIBLE CHANGES

- The GRanges object returned by mapToTranscripts() or
  pmapToTranscripts() takes the transcript lengths as seqlengths.

- pmapToTranscripts() always takes the transcript name as the seqname,
  even when there is no overlap. Before, it used "UNMAPPED" as the
  seqname when there was no overlap.

BUG FIXES

- Fix bug where 'coverageByTranscript(x, transcripts)' was erroring in
  situations where 'transcripts' contains transcripts with an exon that
  receives no coverage and is located on a sequence for which the
  seqlength is not available in 'seqinfo(x)' nor in
  'seqinfo(transcripts)'.

[GenomicRanges](https://bioconductor.org/packages/GenomicRanges)
-------------

Changes in version 1.32.0:

NEW FEATURES

- 2 improvements to the "promoters" method for GenomicRanges objects: -
  The 'upstream' and 'downstream' arguments now can be integer vectors
  parallel to 'x', - The 'use.names' argument now is supported. This is
  for consistency with the other intra range transformations.

SIGNIFICANT USER-VISIBLE CHANGES

- GenomicRanges now is a List subclass. This means that GRanges objects
  and their derivatives are now considered list-like objects (even
  though [[ don't work on them yet, this will be implemented in
  Bioconductor 3.8).

- Add the CompressedGRangesList class as a replacement for the
  GRangesList class. The long term goal is that GRangesList becomes a
  virtual class with CompressedGRangesList as a concrete subclass. Note
  that the GRangesList() constructor now returns a
  CompressedGRangesList instance instead of a GRangesList instance.

- GenomicRangesList is now a virtual class (like IntegerRangesList is).

- GRanges derivatives no longer support the 'x&#91;i, j&#93; <- value' form of
  subassignment. This feature was of very limited usefulness and no
  Bioconductor package was using it.

- Improve performance of nearest(), precede(), and follow() on a
  GRanges object.

- Improve performance of coverage() on a GPos object.

- Improve performance of sort() on a GRangesList object. Also now it
  supports 'ignore.strand'. See
  https://github.com/Bioconductor/GenomicRanges/issues/1 (and note how
  unnicely these changes were requested).

- Improve performance and error handling of coercion from RleList to
  GRanges. This is a 50x speedup or more when the RleList object to
  coerce has thousands of list elements or more.

BUG FIXES

- Fix coercion from RleList to GRanges when some list elements in the
  object to coerce have length 0 (see
  https://support.bioconductor.org/p/105926/ for original report by
  Xiaotong Yao).

- Fix bug in nearest() when an unstranded range in 'query' precedes or
  follows more than one range in 'subject'.

[GenomicScores](https://bioconductor.org/packages/GenomicScores)
-------------

Changes in version 1.4.0:

USER VISIBLE CHANGES

- The function 'scores()' has been deprecated and replaced by the
  function 'gscores()'.

- The argument 'scores.only' in the function 'scores()' has been
  deprecated and replaced by calling the function 'score()'.

- The 'MafDb' class has been deprecated and now the 'GScores' class
  supports former 'MafDb' objects. The 'mafByOverlaps()' and
  'mafById()' functions have been deprecated and replaced by the
  function 'gscores()'. The 'populations()' function from the 'MafDb'
  API has been integrated into the 'GScores' API.

- Added metadata on genomic scores groups, available through the
  function 'gscoresGroups()', on availability of non-single nucleotide
  regions through the function 'gscoresNonSNRs()', and on the default
  population used through the function 'defaultPopulation()'.

- New AnnotationHub resources have been added during this release
  cycle: phyloP60way.UCSC.mm10, LINSIGHT, phastCons46wayPlacental,
  phastcons46wayPrimates.

- Added a BiocSticker at
  https://github.com/Bioconductor/BiocStickers/tree/master/GenomicScores

- Added citation information after package publication has been
  accepted at Bioinformatics.

[genphen](https://bioconductor.org/packages/genphen)
-------

Changes in version 1.7 (2018-01-20):

Introduction

- Additional functions are added sporadically.

- This news file reports changes that have been made as the package has
  been developed.

Changes

- Bayesian inference using stan (rstan package)

- Two bayesian inference methods to support both continuous and
  dichotomous phenotypes

- Additional genotype-phenotype association metrics added

- Posterior predictive checks implemented

- Procedure for data reduction (runDiagnostics) added

- Tutorial updated

- Simple procedure for phylogenetic bias estimation implemented

- Retrospective power analysis module

To do

- Add practical examples where genphen has been used.

- Implement modules for data augmentation

- Update todo after data augmentation

[GEOquery](https://bioconductor.org/packages/GEOquery)
--------

Version: 2.47.1
Text: Bug fixes: * Fixes problems with intermittent connection issues
        (which were not intermittend connection problems, it seems)

[ggtree](https://bioconductor.org/packages/ggtree)
------

Changes in version 1.11.3:

- update msaplot to use DNAbin/AAbin internally and also compatible
  with treedata object <2017-12-14, Thu>

- clean up code <2017-12-13, Thu>

- remove paml_rst, codeml_mlc, codeml and jplace fortify methods
  according to the change of treeio (v = 1.3.3) <2017-12-07, Thu>

Changes in version 1.11.2:

- keep tree order (previously using postorder) <2017-12-06-Wed> +
  https://github.com/GuangchuangYu/ggtree/issues/157

- remove beast object support as read.beast output treedata object in
  treeio <2017-12-05, Tue>

- deprecate subview, annotation_image and phylopic; remove
  theme_transparent <2017-12-04, Mon>

- geom_tiplab now supports geom = "image" or geom = "phylopic"
  <2017-12-04, Mon>

- A new layer geom_nodelab that equivalent to geom_tiplab but works for
  internal node <2017-12-04, Mon>

Changes in version 1.11.1:

- bug fixed in geom_tiplab, now `offset` parameter works with
  `align=TRUE`. <2017-11-20, Mon>

- enable mrsd parameter for treedata object <2017-11-15, Wed>

- set_hilight_legend supports alpha parameter <2017-11-15, Wed> +
  https://github.com/GuangchuangYu/ggtree/issues/149

[GOfuncR](https://bioconductor.org/packages/GOfuncR)
-------

Changes in version 0.99.14:

NEW FEATURES

- allow for custom annotations (alternative to default Bioconductor
  annotation packages)

- allow for custom ontology graph (alternative to default integrated
  GO-graph)

Changes in version 0.99.12:

USER-LEVEL CHANGES

- update GO-graph (version 10-Apr-2018)

[GOpro](https://bioconductor.org/packages/GOpro)
-----

Version: 1.5.1
Text:

[graphite](https://bioconductor.org/packages/graphite)
--------

Changes in version 1.25.9 (2018-04-28):

- Updated all pathway data.

Changes in version 1.25.7 (2018-04-18):

- cytoscapePlot returns an handle to the exported graph.

Changes in version 1.25.6 (2018-04-11):

- Parallel versions of clipper and topologyGSA methods.

Changes in version 1.25.3 (2018-03-23):

- Faster conversion of identifiers.

Changes in version 1.25.2 (2018-03-23):

- Documentation fixes.

- Removed deprecated objects biocarta, humancyc, kegg, nci, panther,
  reactome.

Changes in version 1.25.1 (2017-12-13):

- Export the pathwayURL function.

[GreyListChIP](https://bioconductor.org/packages/GreyListChIP)
------------

Changes in version 1.11.1:

- Allow "karyotype" to be specified as a GRanges object.

- Include Anshul Kundaje's black lists for convenience.

- Support merging of grey lists from multiple input files.

[GSVA](https://bioconductor.org/packages/GSVA)
----

Changes in version 1.28:

USER VISIBLE CHANGES

- Arguments 'rnaseq', 'kernel', 'no.bootstraps' and 'bootstrap.percent'
  have become defunct.

- A Bioconductor sticker has been created and it is available at
  https://github.com/Bioconductor/BiocStickers/tree/master/GSVA

[HDTD](https://bioconductor.org/packages/HDTD)
----

Changes in version 1.13.3 (2018-01-10):

- Updated README.

- Updated vignettes.

Changes in version 1.13.2 (2017-12-10):

- Minor changes.

Changes in version 1.13.1 (2017-11-20):

- Updated CITATION FILE.

- Updated documentation.

- Added RcppArmadillo to improve performance.

[HIBAG](https://bioconductor.org/packages/HIBAG)
-----

Changes in version 1.16.0:

- KIR information in `hlaLociInfo()`

- new functions `hlaGenoSubsetFlank()` and `hlaLDMatrix()`

[hipathia](https://bioconductor.org/packages/hipathia)
--------

Changes in version 0.99.27 (2018-04-08):

- Correcting Bionconductor style

Changes in version 0.99.9 (2018-03-09):

- Adapting hipathia package to SummarizedExperiment and
  MultiAssayExperiment.

Changes in version 0.99.8 (2018-03-07):

- First version using AnnotationHub to load annotation files.

Changes in version 0.99.7 (2018-02-14):

- Last version supporting hpAnnot as annotation package

[HiTC](https://bioconductor.org/packages/HiTC)
----

Changes in version 1.23.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Add use.names and header options in exportC

BUG FIXES

- Fix bug in getExpectedCountsMean for non-symmetrical data

- Deal with NA in getPearson function

- Fix bug in normLGF leading to non symmetrical matrices

[hpar](https://bioconductor.org/packages/hpar)
----

Changes in version 1.21.1:

- use html_document (rather that deprecated html_document2) <2018-01-19
  Fri>

[HTSFilter](https://bioconductor.org/packages/HTSFilter)
---------

Version: 1.19.3
Category: —- Fixed bug for the DESeqDataSet class when there are
        additional columns in the
Text:

Version: 1.19.3
Category: colData of the object beyond those used in the design (thanks
        to Stephanie Durand
Text:

Version: 1.19.3
Category: for finding the bug!)
Text: -- Keep rownames and colnames on filtered data for matrix and
        data.frame class

[ideal](https://bioconductor.org/packages/ideal)
-----

Changes in version 1.4.0:

NEW FEATURES

- Specified single go term selection for generating heatmaps of gene
  signatures

- Added support for logFC shrinkage, following the latest devels of
  DESeq2

BUG FIXES

- Corrected output for the vignette, as html_document2 is now
  deprecated

- Menus are back in the expanded form

- Fixed the behavior with addMLE

OTHER NOTES

- Added further progress indicators to give feedback during lengthy
  steps

- Improved ggplotCounts for better scale display, using exact arg
  matching, defaulting to the transformed counts

[InPAS](https://bioconductor.org/packages/InPAS)
-----

Changes in version 1.11.4:

- fix the warning for help files in windows

Changes in version 1.11.3:

- add utr3.danRer10.

Changes in version 1.11.2:

- fix the error by the change of GRanges.

[InTAD](https://bioconductor.org/packages/InTAD)
-----

Changes in version 0.99.1-4:

- Package accepted in Bioconductor

- Adjustment for R version 3.5

- Minor fixes

[IntEREst](https://bioconductor.org/packages/IntEREst)
--------

Changes in version 1.4.0:

NEW FEATURES

- interest() and interest.sequential() functions now support "IntSpan"
  method, allows counting intron- spanning reads.

- psi() function is added. It calculates Psi values.

- annotateU12() fucntion supports DNAStringSet objects as its refGenome
  input.

- Improved vignette document.

SIGNIFICANT USER-VISIBLE CHANGES

- eBayesInterest() removed.

BUG FIXES

- interest.sequential() and interest() corrections to their object
  output option.

- annotateU12() modified to work correctly with the new changes in
  Biostrings package.

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

[IRanges](https://bioconductor.org/packages/IRanges)
-------

Changes in version 2.14.0:

NEW FEATURES

- Add the windows() generic with various methods. This is a "parallel"
  version of window() for list-like objects i.e. it does
  'mendoapply(window, x, start, end, width)' but uses a fast
  implementation.  Also add heads() and tails() as convenience wrappers
  around windows().  They do 'mendoapply(head, x, n)' and
  'mendoapply(tail, x, n)', respectively, but use a fast
  implementation. They're replacements for S4Vectors::phead() and
  S4Vectors::ptail() which are now deprecated.

- Add equisplit() to split a vector-like object into a specified number
  of partitions with equal (total) width. This is useful for instance
  to ensure balanced loading of workers in parallel evaluation.

- promoters() arguments 'upstream' and 'downstream' now can be integer
  vectors parallel to 'x' (for consistency with the other intra range
  transformations).

- The promoters() generic and methods get the 'use.names' argument.

- Add "resize", "flank", and "restrict" methods for Views objects.

- Add "as.integer" method for Pos objects (equivalent to pos()).

SIGNIFICANT USER-VISIBLE CHANGES

- The Ranges virtual class is now the common parent of the IRanges,
  GRanges, and GAlignments classes (GRanges and GAlignments are defined
  in the GenomicRanges and GenomicAlignments packages, respectively).
  More precisely, Ranges is a virtual class that now serves as the
  parent class for any class that represents a vector of ranges. The
  ranges can be integer ranges (i.e. ranges on the space of integers)
  like in an IRanges object, or genomic ranges (i.e. ranges on a
  genome) like in a GRanges object. Note that because Ranges extends
  List, all Ranges derivatives are considered list-like objects. This
  means that GRanges objects and their derivatives are considered
  list-like objects, which is new (even though [[ don't work on them
  yet, this will be implemented in Bioconductor 3.8).

- Similarly the RangesList virtual class is now the common parent of
  the IRangesList, GRangesList, and GAlignmentsList classes.

- IRanges objects don't support [[, unlist(), as.list(), lapply(), and
  as.integer() anymore. This is a temporary situation only. These
  operations will be re-introduced in Bioconductor 3.8 but with a
  different semantic. The overall goal of all these changes is to bring
  more consitency between IRanges and GRanges objects (GRanges objects
  will also support [[, unlist(), as.list(), and lapply() in
  Bioconductor 3.8).  Non-exported IRanges:::unlist_as_integer() helper
  is a temporary replacement for what unlist() and as.integer() used to
  do a IRanges object.

- Move the pos() generic to BiocGenerics.

- Switch order of breakInChunks() arguments 'chunksize' and 'nchunk' to
  be consistent with tileGenome().

- tile() and slidingWindows() now preserve names.

- Optimize [[<- on a CompressedList object. Was very inefficient. The
  optimized method can be up to 100x faster or more on a long object.

- All the S4Vectors-specific material in the IRangesOverview.Rnw
  vignette has moved to the new S4VectorsOverview.Rnw vignette located
  in the S4Vectors package.

DEPRECATED AND DEFUNCT

- Deprecate the RangesList() constructor. IRangesList() should be used
  instead.

- The "ranges" methods for Hits and HitsList objects are now defunct
  (were deprecated in BioC 3.6).

- The "overlapsAny", "subsetByOverlaps", "coverage" and "range" methods
  for RangedData objects are now defunct (were deprecated in BioC 3.6).

- The universe() getter and setter as well as the 'universe' argument
  of the RangesList(), IRangesList(), RleViewsList(), and RangedData()
  constructor functions are now defunct (were deprecated in BioC 3.6).

[iSEE](https://bioconductor.org/packages/iSEE)
----

Changes in version 0.99.3:

- Custom tours can be restarted via the dropdown menu button,
  overwriting the default tour

- Added functionality to provide a custom title to be displayed in the
  app

- Preserve data points and width ratio upon zoom on discrete variables

Changes in version 0.99.2:

- Added functionality for providing additional custom tours, to be
  launched directly upon starting the app.

Changes in version 0.99.1:

- Added grid-based visual point downsampling for faster plotting,
  including control of resolution.

- Added button "Clear features" for heat maps.

- Reorganized buttons in heat map panels.

- Maintainer badge transferred to Federico

Changes in version 0.99.0:

- Initial submission to _Bioconductor_.

[JASPAR2018](https://bioconductor.org/packages/JASPAR2018)
----------

Changes in version 3.6:

NEW FEATURES

- JASPAR2018 data package

[kissDE](https://bioconductor.org/packages/kissDE)
------

Changes in version 0.99.0:

- Pre-release of the package.

[Logolas](https://bioconductor.org/packages/Logolas)
-------

Changes in version 1.3.1:

- Release

- deprecated the two functions - `logomaker` and `nlogomaker` for
  standard and EDLogo.  All logo plots can be now be generated using
  the same function - `logomaker()`. The type argument in this function
  can be chosen to be Logo or EDLogo.

- trimmed the package down from nearly 60 exported functions to just 7
  exported functions.

- The format of the input data is now made more flexible - it allows
  for a vector of character sequences, along with the PFM or the PWM
  matrix as before (see vignette).

- changed the complicated `color_profile` argument into three separate
  arguments - a `color_type` similar to `color_profile$type` argument
  before, a `colors` argument allowing user to choose a cohort of
  colors, and a `color_seed` argument allowing the user to sample
  different colors from the cohort. We now provide a default cohort of
  `colors` as well as default `color_type` in `per-row` (see vignette).
  The user now can do with not worrying about defining `color_profile`
  at all, and use the defaults instead and change the default cohort by
  `color_seed`.  (see vignette).

- added a `return_heights` option in `logomaker()` function that, when
  set to TRUE, returns the information of the heights of the stacks
  used for both standard and EDLogo (see vignette).

- added a `use_dash` argument that, when set to TRUE, would
  automatically detect if the input is a character sequence of PFM
  matrix and perform adaptive scaling of heights (see vignette).

- updated the vignette completely with major focus on the EDLogo
  representation and the use of the current `logomaker()` functionality

- updated the README - with citation information and a demo example
  added.

- Updated the gallery codes
  (https://kkdey.github.io/Logolas-pages/Gallery.html) here to conform
  to the new system of functions.

- Updated the HTML vignette
  (https://kkdey.github.io/Logolas-pages/workflow.html) to match with
  the pdf version of the vignette attached with the package.

- updated README with examples from String logos (histones and mutation
  signatures).

- moved from having data under `inst/extdata` to the `data` folder.

- added a `demo` folder containing some test gallery examples.

Changes in version 1.2.1:

- Added EDLogo plots highlighting both enrichment and depletion

- Added new fill and border styles for the logos

- Added a Dirichlet Adaptive Shrinkage (dash) for adaptively scaling
  position weights

- Added tutorials in the vignette for multi panel Logos plots and
  combining Logolas plots with ggplot2 graphics.

- Some input arguments are deprecated or passed into control parameters

- Background matrix or vector option has been added for comparative
  logo plot visualization given a prior belief.

- PSSM logo plot function added primarily for protein sequence motif
  visualization

- Functions added to compute the heights of the enrichments and
  depletions of the symbols in logo plot.

- Nomenclature added for calling a base at each position.

- Deprecated depletion weight input for unscaled logos + added unscaled
  log and probKL and wKL approaches to the set of possible logos

[maftools](https://bioconductor.org/packages/maftools)
--------

Changes in version 1.6.00:

NEW FUNCTIONS

- clinicalEnrichment - Performs mutational enrichment analysis for a
  given clinical feature.

- signatureEnrichment - Performs sample stratification based on
  signature exposures and enrichment analysis.

- plotEnrichmentResults - Plots results from clinicalEnrichment and
  signatureEnrichment analysis

- lollipopPlot2 - Compare two lollipop plots

SIGNIFICANT USER-LEVEL IMPROVEMENT

- Forstplot now includes summary table within the plot.

- Included capture_size argument in tcgaCompare.

- annovarToMaf can take multiple annovar annotation files and converts
  them to a single MAF cohort.

[MAGeCKFlute](https://bioconductor.org/packages/MAGeCKFlute)
-----------

Changes in version 0.99.19:

- Beutify figures.

- Remove some unnecessary dependencies.

Changes in version 0.99.18:

- Add HeatmapView and BatchRemove.

- Change view distribution functions which show samples separately.

- Change function ReadBeta to be more friendly.

- Label top ten essential genes in SquareView.

- Change some default parameter values to be better.

Changes in version 0.99.10:

- Change all plot function names ended with 'View'.

- Decrease exported functions

- Decrease package denpendcies

- Revise all documents

Changes in version 0.99.1:

- Remove some new errors, such as error trigered by no GroupA genes.

- Allow users to input their own essential genes to do the cell cycle
  normalization

- Allow users to define the number of genes labeled in rank figure,
  default label top 10 and bottom 10 genes

- Add annotation of other organisms

Changes in version 0.99.0:

- FluteMLE and FluteRRA are two main functions in MAGeCKFlute package.
  FluteMLE run pipeline from gene beta scores caculated by MAGeCK MLE,
  while FluteRRA run pipeline based on MAGeCK RRA results.

[matter](https://bioconductor.org/packages/matter)
------

Changes in version 1.5.4:

NEW FEATURES

- Added 'struct' convenience function for on-disk C-style structs
  (wrapper for 'matter_list')

Changes in version 1.5.3:

NEW FEATURES

- Added remainder of Summary group, including 'range', 'min', 'max',
  'prod', 'any', and 'all' methods

- Added options(matter.cast.warning=FALSE) for turning off C type
  coercion warnings

- Exported low-level utilities 'sizeof', 'make_datamode',
  'convert_datamode', and 'widest_datamode'

- Added 'combiner' generic function and method for setting/getting the
  'combiner' for 'sparse_mat'

- Added 'min' and 'max' combiner functions for 'sparse_mat' matrices
  with tolerance > 0

- Added 'biglm' method for 'matter_df' data frames

BUG FIXES

- Fixed bug when coercing 'matter_list' to 'matter_vec'

Changes in version 1.5.2 (2017-11-12):

NEW FEATURES

- All 'matter' subclasses now support endomorphic subsetting via
  'drop=NULL' wherever appropriate

- Setting the 'dim' slot via 'dim<-' now switches the class between
  'matter_vec' and 'matter_arr'

- Added 'virtual_mat' class for virtual matrices

SIGNIFICANT USER-VISIBLE CHANGES

- Use 'drop=NULL' from now on instead of 'drop=NA' to do endomorphic
  subsetting of matter matrices

- Added 'matter_vt' virtual class for matter objects which may exist
  both on-disk and in-memory

- Added 'matter_tbl' virtual class for data tables

BUG FIXES

- Fixed read/write bug when subsetting across atoms

Changes in version 1.5.1:

NEW FEATURES

- Added 'sparse_mat' class for sparse matrics (potentially on-disk)
  with subclasses 'sparse_matc' for CSC matrices and 'sparse_matr' for
  CSR matrices

- Added 'bsearch' function for fast binary searches

- Added 'uuid' function for generating UUIDs as both 'raw' and and
  'character' vectors

- Added a 'checksum' method for doing sha1 and md5 checksums of all
  files associated with a 'matter' object

[mCSEA](https://bioconductor.org/packages/mCSEA)
-----

Changes in version 0.99.0:

- mCSEA release.

[metagenomeFeatures](https://bioconductor.org/packages/metagenomeFeatures)
------------------

Changes in version 1.99.9 (2018-04-02):

- This version is in prepration for the version 2.0

- The MgDb-class definition was been redefined.

- To reduce memory usage the and sequence data is now stored in the
  SQLite file along with the taxonomy data.

- The mgFeatures-class now extends the DataFrame-class instead of the
  AnnotatedDataFrame-class so that mgFeatures can be used to define the
  rowData slot in a summarizedExperiment-class object.

- The vignettes have been revised and new vignettes were added
  providing examples for working with the new class definitions.

- Along with the new class definitions we have annotation packages for
  the three major 16S rRNA databases, SILVA, RDP, and Greengenes.

- Greengenes version 13.8 85% similarity OTUs database is now included
  in the package.

[metagenomeSeq](https://bioconductor.org/packages/metagenomeSeq)
-------------

Changes in version 1.21:

- Numerous changes. Added greater flexibility to fitFeatureModel

[MetaGxOvarian](https://bioconductor.org/packages/MetaGxOvarian)
-------------

Changes in version 0.99.0:

NEW FEATURES

- This is the first iteration of the MetaGxOvarian package. It is
  essentially an updated version of the curatedOvarianData package with
  5 additional datasets

[metaMS](https://bioconductor.org/packages/metaMS)
------

Changes in version 1.15.1:

- Version bump to get inline with the Bioconductor number

[metaseqR](https://bioconductor.org/packages/metaseqR)
--------

Changes in version 1.19.14 (2018-04-12):

NEW FEATURES

- New option for count.type: when "transcript", the DE analysis is
  performed at the transcript level. Only works for Ensembl annotation.

BUG FIXES

- Fixed a leftover which prevented completing the analysis when
  count.type is "exon".

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

Changes in version 1.9.1:

- no changes, check if package passes R CMD build and R CMD check
  without any error messages and vignette can be run without any errors
  &#91;2018-04-16 Mon&#93;

Changes in version 1.7.0:

- implement calculation for parameter m = 0 in NDP &#91;2017-11-14 Tue&#93;

[methimpute](https://bioconductor.org/packages/methimpute)
----------

Changes in version 1.1.1:

NEW FEATURES

- Vignette describes methylation status calling in windows.

- Vignette describes methylation status calling with a separate-context
  model.

- Data import with data.table::fread for faster performance.

[methylKit](https://bioconductor.org/packages/methylKit)
---------

Changes in version 1.5.3:

IMPROVEMENTS AND BUG FIXES

- resolve absolute paths for dbdir argument

- initial check if output tabix file already exists and if yes rename
  output file

Changes in version 1.5.2:

IMPROVEMENTS AND BUG FIXES

- updated links to vignette and presentations in README

- fixed missspelling in show method of methylRawList object

- fixed getSampleID method for methylDiff object

- correct a not working code snippet in the vignette

[methyvim](https://bioconductor.org/packages/methyvim)
--------

Changes in version 1.1.0:

- An updated release of this package for Bioconductor 3.7, released
  April 2018.

- This release primarily implements minor changes, including the use of
  colors in the plots produced by the visualization methods.

[minfi](https://bioconductor.org/packages/minfi)
-----

Changes in version 1.25:

- Added preliminary support for DelayedArray-backed minfi objects. This
  allows disk-backed minfi objects (e.g., using HDF5). This
  functionality is currently recommended only for developers and
  advanced users. A user-friendly interface is currently in
  development. All existing minfi functionality and serialized objects
  should continue to work as it did in versions prior to 1.25. Please
  report any problems to the GitHub issue tracker.

- Fixing bug in functions readGEORawFile() and
  getGenomicRatioSetFromGEO(). These two functions did not work
  (reported an error). They should work now. Thanks to users who
  reported problems at GitHub issues.

- Updated CITATION and citations in the vignette.

[miRNAmeConverter](https://bioconductor.org/packages/miRNAmeConverter)
----------------

Changes in version 1.7.2:

NEW FEATURES

- getMirnasForMirbaseVersion: Retrieve all mature or precursor miRNA
  names from one or several miRBase version(s) and species

SIGNIFICANT USER-LEVEL CHANGES

- Updated vignette

DEPRECATED AND DEFUNCT

- o

[miRsponge](https://bioconductor.org/packages/miRsponge)
---------

Version: 1.5.2
Category: Update ground truth for validation <2018-04-17, Tues
Text:

Version: 1.5.1
Category: Update netModule and spongeValidate function <2018-04-16, Mon
Text:

Version: 1.5.0
Category: Update Reference Manual <2018-04-15, Sun
Text:

Version: 1.3.3
Category: Update netModule function <2018-04-12, Thur
Text:

Version: 1.3.2
Category: Update DESCRIPTION <2018-03-07, Wed
Text:

Version: 1.3.1
Category: Update Reference Manual <2018-03-03, Sat
Text:

Version: 1.3.0
Category: Update Reference Manual <2018-02-28, Wed
Text:

[MLInterfaces](https://bioconductor.org/packages/MLInterfaces)
------------

Changes in version 1.59.2:

- Fix warning stemming from base::seq_len in
  MLInterfaces:::MLIConverter.svm <2018-04-12 Thu>

- Ammend xvalSpec man page <2018-04-12 Thu>

Changes in version 1.59.1:

- Accommodate new points3d syntax in threejs::scatterplot3d for
  plspinHcube <2018-01-20 Sat>

- remove references to rgl, it is not loading on MacOS and threejs
  seems to suffice

[monocle](https://bioconductor.org/packages/monocle)
-------

Changes in version 2.6.1:

- knn-based density peak clustering is not general for all datasets.
  Rolled back to the previous densityPeak clustering algorithm and set
  it to be the default algorithm.  A new Louvain clustering algorithm
  for dealing with large datasets (> 50 k cells) is added.

- A few bug fixes for importCDS, exportCDS functions.

[motifcounter](https://bioconductor.org/packages/motifcounter)
------------

Changes in version 1.3.2:

- Changed contact information

[motifStack](https://bioconductor.org/packages/motifStack)
----------

Changes in version 1.23.12:

- Silence importMatrix for scan.

Changes in version 1.23.11:

- Fix the bug in importMatrix for meme files.

Changes in version 1.23.10:

- Fix the bug in importMatrix for meme files.

Changes in version 1.23.9:

- Fix the warning in windows.

Changes in version 1.23.8:

- Update the documentation plotMotifStackWithRadialPhylog.

Changes in version 1.23.7:

- Update vignette.

Changes in version 1.23.6:

- Add citation.

Changes in version 1.23.5:

- Update authors email.

Changes in version 1.23.4:

- create global environment to save tmp values.

Changes in version 1.23.3:

- change the colorset name from colorBlindness to blindnessSafe

Changes in version 1.23.2:

- fix some bugs in documentation.

Changes in version 1.23.1:

- add a new color set for color blind safe.

[msa](https://bioconductor.org/packages/msa)
---

Changes in version 1.12.0:

- release as part of Bioconductor 3.7

Changes in version 1.11.2:

- fix of code for using custom substitution matrices in ClustalW

Changes in version 1.11.1:

- minor fix in ClustalW

Changes in version 1.11.0:

- new branch for Bioconductor 3.7 devel

[MSnbase](https://bioconductor.org/packages/MSnbase)
-------

Changes in version 2.5.14:

- Fix changed remote location of mzTab example files <2018-04-19 Thu>

- Fix failing centroided unit test (see issue #338) <2018-04-20 Fri>

Changes in version 2.5.13:

- Reduce unit testing time (see #334) <2018-04-13 Fri>

- Fix bug in write.exprs when only one feature data is passed
  <2018-04-17 Tue>

Changes in version 2.5.12:

- Ensure tic and bpi with initial = TRUE calculate the tic and bpi from
  the data (related to issue #332).  <2018-04-10 Tue>

Changes in version 2.5.11:

- Improve combineFeatures manual to document the effect for missing
  values for different types of aggregation methods <2018-04-07 Sat>

- Update robust summary to hangle missing values (see #330).
  <2018-04-09 Mon>

Changes in version 2.5.10:

- New robust summarisation method in `combineFeatures` contributed by
  Ludger Goeminne, Adriaan Sticker and Lieven Clement <2018-04-03 Tue>

- Adapt `utils.removePeaks` to new `IRanges` implementation; thanks to
  H. Pagès for the implementation (see PR #320 for discussion)
  <2018-03-26>.

- Centroiding information is retrieved from raw files (for mzML/mzXML
  files;.  see issue #325 <2018-03-27>

- Add parameter `timeDomain` to `combineSpectra`,
  `combineSpectraMovingWindow` and `estimateMzScattering` allowing to
  perform the grouping of m/z values from consecutive scans based on
  the square root of the m/z values <2018-03-29>.

- Assure feature CV feature variable names are unique when combining
  feature repeatedly (see issue #303) <2018-04-04 Wed>

Changes in version 2.5.9:

- New combineSpectra, combineSpectraMovingWindow, estimateMzScattering
  and estimateMzResulution functions <2018-03-05>.

- New vignette describing profile mode data centroiding <2018-03-05>.

Changes in version 2.5.8:

- New as(MSnExp, data.frame) method <2018-02-16>

- Speed up readMgfData function - see issue #319 <2018-03-13 Tue>

Changes in version 2.5.7:

- MSmap constructor for OnDiskMSnExp objects (see issue #305)
  <2018-01-31 Wed>

- New filterMsLevel,MSnSet method <2018-02-06 Tue>

Changes in version 2.5.6:

- Fixing processing message when combining two MSnSets (reported by
  kamal.fartiyal84 https://support.bioconductor.org/p/100865/#105206)
  <2018-01-24 Wed>

Changes in version 2.5.5:

- Added TMT11-plex <2018-01-17 Wed>

- Add `phenoData<-` method for `pSet` (issue #299 <2018-01-22>.

- Change readMSnSet2 example <2018-01-23 Tue>

Changes in version 2.5.4:

- Add featureData slot to Chromatograms class and add mz, precursorMz
  and productMz methods for Chromatograms (issue #289) <2017-12-18
  Mon>.

- Add readSRMData function to read chromatographic data from SRM (MRM)
  experiments; issue #286 <2018-01-10 Wed>.

- The MSnSetList class has a new `featureData` slot, accessible with
  `fData` to store metadata for the individial MSnSets of the list.
  `MSnSetList` also not has an `sapply` method. <2018-01-10 Wed>

- combineFeatures now has a new `fcol` argument (see issue #195)
  <2018-01-11 Thu>

Changes in version 2.5.3:

- Add filterPrecursorScan for `MSnExp` and `OnDiskMSnExp`; closes issue
  #282 and PR #287 <2017-12-16 Sat>.

- MSnSet to/from SummarizedExperiment coercion (contributed by Arne
  Smits in PR #284) <2017-12-17 Sun>

- Fix inverted M/Z axis in plot3D,MSmap (reported by Sylvain Dechaumet,
  see issue #292) <2017-12-19 Tue>

Changes in version 2.5.2:

- Use automatic backend detection (based on file name and file content)
  that was introduced in mzR version 2.13.1 (issue #275).

- Fix mzML file writing unit tests to work with recently introduced
  header column "filterString" (issue #278).

- Reduce number of comparsions in in internal `fastquant_max` to get
  little speed improvents for isobaric quantification (see PR #280)
  <2017-11-27 Mon>.

- MIAxE, MSnProcess and AnnotatedDataFrame coercion to list methods (in
  related to PR #280) <2017-12-11 Mon>

- Add support to reduce the featureData for OnDiskMSnExp objects (issue
  #285) <2017-12-11 Mon>.

Changes in version 2.5.1:

- Update dependencies (see issue #271)

- Replace HCD by ETD in TMT10ETD's name/description

Changes in version 2.5.0:

- New version for Bioc 3.7 (devel) # MSnbase 2.4

[msPurity](https://bioconductor.org/packages/msPurity)
--------

Changes in version 1.5.1:

- Updates for database creation (can use CAMERA objects now)

- averageSpectra parameter 'MSFileReader' deprecated MSFileReader.
  Should use csvFile instead, MSFileReader option will still work but a
  warning will be given

Changes in version 1.4.1:

- Updates for Galaxy for Spectral Matching

- Spectral matching ra_thres_t bugfix

- Separation of sqlite database creation. Now can be called on it's own
  or with frag4feature (allows the Galaxy tool to be simplified)

[MSstats](https://bioconductor.org/packages/MSstats)
-------

Version: 3.11.6
Date: 2018-04-23
Text: BUG FIXES - SkylinetoMSstatsFormat : fix the inconsistency of
        column name from Skyline output

Version: 3.11.5
Date: 2018-02-22
Text: BUG FIXES - add the package, stringr, for
        DIAUmpiretoMSstatsFormat function

Version: 3.11.4
Date: 2018-02-19
Text: BUG FIXES - add set.seed for sample size calculation of
        classification

Version: 3.11.3
Date: 2018-02-15
Text: BUG FIXES - nonlinear_quantlim : fix the bug for the resampling
        of the blank sample, increase the default number of bootstrap
        samples - designSampleSize : fix the bug NEW FEATURES - new
        function : designSampleSizeClassification,
        designSampleSizeClassificationPlots - Calculate the optimal
        size of training data for classification problem by simulation.
        - new converter functions : DIAUmpiretoMSstatsFormat,
        OpenMStoMSstatsFormat

Version: 3.10.5
Date: 2018-01-10
Text: BUG FIXES - SpectronauttoMSstatsFormat : TRUE or FALSE are
        allowed for the values of the column,
        F.ExcludedFromQuantification. Check the value for this column.

Version: 3.10.4
Date: 2017-12-22
Text: BUG FIXES - MaxQtoMSstatsFormat : 'fewmeasurements' bug fixed

Version: 3.10.2
Date: 2071-11-27
Text: BUG FIXES - make error messages for QQ plot and residual plot, if
        the protein couldn't be fitted by linear mixed effect model.  -
        ProgenesistoMSstatsFormat : make more generalization for
        different format.

[multiClust](https://bioconductor.org/packages/multiClust)
----------

Changes in version 11-14-17:

- -Package version pushed to 1.8.1 -Fixed package vignette issue with
  loading heat map image

[MutationalPatterns](https://bioconductor.org/packages/MutationalPatterns)
------------------

Version: 3.6
Category: MutationalPatterns v1.3.2 (Release date: 2017-10-24
Text:

Version: 3.6
Category: Bugfixes
Text: Removed deprecated functions from previous release.

Version: 3.6
Category: Bugfixes
Text: Improved examples in documentation.

Version: 3.6
Category: MutationalPatterns v1.3.1 (Release date: 2017-10-24
Text:

Version: 3.6
Category: Bugfixes
Text: Fix running of the code examples.

Version: 3.6
Category: MutationalPatterns v1.3.0 (Release date: 2017-10-22
Text:

Version: 3.6
Category: Bugfixes
Text: To determine the transcriptional strand of mutations in genes,
        all mutations that overlap with multiple genes were excluded.
        When these genes are on different strands, it can indeed not be
        determined whether a mutation is on the transcribed or
        untransribed strand. However, if these overlapping genes are
        all on the same strand, the transcriptional strand can be
        determined, but these were unneccesarily removed from the
        analysis. This bug is now fixed, and as a result more mutations
        are now included in the analysis. This bugfix influences the
        results of: 'mut_strand' (previously 'strand_from_vcf') and
        'mut_matrix_stranded'

Version: 3.6
Category: Renamed functions
Text: 'strand_from_vcf' to 'mut_strand'

Version: 3.6
Category: Renamed functions
Text: 'mutation_types' to 'mut_type'

Version: 3.6
Category: Renamed functions
Text: 'mutation_context' to 'mut_context'

Version: 3.6
Category: New features & parameter changes
Text: Replicative strand bias analyses - 'mut_strand' and
        'mut_matrix_stranded' can now be executed in two modes:
        'transcription' (default) or 'replication' - All downstream
        analyses can be performed for both modes with
        'strand_occurrences', 'strand_bias_test' and 'plot_strand_bias'

Version: 3.6
Category: New features & parameter changes
Text: Condensed plotting option for 'plot_96_profile' and
        'plot_192_profile' condensed = F (default), or condensed = T

Version: 3.6
Category: New functions
Text: 'plot_contribution_heatmap': to visualize the relative
        contribution of mutational signatures in a heatmap. Samples can
        be hierarchically clustered.

Version: 3.6
Category: New functions
Text: 'cos_sim': to calculate the cosine similarity between two
        vectors.

Version: 3.6
Category: New functions
Text: 'cos_sim_matrix': to calculate all pairwise similarities between
        mutational profiles

Version: 3.6
Category: New functions
Text: 'cluster_signatures': to hierarchically cluster signatures based
        on cosine similarity

Version: 3.6
Category: New functions
Text: 'plot_cosine_heatmap': to visualize pairwise cosine similarities
        between mutational profiles in a heatmap Sample can be
        hierarchically clustered.

Version: 3.6
Category: MutationalPatterns v1.1.3 (Release date: 2017-04-20
Text:

Version: 3.6
Category: Fourth preparation release for Bioconductor 3.5
Text:

Version: 3.6
Category: Bugfixes
Text: Add missing package to 'Suggest' field.

Version: 3.6
Category: MutationalPatterns v1.1.3 (Release date: 2017-04-20
Text:

Version: 3.6
Category: Third preparation release for Bioconductor 3.5
Text:

Version: 3.6
Category: Bugfixes
Text: Fix running of a unit test.

Version: 3.6
Category: Bugfixes
Text: Fix another build problem for Windows.

Version: 3.6
Category: MutationalPatterns v1.1.2 (Release date: 2017-04-18
Text:

Version: 3.6
Category: Third preparation release for Bioconductor 3.5
Text:

Version: 3.6
Category: Bugfixes
Text: Properly read external data for tests.

Version: 3.6
Category: Bugfixes
Text: Fix build problems on Windows.

Version: 3.6
Category: MutationalPatterns v1.1.1 (Release date: 2017-04-12
Text:

Version: 3.6
Category: Second preparation release for Bioconductor 3.5
Text:

Version: 3.6
Category: MutationalPatterns v1.1.0 (Release date: 2017-04-06
Text:

Version: 3.6
Category: Preparations for Bioconductor release 3.5
Text:

Version: 3.6
Category: Interface changes
Text: 'read_vcfs_as_granges': The 'genome' parameter must now be the
        name of a BSgenome library, to prevent problems with seqlevels
        style. The function now accepts an optional 'group' parameter
        to use a subset of chromosomes. It also accepts the new
        optional 'check_alleles' parameter to significantly speed up
        the reading of VCF files.

Version: 3.6
Category: Interface changes
Text: 'plot_contribution': This function now accepts an optional
        parameter 'palette' to specify custom colors.

Version: 3.6
Category: Performance updates
Text: Implement parallel execution in 'read_vcfs_as_granges',
        'mut_matrix' and 'mut_matrix_stranded'.

Version: 3.6
Category: Bugfixes
Text: Fix 'mut_type_occurences' to handle missing types.

Version: 3.6
Category: Bugfixes
Text: Fix 'mut_matrix' and 'mut_matrix_stranded' to emit warnings when
        processing empty GRanges.

Version: 3.6
Category: Bugfixes
Text: Fix inconsistencies in the README and the vignette.

Version: 3.6
Category: Other changes
Text: Various vignette updates.

Version: 3.6
Category: Other changes
Text: Added unit tests for 'read_vcfs_as_granges', 'mut_matrix', and
        'mut_matrix_stranded'.

Version: 3.6
Category: MutationalPatterns v1.0.0 (Release date: 2016-10-19
Text:

Version: 3.6
Category: Bioconductor release 3.4
Text:

Version: 3.6
Category: MutationalPatterns v0.99.6 (Release date: 2016-10-14
Text:

Version: 3.6
Category: Changes
Text: Renamed functions: 'mut_type_occurences' to
        'mut_type_occurrences', 'strand_occurences' to
        'strand_occurrences'.

Version: 3.6
Category: MutationalPatterns v0.99.5 (Release date: 2016-10-06
Text:

Version: 3.6
Category: Changes
Text: Added deprecation and defunct messages to functions that have
        changed since the v0.99.0.

Version: 3.6
Category: Changes
Text: Various small vignette and reference manual updates.

Version: 3.6
Category: MutationalPatterns v0.99.4 (Release date: 2016-10-05
Text:

Version: 3.6
Category: Changes
Text: Internal package loading changes.

Version: 3.6
Category: Changes
Text: Removed files that do not belong to the package.

Version: 3.6
Category: MutationalPatterns v0.99.3 (Release date: 2016-09-28
Text:

Version: 3.6
Category: Changes
Text: Renamed functions: 'get_mut_context' to 'mutation_context',
        'get_type_context' to 'type_context', 'get_muts' to
        'mutations_from_vcf', 'get_strand' to 'strand_from_vcf'.

Version: 3.6
Category: Changes
Text: Added an explanation for the difference between SomaticSignatures
        and MutationalPatterns in the vignette.

Version: 3.6
Category: MutationalPatterns v0.99.2 (Release date: 2016-09-23
Text:

Version: 3.6
Category: Changes
Text: Renamed functions: 'vcf_to_granges' to 'read_vcfs_as_granges',
        'get_types' to 'mutation_types'.

Version: 3.6
Category: MutationalPatterns v0.99.1 (Release date: 2016-09-13
Text:

Version: 3.6
Category: Changes
Text: Renamed functions: 'read_vcf' to 'vcf_to_granges'.

Version: 3.6
Category: Changes
Text: Removed functions: 'bed_to_granges', 'estimate_rank',
        'rename_chrom'.

Version: 3.6
Category: Changes
Text: Parameter changes: 'plot_rainfall', 'vcf_to_granges'

Version: 3.6
Category: MutationalPatterns v0.99.0 (Release date: 2016-09-12
Text:

Version: 3.6
Category: Changes
Text: Package created

[mzR](https://bioconductor.org/packages/mzR)
---

Changes in version 2.13.8:

- Document missing chrom argument for chromatogram(s)

Changes in version 2.13.7:

- Add a missing header needed on gcc 6.2.0

Changes in version 2.13.6:

- Add MS CV Term IDs for mzR, MSnbase and CAMERA (issue #151)

- Validate exported mzML files using xsd

Changes in version 2.13.5:

- Fix https://github.com/sneumann/xcms/issues/261

- Fix endian.hpp for new c++ versions (see PR #149)

Changes in version 2.13.4:

- Fix error (see issue #145)

Changes in version 2.13.3:

- Link against Rhdf5lib, allows to read mz5 also on Windows

- Use Rhdf5lib 1.1.4 with c++ headers in /include

- fix BiocStyle related issue in Vignette on Windows

Changes in version 2.13.2:

- Add chromatogramHeader method to read header information for
  chromatograms from an mzML file.

Changes in version 2.13.1:

- Read filter string from mzML files and add it to the data.frame
  returned by the header function (see MSnbase issue #278).

- openMsFile automatically determine the backend to use based on file
  extension and content.

[ndexr](https://bioconductor.org/packages/ndexr)
-----

Changes in version 1.1.7:

- FIX: error in RCX => RCXgraph => RCX conversion

- FIX: NDEx server update for return columns in network list, summary
  and metadata; metadata also not nested anymore

- FIX: tests crashed because of missing network; updated used UUID to
  new version of the public one from ndextutorials

- exclude ..Rcheck from git; added "ndexr" to user agent header

- FIX: build error caused by 'metadata:properties' now being optional

- UPDATE: minor bugfixes due to ndex server update. Added api for ndex
  server version 2.2

Changes in version 1.1.2:

- **Breaking changes of the class and function names!** **NGraph** was
  renamed to **RCXgraph** to avoid naming Disambiguities! **Deprecated
  Functions:**

- *rcx_toRCXgraph*

- *rcxgraph_fromRCX*

- *rcxgraph_toRCX*

- *rcx_fromRCXgraph* Therefore the new funtions are called:

- *rcx_toRCXgraph*

- *rcxgraph_fromRCX*

- *rcxgraph_toRCX*

- *rcx_fromRCXgraph*

[NOISeq](https://bioconductor.org/packages/NOISeq)
------

Changes in version 2.22.1 (2018-02-01):

- Fixed some bugs.

[OmaDB](https://bioconductor.org/packages/OmaDB)
-----

Version: 0.99.9
Date: 2017-10-14
Category: Initial release to Bioconductor
Text:

[omicplotR](https://bioconductor.org/packages/omicplotR)
---------

Changes in version 0.99.4:

- changed vignette to html instead of pdf.

- edited vignette.

Changes in version 0.99.3:

- bug fix. extrainformation.rmd changed to extrainformation.Rmd.

Changes in version 0.99.0:

- release version.

[OncoSimulR](https://bioconductor.org/packages/OncoSimulR)
----------

Changes in version 2.10.0:

- probDetect mechanism changed. This could be a BREAKING CHANGE.  The
  expression divides by the baseline. For fixed initSize, this is
  simply a matter of changing the cPDetect.

- fixation allows exact genotypes, includes tolerance, and checks for a
  successive number of specified periods

- LOD: using only the strict Szendro et al. meaning.

- POM: computed in C++.

- Using fitness landscape directly when given as input (no conversion
  to epistasis) and several improvements in speed when using fitness
  landscapes as input.

Changes in version 2.9.10 (2018-04-19):

- test.Z-fixation: some tests only on Linux because rng is done in C++.

Changes in version 2.9.9 (2018-04-10):

- probDetect mechanism changed. This could be a BREAKING CHANGE.  The
  expression divides by the baseline. For fixed initSize, this is
  simply a matter of changing the cPDetect.

Changes in version 2.9.8 (2018-03-26):

- fixation allows exact genotypes, includes tolerance, and checks for a
  successive number of specified periods

Changes in version 2.9.7 (2018-02-20):

- fixed crash in some conditions when run with stringsAsFactors = FALSE
  as global option

Changes in version 2.9.6 (2017-12-27):

- Updated citation.

- An example (in miscell-files) about using and stopping with modules.

- Prototype for sampling the single larges pop at last period (function
  largest_last_pop, commented out for now).

Changes in version 2.9.5:

- samplePop: new option "single-nowt"

Changes in version 2.9.4 (2017-11-30):

- Deal with the very rare NULL simulations in summary.

Changes in version 2.9.3 (2017-11-27):

- Make clang happy (do not use flandscape as DataFrame)

Changes in version 2.9.2 (2017-11-24):

- LOD: using only the strict Szendro et al. meaning.

- POM: computed in C++.

Changes in version 2.9.1 (2017-11-10):

- Using fitness landscape directly when given as input (no conversion
  to epistasis)

[ORFik](https://bioconductor.org/packages/ORFik)
-----

Changes in version 1.0.0:

- first release of ORFik - find Open Reading Frames, automatic RiboSeq
  footprint shifts, reassignment of Transcription Start Sites with the
  use of CageSeq, plethora of gene identity functions from scientific
  publications

[PanVizGenerator](https://bioconductor.org/packages/PanVizGenerator)
---------------

Changes in version 1.7.0:

- Fix a bug in creation of GO graph where invalid edges were included

[PathoStat](https://bioconductor.org/packages/PathoStat)
---------

Changes in version 1.8.1:

- Update data uploading

- Modify visualization

- Add biomarker tab

[pathview](https://bioconductor.org/packages/pathview)
--------

Changes in version 1.18.2:

- korg now include 5245 KEGG species. To speed up pathview package
  loading, an updated version of korg is now checked out from Pathview
  Web server only when it is acutally used and needed.

[pcaExplorer](https://bioconductor.org/packages/pcaExplorer)
-----------

Changes in version 2.6.0:

NEW FEATURES

- Automatically computing size factors where required

- Added progress indication when compiling the report

BUG FIXES

- Fixed after changes in threejs package

- Edited dropdown menu to remove unused green badge

- Menus start expanded on the side, again

- theme_bw applied when needed, corrected previous behavior

OTHER NOTES

- Updated citation infos

- Slight difference in handling validate/need errors

[PGA](https://bioconductor.org/packages/PGA)
---

Changes in version 1.9.1:

- As the "getTable()" of rtracklayer would produce a Bad Request error.
  We temporarily disabled the code checking of the
  "PrepareAnnotationRefseq2" function until the bug had been fixed.
  Please use the "PrepareAnnotationEnsembl2" to prepare the annotation
  file instead of "PrepareAnnotationRefseq2".

[phantasus](https://bioconductor.org/packages/phantasus)
---------

Changes in version 0.99.34:

- Detecting conditions and replicates in GEO data

Changes in version 0.99.28:

- Submit to Enrichr tool

- Changing shapes in PCA plot

- Better support for OS X and Windows

Changes in version 0.99.24:

- Shiny GAM tool added

- PCA plot: pretty labels, auto-redraw

- Bug fixes

Changes in version 0.99.23:

- Added options "Maximum Mean Probe" and "Maximum Median Probe" to
  Collapse Tool

- Fixed a bug emerging when one applies R-based function after
  collapsing by columns

Changes in version 0.99.22:

- Fixed AdjustTool

- Removed obsolete output produced by servePhantasus

- Updated tutorial

Changes in version 0.99.21:

- Merged changes from morpheus.js repository

- Added function `reparseCachedESs` for updating downloaded
  GEO-datasets, which are saved in cache

Changes in version 0.99.20:

- A bit better code coverage with tests

Changes in version 0.99.19:

- Safer procedure for phenoData parsing

- Even sample-empty datasets deserve to have displayed meta

Changes in version 0.99.18:

- New phenoData parsing strategy for GEO Series

Changes in version 0.99.17:

- Fixed issue with null featureData for GSE37270

Changes in version 0.99.16:

- Fixed issue with incorrect perfomance of LimmaTool when columns are
  passed as args

Changes in version 0.99.15:

- Updated tutorial

Changes in version 0.99.14:

- Correct error messages when dataset is not present in GEO

Changes in version 0.99.13:

- Added possibility to load datasets from preloaded rda-files in
  specialized directory on server

- servePhantasus now opens browser with web-application automatically
  on start-up

- Added tutorial on main site

- Updated biocViews

- Updated tutorial and mans according to changes

- Added experimental support for loading RNA-seq dataset from GEO

Changes in version 0.99.12:

- Added possibility to load GEO-dataset by specifying its id in link
  parameters

[piano](https://bioconductor.org/packages/piano)
-----

Version: 1.20.0
Category: none yet
Text:

Version: 1.18.1
Category: BUG FIXES
Text: Fix parsing of gmt files when gene-set name contains spaces

[Pigengene](https://bioconductor.org/packages/Pigengene)
---------

Changes in version 1.5.22 (2018-04-28):

Bug Fixes

- stats::cor is used in compute.pigengene and draw.cor functions. See
  the NAMESPACE for the important reason.

Changes in version 1.5.9 (2018-03-12):

General

- cor is imported from WGCNA, but not from stats, because WGCNA does
  not call the cor function properly.

Changes in existing functions

- RsquaredCut is added to the arguments of the one.step.pigengene
  function.

Changes in version 1.5.6 (2018-01-19):

Changes in existing functions

- In the get.fitted.leaf function, the function C50:::as.party.C5.0 is
  used, which used to be exported in the previous versions of C50, but
  not in version 0.1.1.

Changes in version 1.5.2 (2017-11-10):

New functions

- The check.nas function is now exported.

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

[PowerExplorer](https://bioconductor.org/packages/PowerExplorer)
-------------

Changes in version 0.99.0:

- The review version.

[powerTCR](https://bioconductor.org/packages/powerTCR)
--------

Version: 0.1.1
Category: Getting compliant with BioConductor policies
Text:

Version: 0.1.1
Category: INITIAL RELEASE
Text:

[pRoloc](https://bioconductor.org/packages/pRoloc)
------

Changes in version 1.19.4:

- Fix regression bug in knntl function <2018-04-12 Thu>

Changes in version 1.19.3:

- Use `dplyr::left_join` without attaching `dplyr` to avoid collision
  between `Biobase::exprs` and `dplyr::exprs` <2018-04-04 Wed>.

- Typo in warning to install rgl <2018-03-27 Tue>

Changes in version 1.19.2:

- Fix bug in QSep that prevented to set non-default fcol <2018-01-29
  Mon>

Changes in version 1.19.1:

- Fix bug in private dimred and set appropriate number of colnames
  <2017-11-07 Tue>

- New `nipals` method in dimensionality reduction for plot2D (closes
  issue #103) <2018-01-16 Tue>

Changes in version 1.19.0:

- Bioconductor devel 3.7

[Prostar](https://bioconductor.org/packages/Prostar)
-------

Changes in version 1.11.13:

BUG FIXES

- Normalization: "Sum by columns" has been modified to provide
  log-abundances compatibles with the other treatments. It can be done
  "for each condition independantly" or "globaly".

NEW FEATURES

- Descriptive statistics: The expression datasets are colored w.r.t the
  nature of missing value (POV or MEC) even when the value has been
  imputed

- Filtering: Manage designs with more than 2 conditions and with
  conditions containing different number of samples

- Filtering: UI more user friendly for the string-based filtering (Tab
  2)

- Normalization: A few modifications in the UI and

- Imputation (protein level): Distinction between missing values on an
  entire condition (Missing on the Entire Condition) and the other ones
  (Partially Observed Value)

- Imputation (protein level): for the POV, it is possible to use SLSA
  which take into account the experimentaldesign experimental

- Imputation (protein level): imputations are all processed condition
  by condition

- Differential analysis: All tests can process datasets with conditions
  of different number of samples

- Differential analysis: Limma takes into account all the hierarchical
  experimental designs

- GO analysis: Add the GeneID nomenclature.

Changes in version 1.11.4:

BUG FIXES

- A bug in the string-based filtering tool was fixed. The case where an
  entity could be both contaminants and reverse was not takien into
  account.This lead to wrong number in the plot.

- Correction of the beahviour of the table in the experimental design
  (convert Data tool). When the user copy-paste some lines it may add
  unneeded rows. These rows can be deleted with an option in the
  contextual menu.

[PureCN](https://bioconductor.org/packages/PureCN)
------

Changes in version 1.10.0:

- New normal database format

- Runtime performance improvements (skip unlikely local optima, support
  for BiocParallel in runAbsoluteCN, pre-calculation of mapping bias)

- Support for replication timing scores in coverage normalization

- More accurate confidence intervals in callMutationBurden

- More accurate copy numbers for high-level amplifications

- Very low or high coverage samples are now by default dropped in
  normal database creation (less than 25% or more than 4 times the
  median sample coverage)

- Improved support for third-party upstream tools like GATK4
  (experimental)

- More checks for wrong or sub-optimal input and providing suggestions
  for fixing those issues

- Gibbs sampling of log tumor/normal coverage error rate

- Better imputation of mapping bias (instead of smoothing over
  neighboring variants in the sample, smooth over neighboring SNPs in
  the pool of normals - only available when pre-calculated)

- Experimental support for indels

- Code cleanups (switch to testthat, removed several obsolete and minor
  features) API CHANGES

- renamed gc.gene.file to interval.file since it now provides more than
  GC-content and gene symbols

- plotAbs ids changed to id (this function now only plots a single
  purity/ploidy solution)

- changed default of runAbsoluteCN max.logr.sdev to 0.6 (from 0.75)

- createTargetWeights does not require tumor coverages anymore

- calculateGCContentByInterval was renamed to preprocessIntervals

- renamed plot.gc.bias to plot.bias in correctCoverageBias since it now
  also includes replication timing

- added calculateMappingBiasVcf to pre-compute mapping bias from a
  panel of normal VCF, thus avoiding time loading and parsing of huge
  VCFs

- max.homozygous.loss now defines the maximum fraction of a chromosome
  lost, not the whole genome, to avoid wrong maximum likelihood
  solutions with completely deleted chromosome arms

[qcmetrics](https://bioconductor.org/packages/qcmetrics)
---------

Changes in version 1.17.1:

- Remove deprecated html_document2 <2018-01-16 Tue>

[qpgraph](https://bioconductor.org/packages/qpgraph)
-------

Changes in version 2.14:

USER VISIBLE CHANGES

- Added argument marginalization argument to qpPathWeight().

BUG FIXES

- Bugfix in the 'rUGgmm()' function when called with an undirected
  graph defined by numeric (integer) vertices, to respect the numeric
  ordering of those vertex lables, i.e., avoid simulating a graph with
  ordered vertices "1", "10", "11", etc., and get instead "1", "2",
  "3", etc.

[RamiGO](https://bioconductor.org/packages/RamiGO)
------

Changes in version 1.26.0:

- Removed RCytoscape code and dependencies since RCytoscape is
  deprecated.

[ramwas](https://bioconductor.org/packages/ramwas)
------

Changes in version 1.3.0:

NEW FEATURES

- Improved error reporting from parallel jobs

- Improved logging behavior

- Added multiple input checks in the pipeline functions

- Improved MLK behavior in Microsoft R

- Added BED and BedGraph export

- Added Manhattan plot function

- Added ROC curve with AUC calculation

- Made fragment size estimation functions public

- Made PCA plotting functions public

- Added rwDataClass for convenient data access.

- Added getLocations, getMWAS*, getData* functions

- Improved QQ-plot function (more parameters)

- Speed up pvalue2qvalue function for sorted vectors

- Added multithreading in ramwas0createArtificialData()

- Improved testPhenotype() to have consistent input (variables by
  columns)

- BAM scanning &#91;ramwas1scanBams&#93; rescans BAMs newer than the rbam files

BUG FIXES

- Use importMethodsFrom(filematrix, as.matrix) in NAMESPACE

[rCGH](https://bioconductor.org/packages/rCGH)
----

Changes in version 1.9.1:

Minor changes

- Vignette formating in order to meet the guidelines of the
  "Bioconductor LaTeX Style 2.0".

[RcisTarget](https://bioconductor.org/packages/RcisTarget)
----------

Changes in version 0.99.7:

- Main update: Re-formatted ranking databases.  They are now loaded
  from .feather files (and therefore, they are transposed: motifs are
  stored as rows, and genes/regions as columns, which allows to load
  only specific genes/regions)

Changes in version 0.99.6:

- Main update: Re-formatted annotation database. For easier subsetting.

Changes in version 0.99.1:

- Main update: New class for the Rankings (rankingRcisTarget)

Changes in version 0.7:

- Converted the ranking databases to an S4 class.

Changes in version 0.6:

- Added function addLogo()

- AUC is now returned as a class

[RCy3](https://bioconductor.org/packages/RCy3)
----

Version: 2.0
Category: For Developers
Text: Reorganized functions into files corresponding to CyREST API,
        e.g., Collections, CytoscapeSystem, Layouts, Networks, etc.
        Normalized all documentation using roxygen2 Streamlined
        interfaces to CyREST and Commands API (see above), greatly
        facilitating the implementation of any new functions matching
        CyREST or Command API additions Reverted all single-instance
        methods to simple functions, replacing class-based signatures
        with simple default values Established handy functions for
        validating network and view SUIDs getNetworkName getNetworkSuid
        getNetworkViewSuid

Version: 2.0
Category: Deprecated
Text: Outdated function names

Version: 2.0
Category: Defunct
Text: CytoscapeConnection and CytoscapeWindow classes, functions and
        parameters

[ReactomePA](https://bioconductor.org/packages/ReactomePA)
----------

Changes in version 1.23.2:

- re-implement viewPathway <2018-03-15, Thu>

- mv site to https://guangchuangyu.github.io/software/ReactomePA

Changes in version 1.23.1:

- import enrichplot

[recount](https://bioconductor.org/packages/recount)
-------

Changes in version 1.5.11:

BUG FIXES

- Change some examples to dontrun and improve the code that cleans up
  after the tests. This should reduce the size of files left in tmp
  although they didn't seem too big to begin with.

Changes in version 1.5.9:

SIGNIFICANT USER-VISIBLE CHANGES

- The functions add_metadata() and add_predictions() now return the
  sample metadata or predictions when the 'rse' argument is missing.

Changes in version 1.5.6:

NEW FEATURES

- Added the function add_metadata() which can be used to append curated
  metadata to a recount rse object. Currently, add_metadata() only
  supports the recount_brain_v1 data available at
  http://lieberinstitute.github.io/recount-brain/ and to be further
  described in Razmara et al, in prep, 2018.

Changes in version 1.5.5:

BUG FIXES

- Fix doc link in geo_characteristics() which affected the Windows
  build machines.

Changes in version 1.5.4:

BUG FIXES

- Fix a unit test for download_study(), add another test for the
  versions, and fix a NOTE in R CMD check.

Changes in version 1.5.3:

NEW FEATURES

- download_study() can now download the transcript counts
  (rse_tx.RData) files. The transcript estimation is described in Fu et
  al, 2018.

SIGNIFICANT USER-VISIBLE CHANGES

- download_study() now has a version parameter (defaults to 2). This
  argument controls which version of the files to download based on the
  change on how exons were defined. Version 1 are reduced exons while
  version 2 are disjoint exons as described in further detail in the
  documentation tab of the recount website
  https://jhubiostatistics.shinyapps.io/recount/.

- recount_url and the example rse_gene_SRP009615 have been updated to
  match the changes in version 2.

[recoup](https://bioconductor.org/packages/recoup)
------

Changes in version 1.7.2 (2018-04-25):

NEW FEATURES

- New faster annotation building system which keeps also versions. May
  break older annotation stores. A rebuild is advised.

- More supported genomes

BUG FIXES

- Fixed a bug resulting in profile bleeding when using a mask of ranges
  from a BAM file.

[RedeR](https://bioconductor.org/packages/RedeR)
-----

Changes in version 1.28.0:

- Regular maintenance, stable release.

[rGREAT](https://bioconductor.org/packages/rGREAT)
------

Changes in version 1.11.1:

- save regions in gzip file then submit to GREAT

- will not change the content of the input regions

[rgsepd](https://bioconductor.org/packages/rgsepd)
------

Changes in version 1.11.2:

BUG FIXES

- updating dependencies caused a crash. Investigating. Biomart seems to
  no longer support "hsapiens_gene_ensembl" (sometimes/bug)
  Intermittent bug fixed in a later biomaRt edition.

- R 3.5.0 support calls for a tweak in GOCatEngine.R

[rhdf5](https://bioconductor.org/packages/rhdf5)
-----

Changes in version 2.24.0:

NEW FEATURES

- Removed bundled HDF5 library - rhdf5 now depends on Rhdf5lib.  This
  updates the version of HDF5 to 1.8.19.

- Functions H5Ldelete() and h5delete() added to provide mechanisms for
  removing items from HDF files.

- Added argument `native` to many functions, which allows data to be
  treated as row-major rather than column-major, improving portability
  with other programming languages.

- Added function H5Sunlimited() allowing creation of extensible
  datasets - thanks to Brad Friedman

BUG FIXES

- Datasets can now be subset using `&#91;` and a range of values e.g.
  did&#91;,1:5&#93;.

- Writing a data.frame that contains factors and setting
  DataFrameAsCompound=FALSE now works.

- Many functions that would leave open file handles after exiting under
  error conditions have been fixed.

- Performance improvements in h5read().

[Rhtslib](https://bioconductor.org/packages/Rhtslib)
-------

Changes in version 1.12.0:

SIGNIFICANT USER-VISIBLE CHANGES

- The HTSlib C library was updated from version 1.1 to version 1.7.

[RiboProfiling](https://bioconductor.org/packages/RiboProfiling)
-------------

Changes in version 1.7.2:

- Corrected Biocstyle bug in author field of the vignette

[RnBeads](https://bioconductor.org/packages/RnBeads)
-------

Changes in version 1.11.9:

- Bugfixes and improvements in gender prediction, GEO import,
  differential variability, installation and others

Changes in version 1.11.8:

- Fixes for bioconductor warnings (combine, etc.)

- docoupled missing value imputation from differential methylation

Changes in version 1.11.7:

- enhanced cross-platform combination methods

- improved installation routines

- enhanced plots in exploratory analysis module

- better LOLA annotation

- improved performance of missing value imputation

- documentation updates

- several minor bugfixes

Changes in version 1.11.6:

- Improved combining methods for RnBSet objects of different data types

- Interpretation of sample mean methylation levels and other statistics

- Improved plots in exploratory analysis module

- Improved missing value imputation

- Several minor bugfixes

Changes in version 1.11.5:

- bugfixes (loading, imputation)

Changes in version 1.11.4:

- Introducing RnBeadsDJ, a shiny-based interface for running RnBeads
  analyses and modules (run with rnb.run.dj())

- Added implementation of the LUMP algorithm for immune cell content
  estimation

- The default normalization method was changed to "wm.dasen"

- Background normalization is now disabled per default.

- Option backwards compatibility

- bugfixes (LOLA dependencies, 1-sample differential variability, QC
  visualization, ...)

Changes in version 1.11.2:

- Fixed a bug concerning failed trackhub exports on Windows

[roar](https://bioconductor.org/packages/roar)
----

Changes in version 1.15.2:

- Updated the vignette with a new url for gtf files

- ranges() for Hits -> overlapsRanges()

[rols](https://bioconductor.org/packages/rols)
----

Changes in version 2.7.2:

- Nothing yet

Changes in version 2.7.1:

- Fix failing partOf unit test <2017-11-27 Mon>

Changes in version 2.7.0:

- Bioconductor devel 3.7

[RPA](https://bioconductor.org/packages/RPA)
---

Changes in version 1.35.2 (2017-11-12):

- phyloseq moved from dependency to import

- vignette format updated

- explicitly imported all functions

[Rsamtools](https://bioconductor.org/packages/Rsamtools)
---------

Changes in version 1.31:

BUG FIXES

- (v.1.31.3) pileup() examples require min_base_quality = 10. See
  https://support.bioconductor.org/p/105515/#105553

[RSeqAn](https://bioconductor.org/packages/RSeqAn)
------

Changes in version 0.99.0:

- First release, only headers files have been included.

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

Version: 1.30.0
Category: NEW FEATURES
Text:

Version: 1.30.0
Category: o New parameters in featureCounts() to give more control on
        the size of overlap between read and feature - 'nonOverlap' and
        'nonOverlapFeature
Text:

Version: 1.30.0
Category: o New parameter in featureCounts() to specify the path where
        files containing counting results for each read are saved -
        'reportReadsPath
Text:

Version: 1.30.0
Category: o New parameter in align() and subjunc() to allow reads in
        the mapping output to have the same order as in the FASTQ file
        - 'keepReadOrder
Text:

Version: 1.30.0
Category: o Reduce running time in outputting BAM files in align() and
        subjunc() via saving data to disk with multi threads
Text:

Version: 1.30.0
Category: o Annotation file (eg. GTF, SAF or VCF annotation) in
        featureCounts(), align(), subjunc() and exactSNP() can be
        provided as a gzipped file
Text:

Version: 1.30.0
Category: o Tilda ('~') is allowed to be included in file names
        provided to functions.
Text:

Version: 1.30.0
Category: o Improve running time of featureCounts for the processing of
        BAM files generated by some other aligners and the Picard tool
Text:

[RTN](https://bioconductor.org/packages/RTN)
---

Changes in version 2.4.0:

- Improved RTN workflows.

[RTNduals](https://bioconductor.org/packages/RTNduals)
--------

Changes in version 1.4.0:

- Improved workflow integration with derivative packages (RTN /
  RTNsurvival).

[rWikiPathways](https://bioconductor.org/packages/rWikiPathways)
-------------

Changes in version 1.0.0:

NEW FEATURES

- Added derivative functions for both listPathways and
  getPathwaysByXref to allow specifying the return of simple lists of
  WPIDs, URLs or names.

- Added support for PNG from getColoredPathway

- Added functions for all remaining web service API methods

- Added a download function for archived pathway sets in various
  formats

- Added Overview and BridgeDbR vignettes

SIGNIFICANT CHANGES UNDER THE HOOD

- Updated all functions to use REST calls instead of Curl

- Simplified function development by abstracting as dedicated
  wikipathwaysGET function

- Added tests for all new functions

- Added enumerated parameters in a number of places for better input
  validation

BUG FIXES

- First official release. No bugs... yet!

[S4Vectors](https://bioconductor.org/packages/S4Vectors)
---------

Changes in version 0.18.0:

NEW FEATURES

- The package gets a new vignette: S4VectorsOverview.Rnw The material
  in this new vignette comes from the IRangesOverview.Rnw vignette
  located in the IRanges package. All the S4Vectors-specific material
  was moved from the IRangesOverview.Rnw vignette to the new
  S4VectorsOverview.Rnw vignette.

- All Vector derivatives now support 'x&#91;i, j&#93;' by default. This allows
  the user to conveniently subset the metadata columns thru 'j'. Note
  that GenomicRanges objects have been supporting this feature for
  years but now all Vector derivatives support it. Developers of Vector
  derivatives with a true 2-D semantic (e.g. SummarizedExperiment) need
  to overwrite this.

- rank() now suports 'by' on Vector derivatives.

- Add concatenateObjects() generic and methods for LLint, vector,
  Vector, Hits, and Rle objects. This is a low-level generic intended
  to facilitate implementation of c() on vector-like objects. The
  "concatenateObjects" method for Vector objects concatenates the
  objects by concatenating all their parallel slots. The method behaves
  like an endomorphism with respect to its first argument 'x'. Note
  that this method will work out-of-the-box and do the right thing on
  most Vector subclasses as long as parallelSlotNames() reports the
  names of all the parallel slots on objects of the subclass (some
  Vector subclasses might require a "parallelSlotNames" method for this
  to happen). For those Vector subclasses on which concatenateObjects()
  does not work out-of-the-box or does not do the right thing, it is
  strongly advised to override the method for Vector objects rather
  than trying to override the (new) "c" method for Vector objects with
  a specialized method. The specialized "concatenateObjects" method
  will typically delegate to the method below via the use of
  callNextMethod(). See "concatenateObjects" methods for Hits and Rle
  objects for some examples. No Vector subclass should need to override
  the "c" method for Vector objects.

- Major refactoring of [[<- for List objects. It's now based on a new
  "setListElement" method for List objects that relies on `[<-` for
  replacement, c() for appending, and `[` for removal, which are the 3
  operations that setListElement() can perform (depending on how it's
  called). As a consequence [[<- now works out-of-the box on any List
  derivative for which `[<-`, c(), and `[` work.

SIGNIFICANT USER-VISIBLE CHANGES

- endoapply() and mendoapply() are now regular functions instead of
  generic functions.

- A couple of minor improvements to how default "showAsCell" method
  handles list-like and non-list like objects.

- Replace strsplitAsListOfIntegerVectors() with
  toListOfIntegerVectors(). (The former is still available but
  deprecated in favor of the latter.) The input of
  toListOfIntegerVectors() now can be a list of raw vectors (in
  addition to be a character vector), in which case it's treated like
  if it was 'sapply(x, rawToChar)'.

- A couple of optimizations to "[<-" method for DataFrame objects (see
  commit e63f4cfd637e3471e4b04015c2938348df17e14a).

DEPRECATED AND DEFUNCT

- phead() and ptail() are deprecated in favor of IRanges::heads() and
  IRanges::tails().

- strsplitAsListOfIntegerVectors() is deprecated in favor of
  toListOfIntegerVectors().

BUG FIXES

- The mcols() setter no more tries to downgrade to DataFrame a supplied
  right value that extends DataFrame (e.g. DelayedDataFrame).

- 'DataFrame(I(x)) and as(I(x), "DataFrame")' now drops the I()
  wrapping before storing 'x' in the returned object. This wrapping was
  ugly, not needed, and breaking S4 objects.

- Fix a couple of long-standing bugs in DataFrame subassignment: - Bug
  in the "[<-" method for DataFrame objects where replacing the 1st
  variable with a rectangular object (e.g. x[1] <-
  DataFrame(aa=I(matrix(1:6, ncol=2)))) was returning a DataFrame with
  the "nrows" slot set incorrectly. - A couple of bugs in the
  "replaceROWS" method for DataFrame objects when used in "rbind mode"
  i.e. when max(i) > nrow(x).

- Fix bug in "cbind" method for DataFrame where it was appending X to
  the column names in some situations (see
  https://github.com/Bioconductor/S4Vectors/issues/8).

- Fix order() on SortedByQueryHits objects (see
  https://github.com/Bioconductor/S4Vectors/issues/6).

- Fix bug in internal new_Hits() constructor where it was not returning
  an object of the class specified via 'Class' in some situations.

- "lapply" for SimpleList objects now calls match.fun(FUN) internally
  to find the function to apply.

[sampleClassifier](https://bioconductor.org/packages/sampleClassifier)
----------------

Version: 22.11.2017
Text:

[scater](https://bioconductor.org/packages/scater)
------

Changes in version 1.7.18:

- Refactored calculateQCMetrics() to ignore potential non-linearity,
  rank genes by highest expression, rename automatically generated
  union sets, allow for output of a compact format.

- Refactored all plotting functions to allow access to nested fields in
  the colData() or rowData(), by supplying a character vector.

- Refactored plotTSNE(), plotPCA(), etc. to dispatch to the calculation
  functions (e.g., runTSNE(), runPCA()), with argument checks.

- Refactored plotColData() and plotRowData() to use the same argument
  types as other functions rather than aes= input.

- Removed all plotting functions that do not operate on
  SingleCellExperiment objects.

- Deprecated read10xResults(), downsampleCounts() in favour of methods
  from the DropletUtils package.

- Deprecated scater_gui() in favour of methods from the iSEE package.

- Deprecated normalizeExprs() as this function made very little sense.

- Added plotHeatmap() function, for easy plotting of heatmaps.

- Added librarySizeFactors() function, to compute size factors from
  library sizes.

- Added by_exprs_values= argument to many plotting functions, to
  distinguish direct plotting of expression values from their use in
  aesthetics.

- Renamed arguments in plotHighestExprs(), plotExprsVsTxLength(),
  plotExprsFreqVsMean() for greater clarity.

- Added centreSizeFactors() function for centralized size factor
  centering.

- Added size_factor_grouping= argument to normalizeSCE(), calcAverage()
  and calculateCPM().

- Added subset_row= argument to calculateCPM().

- Consolidated size_factors= argument into use_size_factors= for
  calcAverage(), calculateCPM().

- Modified normalizeSCE() so that centre_size_factors=FALSE does not
  use centred size factors at all during normalization.

[scDD](https://bioconductor.org/packages/scDD)
----

Changes in version 1.3.4 (2018-03-26):

- The test for differences in proportion of zeroes now uses the Wald
  test p-value instead of the likelihood ratio. This will give very
  similar results, but the Wald test is slightly more conservative.

Changes in version 1.3.3 (2018-03-23):

- An option has been added to skip the categorization step if only
  intereseted significance of difference. This will speed up
  computation.

- The testing zeroes step and the KS test have been parallelized to
  speed up computation.

[scmeth](https://bioconductor.org/packages/scmeth)
------

Changes in version 0.99.38:

PKG FEATURES

- scmeth is a package to analyze methylation data and it generates a
  comprehensive quality control report

- Most functions take bsseq object as the input

[scPipe](https://bioconductor.org/packages/scPipe)
------

Changes in version 1.0.9:

- fix the bug that might misalign the exon mapping reads

- support gff files from gencode or refseq

Changes in version 1.0.8:

- update the gene id annotation code

Changes in version 1.0.6 (2017-12-18):

- Now the id conversion can also be done by using the bioconductor
  annotation package, when biomart fails to connect.

Changes in version 1.0.5 (2017-12-14):

- Fixed bugs in slim report and trimbarcode error message

- Fix incomplete error message

- Documentation updates, new functions and bug fixes

Changes in version 1.0.4 (2017-12-04):

- In `detect_outlier`, give more informative error message when some
  cells or QC metrics have zero values.

Changes in version 1.0.3 (2017-12-03):

- fix a bug in `validObject`. the default value for gene id and
  organism is set to NA

Changes in version 1.0.2 (2017-12-01):

- fix errors in unittest

Changes in version 1.0.1 (2017-11-28):

- Bug Fix: Fixed handling of colData through QC_metrics
  (https://github.com/LuyiTian/scPipe/issues/34)

[scran](https://bioconductor.org/packages/scran)
-----

Changes in version 1.7.28:

- Modified decomposeVar() to return statistics (but not p-values) for
  spike-ins when get.spikes=NA. Added block= argument for mean/variance
  calculations within each level of a blocking factor, followed by
  reporting of weighted averages (using Fisher's method for p-values).
  Automatically record global statistics in the metadata of the output
  for use in combineVar().  Switched output to a DataFrame object for
  consistency with other functions.

- Fixed testVar() to report a p-value of 1 when both the observed and
  null variances are zero.

- Allowed passing of arguments to irlba() in denoisePCA() to assist
  convergence. Reported low-rank approximations for all genes,
  regardless of whether they were used in the SVD. Deprecated design=
  argument in favour of manual external correction of confounding
  effects. Supported use of a vector or DataFrame in technical= instead
  of a function.

- Allowed passing of arguments to prcomp_irlba() in buildSNNGraph() to
  assist convergence. Allowed passing of arguments to get.knn(),
  switched default algorithm back to a kd-tree.

- Added the buildKNNGraph() function to construct a simple
  k-nearest-neighbours graph.

- Fixed a number of bugs in mnnCorrect(), migrated code to C++ and
  parallelized functions. Added variance shift adjustment, calculation
  of angles with the biological subspace.

- Modified trend specification arguments in trendVar() for greater
  flexibility. Switched from ns() to robustSmoothSpline() to avoid bugs
  with unloaded predict.ns(). Added block= argument for mean/variance
  calculations within each level of a blocking factor.

- Added option to avoid normalization in the SingleCellExperiment
  method for improvedCV2(). Switched from ns() to smooth.spline() or
  robustSmoothSpline() to avoid bugs.

- Replaced zoo functions with runmed() for calculating the median trend
  in DM().

- Added block= argument to correlatePairs() to calculate correlations
  within each level of a blocking factor. Deprecated the use of
  residuals=FALSE for one-way layouts in design=. Preserve input order
  of paired genes in the gene1/gene2 output when pairings!=NULL.

- Added block= argument to overlapExprs() to calculate overlaps within
  each level of a blocking factor. Deprecated the use of
  residuals=FALSE for one-way layouts in design=. Switched to automatic
  ranking of genes based on ability to discriminate between groups.
  Added rank.type= and direction= arguments to control ranking of
  genes.

- Modified combineVar() so that it is aware of the global stats
  recorded in decomposeVar(). Absence of global statistics in the input
  DataFrames now results in an error. Added option to method= to use
  Stouffer's method with residual d.f.-weighted Z-scores. Added
  weighted= argument to allow weighting to be turned off for equal
  batch representation.

- Modified the behaviour of min.mean= in computeSumFactors() when
  clusters!=NULL. Abundance filtering is now performed within each
  cluster and for pairs of clusters, rather than globally.

- Switched to pairwise t-tests in findMarkers(), rather than fitting a
  global linear model. Added block= argument for within-block t-tests,
  the results of which are combined across blocks via Stouffer's
  method. Added lfc= argument for testing against a log-fold change
  threshold. Added log.p= argument to return log-transformed
  p-values/FDRs. Removed empirical Bayes shrinkage as well as the
  min.mean= argument.

- Added the makeTechTrend() function for generating a mean-variance
  trend under Poisson technical noise.

- Added the multiBlockVar() function for convenient fitting of multiple
  mean-variance trends per level of a blocking factor.

- Added the clusterModularity() function for assessing the cluster-wise
  modularity after graph-based clustering.

- Added the parallelPCA() function for performing parallel analysis to
  choose the number of PCs.

- Modified convertT() to return raw counts and size factors for
  CellDataSet output.

- Deprecated exploreData(), selectorPlot() in favour of iSEE().

[SDAMS](https://bioconductor.org/packages/SDAMS)
-----

Version: 0.99.7
Category: package submission
Text:

[SeqArray](https://bioconductor.org/packages/SeqArray)
--------

Changes in version 1.20.0:

UTILITIES

- `seqDigest(f, "annotation/filter")` works on a factor variable

- improve the computational efficiency of `seqMerge()` to avoid
  genotype recompression by padding the 2-bit genotype array in bytes

- significantly improve `seqBlockApply()` (its speed is close to
  `seqApply()`)

- reduce the overhead in `seqSetFilter(, variant.sel=...)`

NEW FEATURES

- `seqGDS2VCF()` outputs a bgzip vcf file for tabix indexing

- two more options "Ultra" and "UltraMax" in `seqStorageOption()`

- '@chrom_rle_val' and '@chrom_rle_len' are added to a GDS file for
  faster chromosome indexing

- new function `seqBCF2GDS()` (requiring the software bcftools)

- new function `seqSetFilterPos()`

- new variable "$dosage_alt" in `seqGetData()` and `seqApply()`

- import VCF files with no GT in `seqVCF2GDS()`

Changes in version 1.18.1-1.18.2:

BUG FIXES

- fix an issue: `seqSetFilterChrom()` extends a genomic range upstream
  and downstream 1bp

- use `.onLoad()` instead of `.onAttach()` to fix
  https://support.bioconductor.org/p/104405/#104443

[seqCAT](https://bioconductor.org/packages/seqCAT)
------

Changes in version 1.2.0:

FEATURES

- Add functionality for analysing VCF files containing unannotated
  variants

- Add functionality for listing non-overlapping variants between
  profiles

- Mitochondrial variants can now be optionally skipped when reading SNV
  profiles in the `read_variants` function

- Add the `list_variants` function for listing the genotypes of
  user-specified variants in each provided SNV profile

- Add the `plot_variant_list` function for plotting a genotype grid for
  each variant output by the `list_variants` function

FIXES

- Fix a multi-sample VCF profile creation issue (python only)

- Reading zero-variant profiles now properly returns a GRanges object
  with a dummy-variant profile containing the sample name

- Enable the `plot_impacts` function to properly analyse multi-impact
  SNVs

- Fix reading of SNV profiles containing single-quoted strings

[SeqVarTools](https://bioconductor.org/packages/SeqVarTools)
-----------

Changes in version 1.17.8:

- Add method alleleCount to return count of alleles.

Changes in version 1.17.7:

- Bug fix in alleleFrequency method for SeqVarData where frequency
  calculation was done twice because the first calculation was not
  returned.

Changes in version 1.17.6:

- Add option to return alternate allele dosage in a sparse matrix using
  the Matrix package.

- Improve speed of reading dosages by using seqBlockApply.

Changes in version 1.17.2:

- Change implementation of iterator classes to identify indices of
  selected variants on object creation and store in variantList slot.
  All iterator classes now extend new class SeqVarIterator.

[sevenC](https://bioconductor.org/packages/sevenC)
------

Changes in version 0.99.0:

- New package "sevenC" for predicting chromatin looping interactions
  from ChIP-seq data and DNA-sequence motifs.

[ShortRead](https://bioconductor.org/packages/ShortRead)
---------

Changes in version 1.37:

BUG FIXES

- (v. 1.37.2) FastqQuality() includes the last printable ASCII
  character '~'

- (v. 1.37.3) countLines() returns numeric values, to allow for files
  with more than 2^31-1 lines

[signeR](https://bioconductor.org/packages/signeR)
------

Changes in version 1.5.2:

- support for generating matrices with a DNAStringSet

Changes in version 1.5.1:

- fix GRanges list usage

- fix conflicts of NMF:seed, by adding VariantAnnotation to Depends
  before NMF

[SIMLR](https://bioconductor.org/packages/SIMLR)
-----

Changes in version 1.5.1 (2018-03-09):

- Added CIMLR implementation.

[SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment)
--------------------

Changes in version 1.1.4:

- Added the clearSpikes() function to remove all spike-in information.

- Added the clearSizeFactors() function to remove all size factor
  information.

- Added the sizeFactorNames() function to query the available (named)
  size factor sets.

- isSpike() with an unknown spike-in set in type= will no longer throw
  an error, and will quietly return NULL.

- isSpike<- with type=NULL is deprecated in favour of clearSpikes() for
  removing existing spike-in information. All spike-in sets must also
  be explicitly named during assignment.

- Added the LinearEmbeddingMatrix class for storing dimensionality
  reduction with loadings.

[singleCellTK](https://bioconductor.org/packages/singleCellTK)
------------

Changes in version 0.99.3:

- Consistent use of camel case throughout package

Changes in version 0.6.3:

- Additional links to help documentation

- Example matrices on upload page.

Changes in version 0.4.7:

- Ability to download/reupload annotation data frame and convert
  annotations to factors/numerics

Changes in version 0.4.5:

- Documentation updates to fix NOTES and pass BiocCheck

[singscore](https://bioconductor.org/packages/singscore)
---------

Changes in version 0.99.9:

- remove internal keyword from generateNull

Changes in version 0.99.8:

- vectorize getPvals

- make toy_epxr_se

Changes in version 0.99.5:

- represent expression data in SummarizedExperiement dataset

- function input checkings

- R code re-organised

Changes in version 0.99.3:

- -change argument name bidirectional to knownDirection, default as
  TRUE

Changes in version 0.99.2:

- -bidirection flag added to simpleScore() function -optimise
  generateNull() function -change GSEABase from import to Depends

[SNPRelate](https://bioconductor.org/packages/SNPRelate)
---------

Changes in version 1.14.0:

- the default compression is "LZMA_RA" in `snpgdsBED2GDS()`,
  `snpgdsVCF2GDS()` and `snpgdsVCF2GDS_R()` for annotations

- support Intel C++ compiler with SSE2/AVX2

- allow interrupt requests in the calculation

- new method options in `snpgdsPairScore()`: GVH.major, GVH.minor,
  GVH.major.only, GVH.minor.only

- force to use integers for 'snp.position' in `snpgdsCreateGeno()`

- unit tests for merging GRMs in `snpgdsMergeGRM()`

- the function `snpgdsSNPListStrand()` is merged to
  `snpgdsSNPListIntersect()`, and it is removed from the package

- update `snpgdsSNPListIntersect()` and `snpgdsCombineGeno()` (work
  correctly)

- replace -INF by NaN in the output of `snpgdsIBDKING()`

Changes in version 1.12.2:

- a new option 'method="Jacquard"' in `snpgdsPairIBD()`

- `snpgdsGRM()` can output the GRM matrix to a GDS file

- a new function `snpgdsMergeGRM()` to merge multiple GRMs

Changes in version 1.12.1:

- fix an issue in the C code 'LENGTH or similar applied to NULL object'

[SparseSignatures](https://bioconductor.org/packages/SparseSignatures)
----------------

Changes in version 1.0.0:

- Package released on Bioconductor in May 2018.

[specL](https://bioconductor.org/packages/specL)
-----

Changes in version 1.13.02 (2017-11-20):

- Eliminated C++-11 and Rcpp linking. Passing unit tests.

- use lower_bound_ function exported by <URL:
  https://cran.r-project.org/package=protViz>

[splatter](https://bioconductor.org/packages/splatter)
--------

Changes in version 1.3.5 (2017-04-25):

- Move scater to Imports and add scater version

- Remove lingering references to SCESets

- Add option to use a normal distribution for library sizes in Splat
  simulations

- Allow Splat dropout parameters to be specified by experiment, batch,
  group or cell

- Add SparseDC simulation

- Rename params in metadata slot of simulation to Params for
  consistency

- Improve and colourise Params print output

- Improve test coverage

- Various other minor updates and bug fixes

[SPONGE](https://bioconductor.org/packages/SPONGE)
------

Changes in version 1.0.3:

- Smaller bug fixes of broken unit tests through changes in 1.0.2

Changes in version 1.0.2:

- Big memory regularly breaks the build. Thus it is no longer used by
  default. big memory objects can still be created manually and will be
  used in sponge

Changes in version 1.0.1:

- Release with Bioconductor 3.6

- Fix of bug caused by new version of bigmemory

[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
--------------------

Changes in version 1.10.0:

NEW FEATURES

- Add "subset" method for SummarizedExperiment objects. See
  https://github.com/Bioconductor/SummarizedExperiment/pull/6

- rowRanges() now is supported on a SummarizedExperiment object that is
  not a RangedSummarizedExperiment, and returns NULL. Also doing
  'rowRanges(x) <- NULL' on a RangedSummarizedExperiment object now is
  supported and degrades it to a SummarizedExperiment instance.

- Add 'BACKEND' argument to "realize" method for SummarizedExperiment
  objects.

SIGNIFICANT USER-VISIBLE CHANGES

- saveHDF5SummarizedExperiment() and loadHDF5SummarizedExperiment() are
  now in the HDF5Array package.

- Replace old "updateObject" method for SummarizedExperiment objects
  with a new one. The new method calls updateObject() on all the assays
  of the object. This will update SummarizedExperiment objects (and
  their derivatives like BSseq objects) that have "old" DelayedArray
  objects in their assays. The old method has been around since BioC
  3.2 (released 2.5 years ago) and was used to update objects made
  prior to the change of internals that happened between BioC 3.1 and
  BioC 3.2. All these "old" objects should have been updated by now so
  we don't need this anymore.

BUG FIXES

- Modify the "[<-" method for SummarizedExperiment to leave
  'metadata(x)' intact instead of trying to combine it with
  'metadata(value)'. With this change 'x&#91;i , j&#93; <- x&#91;i , j&#93;' behaves
  like a no-op (as expected) instead of duplicating metadata(x).

- The SummarizedExperiment() constructor does not try to downgrade the
  supplied rowData and/or colData to DataFrame anymore if they derive
  from DataFrame.

[SWATH2stats](https://bioconductor.org/packages/SWATH2stats)
-----------

Changes in version 1.9.3:

BUG FIXES

- removed aLFQ direct link in manual

Changes in version 1.9.2:

NEW FEATURES

- add option check_tranisitions to convert_aLFQ function

- replace align_orig_filename column with filename

Changes in version 1.9.1:

NEW FEATURES

- SWATH2stats in BioC 3.7 development release

Changes in version 1.8.1:

NEW FEATURES

- SWATH2stats in BioC 3.6 release

[synapter](https://bioconductor.org/packages/synapter)
--------

Changes in version 2.3.1:

- Partly revert 95f4094 because MSnbase:::utils.applyColumnwiseByGroup
  is gone.

- Use `html_document` instead of `html_document2` in the vignettes
  &#91;2018-01-17&#93;. # Synapter 2.1

[TargetSearch](https://bioconductor.org/packages/TargetSearch)
------------

Changes in version 1.36.0:

NEW FEATURES

- New dataset object TSExample. This dataset contains data that used to
  be stored in package TargetSearchData.

- New low-level function to search peaks (FindAllPeaks). This allows
  advanced users to refine peak-searches.

- New function to plot peaks across samples (plotPeakRI). Used for
  quality checks of peak annotation and fine-tunning search parameters.

SIGNIFICANT USER-VISIBLE CHANGES

- Add extra checks when manipulating tsLib objects. Extra care needs to
  be taken if changes to the quant/selective/top masses are changed.

- Sample IDs (names) must be unique. These might generate errors when
  loading old TargetSearch workspaces.

BUG FIXES

- tsLib: ensure that every slot in the object contain a library ID.

- Fix warnings during R CMD check.

- Refactor C code for finding peaks. It is possible to return all peaks
  instead of only the most abundant. No visibles changes for the end
  user.

- Refactor C code for NetCDF manipulation and peak finding to reduce
  code duplication. No visible changes for the end user.

- General R code houskeeping: Removal of mixed tabs and spaces, fix
  tabulation, add Rbuildignore.

[TCGAbiolinks](https://bioconductor.org/packages/TCGAbiolinks)
------------

Changes in version 2.7.13:

- Adding new function: PanCancerAtlas_subtypes

- Updating DNA methylation probe information function

- Start to update vignette

- FPPE information is being added in GDCprepare

- Minor issue fixes

[TEQC](https://bioconductor.org/packages/TEQC)
----

Version: 4.0.1
Category: fixed problem with IRanges version >=2.13.28 in
        coverage.target
Text:

Version: 4.0.0
Text: package is now based on 'GRanges' objects instead of deprecated
        'RangedData' objects

Version: 4.0.0
Text: vignette is now created with knitr

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

[ToPASeq](https://bioconductor.org/packages/ToPASeq)
-------

Changes in version 1.13.2:

- Package maintainer change: new maintainer is Ludwig Geistlinger

[topdownr](https://bioconductor.org/packages/topdownr)
--------

Changes in version 1.1.7:

- Add `mz,FragmentViews-method` &#91;2018-02-01&#93;.

- Remove internal `fragmentMass` and `fragmentNames` functions
  &#91;2018-02-22&#93;.

- Parse "spectrumId" column of the mzML header to find the scan number
  (instead of the "acquisitionNum") because ProteomDiscover generates
  non-standard "spectrumId" and proteowizard fails to translated it
  into a valid "acquisitionNum". See #73 for details &#91;2018-02-22&#93;.

- Recalculate TotIonCurrent in the main loop of `.readMzMl`
  &#91;2018-02-22&#93;.

- Add `FragmentCoverage` and `BondCoverage` columns to
  `bestConditions,NCBSet-method` &#91;2018-02-23&#93;.

- Use retention times to test for correct matching between ScanHeadsman
  .txt output and mzML files; closes #74; &#91;2018-02-23&#93;.

Changes in version 1.1.6:

- Rotate fragment labels (vertical orientation) in `plot` &#91;2018-01-17&#93;.

- Replace signature for `updateMedianInjectionTime,TopDownSet-method`
  to `updateMedianInjectionTime,AbstractTopDownSet-method`; closes #69;
  see also #71 &#91;2018-01-27&#93;.

- Fix `.matchFragments` for `length(fmass) == 0` &#91;2018-01-27&#93;.

- Just plot fragments that are present in current `TopDownSet` see #70
  &#91;2018-01-27&#93;.

- Add `combine,FragmentViews,FragmentViews-method` &#91;2018-01-27&#93;.

- Allow to `combine` `TopDownSet` objects with different fragment
  types; closes #71 &#91;2018-01-27&#93;.

- Add `all.equal` for `AbstractTopDownSet` objects &#91;2018-01-27&#93;.

- Allow the user to decide how to handle redundant fragment matching.
  Current default is `redundantFragmentMatch="remove"` and
  `redundantIonMatch="remove"`. This will reduce the number of fragment
  matches. Choose `"closest"` for both to get the old behaviour. See
  also #72 &#91;2018-01-29&#93;.

- `TopDownSet` object store the matching `tolerance` and strategies
  (`redundantIonMatch`, `redundantFragmentMatch`). `AbstractTopDownSet`
  and `NCBSet` lost their `tolerance` slot. Saved objects need to be
  recreated &#91;2018-01-30&#93;.

- `bestConditions,NCBSet-method` returns a 5-column matrix now. Colums
  are: Index, FragmentsAddedToCombination, BondsAddedToCombination,
  FragmentsInCondition, BondsInCondition; see #52 &#91;2018-01-30&#93;.

Changes in version 1.1.5:

- Keep full filename (before `basename` was used) in
  `AbstractTopDownSet` objects &#91;2017-12-28&#93;.

- Add `plot,TopDownSet-method` &#91;2017-12-29&#93;.

- `bestConditions,NCBSet-method` gains a new argument `maximise` that
  allows to optimise for number of fragments or bonds covered (default:
  `"fragments"`); see #52 &#91;2018-01-15&#93;.

Changes in version 1.1.4:

- Add missing export of `combine` and documentation &#91;2017-12-28&#93;.

- Resave `tds` example data set to reflect changes in `colData`
  introduced in version 1.1.2 &#91;2017-12-28&#93;.

Changes in version 1.1.3:

- Add `conditionNames,AbstractTopDownSet-method` to access
  `rownames(colData(tds))` &#91;2017-12-23&#93;.

- Add `updateConditionNames,AbstractTopDownSet-method` (closes #60)
  &#91;2017-12-23&#93;.

- Turn `updateMedianInjectionTime,TopDownSet-method` into
  `updateMedianInjectionTime,AbstractTopDownSet-method` to work with
  `TopDownSet` and `NCBSet` objects &#91;2017-12-27&#93;.

- Add `combine,AbstractTopDownSet-method` to combine multiple
  `TopDownSet`/`NCBSet` objects (closes #69) &#91;2017-12-28&#93;.

Changes in version 1.1.2:

- Add `.rbind` to combine scan and method information with different
  number of colums (could happen when CID/HCD and UVPD scans are taken
  independently with different software versions) &#91;2017-12-22&#93;.

- Don't replace NA values with zeros in the `colData` &#91;2017-12-22&#93;.

- Convert On/Off `character` columns in scan and method information to
  `logical` &#91;2017-12-22&#93;.

- Fix `.camelCase` to avoid "TIC" to "TIc" and
  "UseCalibratedUVPDTimeMs2" to "UseCalibrateduvpdTimems2" conversion
  (now: "Tic" and "UseCalibratedUvpdTimeMs2") &#91;2017-12-22&#93;.

Changes in version 1.1.1:

- Respect assigned intensity in conditions for
  `bestConditions,NCBSet-method` and `fragmentationMap` (closes #62)
  &#91;2017-12-02&#93;.

- Fix explanation of random forest barchart in analysis vignette
  &#91;2017-12-02&#93;.

- Create all fragmentation methods in `.readScanHeadsTable` to avoid
  error if any is missing (fixes #68) &#91;2017-12-20&#93;.

- Never remove Activation column in `colData` (even not if
  `readTopDownFiles(..., dropNonInformativeColumns=TRUE)`)
  &#91;2017-12-20&#93;.

- Allow UVPD in `fragmentationMap,NCBSet-method` &#91;2017-12-20&#93;.

- Add new method: `updateMedianInjectionTime,TopDownSet-method` (closes
  #66) &#91;2017-12-20&#93;.

Changes in version 1.1.0:

- New version for Bioc 3.7 (devel) # topdownr 1.0

[trackViewer](https://bioconductor.org/packages/trackViewer)
-----------

Changes in version 1.15.13:

- add feature.gr parameter to plotGInteractions.

Changes in version 1.15.12:

- add new function plotGInteractions.

Changes in version 1.15.11:

- update the title

Changes in version 1.15.10:

- Change all "class" function to "is" function.

Changes in version 1.15.9:

- Add lollipopData to tracks-class.

Changes in version 1.15.8:

- Plot lollipop plot with tracks.

Changes in version 1.15.7:

- Add new parameter rescale to lolliplot.

Changes in version 1.15.6:

- Add new function geneTrack.

Changes in version 1.15.5:

- fix the warning in windows documentation.

Changes in version 1.15.4:

- fix the warning in windows documentation.

Changes in version 1.15.3:

- fix the position of ylabs in plotOneIdeo.

Changes in version 1.15.2:

- fix the bug in addGuideLine and viewTracks when plot tracks with
  breaks.

Changes in version 1.15.1:

- fix the bug in dandelion.plot

- use Roxygen2 to generate help documents

- add new functions: getLocation, viewGene

[transcriptogramer](https://bioconductor.org/packages/transcriptogramer)
-----------------

Changes in version 1.1.18:

- The plot of the differentiallyExpressed() method shows the number of
  clusters detected and is created by the ggplot2 package.

Changes in version 1.1.17:

- Improved documentation.

Changes in version 1.1.7:

- Added a new argument to the differentiallyExpressed() method,
  supporting the limma-trend approach for RNA-Seq.

Changes in version 1.1.6:

- Updated GPL570 dataset.

Changes in version 1.1.3:

- The datasets of the four species were improved to obtain a better
  clustering.

[treeio](https://bioconductor.org/packages/treeio)
------

Changes in version 1.3.9:

- Exporter.Rmd vignette <2017-12-13, Wed>

Changes in version 1.3.8:

- mv treeio.Rmd vignette to Importer.Rmd and update the contents
  <2017-12-13, Wed>

- write.beast for treedata object <2017-12-12, Tue>

- add "connect" parameter in groupOTU <2017-12-12, Tue> +
  https://groups.google.com/forum/#!msg/bioc-ggtree/Q4LnwoTf1DM/yEe95OFfCwAJ

Changes in version 1.3.7:

- export groupClade.phylo method <2017-12-11, Mon>

Changes in version 1.3.6:

- re-defined groupOTU and groupClade generic using S3 <2017-12-11, Mon>

Changes in version 1.3.5:

- parent, ancestor, child, offspring, rootnode and sibling generic and
  method for phylo <2017-12-11, Mon>

- update mask and merge_tree function according to the treedata object
  <2017-12-11, Mon>

Changes in version 1.3.4:

- support tbl_tree object defined in tidytree <2017-12-08, Fri>

Changes in version 1.3.3:

- read.codeml output treedata, remove codeml class and clean up code
  <2017-12-07, Thu>

Changes in version 1.3.2:

- read.codeml_mlc output treedata object and remove codeml_mlc class
  <2017-12-06, Wed>

- read.paml_rst output treedata and remove paml_rst class <2017-12-06,
  Wed>

- read.phylip.tree and read.phylip.seq

- read.phylip output treedata object and phylip class definition was
  removed

- read.hyphy output treedata object; hyphy class definition was removed

- remove r8s class, read.r8s now output multiPhylo object

- jplace class inherits treedata <2017-12-05, Tue>

- using treedata object to store beast and mrbayes tree

- export read.mrbayes

Changes in version 1.3.1:

- compatible to parse beast output that only contains HPD range
  <2017-11-01, Wed> +
  https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/bioc-ggtree/RF2Ly52U_gc/jEP97nNPAwAJ

[tRNAscanImport](https://bioconductor.org/packages/tRNAscanImport)
--------------

Changes in version 0.99.21 (2018-03-17):

Changes

- Initial version

[TVTB](https://bioconductor.org/packages/TVTB)
----

Changes in version 1.5.4 (2017-09-24):

Bug fix

- min(ranges(GRanges)) does not seem to work anymore, replaced by more
  explicit expression.

- Commented out some unit tests until subsetting of FilterRules by row
  _and_ column throws an error again.

[tximport](https://bioconductor.org/packages/tximport)
--------

Changes in version 1.8.0:

- Added support for StringTie output.

[TxRegInfra](https://bioconductor.org/packages/TxRegInfra)
----------

Version: 0.99.10
Category: NEW FEATURES
Text:

[Uniquorn](https://bioconductor.org/packages/Uniquorn)
--------

Changes in version 1.99.3 (2018-01-03):

Bioconductor compliance

- Minor modifications to comply with Bioconductor regulations

Changes in version 1.99.2 (2018-20-02):

Added RNA-seq and panel-seq capability

- Some more minor bugfixes

Changes in version 1.99.1 (2018-12-02):

Added RNA-seq and panel-seq capability

- Minor bugfixes

Changes in version 1.99.0 (2018-15-01):

Added RNA-seq and panel-seq capability

- Modified the package to work with large scale RNA-seq and small scale
  panel-seq data

- Please note that this might break custom libraries. In case of
  questions please contact the author

[variancePartition](https://bioconductor.org/packages/variancePartition)
-----------------

Changes in version 1.9.9:

- Fix issue when package is autoloaded when starting R

Changes in version 1.9.8:

- Fix issue where if info data.frame contained a column name "gene",
  fitExtractVarPartModel() would not run

Changes in version 1.9.6:

- Fix tximport issue with eval=FALSE

Changes in version 1.9.2:

- Fix vignette

Changes in version 1.9.1:

- Fix formatting of vignette

- add description of canCorPairs() function

[VariantFiltering](https://bioconductor.org/packages/VariantFiltering)
----------------

Changes in version 1.16:

USER VISIBLE CHANGES

- Specific filtering functions have been deprecated in favor of a more
  general filtering mechanism.

- Shiny-app has been rewritten.

[xcms](https://bioconductor.org/packages/xcms)
----

Changes in version 3.1.3:

BUG FIXES

- Fix misplaced parenthesis in the check for multiple spectra in
  findChromPeaks,OnDiskMSnExp,MSWParam. Thanks to @RonanDaly (PR #276).

- Update link to correct metlin page in diffreport result (issue #204).

Changes in version 3.1.2:

NEW FEATURES

- Add filterFeatureDefinitions function.

BUG FIXES

- Fix #273: better error message in case not a single feature could be
  defined by groupChromPeaks.

Changes in version 3.1.1:

NEW FEATURES

- Reading raw files using xcmsSet or xcmsRaw uses now the automatic
  file type detection feature from mzR.

- c function to concatenate XCMSnExp objects.

- groupnames method for XCMSnExp objects (issue #250).

BUG FIXES

- Fix #237: findPeaks.MSW was not throwing an error if applied to
  multi-spectrum MS file.

- Fix #249: quantile call in adjustRtime PeakGroups without na.rm =
  TRUE.

- Fix #259

[xps](https://bioconductor.org/packages/xps)
---

VERSION xps-1.37.2

- configure.in file - unix line endings

VERSION xps-1.37.1

- update INSTALL and README file

Changes in version 3.3:

[zinbwave](https://bioconductor.org/packages/zinbwave)
--------

Changes in version 1.1.6 (2018-04-17):

- `zinbwave` now uses `counts` assay by default.

- Users can now specify which assay to use to fit the zinb model.

Changes in version 1.1.5 (2018-02-15):

- Computational weights are computed in `zinbwave` as saved as assay.

- Modified vignette to include example of Differential Expression.

- Improved documentation for `zinbwave`.

NEWS from new and existing Data Experiment  Packages
===================================

[ASICSdata](https://bioconductor.org/packages/ASICSdata)
---------

Version: 0.99.2
Category: Fixed a few typos
Text:

Version: 0.99.1
Category: Updated R version dependency from 3.4 to 3.5
Text:

Version: 0.99.0
Category: First submitted version to Bioconductor
Text:

[chipenrich.data](https://bioconductor.org/packages/chipenrich.data)
---------------

Changes in version 2.4.0:

- Updated locus definitions based on TxDb 3.4.0 and OrgDb 3.5.0
  packages

- Remove gene expression and EHMN gene sets.

[dsQTL](https://bioconductor.org/packages/dsQTL)
-----

Changes in version 2.17:

USER VISIBLE CHANGES

- ch2locs (retrievable via dsQTL::getSNPlocs) has been changed at about
  1850 locations where rs numbers had been associated with hg19
  addresses; the dsQTL regions are hg18 as are all the chr2... SNP
  addresses.  Previously the discoverable rs numbers used in the
  Chicago distribution from
  http://eqtl.uchicago.edu/dsQTL_data/GENOTYPES/ had be mapped via
  SNPlocs...20111119, but now they come directly from the Chicago text
  file.

[flowPloidyData](https://bioconductor.org/packages/flowPloidyData)
--------------

Changes in version 1.5.3 (2018-01-19):

User Visible Changes

- two additional FCS files have been added: fpVac, which is useful in
  demonstrating gating; and fpBad, which is a very low quality file
  used in unit tests.

[HDCytoData](https://bioconductor.org/packages/HDCytoData)
----------

Version: 1.0.0
Text:

[mCSEAdata](https://bioconductor.org/packages/mCSEAdata)
---------

Changes in version 0.99.0:

- mCSEAdata release.

[MetaGxBreast](https://bioconductor.org/packages/MetaGxBreast)
------------

Changes in version 0.99.0:

NEW FEATURES

- This is the first iteration of the MetaGxBreast package. The package
  contains 39 breast cancer Esets

[MetaGxOvarian](https://bioconductor.org/packages/MetaGxOvarian)
-------------

Changes in version 0.99.0:

NEW FEATURES

- This is the first iteration of the MetaGxOvarian package. It is
  essentially an updated version of the curatedOvarianData package with
  5 additional datasets

[MetaGxPancreas](https://bioconductor.org/packages/MetaGxPancreas)
--------------

Changes in version 0.99.0:

NEW FEATURES

- This is the first iteration of the MetaGxPancreas package.

[microRNAome](https://bioconductor.org/packages/microRNAome)
-----------

Changes in version 1.1.2:

- Fix unit test error.

[mtbls2](https://bioconductor.org/packages/mtbls2)
------

Changes in version 1.9.1:

- Add vignette based on xcms3 new interface, bump minimum xcms version
  to 3.0.0

[pRolocdata](https://bioconductor.org/packages/pRolocdata)
----------

Changes in version 1.17.4:

- Typo in tlopt colnames fixed for f1000 workflow <2018-04-10 Tue>

Changes in version 1.17.3:

- Add Proteome biocView <2018-04-05 Thu>

Changes in version 1.17.2:

- Add data from Beltran et al. 2016 <2018-03-17 Sat>

- Add data from Hirst et al. 2018 <2018-03-18 Sun>

- Add data from Itzhal et al. 2017 <2018-03-20 Tue>

Changes in version 1.17.1:

- Update README with information on how to add data <2018-03-08 Thu>

[PtH2O2lipids](https://bioconductor.org/packages/PtH2O2lipids)
------------

Changes in version 2016-04-21:

- Initial release for Bioconductor

[RcisTarget.hg19.motifDBs.cisbpOnly.500bp](https://bioconductor.org/packages/RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
----------------------------------------

Changes in version 0.99.7:

- Motif rankings are now transposed (motifs in rows, genes in columns),
  and updated to motif collection 9.

- Motif-TF annotation databases have been re-formatted for easier
  selection of the direct/inferred annotations.

[RforProteomics](https://bioconductor.org/packages/RforProteomics)
--------------

Changes in version 1.17.1:

- Fix/update vignettes <2018-01-12 Fri>

- Add mention of DEP package <2018-01-13 Sat>

Changes in version 1.17.0:

- New Bioc devel version

Changes in version 1.16.0:

- New Bioc release version

[TargetSearchData](https://bioconductor.org/packages/TargetSearchData)
----------------

Changes in version 1.17.2:

SIGNIFICANT USER-VISIBLE CHANGES

- The dataset `TargetSearchData` which used to contain `TargetSearch`
  objects has been moved from package `TargetSearchData` to package
  `TargetSearch`. Note that the dataset name has been renamed to
  `TSExample`. To load it use data(TSExample) within `TargetSearch`.

[tissueTreg](https://bioconductor.org/packages/tissueTreg)
----------

Changes in version 0.99.3:

- changing unit tests

Changes in version 0.99.2:

- changing man page example

Changes in version 0.99.1:

- changes to comply with Bioconductor package guide lines

Changes in version 0.99.0:

- First release

[topdownrdata](https://bioconductor.org/packages/topdownrdata)
------------

Changes in version 1.1.1:

- Add UVPD data for CA.

Changes in version 1.1.0:

- New version for Bioc 3.7 (devel) # topdownrdata 1.0

Deprecated and Defunct Packages
===============================

Nine software packages were removed from this release (after being
deprecated in BioC 3.6): BioMedR, ddgraph, EWCE, HCsnip, stepwiseCM,
domainsignatures, iontree, oneChannelGUI, RCytoscape. One software 
package was removed from Bioconductor at user request (mvGST) without
previous deprecation. 

Ten software packages (ontoCat, spliceR, GMRP, MBttest, OperaMate,
DASC, htSeqTools, PAnnBuilder, phenoDist, BrowserVizDemo) are deprecated
in this release and will be removed in BioC 3.8.

One annotation data package (IlluminaHumanMethylation450k.db) was 
removed from this release.

Three experimental data packages (RnaSeqTutorial, cheung2010, MEALData) 
are deprecated in this release and will be removed in Bioc 3.8. 

