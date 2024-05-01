May 1, 2024

**Bioconductor:**

We are pleased to announce Bioconductor 3.19, consisting of
2300 software packages, 430 experiment data packages, 926 annotation
packages, 30 workflows and 5 books.

There are 80 new software packages, 10 new data experiment packages,
5 new annotation packages, no new workflows, 1 new book, and many updates and
improvements to existing packages.

Bioconductor 3.19 is compatible with R 4.4, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher, macOS arm64 and Linux
arm64. This release will also include updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.19](#getting-started-with-bioconductor-319)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [NEWS from existing workflows](#news-from-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.19
======================================

To update to or install Bioconductor 3.19

1. Install R 4.4. Bioconductor 3.19 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 80 new software packages in this release of Bioconductor.

- [AlphaMissenseR](/packages/AlphaMissenseR) The AlphaMissense
  publication
  <https://www.science.org/doi/epdf/10.1126/science.adg7492> outlines
  how a variant of AlphaFold / DeepMind was used to predict missense
  variant pathogenicity. Supporting data on Zenodo
  <https://zenodo.org/record/10813168> include, for instance, 71M
  variants across hg19 and hg38 genome builds. The 'AlphaMissense'
  package allows ready access to the data, downloading individual
  files to DuckDB databases for exploration and integration into *R*
  and *Bioconductor* workflows.

- [Banksy](/packages/Banksy) Banksy is an R package that incorporates
  spatial information to cluster cells in a feature space (e.g. gene
  expression). To incorporate spatial information, BANKSY computes
  the mean neighborhood expression and azimuthal Gabor filters that
  capture gene expression gradients. These features are combined with
  the cell's own expression to embed cells in a neighbor-augmented
  product space which can then be clustered, allowing for accurate
  and spatially-aware cell typing and tissue domain segmentation.

- [BERT](/packages/BERT) Provides efficient batch-effect adjustment
  of data with missing values. BERT orders all batch effect
  correction to a tree of pairwise computations. BERT allows
  parallelization over sub-trees.

- [betaHMM](/packages/betaHMM) A novel approach utilizing a
  homogeneous hidden Markov model. And effectively model
  untransformed beta values. To identify DMCs while considering the
  spatial. Correlation of the adjacent CpG sites.

- [bettr](/packages/bettr) bettr provides a set of interactive
  visualization methods to explore the results of a benchmarking
  study, where typically more than a single performance measures are
  computed. The user can weight the performance measures according to
  their preferences. Performance measures can also be grouped and
  aggregated according to additional annotations.

- [biocroxytest](/packages/biocroxytest) This package provides a
  roclet for roxygen2 that identifies and processes code blocks in
  your documentation marked with `@longtests`. These blocks should
  contain tests that take a long time to run and thus cannot be
  included in the regular test suite of the package. When you run
  `roxygen2::roxygenise` with the `longtests_roclet`, it will extract
  these long tests from your documentation and save them in a
  separate directory. This allows you to run these long tests
  separately from the rest of your tests, for example, on a
  continuous integration server that is set up to run long tests.

- [BREW3R.r](/packages/BREW3R.r) This R package provide functions
  that are used in the BREW3R workflow. This mainly contains a
  function that extend a gtf as GRanges using information from
  another gtf (also as GRanges). The process allows to extend gene
  annotation without increasing the overlap between gene ids.

- [CaMutQC](/packages/CaMutQC) CaMutQC is able to filter false
  positive mutations generated due to technical issues, as well as to
  select candidate cancer mutations through a series of
  well-structured functions by labeling mutations with various flags.
  And a detailed and vivid filter report will be offered after
  completing a whole filtration or selection section. Also, CaMutQC
  integrates serveral methods and gene panels for Tumor Mutational
  Burden (TMB) estimation.

- [ClustAll](/packages/ClustAll) Data driven strategy to find hidden
  groups of patients with complex diseases using clinical data.
  ClustAll facilitates the unsupervised identification of multiple
  robust stratifications. ClustAll, is able to overcome the most
  common limitations found when dealing with clinical data (missing
  values, correlated data, mixed data types).

- [ClusterFoldSimilarity](/packages/ClusterFoldSimilarity) This
  package calculates a similarity coefficient using the fold changes
  of shared features (e.g. genes) among clusters of different
  samples/batches/datasets. The similarity coefficient is calculated
  using the dot-product (Hadamard product) of every pairwise
  combination of Fold Changes between a source cluster i of
  sample/dataset n and all the target clusters j in sample/dataset m

- [CRISPRball](/packages/CRISPRball) A Shiny application for
  visualization, exploration, comparison, and filtering of CRISPR
  screens analyzed with MAGeCK RRA or MLE. Features include
  interactive plots with on-click labeling, full customization of
  plot aesthetics, data upload and/or download, and much more.
  Quickly and easily explore your CRISPR screen results and generate
  publication-quality figures in seconds.

- [crisprShiny](/packages/crisprShiny) Provides means to
  interactively visualize guide RNAs (gRNAs) in GuideSet objects via
  Shiny application. This GUI can be self-contained or as a module
  within a larger Shiny app. The content of the app reflects the
  annotations present in the passed GuideSet object, and includes
  intuitive tools to examine, filter, and export gRNAs, thereby
  making gRNA design more user-friendly.

- [CTexploreR](/packages/CTexploreR) The CTexploreR package
  re-defines the list of Cancer Testis/Germline (CT) genes. It is
  based on publicly available RNAseq databases (GTEx, CCLE and TCGA)
  and summarises CT genes' main characteristics. Several
  visualisation functions allow to explore their expression in
  different types of tissues and cancer cells, or to inspect the
  methylation status of their promoters in normal tissues.

- [cypress](/packages/cypress) CYPRESS is a cell-type-specific power
  tool. This package aims to perform power analysis for the
  cell-type-specific data. It calculates FDR, FDC, and power, under
  various study design parameters, including but not limited to
  sample size, and effect size. It takes the input of a
  SummarizeExperimental(SE) object with observed mixture data
  (feature by sample matrix), and the cell-type mixture proportions
  (sample by cell-type matrix). It can solve the cell-type mixture
  proportions from the reference free panel from TOAST and conduct
  tests to identify cell-type-specific differential expression (csDE)
  genes.

- [CytoMDS](/packages/CytoMDS) This package implements a low
  dimensional visualization of a set of cytometry samples, in order
  to visually assess the 'distances' between them. This, in turn, can
  greatly help the user to identify quality issues like batch effects
  or outlier samples, and/or check the presence of potential sample
  clusters that might align with the exeprimental design. The CytoMDS
  algorithm combines, on the one hand, the concept of Earth Mover's
  Distance (EMD), a.k.a. Wasserstein metric and, on the other hand,
  the Multi Dimensional Scaling (MDS) algorithm for the low
  dimensional projection. Also, the package provides some diagnostic
  tools for both checking the quality of the MDS projection, as well
  as tools to help with the interpretation of the axes of the
  projection.

- [Damsel](/packages/Damsel) Damsel provides an end to end analysis
  of DamID data. Damsel takes bam files from Dam-only control and
  fusion samples and counts the reads matching to each GATC region.
  edgeR is utilised to identify regions of enrichment in the fusion
  relative to the control. Enriched regions are combined into peaks,
  and are associated with nearby genes. Damsel allows for IGV style
  plots to be built as the results build, inspired by ggcoverage, and
  using the functionality and layering ability of ggplot2. Damsel
  also conducts gene ontology testing with bias correction through
  goseq, and future versions of Damsel will also incorporate motif
  enrichment analysis. Overall, Damsel is the first package allowing
  for an end to end analysis with visual capabilities. The goal of
  Damsel was to bring all the analysis into one place, and allow for
  exploratory analysis within R.

- [dar](/packages/dar) Differential abundance testing in microbiome
  data challenges both parametric and non-parametric statistical
  methods, due to its sparsity, high variability and compositional
  nature. Microbiome-specific statistical methods often assume
  classical distribution models or take into account compositional
  specifics. These produce results that range within the specificity
  vs sensitivity space in such a way that type I and type II error
  that are difficult to ascertain in real microbiome data when a
  single method is used. Recently, a consensus approach based on
  multiple differential abundance (DA) methods was recently suggested
  in order to increase robustness. With dar, you can use dplyr-like
  pipeable sequences of DA methods and then apply different consensus
  strategies. In this way we can obtain more reliable results in a
  fast, consistent and reproducible way.

- [DegCre](/packages/DegCre) DegCre generates associations between
  differentially expressed genes (DEGs) and cis-regulatory elements
  (CREs) based on non-parametric concordance between differential
  data. The user provides GRanges of DEG TSS and CRE regions with
  differential p-value and optionally log-fold changes and DegCre
  returns an annotated Hits object with associations and their
  calculated probabilities. Additionally, the package provides
  functionality for visualization and conversion to other formats.

- [DeProViR](/packages/DeProViR) Emerging infectious diseases,
  exemplified by the zoonotic COVID-19 pandemic caused by SARS-CoV-2,
  are grave global threats. Understanding protein-protein
  interactions (PPIs) between host and viral proteins is essential
  for therapeutic targets and insights into pathogen replication and
  immune evasion. While experimental methods like yeast two-hybrid
  screening and mass spectrometry provide valuable insights, they are
  hindered by experimental noise and costs, yielding incomplete
  interaction maps. Computational models, notably DeProViR, predict
  PPIs from amino acid sequences, incorporating semantic information
  with GloVe embeddings. DeProViR employs a Siamese neural network,
  integrating convolutional and Bi-LSTM networks to enhance accuracy.
  It overcomes the limitations of feature engineering, offering an
  efficient means to predict host-virus interactions, which holds
  promise for antiviral therapies and advancing our understanding of
  infectious diseases.

- [dinoR](/packages/dinoR) dinoR tests for significant differences in
  NOMe-seq footprints between two conditions, using genomic regions
  of interest (ROI) centered around a landmark, for example a
  transcription factor (TF) motif. This package takes NOMe-seq data
  (GCH methylation/protection) in the form of a Ranged Summarized
  Experiment as input. dinoR can be used to group sequencing
  fragments into 3 or 5 categories representing characteristic
  footprints (TF bound, nculeosome bound, open chromatin), plot the
  percentage of fragments in each category in a heatmap, or averaged
  across different ROI groups, for example, containing a common TF
  motif. It is designed to compare footprints between two sample
  groups, using edgeR's quasi-likelihood methods on the total
  fragment counts per ROI, sample, and footprint category.

- [epiregulon](/packages/epiregulon) Gene regulatory networks model
  the underlying gene regulation hierarchies that drive gene
  expression and observed phenotypes. Epiregulon infers TF activity
  in single cells by constructing a gene regulatory network
  (regulons). This is achieved through integration of scATAC-seq and
  scRNA-seq data and incorporation of public bulk TF ChIP-seq data.
  Links between regulatory elements and their target genes are
  established by computing correlations between chromatin
  accessibility and gene expressions.

- [epiregulon.extra](/packages/epiregulon.extra) Gene regulatory
  networks model the underlying gene regulation hierarchies that
  drive gene expression and observed phenotypes. Epiregulon infers TF
  activity in single cells by constructing a gene regulatory network
  (regulons). This is achieved through integration of scATAC-seq and
  scRNA-seq data and incorporation of public bulk TF ChIP-seq data.
  Links between regulatory elements and their target genes are
  established by computing correlations between chromatin
  accessibility and gene expressions.

- [faers](/packages/faers) The FDA Adverse Event Reporting System
  (FAERS) is a database used for the spontaneous reporting of adverse
  events and medication errors related to human drugs and therapeutic
  biological products. faers pacakge serves as the interface between
  the FAERS database and R. Furthermore, faers pacakge offers a
  standardized approach for performing pharmacovigilance analysis.

- [findIPs](/packages/findIPs) Feature rankings can be distorted by a
  single case in the context of high-dimensional data. The cases
  exerts abnormal influence on feature rankings are called
  influential points (IPs). The package aims at detecting IPs based
  on case deletion and quantifies their effects by measuring the rank
  changes (DOI:10.48550/arXiv.2303.10516). The package applies a
  novel rank comparing measure using the adaptive weights that stress
  the top-ranked important features and adjust the weights to ranking
  properties.

- [GeDi](/packages/GeDi) The package provides different distances
  measurements to calculate the difference between genesets. Based on
  these scores the genesets are clustered and visualized as graph.
  This is all presented in an interactive Shiny application for easy
  usage.

- [ggtreeSpace](/packages/ggtreeSpace) This package is a
  comprehensive visualization tool specifically designed for
  exploring phylomorphospace. It not only simplifies the process of
  generating phylomorphospace, but also enhances it with the
  capability to add graphic layers to the plot with grammar of
  graphics to create fully annotated phylomorphospaces. It also
  provide some utilities to help interpret evolutionary patterns.

- [ginmappeR](/packages/ginmappeR) Provides functionalities to
  translate gene or protein identifiers between state-of-art
  biological databases: CARD (<https://card.mcmaster.ca/>), NCBI
  Protein, Nucleotide and Gene (<https://www.ncbi.nlm.nih.gov/>),
  UniProt (<https://www.uniprot.org/>) and KEGG
  (<https://www.kegg.jp>). Also offers complementary functionality
  like NCBI identical proteins or UniProt similar genes clusters
  retrieval.

- [gINTomics](/packages/gINTomics) gINTomics is an R package for
  Multi-Omics data integration and visualization. gINTomics is
  designed to detect the association between the expression of a
  target and of its regulators, taking into account also their
  genomics modifications such as Copy Number Variations (CNV) and
  methylation. What is more, gINTomics allows integration results
  visualization via a Shiny-based interactive app.

- [GrafGen](/packages/GrafGen) To classify Helicobacter pylori
  genomes according to genetic distance from nine reference
  populations. The nine reference populations are hpgpAfrica,
  hpgpAfrica-distant, hpgpAfroamerica, hpgpEuroamerica,
  hpgpMediterranea, hpgpEurope, hpgpEurasia, hpgpAsia, and
  hpgpAklavik86-like. The vertex populations are Africa, Europe and
  Asia.

- [gypsum](/packages/gypsum) Client for the gypsum REST API
  (https://gypsum.artifactdb.com), a cloud-based file store in the
  ArtifactDB ecosystem. This package provides functions for uploads,
  downloads, and various adminstrative and management tasks. Check
  out the documentation at
  https://github.com/ArtifactDB/gypsum-worker for more details.

- [hdxmsqc](/packages/hdxmsqc) The hdxmsqc package enables us to
  analyse and visualise the quality of HDX-MS experiments. Either as
  a final quality check before downstream analysis and publication or
  as part of a interative procedure to determine the quality of the
  data. The package builds on the QFeatures and Spectra packages to
  integrate with other mass-spectrometry data.

- [HicAggR](/packages/HicAggR) This package provides a set of
  functions useful in the analysis of 3D genomic interactions. It
  includes the import of standard HiC data formats into R and HiC
  normalisation procedures. The main objective of this package is to
  improve the visualization and quantification of the analysis of HiC
  contacts through aggregation. The package allows to import 1D
  genomics data, such as peaks from ATACSeq, ChIPSeq, to create
  potential couples between features of interest under user-defined
  parameters such as distance between pairs of features of interest.
  It allows then the extraction of contact values from the HiC data
  for these couples and to perform Aggregated Peak Analysis (APA) for
  visualization, but also to compare normalized contact values
  between conditions. Overall the package allows to integrate 1D
  genomics data with 3D genomics data, providing an easy access to
  HiC contact values.

- [HybridExpress](/packages/HybridExpress) HybridExpress can be used
  to perform comparative transcriptomics analysis of hybrids (or
  allopolyploids) relative to their progenitor species. The package
  features functions to perform exploratory analyses of sample
  grouping, identify differentially expressed genes in hybrids
  relative to their progenitors, classify genes in expression
  categories (N = 12) and classes (N = 5), and perform functional
  analyses. We also provide users with graphical functions for the
  seamless creation of publication-ready figures that are commonly
  used in the literature.

- [igvShiny](/packages/igvShiny) This package is a wrapper of
  Integrative Genomics Viewer (IGV). It comprises an htmlwidget
  version of IGV. It can be used as a module in Shiny apps.

- [iSEEfier](/packages/iSEEfier) iSEEfier provides a set of
  functionality to quickly and intuitively create, inspect, and
  combine initial configuration objects. These can be conveniently
  passed in a straightforward manner to the function call to launch
  iSEE() with the specified configuration. This package currently
  works seamlessly with the sets of panels provided by the iSEE and
  iSEEu packages, but can be extended to accommodate the usage of any
  custom panel (e.g. from iSEEde, iSEEpathways, or any panel
  developed independently by the user).

- [knowYourCG](/packages/knowYourCG) knowYourCG automates the
  functional analysis of DNA methylation data. The package tests the
  enrichment of discrete CpG probes across thousands of curated
  biological and technical features. GSEA-like analysis can be
  performed on continuous methylation data query sets. knowYourCG can
  also take beta matrices as input to perform feature aggregation
  over the curated database sets.

- [limpca](/packages/limpca) This package has for objectives to
  provide a method to make Linear Models for high-dimensional
  designed data. limpca applies a GLM (General Linear Model) version
  of ASCA and APCA to analyse multivariate sample profiles generated
  by an experimental design. ASCA/APCA provide powerful visualization
  tools for multivariate structures in the space of each effect of
  the statistical model linked to the experimental design and
  contrarily to MANOVA, it can deal with mutlivariate datasets having
  more variables than observations. This method can handle unbalanced
  design.

- [lute](/packages/lute) Provides a framework for adjustment on cell
  type size when performing bulk transcripomics deconvolution. The
  main framework function provides a means of reference normalization
  using cell size scale factors. It allows for marker selection and
  deconvolution using non-negative least squares (NNLS) by default.
  The framework is extensible for other marker selection and
  deconvolution algorithms, and users may reuse the generics,
  methods, and classes for these when developing new algorithms.

- [MAPFX](/packages/MAPFX) MAPFX is an end-to-end toolbox that
  pre-processes the raw data from MPC experiments (e.g., BioLegend's
  LEGENDScreen and BD Lyoplates assays), and further imputes the
  ‘missing’ infinity markers in the wells without those measurements.
  The pipeline starts by performing background correction on raw
  intensities to remove the noise from electronic baseline
  restoration and fluorescence compensation by adapting a
  normal-exponential convolution model. Unwanted technical variation,
  from sources such as well effects, is then removed using a
  log-normal model with plate, column, and row factors, after which
  infinity markers are imputed using the informative backbone markers
  as predictors. The completed dataset can then be used for
  clustering and other statistical analyses. Additionally, MAPFX can
  be used to normalise data from FFC assays as well.

- [methodical](/packages/methodical) DNA methylation is generally
  considered to be associated with transcriptional silencing.
  However, comprehensive, genome-wide investigation of this
  relationship requires the evaluation of potentially millions of
  correlation values between the methylation of individual genomic
  loci and expression of associated transcripts in a relatively large
  numbers of samples. Methodical makes this process quick and easy
  while keeping a low memory footprint. It also provides a novel
  method for identifying regions where a number of methylation sites
  are consistently strongly associated with transcriptional
  expression. In addition, Methodical enables housing DNA methylation
  data from diverse sources (e.g. WGBS, RRBS and methylation arrays)
  with a common framework, lifting over DNA methylation data between
  different genome builds and creating base-resolution plots of the
  association between DNA methylation and transcriptional activity at
  transcriptional start sites.

- [methyLImp2](/packages/methyLImp2) This package allows to estimate
  missing values in DNA methylation data. methyLImp method is based
  on linear regression since methylation levels show a high degree of
  inter-sample correlation. Implementation is parallelised over
  chromosomes since probes on different chromosomes are usually
  independent. Mini-batch approach to reduce the runtime in case of
  large number of samples is available.

- [MGnifyR](/packages/MGnifyR) Utility package to facilitate
  integration and analysis of EBI MGnify data in R. The package can
  be used to import microbial data for instance into
  TreeSummarizedExperiment (TreeSE). In TreeSE format, the data is
  directly compatible with miaverse framework.

- [MIRit](/packages/MIRit) MIRit is an R package that provides
  several methods for investigating the relationships between miRNAs
  and genes in different biological conditions. In particular, MIRit
  allows to explore the functions of dysregulated miRNAs, and makes
  it possible to identify miRNA-gene regulatory axes that control
  biological pathways, thus enabling the users to unveil the
  complexity of miRNA biology. MIRit is an all-in-one framework that
  aims to help researchers in all the central aspects of an
  integrative miRNA-mRNA analyses, from differential expression
  analysis to network characterization.

- [mobileRNA](/packages/mobileRNA) Genomic analysis can be utilised
  to identify differences between RNA populations in two conditions,
  both in production and abundance. This includes the identification
  of RNAs produced by multiple genomes within a biological system.
  For example, RNA produced by pathogens within a host or mobile RNAs
  in plant graft systems. The mobileRNA package provides methods to
  pre-process, analyse and visualise the sRNA and mRNA populations
  based on the premise of mapping reads to all genotypes at the same
  time.

- [mosdef](/packages/mosdef) This package provides functionality to
  run a number of tasks in the differential expression analysis
  workflow. This encompasses the most widely used steps, from running
  various enrichment analysis tools with a unified interface to
  creating plots and beautifying table components linking to external
  websites and databases. This streamlines the generation of
  comprehensive analysis reports.

- [motifTestR](/packages/motifTestR) Taking a set of sequence motifs
  as PWMs, test a set of sequences for over-representation of these
  motifs, as well as any positional features within the set of
  motifs. Enrichment analysis can be undertaken using multiple
  statistical approaches. The package also contains core functions to
  prepare data for analysis, and to visualise results.

- [multistateQTL](/packages/multistateQTL) A collection of tools for
  doing various analyses of multi-state QTL data, with a focus on
  visualization and interpretation. The package 'multistateQTL'
  contains functions which can remove or impute missing data,
  identify significant associations, as well as categorise features
  into global, multi-state or unique. The analysis results are stored
  in a 'QTLExperiment' object, which is based on the
  'SummarisedExperiment' framework.

- [pathlinkR](/packages/pathlinkR) pathlinkR is an R package designed
  to facilitate analysis of RNA-Seq results. Specifically, our aim
  with pathlinkR was to provide a number of tools which take a list
  of DE genes and perform different analyses on them, aiding with the
  interpretation of results. Functions are included to perform
  pathway enrichment, with muliplte databases supported, and tools
  for visualizing these results. Genes can also be used to create and
  plot protein-protein interaction networks, all from inside of R.

- [Pedixplorer](/packages/Pedixplorer) Routines to handle family data
  with a Pedigree object. The initial purpose was to create
  correlation structures that describe family relationships such as
  kinship and identity-by-descent, which can be used to model family
  data in mixed effects models, such as in the coxme function. Also
  includes a tool for Pedigree drawing which is focused on producing
  compact layouts without intervention. Recent additions include
  utilities to trim the Pedigree object with various criteria, and
  kinship for the X chromosome.

- [pgxRpi](/packages/pgxRpi) The package is an R wrapper for
  Progenetix REST API built upon the Beacon v2 protocol. Its purpose
  is to provide a seamless way for retrieving genomic data from
  Progenetix database—an open resource dedicated to curated
  oncogenomic profiles. Empowered by this package, users can
  effortlessly access and visualize data from Progenetix.

- [PIPETS](/packages/PIPETS) PIPETS provides statistically robust
  analysis for 3'-seq/term-seq data. It utilizes a sliding window
  approach to apply a Poisson Distribution test to identify genomic
  positions with termination read coverage that is significantly
  higher than the surrounding signal. PIPETS then condenses proximal
  signal and produces strand specific results that contain all
  significant termination peaks.

- [PIUMA](/packages/PIUMA) The PIUMA package offers a tidy pipeline
  of Topological Data Analysis frameworks to identify and
  characterize communities in high and heterogeneous dimensional
  data.

- [PLSDAbatch](/packages/PLSDAbatch) A novel framework to correct for
  batch effects prior to any downstream analysis in microbiome data
  based on Projection to Latent Structures Discriminant Analysis. The
  main method is named “PLSDA-batch”. It first estimates treatment
  and batch variation with latent components, then subtracts
  batch-associated components from the data whilst preserving
  biological variation of interest. PLSDA-batch is highly suitable
  for microbiome data as it is non-parametric, multivariate and
  allows for ordination and data visualisation. Combined with
  centered log-ratio transformation for addressing uneven library
  sizes and compositional structure, PLSDA-batch addresses all
  characteristics of microbiome data that existing correction methods
  have ignored so far. Two other variants are proposed for 1/
  unbalanced batch x treatment designs that are commonly encountered
  in studies with small sample sizes, and for 2/ selection of
  discriminative variables amongst treatment groups to avoid
  overfitting in classification problems. These two variants have
  widened the scope of applicability of PLSDA-batch to different data
  settings.

- [pwalign](/packages/pwalign) The two main functions in the package
  are pairwiseAlignment() and stringDist(). The former solves
  (Needleman-Wunsch) global alignment, (Smith-Waterman) local
  alignment, and (ends-free) overlap alignment problems. The latter
  computes the Levenshtein edit distance or pairwise alignment score
  matrix for a set of strings.

- [rawDiag](/packages/rawDiag) Optimizing methods for liquid
  chromatography coupled to mass spectrometry (LC-MS) poses a
  nontrivial challenge. The rawDiag package facilitates rational
  method optimization by generating MS operator-tailored diagnostic
  plots of scan-level metadata. The package is designed for use on
  the R shell or as a Shiny application on the Orbitrap instrument
  PC.

- [rBLAST](/packages/rBLAST) Seamlessly interfaces the Basic Local
  Alignment Search Tool (BLAST) to search genetic sequence data
  bases. This work was partially supported by grant no. R21HG005912
  from the National Human Genome Research Institute.

- [saseR](/packages/saseR) saseR is a highly performant and fast
  framework for aberrant expression and splicing analyses. The main
  functions are: BamtoAspliCounts, convertASpli, calculateOffsets,
  saseRfindEncodingDim, saseRfitFor. information upon how to use these
  functions, check out our [vignette](https://bioconductor.org/packages/release/bioc/vignettes/saseR/inst/doc/saseR-vignette.html)
  and the saseR paper: Segers, A. et al. (2023). Juggling offsets
  unlocks RNA-seq tools for fast scalable differential usage,
  aberrant splicing and expression analyses. bioRxiv.
  https://doi.org/10.1101/2023.06.29.547014.

- [scMitoMut](/packages/scMitoMut) This package is designed for
  analyzing mitochondrial mutations using single-cell sequencing
  data, such as scRNASeq and scATACSeq (preferably the latter due to
  RNA editing issues). It includes functions for mutation filtering
  and visualization. In the future, the visualization tool will
  become an independent package. Mutation filtering is performed by
  fitting a statistical model to account for various sources of
  noise, including PCR error, sequencing error, mtDNA sampling and/or
  heteroplasmy dynamics. The model tests whether the observed allele
  frequency of a locus in a cell can be explained by the noise model.
  If not, we classify it as a mutation. The input for this analysis
  is the allele frequency. The noise model consists of three
  independent models: binomial, binomial-mixture, and beta-binomial
  models.

- [scMultiSim](/packages/scMultiSim) scMultiSim simulates paired
  single cell RNA-seq, single cell ATAC-seq and RNA velocity data,
  while incorporating mechanisms of gene regulatory networks,
  chromatin accessibility and cell-cell interactions. It allows users
  to tune various parameters controlling the amount of each
  biological factor, variation of gene-expression levels, the
  influence of chromatin accessibility on RNA sequence data, and so
  on. It can be used to benchmark various computational methods for
  single cell multi-omics data, and to assist in experimental design
  of wet-lab experiments.

- [shiny.gosling](/packages/shiny.gosling) A Grammar-based Toolkit
  for Scalable and Interactive Genomics Data Visualization.
  http://gosling-lang.org/. This R package is based on gosling.js. It
  uses R functions to create gosling plots that could be embedded
  onto R Shiny apps.

- [simPIC](/packages/simPIC) simPIC is a package for simulating
  single-cell ATAC-seq count data. It provides a user-friendly, well
  documented interface for data simulation. Functions are provided
  for parameter estimation, realistic scATAC-seq data simulation, and
  comparing real and simulated datasets.

- [SingleCellAlleleExperiment](/packages/SingleCellAlleleExperiment)
  Defines a S4 class that is based on SingleCellExperiment. In
  addition to the usual gene layer the object can also store data for
  immune genes such as HLAs, Igs and KIRs at allele and functional
  level. The package is part of a workflow named single-cell
  ImmunoGenomic Diversity (scIGD), that firstly incorporates
  allele-aware quantification data for immune genes. This new data
  can then be used with the here implemented data structure and
  functionalities for further data handling and data analysis.

- [sketchR](/packages/sketchR) Provides an R interface for various
  subsampling algorithms implemented in python packages. Currently,
  interfaces to the geosketch and scSampler python packages are
  implemented. In addition it also provides diagnostic plots to
  evaluate the subsampling.

- [smartid](/packages/smartid) This package enables automated
  selection of group specific signature, especially for rare
  population. The package is developed for generating specifc lists
  of signature genes based on Term Frequency-Inverse Document
  Frequency (TF-IDF) modified methods. It can also be used as a new
  gene-set scoring method or data transformation method. Multiple
  visualization functions are implemented in this package.

- [smoothclust](/packages/smoothclust) Method for segmentation of
  spatial domains and spatially-aware clustering in spatial
  transcriptomics data. The method generates spatial domains with
  smooth boundaries by smoothing gene expression profiles across
  neighboring spatial locations, followed by unsupervised clustering.
  Spatial domains consisting of consistent mixtures of cell types may
  then be further investigated by applying cell type compositional
  analyses or differential analyses.

- [SpaceMarkers](/packages/SpaceMarkers) Spatial transcriptomic
  technologies have helped to resolve the connection between gene
  expression and the 2D orientation of tissues relative to each
  other. However, the limited single-cell resolution makes it
  difficult to highlight the most important molecular interactions in
  these tissues. SpaceMarkers, R/Bioconductor software, can help to
  find molecular interactions, by identifying genes associated with
  latent space interactions in spatial transcriptomics.

- [spillR](/packages/spillR) Channel interference in mass cytometry
  can cause spillover and may result in miscounting of protein
  markers. We develop a nonparametric finite mixture model and use
  the mixture components to estimate the probability of spillover. We
  implement our method using expectation-maximization to fit the
  mixture model.

- [spoon](/packages/spoon) This package addresses the mean-variance
  relationship in spatially resolved transcriptomics data. Precision
  weights are generated for individual observations using Empirical
  Bayes techniques. These weights are used to rescale the data and
  covariates, which are then used as input in spatially variable gene
  detection tools.

- [SpotSweeper](/packages/SpotSweeper) Spatially-aware quality
  control (QC) software for both spot-level and artifact-level QC in
  spot-based spatial transcripomics, such as 10x Visium. These
  methods calculate local (nearest-neighbors) mean and variance of
  standard QC metrics (library size, unique genes, and mitochondrial
  percentage) to identify outliers spot and large technical
  artifacts. Scales linearly with the number of spots and is designed
  to be used with 'SpatialExperiment' objects.

- [SurfR](/packages/SurfR) Identify Surface Protein coding genes from
  a list of candidates. Systematically download data from GEO and
  TCGA or use your own data. Perform DGE on bulk RNAseq data. Perform
  Meta-analysis. Descriptive enrichment analysis and plots.

- [tidyCoverage](/packages/tidyCoverage) `tidyCoverage` framework
  enables tidy manipulation of collections of genomic tracks and
  features using `tidySummarizedExperiment` methods. It facilitates
  the extraction, aggregation and visualization of genomic coverage
  over individual or thousands of genomic loci, relying on
  `CoverageExperiment` and `AggregatedCoverage` classes. This
  accelerates the integration of genomic track data in genomic
  analysis workflows.

- [tidyomics](/packages/tidyomics) The tidyomics ecosystem is a set
  of packages for ’omic data analysis that work together in harmony;
  they share common data representations and API design, consistent
  with the tidyverse ecosystem. The tidyomics package is designed to
  make it easy to install and load core packages from the tidyomics
  ecosystem with a single command.

- [tidySpatialExperiment](/packages/tidySpatialExperiment)
  tidySpatialExperiment provides a bridge between the
  SpatialExperiment package and the tidyverse ecosystem. It creates
  an invisible layer that allows you to interact with a
  SpatialExperiment object as if it were a tibble; enabling the use
  of functions from dplyr, tidyr, ggplot2 and plotly. But,
  underneath, your data remains a SpatialExperiment object.

- [tpSVG](/packages/tpSVG) The goal of `tpSVG` is to detect and
  visualize spatial variation in the gene expression for spatially
  resolved transcriptomics data analysis. Specifically, `tpSVG`
  introduces a family of count-based models, with generalizable
  parametric assumptions such as Poisson distribution or negative
  binomial distribution. In addition, comparing to currently
  available count-based model for spatially resolved data analysis,
  the `tpSVG` models improves computational time, and hence greatly
  improves the applicability of count-based models in SRT data
  analysis.

- [transmogR](/packages/transmogR) transmogR provides the tools
  needed to crate a new reference genome or reference transcriptome,
  using a set of variants. Variants can be any combination of SNPs,
  Insertions and Deletions. The intended use-case is to enable
  creation of variant-modified reference transcriptomes for
  incorporation into transcriptomic pseudo-alignment workflows, such
  as salmon.

- [treeclimbR](/packages/treeclimbR) The arrangement of hypotheses in
  a hierarchical structure appears in many research fields and often
  indicates different resolutions at which data can be viewed. This
  raises the question of which resolution level the signal should
  best be interpreted on. treeclimbR provides a flexible method to
  select optimal resolution levels (potentially different levels in
  different parts of the tree), rather than cutting the tree at an
  arbitrary level. treeclimbR uses a tuning parameter to generate
  candidate resolutions and from these selects the optimal one.

- [txdbmaker](/packages/txdbmaker) A set of tools for making TxDb
  objects from genomic annotations from various sources (e.g. UCSC,
  Ensembl, and GFF files). These tools allow the user to download the
  genomic locations of transcripts, exons, and CDS, for a given
  assembly, and to import them in a TxDb object. TxDb objects are
  implemented in the GenomicFeatures package, together with flexible
  methods for extracting the desired features in convenient formats.

- [UCSC.utils](/packages/UCSC.utils) A set of low-level utilities to
  retrieve data from the UCSC Genome Browser. Most functions in the
  package access the data via the UCSC REST API but some of them
  query the UCSC MySQL server directly. Note that the primary purpose
  of the package is to support higher-level functionalities
  implemented in downstream packages like GenomeInfoDb or txdbmaker.

- [UPDhmm](/packages/UPDhmm) Uniparental disomy (UPD) is a genetic
  condition where an individual inherits both copies of a chromosome
  or part of it from one parent, rather than one copy from each
  parent. This package contains a HMM for detecting UPDs through HTS
  (High Throughput Sequencing) data from trio assays. By analyzing
  the genotypes in the trio, the model infers a hidden state (normal,
  father isodisomy, mother isodisomy, father heterodisomy and mother
  heterodisomy).

- [VisiumIO](/packages/VisiumIO) The package allows users to readily
  import spatial data obtained from either the 10X website or from
  the Space Ranger pipeline. Supported formats include tar.gz, h5,
  and mtx files. Multiple files can be imported at once with *List
  type of functions. The package represents data mainly as
  SpatialExperiment objects.

New Data Experiment Packages
=====================

There are 10 new data experiment packages in this release of Bioconductor.

- [curatedPCaData](/packages/curatedPCaData) The package
  curatedPCaData offers a selection of annotated prostate cancer
  datasets featuring multiple omics, manually curated metadata, and
  derived downstream variables. The studies are offered as
  MultiAssayExperiment (MAE) objects via ExperimentHub, and comprise
  of clinical characteristics tied to gene expression, copy number
  alteration and somatic mutation data. Further, downstream features
  computed from these multi-omics data are offered. Multiple
  vignettes help grasp characteristics of the various studies and
  provide example exploratory and meta-analysis of leveraging the
  multiple studies provided here-in.

- [CytoMethIC](/packages/CytoMethIC) This package provides DNA
  methylation-based prediction of cancer type, molecular signature
  and clinical outcomes. It provides convenience functions for
  missing value imputation, probe ID conversion, model interpretation
  and visualization. The package links to our models on
  ExperimentHub. The package currently supports HM450, EPIC and
  EPICv2.

- [homosapienDEE2CellScore](/packages/homosapienDEE2CellScore) This
  is a data package for normalised homosapien data downloaded from
  DEE2. The package both downloads, normalises, and filters the data,
  and provides a way to access the data from a canonical store
  without needing local processing. This package was built as a way
  to generate and store canonical test data for CellScore.

- [JohnsonKinaseData](/packages/JohnsonKinaseData) The packages
  provides position specific weight matrices (PWMs) for 303 human
  serine/threonine kinases originally published in Johnson et al. 2023. It
  includes gene annotation for each kinase PWM and PWM matching scores for a set
  of 85603 curated human phosphosites which can be used to map a PWM score to
  its percentile rank. The package also includes basic functionality to score
  user provided phosphosites.

- [MouseAgingData](/packages/MouseAgingData) The MouseAgingData
  package provides analysis-ready data resources from different
  studies focused on aging and rejuvenation in mice. The package
  currently provides two 10x Genomics single-cell RNA-seq datasets.
  The first study profiled the aging mouse brain measured across
  37,089 cells (Ximerakis et al., 2019). The second study
  investigated parabiosis by profiling a total of 105,329 cells
  (Ximerakis & Holton et al., 2023). The datasets are provided as
  SingleCellExperiment objects and provide raw UMI counts and cell
  metadata.

- [muleaData](/packages/muleaData) ExperimentHubData package for the
  'mulea' comprehensive overrepresentation and functional enrichment
  analyser R package. Here we provide ontologies (gene sets) in a
  data.frame for 27 different organisms, ranging from Escherichia
  coli to human, all acquired from publicly available data sources.
  Each ontology is provided with multiple gene and protein
  identifiers. Please see the NEWS file for a list of changes in each
  version.

- [scaeData](/packages/scaeData) Contains default datasets used by
  the Bioconductor package SingleCellAlleleExperiment. The raw FASTQ
  files were sourced from publicly accessible datasets provided by
  10x Genomics. Subsequently, our scIGD snakemake workflow was
  employed to process these FASTQ files. The resulting output from
  scIGD constitutes to the contents of this data package.

- [SubcellularSpatialData](/packages/SubcellularSpatialData) This is
  a data package that hosts annotated sub-cellular localised datasets
  from the STOmics, Xenium and CosMx platforms. Specifically, it
  hosts datasets analysed in the publication Bhuva et. al, 2024
  titled "Library size confounds biology in spatial transcriptomics
  data". Raw transcript detections are hosted and functions to
  convert them to SpatialExperiment objects have been implemented.

- [TENxXeniumData](/packages/TENxXeniumData) Collection of Xenium
  spatial transcriptomics datasets provided by 10x Genomics,
  formatted into the Bioconductor classes, the SpatialExperiment or
  SpatialFeatureExperiment (SFE), to facilitate seamless integration
  into various applications, including examples, demonstrations, and
  tutorials. The constructed data objects include gene expression
  profiles, per-transcript location data, centroid, segmentation
  boundaries (e.g., cell or nucleus boundaries), and image.

- [TransOmicsData](/packages/TransOmicsData) Contains a collection of
  trans-omics datasets generated using various sequencing
  technologies such as RNA-seq, Mass spectrometry and ChIP-seq.
  Modalities include the bulk profiling of the phosphoproteome,
  proteome, transcriptome and epigenome. Data reflects the
  timecourses of different developmental systems from the mouse or
  human.


New Annotation Packages
=====================

There are 5 new annotation packages.

- [ath1121501frmavecs](/packages/ath1121501frmavecs) Annotation package for the
  implementation of the frozen Robust Multiarray Analysis procedure for
  Arabidopsis thaliana. This package was generated on the basis of frmaTools
  version 1.52.0.

- [EPICv2manifest](/packages/EPICv2manifest) A data.frame containing
  an extended probe manifest for the Illumina Infinium Methylation
  v2.0 Kit. Contains the complete manifest from the Illumina-provided
  EPIC-8v2-0_EA.csv, plus additional probewise information described
  in Peters et al. (2024).

- [IlluminaHumanMethylationEPICv2anno.20a1.hg38](/packages/IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  An annotation package for Illumina's EPIC v2.0 methylation arrays. The version
  2 covers more than 935K CpG sites in the human genome hg38. It is an update of
  the original EPIC v1.0 array (i.e., the 850K methylation array).
  
- [IlluminaHumanMethylationEPICv2manifest](/packages/IlluminaHumanMethylationEPICv2manifest)
  A manifest package for Illumina's EPIC v2.0 methylation arrays. The version 2
  covers more than 935K CpG sites in the human genome hg38. It is an update of
  the original EPIC v1.0 array (i.e., the 850K methylation array).
  
- [MafH5.gnomAD.v4.0.GRCh38](/packages/MafH5.gnomAD.v4.0.GRCh38) Store minor
  allele frequency data from the Genome Aggregation Database (gnomAD version
  4.0) for the human genome version GRCh38.

New Workflow Packages
=====================

There are no new workflow packages in this release of Bioconductor.

New Online Books
=====================

There is one new online book.

- [OHCA](https://bioconductor.org/books/3.19/OHCA/) The primary aim of this book is to introduce
  the R user to Hi-C analysis. This book starts with key concepts
  important for the analysis of chromatin conformation capture and
  then presents Bioconductor tools that can be leveraged to process,
  analyze, explore and visualize Hi-C data.

NEWS from existing Software Packages
===================================


[a4](/packages/a4)
--

                       Changes in version 1.51.1                        

- add Suggests: hgu95av2.db

[ADaCGH2](/packages/ADaCGH2)
-------

                 Changes in version 2.43.1 (2024-02-06)                 

- Removed dependency from snapCGH (and functionality related to
  snapCGH).
  The snapCGH package was deprecated in Bioconductor 3.18 and removed
  in 3.19.

[adverSCarial](/packages/adverSCarial)
------------

                 Changes in version 1.1.2 (2024-02-04)                  

- Add the "decile+x" modifications

[affxparser](/packages/affxparser)
----------

                 Changes in version 1.75.2 (2024-02-06)                 

Bug Fixes

- Fixed one sprintf-related coercion issue, reported by the GCC
compiler, that would produce an incorrect warning message.

- Fixed two sprintf-related flag issues, reported by the GCC
compiler,
in internal assertions, which would never be reached.

                 Changes in version 1.75.1 (2024-02-05)                 

Documentation

- Fix incorrectly documented arguments in a few functions.

                 Changes in version 1.75.0 (2023-10-24)                 

Notes

- The version number was bumped for the Bioconductor develop version,
which is now Bioconductor 3.19 for R (>= 4.4.0).

[AHMassBank](/packages/AHMassBank)
----------

                         Changes in version 1.3                         

AHMassBank 1.3.1

- Add MassBank releases 2023.06, 2023.09 and 2023.11.

[alevinQC](/packages/alevinQC)
--------

                       Changes in version 1.19.3                        

- Minor fixes to adapt to changes in ggplot2

                       Changes in version 1.19.2                        

- Add 'Close app' button to close the shiny application and return
data

[AlphaMissenseR](/packages/AlphaMissenseR)
--------------

                        Changes in version 1.0.0                        

- (v. 0.99.21) Update Zenodo data source to record 10813168, with
more
accessible 'CC-BY-4.0' license.
- (v. 0.99.16) Housekeeping
- Order vignettes introduction, visualization, issues.
- Use rjsoncons (>= 1.0.1), so that jsonlite is not a hard
dependency.
- Acknowledge additional funding sources; add ImmunoOncology
biocViews term.
- (v. 0.99.15) Respond to Bioconductor reviewer comments
- Use BiocBaseUtils for input assertions.
- Include range join SQL directly rather than via a non-exported
function.
- Report file size when downloading.
- Improve test coverage.
- See GitHub issue comment; thanks @LiNk-NY
- (v. 0.99.11) Use an S4 class for alphamissense_connection,
extending
duckdb_connection.
- (v. 0.99.10) Introduce visualization on AlphaFold predictions.
- (v. 0.99.9) Update maintainer email address, add 'Resource
unavailable' section to 'Issues & Solutions' vignette.
- (v. 0.99.7) Update to revised Zenodo API.
- (v. 0.99.5) BREAKING change. Requires duckdb >= 0.9.1, which cannot
read DuckDB databases created with earlier versions. Also introduces
changes to local cache naming. See the 'Issues & Solutions'
vignette.
- (v. 0.99.4) Change back to Bioconductor 'Software' package, to
allow
for regular testing of data access code.
- (v. 0.99.1) Change to Bioconductor 'Annotation' package.
- (v. 0.99.1) Update to change in Zenodo API.

                       Changes in version 0.99.0                        

- (v. 0.0.18) Rename (including updating existing tables) '#CHROM' to
'CHROM' in hg19 / hg38 tables
- (v. 0.0.17) Initial Bioconductor submission.

[AnnotationHub](/packages/AnnotationHub)
-------------

                       Changes in version 3.11.0                        

BUG CORRECTION

- (3.11.4) Correct subsetting. Didn't account for multiple files
  associated
  with single entry (ie. bam/bai).

[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.16.0                        

USER VISIBLE CHANGES

- (v 1.15.10) Validate API versions against hardcoded variables;
produce warning when discordant (@LiNk-NY, #101).

- (v 1.15.8) Add gcloud_storage() and gcloud_storage_buckets() to
create and manage Google Cloud Storage buckets (@LiNk-NY, #72).

- Gen3 services, avworkflow*_configuration() functions, install(),
repository(), and repositories() are defunct.

- (v 1.15.5) Catch avtable_import_status() errors in the response
object.

- (v 1.15.1) Update vignette with examples for avworkflow_info()
(@mtmorgan, @yubocheng).

BUG FIXES AND MINOR IMPROVEMENTS

- (v 1.15.11) Update Dockstore API file, version, and URL

- (v 1.15.9) Use assertions from BiocBaseUtils

- (v 1.15.7) Use URLencode for table in avtable and direct request to
Rawls endpoint (@LiNk-NY, #98)

- (v 1.15.6) Update the Dockstore API reference URL and use
api_referenc_url instead of API file (@LiNk-NY).

- Update namespace in vignette and examples (@kozo2, #54)

[AnVILWorkflow](/packages/AnVILWorkflow)
-------------

                        Changes in version 1.2.6                        

- New function `getWorkflowConfig` is added and the output from this
  function is now required as an input for the `config` argument of
  `currentInput`, `udpateInput`, and `runWorkflow` functions.
  This change makes the execution of different functions faster because
  they
  don't need to access AnVIL everytime to get the configuration.

- New function `AnVILBrowse` is added

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.27.4                        

- Fix the refresh issue for IGVSnapshot.R

                       Changes in version 1.27.3                        

- Fix the issue if no single end reads in the PE bam files for
shiftGAlignmentsList.R

                       Changes in version 1.27.2                        

- Try to decrease the memory cost for shiftGAlignmentsList.R

                       Changes in version 1.27.1                        

- Fix the issue caused by readGAlignmentsList does not support read
by
chunk

[ATACseqTFEA](/packages/ATACseqTFEA)
-----------

                        Changes in version 1.5.1                        

- add function importFimoBindingSites.

[atena](/packages/atena)
-----

                 Changes in version 1.10.0 (2024-04-27)                 

USER VISIBLE CHANGES

- Annotation metadata is now propagated to parameter and resulting
  SummarizedExperiment objects, which can be retrieved using
  'mcols(features(x))' and 'mcols(rowRanges(x))' methods, respectively.

- The TE subclass annotation getter functions 'getLTRs()',
  'getLINEs()', 'getSINEs()' and 'getDNAtransposons()' work now with
  both, parameter and resulting SummarizedExperiment objects.

- Manual pages and the vignette have been updated to illustrate the
  previous changes.

BUG FIXES

- Several bug fixes, mostly in how annotation metadata was stored and
  propagated.

[AUCell](/packages/AUCell)
------

                        Changes in version 1.25                         

- Remove shiny app to modify thresholds as rbokeh is deprecated.

[bacon](/packages/bacon)
-----

                       Changes in version 1.32.0                        

NEW FEATURES

- Added Random number generator arguments to `bacon`-function,
  `globalSeed` and `parallelSeed` that are passed to `BiocParallel`
  such that proper confidence intervals can be calculated for the
  estimated variables (see the example in the vignette). Thanks to
  Kendra Ferrier who made this contribution.

[BANDITS](/packages/BANDITS)
-------

                       Changes in version 1.18.1                        

- Updated from C++11 to C++17

[Banksy](/packages/Banksy)
------

                       Changes in version 0.99.8                        

- Add feature scaling options for PCA and UMAP

                       Changes in version 0.99.7                        

- SCE/SPE-compatible


[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.15.4 (2024-03-25)                 

- Minor bug fixe in the `subset` method for `BASiCS_Chain` objects

                 Changes in version 2.15.2 (2023-12-22)                 

- Minor bug fixes

[BatchQC](/packages/BatchQC)
-------

                       Changes in version 1.99.00                       

Changes made to dendrogram

- Added dendrogram_color_palette.R for coloring dendrogram
- Updated dendrogram.R allowing batch & category to plot together

[benchdamic](/packages/benchdamic)
----------

                 Changes in version 1.9.4 (2024-03-07)                  

- Re-added corncob (due to its previous removal from CRAN)

- Bug-fix: correct DA_Seurat() due to changes of Seurat::FindMarkers()

- Computed DA direction for ANCOM method

                 Changes in version 1.9.1 (2023-11-04)                  

- Temporarily removed corncob (due to its removal from CRAN)

                 Changes in version 1.9.0 (2023-10-24)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_18 branch

                 Changes in version 1.8.1 (2023-11-04)                  

- Porting the changes of devel version 1.9.1 to release

[BERT](/packages/BERT)
----

                       Changes in version 0.99.0                        

- BERT has been released.

- BERT provides an hierarchical approach to batch effect
  adjustment with tolerance to missing values. BERT has been
  found to outperform similar algorithmms, e.g. HarmonizR with
  regard to execution time, preservation of numerical values in
  the input data and flexibility wrt. covariables and references.

- The authors can be contacted via the BioConductor forum, GitHub
  issues and via email.

[betaHMM](/packages/betaHMM)
-------

                       Changes in version 0.99.0                        

betaHMM VERSION 0.1.x - 0.99.0

- Initial development version of betaHMM.

[bettr](/packages/bettr)
-----

                       Changes in version 0.99.3                        

- Add button to close app and return data

                       Changes in version 0.99.2                        

- Update installation instructions
- Expand vignette

                       Changes in version 0.99.1                        

- Add R (>= 4.4.0) to Depends

                       Changes in version 0.99.0                        

- Ready for submission to Bioconductor
- Add assembleSE() function to create a summary SE


[bioCancer](/packages/bioCancer)
---------

                       Changes in version 1.30.07                       

- fix issues in Networking Tab

                       Changes in version 1.30.06                       

- Merge mutation and Genetic Ptofile tabs

                       Changes in version 1.30.05                       

- rewrite the function getListProfData

- Omit function getMegaProfData user for genlist upper that 500.

                       Changes in version 1.30.04                       

- update ReactomeFI2021.RDS

- update DisGeNet0223.RDS

                       Changes in version 1.30.03                       

- update documentatio with new webapi

- progress in migration

                       Changes in version 1.30.02                       

- Start migration of cBioportal tab with the new cBioPortal webAPI

                       Changes in version 1.30.01                       

- Update api link

[BiocBaseUtils](/packages/BiocBaseUtils)
-------------

                        Changes in version 1.6.0                        

New features

- checkInstalled allows to check if a package is installed and
produces copy-and-paste text for installation.

[BiocCheck](/packages/BiocCheck)
---------

                       Changes in version 1.40.0                        

NEW FEATURES

- Check for duplicate chunk labels in vignettes and produce an error
  when
  present (@hpages, #199)

- Warn when Rnw vignettes are present (@jwokaty, #190)

- Add error when maintainer email not parse-able in DESCRIPTION

- Check for Bioconductor dependencies in the Imports and Depends
  fields;
  if none, produce a warning (@vjcitn).

- Remove redundant `DESCRIPTION` and `NAMESPACE` consistency checks;
  these
  are already included in `R CMD check`

BUG FIXES AND MINOR IMPROVEMENTS

- Produce `NOTE` instead of `WARNING` for infrastructure packages
  without
  Bioconductor dependencies

- `% \VignetteEngine` tags with spaces were not recognized; this has
  been
  fixed

- Skip `\value` checks for package and class documentation; `\formats`
  checked for data documentation

- Delegate system calls to `devtools::build` to build packages

- Update `checkRbuildignore` to flag `test` and `longtests` entries
  appropriately (@hpages, #197)

- Fix issue where BiocCheck would fail when errors caused by duplicate
  chunk
  labels were present (@hpages, #199)

- Avoid format checks e.g., length, indentation, and tabs for Rd files
  that use roxygen2

- Use the `tools:::analyze_licenses` function to check for restrictive
  license use (@hpages, #119)

- Use httr2 over httr for http requests

- Clarify `NOTE` and `WARNING` messages for `CITATION` file inclusion
  (#209).

[BiocFileCache](/packages/BiocFileCache)
-------------

                        Changes in version 2.11                         

BUG FIX

- (2.11.1) Merge PR to fix dbplyr compatibility issue

[BiocHubsShiny](/packages/BiocHubsShiny)
-------------

                        Changes in version 1.4.0                        

- No significant changes.

[BiocIO](/packages/BiocIO)
------

                       Changes in version 1.14.0                        

- No significant changes.

[BioCor](/packages/BioCor)
------

                        Changes in version 1.28                         

- Update tests to latest testthat change.

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.22.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- `biocLastBuildDate` has been removed from the package.

- Added a `recursive = FALSE` argument to `pkgBiocDeps`.

- Add `biocBuildStatusDB` function to download and parse build status
  info
  from BBS artifact.

- Include `hasBiocMaint` to check if a package is maintained by a
  particular
  email address (default to Bioconductor maintainer email).

BUG FIXES AND MINOR IMPROVEMENTS

- Use all `pkgType` options in `pkgBiocDeps` and `pkgBiocRevDeps`.

- Add example to `buildPkgDependencyGraph` using packages from
  `biocMaintained`.

- Include an example using `networkD3` in `buildPkgDependencyGraph`.

[biocroxytest](/packages/biocroxytest)
------------

                       Changes in version 0.99.2                        

Bug fixes

- Remove .Rproj file

                       Changes in version 0.99.1                        

Bug fixes

- Add files to .gitignore

                       Changes in version 0.99.0                        

New features

- Initial Bioconductor submission

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.60.0                        

USER VISIBLE CHANGES

- listEnsemblGenomes() and useEnsemblGenomes() now have a host
  argument,
  allowing you to select an Ensembl Genomes archive site. (Thanks to
  Hervé Pagès
  @hpages for the suggestion:
  https://github.com/grimbough/biomaRt/issues/93)

INTERNAL CHANGES

- Removed dependency on XML package and switched all functionality to
  xml2

- Swiched from httr to httr2 package for submitting queries to BioMart
  servers.

[Biostrings](/packages/Biostrings)
----------

                       Changes in version 2.72.0                        

NEW FEATURES

- Increase size of IO buffer from 20001 to 200000 in read_fasta_files.c
  and read_fastq_files.c. With this change readDNAStringSet() and
  family
  support FASTA/FASTQ files with lines up to 200000 characters.
  See https://github.com/Bioconductor/Biostrings/issues/59

- Add get_XStringSet_width() to Biostrings C interface.

SIGNIFICANT USER-VISIBLE CHANGES

- pairwiseAlignment() and related have moved to the new pwalign
  package.
  List of functions that are now implemented in pwalign:
  - pairwiseAlignment
  - PairwiseAlignments
  - PairwiseAlignmentsSingleSubject
  - writePairwiseAlignments
  - alignedPattern
  - alignedSubject
  - insertion
  - deletion
  - unaligned
  - aligned
  - indel
  - nindel
  - nedit
  - pid
  - mismatchTable
  - mismatchSummary
  - compareStrings
  - stringDist
  - nucleotideSubstitutionMatrix
  - errorSubstitutionMatrices
  - qualitySubstitutionMatrices
  Note that they are still temporarily defined in Biostrings but now
  they just call the corresponding function in pwalign. Since this is a
  temporary redirect, the user also gets a warning that tells them to
  use
  the fully qualified name (e.g. pwalign::pairwiseAlignment()) to call
  the function.

- The BLOSUM and PAM scoring matrices have also moved to the new
  pwalign
  package. List of scoring matrices that are now located in pwalign:
  BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM100, PAM30, PAM40,
  PAM70,
  PAM120, PAM250.

DEPRECATED AND DEFUNCT

- Deprecate matchprobes() and longestConsecutive().

- needwunsQS() is now defunct (after being deprecated for 15+ years).

BUG FIXES

- Detect buffer overflow in writeXStringSet() and raise error instead
  of crash. See https://github.com/Bioconductor/Biostrings/issues/20

- Make sure read*StringSet() closes all input file handles when done.
  See https://support.bioconductor.org/p/9157031/

[bnbc](/packages/bnbc)
----

                        Changes in version 1.25                         

- Fixed an issue with build error because of missing Suggests

- Added support for cooler files using HiCBricks (this happened
  in an earlier verion).

- Various documentation fixes.

[BREW3R.r](/packages/BREW3R.r)
--------

                       Changes in version 0.99.2                        

- Fix some typos and lint code.

                       Changes in version 0.99.1                        

- Internally, a lot of rewritting to remove dplyr usage.

- Use rlang to control verbosity

- Fix cases when 2 genes end at the same base before extension

                       Changes in version 0.99.0                        

- Initial Bioconductor submission.

[BSgenome](/packages/BSgenome)
--------

                       Changes in version 1.72.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette documenting the seed-based forging process has moved
  to the BSgenomeForge package. This vignette used to be called
  the "BSgenomeForge vignette" but it has now be renamed
  the "Advanced BSgenomeForge usage" vignette.

- All the tools related to the seed-based forging process and used by
  the "Advanced BSgenomeForge usage" vignette are now in the
  BSgenomeForge
  package.

[BSgenomeForge](/packages/BSgenomeForge)
-------------

                        Changes in version 1.4.0                        

NEW FEATURES

- Add forgeBSgenomeDataPkgFromTwobitFile()

- fastaTo2bit() improvement: the function now automatically replaces
  IUPAC
  ambiguity letters not supported by the 2bit format with N's and issue
  a warning.

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette documenting the seed-based forging process that used to
  be in the BSgenome software package is now in the BSgenomeForge
  package.
  This vignette used to be called the "BSgenomeForge vignette" but it
  has
  now be renamed the "Advanced BSgenomeForge usage" vignette.

- All the tools related to the seed-based forging process and used by
  the new "Advanced BSgenomeForge usage" vignette are now in the
  BSgenomeForge package.

[bugsigdbr](/packages/bugsigdbr)
---------

                        Changes in version 1.8.2                        

- Uniform first-letter capitalization for Body site and Condition
terms

[CAMERA](/packages/CAMERA)
------

                       Changes in version 1.59.1                        

BUG FIXES

- Add multtest to dependencies to fix testing on daily build

[canceR](/packages/canceR)
------

                       Changes in version 1.36.5                        

- fix issues with indexing Id of Cases and Genetic Profiles

- Update functionnnalities of GSEA

- add an example of GCT and CSL fiiles to use in GSEA

- fix phenoTest analysis

                       Changes in version 1.36.4                        

- modify getCases(), getGenProfs()

- fix get Clinical data

- fix get Mutation data

- fix get specific mutation

- fix get methylation data

- fix get Profile Data for a single gene

                       Changes in version 1.36.3                        

- update welcome panel. get all studies

                       Changes in version 1.36.2                        

- Start migration to the new WebAPI of cBioPortal

[Cardinal](/packages/Cardinal)
--------

                        Changes in version 3.5.6                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add optional 'run' argument to 'slice()'

- Add optional 'plotly' support for all plotting methods

                        Changes in version 3.5.5                        

BUG FIXES

- Fix bug in 'selectROI()' region plot updates

                        Changes in version 3.5.4                        

NEW FEATURES

- Re-add support for 3D coordinates and 'image3D()'

SIGNIFICANT USER-VISIBLE CHANGES

- Export 'vizi_style()' and 'vizi_engine()' functions

- Add 'check' argument to 'readMSIData()' for uuid+checksum

BUG FIXES

- Fix bug in plotting PLS/OPLS coefficients

- Fix removing position columns with 'coord()<-'

                        Changes in version 3.5.3                        

NEW FEATURES

- Add support for multiple instance learning for 'PLS()'
  and 'spatialShrunkenCentroids()' with 'bags' argument

- Add 'estimateReferenceMz()' for more easily estimating
  profile m/z-values from 'MSImagingArrays'

BUG FIXES

- Make sure 'bin()' still respects 'resolution' when 'ref' is
  specified (in which case the range is taken from 'ref')

- Fix 'crossValidate()' for multiple instance learning

                        Changes in version 3.5.2                        

BUG FIXES

- Fix passing 'tolerance' in 'peakProcess()' when not needed

- Fix bugs in 'selectROI()' not updating the plot

- Fix bugs in 'selectROI()' selecting wrong pixels

- Make sure 'peakAlign()' produces similar results when
  re-aligning peaks with the same tolerance

                        Changes in version 3.5.1                        

NEW FEATURES

- Ground-up rewrite to take advantage of matter v2.5 features

- New class 'SpectraArrays' for arrays of spectra

- New class 'SpectralImagingData' for spectral imaging data

- New class 'SpectralImagingArrays' for raw spectra

- New class 'SpectralImagingExperiment' for centroided spectra

- New class 'MSImagingArrays' for raw mass spectra

- Updated class 'MSImagingExperiment'

- Updated class 'XDataFrame'

- Updated classes 'PositionDataFrame' and 'MassDataFrame'

- New apply method 'spectrapply()'

- New processing method 'recalibrate()'

- New processing method 'bin()'

- New processing method 'peakProcess()'

- New classes 'SpatialResults' and 'ResultsList'

- Updated spatial methods 'findNeighbors()' and 'spatialWeights()'

- Updated spatial methods 'colocalized()'

- Updated projection methods 'PCA()' and 'spatialFastmap()'

- New projection method 'NMF()'

- Updated stats methods 'PLS()' and 'OPLS()'

- Updated stats method 'spatialKMeans()'

- Updated stats method 'spatialShrunkenCentroids()'

- Updated stats method 'spatialDGMM()'

- Improved visualization methods 'plot()' and 'image()'

SIGNIFICANT USER-VISIBLE CHANGES

- Deprecated 'smoothSignal()' -- use 'smooth()'

- Deprecated 'mzBin()' -- use 'bin()'

- Deprecated 'mzAlign()' -- use 'recalibrate()'

- Deprecated 'mzFilter()' -- use 'subsetFeatures()'

- Deprecated 'peakFilter()' -- use 'subsetFeatures()'

- Deprecated 'aggregate()' -- use 'summarizeFeatures()'

- Deprecated 'featureApply()' -- use 'summarizeFeatures()'

- Deprecated 'pixelApply()' -- use 'summarizePixels()'

                 Changes in version 3.4.3 (2023-11-22)                  

BUG FIXES

- Fixed default calculation of reference peaks in 'peakAlign()'

                 Changes in version 3.4.2 (2023-11-16)                  

BUG FIXES

- Fixed error in 'readImzML()' if spectrum representation is missing

- Fixed bug in 'writeImzML()' causing overlapping offsets

                 Changes in version 3.4.1 (2023-10-25)                  

BUG FIXES

- Fixed bug in 'writeImzML()' causing malformed cvParam tag

[CardinalIO](/packages/CardinalIO)
----------

                        Changes in version 1.1.6                        

BUG FIXES

- Non-existent 'extraArrays' produce NULL instead of error

                        Changes in version 1.1.5                        

BUG FIXES

- Allow root tags in 'ImzMeta' (e.g., if a tag is unknown)

                        Changes in version 1.1.4                        

SIGNIFICANT USER-VISIBLE CHANGES

- Updated README.md

                        Changes in version 1.1.3                        

NEW FEATURES

- Made 'writeAnalyze()' an S4 generic

- Add optional 'positions' parameter to 'writeAnalyze()'

SIGNIFICANT USER-VISIBLE CHANGES

- Change 'writeAnalyze()' parameter 'mz' to 'domain'

                        Changes in version 1.1.2                        

BUG FIXES

- Remove case sensitivity when checking file extensions

                        Changes in version 1.1.1                        

BUG FIXES

- Fix error when setting 'ImzMeta' elements to NULL

[cBioPortalData](/packages/cBioPortalData)
--------------

                       Changes in version 2.16.0                        

Bug fixes and minor improvements

- Improvements to links in documentation and citation section in the
vignette

[celda](/packages/celda)
-----

                 Changes in version 1.18.2 (2024-04-02)                 

- Updated Makevar files to new CRAN standards
- Fixed unit test causing error

                 Changes in version 1.18.1 (2023-11-05)                 

- Update to match Bioconductor release version
- Removed multipanelfigure as a dependency

[CellBarcode](/packages/CellBarcode)
-----------

                        Changes in version 1.9.1                        

NEW FEATURES

- Adding function bc_create_BarcodeObj enabling import processed barcode data
- a16a9e7 feat: Adding function bc_create_BarcodeObj enabling import process barcode data
- 45d886a test: adding test for processed data import
- 1d5f4de doc: add bc_create_BarcodeObj document

FIX

- Make sure object is data.frame instead of tibble or data.table
- 895e9a9 fix: make sure the data.frame is pure data.frame instead of tibble or data.table
- 31e6989 fix: error while checking data.frame

CHORE

- Error messages and input checks
- f84b2a1 chore: message if no barcodes found in bulk sequencings
- e0cded9 chore: Error message when no barcodes are found
- 3beef87 chore: Check input for sam and bam file

[cellxgenedp](/packages/cellxgenedp)
-----------

                         Changes in version 1.8                         

SIGNIFICANT USER-VISIBLE CHANGES

- (v 1.7.2) Add vignette 'Case studies', include identifying dataset
authors and using ontolgies to identify datasets

- (v 1.7.1) Update vignette section on dataset visualization to

BUG FIXES

- (v 1.7.3) Use all collections, datasets, files to determine
available columns

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.37.4                        

- Handle the biofilecache error.

                       Changes in version 3.37.3                        

- use an updated and unified single vignette to replace the four old
  ones

- remove oligoNucleotideEnrichment.R and incorporate it into
  summarizePatternsInPeaks.R

- update findMotifsInPromoterSeqs.R to match changes in
  summarizePatternsInPeaks.R

                       Changes in version 3.37.2                        

- export function oligoNucleotideEnrichment

                       Changes in version 3.37.1                        

- export disable.logging parameter for venn.diagram called by
  makeVennDiagram

[CircSeqAlignTk](/packages/CircSeqAlignTk)
--------------

                 Changes in version 1.5.2 (2024-04-23)                  

SIGNIFICANT USER-VISIBLE CHANGES

- Slot name `reversed` of `CircSeqAlignTkCoverage` class was changed to
  `reverse`.

NEW FEATURES

- Add GUI mode. See build_app() for details.

[ClassifyR](/packages/ClassifyR)
---------

                        Changes in version 3.8.0                        

- 
  Extraction of seed for random number generator fixed.

[cleaver](/packages/cleaver)
-------

                 Changes in version 1.41.3 (2024-04-24)                 

- Add 3 different types of trypsin (as defined in PeptideMass):
  trypsin-high
  (identical to the current trypsin rule and to PeptideMass' "higher
  specificity"), trypsin-low (cleavage after K/R but not if followed by
  P),
  trypsin-simple (cleavage after K/R even if followed by P);
  The trypsin rule is identical to trypsin-high and will be kept for
  backward compatibility.

                 Changes in version 1.41.2 (2023-12-03)                 

- Add missing documentation for LysargiNase.

                 Changes in version 1.41.1 (2023-11-30)                 

- Add LysargiNase cleavage rule.

[ClustAll](/packages/ClustAll)
--------

                       Changes in version 0.99.0                        

- Initial release of the ClustAll package for multiple robust
  stratifications.

[clusterExperiment](/packages/clusterExperiment)
-----------------

                       Changes in version 2.23.1                        

- Removed option "MB" for argument `mergeMethod` for the function
  `mergeClusters`. This was a method of Meinshausen and Buhlmann
  (2005) to estimate the number of non-null genes, but their supporting
  package `howmany` that implemented the method is no longer supported
  on CRAN.

[clusterProfiler](/packages/clusterProfiler)
---------------

                       Changes in version 4.10.1                        

- bug fixed in parsing KEGG category (2024-03-07, Thu, #664)
- update citation (#656) and wikipedia data URL (2024-01-10, Wed,
#633)

[clustifyr](/packages/clustifyr)
---------

                 Changes in version 1.15.2 (2024-04-03)                 

- Add support for `Seurat` version 5 objects

                 Changes in version 1.15.1 (2023-10-31)                 

- Replace `Seurat` dependency with `SeuratObject`

[ClustIRR](/packages/ClustIRR)
--------

                       Changes in version 1.1.10                        

- New global clustering solution

[CompoundDb](/packages/CompoundDb)
----------

                         Changes in version 1.7                         

Changes in version 1.7.2

- Import method generics from ProtGenerics.

Changes in version 1.7.1

- Adapt script to create CompDb from MassBank to new MassBank
database
format.

[consensusSeekeR](/packages/consensusSeekeR)
---------------

                       Changes in version 1.31.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- This package now depends on R (>= 3.5.0) because serialized objects
  in
  serialize/load version 3 cannot be read in older versions of R.

                       Changes in version 1.31.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- The description field of the package is more detailed.

[COTAN](/packages/COTAN)
-----

                        Changes in version 2.3.6                        

- Refactored DEAOnCluster() to make it run faster.
- Now clustering functions dump the GDI check results for all clusters
- Changed default GDI threshold to 1.43
- Added new input to mergeUniformCellsClusters() to allow proper resume of interrupted merges
- Added possibility to query whether the COEX matrix is available in a COTAN object

                        Changes in version 2.3.5                        

- Made checks more strict when adding a clusterization or condition
- Increased reliability of clustering functions by improved error handling and by allowing retry runs on estimators functions

                        Changes in version 2.3.4                        

- Speed-up of GDI calculation via Rfast package
- Added possibility of using distance between clusters based on Zero-One matrix instead of DEA
- Added average floor to logFoldChangeOnClusters() to dampen extreme results when genes are essentially absent from a cluster.

                        Changes in version 2.3.3                        

- Added method to handle expression levels' change via log-normalized data: logFoldChangeOnClusters()
- Minor fix in the import of operators to align to new version of roxigen2
- Restored default adjustment method of pValueFromDEA() to "none" for backward compatibility reasons

                        Changes in version 2.3.2                        

- Solved issue with cleanPlots() when the number of cells exceeded 65536
- Added methods to calculate the COEX matrix only on a subset of the columns
- Now the function pValueFromDEA() returns the p-value adjusted for multi-test

                        Changes in version 2.3.1                        

- Stopped using explicit PCA via irlba package: using BioConductor PCAtools::pca instead

                        Changes in version 2.3.0                        

- First release in Bioconductor 3.19

[CRISPRball](/packages/CRISPRball)
----------

                 Changes in version 0.99.0 (2023-09-24)                 

- Submitted to Bioconductor.

[crisprBase](/packages/crisprBase)
----------

                        Changes in version 1.7.1                        

- Updated weights in the SaCas9 nuclease.

[crisprDesign](/packages/crisprDesign)
------------

                        Changes in version 1.5.2                        

- Added function addReinitiationFlag.

[crisprShiny](/packages/crisprShiny)
-----------

                        Changes in version 1.0.0                        

- New package crisprShiny, for interactive visualizations of
  GuideSet objects in facilitating the design of CRISPR gRNA.

[CTexploreR](/packages/CTexploreR)
----------

                        Changes in version 0.99                         

CTexploreR 0.99.5

- Address Bioc review comments.

CTexploreR 0.99.4

- Don't run some examples to save time. These functions are
illustrated in the vignette.

CTexploreR 0.99.3

- Reduce check time

CTexploreR 0.99.2

- Updated unit tests

CTexploreR 0.99.1

- Addressing the bioconductor comments

CTexploreR 0.99.0

- Package submission to Biocoductor.

[cTRAP](/packages/cTRAP)
-----

                       Changes in version 1.20.1                        

- When running cTRAP(), raise error if commonPath does not exist
- Fix modular graphical interface functions not showing dropdown
choices
- Update icon names for FontAwesome 6

[cypress](/packages/cypress)
-------

                Changes in version 0.99.24 (2024-04-23)                 

- update new methods for cedar and deseq2

                Changes in version 0.99.21 (2024-04-06)                 

- Coding modified in compliance with BioC requirement.

                 Changes in version 0.99.3 (2024-02-12)                 

- debug for warning

                 Changes in version 0.99.0 (2024-01-15)                 

- Submitted to Bioconductor

[cytomapper](/packages/cytomapper)
----------

                 Changes in version 1.15.4 (2024-03-21)                 

- measureObject fix reversion and test - Bioc 3.19

                 Changes in version 1.15.3 (2024-03-13)                 

- measureObject fixes for Bioc 3.19 release

                 Changes in version 1.15.2 (2023-12-13)                 

- change of package maintenance to Lasse Meyer

                 Changes in version 1.15.1 (2023-12-10)                 

- the unit tests for compImage have now been adjusted to account for
  correct spillover induction

[CytoMDS](/packages/CytoMDS)
-------

                        Changes in version 0.99                         

CytoMDS 0.99.16

- added lineWidth parameter in ggplotSampleMDSShepard()
- running plotly::ggplotly() on ggplotSampleMDSShepard() output now
displays row and column number for each distance point.
- added pointLabelSize and arrowLabelSize in ggplotSample()

CytoMDS 0.99.15

- corrected bug fix (error message) in pwDist() when verbose=TRUE

CytoMDS 0.99.14

- re-factored code portions to replace, as much as possible, for
loops
by apply() family of functions.

CytoMDS 0.99.13

- re-factored code portions to avoid growing lists incrementally

CytoMDS 0.99.12

- removed useBiocParallel parameters from various stats functions
(use
BPPARAM = BiocParallel::SerialParam() as a default)
- implemented MDS class to store MDS projection results
- bi-plots now explicitly discard constant external variables
(+warning) instead of raising an error without producing a plot
- implemented ggplotMarginalDensities()
- updated vignette with Bodenmiller2012 dataset and more biological
interpretation.

CytoMDS 0.99.11

- re-factored package documentation file

CytoMDS 0.99.10

- biplot now handles extVariables with missing values

CytoMDS 0.99.9

- in ggplotSampleMDS() : add label layer after geom_point() (no more
before)

CytoMDS 0.99.8

- renamed getChannelSummaryStats() into channelSummaryStats()
- in channelSummaryStats(), added support for BiocParallel`, and
allowed for not loading the whole flowSet in memory at once.
- replaced NULL defaulted parameters with optional parameters
- added displayPointLabels argument to ggplotSampleMDS()
- added displayLegend argument to ggplotSampleMDSWrapBiplots()
- finalized creating vignette

CytoMDS 0.99.7

- refactored the pairwise distance calculation code, by pre-computing
the unidimensional histograms and store them instead of
recalculating them each time a distance between 2 samples is
calculated. This improves CPU time and memory consumption.

CytoMDS 0.99.6

- added subset argument in ggplotSampleMDS() and
ggplotSampleMDSWrapBiplots

CytoMDS 0.99.5

- renamed getPairwiseEMDDist() into pairwiseEMDDist()
- in pairwiseEMDDist(), added support for BiocParallel, and allowed
for not loading the whole flowSet in memory at once.

CytoMDS 0.99.4

- in getPairwiseEMDDist(), added a second flowSet argument. When the
two flowSet arguments are non-null, distances are calculated for all
sample pairs, where the first element comes from fs, and the second
element comes from fs2.
- renamed ggplotSamplesMDS into ggplotSampleMDS
- renamed ggplotSamplesMDSShepard into ggplotSampleMDSShepard
- renamed getChannelsSummaryStat into getChannelSummaryStats
- new function ggplotSampleMDSWrapBiplots()

CytoMDS 0.99.3

- new version of computeMetricMDS() which automatically sets the
number of dimensions to reach a target pseudo R squared
- added ggplotly() functionality for output MDS plots
- in ggplotSampleMDS(), added flipXAxis, flipYAxis to possibly ease
low dimensional projection comparisons
- in ggplotSampleMDS(), added displayArrowLabels to discard the arrow
labels in biplot. Also added arrowThreshold. Moved arrow labels
toward the end of the arrows.
- in ggplotSampleMDS() and ggplotSampleMDSShepard(): added
displayPseudoRSq parameter.

CytoMDS 0.99.2

- use global Rsquare as an indicator of quality of projection
- use %Var explained per axis

CytoMDS 0.99.1

- in ggplotSamplesMDS(), added parameter pDataForAdditionalLabelling

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

                         Changes in version 1.3                         

CytoPipeline 1.3.6

- execute() now stores the nb of events retained at each
pre-processing step, to speed-up collectNbOfRetainedEvents()

CytoPipeline 1.3.5

- added CITATION file

CytoPipeline 1.3.4

- in execute(), when cache does already exist, make it clean before
executing the pipeline steps (= preventing inconsistent cache upon
crash/forced interruption)

CytoPipeline 1.3.3

- added collectNbOfRetainedEvents() function

CytoPipeline 1.3.2

- systematically override pData in cache upon execute() to allow
running consistently running several times for increasing number of
samples
- sampleFiles<- and pData<-: make sure that order of sample files
follow the one of pData rownames if pData exists.
- added 'verbose' argument in estimateScaleTransforms()

CytoPipeline 1.3.1

- refactored documentation files

[CytoPipelineGUI](/packages/CytoPipelineGUI)
---------------

                         Changes in version 1.1                         

CytoPipelineGUI 1.1.3

- re-factored package documentation file

CytoPipelineGUI 1.1.2

- added CITATION file

CytoPipelineGUI 1.1.1

- re-factored documentation files

[cytoviewer](/packages/cytoviewer)
----------

                 Changes in version 1.3.2 (2024-01-05)                  

- added pixel resolution option

                 Changes in version 1.3.1 (2024-01-04)                  

- updated README and added citation file


[dar](/packages/dar)
---

                       Changes in version 0.99.13                       

Bug Fixes

- Remove humann example form data import vignette

                       Changes in version 0.99.10                       

New Features

- Reimplementing step_corncob after the return of corncob package to
cran

                       Changes in version 0.99.9                        

Bug Fixes

- Set workers parameter to 4 in order to avoid issues with BBS builds

                       Changes in version 0.99.8                        

Improvements

- Reducing examples computation time

                       Changes in version 0.99.7                        

Improvements

- Reducing examples computation time

                       Changes in version 0.99.6                        

Improvements

- Suggest using BiocManager::install() to install dar dependencies

                       Changes in version 0.99.5                        

Improvements

- Reducing vignettes computation time

                       Changes in version 0.99.4                        

Improvements

- Reducing tests and examples computation time

Bug Fixes

- Fixing bug in Github Actions on Linux with rlang installation.

                       Changes in version 0.99.3                        

New Features

- The dar package now accepts both phyloseq class objects and
TreeSummarizedExperiment as inputs.
- The tutorial has been refocused to become a tutorial on how to
import biom, qiime, mothur, metaphlan, and humann into
TreeSummarizedExperiment and phyloseq class objects.
- The Recipe and PrepRecipe classes have been introduced, replacing
the previous recipe and prep_recipe classes.
- The subset and filter operations have been updated to allow all
steps of the recipe to be defined in a chainable manner.
- The functions step_filter_by_abundance, step_filter_by_prevalence,
step_filter_by_rarity, and step_filter_by_variance have been added
to enhance filtering functionality.

Improvements

- The R version dependency has been updated to 4.4.0.
- The dependency on data.table has been removed.
- The re-export of %>% and := has been removed. Now code examples and
vignettes use |>.
- The required_deps function is no longer exported.
- The package now recommends more commonly used installation methods,
such as BiocManager::install() or install.packages().
- A warning message is now displayed whenever the rarefy = TRUE
option
is used, informing users that a fixed seed is being used and how it
could impact their results.
- The package coverage has increased to 82.33%.

Bug Fixes

- Unconventional package installation methods have been avoided, for
example, pak::pkg_install.
- The setting of a seed within a function (run_aldex) has been
addressed.
- The name of the data in the R/data.R documentation has been
corrected from NA.

                       Changes in version 0.99.0                        

- Initial Bioconductor submission.

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

                       Changes in version 0.99.16                       

- Added utils to DESCRIPTION

                       Changes in version 0.99.15                       

- Updated to use R 4.4.0

- Added importFrom for grDevices and utils

- Fixed typos in help files

                       Changes in version 0.99.14                       

BUG FIXES

- Added imported packages to DESCRIPTION

                       Changes in version 0.99.13                       

- Numerous code style changes made to comply with Bioconductor
standards.

- Updated to run require R 4.3.3

- Vignette style changed to Bioconductor

- Selective imports implemented

                       Changes in version 0.99.12                       

BUG FIXES

- Fixed one more bad href in help files.

                       Changes in version 0.99.11                       

BUG FIXES

- Fixed bad hrefs in help files.

                       Changes in version 0.99.10                       

BUG FIXES

- Fixed example for runDegCre.

                       Changes in version 0.99.9                        

BUG FIXES

- Fixed test functions to subset GRanges.

                       Changes in version 0.99.8                        

BUG FIXES

- Fixed runnable examples to load GenomicRanges.

                       Changes in version 0.99.7                        

BUG FIXES

- Changed test functions to run on chr1 only for efficiency.

                       Changes in version 0.99.6                        

BUG FIXES

- Changed examples to run on chr1 for computational efficiency.

                       Changes in version 0.99.4                        

BUG FIXES

- Fixed F to FALSE in example

                       Changes in version 0.99.3                        

BUG FIXES

- Set all non-exported functions to not run examples.
- Fixed vignette code.

NEW FEATURES

- getAssocDistHits is now exported.

                       Changes in version 0.99.2                        

BUG FIXES

- Fixed calcAUC.Rd to not run example.

                       Changes in version 0.99.1                        

BUG FIXES

- Fixed calcBinomFDRperBin to give FDR = 1 for 0 assocProbs in top 10
CRE p-values.

                       Changes in version 0.99.0                        

BUG FIXES

- Added missing unit test for optimizeAlphaDegCre


[DEGreport](/packages/DEGreport)
---------

                       Changes in version 1.39.4                        

- Fix error in check report about stack limit size

                       Changes in version 1.39.3                        

- Add ctb

                       Changes in version 1.39.1                        

- merging changes from 1.38.1

                       Changes in version 1.38.1                        

- remove lasso and fdrtools

- merge https://github.com/lpantano/DEGreport/pull/56

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.30.0                        

NEW FEATURES

- rowsum(&lt;DelayedMatrix&gt;)/colsum(&lt;DelayedMatrix&gt;) now acknowledge the
  current "automatic realization backend" and "automatic BiocParallel
  BPPARAM". See '?DelayedArray::rowsum' for more information.

SIGNIFICANT USER-VISIBLE CHANGES

- Rename supportedRealizationBackends -> registeredRealizationBackends.

- Slightly modify the behavior of 'realize(x, BACKEND=NULL)'.
  See https://github.com/hansenlab/minfi/issues/256

- Two important changes to matrix multiplication of DelayedMatrix
  objects.
  1. Now it returns an ordinary matrix by default (before this change
  an
  ordinary matrix was returned wrapped in a DelayedMatrix object). The
  user can change the default behavior by setting an "automatic
  realization backend". See ?DelayedArray::`%*%` for more information.
  2. Better block processing strategy when only one of the two operands
  is
  a DelayedMatrix object (or derivative). The new strategy acknowledges
  the geometry of the physical chunks of the data in the object. This
  can make a huge difference in some cases. For example, using a subset
  of the "1.3 Million Brain Cell Dataset" from 10x Genomics:
  library(HDF5Array)
  library(ExperimentHub)
  hub <- ExperimentHub()
  tenx <- TENxMatrix(hub&#91;&#91;"EH1039"&#93;&#93;, group="mm10")
  M <- tenx&#91; , 1:25000&#93;
  m <- cbind(runif(ncol(M)), runif(ncol(M)))
  M %*% m
  Doing 'M %*% m' now takes 7.6s and uses 1.1Gb of memory, compared to
  110s / 3.1Gb before this improvement. Furthermore, the new strategy
  operates in linear time and at constant memory:
  with DelayedArray with DelayedArray
  ncol(M) 0.29.2 < 0.29.2
  ------- ----------------- -----------------
  12500 4.3s / 1.1Gb 32s / 2.1Gb
  25000 7.6s / 1.1Gb 110s / 3.1Gb
  50000 13.4s / 1.1Gb 495s / 5.6Gb
  100000 24.0s / 1.2Gg 2409s / 9.1Gb
  Note that the new strategy is implemented in internal helpers
  DelayedArray:::BLOCK_mult_Lgrid() and
  DelayedArray:::BLOCK_mult_Rgrid(). When the two operands are
  DelayedMatrix objects, the old strategy (which is implemented in
  DelayedArray:::.super_BLOCK_mult()) is still used.

[DelayedMatrixStats](/packages/DelayedMatrixStats)
------------------

                        Changes in version 1.25                         

- Scalar center argument for matrix functions are now defunct
  following similar change in matrixStats
  (<https://github.com/HenrikBengtsson/matrixStats/issues/254>).

[DepecheR](/packages/DepecheR)
--------

                 Changes in version 1.19.9 (2024-01-15)                 

- Bug fixes

[DeProViR](/packages/DeProViR)
--------

                       Changes in version 0.99.0                        

- Submitting the package

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.44.0                        

- New method for 'greaterAbs' in results() which has more
  power than the original 2014-2023 method. The old method
  is available as 'greaterAbs2014'.
  Suggested by Nikos Ignatiadis

- Fix from I-Hsuan Lin for the glmGamPoi dispersion estimation
  where the wrong indexing of the fitted mean matrix was used,
  which caused a slowdown.

- Fix from Rasmus Henningsson where the Cook's distances for
  the LRT were computed using a rough estimate of the mean,
  rather than the one from the GLM estimates of the full model.
  Now Cook's distances for Wald and LRT should be consistent.

[DESpace](/packages/DESpace)
-------

                        Changes in version 1.2.1                        

- Update to RELEASE_3_18

[DExMA](/packages/DExMA)
-----

                       Changes in version 1.10.1                        

- Improvement of makeHeatmap function

- Improvement Effects size meta-analysis methods and results

[DifferentialRegulation](/packages/DifferentialRegulation)
----------------------

                        Changes in version 2.0.2                        

- Updated from C++11 to C++17

                        Changes in version 2.0.1                        

- traceplot added, with functions plot_traceplot and
  plot_bulk_traceplot.
  MAJOR UPDATE

[dinoR](/packages/dinoR)
-----

                       Changes in version 0.99.0                        

- initial Bioconductor submission

[distinct](/packages/distinct)
--------

                       Changes in version 1.14.5                        

- Linux Rfast BUILD error solved

                       Changes in version 1.14.2                        

- Updated from C++11 to C++17

                       Changes in version 1.14.1                        

- ERROR solved regarind Rfast dependency

[dittoSeq](/packages/dittoSeq)
--------

                        Changes in version 1.16                         

- Feature Extensions:
1.  Multi-modality functionality: To support visualization of
markers from multiple modalities in the same plot, e.g. gene
expression by RNA and protein capture by ADT, the mechanics of
'assay', 'slot', and 'swap.rownames' inputs have been expanded,
although defaults are unchanged. See the '?GeneTargeting'
documentation page for details. For the standard Seurat CITEseq
case, set 'assay = c("RNA", "ADT")'.
2.  'dittoDotPlot()' vars-categories: Added support for
categorization of markers, as well as for x and y axes swapping.
- Provide 'vars' as a named list to group markers (list
element values) into categories (list element names).
- New input 'vars.dir' controls which axis is used for markers
("x" by default, or "y").
- New boolean inputs 'categories.theme.adjust' or
'categories.split.adjust' can be used to turn off associated
automated additions to the 'theme' input or 'split.adjust'
input as well as faceting mechanics, respectively.
3.  'dittoDotPlot()' 3-color scaling: Added support for injecting a
midpoint color to the color scale via 2 new inputs.
- New input 'mid.color' acts as the switch, and can be set to
the specific strings "ryb", "rwb", or "rgb" (g for gray
here) for a single-point quick update to use corresponding
ColorBrewer inspired scales (effectively updating
'min.color' and 'max.color' as well). 'mid.color' can
alternatively be given a color directly for more fine-grain
control of colors.
- New input 'mid' controls the data value at which 'mid.color'
will be used in the scale.
4.  Additional data representation controls for all
'dittoPlot()'-style plotters, which includes 'dittoFreqPlot()':
- New input 'boxplot.outlier.size' allows control of the
outlier shape's size for "boxplot" representations.
- New input 'vlnplot.quantiles' allows addition of lines at
requested data quantiles for "vlnplot" representations.
5.  Added a new built in 'color.method' style for
'dittoScatterHex()' and 'dittoDimHex()' plotters:
- When 'color.var' targets discrete data, giving 'color.method
= "prop.&lt;value&gt;"', where &lt;value&gt; is an actual data level of
'color.var'-data, will set coloring to represent the
proportion of &lt;value&gt; among the 'color.var'-data of each
bin.
6.  New input 'labels.repel.adjust' allows finer control of the
'do.label' plot addition, via input pass-through to the geom
functions underlying labeling. This affects 'dittoDimPlot()',
'dittoScatterPlot()', 'dittoDimHex()', and 'dittoScatterHex()
functions.
- Bug Fixes:
1.  'dittoHeatmap()': Fixed a bug which blocked provision of
'annotation_row' and 'annotation_colors' inputs to
'dittoHeatmap()' without also generating column annotations via
either 'annot.by' or direct 'annotation_col' provision.
- Deprecation:
1.  Completed deprecation of 'dittoHeatmap()'s 'highlight.genes'
input via removal from the function.
- Dependency Upkeep (generally invisible to users):
1.  ggplot-v3: Replaced all calls to the deprecated 'aes_string()'
function with calls to the standard 'aes()' function. In cases
where mappings are successively built internally to accommodate
customization or flexibility, 'modifyList()' usage replaces the
previous simple 'list' and 'do.call()' management.
2.  Seurat-v5: When the user's Seurat package version is 5.0 or
higher, conditional code switches expression data retrieval from
a call to the reportedly superseded 'GetAssayData()' function to
the newly supported 'SeuratObj&#91;&#91;&lt;assay&gt;&#93;&#93;&#91;&lt;slot&gt;&#93;' syntax.

[DMRcate](/packages/DMRcate)
-------

                        Changes in version 3.0.0                        

- Full utility for EPICv2 implemented. DMRs can now be called from the
  new Illumina Infinium MethylationEPIC v2.0 BeadChip same as the usual
  pipeline, save for replicate probe filtering (mandatory) and
  remapping of cross-hybridising probes (optional).

- Annotation package EPICv2manifest is used as a backend for
  cpg.annotate().

- A new function, rmPosReps(), gives multiple user options for
  filtering replicate probes mapping to the same CpG site. The mean can
  be taken, or, based on Peters et al. (2024) (see documentation), the
  probe that is most precise or sensitive to methylation change may be
  selected.

- Many thanks to Braydon Meyer and Ruth Pidsley for constructive input
  and beta testing.

                       Changes in version 2.16.1                        

- DSS dependency removed

[DOSE](/packages/DOSE)
----

                       Changes in version 3.29.2                        

- fix bugs in get_ont_info() and get_dose_data(): wrong object name
and wrong data type (2023-11-30, Thu)

                       Changes in version 3.29.1                        

- mv 'MPO.db' and 'HPO.db' from 'Imports' to 'Suggests' and fixed
bugs
(2023-11-18, Sat)

[dreamlet](/packages/dreamlet)
--------

                       Changes in version 1.1.18                        

- bug fix

                       Changes in version 1.1.17                        

- new functionality in plotPCA() and outlierByAssay()
- works on any list, not just dreamletProcessedData
- allows outlier analysis on residuals

                       Changes in version 1.1.16                        

- Feb 8, 2024
- fix bug in pbWeights()
- Fix scaling issue in outlierByAssay()

                       Changes in version 1.1.15                        

- Feb 5, 2024
- Fix bug in call to eBayes()
- in processAssays() pass argument scaledByLib to
voomWithDreamWeights()

                       Changes in version 1.1.14                        

- Jan 29, 2024
- fix bug in pbWeights()
- smaller pseudo variance
- limit to only expressed genes by adding getExprGeneNames()

                       Changes in version 1.1.13                        

- Jan 25, 2024
- improve error reporting in seeErrors() and documentation
- update outlier() to compute z-scores. How returns data.frame()
- add outlierByAssay() and plotPCA()

                       Changes in version 1.1.12                        

- Jan 16, 2024
- compositePosteriorTest() allows exclude set to be NULL

                       Changes in version 1.1.11                        

- Jan 16, 2024
- add meta_analysis()

                       Changes in version 1.1.10                        

- Jan 10, 2024
- stackAssays() now includeds metadata()$aggr_means correctly
- add compositePosteriorTest()

                        Changes in version 1.1.9                        

- Jan 3, 2024
- use get_metadata_aggr_means() to extract aggr_means when SCE is
produced by cbind'ing

                        Changes in version 1.1.8                        

- Dec 18, 2023
- fix issue when no genes pass cutoffs
- fix issue with aggr_means in aggregateToPseudoBulk()
- fix bug when rdf is low for all genes

                        Changes in version 1.1.7                        

- Dec 10, 2023
- add plotBeeswarm()
- add rowWeightedVarsMatrix()
- bug fixes
- add isFullRank() check in dreamlet()
- handle exceptions in run_mash()

                        Changes in version 1.1.6                        

- Dec 5, 2023
- pbWeights() add argument maxRatio

                        Changes in version 1.1.5                        

- Nov 28, 2023
- fix edge cases in pbWeights()
- cell weights is not default in dreamlet()

                        Changes in version 1.1.4                        

- Nov 27, 2023
- add pbWeights() to compute precision weights for pseudobulk
counts
- extend extractData() and include it in vignette

                        Changes in version 1.1.3                        

- Nov 22, 2023
- add stackAssays()
- add diffVar()
- fix getVarFromCounts) so zeta is a mean, not a sum

                        Changes in version 1.1.2                        

- Nov 13, 2023
- computeLogCPM() uses augmentPriorCount()

                        Changes in version 1.1.1                        

- bug fix in Bioc 3.18 and devel

[edgeR](/packages/edgeR)
-----

                 Changes in version 4.2.0 (2024-04-28)                  

- 
  The new QL pipeline becomes the default for glmQLFit() by
  setting `legacy=FALSE`.

- 
  Add cameraPR method for DGELRT objects.

- 
  New arguments `prior.n` and `adaptive.span` for voomLmFit().

- 
  New argument `robust` for diffSpliceDGE().

- 
  The NEWS.Rd file has been revised to include the date of each
  version release and to include earlier versions of edgeR.

- 
  catchSalmon() now detects whether resamples are Gibbs or
  bootstrap.

- 
  The catchSalmon help page now explains the columns of the
  `annotation` output data.frame.

- 
  decideTestsDGE() deprecated in favor of decideTests().

[eisaR](/packages/eisaR)
-----

                       Changes in version 1.15.1                        

- Add legacyQLF argument to runEISA (will be passed to
edgeR::glmQLFit)

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.23.2                        

- separate the JC similarity method (2023-12-11, Mon, #265)
- fix the issue in ridgeplot(showCategory) : support a vector of
Description, not ID(2023-12-1, Fri, #193)

                       Changes in version 1.23.1                        

- ridgeplot() supports passing a vector of selected pathways via the
'showCategory' parameter (2023-11-30, Thu, #193)
- fix treeplot() to compatible with the current version of ggtree and
ggtreeExtra. (2023-10-28, Sat)
- add clusterPanel.params&#91;&#91;"colnames_angle"&#93;&#93; parameter to set the
angle of colnames. (2023-10-28, Sat)

[enrichViewNet](/packages/enrichViewNet)
-------------

                        Changes in version 1.1.2                        

NEW FEATURES

- The new function createEnrichMapMultiBasic() enables the creation of
  enrichment maps with groups from simple designs.

- The new function createEnrichMapMultiComplex() enables the creation
  of enrichment maps with groups from complex designs.

SIGNIFICANT USER-VISIBLE CHANGE

- The vignette contains new sections presenting the new functions
  createEnrichMapMultiBasic() and createEnrichMapMultiComplex().

                        Changes in version 1.1.1                        

SIGNIFICANT USER-VISIBLE CHANGE

- The enrichViewNet workflow figure, in the vignette, has been updated.

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.27.1                        

- Functions for coordinate mapping gain support for pre-loaded genomic
  ranges
  hence enabling parallel processing support (contribution from Boyu
  Yu).

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.11.9 (2024-03-21)                 

- uses less memory (due to packed SEQ and XM)

                 Changes in version 1.11.8 (2024-03-12)                 

- stricter filtering in lMHL reports

                 Changes in version 1.11.7 (2024-03-05)                 

- optimised reporting from long-read data

                 Changes in version 1.11.6 (2024-02-29)                 

- float and array tags in simulateBam

- long-read data input (not optimised yet)

                 Changes in version 1.11.5 (2024-02-12)                 

- RRBS-ready

                 Changes in version 1.11.4 (2024-02-11)                 

- methylation calls for bsmap

[epiregulon](/packages/epiregulon)
----------

                       Changes in version 0.99.11                       

- new function aggregateAcrossCells has been added, which is
implemented in c++ and makes run time for calculateP2G much shorter.

                       Changes in version 0.99.5                        

- new function addLogFC which adds log fold changes of gene
expression
to regulons and significance statistics for differential gene
expression
- cellNum argument to calculateP2G set by default to 100 (previously
it was 200).
- in pruneRegulon the check to the uniqueness of gene names has been
added.

                       Changes in version 0.99.0                        

Version number downgraded to 0.99.0 to meet Bioconductor
requirements.


[epiregulon.extra](/packages/epiregulon.extra)
----------------

                       Changes in version 0.99.16                       

- Vignette updated

                       Changes in version 0.99.5                        

- Correct duplication error in plotHeatmapRegulon
- default direction is changed to "any" in getSigGenes

                       Changes in version 0.99.3                        

- default direction is changed to "any" in findDifferentialActivity
- The previous version regulonEnrich function required at least 3
target genes for a transcription factor to perform geneset
enrichment test. Now additional constraint is imposed such that we
require at least 3 genes to be present in the geneset.

                       Changes in version 0.99.2                        

[EpiTxDb](/packages/EpiTxDb)
-------

                 Changes in version 1.15.2 (2024-03-18)                 

- Bugfix due to upstream changes in internal functions

[erccdashboard](/packages/erccdashboard)
-------------

                       Changes in version 1.37.2                        

BUG FIXES AND MINOR IMPROVEMENTS

- Added legacy=TRUE for glmQLFit() function from edgeR to maintain
  behavior equivalent to results with Bioconductor 3.16, rather than
  using the new method of adjustment for small counts from edgeR that
  was introduced in Bioconductor 3.18 (legacy = FALSE by default).
  Changes to the glmQLFit function in edgeR seem to be causing slight
  changes to dispersion results even with legacy=TRUE.

[escape](/packages/escape)
------

                 Changes in version 1.99.1 (2024-02-29)                 

UNDERLYING CHANGES

- ordering by mean values no longer changes the color order
- add explicit BPPARAM argument to runEscape() and escape.matrix()
- added additional details in runEscape() and escape.matrix() for
make.positive.
- removed plotting of splitEnrichment() for group.by = NULL
- separated AUC calculation to rankings and AUC, this was only method
found to get consistent scores.

                 Changes in version 1.99.0 (2024-02-27)                 

NEW FEATURES

- Added runEscape()
- Added geyserEnrichment()
- Added scatterEnrichment()
- Added heatmapEnrichment()
- Changed enrichIt to escape.matrix()
- Changed enrichmentPlot to densityEnrichment()
- performPCA() now works with a matrix or single-cell object
- pcaEnrichment() combines biplot-like functions

UNDERLYING CHANGES

- Updated interaction with gsva package
- Added support for GSVA calculation
- Added support for AUCell calculation
- Added support of visualizations and calculations for single-cell
objects
- Modified getGeneSets() to output a list of gene set objects with
reformatted names following the Seurat "-" convention

DEPRECATED AND DEFUNCT

- Deprecate getSignificance()
- Deprecate masterPCAPlot()

[escheR](/packages/escheR)
------

                        Changes in version 1.3.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add functions add_fill_bin and add_ground_bin to provide hexgon
binning strategy to mitigate the overplotting of data points
- Provide an example to use add_fill_bin and add_ground_bin with
make_escheR to create hexgon binning plot

[esetVis](/packages/esetVis)
-------

                       Changes in version 1.29.2                        

- interactive plot: replace rbokeh (archived) by plotly

- doc: use internal keyword for internal functions

                       Changes in version 1.29.1                        

- vignette: run rbokeh example only if package is available

[EWCE](/packages/EWCE)
----

                       Changes in version 1.11.5                        

Bug fixes

- check_bootstrap_args:
- annotation level for checking cell type names was previous
hardcoded to level 1, now updated to match user input for
annotation level.

                       Changes in version 1.11.4                        

Bug fixes

- drop_uninformative_genes:
- Will now catch cases where expression matrix is a dataframe and
convert to a matrix
- This was causing weird errors, see issue 92:
https://github.com/NathanSkene/EWCE/issues/92
- Fix made in check_sce() function.

                       Changes in version 1.11.2                        

New features

- bootstrap_enrichment_test:
- New args: standardise_sct_data=, standardise_hits=: let users
have more control over data standardisation steps.
- check_ewce_genelist_inputs: updated accordingly.
- New arg: store_gene_data to avoid hitting memory limits.
- Modify test-bootstrap_enrichment_test_2.R to use test new args.

                       Changes in version 1.11.1                        

Bug fixes

- ewce_plot() - Dendrogram not reordering cell types in plot
- see issue
- Occurs when ctd does not contain the plotting info
- Fixed now and unit test added.
- Note that cell type order on the x-axis is based on hierarchical
clustering for both plots if make_dendro = TRUE for ewce_plot().

[ExperimentHub](/packages/ExperimentHub)
-------------

                       Changes in version 2.11.0                        

USER VISIBLE MODIFICATIONS

- (2.11.2) Update to only install package if interactive and approved
  by user

BUG FIXES

- (2.11.3) Herve's fix. No more call to .get_ExperimentHub() in
  createHubAccessors() itself, only in the accessor functions that it
  creates

[ExploreModelMatrix](/packages/ExploreModelMatrix)
------------------

                       Changes in version 1.15.1                        

- Add 'Close app' button to close the application and export the data

[extraChIPs](/packages/extraChIPs)
----------

                        Changes in version 1.7.7                        

- Added merge_within to makeConsensus() for better handling when
method = "coverage"

                        Changes in version 1.7.6                        

- Added the DESeq2 Wald statistic to options for fitAssayDiff()

[faers](/packages/faers)
-----

                       Changes in version 0.99.0                        

- Initial Bioconductor submission.

[fastreeR](/packages/fastreeR)
--------

                 Changes in version 1.7.3 (2024-04-07)                  

- Preparing for next Bioconductor Release.

                 Changes in version 1.7.2 (2024-02-10)                  

- Fix vignette bug.

                 Changes in version 1.7.1 (2024-02-10)                  

- Update R version requirement.

                 Changes in version 1.7.0 (2023-10-29)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_18
  branch.

[fenr](/packages/fenr)
----

                       Changes in version 1.0.10                        

- Due to recurring issues with build and check on Bioconductor's
machines, I have removed all database downloads from the vignette.
Any glitch in the GO server, or simply an internet problem would
cause the vignette build to crash. The GO-term information is now
attached as data and loaded in the vignette.
- Made sure the package passes BUILD and CHECK with no internet
connection.
- Correction in vignette: using yeast genome for topGO, instead of
human (somehow it was not applied in 1.0.5).

                        Changes in version 1.0.9                        

- Changed the way assert_url_path() handles some remote files - it
turns out every time it was called, the entire file was unnecessary
downloaded, leading to duplication. Now we only assert top
directories. Should speed things up!
- Increased default timeout to 30 s.

                        Changes in version 1.0.8                        

- Further improving error handling, making sure assert_url_path()
handles timeouts properly
- Introduced on_error = "ignore" for test purposes

                        Changes in version 1.0.7                        

- Improved error handling with unresponsive servers - timeouts are
now
handled gracefully

                        Changes in version 1.0.6                        

- Changed the Ensembl mapping file downloaded from Reactome to
"Physical entity" mapping, as it contains gene symbols, in addition
to the Ensembl IDs.
- Changed the name of GAF column DB Object Synonym from gene_synonym
to gene_id for consistency with other methods.
- Corrected Reactome test as it failed with multiple gene symbols per
gene id.
- Replaced biomaRt with a single RESTful XML call; as biomaRt is used
only once to obtain GO terms, this replacement reduced dependency
footprint of the package

                        Changes in version 1.0.5                        

- Bug fix: if feature id - term id mapping is not unique (which can
happen), features are duplicated in counting; fixed by
dplyr::distinct() on mapping
- Correction in vignette: using yeast genome for topGO, instead of
human.
- Improving test coverage

                        Changes in version 1.0.4                        

- Reinstated Bioplanet access, this time with graceful fail when the
website is down.
- Minor code changes.

                        Changes in version 1.0.2                        

- Bug fixes, examples need on_error = "warn"

                        Changes in version 1.0.1                        

- First update after Bioconductor release
- Implemented changes to prevent the package from build/check fail,
if
one of the remote servers is not responding
- Moved from httr to httr2
- Tests and examples now generate warnings in case of server failure
- Added tests for behaviour in case of a non-responsive server
- Extended test coverage to 100%, except for the interactive example

[findIPs](/packages/findIPs)
-------

                       Changes in version 0.99.3                        

SIGNIFICANT USER-VISIBLE CHANGES

o Format News

                       Changes in version 0.99.2                        

NEW FEATURES

o Update R version dependency from 4.0.0 to 4.4.0

o Use more importFrom()

BUG FIXES

o Correction on sumRanks()

SIGNIFICANT USER-VISIBLE CHANGES

o Exchange the axis in plotRankScatters()

o Include Installation in Vignette

o Some formatting modification

                       Changes in version 0.99.1                        

SIGNIFICANT USER-VISIBLE CHANGES

o Some formatting modification

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

                       Changes in version 1.11.4                        

- Fix bugs in vignettes

[FRASER](/packages/FRASER)
------

                        Changes in version 2.0.0                        

- Bugfix for contig names containing certain characters

- Update of the plot functions to support colorring aberrant status
  based
  on p values computed on subsets of genes

- Major update to FRASER2:

- Introduction of new & more robust splice metric Intron Jaccard Index

- Only Intron Jaccard Index metric used by default

- Improved gene level pvalue calculation and internal storage

- Introduction of option to limit FDR correction to user-defined
  subsets of genes per sample (e.g. OMIM genes with rare variant)

- Updated internal pseudocount parameter and default delta Jaccard
  cutoff

- Junction filtering adapted to usage of Intron Jaccard Index metric

- Require min expression of N >= 10 in 25% of the samples

- Results table:

- Functionality to flag outliers in blacklist regions of the genome

- Functionality to annotate the predicted type of aberrantSplicing
  (e.g. exon skipping, intron retention etc.)

- Several updates in the plotting functions

- introduction of manhattan plot functionality

- possibility to create sashimi plots to visualize read coverage in
  the bam files for outliers

[gDNAx](/packages/gDNAx)
-----

                 Changes in version 1.2.0 (2024-04-27)                  

USER VISIBLE CHANGES

- Prompt an error when a 0-byte BAM file is passed as an argument.

- Added a new filtering strategy based on stranded genomic windows that
  can be used when input RNA-seq data is stranded.

- gDNAx objects store now also strandedness estimated values.

- identifyStrandMode() has been deprecated in favor of strandedness().

- gDNAdx() can now automatically guess library layout and protocol.

- Added a new function gDNAtx() that performs filtering with the most
  common default strategy.

BUG FIXES

- Bug fix in the show method for gDNAx objects.

- Bug fix when the genome() getter on GAlignment* objects is not NA.

[gDR](/packages/gDR)
---

               Changes in version 2024-02-26 (2024-02-26)               

- improve pkgdown site
- improved references
- valid NEWS.md

               Changes in version 2024-02-12 (2024-02-12)               

- make documentation compatible with pkdgdown

               Changes in version 2024-01-22 (2024-01-22)               

- add new description fields

               Changes in version 2024-01-19 (2024-01-19)               

- update package vignette

               Changes in version 2023-11-22 (2023-11-22)               

- sync main with devel branch

               Changes in version 2023-10-24 (2023-10-24)               

- release Bioc 3.18

- prerelease Bioc 3.18

               Changes in version 2023-10-18 (2023-10-18)               


[gDRcore](/packages/gDRcore)
-------

               Changes in version 2024-04-23 (2024-04-23)               

- add vignette with documentation for data annotation

               Changes in version 2024-04-22 (2024-04-22)               

- add support for external annotation specified in the env var

               Changes in version 2024-03-26 (2024-03-26)               

- fix bug with converting mae into raw_data

               Changes in version 2024-03-15 (2024-03-15)               

- remove unstable tests

               Changes in version 2024-03-14 (2024-03-14)               

- cleanup package

               Changes in version 2024-03-12 (2024-03-12)               

- update function description

               Changes in version 2024-02-26 (2024-02-26)               

- improve pkgdown site
- improved references
- valid NEWS.md

               Changes in version 2024-02-14 (2024-02-14)               

- fix issue with retrieving unique records from mix of control and
treated samples

- make documentation compatible with pkdgdown
- rename 'matches' to 'grr_matches'

               Changes in version 2024-02-12 (2024-02-12)               

- fix unit tests for GRAN

               Changes in version 2024-02-07 (2024-02-07)               

- simplify logic of assays for combination data
- rename matrix into combination

               Changes in version 2024-02-06 (2024-02-06)               

- add support for internal source of annotation

- fix bug with converting standardize MAE into raw data

               Changes in version 2024-02-05 (2024-02-05)               

- add vignette for a data model

               Changes in version 2024-02-01 (2024-02-01)               

- update wrappers for co-dilution data

               Changes in version 2024-01-22 (2024-01-22)               

- add new description fields

               Changes in version 2024-01-04 (2024-01-04)               

- improve logic for normalization and identification of single-agent
and matrix data with Drug3

               Changes in version 2023-12-15 (2023-12-15)               

- fix issue with wrong assignment of untreated records

               Changes in version 2023-11-22 (2023-11-22)               

- sync master with devel branch
- add support for unifying duplicates in combo matrix data
- add "Treatment" as template identifier

               Changes in version 2023-10-24 (2023-10-24)               

- release Bioc 3.18

- prerelease Bioc 3.18

               Changes in version 2023-10-17 (2023-10-17)               


[gDRimport](/packages/gDRimport)
---------

               Changes in version 2024-04-08 (2024-04-08)               

- change output of get_exception_data to data.table

               Changes in version 2024-03-26 (2024-03-26)               

- fix issue with reading tsv files

               Changes in version 2024-03-07 (2024-03-07)               

- clean up the package

               Changes in version 2024-02-26 (2024-02-26)               

- improve pkgdown site
- improved references
- valid NEWS.md

               Changes in version 2024-02-14 (2024-02-14)               

- make documentation compatible with pkdgdown

               Changes in version 2024-01-22 (2024-01-22)               

- add new description fields

               Changes in version 2023-12-19 (2023-12-19)               

- update package vignette

               Changes in version 2023-12-01 (2023-12-01)               

- update validation function used during data submission

               Changes in version 2023-11-22 (2023-11-22)               

- sync main with devel branch
- add "Treatment" as template identifier

               Changes in version 2023-10-24 (2023-10-24)               

- release Bioc 3.18

- prerelease Bioc 3.18

               Changes in version 2023-10-17 (2023-10-17)               


[gDRstyle](/packages/gDRstyle)
--------

               Changes in version 2024-03-25 (2024-03-25)               

- add as_cran flag

               Changes in version 2024-03-20 (2024-03-20)               

- replace remotes:::version_satisfies_criteria function

               Changes in version 2024-03-18 (2024-03-18)               

- temporary skip verifying version test on Bioc

               Changes in version 2024-03-04 (2024-03-04)               

- remove ::: from notes exceptions

               Changes in version 2024-02-26 (2024-02-26)               

- improve pkgdown site
- improved references
- valid NEWS.md

               Changes in version 2024-02-12 (2024-02-12)               

- update documentation - pkdown compatibility

               Changes in version 2024-01-16 (2024-01-16)               

- adjust lintr configs

               Changes in version 2023-11-22 (2023-11-22)               

- sync main with devel branch

               Changes in version 2023-10-24 (2023-10-24)               

- release Bioc 3.18

- prerelease Bioc 3.18

               Changes in version 2023-10-17 (2023-10-17)               


[gDRutils](/packages/gDRutils)
--------

               Changes in version 2024-04-15 (2024-04-15)               

- add get_testdata_combo and get_testdata_codilution

               Changes in version 2024-03-07 (2024-03-07)               

- clean up the package

- simplify keywords

               Changes in version 2024-02-28 (2024-02-28)               

- add fit_source to header list

               Changes in version 2024-02-26 (2024-02-26)               

- improve pkgdown site
- improved references
- valid NEWS.md

               Changes in version 2024-02-22 (2024-02-22)               

- restore tooltips in table

               Changes in version 2024-02-14 (2024-02-14)               

- make documentation compatible with pkdgdown

               Changes in version 2024-01-30 (2024-01-30)               

- rename matrix into combination

               Changes in version 2024-01-22 (2024-01-22)               

- add new description fields

               Changes in version 2023-12-01 (2023-12-01)               

- fix bug with refining rowData
- extend the list of headers

               Changes in version 2023-11-22 (2023-11-22)               

- sync master with devel branch
- update schema to support NA in reference division time
- add minor fix in code styling
- add new function gemoetric_mean
- transform values into numeric in predict_efficacy_from_conc
function
- add "Treatment" as template identifier

               Changes in version 2023-10-24 (2023-10-24)               

- release Bioc 3.18

- prerelease Bioc 3.18

               Changes in version 2023-10-18 (2023-10-18)               


[GDSArray](/packages/GDSArray)
--------

                       Changes in version 1.23.1                        

BUG FIXES

- Removed `gdsfn<-` filename setter for GDSArray and
  GDSArraySeed. Only valid for `GDSFile` class.

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.38.1                        

UTILITIES

- fix the compiler warning: -Wformat-security

[GeDi](/packages/GeDi)
----

                       Changes in version 0.99.5                        

- This version reflects further changes performed upon the
Bioconductor reviewing process.
- Removed icons from the vignette to reduce the overall size of the
package.

                       Changes in version 0.99.4                        

- This version reflects the changes performed upon the Bioconductor
reviewing process
- Added col_name_genesets and col_name_genes as parameter to the
GeDi() main app to allow users to specify the relevant column names
upon executing the command
- Changes in the R code to comply to best practices (replacing single
| with || and similar)
- Loading the example file does not require anymore the setting of
globalVariables()
- All files retrieved do use some form of caching for avoiding
unnecessary re-download operations
- Reworked the allocation of vectors before for loops to avoid
unhealthy growing of vectors/matrices

                       Changes in version 0.99.1                        

- The handling of the parallelization for the distance calculations
is
now unified under the umbrella of BiocParallal, and defaults now to
using SerialParam() to avoid unexpected behaviors

                       Changes in version 0.99.0                        

- Ready for the submission to Bioconductor!

                       Changes in version 0.90.0                        

- Final touches and bug fixes to the main functionality
- Deployment of the package website via pkgdown

                        Changes in version 0.1.0                        

- Officially entering the path of the GeDi!

[gemma.R](/packages/gemma.R)
-------

                        Changes in version 3.0.0                        

- Empty outputs now return data.tables with the corresponding column
names with 0 rows instead of defaulting to lists of length 0
- Column names for the outputs of many functions have changed in this
release to be more standardized. Please refer to the function
documentation. As a general rule names use camelCase separated by .s
to indicate properties of a specific entity (eg.
experiment.sampleCount) with the exception of acronyms which are
always capitalized (experiment.ID)
- update_results function added which allows re-creation of outputs
of
gemma.R functions without relying on the original code.
- get_result_sets function added which allows accessing result sets
directly, filtering them based on certain filterable properties (see
filter_properties()$resultSet).
- gemma_memoise function added which allows setting memoisation
options without manually setting options
- get_child_terms function added which returns child terms of an
ontology term as inferred by Gemma
- gemma_kable is added which returns tables formatted to fit

[GENESIS](/packages/GENESIS)
-------

                       Changes in version 2.33.2                        

- Set extremely small p-values (< Machine\$double.xmin) calculated
  with pchisq to Machine$double.xmin. This change prevents
  GENESIS from returning p-values equal to 0.

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 2.8.0                        

Other notes

- Updated functions used from other packages to reflect changes in
their API/nomenclature
- Adapted the internal code of functions to the latest version of
igraph - no changes happening for the end user

[GenomAutomorphism](/packages/GenomAutomorphism)
-----------------

                        Changes in version 1.5.1                        

- Introducing new functions for DNA and aminoacid sequence
representations with physicochemical properties of DNA and
aminoacids, which would be useful for further downstream statistical
analysis in R.

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.40.0                        

NEW FEATURES

- Register the following NCBI assemblies:
  - GRCr8 assembly (Rat)
  - mCavPor4.1 assembly (domestic guinea pig)
  - 21 Escherichia coli assemblies
  - a few Petromyzon marinus (sea lamprey) assemblies

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.40.0                        

NEW FEATURES

- Add 'seqinfo' argument to GAlignments() constructor function.

BUG FIXES

- mapToAlignments() S at start 0 coordinate fix. By Fedor Bezrukov.
  See https://github.com/Bioconductor/GenomicAlignments/pull/34

- Bugfix in isCompatibleWithSkippedExons() for single-end read
  alignments
  w/ more than 3 junctions. By Robert Castelo.
  See https://github.com/Bioconductor/GenomicAlignments/pull/32

[GenomicDataCommons](/packages/GenomicDataCommons)
------------------

                       Changes in version 1.28.0                        

- Defunct legacy function, methods, endpoints, and arguments


[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.56.0                        

NEW FEATURES

- Add getTerminatorSeq(). Same as getPromoterSeq() but for terminator
  sequences.

SIGNIFICANT USER-VISIBLE CHANGES

- The makeTxDb*() functions and related have moved to the new txdbmaker
  package. Full list:
  - makeTxDb
  - supportedUCSCtables
  - browseUCSCtrack
  - makeTxDbFromUCSC
  - getChromInfoFromBiomart
  - makeTxDbFromBiomart
  - makeTxDbFromEnsembl
  - makeTxDbFromGRanges
  - makeTxDbFromGFF
  - supportedUCSCFeatureDbTracks
  - supportedUCSCFeatureDbTables
  - UCSCFeatureDbTableSchema
  - makeFeatureDbFromUCSC
  - supportedMiRBaseBuildValues
  - makePackageName
  - makeTxDbPackage
  - makeTxDbPackageFromUCSC
  - makeFDbPackageFromUCSC
  - makeTxDbPackageFromBiomart
  Note that they are still temporarily defined in GenomicFeatures but
  now
  they just call the corresponding function in txdbmaker. Since this is
  a
  temporary redirect, the user also gets a warning that tells them to
  use
  the fully qualified name (e.g. txdbmaker::makeTxDbFromUCSC()) to call
  the function.

[GenomicPlot](/packages/GenomicPlot)
-----------

                       Changes in version 1.1.10                        

- Handle NCBI style of seqlevels (seqname)

                        Changes in version 1.1.9                        

- Fixed a bug in plot_region

                        Changes in version 1.1.8                        

NEW FEATURES

- Add handling of .gz files for bed and bedGraph format.
- Add input filter for chromosomes. Only specified chromosomes
will be included in visualization and analysis.

                        Changes in version 1.1.7                        

NEW FEATURES

- Add handle_bedGraph to enable data input of bedGraph format

                        Changes in version 1.1.6                        

- Merged a pull request from Hervé Pagès notifications@github.com

                        Changes in version 1.1.5                        

- Added a function to obtain chromosome info from cached data in the
'circlize' package to produce a Seqinfo object, which is applied to
all GRanges and TxDB object
- Added another function to generate a customized TxDb object from
a genome annotation (GTF or GFF) file.

                        Changes in version 1.1.4                        

Change chromosome size information source from UCSC web service to
cached data in the 'circlize' package, to avoid internet connection
issues and web service issues.

                        Changes in version 1.1.3                        

Increase font size for axis labels, align profile with heatmap

                        Changes in version 1.1.1                        

Removed some suggests in DESCRIPTION to reduce installation time and
size

BUG FIXES

Fixed misalignment between profiles and heatmaps

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.56.0                        

NEW FEATURES

- Add terminators() method, same as promoters() but for terminator
  regions.

BUG FIXES

- Small fix in makeGRangesFromDataFrame(). Fix error when: (1) input
  has
  zero rows, and (2) has no strand field or 'ignore.strand' is TRUE.
  By Marcel Ramos.
  See https://github.com/Bioconductor/GenomicRanges/issues/81

[GenomicScores](/packages/GenomicScores)
-------------

                       Changes in version 2.16.0                        

USER VISIBLE CHANGES

- Added support to the latest version v4.0 of gnomAD MAF data, stored
  in the package MafH5.gnomAD.v4.0.GRCh38.

- Updated vignette to showcase the use of pathogenicity score sets from
  AlphaMissense, as well as the latest gnomAD v4.0 MAF data.

BUG FIXES

- Bugfix in 'getGScores()' to properly handle resource license
  information.

[GeomxTools](/packages/GeomxTools)
----------

                        Changes in version 3.7.3                        

- small bug fix in version comparisons (PKC & Seurat)
- update SpatialExperiment vignette plot coordinates

                        Changes in version 3.7.2                        

- Seurat v5 coercion
- IPA data loading
- improved Proteogenomics error handling
- bug fix on area filtering

[gg4way](/packages/gg4way)
------

                        Changes in version 1.1.3                        

- Legend improvements

                        Changes in version 1.0.2                        

- Improved warning for when genes are not shared by both contrasts

                        Changes in version 1.0.1                        

- Minor fix to add support for custom column names to lists of
data.frames

[ggkegg](/packages/ggkegg)
------

                 Changes in version 1.1.17 (2024-04-06)                 

- stamp function

- Make some colors non-default in overlay_raw_map

[ggsc](/packages/ggsc)
----

                        Changes in version 1.1.4                        

- add background outline of umap plot (2024-04-12, Fri, #22)

                        Changes in version 1.1.3                        

- support plotting pie for spatial data (2023-12-13, Wed, #18)
- extract meta.data of Seurat object (2023-12-12, Tue, #17)

                        Changes in version 1.1.2                        

- add sc_dot() methods (2023-11-29, Wed, #15)
- update vignette to link to the online docs
(https://yulab-smu.top/ggsc)
- add README.Rmd and README.md in github repo
- optimization: retrieve embedding without FetchData (2023-11-27,
Mon,
#14)
- bug fixed for Seurat object (2023-10-31, Tue, #12, #13)

                        Changes in version 1.1.1                        

- ignore the tissue section when image is not exist (2023-10-31, Tue,
#12)
- introduce joint to combine all features with joint.fun and speed up
calculation of kde using RcppParallel (2023-10-25, Wed, #11)

[ggspavis](/packages/ggspavis)
--------

                 Changes in version 1.9.1 (2024-03-16)                  

- major updates to several plotting functions providing improved
  flexibility
  and new features (contributions by Estella Dong)

[ggtree](/packages/ggtree)
------

                       Changes in version 3.11.2                        

- expose 'position' parameter for geom_range() (2024-04-07, Sun,
#611)

                       Changes in version 3.11.1                        

- copy rect_to_poly() from ggplot2 to make it compatible with
ggplot2 3.5.0 (2024-02-13, Tue)

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                       Changes in version 1.13.0                        

- Bioconductor 3.18 released and Bioconductor 3.19 (devel) version
bump. (2023-10-25, Wed)

[ggtreeSpace](/packages/ggtreeSpace)
-----------

                       Changes in version 0.99.4                        

- address BioC review comments (2024-04-09, Tue, #7-#10)
- add vignette (2024-03-23, Sat)
- fixed issues reported by BiocCheck (2024-03-21, Thu)
- phylospm() (2022-11-22)
- geom_tsheatmap() (2022-11-18)
- geom_treeSpace() & ggphylospm() (2022-11-15)
- ggtreeSpace() to plot phylomorphospace (2022-06-25)

[gINTomics](/packages/gINTomics)
---------

                        Changes in version 1.0.0                        

- New package gINTomics, for Multi Omics data integration and
  visualization

[glmGamPoi](/packages/glmGamPoi)
---------

                  Changes in version 1.15 (2023-11-04)                  

- Fix typo in `quasi_gamma_poisson_shrinkage` example (thanks to
  @nlubock)

- Add `sample_fraction` argument to `loc_median_fit` (once again thanks
  @nlubock)

[glmSparseNet](/packages/glmSparseNet)
------------

                       Changes in version 1.22.0                        

Breaking changes

- Adds experiment parameter to different glm* functions.
- Changing in hashing function to use rlang::hash that reduces a
dependency.
- Deprecates hallmarks function as API has been shutdown.

Miscellaneous

- Use of native pipe instead of {magritrr}'s %>%
- Corrects styling and linter issues for better code quality and
readability
- Starts to deprecates parameters using dot.case in favor of
camelCase.
- Increases code coverage to 95%.

[GloScope](/packages/GloScope)
--------

                 Changes in version 1.1.1 (2024-05-01)                  

- Added functionality to compute divergence matrix from cell-type
  proportion vectors

[GNOSIS](/packages/GNOSIS)
------

                 Changes in version 1.99.0 (2023-09-04)                 

- Made the following significant changes
  o added functionality to select and upload cBioPortal study
  o deprecated ability to save R script with executed code

- Submitted to Bioconductor

[GOfuncR](/packages/GOfuncR)
-------

                       Changes in version 1.23.2                        

USER-LEVEL CHANGES

- Bugfix for compilation error on Mac

- Bugfix for refinement "signif" column (see
  https://github.com/sgrote/GOfuncR/pull/7)

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.29.2                        

- update buildGOmap() parameter to consistent with enricher() and
GSEA() (2024-02-06, Tue, #47)

                       Changes in version 2.29.1                        

- extend godata() to support passing a data.frame (can be output of
read.gaf() or read.blast2go()) to 'annoDb' (2023-01-16, Tue)
- deprecate 'OrgDb' and introduce new parameter 'annoDb' in godata()
- standardize the output of read.gaf() and read.blast2go()
- optimize buildGOmap()

[GrafGen](/packages/GrafGen)
-------

                 Changes in version 0.99.0 (2024-03-13)                 

- Submitted to Bioconductor

[GRaNIE](/packages/GRaNIE)
------

              Changes in version 1.7.3-1.7.4 (2024-04-03)               

New features

- We offer an exciting new feature: Integrating Capture Hi-C or more
generally known promoter-enhancer interactions directly into the
GRaNIE framework to complement / guide the peak-gene search. For
more information, see the Package vignette and the R help for
addConnections_peak_gene().
- Added support for the new JASPAR2024 package and TF motives.

Bug fixes

- small bug fixes

                 Changes in version 1.7.2 (2023-12-08)                  

Improvements

- plotPCA_all() now also stores all screeplot and PCA results in the
object within GRN@stats$PCA

                 Changes in version 1.7.1 (2023-10-26)                  

- version jump due to new Bioconductor development cycle

Bug fixes

- fixed an accidentally recently introduced bug that caused an error
in addTFBS when using the JASPAR database
- made the code more error prone related to AnnotationHub and caching
annotation data in cases when the cache directory is corrupted or
deleted

[graphite](/packages/graphite)
--------

                 Changes in version 1.49.1 (2024-04-29)                 

- Updated all pathway data.

[GreyListChIP](/packages/GreyListChIP)
------------

                       Changes in version 1.35.1                        

- Added yieldSize argument to countReads method.

[GSEABase](/packages/GSEABase)
--------

                        Changes in version 1.66                         

SIGNIFICANT USER-VISIBLE CHANGES

- faster implementation of `incidence()` for large gene sets /
  collections
  (https://github.com/Bioconductor/GSEABase/issues/9)

[GSVA](/packages/GSVA)
----

                        Changes in version 1.52                         

USER VISIBLE CHANGES

- Moved old API from deprecated to defunct; see vignette and help pages
  for examples on how to use the new API.

- Documentation fixes.

- Gene/feature filtering is not based anymore on floating point
  arithmetic, but on comparing minimum and maximum values. In the case
  of expression data stored in dgCMatrix objects, this criterion
  applies to non-zero values only.

- Parameter objects have nicer show methods.

- Resulting enrichment score matrices are now always dense.

- Added a first version of the support to SpatialExperiment objects. At
  the moment, only GSVA scores are calculated without using spatial
  coordinates information.

- When expression data is stored in an input sparse matrix of class
  dgCMatrix, now it is handled as such a sparse matrix as much as
  possible to reduce memory consumption.

- Added geneSets() and geneSetSizes() methods that allow one to
  retrieve, respectively, the filtered collection of gene sets and
  their sizes from either the parameter object or the resulting object
  output by the gsva() function. This allows the user to more easily
  give this information to analysis pipelines that exploit it, such as
  limma-trend; see the vignette for an example.

- Added functions readGMT() and deduplicateGeneSets() to read GMT
  files, handling the case when the file contains gene sets with
  duplicated names; see the corresponding help pages for more
  information.

- The underlying code of the ssGSEA method has been optimized after
  discussion on https://github.com/rcastelo/GSVA/issues/71 and now it
  runs one order of magnitude faster and consumes one order of
  magnitude less memory.

BUG FIXES

- Bugfix on accessing the assay names of a SingleCellExperiment object
  from the 'gsva()' function.

- Bugfix on a rare combination of input parameter for gsvaParam().

[hca](/packages/hca)
---

                        Changes in version 1.12                         

Bug fixes

- (v. 1.11.1) unlink cache on Bioconductor build system once every 2
weeks to mitigate cache corruption.

- (v. 1.11.2) manifest() requires modification to accommodate changes
introduced by HCA

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.32.0                        

NEW FEATURES

- Some light refactoring of the HDF5 dump management utilities:
  - All the settings controlled by the get/setHDF5Dump*() functions are
  now formally treated as global options (i.e. they're stored in the
  global .Options vector). The benefit is that the settings will always
  get passed to the workers in the context of parallel evaluation, even
  when using a parallel back-end like BiocParallel::SnowParam.
  In other words, all the workers are now guaranteed to use the same
  settings as the main R process.
  - In addition, getHDF5DumpFile() was further modified to make sure
  that
  it will generate unique "automatique dump files" across workers.

SIGNIFICANT USER-VISIBLE CHANGES

- Change 'with.dimnames' default to TRUE (was FALSE) in
  writeHDF5Array().

BUG FIXES

- Make sure that chunkdim(x) on a TENxRealizationSink,
  CSC_H5SparseMatrixSeed, or CSR_H5SparseMatrixSeed object 'x'
  **always** returns dimensions that are at most dim(x), even
  when 'x' has 0 rows and/or columns.

[hdxmsqc](/packages/hdxmsqc)
-------

                       Changes in version 0.99.3                        

hdxmsqc 0.99.3

- updated documentation in response to bioconductor review

hdxmsqc 0.99.1

- initial version of hdxmsqc package

[hermes](/packages/hermes)
------

                        Changes in version 1.7.1                        

Enhancements

- New plotting function draw_heatmap to produce heatmaps of
(normalized) counts.

Miscellaneous

- The utility function df_cols_to_factor now also converts existing
factors to having explicit missing levels.
- Version bump on forcats dependency.
- Removed ggplot2 deprecation warning ..count...

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.40.0                        

- new option 'all' in `hlaUniqueAllele()`

                       Changes in version 1.38.3                        

- fix compiler warnings: -Wformat & -Wformat-security

- new 'use.matching=TRUE' in `hlaPredMerge()`; to set
  'use.matching=FALSE'
  for backward compatibility

- 'ret.postprob=FALSE' by default in `hlaPredMerge()`

                       Changes in version 1.38.1                        

- fix a compiler warning of "unused-but-set-variable" on Apple ARM
  chips

- fix the failure of package loading on Apple ARM chips in the R
  console

[HicAggR](/packages/HicAggR)
-------

                       Changes in version 0.99.7                        

- caching in vignette

                       Changes in version 0.99.6                        

- default set column names in OrientateMatrix are now removed, this
is
for ggAPA
- ggAPA: "rf" mode has it's own customized axes labels now
- added introduction to the package's vignette

                       Changes in version 0.99.5                        

- added option to remove duplicated submatrices in SearchPairs

                       Changes in version 0.99.4                        

- removed CITATION file

                       Changes in version 0.99.3                        

- NAMESPACE was generated with roxygen2 to define exportable
functions.
- CompareToBackground: to correct the skewedness of o/e values
towards
long distances, computation of z.scores is now calculated using
residuals from a polynomial model that fits the background couples
(log(counts)~distance).
- SearchPairs: added option to remove self interacting bins.
- all internal functions are now in utilities.R.
- Docs were reviewed.

                       Changes in version 0.99.2                        

- Implemented import of corrected matrices for data in .hic and
cool/mcool format.
- Implemented import of O/E matrix for data in .hic format.
- Implemented import of raw data in .h5 format.
- Added GetInfo function to get info on a hic data (.hic,
cool/mcool/h5 formats).
- Removed dependency to BSDA::z.test in CompareToBackground.
- Removed dependency to InteractionSet and added it as package to
import in NAMESPACE to remove the PackageStartUpMessages.
- Removed chatty package start up message and replaced it with nicer
message.
- Some BiocCheck NOTES were also takedn into consideration: changing
sapply to vapply etc.

                       Changes in version 0.99.1                        

- Corrected with Bioconductor's reviews
- Added PrepareMtxList, ImportLoops, plotMultiAPA &
CompareToBackground functions
- Corrected bugs on seqlevels consistancy and name column for GRanges
objects
- Encapsulated small and internal functions in utilities.R
- ExtractSubMatrix has option to remove duplicated submatrices
- Corrected quantilization operations in Aggregation
- Corrected over all code with suggestions from BiocCheck

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[hoodscanR](/packages/hoodscanR)
---------

                        Changes in version 1.1.2                        

Fix bugs in test with R4.4

                        Changes in version 1.1.1                        

Small bug fix:

- when data have cell_id within the colData.

                        Changes in version 1.1.0                        

Published in Bioconductor.

[HybridExpress](/packages/HybridExpress)
-------------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[hypeR](/packages/hypeR)
-----

                        Changes in version 2.0.1                        

- Dots plots are now explicitly sized with size_by=c("genesets",
"significance", "none")

[iCOBRA](/packages/iCOBRA)
------

                       Changes in version 1.31.2                        

- Added 'Close app' button to close the shiny app

[IgGeneUsage](/packages/IgGeneUsage)
-----------

                        Changes in version 1.17                         

- model expanded for designs with biological replicates

- model expanded for paired sample designs (also with replicates)

[igvShiny](/packages/igvShiny)
--------

               Changes in version 2024-04-23 (2024-04-23)               

- change file links from igv-data.systemsbiology.net to
gladki.pl/igvr

               Changes in version 2024-03-16 (2024-03-16)               

- add shinytest2 for igvShinyDemo-GFF3.R

               Changes in version 2024-03-14 (2024-03-14)               

- fix issues with GFF3 data
- make igvShiny demo app for GFF3 working
- update trackName of GFF3 (from URL)
- udpate path to local GFF3

               Changes in version 2024-02-28 (2024-02-28)               

- add pkgdown content

               Changes in version 2024-02-16 (2024-02-16)               

- fix bug in function loadBamTrackFromLocalData
- improvge way of loading BAM files - show mismatches

               Changes in version 2024-02-09 (2024-02-09)               

- fix some Bioconductor NOTEs

               Changes in version 2024-02-05 (2024-02-05)               

- fix some Bioconductor NOTEs

               Changes in version 2024-02-04 (2024-02-04)               

- make the first Bioconductor release

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.9.2 (2023-12-14)                  

- maintainer change

                 Changes in version 1.9.1 (2023-11-20)                  

- Switched from aes_ to .data

- Initialise SPE and SCE objects with colnames

- Additional test within binAcrossPixels

[infercnv](/packages/infercnv)
--------

                 Changes in version 1.18.1 (2023-12-01)                 

- Update to work with Seurat v5 object changes.

- New default setting for leiden_resolution now set to "auto" which
  will very roughly set a resolution value that scales with the number
  of cells to avoid over splitting subclusters when the number of cells
  increases.

- Fix for plotting with dynamic_resolution enabled when the png size
  limit of cairo is hit to not have parts of the plot that try to be
  plotted outside the limit.

- Add useRaster option to plot_per_group method and transfer it
  plot_cnv calls.

- Transfer cluster_by_groups option to HMM when running in samples mode
  to allow running the HMM on "all_observations".


[IRanges](/packages/IRanges)
-------

                       Changes in version 2.38.0                        

NEW FEATURES

- Add terminators(), same as promoters() but for terminator regions.

[iSEE](/packages/iSEE)
----

                       Changes in version 2.15.1                        

- Add button 'Draft out a tour' to navigation bar.
- Add button 'About this data set' to navigation bar.

[iSEEfier](/packages/iSEEfier)
--------

                       Changes in version 0.99.2                        

- Addressing the points raised in the Bioc review
- Better checks of the arguments - more compact and robust
- Explicitly suggesting the instructions to install potentially
missing packages for iSEEnrich()

                       Changes in version 0.99.1                        

- Ready for Bioconductor review

                       Changes in version 0.99.0                        

- Ready for Bioconductor submission!

                        Changes in version 0.3.0                        

- Main functions are equipped with extra parameters determining their
behavior
- Unit test suite is fully in!

                        Changes in version 0.2.0                        

- Functions renamed to the final version, in a more matching and
descriptive manner

                        Changes in version 0.1.0                        

- Initial concept of the package!

[iSEEindex](/packages/iSEEindex)
---------

                        Changes in version 1.1.1                        

- Add possibility to customise the title of the app.

[iSEEu](/packages/iSEEu)
-----

                       Changes in version 1.15.1                        

- Expanded the content of the vignette, to have AggregatedDotPlot()
and the MarkdownBoard() panels highlighted

[ISLET](/packages/ISLET)
-----

                        Changes in version 1.5.1                        

- personalized deconvolution methods added.

[IsoBayes](/packages/IsoBayes)
--------

                        Changes in version 1.0.1                        

- load_data inputs SE objects only

- Updated from C++11 to C++17

[isomiRs](/packages/isomiRs)
-------

                       Changes in version 1.30.1                        

- Add color to PCA plot by `isoTop` function

- Fix depreciated dplyr code

- Many other fixes

[kebabs](/packages/kebabs)
------

                       Changes in version 1.37.2                        

- changed e-mail address of maintainer

- updated README.md and formatting of package vignette

- updated references in documentation

- fixed registration of API calls to C/C++ routines

                       Changes in version 1.37.1                        

- minor changes to C++ source code to avoid warnings related to
  printf-style
  format strings

                       Changes in version 1.37.0                        

- new branch for Bioconductor 3.19 devel

[knowYourCG](/packages/knowYourCG)
----------

                       Changes in version 0.99.0                        

- First submission of knowYourCG package to Bioconductor

[lemur](/packages/lemur)
-----

                         Changes in version 1.1                         

- Make predict function faster and less memory intensive for subset
fits.
- Speed-up internal function get_groups
- Gracefully handle duplicated column names in colData(fit)
- Give better error message in test_de if cond(..) is used for a fit
that was not specified with a design formula (thanks
@MaximilianNuber for reporting)

[limma](/packages/limma)
-----

                       Changes in version 3.60.0                        

- 
  The default settings for the `small.n` and `min.span` arguments
  of chooseLowessSpan() have been increased, increasing the
  chosen span value.

- 
  New argument `adaptive.span` for voom(). If `TRUE`, then
  chooseLowessSpan() is used to select an optimal `span` value
  depending on the number of genes, same as is done for vooma(),
  voomaByGroup() and voomaLmFit().

- 
  The plot information saved by voom() when `save.plot=TRUE` now
  includes `pch` and `cex` plotting character settings.

- 
  The default `span` set by vooma() is increased slightly to the
  value given by `chooseLowessSpan(ngenes, small.n=50,
  min.span=0.3, power=1/3)`. A new argument `legacy.span` had
  been added to optionally restore the old default settings for
  `span` for users who want backward compatibility.

- 
  vooma() argument `covariate` renamed to `predictor`.

- 
  vooma() has been revised in several other ways to match the
  behavior of voom() more closely. The order of arguments has
  been adjusted to match voom() and lmFit(). A new argument
  `save.plot` has been added similar to the same argument for
  voom(). Changes to have also been made to the title and xlab
  for the variance trend plot when `predictor` is non-NULL.

- 
  New function voomaLmFit() with the same functionality as
  vooma() but which automates the estimation of sample weights
  and intrablock correlation from vooma(). It produces an
  MArrayLM fit object instead of an EList object. It is the
  analogous to edgeR::voomLmFit() but for continuous
  microarray-like data.

- 
  Edits to the help pages for decideTests(), voom(),
  voomWithQualityWeights(), vooma() and "11RNAseq". Example added
  to the vooma() help page.

- 
  Add RNA-seq to package description.

- 
  Add Charity Law, Goknur Giner and Mengbo Li to author list in
  package description.

[limpca](/packages/limpca)
------

                       Changes in version 0.0.99                        

- submitted to bioconductor
- added the argument lmpDataList to pcaBySvd(), plotLine(),
plotScatter(), plotScatterM(), plotMeans(), plotDesign()
- added the function data2LmpDataList()

[lute](/packages/lute)
----

                       Changes in version 0.99.0                        

- Add testthat unit tests.

[MACSr](/packages/MACSr)
-----

                 Changes in version 1.11.2 (2023-11-20)                 

- Upgrade to MACS 3.0.0.

[MANOR](/packages/MANOR)
-----

                       Changes in version 1.75.2                        

- Replace calls to 'class' by calls to 'inherits'.

                       Changes in version 1.75.1                        

- Fix minor issues in .Rd files.

[MAPFX](/packages/MAPFX)
-----

                 Changes in version 0.99.7 (2024-04-30)                 

- Fixed issues on WinOS to pass.

- Updated the README.md file.

                 Changes in version 0.99.1 (2024-04-23)                 

User Visible Changes

- Addressed reviewer's comments, including removing repetition and
  adding unit tests.

- Ran examples with 2 cores.

                 Changes in version 0.99.0 (2024-03-18)                 

User Visible Changes

- Created the MAPFX package.

[mastR](/packages/mastR)
-----

                        Changes in version 1.3.6                        

- Add bioRxiv citation for mastR package.

                        Changes in version 1.3.5                        

- Update function remove_bg_exp() using Gaussian distribution
percentiles to replace min-max scaling as relative exppression
within each sample.

                        Changes in version 1.3.4                        

- Update function sig_gseaplot() to allow more custom arguments for
enrichplot::gseaplot2().

                        Changes in version 1.3.3                        

- Fix names mismatch problem when passing user-defined contrast
matrix
to DE analysis functions.

                        Changes in version 1.3.2                        

- Update function voom_fit_treat() to allow pass user-defined
contrast
matrix.

                        Changes in version 1.3.1                        

- Remove Matrix version bound.

                        Changes in version 1.2.1                        

- Specify Matrix (<= 1.6.1.1) to avoid conflicts between SeuratObject
and Matrix 1.6-2. Will fix to update to the latest version in
BiocVersion 3.19.

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.11.7 (2024-04-25)                 

- add dplyr:: in front of pull function in vignette to avoid errors

                 Changes in version 1.11.6 (2024-04-23)                 

- add parameter ... to function updateSE that will be passed
  to SummarizedExperiment::assay within updateSE

                 Changes in version 1.11.5 (2024-04-12)                 

- set parameter multiplyByNormalizationValue in normalizeAssay
  To TRUE in shinyQC

                 Changes in version 1.11.4 (2024-04-02)                 

- add option to display size in dimensionReductionPlot

                 Changes in version 1.11.3 (2024-03-20)                 

- allow quantile normalisation (method = "quantile") for assays with
  columns
  that contain only NA values

                 Changes in version 1.11.2 (2024-03-18)                 

- fix warning in normalizeAssay

                 Changes in version 1.11.1 (2024-03-15)                 

- add parameter multiplyByNormalizationValue in normalizeAssay

[matter](/packages/matter)
------

                       Changes in version 2.5.22                        

BUG FIXES

- Fix y-flipped raster images on non-macOS platforms

                       Changes in version 2.5.21                        

NEW FEATURES

- Add new vizi mark 'image' for plotting pre-rastered images

- Add arguments 'rasterImages' and 'rasterParams' to 'plot_image()'

                       Changes in version 2.5.20                        

SIGNIFICANT USER-VISIBLE CHANGES

- Update 'trans2d()' and 'warp2_trans()' to support arrays

                       Changes in version 2.5.19                        

NEW FEATURES

- Add new vizi mark 'rules' for reference lines

                       Changes in version 2.5.18                        

NEW FEATURES

- Add support for optional 'plotly' graphics output

- Add 'vizi_engine()' for setting the plotting engine

- Add support for 3D images to 'plot_image()'

SIGNIFICANT USER-VISIBLE CHANGES

- Export 'parse_formula()' for developer use

BUG FIXES

- Fix 'plotly' error from non-gridded voxels

                       Changes in version 2.5.17                        

SIGNIFICANT USER-VISIBLE CHANGES

- Allow passing subplots to 'as_facets()' via '...'

BUG FIXES

- Fix 'rowsweep()'/'colsweep()' behavior with NA groups

- Fix 'cv_do()' and 'mi_learn()' behavior with NA labels

                       Changes in version 2.5.16                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add PLS 1-component regressions to 'opls' output

- Add 'fitted()' method for 'opls'

- Update 'predict()' method for 'opls'

- Add multiple instance learning support to 'cv_do()'

BUG FIXES

- Fix 'opls_nipals()' usage with 'mi_learn()'

- Fix error setting bags to negative class in 'mi_learn()'

                       Changes in version 2.5.15                        

BUG FIXES

- Use names "MacroRecall"/"MacroPrecision" in 'cv_do()'

                       Changes in version 2.5.14                        

NEW FEATURES

- Add function 'rocscore()' for calculating ROC AUC

BUG FIXES

- Fix bug in 'mergepeaks()' merging too aggressively

- Fix NAs in probability in 'nscentroids()'

- Fix formatting for 'size_bytes()' when vmem is 0

- Fix bug in 'cv_do()' not processing test set

                       Changes in version 2.5.13                        

BUG FIXES

- Fix error in 'estnoise_filt()' when all peaks are noise

                       Changes in version 2.5.12                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add 'free' argument to 'plot_signal()' and 'plot_image()'

- Add 'style' argument to 'set_par()'

- Return 'fitted.values' from each fold in 'cv_do()'

- Add 'pos' argument to 'mi_learn()' to specify positive class

- Export 'chunkify()' and 'chunk_writer()' utilities

- Export 'size_bytes()' constructor

BUG FIXES

- Improve error messages for invalid plotting options

- Rescale 'alpha' channel for images when 'enhance=TRUE'

- Return 'probability' component for 'sgmixn()'

- Prevent mark 'boxplot' from plotting extra axes

- Fix 'mi_learn()' failing for missing values in response

- Fix 'plot_image()' failing for constant opacity if 'scale=TRUE'

                       Changes in version 2.5.11                        

NEW FEATURES

- Add RNG utility functions 'RNGStreams()',
  'getRNGStream()', and 'setRNGStream()'

- Add 'seeds' argument for parallel-safe RNG
  for 'chunk_lapply()', etc.

- Add 'type' argument to 'drle' constructor
  to allow pure-RLE or sequential encoding

- Add function 'sgmixn()' for fitting multiple
  spatial Gaussian mixture models in parallel

- Add vizi marks 'intervals' and 'boxplot'

- Add 'predscore()' for scoring predictions

- Add 'cv_do()' for performing cross-validation

SIGNIFICANT USER-VISIBLE CHANGES

- Change 'fitted()' argument 'type' to use "response"
  (instead of "probability") for 'nscentroids'

- Print methods for models (e.g., 'pls', etc.) now
  truncate output to 'getOption("matter.show.head.n")'

- Add 'jitter' transformation for vizi mark 'points', etc.

- Export 'avg()' utility function

BUG FIXES

- Fix bug in 'nscentroids()' predictions

- Fix error in 'nscentroids()' for one-class models

- Fix error in 'sgmix()' caused by singleton classes

- Pass 'weights' argument to 'distfun()'
  in 'fastmap()' and 'nscentroids()'

- Fix error in 'plot_image()' for NA-only images

                       Changes in version 2.5.10                        

NEW FEATURES

- Add function 'peakheights()'

- Add nearest shrunken centroids ('nscentroids()')

- Add spatial Gaussian mixture model ('sgmix()')

- Add colocalization coefficients ('coscore()')

- Add multiple instance learning ('mi_learn()')

SIGNIFICANT USER-VISIBLE CHANGES

- Rename 'estres()' parameter 'tol.ref' to 'ref'

- Export 'array_ind()' and 'linear_ind()'

- Add 'nchunks' argument to 'prcomp()' and 'pls()'

- Vectorize 'predict.pls()' over 'k' argument

- Set corresponding attributes to NULL when center/scale
  are FALSE in 'rowscale()' and 'colscale()'

BUG FIXES

- Fix 'tolerance()' returning invalid value for 'tolerance()<-'

- Fix error in 'peakareas()' when peak boundaries are NULL

- Remove duplicate peaks in 'findpeaks_cwt()'

- Fix potential infinite loop in 'binpeaks()'

- Fix potential incorrect vector lengths in 'simspec()'

                        Changes in version 2.5.9                        

BUG FIXES

- When writing to a file in chunk apply functions, preserve
  order in file even if chunks are processed out-of-order

                        Changes in version 2.5.8                        

NEW FEATURES

- Added scaling functions 'rescale_rms()', 'rescale_sum()',
  'rescale_ref()', 'rescale_range()', and 'rescale_iqr()'

                        Changes in version 2.5.7                        

SIGNIFICANT USER-VISIBLE CHANGES

- Moved 'keys' and 'keys<-' generics to Cardinal

- Removed 'chunksize' and 'chunksize<-' generics

BUG FIXES

- Fix error with NAs in 'add_alpha()'

- Fix flipped x/y coordinates in 'inpoly()'

                        Changes in version 2.5.6                        

BUG FIXES

- Fix error in apply functions when 'BPPARAM=NULL'

- Remove zero-height peaks in 'plot_signal()'

                        Changes in version 2.5.5                        

NEW FEATURES

- Added 'plot_signal()' and 'plot_image()'

- Added new mark 'vizi_text'

SIGNIFICANT USER-VISIBLE CHANGES

- Add 'vizi' style 'classic' for transparent background

- Implement 'vm_used()' for 'sparse_arr' objects

- Implement 'combine()' for 'vizi_plot' objects

- Rename 'plot_facets()' to 'as_facets()'

BUG FIXES

- Fixed RStudio plotting bugs due to 'par(bg="transparent")'

- Fixed 'matter_list' subsetting not subsetting type

- Fixed 'vizi_pixels' issues with y-axis

- Fixed plotting functions not passing along parameters

- Handle mixed int/double types in 'ltob()' and 'lttb()'

- Implement 'vm_used()' for 'sparse_arr'

                        Changes in version 2.5.4                        

BUG FIXES

- Improved 'mergepeaks()' efficiency

                        Changes in version 2.5.3                        

BUG FIXES

- Fixed stack overflow in 'sparse_mat' subsetting

- Fixed 'mergepeaks()' failing for missing peaks

                        Changes in version 2.5.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- For 'estres()', allow 'tol=NA' or 'tol=Inf'

                        Changes in version 2.5.1                        

BUG FIXES

- Fixed bug in 'simspec()' not respecting 'units'

- Enabled calling 'vizi()' with no arguments

[metabinR](/packages/metabinR)
--------

                 Changes in version 1.5.1 (2024-04-07)                  

- Preparing for next Bioconductor Release.

                 Changes in version 1.5.0 (2023-10-29)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_18
  branch.

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                         Changes in version 1.7                         

Changes in 1.7.5

- Add parameter addOriginalQueryIndex to matchSpectra() that allows
to
add an additional spectra variable to the query Spectra with the
index in the original object (issue #114).

Changes in 1.7.4

- Import setBackend() generic from ProtGenerics.

Changes in 1.7.3

- Add SingleMatchParam for filterMatches to allow selection of (at
most) a single match to a target element for each query element.
- Add new methods queryVariables and targetVariables to extract the
names of variables (columns) of query and target.

Changes in 1.7.2

- Update the Spectra objects within the package to the new versions.

Changes in 1.7.1

- Add examples and a section to the vignette explaining the use of
createStandardMixes.

[MetaboCoreUtils](/packages/MetaboCoreUtils)
---------------

                        Changes in version 1.11                         

MetaboCoreUtils 1.11.3

- Add examples on isotopes (including deuterium) can be used with
calculateMass (issue #81)

MetaboCoreUtils 1.11.2

- Add functions to compute quality check of the data (issue #77

MetaboCoreUtils 1.11.1

- Add functions to enable linear model-based adjustment of (LC-MS
derived) abundance matrices (issue #75).

[metaseqR2](/packages/metaseqR2)
---------

                 Changes in version 1.15.1 (2024-03-07)                 

NEW FEATURES

- None.

BUG FIXES

- Fixed supported Ensembl annotation versions

[methodical](/packages/methodical)
----------

                       Changes in version 0.99.0                        

- Submission to Bioconductor

[methyLImp2](/packages/methyLImp2)
----------

                 Changes in version 0.99.8 (2024-03-16)                 

- Fixing parallelization set-up.

                 Changes in version 0.99.7 (2024-03-14)                 

- Fixing "exportglobals = FALSE" in the vignette.

                 Changes in version 0.99.6 (2024-03-14)                 

- Fixing "BPPARAM" and "assay" usage.

                 Changes in version 0.99.5 (2024-03-05)                 

- Add "BPPARAM" argument for BiocParallel parallelization.

- Add "which_assay" and "overwrite_res" arguments for better
  manupulation of SummarizedExperiment object.

- Change license.

                 Changes in version 0.99.4 (2024-02-29)                 

- Switch to ChAMPdata package for annotation instead of the internal
  one to save space.

- Add unit tests.

- Put splitting the dataset by chromosomes in a separate function.

- Fix Bioconductor review comments.

                 Changes in version 0.99.3 (2023-07-01)                 

- Fix notes.

                 Changes in version 0.99.2 (2023-06-12)                 

- Change vignette example data to SummarizedExperiment.

                 Changes in version 0.99.1 (2023-06-01)                 

- Change of possible input: either numeric matrix or
  SummarizedExperiment;

- Adding independent imputation for different samples groups.

                 Changes in version 0.99.0 (2023-05-08)                 

- First commit.

[MGnifyR](/packages/MGnifyR)
-------

                       Changes in version 0.99.23                       

Date: 2024-03-04

- getReturn fix: failed constructing MAE if samples in experiments did
  not match

                       Changes in version 0.99.20                       

Date: 2024-02-26

- searchAnalysis returns now a named vector where names are accession
  IDs that was fed as input

                       Changes in version 0.99.19                       

Date: 2024-02-15

- Fix deprecated mgnify_client function

                       Changes in version 0.99.18                       

Date: 2024-02-12

- Last modifications for Biocondutor submission

                       Changes in version 0.99.17                       

- Added getData function for fetching raw data from the database

                       Changes in version 0.99.0                        

- Support for TreeSummarizedExperiment and MultiAssayExperiment

- Submitted to Bioconductor

[mia](/packages/mia)
---

                        Changes in version 1.11                         

- loadFromMetaphlan: support strain rank

- agglomerateByRank: agglomerate tree fix

- Replace taxonomyTree and addTaxonomyTree with getHierarchyTree and
  addHierarchyTree

- splitOn: update rowTree fix

- perSampleDominantFeatures: add new arguments (n, other.name,
  complete)

- loadFromMetaphlan: support "taxonomy" column for specifying taxonomy

- cluster: Overwrite old results instead of failing

- getPrevalence: bugfix, if assay contains NA values, it does not end
  up to NA anymore.

- getExperimentCrossCorrelation fix: enable using of sampleMap in MAE.

- Implemented the setTaxonomyRanks function to specify which ranks are
  recognized as taxonomy ranks.

- Rename cluster to addCluster

- rename importers loadFromBiom, loadFromQIIME2, readQZA,
  loadFromMothur, loadFromMetaphlan, loadFromHumann

- fix typo in loadFromBiom definition (deprecate file)

- deprecate subsetSamples, subsetFeatures and subsetTaxa

- deprecate plotNMDS after moving it to miaViz package

- rename estimateDivergence to addDivergence

- Add details to documentation of function agglomerateByPrevalence

[miaViz](/packages/miaViz)
------

                        Changes in version 1.11                         

- replace addTaxonomyTree with addHierarchyTree after renaming in mia
  package

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.15.1                        

- fix a bug of mp_import_qiime when sample metadata has - character.
(2024-04-12, Fri)
- update mp_plot_diff_cladogram with tidytree and treeio.
(2024-03-26,
Tue)

                       Changes in version 1.15.0                        

- Bioconductor 3.18 released and Bioconductor 3.19 (devel) version
bump. (2023-10-25, Wed)

[midasHLA](/packages/midasHLA)
--------

                       Changes in version 1.11.1                        

- Fix bug on hlaToVariable() function that variable would erroneously
be named NA.

[miloR](/packages/miloR)
-----

                 Changes in version 2.0.1 (2023-11-09)                  

- Introduce NB-GLMM into Milo 2.0 for random effect variables and
modelling dependencies between observations
- Diagnostic function for checking model separation for experimental
variables, i.e. splitting zero from non-zero counts perfectly
- Vignette describing basic usage of GLMM functions in testNhoods

[MIRit](/packages/MIRit)
-----

                       Changes in version 0.99.13                       

A minor fix was made to fix the undefined variables note during R CMD
check.

                       Changes in version 0.99.12                       

This version includes several improvements, including a completely
revised vignette where all chunks are evaluated, minor tweaks to
default
values for differential expression analysis, and some bug fixes to
the
error bars in the plotDE() function. Other issues, such as artifacts
in
show methods, lacks in documentation, and dependence in DESCRIPTION,
have been addressed too.

                       Changes in version 0.99.11                       

This new version introduces the possibility of limiting validated
targets retrieval from miRTarBase to only those interactions
supported
by extensive experimental evidence. Moreover, minor fixes were made
to
the vignette, documentation, and internal data. Finally, some
examples
have been redefined to reduce checking time.

                       Changes in version 0.99.10                       

This patch introduces a minor fix for one unit test.

                       Changes in version 0.99.9                        

With this version, significant p-values originating from functional
enrichment analyses now include extreme values. Further, examples
have
been shortened and unit tests now use smaller datasets.

                       Changes in version 0.99.8                        

This patch further reduces R CMD check time by limiting unnecessary
examples and by using smaller datasets for unit tests.

                       Changes in version 0.99.7                        

This version introduces parallel computing capabilities for the
mirnaIntegration() function. This is particularly useful for
Boschloo's
exact test, whose execution is now faster. Moreover, the test suite
has
been redefined to reduce running times.

                       Changes in version 0.99.6                        

The testing suite has been redefined to allow different results for
different versions of packages employed in differential expression
analysis.

                       Changes in version 0.99.5                        

This update fixes a bug in the batchCorrection() function that
prevented
the correct use of this function with newer versions of the
MultiAssayExperiment package.

                       Changes in version 0.99.4                        

After the implementation of the IS_BIOC_BUILD_MACHINE variable to the
Single Package Builder (SPB), this version bump drives a new build to
fix errors during R CMD check on SPB.

                       Changes in version 0.99.3                        

This version fixes a bug in the topologicalAnalysis() function that
prevented the use of a functional progress bar during permutation
testing. Moreover, the example for the addDifferentialExpression()
function has been updated to reduce its running time.

                       Changes in version 0.99.2                        

MIRit now allows to filter the pathways used for topological analysis
based on the number of nodes.

                       Changes in version 0.99.1                        

Functional enrichment analyses and TAIPA now use cached databases to
reduce running times.

                       Changes in version 0.99.0                        

Initial version for Bioconductor submission.

[miRSM](/packages/miRSM)
-----

                    Changes in version 1.99.7-1.99.8                    

- Improve description <2024-04-05, Fri>

                    Changes in version 1.99.0-1.99.6                    

- Add new methods for identifying miRNA sponge modules <2024-02-04,
  Sun>

                       Changes in version 1.21.1                        

- Add sponge module at single-sample level and internal competition
  sponge module <2024-01-26, Fri>

[mitch](/packages/mitch)
-----

                       Changes in version 1.15.5                        

- The default number of cores has been set to 1.

- An issue with contrasts being identical has been resolved.

                       Changes in version 1.15.4                        

- A bug with tables in the HTML reports has been fixed.

- Minor improvements to the gmt_import() function to prevent empty
  sets.

                       Changes in version 1.15.2                        

- A new vignette has been included to demonstrate how to conduct
  enrichment analysis of profiles
  generated with HM450K and EPIC arrays.

                       Changes in version 1.15.1                        

- A change of the gene aggregation algorithm to allow enrichment
  analysis of infinium methylation analysis.

[mobileRNA](/packages/mobileRNA)
---------

                       Changes in version 0.99.23                       

- Improved RNAfeatures function

                       Changes in version 0.99.22                       

- Fixed issues in RNAimport function

                       Changes in version 0.99.20                       

- Build redo

                       Changes in version 0.99.19                       

- Missing connective in RNAfeatures when using repeats variable.
- Added additional check for RNAsequences methods.

                       Changes in version 0.99.17                       

- Corrected ORCID references for authors

- Amended RNAmobile() to ensure removal of non-zero values
- Improved look of PCA plot
- Corrected documentation issue and code disparage in plotHeatmap(),
and improved styling.
- Included calculation of consensus sequence determination option for
RNAsequence()

                       Changes in version 0.99.15                       

- Alteration of vignette and README
- Included clean FASTQ files as example data sets
- Updated R data objects
- Updated citation & news files
- Addition of new function called mapRNA()
- Deletion of RNAloci() and RNAmean() functions
- Additional parameters to RNAimport() to support mRNA data
importation
- Removed parallel computation in RNAmergeGenomes()
- Improved documentation of functions and removed inconsistencies.

                       Changes in version 0.99.14                       

Previous changes

- Added a NEWS.md file to track changes to the package.
- RNAconsensus() changed to RNAdicercall()
- Alterations to RNAdicercall() algorithm, including tie options,
altered default tidy method.
- RNAdicercall() introduced new column "DicerCount" and improved the
functionality, specifically the exclude parameter.
- RNAmobile() introduced new parameter, "threshold"
- Improved RNAsequences() selection algorithm to consider a threshold
value, and handling ties.
- Improved error calling on functions
- Added RNAdf2e() function
- Improvements to plotSamplePCA() and plotHeatmap()
- Amended RNAmergeAnnotations function to meet requirement
- Removed unnecessary man files for GFF and FASTA files on remote
repo
- Fixed bug in RNAdistribution() plot, when sample specific
- Broadened use of RNAattributes() function.
- Updated vignette
- Updated RNAdicercall() to allow any dicer-classification (not
constricted to 20-24)
- Converted cat() to message() for user information from functions
- Amended example data
- Amended CITATION and NEWs file
- Altered examples in RNAmergeAnnotations/Genomes functions to
prevent
examples saving into users directory
- updates inline with bioconductor requirements
- Amended RNAmergeGenomes() and RNAmergeAnnotations()
- Removed lazy loading of package data
- RNAanalysis() changed name to RNAdifferentialAnalysis()

[monaLisa](/packages/monaLisa)
--------

                        Changes in version 1.9.1                        

- adapt dumpJaspar to also work with JASPAR2024

[Moonlight2R](/packages/Moonlight2R)
-----------

                        Changes in version 1.1.2                        

Summary

- fixed issue with parallelization in URA

                        Changes in version 1.1.1                        

Summary

- added new functionality called GMA (Gene Methylation Analysis)
- added three new visualization functions in connection with GMA:
plotGMA, plotMoonlightMet, plotMetExp
- updated vignette to contain new functionalities related to GMA
- updated class testing in all functions
- resized figures

[mosdef](/packages/mosdef)
------

                       Changes in version 0.99.4                        

- Added persistent location to the script that generates the data
objects

                       Changes in version 0.99.3                        

- Working on the size and time constraints for the vignette/package
to
build and check on the BBS

                       Changes in version 0.99.2                        

- Better specification of data objects format and location in the
installed folder

                       Changes in version 0.99.1                        

- This version contains the newly implemented changes as a response
to
the Bioconductor review. In brief, this includes:

- renaming the functions (and parameters) to a more consistent
style
- reduction of the dependencies and runtime of checks/tests
- more details on the exported data objects, and on the outputs
(detailed in the vignette)
- better structure for the vignette with cross-references among
sections
- modularization of some functions to avoid repetitive code
- implementation of an API which is already framework-agnostic, to
later accommodate e.g. edgeR/limma. Mainly, this implies the
renaming of the dds to a more generic de_container, whereas the
res_de parameter stays constant.

For a full list of all changes implemented, please refer to the PR
on the mosdef repository
https://github.com/imbeimainz/mosdef/pull/11/

                       Changes in version 0.99.0                        

- Ready for submission to Bioconductor!

[MOSim](/packages/MOSim)
-----

                  Changes in version 2.0 (2023-07-21)                   

- Added scMOSim functionality

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.47.1                        

- Update documentation of motifStack function.

[motifTestR](/packages/motifTestR)
----------

                       Changes in version 0.99.0                        

- Submitted to Bioconductor


[msa](/packages/msa)
---

                       Changes in version 1.35.5                        

- major update of package help page man/msa-package.Rd

                       Changes in version 1.35.4                        

- fixes to account for move of substitution matrices from 'Biostrings'
  to 'pwalign'
  package

                       Changes in version 1.35.3                        

- changed e-mail address of maintainer

- updated README.md and formatting of package vignette

- updated references in documentation

                       Changes in version 1.35.2                        

- update of msaMakevars.win in ClustalW to avoid problems arising from
  compiling
  ClustalW with C++ 17: added -std=c++14

                       Changes in version 1.35.1                        

- update of some Makevars and Makefiles to avoid compliation issues on
  FreeBSD
  + minor adaptation in vignette

- minor fix in src/ClustalOmega/src/RClustalOmega.cpp (bug in Rprintf
  arg list)

                       Changes in version 1.35.0                        

- new branch for Bioconductor 3.19 devel

[MSA2dist](/packages/MSA2dist)
--------

                 Changes in version 1.7.5 (2024-03-28)                  

BUG FIXES

- fixed dnastring2kaks to work with local alignments
- fixed indices2kaks to work with local alignments
- fixed cds2codonaln to work with local alignments
- fixed cdsstring2codonaln to work with local alignments

NEW FEATURES

- added return.cds parameter to cds2aa for local alignments

                 Changes in version 1.7.4 (2024-03-18)                  

BUG FIXES

- fixed pal2nal

                 Changes in version 1.7.3 (2024-03-14)                  

BUG FIXES

- fixed dnastring2kaks to return correct orientation
- fixed pal2nal to use gap_pos-n_i_codons_added-1

                 Changes in version 1.7.2 (2024-03-14)                  

BUG FIXES

- fixed pal2nal to cover individual gaps longer than 1

                 Changes in version 1.7.1 (2024-03-13)                  

NEW FEATURES

- added pal2nal
- added cdsstring2codonaln
- added indices2kaks
- added KaKs Calculator 2.0 example to vignette

CHANGES

- changed dnastring2kaks to use cdsstring2codonaln function

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                        Changes in version 1.11                         

Changes in 1.11.2

- Import method generics from ProtGenerics. This requires
ProtGenerics
version 1.35.3.

Changes in 1.11.1

- Remove additional empty line at the end of exported MassBank
records
(issue #49).

[MsBackendMgf](/packages/MsBackendMgf)
------------

                        Changes in version 1.11                         

Changes in 1.11.2

- Import generic methods from ProtGenerics. Requires ProtGenerics
version 1.35.3.

Changes in 1.11.1

- Small runtime improvements in MGF importer.

Changes in 1.11.0

- Bioconductor 3.19 developmental branch.

[MsBackendMsp](/packages/MsBackendMsp)
------------

                         Changes in version 1.7                         

Changes in 1.7.3

- Strip whitespaces in values of comments/header information
- Support also name:value header pairs in addition to name: value
(issue #14).
- Add support for parallel processing for import from a single
(large)
MSP file.

Changes in 1.7.2

- Add additional checks to the format of input MSP files to ensure
proper data import.

Changes in 1.7.1

- Import method generics from ProtGenerics. Requires ProtGenerics
version 1.35.3.

[MsBackendSql](/packages/MsBackendSql)
------------

                         Changes in version 1.3                         

Changes in 1.3.5

- Improve input argument check and error message for
backendInitialize() for MsBackendOfflineSql.
- Update documentation adding () to all function names.

Changes in 1.3.4

- Ensure primary keys from the database are in the correct order for
backendInitialize().

Changes in 1.3.3

- Import method generics from ProtGenerics.

Changes in 1.3.2

- Add a dedicated setBackend method for MsBackendSql and
MsBackendOfflineSql backends (issue #17).

Changes in 1.3.1

- Add description on the use/advantages of different SQL database
systems to the vignette.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                        Changes in version 1.15                         

MsCoreUtils 1.15.7

- Add common_path() function.

MsCoreUtils 1.15.6

- Bump version to force new package build on Bioconductor servers.

MsCoreUtils 1.15.5

- Add function force_sorted() to adjust a numeric vector to ensure
increasing/sorted values.

MsCoreUtils 1.15.4

- Fix partial argument match (see issue #125).

MsCoreUtils 1.15.4

- Fix documentation of ndotproduct.

MsCoreUtils 1.15.3

- Add function breaks_ppm to create a sequence of numbers with
increasing difference between elements (defined by parameter ppm).

MsCoreUtils 1.15.2

- Porting baseline estimation function (see issue 119).

MsCoreUtils 1.15.1

- Remove impute_mle2() since norm2 has been removed from CRAN (see
issue 117).

MsCoreUtils 1.15.0

- New Bioc devel version

[MsDataHub](/packages/MsDataHub)
---------

                         Changes in version 1.3                         

MsDataHub 1.3.4

- Update vignette: load Report.Derks2022.plexDIA.tsv and
benchmarkingDIA.tsv with latest QFeatures.

MsDataHub 1.3.3

- Add Report.Derks2022.plexDIA.tsv, plexDIA DIA-NN output from Derks
et al. (2022).

MsDataHub 1.3.2

- Fix path to benchmarkingDIA.tsv on zenodo.

MsDataHub 1.3.1

- Add benchmarkingDIA.tsv data, contributed by Kristina Gomoryova.

[MsExperiment](/packages/MsExperiment)
------------

                         Changes in version 1.5                         

MsExperiment 1.5.5

- Add spectraSampleIndex() function.

MsExperiment 1.5.4

- Fix missing export of filterSpectra.

MsExperiment 1.5.3

- Add filterSpectra method to allow filtering of Spectra within an
MsExperiment while keeping possibly present relationships between
samples and spectra consistent.

MsExperiment 1.5.2

- Add support to read/write sample data from/to a MsBackendSql
database (issue #39).

MsExperiment 1.5.1

- Fix subset with negative indices (issue #37.)

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.31                         

MSnbase 2.31.1

- Disable nested parallel processing for chromatogram() method.
- Fix Rd notes.

MSnbase 2.31.0

- New Bioconductor devel.

                        Changes in version 2.29                         

MSnbase 2.29.4

- Move XML to suggests.

MSnbase 2.29.3

- Remove parts of XML dependency.

MSnbase 2.29.2

- Use fragmentation functions from PSMatch.
- Mention R for Mass Spectrometry in start-up message.

MSnbase 2.29.1

- Check for identical ion header in readMgfData() (see issue 597).

MSnbase 2.29.0

- New devel

[msqrob2](/packages/msqrob2)
-------

                       Changes in version 1.11.2                        

- Fixed issue related to levels of a ridge variable
- Fixed issue when fitting only one random effect

                       Changes in version 1.11.1                        

- Fixed issues related to reference class changes in the models
- Fixed issue related to colData assay levels
- Refactored internal code

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.30.0                        

Bug fixes and minor improvements

- Updated CITATION information in the main vignette.

[multiGSEA](/packages/multiGSEA)
---------

                 Changes in version 1.13.1 (2024-04-26)                 

- Utilize by default non-adjusted single-omics p-values for
  aggregation into a composite p-value. Adjustment for multiple
  testing will be done on the combined multi-omics p-values.

- Inlcude legacy option for p-value combination to ensure
  reproducibility.

[multistateQTL](/packages/multistateQTL)
-------------

                 Changes in version 0.99.1 (2024-04-21)                 

- Small tweaks e.g. using seq_along()

- Changed some arguments to use camelCase instead of snake_case

- Improved vignette.

                 Changes in version 0.99.0 (2024-03-22)                 

- Submitted to Bioconductor

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.11.10                       

New features

- Can now pass local chain files for liftover (local_chain in
format_sumstats() and liftover()).

                       Changes in version 1.11.9                        

New features

- Can now control what columns are checked for missing data
(drop_na_cols in format_sumstats()). By default, SNP, effect columns
and P/N columns are checked. Set to Null to check all columns or
choose specific columns.

                       Changes in version 1.11.7                        

Bug fix

- Force no tab indexing when writing removed rows of SNPs. This
avoids
any issues where missing data causes sort errors.
- Issue fixed when sorting CHR column based on a format when CHR
column is a factor.

                       Changes in version 1.11.6                        

Bug fix

- Catch for overflow when NA's in SNP col for check_no_rs_snp() check
with imputation_ind=TRUE.

                       Changes in version 1.11.4                        

Bug fix

- Minor fix to get_genome_builds() to help with RAM & CPU usage
during
unit tests. No change in functionality for end user.

                       Changes in version 1.11.3                        

Bug fix

- For LDSC format, rename A1 and A2 as LDSC expects A1 to be the
effect column rather than A2 (the opposite to MSS's default) - see
more here. Although, this didn't seem to make any difference to
results in tests, see more here.

                       Changes in version 1.11.2                        

Bug fix

- Remove unused argument make_ordered from sort_coords()
- Issue fixed with check ldsc format wehn compute_n type chosen

                       Changes in version 1.11.1                        

Bug fix

- Speed up unit test timing for bioc checks (predominately for linux
tests)

[mzR](/packages/mzR)
---

                       Changes in version 2.37.3                        

- Fix compilation on Windows/aarch64 (#290), thanks to Tomas Kalibera
  for the patch

                       Changes in version 2.37.2                        

- remove mentions of mzData in the DESCRIPTION, vignette and manual
  files (support dropped in 2.29.3)

                       Changes in version 2.37.1                        

- fix compilation on centOS 7 / R-4.3.1 reported in #286

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 3.0.0                        

- Breaking change to smoothing strategy in plot_gene(),
plot_region(),
and plot_granges() to use weighted moving mean instead of loess.
This deprecates the span argument in favour of smoothing_window
which is defaulted to 2000 bases.
- Smoothing for various plotting functions was previously
performed using loess smoothing, this performed locally weighted
linear estimation to create a smoothed line, the span argument
controlled the proportion of data used in this smoothing. This
parameter was difficult to tune because under a fixed span, the
smoothed line became flatter as the plot region grew larger.
Internally, NanoMethViz dynamically calculated a span that
changed inversely proportional to the width of the plotting
region, decreasing the span as the plot region grew. However the
calculated span was invisible to users, and it unintuitive to
users how to set a span to change the appearance of the plot.
- Changing the smoothing method to a weighted rolling mean lead to
the new smoothing_window argument which represents the window
size in bases from which data is used for smoothing around each
point. This serves the same purpose as the dynamic calculation
done previously, but is set more explicitly and should be more
intuitive for users. The default is always 2000, and can be
increased to increase smoothness and decreased to decrease
smoothness.
- Breaking change to the appearance of plot_gene() plots, previously
the isoform annotation would be restricted to only the gene of
interest. It is now changed to follow the same behaviour as
plot_region() whereby all isoforms in the region are plotted.
- Breaking change to the default plotting options for plot_region(),
plot_gene() and plot_grange() to plot heatmap by default.
- Possible breaking change to query_methy(simplify = FALSE), it will
now return a list that is the same length as the number of regions
queried, where it previously returned nothing if a particular
sequence was missing from the tabix.
- Added gene_anno argument to plot_region() and plot_granges() to
control whether gene annotation is plotted.
- Added plot_violin() function for creating violin plots for samples
over specific regions.
- Added check to remove hard-clipped reads because they may not have
matching mod strings.
- Changed gene annotation to always put label on visible isoforms,
previously labels are plotted at the center of isoform.
- Changed the order of columns when querying from ModBam to be the
same as when querying from Tabix, with readname at the end instead
of being the second column.
- Fixed memory leak in bam parsing due to out of bounds access.
- Fixed crash when CIGAR doesn't match length of SEQ.

[NanoStringNCTools](/packages/NanoStringNCTools)
-----------------

                 Changes in version 1.11.1 (2024-04-10)                 

- Remove setting font in vignette

[NanoTube](/packages/NanoTube)
--------

                        Changes in version 1.9.1                        

- Corrected a bug in positiveQC(), which caused it to calculate
different positive scale factors from normalize_pos_control(); the
second was confirmed to be correct.

[ngsReports](/packages/ngsReports)
----------

                        Changes in version 2.5.2                        

- Added `summariseOverrep()`

                        Changes in version 2.5.1                        

- Changed method for setting default theme using plotTheme

- Set factor levels for plotSeqContent for FastpDataList

- Added line & cumulative plotTypes for plotInsertSize

[nipalsMCIA](/packages/nipalsMCIA)
----------

                 Changes in version 1.2.0 (2024-04-26)                  

Major changes

Minor improvements and bugfixes

- Clarified the dataset used in the single-cell vignette and
reformatted the workflow diagram
- Removed the final log transformation in the single-cell vignette
when saving for MCIA
- Minor plotting improvements (removed redundancies in legend names,
consistency of axis text, etc.)

                 Changes in version 1.1.0 (2024-03-21)                  

Major changes

- Added BiocFileCache to replace previous method of downloading large
SC datasets.

Minor improvements and bugfixes

- Updated nipals_iter to only use var(gs) to compute the significance
of global scores.
- Added to documentation for predict_gs() warning against use with
CPCA deflation.

[nnSVG](/packages/nnSVG)
-----

                 Changes in version 1.7.1 (2024-03-08)                  

- bug fix for non-SpatialExperiment inputs

[OmaDB](/packages/OmaDB)
-----

                 Changes in version 2.19.1 (2024-04-18)                 

- fix problem with getGenomePairs function after OMA API change

[omicsViewer](/packages/omicsViewer)
-----------

                  Changes in version 1.7 (2024-03-08)                   

new functions added: 

- the enriched terms are clustered

- dose response curve function (test)

[ontoProc](/packages/ontoProc)
--------

                       Changes in version 1.25.0                        

- introducing owl operations, relying on reticulate

[ORFik](/packages/ORFik)
-----

                       Changes in version 1.23.8                        

SIGNIFICANT USER-VISIBLE CHANGES

- Implemented a Ribo-seq ORF detector now included

- Added new export functions for bigwig, covRLE

[Organism.dplyr](/packages/Organism.dplyr)
--------------

                       Changes in version 1.32.0                        

BUG FIXES

- (v. 1.31.1) avoid non-generic arguments in tbl.src_organism,
  closing
  <https://github.com/Bioconductor/Organism.dplyr/issues/19>

[orthos](/packages/orthos)
------

                        Changes in version 1.1.1                        

- Upgrade libwebp to address security concerns

[OUTRIDER](/packages/OUTRIDER)
--------

                       Changes in version 1.22.0                        

- Add Manhattenplot functionallity

- Bugfix in restricting for FDR correction and others (#54 and aa6e56a)

- Ensure correct injections of outliers (#37)

[OutSplice](/packages/OutSplice)
---------

                 Changes in version 1.3.1 (2024-03-06)                  

- Removed Repitools annoGR2DF

- Repitools removed from NAMESPACE and DESCRIPTION files

[pathlinkR](/packages/pathlinkR)
---------

                      Changes in version 0.99.380                       

Lots of updates, including:

- Updated documentation for data and functions
- Cleaner vignette
- Better support for Bioconductor classes
- Better handling of data loading

                      Changes in version 0.99.352                       

- Added value tag to data documentation

                      Changes in version 0.99.350                       

- Rewrite of pathwayPlots to be more efficient

                      Changes in version 0.99.332                       

- Changed methodology for importing functions
- Added two new functions for PPI networks: "ppiEnrichNetwork" and
"ppiExtractSubnetwork"

                      Changes in version 0.99.317                       

- Renamed a number of functions to improve consistency, and make it
clear which ones are related and part of the same "workflow"

                      Changes in version 0.99.310                       

- Big update to the vignette
- Update to the README and added in the package hex logo
- Various function tweaks and minor changes

                      Changes in version 0.99.300                       

- Renamed package to pathlinkR

                       Changes in version 0.99.1                        

- Prepared for bioconductor submission

[pathview](/packages/pathview)
--------

                       Changes in version 1.43.1                        

- add org.EcK12.eg.db to Suggests for BioC build check

- korg now include 9767 KEGG species or 1485 new species beyond 2022.

[peakPantheR](/packages/peakPantheR)
-----------

                 Changes in version 1.17.2 (2024-02-26)                 

- update unittest for `ggplot2` 3.5.0

                 Changes in version 1.17.1 (2023-11-04)                 

- `peakPantheR_quickEIC()` to plot EIC from raw files

- SVG output option for AnnotationDiagnostic summary plot

- Correct Windows-specific file locking issue when calling
  `MSnbase::spectra()` in parallel

[Pedixplorer](/packages/Pedixplorer)
-----------

                       Changes in version 0.99.0                        

- Kinship2 is renamed to Pedixplorer and hosted on Bioconductor.
- Pedigree is now a S4 object, all functions are updated to work with
the new class
- Pedigree constructor now takes a data.frame as input for the
Pedigree informations and for the special relationship. The two
data.fram are normalized before being used.
- plot.pedigree support ggplot generation, mark and label can be
added
to the plot. The plot is now generated in two steps ped_to_plotdf()
and plot_fromdf(). This allows the user to modify the plot before it
is generated.
- All documentation are now generated with Roxygen
- New function available: generate_aff_inds, generate_colors,
is_informative, min_dist_inf, normData, num_child, useful_inds
- All functions renamed to follow the snake_case convention
- All parameters renamed to follow the snake_case convention
- All test now use testthat files
- Vignettes have been updated to reflect the new changes

[pgxRpi](/packages/pgxRpi)
------

                 Changes in version 0.99.7 (2023-10-20)                 

- Add pgxFilter function to expose all available filters.

                 Changes in version 0.99.5 (2023-10-10)                 

- Updated data structure for frequency data, transitioning from a
simple list to Bioconductor containers.
- Removed sections on survival analysis and frequency clustering
analysis from the vignettes to align with the package's scope.

                 Changes in version 0.99.0 (2023-08-25)                 

- Submitted to Bioconductor

[phantasus](/packages/phantasus)
---------

                       Changes in version 1.23.3                        

- Changed a configuration scheme

- setupPhantasus() now can be used to generate a config

- Support for remote count files via HSDS

[phantasusLite](/packages/phantasusLite)
-------------

                        Changes in version 1.2.0                        

- Switch to https://alserglab.wustl.edu/hsds remote for default HSDS
  server

- Depending on rhdf5client >= 1.25.1 to support ARCHS4 v2.3 files

[PhyloProfile](/packages/PhyloProfile)
------------

                       Changes in version 1.18.0                        

- option to specify version for OrthoDB

                       Changes in version 1.16.5                        

- fixed detailed plot when name of a taxon is substring of
  another

- show supertaxon names when rotate the profile plot

                       Changes in version 1.16.3                        

- fixed clustering method not applied

                       Changes in version 1.16.2                        

- check if taxonomy DB exists in cwd instead of r/library

                       Changes in version 1.16.1                        

- fixed bug wrong order in detailed plot

- fixed bug mVar for co-orthologs

[Pigengene](/packages/Pigengene)
---------

                Changes in version 1.29.10 (2023-11-15)                 

Changes in existing functions

- More diagnostic messages from scoreCandidates.

                 Changes in version 1.29.8 (2023-11-13)                 

Bug Fixes

- Using "/" is avoided in determine.modules and learn.bn
  functions to support running them on Windows.  The "toLocal"
  task in learn.bn() and sbatch() still do not support Windows.

General

- inst/script/bn.calculation.job.R is added as an example for the
  bnCalculationJob argument.

                 Changes in version 1.29.6 (2023-11-10)                 

Changes in existing functions

- Setting toCompact to TRUE in make.decision.tree() has the same
  effect as NULL.

                 Changes in version 1.29.4 (2023-11-09)                 

Bug Fixes

- The averaged.network() function now correctly passed test data
  to the module.heatmap() function.

Changes in existing functions

- More details on the effect of test data in the documentation of
  one.step.pigengene() and make.decision.tree() functions.
  make.decision.tree() returns more information on performance.

[PIPETS](/packages/PIPETS)
------

                 Changes in version 0.99.0 (2023-09-26)                 

- Submitted to Bioconductor

[PIUMA](/packages/PIUMA)
-----

                        Changes in version 1.0.0                        

- PIUMA has been released!

[plotgardener](/packages/plotgardener)
------------

                        Changes in version 1.9.5                        

NEW FEATURES

- plotHicSquare has added yaxisDir parameter for flipping the
direction of genomic coordinates on the y-axis.

                        Changes in version 1.9.4                        

BUG FIXES

- plotManhattan accommodates GRanges data input.

                        Changes in version 1.9.3                        

BUG FIXES

- plotGenes appropriately handles tibble input for geneHighlights.

                        Changes in version 1.9.1                        

NEW FEATURES

- plotIdeogram has added flip parameter to allow for flipping the
ideogram so the end can be before/above the start.

BUG FIXES

- plotManhattan pch mapping is compatible when number of data points
is less than the number of levels of a colorby column.

                        Changes in version 1.9.0                        

Version bump for Bioconductor 3.18 release.

[PLSDAbatch](/packages/PLSDAbatch)
----------

                       Changes in version 0.99.3                        

- Date: 2023-12-21
- Text: Debug percentileofscore
- Details: Debug the function percentileofscore(), thus the generated
data frame has rownames and colnames as the input data frame to
ensure the operation of function percentile_norm()

                       Changes in version 0.99.2                        

- Date: 2023-11-14
- Text: Updates for Bioconductor
- Details: Update several documents
- DESCRIPTION: Add a URL field pointing to the GitHub repository
- vignettes: Wrap as many texts as possible to a character limit
of 80; Use library() calls instead of sapply(); Use
TreeSummarizedExperiment package to reshape datasets; Remove
colorize()
- R: Update for() loops; Remove extra ### lines

                       Changes in version 0.99.1                        

- Date: 2023-05-26
- Text: Update Imports
- Details: Move pheatmap, vegan, Biobase, BiocStyle from Suggests to
Imports to ensure the generation of vignette

                       Changes in version 0.99.0                        

- Date: 2023-04-03
- Details: This version is equal to the version 0.2.3 on GitHub
(https://github.com/EvaYiwenWang/PLSDAbatch-source-code)

[pmp](/packages/pmp)
---

                       Changes in version 1.15.1                        

- fix plots due to breaking changes in ggplot2

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.11.1 (2023-05-15)                 

- Vignette spellcheck

[podkat](/packages/podkat)
------

                       Changes in version 1.35.1                        

- changed e-mail address of maintainer

- updated README.md and formatting of package vignette

- updated CITATION

- updated references in documentation

                       Changes in version 1.35.0                        

- new branch for Bioconductor 3.19 devel

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

                       Changes in version 1.13.26                       

- New POMA theme and colorblind-friendly palette
- Available sample normalization (sum and quantile)
- New feature normalization methods
- Extensive review and improvement of all POMA functions
- Major documentation updates
- Rename PomaSummarizedExperiment to PomaCreateObject
- Auto-recognition of variable types and automatic variable
re-labeling in PomaCreateObject
- Available violin plots with PomaBoxplots
- New functions PomaPCA and PomaPLS as stand-alone functions from the
old PomaMultivariate (deprecated)
- Other new functions: PomaLM, PomaLMM
- Post-hoc pairwise comparisons in PomaUnivariate
- Update tests and vignettes

[procoil](/packages/procoil)
-------

                       Changes in version 2.31.1                        

- changed e-mail address of maintainer

- updated README.md and formatting of package vignette

- updated references in documentation

[proDA](/packages/proDA)
-----

                        Changes in version 1.17                         

- Check that sigma2 is larger than zero (see #21. Thanks to
  @elena-krismer)

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.43                         

Changes in version 1.43.2

- Fix/update dunkley2006param object.

Changes in version 1.43.1

- Fix syntax in man pages.

Changes in version 1.43.0

- New devel version

[ProtGenerics](/packages/ProtGenerics)
------------

                       Changes in version 1.35.4                        

- Move `backendParallelFactor` and `backendBpparam` and
  `supportsSetBackend`
  generics from `Spectra`.

                       Changes in version 1.35.3                        

- Move `backendInitialize`, `backendMerge`, `setBackend`, `isReadOnly`,
  `peaksData`, `peaksData<-` and `peaksVariables`, `filterRanges`,
  `filterValues`, `filterPrecursorMzRange`, `filterPrecursorMzValues`,
  `filterProductMzRange`, `filterProductMzValues` generics from the
  Spectra package.

                       Changes in version 1.35.2                        

- Add `filterSpectra` generic.

                       Changes in version 1.35.1                        

- Add `filterFeatures` generic.

[psichomics](/packages/psichomics)
----------

                       Changes in version 1.28.1                        

- Update GTEx data download based on new URL endpoints
- Documentation:
- Update function documentation based on bug fixes below
- Update license copyright years
- Bug fixes:
- Fix deprecated used of size in ggplot2::geom_line
- Fix inconsistent argument in survfit.survTerms()
- Fix unit test issues when testing low coverage using random PSI
values
- Update pkgdown and R CMD check automation in GitHub Actions

[PSMatch](/packages/PSMatch)
-------

                         Changes in version 1.7                         

PSMatch 1.7.2

- Fix connected component dim names in show().

PSMatch 1.7.1

- In addFragments() use ... to pass parameters to
calculateFragments().

PSMatch 1.7.0

- New Bioc devel.

[PureCN](/packages/PureCN)
------

                       Changes in version 2.10.0                        

NEW FEATURES

- adjustLogRatio function for adjusting a tumor vs normal coverage
  ratio for purity and ploidy. Useful for downstream tools that
  expect ratios instead of absolute copy numbers such as GISTIC.
  Thanks @tedtoal (#40).

SIGNIFICANT USER-VISIBLE CHANGES

- Provide interval-level likelihood scores in runAbsoluteCN return
  object. Thanks @tinyheero (#335).

- Documentation updates. Thanks @ddrichel (#325).

BUGFIXES

- Bugfix #296 was not merged into the developer branch and did not make
  it into 2.8.0.

- Log ratios not shiften to median of sample medians as intended
  (#356).
  Thanks @sleyn.

- Fixed crash with small toy examples (fewer than 2000 baits, #363)

[pwalign](/packages/pwalign)
-------

                        Changes in version 1.0.0                        

- First version of the package that is ready for general use.

[qcmetrics](/packages/qcmetrics)
---------

                        Changes in version 1.41                         

qcmetrics 1.41.1

- Fix bug when when setting the env.
- Remove template reporting argument.

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.13                         

QFeatures 1.13.7

- Fix filterFeatures() when the filter variable also exists in the
global environment (issue #208).

QFeatures 1.13.6

- Migrate .splitS&#91;C&#93;E() unit test form scp to QFeatures.

QFeatures 1.13.5

- Reuse readQFeatres() params in readQFeatureFromDIANN().

QFeatures 1.13.4

- Use DIA-NN example data from MsDataHub in readQFeaturesFromDIANN().

QFeatures 1.13.3

- Fix is.vector() (see issue #203)
- readQFeatures() multi-set support (ported from scp::readSCP() - see
issue #199).
- new readQFeaturesFromDIANN() function to import DIA-NN report files
(ported from scp::readSCPfromDIANN() - see issue #199).

QFeatures 1.13.2

- Move the filterFeatures generic method to ProtGenerics.

QFeatures 1.13.1

- Fix typo in vigette.

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

                       Changes in version 1.21.1                        

- Added function `coefVar.R` to calculate and plot coefficient of
  variance
  between all samples in a sample group.

[raer](/packages/raer)
----

                        Changes in version 1.1.3                        

- Fixed bug in smart-seq2 pileup_cells reporting

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.28.0                        

New features

- Added citation information to the vignette and package. See
citation(package = "RaggedExperiment")

[Rarr](/packages/Rarr)
----

                         Changes in version 1.3                         

- Added support for using the zstd compression library for reading
and
writing.

[rawDiag](/packages/rawDiag)
-------

                Changes in version 0.99.25 (2024-02-28)                 

- Test
- add test case for rawDiagServerModule
- Documentation
- add Visualization section in vignette
- add shiny and FAQ section in vignette
- use EH4547(DIA) and EH3222(DDA) data from tatare pkg
- add FAQ section in vignette
- Refactor
- if trailer information is missing add NA
- Feature
- handle "Orbitrap Resolution:" in trailer

                Changes in version 0.99.23 (2024-02-15)                 

- pass R CMD check on bioconductor with Status: OK

                Changes in version 0.99.19 (2024-02-14)                 

- NAMESPACE
- renamed read.raw to readRaw
- Documentation
- add an Introduction and Installation section in vignette
- removed commented code that isn't run
- Code
- removed paste cmd in message
- use the BiocParallel
- refactor rawDiag shiny application, e.g., buildRawDiagShinyApp
returns a shiny::shinyApp object

                 Changes in version 0.99.9 (2024-01-05)                 

- Add NEWS.md file

                 Changes in version 0.99.1 (2023-12-06)                 

- refacor existing code https://github.com/fgcz/rawDiag
- use mono assemplies provides through Bioconductor rawrr package
- camel case for public plot functions, e.g., rename PlotCycleLoad to
plotCycleLoad
- submit to bioconductor #3251

[rawrr](/packages/rawrr)
-----

                 Changes in version 1.10.1 (2023-11-02)                 

- Fix index error #67.

[rBLAST](/packages/rBLAST)
------

                 Changes in version 1.1.1 (2024-04-30)                  

Changes

- Added description of how to deal with multi-part databases.

- Fixed package anchors in man pages.

                        Changes in version 1.1.0                        

Changes

- makeblastDB() gained parameters db_name, hash_index and verbose.

- added has_blast() and made tests, man page code and vignette code
  conditional.

- added blast_db_cache() and blast_db_get() to manage downloading BLAST
  DBs
  using BiocFileCache.

[RCy3](/packages/RCy3)
----

                       Changes in version 2.24.0                        

- Use stringi to replace uchardet

- Add test functions

- Doc fixes:
  - add "c" to anchor choices

- Bug fixes:
  - fix setNodePositionBypass to support network and base.url
  - use viridis color palette for continuous mapping, #210

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.17.3 (2024-02-01)                 

- Fixed a bug in vignette showcasing the reanalysis of public data.

                 Changes in version 1.17.2 (2023-11-28)                 

- Added a set of new functions to support ReactomeGSA's features to
  load public datasets.

- Added a new vignette to showcase how public datasets can easily be
  loaded using the ReactomeGSA package.

                       Changes in version 1.17.1                        

- Fixed vignettes

[ReactomePA](/packages/ReactomePA)
----------

                       Changes in version 1.47.1                        

- Add a disclaimer in package Description to claim that this package
is not affiliated with the Reactome team (2023-11-17, Fri)

[recount](/packages/recount)
-------

                       Changes in version 1.29.1                        

BUG FIXES

- Merged a pull request by @hpages which addresses some changes to
GenomicFeatures. See https://github.com/leekgroup/recount/pull/24
for details.

[recoup](/packages/recoup)
------

                 Changes in version 1.31.1 (2024-03-14)                 

NEW FEATURES

- No new features.

BUG FIXES

- Switched to txdbmaker.

[RedeR](/packages/RedeR)
-----

                        Changes in version 3.0.0                        

- Major code refactoring.

- Improved readability.

- Reduced code complexities.

- Improved compatibility with standard base graphics.

[RegEnrich](/packages/RegEnrich)
---------

                       Changes in version 1.13.4                        

- Fix a bug in GRN.R
- Fix a bug in DEA_LRT_DESeq2
- Move BiocStyle to Imports

                       Changes in version 1.13.1                        

- Fix the bug of re-building the vignette on Windows.

[rgoslin](/packages/rgoslin)
-------

                        Changes in version 1.8.0                        

Please note that this Bioconductor version is based on Goslin
version 2.2.0. See the Goslin repository and Goslin C++ repository
for
more details.

Improvements

- Added more trivial mediators
- Added prostglandins
- Added oxylipins
- Updated functional groups

Bug Fixes

- Fixed oxo handling
- Updated lipid class output for plasmanyl / plasmenyl to allow
distinction

[rGREAT](/packages/rGREAT)
------

                        Changes in version 2.5.2                        

- add `getGeneSetsFromOrgDb()`

                        Changes in version 2.5.1                        

- support extending from the whole genes

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.48.0                        

CHANGES

- R complex types can now be written to HDF5.  These will be
  stored as a compound datatype with two elements (r & i)
  representing the real and imaginary parts.

- Functions H5Screate_simple and H5Sset_extent_simple now accept
  numeric values to the dim and maxdim arguments, allowing the
  creation of HDF5 dataspaces larger than R's maximum integer
  value.  (Thanks to @hpages for reporting this and providing a
  patch https://github.com/grimbough/rhdf5/pull/140).

- Messages about the presence of INT_MIN in datasets created
  outside of rhdf5, and how R will convert them to NA have been
  moved from H5Dread to the high-level h5read function.  (See
  https://github.com/LTLA/scRNAseq/issues/44 for more details).

Bug fixes

- Corrected an issue where the function prototype for
  _h5fileLock() differed from the actual implementation.

- Addressed a bug where fixed length string attributes would be
  one character shorter than they should be. Backported to rhdf5
  2.46.1.  (Thanks to Aaron Lun @LTLA for reporting this
  https://github.com/grimbough/rhdf5/issues/132).

- Fixed an issue where zero length datasets of uint32, int64 or
  uint64 datatypes could not be read.  This would fail with an
  error message saying there was not enough memory. Backported to
  rhdf5 2.46.1.  (Thanks to Aaron Lun @LTLA for reporting this
  https://github.com/grimbough/rhdf5/issues/134).

[rhdf5filters](/packages/rhdf5filters)
------------

                       Changes in version 1.16.0                        

CHANGES

- rhdf5filters no longer sets the HDF5_PLUGIN_PATH environment
variable when it is loaded. Instead this is handled by the
H5PLprepend() function in rhdf5.

BUG FIXES

- Fixed issue compiling VBZ filter when R was installed via conda.
Backported to version 1.14.1. (Reported in
https://github.com/grimbough/rhdf5filters/issues/20)

[Rhdf5lib](/packages/Rhdf5lib)
--------

                        Changes in version 1.26                         

Bug fixes

- Updated parts of the bundled HDF5 configure scripts (config.guess &
  config.sub) which were very old and didn't recognise Apple arm64
  systems
  when running R installed via conda. (Thanks to @eyalbenda for
  reporting,
  https://github.com/grimbough/Rhdf5lib/issues/54) Backported to
  Rhdf5lib version 1.24.2

- The configure script has been updated to build from source on
  Windows/aarch64 architecture, rather than trying to use the
  unsuitable
  pre-compiled Windows x64 binaries.  Thanks to Tomas Kalibera for the
  patch.

[Rhisat2](/packages/Rhisat2)
-------

                       Changes in version 1.19.2                        

- Swap GenomicFeatures dependency for txdbmaker

                       Changes in version 1.19.1                        

- Exclude unsupported flags in Makefile for Linux aarch64 platform

- Make mask2iupac a signed char array

[Rhtslib](/packages/Rhtslib)
-------

                        Changes in version 3.0.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Update htslib to 1.18 (was 1.15.1)

[RNAmodR](/packages/RNAmodR)
-------

                 Changes in version 1.17.1 (2024-03-18)                 

- Bugfix due to upstream changes in internal functions

[RNAseqCovarImpute](/packages/RNAseqCovarImpute)
-----------------

                 Changes in version 1.1.4 (2024-04-29)                  

- bug fix above where dfs were Inf fixed for df AND df_bayes

                 Changes in version 1.1.3 (2024-04-12)                  

- fixed bug in get_gene_bin_intervals where all genes were not included
  in some cases

- fixed bug where dfs were Inf in cases where information loss due to
  missingness (lambda) was 0

- Added MI PCA method where we 1) perform PCA on the log CPM values for
  all genes, 2) create M
  imputed datasets where the imputation predictor matrix includes all
  covariates and the optimum
  number of PCs to retain (e.g., based on Horn’s parallel analysis or
  the number of PCs that account
  for >80% explained variation), 3) conduct the standard limma-voom
  pipeline on each M imputed dataset,
  and 4) pool the results with Rubins’ rules to produce combined
  coefficients, standard errors, and P-values

                 Changes in version 1.0.2 (2023-10-26)                  

- fixed bug in combine_rubins where sigma values were not properly made
  into a numeric vector

                 Changes in version 1.0.1 (2023-10-26)                  

- fixed bug in limmavoom_imputed_data_list_helper function where sigma
  values were not properly output

[rols](/packages/rols)
----

                        Changes in version 2.99                         

CHANGES IN VERSION 2.99.4

- Typo and references in manual pages.

CHANGES IN VERSION 2.99.3

- Add orcid.

CHANGES IN VERSION 2.99.2

- Fix different URI locations (see #42)

CHANGES IN VERSION 2.99.1

- Fix logical syntax in url.

CHANGES IN VERSION 2.99.0

- Refactoring to use REST API for OLS4.
- REST queries now use httr2 instead of superseded httr.
- The term(s) and property constructors are capitalised as Term(),
Terms() and Properties().

[ropls](/packages/ropls)
-----

                       Changes in version 1.35.4                        

- gg_scoreplot: minor update

                       Changes in version 1.35.2                        

- gg_scoreplot: update to handle OPLS models

[rRDP](/packages/rRDP)
----

                 Changes in version 1.37.3 (2024-03-26)                 

Changes

- Updated to the lated RDP Classifier version 2.14 released in August
  2023 which contains the bacterial and archaeal taxonomy training set No.19.

[Rsubread](/packages/Rsubread)
--------

                       Changes in version 2.18.0                        

- Removed the 'maxMOp' parameter from featureCounts.

[RTCGAToolbox](/packages/RTCGAToolbox)
------------

                       Changes in version 2.34.0                        

Significant user-visible changes

- getReport has been removed from the package.

Bug fixes and minor improvements

- Remove unnecessary imports e.g., from RCircos.

[Rvisdiff](/packages/Rvisdiff)
--------

                 Changes in version 1.1.1 (2024-02-12)                  

- Added SVG exportation for line chart, boxplot and heatmap

- PNG styles fixed

- PNG resolution enhanced

- Addded mouse panning to explore boxplot

[S4Arrays](/packages/S4Arrays)
--------

                        Changes in version 1.4.0                        

- No changes in this version.

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.42.0                        

NEW FEATURES

- Add 'make.names' argument to as.data.frame() method for DataFrame.
  See https://github.com/Bioconductor/S4Vectors/issues/121

SIGNIFICANT USER-VISIBLE CHANGES

- Document global options 'showHeadLines' and 'showTailLines'.

[SAIGEgds](/packages/SAIGEgds)
--------

                        Changes in version 2.4.0                        

- update `seqFitLDpruning()`

- reduce the memory usage in `seqAssocGLMM_SPA()` when the genotype
  file
  is not split by chromosomes

                        Changes in version 2.2.1                        

- fix the compiler warning: -Wformat-security

- fix a compiler issue with gcc v7.5

[sangeranalyseR](/packages/sangeranalyseR)
--------------

                        Changes in version 99.1                         

- Base class: SangerReads is designed to store each forward/revers
  reads.

[saseR](/packages/saseR)
-----

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[scanMiR](/packages/scanMiR)
-------

                 Changes in version 1.9.3 (2024-02-09)                  

- since ggseqlogo was removed from CRAN, it was replaced and the logo
  plot layouts are slightly different

[scAnnotatR](/packages/scAnnotatR)
----------

                 Changes in version 1.9.2 (2024-04-01)                  

- Fixed bug in classify_cells_sce.

                 Changes in version 1.9.1 (2024-03-21)                  

- Adapted vignettes to address changes in scRNAseq package.

[SCArray](/packages/SCArray)
-------

                       Changes in version 1.12.0                        

- update according to DelayedArray (>= v0.29.0)

[SCArray.sat](/packages/SCArray.sat)
-----------

                        Changes in version 1.4.0                        

- update for Seurat v5

[scMitoMut](/packages/scMitoMut)
---------

                       Changes in version 0.99.5                        

- Update README add biocondutor installing info
- Make sure the subsetting cells and loc must be unique in
scMitoMutObj
- Add alternative reads count threshold for loc filtering
(alt/_count)

                       Changes in version 0.99.2                        

Initiate package.

[scMultiSim](/packages/scMultiSim)
----------

                       Changes in version 0.99.9                        

Prepare for the Bioconductor release

- Fix build errors

                       Changes in version 0.99.8                        

Prepare for the Bioconductor release

- Tidy up the code, add more comments

                       Changes in version 0.99.7                        

Prepare for the Bioconductor release

- Fix build errors

[scp](/packages/scp)
---

                        Changes in version 1.13                         

scp 1.13.6

(nothing yet)

scp 1.13.5

- fix: fixed small error in degrees of freedom
- fix: break workflow upon infinite values

scp 1.13.4

- fix: first drop variable before centering numerical variables (see
#54)

scp 1.13.3

- Refactor readSCP() and readSCPFromDIANN() to use new QFeatures
implementations.

scp 1.13.2

- New scplainer workflow and citation
- Add addReducedDims() (see #52)
- fix: logFC and associated SE are now correctly computed

scp 1.13.1

- docs: fixed vignette about reporting missing values
- docs: add a QFeatures figure to the nutshell vignette

scp 1.13.0

- New Bioconductor 3.18 (devel) release

[scRepertoire](/packages/scRepertoire)
------------

                 Changes in version 2.0.0 (2024-01-10)                  

NEW FEATURES

- Added percentAA()
- Added percentGenes()
- Added percentVJ()
- Added percentKmer()
- Added exportClones()
- Added positionalEntropy()
- Added positionalProperty()
- Changed compareClonotypes to clonalCompare()
- Changed clonotypeSizeDistribution to clonalSizeDistribution()
- Changed scatterClonotypes to clonalScatter()
- Changed quantContig to clonalQuant()
- Changed highlightClonotypes to highlightClones()
- Changed lengthContigs to clonalLength()
- Changed occupiedscRepertoire to clonalOccupy()
- Changed abundanceContig to clonalAbundance()
- Changed alluvialClonotypes to alluvialClones()
- Added features to clonalCompare() to allow for highlighting
sequences, relabeling clonotypes.

UNDERLYING CHANGES

- Removed internal .quiet() function.
- .theCall() now allows for a custom header/variable and checks the
colnames.
- Replaced data arguments to be more descriptive: df is now
input.data, dir is now input, and sc is now sc.data
- Deep clean on the documentation for each function for increased
consistency and explainability
- StartracDiversity() metric re-implemented to remove startrac-class
object intermediary
- Implemented powerTCR locally to reduce dependencies and continue
support
- Universalized underlying function language and intermediate
variables
- License change to MIT
- group.by and split.by have been consolidated into single group.by
parameter
- Added support for Immcantation pipeline, .json, Omniscope, and
MiXCR
formats for loadContigs()
- Made GitHub.io website for support/vignettes/FAQ
- Restructured NEWS Tracking
- Added testthat for all exported and internal functions
- Fixed issue with clonalQuant() for instance of scale = FALSE and
group.by being set.
- clonalDiversity() no longer automatically orders samples.
- Remove order parameter from clonalQuant(), clonalLength(), and
clonalAbundance()
- x.axis parameter in clonalDiversity() separated from group.by
parameter
- filtering chains will not eliminate none matching chains.

DEPRECATED AND DEFUNCT

- Deprecate stripBarcodes()
- Deprecate expression2List() (now only an internal function).
- Deprecate checkContigs()

[scRNAseqApp](/packages/scRNAseqApp)
-----------

                       Changes in version 1.3.25                        

- Fix the abstract issue for edit module.

                       Changes in version 1.3.24                        

- Update the description.

                       Changes in version 1.3.23                        

- Add a hyperlink to package URL at footer.

                       Changes in version 1.3.22                        

- Fix the duplicated rownames.

                       Changes in version 1.3.21                        

- Add search button in the home page.

                       Changes in version 1.3.20                        

- Change the default waffle plot group to celltype.

                       Changes in version 1.3.19                        

- Update abstract in the admin edit tab.

                       Changes in version 1.3.18                        

- Add abstract parameter for createAppConfig

                       Changes in version 1.3.17                        

- Fix the issue if no ip location information available when update
  from the older version.

- add the possibility of the feature names to replace the rownames.

                       Changes in version 1.3.16                        

- Fix the missmatch UI and ID.

                       Changes in version 1.3.15                        

- Fix the no id column issue when update the database.

                       Changes in version 1.3.14                        

- Update the maxOptions from 6 to 50.

                       Changes in version 1.3.13                        

- Add groupby parameter for correlation for proportion plot.

                       Changes in version 1.3.12                        

- Add correlation for proportion plot.

                       Changes in version 1.3.11                        

- Add Wilcoxon test for violin plot.

                       Changes in version 1.3.10                        

- Add Chisq test for proportion.

                        Changes in version 1.3.9                        

- Fix some typos.

                        Changes in version 1.3.8                        

- Change dependence from hdf5r to rhdf5.

                        Changes in version 1.3.7                        

- Fix the wrong downloader handler in downloader module.

                        Changes in version 1.3.6                        

- Add downloader module.

                        Changes in version 1.3.5                        

- Add gene name map to sql database.

                        Changes in version 1.3.4                        

- Fix the issue of missing stats for cell numbers.

                        Changes in version 1.3.3                        

- Add user visit info and gene symbols into the sql database.

                        Changes in version 1.3.2                        

- Add dataset configuration into the sql database.

                        Changes in version 1.3.1                        

- Make it compatable with SeuratObject V5.

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.44.0                        

UTILITIES

- tweak the display of progress information in `seqVCF2GDS()`

- `seqVCF_Header(, getnum=TRUE, verbose=TRUE)` to show the progress
  information for scanning the VCF file

- new `seqGetData(, "$dosage_alt2")` and `seqGetData(, "$dosage_sp2")`
  for
  sex chromosomes, when the alleles are partially missing (e.g.,
  genotypes
  on chromosome X for males)

- new 'verbose.clean' in `seqExport()` to control how much information
  to
  be displayed

                       Changes in version 1.42.4                        

BUG FIXES

- `seqGetData(, "$dosage_alt")` and `seqGetData(, "$dosage_sp")` work
  correctly when the ploidy is >2 and there are missing alleles

- fix a bug that `seqParallel()` does not call a user-defined
  '.combine'
  when 'parallel=1'

                       Changes in version 1.42.1                        

UTILITIES

- update the help files of `seqBlockApply()` and `seqUnitApply()`

- detect the output filename extension in `seqGDS2VCF()` without
  considering the case of the characters, supporting .gz, .bgz, .bz and
  .xz
  as a filename extension

- fix the compiler warning: -Wformat-security

- new option 'include.pheno=TRUE' in `seqBED2GDS()`

[shiny.gosling](/packages/shiny.gosling)
-------------

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[signeR](/packages/signeR)
------

                        Changes in version 2.5.1                        

- Preserve dimnames from the given input matrices

[signifinder](/packages/signifinder)
-----------

                        Changes in version 1.6.0                        

- Add 15 new signatures.

- Update the evaluationSignPlot function to return an overall
  goodness score.

- Add getSignGenes function to access the signatures gene list.

- When multiple signatures are plotted, values are reported in
  z-scores and not rescaled between 0 and 1 as before.

[SimBu](/packages/SimBu)
-----

                 Changes in version 1.5.3 (2023-11-21)                  

- Added 'seed' parameter to simulation

                 Changes in version 1.5.2 (2023-11-08)                  

- Updated SimBu to be compatibel with Seurat v5

[simPIC](/packages/simPIC)
------

                 Changes in version 0.99.7 (2024-04-14)                 

- Updating title.
- Simplifying estimate sparsity description.

                 Changes in version 0.99.6 (2024-03-16)                 

- Replacing sapply with vapply in simPICsimulate line 238.

                 Changes in version 0.99.5 (2024-03-16)                 

- Addressing notes after reviewer comments.
- Updated simPIC count documentation.
- Updated 1:nCells to seq_len in simPICsimulate line 239.
- Added documentation for testdata.R.

                 Changes in version 0.99.4 (2024-03-10)                 

- Fixed issue with system files found that should not be Git tracked
in BiocCheck::BiocCheckGitClone

                 Changes in version 0.99.3 (2024-03-10)                 

- Updating R dependency to R (>=4.4.0)

                 Changes in version 0.99.2 (2024-03-09)                 

- Fixing issue with system files found in BiocCheck

                 Changes in version 0.99.1 (2024-03-09)                 

- Addressed reviewer suggestions
- Major changes include improved docs, including a package man page,
text file and code to reproduce test data

                 Changes in version 0.99.0 (2024-02-28)                 

- Submitted to Bioconductor

[SingleCellAlleleExperiment](/packages/SingleCellAlleleExperiment)
--------------------------

                       Changes in version 0.99.5                        

Changes regarding the second round of comments during the
Bioconductor
submission review:

- Minor changes to SCAE-show method

- Extended and moved installation-checks for packages used in
optional
functionalities to the SCAE-constructor

- Extended checks for the rowData setter

- Extended unit tests to increase coverage

                       Changes in version 0.99.4                        

Changes regarding the first round of comments during the Bioconductor
submission review:

Removed multiple dependencies, now listed as suggested

- computation of the logcounts assay now optional

- computation of a knee plot and according inflection/knee points now
optional

- computation of new gene-symbols now optional

- minor changes in tests (needs to be extended to test for new
optional code); also extension of the tests to check new code

- major changes in package structure and organisation

- reworked vignette to include optional features

                       Changes in version 0.99.3                        

- retriggering bioc build during package submission.

                       Changes in version 0.99.2                        

- minor changes in the vignette to reduce build-time. Content is the
same.

                       Changes in version 0.99.1                        

- minor changes made regarding the Bioconductor submission review of
scaeData which also affected this package.

                       Changes in version 0.99.0                        

- The package is defining a S4 class that is extending the
SingleCellExperiment class. The multi-layer data structure
integrates data for immune genes at allele and functional level.

- Workflow: To be able to generate a SingleCellAlleleExperiment
object
the data has to be generated by the connected workflow single-cell
ImmunoGenomic Diversity (scIGD) which performs allele typing and
quantification of genes and typed alleles. Use the
read_allele_counts() function to read in the data and generate an
SCAE object. This function also offers parameters to perform
filtering on the data, as well as visualize the filtering step in a
so called knee plot.

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.12.2 (2024-01-28)                 

- Added support for Seurat V5

                 Changes in version 2.12.1 (2024-01-10)                 

- Updates to documentation
- Fixes to runTSCAN and plotSeurat Genes
- Added support for flat file import into SCTK-QC
- Fixed directory issue in importCellRanger
- Added Bubble plot to Shiny GUI
- Updated Dockerfile

[sitadela](/packages/sitadela)
--------

                 Changes in version 1.11.1 (2024-03-07)                 

NEW FEATURES

- Updated Ensembl supported versions based on biomaRt

[sketchR](/packages/sketchR)
-------

                       Changes in version 0.99.2                        

- Version bump to trigger new SPB build, no changes

                       Changes in version 0.99.1                        

- Version bump to trigger new SPB build, no changes

                       Changes in version 0.99.0                        

- Prepare for Bioconductor submission


[smartid](/packages/smartid)
-------

                       Changes in version 0.99.5                        

- Update scale_mgm() function adding pooled SD option, add details
for
scale function.

                       Changes in version 0.99.4                        

- Add details for TF, IDF, IAE functions.

                       Changes in version 0.99.3                        

- Bump R version dependency to >= 4.4 and add details for TF, IDF,
IAE
functions.

                       Changes in version 0.99.2                        

- Added test for gs_score() function.

                       Changes in version 0.99.1                        

- Ready for submission to Bioconductor.

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[smoothclust](/packages/smoothclust)
-----------

                 Changes in version 0.99.0 (2024-03-17)                 

- version for submission to Bioconductor

[SNPRelate](/packages/SNPRelate)
---------

                       Changes in version 1.38.0                        

UTILITIES

- faster `snpgdsLDpruning()` and new multi-threading implementation

- new option 'autosave' in `snpgdsLDpruning()`

- new option `start.pos="random.f500"` in `snpgdsLDpruning()`, allowing
  faster genotype scanning (it is used by default now). To be
  compatible
  with the previous versions, please specify `start.pos="random"`
  manually

- Fix the welcome message when running in a non x86 system

                       Changes in version 1.36.1                        

UTILITIES

- fix the compiler warning: -Wformat-security

- set 'maf=0.005' & 'missing.rate=0.05' by default in
  `snpgdsLDpruning()`
  (for WGS data)

- show a progress bar in `snpgdsLDpruning()` when 'verbose=TRUE'

- update the help files of `snpgdsAdmixProp()` and `snpgdsAdmixPlot()`

[SparseArray](/packages/SparseArray)
-----------

                        Changes in version 1.4.0                        

NEW FEATURES

- The following operation are now multithreaded (via OpenMP):
  - SparseMatrix multiplication and crossprod()
  - the matrixStats methods for SparseMatrix objects
  See ?set_SparseArray_nthread for more information.

- Define some of the basic coercion methods that used to be defined in
  the Matrix package but that the lazy Matrix maintainers have decided
  to deprecate in Matrix 1.7-0 (e.g. coercion from matrix to
  dgRMatrix).

BUG FIXES

- Fix aperm(&lt;SVT_SparseArray&gt;) when inner dimensions are not permuted
  See https://github.com/Bioconductor/SparseArray/issues/6

- Fix long-standing bug in C function
  _dotprod_leaf_vector_and_double_col().
  In some circumstances the function was reading beyond a C array.
  Observed effects varied from no observed effects (on Linux) to
  breaking
  crossprod() unit tests on arm64 Mac (consistently), on Windows
  (sporadically), and also probably on PowerPC (or possibly a crash
  sometimes on this platform.
  See https://github.com/Bioconductor/SparseArray/issues/2

[sparseMatrixStats](/packages/sparseMatrixStats)
-----------------

                        Changes in version 1.15                         

- Throw error if length of 'center' in colVars, colSds, colMads (and
  corrresponding row-functions) does not match size of 'x'.
  This behavior was never supported by `sparseMatrixStats` and
  previously
  could return incorrect results. This change matches the upcoming
  behavior of
  matrixStats
  (https://github.com/HenrikBengtsson/matrixStats/issues/254).

[SpatialDecon](/packages/SpatialDecon)
------------

                 Changes in version 1.13.2 (2023-12-28)                 

- Compatibility with Seurat v5

[SpatialFeatureExperiment](/packages/SpatialFeatureExperiment)
------------------------

                        Changes in version 1.6.0                        

- Changed defaults from sample_id = NULL to sample_id = 1L when
dealing with 1 sample or "all" when dealing with multiple samples
- dim method for BioFormatsImage that doesn't load the image into
memory
- Deal with univariate spatial results in featureData in cbind and
changeSampleID
- Fixed super embarrassing bug in cbind that fails when combining
more
than 2 SFE objects
- Updated readXenium for XOA v2
- Updated BioFormatsImage to store affine transform info rather than
converting to EBImage after transform
- Speed up affine transformation of sf geometries with sfheaders
- Coercion from Seurat to SFE
- SpatRasterImage and EBImage directly inherit from SpatRaster and
Image respectively so the user no longer needs to call imgRaster
every time they plots or operates on the image, which I find really
annoying.
- Changed name EBImage to ExtImage to reduce confusion
- Bug fixes on image affine transformation
- Exporting some util functions: aggBboxes, getPixelSize, and
imageIDs

                        Changes in version 1.5.2                        

- Added readXenium (for XOA v1)
- Added BioFormatsImage and EBImage classes to deal with Xenium
OME-TIFF
- Conversion between SpatRasterImage, BioFormatsImage, and EBImage
- Overhaul of geometry operation functions for images and SFE objects
for the new image classes, including bbox, crop, and affine
transforms
- Don't throw error when there are no rows or columns left after [
subsetting
- cbind for multiple samples that have rowGeometry
- Rewrote df2df with the much faster sfheaders, deprecating the less
efficient BPPARAM argument

                        Changes in version 1.5.1                        

- Added support for rowGeometry and transcript spots
- Reformat transcript spot files from Vizgen and CosMX
- Improved readVizgen for transcript spots
- Added readCosMX

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 2.9.6 (2024-04-16)                  

- The option of pairwise comparision is added in spatial enrichment.

- The layout of static network graph is improved.

                 Changes in version 2.9.4 (2024-02-10)                  

- Added the formal citation.

[SpatialOmicsOverlay](/packages/SpatialOmicsOverlay)
-------------------

                        Changes in version 1.4.0                        

- Compatibility with updated labworksheet

- IPA compatibility

                        Changes in version 1.2.1                        

- Corrected labworksheet attached to package

- Added tibble plottingFactor

- scale bar printing in mm

- Bugs fixes:
  - SpatialPosition printing
  - saving xml files
  - ROI matching on string values
  - overcropping on recolored images

[Spectra](/packages/Spectra)
-------

                        Changes in version 1.13                         

Changes in 1.13.8

- Add estimatePrecursorMz() function to estimate the precursor m/z
for
DDA fragment spectra from previous MS1 spectra issue #315.

Changes in 1.13.7

- Move generics backendBpparam(), backendParallelFactor() and
supportsSetBackend() to ProtGenerics. Required ProtGenerics version
1.35.4 or higher.

Changes in 1.13.6

- Add filterRanges() and filterValues() functions to allow filtering
of a Spectra object based on ranges or similarities of any existing
spectraData variables.

Changes in 1.13.5

- Move generics to ProtGenerics. Requires ProtGenerics version
1.35.3.

Changes in 1.13.4

- Add entropy and nentropy functions to allow to calculate the
(normalized) entropy for each spectrum.

Changes in 1.13.3

- Fix issue in setBackend that might cause chunk-wise processing to
be
not run.

Changes in 1.13.2

- Add possibility to enable and perform chunk-wise (parallel)
processing to Spectra: add functions processingChunkSize,
backendParallelFactor and processingChunkFactor to set or get
definition of chunks for parallel processing. All functions working
on peaks data use this mechanism which is implemented in the
internal .peaksapply function. The Spectra object gains a new slot
"processingChunkSize" that is used to define the size of the
processing chunks for the Spectra. See also issue #304. This ensures
processing also of very large data sets.

Changes in 1.13.1

- Fix issue with bin function (see issue #302). Addition of zero.rm
parameter to prevent creation of empty bins.

[SPIAT](/packages/SPIAT)
-----

                        Changes in version 1.5.3                        

BUG FIXES

- Fixed the wrong layering in plot_cell_categories() and added legend
when layer=TRUE

                        Changes in version 1.5.2                        

NOTES

- plot_cell_percentages() now returns the plot.

                        Changes in version 1.5.1                        

BUG FIXES

- No error message when there are duplicates generated by the "TSNE"
option of dimensionality_reduction_plot() function.
- Deleted duplicated plots returned by plot_cell_categories().

                        Changes in version 1.5.0                        

Development version on Bioconductor 3.19.

[splatter](/packages/splatter)
--------

                 Changes in version 1.28.0 (2024-05-01)                 

- 
  Properly prevent output from BASiCSEstimate() and
  scDDSimulate() when verbose = FALSE

- 
  Minor adjustments to tests to set verbose = FALSE and specify
  expected warnings

- 
  Replace the package man page with an auto-generated page

- 
  Add GitHub usernames and the Bioconductor package URL to
  description

- 
  Add @keywords internal and @noRd to function docs where needed

- 
  Modernise code by restyling using styler and fix other minor
  issues suggested by BiocCheck

[spoon](/packages/spoon)
-----

                 Changes in version 0.99.0 (2024-01-18)                 

- version for submission to Bioconductor

[SpotSweeper](/packages/SpotSweeper)
-----------

                       Changes in version 0.99.1                        

- Added a NEWS.md file to track changes to the package.
- First full version of the package to be submitted to Bioconductor.
See Bioconductor submission here.

[SQLDataFrame](/packages/SQLDataFrame)
------------

                       Changes in version 1.17.1                        

- There is a total overhaul of the package implementation. It now
  defines a 1-dimension DelayedArray called SQLColumnVector to be the
  basic element which can sit inside a normal DataFrame for delayed
  operation on-disk. It also defines SQLDataFrame to represent a SQL
  table (with selected columns), which will be fully file-backed so no
  data is loaded into memory until requested. This allows uses to
  represent large datasets in limited memory. SQLDataFrame inherites
  from DataFrame so it can be used anywhere in Bioconductor's ecosystem
  that accepts DataFrame, such as SummarizedExperiment.

- SQLDataFrame is primarily useful for indicating that the in-memory
  representation is consistent with the underlying SQL file. It
  supports
  all the usual methods for a `DataFrame`, except that the data is kept
  on file and referenced as needed. However, only some operations can
  preserve the SQLDataFrame, such as column subsetting and some
  cbinding
  operations. Most operations that add or change data (e.g., filter for
  row subsetting, mutate for new/modify columns) would collapse to a
  regular DFrame with `SQLColumnVector`s before applying the
  operation. They are still file-backed but lack the guarantee of file
  consistency. The fallback to `DFrame` ensures that a `SQLDataFrame`
  is
  interoperable with other Bioconductor data structures that need to
  perform arbitrary `DataFrame` operations.

                 Changes in version 1.17.0 (2019-10-01)                 

NEW FEATURES

- Supports tidy grammar for 'select', 'filter', 'mutate' and '%>%'
  pipe.

- Supports representation and saving of MySQL database tables.

- Supports lazy cross-MySQL database table aggregations, such as join,
  union, rbind, etc.

BUG FIXES

- Fixed bugs for single square bracket subsetting with key column(s).

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

                       Changes in version 1.15.5                        

- update automated documentation generation

                       Changes in version 1.15.2                        

- add ellipsis to model methods

                       Changes in version 1.15.1                        

- remove duplicates when using param_ids for objects with inherited
  params

- set check.names = FALSE in DatasetExperiment

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.34.0                        

NEW FEATURES

- Add terminators() method, same as promoters() but for terminator
  regions.

- Add "Top-level dimnames vs assay-level dimnames" section to vignette.
  Also fix typo in man/SummarizedExperiment-class.Rd. This is in
  response
  to https://github.com/Bioconductor/SummarizedExperiment/issues/79

BUG FIXES

- Fix typo in man/SummarizedExperiment-class.Rd

[SurfR](/packages/SurfR)
-----

                Changes in version 0.99.54 (2024-01-14)                 

- updated vignette, added references.bib

                Changes in version 0.99.15 (2023-10-09)                 

- added README.md file

                Changes in version 0.99.13 (2023-09-26)                 

- updated ind_deg.Rd file

                Changes in version 0.99.12 (2023-09-26)                 

- removed deprecated rio import

                 Changes in version 0.99.7 (2023-09-26)                 

- Addressed notes

- Added BiocFileCache in Gene2SProtein

                 Changes in version 0.99.6 (2023-09-25)                 

- Addressed notes

                 Changes in version 0.99.5 (2023-09-21)                 

- Addressed notes

                 Changes in version 0.99.4 (2023-09-21)                 

- Addressed notes

                 Changes in version 0.99.3 (2023-09-11)                 

- Removed SurfR.Rproj file and updated gitignore

                 Changes in version 0.99.2 (2023-09-11)                 

- Removed SurfR.Rproj file and updated gitignore

                 Changes in version 0.99.1 (2023-09-11)                 

- Updated gitignore

                 Changes in version 0.99.0 (2023-08-09)                 

- Submitted to Bioconductor

[SVMDO](/packages/SVMDO)
-----

                 Changes in version 1.3.0 (2023-12-27)                  

- Application of 10-fold cross-validation with tune.svm function from
e1071 package

- Addition of BiocStyle in the Suggestion section of Description file

- Updating the package citation

- Updating the package citation

- Fixing merge-related file section problems occured during the
latest
pushing process

- Fixing CITATION-related warning

[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.15.3                        

- Addition of FastLabelOOM function to find communities in
graphs/networks on disk space.

                       Changes in version 1.15.2                        

- Addition of PrepareSeqs function, beginning the process of
deprecating PairSummaries in favor of more cohesive and user
friendly functions.

                       Changes in version 1.15.1                        

- Fixes bug in JRF distance causing scores to be higher than
expected.

[TADCompare](/packages/TADCompare)
----------

                 Changes in version 1.13.3 (2024-03-06)                 

- Remove rGREAT dependency

                 Changes in version 1.13.2 (2024-03-04)                 

- Move rGREAT to Imports

                 Changes in version 1.13.1 (2024-02-28)                 

- Fix build errors

[TargetSearch](/packages/TargetSearch)
------------

                        Changes in version 2.6.0                        

NEW FEATURES

- New options for `ImportFameSettings`: New parameter `standard`
  can take the RI values. The input file can be also a matrix.
  Column names can now be specified.

- A function to plot reference spectra for `tsLib` objects in a
  similar style of other visualization tools.

BUG FIXES

- Check for out-of-bounds in binary search.

DOCUMENTATION

- Add note of arguments passed to `plot` via the `dots` operator.
  The functions `ri_plot_peak` and `ncdf4_plot_peak` accept extra
  arguments passed to `plot` and others via the `...` argument.
  However, not all parameters, in particular `panel.first` and
  `panel.last`, will have an effect.

- Fix misplaced dot man page (`RIcorrrect`).

- Fix wrong braces in man page (`tsSample-class`).

INTERNAL

- Tidy up `NAMESPACE` file.

- Add tests for function `ImportFameSettings`.

- Refactor internal code of `ImportLibrary` and improve the spectrum
  parser algorithm (allow NA and make sure the format is valid).
  Add tests cases.

- Improve `tsLib` spectrum handling by ensuring the internal data
  is consistent. Also allow for empty spectra.

- Improve validity checks for `tsProfile` and `tsMSdata` objects.
  Add tests cases.

                        Changes in version 2.4.2                        

BUG FIXES

- Fix `maybe-uninitialized` warning in some platforms. Set the struct
  elements to `NULL`.

- Fix RI files format incompatibility (affected versions: 2.4.0 and
  .1).
  Starting from 2.4.0, these files allow for empty spectra data, but
  older version raise an error. Therefore, we enforce non empty
  spectra data to keep compatibility.

                        Changes in version 2.4.1                        

BUG FIXES

- An error is thrown if all parameters passed to `FindAllPeaks` are
  integers (a rare occurrence). Make sure that arguments are doubles
  before passing them to the C interface. Also add a test for it.

[TCC](/packages/TCC)
---

                       Changes in version 1.43.4                        

- deleted test code related with deprecated functions.

- deleted 'baySeq' dependency codes from vignette and other code.

- deleted 'baySeq' from Depends.

[TCGAutils](/packages/TCGAutils)
---------

                       Changes in version 1.24.0                        

Significant User-visible changes

- The legacy argument in ID translation functions (UUIDtoBarcode,
UUIDtoUUID, barcodeToUUID, and filenameToBarcode) has been defunct
and removed.

Bug fixes and minor improvements

- UUIDtoBarcode ensures that results are ordered based on the input
UUIDs.
- Include informative error message regarding translation of UUIDs
from legacy files.

[TENxIO](/packages/TENxIO)
------

                        Changes in version 1.6.0                        

Bug fixes and minor improvements

- Use rowData from filtered SingleCellExperiment object
(@MalteThodberg, #5)
- Add examples to TENxFileList documentation and improve constructor
function

[tidySpatialExperiment](/packages/tidySpatialExperiment)
---------------------

                       Changes in version 0.99.20                       

- tidySpatialExperiment is now in Bioconductor. Documentation has
been
updated accordingly.

                       Changes in version 0.99.19                       

- Improve code style.

                       Changes in version 0.99.18                       

- Version bump for Bioconductor.

                       Changes in version 0.99.17                       

- Version bump for Bioconductor.

                       Changes in version 0.99.16                       

- Version bump for Bioconductor.

                       Changes in version 0.99.15                       

- Bug fixes.

                       Changes in version 0.99.14                       

- Bug fixes.

                       Changes in version 0.99.13                       

- Removed unneeded dependencies.
- Updated README and vignette material.

                       Changes in version 0.99.12                       

- Moved pkgdown website deployment from docs directory to gh-pages
branch.
- Bug fixes.

                       Changes in version 0.99.11                       

- Added rworkflow github actions to support continuous and
collaborative development.

                       Changes in version 0.99.10                       

- Rectangular and elliptical gating added by Quentin Clayssen.

                       Changes in version 0.99.9                        

- tidySpatialExperiment now features a NEWS.md file!
- The package also now has its own pkgdown website.
- Updated README to better reflect the development of tidyomics
ecosystem.

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

[topdownr](/packages/topdownr)
--------

                        Changes in version 1.25                         

Changes in version 1.25.1

- New version for Bioc 3.19 (devel)
- Adapt to fragment functions that are moved from MSnbase to PSMatch.

Changes in version 1.19.2

[tpSVG](/packages/tpSVG)
-----

                       Changes in version 0.99.1                        

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- Your main changes to a function foo() or parameter param.

BUG FIXES

- Your bug fixes. See more details at
http://bioconductor.org/developers/package-guidelines/#news.

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.39.17                       

- Fix the legend with alpha for lolliplot.

                       Changes in version 1.39.16                       

- Import txdbmaker.

                       Changes in version 1.39.15                       

- Optimize the reduce function for loopBouquetPlot.

                       Changes in version 1.39.14                       

- reduce the input for loopBouquetPlot.

                       Changes in version 1.39.13                       

- Handle the empty input for loopBouquetPlot.

                       Changes in version 1.39.12                       

- Update the safe text algorithm for loopBouquetPlot.

                       Changes in version 1.39.11                       

- Add parameter safe_text_force for loopBouquetPlot.

                       Changes in version 1.39.10                       

- Add color for geneTrack.

                       Changes in version 1.39.9                        

- Fix multiple bugs in loopBouquetPlot.

                       Changes in version 1.39.8                        

- Change loopBouquetPlot to grid plot.

                       Changes in version 1.39.7                        

- Fix the example for loopBouquetPlot.

                       Changes in version 1.39.6                        

- add loopBouquetPlot.

                       Changes in version 1.39.5                        

- Browse heatmap by js (scratch).

                       Changes in version 1.39.4                        

- Fix the issue if the first label is empty for lolliplot.

                       Changes in version 1.39.3                        

- Fix the mismatch data track when change the ylim in browseTracks.

                       Changes in version 1.39.2                        

- Fix the wrong plot for multiple genes track in browseTracks.

                       Changes in version 1.39.1                        

- Fix the y-axis for lollipop plot.

[treeclimbR](/packages/treeclimbR)
----------

                       Changes in version 0.99.1                        

- Add explicit dependency on R
- Bug fix - make sure that ... are passed on in runDS
- Internal refactoring to reduce repeated code

                       Changes in version 0.99.0                        

- Prepare for Bioconductor submission

                        Changes in version 0.1.5                        

- Add NEWS file

[treeio](/packages/treeio)
------

                       Changes in version 1.27.1                        

- added support for multiple trees / writing phylo objects in
write.beast() (2024-04-08, Mon, #113)
- speed up read.beast() (2023-12-13, Wed, #118)
- optimize write.jtree() (2023-12-13, Wed, #117)
- write.jplace() method to export jplace object to a jplace file
(2023-11-27, Mon, #112, #115)

[tRNA](/packages/tRNA)
----

                 Changes in version 1.21.1 (2024-01-17)                 

- replaced deprecated ggplot2 functions

[txdbmaker](/packages/txdbmaker)
---------

                        Changes in version 1.0.0                        

- The txdbmaker package is a spin-off from the GenomicFeatures package.
  It includes the following functions, all of them originally from the
  GenomicFeatures package:
  - browseUCSCtrack
  - getChromInfoFromBiomart
  - makeFDbPackageFromUCSC
  - makeFeatureDbFromUCSC
  - makePackageName
  - makeTxDb
  - makeTxDbFromBiomart
  - makeTxDbFromEnsembl
  - makeTxDbFromGFF
  - makeTxDbFromGRanges
  - makeTxDbFromUCSC
  - makeTxDbPackage
  - makeTxDbPackageFromBiomart
  - makeTxDbPackageFromUCSC
  - supportedMiRBaseBuildValues
  - supportedUCSCFeatureDbTables
  - supportedUCSCFeatureDbTracks
  - supportedUCSCtables
  - UCSCFeatureDbTableSchema

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.21.4                        

- Changing language in docs to "digest" instead of "checksum".

                       Changes in version 1.21.3                        

- GENCODE 44 (H.s.), M34 (M.m), and Ensembl 111
- RefSeq p13 for human, p6 for mouse

[tximport](/packages/tximport)
--------

                       Changes in version 1.31.1                        

- Experimental support for oarfish. At this point
  we are importing number of reads as both
  'abundance' and 'counts'.

[UCell](/packages/UCell)
-----

                        Changes in version 2.7.6                        

- Changed factor for normalization of UCell scores, to account
  for minimal rank. For typical use cases the behaviour is
  similar to the previous implementation, but has an impact with
  very large signatures (UCell score distributions are now more
  homogeneous). See function UCell:::u_stat()) for details.

                        Changes in version 2.7.4                        

- Add support for Seurat v5 assay and datasets in multiple
  layers. Added dependency on Seurat >= 5.0.

- Change default for chunk.size to 100. On parallelized jobs,
  this can improve up to 2-fold execution time.

- Smaller test dataset (30 cells), to speed up package function
  checks.

[UCSC.utils](/packages/UCSC.utils)
----------

                        Changes in version 1.0.0                        

- First version of the package that is ready for general use.

[UniProt.ws](/packages/UniProt.ws)
----------

                       Changes in version 2.44.0                        

BUG FIXES AND MINOR IMPROVEMENTS

- No significant changes.

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.22.0                        

BUG FIXES

- Address C++ compiler warnings.

- Suppress warnings when using as.data.frame() on DataFrame objects
  in the SequenceSearches.Rmd vignette.

[updateObject](/packages/updateObject)
------------

                        Changes in version 1.8.0                        

- No changes in this version.

[UPDhmm](/packages/UPDhmm)
------

                       Changes in version 0.99.4                        

Final corrections and final editing of vignette

                       Changes in version 0.99.3                        

2nd correction of reviewer's suggestions

                       Changes in version 0.99.2                        

Correction of reviewer's suggestions

                       Changes in version 0.99.1                        

First pre-review corrections

                       Changes in version 0.99.0                        

First release of UPDhmm package

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

                       Changes in version 1.33.11                       

- Feb 7, 2024
- fix bug in dream(...,ddf="Kenward-Roger") that gave false positives
and negatives
- scaled weights properly to fix this issue, since df in
lmerTest::contest() assumes mean of weights is 1.
- old code used unscaled weights, so df values were too large

                       Changes in version 1.33.10                       

- Feb 5, 2024
- in augmentPriorCount() and voomWithDreamWeights(), add argument
scaledByLib=FALSE

                       Changes in version 1.33.9                        

- Jan 25, 2024
- bug fix in BIC() and .fitExtractVarPartModel()

                       Changes in version 1.33.8                        

- Jan 22, 2024
- add BIC() for result of dream() and lmFit()

                       Changes in version 1.33.7                        

- Jan 19, 2024
- add logLik() for result of dream() and lmFit()

                       Changes in version 1.33.6                        

- fix issue with subsetting using logicals

                       Changes in version 1.33.5                        

- Dec 21, 2023
- fix issue with ddf always calling "Satterthwaite"
- enforce subsetting of residuals in assign( "[.MArrayLM2",)
- this comes with when rdf < 1

                       Changes in version 1.33.4                        

- Dec 13, 2023
- fix bug with BiocParallel in Windows
- https://github.com/GabrielHoffman/variancePartition/issues/91
- handle case where a single contrasts fails in makeContrastsDream()

                       Changes in version 1.33.3                        

- Dec 5, 2023
- in residuals() when dividing by sqrt(1-hatvalues) add small offset
to make sure the value is positive

                       Changes in version 1.33.2                        

- Nov 13, 2023
- add augmentPriorCount()
- add prior.count argument to voomWithDreamWeights() and feed it to
augmentPriorCount()

                       Changes in version 1.33.1                        

- Nov 2, 2023
- move omp_set_num_threads() deeper in nesting

[velociraptor](/packages/velociraptor)
------------

                       Changes in version 1.13.1                        

- Robust fallback mechanism using basiliskRun option testload=.

[VisiumIO](/packages/VisiumIO)
--------

                        Changes in version 1.0.0                        

- Package released in Bioconductor!

[Voyager](/packages/Voyager)
-------

                        Changes in version 1.6.0                        

- Option to plot polygon geometry without fill in plotGeometry
- Option to show axes in both light and dark themes
- Plot images of classes BioFormatsImage and ExtImage
- Colorize image plot with different channels
- Assign different grayscale images to channels
- Added plotImage to only plot image without plotting geometry

[xcms](/packages/xcms)
----

                         Changes in version 4.1                         

Changes in version 4.1.13

- Add parameter rtimeDifferenceThreshold to ObiwarpParam allowing to
customize the threshold used by obiwarp to determine whether gaps
are present in the sequence of retention times of a sample. This
addresses/fixes issue #739.

Changes in version 4.1.12

- Implementation of the LamaParama class and method for the
adjustRtime() function. Allowing alignment of a dataset based on
landmarks (lamas) from an external reference dataset.
- Implementation of related user-level function
matchLamasChromPeaks(), summarizeMatchLama() and plot(LamaParama)
which allows for evaluation of matching between lamas and
chromPeaks.

Changes in version 4.1.11

- Clean up of required and suggested packages and namespace imports.
- Re-creation of bundled data objects.

Changes in version 4.1.10

- Ensure backward compatibility for parameter objects that gained
additional slots.

Changes in version 4.1.9

- Fix bug in filterFeatures,PercentMissingFilter.

Changes in version 4.1.8

- Fixing issue #716: edit of .empty_chrom_peaks function so an sn
column is returned. Fixes extracting and plotting of peaks after
using manualChromPeaks

Changes in version 4.1.7

- Implementation of filterFeatures function with filter parameters:
RsdFilter, DratioFilter, PercentMissingFilter, BlankFlag. They can
be used ot filter features from XcmsResult and SummarizedExperiment
objects.
- Addition of a section in the main xcms vignette to describe how to
use it.

Changes in version 4.1.6

- Import filterSpectra from MsExperiment.
- Import breaks_ppm from MsCoreUtils.
- Update featureArea function to consider all chromatographic peaks
per feature, not only the one with the highest intensity. As a
consequence, returned m/z and rt ranges might be higher which has an
influence in featureChromatograms, EIC-based feature grouping and,
to a lesser extent also in gap-filling. Related documentation was
updated.
- Improve performance of the featureArea function (and related of the
PeakAreaParam-based gap filling).
- Add parameter ppm to PeakDensityParam to enable peak-density-based
correspondence throgh m/z-dependent bins along the m/z.

Changes in version 4.1.5

- Improve performance of the chromatogram call for XcmsExperiment
objects.
- Remove internal (not exported) normalization functions. These have
been transferred to the MetaboCoreUtils package.
- Support subsetting of XcmsExperiment with negative indices.

Changes in version 4.1.4

- Rename variable data in the vignette to faahko.
- Fix issue in adjustRtime resulting in corrupt processHistory.
- Add support to perform peakGroups alignment using pre-defined
anchor
peak matrix (i.e., the numeric matrix with retention times of anchor
peaks in the samples that can be used to align these samples).
- Fix errors related to invalid Chromatogram objects extracted from
xcms results: ensure MS level in chromPeaksMatrix is integer.
- Fix definition of anchor peaks for peakGroups alignment with subset
(issue #702).
- Add filterMsLevel method for MsExperiment and XcmsExperiment.
- Ensure chunk-wise processing of Spectra (introduced with
version 1.13.2) is disabled when xcms is using its own chunk-wise
processing.

Changes in version 4.1.3

- Add parameter verboseBetaColumns to CentWaveParam to enable
calculation of additional peak quality metrics comparing the EIC to
an idealized bell curve.

Changes in version 4.1.2

- Add a param = to generic function storeResults: PlainTextParam to
save an XcmsExperiment or MsExperiment object as colleciton of plain
text files.

Changes in version 4.1.1

- Add method storeResults and one of its param =: RDataParam to save
an XcmsExperiment object as an .RData file.

[zellkonverter](/packages/zellkonverter)
-------------

                       Changes in version 1.14.0                        

Major changes

- 
  Add environment for *anndata* v0.10.6. This is now the default
  envrionment for the Python reader/writer.

Minor changes

- 
  Improve warnings when converting matrices fails

- 
  Minor change to writing DelayedArray matrices for compatibility
  with *HDF5Array* >= v1.31.1

Bug fixes

- 
  Correctly handle use_backed = TRUE with newer *anndata*
  versions

- 
  Correctly instantiate the *anndata* v0.10.2 environment

- 
  Minor fixes for typos etc.

[zenith](/packages/zenith)
------

                        Changes in version 1.5.3                        

- March 8, 2024
- add coef to zenithPR_gsa()

                        Changes in version 1.5.2                        

- Feb 2, 2024
- resolve issue with renameing gene sets when multiple coefs are used

                        Changes in version 1.5.1                        

- throw error when zenith is run on multiple coefs

[ZygosityPredictor](/packages/ZygosityPredictor)
-----------------

                        Changes in version 1.3.5                        

- fixed small issues in tool and documentation

                        Changes in version 1.3.2                        

- Updated tool to revised manuscript

- New phasing approaches via predefined haploblocks and allelic
  imbalance

- Confidence measurement for each gene status

NEWS from existing Data Experiment Packages
===================================


[CardinalWorkflows](/packages/CardinalWorkflows)
-----------------

                       Changes in version 1.35.5                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add startup message with available vignettes

                       Changes in version 1.35.4                        

SIGNIFICANT USER-VISIBLE CHANGES

- Added multiple instance learning to
  classification vignette

                       Changes in version 1.35.3                        

BUG FIXES

- Vignette edits

                       Changes in version 1.35.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Updated vignettes for Cardinal v3.6

                       Changes in version 1.35.1                        

NEW FEATURES

- Added 'exampleMSIData()' to load data

SIGNIFICANT USER-VISIBLE CHANGES

- Removed old analysis datasets

[curatedPCaData](/packages/curatedPCaData)
--------------

                 Changes in version 0.99.5 (2023-11-02)                 

- Bioconductor re-revised version

                 Changes in version 0.99.2 (2023-08-10)                 

- Bioconductor revised version

                 Changes in version 0.99.0 (2023-05-18)                 

- Revised whole package according to Bioconductor guidelines

- Package license is now CC-BY 4.0

[CytoMethIC](/packages/CytoMethIC)
----------

                       Changes in version 0.99.26                       

- Submission of CytoMethIC package.

- Fixed typos in documentation.

- Added cmi_checkVersion for user to ensure versions are coordinated

                       Changes in version 0.99.21                       

- Submission of CytoMethIC package.

- Fixed issues mentioned in second review.

[DMRcatedata](/packages/DMRcatedata)
-----------

                       Changes in version 2.20.2                        

- EPICv2 data added for vignette

                       Changes in version 2.20.1                        

- SNP data added for EPICv2

[gDRtestData](/packages/gDRtestData)
-----------

              Changes in version 22023-05-15 (2023-05-15)               

- fix related with data.table

               Changes in version 2024-02-28 (2024-02-28)               

- remove redundant annotation files

               Changes in version 2024-02-26 (2024-02-26)               

- improve pkgdown site

- improved references

- valid NEWS.md

               Changes in version 2024-02-14 (2024-02-14)               

- make documentation compatible with pkdgdown

               Changes in version 2024-02-07 (2024-02-07)               

- update test data as per new combo model

               Changes in version 2024-02-01 (2024-02-01)               

- standardize full response data functions

               Changes in version 2024-01-30 (2024-01-30)               

- update package vignette

               Changes in version 2024-01-22 (2024-01-22)               

- update testdata

- add new description fields

               Changes in version 2024-01-04 (2024-01-04)               

- update duplicated functions

               Changes in version 2023-11-22 (2023-11-22)               

- sync master with devel branch

               Changes in version 2023-10-24 (2023-10-24)               

- release Bioc 3.18

- prerelease Bioc 3.18

               Changes in version 2023-09-18 (2023-09-18)               


[imcdatasets](/packages/imcdatasets)
-----------

                 Changes in version 1.11.1 (2024-03-16)                 

- Removed defunct functions

[JohnsonKinaseData](/packages/JohnsonKinaseData)
-----------------

                       Changes in version 0.99.3                        

- Fix selection of left-most central acceptor

- Add documentation of phospho-peptide processing

                       Changes in version 0.99.2                        

- Fix misspelled reference in vignette

                       Changes in version 0.99.1                        

- Fix typos

                       Changes in version 0.99.0                        

- Initial version

[MouseAgingData](/packages/MouseAgingData)
--------------

                       Changes in version 0.99.13                       

- Replaced Parabiosis10x SingleCellExperiment object to include Cell
  Ontology labels and identifiers

- Submitted AgingBrain10x_2019NN SingleCellExperiment object for
  ExperimentHub

                       Changes in version 0.99.0                        

- Submitted Parabiosis10x SingleCellExperiment object for
  ExperimentHub

[rRDPData](/packages/rRDPData)
--------

                       Changes in version 1.23.1                        

Changes

- Update data to the RDP Classifier 2.14 released in August 2023
  contains the latest bacterial and archaeal taxonomy training set
  No. 19.

- The package now contains: 16srrna, fungalits_unite,
  fungalits_warcup, and fungallsu

[scpdata](/packages/scpdata)
-------

                       Changes in version 1.11                        

scpdata 1.11.4

- fix non-ASCII characters
- split data.R into multiple data specific files

scpdata 1.11.3

- leduc2022_pSCoPE: added cellenONE annotations
- guise2024: added dataset

scpdata 1.11.2

- khan2023: added dataset.

scpdata 1.11.1

- Added contribution vignette.
- gregoire2023_mixCTRL: added dataset.



[scRNAseq](/packages/scRNAseq)
--------

                       Changes in version 2.18.0                        

- Switched to the ArtifactDB representations for the underlying
  files. This uses language-agnostic formats (e.g., HDF5, JSON)
  instead of RDS files to store the various parts of each
  SingleCellExperiment. The user experience should be more or
  less the same as the datasets are indistinguishable once loaded
  into memory.

- Added the fetchDataset, to create SingleCellExperiment objects
  from the ArtifactDB file representations. This uses the
  alabaster.base package to do the loading, with some optional
  realization of the assays into memory. Advanced users can
  achieve faster loading times by keeping the assays as
  file-backed matrices.

- Introduced saveDataset and related functions to facilitate user
  uploads of their own datasets. This is accompanied by some
  step-by-step instructions in the vignette, plus some maintainer
  instructions in the README.

- Added searchDatasets to perform text searches on the metadata
  for each dataset, using the SQLite database compiled from the
  gypsum backend where the files are stored.

- Updated some datasets to reflect upstream changes (e.g., in
  ArrayExpress). Currently, this affects mostly
  SegerstolpePancreasData, where ArrayExpress decided to change
  the names and contents of various column annotations.

- Soft-deprecation of some redundant pieces of information in
  each dataset. Some examples are the column names of the
  Zilionis data, which were not unique and had no meaning; or the
  symbols of the rowData in the Segerstolpe data, which were
  redundant with the row names. These changes will only take
  effect when fetchDataset is used directly; the per-dataset
  getter functions have appropriate back-compatibility patches to
  restore this information.

- All getters now have a legacy= option to pull RDS files from
  ExperimentHub instead of the new formats from gypsum.

[seventyGeneData](/packages/seventyGeneData)
---------------

                       Changes in version 1.39.3                        

- Add data-raw directory for data generation scripts.

                       Changes in version 1.39.2                        

- Initial modernization of the package.

- New contributors added.

- Added a NEWS.md file to track changes to the package.

- Added a README.md file to describe the package.

- Update Vignette to use Rmarkdown and BiocStyle.

[SFEData](/packages/SFEData)
-------

                        Changes in version 1.5.3                        

- Added small subset of output files from Vizgen, CosMX, and Xenium to
  unit test read functions

[SingleCellMultiModal](/packages/SingleCellMultiModal)
--------------------

                       Changes in version 1.16.0                        

New features

- Added citation information to the package; see
  citation("SingleCellMultiModal") and the vignette.


[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.15.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- vis_gene() now has a multi_gene_method argument which provides 3
  methods for combining multiple continuous variables: z_score, pca,
  and sparsity. These options can now be used with run_app() (the
  interactive websites). These methods are further illustrated and
  documented in a new vignette available at
  https://research.libd.org/spatialLIBD/articles/multi_gene_plots.html.
  This work was contributed by @Nick-Eagles.

[SubcellularSpatialData](/packages/SubcellularSpatialData)
----------------------

                       Changes in version 0.99.0                        

- Create package

[TransOmicsData](/packages/TransOmicsData)
--------------

                       Changes in version 0.99.5                        

NEW FEATURES

- Added a `NEWS.md` file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- Made changes to vignette to improve usability.


NEWS from existing Workflows
===================================

No updated NEWS files.

Deprecated and Defunct Packages
===============================

**SOFTWARE:**

Forty six software packages were removed from this release (after being deprecated
in Bioc 3.18):

- BGmix, bigPint, biodbMirbase, BioMM, Clonality, COHCAP, CSSP, deco,
DeepBlueR, DMRforPairs, exomeCopy, fcoex, gaggle, GCSscore, genbankr, GISPA,
GOsummaries, GRridge, HPAStainR, imageHTS, IntOMICS, LineagePulse, logitT, LowMACA,
LPEadj, macat, mAPKL, mbOmic, Metab, MSstatsSampleSize, multiSight,
netbiov, OmicsLonDA, PFP, plethy, pwrEWAS, qrqc, Ringo, SCATE, SEPIRA,
seqbias, seqCNA, SISPA, snapCGH, sscore, Travel, trena

- Please note:  baySeq, biomvRCNS, MEIGOR, and RNAdecay, previously announced as
deprecated in 3.18, has been updated and remain in Bioconductor. Also previously
unannounced packages IntOMICS and COHCAP have been removed without a deprecation
cycle.

Seventy software packages are deprecated in this release and will be removed in Bioc 3.20:

- BDMMAcorrect, beadarraySNP, BHC, biodbLipidmaps, BioNetStat, CancerInSilico,
CancerSubtypes, cellHTS2, cliqueMS, CNVgears, compartmap, contiBAIT, CoRegNet,
CORREP, CoSIA, crisprseekplus, DNABarcodes, dpeak, EBSeqHMM, eegc, enrichTF,
ensemblVEP, exomePeak2, farms, FCBF, flowMap, FoldGO, FScanR, FunChIP, GOSim,
HumanTranscriptomeCompendium, ImmuneSpaceR, InterMineR, IntOMICS, IRISFGM,
iterClust, maigesPack, metagene, MetaVolcanoR, miRmine, MMAPPR2,
MobilityTransformR, multiOmicsViz, NeighborNet, NetPathMiner, oneSENSE,
openPrimeRui, pathVar, pcxn, PERFect, phemd, PloGO2, proteasy, PSEA, pwOmics,
RefPlus, ReQON, restfulSE, RIPAT, RLSeq, SimBindProfiles, SMAP, sparseDOSSA,
SpidermiR, SQUADD, StarBioTrek, STROMA4, TimiRGeN, TNBC.CMS, TnT

**EXPERIMENT DATA:** 

Nine experimental data packages were removed from this release (after being
deprecated in BioC 3.18):

- ccTutorial, ChIC.data, mAPKLData, MAQCsubsetILM, MIGSAdata, pwrEWAS.data,
SCATEData, seqCNA.annot, stjudem

- Please note: DLBCL, previously announced as deprecated in 3.18, has been updated
and remain in Bioconductor.

Five experimental data packages are deprecated in this release and will be
removed in Bioc 3.20:

- CoSIAdata, MMAPPR2data, pcxnData, restfulSEData, RLHub

**ANNOTATION DATA:** 

No annotation packages were removed from this release (after being deprecated
in Bioc 3.18).

One annotation package is deprecated in this release and will be removed in
Bioc 3.20.

- MafH5.gnomAD.v3.1.2.GRCh38

**WORKFLOWS:** 

No workflow packages were removed from this release (after being deprecated in
Bioc 3.18).

No workflow packages were deprecated in this release and will be removed in
3.20.

**BOOKS:**

No books were removed from this release (after being deprecated in
Bioc 3.18).

No books were deprecated in this release and will be removed in
3.20.
