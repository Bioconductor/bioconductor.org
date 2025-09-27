October 30, 2024

**Bioconductor:**

We are pleased to announce Bioconductor 3.20, consisting of
2289 software packages, 431 experiment data packages, 928 annotation
packages, 30 workflows and 5 books.

There are 54 new software packages, 5 new data experiment packages,
4 new annotation packages, no new workflows, no new books, and many updates and
improvements to existing packages.

Bioconductor 3.20 is compatible with R 4.4, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher, macOS arm64 and Linux
arm64. This release will also include updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.20](#getting-started-with-bioconductor-320)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [NEWS from existing workflows](#news-from-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.20
======================================

To update to or install Bioconductor 3.20

1. Install R 4.4. Bioconductor 3.20 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 54 new software packages in this release of Bioconductor.

- [ADAPT](/packages/ADAPT) ADAPT carries out differential abundance
  analysis for microbiome metagenomics data in phyloseq format. It
  has two innovations. One is to treat zero counts as left censored
  and use Tobit models for log count ratios. The other is an
  innovative way to find non-differentially abundant taxa as
  reference, then use the reference taxa to find the differentially
  abundant ones.

- [AnVILAz](/packages/AnVILAz) The AnVIL is a cloud computing
  resource developed in part by the National Human Genome Research
  Institute. The AnVILAz package supports end-users and developers
  using the AnVIL platform in the Azure cloud. The package provides a
  programmatic interface to AnVIL resources, including workspaces,
  notebooks, tables, and workflows. The package also provides
  utilities for managing resources, including copying files to and
  from Azure Blob Storage, and creating shared access signatures
  (SAS) for secure access to Azure resources.

- [AnVILBase](/packages/AnVILBase) Provides generic functions for
  interacting with the AnVIL ecosystem. Packages that use either GCP
  or Azure in AnVIL are built on top of AnVILBase. Extension packages
  will provide methods for interacting with other cloud providers.

- [AnVILGCP](/packages/AnVILGCP) The package provides a set of
  functions to interact with the Google Cloud Platform (GCP) services
  on the AnVIL platform. The package is designed to work with the
  AnVIL package. User-level interaction with this package should be
  minimal.

- [assorthead](/packages/assorthead) Vendors an assortment of useful
  header-only C++ libraries. Bioconductor packages can use these
  libraries in their own C++ code by LinkingTo this package without
  introducing any additional dependencies. The use of a central
  repository avoids duplicate vendoring of libraries across multiple
  R packages, and enables better coordination of version updates
  across cohorts of interdependent C++ libraries.

- [BioGA](/packages/BioGA) Genetic algorithm are a class of
  optimization algorithms inspired by the process of natural
  selection and genetics. This package allows users to analyze and
  optimize high throughput genomic data using genetic algorithms.
  The functions provided are implemented in C++ for improved speed
  and efficiency, with an easy-to-use interface for use within R.

- [broadSeq](/packages/broadSeq) This package helps user to do easily
  RNA-seq data analysis with multiple methods (usually which needs
  many different input formats). Here the user will provid the
  expression data as a SummarizedExperiment object and will get
  results from different methods. It will help user to quickly
  evaluate different methods.

- [CatsCradle](/packages/CatsCradle) This package addresses two broad
  areas.  It allows for in-depth analysis of spatial transcriptomic
  data by identifying tissue neighbourhoods.  These are contiguous
  regions of tissue surrounding individual cells.  'CatsCradle'
  allows for the categorisation of neighbourhoods by the cell types
  contained in them and the genes expressed in them.  In particular,
  it produces Seurat objects whose individual elements are
  neighbourhoods rather than cells.  In addition, it enables the
  categorisation and annotation of genes by producing Seurat objects
  whose elements are genes.

- [CleanUpRNAseq](/packages/CleanUpRNAseq) RNA-seq data generated by
  some library preparation methods, such as rRNA-depletion-based
  method and the SMART-seq method, might be contaminated by genomic
  DNA (gDNA), if DNase I disgestion is not performed properly during
  RNA preparation. CleanUpRNAseq is developed to check if RNA-seq
  data is suffered from gDNA contamination. If so, it can perform
  correction for gDNA contamination and reduce false discovery rate
  of differentially expressed genes.

- [DeepTarget](/packages/DeepTarget) This package predicts a drug’s
  primary target(s) or secondary target(s) by integrating large-scale
  genetic and drug screens from the Cancer Dependency Map project run
  by the Broad Institute. It further investigates whether the drug
  specifically targets the wild-type or mutated target forms. To show
  how to use this package in practice, we provided sample data along
  with step-by-step example.

- [DFplyr](/packages/DFplyr) Provides `dplyr` verbs (`mutate`,
  `select`, `filter`, etc...) supporting `S4Vectors::DataFrame`
  objects. Importantly, this is achieved without conversion to an
  intermediate `tibble`. Adds grouping infrastructure to `DataFrame`
  which is respected by the transformation verbs.

- [dominoSignal](/packages/dominoSignal) dominoSignal is a package
  developed to analyze cell signaling through ligand - receptor -
  transcription factor networks in scRNAseq data. It takes as input
  information transcriptomic data, requiring counts, z-scored counts,
  and cluster labels, as well as information on transcription factor
  activation (such as from SCENIC) and a database of ligand and
  receptor pairings (such as from CellPhoneDB). This package creates
  an object storing ligand - receptor - transcription factor linkages
  by cluster and provides several methods for exploring, summarizing,
  and visualizing the analysis.

- [DuplexDiscovereR](/packages/DuplexDiscovereR) DuplexDiscovereR is
  a package designed for analyzing data from RNA cross-linking and
  proximity ligation protocols such as SPLASH, PARIS, LIGR-seq, and
  others. DuplexDiscovereR accepts input in the form of chimerically
  or split-aligned reads. It includes procedures for alignment
  classification, filtering, and efficient clustering of individual
  chimeric reads into duplex groups (DGs). Once DGs are identified,
  the package predicts RNA duplex formation and their hybridization
  energies. Additional metrics, such as p-values for random ligation
  hypothesis or mean DG alignment scores, can be calculated to rank
  final set of RNA duplexes. Data from multiple experiments or
  replicates can be processed separately and further compared to
  check the reproducibility of the experimental method.

- [EnrichDO](/packages/EnrichDO) To implement disease ontology (DO)
  enrichment analysis, this package is designed and presents a double
  weighted model based on the latest annotations of the human genome
  with DO terms, by integrating the DO graph topology on a global
  scale. This package exhibits high accuracy that it can identify
  more specific DO terms, which alleviates the over enriched problem.
  The package includes various statistical models and visualization
  schemes for discovering the associations between genes and diseases
  from biological big data.

- [EpipwR](/packages/EpipwR) A quasi-simulation based approach to
  performing power analysis for EWAS (Epigenome-wide association
  studies) with continuous or binary outcomes. 'EpipwR' relies on
  empirical EWAS datasets to determine power at specific sample sizes
  while keeping computational cost low. EpipwR can be run with a
  variety of standard statistical tests, controlling for either a
  false discovery rate or a family-wise type I error rate.

- [funOmics](/packages/funOmics) The 'funOmics' package ggregates or
  summarizes omics data into higher level functional representations
  such as GO terms gene sets or KEGG metabolic pathways. The
  aggregated data matrix represents functional activity scores that
  facilitate the analysis of functional molecular sets while allowing
  to reduce dimensionality and provide easier and faster biological
  interpretations. Coordinated functional activity scores can be as
  informative as single molecules!

- [geomeTriD](/packages/geomeTriD) geomeTriD (Three Dimensional
  Geometry Package) create interactive 3D plots using the GL library
  with the 'three.js' visualization library (https://threejs.org) or
  the rgl library. In addition to creating interactive 3D plots, the
  application also generates simplified models in 2D. These 2D models
  provide a more straightforward visual representation, making it
  easier to analyze and interpret the data quickly. This
  functionality ensures that users have access to both detailed
  three-dimensional visualizations and more accessible
  two-dimensional views, catering to various analytical needs.

- [ggseqalign](/packages/ggseqalign) Simple visualizations of
  alignments of DNA or AA sequences as well as arbitrary strings.
  Compatible with Biostrings and ggplot2. The plots are fully
  customizable using ggplot2 modifiers such as theme().

- [HoloFoodR](/packages/HoloFoodR) Utility package to facilitate
  integration and analysis of EBI HoloFood data in R. This package
  streamlines access to the resource, allowing for direct loading of
  data into formats optimized for downstream analytics.

- [HuBMAPR](/packages/HuBMAPR) 'HuBMAP' provides an open, global
  bio-molecular atlas of the human body at the cellular level. The
  `datasets()`, `samples()`, `donors()`, `publications()`, and
  `collections()` functions retrieves the information for each of
  these entity types. `*_details()` are available for individual
  entries of each entity type. `*_derived()` are available for
  retrieving derived datasets or samples for individual entries of
  each entity type. Data files can be accessed using
  `files_globus_url()`.

- [immApex](/packages/immApex) A set of tools to build
  tensorflow/keras-based models in R from amino acid and nucleotide
  sequences focusing on adaptive immune receptors. The package
  includes pre-processing of sequences, unifying gene nomenclature
  usage, encoding sequences, and combining models. This package will
  serve as the basis of future immune receptor sequence
  functions/packages/models compatible with the scRepertoire
  ecosystem.

- [immunogenViewer](/packages/immunogenViewer) Plots protein
  properties and visualizes position of peptide immunogens within
  protein sequence. Allows evaluation of immunogens based on
  structural and functional annotations to infer suitability for
  antibody-based methods aiming to detect native proteins.

- [iSEEtree](/packages/iSEEtree) iSEEtree is an extension of iSEE for
  the TreeSummarizedExperiment. It leverages the functionality from
  the miaViz package for microbiome data visualisation to create
  panels that are specific for TreeSummarizedExperiment objects. Not
  surprisingly, it also depends on the generic panels from iSEE.

- [kmcut](/packages/kmcut) The purpose of the package is to identify
  prognostic biomarkers and an optimal numeric cutoff for each
  biomarker that can be used to stratify a group of test subjects
  (samples) into two sub-groups with significantly different survival
  (better vs. worse). The package was developed for the analysis of
  gene expression data, such as RNA-seq. However, it can be used with
  any quantitative variable that has a sufficiently large proportion
  of unique values.

- [koinar](/packages/koinar) A client to simplify fetching
  predictions from the Koina web service. Koina is a model repository
  enabling the remote execution of models. Predictions are generated
  as a response to HTTP/S requests, the standard protocol used for
  nearly all web traffic.

- [MetMashR](/packages/MetMashR) A package to merge, filter sort,
  organise and otherwise mash together metabolite annotation tables.
  Metabolite annotations can be imported from multiple sources
  (software) and combined using workflow steps based on S4 class
  templates derived from the `struct` package. Other modular workflow
  steps such as filtering, merging, splitting, normalisation and
  rest-api queries are included.

- [MOSClip](/packages/MOSClip) Topological pathway analysis tool able
  to integrate multi-omics data. It finds survival-associated modules
  or significant modules for two-class analysis. This tool have two
  main methods: pathway tests and module tests. The latter method
  allows the user to dig inside the pathways itself.

- [MPAC](/packages/MPAC) Multi-omic Pathway Analysis of Cancer
  (MPAC), integrates multi-omic data for understanding cancer
  mechanisms. It predicts novel patient groups with distinct pathway
  profiles as well as identifying key pathway proteins with potential
  clinical associations. From CNA and RNA-seq data, it determines
  genes’ DNA and RNA states (i.e., repressed, normal, or activated),
  which serve as the input for PARADIGM to calculate Inferred Pathway
  Levels (IPLs). It also permutes DNA and RNA states to create a
  background distribution to filter IPLs as a way to remove events
  observed by chance. It provides multiple methods for downstream
  analysis and visualization.

- [MsBackendMetaboLights](/packages/MsBackendMetaboLights)
  MetaboLights is one of the main public repositories for storage of
  metabolomics experiments, which includes analysis results as well
  as raw data. The MsBackendMetaboLights package provides
  functionality to retrieve and represent mass spectrometry (MS) data
  from MetaboLights. Data files are downloaded and cached locally
  avoiding repetitive downloads. MS data from metabolomics
  experiments can thus be directly and seamlessly integrated into
  R-based analysis workflows with the Spectra and
  MsBackendMetaboLights package.

- [OmicsMLRepoR](/packages/OmicsMLRepoR) This package provides
  functions to browse the harmonized metadata for large omics
  databases. This package also supports data navigation if the
  metadata incorporates ontology.

- [omXplore](/packages/omXplore) This package contains a collection
  of functions (written as shiny modules) for the visualisation and
  the statistical analysis of omics data. These plots can be
  displayed individually or embedded in a global Shiny module.
  Additionaly, it is possible to integrate third party modules to the
  main interface of the package omXplore.

- [PepSetTest](/packages/PepSetTest) Peptide Set Test (PepSetTest) is
  a peptide-centric strategy to infer differentially expressed
  proteins in LC-MS/MS proteomics data. This test detects coordinated
  changes in the expression of peptides originating from the same
  protein and compares these changes against the rest of the
  peptidome. Compared to traditional aggregation-based approaches,
  the peptide set test demonstrates improved statistical power, yet
  controlling the Type I error rate correctly in most cases. This
  test can be valuable for discovering novel biomarkers and
  prioritizing drug targets, especially when the direct application
  of statistical analysis to protein data fails to provide
  substantial insights.

- [Pirat](/packages/Pirat) Pirat enables the imputation of missing
  values (either MNARs or MCARs) in bottom-up LC-MS/MS proteomics
  data using a penalized maximum likelihood strategy. It does not
  require any parameter tuning, it models the instrument censorship
  from the data available. It accounts for sibling peptides
  correlations and it can leverage complementary transcriptomics
  measurements.

- [plyxp](/packages/plyxp) The package provides `rlang` data masks
  for the SummarizedExperiment class. The enables the evaluation of
  unquoted expression in different contexts of the
  SummarizedExperiment object with optional access to other contexts.
  The goal for `plyxp` is for evaluation to feel like a data.frame
  object without ever needing to unwind to a rectangular data.frame.

- [PolySTest](/packages/PolySTest) The complexity of high-throughput
  quantitative omics experiments often leads to low replicates
  numbers and many missing values. We implemented a new test to
  simultaneously consider missing values and quantitative changes,
  which we combined with well-performing statistical tests for high
  confidence detection of differentially regulated features. The
  package contains functions to run the test and to visualize the
  results.

- [PRONE](/packages/PRONE) High-throughput omics data are often
  affected by systematic biases introduced throughout all the steps
  of a clinical study, from sample collection to quantification.
  Normalization methods aim to adjust for these biases to make the
  actual biological signal more prominent. However, selecting an
  appropriate normalization method is challenging due to the wide
  range of available approaches. Therefore, a comparative evaluation
  of unnormalized and normalized data is essential in identifying an
  appropriate normalization strategy for a specific data set. This R
  package provides different functions for preprocessing,
  normalizing, and evaluating different normalization approaches.
  Furthermore, normalization methods can be evaluated on downstream
  steps, such as differential expression analysis and statistical
  enrichment analysis. Spike-in data sets with known ground truth and
  real-world data sets of biological experiments acquired by either
  tandem mass tag (TMT) or label-free quantification (LFQ) can be
  analyzed.

- [rhinotypeR](/packages/rhinotypeR) "rhinotypeR" is designed to
  automate the comparison of sequence data against prototype strains,
  streamlining the genotype assignment process. By implementing
  predefined pairwise distance thresholds, this package makes
  genotype assignment accessible to researchers and public health
  professionals. This tool enhances our epidemiological toolkit by
  enabling more efficient surveillance and analysis of rhinoviruses
  (RVs) and other viral pathogens with complex genomic landscapes.
  Additionally, "rhinotypeR" supports comprehensive visualization and
  analysis of single nucleotide polymorphisms (SNPs) and amino acid
  substitutions, facilitating in-depth genetic and evolutionary
  studies.

- [scDiagnostics](/packages/scDiagnostics) The scDiagnostics package
  provides diagnostic plots to assess the quality of cell type
  assignments from single cell gene expression profiles. The
  implemented functionality allows to assess the reliability of cell
  type annotations, investigate gene expression patterns, and explore
  relationships between different cell types in query and reference
  datasets allowing users to detect potential misalignments between
  reference and query datasets. The package also provides
  visualization capabilities for diagnostics purposes.

- [scDotPlot](/packages/scDotPlot) Dot plots of single-cell RNA-seq
  data allow for an examination of the relationships between cell
  groupings (e.g. clusters) and marker gene expression. The scDotPlot
  package offers a unified approach to perform a hierarchical
  clustering analysis and add annotations to the columns and/or rows
  of a scRNA-seq dot plot. It works with SingleCellExperiment and
  Seurat objects as well as data frames.

- [scoup](/packages/scoup) An elaborate molecular evolutionary
  framework that facilitates straightforward simulation of codon
  genetic sequences subjected to different degrees and/or patterns of
  Darwinian selection. The model was built upon the fitness landscape
  paradigm of Sewall Wright, as popularised by the mutation-selection
  model of Halpern and Bruno. This enabled realistic evolutionary
  process of living organisms to be reproduced seamlessly. For
  example, an Ornstein-Uhlenbeck fitness update algorithm is
  incorporated herein. Consequently, otherwise complex biological
  processes, such as the effect of the interplay between genetic
  drift and mutation on the inference of diversifying selection, may
  now be investigated with minimal effort. Frequency-dependent and
  deterministic fitness landscape update techniques are also
  available.

- [scrapper](/packages/scrapper) Implements R bindings to C++ code
  for analyzing single-cell (expression) data, mostly from various
  libscran libraries. Each function performs an individual step in
  the single-cell analysis workflow, ranging from quality control to
  clustering and marker detection. It is mostly intended for other
  Bioconductor package developers to build more user-friendly
  end-to-end workflows.

- [seahtrue](/packages/seahtrue) Seahtrue organizes oxygen
  consumption and extracellular acidification analysis data from
  experiments performed on an XF analyzer into structured nested
  tibbles.This allows for detailed processing of raw data and
  advanced data visualization and statistics. Seahtrue introduces an
  open and reproducible way to analyze these XF experiments. It uses
  file paths to .xlsx files. These .xlsx files are supplied by the
  userand are generated by the user in the Wave software from Agilent
  from the assay result files (.asyr). The .xlsx file contains
  different sheets of important data for the experiment; 1. Assay
  Information - Details about how the experiment was set up. 2. Rate
  Data - Information about the OCR and ECAR rates. 3. Raw Data - The
  original raw data collected during the experiment. 4. Calibration
  Data - Data related to calibrating the instrument. Seahtrue focuses
  on getting the specific data needed for analysis. Once this data is
  extracted, it is prepared for calculations through preprocessing.
  To make sure everything is accurate, both the initial data and the
  preprocessed data go through thorough checks.

- [SpaNorm](/packages/SpaNorm) This package implements the spatially
  aware library size normalisation algorithm, SpaNorm. SpaNorm
  normalises out library size effects while retaining biology through
  the modelling of smooth functions for each effect. Normalisation is
  performed in a gene- and cell-/spot- specific manner, yielding
  library size adjusted data.

- [spatialSimGP](/packages/spatialSimGP) This packages simulates
  spatial transcriptomics data with the mean- variance relationship
  using a Gaussian Process model per gene.

- [SpectraQL](/packages/SpectraQL) The Mass Spec Query Language
  (MassQL) is a domain-specific language enabling to express a query
  and retrieve mass spectrometry (MS) data in a more natural and
  understandable way for MS users. It is inspired by SQL and is by
  design programming language agnostic. The SpectraQL package adds
  support for the MassQL query language to R, in particular to MS
  data represented by Spectra objects. Users can thus apply MassQL
  expressions to analyze and retrieve specific data from Spectra
  objects.

- [squallms](/packages/squallms) squallms is a Bioconductor R package
  that implements a "semi-labeled" approach to untargeted mass
  spectrometry data. It pulls in raw data from mass-spec files to
  calculate several metrics that are then used to label MS features
  in bulk as high or low quality. These metrics of peak quality are
  then passed to a simple logistic model that produces a
  fully-labeled dataset suitable for downstream analysis.

- [StabMap](/packages/StabMap) StabMap performs single cell mosaic
  data integration by first building a mosaic data topology, and for
  each reference dataset, traverses the topology to project and
  predict data onto a common embedding. Mosaic data should be
  provided in a list format, with all relevant features included in
  the data matrices within each list object. The output of stabMap is
  a joint low-dimensional embedding taking into account all available
  relevant features. Expression imputation can also be performed
  using the StabMap embedding and any of the original data matrices
  for given reference and query cell lists.

- [survClust](/packages/survClust) survClust is an outcome weighted
  integrative clustering algorithm used to classify multi-omic
  samples on their available time to event information. The resulting
  clusters are cross-validated to avoid over overfitting and output
  classification of samples that are molecularly distinct and
  clinically meaningful. It takes in binary (mutation) as well as
  continuous data (other omic types).

- [tidyFlowCore](/packages/tidyFlowCore) tidyFlowCore bridges the gap
  between flow cytometry analysis using the flowCore Bioconductor
  package and the tidy data principles advocated by the tidyverse. It
  provides a suite of dplyr-, ggplot2-, and tidyr-like verbs
  specifically designed for working with flowFrame and flowSet
  objects as if they were tibbles; however, your data remain flowCore
  data structures under this layer of abstraction. tidyFlowCore
  enables intuitive and streamlined analysis workflows that can
  leverage both the Bioconductor and tidyverse ecosystems for
  cytometry data.

- [tidysbml](/packages/tidysbml) Starting from one SBML file, it
  extracts information from each listOfCompartments, listOfSpecies
  and listOfReactions element by saving them into data frames. Each
  table provides one row for each entity (i.e. either compartment,
  species, reaction or speciesReference) and one set of columns for
  the attributes, one column for the content of the 'notes'
  subelement and one set of columns for the content of the
  'annotation' subelement.

- [tidytof](/packages/tidytof) This package implements an
  interactive, scientific analysis pipeline for high-dimensional
  cytometry data built using tidy data principles. It is specifically
  designed to play well with both the tidyverse and Bioconductor
  software ecosystems, with functionality for reading/writing data
  files, data cleaning, preprocessing, clustering, visualization,
  modeling, and other quality-of-life functions. tidytof implements a
  "grammar" of high-dimensional cytometry data analysis.

- [TMSig](/packages/TMSig) The TMSig package contains tools to
  prepare, analyze, and visualize named lists of sets, with an
  emphasis on molecular signatures (such as gene or kinase sets). It
  includes fast, memory efficient functions to construct sparse
  incidence and similarity matrices and filter, cluster, invert, and
  decompose sets. Additionally, bubble heatmaps can be created to
  visualize the results of any differential or molecular signatures
  analysis.

- [xenLite](/packages/xenLite) Define a relatively light class for
  managing Xenium data using Bioconductor.  Address use of parquet
  for coordinates, SpatialExperiment for assay and sample data.
  Address serialization and use of cloud storage.

- [zitools](/packages/zitools) zitools allows for zero inflated count
  data analysis by either using down-weighting of excess zeros or by
  replacing an appropriate proportion of excess zeros with NA.
  Through overloading frequently used statistical functions (such as
  mean, median, standard deviation), plotting functions (such as
  boxplots or heatmap) or differential abundance tests, it allows a
  wide range of downstream analyses for zero-inflated data in a less
  biased manner. This becomes applicable in the context of microbiome
  analyses, where the data is often overdispersed and zero-inflated,
  therefore making data analysis extremly challenging.

New Data Experiment Packages
=====================

There are 5 new data experiment packages in this release of Bioconductor.

- [bugphyzz](/packages/bugphyzz) bugphyzz is an electronic database
  of standardized microbial annotations. It facilitates the creation
  of microbial signatures based on shared attributes, which are
  utilized for bug set enrichment analysis. The data also includes
  annotations imputed with ancestra state reconstruction methods.

- [eoPredData](/packages/eoPredData) Provides access to eoPred
  pretrained model hosted on ExperimentHub. Model was trained on
  placental DNA methylation preeclampsia samples using mixOmics
  splsda. There are two resources: 1. the model object, and 2. a
  testing data set used to demonstrate the function.

- [EpipwR.data](/packages/EpipwR.data) This package provides
  reference data for EpipwR. EpipwR is a fast and efficient power
  analysis for continuous and binary phenotypes of epigenomic-wide
  association studies. This package is only meant to be used in
  conjunction with EpipwR.

- [LegATo](/packages/LegATo) LegATo is a suite of open-source
  software tools for longitudinal microbiome analysis. It is
  extendable to several different study forms with optimal
  ease-of-use for researchers. Microbiome time-series data presents
  distinct challenges including complex covariate dependencies and
  variety of longitudinal study designs. This toolkit will allow
  researchers to determine which microbial taxa are affected over
  time by perturbations such as onset of disease or lifestyle
  choices, and to predict the effects of these perturbations over
  time, including changes in composition or stability of commensal
  bacteria.

- [ProteinGymR](/packages/ProteinGymR) The ProteinGymR package
  provides analysis-ready data resources from ProteinGym, generated
  by Notin et al., 2023. ProteinGym comprises a collection of
  benchmarks for evaluating the performance of models predicting the
  effect of point mutations. This package provides access to 1. Deep
  mutational scanning (DMS) scores from 217 assays measuring the
  impact of all possible amino acid substitutions across 186
  proteins, 2. AlphaMissense pathogenicity scores for ~1.6 M
  substitutions in the ProteinGym DMS data, and 3. five performance
  metrics for 62 variant prediction models in a zero-shot setting.

New Annotation Packages
=====================

There are 4 new annotation packages.

- [HDO.db](/packages/HDO.db) A set of annotation maps describing the entire Human Disease 
    Ontology assembled using data from DO.
    Its annotation data comes from 
    https://github.com/DiseaseOntology/HumanDiseaseOntology/tree/main/src/ontology

- [IlluminaHumanMethylationMSAanno.ilm10a1.hg38](/packages/IlluminaHumanMethylationMSAanno.ilm10a1.hg38)
  An annotation package for Illumina's MSA methylation arrays.

- [IlluminaHumanMethylationMSAmanifest](/packages/IlluminaHumanMethylationMSAmanifest)
  A manifest package for use with Illumina's MSA methylation arrays, compatible
  with minfi.
  
- [TENET.AnnotationHub](/packages/TENET.AnnotationHub) AnnotationHub
  package containing datasets for use in the TENET package. Includes
  GenomicRanges objects representing putative enhancer, promoter, and
  open chromatin regions. All included datasets are aligned to the
  hg38 human genome.
  

New Workflow Packages
=====================

There are no new workflow packages in this release of Bioconductor.

New Online Books
=====================

There are no new books in this release of Bioconductor.

NEWS from existing Software Packages
===================================

[alabaster.base](/packages/alabaster.base)
--------------

                        Changes in version 1.6.0                        

- Distinguish between scalars and length-1 vectors when
  saving/loading lists. This effectively unboxes all length-1
  vectors in a list, by default; this is probably the more
  reasonable expectation for other languages that have a concept
  of scalars. Users can override this by calling I() on elements
  that they want to keep as length-1 vectors, in the same manner
  as jsonlite.

- Streamlined the definition of the Rfc3339 class so that it
  behaves better with I().

- Normalize paths to resolve ~ prior to calling C++ code.

- Open HDF5 files in read-only mode to avoid permission-related
  problems for readObject()-dispatched functions.

- Store numbers at maximum precision when saving lists in the
  JSON format via saveObject().

- Added registerValidateObjectSatisfiesInterface() and
  registerValidateObjectDerivedFrom(), to allow developers to
  declare that custom subclasses satisfy an interface or have an
  inheritance relationship, respectively.

- Updated validateDirectory() so that it works with a directory
  of objects saved via saveObject(). Objects saved under the old
  regime (i.e., stageObject()) are auto-detected but can also be
  explicitly validated by setting legacy=FALSE.

- Added a data.frame method for saveObject(), to avoid fallback
  to the list method.

[alabaster.mae](/packages/alabaster.mae)
-------------

                        Changes in version 1.6.0                        

- Respect application-level overrides when saving child
  components of a MultiAssayExperiment.

[alabaster.matrix](/packages/alabaster.matrix)
----------------

                        Changes in version 1.6.0                        

- Support the SVT_SparseMatrix version 1 class definition in
  saveObject(). However, note that this was not implemented for
  the soft-deprecated writeSparseMatrix(), which now errors if
  such objects are passed in.

- Added a extract_sparse_array() method for the WrapperArraySeed
  class, for some future-proofing when the seeds eventually make
  the switch.

- Bugfix for integer overflow when saving large sparse matrices
  in saveObject().

- Open all HDF5 files in read-only mode for readObject() dispatch
  functions, to avoid permission-related issues.

- Added altReloadDelayedObject(), altStoreDelayedObject(), and
  their associated getters/setters, to allow applications to
  override the delayed operation saving/reading process.

- Added registerReloadDelayedObjectFunction() to allow extension
  developers to register reader functions for new classes.

- Added a ReloadedArray.reuse.files="relsymlink" option in the
  saveObject() method for ReloadedArrays. This creates relative
  symbolic links to the original array files, which is more
  robust to their movement provided the linked files are moved in
  the same manner.

- Enable deduplication of identical seeds across multiple calls
  to storeDelayedObject() within a single "session". This avoids
  making multiple copies of the same seed for different
  DelayedArray instances with the same seeds, e.g., in a
  SummarizedExperiment.

- Added an external.save.args= option to storeDelayedObject() to
  avoid conflicts in the method arguments of saveObject().

[alabaster.ranges](/packages/alabaster.ranges)
----------------

                        Changes in version 1.6.0                        

- Open HDF5 files in read-only mode to avoid permission-related
  problems for readObject()-dispatched functions.

- Use altSaveObject() when saving child components inside
  saveObject() methods.

[alabaster.se](/packages/alabaster.se)
------------

                        Changes in version 1.6.0                        

- Preserve row names when setting the rowRanges in
  readRangedSummarizedExperiment.

[alabaster.spatial](/packages/alabaster.spatial)
-----------------

                        Changes in version 1.6.0                        

- Dispatch to saveObject() for any of the SpatialExperiment's
  images that implement a dedicated method for it.

[AlphaMissenseR](/packages/AlphaMissenseR)
--------------

                        Changes in version 1.2.0                        

- (v. 1.1.9) Manage duckdb connections more completely; all
registered
connections are disconnected (invalidated) whenever a new table is
created.
- (v. 1.1.7) Add 'Benchmarking with ProteinGym' vignette to benchmark
AlphaMissense predictions.
https://github.com/mtmorgan/AlphaMissenseR/pull/8. Thank you
@tram-nguyen-n
- (v. 1.1.6) Add gosling_plot() for visualizing variants as bar or
lollipop plots. Merges
https://github.com/mtmorgan/AlphaMissenseR/pull/6 Thanks @lee-t.
- (v. 1.1.4) Add clinvar_data() and clinvar_plot() for visualizing
ClinVar data. Merges
https://github.com/mtmorgan/AlphaMissenseR/pull/4 Thanks
@tram-nguyen-n
- (v. 1.1.2) af_predictions() returns a tibble with 21 columns,
instead of 20. Merges
https://github.com/mtmorgan/AlphaMissenseR/pull/3.

[AlpsNMR](/packages/AlpsNMR)
-------

                 Changes in version 4.7.2 (2024-08-10)                  

- Disable nested parallellization in nmr_detect_peaks_tune_snr().

                 Changes in version 4.7.1 (2024-06-02)                  

- Added nmr_autophase() for automated phase correction using the
NMRphasing package (#68).
- Added to_ASICS function to export dataset for ASICS quantification
(#68).

[AnnotationForge](/packages/AnnotationForge)
---------------

                       Changes in version 1.48.0                        

BUG FIXES AND MINOR IMPROVEMENTS

- Updated `viableIDs.rda` file in `extdata` directory (@lshep)

- Added `rmarkdown` to `Suggests` in `DESCRIPTION` file

[AnnotationHub](/packages/AnnotationHub)
-------------

                       Changes in version 3.13.0                        

BUG CORRECTION

- (3.13.1) merged @votti PR for not checking connection if localhub is
  selected

NEW FEATURES

- (3.13.2) Add auto clean of cache if corrupt index file

- (3.13.3) Add auto redownload once if resource is not loading in R, if
  localHub
  is not selected and internet available.

[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.18.0                        

USER VISIBLE CHANGES

- (v 1.17.18) Added has_avworkspace function to check for the
existence of an AnVIL workspace environment.

- (v 1.17.10) Internal functions now use AnVILGCP for gcloud
utilities.

- (v 1.17.8) Functions that use gcloud utilities are deprecated and
will be moved to AnVILGCP. See help(package = "AnVIL") for a
complete list. Documentation pages have a *-deprecated suffix.

- (v 1.17.3) Added Terra Data Repository (TDR) service as TDR(). See
service at https://data.terra.bio.

- (v 1.17.1) Gen3 services, avworkflow*_configuration() functions,
install(), repository(), and repositories() have been removed.

- (v 1.17.1) Defunct repository_stats function in favor of
BiocPkgTools::repositoryStats (@LiNk-NY)

BUG FIXES AND MINOR IMPROVEMENTS

- (v 1.17.20) Use lifeCycle from BiocBaseUtils to mark functions as
deprecated or defunct.

- (v 1.17.19) Increase robustness of gcloud_exists by testing gcloud
with the version command.

- (v 1.17.18) Remove mentions of AnVIL::install from the vignette.

- (v 1.17.13) Update to changes in rapiclient and use native pipe
operator.

- (v 1.17.7) Do not evaluate vignette chunks if gcloud_exists() is
FALSE

- (v 1.17.6) Update Dockstore API file, version, and URL

- (v 1.17.2) Use application/json as default Content-Type.

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.29.1                        

- Fix the space inserted at the beginning after formatC for
exportBamFile.

[basilisk](/packages/basilisk)
--------

                       Changes in version 1.18.0                        

- Switch to the latest Miniforge installer (24.3.0-0) by default.
  This is preconfigured to use the conda-forge channel and avoids
  issues with the non-FOSS licensing of the Anaconda
  repositories. Users can switch back to the old Miniconda
  installer by setting the BASILISK_USE_MINIFORGE=0 environment
  variable, but this will likely be deprecated in the next
  release.

- Channel specifications in the user's .condarc are now ignored.
  Only the channels by the developer in setupBasiliskEnv() (e.g.,
  via BasiliskEnvironment) will be respected, to improve the
  consistency of the constructed environments across devices.

[basilisk.utils](/packages/basilisk.utils)
--------------

                       Changes in version 1.18.0                        

- Switch to the latest Miniforge installer (24.3.0-0) by default.
  This is preconfigured to use the conda-forge channel and avoids
  issues with the non-FOSS licensing of the Anaconda
  repositories. Users can switch back to the old Miniconda
  installer by setting the BASILISK_USE_MINIFORGE=0 environment
  variable, but this will likely be deprecated in the next
  release.

- Update the reticulate version in the fallback environment to
  1.38.

[BatchQC](/packages/BatchQC)
-------

                        Changes in version 2.1.6                        

Bug Fixes

- Corrected code for proper division of less than or 20+ samples

#Version 2.1.5

bug Fixes

- Coerce variable of interest to a factor

                        Changes in version 2.1.4                        

Major Changes

- Added negative binomial check for 20+ samples to DESeq2

                        Changes in version 2.1.3                        

Major Changes

- Added negative binomial check for less than 20 samples to DESeq2

                        Changes in version 2.1.2                        

Major Changes

- Added Variation Ratio Statistic to the explained variation tab

Minor Changes

- Removed extra "Samples" column from example data
- Uploaded bladder example data batch variable as a factor

                        Changes in version 2.1.1                        

Bug Fixes

- Updated imports to include shinyjs
- Updated imports to remove dendextend which is no longer utilized
- Corrected typos in Intro vignette

[BayesSpace](/packages/BayesSpace)
----------

                       Changes in version 1.15.3                        

Minor improvements and fixes

- Minor improvement of code style

                       Changes in version 1.15.2                        

Minor improvements and fixes

- Minor bugfixes related to R install, build, check
- Documentation improvements

                       Changes in version 1.15.1                        

Major updates

- Accelerate resolution enhancement with multithreaded spatialEnhance
- Improve mixing of MCMC for spatialEnhance with adaptive MCMC
- Support VisiumHD

Minor improvements and fixes

- Support SpaceRanger v2.0+
- Find the optimal number of cores with coreTune before enhancing
resolution
- Adjust the proportion of samples to be removed during burnin using
adjustClusterLabels
- Faster neighbor finding for subspots
- Customize resolution enhancement for VisiumHD data

                       Changes in version 1.15.0                        

New Bioconductor devel (3.20)

- Version numbering change with Bioconductor version bump

[beachmat](/packages/beachmat)
--------

                       Changes in version 2.22.0                        

- Moved all C++ libraries to the assorthead package. This uses
  the latest versions of the tatami framework.

- Automatically attempt to use beachmat.hdf5 (if installed) when
  encountering the various HDF5Array classes.

- Minor optimization when dealing with some known no-op matrices,
  e.g., WrapperArraySeed objects from alabaster.matrix.

- Map delayed type coercions from type<- to their corresponding
  native representations (or no-ops).

- Bugfix for initializeCpp() to work properly with dense Matrix
  instances.

- Soft-deprecated whichNonZero() in favor of SparseArray's new
  nzwhich() and nzvals() functions.

- Exported tatami.* utilities for manipulating
  already-initialized pointers from initializeCpp(). This
  includes all of the delayed operation wrappers, some previously
  internal functions for extracting matrix data, a new function
  for realizing the matrix into a standard R representation, and
  a new function for matrix multiplication.

[beachmat.hdf5](/packages/beachmat.hdf5)
-------------

                        Changes in version 1.4.0                        

- Migrated C++ libraries to the assorthead package. Also switch
  to the latest version fo the tatami_hdf5 library.

- Use the DelayedArray block size as the cache size limit during
  construction of tatami_hdf5::DenseMatrix, etc.

[BgeeDB](/packages/BgeeDB)
------

                        Changes in version 2.32                         

- Allows to download data from Bgee 15.2 with more bulk RNA-Seq
  and single-cell RNA-Seq

- Possibility to download h5ad single-cell data with the function
  getCellProcessedData

- Possibility to download integrated Bgee gene expression calls
  similar to Bgee gene page results with the function
  getIntegratedCalls

- topAnat has been optimized to run faster

- The getData function has been deprecated and replaced by the
  function getSampleProcessedData

[BiocGenerics](/packages/BiocGenerics)
------------

                       Changes in version 0.52.0                        

NEW FEATURES

- Define the OutOfMemoryObject class (VIRTUAL class with no slots).

- Add S4 generic containsOutOfMemoryData() and implement various
  methods.
  See '?containsOutOfMemoryData' for the details.

- Add S4 generic saveRDS() and a default method that is just a thin
  wrapper
  around base::saveRDS() that issues a warning if the object to
  serialize
  contains out-of-memory data.

[BiocNeighbors](/packages/BiocNeighbors)
-------------

                       Changes in version 1.99.0                        

- Switched to the new knncolle C++ libaries for all
  implementations. This greatly streamlines the internals and
  allows downstream packages to re-use the search indices in
  their own C++ code.

- Parallelization is now performed using the standard &#60;thread&#62;
  library.  This avoids the overhead of forking or starting new
  processes via BiocParallel.

- All *Index classes have been removed.  The output of
  buildIndex() is no longer guaranteed to be serializable, e.g.,
  external pointers to C++-owned objects; attempting to do so
  will raise a "null pointer to prebuilt index" error for all
  methods implemented in this package. The removal of this
  guarantee makes it easier to extend BiocNeighbors to new
  methods where the index structure does not have an obvious R
  representation.

- All functions (findKNN(), queryNeighbors(), etc.) will no
  longer coerce X to a matrix, to avoid the headache of S4
  dispatch ambiguity. Users should coerce their data into matrix
  format before supplying it to these functions.

- The last= option in findKNN() and queryKNN() has been replaced
  by the findDistance() and queryDistance() functions instead.
  This provides a much more intuitive method for the typical use
  of last=, i.e., to obtain the distance to the k-th nearest
  neighbor.

- findNeighbors() no longer reports each point as its own
  neighbor. Also, all neighbors are now sorted by increasing
  distance in findNeighbors() and queryNeighbors().

[BiocSingular](/packages/BiocSingular)
------------

                       Changes in version 1.22.0                        

- Bugfix for scale=TRUE with zero-variance rows.

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.62.0                        

USER VISIBLE CHANGES

- Several deprecated functions are now defunct.

BUG FIXES

- Results returned from BioMart queries will be read using Latin-1
  encounding
  if the default fails. Reported in
  https://support.bioconductor.org/p/9158844/
  (Backported to 2.60.1)

- Fixed issue when only one dataset was listed in a Mart instance,
  causing
  data.frame dimensions to be dropped.  This broke connectivity to
  https://parasite.wormbase.org. (Backported to 2.60.1)

[BioNAR](/packages/BioNAR)
------

                         Changes in version 1.7                         

- Add function to calculate the BowTie decomposition of the graph.

[Biostrings](/packages/Biostrings)
----------

                       Changes in version 2.74.0                        

NEW FEATURES

- Add ellipsis to arg list of matchProbePair() generic and methods.
  Extra arguments passed thru the ellipsis are forwarded to the
  internal
  calls to matchPattern().
  &#91;by Aidan Lakshman &#60;ahl27@pitt.edu&#62;&#93;

SIGNIFICANT USER-VISIBLE CHANGES

- Undeprecate longestConsecutive().

- A warning in translate() has been elevated to an error because
  AA_ALPHABET is now enforced.
  &#91;by Aidan Lakshman &#60;ahl27@pitt.edu&#62;&#93;

- man/XStringSetList-class.Rd previously said that "DNAStringSetList
  and
  AAStringSetList are the only constructors", but the current version
  of
  Biostrings includes BStringSetList and RNAStringSetList constructors.
  This has been updated.
  &#91;by Aidan Lakshman &#60;ahl27@pitt.edu&#62;&#93;

BUG FIXES

- replaceAmbiguities() now checks to ensure input is DNA or RNA. This
  previously caused some odd behavior with AA input.
  &#91;by Aidan Lakshman &#60;ahl27@pitt.edu&#62;&#93;

- Fix bug in consensusMatrix() when input has length zero.
  &#91;by Aidan Lakshman &#60;ahl27@pitt.edu&#62;&#93;

- Get rid of spurious warning in readQualityScaledDNAStringSet().
  &#91;by Aidan Lakshman &#60;ahl27@pitt.edu&#62;&#93;

[bluster](/packages/bluster)
-------

                       Changes in version 1.16.0                        

- Add num.threads= to all neighbor-based methods, to be passed to
  BiocNeighbors.

- Parallelize the SNN graph construction in makeSNNGraph().

[BREW3R.r](/packages/BREW3R.r)
--------

                        Changes in version 1.1.1                        

- Fix case when the input_gr_to_overlap have both stranded intervals
  and unstranded intervals (the users were facing a Warning but the
  results was not the expected one).

[BSgenome](/packages/BSgenome)
--------

                       Changes in version 1.74.0                        

- No significant changes in this version.

[BSgenomeForge](/packages/BSgenomeForge)
-------------

                        Changes in version 1.6.0                        

- No significant changes in this version.

[bsseq](/packages/bsseq)
-----

                        Changes in version 1.41                         

- Changing saveRDS to base::saveRDS following upstream changes.

[bugsigdbr](/packages/bugsigdbr)
---------

                       Changes in version 1.11.4                        

- Update stable release version for data on Zenodo

[CAGEr](/packages/CAGEr)
-----

                       Changes in version 2.12.0                        

BACKWARDS-INCOMPATIBLE CHANGES

- The range argument of the annotateCTSS function is renamed annot.
- The removeSingletons option of clustering methods is removed and
the
default value of keepSingletonsAbove is set to 0, which keeps the
standard behavior.
- In cluster objects, the dominant CTSS score is now stored in the
dominantCTSS object directly.
- The clusterCTSS function is replaced by the new paraclu and distclu
function. CTSS filtering is done beforehand with the new
filterLowExpCTSS function.
- Removed CTSSclusteringMethod() function and stop recording
clustering method name.
- Removed returnInterquantileWidth argument of functions and always
return that width when quantile information is provided. Closes #114
- The CTSStagCountDA function is removed.
- The dominant peak in TagClusters objects is now a GRanges object
like in ConsensusClusters.
- The custom method for tag clustering is removed. It was obsoleted
by
the newer CustomConsensusClusters function.
- The exportToTrack function now exports scores of tag clusters and
consensus clusters instead of setting them to zero.

NEW FEATURES

- Support the use of TxDB objects for annotating clusters.
- The plotReverseCumulatives() function now uses ggplot2. Its main,
legend, xlab, ylab, xlim and ylim arguments were removed as this can
be controlled via ggplot2 functions.
- New TSSlogo function wrapping the ggseqlogo package.
- New distclu and paraclu functions that can run directly on CTSS
objects. You can use them to test parameters before running the
whole CAGEexp object through clusterCTSS.

BUG FIXES

- CTSS filtering now works correctly with threshold = 0,
thresholdIsTpm = TRUE.
- The importPublicData function was repaired for FANTOM samples.
- Remove broken and obsolete customClusters method of clusterCTSS.
Use
CustomConsensusClusters() instead. Fixes #113.
- Ensure consensus clusters can be exported as tracks with
interquantile width information. Fixes #108, #70.

OTHER CHANGES

- Accelerated the computation of cumulative sums ~10×.
- Singleton filtering is now done by the paraclu and distclu
functions
themeselves; .ctss_summary_for_clusters does not change the input
clusters except for adding information.

[Cardinal](/packages/Cardinal)
--------

                 Changes in version 3.7.8 (2024-10-26)                  

SIGNIFICANT USER-VISIBLE CHANGES

- Deprecate 'slice()' in factor of 'sliceImage()'
  to avoid conflicts with users attaching dplyr

                 Changes in version 3.7.7 (2024-10-24)                  

SIGNIFICANT USER-VISIBLE CHANGES

- Update 'peakAlign()' with 'binratio' parameter

                 Changes in version 3.7.6 (2024-10-16)                  

NEW FEATURES

- Reading/writing imzML now preserves 'metadata' slot

- Reading/writing imzML now supports parallel processing

SIGNIFICANT USER-VISIBLE CHANGES

- Update 'estimateDomain()' to select resolution
  based on 'median', 'min', 'max', or 'mean'

- Update 'peakAlign()' to build domain bins based on
  the minimum of peak gaps (instead of the median)

- Write log file when writing imzML/Analyze files

BUG FIXES

- Fix bug in 'process()' resulting in output matrix
  with 'list' elements under certain conditions

                        Changes in version 3.7.5                        

NEW FEATURES

- Add 'setCardinalParallel()' for setting a parallelization
  backend with reasonably-selected defaults

SIGNIFICANT USER-VISIBLE CHANGES

- Add 'SAR=FALSE' option to 'simulateImage()' for faster
  simulation if spatial autoregressive model is not needed

- Change default array order for 'MSImagingArrays' so that
  'intensity' is first array and 'mz' is second array

- Improved logging with pre-processing functions

BUG FIXES

- Fix leaky 'meansTest()' closures

                        Changes in version 3.7.4                        

NEW FEATURES

- Add 'saveCardinalLog()' for saving log file

- Add 'getCardinalLogger()' + 'setCardinalLogger()'

- Add 'fetch()' and 'flash()' methods for moving
  spectra between shared memory and temporary files

SIGNIFICANT USER-VISIBLE CHANGES

- Update compatibility with matter 2.7.6

- Changes to 'simulateSpectra()' and 'simulateImage()'
  in parallelization and spectral noise generation

- RNG in 'simulateSpectra()' and 'simulateImage()'
  now warn if 'RNGkind()' is not "L'Ecuyer-CMRG"

BUG FIXES

- Fix calculation of "adaptive" spatial weights
  in 'spatialWeights()' for accuracy and stability

                        Changes in version 3.7.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add 'setCardinalChunksize()' + 'getCardinalChunksize()'

- Update vignettes due to spectral processing updates

BUG FIXES

- Fix plot() and image() not always respecting dimnames

- Fix plot() and image() not plotting multiple columns

- Fix plot() and image() failing for non-syntactic names

                        Changes in version 3.7.2                        

BUG FIXES

- Merge bug fix from 3.6.4.

                        Changes in version 3.7.1                        

NEW FEATURES

- Update 'spectrapply()' so index array is optional

- Update 'plot()' to support 2 domains (e.g., for ion mobility)

- Add plotting methods for 'XDataFrame' and 'PositionDataFrame'

BUG FIXES

- Fix default 'tolerance' in 'pixels()'

                        Changes in version 3.6.6                        

BUG FIXES

- Fixes for 'simulateImage()' and 'simulateSpectra()'

                        Changes in version 3.6.5                        

BUG FIXES

- Fix 'spatialShrunkenCentroids()' failing when 'r=0'
  or when a pixel has no neighboring pixels

                        Changes in version 3.6.4                        

BUG FIXES

- Bug fixes for 'plot()' on spectra with queued processing

                        Changes in version 3.6.3                        

BUG FIXES

- Bug fixes for 'meansTest()' when random effects are specified

                        Changes in version 3.6.2                        

BUG FIXES

- Version bump for 'matter' 2.6.2 bugfixes

                        Changes in version 3.6.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add 'mass.range' and tolerance' arguments to 'bin()'

BUG FIXES

- Fix 'readImzML' error if imzML fails conversion to 'ImzMeta'

[CardinalIO](/packages/CardinalIO)
----------

                 Changes in version 1.3.6 (2024-10-24)                  

BUG FIXES

- Fix writing scientific notation to imzML output
  (note: only affected parsing, not precision)

                 Changes in version 1.3.5 (2024-10-16)                  

NEW FEATURES

- Support for parallel writing with BiocParallel

- Added 'BPPARAM' argument to 'writeAnalyze()'

- Added 'BPPARAM' argument to 'writeImzML()'

SIGNIFICANT USER-VISIBLE CHANGES

- Return 'outdata' attribute on 'writeAnalyze()'
  results with output 'matter' object

- Return 'outdata' attribute on 'writeImzML()'
  results with output 'matter' object

                        Changes in version 1.3.4                        

SIGNIFICANT USER-VISIBLE CHANGES

- Parse imzML "binary data compression type" tags

- Compressed ibd data arrays are now attached as
  raw byte arrays (_without_ decompression)

                        Changes in version 1.3.3                        

SIGNIFICANT USER-VISIBLE CHANGES

- Changed default 'mz.type' from "float32" to "float64"

                        Changes in version 1.3.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Renaming imzML files when 'asis=TRUE' now signals a warning

                        Changes in version 1.3.1                        

BUG FIXES

- Update broken unit test to reflect recent changes in 'matter'

                        Changes in version 1.2.1                        

SIGNIFICANT USER-VISIBLE CHANGES

- Allow multiple options to 'check' in 'parseImzML'

BUG FIXES

- Check ibd file size against imzML binary data array offsets


[cbaf](/packages/cbaf)
----

                 Changes in version 1.28.0 (2024-10-24)                 

New Features

- Introducing a new two step procedure for processMultipleStudies! If
  you
  intend to download the required data on server, CBAF can download and
  create a database, then it compress the database into a ZIP file.
  User can
  then transfer the zipped file into their local computer and CBAF can
  use
  that ZIP file to recreate the intial database that it needs. Since
  servers
  usually lack any graphical units, CBAF won't be able to generate
  heatmaps
  when running on server.

[cellxgenedp](/packages/cellxgenedp)
-----------

                        Changes in version 1.10                         

SIGNIFICANT USER-VISIBLE CHANGES

- (v.1.9.1) Add progress bar when updated db()

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.39.3                        

- Fix the documentation for Vennerable.

                       Changes in version 3.39.2                        

- Fix the documentation for summarizePatterInPeaks.

                       Changes in version 3.39.1                        

- Unique promoters for genomicElementDistribution.

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.41.3                        

- Better covplot(). Support universal chromosome names, and keep the
default order of multiple peaks when plot a list of GRanges object.
- Robust generate_colors(). Edit the logical of decision, and can
validate color code automatically.
- Extend dplyr verbs (filter(), mutate(), arrange(), rename()) to
peak
(GRanges object or data.frame), see #242.

                       Changes in version 1.41.2                        

- Enhancement of plotDistToTSS(), see #241.

                       Changes in version 1.41.1                        

- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)

[ClassifyR](/packages/ClassifyR)
---------

                       Changes in version 3.10.0                        

- 
  Metric calculation done at fold level rather than permutation
  level.

[clusterProfiler](/packages/clusterProfiler)
---------------

                       Changes in version 4.13.4                        

- re-export DOSE::enrichDO() and DOSE::gseDO() (2024-10-01, Tue)

                       Changes in version 4.13.3                        

- fixed bug in enrichPC() (2024-08-26, Mon)

                       Changes in version 4.13.2                        

- fixed bug of gson_KEGG() (2024-08-19, Mon)

                       Changes in version 4.13.1                        

- update functions to access PathwayCommons data (2024-08-11, Sun,
gson#9)
- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)
- update kegg_category information (7 categories and 572
subcategories) (2024-07-26, Fri)
- Cellular Processes (36)
- Drug Development (75)
- Environmental Information Processing (41)
- Genetic Information Processing (39)
- Human Diseases (99)
- Metabolism (190)
- Organismal Systems (92)

[ClustIRR](/packages/ClustIRR)
--------

                       Changes in version 1.3.34                        

- Stable version but not tested

                       Changes in version 1.3.25                        

- Major changes since this version including
  * community detection
  * Stan model-based differential community occupancy
  * supporting functions

                       Changes in version 1.3.20                        

- Future based parallelization in clustering and graphing

                        Changes in version 1.3.1                        

- Stub functions for integration of data from VDJdb, TCR3d and
  MCPAS-TCR

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

                       Changes in version 2.21.1                        

- `pheatmap()`: `na_col` is passed to `Heatmap()`.

[CompoundDb](/packages/CompoundDb)
----------

                         Changes in version 1.9                         

Changes in version 1.9.5

- Add new extractByIndex() method.

Changes in version 1.9.4

- compound_tbl_lipidblast supports now parallel processing and
extracts more information from MoNA's JSON format (thanks to Prateek
Arora for contribution).

Changes in version 1.9.3

- compound_tbl_lipidblast: ensure exactmass is of type numeric.

Changes in version 1.9.2

- compound_tbl_lipidblast: add parameter n to support reading and
processing MoNA json files in sets (chunks) of lines at a time and
hence reduce memory demand for very large files.

Changes in version 1.9.1

- Allow CompDb to store that database name as alternative to an
active
database connection. This allows to serialize and load an object
to/from disk (serializing an active database connection would not be
possible) . Each call to extract data from the database will however
open (and close) its own connection.

[COTAN](/packages/COTAN)
-----

                       Changes in version 2.5.11                        

Fixed bug in the function cellsUMAPPlot(): restored possibility of
passing a genes vector as genesSel parameter. Also updated the
documentation about the available genes selection methods

                       Changes in version 2.5.10                        

Fixed typo in error message

Fixed bug in function genesCoexSpace(): now primaryMarkers can have
only
a single gene

                        Changes in version 2.5.9                        

Exported utility functions about names arrays: conditionsFromNames(),
niceFactorLevels(), factorToVector()

                        Changes in version 2.5.8                        

feature/add_condition_arguments_to_plot_functions Now the following
plot
functions take in conditions explicitly, instead of just instructions
to
determine them from the cells' names. The changes involved:
cellSizePlot(), ECDPlot(), genesSizePlot(),
mitochondrialPercentagePlot(), scatterPlot()

Added possibility to convert COTAN objects to/from
SingleCellExperiment
objects. SCE objects created by the Seurat package are supported

Hardened arguments' checks for function UMAPPlot()

Solved issue with the function establishGenesClusters(): it was
throwing
an error when one of the sub-lists in the groupMarkers argument did
contain only one element

                        Changes in version 2.5.7                        

Introduced new way to check for the Uniform-Transcript property of
the
clusters based on multiple thresholds calibrated so that the new
method
is more effective at describing really statistically uniform clusters

Functions cellsUniformClustering() and mergeUniformCellsClusters()
have
been re-factored so to support new class hierarchy for UT checkers.
This
allows user to select which method to use for the checks; as of now
the
following methods are supported:

- "SimpleGDIUniformityCheck"
- "AdvancedGDIUniformityCheck"

Avoided issue with pdf file creation: file handle was not closed in
case
of errors

Added possibility of choosing number of features in seuratHVG()

Solved minor issue with with clusterization functions in cases when
only
one cluster was created

                        Changes in version 2.5.6                        

Made function heatmapPlot() more easy to use and in line with the
rest
of the COTAN package

Now the method storeGDI() can take in the output data.frame from the
function calculateGDI()

Solved few minor issues with the vignette and changed a few default
parameters in cellsUMAPPlot(), pValueFromDEA() and
findClustersMarkers()

                        Changes in version 2.5.5                        

Stopped function cellsUniformClustering() from saving the internally
created Seurat object due to possibly long saving times

Split the now deprecated function getNormalizedData() into two
separated
functions: getNuNormData() and getLogNormData()

Re-factored function mergeUniformCellsClusters() to be more precise:
now
it merges clusters starting from the most similar in latest batch and
also runs the merging in multiple steps adjusting gradually the GDI
threshold ranging from a very strict up to the user given ones.

Fixed minor bugs in functions GDIPlot() and
clustersMarkersHeatmapPlot()

                        Changes in version 2.5.4                        

Added possibility to display UMAP plots of cells clusters, using the
function cellsUMAPPlot()

                        Changes in version 2.5.3                        

Updated the vignette to the most recent changes

Allowed user to set the ratio of genes above the threshold allowed in
a
Uniform Transcript cluster

                        Changes in version 2.5.2                        

Solved issue with usage checks about the torch library

Allowed user to explicitly opt-out from the torch library usage:
COTAN
will avoid torch commands when the option "COTAN.UseTorch" is set to
FALSE

                        Changes in version 2.5.1                        

Added support for the torch library to help with the heavy lifting
calculations of the genes' COEX matrix, with consequent substantial
speed-up, especially when a GPU is available on the system

                        Changes in version 2.5.0                        

First release in Bioconductor 3.20

[cqn](/packages/cqn)
---

                        Changes in version 1.51                         

- Various small changes addressing changes in dependencies. For
  example, decideTestsDGE got renamed to decideTests. Likewise, a
  number of imports needed fixing.

[crisprBase](/packages/crisprBase)
----------

                        Changes in version 1.9.1                        

- Added Csm complex RNA-targeting nuclease.

[crisprDesign](/packages/crisprDesign)
------------

                        Changes in version 1.7.1                        

- Fixed getTssObjectFromTxObject for negative strand.

[crisprScore](/packages/crisprScore)
-----------

                        Changes in version 1.9.3                        

- Fixed documentation for enPamGB and DeepCpf1 algorithms.

                        Changes in version 1.9.2                        

- Fixed python environments based on the lastest changes in
  basilisk (using conda-forge installations).

[CTdata](/packages/CTdata)
------

                         Changes in version 1.5                         

CTdata 1.5.3

- Correction of a bug in mean_methylation_in_embryo and
mean_methylation_in_FGC datasets
- Vignette update

CTdata 1.5.2

- Correction of an error in the upload

CTdata 1.5.1

- Used more stringent threshold to define genes activated in CCLE
cell
lines ang in TCGA tumors (TPM > 1 in at least 1% of tumors and
cancer cell lines, and TPM > 5 in at least one tumor and cancer cell
line)
- Changed selection criteria of testis-specific and
testis-preferential genes in GTEX
- Added a slightly more stringent criteria of selection of
testis-specific genes when analysing RNAseq data of normal tissues
with multimapping
- Minor modifications in selection of genes induced by DAC
- Set TCGA_catgeory in TCGA_TPM to "leaky" when the q75 expression in
normal peritumoral TCGA samples is higher than 0.5
- Added HPA_cell_type_specificity (a table giving cell type
specificity of each genes based on Human protein Atlas scRNAseq
data)
- Added FGC_sce, scRNAseq data of human fetal gonads
- Added oocytes_sce, scRNAseq data of human oocytes
- Changed preliminary CT_list into all_genes_prelim list to add the
characterisation to all genes and not only CT in the package
- Changed criteria for testis specificity, including HPA/TCGA/CCLE
category in the definition
- Added CT_gene_type column to differentiate CT genes and CT
preferential genes
- Included all genes in methylation objects (methylation_in_tissues,
mean_methylation_in_tissues and TCGA_methylation)
- Added all_genes containing all the characterisation and analysis
for
all genes
- Changed the criteria for regulation by methylation, needing DAC
induction and methylation in somatic only
- Reordered CT_genes and all_genes
- Changed all documentation
- Added embryo, FGC and hESC data

CTdata 1.5.0

- New devel

[CTexploreR](/packages/CTexploreR)
----------

                         Changes in version 1.1                         

CTexploreR 1.1.3

- Add pkgdown config.
- Correcting typo in fetal_germcells_mean_methylation() function name
- Adding tests for new functions

CTexploreR 1.1.2

- Created functions to visualise expression in fetal germ cells,
oocytes, hESC and early embryos.
- Created functions to visualise methylation in fetal germ cells,
hESC
and embryo.

CTexploreR 1.1.1

- Adaptation due to updates in CTdata : adding all genes to
functions,
choice to include CT preferential.

CTexploreR 1.1.0

- New devel

[cydar](/packages/cydar)
-----

                       Changes in version 1.30.0                        

- Updated to use the latest BiocNeighbors, which means that the
  output of prepareCellData() is no longer serializable. This is
  because the precomputed index is now an external pointer to a
  C++-owned data structure, and cannot be moved between sessions.
  (Arguably, it was a mistake to expose these internals in the
  expected workflow for this package, but that ship has sailed.)
  Hopefully, this should not affect most users as they should be
  primarily interacting with the CyData object.

[cypress](/packages/cypress)
-------

                 Changes in version 1.1.1 (2024-05-27)                  

- Reduce dependent packages

- Revise the rd file

[CytoMDS](/packages/CytoMDS)
-------

                         Changes in version 1.1                         

CytoMDS 1.1.4

- implemented pkgdown customization

CytoMDS 1.1.2 & 1.1.3

- expression matrices as input to EMD calculation

CytoMDS 1.1.1

- added citation to bioRxiv pre-print

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

                         Changes in version 1.5                         

CytoPipeline 1.5.2

- updated unit test to account for flowAI version change
- now suggesting CytoPipelineGUI package

CytoPipeline 1.5.1

- updated processing step argument matching using phenoData

[CytoPipelineGUI](/packages/CytoPipelineGUI)
---------------

                         Changes in version 1.3                         

CytoPipelineGUI 1.3.1

- implemented pkgdown site customization

[dar](/packages/dar)
---

                        Changes in version 1.1.2                        

Bug Fixes

- data_import.Rmd

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

                        Changes in version 1.1.1                        

- Version up.

                     Changes in version 1.1.0.9000                      

NEW FEATURES

- runDegCre now accepts "qvalue" as an input to pAdjMethod which is
now the default. Our testing shows that this method improves the
performance of DegCre predictions. The qvalue calculation does not
require the user to specify a alphaVal.

- Added new function collapseDegCreToGene which converts DegCre
results for gene with multiple TSSs to use only the association that
spans the shortest distance.

- Added new function calcAssocProbOR that calculates the odds-ratio
for DegCre association probabilities. This function can operate on
assocProb or rawAssocProb values by altering the type flag. It is
meant to replace calcRawAssocProbOR.

- Added new function convDegCreResListToCreGeneScoreGR that converts
DegCre results to a simplified GRanges with the predicted gene and
score as metadata.

OTHER CHANGES

- Streamlined runDegCre code to use more subfunctions.

BUG FIXES

- Changed the internal function makePlotGInter to return a sorted
GInteractions to avoid inconsistent returns.

                        Changes in version 1.0.2                        

- Changed citation to Genome Research article

BUG FIXES

- Added parameter minNDegs to optimizeAlphaDegCre to stop it from
trying testedAlphaVals that result in the number of passing DEGs to
be too low for the optimization algorithm to function properly.

[DEGreport](/packages/DEGreport)
---------

                       Changes in version 1.41.1                        

- made that the sequence of clusters in the dendogram matches the
  sequence of clusters in the trajectories @zellerivo

- Fix error in check report about stack limit size

- Add option to avoid running dendextend
  https://github.com/lpantano/DEGreport/issues/62

[DelayedArray](/packages/DelayedArray)
------------

                       Changes in version 0.32.0                        

NEW FEATURES

- Add nzwhich() method for DelayedArray objects (block-processed).

- Support %*%, crossprod(), and tcrossprod() between DelayedMatrix and
  COO_SparseMatrix objects.

SIGNIFICANT USER-VISIBLE CHANGES

- The default realize() method now returns an SVT_SparseArray object
  instead of an SparseArraySeed object when the array-like object to
  realize is sparse and 'BACKEND' is NULL.

- Coercing a DelayedArray object or derivative to SparseArray should be
  much more efficient (thanks to various tweaks that happened in the
  SparseArray and HDF5Array packages).

DEPRECATED AND DEFUNCT

- Deprecate SparseArraySeed objects.

- Deprecate OLD_extract_sparse_array() and read_sparse_block() generics
  and methods.

- Fix bug in coercion from DelayedArray to SparseArray when the object
  to
  coerce has NAs. See
  https://github.com/Bioconductor/HDF5Array/issues/61

[DelayedMatrixStats](/packages/DelayedMatrixStats)
------------------

                        Changes in version 1.27                         

- Remove colAnyMissings() and rowAnyMissings(), which were
  deprecated in Version 1.23 and made defunct in Version 1.25.

[demuxSNP](/packages/demuxSNP)
--------

                        Changes in version 1.4.0                        

Added reassign_centroid() function for improved handling of sparse
and
missing data. Added messaging to previous functions reassign(),
reassign_jaccard() and reassign_balanced() to point towards
reassign_centroid(). Updated add_snps() function to add unfiltered
data.

[DepecheR](/packages/DepecheR)
--------

                 Changes in version 1.21.4 (2024-08-22)                 

- Minor bug fixes in neighSmooth.

                 Changes in version 1.21.3 (2024-06-21)                 

- Minor bugs in depeche affecting the special case of one defining
  marker and two clusters.

                 Changes in version 1.21.2 (2024-06-14)                 

- Correction of the number of used neighbors in neighSmooth function.

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.45.1                        

- Dropping the requirement for DESeqDataSet to contain integers.
  This allows fitType="glmGamPoi" to be used. If non-integer
  counts are present and other fitType values are provided,
  the various functions will give an error, e.g. DESeq,
  estimateDispersions, nbinomWald, and nbinomLRT. To avoid
  integer conversion, use DESeqDataSet() with
  skipIntegerMode=TRUE. For more detail see:
  https://github.com/thelovelab/DESeq2/issues/66

[DNABarcodeCompatibility](/packages/DNABarcodeCompatibility)
-----------------------

                 Changes in version 1.21.1 (2024-09-18)                 

- Addition of the distance() function from the DNABarcodes package

[DOSE](/packages/DOSE)
----

                       Changes in version 3.99.1                        

- return NULL in GSEA if not genes can be mapped (2024-08-26, Mon,
ReactomePA#43)

                       Changes in version 3.31.4                        

- remove mpoSim() and hopSim() (2024-08-21, Wed)

                       Changes in version 3.31.3                        

- remove gseMPO() and gseHPO(), instead using e.g., gseDO(ont="HPO")
(2024-08-15, Thu)
- remove enrichMPO() and enrichHPO(), instead using e.g.,
enrichDO(ont="HPO")
- unify API and removing the usages of MPO.db and HPO.db
- update DO-gene mapping data (2024-08-13, Tue)
- remove HDO.db and use new API in GOSemSim (v>=2.31.1) (2024-08-13,
Tue)
- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)

                       Changes in version 3.31.2                        

- add FoldEnrichment, RichFactor and zScore in ORA result
(2024-06-13,
Thu)

                       Changes in version 3.31.1                        

- fixed bug in options(enrichment_force_universe=TRUE) (2024-05-16,
Thu)
-
https://github.com/YuLab-SMU/clusterProfiler/issues/283#issuecomment-2077869230

[dreamlet](/packages/dreamlet)
--------

                        Changes in version 1.3.3                        

- July 19, 2024
- fix in .read_matrix_block()
- See
https://github.com/GabrielHoffman/dreamlet/pull/23/commits/c453ac98ebc0329279b4dd3ae26a674df0e9b1f2

                        Changes in version 1.3.1                        

- May 31, 2024
- bump Bioc version

[easyRNASeq](/packages/easyRNASeq)
----------

                       Changes in version 2.41.1                        

- Ported changes from 2.40.1

                       Changes in version 2.40.1                        

- Fixed broken dependencies (Biostrings::type is now BiocGenerics::type
  and BiocGenerics::clusterApply is now parallel::clusterApply).

- Adapted to changes for Roxygen

- Changed from using biomaRt::useMart to biomaRt::useEnsembl

[edgeR](/packages/edgeR)
-----

                 Changes in version 4.4.0 (2024-10-30)                  

- 
  New function catchRSEM() to read transcript-level
  quantifications from RSEM output. If the RSEM output includes
  resampling replicates, then catchRSEM() uses them to estimate
  the over-dispersion arising from read-to-transcript-ambuity for
  each transcript. Similar to catchSalmon() and catchKallisto()
  but for RSEM output.

- 
  New function normalizeBetweenArrays.DGEList() to apply
  microarray-style normalization to a DGEList by setting the
  offset matrix appropriately.

- 
  New argument `keep.unit.mat` for glmQLFit(). The unit matrices
  produced when 'legacy=FALSE` are not required for routine
  downstream analysis, so they are now not returned by glmQFit()
  unless `keep.unit.mat=TRUE`. This reduces the size of the
  glmQLFit fitted model object.

- 
  Output components var.prior and var.post from glmQLFit() have
  been renamed to s2.prior and s2.post.

- 
  glmQLFit() with `legacy=FALSE` now uses the improved empirical
  Bayes hyperparameter estimation in limma 3.61.9, which is
  designed especially for scenarios when the residual degrees of
  freedom are unequal between genes. Now that unequal residual df
  are better accounted for, very small df.residual values are no
  longer floored to zero before performing empirical Bayes
  estimation.

- 
  The default value of `top.proportion` in glmQLFit() with
  `legacy=FALSE` now depends on the number of genes and on the
  residual degrees of freedom, with more genes and more df giving
  smaller values. Also, the `DGEList` method for glmQLFit() with
  `legacy=FALSE` will take the dispersion from the mean of the
  right tail values of the trended dispersions, instead of
  re-estimating, if these are found in the `DGEList` object.

- 
  makeCompressedMatrix() now allows `x` to optionally be a row
  vector (matrix with one row) or a column vector (matrix with
  one column) or an ordinary vector. Previously, matrix values
  for `x` were simply returned as output.

- 
  diffSpliceDGE() now passes weights to glmFit.default(), if they
  exist.

- 
  A major revision has been undertaken to the C source code. C++
  has been removed in favor of pure C and, instead, the C
  functions are now wrapped into R in modern style. edgeR no
  longer depends on the Rcpp package. The new C code includes
  careful memory management and a simplified file structure.

- 
  Various edits to help pages. The edgeR publications in
  edgeR-package.Rd are now listed in reverse chronological order.
  Some help page cross-references have been fixed. The term
  "scalar" to indicate a numeric vector of length one has been
  replaced with "single value" in several places. The DGEGLM and
  glmFit help pages now clarify that the output `coefficients`
  component is on the natural log scale.

- 
  Add checks for negative or NA counts to cpm(), cpmByGroup() and
  normLibSizes().

- 
  Minor R code improvements that do not change the user interface
  to various functions, for example replacing any(is.na()) with
  anyNA() and using identical() when testing conditions within
  if() statements.

- 
  Remove decidetestsDGE(), whose functionality is now provided by
  the generic function decideTests().

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.25.6                        

- pretty gene count legend (2024-10-29, Tue, #271)

                       Changes in version 1.25.5                        

- new emaplot(), goplot(), cnetplot() and ssplot(), all power by
'ggtangle' package (2024-10-24, Thu)
- re-export ggtangle::cnetplot() (2024-10-24, Thu)
- remove drag_network() (2024-10-24, Thu)

                       Changes in version 1.25.4                        

- fixed goplot() (2024-10-23, Wed, #297, #732, #718)

                       Changes in version 1.25.3                        

- hplot(): Horizontal plot for GSEA result (2024-08-27, Tue)

                       Changes in version 1.25.2                        

- fixed bug in ridgeplot() (2024-08-19, Mon, clusterProfiler#704)

                       Changes in version 1.25.1                        

- fixed GeneRatio in dotplot as character of fraction issue
(2024-08-16, Fri, clusterProfiler#715)
- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)
- dotplot2 to compare two selected clusters in 'compareClusterResult'
object (2024-06-15, Sat)
- volplot to visualize ORA result using volcano plot (2024-06-13,
Thu)

[enrichViewNet](/packages/enrichViewNet)
-------------

                        Changes in version 1.3.4                        

NEW FEATURES

- The 'force' parameter in createEnrichMapMultiComplex() has been
  removed.

- The 'force' parameter in createEnrichMapMultiBasic() has been
  removed.

- The 'force' parameter in createEnrichMap() has been removed.

                        Changes in version 1.3.3                        

NEW FEATURES

- The createEnrichMapMultiComplex() now has an new parameter ... that
  enables passing extra parameters to imbedded emapplot() function.

- The createEnrichMapMultiBasic() now has an new parameter ... that
  enables passing extra parameters to imbedded emapplot() function.

- The createEnrichMap() now has an new parameter ... that enables
  passing extra parameters to imbedded emapplot() function.

                        Changes in version 1.3.2                        

BUG FIXES

- The package uses the Ensembl list rather than the query list to
  extract information when intersection column is missing, in the
  extractInformationWhenNoIntersection() function.

                        Changes in version 1.3.1                        

BUG FIXES

- The package uses the information about the organism when querying the
  gprofiler database.

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.29.1                        

- Improve documentation of `proteinToGenome()` and require named
  `IRanges` as
  input to `proteinToGenome()`.

[epialleleR](/packages/epialleleR)
----------

                 Changes in version 1.13.4 (2024-10-10)                 

- VariantAnnotation moved to Suggests

                 Changes in version 1.13.3 (2024-10-04)                 

- plotPatterns for pretty plotting

[EpiCompare](/packages/EpiCompare)
----------

                        Changes in version 1.9.5                        

New features

- Remove the soon-to-be-deprecated BRGenomics dependency.
- Port tidyChromosomes function to EpiCompare.

Miscellaneous

- Update maintainer details.

[epiregulon.extra](/packages/epiregulon.extra)
----------------

                        Changes in version 1.0.1                        

fixed documentation for plotHeatmapActivity fixed legend label of
summary.logFC in plotBubble removed duplicate genes from
plotHeatmapRegulon

[erccdashboard](/packages/erccdashboard)
-------------

                       Changes in version 1.39.5                        

- Depends on R > 4.0

- Converted from Sweave PDF vignette to Rmarkdown HTML.

- Updating testDECount.R to handle changes made in edgeR

- Updating saveERCCPlots.R to deal with grob error

[escape](/packages/escape)
------

                 Changes in version 2.1.5 (2024-10-23)                  

- update handling of v5 Seurat versus <v5 Seurat Objects

                 Changes in version 2.1.4 (2024-09-13)                  

- update densityEnrichment() GSVA function pull

                 Changes in version 2.1.3 (2024-09-13)                  

#VERSION BUMP FOR BIOCONDUCTOR

UNDERLYING CHANGES

- update densityEnrichment() for new GSVA function name
- Parallelization of performNormalization()
- Refactor of getGeneSets() to prevent issues with m_df error.

                 Changes in version 2.0.1 (2024-07-26)                  

UNDERLYING CHANGES

- fixed performNormalziation() errors when input.data was a matrix,
now requires single-cell object and enrichment data
- passing parallel processing properly to runEscape() function.

[EWCE](/packages/EWCE)
----

                       Changes in version 1.13.1                        

Bug fixes

- Making ewce_plot functionality more clear - will no fail if
make_dendro=TRUE and ggdendro is not installed or CTD is not
provided rather than issuing a warning.

[ExperimentHub](/packages/ExperimentHub)
-------------

                       Changes in version 2.13.0                        

BUG FIXES

- (2.13.1) merged @votti PR for not checking connection if localhub is
  selected

[extraChIPs](/packages/extraChIPs)
----------

                        Changes in version 1.9.6                        

- Added the function centrePeaks() to recentre peaks using any files
with coverage

[faers](/packages/faers)
-----

                        Changes in version 1.1.6                        

- fda_drugs() now directly use a fixed url to download the data

- faers_meta(internal = TRUE) will always use the cache data in the
package.

- Rename "gndr_cod" into "gender" for periods before 2014q2, and
rename "sex" into "gender" for periods after or equal to 2014q2.

- "sex" was added, which recoded any values other than "F" or "M" as
NA.

                        Changes in version 1.1.4                        

- fix error when download failed

- meddra augment additional argument primary_soc to help save the
full
meddra data

[fenr](/packages/fenr)
----

                        Changes in version 1.2.1                        

- Go term namespace added to the information extracted by fetch_go.

[fgsea](/packages/fgsea)
-----

                       Changes in version 1.31.6                        

- plotCoregulationSpatial supports multiple samples

                       Changes in version 1.31.3                        

- check for empty strings in gene names

- minor fix in leading edge calculation

                       Changes in version 1.31.2                        

- fora() provides fold enrichment scores for more effective
  prioritisation of results

[flowAI](/packages/flowAI)
------

                       Changes in version 1.35.2                        

- Substantial modification of the anomaly detection in the flow rate.
  Added a parameter to remove parts of the flow rate from the mode.
  Added the Loess regression as alternative method. Removed the
  deviationFR parameter.

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

                       Changes in version 1.13.1                        

- Fix bugs in vignettes

[gDRcore](/packages/gDRcore)
-------

               Changes in version 2024-10-24 (2024-10-24)               

- split and refactor annotation functions

               Changes in version 2024-10-21 (2024-10-21)               

- add support for reprocessing data with old data model

               Changes in version 2024-09-30 (2024-09-30)               

- handle properly additional perturbations to get rid of creating
additional columns that can't be hangle properly by the app

               Changes in version 2024-08-21 (2024-08-21)               

- identify additional perturbations hidden in the secondary drug

               Changes in version 2024-08-19 (2024-08-19)               

- utilize calc_sd function

               Changes in version 2024-08-09 (2024-08-09)               

- fix issue with wrong mapping of Day0 data

               Changes in version 2024-08-08 (2024-08-08)               

- fix issue with mapping overrides untreated controls

               Changes in version 2024-08-05 (2024-08-05)               

- fix issue with non-avearaged concentration data

               Changes in version 2024-07-23 (2024-07-23)               

- fix issue with providing empty nested_confounder

               Changes in version 2024-07-17 (2024-07-17)               

- allow using custom functions for calculating HSA and Bliss scores
for combination data

               Changes in version 2024-07-15 (2024-07-15)               

- refactor logic for calculating standard deviation for single values

               Changes in version 2024-07-10 (2024-07-10)               

- update unit tests

               Changes in version 2024-06-04 (2024-06-04)               

- switch to get_supported_experiments

[gDRimport](/packages/gDRimport)
---------

               Changes in version 2024-10-22 (2024-10-22)               

- fix test for load_files

               Changes in version 2024-08-29 (2024-08-29)               

- fix failing github page

               Changes in version 2024-08-12 (2024-08-12)               

- update function for importing public PRISM data

               Changes in version 2024-07-03 (2024-07-03)               

- fix a bug with incorrect recognition of input file types

[gDRstyle](/packages/gDRstyle)
--------

               Changes in version 2024-08-29 (2024-08-29)               

- add GitLab credentials to pkgdown builds

               Changes in version 2024-06-10 (2024-06-10)               

- add debug to the exceptions list

[gDRutils](/packages/gDRutils)
--------

               Changes in version 2024-10-11 (2024-10-11)               

- make duplicates' helpers supporting combo assays as well

               Changes in version 2024-10-07 (2024-10-07)               

- refactor the logic for dealing with duplicates in assay data

               Changes in version 2024-10-03 (2024-10-03)               

- fixed issue in average_biological_replicated (fit_type)

               Changes in version 2024-09-16 (2024-09-16)               

- add functions set_unique_cl_names_dt and set_unique_drug_names_dt

               Changes in version 2024-09-04 (2024-09-04)               

- remove hack with checkDimnames

               Changes in version 2024-08-30 (2024-08-30)               

- remove GDS fit_source from gDRviz

               Changes in version 2024-08-28 (2024-08-28)               

- extend the logic of get_additional_variables to support other
sources of fitting metrics

               Changes in version 2024-08-14 (2024-08-14)               

- extend the logic of average_biological_replicates_dt to calculate
standard deviation

               Changes in version 2024-08-06 (2024-08-06)               

- add functions for setting unique identifiers in the colData and
rowData of SE

               Changes in version 2024-07-30 (2024-07-30)               

- refactor average_biological_replicates_dt and
get_additional_variables to support unprettified identifiers

               Changes in version 2024-07-17 (2024-07-17)               

- update define_matrix_grid_positions

               Changes in version 2024-07-12 (2024-07-12)               

- move get_combo_col_settings and get_iso_colors to gDRplots package

               Changes in version 2024-07-08 (2024-07-08)               

- add residual sum of square and p-value to Metrics assay

               Changes in version 2024-07-03 (2024-07-03)               

- add vignette section about prettifying logic

               Changes in version 2024-06-24 (2024-06-24)               

- fixed issue in vignette

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.40.2                        

UTILITIES

- add "#define STRICT_R_HEADERS 1" to the C++ header file according to
  R_r86984

- update the C codes according to '_R_USE_STRICT_R_HEADERS_=true' &
  '_R_CXX_USE_NO_REMAP_=true'

[GeDi](/packages/GeDi)
----

                        Changes in version 1.2.0                        

- Smaller bug and typo fixes

- Fixed the bug that the getGenes() function would set all gene names
to all caps which lead to the inability to download the correct PPI
information for species like mouse. Also renamed the function to
prepareGenesetData() to reflect more accurately its behavior.

- Updated the default version used in the getId() and getStringDB()
to 12.0, the current version of the String database.

- Fixed the broken zoom feature in the Optional Filtering Step in the
Data Input panel. Additionally added a column Description to the
table of zoomed gene sets to facilitate interpretation.

- Fixed that the clustering will now be reset whenever a new score is
calculated.

- Updated the checkInclusion() function to drastically reduce
runtime.

- Fixed the error that the value of alpha would not be properly pass
to all the sub function used by the getpMMMatrix() function.

- Replaced all occurrences of PMM with pMM to match the notation of
the original publication.

- Replaced the kNN clustering algorithm with PAM (partitioning around
mendoids) as this seems to be more suitable for enrichment data
represented by distance scores.

- Renamed the goSimilarity() function to goDistance() to better
indicated that this is a distance rather than a similarity score.
Also scaled all scores to the &#91;0, 1&#93; interval.

- Fixed the normalization function in the getKappaDistanceMatrix()
function.

- Changed the implementation of the Louvain clustering algorithm to
use a weighted graph now. The graph is weighted by the distance
scores between gene sets.

- Added GeneTonicList as a possible input object for GeDi. Now GeDi
is
directly compatible with GeneTonic.

[GeneNetworkBuilder](/packages/GeneNetworkBuilder)
------------------

                       Changes in version 1.47.4                        

- Add nodeData parameter to polishNetwork function.

                       Changes in version 1.47.3                        

- Fix the duplicated parameters in polishNetwork function.

                       Changes in version 1.47.2                        

- Add cy3Network function.

                       Changes in version 1.47.1                        

- Add edge line-width mapping for polishNetwork function.

[GeneTonic](/packages/GeneTonic)
---------

                       Changes in version 2.99.0                        

Other notes

- The transition to the functions available in the mosdef
Bioconductor
is complete, with the original functions now being deprecated. This
applies to goseqTable() (now replaced by mosdef::run_goseq()), which
has now been made faster and more robust in its functionality and in
the ways it can be executed
- The gene plot widgets now also use the gene_plot() function from
mosdef, instead of the previous ggplotCounts() function -
gene_plot() is more flexible and has more options to control the
behavior of the final plot object
- The deseqresult2tbl() and deseqresult2DEgenes() are now replaced by
the more flexible mosdef::deresult_to_df()
- The internally defined createLinkENS(), createLinkGeneSymbol(), and
createLinkGO() are now replaced by the equivalent functions in
mosdef
- The Roxygen-based documentation now supports markdown. No visible
changes should appear to the user, as the content should have stayed
fairly the same
- Although no visible changes for the end user are expected, the
incoming major version bump will reflect the change in the
dependency graph, ensuring that this is noticed at least at the
version numbering level

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.42.0                        

BUG FIXES

- Use more robust heuristic in internal helper
  get_current_Ensembl_release().

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.42.0                        

- No changes in this version.

[GenomicDataCommons](/packages/GenomicDataCommons)
------------------

                       Changes in version 1.30.0                        

New features

- gdc_clinical includes clinical data from the
cases.follow_ups.other_clinical_attributes entity (&#64;LiNk-NY).

Bug fixes and minor improvements

- Removed legacy function, methods, endpoints, and arguments (&#64;LiNk-NY)
- Use native pipe &#124;> instead of magrittr::%>% (&#64;LiNk-NY)

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.58.0                        

NEW FEATURES

- The TxDb class now extends the new OutOfMemoryObject class defined in
  BiocGenerics (virtual class with no slots).

- Calling saveRDS() on a TxDb object now raises an error with a message
  that redirects the user to AnnotationDbi::saveDb().

[GenomicPlot](/packages/GenomicPlot)
-----------

                        Changes in version 1.3.4                        

- Modify pseudo count for ratio matrix

                        Changes in version 1.3.3                        

- Add gene id to rownames to heatmap

                        Changes in version 1.3.2                        

- Fixed a bug in prepare_3parts_genomic_features
- Add messages to each function

                        Changes in version 1.3.1                        

- Updated dependency on R >= 4.4.0 and genomation >= 1.36.0

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.58.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Use hg38 instead of hg19 in vignettes and examples. See commit
  57af07f
  for the details.

BUG FIXES

- The as.data.frame() methods for GenomicRanges and GPos objects
  now obey the 'optional' argument.
  See https://github.com/Bioconductor/GenomicRanges/issues/86

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

                        Changes in version 1.3.1                        

- Helper functions to check DGEdata and a visual unit test

[ggmanh](/packages/ggmanh)
------

                        Changes in version 1.10                         

New Feature

- Added a visualization that is an extension of manhattan plot called
binned manhattan plot. Instead of plotting each individual variant,
the variants are binned based on position (x-axis) and log10 p-value
(y-axis). More information can be found in the vignette.

- Preprocessing function is binned_manhattan_preprocess and plotting
function is binned_manhattan_plot.

Internal Changes

- Change the way positions of the variants are calculated. MPdata now
saves the unscaled positions and the scaling is done in
manhattan_plot function. (this is for the x-axis)

- To get the actual relative positions used in the manhattan plot,
use
calc_new_pos() on the MPdata object.

- Gap between chromosomes can be rescaled by using chr.gap.scaling=
argument inside manhattan_data_preprocess and manhattan_plot
functions.

- Now accounts for negative positions.

[ggsc](/packages/ggsc)
----

                        Changes in version 1.3.1                        

- add plot_lisa_feature() (2024-09-06, #34, #39)
- add geom_bgpoint() layer (2024-06-18, Tue, #24)

[ggtree](/packages/ggtree)
------

                       Changes in version 3.13.2                        

- mv td_filter(), td_unnest() and td_mutate() to 'ggfun' package
(2024-10-24, Thu)

                       Changes in version 3.13.1                        

- layout argument can be a function to re-calculate the coord of
nodes
(2024-07-27, Sat, #622)
- layout.params = list(as.graph=TRUE) (by default) for converting
the tree to a graph tree (igraph object), so that graph layout
can be directly applied to visualize a tree
- otherwise, the function can assume the input tree as a phylo or
treedata object
- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)
- mv %<+% operator to 'ggfun' (2024-05-28, Tue)

[ggtreeSpace](/packages/ggtreeSpace)
-----------

                        Changes in version 1.1.1                        

- update geom_treespace() (2024-09-07)

[ginmappeR](/packages/ginmappeR)
---------

                        Changes in version 1.0.7                        

ENHANCEMENT

- Manuals' examples revised to avoid failing R CMD CHECK on macOS
  platforms.

                        Changes in version 1.0.6                        

ENHANCEMENT

- Test suite revised and enhanced to avoid random failing builds.

                        Changes in version 1.0.5                        

BUG FIX

- UniProt to KEGG translation tests fixed.

                        Changes in version 1.0.4                        

NEW FEATURES

- Added functions that output databases versions being used:
  getCARDVersion()
  getNCBIVersion()
  getUniProtVersion()
  getKEGGVersion()

ENHANCEMENT

- Mapping progress indicator added.

- Improved CARD download messages.

BUG FIX

- Functions that use UniProt.ws printing empty lines before the result
  fixed.

                        Changes in version 1.0.3                        

BUG FIX

- UniProt tests fixed.

[glmGamPoi](/packages/glmGamPoi)
---------

                  Changes in version 1.17 (2024-05-29)                  

- `test_de` can now compute the standard error of the log2-fold change
  (PR#63, thanks @jackkamm)

- `predict` now considers the `ridge_penalty` when calculating the
  standard errors.

- Detect if smart subset of MatrixGenerics' `rowMeans2` and `rowSums2`
  can be used to speed-up `pseudobulk` calculation.

[GloScope](/packages/GloScope)
--------

                 Changes in version 1.3.1 (2024-08-01)                  

- Fix bug which caused the `gloscope_proportion` function to not be
  properly exported.

[GNOSIS](/packages/GNOSIS)
------

                 Changes in version 1.99.0 (2023-09-04)                 

- Made the following significant changes
  o added functionality to select and upload cBioPortal study
  o deprecated ability to save R script with executed code

- Submitted to Bioconductor

[GOSemSim](/packages/GOSemSim)
--------

                       Changes in version 2.31.1                        

- get_rel_df() to access ontology relation data frame required by the
Wang method (2024-08-13, Tue)
- a set of 'OntDb' methods for accessing sqlite data generated by the
'obolite' package (2024-08-13, Tue)
- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)

[goseq](/packages/goseq)
-----

                       Changes in version 1.58.0                        

- Role of maintainer taken over by Federico Marini
- Refreshing the codebase
- introducing roxygen-based documentation
- applying styler for consistent spacing
- using Github Actions for CI/CD
- refined import statements

[goSorensen](/packages/goSorensen)
----------

                        Changes in version 1.7.0                        

We add the following function:

- enrichedIn: This function generates a binary matrix in which the
rows represent GO terms, and the columns represent lists. The matrix
uses the values TRUE or FALSE to indicate whether each GO term is
enriched or not for each list.

In addition:

- We have updated the buildEnrichedTable and allBuildEnrichedTable
functions. They now uses the outcomes of enrichedIn to build the
contingency table. During the conducted tests, the speed of the
procedure for generating contingency tables increased by a factor of
six compared to the previous version.

Furthermore, we have included the showEnrichedIn argument in these
functions. This is a boolean argument. If the argument is TRUE, in
addition to the enrichment contingency table, the function saves a
matrix in the global environment, which contains the cross table of
the enriched and non-enriched GO terms vs the names of the gene
lists generated with the enrichedIn function. This matrix is stored
under the name of "enrichedIn_" followed by the name of the ontology
and the level being analysed. An object will be produced for every
scenario if there are several levels and ontologies.

- We add the vignette "irrelevance-threshold_Matrix_Dissimilarities."
This vignette illustrates calculating, visualizing, and interpreting
the irrelevance-threshold matrix of dissimilarities D. This matrix
provides dissimilarities between pairs of compared lists. These
dissimilarities are not only a descriptive measure but also based on
the irrelevance threshold determining whether two lists are
equivalent. So, this dissimilarity measure between the two lists is
directly associated with their declaration of equivalence.

The matrix D can be represented in interpretable statistic graphs
such as dendrograms or biplots, which help to visualize the
formation of groups containing equivalent lists. Furthermore,
interpreting the biplot dimensions gives us the biological functions
responsible for the equivalence between lists.

[GRaNIE](/packages/GRaNIE)
------

                 Changes in version 1.9.7 (2024-09-25)                  

Bug fixes

- we made a few adjustments for the still experimental GC correction
feature in addConnections_TF_peak() that threw errors for some edge
cases when percBackground_size = 0

                 Changes in version 1.9.6 (2024-09-01)                  

New features

- we now support the mm39 mouse assembly. Please let us know if
something does not work. We are currently also generating the TFBS
resources for mm39, stay tuned. Once done, the links in the vignette
will work.

Improvements

- we now require patchwork version 1.2.0 and above, due to an
incompatibility reported between ggplot2 3.5 and patchwork versions
< 1.2.0. The error that may be thrown in
plotDiagnosticPlots_peakGene is this: Error in
Ops.data.frame(guide_loc, panel_loc) : ‘==’ only defined for
equally-sized data frames

                 Changes in version 1.9.5 (2024-09-01)                  

Improvements

- single-cell vignette update, made it clearer which settings to
change in the GRaNIE workflow for single-cell data
- warning messages now also contain the name of the function they
were
originally called from
- forcats warning for package versions of 1.0.0 and above fixed

                 Changes in version 1.9.4 (2024-08-01)                  

- changed default values for gene.types in filterGRNAndConnectGenes()
to all and adjusted documentation across functions that have also
this argument - by default, there is now no standard filter for the
gene type (genes were filtered to only protein-coding genes before
using protein_coding). This results in bigger eGRNs and seems more
sensible than to automatically filter genes for most use cases.

                 Changes in version 1.9.3 (2024-07-16)                  

Bug fixes

- fixed a bug that caused the peak-gene diagnostic plots to
incorrectly show background signal in plotDiagnosticPlots_peakGene()
(upper right plot on page 1). This was caused by not correctly
shuffling the RNA cunts data due to a code change that occurred in
version 1.7.4. We apologize for the confusion this may have caused.

              Changes in version 1.9.0-1.9.1 (2024-06-03)               

- version jump due to new Bioconductor development cycle

Improvements

- update for the single-cell vignette
- we now provide a public download for the HOCOMOCO v12 TF database
for both human (hg38) and mouse (mm10) that can be used as input for
GRaNIE in the addTFBS() function. We tested this database
extensively and use it as default TF database for our GRaNIE
networks.

Bug fixes

- explicitly requiring tidyr 1.3.0 or above due to the usage of
separate_wider_delim in the code in a recent upgrade
- fixed an accidental reference to AnnotationHub::getAnnotationHub

[GSAR](/packages/GSAR)
----

                       Changes in version 1.40.0                        

- 
  New function ADtest is introduced. It implements a
  nonparametric multivariate test of means based on the HDP
  ranking of samples in the MST and the Anderson-Darling
  statistic.

- 
  New function RADtest is introduced. It implements a
  nonparametric multivariate test of variance based on the radial
  ranking of samples in the MST and the Anderson-Darling
  statistic.

- 
  New function CVMtest is introduced. It implements a
  nonparametric multivariate test of means based on the HDP
  ranking of samples in the MST and the Cramer-Von Mises
  statistic.

- 
  New function RCVMtest is introduced. It implements a
  nonparametric multivariate test of variance based on the radial
  ranking of samples in the MST and the Cramer-Von Mises
  statistic.

[GSVA](/packages/GSVA)
----

                        Changes in version 2.00                         

USER VISIBLE CHANGES

- Added a missing data policy for the GSVA and the ssGSEA methods by
  which missing values in the input expression data either propagate
  throughout calculations (use="everything"), prompt an error
  (use="all.obs"), or are discarded from calculations (use="na.rm").

- The readGMT() function now returns also a GeneSetCollection object,
  and can automatically identify the type of feature/gene identifier
  used in the GMT file, adding it to the metadata of the return object.

- Mapping between different types of feature/gene identifiers between
  gene sets and expression data has been improved.

- Added a new default for the kcdf parameter in the GSVA method,
  kcdf='auto', which tries to select the best value for this parameter,
  according to the input data.

- Added a new sparse regime of the GSVA algorithm, specifically
  tailored for expression data stored in sparse matrices, such as
  single-cell data, to address some of the shortcomings of the original
  (classical) algorithm with this type of data.

- Added two new methods, gsvaRanks() and gsvaScores(), to run the GSVA
  algorithm in two steps that provide additional flexibility to score
  gene sets and extract enrichment information.

- Added a new method gsvaEnrichment() to be used after gsvaRanks(),
  which allows one to extract enrichment information for a particular
  gene set and sample/column, including the possibility to produce an
  enrichment plot using base R or ggplot2 graphics.

- Package messages are now displayed using the cli package, changing
  and removing some of them, as well as adding other ones, to improve
  the information reported to the user.

BUG FIXES

- Removed API legacy convenience features, which keep causing problems,
  as per request of CRAN maintainers.

- Fixed some bugs in the shiny app.

[gypsum](/packages/gypsum)
------

                        Changes in version 1.2.0                        

- Switch to rappdirs to choose the cache directory. This allows
  re-use of the cache with the equivalent Python client.

- Improve handling of unexpected errors with non-JSON bodies.

- Added the unlockProject() utility for administrators.

- Renamed defineTextQuery() to gsc() (for "gypsum search clause")
  and generalized it to other fields in the SQLite file.
  Specifically, it now can filter based on the path to the
  metadata document, the project, asset or version names, the
  identity of the uploader and the timestamp of the upload.

- Added the translateTextQuery() function to convert a
  human-friendly search string into a gypsum.search.clause.

- Renamed searchMetadataText() to searchMetadata(), to reflect
  the more general nature of the search with the new gsc()
  function. The old name is now soft-deprecated.

- Renamed searchMetadataTextFilter() to searchMetadataFilter(),
  to reflect the more general nature of the search. The old name
  is now soft-deprecated. Added arguments to support custom
  names/aliases for the project, asset, version and path columns.
  Deprecated the capability for query= to accept character
  vectors, callers should supply a gypsum.search.clause directly.

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.34.0                        

NEW FEATURES

- Add 'as.vector' argument to h5mread().

SIGNIFICANT USER-VISIBLE CHANGES

- Improvements to coercions from CSC_H5SparseMatrixSeed,
  H5SparseMatrix,
  TENxMatrix, or H5ADMatrix to SparseArray:
  - should be significantly more efficient, thanks to various tweaks
  that
  happened in the SparseArray and Delayed5Array packages;
  - support coercing an object with more than 2^31 nonzero values.

- Coercion from any of the class above to a sparseMatrix derivative
  now fails early if object to coerce has >= 2^31 nonzero values.

- All *Seed classes in the package now extend the new OutOfMemoryObject
  class defined in BiocGenerics (virtual class with no slots).

BUG FIXES

- Fix long standing bug in t() methods for CSC_H5SparseMatrixSeed and
  CSR_H5SparseMatrixSeed objects.

- Replace internal calls to rhdf5::H5Fopen(), rhdf5::H5Dopen(), and
  rhdf5::H5Gopen(), with calls to new internal helpers .H5Fopen(),
  .H5Dopen(), and .H5Gopen(), respectively.
  See commit 31a7e06 for more information.

[hermes](/packages/hermes)
------

                        Changes in version 1.9.1                        

Miscellaneous

- Change maintainer's email address.

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.42.0                        

- the input of `hlaUniqueAllele()` can be a hlaAttrBagClass object or
  a hlaAttrBagObj object

- `hlaAlleleToVCF()` outputs 0 instead of NaN if the location of HLA
  gene
  is unknown or unspecified

- update the C codes according to '_R_USE_STRICT_R_HEADERS_=true' &
  '_R_CXX_USE_NO_REMAP_=true'

[HicAggR](/packages/HicAggR)
-------

                        Changes in version 1.1.2                        

- minor fix in doc for preparePlotgardener

                        Changes in version 1.1.1                        

- HicAggR is now on bioconductor release 3_19
- added a minor fix to ExtractSubmatrix when using GRanges object to
extract, retain only GRanges with seqnames that are present in
hicList
- added preparePlotgardener function to generate a data.frame that
can
be directly used in plotgardener's plotHicTriangle, plotHicRectangle
and plotHicSquare directly

[HilbertCurve](/packages/HilbertCurve)
------------

                       Changes in version 1.99.0                        

- Add a hilbert curve legend showing the orientation of the curve.

[ideal](/packages/ideal)
-----

                       Changes in version 1.99.0                        

Other notes

- The transition to the functions available in the mosdef
Bioconductor
is complete, with the original functions now being deprecated. This
applies to goseqTable() (now replaced by mosdef::run_goseq()), which
has now been made faster and more robust in its functionality and in
the ways it can be executed
- The gene plot widgets now also use the gene_plot() function from
mosdef, instead of the previous ggplotCounts() function -
gene_plot() is more flexible and has more options to control the
behavior of the final plot object
- The deseqresult2tbl() and deseqresult2DEgenes() are now replaced by
the more flexible mosdef::deresult_to_df()
- The internally defined createLinkENS(), createLinkGeneSymbol(), and
createLinkGO() are now replaced by the equivalent functions in
mosdef
- The Roxygen-based documentation now supports markdown. No visible
changes should appear to the user, as the content should have stayed
fairly the same
- Although no visible changes for the end user are expected, the
incoming major version bump will reflect the change in the
dependency graph, ensuring that this is noticed at least at the
version numbering level

[igvShiny](/packages/igvShiny)
--------

               Changes in version 2024-08-29 (2024-08-29)               

- fix issue with loading bed files when app is run with query strings

               Changes in version 2024-08-25 (2024-08-25)               

- switch from Rcurl::url.exists to httr::http_error (Windows
compatibility)

- stop using Amazon S3 URLs by default

               Changes in version 2024-08-16 (2024-08-16)               

- fix issue with VCF files

               Changes in version 2024-08-10 (2024-08-10)               

- fix issue with custom files not working properly

- sync with Bioconductor (3_19 release)

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.11.4 (2024-10-28)                 

- changed testdetectCommunity function for group_by parameter

                 Changes in version 1.11.3 (2024-10-25)                 

- fixed testdetectCommunity function

                 Changes in version 1.11.2 (2024-10-08)                 

- fixed testInteraction and countInteraction tests due to machine
  precision issues

                 Changes in version 1.11.1 (2024-10-03)                 

- fixed BiocNeighbors error for DFrame type data

[immunoClust](/packages/immunoClust)
-----------

                     Changes in version 1.37.12-14                      

- CHANGES
  * improvements in meta.clustering process
  * minor changes in plot.immunoClust

                      Changes in version 1.37.3-11                      

- CHANGES
  * processing time improvements in meta.SON.clustering and
  meta.SON.normalisation

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.40.0                        

BUG FIXES

- Make sure that internal helper coerceToCompressedList() always
  propagates the mcols.

[iSEE](/packages/iSEE)
----

                       Changes in version 2.17.4                        

- Add explicitly a "Stop app" button to close the application (should
be a means to nicely behave in container-spawned instances).
Addresses #630
- Fix text in 'configure' vignette. Closes #626
- Copy information from rowRanges to rowData Closes #637
- Require at least one valid value for atomic, groupable and numeric
columns. Closes 660
- If heatmap matrix has no finite values, set the color range
arbitrarily. Closes #610
- Limit the number of rows in heatmap annotation legends to 10.
Closes
#591
- Fix license. Closes #661.
- Rename some internal constants.
- Add checkbox to fix aspect ratio to 1. Closes #541.
- Add checkbox to hide violin boundaries. Closes #619.

                       Changes in version 2.17.3                        

- Preserve the existing order of rows when receiving a selection in
RowDataTable, ColumnDataTable panels.

                       Changes in version 2.17.2                        

- Export function .selectInputHidden().

                       Changes in version 2.17.1                        

- Fix typo in setMethod(".getContinuousMetadataChoices",
"RowDotPlot",
...).

[iSEEde](/packages/iSEEde)
------

                        Changes in version 1.3.1                        

- Version bump to rebuild classes derived from DotPlot.

[iSEEfier](/packages/iSEEfier)
--------

                        Changes in version 1.2.0                        

- Adding in a new function iSEEmarker(), more focused on finding
marker genes
- A new iSEEinit() version without the DynamicMarkerTable panel
- Expanding the input of iSEEinit() to accept features as data.frame
in addition to vector
- Adding visualization plots to iSEEnrich()

[iSEEhex](/packages/iSEEhex)
-------

                        Changes in version 1.7.1                        

- Version bump to rebuild classes derived from DotPlot.

[iSEEhub](/packages/iSEEhub)
-------

                        Changes in version 1.7.1                        

- Version bump to rebuild classes derived from DotPlot.

[iSEEindex](/packages/iSEEindex)
---------

                        Changes in version 1.3.2                        

- Added a full implementation of the runr resource class, defining
its
behavior with simple heuristics based on the fact that users can now
also provide not just a path but the call to the command of R to be
run (therefore, the name runr) to obtain an object to explore. A
typical use case would be to deploy a collection of datasets from a
data package.

                        Changes in version 1.3.1                        

- Version bump to rebuild classes derived from DotPlot.

[iSEEpathways](/packages/iSEEpathways)
------------

                        Changes in version 1.3.1                        

- Version bump to rebuild classes derived from DotPlot.

[iSEEu](/packages/iSEEu)
-----

                       Changes in version 1.17.1                        

- Version bump to rebuild classes derived from DotPlot.

[IsoBayes](/packages/IsoBayes)
--------

                        Changes in version 1.2.5                        

- inference function extended (map_iso_gene can be either be a
  data.frame or a .csv)

                        Changes in version 1.2.3                        

- plot_traceplot added

                        Changes in version 1.2.2                        

- bug in "input_data" fixed (only triggered when using PEP with MM
  data)

                        Changes in version 1.2.1                        

- normalization constant added (accounting for N detected peptides)

[isomiRs](/packages/isomiRs)
-------

                       Changes in version 1.33.1                        

- fix isoAnnotate error due to duplicates lines
  https://github.com/lpantano/isomiRs/issues/18

[KEGGREST](/packages/KEGGREST)
--------

                       Changes in version 1.46.0                        

BUG FIXES

- 1.45.1 Fix keggFind URL to use '+' instead of spaces.

[lefser](/packages/lefser)
------

                       Changes in version 1.15.10                       

- &#91;Major&#93; Name of the two arguments for lefser function is changed
from groupCol and blockCol to classCol and subclassCol,
respectively.
- &#91;Major&#93; Defunct expr argument in lefser
- &#91;New function&#93; lefserPlotFeat plots the histogram of relative
abundance (in the (0,1) interval) of the selected features
- &#91;New function&#93; lefserPlotClad draws the cladogram of the
significantly more abundant taxa and their LDA scores
- &#91;New function&#93; lefserClades runs the lefser, returning additional
information (e.g., agglomerates the features abundance at different
taxonomic ranks) required for lefserPlotClad.
- &#91;New feature&#93; Visualization functions are using a color-blind
friendly color palette by default.

                       Changes in version 1.15.7                        

- &#91;Major algorithm update&#93; We remove the step (createUniqueValues) in
the lefser function, which used to add small random numbers to make
all the values unique. Potential issues (e.g., LDA) due to excess 0s
should be managed by filtering out low abundant features from the
input.

                       Changes in version 1.15.3                        

- The column names of lefser output is changed to c("features",
"scores") from c("Names", "scores")
- &#91;New feature&#93; The get_terminal_nodes function to select only the
terminal nodes of the hierarchical features (e.g., taxonomic data).
- &#91;Major algorithm update&#93; Add an option for adjusting the first
Kruskal-Wallis Rank Sum Test and Wilcoxon-Rank Sum test for multiple
hypothesis testing through the new argument method in the lefser
function.

[lemur](/packages/lemur)
-----

                         Changes in version 1.3                         

- lemur now automatically inserts the variables from the design
formula into the group_by argument for find_de_neighborhoods.
(thanks Katha for pushing for this feature)
- The formula parsing automatically detects global variables and adds
them to the colData. This avoids problems with the random test /
training assignment.
- Duplicate column names in colData are now longer allowed.
- Require harmony version >= 1.2.0 (thanks Maija for reporting the
problem)

[limma](/packages/limma)
-----

                       Changes in version 3.62.0                        

- 
  New function fitFDistUnequalDF1() for improved empirical Bayes
  hyperparameter estimation. The new function is an alternative
  to the long-standing limma functions fitFDist() and
  fitFDistRobustly().  The new function gives special attention
  to contexts where the df1 argument is unequal between cases or
  can include small fractional values, such as those resulting
  from edgeR's voomLmFit() or from the edgeR 4.0 quasi-likelihood
  pipeline.  The new function uses profile maximum likelihood
  instead of moment estimation to estimate df.prior, and so may
  be more accurate than fitFDist() or fitFDistRobustly() even
  when the df1 values are all equal. The new hyperpameter
  estimation will be used by default by eBayes(), treat() and
  squeezeVar() whenever the residual degrees of freedom are
  unequal, but users can request legacy behavior for those
  functions by setting `legacy=TRUE`.

- 
  New argument `adaptive.span` for voomWithQualityWeights() and
  normalizeCyclicLoess(), with the same purpose as the same
  argument for voom(). If `TRUE` then an optimal value for the
  span is chosen based on the number of genes in the dataset.
  voom() and voomWithQualityWeights() now return the chosen
  `span` as part of the output object if `adaptive.span=TRUE`.

- 
  contrasts.fit() now converts the coefficient name "(Intercept)"
  created automatically by model.matrix() to "Intercept", same as
  is done by makeContrasts(). This should avoid inconsistent
  coefficient names between the fitted model object and contrast
  matrices created by makeContrasts().

- 
  New argument `keep.EList` for voomaLmFit(), similar to the same
  argument for edgeR::voomLmFit(). The observation weights
  generated by voomaLmFit() are now stored only if
  `keep.EList=TRUE`.

- 
  voomaByGroup() now uses the overall mean-variance trend for
  groups with only one sample, for which a separate trend is not
  estimable.

- 
  Add Ravindra et al (2023) as a reference and Mengbo Li as an
  author to the voomaByGroup() help page.  Emphasize in the help
  details that only simple design matrices are supported by this
  function.

- 
  Revise the eBayes/treat help page to further clarify the
  `legacy` and `fc` arguments.

- 
  Bug fix to voomaByGroup() to save variance plot parameters
  correctly when `save.plot=TRUE`.

- 
  Remove argument `zero.weights` from the MArrayLM method for
  plotMD() and plotMA().

- 
  Fix bug in the MArrayLM methods for plotMD() and plotMA(),
  which were attempting to access a `weights` component of the
  fitted model object, analogously to what plotMD() does with the
  EList or MAList objects, even though MArrayLM objects do not
  contain a `weights` component.

- 
  Use of the stats package digamma() function replaced with
  statmod's logmdigamma() in fitFDist(), fitFDistUnequalDF1() and
  intraspotCorrelation().  Slight improvements to speed, accuracy
  and code simplicity, but results should remain the same to at
  least 13 significant figures.

[Macarron](/packages/Macarron)
--------

                        Changes in version 1.9.1                        

- Synced Version with Bioconductor

- Added citation

- Updated README

- Changed heatmap plotting default to FALSE

- Corrected curly bracket placement in wrapper

[MACSr](/packages/MACSr)
-----

                       Changes in version 1.13.1                        

- Upgrade to MACS 3.0.2.

[maftools](/packages/maftools)
--------

                       Changes in version 2.22.0                        

BUG FIXES

- Bug fix in using keepGeneOrder in coOncoplot(). Issue: 1061
- Bug fix in using selectedPathways in oncoplot(). Issue: 1041
- Add an error message when bai files are missing sampleSwaps().
Issue: 1028
- Bug fix in tmb while handling multiple MAFs. Issue: 1018
- Handle missing NAs while sub-setting for ranges. Issue: 1013
- Better error handling when zero mutated samples are encountered in
clinicalEnrichment. Issue: 1010
- MAJOR: read.maf by default coerces clinical data columns to
character. This bug fix avoids it and is auto detected. Issue: 997

ENHANCEMENTS

- Better handling of color codes for continuous variable annotations
1053
- Add right_mar to gisticOncoPlot1043
- Added PPDPFL to protein domain database Issue: 1025
- Better sorting of oncoplot with collapsePathway
- Changed default background for oncoplot from gray to #ecf0f1
- Changed default signature database to SBS_v3.4 (from legacy)
- Update tmb function

BREAKING CHNAGES

- Column order required for pathways() function changed from
Pathway,Gene to Gene,Pathway. Issue: 1041

NEW FUNCTIONS

- gisticCompare() for comparing two GISTIC objects
- segSummarize() for summarizing DNAcopy segments

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.13.5 (2024-06-27)                 

- add log10 as transformation method in transformAssay

- add "none" as imputation method in imputeAssay

- update shiny application with log10

- update vignette with log10

- update unit tests with log10 and "none"

                 Changes in version 1.13.4 (2024-06-24)                 

- update vignette, add information on ComBat

                 Changes in version 1.13.3 (2024-06-21)                 

- add functionality to perform batch correction using sva::ComBat

- improve visualisation of MAplot in shinyQC

                 Changes in version 1.13.2 (2024-06-11)                 

- use DT::renderDT instead of shiny::renderDataTable

- use DT::DTOutput instead of shiny::dataTableOutput

                 Changes in version 1.13.1 (2024-04-25)                 

- add batch2 and ... arguments in batchCorrectionAssay

[matter](/packages/matter)
------

                 Changes in version 2.7.10 (2024-10-16)                 

NEW FEATURES

- Add 'fetch' and 'flash()' for base R types

SIGNIFICANT USER-VISIBLE CHANGES

- Update 'fetch()' and 'flash()' to pass more
  arguments to underlying 'matter' constructors

                        Changes in version 2.7.9                        

NEW FEATURES

- Add 'pivots' argument to 'fastmap()' to support
  more than 2 pivot candidates per iteration

SIGNIFICANT USER-VISIBLE CHANGES

- Update 'fetch()' and 'flash()' to preserve atoms types

- Update 'fetch()' and 'flash()' to fall back to serial
  execution if BPPARAM is a distributed cluster

BUG FIXES

- Add 'chunkopts' to 'nscentroids()' and 'sgmix()'
  signatures to fix downstream '...' bugs

- Fix graphical parameters not applying to vizi plots

                        Changes in version 2.7.8                        

BUG FIXES

- Fix shared memory segmentation faults on Linux

- Fix voxel plotting when length(unique(z)) == 1L

                        Changes in version 2.7.7                        

NEW FEATURES

- New 'fetch()' and 'flash()' generics for moving
  data between file storage and shared memory

SIGNIFICANT USER-VISIBLE CHANGES

- Support both number of items OR size in bytes when
  setting getOption("matter.default.chunksize")

- Now using '@' prefix as shared memory identifier
  for file system compatibility (e.g., Windows)

- Deprecated 'sgmixn()'; updated 'sgmix()' to support
  multichannel images directly

- Improved RNG behavior for non-L'Ecuyer seeds

BUG FIXES

- Fix potential infinite loop in 'mi_learn()'

- Add warning in chunk-apply functions if 'RNG=TRUE'
  but 'RNGkind()' is _not_ "L'Ecuyer-CMRG"

                        Changes in version 2.7.6                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add 'cpals()' and 'dpals()' to list available palettes

- Moving 'simple_logger' now appends 'sessionInfo()'

- Statistical methods 'fastmap()' and 'nscentroids()'
  now use 'rowDists()' and 'colDists()' directly

- Optional argument 'distfun' has new requirements
  in methods 'fastmap()' and 'nscentroids()'

- Class 'matter_str' now inherits from 'matter_list'
  for consistency between R and C++ interfaces

- Export '.rowDists()' and '.colDists()' functions

BUG FIXES

- Fix 'simple_logger' finalizer on R session exit

                        Changes in version 2.7.5                        

NEW FEATURES

- Shared memory support comes to 'matter' via boost!

- Add 'as.shared()' for coercion to shared memory 'matter' objects

- Use 'path=":memory:"' to construct 'matter' object in shared memory

- Use 'mem()' to monitor shared memory usage

SIGNIFICANT USER-VISIBLE CHANGES

- Add shared memory support to 'mem()' output

- Add cluster memory support to 'memtime()' output

- Add cluster memory monitoring with 'memcl()' function

- New documentation for several previously-undocumented
  utility functions that may be useful for others

- New 'permute' argument in chunk-apply functions

BUG FIXES

- Improved behavior for the general 'matter()' constructor

- Use private RNG for 'uuid()' to avoid collisions
  from 'set.seed()' and to avoid changing '.Random.seed'

- Fix 'is_gridded()' and 'to_raster()' behavior
  for non-integer coordinates and character rasters

                        Changes in version 2.7.4                        

NEW FEATURES

- Add 'SnowfastParam' class to support PSOCK clusters
  for faster startups and faster serialization

- Add getOption("matter.default.chunksize") for setting
  the approximate size of chunks in bytes

- Add getOption("matter.default.serialize") for managing
  serialization behavior in chunk-apply functions

- New 'simple_logger' class for writing log files

- New 'options(matter.temp.dir=tempdir())' replaces
  the previous 'options(matter.dump.dir=tempdir())'

- New 'options(matter.temp.gc=TRUE)' controls whether
  temporary matter files are automatically removed when
  all R objects referencing them are garbage collected

- Temporary matter files are now garbage collected!

SIGNIFICANT USER-VISIBLE CHANGES

- New 'chunkopts' argument replaces 'nchunks' and includes
  the new chunk options ("chunksize" and "serialize")

BUG FIXES

- Fix leaky closures in internal chunk loop functions
  causing slow parallel performance on large datasets

- Fix internal bug in 'eval_exprs()' where 'recursive'
  attribute was always set to 'TRUE'

                        Changes in version 2.7.3                        

NEW FEATURES

- Add new functions 'knnmax()' and 'knnmin()'

- Add new functions 'filtn_ma()', 'filtn_conv()', etc.

- Add support for N-dimensional signals to 'estnoise_*()' functions

- Add 'findpeaks_knn()' for KNN-based peak detection

SIGNIFICANT USER-VISIBLE CHANGES

- Add signal processing vignette

- Faster 'knnsearch()' when the data and query are the same

- Ties in 'knnsearch()' are now broken based on data order

- Updated 'estnoise_sd()', 'estnoise_mad()', and 'estnoise_quant()'

BUG FIXES

- Improvements to 'filt1_adapt()' and 'filt2_adapt()'

                        Changes in version 2.7.2                        

NEW FEATURES

- Support 2D signals in 'plot_signal()' (e.g., for ion mobility)

- New argument 'alphapow' in 'plot_signal()' and 'plot_image()'

SIGNIFICANT USER-VISIBLE CHANGES

- Plot method for 'vizi_bars' now sums duplicate bar counts

BUG FIXES

- Fix chunk-apply functions for BPPARAM=NULL

- Fix 'chunkMapply()' named list errors

- Fix 'is_gridded()' so infinite 'estres()' result gives FALSE

- Give more useful error messages for 'dpal()' and 'cpal()'

- Fix alpha channel in 'plot_image()'

                        Changes in version 2.7.1                        

NEW FEATURES

- Support subsetting of deferred operations

- Add 'chunked' classes for chunked objects

- Add 's_stat()' for grouped streaming statistics

- For 'plot_image()', set 'useRaster' to plot 3d surface/volume vs
  voxels

- Add 'enhance_adj()' function for adjusting extreme values

SIGNIFICANT USER-VISIBLE CHANGES

- Performance improvements for chunk-apply functions

- Support chunk-apply functions on remote cluster workers

- Reduced memory usage for chunk-apply functions on SNOW clusters

- Image smoothing 'filt2_*()' functions support multiple channels

- Contrast enhancement 'enhance_*()' functions support multiple
  channels

- Remove 'biglm' from imports

BUG FIXES

- Fix subplot axis scaling and aspect ratio for 'plotly' graphics

- Fix passing colors to line charts for 'plotly' graphics

                        Changes in version 2.6.3                        

BUG FIXES

- Fix 'filt2_gauss()' behavior with NAs

                        Changes in version 2.6.2                        

BUG FIXES

- Fix 'sgmixn()' on 'matter_mat' and 'sparse_mat' input

                        Changes in version 2.6.1                        

BUG FIXES

- Fix 'matter_list' error on zero-length elements

[memes](/packages/memes)
-----

                       Changes in version 1.13.1                        

- Fixed an error in importAme that prevented import when runAme was
run with sequences = TRUE (Reported by @withermatt on Github. Thank
you!)

[meshes](/packages/meshes)
------

                       Changes in version 1.31.1                        

- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)

[metabCombiner](/packages/metabCombiner)
-------------

                       Changes in version 1.15.2                        

- Bug fixes for duplicate column names

                       Changes in version 1.15.1                        

- New m/z shift modeling functionality added:
  + mzFit(): m/z shift modeling and plotting function
  + mzfitParam(): parameter list function for listing modeling
  parameters

- update to calcScores():
  + new arguments (mzshift and mzfit) for m/z correction before scoring

- new citation

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                         Changes in version 1.9                         

Changes in 1.9.2

- Fix missing ProtGenerics dependency.

Changes in 1.9.1

- Add parameter scalePeaks to plotSpectraMirror to allow scaling peak
intensities before plotting.

[metabolomicsWorkbenchR](/packages/metabolomicsWorkbenchR)
----------------------

                       Changes in version 1.15.1                        

- fix untarg import

[methyLImp2](/packages/methyLImp2)
----------

                 Changes in version 1.1.1 (2024-09-10)                  

- Fixing wrong error/warning messages;

- Adding example of user-provided annotation.

                 Changes in version 1.0.1 (2024-07-03)                  

- Fixing warnings printing.

[mia](/packages/mia)
---

                        Changes in version 1.13                         

- Added new functions getMediation and addMediation

- replace getExperiment* and testExperiment* functions with
  getCrossAssociation

- Replace mergeRows and mergeCols with new function
  agglomerateByVariable

- agglomerateByRank: bugfix related to agglomeration of non-existing
  tree

- getHierarchyTree: bugfix related to empty cells in rowData

- agglomerateByRank: bugfix: trees was not pruned correctly

- rename meltAssay to meltSE

- Rename countDominantFeatures and countDominantTaxa to
  summarizeDominance

- rename subsetByPrevalent* and subsetByRare* to subsetByPrevalent and
  subsetByRare

- Deprecate transformFeatures, transformSamples, transformCounts,
  Ztransform,
  relAbundanceCounts

- rename getRare* functions to getRare, getUnique* functions to
  getUnique,
  getTop* functions to getTop and getPrevalent* functions to
  getPrevalent

- Rename subsampleCounts to rarefyAssay

- Rename perSampleDominant* functions to getDominant and
  addPerSampleDominant*
  functions to addDominant

- Rename splitByRanks to agglomerateByRanks and add option as.list

- Explain that rarefyAssay returns a new SummarizedExperiment object
  that
  includes the newly added subsampled assay.

- Fix bug in mergeFeaturesByPrevalence

- new aliases calculateDPCoA to getDPCoA, calculateNMDS to getNMDS,
  calculateRDA to getRDA,
  calculateCCA to getCCA

- add informative error message in rarefyAssay on assays with
  strictly-negative values

- Use rbiom package in unifrac implementation

- Updated parameter names to follow naming convention "parameter.name"

- rename converters makeTreeSEFrom* to convertFrom* and
  makePhyloseqFromTreeSE to
  convertToPhyloseq

- add rowTree agglomeration and RefSeq agglomeration in
  agglomerateByPrevalence

- Fix tree merging in unsplit and mergeSEs functions

- Added addAlpha; a wrapper for calculating all alpha diversity indices

- Added importTaxpasta

- Changes in default taxonomy ranks; more ranks supported

- Added Tito2024QMP dataset

- Added convertToBIOM

- new methods getLDA and addLDA for LDA ordination with feature
  loadings
  computation

- new methods getNMF and addNMF for NMF ordination with feature
  loadings
  computation

- If missing values, give informative error in *RDA/*CCA functions

- transformAssay can apply transformation to altExp

- Added CSS transformation

- In agglomerateByVariable, splitOn and getDominant, use 'group' to
  specify grouping variable.

[miaViz](/packages/miaViz)
------

                        Changes in version 1.13                         

- plot*Tree: bugfix, ununique nodes

- Added confidence.level parameter to plotCCA

- add plotNMDS to miaViz (plotNMDS deprecated in mia)

- Added getNeatOrder function

- plotAbundance: enable plotting without agglomeration

- Change parameter naming convention from parameter_name to
  parameter.name

- Added plotLoadings function

[MicrobiomeProfiler](/packages/MicrobiomeProfiler)
------------------

                       Changes in version 1.11.1                        

- support gson for PubChem Pathway (2024-08-16, Fri, #6)

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.17.1                        

- update scale_fill_diff_cladogram to be compatible with ggnewscale
>= 0.5.0 (2024-07-24, Wed)

                       Changes in version 1.17.0                        

- Bioconductor 3.19 release, and Biocondutor 3.20 (devel) version
bump. (2024-05-14, Tue)

[miloR](/packages/miloR)
-----

                 Changes in version 2.0.1 (2024-04-30)                  

- Introduce NB-GLMM into Milo 2.0 for random effect variables and
modelling dependencies between observations
- Diagnostic function for checking model separation for experimental
variables, i.e. splitting zero from non-zero counts perfectly
- Vignette describing basic usage of GLMM functions in testNhoods

[MIRit](/packages/MIRit)
-----

                        Changes in version 1.1.1                        

Changes were made to ensure that R CMD check runs without errors or
warnings. Moreover, MIRit has been updated to use the new miRTarBase
version 10 database.

[miRSM](/packages/miRSM)
-----

                        Changes in version 2.1.2                        

- Update CITATION <2024-09-18, Wed>

                        Changes in version 2.1.1                        

- Update module_FA function <2024-08-25, Sun>

[miRspongeR](/packages/miRspongeR)
----------

                        Changes in version 2.9.1                        

- Update moduleDEA function <2024-08-25, Sun>.

[monaLisa](/packages/monaLisa)
--------

                       Changes in version 1.11.3                        

- add show_bin_legend argument to plotMotifHeatmaps (contributed by
@danymukesha, PR #62)

[MOSClip](/packages/MOSClip)
-------

                        Changes in version 1.0.0                        

- Initial BioConductor submission.
- New topological pathway analysis tool able to integrate multi-omics
data. It finds survival-associated gene modules or significant
modules for two-class analysis. This tool have two main methods:
pathway tests and module tests. The latter method allows the user to
dig inside the pathways itself.

[mosdef](/packages/mosdef)
------

                        Changes in version 1.2.0                        

- New functionality added for general scope tasks, such as feature
annotation and summaries on the expression tables
- The gene_plot() function now is (again) able to handle the
specification of more than one grouping variable

[MOSim](/packages/MOSim)
-----

                 Changes in version 2.1.0 (2024-05-09)                  

- Added TF simulation to scMOSim

[motifStack](/packages/motifStack)
----------

                       Changes in version 1.49.1                        

- Fix the issue for importMatrix when 'w' is missing.

[motifTestR](/packages/motifTestR)
----------

                        Changes in version 1.1.5                        

- Added clusterMotifs, testClusterPos and testClusterEnrich
- Enforced strict use of PWMs for all functions

[MouseFM](/packages/MouseFM)
-------

                 Changes in version 1.14.1 (2024-08-26)                 

- Backend server changed

[mpra](/packages/mpra)
----

                       Changes in version 1.27.2                        

- Add endomorphic option to return MPRASet as a final result

- Add library size default option for normalization

                       Changes in version 1.26.1                        

- Pass in block vector in normalize_counts() call in .fit()
  functions

[msa](/packages/msa)
---

                       Changes in version 1.37.4                        

- change of interface of msaPrettyPrint() function in order to allow
  for a
  workaround for texi2dvi() problems (using quiet=TRUE and index=FALSE
  together
  throws an error)

                       Changes in version 1.37.3                        

- fix of major bug in the Muscle interface

                       Changes in version 1.37.2                        

- further fix of ClustalOmega makefile and source code to ensure
  compatibility
  with new Rcpp version (further adaptations for MacOS platform)

                       Changes in version 1.37.1                        

- fix of ClustalOmega makefile and source code to ensure compatibility
  with
  new Rcpp version

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                        Changes in version 1.13                         

Changes in 1.13.1

- Add extractByIndex() method.

[MsBackendMetaboLights](/packages/MsBackendMetaboLights)
---------------------

                        Changes in version 0.99                         

Changes in 0.99.1

- Add mtbls_sync() to synchronize a MsBackendMetaboLights object.
- Add mtbls_sync_data_files() function to cache selected files and
mtbls_cached_data_files() to list all locally cached data files.

Changes in 0.99.0

- Prepare package for submission to Bioconductor.

[MsBackendSql](/packages/MsBackendSql)
------------

                         Changes in version 1.5                         

Changes in 1.5.1

- Implement the new extractByIndex() methods.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                        Changes in version 1.17                         

MsCoreUtils 1.17.3

- Use R_Calloc and R_Free instead of Calloc and Free in
src/lowerConvexHull.c, respectively, to reflect changes in the R API
and fullfil STRICT_R_HEADERS check.

MsCoreUtils 1.17.2

- Fix typo in normalisation methods description.

MsCoreUtils 1.17.1

- Fix common_path() to not return also the file name if the input
parameter contains only identical character strings.

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.31                         

MSnbase 2.31.2

- Suggest pRolocdata (>= 1.43.2.1) (that has some extdata, needed to
other packages' vignettes).

MSnbase 2.31.1

- Disable nested parallel processing for chromatogram() method.
- Fix Rd notes.

                       Changes in version 2.31.0                        

- New Bioconductor devel.

[msPurity](/packages/msPurity)
--------

                       Changes in version 1.30.1                        

- Update example for purityX causing build error

[MSstatsConvert](/packages/MSstatsConvert)
--------------

                       Changes in version 1.16.0                        

- Improved memory management for efficiency.
- Added support for ProteinProspector data conversion to MSstatsTMT
format.
- Updated NA handling for Skyline and Spectronaut converters
- Added MetaMorpheus converter functionality.

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.32.0                        

Bug fixes and minor improvements

- Various documentation improvements to MultiAssayExperimenToMAF,
saveHDF5MultiAssayExperiment, and reexports
- When rownames are numeric characters e.g., "1", ensure they stay
character when converting MultiAssayExperiment to longFormat.

                       Changes in version 1.30.2                        

Bug fixes and minor improvements

- The colData<- replacement method now correctly works with
data.frame
value inputs (@drighelli, #330).
- Updated CITATION information in the main vignette.
- Use reshape2::melt instead of stats::reshape to preserve row names
in longFormat

[multiGSEA](/packages/multiGSEA)
---------

                 Changes in version 1.15.1 (2024-10-22)                 

- Remove proteomics feature with empty identifier to prevent
  errors in vignette building.

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.13.7                        

Bug fix

- infer_eff_direction now includes A0 as an ambiguous case as well as
A1/A2.

New features

- eff_on_minor_alleles parameter added (off by default) - controls
whether MungeSumstats should assume that the effects are
majoritively measured on the minor alleles. Default is FALSE as this
is an assumption that won't be appropriate in all cases. However,
the benefit is that if we know the majority of SNPs have their
effects based on the minor alleles, we can catch cases where the
allele columns have been mislabelled.

                       Changes in version 1.13.6                        

New features

- Mappings added to mapping file for risk and non risk allele.

                       Changes in version 1.13.3                        

Bug fix

- Bug fix for check 3 in infer effect column - previously A1 & A2
were
swapped when there were more matches for the ref genome in A1 rather
than A2 which was incorrect. Corrected now so it will only be
flipped when A2 has more matches to the reference genome.

                       Changes in version 1.13.2                        

New features

- Handling of -log10 p-values (outside of VCFs) added.

                       Changes in version 1.13.1                        

New features

- Mapping for OA (other Alllele) added to A1.

[muscat](/packages/muscat)
------

                       Changes in version 1.19.1                        

- added J Gilis, D Risso, L Clement as authors

- differential detection with 'pbDS(..., method="DD")' or 'pbDD()'
  & stagewise testing &#91;Gilis et al.&#93;, plus corresponding vignette

- replace 'aes_string()' in 'ggplot()' by '.data$.' from 'rlang'

[musicatk](/packages/musicatk)
--------

                 Changes in version 1.14.1 (2024-07-13)                 

- Removed deconstructSigs

[mzR](/packages/mzR)
---

                       Changes in version 2.39.2                        

- Improve openMSfile and openIDfile manual page (contirbuted by
  rdhale92, PR #299)

                       Changes in version 2.39.1                        

- Update CV translator, adds support for Astral

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 3.1.0                        

- Fixed parsing error for BAM files leading to an extre site to be
called whenever skip width is greater than 0.
- Added functions get_cgi_*() to get CpG islands annotation for mm10,
GRCm39, hg19, and hg38.
- Changed colour palette for heatmaps to have a less bright yellow on
the low end.
- Changed default NanoMethViz.site_filter value to 3. This filters
out
any sites with less than 3 coverage.
- Changed BAM file parsing behaviour when parse flag "." or "?" is
not
set to default to "?". This is against SAM spec but follows the
behaviour of IGV and older PacBio datasets.

[nipalsMCIA](/packages/nipalsMCIA)
----------

                 Changes in version 1.2.1 (2024-08-31)                  

Major changes

Minor improvements and bugfixes

- Updated 'Single-Cell-Analysis' vignette to remove potential issue
with BiocFileCache having two copies of data.
- Updated author list and citation.
- Updated references in readme.

[NormalyzerDE](/packages/NormalyzerDE)
------------

                       Changes in version 1.23.2                        

- Fix bug where non-categorical ANOVA setting wasn't applied

                       Changes in version 1.23.1                        

- Remove unnecessary dependency (Biobase) by calling rowMedians with
  matrixStats::rowMedians instead

- Fix deprecation warning by using guide = "none" instead of FALSE

[omicsViewer](/packages/omicsViewer)
-----------

                        Changes in version 1.10                         

new features

- dynamic heatmap - zoomed in heatmap upon subsetting features and
  observations

[ompBAM](/packages/ompBAM)
------

                 Changes in version 1.9.1 (2024-07-27)                  

- Remove dependency on zlibbioc

[OncoSimulR](/packages/OncoSimulR)
----------

                 Changes in version 4.7.8 (2024-10-09)                  

- Remove support for maOS with arm64 (i.e., kjohnson3
  and kjohnson1 machines).

                 Changes in version 4.7.7 (2024-09-18)                  

- Try to work around (and ask for help) for tests that fail
  in the BioC kjohnson3 machine (maOS 13.6.5, arm64). Yes,
  more of this.

                 Changes in version 4.7.6 (2024-09-16)                  

- Try to work around (and ask for help) for tests that fail
  in the BioC kjohnson3 machine (maOS 13.6.5, arm64). Yes,
  more of this.

                 Changes in version 4.7.5 (2024-09-12)                  

- Try to work around (and ask for help) for tests that fail
  in the BioC kjohnson3 machine (maOS 13.6.5, arm64).

                 Changes in version 4.7.4 (2024-07-30)                  

- Try to work around (and ask for help) for a vignette example
  that fails in the BioC kjohnson3 machine (maOS 13.6.5, arm64).
  (Yes, a fourth place where it fails).

                 Changes in version 4.7.3 (2024-07-26)                  

- Try to work around (and ask for help) for a vignette example
  that fails in the BioC kjohnson3 machine (maOS 13.6.5, arm64).
  (Yes, a third place where it fails).

                 Changes in version 4.7.2 (2024-07-25)                  

- Try to work around (and ask for help) for a vignette example
  that fails in the BioC kjohnson3 machine (maOS 13.6.5, arm64).
  (Yes, a second place where it fails).

                 Changes in version 4.7.1 (2024-07-23)                  

- Try to work around (and ask for help) for a vignette example
  that fails in the BioC kjohnson3 machine (maOS 13.6.5, arm64).

- Miscell minor: incorporate 361ccdc to 1bc3cbd (in main github
  repo) to the BioC version.

- Fixed Note "Lost braces in \itemize; meant \describe ?"

- Updated vignette address to IIBM new name.

[ontoProc](/packages/ontoProc)
--------

                       Changes in version 1.99.4                        

- caching the RDS version of ontology_index so setup_entities2() is
  fast after first call

- exploring role for bioregistry from pypi -- added a
  bioregistry_ols_resources function

                       Changes in version 1.27.1                        

- added search_labels(), to use owlready2 search facility on term
  labels

[orthos](/packages/orthos)
------

                        Changes in version 1.3.1                        

- Environment installation addresses Anaconda licencing changes

                        Changes in version 1.3.0                        

- Bioconductor Release update

[pathlinkR](/packages/pathlinkR)
---------

                       Changes in version 1.1.22                        

- Expanded and clarified error messages for pathwayEnrichment input
- Updated vignette

                       Changes in version 1.1.18                        

- Fixed bug for pathwayEnrichment when running fgsea (NA gene names)

                       Changes in version 1.1.17                        

- Fix for default parameters and theme settings in ppiPlotNetwork
- Support for GSEA (via package fgsea), both as an option in
pathwayPlots and pathwayEnrichment
- Global fontSize argument added to pathwayPlots

                       Changes in version 1.1.10                        

- Support for KEGG pathways when running pathwayEnrichment, either
through Sigora or traditional over-representation

                        Changes in version 1.1.8                        

- pathwayPlots now supports the output from fgsea on Reactome
pathways
(109a173)

                        Changes in version 1.1.7                        

- Size of asterisks in pathwayPlots is set based on tringle size
(1f9b751)

                        Changes in version 1.1.6                        

- Updated installation information in README
- "eruption" gained argument "labelCutoffs" to add labels to p value
and fold change cutoff lines
- "ppiPlotNetwork" checks for argument "labelColumn" when
"label=TRUE"
- Fixed some tests

[pcaExplorer](/packages/pcaExplorer)
-----------

                       Changes in version 2.99.0                        

New features

- The pcaplot() function now provides a clever default for the
intgroup parameter, if some content (as it should) is provided in
the colData slot of the main input object

Other notes

- The transition to the functions available in the mosdef
Bioconductor
is complete, with the original functions now being deprecated. This
applies to topGOtable() (now replaced by mosdef::run_topGO())
- The gene plot widgets now also use the gene_plot() function from
mosdef, instead of the previous undocumented internal function
- The Roxygen-based documentation now supports markdown. No visible
changes should appear to the user, as the content should have stayed
fairly the same
- Although no visible changes for the end user are expected, the
incoming major version bump will reflect the change in the
dependency graph, ensuring that this is noticed at least at the
version numbering level

[Pedixplorer](/packages/Pedixplorer)
-----------

                        Changes in version 1.1.5                        

- Change code of ped_to_legdf
- When plotting with the main plot, the legend gets its own space
separate from the plot. This allow better control over the size and
localisation of the legend.
- The graphical parameters are reset after each use of plot_fromdf
- Add tooltips control in Pedigree plots and add it to the app
- Add example of interactivness in vignette
- Fix plot area function and legend creation for better alignment

                        Changes in version 1.1.4                        

- Update website and logo
- Improve ped_shiny() esthetics
- Change plot element order rendering for better looks
- Add more control to line width of box and lines
- Improve legend ordering
- Separate website building workflow from check
- Update function documentation and set to internal all unnecessary
functions for users
- Stabilize unit test
- Standardize the vignettes and add more documentation
- Fix label adjusting position in plot functions

                        Changes in version 1.1.3                        

- Fix github workflows
- Disable ped_shiny() execution in markdown
- Publish with pkgdown

                        Changes in version 1.1.2                        

- Use R version 4.4 and update workflows

                        Changes in version 1.1.1                        

- A shiny application is now available through the ped_shiny()
function.
- Function imports have been cleaned.
- Unit tests have been added as well as more snapshot to increase
package coverage.
- relped dataset allows to easily test special relationship.
- Documentation is enhanced and correctly linted.
- precision parameter has been added to align4() and set_plot_area()
to reduce noise between platform.
- fix_parents() has been fixed and improved.
- More controls over color setting with generate_colors().
- Possibility to force computation of alignement when it fails with
force = TRUE.
- upd_famid_id() to upd_famid().
- Zooming in a pedigree object is now done by subsetting the
dataframe
computed by ped_to_plotdf().
- useful_inds() function has been improved.

[pgxRpi](/packages/pgxRpi)
------

                 Changes in version 1.1.8 (2024-10-11)                  

- Added support for accessing and visualizing level-specific CNV
frequency data.
- Enabled calculation of level-specific CNV frequencies from segment
data.

                 Changes in version 1.1.7 (2024-09-11)                  

- Adapted to Progenetix API change: updated endpoint from
"analyses/?output=cnvstats" to "services/cnvstats/".
- The dataset parameter in pgxLoader is now used to select datasets
directly from the Beacon response, rather than being used
internally.
- Modified pgxSegprocess to support usage with downloaded "pgxseg"
files from Progenetix.

                 Changes in version 1.1.6 (2024-08-05)                  

- Modified extract_general_results function to ensure it adapts
correctly to arrays.
- Moved callset and cnvstats data from the "g_variant" type to
"cnv_fraction" to better align with data types.
- Removed the pgxCount function and integrated its functionality into
pgxLoader with the "sample_count" type, streamlining such query.

                 Changes in version 1.1.5 (2024-07-30)                  

- Added config/datatable_mappings.yaml to define mapping rules
between
Beacon JSON responses and data tables.
- Modified metadata access to retrieve data directly from the Beacon
API instead of using the services/sampletable API.
- Enabled querying of analyses information.
- Updated the type parameter in pgxLoader to align more closely with
Beacon v2 model entities: biosamples, individuals, analyses, and
g_variants.
- Added entry_point parameter to pgxLoader.
- Removed filterLogic parameter from pgxLoader.
- Optimized parallel query for variants.
- Cleaned up code and vignettes.

                 Changes in version 1.1.3 (2024-06-14)                  

- Add pgxMetaplot function to generate survival plots from metadata
- Add num_cores parameter for parallel query of variants

                 Changes in version 1.1.2 (2024-05-03)                  

- Add segtoFreq function to allow CNV frequency calculation from
given
segment data

[phantasus](/packages/phantasus)
---------

                       Changes in version 1.25.5                        

- reworked annotation parsing

[phenomis](/packages/phenomis)
--------

                        Changes in version 1.7.8                        

- hypotesting: new parameter checks

                        Changes in version 1.7.6                        

- hypotesting: bug fix for MultiAssayExperiment objects

                        Changes in version 1.7.4                        

- writing: bug fix when writing MultiAssayExperiment with single
  colData column

                        Changes in version 1.7.2                        

- correcting: bug fix when removing features with 0 or NAs in all
  reference samples

[PhyloProfile](/packages/PhyloProfile)
------------

                       Changes in version 1.20.0                        

- added UMAP clustering

- improved domain plot #110

- fast mode for large profile plot

- option to change font

- option to show gene names

- option to sort genes by defined list

- option to upload preprocessed data from an input folder

[Pigengene](/packages/Pigengene)
---------

                 Changes in version 1.31.2 (2024-06-21)                 

General

- Links to project.eigen were added in related fucntions docs.

Changes in existing functions

- doMinimize was added to compute.pigengene.

[PIPETS](/packages/PIPETS)
------

                 Changes in version 1.1.5 (2024-09-27)                  

- Added strand specific analysis parameters. Can now use different
  threshAdjust and highOutlierTrim values for each strand should they
  choose.

                 Changes in version 1.1.4 (2024-09-11)                  

- Resolved warnings that popped up during Initial Poisson Distribution
  testing. Results will not be greatly changed but now PIPETS correctly
  counts all reads.

- Fixed problems in Genomic Ranges input method so that it now
  correctly mirrors the analysis performed for bed file inputs.

                 Changes in version 1.1.3 (2024-07-12)                  

- Changed read quality score parameter from a set value to a minimum,
  so anything equal to or greater than the input value will be
  considered

                 Changes in version 1.1.2 (2024-07-11)                  

- Changed PIPETS to use read quality score as a cutoff instead of read
  length. Functions as better trimming process

- Input bed files and GRanges objects now both correctly use read
  quality score

- Updated Vignette file wording to properly match recent changes

                 Changes in version 1.1.1 (2024-06-06)                  

- Resolved Git inconsistencies, made minor fixes

                 Changes in version 1.1.0 (2024-05-14)                  

- Added Benjamini-Hochberg multiple testing correction to help reduce
  incidence of Type I error

[plotgardener](/packages/plotgardener)
------------

                       Changes in version 1.11.2                        

NEW FEATURES

- plotgardener now supports .(m)cool files for Hi-C data! Functions
for reading and extracting features from .(m)cool files include:
- readCool
- readCoolBpResolutions
- readCoolNorms
- readCoolChroms .(m)cool files can now be used with any Hi-C
plotting function in place of .hic files or other data types.

                       Changes in version 1.11.1                        

BUG FIXES

- plotPairs and plotPairsArches check and swap order of anchors.
- plotTranscripts NULL label checking bug fix.

                       Changes in version 1.11.0                        

Version bump for Bioconductor 3.19 release.

[PLSDAbatch](/packages/PLSDAbatch)
----------

                        Changes in version 1.1.1                        

- Date: 2024-09-18
- Text: pca function clash
- Details: Avoid the clash of function pca() from vegan and mixOmics
packages

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.13.1 (2023-05-15)                 

- Made PoDCall compatible with multiplexed target assays
  (multiple target channels)

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

                       Changes in version 1.15.17                       

- Analyzing data with random effects in PomaLimma
- Select outcome factor in PomaBoxplots, PomaDensity, and
PomaOutliers
- Documentation improvements
- Introduces PomaEnrichment for gene set analysis
- Allows custom study designs in PomaDESeq

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.45                         

Changes in version 1.45.2

- pRolocmarkers() has a new version argument, to allow for new
markers
versions to be added.
- 14 new marker sets have been added to pRolocmarkers() under version
= 2 (new default).
- Documentation for pRolocmarkers() has been updated to include a
description of version = 2 markers.

Changes in version 1.45.1

- Import 'mclust::estep*()'.

Changes in version 1.45.0

- New devel version

[pRolocGUI](/packages/pRolocGUI)
---------

                        Changes in version 2.15                         

CHANGES IN VERSION 2.15.2

- Fix bug in compare app when selecting marker classes and plotting
the map
- Change organlle profile plots to 25% and 75% quantiles

CHANGES IN VERSION 2.15.1

- Remap in aggregation app.

CHANGES IN VERSION 2.15.0

- New version for Bioc devel

[ProtGenerics](/packages/ProtGenerics)
------------

                       Changes in version 1.37.1                        

- Add `estimatePrecursorIntensity()` generic.

[PSMatch](/packages/PSMatch)
-------

                         Changes in version 1.9                         

PSMatch 1.9.1

- Fix check errors.

PSMatch 1.9.0

- New Bioc devel.

[pwalign](/packages/pwalign)
-------

                        Changes in version 1.2.0                        

- No significant changes in this version.

[QDNAseq](/packages/QDNAseq)
-------

                 Changes in version 1.41.3 (2024-07-24)                 

Miscellaneous

- exportBins(..., format = "vcf") did not add meta-data information
for field END to the VCF file header.

                 Changes in version 1.41.2 (2024-07-23)                 

Bug Fixes

- poolRuns() would give an error "Error in colMeans2(oldphenodata,
cols = numericCols, useNames = FALSE) : Argument 'x' must be a
matrix or a vector."

                 Changes in version 1.41.1 (2024-07-20)                 

Miscellaneous

- Fix markup typo in help("exportBins").

- Add package test for poolRuns().

[QFeatures](/packages/QFeatures)
---------

                        Changes in version 1.15                         

QFeatures 1.15.2

- Fix typo in normalisation methods.

QFeatures 1.15.2

- Fix bug in QFeatures::longFormat() when rownames are numerical
(reported upstream
https://github.com/waldronlab/MultiAssayExperiment/issues/331).
- Starting New API unit test (see issue #214.

QFeatures 1.15.1

- Import reshape2::melt, required for
MultiAssayExperiment::longFormat().

QFeatures 1.15.0

- New Bioc devel

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

                       Changes in version 1.23.1                        

- Add coefCV to NAMESPACE

- Modify pcaPlot to allow user to omit the labelling

- Switch from gridExtra to patchwork in vignette

                       Changes in version 1.23.0                        

- Version bump for Bioconductor 3.20, no changes to package

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.30.0                        

New features

- metadata argument added to RaggedExperiment constructor function
(@LiNk-NY).

Bug fixes and minor improvements

- Ensure order is the same in both sparse and dense outputs of
compactAssay (@csoneson, #30)

[RAIDS](/packages/RAIDS)
-----

                        Changes in version 1.3.3                        

SIGNIFICANT USER-VISIBLE CHANGES

o More comprehensive vignette

                        Changes in version 1.3.2                        

SIGNIFICANT USER-VISIBLE CHANGES

o New functions inferAncestry(), inferAncestryGeneAware() and
getRefSuperPop() to simplify ancestry inference

                        Changes in version 1.3.1                        

SIGNIFICANT USER-VISIBLE CHANGES

o New function createAccuracyGraph() that creates a graphic
representation of the accuracy for different values of PCA dimensions
and K-neighbors through all tested ancestries.

[Rarr](/packages/Rarr)
----

                         Changes in version 1.5                         

- Fixed bug when creating an empty array with a floating datatype.
The
fill value would be interpreted as an integer by read_metadata() and
create and array of the wrong type.
- Fixed bug in update_zarr_array() when NULL was provided to one or
more dimensions in the index argument. This was parsed incorrectly
and the underlying zarr was not modified.
- Fixed bug in reading 64-bit integer arrays compressed with ZLIB or
LZ4.
The calculated decompression buffer size was too small and reading
would fail. (Thanks to Dan Auerbach for the report:
https://github.com/grimbough/Rarr/issues/10)
- Added support for the ZarrArray S4 class and the DelayedArray
framework.
- Improvements to read and write performance.

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

[Rcpi](/packages/Rcpi)
----

                 Changes in version 1.41.3 (2024-09-09)                 

Improvements

- Added early dependency checks for pwalign in functions that use
pairwise alignment. If certain dependency conditions can't be met,
the error is raised immediately, preventing unnecessary computation
(#26).
- Replaced rjson with jsonlite to maintain compatibility with R
< 4.4.0. The recent rjson update (v0.2.22 , 2024-08-20) requires R
(>= 4.4.0), which would break builds on r-oldrel (#28).

                 Changes in version 1.41.2 (2024-08-30)                 

Improvements

- Since Bioconductor 3.19 and Biostrings 2.72.0, the pairwise
sequence
alignment facilities have been moved from Biostrings into the
pwalign package. For maximum compatibility, we now detect the
installed Biostrings version at runtime and decide which package to
use for pairwise alignment, without introducing pwalign as an
additional hard dependency.

When calling calcTwoProtSeqSim() and calcParProtSeqSim(), if users
have Biostrings >= 2.72.0 installed while pwalign is not installed,
expect to see an explicit error in the results saying that pwalign
is required and should be installed from Bioconductor.

                 Changes in version 1.41.1 (2024-07-20)                 

Improvements

- Replaced RCurl with httr2 and curl for retrieving molecular and
sequence data from web APIs. Updated the outdated API endpoint URLs
for DrugBank and RCSB PDB (#18).
- inst/CITATION now uses bibentry() to replace citEntry() (#19).
- Fixed check notes on lost braces when running R CMD check under
R 4.4.x (#21).

[RCy3](/packages/RCy3)
----

                       Changes in version 2.26.0                        

- Added documentation & checks for super-long commands
  - suggest users to use commandsPost() when URI too long

- Bug fixes:
  - check for Inf and -Inf in double columns during loadTableData, #224
  - fix missing base.url in some functions

- Add test functions

[ReactomePA](/packages/ReactomePA)
----------

                       Changes in version 1.49.1                        

- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)

[recount](/packages/recount)
-------

                       Changes in version 1.31.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Remove warnings related to using coverage_matrix() or
expressed_regions() on Windows as rtracklayer::import() does work
with local BigWig files on that operating system. I'm not sure if it
will work with remote BigWig files given that remote BigWig file
access on other operating systems is not working due to
https://github.com/lawremi/rtracklayer/issues/83 and related issues.

[ResidualMatrix](/packages/ResidualMatrix)
--------------

                       Changes in version 1.16.0                        

- Minor change to the recommended way to call
  BiocSingular::runPCA(), to avoid block-processing when
  computing the center.

[Rhtslib](/packages/Rhtslib)
-------

                        Changes in version 3.2.0                        

- No significant changes in this version.

[ribosomeProfilingQC](/packages/ribosomeProfilingQC)
-------------------

                       Changes in version 1.17.1                        

- allow DNAStringSet for genome.

[ropls](/packages/ropls)
-----

                       Changes in version 1.37.6                        

- minor vignette update

                       Changes in version 1.37.4                        

- plot_score new name of the previous gg_scoreplot method

                       Changes in version 1.37.2                        

- gg_scoreplot: now handling the 1D or empty models

[rWikiPathways](/packages/rWikiPathways)
-------------

                       Changes in version 1.25.1                        

- bug fix: downloadPathwayArchive scrapes "File Name" column

[S4Arrays](/packages/S4Arrays)
--------

                        Changes in version 1.6.0                        

NEW FEATURES

- Add fast colsum() method for ordinary matrices.

- Define low-level generics:
  - subset_Array_by_logical_array()
  - subset_Array_by_Lindex()
  - subset_Array_by_Mindex()
  - subset_Array_by_Nindex()
  with default methods. Their purpose is to support subsetting (`[`)
  of Array derivatives and to make it easier for developers of Array
  extensions to implement subsetting for their objects.
  They are not intended to be used directly by the end user.

- Define low-level generics:
  - subassign_Array_by_logical_array()
  - subassign_Array_by_Lindex()
  - subassign_Array_by_Mindex()
  - subassign_Array_by_Nindex()
  with default methods. Their purpose is to support subassignment
  (`[<-`)
  of Array derivatives and to make it easier for developers of Array
  extensions to implement subassignment for their objects.
  They are not intended to be used directly by the end user.

SIGNIFICANT USER-VISIBLE CHANGES

- read_block() now uses SparseArray::extract_sparse_array() instead of
  DelayedArray::OLD_extract_sparse_array() behind the scene to extract
  and load sparse blocks from a sparse array-like object. As a
  consequence,
  sparse blocks are now returned as SparseArray objects (implemented in
  the SparseArray package) instead of SparseArraySeed objects
  (implemented
  in the DelayedArray package).

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.44.0                        

BUG FIXES

- Make sure that internal helper coerceToSimpleList() always propagates
  the mcols.

[SAIGEgds](/packages/SAIGEgds)
--------

                        Changes in version 2.4.1                        

- add an option `load.balancing` to `seqAssocGLMM_SPA()`

[SCArray](/packages/SCArray)
-------

                       Changes in version 1.14.0                        

- update according to DelayedArray(>= 0.31.5)

[SCArray.sat](/packages/SCArray.sat)
-----------

                        Changes in version 1.6.0                        

- update according to DelayedArray(>= 0.31.5)

[scBubbletree](/packages/scBubbletree)
------------

                       Changes in version 1.7.21                        

- new abstract

- new vignette

- new function: get_num_cell_tiles
  Changes after version 1.71.21 (17 Sep, 2024)

- scBubbletree paper accepted in BMC bioinformatics
  Changes after version 1.71.22 (20 Sep, 2024)

- simplified tree visualizations added in outout as object: tree_simple

                        Changes in version 1.7.2                        

- bubbletree comparison function get_bubbletree_comparison

- main functions simplified and defaults changed

- documentation updated

[scDblFinder](/packages/scDblFinder)
-----------

                 Changes in version 1.19.6 (2024-09-19)                 

- added a dbr.per1k parameter to set doublet rates per thousands of
  cells, updated the default from 1 to 0.8\%

- fixed some issues stemming from the cxds score in some corner cases
  (absence of inverse correlation between genes)

- updated documentation

[scMultiSim](/packages/scMultiSim)
----------

                        Changes in version 1.1.3                        

- Added the Shiny app to help users visualize the effect of each
parameter and adjust the simulation options.
- Added the speed.up parameter to enable experimental speed
optimization.
- Bug fixes and improvements.

                        Changes in version 1.1.0                        

Prepare for the Bioconductor release

- Fix build errors

[scp](/packages/scp)
---

                        Changes in version 1.15                         

scp 1.15.2

- fix: fixed x-axis direction annotation for volcano plot on contrast
- fix: solved bug in DA when missing contrast level in modelled
feature (issue #65).
- Add link to Leduc SCP.replication vignette.

scp 1.15.1

- test: added unit tests for scplainer: ScpModel-Class,
scpModelFit-Class, ScpModel-Workflow

scp 1.15.0

- New Bioconductor 3.20 (devel) release

[scran](/packages/scran)
-----

                        Changes in version 1.34                         

- Bugfix for quickCluster() to pass along arguments to the
  internal per-block call.

[scRNAseqApp](/packages/scRNAseqApp)
-----------

                       Changes in version 1.15.22                       

- add warning message to log file when the atac bigwig files are not
  available.

                       Changes in version 1.15.21                       

- fix a bug that counter table is not available.

                       Changes in version 1.15.20                       

- Add csv format for figure data download.

                       Changes in version 1.15.19                       

- Try to handle the error `Timeout was reached: [www.ncbi.nlm.nih.gov]
  SSL/TLS connection timeout`.

                       Changes in version 1.15.18                       

- add protection of column sampleID in cell info editor funcitons.

                       Changes in version 1.15.17                       

- Error handle if metadata contains '|'.

- add cell info editor funcitons.

                       Changes in version 1.15.16                       

- add more download format.

                       Changes in version 1.15.15                       

- fix the low resolution of the plot.

                       Changes in version 1.15.14                       

- fix the dotplot download issue.

                       Changes in version 1.15.13                       

- add detailed fontsize and font family control.

- change the dotplot from ggplot2 to complexheatmap.

                       Changes in version 1.15.12                       

- order the stacked violin plots by input order.

                       Changes in version 1.15.11                       

- add reorder function for splitBy for violin plots.

                       Changes in version 1.15.10                       

- add color picker.

- add splitBy for violin plots.

- fix the issue that refresh explorer will lost the subset info.

                       Changes in version 1.15.9                        

- add help video playlist.

                       Changes in version 1.15.8                        

- add functionality to allow different labels for cell info.

- allow subset groups for explorer.

                       Changes in version 1.15.7                        

- fix the issue of reload everything when edit the cell info.

                       Changes in version 1.15.6                        

- Editable cell info during login mode.

- Hiden the cells when check the hidden box.

                       Changes in version 1.15.5                        

- Response to the download width and height.

                       Changes in version 1.15.4                        

- Access tab via query string.

                       Changes in version 1.15.3                        

- Access data via data title.

                       Changes in version 1.15.2                        

- Add favicon.

                       Changes in version 1.15.1                        

- Fix the unique vistor counter.

[scuttle](/packages/scuttle)
-------

                        Changes in version 1.16                         

- addPerCellQCMetrics() will (optionally) add the feature subset
  identities to the row data of the SingleCellExperiment.

- Added uniquifyDataFrameByGroup() to collapse grouped rows in a
  DataFrame into their unique values.

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.46.0                        

UTILITIES

- `seqGetData()` return NULL, if 'var.name=character()'

                       Changes in version 1.44.3                        

UTILITIES

- update the C codes according to '_R_USE_STRICT_R_HEADERS_=true' &
  '_R_CXX_USE_NO_REMAP_=true'

                       Changes in version 1.44.2                        

BUG FIXES

- fix `seqAddValue(, val=vector("list", NUM_VARIANT))`

- fix the ploidy returned from `seqVCF_Header()`, when there are
  genotypes
  of males and females on Chromosome X

                       Changes in version 1.44.1                        

UTILITIES

- new option 'numvariant' in `seqEmptyFile()`

BUG FIXES

- `seqMerge()` should internally use "chr_position_ref_alt" to
  distinguish
  the variants in different files

- `seqAddValue(, varnm="annotation/filter")` should work with a factor
  variable

- `seqAddValue(, varnm="variant.id")` can reset the variant IDs with a
  different number of the variants

[ShortRead](/packages/ShortRead)
---------

                        Changes in version 1.64                         

BUG FIXES

- (v 163.1) rely on system-provided zlib on all platforms

[simplifyEnrichment](/packages/simplifyEnrichment)
------------------

                        Changes in version 2.0.0                        

- use **simona** for semantic similarity calculation

- removed similarity functions for DO and based on overlaps (the count)

[SingleCellAlleleExperiment](/packages/SingleCellAlleleExperiment)
--------------------------

                        Changes in version 1.0.2                        

Fixed error in tests

                        Changes in version 1.0.1                        

Updated allele names in the vignette after data package "scaeData"
was
updated

[SingleR](/packages/SingleR)
-------

                        Changes in version 2.8.0                        

- Added the test.genes= argument to trainSingleR(), to restrict
  marker detection to only those genes in the test dataset. This
  is also checked against rownames(test) in classifySingleR() to
  ensure that the test's feature space is consistent with the
  space used during training.

- The introduction of test.genes= means that we no longer need to
  explicitly subset the rows of the reference dataset (to match
  the test features) in SingleR(). This saves memory by avoiding
  an unnecessary copy of the reference dataset, but may also
  slightly alter the marker selection as ties are broken in a
  different way. Namely, if the top X genes are used as markers,
  and the X-th and (X+1)-th gene have the same log-fold change,
  tie breaking will be based on the ordering of the rows in the
  reference matrix - which is no longer the same as in the
  previous version of SingleR. This results in some slight
  differences in the markers that propagate down to the
  classification results.

- Restored the BNPARAM= argument in trainSingleR(), to enable
  more fine-grained specification of neighbor search algorithms.
  The approximate= argument is deprecated.

- Soft-deprecated check.missing= in classifySingleR() and
  combineRecomputedResults(). This is because any filtering will
  cause a mismatch between the row names of tests and the
  test.genes in trained. Rather, filtering should be done prior
  to trainSingleR(), as is done in the main SingleR() function.

- combineRecomputedResults() now supports fine-tuning to resolve
  closely-related labels from different references. This is
  similar to the fine-tuning in classifySingleR() where the
  feature space is iterately redefined as the union of markers of
  labels with near-highest scores.

- Added the plotMarkerHeatmap() function to plot a diagnostic
  heatmap of the most interesting markers for each label.

[sketchR](/packages/sketchR)
-------

                        Changes in version 1.1.3                        

- Update environment to avoid setuptools incompatibility

                        Changes in version 1.1.2                        

- Specify Linux aarch environment

                        Changes in version 1.1.1                        

- Update environments to adapt to basilisk conda update

[smartid](/packages/smartid)
-------

                        Changes in version 1.1.2                        

- Update marker selection functions to fix wrong names of marker
list.

                        Changes in version 1.1.1                        

- Update cal_score() function to convert input sparse matrix into
dense matrix.

[SNPRelate](/packages/SNPRelate)
---------

                       Changes in version 1.38.1                        

UTILITIES

- update the C codes according to '_R_USE_STRICT_R_HEADERS_=true' &
  '_R_CXX_USE_NO_REMAP_=true'

[SpaceMarkers](/packages/SpaceMarkers)
------------

                        Changes in version 1.1.3                        

- Optimized the long running row.dunn.test() function
- Corrected sparse -> dense conversions

                        Changes in version 1.1.2                        

- Added getPairwiseInteractingGenes which enables pairwise analysis
of
interacting patterns

                        Changes in version 1.1.1                        

- getSpatialFeatures: add default method to infer the object passed
to
it.

[SparseArray](/packages/SparseArray)
-----------

                        Changes in version 1.6.0                        

NEW FEATURES

- Linear subsetting of a SVT_SparseArray object e.g. svt&#91;11:13&#93;.

- Add the 'dim' and 'dimnames' arguments to SVT_SparseArray()
  constructor
  function. Main use case is to make it easy and efficient to construct
  a zero-only SVT_SparseArray object of arbitrary dimensions.

- Add is_nonzero() generic with a default method and fast methods for
  SVT_SparseArray, COO_SparseArray, and sparseMatrix objects.

- The Sparse Vector Tree of a SVT_SparseArray can have "lacunar
  leaves",
  that is, leaves where the nzvals are missing. These leaves are
  semantically equivalent to leaves where all the values in nzvals are
  ones. This reduces the memory footprint of a "logical"
  SVT_SparseArray
  with no NAs by half.
  Note that a "logical" SVT_SparseMatrix object with no NAs is similar
  to an ngCMatrix object from the Matrix package. ngCMatrix objects
  don't store any non-zero values either (they're implicitly considered
  to be TRUEs), only their offsets.

- col/row summarization methods (a.k.a. "matrixStats methods") now work
  on COO_SparseArray objects (in addition to SVT_SparseArray objects).

- Add row/colSums2() and row/colMeans2() methods for SparseArray
  objects.

- rowsum(), crossprod(), tcrossprod(), and %*% methods now work on
  COO_SparseMatrix objects (in addition to SVT_SparseMatrix objects).

- Add fast colsum() methods for SparseMatrix and dgCMatrix objects.

- Add is.na(), is.nan(), and is.infinite() methods for SVT_SparseArray
  objects.

- Add pmin() and pmax() methods for SparseArray objects.

- aperm(&#60;SVT_SparseArray&#62;) now supports S4Arrays::aperm2() extended
  semantic. See '?S4Arrays::aperm2' for details.

- Add 'dimnames' argument to randomSparseArray(), poissonSparseArray(),
  randomSparseMatrix(), and poissonSparseMatrix().

- Coercions from dgTMatrix, lgTMatrix, or ngTMatrix, to
  COO_SparseMatrix/Array or SVT_SparseMatrix/Array, and vice versa.

- Coercions from ng[C|R|T]Matrix to COO_SparseMatrix, and vice versa.

- Coercion from ngCMatrix to SVT_SparseMatrix, and vice versa.

- NaArray objects: WORK-IN-PROGRESS! New objects that use the same
  internal representation as SVT_SparseArray objects (Sparse Vector
  Tree), but background value is NA instead of zero. See
  https://github.com/fmicompbio/footprintR/issues/7 for the
  motivating use case. They support most operations supported by
  SVT_SparseArray objects: [, [<-, t(), aperm(), cbind(), rbind(),
  abind(), matrixStats operations (col*(), row*()), etc...
  What's missing:
  - a dedicated vignette;
  - some row*() functions are not ready yet.

- Add new generics is_nonna(), nnacount(), nnawhich(), nnavals(), and
  `nnavals<-`(), with default methods and methods for NaArray objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Refactored subassignment of a SVT_SparseArray object, which resulted
  in significant speed improvement and memory footprint reduction.

- Speed up some row summarization methods (e.g. rowSums(x)) by
  implementing
  them natively in C rather than doing a transposition followed by a
  column
  summarization (e.g. colSums(t(x))). This avoids the costly
  transposition
  step which is very time and memory consuming.
  Row summarization methods with native C implementation so far:
  rowAnyNAs() (35x speedup), rowMins() (6x speedup), rowMaxs() (6x
  speedup),
  rowSums() (20x speedup), rowMeans() (18x speedup), rowVars() (5x
  speedup),
  and rowSds() (5x speedup).
  More methods will follow.

- Special-case 'lambda=0' in poissonSparseArray() so the empty
  SVT_SparseArray object gets returned instantaneously.

BUG FIXES

- Fix integer overflow in nzwhich() methods for CsparseMatrix and
  RsparseMatrix objects when the object has a length >= 2^31.

- Fix bug in coercion from COO_SparseMatrix to &#91;d|l&#93;gCMatrix or
  &#91;d|l&#93;gRMatrix when the 'nzcoo' slot of the COO_SparseMatrix object
  contains duplicates.

- Make sure that coercion from CsparseMatrix to SVT_SparseMatrix works
  on any &#91;d|l|n&#93;gCMatrix **derivative** and not just on a
  &#91;d|l|n&#93;gCMatrix
  **instance**.

[sparseMatrixStats](/packages/sparseMatrixStats)
-----------------

                        Changes in version 1.17                         

- Fix handling of missing values in `rowSums2` if `cols` is a boolean
  vector

- Implement optimized code path for `rowSums2` and `rowMeans2` if
  `cols` is provided

[SpatialFeatureExperiment](/packages/SpatialFeatureExperiment)
------------------------

                        Changes in version 1.7.1                        

- Added image setter, Img<-
- Implemented spatial aggregation functions to aggregate directly
from
transcript spot file, from rowGeometry, or from cell geometries in
SFE objects
- Implemented splitByCol to split SFE objects by geometry,
splitSamples to split by sample_id, and splitContiguity to split by
cotiguity of an annotGeometry

                        Changes in version 1.6.1                        

- readRDS converts old style SpatRasterImage to the new style
- readSelectTx and addSelectTx functions to read transcript spots
from
a few select genes from the parquet output of formatTxSpots or add
them to an SFE object
- Added formatTxTech and addTxTech functions, basically thin wrappers
of formatTxSpots and addTxSpots with presets for Vizgen, Xenium, and
CosMX

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 2.11.4 (2024-07-28)                 

- Reduced dependencies.

                 Changes in version 2.11.2 (2024-07-23)                 

- The function cvt_id accepts vectors.

                 Changes in version 2.11.1 (2024-06-14)                 

- Fixed some contrasts in spatial enrichment.

[Spectra](/packages/Spectra)
-------

                        Changes in version 1.15                         

Changes in 1.15.13

- Add precursorMz<- method issue #336.

Changes in 1.15.12

- Add generic backendRequiredSpectraVariables() to allow definition
of
mandatory spectra variables for a backend.

Changes in 1.15.11

- Add reference to MsBackendMetaboLights.

Changes in 1.15.10

- Add new extractSpectra() generic and implementation for MsBackend.
Fixes issue #5.

Changes in 1.15.9

- Restructure and reorganize documentation for Spectra.

Changes in 1.15.8

- Refactor the Spectra() constructor method: better support for
initialization of backends that define their own specific
parameters.

Changes in 1.15.7

- Change estimatePrecursorIntensity() to a method to avoid
overrides/clashes with the same-named implementation in xcms.

Changes in 1.15.6

- Fix in selectSpectraVariables() for MsBackendMzR: ensure peaks
variables "mz" and "intensity" are not by default removed.

Changes in 1.15.5

- Add new filterPeaksRanges() function to filter mass peaks by ranges
on numeric spectra or peak variables.

Changes in 1.15.3

- For evaluation of the Spectra's processing queue: call functions
from the MetaboCoreUtils directly through their namespace
(MsCoreUtils::) to avoid errors if performed in parallel on Windows
machines or if called on a re-loaded object.
- New asDataFrame() function to convert a (small) Spectra object into
a long DataFrame.

Changes in 1.15.2

- Add dataStorageDataPath() and dataStorageDataPath<- methods to
allow
updating/adapting the path of the data storage files of backends
supporting that issue #321.

Changes in 1.15.1

- Improve documentation for combineSpectra() and combinePeaks() issue
#320.

[SPIAT](/packages/SPIAT)
-----

                        Changes in version 1.7.2                        

BUG FIXES

- Fixed cluster assignment issue in function
identify_neighborhoods().

                        Changes in version 1.7.1                        

BUG FIXES

- Added feature_colname parameter to functions
plot_cell_marker_levels() and marker_intensity_boxplot().

                        Changes in version 1.7.0                        

Development version on Bioconductor 3.20.

[splatter](/packages/splatter)
--------

                 Changes in version 1.30.0 (2024-10-30)                 

- 
  Minor maintenance chores:
  
  • Update roxygen2 version
  
  • Remove old CI config files
  
  • Add PR commands GitHub action
  
  • Re-configure pkgdown site
  
  • Add dependabot GitHub action
  
  • Update README

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

                       Changes in version 1.35.1                        

- The long format for output from shiftCor

- The long format for input from statAnalysis

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.36.0                        

NEW FEATURES

- Calling saveRDS() on a SummarizedExperiment object that contains
  out-of-memory data now raises an error with a message that redirects
  the user to HDF5Array::saveHDF5SummarizedExperiment().

SIGNIFICANT USER-VISIBLE CHANGES

- Move documentation of constructor function SummarizedExperiment()
  from RangedSummarizedExperiment-class.Rd to
  SummarizedExperiment-class.Rd.
  See https://github.com/Bioconductor/SummarizedExperiment/issues/80

- Change default value of 'rowData' argument in SummarizedExperiment()
  constructor from GRangesList() to NULL.

BUG FIXES

- Fix typos in error message from assay() and assays() setters.

[SVMDO](/packages/SVMDO)
-----

                 Changes in version 1.5.5 (2024-08-30)                  

- Fixing enrichDO-related error
- Fixing deprecation issues about data table visualization in Shiny

[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.17.7                        

- ExoLabel is much much faster and does a better job cleaning up when
aborted early
- ExoLabel now has fewer arguments
- Updates to man pages

                       Changes in version 1.17.6                        

- Lots of internal improvements to ExoLabel to increase computational
speed and decrease disk usage.
- ExoLabel will no longer crash if given relative paths.
- Adds more internal error checking to prevent some rare bugs.
- Updates man pages to reflect new changes.
- Updates EstimateExoLabel to reflect new changes.

                       Changes in version 1.17.5                        

- ExoLabel will no longer brick R during sorts on large files.
- ExoLabel reports more progress during some lengthy processing
sections when verbose=TRUE
- Known issue: "Copying source file" step is still non-interruptable,
will be fixed in a later update

                       Changes in version 1.17.4                        

- ExoLabel now allows an inflation argument to control application of
inflation

                       Changes in version 1.17.3                        

- predict.EvoWeaver now supports returning p-values separately from
raw score for some algorithms.
- OpenMP implementation for EvoWeaver algorithms that support it has
been fixed

                       Changes in version 1.17.2                        

- RandForest function added to train random forest models
- Associated man pages for RandForest and DecisionTree objects
- New methods for DecisionTree objects to plot and coerce to
dendrogram
- Small bugfix to subset.dendrogram

                       Changes in version 1.17.1                        

- Major updates to EvoWeaver:
- predict.EvoWeaver now returns a data.frame by default
- Method arguments are updated to match their names in the
associated EvoWeaver manuscript
- Above changes have propagated to documentation files
- New Phylogenetic Profiling methods with improved accuracy
- New meta-methods PhylogeneticProfiling, PhylogeneticStructure,
GeneOrganization, SequenceLevel for predict.EvoWeaver
- New pre-trained Ensemble models have been included
- Updates to ExoLabel for better status printing

[TargetSearch](/packages/TargetSearch)
------------

                        Changes in version 2.8.0                        

NEW FEATURES

- Allow incompatible data when combining objects `tsLib` and
  `tsSample`.
  The idea is if the column names of the `data` slots are different,
  the `c` operator do not fail. There are, however, limitations with
  some
  data types (like lists or matrices). If this is the case, then throw
  a more meaningful error message. (commit 6fb2f59)

BUG FIXES

- Fix `FAMEoutliers` manual grouping: The comparison should use
  character
  vectors instead of numeric. Also, an error was thrown due to an `if`
  condition not having an scalar return (due to using `is.na`).
  Note that since this version, the default values of startDay and
  endDay
  have been set as NULL. (commit 69e8a2f).

INTERNAL

- Add tests for `FAMEoutliers`.

                        Changes in version 2.6.2                        

BUG FIXES

- Use \providecommand for missing defines in vignettes instead

                        Changes in version 2.6.1                        

BUG FIXES

- Fix vignettes building by adding missing latex commands due to
  an incompatibility between BiocStyle and knitr

[TBSignatureProfiler](/packages/TBSignatureProfiler)
-------------------

                       Changes in version 1.17.0                        

Minor Changes

- Added Li_3 signature
- Adjusted function to create documentation in Signature Addition
vignette
- Fixed Zhao_Nano_6 signature documentation

Bug Fixes

- Changed out .data pronouns to use proper tidyr variable calling

[TENxIO](/packages/TENxIO)
------

                        Changes in version 1.8.0                        

New features

- The import method for h5 files now includes all rowData which
typically includes "ID", "Symbol", and "Type" columns.

Bug fixes and minor improvements

- rownames are set to first column in the features.tsv.gz data
(rownames for h5 files are determined by HDF5Array::TENxMatrix)
- TENxH5 now tests whether there is a /matrix/features/interval
dataset in the h5 file. It sets ranges to NA_character_ when
interval data is not found.

[tidySpatialExperiment](/packages/tidySpatialExperiment)
---------------------

                        Changes in version 1.0.1                        

- Greatly improved interactive gating, facilitated via tidygate.

[TileDBArray](/packages/TileDBArray)
-----------

                       Changes in version 1.16.0                        

- Minor fix for as.data.frame= deprecation in tiledb_array().

- Support other datatypes for the dimensions and storage when
  configuring a TileDBRealizationSink. This is achieved via the
  new storagetype= and dimtype= arguments. Also added
  getTileDBDimType() and setTileDBDimType() to globally define
  the choice of dimension datatype.

- Bugfix to TileDBArraySeed to correctly handle dimension domains
  that do not start at 1. This requires a modification to the
  class to record the domain offset.

- Added a offset= option to TileDBRealizationSink() to create
  arrays with dimension domains that do not start at 1. This
  requires a modification to the associated class to record the
  domain offset.

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

                       Changes in version 1.41.10                       

- Improve the plot speed for interaction data.

                       Changes in version 1.41.9                        

- Add a simple gene track plot model.

                       Changes in version 1.41.8                        

- Add listNormalizations function.

                       Changes in version 1.41.7                        

- Fix the bug that newpage parameter not used for plotGRanges.

                       Changes in version 1.41.6                        

- Add `legendPosition` parameter for lollipop plot.

                       Changes in version 1.41.5                        

- Export the method `$` for trackStyle class.

                       Changes in version 1.41.4                        

- Move the 3d plots to a different package.

                       Changes in version 1.41.3                        

- plot lolliplot with multiple shapes in one position.

                       Changes in version 1.41.2                        

- export view3dStructure function.

                       Changes in version 1.41.1                        

- Add mdsPlot function.

[transmogR](/packages/transmogR)
---------

                        Changes in version 1.1.1                        

- Added digestSalmon()

[transomics2cytoscape](/packages/transomics2cytoscape)
--------------------

                       Changes in version 1.15.2                        

DATA or DOCUMENT CHANGES

- Update Version Information in vignette

                       Changes in version 1.15.1                        

DATA or DOCUMENT CHANGES

- Add inst/CITATION

[treeclimbR](/packages/treeclimbR)
----------

                        Changes in version 1.1.1                        

- Adapt unit tests to limma updates

[treeio](/packages/treeio)
------

                       Changes in version 1.29.2                        

- speedup read.beast() with multithreading supports (2024-10-27, Sun,
#128)

                       Changes in version 1.29.1                        

- use yulab.utils::yulab_msg() for startup message (2024-07-26, Fri)
- support treetime output (2024-07-18, Thu)
- https://treetime.readthedocs.io/en/latest/tutorials/clock.html

[TVTB](/packages/TVTB)
----

                 Changes in version 1.31.2 (2024-07-05)                 

Bug fix

- Fix tSVE() shiny application.

                 Changes in version 1.31.1 (2024-07-04)                 

Bug fix

- Copy parseCSQToGRanges from ensemblVEP

[txdbmaker](/packages/txdbmaker)
---------

                        Changes in version 1.2.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Adjust makeTxDbFromBiomart() examples to reflect breaking change
  in Ensembl 113. See commit a5a885d for the details.

BUG FIXES

- Small fix to internal helper .Ensembl_getMySQLCoreDir(). See commit
  aa586a6 for the details.

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.23.5                        

- GENCODE 47 (H.s.), M36 (M.m), and Ensembl 113

                       Changes in version 1.23.1                        

- GENCODE 46 (H.s.), M35 (M.m), and Ensembl 112

[UCSC.utils](/packages/UCSC.utils)
----------

                        Changes in version 1.2.0                        

- No significant changes in this version.

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.22.3                        

- Fixed make_DBscores() crash.

                       Changes in version 1.22.2                        

- Added CITATION file.

                       Changes in version 1.22.1                        

BUG FIXES

- Always check object sizes before indexing in universalmotif_cpp().
  Fixes build timeouts on Bioconductor.

[updateObject](/packages/updateObject)
------------

                       Changes in version 1.10.0                        

- No changes in this version.

[variancePartition](/packages/variancePartition)
-----------------
                       Changes in version 1.35.4                        

- June 13, 2024
- fix bug in [.MArrayLM2 dropping residuals to a vector

                       Changes in version 1.35.3                        

- June 11, 2024
- in mvTest() with shrink.cov = TRUE uses lambda = 0.01

                       Changes in version 1.35.2                        

- June 6, 2024
- in mvTest() return beta and se

                       Changes in version 1.35.1                        

- May 31, 2024
- version bump for BioC

[velociraptor](/packages/velociraptor)
------------

                       Changes in version 1.15.13                       

- Add separate environment (untested) for Linux Aarch64. Cannot break
that platform more than it already is.

                       Changes in version 1.15.12                       

- Update Conda environment for Linux to avoid Anaconda channel.

                       Changes in version 1.15.11                       

- Re-enable plotVelocityStream() with a warning.

                       Changes in version 1.15.10                       

- Temporarily disable plotVelocityStream() due to unexplained issue
related to metr::geom_streamline()

                       Changes in version 1.15.9                        

- Update Conda environment for Linux and MacOSX Arm.
- Patch GitHub Action to use GitHub version of remotes.

                       Changes in version 1.15.8                        

- Update Conda environment using micromamba for MacOSX Arm.
- Fix switch between MacOSX and MacOSX Arm environments.

                       Changes in version 1.15.7                        

- Update Conda environment using micromamba for Windows.

                       Changes in version 1.15.6                        

- Update Conda environment to use anaconda channel on Linux. Passed
GitHub Action
https://github.com/kevinrue/velociraptor/actions/runs/10612115572/job/29413105915.

                       Changes in version 1.15.5                        

- Update Conda environment to svelo==0.3.2 on Linux.

                       Changes in version 1.15.4                        

- Set scvelo version triggering deprecation error to 0.3.1.

                       Changes in version 1.15.3                        

- Revert environment for Linux to the one of Bioconductor
release 3.18.

                       Changes in version 1.15.2                        

- Add environment for macOS (Intel); same environment as macos (M1).

                       Changes in version 1.15.1                        

- Fix issue #63.
- Update scvelo to 0.3.2 (conda-forge) for macOS (M1) and Linux.
- Update scvelo to 0.2.5 (bioconda) for Windows.
- Add mechanism to switch Conda environment (and scvelo version)
based
on operating system and architecture.
- Use scanpy.pp.neighbors to calculate neighbors due to deprecation
of
automatic neighbor calculation in scvelo.pp.moments.
- Update vignette to document the change of default value for
n_neighbors from scvelo (30) to scanpy (15).

[VisiumIO](/packages/VisiumIO)
--------

                        Changes in version 1.2.0                        

New features

- Support VisiumHD file formats including parquet.
- Include format argument for h5 file imports (default remains mtx).
- Add example data for Visium and VisiumHD imports.

[Voyager](/packages/Voyager)
-------

                        Changes in version 1.7.2                        

- Fixed bug that caused error when computing Lee's L on DelayedArray

                        Changes in version 1.7.1                        

- Revamped user interface of plotGeometry to allow plotting multiple
col, row, and annot geometries at once
- Allow plotting transcript spots in plotGeometry,
plotSpatialFeature,
plotLocalResult, and spatialReducedDim

[xcms](/packages/xcms)
----

                         Changes in version 4.3                         

Changes in version 4.3.3

- Fix issue #755: chromatogram() with msLevel = 2 fails to extract
chromatographic data if isolationWindowTargetMz is not specified or
available (e.g. for MSe data).
- Support coercing from XcmsExperiment to XCMSnExp with as(object,
"XCMSnExp").
- Change estimatePrecursorIntensity() to a method and add an
implementation for MsExperiment objects.

Changes in version 4.3.2

- Remove data/results import/export functionality as it is being
developed in the MsIO package.

Changes in version 4.3.1

- Support excluding samples or sample groups from defining features
with PeakDensity correspondence analysis (issue #742).
- Add plotPrecursorIons() function.
- Fix in dropFeatureDefinitions() that was not correctly removing
additional metadata from gap-filled chromatographic peaks.


NEWS from existing Data Experiment Packages
===================================

[CardinalWorkflows](/packages/CardinalWorkflows)
-----------------

                       Changes in version 1.37.1                        

BUG FIXES

- Updates for Cardinal v3.8

[curatedPCaData](/packages/curatedPCaData)
--------------

                 Changes in version 1.1.1 (2024-05-28)                  

- Minor fixes for Imports/Suggests for vignettes

[gDRtestData](/packages/gDRtestData)
-----------

              Changes in version 22023-05-15 (2023-05-15)               

- fix related with data.table

               Changes in version 2024-07-09 (2024-07-09)               

- reprocess datasets as per new columns added into Metrics assay

[geneLenDataBase](/packages/geneLenDataBase)
---------------

                       Changes in version 1.42.0                        

- Role of maintainer taken over by Federico Marini

- Refreshing the codebase, introducing roxygen-based documentation and
  applying styler for consistent spacing, as well as using Github
  Actions for CI/CD

[JohnsonKinaseData](/packages/JohnsonKinaseData)
-----------------

                        Changes in version 1.1.1                        

- Version bump

                        Changes in version 1.1.0                        

- Added tyrosine kinase PWMs published in Yaron-Barir et al. 2024

- Added option to mach acceptor specificity

- Added upstream/downstream options to flexibly trim phospho-peptides

[MetaScope](/packages/MetaScope)
---------

                        Changes in version 3.19                         

Bug Fixes

- Identified the Rbowtie2 parameter k as being doubled when specified
  for bowtie filter or align steps.

- Added another call to taxize in convert_animalcules() to catch any
  accessions that were not mapped to a UID in metascope_id(), in
  addition to another call in metascope_id() itself

Major changes

- Altered bt2_params objects to reflect 98% identity (16S), 95%
  identity (metagenomics) and added a parameter for when the origin
  genome is thought to not be present in the reference database.

- Added SILVA species_headers object called internally in
  convert_animalcules_silva

- Added convert_animalcules_silva function

                        Changes in version 3.18                         

Bug Fixes

- Fixed taxonomy table function to output correctly formatted table

- Fixed examples for various functions that were calling genomes with
  download_refseq but genomes were not able to be found.

- Fixed plot generation for metascope_id()

- Fixed premature stopping of download_refseq for strains labeled as
  "no rank" in NCBI.

- Fixed identification of reads as unknown genomes (due to outdated
  reference databases that identify genomes now removed from NCBI).
  The unknown genomes will now be distinctly identified based on NCBI
  accession ID, which will separate them in the final results. IDs can
  also be looked up manually (on NCBI website) to see what the reads
  were aligning to.

- Incorporated unknown taxa (removed from NCBI databases) into output
  of convert_animalcules()

Major changes

- Added ability to identify reads from databases other than NCBI
  (known to work for Silva and Greengenes2)

[MicrobiomeBenchmarkData](/packages/MicrobiomeBenchmarkData)
-----------------------

                 Changes in version 1.7.1 (2021-09-26)                  

- Added a function (scml) for re-calibrating the
  "Stammler_2016_16S_spikein" dataset with SCML: spike-in-based
  calibration to total microbial load.

[microbiomeDataSets](/packages/microbiomeDataSets)
------------------

                 Changes in version 1.13.2 (2024-06-26)                 

- LahtiWAData removed

- HintikkAXOData removed

[pRolocdata](/packages/pRolocdata)
----------

                       Changes in version 1.43.1                        

- move content of inst/extdata to the pRolocdata dataverse
  (https://dataverse.uclouvain.be/dataverse/pRolocdata).

                       Changes in version 1.43.0                        

- new devel version

[scaeData](/packages/scaeData)
--------

                        Changes in version 1.1.2                        

- Lookup tables have been updated to reflect complete allele names

                        Changes in version 1.1.1                        

- Changed filter_mode parameter value in the vignette to avoid check
  error

[scATAC.Explorer](/packages/scATAC.Explorer)
---------------

                       Changes in version 1.12.0                        

- corrected a bug when searching in-between years

- added h5ad support to saveATAC()

- added more examples in the help documentation for queryATAC()

[scMultiome](/packages/scMultiome)
----------

                        Changes in version 1.4.2                        

- replace writeSparseMatrix with alabaster.matrix::writeSparseMatrix

[scRNAseq](/packages/scRNAseq)
--------

                       Changes in version 2.20.0                        

- Support complex queries in searchDatasets() via a
  human-friendly syntax.

[SpatialDatasets](/packages/SpatialDatasets)
---------------

                 Changes in version 1.5.0 (2022-04-27)                  

- reformat datasets to SpatialExperiment version 1.5.3

                 Changes in version 1.3.3 (2022-01-31)                  

- add new datasets ST_mouseOB, SlideSeqV2_mouseHPC

- reformat datasets to SpatialExperiment version 1.5.2

[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.17.6                        

BUG FIXES

- Fixed the bug reported by @lahuuki about vis_grid_clus() not
  handling logical() cluster variables. See
  https://github.com/LieberInstitute/spatialLIBD/issues/80. To resolve
  this, sort_clusters() and get_colors() had to change internally.
  Examples and documentation for both functions have now been updated
  to showcase what happens when you provide a logical() vector as an
  input.

                       Changes in version 1.17.5                        

NEW FEATURES

- Added add_qc_metrics() inspired by
  https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/code/07_spot_qc/01_qc_metrics_and_segmentation.R
  which adds seven new columns to the colData(spe) that can be useful
  when performing quality control of the data. Developed by @lahuuki.

                       Changes in version 1.17.3                        

NEW FEATURES

- Added support for SpatialExperiment objects created with
  visiumStitched::build_spe()
  https://research.libd.org/visiumStitched/reference/build_spe.html
  that stitch together multiple Visium capture areas. Developed by
  @Nick-Eagles.


NEWS from existing Workflows
===================================

[recountWorkflow](/packages/recountWorkflow)
---------------

                       Changes in version 1.29.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Implement a workaround to
https://github.com/lawremi/rtracklayer/issues/83 which currently is
limiting the ability to remotely access BigWig files using
rtracklayer::import(). This affects the functions
recount::expressed_regions(), recount::coverage_matrix(), and
derfinder::getRegionCoverage() used in recountWorkflow.

Deprecated and Defunct Packages
===============================

**SOFTWARE:**

Sixty five software packages were removed from this release (after being deprecated
in Bioc 3.19):

- BDMMAcorrect, beadarraySNP, BHC, biodbLipidmaps, BioNetStat, CancerInSilico, CancerSubtypes, cellHTS2, CNVgears, compartmap, contiBAIT, CoRegNet, CORREP, crisprseekplus, dpeak, EBSeqHMM, eegc, enrichTF, ensemblVEP, exomePeak2, farms, FCBF, flowMap, FoldGO, FScanR, FunChIP, GOSim, HumanTranscriptomeCompendium, ImmuneSpaceR, InterMineR, IntOMICS, IRISFGM, iterClust, maigesPack, metagene, MetaVolcanoR, miRmine, MMAPPR2, MobilityTransformR, multiOmicsViz, NeighborNet, oneSENSE, openPrimeRui, pathVar, pcxn, PERFect, phemd, PloGO2, proteasy, PSEA, pwOmics, RefPlus, ReQON, restfulSE, RIPAT, RLSeq, SimBindProfiles, SMAP, sparseDOSSA, SpidermiR, SQUADD, StarBioTrek, STROMA4, TimiRGeN, TNBC.CMS

- Please note:  DNABarcodes, cliqueMS, CoSIA, NetPathMiner, and TnT, previously announced as
deprecated in 3.19, have been updated and remain in Bioconductor. 

Twenty One software packages are deprecated in this release and will be removed in Bioc 3.21:

- ATACCoGAPS, BiocOncoTK, biodbExpasy, biodbKegg, brainflowprobes, BRGenomics, CellaRepertorium, HTqPCR, microbiomeMarker, MQmetrics, nanotatoR, netOmics, Pi, polyester, psygenet2r, RandomWalkRestartMH, rDGIdb, Risa, RNAinteract, single, SummarizedBenchmark

**EXPERIMENT DATA:** 

Four experimental data packages were removed from this release (after being
deprecated in BioC 3.19):

-  MMAPPR2data, pcxnData, restfulSEData, RLHub

- Please note: CoSIAdata, previously announced as deprecated in 3.19, has been updated
and remain in Bioconductor.

Two experimental data packages are deprecated in this release and will be
removed in Bioc 3.21:

- DmelSGI, RNAinteractMAPK

**ANNOTATION DATA:** 

One annotation package was removed from this release:

- MafH5.gnomAD.v3.1.2.GRCh38

No annotation packages are deprecated in this release.


**WORKFLOWS:** 

No workflow packages were removed from this release.

No workflow packages were deprecated in this release.

**BOOKS:**

No books were removed from this release.

No books were deprecated in this release.
