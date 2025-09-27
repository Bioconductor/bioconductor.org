April 26, 2023

Bioconductors:

We are pleased to announce Bioconductor 3.17, consisting of
2230 software packages, 419 experiment data packages, 912 annotation
packages, 27 workflows and 3 books.

There are 79 new software packages, 7 new data experiment packages,
no new annotation packages, 2 new workflows, no new books, and many updates and
improvements to existing packages.

Bioconductor 3.17 is compatible with R 4.3, and is supported on Linux,
64-bit Windows, Intel 64-bit macOS 11 (Big Sur) or higher and
macOS arm64. This release will also include updated Bioconductor [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3] for details and downloads.

[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.17](#getting-started-with-bioconductor-316)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [New Workflow](#new-workflow-packages)
* [New Books](#new-online-books)
* [NEWS from existing software packages](#news-from-existing-software-packages)
* [NEWS from existing data experiment packages](#news-from-existing-data-experiment-packages)
* [NEWS from existing workflows](#news-from-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)


Getting Started with Bioconductor 3.17
======================================

To update to or install Bioconductor 3.17:

1. Install R 4.3. Bioconductor 3.17 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).


New Software Packages
=====================

There are 79 new software packages in this release of Bioconductor.

- [AHMassBank](/packages/AHMassBank) Supplies AnnotationHub with
  MassBank metabolite/compound annotations bundled in CompDb SQLite
  databases. CompDb SQLite databases contain general compound
  annotation as well as fragment spectra representing fragmentation
  patterns of compounds' ions. MassBank data is retrieved from
  https://massbank.eu/MassBank and processed using helper functions
  from the CompoundDb Bioconductor package into redistributable
  SQLite databases.

- [alabaster](/packages/alabaster) Umbrella for the alabaster suite,
  providing a single-line import for all alabaster.* packages.
  Installing this package ensures that all known alabaster.* packages
  are also installed, avoiding problems with missing packages when a
  staging method or loading function is dynamically requested.
  Obviously, this comes at the cost of needing to install more
  packages, so advanced users and application developers may prefer
  to install the required alabaster.* packages individually.

- [alabaster.base](/packages/alabaster.base) Save Bioconductor data
  structures into file artifacts, and load them back into memory.
  This is a more robust and portable alternative to serialization of
  such objects into RDS files. Each artifact is associated with
  metadata for further interpretation; downstream applications can
  enrich this metadata with context-specific properties.

- [alabaster.bumpy](/packages/alabaster.bumpy) Save BumpyMatrix
  objects into file artifacts, and load them back into memory. This
  is a more portable alternative to serialization of such objects
  into RDS files. Each artifact is associated with metadata for
  further interpretation; downstream applications can enrich this
  metadata with context-specific properties.

- [alabaster.mae](/packages/alabaster.mae) Save MultiAssayExperiments
  into file artifacts, and load them back into memory. This is a more
  portable alternative to serialization of such objects into RDS
  files. Each artifact is associated with metadata for further
  interpretation; downstream applications can enrich this metadata
  with context-specific properties.

- [alabaster.matrix](/packages/alabaster.matrix) Save matrices,
  arrays and similar objects into file artifacts, and load them back
  into memory. This is a more portable alternative to serialization
  of such objects into RDS files. Each artifact is associated with
  metadata for further interpretation; downstream applications can
  enrich this metadata with context-specific properties.

- [alabaster.ranges](/packages/alabaster.ranges) Save GenomicRanges,
  IRanges and related data structures into file artifacts, and load
  them back into memory. This is a more portable alternative to
  serialization of such objects into RDS files. Each artifact is
  associated with metadata for further interpretation; downstream
  applications can enrich this metadata with context-specific
  properties.

- [alabaster.sce](/packages/alabaster.sce) Save SingleCellExperiment
  into file artifacts, and load them back into memory. This is a more
  portable alternative to serialization of such objects into RDS
  files. Each artifact is associated with metadata for further
  interpretation; downstream applications can enrich this metadata
  with context-specific properties.

- [alabaster.schemas](/packages/alabaster.schemas) Stores all schemas
  required by various alabaster.* packages. No computation should be
  performed by this package, as that is handled by alabaster.base. We
  use a separate package instead of storing the schemas in
  alabaster.base itself, to avoid conflating management of the
  schemas with code maintenence.

- [alabaster.se](/packages/alabaster.se) Save SummarizedExperiments
  into file artifacts, and load them back into memory. This is a more
  portable alternative to serialization of such objects into RDS
  files. Each artifact is associated with metadata for further
  interpretation; downstream applications can enrich this metadata
  with context-specific properties.

- [alabaster.spatial](/packages/alabaster.spatial) Save
  SpatialExperiment objects and their images into file artifacts, and
  load them back into memory. This is a more portable alternative to
  serialization of such objects into RDS files. Each artifact is
  associated with metadata for further interpretation; downstream
  applications can enrich this metadata with context-specific
  properties.

- [alabaster.string](/packages/alabaster.string) Save Biostrings
  objects to file artifacts, and load them back into memory. This is
  a more portable alternative to serialization of such objects into
  RDS files. Each artifact is associated with metadata for further
  interpretation; downstream applications can enrich this metadata
  with context-specific properties.

- [alabaster.vcf](/packages/alabaster.vcf) Save variant calling
  SummarizedExperiment to file and load them back as VCF objects.
  This is a more portable alternative to serialization of such
  objects into RDS files. Each artifact is associated with metadata
  for further interpretation; downstream applications can enrich this
  metadata with context-specific properties.

- [AnVILWorkflow](/packages/AnVILWorkflow) The AnVIL is a cloud
  computing resource developed in part by the National Human Genome
  Research Institute. The main cloud-based genomics platform deported
  by the AnVIL project is Terra. The AnVILWorkflow package allows
  remote access to Terra implemented workflows, enabling end-user to
  utilize Terra/ AnVIL provided resources - such as data, workflows,
  and flexible/scalble computing resources - through the conventional
  R functions.

- [BG2](/packages/BG2) This package is built to perform GWAS analysis
  for non-Gaussian data using BG2. The BG2 method uses penalized
  quasi-likelihood along with nonlocal priors in a two step manner to
  identify SNPs in GWAS analysis. The research related to this
  package was supported in part by National Science Foundation awards
  DMS 1853549 and DMS 2054173.

- [BiocHail](/packages/BiocHail) Use hail via basilisk when
  appropriate, or via reticulate. This package can be used in
  terra.bio to interact with UK Biobank resources processed by
  hail.is.

- [BiocHubsShiny](/packages/BiocHubsShiny) A package that allows
  interactive exploration of AnnotationHub and ExperimentHub
  resources. It uses DT / DataTable to display resources for multiple
  organisms. It provides template code for reproducibility and for
  downloading resources via the indicated Hub package.

- [BSgenomeForge](/packages/BSgenomeForge) A set of tools to forge
  BSgenome data packages. Supersedes the old seed-based tools from
  the BSgenome software package. This package allows the user to
  create a BSgenome data package in one function call, simplifying
  the old seed-based process.

- [CBNplot](/packages/CBNplot) This package provides the
  visualization of bayesian network inferred from gene expression
  data. The networks are based on enrichment analysis results
  inferred from packages including clusterProfiler and ReactomePA.
  The networks between pathways and genes inside the pathways can be
  inferred and visualized.

- [cfTools](/packages/cfTools) The cfTools R package provides methods
  for cell-free DNA (cfDNA) methylation data analysis to facilitate
  cfDNA-based studies. Given the methylation sequencing data of a
  cfDNA sample, for each cancer marker or tissue marker, we
  deconvolve the tumor-derived or tissue-specific reads from all
  reads falling in the marker region. Our read-based deconvolution
  algorithm exploits the pervasiveness of DNA methylation for signal
  enhancement, therefore can sensitively identify a trace amount of
  tumor-specific or tissue-specific cfDNA in plasma. cfTools provides
  functions for (1) cancer detection: sensitively detect
  tumor-derived cfDNA and estimate the tumor-derived cfDNA fraction
  (tumor burden); (2) tissue deconvolution: infer the tissue type
  composition and the cfDNA fraction of multiple tissue types for a
  plasma cfDNA sample. These functions can serve as foundations for
  more advanced cfDNA-based studies, including cancer diagnosis and
  disease monitoring.

- [chihaya](/packages/chihaya) Saves the delayed operations of a
  DelayedArray to a HDF5 file. This enables efficient recovery of the
  DelayedArray's contents in other languages and analysis frameworks.

- [clevRvis](/packages/clevRvis) clevRvis provides a set of
  visualization techniques for clonal evolution. These include shark
  plots, dolphin plots and plaice plots. Algorithms for time point
  interpolation as well as therapy effect estimation are provided.
  Phylogeny-aware color coding is implemented. A shiny-app for
  generating plots interactively is additionally provided.

- [concordexR](/packages/concordexR) Many analysis workflows include
  approximation of a nearest neighbors graph followed by clustering
  of the graph structure. The concordex coefficient estimates the
  concordance between the graph and clustering results. The package
  'concordexR' is an R implementation of the original concordex
  Python-based command line tool.

- [CoSIA](/packages/CoSIA) Cross-Species Investigation and Analysis
  (CoSIA) is a package that provides researchers with an alternative
  methodology for comparing across species and tissues using normal
  wild-type RNA-Seq Gene Expression data from Bgee. Using RNA-Seq
  Gene Expression data, CoSIA provides multiple visualization tools
  to explore the transcriptome diversity and variation across genes,
  tissues, and species. CoSIA uses the Coefficient of Variation and
  Shannon Entropy and Specificity to calculate transcriptome
  diversity and variation. CoSIA also provides additional conversion
  tools and utilities to provide a streamlined methodology for
  cross-species comparison.

- [CTdata](/packages/CTdata) Data from publicly available databases
  (GTEx, CCLE, TCGA and ENCODE) that go with CTexploreR in order to
  re-define a comprehensive and thoroughly curated list of CT genes
  and their main characteristics.

- [cytofQC](/packages/cytofQC) cytofQC is a package for initial
  cleaning of CyTOF data. It uses a semi-supervised approach for
  labeling cells with their most likely data type (bead, doublet,
  debris, dead) and the probability that they belong to each label
  type. This package does not remove data from the dataset, but
  provides labels and information to aid the data user in cleaning
  their data. Our algorithm is able to distinguish between doublets
  and large cells.

- [CytoPipeline](/packages/CytoPipeline) This package provides
  support for automation and visualization of flow cytometry data
  analysis pipelines. In the current state, the package focuses on
  the preprocessing and quality control part. The framework is based
  on two main S4 classes, i.e. CytoPipeline and CytoProcessingStep.
  The pipeline steps are linked to corresponding R functions - that
  are either provided in the CytoPipeline package itself, or exported
  from a third party package, or coded by the user her/himself. The
  processing steps need to be specified centrally and explicitly
  using either a json input file or through step by step creation of
  a CytoPipeline object with dedicated methods. After having run the
  pipeline, obtained results at all steps can be retrieved and
  visualized thanks to file caching (the running facility uses a
  BiocFileCache implementation). The package provides also specific
  visualization tools like pipeline workflow summary display, and
  1D/2D comparison plots of obtained flowFrames at various steps of
  the pipeline.

- [cytoviewer](/packages/cytoviewer) This R package supports
  interactive visualization of multi-channel images and segmentation
  masks generated by imaging mass cytometry and other highly
  multiplexed imaging techniques using shiny. The cytoviewer
  interface is divided into image-level (Composite and Channels) and
  cell-level visualization (Masks). It allows users to overlay
  individual images with segmentation masks, integrates well with
  SingleCellExperiment and SpatialExperiment objects for metadata
  visualization and supports image downloads.

- [DELocal](/packages/DELocal) The goal of DELocal is to identify DE
  genes compared to their neighboring genes from the same chromosomal
  location. It has been shown that genes of related functions are
  generally very far from each other in the chromosome. DELocal
  utilzes this information to identify DE genes comparing with their
  neighbouring genes.

- [DESpace](/packages/DESpace) Intuitive framework for identifying
  spatially variable genes (SVGs) via edgeR, a popular method for
  performing differential expression analyses. Based on pre-annotated
  spatial clusters as summarized spatial information, DESpace models
  gene expression using a negative binomial (NB), via edgeR, with
  spatial clusters as covariates. SVGs are then identified by testing
  the significance of spatial clusters. The method is flexible and
  robust, and is faster than the most SV methods. Furthermore, to the
  best of our knowledge, it is the only SV approach that allows: -
  performing a SV test on each individual spatial cluster, hence
  identifying the key regions of the tissue affected by spatial
  variability; - jointly fitting multiple samples, targeting genes
  with consistent spatial patterns across replicates.

- [doubletrouble](/packages/doubletrouble) doubletrouble aims to
  identify duplicated genes from whole-genome protein sequences and
  classify them based on their modes of duplication. The duplication
  modes are: i. whole-genome duplication (WGD); ii. tandem
  duplication (TD); iii. proximal duplication (PD); iv. transposed
  duplication (TRD) and; v. dispersed duplication (DD). If users want
  a simpler classification scheme, duplicates can also be classified
  into WGD- and SSD-derived (small-scale duplication) gene pairs.
  Besides classifying gene pairs, users can also classify genes, so
  that each gene is assigned a unique mode of duplication. Users can
  also calculate substitution rates per substitution site (i.e., Ka
  and Ks) from duplicate pairs, find peaks in Ks distributions with
  Gaussian Mixture Models (GMMs), and classify gene pairs into age
  groups based on Ks peaks.

- [EDIRquery](/packages/EDIRquery) EDIRquery provides a tool to
  search for genes of interest within the Exome Database of
  Interspersed Repeats (EDIR). A gene name is a required input, and
  users can additionally specify repeat sequence lengths, minimum and
  maximum distance between sequences, and whether to allow a 1-bp
  mismatch. Outputs include a summary of results by repeat length, as
  well as a dataframe of query results. Example data provided
  includes a subset of the data for the gene GAA (ENSG00000171298).
  To query the full database requires providing a path to the
  downloaded database files as a parameter.

- [escheR](/packages/escheR) The creation of effective visualizations
  is a fundamental component of data analysis. In biomedical
  research, new challenges are emerging to visualize
  multi-dimensional data in a 2D space, but current data
  visualization tools have limited capabilities. To address this
  problem, we leverage Gestalt principles to improve the design and
  interpretability of multi-dimensional data in 2D data
  visualizations, layering aesthetics to display multiple variables.
  The proposed visualization can be applied to spatially-resolved
  transcriptomics data, but also broadly to data visualized in 2D
  space, such as embedding visualizations. We provide this open
  source R package escheR, which is built off of the state-of-the-art
  ggplot2 visualization framework and can be seamlessly integrated
  into genomics toolboxes and workflows.

- [FeatSeekR](/packages/FeatSeekR) FeatSeekR performs unsupervised
  feature selection using replicated measurements. It iteratively
  selects features with the highest reproducibility across
  replicates, after projecting out those dimensions from the data
  that are spanned by the previously selected features. The selected
  a set of features has a high replicate reproducibility and a high
  degree of uniqueness.

- [flowGate](/packages/flowGate) flowGate adds an interactive Shiny
  app to allow manual GUI-based gating of flow cytometry data in R.
  Using flowGate, you can draw 1D and 2D span/rectangle gates,
  quadrant gates, and polygon gates on flow cytometry data by
  interactively drawing the gates on a plot of your data, rather than
  by specifying gate coordinates. This package is especially geared
  toward wet-lab cytometerists looking to take advantage of R for
  cytometry analysis, without necessarily having a lot of R
  experience.

- [GeoTcgaData](/packages/GeoTcgaData) Gene Expression Omnibus(GEO)
  and The Cancer Genome Atlas (TCGA) provide us with a wealth of
  data, such as RNA-seq, DNA Methylation, SNP and Copy number
  variation data. It's easy to download data from TCGA using the gdc
  tool, but processing these data into a format suitable for
  bioinformatics analysis requires more work. This R package was
  developed to handle these data.

- [HiCExperiment](/packages/HiCExperiment) R generic interface to
  Hi-C contact matrices in `.(m)cool`, `.hic` or HiC-Pro derived
  formats, as well as other Hi-C processed file formats. Contact
  matrices can be partially parsed using a random access method,
  allowing a memory-efficient representation of Hi-C data in R. The
  `HiCExperiment` class stores the Hi-C contacts parsed from local
  contact matrix files. `HiCExperiment` instances can be further
  investigated in R using the `HiContacts` analysis package.

- [HiCool](/packages/HiCool) HiCool provides an R interface to
  process and normalize Hi-C paired-end fastq reads into .(m)cool
  files. .(m)cool is a compact, indexed HDF5 file format specifically
  tailored for efficiently storing HiC-based data. On top of
  processing fastq reads, HiCool provides a convenient reporting
  function to generate shareable reports summarizing Hi-C experiments
  and including quality controls.

- [IFAA](/packages/IFAA) This package offers a robust approach to
  make inference on the association of covariates with the absolute
  abundance (AA) of microbiome in an ecosystem. It can be also
  directly applied to relative abundance (RA) data to make inference
  on AA because the ratio of two RA is equal to the ratio of their
  AA. This algorithm can estimate and test the associations of
  interest while adjusting for potential confounders. The estimates
  of this method have easy interpretation like a typical regression
  analysis. High-dimensional covariates are handled with
  regularization and it is implemented by parallel computing. False
  discovery rate is automatically controlled by this approach. Zeros
  do not need to be imputed by a positive value for the analysis. The
  IFAA package also offers the 'MZILN' function for estimating and
  testing associations of abundance ratios with covariates.

- [INTACT](/packages/INTACT) This package integrates colocalization
  probabilites from colocalization analysis with transcriptome-wide
  association study (TWAS) scan summary statistics to implicate genes
  that may be biologically relevant to a complex trait. The
  probabilistic framework implemented in this package constrains the
  TWAS scan z-score-based likelihood using a gene-level
  colocalization probability. Given gene set annotations, this
  package can estimate gene set enrichment using posterior
  probabilities from the TWAS-colocalization integration step.

- [IntOMICS](/packages/IntOMICS) IntOMICS is an efficient integrative
  framework based on Bayesian networks. IntOMICS systematically
  analyses gene expression (GE), DNA methylation (METH), copy number
  variation (CNV) and biological prior knowledge (B) to infer
  regulatory networks. IntOMICS complements the missing biological
  prior knowledge by so-called empirical biological knowledge (empB),
  estimated from the available experimental data. An automatically
  tuned MCMC algorithm (Yang and Rosenthal, 2017) estimates model
  parameters and the empirical biological knowledge. Conventional
  MCMC algorithm with additional Markov blanket resampling (MBR) step
  (Su and Borsuk, 2016) infers resulting regulatory network structure
  consisting of three types of nodes: GE nodes refer to gene
  expression levels, CNV nodes refer to associated copy number
  variations, and METH nodes refer to associated DNA methylation
  probe(s).

- [magpie](/packages/magpie) This package aims to perform power
  analysis for the MeRIP-seq study. It calculates FDR, FDC, power,
  and precision under various study design parameters, including but
  not limited to sample size, sequencing depth, and testing method.
  It can also output results into .xlsx files or produce
  corresponding figures of choice.

- [mariner](/packages/mariner) Tools for manipulating paired ranges
  and working with Hi-C data in R. Functionality includes
  manipulating/merging paired regions, generating paired ranges,
  extracting/aggregating interactions from `.hic` files, and
  visualizing the results. Designed for compatibility with
  plotgardener for visualization.

- [mastR](/packages/mastR) mastR is an R package designed for
  automated screening of signatures of interest for specific research
  questions. The package is developed for generating refined lists of
  signature genes from multiple group comparisons based on the
  results from edgeR and limma differential expression (DE) analysis
  workflow. It also takes into account the background noise of
  tissue-specificity, which is often ignored by other marker
  generation tools. This package is particularly useful for the
  identification of group markers in various biological and medical
  applications, including cancer research and developmental biology.

- [mbQTL](/packages/mbQTL) mbQTL is a statistical R package for
  simultaneous 16srRNA,16srDNA (microbial) and variant, SNP, SNV
  (host) relationship, correlation, regression studies. We apply
  linear, logistic and correlation based statistics to identify the
  relationships of taxa, genus, species and variant, SNP, SNV in the
  infected host. We produce various statistical significance measures
  such as P values, FDR, BC and probability estimation to show
  significance of these relationships. Further we provide various
  visualization function for ease and clarification of the results of
  these analysis. The package is compatible with dataframe,
  MRexperiment and text formats.

- [microSTASIS](/packages/microSTASIS) The toolkit 'µSTASIS', or
  microSTASIS, has been developed for the stability analysis of
  microbiota in a temporal framework by leveraging on iterative
  clustering. Concretely, the core function uses Hartigan-Wong
  k-means algorithm as many times as possible for stressing out
  paired samples from the same individuals to test if they remain
  together for multiple numbers of clusters over a whole data set of
  individuals. Moreover, the package includes multiple functions to
  subset samples from paired times, validate the results or visualize
  the output.

- [MoleculeExperiment](/packages/MoleculeExperiment)
  MoleculeExperiment contains functions to create and work with
  objects from the new MoleculeExperiment class. We introduce this
  class for analysing molecule-based spatial transcriptomics data
  (e.g., Xenium by 10X, Cosmx SMI by Nanostring, and Merscope by
  Vizgen). This allows researchers to analyse spatial transcriptomics
  data at the molecule level, and to have standardised data formats
  accross vendors.

- [MsBackendSql](/packages/MsBackendSql) SQL-based mass spectrometry
  (MS) data backend supporting also storange and handling of very
  large data sets. Objects from this package are supposed to be used
  with the Spectra Bioconductor package. Through the MsBackendSql
  with its minimal memory footprint, this package thus provides an
  alternative MS data representation for very large or remote MS data
  sets.

- [MsDataHub](/packages/MsDataHub) The MsDataHub package uses the
  ExperimentHub infrastructure to distribute raw mass spectrometry
  data files, peptide spectrum matches or quantitative data from
  proteomics and metabolomics experiments.

- [MsQuality](/packages/MsQuality) The MsQuality provides
  functionality to calculate quality metrics for mass
  spectrometry-derived, spectral data at the per-sample level.
  MsQuality relies on the mzQC framework of quality metrics defined
  by the Human Proteom Organization-Proteomics Standards Initiative
  (HUPO-PSI). These metrics quantify the quality of spectral raw
  files using a controlled vocabulary. The package is especially
  addressed towards users that acquire mass spectrometry data on a
  large scale (e.g. data sets from clinical settings consisting of
  several thousands of samples). The MsQuality package allows to
  calculate low-level quality metrics that require minimum
  information on mass spectrometry data: retention time, m/z values,
  and associated intensities. MsQuality relies on the Spectra
  package, or alternatively the MsExperiment package, and its
  infrastructure to store spectral data.

- [MultimodalExperiment](/packages/MultimodalExperiment)
  MultimodalExperiment is an S4 class that integrates bulk and
  single-cell experiment data; it is optimally storage-efficient, and
  its methods are exceptionally fast. It effortlessly represents
  multimodal data of any nature and features normalized experiment,
  subject, sample, and cell annotations, which are related to
  underlying biological experiments through maps. Its coordination
  methods are opt-in and employ database-like join operations
  internally to deliver fast and flexible management of multimodal
  data.

- [OutSplice](/packages/OutSplice) An easy to use tool that can
  compare splicing events in tumor and normal tissue samples using
  either a user generated matrix, or data from The Cancer Genome
  Atlas (TCGA). This package generates a matrix of splicing outliers
  that are significantly over or underexpressed in tumors samples
  compared to normal denoted by chromosome location. The package also
  will calculate the splicing burden in each tumor and characterize
  the types of splicing events that occur.

- [pairedGSEA](/packages/pairedGSEA) pairedGSEA makes it simple to
  run a paired Differential Gene Expression (DGE) and Differencital
  Gene Splicing (DGS) analysis. The package allows you to store
  intermediate results for further investiation, if desired.
  pairedGSEA comes with a wrapper function for running an
  Over-Representation Analysis (ORA) and functionalities for plotting
  the results.

- [pfamAnalyzeR](/packages/pfamAnalyzeR) Protein domains is one of
  the most import annoation of proteins we have with the Pfam
  database/tool being (by far) the most used tool. This R package
  enables the user to read the pfam prediction from both webserver
  and stand-alone runs into R. We have recently shown most human
  protein domains exist as multiple distinct variants termed domain
  isotypes. Different domain isotypes are used in a cell, tissue, and
  disease-specific manner. Accordingly, we find that domain isotypes,
  compared to each other, modulate, or abolish the functionality of a
  protein domain. This R package enables the identification and
  classification of such domain isotypes from Pfam data.

- [planttfhunter](/packages/planttfhunter) planttfhunter is used to
  identify plant transcription factors (TFs) from protein sequence
  data and classify them into families and subfamilies using the
  classification scheme implemented in PlantTFDB. TFs are identified
  using pre-built hidden Markov model profiles for DNA-binding
  domains. Then, auxiliary and forbidden domains are used with
  DNA-binding domains to classify TFs into families and subfamilies
  (when applicable). Currently, TFs can be classified in 58 different
  TF families/subfamilies.

- [Rarr](/packages/Rarr) The Zarr specification defines a format for
  chunked, compressed, N-dimensional arrays.  It's design allows
  efficient access to subsets of the stored array, and supports both
  local and cloud storage systems. Rarr aims to implement this
  specifcation in R with minimal reliance on an external tools or
  libraries.

- [RBioFormats](/packages/RBioFormats) An R package which interfaces
  the OME Bio-Formats Java library to allow reading of proprietary
  microscopy image data and metadata.

- [Rcollectl](/packages/Rcollectl) Provide functions to obtain
  instrumentation data on processes in a unix environment.  Parse
  output of a collectl run.  Vizualize aspects of system usage over
  time, with annotation.

- [retrofit](/packages/retrofit) RETROFIT is a Bayesian non-negative
  matrix factorization framework to decompose cell type mixtures in
  ST data without using external single-cell expression references.
  RETROFIT outperforms existing reference-based methods in estimating
  cell type proportions and reconstructing gene expressions in
  simulations with varying spot size and sample heterogeneity,
  irrespective of the quality or availability of the single-cell
  reference. RETROFIT recapitulates known cell-type localization
  patterns in a Slide-seq dataset of mouse cerebellum without using
  any single-cell data.

- [ReUseData](/packages/ReUseData) ReUseData is an _R/Bioconductor_
  software tool to provide a systematic and versatile approach for
  standardized and reproducible data management. ReUseData
  facilitates transformation of shell or other ad hoc scripts for
  data preprocessing into workflow-based data recipes. Evaluation of
  data recipes generate curated data files in their generic formats
  (e.g., VCF, bed). Both recipes and data are cached using database
  infrastructure for easy data management and reuse. Prebuilt data
  recipes are available through ReUseData portal
  ("https://rcwl.org/dataRecipes/") with full annotation and user
  instructions. Pregenerated data are available through ReUseData
  cloud bucket that is directly downloadable through
  "getCloudData()".

- [rifiComparative](/packages/rifiComparative) 'rifiComparative' is a
  continuation of rifi package. It compares two conditions output of
  rifi using half-life and mRNA at time 0 segments. As an input for
  the segmentation, the difference between half-life of both
  condtions and log2FC of the mRNA at time 0 are used. The package
  provides segmentation, statistics, summary table, fragments
  visualization and some additional useful plots for further
  anaylsis.

- [S4Arrays](/packages/S4Arrays) The S4Arrays package defines the
  Array virtual class to be extended by other S4 classes that wish to
  implement a container with an array-like semantic. It also
  provides: (1) low-level functionality meant to help the developer
  of such container to implement basic operations like display,
  subsetting, or coercion of their array-like objects to an ordinary
  matrix or array, and (2) a framework that facilitates block
  processing of array-like objects (typically on-disk objects).

- [SCArray.sat](/packages/SCArray.sat) Extends the Seurat classes and
  functions to support Genomic Data Structure (GDS) files as a
  DelayedArray backend for data representation. It relies on the
  implementation of GDS-based DelayedMatrix in the SCArray package to
  represent single cell RNA-seq data. The common optimized algorithms
  leveraging GDS-based and single cell-specific DelayedMatrix
  (SC_GDSMatrix) are implemented in the SCArray package. This package
  introduces a new SCArrayAssay class (derived from the Seurat
  Assay), which wraps raw counts, normalized expressions and scaled
  data matrix based on GDS-specific DelayedMatrix. It is designed to
  integrate seamlessly with the Seurat package to provide common data
  analysis in the SeuratObject-based workflow. Compared with Seurat,
  SCArray.sat significantly reduces the memory usage and can be
  applied to very large datasets.

- [scFeatures](/packages/scFeatures) scFeatures constructs multi-view
  representations of single-cell and spatial data. scFeatures is a
  tool that generates multi-view representations of single-cell and
  spatial data through the construction of a total of 17 feature
  types. These features can then be used for a variety of analyses
  using other software in Biocondutor.

- [screenCounter](/packages/screenCounter) Provides functions for
  counting reads from high-throughput sequencing screen data (e.g.,
  CRISPR, shRNA) to quantify barcode abundance. Currently supports
  single barcodes in single- or paired-end data, and combinatorial
  barcodes in paired-end data.

- [scRNAseqApp](/packages/scRNAseqApp) scRNAseqApp is a Shiny app
  package that allows users to visualize single cell data
  interactively. It was modified from ShinyCell and repackaged to a
  tool to show multiple data. It can visulize the data with multiple
  information side by side.

- [scviR](/packages/scviR) This package defines interfaces from R to
  scvi-tools.  A vignette works through the totalVI tutorial for
  analyzing CITE-seq data.  Another vignette compares outputs of
  Chapter 12 of the OSCA book with analogous outputs based on totalVI
  quantifications. Future work will address other components of
  scvi-tools, with a focus on building understanding of probabilistic
  methods based on variational autoencoders.

- [seq.hotSPOT](/packages/seq.hotSPOT) seq.hotSPOT provides a
  resource for designing effective sequencing panels to help improve
  mutation capture efficacy for ultradeep sequencing projects. Using
  SNV datasets, this package designs custom panels for any tissue of
  interest and identify the genomic regions likely to contain the
  most mutations. Establishing efficient targeted sequencing panels
  can allow researchers to study mutation burden in tissues at high
  depth without the economic burden of whole-exome or whole-genome
  sequencing. This tool was developed to make high-depth sequencing
  panels to study low-frequency clonal mutations in clinically normal
  and cancerous tissues.

- [seqArchRplus](/packages/seqArchRplus) seqArchRplus facilitates
  downstream analyses of promoter sequence architectures/clusters
  identified by seqArchR (or any other tool/method). With additional
  available information such as the TPM values and interquantile
  widths (IQWs) of the CAGE tag clusters, seqArchRplus can order the
  input promoter clusters by their shape (IQWs), and write the
  cluster information as browser/IGV track files. Provided
  visualizations are of two kind: per sample/stage and per cluster
  visualizations. Those of the first kind include: plot panels for
  each sample showing per cluster shape, TPM and other score
  distributions, sequence logos, and peak annotations. The second
  include per cluster chromosome-wise and strand distributions, motif
  occurrence heatmaps and GO term enrichments. Additionally,
  seqArchRplus can also generate HTML reports for easy viewing and
  comparison of promoter architectures between samples/stages
  (future).

- [SGCP](/packages/SGCP) SGC is a semi-supervised pipeline for gene
  clustering in gene co-expression networks. SGC consists of multiple
  novel steps that enable the computation of highly enriched modules
  in an unsupervised manner. But unlike all existing frameworks, it
  further incorporates a novel step that leverages Gene Ontology
  information in a semi-supervised clustering method that further
  improves the quality of the computed modules.

- [SiPSiC](/packages/SiPSiC) Infer biological pathway activity of
  cells from single-cell RNA-sequencing data by calculating a pathway
  score for each cell (pathway genes are specified by the user). It
  is recommended to have the data in Transcripts-Per-Million (TPM) or
  Counts-Per-Million (CPM) units for best results. Scores may change
  when adding cells to or removing cells off the data. SiPSiC stands
  for Single Pathway analysis in Single Cells.

- [SparseArray](/packages/SparseArray) The SparseArray package
  defines the SparseArray virtual class to be extended by other S4
  classes that wish to represent in-memory multidimensional sparse
  arrays. One such extension is the SVT_SparseArray class, also
  defined in the package, that provides an efficient representation
  of the nonzero multidimensional data via a novel layout called the
  "SVT layout". SVT_SparseArray objects mimic the behavior of
  ordinary matrices or arrays in R as much as possible. In
  particular, they suppport most of the "standard array API" defined
  in base R.

- [SpatialOmicsOverlay](/packages/SpatialOmicsOverlay) Tools for
  NanoString Technologies GeoMx Technology. Package to easily graph
  on top of an OME-TIFF image. Plotting annotations can range from
  tissue segment to gene expression.

- [speckle](/packages/speckle) The speckle package contains functions
  for the analysis of single cell RNA-seq data. The speckle package
  currently contains functions to analyse differences in cell type
  proportions. There are also functions to estimate the parameters of
  the Beta distribution based on a given counts matrix, and a
  function to normalise a counts matrix to the median library size.
  There are plotting functions to visualise cell type proportions and
  the mean-variance relationship in cell type proportions and counts.
  As our research into specialised analyses of single cell data
  continues we anticipate that the package will be updated with new
  functions.

- [SVMDO](/packages/SVMDO) It is an easy-to-use GUI using disease
  information for detecting tumor/normal sample discriminating gene
  sets from differentially expressed genes. Our approach is based on
  an iterative algorithm filtering genes with disease ontology
  enrichment analysis and wilk’s lambda criterion connected to SVM
  classification model construction. Along with gene set extraction,
  SVMDO also provides individual prognostic marker detection. The
  algorithm is designed for FPKM and RPKM normalized RNA-Seq
  transcriptome datasets.

- [TDbasedUFE](/packages/TDbasedUFE) This is a comprehensive package
  to perform Tensor decomposition based unsupervised feature
  extraction. It can perform unsupervised feature extraction. It uses
  tensor decomposition. It is applicable to gene expression, DNA
  methylation, and histone modification etc. It can perform
  multiomics analysis. It is also potentially applicable to single
  cell omics data sets.

- [TDbasedUFEadv](/packages/TDbasedUFEadv) This is an advanced
  version of TDbasedUFE, which is a comprehensive package to perform
  Tensor decomposition based unsupervised feature extraction. In
  contrast to TDbasedUFE which can perform simple the feature
  selection and the multiomics analyses, this package can perform
  more complicated and advanced features, but they are not so
  popularly required. Only users who require more specific features
  can make use of its functionality.

- [TOP](/packages/TOP) TOP constructs a transferable model across
  gene expression platforms for prospective experiments. Such a
  transferable model can be trained to make predictions on
  independent validation data with an accuracy that is similar to a
  re-substituted model. The TOP procedure also has the flexibility to
  be adapted to suit the most common clinical response variables,
  including linear response, binomial and Cox PH models.

- [ZygosityPredictor](/packages/ZygosityPredictor) The
  ZygosityPredictor allows to predict how many copies of a gene are
  affected by small variants. In addition to the basic calculations
  of the affected copy number of a variant, the Zygosity-Predictor
  can integrate the influence of several variants on a gene and
  ultimately make a statement if and how many wild-type copies of the
  gene are left. This information proves to be of particular use in
  the context of translational medicine. For example, in cancer
  genomes, the Zygosity-Predictor can address whether unmutated
  copies of tumor-suppressor genes are present. Beyond this, it is
  possible to make this statement for all genes of an organism. The
  Zygosity-Predictor was primarily developed to handle SNVs and
  INDELs (later addressed as small-variants) of somatic and germline
  origin. In order not to overlook severe effects outside of the
  small-variant context, it has been extended with the assessment of
  large scale deletions, which cause losses of whole genes or parts
  of them.

New Data Experiment Packages
=====================

There are 7 new data experiment packages in this release of Bioconductor.

- [CoSIAdata](/packages/CoSIAdata) Variance Stabilized Transformation
  of Read Counts derived from Bgee RNA-Seq Expression Data.
  Expression Data includes annotations and is across 6 species (Homo
  sapiens, Mus musculus, Rattus norvegicus, Danio rerio, Drosophila
  melanogaster, and Caenorhabditis elegans) and across more than 132
  tissues. The data is represented as a RData files and is available
  in ExperimentHub.

- [DNAZooData](/packages/DNAZooData) DNAZooData is a data package
  giving programmatic access to genome assemblies and Hi-C contact
  matrices uniformly processed by the [DNA Zoo
  Consortium](https://www.dnazoo.org/). The matrices are available in
  the multi-resolution `.hic` format. A URL to corrected genome
  assemblies in `.fastq` format is also provided to the end-user.

- [fourDNData](/packages/fourDNData) fourDNData is a data package
  giving programmatic access to Hi-C contact matrices uniformly
  processed by the [4DN consortium](https://www.4dnucleome.org/). The
  matrices are available in the multi-resolution `.mcool` format.

- [gDNAinRNAseqData](/packages/gDNAinRNAseqData) Provides access to
  BAM files generated from RNA-seq data produced with different
  levels of gDNA contamination. It currently allows one to download a
  subset of the data published by Li et al., BMC Genomics, 23:554,
  2022. This subset of data is formed by BAM files with about 100,000
  alignments with three different levels of gDNA contamination.

- [marinerData](/packages/marinerData) Subsampled Hi-C in HEK cells
  expressing the NHA9 fusion with an F to S mutated IDR ("FS") or
  without any mutations to the IDR ("Wildtype" or "WT"). These files
  are used for testing mariner functions and some examples.

- [MetaScope](/packages/MetaScope) This package contains tools and
  methods for preprocessing microbiome data. Functionality includes
  library generation, demultiplexing, alignment, and microbe
  identification.  It is partly an R translation of the PathoScope
  2.0 pipeline.

- [scMultiome](/packages/scMultiome) Single cell multiome data,
  containing chromatin accessibility (scATAC-seq) and gene expression
  (scRNA-seq) information analyzed with the ArchR package and
  presented as MultiAssayExperiment objects.


New Annotation Packages
=====================

There are no new annotation packages.


New Workflow Packages
=====================

There are 2 new workflow packages in this release of Bioconductor.

- [seqpac](/packages/seqpac) Seqpac provides functions and workflows
  for analysis of short sequenced reads. It was originally developed
  for small RNA analysis, but can be implemented on any sequencing
  raw data (provided as a fastq-file), where the unit of measurement
  is counts of unique sequences. The core of the seqpac workflow is
  the generation and subsequence analysis/visualization of a
  standardized object called PAC. Using an innovative targeting
  system, Seqpac process, analyze and visualize sample or sequence
  group differences using the PAC object. A PAC object in its most
  basic form is a list containing three types of data frames. -
  Phenotype table (P): Sample names (rows) with associated metadata
  (columns) e.g. treatment. - Annotation table (A): Unique sequences
  (rows) with annotation (columns), eg. reference alignments. -
  Counts table (C): Counts of unique sequences (rows) for each sample
  (columns). The PAC-object follows the rule: - Row names in P must
  be identical with column names in C. - Row names in A must be
  identical with row names in C. Thus P and A describes the columns
  and rows in C, respectively. The targeting system, will either
  target specific samples in P (pheno_target) or sequences in A
  (anno_target) and group them according to a target column in P and
  A, respectively (see vignettes for more details).

- [spicyWorkflow](/packages/spicyWorkflow) We have developed an
  analytical framework for analysing data from high dimensional in
  situ cytometry assays including CODEX, CycIF, IMC and High
  Definition Spatial Transcriptomics. This framework makes use of
  functionality from our Bioconductor packages spicyR, lisaClust,
  scFeatures, FuseSOM, simpleSeg and ClassifyR and contains most of
  the key steps which are needed to interrogate the comprehensive
  spatial information generated by these exciting new technologies
  including cell segmentation, feature normalisation, cell type
  identification, micro-environment characterisation, spatial
  hypothesis testing and patient classification. Ultimately, our
  modular analysis framework provides a cohesive and accessible entry
  point into spatially resolved single cell data analysis for any
  R-based bioinformatician.

New Online Books
=====================

There are no new online books.


NEWS from existing Software Packages
===================================


[ADaCGH2](/packages/ADaCGH2)
-------

                 Changes in version 2.39.1 (2023-04-14)                 

- Fixed "DLL requires the use of native symbols" with R-4.3.0.

[affxparser](/packages/affxparser)
----------

                 Changes in version 1.73.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor release version,
which is now Bioconductor 3.18 for R (>= 4.4.0).

                 Changes in version 1.72.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor release version,
which is now Bioconductor 3.17 for R (>= 4.3.0).

                 Changes in version 1.71.2 (2023-04-23)                 

Bug Fixes

- fix to src/_mingw.h provided by Tomas Kalibera.

                 Changes in version 1.71.1 (2023-04-04)                 

Bug Fixes

- Fix two instances of "watching polymorphic type 'class Except' by
value &#91;-Wcatch-value=&#93;" compiler warnings.

                 Changes in version 1.71.0 (2022-11-01)                 

Notes

- The version number was bumped for the Bioconductor devel version,
which is now Bioconductor 3.17 for R devel.

[AHMassBank](/packages/AHMassBank)
----------

                        Changes in version 0.99                         

AHMassBank 0.99.3

- Add package documentation.
- Increase required R version.

AHMassBanki 0.99.0

- Add metadata for MassBank 2022.12.1.

[alevinQC](/packages/alevinQC)
--------

                       Changes in version 1.15.1                        

- Support BiocStyle::pdf_document and BiocStyle::html_document output
formats (thanks to Mike Smith for insights regarding BiocStyle)

[AlpsNMR](/packages/AlpsNMR)
-------

                 Changes in version 4.1.6 (2023-02-16)                  

- Download improvements:
- Progress bar
- Detect user interruptions
- Sleep 3 seconds between retrying failed downloads

                 Changes in version 4.1.5 (2023-02-10)                  

- Replace deprecated dplyr::select() calls.
- Remove workaround for mixOmics bug, bump mixOmics dependency
- Simplify implementation (same algorithm) for determining the
optimal
number of latent variables in the plsda models.
- Bump dplyr dependency version.

                 Changes in version 4.1.4 (2022-11-08)                  

- Disable nested parallelization
- Update workaround Biocparallel bpmapply
- Revert bpstop to sleep

                 Changes in version 4.1.3 (2022-11-07)                  

- Closer to the fix

                 Changes in version 4.1.2 (2022-11-04)                  

- Try a more robust fix on palomino4 (bpstop() instead sleep)
- Use register() in an example to avoid further breakage on palomino
- Workaround performance issues on BiocParallel::bpmapply()
(https://github.com/Bioconductor/BiocParallel/pull/228)

                 Changes in version 4.1.1 (2022-11-02)                  

- Remove archive dependency
- Try fixing build on palomino4, due to race condition in R CMD check

[AneuFinder](/packages/AneuFinder)
----------

                       Changes in version 1.27.1                        

SIGNIFICANT USER-LEVEL CHANGES

- Removed parameter "classes" from clusterHMMs(). This is due to a
  missing dependency (the ReorderCluster package). Sorry!

- Removed parameter "reorder.by.class" from heatmapGenomewide(). This
  is due to a missing dependency (the ReorderCluster package). Sorry!

[AnnotationForge](/packages/AnnotationForge)
---------------

                       Changes in version 1.41.0                        

MODIFICATIONS

- 1.41.5 Add RUnit to Suggest

- 1.41.3 Update viableID.rda for upcoming release

- 1.41.2 Use dbExecute() instead of dbGetQuery() or dbSendQuery() to
  avoid warnings

ENHANCEMENT

- 1.41.4 Knoknok Optimizations and allow specification of the desired
  Ensembl release

- 1.41.1 James MacDonald Fixed makeOrgPackageFromNCBI to use existing
  downloaded files

[AnnotationHub](/packages/AnnotationHub)
-------------

                        Changes in version 3.7.0                        

NEW FEATURES

- (3.7.4) Suppress snapshot date unless interactive session

- (3.7.1) Add DispatchClass for keras model weights

DEPRECATED

- (3.7.2) Deprecate display in favor of BiocHubsShiny; see vignette

[AnnotationHubData](/packages/AnnotationHubData)
-----------------

                       Changes in version 1.29.0                        

NEW FEATURES

- 1.29.2 Added HIC as acceptable source type

- 1.29.1 Added CDF as acceptable source type

[AnVIL](/packages/AnVIL)
-----

                       Changes in version 1.12.0                        

USER VISIBLE CHANGES

- (v 1.11.2) update workflow file discovery to use API, rather than
'scraping' google bucket.
https://github.com/Bioconductor/AnVIL/issues/69

- (v 1.11.3) Gen3 services deprecated

- (v 1.11.5) Add na = to handle NA encoding in avtable() /
avtable_import(). Changes default behavior.
https://github.com/Bioconductor/AnVIL/issues/75

BUG FIXES

- (v 1.11.1) consistently URLencode workspace and workflow name, to
allow for spaces. https://github.com/Bioconductor/AnVIL/issues/67

[AnVILWorkflow](/packages/AnVILWorkflow)
-------------

                        Changes in version 1.0.0                        

- Initial release of the 'AnVILWorkflow' package

[aroma.light](/packages/aroma.light)
-----------

                 Changes in version 3.31.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor release version,
which is now Bioconductor 3.18 for R (>= 4.4.0).

                 Changes in version 3.30.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor release version,
which is now Bioconductor 3.17 for R (>= 4.3.0).

                 Changes in version 3.29.0 (2022-11-01)                 

[ArrayExpress](/packages/ArrayExpress)
------------

                       Changes in version 1.58.0                        

- Update to access ArrayExpress Collection at BioStudies.

[ATACseqQC](/packages/ATACseqQC)
---------

                       Changes in version 1.23.1                        

- Add smooth function to TSSEscore.
- Add pseudoPausingIndex function.

[ATACseqTFEA](/packages/ATACseqTFEA)
-----------

                        Changes in version 1.1.1                        

- fix the bugs in sample_scripts.R.

[atena](/packages/atena)
-----

                 Changes in version 1.6.0 (2023-04-20)                  

USER VISIBLE CHANGES

- Improvement of the accuracy of the atena expression quantification
  method.

- Fixed numerical instability in TEtranscripts method.

- Added function (.matchSeqinfo) to harmonize seqinfo() of object with
  alignments and object with feature annotations.

- Implemented 'OneCodeToFindThemAll' annotation parser for RepeatMasker
  annotations.

- Implemented 'atena' annotation parser for RepeatMasker annotations.

- Added examples in the vignettes of TE annotation preprocessing steps
  using the implemented parsers.

- Changed default value of geneFeatures to NULL in the parameters
  objects.

BUG FIXES

- Fixed bug for not properly paired reads.

[BASiCS](/packages/BASiCS)
------

                 Changes in version 2.11.9 (2023-03-07)                 

- Fixes minor bug in test_denoised.R

                 Changes in version 2.11.8 (2023-03-07)                 

- fixes conflict after merge and version bump only

                 Changes in version 2.11.7 (2023-03-07)                 

- Updates NEWS

- version bump to trigger new build

- Fixes small bugs in tests (test_misc.R and test_divide_and_conquer.R)

- Minor changes suggested by BiocCheck

                 Changes in version 2.11.6 (2023-02-07)                 

- Adds `BASiCS_CalculateERCC()` to NAMESPACE

- Creates unit test for `BASiCS_CalculateERCC()`

                 Changes in version 2.11.5 (2023-01-04)                 

- Changes maintainer

- Remove unicode mu in documentation

                 Changes in version 2.11.4 (2022-12-15)                 

- Solves bug in `BASiCS_DiagPlot()`

                 Changes in version 2.11.3 (2022-12-15)                 

- Creates `BASiCS_CalculateSpikeIns()` to perform spike-in calculation
  (depends on the `scRNAseq` package).

- Adds extra option to plotting functions for MCMC diagnostics (e.g.
  `BASiCS_DiagPlot()`)

                 Changes in version 2.11.2 (2022-11-23)                 

- Bugfix in DetectLVG/HVG, see Issue #265; delta column in results was
  actually
  duplicated mu column.

- Change BASiCS_CorrectOffset to also alter the scale of the
  normalisation
  factors when WithSpikes=FALSE.

- Add "rhat" option to BASiCS_DiagHist and BASiCS_DiagPlot

[basilisk.utils](/packages/basilisk.utils)
--------------

                       Changes in version 1.12.0                        

- Added Arm64 support for Linux.

[batchelor](/packages/batchelor)
---------

                       Changes in version 1.16.0                        

- Bugfix to rownames of mnnCorrect() output when using
  correct.all=TRUE with subset.row=.

[benchdamic](/packages/benchdamic)
----------

                 Changes in version 1.5.3 (2023-04-20)                  

- Add CITATION file

- Update README.md file

- Add ORCID for authors

- Bug-fix for DA_MAST() function

                 Changes in version 1.5.2 (2022-11-21)                  

- Correct FDR description in vignette

- Add BiocParallel support for runMocks() and runSplits() functions

- Add example for ANCOMBC based methods in parallel

                 Changes in version 1.5.1 (2022-11-14)                  

- Add FDR computation in Type I Error Control analysis

- Add BiocParallel support for runMocks() and runSplits() functions

- Update vignette

                 Changes in version 1.5.0 (2022-11-01)                  

- New devel version

[BiocCheck](/packages/BiocCheck)
---------

                       Changes in version 1.36.0                        

NEW FEATURES

- Include size limit checks for data files in `data`, `inst/extdata`,
  and
  `data-raw` folders (@lshep, @const-ae, #167, #67)

- Source package directories that include an `inst/doc` folder with
  files
  are now flagged with an error. `doc` folders are generated during
  `R CMD build`.

- The error for packages already hosted on CRAN has been converted to a
  warning (@lshep, #177). Any such incoming packages must be removed
  from CRAN
  before the following Bioconductor release.

BUG FIXES AND MINOR IMPROVEMENTS

- Filter out 'package' docTypes from '\value' documentaiton checks
  (@grimbough, #189)

- Obtain a more complete list of deprecated packages for
  `checkDeprecatedPkgs`

- Fix issue with path seperators on Windows ('\\' vs '/') causing the
  unit
  test for `getBiocCheckDir` to report erroneous mismatches
  (@grimbough, #175)

- Fix bug where the wrong number of functions with length greater than
  50
  was reported (@grimbough, #182)

- biocViews term suggestions should be a scalar character
  (@lcolladotor,
  #184)

- Update email in bioc-devel subscription check (@lshep, #185)

- Handle function length checks when there is no R code (@lshep, #186)

- Edit warning text when empty or missing `value` sections are found in
  a package (@vjcitn)

- `checkForValueSection` re-implemented for clarity; filters out
  comments
  in `value` sections (@LiNk-NY)

- Correctly identify `Rd` comments by updating the regular expression
  to
  identify them (@LiNk-NY)

[BiocFileCache](/packages/BiocFileCache)
-------------

                         Changes in version 2.7                         

USER VISIBLE CHANGE

- (2.7.1) Remove rappdirs officially from package after deprecating
  functionality and moving to using default R tools caching location.

BUG FIX

- (2.7.2) Fix mismatch of arguments from function to generic for
  bfcupdate

[BiocHail](/packages/BiocHail)
--------

                        Changes in version 1.0.0                        

- disabled mac builds in Bioconductor build system, but
  the package should be usable with JDK version <= 11
  

[BiocHubsShiny](/packages/BiocHubsShiny)
-------------

                        Changes in version 1.0.0                        

- Added a NEWS.md file to track changes to the package.

[BiocIO](/packages/BiocIO)
------

                       Changes in version 1.10.0                        

Significant user-visible changes

- All the documentation has been updated to be more user-friendly.

[BiocParallel](/packages/BiocParallel)
------------

                        Changes in version 1.34                         

NEW FEATURES

- (1.33.2) limit worker number via environment variables.
  https://github.com/Bioconductor/BiocParallel/issues/229

- (v1.33.3) bpmapply() does not send the whole list of arguments
  to all workers. Instead, it takes the arguments and slices them,
  passing the corresponding slice to each worker. Thanks Sergio Oller!
  https://github.com/Bioconductor/BiocParallel/issues/229

USER VISIBLE CHANGES

- (1.33.1) Mark BatchJobsParam, bprunMPIslave as defunct.

- (1.33.9) Change default force.GC= to FALSE in MulticoreParam().
  <https://github.com/Bioconductor/BiocParallel/issues/238>

- (1.33.11) change content of 'traceback' on error to include the
  stack from the location of the error up to the invokation of
  FUN. Previously, the traceback was from FUN to the top-level of
  worker code, providing limited insight into nested errors.

- (1.33.12) 'force' function arguments to avoid consequences of
  lazy evaluation discussed in
  <https://github.com/Bioconductor/BiocParallel/issues/241#issuecomment-1445006892>

BUG FIXES

- (1.33.6) Restore 'exported' global variables in SerialParam()
  https://github.com/Bioconductor/BiocParallel/issues/234

- (1.33.7) 'configure.ac' uses C++ compiler and checks for existence
  of required header
  <https://github.com/Bioconductor/BiocParallel/pull/236>

- (1.33.8 / v 1.32.5) set socket idle timeout to a large value, to
  avoid premature worker termination and to be consistent with snow
  / parallel defaults.
  <https://github.com/Bioconductor/BiocParallel/issues/237>

- (1.33.10 / v 1.32.6) be sure to clean up TransientMulticoreParam
  state at start of each job.
  <https://github.com/Bioconductor/BiocParallel/issues/243>

[BiocPkgTools](/packages/BiocPkgTools)
------------

                       Changes in version 1.18.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- `biocLastBuildDate` has been deprecated and moved to `BiocArchive`.

- `pkgBiocDeps` returns only the Bioconductor dependencies for a given
  vector of packages; can be set to all dependencies.

BUG FIXES

- Use `bfcdownload` when a web resource needs updating. Functions that
  download and cache files are affected including `biocBuildReport`,
  `biocDownloadStats`, `biocPkgRanges`, and `CRANstatus`.

- `biocPkgList` throws a descriptive error when using an alternative
  repository to the default CRAN, e.g., `RSPM` where the `packages.rds`
  file does not exist.

[BiocSet](/packages/BiocSet)
-------

                        Changes in version 1.14                         

BUG FIX

- (1.12.1 / 1.13.1) propagate metadata with map_*() functions;
  reexport metadata(), metadata()<- generics and methods for public
  use <https://github.com/Bioconductor/BiocSet/issues/6>

[BiocStyle](/packages/BiocStyle)
---------

                       Changes in version 2.28.0                        

BUG FIXES

- Addressed issue where, when using newer versions of Pandoc, footnotes
  would appear at the bottom of HTML output rather than be moved to
  the margin.

- Fixed problem including references in R Markdown documents when the
  output format was BiocStyle::pdf_document
  (https://github.com/Bioconductor/BiocStyle/issues/94)

[biocthis](/packages/biocthis)
--------

                        Changes in version 1.9.2                        

BUG FIXES

- use_bioc_github_action() now properly works again when docker =
TRUE. Behind the scenes, this function now uses
docker/build-docker-action@v4 instead of the deprecated
docker/build-docker-action@v1. These updates were tested at
https://github.com/lcolladotor/ExampleBiocWorkshop2023.

[biocViews](/packages/biocViews)
---------

                       Changes in version 1.67.0                        

ENHANCEMENT

- (1.67.1) Add biocViews term WorkflowManagement

- (1.67.3) Add biocViews term LongRead

[biodbHmdb](/packages/biodbHmdb)
---------

                 Changes in version 1.5.1 (2023-03-22)                  

- Remove R_front script.

- Define new field secondary.accessions.

[biomaRt](/packages/biomaRt)
-------

                       Changes in version 2.56.0                        

BUG FIXES

- Fix problem when multiple cache entries with the same ID could
  be created.  (Thanks to Hervé Pagès & Henrik Bengtsson for
  independent
  reports of this issue.)

- bmRequest() will now respect the setting in options("timeout")

[BioNAR](/packages/BioNAR)
------

                        Changes in version 1.2.1                        

- Allow analysis of directed graphs and add four new centrality
measures specific for the directed graphs.
- Modify annotation functions to allow annotate nodes not only by its
name propety but by any other arbitrary property if its value
uniquely identify the node.
- Added new framework for the visualisation and analysis of
enrichment
results.

[bluster](/packages/bluster)
-------

                       Changes in version 1.10.0                        

- Corrected calculation of the ratio matrix in pairwiseRand()
  when adjusted=TRUE, so that ratios have an upper bound of 1.

[BRGenomics](/packages/BRGenomics)
----------

                       Changes in version 1.11.1                        

- Updated the import_bam_ATAC() function with slightly revised
options
and operations
- Fixed broken links in documentation

[BridgeDbR](/packages/BridgeDbR)
---------

                        Changes in version 2.9.2                        

NEW FEATURES

- Updated to BridgeDb 3.0.21

                        Changes in version 2.9.1                        

NEW FEATURES

- Updated to BridgeDb 3.0.19

[BSgenomeForge](/packages/BSgenomeForge)
-------------

                        Changes in version 1.0.0                        

- First version of the package that is ready for general use.

[CAGEr](/packages/CAGEr)
-----

                        Changes in version 2.6.0                        

BUG FIXES

- Re-enable the promoter shift functions (Fixes #4).
- Compute dominant TSS information for consensus clusters (Fixes #1).

NEW FEATURES

- Inter-quantile width is also computed when no sample is selected.
- Enhancer detection using a wrapper to the CAGEfighter package's
function quickEnhancers().

[Cardinal](/packages/Cardinal)
--------

                 Changes in version 3.0.1 (2022-11-14)                  

SIGNIFICANT USER-VISIBLE CHANGES

- Overwriting existing MSI files is now a warning instead of an error

BUG FIXES

- Fixed issue in 'peakAlign()' reference m/z's being sort and unique

[cBioPortalData](/packages/cBioPortalData)
--------------

                       Changes in version 2.12.0                        

Bug fixes and minor improvements

- Update successrate thresholds and fix long tests

[CBNplot](/packages/CBNplot)
-------

                 Changes in version 0.99.5 (2023-04-13)                 

- Changed the dependency for R version

                 Changes in version 0.99.4 (2023-04-07)                 

- Fixed the WARNINGS in package builder at bioconductor.org

                 Changes in version 0.99.3 (2023-04-07)                 

- Fixed the problems raised in the review

                 Changes in version 0.99.1 (2022-02-24)                 

- Fixed errors and warnings raised in pre-checking

                 Changes in version 0.99.0 (2021-07-29)                 

- Prepared for submission to Bioconductor

[celda](/packages/celda)
-----

                 Changes in version 1.14.2 (2023-01-19)                 

- Update to match Bioconductor release version

[cellxgenedp](/packages/cellxgenedp)
-----------

                         Changes in version 1.4                         

SIGNIFICANT USER-VISIBLE CHANGES

- (v 1.3.3) add publisher_metadata(), authors(), and links() to make
access to nested 'collections()' data more straight-forward

[CeTF](/packages/CeTF)
----

                       Changes in version 1.11.1                        

- Update for Biocondutor 3.17 release

[cfDNAPro](/packages/cfDNAPro)
--------

                        Changes in version 1.5.4                        

- In addition to "bam" and "picard" files as the input, now we accept
"cfdnapro" as input_type to various functions, this 'cfdnapro' input
is exactly the output of read_bam_insert_metrics function in
cfDNAPro package. It is a tsv file containing two columns, i.e.,
"insert_size" (fragment length) and "All_Reads.fr_count" (the count
of the fragment length).

                        Changes in version 1.5.3                        

- added support for hg38-NCBI version, i.e. GRCh38

[cfTools](/packages/cfTools)
-------

                        Changes in version 1.0.0                        

- 1st version of the package

- Submitted to Bioconductor

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

                       Changes in version 3.33.1                        

- update the documentation for enrichment analysis

- update the disjointExons to exonicParts for ensembldb

[ChIPQC](/packages/ChIPQC)
------

                       Changes in version 1.35.1                        

- rollup bugfixes

[ChIPseeker](/packages/ChIPseeker)
----------

                       Changes in version 1.35.3                        

- fixed R check by removing calling BiocStyle::Biocpkg() in vignette,
instead we use yulab.utils::Biocpkg() (2023-04-11, Tue)

                       Changes in version 1.35.2                        

- fixed R check by adding 'prettydoc' to Suggests (2023-04-04, Tue)

                       Changes in version 1.35.1                        

- use ggplot to plot heatmap (2022-12-30, Fri, #203)
- update startup message to display the 'Current Protocols (2022)'
paper.

[ClassifyR](/packages/ClassifyR)
---------

                        Changes in version 3.4.0                        

- 
  Companion website with more in-depth explanation and examples.

- 
  Default random forest classifier based on ranger now does
  two-step classification; one for variable importance and one
  for model fitting, as recommended by ranger's developer.

- 
  More functions use automatic parameter value selection as their
  defaults.

- 
  randomSelection function to choose random sets of features in
  cross-validation.

- 
  crossValidate function now permits custom parameter tuning via
  extraParams.

- 
  Elastic net GLM and ordinary GLM now calculate class weights be
  default, so as to perform well in class-imbalanced scenarios.

- 
  precisionPathwaysTrain and precisionPathwaysPredict functions
  for building tree-like models of multiple assays and their
  accessory functions calcCostsAndPerformance, bubblePlot,
  flowchart, strataPlot for model performance evaluation.

- 
  crissCrossValidate function that takes a list of data sets with
  the same set of features and the same set of outcomes and does
  all possible pairs of training and prediction to evaluate
  generalisability. crissCrossPlot for visual evaluation.

- 

[clevRvis](/packages/clevRvis)
--------

                 Changes in version 0.99.5 (2023-01-16)                 

- Fixed notes from BiocCheck

                 Changes in version 0.99.4 (2023-01-04)                 

- Included reviewer feedback

                 Changes in version 0.99.3 (2022-12-20)                 

- Removed R project file

                 Changes in version 0.99.2 (2022-12-06)                 

- Information on data added to vignette

- Show method added for seaObject

                 Changes in version 0.99.1 (2022-11-26)                 

- Initial submission to Bioconductor

[clusterProfiler](/packages/clusterProfiler)
---------------

                        Changes in version 4.7.1                        

- update according to the KEGG api changes (2023-03-01, Wed)

[CNVfilteR](/packages/CNVfilteR)
---------

                       Changes in version 1.13.2                        

MINOR

- Fixed bug: sample name column in CNV calling file is always
  interpreted as character, which prevents later errors.

                       Changes in version 1.13.1                        

MINOR

- Vignette: genome is explicitly shown as an option in loadVCFs()

- Vignette: mutiallelic sites are not currently supported and bcftools
  can be used as workaround

[cogeqc](/packages/cogeqc)
------

                        Changes in version 1.3.1                        

CHANGES

- Added functions to explore assembly and annotation statistics in a
context: assembly and annotation stats for NCBI genomes can be
extracted through the Datasets API and compared with user-defined
values. New functions: get_genome_stats(), compare_genome_stats(),
and plot_genome_stats().

[cola](/packages/cola)
----

                        Changes in version 2.5.4                        

- adapt the code with markdown v1.6

                        Changes in version 2.5.3                        

- names are automatically added if `cola_opt$color_set_1` and
  `cola_opt$color_set_2`
  are specific without names.

                        Changes in version 2.5.1                        

- In the rmarkdown, replace `message = FALSE` to `message = NA`.

- use `%dorng%` for paralell computing

[compcodeR](/packages/compcodeR)
---------

                       Changes in version 1.35.1                        

- Add a function for SVA + limma differential expression analysis

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

                       Changes in version 2.15.3                        

- `Legend()`: `legend_gp` also controls line color, width and style.

- `anno_mark()`: labels can be duplicated.

                       Changes in version 2.15.1                        

- `Legend()`: allows `NA` in `pch`.

- `SingleAnnotation()`: correctly calculate the max width/height of a
  vector of texts.

- `to_unit()`: fixed a bug when the unit is negative.

- `Legend()`: add `tick_length` argument.

- `Legend()`: colors are correctly calculated when differences between
  `at` are not equal.

[CompoundDb](/packages/CompoundDb)
----------

                         Changes in version 1.3                         

Changes in version 1.3.3

- Add backendBpparam to define (disable) parallel processing for the
MsBackendCompDb backend.

Changes in version 1.3.2

- Evaluate validity of the MsBackendCompDb using the full test suite
from the Spectra package.

Changes in version 1.3.2

- Add parameter nonStop to compound_tbl_sdf that is passed to
parameter skipErrors of ChemmineR::read.SDFset. Issue #110

Changes in version 1.3.0

- Bioconductor 3.17 developmental version.

[concordexR](/packages/concordexR)
----------

                       Changes in version 0.99.0                        

- Submitted to Biconductor.

[CoSIA](/packages/CoSIA)
-----

                       Changes in version 0.99.0                        

- CoSIA is an R package that provides an alternative framework for
cross-species transcriptomic comparison of non-diseased wild-type RNA
sequencing gene expression data across tissues and species through
visualization of variability, diversity, and specificity metrics.
Check
out the Vignette and Readme files to get more information on how to
load
and use CoSIA.

[COTAN](/packages/COTAN)
-----

                       Changes in version 1.99.3                        

Updated the vignette, README.md and NEWS.md

                       Changes in version 1.99.2                        

Dropped second vignette: will be merged in the other one...

                       Changes in version 1.99.1                        

Minor bug fixes and new function clustersMarkersHeatmapPlot()

                       Changes in version 1.99.0                        

Included new functionalities for Bioc 2.17 release:

- created a new COTAN class to replace the old scCOTAN: this class
provides internal invariants along a wide host of accessors that
allows users to avoid peeking inside the class

- made a new multi-core implementation of the model parameters
estimations and COEX calculations that achieves much higher speeds.

- added new functionality about gene clusters starting from given
markers lists

- added new functionality about uniform cell clustering based on the
maximum GDI level in the cluster

- added function to get a differential expression estimation for each
cluster against background

- added function to get an enrichment score for each cluster given a
list of markers specific for the cells' population

- added plots to asses data-set information at cleaning stage

[crisprBase](/packages/crisprBase)
----------

                        Changes in version 1.3.3                        

- Added BiocStyle to Suggests.

                        Changes in version 1.3.2                        

- Added nuclease MAD.

[crisprBowtie](/packages/crisprBowtie)
------------

                        Changes in version 1.3.2                        

- Chromosomes not found in the BSgenome object (usually extra
  chromosomes) are ignored.

[crisprDesign](/packages/crisprDesign)
------------

                       Changes in version 1.1.25                        

- Added a warning message when having duplicated gRNAs.

                       Changes in version 1.1.21                        

- Fixed intergenic annotation of alignments when non-standard
  chromosomes are found.

                       Changes in version 1.1.20                        

- flattenGuideSet is now GuideSet2DataFrames.

                       Changes in version 1.1.19                        

- Changed argument standard_chr_only to FALSE by default in the
  functions getSpacerAlignments and friends.

                       Changes in version 1.1.18                        

- Updated the object guideSetExampleWithAlignments to contain the
  latest annotations.

                       Changes in version 1.1.17                        

- Updated the object guideSetExampleFullAnnotation to contain the
  latest annotations.

                       Changes in version 1.1.12                        

- Refactored addGeneAnnotation.

                       Changes in version 1.1.10                        

- Added the function getPreMrnaSequence.

                        Changes in version 1.1.8                        

- Added features for non-targeting controls (addNtcs).

[crisprScore](/packages/crisprScore)
-----------

                        Changes in version 1.3.4                        

- Fixed unit testing for MIT scores.

                        Changes in version 1.3.3                        

- One more fix to the MIT scoring. Forgot to save sysdata last
  time.

                        Changes in version 1.3.1                        

- Fixed MIT formula. Previous calculations were erroneous.

[csdR](/packages/csdR)
----

                        Changes in version 1.5.1                        

- Fixed segfault issue which did occur in partial_argsort() when the
n_elements argument was larger than the length of the input vector.
In order to ensure equivalence with order(x, decreasing =
TRUE)&#91;1:n_elements&#93;, the additional elements, if any, are padded at
the end of the answer as NA values.

[CTdata](/packages/CTdata)
------

                        Changes in version 0.99                         

CTdata 0.99.0

- Package submission to Biocoductor.

[cytofQC](/packages/cytofQC)
-------

                 Changes in version 0.99.3 (2022-11-11)                 

- Resubmitted to Bioconductor
  * Added script to inst/script that demonstrates how extdata was
  created
  * Deleted doc folder and contents
  * Added a man level page for cytofQC
  * Added relavence to Bioconductor on vignette
  * Added functions tech, scores, probs, initial, and label to replace
  previous get functions that performed the same function to be more
  consistent with Bioconductor syntax
  * Added match.arg and checks for datatypes to all relavent functions
  * Removed plotInitialGuess2 and added argument to plotInitialGuess to
  choose if one or both graphs are plotted
  * Cleaned up warning and message statements

                 Changes in version 0.99.0 (2022-09-13)                 

- Submitted to Bioconductor

[cytomapper](/packages/cytomapper)
----------

                 Changes in version 1.11.3 (2023-01-23)                 

- loadImages: added option to read in single-channel images to
  multi-channel

                 Changes in version 1.11.2 (2023-02-02)                 

- Bug fix: measureObjects internally sets the correct channel names

                 Changes in version 1.11.1 (2023-01-18)                 

- Bug fix: set default dodge.width = NULL for geom_quasirandom

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

                        Changes in version 0.99                         

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

[cytoviewer](/packages/cytoviewer)
----------

                 Changes in version 0.99.3 (2023-04-18)                 

- add initial color control for cell-level

                 Changes in version 0.99.2 (2023-04-13)                 

- image download fix

                 Changes in version 0.99.1 (2023-04-12)                 

- code adjustments after Bioconductor approval

- included reactive image / mask reading

                 Changes in version 0.99.0 (2023-03-23)                 

- code preparations for Bioconductor submission

[dearseq](/packages/dearseq)
-------

                 Changes in version 1.11.1 (2023-03-16)                 

- add parallel support on Windows

[deconvR](/packages/deconvR)
-------

                        Changes in version 1.5.2                        

- Fixed a bug in deconvolute
- findSignatures now can construct tissue specific methylation
signatures
- findSignatures now can construct tissue specific DMPs

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

                       Changes in version 0.26.0                        

- No changes in this version.

[DELocal](/packages/DELocal)
-------

                       Changes in version 0.99.1                        

- Added a NEWS.md file to track changes to the package.

[DepInfeR](/packages/DepInfeR)
--------

                 Changes in version 1.3.1 (2023-04-04)                  

- Add BiocStyle dependence in the DESCRIPTION file to avoid vingette
  build error on Linux

[DESeq2](/packages/DESeq2)
------

                       Changes in version 1.39.8                        

- Changed `lower=0` to `lower=1e-6` in unmix(), as the
  lower bound of 0 was producing sqrt(negative) errors
  on Linux ARM64.
  https://support.bioconductor.org/p/9150056/

                       Changes in version 1.39.7                        

- Fix bug in independent filtering: with very little variation in
  the curve of number of rejections over threshold, and when the
  maximum was only reached near the end, the default filtering
  wouldn't attain sufficient filtering. This has been addressed
  by also checking for a threshold at which 90%, or 80% of the
  fitted number of rejections is found.
  Note: IHW is the preferred method for filtering, and can easily
  by used by calling `filterFun=ihw`.

                       Changes in version 1.39.6                        

- Fix bug on estimateDispersionsGeneEst when niter is
  larger than 1 (#64 on GitHub).

                       Changes in version 1.39.5                        

- PR from Hendrik Weisser for lfcShrink when the results
  table has additional columns than those produced by
  results().

                       Changes in version 1.39.4                        

- Removing geneplotter dependency.

                       Changes in version 1.39.1                        

- Removing genefilter as dependency, switching to matrixStats.
  This should resolve gfortran issues.

[DESpace](/packages/DESpace)
-------

                        Changes in version 1.0.0                        

- 1st version of the package

- submitted to Bioconductor

[DiffBind](/packages/DiffBind)
--------

                        Changes in version 3.9.4                        

- Rollover bugfixes

- add option to return ggplot in dba.plotVolcano

[DNAfusion](/packages/DNAfusion)
---------

                        Changes in version 1.3.0                        

- Added introns_ALK_EML4()
- Added find_variants()
- Updated break_position_depth()
- Updated break_position()
- Refined EML4_ALK_detection()
- Updated vignette to include new functions
- Cleaned functions for Bioconductor 3.17 release

                        Changes in version 1.1.1                        

- Updated example files to include ALK reads

                        Changes in version 1.1.0                        

- Updated DESCRIPTION

[doubletrouble](/packages/doubletrouble)
-------------

                       Changes in version 0.99.3                        

BUG FIXES

- Updated functions (e.g., get_anchor_list(), collinearity2blocks())
after update in syntenet.

                       Changes in version 0.99.2                        

CHANGES

- Small change in coding style after Bioconductor peer-review (m:n
replaced with c(m, n) and seq(m,n))

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[edgeR](/packages/edgeR)
-----

                       Changes in version 3.42.0                        

- 
  New function Seurat2PB() for creating a pseudo-bulk DGEList
  object from a Seurat object. New case study in User's Guide
  illustrating its use.

- 
  New function normLibSizes() is now a synonym for
  calcNormFactors().

- 
  Rename effectiveLibSizes() to getNormLibSizes().

- 
  DGEList() is now an S3 generic function with a method for
  data.frames. The data.frame method allows users to specify
  which columns contain gene annotation and which contain counts.
  If the annotation columns are not specified, the function will
  check for non-numeric columns and will attempt to set the
  leading columns up to the last non-numeric column as
  annotation. 'y' is now a compulsory argument for DGEList().
  Previously it defaulted to a matrix with zero rows and zero
  columns.

- 
  New case study in User's Guide on a transcript-level different
  expression analysis.

- 
  The case study on alernative splicing in the User's Guide has
  been replaced with a new data example.

[enrichplot](/packages/enrichplot)
----------

                       Changes in version 1.19.2                        

- fix emapplot() for parameter mismatch (2023-02-20, Mon)
- fix ridgeplot for error when x@readable == TRUE and
length(x@gene2Symbol) = 0 (2022-12-5, Mon)
- fix ridgeplot for error when x@readable == TRUE and
length(x@gene2Symbol) = 0 (2022-12-5, Mon, #217)

                       Changes in version 1.19.1                        

- fix cnetplot() for node_label parameter is flipped(2022-12-04, Sun,
#216)
- bug fixed in treeplot() (2022-11-18, Fri)
- enable dotplot() and autofacet() for gseaResultList object

[ensembldb](/packages/ensembldb)
---------

                       Changes in version 2.23.2                        

- Fix SQLite database names to include also the subspecies.

                       Changes in version 2.23.1                        

- Remove disjointExons method.

[EpiCompare](/packages/EpiCompare)
----------

                        Changes in version 1.3.4                        

New features

- Example report:
- Delete report/ folder and upload to Releases instead:
https://github.com/neurogenomics/EpiCompare/releases
- Add Rscript to replicate example report in inst/examples.
- EpiCompare
- New arg add_download_button.
- Always keep download button for post-processed peak files.

Bug fixes

- README.Rmd:
- Fix broken link to example report.

                        Changes in version 1.3.3                        

New features

- download_button:
- Saves and downloads files.
- prepare_blacklist:
- Auto-selects appropriate blacklist, or returns user-specified
option.
- EpiCompare(blacklist=NULL) is now the default.
- prepare_genome_builds:
- Update to handle supplying builds for "peakfiles" and
"reference" but not "blacklist" (so long as the blacklist arg is
not a user-supplied GRanges object)
- Added mm9_blacklist
- Made more plots interactive:
- width_boxplot
- plot_enrichment
- plot_ChIPseeker_annotation
- overlap_stat_plot
- Name elements in output list.
- Change annotation arg to more informative txdb arg, and set default
to NULL, which ChIPseeker functions will automatically handle.
- New function as_interactive:
- Help standardise this.
- New EpiCompare::EpiCompare arguments:
- error: keep knitting even on errors.
- tss_distance: upstream/downstream of TSS.
- quiet: knit quietly
- Rename 'test-EpiCompare_combinations.R' --> 'test-EpiCompare.R'
- Separate test-generalMetrics_functions.R into function-specific
test
files.
- Separate test-peakOverlap_functions.R into function-specific test
files.
- Make fancy header with new func:
- report_header()
- Create EpiCompare command code as text:
- report_command()
- width_boxplot:
- Make more efficient with data.table and lapply
- Update hex sticker to match custom.css palette.
- README.Rmd
- Collapse more detailed sections.

Bug fixes

- tss_plot:
- Fix examples/tests after Sera updated the arguments.
- Pass upstream/downstream to ChIPseeker::getTagMatrix
- Make interactive
- Name plots in list
- Remove unnecessary extra level of list nesting.
- Make documentation width <80 lines where possible.
- EpiCompare.Rmd
- Remove methods::show from all parts
- Name all chunks
- Make explanations more clear
- Add table of contents for main 3 sections.
- Fix header levels
- Set results='asis' globally instead of in each chunk header.
- Automatically number sections with yaml arg: number_sections:
true
- Omit specific headers from numbering system with {-} tags.
- Add custom.css
- plot_chromHMM:
- Error in (function (classes, fdef, mtable) unable to find an
inherited method for function ‘annotateWithFeatures’ for
signature ‘"SimpleGRangesList", "list"’
- Misleading error message; was actually due to
chromHMM_annotation not being converted from a list to a
GRangesList.
- Change yaml arg peakfile --> peakfiles to be consistent with
other variables.

                        Changes in version 1.3.1                        

New features

- Replace badger with rworkfows:
- Use rworkflows::use_badges
- New helper functions:
- precision_recall_matrix
- report_time
- overlap_upset_plot:
- Switched out UpSetR for ComplexUpsetto show percentages.
- Moved up dep checks to beginning of function.
- Handle bug with heatmaply by checking args where it might be used:
- check_heatmap_args
- tss_plot:
- Add unit tests
- Drastically reduce example/test runtime by setting upstream=50
- compute_corr:
- Reduce example runtime by setting bin_size = 200000 (takes
<2s).

Bug fixes

- Fix typo in EpiCompare docs: "hg38 blacklist dataset"
- Avoid explicitly specifying "/" in paths to help cross-platform
testing.
- tss_plot:
- Use parallel::detectCores-1 by default to set workers, but set
to 1 in examples/tests to meet CRAN/Bioc standards.

[EpiDISH](/packages/EpiDISH)
-------

                       Changes in version 2.15.0                        

- Add 12 blood cell-type DNAm reference matrix for 450k and EPIC
  arrays.

[epigraHMM](/packages/epigraHMM)
---------

                        Changes in version 1.7.2                        

- Fix callPeaks writers so that it works on Windows

                        Changes in version 1.7.1                        

- Fix callPeaks function to export .bed, .wig, and .bedGraph files
- Output files are now overwritten, if existing

[epistack](/packages/epistack)
--------

                 Changes in version 1.5.4 (2022-04-05)                  

- CITATION added

                 Changes in version 1.5.3 (2022-04-04)                  

- plotAverageProfile's reversed_z_order is now exposed in
plotEpistack
- 95% confidence interval is now the default in plotEpistack
- plotEpistck's legends are now also displayed as y-axis title in
plotAverageProfile plots

                 Changes in version 1.5.2 (2022-04-03)                  

- changing the default bin colors
- bin_palette in plotEpistack(), plotBoxMetric(), and
plotAverageProfile(), and palette in plotBinning() can now be
vectors of colors instead of palette functions
- tints now accept palette functions and list of palette functions,
in
addition to colors.

                 Changes in version 1.5.1 (2022-03-30)                  

- vignette improvements (thanks to Isabelle Stevant)
- documentation improvments: multiple zlims are possible if provided
as a list

[erccdashboard](/packages/erccdashboard)
-------------

                       Changes in version 1.33.1                        

SIGNIFICANT USER-VISIBLE CHANGES / BUG FIX

- Vignette output showed incorrect mRNA fraction estimates and
  incorrect LODR estimates due to the stringsAsFactors=FALSE switch
  with R > 4.0. Ratio values were not being treated as factors in the
  loadERCCInfo function and this was causing incorrect ratio value
  assignments in the est_r_m.R function. The issue has been resolved by
  an explicit as.factor assignment in the loadERCCInfo.R function.

[escheR](/packages/escheR)
------

                       Changes in version 0.99.8                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add a new argument y_reverse = TRUE to make_escheR to provide a
consistent orientation between spot plot and tissue image (see Issue
#13)

                       Changes in version 0.99.7                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add default color scheme (viridis) to add_fill
- Add explicit reference to spatialLIBD in make_escheR documentation
- Add installation instruction for users whose R version is pre-R4.3

                       Changes in version 0.99.6                        

SIGNIFICANT USER-VISIBLE CHANGES

- Add minimium versions to dependencies and imported packages
- Import individual functions in NAMESPACE from packages
- Clean up comments in code
- Accepted by Bioconductor and will be released in Bioconductor 3.17

                       Changes in version 0.99.1                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.
- First full version of the package to be submitted to Bioconductor.
See Bioconductor submission here.

[EWCE](/packages/EWCE)
----

                        Changes in version 1.7.3                        

New features

- check_ewce_genelist_inputs/bootstrap_enrichment_test
- New arg: sctSpecies_origin lets users clarify that their data
originally came from mouse even when it is currently formatted
as human orthologs. This is necessary for creating the
appropriate background gene lists.
- Remove grDevices as dep entirely.
- fix_celltype_names
- Add new arg make_unique to make this function easily usable for
vectors where the same celltype appears multiple times.
- bootstrap_enrichment_test
- Return gene-level scores based on adaptation of code from
generate_bootstrap_plots. now stored as a list element named
gene_data in data.table format.
- generate_bootstrap_plots
- Revamp wrap code into reusable subfunctions.
- Avoid resampling random genes when the gene data is stored in
the bootstrap results as gene_data. It will also tell you which
of these options it's using.
- Save with ggsave instead if grDevices.
- Facet by celltype instead of generating tons of separate plots.
- Let users decide cutoff threshold with new arg adj_pval_thresh
- Now returns a named list with the plots themselves ("plot") and
the paths to where they're saved ("paths") rather than just a
higher-level directory path in which users had to search for the
right files (and didn't ever have access to the ggplot2
objects).
- Show significance with barplot fill/color instead of asterices.
Much easier to see now.
- Change savePath arg to the more accurate save_dir. Expose
appending BootstrapPlots to the user within the argument.
- generate_bootstrap_plots_for_transcriptome
- Change savePath arg to the more accurate save_dir. Expose
appending BootstrapPlots to the user within the argument.
- Save with ggsave instead if grDevices.
- Standardise hits + hitGenes arg all to hits.
- Update hex:
- Off load large source image from DALLE to Releases instead of
including it within the package.

Bug fixes

- drop_uninformative_genes / generate_celltype_data
- Pass verbose arg to matrix formatting functions.
- generate_controlled_bootstrap_geneset
- Removed combinedGenes arg as it was not being used anywhere
within.
- check_args_for_bootstrap_plot_generation
- Removed unused args: ttSpecies, sctSpecies
- test-bootstrap_enrichment_test_2.R
- "monkey_ctd" tests seems to be running more smoothly than before
(not just getting NAs). This might have to with orthogene
databases improving.
- Reassuringly, "godzilla" tests still fail as expected :)
- Add tess/testthat/Rplots.pdf to .gitignore.

                        Changes in version 1.7.1                        

New features

- Use rworkflows GHA.
- Add rworkflows::use_badges to README.Rmd.
- Remove Dockerfile (no longer necessary).
- Make all 3 platforms (Linux, Mac, Windows) use Bioc dev, as
ewceData (>=1.7.1) is now required, due to a fix made only in
the development version of rtracklayer.
- Remove cowplot dependency.
- Replace all %>% with |>

Bug fixes

- calc_quantiles:
- This function was only used in filter_variance_quantiles
- Compare stats::ecdf vs. dplyr::ntile methods.
- Remove from EWCE as it's not longer used anywhere.
- bin_columns_into_quantiles:
- Rename arg matrixIn --> vec to reflect what the function
actually does.
- filter_variance_quantiles:
- Change to use bin_columns_into_quantiles instead of
calc_quantiles to be consistent with how quantiles are handled
in the rest of EWCE.
- Updated tests in test-get_celltype_table.r to reflect that the
number of genes filtered is unaffected by the normalization
procedure (when quantiles are computed with stats::quantile).
- ewce_plot:
- Celltypes were producing NAs because the names in the results
were not always standardized in the same way as the CTD. Now
this is done internally.
- Celltypes were not ordered factors, meaning the dendrogram
didn't line up correctly.
- Switched from gridArrange/cowplot to patchwork.
- Added dedicated unit tests file: test-ewce_plot.r

New features

- Offline runs enabled with functions using reference datasets (from
ewceData). These functions have the parameter localhub added to
control this.

[ExperimentHub](/packages/ExperimentHub)
-------------

                        Changes in version 2.7.1                        

DEPRECATE

- (2.7.1) Deprecate display in favor of BiocHubsShiny

[extraChIPs](/packages/extraChIPs)
----------

                        Changes in version 1.3.9                        

- Changed labelling strategy for plotPie()

                        Changes in version 1.3.8                        

- Added using coverage to makeConsensus()

                        Changes in version 1.3.7                        

- Added plotAssayHeatmap()
- Added fitAssayDiff() and added coercion of TopTags objects

                        Changes in version 1.3.6                        

- Expanded arguments for plotSpliDonut()
- Added mergeByHMP() for merging overlapping windows using the
harmonic mean p

                        Changes in version 1.3.5                        

- Added plotSplitDonut()
- Bugfix in plotAssayDensities() and plotAssayRle() along with
enabling plotting by group

                        Changes in version 1.3.4                        

- Added mergeBySig()


[FeatSeekR](/packages/FeatSeekR)
---------

                       Changes in version 0.99.1                        

- Submitted to Bioconductor

[fgsea](/packages/fgsea)
-----

                       Changes in version 1.25.1                        

- introduced plotEnrichmentData() function for more flexible plotting

- update GESECA plot behavior

[fishpond](/packages/fishpond)
--------

                        Changes in version 2.5.4                        

- As CellRanger 7 includes both spliced and unspliced counts in their
count matrix, we want to mimic this behavior by adding more
pre-defined output formats in the loadFry function. We added "all"
and "U+S+A" to include all counts in the count matrix. Moreover, now
the "scRNA" output format has an "unspliced" field, which contains
the unspliced count matrix.

                        Changes in version 2.5.3                        

- Fix bug where salmonEC did not correct equivalence class names for
going from 0-indexing to 1-indexing internally. In prior
versions, to correctly link equivalence classes to gene names, users
would have needed to manually add a value of 1 to the equivalence
class names, which was erroneously not mentioned in the man files.
After this bug fix, if the equivalence class identifier reads 1|2|8,
then the equivalence class is immediatly compatible with the
transcripts and their respective genes in rows 1, 2 and 8 of
'tx2gene_matched', without any further user intervention.

                        Changes in version 2.5.1                        

- Fix plotAllelicGene() so that when samples have an allele with no
expression at the gene level, it doesn't throw an error trying to
divide by 0.

[flowGate](/packages/flowGate)
--------

                 Changes in version 0.99.3 (2023-03-15)                 

- Document cleanup as part of bioconductor review

[FlowSOM](/packages/FlowSOM)
-------

                        Changes in version 2.7.9                        

- Adapted documentation of readInput and FlowSOM. It can also use a
  matrix with
  column names

                        Changes in version 2.7.8                        

- Added PlotOutliers and documentation

                        Changes in version 2.7.7                        

- Updates to TestOutliers

                        Changes in version 2.7.5                        

- Edited PlotManualBars so that it shows percentages of cells

- Edited PlotDimRed so that it optionally uses scattermore if it
  installed

                        Changes in version 2.7.4                        

- AggregateFlowFrames can now also resample if more cells are asked for
  than
  available in the fcs files. Default stays FALSE, taking at most the
  number
  of cells in the fcs file

                        Changes in version 2.7.2                        

- Version bump to align with Bioconductor

- Update for aggregateFlowFrames for nices handling of iterative
  aggregation
  (+ according visualisation in FlowSOMmary). AggregateFlowFrames will
  now introduce
  0 values if a channel is not present in one of the fcs files. If
  channels is not
  provided as an argument, the channels of the first file are used.

- Added support for abbrevations in PlotDimRed "colorBy".

- Return 0 instead of NA in case of an empty cluster in
  GetCounts/GetPercentages

[flowSpecs](/packages/flowSpecs)
---------

                 Changes in version 1.13.1 (2023-04-05)                 

- +Updates to specMatCalc to work with BigFoot data.

[fmrs](/packages/fmrs)
----

                        Changes in version 2.0.0                        

IMPROVEMENTS SINCE LAST RELEASE

- The package is rewritten using .Call function.
- The codes for Weibull distribution are improved.

BUG FIXES

- Several bugs are fixed which caused the results to be different for
the same analysis.

[gdsfmt](/packages/gdsfmt)
------

                       Changes in version 1.36.0                        

UTILITIES

- fix the compiler warning: sprintf is deprecated

- LZ4 updated to v1.9.4

- XZ updated to v5.2.9

- update zlib to v1.2.13

NEW FEATURES

- `system.gds()$compiler.flag[1]` is either "64-bit" or "32-bit"
  indicating
  the number of bits of internal data pointer

- new argument 'use.abspath=TRUE' in `openfn.gds()` and
  `createfn.gds()`:
  the behavior before v1.35.4 is the same as 'use.abspath=TRUE'

                       Changes in version 1.34.1                        

UTILITIES

- avoid using `crayon::blurred()` in the display (RStudio blurs the
  screen
  output)

[genbankr](/packages/genbankr)
--------

                       Changes in version 1.27.1                        

BUGFIXES

- parser now correctly reads the "D-loop" feature type (#22)

[GeneTonic](/packages/GeneTonic)
---------

                        Changes in version 2.4.0                        

New features

- enhance_table() has now the possibility to plot the visual
summaries
as ridge lines.
- When plotting the gene expression for the selected features in the
gene-geneset-graph box, it is now possible to disable the labels
from being displayed (could lead to unnecessary clutter sometimes).

Other notes

- Fortified the behavior of gene_plot() to fail early when providing
an invalid value to the intgroup parameter.

[GenomAutomorphism](/packages/GenomAutomorphism)
-----------------

                        Changes in version 1.0.3                        

- Add new functions to work with modular matrix operations

                        Changes in version 1.0.2                        

- Fix error of parallel computation on Windows (12/08/2022)

                        Changes in version 1.0.1                        

- Expanding analyses by including aminoacid similarity based on codon
distances. Three new functions are added: codon_dist,
codon_dist_matrix, and aminoacid_dist. See a tutorial applying these
functions at https://is.gd/oYLDK4.

[GenomeInfoDb](/packages/GenomeInfoDb)
------------

                       Changes in version 1.36.0                        

NEW FEATURES

- Register the following NCBI assemblies:
  - the bCatUst1.alt.v2 and bCatUst1.pri.v2 assemblies
  - the felCat9.1_X and F.catus_Fca126_mat1.0 assemblies
  - the Gossypium_hirsutum_v2.1 assembly
  - the Kamilah_GGO_v0 assembly
  - the Xenopus_laevis_v2 and Xenopus_laevis_v10.1 assemblies

- Register the following UCSC genomes:
  - gorGor6 (linked to Kamilah_GGO_v0)
  - xenLae2 (linked to Xenopus_laevis_v2)

- Add Gossypium_hirsutum.txt to inst/extdata/dataFiles/ (provided
  by Emory Lucas <bararayung123@hotmail.com>)

- Add 'organism' argument to registered_UCSC_genomes() (contributed by
  Kirabo Kakopo).

SIGNIFICANT USER-VISIBLE CHANGES

- UCSC genome hg38 is now based on GRCh38.p14 instead of GRCh38.p13
  (this
  change originated at UCSC). See commit 091b5d2.

- The submitters for NCBI assembly Dog10K_Boxer_Tasha have updated
  the info for the MT sequence, which is reflected in the output of
  getChromInfoFromNCBI("Dog10K_Boxer_Tasha"). See commit 79a066c.

- The Accept-organism-for-GenomeInfoDb vignette was converted from
  Rnw to Rmd (thanks to Haleema Khan and Jen Wokaty for this
  conversion).

- Small improvements to low-level helper find_NCBI_assembly_ftp_dir().

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

                       Changes in version 1.36.0                        

NEW FEATURES

- Add 'strandMode' argument to readGAlignmentsList() (contributed by
  Robert Castelo).

BUG FIXES

- Increase 'cigar_buf' size to reduce risk of buffer overflow in
  cigar-utils C code.

[GenomicDataCommons](/packages/GenomicDataCommons)
------------------

                       Changes in version 1.24.0                        

Bug fixes and minor improvements

- gdc_clinical handles NULL responses when diagnoses are not
available
for all IDs queried (#109, @zx8754).
- Minor updates to somatic mutations vignette and unit tests.

[GenomicFeatures](/packages/GenomicFeatures)
---------------

                       Changes in version 1.52.0                        

NEW FEATURES

- Small improvement to makeTxDbFromEnsembl(): The function can now be
  called on the abbreviated organism, e.g. "hsapiens", in addition
  to "homo sapiens".

DEPRECATED AND DEFUNCT

- Finally remove disjointExons() (got deprecated in BioC 3.13 and
  defunct
  in BioC 3.15).

BUG FIXES

- Fix issue with order of sequences in seqinfo(makeTxDbFromUCSC()).
  See commit e4381bc.

[GenomicRanges](/packages/GenomicRanges)
-------------

                       Changes in version 1.52.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Improve documentation of findOverlaps() argument 'minoverlap':
  The man page now explains how findOverlaps() argument 'minoverlap'
  is interpreted when 'query' or 'subject' is a GRangesList object.

[GeoTcgaData](/packages/GeoTcgaData)
-----------

                       Changes in version 0.99.2                        

- use SummarizedExperiment input (2023-1-29, Sun)
- fix return value of differential_array (2022-10-8, Sat)
- fix gene length bug in countToTpm() and countToFpkm()(2022_9_22,
Tue)
- fix a bug in id_conversion (2022-8-27, Sat)
- fix a bug in differential_RNA(useTopconfects = TRUE) (2022-8-12,
Fir)
- add function methydifferential_ucsc and
methydifferential_limma(2021-10-24, Sun)
- update hgnc_file data(2021-10-24, Sun)
- add function differential_RNA to do difference analysis of RNA-seq
data(2021-7-20, Tue)
- add data hgnc_file
- update function: id_ava()
- add functions: ann_merge(), countToFpkm(), countToTpm()

[ggtree](/packages/ggtree)
------

                        Changes in version 3.7.2                        

- ggtree() supports 'dendro' object (ggdendro::dendro_data() output)
(2023-03-02, Thu)
- update theme_dendrogram() to use ggfun::theme_noxaxis()
(2022-11-21,
Mon)
- using cli::cli_alert_warning() instead of warning_wrap (2022-11-10,
Thu)

                        Changes in version 3.7.1                        

- compatible with ggplot2 v=3.4.0 (2022-11-07, Mon)
- allows setting options(clade_align = TRUE) to align geom_hilight()
layer and allows setting options(clade_width_extend = 0.35) to set
the amount the width extension (in y-axis) of geom_hilight(). These
two features is designed for ggtreeDendro::geom_rect_subtree() layer
(2022-11-06, Sun)

[ggtreeDendro](/packages/ggtreeDendro)
------------

                        Changes in version 1.1.3                        

- autoplot method for 'genoMatriXeR' and 'multiLocalZScore' objects
(both defined in the 'regioneReloaded' package) (2023-02-11, Sat)
- autoplot method for 'ClusterExperiment' object (2022-11-28, Mon)

                        Changes in version 1.1.2                        

- plot_wgcna() for visualizing WGCNA dendrogram (2022-11-06, Sun)
- supports hkmeans object defined in the factoextra package
(2022-11-06, Sun)
- set default options to align geom_rect_subtree() and adjust width
extension.

                        Changes in version 1.1.1                        

- geom_rect_subtree layer to hilight subtrees (2022-11-03, Thu)
- supports bclust object (S3) output by e1071::bclust()
- supports hdbscan object

[ggtreeExtra](/packages/ggtreeExtra)
-----------

                        Changes in version 1.9.2                        

- add limits parameter in axis.params to display the range of x axis.
(2023-03-08, Wen)
- using cli to replace warning or stop. (2022-11-24, Thu)
- using rlang to transfer the geom argument. (2022-11-24, Thu)
- using pwidth to replace width when x only have one unique value
with
geom = geom_tile and width is not provided. (2022-11-25, Fri)
-
https://github.com/YuLab-SMU/ggtreeExtra/issues/11#issuecomment-873335702
- https://github.com/YuLab-SMU/ggtreeExtra/issues/12
- https://github.com/YuLab-SMU/ggtreeExtra/issues/25
- adding the messages to show how to adjust the color scale when the
color aesthetic was use internal with geom=geom_boxplot etc.
(2022-12-01, Thu)
- https://github.com/YuLab-SMU/ggtree/issues/551

                        Changes in version 1.9.1                        

- import cli to avoid warning in checking, which is used by the
latest
ggplot2 (2022-11-08, Tue)

[glmGamPoi](/packages/glmGamPoi)
---------

                  Changes in version 1.11 (2023-01-03)                  

- Breaking change: rename 'pseudobulk_sce' to 'pseudobulk'

- Add a new vignette explaining how and why pseudobulking is
  a powerful concept for single cell data analysis

- Depcreate 'pseudobulk_by' argument in 'test_de'. Use the 'pseudobulk'
  function
  instead.

- Add a new argument 'max_lfc' to to 'test_de' to avoid impractically
  large
  log fold changes for lowly expressed genes.

- Support rlang quosures for the contrast argument in 'test_de'

- Add a helper function called 'fact' that simplifies specification
  of contrast for complex experimental designs

- Add 'use_assay' argument to 'glm_gp'

- Add `vctrs` as dependency. The package is necessary to replicate the
  `group_by` behavior from `dplyr`.

- Add 'size_factors = "ratio"' to emulate the behavior of DESeq2's size
  factor calculation

- Make sure that the 'ignore_degeneracy' argument is propagated to
  'test_de'

[GRaNIE](/packages/GRaNIE)
------

                 Changes in version 1.3.34 (2023-03-29)                 

Bugfixes

- many small bugfixes and other small improvements to homogenize the
user experience due to the usage of systematic unit tests

                 Changes in version 1.3.33 (2023-03-06)                 

New features and vignette updates

- we provide two new functions with this update:

1.  getGRNSummary() that summarizes a GRN object and returns a named
list, which can be used to compare different GRN objects ore
easily among each other, for example.
2.  plotCorrelations() for scatter plots of the underlying data for
either TF-peak, peak-gene or TF-gene pairs. This can be useful
to visualize specific TF-peak, peak-gene or TF-gene pairs to
investigate the underlying data and to judge the reasonability
of the inferred connection.

- methods vignette updates

             Changes in version 1.3.31-1.3.32 (2023-03-06)              

Bugfixes

- various small bugfixes that were accidentally introduced in the
latest change from using the TF.ID instead of TF.name column as
unique TF identifier

             Changes in version 1.3.26-1.3.30 (2023-03-06)              

New features and vignette updates

- added two new supported genomes: rn6/rn7 and dm6 for the rat and
the
Drosophila (fruit fly) genome, respectively
- added preliminary support for a new, alternative way of how to
import TF and TFBS data into GRaNIE. We now additionally offer a
more user-friendly way by making it possible to directly use the
JASPAR2022 database. You do not need any custom files anymore for
this approach! See the Package vignette for more details.

Bugfixes

- fixed a regression bug in addConnections_TF_peak (Column
peak.GC.class doesn't exist.) that was caused due to the recent GC
modifications

                 Changes in version 1.3.25 (2023-02-20)                 

New features and vignette updates

- additional significant methods vignette updates
- updates and clarifications for the workflow vignette
- a new QC plot for plotDiagnosticPlots_TFPeaks (and indirectly in
addConnections_TF_peak when plotDiagnosticPlots = TRUE) on page 1
that shows the total number of connections for real and background
TF-peak links as calculated and stored in the GRN object, stratified
by TF-peak FDR and correlation bin. This is a similar plot as we
show in the paper and helps comparing foreground and background.

Improvements

- speed improvements for plotDiagnosticPlots_TFPeaks (and indirectly
in addConnections_TF_peak when plotDiagnosticPlots = TRUE) when
plotAsPDF = FALSE

Bugfixes

- fixed a bug that only occurred in addConnections_TF_peak when using
useGCCorrection = TRUE

                 Changes in version 1.3.24 (2023-02-20)                 

New features and vignette updates

- significant methods vignette updates that help clarifying methods
details

             Changes in version 1.3.22-1.3.23 (2023-02-15)              

Minor changes

- Small workflow vignette updates

Bug fixes

- we were informed that newer versions of dplyr (1.1.0) changed their
default behavior for the function if_else when NULL is involved,
which caused an error. We changed the implementation to accommodate
for that and now avoid dplyr::if_else and use base R ifelse instead.

             Changes in version 1.3.18-1.3.21 (2023-02-07)              

Minor changes

- Small vignette updates and fixing typos / improved wording

Bug fixes

- due to a change from USCS that affected
GenomeInfoDb::getChromInfoFromUCSC("hg38") (see here for more
details), the minimum required version of GenomeInfoDb had to be
increased to 1.34.8. If you have troubles installing at least this
version, we recommend updating to the newest Bioconductor
version 3.16 or (without warranties) use the following line to
manually install the newest version directly from GitHub outside of
Bioconductor (not recommended):
BiocManager::install("Bioconductor/GenomeInfoDb)"
- small change in addData() so that peak IDs are stored with the same
name in the object in case the user-provided peak IDs have the
format chr:start:end as opposed to the required chr:start-end.
filterData() otherwise incorrectly discarded all peaks because of
the ID mismatch caused by the two different formats.
- fixed a rare edge case in filterGRNAndConnectGenes() that caused an
error when 0 TF-peak connections were found beforehand

                 Changes in version 1.3.17 (2023-01-26)                 

New features

- We are excited to announce that we added a new vignette for how to
use GRaNIE for single-cell data! We plan to update it regularly with
new information. Check it out here!

                 Changes in version 1.3.16 (2023-01-24)                 

New features

- significant updated to the package details vignette
- revisited and improved the internal logging and object history. The
time when a function was called is now added to the list name, which
allows the storage of multiple instances of the same function.
- new parameter in addData(): geneAnnotation_customHost to specify a
custom host and overriding the default and previously hard-coded
hostname when retrieving gene annotation data via biomaRt.
- the function getGRNConnections() can now also include the various
additional metadata for all type parameters and not only the default
type all.filtered.

                 Changes in version 1.3.15 (2023-01-20)                 

Bug fixes

- fixed an error that appeared in rare cases when a chromosome name
from either peak or RNA data could not be found in biomaRt such as
GL000194.1. Peaks from chromosomes with irretrievable lengths are
now automatically discarded.
- significant updates to the package details vignette

             Changes in version 1.3.13-1.3.14 (2023-01-20)              

New features

- the function plotDiagnosticPlots_peakGene() (which is also called
indirectly from addConnections_peak_gene() when setting
plotDiagnosticPlots = TRUE) now stores the plot data for the QC
plots from the first page into the GRN object. It is stored in
GRN@stats$peak_genes
- the columns of the result table from getGRNConnections() are now
explained in detail in the R help, and we reference this from the
Vignette and other places
- various significant Vignette updates

Bug fixes

- optimized the column names for the function getGRNConnections(),
which now does not return duplicate columns for particular cases
anymore
- improved printing in the log for the function filterData() and
addData()
- the loadExampleObject() function has been optimized and should now
force download an example object when requesting it.
- the package version as stored in the GRN object now works
correctly.

Minor changes

- further code cleaning in light of the tidyselect changes in
version 1.2.0 to eliminate deprecated warnings
- the default gene types for addConnections_peak_gene() and
plotDiagnosticPlots_peakGene() have been homogenized and changed to
list(c("all"), c("protein_coding")). Before, the default was
list(c("protein_coding", "lincRNA")), but we decided to now split
this into two separate lists: Once for all genes irrespective of the
gene type and once for only protein-coding genes. As before, lincRNA
or other gene types can of course still be selected and chosen.
- various minor changes

                 Changes in version 1.3.12 (2022-12-22)                 

Bug fixes

- bug fix in plotCommunitiesEnrichment() that was introduced due to
the tidyselect 1.2.0 changes

Minor changes

- further code cleaning in light of the tidyselect changes in
version 1.2.0 to eliminate deprecated warnings

                 Changes in version 1.3.11 (2022-12-16)                 

Major changes

- the default URL for the example GRN object in loadExampleObject()
had to be changed due to changes in the IT infrastructure. The new
stable default URL is now
\url{https://git.embl.de/grp-zaugg/GRaNIE/-/raw/master/data/GRN.rds},
in the same Git repository that provides GRaNIE outside of
Bioconductor.

Bug fixes

- fixing bugs introduced due to the tidyverse 1.2.0 related code
cleaning
- other bugfix accidentally introduced in the previous commits

                 Changes in version 1.3.10 (2022-12-15)                 

Bug fixes

- revisited the import of TADs and made the code more error-prone and
fixed some bugs related to TADs. Importing TADs now works again as
before.

Minor changes

- code cleaning in light of the tidyselect changes in version 1.2.0
to
eliminate deprecated warnings

New features

- new argument for addConnections_peak_gene(): TADs_mergeOverlapping.
See the R help for more details.

                 Changes in version 1.3.9 (2022-12-14)                  

New features

- new argument for addConnections_peak_gene(): shuffleRNACounts. See
the R help for more details.

Minor changes

- first round of code cleaning in light of the tidyselect changes in
version 1.2.0 to eliminate deprecated warnings

              Changes in version 1.3.4-1.3.8 (2022-12-06)               

Major changes

- the topGO package is now required package and not optional anymore.
The reasoning for this is that the standard vignette should run
through with the default arguments, and GO annotation is the default
ontology so topGO is needed for this. Despite this package still
being optional from a strict workflow point of view, we feel this is
a better way and improves user friendliness by not having to install
another package in the middle of the workflow.

Minor changes

- in initializeGRN(), the objectMetadata argument is now checked
whether it contains only atomic elements, and an error is thrown if
this is not the case. As this list is not supposed to contain real
data, checking this prevents the print(GRN) function to
unnecessarily print the whole content of the provided object
metadata, thereby breaking the original purpose.

New features

- addTFBS() got two more arguments to make it more flexible. Now, it
is possible to specify the file name of the translation table to be
used via the argument translationTable, which makes it more flexible
than the previously hard-coded name "translationTable.csv. In
addition, the column separator for this file can now be specified
via the argument translationTable_sep
- Overlapping TFBS data with the peak is now more error-tolerant and
does not error out in case that some chromosome or contig names from
the TFBS BED files contain elements the size of which cannot be
retrieved online. This was the case for some contig names with the
suffix decoy, for example. If such elements are found, a warning is
now thrown and they are ignored as they are usually not wanted
anyway.
- in case a GRN objects contains 0 connections (e..g, because of too
strict filtering), subsequent functions as well as the print
function now give a more user-friendly warning / error message.

[graphite](/packages/graphite)
--------

                 Changes in version 1.45.1 (2023-04-19)                 

- Updated all pathway data.

[hca](/packages/hca)
---

                         Changes in version 1.8                         

New features

- (v. 1.7.1) Implement project_title() and project_information() to
retrieve summary information from project IDs.

[HDF5Array](/packages/HDF5Array)
---------

                       Changes in version 1.28.0                        

- No changes in this version.

[HIBAG](/packages/HIBAG)
-----

                       Changes in version 1.34.1                        

- fix the compiler warning: sprintf is deprecated

- show "64-bit" correctly when run on Windows

[HiCDOC](/packages/HiCDOC)
------

                 Changes in version 1.1.1 (2022-11-30)                  

- Change the computation of selfInteraction (use diagonal only)

- review the distance normalisation (remove bias)

- Change affectation A is the higher compartment in
  selfInteractionRatio

[HilbertCurve](/packages/HilbertCurve)
------------

                       Changes in version 1.29.1                        

- add `hc_which()`.

[hipathia](/packages/hipathia)
--------

                 Changes in version 2.99.0 (2022-12-01)                 

- Adding parameters uni.terms and GO.terms to hipathia, to compute
  functional activity within this function.

- Adding functions DAcomp, DAtop, DAsummary, DAoverview, define_colors,
  plotVG, DAreport.

- Modifyng structure of objects, by creating object DAdata, which
  includes more information than traditional hipathia results object.
  This includes:
  - Activity values of nodes, paths, and selected functions
  - Extra information about the paths, nodes and functions in rowData
  elements of the SummarizedExperiments

[HPAanalyze](/packages/HPAanalyze)
----------

                        Changes in version 1.17                         

- Changes in version 1.17.0
  + Starting devel for Bioconductor 3.18

- Changes in version 1.17.1
  + Update built-in data to version 22.0.

- Changes in version 1.17.2
  + Update HPA download to correctly download new version 22.0 data.
  + Update documentation.

[hpar](/packages/hpar)
----

                        Changes in version 1.41                         

Changes in version 1.41.1

- Drop getHpa() function and use base R or tidyverse (as illustrated
in the vignette).
- Repalce getHpa(type = "details") by browseHPA().
- Serve data through ExperimentHub.
- Update to HPA release 21.1 <2022-05-31 Tue>
- New datasets added: rnaConsensusTissue, rnaHpaTissue,
rnaGtexTissue,
rnaFantomTissue.
- The dataset rnaGeneTissue becomes rnaGeneTissue21.0 as no longer
available in version 21.1

[IFAA](/packages/IFAA)
----

                       Changes in version 0.99.8                        

- Add partition algorithm in phase 1 and phase 2.

                       Changes in version 0.99.4                        

- Update algorithm in phase 1 to improve speed.

[illuminaio](/packages/illuminaio)
----------

                 Changes in version 0.43.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor release version,
which is now Bioconductor 3.18 for R (>= 4.4.0).

                 Changes in version 0.42.0 (2023-04-25)                 

Notes

- The version number was bumped for the Bioconductor release version,
which is now Bioconductor 3.17 for R (>= 4.3.0).

                 Changes in version 0.41.0 (2022-11-01)                 

- The version number was bumped for the Bioconductor devel version,
which is now Bioconductor 3.17 for R-devel.

[imcRtools](/packages/imcRtools)
---------

                 Changes in version 1.5.5 (2023-04-20)                  

- Bug fix: in corner cases, levels of the permuted factors in the
  testInteractions function do not match the baseline.

- Bug fix: testInteractions now works for cells that are not ordered by
  grouping level

- Added explicitly message for some functions that the output is
  ordered by image

                 Changes in version 1.5.4 (2023-04-04)                  

- Added option to get permutation counts returned from testInteractions

                 Changes in version 1.5.3 (2023-03-23)                  

- Replaced merge function by left_join as safety measure

                 Changes in version 1.5.2 (2022-11-23)                  

- Bug fix minDistToCells: return NA when all cells of an image belong
  to a patch

                 Changes in version 1.5.1 (2022-11-22)                  

- fix axis.ratio to 1 in plotSpatial function and set scales = "fixed"

[immunoClust](/packages/immunoClust)
-----------

                       Changes in version 1.31.12                       

- CHANGES
  * introducing single E,M-steps

                       Changes in version 1.31.8                        

- CHANGES
  * bugfix in plot.immunoClust

                       Changes in version 1.31.5                        

- CHANGES
  * bugfix in meta.export

                       Changes in version 1.31.3                        

- NEW FEATURE
  * introduces meta-clustering with SON/ormalization

                       Changes in version 1.31.2                        

- CHANGES
  * bugfix in meta.export functions

                       Changes in version 1.31.1                        

- CHANGES
  * additional option thres to control meta.SubClustering process
  * code cleaning

[infercnv](/packages/infercnv)
--------

                 Changes in version 1.15.3 (2023-03-29)                 

- Now use the first of the 2 color bars on the left side of the
  observations heatmap to display subclusters when they have been
  calculated and k_obs_groups is not used (>1) when cluster_by_groups=F
  is.

- Add export of infercnv subclusters to Seurat object and features file
  generated by add_to_seurat.

- Add write_phylo option to run() and plot_cnv() to control if a file
  with newick strings for the dendrogram is generated.

- Fully transfer subclustering information when running plot_per_group
  to each annotation's object to take advantage of subcluster being
  displayed on color bars now.

- Fix for "meanvar" sim_method in hspike generation when there are no
  references and one of the observation annotations has only a single
  cell.

- Fix plot_subclusters() when using an object that was processed with
  cluster_by_groups=FALSE.

- Change default k_obs_groups in plot_cnv to 1 instead of 3 to use the
  new subclustering coloring by default.

                 Changes in version 1.15.2 (2023-03-08)                 

- "infercnv_subclusters" plot after subclustering step in run() that
  displays the subclusters is now controlled by its own
  "inspect_subclusters" option.

                 Changes in version 1.15.1 (2023-02-24)                 

- Add helper method plot_subclusters() to plot subclusters as the
  annotations to more easily check if the subclustering settings used
  produce good results or not, so settings can be adjusted.

- Change cluster_by_groups default to True.

                 Changes in version 1.14.2 (2023-03-08)                 

- Add "infercnv_subclusters" plot after subclustering step in run()
  that displays the subclusters. Can be disabled with the exisiting
  no_plot argument.

                 Changes in version 1.14.1 (2023-02-24)                 

- Fix per chromosome subclustering to use all the data and not only the
  data from the last annotation group when there are no outliers
  filtered by z_score (which would always happen when no references are
  defined).

- Apply per subcluster consensus on HMM predictions on the object
  directly when running with per chromosome subclustering enabled,
  rather than only running it before writting txt outputs and figure.
  This makes the infercnv_obj also contain the same values plotted.

- Disable per chromosome subclustering by default. The option remains
  available.

- Change how the subclustering is run on the hspike to avoid issues
  with Leiden settings tuning. Now simply keep the structure of how the
  fake cells are generated.

- Fix to order in which Leiden partition is made into an
  hclust/dendrogram to avoid issues with singleton (as long as there
  are non singletons).

- Fix options stored in backup infercnv objects not being updated
  properly updated when changing the settings on a rerun.

- Fix to properly restart past the Bayesian network step if it has
  already been run and only the BayesMaxPNormal filter has been
  changed.

[INTACT](/packages/INTACT)
------

                       Changes in version 0.99.0                        

NEW FEATURES

- Pilot version

SIGNIFICANT USER-VISIBLE CHANGES

- Pilot version

[InteractiveComplexHeatmap](/packages/InteractiveComplexHeatmap)
-------------------------

                        Changes in version 1.7.1                        

- `row_names_max_width` and `column_names_max_height` are all passed to
  sub-heatmaps.

[IntEREst](/packages/IntEREst)
--------

                       Changes in version 1.24.0                        

NEW FEATURES

- strandSpecific is a new parameter added to interest() and
  interest.sequential(). It indicates that IntEREst is now
  strand-aware. All analysis can be run whilst taking into
  account (or ignoring) the strand specificity of the sequencing
  reads.

- referencePrepare can create references (with collapsed exons)
  that does not ignore the strand information.

SIGNIFICANT USER-VISIBLE CHANGES

- GenomicFiles::reduceByYield function is now used instead of
  BiocParallel::bpiterat for improved efficiency and readability.
  This relates to how analysis related to different pieces of the
  alignment data is distributed over the parallel cores and the
  results are sumarized.

BUG FIXES

- In interest() and and interest.sequential() The oredr of the
  results of the different methods in the result matrix is
  corrected. Correct running of analyses using multiple methods
  with a single command is now possible.

[IRanges](/packages/IRanges)
-------

                       Changes in version 2.34.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- Improve error handling in AtomicList constructors when input is too
  big.

[ISAnalytics](/packages/ISAnalytics)
-----------

                 Changes in version 1.9.3 (2023-04-04)                  

ENHANCES AND REFACTORING

- The package no longer depends on magrittr
- All functionality associated with data.table now it’s completely
optional and will be used internally only if the package is
available
- Several packages were moved from Imports to Suggests - functions
will notify when additional packages are requested for the specific
functionality
- All known deprecated or superseded functions from other packages
have been removed or substituted

FIXES AND GENERAL UPDATES

- Added a new tag “barcode_mux” in available_tags()
- The function HSC_population_size_estimate() now better supports the
computation of estimates from different groups of cell types and
tissues at the same time. The tabular output now contains an
additional column “Timepoints_included” that specifies how many time
points the estimate contains
- Function is_sharing() can now handle better limit cases and has the
option of being parallelised provided appropriate packages are
available (better performance)

DEPRECATIONS & BREAKING CHANGES

- Functions import_parallel_Vispa2Matrices_auto() and
import_parallel_Vispa2Matrices_interactive() are officially defunct
and will not be exported anymore starting from the next release
cycle
- The argument mode of import_parallel_Vispa2Matrices() no longer
accepts INTERACTIVE as a valid option and the interactive mode is
considered now defunct, since the usage is very limiting and limited
- The argument association_file of import_parallel_Vispa2Matrices()
no
longer accepts a string representing a path. Association file import
is delegated solely to its dedicated function from now on.
- The function threshold_filter() is deprecated, since its use is
rather complicated instead of using standard filtering with dplyr or
similar tools

MINOR CHANGES

- default_af_transform() now pads time points based on the maximum
number of characters + 1 in the column

                 Changes in version 1.9.2 (2023-01-26)                  

FIXES AND GENERAL UPDATES

- Fixed an issue with ifelse function in top_abund_tableGrob() - now
the function has a new argument transform_by which is useful for
controlling ordering of columns
- Updated CITATION file
- Package DT has been moved (likely temporarily) in Imports - linked
to issue https://github.com/calabrialab/ISAnalytics/issues/2
- Fixed other typos and minor issues

                 Changes in version 1.9.1 (2022-12-01)                  

- Fixed all tidyselect warnings (internal use of .data$ in selection
context)
- Added bonferroni correction in gene_frequency_fisher

[iSEE](/packages/iSEE)
----

                       Changes in version 2.11.4                        

- Fix bug introduced in 2.11.2.

                       Changes in version 2.11.3                        

- Use standard syntax to include empty icon in relevant places.

                       Changes in version 2.11.2                        

- Make it possible to hide the "Visual Parameters" box.

                       Changes in version 2.11.1                        

- Add tooltip.signif to registerAppOptions() to regulate the number
of
significant digits shown in the tooltip.

[iSEEhub](/packages/iSEEhub)
-------

                        Changes in version 1.1.1                        

MINOR UPDATES

- Update screenshots in vignettes.

[iSEEu](/packages/iSEEu)
-----

                       Changes in version 1.11.2                        

- Adjusted code to correctly parse and rename KEGG pathway
identifiers.

                       Changes in version 1.11.1                        

- Update ggplot2 imports for AggregatedDotPlot.

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

                       Changes in version 1.99.16                       

- Update type: Minor.

- Update for Bioconductor

[kebabs](/packages/kebabs)
------

                       Changes in version 1.33.2                        

- fix in src/Kebabs.h in order to avoid memory problems on 32-bit
  architectures

                       Changes in version 1.33.1                        

- fixes in several help pages to avoid warnings

                       Changes in version 1.33.0                        

- new branch for Bioconductor 3.17 devel

[limma](/packages/limma)
-----

                       Changes in version 3.56.0                        

- 
  Rename readSampleInfoFromGEO() to sampleInfoFromGEO().

- 
  Add new argument `p.value` for topGO() and topKEGG().

- 
  Fix a problem with inconsistent pathway names from KEGG. Due to
  a change on the KEGG website, pathway names from
  getGeneKEGGLinks() have "path:" prefixes but pathway names from
  getKEGGPathwayNames() do not. To make the pathway names
  consistent, getGeneKEGGLinks() now removes the "path:" prefix.

- 
  eBayes() now checks whether `fit` is a list object before
  undertaking other tests.

- 
  Revise voom.Rd to specify that y-values in voom plot are sqrt
  residual standard deviations.

- 
  Update out-of-date URLs in the help pages and the User's Guide.
  Convert reference URLs to DOIs.

- 
  Convert CITATION to use the newer bibentry() instead of
  citEntry().

[magpie](/packages/magpie)
------

                       Changes in version 0.99.10                       

- Added a NEWS.md file to track changes to the package.

[mariner](/packages/mariner)
-------

                       Changes in version 0.99.0                        

Bug fixes and improvements:

- Improve dispatch speed of mergePairs() by removing S4 method
dispatch on all arguments to just x and radius.

- Fix bug in mergePairs() where all pairs are altered during mean of
mode transformation. Now original pairs are preserved when accessed
with getPairClusters().

- Set replace method for counts<- accessor for InteractionMatrix
objects. Helpful for converting DelayedMatrix to matrix.

- Update pixelsToMatrix to preserve metadata columns and include some
additional tests.

- Add plotMatrix() function for plotting matrix data as a heatmap.
Useful for visualizing DelayedMatrices from pullHicMatrices() and
aggHicMatrices(). Compatible with plotgardener package.

- Allow plotMatrix() to accept na.color

- Bug fix in mergePairs() that allows columns named "radius" and/or
"method".

- Swap "binSize" and "files" argument order in pullHicPixels and
pullHicMatrices

- Allow pullHicPixels to overwrite existing HDF5 files.

- Validity checks and functions to access/update the HDF5 paths for
InteractionMatrix objects, even when those paths have been broken.

- Add temporary plotBullseye function.

- Selection functions for selecting indices of a matrix:

- selectCenterPixel
- selectRadius
- selectSubmatrix
- selectCoordinates
- selectBlock
- selectTopLeft
- selectTopRight
- selectBottomRight
- selectBottomLeft
- selectCorners
- selectRows
- selectCols
- selectInner
- selectOuter

- calcLoopEnrichment function for flexibly calculating enrichment of
interactions compared to their local background.

- adjustEnrichment and plotEnrichment for adjusting the loop
enrichment to remove the effect of loop size on enrichment and
visualize this correction across chosen parameters.

                        Changes in version 0.2.0                        

Methods for pulling Hi-C pixels and matrices from .hic files and
storing
them on-disk with HDF5Array and DelayedArray.

Overview of functionality

New or updated functions:

- snapToBins()

- "snaps" ranges in GInteractions objects to their nearest bin
boundary. Allows spanning of multiple bins.

- pullHicPixels() extracts contact frequency from .hic files and
returns an InteractionMatrix object containing a matrix of Hi-C
interactions (rows) and samples (columns).

- Includes counts() accessor for matrix.
- Custom show() method.
- rbind() and cbind() methods.

- pullHicMatrices() extracts submatrices of contact frequency from
.hic files and returns an InteractionArray object containing a
4-dimensional array of Hi-C submatrices, rownames, and colnames.

- Includes counts() accessor for submatrices.
- Custom show() method.
- rbind() and cbind() methods.

- pixelsToMatrices() takes GInteractions containing single pixels
(i.e., each range represents one binSize) and expands ranges such
that there is a buffer of pixels around each range.

- changePixelRes() takes a GInteractions object containing pixels of
interest and is resized to the from resolution/binSize (if its not
already). Then count matrices are extracted for each interaction and
.hic file using the new to resolution. Count matrices are aggregated
by interactions with the supplied aggFUN and a new pixel is selected
with the supplied selectFUN. Allows block processing for large
datasets. The object returned is a GInteractions object with the
updated pixel ranges along with a column containing the aggregated
min/max value for that pixel.

- calcLoopEnrichment() pulls Hi-C pixels and calculates the
enrichment
over background returning a DelayedMatrix of enrichment scores where
rows are interactions and columns are Hi-C files.

- Accessors for GInteractions objects such as seqnames1(), start1(),
end1(), seqnames2(), start2(), end2().

                        Changes in version 0.1.0                        

First pre-release of mariner functionality focused on manipulating,
clustering, and merging paired interactions.

Overview of functionality

- Conversion of paired-range data to GInteractions with
as_ginteractions/makeGInteractionsFromDataFrame

- Functions for manipulating GInteractions and GRanges objects with
binPairs, binRanges, shiftRanges.

- Functions for clustering and merging lists of GInteractions objects
with mergePairs.

- Extensions to GInteractions class with MergedGInteractions, and
DelegatingGInteractions.

- Accessor functions for MergedGInteractions:

- aggPairMcols
- Aggregate metadata columns of clustered interactions.
- getPairClusters
- Return interactions for each cluster of interactions.
- selectionMethod
- Method used to select pairs from each cluster of
interactions.
- sources
- List of names (or indices) used as input for clustering and
merging.
- subsetBySource
- Return interactions unique to each source or combination of
sources.

[MassSpecWavelet](/packages/MassSpecWavelet)
---------------

                 Changes in version 1.65.1 (2023-04-07)                 

- Fix .Call() for R-4.3. Thanks to Steffen Neumann. Closes #5

                 Changes in version 1.64.1 (2023-01-30)                 

- Fix undefined variable in MassSpecWavelet.Rmd

[mastR](/packages/mastR)
-----

                       Changes in version 0.99.9                        

- Updated vignette for BiocStyle packages link functions.

                       Changes in version 0.99.8                        

- Linked external packages in vignette using BiocStyle functions,
fixed the note of using paste in conditions.

                       Changes in version 0.99.7                        

- Improved test coverage depth, modified demo vignette and update
Suggests packages.

                       Changes in version 0.99.6                        

- Removed paste in message().

                       Changes in version 0.99.5                        

- Updated R dependency, delete .Rhistory file, move all generic
functions to a separate file "AllGenerics.R".

                       Changes in version 0.99.4                        

- Deleted *.Rproj file to stop git track.

                       Changes in version 0.99.3                        

- Added *.Rproj into .gitignore.

                       Changes in version 0.99.2                        

- Added non-emtpy return value to man files and update pca plot
function.

                       Changes in version 0.99.1                        

- Fixed notes in R CMD check.

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[MatrixQCvis](/packages/MatrixQCvis)
-----------

                 Changes in version 1.7.7 (2023-04-20)                  

- add package statmod to Suggests

                 Changes in version 1.7.6 (2023-04-19)                  

- add package jpeg to Suggests

                 Changes in version 1.7.5 (2023-04-18)                  

- adjust test ERROR messages

                 Changes in version 1.7.4 (2023-01-26)                  

- add ExperimentHub and GEOquery to NAMESPACE

                 Changes in version 1.7.3 (2023-01-23)                  

- use TCGA RNA-seq and cell line proteomics datasets from ExperimentHub
  in
  vignette to showcase the functionality of the package

                 Changes in version 1.7.2 (2022-11-24)                  

- replace aes_string by aes since aes_string is deprecated in newest
  ggplot2
  version

                 Changes in version 1.7.1 (2022-11-08)                  

- adjust errors in unit tests after updating the packages

[matter](/packages/matter)
------

                        Changes in version 2.1.1                        

BUG FIXES

- Fixed NAMESPACE issue by importing 'Matrix::rowSums()' etc.

                 Changes in version 2.0.1 (2022-11-16)                  

SIGNIFICANT USER-VISIBLE CHANGES

- Updated "User guide" vignette

BUG FIXES

- Fixed combining for matter vectors and lists

- Fixed slow 'matter_mat()' constructor for large # of atoms

- Fixed bug with BPPARAM not passed in 'rowStats()'/'colStats()'

[MBECS](/packages/MBECS)
-----

                        Changes in version 1.3.1                        

- Included PLSDA algorithm.

- WiP requires adjustment of vignette and tests, should work
  though.

[mbQTL](/packages/mbQTL)
-----

                 Changes in version 0.99.6 (2023-03-30)                 

- changes made on vignette

                 Changes in version 0.99.5 (2023-02-25)                 

- changes made based on errors of Bioconductor

                 Changes in version 0.99.4 (2023-02-05)                 

- changes made based on errors of Bioconductor

                 Changes in version 0.99.3 (2023-02-05)                 

- changes made based on Warnings of Bioconductor

                 Changes in version 0.99.2 (2023-02-05)                 

- changes made based on Bioconductor team advisory

                 Changes in version 0.99.1 (2023-01-19)                 

- corrected data loading issues in examples

                 Changes in version 0.99.0 (2023-01-15)                 

- initiated Package and NEWS file

[MEB](/packages/MEB)
---

                       Changes in version 1.13.1                        

- 
  Add new function scMEB().

[metabCombiner](/packages/metabCombiner)
-------------

                        Changes in version 1.9.2                        

- featdata renamed to featData

- associated methods also renamed:
  + adductdata() -> adductData()
  + iddata() -> idData()
  + mzdata() -> mzData()
  + Qdata() -> QData()
  + rtdata() -> rtData()

- new method combineData():
  + merges columns from combinedTable() and featData()

                        Changes in version 1.9.1                        

- Changes to metabData():
  + Major changes to duplicate feature handling
  + bug fix to existing duplicate filtering code
  + duplicate argument reworked to accept a list of parameters
  + option to merge duplicate feature rows added
  + duplicate features are filtered or merged before missingness filter
  + new rowID column added to metabData objects

- new function opts.duplicate():
  + lists default parameters for duplicate feature handling in
  metabData()
  + new duplicate feature merging option (opts.duplicate(resolve =
  "merge"))

- Minor Change to evaluateParams():
  + "score" column changed to "totalScore" (for clarity)

[metabinR](/packages/metabinR)
--------

                 Changes in version 1.2.0 (2023-04-21)                  

- Bump x.y.z version to even y prior to creation of RELEASE_3_17
  branch.

                 Changes in version 1.1.0 (2022-11-01)                  

- Bump x.y.z version to odd y following creation of RELEASE_3_16
  branch.

[MetaboAnnotation](/packages/MetaboAnnotation)
----------------

                         Changes in version 1.3                         

Changes in 1.3.2

- Add mzR as suggested package to ensure package vignettes can be
built.

Changes in 1.3.1

- Small changes in matchSpectra to avoid unnecessary object creation.
- Use backendBpparam to disable parallel processing of matchSpectra
if
the backend does not support it.

[MetCirc](/packages/MetCirc)
-------

                 Changes in version 1.29.2 (2023-04-18)                 

- replace the error message in shinyCircos test

                 Changes in version 1.29.1 (2022-11-23)                 

- replace aes_string by aes since aes_string is deprecated in newest
  ggplot2
  version

[MetNet](/packages/MetNet)
------

                 Changes in version 1.17.1 (2022-11-23)                 

- replace aes_string by aes since aes_string is deprecated in newest
  ggplot2
  version

[mia](/packages/mia)
---

                         Changes in version 1.7                         

- Deprecated assay_name arguments, replaced with assay.type

- Removed abund_values argument

- makePhyloseqFromTreeSE: added option for choosing a tree from
  multiple rowTrees

- mergeSEs: match rows based on all available taxonomy level data on
  rowData

- mergeSEs: fix bug related to equally named variables that are
  different class

- mergeSEs: option for merging multiple assays

- calculateUnifrac: option for specifying the tree from TreeSE

- transformCounts: utilize vegan package

- calculateUnifrac: subset tree based on data

- agglomerateByRank: take into account multiple trees

- loadFromBiom: name columns of rowData based on prefixes

- Deprecate transformSamples, *Features, relabundance, ZTransform,
  relAbundanceCounts

- mergeSEs: faster tree merging

- Faith's index: fix bug that occurred when only one taxon is present

[miaViz](/packages/miaViz)
------

                 Changes in version 1.7.2 (2022-02-17)                  

- assay_name argument changed to assay.type

                         Changes in version 1.7                         

- Fixed plotGraph* that was broked due changes in dependencies

[MicrobiotaProcess](/packages/MicrobiotaProcess)
-----------------

                       Changes in version 1.11.5                        

- update ggdiffclade using the geom layer of ggtree. (2023-03-02)
- update mp_cal_dist to support storing the distance between the
features with action='add'. (2023-03-28)
- fix the special symbol in group name with mp_diff_analysis.
(2023-04-11)

                       Changes in version 1.11.4                        

- update mp_plot_diff_res and mp_plot_diff_boxplot to support
visualize the abundance (not relative abundance). (2022-12-02, Fri)
- fix the tip.label and rownames of assays when tree is provided in
mp_import_dada2 (2022-12-02, Fri)
- add the message information when the differential features was
filtered in the first and second test in mp_diff_analysis.
(2022-12-06, Tue)
- fix the dynamic dots issue of left_join. (2022-12-15, Thu)

                       Changes in version 1.11.3                        

- update mp_plot_abundance to be compatible with the latest ggplot2.
(2022-11-08, Tue)
- add mp_plot_diff_manhattan to visualize the different results with
manhattan plot. (2022-11-21, Mon)

[microSTASIS](/packages/microSTASIS)
-----------

                 Changes in version 0.99.0 (2022-09-25)                 

- Submitted to Bioconductor. Previously on CRAN. Made the following
  changes:

- added support for TreeSumarizedExperiment objects

- changed parallelization from future to BiocParallel

- merged iterative_clustering() and stabilitas() into
  iterativeClustering()

- added a wrapper for automatic multiple paired time points in
  pairedTimes() and iterativeClustering()

- split CV_results() into mSerrorCV() and mSlinesCV()

- added apply loops instead of for and general coding implementations

- removed dependencies on other packages

- removed pre_radarPC() and radarPC()

[miQC](/packages/miQC)
----

                 Changes in version 1.7.1 (2023-01-04)                  

- Added new function, get1DCutoff

[mirTarRnaSeq](/packages/mirTarRnaSeq)
------------

                 Changes in version 1.7.1 (2023-04-12)                  

- fixed ubuntu error

[mistyR](/packages/mistyR)
------

                         Changes in version 1.7                         

                        Changes in version 1.6.1                        

- Discontinue the use of .data due to its deprecation since
tidyselect 1.2.0.
- Update minimum required versions for some tidyverse related
packages
due to function deprecation.

[MoleculeExperiment](/packages/MoleculeExperiment)
------------------

                 Changes in version 0.99.5 (2023-04-17)                 

- Implemented review feedback for Bioconductor submission.

                 Changes in version 0.99.0 (2023-03-30)                 

- Submitted to Bioconductor.

[msa](/packages/msa)
---

                       Changes in version 1.31.7                        

- changes of Version 1.31.6 undone (fix did not work on Mac OS)

- update of ClustalW makefile to avoid problems arising from compiling
  ClustalW with C++ 17: added -std=c++14

                       Changes in version 1.31.6                        

- update of ClustalW makefile to avoid problems arising from compiling
  ClustalW with C++ 17: added -D_HAS_AUTO_PTR_ETC=1

                       Changes in version 1.31.5                        

- update of Muscle source code to avoid problems arising from compiling
  Muscle with C++ 17: renamed type 'byte' to 'MByte'

                       Changes in version 1.31.4                        

- updated src/Muscle/subfams.cpp to avoid conflicting definitions of
  INFINITY on some Mac systems

                       Changes in version 1.31.3                        

- updated config.sub and config.guess in source code of ClustalW and
  ClustalOmega to solve compilation issues on aarch64 (thanks to Yikun
  Jiang
  for contributing this fix!)

                       Changes in version 1.31.2                        

- msa() function changed such that it also works if the package
  is not attached to the workspace

                       Changes in version 1.31.1                        

- update of gc

                       Changes in version 1.31.0                        

- new branch for Bioconductor 3.17 devel

[MSA2dist](/packages/MSA2dist)
--------

                 Changes in version 1.3.1 (2022-11-08)                  

Major changes

Minor improvements and bug fixes

- additional Genetic Codes into base.cpp

[MsBackendMassbank](/packages/MsBackendMassbank)
-----------------

                         Changes in version 1.7                         

Changes in 1.7.4

- Add backendBpparam for MsBackendMassbankSql; parallel processing of
Spectra with MsBackendMassbankSql will silently disable parallel
processing.

Changes in 1.7.3

- Add support for EAD parameter.

Changes in 1.7.2

- Avoid export of retention time data if not available.

Changes in 1.7.1

- Run the full test suite from the Spectra package to validate the
MsBackend implementations.

[MsBackendMgf](/packages/MsBackendMgf)
------------

                         Changes in version 1.7                         

Changes in 1.7.5

- Support import of MS levels from MGF files.

Changes in 1.7.4

- Fix bug in import when m/z values are not sorted.

Changes in 1.7.3

- Enforce m/z values to be increasingly ordered while importing.

Changes in 1.7.2

- Run the full unit test suite from the Spectra package to check
validity of the MsBackendMgf.

Changes in 1.7.1

- Use compress = FALSE for NumericList.

Changes in 1.7.0

- Bioconductor 3.17 developmental branch.

[MsBackendMsp](/packages/MsBackendMsp)
------------

                         Changes in version 1.3                         

Changes in 1.3.1

- Use the full unit test suite from the Spectra package to check
validity of the MsBackendMsp.

[MsBackendSql](/packages/MsBackendSql)
------------

                        Changes in version 0.99                         

Changes in 0.99.7

- Decrease required R version to 4.2.

Changes in 0.99.6

- Add mzR to Suggests to ensure package vignettes can be build
properly.

Changes in 0.99.5

- Add MsBackendOfflineSql backend that re-connects to the database
for
each query.

Changes in 0.99.4

- Add backendBpparam method to ensure parallel processing is disabled
for the MsBackendSql backend.

Changes in 0.99.3

- Add backendMerge method.
- Add parameter data to backendInitialize to allow creating a new
MsBackendSql database and store the values from data in it. This
enables the use of Spectra,setBackend to convert any backend to a
MsBackendSql.
- Implement supportsSetBackend to enable
setBackend,Spectra,MsBackendSql.

Changes in 0.99.2

- Evaluate validity of the MsBackendSql using the full unit test
suite
from the Spectra package.

Changes in 0.99.1

- Address Kayla's package review comments.

Changes in 0.99.0

- Prepare the package for submission to Bioconductor.

                        Changes in version 0.98                         

Changes in 0.98.1

- Rename MsqlBackend to MsBackendSql.

Changes in 0.98.0

- Add vignette.

                         Changes in version 0.1                         

Changes in 0.1.0

- Add parameter blob to allow storing m/z and intensity values as
BLOB
data type in the database. MsBackendSql will use different functions
to retrieve data from a database with this type of storage.

[MsCoreUtils](/packages/MsCoreUtils)
-----------

                        Changes in version 1.11                         

MsCoreUtils 1.11.6

- Fix bug in impute_MinDet(MARGIN = 1) and add unit test.
- Check that package is available and namespace is loaded adding
stopifnot() when calling requireNamespace().

MsCoreUtils 1.11.5

- Add function maxi to determine the maximal intensity value. This
function returns NA_real_ instead of -Inf if all values are missing
or if the length of the input parameter is 0.

MsCoreUtils 1.11.4

- Check if parameter y is increasingly ordered in bin: issue #108.

MsCoreUtils 1.11.3

- Add function sumi to sum intensity values with correct NA handling.

MsCoreUtils 1.11.2

- Reimplement between in C (see issue #105).
- Use symbols to call registered C methods for faster lookup (see PR
#106 and Writing R extensions: Converting a package to use
registration).
- Documentation improvement: explicitly mention impute_mle2() in the
MLE imputation paragraph.

MsCoreUtils 1.11.1

- Add a MARGIN argument to (relevant) imputation functions to support
(and make it explicit) along which dimensions (row or columns)
imputation is performed.
- New impute_mle2() function that uses norm2 (see issue #100).

MsCoreUtils 1.11.0

- New Bioconductor 3.17 (devel) release

[MsDataHub](/packages/MsDataHub)
---------

                        Changes in version 0.99                         

MsDataHub 0.99.3

- Updates based on package review (see
https://github.com/Bioconductor/Contributions/issues/2887#issuecomment-1463655683).

MsDataHub 0.99.2

- Fix title names with make.names() to assert that they will
represent
valid function names after running
ExperimentHub::createHubAccessors().
- Execute chunks in vignette.

MsDataHub 0.99.1

- Run createHubAccessors() when loading the package.

MsDataHub 0.99.0

- Package submission to Biocoductor.

[MsExperiment](/packages/MsExperiment)
------------

                         Changes in version 1.1                         

MsExperiment 1.1.4

- Fix and improve show,MSnExperiment.
- New otherData setter and getter functions.

MsExperiment 1.1.3

- Fix problem in unit test.

MsExperiment 1.1.2

- Add readMsExperiment function (issue #32).

MsExperiment 1.1.1

- Use S4Vectors::findMatches in linkSampleData which improves
performance, especially for larger data sets.

MsExperiment 1.1.0

- Bioconductor release 3.17 (devel).

[MSnbase](/packages/MSnbase)
-------

                        Changes in version 2.25                         

MSnbase 2.25.2

- Fix bug in descendPeaks (see #583)

MSnbase 2.25.1

- Fix remote mztab filename.

MSnbase 2.25.0

- New devel (Bioc 3.17)

[MsQuality](/packages/MsQuality)
---------

                        Changes in version 0.99                         

Changes in version 0.99.9 (2023-04-18)

- enable parallel processing in calculateMetricsFromSpectra.
- rename ticQuantileToQuantileLogRatio to
ticQuartileToQuartileLogRatio
- rename rtOverTicQuantile to rtOverTicQuantiles
- return quartiles instead of quantiles in
precursorIntensityQuartiles
- add mzR to Suggests in DESCRIPTION

Changes in version 0.99.8 (2023-09-02)

- adjust behaviour of metrics function when Spectra object of length
0
is presented, return NA values instead of raising an error

Changes in version 0.99.7 (2023-06-02)

- adjust ticQuantileToQuantileLogRatio that it adheres to Quameter
calculation

Changes in version 0.99.6 (2022-11-29):

- extend section Description in DESCRIPTION
- add sections BugReports and URL in DESCRIPTION
- update dependency to R version 4.2.0
- transfer source code in R/Lee2019-data.R to
inst/sources/Lee2019-data-source.R
- use partial_bundle to reduce the file size of plotly graphics

Changes in version 0.99.5 (2022-10-12):

- add section on alternative software in vignette
- simplify the vignette with regard to dealing with the RPLC and
HILIC
example data set and adjust the Lee2019-data.R accordingly to keep
the RPLC/HILIC information in the dataOrigin slot of the Spectra
object

Changes in version 0.99.4 (2022-10-11):

- add missing comma in DESCRIPTION

Changes in version 0.99.3 (2022-10-11):

- delete Maintainer field in DESCRIPTION
- add instructions to install MsQuality in vignette and README.md via
remotes/BiocManager instead of devtools

Changes in version 0.99.2 (2022-09-20):

- create branch msexperiment to store infrastructure and functions
for
MsExperiment objects
- remove MsExperiment objects since this is still in BioC review and
solely rely on Spectra objects
- adjust documentation for the implemented changes (removal of
MsExperiment)

Changes in version 0.99.1 (2021-11-23)

- simplify calculateMetricsFromSpectra:
- the function does not any longer match the arguments by the
formal arguments of the metric functions
- the function does not any longer combine the parameters
- the additional arguments do not take longer the parameter list
of arguments but comma-separated arguments given to ...
- for all metric functions the ... parameter is added
- adjust the vignette and help pages
- rename functions to camel case

Changes in version 0.99.0 (2021-09-10)

- add quality metrics/functions (metrics based on HUPO-mzQC):
- rtDuration (QC:4000053),
- rtOverTICquantile (QC:4000054),
- rtOverMsQuarters (rtOverMSQuarters),
- ticQuantileToQuantileLogRatio (QC:4000057, QC:4000058),
- numberSpectra (QC:4000059, QC:4000060),
- medianPrecursorMz (QC:4000065),
- rtIQR (QC:4000072),
- rtIQRrate (QC:4000073),
- areaUnderTIC (QC:4000077),
- areaUnderTICRTquantiles (QC:4000078),
- extentIdentifiedPrecursorIntensity (QC:4000125),
- medianTICRTIQR (QC:4000130),
- medianTICofRTRange (QC:4000132),
- mzAcquisitionRange (QC:4000138),
- rtAcquisitionRange (QC:4000139),
- precursorIntensityRange (QC:4000144),
- precursorIntensityQuartiles ((QC:4000167, QC:4000228,
QC:4000233),
- precursorIntensityMean (QC:4000168, QC:4000229, QC:4000234),
- precursorIntensitySD (QC:4000169, QC:4000230, QC:4000235),
- msSignal10XChange (QC:4000172, QC:4000173),
- ratioCharge1over2 (QC:4000174, QC:4000179),
- ratioCharge3over2 (QC:4000175, QC:4000180),
- ratioCharge4over2 (QC:4000176, QC:4000181),
- meanCharge (QC:4000177, QC:4000182),
- medianCharge (QC:4000178, QC:4000183)
- .rt_order_spectra (helper function)
- add the functions calculateMetricsFromSpectra,
calculateMetricsFromMsExperiment to calculate the metrics based on
Spectra and MsExperiment objects
- add functions plotMetric, plotMetric_tibble to visualize the
metrics
- add shiny application shinyMsQuality to interactively visualize the
metrics

[MSstatsConvert](/packages/MSstatsConvert)
--------------

                        Changes in version 1.9.1                        

- Added DIA-NN converter.

[MSstatsTMT](/packages/MSstatsTMT)
----------

                 Changes in version 2.6.1 (2023-02-26)                  

- Minor change: update PhilosophertoMSstatsTMTFormat function
  normalization

[MultiAssayExperiment](/packages/MultiAssayExperiment)
--------------------

                       Changes in version 1.26.0                        

New features

- showReplicated displays the actual colnames of technical replicates
by assay and biological unit.
- The bracket replacement method [<- for MultiAssayExperiment now
also
replaces the names with those from the right-hand side of the
operation, if any (@DarioS, #319)

Bug fixes and minor improvements

- During single assay replacement [[<-, the re-ordering of assays
based on the value input was invalid when empty assays present
(@danielinteractive, #322).
- Permuting assays also updates the order of names in the
MultiAssayExperiment and assays in the sampleMap

[MultimodalExperiment](/packages/MultimodalExperiment)
--------------------

                       Changes in version 0.99.0                        

- Added a NEWS.md file to track changes to the package.

[MungeSumstats](/packages/MungeSumstats)
-------------

                       Changes in version 1.7.18                        

New features

- Check added, ensure BP is between 1 - length of chromosome using
reference chromosome.

                       Changes in version 1.7.17                        

New features

- extra mapping for base-pair position (BP) column added

                       Changes in version 1.7.14                        

Bug fix

- Fix ensembl chain file retrieval so works on all environments

                       Changes in version 1.7.13                        

Bug fix

- write_sumstats:
- Fix indexing issues due to incomplete genome coordinates
sorting:
https://github.com/neurogenomics/MungeSumstats/issues/117
- Add default NULL to ref_genome.
- Check ref_genome (only in conditions where its used).
- sort_coord:
- Renamed .R file from sort_coordinates to match current function
name.
- Add multiple sort_methods, including improved/more robust
data.table-native method.
- Added dedicated unit tests within test-index_tabular.R.
- New helper function: check_numeric:
- Ensures relevant sumstats cols are numeric.
- Added internally to: sort_coord, read_header
- rworkflows.yml:
- Omit Windows runner.
- Turn on run_biocheck
- to_GRanges.R / to_VRanges.R:
- Rename files to match current function names.
- Remove extra extdata files (I think these were created by
accident):
- ALSvcf.vcf.bgz
- ALSvcf.vcf.bgz.bgz
- ALSvcf.vcf.bgz.bgz.tbi
- ALSvcf.vcf.bgz.tbi
- ALSvcf.vcf.gz
- Remove .DS_Store files throughout.
- Don't check for duplicates based on RS ID with Indels, remove these
first.

New features

- Implement rworkflows.
- Removed old Dockerfile (not needed anymore) and workflow yaml.
- Add drop_indels parameter so a user can decide to remove indels
from
sumstats.

                       Changes in version 1.7.12                        

Bug fix

- For downloading files use sed -E rather than sed -r as its
compatible with mac which has issues with sed -r

New features

- For instances where a single column contains CHR, BP, A1 and A2.
The
default order has been updated to CHR:BP:A1:A2 to align with
SPDI format. If your format differs and MSS doesn't pick up on it,
update the column name to the true format e.g. CHR:BP:A2:A1

                       Changes in version 1.7.11                        

New features

- Update to where SNP column is given by the four CHR, BP, A1, A2.
Now, if A1 or A2 is also a separate column, these will be used to
infer the order.

                       Changes in version 1.7.10                        

Bug fix

- further fix for Latex issues when rendering PDF of examples.

                        Changes in version 1.7.9                        

Bug fix

- fix for Latex issues when rendering PDF of examples.

                        Changes in version 1.7.3                        

Bug fix

- fix for offline runs and accessing chain files from 1.7.2.

                        Changes in version 1.7.2                        

New features

- New chain files used for lifting over the genome build from Ensembl
have now been added. These will now be set as the default chain file
instead of UCSC due to licensing issues. The choice to use UCSC
files will still be there but the files will not be stored in the
package themselves, they will instead be downloaded for use on the
fly.

                        Changes in version 1.7.1                        

New features

- The use of the log_folder parameter in format_sumstats() has been
updated. It is still used to point to the directory for the log
files and the log of MungeSumstats messages to be stored. And the
default is still a temporary directory. However, now the name of the
log files (log messages and log outputs) are the same as the name of
the file specified in the save_path parameter with the extension
'_log_msg.txt' and '_log_output.txt' respectively.

[muscat](/packages/muscat)
------

                       Changes in version 1.12.1                        

- fixed various typos in both vignettes

- internal fixes to keep up with 'ggplot2' & 'dplyr' updates

- bug fix in 'simDS' computing means when one group is missing

- bug fix in 'resDS' until testing when 'cpm/frq = TRUE'

[mzR](/packages/mzR)
---

                       Changes in version 2.33.1                        

- fix: update to a new PSI-MS OBO for e.g. ZenoTOF CV term
  Closes #278

[NanoMethViz](/packages/NanoMethViz)
-----------

                        Changes in version 2.6.0                        

- Added preliminary modbam file support.
- Changed rug plot to appear under other geoms. This helps with
visibility of data when methylation values are close to 0.
- Changed heatmap alpha from 0.5 to 1, line width from 1.0 to 1.2 and
line colour from black to darkgrey.
- Changed x-axis limits on plots to be controlled using
coord_cartesian instead of scale_x_continuous. Plots should now
accurately represent data around the boundaries.

[NanoStringNCTools](/packages/NanoStringNCTools)
-----------------

                 Changes in version 1.7.1 (2022-01-25)                  

- Parameters added for update in ggiraph calls

[NanoTube](/packages/NanoTube)
--------

                        Changes in version 1.5.1                        

- We're published in Bioinformatics! CITATION file updated.

[ndexr](/packages/ndexr)
-----

                       Changes in version 1.21.1                        

- UPDATE: Using the RCX package for working with networks.

- all function of the old RCX implementation are removed from this
  package

[netZooR](/packages/netZooR)
-------

                        Changes in version 1.3.2                        

- R implementation of SPIDER
- R implementation of DRAGON
- header argument in pandaPy

[ngsReports](/packages/ngsReports)
----------

                        Changes in version 2.1.5                        

- Added FastpData and FastpDataList classes for working with fastp
  reports

                        Changes in version 2.1.4                        

- Added umi_tools dedup to importNgsLogs

                        Changes in version 2.1.3                        

- Bugfix when importing duplicationMetrics

[nullranges](/packages/nullranges)
----------

                       Changes in version 1.5.19                        

- Remove speedglm dependency as it was removed from CRAN (April
2023).

[oncoscanR](/packages/oncoscanR)
---------

                        Changes in version 1.1.0                        

- To simplify the workflow, the gender of the patient is not taken
into account anymore. That implies that in a male sample, a gain
of 3 extra copies on the X or Y chromosome is considered as a gain
and not an amplification anymore. For female samples, nothing
changes.
- The oncoscan coverage has been corrected to reflect only areas
where
there are groups of probes. Isolated probes where causing issues to
identify arm-level alterations as ChAS segments where never extended
to these probes and the 90% threshold could never be met
(particularly on chromosomal arms 9p and Yq).
- Minor corrections in vignette

[OncoSimulR](/packages/OncoSimulR)
----------

                 Changes in version 4.1.4 (2023-03-08)                  

- Updated exprTk (thanks to pull request from Arash
  Partow).

                 Changes in version 4.1.3 (2023-02-20)                  

- no need to specify CXX_STD = CXX14 (the default since
  R-4.1.0).

                 Changes in version 4.1.2 (2023-02-20)                  

- NAMESPACE: do not import .rbind.data.table (no longer exported
  from data.table 1.14.8, and not needed with R >= 4.0.0).

                  Changes in version 4.1 (2022-11-24)                   

- Faster, cleaner, and better tests of interventions.

[OrganismDbi](/packages/OrganismDbi)
-----------

                       Changes in version 1.41.1                        

MODIFICATIONS

- (v. 1.41.1) Convert OrganismDbi.Rnw to OrganismDbi.Rmd

[orthogene](/packages/orthogene)
---------

                        Changes in version 1.5.2                        

New features

- Remove unnecessary Suggests
- map_species
- get_all_orgs: When species=NULL, now returns an extra columns
called "scientific_name_formatted".
- format_species_name: New args remove_parentheses,
remove_subspecies, remove_subspecies_exceptions
- report_orthologs
- Make much more efficient by only querying ref_genes once.
- Added new internal func report_orthologs_i instead of recursion
to make this easier.
- Ensure that map and report get rbound separately and returned
according to the return_report arg.
- format_species
- Export function that was previously named format_species_name.
- all_species:
- New exported function
- Originally implemented a version of this in EWCE::list_species,
but decided to extend it and export it here.
- get_silhouettes
- Previously was internal func: gather_images.
- Now an exported function.
- plot_orthotree
- Add "Invertebrates" to default clades
- Update README to showcase more functions.

Bug fixes

- drop_non121
- New arg symbol_only to ONLY consider gene symbols (not ensemble
IDs) when identifying non-121 orthologs.
- This make a drastic difference in the number of 1:1 orthologs
that get dropped!
- gather_images
- Update to use newly released rphylopic 1.0.0 which uses the new
phylopic API.
- Add another tryCatch for when the SVG is available but not the
png.
- is_human
- Add "9606" and "homo sapiens sapiens" species ID to list of
options.
- all_genes_babelgene
- Don't filter by support when speices is human, because this
column will always be NA since it's irrelevant for humans.
- Fix unit tests:
- report_orthologs
- dMcast
- Fix stats::pass --> stats::na.pass. Weirdly, only a problem on
Linux. Did base R change a fundamental function name?

                        Changes in version 1.5.1                        

New features

- Bumped to 1.5.1 for Bioc devel 3.17
- Merged upstream devel.
- Now using rworkflows for GHA.
- Removed Dockerfile
- Host orthogene data resources on Zenodo:
- https://doi.org/10.5281/zenodo.7315418
- Upgrade TimeTree phylogeny to v5 (2022):
- 50k+ species --> 137k+ species!
- Replace dplyr::%>% usage with |>
- Add CITATION file

Bug fixes

- prepare_tree:
- Ignore species name case and trim "'" when filtering tree.
- map_species:
- Add trimws step to remove flanking " " or "'".

[OUTRIDER](/packages/OUTRIDER)
--------

                       Changes in version 1.16.1                        

- Add option to restrict FDR correction to sets of genes of interest

- Add opttion to retrieve results based on those FDR values on subsets
  of
  genes

[OutSplice](/packages/OutSplice)
---------

                 Changes in version 0.99.9 (2023-04-17)                 

- Fixed issue when assigning gene ids

                 Changes in version 0.99.8 (2023-04-11)                 

- Fixed error that occured when no significant outliers were
  present

                 Changes in version 0.99.7 (2023-04-07)                 

- Add value section to overview man page

                 Changes in version 0.99.6 (2023-04-06)                 

- Add options to specificy row name columns

- Return function output in addition to saving output files

                 Changes in version 0.99.5 (2023-03-29)                 

- Correct package for build report

                 Changes in version 0.99.4 (2023-03-29)                 

- Correct error with BiocGenerics exported functions

                 Changes in version 0.99.3 (2023-03-28)                 

- Added BiocGenerics to Description and import statements

                 Changes in version 0.99.2 (2023-03-22)                 

- Numerous changes to address comments made during the review
  process

                 Changes in version 0.99.1 (2023-02-27)                 

- Fixed Bioc-devel mailing list registration and added package to
  Watched Tags

                 Changes in version 0.99.0 (2023-01-06)                 

- Submitted to Bioconductor

[pairedGSEA](/packages/pairedGSEA)
----------

                       Changes in version 0.99.6                        

- Improved test coverage and depth
- Coerced paired_ora output to DataFrame from data.table
- Added value field to data man pages

                       Changes in version 0.99.5                        

- Remove non-exported man pages

                       Changes in version 0.99.4                        

- Remove filter_gene_sets option, as it is inherent in fgsea::ora
- Added test depth and coverage for paired_diff and paired_ora
- Fixed wrong storing location for splicing intermediate results

                       Changes in version 0.99.3                        

- Rewrote code base to remove tidyverse dependencies
- Added \code{...} and \link{...} where appropriate in documentation
- Added input parameter checks
- Reduced redundant input parameters from aggregate_pvalue
- Modularized paired_ora and plot_ora
- Added filter_gene_sets parameter to help users reduce gene set bias
- Increased test coverage and depth
- Moved data scripts to inst/script

                       Changes in version 0.99.2                        

- Implement limma as alternative analysis method
- Increased test coverage significantly
- Removed usage of deprecated purrr::when

                       Changes in version 0.99.1                        

- Updated vignette with a brief paragraph on the motivation of
pairedGSEA.
- Updated NAMESPACE to include all imported packages. Suggested
packages are added in notes.

                       Changes in version 0.99.0                        

- Submitted to Bioconductor

[peakPantheR](/packages/peakPantheR)
-----------

                Changes in version 1.13.21 (2023-03-28)                 

- Update for next Bioc release

                 Changes in version 1.13.1 (2022-12-11)                 

- Bugfix: plotHistogram changing warning message

- Bugfix: peakPantheR_plotPeakwidth issue due to ggplot2 behaviour
  changes for
  date axis

- Depreciation warning from ggplot2 `..density..` to `after_stat()`

[pfamAnalyzeR](/packages/pfamAnalyzeR)
------------

                 Changes in version 0.99.7 (2023-04-18)                 

- Added license

                 Changes in version 0.99.6 (2023-04-13)                 

- Updates for BioConductor review

                 Changes in version 0.99.5 (2023-04-13)                 

- Updates for BioConductor review

                 Changes in version 0.99.4 (2023-04-13)                 

- Updates for BioConductor review

- Main user related is a extended vignette

                 Changes in version 0.99.3 (2023-03-15)                 

- Redbuild documentation

                 Changes in version 0.99.2 (2023-03-15)                 

- Update of augment_pfam() to use envelope coordiantes

- Update of augment_pfam() documentation

                 Changes in version 0.99.1 (2023-03-15)                 

- Fixed problematic annoation

                 Changes in version 0.99.0 (2023-02-14)                 

- Version bump for Bioconductor

- Better Description

                 Changes in version 0.1.0 (2022-08-12)                  

- Package introduced

[phantasus](/packages/phantasus)
---------

                       Changes in version 1.19.3                        

- Shiny GAM -> Shiny GATOM

[PharmacoGx](/packages/PharmacoGx)
----------

                        Changes in version 3.3.2                        

- Debugging vignette issues on the Bioconductor build system

                        Changes in version 3.3.1                        

- Added new vignette documenting support for drug combination
modelling new drug combination features added in PharmacoGx >=3.0

[phenomis](/packages/phenomis)
--------

                        Changes in version 1.1.2                        

- minor documentation update

[PhosR](/packages/PhosR)
-----

                        Changes in version 1.9.1                        

- Bug fixes

[PhyloProfile](/packages/PhyloProfile)
------------

                       Changes in version 1.14.0                        

- option for uploading sorted taxon list

                       Changes in version 1.12.6                        

- fixed bug ordering gene IDs when using gene categories

                       Changes in version 1.12.5                        

- fixed bug multiple entries for one (super)taxon ID

                       Changes in version 1.12.4                        

- fixed bug download data

                       Changes in version 1.12.3                        

- fixed bug customized profile not showed

                       Changes in version 1.12.1                        

- added number of co-orthologs and number of taxa in each
  supertaxon

- replaced pfam link by interpro url

[Pigengene](/packages/Pigengene)
---------

                Changes in version 1.25.16 (2023-03-22)                 

Changes in existing functions

- Constant genes are now ignored in the
  compute.pigengene(doWgcna=FALSE,...) function to prevent a run
  time error.

                Changes in version 1.25.12 (2023-03-02)                 

Changes in existing functions

- The default value changed in gene.mapping(leaveNA=FALSE, ...).

                Changes in version 1.25.10 (2023-01-04)                 

Changes in existing functions

- The doWgcna option is added to the compute.pigengene function.

                 Changes in version 1.25.4 (2022-12-01)                 

Changes in existing functions

- The get.enriched.pw function adds gene symbols in the excel
  file.

[planttfhunter](/packages/planttfhunter)
-------------

                       Changes in version 0.99.2                        

CHANGES

- Made small changes suggested by Bioconductor reviewer.

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

[plotgardener](/packages/plotgardener)
------------

                        Changes in version 1.5.3                        

BUG FIXES

- Fixed page viewport parsing bug fixes related to R version 4.3.0
updates.

                        Changes in version 1.5.2                        

BUG FIXES

- annoPixels detects and annotates all pixels for plotHicRectangle
plots.

                        Changes in version 1.5.1                        

BUG FIXES

- yscales for plotHicRectangle and plotHicTriangle reflect distance
in
Hi-C bins.

                        Changes in version 1.5.0                        

Version bump for Bioconductor 3.16 release.

[PoDCall](/packages/PoDCall)
-------

                 Changes in version 1.7.1 (2023-04-04)                  

- Suggest user to change reference well when no threshold for ch
  2

[podkat](/packages/podkat)
------

                       Changes in version 1.31.1                        

- changed arguments in qqplot() method for compatibility with new
  version of
  the standard function in the 'stats' package (added dummy arguments
  to
  avoid errors)

- minor adaptations of help pages

                       Changes in version 1.31.0                        

- new branch for Bioconductor 3.17 devel

[polyester](/packages/polyester)
---------

                       Changes in version 1.99.3                        

- NB function now exported

- note that version 1.99.3 on GitHub was version 1.1.0 on Bioconductor.

                       Changes in version 1.99.2                        

- bug fix in fragment generation (last 2 bases of transcript were never
  sequenced)

[pRoloc](/packages/pRoloc)
------

                        Changes in version 1.39                         

Changes in version 1.39.1

- Update transfer learning vignette to use hpar 1.41.

Changes in version 1.39.0

- New devel version (Bioc 3.17)

[pRolocGUI](/packages/pRolocGUI)
---------

                        Changes in version 2.9.0                        

CHANGES IN VERSION 2.9.0

- New version for Bioc 3.17

CHANGES IN VERSION 2.9.1

- Fix bugs in table selection in explore and compare app
- Fix transparency slider in explore app

[protGear](/packages/protGear)
--------

                 Changes in version 1.3.32 (2022-12-11)                 

- Removed kableExtra due to build error
- added the sampleID files

                 Changes in version 1.3.0 (2022-11-15)                  

- Updated shiny app to load
- Added missing paths for launch_protGear_interactive

[PSMatch](/packages/PSMatch)
-------

                         Changes in version 1.3                         

PSMatch 1.3.3

- New fdr variable (default is always NA_character_ for now) that
defines the spectrum FDR (or any similar/relevant metric that can be
used for filtering - see next item).
- New filterPsmFdr() function that filters based on the fdr variable.

PSMatch 1.3.2

- Specific Matrix::rowSums() to fix error in example.

PSMatch 1.3.1

- Fix type in vignette.

[QFeatures](/packages/QFeatures)
---------

                         Changes in version 1.9                         

QFeatures 1.9.4

- New dropEmptyAssays() function (see issue #184).

QFeatures 1.9.3

- Minor rephrasing in vignette and README.

QFeatures 1.9.2

- feat: filterFeatures() now allows to select assays to filter (i
argument)
- feat: aggregateFeatures() can now take multiple assays
- feat: impute() can now take multiple assays
- feat: processing functions (normalize, scaleTransform,
logTransform,
sweep) can now take multiple assays
- refactor: avoid validObject() when possible
- Use |> rather than %>%.

QFeatures 1.9.1

- fix: solved bug in selectRowData()

QFeatures 1.9.0

- New Bioconductor 3.17 (devel) release

[qsea](/packages/qsea)
----

                       Changes in version 1.25.1                        

- Bugfixes:
  - fixed issue with R4.3 ("cannot xtfrm data frames")
  - fixed documentation
  - updated description file

[qsmooth](/packages/qsmooth)
-------

                 Changes in version 1.15.1 (2022-11-04)                 

- Removed CC BY 4.0 license and replaced with GPL-3

[QuasR](/packages/QuasR)
-----

                       Changes in version 1.40.0                        

USER-VISIBLE CHANGES

- count bam alignments with coordinates but flag 'unmapped' (e.g.
  generated by hisat2) as mapped, so that qCount, qQCReport and
  alignmentStats are consistent

[RaggedExperiment](/packages/RaggedExperiment)
----------------

                       Changes in version 1.24.0                        

Bug fixes and minor improvements

- Use full argument names in unit tests.
- Invoke colSums from MatrixGenerics in unit tests.

[Rarr](/packages/Rarr)
----

                       Changes in version 0.99.9                        

- Response it initial package review (thanks @Kayla-Morrell)
- Provided manual page examples for use_* compression filter
functions.
- Add details of how example data in inst/extdata/zarr_examples was
created.
- General code tidying

                       Changes in version 0.99.8                        

- Patch compression libraries to remove R CMD check warnings about C
functions that might crash R or write to something other than the R
console. Working in Linux only.

                       Changes in version 0.99.7                        

- Allow reading and writing chunks with GZIP compression.
- Add compression level arguments to several compression tools.

                       Changes in version 0.99.6                        

- Allow reading and writing chunks with no compression.
- Enable LZ4 compression for writing.
- Fix bug in blosc compression that could result in larger chunks
than
necessary.
- Improve speed of indexing when combining chunks into the final
output array.

                       Changes in version 0.99.5                        

- Fixed bug when specifying nested chunks, where the chunk couldn't
be
written unless the directory already existed.

                       Changes in version 0.99.4                        

- When writing chunks that overlap the array edge, even the undefined
overhang region should be written to disk.

                       Changes in version 0.99.3                        

- Allow choice between column and row ordering when creating a Zarr
array

                       Changes in version 0.99.2                        

- Catch bug when chunk files contain values outside the array extent.
- Add manual page issues identified by BBS

                       Changes in version 0.99.1                        

- Switch from aws.s3 to paws.storage for S3 data retrieval.

                       Changes in version 0.99.0                        

- Initial Bioconductor submission.

                        Changes in version 0.0.1                        

- Added a NEWS.md file to track changes to the package.

[rawrr](/packages/rawrr)
-----

                 Changes in version 1.7.12 (2023-03-17)                 

- Replace rtinseconds by StartTime &#91;min&#93; as provided by the TFS
  assembly #60.

                 Changes in version 1.7.4 (2023-03-01)                  

- Add rawrr::readTrailer.

[Rbwa](/packages/Rbwa)
----

                        Changes in version 1.3.1                        

- Now supporting Linux ARM64.

[Rcollectl](/packages/Rcollectl)
---------

Initial Submission

- This package helps measure CPU/network/storage consumption of a series of R commands.

[RCy3](/packages/RCy3)
----

                       Changes in version 2.20.0                        

- Add a delay in mergeNetworks method
  - Wait for Cytoscape to finish adding annotations column to Network
  table

- Bug fixes:
  - ellipsis args in createNetworkFromDataFrames, #195
  - setNodeColorBypass: List of hex colors, #188

[ReactomeGSA](/packages/ReactomeGSA)
-----------

                 Changes in version 1.13.1 (2023-04-13)                 

- Fixed documentation mismatch in "plot_heatmap"

[RESOLVE](/packages/RESOLVE)
-------

                        Changes in version 1.2.0                        

- Major code refactoring.

- Package released in Bioconductor 3.17.

[retrofit](/packages/retrofit)
--------

                 Changes in version 0.99.0 (2023-03-05)                 

- Submitted to Bioconductor

[ReUseData](/packages/ReUseData)
---------

                       Changes in version 0.99.23                       

- `ReUseData` helps create data recipes for reproducible data
  processing, where any necessary command-line tools are managed using
  conda and docker.

- `ReUseData` has pre-built data recipes for data downloading/curation
  for common biomedical data resources.

- `ReUseData` supports cloud downloading of pre-generated curated data.

- `ReUseData` has recipe landing pages (https://rcwl.org/dataRecipes/)
  with full annotations and instructions.

[rgoslin](/packages/rgoslin)
-------

                        Changes in version 1.4.0                        

Please note that this Bioconductor version is based on Goslin
version 2.0.0. See the Goslin repository and Goslin C++ repository
for
more details.

BioConductor 3.17 - Changes in 1.4.0

Improvements

- Improved handling of Glycosphingolipids and carbohydrates
- Improved headgroup normalization for Glycosphingolipids.
- Added PMeOH.
- Added TG-EST (estolide) Estolides &#91;GL0305&#93;.
- Added more sterol variants.

Bug Fixes

- Fixed classification of SB1a as Sulfoglycosphingolipids
(sulfatides) &#91;SP0602&#93;.
- Fixed classification of SHex2Cer as Sulfoglycosphingolipids
(sulfatides) &#91;SP0602&#93;.
- Fixed classification of SMGDG as Glycosylalkylacylglycerols
&#91;GL0502&#93;, added synonym seminolipid.
- Fixed classification of SQDG Glycosyldiacylglycerols &#91;GL0501&#93;.
- Fixed classification of sterols.

BioConductor 3.16 - Changes in 1.2.0

- No noteworthy changes.

BioConductor 3.15 - Changes in 1.0.0

Improvements

- Reduced memory consumption.
- Added 'ChE' abbreviation.
- Added FG hydroperoxy to mediator nomenclature, refinement of
mediators.
- Added more sphingosine and sphinganine synonyms.
- Added more ether dialects to LipidMaps grammar.
- Improved handling for SP without explicit OH description.
- Added Sa So support.
- Updated old SP shortcuts.
- Added CholE as abbreviation for cholesterol esters.
- Modifications and improvements for Windows.
- Added column of elements to functional group list and class.
- Added 'ChoE'.
- Added functional group butylperoxy -> BOO.

Bug Fixes

- Fixed handling of LIPID MAPS SP notation.
- Fixed critical bug when parsing LIPID MAPS names.
- Fixed implicit hydroxy count.
- Fixed ACer rule for species level.
- Fixed lcb rule in LipidMaps grammar.
- Fixed S1P and Sa1P handling.
- Fixed gangliosides in Goslin grammar.
- Fixed correct handling of dummy FAs during sorting.
- Fixed segmentation fault in FA parser event handler.

Changes in 0.99.1

- The column names within the data frames returned from the parse*
methods now use column names with dots instead of spaces. This makes
it easier to use the column names unquoted within other R
expressions.
- All parse* methods now return data frames.
- The Messages column has been added to capture parser messages. If
parsing succeeds, this will contain NA and Normalized.Name will
contain the normalized lipid shorthand name.
- Parser implementations have been updated to reflect the latest
lipid
shorthand nomenclature changes. Please see the Goslin repository for
more details.
- Exceptions in the C++ part of the library are captured as warnings
in R. However, if you parse multiple lipid names, exceptions will
not stop the parsing process.

[rGREAT](/packages/rGREAT)
------

                        Changes in version 2.1.9                        

- add `getKEGGGenome()`

                        Changes in version 2.1.8                        

- add `getGenomeDateFromNCBI()`

                        Changes in version 2.1.4                        

- online great: the names of input regions are kept.

- region coordinates (1-based or zero-based) are adjusted in
  `getRegionGeneAssociations()`.

[rhdf5](/packages/rhdf5)
-----

                       Changes in version 2.44.0                        

CHANGES

- h5closeAll() now accepts objects as arguments to allow closing
  a set of HDF5 identifiers.

- Functions H5Teunum_create() and H5Tenum_insert() have been
  included.

- h5set_extent() will now test whether a dataset is chunked and
  inform the user if not. This uses the new function
  H5Dis_chunked().

- The function H5Pset_filter() is now exposed to the user.

BUG FIXES

- Modified how the constant H5S_UNLIMITED was being passed to the
  HDF5 library.  The previous strategy was not working on the
  ARM64 architecture, and leading to failures when trying to
  change the size of a dataset.

- Resolved issue when reporting missing filters where R-to-C
  indexing was being applied twice, resulting in the message:
  "'idx' argument is outside the range of filters set on this
  property list"

[rhdf5filters](/packages/rhdf5filters)
------------

                       Changes in version 1.12.0                        

USER VISIBLE CHANGES

- Compression libraries updated:
- lz4: 1.9.2 🠪 1.9.4
- Added the standalone Zstandard filter. The distributed version of
Zstandard is now 1.5.5.
- Added the VBZ filter.

BUG FIXES

- Ensure ranlib is applied to the lzf library after it is compiled.
This cause linking failures on some systems. (Thanks to Sergey
Fedorov for the report
https://github.com/grimbough/rhdf5filters/pull/18)

[Rhdf5lib](/packages/Rhdf5lib)
--------

                        Changes in version 1.22                         

User visible changes

- Removed pre-compiled library for 32-bit Windows as this platform is
  no
  longer supported by R.

[Rhtslib](/packages/Rhtslib)
-------

                        Changes in version 2.2.0                        

- No significant changes in this version.

[ribosomeProfilingQC](/packages/ribosomeProfilingQC)
-------------------

                       Changes in version 1.11.1                        

- add precedence parameter to readsDistribution

[rifiComparative](/packages/rifiComparative)
---------------

                       Changes in version 0.99.5                        

- Update News file

- Usage of back ticks around variable in vignettes

- deleting details e.g. "Default is" in documentation

                       Changes in version 0.99.4                        

- Fixing error while building the package

                       Changes in version 0.99.3                        

- Improve import file by indicating the functions used of soem packages

- making rifiComparative documenation accessible by using
  "?rifiComparative"

- Improving some naming files as in vignette

                       Changes in version 0.99.2                        

- Supress unecessary comments

- Bugs in vignettes are corrected

- Several typo errors are corrected

                       Changes in version 0.99.1                        

- Significant changes in documentation making it easier to read.

- Bugs in vignettes are corrected

- Suppress warnings in visualization is added

- Several typo errors are corrected

[RnBeads](/packages/RnBeads)
-------

                       Changes in version 2.17.1                        

- Temporarily removed the GLAD dependency and CNV functionality

                       Changes in version 2.17.0                        

- Fixed bug in differential variability analysis

[ROTS](/packages/ROTS)
----

                       Changes in version 1.27.1                        

- Added support for survival analysis

[rpx](/packages/rpx)
---

                         Changes in version 2.7                         

rpx 2.7.4

- Delete entries in default cache before running the unit tests. This
assures that the instances matche any updates in the class
definition.

rpx 2.7.3

- New px_file_types() and file_types() functions to infer mass
spectrometry and proteomics file types based on their extensions.
- New pxSubmissionDate() and pxPublicationDate() accessor functions.

rpx 2.7.2

- Deprecation warning: PXDataset class/methods and PXDataset1()
constructor.

rpx 2.7.1

- Don't rely on a PRIDE project's README.txt files anymore, as it has
been discontinued (see issue #21). The files in a project are now
listed from the remote ftp directory.

[rrvgo](/packages/rrvgo)
-----

                       Changes in version 1.10.1                        

- Add citation to microPublication Biology

[Rsamtools](/packages/Rsamtools)
---------

                        Changes in version 2.16                         

NEW FEATURES

- (v 2.15.1) sortBam() gains support for sorting by tag (byTag) and
  using
  multiple threads (nThreads). (See
  https://github.com/Bioconductor/Rsamtools/issues/46. ; kriemo)

[RTCGAToolbox](/packages/RTCGAToolbox)
------------

                       Changes in version 2.30.0                        

New features

- Include additional options for miRNASeqGeneType and RNAseq2Norm
inputs, see ?getFirehoseData.
- The functions getCNGECorrelation, getDiffExpressedGenes, and
getSurvival have been defunct and removed from the package (see
?'RTCGAToolbox-defunct').
- The RNAseq2Norm argument in getFirehoseData allows additional
options: "RSEM_normalized_log2", "raw_counts", "scaled_estimate"
from the 'preprocessed' tarballs in Firehose. The
"normalized_counts" default remains.

Bug fixes and minor improvements

- getFirehoseData when used with the miRNASeqGene argument was
downloading read counts data rather than RPM. This has been fixed
with the miRNASeqGeneType argument. "read_count" and "cross-mapped"
data are still available but must be entered explicitly in
miRNASeqGeneType.

[S4Arrays](/packages/S4Arrays)
--------

                        Changes in version 1.0.0                        

- First version of the package that is ready for general use.

[S4Vectors](/packages/S4Vectors)
---------

                       Changes in version 0.38.0                        

SIGNIFICANT USER-VISIBLE CHANGES

- The RleTricks vignette was converted from Rnw to Rmd (thanks to
  Beryl Kanali and Jen Wokaty for this conversion).

BUG FIXES

- Improve support of DFrame objects with S3-typed columns (commit
  3b83d5f).

- Fix bug in internal helper bindROWS2() (commit 05cfd3c).

[SAIGEgds](/packages/SAIGEgds)
--------

                        Changes in version 2.0.0                        

NEW FEATURES

- New features to use sparse genetic relationship matrix in generalized
  linear mixed models (GLMMs) according to the SAIGE-GENE paper
  (Zhou et al., 2020)

- a new argument 'grm.mat' in `seqFitNullGLMM_SPA()`: the dense or
  sparse
  genetic relationship matrix (GRM) can be specified via 'grm.mat'

- MAC categories for multiple variance ratios in
  `seqFitNullGLMM_SPA()`,
  `seqAssocGLMM_SPA()`, designed for rare variants

- new function `seqAssocGLMM_SKAT()` for the SKAT aggregate method

- new feature `seqAssocGLMM_ACAT_O()` to include the SKAT aggregate
  method
  for full ACAT-O tests

- faster `seqSAIGE_LoadPval()` when merging multiple input files

                       Changes in version 1.12.4                        

- fix the compiling issue on ARM64
  (see
  https://github.com/AbbVie-ComputationalGenomics/SAIGEgds/issues/8)

                       Changes in version 1.12.1                        

- fix the memory issue because of using deprecated
  tbb::task_scheduler_init
  in RcppParallel

[sangeranalyseR](/packages/sangeranalyseR)
--------------

                        Changes in version 99.1                         

- Base class: SangerReads is designed to store each forward/revers
  reads.

[scAnnotatR](/packages/scAnnotatR)
----------

                 Changes in version 1.5.3 (2023-03-19)                  

- Added paramter to "train_classifier" and "test_classifier" to set the
  characters
  that identify cells as ambiguous.

                 Changes in version 1.5.2 (2023-03-05)                  

- Automatically include parent classifiers if only a subset of
  classifires is selected.

                 Changes in version 1.5.1 (2023-01-17)                  

- Added check to ignore any classifiers / cell types for which genes
  are missing in the dataset.

[SCArray](/packages/SCArray)
-------

                        Changes in version 1.8.0                        

- new function `scNumSplit()`

- update the vignettes

- override S4 functions colsum, rowsum, scale, pmin2, pmax2 for
  SC_GDSMatrix

- new S4 generic functions `row_nnzero()`, `col_nnzero()`,
  `scGetFiles()`,
  `scMemory()`, `scRowMeanVar()`, `scColMeanVar()`

- new functions `scSetMax()`, `scSetMin()`, `scSetBounds()`,
  `scReplaceNA()`

- update S4 functions `[`, `[[`, `names<-`, `dimnames<-`, aperm, Ops,
  Math,
  crossprod, tcrossprod for SC_GDSArray

- override S4 functions rowSums, colSums, rowSums2, colSums2, rowProds,
  colProds, rowMeans, colMeans, rowMeans2, colMeans2, rowVars, colVars,
  rowSds, colSds, rowMins, colMins, rowMaxs, colMaxs, rowRanges,
  colRanges,
  rowAnyNAs, colAnyNAs, rowCollapse, colCollapse for SC_GDSMatrix

[SCArray.sat](/packages/SCArray.sat)
-----------

                       Changes in version 0.99.0                        

- initial version of SCArray.sat

[scater](/packages/scater)
------

                       Changes in version 1.28.0                        

- Change exprs_values (and similar) to assay.type.

- Tweak colouring of violin plots.

- Fix use of block arguments in plotGroupedHeatmap.

- Add scattermore and binning support to various plots (eg
  plotReducedDim).

- swap_rownames works in retrieveCellInfo for altExp now as well
  as in the main assay.

- Add point_shape argument to plotDots and plotPlatePosition.

[scBubbletree](/packages/scBubbletree)
------------

                        Changes in version 1.1.8                        

- Manhattan distance metric implemented

- other linkage functions for clustering implemented

- tracking timestamps for clustering, bubbletree construction, end

[scDblFinder](/packages/scDblFinder)
-----------

                Changes in version 1.13.10 (2023-03-23)                 

- fixed serializing error in multithreading large single samples

- computed thresholds now reported in metadata

                 Changes in version 1.13.7 (2023-01-09)                 

- added possibility to provide the genes/features to use, updated docs

                 Changes in version 1.13.4 (2022-11-21)                 

- fixed bug in samples reporting in split mode (doesn't affect doublets
  scores)

                 Changes in version 1.13.3 (2022-11-20)                 

- updated default parameters according to
  https://arxiv.org/abs/2211.00772

                 Changes in version 1.13.2 (2022-11-11)                 

- added two-pass mode for feature aggregation

[scp](/packages/scp)
---

                         Changes in version 1.9                         

scp 1.9.2

- Updated citation

scp 1.9.1

- Fix minor typo in readSCP() man page

scp 1.9.0

- New Bioconductor 3.17 (devel) release

[scran](/packages/scran)
-----

                        Changes in version 1.28                         

- Added a restricted= option to quickSubCluster() to enable
  subclustering on specific clusters.

[scReClassify](/packages/scReClassify)
------------

                        Changes in version 1.5.1                        

- Bug fixes: handling celltypes with 1 cell in multiAdaSampling
function
- Clarification in multiAdaSampling documentation for parameter
label.

                        Changes in version 1.5.0                        

- Bioconductor 3.16 release

[screenCounter](/packages/screenCounter)
-------------

                        Changes in version 1.2.0                        

- Added support for counting dual barcodes via
  countDualBarcodes().

- Refactored countSingleBarcodes() to support arbitrary numbers
  of substitutions, insertions and deletions.

- Simplified countComboBarcodes() at the expense of dropping
  support for edits inside variable regions.

                        Changes in version 1.0.0                        

- New package screenCounter, containing the barcode counting
  utilities previously in gp.sa.screen.

[scRNAseqApp](/packages/scRNAseqApp)
-----------

                        Changes in version 1.1.1                        

- Fix the bug that `coverage` function is not imported.

                       Changes in version 0.99.27                       

- Fix the subtitles for scATACseq data.

                       Changes in version 0.99.25                       

- Fix a typo when fixing the chromsome name style for scATACseq data.

                       Changes in version 0.99.24                       

- Add tracks for scATACseq data.

                       Changes in version 0.99.23                       

- Add title for ATAC plot buttons.

                       Changes in version 0.99.22                       

- Fix the issue if there is no reference.

                       Changes in version 0.99.21                       

- Fix the typo for `createDataset` function.

                       Changes in version 0.99.20                       

- Add zoom in and out for atac track plots.

                       Changes in version 0.99.19                       

- Update the atac track plots.

                       Changes in version 0.99.18                       

- Simplify the about module.

                       Changes in version 0.99.17                       

- Update documentation for cellgene VIP.

                       Changes in version 0.99.16                       

- fix the bug for cell index for waffleplot.

                       Changes in version 0.99.15                       

- fix the bug if no reference can be retrieved.

                       Changes in version 0.99.14                       

- fix the bug if no markers can be detected.

                       Changes in version 0.99.13                       

- remove aes_string.

                       Changes in version 0.99.12                       

- fix the dropdown names issue.

                       Changes in version 0.99.11                       

- add check/uncheck all to subset cells.

                       Changes in version 0.99.10                       

- fix the multiple loading.

                       Changes in version 0.99.9                        

- fix the long table in vignette.

                       Changes in version 0.99.8                        

- add app_path to app to fix the issue of missing files.

                       Changes in version 0.99.6                        

- add ATAC plots.

                       Changes in version 0.99.5                        

- Update R version dependency from 4.2.0 to 4.3.0

- remove 'CellChat', 'monocle3' and 'SeuratWrappers'.

- re-organize the files in extdata to script and add documentation.

- remove the `eval=FALSE` in vignettes and add `publish_folder` to
  scripts.

- change the parameter of `cat` to `category` for `helper1` function.

- remove `suppressWarnings/Messages`.

- re-format the code style to fit 80 line width and 4 spaces for line
  indents.

- add function to save ATAC peaks.

                       Changes in version 0.99.4                        

- Update documentation.

                       Changes in version 0.99.3                        

- Fix the heatmap donwload button.

                       Changes in version 0.99.2                        

- Fix multiple notes.

                       Changes in version 0.99.1                        

- Prepare for release.

[scuttle](/packages/scuttle)
-------

                        Changes in version 1.10                         

- .subset2index now converts factor inputs to character vectors,
  rather than treating them as integers.

[seq.hotSPOT](/packages/seq.hotSPOT)
-----------

                       Changes in version 0.99.6                        

- Submitted to Bioconductor

[seqArchRplus](/packages/seqArchRplus)
------------

                      Changes in version 0.99.0.10                      

New

- (User-facing) An alternative function to plot motif heatmaps using
the R pkg seqPattern. Function name: plot_motif_heatmaps2()

                      Changes in version 0.99.0.7                       

New

- (User-facing) Now possible to perform per cluster GO term
enrichment
analysis via per_cluster_go_term_enrichments() when a orgDb package
is available for the organism

                      Changes in version 0.99.0.6                       

New

- (User-facing) Examples added majority of functions using
example/dummy data

                      Changes in version 0.99.0.5                       

Bug-fixes

- (User-facing) In per_cluster_annotations():
- Default for clusts set to NULL. This enables cleanly specifying
just the BED file as input to tc_gr
- Help pages now have examples except for the curate_clusters
function. Elaborate documentation for this function is available
as part of the vignette

                      Changes in version 0.99.0.4                       

New features

- (User-facing) In write_seqArchR_cluster_track_bed():
- new argument use_q_bound to choose if you wish to use quantiles
as tag cluster boundaries
- new argument use_as_names to specify any column in info_df to be
used as the name column in the output track BED file of clusters
- dominant_ctss information for each tag cluster is presented as
thickStart thickEnd for ease of visualising. See documentation
for more details
- (User-facing) In per_cluster_annotations():
- tc_gr can now accept a bedfile to read records as a GRanges
object
- Details added in documentation for ways to selectively pick
clusters to annotate

                      Changes in version 0.99.0.3                       

Bux-fixes

- Fixed bugs in curate_clusters() function, and touch up its
documentation

                      Changes in version 0.99.0.2                       

New features

- (User-facing) Parallelization support to speed up annotating
genomic
regions

                       Changes in version 0.99.0                        

New

- Package ready for Bioconductor submission

[SeqArray](/packages/SeqArray)
--------

                       Changes in version 1.40.0                        

- fix the compiler warning: sprintf is deprecated

[SGCP](/packages/SGCP)
----

                 Changes in version 0.99.0 (2022-10-06)                 

- The first version 0.99.0 is submitted to Bioconductor

[ShortRead](/packages/ShortRead)
---------

                        Changes in version 1.58                         

BUG FIXES

- (v 1.57.1, 1.56.1) avoid integer overflow in countFastq()
  https://github.com/Bioconductor/ShortRead/issues/10

[signeR](/packages/signeR)
------

                        Changes in version 2.1.1                        

- Fixed doc issues

- Fixed cpp issues

[signifinder](/packages/signifinder)
-----------

                        Changes in version 1.1.1                        

- Add evaluationSignPlot function to show some technical
  information of the signatures computed.

- Add nametype argument to geneHeatmapSignPlot function to allow
  more gene name ID in data.

- The vignette now contains an example with a single-cell dataset
  and an example with a spatial transcriptomics dataset.

[SimBu](/packages/SimBu)
-----

                 Changes in version 1.1.2 (2023-01-12)                  

- Bugfix in the dataset_seurat function, where custom names for the
  parameters `cell_id_col` and `cell_type_col` did not work (thanks to
  @orange-whale for bringing this up)

- Changed default `remove_bias_in_counts` parameter of function
  `simulate_bulk` to FALSE after discussion with @ZheFrench

                 Changes in version 1.1.1 (2022-11-02)                  

- Bugfix in mirror_db scenario, where cell type fractions were not
  correctly annotated in the cell_fraction output (thanks to @arielah
  for helping out here)

[simplifyEnrichment](/packages/simplifyEnrichment)
------------------

                        Changes in version 1.9.1                        

- `term_similarity_from_gmt()`: previously forgot to split rows.

[singleCellTK](/packages/singleCellTK)
------------

                 Changes in version 2.8.1 (2022-03-10)                  

- Added scanpy wrapper functions for use from console
- Added scanpy UI curated workflow
- Integrated scanpy to a la carte workflow
- Fixed a bug in importing fluidigm dataset
- Updated downloading features in the Shiny app
- Added error checking around Enrichr functions
- Minor tweaks to plot defaults

[SiPSiC](/packages/SiPSiC)
------

                       Changes in version 0.99.2                        

- Removing files LICENCE and tar.gz off the package to meet
  Bioconductor package requirements

[SNPRelate](/packages/SNPRelate)
---------

                       Changes in version 1.34.0                        

- fix the compiler warning: sprintf is deprecated

[SparseArray](/packages/SparseArray)
-----------

                        Changes in version 1.0.0                        

- First version of the package that is ready for general use.

[spaSim](/packages/spaSim)
------

                        Changes in version 1.1.2                        

BUG FIXES

- Fixed documentation of argument jitter in TIS() and
simulate_background_cells().
- Fixed incorrectly plotted image in README.

                        Changes in version 1.1.1                        

SIGIFICANT USER CHANGE

- simulate_background_cells() added an option (method = "Even") to
simulate evenly spaced background images. This change is accompanied
with addition of two parameters, method - to choose the background
cell distribution and jitter

- the parameter to simulate evenly spaced background cells. Tutorials
and function TIS() were modified accordingly.

[SpatialExperiment](/packages/SpatialExperiment)
-----------------

                 Changes in version 1.9.5 (2023-03-02)                  

- bugfix for tissue positions read in incorrect order with
  read10xVisium() in
  datasets with multiple samples (bug introduced in version 1.7.1)

[SpatialFeatureExperiment](/packages/SpatialFeatureExperiment)
------------------------

                        Changes in version 1.1.6                        

- Read images as SpatRaster, in read10xVisiumSFE
- read10xVisiumSFE can also convert full resolution image pixels to
microns based on Visium spot spacing
- read10xVisiumSFE no longer transposes output from read10xVisium so
the spots would match the image by default, and to be consistent
with SpatialExperiment
- Read standard Vizgen MERFISH output with readVizgen
- SpatRasterImage class inheriting from VirtualSpatialImage for
SpatialExperiment compatibility
- Methods of addImg, mirrorImg, and transposeImg for SpatRasterImage
and SFE
- Mirror and transpose SFE objects, operating on both geometries and
images
- Images are cropped when the SFE object is cropped
- Images are also shifted when removeEmptySpace is called

                        Changes in version 1.1.4                        

- Store SFE package version in object and added SFE method of
updateObject to pave way for a potential reimplementation of
spatialGraphs.

                        Changes in version 1.1.3                        

- Use BiocNeighbors for k nearest neighbors and distance based
neighbors, preserving distance info to avoid slow step to refind
distances with sf as done in spdep.
- Added swap_rownames argument in localResult(s) getters so gene
symbols from any rowData column can be used to get local results
stored under Ensembl IDs.

                        Changes in version 1.0.3                        

- Correctly move the geometries when there are multiple samples
- Use translate = FALSE when using localResult setter for geometries
- More helpful error messages when geometries, localResult, or
spatial
graphs are absent

                        Changes in version 1.0.2                        

- Correctly move spatialCoords in removeEmptySpace
- Preserve rownames when setting colGeometry for some of all samples

[spatialHeatmap](/packages/spatialHeatmap)
--------------

                 Changes in version 2.5.4 (2023-04-14)                  

- Developed the SPHM class for storing aSVG, bulk data, single-cell
  data, and matching list.

- Spatial Enrichement: outlier spatial featuers in references are
  supported.

- Shiny App: added K-means clustering, functional enrichment, and
  cluster downloading, redesigned Data Mining and Spatial Enrichment.

- Updated data structures in the background (R functions, Shiny App),
  and data pre-processing steps in SHM.

- Co-visualization: developed a new feature of coloring each cell by
  its own value of a chosen biomolecule; developed a new feature of
  co-visualizing bulk and spatially resolved single-cell data;
  developed a new feature of visualizing deconvolution results.

[SpatialOmicsOverlay](/packages/SpatialOmicsOverlay)
-------------------

                 Changes in version 0.99.0 (2022-04-01)                 

- Submitted to Bioconductor

[speckle](/packages/speckle)
-------

                       Changes in version 0.99.7                        

- Added convertDataToList() to allow propeller to work with any
proportions data
- Added unit tests
- Added minor changes for Bioconductor submission
- Update vignette to include example of analysing proportions data
directly

                       Changes in version 0.99.0                        

- remove classifySex functions in preparation for submission to
Bioconductor
- speckle only contains propeller functions

                        Changes in version 0.0.3                        

- Added functions to classify cells as male or female
- Change propeller transform default to logit

                        Changes in version 0.0.2                        

- Added functions to plot the mean variance relationship of the cell
type counts and proportions
- Added functions to estimate parameters of a Beta distribution
- Added logit transformation option to propeller

                        Changes in version 0.0.1                        

- First version of the speckle package contains propeller functions
to
test for differences in cell type composition between groups of
samples in single cell RNA-Seq data

[Spectra](/packages/Spectra)
-------

                         Changes in version 1.9                         

Changes in 1.9.15

- Fix issue in MsBackendMemory failed to return intensity or m/z
values when peaks data is empty.
- Fix bug in filterPrecursorScan() (see #194 and PR #277).

Changes in 1.9.14

- Fix issue with filterMzValues that would only keep (or remove) the
first matching peak instead of all matching peaks (given ppm and
tolerance). Issue #274.
- Add parameter keep to filterMzRange to support keeping or removing
matching peaks.

Changes in 1.9.13

- Add the backendBpparam method that allows to evaluate whether a
MsBackend supports the provided (or the default) BiocParallel-based
parallel processing setup.
- Minor tweaks in the internal .peaksapply function to avoid
splitting/merging of data if not needed (e.g. if no parallel
processing is performed).
- Minor tweaks in spectra comparison functions to avoid repeated
calling of functions in loops.

Changes in 1.9.12

- Extend the list of available MsBackend backends provided by other
packages (in the README and in the package vignette).

Changes in 1.9.11

- Fix headers in MsBackend vignette.

Changes in 1.9.10

- Add supportsSetBackend method for MsBackend to specify whether a
backend supports setBackend,Spectra.
- setBackend checks using supportsSetBackend whether a backend
supports setBackend.

Changes in 1.9.9

- Refactor setBackend to only split and merge backends if necessary
and to not change dataOrigin of the original backend.
- Support setBackend with MsBackendMemory for an empty Spectra object
(issue #268).
- Disable automatic detection of peak variables for MsBackendMemory
(issue #269).
- Fix issue in Spectra with empty character (issue #267).

Changes in 1.9.8

- Address comments from Michele Stravs regarding the MsBackend
vignette.
- Add additional tests checking for MsBackend compliance.

Changes in 1.9.7

- Add a vignette describing how to build a MsBackend from scratch
(issue #262).
- Extend unit test suite to evaluate validity of MsBackend
implementations.

Changes in 1.9.6

- Replace <= with between calls.

Changes in 1.9.5

- Fix bug in containsMz() when mz isn't ordered (see #258).

Changes in 1.9.4

- Fix error when extracting spectra variables from a MsBackendMzR of
length 0.

Changes in 1.9.3

- Add chunkapply function to split a Spectra into chunks and stepwise
apply a function FUN to each.

Changes in 1.9.2

- combineSpectra on Spectra with read-only backends change backend to
an MsBackendMemory instead of an MsBackendDataFrame.

Changes in 1.9.1

- Expand documentation on compareSpectra for GNPS-like similarity
scoring.

Changes in 1.9.0

- Bioconductor 3.17 developmental version.

[SpectralTAD](/packages/SpectralTAD)
-----------

                 Changes in version 1.15.3 (2023-03-07)                 

- SpectralTAD fixed by Kellen: "There was something in the code
  designed to expand the window if no TADs were detected but it was not
  working properly. No idea why it just started now."

- Change example for SpectralTAD_Par(mat_list, chr= chr, labels =
  labels, cores = 2)

- Update authors, make DESCRIPTION current

[SPIAT](/packages/SPIAT)
-----

                        Changes in version 1.1.6                        

BUG FIXES

- Fixed the mixing score and normalised mixing score calculation.
Each
reference-reference interaction is now counted once (was treated
directional and counted twice) and the fraction of normalised mixing
score is fixed.

                        Changes in version 1.1.5                        

SIGNIFICANT USER-VISIBLE CHANGES

- Return message instead of error when there are no cells of interest
present in the image (identify_neighborhoods()).
- Removed the option of manually defining tumour regions in
identify_bordering_cells(). Removed parameters n_of_polygons and
draw.

NOTES

- Moved the following packages from Imports to Suggestions: graphics,
umap, Rtsne, rlang, ComplexHeatmap and elsa. SpatialExperiment
requires version >= 1.8.0. Removed xROI.

                        Changes in version 1.1.4                        

BUG FIXES

- Fixed error when there are only one cell in the clusters.
(identify_neighborhoods()).

                        Changes in version 1.1.3                        

BUG FIXES

- The calculation of cell types of interest to
All_cells_in_the_structure in
calculate_proportions_of_cells_in_structure() was incorrect. Now
fixed.

                        Changes in version 1.1.2                        

SIGNIFICANT USER-VISIBLE CHANGES

- Re-organised the vignettes.

                        Changes in version 1.1.1                        

BUG FIXES

- Fix bug when Cell.ID column is missing from the spe_object in
identify_neighborhood().

                        Changes in version 1.1.0                        

Development version on Bioconductor 3.17.

[splatter](/packages/splatter)
--------

                 Changes in version 1.24.0 (2022-04-26)                 

- 
  Fixed bugs in splatPopSimulate() where conditional group
  assignments were incorrect when batch effects were applied
  (from Christina Azodi)

- 
  Reduced core dependencies by importing scuttle rather than
  scater (scater is suggested) and making ggplot2 a suggested
  dependency.

[SpliceWiz](/packages/SpliceWiz)
---------

                 Changes in version 1.1.8 (2023-04-17)                  

- Users will no longer need to specify separate folders for processBAM
  and
  collateData output. Instead, using the GUI, processBAM will output to
  the
  `pbOutput` subdirectory inside the specified NxtSE folder

- Buttons in the Experiment creation and loading interfaces have been
  streamlined

- Users can select and de-select events using lasso / box / click
  select tools
  in volcano and scatter plots (previously only select was possible for
  all except
  click)

- A unified event filtering interface has been implemented for all
  visualizations

- A new system for creating coverage plots has been implemented.
  Coverage plots
  are now created in a 3-step process: getCoverageData (to get coverage
  data!),
  getPlotObject (customizes data for ASE event normalization and per
  condition),
  and plotView (which actually generates the plot). This system makes
  it easier
  to refine plots without having to fetch data from the disk everytime.

- All GUI visualization can now be exported directly as pdf files

- Added internal functions to NxtSE object - row_gr() fetches
  EventRegion
  GRanges for each ASE

                 Changes in version 1.1.7 (2023-03-27)                  

- Bugfix: t-test track plots as zero any non-finite p-values (arises
  when all
  normalized coverages are the same value)

- Bugfix: BAM2COV uses correct number of threads now

- Bugfix: fixed duplicate junc_* elements on rbind of NxtSE objects

- Feature: Add abs_deltaPSI as a column in differential analysis result
  (absolute value of delta-PSI)

- Feature: GO interactive plot now displays more information

- Feature: faster retrieval of makeMatrix and makeMeanPSI functions

- Feature: makeSE() now gives more verbose loading information

- Feature: improved performance in getCoverageBins()

- Feature: faster retrieval from Ensembl FTP (using rvest instead of
  XML)

- Feature: slight performance optimization of plotCoverage

                 Changes in version 1.1.6 (2023-02-24)                  

- Gene ontology analysis is available! This is implemented via a
  wrapper to
  fgsea's fora() function (over-representation analysis)

- plotCoverage improved - now exons are plotted at higher resolution,
  and can
  be plotted in isolation (i.e., by removing intronic regions) using
  static
  plots (via as_ggplot_cov()). Further plotCoverage improvements:

- better hover-info for plotly-based events

- unstranded coverage now displays unstranded junction counts

- fixed display of novel transcripts to only display those that are
  supported
  by junction counts in the display data. All annotated transcripts are
  still
  shown

- other miscellaneous bugfixes

- collateData's output is improved. Temporary output files are removed,
  the
  reference is compressed, allowing for lower storage footprint. This
  facilitates file transfer among collaborators. Additionally, COV
  files can be
  copied into the NxtSE folder for file-transfer purposes

- Novel splicing - (optionally) tandem junctions can now be
  extrapolated from
  the data. Given known exons and observed junctions, "putative tandem
  junctions"
  are included among observed tandem junctions, during novel splicing
  event
  generation. This allows for better identification, especially for
  novel
  casette exon skipping.

- Introduced StrictAltSS filter - this removes A5SS/A3SS events for
  which the
  two alternate splice sites are separated by an observed intron.

- Integrated GO analysis into GUI. Heatmaps and event lists in COV can
  now be
  subsetted by top gene ontology categories. GO analysis must first be
  performed
  prior to this option being available.

- Incompatibilities with prior versions:

- buildRef in 1.1.6 now generates gene ontology annotation.

- collateData output is incompatible with prior versions

- processBAM output remains largely unchanged compared with 1.1.5

- NxtSE objects are incompatible with that of prior versions

                 Changes in version 1.1.5 (2022-12-20)                  

- Fix vignettes not building

                 Changes in version 1.1.3 (2022-12-18)                  

- Improved performance of SpliceWiz processBAM() in multi-sample
  setting

- Fixed memory leak in processBAM and BAM2COV functions

- Added edgeR-based differential ASE wrappers, including ability to
  construct
  custom model matrices to model complex experimental designs.

- Overhauled STAR wrappers, added functions to allow STAR genome
  reference to
  be generated (without GTF). A temporary STAR genome can be
  subsequently
  derived by supplying a SpliceWiz reference containing the requisite
  GTF file.

- Included more tandem junctions into novel splicing reference (will
  find more
  novel splicing events compared with versions <=1.1.2)

- collateData's lowMemoryMode will now cap usage to 4 threads (instead
  of 1),
  which is expected to limit RAM usage to ~ 16-20 Gb, depending on
  genome size
  and whether novel splicing mode is on/off. To use even less memory,
  consider

- Slightly improved runtimes of buildRef and collateData functions

- collateData is now single-threaded on Windows (as MulticoreParam is
  not
  available)

                 Changes in version 1.1.2 (2022-11-08)                  

- Implemented time series analysis in limma using splines

- Reduced loading time of makeSE's overlapping intron removal

- Optimised H5 database chunking to speed up data loading times

- Added installation instructions for SpliceWiz using conda environment

- Bugfix: resolved mismatched chromosome issues in collateData

- Bugfix: fixed novel splice counts filtering

- Bugfix: Depth calculation in collateData fixed to properly reflect
  maximum
  splicing across junction

- Bugfix: Fixed plotCoverage / plotGenome by coordinates

[SpotClean](/packages/SpotClean)
---------

                        Changes in version 1.1.2                        

- Fixed a bug in vignettes for column sum on sparse matrix.

                        Changes in version 1.0.1                        

- Fixed a bug due to updated slots in SpatialExperiment objects

[standR](/packages/standR)
------

                        Changes in version 1.3.9                        

- Add function prepareSpatialDecon to help using R package
SpatialDecon after using standR to preprocess GeoMx data.

- New RUV-4: now using the RUV-4 from standR allow you to perform
rank-based analysis such as gene-set scoring with the
RUV-4-normalised count.

- New vignette - A quick start guide to the standR package.

- Many bugs fixed.

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

                       Changes in version 1.29.4                        

- Remove the GUI

- Remove the QCRLSC

[stJoincount](/packages/stJoincount)
-----------

                        Changes in version 1.1.1                        

- Update description and introduction

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

                       Changes in version 1.30.0                        

DEPRECATED AND DEFUNCT

- Finally remove readKallisto() (got deprecated in BioC 3.12 and
  defunct
  in BioC 3.15).

[surfaltr](/packages/surfaltr)
--------

                 Changes in version 1.5.1 (2022-08-12)                  

- Fixed warning

[SVMDO](/packages/SVMDO)
-----

                Changes in version 0.99.24 (2023-04-22)                 

- Fixing problems in NEWS file invisibility

                Changes in version 0.99.23 (2023-04-18)                 

- Acceptance of package for nightly devel build by Bioconductor

                 Changes in version 0.99.0 (2022-12-29)                 

- Sending reopen request for bioconductor submission of updated
package

- Reopening of package submission issue

[synapsis](/packages/synapsis)
--------

               Changes in version 2021-07-02 (2021-07-02)               

Changes:

- Added this news file

[SynExtend](/packages/SynExtend)
---------

                       Changes in version 1.11.8                        

- Fixes various small bugs in MoransI
- Adds some multiprocessing support (more will be added in the
future)
- Slight rework to species trees and their interaction with
ProtWeaver
objects
- ProtWeaver has new attribute speciesTree, can be initialized
with a dendrogram object
- New method SpeciesTree to get species tree from a ProtWeaver
object (or compute one, if it doesn't exist)
- Various internal improvements for Bioc style consistency
- Various documentation updates

                       Changes in version 1.11.7                        

- Adds new optimized dendrapply implementation (overloads
stats::dendrapply)
- AddsHungarianAlgorithm for optimal solving of the linear assignment
problem (O(n^3) complexity)
- Adds new C code for fast computation of Pearson's R and p-value
- Adds new Ancestral.ProtWeaver algorithm for calculating coevolution
from correlated residue changes
- Supporting code and documentation for Ancestral.ProtWeaver
- Other new internal methods
- Various updates and optimizations to internal methods and
documentations
- Updates GRF method to be called CI (for Clustering Information
Distance)
- Method="CI" in PhyloDistance now calculates an approximate p-value
using simulated data from Smith (2020)

                       Changes in version 1.11.6                        

- Adds new Residue method NVDT using gene sequence Natural Vector
with
Dinucleotide and Trinucleotide frequency
- Adds some new C methods to speed up calculation of NVDT
- Fixes .Call() not using PACKAGE="SynExtend"
- Updates to documentation

                       Changes in version 1.11.5                        

- Adds new colocalization algorithm ColocMoran, uses Coloc with
MoransI to correct for phylogenetic signal
- Adds new colocalization algorithm TranscripMI, uses mutual
information of transcriptional direction
- Adds new corrections/checks to allow for transcriptional direction
to be in labels
- Various bugfixes to support new four number labelling scheme
- Various documentation updates
- Adds new function MoransI to calculate Moran's I for a set of
spatially distributed signals

                       Changes in version 1.11.4                        

- Internal code refactor
- ShuffleC now supports reproducibility using R's set.seed
- ShuffleC now support sampling with replacement, performance is
around 2.25x faster than sample

                       Changes in version 1.11.3                        

- Internal bugfixes for JRF Distance--previous commit was incorrectly
calculating values
- Adds new TreeDistance predictor for ProtWeaver, incorporating all
tree distance metrics; these metrics are bundled due to some backend
optimizations that improve performance
- Bugfixes for PhyloDistance
- Adds Random Projection for MirrorTree predictor to solve memory
problems and increase accuracy
- New internal random number generator using xorshift, significantly
faster than sample()
- HammingGL changed to CorrGL, now uses Pearson's R weighted by
p-value
- Refactors internal predictors to reduce size of codebase and remove
redundancies
- Internal ShuffleC function to replicate sample functionality
with 2-6x speedup
- Method GainLoss now uses bootstrapping to estimate a p-value
- Updates to documentation files

                       Changes in version 1.11.2                        

- Adds KF Distance for trees
- Adds Jaccard Robinson Foulds Distance for trees
- Reworks tree distances into PhyloDistance function
- Numerous new documentation pages
- Updates internal functions to use rapply instead of dendrapply to
avoid stack overflow issues due to R recursion

                       Changes in version 1.11.1                        

- Minor bugfix to RF distance
- updates gitignore for workflows

                       Changes in version 1.10.1                        

- Memory leak bugfix

[synlet](/packages/synlet)
------

                       Changes in version 1.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

- Drop doBy, reshape2, dplyr dependency.

- Added data.table support.

[syntenet](/packages/syntenet)
--------

                        Changes in version 1.1.6                        

BUG FIXES

- Added a check for which IQ-TREE version is installed. If IQ-TREE2
is
installed (and not IQ-TREE), arguments and call are modified
accordingly, because IQ-TREE developers changed some parameters
(e.g., -bb is now -B).

NEW FEATURES

- Added a parameter as to parse_collinearity() that allows the
extraction on synteny block information from .collinearity files.
The vignette was updated accordingly.

                        Changes in version 1.1.5                        

BUG FIXES

- Fixed species ID retrieval by adding an exported function named
create_species_id_table() that correctly creates unique species IDs
(3-5 characters), even when the first 5 characters are equal.

- intraspecies_synteny() and interspecies_synteny() (originally
unexported) now take the same output of process_input(), which makes
them consistent with the entire package.

NEW FEATURES

- To make it easier for users who want to run DIAMOND from the
command
line, I added the functions export_sequences() and read_diamond(),
which write processed sequences to FASTA files and read the DIAMOND
output, respectively. An example code on how to run DIAMOND from the
command line has been added to the vignette.

- Included a section in the vignette on how to use syntenet as a
synteny detection program (i.e., to find synteny within a single
genome or between two genomes).

                        Changes in version 1.1.4                        

BUG FIXES

- Replaced sprintf calls with snprintf calls in C++ code to address
warnings in the devel branch of Bioc

UPDATES

- Added CITATION file with reference to published paper

                        Changes in version 1.1.3                        

BUG FIXES

- Tidy evaluation with aes_() was deprecated in ggplot 3.0.0, and
testthat now returns warnings for it. Replaced aes_() with aes() and
.data from the rlang package.

                        Changes in version 1.1.2                        

NEW FEATURES

- Added parameters clust_function and clust_params in
cluster_network() to let users pass any igraph::cluster_* function
to cluster the synteny network.

- Added parameters clust_function and clust_params in plot_profiles()
to let users have more control on the method used to cluster the
distance matrix (columns in phylogenomic profiles).

- Updated vignette to reflect the changes mentioned above and
included
an FAQ item with instructions on how to update the R PATH variable.

                        Changes in version 1.1.1                        

NEW FEATURES

- Ward's clustering of synteny clusters is now performed prior to
plotting in plot_profiles(), not in phylogenomic_profile(). As a
consequence, phylogenomic_profile() now returns only a matrix of
profiles, not a list containing the matrix and an hclust object.

- Added an option to handle names in vector cluster_species as new
names for display in the heatmap. This way, species abbreviations
can be easily replaced with species' full names to make plots look
better.

- Added parameters dist_function and dist_params to allow users to
specify function and parameters to calculate the distance matrix
that will be passed to Ward's clustering.

[systemPipeShiny](/packages/systemPipeShiny)
---------------

                       Changes in version 1.9.04                        

Major Change

- Add video tutorials to main page and all modules.

Minor Change

- Replace the major icons on welcome page, now it has 3 circles
instead of 3, reflecting the 3 major functionalities of SPS,
workflow, visualization, and canvas.
- Changed some namings in workflow creation options
- empty: Start from scratch
- existing: upload custom workflows

Bug Fix

- Fix FontAwesome 6.0 introduced name changes

- One would also need to install develop version of spsComps and
drawer.

- Fix not working tab link on welcome page.

- Fix empty icon problems.

- Fix bug where RNAseq module is always disabled

[TargetSearch](/packages/TargetSearch)
------------

                        Changes in version 2.2.0                        

BUG FIXES

- C code: Refactor the function get_line due to failures in files
  with CR line terminators, commonly found in MacOS systems.

[TBSignatureProfiler](/packages/TBSignatureProfiler)
-------------------

                       Changes in version 1.92.0                        

Bug Fixes

- Fixed the tableAUC bootstrapped confidence interval to be the 2.5
and 97.5 percentiles instead of the 5 and 95 percentiles
- Fixed the upper CI value for the pROC/DeLong AUC CI method in the
bootstrapAUC function

Major Changes

- Changed tableAUC confidence interval default to bootstrapped
instead
of DeLong (pROC argument)

Minor Changes

- Updated the github and website introductions.
- Added Natarajan_7, Kaul_3 signatures
- Added Francisco_OD_2, Kwan_186 signatures
- Shortened some example run times

[TCGAutils](/packages/TCGAutils)
---------

                       Changes in version 1.20.0                        

New features

- makeSummarizedExperimentFromGISTIC and splitAssays are now defunct.

[TDbasedUFE](/packages/TDbasedUFE)
----------

                       Changes in version 0.99.7                        

- Bugfix
- selectFeaureVectorLarge,selectFeaureVectorSmall

                        Changes in version 0.1.0                        

- Replacing menu with shiny and spell check
- Added a NEWS.md file to track changes to the package.

[TDbasedUFEadv](/packages/TDbasedUFEadv)
-------------

                Changes in version 0.99.21 (2023-04-14)                 

- Correction to DESCRIPTION:
- Correct wrong format of DESCRIPTION

                Changes in version 0.99.20 (2023-04-12)                 

- Correction to DESCRIPTION:
- modify DESCRIPTION so as to be disrtnct from the previous
packpage, TDbasedUFE

                Changes in version 0.99.19 (2023-03-28)                 

- Addressing reviews comments:
- Removing some "<<-"s in test.

                Changes in version 0.99.18 (2023-03-27)                 

- Addressing reviews comments:
- Fragmenting testthat toward individual functions

                Changes in version 0.99.17 (2023-03-25)                 

- Addressing reviews comments:
- Replace as.integer with L
- Add link to TDbasedUFE in vignettes
- Correct testthat (with using test_that())

                Changes in version 0.99.16 (2023-03-23)                 

- Addressing reviews comments:
- Renaming class
- Removing one "for" loop
- Adding comment to GUI
- Adding explanation of data in inst/extdata

                Changes in version 0.99.15 (2023-03-21)                 

- Addressing reviews comments:
- Adding argument check
- adding test functions
- Rename and merge vignettes
- and so on

[TENxIO](/packages/TENxIO)
------

                        Changes in version 1.2.0                        

New features

- TENxTSV class has been added to handle compressed and uncompressed
TSV files.

Bug fixes and minor improvements

- The import method for TENxFileList was returning a nested
SummarizedExperiment within the SingleCellExperiment. The counts and
assay(..., i="counts") methods should only return the bare Matrix
rather than the embedded SummarizedExperiment (@dgastn, #2).

[terraTCGAdata](/packages/terraTCGAdata)
-------------

                        Changes in version 1.4.0                        

New features

- findTCGAworkspaces is defunct. Use selectTCGAworkspace instead.

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

                 Changes in version 1.3.2 (2022-11-29)                  

- Updated tomo-seq data from Junker et al., 2014.

[TPP](/packages/TPP)
---

                       Changes in version 3.27.1                        

- Fixed bugs and warnings after dplyr update

[trackViewer](/packages/trackViewer)
-----------

                       Changes in version 1.35.5                        

- use straw for hic import.

                       Changes in version 1.35.4                        

- handel label.parameter.draw for labels.

                       Changes in version 1.35.3                        

- Add documentation for selected labels for mutations.

                       Changes in version 1.35.2                        

- Add support for node label parameters.

                       Changes in version 1.35.1                        

- Add support for single resolution cooler file.

[treeio](/packages/treeio)
------

                       Changes in version 1.23.1                        

- supports converting dendro object (output of
ggdendro::dendro_data()) to a phylo object (2023-03-02, Thu, #95)
- add inner_join() method to allows appending data of a variable
- use nested data structure and tidyr::unnest can extract and
convert the data to a simple tibble data frame (#93)
- update full_join method (#92)
- support standard dplyr UI of by = c( 'columnX' = 'cloumnY')
- drop data from external data.frame that are not presented in the
tree
- use nested column if duplicated rows exist

[tRNA](/packages/tRNA)
----

                 Changes in version 1.17.1 (2022-12-19)                 

- bugfix for gettRNAstructureSeqs

[tximeta](/packages/tximeta)
-------

                       Changes in version 1.17.2                        

- Up to GENCODE 43 (H.s.), M32 (M.m), and Ensembl 109

[tximport](/packages/tximport)
--------

                       Changes in version 1.27.1                        

- Changing license from GPL to LGPL.

[UniProt.ws](/packages/UniProt.ws)
----------

                       Changes in version 2.40.0                        

USER VISIBLE CHANGES

- Convert vignette from RSweave to web-based RMarkdown

BUG FIXES AND MINOR IMPROVEMENTS

- Increase fault tolerance of unit tests
- Add examples to `mapUniProt` documentation

[universalmotif](/packages/universalmotif)
--------------

                       Changes in version 1.18.0                        

NEW FEATURES

- read_meme(readsites.meta.tidy): New option to tidy the output of the
  readsites.meta option into a single data.frame.

MINOR CHANGES

- create_motif(): More robust argument checking.

- The rowMeans, colMeans, rowSums, and colSums generics are now
  imported
  from the MatrixGenerics package instead of BiocGenerics.

BUG FIXES

- Clean up output of argument checks internal to exported functions.

- Delete a reference in IntroductionToSequenceMotifs vignette to a
  non-exported function.

- Delete outdated in MotifManipulation vignette regarding
  convert_motifs
  function.

[updateObject](/packages/updateObject)
------------

                        Changes in version 1.4.0                        

- No changes in this version.

[variancePartition](/packages/variancePartition)
-----------------

                       Changes in version 1.28.9                        

- March 14, 2023
- fix rounding error in makeContrastsDream()
- add Pearson residuals to residuals()

                       Changes in version 1.28.8                        

- March 8, 2023
- add mvTest() with features as list

                       Changes in version 1.28.7                        

- March 7, 2023
- Fix bug in makeContrastsDream() by adding droplevels()

                       Changes in version 1.28.6                        

- March 1, 2023
- diffVar() now fits contrasts estimated in first step

                       Changes in version 1.28.5                        

- Feb 24, 2023
- Fix error in vcov() when samples are dropped due to covariate
having
NA value

                       Changes in version 1.28.4                        

- Feb 2, 2023
- Improve documentation for contrasts

                       Changes in version 1.28.3                        

- Jan 19, 2023
- Improve checking and documentation for contrasts

                       Changes in version 1.28.2                        

- Jan 13, 2023
- canCorPairs() now allows random effects in formula
- but won't change results

                       Changes in version 1.28.1                        

- Jan 04, 2023
- bug fixes

[velociraptor](/packages/velociraptor)
------------

                        Changes in version 1.9.2                        

- Remove column names from reduced dimension matrix in gridVectors().

[Voyager](/packages/Voyager)
-------

                       Changes in version 1.1.12                        

- Plot image behind geometries in all functions that plot geometries
- Added dark theme support for functions that plot geometries

                       Changes in version 1.1.11                        

- Added MULTISPATI PCA
- Added multivariate local Geary's C from Anselin 2019
- Added calculateMultivariate as a unified user interface to
multivariate spatial analyses
- Variogram and variogram map with gstat and related plotting
functions
- Allow non-standard names for local results in plotLocalResult

                       Changes in version 1.1.10                        

- Record parameters used to get spatial results
- Force users to use a new name when running the same method with
different parameters

                        Changes in version 1.1.9                        

- Deprecated show_symbol argument, replacing with swap_rownames to be
consistent with scater

                        Changes in version 1.1.7                        

- Added bbox argument to spatial plotting functions to zoom in with a
bounding box

[weitrix](/packages/weitrix)
-------

                       Changes in version 1.11.1                        

- Don't import colSums, rowSums, which have disappeared from
BiocGenerics.

[xcms](/packages/xcms)
----

                       Changes in version 3.21.5                        

- Fix issue in `chromatogram` after filtering a result object (issue
  #511).

                       Changes in version 3.21.4                        

- Move multtest from Suggests to Imports in dependencies

                       Changes in version 3.21.3                        

- Only fixes in the long running tests

                       Changes in version 3.21.2                        

- Re-write the `reconstructChromPeakSpectra` for DIA data analysis to
  fix an
  issue with chromatographic peaks in overlapping SWATH isolation
  windows and
  generally to improve performance.

                       Changes in version 3.21.1                        

- Fix error with `fillChromPeaks` on sparse data (many empty spectra)
  and peak
  detection performed with `MatchedFilterParam` (issue #653).

- Update to newer function names in the `rgl` package (issue #654).

[zellkonverter](/packages/zellkonverter)
-------------

                       Changes in version 1.10.0                        

Major changes

- 
  Add compatibility with the *anndata* v0.8 H5AD format to the
  the native R writer (By @jackkamm and @mtmorgan)

- 
  Add functions for converting *pandas* arrays used by *anndata*
  when arrays have missing values

Minor changes

- 
  Add Robrecht Cannoodt and Jack Kamm as contributors!

- 
  Minor adjustments to tests to match reader changes

[zenith](/packages/zenith)
------

                        Changes in version 1.1.2                        

- April 20, 2023
- fix warnings for limma

                        Changes in version 1.1.1                        

- fix warnings

                        Changes in version 1.0.7                        

- update docs

                        Changes in version 1.0.6                        

- in zenith() set inter.gene.cor=0.01 to be default to be consistent
with limma::camera

                        Changes in version 1.0.5                        

- bug fix in zenith() when progressbar=FALSE

                        Changes in version 1.0.3                        

- fix issue with corInGeneSet() when some residuals are NA

                        Changes in version 1.0.2                        

- add zenithPR_gsa()

- flag to disable correlation in zenith()

                        Changes in version 1.0.1                        

- fixes to Bioconductor

- improve documentation
- get_GeneOntology() uses getGenesets(...,hierarchical=TRUE)

[ZygosityPredictor](/packages/ZygosityPredictor)
-----------------

                 Changes in version 0.99.0 (2023-02-22)                 

- Submitted to Bioconductor

NEWS from existing Data Experiment Packages
===================================


[CoSIAdata](/packages/CoSIAdata)
---------

New Package Release

- CoSIAdata includes Variance Stabilized Transformation of Read Counts from 
Bgee RNA-Seq Expression Data across six species (Homo sapiens, Mus musculus,
Rattus norvegicius, Danio rerio, Drosophila melanogaster, and 
Caenorhabditis elegans) and more than 132 tissues. Each species has its own 
independent data frame with its unique set of tissue and gene specific 
expression data.
 
- CoSIAdata is meant to be integrated into the CoSIA Package, a visualization 
tool for cross species comparison of expression metrics. However, it can be 
used to conduct independent species, tissue, and gene-specific 
expression analysis.


[curatedTCGAData](/packages/curatedTCGAData)
---------------

                       Changes in version 1.22.0                        

Bug fixes and minor improvements

- When the assays argument was RNASeq2Gene, curatedTCGAData would
  incorrectly include RNASeq2GeneNorm assays. Users who want to return
  both assay types should enter RNASeq2Gene* instead (with an
  asterisk).

[gDNAinRNAseqData](/packages/gDNAinRNAseqData)
----------------

                        Changes in version 1.0.0                        

USER VISIBLE CHANGES

- Submission of the first version to the Bioconductor project on March
  30th, 2023.

[imcdatasets](/packages/imcdatasets)
-----------

                 Changes in version 1.7.3 (2023-02-26)                  

- Added HochSchulz_2022_Melanoma dataset.

                 Changes in version 1.7.2 (2023-02-25)                  

- Added IMMUcan_2022_CancerExample dataset.

- Using the IMMUcan_2022_CancerExample dataset in the vignette.

                 Changes in version 1.7.1 (2023-01-31)                  

- Added full dataset (masks and single cell data) for the
  JacksonFischer
  dataset.

- Added "Zurich" cohort for the JacksonFischer dataset.

- Added full datasets (masks and single cell data) for the Damond...
  dataset.

- Make deprecated functions (`DamondPancreas2019Data`,
  `JacksonFischer2020Data`, and `ZanotelliSpheroids2020Data`) defunct.

[marinerData](/packages/marinerData)
-----------

                       Changes in version 0.99.0                        

NEW FEATURES

- Added a NEWS.md file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- Your main changes to a function foo() or parameter param.

BUG FIXES

- Your bug fixes. See more details at
  http://bioconductor.org/developers/package-guidelines/#news.

[MetaScope](/packages/MetaScope)
---------

                       Changes in version 0.99.0                        

- Pre-Release version of MetaScope
  Bug Fixes

- Fixed check error message about data.table::fread for reading .gz
  files by adding R.utils to imports.
  Major Changes

- Submitted to Bioconductor
  Minor Changes

[microRNAome](/packages/microRNAome)
-----------

                       Changes in version 1.21.1                        

- Updated microRNAome data to the latest published version (v3)
  processed with the latest version of miRge.

[pRolocdata](/packages/pRolocdata)
----------

                       Changes in version 1.37.1                        

- add data from Moloney et al. (2023)

- update syntax in the author field of description

- remove redundant alias for itzhak2016

                       Changes in version 1.37.0                        

- new devel version for Bioc

[scMultiome](/packages/scMultiome)
----------

                       Changes in version 0.99.5                        

NEWS.md setup

- added NEWS.md

[SFEData](/packages/SFEData)
-------

                        Changes in version 1.2.0                        

- Added seqFISH mouse gastrulation data

- Corrected Xenium dataset format

[SingleCellMultiModal](/packages/SingleCellMultiModal)
--------------------

                       Changes in version 1.12.0                        

Bug fixes and minor improvements

- Added Ludwig Geistlinger as author (@lgeistlinger) for contributing
  the GTseq dataset.

[spatialLIBD](/packages/spatialLIBD)
-----------

                       Changes in version 1.11.13                       

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette now has a section describing the data from the
  spatialDLFPC, Visium_SPG_AD, and locus-c projects that were done by
  members of the Keri Martinowich, Kristen R. Maynard, and Leonardo
  Collado-Torres LIBD teams as well as our collaborators.

                       Changes in version 1.11.12                       

SIGNIFICANT USER-VISIBLE CHANGES

- fetch_data("Visium_SPG_AD_Visium_wholegenome_spe""),
  fetch_data("Visium_SPG_AD_Visium_targeted_spe"),
  fetch_data("Visium_SPG_AD_Visium_wholegenome_pseudobulk_spe"), and
  fetch_data("Visium_SPG_AD_Visium_wholegenome_modeling_results") have
  been added. Use this to access data from the
  https://github.com/LieberInstitute/Visium_SPG_AD project.

                       Changes in version 1.11.11                       

SIGNIFICANT USER-VISIBLE CHANGES

- fetch_data("spatialDLPFC_snRNAseq") now works if you want to
  download the snRNA-seq data used in
  http://research.libd.org/spatialDLPFC/.

                       Changes in version 1.11.10                       

BUG FIXES

- read10xVisiumAnalysis() now supports spaceranger version 2023.0208.0
  (internal 10x Genomics version) output files that store analysis
  CSVs under the outs/analysis_csv directory instead of outs/analysis
  and also use the gene_expression_ prefix for each of the analysis
  directories. This was tested with @heenadivecha on files from
  https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/code/02_build_spe/01_build_spe.R.

                       Changes in version 1.11.9                        

SIGNIFICANT USER-VISIBLE CHANGES

- gene_set_enrichment() now internally uses fisher.test(alternative =
  "greater") to test for odds ratios greater than 1. Otherwise odds
  ratios of 0 could be significant.

                       Changes in version 1.11.4                        

SIGNIFICANT USER-VISIBLE CHANGES

- Several changes were made to the default plotting aspect of
  vis_gene(), vis_clus() and related plotting functions. This was done
  with input from @lahuuki and @nick-eagles and is described in more
  detail at
  https://github.com/LieberInstitute/spatialLIBD/commit/8fa8459d8fa881d254824d43e52193bf2c3021c0.
  Most noticeably, the aspect ratio is no longer stretched to fill the
  plotting area, the NA values will be shown with a light grey that
  has alpha blending, and the position of the legends has been made
  consistent between the plots.

                       Changes in version 1.11.3                        

NEW FEATURES

- Added the function frame_limits() and introduced the auto_crop
  argument to vis_clus(), vis_gene() and all related functions. This
  new function enables automatically cropping the image and thus
  adjusting the plotting area which is useful in cases where the image
  is not centered and is not a square. This was based on work by
  @lahuuki at
  https://github.com/LieberInstitute/spatialDLPFC/blob/2dfb58db728c86875a86cc7b4999680ba1f34c38/code/analysis/99_spatial_plotting/01_get_frame_limits.R
  and
  https://github.com/LieberInstitute/spatialDLPFC/blob/ef2952a5a0098a36b09488ebd5e36a902bb11b48/code/analysis/99_spatial_plotting/vis_gene_crop.R.

NEWS from existing Workflows
===================================


[seqpac](/packages/seqpac)
------

- Initial Bioconductor submission


Deprecated and Defunct Packages
===============================

Thirty two software packages were removed from this release (after being deprecated
in Bioc 3.16):
AffyCompatible, BAC, BitSeq, BrainSABER, bridge, cellTree, coexnet, conclus,
ctgGEM, CytoTree, DEComplexDisease, flowCL, flowUtils, gaia, gpart, inveRsion,
IsoGeneGUI, iteremoval, MACPET, PoTRA, rama, Rcade, RNASeqR, scAlign, scMAGeCK,
sojourner, TCGAbiolinksGUI, TDARACNE, TimeSeriesExperiment, TraRe, tspair, XCIR

Please note:  coMET, previously announced as deprecated in 3.16, has been
updated and remain in Bioconductor.

Thirty five software packages are deprecated in this release and will be removed in Bioc 3.18:
alpine, ArrayExpressHTS, ASpediaFI, BiocDockerManager, ChIC, chromswitch,
copynumber, CopywriteR, dasper, epihet, GAPGOM, gcatest, GeneAccord,
genotypeeval, lfa, maanova, metavizr, MethCP, MIGSA, MIMOSA, NanoStringQCPro,
NBSplice, netboxr, NxtIRFcore, ODER, pkgDepTools, PrecisionTrialDrawer,
proBatch, proFIA, pulsedSilac, savR, sigPathway, STAN, TarSeqQC, tscR


Two experimental data packages were removed from this release (after being
deprecated in BioC 3.16):
gatingMLData, RNASeqRData

Two experimental data packages are deprecated in this release and will be
removed in Bioc 3.18:
alpineData, plasFIA

No annotation packages were removed from this release (after being deprecated
in Bioc 3.16).

No annotation packages were deprecated in this release and will be removed in
Bioc 3.18.

No workflow package was removed from this release (after being deprecated in
Bioc 3.16).

No workflow packages were deprecated in this release and will be removed in
3.18.
