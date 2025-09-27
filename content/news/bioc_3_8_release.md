Oct 31, 2018

Bioconductors:

We are pleased to announce Bioconductor 3.8, consisting of 1649
software packages, 360 experiment data packages, 941 annotation
packages, and 23 workflows.

There are 95 new software packages, 21 new data experiment packages, 
2 new workflows, and many updates and improvements
to existing packages; Bioconductor 3.8 is compatible with R 3.5.0,
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

* [Getting Started with Bioconductor 3.8](#getting-started-with-bioconductor-35)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Workflows](#new-Workflows)
* [NEWS from new and existing software packages](#news-from-new-and-existing-software-packages)
* [NEWS from new and existing data experiment packages](#news-from-new-and-existing-data-experiment-packages)
* [NEWS from new and existing workflows](#news-from-new-and-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.8
======================================

To update to or install Bioconductor 3.8:

1. Install R >=3.5.0.  Bioconductor 3.8 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 95 new software packages in this release of Bioconductor.

- [abseqR](https://bioconductor.org/packages/abseqR) AbSeq is a
  comprehensive bioinformatic pipeline for the analysis of sequencing
  datasets generated from antibody libraries and abseqR is one of its
  packages. abseqR empowers the users of abseqPy
  (https://github.com/malhamdoosh/abseqPy) with plotting and
  reporting capabilities and allows them to generate interactive HTML
  reports for the convenience of viewing and sharing with other
  researchers. Additionally, abseqR extends abseqPy to compare
  multiple repertoire analyses and perform further downstream
  analysis on its output.

- [ACE](https://bioconductor.org/packages/ACE) Uses segmented copy
  number data to estimate tumor cell percentage and produce copy
  number plots displaying absolute copy numbers.

- [AffiXcan](https://bioconductor.org/packages/AffiXcan) Impute a
  GReX (Genetically Regulated Expression) for a set of genes in a
  sample of individuals, using a method based on the Total Binding
  Affinity (TBA). Statistical models to impute GReX can be trained
  with a training dataset where the real total expression values are
  known.

- [appreci8R](https://bioconductor.org/packages/appreci8R) The
  appreci8R is an R version of our appreci8-algorithm - A Pipeline
  for PREcise variant Calling Integrating 8 tools. Variant calling
  results of our standard appreci8-tools (GATK, Platypus, VarScan,
  FreeBayes, LoFreq, SNVer, samtools and VarDict), as well as up to 5
  additional tools is combined, evaluated and filtered.

- [artMS](https://bioconductor.org/packages/artMS) artMS provides a
  set of tools for the analysis of proteomics label-free datasets. It
  takes as input the MaxQuant search result output (evidence.txt
  file) and performs quality control, relative quantification using
  MSstats, downstream analysis and integration. artMS also provides a
  set of functions to re-format and make it compatible with other
  analytical tools, including, SAINTq, SAINTexpress, Phosfate, and
  PHOTON.

- [AssessORF](https://bioconductor.org/packages/AssessORF) In order
  to assess the quality of a set of predicted genes for a genome,
  evidence must first be mapped to that genome. Next, each gene must
  be categorized based on how strong the evidence is for or against
  that gene. The AssessORF package provides the functions and class
  structures necessary for accomplishing those tasks, using proteomic
  hits and evolutionarily conserved start codons as the forms of
  evidence.

- [bayNorm](https://bioconductor.org/packages/bayNorm) bayNorm is
  used for normalizing single-cell RNA-seq data.

- [BDMMAcorrect](https://bioconductor.org/packages/BDMMAcorrect)
  Metagenomic sequencing techniques enable quantitative analyses of
  the microbiome. However, combining the microbial data from these
  experiments is challenging due to the variations between
  experiments. The existing methods for correcting batch effects do
  not consider the interactions between variables—microbial taxa in
  microbial studies—and the overdispersion of the microbiome data.
  Therefore, they are not applicable to microbiome data. We develop a
  new method, Bayesian Dirichlet-multinomial regression meta-analysis
  (BDMMA), to simultaneously model the batch effects and detect the
  microbial taxa associated with phenotypes. BDMMA automatically
  models the dependence among microbial taxa and is robust to the
  high dimensionality of the microbiome and their association
  sparsity.

- [BiocNeighbors](https://bioconductor.org/packages/BiocNeighbors)
  Implements exact and approximate methods for nearest neighbor
  detection, in a framework that allows them to be easily switched
  within Bioconductor packages or workflows. The exact algorithm is
  implemented using pre-clustering with the k-means algorithm, as
  described by Wang (2012). This is faster than conventional kd-trees
  for neighbor searching in higher (> 20) dimensional data. The
  approximate method uses the Annoy algorithm. Functions are also
  provided to search for all neighbors within a given distance.
  Parallelization is achieved for all methods using the BiocParallel
  framework.

- [BiocPkgTools](https://bioconductor.org/packages/BiocPkgTools)
  Bioconductor has a rich ecosystem of metadata around packages,
  usage, and build status. This package is a simple collection of
  functions to access that metadata from R. The goal is to expose
  metadata for data mining and value-added functionality such as
  package searching, text mining, and analytics on packages.

- [brainImageR](https://bioconductor.org/packages/brainImageR)
  BrainImageR is a package that provides the user with information of
  where in the human brain their gene set corresponds to. This is
  provided both as a continuous variable and as a
  easily-interpretable image. BrainImageR has additional
  functionality of identifying approximately when in developmental
  time that a gene expression dataset corresponds to. Both the
  spatial gene set enrichment and the developmental time point
  prediction are assessed in comparison to the Allen Brain Atlas
  reference data.

- [breakpointR](https://bioconductor.org/packages/breakpointR) This
  package implements functions for finding breakpoints, plotting and
  export of Strand-seq data.

- [BUScorrect](https://bioconductor.org/packages/BUScorrect)
  High-throughput experimental data are accumulating exponentially in
  public databases. However, mining valid scientific discoveries from
  these abundant resources is hampered by technical artifacts and
  inherent biological heterogeneity. The former are usually termed
  "batch effects," and the latter is often modelled by "subtypes."
  The R package BUScorrect fits a Bayesian hierarchical model, the
  Batch-effects-correction-with-Unknown-Subtypes model (BUS), to
  correct batch effects in the presence of unknown subtypes. BUS is
  capable of (a) correcting batch effects explicitly, (b) grouping
  samples that share similar characteristics into subtypes, (c)
  identifying features that distinguish subtypes, and (d) enjoying a
  linear-order computation complexity.

- [CAMTHC](https://bioconductor.org/packages/CAMTHC) An R package for
  tissue heterogeneity characterization by convex analysis of
  mixtures (CAM). It provides basic functions to perform unsupervised
  deconvolution on mixture expression profiles by CAM and some
  auxiliary functions to help understand the subpopulation-specific
  results. It also implements functions to perform supervised
  deconvolution based on prior knowledge of molecular markers, S
  matrix or A matrix. Combining molecular markers from CAM and from
  prior knowledge can achieve semi-supervised deconvolution of
  mixtures.

- [celaref](https://bioconductor.org/packages/celaref) After the
  clustering step of a single-cell RNAseq experiment, this package
  aims to suggest labels/cell types for the clusters, on the basis of
  similarity to a reference dataset. It requires a table of read
  counts per cell per gene, and a list of the cells belonging to each
  of the clusters, (for both test and reference data).

- [CellTrails](https://bioconductor.org/packages/CellTrails)
  CellTrails is an unsupervised algorithm for the de novo
  chronological ordering, visualization and analysis of single-cell
  expression data. CellTrails makes use of a geometrically motivated
  concept of lower-dimensional manifold learning, which exhibits a
  multitude of virtues that counteract intrinsic noise of single cell
  data caused by drop-outs, technical variance, and redundancy of
  predictive variables. CellTrails enables the reconstruction of
  branching trajectories and provides an intuitive graphical
  representation of expression patterns along all branches
  simultaneously. It allows the user to define and infer the
  expression dynamics of individual and multiple pathways towards
  distinct phenotypes.

- [cicero](https://bioconductor.org/packages/cicero) Cicero computes
  putative cis-regulatory maps from single-cell chromatin
  accessibility data. It also extends monocle 2 for use in chromatin
  accessibility data.

- [COCOA](https://bioconductor.org/packages/COCOA) COCOA is a method
  for understanding variation among samples and can be used with data
  that includes genomic coordinates such as DNA methylation. On a
  high level, COCOA uses a database of "region sets" and principal
  component analysis (PCA) of your data to identify sources of
  variation among samples. A region set is a set of genomic regions
  that share a biological annotation, for instance, transcription
  factor binding regions, histone modification regions, or open
  chromatin regions. COCOA works in both supervised (known groups of
  samples) and unsupervised (no groups) situations and can be used as
  a complement to "differential" methods that find discrete
  differences between groups. COCOA can identify biologically
  meaningful sources of variation between samples and increase
  understanding of variation in your data.

- [compartmap](https://bioconductor.org/packages/compartmap)
  Compartmap performs shrunken A/B compartment inference from
  ATAC-seq and methylation arrays.

- [condcomp](https://bioconductor.org/packages/condcomp) For a given
  clustered data, which can also be split into two conditions, this
  package provides a way to perform a condition comparison on said
  clustered data. The comparison is performed on each cluster.
  Several statistics are used and, when analysed in conjunction, they
  might give some insight regarding the heterogeneity of some of the
  clusters.

- [consensus](https://bioconductor.org/packages/consensus) An
  implementation of the American Society for Testing and Materials
  (ASTM) Standard E691 for interlaboratory testing procedures,
  designed for cross-platform genomic measurements. Given three (3)
  or more genomic platforms or laboratory protocols, this package
  provides interlaboratory testing procedures giving per-locus
  comparisons for sensitivity and precision between platforms.

- [consensusDE](https://bioconductor.org/packages/consensusDE) This
  package allows users to perform DE analysis using multiple
  algorithms. It seeks consensus from multiple methods. Currently it
  supports "Voom", "EdgeR" and "DESeq", but can be easily extended.
  It uses RUV-seq (optional) to remove batch effects.

- [coRdon](https://bioconductor.org/packages/coRdon) Tool for
  analysis of codon usage in various unannotated or KEGG/COG
  annotated DNA sequences. Calculates different measures of CU bias
  and CU-based predictors of gene expressivity, and performs gene set
  enrichment analysis for annotated sequences. Implements several
  methods for visualization of CU and enrichment analysis results.

- [countsimQC](https://bioconductor.org/packages/countsimQC)
  countsimQC provides functionality to create a comprehensive report
  comparing a broad range of characteristics across a collection of
  count matrices. One important use case is the comparison of one or
  more synthetic count matrices to a real count matrix, possibly the
  one underlying the simulations. However, any collection of count
  matrices can be compared.

- [cTRAP](https://bioconductor.org/packages/cTRAP) Compare
  differential gene expression results with those from known cellular
  perturbations (such as gene knock-down, overexpression or small
  molecules) derived from the Connectivity Map. Such analyses allow
  not only to infer the molecular causes of the observed difference
  in gene expression but also to identify small molecules that could
  drive or revert specific transcriptomic alterations.

- [DEqMS](https://bioconductor.org/packages/DEqMS) DEqMS is
  developped on top of Limma. However, Limma assumes same prior
  variance for all genes. In proteomics, the accuracy of protein
  abundance estimates varies by the number of peptides/PSMs
  quantified in both label-free and labelled data. Proteins
  quantification by multiple peptides or PSMs are more accurate.
  DEqMS package is able to estimate different prior variances for
  proteins quantified by different number of PSMs/peptides, therefore
  acchieving better accuracy. The package can be applied to analyze
  both label-free and labelled proteomics data.

- [EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano)
  Volcano plots represent a useful way to visualise the results of
  differential expression analyses. Here, we present a
  highly-configurable function that produces publication-ready
  volcano plots. EnhancedVolcano will attempt to fit as many
  transcript names in the plot window as possible, thus avoiding
  'clogging' up the plot with labels that could not otherwise have
  been read.

- [ERSSA](https://bioconductor.org/packages/ERSSA) The ERSSA package
  takes user supplied RNA-seq differential expression dataset and
  calculates the number of differentially expressed genes at varying
  biological replicate levels. This allows the user to determine,
  without relying on any a priori assumptions, whether sufficient
  differential detection has been acheived with their RNA-seq
  dataset.

- [ExCluster](https://bioconductor.org/packages/ExCluster) ExCluster
  flattens Ensembl and GENCODE GTF files into GFF files, which are
  used to count reads per non-overlapping exon bin from BAM files.
  This read counting is done using the function featureCounts from
  the package Rsubread. Library sizes are normalized across all
  biological replicates, and ExCluster then compares two different
  conditions to detect signifcantly differentially spliced genes.
  This process requires at least two independent biological repliates
  per condition, and ExCluster accepts only exactly two conditions at
  a time. ExCluster ultimately produces false discovery rates (FDRs)
  per gene, which are used to detect significance. Exon log2 fold
  change (log2FC) means and variances may be plotted for each
  significantly differentially spliced gene, which helps scientists
  develop hypothesis and target differential splicing events for
  RT-qPCR validation in the wet lab.

- [FastqCleaner](https://bioconductor.org/packages/FastqCleaner) An
  interactive web application for quality control, filtering and
  trimming of FASTQ files. This user-friendly tool combines a
  pipeline for data processing based on Biostrings and ShortRead
  infrastructure, with a cutting-edge visual environment. Single-Read
  and Paired-End files can be locally processed. Diagnostic
  interactive plots (CG content, per-base sequence quality, etc.) are
  provided for both the input and output files.

- [FCBF](https://bioconductor.org/packages/FCBF) This package
  provides a simple R implementation for the Fast Correlation Based
  Filter described in Yu, L. and Liu, H.; Feature Selection for
  High-Dimensional Data: A Fast Correlation Based Filter
  Solution,Proc. 20th Intl. Conf. Mach. Learn. (ICML-2003),
  Washington DC, 2003 The current package is an intent to make easier
  for bioinformaticians to use FCBF for feature selection, especially
  regarding transcriptomic data.This implies discretizing expression
  (function discretize_exprs) before calculating the features that
  explain the class, but are not predictable by other features. The
  functions are implemented based on the algorithm of Yu and Liu,
  2003 and Rajarshi Guha's implementation from 13/05/2005 available
  (as of 26/08/2018) at http://www.rguha.net/code/R/fcbf.R .

- [FoldGO](https://bioconductor.org/packages/FoldGO) FoldGO is a
  package designed to annotate gene sets derived from expression
  experiments and identify fold-change-specific GO terms.

- [GeneAccord](https://bioconductor.org/packages/GeneAccord) A
  statistical framework to examine the combinations of clones that
  co-exist in tumors. More precisely, the algorithm finds pairs of
  genes that are mutated in the same tumor but in different clones,
  i.e. their subclonal mutation profiles are mutually exclusive. We
  refer to this as clonally exclusive. It means that the mutations
  occurred in different branches of the tumor phylogeny, indicating
  parallel evolution of the clones. Our statistical framework
  assesses whether a pattern of clonal exclusivity occurs more often
  than expected by chance alone across a cohort of patients. The
  required input data are the mutated gene-to-clone assignments from
  a cohort of cancer patients, which were obtained by running
  phylogenetic tree inference methods. Reconstructing the
  evolutionary history of a tumor and detecting the clones is
  challenging. For nondeterministic algorithms, repeated tree
  inference runs may lead to slightly different mutation-to-clone
  assignments. Therefore, our algorithm was designed to allow the
  input of multiple gene-to-clone assignments per patient. They may
  have been generated by repeatedly performing the tree inference, or
  by sampling from the posterior distribution of trees. The tree
  inference methods designate the mutations to individual clones. The
  mutations can then be mapped to genes or pathways. Hence our
  statistical framework can be applied on the gene level, or on the
  pathway level to detect clonally exclusive pairs of pathways. If a
  pair is significantly clonally exclusive, it points towards the
  fact that this specific clone configuration confers a selective
  advantage, possibly through synergies between the clones with these
  mutations.

- [GIGSEA](https://bioconductor.org/packages/GIGSEA) We presented the
  Genotype-imputed Gene Set Enrichment Analysis (GIGSEA), a novel
  method that uses GWAS-and-eQTL-imputed trait-associated
  differential gene expression to interrogate gene set enrichment for
  the trait-associated SNPs. By incorporating eQTL from large gene
  expression studies, e.g. GTEx, GIGSEA appropriately addresses such
  challenges for SNP enrichment as gene size, gene boundary, SNP
  distal regulation, and multiple-marker regulation. The weighted
  linear regression model, taking as weights both imputation accuracy
  and model completeness, was used to perform the enrichment test,
  properly adjusting the bias due to redundancy in different gene
  sets. The permutation test, furthermore, is used to evaluate the
  significance of enrichment, whose efficiency can be largely
  elevated by expressing the computational intensive part in terms of
  large matrix operation. We have shown the appropriate type I error
  rates for GIGSEA (<5%), and the preliminary results also
  demonstrate its good performance to uncover the real signal.

- [glmSparseNet](https://bioconductor.org/packages/glmSparseNet)
  glmSparseNet is an R-package that generalizes sparse regression
  models when the features (e.g. genes) have a graph structure (e.g.
  protein-protein interactions), by including network-based
  regularizers. glmSparseNet uses the glmnet R-package, by including
  centrality measures of the network as penalty weights in the
  regularization. The current version implements regularization based
  on node degree, i.e. the strength and/or number of its associated
  edges, either by promoting hubs in the solution or orphan genes in
  the solution. All the glmnet distribution families are supported,
  namely "gaussian", "poisson", "binomial", "multinomial", "cox", and
  "mgaussian".

- [gpart](https://bioconductor.org/packages/gpart) we provide a new
  SNP sequence partitioning method which partitions the whole SNP
  sequence based on not only LD block structures but also gene
  location information. The LD block construction for GPART is
  performed using Big-LD algorithm, with additional improvement from
  previous version reported in Kim et al.(2017). We also add a
  visualization tool to show the LD heatmap with the information of
  LD block boundaries and gene locations in the package.

- [gwasurvivr](https://bioconductor.org/packages/gwasurvivr)
  gwasurvivr is a package to perform survival analysis using Cox
  proportional hazard models on imputed genetic data.

- [HiCBricks](https://bioconductor.org/packages/HiCBricks) A flexible
  framework for storing and accessing high-resolution Hi-C data
  through HDF files. HiCBricks allows import of Hi-C data through
  various formats such as the 2D matrix format or a generalized
  n-column table formats. In terms of access, HiCBricks offers
  functions to retrieve values from genomic loci separated by a
  certain distance, or the ability to fetch matrix subsets using word
  alike terms. HiCBricks will at a later point offer the ability to
  fetch multiple matrix subsets using fewer calls. It offers the
  capacity to store GenomicRanges that may be associated to a
  particular Hi-C experiment, to do basic ranges overlap (any,
  within) with the Hi-C experiment associated Ranges object and also
  to store any metadata that users may think to be relevant for their
  Hi-C experiment. Finally, you can do TAD calls with LSD and create
  pretty heatmaps.

- [hierinf](https://bioconductor.org/packages/hierinf) Tools to
  perform hierarchical inference for one or multiple studies / data
  sets based on high-dimensional multivariate (generalised) linear
  models. A possible application is to perform hierarchical inference
  for GWA studies to find significant groups or single SNPs (if the
  signal is strong) in a data-driven and automated procedure. The
  method is based on an efficient hierarchical multiple testing
  correction and controls the FWER. The functions can easily be run
  in parallel.

- [HIREewas](https://bioconductor.org/packages/HIREewas) In
  epigenome-wide association studies, the measured signals for each
  sample are a mixture of methylation profiles from different cell
  types. The current approaches to the association detection only
  claim whether a cytosine-phosphate-guanine (CpG) site is associated
  with the phenotype or not, but they cannot determine the cell type
  in which the risk-CpG site is affected by the phenotype. We propose
  a solid statistical method, HIgh REsolution (HIRE), which not only
  substantially improves the power of association detection at the
  aggregated level as compared to the existing methods but also
  enables the detection of risk-CpG sites for individual cell types.
  The "HIREewas" R package is to implement HIRE model in R.

- [HPAanalyze](https://bioconductor.org/packages/HPAanalyze) Provide
  functions for retrieving, exploratory analyzing and visualizing the
  Human Protein Atlas data.

- [iasva](https://bioconductor.org/packages/iasva) Iteratively
  Adjusted Surrogate Variable Analysis (IA-SVA) is a statistical
  framework to uncover hidden sources of variation even when these
  sources are correlated. IA-SVA provides a flexible methodology to
  i) identify a hidden factor for unwanted heterogeneity while
  adjusting for all known factors; ii) test the significance of the
  putative hidden factor for explaining the unmodeled variation in
  the data; and iii), if significant, use the estimated factor as an
  additional known factor in the next iteration to uncover further
  hidden factors.

- [icetea](https://bioconductor.org/packages/icetea) icetea
  (Integrating Cap Enrichment with Transcript Expression Analysis)
  provides functions for end-to-end analysis of multiple 5'-profiling
  methods such as CAGE, RAMPAGE and MAPCap, beginning from raw reads
  to detection of transcription start sites using replicates. It also
  allows performing differential TSS detection between group of
  samples, therefore, integrating the mRNA cap enrichment information
  with transcript expression analysis.

- [INDEED](https://bioconductor.org/packages/INDEED) An
  Implementation of Integrated Differential Expression and
  Differential Network Analysis of Omic Data. The differential
  network is obtained based on partial correlation or correlation.

- [ipdDb](https://bioconductor.org/packages/ipdDb) All alleles from
  the IPD IMGT/HLA <https://www.ebi.ac.uk/ipd/imgt/hla/> and IPD KIR
  <https://www.ebi.ac.uk/ipd/kir/> database for Homo sapiens.
  Reference: Robinson J, Maccari G, Marsh SGE, Walter L, Blokhuis J,
  Bimber B, Parham P, De Groot NG, Bontrop RE, Guethlein LA, and
  Hammond JA KIR Nomenclature in non-human species Immunogenetics
  (2018), in preparation.

- [IsoCorrectoR](https://bioconductor.org/packages/IsoCorrectoR)
  IsoCorrectoR is a tool for correcting natural isotope abundance
  contributions in tracing experiments.

- [KinSwingR](https://bioconductor.org/packages/KinSwingR) KinSwingR
  integrates phosphosite data derived from mass-spectrometry data and
  kinase-substrate predictions to predict kinase activity. Several
  functions allow the user to build PWM models of kinase-subtrates,
  statistically infer PWM:substrate matches, and integrate these data
  to infer kinase activity.

- [levi](https://bioconductor.org/packages/levi) The tool integrates
  data from biological networks with transcriptomes, displaying a
  heatmap with surface curves to evidence the altered regions.

- [LoomExperiment](https://bioconductor.org/packages/LoomExperiment)
  The LoomExperiment class provide a means to easily convert
  Bioconductor's "Experiment" classes to loom files and vice versa.

- [LRBaseDbi](https://bioconductor.org/packages/LRBaseDbi) Interface
  to construct LRBase package (LRBase.XXX.eg.db).

- [maser](https://bioconductor.org/packages/maser) This package
  provides functionalities for analysis, annotation and visualizaton
  of alternative splicing events.

- [methylGSA](https://bioconductor.org/packages/methylGSA) The main
  functions for methylGSA are methylglm and methylRRA. methylGSA
  implements logistic regression adjusting number of probes as a
  covariate. methylRRA adjusts multiple p-values of each gene by
  Robust Rank Aggregation. For more detailed help information, please
  see the vignette.

- [MetID](https://bioconductor.org/packages/MetID) This package uses
  an innovative network-based approach that will enhance our ability
  to determine the identities of significant ions detected by LC-MS.

- [MetNet](https://bioconductor.org/packages/MetNet) MetNet contains
  functionality to infer metabolic network topologies from
  quantitative data and high-resolution mass/charge information.
  Using statistical models (including correlation, mutual
  information, regression and Bayes statistics) and quantitative data
  (intensity values of features) adjacency matrices are inferred that
  can be combined to a consensus matrix. Mass differences calculated
  between mass/charge values of features will be matched against a
  data frame of supplied mass/charge differences referring to
  transformations of enzymatic activities. In a third step, the two
  matrices are combined to form a adjacency matrix inferred from both
  quantitative and structure information.

- [miRSM](https://bioconductor.org/packages/miRSM) The package aims
  to identify miRNA sponge modules by integrating expression data and
  miRNA-target binding information. It provides several functions to
  study miRNA sponge modules, including popular methods for inferring
  gene modules (candidate miRNA sponge modules), and a function to
  identify miRNA sponge modules, as well as a function to conduct
  functional analysis of miRNA sponge modules.

- [mixOmics](https://bioconductor.org/packages/mixOmics) Multivariate
  methods are well suited to large omics data sets where the number
  of variables (e.g. genes, proteins, metabolites) is much larger
  than the number of samples (patients, cells, mice). They have the
  appealing properties of reducing the dimension of the data by using
  instrumental variables (components), which are defined as
  combinations of all variables. Those components are then used to
  produce useful graphical outputs that enable better understanding
  of the relationships and correlation structures between the
  different data sets that are integrated. mixOmics offers a wide
  range of multivariate methods for the exploration and integration
  of biological datasets with a particular focus on variable
  selection. The package proposes several sparse multivariate models
  we have developed to identify the key variables that are highly
  correlated, and/or explain the biological outcome of interest. The
  data that can be analysed with mixOmics may come from high
  throughput sequencing technologies, such as omics data
  (transcriptomics, metabolomics, proteomics, metagenomics etc) but
  also beyond the realm of omics (e.g. spectral imaging). The methods
  implemented in mixOmics can also handle missing values without
  having to delete entire rows with missing data. A non exhaustive
  list of methods include variants of generalised Canonical
  Correlation Analysis, sparse Partial Least Squares and sparse
  Discriminant Analysis. Recently we implemented integrative methods
  to combine multiple data sets: N-integration with variants of
  Generalised Canonical Correlation Analysis and P-integration with
  variants of multi-group Partial Least Squares.

- [mlm4omics](https://bioconductor.org/packages/mlm4omics) To conduct
  Bayesian inference regression for responses with multilevel
  explanatory variables and missing values; It uses function from
  'Stan', a software to implement posterior sampling using
  Hamiltonian MC and its variation Non-U-Turn algorithms. It
  implements the posterior sampling of regression coefficients from
  the multilevel regression models. The package has two main
  functions to handle not-missing-at-random missing responses and
  left-censored with not-missing-at random responses. The purpose is
  to provide a similar format as the other R regression functions but
  using 'Stan' models.

- [MPRAnalyze](https://bioconductor.org/packages/MPRAnalyze)
  MPRAnalyze provides statistical framework for the analysis of data
  generated by Massively Parallel Reporter Assays (MPRAs), used to
  directly measure enhancer activity. MPRAnalyze can be used for
  quantification of enhancer activity, classification of active
  enhancers and comparative analyses of enhancer activity between
  conditions. MPRAnalyze construct a nested pair of generalized
  linear models (GLMs) to relate the DNA and RNA observations, easily
  adjustable to various experimental designs and conditions, and
  provides a set of rigorous statistical testig schemes.

- [MSstatsTMT](https://bioconductor.org/packages/MSstatsTMT) Tools
  for protein significance analysis in shotgun mass
  spectrometry-based proteomic experiments with tandem mass tag (TMT)
  labeling.

- [MTseeker](https://bioconductor.org/packages/MTseeker) Variant
  analysis tools for mitochondrial genetics.

- [multiHiCcompare](https://bioconductor.org/packages/multiHiCcompare)
  multiHiCcompare provides functions for joint normalization and
  difference detection in multiple Hi-C datasets. This extension of
  the original HiCcompare package now allows for Hi-C experiments
  with more than 2 groups and multiple samples per group.
  multiHiCcompare operates on processed Hi-C data in the form of
  sparse upper triangular matrices. It accepts four column
  (chromosome, region1, region2, IF) tab-separated text files storing
  chromatin interaction matrices. multiHiCcompare provides cyclic
  loess and fast loess (fastlo) methods adapted to jointly
  normalizing Hi-C data. Additionally, it provides a general linear
  model (GLM) framework adapting the edgeR package to detect
  differences in Hi-C data in a distance dependent manner.

- [NBSplice](https://bioconductor.org/packages/NBSplice) The package
  proposes a differential splicing evaluation method based on isoform
  quantification. It applies generalized linear models with negative
  binomial distribution to infer changes in isoform relative
  expression.

- [NeighborNet](https://bioconductor.org/packages/NeighborNet)
  Identify the putative mechanism explaining the active interactions
  between genes in the investigated phenotype.

- [NormalyzerDE](https://bioconductor.org/packages/NormalyzerDE)
  NormalyzerDE provides screening of normalization methods for LC-MS
  based expression data. It calculates a range of normalized matrices
  using both existing approaches and a novel time-segmented approach,
  calculates performance measures and generates an evaluation report.
  Furthermore, it provides an easy utility for Limma- or ANOVA- based
  differential expression analysis.

- [nuCpos](https://bioconductor.org/packages/nuCpos) nuCpos, a
  derivative of NuPoP, is an R package for prediction of nucleosome
  positions. In nuCpos, a duration hidden Markov model is trained
  with a chemical map of nucleosomes either from budding yeast,
  fission yeast, or mouse embryonic stem cells. nuCpos outputs the
  Viterbi (most probable) path of nucleosome-linker states, predicted
  nucleosome occupancy scores and histone binding affinity (HBA)
  scores as NuPoP does. nuCpos can also calculate local and whole
  nucleosomal HBA scores for a given 147-bp sequence. Furthermore,
  effect of genetic alterations on nucleosome occupancy can be
  predicted with this package. The parental package NuPoP, which is
  based on an MNase-seq-based map of budding yeast nucleosomes, was
  developed by Ji-Ping Wang and Liqun Xi, licensed under GPL-2.

- [OMICsPCA](https://bioconductor.org/packages/OMICsPCA) OMICsPCA is
  an analysis pipeline designed to integrate multi OMICs experiments
  done on various subjects (e.g. Cell lines, individuals), treatments
  (e.g. disease/control) or time points and to analyse such
  integrated data from various various angles and perspectives. In
  it's core OMICsPCA uses Principal Component Analysis (PCA) to
  integrate multiomics experiments from various sources and thus has
  ability to over data insufficiency issues by using the ingegrated
  data as representatives. OMICsPCA can be used in various
  application including analysis of overall distribution of OMICs
  assays across various samples /individuals /time points; grouping
  assays by user-defined conditions; identification of source of
  variation, similarity/dissimilarity between assays, variables or
  individuals.

- [onlineFDR](https://bioconductor.org/packages/onlineFDR) This
  package allows users to control the false discovery rate for online
  hypothesis testing, where hypotheses arrive sequentially in a
  stream, as presented by Javanmard and Montanari (2015, 2018). In
  this framework, a null hypothesis is rejected based only on the
  previous decisions, as the future p-values and the number of
  hypotheses to be tested are unknown.

- [OUTRIDER](https://bioconductor.org/packages/OUTRIDER)
  Identification of aberrent gene expression in RNA-seq data. Read
  count expectations are modeled by an autoencoder to control for
  confounders in the data. Given these expectations, the RNA-seq read
  counts are assumed to follow a negative binomial distribution with
  a gene-specific dispersion. Outliers are then identified as read
  counts that significantly deviate from this distribution. Further
  OUTRIDER provides useful plotting function to analyze and visualize
  the results.

- [PepsNMR](https://bioconductor.org/packages/PepsNMR) This package
  provides R functions for common pre-procssing steps that are
  applied on 1H-NMR data. It also provides a function to read the FID
  signals directly in the Bruker format.

- [plotGrouper](https://bioconductor.org/packages/plotGrouper) A
  shiny app-based GUI wrapper for ggplot with built-in statistical
  analysis. Import data from file and use dropdown menus and
  checkboxes to specify the plotting variables, graph type, and look
  of your plots. Once created, plots can be saved independently or
  stored in a report that can be saved as a pdf. If new data are
  added to the file, the report can be refreshed to include new data.
  Statistical tests can be selected and added to the graphs. Analysis
  of flow cytometry data is especially integrated with plotGrouper.
  Count data can be transformed to return the absolute number of
  cells in a sample (this feature requires inclusion of the number of
  beads per sample and information about any dilution performed).

- [primirTSS](https://bioconductor.org/packages/primirTSS) A fast,
  convenient tool to identify the TSSs of miRNAs by integrating the
  data of H3K4me3 and Pol II as well as combining the conservation
  level and sequence feature, provided within both command-line and
  graphical interfaces, which achieves a better performance than the
  previous non-cell-specific methods on miRNA TSSs.

- [ProteoMM](https://bioconductor.org/packages/ProteoMM) ProteoMM is
  a statistical method to perform model-based peptide-level
  differential expression analysis of single or multiple datasets.
  For multiple datasets ProteoMM produces a single fold change and
  p-value for each protein across multiple datasets. ProteoMM
  provides functionality for normalization, missing value imputation
  and differential expression. Model-based peptide-level imputation
  and differential expression analysis component of package follows
  the analysis described in “A statistical framework for protein
  quantitation in bottom-up MS based proteomics" (Karpievitch et al.
  Bioinformatics 2009). EigenMS normalisation is implemented as
  described in "Normalization of peak intensities in bottom-up
  MS-based proteomics using singular value decomposition."
  (Karpievitch et al. Bioinformatics 2009).

- [qPLEXanalyzer](https://bioconductor.org/packages/qPLEXanalyzer)
  Tools for quantitative proteomics data analysis generated from
  qPLEX-RIME method.

- [QSutils](https://bioconductor.org/packages/QSutils) Set of utility
  functions for viral quasispecies analysis with NGS data. Most
  functions are equally useful for metagenomic studies. There are
  three main types: (1) data manipulation and exploration—functions
  useful for converting reads to haplotypes and frequencies,
  repairing reads, intersecting strand haplotypes, and visualizing
  haplotype alignments. (2) diversity indices—functions to compute
  diversity and entropy, in which incidence, abundance, and
  functional indices are considered. (3) data simulation—functions
  useful for generating random viral quasispecies data.

- [REBET](https://bioconductor.org/packages/REBET) There is an
  increasing focus to investigate the association between rare
  variants and diseases. The REBET package implements the
  subREgion-based BurdEn Test which is a powerful burden test that
  simultaneously identifies susceptibility loci and sub-regions.

- [Rmmquant](https://bioconductor.org/packages/Rmmquant) RNA-Seq is
  currently used routinely, and it provides accurate information on
  gene transcription. However, the method cannot accurately estimate
  duplicated genes expression. Several strategies have been
  previously used, but all of them provide biased results. With
  Rmmquant, if a read maps at different positions, the tool detects
  that the corresponding genes are duplicated; it merges the genes
  and creates a merged gene. The counts of ambiguous reads is then
  based on the input genes and the merged genes. Rmmquant is a
  drop-in replacement of the widely used tools findOverlaps and
  featureCounts that handles multi-mapping reads in an unabiased way.

- [RNASeqR](https://bioconductor.org/packages/RNASeqR) This R package
  is designed for case-control RNA-Seq analysis (two-group). There
  are six steps: "RNASeqRParam S4 Object Creation", "Environment
  Setup", "Quality Assessment", "Reads Alignment & Quantification",
  "Gene-level Differential Analyses" and "Functional Analyses". Each
  step corresponds to a function in this package. After running
  functions in order, a basic RNASeq analysis would be done easily.

- [SCBN](https://bioconductor.org/packages/SCBN) This package
  provides a scale based normalization (SCBN) method to identify
  genes with differential expression between different species. It
  takes into account the available knowledge of conserved orthologous
  genes and the hypothesis testing framework to detect differentially
  expressed orthologous genes. The method on this package are
  described in the article 'A statistical normalization method and
  differential expression analysis for RNA-seq data between different
  species' by Yan Zhou, Jiadi Zhu, Tiejun Tong, Junhui Wang, Bingqing
  Lin, Jun Zhang (2018, pending publication).

- [scruff](https://bioconductor.org/packages/scruff) A pipeline which
  processes single cell RNA-seq (scRNA-seq) reads from CEL-seq and
  CEL-seq2 protocols. Demultiplex scRNA-seq FASTQ files, align reads
  to reference genome using Rsubread, and generate UMI filtered count
  matrix. Also provide visualizations of read alignments and pre- and
  post-alignment QC metrics.

- [sesame](https://bioconductor.org/packages/sesame) Tools For
  analyzing Illumina Infinium DNA methylation arrays.

- [sigFeature](https://bioconductor.org/packages/sigFeature) This
  package provides a novel feature selection algorithm for binary
  classification using support vector machine recursive feature
  elimination SVM-RFE and t-statistic. In this feature selection
  process, the selected features are differentially significant
  between the two classes and also they are good classifier with
  higher degree of classification accuracy.

- [SIMD](https://bioconductor.org/packages/SIMD) This package
  provides a inferential analysis method for detecting differentially
  expressed CpG sites in MeDIP-seq data. It uses statistical
  framework and EM algorithm, to identify differentially expressed
  CpG sites. The methods on this package are described in the article
  'Methylation-level Inferences and Detection of Differential
  Methylation with Medip-seq Data' by Yan Zhou, Jiadi Zhu, Mingtao
  Zhao, Baoxue Zhang, Chunfu Jiang and Xiyan Yang (2018, pending
  publication).

- [slingshot](https://bioconductor.org/packages/slingshot) Provides
  functions for inferring continuous, branching lineage structures in
  low-dimensional data. Slingshot was designed to model developmental
  trajectories in single-cell RNA sequencing data and serve as a
  component in an analysis pipeline after dimensionality reduction
  and clustering. It is flexible enough to handle arbitrarily many
  branching events and allows for the incorporation of prior
  knowledge through supervised graph construction.

- [slinky](https://bioconductor.org/packages/slinky) Wrappers to
  query the L1000 metadata available via the clue.io REST API as well
  as helpers for dealing with LINCS gctx files, extracting data sets
  of interest, converting to SummarizedExperiment objects, and some
  facilities for performing streamlined differential expression
  analysis of these data sets.

- [sparsenetgls](https://bioconductor.org/packages/sparsenetgls) The
  package provides methods of combining the graph structure learning
  and generalized least squares regression to improve the regression
  estimation. The main function sparsenetgls() provides solutions for
  multivariate regression with Gaussian distributed dependant
  variables and explanatory variables utlizing multiple well-known
  graph structure learning approaches to estimating the precision
  matrix, and uses a penalized variance covariance matrix with a
  distance tuning parameter of the graph structure in deriving the
  sandwich estimators in generalized least squares (gls) regression.
  This package also provides functions for assessing a Gaussian
  graphical model which uses the penalized approach. It uses Receiver
  Operative Characteristics curve as a visualization tool in the
  assessment.

- [strandCheckR](https://bioconductor.org/packages/strandCheckR) This
  package aims to quantify and remove putative double strand DNA from
  a strand-specific RNA sample. There are also options and methods to
  plot the positive/negative proportions of all sliding windows,
  which allow users to have an idea of how much the sample was
  contaminated and the appropriate threshold to be used for
  filtering.

- [TimeSeriesExperiment](https://bioconductor.org/packages/TimeSeriesExperiment)
  Visualization and analysis toolbox for short time course data which
  includes dimensionality reduction, clustering, two-sample
  differential expression testing and gene ranking techniques. The
  package also provides methods for retrieving enriched pathways.

- [transite](https://bioconductor.org/packages/transite) transite is
  a computational method that allows comprehensive analysis of the
  regulatory role of RNA-binding proteins in various cellular
  processes by leveraging preexisting gene expression data and
  current knowledge of binding preferences of RNA-binding proteins.

- [tRNA](https://bioconductor.org/packages/tRNA) The tRNA package
  allows tRNA sequences and structures to be accessed and used for
  subsetting. In addition, it provides visualization tools to compare
  feature parameters of multiple tRNA sets and correlate them to
  additional data. The tRNA package uses GRanges objects as inputs
  requiring only few additional column data sets.

- [tRNAdbImport](https://bioconductor.org/packages/tRNAdbImport)
  tRNAdbImport imports the entries of the tRNAdb and mtRNAdb
  (http://trna.bioinf.uni-leipzig.de) as GRanges object.

- [tximeta](https://bioconductor.org/packages/tximeta) Transcript
  quantification import from Salmon with automatic population of
  metadata and transcript ranges. Filtered, combined, or de novo
  transcriptomes can be linked to the appropriate sources with
  linkedTxomes and shared for reproducible analyses.

- [Ularcirc](https://bioconductor.org/packages/Ularcirc) Ularcirc
  reads in STAR aligned splice junction files and provides
  visualisation and analysis tools for splicing analysis. Users can
  assess backsplice junctions and forward canonical junctions.

- [universalmotif](https://bioconductor.org/packages/universalmotif)
  Allows for importing most common motif types into R for use by
  functions provided by other Bioconductor motif-related packages.
  Motifs can be exported into most major motif formats from various
  classes as defined by other Bioconductor packages. A suite of motif
  and sequence manipulation and analysis functions are included,
  including enrichment, comparison, P-value calculation, shuffling,
  trimming, higher-order motifs, and others.

- [Wrench](https://bioconductor.org/packages/Wrench) Wrench is a
  package for normalization sparse genomic count data, like that
  arising from 16s metagenomic surveys.

- [XINA](https://bioconductor.org/packages/XINA) An intuitive R
  package simplifies network analyses output from multiplexed
  high-dimensional proteomics/trascriptomics kinetics data.


New Data Experiment Packages
=====================

There are 21 new data experiment packages in this release of Bioconductor.

- [allenpvc](https://bioconductor.org/packages/allenpvc) Celular
  taxonomy of the primary visual cortex in adult mice based on single
  cell RNA-sequencing from a study performed by the Allen Institute
  for Brain Science. In said study 49 transcriptomic cell types are
  identified.

- [AssessORFData](https://bioconductor.org/packages/AssessORFData)
  This package provides access to mapping and results objects
  generated by the AssessORF package, as well as the genome sequences
  for the strains corresponding to those objects.

- [brainImageRdata](https://bioconductor.org/packages/brainImageRdata)
  brainImageRdata contains image masks for the developing human and
  the adult human brain. These masks can be used in conjunction with
  the gene expression data to generate spatial gene set enrichment
  plots. It also contains the expression data for the 15 pcw human
  brain, the adult human brain, and the developing human brain.

- [breakpointRdata](https://bioconductor.org/packages/breakpointRdata)
  Strand-seq data to demonstrate functionalities of breakpointR
  package.

- [celarefData](https://bioconductor.org/packages/celarefData) This
  experiment data contains some processed data used in the celaref
  package vignette. These are publically available datasets, that
  have been processed by celaref package, and can be manipulated
  further with it.

- [CopyNeutralIMA](https://bioconductor.org/packages/CopyNeutralIMA)
  Provides a set of genomic copy neutral samples hybridized using
  Illumina Methylation arrays (450k and EPIC).

- [DuoClustering2018](https://bioconductor.org/packages/DuoClustering2018)
  Preprocessed experimental and simulated scRNA-seq data sets used
  for evaluation of clustering methods for scRNA-seq data in Duò et
  al (2018). Also contains results from applying several clustering
  methods to each of the data sets, and functions for plotting method
  performance.

- [FlowSorted.Blood.EPIC](https://bioconductor.org/packages/FlowSorted.Blood.EPIC)
  Raw data objects to be used for blood cell proportion estimation in
  minfi and similar packages. The FlowSorted.Blood.EPIC object is
  based in samples assayed by Brock Christensen and colleagues; for
  details see Salas et al. 2018.
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554.

- [GIGSEAdata](https://bioconductor.org/packages/GIGSEAdata) The gene
  set collection used for the GIGSEA package.

- [mcsurvdata](https://bioconductor.org/packages/mcsurvdata) This
  package stores two merged expressionSet objects that contain the
  gene expression profile and clinical information of -a- six breast
  cancer cohorts and -b- four colorectal cancer cohorts. Breast
  cancer data are employed in the vignette of the hrunbiased package
  for survival analysis of gene signatures.

- [MSMB](https://bioconductor.org/packages/MSMB) Data sets for the
  book 'Modern Statistics for Modern Biology', S.P. Holmes and W.
  Huber

- [MTseekerData](https://bioconductor.org/packages/MTseekerData)
  Provides examples for the MTseeker package vignette.

- [OMICsPCAdata](https://bioconductor.org/packages/OMICsPCAdata)
  Supporting data for package OMICsPCA

- [PepsNMRData](https://bioconductor.org/packages/PepsNMRData) This
  package contains all the datasets used in the PepsNMR package.

- [qPLEXdata](https://bioconductor.org/packages/qPLEXdata) qPLEX-RIME
  and Full proteome TMT mass spectrometry datasets.

- [RegParallel](https://bioconductor.org/packages/RegParallel) In
  many analyses, a large amount of variables have to be tested
  independently against the trait/endpoint of interest, and also
  adjusted for covariates and confounding factors at the same time.
  The major bottleneck in these is the amount of time that it takes
  to complete these analyses. With RegParallel, a large number of
  tests can be performed simultaneously. On a 12-core system, 144
  variables can be tested simultaneously, with 1000s of variables
  processed in a matter of seconds via 'nested' parallel processing.
  Works for logistic regression, linear regression, conditional
  logistic regression, Cox proportional hazards and survival models,
  Bayesian logistic regression, and negative binomial regression.

- [RNASeqRData](https://bioconductor.org/packages/RNASeqRData)
  RNASeqRData is a helper experiment package for vignette
  demonstration purpose in RNASeqR software package.

- [sesameData](https://bioconductor.org/packages/sesameData) Provides
  supporting annotation and test data for SeSAMe package.

- [TabulaMurisData](https://bioconductor.org/packages/TabulaMurisData)
  Access to processed 10x (droplet) and SmartSeq2 (on FACS-sorted
  cells) single-cell RNA-seq data from the Tabula Muris consortium
  (http://tabula-muris.ds.czbiohub.org/).

- [tcgaWGBSData.hg19](https://bioconductor.org/packages/tcgaWGBSData.hg19)
  Data package for WGBS Data in TCGA. Data is stored as
  SummarizedExperiment Format. See vignette on how to extract the
  data and perform differential methylation analysis.

- [TENxPBMCData](https://bioconductor.org/packages/TENxPBMCData)
  Single-cell RNA-seq data for on PBMC cells, generated by 10X
  Genomics.

New Workflows
=====================

There are 2 new workflow packages in this release of Bioconductor.

- [maEndToEnd](https://bioconductor.org/packages/maEndToEnd) In this
  article, we walk through an end-to-end Affymetrix microarray
  differential expression workflow using Bioconductor packages. This
  workflow is directly applicable to current "Gene" type arrays, e.g.
  the HuGene or MoGene arrays, but can easily be adapted to similar
  platforms. The data analyzed here is a typical clinical microarray
  data set that compares inflamed and non-inflamed colon tissue in
  two disease subtypes. For each disease, the differential gene
  expression between inflamed- and non-inflamed colon tissue was
  analyzed. We will start from the raw data CEL files, show how to
  import them into a Bioconductor ExpressionSet, perform quality
  control and normalization and finally differential gene expression
  (DE) analysis, followed by some enrichment analysis.

- [rnaseqDTU](https://bioconductor.org/packages/rnaseqDTU) RNA-seq
  workflow for differential transcript usage (DTU) following Salmon
  quantification. This workflow uses Bioconductor packages tximport,
  DRIMSeq, and DEXSeq to perform a DTU analysis on simulated data. It
  also shows how to use stageR to perform two-stage testing of DTU, a
  statistical framework to screen at the gene level and then confirm
  which transcripts within the significant genes show evidence of
  DTU.


NEWS from new and existing Software Packages
===================================


[ABAEnrichment](https://bioconductor.org/packages/ABAEnrichment)
-------------

Changes in version 1.11.7:

USER-LEVEL CHANGES

- simplify functions 'get_expression' and 'plot_expression' (remove
  option to automatically use data from last aba_enrich-call)

- 'plot_expression' now takes matrix from 'get_expression' as input,
  instead of calling 'get_expression' internally

- add color key to 'plot_expression' heatmap when gene-associated
  variables are shown in a colored side bar

Changes in version 1.11.6:

USER-LEVEL CHANGES

- sort output aba_enrich()&#91;&#91;1&#93;&#93; by FWER first and then by age category
  (previously it was sorted by age category first)

Changes in version 1.11.3:

- remove all C++ code and depend on package GOfuncR instead

Changes in version 1.11.2:

USER-LEVEL CHANGES

- when coordinates are used, remove genes with multiple
  gene-coordinates (coordinates are used when the familiy-wise error
  rate is corrected for gene length or spatial clustering of genes)

Changes in version 1.11.1:

NEW FEATURES

- allow for custom gene-coordinates (alternative to integrated
  coordinates; coordinates are used when the familiy-wise error rate is
  corrected for gene length or spatial clustering of genes)

[abseqR](https://bioconductor.org/packages/abseqR)
------

Changes in version 0.99.0:

- Submitted to Bioconductor

[ABSSeq](https://bioconductor.org/packages/ABSSeq)
------

Changes in version 1.35.2 (2018-06-28):

- Updating! Extend aFold to paired samples

[ACE](https://bioconductor.org/packages/ACE)
---

Changes in version 0.99.8 (2018-09-06):

- bugfix in squaremodelsummary: remove standard in singleplot call

Changes in version 0.99.6 (2018-07-26):

- bugfix in linkvariants: no default freqindex

Changes in version 0.99.4 (2018-06-26):

IMPROVEMENTS

- changed linkmutationdata to linkvariants

- linkvariants can estimate copies of heterozygous germline variants

- linkvariants calculates the upper and lower bound of a given
  confidence level when read depth is given

- changed postanalysisloop to accommodate linkvariants

- added argument in runACE to specify genome build

- several minor documentation amendments

- improved coding robustness using seq and seq_along functions

Changes in version 0.99.0 (2018-05-24):

SUBMISSION

- initial submission to Bioconductor

[ADaCGH2](https://bioconductor.org/packages/ADaCGH2)
-------

Changes in version 2.21.1 (2018-07-16):

- Removed all mentions to old installer

[affxparser](https://bioconductor.org/packages/affxparser)
----------

Changes in version 1.53.2 (2018-10-22):

DOCUMENTATION

- Link to Affx Fusion SDK archive on GitHub.

- Spell corrections.

Changes in version 1.53.1 (2018-08-28):

- Updated installation instructions in README.md.

Changes in version 1.53.0 (2018-04-30):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC 3.8 for R (>= 3.5.0).

[AnnotationFilter](https://bioconductor.org/packages/AnnotationFilter)
----------------

Changes in version 1.5.2:

USER VISIBLE CHANGES

- Rename GenenameFilter into GeneNameFilter and deprecate
  GenenameFilter (issue #22).

[AnnotationHubData](https://bioconductor.org/packages/AnnotationHubData)
-----------------

Changes in version 1.11.0:

MODIFICATIONS

- Removed scripts for Pazar DB as website no longer active

- Update from BiocInstaller to BiocManager

NEW FEATURES

- Species and taxonomyId are now validated against GenomeInfoDbData
  object

BUG FIX

- Fix TwoBit resource receipe. Converts DNA that is not A,C,T,G,N to N
  do to design of rtracklayer::export for TwoBit

- Fix bug with assignment of tags in annotationhub

- makeEpigenomeRoadMap recipe updated to account for XML bug that
  cannot handle http urls. updated to https

[anota2seq](https://bioconductor.org/packages/anota2seq)
---------

Changes in version 1.3.1:

- fixed a bug where in anota2seqPlotPvalues and anota2seqPlotFC the
  contrast names were not displayed correctly when selecting only 1
  contrast in case there were multiple

- Using the anota2seqRun function with custom filtering parameters for
  maxP was still based on a maxPAdj of 0.15. This has been fixed, i.e.
  when maxP filtering is applied no maxPAdj filter will be used.

[apeglm](https://bioconductor.org/packages/apeglm)
------

Changes in version 1.3.2:

- Fix simple bug when no sample has only intercept coefficient.

[aroma.light](https://bioconductor.org/packages/aroma.light)
-----------

Changes in version 3.11.2 (2018-09-04):

CODE REFACTORING

- fitPrincipalCurve() now requires princurve (>= 2.1.2) and was updated
  to make use of new principcal_curve class instead of deprecated
  principcal.curve class.  This update "should not" affect the results,
  but see https://github.com/dynverse/princurve/issues/8 for
  information of what has changed in the princurve package in this
  respect.

Changes in version 3.11.1 (2017-08-28):

- Updated installation instructions in README.md.

Changes in version 3.11.0 (2017-04-30):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC 3.8 for R (>= 3.6.0).

[artMS](https://bioconductor.org/packages/artMS)
-----

Changes in version 0.99.25:

- Update Vignette to just output html

Changes in version 0.99.02:

- Add the Bioconductor webhook

Changes in version 0.99.01:

- Submit the package to the Bioconductor project For Developers

- All the functions name must have the prefix 'artms'

- Normalized all documentation using roxygen2 Deprecated

- Nothing yet Defunct

- Nothing yet

[ASICS](https://bioconductor.org/packages/ASICS)
-----

Changes in version 1.1.1:

Improvements

- need to specify number of cores to use a parallel environment

- reference library (improvement of reference spectra by cleaning them)

[ASSIGN](https://bioconductor.org/packages/ASSIGN)
------

Changes in version 1.17.2:

- Bugfix ComBat.step2 by forcing ComBat input to be a matrix

- Added new argument to assign.wrapper to specify direction of
  signature genes.

- parameters file in assign.wrapper is now in yaml format.

[attract](https://bioconductor.org/packages/attract)
-------

Changes in version 1.33.2 (2018-06-29):

- Fix NEWS

- change maintainer email

Changes in version 1.33.1 (2018-05-15):

- Fix NEWS

[BASiCS](https://bioconductor.org/packages/BASiCS)
------

Changes in version 1.3.24 (2018-10-21):

- Version bump to trigger new build

Changes in version 1.3.23 (2018-10-18):

- Updated vignette to include simple SCE example without spike-in
  information

Changes in version 1.3.22 (2018-10-19):

- Update of DESCRIPTION to require C++11

- Once again, replacement of `R::lgammafn` by `std::lgamma` in c++ code

Changes in version 1.3.21 (2018-10-19):

- Restores use of `R::lgammafn` due to errors in Bioconductor build
  report

Changes in version 1.3.20 (2018-10-18):

- Minor updates to README and vignette (e.g. to include BASiCS sticker)

Changes in version 1.3.19 (2018-10-18):

- Replacement of `R::lgammafn` by `std::lgamma` in c++ code (great
  speed-up! with thanks to Shian Su (@Shians))

Changes in version 1.3.18 (2018-10-15):

- `BASiCS_D_TestDE` has been deprecated (replaced by `BASiCS_TestDE`)

- `SingleCellExperiment::` added to all calls of `isSpike`

- Extra input check in `BASiCS_MCMC` to avoid issues due to multiple
  types of spike-ins being present in the data.

Changes in version 1.3.17 (2018-10-11):

- New unit tests for different argument settings for the BASiCS_MCMC
  function

Changes in version 1.3.16 (2018-10-08):

- Reduce checks for Data that was not generated with the newBASiCS_Data
  function

Changes in version 1.3.15 (2018-10-07):

- Fixed bug when storing chains in the non-spikes case

Changes in version 1.3.13 (2018-09-30):

- Some Rd files generated again

Changes in version 1.3.11 (2018-09-30):

- Alan added as contributor

Changes in version 1.3.10 (2018-09-28):

- Prevent BASiCS_Chain from crashing when StoreAdapt=TRUE for the
  regression case

- Prevent plot function errors

- Updated documentation using roxygen2

Changes in version 1.3.9 (2018-09-27):

- Changed raw data accessor to 'counts' instead of 'assay'

- makeExampleBASiCS_Data now generates 30 cells instead of 20

- Updated unit tests to match these 30 cells

Changes in version 1.3.8 (2018-09-05):

- Updated README.md (Travis/codecov badges) + updated citation to Cell
  Systems

Changes in version 1.3.7 (2018-09-05):

- Commit to activate Travis/codecov

Changes in version 1.3.6 (2018-09-04):

- Updated CITATION file

Changes in version 1.3.5 (2018-08-31):

- Updated reference to latest paper

Changes in version 1.3.3 (2018-05-21):

- Missing line added to test_data_examples.R

Changes in version 1.3.2 (2018-05-21):

- Version bump only

Changes in version 1.3.1 (2018-05-21):

- Unit tests created for selected C++ functions (`Hidden_rDirichlet`
  and `Hidden_muUpdate`)

- Updated unit tests in spike-in parameter estimation to have tolerance
  (TEMPORARY FIX)

- Minor updates in `makeExampleBASiCS_Data` and `newBASiCS_Data`
  function to allow different spike-in types (as in
  `SingleCellExperiment`)

[bayNorm](https://bioconductor.org/packages/bayNorm)
-------

Changes in version 0.99.19:

- Replacing foreach with BiocParallel

- Fix some errors in BB_fun

Changes in version 0.99.9:

- Let more examples to be runnable; Updated vignette; Including
  SingleCellExperiment Class.

Changes in version 0.99.8:

- Fix some errors.

Changes in version 0.99.4:

- Apply corrections proposed by the reviewer, for more details please
  see <URL: https://github.com/Bioconductor/Contributions/issues/878>

Changes in version 0.99.1:

- Fix two warnings reported by Bioconductor * WARNING: Use TRUE/FALSE
  instead of T/F, Found in files: tests/testthat/test-bayNorm.r *
  WARNING: Use is() or !is()

Changes in version 0.99.0:

- Pass all checks of both R CMD build and R CMD BiocCheck.

- Future work: * Vignette needs to be further improved.  * Need to
  Improve man pages.  * Submit to Bioconductor.

[bcSeq](https://bioconductor.org/packages/bcSeq)
-----

Changes in version 1.3.3:

- update format for NEWS

Changes in version 1.3.2:

- update citation and fix bugs

Changes in version 1.3.1:

- update citations

[beachmat](https://bioconductor.org/packages/beachmat)
--------

Changes in version 1.4.0:

- Removed native support for RleMatrix and packed symmetric matrices.

- Added multi-row/column getters.

- Added mechanism for native support of arbitrary developer-defined
  matrices.

- Switched to row/colGrid() for defining chunks in unsupported
  matrices.

[BEclear](https://bioconductor.org/packages/BEclear)
-------

Changes in version 1.13.11 (2018-10-25):

- implementation of further tests

Changes in version 1.13.10 (2018-10-24):

- test for the gdepoch function

- simplification of code

- performance improvements

Changes in version 1.13.9 (2018-10-23):

- test for the loss function

- some simplifications

Changes in version 1.13.8 (2018-10-22):

- first implementation of tests with testthat

- simplifications of the code

Changes in version 1.13.7 (2018-10-10):

- small bugfix

Changes in version 1.13.6 (2018-10-05):

- some minor issues with the calculation of p-values are fixed

Changes in version 1.13.5 (2018-10-04):

- checks for the validity of inputs added

- some performance improvements

Changes in version 1.13.4 (2018-10-02):

- major performance improvement

Changes in version 1.13.3 (2018-09-28):

- usage of data.table for various functions for performance improvement

- simplification of the source code

Changes in version 1.13.2 (2018-09-27):

- bug fixed that, where p-values couldn't be calculated for a gene with
  only missing values in a batch (#1)

Changes in version 1.13.1 (2018-09-25):

- BEclear uses now BiocParallel instead of snowfall for parallelisation

- roxygen2 is now used for the generating the documentation

- major code refactoring

- some minor bug fixes

- performance improvements

[BgeeDB](https://bioconductor.org/packages/BgeeDB)
------

Changes in version 2.6.2:

- Fix issue in the formatData() function. It is now possible to format
  using fpkm expression values when using Bgee 14.0.

- Implementation of regression tests

- Update vignette and README

Changes in version 2.6.1:

- Fix issue on Bgee 14.0 tar.gz annotation file management.

- Update README and DESCRIPTION files.

[BiFET](https://bioconductor.org/packages/BiFET)
-----

Changes in version 1.1.5 (2018-07-02):

- Fixed R documentation notes # Updated date field in DESCRIPTION

[bioCancer](https://bioconductor.org/packages/bioCancer)
---------

Changes in version 1.9.07:

- replace DOSE::dotplot by clusterProfiler::dotplot

- replace r_data by r_info in reports

- not need to define genelist in r_data

- update ReactomeFI.RDS file (version 2017)

- update DisGeNet0918.RDS file (version September 2018), move it from
  wiki.ubuntu.com to github/kmezhoud

Changes in version 1.9.06:

- rm warning message for min and max functions

- run getFreqMutData() when getListProfData() instead
  getCoffeeWheel_Mut(). Avoid error when loading x profiles data to
  workspace.

- include Tools panel into Workspace panel.

- update Overview image

Changes in version 1.9.05:

- r_data vs r_info ... https://radiant-rstats.github.io/docs/news.html

- use r_info for dataset list and r_data for genes list

- set progress bar

Changes in version 1.9.04:

- add style.css file: raduce padding-top to 0px

Changes in version 1.9.03:

- modify stop function

- update paste gene list function

Changes in version 1.9.02:

- upload and download using Rstudio file browser remove plot_downloader
  function and replace it by download_link (defined in radiant.R file)

- add radiant_old.R file for needed functions but not longer used by
  radiant.data

Changes in version 1.9.01:

- replace getdata() by get_data()

- replace factorizer() by lapply(.,factor)

[BiocCheck](https://bioconductor.org/packages/BiocCheck)
---------

Changes in version 1.17:

NEW FEATURES

- (1.17.21) Added quit-with-status option to both BiocCheck and
  BiocCheckGitClone for compatibility with travis

- (1.17.18) Update devel to use BiocManager instructions instead of
  BiocInstaller

- (1.17.17) Add a new function that can be run interactive or command
  line BiocCheckGitClone which is only run on a source directory not a
  tarball. This will check for bad system files

- (1.17.17) BiocCheck addition: Checks vignette directory for
  intermediate and end files that should not be included.

- (1.17.16) Checks for Bioconductor package size requirement if
  checking tarball

BUG FIXES

- (1.17.19) Updated internal functions to use BiocManger instead of
  BiocInstaller

[BiocFileCache](https://bioconductor.org/packages/BiocFileCache)
-------------

Changes in version 1.5:

NEW FEATURES

- (v. 1.5.1) Add 'exact = ' for exact matching in bfcquery(),
  bfcrpath(). (v. 1.5.2) defaults to TRUE for bfcrpath()

BUG FIX

- (v. 1.5.2) bfcrpath() more robust when adding regular expression
  rnames.

USER-VISIBLE CHANGES

- (v. 1.5.4) bfcpath() implementation change. Only displays rpath and
  can work with multiple rids. bfcpath only for rpath access while
  bfcrpath is the option to get or add.

[BiocInstaller](https://bioconductor.org/packages/BiocInstaller)
-------------

Changes in version 1.31.2:

NEW FEATURES

- 'BiocInstaller' is currently deprecated. All installation of
  Bioconductor packages will be done via `BiocManager` a CRAN package.

BUG FIXES

- The documentation incorrectly refered to `remotes::install` which is
  not an exported function of the package. This link was fixed.

[BiocNeighbors](https://bioconductor.org/packages/BiocNeighbors)
-------------

Changes in version 1.0.0:

- New package kmknn, for k-means-based k-nearest neighbor detection.

[BioCor](https://bioconductor.org/packages/BioCor)
------

Changes in version 1.5.4:

- Improved spelling

- Reduced the complexity fo combineScoresPar and combineScores.

- Improve efficiency of code in vignettes.

- Change the news section

- Adding a section about GeneOverlap package.

- Changed the License

[BiocParallel](https://bioconductor.org/packages/BiocParallel)
------------

Changes in version 1.16:

NEW FEATURES

- (v 1.15.9) BatchtoolsParam() gains resources=list() for template file
  substitution.

- (v 1.15.12) bpexportglobals() for all BPPARAM exports global options
  (i.e., base::options()) to workers. Default TRUE.

BUG FIXES

- (v 1.15.6) bpiterate,serial-method does not return a list() when
  REDUCE present
  (https://github.com/Bioconductor/BiocParallel/issues/77)

- (v 1.15.7) bpaggregate,formula-method failed to find BPREDO
  (https://support.bioconductor.org/p/110784)

- (v 1.15.13) bplappy,BatchtoolsParam() coerces List to list
  (https://github.com/Bioconductor/BiocParallel/issues/82)

- (v 1.15.14) implicit loading of BiocParallel when loading a third-
  party package failed because reference class `initialize()` methods
  are not installed correctly. This bug fix results in signficant
  revision in the implementation, so that valid objects must be
  constructed through the public constructors, e.g.,
  `BatchtoolsParam()`

[breakpointR](https://bioconductor.org/packages/breakpointR)
-----------

Changes in version 0.99.0:

- This package implements various functions to find template strand
  changepoints in Strand-seq data.

[bsseq](https://bioconductor.org/packages/bsseq)
-----

Changes in version 1.17:

- BSseq() will no longer reorder inputs. Previously, the returned
  _BSseq_ object was ordered by ordering the loci, although this
  behaviour was not documented. BSseq() may still filter out loci if
  rmZeroCov = FALSE or collapse loci if strandCollapse = FALSE or
  duplicate loci are detected, but the relative order of loci in the
  output will match that of the input.

- Fix bug with maxGap argument of BSmooth(). The bug meant that the
  'maximum gap between two methylation loci' was incorrectly set to 2 *
  maxGap + 1 instead of maxGap. This likely did not affect results for
  users who left the default value of maxGap = 10^8 but may have
  affected results for small values of maxGap.

[BUScorrect](https://bioconductor.org/packages/BUScorrect)
----------

Changes in version 0.99.13 (2018-10-09):

- Improved the ``postprob\_DE\_thr\_fun'' function.

Changes in version 0.99.6 (2018-07-27):

- Created a NEWS.Rd file parsable by utils::news within the inst folder

[CAGEfightR](https://bioconductor.org/packages/CAGEfightR)
----------

Changes in version 1.5:

- Added new functions for spatial analysis of clusters: findLinks finds
  nearby pairs of clusters (for example TSSs and enhancers) and
  calculates the correlation of expression between them. findStretches
  find stretches along the genome where clusters are within a certain
  distance of eachother (for example groups of enhancer forming a super
  enhancer) and calculates the average pairwise correlation between
  members.

- Changed the way clustering works: CAGEfightR uses coverage() to
  calculate genome-wide signals and now rounds the resulting signal to
  a certain number of digits (this can be modified via the
  CAGEfightR.round option), to prevent small positive or negative
  values due to floating point errors. This makes clustering more
  stable meaning the tuneTagClustering function is now deprecated. This
  should also increase the speed of most functions.

- CAGEfightR now uses GPos instead of GRanges for storing CTSSs, this
  should result in improved memory performance.

- Several changes to clusterBidirectionality: Balance is now calculated
  using the midpoint as well (preventing some rare cases where the
  midpoint could mask a single highly expressed CTSS), the pooled CTSS
  signal is now prefiltered for bidirectionality to increase speed and
  custom balance function can be provided (Bhattacharyya coefficient
  and Andersson's D are included).

- Added new check-functions to make it easier to check if objects are
  formatted correctly

[CAGEr](https://bioconductor.org/packages/CAGEr)
-----

Changes in version 1.24.0:

NEW FEATURES

- plotAnnot(): the "group" argument now accepts formulas such as "~ a +
  b" to indicate names of metadata column that will be pasted together
  to form a new group factor.

BUG FIXES

- Prevent mergeSamples() from producing colData that cause other
  functions to crash later when coercing to data.frame.

- Repaired paraclu support for CAGEset objects.

- normalizeTagCount() works again on CAGEset objects.

- consensusClustersGR() reports expression score of the selected sample
  (instead of silently ignoring the "sample" argument and reporting
  expression sum on all the samples).

[Cardinal](https://bioconductor.org/packages/Cardinal)
--------

Changes in version 1.99.2 (2018-10-28):

SIGNIFICANT USER-VISIBLE CHANGES

- Updated vignettes for Cardinal 2.0

BUG FIXES

- Removed a unit test broken on Windows

Changes in version 1.99.1 (2018-10-26):

NEW FEATURES

- Added vignettes and documentation for Cardinal 2.0

Changes in version 1.99.0 (2018-10-25):

SIGNIFICANT USER-VISIBLE CHANGES

- Version bump for Cardinal v2 release candidate

Changes in version 1.13.3 (2018-10-24):

NEW FEATURES

- Added 'process' method for queueing delayed processing functions to
  an imaging dataset and applying them

- Added new processing methods for Cardinal v2 including new versions
  of 'normalize', 'smoothSignal', 'reduceBaseline', 'peakPick',
  'peakAlign', and 'peakFilter'

- Added new 'peakBin' function for binning peaks

- Updated 'show' method for new Cardinal v2 classes

- New support for exporting 'processed' imzML files via the
  'writeImzML' function

Changes in version 1.13.2:

NEW FEATURES

- Added new classes for Cardinal v2 including 'XDataFrame',
  'PositionDataFrame', 'MassDataFrame', 'ImagingExperiment',
  'SparseImagingExperiment', and 'MSImagingExperiment'

Changes in version 1.13.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated installation instructions for "CardinalWorkflows"

Changes in version 1.12.1:

BUG FIXES

- Fixed bug in reading Analyze 7.5 files

[cbaf](https://bioconductor.org/packages/cbaf)
----

Changes in version 1.4.0 (2018-10-30):

New Features

- obtainOneStudy() and obtainMultipleStudies() functions can obtain
  data for groups of genes each possess more than 250 genes (Virtually
  unlimited gene number).

- A new argument for xlsxOutput() function to exchange the columns and
  rows.

[ccfindR](https://bioconductor.org/packages/ccfindR)
-------

Changes in version 1.1.2:

Changes

- Added C++ code update step in vb_factorize(..., useC=TRUE)

- Added Singular value decomposition initializer for vb_factorize(...,
  initializer='svd2') random initial condition: initializer='random'

- Changed default: filter_genes(...,rescue.genes=FALSE)

- Changed filter_genes(), vmr.min action from vmr >= vmr.min to vmr >
  vmr.min (removes genes with vmr=0)

- Added parallel run for vb_factorize(..., ncores=10)

- Added feature_map(...)

[celaref](https://bioconductor.org/packages/celaref)
-------

Changes in version 0.99.1:

BUG FIXES

- Code style.

- Do not attempt to multithread on windows (suggest mulithread on
  linux).

Changes in version 0.99.0:

NEW FEATURES

- Initial Version.

[CellTrails](https://bioconductor.org/packages/CellTrails)
----------

Changes in version 0.99.15:

- Bugfixes: - Adjusted NAMESPACE for most recent 'SingleCellExperiment'
  release

Changes in version 0.99.14:

- Bugfixes: - Filtering. Captures invalid user input. - Clustering.
  Parameter settings resulting in a minimum number of samples per
  initial state <= 1 are captured.

Changes in version 0.99.13:

- Bugfixes: - Visualization. Package is compatible with most recent
  'ggplot2' release (3.x.x).

- Updated CITATION

Changes in version 0.99.0:

- Public pre-release for Bioconductor submission

- CellTrails is fully compatible with a 'SingleCellExperiment' object

- Bugfixes

[chimeraviz](https://bioconductor.org/packages/chimeraviz)
----------

Changes in version 1.7.1:

CHANGES

- Introduced support for oncofuse. Ref
  https://github.com/stianlagstad/chimeraviz/issues/47 and
  https://github.com/stianlagstad/chimeraviz/pull/48/files.


[ChIPpeakAnno](https://bioconductor.org/packages/ChIPpeakAnno)
------------

Changes in version 3.15.2:

- fix the bug in binOverRegions and binOverGene

Changes in version 3.15.1:

- Fix the bug in assignChromosomeRegion for the Jaccard index
  calculation.

[chromstaR](https://bioconductor.org/packages/chromstaR)
---------

Changes in version 1.7.2:

BUGFIXES

- Compatibility fixes for the new release of ggplot2 (3.0.0).

- seqlevels() that are smaller than binsize are dropped properly in
  fixedWidthBins() and variableWidthBins().

[CHRONOS](https://bioconductor.org/packages/CHRONOS)
-------

Changes in version 1.9.2:

- Reformatted NEWS file.

Changes in version 1.9.1:

- Reconfigured KEGG pathway downloading via KEGG API.

[ClassifyR](https://bioconductor.org/packages/ClassifyR)
---------

Changes in version 2.2.0:

- getClasses is no longer a slot of PredictParams. Every predictor
  function needs to return either a factor vector of classes, a numeric
  vector of class scores for the second class, or a data frame with a
  column for the predicted classes and another for the second-class
  scores.

- Cross-validations which use folds ensure that samples belonging to
  each class are in approximately the same proportions as they are for
  the entire data set.

- Classification can reuse fitted model from previous classification by
  using previousTrained function.

- Feature selection using gene sets and networks. Classification can
  use meta-features derived from the individual features used for
  feature selection.

- tTestSelection function for feature selection based on ordinary
  t-test statistic ranking. Now the default feature selection function,
  if none is specified.

- Tuning parameter optimisation metric is specified by providing a
  tuneOptimise parameter to TrainParams rather than depending on
  ResubstituteParams being used during feature selection.

[clusterExperiment](https://bioconductor.org/packages/clusterExperiment)
-----------------

Changes in version 2.1.5 (2018-06-28):

Changes

- Add functionality to `getBestFeatures` to allow `edgeR` for DE, as
  well as weights used with `edgeR` for `zinbwave` compatability. As
  part of this change: - Removed `isCount` argument and replaced with
  more fine-grained `DEMethod` argument in `getBestFeatures`,
  `mergeClusters`; or `mergeDEMethod` in `RSEC`. - *Change to class
  definition*: added slot `merge_demethod` to keep track of the DE
  method used in merging

- Change function names (old function name is now depricated): -
  `combineMany` -> `makeConsensus` - `removeUnclustered` ->
  `removeUnassigned`

- Arguments to functions changed: - `combineProportion` ->
  `consensusProportion` in `RSEC` - `combineMinSize` ->
  `consensusMinSize` in `RSEC` - `sampleData` -> `colData` to match
  `SummarizedExperiment` syntax (in many plotting functions).  -
  `alignSampleData` -> `alignColData` in `plotHeatmap` -
  `ignoreUnassignedVar` -> `filterIgnoresUnassigned` in `mergeClusters`
  (and other functions) for clarity. - `removeNegative` ->
  `removeUnassigned` in `getBestFeatures` for uniformity - Removed
  `largeDataset` option to `subsampleClustering` because no longer
  provides advantage. - `nBlank` -> `nBlankFeatures` in `makeBlankData`
  to allow for samples

- Created functions: - `primaryClusterLabel` and `primaryClusterType` -
  `getReducedData` - `assignUnassigned`: assigns unassigned samples to
  nearest cluster - `renameClusters` and `recolorClusters`: assign new
  names/colors to clusters within a particular clustering -
  `clusterMatrixColors`: wrapper to `convertClusterLegend` to return
  matrix like `clusterMatrix` only with colors in place of the internal
  cluster ids (like existing `clusterMatrixNamed`) -
  `plotClustersTable` for plotting a heatmap showing the results of
  `tableClusters` - `subsetByCluster` for subsetting CE object to only
  those samples in a particular cluster(s) of a clustering. -
  `plotFeatureScatter` for a scatter plot of 2+ features (genes)
  colored by cluster - `addToColData` and `colDataClusters` adding
  clustering information to colData of object. `addToColData` returns
  object with `colData` augmented, while `colDataClusters` just returns
  the `DataFrame` with clusterings added. - `updateObject` to update
  historical object created from previous versions to the current class
  definitions.

- Added arguments: - `whichAssay` to all functions to allow the user to
  select the assay on which the operations will be performed. -
  `stopOnErrors` to `RSEC` - `nColLegend` to `plotReducedDims` -
  `subsample` and `sequential` to `RSEC` to allow for opting out of
  those options (but default is `TRUE` unlike `clusterMany`) -
  `nBlankSamples` and `groupsOfSamples` to `makeBlankData` to allow for
  separating samples (columns) - `add` and `location` to
  `plotClusterLegend`

- `makeBlankData` will now allow for making blank columns to separate
  groups of samples.

- `plotDendrogram` now allows for plotting of `colData` (previously
  `sampleData`) like `plotHeatmap` or `plotClusters`

- `clusterMany` now allows user-defined `ClusterFunction` objects to
  argument `clusterFunction`.

- Removed restriction in `plotClustersWorkflow` that only
  `clusterType="clusterMany"` allowed.

- Allow `getClusterManyParams` to search old `clusterMany` runs as
  well.

- Added `table` method to `plotHeatmap` (for plotting heatmap of
  results of `table` function)

- Added error catch if try to give argument `whichCluster` to
  `mergeCluster`.

- Added error catch if give param `whichClusters` to functions that
  only take `whichCluster` (singular) as an argument

- `plotFeatureBoxplot` now returns (invisibly) the colors and
  clusterIds along with the boxplot information.

Bugs

- Add check for merge_nodeMerge table that `mergeClusterId` column
  can't be NA for entries where `isMerged=TRUE`

- Fix internal .makeIntegerClusters so that if given values `1:K` for
  input clustering will retain these same values (Issue #227)

- `mergeClusters` now returned object saves the merge information (and
  deletes old info and updates clusterType/clusterLabel of existing
  merge clusters), even if `mergeMethod="none"`.

- Fix `removeClusterings` so doesn't loose merge info unless deleting
  relevant clusterings.

Changes in version 2.1.4 (2018-06-27):

Bugs

- Fix tests that fail in devel version of Bioconductor (due to changes
  in other packages)

Changes in version 2.1.3 (2018-05-24):

Bugs

- Fix bug in how `clusterMany` and `defaultNDims` dealt with filtering
  choices in `reduceMethod`.

- Fixed bug in how `clusterMany` assigned label names

- Fixed bug so that `clusterMany` labels are increased (iteration
  version added) if `clusterMany` is rerun. Also, user-defined labels
  for functions like `mergeClusters` are now updated with iteration
  value if they are duplicated.

Changes in version 2.1.2 (2018-05-17):

Bugs

- Fix so that estimates of the proportion non-null in `mergeClusters`
  are always positive.

- Fix so that calculation of filter with `ignoreUnassignedVar` doesn't
  delete existing base filter of same type.

- Fix `plotClustersWorkflow` where wrong cluster was grabbed to plot.

- Fix bug in `plotDendrogram` where colors of clusters plotted with
  could be subsumbed by color of previous clustering. (git Issue
  `#220`)

Changes in version 2.1.1 (2018-05-15):

Bugs: 

- Fixed error introduced in 2.0.0 where arguments to `mergeCutoff` and
  `mergeLogFCcutoff` were passed instead to `dendroNDims`.

[cogena](https://bioconductor.org/packages/cogena)
------

Changes in version 1.15.1 (2018-05-30):

- add dotplot in heatmapPEI and vignettes

- multi-group labelling

[COMPASS](https://bioconductor.org/packages/COMPASS)
-------

Changes in version 1.19.4:

- Fixed a bug in the FunctionalityScore and PolyfunctionalityScore APIs
  when `markers=` argument is passed in and not null.

[compcodeR](https://bioconductor.org/packages/compcodeR)
---------

Changes in version 1.17.1:

- Removed support for SAMseq since samr has been removed from CRAN

[compEpiTools](https://bioconductor.org/packages/compEpiTools)
------------

Changes in version 1.15.1:

- The following function is updated + heatmapPlot.R: the function was
  updated removing a useless call to hclust, thus reducing by half the
  time required for clustering.  + topGOres.R: the function was sligtly
  modified to be more efficient in case of GO enrichment analysis of
  multiple gene sets. The indexing of the ontology is now performed
  only once, greatly improving the speed of the function.

[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap)
--------------

Changes in version 1.19.1:

- `Heatmap()`: no column name added if the input matrix is a one-column
  matrix.

- `oncoPrint()`: scales the the row annotations are now the same if
  rows are split.

[consensus](https://bioconductor.org/packages/consensus)
---------

Changes in version 0.99.0:

- Original submission to Bioconductor

[coRdon](https://bioconductor.org/packages/coRdon)
------

Changes in version 0.99.11:

- New package coRdon, for analysis of codon usage.

[countsimQC](https://bioconductor.org/packages/countsimQC)
----------

Changes in version 0.99.0:

- To comply with Bioconductor guidelines, the random seed (for
  subsampling of samples and features for dispersion and correlation
  calculations) is no longer set inside the functions, and the argument
  is removed from `countsimQCReport()`. For reproducible results,
  please set the random seed explicitly in the R session.

Changes in version 0.5.2:

- Add options to silence progress indicators

- Fixes in dispersion visualization

Changes in version 0.5.0:

- Add a number of quantitative evaluation criteria and tests for
  pairwise comparison of data sets

- The argument `subsampleSize` to the `countsimQCReport()` function now
  determines the number of observations for which (time-consuming)
  statistics are calculated

- A new argument `maxNForCorr` is added to the `countsimQCReport()`
  function to indicate the number of observations for which pairwise
  correlations are calculated

- Improvements in documentation

Changes in version 0.4.6:

- Allow code folding when showCode = TRUE

Changes in version 0.4.5:

- Allow data frames or matrices as input (assuming design = ~ 1)

- First implementation of area between ECDFs

- Increase transparency of points in scatter plots

Changes in version 0.4.4:

- Included effective library size evaluation.

Changes in version 0.4.2:

- Allowed using just a subset of the samples for dispersion
  calculations.

Changes in version 0.4.1:

- Added pairwise sample and variable correlation distributions.

- Added box plots and violin plots.

Changes in version 0.4.0:

- Major overhaul of function and argument names for increased
  consistency.

Changes in version 0.3.2:

- Added mean-variance plots.

Changes in version 0.3.0:

- Added K-S statistics.

- Added line density plots.

[CrispRVariants](https://bioconductor.org/packages/CrispRVariants)
--------------

Changes in version 1.9.2:

- Keep metadata columns in targets

- Adds function for splitting insertion sequences

[csaw](https://bioconductor.org/packages/csaw)
----

Changes in version 1.16.0:

- Added normFactors() function to avoid confusion when normOffsets()
  returns factors.

- Deprecated type="scaling" option in normOffsets().

- Added calculateCPM() function for convenient calculation of
  (log-)CPMs.

- Split up consolidateSizes() function into consolidateWindows(),
  consolidateTests() and consolidateOverlaps(). Deprecated
  consolidateSizes() itself.

- Switched output of combineTests() and getBestTest() and related
  functions to a DataFrame.

- Modified mergeWindows() behaviour with specified sign=, for dealing
  with nested windows of opposing sign.

- Altered controlClusterFDR() to take the largest adjusted p-value
  threshold that yields a cluster-level FDR below target=.

- Simplified detailRanges() output so that it no longer returns an
  arbitrary exon number.

[customProDB](https://bioconductor.org/packages/customProDB)
-----------

Changes in version 1.20.3:

UPDATED FUNCTIONS

- Fix bug in function PrepareAnnotationRefseq.R

- Update function aaVariation.R to avoid translate 'TTG' and 'CTG' into
  'M'

[cydar](https://bioconductor.org/packages/cydar)
-----

Changes in version 1.6.0:

- Restructured the CyData class for simplicity and internal fields.

- Deprecated plotCell* functions, renamed them to plotSphere*.

- Added the createColorBar() convenience function.

- Removed the diffIntDist() function.

- Restored option for quantile normalization in normalizeBatch().
  Switched to deterministic algorithm for sampling when mode="warp".

[cytofkit](https://bioconductor.org/packages/cytofkit)
--------

Changes in version 3.6:

- Contains all modifications and fixes from 1.8.2 to 1.9.5

[dada2](https://bioconductor.org/packages/dada2)
-----

Changes in version 1.9.2:

NEW FEATURES

- PacBio CCS reads up to 3 kilobases are now supported. See also
  PacBioErrfun, the new and recommended error-estimation function for
  PacBio CCS data.

SIGNIFICANT USER-VISIBLE CHANGES

- primer.fwd has been replaced by orient.fwd in the filterAndTrim
  function. This option consistently orients mixed-orientation
  single-end or paired-end reads based on matching the provided
  sequence fragment to the start or end of each read (or paired read).
  Intended for use with mixed-orientation reads that included sequenced
  primers. If primers aren't included in the amplicons, an external
  re-orientation solution remains preferable.

- trimRight has been added to the filterAndTrim function. This removes
  the specified number of bases from the end ("right" side) of each
  read.

BUG FIXES

- collapseNoMismatch now collapses identical sequences as well
  (previous behavior is togglable).

Changes in version 1.9.1:

BUG FIXES

- mergePairs now gracefully handles cases when zero reads succesfully
  merge.

- plotQualityProfile now works correclty when given a directory
  containing fastq files.

[DaMiRseq](https://bioconductor.org/packages/DaMiRseq)
--------

Changes in version 1.5.2:

- The DaMiR.normalization function embeds also the 'logcpm'
  normalization.

- Now, DaMiR.EnsembleLearning calculates also the Positive Predicted
  Values (PPV) and the Negative Predicted Values (NPV).

- Three new functions have been implemented for the binary
  classification task: DaMiR.EnsembleLearning2cl_Training,
  DaMiR.EnsembleLearning2cl_Test and DaMiR.EnsembleLearning2cl_Predict.
  The first one allows the user to implement the training task and to
  select the model with the highest accuracy or the average accuracy;
  the second function allows the user to test the selected
  classification model on a test set defined by the user; the last
  function allows the user to predict the class of new samples.

- Removed black dots in the violin plots.

Changes in version 1.4.1:

- Adjusted Sensitivity and Specificity calculations.

[DChIPRep](https://bioconductor.org/packages/DChIPRep)
--------

Changes in version 1.11.1:

- fixed unit test for plotting without replicates as DESeq2 no longer
  supports silentl treating as replicates since version 1.22

[debrowser](https://bioconductor.org/packages/debrowser)
---------

Changes in version 1.9.21:

- Last publication fixes

Changes in version 1.8.6:

- The columns that have strings removed while loading

- Labels can be changed in the main plots

- Plot Information box added to give more information about the plots.

Changes in version 1.8.4:

- Calculated padj and foldchange columns added to all detected genes
  result

- Normalization issue is fixed in comparison table

Changes in version 1.8.3:

- EdgeR and Limma

Changes in version 1.8.2:

- Package installation changed to warning

- Dependancies removed from DESCRIPTION to suppress loading messages

- Biarxiv citation added

Changes in version 1.8.1:

- Interactive heatmap height and width fix

- DEBrowser turned to a modular structure. The modules can be used in
  other shiny applications.

- More interactivity added to Heatmaps and main plots.

- Lasso selection is added to main plots

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

Changes in version 1.17.8:

- Fix: typo in as.DEGSet from DESeq2 object.

Changes in version 1.17.7:

- Fix: re-format NEWS file

Changes in version 1.17.6:

- Fix: Error in degCovariate when no correlation exists.

Changes in version 1.17.5:

- Fix: Fix documentation warning.

Changes in version 1.17.4:

- Feature: Reduce parameter is used to remove outlier points after
  clustering genes.

- Feature: Add maximum log2FoldChange to the significance output when
  multiple comparisons are used as inputs.

- Feature: Remove non-mapped genes in degPlot.

- Feature: Add specific function to plot degPatterns clusters.

- Fix: Support DESeqResults for list of DEGSets.

- Feature: Add variable selection for covariates that correlate with
  PCs in degCovariate function.

- Feature: Add lasso as an option to variable selection in covariate
  analysis.

- Feature: Fill with colors only significant variables by lm or lasso,
  and draw borber for correlated variables by cor.test.

Changes in version 1.17.3:

- Fix: degCovariates works with metadata only with numerical variables

- * Feature: Improve support for DEGSet conversion.

- Fix: Remove theme set up for degPCA plot.

- Feature: Make function to generate colors for metadata variables for
  annotation column in heatmap figure.

- Feature: Improve degCovariates to add effect size of the covariates.
  Thanks to @vbarrera

- Fix: Improve degCovariates man pages.

Changes in version 1.17.1:

- Fix: remove discrete scale color in degPCA.

- Feature: Return same output for degPatterns with single genes. Thanks
  Amir Jassim.

- Feature: Allow custom y-axis lab in degPlot. Thanks @vbarrera.

[DEGseq](https://bioconductor.org/packages/DEGseq)
------

Changes in version 1.34.1:

- remove function samWrapper

[DelayedArray](https://bioconductor.org/packages/DelayedArray)
------

Changes in version 0.8.0:

NEW FEATURES

- Add get/setAutoBlockSize(), getAutoBlockLength(),
      get/setAutoBlockShape() and get/setAutoGridMaker().

- Add rowGrid() and colGrid(), in addition to blockGrid().

- Add get/setAutoBPPARAM() to control the automatic 'BPPARAM' used by
      blockApply().

- Reduce memory usage when realizing a sparse DelayedArray to disk
    
      + On-disk realization of a DelayedArray object that is reported to be sparse
      (by is_sparse()) to a "sparsity-optimized" backend (i.e. to a backend with
      a memory efficient write_sparse_block() like the TENxMatrix backend imple-
      mented in the HDF5Array package) now preserves sparse representation of
      the data all the way. More precisely, each block of data is now kept in
      a sparse form during the 3 steps that it goes thru: read from seed,
      realize in memory, and write to disk.

- showtree() now displays whether a tree node or leaf is considered sparse
      or not.

- Enhance "aperm" method and dim() setter for DelayedArray objects. In
      addition to allowing dropping "ineffective dimensions" (i.e. dimensions
      equal to 1) from a DelayedArray object, aperm() and the dim() setter now
      allow adding "ineffective dimensions" to it.

- Enhance subassignment to a DelayedArray object.
    
      + So far subassignment to a DelayedArray object only supported the **linear
      form** (i.e. x[i] <- value) with strong restrictions (the subscript 'i'
      must be a logical DelayedArray of the same dimensions as 'x', and 'value'
      must be an ordinary vector of length 1).
    
      + In addition to this linear form, subassignment to a DelayedArray object
      now supports the **multi-dimensional form** (e.g. x[3:1, , 6] <- 0). In
      this form, one subscript per dimension is supplied, and each subscript
      can be missing or be anything that multi-dimensional subassignment to
      an ordinary array supports. The replacement value (a.k.a. the right
      value) can be an array-like object (e.g. ordinary array, dgCMatrix object,
      DelayedArray object, etc...) or an ordinary vector of length 1. Like the
      linear form, the multi-dimensional form is also implemented as a delayed
      operation.

- Re-implement internal helper simple_abind() in C and support long arrays.
      simple_abind() is the workhorse behind realization of arbind() and
      acbind() operations on DelayedArray objects.

- Add "table" and (restricted) "unique" methods for DelayedArray objects,
      both block-processed.

- range() (block-processed) now supports the 'finite' argument on a
      DelayedArray object.

- %*% (block-processed) now works between a DelayedMatrix object and an
      ordinary vector.

- Improve support for DelayedArray of type "list".

- Add TENxMatrix to list of supported realization backends.

- Add backend-agnostic RealizationSink() constructor.

- Add linearInd() utility for turning array indices into linear indices.
      Note that linearInd() performs the reverse transformation of
      base::arrayInd().

- Add low-level utilities mapToGrid() and mapToRef() for mapping reference
      array positions to grid positions and vice-versa.

- Add downsample() for reducing the "resolution" of an ArrayGrid object.

- Add maxlength() generic and methods for ArrayGrid objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Multi-dimensional subsetting is no more delayed when drop=TRUE and the
      result has only one dimension. In this case the result now is returned
      as an **ordinary** vector (atomic or list). This is the only case of
      multi-dimensional single bracket subsetting that is not delayed.

- Rename defaultGrid() -> blockGrid(). The 'max.block.length' argument
      is replaced with the 'block.length' argument. 2 new arguments are
      added: 'chunk.grid' and 'block.shape'.

- Major improvements to the block processing mechanism.
      All block-processed operations (except realization by block) now support
      blocks of **arbitrary** geometry instead of column-oriented blocks only.
      'blockGrid(x)', which is called by the block-processed operations to get
      the grid of blocks to use on 'x', has the following new features:
      + It's "chunk aware". This means that, when the chunk grid is known (i.e.
         when 'chunkGrid(x)' is not NULL), 'blockGrid(x)' defines blocks that
         are "compatible" with the chunks i.e. that any chunk is fully contained
         in a block. In other words, blocks are chosen so that chunks don't
         cross their boundaries.
      + When the chunk grid is unknown (i.e. when 'chunkGrid(x)' is NULL),
         blocks are "isotropic", that is, they're as close as possible to an
         hypercube instead of being "column-oriented" (column-oriented blocks,
         also known as "linear blocks", are elongated along the 1st dimension,
         then along the 2nd dimension, etc...)
      + The returned grid has the lowest "resolution" compatible with
         'getAutoBlockSize()', that is, the blocks are made as big as possible
         as long as their size in memory doesn't exceed 'getAutoBlockSize()'.
         Note that this is not a new feature. What is new though is that an
         exception now is made when the chunk grid is known and some chunks
         are >= 'getAutoBlockSize()', in which case 'blockGrid(x)' returns a
         grid that is the same as the chunk grid.
      + These new features are supposed to make the returned grid "optimal" for
      block processing. (Some benchmarks still need to be done to
      confirm/quantify this.)

- The automatic block size now is set to 100 Mb (instead of 4.5 Mb
      previously) at package startup. Use setAutoBlockSize() to change the
      automatic block size.

- No more 'BPREDO' argument to blockApply().

- Replace block_APPLY_and_COMBINE() with blockReduce().

BUG FIXES

- No-op operations on a DelayedArray derivative really act like no-ops.
      Operating on a DelayedArray derivative (e.g. RleArray, HDF5Array or
      GDSArray) will now return an objet of the original class if the result
      is "pristine" (i.e. if it doesn't carry delayed operations) instead of
      degrading the object to a DelayedArray instance. This applies for example
      to 't(t(x))' or 'dimnames(x) <- dimnames(x)' etc...


[derfinder](https://bioconductor.org/packages/derfinder)
---------

Changes in version 1.15.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Switched from 'outfile' to 'log' when invoking
  BiocParalle::SnowParam(). Thus define_cluster() now has a mc.log
  argument instead of mc.outfile.

Changes in version 1.15.2:

BUG FIXES

- Fix a message regarding the deprecated IRanges subset method.

Changes in version 1.15.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Use BiocManager

[derfinderHelper](https://bioconductor.org/packages/derfinderHelper)
---------------

Changes in version 1.15.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Use BiocManager

[derfinderPlot](https://bioconductor.org/packages/derfinderPlot)
-------------

Changes in version 1.15.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Use BiocManager

[DESeq2](https://bioconductor.org/packages/DESeq2)
------

Changes in version 1.22.0:

- No replicate designs no longer supported (previous version began
  deprecation with a warning).

- unmix() now optionally will return the correlation (in the variance
  stabilized space) of the fitted data to the original data, and the
  matrix of fitted data (format="list"). Argument 'loss' was changed to
  'power'. Will give warning if the columns of 'pure' have high
  correlation (in the variance stabilized space).

Changes in version 1.21.21:

- Improved code for 'linearModelMu' (an internal fitting function used
  in dispersion estimation for some models) contributed by Wolfgang
  Huber speeds up an internal step by 2 orders of magnitude.

Changes in version 1.21.15:

- Rows of the weights matrix which would produce a degenerate design
  matrix, instead of giving an error, will produce a warning, and these
  rows will be treated as if they contained all zeros
  (mcols(dds)$allZero and mcols(dds)$weightsFail will be set to TRUE).

Changes in version 1.21.14:

- The nbinom{WaldTest,LRT} functions will not stop if the design
  produces a model matrix that is not full rank and betaPrior=FALSE
  (default). This was assumed by the DESeq2 code, because errors are
  produced at object construction and at dispersion estimation, but it
  was possible to call nbinomLRT() from DEXSeq after dispersion
  estimation and end up with a full model matrix that was not full
  rank. Instead testForDEU() should be called from DEXSeq.

Changes in version 1.21.13:

- Adding back a feature from version 1.15, where contrasts of two
  groups where both had all zero counts would have the LFC zero-ed out,
  rather than output a small but non-zero value. It's preferable for
  the Wald test that the LFC be set to zero for such contrasts.

Changes in version 1.21.9:

- DESeq() now only says one time 'using supplied model matrix',
  previously this was repeated three times from sub-functions.
  Sub-functions therefore no longer print this message.

- Fixed bug when lfcShrink run directly after LRT with supplied model
  matrices.

- Added heuristic to prevent Cook's outlier based filtering when the
  max Cook's sample has lower counts than 3 other samples. Restricted
  to two group comparison datasets.

[DEsingle](https://bioconductor.org/packages/DEsingle)
--------

Changes in version 1.0.5:

- Optimization of speed and memory.

Changes in version 1.0.1:

- Optimization of memory management.

[DEsubs](https://bioconductor.org/packages/DEsubs)
------

Changes in version 1.7.4:

- Reformatted NEWS file.

Changes in version 1.7.3:

- Updated rankedList argument in DEsubs() to receive gene identifiers
  in any of the supported mRNAnomenclatures.

Changes in version 1.7.2:

- Removed 'Significance Analysis of Microarrays' (SAM) from available
  differential expression analysis options due to the removal of its
  package from the ecosystem.

[diffHic](https://bioconductor.org/packages/diffHic)
-------

Changes in version 1.14.0:

- Added the readMTX2IntSet() function to create InteractionSets from
  file.

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

Changes in version 1.7.3:

- Added citation

[DMCHMM](https://bioconductor.org/packages/DMCHMM)
------

Changes in version 1.3.1:

CHANGES SINCE LAST TIME

- A new citation is added.

- Several bugs are corrected.

[dmrseq](https://bioconductor.org/packages/dmrseq)
------

Changes in version 1.1.20 (2018-10-11):

- Minor bug fix for a rare situation that occured in the case of a
  multi-level covariate of interest if the region level model fitting
  procedure did not converge.

Changes in version 1.1.2 (2018-05-09):

- The newly added `chrsPerChunk` argument specifies the number of
  chromosomes to compute at a time (default is 1).

[DominoEffect](https://bioconductor.org/packages/DominoEffect)
------------

Changes in version 1.1.2:

BUG FIX

- add missing space to vignette

Changes in version 1.1.1:

UDPATE

- move from bioclite to BiocManager

Changes in version 1.1.0:

NEW FEATURES

- DominoEffect in BioC 3.8 development release

[drawProteins](https://bioconductor.org/packages/drawProteins)
------------

Changes in version 1.1.1:

FEATURES

- New function called draw_folding(). This function allows the drawing
  of the STRAND, HELIX and TURN types which denote regions of the
  proteins that assemble as beta-strands, alpha helicies or make a turn
  in the 3D structure of the protein.

- New function called parse_gff(). This function imports files or urls
  that link to a GFF3 format and parses the data to allow it to be
  plotted.

[DropletUtils](https://bioconductor.org/packages/DropletUtils)
------------

Changes in version 1.2.0:

- Added removeSwappedDrops() for removing swapping in other types of
  droplet-based data.

- Added alpha= argument to testEmptyDrops() to support overdispersion
  during sampling. Returned arguments and estimates in metadata of
  testEmptyDrops(), emptyDrops().

- Added encodeSequences() for convenient 2-bit encoding of sequences.

- Added get10xMolInfoStats() function to compute per-cell statistics
  from a molecule info file.

- Deprecated read10xMatrix(), as it does not add much practical value
  over Matrix::readMM().

- Support the 10X sparse HDF5 format in read10xCounts().

- Support the 10X sparse HDF5 format in write10xCounts().

[dupRadar](https://bioconductor.org/packages/dupRadar)
--------

Changes in version 1.11.1:

- Added support to count on different feature types rather than exons
  (thanks to John Chuang)

[edgeR](https://bioconductor.org/packages/edgeR)
-----

Changes in version 3.24.0:

- New functions catchKallisto() and catchSalmon() to read outputs from
  kallisto and Salmon and to compute overdispersion factors for each
  transcript from bootstrap samples.

- New function readBismark2DGE() to read coverage files created by
  Bismark for BS-seq methylation data.

- New method "TMMwzp" for calcNormFactors() to better handle samples
  with a large proportion of zero counts.

- The default value for prior.count increased from 0.25 to 2 in cpm()
  and rpkm(). The new value is more generally useful and agrees with
  the default values in aveLogCPM() and with the DGEList method for
  plotMDS().

- zscoreNBinom() now supports non-integer q values.

- The scaleOffset() S3 methods for DGEList and default objects are now
  registered in the NAMESPACE. Previously the functions were exported
  but not registered as S3 methods.

- The rowsum() method for DGEList objects (rowsum.DGEList) now
  automatically removes gene annotation columns that are not
  group-level.

- More specific error messages from DGEList() when invalid (NA,
  negative or infinite) count values are detected.

- Bug fix to glmfit.default() when lib.size is specified.

- Bug fix to column name returned by decideTestsDGE().

[EGSEA](https://bioconductor.org/packages/EGSEA)
-----

Changes in version 1.9.1:

- fixed: a minor bug to allow users generate EGSEA reports for only a
  single base method

- fixed: several minor bugs

- Change: from biocLite to BiocManager

[EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano)
---------------

Changes in version 1.0.0:

- user can now supply her/his/its own colour vector to label the points

- default is to now draw grid lines and only have left and bottom
  borders

- when selectLab is not NULL, even variables that do not pass the
  thresholds are now labelled along with those that do, even when
  DrawConnectors is either TRUE or FALSE.

- correctly catches non-numeric x and / or y variables, and throws
  error.

- the function now tolerates P values of 0 (zero) and replaces these
  with the lowest possible double value, given a user's specific
  computer architecture and R version.

[EnrichedHeatmap](https://bioconductor.org/packages/EnrichedHeatmap)
---------------

Changes in version 1.11.1:

- add `as.normalizedMatrix()` function to convert a matrix to
  `normalizedMatrix` class

[EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser)
-----------------

Changes in version 2.12.0:

- Major refactoring of ID mapping for the rownames of a
  SummarizedExperiment: (functions idmap / probe2gene): - to.ID can
  also be a rowData column to support user-defined mappings - support
  of data-driven strategies for many:1 and 1:many mappings -
  synchronized behavior of microarray probe ID mapping (probe2gene) and
  general gene ID mapping (idmap)

- Alternative representation of gene sets based on GSEABase::GeneSet
  and GSEABase::GeneSetCollection to facilitate gene ID mapping for
  gene sets (function getGenesets)

- Output destination of HTML reports (functions eaBrowse / ebrowser):
  extended control via arguments out.dir and report.name that overwrite
  corresponding config defaults (configEBrowser)

- Separation of nominal and adjusted p-values in DE and EA result
  tables (functions deAna / sbea / nbea)

[ensembldb](https://bioconductor.org/packages/ensembldb)
---------

Changes in version 2.5.9:

- Write the official scientific name into the "Organism" metadata
  field.

Changes in version 2.5.8:

- Further improve MySQL support and performance.

Changes in version 2.5.6:

- Add additional (integer) ID columns to the tables for the MySQL
  backend to improve performance.

- Use integer primary key columns for join queries in MySQL/MariaDB
  EnsDb databases.

Changes in version 2.5.5:

- Switch from RMySQL to RMariaDB.

Changes in version 2.5.2:

- Switch from GenenameFilter to GeneNameFilter for AnnotationHub >=
  1.5.2.

Changes in version 2.5.1:

- Fix bug in getGeneRegionTrackForGviz that throws an error if both
  protein coding and non-coding transcripts are fetched.

[ensemblVEP](https://bioconductor.org/packages/ensemblVEP)
----------

Changes in version 1.24.0:

NEW FEATURES

- add support for Ensembl release 94 (this should also encompass 93)


[ERSSA](https://bioconductor.org/packages/ERSSA)
-----

Changes in version 0.99.8:

- Fix a bug in calling BiocParallel

Changes in version 0.99.7:

- Code improvements

Changes in version 0.99.5:

- Increase font sizes in plots

- Add line of mean DEG discovery in interect plot

- Improve clarity in plots

- Vignette receives major edits to improve clarity

Changes in version 0.99.4:

- Update installation instruction in vignette

Changes in version 0.99.3:

- Add additional citations

- Reduce example runtime with fewer combinations

Changes in version 0.99.2:

- Remove set.seed from package code

Changes in version 0.99.0:

- Codes ready for bioconductor submission

- Added full example dataset

- Cleanup the outputs to be more descriptive and append ERSSA_ to all
  files

- Write vignette

[esetVis](https://bioconductor.org/packages/esetVis)
-------

Changes in version 1.7.3:

- vignette: use print rather 'knit_print.ggvis'

- ggvis: include transparency, fix issue position legend

- rbokeh: use vector instead of column names for ly_hexbin

Changes in version 1.7.2:

- ggplotEsetPlot: enable title/axes labels of type expression

Changes in version 1.7.1:

- fix issue duplicated columns when same column is used for multiple
  aesthetics (reported as issue in ggplot2: 2.2.1.9000)

[EventPointer](https://bioconductor.org/packages/EventPointer)
------------

Changes in version 2.0:

- Minor bugs fixed (ClassifyEvents)

- Relative error is now obtained when obtaining PSI

- New statistical analysis based on PSI and the associated relative
  error

- New pipeline for RNA-Seq based on quantification using a reference
  transcriptome

- Multi-path events detection o Events with more than two alternative
  paths are detected

[ExCluster](https://bioconductor.org/packages/ExCluster)
---------

Changes in version 0.99.13:

- Fixed bugs that would cause ExCluster to either crash or not properly
  plot results in some cases

Changes in version 0.99.12:

- Loosened FDR calculations slightly (ExCluster was a bit too
  stringent)

- Added plot.Type option to plotExonlog2FC function, which accepts
  "bitmap" and "PNG"

- Bug-tested plot.Type so machines with at least Ghostscript or X11
  forwarding will have minimal issues

- Changed how files/folders are written & how write-permissions are
  detected, to avoid bugs

- Updated the vignette to use BiocManager instead of biocLite

Changes in version 0.99.11:

- Removed some duplicated code (stripping ID numbers, computing
  p-values)

- removed several instances of cat() and print() and replaced them with
  message()

- changed apply() code to use matrixStats instead, which is up to 500
  times faster (Thanks Lori!)

- this previous change sped the algorithm up from 1 hour+ to only ~ 20
  minute runtime

- changed the test dataset so one of the genes has a p-value < 0.05 and
  plots results

Changes in version 0.99.10:

- The GFF_convert function now outputs GFF3 formatted annotations

- Amended other functions (processCounts, ExCluster) to work on GFF3
  formatted annotations

- Changed the output of GFF_convert to be a GRanges object of said GFF3
  annotations

- Created a separate, internal, load_ExCluster_functions.R script to
  load helper functions

- Removed many depenencies on variables outside the environment for
  said functions

- Created a library of error messages in the ExCluster_errors.R script
  (internal)

- Separated the large ExCluster function into ExCluster.R,
  ExClust_compute_stats.R, and ExClust_main_function.R

Changes in version 0.99.9:

- Added the function rtracklayerGTFtoGFF, which flattens GTF files
  imported by rtracklayer to GFF format

- Added the function GRangesFromGFF, which converts GFF formatted data
  to GRanges format

- Added the function GRangesFromExClustResults, which converts
  ExCluster function results to GRanges format

- Removed repeated code from the GFF_convert function

- Most functions now have checks to ensure GTF, GFF, and ExClustResults
  data is formatted correctly'

- Made an improvement to the processCounts function, which handles some
  edge cases better (zero reads in some conditions)

Changes in version 0.99.8:

- Fixed build check errors from 0.99.7

Changes in version 0.99.7:

- Fixed a number of bugs resulting from ExCluster() function code
  changes

- Updated the Vignette to correctly reflect changes suggested from the
  Bioconductor review

- Added several package imports to better keep track of global
  variables/functions

- Reduced the lengths of code lines in a number of instances

- Cleaned up error messages for ExCluster() function

Changes in version 0.99.6:

- Completely re-wrote the GFF_convert() function to have helper
  functions and take advantage of GenomicRanges

- Changed processCounts, which was incorrectly counting reads as
  stranded (is now set to unstranded)

- Altered ExCluster null hypothesis simulations to run faster (about 4
  times faster now)

- Made numerous small changes to address the issues with first
  Bioconductor review

Changes in version 0.99.5:

- Fixed incorrectly pushed changes, which caused build error

Changes in version 0.99.4:

- Reformatted the Vignette & DESCRIPTION to process the Vignette
  correctly with Sweave

Changes in version 0.99.3:

- Removed extraneous files causing errors in the build

Changes in version 0.99.2:

- Reformatted the Vignette & DESCRIPTION to process the Vignette
  correctly with Sweave

Changes in version 0.99.1:

- Removed extraneous files causing errors in the build

Changes in version 0.99.0:

- ExCluster package now passes R CMD build and check, and BiocCheck!

- takes about 2-3 hours to run on large datasets

- added a "Toy Dataset" for function manual runnable examples within
  the package

- more code in ExCluster has been turned into functions to reduce
  repeated code

- still some repeated code (such as parsing EnsIDs and exon bins)

[ExperimentHubData](https://bioconductor.org/packages/ExperimentHubData)
-----------------

Changes in version 1.7.0:

NEW FEATURES

- There are now validity checks for tags (valid biocViews in
  DESCRIPTION), RDataPath (must be defined), SourceUrls (must be
  defined and indicate url), genome (valid length), sourcetype (from
  list of valid entries), species and taxonomy id (validated from
  GenomeInfoDbData)

MODIFICATIONS

- Modified to use BiocManager instead of BiocInstaller

[FELLA](https://bioconductor.org/packages/FELLA)
-----

Changes in version 1.1.6:

- Fixed vignette indices

Changes in version 1.1.5:

- Added full vignette on the zebrafish dataset

- Small modifications to the Mus musculus vignette

Changes in version 1.1.4:

- Added a full vignette showing a case study on Mus musculus

- Functions `getCom` and `getGraph` are exported now

- Added `DT` and other packages to `suggests`

- Small fixes for `BiocCheck`

Changes in version 1.1.3:

- Version 1.1.2 did not skip such test

Changes in version 1.1.2:

- Disabled `buildGraphFromKEGGREST` test in 32-bit Windows due to its
  memory usage

[FGNet](https://bioconductor.org/packages/FGNet)
-----

Changes in version 3.15:

- plotKegg() has been deprecated. It will be removed from the package
  in upcoming versions.

[fgsea](https://bioconductor.org/packages/fgsea)
-----

Changes in version 1.7.1:

- Setting colwidth to zero make column not to be drawn

- Changable line width in plotEnrichment

[flowAI](https://bioconductor.org/packages/flowAI)
------

Changes in version 1.10.1:

- add parallel processing compatibility

[flowCL](https://bioconductor.org/packages/flowCL)
------

Changes in version 1.19.1:

- Changed how NEWS is displayed.

Changes in version 1.19.0:

- Same as previous. Bumped by BioConductor.

[FoldGO](https://bioconductor.org/packages/FoldGO)
------

Changes in version 0.99.10 (2018-09-27):

- Support for MgsaSets object (mgsa package) added

Changes in version 0.99.9:

- Documentation update

Changes in version 0.99.8:

- Man pages list shortened to 10 pages

- Documentation update

Changes in version 0.99.7 (2018-07-19):

- Argument specifying column with Gene IDs added to GAFReader class

- Example gaf file truncated (don't use it for real analysis!)

Changes in version 0.99.6 (2018-07-19):

- Malformed description filed problem fixed

- Imports for stats and methods packages added

Changes in version 0.99.5 (2018-07-16):

- R version dependency changed from version 3.4.4 to 3.5

Changes in version 0.99.4 (2018-07-16):

- Vignettes update: less chunks with eval=FALSE

- /data added to .Rbuildignore

Changes in version 0.99.3 (2018-07-13):

- Vignettes update

- Bugs with documentation fixed

- R CMD check with no warnings

Changes in version 0.99.2 (2018-07-12):

- Subscribed to the bioc-devel mailing list.

Changes in version 0.99.1 (2018-07-12):

- .Rproj file removed from repository

Changes in version 0.99.0 (2018-07-12):

- Submitted to Bioconductor

[gdsfmt](https://bioconductor.org/packages/gdsfmt)
------

Changes in version 1.17.1-1.17.6:

NEW FEATURES

- new options 'recursive' and 'include.dirs' in `ls.gdsn()`: the
  listing recurses into child nodes

UTILITIES

- replace BiocInstaller biocLite mentions with BiocManager

- `digest.gdsn()` fails if the digest package is not installed

- SIMD optimization in 2-bit array decoding with a logical vector of
  selection (3x speedup when there are lots of zeros)

BUG FIXES

- bug fixed: `put.attr.gdsn()` fails to update the existing attribute

[GeneAccord](https://bioconductor.org/packages/GeneAccord)
----------

Changes in version 0.99.13 (2018-08-02):

- Improved the code to include more from tidyverse; simplified code

- Using now roxygen2 to create Namespace and importing each function
  indiviually

- Created a NEWS.Rd file parsable by utils::news in folder inst

[genefilter](https://bioconductor.org/packages/genefilter)
----------

Changes in version 1.64.0:

NEW FEATURES

- Add `na.rm =` to row/colttests, requested by
  https://github.com/Bioconductor/genefilter/issues/1

[GeneNetworkBuilder](https://bioconductor.org/packages/GeneNetworkBuilder)
------------------

Changes in version 1.23.2:

- Fix the problem for when the interaction network is breaking.

[GENESIS](https://bioconductor.org/packages/GENESIS)
-------

Changes in version 2.12.0:

- pcair and pcrelate have been completely rewritten for better
  consistency with other methods. Some argument names have changed; see
  the documentation. The output of pcrelate is now a list of
  data.frames instead of a list of matrices.

- pcrelateReadKinship and pcrelateReadInbreed are deprecated, as these
  tables are now returned by pcrelate.

- pcrelateMakeGRM is deprecated; use pcrelateToMatrix with new pcrelate
  output format.

- king2mat is deprecated; use kingToMatrix instead.

- fitNullMM, fitNullReg, assocTestMM, and admixMapMM are deprecated.
  assocTestSeq and assocTestSeqWindow are defunct. Use fitNullModel,
  assocTestSingle, assocTestAggregate, and admixMap instead.

Changes in version 2.11.15:

- Refactor pcrelate.

Changes in version 2.11.14:

- Added assocTestAggregate method for GenotypeData objects.

Changes in version 2.11.11:

- Refactor pcair.

Changes in version 2.11.8:

- Added admixMap function to replace admixMapMM.

Changes in version 2.11.4:

- Added assocTestSingle and fitNullModel methods for GenotypeData
  objects.

[GenomeInfoDb](https://bioconductor.org/packages/GenomeInfoDb)
------------------

Changes in version 1.18.0

NEW FEATURES

- Add checkCompatibleSeqinfo().

SIGNIFICANT USER-VISIBLE CHANGES

- Update genomeMappingTbl.csv, the db used internally by genomeBuilds()
      and family.


[GenomicDataCommons](https://bioconductor.org/packages/GenomicDataCommons)
------------------

Changes in version 1.5.8:

- filters are now chainable with pipes

[GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures)
------------

Changes in version 1.34.0:

NEW FEATURES

- 2 changes to makeTxDbFromGFF() / makeTxDbFromGRanges():
      + Now they support GFF3 files where the CDS parent is an exon instead
        of a transcript. Note that such GFF3 files are rare and not following
        the well established convention documented in the GFF3 specs:
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
      + Now they accept missing/invalid CDS phases (with a warning).

- makeTxDb() now accepts missing CDS phases.


[GenomicFiles](https://bioconductor.org/packages/GenomicFiles)
------------

Changes in version 1.18.0:

NEW FEATURES

- (v. 1.17.3) Add vcfFields,VcfStack-method.

[GenomicRanges](https://bioconductor.org/packages/GenomicRanges)
-------------

Changes in version 1.34.0:

NEW FEATURES

- Add coercions from GenomicRanges to IRangesList and from GenomicRanges
      to CompressedIRangesList. These 2 new coercions are equivalent to
      coercion from GenomicRanges to IntegerRangesList, that is, if 'gr' is
      a GenomicRanges object, the 3 following coercions are equivalent and
      return the same CompressedIRangesList object:
          as(gr, "IntegerRangesList")
          as(gr, "IRangesList")
          as(gr, "CompressedIRangesList")

DEPRECATED AND DEFUNCT

- Deprecate several RangedData methods: seqinfo, seqinfo<-, seqnames, and
      findOverlaps#RangedData#GenomicRanges

      + RangedData objects will be deprecated in BioC 3.9 (their use has been
      discouraged since BioC 2.12, that is, since 2014). Package developers
      that are still using RangedData objects need to migrate their code to
      use GRanges or GRangesList objects instead.

BUG FIXES

- Make [[, as.list(), lapply(), and unlist() fail more graciously on
      a GenomicRanges object.

- Make "show" methods for GenomicRanges and GPos objects robust to
      special metadata column names like "stringsAsFactors".

- Export the "update" method for GRanges objects. This addresses
      https://github.com/Bioconductor/GenomicRanges/issues/7

[GenomicScores](https://bioconductor.org/packages/GenomicScores)
-------------

Changes in version 1.6.0:

USER VISIBLE CHANGES

- Functions and classes deprecated in the previous release (scores,
  MafDb class) have been now removed from the package.

- Added support to latest release 2.1 of gnomAD MAF data, stored in
  packages MafDb.gnomAD.r2.1.hs37d5 and MafDb.gnomADex.r2.1.hs37d5.

[genphen](https://bioconductor.org/packages/genphen)
-------

Changes in version 1.9 (2018-10-20):

Changes

- Bayesian hierarchical model implemented => multiple comparison solved

- Transition from effect size metrics to GLMs

- Posterior predictive checks automatic

- Retrospective power analysis available

- Stan model debugging possible

- Multicore speedup

- Phylogenetic bias quantification

- Data reduction with diagnostics procedure based on random forest

To do

- Include a CONFIG file to change internal parameters

- Add practical examples where genphen has been used.

- Implement modules for data augmentation

- Update todo after data augmentation

[GGBase](https://bioconductor.org/packages/GGBase)
------

Changes in version 3.43.1:

- man page links modernized

- vignette to Rmd

- Warning added to vignette indicating that gQTLBase/gQTLstats are more
  modern


[GlobalAncova](https://bioconductor.org/packages/GlobalAncova)
------------

Changes in version 4.0.0:

- new function 'gGlobalAncova' for generalized linear models: groups of
  tested variables can be quantitative, categorical, ordinal and even
  of mixed types

- new function 'Plot.features' for showing contributions of individual
  variables to global test statistic; equivalent to 'Plot.genes', but
  for 'gGlobalAncova'

- new function 'gGlobalAncova.hierarchical' for hierarchical testing.

- new S4-class 'GAhier' for storing results of
  gGlobalAncova.hierarchical

- new methods 'show', 'results', 'sigEndnodes', and 'Plot.hierarchy' to
  access and visualize results from 'GAhier' objects

- added CITATION and NEWS files

- Vignettes are now created with knitr

[GOfuncR](https://bioconductor.org/packages/GOfuncR)
-------

Changes in version 1.1.5:

USER-LEVEL CHANGES

- add name and domain to categories in get_anno_categories()

- update GO-graph (version 11-Oct-2018)

Changes in version 1.1.2:

NEW FEATURES

- allow for custom gene coordinates (alternative to default
  Bioconductor annotation packages; coordinates are used when the
  familiy-wise error rate is corrected for gene length or spatial
  clustering of genes)

USER-LEVEL CHANGES

- gene coordinates used in correction for gene-length or spatial
  clustering are not restricted to chromosomes 1-22,X,Y,MT anymore

[gpart](https://bioconductor.org/packages/gpart)
-----

Changes in version 1.0.0:

NEW FEATURES

- Add new functions, BigLD, CLQD, GPART, LDblockHeatmap,
  convert2SumExpObj

[graphite](https://bioconductor.org/packages/graphite)
--------

Changes in version 1.27.6 (2018-10-25):

- Updated all pathway data.

- Added PathBank database.

Changes in version 1.27.2 (2018-05-22):

- Vignette describing mixed analysis with genes and metabolites.

[gtrellis](https://bioconductor.org/packages/gtrellis)
--------

Changes in version 1.13.2:

- add_rect_track(): height is not converted to mm anymore

[GWASTools](https://bioconductor.org/packages/GWASTools)
---------

Changes in version 1.27.1:

- Add GenotypeIterator and GenotypeBlockIterator classes. These classes
  allow returning blocks of SNPs with each call to iterateFilter.

[HDF5Array](https://bioconductor.org/packages/HDF5Array)
-----

Changes in version 1.10.0:

NEW FEATURES

- Implement the TENxMatrix container (DelayedArray backend for the
      HDF5-based sparse matrix representation used by 10x Genomics).
      Also add writeTENxMatrix() and coercion to TENxMatrix.

SIGNIFICANT USER-VISIBLE CHANGES

- By default automatic HDF5 datasets (e.g. the dataset that gets written
      to disk when calling 'as(x, "HDF5Array")') now are created with chunks
      of 1 million array elements (revious default was 1/75 of
      'getAutoBlockLength(x)'). This can be controlled with new low-level
      utilities get/setHDF5DumpChunkLength().

- By default automatic HDF5 datasets now are created with chunks of
      shape "scale" instead of "first-dim-grows-first". This can be
      controlled with new low-level utilities get/setHDF5DumpChunkShape().

- getHDF5DumpChunkDim() looses the 'type' and 'ratio' arguments (only 'dim'
      is left).


[HIBAG](https://bioconductor.org/packages/HIBAG)
-----

Changes in version 1.17.1:

- new function `hlaDistance()`

[HilbertCurve](https://bioconductor.org/packages/HilbertCurve)
------------

Changes in version 1.11.1:

- move IRanges and GenomicRanges to imports field in DESCRIPTION

[hipathia](https://bioconductor.org/packages/hipathia)
--------

Changes in version 1.1.2 (2018-09-07):

- Adding gene name to tooltip

- Fixing minor bugs in chart.R, save.R, stats.R

Changes in version 1.0.1 (2018-06-14):

- Fixing minor bug in heatmap_plot about variances in case
  variable_cluster = TRUE.

- Fixing minor bug in visualize_report.

- Adding test_package.R file.

[HPAanalyze](https://bioconductor.org/packages/HPAanalyze)
----------

Changes in version 0.99.19:

- Minor change: - Update CITATION

Changes in version 0.99.18:

- Minor changes: - Add 'journal' to CITATION - Add '+ file LICENSE' to
  DESCRIPTION

Changes in version 0.99.17:

- Minor changes: - Standardize the NEWS file - Add CITATION file

Changes in version 0.99.16:

- Minor changes: - Update License. Now use GPL-3 instead of MIT.

Changes in version 0.99.14:

- Minor changes: - Replace visType == 'all' in hpaVis() with
  identical(visType, 'all') - Replace extractType == 'all' in hpaXml()
  with identical(extractType, 'all')

Changes in version 0.99.13:

- Respond to Bioconductor review on Aug 6, 2018

- Minor changes: - Fix example for hpaXml() (issue #5) - Update
  documentation for hpa_downloaded_histology_v18 dataset (issue #6) -
  Re-name R/data.R to R/hpa_downloaded_histology_v18.R (issue #7) -
  Cross-reference documentations (issues #8, #9) - Combining
  hpaListParam() and hpaSubset() on one man page (issue #10)

Changes in version 0.99.12:

- Minor changes: - Minor edit in hpaXml() function

Changes in version 0.99.11:

- Respond to Bioconductor review on Jun 11, 2018

- Major changes: - Add the hpaVis() function as an umbrella for the
  whole function family - Add the hpaXml() function as an umbrella for
  all XML extraction

- Minor changes: - Remove the ignore/ directory - Import individual
  functions for xml2 - Remove /dontrun{} on man pages - Modifications
  to hpaVis functions for better consistency - Updated the help files
  with details about output of each function

Changes in version 0.99.10:

- Initital version reviewed by Bioconductor

[hpar](https://bioconductor.org/packages/hpar)
----

Changes in version 1.23.3:

- Use BiocManager for installation <2018-07-19 Thu>

Changes in version 1.23.2:

- Rename 16.1 objects <2018-05-23 Wed>

Changes in version 1.23.1:

- New Bioconductor devel

- Update to HPA version 18

[iasva](https://bioconductor.org/packages/iasva)
-----

Changes in version 0.99.3:

MODIFICATIONS

- updated R functions to use stopifnot() for error control

- changed cat() to message() for communicating output

- change use of R CRAN parallel functions to BiocParallel

- replace date field with hard-corded date in vignette

- update variable names to lower snake case in vignette

- added new biocViews (RNASeq, Software, StatisticalMethod,
  FeatureExtraction)

[illuminaio](https://bioconductor.org/packages/illuminaio)
----------

Changes in version 0.23.2 (2018-07-18):

NEW FEATURES

- The readIDAT() function gains a 'what' argument (only for unencrypted
  IDAT files).  This allows the fast return of the number of
  'nSNPsRead' (number of probes in file) and 'IlluminaID' (probenames).
  This is to allow fast handling of IDAT files in the minfi package.

Changes in version 0.23.1 (2018-07-18):

SOFTWARE QUALITY

- ROBUSTNESS: Now registering native routines.

Changes in version 0.23.0 (2018-05-01):

NOTES

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.8 for R (>= 3.5.0).

[immunoClust](https://bioconductor.org/packages/immunoClust)
-----------

Changes in version 1.13.3:

NEW FEATRUES 

- additional variant of t-mixture model fitting respecting truncated observed values

CHANGES 

- minor bugfixes and code cleaning

[InteractionSet](https://bioconductor.org/packages/InteractionSet)
--------------

Changes in version 1.10.0:

- Bug fix to seqinfo<- to support other arguments in the generic.

- Modified behaviour of inflate() with unspecified fill= for
  GInteractions objects.

[IntEREst](https://bioconductor.org/packages/IntEREst)
--------

Changes in version 1.5.2:

NEW FEATURES

- buildSsTypePwms() function supports the possibility to select begin
  and end point of splice sites sequences from which PWMs are built. It
  also supports the possibility to paste splice sites to build PWM.

BUG FIXES

- interest.sequential() and interest() corrections to their object
  output option.

- annotateU12() modified to work correctly with the new changes in
  Biostrings package.

- buildSsTypePwms() corrected and modified to better suit data for all
  species from SpliceRack and U12DB.

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

[IRanges](https://bioconductor.org/packages/IRanges)
----

Changes in version  2.16.0

SIGNIFICANT USER-VISIBLE CHANGES

- Optimize unlist() on Views objects.

- Optimize range(), any() and all() on CompressedRleList objects.

- Optimize start(), end(), width() setters on CompressedRangesList objects.

DEPRECATED AND DEFUNCT

- Deprecate several RangedData methods:
      + score, score<-, lapply, within, countOverlaps;
      + coercions from list, data.frame, DataTable, Rle, RleList, RleViewsList,
        IntegerRanges, or IntegerRangesList to RangedData.
      + RangedData objects will be deprecated in BioC 3.9 (their use has been
      discouraged since BioC 2.12, that is, since 2014). Package developers
      that are still using RangedData objects need to migrate their code to
      use GRanges or GRangesList objects instead.

     + The RangesList() constructor is now defunct (after being deprecated in
      BioC 3.7).

BUG FIXES

- Fix DF&#91;IRanges(...), &#93; on a DataFrame with data.frame columns.

- Make &#91;&#91;, as.list(), lapply(), and unlist() fail more graciously on
      a IRanges object.

- NCList objects now properly support c().


[iSEE](https://bioconductor.org/packages/iSEE)
----

Changes in version 1.1.13:

- Add missing observer for assay type in row data plot panels.

- Add missing observer for colorpicker when colouring by feature name
  in row-based plots, or by sample name in column-based plots.

- Ignore NA values when computing the range of coloring scales.

- Add a size expansion factor (5x) to the selected point when colouring
  by feature name in row-based plots, or by sample name in column-based
  plots.

- Fix redundant coloring of selected point when colouring by feature
  name in row-based plots, or by sample name in column-based plots.

- Update basic vignette.

Changes in version 1.1.12:

- Updated NEWS file.

Changes in version 1.1.11:

- Export list of panel names and codes.

Changes in version 1.1.10:

- Fix colour scale to be invariant when selecting on a different color.

- Protect heat map plot panels against restriction on zero samples.

Changes in version 1.1.9:

- Fix compatibility with DelayedArray assays.

Changes in version 1.1.8:

- Extend unit test coverage.

- Move generics to separate file.

- Minor fix to annotateEnsembl.

- Update list of functionalities in README.

Changes in version 1.1.7:

- Resolved BiocManager message.

Changes in version 1.1.6:

- Minor fix for Windows unit test.

Changes in version 1.1.5:

- New panel colors.

- Control arguments to custom panels through action buttons.

- Distinguish visible from active arguments for custom panels.

Changes in version 1.1.4:

- Split ?defaults help page by panel type.

- Generalized support for custome data plots and statistics tables.

Changes in version 1.1.3:

- Add new _Sample assay plot_ panel type.

- Extend documentation.

- Split vignette into three: basic, advanced, ExperimentColorMap.

- Fix initialization of reduced dimensions with a single plot axis
  choice.

- Substitute discouraged use of sapply.

- Moved roxygen importFrom instructions closer to the relevant code.

- Increase unit test coverage.

- Consistent use of "colormap" through the package.

- Update installation instructions.

- Add CITATION file.

- Add Figure 1 of article in README.

Changes in version 1.1.2:

- Enable faceting by row and column, with appropriate updates to brush
  and lasso.

- Enable shaping on data points.

- Minor fix of jitter for violin and square plots.

- INTERNAL: Enable storage of additional plot.data beyond X and Y in
  all.coordinates.  See constant .allCoordinatesNames.  Necessary for
  correct behaviour of brushes on faceted plots.

Changes in version 1.0.1:

- Rename feature expression plots to feature assay plots, for
  generality.

[IsoCorrectoR](https://bioconductor.org/packages/IsoCorrectoR)
------------

Changes in version 0.1.27:

- prepare for Biocoductor submission

[IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR)
---------------------

Changes in version 1.3.10:

- A problem was fixed with the isoformSwitchTestDEXSeq() which could
  cause continuous co-variables to interpreted as discrete
  co-variables.

- importRdata() was updated to be a bit more versetile with regards to
  accepting isoform_ids as row.names.

- Both isoformSwitchTestDRIMSeq() and isoformSwitchTestDRIMSeq() was
  updated so the the resulting "isoformSwitchAnalysis" entry in the
  switchAnalyzeRlist also contains results with p-values set to NA. (NA
  filter removed).  Furthemore the interpretation of design matrixes
  with regards to continous or discrete variables was improved.

- The vignette was updated all around including the FAQ sections:

- "What Quantification Tool(s) Should I Use?"

- "What constitue an independent biological replicate?"

- An error message was corrected to give the rigth error

- various small updates

Changes in version 1.3.9:

- Update to namespace to fix 1.3.8 update of importCufflinksFiles

- Update to vignette to fix header

Changes in version 1.3.8:

- Update to importCufflinksFiles to make it faster and more robust.

Changes in version 1.3.7:

- isoformSwitchTestDEXSeq() was updated to use testForDEU instead of
  nbinomLRT as now reccomended by the authors.

Changes in version 1.3.5:

- One-line summary: Improved robustness, usability and speed

- Main changes:

- isoformSwitchTestDEXSeq() is introduced as the new default test as it
  is a more robust and much more reliable test for differential isoform
  usage.

- The original isoformSwitchTest() is decommissioned due to it being
  inferior to both isoformSwitchTestDEXSeq() and
  isoformSwitchTestDRIMSeq() in most aspects.

- importIsoformExpression() now also support import of StringTie
  quantifications.

- updates that allows for better handling of Ensemble data.

- updates throughout the R package making IsoformSwitchAnalyzeR (much)
  faster and more reliable.

- Specifically the changes in inlcuded functions are:

- isoformSwitchTestDEXSeq() is introduced as the default switch isoform
  switch test function

- I handles the False Discovery Rate much better

- It allows for batch corrected effect size estimation

- It is a good deal faster (for smaller sample sizes)

- isoformSwitchTestDRIMSeq() was updated to handle continous
  co-variates.

- isoformSwitchTest() has been removed from the package since it is
  obsolte.

- The importRdata() now:

- Allows for import via either replicate abundance or replciate count
  data (or both - which is highly reccomended).

- These changes were reflected in createSwitchAnalyzeRlist()

- Test for full rank of experimental design

- The importCufflinksCummeRbund() and importCufflinksFiles() now also
  extract and replicate isoform abundance estimates.

- The functions importRdata(), importCufflinksCummeRbund() and
  importCufflinksFiles()

- Calculates isoform fractions based on the replicate isoform fraction
  matrix (instead based on average isoform and gene expression)
  providing more accurate estimataes.

- Was uptimized so they are more streamlined and faster.

- importGTF (also used by importRdata() ) was updated to handle the
  problems with version numbering in amongst other Ensembl data.

- To support the batch correction feature in isoformSwitchTestDEXSeq()
  the subsetSwitchAnalyzeRlist() function was modified so when
  subsetting in the the exon entry of the switchAnalyzeRList, as well
  as any replicate matrix entry (counts, abundances or isoform
  fractions), all isoforms from genes where at least one isoform passed
  the filters are kept.

- The isoformToIsoformFraction() - a general purpose function for
  calculateing Isoform Fraction (IFs) from isoform expression - are
  introduced

- The isoformToGeneExp() function was updated to be true general
  purpose (less stringent about data formating) and thanks to a
  tidyverse solution to the central problem is now between 2x-10x
  faster than previously (and becomses faster as the large the datasets
  are)

- createSwitchAnalyzeRlist() was updated to

- handle replicate data

- fix condition name problems

- test for full rank of design

- importIsoformExpression() was updated to:

- Support StringTie data.

- Perform the inter-library normalization after a lenient expression
  cutoff have beeen applied (to remove most very lowly expressed
  isoforms).

- Now uses the "scaledTPM" instead of "lengthScaledTPM" tximport option
  when imporitng with countsFromAbundance=TRUE The ignoreAfterBar
  argument from tximport() is now also supported.

- We introduce the removeAnnoationData() function which eables removal
  of biological sequence and/or the replicate quantification data from
  a switchAnalyzeRlist threby significantly removing the size.

- The default on the IFcutoff in switchPlot() and
  switchPlotTopSwitches() was updated to 0.05 which should result in
  cleaner plots (meaningisoforms only contributing minimally to the
  parent gene expression are now omitted from plot).

- Specifically the package maintenance changes are:

- All around speed improvements mainly due to updates regarding two
  bottelnecks:

- stringr::str_c replaces paste0 since it is up to 10x faster on
  data.frames

- dplyr::inner_join() or dplyr::left_join() have replaced most
  base::merge() opperations since since they are up to 10x faster.

- All documentation and examples are now based on Salmon data.
  Cufflinks is shown as a special case.

- For this switch new example data was included in the package.

- Directy suppor of Cufflinks/Cuffdiff files via the cummeRbund R
  package (via the importCufflinksCummeRbund function) have been
  removed. Use importCufflinksFiles() instead.

- analyzeSignalP, analyzePFAM, analyzeCPAT now better handles empty
  files.

- All documentation regarding PFAM was updated to use EBI's homepage
  (and their restrictions).

- Updated package title to reflect the introduction of the alternative
  splicing module

- A requirement for tximport >= 1.8.0 was introduced (due to problems
  with importing from RSEM in previous versions)

- Highligting that import of GTF files can be done from both unziped
  and gziped gtf files.

- Updated NEWS file to follow bioconductor style guideline

- Genral update to support condition (and covariate) names compatible
  with model building in R.

- All general support functions (potentially) used more than once place
  were moved to tool.R and names were streamlined.

- Various updates in vignette to reflect all changes desribed above as
  well as update of installation instructions.

- Various updates to input testing to catch commonly occuring problems.

- Correction of loads of spelling mistakes kindely pointed out by
  @afonsoguerra - thanks!

[isomiRs](https://bioconductor.org/packages/isomiRs)
-------

Changes in version 1.9.5:

FIX

- Fix NEWS format

Changes in version 1.9.4:

- Support object creation from rawData of a current IsomirDataSeq
  object.

Changes in version 1.9.3:

- Fix warnings in dosc.

Changes in version 1.9.2:

- Fix bug that quantify wrongly the reference sequence.

Changes in version 1.9.1:

MAJOR

- Improve clustering in isoNetwork plot.

- Use varianceStabilizingTransformation to normalize counts.

- Use a less common column as ID for samples in isoPLot fns.

- Add updateIsomirDataSeq to be compatible with previous versions.

- Adapt all functions to new object. Fix documentation.

- Reduce object size and structure

[kmknn](https://bioconductor.org/packages/kmknn)
-----

Changes in version 1.0.0:

- New package kmknn, for k-means-based k-nearest neighbor detection.

[limma](https://bioconductor.org/packages/limma)
-----

Changes in version 3.38.0:

- New function plotExonJunc() to plot results from diffSplice().

- New function logsumexp().

- New argument hl.col for volcanoplot(), allowing users to specify the
  color for the gene names when highlight > 0.

- barcodeplot() no longer assumes that 'statistic' has unique names.
  Previously it returned an error if names(statistic) contained any
  duplicated values.

- The colors "blue", "red" and "yellow" used by coolmap() changed to
  "blue2", "red2" and "yellow2" when used in a color panel with white.

- goana.Rd now explains more explicitly that p-values are unadjusted
  for multiple testing.

- arrayWeights.Rd now mentions minimum dimensions for expression
  object.

- More advice on how to choose 'lfc' added to the treat() help page.

- Minor bug fix to the mixed p-value from roast() and mroast() when
  set.statistic="floormean".

- Bug fix for cumOverlap(), which was under-counting overlaps in some
  cases.

[LRBaseDbi](https://bioconductor.org/packages/LRBaseDbi)
---------

Changes in version 0.99.17:

- Error is occured when the OS is windows?

Changes in version 0.99.16:

- inst/LRBasePkg-template/man/@PKGNAME@.Rd is modified

Changes in version 0.99.15:

- Reviewed

Changes in version 0.99.14:

- vignette is modified

Changes in version 0.99.13:

- vignette is modified

Changes in version 0.99.12:

- vignette is modified

Changes in version 0.99.11:

- vignette is modified

Changes in version 0.99.10:

- Removed inst/doc and added vignettes

Changes in version 0.99.9:

- Removed vignettes

Changes in version 0.99.8:

- Added inst/doc

Changes in version 0.99.7:

- Added inst/doc

Changes in version 0.99.6:

- rm -rf inst/doc

Changes in version 0.99.5:

- Added VignetteBuilder and VignetteEngine

- The vignette set as BiocStyle

Changes in version 0.99.4:

- The vignette-related files (.Rmd/.html) in vignettes directory are
  moved to inst/doc direcotry

Changes in version 0.99.3:

- The vignette-related files (.Rmd/.html) in inst/doc directory are
  removed

Changes in version 0.99.2:

- The vignette-related files (.Rmd/.html) in inst/doc directory are
  removed

Changes in version 0.99.1:

- The vignette-related files (.Rmd/.html) in inst/doc directory are
  removed

Changes in version 0.99.0:

- Package released

[maftools](https://bioconductor.org/packages/maftools)
--------

Changes in version 1.8.0:

NEW FUNCTIONS

- OncogenicPathways - Perform enrichment for known oncogenic pathways
  from TCGA studies.

- PlotOncogenicPathways - Plots OncogenicPathways results

- drugInteractions - Drug gene interactions from DGIB database.

SIGNIFICANT USER-LEVEL IMPROVEMENT

- trinucleotideMatrix functions now works with BSgenomes instead of
  time consuming fasta files

- rainfallPlot now detects Kataegis based on a greedy algorithm

[MAGeCKFlute](https://bioconductor.org/packages/MAGeCKFlute)
-----------

Changes in version 1.1.9:

- Add parameter 'pathway_limit' / 'limit' / 'gmtpath' in FluteRRA,
  FluteMLE, and all enrichment functions which enable users to
  customize gene sets for enrichment analysis.

- Remove DAVID and GOstats, which are not recommended.

- Add enrichment score in enrichment results from all algorithms.

Changes in version 1.1.6:

- Add functions to plot figures in NatureProtocol manuscript.

- Shorten the check time.

- Remove null figures in pathview part.

Changes in version 1.1.1:

- Speed up the enrichment analysis.

- Release memory in time.

[matter](https://bioconductor.org/packages/matter)
------

Changes in version 1.7.6:

NEW FEATURES

- Added 'apply' methods for 'sparse_mat' and 'virtual_mat'

- Added 'vm_used' function to exported internal utilities

- Can infer length of 'matter_vec' from paths when missing

- Can coerce 'matter_list' to 'matter_matc' or 'matter_matr' if all
  elements of the list are the same length

SIGNIFICANT USER-VISIBLE CHANGES

- Initializing a matter object with data to an existing file no longer
  results in a warning if 'filemode' is supplied

BUG FIXES

- Fixed bug where subsetting 'sparse_mat' objects would pull
  'matter_list' key-value pairs into memory

- Subsetting 'sparse_mat' objects with out-of-bounds subscripts when
  'drop=NULL' is now an error

Changes in version 1.7.5:

NEW FEATURES

- Added 'as.matrix' methods for 'sparse_mat' and 'virtual_mat'

Changes in version 1.7.4:

NEW FEATURES

- Added coercion from 'matter_mat' to 'matter_list'

Changes in version 1.7.3:

NEW FEATURES

- Added 'rep_vt' class for virtual replicated vectors

Changes in version 1.7.2:

NEW FEATURES

- Updated installation instructions for BiocManager

- Setting 'sparse_mat' keys also updated nrows/ncols

BUG FIXES

- Fixed mem() to reflect R >= 3.5 gc() function

Changes in version 1.7.1:

BUG FIXES

- Fixed file.exists() bug when length(paths) > 1

[meshes](https://bioconductor.org/packages/meshes)
------

Changes in version 1.7.1:

- add citation <2018-05-25, Fri>

[MetCirc](https://bioconductor.org/packages/MetCirc)
-------

Changes in version 1.11.3 (2018-10-22):

- adjust NEWS file to new format according to ?news: o entries are
  grouped according to version, with version header “Changes in
  version” at the beginning of a line, followed by a version number,
  optionally followed by an ISO 8601 format date, possibly
  parenthesized o entries may be grouped according to category, with a
  category header (different from a version header) starting at the
  beginning of a line o entries are written as itemize-type lists,
  using one of o, *, - or + as item tag. Entries must be indented, and
  ideally use a common indentation for the item texts

Changes in version 1.11.2 (2018-10-15):

- remove Makefile from directory vignettes/

- check if package passes R CMD build and R CMD check without any error
  messages and vignette can be run without any errors

[methInheritSim](https://bioconductor.org/packages/methInheritSim)
--------------

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- New citation referring to the associated published article in Nucleic
  Acids Research.

[methylInheritance](https://bioconductor.org/packages/methylInheritance)
-----------------

Changes in version 1.5.1:

SIGNIFICANT USER-VISIBLE CHANGES

- New citation referring to the associated published article in Nucleic
  Acids Research.

[methylKit](https://bioconductor.org/packages/methylKit)
---------

Changes in version 1.7.10:

IMPROVEMENTS AND BUG FIXES

- update man: update object and function descriptions to resolve the
  'file link treated as topic' warning under windows

Changes in version 1.7.9:

IMPROVEMENTS AND BUG FIXES

- bug fix: make internal function .checkdbdir more error proof

Changes in version 1.7.8:

IMPROVEMENTS AND BUG FIXES

- internal applyTbxByOverlap(): did not work correctly for return.type
  set to "data.frame" or"data.table", due to missing argument.
  https://github.com/al2na/methylKit/issues/131

Changes in version 1.7.7:

IMPROVEMENTS AND BUG FIXES

- unite(): skip checking if destranded tabix file exist

Changes in version 1.7.6:

NEW FUNCTIONS AND FEATURES

- methSeg(): introduce parameter `initialize.on.subset` to subset data
  for initialization of mixture modeling; update description; add tests

Changes in version 1.7.4:

IMPROVEMENTS AND BUG FIXES

- update link to test-file for methSeg() function

Changes in version 1.7.3:

IMPROVEMENTS AND BUG FIXES

- selectByOverlap(): update description; refine method signatures to
  only support GRanges as range argument; update NAMESPACE to import
  subjectHits() method from S4Vectors

Changes in version 1.7.2:

NEW FUNCTIONS AND FEATURES

- New constructor method methylRawList() can be used to combine list of
  methylRaw objects into a methylRawList

IMPROVEMENTS AND BUG FIXES

- fix bug in methSeg: when joinSegments was activated, diagnostic plot
  would always be plotted

- fixes in selectByOverlap() function: update description, to show that
  any methylKit object (tabix or not) can be used with it ; fixed
  broken method after @subjectHits was not available anymore

- fix error in .checkTabixFileExists function, that lead to overwriting
  of files

[MetNet](https://bioconductor.org/packages/MetNet)
------

Changes in version 0.99.24 (2018-10-22):

- adjust NEWS file to new format according to ?news: o entries are
  grouped according to version, with version header “Changes in
  version” at the beginning of a line, followed by a version number,
  optionally followed by an ISO 8601 format date, possibly
  parenthesized o entries may be grouped according to category, with a
  category header (different from a version header) starting at the
  beginning of a line o entries are written as itemize-type lists,
  using one of o, *, - or + as item tag. Entries must be indented, and
  ideally use a common indentation for the item texts

Changes in version 0.99.23 (2018-10-16):

- improve createStructuralAdjacency function

Changes in version 0.99.20 (2018-08-06):

- replace psych::corr.test by WGCNA::corAndPvalue to improve speed

Changes in version 0.99.19 (2018-07-26):

- print message when model calculation in createStatisticalAdjacency is
  finished

Changes in version 0.99.18 (2018-07-19):

- set rfPermute.formula to rfPermute.default in order to use num.cores

Changes in version 0.99.17 (2018-07-19):

- set num.cores to 1 in test_statistical for randomForest

Changes in version 0.99.15 (2018-07-16):

- do not import stabsel from stabs

Changes in version 0.99.14 (2018-07-14):

- use BiocManager instead of BiocLite for installation

Changes in version 0.99.13 (2018-07-12):

- do not export functions threeDotsCall and addToList

Changes in version 0.99.12 (2018-07-12):

- use BiocStyle package for vignette

- remove Makefile

- use BiocParallel instead of parallel, for instance use bplapply
  instead of mclapply

Changes in version 0.99.11 (2018-07-03):

- speed up function rtCorrection by vectorizing

Changes in version 0.99.10 (2018-07-03):

- speed up function createStructuralAdjacencyMatrix by vectorizing

Changes in version 0.99.9 (2018-06-26):

- change 1:... to seq_len()

Changes in version 0.99.8 (2018-06-26):

- implement function rtCorrection

Changes in version 0.99.6 (2018-06-13):

- use camelCaps for functions

- use no spaces between '=' and named arguments

- fix typo in lasso function

Changes in version 0.99.5 (2018-06-12):

- change 1:... to seq_len()

Changes in version 0.99.4 (2018-06-12):

- change 1:... to seq_len()

Changes in version 0.99.3 (2018-06-11):

- remove sum check for correlation

Changes in version 0.99.2 (2018-06-11):

- require R version >= 3.5

Changes in version 0.99.1 (2018-06-11):

- remove bugs that there are no WARNINGS and ERRORs when running R CMD
  check and R CMD BiocCheck

- reduce file size of peaklist_example.RData

- submit to Bioconductor

Changes in version 0.99.0 (2018-05-14):

- implement functionality to calculate statistical models of
  correlation (Pearson, Spearman), LASSO, Random Forest, Context
  likelihood or relatedness network algorithm, algorithm for the
  reconstruction of accurate cellular networks, constraint-based
  structure learning algorithm

- implement the function create_statistical_network to calcululate the
  consensus matrix from the different statistically-infered networks

- implement the function create_structural_network to calculate
  molecular weight differences and create a network

- implement the function combine_structural_statistical to combine the
  structurally-derived and statistically-derived network

- implement model partial and semi-partial pearson/spearman correlation
  using the ppcor package

[minfi](https://bioconductor.org/packages/minfi)
-----

Changes in version 1.27:

- v1.27.2 Fix bug in preprocessQuantile() that arose when checking
  input for previous preprocessing method. Thanks to @DelnazR for the
  report (<URL: https://github.com/hansenlab/minfi/issues/165>).

- v1.27.3 Fixed bug related to switch A and B for SNPs of type I when
  using convertArray / combineArrays.  Reported by Jenny van Dongen.

- v1.27.3 Fixed error in dmpFinder.

- Added preliminary support for HorvathMammalMethylChip40.

[miRSM](https://bioconductor.org/packages/miRSM)
-----

Changes in version 0.99.28:

- Update module_biclust function <2018-10-23, Thue>

Changes in version 0.99.27:

- Update module_NMF function <2018-10-20, Sat>

Changes in version 0.99.26:

- Update module_biclust function <2018-10-18, Thus>

Changes in version 0.99.25:

- Update module_NMF function <2018-10-12, Fri>

Changes in version 0.99.24:

- Update miRSM.R <2018-10-10, Wed>

Changes in version 0.99.23:

- Update Vignettes <2018-08-31, Fri>

Changes in version 0.99.15-0.99.22:

- Update R codes, DESCRIPTION, NEWS, README <2018-08-27, Mon>

Changes in version 0.99.8-0.99.14:

- Update miRSM.Rmd <2018-08-07, Thue>

Changes in version 0.99.4-0.99.7:

- Update runnable examples <2018-08-06, Mon>

Changes in version 0.99.3:

- Update src <2018-08-03, Fri>

Changes in version 0.99.2:

- Update module_biclust function <2018-08-03, Fri>

Changes in version 0.99.1:

- Change the type of input data into a SummarizedExperiment object.
  Polish the file miRSM.R. <2018-08-03, Fri>

Changes in version 0.99.0:

- This is the first version of miRSM package. If any bugs, please let
  me know. Contact Email: zhangjunpeng_411@yahoo.com <2018-08-03, Fri>

[miRsponge](https://bioconductor.org/packages/miRsponge)
---------

Changes in version 1.7.5:

- Update miRsponge.Rmd <2018-08-25, Sat>.

Changes in version 1.7.4:

- Update netModule function <2018-08-10, Fri>.

Changes in version 1.7.3:

- Update netModule function <2018-07-31, Tues>.

Changes in version 1.7.2:

- Update netModule function <2018-07-31, Tues>.

Changes in version 1.7.1:

- Update netModule function <2018-07-29, Sun>.

Changes in version 1.7.0:

- Update netModule function <2018-07-27, Fri>.

[missMethyl](https://bioconductor.org/packages/missMethyl)
----------

Changes in version 1.15.2:

- Updated getMappedEntrezIDs, gometh and gsameth to to speed up
  execution by taking the array annotation in as an optional argument.

- missMethyl now uses the latest
  IlluminaHumanMethylationEPICanno.ilm10b2.hg19 annotation by default
  for EPIC arrays.

Changes in version 1.15.1:

- Added getAdjusted function for extracting RUVm adjusted data for
  visualisation purposes

- Updated vignette to demonstrate use of getAdjusted function

- Vignette now includes an example of how to handle cases with RUVm
  where number of samples is greater than number of Illumina negative
  controls


[motifcounter](https://bioconductor.org/packages/motifcounter)
------------

Changes in version 1.5.4:

- Accept DNAString wherever DNAStringSet was accepted.

Changes in version 1.5.3:

- Fast motif matching Various improvements to improve motif matching
  speed were implemented. 1. Precomputed position weights to avoid the
  invokation of the log function during scanning. 2. threshold argument
  added to avoid having to determine the score distribution for each
  sequence individually. 3. Stop as soon as the motif score cannot or
  is certainly exceeding the threshold.

Changes in version 1.5.2:

- scoreSequence, scoreHistogram, scoreProfile motifHits,
  motifHitProfile can be used with N-containing sequences. Positions at
  which the motif overlaps with N's will be NaN.

- Background can be computed from a DNAString in addition to a
  DNAStringSet object.

Changes in version 1.5.1:

- Background can be computed with sequences that contain N's.

[MotifDb](https://bioconductor.org/packages/MotifDb)
-------

Changes in version 1.24:

NEW FEATURES

- query method now flexible, with "andStrings", "orStrings",
  "notStrings" parameters The previous usage style is still supported.
  See man page.

- associateTranscriptionFactors (with motifs) substantially faster

[motifStack](https://bioconductor.org/packages/motifStack)
----------

Changes in version 1.25.2:

- fix the bug when plot "others" type of motif.

Changes in version 1.25.1:

- Fix a bug for y-axis plot.

[MPRAnalyze](https://bioconductor.org/packages/MPRAnalyze)
----------

Changes in version 0.99.0 (2018-08-07):

- Submitted to Bioconductor

[MSnbase](https://bioconductor.org/packages/MSnbase)
-------

Changes in version 2.7.12:

- Fix warnings on windows (see #371) <2018-10-26 Fri>

- Add parameter ppm to consensusSpectrum and meanMzInts (see #373 for
  details) <2018-10-26 Fri>

Changes in version 2.7.11:

- Change default for `timeDomain` in `combineSpectra` and
  `combineSpectraMovingWindow` to `FALSE` <2018-10-18 Thu>

- Add new spectra combination function `consensusSpectrum` <2018-10-24
  Wed>

- Amend plot,Spectrum 1 and 2 (see #369)

Changes in version 2.7.10:

- Methods for Spectra class <2018-10-15 Mon>

Changes in version 2.7.9:

- Import rather than depend on BiocParallel <2018-10-15 Mon>

- Fix failing test on Windows (requiring normalizePath) <2018-10-15
  Mon>

Changes in version 2.7.8:

- MGF exporter gets a new `addFields` argument (see PR #362)
  <2018-10-12 Fri>

- New `Spectra` (`SimpleList` of `Sepctrum` objects) (see PR #361)
  <2018-10-13 Sat>

Changes in version 2.7.7:

- Fix unit tests (issue #360, wrong MS OBO CV terms for data smoothing
  (MS:1000592) and baseline correction (MS:1000593))

Changes in version 2.7.6:

- Fix wrong MS OBO CV terms for data smoothing (MS:1000592) and
  baseline correction (MS:1000593).

Changes in version 2.7.5:

- Add a note about parallel processing in vignette (see #356 for
  background) <2018-09-04 Tue>

Changes in version 2.7.4:

- Fix filterMz for spectra with no non-NA intensities in m/z range (see
  #355) <2018-08-08 Wed>

Changes in version 2.7.3:

- Fix bug in robust summary (see PR #349) <2018-07-28 Sat>

- Fix failing unit test <2018-07-28 Sat>

Changes in version 2.7.2:

- Handle files without any spectra - see #342 <2018-05-15 Tue>

- New `mergeFeatureVars` and `expandFeatureVars` functions <2018-05-30
  Wed>

- Update plot,Spectrum methods to match the tolerance and relative
  arguments (see #350) <2018-06-29 Fri>

Changes in version 2.7.1:

- Version bump to force new vignette build

Changes in version 2.7.0:

- New devel version for Bioc 3.8 # MSnbase 2.6

[MSnID](https://bioconductor.org/packages/MSnID)
-----

Changes in version 1.15.1:

- mzR now added scan number(s) into the table representation of
  mzIdentML object. As a results it caused an error in my unit test
  checking if the file reads properly. Fixed this check with updated
  hash.

[msPurity](https://bioconductor.org/packages/msPurity)
--------

Changes in version 1.6.1:

- Bug fix. For pos/neg switching acquisition two files are can be
  generated when converting from RAW to mzML (1 for pos, 1 for neg).
  The resulting files retention time scans were not being tracked
  properly in msPurity in these cases. This is now fixed. Thanks to
  Julien (https://github.com/jsaintvanne) for spotting the bug.


[MSstatsTMT](https://bioconductor.org/packages/MSstatsTMT)
----------

Changes in version 0.99.0 (2018-09-21):

- Submitted to Bioconductor

[multiClust](https://bioconductor.org/packages/multiClust)
----------

Changes in version 1.11.2:

- Updated vignette Rmd file to address bug in getGEO to obtain
  expression data

Changes in version 1.11.1:

- Updated NEWS file format to reflect changes by version

[mzR](https://bioconductor.org/packages/mzR)
---

Changes in version 2.15.5:

- Fix bug #181

Changes in version 2.15.4:

- Use new dependency ncdf4 for netCDF reading, removes a lot of build
  hassles with old libnetcdf-dev linking.

- specParams returns a numeric scan.number.s.

Changes in version 2.15.3:

- Adds MS-GF+ information such as Scan Time and a more reliable Scan
  Number, contributed by FarmGeek4Life (see PR #174).

Changes in version 2.15.2:

- Add header column ionMobilityDriftTime to report the corresponding CV
  parameter (issue https://github.com/sneumann/mzR/issues/44).

- Ensure ion injection time is always reported in milliseconds.

- Replace BiocInstaller::biocLite with BiocManager::install (by Bioc
  core)

Changes in version 2.15.1:

- Fix typo (see https://github.com/sneumann/mzR/pull/162)

- New .hasSpectra and .hasChromatograms private function (see
  https://github.com/lgatto/MSnbase/issues/343)

- Fix bug in score when more cvParams than expected are read - see
  https://github.com/sneumann/mzR/issues/136 <2018-05-26 Sat>

Changes in version 2.15.0:

- New Bioc devel version

[NBSplice](https://bioconductor.org/packages/NBSplice)
--------

Changes in version 0.99.4:

- CODE: Modification in the fitModel function to correct the
  estimations with contrasts specified by the user

[NormalyzerDE](https://bioconductor.org/packages/NormalyzerDE)
------------

Changes in version 0.99.24:

- Various input checks added

Changes in version 0.99.22:

- Input parsing corrected so that a clear error message is provided
  when the design matrix samples doesn't match the data matrix header

Changes in version 0.99.21:

- Documentation example fixes

Changes in version 0.99.20:

- Reactivate quantile normalization

Changes in version 0.99.19:

- Empty annotation no longer leads to crash when writing output

Changes in version 0.99.17:

- Very minor fix where function wasn't retrieved properly from
  SummarizedExperiment

Changes in version 0.99.16:

- Crash when feeding SummarizedExperiment as input is fixed

- If full rows contains only NA values they are filtered and a warning
  is generated

Changes in version 0.99.15:

- Corrected expected output from vapply in correlation matrix (where
  previous varying lengths caused crash)

- Clearer error message when providing invalid sample/group conditions
  to statistics module

- Issue where several identical RT-values caused crash in RT-slicing
  fixed


[omicplotR](https://bioconductor.org/packages/omicplotR)
---------

Changes in version 1.1.2:

- bug fix for CZM in RAB plots. if not zeros, plots would fail.

- updated install instructions to use `BiocManager`

Changes in version 1.1.1:

- added slider inputs to change size and opacity of sample names for
  biplot and coloured biplot

- changed UI for filtering page

- added downloadable scripts for PCA biplots, dendrogram/barplots, and
  effect plots (simple effect plot script).


[OncoSimulR](https://bioconductor.org/packages/OncoSimulR)
----------

Changes in version 2.11.1:

- robustify test.fixation.R, Local max, tolerance

[onlineFDR](https://bioconductor.org/packages/onlineFDR)
---------

Changes in version 0.99.7:

MODIFICATIONS

- randomisation of batches now implemented via the randBatch helper
  function

- updated references

- updated vignette

Changes in version 0.99.5:

MODIFICATIONS

- replace date field with hard-coded date in vignette

- updating data.frame rows by group now uses 'split-apply-combine'

- vectorise part of the LOND function

- replace 1:N with seq_len(N)

- applied consistent formatting & indentation

- added new tests

[oposSOM](https://bioconductor.org/packages/oposSOM)
-------

Changes in version 2.0.1:

- Colored edges in sample correlation networks & minor fixes.

Changes in version 2.0.0:

- Great performance boost by implementation of parallel SOM
  calculation. "som" package dependence removed.

[ORFik](https://bioconductor.org/packages/ORFik)
-----

Changes in version 1.1.12:

SIGNIFICANT USER-VISIBLE CHANGES

- The orf finding function now find the longest orf per stop codon if
  you set longestORF = TRUE in findORFS, findMapORFs and findORFsFasta

[PathoStat](https://bioconductor.org/packages/PathoStat)
---------

Changes in version 1.8.1:

- Update data uploading

- Modify visualization

- Add biomarker tab

[pathVar](https://bioconductor.org/packages/pathVar)
-------

Changes in version 1.11.2 (2018-06-29):

- Updated NEWS file of pathVar

- Updated maintainer email address

Changes in version 1.11.1 (2018-05-15):

- Updated NEWS file of pathVar

[pcaExplorer](https://bioconductor.org/packages/pcaExplorer)
-----------

Changes in version 2.8.0:

NEW FEATURES

- PCA plots now are correctly generated with fixed coordinates

[PGA](https://bioconductor.org/packages/PGA)
---

Changes in version 1.11.4:

- Fix bug in function PrepareAnnotationRefseq2.R

Changes in version 1.11.2:

- Add a option 'tabFile' in the 'dbCreator' function for the support of
  Hisat2 "splicesites.tab" in junction construction. As HISAT2 is a
  successor to both HISAT and TopHat2, we recommend that users switch
  to 'tabFile' from HISAT2 instead of 'bedFile' from Tophat2.

[phantasus](https://bioconductor.org/packages/phantasus)
---------

Changes in version 1.1.6:

- Minor bug fixes

Changes in version 1.1.5:

- PCA plot uses annotation color scheme

- Rename annotation column implemented

- log2(1 + x) adjust tool added

- Calculated annotation revamped

- Annotate dataset revamped and moved to file

- Session syncing improved.

- Sweep adjustment implemented

Changes in version 1.1.3:

- Option to get vertical GSEA plot

- Showing heatmap along the GSEA plot

- Changed filter scheme

- Better session synching between fronend and backend

Changes in version 1.1.2:

- Moved to protobuf 3.4 in docker (supports messages up to 2GB)

- Set exact match to be default

- Set row profile as the default chart

Changes in version 1.1.1:

- Fixes in PCA plot

- GSEA plotCHANGES IN VERSION 0.99.34

- Detecting conditions and replicates in GEO data

[piano](https://bioconductor.org/packages/piano)
-----

Changes in version 1.22.0:

DOCUMENTATION

- Update installation instructions in the vignette

Changes in version 1.20.1:

DOCUMENTATION

- Updates to GitHub landing page (README.md)

[Pigengene](https://bioconductor.org/packages/Pigengene)
---------

Changes in version 1.7.2 (2018-05-22):

General

- The version of the package C50 is now required to be at least 0.1.2,
  which exports the as.party.C5.0() function.

[plotGrouper](https://bioconductor.org/packages/plotGrouper)
-----------

Changes in version 0.99.37:

New features

- Added ability to display adjusted p value

Changes in version 0.99.24:

New features

- Added unit testing for the shiny app.

Changes in version 0.99.15:

New features

- Can now load .csv and .tsv files.

Changes in version 0.99.14:

Minor improvements and bug fixes

- Using `vapply()` instead of `sapply()` for safer code.

[polyester](https://bioconductor.org/packages/polyester)
---------

Changes in version 1.99.3:

- NB function now exported

- note that version 1.99.3 on GitHub was version 1.1.0 on Bioconductor.

Changes in version 1.99.2:

- bug fix in fragment generation (last 2 bases of transcript were never
  sequenced)

[PowerExplorer](https://bioconductor.org/packages/PowerExplorer)
-------------

Changes in version 1.1.1:

- Bug fixes for adding missing values in simulations.

Changes in version 1.1.0:

- The release version.

[primirTSS](https://bioconductor.org/packages/primirTSS)
---------

Changes in version 0.99.7:

- Initial release of 'primirTSS' package

[pRoloc](https://bioconductor.org/packages/pRoloc)
------

Changes in version 1.21.9:

- Fix type in vignette header <2018-09-18 Tue>

- Fix bug in plot method for ThetaRegRes object <2018-09-24 Mon>

Changes in version 1.21.8:

- Add an `fcol` argument to `plotDist` to plot and colour all profiles
  <2018-08-09 Thu>

Changes in version 1.21.7:

- Use BiocManager in vignette

- Fix bug in plot2D: pass ... to hexbin <2018-08-02 Thu>

Changes in version 1.21.6:

- Use BiocManager in installation instructions

Changes in version 1.21.5:

- Added new section in Bayesian spatial proteomics vignette detailing
  mcmc output processing <2018-07-07 Sat>

Changes in version 1.21.4:

- Fix bugs in tagmMcmcPredict, where fcol was ignored <2018-06-05 Tue>

- Order vignettes by prefixing the files with numbers <2018-06-05 Tue>

Changes in version 1.21.3:

- New TAGM-MCMC generative model, contributed by Oliver Crook
  <2018-05-18 Fri>

Changes in version 1.21.2:

- Version bump for BiocStyle update: Vignette needed to be rebuilt to
  have bug fixed in BiocStyle footnote rendering.

Changes in version 1.21.1:

- Fix bug in higlightOnPlot with missing fcol (see #105) <2018-05-03
  Thu>

- New TAGM-MAP generative model, contributed by Oliver Crook
  <2018-05-18 Fri>

- New `plotEllipse` function to visualise and assess TAGM models
  <2018-05-18 Fri>

[PureCN](https://bioconductor.org/packages/PureCN)
------

Changes in version 1.12.0:

NEW FEATURES

- normalDB does not need input normal coverage files anymore after
  creation (so the resulting normalDB.rds file can be moved)

- base quality filtering can be turned off by setting min.base.quality
  to 0 or NULL

- possible to change the POP_AF info field name

- possible to change POP_AF cutoff to set a high germline prior

- possible to change min.cosmic.cnt and max.homozygous.loss in PureCN.R

- set number of cores in PureCN.R (thanks Brad)

SIGNIFICANT USER-VISIBLE CHANGES

- renamed reptimingbinsize to reptimingwidth in IntervalFile.R, added
  this feature to preprocessIntervals

- clarified "targets" vs. "intervals"; whenever something affects both
  on-target and off-target, it is now called "intervals". When only
  targets, e.g. in annotateTargets, "targets" was kept.

- made gc.gene.file defunct

- new default for min.cosmic.cnt = 6 (instead of 4)

BUGFIXES

- catch various input problems and provide better error messages
  instead of crashing

- stranded input BED files do not cause problems anymore

- fixed a bug when only a single local optimum was tested (happens only
  when users copy the examples that restrict the search speach to avoid
  long runtimes)

- added missing QC flag to predictSomatic VCF annotation

[qPLEXanalyzer](https://bioconductor.org/packages/qPLEXanalyzer)
-------------

Changes in version 0.99.0:

- initial version with the following functions implemented: +
  convertToMSnset + summarizeIntensities + normalizeQuantiles +
  normalizeScaling + groupScaling + rowScaling + regressIntensity +
  computeDiffStats + getContrastResults + assignColours + intensityPlot
  + intensityBoxplot + peptideIntensityPlot + pcaPlot + maVolPlot +
  corrPlot + rliPlot + hierarchicalPlot + plotMeanVar + coveragePlot

[qsea](https://bioconductor.org/packages/qsea)
----

Changes in version 1.7.4:

- Reformated NEWS file

Changes in version 1.7.3:

- Adapted vignette to BiocManager installation process

Changes in version 1.7.2:

- Bugfixes: - addCoverage fixed fragment size sd for chromosomes with 1
  read - addCoverage fixed missing regions in count_matrix for
  chromosomes without any reads

Changes in version 1.7.1:

- Bugfix: addCoverage crashes when there are no reads for
  chromosome/contig


[RCy3](https://bioconductor.org/packages/RCy3)
----

Changes in version 2.2.0:

- New functions to remove duplicate edges - deleteDuplicateEdges -
  deleteSelfLoops

- New node selection function - selectNodesConnectedBySelectedEdges

- New visual style management functions - importVisualStyles -
  deleteVisualStyle - deleteStyleMapping

- New edge bundling function - bundleEdges

- New custom graphics options for nodes - setNodeCustomBarChart -
  setNodeCustomBoxChart - setNodeCustomHeatMapChart -
  setNodeCustomLineChart - setNodeCustomPieChart -
  setNodeCustomRingChart - setNodeCustomLinearGradient -
  setNodeCustomRadialGradient - setNodeCustomPosition -
  removeNodeCustomGraphics

- New filter functions - applyFilter - createColumnFilter -
  createCompositeFilter - createDegreeFilter - getFilterList -
  exportFilters - importFilters

- Improved speed on bulk node and edge property bypasses

- Bug Fixes - selectEdgesConnectingSelectedNodes -- set default by.col
  = 'name' - setEdgeLineWidthMapping -- fixes input type - getGroupInfo
  -- works without collapsing first - getTableColumns -- work with List
  type columns

- For Developers - Updated many functions to properly pass the base.url
  parameter to functions like getNetworkSuid. Please be aware and
  vigilent about this with future development. - Adopted use of
  seq_len(). Please be aware and vigilent. - Replaced all but one case
  of sapply() with vapply().

- Deprecated - Nothing

- Defunct - Previously deprecated functions in v2.0 from older 1.x
  version of the package

[ReactomePA](https://bioconductor.org/packages/ReactomePA)
----------

Changes in version 1.25.1:

- add keyType parameter in viewPathway (2018-09-04, Tue)

[REBET](https://bioconductor.org/packages/REBET)
-----

Changes in version 0.99.0 (2018-07-12):

- Submitted to Bioconductor

[recount](https://bioconductor.org/packages/recount)
-------

Changes in version 1.7.5:

SIGNIFICANT USER-VISIBLE CHANGES

- add_metadata() can now download the recount_brain_v2 data.

Changes in version 1.7.4:

BUG FIXES

- Fix a NOTE about RefManageR.

Changes in version 1.7.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Use BiocManager

Changes in version 1.7.2:

BUG FIXES

- Fix a unit test.

Changes in version 1.7.1:

SIGNIFICANT USER-VISIBLE CHANGES

- rse_tx URLs now point to v2 to reflect recent changes by Fu et al.

[regioneR](https://bioconductor.org/packages/regioneR)
--------

Changes in version 1.13.2:

BUG FIXES

- createRandomRegions ignored the non.overlapping argument. It does
  work now.

Changes in version 1.13.1:

NEW FEATURES

- Revamped toGRanges now accepts genome region descriptions as used by
  UCSC and IGV ("chr9:23000-25000"). It also may take a genome
  parameter and set the genome information of the GRanges accordingly.

[regionReport](https://bioconductor.org/packages/regionReport)
------------

Changes in version 1.15.4:

BUG FIXES

- Fix a bug in the order that was reported and fixed by @bounlu at
  https://github.com/leekgroup/regionReport/pull/9/

Changes in version 1.15.3:

BUG FIXES

- Fixed an issue with DESeq2Exploration.Rmd that affected both
  DESeq2Report and edgeReport. This should also fix the recount
  bioc-release (3.7) and bioc-devel (3.8) branches.

- Fixed a NAMESPACE issue with rmarkdown::html_document and
  BiocStyle::html_document

Changes in version 1.15.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Use BiocManager

Changes in version 1.15.1:

BUG FIXES

- Fix namespace issue in relation to BiocStyle::html_document2

[rGREAT](https://bioconductor.org/packages/rGREAT)
------

Changes in version 1.13.1:

- `plotRegionGeneAssociationGraphs()`: `par(mfrow)` is automatically
  set according to the length of `type`.

[RNASeqR](https://bioconductor.org/packages/RNASeqR)
-------

Changes in version 0.99.25:

- Update vignette rmd file. <2018-10-26 Thu>

Changes in version 0.99.24:

- Change description in DESCRIPTION file. <2018-10-24 Tue>

Changes in version 0.99.23:

- Add show method for 'RNASeqRParam' S4 object and revise vignette.
  <2018-10-16 Tue>

[RnBeads](https://bioconductor.org/packages/RnBeads)
-------

Changes in version 1.99.0:

- RnBeadsDJ: updated documentation and tooltips

- deactivated intersample plots (exploratory module) as default option
  for BS-seq data

- Roxygen documentation updates

- Extended methods descriptions in reports (imputation, filtering,
  differential, ...)

- Miscellaneous bugfixes and performance improvements in RnBeadsDJ,
  imputation, ...

Changes in version 1.13.3:

- added CNV estimation using the GLAD package to QC module

- some bugfixes

- changed 'gender' to 'sex', affected functions and options are
  'rnb.execute.gender.prediction' and 'import.gender.prediction'

Changes in version 1.13.2:

- added qvalue to suggested packages

Changes in version 1.13.1:

- Updated LOLA DB download links

[rols](https://bioconductor.org/packages/rols)
----

Changes in version 2.9.4:

- Fix unit test <2018-10-06 Sat>

Changes in version 2.9.3:

- Update news and pkgdown

Changes in version 2.9.2:

- replace BiocInstaller biocLite mentions with BiocManager

Changes in version 2.9.1:

- Fix bug ancestors function, reported by Christian Holland (see
  https://github.com/lgatto/rols/issues/26 for details) <2018-06-01
  Fri>

Changes in version 2.9.0:

- Bioconductor devel

[ropls](https://bioconductor.org/packages/ropls)
-----

Changes in version 1.13.8:

INTERNAL MODIFICATION

- minor correction in the documentation

Changes in version 1.13.6:

INTERNAL MODIFICATION

- update of package vignette

Changes in version 1.13.4:

INTERNAL MODIFICATION

- fixed bug in unit test

Changes in version 1.13.2:

BUG CORRECTION

- predict method: naming of the predicted Y matrix output columns in
  case of PLS modeling of multiple responses

INTERNAL MODIFICATION

- ropls.Rproj file added for package management with RStudio

[rpx](https://bioconductor.org/packages/rpx)
---

Changes in version 1.17.2:

- Update NEWS and pkgdown site

Changes in version 1.17.1:

- replace BiocInstaller biocLite mentions with BiocManager

Changes in version 1.17.0:

- New version for Bioconductor devel

[Rsamtools](https://bioconductor.org/packages/Rsamtools)
---------

Changes in version 1.33:

NEW FEATURES

- (v 1.33.4, 1.33.7) scanBamFlag() gains isSupplementaryAlignment
  support.

BUG FIXES

- (v 1.33.1) Do not try to grow NULL (not-yet-encountered) tags
  (https://support.bioconductor.org/p/110609/ ; Robert Bradley)

- (v 1.33.5) Check for corrupt index
  (https://github.com/Bioconductor/Rsamtools/issues/3 ; kjohnsen)


[Rsubread](https://bioconductor.org/packages/Rsubread)
--------

Changes in version 1.32.0:

- New function flattenGTF() that merges overlapping features into a
  single interval.

- New parameter for align() and subjunc(): sortReadsByCoordinates.

- New parameters for featureCounts(): readShiftType, readShiftSize and
  additionalAttributes.

- Specify strand protocol for each library individually in
  featureCounts().

- Much improved speed of align() and subjunc().

- align() and subjunc() return mapping statistics.

- Default setting of buildindex() is changed to building a one-block
  full index.

[RTN](https://bioconductor.org/packages/RTN)
---

Changes in version 2.6.0:

- Make available GSEA-2T and aREA-3T algorithms for single-sample
  analysis.

[RTNduals](https://bioconductor.org/packages/RTNduals)
--------

Changes in version 1.6.0:

- Improved workflow, integration with derivative packages (RTN /
  RTNsurvival).

[S4Vectors](https://bioconductor.org/packages/S4Vectors)
------

Changes in version 0.20.0

NEW FEATURES

- rbind() now supports DataFrame objects with the same column names
      but in different order, even when some of the column names are
      duplicated. How rbind() re-aligns the columns of the various objects
      to bind with those of the first object is consistent with what
      base:::rbind.data.frame() does.

- Add isSequence() low-level helper.

- Add 'nodup' argument to selectHits().

SIGNIFICANT USER-VISIBLE CHANGES

- The rownames of a DataFrame are no more required to be unique.

- Change 'use.names' default from FALSE to TRUE in mcols() getter.

- Coercion to DataFrame now **always** propagates the names.

- Rename low-level generic concatenateObjects() -> bindROWS().

- replaceROWS() now dispatches on 'x' and 'i' instead of 'x' only.

- Speedup row subsetting of DataFrame with many columns.

DEPRECATED AND DEFUNCT

- phead(), ptail(), and strsplitAsListOfIntegerVectors() are now defunct
      (after being deprecated in BioC 3.7).

BUG FIXES

- Fix window() on a DataFrame with data.frame columns.

- 2 fixes to "rbind" method for DataFrame objects:
      + It now properly handles DataFrame objects with duplicated colnames.
        Note that the new behavior is consistent with base::rbind.data.frame().
      + It now properly handles DataFrame objects with columns that are 1D
        arrays.

- Fix showAsCell() on nested data-frame-like objects.

- 2 fixes to "as.data.frame" method for DataFrame objects:
      + It now works if the DataFrame object contains nested data-frame-like
        objects or other complicated S4 objects (as long as these complicated
        objects in turn support as.data.frame()).
      + It now handles 'stringsAsFactors' argument properly. Originally
        reported here: https://github.com/Bioconductor/GenomicRanges/issues/18


[scater](https://bioconductor.org/packages/scater)
------

Changes in version 1.10.0:

- Fixes to all violin plots to ensure scatter matches up with violin
  outlines.

- Rectangle categorical/categorical plots collapse to mirrored bar
  plots when either factor contains only one level.

- Removed scater_gui(), downsampleCounts(), read10xResults(),
  normalizeExprs().

- Simplified plotRLE() to avoid the need for internal faceting.

- Added option for row subsetting in librarySizeFactors().

- Ensured calcAverage() with subset_row= behaves as if the matrix was
  subsetted prior to the function call.  Added support for
  parallelization.

- Ensured calculateCPM() with subset_row= behaves as if the matrix was
  subsetted prior to the function call.

- Added support for parallelization in nexprs().

- Added readSparseCounts() for creating a sparse matrix from a dense
  array on file.

- Added normalizeCounts() for easy division of matrix columns by the
  size factors.  Modified to throw error upon encountering negative, NA
  or zero size factors.

- Added preserve_zeroes= option to normalizeSCE() for preserving
  sparsity with non-unity pseudo-counts.

- Added runUMAP() and plotUMAP() to use the UMAP dimensionality
  reduction method.

- Added plotExplanatoryPCs() and getExplanatoryPCs() to correlate PCs
  with known factors.  Deprecated findImportantPCs().

- Added getVarianceExplained() to get the variance in gene expression
  explained by known factors.

- Removed runKallisto() and runSalmon().

- Switched readTxResults() to use tximport.  Switched
  readSalmonResults() and readKallistoResults() to use readTxResults().

- Removed obsolete fields in calculateQCMetrics().  Moved processing
  into C++ for a single-pass algorithm.  Supported parallelization
  across cells for QC computations.

- Added sumCountsAcrossFeatures() to sum counts across multiple
  redundant features.  Deprecated summariseExprsAcrossFeatures().

- All plotting functions can now access internal fields by using a
  character vector with NA as the first element.

- Returned threshold values in the attributes of the output from
  isOutlier().

- Deprecated the ticks in plotReducedDim().

[SCBN](https://bioconductor.org/packages/SCBN)
----

Changes in version 0.99.1:

- Modify description part in DESCRIPTION file.

- Change seq_len(length(x)) to seq_along(x) in sageTestNew.R function.

[scmeth](https://bioconductor.org/packages/scmeth)
------

Changes in version 1.1.7:

PKG FEATURES

- methylationDist function produces more interpretable mean methylation
  values

- In addition to the report file QC_Summary text file is generated.
  Also mbias table and downsample table are stored as text file

- Added "all" option in functions that utilizes subsampling such that
  given "all" option certain function would be applied to all the CpGs
  without any subsample or offset

[scPipe](https://bioconductor.org/packages/scPipe)
------

Changes in version 1.3.8:

- fix a minor bug in cell barcode demultiplexing

Changes in version 1.3.6:

- support multiple bam file for the same sample pooling as input for
  exon mapping

Changes in version 1.3.5:

- bug fix

Changes in version 1.3.4:

- put `distance_to_end` back.

Changes in version 1.3.2:

- Added gzipped output for `sc_trim_barcode()`, if output filename ends
  with `.gz` then gzipped output will be produced.

- Added `get_read_str()` function for getting common read structures.

Changes in version 1.3.1:

- updated the `sc_exon_mappping` function so it can accept multiple bam
  files now, together with a list of cell id or cell barcode.

[scran](https://bioconductor.org/packages/scran)
-----

Changes in version 1.10.0:

- Removed selectorPlot(), exploreData() functions in favour of iSEE.

- Fixed underflow problem in mnnCorrect() when dealing with the
  Gaussian kernel. Dropped the default sigma= in mnnCorrect() for
  better default performance.

- Supported parallelized block-wise processing in quickCluster().
  Deprecated max.size= in favour of max.cluster.size= in
  computeSumFactors(). Deprecated get.ranks= in favour of
  scaledColRanks().

- Added max.cluster.size= argument to computeSumFactors(). Supported
  parallelized cluster-wise processing. Increased all pool sizes to
  avoid rare failures if number of cells is a multiple of 5. Minor
  improvement to how mean filtering is done for rescaling across
  clusters in computeSumFactors(). Throw errors upon min.mean=NULL,
  which used to be valid. Switched positive=TRUE behaviour to use
  cleanSizeFactors().

- Added simpleSumFactors() as a simplified alternative to
  quickCluster() and computeSumFactors().

- Added the scaledColRanks() function for computing scaled and centred
  column ranks.

- Supported parallelized gene-wise processing in trendVar() and
  decomposeVar(). Support direct use of a factor in design= for
  efficiency.

- Added doubletCluster() to detect clusters that consist of doublets of
  other clusters.

- Added doubletCells() to detect cells that are doublets of other cells
  via simulations.

- Deprecated rand.seed= in buildSNNGraph() in favour of explicit
  set.seed() call. Added type= argument for weighting edges based on
  the number of shared neighbors.

- Deprecated rand.seed= in buildKNNGraph().

- Added multiBlockNorm() function for spike-abundance-preserving
  normalization prior to multi-block variance modelling.

- Added multiBatchNorm() function for consistent downscaling across
  batches prior to batch correction.

- Added cleanSizeFactors() to coerce non-positive size factors to
  positive values based on number of detected genes.

- Added the fastMNN() function to provide a faster, more stable
  alternative for MNN correction.

- Added BPPARAM= option for parallelized execution in makeTechTrend().
  Added approx.npts= option for interpolation-based approximation for
  many cells.

- Added pairwiseTTests() for direct calculation of pairwise
  t-statistics between groups.

- Added pairwiseWilcox() for direct calculation of pairwise Wilcoxon
  rank sum tests between groups.

- Added combineMarkers() to consolidate arbitrary pairwise comparisons
  into a marker list.

- Bugfixes to uses of block=, lfc= and design= arguments in
  findMarkers(). Refactored to use pairwiseTTests() and
  combineMarkers() internally. Added BPPARAM= option for parallelized
  execution.

- Refactored overlapExprs() to sort by p-value based on
  pairwiseWilcox() and combineMarkers(). Removed design= argument as it
  is not compatible with p-value calculations.

- Bugfixes to the use of Stouffer's Z method in combineVar().

- Added combinePValues() as a centralized internal function to combine
  p-values.

[SDAMS](https://bioconductor.org/packages/SDAMS)
-----

Changes in version 1.1.2:

- update NEWS file.

Changes in version 1.1.1:

- allow adjustment of covariates;

- update 'SDA' with additional arguments passed to 'qvalue'.

[SeqArray](https://bioconductor.org/packages/SeqArray)
--------

Changes in version 1.21.1-1.21.7:

UTILITIES

- avoid duplicated meta-information lines in `seqVCF2GDS()` and
  `seqVCF_Header()`

- require >= R_v3.5.0, since reading from connections in text mode is
  buffered

- `seqDigest()` requires the digest package

- optimization in reading genotypes from a subset of samples (according
  to gdsfmt_1.17.5)

NEW FEATURES

- `seqSNP2GDS()` imports dosage GDS files

- `seqVCF_Header()` allows a BCF file as an input

- a new function `seqRecompress()`

- a new function `seqCheck()` for checking the data integrity of a
  SeqArray GDS file

- `seqGDS2SNP()` exports dosage GDS files

BUG FIXES

- `seqVCF2GDS()` and `seqVCF_Header()` are able to import site-only VCF
  files (i.e., VCF with no sample)

- fix `seqVCF2GDS()` and `seqBCF2GDS()` since reading from connections
  in text mode is buffered for R >= v3.5.0

Changes in version 1.20.1:

BUG FIXES

- `seqExport()` fails to export haploid data (e.g., Y chromosome)

- `seqVCF2GDS()` fails to convert INFO variables when Number="R"

[seqCAT](https://bioconductor.org/packages/seqCAT)
------

Changes in version 1.4.0:

FEATURES

- Add convenience functions for creating and reading multiple SNV
  profiles

- Add functionality for reading general COSMIC mutational data, not
  just cell line mutational data

FIXES

- Fix an issue when reading COSMIC data due to new GRanges
  functionality

MISCELLANEOUS: 

- Update the citation info with the now-published seqCAT-specific
  article

[sesame](https://bioconductor.org/packages/sesame)
------

Changes in version 1.0.0:

- First submission of SeSAMe package.


[signeR](https://bioconductor.org/packages/signeR)
------

Changes in version 1.7.1:

- optimization of genCountMatrixFromVcf function

[SIMAT](https://bioconductor.org/packages/SIMAT)
-----

Changes in version 1.13.1:

- fix compatibility issues with new ggplot package &#91;2018-08-02 Wed&#93;

[SIMD](https://bioconductor.org/packages/SIMD)
----

Changes in version 0.99.9:

- put return statement in code at the end.

- instead of writing at each line, you paste output lines in a
  single paste call and then use a single write call.

- use a package BiocManager instead of biocLite in the vignettes.

- remove a few eval=FALSE sections in vignette to ensure vignette to be
  evaluated.

Changes in version 0.99.8:

- Modify vignette output to BiocStyle .

- Add some information in README file.

Changes in version 0.99.6:

- Adjust the code format to follow the coding style.

- Generate NEWS file in .Rd form.

- Add some information in README.md, vignettes and DESCRIPTION file to
  make the package SIMD be more understandable.

Changes in version 0.99.5:

- Modify the 'TIMEOUT' happened in R biocheck.

Changes in version 0.99.4:

- Modify some problems in R biocheck.

Changes in version 0.99.3:

- Modify some problems in R biocheck.

- Recompile .Rd by roxygen2 package in R.

Changes in version 0.99.2:

- Modify some problems in R biocheck.

- Change the NEWS format.

Changes in version 0.99.1:

- Modify the problem in R check.

- To get the vignettes by R Markdown instead of by Sweave.

- Add NEWS.md in package.

[SIMLR](https://bioconductor.org/packages/SIMLR)
-----

Changes in version 1.7.3 (2018-10-13):

- Removed CIMLR implementation.

[simulatorZ](https://bioconductor.org/packages/simulatorZ)
----------

Changes in version 1.15.1:

- Updated getTrueModel, zmatrix, simTime to include additional
  covariates as predictors

- Simplified examples

[SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment)
--------------------

Changes in version 1.4.0:

- Allow ... arguments to be passed to rowData() and colData().

- Added weights() methods for getting/setting observational weights.

- Added reducedDimNames<- method to set the names of reduced dimension
  slots.

- Added withDimnames= argument to reducedDim() and reducedDims().

- Exported getters and setters for internal metadata fields.

- Added developer instructions for making use of internal metadata
  fields.

[singleCellTK](https://bioconductor.org/packages/singleCellTK)
------------

Changes in version 1.1.26 (2018-10-23):

- New UI design for the Differential Expression tab.

- New UI design for the Data Summary & Filtering tab.

- Support for additional assay modification including log transforming
  any assay and renaming assays.

- New function visPlot for creating scatterplots, boxplots, heatmaps,
  and barplots for custom gene sets.

- The Downsample tab now works on a generic counts matrix

- You can upload a SCtkExperiment object or a SingleCellExperiment
  object saved in an RDS file on the Upload tab.

- Differential Expression results can now be saved in the rowData of
  the object and loaded for later analysis.

- Improved ability to save a biomarker based on user options.

- The Differential Expression plot is not automatically created, for
  more user control with large datasets.

Changes in version 1.1.3:

- Improvements to plotting, change text size and hide labels in gsva
  plots.

- MAST violin and linear model plots are now more square when plotting
  less than 49 facets.

- Changed y axis label in plotBatchVariance to "Percent Explained
  Variation"

Changes in version 1.1.2:

- Ability to hide version number in the SCTK GUI.

Changes in version 1.1.1:

- Fixed a bug that would cause the diffex color bar to not display when
  special characters were in the annotation.

[slinky](https://bioconductor.org/packages/slinky)
------

Changes in version 0.99.0:

NEW FEATURES

- Initial submission.  Everything is new and shiny.

[SNPRelate](https://bioconductor.org/packages/SNPRelate)
---------

Changes in version 1.15.1-1.15.5:

- a new option 'useMatrix' to allow for the packed symmetric matrix
  using the Matrix package in `snpgdsIBDMoM()`, `snpgdsIBDKING()` and
  `snpgdsIBS()`

- fix a bug of missing sample and SNP IDs in the output of
  `snpgdsIndInb()`

- new option 'start.pos' in `snpgdsLDpruning()`

- new methods in `snpgdsIndInb()`: gcta1, gcta2, gcta3; progress
  information is shown during running the function

- `snpgdsCombineGeno()` supports dosages

[SparseSignatures](https://bioconductor.org/packages/SparseSignatures)
----------------

Changes in version 1.1.4:

- Move NMF to Depends section

Changes in version 1.1.3:

- Issue with the basis function solved

[statTarget](https://bioconductor.org/packages/statTarget)
----------

Changes in version 2.0:

NEW FEATURES

- New GUI o Mouse Hover for help information o .log file

- New Signal correction o QC-RFSC methods for metabolomics and
  proteomics data

- New feature slection o Random Forest and the Permutation based
  variable importance measures o new MDSplot for Random Forest o
  P-value based importance plot

- New data preprocessing o PQN/SUM/none normalization o center/none
  Scaling method

[strandCheckR](https://bioconductor.org/packages/strandCheckR)
------------

Changes in version 0.99.16:

- Submitted version to Bioconductor

[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
--------------------

Changes in version 1.12.0:

NEW FEATURES

- The package has a new vignette "Extending the SummarizedExperiment class"
      by Aaron Lun intended for developers. It documents in great details the
      process of implementing a SummarizedExperiment extension (a.k.a.
      subclass).

SIGNIFICANT USER-VISIBLE CHANGES

- rowData() gains use.names=TRUE argument; prior behavior was to
  use.names=FALSE. rowData() by default fails when rownames() contains
  NAs.
  
BUG FIXES

- Better error handling in SummarizedExperiment() constructor.
      SummarizedExperiment() now prints an informative error message when
      the supplied assays have insane rownames or colnames. This addresses
      https://github.com/Bioconductor/SummarizedExperiment/issues/7

[SWATH2stats](https://bioconductor.org/packages/SWATH2stats)
-----------

Changes in version 1.11.7:

UPDATE

- avoid sorting when adding gene symbols in add_gene_symbol

Changes in version 1.11.6:

BUG FIXES

- fix manual page and function in convert_protein_ids to copy
  non-converted IDs

Changes in version 1.11.5:

UDPATE

- move from bioclite to BiocManager

Changes in version 1.11.4:

BUG FIXES

- fix manual page to convert_protein_ids

Changes in version 1.11.3:

BUG FIXES

- updates to convert_protein_ids, see that convert4aLFQ outputs
  character vectors

Changes in version 1.11.2:

NEW FEATURES

- add functions convert_protein_ids, load_mart, add_genesymbol

Changes in version 1.11.1:

BUG FIXES

- remove links from manual page SWATH2stats-package

Changes in version 1.11.0:

NEW FEATURES

- SWATH2stats in BioC 3.8 development release

[synapter](https://bioconductor.org/packages/synapter)
--------

Changes in version 2.5.2:

- Use `BiocManager::install` &#91;2018-07-16&#93;. # Synapter 2.3

[TargetSearch](https://bioconductor.org/packages/TargetSearch)
------------

Changes in version 1.38.0:

NEW FEATURES

- New function `checkRimLim` to vizualise a retention index markers
  before the actual time correction. It can be useful to fix the search
  limits.

- Peak detection method (NetCDFPeakFinding) has the option to use a
  gaussian smoothing in addition to usual moving average.

- Detects if CDF files are not found during sample description import.
  In addition, search for column names matching a pattern if the
  expected names are not found.

- Add support for a custom CDF file for faster data retrieval.

SIGNIFICANT USER-VISIBLE CHANGES

- The parameter `massRange` (m/z mass range) which used to be needed in
  some functions is deprected. It was used mostly as a hint and usually
  detected automatically. If it is passed, there would be no effect.

BUG FIXES

- Big refactor of C code to eliminate duplicated code and to separate
  what is R-C code (ie, SEXP structs) out of the C code that actually
  does something.

- General R code refactoring and housekeeping.

[TCC](https://bioconductor.org/packages/TCC)
---

Changes in version 1.21.3:

- fixed bug (order of outputs from devel version of baySeq cannot be
  sorted) in '.testByBayseq'.

- disabled 'samseq' option, since samr package has been removed from
  CRAN.

- change design matrix used in DESeq2

[TFEA.ChIP](https://bioconductor.org/packages/TFEA.ChIP)
---------

Changes in version 1.1.1:

New Features: 

- The function GeneID2entrez now suports translation from mouse ENSEMBL
  gene IDs and MGI symbols to mouse Entrez gene IDs

- 76 new ChIP-Seq experiments added to the database

[tofsims](https://bioconductor.org/packages/tofsims)
-------

Changes in version 099.1:

SIGNIFICANT USER-VISIBLE CHANGES

- changed function behvaiour in the whole package from call-by-ref to
  call-by value. Adjusted accordingly all examples and the vignette.

INTERNALS

- depends now on ProtGenerics from which it uses 'mz'

- exchanged various print() with message()

[topdownr](https://bioconductor.org/packages/topdownr)
--------

Changes in version 1.3.6:

- Add Pavel's and Ole's ORCID to DESCRIPTION &#91;2018-10-23&#93;.

Changes in version 1.3.5:

- Fix format of roxygen links to foreign packages to avoid link warning
  in `R CMD check` &#91;2018-10-10&#93;.

Changes in version 1.3.4:

- Add inst/CITATION file &#91;2018-09-26&#93;.

Changes in version 1.3.3:

- Revert commit c6e8dfd: "Adapt to `MSnbase 2.7.2` with internal
  fragments; see #82 &#91;2018-06-03&#93;."

Changes in version 1.3.2:

- Use `BiocManager::install` &#91;2018-07-16&#93;.

Changes in version 1.3.1:

- Adapt to `MSnbase 2.7.2` with internal fragments; see #82
  &#91;2018-06-03&#93;.

- Fix `FragmentViews` start/end/width and labels for internal fragments
  &#91;2018-06-03&#93;.

- Fix `as(tds, "MSnSet")` unit test &#91;2018-07-06&#93;.

- Use `elementMetadata(..., use.names=FALSE)` in
  `combine,FragmentViews,FragmentViews-method` to avoid duplicated
  rownames in elementMetadata slot &#91;2018-07-06&#93;.

Changes in version 1.3.0:

- New version for Bioc 3.8 (devel) # topdownr 1.2

[trackViewer](https://bioconductor.org/packages/trackViewer)
-----------

Changes in version 1.17.8:

- fix the bug in figure captions of vignette.

Changes in version 1.17.7:

- allow ylim of tracks be not fixed from 0.

Changes in version 1.17.6:

- fix a issue of circle in lollipop plot.

Changes in version 1.17.4:

- add flag type in lollipop plot.

Changes in version 1.17.3:

- adjust plot position for dandelion.plot.

Changes in version 1.17.2:

- add Yscales for dandelion.plot.

Changes in version 1.17.1:

- add smooth function for tracks.

[transcriptogramer](https://bioconductor.org/packages/transcriptogramer)
-----------------

Changes in version 1.3.6:

- Slot genesInTerm added to the Transcriptogram class.

- Argument boundaryConditions from differentiallyExpressed(): default
  changed from FALSE to TRUE.

- Argument onlyGenesInDE from clusterVisualization(): default changed
  from TRUE to FALSE.

- Argument onlyGenesInDE from clusterEnrichment(): default changed from
  TRUE to FALSE.

- Argument colors added to the differentiallyExpressed() method.

- Argument colors added to the clusterVisualization() method.

- Argument colors added to the enrichmentPlot() method.

- Argument alpha added to the enrichmentPlot() method.

Changes in version 1.3.4:

- Argument universe from clusterEnrichment(): default changed from "all
  the proteins present in the transcriptogramS2 slot" to "all the
  proteins present in the ordering slot".

Changes in version 1.3.3:

- Changes on the clusterEnrichment() method return.

- Slots Protein2GO, and Terms added to the Transcriptogram class.

- New methods: enrichmentPlot() and Terms().

Changes in version 1.3.1:

- Argument boundaryConditions added to the differentiallyExpressed()
  method.

- Slots Protein2Symbol, clusters, and pbc added to the Transcriptogram
  class.

- Argument onlyGenesInDE added to the clusterVisualization() method.

- Argument onlyGenesInDE added to the clusterEnrichment() method.

[transite](https://bioconductor.org/packages/transite)
--------

Changes in version 0.99.0:

- initial release

[tRNAscanImport](https://bioconductor.org/packages/tRNAscanImport)
--------------

Changes in version 1.1.9 (2018-10-24):

- moved some functionality to tRNA package

- added dependency for tRNA package

Changes in version 1.1.1 (2018-08-16):

- bugfix for recognizing pseudogene annotation

[tximeta](https://bioconductor.org/packages/tximeta)
-------

Changes in version 0.0.16:

- Added examples to all man pages.

[tximport](https://bioconductor.org/packages/tximport)
--------

Changes in version 1.9.11:

- Exporting simple internal function makeCountsFromAbundance().

Changes in version 1.9.10:

- Added 'infRepStat' argument which offers re-compution of counts and
  abundances using a function applied to the inferential replicates,
  e.g. matrixStats::rowMedian for using the median of posterior samples
  as the point estimate provided in "counts" and "abundance". If
  'countsFromAbundance' is specified, this will compute counts a second
  time from the re-computed abundances.

Changes in version 1.9.9:

- Adding support for gene-level summarization of inferential
  replicates. This takes place by perform row summarization on the
  inferential replicate (counts) in the same manner as the original
  counts (and optionally computing the variance).

Changes in version 1.9.6:

- Added new countsFromAbundance method: "dtuScaledTPM". This is
  designed for DTU analysis and to be used with txOut=TRUE. It provides
  counts that are scaled, with a gene, by the median transcript length
  among isoforms, then later by the sample's sequencing depth, as in
  the other two methods. The transcript lengths are calculated by first
  taking the average across samples. With this new method, all the
  abundances within a gene across all samples are scaled up by the same
  length, preserving isoform proportions calculated from the counts.

Changes in version 1.9.4:

- Made a change to summarizeToGene() that will now provide different
  output with a warning to alert the user. The case is: if tximport()
  is run with countsFromAbundance="scaledTPM" or "lengthScaledTPM" and
  txOut=TRUE, followed by summarizeToGene() with
  countsFromAbundance="no". This is a problematic series of calls, and
  previously it was ignoring the fact that the incoming counts are not
  original counts. Now, summarizeToGene() will throw a warning and
  override countsFromAbundance="no" to instead set it to the value that
  was used when tximport was originally run, either "scaledTPM" or
  "lengthScaledTPM".

Changes in version 1.9.1:

- Fixed edgeR example code in vignette to use scaleOffset after
  recommendation from Aaron Lun (2018-05-25).

[Ularcirc](https://bioconductor.org/packages/Ularcirc)
--------

Changes in version 0.99.0:

- Pre-Bioconductor submission

[Uniquorn](https://bioconductor.org/packages/Uniquorn)
--------

Changes in version 2.1.4 (2018-10-10):

Bioconductor compliance

- Minor bugfixes regarding the removal of CCLs

Changes in version 2.1.3 (2018-10-03):

Bioconductor compliance

- Minor bugfixes regarding the removal of CCLs and reference libraries

[variancePartition](https://bioconductor.org/packages/variancePartition)
-----------------

Changes in version 1.11.13:

- Fix multithreading issue

Changes in version 1.11.11:

- dream can handle multiple contrasts at the same time

Changes in version 1.11.10:

- fix typos in dream vignette

Changes in version 1.11.8:

- Check and stop() if response variable has variance of 0 - in dream(),
  fitExtractVarPartModel(), and fitVarPartModel()

- add standardized_t_stat() implicitly in eBayes() using MArrayLM2
  class - this transforms moderated t-statistics to have same degrees
  of freedom

Changes in version 1.11.7:

- Simplify object return by dream to be more more similar to lmFit -
  now returns MArrayLM instead of MArrayLMM_lmer

- if a fixed effects formula is specified (i.e. not random terms) -
  dream call lmFit in the backend - getContrast() works seamlessly

- dream() now returns gene annotation if passed to function

Changes in version 1.11.6:

- add error checing for L in dream

- fix typoes in dream vignette

- fix typoes in theory_practice_random_effects.Rnw

Changes in version 1.11.5:

- Add dream function for differential expression for repeated measures
  with a linear mixed model

Changes in version 1.11.2:

- Add warnings to canCorPairs for colinear terms

Changes in version 1.11.1:

- Add vignette: theory_practice_random_effects.Rnw

[VariantAnnotation](https://bioconductor.org/packages/VariantAnnotation)
-----------------

Changes in version 1.28.0:

NEW FEATURES

- Update package to support VCF format version 4.3 - SAMPLE field lines
  can now have key 'SAMPLE' or 'META'.  To avoid a name clash, the
  existing 'META' DataFrame has been split by row into separate
  DataFrames. The 'meta(VCFHeader)' getter now returns one DataFrame
  per unique key in the header. - PEDIGREE header line now begins with
  'ID'

- Add vcfFields method for character, VCFHeader, VcfFile and VCF to
  return all available vcf fields in CharacterList().

- Add support for single breakend notation (thanks d-cameron)

BUG FIXES

- .formatInfo() now return a column with all 'NA' for a missing value
  instead of dropping the column.

[xcms](https://bioconductor.org/packages/xcms)
----

Changes in version 3.3.6:

- Add type = "polygon" to highlightChromPeaks allowing to fill the
  actual signal area of identified chromatographic peaks.

Changes in version 3.3.5:

- Performance enhancement of the chromPeakSpectra and featureSpectra
  functions.

Changes in version 3.3.4:

- Add featureChromatograms to extract ion chromatograms for each
  feature.

- Add hasFilledChromPeaks function.

- Add argument skipFilled to the featureSummary function.

Changes in version 3.3.3:

- Add chromPeakSpectra and featureSpectra functions to extract MS2
  spectra for chromatographic peaks and features, respectively (issue
  #321).

- Fix profMat to handle also data files with empty spectra (issue
  #312).

- Add argument ylim to plotAdjustedRtime (issue #314).

- Add imputeRowMin and imputeRowMinRand, two simple missing value
  imputation helper functions.

- Fix additional problem mentioned in issue #301 with obiwarp retention
  time correction if some spectra have m/z values of `NA`.

- Fix issue #300 avoiding chromatographic peaks with rtmin > rtmax.

- Fixes for issues #291, #296.

- Add parameter 'missing' to diffreport allowing to replace NA with
  arbitrary numbers.

- Add exportMetaboAnalyst function to export the feature matrix in
  MetaboAnalyst format.

- Add parameter missing to featureValues allowing to specify how to
  handle/ report missing values.

- The chromPeaks matrix has now rownames to uniquely identify
  chromatographic peaks in an experiment. Chromatographic peak IDs
  start with "CP" followed by a number.

Changes in version 3.3.2:

- Add writeMSData method for XCMSnExp allowing to write mzML/mzXML
  files with adjusted retention times (issue #294).

- Fix profEIC call for single-scan-peak (pull request #287 from
  @trljcl).

- Fix centWave avoiding that the same peak is reported multiple times
  if fitgauss = TRUE is used (issue #284).

- featureSummary reports also RSD (relative standard deviations) of
  features across samples (issue #286).

- Add parameters fixedMz and fixedRt to FillChromPeaksParam that allow
  to increase the features' m/z and rt widths by a constant factor.

- Add option "sum" to featureValues' method parameter allowing to sum
  the intensities of peaks that are assigned to the same feature in a
  file/sample.

Changes in version 3.3.1:

- Add overlappingFeatures function to identify overlapping or close
  features.

- Add support for type = "apex_within" for featureDefinitions.

- Fix a bug in fillChromPeaks that would return the integrated signal
  being Inf.

- Fix for issue #267: error in fillChromPeaks when the retention time
  of the peaks are outside of the retention time range of certain
  files.

- New featureSummary function to calculate basic feature summaries
  (number of samples in which peaks were found etc).

- Parameter 'type' added to plotChromPeakDensity and 'whichPeaks' to
  highlightChromPeaks. Both parameters are passed to the 'type'
  argument of chromPeaks.

- Parameter 'type' in chromPeaks gets additional option "apex_within"
  to return chromatographic peaks that have their apex within the
  defined rt and/or m/z range.

- Add functions rla and rowRla to calculate RLA (relative log
  abundances).

- Add peaksWithMatchedFilter to perform peak detection in
  chromatographic (MRM/SRM) data (issues #277 and #278).

- Add peaksWithCentWave to perform centWave peak detection in
  chromatographic (MRM/SRM) data (issue #279).

- Add findChromPeaks,Chromatogram methods for CentWaveParam and
  MatchedFilterParam (issue #280).

[XINA](https://bioconductor.org/packages/XINA)
----

Changes in version 1.0.0:

NEW FEATURES

- No changes

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

[xps](https://bioconductor.org/packages/xps)
---

Changes in version 1.41:

- update NEWS file

[zinbwave](https://bioconductor.org/packages/zinbwave)
--------

Changes in version 1.3.1 (2018-05-09):

- New `zinbsurf` function implements approximate method for large
  matrices.

- New option `which_genes` in `zinbwave` to specify which genes to use
  to compute `W`.


NEWS from new and existing Data Experiment Packages
===================================

[celarefData](https://bioconductor.org/packages/celarefData)
-----------

Changes in version 0.99.0 (2018-08-05):

- Submitted to Bioconductor

[CopyNeutralIMA](https://bioconductor.org/packages/CopyNeutralIMA)
--------------

Changes in version 1.0:

- Package release

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

[DuoClustering2018](https://bioconductor.org/packages/DuoClustering2018)
-----------------

Changes in version 0.1.0:

- Initial version of the package

[FlowSorted.Blood.EPIC](https://bioconductor.org/packages/FlowSorted.Blood.EPIC)
---------------------

Changes in version 0.99.37 (2018-10-25):

- Made the following significant changes o Initial release to
  Bioconductor o Updated manuscript reference

Changes in version 0.99.36 (2018-03-24):

- Submitted to Bioconductor

- Added NEWS file

[GIGSEAdata](https://bioconductor.org/packages/GIGSEAdata)
----------

Changes in version 0.99.8 (2018-06-30):

- Changed: GO.rda to org.Hs.eg.GO.rda

- Added: vignettes

Changes in version 0.99.1 (2018-06-04):

- Added: GO.rda and Fantom5.TF.rda

[mcsurvdata](https://bioconductor.org/packages/mcsurvdata)
----------

Changes in version 1.0.1:

- Package submission


[pRolocdata](https://bioconductor.org/packages/pRolocdata)
----------

Changes in version 1.19.4:

- New yeast spatial proteome dataset <2018-08-14 Tue>

Changes in version 1.19.3:

- New synechocystis spatial proteome <2018-07-26 Thu>

Changes in version 1.19.2:

- Fix typo in beltran2016 man page <2018-07-24 Tue>

- Added LOPIT-DC and hyperLOPIT U2OS data <2018-07-25 Wed>

Changes in version 1.19.1:

- Adding hyperLOPIT TAGM results <2018-05-21 Mon>

Changes in version 1.19.0:

- New Bioconductor devel version

[PtH2O2lipids](https://bioconductor.org/packages/PtH2O2lipids)
------------

Changes in version 2016-04-21:

- Initial release for Bioconductor

[qPLEXdata](https://bioconductor.org/packages/qPLEXdata)
---------

Changes in version 0.99.0:

- Initial commit with data from the qPLEX-RIME and Full proteome TMT.


[RforProteomics](https://bioconductor.org/packages/RforProteomics)
--------------

Changes in version 1.19.3:

- Use BiocManager

- Fix typo in rmd

Changes in version 1.19.0:

- New Bioconductor devel version

[sesameData](https://bioconductor.org/packages/sesameData)
----------

Changes in version 1.0.0:

- First submission of sesameData package.

[TabulaMurisData](https://bioconductor.org/packages/TabulaMurisData)
---------------

Changes in version 0.99.0:

- Initial version of the package

[tcgaWGBSData.hg19](https://bioconductor.org/packages/tcgaWGBSData.hg19)
-----------------

Changes in version 0.99.3:

PKG FEATURES

- This is a data package containing Whole Genome Bisulfite Sequencing
  (WGBS) data from TCGA.


NEWS from new and existing Workflows
===================================

[maEndToEnd](https://bioconductor.org/packages/maEndToEnd)
----------

Changes in version 1.99.6:

- some minor fixes to RLE section

- minor text edits as suggested by James McDonald at F1000

- modified to code chunk that checks for latex compilation

Changes in version 1.99.2:

- changed installation instructions to BiocManager

Changes in version 0.99.0:

- new submission to Bioc in June 18

[recountWorkflow](https://bioconductor.org/packages/recountWorkflow)
---------------

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Use BiocManager


Deprecated and Defunct Packages
===============================

Seven software packages were removed from this release (after being deprecated
in Bioc 3.7): ontoCat, spliceR, OperaMate, DASC, PAnnBuilder, phenoDist, BrowserVizDemo.

Thirteen software are deprecated in this release and will be removed in Bioc 3.9:
BiocInstaller, GoogleGenomics, IrisSpatialFeatures, facopy, gaucho, nudge,
RamiGO, mQTL.NMR, cytofkit, pbcmc, GeneSelector, ampliQueso, prot2D.

Three experimental data packages were removed in this release (after being
deprecated in BioC 3.7):  RnaSeqTutorial, cheung2010, MEALData.

Two experimental data packages are deprecated in this release and will be
removed in Bioc 3.9: iontreeData, MSBdata.

