October 18, 2016

Bioconductors:

We are pleased to announce Bioconductor 3.4, consisting of 1295
software packages, 309 experiment data packages, and 933
up-to-date annotation packages.

There are 100 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.4 is compatible with R 3.3,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release will include an updated Bioconductor [Amazon Machine Image][1]
and [Docker containers][2].

Visit [http://bioconductor.org][3]
for details and downloads.

[1]: http://bioconductor.org/help/bioconductor-cloud-ami/
[2]: http://bioconductor.org/help/docker/
[3]: http://bioconductor.org

Contents
--------

* [Getting Started with Bioconductor 3.4](#getting-started-with-bioconductor-34)
* [New Software Packages](#new-software-packages)
* [NEWS from new and existing packages](#news-from-new-and-existing-packages)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.4
======================================

To update to or install Bioconductor 3.4:

1. Install R 3.3 (>= 3.3.1 recommended).  Bioconductor 3.4 has been designed expressly for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/) .

New Software Packages
=====================

There are 100 new software packages in this release of Bioconductor.

- [alpine](https://bioconductor.org/packages/alpine) Fragment
  sequence bias modeling and correction for RNA-seq transcript
  abundance estimation.

- [AMOUNTAIN](https://bioconductor.org/packages/AMOUNTAIN) A pure
  data-driven gene network, weighted gene co-expression network
  (WGCN) could be constructed only from expression profile. Different
  layers in such networks may represent different time points,
  multiple conditions or various species. AMOUNTAIN aims to search
  active modules in multi-layer WGCN using a continuous optimization
  approach.

- [anamiR](https://bioconductor.org/packages/anamiR) This package is
  intended to identify potential interactions of miRNA-target gene
  interactions from miRNA and mRNA expression data. It contains
  functions for statistical test, databases of miRNA-target gene
  interaction and functional analysis.

- [Anaquin](https://bioconductor.org/packages/Anaquin) The project is
  intended to support the use of sequins (synthetic sequencing
  spike-in controls) owned and made available by the Garvan Institute
  of Medical Research. The goal is to provide a standard open source
  library for quantitative analysis, modelling and visualization of
  spike-in controls.

- [annotatr](https://bioconductor.org/packages/annotatr) Given a set
  of genomic sites/regions (e.g. ChIP-seq peaks, CpGs, differentially
  methylated CpGs or regions, SNPs, etc.) it is often of interest to
  investigate the intersecting genomic annotations. Such annotations
  include those relating to gene models (promoters, 5'UTRs, exons,
  introns, and 3'UTRs), CpGs (CpG islands, CpG shores, CpG shelves),
  or regulatory sequences such as enhancers. The annotatr package
  provides an easy way to summarize and visualize the intersection of
  genomic sites/regions with genomic annotations.

- [ASAFE](https://bioconductor.org/packages/ASAFE) Given admixed
  individuals' bi-allelic SNP genotypes and ancestry pairs (where
  each ancestry can take one of three values) for multiple SNPs,
  perform an EM algorithm to deal with the fact that SNP genotypes
  are unphased with respect to ancestry pairs, in order to estimate
  ancestry-specific allele frequencies for all SNPs.

- [ASpli](https://bioconductor.org/packages/ASpli) Integrative
  pipeline for the analyisis of alternative splicing using RNAseq.

- [BaalChIP](https://bioconductor.org/packages/BaalChIP) The package
  offers functions to process multiple ChIP-seq BAM files and detect
  allele-specific events. Computes allele counts at individual
  variants (SNPs/SNVs), implements extensive QC steps to remove
  problematic variants, and utilizes a bayesian framework to identify
  statistically significant allele- specific events. BaalChIP is able
  to account for copy number differences between the two alleles, a
  known phenotypical feature of cancer samples.

- [BayesKnockdown](https://bioconductor.org/packages/BayesKnockdown)
  A simple, fast Bayesian method for computing posterior
  probabilities for relationships between a single predictor variable
  and multiple potential outcome variables, incorporating prior
  probabilities of relationships. In the context of knockdown
  experiments, the predictor variable is the knocked-down gene, while
  the other genes are potential targets. Can also be used for
  differential expression/2-class data.

- [bigmelon](https://bioconductor.org/packages/bigmelon) Methods for
  working with Illumina arrays using gdsfmt.

- [bioCancer](https://bioconductor.org/packages/bioCancer) bioCancer
  is a Shiny App to visualize and analyse interactively Multi-Assays
  of Cancer Genomic Data.

- [BiocWorkflowTools](https://bioconductor.org/packages/BiocWorkflowTools)
  Provides functions to ease the transition between Rmarkdown and
  LaTeX documents when authoring a Bioconductor Workflow.

- [CancerInSilico](https://bioconductor.org/packages/CancerInSilico)
  The CancerInSilico package provides an R interface for running
  mathematical models of tumor progresson. This package has the
  underlying models implemented in C++ and the output and analysis
  features implemented in R.

- [CancerSubtypes](https://bioconductor.org/packages/CancerSubtypes)
  CancerSubtypes integrates the current common computational biology
  methods for cancer subtypes identification and provides a
  standardized framework for cancer subtype analysis based on the
  genomic datasets.

- [ccmap](https://bioconductor.org/packages/ccmap) Finds drugs and
  drug combinations that are predicted to reverse or mimic gene
  expression signatures. These drugs might reverse diseases or mimic
  healthy lifestyles.

- [CCPROMISE](https://bioconductor.org/packages/CCPROMISE) Perform
  Canonical correlation between two forms of high demensional genetic
  data, and associate the first compoent of each form of data with a
  specific biologically interesting pattern of associations with
  multiple endpoints. A probe level analysis is also implemented.

- [CellMapper](https://bioconductor.org/packages/CellMapper) Infers
  cell type-specific expression based on co-expression similarity
  with known cell type marker genes. Can make accurate predictions
  using publicly available expression data, even when a cell type has
  not been isolated before.

- [chromstaR](https://bioconductor.org/packages/chromstaR) This
  package implements functions for combinatorial and differential
  analysis of ChIP-seq data. It includes uni- and multivariate
  peak-calling, export to genome browser viewable files, and
  functions for enrichment analyses.

- [clusterExperiment](https://bioconductor.org/packages/clusterExperiment)
  This package provides functions for running and comparing many
  different clusterings of single-cell sequencing data.

- [covEB](https://bioconductor.org/packages/covEB) Using bayesian
  methods to estimate correlation matrices assuming that they can be
  written and estimated as block diagonal matrices. These block
  diagonal matrices are determined using shrinkage parameters that
  values below this parameter to zero.

- [covRNA](https://bioconductor.org/packages/covRNA) This package
  provides the analysis methods fourthcorner and RLQ analysis for
  large-scale transcriptomic data.

- [crisprseekplus](https://bioconductor.org/packages/crisprseekplus)
  Bioinformatics platform containing interface to work with
  offTargetAnalysis and compare2Sequences in the CRISPRseek package,
  and GUIDEseqAnalysis.

- [crossmeta](https://bioconductor.org/packages/crossmeta) Implements
  cross-platform and cross-species meta-analyses of Affymentrix,
  Illumina, and Agilent microarray data. This package automates
  common tasks such as downloading, normalizing, and annotating raw
  GEO data. A user interface makes it easy to select control and
  treatment samples for each contrast and study. This input is used
  for subsequent surrogate variable analysis (models unaccounted
  sources of variation) and differential expression analysis. Final
  meta-analysis of differential expression values can include genes
  measured in only a subset of studies.

- [ctsGE](https://bioconductor.org/packages/ctsGE) Methodology for
  supervised clustering of potentially many predictor variables, such
  as genes etc., in time series datasets Provides functions that help
  the user assigning genes to predefined set of model profiles.

- [CVE](https://bioconductor.org/packages/CVE) Shiny app for
  interactive variant prioritisation in precision cancer medicine.
  The input file for CVE is the output file of the recently released
  Oncotator Variant Annotation tool summarising variant-centric
  information from 14 different publicly available resources relevant
  for cancer researches. Interactive priortisation in CVE is based on
  known germline and cancer variants, DNA repair genes and functional
  prediction scores. An optional feature of CVE is the exploration of
  the tumour-specific pathway context that is facilitated using
  co-expression modules generated from publicly available
  transcriptome data. Finally druggability of prioritised variants is
  assessed using the Drug Gene Interaction Database (DGIdb).

- [CytoML](https://bioconductor.org/packages/CytoML) This package is
  designed to use GatingML2.0 as the standard format to exchange the
  gated data with other software platform.

- [DeepBlueR](https://bioconductor.org/packages/DeepBlueR) Accessing
  the DeepBlue Epigenetics Data Server through R.

- [DEsubs](https://bioconductor.org/packages/DEsubs) DEsubs is a
  network-based systems biology package that extracts
  disease-perturbed subpathways within a pathway network as recorded
  by RNA-seq experiments. It contains an extensive and customizable
  framework covering a broad range of operation modes at all stages
  of the subpathway analysis, enabling a case-specific approach. The
  operation modes refer to the pathway network construction and
  processing, the subpathway extraction, visualization and enrichment
  analysis with regard to various biological and pharmacological
  features. Its capabilities render it a tool-guide for both the
  modeler and experimentalist for the identification of more robust
  systems-level biomarkers for complex diseases.

- [Director](https://bioconductor.org/packages/Director) Director is
  an R package designed to streamline the visualization of molecular
  effects in regulatory cascades. It utilizes the R package htmltools
  and a modified Sankey plugin of the JavaScript library D3 to
  provide a fast and easy, browser-enabled solution to discovering
  potentially interesting downstream effects of regulatory and/or
  co-expressed molecules. The diagrams are robust, interactive, and
  packaged as highly-portable HTML files that eliminate the need for
  third-party software to view. This enables a straightforward
  approach for scientists to interpret the data produced, and
  bioinformatics developers an alternative means to present relevant
  data.

- [dSimer](https://bioconductor.org/packages/dSimer) dSimer is an R
  package which provides computation of nine methods for measuring
  disease-disease similarity, including a standard cosine similarity
  measure and eight function-based methods. The disease similarity
  matrix obtained from these nine methods can be visualized through
  heatmap and network. Biological data widely used in disease-disease
  associations study are also provided by dSimer.

- [eegc](https://bioconductor.org/packages/eegc) This package has
  been developed to evaluate cellular engineering processes for
  direct differentiation of stem cells or conversion
  (transdifferentiation) of somatic cells to primary cells based on
  high throughput gene expression data screened either by DNA
  microarray or RNA sequencing. The package takes gene expression
  profiles as inputs from three types of samples: (i) somatic or stem
  cells to be (trans)differentiated (input of the engineering
  process), (ii) induced cells to be evaluated (output of the
  engineering process) and (iii) target primary cells (reference for
  the output). The package performs differential gene expression
  analysis for each pair-wise sample comparison to identify and
  evaluate the transcriptional differences among the 3 types of
  samples (input, output, reference). The ideal goal is to have
  induced and primary reference cell showing overlapping profiles,
  both very different from the original cells.

- [esetVis](https://bioconductor.org/packages/esetVis) Utility
  functions for visualization of expressionSet (or
  SummarizedExperiment) Bioconductor object, including spectral map,
  tsne and linear discriminant analysis. Static plot via the ggplot2
  package or interactive via the ggvis or rbokeh packages are
  available.

- [ExperimentHub](https://bioconductor.org/packages/ExperimentHub)
  This package provides a client for the Bioconductor ExperimentHub
  web resource. ExperimentHub provides a central location where
  curated data from experiments, publications or training courses can
  be accessed. Each resource has associated metadata, tags and date
  of modification. The client creates and manages a local cache of
  files retrieved enabling quick and reproducible access.

- [ExperimentHubData](https://bioconductor.org/packages/ExperimentHubData)
  Functions to add metadata to ExperimentHub db and resource files to
  AWS S3 buckets.

- [fCCAC](https://bioconductor.org/packages/fCCAC) An application of
  functional canonical correlation analysis to assess covariance of
  nucleic acid sequencing datasets such as chromatin
  immunoprecipitation followed by deep sequencing (ChIP-seq).

- [fgsea](https://bioconductor.org/packages/fgsea) The package
  implements an algorithm for fast gene set enrichment analysis.
  Using the fast algorithm allows to make more permutations and get
  more fine grained p-values, which allows to use accurate stantard
  approaches to multiple hypothesis correction.

- [FitHiC](https://bioconductor.org/packages/FitHiC) Fit-Hi-C is a
  tool for assigning statistical confidence estimates to
  intra-chromosomal contact maps produced by genome-wide genome
  architecture assays such as Hi-C.

- [flowPloidy](https://bioconductor.org/packages/flowPloidy)
  Determine sample ploidy via flow cytometry histogram analysis.
  Reads Flow Cytometry Standard (FCS) files via the flowCore
  bioconductor package, and provides functions for determining the
  DNA ploidy of samples based on internal standards.

- [FunChIP](https://bioconductor.org/packages/FunChIP) Preprocessing
  and smoothing of ChIP-Seq peaks and efficient implementation of the
  k-mean alignment algorithm to classify them.

- [GAprediction](https://bioconductor.org/packages/GAprediction)
  [GAprediction] predicts gestational age using Illumina
  HumanMethylation450 CpG data.

- [gCrisprTools](https://bioconductor.org/packages/gCrisprTools) Set
  of tools for evaluating pooled high-throughput screening
  experiments, typically employing CRISPR/Cas9 or shRNA expression
  cassettes. Contains methods for interrogating library and cassette
  behavior within an experiment, identifying differentially abundant
  cassettes, aggregating signals to identify candidate targets for
  empirical validation, hypothesis testing, and comprehensive
  reporting.

- [GEM](https://bioconductor.org/packages/GEM) Tools for analyzing
  EWAS, methQTL and GxE genome widely.

- [geneAttribution](https://bioconductor.org/packages/geneAttribution)
  Identification of the most likely gene or genes through which
  variation at a given genomic locus in the human genome acts. The
  most basic functionality assumes that the closer gene is to the
  input locus, the more likely the gene is to be causative.
  Additionally, any empirical data that links genomic regions to
  genes (e.g. eQTL or genome conformation data) can be used if it is
  supplied in the UCSC .BED file format.

- [GeneGeneInteR](https://bioconductor.org/packages/GeneGeneInteR)
  The aim of this package is to propose several methods for testing
  gene-gene interaction in case-control association studies. Such a
  test can be done by aggregating SNP-SNP interaction tests performed
  at the SNP level (SSI) or by using gene-gene multidimensionnal
  methods (GGI) methods. The package also proposes tools for a
  graphic display of the results.

- [geneplast](https://bioconductor.org/packages/geneplast) Geneplast
  is designed for evolutionary and plasticity analysis based on
  orthologous groups distribution in a given species tree. It uses
  Shannon information theory and orthologs abundance to estimate the
  Evolutionary Plasticity Index. Additionally, it implements the
  Bridge algorithm to determine the evolutionary root of a given gene
  based on its orthologs distribution.

- [geneXtendeR](https://bioconductor.org/packages/geneXtendeR)
  geneXtendeR is designed to optimally annotate a histone
  modification ChIP-seq peak input file with functionally important
  genomic features (e.g., genes associated with peaks) based on
  optimization calculations.  geneXtendeR optimally extends the
  boundaries of every gene in a genome by some genomic distance (in
  DNA base pairs) for the purpose of flexibly incorporating
  cis-regulatory elements (CREs), such as enhancers and promoters, as
  well as downstream elements that are important to the function of
  the gene relative to an epigenetic histone modification ChIP-seq
  dataset. geneXtender computes optimal gene extensions tailored to
  the broadness of the specific epigenetic mark (e.g., H3K9me1,
  H3K27me3), as determined by a user-supplied ChIP-seq peak input
  file. As such, geneXtender maximizes the signal-to-noise ratio of
  locating genes closest to and directly under peaks. By performing a
  computational expansion of this nature, ChIP-seq reads that would
  initially not map strictly to a specific gene can now be optimally
  mapped to the regulatory regions of the gene, thereby implicating
  the gene as a potential candidate, and thereby making the ChIP-seq
  experiment more successful. Such an approach becomes particularly
  important when working with epigenetic histone modifications that
  have inherently broad peaks.

- [GOpro](https://bioconductor.org/packages/GOpro) Find the most
  characteristic gene ontology terms for groups of human genes. This
  package was created as a part of the thesis which was developed
  under the auspices of MI^2 Group (http://mi2.mini.pw.edu.pl/,
  https://github.com/geneticsMiNIng).

- [GRmetrics](https://bioconductor.org/packages/GRmetrics) Functions
  for calculating and visualizing growth-rate inhibition (GR)
  metrics.

- [HelloRanges](https://bioconductor.org/packages/HelloRanges)
  Translates bedtools command-line invocations to R code calling
  functions from the Bioconductor *Ranges infrastructure. This is
  intended to educate novice Bioconductor users and to compare the
  syntax and semantics of the two frameworks.

- [ImpulseDE](https://bioconductor.org/packages/ImpulseDE) ImpulseDE
  is suited to capture single impulse-like patterns in high
  throughput time series datasets. By fitting a representative
  impulse model to each gene, it reports differentially expressed
  genes whether across time points in a single experiment or between
  two time courses from two experiments. To optimize the running
  time, the code makes use of clustering steps and multi-threading.

- [IPO](https://bioconductor.org/packages/IPO) The outcome of XCMS
  data processing strongly depends on the parameter settings. IPO
  (`Isotopologue Parameter Optimization`) is a parameter optimization
  tool that is applicable for different kinds of samples and liquid
  chromatography coupled to high resolution mass spectrometry
  devices, fast and free of labeling steps. IPO uses natural, stable
  13C isotopes to calculate a peak picking score. Retention time
  correction is optimized by minimizing the relative retention time
  differences within features and grouping parameters are optimized
  by maximizing the number of features showing exactly one peak from
  each injection of a pooled sample. The different parameter settings
  are achieved by design of experiment. The resulting scores are
  evaluated using response surface models.

- [KEGGlincs](https://bioconductor.org/packages/KEGGlincs) See what
  is going on 'under the hood' of KEGG pathways by explicitly
  re-creating the pathway maps from information obtained from KGML
  files.

- [LINC](https://bioconductor.org/packages/LINC) This package
  provides methods to compute co-expression networks of lincRNAs and
  protein-coding genes. Biological terms associated with the sets of
  protein-coding genes predict the biological contexts of lincRNAs
  according to the 'Guilty by Association' approach.

- [LOBSTAHS](https://bioconductor.org/packages/LOBSTAHS) LOBSTAHS is
  a multifunction package for screening, annotation, and putative
  identification of mass spectral features in large, HPLC-MS lipid
  datasets. In silico data for a wide range of lipids, oxidized
  lipids, and oxylipins can be generated from user-supplied
  structural criteria with a database generation function. LOBSTAHS
  then applies these databases to assign putative compound identities
  to features in any high-mass accuracy dataset that has been
  processed using xcms and CAMERA. Users can then apply a series of
  orthogonal screening criteria based on adduct ion formation
  patterns, chromatographic retention time, and other properties, to
  evaluate and assign confidence scores to this list of preliminary
  assignments. During the screening routine, LOBSTAHS rejects
  assignments that do not meet the specified criteria, identifies
  potential isomers and isobars, and assigns a variety of annotation
  codes to assist the user in evaluating the accuracy of each
  assignment.

- [M3Drop](https://bioconductor.org/packages/M3Drop) This package
  fits a Michaelis-Menten model to the pattern of dropouts in
  single-cell RNASeq data. This model is used as a null to identify
  significantly variable (i.e. differentially expressed) genes for
  use in downstream analysis, such as clustering cells.

- [MADSEQ](https://bioconductor.org/packages/MADSEQ) The MADSEQ
  package provides a group of hierarchical Bayeisan models for the
  detection of mosaic aneuploidy, the inference of the type of
  aneuploidy and also for the quantification of the fraction of
  aneuploid cells in the sample.

- [maftools](https://bioconductor.org/packages/maftools) Analyze and
  visualize Mutation Annotation Format (MAF) files from large scale
  sequencing studies. This package provides various functions to
  perform most commonly used analyses in cancer genomics and to
  create feature rich customizable visualzations with minimal effort.

- [MAST](https://bioconductor.org/packages/MAST) Methods and models
  for handling zero-inflated single cell assay data.

- [matter](https://bioconductor.org/packages/matter) Memory-efficient
  reading, writing, and manipulation of structured binary data on
  disk as vectors, matrices, and arrays. This package is designed to
  be used as a back-end for Cardinal for working with high-resolution
  mass spectrometry imaging data.

- [meshes](https://bioconductor.org/packages/meshes) MeSH (Medical
  Subject Headings) is the NLM controlled vocabulary used to manually
  index articles for MEDLINE/PubMed. MeSH terms were associated by
  Entrez Gene ID by three methods, gendoo, gene2pubmed and RBBH. This
  association is fundamental for enrichment and semantic analyses.
  meshes supports enrichment analysis (over-representation and gene
  set enrichment analysis) of gene list or whole expression profile.
  The semantic comparisons of MeSH terms provide quantitative ways to
  compute similarities between genes and gene groups. meshes
  implemented five methods proposed by Resnik, Schlicker, Jiang, Lin
  and Wang respectively and supports more than 70 species.

- [MetaboSignal](https://bioconductor.org/packages/MetaboSignal)
  MetaboSignal is an R package that allows merging, analyzing and
  customizing metabolic and signaling KEGG pathways. It is a
  network-based approach designed to explore the topological
  relationship between genes (signaling- or enzymatic-genes) and
  metabolites, representing a powerful tool to investigate the
  genetic landscape and regulatory networks of metabolic phenotypes.

- [MetCirc](https://bioconductor.org/packages/MetCirc) MetCirc
  comprises a workflow to interactively explore metabolomics data:
  create MSP, bin m/z values, calculate similarity between precursors
  and visualise similarities.

- [methylKit](https://bioconductor.org/packages/methylKit) methylKit
  is an R package for DNA methylation analysis and annotation from
  high-throughput bisulfite sequencing. The package is designed to
  deal with sequencing data from RRBS and its variants, but also
  target-capture methods and whole genome bisulfite sequencing. It
  also has functions to analyze base-pair resolution 5hmC data from
  experimental protocols such as oxBS-Seq and TAB-Seq. Perl is needed
  to read SAM files only.

- [MGFR](https://bioconductor.org/packages/MGFR) The package is
  designed to detect marker genes from RNA-seq data.

- [MODA](https://bioconductor.org/packages/MODA) MODA can be used to
  estimate and construct condition-specific gene co-expression
  networks, and identify differentially expressed subnetworks as
  conserved or condition specific modules which are potentially
  associated with relevant biological processes.

- [MoonlightR](https://bioconductor.org/packages/MoonlightR)
  Motivation: The understanding of cancer mechanism requires the
  identification of genes playing a role in the development of the
  pathology and the characterization of their role (notably oncogenes
  and tumor suppressors). Results: We present an R/bioconductor
  package called MoonlightR which returns a list of candidate driver
  genes for specific cancer types on the basis of TCGA expression
  data. The method first infers gene regulatory networks and then
  carries out a functional enrichment analysis (FEA) (implementing an
  upstream regulator analysis, URA) to score the importance of
  well-known biological processes with respect to the studied cancer
  type. Eventually, by means of random forests, MoonlightR predicts
  two specific roles for the candidate driver genes: i) tumor
  suppressor genes (TSGs) and ii) oncogenes (OCGs). As a consequence,
  this methodology does not only identify genes playing a dual role
  (e.g. TSG in one cancer type and OCG in another) but also helps in
  elucidating the biological processes underlying their specific
  roles. In particular, MoonlightR can be used to discover OCGs and
  TSGs in the same cancer type. This may help in answering the
  question whether some genes change role between early stages (I,
  II) and late stages (III, IV) in breast cancer. In the future, this
  analysis could be useful to determine the causes of different
  resistances to chemotherapeutic treatments.

- [msPurity](https://bioconductor.org/packages/msPurity) Assess the
  contribution of the targeted precursor in fragmentation acquired or
  anticipated isolation windows using a metric called "precursor
  purity". Also provides simple processing steps (averaging,
  filtering, blank subtraction, etc) for DI-MS data. Works for both
  LC-MS(/MS) and DI-MS(/MS) data.

- [MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment)
  Develop an integrative environment where multiple assays are
  managed and preprocessed for genomic data analysis.

- [MutationalPatterns](https://bioconductor.org/packages/MutationalPatterns)
  An extensive toolset for the characterization and visualization of
  a wide range of mutational patterns in base substitution data.

- [netprioR](https://bioconductor.org/packages/netprioR) A model for
  semi-supervised prioritisation of genes integrating network data,
  phenotypes and additional prior knowledge about TP and TN gene
  labels from the literature or experts.

- [normr](https://bioconductor.org/packages/normr) Robust
  normalization and difference calling procedures for ChIP-seq and
  alike data. Read counts are modeled jointly as a binomial mixture
  model with a user-specified number of components. A fitted
  background estimate accounts for the effect of enrichment in
  certain regions and, therefore, represents an appropriate null
  hypothesis. This robust background is used to identify
  significantly enriched or depleted regions.

- [PathoStat](https://bioconductor.org/packages/PathoStat) The
  purpose of this package is to perform Statistical Microbiome
  Analysis on metagenomics results from sequencing data samples. In
  particular, it supports analyses on the PathoScope generated report
  files. PathoStat provides various functionalities including
  Relative Abundance charts, Diversity estimates and plots, tests of
  Differential Abundance, Time Series visualization, and Core OTU
  analysis.

- [PharmacoGx](https://bioconductor.org/packages/PharmacoGx) Contains
  a set of functions to perform large-scale analysis of
  pharmacogenomic data.

- [philr](https://bioconductor.org/packages/philr) PhILR is short for
  Phylogenetic Isometric Log-Ratio Transform. This package provides
  functions for the analysis of compositional data (e.g., data
  representing proportions of different variables/parts).
  Specifically this package allows analysis of compositional data
  where the parts can be related through a phylogenetic tree (as is
  common in microbiota survey data) and makes available the Isometric
  Log Ratio transform built from the phylogenetic tree and utilizing
  a weighted reference measure.

- [Pi](https://bioconductor.org/packages/Pi) Priority index or Pi is
  developed as a genomic-led target prioritisation system, with the
  focus on leveraging human genetic data to prioritise potential drug
  targets at the gene, pathway and network level. The long term goal
  is to use such information to enhance early-stage target
  validation. Based on evidence of disease association from
  genome-wide association studies (GWAS), this prioritisation system
  is able to generate evidence to support identification of the
  specific modulated genes (seed genes) that are responsible for the
  genetic association signal by utilising knowledge of linkage
  disequilibrium (co-inherited genetic variants), distance of
  associated variants from the gene, and evidence of independent
  genetic association with gene expression in disease-relevant
  tissues, cell types and states. Seed genes are scored in an
  integrative way, quantifying the genetic influence. Scored seed
  genes are subsequently used as baits to rank seed genes plus
  additional (non-seed) genes; this is achieved by iteratively
  exploring the global connectivity of a gene interaction network.
  Genes with the highest priority are further used to
  identify/prioritise pathways that are significantly enriched with
  highly prioritised genes. Prioritised genes are also used to
  identify a gene network interconnecting highly prioritised genes
  and a minimal number of less prioritised genes (which act as
  linkers bringing together highly prioritised genes).

- [Pigengene](https://bioconductor.org/packages/Pigengene) Pigengene
  package provides an efficient way to infer biological signatures
  from gene expression profiles. The signatures are independent from
  the underlying platform, e.g., the input can be microarray or RNA
  Seq data. It can even infer the signatures using data from one
  platform, and evaluate them on the other. Pigengene identifies the
  modules (clusters) of highly coexpressed genes using coexpression
  network analysis, summarizes the biological information of each
  module in an eigengene, learns a Bayesian network that models the
  probabilistic dependencies between modules, and builds a decision
  tree based on the expression of eigengenes.

- [proFIA](https://bioconductor.org/packages/proFIA) Flow Injection
  Analysis coupled to High-Resolution Mass Spectrometry is a
  promising approach for high-throughput metabolomics. FIA- HRMS
  data, however, cannot be pre-processed with current software tools
  which rely on liquid chromatography separation, or handle low
  resolution data only. Here we present the proFIA package, which
  implements a new methodology to pre-process FIA-HRMS raw data
  (netCDF, mzData, mzXML, and mzML) including noise modelling and
  injection peak reconstruction, and generate the peak table. The
  workflow includes noise modelling, band detection and filtering
  then signal matching and missing value imputation. The peak table
  can then be exported as a .tsv file for further analysis.
  Visualisations to assess the quality of the data and of the signal
  made are easely produced.

- [psichomics](https://bioconductor.org/packages/psichomics)
  Automatically retrieve data from RNA-Seq sources such as The Cancer
  Genome Atlas or load your own files and process the data. This tool
  allows you to analyse and visualise alternative splicing.

- [qsea](https://bioconductor.org/packages/qsea) qsea (quantitative
  sequencing enrichment analysis) was developed as the successor of
  the MEDIPS package for analyzing data derived from methylated DNA
  immunoprecipitation (MeDIP) experiments followed by sequencing
  (MeDIP-seq). However, qsea provides several functionalities for the
  analysis of other kinds of quantitative sequencing data (e.g.
  ChIP-seq, MBD-seq, CMS-seq and others) including calculation of
  differential enrichment between groups of samples.

- [RCAS](https://bioconductor.org/packages/RCAS) RCAS is an automated
  system that provides dynamic genome annotations for custom input
  files that contain transcriptomic regions. Such transcriptomic
  regions could be, for instance, peak regions detected by CLIP-Seq
  analysis that detect protein-RNA interactions, RNA modifications
  (alias the epitranscriptome), CAGE-tag locations, or any other
  collection of target regions at the level of the transcriptome.
  RCAS is designed as a reporting tool for the functional analysis of
  RNA-binding sites detected by high-throughput experiments. It takes
  as input a BED format file containing the genomic coordinates of
  the RNA binding sites and a GTF file that contains the genomic
  annotation features usually provided by publicly available
  databases such as Ensembl and UCSC. RCAS performs overlap
  operations between the genomic coordinates of the RNA binding sites
  and the genomic annotation features and produces in-depth
  annotation summaries such as the distribution of binding sites with
  respect to gene features (exons, introns, 5'/3' UTR regions,
  exon-intron boundaries, promoter regions, and whole transcripts).
  Moreover, by detecting the collection of targeted transcripts, RCAS
  can carry out functional annotation tables for enriched gene sets
  (annotated by the Molecular Signatures Database) and GO terms. As
  one of the most important questions that arise during protein-RNA
  interaction analysis; RCAS has a module for detecting sequence
  motifs enriched in the targeted regions of the transcriptome. A
  full interactive report in HTML format can be generated that
  contains interactive figures and tables that are ready for
  publication purposes.

- [rDGIdb](https://bioconductor.org/packages/rDGIdb) The rDGIdb
  package provides a wrapper for the Drug Gene Interaction Database
  (DGIdb). For simplicity, the wrapper query function and output
  resembles the user interface and results format provided on the
  DGIdb website (http://dgidb.genome.wustl.edu/).

- [readat](https://bioconductor.org/packages/readat) This package
  contains functionality to import, transform and annotate data from
  ADAT files generated by the SomaLogic SOMAscan platform.

- [recount](https://bioconductor.org/packages/recount) Explore and
  download data from the recount project available at
  https://jhubiostatistics.shinyapps.io/recount/. Using the recount
  package you can download RangedSummarizedExperiment objects at the
  gene, exon or exon-exon junctions level, the raw counts, the
  phenotype metadata used, the urls to the sample coverage bigWig
  files or the mean coverage bigWig file for a particular study. The
  RangedSummarizedExperiment objects can be used by different
  packages for performing differential expression analysis. Using
  http://bioconductor.org/packages/derfinder you can perform
  annotation-agnostic differential expression analyses with the data
  from the recount project as described at
  http://biorxiv.org/content/early/2016/08/08/068478.

- [regsplice](https://bioconductor.org/packages/regsplice)
  Statistical methods for detection of differential exon usage in
  RNA-seq and exon microarray data sets, using L1 regularization
  (lasso) to improve power.

- [sights](https://bioconductor.org/packages/sights) SIGHTS is a
  suite of normalization methods, statistical tests, and diagnostic
  graphical tools for high throughput screening (HTS) assays. HTS
  assays use microtitre plates to screen large libraries of compounds
  for their biological, chemical, or biochemical activity.

- [signeR](https://bioconductor.org/packages/signeR) The signeR
  package provides an empirical Bayesian approach to mutational
  signature discovery. It is designed to analyze single nucleotide
  variaton (SNV) counts in cancer genomes, but can also be applied to
  other features as well. Functionalities to characterize signatures
  or genome samples according to exposure patterns are also provided.

- [SIMLR](https://bioconductor.org/packages/SIMLR) Single-cell
  RNA-seq technologies enable high throughput gene expression
  measurement of individual cells, and allow the discovery of
  heterogeneity within cell populations. Measurement of cell-to-cell
  gene expression similarity is critical to identification,
  visualization and analysis of cell populations. However,
  single-cell data introduce challenges to conventional measures of
  gene expression similarity because of the high level of noise,
  outliers and dropouts. We develop a novel similarity-learning
  framework, SIMLR (Single-cell Interpretation via Multi-kernel
  LeaRning), which learns an appropriate distance metric from the
  data for dimension reduction, clustering and visualization. SIMLR
  is capable of separating known subpopulations more accurately in
  single-cell data sets than do existing dimension reduction methods.
  Additionally, SIMLR demonstrates high sensitivity and accuracy on
  high-throughput peripheral blood mononuclear cells (PBMC) data sets
  generated by the GemCode single-cell technology from 10x Genomics.

- [SNPediaR](https://bioconductor.org/packages/SNPediaR) SNPediaR
  provides some tools for downloading and parsing data from the
  SNPedia web site <http://www.snpedia.com>. The implemented
  functions allow users to import the wiki text available in SNPedia
  pages and to extract the most relevant information out of them. If
  some information in the downloaded pages is not automatically
  processed by the library functions, users can easily implement
  their own parsers to access it in an efficient way.

- [SPLINTER](https://bioconductor.org/packages/SPLINTER) SPLINTER
  provides tools to analyze alternative splicing sites, interpret
  outcomes based on sequence information, select and design primers
  for site validiation and give visual representation of the event to
  guide downstream experiments.

- [SRGnet](https://bioconductor.org/packages/SRGnet) We developed
  SRMnet to analyze synergistic regulatory mechanisms in
  transcriptome profiles that act to enhance the overall cell
  response to combination of mutations, drugs or environmental
  exposure. This package can be used to identify regulatory modules
  downstream of synergistic response genes, prioritize synergistic
  regulatory genes that may be potential intervention targets, and
  contextualize gene perturbation experiments.

- [StarBioTrek](https://bioconductor.org/packages/StarBioTrek) This
  tool StarBioTrek presents some methodologies to measure pathway
  activity and cross-talk among pathways integrating also the
  information of network data.

- [statTarget](https://bioconductor.org/packages/statTarget) An easy
  to use tool provide a graphical user interface for quality control
  based shift signal correction, integration of metabolomic data from
  multi-batch experiments, and the comprehensive statistic analysis
  in non-targeted or targeted metabolomics.

- [SVAPLSseq](https://bioconductor.org/packages/SVAPLSseq) The
  package contains functions that are intended for the identification
  of differentially expressed genes between two groups of samples
  from RNAseq data after adjusting for various hidden biological and
  technical factors of variability.

- [switchde](https://bioconductor.org/packages/switchde) Inference
  and detection of switch-like differential expression across
  single-cell RNA-seq trajectories.

- [synergyfinder](https://bioconductor.org/packages/synergyfinder)
  Efficient implementations for all the popular synergy scoring
  models for drug combinations, including HSA, Loewe, Bliss and ZIP
  and visualization of the synergy scores as either a two-dimensional
  or a three-dimensional interaction surface over the dose matrix.

- [TVTB](https://bioconductor.org/packages/TVTB) The package provides
  S4 classes and methods to filter, summarise and visualise genetic
  variation data stored in VCF files. In particular, the package
  extends the FilterRules class (S4Vectors package) to define news
  classes of filter rules applicable to the various slots of VCF
  objects. Functionalities are integrated and demonstrated in a Shiny
  web-application, the Shiny Variant Explorer (tSVE).

- [uSORT](https://bioconductor.org/packages/uSORT) This package is
  designed to uncover the intrinsic cell progression path from
  single-cell RNA-seq data. It incorporates data pre-processing,
  preliminary PCA gene selection, preliminary cell ordering, feature
  selection, refined cell ordering, and post-analysis interpretation
  and visualization.

- [yamss](https://bioconductor.org/packages/yamss) Tools to analyze
  and visualize high-throughput metabolomics data aquired using
  chromatography-mass spectrometry. These tools preprocess data in a
  way that enables reliable and powerful differential analysis.

- [YAPSA](https://bioconductor.org/packages/YAPSA) This package
  provides functions and routines useful in the analysis of somatic
  signatures (cf. L. Alexandrov et al., Nature 2013). In particular,
  functions to perform a signature analysis with known signatures
  (LCD = linear combination decomposition) and a signature analysis
  on stratified mutational catalogue (SMC = stratify mutational
  catalogue) are provided.

- [yarn](https://bioconductor.org/packages/yarn) Expedite large
  RNA-Seq analyses using a combination of previously developed tools.
  YARN is meant to make it easier for the user in performing basic
  mis-annotation quality control, filtering, and condition-aware
  normalization. YARN leverages many Bioconductor tools and
  statistical techniques to account for the large heterogeneity and
  sparsity found in very large RNA-seq experiments.

NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages since the last release. The following package NEWS is available:


[ABAEnrichment](https://bioconductor.org/packages/ABAEnrichment)
-------------

Changes in version 1.3.7:

NEW FEATURES

- Add figure to vignette describing the effect of the
  "circ_chrom"-option when genomic regions are used as input

Changes in version 1.3.6:

NEW FEATURES

- Add figure to vignette describing the FWER calculation for the
  hypergeometric test.

Changes in version 1.3.5:

NEW FEATURES

- The new option 'gene_len' creates random sets for the FWER of the
  hypergeometric test dependent on the gene lengths.

- In addition to explicitly naming the candidate genes, it is now
  possible to define entire genomic regions as candidate and background
  regions.  The names of the 'genes' vector have to be of the form
  'chr:start-stop' to use this option.

- The background regions by default are independent. But using the
  option 'circ_chrom = TRUE', background regions on the same chromosome
  as the candidate region are used and random candidate regions are
  allowed to overlap multiple background regions.

SIGNIFICANT USER-LEVEL CHANGES

- cutoff-quantiles are now computed using all protein coding genes and
  not just input candidate and background genes.

Changes in version 1.3.4:

NEW FEATURES

- functions get_name, get_superstructures and get_sampled_substructures
  now also accept structure-IDs without the "Allen:"-prefix as input
  (like data are provided in ABAData-package).

Changes in version 1.3.3:

BUG FIXES

- fixed error in get_expression(..., dataset="5_stages")- message
  ("returning log2(RPKM)" to "returning RPKM").

[affxparser](https://bioconductor.org/packages/affxparser)
----------

Changes in version 1.45.1 (2016-09-16):

- CLEANUP: Dropped obsolete src/R_affx_test.*cmdline.cpp files.

- CLEANUP: Using c(x,y) instead of append(x,y) internally.

Changes in version 1.45.0 (2015-05-03):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.4 for R (>= 3.3.0).

[ALDEx2](https://bioconductor.org/packages/ALDEx2)
------

Changes in version 1.5.2:

- added ability to choose the basis for the clr: all, iqlr, zero or
  user-defined. useful when dealing with asymmetric datasets (selex,
  metagenomics, meta-RNA-seq)

- updated vignette to show how the basis affects the analysis

- made BiocParallel the only parallel package for multicore processing

- made zero-replacement a prior probability rather than a pseudocount

[alpine](https://bioconductor.org/packages/alpine)
------

Changes in version 0.99.6:

- Change estimateAbundance() and predictCoverage() interface such that
  model.names is provided instead of models, using a simple character
  vector. The gene expression term '+ gene' is taken care of
  internally.

Changes in version 0.99.5:

- Store readlength, minsize, maxsize in fitpar and so remove as
  arguments to estimateAbundance() and predictCoverage()

Changes in version 0.99.4:

- Renamed the mysterious estimateTheta() to estimateAbundance()

Changes in version 0.99.3:

- Allow custom knots for GC and relative position

Changes in version 0.99.0:

- Package submission!

[AMOUNTAIN](https://bioconductor.org/packages/AMOUNTAIN)
---------

Changes in version 0.99.1:

NEW FEATURES

- Initial version

USER-LEVEL CHANGES

- More comprehensive help pages

- BiocStyle Vignette

[anamiR](https://bioconductor.org/packages/anamiR)
------

Changes in version 1.0.0:

- Initial release

[Anaquin](https://bioconductor.org/packages/Anaquin)
-------

Changes in version 0.1:

- Created Vignette

- Initial submission to Bioconductor

[AneuFinder](https://bioconductor.org/packages/AneuFinder)
----------

Version: 1.1.6
Category: NEW FEATURES
Text: Added DNAcopy algorithm to Strand-seq mode.

Version: 1.1.6
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: Renamed parameter 'most.frequent.state.bivariate' ->
        'most.frequent.state.strandseq'.

Version: 1.1.6
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: Renamed parameter 'most.frequent.state.univariate' ->
        'most.frequent.state'.

Version: 1.1.6
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: New parameter 'strandseq'.

Version: 1.1.6
Category: BUG FIXES
Text: Dendrogram and heatmap are now aligned properly in
        heatmapGenomewide().

Version: 1.1.5
Category: NEW FEATURES
Text: Aneufinder runs DNAcopy algorithm in addition to the Hidden
        Markov Model.

Version: 1.1.5
Category: NEW FEATURES
Text: New function "getQC" to get a data.frame with quality metrics.

Version: 1.1.5
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: Changed folder structure to include DNAcopy method.

Version: 1.1.5
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: Renamed methods from c('univariate','bivariate') to
        c('HMM','biHMM')

Version: 1.1.4
Category: NEW FEATURES
Text: karyotypeMeasures() has new option regions.

Version: 1.1.4
Category: NEW FEATURES
Text: plotHeterogeneity() for easy plotting of karyotype measures.

Version: 1.1.4
Category: NEW FEATURES
Text: BiocStyle vignette.

Version: 1.1.4
Category: NEW FEATURES
Text: New option use.bamsignals=FALSE/TRUE available for the binning
        step.

Version: 1.1.4
Category: NEW FEATURES
Text: getQC() handles NULL entries as NA and is thus more robust.

Version: 1.1.4
Category: NEW FEATURES
Text: complexity estimation via Michaelis-Menten is carried along.

Version: 1.1.4
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: Color scheme for copy number states has been improved for states
        >= 5-somy.

Version: 1.1.4
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: Option format has been removed in all functions. File format is
        determined automatically now.

Version: 1.1.4
Category: SIGNIFICANT USER-LEVEL CHANGES
Text: clusterByQuality() clusters now on complexity as well by default.

Version: 1.1.4
Category: DEPRECATED AND DEFUNCT
Text:

Version: 1.1.4
Category: BUG FIXES
Text: Corrected bug in order of seqlevels after as(..., 'GRanges').

Version: 1.1.4
Category: BUG FIXES
Text: Corrected bug in hotspotter() that caused detection of
        low-abundance hotspots.

[AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi)
-------------

Changes in version 1.36.0:

BUG FIXES

- fix bug in dbschema method for DBIConnection objects

- keys() no longer fails when global variable 'y' exists

[AnnotationForge](https://bioconductor.org/packages/AnnotationForge)
---------------

Changes in version 1.16.0:

NEW FEATURES

- add check that required db0 package is installed in makeDBPackage()

- add script to create viableIDs data file

MODIFICATIONS

- remove outdated unit test and deprecate MAPCOUNTS in ChipDb package
  templates

- update ChipDb package templates to load DBI in unit tests

- allow build without --keep-empty-dirs flag

- misc code cleanup - No duplicate dependencies in DESCRIPTION -
  resolve symbols in name space - re-usable unit tests (new rather than
  reuse of temporary files) - avoid unnecesary paste(), e.g,
  message(paste()) - use of loadNamespace() rather than require() -
  avoid 1:n iterations - formatting, e.g., of SQL and if () {} else {}

- remove unused tests/runalltests.Rout.save file in AnnDbPkg-templates

- use https: rather than http: for NCBI access

- rename R code and man page files consistent with high level functions

- bump OrgDb version for 3.4 release

- add explicit AnnotationDbi:::NCBIORG_DB_SeedGenerator()

BUG FIXES

- getFastaSpeciesDirs() trims '\r' on Windows

[AnnotationHub](https://bioconductor.org/packages/AnnotationHub)
-------------

Changes in version 2.6.0:

NEW FEATURES

- add vignette section on sharing resources on clusters

- add 'preparerclass' to index.rda to allow search by package name for
  ExperimentHub objects

- add GenomicScoresResource class for Robert Castelo

MODIFICATIONS

- return 'tags' metadata as list instead of comma-separated character
  vector

- move AnnotationHubRecipes vignette to AnnotationHubData

- move listResources() and loadResources() from ExperimentHub

- expose additional fields in .DB_RESOURCE_FIELDS()

- modify cache path to avoid creating a '~' directory on Mac

- use https: NCBI rul in documentation

- modify .get1,EpiExpressionTextResource-method to use 'gene_id' column
  as row names

[AnnotationHubData](https://bioconductor.org/packages/AnnotationHubData)
-----------------

Changes in version 1.4.0:

NEW FEATURES

- add script to generate user-contributed resources

- makeEnsemblGtfToGRanges() no longer stores data in S3 but downloads
  and converts to GRanges on the fly

- add EnsemblFastaTwoBitToAHM unit test

- add man page for makeEnsemblTwoBitToAHM and ensemblFastaToTwoBitFile

- add makeAnnotationHubMetadata() helper

MODIFICATIONS

- move GSE62944-related code to ExperimentHub

- move old vignettes to inst/scripts; add 'Introduction to
  AnnotationHubData' vignette

- remove fasta and towbit files on the fly

- add 'uploadToS3' argument to pushResources() and runRecipes()

- move readMetadataFromCsv() from ExperimentHubData to
  AnnotationHubData

- add 'fileName' arg to readMetadataFromCsv(); don't warn when 'Tags'
  are provided

- specify length for args in readMetadataFromCsv()

- makeAnnotationHubMetadata() populates PreparerClass with package name

- add 'fileName' arg to makeAnnotationHubMetadata()

[annotatr](https://bioconductor.org/packages/annotatr)
--------

Changes in version 0.99.13:

PKG FEATURES

- annotatr is a package to quickly and flexibly annotate genomic
  regions to genomic annotations.

- Genomic annotations include CpG features (island, shore, shelves, and
  open sea), genic features (1-5kb upstream of TSS, promoters, 5'UTRs,
  exons, introns, CDS, 3'UTRs, intron/exon boundaries, and exon/ intron
  boundaries), as well as enhancers from the FANTOM5 consortium for
  hg19 and mm9.

- Annotations are built at runtime using the TxDb.*, AnnotationHub, and
  rtracklayer packages. Users can select annotations a la carte, or via
  shortcuts, such as hg19_basicgenes.

- Annotations are currently available for hg19, mm9, mm10, dm3, dm6,
  rn4, rn5, and rn6. Any species is supported through custom
  annotations.

- Genomic regions are read in using the rtracklayer::import() function,
  and the extraCols argument enables users to include an arbitrary
  number of categorical or numerical data with the genomic regions.

- Annotations are determined via GenomicRanges::findOverlaps(), and all
  annotations are returned, rather than imposing a prioritization.

- annotatr provides several helpful summarization (using dplyr) and
  plot functions (using ggplot2) to investigate trends in data
  associated with the genomic regions over annotations.

[aroma.light](https://bioconductor.org/packages/aroma.light)
-----------

Version: 3.3.2
Date: 2016-09-16
Text: BUG FIX: robustSmoothSpline() gave an error since R-devel (>=
        3.4.0 r70682) (Issue #9)

Version: 3.3.2
Date: 2016-09-16
Text: Using NA_real_ (not just NA) everywhere applicable.

Version: 3.3.1
Date: 2016-08-10
Category: CLEANUP: Using seq_len() and seq_along() everywhere (Issue #8
Text:

Version: 3.3.0
Date: 2016-05-03
Text: The version number was bumped for the Bioconductor devel version,
        which is now BioC v3.4 for R (>= 3.3.0).

[ASpli](https://bioconductor.org/packages/ASpli)
-----

Version: 0.98
Category: NEW FEATURES: multiple bins are reclassified using annotated
        junctions

[attract](https://bioconductor.org/packages/attract)
-------

Version: 1.25.2
Date: 2016-10-13
- Update attract to accept MsigDB data sets in both .gmt and
  .gmx file formats. Before it was just .gmt files. An error is
  thrown if custom gene sets are not in .gmt or .gmx format. Also
  update attract to correct bug when you use reactome database
  and microarray data

Version: 1.25.1
Date: 2016-05-23
- Corrected bug found when using MsigDB data sets instead of
  KEGG and reactome. Also added expressionSetGeneFormat to
  findAttractors and calcFuncSynexprs functions. The vignette was
  also changed to reflect the update

[BaalChIP](https://bioconductor.org/packages/BaalChIP)
--------

Changes in version 0.99.0:

NEW FEATURES

- first submission to Bioc-devel

[BgeeDB](https://bioconductor.org/packages/BgeeDB)
------

Changes in version 2.0.0 (2016-10-10):

- Implemented possibility to deal with different Bgee releases.

- Improved storage and versioning of cached files.

- Implemented use of API key to query our servers in order to prevent
  overloading and spamming.

- Improved management of downloading errors.

- Harmonized the use of a Bgee class object by all functions of the
  package. For example, loadTopAnatData() now requires an input Bgee
  class object to specify species, dataType and pathToData arguments.

- Added input Bgee class object to output of loadtopAnatData()
  function.

- Created new getAnnotation(), getData() and formatData() independent
  functions to replace the Bgee class methods get_annotation(),
  get_data() and format_data().

- In formatData() function, when affymetrix data is used, the "stats"
  parameter is automatically set to "intensities".

- Added possibility to reproduce an analysis offline if all data files
  were previously downloaded in cache.

- Fixed data frames headers including spaces to more convenient headers
  with spaces replaced by dots.

- Harmonized use of camelCase in functions arguments.

- Added argument allowing to sort result table in makeTable() function.

- Implemented management of TPMs as expression unit in future Bgee
  releases.

- Updated vignette.

Changes in version 1.0.3 (2016-08-31):

- Update of format_data() function to output an Expression Set object.

- Fixed makeTable "FDR" column which was a factor instead of a numeric.

- Fixed get_data() and format_data() functions, which did not work when
  multiple chip types were available for an experiment.

[Biobase](https://bioconductor.org/packages/Biobase)
-------

Changes in version 2.33:

BUG FIXES

- exprs<- enforces value with correct dim, dimnames.

[bioCancer](https://bioconductor.org/packages/bioCancer)
---------

Changes in version 1.0.03:

- update pivotr from radiant.data

Changes in version 1.0.02:

- replace add_rownames() by rownames_t_column()

Changes in version 1.0.01:

- rm "id=" from navbarMenu - shiny update

Changes in version 1.0.0:

- built vignette with knitr

- modify whichGeneList(GeneListLabel) resolve difference between
  bioCancer server and package FIRST RELEASE

- Package released

- Omit R/code menu: request shinyAce (>=0.2.1) - Omit help_and_report
  function

- Various issues request DT (>=0.1.39) - Omit Show Plot in Pivot
  sidebar menu - Can not download/Store correctly filtered table in
  Handle/View PERSPECTIVE

- Release R 3.4 import(shiny, except= c("dataTableOutput",
  "renderDataTable"))

- use DisGeNet server despite /extdata/disGeNet file (DONE Not Working)

- clusterProfiler GMT file:
  http://www.r-bloggers.com/go-analysis-using-clusterprofiler/

- gene classification using rpart,ggplot2, ggtree:
  http://guangchuangyu.github.io/2016/01/annotate-a-phylogenetic-tree-with-insets/

- Specify the range of mutation frequency in circos plot [Min, Max]

[BiocInstaller](https://bioconductor.org/packages/BiocInstaller)
-------------

Changes in version 1.24.0:

NEW FEATURES

- biocLite() uses lib.loc= to find devtools, reports more informatively
  why devtools fails to load

- biocLite() only offers to update non-masked packages

- biocLite() reports when packages in unwriteable directories are
  out-of-date, but does not try (and fail) to update them.

- isDevel() returns TRUE if the version of BiocInstaller corresponds to
  the development version of Bioconductor.

[BiocStyle](https://bioconductor.org/packages/BiocStyle)
---------

Changes in version 2.2.0:

NEW FEATURES

- New Bioconductor HTML Style. See package vignettes for details.

BUG FIXES AND IMPROVEMENTS to Bioconductor LaTeX Style 2

- Use `\path` for file names to allow long line breaks

- Load 'nowidow' LaTeX package to prevent widows and orphans

- Patch bug in 'titlesec' 2.10.1
  (http://tex.stackexchange.com/q/299969/102422)

- Pass option `multiple` to 'footmisc' for better handling of
  consecutive footnotes

- Load 'marginfix' LaTeX package to prevent margin notes from
  overflowing the bottom margin

- Fix the issue with color spilling out on margin notes
  (https://github.com/Bioconductor/BiocStyle/issues/5)

- Use `fig.asp` to override figure height
  (https://github.com/Bioconductor/BiocStyle/issues/4)

- Fix compatibility with the 'float' package, in particular the `[H]`
  placement specifier

- Load 'marginfix' LaTeX package to prevent margin notes from
  overflowing the bottom margin

- Enclose wide floats in `\blockmargin` and `\unblockmargin` to prevent
  footnotes from entering them

- Move the footnote mark inside margin notes

- Add vertical skip after margin phantoms of wide floats for better
  alignment of margin notes with paragraph text

- Fix concatenation of `includes` reported in
  https://github.com/Bioconductor/BiocStyle/issues/8

- Stratify parnote mark definition depending on package version
  (https://github.com/Bioconductor/BiocStyle/issues/7)

- Capitalize default opening words in `\comment`, `\warning` and
  `\fixme`, and mention the optional argument in the vignette

[biomaRt](https://bioconductor.org/packages/biomaRt)
-------

Changes in version 2.30.0:

SIGNIFICANT USER-LEVEL CHANGES

- Updated vignette to use BiocStyle and execute most code chunks.

[biosigner](https://bioconductor.org/packages/biosigner)
---------

Changes in version 1.1.14:

INTERNAL MODIFICATIONS

- Biobase import restricted to ExpressionSet, exprs, and pData
  function/methods to avoid warning (conflict on 'combine' with
  randomForest)

Changes in version 1.1.12:

INTERNAL MODIFICATIONS

- minor modification in unit tests

Changes in version 1.1.10:

NEW FEATURE

- biosign can now be applied to an ExpressionSet object

- vignette in html format

INTERNAL MODIFICATIONS

- documentation generated with roxygen2 and vignette with rmarkdown

Changes in version 1.1.8:

INTERNAL MODIFICATIONS

- error message fixed when missing values in the accuracy vector
  obtained by the bootstrap step

Changes in version 1.1.6:

INTERNAL MODIFICATIONS

- versioning update

- unit tests silenced on windows platforms because of errors on the
  moscato2 bioc platform running on windows 8

Changes in version 1.1.4:

INTERNAL MODIFICATIONS

- bug fixed (when tierMN contains 0 only)

- 'show' method: better handling of messages when no signature is found

Changes in version 1.1.2:

INTERNAL MODIFICATIONS

- unit tests: test_biosign_diaplasma and test_biosign_sacurine added

- internal renaming of variables (to indicate their type) and functions
  (to facilitate their understanding)

- PLS-DA: to avoid errors during generation of models, the number of
  predictive components is at least 1

Changes in version 1.1.0:

PACKAGE MODIFICATION

- Wellcome to the biosigner package for feature selection from omics
  datasets

- The package implements a new wrapper method detecting the features
  which are important for PLS-DA, Random Forest, or SVM binary
  classification

- The package contains the 'diaplasma' LC-MS metabolomics real dataset
  (plasma samples from diabetic type 1 and 2 patients)

- Please see the vignette for details about the approach and package
  use

- The corresponding publication is currently under review.

[BrowserViz](https://bioconductor.org/packages/BrowserViz)
----------

Changes in version 1.9.8:

BUG FIXES

- Significant (3x) speedup.  A 5000-node, 6000-edge graph transmits to
  Cytoscape from R in about 20 seconds.

Changes in version 1.8.0:

NEW FEATURES

- setNodeOpacityRule, controlling node fill color, border and/or label;
  interpolate & lookup modes both supported

- getNodeSize

- saveImage now supports pdf as well as png and svg formats

- setDefaultEdgeFontSize

- getAdjacentEdgeNames

SIGNIFICANT USER-VISIBLE CHANGES

- changed method names: layout -> layoutNetwork, version ->
  pluginVersion, get/setPosition -> get/setNodePosition

- NAMESPACE now imports four more methods from the graph package,
  helpful for package developers using RCytoscape: edgemode, addNode,
  addEdge, requested by Robert Flight.

BUG FIXES

- Changed getNodePosition node.name.delimiter to eliminate regex token,
  from ':.:' to ':-:' saveLayout now has optional 3rd parameter,
  'timestamp.in.filename'

- Fixed bug in setNodeLabelDirect.  Multiple nodes, one label now
  works.

- setCenter now casts x,y to numeric before sending out to CyRPC

[CAMERA](https://bioconductor.org/packages/CAMERA)
------

Changes in version 1.29.2:

BUG FIXES

- Fix build issues and Vignette encoding (again)

- Fix build issues and Vignette encoding

Changes in version 1.29.1:

NEW FEATURES

- Added function findIsotopesWithValidation

[canceR](https://bioconductor.org/packages/canceR)
------

Version: 1.5.2
- resolve issue 1 about "tl not found
- delete file R/getGeneticProfiles.R
- resolve 'no visible binding for global variable
- add R/aaa.R

[Cardinal](https://bioconductor.org/packages/Cardinal)
--------

Changes in version 1.5.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated 'batchProcess' to support reduceDimension and peakAlign

- Now 'peakAlign' looks for an existing 'mean' column in featureData

- Added 'matter' support to readAnalyze (previously only readImzML)

BUG FIXES

- Fixed bug when indexing into data cube using 'imageData' method

Changes in version 1.5.1:

BUG FIXES

- Corrected author and maintainer contact information

Changes in version 1.5.0:

NEW FEATURES

- Added experimental support for 'matter' on-disk matrices, from
  package 'matter', hosted at https://github.com/kuwisdelu/matter, as a
  replacement for 'Binmat' matrices

BUG FIXES

- Fixed subsetting SImageData objects with variables

[ChAMP](https://bioconductor.org/packages/ChAMP)
-----

Changes in version 2.0.1:

- Updating to ChAMP2.

[Chicago](https://bioconductor.org/packages/Chicago)
-------

Version: 1.1.5

- We fixed a bug in .addTLB() (called by
  estimateTechnicalNoise()) that reduced the number of
  trans-count bins used, though the read counts in those bins
  were still calculated correctly. Since fixing this bug
  increased the number of trans-count bins, we adjusted the
  default values of two settings accordingly

- tlb.minProxOEPerBin: 1,000 changed to 50,000

- tlb.minProxB2BPerBin: 100 changed to 2,500

- These values have been chosen to ensure that results should
  change as little as possible from Chicago 1.1.4. Unless you
  were using custom values, you should not notice any qualitative
  differences

- If you need to re-run chicagoPipeline() (or just
  estimateTechnicalNoise()) on a chicagoData object created using
  version 1.1.4 or earlier, please manually update the parameters
  to their new settings

- cd <- modifySettings(cd, settings=list(tlb.minProxOEPerBin=50000,
        tlb.minProxB2BPerBin=2500))

- Many thanks to Thomas Sexton for bringing the bug to our attention and
  helping us fix it

[chipenrich](https://bioconductor.org/packages/chipenrich)
----------

Changes in version 1.11.4:

PKG FEATURES

- Back-end refactoring of the code base.

- Break out functions in main.R into .R files collecting similar
  functions.

- Transition documentation to roxygen2 blocks.

- Improve commenting in chipenrich() function.

- Assigning peaks using GenomicRanges object rather than than list of
  IRanges.

- Rewrite package vignette in Rmarkdown and render with knitr.

- Improve supported_*() functions to report and check combinations of
  genome, organism, genesets, locusdef, and mappability read length.

- Cleanup DESCRIPTION and NAMESPACE to avoid loading entire packages.

- Follow data() best practices.

[ChIPpeakAnno](https://bioconductor.org/packages/ChIPpeakAnno)
------------

Changes in version 3.7.9:

- fix the problem in the test when there are names for target but not
  for current.

Changes in version 3.7.8:

- remove the dependence of MMDiffBamSubset in documentation becuase it
  is not availble now.

Changes in version 3.7.7:

- fix bugs for toGRanges to import MACS2 broad calls.

Changes in version 3.7.6:

- add new functions tileGRanges, tileCount

Changes in version 3.7.5:

- add parameter select to annoPeaks

Changes in version 3.7.3:

- add new function featureAlignedExtendSignal, estLibSize,
  estFragmentLength

Changes in version 3.7.2:

- Correct typo in documentation for Z-score.

- Improve the efficiency of findOverlapsOfPeaks.

[ChIPQC](https://bioconductor.org/packages/ChIPQC)
------

Changes in version 1.9.2:

- Dependancy cleanup; update handling of how control samples are added

[ChIPseeker](https://bioconductor.org/packages/ChIPseeker)
----------

Changes in version 1.9.8:

- plotAvgProf/plotAvgProf2 order of panel by names of input tagMatrix
  List <2016-09-25, Sun>

- test ENSEMBL ID using '^ENS' instead of '^ENSG' <2016-09-20, Tue> +
  https://github.com/GuangchuangYu/ChIPseeker/issues/41

Changes in version 1.9.7:

- unit test <2016-08-16, Tue>

Changes in version 1.9.6:

- update vignette <2016-08-16, Tue>

Changes in version 1.9.5:

- when TxDb doesn't have gene_id information, converting gene ID
  (ensembl/entrez and symbol) will be omitted instead of throw error.
  <2016-08-02, Tue> + https://www.biostars.org/p/204142

- bug fixed if testing targetPeak is a list of GRanges objects in
  enrichPeakOverlap function <2016-07-20, Wed> +
  https://github.com/GuangchuangYu/ChIPseeker/issues/37 +
  https://github.com/GuangchuangYu/ChIPseeker/issues/36

- fixed typo in determine gene ID type <2016-06-21, Tue> +
  https://github.com/GuangchuangYu/ChIPseeker/issues/28#issuecomment-227212519

- move upsetplot generics to DOSE and import from DOSE to prevent
  function name conflict <2016-06-14, Tue>

Changes in version 1.9.4:

- bug fixed <2016-06-08, Wed> +
  https://github.com/GuangchuangYu/ChIPseeker/issues/17#issuecomment-224407402
  + https://github.com/GuangchuangYu/ChIPseeker/pull/24/files

Changes in version 1.9.3:

- use byte compiler <2016-05-18, Wed>

- https://github.com/Bioconductor-mirror/ChIPseeker/commit/f1ada57b9c66a1a44355bbbbdaf5b0a88e10cf7d

Changes in version 1.9.2:

- name tagMatrix in plotAvgProf automatically if missing <2016-05-12,
  Thu>

- https://github.com/Bioconductor-mirror/ChIPseeker/commit/d5f16b2bc01725e30282c3acb33007ef521a514c

Changes in version 1.9.1:

- bug fixed in getNearestFeatureIndicesAndDistances <2016-05-11, Wed> +
  correct metadata in dummy NA feature

[ClassifyR](https://bioconductor.org/packages/ClassifyR)
---------

Changes in version 1.8.0:

- Ordinary k-fold cross-validation option added.

- Absolute difference of group medians feature selection function
  added.

[clonotypeR](https://bioconductor.org/packages/clonotypeR)
----------

Changes in version 1.12.0:

OTHER CHANGES

- Use Roxygen to generate documentation.

[clustComp](https://bioconductor.org/packages/clustComp)
---------

Changes in version 1.2.0:

- Added new arguments ‘expression’, ‘layout’ and ‘ramp’ to
  flatVShier(), to provide more flexibility to the plot layout. NEW
  FEATURE:

- The expanded version of the flatVShier() function allows adding the
  heatmap of the data, ordered according to the resulting hierarchical
  tree.

[clusterExperiment](https://bioconductor.org/packages/clusterExperiment)
-----------------

Changes in version 0.99.3 (2016-07-26):

Changes

- plot in mergeClusters now uses cluster names and colors from
  clusterLegend

- plotDendrogram now calls plot.phylo

- add 'clusterLabel' argument to `clusterSingle`

- add options 'mad' and 'cv' to the dimensionality reduction. Also made
  option to only use clustered samples for feature reduction for
  relevant functions (e.g. `makeDendrogram`).

- clusterSingle now always returns the D matrix to the slot
  coClustering (previously only did so if D was from subsampling).

- change so that clusterSingle takes dissimilarity matrix, and now
  clusterMany calculates dissimilarities up front (rather than
  recalculating each time)

- add RSEC function for wrapper that leads to RSEC algorithm.

- add test for clusterMany to make sure replicable with past results
  (not unit test because too long to run, so not part of R build)

Bug fixes

- fix bug in .TypeIntoIndices so that handles mix of clusterType and
  clusterLabels in whichClusters

- fixed bug in plotCoClustering so handles clusterSamplesData

- D for clusterD is now distance, not similarity, for 0-1, meaning
  larger values are values that are less similar.

- fix bug in plotClusters that would give clusterLegend entries that
  were vectors, not matrices.

Changes in version 0.99.1 (2016-05-24):

Changes

- changes to pass development version of bioConductor checks.

Changes in version 0.99.0 (2016-05-24):

Changes

- changed number to indicate bioconductor submission

Changes in version 0.2.0 (2016-05-10):

Changes

- Allow 'whichCluster'/'whichClusters' arguments to match to
  clusterLabels, not just clusterTypes

- Added slot 'dendro_index'

- Added 'whichCluster' argument to `makeDendrogram`

- Added 'hierarchicalK' clustering

- Added default distance for 0-1 clustering

- Added ability to define distance for clustering

- Added 'setToCurrent' and 'setToFinal' options to update status of a
  cluster.

- Added unit tests for workflow function (in test_constructor)

- 'getBestFeatures' now calls 'clusterContrasts' internally

- Output for 'clusterContrasts' changed

- Removed 'Index' output for getBestFeatures

- Changed tests for getBestFeatures to run on standard objects (which
  means now have -2 values to test against)

- User can now give clusterLabel for resulting cluster of combineMany
  and mergeClusters

Changes in version 0.1.0 (2016-05-04):

Changes

- Conversion to S4 language for bioConductor submission

- All previous functions have been overhauled, renamed, etc.

Changes in version 0.0.0.9006:

Changes

- fixed so that mergeClusters, clusterHclust, and getBestFeatures will
  appropriately convert if the input of clustering vector is a factor
  rather than numeric (with warning).

- fixed mergeClusters to have option to indicate that input matrix is a
  count matrix (in which case will create dendrogram with log(counts+1)
  and will do getBestFeatures with the voom correction)

- added more tutoral-oriented vignette (old vignette is now the
  documentation vignette with more detail about the internal workings
  of package). Currently is just simulated data, but will be updated to
  real single-cell sequencing dataset.

Changes in version 0.0.0.9005:

Changes

- Changed simulated data so load all with data(simData) rather than
  separate calls for simData and simCount. Also added 'trueCluster'
  vector to give true cluster assignments of simulated data

- added dendro example to getBestFeatures

- added example to clusterHclust

- added single function for converting to phylobase tree (used
  internally by package)

- added functionality to find proportion of significant null hypotheses
  for merging clusters (mergeClusters)

Changes in version 0.0.0.9004:

Changes

- Changed clusterMany.R to only set k<-NA if sequential=FALSE
  (previously for all where findBestK=TRUE)

- Added to vignette

- fixed bug in plotClusters to correctly plot "-1"

[clusterProfiler](https://bioconductor.org/packages/clusterProfiler)
---------------

Changes in version 3.1.9:

- bug fixed of simplify method for compareCluster object <2016-08-17,
  Mon> + use semData parameter according to the change of GOSemSim

Changes in version 3.1.8:

- as.data.frame method for compareClusterResult <2016-09-29, Thu>

- geneID and geneInCategory methods <2016-09-19, Mon>

Changes in version 3.1.7:

- [, [[, head, tail, dim methods for compareClusterResult <2016-08-28,
  Sun>

- fixed R-devel check <2016-08-16, Tue>

Changes in version 3.1.6:

- update vignette <2016-08-16, Tue>

- browseKEGG <2016-08-15, Mon>

- unit test <2016-08-15, Mon>

Changes in version 3.1.5:

- move enrichMeSH & gseMeSH to meshes <2016-08-11, Thu>

Changes in version 3.1.4:

- MeSH term enrichment analysis (enrichMeSH and gseMeSH) <2016-07-28,
  Thu>

Changes in version 3.1.3:

- export download_KEGG <2016-07-25, Mon>

- according to the changes of GOSemSim and DOSE <2016-07-05, Tue>

Changes in version 3.1.2:

- 'by' parameter for GSEA analysis <2016-07-04, Mon>

- duplicated KEGG path id for KEGG orthology, only use ko <2016-06-28,
  Tue> + e.g. map00010 and ko00010 +
  http://www.kegg.jp/dbget-bin/www_bget?ko+ko00010 +
  http://www.kegg.jp/dbget-bin/www_bget?ko+map00010

Changes in version 3.1.1:

- getGOLevel and dropGO support using a vector of GO level <2016-05-20,
  Fri> +
  https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/clusterprofiler/zNYHLI6RpJI/1u4XRKfNAQAJ

- use byte compiler <2016-05-18, Wed>

- https://github.com/Bioconductor-mirror/clusterProfiler/commit/2949fe8db5464c419c9d1665cb24f43c49ac54ec

[ClusterSignificance](https://bioconductor.org/packages/ClusterSignificance)
-------------------

Changes in version 1.2.0:

NEW FEATURES

- Added a slot for group colors for the classes Pcp, Mlp,
  ClassifiedPoints and PermutationResults

[CNEr](https://bioconductor.org/packages/CNEr)
----

Changes in version 3.4:

NEW FEATURES

- Updated CNE class for storing all the information about running the
  pipeline.

- Add read.rmMask.GRanges to read RepeatMasker .out file.

- Add read.rmskFasta to read soft repeat masked fasta file.

- Add the distribution plot of axt alignment matches.

- Add the distribution plot of CNE length.

- Add the syntenicDotplot for axt alignment and GRangePairs.

- When readAxt and readBed, seqinfo is kept when available.

- New fixCoordinates function makes the coordinates of Axt alignments
  always relative to positive strands.

- Add parallel subAxt for Axt alignment.

- Add the plot of genomic distribution of CNE.

- Add the function to make bed and bigwig files of CNEs.

BUG FIXES

- Instead of an error, an empty GRangePairs is returned when no CNEs
  identified.

Changes in version 3.3:

BUG FIXES

- Fix a bug caused by "format" in the blat step.

NEW FEATURES

- Add the pairwise whole genome alignment pipeline

- Add a new class "GRangePairs"

- "Axt" class is now based on "GRangePairs" class.

- readAncora for reading Ancora format CNE files.

[coMET](https://bioconductor.org/packages/coMET)
-----

Changes in version 1.5.5 (2016-06-12):

- Update the script for coMET website

- remove big files

- Add examples data that was not uploaded correctly the first time

Changes in version 1.5.3:

- -Change the name of functions: xenorefGenesUCSC into
  xenorefGenes_UCSC transcriptENSEMBL into transcript_ENSEMBL
  structureBiomart into structureBiomart_ENSEMBL snpLocationsUCSC into
  snpLocations_UCSC snpBiomart into snpBiomart_ENSEMBL
  RepeatMaskerTrack into repeatMasker_UCSC regulationBiomart into
  regulationBiomart_ENSEMBL knownGenesUCSC into knownGenes_UCSC
  ISCATrack into ISCA_UCSC HistoneOne into HistoneOne_UCSC HistoneAll
  to HistoneAll_UCSC GWASTrack to GWAScatalog_UCSC genesNameENSEMBL to
  genesName_ENSEMBL genesENSEMBL to genes_ENSEMBL GeneReviewsTrack to
  GeneReviews_UCSC gcContent to gcContent_UCSC GADTrack to GAD_UCSC
  DNAseUCSC to DNAse_UCSC COSMICTrack to COSMIC_UCSC cpgIslandsUCSC to
  cpgIslands_UCSC CoreillCNVTrack to CoreillCNV_UCSC chromatinHMMAll to
  chromatinHMMAll_UCSC chromatinHMMOne to chromatinHMMOne_UCSC
  ClinVarCnvTrack to ClinVarCnv_UCSC ClinVarMainTrack to
  ClinVarMain_UCSC

Changes in version 1.4.3 (2016-04-29):

- Update the management of multiple secondary info data

- Add the update done in the devel version that was not updated in the
  new release for x reasons

[compEpiTools](https://bioconductor.org/packages/compEpiTools)
------------

Changes in version 1.7.4:

- The following function is updated + GRcoverageInbins: A small bug in
  the GRcoverageInbins method was fixed: the problem resulted in
  incorrect coverage (reported as NA) in the last bin(s) of short
  intervals, those for which the number of bins is similar to their
  width.

Changes in version 1.7.3:

- The following function is updated + topGOres: Updated the bug in
  topGOres that was returning very low p-values as characters instead
  of numeric values.

[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap)
--------------

Changes in version 1.11.8:

- `anno_barplot()`: accept a matrix as input to plot stacked barplots

Changes in version 1.11.7:

- SingleAnnotation: if `col` is a vector with no names, it will be
  assigned as `level(value)` or `unique(value)`

- HeatmapAnnotation: give warnings if color is defined while with no
  annotations

- HeatmapAnnotation: check `col`, if it is not valid, give warnings

- catch error when making annotation graphics

Changes in version 1.11.6:

- simply bump the verison number

- `gap` in `Heatmap()` now can be a vector

Changes in version 1.11.5:

- `gap` in `HeatmapAnnotation` has been adjusted

- annotations support drawing names of either sides

- `densityHeatmap()`: quantile lines are also reordered

- export `anno_oncoprint_barplot`

- `Heatmap()`: if `col` is a unnamed vector and the number of colors is
  same as unique itemes in `mat`, the name of `col` vector is set to
  `sort(unique(mat))

- adjusted the order of annotation legends

- discreat legend: if a level is not in the data while defined by
  `col`, it will be removed.

Changes in version 1.11.2:

- `grid.dendrogram()`: do not draw dendrogram if the height is zero

- `densityHeatmap()`: support clustering on columns and more controls
  on column settings

Changes in version 1.11.1:

- `draw,HeatmapList-method` can control row order and clustering of the
  main heatmap

[CountClust](https://bioconductor.org/packages/CountClust)
----------

Changes in version 1.1.3:

- Release We have added a new function `ExtractHighCorFeatures.R` to
  track the top correlated features with expression data per topic. We
  added an option for making StructureGGplot without phenotype
  information. For the FitGoM function that fits the Grade of
  Membership model, added a maximum number of iterations input, set to
  10,000 as default but flexible to user change.

[crisprseekplus](https://bioconductor.org/packages/crisprseekplus)
--------------

Version: 0.99.0
Text:

[CrispRVariants](https://bioconductor.org/packages/CrispRVariants)
--------------

Changes in version 1.1.6:

- plotAlignments can now mark codon boundaries if codon frame is
  specified.

- Added citation

Changes in version 1.1.5:

- More flexible specification of strand with new readsToTarget
  parameter 'orientation'

- Fixed warning caused by implicit embedding of S4 objects

- Added tests for 'plotAlignments' and 'annotateGenePlot'

- Minor speedup and internal restructuring of 'annotateGenePlot'

- Added CRISPR biocView

- Changed NEWS to rd format

Changes in version 1.1.4:

- Fixed a bug that prevented SNV settings being used in some
  circumstances

Changes in version 1.1.2:

- new function consensusSeqs returns the consensus sequences of the
  variant alleles

[csaw](https://bioconductor.org/packages/csaw)
----

Changes in version 1.7.4:

- Added protection against NA values in filterWindows().

- Deprecated the use of parameter lists in any param= arguments.

- Tightened up allowable values of ext= arguments in various functions.

- Added the BPPARAM slot in the readParam class to store
  BiocParallelParam objects.

- Added support for parallelization in windowCounts(), regionCounts()
  and others.

- Updated documentation, user's guide.

[customProDB](https://bioconductor.org/packages/customProDB)
-----------

Changes in version 1.13.2:

UPDATED FUNCTIONS

- Update functions PrepareAnnotationEnsembl.R and
  PrepareAnnotationRefseq.R due to updates in depending packages

Changes in version 1.13.1:

BUG FIXES

- Fix a small bug in function JunctionType.R

NEW FEATURES

- Update the function OutputVarproseq.R to output a data frame
  'snvproseq'

- Add a function OutputVarprocodingseq.R to output a data frame
  'snvprocoding'

[cytofkit](https://bioconductor.org/packages/cytofkit)
--------

Changes in version 1.5.10 (2016-10-05):

MODIFICATION

- corrected the citation title

- added cluster filter in rateChange line plot in ShinyAPP

- debugged the error of loading back exported RData file

- updated the saving button, make it more robust. If cannot find FCS
  path, then doesn't save new FCS files

- added progression indicator in shinyAPP

Changes in version 1.5.9 (2016-10-03):

MODIFICATION

- debug when w is negative in autoLgcl transformation, convert to
  logicle transformation in this case

- modify the cytof_writeResults function to make it more robust.

- set autoLgcl as the default transformation method.

- add cytofkit plos computational biology paper in CITATION

Changes in version 1.5.8 (2016-08-19):

MODIFICATION

- add none option for transformation method to support FCS files with
  data already transformed

- big updates on the layout of shinyAPP, functions categorized into
  four panels ("cluster", "marker", "sample", "progression")

- add "Seperate Plot by Samples" on side panel, remove option "Label
  Samples by Shapes"

- make sample filter works on subset progression panel, support plot
  labeld on grid plot

- added "group samples" function to relabel and group samples in sample
  panel

- added subset precentage change plot in sample panel

- added case checking of Nan of w in autoLgcl function

Changes in version 1.5.7 (2016-07-13):

MODIFICATION

- Smart handle of non-transformation markers(FSC,SSC), and exclude of
  channels(like Event, Time) when markers==NULL specified in
  cytof_exprsExtract function

NEW FEATURES

- Redesigned the shinyAPP tab panels: Cluster Plot, Marker Plot, Subset
  Progression

- Added cluster filtering and cluster table on shinyAPP diffusionmap
  set up page

- Added combined view of marker expression patten on scatter plot and
  marker expression trend on subset progression

- Added stack density plot in Marker Plot

- Added cluster annotation (Label Clusters) in Marker Plot

Changes in version 1.5.6 (2016-07-08):

MODIFICATION

- Adjusted the windown width and height of getParameters_GUI() to fit
  long name of FCM data

- debugged, set full.names = TRUE in getParameters_GUI() when fcsFile
  == NULL

Changes in version 1.5.5 (2016-07-04):

MODIFICATION

- Concised the autoLgcl function

Changes in version 1.5.4 (2016-06-14):

MODIFICATION

- concised the title in vignette, tiny modification

- added missed halo variable in ClusterX

- added transformation=FALSE in calling read.FCS function to avoid
  unexpected transformation for flow cytometry data.

- changed cast from reshap to dcast from reshap2 in function
  cytof_writeResults

- added projectName, rawFCSdir, resultDir entries in analysis_results
  object

- modified parameters for cytof_writeResults, only need
  analysis_results object

- replaced parameter uniformClusterSize in cytofkit to
  clusterSampleSize

NEW FEATURES

- redesigned shiny APP, big updates on subset progression tab, support
  FlowSOM and Diffusion map running. Added save buttion.

- rewrited most of the codes in function cytof_progression

- added diffusionmap in cytof_progression, updated on GUI

- added reverseOrder option in cytof_progressionPlot function

- added clusterLabelSize option in cytof_progressionPlot function

- added segmentSize option in cytof_progressionPlot function

- added cluster filetering and addClusterLabel option in
  cytof_progressionPlot

- added fixCoord option to function cytof_clusterPlot

- added distMethod option in cytof_progression function

- added distance calculation options in cytof_dimReduction

- added tsneSeed in cytof_dimReduction for reproducible t-SNE results

- added cytof_clusterStat and cytof_colorPlot function in
  cytof_postProcess

- added a button on GUI to open the resultDir once the
  cytof_writeResults was done, cross platform.

[CytoML](https://bioconductor.org/packages/CytoML)
------

Version: 1.0.0
Category: First submission
Text:

[DChIPRep](https://bioconductor.org/packages/DChIPRep)
--------

Changes in version 1.3.3:

- updated CITATION file to point to the published paper

Changes in version 1.3.2:

- Bioc 3.4 devel version

[debrowser](https://bioconductor.org/packages/debrowser)
---------

Changes in version 1.1.11:

- Bug fixes

Changes in version 1.1.10:

- Normalization method selection added to main plots and tables too

- Batch effect correction

Changes in version 1.1.9:

- edgeR and d3heatmap dependancy fix

- undocumented functions fix

Changes in version 1.1.8:

- Various bug fixes

- Data preperation is reachable from everywhere

Changes in version 1.1.7:

- Limma, EdgeR and DESeq2 support

- Interactive Heatmap

- Interactive PCA

- IQR and Density plots

[DeepBlueR](https://bioconductor.org/packages/DeepBlueR)
---------

Version: 0.99.0
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

[DEGreport](https://bioconductor.org/packages/DEGreport)
---------

Version: 1.9.9
Text: 10-14-2016 Lorena Pantano <lorena.pantano@gmail.com> Fixes: Fix
        the reduccion of clusters to use correlation values directly.

Version: 1.9.9
Text: 10-06-2016 Lorena Pantano <lorena.pantano@gmail.com> Fixes: GO
        enrichment function. Input genes was wrong.

Version: 1.9.7
Text: 07-29-2015 Lorena Pantano <lorena.pantano@gmail.com> Features:
        Add functions to return markdown report from DESeq2 object, and
        clustering function for time course data. Fixes: small bugs
        related to Nan values or small diversity values in QC plots.

[derfinder](https://bioconductor.org/packages/derfinder)
---------

Changes in version 1.7.16:

SIGNIFICANT USER-VISIBLE CHANGES

- Help pages now document advanced arguments.

- Deprecated advancedArg().

Changes in version 1.7.14:

NEW FEATURES

- Added the function getTotalMapped() for calculating the total number
  of mapped reads for a BAM file or the area under the curve (AUC) for
  a BigWig file. This information can then be used with fullCoverage(),
  filterData() and other functions. Note that if you 'totalMapped' in
  fullCoverage() you should not use 'totalMapped' again in
  filterData().

Changes in version 1.7.12:

BUG FIXES

- Updated links to BrainSpan. Issue reported by Steve Semick
  https://github.com/SteveSemick.

Changes in version 1.7.2:

BUG FIXES

- Now derfinder uses DataFrame(check.names = FALSE) to avoid naming
  issues.

Changes in version 1.7.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Dropped defunct functions.

[derfinderPlot](https://bioconductor.org/packages/derfinderPlot)
-------------

Changes in version 1.7.10:

SIGNIFICANT USER-VISIBLE CHANGES

- Help pages now document advanced arguments.

Changes in version 1.7.8:

BUG FIXES

- Updated links to BrainSpan. Issue reported by Steve Semick
  https://github.com/SteveSemick.

Changes in version 1.7.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Dropped defunct functions.

[DESeq2](https://bioconductor.org/packages/DESeq2)
------

Changes in version 1.13.8:

- Use a linear model to estimate the expected counts for dispersion
  estimation in estDispGeneEst() if the number of groups in the model
  matrix is equal to the number of columns of the model matrix. Should
  provide a speed-up for dispersion estimation for model matrices with
  many samples.

Changes in version 1.13.3:

- Fixed bug: fpm() and fpkm() for tximport.

- Fixed bug: normalization factors and VST.

- Added an error if tximport lengths have 0.

- Added an error if user matrices are not full rank.

- More helpful error for constant factor in design.

[DiffBind](https://bioconductor.org/packages/DiffBind)
--------

Changes in version 2.2.0:

- Feature: Control which principal components are plotted using
  components parameter in dba.plotPCA

- Feature: Control axis range using xrange and yrange parameters in
  dba.plotMA

- Feature: Filtering per-contrast using filter and filterFun parameters
  in dba.analyze)

- Feature: Flip which group in contrast shows gain/loss (sign of fold
  change) using bFlip parameter in dba.report

- Feature: Flip which group in contrast shows gain/loss (sign of fold
  change) using bFlip parameter in dba.plotMA

[diffHic](https://bioconductor.org/packages/diffHic)
-------

Changes in version 1.5.6:

- Relaxed checks in preparePairs(), prepPseudoPairs() when more
  chromosomes are supplied than needed.

- Relaxed checks in connectCounts() when extra chromosomes are in the
  input ranges.

- Fixed an issue with preparePairs() when too many file handles are
  open.

- Fixed clash between BiocGenerics and Matrix which().

- Modified domainDirections() to return a RangedSummarizedExperiment.

- Removed defunct DIList class and methods.

- Switched from seqlevels() to seqlevelsInUse() for fragment intervals.

- Updated user's guide, documentation.

[Director](https://bioconductor.org/packages/Director)
--------

Changes in version 0.99.6:

- Fixed example of function initSankey.

Changes in version 0.99.5:

- Fixes to compilation to pass BiocCheck.

Changes in version 0.99.4:

- Changed R package dependency to htmltools, resulting in modifications
  to initSankey(), drawSankey(), and writeSankey() functions.

Changes in version 0.99.3:

- Modified colour assignment rules to be more intuitively scaled.

- Added additional makeSankey parameters to allow for path colouring
  with a defined threshold/cutoff.

- Updated examples to use package data.

Changes in version 0.99.2:

- Added 'ovca' dataset

Changes in version 0.99.1:

- Remove a package dependency

Changes in version 0.99.0:

- First submission to Bioc-devel

[DOSE](https://bioconductor.org/packages/DOSE)
----

Changes in version 2.11.12:

- deprecate summary method and use as.data.frame instead <2016-09-29,
  Thu>

Changes in version 2.11.11:

- export geneID and geneInCategory <2016-09-19, Mon>

Changes in version 2.11.10:

- geneID and geneInCategory accessor functions <2016-09-19, Mon>

Changes in version 2.11.9:

- enrichMap works with result contains only one term <2016-08-15, Wed>
  + https://github.com/GuangchuangYu/DOSE/issues/15

- build_Anno now works with tibble <2016-08-18, Tue>

Changes in version 2.11.8:

- add unit test <2016-08-15, Mon>

Changes in version 2.11.7:

- change for meshes packages <2016-08-11, Thu>

Changes in version 2.11.6:

- user can use options(DOSE_workers = x) to set using x cores for GSEA
  analysis <2016-08-02, Tue>

- support DisGeNET enrichment analyses <2016-08-01, Mon> + enrichDGN,
  enrichDGNv, gseDGN

- update vignettes <2016-07-29, Fri>

Changes in version 2.11.5:

- enrichMap now output igraph object that can be viewed using other
  software like networkD3 <2016-07-25, Mon>

- dim methods for enrichResult and gseaResult <2016-07-22, Fri>

- $ methods for enrichResult and gseaResult <2016-07-20, Wed>

- switch from parallel to BiocParallel <2017-07-07, Thu>

- [, head and tail methods for enrichResult and gseaResult <2016-07-06,
  Wed>

- change according to GOSemSim <2016-07-05, Tue>

Changes in version 2.11.4:

- 'title' parameter for gseaplot <2016-07-04, Mon> + contributed by
  https://github.com/pedrostrusso +
  https://github.com/GuangchuangYu/DOSE/pull/13

- 'by' parameter in GSEA_internal, by default by = 'fgsea' <2016-07-04,
  Mon> + by = 'fgsea', use GSEA algorithm implemented in fgsea + by =
  'DOSE', use GSEA algorithm implemented in DOSE

- leading edge analysis for GSEA <2016-07-04, Mon>

Changes in version 2.11.3:

- output igraph object in cnetplot <2016-06-21, Tue>

- upsetplot generics <2016-06-14, Tue>

- [[ methods for enrichResult and gseaResult for accessing gene set
  <2016-06-14, Tue>

Changes in version 2.11.2:

- use byte compiler <2016-05-18, Wed>

- https://github.com/Bioconductor-mirror/DOSE/commit/6c508c6a6816f465bb372f30f4ab99c839d81767

Changes in version 2.11.1:

- https://github.com/Bioconductor-mirror/DOSE/commit/7e87d01e671ce1b5fbe974c06b796b1a2970f11c

[EBImage](https://bioconductor.org/packages/EBImage)
-------

Changes in version 4.16.0:

SIGNIFICANT USER-VISIBLE CHANGES

- made defunct deprecated '...GreyScale' family morphological
  functions; use common functions 'dilate', 'erode', 'opening',
  'closing', 'whiteTopHat', 'blackTopHat' and 'selfComplementaryTopHat'
  for filtering both binary and grayscale images

- removed defunct 'getNumberOfFrames' function

PERFORMANCE IMPROVEMENTS

- 'readImage': use 'vapply' instead of 'abind' to reduce memory
  footprint

[edgeR](https://bioconductor.org/packages/edgeR)
-----

Changes in version 3.16.0:

- estimateDisp() now respects weights in calculating the APLs.

- Added design matrix to the output of estimateDisp().

- glmFit() constructs design matrix, if design=NULL, from
  y$samples$group.

- New argument 'null' in glmTreat(), and a change in how p-values are
  calculated by default.

- Modified the default 'main' in plotMD().

- Created a new S3 class, compressedMatrix, to store offsets and
  weights efficiently.

- Added the makeCompressedMatrix() function to make a compressedMatrix
  object.

- Switched storage of offsets in DGEGLM objects to use the
  compressedMatrix class.

- Added the addPriorCount() function for adding prior counts.

- Modified spliceVariants() calculation of the average log-CPM.

- Migrated some internal calculations and checks to C++ for greater
  efficiency.

[EGSEA](https://bioconductor.org/packages/EGSEA)
-----

Changes in version 1.1.10:

- Added: two slots to the GSCollectionIndex: version and date

- Added: citations of the base methods to the documentation of
  egsea.base()

- Removed: rdata.dir from buildIdx functions

- Added: the gene set collection version/update date to the
  GSCollectionIndex class

- Added: an argument to egsea() and egsea.cnt() to return the analysis
  of limma results, which is keep.imma

- Added: an argument to egsea() and egsea.cnt() to return the set
  scores of ssgsea, keep.set.scores

- Added: a slot to the EGSEAResults, which is limmaResults

- Added: a slot to the EGSEAResults, baseInfo

- Added: limmaTopTable, getlimmaResults and getSetScores to the class
  EGSEAResults

- Added: plotSummaryHeatmap to the class EGSEAResults

- Improved: documentation across several functions.

Changes in version 1.1.9 (2016-08-19):

- Fixed: a bug in buildIdx of mouse H gene set

Changes in version 1.1.8 (2016-07-12):

- Removed: EGSEAResults of IL13 from EGSEA and moved it to EGSEAdata

- Updated: EGSEAdata object names in idxAnno

Changes in version 1.1.7 (2016-06-30):

- Improved: the documentation of the methods in the vignette

- Improved: the interpretation of the results in the vignette

- Improved: the ranking when ties occur

- Added: useDingbats = FALSE to pdf() when generating summary plots

- Added: S4 class named GSCollectionIndex to store indexed gene set
  collections

- Added: showSetByName() and showSetByID() to EGSEAResults and
  GSCollectionIndex

- Added: plotGOGraph() to EGSEAResults

- Updated: the GO graphs page of the comparative analysis

- Added: GO graphs to the GO collection of the GeneSetDB

- Fixed: minor bugs

Changes in version 1.1.6 (2016-05-31):

- Fixed: a minor bug in the calculation of the comparative analysis
  p-value.

Changes in version 1.1.5 (2016-05-27):

- Fixed: a minor bug in EGSEAResults when symbolsMap = NULL

- Added: NA Gene Symbols are replaced with Feature IDs in the
  symbolsMap

- Added: FRY to egsea.base()

- Added: several sanity checks on the input parameters

Changes in version 1.1.4 (2016-05-23):

- Improved: ORA to adapt a cut-off threshold logFC=0 if no DE genes
  were found at logFC=1. In both cases, the cut-off threshold of
  adjusted p-value = 0.05.

- Improved: the robustness of the package and allows for single GSE
  analysis to be carried out using EGSEA.

- Added: Ensemble mode is disabled if one base GSE method is provided.

- Added: S4 class for the egsea() output, named EGSEAResults.

- Added: generic methods: show(), summary(), plotHeatmap(),
  plotPathway(), plotMDS(), and plotSummary().

Changes in version 1.1.1 (2016-05-10):

- Improved: topSets(...) and the functionality of the "report" argument
  in the egsea(...) function.

- Improved: verbosity and "print.base" usability in egsea(). The
  statistics of individual methods can now be exported in the output of
  egsea when print.base = TRUE.

- Added: the fry(...) gene set test from the limma package.

- Added: multiple methods to combine the p-values of multiple methods.
  See egsea.combine().

- Fixed: various minor bugs.

[ENmix](https://bioconductor.org/packages/ENmix)
-----

Changes in version 1.9.4:

- bug fix

Changes in version 1.9.3:

- added function relic

Changes in version 1.9.2:

- added function ctrlsva

Changes in version 1.9.1:

- added function oxBS.MLE

[EnrichedHeatmap](https://bioconductor.org/packages/EnrichedHeatmap)
---------------

Changes in version 1.3.5:

- anno_enriched(): add `abs_mean` and `abs_sum` for `value` argument

Changes in version 1.3.4:

- add a new section in the vignette ("Use your own matrix")

- standard error is used instead of standard deviation for the
  annotation chagnes in version 1.3.2

- makeMatrix(): consider when no target is overlaped to any signal.

[EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser)
-----------------

Changes in version 2.3.2:

- Additional sbea methods: gsa, mgsa, padog, globaltest, roast, camera,
  gsva

- Additional nbea methods: netgsa, degraph, topologygsa, ganpa, cepa

[ensembldb](https://bioconductor.org/packages/ensembldb)
---------

Changes in version 1.5.14:

NEW FEATURES

- listEnsDbs function to list EnsDb databases in a MySQL server.

- EnsDb constructor function allows to directly connect to a EnsDb
  database in a MySQL server.

- useMySQL compares the creation date between database and SQLite
  version and proposes to update database if different.

Changes in version 1.5.13:

NEW FEATURES

- useMySQL method to insert the data into a MySQL database and switch
  backend from SQLite to MySQL.

Changes in version 1.5.12:

USER VISIBLE CHANGES

- Add additional indices on newly created database which improves
  performance considerably.

BUG FIXES

- Fix issue #11: performance problems with RSQLite 1.0.9011. Ordering
  for cdsBy, transcriptsBy, UTRs by is performed in R and not in SQL.

- Fix ordering bug: results were sorted by columns in alphabetical
  order (e.g. if order.by = "seq_name, gene_seq_start" was provided
  they were sorted by gene_seq_start and then by seq_name

Changes in version 1.5.11:

BUG FIXES

- makeEnsemblSQLiteFromTables and ensDbFromGRanges perform sanity
  checks on the input tables.

Changes in version 1.5.10:

USER VISIBLE CHANGES

- Using html_document2 style for the vignette.

Changes in version 1.5.9:

NEW FEATURES

- New SymbolFilter.

- returnFilterColumns method to enable/disable that filter columns are
  also returned by the methods (which is the default).

- select method support for SYMBOL keys, columns and filter.

- Select method does ensure result ordering matches the input keys if a
  single filter or only keys are provided.

Changes in version 1.5.8:

BUG FIXES

- Fix problem with white space separated species name in
  ensDbFromGRanges.

Changes in version 1.5.7:

OTHER CHANGES

- Fixed typos in documentation

Changes in version 1.5.6:

BUG FIXES

- Fix warning fo validation of numeric BasicFilter.

Changes in version 1.5.5:

BUG FIXES

- exonsBy: did always return tx_id, even if not present in columns
  argument.

Changes in version 1.5.4:

Bug fixes

- Column tx_id was always removed from exonsBy result even if in the
  columns argument.

- exon_idx was of type character if database generated from a GTF file.

Changes in version 1.5.2:

NEW FEATURES

- Added support for column tx_name in all methods and in the keys and
  select methods. Values in the returned tx_name columns correspond to
  the tx_id.

- Update documentation.

Changes in version 1.5.1:

BUG FIXES

- tx_id was removed from metadata columns in txBy.

- Fixed a bug that caused exon_idx column to be character if database
  created from a GTF.

[ensemblVEP](https://bioconductor.org/packages/ensemblVEP)
----------

Changes in version 1.14.0:

NEW FEATURES

- add support for Ensembl release 85 and 86

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

Changes in version 1.7.2:

SIGNIFICANT USER-VISIBLE CHANGES

- 

BUG FIXES AND MINOR IMPROVEMENTS

- Fixed minor syntax issues in Rd files

- Added more importFrom calls to Namespace file

Changes in version 1.7.1:

SIGNIFICANT USER-VISIBLE CHANGES

- 

BUG FIXES AND MINOR IMPROVEMENTS

- Fixed conditional statement for normalized microarray data in
  testDEArray.

- Added dev.off to close the last plot in an attempt to prevent
  Rplots.pdf when run via Rscripts on commandline.

[esetVis](https://bioconductor.org/packages/esetVis)
-------

Changes in version 0.99.9:

- update documentation

Changes in version 0.99.8:

- fix issue for ggplot2 plot when variable names contain space

Changes in version 0.99.7:

- fix issue: no annotated samples/genes and returnAnalysis set to TRUE

Changes in version 0.99.6:

- fix issue when returnAnalysis and no top elements

Changes in version 0.99.5:

- annotated one gene/sample when topGenes/topSamples is set to 1

Changes in version 0.99.4:

- upgrade package version

Changes in version 0.99.3:

- support object of class 'SummarizedExperiment'

Changes in version 0.99.2:

- increase version number after subscribed to bioc-devel list

Changes in version 0.99.1:

- add extra imports (packages base R) + vignette with clean = FALSE

Changes in version 0.99.0:

- first version of the package submitted to Bioconductor

[ExperimentHub](https://bioconductor.org/packages/ExperimentHub)
-------------

Version: 1.0.0
Category: Package added to release BioC 3.4
Text:

[ExperimentHubData](https://bioconductor.org/packages/ExperimentHubData)
-----------------

Version: 1.0.0
Category: Package added to release BioC 3.4
Text:

[fgsea](https://bioconductor.org/packages/fgsea)
-----

Changes in version 0.99.8:

- Results are now reproducible with set.seed()

Changes in version 0.99.7:

- Performance improvement: slightly better sampling and better segment
  tree

Changes in version 0.99.6:

- Fixed bug with failing when zero pathways are analyzed

Changes in version 0.99.5:

- Using `BiocParallel` package instead of `parallel`

- Implemented leading edge analysis (`leadingEdge` column of fgsea
  result)

[FindMyFriends](https://bioconductor.org/packages/FindMyFriends)
-------------

Changes in version 1.3.3:

- Use homegrown Rcpp based graph class and algorithm for clique
  extraction

[FitHiC](https://bioconductor.org/packages/FitHiC)
------

Version: 0.99.0
Text:

Version: 0.99.1
Text: grammar in FitHiC.R

Version: 0.99.2
Text: benjamini_hochberg_correction in R

Version: 0.99.3
Text:

Version: 0.99.4
Text:

Version: 0.99.5
Text: vignettes

[flowCL](https://bioconductor.org/packages/flowCL)
------

Version: 1.11.1
Text:

Version: 1.11.2
Text:

[flowDensity](https://bioconductor.org/packages/flowDensity)
-----------

Version: 1.7.1
Text:

[flowPloidy](https://bioconductor.org/packages/flowPloidy)
----------

Changes in version 0.99.3 (2016-10-17):

Internal Changes

- Updated DESCRIPTION to include URL and BugReports fields, mostly as
  an exercise to test SVN access and syncing SVN <-> git

Changes in version 0.99.2 (2016-10-12):

Internal Changes

- added a new slot to ModelComponents, `paramLimits`, which allows
  lower and upper limits to be set for each model parameter. (corrects
  bug where linearity wanders below 1, giving nonsense results).

- rationalized the bounds of the data fit in the NLS procedure. Model
  fitting, and RCS calculation, are now all tied to the bin identified
  by fhStart. This selects the highest intensity (peak) in the first 20
  non-zero channels, and ignores all channels below this point. Prior
  to this, the number of observations and associated degrees of freedom
  was calculated in an ad-hoc manner, making the RCS values (even-more)
  difficult to interpret; in addition, the single- and multi- cut
  values started at one channel, but the RCS calculations started on
  another channel, which didn't make sense.

Changes in version 0.99.1 (2016-10-11):

User Visible Changes

- Vignette updated to include installation instructions for
  BioConductor

Internal Changes

- Added accessor functions for ModelComponent and FlowHist classes, so
  direct access of slots via the `@` operator is no longer used outside
  of the initialization functions.

- replaced some loops with vectorized calculations

- replaced call to `eval` with a normal function call to `nlsLM` in
  flowAnalyze.R.

- formatted NEWS file

Changes in version 0.99.0 (2016-08-25):

- What's NEW? Everything so far!

[flowQB](https://bioconductor.org/packages/flowQB)
------

Changes in version 2.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- All new code and new functions, properly documented and covered by
  unit tests.

Changes in version 1.21.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Example files and data have been moved to the flowQBData package.

NEW FEATURES

- Work in progress on functionality that allows users to calculate Q
  and B of an instrument.

Changes in version 1.20.1:

SIGNIFICANT USER-VISIBLE CHANGES

- This is the beginning of major refactoring.  The existing functions
  have been removed and the functionality will be re-implemented. The
  current version is not useful to any users, sorry.

NEW FEATURES

- None, all old features removed for now.

[FunChIP](https://bioconductor.org/packages/FunChIP)
-------

Version: 0.99.0
Text:

Version: 0.99.3
Text:

Version: 0.99.4
Text:

[gCMAP](https://bioconductor.org/packages/gCMAP)
-----

Changes in version 1.17.3:

- BUGFIX: Simpified the mergeCMAPs function.

- BUGFIX: Added additional imports to NAMESPACE, as recommended by R
  CMD check.

Changes in version 1.17.2:

- BUGFIX: Added import of the Matrix package to NAMESPACE file.

Changes in version 1.17.0:

- BUGFIX: Added camera's newly introduced inter.gene.cor parameter in
  camera_score-methods.

[gCMAPWeb](https://bioconductor.org/packages/gCMAPWeb)
--------

Changes in version 1.13.3:

- BUGFIXUpdated NAMESPACE by adding imports recommended by R CMD check.

Changes in version 1.13.1:

- UPDATEAdded note about incompatibility with RStudio.

[gCrisprTools](https://bioconductor.org/packages/gCrisprTools)
------------

Version: 0.99.0
Category: Initial bioconductor submission.
Text:

[gdsfmt](https://bioconductor.org/packages/gdsfmt)
------

Changes in version 1.10.0:

- the version number was bumped for the Bioconductor release version
  3.4

Changes in version 1.8.0-1.8.3:

- the version number was bumped for the Bioconductor release version
  3.3

- define C MACRO 'COREARRAY_ATTR_PACKED' and
  'COREARRAY_SIMD_ATTR_ALIGN' in CoreDEF.h

- SIMD optimization for 1-bit and 2-bit array encode/decode (e.g.,
  decode, RAW output: +20% for 2-bit, +50% for 1-bit)

[geecc](https://bioconductor.org/packages/geecc)
-----

Changes in version 1.7.7 (2016-09-19):

- major revision of init and runConCub resulting in significant
  reduction of memory allocation and runtime

- speed-up of set-operations

- corrected typos in some messages

[GEM](https://bioconductor.org/packages/GEM)
---

Changes in version 0.99.3 (2016-05-24):

MODIFICATION

- add savePlot to GEM_GxEmodel and GEM_Emodel

- modify vignette to plot the figures from codes

Changes in version 0.99.2 (2016-05-21):

MODIFICATION

- change F to FALSE, T to TRUE in codes

- add unit test to inst

- modify man page for SlicedData-class

- modify NEWS to correct format

Changes in version 0.99.1 (2016-05-16):

MODIFICATION

- modify documentations to pass the R CMD check

[genbankr](https://bioconductor.org/packages/genbankr)
--------

Changes in version 1.1.5:

MINOR CHANGES

- Change vignette to use system file to retrieve sample file in
  response to ROpenSci reviewer comment.

- Change vignette builder from knitr to rmarkdown in response to
  ROpenSci reviewer comment.

Changes in version 1.1.4:

BUGFIXES

- Fix bug when mRNA features are present and not annototated with gene.
  Addresses half of https://github.com/gmbecker/genbankr/issues/1 Will
  be backported to 1.0.4

- Fix bug when joins include incomplete range (< or >).  Addresses
  second half of https://github.com/gmbecker/genbankr/issues/1 Will be
  backported to 1.0.4

Changes in version 1.1.3:

BUGFIXES

- Fix bug in .seqTypeFromLOCUS when accession has non-alphanumeric
  characters (will be backported to version 1.0.3.

[GeneNetworkBuilder](https://bioconductor.org/packages/GeneNetworkBuilder)
------------------

Changes in version 1.15.2:

- Fix the bug for error message when missing miRNAlist in
  filterNetwork.

[geneplast](https://bioconductor.org/packages/geneplast)
---------

Changes in version 1.0.0:

- 1st Bioconductor release of geneplast [2015-?-?].

[GeneticsPed](https://bioconductor.org/packages/GeneticsPed)
-----------

Changes in version 2007-09-12:

- Added prune function for pruning/trimming the pedigree. It does not
  work on pedigree, but assumes a data.frame with defined structure.
  Adaption needed.

- We depend on genetics package since gpLong2Wide and hwp assume that
  input is of genotype class.

- Added two utility functions (gpLong2Wide and hwp) for work with
  gpi(). There is also a separate help page for them.

- Internal fixes in gpi - no need to transpose inputs for Fortran call
  anymore - check that there are no NA values in gp and hwp

Changes in version 2007-04-25:

- Version bump to follow BioC releases. # 0.1.4

Changes in version 2007-04-19:

- Added internal .get* functions for retrieving slot values and use of
  it.

Changes in version 2007-04-18:

- R CMD check should now fail also when R error (usually call to
  stop()) occurs in unit testing.

- Added unit tests for sort, examples in help page and clarified sort
  help page.

Changes in version 2007-04-06:

- New small dataset Falconer5.1.

- Added sex(), sex<-(), ascendantSex() and ascendantSex<-() functions.

- Added some more tests for relationshipAdditive, inverseAdditive and
  inbreeding. # 0.1.2

Changes in version 2007-04-01:

- Now we are more rigorous for value of ascendantSex argument in
  Pedigree(). It must accord to values in sex column, if that one is
  passed of course.

- MASS added to depends due to use of fractions() in many places in
  documnetation.

- Added vignette on quantitative genetic (animal) model and
  model.matrix.Pedigree() functions for educational purposes.

- Created data directory and added pedigree example Mrode2.1 and
  Mrode3.1.

- Added unit tests for genetic relationship matrices that should test
  all related functions - mainly against Mrode's book examples -
  runit.genRelMatrix.R.

- Added arguments sort and names to inverseAdditive(),
  relationshipAdditive() and inbreeding().

- Reworked core for relationshipAdditive() into vectorized form.

- Added geneFlowT(), geneFlowTinv(), geneFlowM() and
  mendelianSamplingD() functions.

Changes in version 2007-04:

- Implemented sort by pedigree information only. # 0.1.3

Changes in version 2007-03-02:

- Registration of native routines - src/register.cc.

- Added gpi() function. # 0.1.1

- Temporarily removed unknown funcs, due to planned move to BioC and
  change to S4.

- All depends are now in gdata --> removing ggmisc.

Changes in version 2006-03-29:

- codeUnit to get internal codes for subject and ascendants if factors
  are used.

- Handled unused levels in nlevels.pedigree and summary.pedigree now
  produces a simple summary

- Started to implement checks in check.pedigree and friends

- Proper factor handling

Changes in version 2006-03-16:

- Playing around with NA/unknown representation - it is very likely
  that I messed up some things with factors --> subject to changes. #
  Version 0.1

- Initial version # NEWS ends here

[genomation](https://bioconductor.org/packages/genomation)
----------

Changes in version 1.5.6:

NEW FUNCTIONS AND FEATURES

- annotateWithFeatures() function added for a diversity of annotation
  operations bases on genomic interval overlap.

- heatTargetAnnotation() returns a heatmap for percentage of genomic
  intervals overlapping with annotation features. These functions are
  more general than the old annotation functions and could be used in
  variety of settings.

- plotGeneAnnotation() is depracated, use heatTargetAnnotation()
  instead. It does the same job and has more features.

Changes in version 1.5.5:

IMPROVEMENTS AND BUG FIXES

- fixed warnings related to "[" function of ScoreMatrix and
  ScoreMatrixList

- fixed issue "scoreMatrixList drops names after subsetting #141"

- fixed issue "heatMeta returns a list not a matrix #142"

- fixed vignette issue, now things do not depend on knitrBootstrap #143

Changes in version 1.5.4:

IMPROVEMENTS AND BUG FIXES

- If strand.aware=TRUE and the GRanges object provided to
  ScoreMatrix/ScoreMatrixBin for windows was not sorted by coordinates,
  then the results were distorted.

Changes in version 1.5.3:

IMPROVEMENTS AND BUG FIXES

- Count bam faster in cpp. Fast retrieval of library size in c++ using
  Samtools idxstats-like code implemented by alexg9010.

Changes in version 1.5.2:

IMPROVEMENTS AND BUG FIXES

- merged updates from BioC SVN repository: - elementLengths was renamed
  -> elementNROWS in S4Vectors (new name reflects TRUE) - use new
  invertStrand() instead of GenomicAlignments:::invertRleStrand() -
  bump version prior to creation of 3.3 branch, bump version of all
  packages that use knitr for vignettes, and bump version to trigger
  package rebuilding now that purl()'ing issue.

[GenomeInfoDb](https://bioconductor.org/packages/GenomeInfoDb)
------------

Changes in version 1.10.0:

NEW FEATURES

- Add function mapGenomeBuilds() that maps between UCSC and Ensembl
  builds.

- Add function genomeBuilds() that list all the available UCSC or
  Ensembl builds for a given organism[s] that can be used in
  mapGenomeBuilds()

- Add listOrganism() that list all currently available organism[s]
  included for use in genomeBuilds()

DEPRECATED AND DEFUNCT

- After being deprecated, the species() method for GenomeDescription
  objects is now defunct

MODIFICATIONS

- Zebra finch is removed as option for fetchExtendedChromInfoFromUCSC()
  as it is not support yet

- keepStandardChromosomes() chooses first style when multiple are
  matched

BUG FIXES

- Fix WARNING occuring when determining style in
  keepStandardChromosomes()

[GenomicAlignments](https://bioconductor.org/packages/GenomicAlignments)
-----------------

Changes in version 1.10.0:

NEW FEATURES

- The GAlignmentPairs container now supports pairs with discordant
  strand and/or seqnames. The "granges" and "ranges" methods for
  GAlignmentPairs objects get new argument 'on.discordant.seqnames' to
  let the user control how to handle pairs with discordant seqnames.
  See ?GAlignmentPairs for more information.

- Add "invertStrand" method for GAlignmentPairs objects.

- Add 'use.names' argument to the "ranges", "granges", "grglist" and
  "rglist" methods for GAlignments and GAlignmentsList objects.

- Add 'use.names' argument to the "granges" and "grglist" methods for
  GAlignmentPairs objects.

- Add "ranges" method for GAlignmentPairs objects.

SIGNIFICANT USER-LEVEL CHANGES

- The 'at' argument of pileLettersAt() is now expected to be a GPos
  object (GRanges still accepted).

- 50x speed-up of the granges() extractor for GAlignmentPairs object.
  The improvement is based on a suggestion by Arne Muller.

DEPRECATED AND DEFUNCT

- Remove left() and right() generics and methods (were defunct in BioC
  3.3).

- Remove 'invert.strand' argument from "first" and "last" methods for
  GAlignmentPairs objects (was defunct in BioC 3.3).

- Remove strand() setter for GAlignmentPairs objects (was defunct in
  BioC 3.3).

- Remove 'order.as.in.query' argument from "grglist" method for
  GAlignmentPairs objects and from "grglist" and "rglist" methods for
  GAlignmentsList objects (was defunct in BioC 3.3).

BUG FIXES

- Fix 'use.names=FALSE' in "grglist" and "rglist" methods for
  GAlignmentsList objects.

[GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures)
---------------

Version: 1.26
Category: NEW FEATURES
Text: makeTxDbFromGRanges() now recognizes features of type lnc_RNA,
        antisense_lncRNA, transcript_region, and pseudogenic_tRNA, as
        transcripts.

Version: 1.26
Category: NEW FEATURES
Text: Add 'intronJunctions' argument to mapToTranscripts().

Version: 1.26
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text:

Version: 1.26
Category: DEPRECATED AND DEFUNCT
Text: The 'vals' argument of the "transcripts", "exons", "cds", and
        "genes" methods for TxDb objects is now defunct (was deprecated
        in BioC 3.3).

Version: 1.26
Category: DEPRECATED AND DEFUNCT
Text: The "species" method for TxDb object is now defunct (was
        deprecated in BioC 3.3).

Version: 1.26
Category: BUG FIXES
Text:

[GenomicRanges](https://bioconductor.org/packages/GenomicRanges)
-------------

Version: 1.26.0
Category: NEW FEATURES
Text: Add 'with.revmap' argument to "reduce" method for GRangesList
        objects.

Version: 1.26.0
Category: NEW FEATURES
Text: Add 'with.revmap' argument to various "disjoin" methods.

Version: 1.26.0
Category: NEW FEATURES
Text: makeGRangesFromDataFrame() now tries to turn the "start" and
        "end" columns of the input data frame into numeric vectors if
        they are not already.

Version: 1.26.0
Category: NEW FEATURES
Text: Add makeGRangesListFromDataFrame() function.

Version: 1.26.0
Category: NEW FEATURES
Text: Add "summary" method for GenomicRanges objects.

Version: 1.26.0
Category: NEW FEATURES
Text: Add 'use.names' argument to the granges(), grglist(), and
        rglist() generics and methods, as well as to a bunch of
        "ranges" methods (for GRanges, GPos, GNCList, GRangesList, and
        DelegatingGenomicRanges). Default is TRUE to preserve existing
        behavior.

Version: 1.26.0
Category: NEW FEATURES
Text: Add 'use.mcols' arguments to the "ranges" methods for GPos
        objects.

Version: 1.26.0
Category: SIGNIFICANT USER-LEVEL CHANGES
Text:

Version: 1.26.0
Category: DEPRECATED AND DEFUNCT
Text:

Version: 1.26.0
Category: BUG FIXES
Text: Fix bug in distanceToNearest() related to ranges starting at
        zero.

Version: 1.26.0
Category: BUG FIXES
Text: Fix GRanges(Seqinfo()).

[GenRank](https://bioconductor.org/packages/GenRank)
-------

Version: 1.0.3
Category: Made slight modifications to the vignette and added a couple
        of figures. Changed the way RankProduct method handles the
        missing evidence across evidence layers
Text:

[GenVisR](https://bioconductor.org/packages/GenVisR)
-------

Changes in version 1.1.5:

- added clarification for use of cnFreq (cnFreq requires consistent
  windows across samples!)

- minor documentation changes

- bug fix for lohView VAF values are no longer assumed to be from 0-100
  and may range from 0-1

[ggtree](https://bioconductor.org/packages/ggtree)
------

Changes in version 1.5.17:

- read.nhx support newick file <2016-10-17, Mon> +
  https://github.com/GuangchuangYu/ggtree/issues/79

Changes in version 1.5.16:

- read.phyloT for parsing newick format of phyloT output <2016-10-11,
  Tue> + https://www.biostars.org/p/210401/#216128

- fixed aes mapping in geom_strip <2016-10-11, Tue>

- fixed R check <2016-10-10, Mon> + check.aes parameter is not
  available in release version of ggplot2 yet

Changes in version 1.5.15:

- check.aes for layers defined in ggtree <2016-10-07, Fri>

- recalculate 'angle' when collapse, expand and rotate clade
  <2016-10-06, Thu> + https://github.com/GuangchuangYu/ggtree/issues/78

Changes in version 1.5.14:

- subset tip in geom_tiplab2 <2016-10-05, Wed>

- add `compute_group` according to ggplot (v2.1.0) <2016-09-29, Thu> +
  https://github.com/hadley/ggplot2/issues/1797

- unit test for groupOTU and groupClade <2016-09-22, Thu>

- groupOTU label groups by input group names (when input is a named
  list) <2016-09-22, Thu>

- update angle calculation for geom_tiplab <2016-09-13, Thu>

- as.polytomy to collapse binary tree to polytomy by applying 'fun' to
  selected 'feature' (e.g. bootstrap value less than 70). <2016-09-13,
  Tue> + currently only phylo object supported. + add test for
  as.polytomy

Changes in version 1.5.13:

- facet_plot for plotting data with tree <2016-09-06, Tue>

- more parameters for column names in gheatmap <2016-09-06, Tue> +
  colnames_angle + colnames_offset_x + colnames_offset_y + hjust

- offset parameter in geom_tiplab and geom_tiplab2 <2016-09-05, Mon>

Changes in version 1.5.12:

- use data in all layers instead of the base layer for coordination
  calculation in subview <2016-09-01, Thu>

- bug fixed in subview, width & height should be width/2 & height/2
  <2016-09-01, Thu>

Changes in version 1.5.11:

- gheatmap works with matrix <2016-08-28, Sun> +
  https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/bioc-ggtree/2YLvXHMJJ6U/c4zS7yfGCAAJ

- support parsing expression in geom_strip <2016-08-18, Thu>

- bug fixed in geom_tiplab <2016-08-17, Wed> +
  https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/bioc-ggtree/Tm9ULK7hd9E/HviXEh3CBwAJ

- update citation info, add doi. <2016-08-16, Tue>

Changes in version 1.5.10:

- fixed issue #72 for label of geom_treescale not displayed
  <2016-08-16, Tue> + https://github.com/GuangchuangYu/ggtree/issues/72

Changes in version 1.5.9:

- update citation info <2016-08-12, Fri>

Changes in version 1.5.8:

- add color parameter in geom_cladelabel, color should be of length 1
  or 2 <2016-08-11, Thu>

- geom_cladelabel support parsing expression <2016-08-11, Thu>

Changes in version 1.5.7:

- geom_strip can accept taxa name as input but labeling strip will not
  supported. To support labeling strip, user need to input node id
  <2016-07-27, Wed>

- nodeid function for converting node label(s) to node id(s)
  <2016-07-27, Wed>

Changes in version 1.5.6:

- remove dependency of Biostring for installing ggtree <2016-07-21,
  Thu> + still needed for building vignette and for processing FASTA
  file

- remove dependency of EBImage for building & installing ggtree
  <2016-07-21, Thu> + the package is still needed if user want to
  annotate tree with image file

- `%<+%` now works with tbl_df <2016-07-21, Thu> +
  https://github.com/GuangchuangYu/ggtree/issues/66

- identify method for ggtree <2016-06-28, Tue> + see
  https://guangchuangyu.github.io/2016/06/identify-method-for-ggtree

- geom_balance contributed by Justin Silverman <2016-06-22, Wed> + see
  https://github.com/GuangchuangYu/ggtree/pull/64

Changes in version 1.5.5:

- update geom_tiplab2 according to angle change introduced by open_tree
  <2016-06-20, Mon>

- bug fixed in collapse, now work with collapse a clade that contain a
  subclade that was already collapsed <2016-06-02-Thu>

- bug fixed if time-scaled tree extend into the BCE. <2016-06-02, Thu>
  + as.Date won't work for BCE time. + as.Date=FALSE by default in
  fortify method, just use the time in decimal format (real number, not
  Date object).

Changes in version 1.5.4:

- reroot method for raxml object <2016-05-22, Sun>

- bug fixed in scaleClade, now y positions are (hopefully) always
  correct. <2016-05-20, Fri>

- bug fixed in collapse <2016-05-20, Fri> + if user collapse a node
  that is an offspring of a collapsed node, print warning msg and
  return the tree directly

- use byte compiler <2016-05-18, Wed>

- change any(is.na()) to anyNA() which is more efficient <2016-05-18,
  Wed>

- https://github.com/Bioconductor-mirror/ggtree/commit/559548c66b51253e8ccb983d353385838a81f106

Changes in version 1.5.3:

- add examples in vignettes <2016-05-13, Fri> + add fan layout example
  in treeVisualization vignette + add open_tree and rotate_tree example
  in treeManipulation vignette

- add angle in ggtree function, fan layout supported <2016-05-12, Thu>

- rotate_tree and open_tree function <2016-05-12, Thu>

- support reading BEAST MCC trees (multiple trees in one file) via the
  read.beast function <2016-05-12, Thu>

- https://github.com/Bioconductor-mirror/ggtree/commit/51eec4721595c274c24dc4df2f1fdf40700cb1a5

Changes in version 1.5.2:

- add multiplot in ggtreeUtilities vignette <2016-05-12, Thu>

- add example of integrate user's data using phylo4d in treeAnnotation
  vignette <2016-05-11, Wed>

- add extend, extendto parameter in geom_hilight <2016-05-10, Tue>

- geom_hilight now supports hilight tips <2016-05-10, Tue> +
  https://github.com/GuangchuangYu/ggtree/issues/53

- more accurate ylim & angle for circular layout <2016-05-10, Tue> +
  https://github.com/GuangchuangYu/ggtree/issues/40

- supports phylo4d object <2016-05-10, Tue> +
  https://github.com/GuangchuangYu/ggtree/issues/47

Changes in version 1.5.1:

- update vignettes <2016-05-10, Tue> + add geom_range example in
  treeImport + add geom_strip and geom_taxalink example in
  treeAnnotation + add ggtreeUtilities vignette

- gheatmap now works with data.frame of only one column <2016-05-09,
  Mon> + contributed by Justin Silverman <jsilve24@gmail.com> +
  https://github.com/GuangchuangYu/ggtree/pull/57

- geom_strip for associated taxa <2016-05-09, Mon> +
  https://github.com/GuangchuangYu/ggtree/issues/52

[Glimma](https://bioconductor.org/packages/Glimma)
------

Changes in version 1.2.0:

- Added option to turn off internal cpm transform (transform=FALSE) in
  MD Plot.

- Added gridlines to MD Plot.

- Added enter to search on MD plot table search box.

- Added sample.cols argument to MD Plot.

- Added data.frame handling for groups argument in MDS plot.

- Added option to not use counts argument.

- Added option to not use anno argument.

- Added glXYPloy for more general plotting.

- Added option to use numerical value for "groups" in MD Plot.

- Fixed multiple plots in same directory overwriting each other's data.

- Fixed logical values in annotation breaking graphs.

- Fixed numeric values not working for colours.

Changes in version 1.1.0:

- Added tables to MD Plot.

[globalSeq](https://bioconductor.org/packages/globalSeq)
---------

Changes in version 1.1.0 (2016-10-18):

- Minor improvements

[GOexpress](https://bioconductor.org/packages/GOexpress)
---------

Changes in version 1.7.1:

GENERAL UPDATES

- Updated email address.

[GOpro](https://bioconductor.org/packages/GOpro)
-----

Version: 0.99.2
Text:

[GOSemSim](https://bioconductor.org/packages/GOSemSim)
--------

Changes in version 1.99.4:

- fixed NOTE in R check <2016-08-12, Fri>

- add unit test using testthat <2016-08-11, Thu>

Changes in version 1.99.3:

- changes to satisfy meshsim package <2016-08-05, Fri>

Changes in version 1.99.2:

- fixed Rcpp issue <2016-07-19, Tue> +
  https://github.com/GuangchuangYu/GOSemSim/issues/6

Changes in version 1.99.1:

- update vignette <2016-07-14, Thu>

Changes in version 1.99.0:

- support all organisms that have OrgDb object <2016-07-05, Tue>

- optimize Wang method <2016-07-04, Mon>

Changes in version 1.31.2:

- use byte compiler <2016-05-18, Wed>

- https://github.com/Bioconductor-mirror/GOSemSim/commit/71c29280c560e0293569121aeeecb0ed7b37055a

Changes in version 1.31.1:

- https://github.com/Bioconductor-mirror/GOSemSim/commit/a829a50a017b90f08c41b5955df176dfad333d06

[gprege](https://bioconductor.org/packages/gprege)
------

Changes in version 1.17.1 (2016-10-07):

- Replaced download urls for downloading DellaGattaData.RData while
  building the vignette.

[GSALightning](https://bioconductor.org/packages/GSALightning)
------------

Version: 1.1.1
Category: Codes
Text:

Version: 1.1.1
Category: Added “maxmean” statistics
Text:

Version: 1.1.1
Category: Added “wilcoxTest” for Mann-Whitney-U tests for single gene
        analysis
Text:

Version: 1.1.1
Category: Default for GSALight() is now “maxmean” with
        restandardization
Text:

Version: 1.1.1
Category: Fix some issues to make GSALight() equivalent to GSA of Efron
        and Tibshirani
Text:

Version: 1.1.1
Category: The target gene list is now a list instead of a data table
Text:

Version: 1.1.1
Category: The Vignette now contains a comprehensive user guide
Text:

Version: 1.1.1
Category: All documentations are updated accordingly
Text:

Version: 1.1.2
Category: Added “importFrom("stats", "p.adjust", "rbinom", "var",
        "wilcox.test”)” to NAMESPACE
Text:

Version: 1.1.2
Category: Import “stats” in description
Text:

Version: 1.1.2
Category: Added back manual for “wilcoxTest
Text:

Version: 1.1.3
Category: Bioconductor version updates
Text:

Version: 1.1.5
Category: Updated references and manual and vignette
Text:

Version: 1.1.5
Category: Added Github URL
Text:

Version: 1.1.6
Category: Fixed default settings of permTestLight() by adding
        match.args and setting nperm=NULL
Text:

[GSAR](https://bioconductor.org/packages/GSAR)
----

Changes in version 1.8.0:

- New function MDtest is introduced. It implements a nonparametric
  multivariate test of means based on sample ranking in the MST similar
  to function KStest, but the test statistic is the mean deviation
  between the CDFs of two conditions.

- New function RMDtest is introduced. It implements a nonparametric
  multivariate test of variance based on sample ranking in the MST
  similar to function RKStest, but the test statistic is the mean
  deviation between the CDFs of two conditions.

- New function AggrFtest is introduced. It implements a nonparametric
  test of variance by aggregating the univariate p-values obtained by
  the F-test using Fisher's probability combining method. It test the
  hypothesis that all genes in a gene set show no significant
  difference in variance between two conditions against the alternative
  hypothesis that at least one gene in the gene set shows significant
  difference in variance between two conditions.

- New function findMST2.PPI is introduced. It finds the union of the
  first and second MSTs similar to function findMST2, but it accepts an
  object of class igraph as input rather that a matrix of gene
  expression data. The input igraph object represents a protein-protein
  interaction (PPI) network that can be binary or weighted, directed or
  undirected.

- New wrapper function TestGeneSets is introduced. It performs a
  specific statistical method from the ones available in package GSAR
  for multiple gene sets. The gene sets are provided as a list of
  character vectors where each entry has the feature (gene) identifiers
  in a single gene set.

- New argument pvalue.only added to all available statistical methods
  in the package. When pvalue.only=TRUE (default), each statistical
  method returns the p-value only. When pvalue.only=FALSE, each
  statistical method returns a list of length 3 consisting of the
  observed statistic, vector of permuted statistics, and p-value.

- New arguments leg.x, leg.y, group1.name, group2.name, label.color,
  label.dist, vertex.size, vertex.label.font, and edge.width added to
  function plotMST2.pathway to allow more flexibility in generating
  plots. The values of most of these arguments are passed to function
  plot.igraph.

[GSEABase](https://bioconductor.org/packages/GSEABase)
--------

Changes in version 1.35:

BUG FIXES

- some multi-line warnings and errors would fail without reporting the
  error message.

- getGmt(), GO GeneSetCollection methods much faster.

[Guitar](https://bioconductor.org/packages/Guitar)
------

Changes in version 1.11.9:

- Added "combinedGuitarPlot" function to support comparison of
  different species within a single figure.

- Allow rescale

[GWASTools](https://bioconductor.org/packages/GWASTools)
---------

Changes in version 1.19.1:

- Speed up corr.by.snp in duplicateDiscordance.

[Harman](https://bioconductor.org/packages/Harman)
------

Changes in version 1.2.0:

- Dynamically resizes legends depending upon the number of batches in
  pcaPlot

- A custom prcomp function to get the appropriate scores. The standard
  R prcomp function did not work in instances where the number of
  samples (matrix columns) was greater than the number of assays
  (matrix rows), so a special case is needed for less assays than
  samples. We need to use u' instead of v from the SVD. Presently, this
  is under developement, so an error is thrown if rows < cols.

- prcompPlot now has an argument for scaling, which defaults to FALSE.
  Previously, scaling was always TRUE. This new default makes
  prcompPlots plots agree with the 'original' plots of plotting
  harmanresults objects.

- Extensive updating of the vignette with a new comparison to ComBat
  from sva.

Changes in version 1.0.2:

- First public version on Bioconductor.

[HDTD](https://bioconductor.org/packages/HDTD)
----

Changes in version 1.7.1 (2016-09-15):

- Updated CITATION FILE.

[hiAnnotator](https://bioconductor.org/packages/hiAnnotator)
-----------

Changes in version 1.7.1:

- code spacing edits & namespace conflict fixes

[HIBAG](https://bioconductor.org/packages/HIBAG)
-----

Changes in version 1.10.0:

- the version number was bumped for the Bioconductor release version
  3.4

Changes in version 1.9.0-1.9.3:

- the development version

Changes in version 1.8.0-1.8.3:

- the version number was bumped for the Bioconductor release version
  3.3

- new arguments 'pos.start' and 'pos.end' in `hlaFlankingSNP()`

[hiReadsProcessor](https://bioconductor.org/packages/hiReadsProcessor)
----------------

Changes in version 1.9.3:

- Namespace fixes for base R packages to alleviate build errors

- Updated findOverlap calls to use the new drop.self & drop.redundant
  args

Changes in version 1.9.1:

- Using readxl to import sampleInfo excel file

- Minor code improvements for merging internal objects and cleaning
  data

- Removed parameter interactive in read.SeqFolder in favor of
  base::interactive()

- findOverlap calls updated to use ignoreRedundant & ignoreSelf

- splitSeqsToFiles receives outDir parameter

[HiTC](https://bioconductor.org/packages/HiTC)
----

Changes in version 1.17.1:

NEW FEATURES

- New getPearsonMap function. Will generate the correlation map used by
  the pca analysis

- The pca.hic function is now able to detect and to assign the A/B
  compartment if a gene annotation is provided

SIGNIFICANT USER-VISIBLE CHANGES

- update of reduce method for HTClist object

- update of the getExpectedCounts function with two methods ; loess and
  mean. The mean method allows to estimated the expected counts using
  the mean of diagonal matrices. The method is adviced in case of high
  resolution maps when the loess smoothing can take time and might not
  give good results

- By default, obs/exp maps are now centered before calculated the
  pearson correlation map. This allow correlation of small values to be
  as valuable as correlation of big values

- asRangedData from importC is now deprecated due to rtracklayer change

- normPerExpected function now reports 0 instead of NA to avoid error
  in correlation calculation

BUG FIXES

- Fix bug in directionalityIndex. The input count matrix is converted
  into a dense matrix before calculation

[hpar](https://bioconductor.org/packages/hpar)
----

Changes in version 1.15.3:

- Updating to HPA version 15, with new datasets <2016-09-14 Wed>

Changes in version 1.15.2:

- Version bump to trigger package rebuilding now that purl()'ing issue
  has been correctly identified. knitr does not create purl()'ed
  (Stangle equivalent) .R files if _R_CHECK_TIMINGS_ is set, which the
  build system was setting. Now it's not set, so these .R files are now
  created. See https://github.com/yihui/knitr/issues/1212 for more.
  r117512 | d.tenenbaum | 2016-05-15 21:14:22 +0100 (Sun, 15 May 2016)
  | 6 lines

Changes in version 1.15.1:

- Bump version of all packages that use knitr for vignettes. This is
  because of an issue (now fixed) in knitr which failed to create
  purl()'ed R files from vignette sources and include them in the
  package. This version bump will cause these packages to propagate
  with those R files included. r117081 | d.tenenbaum | 2016-05-03
  22:30:44 +0100 (Tue, 03 May 2016) | 2 lines

Changes in version 1.15.0:

- Bioc version 3.3

[HTSFilter](https://bioconductor.org/packages/HTSFilter)
---------

Changes in version 1.13.1:

- -- Fixed bug in vignette causing error message on dev (related to new
  estimateDisp function in edgeR package) -- All documentation now
  automatically done through Roxygen -- Functionality for DESeq
  pipeline (CountDataSet objects) now deprecated -- Add possibility for
  parallel computing via BiocParallel -- Recorrect multiple testing
  correction of p-values in DESeq2 framework (only correct p-values for
  genes that pass filter) using pAdjustMethod argument in HTSFilter and
  HTSBasicFilter functions -- Started adding unit testing framework

[illuminaio](https://bioconductor.org/packages/illuminaio)
----------

Changes in version 0.15.1 (2016-08-27):

- Now the package DLL is unloaded when the package is unloaded.

Changes in version 0.15.0 (2015-05-03):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.4 for R (>= 3.3.0).

[InPAS](https://bioconductor.org/packages/InPAS)
-----

Changes in version 1.5.7:

NEW FEATURE

- decrease the memory to run InPAS.

Changes in version 1.5.6:

NEW FEATURE

- Add a parameter to control search the proximal site from both
  direction or not.

Changes in version 1.5.5:

NEW FEATURE

- search the proximal site from both direction.

Changes in version 1.5.4:

NEW FEATURE

- Categorize novel distal events into "novel distal" and "truncated
  distal"

BUG FIXES

- fix the bug in proximalAdjByCleanUpdTSeq.

Changes in version 1.5.1:

NEW FEATURE

- new function coverageRate for quality control.

BUG FIXES

- fix the error when not all the seqnames is in bedgraphs for CPsites.

[InteractionSet](https://bioconductor.org/packages/InteractionSet)
--------------

Changes in version 1.1.6:

- Allowed specification of NULL in row/column arguments to inflate().

- Fixes to tests and code in response to updates to BiocGenerics,
  S4Vectors.

- Added CITATION to the F1000Res article.

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

Changes in version 1.6:

- Added support of CAMERA isotope identification (findIsotopes.CAMERA)

- selectivity of findIsotopes.IPO may be increased if
  checkBorderIntensity is set to TRUE: 'maxo' of each peak of an
  isotopologue must be three times higher than the intensities at
  'rtmin' and 'rtmax'

- simplified return value of calcPPS() to vector with meaningful names

- changed getDefaultXcmsSetStartingParams() for min_peakdwidth = c(12,
  28) and for ppm to c(17, 32)

Changes in version 1.5.7.1:

- using predict() to identify best levels and expand.grid to generate
  testdata for model

Changes in version 1.5.7:

- supporting single parameter optimization. Only basic version with
  redundant levels in consecutive DoEs

- removed integer-rounding in maximum focusing for all findPeaks
  parameters except prefilter(I) and steps

- Update documentation to match code

- Remove use of getwd() preventing absolute subdir paths

Changes in version 1.5.6:

- updated vignette

- generally using Central-Composite design instead of Box-Behnken
  design

- fixed bug when defining subdir=NULL in functions optimizeXcmsSet and
  optimizeRetGroup

- modified unit tests to handle versions > 1.5.6

- added function writeRScript to NAMESPACE export

- updated man for optimizeXcmsSet and optimizeRetGroup

Changes in version 1.5.5:

- increased recall of reliable peaks in calcPPS. Changed exponent for
  PPS calculation from 1.5 to 2

Changes in version 1.5.4.8:

- remove file lookup in optimizeXcmsSet, and leave that to xcmsSet()

Changes in version 1.5.4.7:

- explicitely use serial evaluation if nSlaves=1, to help debugging

Changes in version 1.5.4.6:

- remove dependency on Rmpi

Changes in version 1.5.4.5:

- added examples from msdata

- fix Depends, imports and library() and require() calls

Changes in version 1.5.4.4:

USER VISIBLE CHANGES

- packaged script

- changed method name attachparams to attachList

- changed method name calculateRGTV to getRGTVValues

- changed method name getDefaultStartingXcmsParams to
  getDefaultXcmsSetStartingParams

- changed method name typeCastFactor to typeCastParams

- changed method name writeRSkript to writeRScript

- changed the parameter name n_slaves to nSlaves

- resultIncreased: if last optimization score was 0, no isotopes have
  been found hence the dataset is not optimizable with IPO.

- added man files for attachList, calcPPS, combineParams, createModel,
  decode, decodeAll, encode, getBbdParameter, getCcdParameter,
  getDefaultRetCorCenterSample, getDefaultRetGroupStartingParams,
  getDefaultXcmsSetStartingParams, getNormalizedResponse,
  getRGTVValues, IPO-package, optimizeRetGroup, optimizeXcmsSet,
  startSlaves, toMatrix, typeCastParams

- removed xcmsSetsettingsAsString.R

- getResponses: now able to handle NULL value for slices parameter

- calcPPS: peaks with NA values are removed before isotopes
  identification

- if subdir is NULL, no rsm's are saved

- optimizeXcmsSet: lowere minimum value for min_peakwith from 5 to 3
  #IPO_V1.5.4.3: * LIP calculation in calcPPS fixed #IPO_V1.5.4.2: *
  added initial parameter check # * renamed all factor-variables to
  params # * in optimizeXcmsSet: - also look for mzML-files # - check
  if files were found # * bug in optimization for matchedFilter fixed;
  sigma and mzdiff have to be # definded later (combineFactors()) when
  sigma and step as well as steps are already known #IPO_V1.5.4.1:
  changes in calcPPS: # rt_window <- rt * 0.005 # rt_lower <-
  part_peaks[,"rt"] - rt_window # rt_upper <- part_peaks[,"rt"] +
  rt_window #IPO_V1.5.4: if bad_group == 0; bad_group = 1 && good_group
  += 1 #IPO_V1.5.3: no parameter for isotope detection.  #
  c13_peak[,"mz"] has to be within (mzmin + isotope_mass) and (mzmax +
  isotope_mass) # c13_peak[,"rt"] has to be within (rtmin +
  isotope_mass) and (rtmax + isotope_mass) #IPO_V1.5.: in
  RCSandGSIncreased: also used good_groups ^ 2 #IPO_V1.4.: vectorized
  isotope identification; # no intensity window, between intensity of
  max carbon and 1 #IPO_V1.3.: good_groups ^ 2 to increase recall

[IRanges](https://bioconductor.org/packages/IRanges)
-------

Changes in version 2.8.0:

NEW FEATURES

- "disjoin" methods now support 'with.revmap' argument.

- Add 'invert' argument to subsetByOverlaps(), like grep()'s invert.

- Add "unstrsplit" method for RleList objects.

- findOverlapPairs() allows 'subject' to be missing for self pairing.

- Add "union", "intersect" and "setdiff" methods for Pairs.

- Add distance,Pairs,missing method.

- Add ManyToManyGrouping, with coercion targets from FactorList and
  DataFrame.

- Add Hits->List and Hits->(ManyToMany)Grouping coercions.

- Add "as.matrix" method for AtomicList objects.

- Add "selfmatch", "duplicated", "order", "rank", and "median" methods
  for CompressedAtomicList objects.

- Add "anyNA" method for CompressedAtomicList objects that ensures
  recursive=FALSE.

- Add "mean" method for CompressedRleList objects.

- Support 'global' argument on "which.min" and "which.max" methods for
  CompressedAtomicList objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Make mstack,Vector method more consistent with stack,List method.

- Optimize and document coercion from AtomicList to RleViews objects.

DEPRECATED AND DEFUNCT

- Are now defunct (were deprecated in BioC 3.3): - RangedDataList
  objects.  - RDApplyParams objects and rdapply().  - The "split" and
  "reduce" methods for RangedData objects.  - The 'ignoreSelf' and/or
  'ignoreRedundant' arguments of the findOverlaps,Vector,missing method
  (a.k.a. "self findOverlaps" method).  - grouplength() - GappedRanges
  objects.

BUG FIXES

- Fix special meaning of findOverlaps's maxgap argument when
  type="within".

- isDisjoint(IRangesList()) now returns logical(0) instead of NULL.

- Fixes to regroup() and Grouping construction.

- Fix rank,CompressedAtomicList method.

- Fix fromLast=TRUE for duplicated,CompressedAtomicList method.

[isobar](https://bioconductor.org/packages/isobar)
------

Changes in version 1.19.1:

- Fix error in subsetIBSpectra and exclude and namespace imports
  (Laurent Gatto <2016-10-07 Fri>)

[isomiRs](https://bioconductor.org/packages/isomiRs)
-------

Changes in version 1.1.5:

OTHERS

- Suppress message from readr when reading

- Add instruction of how to create DESeq2 object from ids one.

Changes in version 1.1.4:

OTHERS

- Add mRNA-miRNA integration

Changes in version 1.1.2:

OTHERS

- Update version dependecy for packages

[IVAS](https://bioconductor.org/packages/IVAS)
----

Changes in version 1.5.1:

- Debug error when result is NULL in the MsqtlFinder.

- Adjust Chr names between GTF and SNP locus data.

- Change test of UTR region.

- Debug error in Boxplot, calSignificant and sqtlfinder.

[JunctionSeq](https://bioconductor.org/packages/JunctionSeq)
-----------

Changes in version 1.3.4:

Minor bugfixes

- Manual updated to reflect changes in the most recent version of
  QoRTs: the mergeNovelSplices function now automatically calculates
  size factors, if size factors are not set explicitly.

- Fixed typos in the manual.

- Bugfix: Basic functionality now works without a GFF file being
  specified (as per specification)

- Iterated version number to match Bioconductor versioning.

- Added the ability to suppress the creation of a DESeq2 count
  container object (for advanced users).

[kebabs](https://bioconductor.org/packages/kebabs)
------

Changes in version 1.7.3:

- fix of citation info file

- fix of inst/NEWS

Changes in version 1.7.2:

- version number bump for technical reasons related to Bioconductor
  build servers

Changes in version 1.7.1:

- version number bump for technical reasons related to Bioconductor
  build servers

Changes in version 1.7.0:

- new branch for Bioconductor 3.4 devel

[kimod](https://bioconductor.org/packages/kimod)
-----

Version: 1.1.1
Text:

[LedPred](https://bioconductor.org/packages/LedPred)
-------

Version: 1.7.4
Text: Deals with error when at least one feature in some training fold
        of the ModelPerformance CV algo has all value zero -> removes
        it

Version: 1.7.3
Text: Fixed opposite sign of feature weigths in some situations

Version: 1.7.3
Text: LedPred function returns object with many more informations

Version: 1.7.1
Category: Hotfix: Fixed sign of prediction
Text:

[limma](https://bioconductor.org/packages/limma)
-----

Changes in version 3.30.0:

- New function cumOverlap() to analyse the overlap between two ordered
  lists.

- New function detectionPValues() to compute detection p-values from
  negative control probes.

- New function fitmixture() to estimate genewise fold changes and
  expression values from mixture experiments.  Previously this funtion
  was only available as part of the illumina package available from
  http://bioinf.wehi.edu.au/illumina

- New function logcosh() to compute log(cosh(x)) accurately without
  floating underflow or overflow.

- The default settings for the 'inter.gene.cor' and 'use.neg.cor'
  arguments of camera() have been changed.  camera() now uses by
  default a preset value for the inter-gene correlation. This has the
  effect that it tends to rank highly co-regulated,
  biologically-interpretable more highly than before.

- New flexibility for the roast() and mroast() functions, similar to
  that previously implemented for fry().  The index vector for each
  gene set can now be a data.frame, allowing each gene set to have its
  own set of gene weights.  The indices can now optionally be a vector
  of gene names instead of a vector of indices.  roast() and mroast()
  now support the robust empirical Bayes option of squeezeVar().
  roast() and mroast() can now accept, via '...', any argument that
  would be normally passed to lmFit() or eBayes().

- Slight change to the standardize="posterior.sd" method for fry().

- goana(), alias2Symbol() and alias2SymbolTable() now work for any
  species for which an Entrez Gene based organism package exists.

- The 'species' argument of kegga() can now accept any Bioconductor
  species abbreviation.

- topGO() now breaks ties by the number of genes in the GO Term and by
  the name of the Term if the p-values are equal. (This is the same
  behavior as topKEGG.)

- New argument 'plot' for plotMDS(), to optionally allow an MDS object
  to be returned without making a plot.

- New arguments 'annotation' and 'verbose' for read.idat().  The first
  of these to allows users to read any required columns from the
  manifest file.

- New arguments 'pch.y' and 'pch.z' for plotRLDF.

- Unnecessary argument 'design' removed from fitted.MArrayLM.

- normalizeBetweenArrays() now checks whether the input 'object' is a
  data.frame, and converts to a matrix if possible.

- duplicateCorrelation() now expands weights using expandAsWeights(),
  making it consistent with lmFit().

- The lowess line drawn by plotSA() is now more robust with respect to
  NA variances.

- More informative error message from voom() when there is only one row
  of data.

- More informative error message from getEAWP() and lmFit() when
  'object' is a non-normalized data object.

- Update Phipson et al (2016) reference for robust empirical Bayes.

- Bug fix to fitFDistRobustly(), which in some circumstances was trying
  to save extra (undocumented) results that had not been computed.

- Bug fixes for fry() with robust=TRUE or when 'index' has NULL names.

- Bug fix to propTrueNull() when method="hist" and all the p-values are
  less than 1/nbins.

- Bug fix to alias2Symbol() with expand.symbol=TRUE.

[LINC](https://bioconductor.org/packages/LINC)
----

Changes in version 1.1.6:

NEW FEATURES

- There is now a link to a GitHub repository included in the package
  vignette The example data for 'justlinc' can be downloaded from this
  website.

SIGNIFICANT USER-LEVEL CHANGES

- Co-expressed genes can now be returned as gene symbols or other
  keytypes using the function 'getcoexpr'.

BUG FIXES

- Standard values of arguments in 'singlelinc' are now consistent
  across the definitions of generic function and methods.

[LOBSTAHS](https://bioconductor.org/packages/LOBSTAHS)
--------

Changes in version 2016-04-21:

- Initial release for Bioconductor

[M3D](https://bioconductor.org/packages/M3D)
---

Version: 1.7.4
Text:

Version: 1.7.5
Text:

Version: 1.7.11
Text:

Version: 1.7.12
Text:

[maftools](https://bioconductor.org/packages/maftools)
--------

Changes in version 0.99.50:

NON SIGNIFICANT CHANGES

- Bux fix in annovartomaf.

- Rainfall plot improvement - Better handling of empty contigs.

- Better parsing of HGVSp annotations.

- Default value for useAll in read.maf changed from FALSE to TRUE.

- FLAG notes

- Added CITATION file.

Changes in version 0.99.45:

NON SIGNIFICANT CHANGES

- Bux fixes in oncoplot while drawing rowbar.

- Bux fixes in subsetMaf.

- Better error messages in read.maf.

Changes in version 0.99.40:

NON SIGNIFICANT CHANGES

- Bux fixes in plotMAFsummary.

- Support new version of ComplexHeatmap.

Changes in version 0.99.34:

SIGNIFICANT USER-LEVEL CHANGES

- rainfall plot can detect katageis when detectChangePoints is TRUE.

NON SIGNIFICANT CHANGES

- Bux fixes.

Changes in version 0.99.33:

NEW FEATURES

- icgcSimpleMutationToMAF - convert ICGC Simpale Somatic Mutation
  Format to MAF

SIGNIFICANT USER-LEVEL CHANGES

- Now read.maf validates MAF for repeated variants mapping to multiple
  transcripts of same gene.

NON SIGNIFICANT CHANGES

- Fixed vignette TOC

Changes in version 0.99.32:

NON SIGNIFICANT CHANGES

- Fixed warnings for bioc-check.

Changes in version 0.99.31:

NEW FEATURES

- CoOncoplot - Plot two oncoplots side by side

- rainfallPlot

- oncoplot now supports copy number data

- New function readGistic to read summarize output files from GISTIC.

- dashboard style plot for mafsummary

- Plotting functions for GISTIC results: plotGisticResults, gisticPlot

- geneCloud to plot wordcloud of frequently mutated genes and copy
  number alteartions.

SIGNIFICANT USER-LEVEL CHANGES

- More plotting parameters to plotCBSsegments.

BUG FIXES

- Bug fix while performing fisher exact for differences between
  cohorts.

DEPRECATED AND DEFUNCT

- Remove AddReadcounts function, which relied on extrenal patform
  specific parogram bam-readcount

Changes in version 0.99.30:

NEW FEATURES

- Plotting and mapping variants on copy number data (require cbs
  segments)

- Compare two maf files for differentially mutated genes

- Forestplots

- Better clustering with copynumber and outlier detection for
  heterogeneity.

Changes in version 0.99.25:

NEW FEATURES

- inferHeterogeneity: Clustering using nonparametric dirichlet process;
  more plotting options.

- New subsetMaf function to subset MAF

SIGNIFICANT USER-LEVEL CHANGES

- Fixed a warning in oncodrive function

- Upto 5X faster oncomatrix creation and sorting using data.table

- Examples won't take >5 seconds

- Updated vignette

Changes in version 0.99.20:

CHANGES

- Resubmission to bioconductor using r-devel 3.3

Changes in version 0.99.15:

CHANGES

- Resubmission to bioconductor using r-devel 3.4

- added plotBestFitRes argument to extractSignatures

- added DriverMutation to biocviews

Changes in version 0.99.1:

CHANGES

- Modified vignettes to fix build vignettes error

[MANOR](https://bioconductor.org/packages/MANOR)
-----

Changes in version 2014-10-01:

- Removed unused C functions for I/O in 'lib_io'.

- Minor changes to C code to avoid WARNINGS upon R CMD check.

- Changed 'par' settings in vignette to avoid 'Figure margins too
  large' errors.

Changes in version 2013-07-03 (2013-07-03):

- Replaced calls to 'exit', 'fprintf' which now raise a WARNING upon
  check.

- Iteration index not printed anymore in 'nem' (raised an error when
  compiling the vignette).

Changes in version 2011-03-18 (2011-03-18):

- Updated calls to 'GLAD:::daglad' in the vignette to fix an error
  caused by changes in the defaults of 'daglad'.

Changes in version 2010-10-01 (2010-10-01):

- Cleaned 'data/flags.RData" which contained objects from an old
  'globalenv'.

- Updated maintainer's email address.

Changes in version 2010-01-24 (2010-01-24):

- Added (back) 'intensity.flag' to 'data/flags.RData' (had been removed
  since v 1.12.0 for an unknown reason).

Changes in version 2009-01-15 (2009-01-15):

- misplaced alignment tab in man/spatial.Rd.

Changes in version 2009-01-13 (2009-01-13):

- updated references in .Rd files.

- fixed warnings due to incorrect use of \item in .Rd files.

Changes in version 2009-01-06 (2009-01-06):

- 'norm' and 'sort' are now S3 methods as well

Changes in version 2009-01-04 (2009-01-04):

- (almost) one file per function in R/

- removed empty section \details in man/qscore.Rd

- added a NAMESPACE

- removed inst/doc/Makefile (not needed anymore because no html output
  required)

Changes in version 2009-01-02 (2009-01-02):

- removed another non-standard keyword

Changes in version 2009-01-01 (2009-01-01):

- only one keyword per \keyword entry...

Changes in version 2008-12-31 (2008-12-31):

- now use standard "keyword"s

- changed \link{\code{stuff}} into \code{\link{stuff}}

Changes in version 2008-11-26 (2008-11-26):

- filled in "keyword" sections in .Rd files.

- removed empty "examples" sections from .Rd files.

- initialized a few variables upon declaration in C code to prevent
  warnings in R CMD CHECK.

Changes in version 2008-09-23 (2008-09-23):

- modification de la fonction cv pour retourner NA lorsque toutes la
  valeurs du vecteur sont <e0> NA

- modification de la function getChromosomeArm pour que cytoband ne
  soit pas positionn<e9>e <e0> NULL

Changes in version 2008-09-04 (2008-09-04):

- added a CHANGELOG

- updated outdated reference in the .bib file

- changed the definition of flag "rep.flag" to avoid the error now
  caused by sd(NA, na.rm=TRUE)

[matter](https://bioconductor.org/packages/matter)
------

Changes in version 0.99.11 (2016-10-11):

BUG FIXES

- Updated documentation details.

Changes in version 0.99.10 (2016-10-11):

SIGNIFICANT USER-VISIBLE CHANGES

- Updated vignettes with installation instructions, faster build time

- Added 'adata' method for accessing 'matter' class @data slot

Changes in version 0.99.9:

NEW FEATURES

- Added delayed scaling and centering via 'scale' method

- Added 'prcomp' method for principal components analysis

SIGNIFICANT USER-VISIBLE CHANGES

- Renamed 'colSd' -> 'colSds', 'colVar' -> 'colVars', etc.

- Renamed 'filepaths' -> 'paths' and 'file_id' -> 'source_id'

- Moved 'irlba' from Suggests to Imports to support new 'prcomp' method

- Updated vignette to use new 'prcomp' method in the PCA example

BUG FIXES

- Fixed bug when combining 'matter' objects with multiple data sources

Changes in version 0.99.8:

SIGNIFICANT USER-VISIBLE CHANGES

- In S4 class 'atoms', slot 'file_id' is now type 'integer' to save
  space

- In S4 class 'atoms', slot 'datamode' is now type 'integer' to save
  space

- More comprehensive error messages in constructors for S4 classes

Changes in version 0.99.7:

BUG FIXES

- Try to fix namespace std::isnan scoping issues on Windows

Changes in version 0.99.6:

BUG FIXES

- Fixed handling of NA, NaN, Inf, and -Inf during C type coercion

- Improved handling of NA, NaN, Inf, and -Inf in summary stats

- Fixed handling of NAs in matrix multiplication for integers

Changes in version 0.99.5:

BUG FIXES

- Fixed .Call native routine registration for C++ methods

- Added "C_" prefix for C++ methods called through .Call

Changes in version 0.99.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Added .Call native routine registration for C++ methods

Changes in version 0.99.3:

BUG FIXES

- Import generics from S4Vectors for 'colMeans', 'colSums', 'rowMeans',
  and 'rowSums'

- Cleaned up method signatures and class unions

Changes in version 0.99.2:

BUG FIXES

- Cleaned up double assignments in C++ code

Changes in version 0.99.1:

BUG FIXES

- Version bump for Bioconductor build system

Changes in version 0.99.0:

BUG FIXES

- Updated PCA example in vignette (irlba now requires 'mult' argument
  to be non-missing for non-C execution)

- Added irlba unit test

Changes in version 0.6:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated maintainer and author email address

BUG FIXES

- Fixed bug in rowVar for matter_matc

Changes in version 0.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated documentation with examples and added unit tests

BUG FIXES

- Fixed matrix multiplication on mixed data types (int x double)

Changes in version 0.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Added support for class-preserving subsetting of 'matter' matrices
  with 'drop=NA' argument for Cardinal compatibility.

- Added new 'show' method for 'matter' vectors and matrices showing
  their size in in-memory and size on disk.

Changes in version 0.3:

NEW FEATURES

- Added C++ class 'MatterAccessor' for iterating through.  a buffered
  version of a 'matter' vector or matrix

- Added summary statistics including 'sum', 'mean', 'var', 'sd',
  'colSums', 'colMeans', 'colVar', 'colSd', 'rowSums', 'rowMeans',
  'rowVar', and 'rowSd'.

- Added support for 'apply' method for 'matter' matrices

- Added support for 'bigglm' linear regression.

- Added basic matrix multiplication for 'matter' matrices with an
  in-memory R matrix or vector.

Changes in version 0.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Overhauled backend to use C++ classes 'Matter' and 'Atoms', to
  maximizes use of sequential reads versus random reads.

Changes in version 0.1:

SIGNIFICANT USER-VISIBLE CHANGES

- First rough implementation of matter including the classes 'atoms'
  and 'matter', and subclasses 'matter_vec', 'matter_matc',
  'matter_matr', with a C backend.

[meshes](https://bioconductor.org/packages/meshes)
------

Changes in version 0.99.7:

- Depends DOSE <2016-09-29, Thu>

Changes in version 0.99.6:

- update docs <2016-09-19, Mon>

Changes in version 0.99.3:

- remove non-ASCII character in vignettes/meshes.bib <2016-08-16, Tue>

Changes in version 0.99.2:

- all exported functions have runnable examples required by BiocCheck
  <2016-08-16, Tue>

Changes in version 0.99.1:

- R check passed <2016-08-16, Tue>

Changes in version 0.99.0:

- add vignette <2016-08-11, Thu>

- enrichMeSH and gseMeSH <2016-08-11, Thu> + move from clusterProfiler

Changes in version 0.0.1:

- geneSim <2016-08-10, Wed>

- meshSim <2016-08-05, Fri>

[MeSHSim](https://bioconductor.org/packages/MeSHSim)
-------

Version: 2015-01-12
Text:

Version: 2015-01-12
Text:

Version: 2015-01-12
Text:

Version: 2015-01-12
Text:

Version: 2015-01-12
Text:

Version: 2015-01-11
Text:

Version: 2015-01-08
Text:

Version: 2014-12-14
Text:

Version: 2014-11-16
Text:

[metabomxtr](https://bioconductor.org/packages/metabomxtr)
----------

Changes in version 1.7.1:

- *modified package functions to handle situations where metabolites
  are entirely missing for a particular level of a categorical
  predictor in the mixture model.

[MetaboSignal](https://bioconductor.org/packages/MetaboSignal)
------------

Changes in version 1.0.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Package introduced.

NEW FEATURES

- Package introduced.

[metagenomeFeatures](https://bioconductor.org/packages/metagenomeFeatures)
------------------

Changes in version 1.1.3 (2016-10-01):

- replaced msd16S with mock community dataset for examples

- changed tree mgFeatures and MgDb slot from class phylo to phyloOrNULL

Changes in version 1.1.2 (2016-03-26):

- Added mgFeatures class - this class replaces metagenomeAnnotation in
  the 16S workflow, and contains database information for a user
  provided list of database sequence ids. A new
  metagenomeAnnotation-class will be added to the package when a
  suitable R native 16S taxonomic classificaiton method is available.

- new `aggregate_taxa` function for aggregating MRexperiment objects to
  user defined taxonomic level, aggretation of count data is performed
  using sums by dafault, users can pass any column wise matrix
  operation.

[metagenomeSeq](https://bioconductor.org/packages/metagenomeSeq)
-------------

Changes in version 1.15:

- Added 'mergeMRexperiment' function

- Added 'normFactors' and 'libSize' generics

- Added 'fitMultipleTimeSeries' function

- Replaced RUnit with testthat library for unit testing

- Adding multiple upgrades and changes throughout

- Deprecated the load_* functions and created load* function.

[MetCirc](https://bioconductor.org/packages/MetCirc)
-------

Changes in version 0.99.6:

- improve allocation of clicked features in shinyCircos [2016-10-10
  Mon]

- change line width and colour in plotCircos [2016-10-10 Mon]

Changes in version 0.99.5:

- add arguments colour and transparency for functions plotCircos and
  highlight [2016-09-15 Thu]

Changes in version 0.99.4:

- add arguments splitIndMZ and splitIndRT in function convert2MSP and
  fix bug (colnames classes instead of class) in convert2MSP
  [2016-08-17 Wed]

- update vignette, unit tests and manuals [2016-08-17 Wed]

Changes in version 0.99.3:

- remove verbatim Textouput(help) in shinyApp [2016-08-01 Mon]

- change function createSimilarityMatrix in order that it does not
  change ordering of row names [2016-08-01 Mon]

- update unit tests and manuals [2016-08-01 Mon]

Changes in version 0.99.2:

- do not reorder again in createOrderedSimMat function [2016-07-26 Tue]

- change allocatePrecursor2mz function such, that it is compatible with
  all sd01 and sd02 objects [2016-07-26 Tue]

Changes in version 0.99.1:

- rewrite functions that they do not require the argument dfNameGroup
  (data.frame containing group and unique identifier) any longer
  [2016-06-14 Tue]

- the function binning uses now the function cut to create bins binning
  has two methods implemented to calculate from these bins [2016-06-10
  Fri] median and mean m/z values from the fragment m/z values
  [2016-06-10 Fri]

- documentation about data sets is extended [2016-06-10 Fri]

- include 'Suggests: BiocGenerics' in the DESCRIPTION file [2016-06-10
  Fri]

Changes in version 0.99.0:

- submit package to Bioconductor file tracker [2016-05-14 Sat]

Changes in version 0.98.0:

- add information of hovered objects in shinyCircos [2016-04-16 Tue]

- allow for subsetting of MSP objects [2016-04-12 Tue]

- add unit tests for exported functions [2016-04-12 Tue]

[methylKit](https://bioconductor.org/packages/methylKit)
---------

Version: 0.99.4
Category: IMPROVEMENTS AND BUG FIXES
Text: all annotation functions are now removed and are available
        through genomation package. See the vignette for examples. Most
        functions have similar names and functionality in genomation.
        The users need to convert methylKit objects to GRanges using
        as(methylKit.obj,"GRanges") before they can use genomation
        functions.

Version: 0.99.3
Category: IMPROVEMENTS AND BUG FIXES
Text: Fixed a bug in processBismark() C++ function where BAM files are
        should be treated as 0-based.  Now this is fixed.

Version: 0.99.3
Category: IMPROVEMENTS AND BUG FIXES
Text: Bug in calculateDiffMethDSS is fixed
        (https://github.com/al2na/methylKit/issues/49)

Version: 0.99.2
Category: IMPROVEMENTS AND BUG FIXES
Text: Fixed a bug in methRead() introduced after the addition of mincov
        argument.  the bug occurred only reading files legacy text
        files that have CpGs with coverage below 10

Version: 0.99.2
Category: IMPROVEMENTS AND BUG FIXES
Text: Changes to vignette for better description of the tests.

Version: 0.99.2
Category: IMPROVEMENTS AND BUG FIXES
Text: Compiler error that occurs in older compilers are fixed via this
        PR https://github.com/al2na/methylKit/pull/43

Version: 0.99.1
Category: IMPROVEMENTS AND BUG FIXES
Text: mostly changes to meet BioCcheck() requirements and
        reccomendations

Version: 0.99.1
Category: IMPROVEMENTS AND BUG FIXES
Text: C++ code compiles on windows now, regex requirement is no longer
        there.

Version: 0.9.6
Category: IMPROVEMENTS AND BUG FIXES
Text: changes to following function names: read() to methRead()
        read.bismark to processBismarkAln() adjust.methylC() to
        adjustMethylC() get.methylDiff() to getMethylDiff()
        annotate.WithFeature() to annotateWithFeature()
        annotate.WithFeature.Flank() to annotateWithFeatureFlank()
        annotate.WithGenicParts() to annotateWithGenicParts()
        read.bed() to readBed() read.feature.flank() to
        readFeatureFlank() read.transcript.features() to
        readTranscriptFeatures()

Version: 0.9.6
Category: IMPROVEMENTS AND BUG FIXES
Text: Improved documentation for methRead() (old read())

Version: 0.9.6
Category: IMPROVEMENTS AND BUG FIXES
Text: Now, bismark cytosine report and coverage files can be read using
        methRead() pipeline argument. see ?methRead

Version: 0.9.6
Category: IMPROVEMENTS AND BUG FIXES
Text: Ported the Perl script for methylation base calling to C/C++ via
        Rcpp.  Contributed by Alexander Gosdschan.

Version: 0.9.6
Category: IMPROVEMENTS AND BUG FIXES
Text: methRead() uses data.table::fread() to read files faster.

Version: 0.9.6
Category: IMPROVEMENTS AND BUG FIXES
Text: methRead() has a new argument mincov, which sets the minimum
        number of reads that needs to cover a base. Positions with
        coverage below this number are discarded.

Version: 0.9.6
Category: NEW FUNCTIONS AND FEATURES
Text: new function methSeg() can segment methylation (methylRaw
        objects) and differential methylation (methylDiff objects)
        profiles to segments.  Associated function methSeg2bed()
        creates BED files from segments.  see ?methSeg. A test is added
        to check this in R CMD check.  Contributed by Arsene Wabo and
        Alexander Gosdschan.

Version: 0.9.6
Category: NEW FUNCTIONS AND FEATURES
Text: new tabix based classes methylRawDB, methylRawListDB,
        methylBaseDB, methylDiffDB and respective methods implemented.
        Tests are updated to check proper function in R CMD check.
        Contributed by Alexander Gosdschan.

Version: 0.9.6
Category: NEW FUNCTIONS AND FEATURES
Text: calculateDiffMeth() now supports basic overdispersion correction
        and multiple methods for pvalue correction. The function also
        now handles covariates such as age,sex etc. A test is added to
        check this in R CMD check.  Contributed by Adrian Bierling.

Version: 0.9.6
Category: NEW FUNCTIONS AND FEATURES
Text: New function calculateDiffMethDSS() is using beta-binomial model
        from DSS package to calculate differential methylation.
        Contributed by Dhruva Chandramohan with modifications from
        Altuna Akalin. This is a modified version of the function from
        DSS package so that it can work with methylKit objects.

Version: 0.9.6
Category: NEW FUNCTIONS AND FEATURES
Text: dataSim creates a methylBase object with simulated methylation
        data. A test is added to check this in R CMD check. Contributed
        by Adrian Bierling.

Version: 0.9.5
Category: IMPROVEMENTS AND BUG FIXES
Text: travis CI build shield added

Version: 0.9.4
Category: IMPROVEMENTS AND BUG FIXES
Text: tileMethylCounts now works on methylBase objects, affected by
        BioC 3.0 upgrade.  A test is added to check this in R CMD
        check.

Version: 0.9.3
Category: IMPROVEMENTS AND BUG FIXES
Text: compatibility with BioC 3.0. Multiple BioC functions moved to
        other packages within BioC,which broke some code and caused
        installation issues.  Now this is fixed.

Version: 0.9.3
Category: IMPROVEMENTS AND BUG FIXES
Text: data.table::merge is now stable when all=TRUE. Removed my(altuna)
        version of data.table::merge the code base.

Version: 0.9.2.5
Category: IMPROVEMENTS AND BUG FIXES
Text: calculateDiffMeth slim=FALSE argument works correctly when there
        are multiple samples per group.

Version: 0.9.2.4
Category: IMPROVEMENTS AND BUG FIXES
Text: install_github() now works correctly. Removed blank lines at the
        end of DESCRIPTION file.

Version: 0.9.2.4
Category: IMPROVEMENTS AND BUG FIXES
Text: calculateDiffMeth and tileMethylCounts now works correctly, typos
        were introduced in code with 0.9.2.2. Now these are fixed

Version: 0.9.2.3
Category: IMPROVEMENTS AND BUG FIXES
Text: regionCounts() bug occurring when the first argument is a class
        of methylBase is fixed.  The bug introduced an additional
        column to the resulting methylBase object.
        https://groups.google.com/forum/#!topic/methylkit_discussion/p19K-pgavAI

Version: 0.9.2.2
Category: IMPROVEMENTS AND BUG FIXES
Text: unite() issues when destrand=FALSE is resolved. The issue
        appeared due to numeric vs.  integer conflict when merging data
        sets using chr,start and end locations.

Version: 0.9.2.2
Category: IMPROVEMENTS AND BUG FIXES
Text: typo in calculateDiffMeth() argument "weighted.mean" is fixed.

Version: 0.9.2.2
Category: IMPROVEMENTS AND BUG FIXES
Text: unite() is faster when destrand=TRUE, due to improvements in
        internal function .CpG.dinuc.unify()

Version: 0.9.2.2
Category: IMPROVEMENTS AND BUG FIXES
Text: regionCounts() strand.aware argument now works correctly. This
        bug had no effect on tileMethylCounts. With the default
        arguments, strand aware setting was not in effect at all, so
        every region was treated as strandless.

Version: 0.9.2.2
Category: IMPROVEMENTS AND BUG FIXES
Text: Conversion to GRanges objects now removes seqlevels(chromosome
        names) that are not used. Not removing these could cause
        problems in regionCounts() function and functions depending on
        regionCounts().

Version: 0.9.2.2
Category: IMPROVEMENTS AND BUG FIXES
Text: read() function in some cases was treating strand columns of the
        flat CpG files as logical if strand had F or R letters and most
        of the first CpGs were having F strand (Forward strand). Now
        this is fixed.

Version: 0.9.2.1
Category: IMPROVEMENTS AND BUG FIXES
Text: Fixes a bug in getMethylationStats() function that generates
        incorrect label numbers when the user overrides the default
        number of breaks for the histogram. Contributed by Bonnie
        Barrilleaux.

Version: 0.9.2
Category: IMPROVEMENTS AND BUG FIXES
Text: A bug introduced with 0.9.1 is fixed. The bug occured when
        unite() function

Version: 0.9.2
Category: is used with destrand=TRUE argument. It returned less number
        of CpGs then it was
Text:

Version: 0.9.2
Category: supposed to, although returned CpGs had correct methylation
        and coverage values
Text:

Version: 0.9.1
Category: NEW FUNCTIONS AND FEATURES
Text: objects takes less memory and they are reorganized.
        updateMethObject() updates the objects from previous versions
        the latest version.

Version: 0.9.1
Category: NEW FUNCTIONS AND FEATURES
Text: New batch effect control functions are implemented. You can
        control if certain principal components are associated with
        batch effects and remove those components from your data. see
        ?assocComp amd ?removeComp. In addition, if you have corrected
        for the batch effects via other methods, you can reconstruct a
        corrected methylBase object ( see ?reconstruct). Check the
        "batch effects" section in the vignette.

Version: 0.9.1
Category: IMPROVEMENTS AND BUG FIXES
Text: unite() function takes less time due to use of data.table::merge

Version: 0.9.1
Category: IMPROVEMENTS AND BUG FIXES
Text: fixes a bug appeared with R 3.02, where getting exon and intron
        coordinates from BED12 files produced an error.

Version: 0.9.1
Category: IMPROVEMENTS AND BUG FIXES
Text: tested with R 3.2 and matching bioconductor packages

Version: 0.9.1
Category: IMPROVEMENTS AND BUG FIXES
Text: data(methylKit) has the new objects

Version: 0.5.7
Category: IMPROVEMENTS AND BUG FIXES
Text: tested with R 3.0 and matching bioconductor packages

Version: 0.5.7
Category: IMPROVEMENTS AND BUG FIXES
Text: deprecated matchMatrix import from IRanges was causing a problem
        with package installation. v0.5.7 fixes that problem.

Version: 0.5.7
Category: IMPROVEMENTS AND BUG FIXES
Text: Now there is no need for "chr" string in BED files when reading
        them as annotation files. Some assemblies do not have the
        "chr"" string in their chromosome names.

Version: 0.5.6
Category: NEW FUNCTIONS AND FEATURES
Text: new arguments for clusterSamples() and PCASamples() functions are
        added. With the new options
        "sd.filter","sd.threshold","filterByQuantile" are added. These
        options help finetune how low variation bases/regions are
        discarded prior to clustering or PCA. See ?PCASamples and
        ?clusterSamples() for details on the new options

Version: 0.5.6
Category: NEW FUNCTIONS AND FEATURES
Text: FAQ section added to the vignette

Version: 0.5.6
Category: NEW FUNCTIONS AND FEATURES
Text: show methods added for each class. Now, typing the variable name
        containing the object will display concise information about
        the contents of the object.

Version: 0.5.6
Category: NEW FUNCTIONS AND FEATURES
Text: Subsetting objects via "[" notation is now enabled. You can
        subset rows of the objects and it will return a new object
        rather than just a data frame.

Version: 0.5.6
Category: IMPROVEMENTS AND BUG FIXES
Text: tileMethylCounts() error is fixed. Error occurred when tilling
        sparsely covered small chromosomes like chrM in human RRBS
        data.

Version: 0.5.6
Category: IMPROVEMENTS AND BUG FIXES
Text: read.bismark() can deal with Bismark output with Bowtie2. Bowtie2
        can put gaps in the alignment, now read.bismark() can deal with
        those gaps when parsing the SAM file.

Version: 0.5.6
Category: IMPROVEMENTS AND BUG FIXES
Text: Coverage columns are coerced to integer when reading generic
        methylation per base files. BSMAP scripts can produce a
        methylation ratio file where coverages (or "effective CT
        counts" as they are called) are not always integers, which
        causes a problem in the downstream analysis. Now, these
        non-integer columns are rounded to nearest integer while
        reading.  See
        http://zvfak.blogspot.com/2012/10/how-to-read-bsmap-methylation-ratio.html
        for example usage of this functionality.

Version: 0.5.5
Category: IMPROVEMENTS AND BUG FIXES
Text: Differential methylation percentage calculation bug fixed. The
        bug occurred when "min.per.group" argument used in unite()
        function and when "weighted.mean=TRUE" in calculateDiffMeth()
        function.

Version: 0.5.5
Category: IMPROVEMENTS AND BUG FIXES
Text: plotTargetAnnotation bug is fixed. Bug occured when "precedence"
        argument set to FALSE.

Version: 0.5.4
Category: IMPROVEMENTS AND BUG FIXES
Text: Examples added to help pages

Version: 0.5.4
Category: IMPROVEMENTS AND BUG FIXES
Text: Changes to DESCRIPTION for complying with the bioconductor
        guidelines

Version: 0.5.4
Category: IMPROVEMENTS AND BUG FIXES
Text: unused "cor" option removed from PCASamples() function

Version: 0.5.4
Category: IMPROVEMENTS AND BUG FIXES
Text: some irrelevant functions are not exported (i.e they are not
        public) anymore.

Version: 0.5.4
Category: IMPROVEMENTS AND BUG FIXES
Text: getContext() looks for the correct slot name now

Version: 0.5.3
Category: NEW FUNCTIONS AND FEATURES
Text: new function adjust.methylC() can be used to adjust measured 5mC
        levels by measured 5hmC levels

Version: 0.5.3
Category: IMPROVEMENTS AND BUG FIXES
Text: pool() function bug fixed where one of the groups have one sample
        the pool function was not returning correct values

Version: 0.5.3
Category: IMPROVEMENTS AND BUG FIXES
Text: calculateMethDiff() function option "SLIM" is now working. If set
        to TRUE SLIM method for q-value calculation will be used. If
        FALSE, p.adjust with method="BH" option will be used for
        P-value correction.

Version: 0.5.3
Category: IMPROVEMENTS AND BUG FIXES
Text: read.bismark() function now works correctly under Windows.

Version: 0.5.3
Category: IMPROVEMENTS AND BUG FIXES
Text: read.bismark() bug occurring when the reads are paired and
        overlapping is fixed now.

Version: 0.5.3
Category: IMPROVEMENTS AND BUG FIXES
Text: read.bismark() bug occurring when the alignment is done in
        non-directional manner is fixed now. However, illumina
        sequencing protocols are directional and you are unlikely to
        have encountered this error if you aligned your sequences in a
        directional manner.

Version: 0.5.3
Category: IMPROVEMENTS AND BUG FIXES
Text: getCorrelation() function has a new option called "method" can
        take the value of "spearman","person" or "kendall"

Version: 0.5.2
Category: NEW FUNCTIONS AND FEATURES
Text: new function pool() sums up coverage, numCs and numTs values
        within each group so one representative sample for each group
        will be created in a new methylBase object.

Version: 0.5.2
Category: NEW FUNCTIONS AND FEATURES
Text: new function normalizeCoverage() normalizes coverage and
        associated number of Cs and number of Ts values between samples
        using a scaling factor derived from the ratio between mean or
        median of coverage distributions of samples.

Version: 0.5.1
Category: IMPROVEMENTS AND BUG FIXES
Text: calculateDiffMeth() can now deal with differential methylation
        calculations where one group has multiple samples but the other
        one has only one.

Version: 0.5
Category: NEW FUNCTIONS AND FEATURES
Text: new function reorganize() can be used to reorganize methylRawList
        and methylBase objects to create new objects from subset of
        samples or can be used to change the order of samples and
        treatment vector.

Version: 0.5
Category: NEW FUNCTIONS AND FEATURES
Text: new function bedgraph() can output UCSC bedgraph files for
        methylRaw and methylDiff objects

Version: 0.5
Category: NEW FUNCTIONS AND FEATURES
Text: new function percMethylation() extracts percent methylation
        values from a methylBase object and returns a matrix

Version: 0.5
Category: NEW FUNCTIONS AND FEATURES
Text: unite() now can merge bases/regions that are not covered by all
        samples by setting "min.per.group" option

Version: 0.5
Category: NEW FUNCTIONS AND FEATURES
Text: PCASamples() has new options "transpose" and "sd.threshold". Now,
        one can do PCA analysis on the transposed % methylation values.
        "sd.threshold" is for removing rows with small variation prior
        to PCA analysis .see ?PCASamples for details

Version: 0.5
Category: IMPROVEMENTS AND BUG FIXES
Text: regionCounts() bug that appeared with data.table 1.8.0 is fixed

Version: 0.5
Category: IMPROVEMENTS AND BUG FIXES
Text: unite() bug that appeared with R 2.15 is fixed

Version: 0.5
Category: IMPROVEMENTS AND BUG FIXES
Text: calculateDiffMeth() can deal with NA values introduced by setting
        "min.per.replicate" option of unite() function

Version: 0.5
Category: IMPROVEMENTS AND BUG FIXES
Text: PCASamples() now uses "prcomp" function to do the PCA analysis

Version: 0.5
Category: IMPROVEMENTS AND BUG FIXES
Text: external annotation data (cpgi.hg18.bed.txt and
        refseq.hg18.bed.txt) that comes with the package is now only a
        subset of full datasets. Do not use them for your own analysis,
        download the complete annotation data from UCSC or other
        sources in BED format.

Version: 0.4.1
Category: NEW FUNCTIONS AND FEATURES
Text: new function select() can be used to subset methylRaw, methylBase
        and methylDiff objects to create new objects with a subset of
        methylation data useful if you want to use only a particular
        portion of methylation events on the genome.

Version: 0.4.1
Category: IMPROVEMENTS AND BUG FIXES
Text: read.bismark bug fixed, now it should be able to save methylation
        call files with no problem

Version: 0.4
Category: NEW FUNCTIONS AND FEATURES
Text: New read.bismark() function can read directly from a sorted SAM
        file output by Bismark aligner, the function can save
        methylation calls as text files or read them as methylRaw
        object to be used in analysis

Version: 0.4
Category: NEW FUNCTIONS AND FEATURES
Text: calculateDiffMeth() can do all differential methylation
        calculations using multiple-cores if 'num.cores' is set to an
        integer denoting how many cores should be used.

Version: 0.4
Category: NEW FUNCTIONS AND FEATURES
Text: getCoverageStats() & getMethylationStats() have a new option
        called 'labels', if set to FALSE, no labels will be drawn on
        top of the histogram bars.

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: cov.bases option in tileMethylCounts() now works

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: methylRaw,methylBase and methylDiff objects have a new slot
        'resolution', which designates whether methylation information
        is base-pair resolution or regional resolution. allowed values
        'base' or 'region'

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: getCoverageStats() & getMethylationStats() now print the sample
        ids and methylation context automatically in the title when
        plotting

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: getCoverageStats() & getMethylationStats() takes extra options to
        be passed to hist() function

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: destrand option of unite() function will be over-ridden when
        methylRawList to be united contains regions rather than bases

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: unite() function checks if supplied elements of methylRawList
        object have the same context,assembly and resolution before
        uniting data sets.

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: validity function checking the format of the data for a methylRaw
        object is implemented

Version: 0.4
Category: IMPROVEMENTS AND BUG FIXES
Text: clusterSamples() and PCASamples() take methylation context
        information automatically and use it in plot titles

Version: 0.3.1
Category: IMPROVEMENTS AND BUG FIXES
Text: syntax error when fisher.test applied at calculateMethDiff() is
        fixed.

[MGFR](https://bioconductor.org/packages/MGFR)
----

Changes in version 0.99.2:

- First submission to BioConductor

[minfi](https://bioconductor.org/packages/minfi)
-----

Changes in version 1.19:

- preprocessNoob gets a dyeMethod argument which now allows for true
  single sample processing.

- combineArrayTypes is added; the intention is to be able to combine
  450k and EPIC array data at the RGChannelSet level.

- Support for early access IDAT files form the EPIC array.

- message() is now used instead of cat().

- Some functions moved from deprecated to defunct.

- Addressing a bug in preprocessQuantile which led to reduced
  performance for Type I probes when run with default paramters
  (stratified=TRUE).  Users are strongly encouraged to update to the
  latest version (1.19.7 or greater) and rerun the function.

- Extended combineArrayTypes to deal with control probes with the same
  address, but different characteristics (Color, Type, ExtendedType).
  Discussions with Illumina support reveals that, for control probes,
  same address is same probe.

- Extended combineArrayTypes to support [Genomic](Methyl|Ratio)Set.

- Fixed a bug that made detectionP fail with an error if used on only 1
  sample.

- Fixed a bug in read.metharray where we assumed a certain ordering is
  consistent in IDAT files from different samples. This is no longer
  assumed, but as a consequence the function is a bit slower.  Bug
  (indirectly) observed by Giovanni Calice <giovcalice@gmail.com>.

- Changing internals of MethylSet() to follow recent changes in
  assayDataElement<- in Biobase 2.33.1.

- Changing internals of MethylSet() (again) to follow recent changes in
  assayDataElement<- in Biobase 2.33.2.

- Fixing issues with combine() on various classes where pData columns
  doesn't have the same class. This translates to fixes for
  combineArrays and estimateCellType.

- estimateCellCounts gets a referencePlatform array (defaulting to
  450k) and now silently converts the input data to the desired
  platform using convertArray.

- Major refactoring of the annotation packages to reduce memory
  consumption.

[miRcomp](https://bioconductor.org/packages/miRcomp)
-------

Changes in version 1.3.3:

- Added miRcomp Shiny app.

- Added Lauren Kemperman as an author.

- Added CITATION file.

- Added NEWS.Rd file.

[mirIntegrator](https://bioconductor.org/packages/mirIntegrator)
-------------

Version: 1.3.1
Category: The warning message fixed: "Removed 3 rows containing missing
        values (geom_path)."
Text:

Version: 1.3.1
Category: microRNAadded <- microRNAadded[complete.cases(microRNAadded
Text:

Version: 1.3.1
Category: Plots were improved on this version
Text:

[miRNAmeConverter](https://bioconductor.org/packages/miRNAmeConverter)
----------------

Changes in version 1.1.6:

NEW FEATURES

- Added new attribute 'sequence', which is attached to the result data
  frame of the 'translateMiRNAName'-function.

SIGNIFICANT USER-LEVEL CHANGES

- Removed parameter 'correct' from 'checkMiRNAName' function.

- Updated citation

- Updated vignette

BUG FIXES

- Fixed bug regarding filter that neglected 'let' miRNAs, such as
  'hsa-let-7a'.

[miRNAtap](https://bioconductor.org/packages/miRNAtap)
--------

Changes in version 1.7.2:

- added 5th data source - miRDB, updated TargetScan and DIANA sources

[missMethyl](https://bioconductor.org/packages/missMethyl)
----------

Changes in version 1.7.3:

- modified gene set testing functions gometh and gsameth to accommodate
  EPIC arrays

[MLP](https://bioconductor.org/packages/MLP)
---

Changes in version 1.21.1:

- mlpBarplot: add ylab, cex, use title instead of mtext to set title

- plotGOgraph: add nCutDescPath, replace '\\\n' by \n', add error when
  graphs without edges

[MODA](https://bioconductor.org/packages/MODA)
----

Changes in version 0.99.4:

NEW FEATURES

- Initial version

USER-LEVEL CHANGES

- More comprehensive help pages

- BiocStyle Vignette

BUG FIXES

- use `seq_len(numP)` instead of `1:numP`, to avoid surprises when numP
  == 0

[monocle](https://bioconductor.org/packages/monocle)
-------

Version: 2.0.0
Text:

Version: 2.1.1
Category: BUGFIXES
Text: Fixed an problem where reduceDimension would return different
        results on repeated runs given the same inputs.The problem was
        actually in DDRTree in two places: kmeans and irlba. We now
        call irlba with deterministically initialized eigenvectors and
        kmeans with deterministically selected rows of the input.

Version: 2.1.1
Category: BUGFIXES
Text: Fixed a problem in classifyCells related to joining factors and
        levels. This would generate annoying warnings.

Version: 2.1.1
Category: BUGFIXES
Text: Fixed the check for valid sizeFactors prior to
        differentialGeneTest. Without this check, differentialGeneTest
        would report FAIL on all genes because of factors not having
        enough levels.

Version: 2.1.0
Category: FEATURES
Text: Monocle's algorithm for converting relative expression values
        (e.g. TPMs) into absolute transcript counts, called "Census"
        has been re-designed. The new version is much more accurate.
        Census is available through the relative2abs() function. The
        interface to this function has changed, and the output values
        will be quite different from the old version. The main change
        is that the version in 1.99.0 reported cDNA counts, as
        estimated in each cell's cDNA library. Now, relative2abs
        reports estimates of mRNA counts in the lysate. These are
        easier to compare with experiments that estimate mRNA counts
        via the use of spike-in controls.

Version: 2.1.0
Category: FEATURES
Text: A new heatmap function plot_pseudotime_heatmap() replaces the old
        plot_genes_heatmap().

Version: 2.1.0
Category: FEATURES
Text: Lots of new documentation. More to follow in upcoming releases.

Version: 2.1.0
Category: BUGFIXES
Text: Ordering cells with DDRTree no longer compresses cells at the
        tips of trajectories.

Version: 2.1.0
Category: BUGFIXES
Text: Pseudo counts are now applied correctly no matter what the
        underlying distribution used to model expression is.

Version: 2.1.0
Category: BUGFIXES
Text: Variance stabilization is applied correctly when a CellDataSet
        object's expressionFamily is negbinomial.size().

Version: 2.1.0
Category: BUGFIXES
Text: The gene_id field in results from calculateMarkerSpecificity is
        now a character instead of a factor. This was leading to
        indexing errors and nonsensical semi-supervised clustering and
        ordering results.

Version: 2.1.0
Category: BUGFIXES
Text: BEAM() and plot_branched_heatmap() now use cBind instead of
        cbind, fixing an issue with sparse matrices.

Version: 1.99.0
Text: The first public release of the Monocle 2 series. For a summary
        of new features and changes, please see:
        http://cole-trapnell-lab.github.io/monocle-release/features/

[msa](https://bioconductor.org/packages/msa)
---

Changes in version 1.5.5:

- fixes in ClustalOmega source code to ensure Windows compatibility of
  GCC6 compatibility fix

Changes in version 1.5.4:

- bug fix in msaClustalW(): unsupported parameter 'tree' deactivated

- fixes in ClustalOmega source code to ensure GCC6 compatibility

- fix in msaConvert() function to improve safety of call to suggested
  package 'phangorn'

Changes in version 1.5.3:

- additional conversions implemented for msaConvert() function

- corresponding changes in documentation

Changes in version 1.5.1:

- version number bumps for technical reasons related to Bioconductor
  build servers

Changes in version 1.5.0:

- new branch for Bioconductor 3.4 devel

[MSnbase](https://bioconductor.org/packages/MSnbase)
-------

Changes in version 1.99.6:

- Reverting to old initialize,Spectrum (see issue #163) <2016-10-07
  Fri>

- Setting Spectrum class versions outside of prototype (see issue
  #163). For this, there is now a vector of class version in
  .MSnbaseEnv <2016-10-10 Mon>

Changes in version 1.99.5:

- Add removeReporters,OnDiskMSnExp (see issue #161 for details)
  <2016-10-07 Fri>

Changes in version 1.99.4:

- Fix bug in readMzTabData_v0.9 <2016-10-07 Fri>

Changes in version 1.99.3:

- Injection time is now added to the header when reading mzML files
  using readMSData2 (see issue #159) <2016-10-04 Tue>

Changes in version 1.99.2:

- Added isolationWindow,MSnExp method <2016-09-23 Fri>

- Updated readMSData Spectrum2 class to support MS levels > 2
  <2016-09-30 Fri>

Changes in version 1.99.1:

- Fix MS level test in quantify,OnDiskMSnExp to support MS2+ level
  quantitation <2016-09-14 Wed>

Changes in version 1.99.0:

- Add normalize method for OnDiskMSnExp <2016-06-07 Tue>.

- Fix bug in readMSData <2016-06-07 Tue>.

- Added documentation and unit test for trimMz <2016-06-01 Wed>.

- Added trimMz method for OnDiskMSnExp objects <2016-05-31 Tue>.

- Added (internal) method spectrapply to apply a function to all
  spectra of an OnDiskMSnExp object; does the data import, subsetting,
  application of lazy processing steps etc. <2016-05-27 Fri>.

- Added a section to the MSnbase-development vignette <2016-05-27 Fri>.

- Finished with all pSet inherited methods and docs <2016-05-26 Thu>.

- Added rtlim argument to spectra method for OnDiskMSnExp <2016-05-26
  Thu>.

- Methods rtime, tic, ionCount, polarity and acquisitionNum implemented
  <2016-05-25 Wed>.

- Documentation added for most methods <2016-05-25 Wed>.

- Methods peaksCount and spectra for OnDiskMSnExp objects implemented
  <2016-05-24 Tue>.

- Performance and validation tests for the Spectrum1 C-constructor
  <2016-05-23 Mon>.

- Spectrum1 constructor implemented in C (presently not exported)
  <2016-05-20>

- length, scanIndex, acquisitionNum and centroided implemented
  <2016-05-20>

- Implemented fromFile and msLevel for OnDiskMSnExp <2016-05-19 Thu>

- Implemented an OnDiskMSnExp class for on-the-fly data import
  <2016-05-19 Thu>

- Fixed the validate method for pSet to play nicely with OnDiskMSnExp
  objects <2016-05-19 Thu>

- Added slot onDisk to pSet objects (TRUE for OnDiskMSnExp objects,
  FALSE otherwise).  The getter method isOnDisk checks for the presence
  of the slot in the object (backward compatibility) <2016-05-19 Thu>

- Implemented a ProcessingStep class that helps to keep track how
  (spectra) data should be processed on the fly <2016-05-18 Wed>

- In OnDiskMSnExp validity, check that the assaydata is empty
  <2016-06-28 Tue>

- Pass neutralLoss in plot,Spectrum,Spectrum-method to
  .calculateFragments; fixes #146 <2016-08-12 Fri>

- Allow the user to specify the `cex`, `lwd`, `pch` for peaks and
  fragments in plot,Spectrum,Spectrum-method; closes #148 <2016-08-12
  Fri>

- Update centroided with an na.fail argument (see issue #150 for
  details) <2016-08-12 Fri>

- Fix warning in readMgfData if TITLE contains multiple "=" <2016-08-24
  Wed>

Changes in version 1.21.8:

- Update MSnSet validity method to guard agains empty string feature
  names <2016-06-20 Mon>

- Simplify show,MSnExp method to work for various MS level cases and
  OnDiskMSnExp - addressed issue #98 <2016-06-28 Tue>

- removeMultipleAssignment also removes features that were not assigned
  (i.e. that have fcol (nprots) NA) <2016-07-02 Sat>

- New smoothed slot/accessor/replacement methods <2016-07-08 Fri>

- Update reporter masses and add TMT10 ETD/HCD <2016-07-20 Wed>

- returning empty spectrum when fliterMz has empty range - see issue
  #134 <2016-07-22 Fri>

- (mz, intensity) values are reordered based on order(mz) - see issue
  #135 <2016-07-26 Tue>

- fix bug in bin,Spectrum - see issue #137 <2016-08-05 Fri>

Changes in version 1.21.7:

- Update iPQF reference <2016-06-01 Wed>

- Fix a bug in normalize method for MSnExp objects: assigning
  normalized spectra directly to assayData is not possible, as the
  environment is locked. See PR #91.

- readMSData: if no phenodata is provided it creates an empty one with
  rownames corresponding to the file names. See PR #91.

- Lock itraqdata's assaydata bindings <2016-06-08 Wed>

Changes in version 1.21.6:

- MSmap unit test <2016-05-23 Mon>

- Fix bug in as(MSmap, "data.frame") <2016-05-23 Mon>

Changes in version 1.21.5:

- Added Johannes as contributor <2016-05-12 Thu>

- Deprecate MzTab v0.9 <2016-05-19 Thu>

- Fix old googlecode URLs to old MzTab <2016-05-19 Thu>

Changes in version 1.21.4:

- More MzTab and Spectrum1 unit testing <2016-05-08 Sun>

- Speed up readMSData (PR #86 by jotsetung) <2016-05-12 Thu>

- Replace example file URL to use github instead of googlecode
  <2016-05-12 Thu>

Changes in version 1.21.3:

- No fileNames replacement method <2016-05-07 Sat>

- fileNames unit tests <2016-05-07 Sat>

- add fileNames to class that had fileName accessor (MSmap, MzTab)
  <2016-05-07 Sat>

Changes in version 1.21.2:

- Check for rownames/fnames in readMSnSet2 and unit test <2016-05-05
  Thu>

Changes in version 1.21.1:

- Fix wrong indexing in readMSdata, msLevel==1 (PR #85 by jotsetung)
  <2016-05-04 Wed>

- grep/getEcols have a 'n' param specifying which line to grep/get
  <2016-05-04 Wed>

Changes in version 1.21.0:

- Version bump for new Bioc devel

[MSnID](https://bioconductor.org/packages/MSnID)
-----

Changes in version 1.7.3:

- Logical error fix in PSM FDR calculation. Prior it ignored redundancy
  in peptide to protein mapping. Now protein accession are discarded
  and unique rows considered before compuring the FDR.

- Added unit test to cover new PSM FDR calculation.

Changes in version 1.7.2:

- added high performance parser based on mzR::openIDFile

Changes in version 1.7.1:

- added infer_parsimonious_accessions method

- version leap from 1.3.1 to 1.7.1 to catch up with Bioconductor

- fixed conflict between reshape2 and data.table

[msPurity](https://bioconductor.org/packages/msPurity)
--------

Changes in version 0.99.10:

- Troubleshooting windows build failure

Changes in version 0.99.9:

- Offset bug fixed (previous only using extracting lower offset from
  mzML file)

- Updated handling of RT corrected xcmsSet objects for frag4feature
  function

- Additional column added for tracking ms/ms spectra

Changes in version 0.99.8:

- User option to change the mzR backend library

Changes in version 0.99.4:

- Troubleshooting mac build failure

Changes in version 0.99.3:

- Grouping multiple peaklist into one wide dataframe ** Peaklists can
  now be averaged across each class using the function groupPeaks() for
  the class purityD ** A list of dataframes can also be grouped togther
  using the function groupPeakEx()

Changes in version 0.99.2:

- Updated class names purityPD to purityD

- Updated class names purityLC to purityX

- Updated vignette to reflect slightly different terminology

- Added normalised TIC option for purityD msPurity v0.99.0 (Release
  date: 2016-04-08) Initial release!

[MSstats](https://bioconductor.org/packages/MSstats)
-------

Version: 3.5.5
Date: 2016-09-30
Text: BUG FIXES - dataProcess with fractionation sample when filling
        incomplete rows. Especially, not balanced fractionation for
        heavy and light, (heavy in one fractionation, no heavy in other
        fractionation) - groupComparison function : fix the issue with
        different columns from different summary Methods.  -
        MaxQtoMSstats function : option removeMpeptides=FALSE are now
        available. (Thanks, Danielle) - In case of multiple injections
        per sample (from fractionation or multiple injections with
        different range of m/x), normalization is performed separately
        and multiple injections are merged in each sample. NEW FEATURES
        - Add ‘originalRUN’ column in xx$ProcessedData after
        dataProcess function.  - Profile plot from dataProcessPlot
        distinguish censored missing data or not with different symbol.
        SIGNIFICANT CHANGES FOR METHOD - applied the algorithm for
        deciding the threshold of censoring.  - Method for calculation
        of the LOB/LOD is changed. LOQ is not calculated anymore.
        Please check help files for details.  -
        summaryMethod='logOfSum' option in dataProcess is retired.  -
        modelBasedQCPlots work with output from groupComparison in
        order to check the normality assumption of linear mixed effect
        model for whole plot level inference.

Version: 3.5.1
Text: Fix bug : summaryMethod=‘logOfSum’, redesign for result table.

[multiClust](https://bioconductor.org/packages/multiClust)
----------

Changes in version 6-19-16:

- -Package version pushed to 1.0.3 -Updated package title and abstract
  to reflect that of publication in Cancer Informatics DOI:
  10.4137/CIN.S38000.

Changes in version 3-2-16:

- -Package version pushed to 0.99.6 -Added Biobase, GEOquery, and
  preprocessCore to package suggests -Minor revisions in the package
  vignette

Changes in version 2-13-16:

- -Package version changed to 0.99.5 -Added option to specify FDR
  cutoff when using the Adaptive GMM method in the number_probes
  function. -Updated documentation in the vignette

Changes in version 2-11-16:

- -Packaged changed to version 0.99.4 -Revised the code in the vignette
  -Added explanation of using RNA-seq data with package in vignette
  -Revised code documentation for probe_ranking function

[MultiDataSet](https://bioconductor.org/packages/MultiDataSet)
------------

Version: 1.2.0
Category: o Add wrappers for two functions to integrate omic data: mcia
        (from omicade4 package) and iClusterPlus
Text:

Version: 1.2.0
Category: o Add advanced subsetting by phenotype and feature
Text:

Version: 1.2.0
Category: BUG FIXES
Text:

Version: 1.2.0
Category: o Solve subsetting issues when sampleNames is different from
        ID column.
Text:

[mzR](https://bioconductor.org/packages/mzR)
---

Changes in version 2.7.13:

- Fix waring on OSX (see issue #60)

Changes in version 2.7.12:

- Fixing last warning on MacOS

Changes in version 2.7.11:

- Fix failing example and test due to updated msdata::protemics data

Changes in version 2.7.10:

- compiling and loading on Windows (hopefully)

Changes in version 2.7.9:

- Updating Makevars.win

Changes in version 2.7.8:

- compile pwiz when installing on Windows <2016-09-26 Mon>

Changes in version 2.7.7:

- add netcdf static lib; potential incompatibility <2016-09-22 Thu>

Changes in version 2.7.6:

- new isolationWindow accessor <2016-09-23 Fri>

Changes in version 2.7.5:

- Apply Martin's free/delete patch - see
  https://github.com/sneumann/mzR/issues/52 <2016-09-22 Thu>

Changes in version 2.7.4:

- upgrade pwiz

- ensure file connections are cleaned and closed for pwiz.

Changes in version 2.7.1:

- revert added compiler switch, and fixed boost code instead to get rid
  of warnings

[netprioR](https://bioconductor.org/packages/netprioR)
--------

Changes in version 0.99.1:

- added examples and vignette

- Re-structured to S4

[normr](https://bioconductor.org/packages/normr)
-----

Changes in version 0.99.6:

FEATURES

- summary() for NormRFit instances shows now percentage of bins at each
  significance level

- getEnrichment() takes argument for computing a non-standardized
  enrichment

PERFORMANCE

- SIMD optimization for OpenMP >= ver4

Changes in version 0.99.5:

FEATURES

- summary() for NormRFit instances of types "diffR" and "regimeR" now
  shows componentwise test results

BUGFIXES

- Removing MAKEFLAG "-D_GLIBCXX_PARALLEL" because it is experimental

Changes in version 0.99.4:

BUGFIXES

- Removing PKG_CPPFLAGS and PKG_LIBS variables from src/Makevars

Changes in version 0.99.3:

BUGFIXES

- Fixing imports in roxygen2

Changes in version 0.99.2:

FEATURES

- Vignette finished

- Documentation for NormRFit-class & NormRCountConfig-class updated

BUGFIXES

- deploy.log removed

Changes in version 0.99.1:

FEATURES

- Added examples for normr-methods, BamCountConfig-class &
  NormRFit-class

- Moved all documentation to one manpage (normr.Rd)

- Vignette: knitr-based vignette for normr finished

- exportR(): Export NormRFit objects to bed, bedGraph and bigWig

- plot(): Produce a small diagnostic plot for NormRFit objects

- diffR(): Intersection of results for switched labels increases
  specificity

BUGFIXES

- Removed deploy.log

- Fixed build errors and warnings

- Fixed failing test

Changes in version 0.99.0:

FEATURES

- huge performance speedup by data compression to unique tupels only

- Kahan summation for log posterior summation to reduce numerical error

- Standardized enrichment computation in C++. This enrichment can be
  well compared between various ChIP-seq experiments

- exportR() function provides methods to export results of fits as a
  region bed file or bigWig track of standardized enrichment

[OncoScore](https://bioconductor.org/packages/OncoScore)
---------

Changes in version 1.1.2 (2016-09-28):

- Various bug fixes to the queries to PubMed.

[OncoSimulR](https://bioconductor.org/packages/OncoSimulR)
----------

Changes in version 2.4.0:

- Mutator phenotype and gene-specific mutation rates.

- End simulations stochastically as a function of size.

- Specify fitness by giving genotype-> fitness mapping.

- Random fitness landscape generation.

- Plots of fitness landscapes.

- Vignette: using Rmd.

- Several improvements in help and vignette.

- Improved test coverage.

- samplePop: sample at arbitrary sizes.

- evalAllGenotypes: order = FALSE by default.

Changes in version 2.3.17 (2017-09-22):

- random2 for rfitness.

- Vignette: decrease size and running time.

Changes in version 2.3.16 (2017-09-19):

- Help was not accurate for some probSize args

Changes in version 2.3.15 (2017-09-14):

- Typo in name (progresion)

Changes in version 2.3.14 (2017-08-26):

- Documentation improvements.

Changes in version 2.3.13 (2017-08-19):

- bioRxiv citation.

- Using Rmd for vignette.

- Improvements in vignette.

Changes in version 2.3.12 (2017-08-10):

- Increase N in some tests.

Changes in version 2.3.11 (2017-08-10):

- evalAllGenotypes: order = FALSE by default.

- Clarified difference plotFitnessEffects and plotFitnessLandscape.

Changes in version 2.3.10 (2017-07-14):

- The windows check failure with mc.cores.

Changes in version 2.3.9 (2017-07-09):

- Accessible genotypes in rfitness and plotFitnessLandscape

Changes in version 2.3.8 (2017-07-08):

- PDBasline default is now 1.2.

Changes in version 2.3.7 (2016-07-05):

- Unused C++ code reorganiz.

- Added tests

- Vignette and documentation improvements.

- End simulations stochastically as a function of size.

Changes in version 2.3.6 (2016-06-25):

- Improved test coverage

Changes in version 2.3.5 (2016-06-25):

- to_Magellan.

Changes in version 2.3.4 (2016-06-24):

- Failing some tests in Win 32-bits

Changes in version 2.3.3 (2016-06-23):

- Vignette improvements and typo fixes.

- rfitness: generate fitness landscapes.

- Plot of fitness landscapes.

- Specify fitness by giving genotype-> fitness mapping.

- Tests showing same gene in epist./DAG/order.

- Clarified internal C++ unique/sorted in genotypes.

- Checks initMutant correct and bug initMutant mutable pos.

- samplePop: sample at arbitrary sizes.

- plot.oncosimulpop using auto for color.

- Mutator phenotype and gene-specific mutation rates.

- Bug fixed: to_update set at 2 when mutating to pre-existing.

- Lots and lots of new tests.

[oppar](https://bioconductor.org/packages/oppar)
-----

Changes in version 1.2.0:

BUG FIXES

- Enhance opa() function: - Improve speed and simplicity of the code -
  Improve grouping normal/cancer samples. Previously a factor vector
  passed to group argument would have failed to identify normal/cancer
  samples correctly.

[OrganismDbi](https://bioconductor.org/packages/OrganismDbi)
-----------

Changes in version 1.16.0:

NEW FEATURES

- add check for missing OrgDb package in .taxIdToOrgDb()

- add 'orgdb' argument to makeOrganismDbFromBiomart()

MODIFICATIONS

- modify error message in .taxIdToOrgDb()

[PAA](https://bioconductor.org/packages/PAA)
---

Changes in version 1.7.1 (2016-05-13):

GENERAL

- Built with the latest R version (R-3.2.4 Revised).

- Built with the latest Rtools version (Rtools 3.3).

- Update of the PAA citation: PAA applications note, Turewicz et al.,
  Bioinformatics, 2016, PMID: 26803161, added to the CITATION file.

- Update of the URL in the DESCTIPTION file.

- Correction of some typos in the documentation.

- Update of the vignette to demonstrate the new features.

NEW FEATURES

- New function: batchFilter.anova() for multi-batch scenarios.

- New function: plotFeaturesHeatmap.2() as an alternative to
  plotFeaturesHeatmap().

- New arguments for loadGPR(): "description", "description.features"
  and "description.discard" for data import when gpr files don't
  provide the column "Description". Furthermore, new dummy data has
  been added to the package in order to demonstrate this new feature.

- Now the customized object classes EListRaw and EList contain the new
  component "array.type" added by loadGPR when the data is imported.

IMPROVEMENTS

- Additional and more informative error messages added.

- The function plotArray() is more flexible. Now, not only ProtoArrays
  can be plotted. Moreover, now the plot can be saved also as png file.

MODIFICATIONS

- The loadGPR() and plotArray()-argument "protoarray.aggregation" has
  been renamed to "aggregation" and the default value has been changed
  to "none".

- The loadGPR() argument "array.type" is mandatory now.

BUG FIXES

- Plots created during rlm normalization by normalizeArrays() were in
  logE scale and not in log2 scale. This has been fixed.

- Correction of regular expressions for rlm normalization in
  normalizeArrays() and for data import in loadGPR().

[PanVizGenerator](https://bioconductor.org/packages/PanVizGenerator)
---------------

Changes in version 1.1.3:

- Updating to PanViz 0.3.1

[PathoStat](https://bioconductor.org/packages/PathoStat)
---------

Changes in version 1.0.0:

- Relative Abundance plots (Stacked Bar Plot, Heatmap)

- Diversity plots (Alpha and Beta diversity, Exploratory Tree, BiPlot,
  Co-Occurrence)

- Differential Expression (Expression Plots, Limma)

- Confidence Region Plots

- PCA plots

- PCoA plots

- Alluvial Plots for longitudinal data

- Core OTU analysis

[pathVar](https://bioconductor.org/packages/pathVar)
-------

Version: 1.3.1
Date: 2016-10-13
Category: Updated pathVar to remove any warnings for Bioconductor
        release
Text:

[Pbase](https://bioconductor.org/packages/Pbase)
-----

Changes in version 0.13.3:

- Update mapping vignette index entry, as suggested by Mike Love
  <2016-09-13 Tue>

Changes in version 0.13.0:

- Bioc devel

[pbcmc](https://bioconductor.org/packages/pbcmc)
-----

Changes in version 1.1.1:

CODE

- `PAM50Filtrate` now supports a single sample i.e.  `drop=FALSE`
  (Thanks to anonymous reviewer).

- `PAM50Classify` now Force filtrate if not already performed. In
  addition, for `std="median"` checks for appropriate annotation i. e.
  EntrezGene.ID-Symbol tuples for the same EntrezGene.ID have the same
  Symbol (Thanks to anonymous reviewer).

- `PAM50Permutate` now checks for previous call to `classify` (Thanks
  to anonymous reviewer).

DOCUMENTATION

- The vignette has been update in order to include two additional
  examples (single subject and custom data) and computational
  requirements (e.g. memory, time) (Thanks to anonymous reviewer).

Changes in version 1.1.0:

VERSION

- Bump version after creating 3.4 devel branch

[pcaExplorer](https://bioconductor.org/packages/pcaExplorer)
-----------

Changes in version 1.99.0:

OTHER NOTES

- Reflecting the major feature added, will trigger a major version
  number bump. Welcome soon, pcaExplorer 2.0.0!

Changes in version 1.1.5:

NEW FEATURES

- Automated report generation - template available + editor in the app
  tab for advance user customization

- Support for state saving, in the global environment as well as with
  binary data

- All plots generated can be now exported with the dedicated button

- Added confidence ellipse for PCA plot

- Added 3d pca plot

- Added functions to automatically retrieve the annotation in format
  ready to use for the app

- Added profile explorer function, for plotting behaviour across all
  samples for subset of genes

- Added distribution plots

- Added pairwise correlation plot

- Added table to enhance readability of the gene finder plot, also by
  annotating sample names

BUG FIXES

- Minor typos fixed in the tabs

- Added option row.names to read.delim for allowing row names when
  uploading the data

OTHER NOTES

- Added extra info in the about section

- Instructions and vignette rewritten to reflect new design of the app

Changes in version 1.1.3:

BUG FIXES

- Remove y axis limits to gene boxplots

- Fixed: correct labels and colors assignements for genespca

[PGA](https://bioconductor.org/packages/PGA)
---

Changes in version 1.3.15:

- Fix some bugs

Changes in version 1.3.10:

- Fix the bug of citation file

- Update site of ensembl when using useMart

Changes in version 1.3.6:

- Add citation information

Changes in version 1.3.5:

- Fixed some bugs

- Add the funcion of constructing the Trinity protein DB

- and the "addGeneName4Ensembl" function

[philr](https://bioconductor.org/packages/philr)
-----

Changes in version 0.99.0:

USER-VISIBLE CHANGES

- Removed deprecated functions (c.to.nn and t.to.nn)

- Removed options for parallel processing in phylo2sbp, with the
  algorithmic speedups to this function parallel processing is
  superfluous and not used (even for trees of >45,000 leaves).

- Updated introduction to philr-intro vignette.

- Added Install instructions (from source) to readme

- Added citation info to readme (paper not on bioRxiv)
  http://biorxiv.org/content/early/2016/08/31/072413 Silverman JS,
  Washburne A, Mukherjee S, David LA. 2016.  A phylogenetic transform
  enhances analysis of compositional microbiota data.  bioRxiv doi:
  10.1101/072413

- News file is now parsed by news() function.

Changes in version 0.3.0:

INTERNAL CHANGES

- Fixed confusing difference between gp.rowMeans and g.colMeans (now
  gp.rowMeans -> g.rowMeans)

- Various other bug fixes

Changes in version 0.1.4:

USER-VISIBLE CHANGES

- Introduced new vignette ('philr-intro') based on Global Patterns
  dataset from phyloseq

- name.to.nn which (as well as nn.to.name) are now exported!

- Internal plotting functions replaced with annotate_balance and new
  geom_balance which was implemented in the package ggtree.

- Resolved anorm vs. enorm bug (previously anorm was calculating the
  euclidean norm due to subsetting behavior in compositions package).
  With this, also removed dependency on compositions package and
  reimplemented closure in philr.

INTERNAL

- gp.rowMeans and g.colMeans now handle calculation of geometric means
  for rows and columns respectively. Note 'gp' because rowMeans for
  rows needs to be calculated with weights. (See reference in the
  documentation for that function).

Changes in version 0.1.3:

USER-VISIBLE CHANGES

- Weighted / Genralized ILR functions now exported (shiftp, clrp, ilrp,
  buildilrBasep)

- Renamed function blw.mean.descendants to mean_dist_to_tips

- t.to.nn and c.to.nn replaced by name.to.nn which (as well as
  nn.to.name) are now vectorized.

Changes in version 0.1.0:

USER-VISIBLE CHANGES

- Basic functions from paper in added

- name.balance can now show vote tallies

- Philr function now warns if zeroes present

- Added convert_to_long function

[Pi](https://bioconductor.org/packages/Pi)
--

Changes in version 0.99.5:

NEW FEATURES

- A genomic-led target prioritisation system

- Leveraging genetic data to prioritise targets at the gene, pathway
  and network level

[piano](https://bioconductor.org/packages/piano)
-----

Changes in version 1.14:

NEW FEATURES

- Added new GSA method fgsea (fast GSEA) from the fgsea package.

- Added max method as an option in consensusScores().

IMPROVEMENTS [credits to Alexey Sergushichev

- Updated the p-value calculation (while using gene permutation) from:
  sum(background >= value) / length(background) to: (sum(background >=
  value) + 1) / (length(background) + 1) where bakground is a vector of
  permuted gene-set statistics and value is the actual calculated
  gene-set statistic. This change avoids returning p-values = 0,
  instead the theoretically smallest possible p-value will be
  determined by the number of permutations used. E.g. for nPerm =
  10000, the smallest p-value will be 1 / 10001 = 9.999e-05.

- Improved speed of checkLoadArg().

- Improved speed of gene-set statistic calculation for GSEA in
  calcGeneSetStat().

- Improved speed of FDR calculation for GSEA in fdrGSEA().

BUG FIXES

- Updated NAMESPACE to load dependencies in a better way.

- Updated the parsing of sbml files in loadGSC().

- Updated networkPlot() to conform to changes in the igraph package.

[Pigengene](https://bioconductor.org/packages/Pigengene)
---------

Changes in version 0.99.25 (2016-10-02):

Changes in existing functions

- The compute.pigengene() function now uses welch.pvalue() instead of
  pvalues.manov().

Changes in version 0.99.8 (2016-05-18):

General

- Under review by Bioconductor.

- Created.

[podkat](https://bioconductor.org/packages/podkat)
------

Changes in version 1.5.4:

- fix in vignette code that became necessary because of a change in the
  hg38 transcript annotation

Changes in version 1.5.3:

- readRegionsFromBedFile() now allows for including metadata columns

- corresponding changes in documentation

- fix of call to method rowRanges() into suggested package

- fix of citation info file

- fix of inst/NEWS

Changes in version 1.5.2:

- version number bump for technical reasons related to Bioconductor
  build servers

Changes in version 1.5.1:

- version number bump for technical reasons related to Bioconductor
  build servers

Changes in version 1.5.0:

- new branch for Bioconductor 3.4 devel

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

[procoil](https://bioconductor.org/packages/procoil)
-------

Changes in version 2.1.3:

- fix of citation info file

- fix of inst/NEWS

Changes in version 2.1.2:

- version number bump for technical reasons related to Bioconductor
  build servers

Changes in version 2.1.1:

- version number bump for technical reasons related to Bioconductor
  build servers

Changes in version 2.1.0:

- new branch for Bioconductor 3.4 devel

[pRoloc](https://bioconductor.org/packages/pRoloc)
------

Changes in version 1.13.16:

- removing visualTest from suggest, as plotting tests currently
  disabled (see https://github.com/MangoTheCat/visualTest/issues/19 for
  details) <2016-10-07 Fri>

Changes in version 1.13.15:

- Fix bug in plot,QSep with NAs <2016-09-09 Fri>

- Adapt tl vignette to new hpar <2016-09-15 Thu>

Changes in version 1.13.14:

- Pass ... in levelPlot,QSep <2016-09-06 Tue>

- Various improvements to the tutorial vignette <2016-09-07 Wed>

- Fix bug in QSep to honour fcol <2016-09-09 Fri>

Changes in version 1.13.13:

- Added plot3D, equivalent to plot2D, but using rgl in 3 dimensions
  <2016-09-03 Sat>

- Parametrise logging in knntlOptimisation to reduce object size
  <2016-09-06 Tue>

- New mrkHClust function to plot a dendrogram of subcellular markers
  <2016-09-06 Tue>

Changes in version 1.13.12:

- Remove verbose argument from getMarkerClasses function <2016-08-23
  Tue>

- plot2D has new addLegend argument <2016-08-23 Tue>

- addLegend default position is now bottomleft (requested by Lisa)
  <2016-08-23 Tue>

- new hexbin plot2D method <2016-08-28 Sun>

Changes in version 1.13.11:

- Add a norm agument to QSep's plotting functions to visualise
  normalised and raw distances <2016-08-22 Mon>

- New zeroInBinMSnSet function to visualise rowSums <2016-08-23 Tue>

Changes in version 1.13.10:

- Send biomart queries in getGOFromFeatues in chunks - see issue #85
  for details <2016-08-20 Sat>

Changes in version 1.13.9:

- new QSep infrastructure to assess spatial resolution <2016-08-17 Wed>

- Use %age of total variance in plot2D's scree plot <2016-08-17 Wed>

- fixed bug in combineThetaRegRes <2016-08-19 Fri>

Changes in version 1.13.8:

- fDataToUnknown accepts from = NA <2016-08-10 Wed>

Changes in version 1.13.7:

- Added an lda method to plot2D <2016-08-08 Mon>

Changes in version 1.13.6:

- Fixed bug in addGoAnnotations <2016-07-29 Fri>

Changes in version 1.13.5:

- Fixed bug in show method for ThetaRegRes <2016-06-08 Wed>

- Profiled knntl code so it's now much faster. Added more unit tests
  for knntl <2016-06-10 Fri>

- Fixed bug in plotDist x-axis labelling <2016-06-14 Tue>

Changes in version 1.13.4:

- Updated mouse pRolocmarkers <2016-06-02 Thu>

Changes in version 1.13.3:

- Fix bug in mirrorX/Y highlightOnPlot - see problem 1 in issue #79
  <2016-05-31 Tue>

Changes in version 1.13.2:

- Version bump to trigger package rebuilding now that purl()'ing issue
  has been correctly identified. knitr does not create purl()'ed
  (Stangle equivalent) .R files if _R_CHECK_TIMINGS_ is set, which the
  build system was setting. Now it's not set, so these .R files are now
  created. See https://github.com/yihui/knitr/issues/1212 for more.
  d.tenenbaum

Changes in version 1.13.1:

- Added keepNA argument to goTermToId so that if a GO term becomes
  obsolete and you cannot replace it with the ID name, you have the
  option to either replace it with a NA (previous and current default
  option) or now with keepNA = FALSE the term name will be replaced
  with the id name <2016-05-25 Wed>

- Bump version of all packages that use knitr for vignettes. This is
  because of an issue (now fixed) in knitr which failed to create
  purl()'ed R files from vignette sources and include them in the
  package. This version bump will cause these packages to propagate
  with those R files included. d.tenenbaum

Changes in version 1.13.0:

- Version bump for new Bioc devel

[pRolocGUI](https://bioconductor.org/packages/pRolocGUI)
---------

Changes in version 1.7.5:

- Various bug fixes <2016-10-07 Fri>

Changes in version 1.7.4:

- Add DT version dependency (see issue #71) <2016-08-09 Tue>

Changes in version 1.7.3:

- Fixed bug with matrix of markers in fcol. Now if one of the columns
  in the feature data is a matrix, it is converted to a vector using
  mrkMatToVec. This helps clarity and reduces wasted table space of 1's
  and 0's. <2016-07-29 Fri>

- Fixed bug in pRolocVis_compare application. Same issue as above with
  passing matrix as a fcol which manifested as an extra issue with
  zooming. Now fixed. <2016-07-29 Fri>

[Prostar](https://bioconductor.org/packages/Prostar)
-------

Version: 1.5.13
Category: NEW FEATURES
Text: In this version, the methods for missing values impuation from
        the package imp4p have been disabled because they are still
        under development. They will be available later.

Version: 1.5.11
Category: NEW FEATURES
Text: When the user wants to convert a text file into a MSnset dataset,
        the options are available only when a file has been loaded.

Version: 1.5.9
Category: NEW FEATURES
Text: This version is compliant with the version 1.5.9 of DAPAR

Version: 1.5.7
Category: NEW FEATURES
Text: The R commands of all the operations made with prostar are
        available in the logSession tab. This allows the user to create
        a script to automatize the analysis of a dataset.

Version: 1.5.5
Category: NEW FEATURES
Text: Prostar integrates the imp4p package for the imputation of
        missing values, w.r.t MNAR / MCAR

Version: 1.5.5
Category: NEW FEATURES
Text: In the export tool, the user can select the columns of metadata
        he wants to keep

Version: 1.5.3
Category: NEW FEATURES
Text: Prostar handles the format xlsx for Excel files

Version: 1.5.1
Category: NEW FEATURES
Text: Add modules

Version: 1.5.1
Category: BUG FIXES
Text:

[ProteomicsAnnotationHubData](https://bioconductor.org/packages/ProteomicsAnnotationHubData)
---------------------------

Changes in version 1.3.3:

- Fix formatting in vignette <2016-06-26 Sun>

[proteoQC](https://bioconductor.org/packages/proteoQC)
--------

Changes in version 1.9.4:

- Update the parameters for new version of ggplot2

[ProtGenerics](https://bioconductor.org/packages/ProtGenerics)
------------

Changes in version 1.5.1:

- added isCentroided generic <2016-07-23 Sat>

Changes in version 1.5.0:

- Version bump Bioc devel

[PureCN](https://bioconductor.org/packages/PureCN)
------

Changes in version 1.2.0:

Focus on PureCN 1.2 was to dramatically improve user-friendliness

- GATK requirement was dropped: PureCN now comes with functionality to
  generate coverage data when GATK is not available.

- Experimental support for VCFs generated by other callers than MuTect
  1.1.7

- Experimental automatic curation of results.

- Output plots were polished.

- Both numerical (NCBI-style) and non-numerical (UCSC-style) chromosome
  names (1 vs. chr1) supported.

- Better integration into existing copy number pipelines.

- Automatic COSMIC annotation.

- Thorough checks of input data, resulting in fewer crashes with
  unhelpful error messages.

- Documentation improvements.

OTHER MAJOR ENHANCEMENTS

- More thorough initial grid search should minimize cases where
  purity/ploidy is wrong because the solution was not considered.

- More robust and accurate likelihood model that provides robust
  estimates of minor segment copy numbers and more accurate LOH.

- Support for SNVs outside the target interval file
  (remove.off.target.snvs=FALSE).

- Improved segmentation and log-ratio normalization.

- Support for non-human samples.

- Support for 100% pure samples (matched normal mode only).

- Faster SNV fitting.

- Experimental support for correcting non-reference mapping bias.

API CHANGES

- New functions: autoCurateResults, bootstrapResults,
  calculateBamCoverageByInterval, calculateGCContentByInterval,
  calculateLogRatio, calculatePowerDetectSomatic, callAlterations,
  callAlterationsFromSegmentation, callLOH, filterTargets, getDiploid,
  getSexFromVcf, setMappingBiasVcf

- Deprecated functions: segmentationPSCBS

- Renamed functions: createExonWeightFile to createTargetWeights

- Renamed function arguments: exon.weight.file to target.weight.file
  gatk.normal.file to normal.coverage.file gatk.tumor.file to
  tumor.coverage.file coverage.cutoff to min.coverage

- Changed defaults: findFocal: made defaults less stringent, now 3Mb
  (instead of 2Mb) and minimum copy number of 5 (callAlterations still
  uses 6 as default).  runAbsoluteCN(genome): no default anymore
  (instead of 'hg19')

OTHER NEW FEATURES

- Functions to calculate power to detect mono-clonal and sub-clonal
  somatic mutations (Carter et al., Nature Biotech, 2012) were added.

PLANNED FEATURES FOR PURECN 1.4 (BIOCONDUCTOR 3.5

- Support for indels.

- Better runtime performance by ignoring unlikely solutions early.

- Better support for known, small deletions and amplifications (e.g.
  EGFRvIII, MYC)

- Better support for pool of normals.

- Switch to S4 data structures (maybe).

- Whole dataset visualizations (maybe).

[pwOmics](https://bioconductor.org/packages/pwOmics)
-------

Changes in version 3.1.1:

- include static consensus profiles function

- include activating/inhibiting edges in dynamic consensus net

[qcmetrics](https://bioconductor.org/packages/qcmetrics)
---------

Changes in version 1.11.3:

- Ammend biocViews <2016-09-23 Fri>

[qpgraph](https://bioconductor.org/packages/qpgraph)
-------

Changes in version 2.80:

USER VISIBLE CHANGES

- Updated documentation about qpPathWeight().

Changes in version 2.60:

USER VISIBLE CHANGES

- Added a first version of the function qpPathWeight() implementing the
  path weight calculations described in Roverato and Castelo, J. R.
  Soc. Ser.-C Appl. Stat., accepted.

Changes in version 2.40:

BUG FIXES

- Bugfix on qpPrecisionRecall() when argument 'refGraph' is a graphBAM
  object.

Changes in version 2.20:

USER VISIBLE CHANGES

- Updated the vignette "Estimate eQTL networks using qpgraph". It
  includes more detailed simulations illustrating the steps involved in
  the estimation of eQTL networks with qpgraph.

BUG FIXES

- Bugfix on the display of eQTL networks with hive plots

[qsea](https://bioconductor.org/packages/qsea)
----

Changes in version 0.99.6:

Bugfixes

- addCNV: check for sample table now checks specified files for CNV

Changes in version 0.99.5:

New features

- user level functions now have sample id check

- more detailed tutorial, including case study

Bugfixes

- makeTable: fixed groupMeans without individual samples

Changes in version 0.99.4:

New functions

- multicore support for bam file processing

- memory and runtime optimizations

Bugfixes

- bug in consideration of zygosity of chromosomes

Changes in version 0.99.3:

Bugfixes

- updated imported NAMESPACE

Changes in version 0.99.2:

New functions

- option to set ploidity for chromosomes, i.e. to account for sex
  chromosomes

Bugfixes

- fixed bug resulting in NA values at quantile calculation for beta
  estimates

- minor fixes

Documentation

- changed vignette builder to knitr/BiocStyle/html

Changes in version 0.99.1:

Bugfixes

- fixed crashing plotCoverage with one sample

- qseaPCA-show

- typo in error-message of getOffset (rpw)

- resolved namespace issues

Changes in version 0.99.0:

- Submission to Bioconductor

[QUBIC](https://bioconductor.org/packages/QUBIC)
-----

Changes in version 1.0.3:

- OpenMP is enabled

[rBiopaxParser](https://bioconductor.org/packages/rBiopaxParser)
-------------

Changes in version 2.12:

- Integrated functionality writen by Nirupama Benis
  (nirupama.benis@wur.nl). Added function pathway2Graph which creates a
  graph object containing all the interactions within a pathway. This
  is in contrast to pathway2RegulatoryGraph which only honors
  regulatory controls like activations and inhibitions.

- BUGFIX: fixed download link for Reactome Homo Sapiens data

- BUGFIX: fixed DESCRIPTION file. URL of the GitHub Repo is now:
  https://github.com/frankkramer-lab/rBiopaxParser

[RCyjs](https://bioconductor.org/packages/RCyjs)
-----

Changes in version 1.9.8:

BUG FIXES

- Significant (3x) speedup.  A 5000-node, 6000-edge graph transmits to
  Cytoscape from R in about 20 seconds.

Changes in version 1.8.0:

NEW FEATURES

- setNodeOpacityRule, controlling node fill color, border and/or label;
  interpolate & lookup modes both supported

- getNodeSize

- saveImage now supports pdf as well as png and svg formats

- setDefaultEdgeFontSize

- getAdjacentEdgeNames

SIGNIFICANT USER-VISIBLE CHANGES

- changed method names: layout -> layoutNetwork, version ->
  pluginVersion, get/setPosition -> get/setNodePosition

- NAMESPACE now imports four more methods from the graph package,
  helpful for package developers using RCytoscape: edgemode, addNode,
  addEdge, requested by Robert Flight.

BUG FIXES

- Changed getNodePosition node.name.delimiter to eliminate regex token,
  from ':.:' to ':-:' saveLayout now has optional 3rd parameter,
  'timestamp.in.filename'

- Fixed bug in setNodeLabelDirect.  Multiple nodes, one label now
  works.

- setCenter now casts x,y to numeric before sending out to CyRPC

[ReactomePA](https://bioconductor.org/packages/ReactomePA)
----------

Changes in version 1.17.4:

- unit test <2016-08-15, Mon>

Changes in version 1.17.3:

- move getDb from GOSemSim (no longer need this function) to ReactomePA
  and remove dependency of GOSemSim <2016-07-06, Wed>

- 'by' parameter for GSEA analysis <2016-07-04, Mon>

Changes in version 1.17.2:

- use byte compiler <2016-05-18, Wed>

- https://github.com/Bioconductor-mirror/ReactomePA/commit/6ce32c8e03e1b252662a07901cce022fab038086

Changes in version 1.17.1:

- https://github.com/Bioconductor-mirror/ReactomePA/commit/5d150f5fe545cfa3983872bf5485af1b9ba3129d

[readat](https://bioconductor.org/packages/readat)
------

Version: 0.99.0
Text:

[recount](https://bioconductor.org/packages/recount)
-------

Changes in version 0.99.30:

NEW FEATURES

- Added the function snaptron_query() which queries Intropolis via
  Snaptron to find if an exon-exon junction is present in the data.

Changes in version 0.99.29:

BUF FIXES

- Fixed an bug in the vignette. Thanks to Michael Love for noticing it!

Changes in version 0.99.0:

NEW FEATURES

- Created the package skeleton for recount

- Added the function reproduce_ranges() for re-creating the gene or
  exon level information used in the recount project.

- Added the function abstract_search() for identifying SRA projects of
  interest by searching the abstracts.

- Added the function browse_study() for opening a browser tab for
  further exploring a project.

- Added the function download_study() for downloading the data from the
  recount project.

- Added the function scale_counts() for properly scaling the counts
  before performing a differential expression analysis with the
  RangedSummarizedExperiment objects hosted in the recount project.

- Added the function expressed_regions() for defining the expressed
  regions in a chromosome for a given SRA study.

- Added the function coverage_matrix() for computing the coverage
  matrix based on the regions of interest for a given SRA study.

- Added the function geo_info() for obtaining sample information from
  GEO.

- Added the function find_geo() for finding the GEO accession id given
  a SRA run accession (id). This function will be useful for SRA
  projects that did not have GEO entries at the time recount's data was
  created.

- Added the function geo_characteristics() for building a data.frame
  from geo_info()'s results for the characteristics.

- Added the function all_metadata() which downloads all the phenotype
  data for all projects. This function can be useful for identifying
  projects and/or samples of interests.

[regionReport](https://bioconductor.org/packages/regionReport)
------------

Changes in version 1.7.10:

SIGNIFICANT USER-VISIBLE CHANGES

- Help pages now document advanced arguments.

Changes in version 1.7.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Dropped defunct functions.

[regsplice](https://bioconductor.org/packages/regsplice)
---------

Version: 1.0.0
Text:

[ReportingTools](https://bioconductor.org/packages/ReportingTools)
--------------

Changes in version 2015-3-27:

- Updated email for Jessica L. Larson

Changes in version 2013-4-1:

- Minor bug fixes in NAMESPACE and DESCRIPTION and css

- Enhanced vignettes to include information on our website and
  publication

- Enhanced vignettes to clarify use of .modifyDF()

- HTML reports are now represented by the HTMLReportRef referenceClass

- HTML output now fully customizable via .toHTML, .toDF and .modifyDF
  arguments to publish (see vignette)

- Publication mechanism is abstracted and customizable via
  ReportHandlers class

- ReportingTools output can be used within knitr documents and shiny
  Web applications (see vignettes knitr.Rmd and shiny.Rnw)

- Persistent representation of the HTML report being created is stored
  and accessible in the .reportDOM field of HTMLReportRef objects

- [[ and [[<- methods created for HTMLReportRef objects which allow
  selection, replacement and insertion of objects directly into reports

- Publish generic now accepts a 'name' argument.

- Existing reports can be read in via readReport, modified (via
  publish, [[<-, or direct manipulation of .reportDOM), and rewritten
  to file

- Path generic now returns a list/vector of the location slot values of
  the attached ReportHandlers object(s). These can be paths,
  connections, or other indications of report destination.

- Link generic function provided to build tables/sets of HTML links

- Added support for publishing ggbio and recordedplot objects

- CSS changed to Twitter Bootstrap

- Bugfixes to how NAs are handled when filtering and sorting columns

- New methods to handle output from running a glmLRT test in edgeR or
  nbinomTest in DESeq

- DEPRECATED: HTMLReport class is superseded by HTMLReportRef

- DEPRECATED: publication of HTMLReportRefs directly to a report (in
  order to make an index page) is no longer supported. Use the Link
  function.

- DEPRECATED: the page generic is not meaningful for HTMLReportRef
  objects (not all of which have a corresponding connection) and is
  deprecated. Use path instead.

[RGraph2js](https://bioconductor.org/packages/RGraph2js)
---------

Changes in version 1.0.1:

Improvement

- Remove the version() function from the Javascript code

[rGREAT](https://bioconductor.org/packages/rGREAT)
------

Changes in version 1.5.3:

- if background is provided, input regions are automaticially adjusted
  to make sure all input regions are subset of backgrounds.

[rhdf5](https://bioconductor.org/packages/rhdf5)
-----

Changes in version 2.18.0:

NEW FEATURES

- The low-level functions H5Pset_libver_bounds and H5Pget_libver_bounds
  is implemented. Creating files that can only be read by library
  versions 1.8 or later allows the usage of large attributes and
  improves performance.

USER VISIBLE CHANGES

- Per default all HDF5 files will be created with version 1.8 as lower
  bound. That means the created files can only be read with HDF5
  library versions >= 1.8. This changes allows the usage of large
  attributes and leads to performance improvements. If one wants to
  create a file that is readable with the earliest version of HDF5, one
  has to call H5Fcreate with fapl=h5default("H5P").

- Warning messages from the package C code can now be suppressed by the
  R-function suppressWarnings().

[RnBeads](https://bioconductor.org/packages/RnBeads)
-------

Changes in version 1.5.1:

- Reference-based cell type composition estimation now uses the 50,000
  most variable sites by default.

- Minor bug fixes

[roar](https://bioconductor.org/packages/roar)
----

Changes in version 1.9.1:

- Updated the vignette with some details about overlapping gene
  features

[rols](https://bioconductor.org/packages/rols)
----

Changes in version 2.1.3:

- Fix failing tests <2016-09-16 Fri>

[ropls](https://bioconductor.org/packages/ropls)
-----

Changes in version 1.5.22:

INTERNAL MODIFICATION

- minor code formatting

Changes in version 1.5.20:

BUG CORRECTION

- PCA: correction in the computation of variance and R2X (svd), and the
  determination of the optimal number of components

Changes in version 1.5.18:

INTERNAL MODIFICATION

- correction in the vignette (determination of optimal number of
  components in PCA)

Changes in version 1.5.16:

INTERNAL MODIFICATION

- minor vignette formatting

Changes in version 1.5.14:

NEW FEATURE

- 'fromW4M' function and 'toW4M' method to import preprocessed tables
  generated by the online Workflow4metabolomics infrastructure into an
  ExpressionSet bioconductor object (and vice versa)

Changes in version 1.5.12:

INTERNAL MODIFICATION

- minor vignette formatting

Changes in version 1.5.10:

INTERNAL MODIFICATION

- checking automated documentation with the 'roxygen2' package;
  automated vignette with rmarkdown and knitr

Changes in version 1.5.8:

INTERNAL MODIFICATION

- deleting missing files from svn repository

Changes in version 1.5.6:

NEW FEATURE

- opls can now be applied to ExpressionSet instances

- vignette is now in html format

- automated documentation with the 'roxygen2' package; automated
  vignette with rmarkdown and knitr

Changes in version 1.5.4:

NEW FEATURE

- package documentation now available with ?ropls, help("ropls"), or
  help("ropls-package")

Changes in version 1.5.2:

BUG FIXED

- stop messages generated when non-significant predictive or first
  orthogonal components are obtained

Changes in version 1.5.0:

NEW FEATURES

- 'opls' is now an S4 class: please use the accessors e.g.
  'getScoreMN(oplsModel)' instead of 'oplsModel$scoreMN', or
  'getScoreMN(oplsModel, orthoL = TRUE)' instead of
  'oplsModel$orthoScoreMN'

- The full list of accessors is: getLoadingMN, getScoreMN, getVipVn,
  getWeightMN (all with the orthoL argument), in addition to
  getSummaryDF, getPcaVarVn, getSubsetVi

- Please see the vignette and documentation for examples

[rpx](https://bioconductor.org/packages/rpx)
---

Changes in version 1.9.4:

- Updating unit tests <2016-10-04 Tue>

Changes in version 1.9.3:

- Update vignette to use readMzTabData v0.9 <2016-07-24 Sun>

[Rsamtools](https://bioconductor.org/packages/Rsamtools)
---------

Changes in version 1.25:

NEW FEATURES

- idxstatsBam() quickly summarizes the number of mapped and unmapped
  reads on each sequence in a BAM file.

BUG FIXES

- *File constructors now check that the file argument is length 1, and
  that the index argument is length 0 or 1.

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

Version: 1.24.0
Category: NEW FEATURES
Text:

Version: 1.24.0
Category: o New parameters in featureCounts() function: `fracOverlap'
        and `tmpDir
Text:

Version: 1.24.0
Category: o FeatureCounts can now perfrom fractional counting for both
        multi-mapping reads and multi-overlapping reads via its
        `fraction' parameter
Text:

Version: 1.24.0
Category: o Depreciate `PE_orientation' parameter in featureCounts
Text:

Version: 1.24.0
Category: o In-built annotation for hg38(GRCh38) is added to the
        package
Text:

Version: 1.24.0
Category: o Improved reporting of mapping quality score (MQS) for
        align() and subjunc() functions
Text:

Version: 1.24.0
Category: o Improved efficiency of exactSNP() function when calling
        SNPs from data with very high sequencing depth
Text:

Version: 1.24.0
Category: o Improved documentation and screen output
Text:

[RTCA](https://bioconductor.org/packages/RTCA)
----

Changes in version 2009-07-13:

- combineRTCA(list): Additional column is renamed into Plate. The vlues
  is evaluated from list item names. When the list has no name, an
  integer index beginning from 1 is used. Special attentions to list
  partially with names is noted in the documentation.

- parseRTCA(file, dec=".",phenoData, skipWell,...): Example is added in
  the documentation how to import pre-configured phenoData. Details
  section in the documentation is re-written to describe the process of
  parsing.

- RTCA-class: Experiment ID added to RTCA class

- Makefile: add Makefile to simplify common tasks like check and
  install

- plotGridEffect: takes 'column' instead of 'col' as mode parameter,
  and renders the mode as the title of the legend. Documentation
  updated.

- plotRTCA: is removed from the package and is substituted by the plot
  function.

[RUVSeq](https://bioconductor.org/packages/RUVSeq)
------

Changes in version 1.7:

- Fixed an error with NA and .isWholeNumber

[S4Vectors](https://bioconductor.org/packages/S4Vectors)
---------

Changes in version 0.12.0:

NEW FEATURES

- Add n-ary "merge" method for Vector objects.

- "extractROWS" methods for atomic vectors and DataFrame objects now
  support NAs in the subscript. As a consequence a DataFrame can now be
  subsetted by row with a subscript that contains NAs. However that
  will only succeed if all the columns in the DataFrame can also be
  subsetted with a subscript that contains NAs (e.g. it would fail at
  the moment if some columns are Rle's but we have plans to make this
  work in the future).

- Add "union", "intersect", "setdiff", and "setequal" methods for
  Vector objects.

- Add coercion from data.table to DataFrame.

- Add t() S3 methods for Hits and HitsList.

- Add "c" method for Pairs objects.

- Add rbind/cbind methods for List, returning a list matrix.

- aggregate() now supports named aggregator expressions when 'FUN' is
  missing.

SIGNIFICANT USER-VISIBLE CHANGES

- "c" method for Rle objects handles factor data more gracefully.

- "eval" method for FilterRules objects now excludes NA results, like
  subset(), instead of failing on NAs.

- Drop "as.env" method for List objects so that as.env() behaves more
  like as.data.frame() on these objects.

- Speed up "replaceROWS" method for Vector objects when 'x' has names.

- Optimize selfmatch for factors.

DOCUMENTATION IMPROVEMENTS

- Add S4QuickOverview vignette.

DEPRECATED AND DEFUNCT

- elementLengths() and compare() are now defunct (were deprecated in
  BioC 3.3).

- Remove "ifelse" methods for Rle objects (were defunct in BioC 3.3),

BUG FIXES

- Fix bug in showAsCell(x) when 'x' is an AsIs object.

- DataFrame() avoids NULL names when there are no columns.

- DataFrame with NULL colnames are now considered invalid.

[SC3](https://bioconductor.org/packages/SC3)
---

Changes in version 1.1.7:

- Added integration with scater package

- Added full non-interactive support 04/08/2016

Changes in version 1.1.6:

- Optimised distance and consensus calculations using Rccp

- Optimised some other bits and pieces

- SC3 is now up to 2 times faster for 1000 cells datasets 28/06/2016

Changes in version 1.1.5:

- Added svm.train.inds parameter - it allows a user to train SVM on a
  selected subset of cells 15/04/2016 The latest version 0.99.37
  contains the following updates

- tSNE panel is added to the interactive session

- additional table has been added to the output - it contains original
  and new cell labels in the order corresponding to the order of cells
  in the input expression matrix

- three functions corresponding to biological interpretation are now
  exported and can be used manually: get_de_genes, get_marker_genes,
  get_outl_cells. See their documentation for more details

- several bugs have been fixed 01/03/2016 The latest version 0.99.23
  contains the following major updates

- Major redesign

- Added description panels

- Fixed a bug with RSelenium

- Added a proper Excel export using WriteXLS library 10/12/2015 The
  first version 0.99.0 is submitted to Bioconductor

[scater](https://bioconductor.org/packages/scater)
------

Changes in version 1.1.26:

- Key code ported to C++ for greater computational and memory
  efficiency

- Added support/integration for SC3 package for single-cell clustering

- writeSCESet() function added to write SCESets to HDF5 format on disk

- mergeSCESet() function to merge two SCESet objects that share
  features

- plotPlatePosition() function to visualise gene expression and cell
  metadata for cells in their position on a plate

- Adding plotExprsVsTxLength() to plot expression against transcript
  length

- Added fitted line and some tweaks to plotExprsFreqVsMean().

- Adding support for scaled TPM counts when summarising expression at
  feature level.

- Adding NULL method for set_exprs(). Adding tests.

- Adding import of effective feature lengths with readKallistoResults()

- runSalmon() function for running Salmon from within R, following
  suggestions from Rob Patro.

- Added cellNames<- assignment function

- Added extra QC metrics

- Numerous other bug fixes and minor improvements

[scran](https://bioconductor.org/packages/scran)
-----

Changes in version 1.1.10:

- Transformed correlations to a metric distance in quickCluster().

- Removed normalize() in favour of scater's normalize().

- Switched isSpike()<- to accept a character vector rather than a
  logical vector, to enforce naming of spike-in sets. Also added
  warning code when the specified spike-in sets overlap.

- Allowed compute*Factors() functions to directly return the size
  factors.

- Added selectorPlot() function for interactive plotting.

- Switched to a group-based weighted correlation for one-way layouts in
  correlatePairs() and correlateNull(), and to a correlation of
  residuals for more complex design matrices.

- Added phase assignments to the cyclone() output.

- Implemented Brennecke et al.'s method in the technicalCV2() function.

- Updated convertTo() to store spike-in-specific size factors as
  offsets.

- Moved code and subsetting into C++ to improve memory efficiency.

- Switched to loess-based trend fitting as the default in trendVar(),
  replaced polynomial with semi-loess fitting.

- Added significance statistics to output of decomposeVar(), with only
  the p-values replaced by NAs for spike-ins.

- Updated documentation and tests.

[SeqArray](https://bioconductor.org/packages/SeqArray)
--------

Changes in version 1.14.0:

- the version number was bumped for the Bioconductor release version
  3.4

Changes in version 1.12.0-1.12.9:

- the version number was bumped for the Bioconductor release version
  3.3

- `seqVCF_SampID()`, `seqVCF_Header()` and `seqVCF2GDS()` allow a
  connection object instead of a file name

- "$num_allele" is allowed in `seqGetData()` and `seqApply()` (the
  numbers of distinct alleles)

- a new option '.progress' in `seqAlleleFreq()`, `seqMissing()` and
  `seqAlleleCount()`

- 'as.is' can be a `gdsn.class` object in `seqApply()`

- v1.12.7: a new argument 'parallel' in `seqApply()`, BiocParallel
  integration in `seqParallel()` and a new function `seqBlockApply()`

- v1.12.8: a new function `seqGetParallel()`

[SeqVarTools](https://bioconductor.org/packages/SeqVarTools)
-----------

Changes in version 1.11.3:

- Add hethom method to calculate heterozygosity / non-reference
  homozygosity in one step

- Add countSingletons method

Changes in version 1.11.1:

- Add variantData slot to SeqVarData class

[SGSeq](https://bioconductor.org/packages/SGSeq)
-----

Changes in version 1.8.0:

- Bug fixes and other improvements

[signeR](https://bioconductor.org/packages/signeR)
------

Changes in version 1.0.0:

- First Bioconductor release.

Changes in version 0.99.17:

- Reduced memory usage.

[sincell](https://bioconductor.org/packages/sincell)
-------

Changes in version 1.5.3:

- added extra information and examples about how to work with the
  clusters

[SNPhood](https://bioconductor.org/packages/SNPhood)
-------

Changes in version 1.3.4:

BUG FIXES

- updated procedure to parse BAM headers. Fixes error when performing
  consistency checks of chromosome sizes from BAM files. Fixes the
  error: Error in if (as.numeric(chrSizes[names(chrSizes)[i]]) !=
  chrSizes.df$size[chrSizes.df$chr == : argument is of length zero

[SNPRelate](https://bioconductor.org/packages/SNPRelate)
---------

Changes in version 1.8.0:

- add a new function `snpgdsIndivBeta()`

Changes in version 1.6.0-1.6.6:

- the version number was bumped for the Bioconductor release version
  3.3

- new implement of thread pool

- bitwise intrinsics (SSE2/AVX2) to accelerate `snpgdsIBSNum()`,
  `snpgdsIBS()`, `snpgdsIBDMoM()`, `snpgdsIBDKing()` (+50% to +300%)

- v1.6.4: bug fix in v1.6.3 (allele counting error with SSE2 implement)

- v1.6.5: `snpgdsGRM()`, renames the option "Visscher" to "GCTA", new
  option 'dosage' in `snpgdsPairScore()`, new function
  `plot.snpgdsPCAClass()`

[specL](https://bioconductor.org/packages/specL)
-----

Changes in version 1.7.4 (2016-05-19):

- USER VISIBLE CHANGES
  
  * added to specLSet summary "which std peptides (iRTs) where found in
  which raw files"
  
  * one plot per raw file in plot methode of specLSet object

Changes in version 1.7.1 (2016-05-13):

- USER VISIBLE CHANGES
  
  * replaced NEWS by NEWS.Rd file
  
  * modified specL object /replace decoy by score attribute #1,#4
  
  * specL on BioC 3.3

[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
--------------------

Version: 1.4.0
Category: NEW FEATURES
Text: Add makeSummarizedExperimentFromDataFrame() function.

Version: 1.4.0
Category: NEW FEATURES
Text: Add "acbind" and "arbind" methods for Matrix objects.

Version: 1.4.0
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text: Speed up "cbind" method for SummarizedExperiment objects based on
        a suggestion by Peter Hickey.

Version: 1.4.0
Category: DEPRECATED AND DEFUNCT
Text: Remove exptData() getter and setter (were defunct in BioC 3.3).

Version: 1.4.0
Category: BUG FIXES
Text:

[SWATH2stats](https://bioconductor.org/packages/SWATH2stats)
-----------

Changes in version 1.3.10:

BUG FIXES

- sample_annotation: add stop if column for files is not present,
  instead of printing a warning

- convert4aLFQ: improve warning message

- tests: add tests

Changes in version 1.3.9:

NEW FEATURES

- sample_annotation: improved architecture of function

- tests: add tests

- MSstats_data: add data in MSstats format to use for testing

- assess_decoy_rate: add stop instead of warnings if no decoy column is
  present

- mscore4pepfdr, mscore4assayfdr: add default value in function
  description

Changes in version 1.3.8:

BUG FIXES

- DESCRIPTION: move MSstats from suggests to enhances

Changes in version 1.3.7:

BUG FIXES

- test_filtering: tests for count_analytes

Changes in version 1.3.6:

NEW FEATURES

- filter_all_peptides: Improvement in display of protein Names using
  collapse

DOCUMENTATION

- Vignette: make MSstats and aLFQ eval=FALSE to prevent error when
  these are not build properly

Changes in version 1.3.5:

NEW FEATURES

- filter_mscore_fdr: add fail-safe in case target protein FDR cannot be
  reached.

- plot_variation, plot_variation_vs_total,
  plot_correlation_between_samples, count_analytes: add tests

Changes in version 1.3.4:

NEW FEATURES

- convert4mapDIA: improve error message in the case that multiple
  values exist per condition

- plot_correlation_between_samples: add option to not write any labels

BUG FIXES

- plot_variation: improved function to also deal with data that did not
  have the same size in the different conditions

DOCUMENTATION

- DESCRIPTION: update packages to import

- convert4MSstats: updated warning message

Changes in version 1.3.3:

BUG FIXES

- Add fail-safe functionality to sample_annotation function for
  non-unique file names.

Changes in version 1.3.0:

NEW FEATURES

- Development version of SWATH2stats in BioC 3.4

Changes in version 1.2.3:

BUG FIXES

- Add fail-safe functionality to sample_annotation function for
  non-unique file names.

[TarSeqQC](https://bioconductor.org/packages/TarSeqQC)
--------

Changes in version 1.3.13:

CODE

- Modification in readFrequencies method to fix countOverlaps error

- Modification in plot methods to capitalize axis names and titles

Changes in version 1.3.12:

CODE

- Modification in readFrequencies method in order to count only those
  reads falling out the features and non the overlapped.

- Modification in plotRegion method to correclty plot overlapped
  features

Changes in version 1.3.1:

CODE

- Addition of a new class called TargetExperimentList that allows the
  the comparisson of several TargetExperiment objects obtained using
  the same bed file

[TCGAbiolinks](https://bioconductor.org/packages/TCGAbiolinks)
------------

Changes in version 2.2.0:

NEW FEATURES

- Add GISTIC2 information and Mutation information in the
  summarizedExperiment while preparing GDC data

SIGNIFICANT USER-LEVEL CHANGES

- Speed up GDCquery

[TEQC](https://bioconductor.org/packages/TEQC)
----

Changes in version 3.13.2:

- since 'reduce' method for 'RangedData' objects was deprecated, copied
  its content to function within package. Transition from 'RangedData'
  to 'GRanges' objects all over the package should be pursued...

Changes in version 3.13.1:

- added flag 'isSecondaryAlignment=FALSE' for reading bam files in
  'get.reads', such that only unique (primary) alignments are
  considered for TEQC analysis

[TFBSTools](https://bioconductor.org/packages/TFBSTools)
---------

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

[TPP](https://bioconductor.org/packages/TPP)
---

Version: 2.3.1
Category: Update version number after the new Bioconductor Release 3.3
        was launched on
Text:

Version: 2.3.1
Category:
        https://www.bioconductor.org/packages/release/bioc/html/TPP.html
Text:

Version: 2.3.1
Category: The version number stays in accordance with Devel version in
        the new Release
Text:

Version: 2.3.2
Category: Removed unit test that causes R CMD check to crash since the
        latest update of package 'testthat
Text:

Version: 2.3.3
Text: Bugfixes in the functions responsible for filtering by
        user-specified QC columns (e.g. 'qssm') for normalization set
        creation. It is now possible to leave the 'otherRequirements'
        slot in the normalization criteria empty.

Version: 2.3.3
Text: Fixed typos in function documentation of 'analyzeTPPTR()'

Version: 2.3.3
Text: Fixed bug in plotColors for the case when no comparisons are
        specified

Version: 2.3.4
Text: Improved documentation of non-exported functions by ensuring that
        there is at least a one-line comment explaining what the
        function does.

Version: 2.99.0
Text: Major update: Extended package to analyse 2D-TPP experiments

Version: 2.99.0
Text: Major update: Extended package to analyse TPP-TR experiments by
        "non-parametric analysis of response curves" (NPARC), a
        non-parametrc spline-based test for treatment effects.

Version: 2.99.0
Text: Bugfix in 'tpptrCurveFit': restrict file names to 255 characters
        when creating melting curve plots in order to prevent file
        system crashes.

Version: 2.99.0
Text: New packages dependencies: dplyr, margrittr

Version: 2.99.0
Text: Bugfix in CCR workflow: avoid division by 0 when transforming
        foldchanges before dose-response curve fitting.

Version: 2.99.1
Text: Remove the complex reference data object from the data folder.
        Will be re- introduced lalter in a leaner format.

Version: 2.99.2
Text: Reduce data by "xz" compression

Version: 2.99.3-5
Text: Adapt notation of function names and input arguments for
        consistency with 1D (TPP-TR and TPP-CCR) part.

[trackViewer](https://bioconductor.org/packages/trackViewer)
-----------

Changes in version 1.9.14:

- add controls for lolliplot yaxis.

Changes in version 1.9.13:

- the height of features could be unit object.

Changes in version 1.9.12:

- Update the vignette.

- Fix the position of ylab.

Changes in version 1.9.11:

- Add caterpillar layout to lolliplot.

Changes in version 1.9.10:

- Fix the bug that if the features is outside of granges for lolliplot.

Changes in version 1.9.9:

- Fix the bug that no xscale is drawn when there is no positive value
  for track 1.

Changes in version 1.9.8:

- lolliplot could accept multiple features for one layer.

Changes in version 1.9.7:

- lolliplot could accept multiple types for multiple layers.

Changes in version 1.9.6:

- fix the bugs for the lolliplot when the last snp need to be jittered.

Changes in version 1.9.5:

- fix the bugs in the documentation.

Changes in version 1.9.4:

- add pie.stack layout for lolliplot.

Changes in version 1.9.3:

- add controls for lolliplot labels.

Changes in version 1.9.2:

- update documentation to include the how to set ylim.

[TRONCO](https://bioconductor.org/packages/TRONCO)
------

Changes in version 2.4.3:

- New algorithms: Prim, Chow-Liu, Prim and Edmonds.

- New scores: PMI, CPMI, MI

- A better bootstrap implementation

Changes in version 2.4.2:

- Implementation of maximum spanning tree algorithms.

- Implementation of a noise model for cumulative phenomena.

Changes in version 2.4.1:

- A better implementation of remove.cycles.

[TVTB](https://bioconductor.org/packages/TVTB)
----

Changes in version 0.99.8 (2016-10-12):

Minor changes

- Set vignette output format to BiocStyle::html_document2.

Changes in version 0.99.7 (2016-10-12):

Minor changes

- Fixed a closign bracket in NEWS file.

Changes in version 0.99.6 (2016-10-12):

Experimental changes

- Replaced BiocStyle::pdf_document2() by BiocStyle::pdf_document(); the
  former fails at the pandoc step of vignette production, the latter
  does not.

Changes in version 0.99.5 (2016-10-12):

Bug fix

- Fixed Collate: field of the DESCRIPTION file.

Changes in version 0.99.4 (2016-10-12):

Experimental changes

- Reverted changes applied in version 0.99.3. BiocParallel does not
  seem to be causing the build error on Windows Server.

Changes in version 0.99.3 (2016-10-12):

Experimental changes

- Disabled BiocParallel code to see if it resolves build errors on the
  _Bioconductor_ Windows Server. Note that the documentation was not
  yet updated to reflect this experimental change; this is intended to
  facilitate code reversion-or alternatively document the change-in the
  next commit.

Changes in version 0.99.2 (2016-10-11):

Bug fixes

- Use suffix accessor in add*Frequencies methods.

Changes in version 0.99.1 (2016-10-10):

Major changes

- New dedicated Genotypes class to store homozygote reference,
  heterozygote, and homozygote alternate genotype codes, along with the
  suffixes that define the INFO keys used to store their respective
  data in the VCF object.

- Removed families of methods tabulate* and density* from the
  NAMESPACE. The features may be revisited in the future.  The
  associated code and documentation was saved in the inst/sandbox
  subfolder for future reference.

- New slot svp in TVTBparam class to store ScanVcfParam objects. New
  associated accessor methods.  Moreover, TVTBparam may be coerced to
  ScanVcfParam.

- New signatures for method readVcf that supports param=TVTBparam, and
  optional phenotypes.  The method stores TVTBparam in the metadata
  slot of the VCF object, and phenotypes using the colData accessor.

- TVTBparam are no longer an argument of downstream methods; instead,
  they must be stored in metadata(vcf)[["TVTBparam"]]

Minor changes

- hRef and hAlt accessors renamed to refand alt, respectively.

- suffix accessor to returned named character vector for classes
  Genotypes and TVTBparam.

- Removed functions relevant only to the _Shiny_ application from the
  NAMESPACE (getEdb, EnsDbFilter, chr2file).

- VcfFilterRules can also store instances of the parent FilterRules
  class.

- Defined default return value for accessors vep and type, to avoid
  unnecessary switch statements.

- Simplified code of inherited methods following updates to the
  relevant packages (_e.g._ S4Vectors).

- Better respect of coding standards: removed superfluous usage of
  explicit argument naming, removed superfluous initialize methods.

- Man pages, vignettes, unit tests and _Shiny_ application updated to
  reflect changes to the package.

Changes in version 0.99.0 (2015-09-15):

New features

- First release submitted to the Bioconductor review process See
  DESCRIPTION file for details.

[Uniquorn](https://bioconductor.org/packages/Uniquorn)
--------

Changes in version 1.0.8 (2016-06-17):

Optimization confidence score

- Fixed the installation of the package

Changes in version 1.0.7 (2016-06-13):

Optimization confidence score

- Optimized the way confidence scores are calculated to increase the
  sensitivity and specificity of the method

Changes in version 1.0.6 (2016-06-05):

Optimization of default confidence score

- Adjusted default confidence score to optimal threshold

Changes in version 1.0.5 (2016-05-30):

Minor Bugfixes

- Fixed errors connected to adding CLs

Changes in version 1.0.4 (2016-05-27):

Introduction of confidence score

- Confidence score is the central threshold now

- Default is 3.0

- The score is the negative log e of the q-value

Changes in version 1.0.3 (2016-05-18):

Changes to treshold calculation

- Fixed error in threshold calculation

Changes in version 1.0.1 (2016-05-09):

Minor update of BED files

- Extended length of BED file mutations shown in the IGV

[uSORT](https://bioconductor.org/packages/uSORT)
-----

Changes in version 0.99.4 (2016-10-12):

MODIFICATION

- commented some example codes in monocle_wrapper, due to some unknow
  error.

Changes in version 0.99.2 (2016-09-13):

MODIFICATION

- updated according to bioCheck results for version 0.99.1.

Changes in version 0.99.1:

MODIFICATION

- updated according to first bioCheck results.

Changes in version 0.99.0 (2016-08-22):

MODIFICATION

- added the main function for uSORT named uSORT in file uSORT_main.R

- added function uSORT_sorting_wrapper

- added GUI function uSORT_GUI

- structured all codes for preprocessing to file uSORT_preProcess.R

- structured all codes for postprocessing to file uSORT_postProcess.R

- structured all codes for gene selection to file uSORT_GeneSelection.R

- added documentations for exported functions

- added one vignettes for the package

- added test code for the pacakge

- added namespace and dependencies declarition

[variancePartition](https://bioconductor.org/packages/variancePartition)
-----------------

Changes in version 1.3.11:

- in canCorPairs() and other functions, convert formula with
  as.formula()

- improve error messages for canCorPairs()

Changes in version 1.3.10:

- Add plotStratify()

- Update documentation

Changes in version 1.3.8:

- Add additional examples to vignette

- show projected memory usage of fitVarPartModel()

Changes in version 1.3.7:

- fitVarPartModel warns if names in exprObj and data are not identical

- residuals() and other functions deal with missing values properly

Changes in version 1.3.6:

- Small changes to vignette

Changes in version 1.3.5:

- Fix Bioconductor error

Changes in version 1.3.4:

- Fix typos

Changes in version 1.3.3:

- Improve documentation

[VariantAnnotation](https://bioconductor.org/packages/VariantAnnotation)
-----------------

Changes in version 1.20.0:

NEW FEATURES

- add import() wrapper for VCF files

MODIFICATIONS

- use now-public R_GetConnection

- remove defunct readVcfLongForm() generic

- remove 'genome' argument from readVcf()

- improvements to VCF to VRanges coercion

- support Varscan2 AD/RD convention when coercing VCF to VRanges

- use [["FT"]] to avoid picking up FTZ field

- summarizeVariants() recognize '.' as missing GT field

- document scanVcfheader() behavior for duplicate row names

BUG FIXES

- ensure only 1 matching hub resource selected in filterVcf vignette

- fix check for FILT == "PASS"

- correct column alignment in makeVRangesFromGRanges()

- fix check for AD conformance

[VariantFiltering](https://bioconductor.org/packages/VariantFiltering)
----------------

Changes in version 1.10:

USER VISIBLE CHANGES

- Updated mafById() function to speed up rs identifier look up.

- Updated scoring with weight matrices to enable using any type of
  weight matrix, and not only splice site matrices.

[wavClusteR](https://bioconductor.org/packages/wavClusteR)
----------

Changes in version 2.7.0:

- removed use of labeller (ggplot2) in annotateClusters

[xcms](https://bioconductor.org/packages/xcms)
----

Changes in version 1.49.7:

BUG FIXES

- Fix documentation warnings

Changes in version 1.49.6:

USER VISIBLE CHANGES

- Peak Picking function findPeaks.centWaveWithPredictedIsotopeROIs()
  and findPeaks.addPredictedIsotopeFeatures(), which allow more
  sensitive detection of isotope features.

Changes in version 1.49.5:

USER VISIBLE CHANGES

- Some documentation updates.

- Preparation for a new binning function

Changes in version 1.49.4:

BUG FIXES

- Fix getXcmsRaw that would prevent retention time correction to be
  applied (issue #44 reported by Aleksandr).

Changes in version 1.49.3:

NEW FEATURE

- updateObject method for xcmsSet.

USER VISIBLE CHANGES

- xcms uses now BiocParallel for parallel processing. All other
  parallel processing functions have been deprecated.

BUG FIXES

- Added missing package imports.

- Fix bug in fillPeaksChromPar referencing a non-existing variables i
  and object.

- Fix bug in group.nearest: variable scoreList was mis-spelled
  (coreList).

- Remove all DUP = FALSE from the .C calls as they are ignored anyways.

OTHER CHANGES

- Re-organization of class, function and method definitions in R-files.

- Use roxygen2 to manage the DESCRIPTION's collate field.

Changes in version 1.49.2:

NEW FEATURE

- Initial support for exporint mzTab format. Since Changes are still to
  be expected, xcms:::writeMzTab() is not yet exported.

Changes in version 1.49.1:

NEW FEATURE

- The raw CDF/mzXML/mzData/mzML is assumed to have scans sorted by m/z.
  Instead of throwing an "m/z sort assumption violated !" error, the
  data is re-read and on-demand sorted by m/z.

[xps](https://bioconductor.org/packages/xps)
---

Changes in version 3.2:

VERSION xps-1.29.1

- update README file

[YAPSA](https://bioconductor.org/packages/YAPSA)
-----

Changes in version 0.99.0:

- first version to be submitted to Bioconductor

Deprecated and Defunct Packages
===============================

1 software package (betr) was marked as deprecated, to be removed in the next release.

16 previously deprecated software packages were removed from the release.

