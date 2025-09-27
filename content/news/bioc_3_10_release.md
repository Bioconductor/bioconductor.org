October 30, 2019

Bioconductors:

We are pleased to announce Bioconductor 3.10, consisting of 1823
software packages, 384 experiment data packages, 953 annotation
packages, and 27 workflows.

There are 94 new software packages, 15 new data experiment packages,
3 new annotation packages, and many updates and improvements
to existing packages; Bioconductor 3.10 is compatible with R 3.6.1,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release will include an updated Bioconductor [Amazon Machine Image][1]
and [Docker containers][2].

Thank you to everyone for your contribution to Bioconductor

Visit [Bioconductor BiocViews][3]
for details and downloads.

[1]: /help/bioconductor-cloud-ami/
[2]: /help/docker/
[3]: /packages/release/BiocViews.html

Contents
--------

* [Getting Started with Bioconductor 3.10](#getting-started-with-bioconductor-35)
* [New Software Packages](#new-software-packages)
* [New Data Experiment Packages](#new-data-experiment-packages)
* [New Annotation Packages](#new-annotation-packages)
* [NEWS from new and existing software packages](#news-from-new-and-existing-software-packages)
* [NEWS from new and existing data experiment packages](#news-from-new-and-existing-data-experiment-packages)
* [NEWS from new and existing workflows](#news-from-new-and-existing-workflows)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.10
======================================

To update to or install Bioconductor 3.10:

1. Install R >=3.6.1.  Bioconductor 3.10 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [Installing Bioconductor](/install/).

New Software Packages
=====================

There are 93 new software packages in this release of Bioconductor.

- [AlphaBeta](/packages/AlphaBeta) AlphaBeta
  is a computational method for estimating epimutation rates and
  spectra from high-throughput DNA methylation data in plants. The
  method has been specifically designed to: 1. analyze 'germline'
  epimutations in the context of multi-generational mutation
  accumulation lines (MA-lines). 2. analyze 'somatic' epimutations in
  the context of plant development and aging.

- [ALPS](/packages/ALPS) The package provides
  analysis and publication quality visualization routines for
  genome-wide epigenomics data such as histone modification or
  transcription factor ChIP-seq, ATAC-seq, DNase-seq etc. The
  functions in the package can be used with any type of data that can
  be represented with bigwig files at any resolution. The goal of the
  ALPS is to provide analysis tools for most downstream analysis
  without leaving the R environment and most tools in the package
  require a minimal input that can be prepared with basic R, unix or
  excel skills.

- [APAlyzer](/packages/APAlyzer) Perform
  3'UTR APA, Intronic APA and gene expression analysis using RNA-seq
  data.

- [ASpediaFI](/packages/ASpediaFI) This
  package provides functionalities for a systematic and integrative
  analysis of alternative splicing events and their functional
  interactions.

- [Autotuner](/packages/Autotuner) This
  package is designed to help faciliate data processing in untargeted
  metabolomics. To do this, the algorithm contained within the
  package performs statistical inference on raw data to come up with
  the best set of parameters to process the raw data.

- [AWFisher](/packages/AWFisher)
  Implementation of the adaptively weighted fisher's method,
  including fast p-value computing, variability index, and
  meta-pattern.

- [BiocSet](/packages/BiocSet) BiocSet
  displays different biological sets in a triple tibble format. These
  three tibbles are `element`, `set`, and `elementset`. The user has
  the abilty to activate one of these three tibbles to perform common
  functions from the dplyr package. Mapping functionality and
  accessing web references for elements/sets are also available in
  BiocSet.

- [BioTIP](/packages/BioTIP) Adopting
  tipping-point theory to transcriptome profiles to unravel disease
  regulatory trajectory.

- [biscuiteer](/packages/biscuiteer) A test
  harness for bsseq loading of Biscuit output, summarization of WGBS
  data over defined regions and in mappable samples, with or without
  imputation, dropping of mostly-NA rows, age estimates, etc.

- [blacksheepr](/packages/blacksheepr)
  Blacksheep is a tool designed for outlier analysis in the context
  of pairwise comparisons in an effort to find distinguishing
  characteristics from two groups. This tool was designed to be
  applied for biological applications such as phosphoproteomics or
  transcriptomics, but it can be used for any data that can be
  represented by a 2D table, and has two sub populations within the
  table to compare.

- [brainflowprobes](/packages/brainflowprobes)
  Use these functions to characterize genomic regions for BrainFlow
  target probe design.

- [brendaDb](/packages/brendaDb) R interface
  for importing and analyzing enzyme information from the BRENDA
  database.

- [BUSpaRse](/packages/BUSpaRse) The kallisto
  | bustools pipeline is a fast and modular set of tools to convert
  single cell RNA-seq reads in fastq files into gene count or
  transcript compatibility counts (TCC) matrices for downstream
  analysis. Central to this pipeline is the barcode, UMI, and set
  (BUS) file format. This package serves the following purposes:
  First, this package allows users to manipulate BUS format files as
  data frames in R and then convert them into gene count or TCC
  matrices. Furthermore, since R and Rcpp code is easier to handle
  than pure C++ code, users are encouraged to tweak the source code
  of this package to experiment with new uses of BUS format and
  different ways to convert the BUS file into gene count matrix.
  Second, this package can conveniently generate files required to
  generate gene count matrices for spliced and unspliced transcripts
  for RNA velocity. Third, this package implements utility functions
  to get transcripts and associated genes required to convert BUS
  files to gene count matrices, to write the transcript to gene
  information in the format required by bustools, and to read output
  of bustools into R as sparses matrices.

- [calm](/packages/calm) Statistical methods
  for multiple testing with covariate information. Traditional
  multiple testing methods only consider a list of test statistics,
  such as p-values. Our methods incorporate the auxiliary
  information, such as the lengths of gene coding regions or the
  minor allele frequencies of SNPs, to improve power.

- [circRNAprofiler](/packages/circRNAprofiler)
  R-based computational framework for a comprehensive in silico
  analysis of circRNAs. This computational framework allows to
  combine and analyze circRNAs previously detected by multiple
  publicly available annotation-based circRNA detection tools. It
  covers different aspects of circRNAs analysis from differential
  expression analysis, evolutionary conservation, biogenesis to
  functional analysis.

- [cliqueMS](/packages/cliqueMS) Annotates data from
  liquid chromatography coupled to mass spectrometry (LC/MS) metabolomics
  experiments. Based on a network algorithm (O.Senan, A. Aguilar- Mogas,
  M. Navarro, O. Yanes, R.Guimerà and M. Sales-Pardo, Metabolomics Conference
  (2016), Dublin), 'CliqueMS' builds a weighted similarity network where nodes
  are features and edges are weighted according to the similarity of this
  features. Then it searches for the most plausible division of the similarity
  network into cliques (fully connected components). Finally it annotates
  metabolites within each clique, obtaining for each annotated metabolite the
  neutral mass and their features, corresponding to isotopes, ionization adducts
  and fragmentation adducts of that metabolite.

- [CNVfilteR](/packages/CNVfilteR) CNVfilteR
  identifies those CNVs that can be discarded by using the single
  nucleotide variant (SNV) calls that are usually obtained in common
  NGS pipelines.

- [CrossICC](/packages/CrossICC) CrossICC
  utilizes an iterative strategy to derive the optimal gene set and
  cluster number from consensus similarity matrix generated by
  consensus clustering and it is able to deal with multiple cross
  platform datasets so that requires no between-dataset
  normalizations. This package also provides abundant functions for
  visualization and identifying subtypes of cancer. Specially, many
  cancer-related analysis methods are embedded to facilitate the
  clinical translation of the identified cancer subtypes.

- [debCAM](/packages/debCAM) An R package for
  fully unsupervised deconvolution of complex tissues. It provides
  basic functions to perform unsupervised deconvolution on mixture
  expression profiles by Convex Analysis of Mixtures (CAM) and some
  auxiliary functions to help understand the subpopulation-specific
  results. It also implements functions to perform supervised
  deconvolution based on prior knowledge of molecular markers, S
  matrix or A matrix. Combining molecular markers from CAM and from
  prior knowledge can achieve semi-supervised deconvolution of
  mixtures.

- [deltaCaptureC](/packages/deltaCaptureC)
  This package discovers meso-scale chromatin remodelling from 3C
  data.  3C data is local in nature.  It givens interaction counts
  between restriction enzyme digestion fragments and a preferred
  'viewpoint' region.  By binning this data and using permutation
  testing, this package can test whether there are statistically
  significant changes in the interaction counts between the data from
  two cell types or two treatments.

- [DEWSeq](/packages/DEWSeq) Differential
  expression analysis of windows for next-generation sequencing data
  like eCLIP or iCLIP data.

- [DMCFB](/packages/DMCFB) DMCFB is a
  pipeline for identifying differentially methylated cytosines using
  a Bayesian functional regression model in bisulfite sequencing
  data. By using a functional regression data model, it tries to
  capture position-specific, group-specific and other
  covariates-specific methylation patterns as well as spatial
  correlation patterns and unknown underlying models of methylation
  data. It is robust and flexible with respect to the true underlying
  models and inclusion of any covariates, and the missing values are
  imputed using spatial correlation between positions and samples. A
  Bayesian approach is adopted for estimation and inference in the
  proposed method.

- [fcoex](/packages/fcoex) The fcoex package
  implements an easy-to use interface to co-expression analysisbased
  on the FCBF (Fast Correlation-Based Filter) algorithm. it was
  implemented especifically to deal with single-cell data. The
  modules found can be used to redefine cell populations, unrevel
  novel gene associations and predict gene function by
  guilt-by-association. The package structure is based on the
  CEMiTool package.

- [fcScan](/packages/fcScan) This package is
  used to detect combination of genomic coordinates falling within a
  user defined window size along with user defined overlap between
  identified neighboring clusters. It can be used for genomic data
  where the clusters are built on a specific chromosome or specific
  strand. Clustering can be performed with a "greedy" option allowing
  thus the presence of additional sites within the allowed window
  size.

- [flowSpecs](/packages/flowSpecs) This
  package is intended to fill the role of conventional cytometry
  pre-processing software, for spectral decomposition,
  transformation, visualization and cleanup, and to aid further
  downstream analyses, such as with DepecheR, by enabling
  transformation of flowFrames and flowSets to dataframes. Functions
  for flowCore-compliant automatic 1D-gating/filtering are in the
  pipe line. The package name has been chosen both as it will deal
  with spectral cytometry and as it will hopefully give the user a
  nice pair of spectacles through which to view their data.

- [flowSpy](/packages/flowSpy) A trajectory
  inference and visualization toolkit for flow and mass cytometry
  data. flowSpy offers complete analyzing workflow for flow and mass
  cytometry data. flowSpy can be a valuable tool for application
  ranging from clustering and dimensionality reduction to trajectory
  reconstruction and pseudotime estimation for flow and mass
  cytometry data.

- [GCSscore](/packages/GCSscore) For
  differential expression analysis of 3'IVT and WT-style microarrays
  from Affymetrix/Thermo-Fisher.  Based on S-score algorithm
  originally described by Zhang et al 2002.

- [gemini](/packages/gemini) GEMINI uses
  log-fold changes to model sample-dependent and independent effects,
  and uses a variational Bayes approach to infer these effects. The
  inferred effects are used to score and identify genetic
  interactions, such as lethality and recovery. More details can be
  found in Zamanighomi et al. 2019 (in press).

- [GenomicOZone](/packages/GenomicOZone) The
  package clusters gene activity along chromosome into zones, detects
  differential zones as outstanding, and visualizes maps of
  outstanding zones across the genome. The method guarantees cluster
  optimality, linear runtime to sample size, and reproducibility. It
  enables new characterization of effects due to genome
  reorganization, structural variation, and epigenome alteration.

- [GmicR](/packages/GmicR) This package uses
  bayesian network learning to detect relationships between Gene
  Modules detected by WGCNA and immune cell signatures defined by
  xCell. It is a hypothesis generating tool.

- [gramm4R](/packages/gramm4R) Generalized
  Correlation Analysis for Metabolome and Microbiome (GRaMM), for
  inter-correlation pairs discovery among metabolome and microbiome.

- [gscreend](/packages/gscreend) Package for
  the analysis of pooled genetic screens (e.g. CRISPR-KO). The
  analysis of such screens is based on the comparison of gRNA
  abundances before and after a cell proliferation phase. The
  gscreend packages takes gRNA counts as input and allows detection
  of genes whose knockout decreases or increases cell proliferation.

- [HCAExplorer](/packages/HCAExplorer)
  Search, browse, reference, and download resources from the Human
  Cell Atlas data portal. Development of this package is supported
  through funds from the Chan / Zuckerberg initiative.

- [HiLDA](/packages/HiLDA) A package built
  under the Bayesian framework of applying hierarchical latent
  Dirichlet allocation to statistically test whether the mutational
  exposures of mutational signatures (Shiraishi-model signatures) are
  different between two groups.

- [idr2d](/packages/idr2d) A tool to measure
  reproducibility between genomic experiments that produce
  two-dimensional peaks (interactions between peaks), such as
  ChIA-PET, HiChIP, and HiC. idr2d is an extension of the original
  idr package, which is intended for (one-dimensional) ChIP-seq
  peaks.

- [IgGeneUsage](/packages/IgGeneUsage)
  Decoding the properties of immune repertoires is key in
  understanding the response of adaptive immunity to challenges such
  as viral infection. One important task in immune repertoire
  profiling is the detection of biases in Ig gene usage between
  biological conditions. IgGeneUsage is a computational tool for the
  analysis of differential gene usage in immune repertoires. It
  employs Bayesian hierarchical models to fit complex gene usage data
  from immune repertoire sequencing experiments and quantifies Ig
  gene usage biases as probabilities.

- [KnowSeq](/packages/KnowSeq) KnowSeq
  proposes a whole pipeline that comprises the most relevant steps in
  the RNA-seq gene expression analysis, with the main goal of
  extracting biological knowledge from raw data (Differential
  Expressed Genes, Gene Ontology enrichment, pathway visualization
  and related diseases). In this sense, KnowSeq allows aligning raw
  data from the original fastq or sra files, by using the most
  renowned aligners such as tophat2, hisat2, salmon and kallisto.
  Nowadays, there is no package that only from the information of the
  samples to align -included in a text file-, automatically performs
  the download and alignment of all of the samples. Furthermore, the
  package includes functions to: calculate the gene expression
  values; remove batch effect; calculate the Differentially Expressed
  Genes (DEGs); plot different graphs; and perform the DEGs
  enrichment with the GO information, pathways visualization and
  related diseases information retrieval. Moreover, KnowSeq is the
  only package that allows applying both a machine learning and DEGs
  enrichment processes just after the DEGs extraction. This idea
  emerged with the aim of proposing a complete tool to the research
  community containing all the necessary steps to carry out complete
  studies in a simple and fast way.

- [LinkHD](/packages/LinkHD) Here we present
  Link-HD, an approach to integrate heterogeneous datasets, as a
  generalization of STATIS-ACT (“Structuration des Tableaux A Trois
  Indices de la Statistique–Analyse Conjointe de Tableaux”), a family
  of methods to join and compare information from multiple subspaces.
  However, STATIS-ACT has some drawbacks since it only allows
  continuous data and it is unable to establish relationships between
  samples and features. In order to tackle these constraints, we
  incorporate multiple distance options and a linear regression based
  Biplot model in order to stablish relationships between
  observations and variable and perform variable selection.

- [lionessR](/packages/lionessR) LIONESS, or
  Linear Interpolation to Obtain Network Estimates for Single
  Samples, can be used to reconstruct single-sample networks
  (https://arxiv.org/abs/1505.06440). This code implements the
  LIONESS equation in the lioness function in R to reconstruct
  single-sample networks. The default network reconstruction method
  we use is based on Pearson correlation. However, lionessR can run
  on any network reconstruction algorithms that returns a complete,
  weighted adjacency matrix. lionessR works for both unipartite and
  bipartite networks.

- [Maaslin2](/packages/Maaslin2) MaAsLin2 is
  comprehensive R package for efficiently determining multivariable
  association between clinical metadata and microbial meta'omic
  features. MaAsLin2 relies on general linear models to accommodate
  most modern epidemiological study designs, including
  cross-sectional and longitudinal, and offers a variety of data
  exploration, normalization, and transformation methods. MaAsLin2 is
  the next generation of MaAsLin.

- [MACSQuantifyR](/packages/MACSQuantifyR)
  Automatically process the metadata of MACSQuantify FACS sorter. It
  runs multiple modules: i) imports of raw file and graphical
  selection of duplicates in well plate, ii) computates statistics on
  data and iii) can compute combination index.

- [MBQN](/packages/MBQN) Modified quantile
  normalization for omics or other matrix-like data distorted in
  location and scale.

- [MEB](/packages/MEB) Identifying
  differential expression genes between the same or different species
  is an urgent demand for biological research. In most of cases,
  normalization is the first step to solve this problem, then by
  employing the hypothesis testing, we could detect statistically
  significant genes. With the development of machine learning, it
  gives us a new perspective on discrimination between differential
  expression (DE) and non-differential expression (non-DE) genes.
  Provided a set of training data, the procedure of distinguishing
  genes could be simplified as a classification problem. However, in
  reality, it is hard for us to get the information from both DE and
  non-DE genes. To solve this problem, we try to identify
  differential cases only in the domain of non-DE genes, and
  transform the problem to an outlier detection in machine learning.
  Given that non-DE genes have some similarities in features, we
  build a Minimum Enclosing Ball (MEB) to cover those non-DE genes in
  feature space, then those DE genes, which are enormously different
  from non-DE genes, being regarded as outliers and rejected outside
  the ball. Compared with existing methods, it is no need for the MEB
  method to normalize data in advance. Besides, the MEB method could
  be easily applied to the same or different species data and without
  changing too much.

- [MetaVolcanoR](/packages/MetaVolcanoR) This
  package combines differential gene expression results. It
  implements three strategies to summarize differential gene
  expression from different studies. i) Random Effects Model (REM)
  approach, ii) a p-value combining-approach, and iii) a naive
  vote-counting approach. In all cases, MetaVolcano exploits the
  Volcano plot reasoning to visualize the gene expression
  meta-analysis results.

- [MethCP](/packages/MethCP) MethCP is a
  differentially methylated region (DMR) detecting method for
  whole-genome bisulfite sequencing (WGBS) data, which is applicable
  for a wide range of experimental designs beyond the two-group
  comparisons, such as time-course data. MethCP identifies DMRs based
  on change point detection, which naturally segments the genome and
  provides region-level differential analysis.

- [methrix](/packages/methrix) Bedgraph files
  generated by Bisulfite pipelines often come in various flavors.
  Critical downstream step requires summarization of these files into
  methylation/coverage matrices. This step of data aggregation is
  done by Methrix, including many other useful downstream functions.

- [methylCC](/packages/methylCC) A tool to
  estimate the cell composition of DNA methylation whole blood sample
  measured on any platform technology (microarray and sequencing).

- [microbiomeDASim](/packages/microbiomeDASim)
  A toolkit for simulating differential microbiome data designed for
  longitudinal analyses. Several functional forms may be specified
  for the mean trend. Observations are drawn from a multivariate
  normal model. The objective of this package is to be able to
  simulate data in order to accurately compare different longitudinal
  methods for differential abundance.

- [MMAPPR2](/packages/MMAPPR2) MMAPPR2 maps
  mutations resulting from pooled RNA-seq data from the F2 cross of
  forward genetic screens. Its predecessor is described in a paper
  published in Genome Research (Hill et al. 2013). MMAPPR2 accepts
  aligned BAM files as well as a reference genome as input,
  identifies loci of high sequence disparity between the control and
  mutant RNA sequences, predicts variant effects using Ensembl's
  Variant Effect Predictor, and outputs a ranked list of candidate
  mutations.

- [MMUPHin](/packages/MMUPHin) MMUPHin is an
  R package for meta-analysis tasks of microbiome cohorts. It has
  function interfaces for: a) covariate-controlled batch- and cohort
  effect adjustment, b) meta-analysis differential abundance testing,
  c) meta-analysis unsupervised discrete structure (clustering)
  discovery, and d) meta-analysis unsupervised continuous structure
  discovery.

- [MOSim](/packages/MOSim) MOSim package
  simulates multi-omic experiments that mimic regulatory mechanisms
  within the cell, allowing flexible experimental design including
  time course and multiple groups.

- [MSstatsSampleSize](/packages/MSstatsSampleSize)
  The packages estimates the variance in the input protein abundance
  data and simulates data with predefined number of biological
  replicates based on the variance estimation. It reports the mean
  predictive accuracy of the classifier and mean protein importance
  over multiple iterations of the simulation.

- [muscat](/packages/muscat) `muscat`
  provides various methods and visualization tools for DS analysis in
  multi-sample, multi-group, multi-(cell-)subpopulation scRNA-seq
  data, including cell-level mixed models and methods based on
  aggregated “pseudobulk” data, as well as a flexible simulation
  platform that mimics both single and multi-sample scRNA-seq data.

- [ncGTW](/packages/ncGTW) The purpose of
  ncGTW is to help XCMS for LC-MS data alignment. Currently, ncGTW
  can detect the misaligned feature groups by XCMS, and the user can
  choose to realign these feature groups by ncGTW or not.

- [OmnipathR](/packages/OmnipathR) Import
  data from https://www.omnipathdb.org webservice. It also includes
  functions to transform and print this data.

- [oppti](/packages/oppti) The aim of oppti
  is to analyze protein (and phosphosite) expressions to find
  outlying markers for each sample in the given cohort(s) for the
  discovery of personalized actionable targets.

- [peakPantheR](/packages/peakPantheR) An
  automated pipeline for the detection, integration and reporting of
  predefined features across a large number of mass spectrometry data
  files.

- [PERFect](/packages/PERFect) PERFect is a
  novel permutation filtering approach designed to address two
  unsolved problems in microbiome data processing: (i) define and
  quantify loss due to filtering by implementing thresholds, and (ii)
  introduce and evaluate a permutation test for filtering loss to
  provide a measure of excessive filtering.

- [PhyloProfile](/packages/PhyloProfile)
  PhyloProfile is a tool for exploring complex phylogenetic profiles.
  Phylogenetic profiles, presence/absence patterns of genes over a
  set of species, are commonly used to trace the functional and
  evolutionary history of genes across species and time. With
  PhyloProfile we can enrich regular phylogenetic profiles with
  further data like sequence/structure similarity, to make
  phylogenetic profiling more meaningful. Besides the interactive
  visualisation powered by R-Shiny, the package offers a set of
  further analysis features to gain insights like the gene age
  estimation or core gene identification.

- [proDA](/packages/proDA) Account for
  missing values in label-free mass spectrometry data without
  imputation. The package implements a probabilistic dropout model
  that ensures that the information from observed and missing values
  are properly combined. It adds empirical Bayesian priors to
  increase power to detect differentially abundant proteins.

- [pulsedSilac](/packages/pulsedSilac) This
  package provides several tools for pulsed-SILAC data analysis.
  Functions are provided to organize the data, calculate isotope
  ratios, isotope fractions, model protein turnover, compare turnover
  models, estimate cell growth and estimate isotope recycling.
  Several visualization tools are also included to do basic data
  exploration, quality control, condition comparison, individual
  model inspection and model comparison.

- [pwrEWAS](/packages/pwrEWAS) pwrEWAS is a
  user-friendly tool to assists researchers in the design and
  planning of EWAS to help circumvent under- and overpowered studies.

- [Qtlizer](/packages/Qtlizer) This R package
  provides access to the Qtlizer web server. Qtlizer annotates lists
  of common small variants (mainly SNPs) and genes in humans with
  associated changes in gene expression using the most comprehensive
  database of published quantitative trait loci (QTLs).

- [ReactomeGSA](/packages/ReactomeGSA) The
  ReactomeGSA packages uses Reactome's online analysis service to
  perform a multi-omics gene set analysis. The main advantage of this
  package is, that the retrieved results can be visualized using
  REACTOME's powerful webapplication. Since Reactome's analysis
  service also uses R to perfrom the actual gene set analysis you
  will get similar results when using the same packages (such as
  limma and edgeR) locally. Therefore, if you only require a gene set
  analysis, different packages are more suited.

- [RNAmodR](/packages/RNAmodR) RNAmodR
  provides classes and workflows for loading/aggregation data from
  high througput sequencing aimed at detecting post-transcriptional
  modifications through analysis of specific patterns. In addition,
  utilities are provided to validate and visualize the results. The
  RNAmodR package provides a core functionality from which specific
  analysis strategies can be easily implemented as a seperate
  package.

- [RNAmodR.AlkAnilineSeq](/packages/RNAmodR.AlkAnilineSeq)
  RNAmodR.AlkAnilineSeq implements the detection of m7G, m3C and D
  modifications on RNA from experimental data generated with the
  AlkAnilineSeq protocol. The package builds on the core
  functionality of the RNAmodR package to detect specific patterns of
  the modifications in high throughput sequencing data.

- [RNAmodR.ML](/packages/RNAmodR.ML)
  RNAmodR.ML extend the functionality of the RNAmodR package and
  classical detection strategies towards detection through machine
  learning models. RNAmodR.ML provides classes, functions and an
  example workflow to establish a detection stratedy, which can be
  packaged.

- [RNAmodR.RiboMethSeq](/packages/RNAmodR.RiboMethSeq)
  RNAmodR.RiboMethSeq implements the detection of 2'-O methylations
  on RNA from experimental data generated with the RiboMethSeq
  protocol. The package builds on the core functionality of the
  RNAmodR package to detect specific patterns of the modifications in
  high throughput sequencing data.

- [RNAsense](/packages/RNAsense) RNA-sense
  tool compares RNA-seq time curves in two experimental conditions,
  i.e. wild-type and mutant, and works in three steps. At Step 1, it
  builds expression profile for each transcript in one condition
  (i.e. wild-type) and tests if the transcript abundance grows or
  decays significantly.  Dynamic transcripts are then sorted to
  non-overlapping groups (time profiles) by the time point of switch
  up or down. At Step 2, RNA-sense outputs the groups of
  differentially expressed transcripts, which are up- or
  downregulated in the mutant compared to the wild-type at each time
  point. At Step 3, Correlations (Fisher's exact test) between the
  outputs of Step 1 (switch up- and switch down- time profile groups)
  and the outputs of Step2 (differentially expressed transcript
  groups) are calculated. The results of the correlation analysis are
  printed as two-dimensional color plot, with time profiles and
  differential expression groups at y- and x-axis, respectively, and
  facilitates the biological interpretation of the data.

- [SAIGEgds](/packages/SAIGEgds) Scalable
  implementation of generalized mixed models with highly optimized
  C++ implementation and integration with Genomic Data Structure
  (GDS) files. It is designed for single variant tests in large-scale
  phenome-wide association studies (PheWAS) with millions of variants
  and samples, controlling for sample structure and case-control
  imbalance. The implementation is based on the original SAIGE R
  package (v0.29.4.4). SAIGEgds also implements some of the SPAtest
  functions in C to speed up the calculation of Saddlepoint
  approximation. Benchmarks show that SAIGEgds is 5 to 6 times faster
  than the original SAIGE R package.

- [SBGNview](/packages/SBGNview) SBGNview is
  an R package for visualizing omics data on SBGN pathway maps. Given
  omics data and a SBGN-ML file with layout information, SBGNview can
  display omics data as colors on glyphs and output image files.
  SBGNview provides extensive options to control glyph and edge
  features (e.g. color, line width etc.). To facilitate pathway based
  analysis, SBGNview also provides functions to extract molecule sets
  from SBGN-ML files. SBGNview can map a large collection of gene,
  protein and compound ID typs to glyphs.

- [SCANVIS](/packages/SCANVIS) SCANVIS is a
  set of annotation-dependent tools for analyzing splice junctions
  and their read support as predetermined by an alignment tool of
  choice (for example, STAR aligner). SCANVIS assesses each
  junction's relative read support (RRS) by relating to the context
  of local split reads aligning to annotated transcripts. SCANVIS
  also annotates each splice junction by indicating whether the
  junction is supported by annotation or not, and if not, what type
  of junction it is (e.g. exon skipping, alternative 5' or 3' events,
  Novel Exons). Unannotated junctions are also futher annotated by
  indicating whether it induces a frame shift or not. SCANVIS
  includes a visualization function to generate static sashimi-style
  plots depicting relative read support and number of split reads
  using arc thickness and arc heights, making it easy for users to
  spot well-supported junctions. These plots also clearly delineate
  unannotated junctions from annotated ones using designated color
  schemes, and users can also highlight splice junctions of choice.
  Variants and/or a read profile are also incoroporated into the plot
  if the user supplies variants in bed format and/or the BAM file.
  One further feature of the visualization function is that users can
  submit multiple samples of a certain disease or cohort to generate
  a single plot - this occurs via a "merge" function wherein junction
  details over multiple samples are merged to generate a single
  sashimi plot, which is useful when contrasting cohorots (eg.
  disease vs control).

- [scBFA](/packages/scBFA) This package is
  designed to model gene detection pattern of scRNA-seq through a
  binary factor analysis model. This model allows user to pass into a
  cell level covariate matrix X and gene level covariate matrix Q to
  account for nuisance variance(e.g batch effect), and it will output
  a low dimensional embedding matrix for downstream analysis.

- [scDblFinder](/packages/scDblFinder)
  Efficient identification of doublets in single-cell RNAseq directly
  from counts using overclustering-based generation of artifical
  doublets.

- [scGPS](/packages/scGPS) The package
  implements two main algorithms to answer two key questions: a SCORE
  (Stable Clustering at Optimal REsolution) to find subpopulations,
  followed by scGPS to investigate the relationships between
  subpopulations.

- [schex](/packages/schex) Builds hexbin
  plots for variables and dimension reduction stored in single cell
  omics data such as SingleCellExperiment and SeuratObject. The ideas
  used in this package are based on the excellent work of Dan Carr,
  Nicholas Lewin-Koh, Martin Maechler and Thomas Lumley.

- [scPCA](/packages/scPCA) A toolbox for
  sparse contrastive principal component analysis (scPCA) of
  high-dimensional biological data. scPCA combines the stability and
  interpretability of sparse PCA with contrastive PCA's ability to
  disentangle biological signal from techical noise through the use
  of control data. Also implements and extends cPCA.

- [scTGIF](/packages/scTGIF) scTGIF connects
  the cells and the related gene functions without cell type label.

- [SEtools](/packages/SEtools) This includes
  a set of tools for working with the SummarizedExperiment class,
  including handy merging and plotting functions.

- [SharedObject](/packages/SharedObject) This
  package is developed for facilitating parallel computing in R. It
  is capable to create an R object in the shared memory space and
  share the data across multiple R processes. It avoids the overhead
  of memory dulplication and data transfer, which make sharing big
  data object across many clusters possible.

- [signatureSearch](/packages/signatureSearch)
  This package implements algorithms and data structures for
  performing gene expression signature (GES) searches, and
  subsequently interpreting the results functionally with specialized
  enrichment methods.

- [SigsPack](/packages/SigsPack) Single
  sample estimation of exposure to mutational signatures. Exposures
  to known mutational signatures are estimated for single samples,
  based on quadratic programming algorithms. Bootstrapping the input
  mutational catalogues provides estimations on the stability of
  these exposures. The effect of the sequence composition of
  mutational context can be taken into account by normalising the
  catalogues.

- [SingleR](/packages/SingleR) Performs
  unbiased cell type recognition from single-cell RNA sequencing
  data, by leveraging reference transcriptomic datasets of pure cell
  types to infer the cell of origin of each single cell
  independently.

- [sojourner](/packages/sojourner) Single
  molecule tracking has evolved as a novel new approach complementing
  genomic sequencing, it reports live biophysical properties of
  molecules being investigated besides properties relating their
  coding sequence; here we provided "sojourner" package, to address
  statistical and bioinformatic needs related to the analysis and
  comprehension of high throughput single molecule tracking data.

- [Spaniel](/packages/Spaniel) Spaniel
  includes a series of tools to aid the quality control and analysis
  of Spatial Transcriptomics data. The package contains functions to
  create either a Seurat object or SingleCellExperiment from a count
  matrix and spatial barcode file and provides a method of loading a
  histologial image into R. The spanielPlot function allows
  visualisation of metrics contained within the S4 object overlaid
  onto the image of the tissue.

- [SQLDataFrame](/packages/SQLDataFrame)
  SQLDataFrame is developed to lazily represent and efficiently
  analyze SQL-based tables in _R_. SQLDataFrame supports common and
  familiar 'DataFrame' operations such as '[' subsetting, rbind,
  cbind, etc.. The internal implementation is based on the widely
  adopted dplyr grammar and SQL commands. In-memory datasets or plain
  text files (.txt, .csv, etc.) could also be easily converted into
  SQLDataFrames objects (which generates a new database on-disk).

- [ssPATHS](/packages/ssPATHS) This package
  generates pathway scores from expression data for single samples
  after training on a reference cohort. The score is generated by
  taking the expression of a gene set (pathway) from a reference
  cohort and performing linear discriminant analysis to distinguish
  samples in the cohort that have the pathway augmented and not. The
  separating hyperplane is then used to score new samples.

- [target](/packages/target) Implement the
  BETA algorithm for infering direct target genes from DNA-binding
  and perturbation expression data Wang et al. (2013) <doi:
  10.1038/nprot.2013.150>. Extend the algorithm to predict the
  combined function of two DNA-binding elements from comprable
  binding and expression data.

- [TOAST](/packages/TOAST) This package is
  devoted to analyzing high-throughput data (e.g. gene expression
  microarray, DNA methylation microarray, RNA-seq) from complex
  tissues. Current functionalities include 1. detect cell-type
  specific or cross-cell type differential signals 2. improve
  variable selection in reference-free deconvolution.

- [tradeSeq](/packages/tradeSeq) tradeSeq
  provides a flexible method for finding genes that are
  differentially expressed along one or multiple trajectories, using
  a variety of tests suited to answer questions of interest, e.g. the
  discovery of genes that whose expression is associated with
  pseudotime, or who are differentially expressed (in a specific
  region) along the trajectory. It fits a generalized additive model
  (GAM) for each gene, and performs inference on the parameters of
  the GAM.

- [VariantExperiment](/packages/VariantExperiment)
  VariantExperiment is a Bioconductor package for saving data in
  VCF/GDS format into RangedSummarizedExperiment object. The
  high-throughput genetic/genomic data are saved in GDSArray objects.
  The annotation data for features/samples are saved in
  DelayedDataFrame format with mono-dimensional GDSArray in each
  column. The on-disk representation of both assay data and
  annotation data achieves on-disk reading and processing and saves
  memory space significantly. The interface of
  RangedSummarizedExperiment data format enables easy and common
  manipulations for high-throughput genetic/genomic data with common
  SummarizedExperiment metaphor in R and Bioconductor.

- [ViSEAGO](/packages/ViSEAGO) The main
  objective of ViSEAGO package is to carry out a data mining of
  biological functions and establish links between genes involved in
  the study. We developed ViSEAGO in R to facilitate functional Gene
  Ontology (GO) analysis of complex experimental design with multiple
  comparisons of interest. It allows to study large-scale datasets
  together and visualize GO profiles to capture biological knowledge.
  The acronym stands for three major concepts of the analysis:
  Visualization, Semantic similarity and Enrichment Analysis of Gene
  Ontology. It provides access to the last current GO annotations,
  which are retrieved from one of NCBI EntrezGene, Ensembl or Uniprot
  databases for several species. Using available R packages and novel
  developments, ViSEAGO extends classical functional GO analysis to
  focus on functional coherence by aggregating closely related
  biological themes while studying multiple datasets at once. It
  provides both a synthetic and detailed view using interactive
  functionalities respecting the GO graph structure and ensuring
  functional coherence supplied by semantic similarity. ViSEAGO has
  been successfully applied on several datasets from different
  species with a variety of biological questions. Results can be
  easily shared between bioinformaticians and biologists, enhancing
  reporting capabilities while maintaining reproducibility.

- [waddR](/packages/waddR) Wasserstein
  distance based statistical test for detecting and describing
  differential distributions in one-dimensional data. Functions for
  wasserstein distance calculation, differential distribution
  testing, and a specialized test for differential expression in
  scRNA data are provided.

- [XCIR](/packages/XCIR) Models and tools for
  subject level analysis of X chromosome inactivation (XCI) and
  XCI-escape inference.

New Data Experiment Packages
=====================

There are 15 new data experiment packages in this release of Bioconductor.

- [benchmarkfdrData2019](/packages/benchmarkfdrData2019)
  Benchmarking results for experimental and simulated data sets used
  in Korthauer and Kimes et al. (2019) to compare methods for
  controlling the false discovery rate.

- [biscuiteerData](/packages/biscuiteerData)
  Contains default datasets used by the Bioconductor package
  biscuiteer.

- [depmap](/packages/depmap) The depmap
  package is a data package that accesses datsets from the Broad
  Institute DepMap cancer dependency study using ExperimentHub.
  Datasets from the most current release are available, including
  RNAI and CRISPR-Cas9 gene knockout screens quantifying the genetic
  dependency for select cancer cell lines. Additional datasets are
  also available pertaining to the log copy number of genes for
  select cell lines, protein expression of cell lines as measured by
  reverse phase protein lysate microarray (RPPA), 'Transcript Per
  Million' (TPM) data, as well as supplementary datasets which
  contain metadata and mutation calls for the other datasets found in
  the current release. The 19Q3 release adds the drug_dependency
  dataset, that contains cancer cell line dependency data with
  respect to drug and drug-candidate compounds. This package will be
  updated on a quarterly basis to incorporate the latest Broad
  Institute DepMap Public cancer dependency datasets. All data made
  available in this package was generated by the Broad Institute
  DepMap for research purposes and not intended for clinical use.
  This data is distributed under the Creative Commons license
  (Attribution 4.0 International (CC BY 4.0)).

- [HMP2Data](/packages/HMP2Data) HMP2Data is
  a Bioconductor package of the Human Microbiome Project 2 (HMP2) 16S
  rRNA sequencing data. Processed data is provided as phyloseq,
  SummarizedExperiment, and MultiAssayExperiment class objects.
  Individual matrices and data.frames used for building these S4
  class objects are also provided in the package.

- [MMAPPR2data](/packages/MMAPPR2data)
  Contains data for illustration purposes in the MMAPPR2 package,
  namely simulated BAM files containing RNA-Seq data for a mutation
  in the slc24a5 gene, taken from the GRCz11 genome. Also contains
  reference sequence and annotation files for the region.

- [MouseGastrulationData](/packages/MouseGastrulationData)
  Provides processed and raw count matrices for single-cell RNA
  sequencing data from a timecourse of mouse gastrulation and early
  organogenesis.

- [muscData](/packages/muscData) Data package
  containing a collection of multi-sample multi-group scRNA-seq
  datasets in SingleCellExperiment Bioconductor object format.

- [PhyloProfileData](/packages/PhyloProfileData)
  Two experimental datasets to illustrate running and analysing
  phylogenetic profiles with PhyloProfile package.

- [pwrEWAS.data](/packages/pwrEWAS.data) This
  package provides reference data required for pwrEWAS. pwrEWAS is a
  user-friendly tool to estimate power in EWAS as a function of
  sample and effect size for two-group comparisons of DNAm (e.g.,
  case vs control, exposed vs non-exposed, etc.).

- [ReactomeGSA.data](/packages/ReactomeGSA.data)
  Companion data sets to showcase the functionality of the
  ReactomeGSA package. This package contains proteomics and RNA-seq
  data of the melanoma B-cell induction study by Griss et al.

- [RNAmodR.Data](/packages/RNAmodR.Data)
  RNAmodR.Data contains example data, which is used for vignettes and
  example workflows in the RNAmodR and dependent packages.

- [SBGNview.data](/packages/SBGNview.data)
  This package contains: 1. A microarray gene expression dataset from
  a human breast cancer study. 2. A RNA-Seq gene expression dataset
  from a mouse study on IFNG knockout. 3. ID mapping tables between
  gene IDs and SBGN-ML file glyph IDs. 4. Percent of orthologs
  detected in other species of the genes in a pathway. Cutoffs of
  this percentage for defining if a pathway exists in another
  species. 5. XML text of SBGN-ML files for all pre-collected
  pathways.

- [signatureSearchData](/packages/signatureSearchData)
  CMAP/LINCS hdf5 databases and other annotations used for
  signatureSearch software package.

- [tartare](/packages/tartare) provides raw
  files (size=278MBytes) recorded on different Liquid Chromatography
  Mass Spectrometry (LC-MS) instruments. All included MS instruments
  are manufactured by Thermo Fisher Scientific and belong to the
  Orbitrap Tribrid or Q Exactive Orbitrap family of instruments.
  Despite their common origin and shared hardware components (e.g.
  Orbitrap mass analyser), the above instruments tend to write data
  in different "dialects" in a shared binary file format (.raw). The
  intention behind tartare is to provide complex but slim real-world
  files that can be used to make code robust with respect to this
  diversity. In other words, it is intended for enhanced unit
  testing. The package is considered to be used with the rawDiag
  package (Trachsel, 2018 <doi:10.1021/acs.jproteome.8b00173>) and
  the Spectra MsBackends.

- [TENxBUSData](/packages/TENxBUSData)
  Download Barcode, UMI, and Set (BUS) format of 10x datasets from
  within R. This package accompanies the package BUSpaRse, which can
  load BUS format into R as a sparse matrix, and which has utility
  functions related to using the C++ command line package bustools.


New Annotation Packages
=====================

There are 3 new annotation packages in this release of Bioconductor.

- [MafDb.gnomAD.r3.0.GRCh38](/packages/MafDb.gnomAD.r3.0.GRCh38)
- [org.Mxanthus.db](/packages/org.Mxanthus.db)
- [GenomicState](/packages/GenomicState)


NEWS from new and existing Software Packages
===================================


[AffiXcan](/packages/AffiXcan)
--------

Changes in version 1.3.7

UPDATED VIGNETTE

- Package vignette has been updated to address the new features

ALSO MODELS' P-VALUES AND RHO ARE RETURNED IN CROSS-VALIDATION

- For each fold, for each gene for which a GReX could be imputed, the
  following values are now returned

- Squared correlation of the model's predictions with training data

- P-value of the model

- Corrected p-value of the model with benjamini-hochberg procedure

- Correlation of the model's predictions with validation data

- Squared correlation of the model's predictions with validation data

- P-value of the correlation test of the model's predictions with
  validation data

Changes in version 1.3.5

P-VALUE OF THE CORRELATION TEST

- P-value of the cor.test() between predicted GReX and real expression
  values of genes is returned in when performing cross-validation

- When using affiXcanTrain in cross-validation mode, three main values
  for each gene are therefore returned: rho and rho squared (see
  changes in v 1.3.1), and the p-value of cor.test()

Changes in version 1.3.3

UPDATED DOCUMENTATION

- Formatting of functions documentation has been improved

- Important: vignette is still outdated (AffiXcan 1.2.0)

Changes in version 1.3.2

POPULATION STRUCTURE COVARIATES ARE OPTIONAL

- Providing population structure covariates is now not mandatory to
  perform models training

Changes in version 1.3.1

K-FOLD CROSS-VALIDATION IS SUPPORTED

- ANOVA p-value < 0.05 to assess prediction significance is not used
  anymore; instead:

- A k-fold cross-validation on the training dataset may be performed; k
  can be defined by the user

- Pearson correlation coefficients (R) and determination coefficients
  (R^2) between predicted GReX and real expression values of genes are
  computed

- In literature, GReX of genes for which the mean of the R^2 is above
  0.01 are considered non-randomly predicted, according to the new
  benchmarking standards &#91;ref&#93;

[alevinQC](/packages/alevinQC)
--------

Changes in version 1.1

- Added ability to process output from Salmon v0.14 or later

[AlphaBeta](/packages/AlphaBeta)
---------

Changes in version 0.99.2 (2019-08-15)

- Fixed bugs.

- Made the following significant changes o Using BiocParallel for
  parallel evaulation. o Updated documentation file.

[ALPS](/packages/ALPS)
----

Changes in version 0.99.8

- Bioc review fixes

Changes in version 0.99.7

- Bioc review fixes

Changes in version 0.99.0

- FIRST RELEASE

[AneuFinder](/packages/AneuFinder)
----------

Changes in version 1.12.1

BUG FIXES

- Moved ggplot2 from Imports to Depends. This was necessary because
  cowplot doesn't call ggplot2 explicitly anymore.

[AnnotationForge](/packages/AnnotationForge)
---------------

Changes in version 1.28.0

NEW FEATURES

- (version 1.27.1) `makeOrgPackage()` supports GO ontologies.

[AnnotationHub](/packages/AnnotationHub)
-------------

Changes in version 2.17.0

NEW FEATURES

- (2.17.6) remove debugging message of loading resource (`AH: 1`)

- (2.17.5) system environment variable to control localHub option for
  creating hub based only on previously downloaded resources

- (2.17.9) Allow force redownload of Hub sqlite file with refreshHub

- (2.17.12) Only display download message when something to download.

- (2.17.13) The output list of files is the AH/EH id not
  AHid:resourceid

BUG FIXES

- (2.17.4) Fix localHub when no internet connection.  The internal use
  of isDevel was preventing Hub creation when no internet connection.
  Fixed by checking connection. This code pretained to orgDb filters

- (2.17.8) On chance of very first download of hub failure, next call
  to construtor will redownload

- (2.17.10) Fix ability to use hubs when offline

- (2.17.11) Add BiocVersion to Imports. Fixes bug with R CMD check when
  testing if library can be loaded off search path. BiocManager doesn't
  Import BiocVersions and this is needed to get the correct BiocManager
  version of the snapshot date.

[AnnotationHubData](/packages/AnnotationHubData)
-----------------

Changes in version 1.15.0

MODIFICATIONS

- 1.15.13 Added "BLOB" as a valid source type

- 1.15.7 Added "MTX" as a valid source type

- 1.15.6 Expanded documentation to clarify that data can be hosted
  publically not strictly Bioconductor AWS

- 1.15.4 Added "XLS/XLSX" as valid source type

INTERNAL BUG CORRECTION

- 1.15.11 updated GencodeGFF recipes for potential future use (still
  would revisit this with another update to do like ensembl on the fly)

- 1.15.5 remove validity check that is wrong/outdated

- 1.15.1 needToRerunNonStandardOrgDb added as helper function for when
  generating non standard org dbs. 1.15.3 added try catch in case aws
  buckets unreachable.

[APAlyzer](/packages/APAlyzer)
--------

Changes in version 0.99.2 (2019-08-18)

- Revise formats to consistent with Bioconductor coding styles.

- Added unit tests

Changes in version 0.99.1 (2019-07-14)

- Submitted to Bioconductor

[apeglm](/packages/apeglm)
------

Changes in version 1.8.0

- Bringing in Josh Zitovsky's work on fast beta binomial GLM fitting in
  C++.

[aroma.light](/packages/aroma.light)
-----------

Changes in version 3.15.1 (2019-08-28)

BUG FIXES

- wpca() for matrices had a 'length > 1 in coercion to logical' bug.

Changes in version 3.15.0 (2019-05-02)

- The version number was bumped for the Bioconductor devel version,
  which is now BioC 3.10 for R (>= 3.6.0).

[ASpediaFI](/packages/ASpediaFI)
---------

Changes in version 0.99.10 (2019-10-25)

- Fixed issue from quantifying RI and ALSS

Changes in version 0.99.9 (2019-10-22)

- Fixed issue from quantifying single-end reads

Changes in version 0.99.8 (2019-10-17)

- Made additional changes to conform to Bioconductor guidelines

Changes in version 0.99.7 (2019-10-17)

- Made changes to conform to Bioconductor guidelines

Changes in version 0.99.6 (2019-10-16)

- Replaced "1:" with "seq_len".

Changes in version 0.99.5 (2019-10-15)

- Implemented S4 class instead of reference class.

Changes in version 0.99.4 (2019-10-01)

- Fixed vignette issue.

Changes in version 0.99.3 (2019-10-01)

- Fixed package loading issue.

Changes in version 0.99.2 (2019-10-01)

- Removed .Rproj file

Changes in version 0.99.1 (2019-10-01)

- Added NEWS file.

Changes in version 0.99.0 (2019-10-01)

- Submitted to Bioconductor.

[AssessORF](/packages/AssessORF)
---------

Changes in version 1.3.1 (2019-10-22)

- Data for six new strains, S. pyogenes AP1, E. coli BW25113, S. aureus
  HG001, L. lactis MG1363, B. subtilis NCIB 3610, and L. monocytogenes
  10403S, have been added to AssessORFData

- CITATION file added as corresponding paper has now been published

- Updated MapAssessmentData such that verbose output is cleaner and
  more informative # AssessORF 1.1

[ATACseqQC](/packages/ATACseqQC)
---------

Changes in version 1.9.9

- export prepareBindingSitesList function.  - Add rownames for
footprintsScanner counts data.

Changes in version 1.9.8

- Add error message for vPlot when no paired reads in bam file.

Changes in version 1.9.7

- Fix the bug that gscore changed the output for
splitGAlignmentsByCut.

Changes in version 1.9.6

- Try to decrease the memory cost for splitGAlignmentsByCut.

Changes in version 1.9.5

- Try to decrease the memory cost for splitGAlignmentsByCut.

Changes in version 1.9.4

- Add the error handle if not enough mononucleosome reads for
splitGAlignmentsByCut.

Changes in version 1.9.3

- Try to decrease the memory cost for splitGAlignmentsByCut.

Changes in version 1.9.2

- Fix the bug if the bam file containsupplementary alignments.

Changes in version 1.9.1

- Fix the bug if the bam file contain mix of single ends and paired
ends.

[AUCell](/packages/AUCell)
------

Changes in version 1.7

- Added function cbind()

[Autotuner](/packages/Autotuner)
---------

Changes in version 1.0.1 (2019-07-06)

- Added DOI to cite code through zenodo.

Changes in version 1.0.0 (2019-07-06)

- Introduced first official release of AutoTuner.

Changes in version 0.99.9 (2019-08-06)

- Added changes to satisfy second round of bioc review o Fixed source
  path for mzDb data object o Removed redundant table of contents o
  Removed non used redundant data object eicParamEsts.rds

Changes in version 0.99.8 (2019-08-02)

- Added changes to fix warnings o Replaced class == to is(). o Added
  corrections to dependencies regarding xcms.

Changes in version 0.99.7 (2019-08-02)

- Introduced changes as part of bioconductor review process. o Added
  Author info to vignette. o Added installation instructions for
  package. o Removed mention of mmetspData in favor of test data
  package mtbls2. o Added table of contents to vignette. o Removed png
  within vignette folder that was not being used.  o Changed dimensions
  of image returned from TIC plot in vignette. o Removed "-" from
  documentation. o Fixed source paths on files using internal data for
  examples o Removed peak_table and peak_difference from data/ dir o
  Added documentation for eicParamsEsts object used in testing o
  swapped sapplys for vapplys o switched 1:.. to seq_len or seq_along o
  Added and implemented accessor and setter functions for Autotuner
  slots o Added unit tests for accessor functions.

[AWFisher](/packages/AWFisher)
--------

Changes in version 0.99.2 (2019-06-14)

- Submitted to Bioconductor

[bamsignals](/packages/bamsignals)
----------

Changes in version 1.17.3

- DESCRIPTION updated

- Fixed WARNING: bamsignals.cpp:521:9: warning: ignoring return value
  of function declared with 'warn_unused_result' attribute
  &#91;-Wunused-result&#93; bamsignals.cpp:516:5: warning: ignoring return
  value of function declared with 'warn_unused_result' attribute
  &#91;-Wunused-result&#93; bamsignals.cpp:530:5: warning: ignoring return
  value of function declared with 'warn_unused_result' attribute
  &#91;-Wunused-result&#93;

Changes in version 1.17.2

- Maintainer E-mail adress updated

[BANDITS](/packages/BANDITS)
-------

Changes in version 1.2.0

- Added mean and standard deviation estimates of the precision
  parameter

- Allow estimation of parameters (without testing) when 1 group only is
  provided

- Added 1 section of the vignette for inference with 1 group only

- Added reference to BANDITS manuscript

[BASiCS](/packages/BASiCS)
------

Changes in version 1.7.20 (2019-10-17)

- Bumps version number to trigger new build.

Changes in version 1.7.19 (2019-10-08)

- BASiCS_TestDE now checks to ensure that both input chains have been
  run with Regression = FALSE or both with Regression = TRUE.

- Remove duplicated text from `BASiCS_Chain`'s show method.

Changes in version 1.7.18 (2019-10-06)

- Refactor HVG/LVG plots code to use `ggplot2`.  When calling
  `BASiCS_DetectHVG` or `BASiCS_DetectLVG`, the plots are now stored in
  the returned list in the named element `Plots`. The vignette has been
  updated to show the usage of this functionality.

- In all unit tests: use `expect_equal(foo, bar)` instead of
  `expect_true(all.equal(foo, bar))`. This gives better printing (ie,
  if it's not true it tells by what margin)

- Rename some functions for consistency in capitalisation. For example,
  `BASiCS_showFit` has been deprecated and renamed to `BASiCS_ShowFit`.

- Add a `Smooth` argument to `BASiCS_DiagPlot`.

- BASiCS_TestDE is now explicit about requiring the same number of
  samples in both chains. Previously it would fail due to arrays of
  non-conformable dimensions.

- Add some internal utility functions.

- Remove exportPattern("^[^\\Hidden]") meaning functions must be
  explicitly exported using roxygen tags.

Changes in version 1.7.17 (2019-10-04)

- Minor change in `newBASiCS_Data` to avoid missing `colnames` when
  adding `colData`

- New unit test to verify that `colnames` are not lost

Changes in version 1.7.16 (2019-10-04)

- Minor changes to NAMESPACE to pass R CMD Check

Changes in version 1.7.15 (2019-09-30)

- Swtich from `matrixStats` dependency to `Matrix` as it supports more
  general input classes (including sparse matrices)

- `newBASiCS_Data` now requires input counts to be a `matrix`.

Changes in version 1.7.14 (2019-09-27)

- Changed internal functions from `matrixStats` to `Matrix` to support
  more classes of matrix (eg, `dgCMatrix`). BASiCS_MCMC now supports
  DelayedArray, dgEMatrix, dgCMatrix objects and likely more.

Changes in version 1.7.13 (2019-09-25)

- Updated unit tests to account for new default choice for `min.mean`
  parameter in `scran::computeSumFactors`

- New unit test to check for changes in `scran::computeSumFactors`

Changes in version 1.7.12 (2019-09-15)

- changes maintainer's email address

Changes in version 1.7.11 (2019-09-15)

- `Uncertainty` parameter added to `BASiCS_showFit`. This enables
  optional inclusion of uncertainty measure around the regression
  trend.

- Re-ordering of parameter in c++ MCMC samplers (does not affect
  output)

- Avoid unnecessary transformations between NumericMatrix/Vector and
  arma::

- Removes `isSpike` usage in `BASiCS_Sim` documentation

- Removes `isSpike` call from the vignette

- Adds missing parameters to the documentation of
  `BASiCS_CorrectOffset`

Changes in version 1.7.10

- `metadata(Data)$SpikeInput` is now required to be a `data.frame`.
  This allows us to verify the correct order in spike-in inputs when
  the user manually modifies an existing `SingleCellExperiment` object.

- Additional changes in `newBASiCS_Data`,
  `HiddenBASiCS_MCMC_InputCheck`, `HiddenChecksBASiCS_Data` to avoid
  errors.

- `HiddenBASiCS_MCMC_Start`, `HiddenBASiCS_MCMC_GlobalParams`,
  `BASiCS_DenoisedCounts` and `BASiCS_DenoisedRates` modified to
  replace `isSpike` by `altExp`

- Some unit tests adapted accordingly

Changes in version 1.7.9 (2019-09-05)

- `newBASiCS_Data`, `HiddenBASiCS_MCMC_InputCheck`,
  `HiddenChecksBASiCS_Data` and `BASiCS_MCMC` modified to replace
  `isSpike` by `altExp`

- WIP - unit test do not pass

Changes in version 1.7.8 (2019-08-20)

- `min.mean` parameter exposed in `BASiCS_CorrectOffset` and
  `BASiCS_TestDE`

Changes in version 1.7.7 (2019-08-20)

- Preliminary version of `BASiCS_CorrectOffset` has been added. This
  includes a trimmed option for the offset calculation that excludes
  lowly expressed genes. This is similar to what is implemented in
  `scran:::.rescale_clusters

Changes in version 1.7.6 (2019-08-19)

- Offset correction within `BASiCS_TestDE` modified to use `rowMedians`
  instead of `rowSums2`. This makes it more robust to outlier genes.

- Unit tests updated accordingly.

Changes in version 1.7.5 (2019-07-30)

- Returns to original updates for X (fixed after burn-in)

Changes in version 1.7.4 (2019-07-30)

- Show method for regression objects now correctly shows number of
  cells.

- `HiddenBASiCS_MCMCcppReg` and `HiddenBASiCS_MCMCcppRegNoSpikes` (C++
  code) modified to update design matrix throughout MCMC (GRBF
  locations remain fixed after the burn-in period is over)

- Unit tests updated accordingly.

Changes in version 1.7.3 (2019-07-29)

- Minor typo fixed in `BASiCS_MCMC`.

Changes in version 1.7.2 (2019-07-28)

- Specific minimum tolerance thresholds (e.g. 1e-3 for mu updates)
  replaced by global parameters (e.g. `mintol_mu`). Optional
  parameters; internal use only.  Incorporated within
  `HiddenBASiCS_MCMC_ExtraArgs`.

- General clean up to remove code redundancy in `BASiCS_MCMC`

- `HiddenBASiCS_MCMC_GlobalParams` created to facilitate clean-up above

Changes in version 1.7.1 (2019-07-27)

- `is_true` deprecated in `testthat`. Unit tests updated to use
  `expect_true`

[batchelor](/packages/batchelor)
---------

Changes in version 1.2.0

- Deprecated rotate.all= in favour of get.all.genes= in
  multiBatchPCA().

- Switched BSPARAM= to use IrlbaParam(deferred=TRUE) by default in
  fastMNN(), so that the default behaviour is actually fast.

- Deprecated auto.order= in favor of merge.order= and auto.merge= in
  fastMNN() and mnnCorrect(). Automatic merging now detects potential
  tree-based merges. Merge trees can also be specified as input.

- Added the correctExperiments() function to cbind the original assays
  alongside the merged values.

- Added the subset.row= argument to cosineNorm() for in-place
  subsetting.

- Added batch= and preserve.single= arguments to multiBatchNorm().
  Standardized behavior of subset.row= by adding a normalize.all=
  argument.

- Added the regressBatches() function for correction via standard
  linear regression.

- Added the prop.k= argument in all MNN-related functions, to allow the
  value of k to adapt asymmetrically to the size of each batch.

[bcSeq](/packages/bcSeq)
-----

Changes in version 1.7.1

- using array back structure for the library to speed up the mapping
  process

- depending on the configuration of the mechine, the alignment may
  speed up

- to a factor as much as 3.

[BEclear](/packages/BEclear)
-------

Changes in version 2.1.4 (2019-05-17)

- minor changes

Changes in version 2.1.3 (2019-05-13)

- Performance improvments o regarding function calcSummary o and
  calcScore

Changes in version 2.1.2 (2019-05-10)

- New feature o y-axis in boxplot function can now be logged

Changes in version 2.1.1 (2019-05-09)

- New feature o Implementation of outlier detection

[BiocNeighbors](/packages/BiocNeighbors)
-------------

Changes in version 1.4.0

- Allow memory-efficient retrieval of the distance to the furthest
  neighbors.

- Added a warn.ties= argument to turn off tie-related warnings in the
  KMKNN and VP tree algorithms.

- Return neighbor counts in rangeFind*() and rangeQuery*() functions
  when get.index=FALSE and get.distance=FALSE.

[BioCor](/packages/BioCor)
------

Changes in version 1.10.1

- Added a pkgdown website

Changes in version 1.8.1

- Exported inverseList function

- Improved spelling

- Corrected a bug related to logical coercion of length greater than 1.

[BiocParallel](/packages/BiocParallel)
------------

Changes in version 1.20

BUG FIXES

- (v 1.19.2) Improve efficiency of MulticoreParam() when state does not
  persist across calls to bplapply().

[BiocSingular](/packages/BiocSingular)
------------

Changes in version 1.2.0

- Added the ResidualMatrix class for computing PCA on residuals
  efficiently.

- Fixed runIrlba() to avoid errors at the limit of available PCs.

- Added the FastAutoParam class to automatically choose a fast SVD
  depending on the matrix representation.

- Added the bsparam() function to quickly set or get a global default
  algorithm choice.

[biomaRt](/packages/biomaRt)
-------

Changes in version 2.42.0

NEW FEATURES

- The results of queries will now be cached, and if repeated queries
  are detected the results are loaded from disk.

MINOR CHANGES

- Ensembl users will be redirected to their closest mirror unless the
  host argument is explicitly provided.  In this case the defined value
  will be enforced.

- Unused argument 'ssl.verifypeer' removed from listMarts() and
  useMarts().

- RCurl removed from package dependecies.

BUG FIXES

- Improvements made to selecting the correct port when using http vs
  https

- Results that contain unescaped new line characters are now returned
  successfully.

[BioMM](/packages/BioMM)
-----

Changes in version 1.1.10

- updated tutorial: added circular plot

Changes in version 1.1.9

- updated tutorial: updated description for two parallel computing;
  cirPlot4pathway() example added; seeds added

- updated plotRankedFeature().

Changes in version 1.1.8

- updated tutorial: fixed typo 'param1' to 'param2'

- added a new function cirPlot4pathway()

Changes in version 1.1.7

- updated tutorial

Changes in version 1.1.6

- improved plotVarExplained() and plotRankedFeature()

- renamed R functions for getDataAfterFS, BioMMreconData,
  BioMMstage1pca

- updated installation approaches (including R 3.5 from Github)

- updated tutorial to adapt the usage of the parallel package installed
  from Github

Changes in version 1.1.5

- updated BioMM(); 'dataMode' added.

- fixed roc() in getMetrics()

Changes in version 1.1.4

- removed the examples with omcis2chrlist()

- updated the omics2pathlist()

Changes in version 1.1.3

- added required package pROC for getMetrics()

- updated NAMESPACE and Rd file.

Changes in version 1.1.2

- added the library pROC in the tutorial.

Changes in version 1.1.1

- updated functions getMetrics(), plotVarExplained() and
  plotRankedFeature(); to focus on pathway based result report.

- updated BioMMtutorial.Rmd

Changes in version 1.1.0

- removed gene and chromosome based stratification methods.

- updated feature selection method based on filtering getDataAfterFS().

- updated BioMM() function to focus on pathway based machine learning.

- updated BioMMtutorial.Rmd

[BioQC](/packages/BioQC)
-----

Changes in version 1.13.1

- Functions to manipulate GmtList objects are considerably expanded

- All documents and namespaces are now managed by roxygen2

- readGmt by default read unique genes from GMT files

[biosigner](/packages/biosigner)
---------

Changes in version 1.13.20

BUG FIXED

- plot.biosignMultiDataSet correction to include all plots in the .pdf
  file

Changes in version 1.13.18

MINOR MODIFICATION

- Minor correction in the plot.biosignMultiDataSet documentation

Changes in version 1.13.16

MINOR MODIFICATION

- Minor correction in the biosignMultiDataSet class documentation

Changes in version 1.13.14

MINOR MODIFICATION

- Minor correction in the biosignMultiDataSet class documentation

Changes in version 1.13.12

INTERNAL MODIFICATION

- Minor internal modification

Changes in version 1.13.10

INTERNAL MODIFICATION

- Minor internal modification

Changes in version 1.13.8

NEW FEATURE

- final call to warnings() omitted

Changes in version 1.13.6

NEW FEATURES

- seedI argument now available (set to 123 by default)

- biosign can now be applied to MultiDataSet objects (getMset method to
  get the updated MultiDataSet back)

Changes in version 1.13.4

MINOR MODIFICATION

- minor internal modification

Changes in version 1.13.2

NEW FEATURE

- info.txtC and fig.pdfC argument values NULL and NA replaced by 'none'
  and 'interactive', respectively

[blacksheepr](/packages/blacksheepr)
-----------

Changes in version 0.99.13 (2019-10-11)

- Final Clean up of code

- Added normalization helper function

Changes in version 0.99.9 (2019-09-30)

- Major rewrite to functionalize outlier_analysis and reduce redundant
  code

- Fixed Reviewer suggestions o changed assigner symbol o renamed
  vignette o changed getwd() -> tempdir()

Changes in version 0.99.7 (2019-09-10)

- Submitted to Bioconductor

- Fixed Reviewer suggestions o Got rid of lazydata loading in
  DESCRIPTION o properly formetted NEWS object

[breakpointR](/packages/breakpointR)
-----------

Changes in version 1.2.1

- Added new genotyping method based on binomial probabilities see
  parameter genoT = 'fisher' or 'binom'; default: 'fisher'

[brendaDb](/packages/brendaDb)
--------

Changes in version 0.99.21 (2019-10-17)

- Fix: BiocycPathwayGenes now deals with finding multiple Ensembl IDs

Changes in version 0.99.20 (2019-08-16)

- Feature: function to extract field information from brenda.entries
  objects

- Feature: DownloadBrenda now utilizes BiocFileCache

- Performance: ReadBrenda is now 50% faster

Changes in version 0.99.10 (2019-08-12)

- Doc: updated the vignette, readme and package help page

- Doc: added documentation for data file acronyms.RData

Changes in version 0.99.0 (2019-07-22)

- Submitted to Bioconductor

[BUSpaRse](/packages/BUSpaRse)
--------

Changes in version 0.99.25 (2019-09-11)

- Added message to indicate when get_velocity_files is extracting
exon-exon junctions.  - Restored "separate" to be default
isoform_action in get_velocity_files.

Changes in version 0.99.24 (2019-09-06)

- When sorting tr2g from file, now the file must be formatted in ways
required by bustools.

Changes in version 0.99.23 (2019-09-06)

- Previous two version bumps did not accompany changes; those were
used to trigger rebuilds on Bioconductor.  - Added the functionality
to use L-1 (L is read length) bases around exon-exon junction to
better distinguish between spliced and unspliced transcripts for RNA
velocity.  - Fixed serious problem with get_velocity_file that
counted reads mapping to exons of length between L-1 and 2(L-1) as
from unspliced transcript. This was done by an reimplementation of
the method to get flanked intronic ranges.  - Changed default
isoform_action to "collapse".  - Make sure that all transcripts in
tr2g.tsv from get_velocity_files are in the transcriptome.

Changes in version 0.99.20 (2019-08-26)

- Addressed Bioconductor review.

Changes in version 0.99.19 (2019-07-23)

- Finished get_velocity_files to generate files required for kallisto
| bustools RNA velocity.

Changes in version 0.99.0 (2019-06-21)

- Submitted package to Bioconductor for review

[CAGEfightR](/packages/CAGEfightR)
----------

Changes in version 1.5.3

- Added citation to BMC Bioinformatics paper.

Changes in version 1.5

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

[CAGEr](/packages/CAGEr)
-----

Changes in version 1.28.0

BACKWARDS-INCOMPATIBLE CHANGES

- The CTSS objects are now based on the UnstitchedGPos class instead of
  the parent GPos class.  Existing CAGEexp objects might give errors
  after package upgrade.

NEW FEATURES

- Depends on MultiAssayExperiemnt instead of previously Importing it.
  No need to load it by hand (nor SummaryExperiment) anymore.

BUG FIXES

- Correct strand in .remove.added.G private function (PR#26).

[calm](/packages/calm)
----

Changes in version 0.99.0 (2019-07-18)

- Submitted to Bioconductor

[CAMERA](/packages/CAMERA)
------

Changes in version 1.41.1

NEW FEATURES

- allow findAdducts alternative intval (thanks to Sunil Dhakad)

[canceR](/packages/canceR)
------

Changes in version 1.19.05

- move vignette to inst/doc

Changes in version 1.19.04

- omit GSEAlm dependency

- include used functions from gsealm to canceR package in gsealm.R
  file.

Changes in version 1.19.03

- compress RData files to RDS and move them from /data to
  /extdata/rdata

- update running examples

Changes in version 1.19.02

- remove /extdata/gct_cls

- do not import grDevices::quartz

- add empty line in test functions

Changes in version 1.19.01

- remove .txt file from /data.

- import grDevices package

- Adjust the vignetteEngine package

[Cardinal](/packages/Cardinal)
--------

Changes in version 2.3.18 (2019-10-27)

SIGNIFICANT USER-VISIBLE CHANGES

- Processing in 'crossValidate()' now allows processing unprocessed
  data by performing peak-picking on the mean spectra of the training
  sets

BUG FIXES

- Fixed some errors in user messages during peak processing

Changes in version 2.3.17 (2019-10-25)

SIGNIFICANT USER-VISIBLE CHANGES

- Default for 'peakBin()' argument 'type' is now "area"

- In 'peakBin()', peak boundaries should be calculated more accurately
  now, and general speed improvements

- In 'peakAlign()', peak centers are now calculated as weighted average
  mass rather than the highest point

Changes in version 2.3.16 (2019-10-14)

NEW FEATURES

- New class 'ImagingSummary' with sub-classes including
  'SparseImagingSummary' and 'MSImagingSummary' with appropriate
  'plot()' and 'image()' methods

SIGNIFICANT USER-VISIBLE CHANGES

- The 'summarize()' method for 'SparseImagingExperiment' now returns a
  'SparseImagingSummary', to more closely reflect "tidy" data
  principles by returning an object of a similar class; the previous
  behavior can be reproduced by specifying '.as="DataFrame"'

Changes in version 2.3.15 (2019-10-13)

NEW FEATURES

- For methods requiring 'resolution' or 'tolerance', the default
  arguments have been updated to automatically guess based on the data

Changes in version 2.3.14

NEW FEATURES

- Add spectraData() as an alias for 'imageData()' for
  'MSImagingExperiment' sub-classes

- Formalize 'mzData()' and 'intensityData()' getters and setters for
  'MSProcessedImagingExperiment'

- Add 'peaks()' and 'peakData()' methods for extracting peak matrices
  and/or peak information

- Add 'isCentroided()' method for guessing whether spectra are
  centroided (without using the @centroided slot)

SIGNIFICANT USER-VISIBLE CHANGES

- Allow 'NA' for @centroided slot for 'MSImagingExperiment'

- 'mzBin()' method now sets centroided = NA

- Update 'mzFilter()' with parameter defaults so that 'thresh.max = NA'
  and new arg 'rm.zero = TRUE'

- Log more pre-processing information (e.g., method name)

Changes in version 2.3.13

BUG FIXES

- Try using 'parent.frame(1)' instead of 'parent.frame(2)' to fix NSE
  methods when used in LHS of a maggritr pipe

- Fix weird 'iData()<-' missing argument 'i' bug

Changes in version 2.3.12

SIGNIFICANT USER-VISIBLE CHANGES

- Changed default 'peakPick()' method to 'mad'

- In 'peakPick()' method 'mad', change the default number of blocks to
  1 (no adaptive smoothing)

- In 'peakPick()' method 'mad', update w/ new arguments w/ new defaults
  'fun=median' and 'tform=diff'

BUG FIXES

- In 'peakPick()' methods 'simple' and 'adaptive', warn if kurtosis
  cannot be estimated and try to recover

- In 'normalize()' method 'reference', provide a warning if the
  reference value is 0 for a pixel

Changes in version 2.3.11

SIGNIFICANT USER-VISIBLE CHANGES

- Improved speed in 'spatialFastmap()'

- Improved speed in 'spatialShrunkenCentroids()'

- New dissimilarity metrics for 'spatialFastmap()' including a new
  default metric='average'

BUG FIXES

- Fix error in 'print()' for facet plots where lims=NULL

Changes in version 2.3.10

SIGNIFICANT USER-VISIBLE CHANGES

- Improved speed in 'spatialDGMM()', by moving spatial filtering of
  probabilities to C code, up to 10x faster

- Linesearch in 'spatialDGMM()' now uses 'optimize()' rather than
  'optim()' -- results may differ slightly

Changes in version 2.3.9

NEW FEATURES

- Added 'Cardinal.history()' and 'Cardinal.version()'

SIGNIFICANT USER-VISIBLE CHANGES

- Removed under-used setter generic function definitions

BUG FIXES

- Cleaned up generics to reflect ProtGenerics >= 1.17.2

Changes in version 2.3.8

NEW FEATURES

- Added boxplot, histogram, and bar chart functionality to the 'plot()'
  method for 'XDataFrame'

- Added 'plot()' plotting for 'AnnotatedImageList'

- Added 'plot()' methods for 'SpatialDGMM', 'MeansTest', and
  'SegmentationTest' result classes

- Added 'image()' method for 'MeansTest' result class

SIGNIFICANT USER-VISIBLE CHANGES

- Updated 'plot()' and 'image()' methods for 'SpatialDGMM',
  'MeansTest', and 'SegmentationTest' result classes

BUG FIXES

- Various bug fixes in object printing and plot auto-layout

Changes in version 2.3.7

NEW FEATURES

- Added 'image()' plotting for 'AnnotatedImageList'

BUG FIXES

- Plotting with 'add=TRUE' now respects 'par('usr')' coordinates

Changes in version 2.3.6

NEW FEATURES

- Add 'AnnotatedImageList' class for list of 'AnnotatedImage' objects

- Add 'AnnotatedImagingExperiment' class for containing data for an
  optical imaging experiment (e.g., a microscopy experiments)

- Add 'image()' plotting for 'AnnotatedImagingExperiment'

SIGNIFICANT USER-VISIBLE CHANGES

- Redefine '@featureData' slot of a 'SparseImagingExperiment' to be a
  'DataFrame' rather than requiring an 'XDataFrame'

BUG FIXES

- Respect 'layout' and 'byrow' passed through ... args to 'print()'
  method on facet plot objects

Changes in version 2.3.5

SIGNIFICANT USER-VISIBLE CHANGES

- Moved some S4 method definitions from 'ImagingExperiment' to
  'SparseImagingExperiment' so that the former can be more flexible for
  a wider variety of imaging modalities

BUG FIXES

- Pass more ... args through to 'par()' in plotting functions

Changes in version 2.3.4

NEW FEATURES

- Added 'AnnotatedImage' class for optical images

Changes in version 2.3.3

SIGNIFICANT USER-VISIBLE CHANGES

- Improved facet plotting when 'add=TRUE'

BUG FIXES

- Better 'cex.axis' defaults and user setting for colorkeys

Changes in version 2.3.2

SIGNIFICANT USER-VISIBLE CHANGES

- Improved 'writeMSIData()' for 3D and non-gridded data

Changes in version 2.3.1

NEW FEATURES

- Output directly to imzML while processing with 'process()'

SIGNIFICANT USER-VISIBLE CHANGES

- Improved auto-layout for visualization with multiple runs

- Added 'parse.only' option to 'readImzML()' for parsing only

Changes in version 2.2.4

BUG FIXES

- Fix large external array offsets in 'writeImzML'

Changes in version 2.2.3

BUG FIXES

- Cleaned up some 'writeImzML' mapping validity issues

Changes in version 2.2.2

BUG FIXES

- Removed curly braces around UUID when writing imzML

Changes in version 2.2.1

BUG FIXES

- Fixed bug in plotting results where 'column' argument would get
  matched before the 'col' argument

[cbaf](/packages/cbaf)
----

Changes in version 1.8.0 (2019-10-28)

New Features

- Optimization has led obtainOneStudy() and obtainMultipleStudies()
  functions to work faster.

- If an entered cancer has corrupted data or lacks the requested data
  type, obtainMultipleStudies() doesn't return error. That study is
  automatically omitted from results and its name is printed on
  console.

- AvailableData() function now works more accurately but unfortunately,
  it is generally slower: Due to inconsistancy in the terms that cgdsr
  uses, AvailableData() has to check the availability of the data at
  two different levels.

[ccfindR](/packages/ccfindR)
-------

Changes in version 1.5.0

Changes

- Revised default meta_gene.cv(..., cv.max=Inf)

- Added assignCelltype(...)

Changes in version 1.4.1

Changes

- Added URL to published article.

- Revised optimal_rank(...) such that Bayes factor criterion is used.

- Fixed filter_genes(...) for cases with non-expressed genes

[celda](/packages/celda)
-----

Changes in version 1.1.6 (2019-07-16)

- Add multiclass decision tree

Changes in version 1.1.4 (2019-05-28)

- Add Alternate headings support for plotDimReduceFeature

Changes in version 1.1.3 (2019-05-14)

- Add multiclass decision tree (MCDT) cell cluster annotation

Changes in version 1.1.2 (2019-05-14)

- Fix a bug in celdaHeatmap

Changes in version 1.1.1 (2019-05-09)

- Default seed setting to maintain reproducibility

[CellBench](/packages/CellBench)
---------

Changes in version 1.1.3

Bug Fixes

- Data loading functions now appear in package index and
documentation

Changes in version 1.1.2

Bug Fixes

- Updated make_combinations to work with tidyr 1.0.0

Modifications

- Updated the WritingWrappers vignette.  - Added a case study
precompiled vignette.

Changes in version 1.1.1

New Features

- Added any_task_errors() function to check if any tasks failed in
benchmark tibble.

[CellMixS](/packages/CellMixS)
--------

Changes in version 1.1.3

- Add entropy and weighted isi function

- Add evalIntegration function

- Replace `metric_prefix` from plotting functions by `metric`
  parameter.

[ChAMP](/packages/ChAMP)
-----

Changes in version 2.15.1

- Added scree plot into champ.SVD() function.

[Chicago](/packages/Chicago)
-------

Version: 1.13
Category: Defaults are now updated based on .npb file header rather
        than just checking them for consistency with it
Text:

Version: 1.13
Category: When using Chicago with four-cutter enzymes, make sure you
        reduce minFragLen parameter. However, this can now only be done
        when generating the design files, and not in the R package
        itself
Text:

[chipenrich](/packages/chipenrich)
----------

Changes in version 2.10.0

NEW FEATURES

- A new test, proxReg(), tests for genomic region binding proximity to
  either gene transcription start sites or enhancer regions within gene
  sets. Used as an addendum to any gene set enrichment test, not
  exclusive to those in this package.

IMPROVEMENTS

- Poly-Enrich now uses the likelihood ratio test instead of the Wald
  test, as LRT is more robust when using a negative binomial GLM.

BUG FIXES

- Poly-Enrich Approximate method that uses the score test now uses the
  correct formula.

[ChIPpeakAnno](/packages/ChIPpeakAnno)
------------

Changes in version 3.19.5

- fix the issue that findOverlapsOfPeaks will connect the peaks in same
  peak list.

Changes in version 3.19.4

- add last choice for xget.

Changes in version 3.19.3

- fix the colnames of addMetadata.

Changes in version 3.19.2

- fix the issue that seqlevelsStyle(peak) == seqlevelsStyle(annotation)
  are not all TRUE

Changes in version 3.19.1

- remove RangedData

[ChIPseeker](/packages/ChIPseeker)
----------

Changes in version 1.21.1

- new implementation of upsetplot (2019-08-29, Thu) - use ggupset,
ggimage and ggplotify - subset method for csAnno object (2019-08-27,
Tue)

[chromstaR](/packages/chromstaR)
---------

Changes in version 1.11.1

BUGFIXES

- Bugfix for error when exporting empty peaks in univariateHMM.

[circRNAprofiler](/packages/circRNAprofiler)
---------------

Changes in version 0.1.1 (2019-05-30)

- First releaase

[cleanUpdTSeq](/packages/cleanUpdTSeq)
------------

Changes in version 1.23.1

- update citation.

[clusterExperiment](/packages/clusterExperiment)
-----------------

Changes in version 2.5.7

Bugs:

- Fix logical precedence error in C++ code for subsample loop

Changes in version 2.5.3

Changes

- Made "kmeans" the default in subsampling if data type is "X"

- Set checkDiss=FALSE by default for most all functions (exception is
  when user defines distance function)

- Add warnings when forced to calculate a nxn distance matrix

- Improve generic classification to centroids (used by "pam") so not
  calculate unnecessary distances.

- Added mbkmeans as built-in cluster function

Bugs

- Remove bug in built-in cluster function "kmeans" so as to not make
  unnecessary nxn distance matrix.

[CNVRanger](/packages/CNVRanger)
---------

Changes in version 1.2.0

- New function `plotRecurrentRegions` to visualize the landscape of
  recurrent CNV regions

- New function `plotEQTL` to explore differential expression of genes
  in the neighborhood of a CNV region

- Reworked vignette - dedicated section on applicability and scope, -
  overview of key functions, - extended input data format description,
  - visualizations and dedicated Gviz plots to illustrate key concepts

[cola](/packages/cola)
----

Changes in version 1.1.2

- improve documentations

Changes in version 1.1.1

- add `GO_enrichment()` and `map_to_entrez_id()`.

- add `ncol()`/`nrow()`/`colnames()`/`rownames()`/`dim()` helper
  functions.

- use `eulerr::euler()` to make the Euler diagram.

- simplified the rules for deciding the best k

[ComplexHeatmap](/packages/ComplexHeatmap)
--------------

Changes in version 2.1.1

- `Heatmap()`: give error when heatmap has empty string as its name.

- `anno_mark()`: text positions are correctly calculated now with
  rotations.

- The order of legend labels are ordered by either `sort` or `levels`.

Changes in version 2.1.0

- check the length of the clustering objects and the matrix
  rows/columns

- `anno_oncoprint_barplot()`: add `ylim` argumnet

- `anno_mark()`: add `labels_rot` argument

- `draw_legend()`: legends for annotations with the same names are
  merged

- `densityHeatmap()`: `ylim` works as it is expected.

- add `cluster_row_slices` and `cluster_column_slices` to
  `draw,HeatmapList-method()`

- `densityHeatmap()`: `col` can be set as a function

- add `cluster_rows`/`cluster_columns` in `oncoPrint()`

- legend labels support symbols

- `Heatmap()`: add `jitter` argument to add tiny random shift to
  original matrix.  It is mainly to solve the problem of "Error: node
  stack overflow" when there are too many identical rows/columns for
  plotting the dendrograms.

[consensusDE](/packages/consensusDE)
-----------

Changes in version 1.3.4 (2019-09-16)

- table_means function of multi_de_pairs updated to be n=1 aware

- add warning when less than two biological replicates

Changes in version 1.3.3 (2019-06-27)

- norm_method paramter in multi_de_pairs to allow only EDASeq
  normalisation

Changes in version 1.3.2 (2019-06-21)

- make character output for annotation

- remove prior normalisation for RUVr

Changes in version 1.3.1 (2019-06-18)

- fix buildSummarized bug in detection of minimum paired numbers

- fix buildSummarized htseq with out of order rownames

- add parameter to buildSummarized for technical_reps

- technical_reps merges reads from technical replicates

- change merged column name from p_max to p_intersect

- change LogFC of merged results to mean(LogFC) of all methods

- change AveExpr of merged results to mean(AveExpr) of all methods

- addition of merged column name p_union (Union)

- addition of merged column name LogFC_sd (Standard Deviation of FC)

- updating plotting functions to add weight for LogFC_sd

- add numbers for each category to legend of plots

- disabled cooksCutoff for DESeq2 for comparability of all p-value
  rankings

- add option gtf_annotate multi_de_pairs for annotation of gene symbols
  from gtf and combine with tx_db

- update vignette

- version 1.3.0 onwards is BioC 3.9

[consensusOV](/packages/consensusOV)
-----------

Changes in version 1.8.0

- New function `get.hao.subtypes` to predict the tissue of origin of
  ovarian tumors as either fallopian tube (FT) or ovarian surface
  epithelium (OSE) based on Hao et al., Clin Cancer Res, 2017

[CRISPRseek](/packages/CRISPRseek)
----------

Changes in version 1.25.6

- added parameter ignore.strand to indidate whether gene annotation
  should be strand-specific

Changes in version 1.25.5

- efficacy is calculated only for on-target

Changes in version 1.25.4

- annotation is now strand-specific

Changes in version 1.25.1

NEW FEATURES

- added paired.orientation parameteter

[CrossICC](/packages/CrossICC)
--------

Changes in version 0.99.27 (2019-10-22)

- Use tempdir() in test to avoid using home dir of bioconductor build
  machine

Changes in version 0.99.26 (2019-09-23)

- Use tempdir() in examples and vignette to prevent from contaminating
  Bioc and users' machine

Changes in version 0.99.24 (2019-09-23)

- Update source code for robustness

Changes in version 0.99.23 (2019-09-14)

- Add support for type SummarizedExperiment in function predictor()

- Fix output directory problem in main function

- Remove some unused functions

Changes in version 0.99.22 (2019-09-06)

- Fix bugs by removing lazy load

Changes in version 0.99.21 (2019-09-06)

- Add overwrite function for main function

- Optimize packages suggested checking

Changes in version 0.99.20 (2019-07-30)

- Fix predictor() error by complete cases

Changes in version 0.99.19 (2019-07-29)

- Fix m.f.s bug when calculating centroid of new matrix

Changes in version 0.99.18 (2019-07-29)

- Add a choice for using kmeans() for super-clustering

Changes in version 0.99.17 (2019-07-29)

- Add a parameter allowing keep rows with no variance

- List for centroid2exp() should not be filtered by variance

Changes in version 0.99.16 (2019-07-29)

- CrossICC now need R version >= 3.5

- Update ssGSEA function

Changes in version 0.99.15 (2019-07-28)

- CrossICC now need R version >= 3.6

- Fix NAMESPACE bug

Changes in version 0.99.14 (2019-07-28)

- Update vignette and man pages

Changes in version 0.99.13 (2019-07-28)

- Update vignette and man pages

Changes in version 0.99.12 (2019-07-28)

- Fix ssGSEA() bug

- Remove external shiny calling function

- CrossICC() can reset working directory to previous one now

Changes in version 0.99.11 (2019-07-28)

- Add sankey plot to shiny app

- Optimize the main function

Changes in version 0.99.10 (2019-07-27)

- Update ssGSEA man page

Changes in version 0.99.9 (2019-07-27)

- Label ssGSEA examples as donttest

Changes in version 0.99.8 (2019-07-27)

- Update help page of function CrossICCInput(): example should contain
  file name pattern

Changes in version 0.99.7 (2019-07-26)

- Update man pages

Changes in version 0.99.6 (2019-07-26)

- Use all functions in reloaded MergeMaid code

- To test all examples in exported function

Changes in version 0.99.5 (2019-07-26)

- Turn off default shiny calling

Changes in version 0.99.4 (2019-07-26)

- Update some man pages

Changes in version 0.99.3 (2019-07-26)

- Update examples for most functions

- Remove random seed of ConcensusClusterPlus function

Changes in version 0.99.2 (2019-07-25)

- Fix some mistakes in Rd files

Changes in version 0.99.1 (2019-07-25)

- Main function CrossICC will not use shiny app by default

- Add sankey plot

- Remove some default dependencies (only needed when use shiny)

- Fix some check error and warnings

Changes in version 0.1.1 (2019-06-25)

- Update vignette and some man pages for functions

Changes in version 0.1.0 (2019-06-25)

- The first version of CrossICC

- Submitted to Bioconductor

[csaw](/packages/csaw)
----

Changes in version 1.20.0

- Removed deprecated functionality in normOffsets(), readParam(),
  scaledAverage().

- Added mergeResults(), overlapResults() wrapper functions to simplify
  getting region-level results.

- Added mergeWindowsList(), findOverlapsList() functions for
  consolidating windows from multiple analyses. Deprecated
  consolidateWindows().

- Added mergeResultsList(), overlapResultsList() wrapper functions to
  obtain consolidated region-level results. Deprecated
  consolidateTests(), consolidateOverlaps().

- Renamed regions= to ranges= in mergeWindows() for consistency.

- Added clusterWindowsList() to replace consolidateClusters().

- Added filterWindowsGlobal(), filterWindowsLocal(),
  filterWindowsProportion() and filterWindowsControl(). Deprecated
  filterWindows().

- Minor renaming of scaleControlInfo() arguments, added assay.data= and
  assay.back= arguments.

[cTRAP](/packages/cTRAP)
-----

Changes in version 1.4

New features

- Predict targeting drugs (predictTargetingDrug()): - Based on
expression and drug sensitivity associations derived from NCI60, CTRP
and GDSC data (see loadExpressionDrugSensitivityAssociation()) -
Compare user-provided differential expression profile with gene
expression and drug sensitivity associations to predict targeting
drugs and their targeted genes - Compounds are ranked based on their
relative targeting potential - Plot candidate targeting drugs against
ranked compound perturbations using
plotTargetingDrugsVSsimilarPerturbations(), highlighting compounds
that selectively select against cells with a similar differential
gene expression profile - Analyse drug set enrichment
(performDSEA()): - Prepare drug sets based on a table with compound
identifiers and respective 2D and 3D molecular descriptors using
prepareDrugSets() - Test drug set enrichment on results from
rankSimilarPerturbations() (when ranking against compound
perturbations) and predictTargetingDrugs() - Convert ENSEMBL
identifiers to gene symbols using convertENSEMBLtoGeneSymbols()

Major changes

- Update the tutorial and function documentation - Remove most L1000
instances, including in function names: - getL1000perturbationTypes()
-> getCMapPerturbationTypes() - getL1000conditions() ->
getCMapConditions() - downloadL1000data() -> loadCMapData() -
filterL1000metadata() -> filterCMapMetadata() -
loadL1000perturbations() -> prepareCMapPerturbations() -
compareAgainstL1000() -> rankSimilarPerturbations() -
plotL1000comparison() -> plot() - Improve loading of ENCODE samples
(loadENCODEsamples()): - Rename function from downloadENCODEsamples()
to loadENCODEsamples() - Load ENCODE samples regarding multiple cell
lines and experiment targets using loadENCODEsamples() - Improve CMap
data and metadata retrieval: - By default, do not return control
perturbation types when using getCMapPerturbationTypes() (unless if
using argument control = TRUE) - Parse CMap identifiers using
parseCMapID() - Load CMap's compound metadata using loadCMapData() -
Ask to download CMap perturbations z-scores file for differential
expression if not found (avoiding downloading a huge file without
user consent) - Improve preparation of CMap perturbations
(prepareCMapPerturbations()): - Allow to load CMap metadata directly
from files when using file paths as arguments of
prepareCMapPerturbations() - Significantly decrease memory required
to use cTRAP by loading chunks of z-scores from CMap perturbations
on-demand (a slight decrease in time performance is expected), unless
prepareCMapPerturbations() is run with argument loadZscores = TRUE -
Display summary of loaded perturbations after running
prepareCMapPerturbations() - Improve ranking of similar perturbations
(rankSimilarPerturbation()): - Redesigned output: long (instead of
wide) table - By default, calculate mean across cell lines if there
is more than one cell line available; disabled if argument
cellLineMean = FALSE - Allow to rank (or not) individual cell line
perturbations (argument rankIndividualCellLinePerturbations) when the
mean is calculated - Allow to perform multiple comparison methods if
desired (by providing a vector of supported methods via the method
argument) - Calculate the rank product's rank to assess ranks across
multiple methods - Sort results based on rank product's rank (or the
rank of the only comparison method performed, otherwise) - Include
information for calculated means across cell lines in metadata -
Include run time as an attribute - Improve metadata display for a
similarPerturbations object, obtained after running
rankSimilarPerturbations(): - Show further metadata information
(including compound data, if available) related with a given
perturbation by calling print() with a similarPerturbations object
and a specific perturbation identifier - Show a complete table with
metadata (and compound information, if available) when calling
as.table() with a similarPerturbations object - Improve plotting
(plot()): - Plot comparison results against all compared data by
calling plot() with the results obtained after running
rankSimilarPerturbations() or predictTargetingDrugs(); non-ranked
compared data can also be plotted with argument
plotNonRankedPerturbations = TRUE - Render scatter and Gene Set
Enrichment Analysis (GSEA) plots between differential expression
results and a single perturbation by calling plot() with a
perturbationChanges object (if an identifier regarding the summary of
multiple perturbations scores across cell lines is given, the plots
are coloured by cell line) - When displaying GSEA plots, plot results
for most up- and down-regulated user-provided differentially
expressed genes (by default) - Improve GSEA plot style, including rug
plot in enrichment score plot (replacing the gene hit plot)

Bug fixes and minor changes

- CMap metadata minor improvements: - Improve list returned by
getCMapConditions(), including sorting of dose and time points -
Correctly set instances of -666 in CMap metadata as missing values
and fix specific issues with metadata (such as doses displayed as 300
ng|300 ng) - In compound metadata, fix missing values showing as
literal "NA" values - CMap perturbation minor improvements: - Fix
error when subsetting a perturbationChanges object with only one row
- Improve performance when subsetting perturbationChanges objects -
Minor improvements to rankSimilarPerturbations(): - Correctly set
name of perturbations depending on their type (genes, biological
agents or compounds) - Improve performance when correlating against
multiple cell lines - Remove cellLine argument (please filter
conditions with upstream functions such as filterCMapMetadata()) -
Fix incorrect label of first column identifiers - Report run time and
settings used - Perform comparisons against perturbations
disregarding their cell lines (faster runtime) - Fix error when
trying to calculate the mean for cell lines with no intersecting
conditions available - Clearly state to the user when no intersecting
genes were found between input dataset and CMap data - Minor
improvements to plot(): - Improve rendering performance of the GSEA
plot - Fix disproportionate height between top and bottom enrichment
score panels in GSEA plots - Update demo datasets: - Update the
cmapPerturbationsCompounds and cmapPerturbationsKD datasets according
to new internal changes and fix their respective code in the
documentation - Include license and copyright text for cmapR code

[dagLogo](/packages/dagLogo)
-------

Changes in version 1.23.5

- add dontrun for prepareProteomeByFTP in case there is net connection
  issue.

Changes in version 1.23.4

- return errors when type of getSequence is not correct.

Changes in version 1.23.3

- Update documentation.

Changes in version 1.23.1

- Add function prepareProteomeByFTP

- fix the bug in addScheme

[debCAM](/packages/debCAM)
------

Changes in version 1.3.3 (2019-10-28)

- Change package name from CAMTHC to debCAM

Changes in version 1.1.1 (2018-12-28)

- Add reselectMG() to help select markers from all probes

- Add redoASest() to re-estimate A and S matrix and optionally apply
  ALS

- Add quick.select option for greedy search by sffsHull() function

- Add sample.wight option for CAM(), CAMPrep() and add SW for
  CAMPrepObj class

- Add generalNMF option which has no sum-to-one constraint

- Fix bug caused by NMF::.fcnnls() and import more robust function
  nnls::nnls()

- Fix bug in space median when dimenion is 2

- Decrease Kmeans repetition times when input data has too many data
  points

- Enhance simplex plot

[debrowser](/packages/debrowser)
---------

Changes in version 1.12.2

- Bar bax plot name fix

[deco](/packages/deco)
----

Changes in version 1.0.1

- Release 3.9 Bioconductor

- Changes in vignette.

[decompTumor2Sig](/packages/decompTumor2Sig)
---------------

Changes in version 2.1.0 (2019-08-18)

- adapted readAlexandrovSignatures to read the file format used by the
  COSMIC mutational signatures version 3 (May 2019; Single Base
  Substitution/SBS signatures only)

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

Changes in version 1.21.1

- Fix: call summary from DESeq2

- Fix: degPlot to avoid reordering `ann` vector when checking if they
  exists.

- Fix: degSignature to work with new version of melt.

- Fix: degVolcano doesn't uses ranges of 0.5 in the x axis anymore

- Fix: degMA used always raw table, now fixed to use the right one.

[DelayedArray](/packages/DelayedArray)
------------

Changes in version 0.12.0

NEW FEATURES

- Add isPristine()

- Delayed subassignment now accepts a right value with dimensions that
  are not strictly the same as the dimensions of the selection as long
  as the "effective dimensions" are the same

- Small improvement to delayed dimnames setter: atomic vectors or
  factors in the supplied 'dimnames' list are now accepted and passed
  thru as.character()

SIGNIFICANT USER-VISIBLE CHANGES

- Improve show() method for DelayedArray objects (see commit 54540856)

BUG FIXES

- Setting and getting the dimnames of a DelayedArray object or
  derivative now preserves the names on the dimnames

- Some fixes related to DelayedArray objects with list array seeds (see
  commit 6c94eac7)

[DepecheR](/packages/DepecheR)
--------

Changes in version 1.1.9 (2019-10-06)

- Two changes making package more BioConductor compliant: o
  Introduction of a "Installation" and a "Session info" slot in the
  vinjette o Formatting change of this NEWS file.

Changes in version 1.1.8 (2019-09-12)

- Bug fixes. True also for last version.

Changes in version 1.1.6 (2019-09-12)

- Adding export option for the vector created by groupProbPlot

- Correcting dColorPlot so that it can take a hex color vector as
  input, and

- keep the colors.

- Adding a test for groupProbPlot.

- Minor bug fixes for groupProbPlot.

- Small edits to the vinjette, to clarify that umap is preferred over
  tSNE and

- that dViolins currently does not fully support plotting of CyTOF data
  due to

- its sparsity.

Changes in version 1.1.5 (2019-08-07)

- 4x speedup of the new groupProbPlot.

- Proper citations of the MixOmics package, that the dSplsda function
  uses

- heavily.

- Bug fixes of the groupProbPlot.

Changes in version 1.1.4 (2019-07-30)

- A brand new function is added! Here, the user can get single-cell
  resolution

- on which group of two that a cell is more probable to belong to. See
  docs for

- further information. Highly recommended to test!

- Internal c++ test corrected by external contribution by Zimoun, which
  we are

- most grateful for.

Changes in version 1.1.3 (2019-06-28)

- dResidualPlot bug making colors sometimes represent the wrong group
  fixed.

Changes in version 1.1.2 (2019-06-28)

- Fixing the dColorPlot and the dDensityPlot functions, so that
  original names

- in character vectors and factors are correctly displayed.

[derfinder](/packages/derfinder)
---------

Changes in version 1.19.9

BUG FIXES

- Fixed an important bug in findRegions() that affected the end
positions of the regions when maxRegionGap was supplied with a value
greater than the default of 0 and the data was filtered (so position
was not all TRUE in the findRegions() call). To check this scenario
now there is a new unit test under
tests/testthat/test-maxRegionGap.R.  - The above bug went unnoticed
in getRegionCoverage() and thus for regionMatrix() for the same type
of situations (filtered data with a non-zero maxRegionGap). The
coverage values are ok, it's just the end positions of the regions
returned by findRegions() that were incorrect and that would need to
be re-computed with the fixed version.  - Changed some internal tests
to check bumphunter::loessByCluster() instead of
bumphunter::runmedByCluster() given some issues with the second one.

Changes in version 1.19.8

- Added a NEWS.md file to track changes to the package.

Changes in version 1.19.4

BUG FIXES

- railMatrix() and loadCoverage() helper functions had an issue when
the input set of regions was duplicated. This could be reproduced
with

sampleFile <- c('SRR387777' =
'http://duffel.rail.bio/recount/SRP009615/bw/SRR387777.bw') regs <-
GenomicRanges::GRanges('chrY', IRanges(start = c(1, 1), width = 10),
strand = '-') names(regs) <- c(1:2) result <-
rtracklayer::import(sampleFile, selection = regs, as = 'RleList')

This error affected recount and other reverse dependencies that use
derfinder for processing BigWig files.

Changes in version 1.19.2

BUG FIXES

- railMatrix() and loadCoverage() helper functions now attempt to
import a BigWig file 3 times before giving up. Based on
http://bioconductor.org/developers/how-to/web-query/ and
https://github.com/leekgroup/recount/commit/8da982b309e2d19638166f263057d9f85bb64e3f
which will make these functions more robust to occasional web access
issues.

[derfinderHelper](/packages/derfinderHelper)
---------------

Changes in version 1.19.3

- Added a NEWS.md file to track changes to the package.

[derfinderPlot](/packages/derfinderPlot)
-------------

Changes in version 1.19.3

- Added a NEWS.md file to track changes to the package.

[DESeq2](/packages/DESeq2)
------

Changes in version 1.26.0

- Incorporation of fast code from Constantin Ahlmann-Eltze which speeds
  up DESeq2 for large sample sizes (n > 100) by at least an order of
  magnitude. In fact the speed is now linear with number of samples
  whereas previously DESeq2 would scale quadratically. The critical
  merge commits were: c96c1c0ad43280c82403d3e6bc3501332a62e7b8
  (2019-07-16) 0a47a0c750aa5c31df759a171c737d6ed782d6c2 (2019-07-30)

- Fixed a bug where rbind() in parallel=TRUE would proliferate metadata
  items.

- Updated vignette to discuss tximeta (workflow also updated to show
  use of tximeta instead of read counting).

[DEWSeq](/packages/DEWSeq)
------

Changes in version 0.99.0 (2019-09-23)

- Submitted to Bioconductor

[DiffBind](/packages/DiffBind)
--------

Changes in version 3.13.1

- Change dba.plotPCA to use proper loading

[diffcyt](/packages/diffcyt)
-------

Changes in version 1.5.6

- Update 'diffcyt()' wrapper function to handle new
  SingleCellExperiment object structure from CATALYST.

[DiscoRhythm](/packages/DiscoRhythm)
-----------

Changes in version 1.1.1 (2019-09-27)

- Added report execution mode with zip archive of results. Results can
  now be emailed from the public server.

- Automatically fix duplicate row names.

- Allow direct upload of metadata table to web application.

- Cosinor now works with NA values (using discoODAs only).

- Bug fixes.

[DMCFB](/packages/DMCFB)
-----

Changes in version 0.99.1

New Features

- None.

CHANGES IN VERSION 0.99.1 (2019-06-02)

Added

- None.

Changed

- None.

Removed

- None.

Fixed

- None.

[DMRcate](/packages/DMRcate)
-------

Changes in version 2.0.0

- Full utility for WGBS and RRBS assays implemented using
  sequencing.annotate(): Users can either input a) A BSseq object and
  model matrix from edgeR::modelMatrixMeth, or b) Output from
  DSS::DMLtest() or DSS::DMLtest.multiFactor(),

- Major reconstruction of class types in S4: a S3 "annot" object is now
  a S4 "CpGannotated" object and a S3 "dmrcate.output" object has had
  its "input" and "pcutoff" slots removed, and "results" are now
  represented in an S4

- Improved DMR.plot() using more detailed transcript annotation from
  hg19, hg38 and mm10 GeneRegionTracks from updated DMRcatedata, as
  well as smoothed group means (group specified via "phen.col"
  argument). For bisulfite sequencing assays, the CpGs argument now
  takes a BSseq object instead of a GRanges object

- Extra DMR-level summary statistics including Fisher's multiple
  comparison test and harmonic mean of individual CpG FDRs

- Addition of the changeFDR() utility function that allows the
  re-thresholding of a "CpGannotated" object without fitting the entire
  model again

- Simplification of the rmSNPandCH() function

- Overlapping.promoters in extractRanges() are now overlapping.genes

- All vignette examples use ExperimentHub data

- Extra data object from DMRcatedata needed for rmSNPandCH(),
  extractRanges() and DMR.plot() are now in ExperimentHub

- Removal of the "p.adjust.method" argument to dmrcate() - it is
  confusing since thresholding should be performed at the
  (cpg|sequencing).annotate level

- Removal of the "samps" argument to DMR.plot() - it is redundant and
  usage can be specified by subsetting "CpGs" and "phen.col"

- Multicore processing removed since WGBS DMRs should be able to be
  produced in serial in < 1 hour

[DominoEffect](/packages/DominoEffect)
------------

Changes in version 1.5.2

UDPATE

- adding importFrom(methods, is)

- Change class() == to is()

Changes in version 1.5.1

UDPATE

- Adapt to new GenomicRanges::GPo function

Changes in version 1.5.0

UDPATE

- DominoEffect BioC 3.10 development version

[DOSE](/packages/DOSE)
----

Changes in version 3.11.2

- ignore universe and print a waring message if users passing
accidentally passing wrong input (2019-10-24, Thu) -
https://github.com/YuLab-SMU/clusterProfiler/issues/217 - gene with
minimal ES value (NES < 0) will be reported in core_enrichment
(2019-07-31, Wed)

Changes in version 3.11.1

- build_Anno now compatible with tibble (2019-05-28, Tue)

[DropletUtils](/packages/DropletUtils)
------------

Changes in version 1.8.0

- Added the downsampleBatches() function for convenient downsampling of
  batches.

Changes in version 1.6.0

- Switched emptyDrops() to use Boost's discrete_distribution for
  weighted sampling.  This results in some minor stochastic changes to
  the Monte Carlo p-values. Automatically round non-integer count
  matrices.

[edgeR](/packages/edgeR)
-----

Changes in version 3.28.0

- Add head() and tail() methods for edgeR classes.

- Remove the 'mixed.df' argument and add a 'locfit.mixed' option to
  'trend.method' in estimateDisp() and WLEB().

- Add two new arguments 'large.n' and 'min.prop' to filterByExpr() to
  allow users to change parameters previously hard-wired.

- Remove 'values' and 'col' arguments to plotMD.DGELRT() and
  plotMD.ExactTest() as no longer needed because of changes to
  plotWithHighlights().

- roast.DGEList() and mroast.DGEList() now pass the 'nrot' argument to
  roast.default().

- Rename dglmStdResid() to plotMeanVar2().

- getDispersions() is no longer exported.

- Estimated dispersions are now numeric even if NA.

- Bug fix to goana.DGELRT() and kegga.DGELRT() when the LRT was on more
  than 1 df.

[EnhancedVolcano](/packages/EnhancedVolcano)
---------------

Changes in version 1.4

- modified behaviour where a p-value of 0 is found: now converts these
  to 10^-1 * lowest non-zero p-value

- added new paremeter 'legendLabels', which allows user to use
  expressions in the legend label

- added support for tibbles

- transcriptPointSize, transcriptLabSize, transcriptLabCol,
  transcriptLabFace, transcriptLabhjust, and transcriptLabvjust now
  deprecated. Use pointSize, labSize, labCol, labFace, labhjust, and
  labvjust, respectively, instead

- pointSize default changed to 2.0

- boxedlabels now deprecated. Use boxedLabels

[enrichplot](/packages/enrichplot)
----------

Changes in version 1.5.2

- update node_label parameter in cnetplot to support selection of
subset to be labeled (2019-09-27, Fri) -
https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#fig:cnetNodeLabel
- upsetplot for gseaResult (2019-09-25, Wed) - reimplement upsetplot
based on ggupset

Changes in version 1.5.1

- gseadist for plotting logFC distribution of selected gene sets.
(2019-06-25, Tue)

[ensemblVEP](/packages/ensemblVEP)
----------

Changes in version 1.28.0

- add support for Ensembl release 97/98

[epivizrServer](/packages/epivizrServer)
-------------

Changes in version 999.999

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrServer>

[ExperimentHub](/packages/ExperimentHub)
-------------

Changes in version 1.11.0

NEW FEATURES

- (1.11.2) system environment variable to control localHub option for
  creating hub based only on previously downloaded resources

- (1.11.5) With change in AnnotationHub. All force redownload of Hub
  sqlite file with refreshHub

BUG FIXES

- (1.11.2) Fix localHub when no internet connection.  The internal use
  of isDevel was preventing Hub creation when no internet connection.
  Fixed by checking connection. This code pretained to orgDb filters

[FamAgg](/packages/FamAgg)
------

Changes in version 1.13.1

- Parameters id.col, father.col, mother.col, family.col and sex.col for
  FAData allow now to specify the column names for the respective
  content also if a data.frame is provided.

[fcScan](/packages/fcScan)
------

Changes in version 0.99

- Fixing bugs and preparing for submission

[FELLA](/packages/FELLA)
-----

Changes in version 1.5.3

- Changed `buildDataFromGraph()` so that it looks for organismal
  annotations in `keggInfo()`.

- If `ncbi-geneid` was not available, `buildDataFromGraph()` would
  crash with 404. Now, it can use `ncbi-proteinid` if `ncbi-geneid` is
  missing.

- Discovered by LY by building db for `"cvi"`

- Small update on `sanitise()`

Changes in version 1.4.2

- **The REST KEGG service changed and modules are no longer listable by
  organism**. `FELLA` now chooses modules that have at least one
  organismal gene. This seems equivalent to picking the modules from
  `keggLink("genome", "module")`, but the latter is slow (90s).

Changes in version 1.4.1

- Fixed bug in vignette due to changes in `biomaRt`

[fgsea](/packages/fgsea)
-----

Changes in version 1.11.2

- Simpler handling of conditional probabilities

- Added the exact algorithm to inst folder

Changes in version 1.11.1

- Proper absEps handling

[fishpond](/packages/fishpond)
--------

Changes in version 1.2.0

- Switching to a faster version of Swish which only computes the ranks
  of the data once, and then re-uses this for the permutation
  distribution. This bypasses the addition of uniform noise per
  permutation and is 10x faster. Two designs which still require
  re-computation of ranks per permutation are the paired analysis and
  the general interaction analysis. Two-group, stratified two-group,
  and the paired interaction analysis now default to the new fast
  method, but the original, slower method can be used by setting fast=0
  in the call to swish().

- Adding Rcpp-based function readEDS() written by Avi Srivastava which
  imports the sparse counts stored in Alevin's Efficient Data Storage
  (EDS) format.

- Changed the vignette so that it (will) use a linkedTxome, as sometime
  the build would break if the Bioc build machine couldn't access
  ftp.ebi.ac.uk.

- Add 'computeInfRV' function. InfRV is not used in the Swish methods,
  only for visualization purposes in the Swish paper.

- removed 'samr' from Imports, as it required source installation,
  moved to Suggests, for optional qvalue calculation

[flowcatchR](/packages/flowcatchR)
----------

Changes in version 1.20.0

Bug fixes

- In the shinyFlow app, include.area now defaults to FALSE

[flowDensity](/packages/flowDensity)
-----------

Changes in version 1.19.8

- -Adding densityoverlay option to plotDens, and changing sweave to
  markdown

[flowSpecs](/packages/flowSpecs)
---------

Changes in version 0.99.4 (2019-09-30)

- responses to BioConductor review #1

Changes in version 0.99.3 (2019-09-18)

- Submitted to BioConductor

Changes in version 0.99.2 (2019-09-18)

- Diminishing the size of the example data, to fit all size
  requirements.

Changes in version 0.99.1 (2019-09-18)

- Version change, to suit BioConductor.

Changes in version 0.9.2 (2019-09-17)

- Formal tests included for all user functions.

- Using specMatCalc with one color is deprecated.

Changes in version 0.9.1 (2019-09-16)

- The name "flowSpecs" is introduced.

- The older "theFlowSpec" is from now on deprecated.

[flowSpy](/packages/flowSpy)
-------

Version: 2018-12-19
Text:

Version: 2019-03-31
Text:

Version: 2019-05-04
Text:

Version: 2019-05-07
Text:

Version: 2019-05-29
Text:

Version: 2019-08-08
Text:

[FoldGO](/packages/FoldGO)
------

Changes in version 1.2.3 (2019-07-17)

- issue with wrong amount of tests fixed

Changes in version 1.2.2 (2019-07-16)

- Multiple testing procedure corrected

Changes in version 1.2.1

- Bug with wrong arguments order in plotting functions fixed

[GCSscore](/packages/GCSscore)
--------

Changes in version 1.0.0

- All GCSscore probe packages are now automatically generated from
  Bioconductor sources (platform design (pd) and annotation (.db)
  packages) by using the makeProbePackage() function from the
  'AnnotationForge' package.  These are generated and installed based
  on the chip-type being analyzed by the end user on an as-needed
  basis.  NEW FEATURES

- Initial release of R package.

- All GCSscore probe packages are now automatically generated from
  Bioconductor sources (platform design (pd) and annotation (.db)
  packages) by using the makeProbePackage() function from the
  'AnnotationForge' package.  These are generated and installed based
  on the chip-type being analyzed by the end user on an as-needed
  basis.

[GDSArray](/packages/GDSArray)
--------

Changes in version 1.5.3

NEW FEATURES

- Dollar completion for 'GDSFile' class.

BUG FIXES

- dimnames() of GDSArray and GDSArraySeed is now list of all character
  vectors.

[gdsfmt](/packages/gdsfmt)
------

Changes in version 1.22.0

NEW FEATURES

- a new function `unload.gdsn()` to unload a GDS node from memory

UTILITIES

- add '#pragma GCC optimize("O3")' to some of C++ files when GCC is
  used

- add the compiler information in `system.gds()`

- change the file name "vignettes/gdsfmt_vignette.Rmd" to
  "vignettes/gdsfmt.Rmd", so `vignette("gdsfmt")` can work directly

BUG FIXES

- avoid the segfault if the data type is not registered internally

- use O_CLOEXEC (the close-on-exec flag) when open and create files to
  avoid potentially leaking file descriptors in forked processes

[gemini](/packages/gemini)
------

Changes in version 0.3.0

- 6-01-19: Introducing GEMINI, a variational Bayesian approach to
analyze pairwise CRISPR screens.

- Note: This is a pre-release version of GEMINI.

- Added a NEWS.md file to track changes to the package.

- Minor modifications were made to prepare for repository submission.

- 6-13-19: Incremented version number from 0.2.0 to 0.99.0 for
pre-release.

- 6-24-19: Bioconductor revision in progress. Version 0.99.9

- Note: This is a pre-release version of GEMINI.

- Modified documentation for all functions and added a workable
vignette

[geneClassifiers](/packages/geneClassifiers)
---------------

Changes in version 1.9.0

- Bugfix: Added signature( "FixedExpressionData", i=ANY, j=missing )
  for the "[" and "[[" functions

- More indepth checking of function arguments

- Additional unit tests

[GENESIS](/packages/GENESIS)
-------

Changes in version 2.15.3

- Add option "fastSKAT" to assocTestAggregate. Some other arguments and
  names of output columns for SKAT tests have also changed. This update
  includes code from the bigQF package
  (https://github.com/tslumley/bigQF). With the addition of C code,
  GENESIS now requires compliation.

[GenomicAlignments](/packages/GenomicAlignments)
-----------------

Changes in version 1.22.0

NEW FEATURES

- Add stackStringsFromGAlignments(). Analog to stackStringsFromBam()
  except that it stacks the read sequences stored in a GAlignments
  object instead of a BAM file.

BUG FIXES

- Fix summarizeJunctions() error when no junctions are found and
  'genome' is specified.

[GenomicFeatures](/packages/GenomicFeatures)
---------------

Changes in version 1.38.0

NEW FEATURES

- Small improvement to exonicParts() and intronicParts(): - When
  'linked.to.single.gene.only' is set to TRUE, exonicParts() now
  returns an additional "exonic_part" metadata column that indicates
  the rank of each exonic part within all the exonic parts linked to
  the same gene. This is for compatibility with old disjointExons().  -
  intronicParts() does something similar except that the additional
  metadata column is named "intronic_part".

SIGNIFICANT USER-VISIBLE CHANGES

- makeTxDbFromGRanges() and makeTxDbFromGFF() now recognize more
  features as transcripts, exons, or CDS (see commits 822665f8 and
  6281f856)

BUG FIXES

- Fix makeTxDbFromUCSC(..., "refGene") for
  bosTau9/galGal6/panTro6/rheMac10

[GenomicOZone](/packages/GenomicOZone)
------------

Changes in version 0.99.9

- Added clutering method: Multi-channel weighted univariate clustering.

- In function MD.Chr.zoning.Granges(), if only using one core, stop
  calling parallel computing.

- Correct a bug in function MD.rank.statistic(). ANOVA requires every
  group has no less than 2 elements. If not, report p-value 1 directly.

- Fixed the mal-formated NEWS file.

Changes in version 0.99.8

- Solved problem connecting to GitHub and bioconductor git server.

Changes in version 0.99.7

- Solved problem connecting to GitHub and bioconductor git server.

Changes in version 0.99.6

- Solved problem connecting to GitHub and bioconductor git server.

Changes in version 0.99.5

- Solved problem connecting to GitHub and bioconductor git server.

Changes in version 0.99.4

- Modified R document files. + Corrected typos. + Replaceded
  "expression" with "activity", because not only expression can be
  analyzed, but also all kinds of gene activity. + Replaced the
  "GenomicOZone list" with "GenomicOZone object".

Changes in version 0.99.3

- Replaced c(1:nrow()) into seq_len(nrow()) to avoid potential issues.

- Modified the vignettes to avoid long lines in code chunks.

- Renamed parameter 'p.value.cutoff' with 'alpha'.

- Renamed parameter 'effect.size.rate' with 'min.effect.size'.

- Updated .bib files. Added missing citations. Removed duplicated
  citations.

- Replace all c(1:...) with seq_len().

- Replaced the \texttt in vignettes Rmd file with a correct format
  using "``".

- Replaced the \textit in vignettes Rmd file with a correct format
  using "**".

- Adjusted the image size in vignettes Rmd file.

- Updated R documents.

Changes in version 0.99.2

- Submitted on 2019-08-18

- Added joemsong as a coresponding auther to the submission.

Changes in version 0.99.1

- Submitted on 2019-08-17

- Reformatted the package for Bioconductor submission.

Changes in version 0.99.0

- Submitted on 2019-08-17

- Packed the completed package.

[GenomicRanges](/packages/GenomicRanges)
-------------

Changes in version 1.38.0

NEW FEATURES

- GPos objects now exist in 2 flavors: UnstitchedGPos and StitchedGPos
  GPos is now a virtual class with 2 concrete subclasses:
  UnstitchedGPos and StitchedGPos. In an UnstitchedGPos instance the
  positions are stored as an integer vector. In a StitchedGPos
  instance, like with old GPos instances, the positions are stored as
  an IRanges object where each range represents a run of consecutive
  positions. This is analog to the IPos/UnstitchedIPos/StitchedIPos
  situation. See ?GPos for more information. Old serialized GPos
  instances can be converted to StitchedGPos instances with
  updateObject().

- GPos objects now can hold names

- Coercion to GPos now propagates the names

- Add GRangesFactor class (Factor derivative). See ?GRangesFactor

SIGNIFICANT USER-VISIBLE CHANGES

- Export from_GPos_to_GRanges()

- Some reorgnization of the GenomicRangesList hierarchy (see commit
  f988a5a9).

- Swap order of arguments 'seqlengths' and 'seqinfo' of the GRanges()
  constructor so now the latter comes before the former.

DEPRECATED AND DEFUNCT

- Remove findOverlaps, seqnames, and seqinfo<- methods for RangedData
  objects. These methods were deprecated in BioC 3.8 and defunct in
  BioC 3.9.

BUG FIXES

- Coercion from RangesList to GRanges is more robust to seqlevel
  differences

- Fix bug in isSmallGenome() (introduced by change in sum() in R >=
  3.5)

[GenomicScores](/packages/GenomicScores)
-------------

Changes in version 1.10.0

USER VISIBLE CHANGES

- Added support to latest release 3.0 of gnomAD MAF data, stored in the
  package MafDb.gnomAD.r3.0.GRCh38.

- Individual allele frequencies can be now retrieved from MafDb.*
  packages when 'ref' and 'alt' arguments are given to the functions
  'gscores()' and 'score()'. See manual pages and vignette for further
  details.

- NonSNRs are now searched giving the argument type="equal" to
  findOverlaps(). This means that only scores from exact matches to
  nonSNRs are returned.

BUG FIXES

- Bugfix on the 'getGScores(') function that precluded accessing the
  files downloaded by the AnnotationHub

- Bugfix in accessing MAF values from nonSNVs when multiallelic
  variants are stored in different records from the VCF file.

[ggtree](/packages/ggtree)
------

Changes in version 1.99.1

- bug fixed of geom_hilight for tree$edge.length = NULL (2019-10-16,
Wed) -
https://groups.google.com/d/msg/bioc-ggtree/GULj-eoAluI/Llpm-HbfCwAJ
- fortify method for igraph (only work with tree graph) (2019-09-28,
Sat) - ggdensitree (2019-09-11, Wed) -
https://github.com/YuLab-SMU/ggtree/pull/253 -
https://github.com/YuLab-SMU/ggtree/pull/255 -
https://yulab-smu.github.io/treedata-book/chapter4.html#visualize-a-list-of-trees

Changes in version 1.99.0

- prepare for ggtree v=2.0.0

Changes in version 1.17.5

- fortify methods for hierarchical clustering objects, including
agnes, diana and twins (2019-08-30, Fri) - now geom_hilight supports
unrooted and daylight layouts (2019-08-28, Wed) - by calling
geom_hilight_encircle - update geom_motif according to the change of
gggenes and allow labeling genes (2019-08-27, Tue) -
https://yulab-smu.github.io/treedata-book/chapter11.html#genome-locus
- re-implement geom_strip with more robust support of labelling
strip, either input taxa using name or id.  - support phylog defined
in ade4 package (2019-08-14, Wed) -
https://yulab-smu.github.io/treedata-book/chapter9.html#phylog

Changes in version 1.17.4

- now geom_cladelabel supports unrooted and daylight layouts
(2019-08-14, Wed) - by integrating geom_cladelabel2 - defined nodelab
method for ggtree to convert node number to label (2019-08-09, Fir) -
redefined nodeid as S3 generic in tidytree v=0.2.6 - change the
original function as a method for ggtree - move the nodeid function
for tree object to treeio - defunct gzoom function - introduce
rootnode parameter in geom_tree with default = TRUE and behave as
previous version (2019-08-08, Thu) - the invisible root to itself
line segment have advantage for the number of line segments is
consistent with the number of nodes.  - if rootnode = FALSE, there
will be no line segment of root to itself.  - extend gheatmap to
support collapsed node (2019-08-06, Tue) -
https://github.com/GuangchuangYu/ggtree/pull/243 - support hclust and
dendrogram (2019-07-31, Wed)

Changes in version 1.17.3

- remove re-export treeio parser function, user now need to load
treeio explictly (2019-07-24, Wed) - export layout_circular,
layout_fan and layout_rectangular - layout_dendrogram and
theme_dendrogram -
https://yulab-smu.github.io/treedata-book/chapter10.html#dendrogram -
scale_x_range for adding second x-axis for geom_range (2019-07-23,
Tue) - change branch.length parameter to center for geom_range

Changes in version 1.17.2

- extend expand according to the change of collapse (2019-07-11, Thu)
- mode parameter in collapse - geom_tiplab now works with 'circular'
and 'fan' layouts (2019-07-05, Fri) - geom_inset for adding subplots
to specific nodes (see also the inset function introduced in v=1.3.8)

Changes in version 1.17.1

- facet_data to extract data used in facet_plot or geom_facet
(2019-07-02, Tue) - continuous parameter in geom_tree to to
continuous color edge from parent to child (2019-09-25, Tue) -
https://yulab-smu.github.io/treedata-book/chapter4.html#color-tree -
root.position parameter for fortify and ggtree (2019-05-27, Mon) -
geom_facet, a geom layer version of facet_plot (2019-05-23, Thu) -
update scale_x_ggtree, now we can use gheatmap() + scale_x_ggtree()
(2019-05-22, Wed) - extend xlim_expand to work with ggplot2
(2019-05-20, Tue) -
https://yulab-smu.github.io/treedata-book/chapter9.html#xlim_expand -
add legend_title variable in gheatmap (2019-05-16, Thu)

[GNET2](/packages/GNET2)
-----

Changes in version 1.1.3

- Fix module capatiblality issue.

Changes in version 1.1.2

- Update plot format.

Changes in version 1.1.1

- Fix several conditions that may cause error in tree construction.

[GOfuncR](/packages/GOfuncR)
-------

Changes in version 1.5.2

USER-LEVEL CHANGES

- update GO-graph (version 07-Oct-2019)

Changes in version 1.5.1

NEW FEATURES

- add refine() function to restrict results from enrichment analysis to
  more specific categories

[graphite](/packages/graphite)
--------

Changes in version 1.31.1 (2019-10-24)

- Updated all pathway data.

- Removed HumanCyc pathways (database now requires subscription).

[gscreend](/packages/gscreend)
--------

Changes in version 0.99.4

- Replaced parallel::mclapply() with BiocParallel::bplapply()

- Updated stop/warning/message functions. Testing if sampling
  timepoints are named correctly.

[GSEABenchmarkeR](/packages/GSEABenchmarkeR)
---------------

Changes in version 1.6.0

- New function `evalTypeIError` for type I error rate evalution by
  sample permutation: - evaluation of one or more enrichment methods on
  one or more expression datasets - support for splitting permutations
  into blocks of defined size, and invoking parallel evaluation of the
  partitions

- New function `evalRandomGS` for evaluation of random gene sets: -
  estimates proportion of rejected null hypotheses (= fraction of
  significant gene sets) of an enrichment method when applied to random
  gene sets of defined size - evaluation of one or more enrichment
  methods on an expression dataset of choice

- New argument `method` to the `evalRelevance` function for the
  evaluation of phenotype relevance of gene set rankings, choices
  include: - "wsum": computes a weighted sum of the relevance scores
  (default), - "auc": performs a ROC/AUC analysis based on the ROCR
  package, - "cor": computes a standard correlation measure such as
  Spearman's rank correlation, - a user-defined function for customized
  behaviors.

- New function `metaFC` for summarizing fold changes of individual
  datasets across a compendium of expression datasets

- New functions `plotDEDistribution` and `plotNrSamples` for exploring
  differential expression and sample size across a compendium of
  expression datasets

- Extended support for user-defined benchmarking inputs including
  simplified plug-in of user-defined enrichment methods (thanks to
  Marcel Ramos @LiNk-NY)

[GSVA](/packages/GSVA)
----

Changes in version 1.34

BUG FIXES

- Bugfix to handle when parallel::detectCores() returns NA instead of
  an integer number of cores, which may happen when running GSVA in a
  docker container. Bug reporting and pull request fix thanks to Aaron
  (https://github.com/rcastelo/GSVA/pull/10).

- Bugfix to handle when arguments 'method="ssgsea"' and 'tau=0'. Bug
  reporting thanks to Lena Morill
  (https://github.com/rcastelo/GSVA/issues/4).

[gtrellis](/packages/gtrellis)
--------

Changes in version 1.17.1

- fixed a bug of selecting chromosomes.

[Gviz](/packages/Gviz)
----

Changes in version 1.29.1

NEW FEATURES

- Exposed the possibility to specify width and distance for feather
  (arrowBar) in AnnotationTrack class

BUG FIXES

- Fixed issue with `subseq` function for `ReferenceSequenceTrack`

- Changed the visualized position of tickmarks, and values in
  `DataTrack` => align all to + 0.5 position, which matches the
  `SequenceTrack` and `AnnotationTrack` visualization

- Changed the check for transparency support in currently opened device
  `supportsAlpha`, point moved from center to the left bottom corner

[HDF5Array](/packages/HDF5Array)
---------

Changes in version 1.14.0

NEW FEATURES

- Add coercions from TENxMatrix (or TENxMatrixSeed) to dgCMatrix

SIGNIFICANT USER-VISIBLE CHANGES

- h5mread() argument 'starts' now defaults to NULL

BUG FIXES

- h5mread() now supports datasets with contiguous layout (i.e. not
  chunked)

[HDTD](/packages/HDTD)
----

Changes in version 1.19.1 (2019-10-24)

- Updated CITATION FILES.

- Fixed minor bug.

[hiAnnotator](/packages/hiAnnotator)
-----------

Changes in version 1.19.0

- Drop forced strand conversion to '+' and makeGRanges with add strand
  as '*' if none was found.

[HIBAG](/packages/HIBAG)
-----

Changes in version 1.22.0

- change the file name "vignettes/HIBAG_Tutorial.Rmd" to
  "vignettes/HIBAG.Rmd", so `vignette("HIBAG")` can work directly

[HiLDA](/packages/HiLDA)
-----

Changes in version 0.99.11 (2019-09-10)

- Update the CITATION after the paper has been accepted by PeerJ

Changes in version 0.99.10 (2019-07-24)

- Fix the issues commented by the Bioconductor reviewer

Changes in version 0.99.9 (2019-07-24)

- Fix the issues commented by the Bioconductor reviewer

Changes in version 0.99.8 (2019-06-22)

- Remove git tracking .o and .dll files

Changes in version 0.99.7 (2019-06-22)

- Remove git tracking .o and .dll files

Changes in version 0.99.6 (2019-06-22)

- Remove git tracking .o and .dll files

Changes in version 0.99.5 (2019-06-22)

- Remove git tracking .o and .dll files

Changes in version 0.99.4 (2019-06-22)

- Remove git tracking .o and .dll files

Changes in version 0.99.3 (2019-06-21)

- Fix the issues commented by the Bioconductor reviewer

Changes in version 0.99.2 (2019-06-05)

- Remove pmsignature dependency

Changes in version 0.99.1 (2019-06-04)

- Submitted to Bioconductor

[hipathia](/packages/hipathia)
--------

Changes in version 2.1.1 (2019-05-17)

- Adding function mgi_from_sif, which allows to create a pathways
  object from SIF + ATT files.

[HIREewas](/packages/HIREewas)
--------

Changes in version 1.3.1

- Added the citation file.

[HPAanalyze](/packages/HPAanalyze)
----------

Changes in version 1.3

- Changes in version 1.3.3
    + Fix sn error introduced in version 1.3.2 where hpaVisPatho plotted incorrectly.
    + Added a vignette with codes for figures in the HPAanalyze manuscript.

- Changes in version 1.3.2
    + To reduce dependency, hpaVis now use gridExtra for multiple plots instead of cowplot
    + Removed dependency on reshape2

- Changes in version 1.3.1
    + hpaExport is now supporting csv and tsv.
    + Support for JSON export is added via vignette.
    + Fix the error where hpaXmlGet give an error when curl is not installed.
    + Most functions with now take both HGNC gene names and ensemble ids.
    + hpaVisSubcell now allows user to choose reliability scores.
    + Removed dependency on magrittr

- Changes in version 1.3.0
    + Starting devel for Bioconductor 3.10

[ideal](/packages/ideal)
-----

Changes in version 1.10.0

New features

- ggplotCounts gains a new parameter, labels_repel, to control the
placement of the different labels in the plot - this can be useful
when a large number of samples is available

Bug fixes

- Fixed an error in the initialization of the app due to a new
behavior introduced by shinyAce in version >= 0.4.0 - occurred in the
same way as for pcaExplorer

Other notes

- Better class checks via is(...) as per BiocCheck suggestion

[idr2d](/packages/idr2d)
-----

Changes in version 0.99.2

- added consistent ticks and limits to IDR plot functions

- updated function documentation

Changes in version 0.99.1

- added reference to Li et al. paper

- added diagnostics plots

- changed rank order

- added IDR1D functionality

- fixed local / global IDR issue

Changes in version 0.99.0

- initial release

[IgGeneUsage](/packages/IgGeneUsage)
-----------

Changes in version 0.99.0 (2019-07-15)

- Submitted to Bioconductor

[igvR](/packages/igvR)
----

Changes in version 1.6

NEW FEATURES

- bam (Alignment) tracks now supported

- Motif logos now displayed by clicking on appropriately configured
  tracks

- pkgdown website (vignettes and man pages) at
  https://paul-shannon.github.io/igvR/index.html

- more genomes supported: hg38, hg19, hg18, mm10, bostau8, canfam3,
  ce11, danrer10, danrer11, dm6, gorgor4, panpan2, pantro4, pfal3d7,
  rn6, saccer3, susscr11, tair10

- R commands sent to the browser now return only when the Javascript
  command completes

- built with version 2.3.2 of igv.js

[immunoClust](/packages/immunoClust)
-----------

Changes in version 1.17.3

- CHANGES * bugfix in subset.immunoMeta

Changes in version 1.17.2

- NEW FEATURES * introducing immunoMeta-class on meta-clustering
  results to buildup and annotate a hierachical
  population/sub-populaton tree

- CHANGES * trail of automated annotation of meta-clusters using
  scatter-clustering is removed. The approach does not work in a
  usefull manner. Instead the immunoMeta-class is introduced providing
  methods for a manual annotation of meta-cluster. See man-pages and
  vignette for more details. * set.seed is removed from clustering
  routines. To obtain reproducable results with cell.process function
  set.seed has to set explicit before.

Changes in version 1.17.1

- CHANGES * minor bugfixes and code cleaning * minor additional options
  for immunoClust.plot/splom

[infercnv](/packages/infercnv)
--------

Changes in version 1.1.4 (2019-10-29)

- Fix reading of input annotations when some are only digits to be
  properly read as characters.

- Added checks that HMM_report_by option is compatible with
  analysis_mode option. If not, change it automatically.

- Fix reading of bayesian filtered HMM results in add_to_seurat after
  previous version changes to keep CNV ids and states scale.

Changes in version 1.1.3 (2019-09-16)

- Fix to reload checks on HMM steps. +Added new smoothing method,
  'coordinates", that smooths the per cell data using a window based on
  a base pairs distance (around 10.000.000 seems to be a good start for
  the window size) to the current gene. As the hspike does not model
  gene distances/positions on a chromosome at this time, the HMM i6
  mode is not compatible with this smoothing method and an error will
  be returned if both try to be used together.

- Added a bp distance tolerance to "merge" top CNVs that are actually
  the same CNV in different subclusters in add_to_seurat method.

- Removed the top any type of CNV field, as it is redundant to top loss
  and top dupli.

- Added text output of identified top CNVs as they the base pair
  tolerance aggregates some compared to the original HMM output.

- Added an argument "up_to_step" to stop infercnv::run() after a given
  step.

- Fix contents of @options field in infercnv_obj to store non default
  run time arguments in the same form as object creation arguments.

- Update so that HMM predictions outputs have the analysis_mode in
  their name and do not get overwritten when using different modes. Can
  also know which pairs of file to reload together now.

- Update add_to_seurat so that it checks what analysis_mode was used in
  the run based on the @options field to reload the matching HMM
  predictions or Bayesian Network filtered predictions.

- Updated "subclustering" method for sample mode so that
  @tumor_subclusters$subclusters indices have the cell names attached
  to be able to map back in add_to_seurat.

- Allow HMM steps to be resumed if needed even if steps 20/21 are done.

Changes in version 1.1.2 (2019-07-08)

- Added method to write table of wide array of predicted features from
  HMM results to file or add them as meta.data to a Seurat object if
  one is provided.

- Overhaul of save/reload system to store non default arguments and
  keep track of relevant options at each step when trying to reload
  backups. Also check for input counts matrix identity with reloaded
  one (hash at object creation time).

- Added linking of image() option useRaster to run() and plot_cnv() to
  be able to enable by default, speeding up plotting significantly.

Changes in version 1.1.1 (2019-05-20)

- Added method to sample an infercnv object to a given number of cells,
  or at a given frequency, per annotation group. This is to make it
  easy to plot figures where all annotations groups have the same
  overall height, as well as downsample very large datasets that would
  otherwise take too long to plot (while still running the analysis on
  the full data)

- Added method to plot each annotation group to a different figure and
  combine with sampling. Mostly intended to split data for larger
  datasets.

- Added support for output_format option within run() to link to
  plot_cnv() to support only writting text outputs during the analysis.

Changes in version 1.0.4 (2019-09-16)

- Fix check that contig to cluster by was found when specified.

- Added support to plot_cnv for cell groups with exactly 2 cells.

- Fix which input file type is checked.

- Made it so that plot_cnv recalculates clustering automatically if non
  null ref_contig argument is provided.

- Fix for plot_cnv() when providing multiple ref_contigs and
  cluster_by_group is False.

- Fix only 1/n genes being taken into account when using n ref_contig
  in plot_cnv.

- Fix error in file creation when using multiple ref_contig and
  cluster_by_groups=FALSE in plot_cnv.

- Bayesian filtering now preserves CNV ids in outputs

Changes in version 1.0.3 (2019-07-05)

- Fix missing dendrograms in text output when drawing figures.

- Fix path to save object to when splitting references.

- Fix file name creation when using num_ref_groups option.

- Fix reference cells indices returned from method that splits
  references in num_ref_groups when references are not sorted and at
  the beginning of the matrix.

- Fix to support of data.frame as input type for counts matrix.

Changes in version 1.0.2 (2019-05-21)

- Reduce peak memory usage.

- Fix to subclusters definition when using a sparse matrix and a non
  random trees method with no references.

Changes in version 1.0.1 (2019-05-20)

- Improved when the clustering is defined for groups when running in
  sample mode.

- Fixed support for NA to be understood as an output_format value to
  plot_cnv() in case a user only wants to generate the text outputs and
  not the plot.

- Fix ordering of cells and color bars on the heatmap and text outputs
  when cells are not sorted in the same order in the input matrix and
  the annotation file.

- Fix to (sub)cluster definition when a group only has 1 cell.

- Fix plot_cnv() to handle observations groups with only 1 cell (that
  can't be hierarchically clustered), a single reference group when no
  references are not clusterd, and a single reference group with a
  single cell.

[InPAS](/packages/InPAS)
-----

Changes in version 1.17.3

- update documentation.

Changes in version 1.17.2

- add citation of cleanUpdTSeq.

[IRanges](/packages/IRanges)
-------

Changes in version 2.20.0

NEW FEATURES

- IPos objects now exist in 2 flavors: UnstitchedIPos and StitchedIPos
  IPos is now a virtual class with 2 concrete subclasses:
  UnstitchedIPos and StitchedIPos. In an UnstitchedIPos instance the
  positions are stored as an integer vector. In a StitchedIPos
  instance, like with old IPos instances, the positions are stored as
  an IRanges object where each range represents a run of consecutive
  positions. See ?IPos for more information.  Old serialized IPos
  instances need to be converted to StitchedIPos instances with
  updateObject().

- IPos objects now can hold names

- The IRanges() and IPos() constructors now accept user-supplied
  metadata columns

- Add grep(), startsWith() and endsWith() methods for CharacterList
  objects

SIGNIFICANT USER-VISIBLE CHANGES

- as.data.frame(IRanges) now propagates the metadata columns

- Move splitAsList() to the S4Vector package

- Move S4 class "atomic" from the S4Vector package

- No longer export %in% (was a leftover from an older time when the
  package was defining an %in% method)

DEPRECATED AND DEFUNCT

- After being deprecated in BioC 3.9, the following RangedData methods
  are now defunct: findOverlaps, rownames<-, colnames<-,
  columnMetadata, columnMetadata<-, c, rbind, as.env, as.data.frame,
  and coercion from RangedData to DataFrame.

- Remove the following RangedData methods: - score, score<-, lapply,
  within, countOverlaps; - coercions from list, data.frame, DataTable,
  Rle, RleList, RleViewsList, IntegerRanges, or IntegerRangesList to
  RangedData.  These methods were deprecated in BioC 3.8 and defunct in
  BioC 3.9.

BUG FIXES

- Fix integer overflow issue in end() setter for IRanges objects.

[iSEE](/packages/iSEE)
----

Changes in version 1.5.13

- Order features selected in heat map selectize from top to bottom.

Changes in version 1.5.12

- Support gene list input from aceEditor and fileInput.

Changes in version 1.5.11

- Rename isColorMapCompatible to checkColormapCompatibility.

- Fix graceful server side handling of checkColormapCompatibility.

- Update documentation about panel organisation in vignette.

Changes in version 1.5.10

- Fix test to provide a non-empty selection to custom plot function.

Changes in version 1.5.9

- Introduce Bugs Easter egg.

Changes in version 1.5.8

- Substitute deprecated scater::normalize by logNormCounts.

Changes in version 1.5.7

- Simplify protection of redDimPlotDefaults against empty reducedDims.

- Fix to declare all panel types not available.

Changes in version 1.5.6

- Updates following deprecation of isSpike and sizeFactorNames.

Changes in version 1.5.5

- Add modeEmpty().

- Support zero-row initialPanels argument.

Changes in version 1.5.4

- Added support for file upload with server re-initialization.

- Moved observers to separate file. Exclude from code coverage.

- Updating calls to ReprocessedAllenData() to load only tophat_counts
  assay.

Changes in version 1.5.3

- Use ReprocessedAllenData() following the deprecation of data(allen).

Changes in version 1.5.2

- Minor doc fix..

- Do not allow duplicated values in Name field of initialPanels.

- Downsample points randomly.

Changes in version 1.5.1

- Fix report of table links.

Changes in version 1.5.0

- Bioconductor release.

[IsoformSwitchAnalyzeR](/packages/IsoformSwitchAnalyzeR)
---------------------

Changes in version 1.7.2 (2019-10-18)

- Update type: Major.

- analyzeIUPred2A() for analyzing intrincially disordered regions (and
  binding sites therein) was introduced. To enable this the following
  changes were also made: * analyzeNetSurfP2() was extended to also
  create the idr_type column in the result *
  analyzeSwitchConsequences() was extended to handle idr_type. Also it
  was upgrated to handle large differences in IDR lengths. * The data
  included in the "exampleSwitchListAnalyzed" object was updated to
  include the result of an IUPred2A analysis (instead of the NetSurfP2
  analysis) * The build in data file for analysis of NetSurfP-2 in
  relation to exampleSwitchListIntermediary was replaced by the
  corresponding data for the IUPred2A analysis. *
  switchPlotTranscript() (which is used by switchPlot() internally) was
  extended to also handle IDR types * the switchPlot() layout was
  re-optimied for the new annotation. * isoformSwitchAnalysisPart2()
  was updated to also handle IUPred2A input. * The vignette was updated
  accordingly.

- switchPlotTranscript() (and thereby also switchPlot) now use the
  annotationImportance in a much nicer way. Instead of removing the
  annotation (which could cause problems when comparing computational
  analysis to visual output) it now uses annotationImportance to plot
  the data as layers with the most important on top - meaning no
  annotation is skipped.

- switchPlotGeneExp() was updated to follow the condition coloring used
  by switchPlotIsoExp() and switchPlotIsoUsage() when used by the
  switchPlot() function.

- Corrected a bug in switPlot() which caused the interpretation of the
  "increased/decreased usage" added to the plot to be the min instead
  of max of the supplied alphas.

- Corrected a bug in analyzeSwitchConsequences() that could cause the
  "domain_length" consequnce type to give wrong results. Now the
  'domain_length' test transcripts for differences in the length of
  overlapping domains of the same type (same hmm_name)

- isoformSwitchAnalysisPart2() now also uses n=Inf to create all plots
  (NA have same function for backward compatability).

- importRdata now ensures the order of columns in the designmatrix is
  always.

- All functions for importing external analysis, which supports
  multiple files, now automatically remove duplicated interies.

- The extractExpressionMatrix function was depreciated.

- All function documentation was spell-checked.

- Various documentation improvements.

- Error message improvements.

Changes in version 1.7.1 (2019-07-19)

- Update type: Minor.

- Version bump due to Bioconductor release.

- Updated NEWS layout in accordance with Bioconductor guidelines.

- switchPlotTranscript() was extended to also indicate
  increased/decreased/unchanged isoform usage making interpretation
  easier. This also required switchPlotTranscript() and switchPlot()
  was updated with extra arguments to control this behaviour. The
  switchPlotTranscript() function was furthermore updated to also
  indicate significance (indated by asterisks) and size (dIF) when used
  alone (aka not from within switchPlots) making it a good alternative
  to the switchPlot.

- switchPlotGeneExp(), switchPlotIsoExp(), switchPlotIsoUsage() was
  prettyfied and now also show the name of the gene plotted.

- importRdata() and importGTF() now also supports import of RefSeq GFF
  files (downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/, see FAQ
  in vignette). This should increase ease of usage for a long range of
  species not in the Ensembl catalogue.

- importRdata() * Now also removes non-exsisting introns from
  annotation even when it is supplied as a GRange (previously only done
  for GTF files). * No longer removes isoforms with NA as biotypes when
  removeTECgenes = TRUE. * Was extended to better handle gene_names
  when novelt transcripts are predicted. Specificallt if there are NA
  in the gene_name column (e.g. like done by StringTie) these are
  automatically assigned the same gene name as the other isoforms from
  the same gene_id (only for cases where a single gene_name is
  associated to the gene_id).

- isoformSwitchTestDEXSeq() was updated to: * Better handle rare design
  setups that could cause an error to occure. * Now handle analysis of
  data with some isoforms only analyzed in a subset of comparisons

- The extractSwitchSummary() was extended to also print number of
  switches.

- A bug was fixed in extractSequence() which cased a fail when CDS
  sequences with multiple stop codons where annotated

- A bug was fiexed which caused extractSplicingSummary() to only return
  the summary of splicing types with more than "minEventsForPlotting"
  events.

- analyzeNetSurfP2() was updated to handle multiple files due to recent
  restrictions on the number of sequences one can upload to the
  webserver.

- switchPlotTopSwitches() and extractTopSwitches() now uses n = Inf to
  to output all (although internally NA is converted to Inf for
  backward compatability).

- subsetSwitchAnalyzeRlist() was improved to be more stable to edgecase
  sitiuations.

- isoformToGeneExp() was improved * To be more userfriendly. * To
  directly support annotation stored in a GTF file (which it itself
  imports into R). * To directly support switchAnalyzeRlists.

- analyzePFAM() was updated to be more robust to edge usecases

- Improved error messages in mutliple function.

- Various documentation updates.

- Various stability updates.

[kebabs](/packages/kebabs)
------

Changes in version 1.20.0

- release as part of Bioconductor 3.10

Changes in version 1.19.1

- removed change history from package vignette for easier maintenance

Changes in version 1.19.0

- new branch for Bioconductor 3.10 devel

[KEGGprofile](/packages/KEGGprofile)
-----------

Changes in version 1.27.4

- Fix bugs in vignette

- Improvements for notes and warnings in R cmd check

Changes in version 1.27.3

- Improvement for plot_pathway_overall function

Changes in version 1.27.2

- Fix bugs in convertId function, which was caused by the updates of
  biomaRt package

Changes in version 1.27.1

- Fix bugs in download_KEGGfile function, which was caused by the
  updates of KEGG web site.

- New function: plot_pathway_overall function to plot gene expression
  in pathway level.

- Improvement for documents and examples.

[KnowSeq](/packages/KnowSeq)
-------

Changes in version 0.99.55 (2019-10-09)

- Several minor bug fixes and enhancements and xgrid/ygrid added to
  classification plots Further versions

- MAC_OS alignment tools support

- Incorporation of RUV to batch effect methods

Changes in version 0.99.51 (2019-08-19)

- New method web platform added (targetValidation) to retrieve the
  diseases in the DEGsToDiseases function

Changes in version 0.99.42

- dataPlots genesBoxplot mode improvements

Changes in version 0.99.31 (2019-07-17)

- Best tunning parameters for SVM functions

Changes in version 0.99.30 (2019-06-05)

- Initial release with Bioconductor

- Alignment functions only works on Unix system (For now)

- The rest of the pipeline is completely functional for all the OS
  including windows and MAC_OS

[limma](/packages/limma)
-----

Changes in version 3.42.0

- New head() and tail() methods for all limma data classes.

- New unique() and subsetting methods for TestResults objects.
  Previously a subsetted TestResults object became an ordinary numeric
  matrix, but now the TestResults class is preserved. Single index
  subsetting is also allowed but produces a numeric vector.

- roast() and mroast() now use less memory. Memory usage now remains
  bounded regardless of the number of rotations used. Both functions
  are faster than before when approx.zscore=TRUE. A new argument
  `legacy` is added to allow users to turn these improvements off if
  they wish to reproduce numerical results from limma 3.40.0 exactly.
  The default number of rotations `nrot` is increased from 999 to 1999.
  A bug has been fixed for roast() and mroast() when `index` is a
  character vector of geneids or rownames. Previously this usage
  failed.

- Complete rewrite of zscoreT(), which is now somewhat faster and offer
  two new options for the normalizing transformation. The new argument
  `method` indicates which transformation is used when `approx=TRUE`.
  The default approximation method is changed from "hill" to "bailey".
  zscoreT() now allows NA values when `approx=TRUE` and works correctly
  for all valid df values. Previously `approx=TRUE` returned NA results
  for `df` infinite or <= 0.05. A simple robust approximation from
  Wallace (1959) has been introduced to handle very large or very small
  `df` values. Previously `approx=FALSE` treated any `df` greater than
  10000 as infinite; this threshold is now raised to 1e300.

- tZScore() is slightly faster.

- New argument `gene.weights` for fry() so as to match the arguments
  and behavior of mroast(). Previously gene weights could be input to
  fry() only through the `index` argument. fry's mixed p-values now
  take account of gene weights. Previously gene weights were not used
  in fry's mixed p-value calculation.

- Add arguments `block`, `correlation` and `weights` to voom() and
  remove the `...` argument.

- Update the arguments of voomWithQualityWeights() to match the
  revisions to arrayWeights() in limma 3.40.0.  Add new argument
  `var.group`.

- New arguments `xlab`, `ylab` and `main` for plotFB(). The default
  plot title now includes the column name from the data object.
  Default for `pch` increased to 0.3.

- plotWithHighlights and plotMD.MArrayLM now gives special treatment to
  the case that `status` is a TestResults object.

- Update Yoruba case study in User's Guide to call voom() and
  duplicateCorrelation() twice each instead of once.

- Update Users Guide: update Rsubread reference, correct remarks about
  prior distribution for 1/sigma^2 in Section 13.2 and rerun the Yoruba
  and Pasilla case studies.

- Use "RNA-seq" and "ChIP-seq" consistently in the documentation
  instead of "RNA-Seq" and "ChIP-Seq".

- All uses of the approx() function now set
  `ties=list("ordered",mean)`, except for fitFDistRobustly() which uses
  `ties=mean`.  The change was prompted by a change in R 3.6.0 whereby
  using the default for `ties` generates a warning message whenever
  ties are present in `x`. The new code gives the same results but is
  very slightly faster and avoids the warning message.  limma now
  depends on R >= 3.6.0 because of this change.

- Fix bug in barcodeplot() when `index` or `index2` are character
  vectors.

[lionessR](/packages/lionessR)
--------

Changes in version 0.99.2 (2019-09-30)

- Modify behaviour when SummarizedExperiment object used as input

- Use SummarizedExperiment object in vignettes

- Add description for the vignettes figure

- Change examples in manual for lioness function

Changes in version 0.99.1 (2019-09-05)

- Bug fixed for matrix without column names

- Accept SummarizedExperiment object as input

Changes in version 0.99.0 (2019-08-19)

- Submitted to Bioconductor

[lipidr](/packages/lipidr)
------

Changes in version 1.99.2

- Breaking changes: changed the core object to LipidomicsExperiment.

- Added support for both targeted and untargeted lipidomics analysis.

- lipidr can accept numerical matrix as input.

- Added integration with Metabolomics Workbench API for enable data
  mining.

[Maaslin2](/packages/Maaslin2)
--------

Changes in version 0.99.15 (2019-08-05)

- Add value to categorical plots.

- Adding another dependency required by bioconductor MacOS automated
  build/test

Changes in version 0.99.14 (2019-07-31)

- Adding two more dependencies required by bioconductor MacOS automated
  build/test

Changes in version 0.99.13 (2019-07-31)

- Iterations for addition to bioconductor: Add back in license file and
  small changes to coding sections of vignette.

- Modifications to man page to include new option.

Changes in version 0.99.12 (2019-07-26)

- Add new option to set the max number of features shown in heatmap.

- Fix heatmap to include all rows of significant values for the top N
  features instead of only including the rows after finding the top N
  features.

Changes in version 0.99.11 (2019-07-24)

- Update demo to data from HMP2 (provided by Himel).

Changes in version 0.99.10 (2019-07-19)

- Small modifications to documentation to update dependency install
  notes to match bioconductor.

Changes in version 0.99.9 (2019-07-18)

- Add dependency to namespace for automated tests.

Changes in version 0.99.8 (2019-07-18)

- Add one more dependency for bioconductor MacOS build tests

Changes in version 0.99.7 (2019-07-17)

- Change file paths to relative to package for windows tests.

Changes in version 0.99.6 (2019-07-17)

- Fix format of R sections in vignette to pass tests.

Changes in version 0.99.5 (2019-07-17)

- Update R sections of vignette format.

- Change test paths for windows.

Changes in version 0.99.4 (2019-07-17)

- Update required R version.

- Modifications in vignette format for bioconductor build.

Changes in version 0.99.3 (2019-07-17)

- Modifications to package based on feedback from bioconductor review
  (additions to description, notes on data files, update testing to
  testthat, condensing vignette sections and updating install to
  bioconductor method, remove tests for packages, use seq_len and
  lapply)

Changes in version 0.99.2 (2019-06-27)

- Show top N features in heatmap instead of top N associations

Changes in version 0.99.1 (2019-06-05)

- Fix plots to allow for NAs in values

Changes in version 0.99.0 (2019-05-24)

- Only show the top 50 associations in the heatmaps

- Use static heatmap plot colors

- In boxplots, use angle for x axis text for long text strings (if any
  in set is over 5 chars)

- For larger y axis labels reduce the font size (if over 15 chars)

- Add Ns to plots in annotation for continuous and x axis label for
  categorical

Changes in version 0.3.0 (2019-05-20)

- Plots now show normalized/filtered/transformed data

- Package modifications for submission to bioconductor

Changes in version 0.2.3 (2018-12-20)

- Move filtering to after normalization

- Updates to barplots

Changes in version 0.2.2 (2018-11-15)

- Fix issue with single column in visualizations (Thanks, sma!)

- Add hash to dependencies

- Change output column names to match data.frame names

- Add options to bypass plotting

- Add crossed random effects for LM using lme4 and lmerTest

- Fix ZICP fitting errors

- Add stderr to results

- Rotate heatmap column names by 45 degrees

Changes in version 0.2.1 (2018-10-10)

- Update read/slicing to support input files with a single feature.

Changes in version 0.2.0 (2018-10-09)

- Group boxplots/scatter plots by metadata name.

- Replace ggsave with pdf to print heatmap/plots to resolve ggsave
  Rplot.pdf issue.

- Add tryCatch to allow for error in heatmap but still print other
  plots.

- Allow data/metadata inputs to be paths to files or data.frame.

- Return fit data from maaslin2 function.

- Set na.action default in model fit to na.exclude.

Changes in version 0.1.0 (2018-09-27)

- Initial tagged release.

[MACSQuantifyR](/packages/MACSQuantifyR)
-------------

Changes in version 0.99.0 (2018-11-14)

- Not submitted o GUI interface o Sort o Statistics o Graphics

[maftools](/packages/maftools)
--------

Changes in version 2.2.0

NEW FUNCTIONS AND FEATURES

- `survGroup`, `mafSurvGroup` - Predict genes/genesets associated with
  survival. Issue: #396

- New argument `altered` in `oncoplot` to plot top genes based on CNV
  or mutation. Default is `FALSE`. Issue: #405

DEPRECATED

- `oncotate` - Oncotator is no longer supported by Broad. Issue: #403
  #384 #381

- `findPathways` argument in `somaticInteractions` function has been
  deprecated

SIGNIFICANT USER-LEVEL IMPROVEMENT

- Signature analysis has been modified to be more flexible. It now
  consits of three functions namely:

- `estimateSignatures` - which measures cophenetic correlation metric
  (a measure of goodness of fit) for a range of values

- `plotCophenetic` - which draws an elbow plot of cophenetic
  correlation metric from `estimateSignatures` and helps to decide an
  optimal number of signature `(n)`.

- `extractSignatures` - which then extracts final `n` signatures

- `compareSignatures` now has two databases of known signatures.
  Classic 30 signatures and newer 67 SBS signatures from COSMIC.

BUG FIX

- `gisticChromPlot` bug fix. Issue: #392

- `annovarToMaf` bug fix for results without exonic variants. Issue:
  #388

- Avoid conflict with `plotly::layout` with `graphcis::layout`. Issue:
  #387 #202

- Fix side-bar heights in Oncoplot with copy-number data. Issue: #383

- Update y-axis label for `tcgaCompare` plot. Issue: #366

- `annovarToMaf` annotating variants with MNPs. Issue: #335

[MAGeCKFlute](/packages/MAGeCKFlute)
-----------

Changes in version 1.4.3

- Remove bugs when perform enrichment analysis based on use-defined
  gene sets.

- Prioritize NormalizeBeta.

- Customize KEGG pathways (04***).

Changes in version 1.4.1

- Prioritize many functions.

[matter](/packages/matter)
------

Changes in version 1.11.7 (2019-10-25)

BUG FIXES

- Fixed errors in 'rowStats()' and 'colStats()' when calculating group
  statistics over minor dimensions

Changes in version 1.11.6 (2019-10-23)

NEW FEATURES

- Added 'locmax()' function for finding local maxima

- Added 'binvec()' function for binning vectors

Changes in version 1.11.5 (2019-10-13)

NEW FEATURES

- Exposed 'chunk_apply()' for applying functions over parallelized
  chunks of vectors and matrices

Changes in version 1.11.4 (2019-10-13)

NEW FEATURES

- Added 'rowStats()' and 'colStats()' methods for chunk-applying
  grouped summary statistics over the rows and columns of matrices

SIGNIFICANT USER-VISIBLE CHANGES

- Family of 'apply()' and 'lapply()' methods now attempt to respect the
  object's 'chunksize()'

- Added getOption('matter.default.chunksize') for setting default
  chunksize for matter objects

Changes in version 1.11.3

SIGNIFICANT USER-VISIBLE CHANGES

- Faster subsetting of matter-backed ALTREP objects

- Updated 'bsearch()' to return the closest match when tol > 0

- Added getOption('matter.dump.dir') to control where temporary files
  are stored for matter objects

Changes in version 1.11.2

NEW FEATURES

- Coercing to native R types now returns an ALTREP representation for
  most 'matter' objects

- Use 'as.altrep()' method to coerce an existing 'matter' object to
  ALTREP or to coerce native R types to 'matter'-backed out-of-memory
  ALTREP objects

- ALTREP coercion can be controlled by new options:
  getOption('matter.coerce.altrep')
  getOption('matter.coerce.altrep.list')
  getOption('matter.wrap.altrep')

- See ?`matter-options` for details

SIGNIFICANT USER-VISIBLE CHANGES

- The 'show()' method for 'matter' objects now prints a preview of the
  data head

- Printing of data can be controlled with new options
  getOption('matter.show.head') and getOption('matter.show.head.n')

Changes in version 1.11.1

NEW FEATURES

- Added 'stream_stat' class for streaming statistics with functions
  s_mean(), s_var(), s_sd(), etc.

[MEB](/packages/MEB)
---

Changes in version 0.99.3

- The latest one.

Changes in version 0.99.2

- Update package to build NAMESPACE with roxygen and use seq_len() in
  MEB function.

Changes in version 0.99.1

- Update package.

Changes in version 0.99.0

- Update package to transform data format so that the data inputs and
  outputs are consistent with standard Bioconductor representations
  such as SummarizedExperiment.

[MetaboSignal](/packages/MetaboSignal)
------------

Changes in version 1.14.1

- MS_keggNetwork returns interaction subtype (e.g. activation) instead
  of type (e.g. PPrel). Thank you Shilpa Harshan for noticing this!

[metagene2](/packages/metagene2)
---------

Changes in version 1.1.4 (2019-07-12)

- Improved documentation of control samples in design.

Changes in version 1.1.3 (2019-05-13)

- When parameter validation fails, the metagene2 object will now fall
  back to the previous set of parameters, instead of ending up in an
  invalid state.

Changes in version 1.1.2 (2019-05-09)

- rnaseq mode now sets strand_specific to TRUE by default.

Changes in version 1.1.1 (2019-05-03)

- Added log2_ratio normalization option.

- Fix a bug that prevented strand_specific and stitch mode to work
  together.

Changes in version 1.1.0 (2019-04-05)

- Submitted package to Bioconductor.

[metaMS](/packages/metaMS)
------

Changes in version 1.21.4

- Pull request accepted bug correction in annotation pipeline

Changes in version 1.21.2

- change Maintainer

[MetCirc](/packages/MetCirc)
-------

Changes in version 1.15.2 (2019-09-09)

- add codecov

- add Travis-CI for continuous integration

Changes in version 1.15.1 (2019-08-29)

- add ggplot2 in dependencies

- change GPL-2 to GPL-3 license

Changes in version 1.15.0 (2019-04-23)

- implement the MSnbase Spectra/Spectrum2 as the container for MS2
  spectra, change all functions that they accept Spectra objects

- write function convertMsp2Spectra that converts MSP files to Spectra
  files

- change data files: change convertMSP2MSP.RData to convertMsp2Spectra,
  create spectra.RData, update similarityMat.RData

- remove data files: binnedMSP.RData, idMSMStoMSP.RData

[MethCP](/packages/MethCP)
------

Changes in version 0.99.0

- Submitted to Bioconductor

[methrix](/packages/methrix)
-------

Changes in version 0.99.0

o Submission to Bioconductor

[methylGSA](/packages/methylGSA)
---------

Changes in version 1.3.5

- Using IlluminaHumanMethylationEPICanno.ilm10b4.hg19

Changes in version 1.3.2

- Shiny app is available within the package

Changes in version 1.3.1

- Bug fixes in methylgometh

[MGFR](/packages/MGFR)
----

Changes in version 1.10.1

- fixed the error that caused failed build of the package as follows:
  Changed the argument "entrezgene" used in getBM() to "entrezgene_id",
  since this was changed in the Ensembl database and caused an error in
  the internal function '.get.genes.rnaseq()'.

[microbiome](/packages/microbiome)
----------

Changes in version 2.0 (2019-08-20)

- spreadplot function added

- removed ready made themes from functions

- Added is.compositional

- Fixed a bug in core_members (also non-compositional detection now
  allowed)

- removed rm.na option from aggregate_taxa

- Deprecating noncore_* functions (replacing with rare_* functions
  everywhere)

- Removed variable_members function

- Support removed from R-3.3.3 and lower

[microbiomeDASim](/packages/microbiomeDASim)
---------------

Changes in version 0.99.0 (2019-08-26)

- Submitted to Bioconductor

[miRSM](/packages/miRSM)
-----

Changes in version 1.3.2

- Update miRSM function <2019-09-11, Wed>

Changes in version 1.3.1

- Add modular analysis <2019-08-12, Mon>

[miRspongeR](/packages/miRspongeR)
----------

Changes in version 1.11.1

- Add citation <2019-05-10, Fri>.

[mixOmics](/packages/mixOmics)
--------

Changes in version 6.8.6

new features / enhancements

-

bug fixes

- predict function bug for single sample prediction fixed -
plotLoadings bug for long variable names fixed

minor improvements

- missing values in plotIndiv's group argument no more throws error -
mixOmics::predict function documentation now more accessible

Changes in version 6.8.5

bug fixes

- names of linnerud datasets fixed.

Changes in version 6.8.4

minor improvements

- package startup message with direct liks to useful resources -
mixOmics function documentation disambiguated with instruction on how
to get package help.

Changes in version 6.8.3

new features / enhancements

- You can now customise auroc plots. Refer to documentation for more
info.

bug fixes

- Fixed tune.spls and pef.plsda bugs when using cpus argument for
parallel processing

minor improvements

- auroc help files now updated with latest changes

Changes in version 6.8.2

minor improvements

- Updated onLoad message with discussion forum info, bug reports, and
more - Dropped legacy comp.tol argument from pca

Changes in version 6.8.1

bug fixes

- perf.plot bug in extracting names fixed - Few fixes for tune.splsda
with AUC

minor improvements

- plot.perf now respects ylim arguments for custom y range - Added
code of conduct - Updated DESCRIPTION with bug reports and biocViews
- Updated README

[MMAPPR2](/packages/MMAPPR2)
-------

Changes in version 0.99.0

- Submitted to Bioconductor

[MMUPHin](/packages/MMUPHin)
-------

Changes in version 0.99.2 (2019-10-21)

- Changed default output of Maaslin2_wrapper and rma_wrapper to
  tempdir()

Changes in version 0.99.1 (2019-10-15)

- Switched reference links for igraph and fpc in Rd documents.

Changes in version 0.99.0 (2019-10-15)

- Initial submission to Bioconductor.

[MOSim](/packages/MOSim)
-----

Changes in version 0.99.0 (2019-06-13)

- Submitted to Bioconductor

[motifStack](/packages/motifStack)
----------

Changes in version 1.29.8

- remove google scholar from vignette.

Changes in version 1.29.7

- Add google scholar in vignette.

Changes in version 1.29.6

- fix a bug that 'grid' and 'graphics' output mixed for
  plotMotifStackWithPhylog

Changes in version 1.29.5

- fix a bug that importMatrix can not handle empty matrix.

Changes in version 1.29.4

- fix a bug that markers shift position when alignment.

Changes in version 1.29.3

- accept pcm for plotMotifLogoA

- geom_motif accept x,y,width,height

Changes in version 1.29.2

- add markers

Changes in version 1.29.1

- add function geom_motif

[msa](/packages/msa)
---

Changes in version 1.18.0

- release as part of Bioconductor 3.10

Changes in version 1.17.2

- removed change history from package vignette for easier maintenance

Changes in version 1.17.1

- fixed regular expression to comply with PCRE2

- fixed Windows makefile for gc lib

- fixed Windows cleanup script

- fixed src/Makevars.win

Changes in version 1.17.0

- new branch for Bioconductor 3.10 devel

[MSnbase](/packages/MSnbase)
-------

Changes in version 2.11

Changes in 2.11.13

- Fix bug in readSRMData (issue #486) which caused the function to
not discriminate between chromatograms with different precursor
collision energy.

Changes in 2.11.12

- Fix bug in pickPeaks: set variable centroided only to TRUE for
spectra matching the provided msLevel..  - Add parameter msLevel. to
smooth.

Changes in 2.11.11

- Update news files titles to be consistent with other Bioc packages.

Changes in 2.11.10

- Updates to consensusSpectrum: add parameter mzFun allowing to
define the aggregation function for m/z values, change default for
mzd to 0 (hence avoid estimating m/z differences by default) and
change the default for parameter intensityFun to median.

Changes in 2.11.9

- Use latest mzR

Changes in 2.11.8

- Update URL to point to pkgdown site

Changes in 2.11.7

- Add parameter weighted to consensusSpectrum and change the default
from reporting the intensity-weighted mean of m/z values for
consensus peaks to reporting the m/z of the largest peak <2019-09-12
Thu>.

Changes in 2.11.6

- Update to match changes in mzR version 2.17.4 <2019-09-04 Wed>.  -
Add parameter msLevel to pickPeaks to allow peak picking in specific
MS levels. See #478 <2019-09-04 Wed>.

Changes in 2.11.5

- Use filterMsLevel, filterMz, filterPolarity, filterRt,
filterAcquisitionNum, filterEmptySpectra, filterPrecursorScan,
productMz, filterPrecursorMzandfilterIsolationWindowgenerics
fromProtGenerics`.

Changes in 2.11.4

- plot,Spectrum,Spectrum now also supports MS1 spectra (see #477)
<2019-07-23 Tue>

Changes in 2.11.3

- Make combineFeatures a method <2019-05-31 Fri> - Remove message
about changed meaning of the "modifications" argument in
calculateFragments' that was introduced in MSnbase 1.17.6
(2015-06-21). <2019-06-01 Sat> - Implement combineSpectra,MSnExp (see
#474) <2019-06-02 Tue>.

Changes in 2.11.2

- Fix bug in calculateFragments for neutral loss calculation. For the
"loss of water" the mass of HO~2~ instead of H~2~O was removed (see
#462). Thanks to Max Helf (@mjhelf) for the fix (see #463)
<2019-05-31 Fri>.

Changes in 2.11.1

- Migrate generics to ProtGenerics

Changes in 2.11.0

- Bioconductor 3.10 (devel)

[msPurity](/packages/msPurity)
--------

Changes in version 1.11.5

- frag4feature fileid fix for conversion from factor to character

- Add missing plyr:: reference (thanks jsaintvanne)

Changes in version 1.11.3

- Overhaul of combineAnnotation function. Uses local database now as
  previously API calls would take too much time to finish and was not
  usable

- Various updates of createMSP to make compatible with Galaxy workflows

- Parameter added to purityA to allow user to change the PPM tolerance
  for MZ values between scans when calculated the interpolated
  precursor ion purity

- Update of spectralMatching results columns to include additional
  details (e.g. retention time)

- Update of spectralMatching so that either PostgreSQL or MySQL
  database can be used as input to either query or library

Changes in version 1.11.2

- Bug fix for EIC with MSMS data

Changes in version 1.11.1

- Bug fix for duplicate MSP spectra when not using metadata table

- Added xcms3 to xcmsSet conversion for "create database" code

- Fix for sirius combine annotations (incorrect column format)

Changes in version 1.11.0

- Bioconductor dev (automatic version bump)


[mzR](/packages/mzR)
---

Changes in version 2.19.6

- header for the pwiz backend returns NA instead of 0 for not defined
  or missing information <2019-09-24 Tue>.

- peaks for pwiz backend rewritten (small performance improvement)
  <2019-09-26 Thu>.

Changes in version 2.19.5

- version bump to force build with latest Rcpp

Changes in version 2.19.4

- Add header columns scanWindowLowerLimit and scanWindowUpperLimit

Changes in version 2.19.3

- use ProtGenerics::tolerance generic <2019-08-16 Fri>

Changes in version 2.19.2

- Fix issue 190, compiles on clang-8.0

Changes in version 2.19.1

- Remove analyzer generics, now in ProtGenerics <2019-05-13 Mon>

[NADfinder](/packages/NADfinder)
---------

Changes in version 1.9.2

- Update citation.  - Update trimPeaks.R to filter peaks using
adjusted pvalues and trimmed peaks by removing windows with zscore <
1.696

Changes in version 1.9.1

- Add citation.

[NBAMSeq](/packages/NBAMSeq)
-------

Changes in version 1.1.1 (2019-08-16)

- Added a new function makeplot

- Visualization part in vignette is modified

[ncGTW](/packages/ncGTW)
-----

Changes in version 0.99.7 (2019-08-22)

- Added citaion

Changes in version 0.99.2 (2019-06-19)

- Submitted to Bioconductor

Changes in version 0.5.0 (2019-04-02)

- Uploaded to GitHub

[netboost](/packages/netboost)
--------

Changes in version 1.1.3 (2019-08-01)

- Introduction of the fully rank based extension
  (netboost(...,robust_PCs=TRUE,filter_method="spearman",method="spearman")).

Changes in version 1.1.1 (2019-05-06)

- Introduction of Pearson-, Spearman- and Kendall-based filtering.

[ngsReports](/packages/ngsReports)
----------

Changes in version 1.1.1

- Added plotAlignmentSummary()

- Added plotFastqcPCA()

- Added quast, busco, cutadapt, featureCounts, trimmomatic, flagstats &
  AdapterRemoval support to importNgsLogs()

- Enabled auto detection for report type for importNgsLogs()

Changes in version 1.0.2

- Added Transcriptomic GC Content for A.thaliana to default
  gcTheoretical object

Changes in version 1.0.1

- Table in default FastQC template now scroll for larger datasets

- Kmers removed from default FastQC template

- Typos in vignette corrected, seperate LICENSE file added & dplyr
  updates corrected

- Corrected dependencies for writeHtmlReport

[NormalyzerDE](/packages/NormalyzerDE)
------------

Changes in version 1.3.3

- Extended input validation when executing NormalyzerDE using
  experimentObj

Changes in version 1.3.2

- Correctly showing the number of values replaced with NA in status
  text

Changes in version 1.3.0

- Sync with Bioconductor changes

[normr](/packages/normr)
-----

Changes in version 1.11.2

- Added correct CITATION

Changes in version 1.11.1

- Maintainer E-mail adress updated

[omicplotR](/packages/omicplotR)
---------

Changes in version 1.5.4

- explicitly calling some functions to prevent conflicts
  (jsonlite::fromJSON)

- corrected namespace and description

- removed "require(libraries") from ui.

- bug fix for calculating aldex object. sometimes would need to click
  generate effect plot twice.

- version bump

- changed T to TRUE in rab_script and server.R

Changes in version 1.5.1

- choice of pseudocount or CZM

- changed how conditions are selected manually for effect size. now,
  you input the column number

- added ability to download GO slim annotated feature tables directly
  from MGNify database by inputting a Study ID

- added density plots for interactive effect sizes

- currently depends on specific version of ALDEx2 (temporary)

[OmnipathR](/packages/OmnipathR)
---------

Changes in version 0.99.12 (2019-10-21)

- Modification in the separation between genes within a complex (From
dash) to Underscore

Changes in version 0.99.0 (2019-10-10)

- Submitted to Bioconductor

[OncoSimulR](/packages/OncoSimulR)
----------

Changes in version 2.15.2 (2019-08-14)

- Trying to prevent fscanf warning in FitnessLandscape/input.c

Changes in version 2.15.1 (2019-06-06)

- Added MAGELLAN's sources and functionality from MAGELLAN.

Changes in version 2.15.0 (2019-06-06)

- Bumped version to match current Biocdevel.

[onlineFDR](/packages/onlineFDR)
---------

Changes in version 1.4.0

MODIFICATIONS

- added online FWER algorithms of Tian and Ramdas &#91;2019b&#93;

- added the ADDIS algorithms of Tian and Ramdas &#91;2019a&#93;

- added asynchronous online testing algorithms of Zrnic et al. &#91;2018&#93;

- added the SAFFRON procedure for online FDR control &#91;Ramdas et al.,
  2018&#93;

- added the Alpha-investing procedure of Ramdas et al. &#91;2018&#93;

- updated vignette

- deprecated LORDdep and added functionality to LORD

- deprecated Bonfinifinite, which is replaced by AlphaSpending

- removed LORD versions 1 and 2

- added unit tests

- updated references

- updated authors

[oppti](/packages/oppti)
-----

Changes in version 0.99.13

Features in the first version, Bioconductor 3.10 Release (September 2019

- artImpute Artificially miss and impute each data entry individually
  by ignoring outlying values

- clusterData Hierarchical cluster analysis

- dropMarkers Filter out markers

- dysReg Analyze dysregulated (protruding) events

- markOut Display outlying expressions

- oppti Outlier protein and phosphosite target identification

- outScores Analyze putative outliers

- plotDen Draw densities

- rankPerOut Rank markers by the percentage of outlying events

- statTest Analyze dysregulation significance

[Organism.dplyr](/packages/Organism.dplyr)
--------------

Changes in version 1.14.0

NEW FEATURES

- src_organism() supports an option overwrite=FALSE to optionally
  over-write exisiting (cached) resources created from a previous txdb
  version.

- src_organism() supports construction from a TxDb object.

[OUTRIDER](/packages/OUTRIDER)
--------

Changes in version 1.3.2

- Documentation

- New plot functionality - plotExpectedVsObservedCounts() -
  plotCountGeneSampleHeatmap() - plotExpressedGenes()

- Minor bug fixes

- Linking to publication REVISED VERSION 0.99.29

- Major changes

- Improved autoencoder model.

- Updated API

- Improved default parameters. SECOND BIOCONDUCTOR SUBMISSION 0.99.10

- Smaller and faster example data set

- Loss and gradient of loss in C++

- Bugfixes INITIAL BIOCONDUCTOR SUBMISSION 0.99.8

- Better documentation

- Little changes in default values

- Code cleanup and improvements

- Bugfixes R IMPLEMENTATION OF THE AUTOENCODER 0.99.7

- Adding the R implementation of the autoencoder

- Little changes in default values

- Updating and completing documentation

- Bugfixes MINOR BUGFIXES AND UPDATE OF AUTOCORRECTION 0.99.6

- Minor changes in plotting functions

- Update interface to autoCorrection

- Bugfixes GITHUB RELEASE OF OUTRIDER 0.99.5

- Pre-release of OUTRIDER on GitHub INITIAL SETUP OF OUTRIDER VERSION
  0.99.2

- Initial setup of OUTRIDER package

- Added autoCorrect (Auto Encoder) as normalization function

[pathwayPCA](/packages/pathwayPCA)
----------

Changes in version 1.1.1

2019-06-06

Our build on Bioconductor 3.9 devel fails for the second vignette.
This patch resolves this issue.

[Pbase](/packages/Pbase)
-----

Changes in version 0.25.1

- Removing mapping functions, moved to the ensembldb package.
  <2019-08-08 Thu>

- Marking package for deprecation <2019-08-08 Thu>

[pcaExplorer](/packages/pcaExplorer)
-----------

Changes in version 2.12.0

Bug fixes

- Fixed an error in the initialization of the app due to a new
behavior introduced by shinyAce in version >= 0.4.0 - topGOtable does
not generate rows with NAs if providing a too high number for the
categories to report

Other notes

- The type of the columns in the data.frame returned by topGOtable
are now correctly referring to the type they contain - e.g. the p
values are now stored as numeric values - Citation now refers to the
published manuscript - https://doi.org/10.1186/s12859-019-2879-1

[PCAtools](/packages/PCAtools)
--------

Changes in version 2.0.0

- added parallelPCA function to perform Horn's parallel analysis, which
  chooses an ideal number of principal components to retain (courtesy
  Aaron Lun)

- added findElbowPoint function, which finds the elbow point in the
  curve of variance explained and which can also be used to determine
  the number of principal components to retain (courtesy Aaron Lun)

- user can now specify custom labels for points

- fixed bug with singlecol parameter for biplot colouring everything
  black

[peakPantheR](/packages/peakPantheR)
-----------

Changes in version 0.99.3 (2019-10-01)

- Revisions for Bioconductor submission

Changes in version 0.99.2 (2019-09-10)

- Revisions for Bioconductor submission

Changes in version 0.99.0 (2019-06-09)

- Submitted to Bioconductor

[PERFect](/packages/PERFect)
-------

Changes in version 0.99.11 (2019-09-19)

- Fixed bugs

- Made the following significant changes o Removed Knight dataset

Changes in version 0.99.10 (2019-09-17)

- Submitted to Bioconductor

[PGA](/packages/PGA)
---

Changes in version 1.15.1

- Implementation for the protein database construction from fusion
  events(buildFusionProteinDB). Currently, this function only supports
  the calling result from STAR-fusion.

Changes in version 1.14.1

- Add a parameter to set FDR

- Update parser program V

[PhyloProfile](/packages/PhyloProfile)
------------

Changes in version 0.99.31

- Submitted to Bioconductor

[piano](/packages/piano)
-----

Changes in version 2.2.0

BUG FIXES

- Fix import exprs<- bug in loadMAdata.

Changes in version 2.0.1

BUG FIXES

- Fix error message handling in network plot in exploreGSAres().

[Pigengene](/packages/Pigengene)
---------

Changes in version 1.11.34 (2019-10-21)

General

- The pipeline is now explained step by step in the vignette.

Changes in version 1.11.32 (2019-10-18)

Changes in existing functions

- The Data argument of compute.pigengene can now be a matrix with only
  1 column.

Changes in version 1.11.30 (2019-10-02)

Changes in existing functions

- The doRetuNetworks argument is now added to the combine.networks
  function.

Changes in version 1.11.28 (2019-10-01)

New functions

- message.if() is now exported.

Changes in version 1.11.26 (2019-09-26)

Changes in existing functions

- Better QC in the gene.mapping() function, the possible keys will be
  printed if the input is not appropriate.

New functions

- save.if() is now exported.

Changes in version 1.11.24 (2019-09-03)

Bug Fixes

- In the combine.network() function, selectedModules does not need to
  be "All". Also, if saveFile=NULL, nothing will be saved without any
  error.

Changes in version 1.11.20 (2019-05-15)

Changes in existing functions

- Data and Labels can now be lists, which will be combined using
  combine.network() before analysis.

Changes in version 1.11.4 (2019-05-02)

Bug Fixes

- repeat.data(times=1,...) now produces valid output.

[podkat](/packages/podkat)
------

Changes in version 1.18.0

- release as part of Bioconductor 3.10

Changes in version 1.17.3

- removed change history from package vignette for easier maintenance

Changes in version 1.17.2

- minor changes to DESCRIPTION file (system requirement GNU make) and
  src/Makevars

Changes in version 1.17.1

- changed summary() method for VariantAnnotation class in order to stay
  compatible with print() method in GenomicRanges package

- corresponding minor adaptations in documentation and package vignette

Changes in version 1.17.0

- new branch for Bioconductor 3.10 devel

[polyester](/packages/polyester)
---------

Changes in version 1.99.3

- NB function now exported

- note that version 1.99.3 on GitHub was version 1.1.0 on Bioconductor.

Changes in version 1.99.2

- bug fix in fragment generation (last 2 bases of transcript were never
  sequenced)

[pqsfinder](/packages/pqsfinder)
---------

Changes in version 2.2

SIGNIFICANT USER-VISIBLE CHANGES

- Default minimal PQS score was decreased from 52 to 47. The score 47
  shows the best balanced accuracy on new G4 sequencing data provided
  by Marsico et al. 2019.

BUG FIXES

- Fixed bug allowing unlimited length of third loop leading to invalid
  memory access and random scores.

[proDA](/packages/proDA)
-----

Changes in version 0.99.0 (2019-06-02)

- Initial submission

[pRoloc](/packages/pRoloc)
------

Changes in version 1.25

Changes in version 1.25.2

- Fix new biomart attribute <2019-08-09 Fri> - Bug fix: pass fcol to
helper function (see https://support.bioconductor.org/p/123614/)
<2019-08-09 Fri>

Changes in version 1.25.1

- Always use mvtnorm::dmvnorm <2019-06-21 Fri>

Changes in version 1.25.0

- Version bump for Bioc 3.10 (devel)

[ProtGenerics](/packages/ProtGenerics)
------------

Changes in version 1.17.4

- New generics needed in Spectra and Chromatograms packages <2019-08-20
  Tue>

Changes in version 1.17.3

- New tolerance generic <2019-08-16 Fri>

Changes in version 1.17.2

- New generics (see issue #8) <2019-05-12 Sun>

Changes in version 1.17.1

- Move many generics from MSnbase <2019-05-11 Sat>

[PureCN](/packages/PureCN)
------

Changes in version 1.16.0

NEW FEATURES

- Flag segments in poor quality regions

- predictSomatic now provides log-likelihood of allelic balance
  (ALLELIC.IMBALANCE column) for each variant

- Added readLogRatioFile function to read GATK4 DenoiseReadCounts
  output files containing log2 tumor/normal ratios

- Added readSegmentationFile function to read GATK4 ModelSegment output
  files containing segmented log2 tumor/normal ratios

- Added callAmplificationsInLowPurity to call gene-level amplifications
  in samples < 10% purity

- Dx.R now reports chromosomal instability scores (available also via
  callCIN function)

- Dx.R supports deconstructSigs 1.9.0 and COSMIC signatures v3. To run
  both v2 and v3, simply add --signature_databases
  signatures.exome.cosmic.v3.may2019:signatures.cosmic to Dx.R

SIGNIFICANT USER-VISIBLE CHANGES

- Made filterTargets and createTargetWeights defunct

- setMappingBiasVcf now returns a data.frame

- Best practices vignette now HTML-based

- Renamed normal.panel.vcf.file in setMappingBiasVcf to
  mapping.bias.file; in 1.18, setMappingBiasVcf will not accept a VCF
  anymore but requires a precomputed mapping bias RDS file.

- calculateIntervalWeights now directly called by createNormalDatabase
  and information included in the normalDB RDS object. This function is
  thus deprecated.

- Column gene.mean in callAlterations output now weighted by interval
  weights when available

- Changed default of min.target.width in preprocessIntervals from 10 to
  100 (#73)

- replaced write.table with data.table::fwrite to automatically support
  producing gzipped output (requires data.table 1.12.4, #106)

- Coverage.R now gzips BAM file coverage (requires data.table 1.12.4,
  #106)

- Output coverage files now code FALSE as 0 and TRUE as 1

- PureCN.R now bgzips and tabix indexes VCFs when --vcf is provided

BUGFIXES

- Fix for bug in CCF calculation resulting in NAs (happens in high
  coverage samples, early mutations with > 1 allele copy number)

- Fix for a bug in preprocessIntervals when small targets (<
  min.target.width) were present

- Fix for a bug in callMutationBurden when VCF contained indels (#82)

- Die with helpful error message when snp.blacklist import failed

- Check input segmentation files for missing values resulting in crash

- Fixed a crash in Varscan2 produced VCFs when ALT field missed ref
  counts (#109)

[pwrEWAS](/packages/pwrEWAS)
-------

Changes in version 0.1.0

- Initial version

[QDNAseq](/packages/QDNAseq)
-------

Changes in version 1.21.6 (2019-09-25)

BUG FIXES

- Link to the 'GEM library' tool was broken.

Changes in version 1.21.5 (2019-09-09)

BUG FIXES

- plot(..., logTransform=TRUE, doSegments=TRUE) on QDNAseqSignals would
  position segments that were out of range incorrectly, because it
  forgot to take the log transform on those outliers.

Changes in version 1.21.4 (2019-09-06)

SIGNIFICANT CHANGES

- Bin annotation data files are no longer downloaded automatically.
  They are instead part of Bioconductor annotation packages
  QDNAseq.hg19 and QDNAseq.mm10, which needs to be installed by the
  user.  If not installed, an error is now produced.  The reason for
  this change, is that the QDNAseq maintainers will no longer host
  QDNAseq bin annotation files online (in the cloud).

Changes in version 1.21.3 (2019-09-04)

NEW FEATURES

- All functions that produce verbose output gained argument 'verbose'.
  For backward compatibility, it currently defaults to 'verbose=TRUE'
  but that may be changed to 'verbose=FALSE' in a future release.

IMPROVEMENTS

- callBins() now respects option 'QDNAseq::verbose' for controlling
  whether output from the CGHcall package should be relayed or not.

- MEMORY: Utilize more memory-efficient matrixStats functions
  colSums2(), colMeans2(), etc.

DEPRECATION AND DEFUNCT

- Argument 'seeds' of segmentBins() is deprecated because it did not
  use proper parallel random number generation (RNG).  We now instead
  rely on future.apply::future_lapply(..., future.seed=TRUE) for this.

Changes in version 1.21.2 (2019-09-03)

SIGNIFICANT CHANGES

- Package now imports the 'future' and 'future.apply' packages;
  previously the 'future' was listed as a suggested package.

- Package no longer depends on BiocParallel.

- binReadCounts() now uses the future framework instead of BiocParallel
  for parallelization.

IMPROVEMENTS

- MEMORY: Avoiding data type coercions in more places by for instance
  making sure that vectors and matrices are initated with the values of
  the correct data type (instead of the default NA, which is a logical
  value).

SOFTWARE QUALITY

- Using future_lapply() and future_apply() of the well-tested
  future.apply package instead of internal analogue implementations.

- TESTS: Now testing numerical reproducibility also for parallel
  processing (using future strategies 'multisession' and 'multicore').

- TESTS: Now asserting numerical reproducibility of also segmentBins()
  and callBins().

Changes in version 1.21.1 (2019-08-30)

SIGNIFICANT CHANGES

- exportVCF() is no longer exported. Use exportBins(..., format="vcf")
  instead.

[qPLEXanalyzer](/packages/qPLEXanalyzer)
-------------

Changes in version 1.3.2

- Added mergePeptides to NAMESPACE

Changes in version 1.3.1

- Added mergePeptides function to summarize identical peptides to
  single peptide for specific protein.

[Qtlizer](/packages/Qtlizer)
-------

Changes in version 0.99.12 (2019-09-06)

- eval=false for package installation in vignette

- Title change in vignette

Changes in version 0.99.11 (2019-09-03)

- Incorporated Bioconductor reviewer comments

Changes in version 0.99.10 (2019-08-27)

- Minor

Changes in version 0.99.9 (2019-08-27)

- Bugfixes, error handling

Changes in version 0.99.8 (2019-08-27)

- Additional HTTP error handling

Changes in version 0.99.7 (2019-08-27)

- Minor

Changes in version 0.99.6 (2019-08-27)

- Code revised

Changes in version 0.99.5 (2019-08-27)

- Added parameter return_obj to have the possibility to get a
  GenomicRanges::GRanges object

Changes in version 0.99.4 (2019-08-14)

- Added parameter ld_method and corr

Changes in version 0.99.0 (2019-08-12)

- Submitted to Bioconductor

[RaggedExperiment](/packages/RaggedExperiment)
----------------

Changes in version 1.10.0

Bug fixes and minor improvements

- Include reference to TCGAutils functions for qreduceAssay examples
- Add robustness to RaggedExperiment constructor including unit tests
- Include class and assay operations overview schematic in the
vignette

[Rcpi](/packages/Rcpi)
----

Changes in version 1.21.1 (2019-05-17)

Improvements

- Removed AppVeyor CI due to the frequent Bioconductor installation
and dependency issues which are not related to the package itself.  -
Updated GitHub repository links due to the recent handle change.  -
Updated the vignette style.

[Rcwl](/packages/Rcwl)
----

Changes in version 1.1.12 (2019-09-09)

- Add ExpressionTool, Extensions and Metadata

Changes in version 1.1.7 (2019-07-25)

- New feature, baseCommand works with R function

[RcwlPipelines](/packages/RcwlPipelines)
-------------

Changes in version 1.1.14 (2019-10-21)

- Add pvactools neoantigen prediction pipeline

Changes in version 1.1.8 (2019-08-22)

- Add more Somatic Variant Callers: VarScan2, LoFreq, MuSE,
  SomaticSniper, VarDict, lancet, manta, strelka2

Changes in version 1.1.5 (2019-07-19)

- Structure updates: export all tools and pipelines.

- Added GATK4 Mutect2 tools and pipelines.

- Added Command and Container to metadata of cwlTools.

- Vignette updated.

Changes in version 1.1.2 (2019-05-22)

- fix samtools index bug for cwltool 2019

[RCy3](/packages/RCy3)
----

Changes in version 2.6.0

- New functions: - createGroupByColumn - clearEdgeBends -
  getNodePosition

- New parameter to return SUIDs for - getSelectedNodes -
  getSelectedEdges

- Node and edge property values returned as named lists

- Faster results for getting all node and edge property values, #78

- More robust handling of file type in export functions

- More robust handling of dataframes in createNetworkFromDataFrames

- New support for loading list data

- Doc Fixes - added Filters to Overview vignette - improved file type
  handling descriptions

Changes in version 2.4.4

- Bug Fixes - filter functions -- #73 wrong params

Changes in version 2.4.3

- Bug Fixes - import functions -- #62 fixed default directory

Changes in version 2.4.2

- Bug Fixes - getEdgeInfo -- #61 missing function

Changes in version 2.4.1

- Bug Fixes - getLayoutPropertyNames -- #59 fixed returned values -
  createNetworkFromIgraph -- #58 flatten list attributes to strings

[ReactomeGSA](/packages/ReactomeGSA)
-----------

Changes in version 0.99.8 (2019-08-13)

- Removed package startup message

- Changed spelling of name to ReactomeGSA

- Removed all default values from the method signatures

- Added constructor function for the ReactomeAnalysisRequest class.

Changes in version 0.99.7 (2019-07-22)

- Adapted the get_reactome_methods function to provide a more readable
  overview

- Changed to reactome analysis service API URL to the new domain name.

Changes in version 0.99.6 (2019-07-09)

- Minor bugfixes in examples

Changes in version 0.99.5 (2019-07-09)

- Updated vignette to new API data types

- Added new function remove_dataset

- Fixed bugs when overwriting existing datasets in AnalysisRequests

- Fixed bug that pathways function did not sort the pathway table

Changes in version 0.99.0 (2019-07-01)

- Submitted to Bioconductor

[recount](/packages/recount)
-------

Changes in version 1.11.14

- Added a NEWS.md file to track changes to the package.

Changes in version 1.11.13

NEW FEATURES

- Added the function getTPM() as discussed in
https://support.bioconductor.org/p/124265 and based on Sonali Arora
et al https://www.biorxiv.org/content/10.1101/445601v2.

Changes in version 1.11.12

BUG FIXES

- Now geo_characteristics() can deal with the scenario reported at
https://support.bioconductor.org/p/116480/ by @Jacques.van-Helden.

Changes in version 1.11.7

SIGNIFICANT USER-VISIBLE CHANGES

- Renamed .load_install() as .load_check() as this function now only
checks that the package(s) was installed and returns an error if
missing. The error shows the user how to install the package(s) they
are missing instead of installing them automatically. This complies
with Marcel Ramos' request at
https://github.com/leekgroup/recount/issues/14.

Changes in version 1.11.4

NEW FEATURES

- Added the function download_retry() based on
http://bioconductor.org/developers/how-to/web-query/ such that
download_file() and other recount functions will re-try to download a
file 3 times before giving up. This should help reduce the number of
occasional failed Bioconductor nightly checks.

[regioneR](/packages/regioneR)
--------

Changes in version 1.18.0

NEW FEATURES

- Expanded toGRanges support. It is now possible to transform coverage
  objects (i.e. toGRanges(coverage(A))) into GRanges. It also supports
  ".assoc" files produced by PLINK.

- overlapPermTest now supports multiple region sets in B and will
  perform a multi- permutation test against each one much faster than
  testing them independently.

BUG FIXES

- Multiple bug fixes

[regionReport](/packages/regionReport)
------------

Changes in version 1.19.2

- Added a NEWS.md file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

- Renamed load_install() as load_check() as this function now only
checks that the package(s) was installed and returns an error if
missing. The error shows the user how to install the package(s) they
are missing instead of installing them automatically. This complies
with Marcel Ramos' request at
https://github.com/leekgroup/recount/issues/14.

[rGREAT](/packages/rGREAT)
------

Changes in version 1.17.1

- assign to foreName and backName when querying GREAT

- use rmarkdown for vignette

- support GREAT version 4.0.4

- add startup messages

[rhdf5](/packages/rhdf5)
-----

Changes in version 2.30.0

NEW FEATURES

- Functions H5Lmove & H5Lcopy are now exported and accessible.

BUG FIXES

- Source file names are no longer mangled when printing error messages.

- NA values in a character() vector can now be written to an HDF5
  dataset.

[Rhdf5lib](/packages/Rhdf5lib)
--------

Changes in version 1.8

New features

- Updated internal version of HDF5 to 1.10.5

Bug fixes

- Quote paths reported by pkgconfig() to allow installation in librarys
  with whitespace in the path.

[Rhtslib](/packages/Rhtslib)
-------

Changes in version 1.18.0

NEW FEATURES

- Compile HTSlib with libcurl enabled

- Support an installation path that contains whitespaces

SIGNIFICANT USER-VISIBLE CHANGES

- Switch from dynamic to static linking on all Unix-like systems (see
  commit db1d8e17)

- Package now requires libbz2 & liblzma & libcurl (with header files),
  and GNU make. This is declared in new SystemRequirements field.

BUG FIXES

- Use preprocessor flag -D_FILE_OFFSET_BITS=64. This addresses nasty
  problem with big files that get truncated on Windows. See
  https://support.bioconductor.org/p/124568/

- Don't overwrite CPPFLAGS, CFLAGS, or LDFLAGS values set in
  ${R_HOME}/etc/Makeconf on Linux or Mac

[RNAmodR](/packages/RNAmodR)
-------

Changes in version 0.99.0 (2019-04-29)

- Submitted to Bioconductor

[RNAmodR.AlkAnilineSeq](/packages/RNAmodR.AlkAnilineSeq)
---------------------

Changes in version 0.99.0 (2019-04-29)

- Submitted to Bioconductor

[RNAmodR.ML](/packages/RNAmodR.ML)
----------

Changes in version 0.99.0 (2019-04-29)

- Submitted to Bioconductor

[RNAmodR.RiboMethSeq](/packages/RNAmodR.RiboMethSeq)
-------------------

Changes in version 0.99.0 (2019-04-29)

- Submitted to Bioconductor

[RnBeads](/packages/RnBeads)
-------

Changes in version 2.3.3

- Fixed non ASCII column names in age prediction and immune estimation

Changes in version 2.3.2

- Added stratification plots for inferred covariates in different
  sample groups

[rols](/packages/rols)
----

Changes in version 2.13

CHANGES IN VERSION 2.13.2

- Remove failing test - term is now obsolete <2019-10-12 sam.>

CHANGES IN VERSION 2.13.1

- Fix failing test <2019-08-12 Mon>

[ropls](/packages/ropls)
-----

Changes in version 1.17.34

NEW FEATURE

- 'view' method (wrapper of the imageF and strF functions)

Changes in version 1.17.32

BUG FIXED

- plot.oplsMultiDataSet correction to include all plots in the .pdf
  file

Changes in version 1.17.30

MINOR MODIFICATION

- Additional correction in getMset documentation (example)

Changes in version 1.17.28

MINOR MODIFICATION

- Correcting typo in getMset documentation (example)

Changes in version 1.17.26

MINOR MODIFICATION

- Including the set name in the figure of MultiDataSet models

Changes in version 1.17.24

MINOR MODIFICATION

- Correction of the documentation for 'residuals' and
  'plot-oplsMultidDataSet'

Changes in version 1.17.22

INTERNAL MODIFICATION

- Correction of the example in the oplsMultiDataSet class documentation

Changes in version 1.17.20

NEW FEATURE

- plotPhenoDataC parameter for coloring (score plots) according to one
  column from the pData data frame when the opls method has been
  applied to an ExpresssionSet

Changes in version 1.17.18

INTERNAL MODIFICATION

- Minor internal modification

Changes in version 1.17.16

INTERNAL MODIFICATION

- Minor internal modification

Changes in version 1.17.14

NEW FEATURE

- Example of the analyis of a MultiDataSet ('NCI60_4arrays' from the
  'omicade4' package)

Changes in version 1.17.12

NEW FEATURE

- Multi data set names display on graphics

Changes in version 1.17.10

NEW FEATURE

- plot method for oplsMultiDataSet objects

Changes in version 1.17.8

NEW FEATURE

- opls can now be applied to MultiDataSet objects; getMset method to
  extract the complemented MultiDataSet

Changes in version 1.17.6

INTERNAL MODIFICATION

- minor internal modification

Changes in version 1.17.4

INTERNAL MODIFICATION

- minor internal modification

Changes in version 1.17.2

NEW FEATURE

- info.txtC and fig.pdfC argument values NULL and NA replaced by 'none'
  and 'interactive', respectively

[rpx](/packages/rpx)
---

Changes in version 1.21

Changes in version 1.21.3

- Update NEWS titles

Changes in version 1.21.2

- Temporary fix for issue #5 (wrong URL from PRIDE server)
<2019-10-02 Wed>

Changes in version 1.21.1

- Don't set old class <2019-08-09 Fri>

[Rsubread](/packages/Rsubread)
--------

Changes in version 2.0.0

- Rsubread package is ported to Windows OS.

- New function cellCounts(): generate UMI counts for Chromium 10X
  single-cell RNA-seq data.

- flattenGTF() function can merge or chop overlap features.

- Check and display the amount of memory available on the computer
  before starting read mapping.

- Optimize the data structure used in buildindex() function to reduce
  its memory use.

- qualityScores() function can optionally retrieve quality scores from
  all the reads.

- File paths included in column names of objects returned by
  featureCounts(), align(), subjunc() and propmapped() functions are
  removed or shortened where appropriate.

- featureCounts() will be terminated if both single-end and paired-end
  reads are found in the same input file.

- Limit on the length of input file names is increased to 1000 bytes
  for all functions.

[RTCA](/packages/RTCA)
----

Changes in version 2009-07-13

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

[RTCGAToolbox](/packages/RTCGAToolbox)
------------

Changes in version 2.16.0

New features

- RNASeq2GeneNorm slot in the FirehoseData class is a list now (from
matrix) - Use tempdir() as the default directory for downloading data
in getFirehoseData

Bug fixes and minor improvements

- Save all RNASeq2GeneNorm datasets within the output object as a
list. Previously, only the last dataset would get returned (#30) -
Read files from the appropriate download location in getFirehoseData
- Move static text file references from 'canevolve.org' to GitHub
hosted locations - Check file sizes using httr instead of
'canevolve.org' query (@mksamur, #32)

[rWikiPathways](/packages/rWikiPathways)
-------------

Changes in version 1.6.0

- Minor updates to the Pathway Analysis vignette

Changes in version 1.4.1

- Bug fixes - Updated the findPathwaysByLiterature test case (fixing a
  false positive fail) - Removed tests for unavailable service:
  getColoredPathway

[S4Vectors](/packages/S4Vectors)
---------

Changes in version 0.24.0

NEW FEATURES

- Add Factor class. Serves a similar role as factor in base R except
  that the levels of a Factor object can be any Vector derivative.

- New methods for DataFrame comparisons (by Aaron Lun)

- Add sameAsPreviousROW() generic and methods for ANY, atomic, integer,
  numeric, complex, Rle, DataFrame, and Pairs (by Aaron Lun)

- Support more comparison methods for Pairs objects

- Add methods for coercing back and forth between HitsList and
  SortedByQueryHitsList.

- Add anyDuplicated() method for Vector derivatives.

- Support 'by=' argument on sort,List

- Add is.finite() method for Rle objects

- Add add "&" method for FilterRules objects as a convenience for
  concatenation

SIGNIFICANT USER-VISIBLE CHANGES

- Add DFrame class (commit 36837bdf). DataFrame() now returns a DFrame
  instance (commit 83b09b19).

- Now 'stringsAsFactors' is set to FALSE when coercing something to a
  DataFrame.

- Move splitAsList() from the IRanges package

- Move S4 class "atomic" from the IRanges package

- Improve handling of user-supplied metadata columns

DEPRECATED AND DEFUNCT

- Remove phead(), ptail(), and strsplitAsListOfIntegerVectors(). These
  functions were deprecated in BioC 3.7 and defunct in BioC 3.8.

BUG FIXES

- Fix split() on a SortedByQueryHits object (issue #39)

- Fix the following coercions: - Hits -> SelfHits - SortedByQueryHits
  -> SortedByQuerySelfHits - SelfHits -> SortedByQuerySelfHits - Hits
  -> SortedByQuerySelfHits Before this fix all these coercions
  **seemed** to work but they were in fact silently producing invalid
  objects.

- A fix to anyDuplicated() method for Rle objects (commit 63495d6)

- A fix related to replacing DataFrame columns with matrix columns
  (commit 00169dd6)

- All show() methods now return an invisible NULL (commit f4b4ee76)

[SAIGEgds](/packages/SAIGEgds)
--------

Changes in version 1.0.0

- first Bioconductor release

Changes in version 0.99.0

- package submitted to Bioconductor

Changes in version 0.9.10

- SAIGE algorithm implementation for quantitative outcomes

Changes in version 0.9.9

- add a vignette to the package

- the default of random number generator changes in R: "Rounding" was
  the default in `RNGkind()` prior to R_3.6.0, but "Rejection" is used
  in R (>= v3.6.0). For reproduction of the results created by R (<
  v3.6.0), please use `RNGkind("Mersenne-Twister", "Inversion",
  "Rounding")` in R (>= v3.6.0)

Changes in version 0.9.7

- seqAssocGLMM_SPA(): load balancing in parallel

Changes in version 0.9.0

- first version of SAIGEgds

[scater](/packages/scater)
------

Changes in version 1.14.0

- Removed deprecated dplyr verbs.

- Removed deprecated method= option in runPCA().  Increased
  ncomponents= default to 50.  Deprecated use_coldata= and related
  options in favour of runColDataPCA().  Switched BSPARAM= default to
  bsparam().

- Added runColDataPCA() function for running PCA on colData(). Switch
  outlier detection strategy to avoid mvoutlier's dependency tree.

- Added the annotateBMFeatures() function to perform annotation without
  modifying the input.

- Pass all ... options to biomaRt::useMart() in getBMFeatureAnnos().

- Added name= arguments to runPCA(), etc. to change the name of the
  output reducedDim.

- Added the logNormCounts() function to compute log-normalized counts
  in an alternative experiment-aware manner.  Added a
  normalization-by-downsampling option via DropletUtils.

- Added the perCellQCMetrics() function to compute per-cell QC metrics
  in an alternative experiment-aware manner.

- Deprecated the normalize() method, which was considered too vague to
  describe what the function was actually doing.

- Added the perFeatureQCMetrics() function to compute per-feature QC
  metrics.

- Deprecated the calculateQCMetrics() function, to be replaced by the
  streamlined addQCPerCell() and addQCPerFeature().

- Generalized all functions, where possible, to operate on
  SummarizedExperiment and numeric matrices.  This involved converting
  a number of them to S4 methods to take advantage of dispatch.
  Affected functions include normalizeCounts(), calculateCPM(),
  librarySizeFactors() and so on.

- Added calculateTSNE() and related methods to operate directly on an
  input matrix.

- Renamed the use_dimred= argument to dimred=, along with similar
  renamings of other arguments for consistency.

- Report all percentages of variances explained as actual variances in
  runPCA() and getVarianceExplained().

- Added aggregateAcrossCells() and aggregateAcrossFeatures() to create
  a summed SingleCellExperiment object.

- Added the mockSCE() function to generate example objects for the
  documentation.

- Support multiple factors for grouping cells in
  sumCountsAcrossCells().

- Support list of grouping vectors in sumCountsAcrossFeatures().

- Added the order_columns_by= argument to plotHeatmap() for easy
  plotting by a given factor.  Changed defaults to more common values.

- Added a plotDots() function to create a Seurat-style dot plot.

- Dropped default nmads= to 3 in isOutlier().

[scds](/packages/scds)
----

Changes in version 1.2.0

- cxds performance improvement

- Added heuristics to estimate number of doublets

[schex](/packages/schex)
-----

Changes in version 0.99.7

- Changed sapply to vapply in utility-function.

Changes in version 0.99.6

- Changed the link from Seurat to Seurat-class.

Changes in version 0.99.5

- Fixed error introduced through merge conflict.

Changes in version 0.99.4

- Fixed typos in vignettes.  - Added explanation of ggplot connection
to vignettes.  - Added explanation for interactivity in iSEE to
vignette.  - Fixed global variables.  - Fixed documentation in
plot_hexbin_gene.  - Reformatted news.  - Added packages to suggests
in description.

Changes in version 0.99.3

- Minimal functioning package created.  - Mostly compliant with
BiocCheck::BiocCheck() and goodpractice::goodpractice().

[scMerge](/packages/scMerge)
-------

Changes in version 1.1.6

- Accepts DelayedArray, HDF5Array and dgCMatrix inputs in the slots
of input SCE objects.  - Significant speed optimisation on scSEGIndex
and add BiocParallel support.  - Updated scSEGIndex references after
publication.  - scMerge now has the svd_k input that controls the
number eigenvectors needed in the RUV step to allow fast
approximation for large dataset.  - Now using BiocSingular to manage
all SVD components.  - Now automatically remove zeroes in the rows
and columns of the SCE.

Changes in version 1.1.5

- Adding version restrictions on S4Vectors and SingleCellExperiment
dependent packages.

Changes in version 1.1.4

- plot_igraph would allow suppression of igraph output during
unsupervised scMerge

Changes in version 1.1.3

- Column name must be non-NULL and without duplicates

Changes in version 1.1.2

- Resolved problems with only a single linking cell-type across
multiple batches

Changes in version 1.1.0

- Accepted by Bioconductor

[scPCA](/packages/scPCA)
-----

Changes in version 0.99.0 (2019-09-13)

- Submitted to Bioconductor

[scran](/packages/scran)
-----

Changes in version 1.14.0

- Removed deprecated approximate= and pc.approx= arguments.

- Removed deprecated batch correction functions.

- Added option to pairwiseTTests() for standardization of log-fold
  changes.

- Changed default BSPARAM= to bsparam() in quickCluster(),
  denoisePCA(), doubletCells() and build*NNGraph().

- Added the pairwiseBinom() function for pairwise binomial tests of
  gene expression.

- Renamed output fields of pairwiseWilcox() to use AUC for less
  confusion. Added the lfc= argument to test against a log-fold change.

- Added the fitTrendVar(), fitTrendCV2(), modelGeneVar(),
  modelGeneVarWithSpikes(), modelGeneCV2(), modelCV2WithSpikes(),
  fitTrendPoisson() and modelGeneVarByPoisson() functions to model
  variability.

- Deprecated the trendVar(), technicalCV2(), improvedCV2(),
  decomposeVar(), trendVar(), testVar(), makeTechTrend(),
  multiBlockVar() and multiBlockNorm() functions.

- Modified combineVar() to not weight by residual d.f. unless
  specifically instructed.

- Added the combineCV2() function to combine separate CV2 modelling
  results.

- Added the test.type= argument in findMarkers() to switch between
  pairwise DE tests. Added the row.data= argument to easily include row
  metadata in reordered tables. Deprecated overlapExprs(), which is
  replaced by type="wilcox" in findMarkers().

- Added the getTopMarkers() function to easily retrieve marker lists
  from pairwise DE results.

- Added the getTopHVGs() function to easily retrieve HVG sets from
  variance modelling results.

- In all functions that accept a block= argument, any level of the
  blocking factor that cannot yield a result (e.g., due to insufficient
  degrees of freedom) will now be completely ignored and not contribute
  to any statistic.

- Added the getDenoisedPCs() function for general-purpose PCA-based
  denoising on non-SingleCellExperiment inputs. Converted denoisePCA()
  to a normal function, removed the method for ANY matrix. Dropped
  max.rank= default to 50 for greater speed in most cases.

- Added the calculateSumFactors() function for general-purpose
  calculation of deconvolution factors on non-SingleCellExperiment
  inputs. Converted computeSumFactors() to a normal function, removed
  the method for ANY input. Auto-guess min.mean= based on the average
  library size.

- Deprecated all special handling of spike-in rows, which are no longer
  necessary when spike-ins are stored as alternative experiments.

- Deprecated general.use= in computeSpikeFactors(), which is no longer
  necessary when spike-ins are stored as alternative experiments.

- Deprecated parallelPCA(), which has been moved to the PCAtools
  package.

- Modified clusterModularity() to return upper-triangular matrices,
  fixing a bug where the off-diagonal weights were split into two
  entries across the diagonal. Added the as.ratio= argument to return a
  matrix of log-ratios. Renamed the get.values= argument to
  get.weights=.

- Simplified density calculation in doubletCells() for greater
  robustness.

- Added a method="holm-middle" option to combinePValues(), to test if
  most individual nulls are true. Added a min.prop= option to control
  the definition of "most".

- Added a pval.type="some" option to combineMarkers(), as a compromise
  between the two other modes. Added a min.prop= option to tune
  stringency for pval.type="some" and "any".

- Added the getClusteredPCs() function to provide a cluster-based
  heuristic for choosing the number of PCs.

- Added the neighborsTo*NNGraph() functions to generate (shared)
  nearest neighbor graphs from pre-computed NN results.

- Switched to using only the top 10% of HVGs for the internal PCA in
  quickCluster().

[scTensor](/packages/scTensor)
--------

Changes in version 1.2.0

- goenrich, meshenrich, reactomeenrich, doenrich, ncgenrich, and
  dgnenrich in cellCellReport are added

- A bug related in sparse matrix in cellCellSetting is fixed

- All the vignettes are updated

- A vignette for reanalysis of the results of scTensor is added

- Some bugs are fixed

[scTGIF](/packages/scTGIF)
------

Changes in version 0.99.0

- Package released

[SeqArray](/packages/SeqArray)
--------

Changes in version 1.26.0

NEW FEATURES

- new function `seqAddValue()`

UTILITIES

- RLE chromosome coding in `seqBED2GDS()`

- change the file name "vignettes/R_Integration.Rmd" to
  "vignettes/SeqArray.Rmd", so `vignette("SeqArray")` can work directly

- correct Estimated remaining Time to Complete (ETC) for load balancing
  in `seqParallel()`

BUG FIXES

- `seqBED2GDS(, verbose=FALSE)` should have no display

CHANGES

- use a svg file instead of png in vignettes

Changes in version 1.24.2

NEW FEATURES

- add the compiler information in `seqSystem()`

- new arguments '.balancing', '.bl_size' and '.bl_progress' in
  `seqParallel()` for load balancing

UTILITIES

- improve unix forking processes for load balancing in `seqParallel()`

BUG FIXES

- fix `seqSummary()` when no phase data

[seqCAT](/packages/seqCAT)
------

Changes in version 1.8.0

Features

- The Python implementation of profile creation has been deprecated

- SNV profiles are now stored as data frames, rather than GRanges
  objects

- The `create_profile` function now now reads profiles into memory,
  instead of storing them on disk. De-duplication and removal of
  mitochondrial variants is now also performed at this stage

- The `create_profiles` function now returns a list of data frames

- A new function, `write_profile`, can write profiles to disk

- The `read_profile` function now reads profiles without performing any
  de-duplication or removing mitochondrial chromosomes

- The `compare_profiles`, `compare_many` and `list_variants` functions
  now converts input profiles to GRanges internally

- Keep the FILTER column in created SNV profiles

- Add a check for the creation of zero-variant profiles

- Add a check for the existance of the specified input samples

- Removal of non-standard chromosomes and variant de-duplication is now
  optional, and filtration documentation has been extended

[sevenbridges](/packages/sevenbridges)
------------

Changes in version 1.15.1

Improvements

- Added new fields created_by, created_on, and modified_on to the
Project class following the recent API improvements. This enables
better project filtering when querying projects. See the vignette for
details.

[SharedObject](/packages/SharedObject)
------------

Changes in version 0.99.2

- Clean some formating issues

Changes in version 0.99.0 (2019-05-16)

- Submitted to Bioconductor

[signatureSearch](/packages/signatureSearch)
---------------

Changes in version 1.1.0 (2019-10-23)

- Initial version

Changes in version 0.99.20 (2019-10-22)

- Submitted to Bioconductor - Major changes - used HDF5 file to read
and write the matrix in batches - data stored in ExperimentHub

[signeR](/packages/signeR)
------

Changes in version 1.11.1

- ignore mutations on edges of chromosomes

[SIMLR](/packages/SIMLR)
-----

Changes in version 1.11.1 (2019-10-13)

- Fix C++ bug.

[SingleCellExperiment](/packages/SingleCellExperiment)
--------------------

Changes in version 1.8.0

- Added altExp() and related methods to get and set alternative
  Experiments.

- Added the splitAltExps() utility to create many alternative
  Experiments at once.

- Added the swapAltExp() utility to swap between main and alternative
  Experiments.

- Deprecated isSpike(), spikeNames() and related arguments for handling
  spike-ins, in favor of representing spike-ins as alternative
  Experiments.

- Deprecated type= in sizeFactors() and sizeFactorNames(), which were
  previously only required to store size factors for spike-ins.

- Internal change to the representation of reducedDims() to streamline
  subsetting and combining.

[singscore](/packages/singscore)
---------

Changes in version 1.5.1

- added link to the F1000research published workflow that
demonstrates usage of singscore on a real dataset
(https://f1000research.com/articles/8-776/v2) - allow labelling of
samples on dispersion plots (plotDispersion) in a similar manner to
landscape plots (projectScoreLandscape)

[sitePath](/packages/sitePath)
--------

Changes in version 1.1.10

- Allow user to choose whether to show tip labels in the plot functions

Changes in version 1.1.9

- Add progress bar for the resampling and summarizing step of the
  function 'multiFixationSites'

Changes in version 1.1.8

- Add 'plot' function for directly plotting the return of 'extractSite'

- Apply resampling method for 'multiFixationSites'

- The function 'fixationSites' applys the old 'multiFixationSites'

Changes in version 1.1.7

- Add functionality 'extractSite' to allow accessing a single site from
  the result of 'fixationSites' and 'multiFixationSites'

Changes in version 1.1.6

- Add functionality 'setSiteNumbering' to allow manipulating the
  reference of site numbering

Changes in version 1.1.5

- Move similarity calculation to 'addMSA'. This will slow the function

Changes in version 1.1.4

- Use total number of tips divided by number of nodes as
  'minEffectiveSize'

- Ignore invariant sites when search for fixation sites

Changes in version 1.1.3

- Expose 'searchDepath' for 'multiFixationSites'

Changes in version 1.1.2

- Use 'multiFixationSites' for both single and multiple fixation sites

Changes in version 1.1.1

- Bug fix: Error when adding new result in 'fixationSites'

[SNPRelate](/packages/SNPRelate)
---------

Changes in version 1.20.0

- a leading tilde in the file path is allowed in `snpgdsGDS2BED()`

- change the file name "vignettes/SNPRelateTutorial.Rmd" to
  "vignettes/SNPRelate.Rmd", so `vignette("SNPRelate")` can work
  directly

Changes in version 1.18.1

- support long vector in `snpgdsIBDSelection()`

[sojourner](/packages/sojourner)
---------

Changes in version 0.99.8

- added unit test - createTrackll() and msd() in this version

Changes in version 0.99.7

- fixed BiocCheck error checkForSupportSiteRegistration - unblocked
maintainer email account

Changes in version 0.99.6

- fixed BiocCheck warning - revised sojourner.r

Changes in version 0.99.5

- fixed various issues raised in review for Bioconductor submission.

Changes in version 0.99.4

- fixed Linux build error and removed extraneous .pdf files

Changes in version 0.99.3

- fixed git clone issues for Windows

Changes in version 0.99.2

- reduced time runing for examples

Changes in version 0.99.1

- changed system.time() format

Changes in version 0.99.0

- Submitted to Bioconductor

[Spaniel](/packages/Spaniel)
-------

Changes in version 0.99.0

- Submitted to Bioconductor

[SpatialCPie](/packages/SpatialCPie)
-----------

Changes in version 1.1

- (1.1.9) Added example on how to customize resolutions and assignment
  function to vignette

- (1.1.8) Changed UI to make better use of screen real estate

- (1.1.7) Added "top genes" barplot in array plot tooltip

- (1.1.6) Scoring is now based on normalized distances

- (1.1.5) Added interactivity to array plots

- (1.1.4) Increased numerical stability

- (1.1.3) Added log messages

- (1.1.2) Performance improvements

- (1.1.1) Added "top features" barplot in cluster tree tooltip

[splatter](/packages/splatter)
--------

Changes in version 1.10.0 (2019-10-20)

- Add the (experimental) Kersplat simulation model. This model
  incorporates a gene network and other useful features.

- Refactor the summariseDiff function and add the KS statistic.

- Add variable gene correlation plot to compareSCEs and violins to
  other comparison plots.

- Check for counts assay when estimating from SingleCellExperiment
  objects.

- Fix where simpleSimulate stores parameters.

- Fix bugs where parameters were not being passed correctly in
  BASiCSEstimate and sparseDCEstimate.

- Replace the sc_example_counts dataset from scater with the mockSCE
  function.

- Tidy and improve estimation function examples and add checks for
  suggested packages.

- Various fixes for compatibility with updates to other packages.

[SQLDataFrame](/packages/SQLDataFrame)
------------

Changes in version 1.0.0 (2019-10-01)

NEW FEATURES

- Supports tidy grammar for 'select', 'filter', 'mutate' and '%>%'
  pipe.

- Supports representation and saving of MySQL database tables.

- Supports lazy cross-MySQL database table aggregations, such as join,
  union, rbind, etc.

BUG FIXES

- Fixed bugs for single square bracket subsetting with key column(s).

Changes in version 0.99.0 (2019-04-05)

- Submitted to Bioconductor

[sRACIPE](/packages/sRACIPE)
-------

Changes in version 1.1.3 (2019-10-17)

- GeneEx added to inst

Changes in version 1.1.2 (2019-10-17)

- Time series option added to sracipeSimulate

- Minor bug fixes

[ssrch](/packages/ssrch)
-----

Changes in version 1.1.1

- parseDoc did not handle doctitles correctly, fixed

- DocSet() did not succeed in 1.0+, new defaults added to allow this

- parseDoc had presumptive elimination of rownames-associated columns
  generated by read.csv, which was removed in 1.1.1

- increased testing for DocSet updating via parseDoc

- added cautions about parseDoc updating to man page

- for 1.1.2, adding title as a string into searchable 'token' set

[statTarget](/packages/statTarget)
----------

Changes in version 2.0

NEW FEATURES

- New GUI o Mouse Hover for help information o .log file

- New Signal correction o Combat for QC-free Signal correction o
  QC-RFSC methods for metabolomics and proteomics data

- New feature slection o Random Forest and the Permutation based
  variable importance measures o new MDSplot for Random Forest o
  P-value based importance plot

- New data preprocessing o PQN/SUM/none normalization o center/none
  Scaling method

Changes in version 1.15.4

- edit colum name for the result table

- skip the roc analysis once the number of samples in any groups was
  less than 5.

Changes in version 1.15.3

- odds.ratio.  Add the ratio of the odds of an event occurring in one
  group (A) to the odds of it occurring in reference group (B), such as
  odds_of_AtoB.

[Structstrings](/packages/Structstrings)
-------------

Changes in version 1.1.6 (2019-09-23)

- DotBracketDataFrame refactored to split concept from data
  implementation

- fixed sequence column not returned when calling getBasePairing with a
  StructuredXStringSet

[SummarizedBenchmark](/packages/SummarizedBenchmark)
-------------------

Changes in version 2.3.7 (2019-09-18)

- Bug fix: update call on `tidyr::gather` in `tidyUpMetrics` to work
  after bump in `tidyr` version to 1.0.0 on CRAN

Changes in version 2.3.6 (2019-09-06)

- Replace `sc_example_counts` dataset from `scater` with call to
  `scater::mockSCE`

- Fix error thrown by `show.BDData` when dataset is a list

Changes in version 2.3.5 (2019-06-22)

- Fix index title for "Feature: Error Handling" vignette

- Replace Rd references to `rlang::quos` with `rlang:quotation` topic
  page since aliasing appears to be causing warnings on Windows build

- Update pkgdown site with latest devel docs

Changes in version 2.3.4 (2019-06-18)

- Fix scRNAseq simulation case study vignette to work with changes in
  `scRNAseq` package

- Add code examples and import statements to pass BiocCheck

- Update pkgdown site with latest devel docs

Changes in version 2.3.3 (2019-06-14)

- Added NEWS file to track changes

- Incorporated major updates to package documentation and vignettes

- Added pkgdown site

[SummarizedExperiment](/packages/SummarizedExperiment)
--------------------

Changes in version 1.16.0

NEW FEATURES

- Some improvements to the SummarizedExperiment() constructor (see
  commit 0d74843c)

- Support 'colData(SummarizedExperiment) <- NULL' to clear colData

SIGNIFICANT USER-VISIBLE CHANGES

- All the arguments of the SummarizedExperiment() constructor are now
  visible (no more ellipsis) and have default values. So tab completion
  works. See commit 0d74843c

- The dimnames on the individual assays of a SummarizedExperiment
  derivative now can be anything (see issue #25 for the details)

BUG FIXES

- Some fixes to the SummarizedExperiment() constructor (see commit
  0d74843c)

- Address all.equal() false positives on SummarizedExperiment objects
  (see issue #16 for the details)

[SynMut](/packages/SynMut)
------

Changes in version 1.1.3 (2019-05-08)

- Revise dinu_to.keep algorithm, enhance performance.

Changes in version 1.1.2 (2019-05-08)

- Bug fix: "dinu_to" ifelse issue in get_optimal_codon.

Changes in version 1.1.1 (2019-05-06)

- Bug fix: "dinu_to" fix wrong result with "keep == TRUE" parameter.

Changes in version 1.1.0 (2019-05-06)

- Bump

[target](/packages/target)
------

Changes in version 0.1.0

- Submit to Bioconductor

[TargetSearch](/packages/TargetSearch)
------------

Changes in version 1.42.0

SIGNIFICANT USER-VISIBLE CHANGES

- As announced in version 1.40.0, the graphical user interface is gone
  for good (i.e., defunct). The source code is available in my github
  repository.

NEW FEATURES

- The most interesting feature is the introduction of a custom CDF-4
  format which hold the same data as a normal CDF-3 (as exported by the
  software vendors), but allows faster read-access (specially for
  plotting) and compression (among other features). This is at the cost
  of compatibility as the CDF-4 files are unlikely to be used outsied
  TargetSearch

- A new baseline correction method based on quantiles around a
  retention time window. In addition, the new CDF-4 file format allows
  storing of baseline-corrected values so it is not needed to recompute
  the baseline each time like in older TargetSearch versions.

- New function to transform to nominal mass. Some GC instruments export
  CDF not in nominal mass format (some even export high mass accuracy).
  Formely, this type of files were not supported and TargetSearch would
  refuse to process them. Now, all types of mass accuracy are allowed,
  obviously at the cost of losing that accuracy.

BUG FIXES

- Mostly code refactoring and house-keeping.

[TCGAutils](/packages/TCGAutils)
---------

Changes in version 1.6.0

New features

- oncoPrintTCGA: Create an oncoPrint visualization for mutation data
- Support aliquot_ids as input to UUIDtoBarcode function - Additional
sections in the vignette: CpGtoRanges, UUIDtoBarcode for aliquot_ids
- TCGAprimaryTumors allows users to select all primary tumors for a
given curatedTCGAData MultiAssayExperiment object (@vjcitn)

Minor changes and bug fixes

- Now merging clinical data using both rows and columns in
mergeColData - Added informative error when query results are empty
in UUIDtoBarcode - Updates to makeGRangesListFromExonFiles to use
S4Vectors::splitAsList (@hpages)

[TOAST](/packages/TOAST)
-----

Changes in version 0.99.7

- Add options to improve feature selections in reference-free
  deconvolution.

Changes in version 0.99

NEW FEATURES

- Initial release.

[tofsims](/packages/tofsims)
-------

Changes in version 099.1

SIGNIFICANT USER-VISIBLE CHANGES

- changed function behvaiour in the whole package from call-by-ref to
  call-by value. Adjusted accordingly all examples and the vignette.

INTERNALS

- depends now on ProtGenerics from which it uses 'mz'

- exchanged various print() with message()

[topconfects](/packages/topconfects)
-----------

Changes in version 1.1.4 (2018-12-29)

- Prepare for Bioconductor submission.

Changes in version 1.1.3 (2018-11-14)

- Add "full" output option for limma_confects and normal_confects,
  which adds se, df, and fdr_zero columns.

Changes in version 1.1.2 (2018-06-30)

- Add DESeq2 support.

- Code to do with non-linear effect sizes has been move to the "ql"
  branch on github, pending publication of the basic method and a
  possible rethink and use of a simpler method.

Changes in version 1.1.1 (2019-09-20)

- Update citation information.

Changes in version 1.0.1 (2018-02-04)

- Initial release.

[topdownr](/packages/topdownr)
--------

Changes in version 1.7

- New version for Bioc 3.10 (devel)

Changes in version 1.7.1

- Remove NEWS file (just keep NEWS.md).  - Never remove "AgcTarget"
column from colData DataFrame.  - Strip white spaces from
ScanHeadsman output.  - Defunct defaultMs1Settings and
defaultMs2Settings. They will be removed in 3.11 &#91;2019-06-19&#93;.

Changes in version 1.7.2

- Add readTopDownSet(..., conditions="ScanDescription") as a new way
to read scan conditions (see #80/#81) &#91;2019-08-08&#93;.

[TPP](/packages/TPP)
---

Changes in version 3.13.3

- defunct the following functions: tpp2dPlotCCRGoodCurves,
  tpp2dPlotCCRAllCurves, tpp2dPlotCCRSingleCurves,
  tpp2dEvalConfigTable, tpp2dRemoveZeroSias, tpp2dReplaceColNames,
  tpp2dCreateCCRConfigFile

Changes in version 3.13.2

- Revise vignette

Changes in version 3.13.1

- Avoid warnings due to factor/character conversions in vignette

[TPP2D](/packages/TPP2D)
-----

Changes in version 1.1.0 (2019-05-06)

- Carry-over detection now works via a LDA model

[trackViewer](/packages/trackViewer)
-----------

Changes in version 1.21.18

- update NEWS file from "CHANGES IN VERSION" to "Changes in version".

Changes in version 1.21.17

- remove google scholar in vignette.

Changes in version 1.21.16

- Try to avoid the error introduced by BiocStyle: Invalid Parameter -
  /figure-html.

Changes in version 1.21.15

- Add google scholar in vignette.

Changes in version 1.21.14

- Update the citation in readme file.

Changes in version 1.21.13

- Fix a typo in plotLollipops.

Changes in version 1.21.12

- setTrackXscaleParam accept position attribute.

Changes in version 1.21.11

- set rescale of lolliplot by precentage.

Changes in version 1.21.10

- fix the bug of score is NumericList.

Changes in version 1.21.9

- accept negative values for bigwig files.

Changes in version 1.21.8

- update legend method.

Changes in version 1.21.7

- fix the label label.parameters.gp.

Changes in version 1.21.6

- add sample code of Proteins API.

Changes in version 1.21.5

- fix the bug that lollipop stem is too long.

Changes in version 1.21.4

- add shape for lolliplot.

Changes in version 1.21.3

- add shinyApp video.

Changes in version 1.21.2

- add citation.

Changes in version 1.21.1

- fix the bug if condiction with multiple logical values.

[tradeSeq](/packages/tradeSeq)
--------

Changes in version 0.99.9902 (2019-10-23)

- tradeSeq's fitGam and evaluateK functions do not have a seed argument
  anymore, as per Bioconductor's guidelines. Users are encouraged to
  set the seed manually with the set.seed function before running those
  functions for reproducibility purposes.

Changes in version 0.99.47 (2019-09-02)

- tradeSeq now provides `singleCellExperiment` output o fitGam now
  accepts a `slingshotDataSet` object as input o All tests and plotting
  functions accept a `singleCellExperiment` object that contain
  tradeSeq output

Changes in version 0.99.0 (2019-06-22)

- Submitted to Bioconductor

Changes in version 0.9.0 (2019-03-15)

- Reformatted to fulfill Bioconductor guidelines

[transite](/packages/transite)
--------

Changes in version 1.2.1

- updated author list

[treeio](/packages/treeio)
------

Changes in version 1.9.3

- add citation information (2019-10-05, Sta) - rename phyPML to
as.treedata.pml (2019-10-01, Tue) - as.phylo method for igraph (only
work with tree graph) (2019-09-28, Sat)

Changes in version 1.9.2

- nodeid and nodelab methods for converting between node number and
labels (2019-08-09, Fri) - parent, 'ancestor, child,
offspringandrootnodemethods fortreedata` (2019-08-07, Wed) -
read.mega_tabular to parse MEGA Tabular output (2019-07-16, Tue) -
read.mega to parse MEGA NEXUS (actually BEAST compatible)

Changes in version 1.9.1

- rename_taxa now use 1st column as key and 2nd column as value by
default (2019-05-28, Tue) - enable tree_subset to specify group_name
and enable to incorporate root.edge by setting root_edge = TRUE
(2019-05-27, Mon) - full_join method for phylo object (2019-05-22,
Wed) - redefined root method to wrape ape::root.phylo for
compatibility (2019-05-20, Mon) -
https://github.com/GuangchuangYu/treeio/issues/18

[tRNAscanImport](/packages/tRNAscanImport)
--------------

Changes in version 1.5.3 (2019-08-25)

- added get.tRNAprecursor function to retrieve tRNA precursor sequences
  in combination with genomic sequences

- fixed typos in the NEWS file

- updated tRNAscan example file for human (the high confidence set is
  now included)

[TSRchitect](/packages/TSRchitect)
----------

Changes in version 1.11.12

- Updates to include 3'-read capturing from bam input files and data
  structures to allow ignoring spurious TSSs.

- Added capabilities for sorting by strand (and seq/TSS) for TSS
  merging.

- Made a change leading to a substantial speed-up of mergeSampleData().

- Added qname to saved columns from .bam input.

- Minor changes to writeTSR() and writeTSS().

- Rewrite to mergeSampleData).

- Improvements to loadTSSobj().

- Necessary updates to the Singularity recipe to reflect the above
  changes.

[tximeta](/packages/tximeta)
-------

Changes in version 1.4.0

- tximeta will now pull down RefSeq seqinfo, using the dirname() of the
  GTF location, and assuming some consistency in the structure of the
  assembly_report.txt that is located in the same directory. Needs more
  testing though across releases and organisms.

- expanded caching of ranges to exons and genes as well. Exons in
  particular take a long time to build from TxDb, so this saves quite a
  lot of time.

- new 'addExons' function will add exons to trancript-level summarized
  experiments, by replacing transcript GRanges with exon-by-transcript
  GRangesList. Purposely designed only for transcript-level, see note
  in ?addExons

- tximeta now also caches the transcript ranges themselves, rather than
  just the TxDb. This shaves extra seconds off the tximeta() call!

- add 'skipSeqinfo' argument, which avoids attempting to fetch
  chromosome information (from UCSC) if set to TRUE.

[tximport](/packages/tximport)
--------

Changes in version 1.14.0

- Alevin count and inferential variance can be imported now ~40x faster
  for large number of cells, leveraging C++ code from the fishpond
  package (>= 1.1.18).

- Alevin inferential replicates can be imported (also sparse). To not
  import the inferential replicates, set dropInfReps=TRUE.

[Ularcirc](/packages/Ularcirc)
--------

Changes in version 1.3.22

NEW FEATURES

- New stand alone functions that don't require shiny app :
  bsj_to_circRNA_sequence() and bsj_fastq_generate()

BUG FIXES

- BUG fix: All warnings in bioconductor now resolved

- BUG fix: Added dependencies to DESCRIPTION to allow automatic install

Changes in version 1.3.1

BUG FIXES

- BUG fix: Fixed code in Compatible_Annotation_DBs examples that was
  causing error in bioconductor builds.

[universalmotif](/packages/universalmotif)
--------------

Changes in version 1.4.0

NEW FEATURES

- scan_sequences(..., threshold.type) option: 'logodds.abs'. Allows the
  exact threshold scores to be provided.

- compare_motifs() option: 'min.position.ic'. Prevent low-IC positions
  in an alignment from contributing to the final alignment score.

- compare_motifs() option: 'score.strat'. Instruct the function how to
  deal with individual column scores in an alignment. This is also
  replaced the old way of choosing between sum and mean via prepending
  an 'M' to the metric name. Strategies for combining column scores
  include: sum, arithmetic mean, geometric mean, median, Fisher
  Z-transform, and weighted means.

- Motif comparison metrics: average log-likelihood ratio, squared
  Euclidean distance, Hellinger distance, Bhattacharyya coefficient,
  Manhattan distance, lower limit average log-likelihood ratio,
  weighted Euclidean distance, weighted Pearson correlation
  coefficient.

- compare_columns() utility: Compare two 1d numeric vectors using the
  comparison metrics from compare_motifs().

- compare_motifs() option: 'output.report'. Generate an output report
  when 'compare.to' is provided, showing motif alignments of top
  matches.

- get_scores() utility: Extract all possible scores from a motif.

- filter_motifs(): Filter using the 'extrainfo' slot.

- MotifComparisonAndPvalues.pdf vignette: the comparisons and P-values
  sections have been moved from AdvancedUsage.pdf to their own
  vignette. Higher order motifs, enrichment and run_meme() usage
  sections have been moved to SequenceSearches.pdf.

MINOR CHANGES

- Removed 'random' shuffling method.

- Using RcppThread instead of BiocParallel in several functions:
  compare_motifs(), create_sequences(), get_bkg(), motif_pvalue(),
  scan_sequences(), shuffle_sequences(). This means parallelization can
  occur within C++ code which is much faster than having to jump
  between R and C++. Currently motif_peaks(), read_motifs() and
  write_motifs() are the only remaining functions which offer optional
  BiocParallel usage.

- Many performance improvements to functions relying on internal C++
  code. Several internal R functions have been replaced with C++
  versions.

- Changed behaviour of make_DBscores() and motif comparison P-values.
  Re-calculated internal P-value databases.

- For merge_motifs(..., use.type): now only accepts 'PPM'.

- When comparing all motifs to all motifs with any method in
  compare_motifs(), the diagonal entries now properly show the max/min
  possible similarity/distance scores.

- New internal merge_motifs() implementation. This also fixes a
  previous bug with incorrect PPM averaging.

- read_homer(): the logodds score is converted to a P-value.

- motif_pvalue(): New score calculator. Exact scores are still
  calculated the same (but with a faster C++ function), but approximate
  scores are now calculated by randomly generating score distributions
  from size 'k' motif score blocks.

- motif_pvalue(): Added a safety check when trying to use this function
  with large motifs. Will throw a warning when nrow(matrix)^k > 1e8 and
  reduce k accordingly before continuing.

- Adjusted P-value calculation in motif_peaks() to not display Pval = 0
  so easily by instead estimating a normal distribution from random
  peaks.

- convert_type(): make sure not to leave any zeros in bkg vector when a
  pseudocount greater than zero is used.

- enrich_motifs(): split up 'hits' and 'positional' resuts into their
  own data.frames.

- Replaced several instances of cat() with message() for printing
  progress updates.

- Positional tests have been removed from enrich_motifs(). See
  motif_peaks() for testing motif-sequence positional preferences.

- In read_meme(): E-values are now additionally stored in the extrainfo
  slot. This is to preserve E-values smaller than the R double
  precision limit.

- In read_transfac(): Matrix values are rounded, to prevent errors when
  reading in matrices with non-integers.

- Update JASPAR2018_CORE_DBSCORES with new compare_motifs() methods and
  params.

- universalmotif print() method now returns the object invisibly,
  instead of NULL.

BUG FIXES

- read_meme() will now properly parse background letter frequencies
  which span more than one line.

- convert_motifs() will not error-out when trying to convert a PFMatrix
  with a family character vector longer than one.

- Fixed P-value calculation when importing HOMER motifs. Peviously it
  would simply assume the log threshold value was the P-value. Now
  motif_pvalue() is used to properly calculate a P-value.

Changes in version 1.2.1

BUG FIXES

- Mispelled variable in enrich_motifs()

[variancePartition](/packages/variancePartition)
-----------------

Changes in version 1.15.8

- Replace cat() with message()

- add quiet option to a few functions

- dream() does not call eBayes() when lmFit is used

Changes in version 1.15.7

- Official release to development branch

Changes in version 1.15.6

- fix convergence errror when recycling parameters values from first
  gene

- add column z.std and F.std to topTable

Changes in version 1.15.4

- Add error message when scale of fixed effects causes a problem

Changes in version 1.15.3

- Try changing order of eBayes

Changes in version 1.15.2

- Update vignette

Changes in version 1.15.0

- Push changes to Bioconductor devel

[VariantExperiment](/packages/VariantExperiment)
-----------------

Changes in version 0.99.8

- `saveVariantExperiment` returns newly generated VariantExperiment
  object by calling `loadVariantExperiment`

[ViSEAGO](/packages/ViSEAGO)
-------

Changes in version 0.99.42

- remove environment signature in merge_enrich_terms

Changes in version 0.99.41

- keep require("topGO") in create_topGOdata method

Changes in version 0.99.40

- add attachNamespace("topGO") instead require("topGO") in
  create_topGOdata method for prevent build check warnings

Changes in version 0.99.39

- check gene background in merge_enrich_terms upgrade

Changes in version 0.99.38

- show method for GO_cluster print bug correction

Changes in version 0.99.37

- add environment option for merge_enrich_terms

Changes in version 0.99.36

- annotate id argument signature with "NULL"

Changes in version 0.99.35

- README update

Changes in version 0.99.34

- README update (download from Bioconductor and forgeMIA gitlab)

- add require("topGO") in create_topGOdata method

Changes in version 0.99.33

- merge_enrich_terms exemple update.

Changes in version 0.99.32

- GOterms_heatmap GO.clusters column convert to character in order to
  prevent data.table merge error.

Changes in version 0.99.30

- ViSEAGO citation update

[waddR](/packages/waddR)
-----

Changes in version 0.99.6 (2019-10-28)

- P-value reproducibility feature: o Now using nextRNGstream to
  generate independent seeds for each gene o Setting the seed (with
  set.seed()) only in .testWass o Removed use of seeds and the seed
  argument everywhere in WassersteinTest.R

Changes in version 0.99.5 (2019-10-14)

- P-value reproducibility feature for the permutation procedure: o as
  described
  [here](https://github.com/Bioconductor/Contributions/issues/1218#issuecomment-53710763),
  a seed argument has been added to the functions wasserstein.test and
  wasserstein.sc.

- Fix: o Checks for the validity of the "method" argument in
  wasserstein.sc and wasserstein.test that have become unnecessary with
  the use of match.arg have been removed

Changes in version 0.99.4 (2019-10-01)

- Fixes: o Fixed a bug in .wassersteinTestSp where the names in the
  output vector were changed unexpectedly and added a test for this bug

- Bioc Review I: o vignette: Added SessionInfo() to each vignettes o
  vignette/README: Changed the install instructions o unit tests:
  removed unused and commented-out code o R: Changed to switch
  statements to dispatch different methods in wasserstein.test and
  wasserstein.sc o R: Changed the order of arguments in
  wasserstein.test and wasserstein.sc and added default methods o R:
  wasserstein.test.sp has been renamed to .wassersteinTestSp;
  wassersetin.test.asy has been renamed to .wassersteinTestAsy -> both
  are now private o R/NAMESPACE: removed the previously private
  functions .fishersCombinedPval and .combinePVal from NAMESPACE by
  removing @export decorators

Changes in version 0.99.3 (2019-08-30)

- Fixes: o Fixed a bug in wasserstein.test that led to NAs during gpd
  fitting o Fixed a bug in .gpdFittedPValue that led to NAs during gpd
  fitting

- Modified wasserstein.sc tests and added new tests to reproduce the
  bugs and challenge the fixes

- Change in squared_wass_decomp: If the standard deviation of one
  condition is 0, quantile-quantile corelation is not computed, since
  the shape term would be zero anyway. Previously, NAs were produced in
  some cases.

- Swapped 'true' and 'test' values in call to .relativeError

Changes in version 0.99.2 (2019-08-27)

- Changes to DESCRIPTION: o removed the "Maintainer" field in favor of
  "Authors@R" to address build warnings o Changed roles: 'cre' =
  Maintainer, 'aut' = author

Changes in version 0.99.0 (2019-08-27)

- Version Bump for BioC Submission

- Vignettes: o Added introduction section and sessionInfo() to main
  vignette o Added link to main vignette

- Code style: o renaming where possible to conform with BioC convention
  o file renaming to uniform upper camel-case o comments edited o
  additional helper functions introduced

Changes in version 0.2.8 (2019-08-19)

- R Code redesign: o new unexported functions to help avoiding repeats
  o Redesign of all single-cell methods as S4 methods that are also
  capable to take SingleCellExperiment objects as input o Bioc-style
  function naming implemented o Now using the cpp decomp functions
  instead computations in R

- Output names changed

- Descriptions changed to match altered code

Changes in version 0.2.7 (2019-08-13)

- Bug fix in interval table / wasserstein_metric causing wasserstein.sc
  runs to fail

- Tests for that bug added

Changes in version 0.2.6 (2019-08-08)

- Fix in wasserstein_metric where a result was squared

- Work on CPP implementations: o Rework on the quantile computation in
  CPP that now produces more accurate approximations o Removal of
  obsolete parameter "p" in approximation functions

- Now using wasserstein_metric in two sample testing procedures

Changes in version 0.2.5 (2019-08-06)

- Fixing R CMD check size error by reducing package size: Brownian
  bridge distribution, used in the asymptotic implementation of
  wasserstein.test now is downloaded into a local cache during the
  first run of wasserstein.test. From there it is loaded in all
  subsequent runs.

- Fixed Bug in wasserstein.test

Changes in version 0.2.4 (2019-07-31)

- Added Vignettes for wasserstein distance, wasserstein.test, and
  wasserstein.sc

- Fixed examples and unparsable comments

- Adressing Notes of R CMD BiocCheck ...: o Removed use of 1:... in
  favor of seq() or seq_len() o Removed use of set.seed(), which
  changed the signature of testing functions o Using 4 spaces instead
  of tabs and multiples of 4 spaces for indentation

Changes in version 0.2.3 (2019-07-29)

- Adressing Notes and Warnings of R CMD check ...: o Reverted changes
  of BiocParallel from Import to Enhancement o Removed redundant
  imports of arm, eva, BiocParallel o Removed CITATION, as it could be
  formatted to satisfy the checks without naming a paper

Changes in version 0.2.2 (2019-07-29)

- Added this NEWS file for change announcements

- Added a file inst/CITATION that is supposed to hold the citation
  (after publication)

- DESCRIPTION file: o Title and Description improved and shortened o
  Added bioView Categories: StatisticalMethod, SingleCell,
  DifferentialExpression o Added BugReports and URL in DESCRIPTION o
  Changed the former Import BiocParallel to an Enhancement

- Changes to NAMESPACE: o Explicitly declare the functions that should
  be imported from the packages arm, eva, and BiocParallel

- Changed all code files in the package to have max. 80 Character per
  line

[xcms](/packages/xcms)
----

Changes in version 3.7.5

- Remove xcmsMSn vignette (based on old xcms).

Changes in version 3.7.4

- mzClust correspondence analysis: check and fix missing values in
  column mz of the peaks matrix (issue #416).

Changes in version 3.7.3

- plot type = "XIC" on an XCMSnExp object will draw rectangles
  indicating the identified chromatographic peaks.

- Add a vignette describing LC-MS/MS data analysis with xcms.

Changes in version 3.7.2

- Fix documentation (issue #401).

- Add support for SWATH data analysis.

Changes in version 3.7.1

- Add correlate method for Chromatogram objects.

- Add parameter lwd to plotAdjustedRtime.

- Add align method for Chromatogram objects.

- Add findChromPeaksIsolationWindow to enable chromatographic peak
  detection in isolation windows.

- Fix issue in chromPeakSpectra with method = "signal".

- chromPeakSpectra and featureSpectra return now MS2 spectra with an
  precursor m/z >= mzmin, <= mzmax and retention time >= rtmin, <=
  rtmax.

- Improve performance of chromPeakSpectra and featureSpectra.

[zinbwave](/packages/zinbwave)
--------

Changes in version 1.7.5 (2019-10-08)

- Changed default of `zinbwave` to `observationalWeights=FALSE` to
  speed up computations when weights are not needed.

- Added argument `zeroinflation = TRUE`: when set to FALSE a negative
  binomial model is fit.

- Removed dependence on the `copula` package to avoid depending on
  `gsl`.


NEWS from new and existing Data Experiment Packages
===================================


[benchmarkfdrData2019](/packages/benchmarkfdrData2019)
--------------------

Changes in version 0.99.15 (2019-06-06)

- Updated vignette to remove excessive printing of data objects

Changes in version 0.99.14 (2019-05-08)

- bioc4: Added non-eval code chuck to vignette for installing package

Changes in version 0.99.13 (2019-05-08)

- Changes based on feedback from Bioconductor review

- bioc1: Added 'BugReports:' field to DESCRIPTION file

- bioc2: Removed reference to 'remotes' pkg from vignette

- bioc3: Added details to vignette for alt approach to load resource
  from ExperimentHub

Changes in version 0.99.12 (2019-05-07)

- Explicitly add metadata param documentation to each Rd

- Update description in DESCRIPTION file

Changes in version 0.99.11 (2019-05-06)

- Fixed typo in one doc reference again

Changes in version 0.99.10 (2019-05-06)

- Added `metadata` parameter to docs

- Changed examples to just load metadata

Changes in version 0.99.9 (2019-05-06)

- Added running examples to all docs

- Fixed typo in one doc reference

Changes in version 0.99.8 (2019-05-04)

- Added documentation for all data objects with roxygen2 code in R
  files

Changes in version 0.99.7 (2019-04-24)

- Changed NEWS message for 0.99.6 to prevent BiocCheck warning

Changes in version 0.99.6 (2019-04-24)

- Changed T/F to TRUE/FALSE in make-data Rmd files

- Changed how object class is checked in plotMethodRank code

Changes in version 0.99.5 (2019-04-23)

- Deleted and untracked .Rproj file

Changes in version 0.99.4 (2019-04-23)

- Removed .Rproj file from git tracking

Changes in version 0.99.3 (2019-04-23)

- Fixed typo in DESCRIPTION file

- Added NAMESPACE file

Changes in version 0.99.2 (2019-04-23)

- Triggering new build after adding webhook

Changes in version 0.99.1 (2019-04-23)

- Updating vignette and triggering new build after data have been moved
  into ExperimentHub

- Adding CITATION and NEWS files

- Adding ORCID for authors

Changes in version 0.99.0 (2019-04-09)

- Submitted to Bioconductor

[chipenrich.data](/packages/chipenrich.data)
---------------

Changes in version 2.10.0

- Added data required for proxReg: enhancers.dnase_thurman.0,
  gene.enh.desc, spline.log_dtss.90ENCODE

[depmap](/packages/depmap)
------


Changes in version 0.99.6

- 19Q3 data
- New data loading functions



Changes in version 0.99.0

- Preparing for Bioconductor submission


[derfinderData](/packages/derfinderData)
-------------

Changes in version 2.3.4

- Added a NEWS.md file to track changes to the package.

[DMRcatedata](/packages/DMRcatedata)
-----------

Changes in version 2.0.0

- Data package now uses ExperimentHub

- Objects myBetas, CpGs, tx.hg19, tx.hg38 and tx.mm10 have been retired

- Objects hg19.generanges, hg19.grt, hg38.generanges, hg38.grt,
  mm10.generanges and mm10.grt have been added as annotation for
  extractRanges() and DMR.plot()

[dsQTL](/packages/dsQTL)
-----

Changes in version 2.17

USER VISIBLE CHANGES

- ch2locs (retrievable via dsQTL::getSNPlocs) has been changed at about
  1850 locations where rs numbers had been associated with hg19
  addresses; the dsQTL regions are hg18 as are all the chr2... SNP
  addresses.  Previously the discoverable rs numbers used in the
  Chicago distribution from
  http://eqtl.uchicago.edu/dsQTL_data/GENOTYPES/ had be mapped via
  SNPlocs...20111119, but now they come directly from the Chicago text
  file.

[HDCytoData](/packages/HDCytoData)
----------

Changes in version 1.5.4

- Add Weber_AML_sim and Weber_BCR_XL_sim datasets

- Update documentation

[MMAPPR2data](/packages/MMAPPR2data)
-----------

Changes in version 0.99.20 (2019-05-30)

- Submitted to Bioconductor, nearing end of review process

- See NEWS in MMAPPR2 package for relevant updates

[MouseGastrulationData](/packages/MouseGastrulationData)
---------------------

Changes in version 0.99.0 (2019-06-17)

- Submitted to Bioconductor

[PhyloProfileData](/packages/PhyloProfileData)
----------------

Changes in version 0.99.0 (2019-05-30)

- Submitted to Bioconductor

[PtH2O2lipids](/packages/PtH2O2lipids)
------------

Changes in version 2016-04-21

- Initial release for Bioconductor

[pwrEWAS.data](/packages/pwrEWAS.data)
------------

Changes in version 0.1.0

- Initial version

[RforProteomics](/packages/RforProteomics)
--------------

Changes in version 1.23.1

- remove getThermoHelaPRTC function (server down) <2019-10-24 Thu>

[RNAmodR.Data](/packages/RNAmodR.Data)
------------

Changes in version 0.99.0 (2019-04-29)

- Submitted to Bioconductor

[SBGNview.data](/packages/SBGNview.data)
-------------

Changes in version 0.99.10

- Added gene expression demo datasets.

[scRNAseq](/packages/scRNAseq)
--------

Changes in version 2.0.0

- Added lots of new ExperimentHub datasets, inspired by
  simpleSingleCell use cases and Martin Hemberg's website.

- All outputs are now SingleCellExperiment instances with spike-ins
  stored as alternative experiments.

- Deprecated inbuilt datasets in favor of ExperimentHub equivalents.

[signatureSearchData](/packages/signatureSearchData)
-------------------

Changes in version 1.1.0 (2019-10-23)

- Initial version

Changes in version 0.99.10 (2019-10-22)

- Stored data in ExperimentHub

- Saved LINCS gctx file to hdf5 file in batches

- Saved cmap databases (cmap, cmap_rank, cmap_expr) to hdf5 files

- Loaded hdf5 file into R as SummarizedExperiment object

- hdf5 file includes matrix, rownames and colnames

Changes in version 0.99.0 (2019-04-02)

- Needed 50Gb memory to load matrix in LINCS gctx file and save as HDF5
  backed SE object

[TENxBUSData](/packages/TENxBUSData)
-----------

Changes in version 0.99.0 (2019-01-08)

- Submitted to Bioconductor for review


NEWS from new and existing Workflows
===================================


[BgeeCall](/packages/BgeeCall)
--------

Changes in version 1.1.1

- Better manage transcript version

- Manage output directory when calls are generated for several
  libraries

- Manage transcript version when calls are generated for several
  libraries

Changes in version 1.1.0

- Change kmer size of kallisto index from 21 to 15 for libraries with
  readLength <= 50bp

- Can choose the output directory

- By default use a simpler arborescence for the output directory

- Can use reference intergenic sequences generated by the community

- Can use custom reference intergenic sequences from a local fastq file

[CAGEWorkflow](/packages/CAGEWorkflow)
------------

Changes in version 1.1.5

- Updated main text with suggestions from F1000 reviewers.

- Added F1000 paper as citation.

- Slightly updated code examples to be compliant with new Bioconductor
  release.


Deprecated and Defunct Packages
===============================

Twelve software packages were removed from this release (after being deprecated
in Bioc 3.9): flowQ, rMAT, TSSi, flowQB, rSFFreader, ProCoNA, spliceSites, DOQTL,
NGScopy, SVM2CRM, miRsponge, htSeqTools

Nineteen software are deprecated in this release and will be removed in Bioc 3.11:
SNPchip, rHVDM, GenomeGraphs, plateCore, charm, HTSanalyzeR, PathNet, Rchemcpp,
exomePeak, flipflop, Pbase, RnaSeqSampleSize, birte, SEPA, CNPBayes, dSimer,
mlm4omics, condcomp, brainImageR

Two experimental data packages were removed in this release (after being
deprecated in BioC 3.9): PGPC, flowQBData.

Three experimental data packages are deprecated in this release and will be
removed in Bioc 3.11: charmData, facopy.annot, allenpvc.

Two annotation packages were removed this release: MafDb.gnomADex.r2.0.1.GRCh38,
MafDb.gnomAD.r2.0.1.GRCh38 (they have been replaced with
MafDb.gnomADex.r2.1.GRCh38,
MafDb.gnomAD.r2.1.GRCh38)

Two annotation packages are deprecated in this release and will be removed
in Bioc 3.11: MafDb.ESP6500SI.V2.SSA137.hs37d5, MafDb.ESP6500SI.V2.SSA137.GRCh38
