Bioconductors:

We are pleased to announce Bioconductor 3.0, consisting of 929
software packages, 220 experiment data packages, and 870
up-to-date annotation packages. 

There are 109 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.0 is compatible with R 3.1.1,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release includes an updated Bioconductor [Amazon Machine Image]
(http://bioconductor.org/help/bioconductor-cloud-ami/).

Visit [http://bioconductor.org](http://bioconductor.org)
for details and downloads.

Contents
--------

* Getting Started with Bioconductor 3.0
* New Software Packages
* NEWS from new and existing packages
* Packages removed from the release

Getting Started with Bioconductor 3.0
======================================

To update to or install Bioconductor 3.0:

1. Install R 3.1.1.  Bioconductor 3.0 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 109 new packages in this release of Bioconductor.

ALDEx2 - A differential abundance analysis for the comparison of two or more conditions. For example, single-organism and meta-RNA-seq high-throughput sequencing assays, or of selected and unselected values from in-vitro sequence selections. Uses a Dirichlet-multinomial model to infer abundance from counts, that has been optimized for three or more experimental replicates. Infers sampling variation and calculates the expected false discovery rate given the biological and sampling variation using the Wilcox rank test or Welches t-test (aldex.ttest) or the glm and Kruskal Wallis tests (aldex.glm). Reports both P and fdr values calculated by the Benjamini Hochberg correction.

ASGSCA - The package provides tools to model and test the association between multiple genotypes and multiple traits, taking into account the prior biological knowledge. Genes, and clinical pathways are incorporated in the model as latent variables. The method is based on Generalized Structured Component Analysis (GSCA).

ballgown - Tools for statistical analysis of assembled transcriptomes, including flexible differential expression analysis, visualization of transcript structures, and matching of assembled transcripts to annotation.

blima - Package blima includes several algorithms for the preprocessing of Illumina microarray data. It focuses to the bead level analysis and provides novel approach to the quantile normalization of the vectors of unequal lengths. It provides variety of the methods for background correction including background subtraction, RMA like convolution and background outlier removal. It also implements variance stabilizing transformation on the bead level. There are also implemented methods for data summarization. It also provides the methods for performing T-tests on the detector (bead) level and on the probe level for differential expression testing.

BridgeDbR - Use BridgeDb functions and load identifier mapping databases in R

CFAssay - The package provides functions for calculation of linear-quadratic cell survival curves and for ANOVA of experimental 2-way designs along with the colony formation assay.

ClassifyR - The software formalises a framework for classification in R. There are four stages. Data transformation, feature selection, and prediction. The requirements of variable types and names are fixed, but specialised variables for functions can also be provided. The classification framework is wrapped in a driver loop, that reproducibly does a couple of cross-validation schemes. Functions for differential expression, differential variability, and differential distribution are included. Additional functions may be developed by the user, if they have better performing methods.

compEpiTools - Tools for computational epigenomics developed for the analysis, integration and simultaneous visualization of various (epi)genomics data types across multiple genomic regions in multiple samples.

CoRegNet - This package provides methods to identify active transcriptional programs. Methods and classes are provided to import or infer large scale co-regulatory network from transcriptomic data. The specificity of the encoded networks is to model Transcription Factor cooperation. External regulation evidences (TFBS, ChIP,...) can be integrated to assess the inferred network and refine it if necessary. Transcriptional activity of the regulators in the network can be estimated using an measure of their influence in a given sample. Finally, an interactive UI can be used to navigate through the network of cooperative regulators and to visualize their activity in a specific sample or subgroup sample. The proposed visualization tool can be used to integrate gene expression, transcriptional activity, copy number status, sample classification and a transcriptional network including co-regulation information.

cosmiq - cosmiq is a tool for the preprocessing of liquid- or gas - chromatography mass spectrometry (LCMS/GCMS) data with a focus on metabolomics or lipidomics applications. To improve the detection of low abundant signals, cosmiq generates master maps of the mZ/RT space from all acquired runs before a peak detection algorithm is applied. The result is a more robust identification and quantification of low-intensity MS signals compared to conventional approaches where peak picking is performed in each LCMS/GCMS file separately. The cosmiq package builds on the xcmsSet object structure and can be therefore integrated well with the package xcms as an alternative preprocessing step.

COSNet - Package that implements the COSNet classification algorithm. The algorithm predicts node labels in partially labeled graphs.

csaw - Detection of differentially bound regions in ChIP-seq data with sliding windows, with methods for normalization and proper FDR control.

DEGreport - Creation of a HTML report of differential expression analyses of count data. It integrates some of the code mentioned in DESeq2 and edgeR vignettes, and report a ranked list of genes according to the fold changes mean and variability for each selected gene.

derfinder - Fast differential expression analysis of RNA-seq data at base-pair resolution

derfinderHelper - Helper package for speeding up the derfinder package when using multiple cores.

derfinderPlot - Plotting functions for derfinder

DOQTL - DOQTL is a quantitative trait locus (QTL) mapping pipeline designed for Diversity Outbred mice and other multi-parent outbred populations. The package reads in data from genotyping arrays and perform haplotype reconstruction using a hidden Markov model (HMM). The haplotype probabilities from the HMM are then used to perform linkage mapping. When founder sequences are available, DOQTL can use the haplotype reconstructions to impute the founder sequences onto DO genomes and perform association mapping.

DupChecker - Meta-analysis has become a popular approach for high-throughput genomic data analysis because it often can significantly increase power to detect biological signals or patterns in datasets. However, when using public-available databases for meta-analysis, duplication of samples is an often encountered problem, especially for gene expression data. Not removing duplicates would make study results questionable. We developed a Bioconductor package DupChecker that efficiently identifies duplicated samples by generating MD5 fingerprints for raw data.

EnrichmentBrowser - The EnrichmentBrowser package implements essential functionality for the enrichment analysis of gene expression data. The analysis combines the advantages of set-based and network-based enrichment analysis in order to derive high-confidence gene sets and biological pathways that are differentially regulated in the expression data under investigation. Besides, the package facilitates the visualization and exploration of such sets and pathways.

erccdashboard - Technical performance metrics for differential gene expression experiments using External RNA Controls Consortium (ERCC) spike-in ratio mixtures.

facopy - facopy is an R package for fine-tuned cancer CNA association modeling. Association is measured directly at the genomic features of interest and, in the case of genes, downstream gene-set enrichment analysis can be performed thanks to novel internal processing of the data. The software opens a way to systematically scrutinize the differences in CNA distribution across tumoral phenotypes, such as those that relate to tumor type, location and progression. Currently, the output format from 11 different methods that analyze data from whole-genome/exome sequencing and SNP microarrays, is supported. Multiple genomes, alteration types and variable types are also supported.

FEM - FEM can dentify interactome hotspots of differential promoter methylation and differential ex-pression, where an inverse association between promoter methylation and gene expression is assumed.

flowcatchR - flowcatchR is a set of tools to analyze in vivo microscopy imaging data, focused on tracking flowing blood cells. It guides the steps from segmentation to calculation of features, filtering out particles not of interest, providing also a set of utilities to help checking the quality of the performed operations (e.g. how good the segmentation was). The main novel contribution investigates the issue of tracking flowing cells such as in blood vessels, to categorize the particles in flowing, rolling and adherent. This classification is applied in the study of phenomena such as hemostasis and study of thrombosis development.

flowCHIC - A package to analyze flow cytometric data of complex microbial communities based on histogram images

flowClean - A quality control tool for flow cytometry data based on compositional data analysis.

flowDensity - This package provides tools for automated sequential gating analogous to the manual gating strategy based on the density of the data.

focalCall - Detection of genomic focal aberrations in high-resolution DNA copy number data

FourCSeq - FourCSeq is an R package dedicated to the analysis of (multiplexed) 4C sequencing data. The package provides a pipeline to detect specific interactions between DNA elements and identify differential interactions between conditions. The statistical analysis in R starts with individual bam files for each sample as inputs. To obtain these files, the package contains a python script (extdata/python/demultiplex.py) to demultiplex libraries and trim off primer sequences. With a standard alignment software the required bam files can be then be generated.

geecc - Use log-linear models to perform hypergeometric and chi-squared tests for gene set enrichments for two (based on contingency tables) or three categories (contingency cubes). Categories can be differentially expressed genes, GO terms, sequence length, GC content, chromosmal position, phylostrata, ....

GenomicInteractions - R package for handling Genomic interaction data, such as ChIA-PET/Hi-C, annotating genomic features with interaction information and producing various plots / statistics

GenomicTuples - GenomicTuples defines general purpose containers for storing genomic tuples. It aims to provide functionality for tuples of genomic co-ordinates that are analogous to those available for genomic ranges in the GenomicRanges Bioconductor package.

GenoView - Superimposing input data over existing genomic references allows for fast, accurate visual comparisons. The GenoView package is a novel bioinformatics package which condenses genomic data tracks to offer a comprehensive view of genetic variants. Its main function is to display mutation data over exons and protein domains, which easily identifies potential genomic locations of interest.

GOexpress - The package contains methods to visualise the expression levels of genes from a microarray or RNA-seq experiment and offers a clustering analysis to identify GO terms enriched in genes with expression levels best clustering predefined two or more groups of samples. Annotations for the genes present in the expression dataset are obtained from Ensembl through the biomaRt package. The random forest framework is used to evaluate the ability of each gene to cluster samples according to the factor of interest. Finally, GO terms are scored by averaging the rank (alternatively, score) of their respective gene sets to cluster the samples. An ANOVA approach is also available as an alternative statistical framework.

GOsummaries - A package to visualise Gene Ontology (GO) enrichment analysis results on gene lists arising from different analyses such clustering or PCA. The significant GO categories are visualised as word clouds that can be combined with different plots summarising the underlying data.

groHMM - A pipeline for the analysis of GRO-seq data.

GSAR - Gene set analysis using specific alternative hypotheses. Tests for differential expression, scale and net correlation structure.

GSReg - A package for gene set analysis based on the variability of expressions. It implements DIfferential RAnk Conservation (DIRAC) and gene set Expression Variation Analysis (EVA) methods.

HDTD - Characterization of intra-individual variability using physiologically relevant measurements provides important insights into fundamental biological questions ranging from cell type identity to tumor development. For each individual, the data measurements can be written as a matrix with the different subsamples of the individual recorded in the columns and the different phenotypic units recorded in the rows. Datasets of this type are called high-dimensional transposable data. The HDTD package provides functions for conducting statistical inference for the mean relationship between the row and column variables and for the covariance structure within and between the row and column variables.

hiAnnotator - hiAnnotator contains set of functions which allow users to annotate a GRanges object with custom set of annotations. The basic philosophy of this package is to take two GRanges objects (query & subject) with common set of seqnames (i.e. chromosomes) and return associated annotation per seqnames and rows from the query matching seqnames and rows from the subject (i.e. genes or cpg islands). The package comes with three types of annotation functions which calculates if a position from query is: within a feature, near a feature, or count features in defined window sizes. Moreover, each function is equipped with parallel backend to utilize the foreach package. In addition, the package is equipped with wrapper functions, which finds appropriate columns needed to make a GRanges object from a common data frame.

hiReadsProcessor - hiReadsProcessor contains set of functions which allow users to process LM-PCR products sequenced using any platform. Given an excel/txt file containing parameters for demultiplexing and sample metadata, the functions automate trimming of adaptors and identification of the genomic product. Genomic products are further processed for QC and abundance quantification.

IdeoViz - Plots data associated with arbitrary genomic intervals along chromosomal ideogram.

IMPCdata - Package contains methods for data retrieval from IMPC Database.

interactiveDisplayBase - The interactiveDisplayBase package contains the the basic methods needed to generate interactive Shiny based display methods for Bioconductor objects.

kebabs - The package provides functionality for kernel-based analysis of DNA, RNA, and amino acid sequences via SVM-based methods. As core functionality, kebabs implements following sequence kernels: spectrum kernel, mismatch kernel, gappy pair kernel, and motif kernel. Apart from an efficient implementation of standard position-independent functionality, the kernels are extended in a novel way to take the position of patterns into account for the similarity measure. Because of the flexibility of the kernel formulation, other kernels like the weighted degree kernel or the shifted weighted degree kernel with constant weighting of positions are included as special cases. An annotation-specific variant of the kernels uses annotation information placed along the sequence together with the patterns in the sequence. The package allows for the generation of a kernel matrix or an explicit feature representation in dense or sparse format for all available kernels which can be used with methods implemented in other R packages. With focus on SVM-based methods, kebabs provides a framework which simplifies the usage of existing SVM implementations in kernlab, e1071, and LiblineaR. Binary and multi-class classification as well as regression tasks can be used in a unified way without having to deal with the different functions, parameters, and formats of the selected SVM. As support for choosing hyperparameters, the package provides cross validation - including grouped cross validation, grid search and model selection functions. For easier biological interpretation of the results, the package computes feature weights for all SVMs and prediction profiles which show the contribution of individual sequence positions to the prediction result and indicate the relevance of sequence sections for the learning result and the underlying biological functions.

M3D - This package identifies statistically significantly differentially methylated regions of CpGs. It uses kernel methods (the Maximum Mean Discrepancy) to measure differences in methylation profiles, and relates these to inter-replicate changes, whilst accounting for variation in coverage profiles.

MAIT - The MAIT package contains functions to perform end-to-end statistical analysis of LC/MS Metabolomic Data. Special emphasis is put on peak annotation and in modular function design of the functions.

MBAmethyl - This package provides a function for reconstructing DNA methylation values from raw measurements. It iteratively implements the group fused lars to smooth related-by-location methylation values and the constrained least squares to remove probe affinity effect across multiple sequences.

MBASED - The package implements MBASED algorithm for detecting allele-specific gene expression from RNA count data, where allele counts at individual loci (SNVs) are integrated into a gene-specific measure of ASE, and utilizes simulations to appropriately assess the statistical significance of observed ASE.

MEIGOR - Global Optimization

metabomxtr - The functions in this package return optimized parameter estimates and log likelihoods for mixture models of truncated data with normal or lognormal distributions.

metagene - This package produces metagene plots to compare the behavior of DNA-interacting proteins at selected groups of genes/features. Pre-calculated features (such as transcription start sites of protein coding gene or enhancer) are available. Bam files are used to increase the resolution. Multiple combination of group of features and or group of bam files can be compared in a single analysis. Bootstraping analysis is used to compare the groups and locate regions with statistically different enrichment profiles.

MethylAid - A visual and interactive web application using RStudio's shiny package. Bad quality samples are detected using sample-dependent and sample-independent controls present on the array and user adjustable thresholds. In depth exploration of bad quality samples can be performed using several interactive diagnostic plots of the quality control probes present on the array. Furthermore, the impact of any batch effect provided by the user can be explored.

MethylMix - MethylMix is an algorithm implemented to identify hyper and hypomethylated genes for a disease. MethylMix is based on a beta mixture model to identify methylation states and compares them with the normal DNA methylation state. MethylMix uses a novel statistic, the Differential Methylation value or DM-value defined as the difference of a methylation state with the normal methylation state. Finally, matched gene expression data is used to identify, besides differential, functional methylation states by focusing on methylation changes that effect gene expression.

methylPipe - Memory efficient analysis of base resolution DNA methylation data in both the CpG and non-CpG sequence context. Integration of DNA methylation data derived from any methodology providing base- or low-resolution data.

MGFM - The package is designed to detect marker genes from Microarray gene expression data sets

miRNAtap - The package facilitates implementation of workflows requiring miRNA predictions, it allows to integrate ranked miRNA target predictions from multiple sources available online and aggregate them with various methods which improves quality of predictions above any of the single sources. Currently predictions are available for Homo sapiens, Mus musculus and Rattus norvegicus (the last one through homology translation).

missMethyl - Normalisation and testing for differential variability for data from Illumina's Infinium HumanMethylation450 array. The normalisation procedure is subset-quantile within-array normalisation (SWAN), which allows Infinium I and II type probes on a single array to be normalised together. The test for differential variability is based on an empirical Bayes version of Levene's test.

monocle - Monocle performs differential expression and time-series analysis for single-cell expression experiments. It orders individual cells according to progress through a biological process, without knowing ahead of time which genes define progress through that process. Monocle also performs differential expression analysis, clustering, visualization, and other useful tasks on single cell expression data.  It is designed to work with RNA-Seq and qPCR data, but could be used with other types as well.

MoPS - Identification and characterization of periodic fluctuations in time-series data.

MPFE - Estimate distribution of methylation patterns from a table of counts from a bisulphite sequencing experiment given a non-conversion rate and read error rate.

mQTL.NMR - mQTL.NMR provides a complete mQTL analysis pipeline for 1H NMR data. Distinctive features include normalisation using most-used approaches, peak alignment using RSPA approach, dimensionality reduction using SRV and binning approaches, and mQTL analysis for animal and human cohorts.

MSGFgui - This package makes it possible to perform analyses using the MSGFplus package in a GUI environment. Furthermore it enables the user to investigate the results using interactive plots, summary statistics and filtering. Lastly it exposes the current results to another R session so the user can seamlessly integrate the gui into other workflows.

MSGFplus - This package contains function to perform peptide identification using MS-GF+

MSnID - Extracts MS/MS ID data from mzIdentML (leveraging mzID package) or text files. After collating the search results from multiple datasets it assesses their identification quality and optimize filtering criteria to achieve the maximum number of identifications while not exceeding a specified false discovery rate. Also contains a number of utilities to explore the MS/MS results and assess missed and irregular enzymatic cleavages, mass measurement accuracy, etc.

MultiMed - Implements permutation method with joint correction for testing multiple mediators

mvGST - mvGST provides platform-independent tools to identify GO terms (gene sets) that are differentially active (up or down) in multiple contrasts of interest.  Given a matrix of one-sided p-values (rows for genes, columns for contrasts), mvGST uses meta-analytic methods to combine p-values for all genes annotated to each gene set, and then classify each gene set as being significantly more active (1), less active (-1), or not significantly differentially active (0) in each contrast of interest.  With multiple contrasts of interest, each gene set is assigned to a profile (across contrasts) of differential activity.  Tools are also provided for visualizing (in a GO graph) the gene sets classified to a given profile.

netbiov - A package that provides an effective visualization of large biological networks

NGScopy - NGScopy provides a quantitative caller for detecting copy number variations in next generation sequencing (NGS), including whole genome sequencing (WGS), whole exome sequencing (WES) and targeted panel sequencing (TPS). The caller can be parallelized by chromosomes to use multiple processors/cores on one computer.

OncoSimulR - Functions for simulating and plotting cancer progression data, including drivers and passengers, and allowing for order restrictions. Simulations use continuous-time models (based on Bozic et al., 2010 and McFarland et al., 2013) and fitness functions account for possible restrictions in the order of accumulation of mutations.

oposSOM - This package translates microarray expression data into metadata of reduced dimension. It provides various sample-centered and group-centered visualizations, sample similarity analyses and functional enrichment analyses. The underlying SOM algorithm combines feature clustering, multidimensional scaling and dimension reduction, along with strong visualization capabilities. It enables extraction and description of functional expression modules inherent in the data.

PAA - PAA imports single color (protein) microarray data that has been saved in gpr file format - esp. ProtoArray data. After pre-processing (background correction, batch filtering, normalization) univariate feature pre-selection is performed (e.g., using the "minimum M statistic" approach - hereinafter referred to as "mMs"). Subsequently, a multivariate feature selection is conducted to discover biomarker candidates. Therefore, either a frequency-based backwards elimination aproach or ensemble feature selection can be used. PAA provides a complete toolbox of analysis tools including several different plots for results examination and evaluation.

paxtoolsr - The package provides a basic set of R functions for interacting with BioPAX OWL files and the querying Pathway Commons (PC) molecular interaction data server, hosted by the Computational Biology Center at Memorial-Sloan-Kettering Cancer Center (MSKCC).

Pbase - A set of classes and functions to investigate and understand protein sequence data in the context of a proteomics experiment.

pepStat - Statistical analysis of peptide microarrays

polyester - This package can be used to simulate RNA-seq reads from differential expression experiments with replicates. The reads can then be aligned and used to perform comparisons of methods for differential expression.

Polyfit - Polyfit is an add-on to the packages DESeq which ensures the p-value distribution is uniform over the interval [0, 1] for data satisfying the null hypothesis of no differential expression, and uses an adpated Storey-Tibshiran method to calculate q-values.

proBAMr - Mapping PSMs back to genome. The package builds SAM file from shotgun proteomics data The package also provides function to prepare annotation from GTF file.

pRolocGUI - The package pRolocGUI comprises functions to interactively visualise organelle (spatial) proteomics data on the basis of pRoloc, pRolocdata and shiny.

proteoQC - This package creates a HTML format QC report for MS/MS-based proteomics data. The report is intended to allow the user to quickly assess the quality of proteomics data.

PSEA - Deconvolution of gene expression data by Population-Specific Expression Analysis (PSEA).

Pviz - Pviz adapts the Gviz package for protein sequences and data.

quantro - A data-driven test for the assumptions of quantile normalization using raw data such as objects that inherit eSets (e.g. ExpressionSet, MethylSet). Group level information about each sample (such as Tumor / Normal status) must also be provided because the test assesses if there are global differences in the distributions between the user-defined groups.

rain - This package uses non-parametric methods to detect rhythms in time series. It deals with outliers, missing values and is optimized for time series comprising 10-100 measurements. As it does not assume expect any distinct waveform it is optimal or detecting oscillating behavior (e.g. circadian or cell cycle) in e.g. genome- or proteome-wide biological measurements such as: micro arrays, proteome mass spectrometry, or metabolome measurements.

regionReport - Generate HTML reports to explore a set of regions such as the results from annotation-agnostic expression analysis of RNA-seq data at base-pair resolution performed by derfinder.

RGSEA - Combining bootstrap aggregating and Gene set enrichment analysis (GSEA), RGSEA is a classfication algorithm with high robustness and no over-fitting problem. It performs well especially for the data generated from different exprements.

riboSeqR - Plotting functions, frameshift detection and parsing of sequencing data from ribosome profiling experiments.

Rnits - R/Bioconductor package for normalization, curve registration and inference in time course gene expression data

Rqc - Rqc is an optimised tool designed for quality control and assessment of high-throughput sequencing data. It performs parallel processing of entire files and produces a report which contains a set of high-resolution graphics.

rRDP - Seamlessly interfaces RDP classifier (version 2.9).

RUVnormalize - RUVnormalize is meant to remove unwanted variation from gene expression data when the factor of interest is not defined, e.g., to clean up a dataset for general use or to do any kind of unsupervised analysis.

RUVSeq - This package implements the remove unwanted variation (RUV) methods of Risso et al. (2014) for the normalization of RNA-Seq read counts between samples.

S4Vectors - The S4Vectors package defines the Vector and List virtual classes and a set of generic functions that extend the semantic of ordinary vectors and lists in R. Package developers can easily implement vector-like or list-like objects as concrete subclasses of Vector or List. In addition, a few low-level concrete subclasses of general interest (e.g. DataFrame, Rle, and Hits) are implemented in the S4Vectors package itself (many more are implemented in the IRanges package and in other Bioconductor infrastructure packages).

SemDist - This package implements methods to calculate information accretion for a given version of the gene ontology and uses this data to calculate remaining uncertainty, misinformation, and semantic similarity for given sets of predicted annotations and true annotations from a protein function predictor.

seqplots - SeqPlots is a tool for plotting next generation sequencing (NGS) based experiments' signal tracks, e.g. reads coverage from ChIP-seq, RNA-seq and DNA accessibility assays like DNase-seq and MNase-seq, over user specified genomic features, e.g. promoters, gene bodies, etc. It can also calculate sequence motif density profiles from reference genome. The data are visualized as average signal profile plot, with error estimates (standard error and 95% confidence interval) shown as fields, or as series of heatmaps that can be sorted and clustered using hierarchical clustering, k-means algorithm and self organising maps. Plots can be prepared using R programming language or web browser based graphical user interface (GUI) implemented using Shiny framework. The dual-purpose implementation allows running the software locally on desktop or deploying it on server. SeqPlots is useful for both for exploratory data analyses and preparing replicable, publication quality plots. Other features of the software include collaboration and data sharing capabilities, as well as ability to store pre-calculated result matrixes, that combine many sequencing experiments and in-silico generated tracks with multiple different features. These binaries can be further used to generate new combination plots on fly, run automated batch operations or share with colleagues, who can adjust their plotting parameters without loading actual tracks and recalculating numeric values. SeqPlots relays on Bioconductor packages, mainly on rtracklayer for data input and BSgenome packages for reference genome sequence and annotations.

SGSeq - RNA-seq data are informative for the analysis of known and novel transcript isoforms. While the short length of RNA-seq reads limits the ability to predict and quantify full-length transcripts, short read data are well suited for the analysis of individual alternative transcripts events (e.g. inclusion or skipping of a cassette exon). The SGSeq package enables the prediction, quantification and visualization of alternative transcript events from BAM files.

shinyMethyl - Interactive tool for visualizing Illumina's 450k array data

SigCheck - While gene signatures are frequently used to classify data (e.g. predict prognosis of cancer patients), it it not always clear how optimal or meaningful they are (cf David Venet, Jacques E. Dumont, and Vincent Detours' paper "Most Random Gene Expression Signatures Are Significantly Associated with Breast Cancer Outcome"). Based partly on suggestions in that paper, SigCheck accepts a data set (as an ExpressionSet) and a gene signature, and compares its classification performance (using the MLInterfaces package) against a) random gene signatures of the same length; b) known, (related and unrelated) gene signatures; and c) permuted data.

simulatorZ - simulatorZ is a package intended primarily to simulate collections of independent genomic data sets, as well as performing training and validation with predicting algorithms. It supports ExpressionSets and SummarizedExperiment objects.

SNPRelate - Genome-wide association studies (GWAS) are widely used to investigate the genetic basis of diseases and traits, but they pose many computational challenges. We developed an R package SNPRelate to provide a binary format for single-nucleotide polymorphism (SNP) data in GWAS utilizing CoreArray Genomic Data Structure (GDS) data files. The GDS format offers the efficient operations specifically designed for integers with two bits, since a SNP could occupy only two bits. SNPRelate is also designed to accelerate two key computations on SNP data using parallel computing for multi-core symmetric multiprocessing computer architectures: Principal Component Analysis (PCA) and relatedness analysis using Identity-By-Descent measures. The SNP format in this package is also being used by the GWASTools package with the support of S4 classes and generic functions.

specL - specL provides a function for generating spectra libraries which can be used for MRM SRM MS workflows in proteomics. The package provides a BiblioSpec reader, a function which can add the protein information using a FASTA formatted amino acid file, and an export method for using the created library in the Spectronaut software.

ssviz - Small RNA sequencing viewer

STAN - STAN (STrand-specic ANnotation of genomic data) implements bidirectional Hidden Markov Models (bdHMM), which are designed for studying directed genomic processes, such as gene transcription, DNA replication, recombination or DNA repair by integrating genomic data. bdHMMs model a sequence of successive observations (e.g. ChIP or RNA measurements along the genome) by a discrete number of 'directed genomic states', which e.g. reflect distinct genome-associated complexes. Unlike standard HMM approaches, bdHMMs allow the integration of strand-specific (e.g. RNA) and non strand-specific data (e.g. ChIP).

STATegRa - Classes and tools for multi-omics data integration.

switchBox - The package offer different classifiers based on comparisons of pair of features (TSP), using various decision rules (e.g., majority wins principle).

systemPipeR - R package for building end-to-end analysis pipelines with automated report generation for next generation sequence (NGS) applications such as RNA-Seq, ChIP-Seq, VAR-Seq and many others. An important feature is support for running command-line software, such as NGS aligners, on both single machines or compute clusters. This includes both interactive job submissions or batch submissions to queuing systems of clusters.

ToPASeq - Implementation of seven methods for topology-based pathway analysis of both RNASeq and microarray data: SPIA, DEGraph, TopologyGSA, TAPPA, TBS, PWEA and a visualization tool for a single pathway.

tracktables - Methods to create complex IGV genome browser sessions and dynamic IGV reports in HTML pages.

TSCAN - TSCAN enables users to easily construct and tune pseudotemporal cell ordering as well as analyzing differentially expressed genes. TSCAN comes with a user-friendly GUI written in shiny. More features will come in the future.

wavClusteR - Infer PAR-CLIP induced transitions and discriminate them from sequencing error, SNPs, contaminants and additional non-experimental causes, using a non-parametric mixture model. wavClusteR resolves cluster boundaries at high resolution and provides robust estimation of cluster statistics. In addition, the package allows to integrate RNA-Seq data to estimate FDR over the entire range of relative substitution frequencies. Furthermore, the package provides post-processing of results and functions to export results for UCSC genome browser visualization and motif search analysis. Key functions support parallel multicore computing. While wavClusteR was designed for PAR-CLIP data analysis, it can be applied to the analysis of other Next-Generation Sequencing data obtained from substitution inducing experimental procedures (e.g. BisSeq)



NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages. The following package NEWS is available:




Packages removed from the release
=================================

The following packages are no longer in the release:

Clonality, Resourcerer, RMAPPER, virtualArray


