Bioconductors:

We are pleased to announce Bioconductor 3.0, consisting of 934
software packages, 219 experiment data packages, and 870
up-to-date annotation packages. 

There are 114 new software packages, and many updates and improvements
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

There are 114 new packages in this release of Bioconductor.

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

EBSeqHMM - The EBSeqHMM package implements an auto-regressive hidden
Markov model for statistical analysis in ordered RNA-seq experiments
(e.g. time course or spatial course data). The EBSeqHMM package
provides functions to identify genes and isoforms that have non-
constant expression profile over the time points/positions, and
cluster them into expression paths.

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

Metab - Metab is an R package for high-throughput processing of
metabolomics data analysed by the Automated Mass Spectral
Deconvolution and Identification System (AMDIS)
(http://chemdata.nist.gov/mass-spc/amdis/downloads/). In addition, it
performs statistical hypothesis test (t-test) and analysis of variance
(ANOVA). Doing so, Metab considerably speed up the data mining process
in metabolomics and produces better quality results. Metab was
developed using interactive features, allowing users with lack of R
knowledge to appreciate its functionalities.

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

mygene - MyGene.Info_ provides simple-to-use REST web services   to
query/retrieve gene annotation data. It's designed   with simplicity
and performance emphasized. *mygene*,   is an easy-to-use R wrapper to
access MyGene.Info_ services.

netbiov - A package that provides an effective visualization of large biological networks

NGScopy - NGScopy provides a quantitative caller for detecting copy number variations in next generation sequencing (NGS), including whole genome sequencing (WGS), whole exome sequencing (WES) and targeted panel sequencing (TPS). The caller can be parallelized by chromosomes to use multiple processors/cores on one computer.

OncoSimulR - Functions for simulating and plotting cancer progression data, including drivers and passengers, and allowing for order restrictions. Simulations use continuous-time models (based on Bozic et al., 2010 and McFarland et al., 2013) and fitness functions account for possible restrictions in the order of accumulation of mutations.

oposSOM - This package translates microarray expression data into metadata of reduced dimension. It provides various sample-centered and group-centered visualizations, sample similarity analyses and functional enrichment analyses. The underlying SOM algorithm combines feature clustering, multidimensional scaling and dimension reduction, along with strong visualization capabilities. It enables extraction and description of functional expression modules inherent in the data.

PAA - PAA imports single color (protein) microarray data that has been saved in gpr file format - esp. ProtoArray data. After pre-processing (background correction, batch filtering, normalization) univariate feature pre-selection is performed (e.g., using the "minimum M statistic" approach - hereinafter referred to as "mMs"). Subsequently, a multivariate feature selection is conducted to discover biomarker candidates. Therefore, either a frequency-based backwards elimination aproach or ensemble feature selection can be used. PAA provides a complete toolbox of analysis tools including several different plots for results examination and evaluation.

paxtoolsr - The package provides a basic set of R functions for interacting with BioPAX OWL files and the querying Pathway Commons (PC) molecular interaction data server, hosted by the Computational Biology Center at Memorial-Sloan-Kettering Cancer Center (MSKCC).

Pbase - A set of classes and functions to investigate and understand protein sequence data in the context of a proteomics experiment.

pepStat - Statistical analysis of peptide microarrays

pepXMLTab - Parsing pepXML files based one XML package.  The package
tries to handle pepXML files generated from different softwares. The
output will be a peptide-spectrum-matching tabular file. The package
also provide function to filter the PSMs based on FDR.

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

seqTools - Analyze read length, phred scores and alphabet frequency
and DNA k-mers on uncompressed and compressed fastq files.

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


ADaCGH2
-------

Changes in version 2.5.2 (2014-10-03):

- Corrected help for "inputToADaCGH"

Changes in version 2.5.1 (2014-10-03):

- Fixed bug reading files when R object was called "inputData"

Changes in version 2.5.0:

- Bumped BioC version.

affxparser
----------

Changes in version 1.37.2 (2014-09-28):

- Minor modifications due to the move to GitHub.

Changes in version 1.37.1 (2014-08-25):

- CLEANUP: Removed R CMD check NOTEs that appeared in recent R versions.

Changes in version 1.37.0 (2014-04-11):

- The version number was bumped for the Bioconductor devel version, which is now BioC v2.15 for R (>= 3.1.0).

ALDEx2
------

Changes in version 0.99.2:

NEW FEATURES

- made aldex.clr into a class

- allowed input of SummarizedExperiment object instead of a reads data frame

- prioritized use of the BiocParallel package for multicore processing. If BiocParallel is not installed then the parallel
  package used, if neither packages are installed, then serial processing is used

Changes in version 0.99.1:

NEW FEATURES

- changed conditional tests for multiprocessor use, defaults to FALSE

Changes in version 0.99.0:

NEW FEATURES

- first submission to Bioc-devel

annotate
--------

Changes in version 1.43:

NEW FEATURES

- blastSequences accepts an argument timeout limiting waiting time for a response; in an interactive session and after the
  timeout is reached, the user may opt to retry the query.

- blastSequences accepts an argument as controlling the representation of the return value, either a DNAMultipleAlignment, a
  data.frame, or the XML.

AnnotationHub
-------------

Changes in version 1.6.0:

NEW FEATURES

- add query and subset,AnnotationHub=methods

aroma.light
-----------

Changes in version 2.1.2 (2014-09-23):

- Minor tweaks due to the move to GitHub.

Changes in version 2.1.1 (2014-09-16):

- IMPORTANT/CLEANUP: The matrixStats package is no longer attached with this package.  In other words, you now might have to
  add library('matrixStats') to your scripts.

- CLEANUP: Now importing R.utils (instead of only suggesting it).

- Fixed some new R CMD check NOTEs.

Changes in version 2.1.0 (2014-04-11):

- The version number was bumped for the Bioconductor devel version, which is now BioC v2.15 for R (>= 3.1.0).

ASGSCA
------

Changes in version 0.99.2:

- Minor bioconductor submission issues resolved

Changes in version 0.99.1:

- Added real data in examples

- Changed biocViews tag

Changes in version 0.99.0:

- Initial release to Bioconductor.

- Added NEWS file.

ballgown
--------

9.7: new methods for accessing gene and transcript names and IDs

9.7: development package (on GitHub) now has Travis CI integration, so the package's is built/tested with every new push
             and a status image is visible on the package README

9.6: ballgownrsem can now handle gzipped input files

9.6: small fixes to eliminate warnings/errors on R CMD CHECK

9.5: "subset" function now behaves appropriately if sampleNames contain "dot" characters

9.4: Tablemaker source code moved out of Ballgown repository (both GitHub and BioC svn repos)

9.3: `exprfilter` function added

9.2: `ballgownrsem` function added (now compatible with RSEM output)

9.2: RSEM slot added to ballgown S4 class

9.2: bug fixes in `subset` function

9.2: library size adjustment in `stattest` is now CSS normalization (using *log* FPKM) by default and is more
             customizable

0.99.1:

9.0: "meas" argument added to constructor function

9.0: unit tests added

9.0: vignette updated

9.0: several bug fixes in plotting functions and statistical testing function Alpha release, 30 March 2014:

9.0: package announced


BiocParallel
------------

CHANGES IN VERSION 1.0.0
------------------------

NEW FEATURES

    o Add vignette sections for cluster managers, AMI

    o Add bpiterate generic and methods 

    o Add REDUCE to bpiterate()

    o Add 'reduce.in.order' to bpiterate()


MODIFICATIONS

    o Update vignette examples, reorganize sections 

    o Allow 'workers' in BiocParallelParam to be character or integer

    o Enhance bpresume() man page; add examples

    o Enhance register() man page; add examples

    o Improve default registration for SnowParam:
      - max 8 cores
      - use detectcores() / mc.cores if available

    o Modify .convertToSimpleError() to convert NULL to NA_character_


BUG FIXES

    o Fix recursion problem for BPPARAM as list 

    o Modify bpaggregate() to run in parallel


BiocStyle
---------

Changes in version 1.4.0:

USER VISIBLE CHANGES

- Support for markdown documents

- Add \Githubpkg markup command

biosvd
------

Changes in version 2.0.4:

Category Class

- Class EigensystemPlotParam has been added, to allow all plot-related settings to be specified in an object of this class.

Category Functionality

- The plots function now provides up to ten visualizations. This function generates a heatmap of the eigenfeatures by assays
  with use of the given contrast factor (eigenfeatureHeatmap), a heatmap of the features by eigenassays with use of the given
  contrast factor (eigenassayHeatmap), a heatmap of the features by assays with features sorted according to their correlation
  with two eigenfeatures (sortedHeatmap), a bar plot with the eigenexpression fractions of all eigenfeatures (fraction), a
  screeplot for the eigenexpression fractions (scree), a bar plot with the eigenexpression fractions of the eigenfeatures
  without the dominant eigenfeature(s) (zoomedFraction), the intensity levels of selected eigenfeatures across the assays (by
  default eigenfeatures 1-4) (lines), the intensity levels of all eigenfeatures across the assays (allLines), polar plot for
  the assays according to their correlation with two eigenassays (eigenassayPolar), and polar plot for the features according
  to their correlation with two eigenfeatures (eigenfeaturePolar).

Category Dependencies

- The package now longer depends on ReportingTools.

Category Deprecation

- The sort function has been deprecated, and replaced by the function project. This function returns the rectangular and polar
  coordinates of the eigensystem projection onto two specified eigenfeatures and eigenassays.

BitSeq
------

Changes in version 1.10.0 (2014-10-01):

NEW FEATURES

- This release corresponds to C++ BitSeq version 0.7.5.

- adding "unstranded" option to parseAlignment, to allow read pairs with various directions to be used

- enable excluding singletons (single-mate alignments of paired reads) in parseAlignment

BridgeDbR
---------

Changes in version 0.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- No changes

NEW FEATURES

- New package

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)


CGEN
----

Changes in version 3.0:

USER VISIBLE CHANGES

- The functions score.test, GxE.scan, GxE.scan.partition and GxE.scan.combine have been added.

- The functions snp.logistic and snp.scan.logistic can now handle imputed genotypes.

BUG FIXES

- The function snp.logistic checks that the snp variable is used in main.vars or int.vars.

PLANS

- The next release will have a new Wald test function.

chimera
-------

Changes in version 1.7.10:

NEW FEATURES

- Fusions data generated by fusionCatcher can be imported.

Changes in version 1.7.1:

NEW FEATURES

- Fusions data generated by Rsubread can be imported.

chipenrich
----------

Changes in version 1.3.4:

PKG FEATURES

- A new method, broadenrich, is available in the chipenrich function which is designed for gene set enrichment on broad
  genomic regions, such as peaks resulting from histone modificaiton based ChIP-seq experiments.

- Methods chipenrich and broadenrich are available in multicore versions (on every platform except Windows). The user selects
  the number of cores when calling the chipenrich function.

- Peaks downloaded from the ENCODE Consortium as .broadPeak or .narrowPeak files are supported directly.

- Peaks downloaded from the modENCODE Consortium as .bed.gff or .bed.gff3 files are also supported directly.

- Support for D. melanogaster (dm3) genome and enrichment testing for GO terms from all three branches (GOBP, GOCC, and GOMF).

- New gene sets from Reactome (http://www.reactome.org) for human, mouse, and rat.

- New example histone data set, peaks_H3K4me3_GM12878, based on hg19.

ChIPQC
------

Changes in version 1.1.2:

- Big fixes/man page fixes

ChIPseeker
----------

Changes in version 1.1.21:

- use data.table instead of data.frame to optimize covplot <2014-10-06, Mon>

Changes in version 1.1.20:

- update annotatePeak to store the seqinfo information <2014-09-30, Tue>

- modified runValue(x) to sapply(x, runValue) <2014-09-30, Tue>

Changes in version 1.1.19:

- implement csAnno S4 object <2014-09-28, Sun>

- modify plot function for csAnno instance <2014-09-28, Sun>

- implement vennpie function <2014-09-28, Sun>

Changes in version 1.1.17:

- deprecate plotChrCov to new function covplot <2014-08-18, Mon>

- add new paramter chrs and xlim to covplot <2014-08-18, Mon>

Changes in version 1.1.16:

- optimize plotChrCov, running time reduce drastically <2014-08-15, Fri>

Changes in version 1.1.15:

- remove un-mappable peak to prevent fail in peak annotation <2014-08-14, Thu>

Changes in version 1.1.14:

- bug fixed in plotDistToTSS <2014-08-14, Thu>

Changes in version 1.1.13:

- change TranscriptDb to TxDb according to GenomicFeatures <2014-07-29, Tue>

Changes in version 1.1.12:

- bug fixed in plotChrCov <2014-07-21, Mon>

Changes in version 1.1.10:

- bug fixed in calculating distances from peak end <2014-06-18, Wed>

Changes in version 1.1.9:

- add level parameter to annotatePeak, and set it to "transcript" by default.  Now annotatePeak will annotate peaks in
  transcript level except user specify level = "gene" <2014-06-16, Mon>

- add addFlankGeneInfo parameter to annotatePeak.  If it set to true, all features within the flankDistance will be annotated.
  <2014-06-16, Mon>

Changes in version 1.1.8:

- bug fixed when peak overlap with feature <2014-06-11, Wed>

- optimize for getting overlap features of peaks <2014-06-11, Wed>

- update plotAnnoPie, separate the pie and legend to prevent label overlap <2014-06-12, Thu>

Changes in version 1.1.7:

- bug fixed in calculating distanceToTSS <2014-06-03, Tue>

Changes in version 1.1.6:

- add chainFile parameter in enrichAnnoOverlap and enrichPeakOverlap to support different genome version comparision
  <2014-06-01, Sun>

- fixed color bug in peakHeatmap.internal2 and plotAnnoBar <2014-06-02, Mon>

- update vignettes <2014-06-02, Mon>

Changes in version 1.1.5:

- export getPromoters and getTagMatrix <2014-05-31, Sat>

- rename plotAvgProf to plotAvgProf2 and implement plotAvgProf based on tagMatrix <2014-05-31, Sat>

- implement tagHeatmap for visualize heatmap of the tagMatrix or a list of tagMatrix <2014-05-31, Sat>

- implement shuffle function to generate a random ChIP data based on a real one <2014-05-31, Sat>

- implement enrichPeakOverlap to calcuate significant of ChIP experiments based on the genome coordinations <2014-05-31, Sat>

- implement enrichAnnoOverlap to calculate significant of ChIP experiments based on their nearest gene annotation <2014-05-31,
  Sat>

- incorporate GEO database for mining significant overlap of ChIP data <2014-05-31, Sat> + getGEOspecies summarize the
  collected data by species + getGEOgenomeVersion summarize the colleted data by genome version + getGEOInfo extract the
  information by genome version query + downloadGEObedFiles download all bed files of a particular genome version +
  downloadGSMbedFiles download the bed files of the input GSM list.

Changes in version 1.1.4:

- in the annotation column of output of annotatePeak function, if Exon/Intron, the output change to 'Transcript_Name/GeneID,
  Exon no. of total_no.' <2014-05-14, Wed>

Changes in version 1.1.3:

- bug fixed when metadata(TranscriptDb) contained NA <2014-04-30, Wed>

- support ID type of Ensembl in annotatePeak (Entrez was supported) <2014-04-30, Wed>

Changes in version 1.1.2:

- implemented plotChrCov <2014-04-25, Fri>

- implemented plotAvgProf and peakHeatmap <2014-04-24, Thu>

Changes in version 1.1.1:

- output of annotatePeak now contain chromosome length information <2014-04-22, Tue>

- re-implement plotAnnoPie to use ordinary pie plot instead of pie3D <2014-04-21, Mon>

cisPath
-------

Changes in version 1.5.8:

- documentation improvements

Changes in version 1.5.7:

- Fix a bug

Changes in version 1.5.6:

- Fix a bug

Changes in version 1.5.5:

- improvements

Changes in version 1.5.4:

- improvements

Changes in version 1.5.3:

- A method has been added to format PPI data (SIF format) downloaded from PINA2.

Changes in version 1.5.2:

- Fix a bug

Changes in version 1.5.1:

- Fix a bug

cleaver
-------

Changes in version 1.3.8 (2014-09-28):

- Use title case in the "Title:" field of the DESCRIPTION file.

Changes in version 1.3.7 (2014-05-08):

- Create an IRangesList using `split` instead of `lapply` in the cleavageRanges,AAStringSet-method is more than two times
  faster.

Changes in version 1.3.6 (2014-05-03):

- Add cleavageRanges method for character, AAString and AAStringSet.

- cleave,AAString returns an AAStringSet instead of an AAStringSetList object.

- Fix return value of cleavageSites,AAStringSet.

Changes in version 1.3.5 (2014-04-30):

- Fix defintion of cleavageSites,AAStringSet.

Changes in version 1.3.4 (2014-04-30):

- Add cleavageSites method for character, AAString and AAStringSet.

Changes in version 1.3.3 (2014-04-28):

- Avoid duplicated digest of peptides results in a hugh speed improvement and a hugh memory reduction (removes fix from 1.1.8
  and partly reintroduces original algorithm).

- Remove memory test and "memoryThreshold" argument (fails on different platforms and is not important anymore using the "new"
  cleavage algorithm).

- Change default of "unique" argument to "unique=TRUE".

Changes in version 1.3.2 (2014-04-25):

- Introduce simple memory test and "memoryThreshold" argument to avoid the calculation of ridiculous high numbers of
  "missedCleavages" and peptide combinations.

Changes in version 1.3.1 (2014-04-25):

- Add "custom" argument to allow the user defining own cleavage rules.

Clomial
-------

Changes in version 1.1.9 (2014-07-26):

- Adding QC for the situation when only 1 sample is provided.

Changes in version 1.1.7 (2014-04-24):

- The function Clomial.generate.data() can now model sequencing error.

Changes in version 1.1.2 (2014-04-18):

- Some bugs fixed regarding ignoring samples.

- Some bugs fixed regarding loading saved results when running in parallel.

- Some typos fixed in the man files.

- More efficient job submissions.

- The function compute.errors() is now exported, and can be accessed by the user.

Changes in version 1.1.0 (2014-03-11):

- Approved by Bioconductor.

clonotypeR
----------

Changes in version 1.4.0:

BUG FIXES

- Minor packaging corrections (NAMESPACE, etc.).

clusterProfiler
---------------

Changes in version 1.99.0:

- set version to 1.99.0 and will release version 2.0.0 in next BioConductor release. <2014-09-03, Wed> + In version 2.0.0, the
  package was extended to support about 20 species for both of GO and KEGG analyses. + support gene set enrichment analysis
  algorithm for both of GO and KEGG. + support enrichMap visualization.

- support about 20 species of KEGG analyses. <2014-09-03, Wed>

- update man file to indicate that enrichGO is now support more than 20 species <2014-09-03, Wed>

Changes in version 1.13.3:

- add keepGOlevel, keepGOterm, excludeGOlevel, excludeGOterm, not exported yet <2014-08-28, Thu>

Changes in version 1.13.2:

- implement gseGO and gseKEGG for gene set enrichment analysis <2014-07-31, Thu>

- import enrichMap from DOSE <2014-07-31, Thu>

- add corresponding section in vignettes <2014-07-31, Thu>

Changes in version 1.13.1:

- update man files according to the change of roxygen2 (ver 4.0.0) <2014-05-16, Fri>

compEpiTools
------------

Changes in version 1.0.0:

- initial version with the following functions implemented: + GRbaseCoverage + GRcoverage + GRcoverageInbins +
  GRcoverageSummit + GRenrichment + countOverlapsInBins + stallingIndex + TSS + distanceFromTSS + GRangesInPromoters +
  GRmidpoint + GRannotate + GRannotateSimple + makeGtfFromDb + enhancers + matchEnhancers + topGOres + simplifyGOterms +
  findLncRNA + getPromoterClass + heatmapData + palette2d + heatmapPlot + plotStallingIndex + GR2fasta + overlapOfGRanges +
  GRsetwidth + unionMaxScore + GRanges2ucsc + ucsc2GRanges

CompGO
------

Changes in version 1.1.3:

- Added ability to plot absolute values of Z-scores in plotZScores()

- Fixed vignette to include Abs() plotting

Changes in version 1.1.2:

- Included static vignette for illustration purposes. Some steps in our pipeline require either writing to disk or the use of
  a registered email, which are difficult to include in a dynamically-compiled vignette. Using a vignette that is able to show
  each of the steps in the pipeline with a helpful degree of detail is better than having a dynamic but less expressive one.

- The vignette included uses the example datasets from He (2011) to run through the entire pipeline, illustrating the major
  steps and providing a useful template for analysis.

Changes in version 1.1.1:

- Added functionality to plot differentially-enriched KEGG pathways using software package pathview

- Added interactive plotting function to generate comparisons of many lists at once using PCA and dendrogram plots

- Added several example datasets of .bed files from He, A. (2011).


cosmiq
------

Changes in version 0.99.1:

USER VISIBLE CHANGES

- reduced rt window to decreasse the overall vignette compilation time

COSNet
------

Changes in version 0.99.5:

- Modified the Section 1 of vignettes. USER VISIBLE CHANGES BUG FIXES PLANS

CRISPRseek
----------

Changes in version 1.3.10:

NEW FEATURE

- gRNAs are automatically output to a bed file for view in the UCSC genome browser

Changes in version 1.3.8:

NEW FEATURE

- Added calculategRNAefficiency to calculate gRNA cleavage efficiency using DoenchNBT2014 predictive logistic model

Changes in version 1.3.7:

NEW FEATURE

- Added searchDirection to compare2Sequences to allow search one against the other and many to many sequence search

- Added exception handling to catch no gRNA found error in compare2Sequences

Changes in version 1.3.5:

BUG FIX

- TopN offtarget total score was sometimes missing for sequences containing gRNAs with less than 6 offtargets in version 1.3.3

Changes in version 1.3.3:

NEW FEATURE

- Search for off-targets is much faster when more than 10 gRNAs are searched

- Added new optional parameter orgAnn in offTargetAnalysis

- Added gene ID and optional gene symbol in off-target output file

- Added gRNA target region, GC content of gRNA and number of Ts in the last 4 postion of gRNA (not including PAM sequence) in
  the summary output file

csaw
----

Changes in version 1.0.0:

- New package csaw, for de novo detection of differential binding in ChIP-seq data.

cummeRbund
----------

2.7.1: Bugfixes: - Fixed 'fullnames' argument to cuffData::*Matrix() methods so that it does what it's supposed to do.
             - Added 'showPool' argument to fpkmSCVPlot.  When TRUE, empirical mean and standard deviation are
             determined across all conditions as opposed to cross-replicate. This is set to TRUE anytime you have n<2
             replicates per condition.  - Added stat="identity" to expressionBarplot to comply with ggplot 0.9.3
             enforcement.  - 'labels' argument to csScatter is now working as it's supposed to.  You can pass a vector
             of 'gene_short_name' identifiers to labels and these will be specifically called out in red text on
             scatterplot.  - Added repFpkmMatrix() and replicates() methods to CuffFeature objects.  - Removed
             unnecessary Joins to optimize retrieval speed for several key queries.  - Fixed bug in csVolcano matrix
             that forced ylimits to be c(0,15) New Features: - Added csNMF() method for CuffData and CuffFeatureSet
             objects to perform non-negative matrix factorization.  As of now, it's merely a wrapper around the default
             settings for NMFN::nnmf(), but hope to expand in the future.  * Does not adjust sparsity of matrices after
             output, must be done by user as needed.  - Added csPie() method for CuffGene objects. Allows for
             visualization of relative isoform, CDS, and promoter usage proportions as a pie chart by condition (or
             optionally as stacked bar charts by adding + coord_cartesian() ).  - Added 'method' argument to csCluster
             and csHeatmap to allow custom distance functions for clustering. Default = "none" = JSdist(). You can now
             provide a function that returns a 'dist' object on rows of a matrix.  - Added varModel.info tracking for
             compatibility with cuffdiff >=2.1. Will now find varModel.info file if exists, and incorporate into
             database.  - dispersionPlot() method added for CuffSet object.  This now appropriately draws from
             varModel.info and is the preferred visualization for dispersion of RNA-Seq data with cummeRbund.  - Added
             diffTable() method to CuffData and CuffFeatureSet objects to allow a 'one-table' snapshot of results for
             all Features (CuffData) or a set of Features (CuffFeatureSet). This table outputs key values including gene
             name, gene short name, expression estimates and per-comparison fold-change, p-value, q-value, and
             significance values (yes/no). A convenient 'data-dump' function to merge across several tables.  - Added
             coercion methods for CuffGene objects to create GRanges and GRangeslist objects (more BioC friendly!). Will
             work on making this possible on CuffFeatureSet and CuffFeature objects as well.  - Added pass-through to
             select p.adjust method for getSig (method argument to getSig) - Added ability to revert to cuffdiff
             q-values for specific paired-wise interrogations with getSig as opposed to re-calculating new ones
             (useCuffMTC; default=FALSE) Notes: - Removed generic for 'featureNames'.  Now appropriately uses
             featureNames generic from Biobase.  As a consequence, Biobase is now a dependency.  - Added passthrough to
             as.dist(...) in JSdist(...)  - Added 'logMode' argument to csClusterPlot.  - Added 'showPoints' argument to
             PCAplot to allow disabling of gene values in PCA plot. If false, only sample projections are plotted.  -
             Added 'facet' argument to expressionPlot to disable faceting by feature_id.  - shannon.entropy now uses
             log2 instead of log10 to constrain specificity scores between 0 and 1.

customProDB
-----------

Changes in version 1.5.5:

BUG FIXES

- Fix a bug in function OutputNovelJun.R

Changes in version 1.5.4:

BUG FIXES

- Fix a bug, novel junction coding sequence on '-' strand, in function OutputNovelJun.R

Changes in version 1.5.3:

BUG FIXES

- Fix a bug, ID mapping, in function OutputNovelJun.R

dagLogo
-------

Changes in version 1.3.6:

NEW FEATURES

- add test files for fetchSequence function

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.3.4:

NEW FEATURES

- No changes classified as 'bug fixes' (package under active development)

BUG FIXES

- fetchSequences(): when using proteome object created from AAStringSet, the columns 2 and 3 in the sequences data frame do
  not contain the expected values.

- testDAU(): the counts() function returned NA for letters that were missing in a given matrix; that lead to NaNs in the
  background matrix and to errors during dagLogo() plotting. changed counts() to return 0 for missing letters.

- fix the typo in dagLogo()

Changes in version 1.3.3:

NEW FEATURES

- add background Noise into test

- column label relative to the anchor when draw logos

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

deepSNV
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

Changes in version 1.11.1:

Bugfixes

- Report correct table in summary(..., value="VCF"); Thanks to David Gacquer.

DEGreport
---------

0.99.12: 09-04-2014 Lorena Pantano <lorena.pantano@gmail.com> CORRECT DOC OF createReport

0.99.11: 09-02-2014 Lorena Pantano <lorena.pantano@gmail.com> ADD ncores TO REPORT CREATION

0.99.10: 08-16-2014 Lorena Pantano <lorena.pantano@gmail.com> REMOVE BIOMART FROM VIGNETTES DUE TO WINDOWS UNKNOWN BUILD
             ISSUES

0.99.9: 08-02-2014 Lorena Pantano <lorena.pantano@gmail.com> ADDING PARALLELIZATION TO BAYESIAN INFERENCE

0.99.4: 05-19-2014 Lorena Pantano <lorena.pantano@gmail.com> CODING STYLE * replacing tab by 4 spaces * cleaning up
             function * adding unit tests

derfinder
---------

Changes in version 0.99.0:

NEW FEATURES

- Preparing to submit to Bioconductor.


derfinderHelper
---------------

Changes in version 0.99.0:

NEW FEATURES

- Preparing to submit to Bioconductor.

- Added tests and vignette.

derfinderPlot
-------------

Changes in version 0.99.0:

NEW FEATURES

- Preparing to submit to Bioconductor.

- Added tests and vignette.

SIGNIFICANT USER-VISIBLE CHANGES

- plotOverview() and plotCluster() can now plot FWER adjusted p-values if calculated with derfinder::mergeResults()

DESeq2
------

Changes in version 1.6.0:

- DESeq() and results() gets a 'parallel' argument.

- results() gets an 'addMLE' argument.

- results() gets a 'test' argument, for constructing Wald tests after DESeq() was run using the likelihood ratio test.

- results() argument 'format' for GRanges or GRangesList results.

- new plotCounts() function.

- Less outlier calling from Cook's distance for analyses with many samples and many conditions.

- More robust beta prior variance and log fold change shrinkage.



DiffBind
--------

Changes in version 1.12.0:

- Mostly bug fixes!

DMRforPairs
-----------

Changes in version 1.1.1 (2014-05-06):

- The publication describing DMRforPairs was accepted in BMC Bioinformatics. The citation info of the package has been updated
  accordingly. Please use this updated reference when citing DMRforPairs: Rijlaarsdam MA, Zwan YGvd, Dorssers LC, Looijenga
  LH. DMRforPairs: identifying Differentially Methylated Regions between unique samples using array based methylation
  profiles. BMC Bioinformatics. 2014(in press).

DOQTL
-----

Changes in version 0.99.1:

NEW FEATURES

- Added limited support for Heterogeneous Stock mice.

- assoc.plot Added strain distribution patterns above the association mapping plot.

CHANGES

- rankZ Fixed bug relating to NA values.

Changes in version 0.99.0:

NEW FEATURES

- read.vcf reads Sanger SNP VCF files.

- assoc.map imputes Sanger SNPs onto DO genomes and performs association mapping.

- Fixed bug in kinship.probs in which kinship per chromosome was not calculated correctly.

- Improved gene layout algorithm in gene.plot.

CHANGES

- scanone returns p-values and -log10(p-values).

- doqtl.plot plots either LOD or -log10(p-values).

DOSE
----

Changes in version 2.3.6:

- add readable parameter in simplot <2014-10-09, Thu>

Changes in version 2.3.5:

- implement enrichMap <2014-07-28, Wed>

Changes in version 2.3.4:

- return ggplot2 objects in gseaplot <2014-07-21, Mon>

Changes in version 2.3.3:

- geneSim can only accept one gene ID vector and perform as mgeneSim in GOSemSim <2014-06-23, Mon>

- update man files <2014-06-23, Mon>

Changes in version 2.3.2:

- bug fixed in scaleNodeColor <2014-06-08, Sun>

Changes in version 2.3.1:

- upgrade man file according to roxygen2 (ver 4.0.0) . <2014-05-16, Fri>

easyRNASeq
----------

Changes in version 2.1.15:

- Removed an apparently unnecessary require call.

Changes in version 2.1.14:

- Fixed a missing documentation link and a missing object documentation as well as the move of two objects' (DataFrame,
  SimpleList) generic from IRanges to S4Vectors

Changes in version 2.1.13:

- Bioconductor Core Team changes to underlying packages

Changes in version 2.1.12:

- Bioconductor Core Team changes to underlying packages

Changes in version 2.1.11:

- Fixed more IRanges -> S4Vectors import changes.

- Relaxed the BAM header validation of the SO field.  Thanks to John (Zang Jianhua) for finding the issue and providing the
  fix. Same changes as in release 2.0.9.

Changes in version 2.1.10:

- Bioconductor core changes Iranges -> S4Vectors

Changes in version 2.1.9:

- Ported changes from version 2.0.7

Changes in version 2.1.8:

- Ported changes from the release version 2.0.5 and 2.0.6

Changes in version 2.1.7:

- Changed to use roxygen 4.0.0

- Corrected a dependency mismatch.  Beats me why it did not work in v2.1.6...

Changes in version 2.1.6:

- Corrected a dependency mismatch.

Changes in version 2.1.5:

- Same changes as 2.0.3 and adapted dependencies from IRanges to S4Vectors

Changes in version 2.1.4:

- Some haphazard modification from Bioc.

Changes in version 2.1.3:

- Same changes as 2.0.2

Changes in version 2.1.2:

- Completed the previous commit

Changes in version 2.1.1:

- Same changes as 2.0.1

Changes in version 2.1.0:

- None, Bioc new devel branch


EBImage
-------

Changes in version 4.8.0:

NEW FEATURES

- 'otsu' thresholding method (contributed by Philip A. Marais, University of Pretoria, South Africa)

- Support for dimnames in Image objects

- 'bg.col' argument to 'affine' transformations

- 'reenumerate' argument to 'rmObjects'

- 'names' argument to 'readImage'

- 'as.array' method for Image objects

- 'as.nativeRaster' function

SIGNIFICANT USER-VISIBLE CHANGES

- Performance improvements to 'Image', 'selectChannel', 'combine ' and 'reenumerate'

- Use a more efficient 'nativeRaster' representation in 'displayRaster'

- Cleaner output of the 'show-Image' method; print true object class name and dimorder (if set)

- 'readImage' sets Image dimnames to corresponding file names

- 'filter2' and 'affine' return object of the same class as input

- Renamed 'getNumberOfFrames' to 'numberOfFrames'

BUG FIXES

- Handling of dimensions of character arrays

- Drawing of grid lines in 'displayRaster'

- Passing of '...' arguments in 'readImage'

EBSeq
-----

Changes in version 1.5.4:

- An extra numerical approximation step is implemented in EBMultiTest() function to avoid underflow. The underflow is likely
  due to large number of samples. A bug in EBMultiTest() is fixed. The bug will cause error when there is exactly 1
  gene/isoform that needs numerical approximation.

Changes in version 1.5.3:

BUG FIXES

- Fixed a bug that may generate NA FC estimates when there are no replicates.

Changes in version 1.5.2:

NEW FEATURES

- An extra numerical approximation step is implemented in EBTest() function to avoid underflow. The underflow is likely due to
  large number of samples.

EDASeq
------

Changes in version 1.99:

- "exprs" and "exprs<-" are now deprecated: use "counts" and "counts<-" instead (for compatibility with the DESeq class).

- "counts" now accesses the original counts (even after normalization) and "normCounts" accesses the normalized counts.

- Added the slot "normalizedCounts" to the SeqExpressionSet class: now the normalization functions will save the normalized
  counts in this new slot while keeping the original counts in "counts."

- Added option to "MDPlot" to visualize control genes in red.

- withinLaneNormalization and betweenLaneNormalization now always store the normalized counts in the normalizedCounts slot,
  even when offset=TRUE is used.

- Added a new method plotPCA for Principal Components Analysis (PCA) plots.

- Updated the vignette to reflect these changes.

- DESCRIPTION and NAMESPACE cleaned up.

edgeR
-----

Changes in version 3.8.0:

- New goana() methods for DGEExact and DGELRT objects to perform Gene Ontology analysis of differentially expressed genes
  using Entrez Gene IDs.

- New functions diffSpliceDGE(), topSpliceDGE() and plotSpliceDGE() for detecting differential exon usage and displaying
  results.

- New function treatDGE() that tests for DE relative to a specified log2-FC threshold.

- glmQLFTest() is split into three functions: glmQLFit() for fitting quasi-likelihood GLMs, glmQLFTest() for performing
  quasi-likelihood F-tests and plotQLDisp() for plotting quasi-likelihood dispersions.

- processHairpinReads() renamed to processAmplicons() and allows for paired end data.

- glmFit() now stores unshrunk.coefficients from prior.count=0 as well as shrunk coefficients.

- estimateDisp() now has a min.row.sum argument to protect against all zero counts.

- APL calculations in estimateDisp() are hot-started using fitted values from previous dispersions, to avoid discontinuous APL
  landscapes.

- adjustedProfileLik() is modified to accept starting coefficients. glmFit() now passes starting coefficients to
  mglmOneGroup().

- calcNormFactors() is now a S3 generic function.

- The SAGE datasets from Zhang et al (1997) are no longer included with the edgeR package.

EnrichmentBrowser
-----------------

Changes in version 0.99.8:

NEW FEATURES

- Initial release of 'EnrichmentBrowser' package

ensemblVEP
----------

Changes in version 1.6.0:

NEW FEATURES

- Add VEPParam77 class to support Ensembl 77

epivizr
-------

Changes in version 1.3.20:

- Standalone mode introduced, a version of epiviz with reduced capabilities is now included as part of epivizr. The epiviz web
  app is run locally using 'httpuv's http server

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

- Changes slots using lists in 'EpivizDeviceMgr' to environments to avoid crashing RStudio due to inspection of manager
  objects

Changes in version 1.3.5:

- Fails gracefully on daemonization request on Windows

- Deprecates the 'proxy' argument to 'startEpiviz'

Changes in version 1.3.4:

- Upgrading to Epiviz v2 webapp

erccdashboard
-------------

Changes in version 0.99.4:

BUG FIXES AND MINOR IMPROVEMENTS

- Updated saveERCCPlots, description file, and citation file.

Changes in version 0.99.1:

BUG FIXES AND MINOR IMPROVEMENTS

- Fixed bug in erccROC and estLODR code to only use P-values from ERCCs by coding "ERCC-00" prefix instead of "ERCC-" which in
  some cases might find endogenous transcripts named with this prefix.

Changes in version 0.99.0:

NEW FEATURES

- Version set to 0.99.0 for submission to Bioconductor

ISSUES

- CMD check still takes too long (> 5 minutes), probably due to differential expression testing in underlying package. This
  will be difficult to change.

Changes in version 0.9.12:

BUG FIXES AND MINOR IMPROVEMENTS

- Shortened vignette to reduce installation time

- Fixed userMixFile as.factor bug

- Updated namespace to pull in mgcv and nlme package functionality

ISSUES

- CMD check still takes too long (> 5 minutes), probably due to differential expression testing in underlying package. This
  will be difficult to change.

Changes in version 0.9.11:

BUG FIXES AND MINOR IMPROVEMENTS

- Added examples to several functions including "dontrun" to minimize installation time

- Updated Namespace to include gplots

exomePeak
---------

Changes in version 1.4.1 (2014-07-02):

- fixed bug: sometimes bed file is not displayed properly

- added "NEWS" section

- added reference information

FGNet
-----

Changes in version 3.0:

NEW FEATURES

- Graphical User Interface (GUI)

- New FEA tools: topGO and gage (GSA)

- Plot gene expression (up/down) in the functional networks (Available since devel. version 2.1)

- plotGoAncestors, plotKegg

INTERFACE CHANGES 

- Functions renamed or replaced:

- toMatrix/adjMatrix: fea2incidMat

- plotMetagroupsDistance: clusterDistance

- query_gtLinker & getResults_gtLinker: fea_gtLinker & fea_gtLinker_getResults

- query_david & getResults_david: fea_david & format_david

- report_gtLinker & report_david: FGNet_report

- intersectionNetwork: is now included in functionalNetwork

flipflop
--------

Changes in version 1.3.9:

- [pre-processing] redefine extreme boundary of segments (first and last), see calculateBound in readgroup.cpp

Changes in version 1.3.7:

- solve non-portable flags issues, see Makevars

Changes in version 1.3.5:

- solve compilation isue on Mavericks (using clang), see Makevars.

Changes in version 1.3.4:

- properly take into account 'S', 'X' and '=' characters (soft clipping, strict mismatch and match) in the CIGAR of the sam
  file

Changes in version 1.3.2:

- Solve small bugs in the graph construction

- Update the optimization solver for refit. More accurate and faster!

- Remove warnings from spams (to continue ...)

flowBin
-------

Changes in version 1.1.1:

DOCUMENTATION

- Minor fixes to errors or incompleteness in man files.

flowcatchR
----------

Changes in version 0.99.4 (2014-10-03):

NEW FEATURES

- Final fixes. Accepted in the current status at Bioconductor.

Changes in version 0.99.3:

NEW FEATURES

- Changes in the code to be adapted to the latest devel version of EBImage (new functions were added for the Frames class).

Changes in version 0.99.2:

NEW FEATURES

- Updated/added documentation.

Changes in version 0.99.1:

NEW FEATURES

- Adapted after first round of review.

Changes in version 0.99.0:

NEW FEATURES

- Initial submit to Bioconductor.

flowCore
--------

Changes in version 1.31.15:

Deprecation

- 'filterSet' and 'workFlow' are deprecated by 'flowWorkspace::GatingSet'

- update vignette by using 'GatingSet'

Enhancement

- support 'SPILLOVER` in `read.FCS` besides the existing keywords ("SPILL", "spillover")

- support reading(`read.FCS` ) bad FCS files exported by flowJo that do not follow standards strictly

- support 'ncdfFlowList' in 'findTimeChannel' function


flowType
--------

Changes in version 2.3.0:

- To maintain some backwards compatibility, allow users to specify method='thresholds' as either capital or lower case.

flowViz
-------

Changes in version 1.29.30:

bug fixes

- fix the error when xyplot on rare cells

- error when xyplot on rare cell

enhancements

- speed up xyplot by subsampling the data

- support multiple overlays in xyplot

flowWorkspace
-------------

Changes in version 3.11.32:

Enhancements

- stores the axis information in 'axis' slot of 'GatingSet' class to be used by plotGate

- allow different orders of colnames of flow data when merging multiple 'GatingSet's into 'GatingSetList'

- add hidden nodes support

- annotate hidden and boolean nodes in 'gating tree'

- wrap 'getNodes' logic into c++ to speed it up

- add 'flowJoTrans' function to construct the flowJo-type biexponentioal transformation function

- add multiple overlays support to plotGate method

- add raw.scale argument to plotGate

- add 'subset' S3 function to subset the GatingSet/GatingSetList based on 'pData'

- add 'long' format output from 'getPopStats` method

- add 'transform' method for 'GatingSet' to transform the flow data associated with the GatingSet and save the transformation
  functions within 'GatingSet'

- add some internal functions ('merge-GatingSet.R') to handle merging 'GatingSets'

focalCall
---------

Changes in version 0.99.3 (2014-10-06):

IMPROVEMENTS

- Fixed unit tests

- Changed column names in CNVset in example data

Changes in version 0.99.2 (2014-10-06):

IMPROVEMENTS

- Added initial unit testing

- Fixed indenting in .Rd files

- Added this NEWS file

- Changed all hard coded column numbers to name of column

- Added example to BierkensCNA.Rd

- Changed names example data ("BierkensCNA")

- Fixed NAMESPACE file

Changes in version 0.99.1 (2014-10-02):

IMPROVEMENTS

- Changed F to FALSE

- Fixed /value section in man pages

- Added runnable examples for all man pages

- Updated DESCRIPTION file

Changes in version 0.99.0 (2014-04-14):

RELEASE

- initial submitted version

gCMAP
-----

Changes in version 1.9.2:

- BUGFIX: Updated geneIndex method to accept both symbols2indices and id2indices functions from limma, depending on the limma
  version.


GeneAnswers
-----------

Changes in version 2.8.0:

NEW FEATURES

- New function getListGIF is added to interact with ListGIF web server (http://listgif.nubic.northwestern.edu). ListGIF
  provides a quick analysis solution to identify the overrepresented biomedical concepts from a list of query genes. Result is
  rendered as a wordle.

BUG FIXES

- Reactome hyperlink has been fixed.



GeneNetworkBuilder
------------------

Changes in version 1.7.1:

NEW FEATURES

- No changes classified as 'new features' (package under active development)

BUG FIXES

- Fixed no matching function bug under Mac OS X Mavericks


genomeIntervals
---------------

Changes in version 1.21.1:

- Ported version 1.20.1 changes

Changes in version 1.20.1:

- Added an argument to the readGff3 function to enable quietness.

- Corrected some R CMD check warnings

GenomicFiles
------------

Changes in version 1.2.0:

NEW FEATURES

- Add pack / unpack generics and methods

- Add GenomicFiles class

- Add reduceByFile / reduceByRange methods for GenomicFiles class that expect 'file' to be character and 'ranges' a GRanges

- Move yieldReduce() from Rsamtools to GenomicFiles and rename as reduceByYield()

- Allow GRange or GRangesList as @rowData in GenomicFiles class

MODIFICATIONS

- Remove unused .FileList, VCFFileViews and FaFileViews class

- Add checks for 'summarize=FALSE' when REDUCER is used

- Clean up vignette introduction

- Change REDUCER() signature to single argument reguardless of the value of 'iterate'

- Rework reduceByYield() arguments for consistency with other reduceBy* functions

GenomicRanges
-------------

CHANGES IN VERSION 1.18.0
-------------------------

NEW FEATURES

    o Add 'use.mcols' arg to "ranges" method for GRangesList objects. 

    o "assays<-" methods may be invoked with 'withDimnames' arg.

    o Add mapCoords() generic and methods (replacing map()).

    o Add granges,GenomicRanges method. 

    o Add strand<-,GRangesList,character method for global replacement
      (i.e., all strands become 'value'). 

    o Add resize,GRangesList-method.

    o Add DelegatingGenomicRanges class and vignette on how to extend
      GenomicRanges.

    o Document subsetting a named list-like object by a GRanges subscript.

SIGNIFICANT USER-LEVEL CHANGES

    o Modify "show" methods for GRanges and GRangesList objects so
      they print a 1-line summary of the seqinfo component.
      
    o Remove as.data.frame,GRangesList-method; use as.data.frame,List.

    o "trim" method for GenomicRanges only trims out-of-bound
      ranges on non-circular sequences whose length is not NA.
      This behavior is consistent with the GenomicRanges validity method.

    o Changes to flank(), resize() and start/end/width setters:
      - no longer trim the result ranges when called on a GRanges
      - warning is issued by GenomicRanges validity method when
        out-of-bound ranges are on non-circular sequences whose
        length is not NA
      Note this behavior is now consistent with that of shift().

    o Speed up validation of GenomicRanges objects by 1.2x. 

    o Speed up trim() on GenomicRanges objects by 1.2x. 

    o Improve warning when GenomicRanges object contains out-of-bound ranges. 

    o Work on vignette HOWTOs:
      - split 'How to read BAM files into R' into 3 HOWTOs
      - split 'How to prepare a table of read counts for RNA-Seq
        differential gene expression' into 3 HOWTOs
      - split 'How to extract DNA sequences of gene regions' into 2 HOWTOs 
      - make individual HOWTOs subsections of single HOWTO section

    o Follow renaming of TranscriptDb class to TxDb. 

    o Replace references to plantsmart21 with plantsmart22.

DEPRECATED AND DEFUNCT

    o Defunct map() (skip deprecation). Replace with mapCoords().

BUG FIXES

    o [cr]bind,SummarizedExperiment methods respect derived classes.

    o assays(se, withDimnames=TRUE) <- value no longer tries to access
      a slot 'withDimnames'.

    o cbind and rbind,SummarizedExperiment-methods respect derived
      classes

    o "ranges" method for GRangesList objects should not propagate
      inner metadata columns by default.

    o GRanges() constructor now preserves the seqlevels in the order
      supplied by the user.

    o Ensure tileGenome() breakpoints do not extend past end of genome. 

    o Fix "show" method for GenomicRanges objects when 'showHeadLines'
      global option is set to Inf.

    o [rc]bind,SummarizeExperiment-methods now compare all elements.
    
    o Remove "==" and "<=" methods for GenomicRanges objects (not needed). 


genoset
-------

Changes in version 1.19.27:

NEW FEATURES

- rangeMeans for numeric gets a na.rm argument and uses 50% less time for na.rm = FALSE ands 25% less time for na.rm=TRUE.

Changes in version 1.19.7:

NEW FEATURES

- The version introduces the rangeMeans family of functions, which are closely related to IRanges::viewMeans family of
  functions. These new functions take an RleDataFrame and an IRanges (or matrix) of row ranges and compute summary stats.
  These are useful for summaries like the average DNA copy number or coverage in a gene, by sample. These functions differ
  from view* in that they are more time efficient and they offer a 'simplify' argument to return a matrix, rather than a list.

Changes in version 1.19.0:

NEW FEATURES

- This version introduces the RleDataFrame class. This class extends SimpleRleList and DataFrame. It stores a collection of
  like-length Rle objects.  With the AtomicList features it behaves much like a matrix. For example, log2(x) - 1 works. We
  also have row and column sums and means. rangeSampleMeans functions like viewMeans. This object can be used to store runs of
  data along the genome for multiple samples, like coverage or DNA copy number. It can be used as an assayDataElement in a
  genoset (or any eSet).

- segTable and segPairTable are dramatically faster. Stacking data.frames with segTable(RleDataFrame,stack=TRUE) argument is
  now instantaneous.

DEPRECATED AND DEFUNCT

- Previously we used a DataFrame for storing collections of Rle objects.  This practice is now deprecated. Similarly, colMeans
  on a DataFrame is now deprecated.

GGBase
------

Changes in version 3.27.5:

- clipPCs and regressOut now work for SummarizedExperiment instances


gmapR
-----

CHANGES IN VERSION 1.8.0
-----------------------

NEW FEATURES

    o GmapGenomes can be built from any file supported by rtracklayer
      (so 2bit now works, as well as fasta).

    o Tally BAM files by codon given a set of transcript
      structures. This happens at the read level, i.e., a codon is
      observed within an individual read.

    o Tally BAM files strand by XS tag (inferred strand of
      transcription, instead of strand of alignment).



GOexpress
---------

Changes in version 0.99.17:

BUG FIXES

- NAs were introduced in the average score and rank of GO terms following analysis of microarray data (ANOVA and randomForest)
  for GO terms without associated features in the ExpressionSet. The problem was not found in any case for analysis of RNA-Seq
  data. This was causing issues during the subset_scores function and subsequent plots, such as main titles made of multiple
  lines and NAs. GO terms without associated probesets are now given average score of 0 and average rank equal to the maximum
  in the dataset plus 1.

Changes in version 0.99.16:

GENERAL UPDATES

- Included co-authors who participated in the generation and analysis of data used to test the package.

Changes in version 0.99.15:

GENERAL UPDATES

- Updated man pages GO_analyse to describe the recently added subset slot in the result variable.

Changes in version 0.99.14:

BUG FIXES

- The anova method of the GO_analyse method was broken since the introduction of ExpressionSet, release 0.99.4.

UPDATED FEATURES

- The subset argument was also added to the GO_analyse method. This allows the identification of genes clustering a subset of
  groups, at a subset of time-points, etc. while plotting the expression profile of the entire dataset, if desired.

GENERAL UPDATES

- Updated man pages GO_analyse to allow only one example to be run. This saves time during CMD check, while making the man
  page more easily readable.

Changes in version 0.99.13:

UPDATED FEATURES

- All expression plots were updated to a default ylim range corresponding to the minimum and maximum expression values found
  in the entire ExpressionSet. This is meant to avoid mis-interpretation of the amplitude of variation between sample groups,
  as suggested in the paper: Rougier, N.P., Droettboom, M., and Bourne, P.E. (2014). Ten simple rules for better figures. PLoS
  computational biology 10, e1003833.

Changes in version 0.99.12:

UPDATED FEATURES

- All expression plots were given new arguments for more control: xlab allows users to change the default title for the
  X-axis, ylim allows users to override the lower and upper boundaries of the Y axis.

GENERAL UPDATES

- Updated man pages AlvMac, AlvMac_results, microarray2dataset, and prefix2dataset. Replaced the LaTeX describe statement by
  an itemise statement to make the man page more readable in a terminal window.

- Moved functions between the post_analysis script and the toolkit script. From now on, only functions accessible to the users
  should be present in the post_analysis script, while toolkit should contain method called internally.

- The User's Guide was updated to redirect the users to the new support site for Bioconductor rather than the bioc-devel
  mailing list.

Changes in version 0.99.11:

NEW FEATURES

- The new subEset() method subsets an ExpressionSet given a list where item names are column names of the phenoData slots and
  item values are vectors of values corresponding to sample to retain (e.g. list(Time=c("2H", "6H")) will retain samples with
  value "2H" or "6H" in the "Time" column of the phenoData slot).

UPDATED FEATURES

- All relevant visualisation methods have been added an argument to subset samples to plot to those with a given set of values
  for a given column in the phenoData (uses the new function subEset described above).

- The example analysis results was renamed from "raw_results" to "AlvMac_results", so that it can now loaded running
  data(AlvMac_results). The new name is meant to be more specific to the package.

Changes in version 0.99.10:

GENERAL UPDATES

- Replaced \format by \value sections to the following man pages: AlvMac.Rd, microarray2dataset.Rd, prefix2dataset.Rd,
  raw_results.Rd to avoid having an empty value section, as recommended by the Bioconductor package tracker. Hopefully, it
  does not require a non-empty \format section as well...

- Restricted lines to less than 80 characters and indentation by multiple of four space characters.

Changes in version 0.99.9:

UPDATED FEATURES

- all post analysis functions were given more sanity checks to verify that the "result" argument contains the required slots
  of a GO_analyse() output and the the arguments pointing at a phenotypic data column are valid column names.

BUG FIXES

- expression_plot_symbol and expression_profiles symbol used the default value of col.palette and colourF instead of
  forwarding the user-defined one to the expression_plot and expression_profiles functions.

GENERAL UPDATES

- Added cross-reference in UsersGuide to point at examples of usages of factor and numeric values for the expression plots

Changes in version 0.99.8:

GENERAL UPDATES

- Resaved R data files to reduced package disk size. No more WARNING in R CMD check.

- Removed reference to GitHub in the README file. The weblink is given in the DESCRIPTION file anyway if users are interested.

Changes in version 0.99.7:

UPDATED FEATURES

- fixed a typo in the code of heatmap_GO which made it crash for any other dataset than the example dataset.

Changes in version 0.99.6:

UPDATED FEATURES

- expression_profiles_symbol() method was missing the "index" argument to select the feature identifier to plot alone.

Changes in version 0.99.5:

NEW FEATURES

- expression_profiles() method plots the expression profiles of individual samples series, as opposed to grouped samples
  series handled by expression_plot functions.

- expression_profiles_symbol() method plots the expression profiles of individual samples series using a gene name instead of
  an Ensembl gene identifier.

UPDATED FEATURES

- overlap_GO can print to screen, if filename argument is set to NULL (Default).

- heatmap_GO, cluter_GO and plot_design can resize title font and wrap the text on multiple lines.

- expression_plot and expression_plot_symbol can orient X axis labels at a given angle.

- replaced return(NULL) statement by stop() when no close match is found to a gene name in the family of expression_plot
  functions.

GENERAL UPDATES

- User's Guide updated.

- List of contributors updated in User's Guide and DESCRIPTION.

Changes in version 0.99.4:

UPDATED FEATURES

- Use of ExpressionSet instead of numeric named matrix and AnnotatedDataFrame. Better consistency with other Bioconductor
  packages.

GENERAL UPDATES

- Implemented corrections requested following the Bioconductor review. Includes typos, consistent terminology through the
  package code and metafiles, additional information in help files, no reference to GitHub as an alternate installation
  option, use of arrow signs instead of equal signs for value assignment.

- Restricted lines to 80 characters, and used 4-space tabulations.

- Corrected out-of-date documentation.

Changes in version 0.99.3:

UPDATED FEATURES

- Control the size of the legend text in the two expression plot figures. Updated help files accordingly.

- Updated vignette with new section "Statistics".

- Complete cleaning of code files for lines shorter than 80 columns.

- Cleanup of help files for lines shorter than 80 columns.

- Enabled filtering of raw results on the average score of a GO term.

Changes in version 0.99.2:

UPDATED FEATURES

- Metadata lines in the preamble of the Sweave file

Changes in version 0.99.1:

UPDATED FEATURES

- Sweave vignette implemented.

- Replaced all message() statements by cat() to make Sweave output the full message in the vignette.

- Updated a missed F into FALSE

- Updated an invalid biocViews (typo)

- Date field added for a proper citation() method.

Changes in version 0.99.0:

UPDATED FEATURES

- Replaced all cat() statements by message() to match the Biocondcutor guidelines.


GOSemSim
--------

Changes in version 1.23.2:

- fast IC-method implemented, contributed by Alexey Stukalov <2014-09-27, Sat>

Changes in version 1.23.1:

- add support of species coelicolor and gondii <2014-09-03, Wed>

GOsummaries
-----------

Changes in version 1.99.3:

CHANGES

- Bumped the version number to 1.99.3

Changes in version 1.1:

- Added support for custom data in word clouds (see gosummaries functions)

- Now it is possible to display genes instead of GO categories (see show_genes parameter in gosummaries.MArrayLM,
  gosummaries.prcomp and gosummaries.matrix)

- Added function gosummaries.matrix that takes in a matrix that is a MDS representation of data and expression matrix and then
  finds most correlated features for each MDS component.

Changes in version 1.0:

NEW FEATURES

- First version of the package

Changes in version 0.99.3:

NEW FEATURES

- Updated Vignette to biocStyle

- Form now GOsummaries is going to live in Bioconductor and this version supersedes the CRAN version 1.1

Changes in version 0.99.2:

NEW FEATURES

- Reformatted code for Bioconductor

- Improved Vignette

graphite
--------

Changes in version 1.11.4 (2014-10-01):

- Updated Biocarta, KEGG, NCI, HumanCyc, Panther and Reactome data.

GSAR
----

Changes in version 1.0.0:

- This is the first version of the R package GSAR.

- The package provides two-sample nonparametric multivariate statistical methods to test specific alternative hypotheses
  against a null hypothesis.

- GSAR depends on package igraph to handle graphs in objects of class igraph and uses some functions too.

- New capabilities and future changes will be reported in subsequent versions.

Gviz
----

Changes in version 1.9.0:

NEW FEATURES

- The new CustomTrack class to allow for user-defined plotting functions.

SIGNIFICANT USER-VISIBLE CHANGES

- The collapseTranscripts parameter now offers more control over the type of collapsing.

gwascat
-------

Changes in version 1.9.8:

USER VISIBLE CHANGES

- makeCurrentGwascat has new arguments, useHg38seqinfo and altseqinfo. These address the fact that the textual version of the
  catalog served by NHGRI has hg38 addresses.  Two snapshots of the data are available, created Sept 8 2014, one direct in
  hg38 (gwrngs38 in data) and the other by liftOver to hg19 addresses (gwrngs19)

GWASTools
---------

Changes in version 1.11.33:

- Allow getting variables from sub-nodes in a GDS file (e.g., getVariable(GdsReader, "snp.annot/qual")).

- Add getNodeDescription method to GdsReader.

- Added examples of converting from PLINK and VCF in Formats vignette.

Changes in version 1.11.32:

- imputedDosageFile replaces ncdfImputedDosage and gdsImputedDosage

Changes in version 1.11.31:

- setMissingGenotypes replaces ncdfSetMissingGenotypes and gdsSetMissingGenotypes

Changes in version 1.11.22:

- convertNcdfGds and convertGdsNcdf will convert files with any variable names (not just genotype)

Changes in version 1.11.21:

- Fixed bug in vcfWrite to output missing data code for ID column

- Data cleaning vignette uses createDataFile instead of ncdfCreate and ncdfAddData

- Data cleaning vignette uses snpgdsOpen and snpgdsClose

Changes in version 1.11.20:

- Fixed bug in assocTestCPH when there is no Y chromosome in the data

Changes in version 1.11.19:

- convertNcdfGds will not write entire snp and sample annotations to file

- createDataFile replaces ncdfCreate and ncdfAddData

Changes in version 1.11.18:

- patch from Karl Forner to allow use of open gds objects in constructors for GdsReader and GdsGenotypeReader

Changes in version 1.11.17:

- removed duplicated .probToDosage function from ncdfImputedDosage.R source file

Changes in version 1.11.16:

- expanded matching options in duplicateDiscordanceAcrossDatasets

Changes in version 1.11.15:

- allowed truncate to be a numeric value or TRUE in qqPlot

Changes in version 1.11.14:

- added pasteSorted function

Changes in version 1.11.13:

- in case of missing allele code, return character genotype as NA

Changes in version 1.11.12:

- bug fix in assocTestRegression when a block contains only 1 SNP

Changes in version 1.11.11:

- added vcfCheck function to compare VCF file to GenotypeData object

Changes in version 1.11.10:

- bug fix in gwasExactHW when a block contains only 1 genotype

Changes in version 1.11.9:

- changed colors of BAF plots so points can be more easily distinguished

Changes in version 1.11.8:

- added ref.allele option to vcfWrite to select either A or B as the reference allele for each SNP

Changes in version 1.11.7:

- added vcfWrite function to write VCF file from GenotypeData object

Changes in version 1.11.6:

- bug fix in qqPlot, manhattanPlot when requesting thinning when bins only have 1 point

Changes in version 1.11.5:

- added pointsPerBin argument to manhattanPlot

Changes in version 1.11.4:

- added optional thinThreshold argument to manhattanPlot and qqPlot functions

Changes in version 1.11.3:

- updated gdsSubset for new gdsfmt read.gdsn syntax (also changed in release version)

Changes in version 1.11.1:

- Added ylim argument to qqPlot.

HDTD
----

Changes in version 0.99.4 (2014-10-01):

- Updated CITATION file.

Changes in version 0.99.3 (2014-09-25):

- Introduced warning messages regarding the sample size.

- Updated the vignette file.

Changes in version 0.99.2 (2014-08-12):

- Updated the output of the core functions.

Changes in version 0.99.1 (2014-05-14):

- Accepted to Bioconductor.

Changes in version 0.99.0 (2014-04-17):

- Submitted to Bioconductor.

hiAnnotator
-----------

Changes in version 0.99.9:

- metadatacol bug fix

Changes in version 0.99.8:

- plotdisFeature() obtains a geom parameter.

Changes in version 0.99.7:

- Documentation fixes

Changes in version 0.99.6:

- Using GRanges::sort in makeGRanges().

Changes in version 0.99.5:

- Documentation changes.

Changes in version 0.99.4:

NEW FEATURES

- plotdisFeature(), function to plot distance distribution to feature boundary.

- Includes sites.ctrl dataset to compliment sites dataset.

Changes in version 0.99.3:

NEW FEATURES

- makeGRanges() inherits total functionality of makeRangedData() and an option to turn off factor to character conversion.

DEPRECATED AND DEFUNCT

- Functions resizeRangedData() & makeRangedData() has been removed in favor of using only GenomicRanges-derived objects.

BUG FIXES

- Improved seqinfo slot population method in makeGRanges()

hiReadsProcessor
----------------

0.99.0: Added examples for almost all functions

0.99.0: Removed Subread support

HiTC
----

Changes in version 1.9.5:

NEW FEATURES

- The normLGF function can now be applied both on intra and inter-chromosomal maps

SIGNIFICANT USER-VISIBLE CHANGES

- Access the contact map is now ALWAYS performed using chromosomes' name such as ygi/xgi, i.e. rownames/colnames of contact
  maps

- Update of the setGenomicFeatures method and speed improvment

BUG FIXES

- Update of isComplete method

- Update of import functions

Changes in version 1.9.4:

NEW FEATURES

- ICE normalization can now be applied on both HTCexp and HTClist objects

- New methods for HTClist-class : isComplete, isPairwise, forcePairwise, forceSymmetric

- New methods for HTCexp-class : forceSymmetric

- New methods (not exported) for HTClist-class : getCombinedIntervals, getCombinedContacts

SIGNIFICANT USER-VISIBLE CHANGES

- Update of the binningC method. Improvement of speed and memory usage. The binningC function can now be applied on Hi-C data
  (or any already binned data). The goal is to move from a very high resolution map (for instance 40kb) to a lower resolution
  ( for instance 1Mb) by aggregationg and summing the bins.

- Update of the importC/exportC functions which are now based on a new format (list with interactor1/interactor2/count + BED
  files). This format is recommanded to store sparse data because only the non null values are exported/imported.

- Update of the import.my5C/export.m5C functions. Only the matrix format is now supported. For the list format, see the
  importC/exportC functions.

BUG FIXES

- Bug fixed in divide method

Changes in version 1.9.3:

NEW FEATURES

- getRestrictionFragmentsPerChromosome

SIGNIFICANT USER-VISIBLE CHANGES

- Update output of summary method to add the seqlevels of both interactors

- Speed improvement; seqlevels

- Update import.my5C function

BUG FIXES

- Fix bug in HTCexp constructor. Only intrachromosomal data can be force to be symmetrical

- Fix bug in summary function for HTClist object

- Fix bug in mapC function, in case of empty matrices (only zero values)

HTSeqGenie
----------

Changes in version 3.15.3:

- parallelize Indel realignment by chromosome

Changes in version 3.15.2:

- use a binary tree-reduce like algorithm to merge coverages.  For a big WGS this reduced runtime from 12h to 2 h.

Changes in version 3.15.1:

- start next dev cycle

- remove count_transcripts and count_ncRNA_nongenic from QA

Changes in version 3.14.1:

- 3.14.1 release

- fixed the fragmentLength NA bug spotted by Jinfeng

illuminaio
----------

Changes in version 0.7.2 (2014-10-02):

- Now on GitHub.

Changes in version 0.7.1 (2014-09-21):

- The vignette now reads the GenomeStudio example file from the IlluminaDataTestFiles package instead of the internet.

- Updated CITATION.

- Added citation to vignette.

Changes in version 0.7.0 (2014-04-11):

- The version number was bumped for the Bioconductor devel version, which is now BioC v2.15 for R (>= 3.1.0).

IMPCdata
--------

Changes in version 1.0.1:

FEATURES

- Package allows systematically explore the IMPC dataset's multiple dimensions until the correct combination of filters has
  been selected.

- Package helps to obtain datasets from IMPC database.

- The data retrieved from IMPCdata package can be directly used by PhenStat -- an R package that encapsulates the IMPC
  statistical pipeline, available at <URL: http://www.bioconductor.org/packages/release/bioc/html/PhenStat.html>.

inSilicoMerging
---------------

Changes in version 1.10.0:

- Method DWD removed because the DWD package has been archived by CRAN

IRanges
-------

CHANGES IN VERSION 1.20.0
-------------------------

NEW FEATURES

    o Add IntervalForest class from Hector Corrada Bravo.

    o Add a FilterMatrix class, for holding the results of multiple filters.

    o Add selfmatch() as a faster equivalent of 'match(x, x)'.

    o Add "c" method for Views objects (only combine objects with same
      subject).

    o Add coercion from SimpleRangesList to SimpleIRangesList.

    o Add an `%outside%` that is the opposite of `%over%`.

    o Add validation of length() and names() of Vector objects.

    o Add "duplicated" and "table" methods for Vector objects.

    o Add some split methods that dispatch to splitAsList() even when only
      'f' is a Vector.

    o Add set methods (setdiff, intersect, union) for Rle.

    o Add anyNA methods for Rle and Vector.

    o Add support for subset(), with(), etc on Vector objects,
      where the expressions are evaluated in the scope of the
      mcols and fixed columns. For symbols that should resolve
      in the calling frame, it is supported and encouraged to escape
      them with bquote-style ".(x)".

    o Add "tile" generic and methods for partitioning a ranges object
      into tiles; useful for iterating over subregions.

SIGNIFICANT USER-VISIBLE CHANGES

    o All functionalities related to XVector objects have been moved to the
      new XVector package.

    o Refine how isDisjoint() handles empty ranges.

    o Remove 'keepLength' argument from "window<-" methods.

    o unlist( , use.names=FALSE) on a CompressedSplitDataFrameList object
      now preserves the rownames of the list elements, which is more
      consistent with what unlist() does on other CompressedList objects.

    o Splitting a list by a Vector just yields a list, not a List.

    o The rbind,DataFrame method now handles the case where Rle and vector
      columns need to be combined (assuming an equivalence between Rle and
      vector). Also the way the result DataFrame is constructed was changed
      (avoids undesirable coercions and should be faster).

    o as.data.frame.DataFrame now passes 'stringsAsFactors=FALSE' and
      'check.names=!optional' to the underlying data.frame() call.
      as(x,"DataFrame") sets 'optional=TRUE' when delegating. Most places
      where we called as.data.frame(), we now call 'as(x,"data.frame")'.

    o The [<-,DataFrame method now coerces column sub-replacement value to
      class of column when the column already exists.

    o DataFrame() now automatically derives rownames (from the first argument
      that has some). This is a fairly significant change in behavior, but it
      probably does better match user behavior.

    o Make sure that SimpleList objects are coerced to a DataFrame with a
      single column. The automatic coecion methods created by the methods
      package were trying to create a DataFrame with one column per element,
      because DataFrame extends SimpleList.

    o Change default to 'compress=TRUE' for RleList() constructor.

    o tapply() now handles the case where only INDEX is a Vector (e.g.
      an Rle object).

    o Speedup coverage() in the "tiling case" (i.e. when 'x' is a tiling
      of the [1, width] interval). This makes it much faster to turn into an
      Rle a coverage loaded from a BigWig, WIG or BED as a GRanges object.

    o Allow logical Rle return values from filter rules.

    o FilterRules no longer requires its elements to be named.

    o The select,Vector method now returns a DataFrame even when a single
      column is selected.

    o Move is.unsorted() generic to BiocGenerics.

DEPRECATED AND DEFUNCT

    o Deprecate seqselect() and subsetByRanges().

    o Deprecate 'match.if.overlap' arg of "match" method for Ranges objects.

    o "match" and "%in%" methods that operate on Views, ViewsList, RangesList,
      or RangedData objects (20 methods in total) are now defunct.

    o Remove previously defunct tofactor().

BUG FIXES

    o The subsetting code for Vector derivatives was substancially refactored.
      As a consequence, it's now cleaner, simpler, and [ and [[ behave more
      consistently across Vector derivatives. Some obscure long-standing bugs
      have been eliminated and the code can be slightly faster in some
      circumstances.

    o Fix bug in findOverlaps(); zero-width ranges in the query no longer
      produce hits ever (regardless of 'maxgap' and 'minoverlap' values).

    o Correctly free memory allocated for linked list of results compiled for
      findOverlap(select="all").

    o Various fixes for AsIs and DataFrames.

    o Allow zero-row replacement values in [<-,DataFrame.

    o Fix long standing segfault in "[" method for Rle objects (when doing
      Rle()[0]).

    o "show" methods now display its most specific class when a column or
      slot is an S3 object for which class() returns more than one class.

    o "show" methods now display properly cells that are arrays.

    o Fix the [<-,DataFrame method for when a value DataFrame has matrix
      columns.

    o Fix ifelse() for when one or more of the arguments are Rle objects.

    o Fix coercion from SimpleList to CompressedList via AtomicList
      constructors.

    o Make "show" methods robust to "showHeadLines" and "showTailLines" global
      options set to NA, Inf or non-integer values.

    o Fix error condition in eval,FilterRules method.

    o Corrected an error formatting in eval,FilterRules,ANY method.






kebabs
------

Changes in version 1.0.0:

- first official release as part of Bioconductor 3.0

KEGGprofile
-----------

Changes in version 1.7.6:

NEW FEATURES

- New function: Downloading the latest pathway/gene link from KEGG website.

Changes in version 1.7.5:

NEW FEATURES

- Changes in examples.

Changes in version 1.7.4:

NEW FEATURES

- Improvement for documents and examples.

Changes in version 1.7.3:

NEW FEATURES

- Improvement for the web interface.

Changes in version 1.7.2:

NEW FEATURES

- New function plot_pathway_cor: Caculating correlations between groups in each pathway and comparing with the correlations
  for other genes.

Changes in version 1.7.1:

NEW FEATURES

- Compound data was also supported; The parameters were modified for web interface.

limma
-----

Changes in version 3.22.0:

- New functions goana() and topGO() provide gene ontology analyses of differentially genes from a linear model fit. The tests
  include the ability to adjust for gene length or abundance biases in differential expression detection, similar to the goseq
  package.

- Improvements to diffSplice. diffSplice() now calculates Simes adjusted p-values for gene level inferences, in addition to
  the exon level t-tests and gene level F-tests. topSplice() now has three ranking methods ("simes", "F" or "t"), with "simes"
  now becoming the default. diffSplice() also has a new argument 'robust' giving access to robust empirical Bayes variance
  moderation.

- New function plotExons() to plot log-fold-changes by exon for a given gene.

- New function voomWithQualityWeights() allows users to estimate sample quality weights or allow for heteroscedasticity
  between treatment groups when doing an RNA-seq analysis.

- Improvement to arrayQualightyWeights(). It now has a new argument 'var.design' which allows users to model variability by
  treatment group or other covariates.

- Improved plotting for voomaByGroup().

- barcodeplot() can now plot different weights for different genes in the set.

- Improvements to roast() and mroast(). The directional (up and down) tests done by roast() now use both the original
  rotations and their opposite signs, effectively doubling the number of effective rotations for no additional computational
  cost. The two-sided tests are now done explicitly by rotation instead of doubling the smallest one-sided p-value. The
  two-sided p-value is now called "UpOrDown" in the roast() output. Both functions now use a fast approximation to convert
  t-statistics into z-scores, making the functions much faster when the number of rotations or the number of genes is large.
  The contrast argument can now optionally be a character string giving a column name of the design matrix.

- zscoreT() can optionally use a fast approximation instead of the slower exact calculation.

- symbols2indices() renamed to ids2indices().

- Improvements to removeBatchEffect(). It can now take into account weights and other arguments that will affect the linear
  model fit. It can now accept any arguments that would be acceptable for lmFit(). The behavior of removeBatchEffect() with
  design supplied has also changed so that it is now consistent with that of lmFit() when modelling batches as additive
  effects. Previously batch adjustments were made only within the treatment levels defined by the design matrix.

- New function plotWithHighlights(), which is now used as the low-level function for plotMA() and plot() methods for limma
  data objects.

- The definition of the M and A axes for an MA-plot of single channel data is changed slightly.  Previously the A-axis was the
  average of all arrays in the dataset - this has been definition since MA-plots were introduced for single channel data in
  April 2003.  Now an artificial array is formed by averaging all arrays other than the one to be plotted.  Then a
  mean-difference plot is formed from the specified array and the artificial array.  This change ensures the specified and
  artificial arrays are computed from independent data, and ensures the MA-plot will reduce to a correct mean-difference plot
  when there are just two arrays in the dataset.

- plotMDS() can now optionally plot samples using symbols instead of text labels. It no longer has a 'col' argument, which
  instead is handled by ....

- vennDiagram() now supports circles of different colors for any number of circles.  Previously this was supported only up to
  three sets.

- getEAWP() will now find a weights matrix in an ExpressionSet object if it exists.

- update to helpMethods().

- Substantial updates to the two RNA-seq case studies in the User's Guide. In both cases, the short read data has been
  realigned and resummarized.

- Improvements to many Rd files. Many keyword entries have been revised. Many usage and example lines been reformated to avoid
  over long lines.

- biocViews keywords updated.

- Subsetting columns of a MArrayLM object no longer subsets the design matrix.

- Bug fix for read.maimages: default value for 'quote' was not being set correctly for source="agilent.mean" or
  source="agilent.median".

- Bug fix to topTableF() and topTable(). The ordering of Amean values was sometimes incorrect when sorting by F-statistic and
  a lfc or p.value filter had been set.

- Bug fix to read.ilmn() when sep=",".

MBAmethyl
---------

Changes in version 0.99.1:

- A package ready for resubmission.


MEIGOR
------


0.9.4: New functions: vns_default function, MEIGO, vns_optset, get_VNS_settings, rosen10

0.9.4: Lots of clean-up in other functions with some new prototypes: e.g., CeSSR, rvnds_hamming, rvnds_local.R,
             ssm_localsolver.R


meshr
-----

Changes in version 1.0.4:

- Bug Fix (<NA> rows in using with MeSH.db and org.MeSH.XXX.db)

Changes in version 1.0.3:

- documentation improvements

Changes in version 1.0.2:

- documentation improvements

Changes in version 1.0.1:

- documentation improvements

Metab
-----

CHANGES IN VERSION 0.99
-----------------------

 o Function clean.fix renamed to MetReport

 o The function raw.peaks has been removed.

 o Improved speed and error messages.

 o MetReport is now able to extract directly the 
   area and/or the base peak of an AMDIS report.

 o The new function MetReportNames allows users to 
   extract sample's abundances from an AMDIS report 
   based only in sample's names.

 o The new function buildLib allows users to convert an 
   AMDIS library to a csv file in the format required by 
   Metab.

 o Function Htest allows users to adjust p-values for multiple comparisons.




metagenomeSeq
-------------

Changes in version 1.7 (2014-05-07):

- Added function plotBubble

- Added parallel (multi-core) options to fitPA, fitDO

- Fixed bug for fitMeta when useCSSoffset=FALSE and model matrix ncol==2

- (1.7.10) Updated default quantile estimate (.5) for low estimates

- (1.7.10) Added short description on how to do multiple group comparisons

- (1.7.15) Output of fitZig (eb) is now a result of limma::eBayes instead of limma::ebayes

- (1.7.16) plotMRheatmap allows for sorting by any stat (not just sd)

- (1.7.18) fitTimeSeries Including times series method for differentially abundant time intervals

- (1.7.20) Fixed minor bug for OTU level time series analyses and added plotClassTimeSeries

- (1.7.26) Added warning / fix if any samples are empty in cumNormStat

- (1.7.27) Added a few unit tests

- (1.7.29) Added interactiveDisplay to namespace (display function allows interactive exploration / plots through browser)

metaseqR
--------

Changes in version 1.3.5 (2014-09-30):

NEW FEATURES

- Re-analysis based on saved gene models is now faster.

- Added the ability to save time-consuming part of analyses also when count.type="gene".

BUG FIXES

- Fixed problem with RMySQL dependency when using annotation from UCSC. Solved with the usage of RSQLite, however, the process
  is a bit longer.

Changes in version 1.3.4 (2014-09-03):

NEW FEATURES

- Added support for hg38

BUG FIXES

- Fixed more problems occured with the change of biomaRt attributes for newer organisms.

- Fixed bug causing report crash when no genes are passing FDR threshold in paired comparisons.

Changes in version 1.3.3 (2014-08-21):

NEW FEATURES

- Added the ability to retrieve annotation for genes/exons from UCSC or RefSeq through connection to UCSC Genome Browser
  public SQL database. The GC-content for each region is retrieved through BSgenome packages.

- Added details regarding the number of genes returned by each algorithm when conducting combined analysis.

BUG FIXES

- Fixed minor bug with warning level logging.

- Fixed problem occured with the change of biomaRt attributes for newer organisms.

Changes in version 1.3.2 (2014-05-05):

NEW FEATURES

- Function to merge exons belonging to different isoforms to a set of "virtual" exons to help construct a single gene model
  with unique exons

- Simplified the usage of read2count()

BUG FIXES

- Major bug in "exon" mode that inflated the number of reads for certain genes with many isoforms

- Minor bug in argument checking, not allowing to not save the gene model in "exon" mode

MethylAid
---------

Changes in version 0.99.9:

BUG FIXES

- fixed outlier detection using all quality control plots using an initialization method

MethylMix
---------

Changes in version 0.99.1 (2014-08-08):

- This marks the first release to BioConductor of MethylMix. We included two data sets to test MethylMix: a glioblastoma and
  breast cancer data set. We have included one vignettes which uses a breast cancer example to show the functionality of
  MethylMix.

methylPipe
----------

Changes in version 1.0.0:

- initial version with the following functions implemented: + mCsmoothing + findPMDs + plotAnnoBar + mapBSdata2GRanges +
  findDMR + consolidateDMRs + getCpos + profileDNAmetBin + plotMeth + BSprepare + meth.call + GElist-class +
  GEcollection-class + BSdataSet-class + BSdata-class

MGFM
----

Changes in version 0.99.2:

- replaced lookUp by select

Changes in version 0.99.1:

- Some Improvements

Changes in version 0.99.0:

- Finalizing for BioConductor release

minfi
-----

Changes in version 1.11:

- Updated CITATION.

- Added dropLociWithSnps for easy exclusion of certain methylation loci.

- Add getAnnotationObject for easy printing of contents of the annotation object.

- Changes in 1.10 imported into 1.11.

- Fixed an issue with bumphunter calling the bumphunter package in a wrong way.

- Added getOOB and getSnpBeta convenience functions for accessing the OOB probes and the SNP probes.

- read.450k.sheet now forces a column named Slide to be character.

- The NOOB background correction method is now available throguh preprocessNoob.

- One can now supply the permutations to be used in permutation analysis. This is useful for cases in which the total number
  of possilbe permutations is small and one wants to use them all or in cases in which one wants to assure balance, for
  example, between cases and controls.

- The bumphunter method now has the option to create null distributions using a bootstrap approach.

- Fixed a man page issue.

- Added GitHub URL to DESCRIPTION.

- Functional normalization now supports background correction by NOOB (see preprocessNoob); this is recommended (and the new
  default).

missMethyl
----------

Changes in version 0.99.9:

- This version of the package has all calculations in matrix formulation and handles nuisance parameters appropriately.

Changes in version 0.99.7:

- Added new function, contrasts.varFit.

Changes in version 0.99.6:

- Modified getLeveneResiduals

Changes in version 0.99.5:

- Modified the varFit function to accept DGEList objects for differential variability testing for RNA-Seq data.

Changes in version 0.99.4:

- First version of package contains functions to perform SWAN normalisation and differential variability analysis for DNA
  methylation data from Illumina's Infinium HumanMethylation450 beadchip.

monocle
-------

0.99.6: plot_genes_jitter, plot_genes_in_pseudotime, and plot_genes_positive_cells now accept a new parameter,
             labe_by_short_name, to help control faceting

0.99.6: orderCells() now accepts the name of the root cell, so you can fix the beginning of the pseudotime trajectory

0.99.6: CellDataSet objects accept an argument for the VGAM family to be used for the distribution of expression values
             (e.g. negbinomial())

0.99.6: Dimensionality reduction, plotting, model fitting, and differential analyis routines now recognize CellDataSet
             expression families and alter their behavior when expression is count based.  Allows analysis of absolute,
             as opposed to relative, single cell expression data.

0.99.5: Final changes from BioC for inclusion in development branch

0.99.4: A number of changes to the vignette

0.99.3: orderCells() now accepts the root_cell_name argument to specify the root of the ordering tree.

0.99.3: various fixes to accomodate the BioConductor build system and coding standards. Thanks to Sonali Arora for help
             with this.

0.99.0: INITIAL RELEASE

mosaics
-------

Changes in version 1.99.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Typos in the vignette are fixed.

Changes in version 1.99.3:

BUG FIXES

- mosaicsFit(): mismatch between function and help fixed.

Changes in version 1.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- mosaicsFit(): Sets bgEst="rMOM" as default.

- Help documents are polished and updated.

- Vignette is updated.

Changes in version 1.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- mosaicsFitHMM() & mosaicsPeakHMM(): Hidden-Markov-Model-based MOSAiCS model fitting & peak calling, respectively, to
  identify broad peaks such as histone modifications.

- Add new class 'MosaicsHMM' with methods show(), plot(), & estimates().

- mosaicsFit(): Introduces a new argument 'trans'.

- mosaicsFit(): Stability & robustness of model fitting were improved.

- Polish help documents of constructBins(), generateWig(), and mosaicsRunAll().

- Tested to work with >= R 3.0 properly.

BUG FIXES

- constructBins() & export(): Use correct base for BED file (one base shift).

- Reflect the changes in Rcpp packages that mosaics package depends on.

motifStack
----------

Changes in version 1.9.9:

NEW FEATURES

- Add distance axis for plotMotifStackWithRadiaPhylog

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.9.8:

NEW FEATURES

- No changes classified as 'new features' (package under active development)

BUG FIXES

- fix the bug draw motif outside of canvas when plotMotifStackWithRadiaPhylog

Changes in version 1.9.7:

NEW FEATURES

- Change the position of motifs index when plotMotifStackWithRadialPhylog.

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.9.6:

NEW FEATURES

- Increase the size of motifs when plotMotifStackWithRadialPhylog.

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.9.5:

NEW FEATURES

- Change default of plotIndex to FALSE when plotMotifStackWithRadialPhylog.

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.9.4:

NEW FEATURES

- Restyle index number when plotMotifStackWithRadialPhylog.

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.9.3:

NEW FEATURES

- Add index number when plotMotifStackWithRadialPhylog.

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.9.2:

NEW FEATURES

- No changes classified as 'new features' (package under active development)

BUG FIXES

- Fix bugs using motifDB in tests.

MSGFgui
-------

Changes in version 0.99.6:

- Fixed a visual glitch in the settings modal box

Changes in version 0.99.5:

- Yet another refactoring of shiny server functionality following Dan's (Tenenbaum) suggestions

Changes in version 0.99.4:

- Fix for namespace clash with new version of mzR

Changes in version 0.99.3:

- Refactoring of shiny part of code

- Change how mzR is used. Now connections are opened and closed whenever raw data is needed

- Fix to namespace

- Added warning to RStudio users about potential problems due to mzR problems

Changes in version 0.99.1:

- Small changes to documentation, vignette and tests

Changes in version 0.99.0:

- Submission to Bioconductor

MSGFplus
--------

Changes in version 0.99.3:

- Fix in documentation

Changes in version 0.99.2:

- Added introduction to vignette

- Added sessionInfo to vignette

- Removed calls to library(MSGFplus) in tests

Changes in version 0.99.1:

- Update vignette styling and code

Changes in version 0.99.0:

- Submission to Bioconductor

- Update MS-GF+ to v10072

MSnbase
-------

Changes in version 1.13.16:

- width generic if non-existant [2014-09-03 Wed]

- fix undocumented S4 methods warnings [2014-09-27 Sat]

Changes in version 1.13.15:

- msMap(MSmap)<- method [2014-08-12 Tue]

- Typo in MIAPE man [2014-08-29 Fri]

Changes in version 1.13.14:

- Fix xic example [2014-08-08 Fri]

- importing lattice, ggplot2 (was depends) and suggesting rgl [2014-08-09 Sat]

- MSmap infrastructure [2014-08-09 Sat]

Changes in version 1.13.13:

- fix issue with readMSnSet2 without fdata [2014-07-28 Mon]

Changes in version 1.13.12:

- Don't import width from IRanges [2014-07-23 Wed]

- qual slot is populated again [2014-07-23 Wed]

Changes in version 1.13.11:

- remove Vennerable::Venn from example and DESCRIPTION [2014-06-19 Thu]

Changes in version 1.13.10:

- new averageMSnSet function to generate an average over a list of MSnSets. [2014-06-17 Tue]

- new non-parametric coefficent of variation function [2014-06-17 Tue]

- Using/importing IRanges::width [2014-06-18 Wed]

Changes in version 1.13.9:

- new listOf helper function [2014-06-16 Mon]

- compfnames to compare and document differences in MSnSet feature names [2014-06-16 Mon]

- suggesting Vennerable (compfnames example) and roxygen2 (generate rd) [2014-06-16 Mon]

Changes in version 1.13.8:

- Add recommended biocView [2014-06-05 Thu]

Changes in version 1.13.7:

- Bug tracking link [2014-05-26 Mon]

Changes in version 1.13.6:

- new isEmpty method for Spectrum instances [2014-05-14 Wed]

- removePeaks does not try for empty Spectra [2014-05-14 Wed]

- isEmpty unit test [2014-05-14 Wed]

Changes in version 1.13.5:

- export .get.amino.acids function [2014-05-08 Thu]

- removePeaks for centroided data [2014-05-08 Thu]

- add new exported function get.atomic.mass [2014-05-08 Thu]

- Spectrum[2] prototype sets centroided=FALSE by default (instead of logical()) [2014-05-08 Thu]

- fnamesIn also supports y = "data.frame" [2014-05-14 Wed]

- Using BiocParallel for parallel support; replaced parallel argument with BPPARAM [2014-05-14 Wed]

Changes in version 1.13.4:

- Document difference between traceable and non-traceable FeaturesOfInterest instances [2014-04-30 Wed]

- ignoring desciptions with length > 1 [2014-04-30 Wed]

Changes in version 1.13.3:

- FeaturesOfInterest [2014-04-29 Tue]

Changes in version 1.13.2:

- add compareSpectra to compare Spectrum objects [2014-04-09 Wed]

- add bin method for Spectrum objects [2014-04-09 Wed]

- recreate inst/extdata/msx.rda for R 3.1 and Biobase 2.24 [2014-04-13 Sun]

- add plot,Spectrum,Spectrum method [2014-04-14 Mon]

- add calculateFragments,character,Spectrum and calculateFragments,character,missing methods [2014-04-15 Tue]

- updated readMSData test unit to ignore object@.__classVersion__ that fails with latest R/Bioc versions [2014-04-15 Tue]

- formatRt conversion from 'mm:sec' to sec and unit test [2014-04-17 Thu]

- Update nprot/npsm.prot/npsm.pep/npep.prot feature variables and relevant man/tests/argument defaults [2014-04-17 Thu]

Changes in version 1.13.1:

- add precursor method to Spectrum2 normalisation method [2014-04-08 Tue]

- add smooth for MSnExp and Spectrum classes [2014-04-10 Thu]

- add pickPeaks for MSnExp and Spectrum classes [2014-04-10 Thu]

- updated show,MSnExp to display only first/last files when > 2 [2014-04-11 Fri]

Changes in version 1.13.0:

- New devel version for Bioc 3.0

MSnID
-----

Changes in version 0.99.3:

- Rcpp back into "Depends"

Changes in version 0.99.2:

- Updated Rcpp dependency

- Updated contact info to match the bioc-devel mailing list

Changes in version 0.99.1:

- Updated dependency on mzID version

Changes in version 0.99.0:

- Minor bug fixes

- Added unit tests for filter optimization and assessment of non-tryptic and missed cleavages

- Version bump for Bioconductor submission

Changes in version 0.1.0:

- First release version

MultiMed
--------

Changes in version 0.99.1:

- MultiMed - a package for testing multiple biological mediators simultaneosly - is now available!

mvGST
-----

Changes in version 0.99.3:

OTHER CHANGES

- Removed suppressWarnings() from plot calls of interactiveGraph(), and made arrowheads be 'none' -- to avoid trivial warnings
  about zero-length arrows of indeterminate angle

- Cleaned up usage of switch() function calls in graphCell() and generateGeneSets() functions

Changes in version 0.99.2:

BUG FIXES

- Added section to p.adjust.SFL to have better error catching and restrict GO terms to names in rawp (to match Focus Level
  method)

OTHER CHANGES

- Added suppressWarnings() to plot calls of interactiveGraph (to avoid trivial warnings about edges that are too short to be
  plotted, which can happen when the number of GO terms is large in the call to graphCell)

- Updated nested if/else statements to 'switch' function calls

- Changed many 'cat' function calls to 'message' or 'warning'

- Made 'for' loops more robust through use of 'seq_along', 'seq_len', and 'seq.int'

Changes in version 0.99.1:

OTHER CHANGES

- Modified Description field in DESCRIPTION file

- Reordered import() lines in NAMESPACE file in an attempt to avoid the following two Bioconductor warnings: Warning: multiple
  methods tables found for 'unsplit' Warning: replacing previous import by 'IRanges::unsplit' when loading 'GenomeInfoDb'
  (mvGST does not call 'unsplit' in any way)

Changes in version 0.99.0:

OTHER CHANGES

- Submitted to Bioconductor 8 August 2014

Changes in version 0.1.3:

NEW FEATURES

- Added package vignette

- Made graphCell more flexible with 'background' colors for graph elements of lesser interest.

- Made p.adjust.SFL more flexible to control FWER at a user-supplied FWER (rather than automatic 0.05).

- Implemented minsize and maxsize arguments to restrict sizes of gene sets of interest.

BUG FIXES

- N/A

OTHER CHANGES

- N/A

Changes in version 0.1.2:

NEW FEATURES

- Added sample objects in data(mvGSTsamples)

BUG FIXES

- Edited pickOut function to appropriately handle non-multivariate profiles (i.e., where results.table has only one column, or
  where only one contrast is of interest)

OTHER CHANGES

- Cleaned up package dependencies in DESCRIPTION and NAMESPACE files

- Cleaned up some .Rd files (especially examples)

Changes in version 0.1.1:

NEW FEATURES

- The default method for accounting for gene name translation issues is now 'method = 2'.

- The arguments 'contrasts' and 'gene.names' in 'profileTable' are now optional. The information may be passed to
  'profileTable' as the row and column names of the argument 'pvals'.

- The function 'pickOut' now returns a data frame with the GO descriptions and p-values for each contrasts tested as well as
  the GO ID.

- The function 'graphCell' now has an optional argument that allows the user to enter the returns from 'pickOut' as the only
  argument.

BUG FIXES

- The output from 'profileTable' for 1 dimensional profiles now displays correctly.

OTHER CHANGES

- Minor clarifications in help files.

- 'pvals' is now the first argument in 'profileTable'.

mzID
----

Changes in version 1.3.6:

- Huge restructuring in documentation

- Fixed a bug that would throw an error when an empty mzID object was flattened

- Restructuring of unit tests so they conform with 'new' guidelines and gets called during check

Changes in version 1.3.5:

- Bug fixes

Changes in version 1.3.4:

- Bug fixes

Changes in version 1.3.3:

- Added method removeDecoy, to remove decoy specific information from an mzID or mzIDCollection object thus decreasing its
  size.

- Added formal getter methods to extract specific information from mzID and mzIDCollection objects

- All attribute names are now parsed as-is despite inconsitency between schema versions. Controlling compatability is now done
  using the safeName parameter in flatten and getter methods

- Substantial cleanup in documentation

Changes in version 1.3.2:

- Fixed bug due to R's change in type.convert() that would result in precise numbers being reported as strings

Changes in version 1.3.1:

- Fixed bug that would cause error on files with no modifications in the search parameters

mzR
---

Changes in version 1.99.4:

- don't run pwiz example [2014-10-05 Sun]

Changes in version 1.99.3:

- revising unit testing and using suggested convention as described in Bioc unitTesting guidelines [2014-10-02 Thu]

Changes in version 1.99.2:

- bump to new Rcpp 0.11.3 version

Changes in version 1.99.1:

- annoucning pwiz in vignette and that it will become the default backend in Bioc 3.1 [2014-09-25 Thu]

- adding a dummy close for pwiz backend to avoid breaking code that properly closes the previous ramp backend [2014-09-21 Sun]

- Not using ListBuilder to make psms data.frame and fix segfautls, via KK [2014-09-27 Sat]

Changes in version 1.99.0:

- New pwiz backend and support for mzid, contributed by Qiang Kou as part of GSoC 2014.

- Using BiocStyle for vignette [2014-08-26 Tue]


NarrowPeaks
-----------

Changes in version 1.9.2:

NEW FEATURES

- Corrected bug in getBindingProfiles.R concerning chromosome-name processing.

- Modified wig2CSARScore.c increase the max. number of chromosomes allowed in the WIG file.

- Package title modified from "Analysis of Variation in ChIP-Seq using Functional PCA Statistics" to "Shape-based Analysis of
  Variation in ChIP-Seq using Functional PCA Statistics"

ncdfFlow
--------

Changes in version 2.11.39:

Enhancement

- update 'samples' slot of 'ncdfFlowList' class from 'character' to 'named integer' to speed up look up

- replace 'coerce' method with 'ncdfFlowList' constructor function

- update validity check of 'ncdfFlowList' class

- add 'subset.ncdfFlowList' and 'subset.ncdfFlowSet' S3 functions to subset the ncdfFlowSet/ncdfFlowList based on 'pData'

- wrap "[[" logic into c++ to speed it up

- add 'save_ncfs' and 'load_ncfs' functions to save/load a ncdfFlowSet object to/from disk.

NetPathMiner
------------

Changes in version 1.1.6:

- plotCytoscapeGML: a new function to export network plots as GML files, for Cytoscape 3.0 compatibility.

- getGetsetNetworks now can return pathway-class objects, used by graphite package, which allows fast topology-based geneset
  analyses.

- getGenesets now can export the results in GMT file format, readily parsed by most GSEA packages.

- Bug fixes on KGML signaling network construction.

Changes in version 1.1.3:

- Bioconductor stable release on the development branch

npGSEA
------

Changes in version 1.1.2:

- Added x-axis labels to plots

- Changed 'z' to 'covars' in the npGSEA function

- Added a new function 'pValues' which returns all appropriate p-values for a given npGSEAResult object

- Added a new function slot to the npGSEAResult objects, betaHats.  This vector contains the betaHats for each gene in the
  gene set.  Users can thus now see each gene's individual contribution to the test statistics.

Changes in version 1.1.1:

- Fixed a typo in stat accessor for Beta approximation

OncoSimulR
----------

Changes in version 0.99.2 (2014-07-14):

- Consistently using indentation in .Rd files.

Changes in version 0.99.1 (2014-07-14):

- Minor changes for BioConductor submission:

- extended description

- minor changes to vignette

- improved documentation of posets

Changes in version 0.99.0 (2014-06-26):

- First version submitted to BioConductor.

openCyto
--------

Changes in version 1.3.17:

Enhancements

- more robust csv parsing

- register wrapper function for quadGate.seq

- parse 'subSample' argument from 'gating_args' column in csv template to allow subsampling the data prior to gating.

- pass ... from quantileGate to quantile function

- add method for multipleFitlerResult

- support multi-dimensions in framework by deprecating x,y with channels.

bug fixes

- fix the bug that stop.at argument to gating doesn't recognize partial or full path

PAA
---

Changes in version 1.0.0:

GENERAL

- First release version.

pathview
--------

Changes in version 1.5.4:

- adjust run-time messages into 3 consistent classes: Info (on progress), Note and Warning.

- include paths.hsa data, the full list of human pathway ID/names from KEGG, as to help user specify target pathways when
  calling pathview.

- updated korg to included over 80 newly added species, such as sheep, apple, mandarin orange etc. Pathview can work with 3050
  species now.

- adjust the definitions of 7 arguments for pathview function: discrete, limit,bins, both.dirs, trans.fun, low, mid, high.
  These used to be a list of two logical elements with "gene" and "cpd" as the names. They can be vectors of two or one
  element(s) now. This makes pathview easier to use now.

- Vigette has been reformatted: add a Citation section, and some example on reading user data into R, fix a few typos.


PhenStat
--------

Changes in version 2.0.0:

NEW FEATURES

- Two new statistical methods implemented RR (Reference Range Plus) and TF (Time Fixed Effect) that can be called from
  testDataset() function with argument method set to RR and TF correspondingly.

- Additional measure of biological effect added: Percentage Change to summaryOutput() and vectorOutput functions.

- New function is added to suggest analysis paths for dataset: recommendMethod()

- Function to implement the RR method: RRTest().

- Functions to implement the TF method: TFDataset() creates dataset suitable for TF, startTFModel(), finalTFModel creates
  model and fits it.

- Changes in function summaryOutput() - additional information and clearer layout.

COMPTABILITY ISSUES

- The vectorOutput() had additional elements which have increased its length.

- The function boxplotSexGenotypeBatch() has been deprecated and replace with scatterplotSexGenotypeBatch().

- Additional argument (phenotypeThreshold) with default value 0.01 has been added to the summaryOutput() function

BUG FIXES

- Number of critical bugs fixed. These are detailed at <URL:
  https://github.com/mpi2/stats_working_group/issues?q=is%3Aissue+is%3Aclosed} GitHub software tracking tools.  >

phyloseq
--------

Changes in version 1.9.15:

BUG FIXES

- `phyloseq_to_deseq2` was adding an unnecessary pseudocount of `1` to the count matrix. No longer.

- Originally described at https://github.com/joey711/phyloseq/issues/387

Changes in version 1.9.14:

BUG FIXES

- `distance` erroneously transformed Rao distance results for method DPCoA. Now Fixed.

- Originally described at https://github.com/joey711/phyloseq/issues/390

Changes in version 1.9.13:

USER-VISIBLE CHANGES

- `distance` function now supports "wunifrac" option, for `UniFrac(..., weighted=TRUE)`

- `distance` function regexpr-matching for range of variants for weighted-UniFrac, unweighted-UniFrac method option

- `distance` function regexpr-matching for `type` argument range of alternatives

- Proposed in https://github.com/joey711/phyloseq/pull/384

Changes in version 1.9.12:

BUG FIXES

- `psmelt` function now properly handles single-OTU data

- Related to https://github.com/joey711/phyloseq/issues/338

- Builds on https://github.com/joey711/phyloseq/pull/373

Changes in version 1.9.11:

USER-VISIBLE CHANGES

- More robust `plot_ordination` behavior with clearer warning/error messages.

- Coordinates automatically checked/assigned to OTU or samples

- Attempt to calculate OTU or sample weighted-average coordinates via `vegan::wascores`, if-needed

- Species weighted-average via `vegan::wascores` supported now in phyloseq::scores.pcoa

- Related to https://github.com/joey711/phyloseq/pull/364

Changes in version 1.9.10:

USER-VISIBLE CHANGES

- Massive speed/memory improvement for UniFrac calculations (via `UniFrac` or `distance`)

- Added unit-tests for the correctness of UniFrac results (no bugs detected. results from pycogent)

- Moved all unit tests to tests/testthat as recommended by CRAN maintainers and testthat doc

Changes in version 1.9.9:

USER-VISIBLE CHANGES

- `plot_net` faster, more-flexible network plot with improved defaults

- https://github.com/joey711/phyloseq/pull/353

Changes in version 1.9.8:

USER-VISIBLE CHANGES

- `mt` includes other corrections, FDR by default.

- Resolves [Issue 59](https://github.com/joey711/phyloseq/issues/59)

Changes in version 1.9.7:

BUG FIXES

- Now requires ggplot2 version 1.0.0

- Fixes bug in which ggplot 1.0 breaks in a phyloseq vignette

- Resolves [Issue 347](https://github.com/joey711/phyloseq/issues/347)

Changes in version 1.9.6:

USER-VISIBLE CHANGES

- New `sortby` argument in `plot_richness` function.

- Sort discrete x by one or more alpha-diversity measures

- Solves [Issue 342](https://github.com/joey711/phyloseq/issues/342)

- Resolves/merges [Pull 343](https://github.com/joey711/phyloseq/pull/343)

Changes in version 1.9.5:

USER-VISIBLE CHANGES

- `microbio_me_qiime` function now handles string study number for first argument.  This is in addition to numeric study
  number, already supported.

Changes in version 1.9.4:

BUG FIXES

- `rarefy_even_depth()` function no longer enforces an orientation.

- It used to always coerce to OTU-by-sample orientation.

- Solves Issue 320 https://github.com/joey711/phyloseq/issues/320

Changes in version 1.9.3:

USER-VISIBLE CHANGES

- Massive Revision to plot_tree()

- `plot_tree` now uses the `psmelt` function

- All covariates are available for aesthetic mapping.

- `plot_tree` substantial speed improvement

- Uses native ape-package C code for tree computation

- Efficient `data.table` consolidated graphic data passed to ggplot2

- Additional arguments: `treetheme` and `justify`

- `tree_layout` - new, user-accessible function

- for building alternative trees from phyloseq data.

- Foundation for solving Issue 313 and Issue 331

- https://github.com/joey711/phyloseq/issues/313

- https://github.com/joey711/phyloseq/issues/331

Changes in version 1.9.2:

BUG FIXES

- Large files cause import_usearch_uc() to have error.  Error in paste0(readLines(ucfile), collapse = "\n") : result would
  exceed 2^31-1 bytes

- Solves Issue 327: https://github.com/joey711/phyloseq/issues/327

Changes in version 1.9.1:

BUG FIXES

- Bug in `psmelt` causing unnecessary error for phyloseq datasets with empty components.

- Solves Issue 319: https://github.com/joey711/phyloseq/issues/319

piano
-----

Changes in version 1.6.0:

NEW FEATURES

- Added a function writeFilesForKiwi() that enables a seamless integration of the output from a gene set analysis with piano
  to the network-based visualization offered by the python tool Kiwi.


plethy
------

Changes in version 1.3.1:

NEW FEATURES

- Added the 'retrieveMatrix' method

- Added the 'tsplot' method

- Added the 'mvtsplot' method

SIGNIFICANT USER VISIBLE CHANGES

- Removed dependency on the 'batch' package

- Exported 'proc.sanity'

- Added dependency on ggplot2

plgem
-----

Changes in version 1.37.1:

- minor changes: --added biocViews `GeneExpression' and `MassSpectrometry' --moved Biobase and MASS from `Depends' to
  `Imports'

polyester
---------

Changes in version 9.1:

- fixed bug where read names were incorrect in output fasta files

- updated vignette

- added biocViews

- fixed coding style Package submission, version 0.99.0, 18 July 2014

prebs
-----

Changes in version 1.5.3:

- Introduced new parameter paired_ended_reads. When the data contains paired-ended reads, this parameter should be set to
  TRUE. Otherwise, the two read mates will be treated as independent units.

- Introduced new parameter ignore_strand. If set to TRUE then the strand from which read comes is ignored when counting
  overlaps. If you use strand-specific RNA-seq protocol, you should set it to FALSE, otherwise set it to TRUE.


pRoloc
------

Changes in version 1.5.19:

- added Video tag in DESCRIPTION [2014-10-07 Tue]

Changes in version 1.5.18:

- HUPO 2014 poster [2014-01-02 Thu]

Changes in version 1.5.17:

- fix 'replacing previous import by MLInterfaces::plot when loading pRoloc' warning by using specific importFrom [2014-09-27
  Sat]

Changes in version 1.5.16:

- new pRolocGUI section [2014-08-15 Fri]

- new foi section [2014-08-16 Sat]

Changes in version 1.5.15:

- svmOpt sigma defaults changed from 10^(-2:3) to 10^(-3:2) [2014-08-15 Fri]

- in xxxOptimisation, the best parameter(s) for the validation classification runs are now chosen at random instead of using
  the first best param (see change in pRoloc:::getBestParam that got a sample argument defaulted to TRUE) [2014-08-15 Fri]

- When calculating macroF1 scores (xval and validation), NAs are set to 0 (via MLInterfaces:::.macroF1(..., naAs0 = TRUE)).
  The macro F1 will not be NA (when mean of F1s is calculated) but lowered. This avoids having an NA macro F1 when 1 (or more)
  classe(s) end(s) up with NA (also set to 0) precision(s) or recall(s) [2014-08-15 Fri]

Changes in version 1.5.14:

- add title to plotDist figures [2014-08-13 Wed]

Changes in version 1.5.13:

- none

Changes in version 1.5.12:

- support mirrorX/mirrorY in highlightOnPlot [2014-07-22 Tue]

Changes in version 1.5.11:

- Remove plot2D outliers param [2014-07-10 Thu]

- fix subsetAsDataFrame for keepColNames and write unit test [2014-07-11 Fri]

Changes in version 1.5.10:

- alias lopims4 function [2014-07-01 Tue]

Changes in version 1.5.9:

- Export single steps of lopims [2014-06-23 Mon]

- Rephrase classifier parameter optimisation proceudres [2014-06-23 Mon]

Changes in version 1.5.8:

- nndistx_matrix function added to allow use of query matrix when calculating knn distances [2014-06-18 Wed]

- nndist[x]_[matrix|msnset] are now available using the experted nndist method [2014-06-18 Wed]

- markerSet and unknownSet renamed to markerMSnSet and unknownMSnSet [2014-06-19 Thu]

- functions sampleMSnSet and testMSnSet added [2014-06-19 Thu]

- fix keepColNames in pRoloc:::subsetAsDataFrame - fcol was always renamed to "markers" [2014-06-19 Thu]

Changes in version 1.5.7:

- add recommended biocView [2014-06-05 Thu]

Changes in version 1.5.6:

- addMarkers has a new mcol argument to set the markers feature variable label [2014-05-29 Thu]

Changes in version 1.5.5:

- replaced MSVBAR::rmultnorm with mvtnorm::rmvnorm since the former has been removed from CRAN and don't import [2014-05-21
  Wed]

- Bug tracking [2014-05-26 Mon]

Changes in version 1.5.4:

- testMarkers gets an error argument [2014-05-14 Wed]

- plotDist now has a ylim argument [2014-05-21 Wed]

Changes in version 1.5.3:

- import all MLInterfaces [2014-04-30 Wed]

- new param optim secion in ml vignette [2014-05-05 Mon]

- various ml typos and pRolocmakers man update [2014-05-06 Tue]

Changes in version 1.5.2:

- In plotDist, ... is now passed to matlines instead of plot and has a new lty parameter [2014-04-17 Thu]

- new highlightOnPlot function, using the new features of interest infrastructure [2014-04-29 Tue]

Changes in version 1.5.1:

- new dunkley2006 pdunit object created with mclust 4.3 [2014-04-08 Tue]

- addMarker also accepts fcol and addMarkers unit test [2014-04-14 Mon]

Changes in version 1.5.0:

- Bioc devel 3.0

pRolocGUI
---------

Changes in version 0.99.12:

- added Video tag in DESCRIPTION [2014-10-07 Tue]

Changes in version 0.99.11:

- Add screenshot to README [2014-09-05 Fri]

- fix bug when features of fois are not present

- move data tab to the end

Changes in version 0.99.10:

- Updated README [2014-09-04 Thu]

- Remove old R code file [2014-09-04]

Changes in version 0.99.9:

- fix vignette error on Windows (Dan Tenenbaum) [2014-08-26 Tue]

Changes in version 0.99.8:

- selection and display of multiple features of interest in pRolocVis and pRolocComp [2014-08-13 Wed]

Changes in version 0.99.7:

- display feature meta-data instead of protein name when hovering [2014-07-27 Tue]

- better feature highlighting [2014-07-27 Tue]

- only 1 vignette [2014-07-27 Tue]

Changes in version 0.99.6:

- support mirroring of PCA plots in pRolocComp [2014-07-22 Tue]

Changes in version 0.99.5:

- Update to latest knitcitations version and fix vignette [2014-07-15 Tue]

Changes in version 0.99.4:

- add function pRolocComp [2014-06-24 Tue]

Changes in version 0.99.3:

- multiple objects can be passed to pRolocVis by using a list [2014-06-03 Tue]

- improve query search (submit and select check box) [2014-06-02 Mon]

- add unit tests for helper functions and add manual unit test [2014-06-03 Tue]

Changes in version 0.99.2:

- change access and assignment to object pRolocGUI_SearchResults in .GlobalEnv [2014-05-30 Fri]

Changes in version 0.99.1:

- fix biocViews [2014-05-27 Tue]

- misc refactoring [2014-05-27 Tue]

Changes in version 0.99.0:

- Bioc submission

proteoQC
--------

Changes in version 1.1.10:

- fixed typos in vignette

Changes in version 1.1.9:

- add function "proteinGroup"

Changes in version 1.1.8:

- update function "chargeStat"

- change file name

Changes in version 1.1.7:

- add new function for calculation of labeling efficiency

Changes in version 1.1.6:

- update vignette

Changes in version 1.1.5:

- add function "chargeStat"

Changes in version 1.1.4:

- fix several typos in vignette

Changes in version 1.1.3:

- support MS/MS searching without modification

Changes in version 1.1.2:

- change the figure of missed cleavages

- fix several typos in vignette

Changes in version 1.1.1:

- fix a minor bug

Changes in version 1.1.0:

- support mz[X]ML format file

- add identification-independent metrics

PSEA
----

Changes in version 1.0.0:

NEW FEATURES

- Initial release

qcmetrics
---------

Changes in version 1.3.1:

- using S4Vectors' metadata and metadata<- [2014-06-16 Mon]

Changes in version 1.3.0:

- new devel version, Bioc 3.0

QDNAseq
-------

Changes in version 1.2.0 (2014-10-14):

RELEASE

- Bioconductor 3.0

IMPROVEMENTS

- segmentBins() supports another transformation option besides log2: sqrt(x + 3/8), which stabilizes the variance

- plot() can skip plotting of segments and calls by specifying doSegments=FALSE and doCalls=FALSE

- exportBins() also supports BED files

BUG FIXES

- exportBins() saves base pair positions in fixed notation instead of scientific

Changes in version 1.0.5 (2014-06-13):

BUG FIXES

- fix a bug caused by package matrixStats changing madDiff() from an S4 to an S3 method in version 0.9.4, released on
  2014-05-23

Changes in version 1.0.4 (2014-05-23):

BUG FIXES

- getBinAnnotations() fixed after being broken by a change in Bitbucket

Changes in version 1.0.3 (2014-05-23):

IMPROVEMENTS

- added exportBins() for exporting to a file

- switch graphics in the vignette to PNG

Changes in version 1.0.2 (2014-05-15):

IMPROVEMENTS

- plot() honors user-specified values for xlab and xaxt

- plot() allows omission of labels for the standard deviation and the number of data points

- improve diagnostic messages

Changes in version 1.0.1 (2014-04-17):

BUG FIXES

- smoothOutlierBins() correctly ignores bins filtered out

qpgraph
-------

Changes in version 2.00:

USER VISIBLE CHANGES

- Function qpGraph() has been replaced by a class of objects called qpGraph whose constructor implements the functionality
  that the qpGraph() function was providing. However, a couple of arguments have changed and it returns an object of the new
  defined class qpGraph. Please consult its manual page for more information.

- Removed methods qpCItest() and qpEdgeNrr() taking an 'smlSet' object as input.

- Added a new (first version) vignette showing the new functionality to estimate eQTL networks from genetical genomics data.

NEW FEATURES

- New functions, eQTLnetworkEstimate(), object classes, eQTLnetwork, and corresponding methods to ease the estimation of eQTL
  networks from genetical genomics data using higher-order conditional independence as described in the recent paper by Tur,
  Roverato and Castelo (2014).

BUG FIXES

- Bugfix in qpEdgeNrr() when using arguments restrict.Q and fix.Q and input was 'data.frame' or 'ExpressionSet'.

- Bugfix in the calculation of marginal covariances with missing observations.

- Bugfix in qpAnyGraph()

Rbowtie
-------

Changes in version 1.4.1:

NEW FEATURES

- updated bowtie to version 1.0.1 (patched in 1.5.2/1.5.3 to compile on OS X 10.9)

RDAVIDWebService
----------------

Changes in version 1.3.1:

BUG FIXED

- `DAVIDGODag` bug was fixed. Now in case that pvalueCutoff value is lower than every value present an empty goDag is returned
  (Thanks to Ulrik Stervbo)

Changes in version 1.3.0:

MINOR CHANGES

- Second bump y in version x.y.z after creating 2.14 release branch.

ReactomePA
----------

Changes in version 1.9.4:

- add gseaplot function <2014-07-31, Thu>

Changes in version 1.9.3:

- import enrichMap from DOSE <2014-07-31, Thu>

Changes in version 1.9.2:

- update roxygen to version 4 <2014-06-09, Mon>

Changes in version 1.9.1:

- bug fixed of TERM2NAME.Reactome <2014-04-21, Mon>

RedeR
-----

Changes in version 1.14.0:

- Improved interface compatibility (call-back functions) for some Java platforms.

RefNet
------

Changes in version 1.1.8:

NEW FEATURES

- a simple shiny app in the inst/apps directory

regionReport
------------

Changes in version 0.99.0:

NEW FEATURES

- Preparing to submit to Bioconductor.

SIGNIFICANT USER-VISIBLE CHANGES

- Updated the vignette and the package to work with recent versions of the packages this package depends on.

- Renamed the package from derfinderReport to regionReport and generateReport() to derfinderReport(). In the future we will
  add another report for a general GRanges object.

- Simplified derfinderReport()'s call by using advanced arguments.

- Added Travis integration.



Rgraphviz
---------

Changes in version 2.9:

- Fixed a bug in the C code of Graphviz which manifests itself when the following options are given to GCC "-Wformat
  -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2", this is the case for Ubuntu 12.04.  Reported by Vladimir
  Zhurov <vzhurov2@uwo.ca>, Venkat Seshan <veseshan@gmail.com> and Chong Tang <tangc3@unr.edu>.

RGSEA
-----

Changes in version 0.99.2:

BUG FIXES

- correct NEWS file format

Changes in version 0.99.1:

NEW FEATURES

- added new functions: RGSEAsd() RGSEApredict() RGSEAfix()

- added a vignette to explain functionality of RGSEA

- add unit tests for RGSEA.

rhdf5
-----

Changes in version 2.10.0:

NEW FEATURES

- Added support for HDF5 property lists.

- Added property list arguments to H5Dcreate and H5Dopen.

- New function h5readAttributes implemented that reads all HDF5 attributes of one object.

- New function h5version implemented.

- fillValue parameter added to h5createDataset.

- New low level general library functions H5Lcreate_external, H5Fis_hdf5, H5Fget_filesize, H5Fget_name, H5Pcreate, H5Pcopy,
  H5Pget_class, H5Pclose, H5Pclose_class, H5Pset_char_encoding, H5Pset_create_intermediate_group, H5Pset_chunk_cache,
  H5Pset_layout, H5Pset_chunk, H5Pget_chunk, H5Pset_deflate, H5Pset_fill_value, H5Pset_fill_time, H5Pset_alloc_time, H5Pequal
  implemented.

- Support for parallel Make (make -j)

USER VISIBLE CHANGES

- A warning is shown in high level function (h5read, h5write and others), if an open HDF5 handle already exists for the
  specified filename.

BUG FIXES

- Error in h5write for 0-length objects, as a consequence of automatic determining chunk size

- missing size parameter message in h5createDataset now correctly display

- checking for open file identifiers in h5read and h5ls now only searches for file names in open files, groups and datasets.

- assignment has now correct pointer target type (void *) in H5Pset_fill_value

roar
----

Changes in version 1.1.2:

NEW FEATURES

- New code to deal with experimental setups with a natural pairing between treatment and control samples. The user have to
  define the pairings and the resulting p-values of the independent tests are combined with the Fisher method. A method that
  corrects these p-values (and also those obtained in single samples settings) considering the multiple testing issue has also
  been added.

BUG FIXES

- No changes classified as 'bug fixes'

rols
----

Changes in version 1.7.2:

- Removing (temporarily) allIds("GO")) example as currently returns illegal message to fix testing error [2014-10-07 Tue]

- Updating roxygem inline docs to fix errors in generated Rds [2014-10-07 Tue]

Changes in version 1.7.1:

- add utils to Imports [2014-04-28 Mon]

Changes in version 1.7.0:

- new devel, Bioc 3.0

RRHO
----

1.6.0: () Major documentation cleanup. () Added two sided hypotheses to RRHO() with signed pvalue plotting. () Added
             comparison of three lists using RRHOCOmparison(). () Added an option for log10 pvalues in RRHO. () Added
             error handling to RRHO so that failed plotting does not crash the computation.

Rsamtools
---------

Changes in version 1.17.0:

NEW FEATURES

- pileup visits entire file if no 'which' argument specified for 'ScanBamParam' parameter of pileup. Buffered functionality
  with 'yieldSize' available to manage memory consumption when working with large BAM files

- pileup 'read_pos_breaks' parameter renamed to 'cycle_bins': cycle_bins allows users to differentiate pileup counts based on
  user-defined regions within a read.

- pileup uses PileupParam and ScanBamParam instances to calculate pileup statistics for a BAM file; returns a data.frame with
  columns summarizing information extracted from alignments overlapping each genomic position

- scanBam,BamSampler-method returns requested and actual yieldSize, and total reads

- seqinfo,BamFileList-method returns the merged seqinfo of each BamFile; seqlevels and seqlengths behave similarly.

- scanBamHeader accepts a 'what' argument to control input of the targets and / or text portion of the header, and is much
  faster for BAM files with many rnames.

SIGNIFICANT USER-VISIBLE CHANGES

- rename PileupParam class and constructor -> ApplyPileupsParam

- seqinfo,BamFile-method orders levels as they occur in the file, reverting a change introduced in Rsamtools version 1.15.28
  (version 1.17.16).

BUG FIXES

- scanBam(BamSampler(), param=param) with a 'which' argument no longer mangles element names, and respects yield size

- applyPileups checks that seqlevels are identical across files

- scanFa documentation incorrectly indicated that end coordinates beyond the range of the sequence would be truncated; they
  are an error.

- applyPileups would fail on cigars with insertion followed by reference skip, e.g., 2I1024N98M (bug report of Dan Gatti).

Rsubread
--------

Changes in version 1.16.0:

NEW FEATURES

- Subread aligner (align function) can accurately map micro RNA (miRNA) sequencing reads. A full index without gaps should be
  built for the reference genome before mapping miRNA-seq reads.

- Subjunc requires the number of consensus subreads to be at least 30 percent of the total number of extracted subreads when
  reporting hits for exonic reads. This improves the mapping accuracy of such reads.

- Reads are allowed to be extended in featureCounts.

- Minimum required number of overlapped bases can be specified when assigning reads to features in featureCounts.

- By default, Subread and Subjunc aligners do not allow more than three mismatches in the reported alignments. This however
  can be tuned via the `maxMismatches' parameter.

- A number of bug fixes.

rtracklayer
-----------

Changes in version 1.26:

NEW FEATURES

- ucscGenomes() retrieves organism information

- New function exportToTabix() exports a GRanges to a tabix-indexed tab separated file that contains all of the metadata
  columns. Use import,TabixFile to load specific ranges of data back into a GRanges.

- BigWig import/export to/from Integer/Numeric/RleList is now much more efficient, and uses a more efficient storage format
  within the BigWig file, when possible.

SIGNIFICANT USER-VISIBLE CHANGES

- BSgenome export methods are now in BSgenome.

BUG FIXES

- Handling of quotes in GFF3 is now consistent with the spec.


rTRMui
------

Changes in version 1.4:

- Improvements to user interface and general look and feel. Fontawesome icons are used now. Fixed behavior of some web
  controls (sliders). Updated tutorial. Data tables are paginated (and the number of entries per page customizable).

- Target transcription factor selection process improved. Now all possible transcription factors are added to a searchable
  dropdown menu.

RUVSeq
------

Changes in version 0.1:

- Created S4 methods for RUV* functions

- Created Vignette

- Initial submission to Bioconductor

S4Vectors
---------

CHANGES IN VERSION 0.4.0
------------------------

NEW FEATURES

    o Add isSorted() and isStrictlySorted() generics, plus some methods.

    o Add low-level wmsg() helper for formatting error/warning messages.

    o Add pc() function for parallel c() of list-like objects.

    o Add coerce,Vector,DataFrame; just adds any mcols as columns on top of the
      coerce,ANY,DataFrame behavior.

    o [[ on a List object now accepts a numeric- or character-Rle of length 1.

    o Add "droplevels" methods for Rle, List, and DataFrame objects.

    o Add table,DataTable and transform,DataTable methods.

    o Add prototype of a better all.equals() for S4 objects.

SIGNIFICANT USER-VISIBLE CHANGES

    o Move Annotated, DataTable, Vector, Hits, Rle, List, SimpleList, and
      DataFrame classes from the IRanges package.

    o Move isConstant(), classNameForDisplay(), and low-level argument
      checking helpers isSingleNumber(), isSingleString(), etc... from the
      IRanges package.

    o Add as.data.frame,List method and remove other inconsistent and not
      needed anymore "as.data.frame" methods for List subclasses.

    o Remove useless and thus probably never used aggregate,DataTable method
      that followed the time-series API.

    o coerce,ANY,List method now propagates the names.

BUG FIXES

    o Fix bug in coercion from list to SimpleList when the list contains
      matrices and arrays.

    o Fix subset() on a zero column DataFrame.

    o Fix rendering of Date/time classes as DataFrame columns.




SANTA
-----

Changes in version 2.1.0 (2014-09-12):

NOTES: 

- CITATION file added.

sapFinder
---------

Changes in version 1.3.0:

NEW FEATURES

- *Add the support of Mascot.


SCAN.UPC
--------

Changes in version 2.7.1:

FIXES

- ParseMetaFromGtfFile improvements (ambiguous function call, spaces in FASTA descriptors)

- BrainArray summarization when probe maps to multiple genes/probesets

NEW FEATURES

- This package provides support for the Affymetrix HTA 2.0 arrays

SemDist
-------

Changes in version 0.99.0:

- This is the first version of the package. Currently there are no updates to report. Changes to future versions will be
  listed here.

SeqArray
--------

Changes in version 1.5.2:

- fix the error in haploid genotypes (Y chromosome)

Changes in version 1.5.1:

- fix a bug in 'seqVCF2GDS' when the values in the FILTER column are all missing

- enhance 'seqVCF.Header'

- support the LinkingTo mechanism

seqplots
--------

Changes in version 0.99.4:

PACKAGE

- plotHeatmap function returns cluster report as GRanges structure

- redundant parameters removed from plotting functions

- plotHeatmap function have "embed" parameter for plotHeatmap - allows to plot 1st heatmap without using grid system, intended
  to use with complex plots

BUGFIX

- motif plot orientation properly dependents on strand

- GUI - reordering the heatmap respects previously set include/exclude parameters

Changes in version 0.99.1:

PACKAGE

- heatmap plotting function returns cluster report ad data.table

- getPlotSetArray function have "verbose" parameter that controls messages and warnings output

- references added to documentation

BUGFIX

- plotHeatmap and plotAverage generic methods for SeqPlots-classes respect the parameters

- package passes tests and check on 32bit Windows (plotting only, because no rtrackalyer::BigWigFile support for Win32)

Changes in version 0.99:

GENERAL

- Anchored plots and heatmaps uses [downstream]--0--0--[upstream] X-axis coordinate system instead
  [downstream]--0--[anchored]-[upstream+anchored]

PACKAGE

- package really on reference class system including MotifSetup, PlotSetArray, PlotSetList and PlotSetPair

- generic subset and data manipulation methods for SeqPlots-classes including '[', "[[" and "unlist", which allows to switch
  between classes

- automatic tests for class system, calculations and plotting functions

- documentation for all functions and classes

- PDF vignette engine replaced by HTML one

- QuickStart vignette added

GUI

- automated GUI tests using Rselenium package

BUGFIX

- issue #1: some server instances loads empty .Rdata file on startup

seqTools
--------

CHANGES IN VERSION 0.99.42
-------------------------

BUG FIXES

    o trimFastq also includes record identifiers in output files
        (=read title), e.g. SRR014849.1 EIXKN4201CFU84 length=93.


CHANGES IN VERSION 0.99.40
-------------------------

NEW FEATURES

    o Included this NEWS file



ShortRead
---------

Changes in version 1.23:

NEW FEATURES

- alphabetScore,PhredQuality-method implemented

- reverse, reverseComplement methods for ShortReadQ objects

- srlist, to access SRList data as a base R list.

SIGNIFICANT USER-VISIBLE CHANGES

- readFastq qualityType="Auto" chooses base-64 encoding when no characters are encoded at less than 59, and some are encoded
  at greater than 74.

BUG FIXES

- report() prints adapter contaminants correctly when user has stringsAsFactors=FALSE

- qa(..., sample=FALSE) no longer tries to re-match 'pattern' argument

SNPRelate
---------

Changes in version 0.99.2:

- an option to create an integer snp.id when converting from PLINK

- a new function 'snpgdsFst' to estimate Fst

Changes in version 0.99.0:

- initial submission to Bioconductor

- moving from CRAN to Bioconductor

SomaticSignatures
-----------------

Changes in version 2.0.0:

- Major refractioning

- Reduce external dependenies: stringr, h5vc, h5vcData

specL
-----

Changes in version 0.99.23:

USER VISIBLE CHANGES

- added methods for specLSet class: ionlibrary, rt.input, rt.normalized

- fixed Sys.time() units in message.

USER UNVISIBLE CHANGES

- genSwathIonLib using bpmapply

Changes in version 0.99.22:

USER VISIBLE CHANGES

- specLSet plot method

Changes in version 0.99.21:

USER VISIBLE CHANGES

- specLSet class

- replace print by show and write.Spectronaut method in specL and specLSet classes

Changes in version 0.99.20:

USER VISIBLE CHANGES

- replaced mclapply by BiocParalle::bplapply

- print method

- include iRTpeptide data

ssviz
-----

Changes in version 0.99.3:

- changed pingpong to be more generic (instead of for just piRNAs)

Changes in version 0.99.2:

- added pseudo option to getCountMatrix to handle bam files that do not have count information

- changed outputs to "message" instead of "print"

Changes in version 0.99.1:

- Initial release of ssviz : small RNA-seq visualizer and analysis toolkit

STATegRa
--------

Changes in version 1.0.0:

- Initial release of the STATegRa package for multi-omics data analysis.

TargetSearch
------------

Changes in version 1.22.0:

BUG FIXES

- Fix bug in quantMatrix. make sure that selection mass IDs match the library IDs. Add an attribute to indicate that the
  quantification mass is also a correlation mass.

TEQC
----

Changes in version 3.5.4:

- bug fix in '.coverage.hist' regarding calculation of cumulative coverage fractions (affects values in the 'sensitivity.txt'
  output table and the sensitivity barplots in the multiTEQCreport)

- bug fix in 'htmlDuplicatesBarplot': also working now for paired-end data

Changes in version 3.5.3:

- bug fix in 'coverage.hist': coverage outlier values are removed from the histogram if 'outline=FALSE', but cumulative base
  fractions are now still calculated based on the complete coverage data (affects the orange line in the coverage histogram,
  the values in the 'sensitivity.txt' output table and the sensitivity barplots in the multiTEQCreport)

Changes in version 3.5.2:

- new option 'covthreshold' in 'TEQCreport' (same as in 'coverage.hist') to manipulate which coverage value should be
  highlighted by dashed lines in the coverage histogram

Changes in version 3.5.1:

- choice between jpeg, png and tiff figures in 'TEQCreport' and 'multiTEQCreport'


trackViewer
-----------

Changes in version 1.1.1:

NEW FEATURES

- addArrowMark with label

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

VariantAnnotation
-----------------

Changes in version 1.12.0:

NEW FEATURES

- allow GRanges in 'rowData' to hold user-defined metadata cols (i.e., cols other than paramRangeID, REF, ALT, etc.)

- add isSNV() family of functions

- add faster method for converting a list matrix to an array

- add 'c' method for typed Rle classes so class is preserved

- add CITATION file

- rework writeVCF(): - FORMAT and genotype fields are parsed in C - output file is written from C - chunking added for large
  VCFs

MODIFICATIONS

- add 'row.names' to readVcf()

- deprecate restrictToSNV(); replaced by isSNV() family

- remove use of seqapply()

- show info / geno headers without splitting across blocks

- use mapCoords() in predictCoding() and locateVariants()

- deprecate refLocsToLocalLocs()

- propagate strand in predictCoding()

- replace deprecated seqsplit() with splitAsList()

- ensure GT field, if present, comes first in VCF output

- modify DESCRIPTION Author and Maintainer fileds with @R

- add 'row.names' to info,VCF-method

BUG FIXES

- modify expand() to work with no 'info' fields are imported

- remove duplicate rows from .splicesites()

- fix handling of real-valued NAs in geno omatrix construction in writeVcf()

VariantFiltering
----------------

Changes in version 1.2:

USER VISIBLE CHANGES

- Added the calculation of the position of each variant within the cDNA (when it applies), adding also a new 'Transcript' tab
  in the shiny app

- Showing the number of individuals in the 'VariantFilteringParam' object, taken from the input VCF file

BUG FIXES

- Some bugfixes in the shiny app

- Replaced a .gz file in the 'extdata'f older by the corresponding .bgz file

VariantTools
------------

Changes in version 1.8.0:

NEW FEATURES

- Add callGenotypes function for annotating a set of tallies with diploid genotype quality, likelihoods, etc in a way that is
  conformant with the gVCF spec. The methods are somewhat similar to those used by the GATK UnifiedGenotyper.


xcms
----

Changes in version 1.41.1:

BUG FIXES

- fix sampclass generation from phenoData if some combinations of factors don't exist

- disable parallel code in manpages to avoid issues on BioC windows build farm machines

xps
---

VERSION xps-1.25.2

- update XPSPreProcessing.cxx to reset mask to initial values for each array (necessary for whole genome/exon arrays)

VERSION xps-1.25.1

- update XPSPreProcessing.cxx to remove -Wsequence-point warning for numsels = ++numsels

Changes in version 3.00:


yaqcaffy
--------

Changes in version 1.25.1:

- fixed typo in vignette code chunk [2014-05-08 Thu]



Packages removed from the release
=================================

The following packages are no longer in the release:

Clonality, Resourcerer, RMAPPER, virtualArray


