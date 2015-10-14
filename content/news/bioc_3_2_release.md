Bioconductors:

We are pleased to announce Bioconductor 3.2, consisting of 1104
software packages, 257 experiment data packages, and 917
up-to-date annotation packages. 

There are 80 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.2 is compatible with R 3.2,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release includes an updated Bioconductor [Amazon Machine Image][1]
and [Docker containers][2].

Visit [http://bioconductor.org][3]
for details and downloads.

[1]: http://bioconductor.org/help/bioconductor-cloud-ami/
[2]: http://bioconductor.org/help/docker/
[3]: http://bioconductor.org

Contents
--------

* Getting Started with Bioconductor 3.2
* New Software Packages
* NEWS from new and existing packages
* Packages removed from Bioconductor since the last release

Getting Started with Bioconductor 3.2
======================================

To update to or install Bioconductor 3.2:

1. Install R 3.2.  Bioconductor 3.2 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/) .

New Software Packages
=====================

There are 80 new packages in this release of Bioconductor.



ABAEnrichment - The package ABAEnrichment is designed to test for enrichment of user defined candidate genes in the set of expressed genes in different human brain regions. The core function 'aba_enrich' integrates the expression of the candidate gene set (averaged across donors) and the structural information of the brain using an ontology, both provided by the Allen Brain Atlas project. 'aba_enrich' interfaces the ontology enrichment software FUNC to perform the statistical analyses. Additional functions provided in this package like 'get_expression' and 'plot_expression' facilitate exploring the expression data.

acde - This package provides a multivariate inferential analysis method for detecting differentially expressed genes in gene expression data. It uses artificial components, close to the data's principal components but with an exact interpretation in terms of differential genetic expression, to identify differentially expressed genes while controlling the false discovery rate (FDR). The methods on this package are described in the vignette or in the article 'Multivariate Method for Inferential Identification of Differentially Expressed Genes in Gene Expression Experiments' by J. P. Acosta, L. Lopez-Kleine and S. Restrepo (2015, pending publication).

AnnotationHubData - These recipes convert a wide variety and a growing number of public bioinformatic data sets into easily-used standard Bioconductor data structures.

BBCAnalyzer - BBCAnalyzer is a package for visualizing the relative or absolute number of bases, deletions and insertions at defined positions in sequence alignment data available as bam files in comparison to the reference bases. Markers for the relative base frequencies, the mean quality of the detected bases, known mutations or polymorphisms and variants called in the data may additionally be included in the plots.

biobroom - This package contains methods for converting standard objects constructed by bioinformatics packages, especially those in Bioconductor, and converting them to tidy data. It thus serves as a complement to the broom package, and follows the same the tidy, augment, glance division of tidying methods. Tidying data makes it easy to recombine, reshape and visualize bioinformatics analyses.

caOmicsV - caOmicsV package provides methods to visualize multi-dimentional cancer genomics data including of patient information, gene expressions, DNA methylations, DNA copy number variations, and SNP/mutations in matrix layout or network layout.

CausalR - Causal Reasoning algorithms for biological networks, including predictions, scoring, p-value calculation and ranking

ChIPComp - ChIPComp detects differentially bound sharp binding sites across multiple conditions considering matching control.

CNPBayes - Bayesian hierarchical mixture models for batch effects and copy number.

CNVPanelizer - A method that allows for the use of a collection of non-matched normal tissue samples. Our approach uses a non-parametric bootstrap subsampling of the available reference samples to estimate the distribution of read counts from targeted sequencing. As inspired by random forest, this is combined with a procedure that subsamples the amplicons associated with each of the targeted genes. The obtained information allows us to reliably classify the copy number aberrations on the gene level.

DAPAR - This package contains a collection of functions for the visualisation and the statistical analysis of proteomic data.

DChIPRep - The DChIPRep package implements a methodology to assess differences between chromatin modification profiles in replicated ChIP-Seq studies as described in Chabbert et. al - http://www.dx.doi.org/10.15252/msb.20145776.

DeMAND - DEMAND predicts Drug MoA by interrogating a cell context specific regulatory network with a small number (N >= 6) of compound-induced gene expression signatures, to elucidate specific proteins whose interactions in the network is dysregulated by the compound.

destiny - Create and plot diffusion maps

DiffLogo - DiffLogo is an easy-to-use tool to visualize motif differences.

DNABarcodes - The package offers a function to create DNA barcode sets capable of correcting insertion, deletion, and substitution errors. Existing barcodes can be analysed regarding their minimal, maximal and average distances between barcodes. Finally, reads that start with a (possibly mutated) barcode can be demultiplexed, i.e., assigned to their original reference barcode.

dupRadar - Duplication rate quality control for RNA-Seq datasets.

ELMER -  ELMER is designed to use DNA methylation and gene expression from a large number of samples to infere regulatory element landscape and transcription factor network in primary tissue.

EnrichedHeatmap - Enriched heatmap is a special type of heatmap which visualizes the enrichment of genomic signals on specific target regions. Here we implement Enriched heatmap by ComplexHeatmap package. Since this type of heatmap is just a normal heatmap but with some special settings, with the functionality of ComplexHeatmap, it would be much easier to customize the heatmap as well as concatenating to a list of heatmaps to show correspondance between different data sources.

erma - Software and data to support epigenomic road map adventures.

eudysbiome - eudysbiome a package that permits to annotate the differential genera as harmful/harmless based on their ability to contribute to host diseases (as indicated in literature) or as unknown based on their ambiguous genus classification. Further, the package statistically measures the eubiotic (harmless genera increase or harmful genera decrease) or dysbiotic(harmless genera decrease or harmful genera increase) impact of a given treatment or environmental change on the (gut-intestinal, GI) microbiome in comparison to the microbiome of the reference condition.

fCI - (f-divergence Cutoff Index), is to find DEGs in the transcriptomic & proteomic data, and identify DEGs by computing the difference between the distribution of fold-changes for the control-control and remaining (non-differential) case-control gene expression ratio data. fCI provides several advantages compared to existing methods.

FindMyFriends - A framework for doing microbial comparative genomics in R. The main purpose of the package is assisting in the creation of pangenome matrices where genes from related organisms are grouped by similarity, as well as the analysis of these data. FindMyFriends provides many novel approaches to doing pangenome analysis and supports a gene grouping algorithm that scales linearly, thus making the creation of huge pangenomes feasible.

gcatest - GCAT is an association test for genome wide association studies that controls for population structure under a general class of trait. models.

GeneBreak - Recurrent breakpoint gene detection on copy number aberration profiles.

genotypeeval - Takes in a gVCF or VCF and reports metrics to assess quality of calls.

GEOsearch - GEOsearch is an extendable search engine for NCBI GEO (Gene Expression Omnibus). Instead of directly searching the term, GEOsearch can find all the gene names contained in the search term and search all the alias of the gene names simultaneously in GEO database. GEOsearch also provides other functions such as summarizing common biology keywords in the search results.

GUIDEseq - The package implements GUIDE-seq analysis workflow including functions for obtaining unique cleavage events, estimating the locations of the cleavage sites, aka, peaks, merging estimated cleavage sites from plus and minus strand, and performing off target search of the extended regions around cleavage sites.

Guitar - The package is designed for visualization of RNA-related genomic features with respect to the landmarks of RNA transcripts, i.e., transcription starting site, start codon, stop codon and transcription ending site.

hierGWAS - Testing individual SNPs, as well as arbitrarily large groups of SNPs in GWA studies, using a joint model of all SNPs. The method controls the FWER, and provides an automatic, data-driven refinement of the SNP clusters to smaller groups or single markers.

HilbertCurve - Hilbert curve is a type of space-filling curves that fold one dimensional axis into a two dimensional space, but with still keep the locality. This package aims to provide a easy and flexible way to visualize data through Hilbert curve.

iCheck - QC pipeline and data analysis tools for high-dimensional Illumina mRNA expression data.

iGC - This package is intended to identify differentially expressed genes driven by Copy Number Alterations from samples with both gene expression and CNA data.

Imetagene - This package provide a graphical user interface to the metagene package. This will allow people with minimal R experience to easily complete metagene analysis.

INSPEcT - INSPEcT (INference of Synthesis, Processing and dEgradation rates in Time-Course experiments) analyses 4sU-seq and RNA-seq time-course data in order to evaluate synthesis, processing and degradation rates and asses via modeling the rates that determines changes in mature mRNA levels.

IONiseR - IONiseR provides tools for the quality assessment of Oxford Nanopore MinION data. It extracts summary statistics from a set of fast5 files and can be used either before or after base calling.  In addition to standard summaries of the read-types produced, it provides a number of plots for visualising metrics relative to experiment run time or spatially over the surface of a flowcell.

ldblock - Define data structures for linkage disequilibrium measures in populations.

LedPred - This package aims at creating a predictive model of regulatory sequences used to score unknown sequences based on the content of DNA motifs, next-generation sequencing (NGS) peaks and signals and other numerical scores of the sequences using supervised classification. The package contains a workflow based on the support vector machine (SVM) algorithm that maps features to sequences, optimize SVM parameters and feature number and creates a model that can be stored and used to score the regulatory potential of unknown sequences.

lfa - LFA is a method for a PCA analogue on Binomial data via estimation of latent structure in the natural parameter.

LOLA - Provides functions for testing overlap of sets of genomic regions with public and custom region set (genomic ranges) databases. This make is possible to do automated enrichment analysis for genomic region sets, thus facilitating interpretation of functional genomics and epigenomics data.

MEAL - Package to integrate methylation and expression data. It can also perform methylation or expression analysis alone. Several plotting functionalities are included as well as a new region analysis based on redundancy analysis. Effect of SNPs on a region can also be estimated.

metagenomeFeatures - metagenomeFeatures was developed for use in exploring the taxonomic annotations for a marker-gene metagenomic sequence dataset. The package can be used to explore the taxonomic composition of a marker-gene database or annotated sequences from a marker-gene metagenome experiment.

metaX - The package provides a integrated pipeline for mass spectrometry-based metabolomic data analysis. It includes the stages peak detection, data preprocessing, normalization, missing value imputation, univariate statistical analysis, multivariate statistical analysis such as PCA and PLS-DA, metabolite identification, pathway analysis, power analysis, feature selection and modeling, data quality assessment.

miRcomp - Based on a large miRNA dilution study, this package provides tools to read in the raw amplification data and use these data to assess the performance of methods that estimate expression from the amplification curves.

mirIntegrator - Tools for augmenting signaling pathways to perform pathway analysis of microRNA and mRNA expression levels.

miRLAB - Provide tools exploring miRNA-mRNA relationships, including popular miRNA target prediction methods, ensemble methods that integrate individual methods, functions to get data from online resources, functions to validate the results, and functions to conduct enrichment analyses.

motifbreakR - We introduce motifbreakR, which allows the biologist to judge in the first place whether the sequence surrounding the polymorphism is a good match, and in the second place how much information is gained or lost in one allele of the polymorphism relative to another. MotifbreakR is both flexible and extensible over previous offerings; giving a choice of algorithms for interrogation of genomes with motifs from public sources that users can choose from; these are 1) a weighted-sum probability matrix, 2) log-probabilities, and 3) weighted by relative entropy. MotifbreakR can predict effects for novel or previously described variants in public databases, making it suitable for tasks beyond the scope of its original design. Lastly, it can be used to interrogate any genome curated within Bioconductor (currently there are 22).

myvariant - MyVariant.info is a comprehensive aggregation of variant annotation resources. myvariant is a wrapper for querying MyVariant.info services

NanoStringDiff - This Package utilizes a generalized linear model(GLM) of the negative binomial family to characterize count data and allows for multi-factor design. NanoStrongDiff incorporate size factors, calculated from positive controls and housekeeping controls, and background level, obtained from negative controls, in the model framework so that all the normalization information provided by NanoString nCounter Analyzer is fully utilized.

OGSA - OGSA provides a global estimate of pathway deregulation in cancer subtypes by integrating the estimates of significance for individual pathway members that have been identified by outlier analysis.

OperaMate - OperaMate is a flexible R package dealing with the data generated by PerkinElmer's Opera High Content Screening System. The functions include the data importing, normalization and quality control, hit detection and function analysis.

Oscope - Oscope is a statistical pipeline developed to identifying and recovering the base cycle profiles of oscillating genes in an unsynchronized single cell RNA-seq experiment. The Oscope pipeline includes three modules: a sine model module to search for candidate oscillator pairs; a K-medoids clustering module to cluster candidate oscillators into groups; and an extended nearest insertion module to recover the base cycle order for each oscillator group.

Path2PPI - Package to predict pathway specific protein-protein interaction (PPI) networks in target organisms for which only a view information about PPIs is available. Path2PPI uses PPIs of the pathway of interest from other well established model organisms to predict a certain pathway in the target organism. Path2PPI only depends on the sequence similarity of the involved proteins.

pathVar - This package contains the functions to find the pathways that have significantly different variability than a reference gene set. It also finds the categories from this pathway that are significant where each category is a cluster of genes. The genes are separated into clusters by their level of variability.

PGA - This package provides functions for construction of customized protein databases based on RNA-Seq data, database searching, post-processing and report generation. This kind of customized protein database includes both the reference database (such as Refseq or ENSEMBL) and the novel peptide sequences form RNA-Seq data.

Prize - The high throughput studies often produce large amounts of numerous genes and proteins of interest. While it is difficult to study and validate all of them. Analytic Hierarchy Process (AHP) offers a novel approach to narrowing down long lists of candidates by prioritizing them based on how well they meet the research goal. AHP is a mathematical technique for organizing and analyzing complex decisions where multiple criteria are involved. The technique structures problems into a hierarchy of elements, and helps to specify numerical weights representing the relative importance of each element. Numerical weight or priority derived from each element allows users to find alternatives that best suit their goal and their understanding of the problem.

Prostar - This package provides a GUI interface for DAPAR.

ProteomicsAnnotationHubData - These recipes convert a variety and a growing number of public proteomics data sets into easily-used standard Bioconductor data structures.

RareVariantVis - Genomic variants can be analyzed and visualized using many tools. Unfortunately, number of tools for global interrogation of variants is limited. Package RareVariantVis aims to present genomic variants (especially rare ones) in a global, per chromosome way. Visualization is performed in two ways - standard that outputs png figures and interactive that uses JavaScript d3 package. Interactive visualization allows to analyze trio/family data, for example in search for causative variants in rare Mendelian diseases.

rCGH - A comprehensive pipeline for analyzing and interactively visualizing genomic profiles generated through Agilent and Affymetrix microarrays. As inputs, rCGH supports Agilent dual-color Feature Extraction files (.txt), from 44 to 400K, and Affymetrix SNP6.0 and cytoScan probeset.txt, cychp.txt, and cnchp.txt files, exported from ChAS or Affymetrix Power Tools. This package takes over all the steps required for a genomic profile analysis, from reading the files to the segmentation and genes annotations, and provides several visualization functions (static or interactive) which facilitate profiles interpretation. Input files can be in compressed format, e.g. .bz2 or .gz.

RCy3 - Vizualize, analyze and explore graphs, connecting R to Cytoscape >= 3.2.1.

RiboProfiling - Starting with a BAM file, this package provides the necessary functions for quality assessment, read start position recalibration, the counting of reads on CDS, 3'UTR, and 5'UTR, plotting of count data: pairs, log fold-change, codon frequency and coverage assessment, principal component analysis on codon coverage.

rnaseqcomp - Several quantitative and visualized benchmarks for RNA-seq quantification pipelines. Two-replicate quantifications for genes, transcripts, junctions or exons by each pipeline with nessasery meta information should be organizd into numeric matrix in order to proceed the evaluation.

ropls - Latent variable modeling with Principal Component Analysis (PCA) and Partial Least Squares (PLS) are powerful methods for visualization, regression, classification, and feature selection of omics data where the number of variables exceeds the number of samples and with multicollinearity among variables. Orthogonal Partial Least Squares (OPLS) enables to separately model the variation correlated (predictive) to the factor of interest and the uncorrelated (orthogonal) variation. While performing similarly to PLS, OPLS facilitates interpretation. Successful applications of these chemometrics techniques include spectroscopic data such as Raman spectroscopy, nuclear magnetic resonance (NMR), mass spectrometry (MS) in metabolomics and proteomics, but also transcriptomics data. In addition to scores, loadings and weights plots, the package provides metrics and graphics to determine the optimal number of components (e.g. with the R2 and Q2 coefficients), check the validity of the model by permutation testing, detect outliers, and perform feature selection (e.g. with Variable Importance in Projection or regression coefficients). The package can be accessed via a user interface on the Workflow4Metabolomics.org online resource for computational metabolomics (built upon the Galaxy environment).

RTCGA - The Cancer Genome Atlas (TCGA) Data Portal provides a platform for researchers to search, download, and analyze data sets generated by TCGA. It contains clinical information, genomic characterization data, and high level sequence analysis of the tumor genomes. The key is to understand genomics to improve cancer care. RTCGA package offers download and integration of the variety and volume of TCGA data using patient barcode key, what enables easier data possession. This may have an benefcial infuence on impact on development of science and improvement of patients' treatment. Furthermore, RTCGA package transforms TCGA data to tidy form which is convenient to use.

RTCGAToolbox - Managing data from large scale projects such as The Cancer Genome Atlas (TCGA) for further analysis is an important and time consuming step for research projects. Several efforts, such as Firehose project, make TCGA pre-processed data publicly available via web services and data portals but it requires managing, downloading and preparing the data for following steps. We developed an open source and extensible R based data client for Firehose pre-processed data and demonstrated its use with sample case studies. Results showed that RTCGAToolbox could improve data management for researchers who are interested with TCGA data. In addition, it can be integrated with other analysis pipelines for following data analysis.

sbgr - R client for Seven Bridges Genomics API.

SEPA - Given single-cell RNA-seq data and true experiment time of cells or pseudo-time cell ordering, SEPA provides convenient functions for users to assign genes into different gene expression patterns such as constant, monotone increasing and increasing then decreasing. SEPA then performs GO enrichment analysis to analysis the functional roles of genes with same or similar patterns.

SICtools - This package is to find SNV/Indel differences between two bam files with near relationship in a way of pairwise comparison thourgh each base position across the genome region of interest. The difference is inferred by fisher test and euclidean distance, the input of which is the base count (A,T,G,C) in a given position and read counts for indels that span no less than 2bp on both sides of indel region.

SISPA - Sample Integrated Gene Set Analysis (SISPA) is a method designed to define sample groups with similar gene set enrichment profiles.

SNPhood - To date, thousands of single nucleotide polymorphisms (SNPs) have been found to be associated with complex traits and diseases. However, the vast majority of these disease-associated SNPs lie in the non-coding part of the genome, and are likely to affect regulatory elements, such as enhancers and promoters, rather than function of a protein. Thus, to understand the molecular mechanisms underlying genetic traits and diseases, it becomes increasingly important to study the effect of a SNP on nearby molecular traits such as chromatin environment or transcription factor (TF) binding. Towards this aim, we developed SNPhood, a user-friendly *Bioconductor* R package to investigate and visualize the local neighborhood of a set of SNPs of interest for NGS data such as chromatin marks or transcription factor binding sites from ChIP-Seq or RNA-Seq experiments. SNPhood comprises a set of easy-to-use functions to extract, normalize and summarize reads for a genomic region, perform various data quality checks, normalize read counts using additional input files, and to cluster and visualize the regions according to the binding pattern. The regions around each SNP can be binned in a user-defined fashion to allow for analysis of very broad patterns as well as a detailed investigation of specific binding shapes. Furthermore, SNPhood supports the integration with genotype information to investigate and visualize genotype-specific binding patterns. Finally, SNPhood can be employed for determining, investigating, and visualizing allele-specific binding patterns around the SNPs of interest.

subSeq - Subsampling of high throughput sequencing count data for use in experiment design and analysis.

SummarizedExperiment - The SummarizedExperiment container contains one or more assays, each represented by a matrix-like object of numeric or other mode. The rows typically represent genomic ranges of interest and the columns represent samples.

SWATH2stats - This package is intended to transform SWATH data from the OpenSWATH software into a format readable by other statistics packages while performing filtering, annotation and FDR estimation.

synlet - Select hits from synthetic lethal RNAi screen data. For example, there are two identical celllines except one gene is knocked-down in one cellline. The interest is to find genes that lead to stronger lethal effect when they are knocked-down further by siRNA. Quality control and various visualisation tools are implemented. Four different algorithms could be used to pick up the interesting hits. This package is designed based on 384 wells plates, but may apply to other platforms with proper configuration.

TarSeqQC - The package allows the representation of targeted experiment in R. This is based on current packages and incorporates functions to do a quality control over this kind of experiments and a fast exploration of the sequenced regions. An xlsx file is generated as output.

TCGAbiolinks -  The aim of TCGAbiolinks is : i) facilitate the TCGA open-access data retrieval, ii) prepare the data using the appropriate pre-processing strategies, iii) provide the means to carry out different standard analyses and iv) allow the user to download a specific version of the data and thus to easily reproduce earlier research results. In more detail, the package provides multiple methods for analysis (e.g., differential expression analysis, identifying differentially methylated regions) and methods for visualization (e.g., survival plots, volcano plots, starburst plots) in order to easily develop complete analysis pipelines.

traseR - traseR performs GWAS trait-associated SNP enrichment analyses in genomic intervals using different hypothesis testing approaches, also provides various functionalities to explore and visualize the results.

variancePartition - Quantify and interpret multiple sources and biological and technical variation in gene expression experiments.  Uses linear mixed model to quantify variation in gene expression attributable to individual, tissue, time point, or technical variables.

XBSeq - We developed a novel algorithm, XBSeq, where a statistical model was established based on the assumption that observed signals are the convolution of true expression signals and sequencing noises. The mapped reads in non-exonic regions are considered as sequencing noises, which follows a Poisson distribution. Given measureable observed and noise signals from RNA-seq data, true expression signals, assuming governed by the negative binomial distribution, can be delineated and thus the accurate detection of differential expressed genes.


NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages since the last release. The following package NEWS is available:



acde
----

Changes in version 0.99.0:

NEW FEATURES

- first submission to Bioc-devel

ADaCGH2
-------

Changes in version 2.9.3 (2015-05-30):

- Fixed Note about require(Cairo).

- Fixed failure with changes in ffbase and not exporting min.ff and
  max.ff

Changes in version 2.9.2 (2015-05-07):

- Remove wrong {} in CITATION.

Changes in version 2.9.1 (2015-04-30):

- Added the CITATION file (which was never uploaded to the svn repos)

- Added in one reference in pSegement.Rd

affxparser
----------

Changes in version 1.41.7 (2015-09-14):

- ROBUSTNESS: Explicitly importing core R functions.

Changes in version 1.41.6 (2015-07-29):

- Updated the BiocViews field of DESCRIPTION.

Changes in version 1.41.5 (2015-06-17):

- New maintainer address (in all fields).

Changes in version 1.41.4 (2015-05-26):

- New maintainer address.

Changes in version 1.41.3 (2015-05-13):

- AVAILABILITY: Removed requirement for 'GNU make'.

Changes in version 1.41.2 (2015-05-05):

- BUG FIX/ROBUSTNESS: readCelHeader() and readCel() would core dump
  R/affxparser if trying to read multi-channel CEL files (Issue #16).
  Now an error is generated instead.  Multi-channel CEL files (e.g.
  Axiom) are not supported by affxparser. Thanks to Kevin McLoughlin
  (Lawrence Livermore National Laboratory, USA) for reporting on this.

- BUG FIX/ROBUSTNESS: readCelHeader() and readCel() on corrupt CEL
  files could core dump R/affparser (Issues #13 & #15).  Now an error
  is generated instead.  Thanks to Benilton Carvalho (Universidade
  Estadual de Campinas, Sao Paulo, Brazil) and Malte Bismarck (Martin
  Luther University of Halle-Wittenberg) for reports.

Changes in version 1.41.1 (2015-04-25):

- MEMORY ERRORS: Native functions R_affx_GetCHPEntries() and
  R_affx_ReadCHP() had unbalanced PROTECT()/UNPROTECT().  Also, native
  R_affx_GetCHPGenotypingResults() had two non-PROTECT():ed usages of
  mkString().  Thanks to Tomas Kalibera at Northeastern University for
  reporting on this.

Changes in version 1.41.0 (2015-04-16):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.2 for R (>= 3.3.0).

AllelicImbalance
----------------

Changes in version 1.8.0:

NEW FEATURES

- genotype accessor requires reference and alternative allele
  information

- genotype setter requires reference allele information

- reference fraction is now accessed through fraction(ASEset,
  top.allele.criteria="ref")

- for a more robust performance unit tests are now covering the most
  crucial calculations, such as e.g. fraction, frequency, summary,
  mapbias and phase specific calculations.

annotate
--------

Changes in version 1.47:

DEFUNCT

- probesByLL is now defunct; use AnnotationDbi::select() instead.

- blastSequences supports multiple sequence queries; use
  as="data.frame" for output.

- Improve blastSequences strategy for result retrieval, querying the
  appropriate API for status every 10 seconds after initial estimated
  processing time.

AnnotationDbi
-------------

Changes in version 1.31:

NEW FEATURES and API changes

- columns() and keytypes() sort their return values.

- ls() on a Bimap option returns keys in sort()ed order, by default * *
  * 1.30.x SERIES NEWS * * *

AnnotationHub
-------------

Changes in version 2.1.26:

SIGNIFICANT USER-VISIBLE CHANGES

- seqinfo(GRanges) for all genomes supported by GenomeInfoDb now
  contain seqlengths.


AnnotationHubData
-----------------

Changes in version 1.0.0:

BUG FIXES

- ENSEMBL recipes discover gtf files on Windows.

Changes in version 0.0.214:

NEW FEATURES

- Have added vcf files from the following genome builds for humans
  "human_9606/VCF/clinical_vcf_set/", "human_9606_b141_GRCh37p13/VCF/",
  "human_9606_b142_GRCh37p13/VCF/",
  "human_9606_b142_GRCh37p13/VCF/clinical_vcf_set/"

- For each genome build, where available, the following VCF file
  formats are available a) all.vcf.gz b) all_papu.vcf.gz c)
  common_all.vcf.gz d) clinvar.vcf.gz e) clinvar_papu f)
  common_and_clinical g) common_no_known_medical_impact

- The user can refer to
  http://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/ for
  VCF file type formats

antiProfiles
------------

Changes in version 1.9.1:

- Add NEWS file

- Add code to compute expression variability measure from Alemu, et al.
  NAR 2014.

aroma.light
-----------

Changes in version 2.99.0 (2015-10-06):

- No changes.


ballgown
--------

Changes in version 2.1.1:

- bug fix in texpr, eexpr, iexpr, and gexpr methods; they no longer
  crash with single-replicate objects

- Small documentation clarifications

BaseSpaceR
----------


biobroom
--------

Changes in version 1.0.0:

- Tidying methods for DESeq2, limma, edgeR, qvalue and MSnSet objects.

BiocInstaller
-------------

Changes in version 1.20.0:

BUG FIXES

- biocLite() uses lib.loc when calling update.packages

BiocStyle
---------

Changes in version 1.8.0:

- R Markdown templates for Bioconductor HTML and PDF documents

- Suggest 'rmarkdown' as the default engine for .Rmd documents

- Simplified use with 'rmarkdown' - no need to include a separate code
  chunk calling 'BiocStyle::markdown' anymore

- Functions facilitating the inclusion of document compilation date and
  package version in the .Rmd document header

biocViews
---------

Changes in version 1.37.2:

NEW FEATURES

- Add new function recommendPackages for finding packages tagged with a
  set of biocViews.

BrowserViz
----------

Changes in version 1.9.8:

BUG FIXES

- Significant (3x) speedup.  A 5000-node, 6000-edge graph transmits to
  Cytoscape from R in about 20 seconds.


bsseq
-----

Changes in version 1.5:

- new function strandCollapse for collapsing forward and reverse strand
  data to be unstranded.

- Updated read.bismark() to support the cytosine report files; both
  formats are supported. Other minor updates (mostly internal) to
  read.bismark(). Greatly improved documentation of this function,
  paying particular attention to differences in file formats between
  versions of Bismark.

CAMERA
------

Changes in version 1.25.3:

BUG FIXES

- Fix parallel processing with SNOW

Changes in version 1.25.2:

USER VISIBLE CHANGES

- Disabled Rmpi support and usage on Windows

Changes in version 1.25.1:

BUG FIXES

- Fix the import of the igraph methods, because it also exported a
  groups function. Solves the error "Error in x$membership : $ operator
  not defined for this S4 class"


Cardinal
--------

Changes in version 1.1.0:

BUG FIXES

- Fixed bug in formatting m/z labels affecting R 3.2.2

- Removed dependency on 'fields' because 'maps' is broken on Windows


ChIPpeakAnno
------------

Changes in version 3.3.8:

- update the documentation

Changes in version 3.3.7:

FIX BUGS

- fix the problem in creating vignettes in linux platform.

Changes in version 3.3.6:

NEW FEATURE

- add new function featureAlignedSignal, featureAlignedDistribution,
  featureAlignedHeatmap, pie1

- add new dataset HOT.spots, wgEncodeTfbsV3

- update annoGR class

- update vignettes

Changes in version 3.3.5:

NEW FEATURE

- remove all the RangedData

- add annoGR class

- update vignettes

Changes in version 3.3.4:

NEW FEATURE

- toGRanges from MACS2, narrowPeak.

- calculate the p value of overlapping peaks by reagioneR

- update documentation

Changes in version 3.3.3:

NEW FEATURE

- mergePlusMinusPeaks

Changes in version 3.3.2:

NEW FEATURE

- update citation

Changes in version 3.3.1:

NEW FEATURE

- add new citation

ChIPseeker
----------

Changes in version 1.5.11:

- remove ellipsis parameter in enrichPeakOverlap function and extend it
  to support GRanges objects <2015-10-08, Thu> + see
  https://support.bioconductor.org/p/73069/

- fixed the issue,
  https://github.com/GuangchuangYu/ChIPseeker/issues/13 <2015-10-05,
  Mon>

- update GEO info, now contains >18,000 bed file information
  <2015-09-24, Thu>

Changes in version 1.5.10:

- dropAnno function, eg. drop nearest gene annotation that far from TSS
  (>10k). <2015-09-17, Thu> + see
  https://github.com/GuangchuangYu/ChIPseeker/issues/9 + add parameter
  distanceToTSS_cutoff to enrichAnnoOverlap

- use base::subset in plotDistToTSS instead of subsetting data within
  geom_bar <2015-09-17, Thu> + see
  https://github.com/hadley/ggplot2/issues/1295 + subset parameter in
  layer will be removed in next release of ggplot2.

Changes in version 1.5.9:

- bug fixed of enrichAnnoOverlap <2015-08-26, Wed>

- change parameter order.matrix to order.by in upsetplot to meet the
  change of UpSetR pkg <2015-08-26, Wed>

Changes in version 1.5.8:

- better implementation of getFirstHitIndex.  <2015-07-29, Wed> +
  contributed by Herve Pages. + see
  https://support.bioconductor.org/p/70432/#70545.

Changes in version 1.5.7:

- add vennpie parameter in upsetplot <2015-07-20, Mon>

- upsetplot function for csAnno object <2015-07-20, Mon>

Changes in version 1.5.6:

- update citation info <2015-07-09, Thu>

- BED file +1 shift for BED coordinate system start at 0 <2015-07-07,
  Tue>

Changes in version 1.5.5:

- seq2gene for linking genomic regions to genes by many-to-many
  mapping. <2015-06-29, Mon>

Changes in version 1.5.4:

- add pseudocount in enrichPeakOverlap to prevent 0 pvalue <2015-05-22,
  Fri>

Changes in version 1.5.3:

- convert the vignette from Rnw to Rmd format <2015-05-17, Sun>

Changes in version 1.5.1:

- minor bug fixed in getChrCov <2015-04-27, Mon>

ClassifyR
---------

Changes in version 1.4.0:

- Weighted voting mode that uses the distance from an observation to
  the nearest crossover point of the class densities added.

- Bartlett Test selection function included.

- New class SelectResult. rankPlot and selectionPlot can additionally
  work with lists of SelectResult objects. All feature selection
  functions now return a SelectResult object or a list of them.

- priorSelection is a new selection function for using features
  selected in a prior cross validation for a new dataset
  classification.

- New weighted voting mode, where the weight is the distance of the x
  value from the nearest crossover point of the two densities. Useful
  for predictions with skewed features.

cleanUpdTSeq
------------

Changes in version 1.7.2:

- update documentation

Changes in version 1.7.1:

BUG FIXES

- Update the author list.

clonotypeR
----------

Changes in version 1.8.0:

NEW FEATURES

- Added a “private_clonotypes” function.

clusterProfiler
---------------

Changes in version 2.3.8:

- dropGO function <2015-09-24, Thu> + see
  https://github.com/GuangchuangYu/clusterProfiler/issues/26

- use_internal_data = TRUE in enrichKEGG example to speedup compilation
  of vignette and prevent error when online is not available
  <2015-09-23, Wed>

Changes in version 2.3.7:

- bug fixed in sorting pvalue of compareClusterResult. <2015-09-16,
  Wed> For compareCluster(fun=groupGO), there is no pvalue, use Count
  for sorting

- use_internal_data= TRUE in enrichKEGG and gseKEGG demonstrated in
  vignette due to the issue
  https://github.com/GuangchuangYu/clusterProfiler/issues/20
  <2015-08-26, Wed>

Changes in version 2.3.6:

- merge_result function <2015-07-15, Wed>

- add citation of ChIPseeker <2015-07-09, Thu>

- add 'Functional analysis of NGS data' section in vignette
  <2015-06-29, Mon>

- update vignette <2015-06-24, Wed>

Changes in version 2.3.5:

- import dotplot from DOSE <2015-06-21, Sun>

Changes in version 2.3.4:

- bug fixed in getGeneSet.GO <2015-06-16, Tue>

Changes in version 2.3.3:

- convert vignette from Rnw to Rmd format <2015-05-17, Sun>

Changes in version 2.3.2:

- bug fixed of build_Anno <2015-05-07, Thu>

- add plotGOgraph function <2015-05-05, Tue>

Changes in version 2.3.1:

- remove import RDAVIDWebService <2015-04-29, Wed>

- update buildGOmap <2015-04-29, Wed>

- remove import KEGG.db <2015-04-29, Wed>

CNVPanelizer
------------

Changes in version 0.99.12:

- Add robust as option to calculate Background calculation, then
  replacing mean with median and sd with mad

Changes in version 0.99.11:

- Calculate statistics base on ratio log transformation

Changes in version 0.99.10:

- Fix plotting and single sample bugs

Changes in version 0.99.9:

- Add fix for inconsistencies

- Remove CheckSignificance documentation for function not available

- Fix for vignette plot inconsistency bootstrap count

- Increment version number

Changes in version 0.99.8:

- Small fix for bioconductor submission

Changes in version 0.99.7:

- Specificity improvements BUG FIXES

- Improve code Styling

- Add documentation

Changes in version 0.99.6:

- Initial Bioconductor release BUG FIXES

- Add required imports at NAMESPACE

- Improve code Styling

- Add documentation

- Remove unnecessary files

codelink
--------

Changes in version 1.38.0:

- Fixed error reading Codelink files with option type="Raw" or
  type="Norm".

- Now readCodelinkSet() accepts "path=" as argument to enable reading
  files from a target directory.

cogena
------

Changes in version 1.2.0 (2015-09-01):

- update exmaple dataset.

- add drug repositioning gene sets and anlaysis.



compEpiTools
------------

Changes in version 1.3.4:

- The following function is updated + GR2fasta: The function is updated
  to handle out of limit GRanges which will not be trimmed.

Changes in version 1.3.3:

- The following function is updated + GRannotate, GRangesInPromoters,
  getPromoterClass: The functions are updated to create consistency in
  how the transcripts are retreived from txdb.

Changes in version 1.3.2:

- The following function is updated + GRcoverageSummit: the BAM file
  base-coverage can now be optionally corrected based on a control BAM
  file (subtracting the control base-coverage, after normalizing both
  for library size), before looking for the maximum coverage position
  (summit).

ComplexHeatmap
--------------

Changes in version 1.5.1:

- `oncoPrint`: there are default graphics if type of alterations is
  less than two.

- `anno_*`: get rid of lazy loading


CopywriteR
----------

Changes in version 2.0.6 (2015-06-15):

- Fixed another bug resulting from checking whether all co-analyzed
  bams were mapped to the same reference genome

Changes in version 2.0.5 (2015-06-10):

- Fixed an error resulting from checking whether all co-analyzed bams
  were mapped to the same reference genome

Changes in version 2.0.4 (2015-06-05):

- Fixed an error in plotCNA leading to mixing of chromosomes during
  plotting

Changes in version 2.0.3 (2015-05-06):

- Included MAD-values in CopywriteR log file

Changes in version 2.0.2 (2015-04-24):

- Fixed the bug that under certain circumstances error messages are
  thrown that .bam files are not coordinate sorted

Changes in version 2.0.1 (2015-04-18):

- Implemented resetting the original work directory upon exit of
  CopywriteR functions

- Corrected R dependency to version 3.2

- Updated DESCRIPTION file

- Fixed bug resulting in failure to calculate loesses RELEASE (version
  2.0)

COSNet
------

Changes in version 1.3.2:

- Modified the fields title and description of the file DESCRIPTION to
  cerrect the title in the citation of the package.

- Modified the field title of the file COSNet-package.Rd USER VISIBLE
  CHANGES BUG FIXES PLANS


csaw
----

Changes in version 1.3.14:

- Added clusterFDR() function to compute the FDR for clusters of DB
  windows.

- Added checkBimodality() function to compute bimodality scores for
  regions.

- Modified normalize(), asDGEList() to allow manual specification of
  library sizes.

- Switched from normalizeCounts(), normalize() to S4 method
  normOffsets().

- Modified default parameter specification in strandedCounts(), to
  avoid errors.

- Switched to warning from error when a restricted chromosome is
  specified in extractReads().

- Modified extractReads() interface with improved support for extended
  read and paired read extraction.

- Added normalization options to filterWindows() when using control
  samples.

- Fixed bug for proportional filtering in filterWindows().

- Allowed correlateReads() to accept paired-end specification when
  extracting data.

- Added maximizeCcf() function to estimate the average fragment length.

- Added support for strand-specific overlapping in detailRanges().

- Increased the fidelity of retained information in dumped BAM file
  from dumpPE().

- Modified strand specification arguments for profileSites(), allowed
  reporting of individual profiles.

- Removed param= specification from wwhm().

- Switched to RangedSummarizedExperiment conventions for all relevant
  functions.

- Switched to mapqFilter for scanBam() when filtering on mapping
  quality.

- Added tests for previously untested functions.

- Slight updates to documentation, user's guide.

customProDB
-----------

Changes in version 1.8.2:

UPDATED FUNCTIONS

- Update functions PrepareAnnotationEnsembl.R and
  PrepareAnnotationRefseq.R due to updates in depending packages

dagLogo
-------

Changes in version 1.7.2:

NEW FEATURES

- add html documentation

- customize the x axis as the amino acid physical position

- update the documentation

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.7.1:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- better layout of the sequence logo.

DART
----

Changes in version 1.17.1:

- Some minor changes in value names in function DoDART and DoDARTCLQ.

DChIPRep
--------

Changes in version 0.99.5:

- changed qplot() calls to standard evaluation based calls using
  ggplot()

Changes in version 0.99.4:

- python script help in the vignette has been switched to static again

Changes in version 0.99.3:

- updated the python script to be executable

- python script help in the vignette is now "live"

- added negative tests for empty and flawed input in test_dataInput.R

- small fixes in the tests and the source code

Changes in version 0.99.2:

- updated general package help

- small fixes in the vignette

- updated newsfile

Changes in version 0.99.1:

- fixed typos in the documentation

- clarified dimensions of the count tables in the vignette

- added additional checks to the data import functions

Changes in version 0.99.0:

- initial devel version

ddCt
----

Changes in version 1.24.0:

- QuantStudio (LifeTechnologies) output files are now supported

- All packages that were 'depended on' are not explicitly 'imported
  from'. This has the effect that scripts depending on ddCt may
  explicitly load libraries (such as RColorBrewer and Biobase) to use
  functions therein.

- ColMap is now simplied to contain only three slots: Sample, Detector,
  and Ct.

- Class 'CSVReader' is now renamed into 'TSVReader' since it supports
  rather tab-delimited file, instead of comma-separated file

- A direct way to convert a data.frame to a InputFrame object is
  documented. See 'example(InputFrame)'.

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

DEGreport
---------

1.4.0: 07-03-2015 Lorena Pantano <lorena.pantano@gmail.com> FIX SOME
        TEXT IN VIGNETTE, AND CLEAN DEPENDS FLAG

deltaGseg
---------

Changes in version 1.9.0:

- function name change in package changepoint (multiple.mean.norm -->
  cpt.mean)

derfinder
---------

Changes in version 1.3.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Brought back the 'mc.outfile' argument for specifying the 'outfile'
  argument in SnowParam(). See more details at
  https://stat.ethz.ch/pipermail/bioc-devel/2015-May/007531.html

Changes in version 1.3.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Deprecated functions with underscores in their names in favor of
  camelCase functions. This was done to simplify the package.

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Greatly increased the speed of the p-values calculation step. See
  https://github.com/lcolladotor/derfinder/issues/29 for details.

derfinderHelper
---------------

Changes in version 1.3.1:

NEW FEATURES

- https://github.com/leekgroup/derfinderHelper is now in synch with the
  new Bioconductor mirror for it. That is:
  https://github.com/Bioconductor-mirror/derfinderHelper

derfinderPlot
-------------

Changes in version 1.3.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Dropped tMatrix() because it was a confusing plot and also lead to
  build errors.

Changes in version 1.3.3:

NEW FEATURES

- Added the function tMatrix() which uses a GRanges object that has a
  variable of interest to compute t-Statistics between bins of the
  genome. For each bin, the values of the variable of interest from the
  regions overlapping the bin are used. We use t-Statistics instead of
  correlation because not all bins will have the same number of
  regions. This type of plot is similar to interaction plots made for
  HiC data.

Changes in version 1.3.2:

NEW FEATURES

- Added the vennRegions() function to visualize how many regions
  overlap known exons, introns, intergenic regions, none of them or
  several of these groups.

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Deprecated functions with underscores in their names in favor of
  camelCase functions. This was done to simplify the package.

DESeq2
------

Changes in version 1.10.0:

- Added MLE argument to plotMA().

- Added normTransform() for simple log2(K/s + 1) transformation, where
  K is a count and s is a size factor.

- When the design contains an interaction, DESeq() will use
  betaPrior=FALSE. This makes coefficients easier to interpret.

- Independent filtering will be less greedy, using as a threshold the
  lowest quantile of the filter such that the number of rejections is
  within 1 SD from the maximum. See ?results.

- summary() and plotMA() will use 'alpha' from results().


DiffBind
--------

Changes in version 1.16.0:

- Roll up bugfixes

- dba.plotHeatmap returns binding sites in row order

diffHic
-------

Changes in version 1.1.17:

- Renamed normalize() to normOffsets().

- Added library size specification to DIList methods normOffsets(),
  asDGEList().

- Fixed bugs under pathological settings in plotPlaid(), plotDI(),
  rotPlaid(), rotDI().

- Optimized C++ code for connectCounts(), squareCounts().

- Streamlined various R utilities used throughout all functions.

- Added iter_map.py to inst/python, for iterative mapping of DNase Hi-C
  data.

- Added the neighborCounts() function, for simultaneous read counting
  and enrichment calculation.

- Added exclude for enrichedPairs(), to provide an exclusion zone in
  the local neighborhood.

- Switched default colour in rotPlaid(), plotPlaid() to black.

- Added compartmentalize() function to identify genomic compartments.

- Added domainDirections() function to help identify domains.

- Modified correctedContact() to allow distance correction and report
  factorized probabilities directly.

- Modified marginCounts() function for proper single-end-like treatment
  of Hi-C data.

- Extended clusterPairs() to merge bin pairs from multiple DILists.

- Switched to reporting ranges directly from boxPairs(), added support
  for minimum bounding box output.

- Modified consolidatePairs() to accept index vectors for greater
  modularity.

- Added reference argument for large bin pairs, in filterDirect() and
  filterTrended().

- Added filterDiag() convenience function for filtering of
  (near-)diagonal bin pairs.

- Slight change to preparePairs() diagnostic reports when dedup=FALSE,
  and for unpaired reads.

- Added option for a distance-based threshold to define invalid
  chimeras in preparePairs().

- Updated documentation, tests and user's guide.

- Added diffHic paper entry to CITATION.


DOSE
----

Changes in version 2.7.12:

- output a list of data frames by enrichMap via invisible, so that the
  graph info can be viewed by other tools, eg. Cytoscape. <2015-09-30,
  Wed>

Changes in version 2.7.11:

- check whether input geneList is sorted <2015-09-22, Tue>

- order first followed by showCategory subsetting in barplot
  <2015-09-08, Tue> + https://support.bioconductor.org/p/71917/#71939

- bug fixed in EXTID2NAME, keytype is TAIR for arabidopsis and ORF for
  malaria <2015-08-26, Wed>

Changes in version 2.7.10:

- add Giovanni Dall'Olio as contributor in author list. <2015-07-21,
  Tue>

- update NCG data to version cancergenes_4.9.0_20150720 contributed by
  Giovanni Dall'Olio. https://github.com/GuangchuangYu/DOSE/pull/8
  <2015-07-21, Tue>

- geneInCategory may simplify to vector by R in very rare case, which
  violate the assumption of list in S4 validation checking. This issue
  was fixed. <2015-07-19, Sun>

Changes in version 2.7.9:

- add citation of ChIPseeker <2015-07-09, Thu>

- add 'Disease analysis of NGS data' section in vignette <2015-06-29,
  Mon>

- convert vignette from Rnw to Rmd <2015-06-24, Wed>

Changes in version 2.7.8:

- dotplot for enrichResult <2015-06-21, Sun>

Changes in version 2.7.7:

- speed up by using sample.int instead of sample <2015-06-04, Thu>

- add seed = FALSE to control reproduciblility of gsea function. to
  make result reproducible, explicitly set seed=TRUE <2015-06-04, Thu>
  - contributed by Vlad, see
  https://github.com/GuangchuangYu/DOSE/pull/5/ - modified by
  Guangchuang

Changes in version 2.7.6:

- bug fixed in cnetplot legend contributed by Vladislav Petyuk
  <2015-05-28, Thu>

Changes in version 2.7.5:

- update vignette <2015-05-15, Fri>

Changes in version 2.7.4:

- update permutation test with pvalue = (K+1)/(N+1) instead of K/N, so
  that p value will never be 0 <2015-05-12, Tue>

Changes in version 2.7.3:

- add setType slot in gseaResult <2015-05-15, Tue>

- add universe and geneSets slots in enrichResult <2015-05-05, Tue>

Changes in version 2.7.2:

- add vertex.label.font parameter in enrichMap <2015-04-27, Mon>

Changes in version 2.7.1:

- update NCG description in enrichNCG <2015-04-22, Wed>

DSS
---

Changes in version 2.8.0:

- Add function to detect DML for general experimental design.

dupRadar
--------

Changes in version 0.99.3:

DOCUMENTATION

- First release to Bioconductor

- `NEWS` file was added.

easyRNASeq
----------

Changes in version 2.5.6:

- Ported changes from version 2.4.5 - 2.4.7

- Added a function to create the synthetic transcripts

- Deprecated functions fetchAnnotation and knowOrganisms are now
  defunct

- Export 'basename', 'seqlevels', 'seqlevels<-' and 'seqnames<-'

Changes in version 2.5.5:

- Ported changed from release version 2.4.4

Changes in version 2.5.4:

- Ported changed from release version 2.4.3

Changes in version 2.5.3:

- Ported changed from release version 2.4.1

- Adapted to the genomeIntervals API changes (change from seq_name to
  seqnames and addition of the coercion methods to GRanges and
  consort).

Changes in version 2.5.2:

- Bioc. API changes

Changes in version 2.5.1:

- Bioc. API changes

Changes in version 2.5.0:

- Bioc. Devel Version 3.2


EBImage
-------

Changes in version 4.12.0:

NEW FEATURES

- 'antialias' argument to the function 'affine' allowing to switch off
  bilinear filtering around image borders

SIGNIFICANT USER-VISIBLE CHANGES

- deprecated '...GreyScale' morphological functions; use common
  functions 'dilate', 'erode', 'opening', 'closing', 'whiteTopHat',
  'blackTopHat' and 'selfComplementaryTopHat' for filtering both binary
  and grayscale images

- 'resize' doesn't perform bilinear filtering at image borders anymore
  in order to prevent the blending of image edges with the background
  when the image is upscaled; to switch on bilinear sampling at image
  borders use the function argument 'antialias = TRUE'

- 'floodFill' and 'fillHull' preserve storage mode

PERFORMANCE IMPROVEMENTS

- all morphological operations use the efficient Urbach-Wilkinson
  algorithm (up to 3x faster compared to the previous implementation)

- 'rotate': perform lossless 90/180/270 degree rotations by disabling
  bilinear filtering

BUG FIXES

- reimplemented the Urbach-Wilkinson algorithm used to perform
  grayscale morphological transformations

- improved pixel-level accuracy of spatial linear transformations:
  'affine', 'resize', 'rotate' and 'translate'

- 'display(..., method = "raster")': displaying of single-channel color
  images

- 'drawCircle': corrected x-y offset

- 'equalize': in case of a single-valued histogram issue a warning and
  return its argument

- 'hist': accept images with 'colorMode = Color' containing less than
  three color channels

- 'image': corrected handling of image frames

- 'medianFilter': filter size check

- 'normalize': normalization of a flat image when the argument
  'separate' is set to 'FALSE'

- 'reenumerate': corrected handling of images without background

- 'stackObjects': corrected handling of blank images without any
  objects

- 'tile': reset dimnames


EBSeq
-----

Changes in version 1.9.3:

- Correct typos in GetDEResults help file.

- Include an additional method for normalization.

Changes in version 1.9.2:

- Fixed a bug which may cause error when input a matrix to the
  sizeFactors parameter

Changes in version 1.9.1:

- Added Q&A seqction in vignette to address common questions

EDASeq
------

Changes in version 2.3:

- Added function getGeneLengthAndGCContent to compute gene length and
  GC-content.

- Updated vignette.

edge
----

Changes in version 2.1.1:

- Moderated F-test has been added for likelihood ratio test

- Weights can be inputted into odp/lrt which allows it to work for
  RNA-Seq experiments with low samples

- added function apply_jackstraw

- fixed bug in build_study

edgeR
-----

Changes in version 3.12.0:

- New argument tagwise for estimateDisp(), allowing users not to
  estimate tagwise dispersions.

- estimateTrendedDisp() has more stable performance and does not return
  negative trended dispersion estimates.

- New plotMD methods for DGEList, DGEGLM, DGEExact and DGELRT objects
  to make a mean-difference plot (aka MA plot).

- readDGE() now recognizes HTSeq style meta genes.

- Remove the F-test in glmLRT().

- New argument contrast for diffSpliceDGE(), allowing users to specify
  the testing contrast.

- glmTreat() returns both logFC and unshrunk.logFC in the output table.

- New method implemented in glmTreat() to increase the power of the
  test.

- New kegga methods for DGEExact and DGELRT objects to perform KEGG
  pathway analysis of differentially expressed genes using Entrez Gene
  IDs.

- New dimnames<- methods for DGEExact and DGELRT objects.

- Bug fix to dimnames<- method for DGEGLM objects.

- User's Guide updated. Three old case studies are replaced by two new
  comprehensive case studies.

ENmix
-----

Changes in version 1.2.0:

- Removed function QCfilter

- Heavily modified function QCinfo

- Add an argument exSample to preprocessENmix


EnrichedHeatmap
---------------

Changes in version 0.99.3:

- `anno_enrich`: get rid of lazy loading

- smoothing by locfit

Changes in version 0.99.2:

- NULL can be added to the heatmap list

Changes in version 0.99.1:

- correctly use `system.file()` now

Changes in version 0.99.0:

- first release


ensembldb
---------

Changes in version 1.1.9:

BUG FIXES

- Fixed a figure placement problem that can result in an error on
  certain systems using a recent TexLive distribution.

Changes in version 1.1.6:

BUG FIXES

- Fix a bug in lengthOf that caused an error if no filter was supplied.

Changes in version 1.1.5:

NEW FEATURES

- Implemented a shiny web app to search for genes/transcripts/exons
  using annotation of an EnsDb annotation package (function
  runEnsDbApp).

Changes in version 1.1.4:

NEW FEATURES

- Added promoters method.

Changes in version 1.1.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Added method ensemblVersion that returns the Ensembl version the
  package bases on.

- Added method getGenomeFaFile that queries AnnotationHub to retrieve
  the Genome FaFile matching the Ensembl version of the EnsDb object.

Changes in version 1.1.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Added examples to the vignette for building an EnsDb using
  AnnotationHub along with the matching genomic sequence.

- Added an example for fetching the sequences of genes, transcripts and
  exons to the vignette.

BUG FIXES

- Fixed a bug in ensDbFromGRanges and ensDbFromGtf in which the genome
  build version was not set even if provided.

Changes in version 1.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- The filter argument in all functions supports now also submission of
  a filter object, not only of a list of filter objects.

ensemblVEP
----------

Changes in version 1.10.0:

NEW FEATURES

- VEPParam78 class supports Ensembl 80

epivizr
-------

Changes in version 1.6.1:

- Transition cached ranges for querying to 'GNCLists' from
  'GIntervalTree'

erccdashboard
-------------

Changes in version 1.3.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated help files

BUG FIXES AND MINOR IMPROVEMENTS

- add warning if sample1Name or sample2Name include spaces or special
  characters.

Changes in version 1.3.4:

SIGNIFICANT USER-VISIBLE CHANGES

- None

BUG FIXES AND MINOR IMPROVEMENTS

- Fixed allReps and labelReps value assignment

- Fixed bug for duplicate row names in array data

- Fixed warning for incorrect sample names

Changes in version 1.3.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated dynRangePlot - function can now plot all replicates with or
  without replicate labels

BUG FIXES AND MINOR IMPROVEMENTS

- Added warnings to initDat function for incorrect sample names and
  missing 1:1 controls from userMixFileß

- Set minimum version requirements for ggplot2 and gridExtra

Changes in version 1.3.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette is updated to include use of grid.arrange for viewing
  figures

- multiplot function is deprecated and removed from the package,
  function is replaced with grid.arrange

BUG FIXES AND MINOR IMPROVEMENTS

- updated erccROC, estLODR, and maSignal to address changes to
  tableGrob in gridExtra package

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Version incremented to 1.3.1 for development.

BUG FIXES AND MINOR IMPROVEMENTS

- Updated DESCRIPTION file, NEWS file and .gitignore file

exomePeak
---------

Changes in version 1.9.2:

- corrected fold change calculation in bltest, by Lin Zhang

Changes in version 1.9.1:

- integrated the RMT tool for extracting the combinatorial RNA
  methylome from multiple MeRIP-Seq datasets, see ?RMT

- 3 authors: Jia Meng, Lin Zhang and Lian Liu

- added another citation (Liu 2015)


FGNet
-----

Changes in version 3.4:

- Added argument 'downloadFile' to fea_david() to choose whether to
  save the analysis results to the current directory

- DAVID now requires https. This causes errors in some systems. A
  (hopefully) temporary solution is to install some certificates
  locally.  See RDAVIDWebService help:
  https://support.bioconductor.org/p/70090/#72226

FindMyFriends
-------------

Changes in version 0.99.0:

- Submission to Bioconductor

flipflop
--------

Changes in version 1.7.3:

- add the header ie the sample names to the output tables (created when
  using output.type='table')

- correct little bug concerning the SAM file when header is present
  (see samio.cpp first test over linecount replace by totalnumread)

flowcatchR
----------

Changes in version 1.4.0:

NEW FEATURES

- added details in the vignette on the dockerflow proposal, where
  flowcatchR is made available in preconfigured Docker containers

BUG FIXES

- penaltyFunctionGenerator: simmetric penalty function is now correctly
  computed

OTHER NOTES

- added LazyData: true in the DESCRIPTION

- updated travis.yml to use native R integration

- vignette: improved structure, to enhance readability and usability;
  removed rgl hook which was not used anyway

- documentation: slight updates and enhancements

FlowRepositoryR
---------------

Changes in version 1.1.1:

NEW FEATURES

- Support for dataset searching (new function: flowRep.search).

flowType
--------

Changes in version 2.7.0:

- Minor bugfix to make partitions matrix come out correctly when more
  than one threshold is used.

gdsfmt
------

Changes in version 1.6.0:

- the version number was bumped for the Bioconductor release version
  3.2


GeneNetworkBuilder
------------------

Changes in version 1.11.1:

- update the documentation

GENESIS
-------

Changes in version 2.0.0:

- Added functions for PC-Relate.  PC-Relate provides model-free
  estimation of recent genetic relatedness in general samples.  It can
  be used to estimate kinship coefficients, IBD sharing probabilities,
  and inbreeding coefficients using genome-wide SNP data.  PC-Relate
  accounts for population structure (ancestry) among sample individuals
  through the use of ancestry representative principal components (PCs)
  to provide accurate relatedness estimates due only to recent family
  (pedigree) structure.

- GENESIS now imports the package gdsfmt.


geNetClassifier
---------------

Changes in version 1.9.3:

- Updated for compatibility with igraph 1.0

Changes in version 1.9:

- Added CITATION

- GenesRanking help: How to create a GenesRanking object based on
  scores provided by another algorithm

- Several message and warning improvements

- Bugfix: 0 IQR filter now returns the original matrix


genomation
----------

Changes in version 1.1.27:

IMPROVEMENTS AND BUG FIXES

- fixed extimating ylim and plotting notches in plotMeta

Changes in version 1.1.26:

NEW FUNCTIONS AND FEATURES

- new function patternMatrix looks for k-mer or PWM matrix occurence
  over specified DNA sequences or windows in a given genome. It depends
  on BioString and seqPattern R packages.

Changes in version 1.1.25:

IMPROVEMENTS AND BUG FIXES

- faster assigning colors to the heatmap in heatMatrix and
  multiHeatMatrix functions

Changes in version 1.1.24:

IMPROVEMENTS AND BUG FIXES

- Fixed a bug in multiHeatMatrix that occured after version 1.1.19

Changes in version 1.1.23:

IMPROVEMENTS AND BUG FIXES

- Fixed a bug in plotMeta that occured after version 1.1.19

Changes in version 1.1.22:

NEW FUNCTIONS AND FEATURES

- New argument of multiHeatMatrix clust.matrix defines which matrices
  are used for clustering

Changes in version 1.1.21:

IMPROVEMENTS AND BUG FIXES

- Improve reading paired-end alignment in score matrix functions.  In
  order to avoid duplicated reads from
  GenomicAlignments:readGAlignmentPairs reads with repeated ids are
  removed (with assumption that reads have unique ids).

Changes in version 1.1.20:

IMPROVEMENTS AND BUG FIXES

- test more examples (remove \donttest{} from examples responsible for
  reading data, e.g. readBed, include into \donttest{} only those parts
  of examples that are slow, e.g. plotMeta(.,.))

Changes in version 1.1.19:

NEW FUNCTIONS AND FEATURES

- Arithmetic, indicator and logic operations as well as subsetting work
  on ScoreMatrix, ScoreMatrixBin and ScoreMatrixList objects. New
  functionality at "ScoreMatrix-class" and "ScoreMatrixList-class" are
  documented in help pages.

- Commented examples of functions are uncommented or \donttest{} is
  used.

Changes in version 1.1.18:

NEW FUNCTIONS AND FEATURES

- Integration of Travis CI

Changes in version 1.1.17:

IMPROVEMENTS AND BUG FIXES

- Bioconductor Git-SVN bridge improved rtracklayer::import() function
  that has no more 'asRangedData' arg.

Changes in version 1.1.16:

IMPROVEMENTS AND BUG FIXES

- issue with readBam() used in scoreMatrix() to read bam files. Reading
  paired-end reads was extremely slow, because of lines of code that
  were responsible for checking if both mates from pair are not counted
  twice, e.g. when mates map into two different windows of interest.
  Right now we allow it. Fixed bug in reading single-end reads.

Changes in version 1.1.15:

IMPROVEMENTS AND BUG FIXES

- fix bug of plotMeta() when dispersion=NULL and smoothfun is not NULL

Changes in version 1.1.14:

IMPROVEMENTS AND BUG FIXES

- improve checking if the ncols of matrices match in heatMeta() and
  plotMeta()

Changes in version 1.1.13:

IMPROVEMENTS AND BUG FIXES

- add new argument cex.legend to plotTargetAnnotation() to specify the
  size of the legend.

- changed ScoreMatrixBin() to run faster when noCovNA=TRUE

- check not only for .bw but also .bigWig and .bigwig extensions of
  BigWig file in ScoreMatrix() and ScoreMatrixBin()

Changes in version 1.1.12:

IMPROVEMENTS AND BUG FIXES

- gffToGRanges is now a wrapper for import from rtracklayer

Changes in version 1.1.11:

IMPROVEMENTS AND BUG FIXES

- plotMeta function bug fixed. The error occured during plotting
  meta-profile when dispersion=NULL. Implemented by Katarzyna
  Wreczycka.

Changes in version 1.1.10:

IMPROVEMENTS AND BUG FIXES

- changes in heatMatrix() and plotMeta() - dispersion, smoothfun and
  clustfun arguments changed from be default FALSE to NULL.

Changes in version 1.1.9:

IMPROVEMENTS AND BUG FIXES

- ScoreMatrix function bug fixed. The bug occured when RleList object
  ("target" argument) and GRanges object ("windows" argument) did not
  have the same chromosome ordering. In that case, the returned score
  matrix was correct but the row ordering did not correspond to the row
  ordering supplied in the windows argument.

Changes in version 1.1.8:

IMPROVEMENTS AND BUG FIXES

- changes in the vignette, updated the AnnotationHub example, thanks to
  Christopher Bottoms Implemented by Katarzyna Wręczycka.

Changes in version 1.1.7:

IMPROVEMENTS AND BUG FIXES

- extended capabilities for meta-plots; added new parameters to
  plotMeta() that enable to show dispersion interval bands around a
  central tendency (mean or median) and smoothing them.  Implemented by
  Katarzyna Wręczycka.

Changes in version 1.1.6:

IMPROVEMENTS AND BUG FIXES

- heatMatrix() and multiHeatMatrix() can visualizate score
  matrix/matrices clustered by using not only kmeans, but custom
  clustering function (clustfun argument). Implemented by Katarzyna
  Wreczycka.

Changes in version 1.1.5:

IMPROVEMENTS AND BUG FIXES

- changes in the vignette,explicit call to
  GenomicRanges::countOverlaps, thanks to Christopher Bottoms

Changes in version 1.1.4:

IMPROVEMENTS AND BUG FIXES

- scoreMatrix, scoreMatrixBin and ScoreMatrixList read paired-end BAM
  files.  Paired reads are treated as fragments. Implemented by
  Katarzyna Wręczycka.

Changes in version 1.1.3:

IMPROVEMENTS AND BUG FIXES

- Functions that read data from text files re-written from
  data.table::fread() to readr::read_delim() and now they can read
  compressed files from a URL

- Changes in readBed(), readNarrowPeak(), readBroadPeak() and
  gffToGRanges() arguments; All of them have now track.line argument
  that can be FALSE, "auto" or an integer indicating number of first
  lines to skip. The zero.based argument was added to readBed() that
  tells whether ranges in the bed file are in 0 or 1-based coordinate
  system. Implemented by Katarzyna Wręczycka.

Changes in version 1.1.2:

IMPROVEMENTS AND BUG FIXES

- scoreMatrixList() runs in parallel by using parallel::mclapply().
  Implemented by Katarzyna Wręczycka.

Changes in version 1.1.1:

NEW FUNCTIONS AND FEATURES

- We are in BioC

IMPROVEMENTS AND BUG FIXES

- readGeneric() and all the functions that read data from text files
  are faster due to use of data.table::fread(). Implemented by
  Katarzyna Wręczycka

genomeIntervals
---------------

Changes in version 1.25.3:

- Changed the "c" function from an S3method to use the S4 generic.

Changes in version 1.25.2:

- Fixed the width function to handle correctly the left end of the
  intervals

Changes in version 1.25.1:

- Changed readGff3 to use closed intervals by default. Implemented two
  sub-functions that implement reading a gff3 as base-pair features
  only (no zero length intervals, i.e. right-closed intervals) or which
  allows for zero length intervals, i.e. right-open intervals, when
  start equals end)

- Deprecated the seq_name accessors in favour of the BiocGenerics
  seqnames

- Added a width accessor - similar to the IRanges functionality

- Added coercion to GRangesList and RangedData

- Edited some of the documentation (man page) and NAMESPACE generation
  to use roxygen2

Changes in version 1.25.0:

- Bioc Version 3.2 devel

GenomicAlignments
-----------------

Changes in version 1.6.0:

NEW FEATURES

- Add strandMode() getter and setter for GAlignmentPairs objects in
  response to the following post:
  https://support.bioconductor.org/p/65844/ See ?strandMode for more
  information.

- The readGAlignment*() functions now allow repeated seqnames in the
  BAM header.

- Add "coverage" method for GAlignmentsList objects.

- The strand setter now works on a GAlignmentsList object in a
  restricted way (only strand(x) <- "+" or "-" or "*" is supported).

SIGNIFICANT USER-LEVEL CHANGES

- summarizeOverlaps() now returns a RangedSummarizedExperiment object
  (defined in the new SummarizedExperiment package) instead of an "old"
  SummarizedExperiment object (defined in the GenomicRanges package).

- Slightly modify the behavior of junctions() on a GAlignmentPairs
  object so that the returned ranges now have the "real strand" set on
  them. See ?junctions and the documentation of the 'real.strand'
  argument in the man page of GAlignmentPairs objects for more
  information.

- Add 'real.strand' argument to first() and last() getters for
  GAlignmentPairs objects.

DEPRECATED AND DEFUNCT

- Deprecate left() and right() getters and strand() setter for
  GAlignmentPairs objects.

- Deprecate 'invert.strand' argument of first() and last() getters for
  GAlignmentPairs objects.

- Deprecate 'order.as.in.query' argument of "grglist" method for
  GAlignmentPairs objects.

- Deprecate 'order.as.in.query' argument in "rglist" method for
  GAlignmentsList objects (this concept is not defined for these
  objects in general and the argument was ignored anyway).

- After being deprecated in BioC 3.1, the "mapCoords" and "pmapCoords"
  methods are now defunct. mapToAlignments() should be used instead.

- After being deprecated in BioC 3.1, the readGAlignment*FromBam()
  functions are now defunct. Everybody says "Let's all use the
  readGAlignment*() functions instead! (no FromBam suffix). Yeah!"

BUG FIXES

- Various fixes to grglist/granges/rglist/ranges methods for
  GAlignmentsList objects: - Respect cigar information (as claimed in
  man page). - Restore 'drop.D.ranges' argument in "grglist" method
  (mistakenly got deprecated at the beginning of BioC 3.2 devel cycle).
  - The 'drop.D.ranges' argument in "rglist" method now works (was
  ignored). - Handle empty list elements.

GenomicFeatures
---------------

Changes in version 1.22:

NEW FEATURES

- Add coverageByTranscript() and pcoverageByTranscript(). See
  ?coverageByTranscript for more information.

- Various improvements to makeTxDbFromGFF(): - Now supports
  'format="auto"' for auto-detection of the file format.  - Now
  supports naming features by dbxref tag (like GeneID). This has proven
  useful when importing GFFs from NCBI.

- Improvements to the coordinate mapping methods: - Support recycling
  when length(transcripts) == 1 for parallel mapping functions.  - Add
  pmapToTranscripts,Ranges,GRangesList and
  pmapFromTranscripts,Ranges,GRangesList methods.

- Adds 'taxonomyId' argument to the makeTxDbFrom*() functions.

- Improvements to makeTxDbPackage(): - Add 'pkgname' argument to
  makeTxDbPackage() to let the user override the automatic naming of
  the package to be generated.  - Support person objects for
  'maintainer' and 'author' arguments to makeTxDbPackage().

- The 'chrominfo' vector passed to makeTxDb() can now mix NAs and
  non-NAs.

SIGNIFICANT USER-VISIBLE CHANGES

- 2 significant changes to makeTxDbFromGRanges() and makeTxDbFromGFF():
  - They now also import transcripts of type pseudogenic_transcript and
  pseudogenic_exon.  - They normally get the "gene_id" and
  "[tx|exon|CDS]_name" columns from the Name tag. Now they will also
  infer these columns from the ID tag when the Name tag is missing.

- Improve handling of 'circ_seqs' argument by makeTxDbFromUCSC(),
  makeTxDbFromGFF(), and makeTxDbFromBiomart(): no more annoying
  warning when none of the strings in DEFAULT_CIRC_SEQS matches the
  seqlevels of the TxDb object to be made.

- 2 minor changes to makeTxDbFromBiomart(): - Now drops unneeded
  chromosome info when importing an incomplete transcript dataset.  -
  Now returns a TxDb object with 'Full dataset' field set to 'no' when
  makeTxDbFromBiomart() is called with user-supplied 'filters'.

- makeTxDbPackage() now includes data source in the package name by
  default (for non UCSC and BioMart databases).

- The following changes were made to the coordinate mapping methods: -
  mapToTranscripts() now reports mapped position with respect to the
  transcription start site regardless of strand.  - Change
  'ignore.strand' default from TRUE to FALSE in all coordinate mapping
  methods for consistency with other functions that already have the
  'ignore.strand' argument.  - Name matching in mapFromTranscripts() is
  now done with seqnames(x) and names(transcripts).  - The
  pmapFromTranscripts,*,GRangesList methods now return a GRangesList
  object. Also they no longer use 'UNMAPPED' seqname for unmapped
  ranges.  - Remove uneeded ellipisis from the argument list of various
  coordinate mapping methods.

- Change behavior of seqlevels0() getter so it does what it was always
  intended to do.

- The order of the transcripts returned by transcripts() has changed:
  now they are guaranteed to be in the same order as in the GRangesList
  object returned by exonsBy().

- Code improvements and speedup to the transcripts(), exons(), cds(),
  exonsBy(), and cdsBy() extractors.

- In order to avoid loss of information (and make it reversible with
  makeTxDbFromGRanges()), the "asGFF" method for TxDb objects now
  propagates the "exon_name" and "cds_name" columns to the GRanges
  object.

DEPRECATED AND DEFUNCT

- After being deprecated in BioC 3.1, the makeTranscriptDb*() functions
  are now defunct.

- After being deprecated in BioC 3.1, the 'exonRankAttributeName',
  'gffGeneIdAttributeName', 'useGenesAsTranscripts', 'gffTxName', and
  'species' arguments of makeTxDbFromGFF() are now defunct.

- Remove sortExonsByRank() (was defunct in BioC 3.1).

BUG FIXES

- Fix bug in fiveUTRsByTranscript() and threeUTRsByTranscript()
  extractors when the TxDb object had "user defined" seqlevels and/or a
  set of "active/inactive" seqlevels defined on it.

- Fix bug in isActiveSeq() setter when the TxDb object had "user
  defined" seqlevels on it.

- Fix many issues with seqlevels() setter for TxDb objects. In
  particular make the 'seqlevels(x) <- seqlevels0(x)' idiom work on
  TxDb objects.

- Fix bug in makeTxDbFromBiomart() when using it to retrieve a dataset
  that doesn't provide the cds_length attribute (e.g. sitalica_eg_gene
  dataset in plants_mart_26).

GenomicFiles
------------

Changes in version 1.6.0:

BUG FIXES

- dimnames<- correctly updates dim names

GenomicInteractions
-------------------

Changes in version 1.3.9:

NEW FEATURES

- InteractionTrack class for plotting interactions with Gviz

- plotAvgViewpoint for plotting summarised interactions around a set of
  features

- summariseByFeaturePairs: to summarise interactions between all pairs
  of two feature sets

SIGNIFICANT USER-LEVEL CHANGES

- annotateInteractions is significantly faster

- Data import has been made stricter and more consistent across
  different file formats.  Homer interaction files now have data
  imported that was previously discarded.

BUG FIXES

- single viewpoint plotting (plotViewpoint)

- is.dt, is.pt

- calculateDistances

- some import/export issues

GenomicRanges
-------------

Changes in version 1.22.0:

NEW FEATURES

- Support coercions back and forth between a GRanges object and a
  character vector (or factor) with elements in the format
  'chr1:2501-2800' or 'chr1:2501-2800:+'.

- Add facilities for manipulating "genomic variables": bindAsGRanges(),
  mcolAsRleList(), and binnedAverage(). See ?genomicvars for more
  information.

- Add "narrow" method for GRangesList objects.

- Enhancement to the GRanges() constructor. If the 'ranges' argument is
  not supplied then the constructor proceeds in 2 steps: 1. An initial
  GRanges object is created with 'as(seqnames, "GRanges")'.  2. Then
  this GRanges object is updated according to whatever other arguments
  were supplied to the call to GRanges(). Because of this enhancement,
  GRanges(x) is now equivalent to 'as(x, "GRanges")' e.g. GRanges() can
  be called directly on a character vector representing ranges, or on a
  data.frame, or on any object for which coercion to GRanges is
  supported.

- Add 'ignore.strand' argument to "range" and "reduce" methods for
  GRangesList objects.

- Add coercion from SummarizedExperiment to RangedSummarizedExperiment
  (also available via updateObject()). See 1st item in DEPRECATED AND
  DEFUNCT section below for more information about this.

- GNCList objects are now subsettable.

- "coverage" methods now accept 'shift' and 'weight' supplied as an
  Rle.

SIGNIFICANT USER-LEVEL CHANGES

- Modify behavior of "*" strand in precede() / follow() to mimic
  'ignore.strand=TRUE'.

- Revisit "pintersect" methods for GRanges#GRanges,
  GRangesList#GRanges, and GRanges#GRangesList: - Sanitize their
  semantic. - Add 'drop.nohit.ranges' argument (FALSE by default). - If
  'drop.nohit.ranges' is FALSE, the returned object now has a "hit"
  metadata column added to it to indicate the elements in 'x' that
  intersect with the corresponding element in 'y'.

- binnedAverage() now treats 'numvar' as if it was set to zero on
  genomic positions where it's not set (typically happens when 'numvar'
  doesn't span the entire chromosomes because it's missing the trailing
  zeros).

- GRanges() constructor no more mangles the names of the supplied
  metadata columns (e.g. if the column is "_tx_id").

- makeGRangesFromDataFrame() now accepts "." in strand column (treated
  as

- GNCList() constructor now propagates the metadata columns.

- Remove "seqnames" method for RangesList objects.

DEPRECATED AND DEFUNCT

- The SummarizedExperiment class defined in GenomicRanges is deprecated
  and replaced by 2 new classes defined in the new SummarizedExperiment
  package: SummarizedExperiment0 and RangedSummarizedExperiment. In
  BioC 2.3, the SummarizedExperiment class will be removed from the
  GenomicRanges package and the SummarizedExperiment0 class will be
  renamed SummarizedExperiment. To facilitate this transition, a
  coercion method was added to coerce from old SummarizedExperiment to
  new RangedSummarizedExperiment (this coercion is performed when
  calling updateObject() on an old SummarizedExperiment object).

- makeSummarizedExperimentFromExpressionSet() and related stuff was
  moved to the new SummarizedExperiment package.

- After being deprecated in BioC 3.1, the rowData accessor is now
  defunct (replaced with the rowRanges accessor).

- After being deprecated in BioC 3.1, GIntervalTree objects and the
  "intervaltree" algorithm in findOverlaps() are now defunct.

- After being deprecated in BioC 3.1, mapCoords() and pmapCoords() are
  now defunct.

BUG FIXES

- 2 tweaks to subsetting *by* an GenomicRanges: - Improve speed when
  the object to subset is a SimpleList (e.g.  SimpleRleList). - Fix
  issue when the GenomicRanges subscript is empty.

genotypeeval
------------

Changes in version 0.99.0:

- genotypeeval in use at Genentech and HLI.

- genotypeeval submitted to Bioconductor.

GEOquery
--------

2.36: New Features * New, faster SOFT format parsing (Leonardo Gama) *
        Turned on unit tests in Travis CI * Test coverage metrics added
        Bug fixes * default download method no longer assumes that curl
        is installed on linux * GSEMatrix parsing from file now finds
        cached GPLs


ggtree
------

Changes in version 1.1.20:

- fixed bug in geom_tiplab when x contains NA (eg, removing by collapse
  function) <2015-10-01, Thu>

- bug fixed in %add2%, if node available use node, otherwise use label
  <2015-09-04, Fri>

- bug fixed of subview for considering aes mapping of x and y
  <2015-09-03, Thu>

- update vignette by adding r8s example <2015-09-02, Wed>

- defined r8s class, see http://loco.biosci.arizona.edu/r8s/
  <2015-09-02, Wed> + add r8s sample files + read.r8s, parser function
  + fortify method + plot, get.tree, get.fields, groupOTU, groupClade,
  scale_color, gzoom and show methods

- bug fixed in fortify.multiPhylo, convert df$.id to factor of
  levels=names(multiPhylo_object) <2015-09-02, Wed>

- update scale_x_ggtree to support Date as x-axis <2015-09-01, Tue>

- add mrsd parameter for user to specify 'most recent sampling date'
  for time tree <2015-09-01, Tue> - remove 'time_scale' parameter.

- defined 'raxml' class for RAxML bootstrapping analysis result
  <2015-09-01, Tue> + see
  http://sco.h-its.org/exelixis/web/software/raxml/hands_on.html +
  read.raxml, parser function + plot, get.tree, get.fields, groupOTU,
  groupClade, scale_color, gzoom and show methods + fortify.raxml
  method

- bug fixed in edgeNum2nodeNum for jplace parsing jplace file
  <2015-09-01, Tue>

Changes in version 1.1.19:

- use fortify instead of fortify.phylo in fortify.multiPhylo, so that
  multiPhylo can be a list of beast/codeml or other supported objects.
  <2015-08-31, Mon>

- support multiPhylo object, should use + facet_wrap or + facet_grid
  <2015-08-31, Mon>

- remove dependency of EBImage and phytools to speedup the installation
  process of ggtree <2015-08-31, Mon> + these two packages is not
  commonly used, and will be loaded automatically when needed.

Changes in version 1.1.18:

- layout name change to 'rectangular', 'slanted', 'circular'/'fan' for
  phylogram and cladogram (if branch.length = 'none') 'unroot' is not
  changed. <2015-08-28. Fri>

- implement geom_point2, geom_text2, geom_segment2 to support
  subsetting <2015-08-28, Fri> see
  https://github.com/hadley/ggplot2/issues/1295

- update geom_tiplab according to geom_text2 and geom_segment2
  <2015-08-28, Fri>

- add geom_tippoint, geom_nodepoint and geom_rootpoint <2015-08-28,
  Fri>

Changes in version 1.1.17:

- bug fixed in rm.singleton.newick by adding support of scientific
  notation in branch length <2015-08-27, Thu>

- bug fixed in gheatmap, remove inherit aes from ggtree <2015-08-27,
  Thu>

- add 'width' parameter to add_legend, now user can specify the width
  of legend bar <2015-08-27, Thu>

- add 'colnames_position' parameter to gheatmap, now colnames can be
  display on the top of heatmap <2015-08-27, Thu>

- theme_transparent to make background transparent <2015-08-27, Thu>

- subview for adding ggplot object (subview) to another ggplot object
  (mainview) <2015-08-27, Thu>

Changes in version 1.1.16:

- update citation <2015-08-17, Mon>

Changes in version 1.1.15:

- open text angle parameter for annotation_clade/annotation_clade2
  <2015-08-13, Thu>

- support changing size of add_legend <2015-08-13, Thu>

- reroot methods for phylo and beast <2015-08-07, Fri>

Changes in version 1.1.14:

- update paml_rst to compatible with only marginal ancestral sequence
  or joint ancestral sequence available <2015-08-07, Fri>

Changes in version 1.1.13:

- implement annotation_image <2015-08-01, Sat>

- better implementation of geom_tiplab for accepting aes mapping and
  auto add align dotted line <2015-08-01, Sat>

- open group_name parameter of groupOTU/groupClade to user <2015-08-01,
  Sat>

Changes in version 1.1.12:

- update vignette according to the changes <2015-07-31, Fri>

- add mapping parameter in ggtree function <2015-07-31, Fri>

- extend groupClade to support operating on tree view <2015-07-31, Fri>

- extend groupOTU to support operating on tree view <2015-07-31, Fri>

- new implementation of groupClade & groupOTU <2015-07-31, Fri>

Changes in version 1.1.11:

- annotation_clade and annotation_clade2 functions. <2015-07-30, Thu>

- better add_legend implementation. <2015-07-30, Thu>

- add ... in theme_tree & theme_tree2 for accepting additional
  parameter. <2015-07-30, Thu>

- better geom_tree implementation. Now we can scale the tree with
  aes(color=numVar). <2015-07-30, Thu>

Changes in version 1.1.10:

- solve overlapping branches for layout = "fan" || "radial", that are
  coord_polar-ized layouts. see
  https://github.com/GuangchuangYu/ggtree/issues/6, contributed by
  Vincent Bonhomme. <2015-07-22, Wed>

Changes in version 1.1.9:

- update add_legend to align legend text <2015-07-06, Mon>

- bug fixed in internal function, getChild.df, which should not include
  root node if selected node is root <2015-07-01, Wed>

- rotate function for ratating a clade by 180 degree and update
  vignette <2015-07-01, Wed>

- get_taxa_name function will return taxa name vector of a selected
  clade <2015-06-30, Tue>

- add example of flip function in vignette <2015-06-30, Tue>

- flip function for exchanging positions of two selected branches
  <2015-06-30, Tue>

Changes in version 1.1.8:

- update get.placement <2015-06-05, Fri>

- edgeNum2nodeNum for converting edge number to node number for
  EPA/pplacer output <2015-06-04, Thu>

- mv scale_x_gheatmap to scale_x_ggtree, which also support msaplot
  <2015-06-02, Tue>

- add mask function <2015-06-02, Tue>

Changes in version 1.1.7:

- add example of msaplot in vignette <2015-05-22, Fri>

- msaplot for adding multiple sequence alignment <2015-05-22, Fri>

Changes in version 1.1.6:

- add vertical_only parameter to scaleClade and set to TRUE by default.
  only vertical will be scaled by default. <2015-05-22, Fri>

- update add_colorbar & add_legend <2015-05-21, Thu>

- add example of add_legend and gheatmap in vignette <2015-05-18, Mon>

- gheatmap implementation of gplot <2015-05-18, Mon>

- add_legend for adding evolution distance legend <2015-05-18, Mon>

Changes in version 1.1.5:

- implement scaleClade <2015-05-12, Tue>

Changes in version 1.1.4:

- better performance of parsing beast tree <2015-05-11, Mon> + support
  beast tree begin with 'tree tree_1 = ' and other forms. + support
  file that only contains one evidence for some of the nodes/tips

- update add_colorbar to auto determine the position <2015-05-04, Mon>

- add_colorbar function <2015-04-30, Thu>

Changes in version 1.1.3:

- add space between residue substitution (e.g. K123R / E155D)
  <2015-04-30, Thu>

- remove slash line in heatmap legend <2015-04-30, Thu>

- update vignette to add example of merge_tree <2015-04-29, Wed>

Changes in version 1.1.2:

- in addition to parsing beast time scale tree in XXX_year[\\.\\d]*,
  now supports XXX/year[\\.\\d]* <2015-04-29, Wed>

- add examples folder in inst that contains sample data <2015-04-29,
  Wed>

- update gplot, now rowname of heatmap will not be displayed
  <2015-04-28, Tue>

- add line break if substitution longer than 50 character <2015-04-28,
  Tue>

- support calculating branch for time scale tree <2015-04-28, Tue>

- remove parsing tip sequence from mlb and mlc file <2015-04-28, Tue>

- remove tip.fasfile in read.paml_rst for rstfile already contains tip
  sequence <2015-04-28, Tue>

- scale_color accepts user specific interval and output contains
  'scale' attribute that can be used for adding legend <2015-04-28,
  Tue>

- extend fortify methods to support additional fields <2015-04-28, Tue>

- extend get.fields methods to support additional fields <2015-04-28,
  Tue>

- extend tree class to support additional info by merging two tree
  <2015-04-28, Tue>

- implement merge_tree function to merge two tree objects into one
  <2015-04-28, Tue>

Changes in version 1.1.1:

- minor bug fixed in extracting node ID of rst file <2015-04-27, Mon>

- update parsing beast time scale tree to support _year (originally
  supports _year.\\d+) <2015-04-27, Mon>

- add Tommy in author <2015-04-27, Mon>

girafe
------

Changes in version 1.21.1:

- replaced seq_name method by seqnames (due to deprecation of seq_name
  generic)

GOexpress
---------

Changes in version 1.3.2:

BUG FIX

- expression*plot() functions return error if an empty gene name is
  given

Changes in version 1.3.1:

BUG FIX

- heatmap_GO() function does not automatically resize margins when
  using either gene names or gene identifiers. This caused the
  user-defined margins to be systematically ignored.


GOSemSim
--------

Changes in version 1.27.4:

- convert vignette from Rnw to Rmd <2015-06-23, Tue>

Changes in version 1.27.3:

- bug fixed in getSupported_Org <2015-05-31, Sun>

Changes in version 1.27.2:

- deprecate 'worm' and use 'celegans' instead.

Changes in version 1.27.1:

- add external documents <2015-05-07, Thu>

graphite
--------

Changes in version 1.15.3 (2015-10-07):

- Updated all pathway data.

- More descriptive edge types.

- Use the package "rappdirs" to select a directory for cached files.

GreyListChIP
------------

Changes in version 1.1.1:

- Change maintainer email.

GSAR
----

Changes in version 1.4.0:

- New argument return.weights added to function plotMST2.pathway. If
  return.weights=TRUE, the weight vectors found by GSNCA for the genes
  under two classes are returned as a 2-column matrix.


GSEABase
--------

Changes in version 1.31:

SIGNIFICANT USER-VISIBLE CHANGES

- Support 'h' Broad set

BUG FIXES

- as(OBOCollection, "graphNEL") failed to include multiple children

gtrellis
--------

Changes in version 1.1.1:

- print message if input chromosome names are without 'chr'

- call `start()` and `end()` explictely with `GenomicRanges::`

- unit of x-axis is correct now

Guitar
------

Changes in version 0.99.9:

- package renamed as "guitar"

gwascat
-------

Changes in version 2.1.0:

USER VISIBLE CHANGES

- GWAS catalog data curated at EMBL/EBI will now be used with
  makeCurrentGwascat: Updated 3 August 2015, as comma-space delim. for
  genes and EFO tags has been introduced

- Interfaces to the Experimental Factor Ontology (www.ebi.ac.uk/efo/)
  are provided: efo.obo.g is an annotated graphNEL with the ontology

- Vignette gwascatOnt.Rmd introduced

GWASTools
---------

Changes in version 1.15.16:

- added permute argument to exactHWE

Changes in version 1.15.15:

- allow multiple color schemes for plots color-coded by genotype

Changes in version 1.15.13:

- pedigreePairwiseRelatedness identifies great grandparent/great
  grandchild (GGp) and grand avuncular (GAv)

Changes in version 1.15.12:

- allow character scanID in createDataFile

Changes in version 1.15.11:

- added col argument to manhattanPlot

Changes in version 1.15.10:

- createDataFile converts non-finite values to NA.

Changes in version 1.15.8:

- alleleFrequency includes scans with missing sex.

Changes in version 1.15.7:

- Added option to reorder samples in vcfWrite.

Changes in version 1.15.5:

- Added option to read genotypes coded with nucleotides in
  createDataFile.

Changes in version 1.15.3:

- Added beta and standard error for GxE term to assocRegression output.

Changes in version 1.15.2:

- Added number of cases and controls to assocRegression output.

Changes in version 1.15.1:

- Added "ci" argument to qqPlot function.

HDTD
----

Changes in version 1.3.4 (2015-09-08):

- Updated email address of maintainer.

Changes in version 1.3.3 (2015-07-20):

- Updated output format of core functions.

- Fixed minor bugs.

Changes in version 1.3.2 (2015-05-26):

- Updated DESCRIPTION FILE.

- Updated the CITATION file.

Changes in version 1.3.1 (2015-05-12):

- Updated DESCRIPTION FILE.

- Updated the CITATION file.

- Updated documentation.

Heatplus
--------

Changes in version 2.15.2:

- New color scheme for clusters to make neighboring clusters distinct

hiAnnotator
-----------

Changes in version 1.3.2:

- fixed duplicated featureName column in getSitesInFeature when
  allSubjectCols is TRUE.

Changes in version 1.3.1:

- replaced plyr with dplyr

HIBAG
-----

Changes in version 1.6.0:

- the version number was bumped for the Bioconductor release version
  3.2

hierGWAS
--------

Changes in version 0.99.3:

- Removed hlr dependency and added the MEL function

- Fixed computation of r2 when the SNP_index = NULL

HilbertCurve
------------

Changes in version 0.99.3:

- transparency is also averaged

Changes in version 0.99.2:

- jump to top viewport for each low-level functions

Changes in version 0.99.1:

- add dependency to base packages

Changes in version 0.99.0:

- new release

hiReadsProcessor
----------------

Changes in version 1.4.0:

- Using dplyr for small computations

HiTC
----

Changes in version 1.13.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Change in the way the interaction matrix is loaded. During the
  import, the lazyload option will force the data to be stored as
  triangular matrix. Note that a symmetric matrix is not triangular by
  construction. Then, to avoid any error in the data processing, the
  triangular matrix is always converted into symmetrical matrix before
  being returned by the intdata method.

- Update of track display in mapC function. Off-set plot of adjacent
  features

BUG FIXES

- Bug fixed in quality control function when empty maps are used

- Bug fixed when plotting empty matrix

Changes in version 1.13.1:

NEW FEATURES

- Export methods getCombinedContacts, getCombinedIntervals for HTClist
  objects

SIGNIFICANT USER-VISIBLE CHANGES

- getCombinedContact is now able to merge HTCexp objects for non
  complete HiTClist objects. Missing maps are replaced by NA matrices

- Update of isComplete, isPairwise, getCombinedIntervals methods for
  Hi-C data with no intrachromosomal maps

- When the maxrange argument is set in the mapC function, all maps are
  displayed on the same scale so that they can be compared to each
  other

BUG FIXES

- Bug fixed in mapC with the contact map is empty

- Bug fixed in seqlevels(HTClist) method

hpar
----

Changes in version 1.11.1:

- unit tests <2015-06-30 Tue>

Changes in version 1.11.0:

- Bio version 3.2 (devel)

iGC
---

Changes in version 1.0.0:

- Initial release

illuminaio
----------

Changes in version 0.11.2 (2015-09-11):

- BUG FIX: readIDAT() can now read non-encrypted IDAT files with
  strings longer than 127 characters.

- BUG FIX: readIDAT() incorrectly assumed that there were exactly two
  blocks in RunInfo fields of non-encrypted (v3) IDAT files. Thanks to
  Gordon Bean (GitHub @brazilbean) for reporting on and contributing
  with code for the above two bugs.

Changes in version 0.11.1 (2015-07-29):

- Updated the BiocViews field of DESCRIPTION.

Changes in version 0.11.0 (2015-04-16):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.2 for R (>= 3.3.0).

immunoClust
-----------

1.1.1: NEW FEATRUES: * A normalization step which the sample
        cell-clusters to the common meta-cluster model is included an
        optionally activated during the major meta-clustering process.
        CHANGES: * The meta.ME C-binding and return value was modified
        in a way that the A-Posterior probability matrix Z for a
        cell-cluster belonging to a meta-cluster is also calculated and
        returned. BUFIXES: * Ellipse position were not correct when
        ploting a parameter subset

InPAS
-----

Changes in version 1.1.9:

- update documentation

Changes in version 1.1.8:

BUG FIXES

- improve the distCP site search.

Changes in version 1.1.7:

BUG FIXES

- typo "novel distal".

Changes in version 1.1.6:

BUG FIXES

- fix the seqlevels style issue.

Changes in version 1.1.5:

BUG FIXES

- improve limmaAnalyze algorithm

Changes in version 1.1.4:

BUG FIXES

- improve the CP sites adjustment algorithm.

Changes in version 1.1.3:

BUG FIXES

- shorten the CPU time of CP sites prediction

Changes in version 1.1.2:

BUG FIXES

- update description

Changes in version 1.1.1:

BUG FIXES

- fix the bug in FFT smooth when length of 3UTR is shorter than smooth
  power

interactiveDisplayBase
----------------------

Changes in version 1.7:

NEW FEATURES

- .runApp runs the app in RStudio's 'viewer' pane (if the app is
  launched under RStudio), or in the browser.

BUG FIXES

- Applications would only start under some versions of shiny, due to a
  reference to either 'rstudio' or 'rstudioapp' search path element.


IRanges
-------

Changes in version 2.4.0:

NEW FEATURES

- Add "cbind" methods for binding Rle or RleList objects together.

- Add coercion from Ranges to RangesList.

- Add "paste" method for CompressedAtomicList objects.

- Add "expand" method for Vector objects for expanding a Vector object
  'x' based on a column in mcols(x).

- Add overlapsAny,integer,Ranges method.

- coverage" methods now accept 'shift' and 'weight' supplied as an Rle.

SIGNIFICANT USER-VISIBLE CHANGES

- The following was moved to S4Vectors: - The FilterRules stuff.  - The
  "aggregate" methods.  - The "split" methods.

- The "sum", "min", "max", "mean", "any", and "all" methods on
  CompressedAtomicList objects are 100X faster on lists with 500k
  elements, 80X faster for 50k elements.

- Tweak "c" method for CompressedList objects to make sure it always
  returns an object of the same class as its 1st argument.

- NCList() constructor now propagates the metadata columns.

DEPRECATED AND DEFUNCT

- RangedData/RangedDataList are not formally deprecated yet but the
  documentation now officially declares them as superseded by
  GRanges/GRangesList and discourage their use.

- After being deprecated in BioC 3.1, IntervalTree and IntervalForest
  objects and the "intervaltree" algorithm in findOverlaps() are now
  defunct.

- After being deprecated in BioC 3.1, mapCoords() and pmapCoords() are
  now defunct.

- Remove seqapply(), mseqapply(), tseqapply(), seqsplit(), and seqby()
  (were defunct in BioC 3.1).

BUG FIXES

- Fix FactorList() constructor when 'compress=TRUE' (note that the
  levels are combined during compression).

- Fix c() on CompressedFactorList objects (was returning a
  CompressedIntegerList object).

kebabs
------

Changes in version 1.3.4:

- correction of Ubuntu problem with realloc for 0 elements in
  linearKernel generating a sparse empty kernel matrix

- correction of problem with feature weights and prediction profiles
  for position specific gappy pair kernel

- correction of problem with feature weights and prediction profiles
  for position specific motif kernel

- corrections for feature weights, prediction via feature weights and
  prediction profile for distance weighted kernels

- update of KeBABS citation

Changes in version 1.3.3:

- new export kebabsCollectInfo for collection of package info

- update of version dependency to Biostrings, XVector, S4Vector

- correction for leading + or - in factor label

- change of bibtex style sheet in vignette to plainnat.bst

Changes in version 1.3.2:

- correction of error in kernel lists

- user defined sequence kernel example SpectrumKernlabKernel moved to
  separate directory

Changes in version 1.3.1:

- correction of error in model selection for processing via dense
  LIBSVM

- remove problem in check for loading of SparseM

Changes in version 1.3.0:

- first devel version created from release version 1.2.0

KEGGgraph
---------

Changes in version 1.27.2 (2015-05-25):

- Only a subset of functions are imported from graph. It should solve
  the problem reported here:
  http://stackoverflow.com/questions/30428860/proper-use-of-optional-package-features-and-dependencies.

Changes in version 1.27.1 (2015-05-04):

- Thanks to input of Anders Ellern Bilgrau, a bug is fixed which caused
  multiple edges of the graphNEL object



mAPKL
-----

Changes in version 1.1.1:

- A new function, the "probes2pathways" has been added.

- I have added to more network reconstruction approaches. The
  "aracne.a" and the "aracne.m" as included in the "parmigene"
  r-package.

- I have replaced in the "preprocess" function the type of the produced
  image from jpeg to tiff format, with resolution=300dpi.


metabomxtr
----------

Changes in version 1.3.6:

- *made small modifications to arguments for function mixnorm *added
  vignette for mixnorm

Changes in version 1.3.1:

- *added mixnorm function for normalization of non-targeted
  metabolomics studies

metagene
--------

Changes in version 2.2.0:

- Added checks to avoid producing identical matrices or data frame when
  the parameters are still the same after first function call.

- Splitted the analysis in multiple (optionnal) intermediate steps
  (add_design, produce_matrices and produce_data_frame).

- narrowPeak and broadPeak format is now supported.

- Added multiple getter to access metagene members that are all now
  private (get_params, get_design, get_regions, get_matrices,
  get_data_frame, get_plot, get_raw_coverages and
  get_normalized_coverages.

- Added the NCIS algorithm for noise removal.

- Replaced the old datasets with promoters_hg19, promoters_hg18,
  promoters_mm10 and promoters_mm9 that can be accessed with
  data(promoters_????).

- Added flip_regions and unflip_regions to switch regions orientation
  based on the strand.

metagenomeFeatures
------------------

Changes in version 0.0.0.9 (2015-09-14):

- Pre-Release Bioconductor

metagenomeSeq
-------------

Changes in version 1.11:

- Adding fitFeatureModel - a feature based zero-inflated log-normal
  model.

- Added MRcoefs,MRtable,MRfulltable support for fitFeatureModel output.

- Added mention in vignette.

- Added support for normalizing matrices instead of just MRexperiment
  objects.

- Fixed cumNormStat's non-default qFlag option

metaseqR
--------

Changes in version 1.9.21 (2015-10-08):

NEW FEATURES

- Added one more output list member of the main pipeline (metaseqr)
  containing all the processed data in a way that can be easily used
  for further downstream analysis in R. Contains also the metaseqr
  pipeline call and the parameters.

BUG FIXES

- Fixed broken dependency with removed CRAN package MADAM which was
  used for the Fisher p-value combination. The fix was done by copying
  the two required functions from the last MADAM archived version
  (1.2).

- Fixed a bug in make.sample.list which rendered the function unusable
  (credits to Marina Adamou-Tzani).

- Fixed a very special-case small bug in the report generation only a
  single gene passes the multiple testing correction cutoff (thans to
  Martic Reszcko, BSRC 'Alexander Fleming').

- Fixed problem with Fisher p-value combination method (returning all
  NA due to the inexistence of names in the p-value vecror).

- Fixed bug with flags of biotype filtered genes.

Changes in version 1.9.0 (2015-05-27):

NEW FEATURES

- Added ability to export JSON for a few graphs to be used with
  Highcharts library for interactive visualization.

BUG FIXES

- None.

metaX
-----

Changes in version 0.1.0:

- First release version

methylPipe
----------

Changes in version 1.3.3:

- updated or fixed the following functions: + profileDNAmetBinParallel:
  fixed the error when there is single minus strand in GRanges for
  profiling

Changes in version 1.3.2:

- updated or fixed the following functions: + meth.call: coverage is
  outputted now as an GRanges object + findDMR: coverage threshold is
  used even for unmethylated cytosines + getCposChr: fixed bug which
  was ignoring Cytosines at the boundary of genomic ranges on both
  sides

minfi
-----

Changes in version 1.15:

- Adding testing for preprocessNoob, preprocessFunnorm.

- Fxing some verbose output of preprocessNoob.

- Adding non-exported function .digestVector for testing.

MLInterfaces
------------

Changes in version 1.49.8:

- Enable learnerIn3D customisation <2015-07-28 Tue>

Changes in version 1.49.7:

- fix proper asNA0 propagation in performance analytics functions
  (implemented 2014-10-16 Thu, commited on 2015-07-28 Tue)

Changes in version 1.49.1:

- add hclustWidget

mogsa
-----

Changes in version 1.1.1:

NEW FEATURES

- updated the moCluster functions, including the following exported
  function: "bootMbpca", "distMoa", "moGap", "softK", "mbpca",
  "moaScore", "moaCoef"

motifStack
----------

Changes in version 1.13.8:

NEW FEATURES

- improve plot plotMotifStackWithPhylog.

Changes in version 1.13.7:

- update the documentation

Changes in version 1.13.6:

BUG FIXES

- fix bug when plot 2 motifs with motifStack function.

Changes in version 1.13.5:

NEW FEATURES

- improve the algorithm of DNAmotifAlignment. Replaced the ALLR by
  average information content.

Changes in version 1.13.4:

NEW FEATURES

- improve the algorithm of DNAmotifAlignment

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.13.3:

NEW FEATURES

- add reorderUPGMAtree function

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.13.2:

NEW FEATURES

- add mergeMotifs function

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.13.1:

NEW FEATURES

- add plotMotifOverMotif function

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

msa
---

Changes in version 1.2.0:

- new branch for Bioconductor 3.2 release


MSnbase
-------

Changes in version 1.17.16:

- new default parameter 'feature weight' in iPQF by Martina Fisher (see
  PR#65) <2015-10-07 Wed>

Changes in version 1.17.15:

- Fix typo in plot man <2015-08-23 Sun>

Changes in version 1.17.14:

- Fix typo in impute man <2015-08-17 Mon>

Changes in version 1.17.13:

- partly rewrite writeMgfData <2015-05-16 Thu>

- initial hmap function <2015-07-16 Thu>

- fix bug in plotting MS1 spectra (closes issue #59) <2015-07-16 Thu>

- new image implementation, based on @vladpetyuk's
  vp.misc::image_msnset <2015-07-25 Sat>

- Changed the deprecated warning to a message when reading MzTab data
  version 0.9, as using the old reader can not only be achived by
  accident and will be kept for backwards file format compatibility
  <2015-07-30 Thu>

Changes in version 1.17.12:

- fix show MIAPE when some otherInfo at NA <2015-07-15 Wed>

Changes in version 1.17.11:

- adding unit tests <2015-07-01 Wed>

- fix abundance column selection when creating MSnSet form MzTab
  <2015-07-01 Wed>

- new mzTabMode and mzTabType shortcut accessors for mode and type of
  an mzTab data <2015-07-01 Wed>

Changes in version 1.17.10:

- Fix URL <2015-06-30 Tue>

Changes in version 1.17.9:

- calculateFragments' "neutralLoss" argument is now a list (was a
  logical before), see #47. <2015-06-24 Wed>

- add defaultNeutralLoss() function to fill calculateFragments'
  modified "neutralLoss" argument, see #47. <2015-06-24 Wed>

Changes in version 1.17.8:

- coercion from IBSpectra to MSnSet, as per user request <2015-06-23
  Tue>

- new iPQF combineFeature method <2015-06-24 Wed>

Changes in version 1.17.7:

- Fix support of identification-only mzTab files <2015-06-22 Mon>

Changes in version 1.17.6:

- Export metadata,MzTab-method <2015-06-19 Fri>

- Replace spectra,MzTab-method by psms,MzTab-method <2015-06-20 Sat>

- Change the meaning of calculateFragments' "modifications" argument.
  Now the modification is added to the mass of the amino acid/peptide.
  Before it was replaced. <2015-06-21 Sun>

- calculateFragments gains the feature to handle N-/C-terminal
  modifications, see #47. <2015-06-21 Sun>

- update readMzTabData example <2015-06-22 Mon>

Changes in version 1.17.5:

- new MzTab class to store a simple parsing of an mzTab version 1 data
  file. See ?MzTab for details. <2015-06-16 Tue>

Changes in version 1.17.4:

- New lengths method for FoICollection instances <2015-06-06 Sat>

- New image2 function for matrix object, that behaves like the method
  for MSnSets <2015-06-10 Wed>

- image,MSnSet labels x and y axis as Samples and Features <2015-06-10
  Wed>

- fixed bug in purityCorrect, reported by Dario Strbenac <2015-06-11
  Thu>

Changes in version 1.17.3:

- iTRAQ 8-plex correction factors and impurity matrix <2015-05-22 Fri>

Changes in version 1.17.2:

- new filterZero function <2015-05-01 Fri>

Changes in version 1.17.1:

- new MSnSetList class <2015-04-19 Sun>

- new commonFeatureNames function <2015-04-14 Tue>

- new compareMSnSets function <2015-04-19 Sun>

- splitting and unsplitting MSnSets/MSnSetLists <2015-04-19 Sun>

Changes in version 1.17.0:

- new devel version for Bioc 3.2

MSnID
-----

Changes in version 1.3.1:

- S4 methods "accessions" and "proteins" outsourced to ProtGenerics

MSstats
-------

Changes in version 3.0.12:

- remove vignetter folder to remove install and build error in
  Bioconductor

Changes in version 3.0.9:

- dataProcess

- add options for ‘cutoffCensored=“minFeatureNRun”’.

- summaryMethods=“TMP” : output will have ‘more50missing’column.

- remove50missing=FALSE option : remove runs which has more than 50% of
  missing measurement. It will be affected for TMP, with censored
  option.

- MBimpute : impute censored by survival model (AFT) with cutoff
  censored value

- featureSubset option : “all”,”top3”, “highQuality”

- change the default. -groupComparisonPlots

- heatmap, for logBase=10, fix the bug for setting breaks.

Changes in version 3.0.8:

- dataProcess : when censoredInt=“0”, intensity=0 works even though
  skylineReport=FALSE.

- dataProcess, with censored=“0” or “NA” : fix the bug for certain run
  has completely missing.

- cutoffCensored=“minRun” or “minFeature” : cutoff for each Run or each
  feature is little less (99%) than minimum abundance.
  -summaryMethod=“TMP”, censored works. censoredInt=NA or 0, and
  cutoffCensored=0, minFeature, minRun

Changes in version 3.0.3:

- dataProcess : new option, skylineReport. for skyline MSstats report,
  there is ‘Truncated’ column. If Truncated==True, remove them. and
  keep zero value for summaryMethod=“skyline”.

- groupComparison : for skyline option, t.test, val.equal=TRUE, which
  is no adjustment for degree of freedom, just pooled variance.

mzID
----

Changes in version 1.7.1:

- Imports generics from ProtGenerics

- Export base classes

mzR
---

Changes in version 2.3.3:

- bump to new Rcpp 0.12.1 version <2015-09-29 Tue>

Changes in version 2.3.2:

- Fix typo/bug in peaks,mzRpwiz <2015-07-31 Fri>

Changes in version 2.3.1:

- update Rcpp version message to point to support site <2015-05-08 Fri>

NanoStringDiff
--------------

Changes in version 0.99.0:

NEW FEATURES

- Initial release.

nethet
------

1.1.1: package. Added covariance matrix as output of screen_cvglasso.

nondetects
----------

Changes in version 1.99.0:

- Added option for the output of model parameters.

- Added option for multiple imputation.

- Added option for the model function.

npGSEA
------

Changes in version 1.5.1:

- Added paper citation

omicade4
--------

Changes in version 1.9.1:

NEW FEATURES

- new function topVar

OmicsMarkeR
-----------

Changes in version 1.1.0:

- Changes:
  
  * the modelList has been added to provide the user a list of
  currently implemented methods.
  
  * The 'verbose' argument in fs.stability, fs.ensembl.stability, &
  fit.only.model has been changed to a character option indicating the
  extent of verbose output.
  
  * Canberra stability has been added as the previous implementation
  was not compatible with RPT. See ?canberra.stability for more
  details.
  
  * Many more tests have been implemented to make the package more
  stable.

OncoSimulR
----------

Changes in version 1.99.9 (2015-10-08):

- Fixed NEWS file.

- Removed empty file.

Changes in version 1.99.8 (2015-10-01):

- Test "initMutant with oncoSimulSample, 2" occasionally failed.

Changes in version 1.99.7 (2015-09-27):

- initMutant available in oncoSimulPop and oncoSimulSample.

Changes in version 1.99.6 (2015-09-26):

- Improved test coverage and removed stringsAsfactors from tests.

- Consistent handling of corner cases in Bozic.

- Miscell minor documentation improvements.

Changes in version 1.99.5 (2015-06-25):

- Fixed bug in initMutant, added tests, and vignette section.

Changes in version 1.99.4 (2015-06-22):

- Plotting true phylogenies.

- Tried randutils, from O'Neill. Won't work with gcc-4.6.

- bool issue in Windows/gcc-4.6.

- Most all to all.equal in tests.

Changes in version 1.99.3 (2015-06-19):

- More examples to vignette

- Using Makevars

- More functionality to plot.fitnessEffects

- Will Windoze work now?

Changes in version 1.99.2 (2015-06-19):

- Fixed typos and other minor in vignett.

Changes in version 1.99.01 (2015-06-17):

- Many MAJOR changes: we are done moving to v.2

- New way of specifying restrictions (v.2) that allows arbitrary
  epistatic interactions and order effects, and very large (larger than
  50000 genes) genomes.

- When onlyCancer = TRUE, all iterations now in C++.

- Many tests added.

- Random DAG generation.

- Some defaults for v.1 changed.

Changes in version 1.99.1 (2015-06-18):

- Try to compile in Windoze with the SSTR again.

- Reduce size of RData objects with resaveRdaFiles.

- Try to compile in Mac: mt RNG must include random in all files.

Changes in version 1.99.00 (2015-04-23):

- Accumulated changes of former 99.1.2 to 99.1.14:

- changes in intermediate version 1.99.1.14 (2015-04-23):

- Now are things OK (I messed up the repos)

- changes in intermediate version 99.1.13 (2015-04-23)

- Added a couple of drop = FALSE. Their absence lead to crashes in some
  strange, borderline cases.

- Increased version to make unambiguous version used for anal. CBN.

- changes in intermediate version 99.1.12 (2015-04-22)

- Removed lots of unused conversion helpers and added more strict
  checks and tests of those checks.

- changes in intermediate version 99.1.11 (2015-04-18)

- Tests of conversion helpers now really working.

- changes in intermediate version 99.1.10 (2015-04-16)

- Added conversion helpers as separate file.

- Added tests of conversion helpers.

- Added generate-random-trees code (separate file).

- More strict now on the poset format and conversions.

- changes in intermediate version 99.1.9 (2015-04-09)

- added null mutation for when we run out of mutable positions, and
  since not clear how to use BNB then.

- changes in intermediate version 99.1.8 (2015-04-03)

- added extraTime.

- changes in intermediate version 99.1.7 (2015-04-03)

- endTimeEvery removed. Now using minDDrPopSize if needed.

- changes in intermediate version 99.1.6 (2015-03-20)

- untilcancer and oncoSimulSample working together

- changes in intermediate version 99.1.5 (2015-03-20)

- Using the untilcancer branch

- changes in intermediate version 99.1.4 (2014-12-24)

- Added computation of min. of ratio birth/mutation and death/mutation.

- changes in intermediate version 99.1.3 (2014-12-23)

- Fixed segfault when hitting wall time and sampling only once.

- changes in intermediate version 99.1.2 (2014-12-16)

- Sampling only once

openCyto
--------

Changes in version 1.7.1:

Enhancements

- New API: add_pop function allows users to apply a single gating
  method to GatingSet without writing the compelete csv template


PAA
---

Changes in version 1.3.3 (2015-07-14):

GENERAL

- No changes.

NEW FEATURES

- No changes.

IMPROVEMENTS

- No changes.

MODIFICATIONS

- No changes.

BUG FIXES

- Bug in the function loadGPR() fixed (due to a wrong regular
  expression some data rows were not imported from gpe files in PAA
  versions 1.3.2 and 1.3.1 -> the bug was not relevant for the versions
  1.3.0 and older).

Changes in version 1.3.2 (2015-06-22):

GENERAL

- Update to the latest R version which has been released a few days ago
  (2015-06-18).

- Built with R-3.2.1

NEW FEATURES

- No changes.

IMPROVEMENTS

- No changes.

MODIFICATIONS

- NEWS file of version 1.3.1 corrected (wrong version number for the
  last update in the NEWS file).

BUG FIXES

- No changes.

Changes in version 1.3.1 (2015-06-17):

GENERAL

- The optional rlm normalization (for ProtoArrays) used by the
  functions normalizeArrays(), plotNormMethods() and plotMAPlots() has
  been completely reimplemented in order to fix some bugs and simplify
  the usage of these functions (details: see below).

- Built with R-3.2.0

NEW FEATURES

- No changes.

IMPROVEMENTS

- Due to the modification of the EListRaw objects for ProtoArrays the
  usage of the functions normalizeArrays(), plotNormMethods() and
  plotMAPlots() is simplified since some obsolete arguments have been
  removed (see 'MODIFICATIONS'). Esp., the provision of a second
  controls-specific EListRaw object for controls data and information
  was too complex and has confused some users. Now all control
  spot-specific data and information are stored together with
  probe-specific data and information in one extended EListRaw object.

MODIFICATIONS

- On the one hand, the arguments 'controls.elist', 'gpr.path',
  'targets.path' and 'contr.names' have been removed from the functions
  normalizeArrays(), plotNormMethods() and plotMAPlots(). On the other
  hand, the argument 'controls' has been added to the functions
  normalizeArrays(), plotNormMethods() and plotMAPlots().

BUG FIXES

- Due to the re-implementation of the optional rlm normalization (for
  ProtoArrays) used by the functions normalizeArrays(),
  plotNormMethods(), plotMAPlots(), all bugs that have been reported by
  users since version 1.0.0 are fixed.

pandaR
------

Changes in version 3.2.0:

- Initial version of pandaR implements the PANDA algorithm for gene
  regulatory network inference from gene expression, sequence motif and
  protein-protein interaction data.

- Includes functions for subsetting bipartite network and discretizing
  TF-gene edgeweights.

- Utilizes igraph library for displaying a subsetted gene regulatory
  network of interest.


pathview
--------

Changes in version 1.9.3:

- pathview can accept a vector of multiple pathway ids, and map/render
  the user data onto all these pathways in one call.

- one extra column "all.mapped" was added to pathview output
  data.frames as to show all the gene/compound IDs mapped to each node.

- add geneannot.map as a generic function for gene ID or annotation
  mapping.

- sim.mol.data now generate data with all major gene ID types for all
  19 species in bods, not just human.

- download.kegg now let the user to choose from xml, png or both file
  types to download for each input pathway. In the meantime, it uses
  the KEGG REST API instead of the classical KEGG download links. All
  potential pathways including the general pathways can be downloaded
  this way.

- solve the redundant import from graph package.

- import specific instead of all functions from XML package.

Changes in version 1.9.1:

- solve the install error due to the recent change in KEGGgraph
  package.

Pbase
-----

Changes in version 0.9.1:

- update cleave method to address build error <2015-10-13 Tue>

Changes in version 0.9.0:

- Bioc devel 3.2

PECA
----

Changes in version 1.5.2 (2015-08-07):

ADDED FUNCTIONS

- PECASI - Differential splicing between two groups of samples in
  Affymetrix exon array studies.


PGA
---

Changes in version 0.1.0:

- First release version

PhenStat
--------

Changes in version 2.3.1:

NEW FUNCTIONALITY

- New output format is added: JSON format. In odred to get output in
  JSON format apply JSONOutput with the PhenTestResult object as an
  argument.

COMPATIBILITY ISSUES

- The Box-Cox transformation is switched off by default. In order to
  swithc in on the testDataset() function's argument “transformValues"
  has to be set to TRUE.

phyloseq
--------

Changes in version 1.13.6:

BUG FIXES

- droplevels suggestion for sample-data
  https://github.com/joey711/phyloseq/pull/476

- DESeq2 migrated to suggests
  https://github.com/joey711/phyloseq/pull/533

- `extend_metagenomeSeq` functionality
  https://github.com/joey711/phyloseq/pull/533

- bugs related to previous version distance uptick, mostly in tests and
  vignette

Changes in version 1.13.5:

USER-VISIBLE CHANGES

- Help avoid cryptic errors due to name collision of `distance` with
  external loaded packages by making `distance` a formal S4 method in
  phyloseq.

- Improve documentation of `distance` function and the downstream
  procedures on which it depends

- Migrate the list of supported methods to a documented, exported list
  object, called `distanceMethodList`.

- Improved distance unit tests with detailed checks that dispatch works
  and gives exactly expected distance matrices for all methods defined
  in distanceMethodList.

- Improved JSD doc, performance, code, deprecated unnecessary
  `parallel` argument in JSD

Changes in version 1.13.4:

BUG FIXES

- `psmelt` bug if user has also loaded the original "reshape" package,
  due to name collision on the function called `melt`. `psmelt` now
  explicitly calls `reshape2::melt` to avoid confusion.
  https://github.com/joey711/phyloseq/pull/489

- Fix following note... There are ::: calls to the package's namespace
  in its code. A package almost never needs to use ::: for its own
  objects: ‘JSD.pair’

piano
-----

Changes in version 1.10:

- No changes yet


podkat
------

Changes in version 1.1.3:

- further fix of weights() method for signature 'AssocTestResultRanges'

Changes in version 1.1.2:

- fix of weights() method for signature 'AssocTestResultRanges'

Changes in version 1.1.1:

- fix of filterResults() method for signature 'GRanges'

Changes in version 1.1.0:

- new branch for Bioconductor 3.2 devel

polyester
---------

1.99.3: NB function now exported

1.99.3: note that version 1.99.3 on GitHub was version 1.1.0 on
        Bioconductor.

1.99.2: bug fix in fragment generation (last 2 bases of transcript were
        never sequenced)

pRoloc
------

Changes in version 1.9.7:

- New SpatProtVis visualisation class <2015-08-13 Thu>

- add link to explanation of supportive/uncertain reliability scores in
  tl vignette <2015-09-02 Wed>

Changes in version 1.9.6:

- Update REAMDE with TL ref

- Update refs in lopims documentation <2015-07-30 Thu>

Changes in version 1.9.5:

- update doc <2015-07-15 Wed>

Changes in version 1.9.4:

- Add reference to TL paper and link to lpSVM code <2015-07-06 Mon>

- highlightOnPlot throws a warning and invisibly returns NULL instead
  of an error when no features are in the object <2015-07-08 Wed>

- highlightOnPlot has a new labels argument <2015-07-10 Fri>

Changes in version 1.9.3:

- Clarify error when no annotation params are provided <2015-05-11 Mon>

- support for matrix-encoded markers <2015-05-19 Tue>

- New default in addLegend: bty = "n" <2015-05-20 Wed>

- getMarkers now supports matrix markers <2015-05-20 Wed>

- getMarkerClasses now supports matrix markers <2015-05-20 Wed>

- markerMSnSet and unknownMSnSet now support matrix markers <2015-05-20
  Wed>

- sampleMSnSet now supports matrix markers <2015-05-23 Sat>

- updated yeast markers and added uniprot ids <2015-05-27 Wed>

- plot2D support a pre-calculated dim-reduced data matrix as method
  parameter to avoid recalculation <2015-05-27 Wed>

Changes in version 1.9.2:

- Lisa's colour palette <2015-05-08 Fri>

Changes in version 1.9.1:

- new plot2Ds function to overlay two data sets on the same PCA plot
  [2015-04-17 Fri]

- regenerate biomart data used by setAnnotationParams [2015-04-24 Fri]

- new setStockcolGui function to set the default colours manually via a
  simple interface [2015-04-29 Wed]

- new move2Ds function to produce an transition movie between two
  MSnSets [2015-04-29 Wed]

- functions to convert GO ids to/from terms. See ?goTermToId for
  details <2015-05-08 Fri>

Changes in version 1.9.0:

- new devel version for Bioc 3.2

pRolocGUI
---------

Changes in version 1.3.2:

- Fix check warnings and errors <2015-07-20 Mon>

Changes in version 1.3.1:

- new plotMat2D function <2015-05-20 Wed>

- Fix query search in pRolocVis, contributed by pierremj <2015-05-27
  Wed>

- pRolocVis has new method arg <2015-05-29 Fri>

Changes in version 1.3.0:

- Devel for Bioc 3.2

Prostar
-------

Changes in version 0.99.0:

- First submission of Prostar to Bioconductor

ProteomicsAnnotationHubData
---------------------------

Changes in version 0.99.3:

- removing .get1 methods and pointing to relevant AnnotationHub file
  <2015-10-09 Fri>

Changes in version 0.99.2:

- Use AnnotationHub::query() instead of manual subsetting <2015-09-02
  Wed>

Changes in version 0.99.1:

- Submission of Bioc (Issue 1292) <2015-08-25 Tue>

Changes in version 0.99.0:

- First release with data PXD000001

ProtGenerics
------------

Changes in version 1.1.1:

- new mz<- generic function <2015-10-01 Thu>

Changes in version 1.1.0:

- Bioc devel version 3.2

PWMEnrich
---------

Changes in version 4.5.1:

- Convert log(P-values) back to P-values for human using a chi-sq
  distribution Version 4:

- New algorithm for human backgrounds

- New function: toPWM() that takes both PFMs and PPMs

pwOmics
-------

Changes in version 1.1.10:

- include threshold parameter into 'readTFdata' function

- change STRINGdb to version 10

Changes in version 1.1.1:

- rename consensus-based dynamic analysis

- adjust vignette workflow figure

QDNAseq
-------

Changes in version 1.6.0:

RELEASE

- Bioconductor 3.2

IMPROVEMENTS

- chromosomes() no longer tries to return chromosome names as an
  integer vector, but returns a character vector instead

BUG FIXES

- plot() and calculateMappability() now work also for other chromosome
  names besides numbers and "X", "Y", and "MT"

Changes in version 1.4.2 (2015-08-20):

BUG FIXES

- createBins() now properly selects chromosomes when ignoring
  mitochondrial DNA. Please note that the mitochondrial DNA is only
  ignored when it is called either "chrMT", "chrM", "MT", or "M"

Changes in version 1.4.1 (2015-06-30):

IMPROVEMENTS

- correctBins() now filters out bins with missing loess correction
  estimate

qpgraph
-------

Changes in version 2.40:

BUG FIXES

- Bugfix on qpPrecisionRecall() when argument 'refGraph' is a graphBAM
  object.

QuasR
-----

Changes in version 1.10.0:

NEW FEATURES

- qExportWig gained createBigWig argument and can now create bigWig
  files directly

- qQCReport now also produces base quality plots for bam-file projects
  by sampling reads from the bam files

- qCount, qProfile and qExportWig have gained a includeSecondary
  argument to include/exclude secondary alignments while counting

qvalue
------

Changes in version 2.1.1:

- handles NA values

- upgraded plotting functions

R3CPET
------

1.2.0: Updates: * Removed dependency from DAVIDQuery package as it will
        deprected.  * Fixed some bugs in DAVIDQuery functions and
        integrated them into R3CPET.  * Updated the Readme.Rd file.  *
        Updated the HPRD.RData and Biogrid.RData to the new igraph
        class.  * some small changes.

RCyjs
-----

Changes in version 1.9.8:

BUG FIXES

- Significant (3x) speedup.  A 5000-node, 6000-edge graph transmits to
  Cytoscape from R in about 20 seconds.


ReactomePA
----------

Changes in version 1.13.5:

- add support of fly (Drosophila melanogaster) <2015-09-22, Tue>

Changes in version 1.13.4:

- update vignette <2015-08-19, Wed>

Changes in version 1.13.3:

- add citation of ChIPseeker <2015-07-09, Thu>

- add 'Pathway analysis of NGS data' section in vignette <2015-06-29,
  Mon>

- convert vignette from Rnw to Rmd <2015-06-29, Mon>

Changes in version 1.13.2:

- enable other organisms for viewPathway by Vladislav Petyuk
  <2015-05-28, Thu> see
  https://github.com/GuangchuangYu/ReactomePA/pull/1

Changes in version 1.13.1:

- update enrichMap and viewPathway <2015-05-13, Wed>

regioneR
--------

Changes in version 1.1.8:

NEW FEATURES

- Added new functionality to permTest to use multiple evaluation
  functions with a single randomization procedure. This gives a
  significant speedup when comparing a single region set with multiple
  other features

- Created a new function createFunctionsList() that given a function
  and a list of values, creates a list of curried functions (e.g with
  one parameter preassigned to each of the given values)

PERFORMANCE IMPROVEMENTS

- Complete rewrite of randomizeRegions() resulting in a 10 to 100 fold
  speedup

BUG FIXES

- Multiple minor bug fixes

regionReport
------------

Changes in version 1.3.8:

SIGNIFICANT USER-VISIBLE CHANGES

- renderReport() and derfinderReport() now show Manhattan plots for
  p-value variables (p-value, q-value, FWER adjusted p-value).

Changes in version 1.3.7:

NEW FEATURES

- renderReport() now has the 'densityTemplates' argument via which
  users can customize the density plots for the p-value variables and
  the continuous variables. This addresses one of David Robinson's
  requests at http://f1000research.com/articles/4-105/v1

Changes in version 1.3.6:

NEW FEATURES

- Added a vignette with an example report from bumphunter results.

Changes in version 1.3.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Merged pull request https://github.com/leekgroup/regionReport/pull/7

- Added "template" argument to renderReport and derfinderReport to
  customize the knitr template used

- Wrapped code that works in a temporary directory in with_wd function,
  which evaluates in the directory but returns to the original
  directory in the case of a user interrupt or error (with on.exit)

Changes in version 1.3.4:

NEW FEATURES

- Reports now have a link to the BibTeX file used for the references.
  This addresses http://f1000research.com/articles/4-105/v1#reflist
  Karthik Ram's bullet point number 4.

Changes in version 1.3.3:

NEW FEATURES

- Now uses derfinderPlot::vennRegions() to show venn diagram of genomic
  states. Requires derfinderPlot version 1.3.2 or greater.

- derfinderReport() now has a 'significantVar' argument that allows
  users to choose between determining significant regions by P-values,
  FDR adjusted P-values, or FWER adjusted P-values (if FWER adjusted
  P-values are absent, then FDR adjusted P-values are used instead,
  with a warning).

Changes in version 1.3.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Deprecated functions with underscores in their names in favor of
  camelCase functions. This was done to simplify the package.

Changes in version 1.3.1:

BUG FIXES

- Fixed renderReport() and derfinderReport() so they'll open the
  correct URL when interactive() == TRUE and the user has
  knitrBootstrap version 0.9.0 installed instead of the latest GitHub
  version.


rGREAT
------

Changes in version 1.1.6:

- versions of GREAT can be selected

Changes in version 1.1.5:

- minor changes in documentations

Changes in version 1.1.4:

- minor changes to adjust to github mirror

Changes in version 1.1.3:

- error message is extract from GREAT now.

- default value of `bgChoice` depeneds on `gr`

Changes in version 1.1.2:

- add dependency to base packages

Changes in version 1.1.1:

- enrivonment for the slots in `GreatJob` object are initialized inside
  `submitGreatJob`

rhdf5
-----

Changes in version 2.14.0:

NEW FEATURES

- improved handling of error messages: HDF5 error messages are
  simplified and forwarded to R.

- When reading integer valued data, especially 64-integers and unsigned
  32-bit integers, overflow values are now replaced by NA's and a
  warning is thrown in this case.

- When coercing HDF5-integers to R-double, a warning is displayed when
  integer precision is lost.

- New low level general library function H5Dget_storage_size
  implemented.

BUG FIXES

- Memory allocation on heap instead of stack for reading large datasets
  (Thanks to a patch from Jimmy Jia).

- Some bugs have been fixed for reading large 64-bit integers and
  unsigned 32-bit integers.

- A bug was fixed for reading HDF5 files containing soft links.

RiboProfiling
-------------

Changes in version 0.99.7:

- In riboSeqFromBAM function capture.output called with utils::

Changes in version 0.99.6:

- In codonPCA functions prcomp and kmeans called with stats::

Changes in version 0.99.5:

- Added @return for roxygen comment in printPCA.R

- Updated packages IRanges and BSgenome.Mmusculus.UCSC.mm10 before
  check.

Changes in version 0.99.4:

- Non-null coverage values define the percentage of best expressed
  CDSs.

- New function: readsToReadStart - it builds the GRanges object of the
  read start genomic positions

- Acronym BAM used instead of bam

- Correction of read start coverage (readStartCov) for the reverse
  strand

- Title in Description in Title Case

- Typos corrections in vignette.

- Style inconsistencies solved: = replaced by <- outside named
  arguments No space around “=” when using named arguments to
  functions. This: somefunc(a=1, b=2) spaces around binary operators a
  space after all commas use of camelCase for both variable and
  function names ORFrelativePos -> orfRelativePos

- Replaced 1:length(x) by seq_len(length(x))

- Replaced 1:nrow(x) by seq_len(NROW(x))

- Replaced trailing white spaces with this command: 'find . -type f
  -path './*R' -exec perl -i -pe 's/ +$//' {} \;'

- countsPlot, histMatchLength, plotSummarizedCov: no longer print
  directly the graphs. Instead they return a list of graphs.

- Replaced 'class()' tests by 'inherits' or 'is'.

- The codonPCA function no longer prints the PCA graphs sequantially.
  The 5 PCA graphs are returned, together with the PCA scores.

- New function printPCA prints the 5 PCA plots produced by codonPCA.

- A BAM file is now available in the inst/extdata ctrl_sample.bam.  It
  is used in the testriboSeqFromBAM testthat

Changes in version 0.99.3:

- Correction of a infinite loop in function riboSeq_fromBam

Changes in version 0.99.2:

- Modified the vignette: small corrections of the explanatory text.

- Added testthat tests for the following functions:
  test-aroundPromoter.R, test-riboSeq_fromBAM.R, test-readStartCov.R

- Introduced the 4 spaces tabulation

- Reduced the percentage of lines > 80 characters to 2%

- Added imports in the Namespace for the proposed packages or methods

- New R version: 3.2.2 and bioconductor packages update

Changes in version 0.99.1:

- Added biocViews: Sequencing, Coverage, Alignment, QualityControl,
  Software, PrincipalComponent

- Added testthat tests for the following functions:
  test-aroundPromoter.R, test-riboSeq_fromBAM.R, test-readStartCov.R

- Introduced the 4 spaces tabulation

- Reduced the percentage of lines > 80 characters to 2%

- Added imports in the Namespace for the proposed packages or methods
  v.0.99.0 Initial release.

RNAprobR
--------

Changes in version 1.1.2:

BUG FIXES

- Fixed a bug in bedgraph2norm() by adding 'ignore.strand = "FALSE"'
  the absence of which could produce erroneous results.

MISC

- Added CITATION

RnBeads
-------

Changes in version 1.1.9:

- Filtering report fix when no normalization is conducted

- Bugfixes for combine(RnBSet, RnBSet) and BigFf matrices

Changes in version 1.1.8:

- Corrected coverage statistics in sample summary table: Sites with NA
  methylation values are no longer considered in the coverage
  statistics (makes a difference if some coverage threshold is applied)

- Improved method for gender prediction. Predicted genders are also
  included in the exported annotation table

Changes in version 1.1.7:

- Improvements to mergeSamples function for RnBiseqSets

- Some more memory clean-up

Changes in version 1.1.6:

- Differential methylation based on region level only is now supported

- Minor updates to the differential methylation report generation

- Performance improvements and minor bugfixes for using disk.dump.bigff

- Performance improvements (more memory clean-up)

Changes in version 1.1.5:

- New option: disk.dump.bigff for using a wrapper class of the ff
  package to avoid their INT_MAX issue

Changes in version 1.1.4:

- Fixes in parallel environment setup from bioconcuctor

Changes in version 1.1.3:

- Some fixes in data loading

- Fixes in parallel environment setup

Changes in version 1.1.2:

- New annotation package format

- Support for filtering out cross-reactive probes in Infinium 450k
  dataset

- Improved logging on a Mac

rols
----

Changes in version 1.11.6:

- fix failing unit tests <2015-07-02 Thu>

Changes in version 1.11.5:

- add test script <2015-06-30 Tue>

- more unit tests <2015-06-30 Tue>

Changes in version 1.11.4:

- Fix bug in as("[MS, MS:123, ]", "CVParam") and add unit test
  <2015-06-19 Fri>

Changes in version 1.11.3:

- New charIsCVParam function to check if a character represents a valid
  CV param <2015-06-16 Tue>

Changes in version 1.11.2:

- Unit test <2015-06-15 Mon>

Changes in version 1.11.1:

- New char to CVParam coerce method <2015-06-15 Mon>

Changes in version 1.11.0:

- Bioc devel 3.2

ropls
-----

Changes in version 1.1.11:

NEW FEATURES

- NEWS format updated

Changes in version 1.1.10:

NEW FEATURES

- opls: PLS models can now be built for 'x' data with a single variable

Changes in version 1.1.9:

NEW FEATURES

- opls: default number of permutations set to 10 (instead of 100) as a
  compromise to enable both quick computation and a first hint at model
  significance

- opls: maximum number of components in automated mode (predI = NA or
  orthoI = NA) set to 10 (instead of 15)

- plot.opls: bug corrected in case of single component model without
  permutation testing

Changes in version 1.1.8:

SIGNIFICANT USER-VISIBLE CHANGES

- sacurine dataset name simplification: names of sample metadata are
  now 'age', 'bmi', and 'gender' instead of 'ageVn', 'bmiVn', and
  'genderFc'; names of variableMetadata are now 'msiLevel', 'hmdb',
  'chemicalClass' instead of 'msiLevelVn', 'hmdbVc', 'chemicalClassVc'

Changes in version 1.1.7:

NEW FEATURES

- Control to avoid overfitting was strenghtened by: i) setting the
  default number of permutations to 100 (instead of 0) and ii) changing
  the default plot ("summary") to include both the "permutation" and
  the "overview" graphics

Changes in version 1.1.6:

NEW FEATURES

- .sinkC argument added in the opls and plot.opls methods: Diversion of
  messages is required for Galaxy integration

Changes in version 1.1.5:

NEW FEATURES

- current implementation supports two-class only classification

Changes in version 1.1.4:

NEW FEATURES

- renamed slot of opls object: rotationMN -> weightStarMN

Changes in version 1.1.3:

NEW FEATURES

- Minor internal changes regarding default ellipse plotting options and
  checking the compatibility between the subset size and the
  cross-validation fold

Changes in version 1.1.2:

BUG FIXES

- Minor internal changes

Changes in version 1.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- The packaging was modified (but not the algorithms) to be consistent
  with other machine learning packages: 'opls' is now a class and the
  'print', 'plot', 'predict', 'summary', 'fitted', 'coefficients' and
  'residuals' methods are available (see the vignette)

- renamed method: roplsF -> opls

- renamed arguments testVi -> now 'subset' which indicates the indices
  of the training (instead of the testing) observations

- values: tMN -> scoreMN pMN -> loadingMN wMN -> weightMN bMN ->
  coefficients rMN -> rotationMN varVn -> pcaVarVn tOrthoMN ->
  orthoScoreMN pOrthoMN -> orthoLoadingMN wOrthoMN -> orthoWeightMN

- new (S3) methods for objects of class 'opls' print summary plot
  predict fitted coefficients residuals


rpx
---

Changes in version 1.5.1:

- unit tests <2015-06-30 Tue>

- export and document pxnodes <2015-06-30 Tue>

Changes in version 1.5.0:

- Bioc version 3.2 (devel)

Rqc
---

Changes in version 1.4:

NEW FEATURES

- Custom template support added

- Read frequency result and plot (rqcReadFrequencyPlot) added

- Per file top represented reads added

- Top representated reads added

- Per file heatmap plot (rqcFileHeatmap) added

- checkpoint function added (experimental)

USER VISIBLE CHANGES

- BPPARAM argument replaced by workers argument

- File information table added to default report template

- Function rqcReadWidthPlot, y-axis changed to proportion (%)

- Almost all plots use colorblid scheme

Rsamtools
---------

Changes in version 1.21:

SIGNIFICANT USER-VISIBLE CHANGES

- pileup adds query_bins arg to give strand-sensitive cycle bin
  behavior; cycle_bins renamed left_bins; negative values allowed
  (including -Inf) to specify bins based on distance from end-of-read.

- mapqFilter allows specification of a mapping quality filter threshold

- PileupParam() now correctly follows samtools with
  min_base_quality=13, min_map_quality=0 (previously, values were
  assigned as 0 and 13, respectively)

- Support parsing 'B' tags in bam file headers.

BUG FIXES

- segfault on range iteration introduced 1.19.35, fixed in 1.21.1

- BamViews parallel evaluation with BatchJobs back-end requires named
  arguments


Rsubread
--------

Changes in version 1.20.0:

NEW FEATURES

- Fast sorting of input bam files in featureCounts.

- Fractional counting of multi-mapping reads in featureCounts.

- Detection of complex indels in Subread and Subjunc aligners.

- Including more candidate locations in read re-alignment step to
  improve mapping performance.

- New formula for mapping paired-end reads that takes into account
  paired-end distance, number of subread votes and number of mismatched
  bases.


rtracklayer
-----------

Changes in version 1.30:

NEW FEATURES

- Add readGFF(), a fast and flexible GFF/GTF parser implemented in C.
  See ?readGFF for more information.

SIGNIFICANT USER-VISIBLE CHANGES

- import.gff() now uses readGFF() internally which makes it at least 5x
  faster in most use cases and dramatically reduces its memory
  footprint.

RUVSeq
------

Changes in version 1.3:

- New argument 'isLog' to RUV* methods: if counts provided on the log
  scale, normalizedCounts will also be on the log scale.

- Fixed a bug: epsilon is now removed from the corrected counts.

- New function makeGroups to make scIdx matrix for RUVs.

sbgr
----

Changes in version 1.0.0 (2015-05-01):

- Initial version

- All the APIs of the SBG platform are supported

- First vignette added

seq2pathway
-----------

Changes in version 1.1.8:

- Used Sys.info() rather than sessionInfor() to extract platform
  information

- Avoided the use of "\\" in R

- Added a parameter "Ontology" to the function FisherTest_GO_BP_MF_CC()

Changes in version 1.1.6:

- Corrected the python path for Linux users

- Activated the demo code in R help files and vignettes

- Changed the maintainer

Changes in version 1.1.4:

- Update the Fisher's exact test used mouse gene background.

- Replace the mm10 GENCODE database V3 with V4 in the seq2pathway.data

Changes in version 1.1.2:

- add citation

Changes in version 1.0.2:

- add parameter UTR3 in function runseq2gene

SeqArray
--------

Changes in version 1.10.0:

- The version number was bumped for the Bioconductor release version
  3.2


SeqVarTools
-----------

Changes in version 1.7.9:

- added more options to return dosage of different alleles

Changes in version 1.7.7:

- added methods to calculate and plot reference allele fraction

Changes in version 1.7.6:

- HWE method returns additional columns and allows permutation of
  genotypes

Changes in version 1.7.5:

- Added SeqVarData class to combine sample annotation with GDS object

Changes in version 1.7.4:

- HWE works on biallelic INDELs as well as SNVs

Changes in version 1.7.3:

- Add methods for duplicateDiscordance with two datasets

- Add alternateAlleleDetection

SGSeq
-----

Changes in version 1.4.0:

- Added importTranscripts() for importing annotation from GFF format

- Added plotCoverage() for visualization of per-base read coverage and
  junction read counts

- Added predictVariantEffects() for predicting the effect of splice
  variants on annotated protein-coding transcripts

- findSGVariants() is now able to deal with more complex gene models

- SGVariants columns closed5p and closed3p now refer to individual
  splice variants rather than the splice event they belong to

- Bug fixes and other improvements

ShortRead
---------

Changes in version 1.27:

SIGNIFICANT USER-VISIBLE CHANGES

- fastqFilter allows several input 'files' to be written to a single
  'destinations'.

- readAligned() for BAM files is defunct. QA and associated methods
  removed.

- srapply removed


sincell
-------

Changes in version 1.1.01:

- paramether tsne.theta has been added to function
  sc_DimensionalityReductionObj

- function sc_InSilicoCellsReplicatesObj has been expanded to
  incorporate a forth model of noise generation following the negative
  binomial (NB) distribution. For each gene, dispersion parameters for
  the NB can be estimated in three alternative ways allowing different
  types of estimates (e.g. technical noise provided by the user). See
  documentation of sc_InSilicoCellsReplicatesObj() function.

SNPhood
-------

Changes in version 0.99.0 (2015-08-07):

- various bugfixing


SNPRelate
---------

Changes in version 1.4.0:

- the version number was bumped for the Bioconductor release version
  3.2


specL
-----

Changes in version 1.3.7:

USER UNVISIBLE CHANGES

- getProteinPeptideTable - added

Changes in version 1.3.5:

USER UNVISIBLE CHANGES

- read.bibliospec - bugfixes

Changes in version 1.3.4:

USER VISIBLE CHANGES

- added Witold Wolski as maintainer

USER UNVISIBLE CHANGES

- read.bibliospec - replaced old code (for loop) by using mcmapply

- added time meassurements to read.bibliospec

Changes in version 1.3.3:

USER VISIBLE CHANGES

- plot::specLSet draws alpha circles iff plot(..., art=TRUE)

USER UNVISIBLE CHANGES

- .mascot2psmSet buxfix

- renamed column name in spectronaut outpu from irt to irt_or_rt

Changes in version 1.3.2:

USER VISIBLE CHANGES

- added ssrc (Sequence Specific Retention Calculator) function

- added a CITATION file

Changes in version 1.3.1:

USER VISIBLE CHANGES

- added fucntion cdsw

USER UNVISIBLE CHANGES

- modified unit test for genSwathIonLib

subSeq
------

Changes in version 1.0.0:

- Initial release of subSeq package to subsample (or downsample)
  RNA-Seq experiments.

supraHex
--------

Changes in version 1.7.3:

NEW FEATURES

- Allow the function "visHexMulComp" to have user-input the number of
  rows and columns for a rectangle grid wherein the planes are placed

Changes in version 1.7.2:

NEW FEATURES

- Add a new function "visHexAnimate" to animate multiple component
  planes of a supra-hexagonal grid


synapter
--------

Changes in version 1.11.2:

- fixing bug introduced in 1.11.1 when coercing from Synapter to MSnSet
  object - see this post on the support forum for details
  https://support.bioconductor.org/p/71087/ <2015-09-27 Sun>

Changes in version 1.11.1:

- avoiding error when some fcols are missing when coercing from
  Synapter to MSnSet object - see this post on the support forum for
  details https://support.bioconductor.org/p/71087/ <2015-09-07 Mon>

systemPipeR
-----------

Changes in version 1.3:

OVERVIEW

- systemPipeR is an R/Bioconductor package for building and running
  automated analysis workflows for a wide range of next generation
  sequence (NGS) applications. Important features include a uniform
  workflow interface across different NGS applications, automated
  report generation, and support for running both R and command-line
  software, such as NGS aligners or peak/variant callers, on local
  computers or compute clusters. Efficient handling of complex sample
  sets and experimental designs is facilitated by a consistently
  implemented sample annotation infrastructure.

- The most important enhancements in the upcoming release of the
  package are outlined below.

NGS WORKFLOWS

- Added new end-to-end workflows for 3 additional NGS application
  areas: - Ribo-Seq and polyRibo-Seq - ChIP-Seq - VAR-Seq The previous
  version of systemPipeR included only a complete workflow for RNA-Seq.

- Added the data package 'systemPipeRdata' to generate systemPipeR
  workflow environments with a single command (genWorkenvir) containing
  all parameter files and sample data required to quickly test and run
  workflows. This change will also allow evaluation of much more code
  examples in the vignettes during the package build/test process than
  this was possible in the past.

- About 20 new functions have been added to the package. Some examples
  are: - Read pre-processor function with support for SE and PE reads -
  Parallelization option of detailed FASTQ quality reports - Read
  distribution plots across all features available in a genome
  annotation (see ?featuretypeCounts) - Visualization of coverage
  trends along transcripts summarized for any number of transcripts
  (see ?featureCoverage) - Functionalities to predict uORFs/sORFs and
  to use them for expression profiling - Differential
  expression/binding analysis includes now DESeq2 as well as edgeR

- Added param templates for additional command-line software including,
  but not limited to: BWA-MEM, GATK, BCFtools, MACS2

- Adoption of R Markdown for main vignette. Future plans are to provide
  for all workflows the report templates in both formats: Latex/PDF and
  R_Markdown/HTML.

WORFLOW FRAMEWORK

- Simplified design of complex analysis workflows. Workflows can now
  include any number or combination of R and/or command-line steps

- Improvements to workflow automation and parallelization on single
  machines and computer clusters. This also includes now many
  additional parallelization examples in the workflow vignettes.

TargetSearch
------------

Changes in version 1.26.0:

SIGNIFICANT USER-VISIBLE CHANGES

- TargetSearch now depends on "ncdf" rather than mzR.

BUG FIXES

- Fix potential potencial pearson correlation errors that occur when
  the standard is deviation zero. If that is the case, replace the
  resulting NA/NaN by zero.

- Fix R check warnings.

TarSeqQC
--------

Changes in version 0.99.9:

CODE

- Modifications in pileupCounts in order to consider stranded counts

- Changes in plot methods. arrangeGrob from gridExtra was changed by
  plot_grid form cowplot package.

- Update the package vignette.

Changes in version 0.99.8:

CODE

- Modifications in pileupCounts and buildFeaturePanel in order to
  process overlapped features.

- Wrap plotting examples.

Changes in version 0.99.7:

CODE

- Modifications in pileupCounts and buildFeaturePanel in order to
  process overlapped features.

- Modifications in bedFile building in order to remove duplicated
  features based on duplicated start, end and chromosome definitions.

Changes in version 0.99.6:

CODE

- Definition of pileupCounts as a function instead a TargetExperiment
  S4 method.

- Adaptation and optimization of pileupCounts and bplapply usage in the
  buildFeaturePanel method.

Changes in version 0.99.4:

DESCRIPTION file

- Modification of the Rsatmools version dependency.

Changes in version 0.99.3:

CODE

- Modifications in buildFeaturePanel and summarizePanel for a correct
  implementation of parallel computing.

Changes in version 0.99.2:

CODE

- BiocParallel was used instead of parallel package for parallel
  computing.

Changes in version 0.99.1:

DOCUMENTATION

- The example data has been compiled using Bioconductor 3.2 packages.

Changes in version 0.99.0:

DOCUMENTATION

- `NEWS` file was added.

TCC
---

Changes in version 1.9.3:

- changed default value. Set 'test.method = "DESeq2" as default value
  when analyze multi-group without replicate data.

- add 'makeFCMatrix' function for generating the foldchange matrix that
  is used in 'simulateReadCounts' funtion.

- add function to simulate DEGs using fondchange matrix into
  'simulateReadCounts' function.

TEQC
----

Changes in version 3.9.1:

- new parameter 'plotchroms' in function 'chrom.barplot' that allows to
  specify the chromosomes (and their desired order) that shall be
  included in the plot

- bug fix regarding use of 'Offset' bases in 'TEQCreport'

TFBSTools
---------

Changes in version 1.7.2:

NEW FEATURES

- New class TFFMFirst and TFFMDetail for next generation TFBSs.

- Novel TFFM sequence logo.

BUG FIXES

- Fix the bug in runMEME, which always return positive strand for site
  sequence.

TPP
---

Changes in version 1.9.9:

- Major update to the CCR-part: It is now possible to fit and plot
  multiple experiments simultaneously.

- It is now possible to perform user-specified comparisons of different
  experiments. They are specified in the 'comparison' column of the
  config table.

- TR-part: Hypothesis testing is now separated from result table
  creation. Therefore, a new function was introduced
  (tpptrAnalyzeMeltCurves) -> see vignette.

- CCR-part: Curve fitting and plotting are now conducted by separate
  functions -> see vignette.

- Bug fixes

- Introduced color coding of the columns belonging to different
  experiments in the excel output. This requires openxlsx version >=
  2.4.0.

- The CCR workflow now only returns normalized measurements if
  normalization was actually performed. Unmodified measurements are
  always returned and indicated by the suffix 'unmodified'.

- Data import from tab-delimited files now ignores quotes so that
  protein annotation fields can contain single ' or " characters.

- Now enabling arbitrary numbers of plot colors for melting curves or
  dose response curves.

Changes in version 1.2.5:

- CCR-part: fold changes can now be normalized to their lowest
  temperature during import by the function tppccrNormalizeToReference.
  This ensures that the transformed values will always be between 0 and
  1.  When normalizing to the lowest concentration, and this value is
  0, a small constant (1e-15) is added to prevent division by zero.

- CCR-import: Argument nonZeroCols can be NULL if no additional
  filtering is needed.

- CCR output: To make filtering easier, the 'passed_filter' column now
  shows FALSE even if proteins could not be used for fitting (instead
  of NAs).

- CCR output: To make filtering easier, the 'passed_filter' column now
  shows FALSE even if proteins could not be used for fitting (instead
  of NAs).

- CCR output: column with normalization results was renamed from
  'normalized' to 'median_normalized' to distinguish from the newly
  introduced normalization to value at lowest concentration.

- Excel export: Columns are only color-coded by experiment, when the
  number of experiments is > 1.

- Excel export: Columns are only color-coded by experiment when the
  number of experiments is > 1.

- Excel export: Relative paths to the plots work now for TR- and CCR
  part

- Bug fix for data import: Unique identifiers now treated correctly
  again.

- Bug fix for Excel output: Boolean column entries are only transformed
  to "yes"/"no" for non-missing values.

- Bugfix in TR-QC plots: do not attempt to create Tm difference
  histograms if only one experiment is provided.

- Bugfix in TR normalization: fixedReference argument works again

- Bugfix in TR-QC plots: do not attempt to plot minSlope vs. Tm-diffs
  if no valid Tm-diff values are available.

Changes in version 1.1.4:

- Curve parameters "a" and "b" are now reported in the output of
  analyzeTPPTR and tpptrCurveFit.

Changes in version 1.1.3:

- fixed problem in QC-plot generation due to the gridExtra package
  update: adapted height parameter in grid.arrange

Changes in version 1.1.2:

- the recent update of the gridExtra package to Version 2.0.0 required
  a small fix to ensure sucessful plotting of the melting curves and
  parameter tables.

Changes in version 1.1.1:

- bugfixes in config table import from .csv files: 1. can now handle
  different types of delimitors (";", "\t", ",") 2. automatically
  removes empty rows after import

Changes in version 1.1.0:

- Changes in devel will be characterized by version number 1.1.x

trackViewer
-----------

Changes in version 1.5.6:

- update documentation.

- add new feature for optimizing the styles with theme.

Changes in version 1.5.5:

NEW FEATURES

- remove interactiveViewer from the package.

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.5.4:

NEW FEATURES

- parse2GRanges function to parse text into GRanges

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.5.3:

NEW FEATURES

- export importData function to import data to RleList.

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.5.2:

NEW FEATURES

- Add GRanges operators: +, -, *, /

- export parseWIG function.

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

trigger
-------

Changes in version 1.15.1:

- Minor changes to fixed check error

TRONCO
------

Changes in version 2.0.0-16:

- Release

variancePartition
-----------------

Changes in version 0.99.8:

- improve warnings / errors when design matrix is close to or exactly
  singular

Changes in version 0.99.7:

- added new class varPartResults to store results of
  fitExtractVarPartModel() and fitExtractVarPartModel() - the user will
  not notice any change, only the backend is different

- Allow computation of adjusted ICC in addition to ICC.

- add warning when categorical variables are modeled as fixed effects

- fix computation of variance fractions for varying coefficient models

- add getVarianceComponents() to return variances from lmer() or lm()
  model fit

- showWarnings=FALSE suppresses warning messages

- add fxn argument to fitVarPartModel to evaluate any function on the
  model fit

Changes in version 0.99.6:

- Update DESCRIPTION information

Changes in version 0.99.5:

- residuals deals with missing data gracefully and returns a matrix

Changes in version 0.99.4:

- add documentation for example datasets

- convert calcVarPart() to S4 from S3 function call

- fix typos in vignette

Changes in version 0.99.3:

- fitVarPartModel() and fitExractVarPartModel() use S4 instead of S3
  calls

Changes in version 0.99.2:

- rename sort.varParFrac to sortCols

- support ExpressionSet

- change options for plotStratifyBy() # Before Bioconductor submission

Changes in version 0.99.0:

- Initial version


VariantAnnotation
-----------------

Changes in version 1.16.0:

NEW FEATURES

- support REF and ALT values ".", "+" and "-" in predictCoding()

- return non-translated characters in VARCODON in predictCoding()
  output

- add 'verbose' option to readVcf() and friends

- writeVcf() writes 'fileformat' header line always

- readVcf() converts REF and ALT values "*" and "I" to '' and '.'

MODIFICATIONS

- VRanges uses '*' strand by default

- coerce 'alt' to DNStringSet for predictCoding,VRanges-method

- add detail to documentation for 'ignore.strand' in predictCoding()

- be robust to single requrested INFO column not present in vcf file

- replace old SummarizedExperiment class from GenomicRanges with the
  new new RangedSummarizedExperiment from SummarizedExperiment package

- return strand of 'subject' for intronic variants in locateVariants()

BUG FIXES

- writeVcf() does not duplicate header lines when chunking

- remove extra tab after INFO when no FORMAT data are present

- filteVcf() supports 'param' with ranges

VariantFiltering
----------------

Changes in version 1.6:

USER VISIBLE CHANGES

- Update on the scores() method for PhastConsDb objects that enables a
  10-fold faster retrieval of mean phastCons scores over genomic
  intervals.

- Added two new annotated regions coded as fiveSpliceSite and
  threeSpliceSite and the scoring of their binding affyinity, if
  scoring matrices are provided.

BUG FIXES

- Fix on the scores() method for PhastConsDb objects, that was
  affecting multiple-nucleotide ranges from an input GRanges object
  with unordered sequence names.

- Fix on snpid2maf() method when decoding MafDb variants with AF values
  whose significant digits start with 95.

xcms
----

Changes in version 1.45.7:

USER VISIBLE CHANGES

- Disabled Rmpi support and usage on Windows

Changes in version 1.45.6:

NEW FEATURE

- J. Rainer implemented a [ method that allows to subset an xcmsSet.

BUG FIXES

- Fixed a problem in split.xcmsSet that did not split the phenoData
  properly. Added some details to the documentation of xcmsSet-class.

Changes in version 1.45.5:

USER VISIBLE CHANGES

- The sampclass method for xcmsSet will now return the content of the
  column "class" from the data.frame in the phenoData slot, or if not
  present, the interaction of all factors (columns) of that data.frame.

- The sampclass<- method replaces the content of the "class" column in
  the phenoData data.frame. If a data.frame is submitted, the
  interaction of its columns is calculated and stored into the "class"
  column.

BUG FIXES

- Fixed a bug that resulted in a cryptic error message when no input
  files are available to the xcmsSet function.

Changes in version 1.45.4:

BUG FIXES

- Fixed a bug in the levelplot method for xcmsSet.

Changes in version 1.45.3:

NEW FEATURE

- xcmsSet now allows phenoData to be an AnnotatedDataFrame.

- new slots for xcmsRaw: - mslevel: store the mslevel parameter
  submitted to xcmsRaw. - scanrange: store the scanrange parameter
  submitted to xcmsRaw.

- new slots for xcmsSet: - mslevel: stores the mslevel argument from
  the xcmsSet method. - scanrange: to keep track of the scanrange
  argument of the xcmsSet method.

- new methods for xcmsRaw: - levelplot: similar to the image method,
  plots m/z vs RT with color coded intensities. - mslevel: returns the
  value for the .mslevel slot. For downstream compatibility, this
  method returns NULL if the object does not have the same named slot.
  - profinfo: same functionality as the profinfo method for xcmsSet. -
  scanrange: returns the value for the scanrange slot. For downstream
  compatibility, this method returns NULL if the object does not have
  the same named slot.

- new methods for xcmsSet: - getXcmsRaw: returns a xcmsRaw object for
  one or more files in the xcmsSet, eventually applying retention time
  correction etc. - levelplot: similar to the image method, plots m/z
  vs RT with color coded intensities. Allows in addition to highlight
  identified peaks. - mslevel: returns the value for the mslevel slot.
  For downstream compatibility, this method returns NULL if the object
  does not have the same named slot. - profMethod: same functionality
  as the profMethod method of xcmsRaw. - profStep: same functionality
  as the profStep method of xcmsRaw. - scanrange: returns the value for
  the scanrange slot. For downstream compatibility, this method returns
  NULL if the object does not have the same named slot.

USER VISIBLE CHANGES

- show method for xcmsSet updated to display also informations about
  the mslevel and scanrange.

- Elaborated some documentation entries.

- rtrange and mzrange for xcmsRaw method plotEIC use by default the
  full RT and m/z range.

- Added arguments "lty" and "add" to plotEIC method for xcmsRaw.

- getEIC without specifying mzrange returns the ion chromatogram for
  the full m/z range (i.e. the base peak chromatogram).

BUG FIXES

- Checking if phenoData is a data.frame or AnnotatedDataFrame and throw
  an error otherwise.

- xcmsSet getEIC method for water Lock mass corrected files for a
  subset of files did not evaluate whether the specified files were
  corrected.

Changes in version 1.45.2:

BUG FIXES

- The xcms split() function now accepts factors that are shorter than
  the number of samples in the xcmsSet, following more closely the
  standard split() behaviour

Changes in version 1.45.1:

NEW FEATURE

- plotrt now allows col to be a vector of color definition, same as the
  plots for retcor methods.

- Added $ method to access phenoData columns in a eSet/ExpressionSet
  like manner.

- Allow to use the "parallel" package for parallel processing of the
  functions xcmsSet and fillPeaks.chrom.

- Thanks to J. Rainer!

xps
---

Changes in version 3.2:

VERSION xps-1.29.1

- update README file


yaqcaffy
--------

Changes in version 1.29.1:

- Update email address <2015-05-26 Tue>


Packages removed since the last release
=================================

No packages were removed in this release.