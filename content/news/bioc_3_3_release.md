April 4, 2016

Bioconductors:

We are pleased to announce Bioconductor 3.3, consisting of 1211
software packages, 293 experiment data packages, and 916
up-to-date annotation packages.

There are 107 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.3 is compatible with R 3.3,
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

* Getting Started with Bioconductor 3.3
* New Software Packages
* NEWS from new and existing packages
* Packages removed from Bioconductor since the last release

Getting Started with Bioconductor 3.3
======================================

To update to or install Bioconductor 3.3:

1. Install R 3.3.  Bioconductor 3.3 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/) .

New Software Packages
=====================

There are 107 new packages in this release of Bioconductor.

AneuFinder - This package implements functions for CNV calling, plotting, export and analysis from whole-genome single cell sequencing data.

bacon - Bacon can be used to remove inflation and bias often observed in epigenome- and transcriptome-wide association studies. To this end bacon constructs an empirical null distribution using a Gibbs Sampling algorithm by fitting a three-component normal mixture on z-scores.

BadRegionFinder - BadRegionFinder is a package for identifying regions with a bad, acceptable and good coverage in sequence alignment data available as bam files. The whole genome may be considered as well as a set of target regions. Various visual and textual types of output are available.

BasicSTARRseq - Basic peak calling on STARR-seq data based on a method introduced in "Genome-Wide Quantitative Enhancer Activity Maps Identified by STARR-seq" Arnold et al. Science. 2013 Mar 1;339(6123):1074-7. doi: 10.1126/science. 1232542. Epub 2013 Jan 17.

BatchQC - Sequencing and microarray samples often are collected or processed in multiple batches or at different times. This often produces technical biases that can lead to incorrect results in the downstream analysis. BatchQC is a software tool that streamlines batch preprocessing and evaluation by providing interactive diagnostics, visualizations, and statistical analyses to explore the extent to which batch variation impacts the data. BatchQC diagnostics help determine whether batch adjustment needs to be done, and how correction should be applied before proceeding with a downstream analysis. Moreover, BatchQC interactively applies multiple common batch effect approaches to the data, and the user can quickly see the benefits of each method. BatchQC is developed as a Shiny App. The output is organized into multiple tabs, and each tab features an important part of the batch effect analysis and visualization of the data. The BatchQC interface has the following analysis groups: Summary, Differential Expression, Median Correlations, Heatmaps, Circular Dendrogram, PCA Analysis, Shape, ComBat and SVA.

BgeeDB - A package for the annotation and gene expression data download from Bgee database, and TopAnat analysis: GO-like enrichment of anatomical terms, mapped to genes by expression patterns.

biomformat - This is an R package for interfacing with the BIOM format. This package includes basic tools for reading biom-format files, accessing and subsetting data tables from a biom object (which is more complex than a single table), as well as limited support for writing a biom-object back to a biom-format file. The design of this API is intended to match the python API and other tools included with the biom-format project, but with a decidedly "R flavor" that should be familiar to R users. This includes S4 classes and methods, as well as extensions of common core functions/methods.

BioQC - BioQC performs quality control of high-throughput expression data based on tissue gene signatures

biosigner - Feature selection is critical in omics data analysis to extract restricted and meaningful molecular signatures from complex and high-dimension data, and to build robust classifiers. This package implements a new method to assess the relevance of the variables for the prediction performances of the classifier. The approach can be run in parallel with the PLS-DA, Random Forest, and SVM binary classifiers. The signatures and the corresponding 'restricted' models are returned, enabling future predictions on new datasets. A Galaxy implementation of the package is available within the Workflow4metabolomics.org online infrastructure for computational metabolomics.

cellity - A support vector machine approach to identifying and filtering low quality cells from single-cell RNA-seq datasets.

cellTree - This packages computes a Latent Dirichlet Allocation (LDA) model of single-cell RNA-seq data and builds a compact tree modelling the relationship between individual cells over time or space.

Chicago - A pipeline for analysing Capture Hi-C data.

chromPlot - Package designed to visualize genomic data along the chromosomes, where the vertical chromosomes are sorted by number, with sex chromosomes at the end.

CHRONOS - A package used for efficient unraveling of the inherent dynamic properties of pathways. MicroRNA-mediated subpathway topologies are extracted and evaluated by exploiting the temporal transition and the fold change activity of the linked genes/microRNAs.

CINdex - The CINdex package addresses important area of high-throughput genomic analysis. It allows the automated processing and analysis of the experimental DNA copy number data generated by Affymetrix SNP 6.0 arrays or similar high throughput technologies. It calculates the chromosome instability (CIN) index that allows to quantitatively characterize genome-wide DNA copy number alterations as a measure of chromosomal instability. This package calculates not only overall genomic instability, but also instability in terms of copy number gains and losses separately at the chromosome and cytoband level.

clustComp - clustComp is a package that implements several techniques for the comparison and visualisation of relationships between different clustering results, either flat versus flat or hierarchical versus flat. These relationships among clusters are displayed using a weighted bi-graph, in which the nodes represent the clusters and the edges connect pairs of nodes with non-empty intersection; the weight of each edge is the number of elements in that intersection and is displayed through the edge thickness. The best layout of the bi-graph is provided by the barycentre algorithm, which minimises the weighted number of crossings. In the case of comparing a hierarchical and a non-hierarchical clustering, the dendrogram is pruned at different heights, selected by exploring the tree by depth-first search, starting at the root. Branches are decided to be split according to the value of a scoring function, that can be based either on the aesthetics of the bi-graph or on the mutual information between the hierarchical and the flat clusterings. A mapping between groups of clusters from each side is constructed with a greedy algorithm, and can be additionally visualised.

ClusterSignificance - The ClusterSignificance package provides tools to assess if clusters have a separation different from random or permuted data. ClusterSignificance investigates clusters of two or more groups by first, projecting all points onto a one dimensional line. Cluster separations are then scored and the probability of the seen separation being due to chance is evaluated using a permutation method.

CONFESS - Single Cell Fluidigm Spot Detector.

consensusSeekeR - This package compares genomic positions and genomic ranges from multiple experiments to extract common regions. The size of the analyzed region is adjustable as well as the number of experiences in which a feature must be present in a potential region to tag this region as a consensus region.

contiBAIT - Using strand inheritance data from multiple single cells from the organism whose genome is to be assembled, contiBAIT can cluster unbridged contigs together into putative chromosomes, and order the contigs within those chromosomes.

CountClust - Fits grade of membership models (GoM, also known as admixture models) to cluster RNA-seq gene expression count data, identifies characteristic genes driving cluster memberships, and provides a visual summary of the cluster memberships.

CrispRVariants - CrispRVariants provides tools for analysing the results of a CRISPR-Cas9 mutagenesis sequencing experiment, or other sequencing experiments where variants within a given region are of interest. These tools allow users to localize variant allele combinations with respect to any genomic location (e.g. the Cas9 cut site), plot allele combinations and calculate mutation rates with flexible filtering of unrelated variants.

dada2 - The dada2 package provides "OTU picking" functionality, but instead of picking OTUs the DADA2 algorithm exactly infers samples sequences. The dada2 pipeline starts from demultiplexed fastq files, and outputs inferred sample sequences and associated abundances after removing substitution and chimeric errors. Taxonomic classification is also available via a native implementation of the RDP classifier method.

dcGSA - Distance-correlation based Gene Set Analysis for longitudinal gene expression profiles. In longitudinal studies, the gene expression profiles were collected at each visit from each subject and hence there are multiple measurements of the gene expression profiles for each subject. The dcGSA package could be used to assess the associations between gene sets and clinical outcomes of interest by fully taking advantage of the longitudinal nature of both the gene expression profiles and clinical outcomes.

debrowser - Bioinformatics platform containing interactive plots and tables for differential gene and region expression studies. Allows visualizing expression data much more deeply in an interactive and faster way. By changing the parameters, user can easily discover different parts of the data that like never have been done before. Manually creating and looking these plots takes time. With this system users can prepare plots without writing any code. Differential expression, PCA and clustering analysis are made on site and the results are shown in various plots such as scatter, bar, box, volcano, ma plots and Heatmaps.

DEFormats - Covert between different data formats used by differential gene expression analysis tools.

diffloop - A suite of tools for subsetting, visualizing, annotating, and statistically analyzing the results of one or more ChIA-PET experiments.

DNAshapeR - DNAhapeR is an R/BioConductor package for ultra-fast, high-throughput predictions of DNA shape features. The package allows to predict, visualize and encode DNA shape features for statistical learning.

doppelgangR - The main function is doppelgangR(), which takes as minimal input a list of ExpressionSet object, and searches all list pairs for duplicated samples.  The search is based on the genomic data (exprs(eset)), phenotype/clinical data (pData(eset)), and "smoking guns" - supposedly unique identifiers found in pData(eset).

DRIMSeq - The package provides two frameworks. One for the differential splicing analysis between different conditions and one for the sQTL analysis. Both are based on modeling the counts of genomic features (i.e., transcripts, exons or exonic bins) with Dirichlet-multinomial distribution. The package also makes available functions for visualization and exploration of the data and results.

EBSEA - Calculates differential expression of genes based on exon counts of genes obtained from RNA-seq sequencing data.

EGAD - The package implements a series of highly efficient tools to calculate functional properties of networks based on guilt by association methods.

EGSEA - This package implements the Ensemble of Gene Set Enrichment Analyses (EGSEA) method for gene set testing.

EmpiricalBrownsMethod - Combining P-values from multiple statistical tests is common in bioinformatics. However, this procedure is non-trivial for dependent P-values. This package implements an empirical adaptation of Brown’s Method (an extension of Fisher’s Method) for combining dependent P-values which is appropriate for highly correlated data sets found in high-throughput biological experiments.

epivizrData - Serve data from Bioconductor Objects through a WebSocket connection.

epivizrServer - This package provides objects to manage WebSocket connections to epiviz apps. Other epivizr package use this infrastructure.

epivizrStandalone - This package imports the epiviz visualization JavaScript app for genomic data interactive visualization. The 'epivizrServer' package is used to provide a web server running completely within R. This standalone version allows to browse arbitrary genomes through genome annotations provided by Bioconductor packages.

EWCE - Used to determine which cell types are enriched within gene lists. The package provides tools for testing enrichments within simple gene lists (such as human disease associated genes) and those resulting from differential expression studies. The package does not depend upon any particular Single Cell Transcriptome dataset and user defined datasets can be loaded in and used in the analyses.

ExpressionAtlas - This package is for searching for datasets in EMBL-EBI Expression Atlas, and downloading them into R for further analysis. Each Expression Atlas dataset is represented as a SimpleList object with one element per platform. Sequencing data is contained in a SummarizedExperiment object, while microarray data is contained in an ExpressionSet or MAList object.

FamAgg - Framework providing basic pedigree analysis and plotting utilities as well as a variety of methods to evaluate familial aggregation of traits in large pedigrees.

flowAI - The package is able to perform an automatic or interactive quality control on FCS data acquired using flow cytometry instruments. By evaluating three different properties: 1) flow rate, 2) signal acquisition, 3) dynamic range, the quality control enables the detection and removal of anomalies.

garfield - GARFIELD is a non-parametric functional enrichment analysis approach described in the paper GARFIELD: GWAS analysis of regulatory or functional information enrichment with LD correction. Briefly, it is a method that leverages GWAS findings with regulatory or functional annotations (primarily from ENCODE and Roadmap epigenomics data) to find features relevant to a phenotype of interest. It performs greedy pruning of GWAS SNPs (LD r2 > 0.1) and then annotates them based on functional information overlap. Next, it quantifies Fold Enrichment (FE) at various GWAS significance cutoffs and assesses them by permutation testing, while matching for minor allele frequency, distance to nearest transcription start site and number of LD proxies (r2 > 0.8).

genbankr - Reads Genbank files.

GenoGAM - This package allows statistical analysis of genome-wide data with smooth functions using generalized additive models based on the implementation from the R-package 'mgcv'. It provides methods for the statistical analysis of ChIP-Seq data including inference of protein occupancy, and pointwise and region-wise differential analysis. Estimation of dispersion and smoothing parameters is performed by cross-validation. Scaling of generalized additive model fitting to whole chromosomes is achieved by parallelization over overlapping genomic intervals.

genphen - Given a set of genetic polymorphisms in the form of single nucleotide poylmorphisms or single amino acid polymorphisms and a corresponding phenotype data, often we are interested to quantify their association such that we can identify the causal polymorphisms. Using statistical learning techniques such as random forests and support vector machines, this tool provides the means to estimate genotype-phenotype associations. It also provides visualization functions which enable the user to visually inspect the results of such genetic association study and conveniently select the genotypes which have the highest strenght ofassociation with the phenotype.

GenRank - Methods for ranking genes based on convergent evidence obtained from multiple independent evidence layers. This package adapts three methods that are popular for meta-analysis.

GenVisR - Produce highly customizable publication quality graphics for genomic data primarily at the cohort level.

ggcyto - With the dedicated fority method implemented for flowSet, ncdfFlowSet and GatingSet classes, both raw and gated flow cytometry data can be plotted directly with ggplot. ggcyto wrapper and some customed layers also make it easy to add gates and population statistics to the plot.

Glimma - This package generates interactive visualisations of RNA-sequencing data based on output from limma, edgeR or DESeq2. Interactions are built on top of popular static displays from the limma package, providing users with access to gene IDs and sample information. Plots are generated using d3.js and displayed in HTML pages.

globalSeq - The method may be conceptualised as a test of overall significance in regression analysis, where the response variable is overdispersed and the number of explanatory variables exceeds the sample size.

GMRP - Perform Mendelian randomization analysis of multiple SNPs to determine risk factors causing disease of study and to exclude confounding variabels and perform path analysis to construct path of risk factors to the disease.

GSALightning - GSALightning provides a fast implementation of permutation-based gene set analysis for two-sample problem. This package is particularly useful when testing simultaneously a large number of gene sets, or when a large number of permutations is necessary for more accurate p-values estimation.

Harman - Harman is a PCA and constrained optimisation based technique that maximises the removal of batch effects from datasets, with the constraint that the probability of overcorrection (i.e. removing genuine biological signal along with batch noise) is kept to a fraction which is set by the end-user.

HDF5Array - This package implements the HDF5Array class for convenient access and manipulation of HDF5 datasets. In order to reduce memory usage and optimize performance, operations on an HDF5Array object are either delayed or executed using a block processing mechanism. The delaying and block processing mechanisms are independent of the on-disk backend and implemented via the DelayedArray class. They even work on ordinary arrays where they can sometimes improve performance.

iCARE - An R package to compute Individualized Coherent Absolute Risk Estimators.

iCOBRA - This package provides functions for calculation and visualization of performance metrics for evaluation of ranking and binary classification (assignment) methods. It also contains a shiny application for interactive exploration of results.

IHW - Independent hypothesis weighting (IHW) is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg by assigning data-driven weights to each hypothesis. The input to IHW is a two-column table of p-values and covariates. The covariate can be any continuous-valued or categorical variable that is thought to be informative on the statistical properties of each hypothesis test, while it is independent of the p-value under the null hypothesis.

ImmuneSpaceR - Provides a convenient API for accessing data sets within ImmuneSpace (www.immunespace.org), the data repository and analysis platform of the Human Immunology Project Consortium (HIPC).

InteractionSet - Provides the GInteractions, InteractionSet and ContactMatrix objects and associated methods for storing and manipulating genomic interaction data from Hi-C and ChIA-PET experiments.

ISoLDE - This package provides ISoLDE a new method for identifying imprinted genes. This method is dedicated to data arising from RNA sequencing technologies. The ISoLDE package implements original statistical methodology described in the publication below.

isomiRs - Characterization of miRNAs and isomiRs, clustering and differential expression.

JunctionSeq - A Utility for Detection and Visualization of Differential Exon or Splice-Junction Usage in RNA-Seq data.

kimod - This package allows to work with mixed omics data (transcriptomics, proteomics, microarray-chips, rna-seq data), introducing the following improvements: distance options (for numeric and/or categorical variables) for each of the tables, bootstrap resampling techniques on the residuals matrices for all methods, that enable perform confidence ellipses for the projection of individuals, variables and biplot methodology to project variables (gene expression) on the compromise. Since the main purpose of the package is to use these techniques to omic data analysis, it includes an example data from four different microarray platforms (i.e.,Agilent, Affymetrix HGU 95, Affymetrix HGU 133 and Affymetrix HGU 133plus 2.0) on the NCI-60 cell lines.NCI60_4arrays is a list containing the NCI-60 microarray data with only few hundreds of genes randomly selected in each platform to keep the size of the package small. The data are the same that the package omicade4 used to implement the co-inertia analysis. The references in packages follow the style of the APA-6th norm.

Linnorm - Linnorm is an R package for the analysis of RNA-seq, scRNA-seq, ChIP-seq count data or any large scale count data. Its main function is to normalize and transform these datasets for parametric tests. Examples of parametric tests include using limma for differential expression analysis or differential peak detection, or calculating Pearson correlation coefficient for gene correlation study. Linnorm can work with raw count, CPM, RPKM, FPKM and TPM. Additionally, Linnorm provides the RnaXSim function for the simulation of RNA-seq raw counts for the evaluation of differential expression analysis methods. RnaXSim can simulate RNA-seq dataset in Gamma, Log Normal, Negative Binomial or Poisson distributions.

lpsymphony - This package was derived from Rsymphony_0.1-17 from CRAN. These packages provide an R interface to SYMPHONY, an open-source linear programming solver written in C++. The main difference between this package and Rsymphony is that it includes the solver source code (SYMPHONY version 5.6), while Rsymphony expects to find header and library files on the users' system. Thus the intention of lpsymphony is to provide an easy to install interface to SYMPHONY. For Windows, precompiled DLLs are included in this package.

LymphoSeq - This R package analyzes high-throughput sequencing of T and B cell receptor complementarity determining region 3 (CDR3) sequences generated by Adaptive Biotechnologies' ImmunoSEQ assay.  Its input comes from tab-separated value (.tsv) files exported from the ImmunoSEQ analyzer.

MBttest - MBttest method was developed from beta t-test method of Baggerly et al(2003). Compared to baySeq (Hard castle and Kelly 2010), DESeq (Anders and Huber 2010) and exact test (Robinson and Smyth 2007, 2008) and the GLM of McCarthy et al(2012), MBttest is of high work efficiency,that is, it has high power, high conservativeness of FDR estimation and high stability. MBttest is suit- able to transcriptomic data, tag data, SAGE data (count data) from small samples or a few replicate libraries. It can be used to identify genes, mRNA isoforms or tags differentially expressed between two conditions.

Mergeomics - The Mergeomics pipeline serves as a flexible framework for integrating multidimensional omics-disease associations, functional genomics, canonical pathways and gene-gene interaction networks to generate mechanistic hypotheses. It includes two main parts, 1) Marker set enrichment analysis (MSEA); 2) Weighted Key Driver Analysis (wKDA).

metaCCA - metaCCA performs multivariate analysis of a single or multiple GWAS based on univariate regression coefficients. It allows multivariate representation of both phenotype and genotype. metaCCA extends the statistical technique of canonical correlation analysis to the setting where original individual-level records are not available, and employs a covariance shrinkage algorithm to achieve robustness.

MethPed - Classification of pediatric tumors into biologically defined subtypes is challenging and multifaceted approaches are needed. For this aim, we developed a diagnostic classifier based on DNA methylation profiles. We offer MethPed as an easy-to-use toolbox that allows researchers and clinical diagnosticians to test single samples as well as large cohorts for subclass prediction of pediatric brain tumors.  The current version of MethPed can classify the following tumor diagnoses/subgroups: Diffuse Intrinsic Pontine Glioma (DIPG), Ependymoma, Embryonal tumors with multilayered rosettes (ETMR), Glioblastoma (GBM), Medulloblastoma (MB) - Group 3 (MB_Gr3), Group 4 (MB_Gr3), Group WNT (MB_WNT), Group SHH (MB_SHH) and Pilocytic Astrocytoma (PiloAstro).

miRNAmeConverter - Package containing an S4 class for translating mature miRNA names to different miRBase versions, checking names for validity and detecting miRBase version of a given set of names (data from http://www.mirbase.org/).

MMDiff2 - This package detects statistically significant differences between read enrichment profiles in different ChIP-Seq samples. To take advantage of shape differences it uses Kernel methods (Maximum Mean Discrepancy, MMD).

multiClust - Whole transcriptomic profiles are useful for studying the expression levels of thousands of genes across samples. Clustering algorithms are used to identify patterns in these profiles to determine clinically relevant subgroups. Feature selection is a critical integral part of the process. Currently, there are many feature selection and clustering methods to identify the relevant genes and perform clustering of samples. However, choosing the appropriate methods is difficult as recent work demonstrates that no method is the clear winner. Hence, we present an R-package called `multiClust` that allows researchers to experiment with the choice of combination of methods for gene selection and clustering with ease. In addition, using multiClust, we present the merit of gene selection and clustering methods in the context of clinical relevance of clustering, specifically clinical outcome. Our integrative R- package contains: 1. A function to read in gene expression data and format appropriately for analysis in R. 2. Four different ways to select the number of genes a. Fixed b. Percent c. Poly d. GMM 3. Four gene ranking options that order genes based on different statistical criteria a. CV_Rank b. CV_Guided c. SD_Rank d. Poly 4. Two ways to determine the cluster number a. Fixed b. Gap Statistic 5. Two clustering algorithms a. Hierarchical clustering b. K-means clustering 6. A function to calculate average gene expression in each sample cluster 7. A function to correlate sample clusters with clinical outcome Order of Function use: 1. input_file, a function to read-in the gene expression file and assign gene probe names as the rownames. 2. number_probes, a function to determine the number of probes to select for in the gene feature selection process. 3. probe_ranking, a function to select for gene probes using one of the available gene probe ranking options. 4. number_clusters, a function to determine the number of clusters to be used to cluster genes and samples. 5. cluster_analysis, a function to perform Kmeans or Hierarchical clustering analysis of the selected gene expression data. 6. avg_probe_exp, a function to produce a matrix containing the average expression of each gene probe within each sample cluster. 7. surv_analysis, a function to produce Kaplan-Meier Survival Plots of selected gene expression data.

MultiDataSet - Implementation of the BRGE's (Bioinformatic Research Group in Epidemiology from Center for Research in Environmental Epidemiology) MultiDataSet and MethylationSet. MultiDataSet is designed for integrating multi omics data sets and MethylationSet to contain normalized methylation data. These package contains base classes for MEAL and rexposome packages.

normalize450K - Precise measurements are important for epigenome-wide studies investigating DNA methylation in whole blood samples, where effect sizes are expected to be small in magnitude. The 450K platform is often affected by batch effects and proper preprocessing is recommended. This package provides functions to read and normalize 450K '.idat' files. The normalization corrects for dye bias and biases related to signal intensity and methylation of probes using local regression. No adjustment for probe type bias is performed to avoid the trade-off of precision for accuracy of beta-values.

nucleoSim - This package can generate a synthetic map with reads covering the nucleosome regions as well as a synthetic map with forward and reverse reads emulating next-generation sequencing. The user has choice between three different distributions for the read positioning: Normal, Student and Uniform.

odseq - Performs outlier detection of sequences in a multiple sequence alignment using bootstrap of predefined distance metrics. Outlier sequences can make downstream analyses unreliable or make the alignments less accurate while they are being constructed. This package implements the OD-seq algorithm proposed by Jehl et al (doi 10.1186/s12859-015-0702-1) for aligned sequences and a variant using string kernels for unaligned sequences.

OncoScore - OncoScore is a tool to measure the association of genes to cancer based on citation frequency in biomedical literature. The score is evaluated from PubMed literature by dynamically updatable web queries.

oppar - The R implementation of mCOPA package published by Wang et al. (2012). Oppar provides methods for Cancer Outlier profile Analysis. Although initially developed to detect outlier genes in cancer studies, methods presented in oppar can be used for outlier profile analysis in general. In addition, tools are provided for gene set enrichment and pathway analysis.

PanVizGenerator - PanViz is a JavaScript based visualisation tool for functionaly annotated pangenomes. PanVizGenerator is a companion for PanViz that facilitates the necessary data preprocessing step necessary to create a working PanViz visualization. The output is fully self-contained so the recipient of the visualization does not need R or PanVizGenerator installed.

pbcmc - The pbcmc package characterizes uncertainty assessment on gene expression classifiers, a. k. a. molecular signatures, based on a permutation test. In order to achieve this goal, synthetic simulated subjects are obtained by permutations of gene labels. Then, each synthetic subject is tested against the corresponding subtype classifier to build the null distribution. Thus, classification confidence measurement can be provided for each subject, to assist physician therapy choice. At present, it is only available for PAM50 implementation in genefu package but it can easily be extend to other molecular signatures.

pcaExplorer - This package provides functionality for interactive visualization of RNA-seq datasets based on Principal Components Analysis. The methods provided allow for quick information extraction and effective data exploration. A Shiny application encapsulates the whole analysis.

PCAN - Phenotypes comparison based on a pathway consensus approach. Assess the relationship between candidate genes and a set of phenotypes based on additional genes related to the candidate (e.g. Pathways or network neighbors).

pqsfinder - The main functionality of the this package is to detect DNA sequence patterns that are likely to fold into an intramolecular G-quadruplex (G4). Unlike many other approaches, this package is able to detect sequences responsible for G4s folded from imperfect G-runs containing bulges or mismatches and as such is more sensitive than competing algorithms.

profileScoreDist - Regularization and score distributions for position count matrices.

psygenet2r - Package to retrieve data from PsyGeNET database (www.psygenet.org) and to perform comorbidity studies with PsyGeNET's and user's data.

PureCN - This package estimates tumor purity, copy number, loss of heterozygosity (LOH), and status of short nucleotide variants (SNVs). PureCN is designed for hybrid capture next generation sequencing (NGS) data, integrates well with standard somatic variant detection pipelines, and has support for tumor samples without matching normal samples.

QuaternaryProd - QuaternaryProd is an R package that performs causal reasoning on biological networks, including publicly available networks such as String-db. QuaternaryProd is a free alternative to commercial products such as Quiagen and Inginuity pathway analysis. For a given a set of differentially expressed genes, QuaternaryProd computes the significance of upstream regulators in the network by performing causal reasoning using the Quaternary Dot Product Scoring Statistic (Quaternary Statistic), Ternary Dot product Scoring Statistic (Ternary Statistic) and Fisher's exact test. The Quaternary Statistic handles signed, unsigned and ambiguous edges in the network. Ambiguity arises when the direction of causality is unknown, or when the source node (e.g., a protein) has edges with conflicting signs for the same target gene. On the other hand, the Ternary Statistic provides causal reasoning using the signed and unambiguous edges only. The Vignette provides more details on the Quaternary Statistic and illustrates an example of how to perform causal reasoning using String-db.

QUBIC - The core function of this R package is to provide the implementation of the well-cited and well-reviewed QUBIC algorithm, aiming to deliver an effective and efficient biclustering capability. This package also includes the following related functions: (i) a qualitative representation of the input gene expression data, through a well-designed discretization way considering the underlying data property, which can be directly used in other biclustering programs; (ii) visualization of identified biclusters using heatmap in support of overall expression pattern analysis; (iii) bicluster-based co-expression network elucidation and visualization, where different correlation coefficient scores between a pair of genes are provided; and (iv) a generalize output format of biclusters and corresponding network can be freely downloaded so that a user can easily do following comprehensive functional enrichment analysis (e.g. DAVID) and advanced network visualization (e.g. Cytoscape).

R4RNA - A package for RNA basepair analysis, including the visualization of basepairs as arc diagrams for easy comparison and annotation of sequence and structure.  Arc diagrams can additionally be projected onto multiple sequence alignments to assess basepair conservation and covariation, with numerical methods for computing statistics for each.

recoup - recoup calculates and plots signal profiles created from short sequence reads derived from Next Generation Sequencing technologies. The profiles provided are either sumarized curve profiles or heatmap profiles. Currently, recoup supports genomic profile plots for reads derived from ChIP-Seq and RNA-Seq experiments. The package uses ggplot2 and ComplexHeatmap graphics facilities for curve and heatmap coverage profiles respectively.

RGraph2js - Generator of web pages which display interactive network/graph visualizations with D3js, jQuery and Raphael.

RImmPort - The RImmPort package simplifies access to ImmPort data for analysis in the R environment. It provides a standards-based interface to the ImmPort study data that is in a proprietary format.

ROTS - Calculates the Reproducibility-Optimized Test Statistic (ROTS) for differential testing in omics data.

SC3 - Interactive tool for clustering and analysis of single cell RNA-Seq data.

scater - A collection of tools for doing various analyses of single-cell RNA-seq gene expression data, with a focus on quality control.

scde - The scde package implements a set of statistical methods for analyzing single-cell RNA-seq data. scde fits individual error models for single-cell RNA-seq measurements. These models can then be used for assessment of differential expression between groups of cells, as well as other types of analysis. The scde package also contains the pagoda framework which applies pathway and gene set overdispersion analysis to identify and characterize putative cell subpopulations based on transcriptional signatures. The overall approach to the differential expression analysis is detailed in the following publication: "Bayesian approach to single-cell differential expression analysis" (Kharchenko PV, Silberstein L, Scadden DT, Nature Methods, doi: 10.1038/nmeth.2967). The overall approach to subpopulation identification and characterization is detailed in the following pre-print: "Characterizing transcriptional heterogeneity through pathway and gene set overdispersion analysis" (Fan J, Salathia N, Liu R, Kaeser G, Yung Y, Herman J, Kaper F, Fan JB, Zhang K, Chun J, and Kharchenko PV, Nature Methods, doi:10.1038/nmeth.3734).

scran - This package implements a variety of low-level analyses of single-cell RNA-seq data. Methods are provided for normalization of cell-specific biases, assignment of cell cycle phase, and detection of highly variable and significantly correlated genes.

SMITE - This package builds on the Epimods framework which facilitates finding weighted subnetworks ("modules") on Illumina Infinium 27k arrays using the SpinGlass algorithm, as implemented in the iGraph package. We have created a class of gene centric annotations associated with p-values and effect sizes and scores from any researchers prior statistical results to find functional modules.

SpidermiR - The aims of SpidermiR are : i) facilitate the network open-access data retrieval from GeneMania data, ii) prepare the data using the appropriate gene nomenclature, iii) integration of miRNA data in a specific network, iv) provide different standard analyses and v) allow the user to visualize the results. In more detail, the package provides multiple methods for query, prepare and download network data (GeneMania), and the integration with  validated and predicted miRNA data (mirWalk, miR2Disease,miRTar, miRandola,Pharmaco-miR,DIANA, Miranda, PicTar and TargetScan) and the use of standard analysis (igraph) and visualization methods (networkD3).

splineTimeR - This package provides functions for differential gene expression analysis of gene expression time-course data. Natural cubic spline regression models are used. Identified genes may further be used for pathway enrichment analysis and/or the reconstruction of time dependent gene regulatory association networks.

sscu - The package can calculate the selection in codon usage in bacteria species. First and most important, the package can calculate the strength of selected codon usage bias (sscu) based on Paul Sharp's method. The method take into account of background mutation rate, and focus only on codons with universal translational advantages in all bacterial species. Thus the sscu index is comparable among different species. In addition, detainled optimal codons (selected codons) information can be calculated by optimal_codons function, so the users will have a more accurate selective scheme for each codons. Furthermore, we added one more function optimal_index in the package. The function has similar mathematical formula as s index, but focus on the estimates the amount of GC-ending optimal codon for the highly expressed genes in the four and six codon boxes. The function takes into account of background mutation rate, and it is comparable with the s index. However, since the set of GC-ending optimal codons are likely to be different among different species, the index can not be compared among different species.

SwathXtend - It contains utility functions for integrating spectral libraries for SWATH and statistical data analysis for SWATH generated data.

tofsims - This packages offers a pipeline for import, processing and analysis of ToF-SIMS 2D image data. Import of Iontof and Ulvac-Phi raw or preprocessed data is supported. For rawdata, mass calibration, peak picking and peak integration exist. General funcionality includes data binning, scaling, image subsetting and visualization. A range of multivariate tools common in the ToF-SIMS community are implemented (PCA, MCR, MAF, MNF). An interface to the bioconductor image processing package EBImage offers image segmentation functionality.

transcriptR - The differences in the RNA types being sequenced have an impact on the resulting sequencing profiles. mRNA-seq data is enriched with reads derived from exons, while GRO-, nucRNA- and chrRNA-seq demonstrate a substantial broader coverage of both exonic and intronic regions. The presence of intronic reads in GRO-seq type of data makes it possible to use it to computationally identify and quantify all de novo continuous regions of transcription distributed across the genome. This type of data, however, is more challenging to interpret and less common practice compared to mRNA-seq. One of the challenges for primary transcript detection concerns the simultaneous transcription of closely spaced genes, which needs to be properly divided into individually transcribed units. The R package transcriptR combines RNA-seq data with ChIP-seq data of histone modifications that mark active Transcription Start Sites (TSSs), such as, H3K4me3 or H3K9/14Ac to overcome this challenge. The advantage of this approach over the use of, for example, gene annotations is that this approach is data driven and therefore able to deal also with novel and case specific events. Furthermore, the integration of ChIP- and RNA-seq data allows the identification all known and novel active transcription start sites within a given sample.

tximport - Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. Average transcript length, weighted by sample-specific transcript abundance estimates, is provided as a matrix which can be used as an offset for different expression of gene-level counts.

Uniquorn - This packages enables users to identify cancer cell lines. Cancer cell line misidentification and cross-contamination reprents a significant challenge for cancer researchers. The identification is vital and in the frame of this package based on the locations/ loci of somatic and germline mutations/ variations. The input format is vcf/ vcf.gz and the files have to contain a single cancer cell line sample (i.e. a single member/genotype/gt column in the vcf file). The implemented method is optimized for the Next-generation whole exome and whole genome DNA-sequencing technology.

NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages since the last release. The following package NEWS is available:

ADaCGH2
-------

Changes in version 2.11.2 (2016-04-29):

- GLAD is now in Depends, not Imports, as it gives warnigns and notes.

Changes in version 2.11.1 (2016-03-29):

- Vignette changes to avoid Windoze build issues?

- Changes in NAMESPACE

affxparser
----------

Changes in version 1.43.2:

- Version bump for Bioconductor

Changes in version 1.43.1-9000 (2016-04-05):

- WINDOWS: Package now compiles with both the old gcc-4.6.3 toolchain
  as well as the new gcc-4.9.3 toolchain - introduced in R (>= 3.3.0).
  Thanks to Jim Hester and Dan Tenenbaum for help with this.

Changes in version 1.43.1 (2016-02-28):

- Fixed a bug related to including <R.h> and extern C, reported by
  Brian Ripley.

Changes in version 1.43.0 (2015-10-23):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.3 for R (>= 3.3.0).

AllelicImbalance
----------------

Changes in version 1.10.0:

NEW FEATURES

- minFreqFilt and minCountFilt

SIGNIFICANT USER-VISIBLE CHANGES

- The Phase2Genotype function now treats phase 0 as the ref allele and
  phase 1 as the alternative allele.

AnnotationDbi
-------------

Changes in version 1.34.0:

NEW FEATURES

- export ALIAS2EG symbol in NAMESPACE for frog, mosquito, chimp and
  rhesus OrgDbs

- add how to use EnsDb section to vignette

MODIFICATIONS

- work on code base and exported functions - break up code in
  geneCentricDbs file - re-name and export functions used by select()
  in AnnotationForge and GenomicFeatures - re-name and export functions
  used by annotation package builders in AnnotationForge

- remove library(RSQLite) from dbFileConnect()

- modify unit tests for new PFAM identifiers in Bioconductor 3.3

- use elementNROWS() instead of elementLengths()

- modify mapIds() to preserve data type returned from select()

- reduce time of AnnDbPkg-checker.Rd example

BUG FIXES

- bugfix for select() and mapIds() when there are many:many mappings

- bugfix in test_generateExtraRows unit test.

AnnotationForge
---------------

Changes in version 1.14.0:

NEW FEATURES

- build 'alias' table in OrgDbs sqlite db for frog, chimp, rhesus and
  mosquito

MODIFICATIONS

- update AnnDbPkg template man to reference select() interface

- work on makeOrgPackageFromNCBI: - error early when tax id is not
  found in gene_info.gz - add 'rebuildCache' arg to control repeated
  downloads - remove old code and comments - update man page

- add PFAM and PROSITE man pages to NCBICHIP and NCBIORG package
  templates

- allow passing of directory location in wrapBaseDBPackages()

- change format of licence; report current version of AnnotationDbi

- modify appendArabidopsisGenes() to check for null 'gene_id'

- add DBI to 'Suggests'; load DBI in `_dbconn` man page

- load SQLite in vignettes; no longer free from
  AnnotationDbi::dbFileConnect()

AnnotationHub
-------------

Changes in version 2.4.0:

NEW FEATURES

- add new status codes '4' and '5' to 'statuses' mysql table; change
  'status_id' field to '4' for all removed records to date

- add getRecordStatus() generic

- add package() generic

- create 'Hub' VIRTUAL class - add new .Hub() base constructor for all
  hubs - add getAnnotationHubOption() and setAnnotationHubOption() -
  promote cache() to generic - add getHub() getter for
  AnnotationHubResource class - add getUrl(), getCache(), getDate()
  getters - export as few db helpers as possible

- add 'EpigenomeRoadmapNarrowAllPeaks' and 'EpigenomeRoadmapNarrowFDR'
  classes

MODIFICATIONS

- distinguish between broad and narrow peak files in
  EpigenomeRoadmapFileResource dispatch class

- don't use cache for AnnotationHub SQLite connection - originally
  introduced so could be closed if needed, but creates complexity -
  instead, open / close connection around individual queries (not a
  performance concern) - expose hub, cache, proxy in AnnotationHub
  constructor - document dbconn,Hub-method, dbfile,Hub-method,
  .db_close

- snapshotDate now uses timestamp (last date any row was modified)
  instead of rdatadateadded

- .require fails rather than emits warning - unit test on .require() -
  also, cache(hub[FALSE]) does not create spurious error

- work on removed records and biocVersion - .uid0() was reorganized and
  no longer groups by record_id - metadata is returned for records with
  biocversion field <= current biocVersion instead of an exact match
  with the current version - metadata is not returned for removed
  records

BUG FIXES

- Work around httr() progress() bug by disabling progress bar

AnnotationHubData
-----------------

Changes in version 1.2.0:

NEW FEATURES

- add makeEnsemblTwoBit()

- add hubError(), hubError<- generics and methods

- create 'HubMetadata' class which 'AnnotationHubMetadata' inherits
  from

MODIFICATIONS

- export ensemblFastaToTwoBitFile()

- modifications due to changes in httr::HEAD(): - AFAICT httr::HEAD()
  1.1.0 and higher accepts https only, not ftp - use xml2 instead of XML for
  parsing (httr >= 1.1.0 dependency change)

- work on recipes: - clean up ChEA and Gencode - don't export
  tracksToUpdate(); was broken and not used - reorg man pages; combine
  Ensembl Fasta and TwoBit on single man page

- work on updateResources(): - push data to S3 before inserting
  metadata in db - isolate pushResources() and pushMetadata() from
  updateResources() - NOTE: Epigenome unit test is failing due to bad
  url. If not fixed by the host the recipe will need to change.

- update makedbSNPVCF() to look in new clinvar location

BUG FIXES

- fix bugs in makedbSNPVCF() recipe related to genome and tags

antiProfiles
------------

Changes in version 1.11.1:

- Use subsetted computation from matrixStats package.

aroma.light
-----------

Changes in version 3.1.1 (2016-01-06):

- Package requires R (>= 2.15.2).

- CLEANUP: robustSmoothSpline() no longer generates messages that
  ".nknots.smspl() is now exported; use it instead of n.knots()" for R
  (>= 3.1.1).

Changes in version 3.1.0 (2015-10-23):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.3 for R (>= 3.3.0).

bamsignals
----------

Changes in version 1.3:

- Bioconductor Release (3.3) Version NEW FEATURES

- Fragment Length Filter: Filtering of paired end read pairs by TLEN
  field with a minimum and maximum TLEN

- SAMFLAG Filtering: Filter out reads with a certain SAMFLAG set, e.g.
  1024 for marked optical duplicates BUG FIXES

- Path extension for relative paths, e.g. resolve "~" to /home/$USER in
  UNIX

BasicSTARRseq
-------------

Changes in version 1.0.0:

- Provided function getPeaks for detecting significant peaks in
  STARR-seq data with two different binomial models.

BatchQC
-------

Changes in version 1.0.0:

- Summary and Sample Diagnostics

- Differential Expression Plots and Analysis using LIMMA

- Principal Component Analysis and plots to check batch effects

- Heatmap plot of gene expressions

- Median Correlation Plot

- Circular Dendrogram clustered and colored by batch and condition

- Shape Analysis for the distribution curve based on HTShape package

- Batch Adjustment using ComBat

- Surrogate Variable Analysis using sva package

- Function to generate simulated RNA-Seq data

BgeeDB
------

Changes in version 0.99.6 (2016-04-28):

- Bug fix in format_data function.

Changes in version 0.99.4 (2016-04-26):

- Bug fix for Affymetrix data.

Changes in version 0.99.3 (2016-04-22):

- Added download for Affymetrix data.

Changes in version 0.0.4 (2016-03-20):

- Added loadTopAnatData function.

Changes in version 0.0.1 (2016-03-13):

- Added topAnat script.  This function produces a topAnatObject, ready
  to use for gene set enrichment testing using functions from the topGO
  package

bigmemoryExtras
---------------

Changes in version 1.17.0:

- All Biobase-related content has been removed including
  annotatedDataFrameFrom. This means that BigMatrix objects can no
  longer be added to eSet/ExpressionSet/GenoSet objects. Try the new
  SummarizedExperiment, which is the eSet replacement.

Biobase
-------

Changes in version 2.31:

NEW FEATURES

- tab completion implemented for eSet classes

- head and tail.AnnotatedDataFrame methods introduced

BiocStyle
---------

Changes in version 2.0.0:

- New Bioconductor LaTeX Style. See package vignettes for details.

biomformat
----------

Changes in version 0.3.13:

USER-VISIBLE CHANGES

- Added make_biom function. Creates biom object from standard R data
  table(s).

Changes in version 0.3.12:

USER-VISIBLE CHANGES

- No user-visible changes. All future compatibility changes.

BUG FIXES

- Unit test changes to work with upcoming R release and new testthat
  version.

- This solves Issue 4: https://github.com/joey711/biom/issues/4

Changes in version 0.3.11:

USER-VISIBLE CHANGES

- No user-visible changes. All future CRAN compatibility changes.

BUG FIXES

- Clarified license and project in the README.md

- Added TODO.html, README.html, and TODO.md to .Rbuildignore (requested
  by CRAN)

- Moved `biom-demo.Rmd` to `vignettes/`

- Updated `inst/NEWS` (this) file to official format

- Removed pre-built vignette HTML so that it is re-built during package
  build.  This updates things like the build-date in the vignette, but
  also ensures that the user sees in the vignette the results of code
  that just worked with their copy of the package.

Changes in version 0.3.10:

USER-VISIBLE CHANGES

- These changes should not affect any package behavior.

- Some of the top-level documentation has been changed to reflect new
  development location on GitHub.

BUG FIXES

- Minor fixes for CRAN compatibility

- This addresses Issue 1: https://github.com/joey711/biom/issues/1

Changes in version 0.3.9:

SIGNIFICANT USER-VISIBLE CHANGES

- speed improvement for sparse matrices

NEW FEATURES

- The `biom_data` parsing function now uses a vectorized
  (matrix-indexed) assignment while parsing sparse matrices.

- Unofficial benchmarks estimate a few 100X speedup.

Changes in version 0.3.8:

SIGNIFICANT USER-VISIBLE CHANGES

- First release version released on CRAN

BioQC
-----

Changes in version 1.1-5:

- Add C-level implementation of Wilcoxon-Mann-Whitney rank sum test

- Documentation and vignettes updated to be ready for Bioconductor
  submission.


biosigner
---------

Changes in version 0.99.12:

PACKAGE MODIFICATION

- Welcome to the biosigner package for feature selection from omics
  datasets

- The package implements a new wrapper method detecting the features
  which are important for PLS-DA, Random Forest, or SVM binary
  classification

- The package contains the 'diaplasma' LC-MS metabolomics real dataset
  (plasma samples from diabetic type 1 and 2 patients)

- Please see the vignette for details about the approach and package
  use

- The corresponding publication is currently under review.

Changes in version 0.99.11:

BUG FIXED

- in vignette (due to switch in S4 methods for ropls)

Changes in version 0.99.10:

PACKAGE MODIFICATION

- grDevices, graphics, stats, utils imported in DESCRIPTION

Changes in version 0.99.9:

PACKAGE MODIFICATION

- vignette update

Changes in version 0.99.8:

PACKAGE MODIFICATION

- adding the import of the following function in NAMESPACE: abline
  arrows axis box boxplot dev.new dev.off head image layout median
  mtext par pdf rect tail title var

- defining the getAccuracyMN and getSignatureLs accessors

Changes in version 0.99.7:

PACKAGE MODIFICATION

- additional unit test silenced in test_biosign_randomforest (because
  of errors on the moscato2 Windows 2008 platform)

Changes in version 0.99.6:

PACKAGE MODIFICATION

- unit tests silenced in test_biosign_randomforest and
  test_biosign_predict (because of errors on the moscato2 Windows 2008
  platform)

Changes in version 0.99.5:

PACKAGE MODIFICATION

- importing of packages in NAMESPACE fixed

- use of S4 methods (instead of S3)

Changes in version 0.99.4:

PACKAGE MODIFICATION

- correction of a bug in the test_biosign_plsdaF test function

Changes in version 0.99.0:

PACKAGE MODIFICATION

- set version number to 0.99.0 for Bioconductor


biovizBase
----------

Changes in version 1.19.1:

NEW FEATURES

- crunch method supports now EnsDb objects and related filter objects.

BiRewire
--------

Changes in version 3.1:

NEW FEATURES

- birewire.rewire.dsg has been introduced in order to rewire directed
  signed network (dsg). Moreover the sampler routine, analysis, visual
  monitoring and similarity has been written in order to have the same
  functionalities the package have for bipartite and undirected graphs
  also for the dsg.

BUG FIXES

- birewire.sampler.bipartite creates more files than required in the
  case of non-sparse writing procedure: fixed in 3.1

BrowserViz
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

bsseq
-----

Changes in version 1.7:

- Fixing an error with reading multiple samples aligned with bismark in
  the format "cytosineReport".

CAMERA
------

Changes in version 1.27.2:

BUG FIXES

- Fix imports in DESCRIPTION to avoid igraph overwriting xcms::group()

Cardinal
--------

Changes in version 1.3.3:

BUG FIXES

- Subsetting the S4 part of Binmat objects by row is now an error

- Providing non-positive m/z values to 'readImzML' is now an error

- Elements of 'imageData' that fail to 'combine' or which are missing
  from one or more of the objects are now dropped from the result with
  warning rather than failing

- Moved unit tests in 'ints/tests' to 'tests/testthat'

Changes in version 1.3.2:

NEW FEATURES

- Added 'image3D' method for plotting 3D images

- Added 'batchProcess' method for batch pre-processing

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Added 'mass.accuracy' and 'units.accuracy' arguments for controlling
  the m/z accuracy when reading 'processed' imzML

- Function 'reduceDimension.bin' now takes argument 'units' with value
  'ppm' or 'mz', and new corresponding defaults

BUG FIXES

- Fixed bug in reading 'processed' imzML format that caused mass
  spectra to be reconstructed in the wrong order

- Improved speed accessing columns of Hashmat sparse matrices

- In 'pixelApply' and 'featureApply', zero-length return values are no
  longer returned as a list when '.simplify=FALSE'

- Function 'peakAlign.diff' should be more memory efficient now

Changes in version 1.3.0 (2015-12-16):

NEW FEATURES

- Added experimental Binmat class for working with on-disk matrices

- Added experimental support for 3D files from benchmark datasets

- Added experimental support for plotting 3D images

- Added experimental support for 'processed' imzML format

SIGNIFICANT USER-VISIBLE CHANGES

- Added 'attach.only' argument to readImzML and readAnalyze

BUG FIXES

- Fixed bug in plotting 3D image slices in the z dimension

- Fixed bug where large imzML files could not be read due to byte
  offsets being stored as ints; they are now stored as doubles.

- Fixed bug with strip labels in 3D plotting and with mixed labels

- Fixed bug with unique m/z feature names for high mass resolutions

cellity
-------

Changes in version 0.99.2 (2016-02-22):

- Package prepared for Bioconductor submission.

ChAMP
-----

Changes in version 1.9.8:

- Make ChAMP package suitable for both 450K array and EPIC array.

- Updated BMIQ normalization to newest version.

- Add champ.refbase function to do reference-based cell proportion
  detection and correction.

- Add champ reffree function to do reference-free EWAS.

- removed champ.lasso function, but added champ.DMR function, which
  combined bumphunter algorithm and probe lasso algorithm to detect
  DMRs.

- Updated a new version of vignette.

ChIPpeakAnno
------------

Changes in version 3.5.18:

- Fix the bugs for featureAlignedSignal when features contain strand
  info.

Changes in version 3.5.17:

- Fix the bugs for getEnrichedGO/PATH when recodes in database such as
  reactome database is not unique.

Changes in version 3.5.16:

- fix the warning in windows for replacing previous import

Changes in version 3.5.15:

- not do log2 transform for featureAlignedSingal

Changes in version 3.5.14:

- remove NA for featureAlignedDistribution

- handle NA and infinite value for featureAlignedSingal

Changes in version 3.5.13:

- throw message when seqlevels are not identical for
  featureAlignedSignal

Changes in version 3.5.12:

- add new function summarizeOverlapsByBins.

Changes in version 3.5.11:

- clean output of findEnhancers

Changes in version 3.5.10:

- add new function findEnhancers

- normalize the output of gene region for binOverFeature

- change the documentation to avoid time-out error

Changes in version 3.5.9:

- merge peaksNearBDP and bdp function.

- improve oligoSummary.

- update the documentations.

Changes in version 3.5.8:

- change toGRanges from function to method

- update the documentations.

Changes in version 3.5.7:

- change test from RUnit to testthat

- add new function addMetadata

- change the output and parameters of annoPeaks

- simple the parameter output of annotatePeakInBatch

- allow bdp function to accept GRanges annotation

- add error bar function for binOverFeature function

- remove the log file after plot for makeVennDiagram function

- add private function trimPeakList

- update the documentation of annoPeaks and annotatePeakInBatch

Changes in version 3.5.5:

- update the documentation to fix the typo in quickStart.

- change the default value of annoPeaks.

- update annoGR class to fix the error: identical(colnames(classinfo),
  colnames(out)) is not TRUE

Changes in version 3.5.2:

- update the documentation to fix the error on windows (import bigwig
  error)

- avoid the output of addGeneIDs as factors

Changes in version 3.5.1:

- update the documentation NEW FEATURE

- toGRanges can accept connection object

- add annoPeaks function

- add xget function

- update the peakPermTest algorithm to make more reasonable result.

- add oligoSummary function

- add IDRfilter function

- add reCenterPeaks function

ChIPQC
------

Changes in version 1.7.2:

- Convert data objects for DiffBind compatibility

Changes in version 1.7.1:

- Maintain compatibility with changes in DiffBind

ChIPseeker
----------

Changes in version 1.7.15:

- update GEO data <2016-03-21, Mon>

Changes in version 1.7.14:

- list_to_dataframe works with data frames that have different colnames
  <2016-03-10, Thu>

Changes in version 1.7.13:

- support annotate peaks with custom regions via passing
  TxDb=user_defined_GRanges to annotatePeak <2016-03-06, Sun>

Changes in version 1.7.12:

- fixed R check <2016-03-05, Sat>

- implement list_to_dataframe that mimic ldply and remove ldply
  dependency <2016-03-05, Sat>

Changes in version 1.7.11:

- fixed issue in testing list in covplot introduced in 1.7.9
  <2016-03-02, Wed>

Changes in version 1.7.10:

- determined gene ID type if TxDb doesn't contain corresponding
  metadata <2016-03-01, Tue> + fixed
  https://github.com/GuangchuangYu/ChIPseeker/issues/28

Changes in version 1.7.9:

- covplot support GRangesList <2016-02-24, Wed>

- update ReactomePA citation info <2016-02-17, Wed>

Changes in version 1.7.8:

- fixed BUG of Peaks upstream first or downstream last gene not
  annotated <2016-01-20, Wed> + contributed by Michael Kluge + see
  https://github.com/GuangchuangYu/ChIPseeker/pull/24

Changes in version 1.7.7:

- bug fixed in newly introduced parameter 'overlap'. solve NA issue.
  <2016-01-13, Wed>

Changes in version 1.7.6:

- introduce 'overlap' parameter in annotatePeak, by default
  overlap="TSS" and only overlap with TSS will be reported as the
  nearest gene. if overlap="all", then gene overlap with peak will be
  reported as nearest gene, no matter the overlap is at TSS region or
  not. <2016-01-12, Tue>

- bug fixed in find overlap with peaks have strand info. <2016-01-12,
  Tue> + see https://github.com/GuangchuangYu/ChIPseeker/issues/23

Changes in version 1.7.5:

- add paramters, sameStrand,ignoreOverlap, ignoreUpstream and
  ignoreDownstream in annotatePeak <2016-01-10, Sun> + see
  https://github.com/GuangchuangYu/ChIPseeker/issues/17

- bug fixed in peak orientation <2016-01-10, Sun> + see
  https://github.com/GuangchuangYu/ChIPseeker/issues/22

Changes in version 1.7.4:

- stop if input list of csAnno object has no name attribute + see
  https://github.com/GuangchuangYu/ChIPseeker/issues/21 + plotAnnoBar +
  plotDistToTSS

- [covplot] xlim now not only restrict the window of data but also set
  the limit of the graphic object <2015-12-30, Wed> + see
  https://github.com/GuangchuangYu/ChIPseeker/issues/20

Changes in version 1.7.3:

- fixed R check <2015-12-29, Tue>

Changes in version 1.7.2:

- use geom_rect instead of geom_segment in covplot <2015-11-30, Mon>

- open lower parameter (by default =1) to specific lower cutoff of
  coverage signal <2015-11-29, Sun>

- fixed covplot to work with None RleViews of specific chromosome
  <2015-11-29, Sun>

- addFlankGeneInfo now works with level="gene" <2015-11-19, Thu> + see
  https://github.com/GuangchuangYu/ChIPseeker/issues/18

Changes in version 1.7.1:

- fixed extracting ID type from TxDb object, since the change of
  metadata(TxDb). now using grep to extract. <2015-10-27, Tue>

- add vp parameter to set viewport of vennpie on top of upsetplot by
  user request <2015-10-26, Mon> + see
  http://ygc.name/2015/07/28/upsetplot-in-chipseeker/#comment-19470

- getBioRegion function <2015-10-20, Tue> + see
  https://github.com/GuangchuangYu/ChIPseeker/issues/16

Changes in version 1.7.0:

- BioC 3.3 branch

clonotypeR
----------

Changes in version 1.10.0:

OTHER CHANGES

- Minor metadata updates.

clusterProfiler
---------------

Changes in version 2.99.2:

- add keyType parameter in enrichKEGG, enrichMKEGG, gseKEGG and
  gseMKEGG <2016-05-03, Tue>

- search_kegg_species function <2016-05-03, Tue>

- bitr_kegg function <2016-05-03, Tue>

Changes in version 2.99.1:

- go2ont function <2016-04-28, Thu> +
  https://www.biostars.org/p/188863/

- go2term function <2016-04-21, Thu>

- export buildGOmap <2016-04-08, Fri>

- fixed enrichDAVID according to
  https://support.bioconductor.org/p/72188/ <2016-04-08, Fri>

Changes in version 2.99.0:

- version bump, ready to release version 3. <2016-03-29, Tue> + KEGG
  pathway and module support all species listed in KEGG website. + GO
  supports all species that have OrgDb which can be accessed online by
  AnnotationHub + User's defined/customized annotation is supported via
  enricher & GSEA

Changes in version 2.5.6:

- write my own code to parse KEGG REST instead of using KEGGREST
  package <2016-03-10, Thu> + due to the issue of KEGGREST can't work
  with network of HTTP proxy with password + see:
  http://guangchuangyu.github.io/2015/02/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/#comment-2561443364

Changes in version 2.5.5:

- maxGSSize paramter <2016-03-9, Wed>

- update show method of compareClusterResult <2016-03-06, Sun>

- fixed R check <2016-03-06, Sun>

- update ReactomePA citation info <2016-02-17, Wed>

Changes in version 2.5.4:

- add use_internal_data=FALSE parameter in gseKEGG <2016-01-19, Tue>

- fixed bug in simplify for organism extracted from OrgDb is not
  consistent to GOSemSim. <2016-01-05, Tue>

- update vignette <2015-12-29, Tue>

- re-designed internal function <2015-12-20, Sun>

Changes in version 2.5.3:

- enrichMKEGG and gseMKEGG functions to support enrichment analysis of
  KEGG Module <2015-12-03, Thu> + see
  https://github.com/GuangchuangYu/clusterProfiler/issues/37

Changes in version 2.5.2:

- facet with compareCluster result <2015-11-03, Tue> + see
  https://github.com/GuangchuangYu/clusterProfiler/issues/32

Changes in version 2.5.1:

- read.gmt function for parsing gmt file for enricher and GSEA
  <2015-10-28, Wed>

- gofilter to filt result at specific GO level <2015-10-23, Fri>

- simplify function <2015-10-21, Wed> + see
  https://github.com/GuangchuangYu/clusterProfiler/issues/28

Changes in version 2.5.0:

- BioC 3.3 branch


CNEr
----

Changes in version 3.3:

BUG FIXES

- Fix a bug caused by "format" in the blat step

NEW FEATURES

- Add the pairwise whole genome alignment pipeline

- Add a new class "GRangePairs"

cogena
------

Changes in version 1.5.2 (2016-03-31):

- update vignette

Changes in version 1.5.1 (2015-10-30):

- add gmtlist2file function

Changes in version 1.5.0 (2015-10-14):

- Bioconductor 3.2 Dev

ComplexHeatmap
--------------

Changes in version 1.9.7:

- add `Legend()` function which is more flexible to generate different
  types of legends.

Changes in version 1.9.6:

- `color_mapping_legend()`, there are ticks on continuous color bar

Changes in version 1.9.5:

- add a section in the vignette to show how to adjust positions of
  column names when there are bottom annotations.

- fixed a bug that character NA values can not to assigned with na_col

- extra character 'at' and 'labels' in legends will be removed
  automatically

- all arguments which are passed to `make_layout()` are all explicitly
  put in `draw()` instead of using ...

Changes in version 1.9.4:

- heatmap bodied can be set as raster images if number of rows are too
  huge

Changes in version 1.9.3:

- graphic parameters are correctly recycled in row annotations

- if there is only one row after splitting, there will be no dendrogram

- add `range` option in `densityHeatmap()`

- when `gap` is set for the main heatmap, other heatmps also adjust
  their `gap` values to it

- fixed a bug that when rownames/colnames are not complete, dendrograms
  are corrupted

- `alter_fun` now supports adding graphics grid by grid

- add `show_pct` option in `oncoPrint()`

- add `column_order` in `densityHeatmap()`

Changes in version 1.9.2:

- imporved example of `anno_link()`

Changes in version 1.9.1:

- width of the heatmap body are calculated correctly if it is set as a
  fixed unit

- there is no dendrogram is nrows in a row-slice is 1

- add `anno_link()` annotation function

- bottom annotations are attached to the bottom edge of the heatmap if
  there are additional blank space

- colors for NA can be set by "_NA_" in annotations

- `row_dend_reorder` and `column_dend_reorder` are set to `TRUE` by
  default again

- optimize the way to specify na_col in heatmap annotations

- correct wrong viewport names in decorate_* functions

CONFESS
-------

Changes in version 1.0.0 (2016-04-28):

NEW FEATURES

- Cell detection and signal estimation model for images produced using
  the Fluidigm C1 system.

- We applied this model to a set of HeLa cell expressing fluorescence
  cell cycle reporters.

- Accompanying dataset available in the CONFESSdata package.


cosmiq
------

Changes in version 1.5.1:

USER VISIBLE CHANGES

- added Dockerfile

COSNet
------

Changes in version 1.5.1:

- Added the URL field in description file to link the corresponfing
  GitHub data USER VISIBLE CHANGES BUG FIXES PLANS


csaw
----

Changes in version 1.5.10:

- Restored normalize() as a S4 method returning a
  RangedSummarizedExperiment object.

- Modified asDGEList() to use any available normalization data in the
  input object.

- Generalized S4 methods to apply on SummarizedExperiment objects.

- Removed the rescue.ext option for PE handling, to maintain consistent
  totals calculations.

- Removed the fast.pe option for PE data handling, in favour of
  improved default processing.

- Removed dumpPE(), which is not required without the fast.pe option.

- Removed makeExtVector() in favour of list/DataFrame specification.

- windowCounts() and regionCounts() now compute and store the mean PE
  size and read length.

- Minor fix in correlateReads() for end-of-chromosome behaviour.

- Modified checkBimodality() so that the width argument behaves like
  ext in windowCounts().

- extractReads() with as.reads=TRUE for PE data now returns a
  GRangesList.

- Added the controlClusterFDR(), clusterWindows() and
  consolidateClusters() functions to automate control of the
  cluster-level FDR.

- Added protection against NA values in the cluster IDs for
  combineTests(), getBestTest(), upweightSummits().

- All read extraction methods are now CIGAR-aware and will ignore
  soft-clipped parts of the alignment.

dagLogo
-------

Changes in version 1.9.3:

- To catch the error when biomaRt server does not response correctly in
  documentations.

Changes in version 1.9.2:

- change fetchSequence to accept sequences.

Changes in version 1.9.1:

- update the test to testthat

dcGSA
-----

Changes in version 0.99.0:

- The paper on dcGSA is under review process and will update here when
  published.

DChIPRep
--------

Changes in version 1.1.8:

- summarizeCountsPerPosition has been improved and erroneously
  hard-coded options can now be actually set. This improved version was
  developed by Sebastian Gibb (http://sebastiangibb.de/).

Changes in version 1.1.7:

- added new function importData_soGGi to import .bam files in R
  directly

Changes in version 1.1.6:

- added CITATION file containing the information on the PeerJ preprint
  BUG FIXES

- added require(mgcv) to the examplse from from plotProfiles to pass
  the check() from devtools

- also added mgcv to suggests

Changes in version 1.1.5:

BUG FIXES

- bug fixes by Hervé Pages

Changes in version 1.1.4:

BUG FIXES

- bug fixes by Hervé Pages

Changes in version 1.1.3:

BUG FIXES

- fixed a wrong preview of the sample annotation table in the vignette

- fixed a wrong reference in the runTesting function help

- made small changes to the export statement in the roxygen docs

Changes in version 1.1.2:

BUG FIXES

- changed the bibliography to a standard rmarkdown one. This should fix
  the vignette compilation problem.

Changes in version 1.1.1:

BUG FIXES

- fixed a literature query in the vignette that tried to resolve a DOI
  via knitcitations::citep() and would often lead to a compilation
  error on the Bioconductor servers.

Changes in version 1.1.0:

- initial Bioc devel release


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

DEFormats
---------

Changes in version 1.0.0:

INITIAL RELEASE OF THE PACKAGE

- Converter functions between the ‘DESeqDataSet’ and ‘DGEList’ objects
  from the DESeq2 and edgeR packages, respectively.

- S4 generic for the ‘DGEList’ constructor and a corresponding method
  for ‘RangedSummarizedExperiment’ objects.

derfinder
---------

Changes in version 1.5.39:

BUG FIXES

- annotateRegions() now ignores strand by default

Changes in version 1.5.37:

SIGNIFICANT USER-VISIBLE CHANGES

- The users guide vignette now has a short section explaining how to
  use derfinder for differential binding analysis with ChIP-seq data.

Changes in version 1.5.27:

SIGNIFICANT USER-VISIBLE CHANGES

- Added smoothing arguments for the single base-level approach based on
  functions from the bumphunter package. These arguments are useful for
  identifying differentially bounded ChIP-seq peaks.

Changes in version 1.5.13:

BUG FIXES

- Fixed railMatrix()'s flexibility for defining the cluster used for
  loading the BigWig files. You can now use 'BPPARAM.railChr' which
  will take priority over 'file.cores'. Also, if 'file.cores = 1L',
  then the default will be to use SerialParam(), which was the
  implementation available prior to 1.5.11.

Changes in version 1.5.11:

SIGNIFICANT USER-VISIBLE CHANGES

- Now coverageToExon(), regionMatrix() and railMatrix() can take an 'L'
  argument of length equal to the number of samples in case not all
  samples have the same read length.

- railMatrix() has a new argument called 'file.cores' for controlling
  how many cores are used for loading the BigWig files. In theory this
  allows using railMatrix() with 'BPPARAM.custom' equal to a
  BatchJobsParam() to submit 1 job per chromosome, then 'file.cores'
  determines the number of cores for reading the files. This is a
  highly experimental feature.

Changes in version 1.5.9:

SIGNIFICANT USER-VISIBLE CHANGES

- Dropped the old introductory and advanced vignettes. We think that
  the new vignettes are clearer. In particular, they do a better job at
  highlighting the differences between the expressed regions-level and
  single base-level F-statistics implementations of the DER Finder
  approach to RNA-seq data.

Changes in version 1.5.8:

SIGNIFICANT USER-VISIBLE CHANGES

- Added a users guide vignette which explains nearly every detail you
  would want to know as a user.

Changes in version 1.5.7:

SIGNIFICANT USER-VISIBLE CHANGES

- Added a quick start vignette.

Changes in version 1.5.6:

NEW FEATURES

- Introduced railMatrix() which generates similar output to
  regionMatrix() but is much faster and less memory intensive. It
  achieves this by extracting the required information from BigWig
  files.

derfinderPlot
-------------

Changes in version 1.5.8:

BUG FIXES

- plotRegionCoverage() used to take into account the strand of the
  regions for finding transcripts that overlapped the regions. This was
  not a problem with DERs from derfinder since they have strand * by
  default but it is a problem when using it with stranded regions.

- plotCluster() will also now ignore strand for finding neighboring
  regions.

Changes in version 1.5.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Only use distance if it's not NA in plotRegionCoverage()

DESeq2
------

Changes in version 1.12.0:

- Added DESeqDataSetFromTximport() to import counts using tximport.

- Added vst() a fast wrapper for the VST.

- Added support for IHW p-value adjustment.

Changes in version 1.11.42:

- Update summary() to be IHW-results-aware.

- Small change to fitted mu values to improve fit stability when counts
  are very low. Inference for high count genes is not affected.

- Galaxy script inst/script/deseq2.R moves to Galaxy repo.

Changes in version 1.11.33:

- Changed 'filterFun' API to accommodate IHW: independent hypothesis
  weighting in results(), see vignette for example code.  Thanks to
  Nikolaos Ignatiadis, maintainer of IHW package.

Changes in version 1.11.18:

- Added a function vst(), which is a fast wrapper for
  varianceStabilizingTransformation(). The speed-up is accomplished by
  subsetting to a smaller number of genes for the estimation of the
  dispersion trend.

Changes in version 1.11.5:

- Adding in functionality to import estimated counts and average
  transcript length offsets from tximport, using
  DESeqDataSetFromTximport().

destiny
-------

Changes in version 1.1 (2015-11-27):

- Implemented rank correlation distance and cosine distance

- Added updateObject method

DEXSeq
------

Changes in version 1.17.28:

- Added support for estimating fold change rates with continous
  variables.

- Reduced RAM usage when parallelizing

DiffBind
--------

Changes in version 2.0.0:

- Feature changes

  * Change default analysis method to DESeq2 for its more conservative
  normalization

  * Designate DESeq method as obsolete (in favor of DESeq2); alter
  documentation and vignette accordingly.

  * Change default FDR threshold to 0.05

  * Add bNot paramater to dba.contrast to remove ! contrasts by default

  * Remove bReturnPeaksets parameter from dba.plotVenn (does this by
  default)

  * Change bCorPlot default to FALSE (no more automatic clustering
  heatmaps)

- Internal changes

  * Bump version number to 2.0

  * Update vignette

  * Remove $allvectors and $vectors; replace with $merged (without
  score matrix) and $binding

  * Upgrade peaksort to use inplace peakOrder

  * Optimize peak merging memory usage

  * Change PCA method from princomp() to prcomp()

  * maxGap implemented in merge

  * Include the beginnings of some unit tests using the testthat
  package.

- Bug fixes

  * Fix bug in retrieving SummarizedExperiment

  * Fix bug when no peaks

  * Fix bugs in non-standard chromosome names and chromsome order

  * Fix bugs in Called logical vectors

  * Ensure loading sample sheet doesn't treat hash as comment char.

  * Tildes in file paths now accepted.

  * Spaces trimmed from entries in sample sheets (with warning).

  * Functions added to importFrom stats to satisfy BiocCheck.

diffHic
-------

Changes in version 1.3.5:

- Deprecated DIList objects and methods in favour of InteractionSet
  objects.

- marginCounts() now returns a RangedSummarizedExperiment for all bins.

- Added the max.height argument to the rotPlaid() and rotDI()
  functions.

- Added the diClusters() function for post-hoc cluster-level FDR
  control.

- Added the annotatePairs() function for convenient annotation of
  (clusters of) interactions.

- Fixed a bug in plotPlaid() when the interaction space was empty.

- Fixed a bug in preparePairs() where unmapped chimeric segments led to
  the loss of the entire pair.

- Updated user's guide, documentation and tests.

diffloop
--------

Changes in version 1.0.0:

- Initial release to Bioconductor.

- Added NEWS file.

- Fixes to documentation.

- Improved automated ploting including different colors

- Enhanced loop annotation

DMRcate
-------

Changes in version 1.7.2:

- Major changes. WGBS pipeline is now implemented with DSS as a
  regression step instead of limma. 450K pipeline is the same, but with
  slight cosmetic changes in anticipation of the transition to the EPIC
  array.

- DMR.plot() has been completely rewritten, now with Gviz and inbuilt
  transcript annotation for hg19, hg38 and mm10.

- DMRs are now ranked by the Stouffer transformations of the limma- and
  DSS- derived FDRs of their constituent CpG sites.

DNAshapeR
---------

Changes in version 0.99.7:

- improved the performace of sequence encoding process.

Changes in version 0.99.6:

- corrected the terms of BioView in DESCRIPTION.

Changes in version 0.99.5:

- Upgraded R dependency version to 3.3.

- Added importFrom methods in NAMESPACE.

Changes in version 0.99.4:

- Added caret package to the SUGGESTIONS in DESCRIPTION.

Changes in version 0.99.3:

- Fixed warnings from R CMD Check.

- Evaluated more vignette chunks.

Changes in version 0.99.1:

- Added documentation for the normalizeShape function. A manual entry
  is now also available for most non-exported functions.

- Added citation entry.

Changes in version 0.99.0:

- normalization parameter in encodeSeqShape.R changed to normalize.

- The y-axis plot range in plotShape can now be defined by the user.

- Completed package documentation.

- DNA shape feature matrix returned by encodeSeqShape is now normalized
  by default.

DOQTL
-----

Changes in version 1.07:

NEW FEATURES

- Added scanone.assoc to perform genome wide association mapping.

- Added calc.genoprob2 to run forward/backward algorithm once cluster
  parameters have been estimated.

CHANGES

- Fixed bug in assoc.map where the SDPs were not in the correct order.

- Change assoc.map to use the Sanger VCF file rather than our own
  custom SNP file.

DOSE
----

Changes in version 2.9.7:

- barplot accepts x and colorBy parameters as in dotplot <2016-04-13,
  Wed>

- gsfilter function for restricting result with minimal and maximal
  gene set size <2016-03-31, Thu> + see
  https://github.com/GuangchuangYu/clusterProfiler/issues/46

Changes in version 2.9.6:

- add maxGSSize parameter for hypergeometric test <2016-03-09, Wed>

- add maxGSSize parameter in GSEA analysis. <2016-03-09, Wed> + Default
  of 500 is fairly standard for GSEA analysis. + Usually if the geneset
  > 500, its probability of being called significant by GSEA rises
  quite dramatically.

- fixed R check <2016-03-05, Sat>

- update ReactomePA citation info <2016-02-17, Wed>

Changes in version 2.9.5:

- upset plot for enrichResult object <2016-01-22, Fri>

Changes in version 2.9.4:

- bug fixed in scaling category sizes in enrichMap <2016-01-10, Sun> +
  use setSize in gseResult

Changes in version 2.9.3:

- update enrichMap to scale category sizes <2016-01-04, Mon>

- update 'show' methods of enrichResult and gseaResult <2015-12-29,
  Tue>

Changes in version 2.9.2:

- re-design internal function <2015-12-20, Sun>

Changes in version 2.9.1:

- GSEA: test bimodal separately <2015-10-28, Wed>

- add NES column in GSEA result <2015-10-28, Wed>

- use NES instead of ES in calculating p-values. <2015-10-28, Wed>

- duplicated gene IDs in enrich.internal is not allow. add `unique` to
  remove duplicated ID. <2015-10-20, Tue> + see
  https://github.com/GuangchuangYu/clusterProfiler/issues/27

Changes in version 2.9.0:

- BioC 3.3 branch


easyRNASeq
----------

Changes in version 2.7.7:

- Ported version 2.6.3 changes

Changes in version 2.7.6:

- Ported version 2.6.2 changes

Changes in version 2.7.5:

- Changes by Bioc Core Dev

Changes in version 2.7.4:

- Changes by Bioc Core Dev

Changes in version 2.7.3:

- Changes by Bioc Core Dev

Changes in version 2.7.2:

- Ported version 2.6.1 changes

Changes in version 2.7.1:

- Changes by Hervé to update the NAMESPACE

Changes in version 2.7.0:

- Original Bioconductor version 3.2

Changes in version 2.6.3:

- Updated the maintainer address

Changes in version 2.6.2:

- Fixed issues with the gff3 synthetic transcript generation. Several
  mRNA lines per mRNA were kept and the feature selection was failing
  for features other than mRNA (e.g. tRNA or miRNA)

- Extended the SimpleRNASeq vignette

Changes in version 2.6.1:

- Upgraded the dependencies

- Introduced the new vignette (SimpleRNASeq) structure

- Fixed a cosmetic issue

- Corrected man pages

- Fixed issues with the synthetic transcript generation from gtf file.
  Thanks to Sylvain Foisy for reporting this one.

EBImage
-------

Changes in version 4.14.0:

NEW FEATURES

- ‘boundary’ argument to ‘filter2()’ for specifying behaviour at image
  boundaries

- the ‘hist’ method now returns a (list of) “histogram” object(s)
  (http://stackoverflow.com/q/33831331/2792099 and
  http://stackoverflow.com/a/35838861/2792099)

- ‘colormap()’ function for mapping a greyscale image to color using a
  color palette

PERFORMANCE IMPROVEMENTS

- ‘as.raster.Image’: increased performance by not transposing the image
  twice and improved support for Color images with less than 3 channels

BUG FIXES

- fixed the ‘log’ method for Image objects
  (https://support.bioconductor.org/p/74236/)

- ‘affine’: fixed handling of images containing an alpha channel
  (https://support.bioconductor.org/p/74876/)

- retain NA's in morphological operations: ‘dilate’, ‘erode’,
  ‘opening’, ‘closing’, ‘whiteTopHat’, ‘blackTopHat’,
  ‘selfComplementaryTopHat’ (https://support.bioconductor.org/p/77295/)

- fix to potential unsafe code in C function ‘affine()’ (thanks Tomas
  Kalibera)

- medianFilter.c: use proper rounding rather than truncation during
  float to int coercion

EBSeq
-----

Changes in version 1.11.1:

- Fixed a bug in EBTest() which may cause error when performing isoform
  DE testing 1 sample vs. multiple samples.

edgeR
-----

Changes in version 3.14.0:

- estimateDisp(), estimateCommonDisp(), estimateTrendedDisp(),
  estimateTagwiseDisp(), splitIntoGroups() and equalizeLibSizes() are
  now S3 generic functions.

- The default method of estimateGLMTrendedDisp() and
  estimateGLMTagwiseDisp() now only return dispersion estimates instead
  of a list.

- Add fry method for DGEList objects.

- Import R core packages explicitly.

- New function gini() to compute Gini coefficients.

- New argument poisson.bound for glmQLFTest(). If TRUE (default), the
  p-value returned by glmQLFTest() will never be less than what would
  be obtained for a likelihood ratio test with NB dispersion equal to
  zero.

- New argument samples for DGEList(). It takes a data frame containing
  information for each sample.

- glmFit() now protects against zero library sizes and infinite offset
  values.

- glmQLFit.default() now avoids passing a NULL design to .residDF().

- cpm.default() now outputs a matrix of the same dimensions as the
  input even when the input has 0 row or 0 column.

- DGEList() pops up a warning message when zero lib.size is detected.

- Bug fix to calcNormFactors(method="TMM") when two libraries have
  identical counts but the lib.sizes have been set unequal.

- Add a CRISPR-Cas9 screen case study to the users' guide and rename
  Nigerian case study to Yoruba.

EmpiricalBrownsMethod
---------------------

Changes in version 0.99.2:

- Added Kost's Method

Changes in version 0.99.1:

- Spelling corrected

Changes in version 0.99.0:

- Package released

ENCODExplorer
-------------

Changes in version 1.3.0:

NEW FEATURES

- Added quiet param to the queryEncode and searchEncode functions.

- Added the Roadmap datasets.

- Updated the encode_df object with the latest changes in the ENCODE
  database.

BUG FIXES

- Solved a bug that caused get_schemas function to fail because a
  directory was added in the ENCODE database schemas.

- Solved a bug with the searchEncode function caused by changes in the
  ENCODE REST API.

ENmix
-----

Changes in version 1.7.5:

- added function rcp for probe type bias correction

Changes in version 1.7.3:

- Bugfix in function rm.outlier

Changes in version 1.6.0:

- Improved code for parallel computing

- Added a function freqpoly

EnrichedHeatmap
---------------

Changes in version 1.1.7:

- add legends by Legend() function in the vignette

Changes in version 1.1.5:

- `normalizeToMatrix`: modified default of `trim` to 0

- `makeMatrix`: w0 mode can deal with overlapping signals

Changes in version 1.1.4:

- the color of error areas is set 25% light of the color of
  corresponding line.

- support raster image for heatmap body

Changes in version 1.1.3:

- parameter name was wrong when constructin
  `ht@layout$graphic_fun_list`

- add `row_order` option in `EnrichedHeatmap()`

- `normalizeToMatrix()` supports self-defined smoothing function

- there can be no upstream or/and downstream in the heatmap

- smoothing function can be self-defined

Changes in version 1.1.2:

- `normalizeToMatrix()`: more options can be passed to
  `locfit::locfit()`

- `anno_enrich()`: exclude NA values

Changes in version 1.1.1:

- add `getSignalsFromList()` which summarize signals from a list of
  matrix

- export `copyAttr()`

EnrichmentBrowser
-----------------

Changes in version 2.1.15:

- Slight modification on ORA's hypergeometric p-value calculation
  according to
  http://mengnote.blogspot.de/2012/12/calculate-correct-hypergeometric-p.html

Changes in version 2.1.14:

- Adapting GGEA to also deal with 2-column GRNs

Changes in version 2.1.12:

- Adapting SPIA to deal with non-kegg gene sets

Changes in version 2.1.10:

- Including EmpiricalBrownsMethod (ebm) among sbea methods

ensembldb
---------

Changes in version 1.3.20:

BUG FIXES

- methods transcripts, genes etc don't result in an error when columns
  are specified which are not present in the database and the
  return.type is GRanges.

- Removed the transcriptLengths method implemented in ensembldb in
  favor of using the one from GenomicFeatures.

Changes in version 1.3.19:

BUG FIXES

- ensDbFromGRanges (and thus ensDbFromGtf, ensDbFromGff and
  ensDbFromAH) support now Ensembl GTF file formats from version 74 and
  before.

Changes in version 1.3.18:

NEW FEATURES

- New ExonrankFilter to filter based on exon index/rank.

Changes in version 1.3.17:

BUG FIXES

- Use setdiff/intersect instead of psetdiff/pintersect.

Changes in version 1.3.16:

BUG FIXES

- Fixed failing test.

Changes in version 1.3.15:

NEW FEATURES

- GRangesFilter now supports GRanges of length > 1.

- seqlevels method for GRangesFilter.

- New methods exonsByOverlaps and transcriptsByOverlaps.

Changes in version 1.3.14:

NEW FEATURES

- seqlevelsStyle getter and setter method to change the enable easier
  integration of EnsDb objects with UCSC based packages.
  supportedSeqlevelsStyle method to list possible values. Global option
  "ensembldb.seqnameNotFound" allows to adapt the behaviour of the
  mapping functions when a seqname can not be mapped properly.

- Added a seqlevels method for EnsDb objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Add an example to extract transcript sequences directly from an EnsDb
  object to the vignette.

- Add examples to use EnsDb objects with UCSC chromosome names to the
  vignette.

BUG FIXES

- Seqinfo for genes, transcripts and exons contain now only the
  seqnames returned in the GRanges, not all that are in the database.

Changes in version 1.3.13:

NEW FEATURES

- EnsDb: new "hidden" slot to store additional properties and a method
  updateEnsDb to update objects to the new implementation.

- New method "transcriptLengths" for EnsDb that creates a similar
  data.frame than the same named function in the GenomicFeatures
  package.

BUG FIXES

- fiveUTRsByTranscript and threeUTRsByTranscript returned wrong UTRs
  for some special cases in which the CDS start and end were in the
  same exon. This has been fixed.

Changes in version 1.3.12:

NEW FEATURES

- ensDbFromGff and ensDbFromAH functions to build EnsDb objects from
  GFF3 files or directly from AnnotationHub ressources.

- getGenomeFaFile does now also retrieve Fasta files for the "closest"
  Ensembl release if none is available for the matching version.

SIGNIFICANT USER-VISIBLE CHANGES

- Removed argument 'verbose' in ensDbFromGRanges and ensDbFromGtf.

- Updated parts of the vignette.

- Removed method extractTranscriptSeqs again due to some compatibility
  problems with GenomicRanges.

BUG FIXES

- Avoid wrong CDS start/end position definition for Ensembl gtf files
  in which the start or end codon is outside the CDS.

Changes in version 1.3.11:

BUG FIXES

- "select" method returns now also the keytype as a column from the
  database.

Changes in version 1.3.10:

NEW FEATURES

- Implemented methods columns, keys, keytypes, mapIds and select from
  AnnotationDbi.

- Methods condition<- and value<- for BasicFilter.

Changes in version 1.3.9:

SIGNIFICANT USER-VISIBLE CHANGES

- The shiny app now allows to return the search results.

Changes in version 1.3.7:

SIGNIFICANT USER-VISIBLE CHANGES

- Some small changes to the vignette.

BUG FIXES

- Fixed a problem in an unit test.

Changes in version 1.3.6:

BUG FIXES

- Fixed a bug in ensDbFromGRanges.

Changes in version 1.3.5:

NEW FEATURES

- Added GRangesFilter enabling filtering using a (single!) GRanges
  object.

- Better usability and compatibility with chromosome names:
  SeqnameFilter and GRangesFilter support both Ensembl and UCSC
  chromosome names, if option ucscChromosomeNames is set to TRUE
  returned chromosome/seqnames are in UCSC format.

SIGNIFICANT USER-VISIBLE CHANGES

- Added method "value" for BasicFilter objects.

BUG FIXES

- transcripts, genes, exons return now results sorted by seq name and
  start coordinate.

Changes in version 1.3.4:

NEW FEATURES

- Added extractTranscriptSeqs method for EnsDb objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Added a section to the vignette describing the use of ensembldb in
  Gviz.

- Fixed the vignette to conform the "Bioconductor style".

- Added argument use.names to exonsBy.

BUG FIXES

- Fixed bug with getGeneRegionTrackForGviz with only chromosome
  specified.

- Fixed an internal problem subsetting a seqinfo.

Changes in version 1.3.3:

NEW FEATURES

- Add method getGeneRegionTrackForGviz to enable using EnsDb databases
  for Gviz.

BUG FIXES

- cdsBy, fiveUTRsForTranscript and threeUTRsForTranscript do no longer
  throw an error if nothing was found but return NULL and produce a
  warning.

Changes in version 1.3.2:

NEW FEATURES

- Implemented methods cdsBy, fiveUTRsForTranscript and
  threeUTRsForTranscript for EnsDb.

Changes in version 1.3.1:

BUG FIXES

- Ensuring that methods exons, genes and transcripts return columns in
  the same order than provided with argument 'columns' for return.type
  'data.frame' or 'DataFrame'.

ensemblVEP
----------

Changes in version 1.12.0:

NEW FEATURES

- add support for Ensembl release 82

- add support for Ensembl relase 84

MODIFICATIONS

- elementLengths was renamed -> elementNROWS in S4Vectors

- mark as unsupported on windows; Ensembl release 84 requires tabix

epivizr
-------

Changes in version 2.0.0:

- Move socket connection and data serving code outside of package to
  new packages.

- Use new 'epivizrServer' and 'epivizrData' packages.

- Move standalone to package 'epivizrStandalone'.

- Use simplified 'plot' and 'visualize' interface to add charts.

Changes in version 1.9.4:

- Set 'useCookie' URL parameter to FALSE so that empty workspaces are
  started

epivizrData
-----------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrData>

epivizrServer
-------------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrServer>

FamAgg
------

Changes in version 0.99.9:

BUG FIXES

- Fixes to the Vignette and addition of citation.

- Small fixes mainly to the Vignette and one left over problem from the
  git/svn conflict merge.

Changes in version 0.99.8:

NEW FEATURES

- FAData constructor recognizes \*.ped and \*.fam files and imports their
  pedigree information correctly.

- export method to export pedigree information from a FAData to a ped
  or fam file.

- Add methods getFounders and getSingletons.

SIGNIFICANT USER-VISIBLE CHANGES

- Founders are now represented by NA in columns 'father' and 'mother'.
  This fixed potential problems when IDs are character strings and not
  numeric.

- Added column 'family' to the results of the probability test.

- genealogicalIndexTest: renamed argument prune into rm.singletons.

- Removed prune argument for methods calculating per-individual
  statistics.

- Re-formated and re-structured the vignette.

BUG FIXES

- Validation of pedigree information in FAData improved.

- [ subsetting now ensures that father or mother IDs which are not
  available in column 'id' are set to NA.

- clique names are no longer dropped when setting cliques(object) <-
  value.

Changes in version 0.99.6:

NEW FEATURES

- Monte Carlo simulation to assess significance of familial incidence
  rate and familial standardized incidence rates.

- New FAIncidenceRateResults object along with its methods.

- New FAStdIncidenceRateResults object along with its methods.

- New function factor2matrix to convert a factor into a matrix.

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette: extended content and adapted to use Bioconductor style.

- Results from the kinship sum test are now, in addition to the
  p-value, also sorted by the kinship sum.

BUG FIXES

- Fixed some problems in the documentation.

Changes in version 0.99.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Some changes related to the github repository.

- Added a readme.org file.

Changes in version 0.99.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Renamed genealogicalIndex method into genealogicalIndexTest to
  conform the naming convention that all methods ending with Test use
  simulations to assess signifcance levels.

Changes in version 0.99.3:

NEW FEATURES

- Stratified sampling for kinship group test.

BUG FIXES

- Fixed a bug in plotPed related to optional labels.

- Fixed a bug in familialIncidenceRate: self-self kinship was not
  excluded.

Changes in version 0.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Vignette streamlining.

Changes in version 0.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Fixed pdf export of vignette.

BUG FIXES

- Removed the Haplopaint requirement checking during package loading as
  this resulted in errors on some systems.

Changes in version 0.99.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Improved the vignette.

- Fixed several issues in the documentation.


FindMyFriends
-------------

Changes in version 1.1.11:

- Improvement: Various speed improvements

- Fix namespace clash between Matrix and S4Vectors

Changes in version 1.1.10:

- Feature: New class pgSlim for handling pangenomes with no ref to
  sequence data

- Bug fix: safeAAread/safeDNAread would return wrong sequences when
  number of fasta files exceeded 2000

Changes in version 1.1.9:

- Feature: Threshold for core group classificaton can now be set
  (defaults to 1)

- Improvement: Only investigate the neighbors to groups that have
  changed during iteration in neighborhoodMerge.

Changes in version 1.1.8:

- Bug fix: Vignette error due to problems with reutils

Changes in version 1.1.7:

- Feature: cdhitGrouping creates initial grouping based on cdhit
  algorithm

- Feature: neighborhoodSplit now refines the splitting as a final step
  by merging highly similar groups sharing gene group up- or downstream

- Feature: gpcGrouping can now precluster using CD-Hit

- Feature: Key algorithms now reports progress and timing information

- Feature: Custom linearKernel function that takes an upper similarity
  threshold to speed up comparisons.

- Feature: Updated vignette, focusing on recommended workflow

- Improvement: More performant pcGraph, neighborhoodSplit, pgMatrix
  methods

- Improvement: pgMatrix now returns a sparseMatrix for lower memory
  footprint

- Improvement: pangenome matrix no longer stored in pgInMem

- Bug fix: Remove zerolength genes upon pangenome creation.

- Bug fix: Batch accessing fasta files to avoid "too many open
  connections" error

Changes in version 1.1.6:

- Skips some long running examples and tests to cut down on check time

Changes in version 1.1.5:

- Bug fix: .fillDefaults now correctly extracts the frame in advance

Changes in version 1.1.4:

- Fixing imports to import BiocGenerics, S4Vectors and IRanges in full

Changes in version 1.1.3:

- Minor optimization of code

- getRep now names genes by group name

- transformSim now works on sparseMatrix rather than matrix objects -
  avoids coercing huge sparse matrices down to matrix format

Changes in version 1.1.1:

- Fix bug that would affect pangenomes above 100 genomes when GPC was
  used with automatic tree creation.

Changes in version 1.1.0:

- Added kmerSplit method


flowCore
--------

Changes in version 1.37.6:

read.FCS

- supports FCS that has diverse bit widths across parameters/channels

- supports FCS that uses big integer (i.e. uint32 > R's integer.max)

write.FCS

- Fixes several bugs so that this API is now more usable


gdsfmt
------

Changes in version 1.8.0:

- the version number was bumped for the Bioconductor release version
  3.3

Changes in version 1.7.18-1.7.22:

NEW FEATURES

- LZMA compression algorithm is available in the GDS system (LZMA,
  LZMA_RA)

- faster implementation of variable-length string: the default string
  becomes string with the length stored in the file instead of
  null-terminated string (new GDS data types: dStr8, dStr16 and dStr32)

UTILITIES

- improve the read speed of characters (+18%)

- significantly improve random access of characters

- correctly interpret factor variable in `digest.gdsn()` when
  `digest.gdsn(..., action="Robject")`, since factors are not integers

Changes in version 1.7.0-1.7.17:

NEW FEATURES

- `digest.gdsn()` to create hash function digests (e.g., md5, sha1,
  sha256, sha384, sha512), requiring the package digest

- new function `summarize.gdsn()`

- `show()` displays the content preview

- define C MACRO 'COREARRAY_REGISTER_BIT32' and
  'COREARRAY_REGISTER_BIT64' in CoreDEF.h

- new C functions `GDS_R_Append()` and `GDS_R_AppendEx()` in R_GDS.h

- allows efficiently concatenating compressed blocks (i.e., ZIP_RA and
  LZ4_RA)

- v1.7.12: add a new data type: packedreal24

- define C MACRO 'COREARRAY_SIMD_SSSE3' in CoreDEF.h

- v1.7.13: `GDS_Array_ReadData()`, `GDS_Array_ReadDataEx()`,
  `GDS_Array_WriteData()` and `GDS_Array_AppendData()` return `void*`
  instead of `void` in R_GDS.h

- v1.7.15: `GDS_Array_ReadData()` and `GDS_Array_ReadDataEx()` allow
  Start=NULL or Length=NULL

- v1.7.16: new C function `GDS_Array_AppendStrLen()` in R_GDS.h

UTILITIES

- `paste(..., sep="")` is replaced by `paste0(...)` (requiring
  >=R_v2.15.0)

- improve random access for ZIP_RA and LZ4_RA, e.g., in the example of
  vignette, +12% for ZIP_RA and 1.7-fold speedup in LZ4_RA

DEPRECATED AND DEFUNCT

- bit17, ..., bit23, bit25, ..., bit31, sbit17, ..., sbit23, sbit25,
  ..., sbit31 are deprecated, and instead it is suggested to use data
  compression

BUG FIXES

- v1.7.7: fix a potential issue of uninitialized value in the first
  parameter passed to 'LZ4_decompress_safe_continue' (detected by
  valgrind)

- v1.7.14: fix an issue of 'seldim' in `assign.gdsn()`: 'seldim' should
  allow NULL in a vector

Changes in version 1.6.0-1.6.2:

- the version number was bumped for the Bioconductor release version
  3.2

- 'attribute.trim=FALSE' in `print.gdsn.class()` by default NEW
  FEATURES

- `diagnosis.gds()` returns detailed data block information BUG FIXES

- v1.6.2: it might be a rare bug (i.e., stop the program when getting
  Z_BUF_ERROR), now the GDS kernel ignores Z_BUF_ERROR in deflate() and
  inflate(); see http://zlib.net/zlib_how.html for further explanation

genbankr
--------

Changes in version 0.99.9:

USER FACING CHANGES

- use of genbankr and genbankr version are now listed in TxDb objects
  created via makeTxDbFromGenBank.

BUGFIXES

- NEWS file entries are now in most-recent-first order.

Changes in version 0.99.7:

MAJOR USER FACING CHANGES

- gene_id is now taken from locus_tag when present (in any features of
  a given type) and gene is now only used as a backup when locus_tag is
  absent for all features of that type.

USER VISIBLE CHANGES

- gene_synonym added to set of attribute that are assumed multi-valued
  and always returned as a CharacterList column (without a warning).

BUGFIXES

- Added Michael Lawrence as second author.

Changes in version 0.99.6:

USER VISIBLE CHANGES

- Added support for creating of TxDb objects from GenBankRecord objects

- Added support for genpept files to readGenBank and parseGenBank

BUGFIXES

- Fix bug where partial was not being passed down correctly when
  specified in readGenBank call.

Changes in version 0.99.5:

MAJOR USER VISIBLE CHANGES

- Move to single-class model with nullable sequence slot, class name
  changed to GenBankRecord. Removed methods, etc for old classes

- parseGenBank now accepts ret.anno and ret.seq to control whether the
  the sequence and annotations are parsed.

- Change how repeated annotation fields within a single feature entry
  is handled. Known multivalue fields (currently db_xref and EC_number)
  *always* return CharacterList. Others return CharacterList in case of
  duplicate entries only, with a warning.

BUGFIXES

- correctly support cases where arbitrary annotations appear more than
  once within a single feature.

- fix regression problem with handling GenBank files with no sequence
  information when ret.seq is TRUE

- fix problem where tests weren't being invoked properly during check

- More informative error message when non-existent filename is passed
  to readGenBank or parseGenBank.

- Fix bug when *all* variation features in a file lack /replace (for
  VRanges alt is a mandatory field).

Changes in version 0.99.4:

BUGFIXES

- Fix small documentation bugs to squash check warnings

Changes in version 0.99.2:

USER VISIBLE CHANGES

- Version number jump as part of Bioconductor submission process

- Added runnable examples to (most) help files

- Hook existing unit tests into Bioconductor testing harness

- Add fastpass to extract sequence, in the form of seq.only argument to
  parseGenBank, and getSeq methods for GBAccession and GenBankFile
  classes

BUGFIXES

- Fix some accessor methods for GenBankFull which didn't call
  annotations() to grab GenBankAnnot object before accessing slots.


genefilter
----------

Changes in version 1.54.0:

DEPRECATED AND DEFUNCT

- remove deprecated anyNA(); contradicted base::anyNA

- remove deprecated allNA()

GeneNetworkBuilder
------------------

Changes in version 1.13.3:

- remove destructor

Changes in version 1.13.2:

- debug for crash in windows

Changes in version 1.13.1:

- change int i to unsigned int i to fix the warning


genomation
----------

Changes in version 1.3.3:

IMPROVEMENTS AND BUG FIXES

- fixed reading genomic files with numeric metadata columns
  (readTableFast converted their numeric values to positive integers)

Changes in version 1.3.2:

IMPROVEMENTS AND BUG FIXES

- strands of reads in paired-end BAM files are inferred depending on
  strand of first alignment from the pair. It's a default setting of
  the strandMode argument in the readGAlignmentPairs function.

- added new argument to ScoreMatrix, ScoreMatrixBin and ScoreMatrixList
  functions library.size indicating total number of aligned reads of a
  BAM file for normalization.

- removed argument stranded from ScoreMatrix, ScoreMatrixBin and
  ScoreMatrixList functions

Changes in version 1.3.1:

IMPROVEMENTS AND BUG FIXES

- according to changes in behaviour of readr::read_delim() that e.g.
  converts "." to 0 readTableFast() was changed

GenomicAlignments
-----------------

Changes in version 1.8.0:

NEW FEATURES

- Add coercion from GAlignments or GAlignmentPairs to DataFrame, and
  from GAlignmentsList to GAlignmentPairs.

SIGNIFICANT USER-LEVEL CHANGES

- Use DESeq2 instead of DESeq in summarizeOverlaps examples (better
  late than never).

DEPRECATED AND DEFUNCT

- After being deprecated in BioC 3.2, the left() and right() getters
  and strand() setter for GAlignmentPairs objects are now defunct.

- After being deprecated in BioC 3.2, the 'invert.strand' argument of
  the first() and last() getters for GAlignmentPairs objects are now
  defunct.

- After being deprecated in BioC 3.2, the 'order.as.in.query' argument
  of the "grglist" method for GAlignmentPairs objects is now defunct.

- After being deprecated in BioC 3.2, the 'order.as.in.query' argument
  of the "rglist" and "grglist" methods for GAlignmentsList objects are
  now defunct.

- Remove the "mapCoords" and "pmapCoords" methods (were defunct in BioC
  3.2).

- Remove the readGAlignment*FromBam() functions (were defunct in BioC
  3.2).

BUG FIXES

- seqnames() setter for GAlignments objects is now consistent with
  seqnames() setter for GRanges objects.

GenomicFeatures
---------------

Changes in version 1.24:

NEW FEATURES

- Add mapRangesToIds() and mapIdsToRanges() for mapping genomic ranges
  to IDs and vice-versa.

- Support makeTxDbFromUCSC("hg38", "knownGene") (gets "GENCODE v22"
  track).

- Add pmapToTranscripts,GRangesList,GRangesList method.

SIGNIFICANT USER-VISIBLE CHANGES

- Rename the 'vals' argument of the transcripts(), exons(), cds(), and
  genes() extractors -> 'filter'. The 'vals' argument is still
  available but deprecated.

- Rename the 'filters' argument of makeTxDbFromBiomart() and
  makeTxDbPackage() -> 'filter'.

- When grouping the transcripts by exon or CDS, transcriptsBy() now
  returns a GRangesList object with the "exon_rank" information (as an
  inner metadata column).

- For transcripts with no exons (like in the GFF3 files from GeneDB),
  makeTxDbFromGRanges() now infers the exons from the CDS.

- For transcripts with no exons and no CDS (like in the GFF3 files from
  miRBase), makeTxDbFromGRanges() now infers the exon from the
  transcript.

- makeTxDbFromGRanges() and makeTxDbFromGFF() now support GFF/GTF files
  with one (or both) of the following peculiarities: - The file is GTF
  and contains only lines of type transcript but no transcript_id tag
  (not clear this is valid GTF but some users are working with this
  kind of file).  - Each transcript in the file is reported to be on
  its own contig and spans it (start=1) but no strand is reported for
  the transcript.  makeTxDbFromGRanges() now sets the strand to "+" for
  all these transcripts.

- makeTxDbFromGRanges() now recognizes features of type miRNA,
  miRNA_primary_transcript, SRP_RNA, RNase_P_RNA, RNase_MRP_RNA,
  misc_RNA, antisense_RNA, and antisense as transcripts. It also now
  recognizes features of type transposable_element_gene as genes.

- makeTxDbFromBiomart() now points to the Ensembl mart by default
  instead of the central mart service.

- Add some commonly used alternative names (Mito, mitochondrion,
  dmel_mitochondrion_genome, Pltd, ChrC, Pt, chloroplast, Chloro, 2uM)
  for the mitochondrial and chloroplast genomes to DEFAULT_CIRC_SEQS.

DEPRECATED AND DEFUNCT

- Remove the makeTranscriptDb*() functions (were defunct in BioC 3.2).

- Remove the 'exonRankAttributeName', 'gffGeneIdAttributeName',
  'useGenesAsTranscripts', 'gffTxName', and 'species' arguments from
  the makeTxDbFromGFF() function (were defunct in BioC 3.2).

BUG FIXES

- Try to improve heuristic used in makeTxDbFromGRanges() for detecting
  the format (GFF3 or GTF) of input GRanges object 'x'.

GenomicInteractions
-------------------

Changes in version 1.5.3:

NEW FEATURES

- Integration with the new InteractionSet package provides new features
  including more overlaps methods and conversion to other classes for
  storing genomic interaction data. See InteractionSet documentation
  for more details.

SIGNIFICANT USER-LEVEL CHANGES

- Following InteractionSet integration, distance calculation is done
  differently and some distances may differ by 1bp.

DEPRECATED AND DEFUNCT

- 'annotateAnchors' is deprecated and replaced by 'annotateRegions'.

GenomicRanges
-------------

Changes in version 1.24.0:

NEW FEATURES

- Add the GPos class, a container for storing a set of "genomic
  positions" (i.e. genomic ranges of width 1). Even though a GRanges
  object can be used for that, using a GPos object can be much more
  memory-efficient, especially when the object contains long runs of
  adjacent positions.

- Add a bunch of "invertStrand" methods to support strand inversion of
  any "stranded" object (i.e. any object with a strand() getter and
  setter). E.g. invertStrand() works on GRanges, GRangesList,
  GAlignments, GAlignmentPairs, GAlignmentsList, and
  RangedSummarizedExperiment objects.

- Add "is.unsorted" method for GenomicRanges objects (contributed by
  Pete Hickey).

- base::rank() gained a new 'ties.method="last"' option and
  base::order() a new argument ('method') in R 3.3. Thus so do the
  "rank" and "order" methods for GenomicRanges objects.

- Add "selfmatch" method for GenomicRanges objects.

- Add "union" method for GRangesList objects.

SIGNIFICANT USER-LEVEL CHANGES

- Remove old SummarizedExperiment class from the GenomicRanges package
  (this class is now defined in the SummarizedExperiment package).

- Move the following generic functions from the GenomicRanges package
  to the SummarizedExperiment package: - SummarizedExperiment -
  exptData, "exptData<-" - rowRanges, "rowRanges<-" - colData,
  "colData<-" - assayNames, "assayNames<-" - assays, "assays<-" -
  assay, "assay<-"

- Rename "pintersect" and "psetdiff" methods for GRangesList objects ->
  "intersect" and "setdiff" without changing their behavior (they still
  do mendoapply(intersect, x, y) and mendoapply(setdiff, x, y),
  respectively). The old names were misnomers (see svn commit message
  for commit 113793 for more information).

- Remove the ellipsis (...) from all the setops methods, except from: -
  "punion" method for signature GRanges#GRangesList; - "pintersect" and
  "psetdiff" methods for signature GRangesList#GRangesList; - "pgap"
  method for GRanges objects.

- Use DESeq2 instead of DESeq in the vignettes (better late than
  never).

DEPRECATED AND DEFUNCT

- Remove GIntervalTree class and methods (were defunct in BioC 3.2).

- Remove mapCoords() and pmapCoords() (were defunct in BioC 3.2).

GenomicTuples
-------------

Changes in version 1.5:

- New findOverlaps,GTuples,GTuples-method when type = "equal" gives
  10-100x speedup by using data.table.

- After being deprecated from GenomicRanges in BioC 3.1, mapCoords()
  and pmapCoords() are now defunct.

- After being deprecated in BioC 3.1, the "intervaltree" algorithm in
  findOverlaps() is now defunct.

genoset
-------

Changes in version 1.27.0:

- GenoSet objects now inherit from RangedSummarizedExperiment (RSE)
  rather than eSet. Use the RSE API rather than that for eSet (e.g.
  colnames rather than sampleNames and assays rather than assayData).
  All Biobase-related content has been removed including
  annotatedDataFrameFrom. This means that BigMatrix objects can no
  longer be added to eSet/ExpressionSet/GenoSet objects.

genphen
-------

Changes in version 1.0 (2016-01-01):

Introduction

- Additional functions are added sporadically.

- This news file reports changes that have been made as the package has
  been developed.

To do

- Add practical examples where genphen has been used.

- Implement multi-core execution.

- Implement genphen for categorical phenotypes.


GenVisR
-------

Changes in version 0.99.20:

- Bug fix for new code introduced in 0.99.19

Changes in version 0.99.19:

- Added option for custom sort of gene/sample in waterfall plot

- Added waterfall specific vignette

Changes in version 0.99.18:

- Minor bug fixes in cnSpec

Changes in version 0.99.17:

- Adding new feature allowing functions to export data, grob, or plots

Changes in version 0.99.15:

- Added documentation fo new function lohView

Changes in version 0.99.10:

- Updated Documentation

Changes in version 0.99.0:

- Develpment Version - Preparing for Release

ggbio
-----

Changes in version 1.19.1:

NEW FEATURES

- autoplot method supports now EnsDb objects and related filter
  objects.


ggtree
------

Changes in version 1.3.16:

- geom_treescale() supports family argument <2016-04-27, Wed> +
  https://github.com/GuangchuangYu/ggtree/issues/56

- update fortify.phylo to work with phylo that has missing value of
  edge length <2016-04-21, Thu> +
  https://github.com/GuangchuangYu/ggtree/issues/54

- support passing textConnection(text_string) as a file <2016-04-21,
  Thu> + contributed by Casey Dunn <casey_dunn@brown.edu> +
  https://github.com/GuangchuangYu/ggtree/pull/55#issuecomment-212859693

Changes in version 1.3.15:

- geom_tiplab2 supports parameter hjust <2016-04-18, Mon>

- geom_tiplab and geom_tiplab2 support using geom_label2 by passing
  geom="label" <2016-04-07, Thu>

- geom_label2 that support subsetting <2016-04-07, Thu>

- geom_tiplab2 for adding tip label of circular layout <2016-04-06,
  Wed>

- use plot$plot_env to access ggplot2 parameter <2016-04-06, Wed>

- geom_taxalink for connecting related taxa <2016-04-01, Fri>

- geom_range for adding range of HPD to present uncertainty of
  evolutionary inference <2016-04-01, Fri>

Changes in version 1.3.14:

- geom_tiplab works with NA values, compatible with collapse
  <2016-03-05, Sat>

- update theme_tree2 due to the issue of
  https://github.com/hadley/ggplot2/issues/1567 <2016-03-05, Sat>

- offset works in `align=FFALSE` with `annotation_image` function
  <2016-02-23, Tue> + see
  https://github.com/GuangchuangYu/ggtree/issues/46

- subview and inset now supports annotating with img files <2016-02-23,
  Tue>

Changes in version 1.3.13:

- add example of rescale_tree function in treeAnnotation.Rmd
  <2016-02-07, Sun>

- geom_cladelabel works with collapse <2016-02-07, Sun> + see
  https://github.com/GuangchuangYu/ggtree/issues/38

Changes in version 1.3.12:

- exchange function name of geom_tree and geom_tree2 <2016-01-25, Mon>

- solved issues of geom_tree2 <2016-01-25, Mon> +
  https://github.com/hadley/ggplot2/issues/1512

- colnames_level parameter in gheatmap <2016-01-25, Mon>

- raxml2nwk function for converting raxml bootstrap tree to newick
  format <2016-01-25, Mon>

Changes in version 1.3.11:

- solved issues of geom_tree2 <2016-01-25, Mon> +
  https://github.com/GuangchuangYu/ggtree/issues/36

- change compute_group() to compute_panel in geom_tree2() <2016-01-21,
  Thu> + fixed issue, https://github.com/GuangchuangYu/ggtree/issues/36

- support phyloseq object <2016-01-21, Thu>

- update geom_point2, geom_text2 and geom_segment2 to support
  setup_tree_data <2016-01-21, Thu>

- implement geom_tree2 layer that support duplicated node records via
  the setup_tree_data function <2016-01-21, Thu>

- rescale_tree function for rescaling branch length of tree object
  <2016-01-20, Wed>

- upgrade set_branch_length, now branch can be rescaled using feature
  in extraInfo slot <2016-01-20, Wed>

Changes in version 1.3.10:

- remove dependency of gridExtra by implementing multiplot function
  instead of using grid.arrange <2016-01-20, Wed>

- remove dependency of colorspace <2016-01-20, Wed>

- support phylip tree format and update vignette of phylip example
  <2016-01-15, Fri>

Changes in version 1.3.9:

- optimize getYcoord <2016-01-14, Thu>

- add 'multiPhylo' example in 'Tree Visualization' vignette
  <2016-01-13, Wed>

- viewClade, scaleClade, collapse, expand, rotate, flip, get_taxa_name
  and scale_x_ggtree accepts input tree_view=NULL. these function will
  access the last plot if tree_view=NULL. <2016-01-13, Wed> + >
  ggtree(rtree(30)); viewClade(node=35) works. no need to pipe.

Changes in version 1.3.8:

- add example of viewClade in 'Tree Manipulation' vignette <2016-01-13,
  Wed>

- add viewClade function <2016-01-12, Tue>

- support obkData object defined by OutbreakTools <2016-01-12, Tue>

- update vignettes <2016-01-07, Thu>

- 05 advance tree annotation vignette <2016-01-04, Mon>

- export theme_inset <2016-01-04, Mon>

- inset, nodebar, nodepie functions <2015-12-31, Thu>

Changes in version 1.3.7:

- split the long vignette to several vignettes + 00 ggtree <2015-12-29,
  Tue> + 01 tree data import <2015-12-28, Mon> + 02 tree visualization
  <2015-12-28, Mon> + 03 tree manipulation <2015-12-28, Mon> + 04 tree
  annotation <2015-12-29, Tue>

Changes in version 1.3.6:

- MRCA function for finding Most Recent Common Ancestor among a vector
  of tips <2015-12-22, Tue>

- geom_cladelabel: add bar and label to annotate a clade <2015-12-21,
  Mon> - remove annotation_clade and annotation_clade2 functions.

- geom_treescale: tree scale layer. (add_legend was removed)
  <2015-12-21, Mon>

Changes in version 1.3.5:

- bug fixed, read.nhx now works with scientific notation <2015-11-30,
  Mon> + see https://github.com/GuangchuangYu/ggtree/issues/30

Changes in version 1.3.4:

- rename beast feature when name conflict with reserve keywords (label,
  branch, etc) <2015-11-27, Fri>

- get_clade_position function <2015-11-26, Thu> +
  https://github.com/GuangchuangYu/ggtree/issues/28

- get_heatmap_column_position function <2015-11-25, Wed> + see
  https://github.com/GuangchuangYu/ggtree/issues/26

- support NHX (New Hampshire X) format via read.nhx function
  <2015-11-17, Tue>

- bug fixed in extract.treeinfo.jplace <2015-11-17, Thu>

Changes in version 1.3.3:

- support color=NULL in gheatmap, then no colored line will draw within
  the heatmap <2015-10-30, Fri>

- add `angle` for also rectangular, so that it will be available for
  layout='rectangular' following by coord_polar() <2015-10-27, Tue>

Changes in version 1.3.2:

- update vignette, add example of ape bootstrap and phangorn ancestral
  sequences <2015-10-26, Mon>

- add support of ape bootstrap analysis <2015-10-26, Mon> see
  https://github.com/GuangchuangYu/ggtree/issues/20

- add support of ancestral sequences inferred by phangorn <2015-10-26,
  Mon> see https://github.com/GuangchuangYu/ggtree/issues/21

Changes in version 1.3.1:

- change angle to angle + 90, so that label will in radial direction
  <2015-10-22, Thu> + see
  https://github.com/GuangchuangYu/ggtree/issues/17

- na.rm should be always passed to layer(), fixed it in geom_hilight
  and geom_text2 <2015-10-21, Wed> + see
  https://github.com/hadley/ggplot2/issues/1380

- matching beast stats with tree using internal node number instead of
  label <2015-10-20, Tue>

Glimma
------

Changes in version 0.99.4:

- Added sample colours for MD plot.

- Added DESeqResults method for glMDPlot

Changes in version 0.99.3:

- Imported p.adjust for relevant functions.

Changes in version 0.99.2:

- Added don't test cases for non-exported functions to pass checks.

Changes in version 0.99.1:

- Initial package creation.


GOexpress
---------

Changes in version 1.5.5:

BUG FIX

- Function subEset also removes empty levels in factors not mentionned
  in the subset argument. Necessary for randomForest as the factor of
  interest cannot have empty levels.

Changes in version 1.5.4:

NEW FEATURES

- plot_design supports subset argument to restrict the plot to certain
  subsets of samples.

Changes in version 1.5.3:

GENERAL UPDATES

- Clarify that GOexpress in not a conventional gene set enrichment
  analysis tool. It is based on ranking, not enrichment.

Changes in version 1.5.2:

BUG FIX

- Rename "namespace" column to "namespace_1003"", not "namespace_1006".

Changes in version 1.5.1:

GENERAL UPDATES

- Updated maintainer email address.

GOSemSim
--------

Changes in version 1.29.2:

- update IC data <2016-04-21, Mon>

Changes in version 1.29.1:

- fixed R check <2016-03-05, Sat>

graphite
--------

Changes in version 1.17.6 (2016-04-26):

- Added celegans pathways from KEGG.

Changes in version 1.17.5 (2016-04-25):

- Updated all pathway data.

Changes in version 1.17.4 (2016-04-04):

- Cytoscape 3 support.

Changes in version 1.17.3 (2016-04-04):

- Added the buildPathway function.


gtrellis
--------

Changes in version 1.3.6:

- change format of NEWS

Changes in version 1.3.5:

- add citation

Changes in version 1.3.4:

- identify normal chromosomes more wisely

Changes in version 1.3.3:

- revised vignette

Changes in version 1.3.2:

- support 'compact' mode to arrange chromosomes

- add customized functions to add points/lines/...

Changes in version 1.3.1:

- legend is supported

Guitar
------

Changes in version 1.9.4:

- Bug Fix

Changes in version 1.9.3:

- Quality Check; Revised Package Overview

gwascat
-------

Changes in version 2.3.6:

USER VISIBLE CHANGES

- traitsManh now emits an xlab concatenating genome tag to seqnames tag

GWASTools
---------

Changes in version 1.17.9:

- Replace ncdf with ncdf4

- Deprecate plinkToNcdf and convertVcfGds (use SNPRelate functions
  instead)

- Add function kingIBS0FSCI to define expected IBS0 spread of full
  siblings based on allele frequency.

Changes in version 1.17.8:

- do not compute qbeta for all points in qqPlot if thinning

Changes in version 1.17.7:

- add error handling to close GdsGenotypeReader and GdsIntensityReader
  gds files if they fail the validity method check

Changes in version 1.17.6:

- Use ZIP_RA as default compression in GDS files for faster access to
  compressed data

Changes in version 1.17.5:

- bug fix in checkImputedDosageFile if not writing a log file of
  missing values and an entire sample is missing from the file

Changes in version 1.17.4:

- bug fix for coloring truncated points in manhattanPlot

Changes in version 1.17.3:

- added support for hard-calling genotypes from imputed genotype
  probabilities in imputedDosageFile

Changes in version 1.17.1:

- changed colors for ibdPlot


HDTD
----

Changes in version 1.5.2 (2016-04-18):

- Updating CITATION FILE.

Changes in version 1.5.1 (2016-03-02):

- Minor fixes to vignette.

- Updating R code.

hiAnnotator
-----------

Changes in version 1.5.1:

- ggplot2 update bug fixes.

HIBAG
-----

Changes in version 1.8.0:

- the version number was bumped for the Bioconductor release version
  3.3

Changes in version 1.7.0-1.7.7:

- new functions `hlaCheckAllele()`, `hlaAssocTest()`,
  `hlaConvSequence()` and `summary.hlaAASeqClass()`

HilbertCurve
------------

Changes in version 1.1.6:

- change format of NEWS

Changes in version 1.1.5:

- use new version of testthat

Changes in version 1.1.4:

- add citation

Changes in version 1.1.3:

- revised vignette

Changes in version 1.1.2:

- add examples in vignette

- add `hc_polygon()`

- `hc_map()`: support borders

- `HilbertCurve()`: start position can be selected, also the
  orientation of the first segment can be selected

Changes in version 1.1.1:

- add functions specific for genomic data visualization

HiTC
----

Changes in version 1.15.1:

NEW FEATURES

- new saveContactMaps function

SIGNIFICANT USER-VISIBLE CHANGES

- Change default behavior of binningC function to use sum of intervals
  instead of median

- Hi-C color are now defined as a vector, so that more than 3 colors
  can be used to for the gradient

BUG FIXES

- Fix bug in import.my5C in reading file

HTSeqGenie
----------

Changes in version 4.1.15:

- investigating a memory fault produced by tallyVariants() in the
  test.callVariantsVariantTools.genotype() unit test

Changes in version 4.1.13:

- fixing NGSPIPE-191: removed mergeReads and peared references ;
  changed BioC imports ; should build on the BioC server (but cannot
  right now due to VariantTools::tallyVariants memory issues)

Changes in version 4.1.11:

- fixed NGSPIPE-182: updated the documentation of
  generateSingleGeneDERs(), highlighting that the function generates
  DEXSeq-ready exons by disjoining the whole exon set, keeping only the
  exons of coding regions and keeping only the exons that belong to
  unique genes.

Changes in version 4.1.10:

- add simple check for duplicate ids to prepocessReads. Note this will
  only catch dups within a chunk.

Changes in version 4.1.9:

- fixed a small glitch that prevented truncateReads to work with
  trimReads.trim5 config parameter only (without length truncation)

Changes in version 4.1.8:

- fixed NGSPIPE-167: added 5'-read trimming with the configuration
  parameter trimReads.trim5 (code, unit test, configuration check code,
  handling corner cases)

Changes in version 4.1.7:

- fixed the NGSPIPE-150 ticket: cleaned up readRNASeqEnds and removed
  consolidateByRead

Changes in version 4.1.6:

- fixed a buglet (absence of summary_alignment.tab file) that prevented
  test.callVariantsVariantTools.genotype() to run

- moved the loading of genomic_features in countGenomicFeatures and not
  countGenomicFeaturesChunk, for speedup

- R CMD check ok

Changes in version 4.1.5:

- fixed the NGSPIPE-172 ticket on "variant calling slow on GRCH38". The
  selftest took about 1.5 h to complete, due to the analyzeVariants()
  phase of the two test.runPipeline.Exome() tests, each lasting about
  0.6 h.  This was due to (1) a wrong loading of the huge dbSNP
  database for avgNborCount filtering (needed only if the avgNborCount
  keyword is present in the config parameter
  "analyzeVariants.postFilters". And due to (2) a costly tiling of the
  genome, which is not required for small datasets (with a number of
  reads lower than 100e6). The selftest now takes about 34 minutes to
  complete.

Changes in version 4.1.4:

- added test.readRNASeqEnds.dupmark() to test if readRNASeqEnds()
  handles properly dups (which is the case; the bug NGSPIPE-180 on
  "Ignore Dup marked reads during feature counting" is due to the fact
  that markdup is done on whole bams, while featuring counting is done
  on chunks). This revision closes NGSPIPE-180.

Changes in version 4.1.3:

- filter out alignments that are marked as dups and QC failed using the
  proper bam flag before feature counting and coverage

Changes in version 4.1.2:

- activate read quality trimming for GATK-rescaled

- move sanger quality at the end of the list of possible qualities, as
  we do not do quality trimming for real sanger.  Current and Recent
  data from illuma should always come with illumina 1.5 or 1.8 range,
  so we want to make sure the read triming is triggered.

Changes in version 4.1.1:

- createTmpDir() is now more robust and retries if dir already exists.
  this is needed as the indel realigner step uses this function to
  parallelize by chromosomes.  For genomes with lots of
  chromomes/cotigs such as GRC38 (>200) this fixed a rare bug.

Changes in version 4.1.0:

- starting new dev cycle


illuminaio
----------

Changes in version 0.13.1 (2016-01-12):

- DEPRECATED: readBPM() is deprecated because it was only a stub that
  never really worked. (Issues #5 and #6)

Changes in version 0.13.0 (2015-10-23):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.3 for R (>= 3.3.0).

ImmuneSpaceR
------------

Changes in version 0.99.0:

- Submission to BioConductor

immunoClust
-----------

1.3.3: CHANGES * Optional plotting of additional parameter in FCS files
        * Minor changes in export meta-cluster features BUGFIXES *
        Problems with clustering of 1-D data sets

1.3.1: CHANGES * The normalization step within the meta.clusering
        proces is improved and extended which is also effects the
        related parameter settings

InPAS
-----

Changes in version 1.3.2:

BUG FIXES

- re-run coverage when worker is died for bpapply

Changes in version 1.3.1:

- suppress the error not needed in phmm

INSPEcT
-------

Changes in version 1.1.4:

- Added capabilities for the comparison of two steady-state conditions

Changes in version 1.1.3:

- INSPEcT is now independent from pyhton and HTSeq for counting reads
  in intronic and exonic regions and now uses the GenomicAlignments and
  Rsamtools Bioconductor packages functionalities. Thanks to this
  modification INSPEcT is now fully contained into the R/Bioconductor
  framework, being therefore more complant with Bioconductor
  guidelines. Additionally BAM files can be directly used, without the
  time consuming and memory consuming step of BAM->SAM conversion.
  Methods related to the creation of the annotation of intronic and
  exonic features in GTF format to provide to pyhton are now dismissed
  consequently and will be probably introduced as methods in the
  Bioconductor package 'compEpiTools'

Changes in version 1.1.2:

- Parallel computation is now managed completely within the
  BiocParallel framework. "newINSPEcT" function and "modelRates" method
  take as input the argument BPPARAM which is used to handle the
  parallelization. By default, BPPARAM is assigned with bpparam()
  function of BiocParallel package, which guarantee the maximum number
  of available cores used and the usage of forking in Linux and MacOS-X
  and the usage of the package Snow for Windows machines.

- nCores methods and arguments are now deprecated.

Changes in version 1.1.1:

- Re-introduced inferKBetaFromIntegralWithPre, which disappeared in the
  devel version following 1.0.1 (excluded)

- selection of best model is now done applying brown test on the pairs
  of model where at least one of the two has a chi-squared test lower
  than the threshold. This is done because in case only one rate leads
  the dynamics, all the model which don't involve that rate won't have
  low chi-squared and no comparison will be made. This leads to brown
  p-values of 1 on that specific rates (change in method "ratePvals")

- in newINSPEcT, the guess of new rates can be done without assuming
  that degradation does not occur during the pulse

- Solved two problems. One that occurred during modeling for genes with
  estimated variance within replicates equal to zero: in these cases
  the variance is estimated within the time course. A second problem
  was encountered in the parameter initialization of impulse model: h1
  cannot be zero in order to evaluate a finite value."

- Evaluate 'modelRates' within the vignette in parallel only in Linux
  and Dawin environments. This is done to avoid timeout in the build
  process on Bioc servers.

- Better estimation of rates in case degDuringPulse=TRUE (newINSPEcT
  function). Added controls on input arguments in newINSPEcT function.
  Fixed a bug in the saturation of values out of breaks in inHeatmap
  method (this change could cause a different clustering order of genes
  in the heatmap). Added the palette argument to inHeatmap method.

- Fixed a bug in '[' method and unpdated the NAMESAPACE and DESCRIPTION
  files according to the update of the 'unlist' method that is exported
  from BiocGenerics and not from GenomicRanges anymore.

InteractionSet
--------------

Changes in version 0.99.0:

- New package InteractionSet, containing base classes for genomic
  interaction data.

IRanges
-------

Changes in version 2.6.0:

NEW FEATURES

- Add regroup() function.

SIGNIFICANT USER-VISIBLE CHANGES

- Remove 'algorithm' argument from findOverlaps(), countOverlaps(),
  overlapsAny(), subsetByOverlaps(), nearest(), distanceToNearest(),
  findCompatibleOverlaps(), countCompatibleOverlaps(),
  findSpliceOverlaps(), summarizeOverlaps(), Union(),
  IntersectionStrict(), and IntersectionNotEmpty(). The argument was
  added in BioC 3.1 to facilitate the transition from an Interval Tree
  to a Nested Containment Lists implementation of findOverlaps() and
  family. The transition is over.

- Restore 'maxgap' special meaning (from BioC < 3.1) when calling
  findOverlaps() (or other member of the family) with 'type' set to
  "within".

- No more limit on the max depth of *on-the-fly* NCList objects. Note
  that the limit remains and is still 100000 when the user explicitely
  calls the NCList() or GNCList() constructor.

- Rename 'ignoreSelf' and 'ignoreRedundant' argument of the
  findOverlaps,Vector,missing method -> 'drop.self' and
  'drop.redundant'.  The old names are still working but deprecated.

- Rename grouplength() -> grouplengths() (old name still available but
  deprecated).

- Modify "replaceROWS" method for IRanges objects so that the replaced
  elements in 'x' get their metadata columns from 'value'. See this
  thread on bioc-devel:
  https://stat.ethz.ch/pipermail/bioc-devel/2015-November/008319.html

- Optimized which.min() and which.max() for atomic lists.

- Remove the ellipsis (...) from all the setops methods, except the
  methods for Pairs objects.

- Add "togroup" method for ManyToOneGrouping objects and deprecate
  default method.

- Modernize "show" method for Ranges objects: now they're displayed
  more like GRanges objects.

- Coercion from IRanges to NormalIRanges now propagates the metadata
  columns when the object to coerce is already normal.

- Don't export CompressedHitsList anymore from the IRanges package.
  This doesn't seem to be used at all and it's not clear that we need
  it.

DEPRECATED AND DEFUNCT

- Deprecate RDApplyParams objects and rdapply().

- Deprecate RangedDataList objects.

- Deprecate the "reduce" method for RangedData objects.

- Deprecate GappedRanges objects.

- Deprecate the 'ignoreSelf' and 'ignoreRedundant' arguments of the
  findOverlaps,Vector,missing method in favor of the new 'drop.self'
  and 'drop.redundant' arguments.

- Deprecate grouplength() in favor of grouplengths().

- Default "togroup" method is deprecated.

- Remove IntervalTree and IntervalForest classes and methods (were
  defunct in BioC 3.2).

- Remove mapCoords() and pmapCoords() generics (were defunct in BioC
  3.2).

- Remove all "updateObject" methods (they were all obsolete).

BUG FIXES

- Fix segfault when calling window() on an Rle object of length 0.

- Fix "which.min" and "which.max" methods for IntegerList, NumericList,
  and RleList objects when 'x' is empty or contains empty list
  elements.

- Fix mishandling of zero-width ranges when calling findOverlaps() (or
  other member of the family) with 'type' set to "within".

- Various fixes to "countOverlaps" method for Vector#missing. See svn
  commit message for commit 116112 for the details.

- Fix validity method for NormalIRanges objects (was not checking
  anything).

isobar
------

Changes in version 1.17.1:

- Fixed critical bug in isotope impurity correction method.  Isobar
  used a transposed matrix for isotope impurity correction prior to
  this fix. Quantification results will change when isotope impurity
  correction is performed Bug discovered by Dario Strbenac
  https://support.bioconductor.org/p/74301/#79900).


kebabs
------

Changes in version 1.6.0:

- release as part of Bioconductor 3.3

Changes in version 1.5.4:

- importing apcluster package for avoiding method clashes

- improved and completed change history in inst/NEWS and package
  vignette

Changes in version 1.5.3:

- correction in prediction via feature weights for very large sparse
  explicit representation

- adaption of vignette template

- vignette engine changed from Sweave to knitr

Changes in version 1.5.2:

- correction in distance weights for mixed distance weighted spectrum
  and gappy pair kernel

- allow featureWeights as numeric vector for method
  getPredictionProfile

- correction for plot of single prediction profile without legend

- change of copyright note

- namespace fixes

Changes in version 1.5.1:

- new method to compute prediction profiles from models trained with
  mixture kernels

- correction for position specific kernel with offsets

- corrections for prediction profile of motif kernel

- additional hint on help page of kbsvm

Changes in version 1.5.0:

- devel branch created from version 1.4.0


limma
-----

Changes in version 3.28.0:

- Improved capabilities and performance for fry().  fry() has two new
  arguments, 'geneid' and 'standardize'.  The index argument of fry()
  can now be a list of data.frames containing identifiers and weights
  for each set.  The options introduced by the standardize argument
  allows fry() to be more statistically powerful when genes have
  unequal variances.  fry() now accepts any arguments that would be
  suitable for lmFit() or eBayes().  The 'sort' argument of fry() is
  now the same as for mroast().

- roast(), mroast(), fry() and camera() now require the 'design' matrix
  to be set.  Previously the design matrix defaulted to a single
  intercept column.

- Two changes to barcodeplot(): new argument 'alpha' to set
  semitransparency of positional bars in some circumstances; the
  default value for 'quantiles' is increased slightly.

- kegga() has several new arguments and now supports any species
  supported by KEGG.  kegga() can now accept annotation as a
  data.frame, meaning that it can work from user-supplied annotation.

- New functions getGeneKEGGLinks() and getKEGGPathwayNames() to get
  pathway annotation from the rest.kegg website.

- topKEGG() now breaks tied p-values by number of genes in pathway and
  by name of pathway.

- goana() now supports species="Pt" (chimpanzee).

- plotRLDF() now produces more complete output and the output is more
  completely documented.  It can now be used as a complete LDF analysis
  for classification based on training data.  The argument 'main' has
  been removed as an explicit argument because it and other plotting
  arguments are better passed using the ... facility.  New arguments
  'ndim' and 'plot' have been added.  The first allows all possible
  discriminant functions to be completed even if only two are plotted.
  The second allows dicriminant functions to be computed without a
  plot.

  plotRLDF() also now uses squeezeVar() to estimate how by how much the
  within-group covariance matrix is moderated.  It has new arguments
  'trend' and 'robust' as options to the empirical Bayes estimation
  step.  It also has new argument 'var.prior' to allow the prior
  variances to be optionally supplied.  Argument 'df.prior' can now be
  a vector of genewise values.

- new argument 'save.plot' for voom().

- diffSplice() now returns genewise residual variances.

- removeExt() has new 'sep' argument.

- Slightly improved algorithm for producing the df2.shrunk values
  returned by fitFDistRobustly().  fitFDistRobustly() now returns
  tail.p.value, prob.outlier and df2.outlier as well as previously
  returned values.  The minimum df.prior value returned by eBayes(fit,
  robust=TRUE) may be slightly smaller than previously.

- tmixture.vector() now handles unequal df in a better way.  This will
  be seen in slightly changed B-statistics from topTable when the
  residual or prior df are not all equal.

- If a targets data.frame is provided to read.maimages() through the
  'files' argument, then read.maimages() now overwrites its rownames on
  output to agree with the column names of the RGList object.

- More explicit handling of namespaces.  Functions needed from the
  grDevices, graphics, stats and utils packages are now imported
  explicitly into the NAMESPACE, avoiding any possible masking by other
  packages.  goana(), alias2Symbol(), alias2SymbolTable() now load the
  name spaces of GO.db the relevant organism package org.XX.eg.db
  instead of loading the packages.  This keeps the user search path
  cleaner.

- Various minor bug fixes and improvements to error messages.


MEAL
----

Changes in version 1.1:

- Adding plotRDAMulti and topRDAhits

- Improve plotRegion

- Improve caseExample and MEAL vignettes

- Return a sorted data.frame in correlationMethExprs

Mergeomics
----------

Changes in version 0.0-1:

- first introduction as separate package

meshr
-----

Changes in version 1.6.1:

- examples in vignette is simplified because of TIMEOUT in R CMD CHECK.

Metab
-----

Changes in version 1.5.1:

- Fixed bug in htest and fixed normalizeByBiomass in order to manually
  input biomasses.

metagene
--------

Changes in version 2.3.0:

NEW FEATURES

- Refactored the bootstrap approach to improve memory usage and reduce
  calculation time.

- Added a plot_metagene function to produce a metagene plot from a
  data.frame to avoid always having to rely on the metagene object.

- Deprecated the range and bin_size params.

- Added new sections in the vignettes: "Managing large datasets" and
  "Comparing profiles with permutations".

BUG FIXES

- Removed params that no longer works with ggplot2 > 2.2.0.

- Changed the seqlevels checks to match changes in GenomicAlignments.

- Added a check to stop early and output clear error message when
  position in a GRanges is greater than the size of a chromosome.

metagenomeFeatures
------------------

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

metagenomeSeq
-------------

Changes in version 1.13:

- Upgrade support for biom-format vs. 2.0

- Fixed issue - "MRtable, etc will report NA rows when user requests
  more features than available"

- Fixed s2 miscalculation in calcZeroComponent

MethylAid
---------

Changes in version 1.5.2:

- adapted to minfi version 1.17.8 (almost EPIC ready)

methylPipe
----------

Changes in version 1.4.3:

- updated or fixed the following functions: + findDMR: updated the
  function to allow it to run without number of processor specification

Changes in version 1.4.2:

- updated or fixed the following functions: + meth.stats: updated the
  coverage object error in this function

Changes in version 1.4.1:

- updated or fixed the following functions: + pool.reads: added this
  new function to pool reads of replicates into a single file for each
  condition

missMethyl
----------

Changes in version 1.5.1:

- New gene set testing function gsameth() added, kegg testing
  functionality added to gometh().

mosaics
-------

Changes in version 2.9.9:

BUG FIXES

- generateWig(): fix the error specific to the BAM file format, which
  generated unnecessary lines in the output WIG file.

Changes in version 2.9.8:

BUG FIXES

- construtBins(), generateWig(), extractReads(): 'readGAlignments'
  replaces the defunct 'readGAlignmentsFromBam'.

- extractReads(), findSummit(), adjustBoundary(), filterPeak(): more
  safeguards added.

- constructBins(): fix typo.

Changes in version 2.9.7:

BUG FIXES

- mosaicsPeak(): fix the error that chrID is missing.

- adjustBoundary(): fix the error that chrID is incorrectly processed.

- constructBins(), generateWig(), extractReads(): fix the error for the
  processing of BAM files for PET data.

Changes in version 2.9.6:

SIGNIFICANT USER-VISIBLE CHANGES

- export(): aveLogP is report as column 9 for 'narrowPeak' and
  'broadPeak' file formats.

Changes in version 2.9.5:

SIGNIFICANT USER-VISIBLE CHANGES

- extractReads(): Users can now choose whether to keep read-level data
  with the argument 'keepReads'.

- Object sizes are significantly decreased for the output of
  extractReads(), findSummit(), adjustBoundary(), and filterPeak().

- MosaicsPeak class definition is modified to reflect the changes
  above.

- In the peak lists, now, logMinP and logAveP (i.e., -log10
  transformation of minP and aveP, respectively) are reported instead
  of minP and aveP, respectively.

- show() method becomes significantly faster.

Changes in version 2.9.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Peak list now incorporates mean(-log10(PP)), summitSignal, and
  summit.

- In the peak list, the counts of control samples and the log ratio of
  ChIP over control counts are adjusted by the ratio of sequencing
  depth, instead of the ratio of sum of ChIP and control counts.

- postProb(): Return posterior probabilities for arbitrary peak
  regions.

- export() becomes significantly faster.

- construtBins(): calculate sequencing depth and keep this information
  in the first line (commented) of bin-level files.

- seqDepth(): returns sequencing depth information, which can be
  applied to all of BinData, MosaicsFit, MosaicsHMM, MosaicsPeak class
  objects.

- Name of method coverage() is changed to readCoverage().

BUG FIXES

- findSummit() & adjustBoundary(): fix the error that an average point
  of multiple apart summit ties is reported as a summit. Now, the first
  summit block is chosen first and then an average point of the first
  summit block is reported as a summit. Also, fix some minor numerical
  issues regarding the calculation of summit locations.

- filterPeak(): fix the error that the improvement of ChIP over control
  samples is set to zero when there is no control signal at the
  position. Now, in this case, control signal is set to zero.

- adjustBoundary(): fix the error "multiple methods tables found for
  ‘coverage’" in R CMD check.

Changes in version 2.9.3:

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette and the help documents are updated and polished.

BUG FIXES

- generateWig(): fix the error that values in the exported files are
  written in scientific notation.

- constructBins(): fix the error that values in the exported files are
  written in scientific notation.

Changes in version 2.9.2:

BUG FIXES

- extractReads(): fix the error that strands are incorrectly handle
  when loading read-level data.

- export(): fix the error to incorrectly ask to run exportReads() when
  it is not needed.

Changes in version 2.9.1:

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette and the help documents are updated and polished.

BUG FIXES

- adjustBoundary(): fix the error that boundaries are incorrectly
  adjusted.

Changes in version 2.9.0:

SIGNIFICANT USER-VISIBLE CHANGES

- extractReads(): Load read-level data and extract reads corresponding
  to each peak region.

- findSummit(): Find a summit within each peak, based on local ChIP
  profile.

- adjustBoundary(): Adjust peak boundaries (designed for histone
  modification peaks).

- filterPeak(): Filter peaks based on their peak lengths and signal
  strengths.

- mosaicsPeakHMM: Posterior decoding is set to default
  (decoding="posterior").

- mosaics package now additionally depends on GenomicRanges,
  GenomicAlignments, Rsamtools, GenomeInfoDb, and S4Vectors packages.

BUG FIXES

- export(): fix the error that values in the exported files are written
  in scientific notation.

msa
---

Changes in version 1.4.0:

- release as part of Bioconductor 3.3

Changes in version 1.3.7:

- fixes in msaPrettyPrint() function

Changes in version 1.3.6:

- msaPrettyPrint() now also accepts dashes in file names

- added section about pretty-printing wide alignments to package
  vignette

Changes in version 1.3.5:

- adaptation of displaying help text by msa() function

Changes in version 1.3.4:

- added function for checking and fixing sequence names for possibly
  problematic characters that could lead to LaTeX errors when using
  msaPrettyPrint()

- corresponding changes in documentation

- minor namespace fix

Changes in version 1.3.3:

- added function for converting multiple sequence alignments for use
  with other sequence alignment packages

- corresponding changes in documentation

Changes in version 1.3.2:

- further fixes in Makefiles and Makevars files to account for changes
  in build system

- update of citation information

Changes in version 1.3.1:

- fixes in Makefiles and Makevars files to account for changes in build
  system

Changes in version 1.3.0:

- new branch for Bioconductor 3.3 devel

MSnbase
-------

Changes in version 1.19.24:

- more unit tests and bug fixes <2016-04-30 Sat> <2016-05-02 Mon>

Changes in version 1.19.23:

- more unit tests and bug fixes <2016-04-27 Wed> <2016-04-28 Thu>

Changes in version 1.19.22:

- more unit tests <2016-04-23 Sat> <2016-04-24 Sun> <2016-04-26 Tue>

- remove readIspyData functions <2016-04-23 Sat>

- Fixed bug in error catching in
  utils.mergeSpectraAndIdentificationData <2016-04-23 Sat>

Changes in version 1.19.21:

- nadata.R unit tests and bugs fixed <2016-04-22 Fri>

Changes in version 1.19.20:

- Moved makeNaData[2] and whichNA (was in pRoloc) <2016-04-21 Thu>

Changes in version 1.19.19:

- new naplot function to visualise missing values as a heatmap and
  barplots along the samples and features. <2016-04-04 Mon>

Changes in version 1.19.18:

- fix typo in man <2016-04-02 Sat>

Changes in version 1.19.17:

- Write support of mzTab has been dropped (writeMzTabData, makeMTD,
  makePEP and makePRT are now defunct) <2016-03-18 Fri>

Changes in version 1.19.16:

- add estimateNoise,[Spectrum|MSnExp]-method; closes #78 <2016-03-10>

- import a lot of functions from recommended packages, namely graphics,
  stats and utils to avoid many "Undefined global functions or
  variables" NOTEs in R CMD check <2016-03-10>

Changes in version 1.19.15:

- limit readMSData unit test due to Windows-only error <2016-03-09 Wed>

- fix unit test for utils.colSd(x, na.rm=TRUE)

Changes in version 1.19.14:

- Document change in nQuants param fcol to groupBy <2016-03-02 Wed>

Changes in version 1.19.13:

- Fixed bug in bin_Spectrum, reported by Weibo Xie <2016-02-16 Tue>

- Added unit test for bug above <2016-02-16 Tue>

- Merged 'apply functions columnwise by group' <2016-02-16 Tue>

- In nQuants, the fcol argument has been replaced with groupBy to make
  the signature consistent with featureCV <2016-02-16 Tue>

Changes in version 1.19.12:

- Moved polarity slot from Spectrum1 (0.1.0 -> 0.2.0) to Spectrum
  (0.2.0 -> 0.3.0) superclass; also bumped Spectrum2 class (0.1.0 ->
  0.2.0) and MSnExp (0.3.0 -> 0.3.1, to trace changes in Spectrum2).
  Wrote MSnExp and Spectrum2 updateObject methods. <2016-02-11 Thu>

Changes in version 1.19.11:

- Fix trimws generic/methods with useAsDefault (see issue #75)
  <2016-02-02 Tue>

- add exprs,MSnSet-method alias (since exprs is now exported)
  <2016-02-02 Tue>

Changes in version 1.19.10:

- export exprs, fvarLabels, and sampleNames[<-] <2016-01-30 Sat>

Changes in version 1.19.9:

- new trimws method for data.frames and MSnSets <2016-01-29 Fri>

Changes in version 1.19.8:

- readMSnSet2 now also accepts a data.frame as input <2016-01-29 Fri>

- selective ggplot2 import <2016-01-29 Fri>

Changes in version 1.19.7:

- new sampleNames<- for pSet and MSnExp objects <2015-12-15 Tue>

- Fix bug preventing to write MS1 to mgf (fixes issue #73 reported by
  meowcat) <2015-12-18 Fri>

Changes in version 1.19.6:

- MSnExp feautreNames are now X01, X02 (0 after X) to maintain
  numerical sorting of ASCII characters; contributed by sgibb
  <2015-12-14 Mon>

- Update MSnbase:::subsetBy to use split instead of lapply, which makes
  topN faster. This nevertheless changes the order of the resulting
  MSnSet (see issue #63 for details and background); contributed by
  sgibb <2015-12-14 Mon>

Changes in version 1.19.5:

- Merged pull request #67 from lgatto/featureCV by sgibb: featureCV
  ignores its na.rm argument <2015-12-12 Sat>

Changes in version 1.19.4:

- Replacement method for MSnSetList names <2015-11-24 Tue>

Changes in version 1.19.3:

- New argument fcol to selectFeatureData to select feature variables
  using a vector <2015-10-28 Wed>

Changes in version 1.19.2:

- new selectFeatureData function to subset the feature data <2015-10-24
  Sat>

Changes in version 1.19.1:

- Remove generics that are defined in ProtGenerics <2015-10-15 Thu>

Changes in version 1.19.0:

- Bioc devel 3.3

MSstats
-------

Changes in version 3.3.11:

- New functionalities : calculation of the LOD and LOQ, 1)
  linear_quantlim, 2) nonlinear_quantlim, 3) plot_quantlim, and two
  example datasets, SpikeInDataLinear, SpikeInDataNonLinear are
  available.

- Update for featureSelection =‘highQuality’ in dataProcess

- allow colon(“:”) in the peptide sequence

- fix the bug for ‘fill in incomplete rows’ in dataProcess. If there
  are only one feature has incomplete rows, the issue for getting run
  and feature ID in dataProcess and not show the list. Now, it works.

- change the default for blimp in dataProcessPlots for profile plot and
  QC plot. The upper limit of y-axis with ylimUp=FALSE is calculated by
  maximum log2(intensity) across all proteins after normalization + 3
  and then rounding off to the nearest integer.

Changes in version 3.3.10:

- fix the bug for dataProcess

- When the number of proteins for $ProcessedData and $RunlevelData are
  different, the bug happened for calculating %missing and imputation.

- fix the bug for groupComparison

- when one of condition is completely missing or other special case,
  .fit.model.single can handle and output of .fit.model.single is not
  try-error. Then output for fitted and residual should be updated.

Changes in version 3.3.9:

- Condition plot from dataProcessPlots : Now condition plots are drawn
  with run-level summarized intensities per condition.

- ComparisonResult from groupComparison

- flag about missingness and imputation : Calculation for
  MissingPercentage and ImputationPercentage columns is changed 1)
  MissingPercentage : number of measured intensities/ total number of
  intensities (which is the number of features * the number of runs in
  a protein) in the conditions used for comparison (from ‘Label’
  column) by protein. Therefore different comparisons(Label in the
  output) from the same protein can have the different percentage of
  missingness.  2) ImputationPercentage : number of imputed
  intensities/total number of intensities in the conditions used for
  comparison (from ‘Label’ column) by protein. Therefore different
  comparisons(Label in the output) from the same protein can have the
  different percentage of imputation.

- new column, ‘issue’, shows special cases, such as completely missing
  in a condition or all conditions for comparisons.

- VolcanoPlot

- flag the proteins which have the condition with completely missing.
  On the left of protein name, ‘\*’ will appear in Volcano plot

Changes in version 3.3.8:

- normalization : overall median -> median of medians. For all workflow
  for MSstats, the result should not be changed. But, log(sum) will
  have slightly different result.

- flag about missingness and imputation

- RunlevelData from dataProcess include two or three more columns 1)
  NumMeasuredFeature : number of measured features in a run 2) Missing
  percentage : number of measured features/total number of features by
  run 3) NumImputedFeature : number of imputed intensities in a run.
  This column is shown only if users allow to impute the missing value.

- ComparisonResult from groupComparison : one or two columns will be
  added.  1) MissingPercentage : number of measured intensities/ total
  number of intensities (which is the number of features * the number
  of runs in a protein) by protein 2) ImputationPercentage : number of
  imputed intensities/total number of intensities by protein

Changes in version 3.3.4:

- fix the bug for featureSubset=‘highQuality’ with label-based
  experiment.

Changes in version 3.3.3:

- add new option, remove_proteins_with_interference=FALSE (default), in
  dataProcess. whether it allows to remove proteins if deem interfered
  or not.

Changes in version 3.3.2:

- ProteinName=TRUE in groupComparisonPlots shows only the name of
  significant proteins, adjusting location. ggrepel package is used.

- change featureSubset option in ‘dataProcess’ for selecting high
  quality features. featureSubset=‘highQuality’

- Fix the bug for ‘fillIncompleteRows=TRUE’ for label-based
  experiments.

- change ‘quantification’ function. use run summarization from
  dataProcess. If there are technical replicates, use median run
  summarization for each subject.

Changes in version 3.3.1:

- fix the bug for volcano plot in groupComparisonPlots, with logbase=2.

- update all plots for ggplot2

- Change the default for ‘cutoffCensored’. Now the default is
  “minFeature”.

- for imputing the censored peak intensities, remove the features which
  has only 1 measurement for survreg function.


multiClust
----------

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

mzID
----

Changes in version 1.9.1:

- Move example files into package to avoid spurious check errors in
  windows

mzR
---

Changes in version 2.5.8:

- fix compiler Warning with clang on MacOS

Changes in version 2.5.7:

- fix compilation with clang on MacOS

Changes in version 2.5.6:

- fix compilation on new windows toolchain, thanks to Kasper Daniel
  Hansen and Dan Tenenbaum

Changes in version 2.5.3:

- new pwiz.version() function returning the pwiz backend version (KK)

Changes in version 2.5.2:

- Provide pre-compiled windows libraries, again thanks to Qiang Kou
  (KK)

Changes in version 2.5.1:

- Import Boost 1.59, thanks to Qiang Kou (KK)

NOISeq
------

Changes in version 2.16.0 (2016-02-11):

- NOISeqBIO has been modified when few replicates are available and the
  computation time has been drastically removed.

- Gene clustering in NOISeqBIO when few replicates are available: It
  will be done when total number of samples is 9 or less (instead of 10
  or less).

- Fixed a bug in "biotype detection" plot. It failed when none of the
  genes in the sample had values = 0.

- Corrected an error in the calculation of standard deviation of D
  statistic in NOISeqBIO.

OncoScore
---------

Changes in version 0.99.0 (2016-03-16):

- New package OncoScore released on BioCondictor

OncoSimulR
----------

Changes in version 2.2.0:

- Plots of genotypes.

- Stacked area and stream plots (code from Marc Taylor).

- Example of modules and no epistasis.

- Removed requirement of Root in geneToModule.

- More tests (and reorganized them)

- Miscell. improvements and typos fixed in documentation and vignette.

- Added mutationPropGrowth as argument.

- Some minor bug fixes and additional checks for user errors.

Changes in version 2.1.6 (2016-04-14):

- Adapt to changes in today's release of testthat (1.0.0)

Changes in version 2.1.5 (2016-04-09):

- Added a test of driverCounts, that does not depend on OS/compiler.

Changes in version 2.1.4 (2016-04-09):

- Moved to manual tests that depend on OS/compiler (for reproduction of
  random number stream in C++).

Changes in version 2.1.3 (2016-04-04):

- Fixed sporadic bug in countDrivers

Changes in version 2.1.2 (2016-03-27):

- Arguments to BNB_Algo5 explicit.

- Example of modules and no epistasis.

- Removed requirement of Root in geneToModule.

- More tests (and reorganized them)

- Miscell. improvements in documentation and vignette.

Changes in version 2.1.1 (2016-03-07):

- Added mutationPropGrowth as argument.

- Stacked area and stream plots (code from Marc Taylor).

- Plots of genotypes.

- Expanded vignette.


OrganismDbi
-----------

Changes in version 1.14.0:

MODIFICATIONS

- replace www.biomart.org with www.ensembl.org

- import 'mcols', 'mcols<-' from S4Vectors

- follow name change for GenomicFeatures:::.set_group_names()

- add biomaRt, rtracklayer to 'Suggests'; used in unit tests/man pages

- elementLengths was renamed -> elementNROWS in S4Vectors

- replace require() with requireNamespace()

- adjustments in response to the 'vals' -> 'filter' renaming in
  GenomicFeatures

- update unit tests to reflect new PFAM data

- load RSQLite in unit tests; no longer free from
  AnnotationDbi::dbFileConnect

- use newly exported functions from AnnotationDbi related to select()
  and building annotation packages

PAA
---

Changes in version 1.5.1 (2015-10-20):

GENERAL

- No changes.

NEW FEATURES

- The new function plotArray() is available. It can be used for visual
  inspection of ProtoArrays as well as for the monitoring of the impact
  of any pre-processing method applied to the data. In order to support
  the plotting of ProtoArrays before duplicate aggregation and mimic
  the original scan image of ProtoArrays, now loadGPR() supports the
  option protoarray.aggregation="none". Finally, the vignette has been
  updated in order to describe the new function via an exemplary
  ProtoArray showing a spatial bias.

IMPROVEMENTS

- No changes.

MODIFICATIONS

- No changes.

BUG FIXES

- Since the exemplary data have been imported with the bug which has
  been fixed in PAA version 1.3.3, all exemplary data sets have been
  reimported with the latest version of loadGPR() and saved as RData
  files. Now all exemplary data sets are complete.

PanVizGenerator
---------------

Changes in version 0.99.0:

- Submission for Bioconductor

- Added vignette

- Added panviz method for use of functionality in R (instead of just
  shiny)

- Make gene ontology a subsequent download instead of bundle with
  package

- Check usability of current GO before building PanViz

Path2PPI
--------

Changes in version 1.1.3:

- Uploaded citation file

Changes in version 1.1.2:

- Change of maintainers' email adress

- Minor changes in the package vignette (e.g. section arrangement)

- Updated citation

Changes in version 1.1.1:

- Changes in package titel and description

- Changes and additional explanations in the package vignette

pathview
--------

Changes in version 1.10.2:

- fixed problem that no node mapped when one gene/protein ID maps to
  multiple Entrez Gene ID (like Enzyme IDs): id2eg call with
  unique.map=F.

Changes in version 1.10.1:

- fix bug in reaction2edge function, which throw error for a few
  metabolic pathways with no multi-node reaction group (examples
  including 00072, 00100).

Pbase
-----

Changes in version 0.11.3:

- acols and pranges replacement method <2016-02-22 Mon>

Changes in version 0.11.2:

- add addPeptideFragments; see #24 for details [2016-02-07 Sun]

Changes in version 0.11.1:

- elementLengths was renamed -> elementNROWS in S4Vectors (new name
  reflects TRUE semantic, old name will be deprecated soon) [r113044 |
  hpages@fhcrc.org | 2016-01-29 01:22:03 +0000 (Fri, 29 Jan 2016)]

Changes in version 0.11.0:

- New Bioc devel version

pbcmc
-----

Changes in version 0.99.5:

DEPENDENCIES

- R (>= 3.3.0) now added (Thanks to Valerie Oberchain).

Changes in version 0.99.4:

DEPENDENCIES

- Depends, Imports and NAMESPACES libraries were modified following
  codetoolsBioC suggestions (Thanks to Valerie Oberchain).

- `cowplot` was used in subjectReport function as textGrob update broke
  the code.

DOCUMENTATION

- Documentation updated (Thanks to Valerie Oberchain).

CODE

- Modification according to Bioconductor style: -@ only used within
  accessor functions.  -No requiere(.) inside functions -prototype and
  validity for each class definition.  (Thanks to Valerie Oberchain).

Changes in version 0.99.3:

DEPENDENCIES

- `BiocParallel (>= 1.3.13)` has been updated to feedback the user with
  progressbar if verbose=TRUE in permutate (Thanks to Valerie
  Oberchain).

DESCRIPTION

- `StatisticalMethod` has been removed from biocViews.

Changes in version 0.99.2:

DEPENDENCIES

- `R (>= 3.x.y)` has been updated to 3.2

- Imports: `BiocParallel` has replaced `parallel` (Thanks to Valerie
  Obenchain).

Changes in version 0.99.1:

DOCUMENTATION & CODE  

- Minor modifications to cope with BiocCheck policies.

- RUnit tests were added.

Changes in version 0.99.0:

DOCUMENTATION

- `NEWS` file was added.

- First functional version

pcaExplorer
-----------

Changes in version 0.99.1:

OTHER NOTES

- Changed format of the NEWS file

Changes in version 0.99.0:

OTHER NOTES

- Ready for submission to Bioconductor


piano
-----

Changes in version 1.12:

- No changes yet

Changes in version 1.10.2:

DOCUMENTATION

- Updated the vignette R code to avoid an error in the call to biomart
  throught the getBM function.

Changes in version 1.10.1:

BUG FIXES

- Fixed a bug in runGSA which returned wrong numbers for the up- and
  down- regulated genes, for the GSEA method.

podkat
------

Changes in version 1.3.1:

- added missing method readVariantInfo() for signature 'character',
  'GRanges'

- minor streamlining of source code of readGenotypeMatrix()

- corrections of namespace imports

Changes in version 1.3.0:

- new branch for Bioconductor 3.3 devel

polyester
---------

1.99.3: NB function now exported

1.99.3: note that version 1.99.3 on GitHub was version 1.1.0 on
        Bioconductor.

1.99.2: bug fix in fragment generation (last 2 bases of transcript were
        never sequenced)


pqsfinder
---------

Changes in version 1.0.0:

- Novel algorithm for identification of potential intramolecular
  G-quadruplex (G4) patterns in DNA sequence.

- Supports multiple defects in G-runs like bulges or mismatches.

- Provides the most accurate results currently available.

- Highly customizable to detect even novel G4 types that might be
  discovered in the future.

procoil
-------

Changes in version 2.0.0:

- updated models PrOCoilModel and PrOCoilModelBA that have been trained
  with newer data and up-to-date methods

- general re-design of classes and functions to allow for multiple
  predictions per run; that comes along with a more streamlined and
  more versatile interface how to supply sequences and registers to
  predict().

- predictions and plots are now performed by the 'kebabs' package; this
  led to a major performance increase.

- the integration of the 'kebabs' package also allowed for inheriting
  functions like heatmap() and accessors, such as, sequences(),
  baselines(), and profiles()

- added a fitted() method to allow for easy extraction of predictions

- addition of small example model file inst/examples/testModel.CCModel

- streamlining/simplification of some man pages

- several corrections and updates of man pages and package vignette

- changed vignette building engine from Sweave to knitr

- removal of reference to Git-SVN bridge


pRoloc
------

Changes in version 1.11.23:

- New gomarkers functionality for adding annotation information to
  spatial proteomics data and accompanying new vignette <2016-04-18
  Mon>

- Added unit tests <2016-04-20 Wed> <2016-04-21 Thu>

- Moved makeNaData[2] and whichNA to MSnbase <2016-04-21 Thu>

- Renamed addGoMarkers and orderGoMarkers to addGoAnnotations and
  orderGoAnnotations and all associated documentation.  Vignette also
  renamed to pRoloc-goannotations <2016-04-21 Thu>

Changes in version 1.11.22:

- fix bug in knntlOptimisation to allow passing of th matrix with 1
  column <2016-04-08 Fri>

Changes in version 1.11.21:

- fix bug in getParams method = 'max' <2016-04-07 Thu>

Changes in version 1.11.20:

- Update plotDist signature to support different types and pch
  <2016-04-05 Tue>

Changes in version 1.11.19:

- Update dunkley2006params <2016-04-01 Fri>

- Update dunkley2006params <2016-04-01 Fri>

Changes in version 1.11.18:

- Update dunkley2006params <2016-03-30 Wed>

- Update dunkley2006params <2016-03-30 Wed>

Changes in version 1.11.17:

- Selective imports <2016-03-20 Sun>

- Selective imports <2016-03-20 Sun>

Changes in version 1.11.16:

- Update package startup msg <2016-03-11 Fri>

Changes in version 1.11.15:

- Clarify that score during optim are not a reflection of final
  assignment accuracy <2016-03-01 Tue>

Changes in version 1.11.14:

- fix build error due to doi/url confusion <2016-03-01 Tue>

Changes in version 1.11.13:

- new method argument added to knntlOptimisation that allows
  optimisation of class weights as per Wu and Dietterich's original
  k-NN TL method <2016-02-19 Fri>

- seed argument added to knntlOptimisation for reproducibility
  <2016-02-22 Mon>

- New section in tl vignette describing preparation of auxiliary PPI
  data <2016-02-29 Mon>

Changes in version 1.11.12:

- the colour and pch setters now invisibly return the old values
  <2016-02-17 Wed>

Changes in version 1.11.11:

- added error message when sampleNames differ between datasets in an
  MSnSetList when using remap function <2016-02-11 Thu>

Changes in version 1.11.10:

- Update colours man page to document change in default colour palette
  <2016-02-08 Mon>

Changes in version 1.11.9:

- Lisa's colour palette is now default. Old colours can be accessed and
  set with get/setOldcol <2016-02-04 Thu>

Changes in version 1.11.8:

- New Lisa cols and changed default unknown col <2016-02-03 Wed>

- mrkVecToMat has been updated so that the column order reflects the
  factor levels of fcol, rather than calling unique on fcol.  This
  change means that the order of the classes in fcol are now consistent
  between plot2D and new visualisation apps that rely on mrkVecToMat.
  <2016-02-03 Wed>

Changes in version 1.11.7:

- Various non-visible changes. <2016-01-19 Tue>

Changes in version 1.11.6:

- new classWeights function <2015-12-26 Sat>

Changes in version 1.11.5:

- highlightOnPlot support labels = TRUE to use featureNames as labels
  <2015-12-21 Mon>

- selective ggplot2 import <2015-12-21 Mon>

- highlightOnPlot also support a vector of feature names in addition to
  an instance of class FeaturesOfInterest <2015-12-21 Mon>

Changes in version 1.11.4:

- plot2D: Mirror PCs even when not plotting, addressing issue #68
  <2015-12-18 Fri>

Changes in version 1.11.3:

- Update dunkley2006params to use plant_mart_30 <2015-12-16 Wed>

- API change in plot2D: to plot data as is, i.e. without any
  transformation, method can be set to "none" (as opposed to passing
  pre-computed values to method as in previous versions). If object is
  an MSnSet, the untransformed values in the assay data will be
  plotted. If object is a matrix with coordinates, then a matching
  MSnSet must be passed to methargs. <2015-12-16 Wed>

Changes in version 1.11.2:

- Internally using MartInterface to query individual mart servers to
  bypass the biomart.org downtime <2015-12-09 Wed>

Changes in version 1.11.1:

- New orgQuants function and update to getPredictions <2015-10-13 Tue>

- Deprecate minClassScore replaced by getPredictions <2015-10-19 Mon>

- Add pRolocVisMethods and check for new apps in pRolocGUI <2015-10-19
  Mon>

- new fDataToUnknown function <2015-10-23 Fri>

- New section in vignette describing readMSnSet2 <2015-11-30 Mon>

Changes in version 1.11.0:

- Bioc devel version 3.3

pRolocGUI
---------

Changes in version 1.5.6:

- Removing plotMat2D app (closes issue #69) <2016-03-11 Fri>

- add package startup msg <2016-03-11 Fri>

- instruct users to install latest version from github

Changes in version 1.5.5:

- Depend on DT >= 0.1.40

Changes in version 1.5.4:

- replace getLisacol by getStockcol (which are now Lisa's colours;
  since pRoloc version 1.11.9) <2016-02-16 Tue>

Changes in version 1.5.3:

- update refs in vignette <2016-02-09 Tue>

Changes in version 1.5.2:

- Updated pca app <2016-01-11 Mon>

- Updated vignette <2016-01-12 Tue>

- Fixed bugs, pca app renamed main app, removed profiles app
  <2016-01-14 Thu>

- new compare app <2016-01-30 Sat>

- updated vignette <2016-02-03 Wed>

Changes in version 1.5.1:

- New shiny apps <2015-10-12 Mon>

- New vignette <2015-10-29 Thu>

- Fixed bugs and updated examples in classify app <2015-11-09 Mon>

Changes in version 1.5.0:

- Bioc devel 3.3

ProteomicsAnnotationHubData
---------------------------

Changes in version 1.1.2:

- Major rewrite of the data preparation, now relying on simple dcf
  files as input and intermediate PAHD objects. <2015-12-17 Thu>

- import read.table <2016-03-29 Tue>

- replace curl code by RCurl::getURL <2016-03-29 Tue>

Changes in version 1.1.1:

- Manually setting flInfo in .ftpfileinfo2 for files located on S3
  <2015-12-15 Tue>

Changes in version 1.1.0:

- Bioc devel 3.3

ProtGenerics
------------

Changes in version 1.3.3:

- exporting mz<- generic <2015-10-19 Mon>

Changes in version 1.3.2:

- Not depending on BiocGenerics <2015-10-16 Fri>

Changes in version 1.3.1:

- added chomatogram generic (previously defined in MSnbase) <2015-10-15
  Thu>

Changes in version 1.3.0:

- Bioc release 3.3

pwOmics
-------

Changes in version 3.1.1:

- include static consensus profiles function

- include activating/inhibiting edges in dynamic consensus net

qcmetrics
---------

Changes in version 1.9.1:

- QcMetric objects now have a description slots, and associated
  accessor and replacement methods <2015-11-18 Wed>

QDNAseq
-------

Changes in version 1.8.0 (2016-04-15):

RELEASE

- Bioconductor 3.3

IMPROVEMENTS

- estimateCorrection(), segmentBins(), createBins(), and
  calculateBlacklist() now support parallel computing (see vignette for
  more details)

- callBins() can now also use cutoffs instead of CGHcall

- binReadCounts() now contains parameter pairedEnds to specify when
  using paired-end data, so that expected variance can be calculated
  correctly

- segmentBins() now allows seeds for random number generation to be
  specified for reproducibility

- binReadCounts() supports chunked processing of bam files

- estimateCorrection() now also allows correcting for only GC content
  or mappability

BUG FIXES

- applyFilters() and highlightFilters() now work properly when using a
  numerical value for parameter residual

- highlightFilters() no longer highlights entire chromosomes for which
  the residual filter is missing altogether, which matches the behavior
  of applyFilters()

- getBinAnnotations() now allows custom bin annotations to be loaded
  via the path parameter even when an annotation package has been
  installed

- phenodata files with a single variable are now handled correctly

- calculateMappability() now retains correct chromosome order even when
  bigWigAverageOverBed reorders them

- calculateBlacklist() now correctly handles non-integer chromosome
  names

qpgraph
-------

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

QUBIC
-----

Changes in version 1.0.0:

- First Bioconductor upload

R3CPET
------

1.4.0: Updates: * Fixed some import issues

RCyjs
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

ReactomePA
----------

Changes in version 1.15.6:

- fixed issue of duplicated name in PATHID2NAME <2016-04-21, Thu>

Changes in version 1.15.5:

- add maxGSSize parameter <2016-03-10, Thu>

- update ReactomePA citation info <2016-02-17, Wed>

Changes in version 1.15.4:

- according to the update of DOSE <2015-12-20, Sun>

Changes in version 1.15.3:

- update citation info <2015-11-24, Tue>

Changes in version 1.15.2:

- update vignette <2015-11-16, Mon>

Changes in version 1.15.1:

- update vignette <2015-10-28, Wed>

recoup
------

Changes in version 0.99.4 (2016-04-02):

NEW FEATURES

- Added support for BigWig files which greatly increases coverage
  reading speed.

- Added plyr import to reduce memory footprints in some averaging
  calculations.

- Added another plot type which shows correlation between average
  coverage in summarized genomic regions and respective plot control
  parameters.

- Plotting of confidence intervals for profile (and the newly added
  correlation) plots (geom_ribbon) is now an option.

- The setter function now supports mulitple argument setting at once.

- Object slicing is now also performed on genomic position instead of
  only reference regions and samples.

- Stopped automatic width in heatmaps and passed control to the user
  through the use of ... and also using ... for plots rendered with
  ggplot.

- Moved sumStat and smooth options from binParams to plotParams as
  smoothing should be available for reusing/replotting recoup objects.
  sumStat remained in binParams to be used for region binning over e.g.
  gene bodies.

- Added documentation for recoup_test_data

- Added small BAM files for testing of the preprocessRanges())
  function.

- Updated vignettes.

BUG FIXES

- Fixed bug when reading reads from bed files. GenomeInfoDb is used to
  fill the seqinfo slot of the produced GenomicRanges. Credits to
  Orsalia Hazapis, BSRC 'Alexander Fleming'.

- Fixed bug when region is "custom" and the intervals were not of fixed
  width. Credits to Orsalia Hazapis, BSRC 'Alexander Fleming'.

- Fixed bug in custom heatmap ordering.

- Fixed bug in calculation of average profiles in recoupCorrelation
  when using mean/median instead of splines.

Changes in version 0.4.0 (2016-02-02):

NEW FEATURES

- Removed bigmemory/biganalytics as storage is a trouble.

- Added (almost) full reusability of recoup objects. Now the most
  serious, memory and time consuming calculations need to perform just
  once.

- Added slicing/subsetting or recoup list object.

- Broke the code in more collated files.

- k-means design function is now much more flexible.

BUG FIXES

- Fixed bug with colors in curve profile single-design plots.

Changes in version 0.2.0 (2016-01-28):

NEW FEATURES

- Added a global sampling factor to reduce the size of the total reads
  and genomic areas for fast overviews.

- Coverages and profiles are not recalculated when only ordering
  (k-means or other) changes.

- Exported kmeansDesign and removeData functions.

- Added full reusability of the output list object of recoup. Upon a
  call using this object, only the required elements are recalculated
  according to the change in input parameters, saving a lot of
  computation time for simple plotting/profile ordering/binning
  changes.

- Reduced binned coverage calculation time by two

- Switched to bigmemory and biganalytics packages to handle coverage
  profile matrices.

- Added sample subsetting when plotting. In this way, profiles are
  computed once and the user may choose what to plot later.

- Added a simple setter and getter for easier manipulation and
  reusability of recoup objects.

BUG FIXES

- Fixed problem with title and curve sizes in some plots

Changes in version 0.1.0 (2016-01-18):

NEW FEATURES

- First release

regionReport
------------

Changes in version 1.5.48:

BUG FIXES

- Fixed a bug in derfinderReport() for a case when there are
  significant regions but not all regions have finite areas.

Changes in version 1.5.43:

SIGNIFICANT USER-VISIBLE CHANGES

- edgeReport() now includes two edgeR specific plots: one showing the
  BCV and another showing a 2-dim MDS. Also added more edgeR citations
  that I missed earlier: thank you Gordon Smyth!

Changes in version 1.5.33:

NEW FEATURES

- Added the function edgeReport() for creating HTML or PDF reports
  based on edgeR results. Together with DESeq2Report() now regionReport
  supports the two most used packages for RNA-seq feature-level
  analysis.

Changes in version 1.5.19:

NEW FEATURES

- Added the templates 'templatePvalueHistogram' and 'templateHistogram'
  to be used with renderReport() if you prefer histogram plots instead
  of density plots.

Changes in version 1.5.12:

NEW FEATURES

- Added the function DESeq2Report() for creating HTML or PDF reports
  based on DESeq2 results. This should also be useful to explore
  derfinder results created from the expressed regions-level approach.

SIGNIFICANT USER-VISIBLE CHANGES

- Added a 'digits' argument to control how to round some numerical
  variables in all type of reports.

- Added a 'theme' argument to allow setting the ggplot2 theme for the
  plots.

BUG FIXES

- Improved the PDF versions of all reports by hiding code and
  shortening tables. Also added a warning to switch the device to 'pdf'
  for PDF output since it looks better than the default 'png'.

Changes in version 1.5.6:

SIGNIFICANT USER-VISIBLE CHANGES

- Switched to using rmarkdown instead of knitrBootstrap as the default
  engine for creating the reports.

rGREAT
------

Changes in version 1.3.3:

- documentations were improved.

- test files are moved to tests/testthat

Changes in version 1.3.2:

- bug fixed when background is self-provided.

Changes in version 1.3.1:

- fixed a bug when gene symbol annotated by GREAT is ''

- print warnings when GREAT gives warnings


rhdf5
-----

Changes in version 2.16.0:

NEW FEATURES

- New access of HDF5 files by file, group and dataset handles. HDF5
  groups and datasets can be read and written by the $-operator (e.g.
  h5f$A) and the [-operator can be used for partial reading and writing
  of datasets (e.g. h5d[3,,]).

- New low level general library function H5Dget_create_plist
  implemented.

- Removed #include <R.H> from external C code. To be compatible with
  newest C-compilers and current R-devel

RiboProfiling
-------------

Changes in version 1.1.6:

- replaced IRanges::as.table since it is not exported by IRanges
  anymore, with as.table

- motifSize in countShiftReads and codonInfo can be only 3, 6 or 9

- in countShiftReads and codonInfo, motifs of size 6 or 9 correspond to
  overalpping windows, sliding 1 codon at a time

- in countShiftReads, for motifSize = 6, the P-site is the first codon
  of the 2 it means that only the reads of the 1st codon in the 2
  codon-motif are counted

- in countShiftReads, for motifSize = 9, the P-site is the second codon
  of the 3 it means that only the reads of the 2nd codon in the 3
  codon-motif are counted

- adapted the RiboProfiling.Rnw to explain the treatment of motifs of
  codons (1, 2, or 3) as opposed to the previous treatment of unique
  codons.

- changed S4Vectors::elementLengths to S4Vectors::elementNROWS in
  riboSeqFromBAM, aroundPromoter

Changes in version 1.1.5:

- replaced GenomicRanges::assays since it is not exported by
  GenomicRanges anymore, by SummarizedExperiment::assays, and
  importFrom

Changes in version 1.1.4:

- replaced GenomicRanges::assays since it is not exported by
  GenomicRanges anymore, by assays

Changes in version 1.1.3:

- function "codonInfo" modified optimized for motifs of 1 or 2
  consecutive codons

- getCodons (utils), countShiftReads, and codonInfo functions modified
  it allows to look at sequence motifs not only for codons (3
  nucleotides) but for other sizes as well motifSize parameter has been
  added to these functions to specify the number of nucleotides on
  which to compute usage and coverage.

- function "readsToReadStart" replaced by "readsToStartOrEnd" this new
  function focuses not only the read 5' (start) but also the read 3'
  (end)

- In riboSeqFromBAM added a parameter offsetStartEnd. It specifies if
  the offset is to be computed on the 5' or 3' of the read.

Changes in version 1.1.2:

- added slope to geom_abline function

Changes in version 1.1.1:

- Added RiboSeq biocViews

RnBeads
-------

Changes in version 1.3.8:

- Added age prediction from 450K and bisulfite data to inference module

Changes in version 1.3.7:

- Subsampling for reporting and plotting methylation values
  distributions in filtering report to reduce memory footprint.

Changes in version 1.3.6:

- Support for Illumina EPIC array and related bugfixes

- Added a section on whitelist and blacklist in the Preprocessing
  module

Changes in version 1.3.5:

- Fixes with ggplot in 450K array QC plotting

- Added an option 'differential.report.sites' to decrease runtime by
  skipping the differential site methylation section in the report

Changes in version 1.3.4:

- Added support for Illumina EPIC array

Changes in version 1.3.3:

- Added 'smooth.profile' parameter to rnb.plot.locus.profile function

Changes in version 1.3.2:

- Performance and memory improvements during loading, qc and
  preprocessing of large bisulfite-derived RnBSets

- Added optional arguments to meth() and covg() allowing for partial
  retrieval of methylation and coverage information from RnBSets.
  Particularly makes sense for large, disk-based datasets.

- Added option "disk.dump.bigff.finalizer" for setting the finalizer
  for temporary BigFfMat objects

Changes in version 1.3.1:

- Added "remove.regions" function for RnBSet class

- summarize.regions now summarizes for each sample individually in
  order to reduce memory consumption for large WGBS datasets

- import.skip.object.check option for keeping the memory profile low
  while loading huge datasets

- NA sites (filtering.missing.value.quantile) are now removed even if
  they were previously masked (filtering.low.coverage.masking)

- Added "nsites" method for quickly extracting the number of sites for
  large RnBSet objects without having to retrieve the full methylation
  matrix

- Added "hasCovg" method for quickly determining whether coverage
  information is present in large RnBSet objects without having to
  retrieve the full coverage matrix

roar
----

Changes in version 1.7.6:

BUG FIXES

- Now correcly using only reads falling on PRE portions in methods
  countResults and fpkmResults also for RoarDatasetMultipleAPA objects.

Changes in version 1.7.2:

- New code to efficiently consider more than a single APA foreach gene
  - under active development (man pages, vignettes and UnitTests are
  being updated)

- Tests now performed with testthat BUG FIXES

- No changes classified as 'bug fixes'

Changes in version 1.7.1:

NEW FEATURES

- Now the gtf can contain an attribute that represents the lengths of
  PRE and POST portions on the transcriptome.  If this is omitted the
  lengths on the genome are used (like before).  Note that right now
  every gtf entry (or none of them) should have it.

BUG FIXES

- No changes classified as 'bug fixes'

rols
----

Changes in version 1.99.1:

- Explicitly extract raw content to fix error <2016-04-02 Sat>

Changes in version 1.99.0:

- Using the REST API

Changes in version 1.13.1:

- Package deprecation message <2015-12-30 Wed>

Changes in version 1.13.0:

- Bioc devel 3.3

ROntoTools
----------

Changes in version 2.0.0:

- added a novel pathway analysis method based on the primary
  dis-regulation of the genes in each pathway

ropls
-----

Changes in version 1.3.20:

NEW FEATURES

- 'opls' is now an S4 class: please use the accessors e.g.
  'getScoreMN(oplsModel)' instead of 'oplsModel$scoreMN', or
  'getScoreMN(oplsModel, orthoL = TRUE)' instead of
  'oplsModel$orthoScoreMN'

- The full list of accessors is: getLoadingMN, getScoreMN, getVipVn,
  getWeightMN (all with the orthoL argument), in addition to
  getSummaryDF, getPcaVarVn, getSubsetVi

- Please see the vignette and documentation for examples

Changes in version 1.3.19:

BUG FIXED

- Committing the new files

Changes in version 1.3.18:

NEW FEATURE

- 'coefficients' method renamed as 'coef'

BUG FIXED

- 'methods-opls_class.R' file renamed as 'opls-methods.R' to ensure
  'opls' class is defined before the methods

Changes in version 1.3.17:

BUG FIXED

- Committing the new files

Changes in version 1.3.16:

BUG FIXED

- Committing the new files

Changes in version 1.3.15:

NEW FEATURES

- Switching to S4 class

- Evaluating the metabolomics workflow in the vignette

Changes in version 1.3.14:

NEW FEATURES

- plot.opls: when 'parAsColFcVn' is a character vector, it is converted
  to a factor (thus drawing ellipses by default)

BUG FIXED

- plot.opls: handling of ellipse display when parEllipseL parameter set
  to NA (default)

Changes in version 1.3.13:

NEW FEATURES

- opls: no scaling option now available

BUG FIXED

- opls: computation of the matrix of correlations (when the number of
  features is <= 100)

Changes in version 1.3.12:

NEW FEATURES

- opls: multiclass PLS-DA implemented (PLS2 approach; comments added in
  the vignette) default number of permutations set to 20 (instead of
  10) predictive components denoted in the tables by 'p' (instead of
  'h' previously) OPLS(-DA): simplified "modelDF" data frame
  ('rotation' row has been suppressed)

BUG FIXED

- opls: PCA model with 'svd' and a single component

Changes in version 1.3.11:

NEW FEATURES

- opls: for OPLS(-DA), vipVn and orthoVipVn are now computed as the
  VIP4,p and VIP4,o described in Galindo-Prieto et al (2014)

- plot.opls: changes in palette: black/grey colors for diagnostics and
  other colors for scores

Changes in version 1.3.10:

BUG FIXED

- plot.opls: minor bug fixed ('x-score' label and color display of test
  samples when 'subset' is not NULL)

Changes in version 1.3.9:

NEW FEATURES

- vignette: minor update

Changes in version 1.3.8:

NEW FEATURES

- vignette: "some words of warning" section added

Changes in version 1.3.7:

NEW FEATURES

- vignette: default figures from plot.opls included

Changes in version 1.3.6:

NEW FEATURES

- strF: object size now displayed in Mb (instead of bytes); minor
  corrections to handle all matrices and data frames (whatever the
  dimensions, mode, row and column names)

- unit tests: new tests added to increase test coverage

Changes in version 1.3.5:

NEW FEATURES

- unit tests: warning message corrected in order to allow test coverage
  display on the Bioconductor page

Changes in version 1.3.4:

NEW FEATURES

- vignette: minor correction in the default number of permutations for
  single response (O)PLS(-DA) models

Changes in version 1.3.3:

NEW FEATURES

- vignette: detailed explanations of Q2Y computation

Changes in version 1.3.2:

NEW FEATURES

- plot.opls: 'plotVc' argument renamed 'typeVc' similarly to the
  default 'plot' function

Changes in version 1.3.1:

NEW FEATURES

- opls: now takes either a (numeric) data frame or matrix as 'x' input
  (instead of matrix only)

- predict: now takes either a (numeric) data frame or matrix as
  'newdata' input (instead of matrix only)

ROTS
----

Changes in version 1.0.0:

- Added log fold change in output.

- Added default plotting for ROTS class objects

- Added support for paired tests

- Updated vignette

- Bug fixes

- Modified for Bioconductor release

RPA
---

Changes in version 1.27.41 (2016-04-14):

- Functions for phylogenetic microarray data added

Rqc
---

Changes in version 1.6:

NEW FEATURES

- Support to paired-end reads added

USER VISIBLE CHANGES

- Per file heatmap plot improved

Rsamtools
---------

Changes in version 1.23:

NEW FEATURES

- filterBam can filter one source file into multiple destinations by
  providing a vector of destination files and a list of FilterRules.

- phred2ASCIIOffset() helps translate PHRED encodings (integer or
  character) to ASCII offsets for use in pileup()

BUG FIXES

- scanBam() fails early when param seqlevels not present in file.

- Rsamtools.mk for Windows avoids spaces in file paths

Rsubread
--------

Changes in version 1.22.0:

NEW FEATURES

- featureCounts() is able to count reads of up to 250kb long.

- A parameter `juncCounts` is added to featureCounts() function to
  report counts for exon-exon junctions in RNA-seq data.

- A parameter `nonSplitOnly` is added to featureCounts() function to
  count non-split alignments only.

- Improved parsing of gzipped fastq files in align() and subjunc()
  aligners.

- Improved screen output and error reporting for align(), subjunc() and
  featureCounts().


RUVSeq
------

Changes in version 1.5:

- When data are not integer, only a warning is thrown and not an error.

- Added an example on how to use RUV with DESeq2.

- Added volume and pages to citation.

S4Vectors
---------

Changes in version 0.10.0:

NEW FEATURES

- Add SelfHits class, a subclass of Hits for representing objects where
  the left and right nodes are identical.

- Add utilities isSelfHit() and isRedundantHit() to operate on SelfHits
  objects.

- Add new Pairs class that couples two parallel vectors.

- head() and tail() now work on a DataTable object and behave like on
  an ordinary matrix.

- Add as.matrix.Vector().

- Add "append" methods for Rle/vector (they promote to Rle).

SIGNIFICANT USER-VISIBLE CHANGES

- Many changes to the Hits class: - Replace the old Hits class (where
  the hits had to be sorted by query) with the SortedByQueryHits class.
  - A new Hits class where the hits can be in any order is
  re-introduced as the parent of the SortedByQueryHits class. - The
  Hits() constructor gets the new 'sort.by.query' argument that is
  FALSE by default. When 'sort.by.query' is set to TRUE, the
  constructor returns a SortedByQueryHits instance instead of a Hits
  instance. - Bidirectional coercion is supported between Hits and
  SortedByQueryHits.  When going from Hits to SortedByQueryHits, the
  hits are sorted by query. - Add "c" method for Hits objects. - Rename
  Hits slots: queryHits -> from subjectHits -> to queryLength -> nLnode
  (nb of left nodes) subjectLength -> nRnode (nb of right nodes) - Add
  updateObject() method to update serialized Hits objects from old
  (queryHits/subjectHits) to new (from/to) internal representation. -
  The "show" method for Hits objects now labels columns with from/to by
  default and switches to queryHits/subjectHits labels only when the
  object is a SortedByQueryHits object. - New accessors are provided
  that match the new slot names: from(), to(), nLnode(), nRnode(). The
  old accessors (queryHits(), subjectHits(), queryLength(), and
  subjectLength()) are just aliases for the new accessors. Also
  countQueryHits() and countSubjectHits() are now aliases for new
  countLnodeHits() and countRnodeHits().

- Transposition of Hits objects now propagates the metadata columns.

- Rename elementLengths() -> elementNROWS() (the old name was clearly a
  misnomer). For backward compatibility the old name still works but is
  deprecated (now it's just an "alias" for elementNROWS()).

- Rename compare() -> pcompare(). For backward compatibility the old
  name still works but is just an "alias" for pcompare() and is
  deprecated.

- Some refactoring of the Rle() generic and methods: - Remove ellipsis
  from the argument list of the generic. - Dispatch on 'values' only. -
  The 'values' and 'lengths' arguments now have explicit default values
  logical(0) and integer(0) respectively. - Methods have no more
  'check' argument but new low-level (non-exported) constructor
  new_Rle() does and is what should now be used by code that needs this
  feature.

- Optimize subsetting of an Rle object by an Rle subscript: the
  subscript is no longer decoded (i.e. expanded into an ordinary
  vector). This reduces memory usage and makes the subsetting much
  faster e.g. it can be 100x times faster or more if the subscript has
  many (e.g. thousands) of long runs.

- Modify "replaceROWS" methods so that the replaced elements in 'x' get
  their metadata columns from 'value'. See this thread on bioc-devel:
  https://stat.ethz.ch/pipermail/bioc-devel/2015-November/008319.html

- Remove ellipsis from the argument list of the "head" and "tail"
  methods for Vector objects.

- pc() (parallel combine) now returns a List object only if one of the
  supplied objects is a List object, otherwise it returns an ordinary
  list.

- The "as.data.frame" method for Vector objects now forwards the
  'row.names' argument.

- Export the "parallelSlotNames" methods.

DEPRECATED AND DEFUNCT

- Deprecate elementLengths() in favor of elementNROWS(). New name
  reflects TRUE semantic.

- Deprecate compare() in favor of pcompare().

- After being deprecated in BioC 3.2, the "ifelse" methods for Rle
  objects are now defunct.

- Remove "aggregate" method for vector objects which was an
  undocumented bad idea from the start.

BUG FIXES

- Fix 2 long-standing bugs in "as.data.frame" method for List objects:
  - must always return an ordinary data.frame (was returning a
  DataFrame when 'use.outer.mcols' was TRUE), - when 'x' has names and
  'group_name.as.factor' is TRUE, the levels of the returned group_name
  col must be identical to 'unique(names(x))' (names of empty list
  elements in 'x' was not showing up in 'levels(group_name)').

- Fix and improve the elementMetadata/mcols setter method for Vector
  objects so that the specific methods for GenomicRanges, GAlignments,
  and GAlignmentPairs objects are not needed anymore and were removed.
  Note that this change also fixes setting the elementMetadata/mcols of
  a SummarizedExperiment object with NULL or an ordinary data frame,
  which was broken until now.

- Fix bug in match,ANY,Rle method when supplied 'nomatch' is not NA.

- Fix findMatches() for Rle table.

- Fix show,DataTable-method to display all rows if <= nhead + ntail + 1

scater
------

Changes in version 0.99.3 (2016-02-29):

- Package added to Bioconductor

- Bioc-submission branch merged with master

Changes in version 0.99.2 (2016-02-21):

- Package prepared for Bioconductor submission.

SeqArray
--------

Changes in version 1.11.19-1.11.22:

- utilizes the official C API `R_GetConnection()` to accelerate text
  import and export, requiring R (>=v3.3.0); alternative version
  (backward compatible with R_v2.15.0) is also available on github
  (https://github.com/zhengxwen/SeqArray/releases/tag/v1.11.18)

- ~4x speedup in the sequential version of `seqVCF2GDS()`, and
  `seqVCF2GDS()` can run in parallel

- variables in "annotation/format/" should be two-dimensional as what
  mentioned in the vignette.

Changes in version 1.11.0-1.11.18:

- rewrite `seqSummary()`

- a new vignette file with Rmarkdown format (replacing
  SeqArray-JSM2013.pdf)

- bug fix in `seqBED2GDS()` if the total number of genotypes > 2^31
  (integer overflow)

- bug fixes in `seqMerge()` if chromosome and positions are not unique

- `seqStorage.Option()` is renamed to `seqStorageOption()`

- new function `seqDigest()`

- `seqVCF.Header()` is renamed to `seqVCF_Header()`, `seqVCF.SampID()`
  is renamed to `seqVCF_SampID()`

- seqSetFilter(): 'samp.sel' is deprecated since v1.11.12, please use
  'sample.sel' instead

- accelerate reading genotypes with SSE2(+13%) and AVX2(+23%)

- new function `seqSystem()`

- allow "$dosage" in `seqGetData()` and `seqApply()` for the dosages of
  reference allele

- accelerate `seqSetFilterChrom()` and allow a selection with multiple
  regions

- new methods `\S4method{seqSetFilter}{SeqVarGDSClass, GRanges}()` and
  `\S4method{seqSetFilter}{SeqVarGDSClass, GRangesList}()`

- 'as.is' in `seqApply()` allows a 'connection' object (created by
  file, gzfile, etc)

- `seqSummary(f, "genotype")$seldim` returns a vector with 3 integers
  (ploidy, # of selected samples, # of selected variants) instead of 2
  integers

Changes in version 1.10.0-1.10.6:

- the version number was bumped for the Bioconductor release version
  3.2

- fix a memory issue in `seqAlleleFreq()` when 'ref.allele' is a vector

- `seqSetFilter()` allows numeric vectors in 'samp.sel' and
  'variant.sel'

- `seqSummary()` returns ploidy and reference

- `seqStorage.Option()` controls the compression level of FORMAT/DATA

- `seqVCF2GDS()` allows extract part of VCF files via 'start' and
  'count'

- `seqMerge()` combines multiple GDS files with the same samples

- export methods for compatibility with VariantAnnotation

- a new argument '.useraw' in `seqGetFilter()`

- a new argument 'allow.duplicate' in `seqOpen()`

- fix a bug in `seqParallel()`
  (https://github.com/zhengxwen/SeqArray/issues/11) and optimize its
  performance

- 'gdsfile' could be NULL in `seqParallel()`

SeqVarTools
-----------

Changes in version 1.9.11:

- Add Firth test option to regression

- Bug fix for refFracPlot: hets significantly different from 0.5
  plotted as triangles, median line shown

Changes in version 1.9.10:

- duplicateDiscordance across two GDS files can calculate discordance
  based on heterozygote/homozygote status instead of genotype

Changes in version 1.9.9:

- duplicateDiscordance across two GDS files can match on either
  position or position and alleles, with the ability to recode
  genotypes if th reference allele in one dataset is the alternate
  allele in the other dataset

Changes in version 1.9.8:

- duplicateDiscordance and alternateAlleleDetection require SeqVarData
  objects; both can match on a subject.id instead of sample.id

Changes in version 1.9.4:

- alleleDosage returns list with dosage of each allele separately

Changes in version 1.9.2:

- added by.variant option to duplicateDiscordance for two gds files

sevenbridges
------------

Changes in version 1.1.16:

NEW FEATURES

- Full support for API V2, user-friendly call from R

- CWL Draft 2+ generator in R, create JSON/YAML tool directly

- 5 Vignettes added for comprehensive tutorials and reference

- Three examples inlcuded under inst/docker for cwl app examples

- Auth configuration file to maintain multiple platforms and user
  account

- Works for multiple Seven Bridges supported platforms

- More features like task hook function to ease the automation

Changes in version 1.0.0:

- Initial version

- All the APIs of the SBG platform are supported

- First vignette added

SGSeq
-----

Changes in version 1.6.0:

- New vignette

- Added predictVariantEffects() for predicting the effect of a splice
  variant on annotated protein-coding transcripts

- Changes in the SGVariants and SGVariantCounts class. Instances
  created with previous versions of SGSeq have to be updated.

- Replaced functions for accessing assay data with two generic
  functions counts() and FPKM()

- Support BamFileLists in sample info

- Changed behavior of the annotate() function when assigning gene names

- Changed behavior of the min_denominator argument in analyzeVariants()
  and getSGVariantCounts(). The specfied minimum count now has to be
  achieved at either the start or end of the event.

- Adjacent exons no longer cause a warning in convertToTxFeatures()

- Deprecated legacy classes TxVariants, TxVariantCounts

- Bug fixes and other improvements

SigCheck
--------

Changes in version 2.4:

- Changed main plots to use -log10(p) as X axis

- Added nolegend parameter to sigCheckPlot

- Added title parameter to sigCheckPlot

- Updated nkiResults data object.

- Section added to vignette explaining how to get signatures from
  MSigDB.

sincell
-------

Changes in version 1.3.1:

- fixed bug in IMC algorithm that allowed cells in the same cluster
  connnect between each others instead of keep growing


SNPRelate
---------

Changes in version 1.6.0:

- the version number was bumped for the Bioconductor release version
  3.3

Changes in version 1.5.0-1.5.2:

- fix an issue in `snpgdsVCF2GDS()` if sample.id has white space

- bug fix in `snpgdsPCASampLoading()` when the input is SeqArray GDS
  file

- improve `snpgdsGetGeno()`

specL
-----

Changes in version 1.5.10-13:

USER UNVISIBLE CHANGES

- added specLSet class support for cdsw methode

- changed Rmd5 vignette style

- added cdsw test case

- intro new vignette for cdsw method

Changes in version 1.5.9:

USER UNVISIBLE CHANGES

- added test case for read.biliospec

Changes in version 1.5.5:

USER VISIBLE CHANGES

- added RT prediction vignette file

Changes in version 1.5.4:

USER UNVISIBLE CHANGES

- changed NAMESPACES and read.bibliospec docu to avoid warnings in 3.3
  check

Changes in version 1.5.3:

USER UNVISIBLE CHANGES

- added sqlite files for peptideStd RData

Changes in version 1.5.2:

USER UNVISIBLE CHANGES

- find all signals having two or more in-silico fragment ions.

- keep only the nearest fragment ion; if there are more take the first
  in line

spliceSites
-----------

Changes in version 1.15.0:

- Changed signature for functions dnaRanges and write.annDNA

- Added further explanations for calculations on HBond in vignette BUG
  FIXES


SummarizedExperiment
--------------------

Changes in version 1.2.0:

NEW FEATURES

- Add 'rowData' argument to SummarizedExperiment() constructor. This
  allows the user to supply the row data at construction time.

- The SummarizedExperiment() constructor function and the assay()
  setter now both take any matrix-like object as long as the resulting
  SummarizedExperiment object is valid.

- Support r/cbind'ing of SummarizedExperiment objects with assays of
  arbitrary dimensions (based on a patch by Pete Hickey).

- Add "is.unsorted" method for RangedSummarizedExperiment objects.

- NULL colnames() supported during SummarizedExperiment construction.

- readKallisto() warns early when files need names.

- base::rank() gained a new 'ties.method="last"' option and
  base::order() a new argument ('method') in R 3.3. Thus so do the
  "rank" and "order" methods for RangedSummarizedExperiment objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Re-introduce the rowData() accessor (was defunt in BioC 3.2) as an
  alias for mcols() and make it the preferred way to access the row
  data. There is now a pleasant symmetry between rowData and colData.

- Rename SummarizedExperiment0 class -> SummarizedExperiment.

- Improved vignette.

- Remove updateObject() method for "old" SummarizedExperiment objects.

DEPRECATED AND DEFUNCT

- exptData() is now defunct, metadata() should be used instead.

BUG FIXES

- Fix bug in "sort" method for RangedSummarizedExperiment objects when
  'ignore.strand=TRUE' (the argument was ignored).

- Fix 2 bugs when r/cbind'ing SummarizedExperiment objects: -
  r/cbind'ing assays without names would return only the first element.
  See
  https://stat.ethz.ch/pipermail/bioc-devel/2015-November/008318.html -
  r/cbind'ing assays with names in different order would stop() with
  'Assays must have the same names()"

- Fix validity method for SummarizedExperiment objects reporting
  incorrect numbers when the nb of cols in assay(x) doesn't match the
  nb of rows in colData(x).

- assay colnames() must agree with colData rownames()

- Fix bug where assays(se, withDimnames=TRUE) was dropping the dimnames
  of the 3rd and higher-order dimensions of the assays. Thanks to Pete
  Hickey for catching this and providing a patch.

- A couple of minor tweaks to the rowData() setter to make it behave
  consistently with mcols()/elementMetadata() setters for Vector
  objects in general.

SWATH2stats
-----------

Changes in version 1.1.17:

DOCUMENTATION

- updated NEWS file

Changes in version 1.1.16:

BUG FIXES

- Error fix in test_convert.R

Changes in version 1.1.15:

DOCUMENTATION

- Add Citation of PLoS ONE publication of SWATH2stats

Changes in version 1.1.14:

DOCUMENTATION

- import_data: minor changes in documentation

- defining several functions at the beginning of function script

Changes in version 1.1.13:

DOCUMENTATION

- assess_fdr_byrun, assess_fdr_overall, filter_mscore_fdr: add sentence
  to manual page about FFT

- vignettes: Add/change title

BUG FIXES

- DESCRIPTION: add knitr to Suggests and VignetteBuilder, add R>=2.10.0
  to Depends

Changes in version 1.1.12:

NEW FEATURES

- sample_annotation: remove option column.runid.

- assess_fdr_byrun and assess_fdr_overall: add option to set range of
  plotting with n.range

DOCUMENTATION

- SWATH2stats_example_script, Spyogenes data: added Example script and
  S.pyogenes data

Changes in version 1.1.11:

NEW FEATURES

- plot_correlation_between_samples: import cor function from stats
  package and reorder y axis.

Changes in version 1.1.10:

BUG FIXES

- test for FDR filtering: bug fix due to fix in mscore4assayfdr

Changes in version 1.1.9:

DOCUMENTATION

- Update vignette

DEPRECATED AND DEFUNCT

- remove function: filter_all_peptides

BUG FIXES

- mscore4assayfdr: bug fix

Changes in version 1.1.8:

BUG FIXES

- add functions that were not updated

Changes in version 1.1.7:

NEW FEATURES

- added functions: count_analytes, plot_correlation_between_samples,
  plot_variation, plot_variation_vs_total, transform_MSstats_OpenSWATH.

- filter_mscore_requant renamed to filter_mscore_freqobs

BUG FIXES

- sample_annotation: Bug fix - it reported upon error different
  conditions instead of filenames.

- assess_fdr_overall: Correction of transition level column name from
  "id" to "transition_group_id"

- assess_fdr_byrun: Correction of transition level column name from
  "id" to "transition_group_id"

- plot.fdr_cube: added na.rm=TRUE to plotting functions

Changes in version 1.1.6:

BUG FIXES

- Correction of typographical error in Vignette for Wolski et al.

Changes in version 1.1.5:

BUG FIXES

- Bug fixes for assess_fdr_byrun that was also introduced into version
  v1.0.2 at the same time. It fixes problems with labelling the runs
  correctly in some cases.

Changes in version 1.1.4:

BUG FIXES

- Bug fixes in test_convert.

Changes in version 1.1.1:

NEW FEATURES

- Improved the function disaggregate() that also data with different
  number of transitions per precursor than 6 can be used

- Added tests for disaggregate.R and convert4pythonscript.R

Changes in version 1.1.0:

NEW FEATURES

- Development version of SWATH2stats in BioC 3.3

synapter
--------

Changes in version 1.13.1:

- Update call to nQuants to accomodate changes in MSnbase

- Defunct synapterGUI <2016-02-29 Mon>

Changes in version 1.13.0:

- Bioc devel 3.2

TarSeqQC
--------

Changes in version 1.1.6:

CODE

- Changes in summaryIntervals in order to allow the exploration at pool
  levels, useful when a targeted sequencing involving several PCR pools
  was performed.

- plotAttrPerform method was added. This function produces a ggplot
  graph illustrating relative and cumulative frequencies of features in
  attribute intervals. If the panel has several pools, then the graph
  shows the mentioned results for each pool.

Changes in version 1.1.2:

CODE

- Changes in buildFeaturePanel in order to reduce the run time. Now,
  the pileupCounts is not called, thus the pileup matrix is not built.
  Instead it, coverage and others Rsamtools and IRanges methods are
  used.

- readPercentages and plotInOutFeatures methods were added in order to
  explore experiment efficiency.

- biasExploration and plotMetaDataExpl were added in order to perform
  bias sources exploration. The first allows the gc content, feature
  length or other source distribution exploration. The second implement
  a plot in which attribute distribution for each source bias quartile
  or group is explored.

VIGNETTE

- English correction and new methods incorporation

Changes in version 1.1.1:

CODE

- In setters and initialize functions, examples were updated to shorten
  elapsed times (only interactive). Thanks to V. Obenchain

VIGNETTE

- Elapsed times reduction using TarSeqQC data example instead of
  building myPanel from scratch. Thanks to V. Obenchain

Changes in version 1.1.0:

DOCUMENTATION

- Bumped version number of all packages after creation of 3.2 branch

TCGAbiolinks
------------

Changes in version 1.1.26:

- Bug fix: TCGAprepare for IlluminaHiSeq_miRNASeq platform was not
  considering all file types

Changes in version 1.1.25:

- Bug fix: TCGAquery_maf had a bug if only one table was found.

Changes in version 1.1.24:

- Bug fix: types argument for Genome_Wide_SNP_6 was not working.

Changes in version 1.1.23:

- Bug fix: logFC result from TCGAanalyze_DEA for method "glmLRT" was
  considerent the alphabetical order of groups, which migth induce the
  user to error

Changes in version 1.1.22:

- Bug fix: when checking the data integrity, manifest was being read
  with argument header = T,but it should be false.

Changes in version 1.1.21:

- Improvement: TCGADownload checks for data integrity.

Changes in version 1.1.20:

- Bug fix: TCGAanalyze_Preprocessing subseting was incorrect in case
  there were outliers.

Changes in version 1.1.19:

- Bug fix: preparing the data for "IlluminaGA_DNASeq" platform and
  TCGAquery_maf had a problem for some files.

Changes in version 1.1.18:

- Update lgg and gbm subtype information. Source:
  http://dx.doi.org/10.1016/j.cell.2015.12.028

Changes in version 1.1.17:

- Adding batch information to the package. batch.info object is
  available for user and TCGAprepare adds automtically info for the
  summarizedExperiment object

- TCGAvisualize_starburst new parameter: circle, to draw or not the
  circles in the plot

- Database update

Changes in version 1.1.16:

- Bug fix: subsetByOverlaps was removed from SummarizedExperiment
  package TCGAbiolinks should not import it

- Small fixes in documentation

- TCGAanalyze_DMR is now saving the results in a cvs file

Changes in version 1.1.15:

- ggplot2 updat broke the package. Some small fixes were made, but the
  function TCGAvisualize_profilePlot is not working as sjPlot is not
  updated yet.

- small fixes in documentation

Changes in version 1.1.14:

- Change in TCGAvisualize_Heatmap: coloring the columns by patient
  might give wrong results if patient has more than one sample in the
  hetamap. To solve that we added the possibility to use the sample
  columns and we added a warning if there is more than one sample for
  the patients and a patient column is being used.

Changes in version 1.1.11:

- Update in the citation

Changes in version 1.1.10:

- TCGAPrepare: bug fix for bt.exon_quantification files from
  IlluminaHiSeq_RNASeqV2 platform

- Database update TCGAbiolinks 1.1.8

- TCGAvisualize_Heatmap Now it is using Heatmap plus package and is
  calculating z-cores

- TCGAvisualize_profilePlot Visualize the distribution of subgroups in
  the groups distributions

- Database update

- From version 1.0: small bugs corrections in some plots and
  TCGAprepare_elmer, documentation improvement. TCGAbiolinks 0.99.2
  FIRST VERSION - FEATURES

- TCGAanalyze_DEA Differentially expression analysis (DEA) using edgeR
  package.

- TCGAanalyze_DMR Differentially methylated regions Analysis

- TCGAanalyze_EA Enrichment analysis of a gene-set with GO [BP,MF,CC]
  and pathways.

- TCGAanalyze_EAcomplete Enrichment analysis for Gene Ontology (GO)
  [BP,MF,CC] and Pathways

- TCGAanalyze_Filtering Filtering mRNA transcripts and miRNA selecting
  a threshold.

- TCGAanalyze_LevelTab Adding information related to DEGs genes from
  DEA as mean values in two conditions.

- TCGAanalyze_Normalization normalization mRNA transcripts and miRNA
  using EDASeq package.

- TCGAanalyze_Preprocessing Array Array Intensity correlation (AAIC)
  and correlation boxplot to define outlier

- TCGAanalyze_survival Creates survival analysis

- TCGAanalyze_SurvivalKM survival analysis (SA) univariate with
  Kaplan-Meier (KM) method.

- TCGAbiolinks Download data of samples from * TCGA

- TCGAdownload Download the data from * TCGA using as reference the
  output from * TCGAquery

- TCGAintegrate Filtering common samples among platforms from *
  TCGAquery for the same tumor

- TCGAinvestigate Find most studied TF in pubmed related to a specific
  cancer, disease, or tissue

- TCGAprepare Read the data from level 3 the experiments and prepare it
  for downstream analysis into a SummarizedExperiment object.

- TCGAquery Searches * TCGA open-access data providing also latest
  version of the files.

- TCGAquery_clinic Get the clinical information

- TCGAquery_clinicFilt Filter samples using clinical data

- TCGAquery_MatchedCoupledSampleTypes Retrieve multiple tissue types
  from the same patients.

- TCGAquery_samplesfilter Filtering sample output from * TCGAquery

- TCGAquery_SampleTypes Retrieve multiple tissue types not from the
  same patients.

- TCGAquery_Version Shows a summary (version, date, number of samples,
  size of the data) of all versions of data for a given tumor and
  platform.

- TCGAsocial Finds the number of downloads of a package on CRAN or BIOC
  and find questions in website ("bioconductor.org", "biostars.org",
  "stackoverflow).

- TCGAvisualize_EAbarplot barPlot for a complete Enrichment Analysis

- TCGAvisualize_meanMethylation Mean methylation boxplot

- TCGAvisualize_PCA Principal components analysis (PCA) plot

- TCGAvisualize_starburst Create starburst plot

- TCGAvisualize_SurvivalCoxNET Survival analysis with univariate Cox
  regression package (dnet)

TEQC
----

Changes in version 3.11.1:

- bug fix in 'reads2pairs': when there are chromosomes/contigs without
  any mapped reads, those will be removed automatically

TFBSTools
---------

Changes in version 3.3:

BUG FIXES

- Adapt the runMEME to work with meme 4.10.x version.

- Fix the scientific notation in run_MEME

- Better error handling of MEME wrappe

Changes in version 1.9.4:

NEW FEATURES

- Add conversion from IUPAC string to matrix

- toPWM, toICM work on PFMatrixList

BUG FIXES

- Fix the seqLogo error when the frequency of each base is same at
  certain site. Thanks to Liz.

- Fix the database interface to deal with the duplicated/missing
  records for certain TFBS in JASPAR2016.

- Fix coersion method failure on certain cases.

tigre
-----

Changes in version 1.24.2:

BUG FIXES

- R-3.2.4 compatibility update

Changes in version 1.24.1:

BUG FIXES

- Fix target ranking when data have no variances


trackViewer
-----------

Changes in version 1.7.9:

- fix the bug that legend write to outside.

Changes in version 1.7.8:

- fix the typo in documentation.

- update the lollipop plot.

Changes in version 1.7.7:

- add dandelion plot function.

Changes in version 1.7.6:

- add legend for lollipop plot.

Changes in version 1.7.5:

- add lables in the nodes of lollipop plot.

Changes in version 1.7.4:

- Update documentation.

Changes in version 1.7.3:

- add Y axis to lollipop plot.

Changes in version 1.7.2:

- add lollipop plot!

Changes in version 1.7.1:

- adjust the fontsize for optimizing styles with theme.

- add gene symbols if possible for geneModelFromTxdb.

TransView
---------

Changes in version 1.15.1:

BUG FIXES

- Importing from S4Vectors instead of IRanges

TRONCO
------

Changes in version 2.4.0:

- New statistics available for model confidence via cross-validation
  routines. New algorithms based on Minimum Spanning Tree extraction.

tximport
--------

Changes in version 0.99.0:

- Preparing package for Bioconductor submission.

Changes in version 0.0.19:

- Added `summarizeToGene` which breaks out the gene-level summary step,
  so it can be run by users on lists of transcript-level matrices
  produced by `tximport` with `txOut=TRUE`.

Changes in version 0.0.18:

- Changed argument `gene2tx` to `tx2gene`. This order is more
  intuitive: linking transcripts to genes, and matches the `geneMap`
  argument of Salmon and Sailfish.

Uniquorn
--------

Changes in version 0.99.15 (2016-05-03):

Minor Update of statistics

- Minor fixes to statistics

Changes in version 0.99.14 (2016-05-03):

Major Update of statistics

- Replaced the absolute and relative cutoff with a p an q-value based
  on a binomial distribution

Changes in version 0.99.13 (2016-04-08):

Fixed unit tests

- Further update test suite

Changes in version 0.99.12 (2016-04-08):

Fixed unit tests

- Update testing suite

Changes in version 0.99.11 (2016-04-01):

Fixed unit tests

- Fixed error in unit tests

Changes in version 0.99.10 (2016-03-30):

Added Unit tests

- Added Unit test feature

Changes in version 0.99.9 (2016-03-29):

Improved documentation

- Minor documentation bugfixes.

Changes in version 0.99.8 (2016-03-25):

Fixed Bug

- General bugfixing.

Changes in version 0.99.7 (2016-03-24):

Fixed Bug

- Minor bugfix concerning output of CL identifications.

Changes in version 0.99.6 (2016-03-23):

Fixed Bug

- Minor bugfix concerning data parsing of CellMiner data.

Changes in version 0.99.5 (2016-03-22):

Fixed Bug

- Fixed a bug that lead to incorrect calculation of CoSMIC CLP CL
  mutation weights.

Changes in version 0.99.4 (2016-03-16):

Added features

- Minimum amount of matching mutations required for a positive
  prediction between query and training sample can now be manually
  adjusted in the 'identify_vcf_file' function using the parameter
  'miminum_matching_mutations'

Changes in version 0.99.3 (2016-03-02):

Fixed bugs

- Fixed vigniette

Changes in version 0.99.2 (2016-03-02):

Fixed bugs

- Moved NEWS.Rd from 'Uniquorn' to 'Uniquorn/inst'

Changes in version 0.99.1 (2016-03-02):

Fixed bugs

- Fixed NEWS section

- Fixed unit testing

- Fixed the Vignette: replaced the github-dependent installation with
  the biocLite installation

Changes in version 0.99.0 (2015-01-27):

New features

- Identify cancer cell lines

- Show which cancer cell line are contained

- Show which mutations are annotated for selected cancer cell lines

- Show which mutaitons are overall included

- Parse cancer cell line custom data -> add your own samples and
  identify these

variancePartition
-----------------

Changes in version 1.1.7:

- GPL License

Changes in version 1.1.6:

- Move packages from Depends to Imports

- For clarity, replace = with <- in parts of examples and vignette

- Stop cluster in examples to solve error on Windows machines

Changes in version 1.1.5:

- Stop cluster in vignette to solve error on Windows machines

Changes in version 1.1.4:

- Fix Bioconductor check error

Changes in version 1.1.3:

- Add details to vignette

- Fix ggplot2 compatibility issues

Changes in version 1.1.2:

- Add details to vignette

Changes in version 1.1.1:

- add plotPercentBars() to vizualize variance fractions for a subset of
  genes

- add ESS() to compute effective sample size

- fix x.labels argument in plotStratifyBy().  Previously, this argument
  was not used correctly

VariantAnnotation
-----------------

Changes in version 1.20.0:

NEW FEATURES

- add SnpMatrixToVCF()

- add patch from Stephanie Gogarten to support 'PL' in
  genotypeToSnpMatrix()

MODIFICATIONS

- move getSeq,XStringSet-method from VariantAnnotation to BSgenome

- update filterVcf vignette

- remove 'pivot' export

- work on readVcf(): - 5X speedup for readVcf (at least in one case) by
  not using "==" to compare a list to a character (the list gets
  coerced to character, which is expensive for huge VCFs) - avoiding
  relist.list()

- update summarizeVariants() to comply with new SummarizedExperiment
  rownames requirement

- defunct VRangesScanVcfParam() and restrictToSNV()

- use elementNROWS() instead of elementLengths()

- togroup(x) now only works on a ManyToOneGrouping object so replace
  togroup(x, ...) calls with togroup(PartitioningByWidth(x), ...) when
  'x' is a list-like object that is not a ManyToOneGrouping object.

- drop validity assertion that altDepth must be NA when alt is NA there
  are VCFs in the wild that use e.g. "\*" for alt, but include depth

- export PLtoGP()

- VariantAnnotation 100% RangedData-free

BUG FIXES

- use short path names in src/Makevars.win

Changes in version 1.18.0:

MODIFICATIONS

- defunct VRangesScanVcfParam()

- defunct restrictToSNV()

BUG FIXES

- scanVcf,character,missing-method ignores blank data lines.

- Build path for C code made robust on Windows.

VariantFiltering
----------------

Changes in version 1.8:

USER VISIBLE CHANGES

- Changed human defaults of VariantFilteringParam() and updated
  documentation on MafDb packages to reflect newer package (shorter)
  names of the 1000 Genomes Project.

- Sequence Ontology (SO) annotations have been updated from
  https://github.com/The-Sequence-Ontology/SO-Ontologies to the latest
  version of the data available on April, 2016.

- Analysis methods (autosomalDominant(), etc.) now accept an argument
  called 'svparam' that takes an object produced by the
  'ScanVcfParam()' function from the VariantAnnotation package. This
  allows one to parametrize the way in which VCF files are read. For
  instance, if one wishes to analyze a specific set of genomic ranges.

- Analysis methods (autosomalDominant(), etc.) now accept an argument
  called 'use' that allows the user to select among three simple
  strategies to handle missing genotypes. See the corresponding help
  page for further information.

- No messages from the AnnotationDbi::select() method are given anymore
  about the 1:1, 1:many or many:1 results obtained when fetching
  annotations.

BUG FIXES

- Several bugfixes including dealing with transcript-centric
  annotations anchored at Ensembl gene identifiers, avoid querying for
  an OMIM column when working with organisms other than human,
  correctly identifying unaffected individuals in the autosomal
  recessive heterozygous analysis.

- Added methods to deal with the new ExAC MafDb packages that enable
  the user to query allele frequencies by position.

wavClusteR
----------

Changes in version 2.5.0:

- fixed bug in plotSubstitutions. The transition type of interest was
  not correctly highlighted in some instances. Thanks to Charlotte
  Sonenson for pointing this out.

- Vignette migrated to Rmarkdown.

xcms
----

Changes in version 1.47.3:

- Disable parallel processing in unit tests causing a timeout on BioC
  build machines

Changes in version 1.47.2:

BUG FIXES

- Fix problem in getEIC on xcmsSet objects reported by Alan Smith in
  issue #7 and add a RUnit test case to test for this (test.issue7 in
  runit.getEIC.R).

- Changed some unnecessary warnings into messages.

USER VISIBLE CHANGES

- Disabled parallel processing in unit tests * migrate dependencies
  from ncdf -> ncdf4



Packages removed since the last release
=================================

No packages were removed from the release.

17 packages were marked as deprecated, to be removed in the next release.

One package, sbgr, was renamed to sevenbridges.

Deprecated packages:

* AffyTiling
* caFlowQ
* cellHTS
* DASiR
* DAVIDQuery
* GenoView
* inSilicoDb
* inSilicoMerging
* jmosaics
* metaX
* MMDiff
* neaGUI
* Rolexa
* RWebServices
* SJava
* SomatiCA
* spade
