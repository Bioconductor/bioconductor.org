Bioconductors:

We are pleased to announce Bioconductor 3.1, consisting of 1024
software packages, 241 experiment data packages, and 917
up-to-date annotation packages. 

There are 95 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.1 is compatible with R 3.2,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release includes an updated Bioconductor [Amazon Machine Image][1]
and [Docker containers][2].

Visit [http://bioconductor.org][3]
for details and downloads.

[1]: /help/bioconductor-cloud-ami/
[2]: /help/docker/
[3]: http://bioconductor.org

Contents
--------

* Getting Started with Bioconductor 3.1
* New Software Packages
* NEWS from new and existing packages
* Packages removed from Bioconductor since the last release.

Getting Started with Bioconductor 3.1
======================================

To update to or install Bioconductor 3.1:

1. Install R 3.2.  Bioconductor 3.1 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 95 new packages in this release of Bioconductor.

AIMS - This package contains the AIMS implementation. It
contains necessary functions to assign the five intrinsic
molecular subtypes (Luminal A, Luminal B, Her2-enriched,
Basal-like, Normal-like). Assignments could be done on
individual samples as well as on dataset of gene expression
data.

AnalysisPageServer - AnalysisPageServer is a modular system
that enables sharing of customizable R analyses via the
web.

bamsignals - This package allows to efficiently obtain
count vectors from indexed bam files. It counts the number
of reads in given genomic ranges and it computes reads
profiles and coverage profiles. It also handles paired-end
data.

BEclear - Provides some functions to detect and correct for
batch effects in DNA methylation data. The core function
"BEclear" is based on latent factor models and can also be
used to predict missing values in any other matrix
containing real numbers.

birte - Expression levels of mRNA molecules are regulated
by different processes, comprising inhibition or activation
by transcription factors and post-transcriptional
degradation by microRNAs. biRte uses regulatory networks of
TFs, miRNAs and possibly other factors, together with mRNA,
miRNA and other available expression data to predict the
relative influence of a regulator on the expression of its
target genes. Inference is done in a Bayesian modeling
framework using Markov-Chain-Monte-Carlo. A special feature
is the possibility for follow-up network reverse
engineering between active regulators.

BrowserViz - Interactvive graphics in a web browser from R,
using websockets and JSON

BrowserVizDemo - A BrowserViz subclassing example, xy
plotting in the browser using d3

BubbleTree - BubbleTree utilizes homogenous pertinent
somatic copy number alterations (SCNAs) as markers of tumor
clones to extract estimates of tumor ploidy, purity and
clonality.

canceR - The package is user friendly interface based on
the cgdsr and other modeling packages to explore, compare,
and analyse all available Cancer Data (Clinical data, Gene
Mutation, Gene Methylation, Gene Expression, Protein
Phosphorylation, Copy Number Alteration) hosted by the
Computational Biology Center at Memorial-Sloan-Kettering
Cancer Center (MSKCC).

CAnD - Functions to perform the non-parametric and
parametric CAnD tests on a set of ancestry proportions. For
a particular ancestral subpopulation, a user will supply
the estimated ancestry proportion for each sample, and each
chromosome or chromosomal segment of interest. A p-value
for each chromosome as well as an overall CAnD p-value will
be returned for each test. Plotting functions are also
available.

Cardinal - Implements statistical & computational tools for
analyzing mass spectrometry imaging datasets, including
methods for efficient pre-processing, spatial segmentation,
and classification.

chromDraw - Package chromDraw is a simple package for
linear and circular type of karyotype visualization. The
linear type of visualization is usually used for
demonstrating chromosomes structures in karyotype and the
circular type of visualization is used for comparing of
karyotypes between each other. This tool has own input data
format or genomicRanges structure can be used as input.
Each chromosome containing definition of blocks and
centromere position.  Output file formats are *.eps and
*.svg.

Clonality - Statistical tests for clonality versus
independence of tumors from the same patient based on their
LOH or genomewide copy number profiles

CODEX - A normalization and copy number variation calling
procedure for whole exome DNA sequencing data. CODEX relies
on the availability of multiple samples processed using the
same sequencing pipeline for normalization, and does not
require matched controls. The normalization model in CODEX
includes terms that specifically remove biases due to GC
content, exon length and targeting and amplification
efficiency, and latent systemic artifacts. CODEX also
includes a Poisson likelihood-based recursive segmentation
procedure that explicitly models the count-based exome
sequencing data.

cogena - Description: Gene set enrichment analysis is a
valuable tool for the study of molecular mechanisms that
underpin complex biological traits. As the method is
conventionally used on entire omic datasets, such as
transcriptomes, it may be dominated by pathways and
processes that are substantially represented in a dataset,
however the approach may overlook smaller scale, but highly
correlated cellular events that may be of great biological
relevance. In order to detect these discrete molecular
triggers, we developed a tool, co-expressed gene-set
enrichment analysis (cogena), for clustering differentially
expressed genes and identification of highly correlated
molecular expression clusters. Cogena offers the user a
range of clustering methods, including hierarchical
clustering, model based clustering and self-organised
mapping, based on different distance metrics like
correlation and mutual information. After obtaining and
visualising clusters, cogena performs gene set enrichment.
These gene sets can be sourced from the Molecular
Signatures Database (MSigDB) or user-defined gene sets. By
performing gene set enrichment across expression clusters,
we find considerable enhancement in the resolution of
molecular signatures in omic data at the cluster level
compared to the whole.

coMET - Visualisation of EWAS results in a genomic region.
In addition to phenotype-association P-values, coMET also
generates plots of co-methylation patterns and provides a
series of annotation tracks. It can be used to other
omic-wide association scans as long as the data can be
translated to genomic level and for any species.

ComplexHeatmap - Complex heatmaps are efficient to
visualize associations between different sources of data
sets and reveal potential features. Here the ComplexHeatmap
package provides a highly flexible way to arrange multiple
heatmaps and supports self-defined annotation graphics.

conumee - This package contains a set of processing and
plotting methods for performing copy-number variation (CNV)
analysis using Illumina 450k methylation arrays.

CopywriteR - CopywriteR generate DNA copy number profiles
from whole exome sequencing by analysing the offtarget
sequence reads. By exploiting the offtarget sequence reads,
it allows for creation of robust copy number profiles from
WES with a uniform read depth and evenly distributed data
points over the genome.

cpvSNP - Gene set analysis methods exist to combine
SNP-level association p-values into gene sets, calculating
a single association p-value for each gene set. This
package implements two such methods that require only the
calculated SNP p-values, the gene set(s) of interest, and a
correlation matrix (if desired). One method (GLOSSI)
requires independent SNPs and the other (VEGAS) can take
into account correlation (LD) among the SNPs. Built-in
plotting functions are available to help users visualize
results.

cytofkit - An integrated mass cytometry data analysis
pipeline that enables simultaneous illustration of cellular
diversity and progression.

diffHic - Detects differential interactions across
biological conditions in a Hi-C experiment. Methods are
provided for read alignment and data pre-processing into
interaction counts. Statistical analysis is based on edgeR
and supports normalization and filtering. Several
visualization options are also available.

diggit - Inference of Genetic Variants Driving Cellullar
Phenotypes by the DIGGIT algorithm

DMRcaller - Uses Bisulfite sequencing data in two
conditions and identifies differentially methylated regions
between the conditions in CG and non-CG context. The input
is the CX report files produced by Bismark and the output
is a list of DMRs stored as GRanges objects.

edge - The edge package implements methods for carrying out
differential expression analyses of genome-wide gene
expression studies. Significance testing using the optimal
discovery procedure and generalized likelihood ratio tests
(equivalent to F-tests and t-tests) are implemented for
general study designs. Special functions are available to
facilitate the analysis of common study designs, including
time course experiments. Other packages such as snm, sva,
and qvalue are integrated in edge to provide a wide range
of tools for gene expression analysis.

EMDomics - The EMDomics algorithm is used to perform a
supervised two-class analysis to measure the magnitude and
statistical significance of observed continuous genomics
data between two groups. Usually the data will be gene
expression values from array-based or sequence-based
experiments, but data from other types of experiments can
also be analyzed (e.g. copy number variation). Traditional
methods like Significance Analysis of Microarrays (SAM) and
Linear Models for Microarray Data (LIMMA) use significance
tests based on summary statistics (mean and standard
deviation) of the two distributions. This approach lacks
power to identify expression differences between groups
that show high levels of intra-group heterogeneity. The
Earth Mover's Distance (EMD) algorithm instead computes the
"work" needed to transform one distribution into the other,
thus providing a metric of the overall difference in shape
between two distributions. Permutation of sample labels is
used to generate q-values for the observed EMD scores.

ENCODExplorer - This package allows user to quickly access
ENCODE project files metadata and give access to helper
functions to query the ENCODE rest api, download ENCODE
datasets and save the database in SQLite format.

ENmix - Illumina HumanMethylation450 BeadChip has a complex
array design, and the measurement is subject to
experimental variations. The ENmix R package provides tools
for low level data preprocessing to improve data quality.
It incorporates a model based background correction method
ENmix, and functions for inter-array quantile
normalization, data quality checking, exploration of
multimodally distributed CpGs and source of data variation.
To support large scale data analysis, the package also
provides multi-processor parallel computing wrappers for
some commonly used data preprocessing methods, such as BMIQ
probe design type bias correction and ComBat batch effect
correction.

ensembldb - The package provides functions to create and
use transcript centric annotation databases/packages. The
annotation for the databases are directly fetched from
Ensembl using their Perl API. The functionality and data is
similar to that of the TxDb packages from the
GenomicFeatures package, but, in addition to retrieve all
gene/transcript models and annotations from the database,
the ensembldb package provides also a filter framework
allowing to retrieve annotations for specific entries like
genes encoded on a chromosome region or transcript models
of lincRNA genes.

FISHalyseR - FISHalyseR provides functionality to process
and analyse digital cell culture images, in particular to
quantify FISH probes within nuclei. Furthermore, it extract
the spatial location of each nucleus as well as each probe
enabling spatial co-localisation analysis.

FlowRepositoryR - This package provides an interface to
search and download data and annotations from
FlowRepository (flowrepository.org). It uses the
FlowRepository programming interface to communicate with a
FlowRepository server.

FlowSOM - FlowSOM offers visualization options for
cytometry data, by using Self-Organizing Map clustering and
Minimal Spanning Trees

flowVS - Per-channel variance stabilization from a
collection of flow cytometry samples by Bertlett test for
homogeneity of variances. The approach is applicable to
microarrays data as well.

gdsfmt - This package provides a high-level R interface to
CoreArray Genomic Data Structure (GDS) data files, which
are portable across platforms and include hierarchical
structure to store multiple scalable array-oriented data
sets with metadata information. It is suited for
large-scale datasets, especially for data which are much
larger than the available random-access memory. The gdsfmt
package offers the efficient operations specifically
designed for integers with less than 8 bits, since a single
genetic/genomic variant, like single-nucleotide
polymorphism (SNP), usually occupies fewer bits than a
byte. Data compression and decompression are also supported
with relatively efficient random access. It is allowed to
read a GDS file in parallel with multiple R processes
supported by the package parallel.

GENESIS - The GENESIS package provides methodology for
estimating, inferring, and accounting for population and
pedigree structure in genetic analyses.  The current
implementation provides functions to perform PC-AiR
(Conomos et al., 2015): a Principal Components Analysis
with genome-wide SNP genotype data for robust population
structure inference in samples with related individuals
(known or cryptic).

genomation - A package for summary and annotation of
genomic intervals. Users can visualize and quantify genomic
intervals over pre-defined functional regions, such as
promoters, exons, introns, etc. The genomic intervals
represent regions with a defined chromosome position, which
may be associated with a score, such as aligned reads from
HT-seq experiments, TF binding sites, methylation scores,
etc. The package can use any tabular genomic feature data
as long as it has minimal information on the locations of
genomic intervals. In addition, It can use BAM or BigWig
files as input.

gespeR - Estimates gene-specific phenotypes from off-target
confounded RNAi screens. The phenotype of each siRNA is
modeled based on on-targeted and off-targeted genes, using
a regularized linear regression model.

ggtree - ggtree extends the ggplot2 plotting system which
implemented the grammar of graphics. ggtree is designed for
visualizing phylogenetic tree and different types of
associated annotation data.

GoogleGenomics - Provides an R package to interact with the
Google Genomics API.

gQTLBase - Infrastructure for eQTL, mQTL and similar
studies.

gQTLstats - computationally efficient analysis of eQTL,
mQTL, dsQTL, etc.

GreyListChIP - Identify regions of ChIP experiments with
high signal in the input, that lead to spurious peaks
during peak calling. Remove reads aligning to these regions
prior to peak calling, for cleaner ChIP analysis.

gtrellis - Genome level Trellis graph visualizes genomic
data conditioned by genomic categories (e.g. chromosomes).
For each genomic category, multiple dimensional data which
are represented as tracks describe different features from
different aspects. This package provides high flexibility
to arrange genomic categories and add self-defined graphics
in the plot.

HIBAG - It is a software package for imputing HLA types
using SNP data, and relies on a training set of HLA and SNP
genotypes. HIBAG can be used by researchers with published
parameter estimates instead of requiring access to large
training sample datasets. It combines the concepts of
attribute bagging, an ensemble classifier method, with
haplotype inference for SNPs and HLA types. Attribute
bagging is a technique which improves the accuracy and
stability of classifier ensembles using bootstrap
aggregating and random variable selection.

immunoClust - Model based clustering and meta-clustering of
Flow Cytometry Data

InPAS - Alternative polyadenylation (APA) is one of the
important post-transcriptional regulation mechanism which
occurs in most human genes. InPAS, developed form DaPars
algorithm, predicts and estimates APA and cleavage sites
for mRNA-seq data. It uses the power of cleanUpdTSeq to
adjust cleavage sites.

IVAS - Identification of genetic variants affecting
alternative splicing.

LEA - LEA is an R package dedicated to landscape genomics
and ecological association tests. LEA can run analyses of
population structure and genome scans for local adaptation.
It includes statistical methods for estimating ancestry
coefficients from large genotypic matrices and evaluating
the number of ancestral populations (snmf, pca); and
identifying genetic polymorphisms that exhibit high
correlation with some environmental gradient or with the
variables used as proxies for ecological pressures (lfmm),
and controlling the false discovery rate. LEA is mainly
based on optimized C programs that can scale with the
dimension of very large data sets.

LowMACA - The LowMACA package is a simple suite of tools to
investigate and analyze the mutation profile of several
proteins or pfam domains via consensus alignment. You can
conduct an hypothesis driven exploratory analysis using our
package simply providing a set of genes or pfam domains of
your interest.

mAPKL - We propose a hybrid FS method (mAP-KL), which
combines multiple hypothesis testing and affinity
propagation (AP)-clustering algorithm along with the
Krzanowski & Lai cluster quality index, to select a small
yet informative subset of genes.

MatrixRider - Calculates a single number for a whole
sequence that reflects the propensity of a DNA binding
protein to interact with it. The DNA binding protein has to
be described with a PFM matrix, for example gotten from
Jaspar.

mdgsa - Functions to preform a Gene Set Analysis in several
genomic dimensions. Including methods for miRNAs.

MeSHSim - Provide for measuring semantic similarity over
MeSH headings and MEDLINE documents

MethTargetedNGS - Perform step by step methylation analysis
of Next Generation Sequencing data.

mogsa - This package provide a method for doing gene set
analysis based on multiple omics data.

msa - This package provides a unified R/Bioconductor
interface to the multiple sequence alignment algorithms
ClustalW, ClustalOmega, and Muscle. All three algorithms
are integrated in the package, therefore, they do not
depend on any external software tools and are available for
all major platforms. The multiple sequence alignment
algorithms are complemented by a function for
pretty-printing multiple sequence alignments using the
LaTeX package TeXshade.

muscle - MUSCLE performs multiple sequence alignments of
nucleotide or amino acid sequences.

NanoStringQCPro - NanoStringQCPro provides a set of quality
metrics that can be used to assess the quality of
NanoString mRNA gene expression data -- i.e. to identify
outlier probes and outlier samples. It also provides
different background subtraction and normalization
approaches for this data. It outputs suggestions for
flagging samples/probes and an easily sharable html quality
control output.

netbenchmark - This package implements a benchmarking of
several gene network inference algorithms from gene
expression data.

nethet - Package nethet is an implementation of statistical
solid methodology enabling the analysis of network
heterogeneity from high-dimensional data. It combines
several implementations of recent statistical innovations
useful for estimation and comparison of networks in a
heterogeneous, high-dimensional setting. In particular, we
provide code for formal two-sample testing in Gaussian
graphical models (differential network and GGM-GSA; Stadler
and Mukherjee, 2013, 2014) and make a novel network-based
clustering algorithm available (mixed graphical lasso,
Stadler and Mukherjee, 2013).

OmicsMarkeR - Tools for classification and feature
selection for 'omics' level datasets.  It is a tool to
provide multiple multivariate classification and feature
selection techniques complete with multiple stability
metrics and aggregation techniques.  It is primarily
designed for analysis of metabolomics datasets but
potentially extendable to proteomics and transcriptomics
applications.

pandaR - Runs PANDA, an algorithm for discovering novel
network structure by combining information from multiple
complimentary data sources.

parglms - support for parallelized estimation of GLMs/GEEs,
catering for dispersed data

pmm - The Parallel Mixed Model (PMM) approach is suitable
for hit selection and cross-comparison of RNAi screens
generated in experiments that are performed in parallel
under several conditions. For example, we could think of
the measurements or readouts from cells under RNAi
knock-down, which are infected with several pathogens or
which are grown from different cell lines.

podkat - This package provides an association test that is
capable of dealing with very rare and even private
variants. This is accomplished by a kernel-based approach
that takes the positions of the variants into account. The
test can be used for pre-processed matrix data, but also
directly for variant data stored in VCF files. Association
testing can be performed whole-genome, whole-exome, or
restricted to pre-defined regions of interest. The test is
complemented by tools for analyzing and visualizing the
results.

PROPER - This package provide simulation based methods for
evaluating the statistical power in differential expression
analysis from RNA-seq data.

ProtGenerics - S4 generic functions needed by Bioconductor
proteomics packages.

pwOmics - pwOmics performs pathway-based level-specific
data comparison of matching omics data sets based on
pre-analysed user-specified lists of differential
genes/transcripts and proteins. A separate downstream
analysis of proteomic data including pathway identification
and enrichment analysis, transcription factor
identification and target gene identification is opposed to
the upstream analysis starting with gene or transcript
information as basis for identification of upstream
transcription factors and regulators. The cross-platform
comparative analysis allows for comprehensive analysis of
single time point experiments and time-series experiments
by providing static and dynamic analysis tools for data
integration.

QuartPAC - Identifies clustering of somatic mutations in
proteins over the entire quaternary structure.

R3CPET - The package provides a method to infer the set of
proteins that are more probably to work together to
maintain chormatin interaction given a ChIA-PET experiment
results.

RBM - Use A Resampling-Based Empirical Bayes Approach to
Assess Differential Expression in Two-Color Microarrays and
RNA-Seq data sets.

rcellminer - The NCI-60 cancer cell line panel has been
used over the course of several decades as an anti-cancer
drug screen. This panel was developed as part of the
Developmental Therapeutics Program (DTP,
http://dtp.nci.nih.gov/) of the U.S. National Cancer
Institute (NCI). Thousands of compounds have been tested on
the NCI-60, which have been extensively characterized by
many platforms for gene and protein expression, copy
number, mutation, and others (Reinhold, et al., 2012). The
purpose of the CellMiner project
(http://discover.nci.nih.gov/cellminer) has been to
integrate data from multiple platforms used to analyze the
NCI-60 and to provide a powerful suite of tools for
exploration of NCI-60 data.

RCyjs - Interactvive viewing and exploration of graphs,
connecting R to Cytoscape.js

regioneR - regioneR offers a statistical framework based on
customizable permutation tests to assess the association
between genomic region sets and other genomic features.

rGREAT - This package makes GREAT (Genomic Regions
Enrichment of Annotations Tool) analysis automatic by
constructing a HTTP POST request according to user's input
and automatically retrieving results from GREAT web server.

rgsepd - R/GSEPD is a bioinformatics package for R to help
disambiguate transcriptome samples (a matrix of RNA-Seq
counts at RefSeq IDs) by automating differential expression
(with DESeq2), then gene set enrichment (with GOSeq), and
finally a N-dimensional projection to quantify in which
ways each sample is like either treatment group.

Rhtslib - This package provides version 1.1 of the 'HTSlib'
C library for high-throughput sequence analysis. The
package is primarily useful to developers of other R
packages who wish to make use of HTSlib. Motivation and
instructions for use of this package are in the vignette,
vignette(package="Rhtslib", "Rhtslib").

RNAprobR - This package facilitates analysis of Next
Generation Sequencing data for which positional information
with a single nucleotide resolution is a key. It allows for
applying different types of relevant normalizations, data
visualization and export in a table or UCSC compatible
bedgraph file.

RnaSeqSampleSize - RnaSeqSampleSize package provides a
sample size calculation method based on negative binomial
model and the exact test for assessing differential
expression analysis of RNA-seq data

RnBeads - RnBeads facilitates comprehensive analysis of
various types of DNA methylation data at the genome scale.

RUVcorr - RUVcorr allows to apply global removal of
unwanted variation (ridged version) to real and simulated
gene expression data.

saps - Functions implementing the Significance Analysis of
Prognostic Signatures method (SAPS). SAPS provides a robust
method for identifying biologically significant gene sets
associated with patient survival. Three basic statistics
are computed. First, patients are clustered into two
survival groups based on differential expression of a
candidate gene set. P_pure is calculated as the probability
of no survival difference between the two groups. Next, the
same procedure is applied to randomly generated gene sets,
and P_random is calculated as the proportion achieving a
P_pure as significant as the candidate gene set. Finally, a
pre-ranked Gene Set Enrichment Analysis (GSEA) is performed
by ranking all genes by concordance index, and P_enrich is
computed to indicate the degree to which the candidate gene
set is enriched for genes with univariate prognostic
significance. A SAPS_score is calculated to summarize the
three statistics, and optionally a Q-value is computed to
estimate the significance of the SAPS_score by calculating
SAPS_scores for random gene sets.

SELEX - Tools for quantifying DNA binding specificities
based on SELEX-seq data

seq2pathway - Seq2pathway is a novel tool for functional
gene-set (or termed as pathway) analysis of next-generation
sequencing data, consisting of "seq2gene" and "gene2path"
components. The seq2gene links sequence-level measurements
of genomic regions (including SNPs or point mutation
coordinates) to gene-level scores, and the gene2pathway
summarizes gene scores to pathway-scores for each sample.
The seq2gene has the feasibility to assign both coding and
non-exon regions to a broader range of neighboring genes
than only the nearest one, thus facilitating the study of
functional non-coding regions. The gene2pathway takes into
account the quantity of significance for gene members
within a pathway compared those outside a pathway. The
output of seq2pathway is a general structure of
quantitative pathway-level scores, thus allowing one to
functional interpret such datasets as RNA-seq, ChIP-seq,
GWAS, and derived from other next generational sequencing
experiments.

seqPattern - Visualising oligonucleotide patterns and
sequence motifs occurrences across a large set of sequences
centred at a common reference point and sorted by a user
defined feature.

sigsquared - By leveraging statistical properties (log-rank
test for survival) of patient cohorts defined by binary
thresholds, poor-prognosis patients are identified by the
sigsquared package via optimization over a cost function
reducing type I and II error.

SIMAT - This package provides a pipeline for analysis of
GC-MS data acquired in selected ion monitoring (SIM) mode.
The tool also provides a guidance in choosing appropriate
fragments for the targets of interest by using an
optimization algorithm. This is done by considering
overlapping peaks from a provided library by the user.

similaRpeak - This package calculates metrics which assign
a level of similarity between ChIP-Seq profiles.

sincell - Cell differentiation processes are achieved
through a continuum of hierarchical intermediate
cell-states that might be captured by single-cell RNA seq.
Existing computational approaches for the assessment of
cell-state hierarchies from single-cell data might be
formalized under a general workflow composed of i) a metric
to assess cell-to-cell similarities (combined or not with a
dimensionality reduction step), and ii) a graph-building
algorithm (optionally making use of a cells-clustering
step). Sincell R package implements a methodological
toolbox allowing flexible workflows under such framework.
Furthermore, Sincell contributes new algorithms to provide
cell-state hierarchies with statistical support while
accounting for stochastic factors in single-cell RNA seq.
Graphical representations and functional association tests
are provided to interpret hierarchies.

skewr - The skewr package is a tool for visualizing the
output of the Illumina Human Methylation 450k BeadChip to
aid in quality control. It creates a panel of nine plots.
Six of the plots represent the density of either the
methylated intensity or the unmethylated intensity given by
one of three subsets of the 485,577 total probes. These
subsets include Type I-red, Type I-green, and Type II.The
remaining three distributions give the density of the
Beta-values for these same three subsets. Each of the nine
plots optionally displays the distributions of the "rs" SNP
probes and the probes associated with imprinted genes as
series of 'tick' marks located above the x-axis.

soGGi - The soGGi package provides a toolset to create
genomic interval aggregate/summary plots of signal or motif
occurence from BAM and bigWig files as well as PWM,
rlelist, GRanges and GAlignments Bioconductor objects.
soGGi allows for normalisation, transformation and
arithmetic operation on and between summary plot objects as
well as grouping and subsetting of plots by GRanges objects
and user supplied metadata. Plots are created using the
GGplot2 libary to allow user defined manipulation of the
returned plot object. Coupled together, soGGi features a
broad set of methods to visualise genomics data in the
context of groups of genomic intervals such as genes,
superenhancers and transcription factor binding events.

SVM2CRM - Detection of cis-regulatory elements using svm
implemented in LiblineaR.

TIN - The TIN package implements a set of tools for
transcriptome instability analysis based on exon expression
profiles. Deviating exon usage is studied in the context of
splicing factors to analyse to what degree transcriptome
instability is correlated to splicing factor expression. In
the transcriptome instability correlation analysis, the
data is compared to both random permutations of alternative
splicing scores and expression of random gene sets.

TPP - Analyze thermal proteome profiling (TPP) experiments
with varying temperatures (TR) or compound concentrations
(CCR).

TRONCO - Genotype-level cancer progression models describe
the ordering of accumulating mutations, e.g., somatic
mutations / copy number variations, during cancer
development. These graphical models help understand the
causal structure involving events promoting cancer
progression, possibly predicting complex patterns
characterising genomic progression of a cancer.
Reconstructed models can be used to better characterise
genotype-phenotype relation, and suggest novel targets for
therapy design. TRONCO (TRanslational ONCOlogy) is a R
package aimed at collecting state-of-the-art algorithms to
infer progression models from cross-sectional data, i.e.,
data collected from independent patients which does not
necessarily incorporate any evident temporal information.
These algorithms require a binary input matrix where: (i)
each row represents a patient genome, (ii) each column an
event relevant to the progression (a priori selected) and a
0/1 value models the absence/presence of a certain mutation
in a certain patient. The current first version of TRONCO
implements the CAPRESE algorithm (Cancer PRogression
Extraction with Single Edges) to infer possible progression
models arranged as trees; cfr. Inferring tree causal models
of cancer progression with probability raising, L. Olde
Loohuis, G. Caravagna, A. Graudenzi, D. Ramazzotti, G.
Mauri, M. Antoniotti and B. Mishra. PLoS One, to appear.
This vignette shows how to use TRONCO to infer a tree model
of ovarian cancer progression from CGH data of copy number
alterations (classified as gains or losses over
chromosome's arms). The dataset used is available in the
SKY/M-FISH database.



NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages since the last release. The following package NEWS is available:


affxparser
----------

Changes in version 1.39.4 (2015-01-18):

- ROBUSTNESS: 'GNU make' is a SystemRequirements (for now).

- ROBUSTNESS: Did not seem to be needed, but package is now a good
  citizen and do library.dynlib.unload() when unloaded.

- CLEANUP: Now using requireNamespace() instead of require().

- CLEANUP: Internal cleanup of native code.

Changes in version 1.39.3 (2014-11-26):

- BUG FIX: readPgf() and readPgfEnv() failed to read all units
  (probesets) on some systems.  Extensive package tests have been added
  to test this and other cases.  Thanks to Grischa Toedt at EMBL
  Germany for reporting on, troubleshooting, and helping out with
  patches for this bug.

Changes in version 1.39.2 (2014-10-28):

- BUG FIX: The range test of argument 'units' to readCdf() and
  readCdfQc() was never performed due to a typo, meaning it was
  possible to request units out of range.  Depending on system this
  could result in either a core dump or random garbage read for the out
  of range units.

- ROBUSTNESS: Added package system tests for out of range 'units' and
  'indices' arguments for most read functions.

Changes in version 1.39.1 (2014-10-26):

- ROBUSTNESS: Now all methods gives an informative error message if
  zero elements are requested, i.e. via zero-length argument 'indices'
  or 'units' that is not NULL. Previously this case would access all
  values just like NULL does.

- ROBUSTNESS: Now readCelRectangle() gives an informative error message
  if argument 'xrange' or 'yrange' is not of length two.

- BUG FIX: readPgf() and readPgfEnv() would give an error if argument
  'indices' was specifies as a double rather than as an integer vector.

Changes in version 1.39.0 (2014-10-13):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.1 for R (>= 3.2.0).

AllelicImbalance
----------------

Changes in version 1.6.0:

NEW FEATURES

- ASEset slot for Phase information

- generation of masked reference genomes

- calculation of reference fraction 'refFraction'

AnalysisPageServer
------------------

Changes in version 0.99.5:

- First release to Bioconductor

- See AnalysisPageServer vignette for introduction.

AnnotationDbi
-------------

Changes in version 1.30:

NEW FEATURES and API changes

- Adds mapIds() which allows users who miss mget() to extract data from
  AnnnotationDb objects: but without the dangers of using mget() (which
  fails if you pass in a bad key)

- Adds dbconn() and dbfile() methods to the list of things that
  AnnotationDb derived objects shoould always be expected to support

BUG FIXES AND CODE MAINTENANCE

- loadDb() now is smarter about whether or not it attaches a supporting
  package * * * 1.24.x SERIES NEWS * * *

AnnotationHub
-------------

Changes in version 2.0.0:

NEW FEATURES and API changes

- AnnotationHub is all new.  We basically rewrote the entire thing.

- The back end is new (new database, new way of tracking/routing data
  etc.)

- The front end is new (new AnnotationHub object, new methods, new
  behaviors, new ways of finding and downloading data)

- The metadata has also been cleaned up and made more
  consistent/searchable

- The recipes that are used to populate these data have also been
  cleaned up.

- There is also a new vignette to explain how to use the new
  AnnotationHub in detail

Improvements since last time

- The old way of finding data (an enormous tree of paths), was not
  really scalable to the amount of data we have to provide access to.
  So we junked it.  Now you have a number of better methods to allow
  you to search for terms instead.

- The new hub interface can be searched using a new display method, but
  it can *also* be searched entirely from the command line.  This
  allows you to use it in examples and scripts in a way that is
  friendlier for reproducible research.

- For users who want to contribute valuable new annotation resources to
  the AnnotationHub, it is now possible to write a recipe and test that
  it works for yourself.  Then once you are happy with it, you can
  contact us and we can add data to the AnnotationHub.


aroma.light
-----------

Changes in version 2.3.3 (2015-02-18):

- If a value of argument 'xlim' or 'ylim' for plotDensity() is NA, then
  it defaults to the corresponding extreme value of the data, e.g.
  plotDensity(x, xlim=c(0,NA)).

Changes in version 2.3.2 (2015-02-17):

- ROBUSTNESS: Added package tests. Code coverage is 76%.

- CLEANUP: Using requestNamespace() instead of request().

Changes in version 2.3.1 (2014-12-08):

- Same updates as in 2.2.1.

Changes in version 2.3.0 (2014-10-13):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.1 for R (>= 3.2.0).

Changes in version 2.2.1 (2014-12-08):

- Minor code cleanup.

ballgown
--------

Changes in version 1.99.6:

- bug fix in constructor: objects are now constructed correctly when
  only one sample is included in the object.

Changes in version 1.99.5:

- clarification of docs/warning messages when "time" variable has too
  few timepoints to fit natural splines

- bug fix in subset function: factor levels that don't appear in subset
  pData are now dropped

Changes in version 1.99.4:

- bug fix in plotting of single-transcript genes

Changes in version 1.99.3:

- unit test update due to changes in findOverlaps

Changes in version 1.99.2:

- bug fix: "ballgown" function now handles case when input feature
  names include quotes (').

Changes in version 1.99.1:

- Bioconductor 3.0 release. Version 0.99.7 is released as version
  1.0.0.

bamsignals
----------

Changes in version 0.99.0:

NEW FEATURES

- Creation of the bamsignals package, a package to load count data from
  indexed bam files efficiently.

BiGGR
-----

Changes in version 1.3.4:

- BiGGR publication is now available at PLoS ONE! Use citation("BiGGR")
  to see how to correctly cite BiGGR.

Biobase
-------

Changes in version 2.27:

NEW FEATURES

- Add write.AnnotatedDataFrame function; request of samuel.granjeaud

BiocInstaller
-------------

Changes in version 1.18.0:

NEW FEATURES

- biocLite() supports github repositories (implicitly, 'packages'
  following the 'maintainer/package' convention)



CAMERA
------

Changes in version 1.23.2:

BUG FIXES

- Fix bug in generateRules() if no ions matching the polarity are
  provided explicitly

Changes in version 1.23.1:

BUG FIXES

- Fix bug in generateRules() with empty neutraladditions

canceR
------

Changes in version 0.99.3:

- getProfilesDataMultipleGenes: line 78: Test only Genetic Profiles
  having mRNA expression to get Profile Data.

- getMutData(): line 96 change .GlobalEnv to myGlobalEnv)

- getProfDataMultipleGenes(): line 7 add testCheckedCaseGenProf()

Changes in version 0.99.2:

- Add documentation for RUN.GSEA() function

- Remove dependency of RSvgDevice. This package is not available for
  Windows OS.

Changes in version 0.99.1:

- Add examples in documents SIGNIFICANT USER-VISIBLE CHANGES

- Add documentation for RUN.GSEA() function

- Remove dependency of RSvgDevice. This package is not available for
  Windows OS.

- Add examples in documents

- Package released


Cardinal
--------

Changes in version 0.99.6:

SIGNIFICANT USER-VISIBLE CHANGES

- Bioconductor Release Candidate 7

- Added CITATION file for citing Cardinal

BUG FIXES

- In 'readAnalyze' and 'readMSIData', removed endianness check in
  Analyze 7.5 headers because some ABSciex data files specify an
  incorrect header size, thereby fixing a bug where bits would be
  swapped wrongly and file read incorrectly.

Changes in version 0.99.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Bioconductor Release Candidate 6

- Added new vignette for Cardinal design and development

- Now using ProtGenerics generics for 'spectra', 'peaks', and 'mz'

Changes in version 0.99.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Bioconductor Release Candidate 5

- Updated vignette with biological examples

- Added new citations to vignette

Changes in version 0.99.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Bioconductor Release Candidate 4

- In plot and image methods, 'groups' arg now coerced to factor

BUG FIXES

- Fixed bug in subset arg in select method

- Fixed bug in plot and image methods with NA in 'groups' arg

Changes in version 0.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Bioconductor Release Candidate 3

- Adjusts NIPALS unit tests for Windows build

Changes in version 0.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Bioconductor Release Candidate 2

- Cleaned up biocViews

BUG FIXES

- Fixed bug in SImageData coord factor levels not being properly
  updated when SImageSet is subsetted

Changes in version 0.99.0 (2014-12-22):

SIGNIFICANT USER-VISIBLE CHANGES

- Bioconductor Release Candidate 1


ChemmineOB
----------

Changes in version 1.6.0:

NEW FEATURES

- Compiled against OpenBabel 2.3.2

ChemmineR
---------

Changes in version 2.20.0:

NEW FEATURES

- Fingerprint search now provides E-values

- New SDF plotting function, draw_sdf

chimera
-------

Changes in version 1.9.2:

NEW FEATURES

- Fusions data generated by tophat-fusion-post can be imported.

- Star import function was improved.

- Annotation of fusion is more robust as it it is clearly linked to
  specific genomic releases. e.g. h19 uses BSgenome.Hsapiens.UCSC.hg19
  and TxDb.Hsapiens.UCSC.hg19.knownGene

ChIPpeakAnno
------------

Changes in version 3.0.1:

BUG FIXES

- give errors when inputs of addGeneIDs is incorrect.

ChIPseeker
----------

Changes in version 1.3.15:

- update vignette <2015-03-31, Tue>

Changes in version 1.3.14:

- add pool parameter in enrichPeakOverlap <2015-03-30, Mon>

Changes in version 1.3.13:

- update enrichPeakOverlap to support nShuffle = 0, which now will
  report only overlay with pvalue = NA <2015-03-29, Sun>

- add facet and free_y parameter for plotAvgProf and plotAvgProf2
  <2015-03-29, Sun>

- update docs <2015-03-29, Sun>

- update plotAvgProf and plotAvgProf2 to fully supporting confidence
  interval, see https://github.com/GuangchuangYu/ChIPseeker/pull/6
  <2015-03-29, Sun>

Changes in version 1.3.12:

- add confidence interval for plotAvgProf, see
  https://github.com/GuangchuangYu/ChIPseeker/issues/3 <2015-03-26,
  Thu>

Changes in version 1.3.11:

- add citation <2015-03-16, Mon>

Changes in version 1.3.10:

- update GEO data <2015-03-03, Tue>

Changes in version 1.3.9:

- add parameter *genomicAnnotationPriority* for annotatePeak function
  <2015-02-27, Fri>

Changes in version 1.3.8:

- add DOSE citation <2015-02-13, Fri>

Changes in version 1.3.7:

- bug fixed in plotDistToTSS <2015-02-06, Fri>

Changes in version 1.3.6:

- when peak is exactly located at gene end and near the end of
  chromosome, NA will be generated and throw error when assigning
  downstream of gene end.  This bug has been fixed <2015-02-03, Tue>

Changes in version 1.3.5:

- bug fixed in getNearestFeatureIndicesAndDistances when peak in the
  very begining or end of the chromosome <2015-01-30, Fri>

Changes in version 1.3.4:

- bug fixed for introducing dplyr in plotDistToTSS <2015-01-28, Wed>

Changes in version 1.3.3:

- update vignette to use BiocStyle::latex() <2015-01-26, Mon>

Changes in version 1.3.2:

- fixed import issue to meet the changes of AnnotationDbi and S4Vectors
  <2015-01-22, Thu>

cisPath
-------

Changes in version 1.7.4:

- improvements

Changes in version 1.7.3:

- Fix a bug

Changes in version 1.7.2:

- documentation improvements

Changes in version 1.7.1:

- several improvements

ClassifyR
---------

Changes in version 1.2.0:

- More classification flexibility, now with parameter tuning integrated
  into the process.

- New performance evaluation functions, such as a ROC curve and a
  performance plot.

- Some existing predictor functions are able to return class scores,
  not just class labels.

cleanUpdTSeq
------------

Changes in version 1.5.4:

NEW FEATURES

- export $ and $<- methods for class PASclassifier, modelInfo and
  featureVector.

BUG FIXES

- No changes classified as "bug fixes" (package under active
  development)

Changes in version 1.5.3:

NEW FEATURES

- No changes classified as "new features" (package under active
  development)

BUG FIXES

- Fix the bugs in test files for different output of 64bit vs 32bit
  version.

Changes in version 1.5.2:

NEW FEATURES

- No changes classified as "new features" (package under active
  development)

BUG FIXES

- Fix the bugs in documentation of predictTestSet.

Changes in version 1.5.1:

NEW FEATURES

- Add buildClassifier function to save the classifier.

- Add test file

- predictTestSet also give a invisible output.

- Add class "modelInfo", "PASclassifier" and "featureVector"

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

cleaver
-------

Changes in version 1.5.3 (2015-02-02):

- Set taxId=9606 for UniProt.ws in vignette.

Changes in version 1.5.1 (2015-01-24):

- Add test case if number of cleavage sites is smaller than the number
  of allowed missed cleavages.

Clomial
-------

Changes in version 1.3.2 (2015-04-04):

- Minor error checking for when the total counts are zero.

clusterProfiler
---------------

Changes in version 2.1.14:

- [enrichDAVID] report NA in qvalue column when it fail to calculate
  instead of throw error <2015-03-10, Fri>

Changes in version 2.1.13:

- update vignette <2015-03-07, Tue>

Changes in version 2.1.12:

- user input is now supported with enricher and GSEA function
  <2015-04-01, Wed>

- change use_internal_data parameter to ... according to the change of
  DOSE <2015-04-01, Wed>

Changes in version 2.1.11:

- add example of bitr in vignette <2015-03-31, Tue>

Changes in version 2.1.10:

- idType and bitr function for translating biological ID. <2015-03-20,
  Fri>

Changes in version 2.1.9:

- throw stop msg when gene can not be mapped by DAVID <2015-03-12, Thu>

- add return NULL when DAVID fail (no enrichment found or wrong ID
  type) <2015-03-11, Wed>

Changes in version 2.1.8:

- implement enrichDAVID <2015-03-07, Fri>

Changes in version 2.1.7:

- set includeAll = TRUE by default <2015-03-03, Tue>

- added an includeAll option in plot(compareCluster) <2015-02-26, Thu>
  refer to: https://github.com/GuangchuangYu/clusterProfiler/issues/9

- fixed typo in groupGO <2015-02-26, Thu>

Changes in version 2.1.6:

- add DOSE citation <2015-02-13, Fri>

Changes in version 2.1.5:

- gseKEGG is also use online KEGG data by default. <2015-02-11, Wed>

- change 'use.KEGG.db' parameter to 'use_internal_data', that can be
  used when I implement GO analysis with user specific annotation data
  <2015-02-11, Wed>

Changes in version 2.1.4:

- update vignette <2015-02-11, Wed>

- support formula like x ~ y + z, contributed by Giovanni Dall'Olio
  <2015-02-11, Wed>

- add example of formula interface <2015-02-11, Wed>

- formula interface for compareCluster contributed by Giovanni
  Dall'Olio <2015-02-11, Wed>

Changes in version 2.1.3:

- update vignettes <2015-01-28, Wed>

- fully support of using online version of KEGG in enrichKEGG, using
  online version is now by default of enrichKEGG <2015-01-28, Wed>

- update CITATION <2015-01-28, Wed>

Changes in version 2.1.2:

- support downloading latest KEGG annotation data <2015-01-27, Tue>

- update vignette to use BiocStyle <2015-01-27, Tue>

Changes in version 2.1.1:

- add BugReports URL <2014-12-17, Wed>

CNAnorm
-------

1.13.1: chromosome does not start from 1

cogena
------

Changes in version 0.99.3 (2015-04-02):

- use makeCluster to creat cluster to avoid checking error from
  Windows.

Changes in version 0.99.2 (2015-04-01):

- Update indent based on the coding style of Bioconductor.

- Trim long line.

Changes in version 0.99.1 (2015-03-31):

- Add runnable examples for functions.

- Update R version dependency from 2.10 to 3.2.

- Add non-empty value sections for some functions.

Changes in version 0.99.0 (2015-03-27):

- Submit cogena to Bioconductor.

coMET
-----

Changes in version 0.99.10 (2015-04-10):

- Update regulationBiomart function and the pre-computed data
  associated because ENSEMBL has recently updated their schema in
  Biomart for these data.

- Update the documentations (man and vignette)

Changes in version 0.99.9 (2015-03-10):

- Update snpBiomart and structureBiomart functions and the pre-computed
  data associated because ENSEMBL has recently updated their schema in
  Biomart for these data.

- Update the documentations (man and vignette)

- Add the function compute.pvalue.cormatrix and allow the visualisation
  of only correlation under a significant level. (use the library
  psych)

- Add the function comet.list to list the correlations between omic
  features

- Replace capital letter by low letter

- Update Shiny-based web application according to new functions

Changes in version 0.99.8 (2015-02-15):

- Correction genesNameENSEMBL in order not to run if no ENSEMBL gene in
  the region of interest

- Correction the error related to the loading of correlation matrice

- Update the visualisation of Chromosome with the name of band (related
  to an update of Gviz)

- Remove the visualisation of border for genes and transcripts to help
  the visualisation of different exons

- Update the vignette about comet's Shiny website and colors of tracks

- Update the annotation tracks as Gviz was updated and now coMET can
  work on Gviz (version 1.10.9, current version)

Changes in version 0.99.7 (2014-12-19):

- Update the function "genesENSEMBL", "genesNameENSEMBL",
  "transcriptENSEMBL", and "regulationBiomart" because ENSEMBL mart
  changed the names of field in GChr37

Changes in version 0.99.6 (2014-11-25):

- add in "import IRanges and S4Vectors" in NAMESPACE

- Fix mistakes in the vignette

Changes in version 0.99.5 (2014-11-06):

- Update the annotation tracks using Gviz because of the update of Gviz
  (version 1.11.2, development version)

Changes in version 0.99.4 (2014-11-05):

- Update the manual and the vignette

- Fix little bugs

- Update the annotation tracks using BiomartGeneRegionTrack of Gviz

Changes in version 0.99.3 (2014-10-25):

- Update the manuel and the vignette

- Update the function to create annotation tracks connecting to ENSEMBL
  mart (change of schema)

- Fix little bug

- Add the function to define the reference genome (cf the change of
  schema of ENSEMBL mart)

Changes in version 0.99.2 (2014-10-16):

- Update the manual comet.web file

Changes in version 0.99.1 (2014-10-16):

- Add the package "BiocStyle" in Suggests of DESCRIPTION file

- Change the absolute path to relative path of files (info file,
  expression file, and correlation file) in vignette and manual

Changes in version 0.99.0 (2014-09-24):

- Version draft of this script

compcodeR
---------

1.3.1: Added citation file

ComplexHeatmap
--------------

Changes in version 0.99.4:

- fixed a bug when setting `cluster_rows` to FALSE but still cluster on
  rows.

Changes in version 0.99.2:

- add two examples in vignette

- add chunk labels in the vignette

Changes in version 0.99.1:

- x and y in `cell_fun` are now `unit` objects.

Changes in version 0.99.0:

- First release

CopywriteR
----------

Changes in version 1.99.4 (2015-04-08):

- Fixed bug related to logging when .bam are not indexed

Changes in version 1.99.3 (2015-04-06):

- Fixed bug related to testing of sort mode of .bam files

Changes in version 1.99.2 (2015-03-21):

- Provided vignette that allows R CMD CHECK --no-build-vignettes in <5
  min

- Fixed minor bug affecting runs in which not all samples are analyzed
  by CopywriteR function

- Included a test for sorting of .bam files

- Improved error messages

- Fixed minor bug affecting plotting of samples

- Fixed minor bug causing error when processing >100 samples
  simultaneously

- Added support for capture regions bed files lacking the correct
  chromosome name prefix

- Fixed minor bug affecting read count statistics on single-sample
  analysis

Changes in version 1.99.0 (2015-03-05):

- Released new version

- Improved the speed of the package

- Implemented logging with the futile.logger package

- Replaced direct slot or field access with the use of accessor
  functions

- Implemented text wrapping

- Replaced the use of bed-file format with GRanges objects

- Created annotation and experiment packages

- Addressed all NOTES upon R CMD CHECK

- Restore options upon exit

- Adjusted version numbering to 1.99.0

- Implemented BiocParallel for parallel computing

- Changed name from ENCODER to CopywriteR

- Shortened lines to 80 characters

- Added NEWS file

- Changed F to FALSE

- Changed T to TRUE

- Added biocViews

- Changed version to 0.99.0

cpvSNP
------

Changes in version 0.99.0:

- Package submitted to Bioconductor

CRISPRseek
----------

Changes in version 1.7.6:

NEW FEATURES

- added parameter foldgRNA to calculate the minimum free energy (mfe)
  of the secondary structure of the gRNA with gRNA backbone sequence.
  In addition, summary file also includes the difference of mfe between
  the secondary structure of gRNA backbone alone and the secondary
  structure of gRNA backbone with the variable region of the gRNA, and
  the bracket notation of the secondary structure of gRNA backbone with
  the variable region of the gRNA to facilitate gRNA selection.

Changes in version 1.7.3:

NEW FEATURES

- added parameter chromToExclude in offTargetAnalysis to specify
  chromosomes not to search for offtargets, which can be used to
  exclude haplotype blocks

Changes in version 1.7.1:

BUG FIXES

- only search for gRNAs for input sequences longer than gRNA.size plus
  PAM size

csaw
----

Changes in version 1.1.28:

- Added getWidths(), scaledAverage() and filterWindows(), to facilitate
  comparison of abundances during filtering.

- Added findMaxima() to identify locally maximal windows from range
  data.

- Added profileSites() to examine the coverage profile around specified
  regions, with wwhm() to guess the ideal window size from the profile.

- Changed default window width in windowCounts() to 50 bp, default
  filter to a fixed count of 10. Also, filter=0 is honored when
  bin=TRUE.

- Switched from the depracated rowData to rowRanges for all
  manipulations of SummarizedExperiment.

- Changed all instances of `pet' to `pe' in read parameter
  specification, and renamed getPETSizes() to getPESizes().

- Removed the redundant rescue.pairs parameter in readParam().

- Added fast.pe parameter in readParam(), for fast paired-end data
  extraction. Added dumpPE() to pre-process paired-end BAM files for
  fast downstream extraction.

- Added support for custom column specification in getBestTest(),
  combineTests().

- Switched from reporting average log-FC to numbers of up/down windows
  in combineTests().

- Allowed getBestTest() to return all fields associated with the best
  window in the output table.

- Added upweightSummits() to compute weights favouring high-abundance
  windows.

- Added combineOverlaps(), getBestOverlaps() and summitOverlaps()
  wrapper functions for processing of Hits objects.

- Added consolidateSizes(), to consolidate DB results from multiple
  window sizes.

- Added support for custom key/name specification in detailRanges() for
  non-human/mouse systems.

- Added support for strand-specific read extraction in readParam(),
  strand-specific counting via strandedCounts().

- Added strand-awareness to mergeWindows(). Added protection against
  stranded input regions in extractReads(), detailRanges().

- Changed algorithm for splitting of large peaks in mergeWindows().

- Stored counting parameters in exptData for windowCounts(),
  regionCounts().

- Fixed small inaccuracies with continuity correction addition in
  normalizeCounts() for NB-loess.

- Switched to fragment midpoint for binning of paired-end data in
  windowCounts().

- Added support for lists of library-specific readParam objects in
  windowCounts(), regionCounts(), correlateReads().

- Added makeExtVector(), to support variable read extension lengths
  between libraries in windowCounts(), regionCounts().

- Added support for read extension within extractReads().

- Updated the user's guide to reflect new and modified functions.

- Added sra2bam.sh in inst/doc to reproducibly generate BAM files to
  run UG examples.

- Cleaned up code in inst/tests for modified functions, added new tests
  for new functions.

cummeRbund
----------

2.9.3: Bugfix: - Introduced CHECK error by adding to
        .Rbuildignore...this is now fixed.

2.9.2: version bump to let BioC nightly build grab commit.

2.9.1: version bump for BioC devel release 3.1

2.8.2: Bugfixes: - removed reference to sqliteCloseConnection() (not
        exported by RSQLite 1.0.0) in vignette.

2.8.1: Bugfixes: - Made minimal changes for compatibility with RSQLite
        1.0.0

dagLogo
-------

Changes in version 1.5.1:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- try to catch the error when mart server does not work.

DART
----

Changes in version 1.15.3:

- Function DoDARTCLQ is added. And some arguments and values changes in
  previous functions.

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

1.01.01: 10-17-2014 Lorena Pantano <lorena.pantano@gmail.com> FIX BUG
        WHEN ONLY ONE GENE IS DEG

derfinder
---------

Changes in version 1.1.18:

BUG FIXES

- Updated to work with qvalue 1.99.0

Changes in version 1.1.17:

SIGNIFICANT USER-VISIBLE CHANGES

- Changed citation information to reference the bioRxiv pre-print

BUG FIXES

- Polished the interaction with bumphunter >= 1.7.3

- Updated the default cluster option now that BiocParallel::SnowParam()
  no longer has an 'outfile' argument.

Changes in version 1.1.16:

SIGNIFICANT USER-VISIBLE CHANGES

- analyzeChr() now uses annotateTranscripts() and matchGenes() from
  bumphunter version 1.7.3 (or greater). As announced at
  https://support.bioconductor.org/p/63568/ these changes in bumphunter
  allow straight forward use of non-human annotation. In analyzeChr()
  using a different organism can be used by changing the 'txdb'
  argument: finer control can be achieved through '...'. For example,
  by specifying the 'annotationPackage' argument used in
  annotateTranscripts().

Changes in version 1.1.15:

BUG FIXES

- makeGenomicState() incorrectly labeled regions as intragenic. The
  correct name is intergenic.

Changes in version 1.1.14:

BUG FIXES

- Fixed an important bug on calculatePvalues()! Basically, internally
  `maxRegionGap` was set to 300 instead of 0 in one step by default.
  Thus the process of mapping regions to genomic coordinates was messed
  up. If you have results prior to this fix you can try using
  https://gist.github.com/bf85e2c7d5d1f8197707 to fix the results as
  much as possible. Basically, regions will be correct but the p-values
  will be approximated with the available information from the null
  regions. Truly fixing the p-values can only be done by re-running
  derfinder.

Changes in version 1.1.5:

NEW FEATURES

- Introduced function extendedMapSeqlevels() for using GenomeInfoDb
  when there is information regarding the species and naming style of
  interest. Otherwise sequence names are left unchanged. If used with
  'verbose = TRUE', a message is printed whenever GenomeInfoDb could
  not be used or if some information had to be guessed.

BUG FIXES

- Fixes https://support.bioconductor.org/p/62136

Changes in version 1.1.3:

NEW FEATURES

- loadCoverage() and fullCoverage() now support BamFile and BigWigFile
  objects

BUG FIXES

- Fixed a bug in loadCoverage() when the input was a BamFileList.
  Implemented tests based on the bug. Bug reported at
  https://support.bioconductor.org/p/62073

derfinderPlot
-------------

Changes in version 1.1.6:

SIGNIFICANT USER-VISIBLE CHANGES

- Adapted to work with bumphunter version >= 1.7.6

Changes in version 1.1.3:

BUG FIXES

- Adapted plotCluster() and plotOverview() to derfinder 1.1.5

DESeq2
------

Changes in version 1.8.0:

- Added support for user-supplied model matrices to DESeq(),
  estimateDispersions() and nbinomWaldTest(). This helps when the model
  matrix needs to be edited by the user.

Changes in version 1.7.45:

- Added a test in rlog for sparse data, mostly zero and some very large
  counts, which will give a warning and suggestion for alternate
  transformations.

- Added plotSparsity() which will help diagnose issues for using rlog:
  data which do not resemble negative binomial due to many genes with
  mostly zeros and a few very large counts.

Changes in version 1.7.43:

- Added 'replaced' argument to counts() and plotCounts() such that the
  assay in "replaceCounts" will be used if it exists. Raised a minimum
  dispersion value used in Cook's calculation, so that other counts in
  a group with an outlier won't get extreme Cook's distances
  themselves.

Changes in version 1.7.32:

- Added logic to results() which will zero out the LFC, Wald statistic
  and set p-value to 1, for 'contrast' argument results tables where
  the contrasted groups all have zero count. Non-zero LFCs were
  otherwise occuring due to large differences in the size factors.

Changes in version 1.7.11:

- Added support for user-supplied model matrices to DESeq(),
  estimateDispersions() and nbinomWaldTest().

Changes in version 1.7.9:

- Added Genome Biology citation for the DESeq2 methods.

- Introduced type="iterate" for estimateSizeFactors, an alternative
  estimator for the size factors, which can be used even when all genes
  have a sample with a count of zero. See man page for details.

Changes in version 1.7.3:

- Fixed two minor bugs: DESeq() with parallel=TRUE was dropping rows
  with all zero counts, instead of propogating NAs. nbinomLRT() with
  matrices provided to 'full' and 'reduced' and a design of ~ 1, the
  matrices were being ignored.

DiffBind
--------

Changes in version 1.13.5:

- Add support for reading Excel-format sample sheets (.xls, .xlsx
  extensions

Changes in version 1.13.4:

- Roll up all updates

- Update DESeq2 reference in vignette; fix vignette samplesheet

- use vennPlot from systemPiper

Changes in version 1.13.3:

- Fix Makevars to avoid gnu-specific extensions

- Replace 'require' with 'requireNamespace' to eliminate NOTEs
  regarding misuse of 'require'

- Remove non-ASCII characters from a couple of comments

Changes in version 1.13.2:

- Change Gord's email address

Changes in version 1.13.1:

- New: color vector lists for dba.plotHeatmap and colos for dba.plotPCA
  labels

- Fix: bug causing two plots when changing score in dba.plotHeatmap and
  dba.plotPCA

diffHic
-------

Changes in version 0.99.0:

- New package diffHic, for detection of differential interactions from
  Hi-C data.

DOSE
----

Changes in version 2.5.12:

- report NA in qvalue column when it fail to calculate instead of throw
  error <2015-03-09, Thu>

- update IC data <2015-03-08, Wed>

Changes in version 2.5.11:

- implement clusterSim and mclusterSim <2015-03-07, Tue>

Changes in version 2.5.10:

- implement enrichNCG and now gseAnalyzer accept setType = "NCG"
  <2015-04-01, Wed> see http://ncg.kcl.ac.uk

- change use_internal_data parameter to ... in all S3 generics and
  methods <2015-04-01, Wed>

Changes in version 2.5.9:

- extend gseAnalyzer to support use_internal_data parameter
  <2015-03-29, Sun>

Changes in version 2.5.8:

- keep order of barplot, sorted by pvalue by defalt <2015-02-11, Wed>

Changes in version 2.5.7:

- fixed typo in vignette <2015-02-10, Tue>

Changes in version 2.5.6:

- add organismMapper to satisfy the change of enrichKEGG <2015-01-28,
  Wed>

Changes in version 2.5.5:

- introduce use_internal_data parameter for enrichKEGG of
  clusterProfiler <2015-01-27, Tue>

Changes in version 2.5.3:

- update vignette <2015-01-27, Tue>

Changes in version 2.5.1:

- add CITATION <2015-01-19, Mon>

DSS
---

Changes in version 2.5.0:

- Implemented methods for detecting differentially methylated regions
  (DMR).

easyRNASeq
----------

Changes in version 2.3.4:

- Bioc. API changes: added missing imports, removed unneeded imports,
  added a missing package dependency

- Removed generics now defined in BiocGenerics

Changes in version 2.3.3:

- Bioc. API changes

Changes in version 2.3.2:

- Bioc. API changes

Changes in version 2.3.1:

- Same as release 2.2.1

Changes in version 2.3.0:

- Bioconductor devel 3.1


EBImage
-------

Changes in version 4.10.0:

NEW FEATURES

- 'paintObjects' allows control over behavior at image edges

- 'getFrames' function returning a list of Image objects

- default 'display' method can be set by options("EBImage.display")

- 'short' argument to the Image 'print' method

- 'equalize' function performing histogram equalization

SIGNIFICANT USER-VISIBLE CHANGES

- 'translate' moves the image in the opposite direction as before

- 'rotate' returns an image centered in a recalculated bounding box

- 'as.Image' coerces subclasses to Image objects

- 'getFrame' returns individual channels of a Color image as Grayscale
  frames

- 'display' method defaults to 'raster' unless used interactively

- 'combine' allows to combine images of the same color mode only

PERFORMANCE IMPROVEMENTS

- 'filter2': calculate FFT using methods from the 'fftwtools' package

- 'as.nativeRaster': now implemented in C

- 'medianFilter': double to integer conversion moved to C code

- 'hist': fixed binning issues and improved performance

BUG FIXES

- 'medianFilter': fixed multi-channel image handling; preserve original
  object class

- 'as.nativeRaster': allow for an arbitrary number of color channels

- 'display': set missing color channels to blank

EBSeq
-----

Changes in version 1.7.1:

- In EBSeq 1.7.1, EBSeq incorporates a new function GetDEResults()
  which may be used to obtain a list of transcripts under a target FDR
  in a two-condition experiment.  The results obtained by applying this
  function with its default setting will be more robust to transcripts
  with low variance and potential outliers.  By using the default
  settings in this function, the number of genes identified in any
  given analysis may differ slightly from the previous version (1.7.0
  or order).  To obtain results that are comparable to results from
  earlier versions of EBSeq (1.7.0 or older), a user may set
  Method="classic" in GetDEResults() function, or use the original
  GetPPMat() function.  The GeneDEResults() function also allows a user
  to modify thresholds to target genes/isoforms with a pre-specified
  posterior fold change.

- Also, in EBSeq 1.7.1, the default settings in EBTest() and
  EBMultiTest() function will only remove transcripts with all 0's
  (instead of removing transcripts with 75th quantile less than 10 in
  version 1.3.3-1.7.0).  To obtain a list of transcripts comparable to
  the results generated by EBSeq version 1.3.3-1.7.0, a user may change
  Qtrm = 0.75 and QtrmCut = 10 when applying EBTest() or EBMultiTest()
  function.

EDASeq
------

Changes in version 2.1:

- Fixed bug in plotPCA: parameter k was not passed properly in
  SeqExpressionSet method.

- Fixed bug in plotQuality: trimming when reads of different lengths in
  same BAM/FASTQ files.

edgeR
-----

Changes in version 3.10.0:

- An DGEList method for romer() has been added, allowing access to
  rotation gene set enrichment analysis.

- New function dropEmptyLevels() to remove unused levels from a factor.

- New argument p.value for topTags(), allowing users to apply a p-value
  or FDR cutoff for the results.

- New argument prior.count for aveLogCPM().

- New argument pch for the plotMDS method for DGEList objects. Old
  argument col is now removed, but can be passed using .... Various
  other improvements to the plotMDS method for DGEList objects, better
  labelling of the axes and protection against degenerate dimensions.

- treatDGE() is renamed as glmTreat(). It can now optionally work with
  either likelihood ratio tests or with quasi-likelihood F-tests.

- glmQLFit() is now an S3 generic function.

- glmQLFit() now breaks the output component s2.fit into three separate
  components: df.prior, var.post and var.prior.

- estimateDisp() now protects against fitted values of zeros, giving a
  more accurate estimate of prior.df.

- DGEList() now gives a message rather than an error when the count
  matrix has non-unique column names.

- Minor corrections to User's Guide.

- requireNamespace() is now used internally instead of require() to
  access functions in suggested packages.


EnrichmentBrowser
-----------------

Changes in version 1.1.1:

NEW FEATURES

- de.ana: now based on limma functionality

- ggea.graph: extended control of layout

- ggea: * permutation p-value now based on fast edge resampling *
  p-value approximation now based on gaussian mixture * extended
  control for edge selection (consistency threshold, edge type, ...)

ensembldb
---------

Changes in version 0.99.17:

NEW FEATURES

- Added new function ensDbFromGRanges that builds an EnsDB database
  from information provided in a GRanges object (e.g. retrieved from
  the AnnotationHub).

Changes in version 0.99.16:

SIGNIFICANT USER-VISIBLE CHANGES

- Added argument outfile to ensDbFromGtf that allows to manually
  specify the file name of the database file.

- ensDbFromGtf tries now to automatically fetch the sequence lengths
  from Ensembl.

BUG FIXES

- Fixed the function that extracts the genome build version from the
  gtf file name.

Changes in version 0.99.15:

NEW FEATURES

- metadata method to extract the information from the metadata database
  table.

- ensDbFromGtf function to generate a EnsDb SQLite file from an
  (Ensembl) GTF file.

Changes in version 0.99.14:

BUG FIXES

- Fixed a problem when reading tables fetched from Ensembl that
  contained ' or #.

Changes in version 0.99.13:

SIGNIFICANT USER-VISIBLE CHANGES

- Added argument "port" to the fetchTablesFromEnsembl to allow
  specifying the MySQL port of the database.

Changes in version 0.99.12:

BUG FIXES

- argument "x" for method organism changed to "object".

ensemblVEP
----------

Changes in version 1.8.0:

NEW FEATURES

- Add VEPParam78 class to support Ensembl 78

- Add 'flag_pick', 'flag_pick_allele', 'pick_order' and 'tsl' flags

MODIFICATIONS

- Add new biocViews terms

erccdashboard
-------------

Changes in version 1.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Version incremented to 1.1.0 for Bioconductor 3.1 release.

BUG FIXES AND MINOR IMPROVEMENTS

- Updated NEWS file, DESCRIPTION file,and improved documentation in
  vignette.

FGNet
-----

Changes in version 3.2:

NEW FEATURES

- HTML vignette

- Bipartite network: Nodes are now circles and squares (only available
  for "static" plot)

- PlotGoAncestors: added argument nCharTerm

- clustersDistance: added argument 'clustMethod' and set "average" as
  default clustering method.

- GO evidence: Added argument to fea_topGO (new data object available:
  GOEvidenceCodes)

- Added argument vExprColors to functionalNetwork.

- Vertex size now allows to set a value for each gene.

BUGFIXES

- fea2incidMat: filter negative values

FISHalyseR
----------

Changes in version 0.99:

- Multi probe support

- Improved speed

- Parameters are now stored in arguments.txt file

flipflop
--------

Changes in version 1.5.15:

- correct a bad memory managment in readgroup.cpp (do not use erase on
  a vector!)

- correct a bug that may happen when the pre-processing encountered a *
  in the CIGAR field.  now it simply ignore the read (align.cpp).

Changes in version 1.5.14:

- correct a bug which might happened when reading the RG:Z tag in the
  SAM file (align.cpp)

- allow H characters in the CIGAR field

- correction in spams (prox/project.h) related to topological order

Changes in version 1.5.12:

MAJOR CHANGES

- parallelization -----> FlipFlop can be run on multiple cores.

NEW FEATURES

- 2 novel options related to the parallelization: - 'sliceCount' number
  of slices in which the pre-processing '.instance' file is cut.  when
  OnlyPreprocess is set to TRUE, it creates slices and you can run
  FlipFlop independently afterwards on each one of the slice.  -
  mc.cores it automatically distributes the slices on several cores
  when sliceCount>1.  when using a given preprocessing.instance and
  several samples, it uses several cores for some operations concerning
  multiple samples.

Changes in version 1.5.11:

- added the fraction of each predicted transcript (among one gene) in
  the output GTF file.

Changes in version 1.5.10:

MAJOR CHANGES

- implementation of a MULTI-SAMPLE version of FlipFlop -----> the
  multi-sample procedure uses several samples simultaneously in order
  to increase statistical power.  by default the multi-sample procedure
  is a penalized likelihood with the group-lasso penalty.  you can also
  choose to perform a simple refit after using the pooled data.

NEW FEATURES

- new features as a result of the multi-sample extension. Novel options
  are - 'samples.file' optional samples file (one line per sample name)
  - 'output.type' type of output when using several samples
  simultaneously.  when equal to "gtf" the output corresponds to a gtf
  file per sample with specific abundances.  when equal to "table" the
  output corresponds to a single gtf file storing the structure of the
  transcripts, and an associated table with the abundances of all
  samples.  - 'mergerefit' if TRUE use a simple refit strategy instead
  of the group-lasso. Default FALSE.

MINOR CHANGES

- change a bit the paired-end graph paradigm. Ensure that one isoform
  is one unique path in the graph, as for single-end graph.

- [parameter] add the delta (poisson offset in the loss) parameter as
  an option, in case of fine-tuning.

- [parameter] the cutoff option is now active for both single and
  paired, and the default value is 1%.

Changes in version 1.5.8:

- change default MinCvgCut value from 0.25 to 0.05

Changes in version 1.5.7:

- correct a stupid error (in flipflop.R) for calculating the total
  number of mapped fragments in the paired-end case

Changes in version 1.5.6:

MAJOR CHANGES

- faster single-end and paired-end graph implementation

- use_TSSPAS: if yes then start/end bins are the ones with no
  entering/outgoing junctions (+ annotated ones if given annotations)

MINOR CHANGES

- remove the sanity-check for possibly duplicated types, as it was
  solved into readgroup.cpp/getType

- [pre-processing] back to default 1 for minJuncCount (I changed my
  mind !!)

- [pre-processing] cancel change of version 1.3.9 (extreme boundary of
  segments calculateBound in readgroup.cpp) (I changed my mind !!)

Changes in version 1.5.5:

- allow more easy 2 steps run (preprocessing+prediction). The option
  "OnlyPreprocess" allows to perform preprocessing only. It creates two
  files 'prefix.instance' and 'prefix.totalnumread'. Running then
  flipflop with the option preprocess.instance='prefix.instance'
  performs the prediction (it reads automatically the total number of
  reads in prefix.totalnumread)

- [pre-processing] push the default minJuncCount from 1 to 2 (high
  enough for a default?)

- [pre-processing] slight change for paired-end preprocessing, see
  option to processsam (--single-only not automatic anymore)

- start using Rprintf in c/cpp codes for user messages (to continue)

Changes in version 1.5.4:

- [pre-processing] solve small issue on CIGAR character 'S,X,='
  (again!), see align.cpp

Changes in version 1.5.3:

NEW FEATURE

- [pre-processing] add minJuncCount as an option, ie number of required
  reads to consider a junction/boundary as valid.

Changes in version 1.5.1:

- [pre-processing] solve small issue of possibly duplicated read type,
  see getType in readgroup.cpp

flowcatchR
----------

Changes in version 1.2.0:

NEW FEATURES

- addTrajectory: changed implementation in the addTrajectory function,
  notable gain in computing time

- particles: implementation of particles() now used BiocParallel -
  massive gain in performances!

- snap() function to allow interaction of the user with the image, to
  display additional information on the trajectories identified on the
  Frames/ParticleSet

- shinyFlow() released! A Shiny App is delivered along with the
  package, documentation is available to guide exploration of data and
  parameters

- Jupyter/IPython notebooks available in the installed folder of the
  package - provided as template for further analysis

- Vignette: updated to reflect the current implementation and newest
  features

- preprocess.Frames: extended documentation

- BiocViews: added a couple of terms that fit better with current
  status of the package

- Added url for development version on github.  Active development is
  performed there, therefore feel free to send pull requests/ feature
  wishes (https://github.com/federicomarini/flowcatchR/issues)!

- Added integration of the repository with Travis-CI

- add.contours: if image is in grayscale, use false colours to enhance
  detected particles

- plot2d.TrajectorySet is now able to draw a grid and help backtracking
  points from the trajectories

BUG FIXES

- select.Frames: fixed behaviour in select.Frames when selecting a
  single frame

- add.contours: corrected bug in add.contours for grayscale images, was
  preventing correct combining of images. now uses correctly false
  colours and aids identification

- particles: for computing particles,added check whether dimnames is
  set or not. previously it raised an error, but it was kind of
  overkilling

- fixed behaviour in export.Frames

- typos fixed

- Vignette: fixed typos, fixed parameter names in the preprocessing
  steps

- read.Frames: additional check on the existence of all files

- Grayscale images can now be directly preprocessed by selecting
  channel = "all"


flowMap
-------

Changes in version 1.5.1:

USER VISIBLE CHANGES

- Remove functions that support computation of permutation-based
  p-values of the Friedman- Rafsky (FR) statistic,including
  permutePops() and getFRvalsPerm(). Current release computes p-values
  of the FR statistic directly from the standard normal curve.

- Remove functions for building a reference sample for multiple flow
  cytometry sample comparisons, including makeRefMap() and
  makeRefSample(). Instead, curent release encourages the user to
  generate a multi-sample similar matrix of cell populations which can
  in turn be used to decide on matched or mismatched cell populations

- Remove getSomePops(), which was used to compute selected cell
  population pairs in a cross-sample comparison. Current release
  computes the FR statistics for all possible cell population pairs in
  a cross-sample comparison.

FlowRepositoryR
---------------

Changes in version 0.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Added various accessor methods.

Changes in version 0.99.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Included addition unit tests.

- Version first submitted to BioConductor.


FlowSOM
-------

Changes in version 0.99.4:

NEW FEATURES

- Extra functionality to differentiate between groups of samples

Changes in version 0.99.0:

NEW FEATURES

- First version of the package

flowType
--------

Changes in version 2.7.0:

- Minor bugfix to make partitions matrix come out correctly when more
  than one threshold is used.

frma
----

Changes in version 1.19:

- Correct citation for use of frma with Gene and Exon ST arrays added.

- The internal workings of the frma function have been modified to ease
  the eventual transition AffyBatch (affy) objects to FeatureSet
  (oligo) objects as input.

- Changed the way in which the user can pass their own fRMA vectors to
  the frma function via the "input.vecs" argument. Previously, the user
  could supply a partial set of vectors and the remaining ones would be
  obtained from the default frmavecs packages. In retrospect, this
  seems very likely to introduce a mistake. Now, the user has to either
  supply all of their own vectors or use all of the default vectors. If
  anyone was mixing vectors before, this can still be accomplished by
  loading the respective frmavecs package and manually creating your
  own mixture of vectors.

- Bug fix: on certain operating systems, due to differences in the
  locale, the probe ordering returned by the pm function (affy) did not
  match the order stored in the frmavecs packages. Updated frmavecs
  packages will now contain the pmindex values for all probes and affy
  ids for all probesets. Also the frma and barcode functions now check
  that the order matches. If the order is not the same, these functions
  attempt to reorder the probes to match.

- Added NEWS.Rd file.

gCMAP
-----

Changes in version 1.11.3:

- BUGFIX: The wiki2cmap function now excludes empty files downloaded
  from wikipathways from parsing.

- BUGFIX: The go2cmap function is now exported.

- BUGFIX: The reactom2cmap function is currently defunct, because the
  reactome.db package is broken.

Changes in version 1.11.2:

- CHANGE: Fixed typo in mroast_score method.

Changes in version 1.11.1:

- BUGFIX: Fixed the romer_score method to accomodate the change in
  argument order in limma's romer function.

Changes in version 1.11.0:

- CHANGE: Version bump for BioC update.

gCMAPWeb
--------

Changes in version 1.7.2:

- BUGFIXRemoved redundant Rd file.

Changes in version 1.7.1:

- BUGFIXUpdated roxygen tags and specified imports in a more granular
  fashion.

gdsfmt
------

Changes in version 1.3.0-1.3.10:

NEW FEATURES

- a new argument 'visible' is added to the functions `add.gdsn()`,
  `addfile.gdsn()` and `addfolder.gdsn()`

- `objdesp.gdsn()` returns 'encoder' to indicate the compression
  algorithm

- add a new function `system.gds()` showing system configuration

- support efficient random access of zlib compressed data, which are
  composed of independent compressed blocks

- support LZ4 compression format (http://code.google.com/p/lz4/), based
  on "lz4frame API" of r128

- allow R RAW data (interpreted as 8-bit signed integer) to replace
  32-bit integer with `read.gdsn()`, `readex.gdsn()`, `apply.gdsn()`,
  `clusterApply.gdsn()`, `write.gdsn()`, `append.gdsn()`

- a new argument 'target.node' is added to `apply.gdsn()`, which allows
  appending data to a target GDS variable

- `apply.gdsn()`, `clusterApply.gdsn()`: the argument 'as.is' allows
  "logical" and "raw"

- more argument checking in `write.gdsn()`

- new components 'trait' and 'param' in the return value of
  `objdesp.gdsn()`

- add new data types: packedreal8, packedreal16, packedreal32

- a new argument 'permute' in the function `setdim.gdsn()`

SIGNIFICANT USER-VISIBLE CHANGES

- v1.3.5: add a R Markdown vignette

BUG FIXES

- minor fixes

- v1.3.2: fixes https://github.com/zhengxwen/gdsfmt/issues/7

- v1.3.3: minor fixes, 'sync.gds' synchronizes the GDS file

- v1.3.4: fixes https://github.com/zhengxwen/gdsfmt/issues/6

- v1.3.5: fixes https://github.com/zhengxwen/gdsfmt/issues/8

GeneNetworkBuilder
------------------

Changes in version 1.9.1:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- Update author names.

GENESIS
-------

Changes in version 0.99.4:

- Fixed a minor bug to zero out diagonal of divMat in pcairPartition()

Changes in version 0.99.0:

- Initial version of GENESIS.  Contains functionality to perform
  PC-AiR.

geNetClassifier
---------------

Changes in version 1.6.1:

- Vignette update

- Added Spearman and Kendall correlations (argument correlationMethod)


genomation
----------

Changes in version 0.99.9:

IMPROVEMENTS AND BUG FIXES

- vignette built with knitrBootstrap

Changes in version 0.99.8:

IMPROVEMENTS AND BUG FIXES

- added the Bioinformatics citation

- vignette is converted to html format

- changed tests form test_that to RUnit

Changes in version 0.99.0.2:

IMPROVEMENTS AND BUG FIXES

- multiHeatMatrix() works correctly when common.scale=TRUE

Changes in version 0.99.0.1:

NEW FUNCTIONS AND FEATURES

- readBed() function to read bed files in to R as GRanges objects

IMPROVEMENTS AND BUG FIXES

- typo in readGeneric,readFeatureFlank,read* functions argument
  "remove.unusual" is fixed

Changes in version 0.99:

NEW FUNCTIONS AND FEATURES

- new plotting functions for visualization of ScoreMatrix and
  ScoreMatrixList multiHeatMatrix, heatMatrix, metaHeat and metaPlot

- ScoreMatrix constructor can deal with a variety of inputs, including
  BAM and wig


genomeIntervals
---------------

Changes in version 1.23.2:

- Ported version 1.22.2-1.22.3 changes

Changes in version 1.23.1:

- Ported version 1.22.1 changes

Changes in version 1.23.0:

- Bioc Version 3.1

Changes in version 1.22.3:

- Fixed a documentation cosmetic issue

Changes in version 1.22.2:

- Fixed a pretty printing issue if the formatting (inclusion of space)

Changes in version 1.22.1:

- Fixed an issue in the formatting of the coordinates in the writeGff3
  function

GenomicAlignments
-----------------

Changes in version 1.4.0:

NEW FEATURES

- All "findOverlaps" methods now support 'select' equal "last" or
  "arbitrary" (in addition to "all" and "first").

SIGNIFICANT USER-LEVEL CHANGES

- Add mapToAlignments(), pmapToAlignments(), mapFromAlignments(), and
  pmapFromAlignments() as replacements for the "mapCoords" and
  "pmapCoords" methods for GAlignments objects.

- Clarify use of 'fragments' in summarizeOverlaps() man page.

- Tweak "show" method for GAlignments objects to display a shorter
  version of long CIGARs.

- Add checks and more helpful error message for summarizeOverlaps()
  when "file does not exist"

DEPRECATED AND DEFUNCT

- Deprecated readGAlignment*FromBam() functions in favor of
  readGAlignments(), readGAlignmentPairs(), readGAlignmentsList() and
  readGappedReads().

- Deprecated "mapCoords" and "pmapCoords" methods.

- Removed Lngap(), Rngap(), introns(), and
  makeGAlignmentsListFromFeatureFragments() functions, and "ngap",
  "map", "pmap", and "splitAsListReturnedClass" methods (were defunct
  in GenomicAlignments 1.2.0).

BUG FIXES

- Fix off-by-one error when processing 'S' in query_locs_to_ref_locs().

GenomicFeatures
---------------

Changes in version 1.20:

NEW FEATURES

- Add makeTxDbFromGRanges() for extracting gene, transcript, exon, and
  CDS information from a GRanges object structured as GFF3 or GTF, and
  returning that information as a TxDb object.

- TxDb objects have a new column ("tx_type" in the "transcripts" table)
  that the user can request thru the 'columns' arg of the transcripts()
  extractor. This column is populated when the user makes a TxDb object
  from Ensembl (using makeTxDbFromBiomart()) or from a GFF3/GTF file
  (using makeTxDbFromGFF()), but not yet (i.e. it's set to NA) when
  s/he makes it from a UCSC track (using makeTxDbFromUCSC()). However
  it seems that UCSC is also providing that information for some tracks
  so we're planning to have makeTxDbFromUCSC() get it from these tracks
  in the near future.  Also low-level makeTxDb() now imports the
  "tx_type" column if supplied.

- Add transcriptLengths() for extracting the transcript lengths from a
  TxDb object. It also returns the CDS and UTR lengths for each
  transcript if the user requests them.

- extractTranscriptSeqs() now works on a FaFile or GmapGenome object
  (or, more generally, on any object that supports seqinfo() and
  getSeq()).

SIGNIFICANT USER-VISIBLE CHANGES

- Renamed makeTranscriptDbFromUCSC(), makeTranscriptDbFromBiomart(),
  makeTranscriptDbFromGFF(), and makeTranscriptDb() ->
  makeTxDbFromUCSC(), makeTxDbFromBiomart(), makeTxDbFromGFF(), and
  makeTxDb(). Old names still work but are deprecated.

- Many changes and improvements to makeTxDbFromGFF(): - Re-implemented
  it on top of makeTxDbFromGRanges().  - The geneID tag, if present, is
  now used to assign an external gene id to transcripts that couldn't
  otherwise be linked to a gene. This is for compatibility with some
  GFF3 files from FlyBase (see for example dmel-1000-r5.11.filtered.gff
  included in this package).  - Arguments 'exonRankAttributeName',
  'gffGeneIdAttributeName', 'useGenesAsTranscripts', and 'gffTxName',
  are not needed anymore so they are now ignored and deprecated.  -
  Deprecated 'species' arg in favor of new 'organism' arg.

- Some tweaks to makeTxDbFromBiomart(): - Drop transcripts with UTR
  anomalies with a warning instead of failing.  We've started to see
  these hopeless transcripts with the release 79 of Ensembl in the
  dmelanogaster_gene_ensembl dataset (based on FlyBase r6.02 /
  FB2014_05). With this change, the user can still make a TxDb for
  dmelanogaster_gene_ensembl but some transcripts will be dropped with
  a warning.  - BioMart data anomaly warnings and errors now show the
  first 3 problematic transcripts instead of 6.

- 'gene_id' metadata column returned by genes() is now a character
  vector instead of a CharacterList object.

- Use # prefix instead of | in "show" method for TxDb objects.

DEPRECATED AND DEFUNCT

- Deprecated makeTranscriptDbFromUCSC(), makeTranscriptDbFromBiomart(),
  makeTranscriptDbFromGFF(), and makeTranscriptDb(), in favor of
  makeTxDbFromUCSC(), makeTxDbFromBiomart(), and makeTxDbFromGFF(), and
  makeTxDb().

- Deprecated species() accessor in favor of organism() on TxDb objects.

- sortExonsByRank() is now defunct (was deprecated in GenomicFeatures
  1.18)

- Removed extractTranscriptsFromGenome(), extractTranscripts(),
  determineDefaultSeqnameStyle() (were defunct in GenomicFeatures
  1.18).

BUG FIXES

- makeTxDbFromBiomart(): - Fix issue causing the download of
  'chrominfo' data frame to fail when querying the primary Ensembl mart
  (with host="ensembl.org" and biomart="ENSEMBL_MART_ENSEMBL").  - Fix
  issue with error reporting code: when some transcripts failed to pass
  the sanity checks, the error message was displaying the wrong
  transcripts. More precisely, many good transcripts were mistakenly
  added to the set of bad transcripts and included in the error
  message.

- extractTranscriptSeqs(): fix error message when internal call to
  exonsBy() fails on 'transcripts'.

GenomicFiles
------------

Changes in version 1.4.0:

NEW FEATURES

- Add reduceFiles() and reduceRanges()

- Add 'algorithm' argument to summarizeOverlaps methods

- Add REDUCEsampler() from Martin

MODIFICATIONS

- Deprecate *FileViews classes

- Modify show() for GenomicFiles class

- Add 'Chunking' section to vignette

- Update vignette figures

- Change REDUCE default from `+` to `c` for reduceByYield()

BUG FIXES

- Bug fix for reduceByRange,GenomicFiles-method

GenomicRanges
-------------

Changes in version 1.20.0:

NEW FEATURES

- Add coercion methods to go back and forth between ExpressionSet and
  SummarizedExperiment.

- Add 'assayNames', 'assayNames<-' for SummarizedExperiment

- assays() supports arrays of up to 4 dimensions.

- Add GNCList() for preprocessing a GenomicRanges object into a GNCList
  object that can be used for fast overlap seach with findOverlaps().
  GNCList() is a replacement for GIntervalTree() that uses Nested
  Containment Lists instead of interval trees. Unlike GIntervalTree(),
  GNCList() supports preprocessing of a GenomicRanges object with
  ranges located on a circular sequence. For a one time use, it's not
  advised to explicitely preprocess the input. This is because
  findOverlaps() or countOverlaps() will take care of it and do a
  better job at it (that is, they preprocess only what's needed when
  it's needed and release memory as they go).

- All "findOverlaps" methods now support 'select' equal "last" or
  "arbitrary" (in addition to "all" and "first").

- Add absoluteRanges() and relativeRanges() to transform back and forth
  between absolute and relative genomic ranges.

SIGNIFICANT USER-LEVEL CHANGES

- Renamed 'rowData' and 'rowData<-' -> 'rowRanges', 'rowRanges<-'. Old
  names still work but are deprecated.

- Some improvements to makeGRangesFromDataFrame(): - Improve internal
  logic used for finding the GRanges columns in the input. - If
  'seqinfo' is not supplied, the seqlevels are now ordered according to
  the output of GenomeInfoDb::rankSeqlevels(). - Now an attempt is made
  to turn 'df' into a data frame (with 'as.data.frame(df)') if it's not
  a data frame or a DataFrame object.

- The GRanges() constructor now propagates the metadata cols that are
  on 'ranges' if no metadata cols are explicitly passed to the
  constructor.

DEPRECATED AND DEFUNCT

- Deprecated 'rowData' and 'rowData<-' in favor of 'rowRanges' and
  'rowRanges<-'.

- Deprecated mapCoords() and pmapCoords(). They're replaced by
  mapToTranscripts() and pmapToTranscripts() from the GenomicFeatures
  package and mapToAlignments() and pmapToAlignments() from the
  GenomicAlignments package.

- Deprecated GIntervalTree objects.

- Removed "map" and "splitAsListReturnedClass" methods (were defunct in
  GenomicRanges 1.18.0).

- Removed makeSeqnameIds() (was defunct in GenomicRanges 1.18.0).

- Removed 'with.mapping' argunment from "reduce" methods (was defunct
  in GenomicRanges 1.18.0).

BUG FIXES

- Fix 'findOverlaps(..., type="start")' on GRangesList objects which
  has been broken for years.

- Fix self overlap search on a GRanges object when 'ignore.strand=TRUE'
  (i.e. 'findOverlaps(gr, ignore.strand=TRUE)').

GenomicTuples
-------------

Changes in version 1.1.14:

- Switched to a html vignette from a pdf one.

genoset
-------

Changes in version 1.21.10:

NEW FEATURES

- calcGC gets a 'bases' argument and handles the presence or absence of
  'chr' prefixes in a safer manner.

GGtools
-------

Changes in version 5.3:

- new VCF extract for demonstrating cisAssoc VCF-SummarizedExperiment
  analysis

- GGtools kept for legacy/migration purposes only; use gQTLBase and
  gQTLstats instead

ggtree
------

Changes in version 0.99.28:

- update vignette with floating table of content <2015-04-08, Wed>

Changes in version 0.99.27:

- bug fixed, see https://github.com/GuangchuangYu/ggtree/issues/4
  <2015-03-07, Tue>

Changes in version 0.99.26:

- update geom_tiplab <2015-03-31, Tue>

- update plot method of beast <2015-03-17, Tue>

Changes in version 0.99.25:

- implement groupClade <2015-03-13, Fri>

Changes in version 0.99.24:

- use "round" segment end, look very better <2015-03-12, Thu>

- update vignett <2015-03-11, Wed>

Changes in version 0.99.23:

- mv geom_hilight to hilight <2015-03-11, Wed>

- mv geom_phylopic to phylopic <2015-03-11, Wed>

- implement collapse and expand for collapse and expand a selected
  clade <2015-03-11, Wed>

Changes in version 0.99.22:

- remove quote in beast tip/node labels <2015-03-10, Tue>

Changes in version 0.99.21:

- fixed downloading png file in Windows platform, should explicitly
  setting mode="wb". <2015-03-03, Tue>

Changes in version 0.99.19:

- for time scale tree inferred by BEAST, now user can use
  time_scale=TRUE parameter in ggtree function <2015-02-12, Thu>

Changes in version 0.99.18:

- bug fixed in reorder the labels in gplot.heatmap <2015-02-12, Thu>

Changes in version 0.99.17:

- add angle and branch.y variable in cladogram layout <2015-02-10, Tue>

Changes in version 0.99.16:

- correct typo in vignette <2015-02-10, Tue>

Changes in version 0.99.15:

- fully support of replace operator, %<% <2015-02-09, Mon>

Changes in version 0.99.14:

- add example in groupOTU for adding legend manually <2015-02-09, Mon>.

Changes in version 0.99.13:

- two dimensional tree <2015-02-06, Fri>

Changes in version 0.99.12:

- update vignette <2015-02-04, Wed>

- gzoom methods that supports all tree objects <2015-02-04, Wed>

- geom_hilight layer for highlighting clade <2015-02-04, Wed>

Changes in version 0.99.11:

- add scale_color to support colored lines and text based on numerical
  values and update vignette <2015-02-04, Wed>

- revised groupOTU that output index can be used in geom_text and
  update vignette <2015-02-04, Wed>

Changes in version 0.99.10:

- support y scale by category variable <2015-02-03, Tue>

- support order nodes by yscale <2015-02-03, Tue>

Changes in version 0.99.9:

- update vignette <2015-02-02, Mon>

Changes in version 0.99.8:

- add get.phylopic function to read the online phylo pic and convert it
  to grob object, which can be use to annotate ggplot figure using
  annotation_custom <2015-01-30, Fri>

Changes in version 0.99.7:

- add angle information for 'fan' & 'unrooted' layout <2015-01-29, Thu>

Changes in version 0.99.6:

- read.beast now supports support values of sets such as {x, y, z}
  <2015-01-19, Mon>

- now read.beast supports characters in support values <2015-01-18,
  Sun>

- add example of gzoom and groupOTU in vignette <2015-01-14, Wed>

- implement groupOTU methods <2015-01-14, Wed>

- export get.offspring.tip <2015-01-14, Wed>

Changes in version 0.99.5:

- move ape and ggplot2 from Depends to Imports <2015-01-12, Mon>

- get.tipseq method for paml_rst and codeml object <2015-01-08, Thu>

- add gzoom function, similar to zoom function in ape <2015-01-07, Wed>

- add examples in man pages of %<% and %<+% operators <2015-01-06, Tue>

- remove <<- and update vignette <2015-01-06, Tue>

- update vignette and use BibTex and CSL for references <2015-01-05,
  Mon>

- update cladogram layout <2015-01-05, Mon>

- read.baseml function and update vignette with baseml example
  <2015-01-04, Sun>

- plot method for hyphy and hyphy example in vignette <2015-01-04, Sun>

- merge all vignettes to ggtree vignette <2015-01-04, Sun>

Changes in version 0.99.4:

- ggtree now support branch.length = "none" to only draw tree topology
  <2015-01-03, Sat>

- get.subs method for hyphy object <2015-01-03, Sat>

- show, get.tree and get.subs methods of hyphy <2015-01-02, Fri>

- export read.hyphy <2015-01-02, Fri>

- export hyphy class <2015-01-01, Thu>

- plot method for beast class and get.tree method for codeml class
  <2014-12-31, Wed>

- show, get.fields, get.subs and plot methods for codeml class
  <2014-12-30, Tue>

- plot method for paml_rst class <2014-12-30, Tue>

- get.subs, method for paml_rst class <2014-12-30, Tue>

- show, plot, get.tree and fority methods for codeml_mlc class
  <2014-12-29, Mon>

- codeml class <2014-12-29, Mon>

- hyphy class and read.hyphy prototype <2014-12-27, Sat>

- update man file and add example file of beast output <2014-12-26,
  Fri>

- get.tree and get.fileds methods of beast class <2014-12-26, Fri>

- read.beast <2014-12-26, Fri>

- beast class and show method <2014-12-26, Fri>

- coplot prototype<2014-12-24, Wed>

- parse translation matrix in beast nexus <2014-12-24, Wed>

- extract beast stats info <2014-12-23, Tue>

Changes in version 0.99.3:

- gplot function that can view a tree and an associated matrix
  simultaneously <2014-12-22, Mon>

- modified vignette to show based on branch position and break the
  evolution distance scale <2014-12-22, Mon>

- label and annotation can be put based on branch. <2014-12-22, Mon>

- write.jplace and fully supports of jplace by ggtree. <2014-12-21,
  Sun>

- support unrooted layout in ggplot. <2014-12-21, Sun>

- support fan, radial, dendrogram layout in geom_tree. <2014-12-21,
  Sun>

Changes in version 0.99.2:

- layout of unrooted tree, implemented equal-angle algorithm that
  described in Chapter 34 of 'Inferring Phylogenies' (page 578-580)
  <2014-12-20, Sat>

- add layout parameter in ggtree and geom_tree, now supports phylogram
  and cladogram <2014-12-20, Sat>

- %<+% function for inserting annotation data to a tree view
  <2014-12-20, Sat>

- update ggtree-treeAnnotation vignette <2014-12-20, Sat>

Changes in version 0.99.1:

- %<% function for updating tree view with a new tree <2014-12-19, Fri>

- add examples in man files <2014-12-19, Fri>


GOexpress
---------

Changes in version 1.1.12:

NEW FEATURES

- Summarisation function of scores and rank from feature-level to
  ontology-level also affects the calculation of P-values by pValue_GO.
  The default behaviour is to use the same function as used in the call
  to GO_analyse, stored in the result object.

Changes in version 1.1.11:

NEW FEATURES

- Summarisation function of scores and rank from feature-level to
  ontology-level can be overriden from the default "mean" to any
  function specified by the user. The function need to support a list
  of numeric values as an input, and return a unique numeric value as
  an output.

Changes in version 1.1.10:

BUG FIX

- margins argument in heatmap_GO() was not used in the call to
  heatmap.2()

Changes in version 1.1.9:

TYPOS

- Missing space character in error message.

Changes in version 1.1.8:

TYPOS

- Simone Coughlan, not Simone Coughland. Apologies!

Changes in version 1.1.7:

BUG FIXES

- Changed the rank of annotated features absent from the given
  ExpressionSet to (number of features in ExpressionSet + 1) instead of
  max(rank) + 1 to match the documentation, and desired behaviour. This
  follows the logic of giving the same rank R to N tied-scored
  features, while the next feature receives a rank of R + N.

Changes in version 1.1.6:

BUG FIXES

- overlap_GO() was crashing for 3-group Venn diagrams, except if the
  VennDiagram was loaded manually loaded in the workspace using
  libray(VennDiagram). The function will now run seemlessly without
  that manual step, as loading GOexpress will immediately load
  VennDiagram in the workspace (stated as a dependency in the
  DESCRIPTION file).

NEW FEATURES

- New function pValue_GO() allows calculation of P-value for each
  ontology using permutation of genes labels. This allows users to
  estimate the chance of seeing a GO term reach a particular rank (or
  score). Features a fancy progres bar shamelessly adapted from
  StackOverflow.

- heatmap_GO now semi-autmoatically resizes the bottom and right
  margins to accomodate large gene and sample labels, respectively. The
  user may control those margins using the "margins" argument of the
  function.

- heatmap_GO default call now shows the gene feature identifier for
  those missing an annotated gene name, when gene names are requested
  (also the default).

- A rank.by slot is now created by the GO_analyse() function in the
  result object to state the metric used to order the result tables.

- a filters.GO slot stating the filters and cutoffs applied to the
  result object is now created or updated by successive uses of the
  subset_scores() function. Warnings and notes are displayed if
  conflicting filters and cutoffs are applied on a previously filtered
  result object.

- rerank() function now supports re-ordering by P-value. Note that this
  is only applicable to the output of the pValue_GO() function
  mentioned above.

- rerank() function now updates the rank.by slot of the result object
  to state the current ordering metric.

- subset_scores() function now allows filtering by P-value. Note that
  this is only applicable to the output of the pValue_GO() function
  mentioned above.

- Backward compatibility with Ensembl annotation releases 75 and
  earlier, which used 'external_gene_id', which was renamed to
  'external_gene_name' in releases 76 and later.

- table_genes() function defaults to sorting genes by decreasing score
  (equivalent to increasing rank). Gene feature name or identifier are
  supported alternative filters for sorting.

UPDATED FEATURES

- Allow user to override row_labels in heatmap_GO. This way, the
  color-coding of the sample can be kept, while better description of
  the samples can be used to label them, instead of the phenodata
  values.

- In heatmap_GO(), if the labRow argument is of length 1, it is assumed
  to be the name of a column in the phenoData slot. Useful to re-label
  subsetted ExpressionSet objects.

GENERAL UPDATES

- Updated the AlvMac training dataset to include 'RPL36A' an example of
  multiple Ensembl gene identifier annotated to the same gene name.

- Updated the AlvMac example custom annotations to match the updated
  dataset.

- Updated the example AlvMac_results to match the updated dataset.

- Set the random seed prior to running the GO_analyse() example in the
  vignette. Hopefully, this should allow reproducible testing by the
  users.

- In User Guide, load package before loading the attached data.

- In User Guide, new sections and examples dealing with the
  re-labelling of heatmap samples, the use of P-values, the re-ranking
  and subsetting of results using P-values. New sub-sections for
  clarity. Emphasis on the use and generation of local annotation,
  rather than use of current online Ensembl annotation release.

- No more code connecting to the Ensembl server in any the help files
  and User Guide.

- Help pages examples with more consistent indentation of code.

Changes in version 1.1.5:

BUG FIXES

- Forgot to commit image file of shiny screenshot in release 1.1.5

Changes in version 1.1.4:

NEW FEATURES

- Custom annotations can be provided to the GO_analyse function using
  three new arguments: "GO_genes", "all_GO", and "all_genes". See below
  for individual description.

- The GO_genes argument allows the user to provide associations between
  feature identifiers in the ExpressionSet and gene ontology
  identifiers. This will skip all calls to the Ensembl BioMart server.
  Consequently, we recommend the use of the "all_GO" and "all_genes"
  arguments whenever "GO_genes" is used. Otherwise, some downstream
  functions may not work. For instance, the expression_plot,
  expression_profiles, and heatmap_GO function can be used to generate
  plots,although lacking gene and gene ontology names. See below for
  more details.

- The all_GO argument allows the user to provide the name and namespace
  ("biological_process", "molecular_function", or "cellular_component")
  corresponding to gene ontology identifiers. This enables subsequent
  filtering of result tables by namespace, and annotation of heatmaps
  with the name of the gene ontology.

- The all_genes argument allows the user to provide the name and
  description corresponding to feature identifiers in the
  ExpressionSet. This enables annotation of expression plots with the
  gene name associated with the feature.

UPDATED FEATURES

- In the GO_analyse function, the eSet argument is now formally checked
  to be of class ExpressionSet prior to any calculation. If not, the
  function returns an appropriate error message.

- The error messages caused when the user gives a name that does not
  exist in the phenoData slot of the ExpressionSet were updated to use
  the word "column", rather than "factor".

GENERAL UPDATES

- Updated the help page for function GO_analyse to describe the new
  features described above and provide a code example.

- Used the BiocStyle package to format the vignette.

- Added a new section in the vignette to document the new features
  described above.

- Added a new section in the vignette to mention the creation of shiny
  applications using the output of GOexpress. Included a screenshot of
  a shiny application developed from the original AlvMac dataset.

- Added a new section in the vignette to highlight the availability of
  a "subset" argument to avoid the need for additional ExpressionSet
  objects to analyse or visualise subsets of samples.

Changes in version 1.1.3:

NEW FEATURES

- Included private function used to generate the prefix2dataset table.

- Included private function used to generate the microarray2dataset
  table.

- GO_analyse states the number of gene features in the given dataset
  that were found in the BioMart dataset. This allows the user to
  interrupt the script if a suspiciously low number of mapped features
  suggests an incorrect BioMart dataset was used.

UPDATED FEATURES

- Updated the prefix2dataset table. Three more species ("Chlorocebus
  sabaeus", "Papio anubis" , and "Poecilia formosa"), and two more
  columns ("species" and "sample").

- Updated the microarray2dataset table. Probeset patterns are now hard
  coded in the package, but dynamically defined as unique to a platform
  (or not) by the function building the microarray2dataset table.
  Species without microarray platforms were also included in the table
  with NA values. Code was added in the GO_analyse method to ignore
  those before trying to resolve the origin of the expression data, if
  not specified by the user.

- Renamed column "prefix" to "pattern" in microarray2dataset table.

- Disabled manual check of C. elegans and D. melanogaster, as the
  automated detection is doing the exact same thing. Only S. cerevisiae
  requires a manual check of the "Y" prefix instead of the
  automatically extracted "Y[:LETTERS:]{2}" pattern.

- Methods expression_plot_symbol and expression_profiles_symbol now
  return the list of gene feature identifiers instead of NULL when
  multiple plots are produced. That may be confusing as they return the
  list of close matches when an invalid gene symbol was given, and they
  return the ggplot when obly one plot is produced.

GENERAL UPDATES

- Code layout. Lines of code over 80 characters were split around
  brackets and commas to use the built-in indentation defaulted to 4
  space characters. Plus, this made the code more readable in many
  cases.

Changes in version 1.1.2:

GENERAL UPDATES

- DESCRIPTION file minor correction: Inappropriate "Metagenomics"
  biocView removed.

Changes in version 1.1.1:

BUG FIXES

- Ensembl BioMart has changed the column 'external_gene_id' to
  'external_gene_name'. Renamed my biomaRt queries accordingly to
  prevent GO_analyse() from crashing during the analysis step.

- Updated the AlvMac_results variable containing sample results
  annotated with the deprecated 'external_gene_id'. It now includes the
  new 'external_gene_name'.

GENERAL UPDATES

- DESCRIPTION file minor correction.

- User Guide's chunk about installation using biocLite is not evaluated
  anymore, as advised by Dan Tenenbaum.

GoogleGenomics
--------------

Changes in version 1.0.0:

NEW FEATURES

- Authentication support for native apps, service accounts, and API
  keys.

- Support for retreiving Reads and Variants.

- Simple converters to GAlignments for Reads, and GRanges and VRanges
  for Variants.

GOSemSim
--------

Changes in version 1.25.5:

- update information content files <2015-03-12, Thu>

Changes in version 1.25.4:

- update vignette and add DOSE citation <2015-02-13, Fri>

Changes in version 1.25.3:

- add doi in CITATION <2015-01-28, Wed>

- update vignette using BiocStyle <2015-01-26, Mon>

Changes in version 1.25.2:

- add BugReports URL <2014-12-17, Wed>

Changes in version 1.25.1:

- import Rcpp <2014-10-23, Thu>

graphite
--------

Changes in version 1.13.2 (2015-04-02):

- Pathways for multiple species: see pathwayDatabases() output.

- Native pathway IDs.

- Removed SPIKE database.

- Deprecated objects: biocarta, humancyc, kegg, nci, panther, reactome.

- Deprecated functions: runClipperMulti, runDEGraphMulti,
  runTopologyGSAMulti.

gtrellis
--------

Changes in version 0.99.3:

- modified infinite loop when making asistant ticks

Changes in version 0.99.2:

- change `gtrellis_info` to `gtrellis_show_index`

Gviz
----

Changes in version 1.11.0:

NEW FEATURES

- BiomartGeneRegionTracks can now be created based on a gene symbol,
  Ensembl trancript id, Ensebml gene id or ENTREZ gene id.

SIGNIFICANT USER-VISIBLE CHANGES

- Streaming behaviour for BiomartGeneRegionTracks. Already queried
  regions are cached, and new data is fetched from Biomart on demand.

- Proper handling of UCSC genome identifiers in
  BiomartGeneRegionTracks. Automated mapping to Ensembl genome versions
  and Biomart archives.

BUG FIXES

- A number of significant fixes.

GWASTools
---------

Changes in version 1.13.27:

- Added "snpID" and "scanID" arguments to getGenotypeSelection.

Changes in version 1.13.26:

- Added getScanAnnotation, getSnpAnnotation accessors for GenotypeData
  objects.

Changes in version 1.13.24:

- Added data tables for genome build 38: centromeres.hg38.RData,
  pseudoautosomal.hg38.RData, HLA.hg38.RData, pcaSnpFilters.hg38.RData.

Changes in version 1.13.23:

- convertGdsNcdf works for transposed (sample x snp) genotype files.

Changes in version 1.13.22:

- Removed "outfile" arguments from batchChisqTest, batchFisherTest, and
  mendelErr. Saving output to a file should happen outside the function
  calls.

- batchChisqTest and batchFisherTest have snp.include arguments to run
  on individual SNPs. Using batchFisherTest with this argument is
  recommended to replace the deprecated assocTestFisherExact.

Changes in version 1.13.21:

- assocRegression replaces assocTestRegression. Only one model is
  allowed per function call.

- assocCoxPH replaces assocTestCPH. Output format is now similar to
  assocRegression.

- exactHWE replaces gwasExactHWE.

- assocRegression, assocCoxPH, and exactHWE include the option to
  select blocks of SNPs by index for easier parallelization.

- scan.chromosome.filter is no longer an option; use
  setMissingGenotypes to filter data prior to running other functions.

Changes in version 1.13.9:

- Add use.names argument to getGenotype and getGenotype selection

- Add order=c("file", "selection") argument to getGenotypeSelection

- duplicateDiscordanceAcrossDatasets and dupDosageCorAcrossDatasets
  will not match on unmapped SNPs

Changes in version 1.13.8:

- Add drop=TRUE argument to getVariable, etc.

Changes in version 1.13.7:

- Add dupDosageCorAcrossDatasets.

- Add getGenotypeSelection method to MatrixGenotypeReader.

Changes in version 1.13.6:

- Add getGenotypeSelection to select non-continguous SNPs and scans
  from GDS files.

Changes in version 1.13.3:

- Bug fix in imputedDosageFile - for IMPUTE2, include columns from
  .samples file in output scan annotation.

HDTD
----

Changes in version 1.1.2 (2015-01-09):

- Updated DESCRIPTION FILE.

Changes in version 1.1.1 (2014-11-04):

- Added the function transposedata.

- Updated the vignette file.

- Updated the CITATION file.

- Updated output of the core fuctions.

- Updated documentation.

hiAnnotator
-----------

Changes in version 1.1.1:

- added allSubjectCols & overlapType params to getSitesInFeature

HIBAG
-----

Changes in version 1.4.0:

- The version number was bumped for the Bioconductor release version

Changes in version 1.3.0-1.3.2:

NEW FEATURES

- support the human genome "hg20"

- add a new function `hlaGDS2Geno()` to support SNP GDS files (in the
  R/Bioconductor package SNPRelate)

- `hlaReport()` outputs text with markdown format

SIGNIFICANT USER-VISIBLE CHANGES

- optimize the calculation of hamming distance using SSE2 and hardware
  POPCNT instructions if available

- hardware POPCNT: 2.4x speedup for large-scale data, compared to the
  implementation in v1.2.4

- SSE2 popcount implementation: 1.5x speedup for large-scale data,
  compared to the implementation in v1.2.4

BUG FIXES

- bug fixes on big-endian machines (like Solaris SPARC, Apple PowerPC)

- minor fix on random sampling from discrete uniform distribution

- bug fix if 'requireNamespace("HIBAG")' instead of 'require(HIBAG)' is
  called from other packages


hiReadsProcessor
----------------

Changes in version 1.1.3:

- Modified all bplapply calls

- Moved BiocGenerics and rSFFreader to Imports

HiTC
----

Changes in version 1.11.4:

NEW FEATURES

- New directionalyIndex function as a first step of TADs detection

SIGNIFICANT USER-VISIBLE CHANGES

- new forcePairwise method for HTClist or HTCexp objects. The methods
  reverse the forceSymmetrical method. Useful for plotting function.

BUG FIXES

- Bug in mapC for Pairwise Hi-C map

Changes in version 1.11.3:

BUG FIXES

- Bug in setEnvDisplay report by S. Thomas. Change decimal to 5 to
  avoid layout of size 0 for small chromosomes

Changes in version 1.11.1:

NEW FEATURES

- Changes in forceSymmetric method on HTCexp object. The default is now
  to sum up both upper and lower maps.

BUG FIXES

- Bug in import.my5C function

hpar
----

Changes in version 1.9.1:

- Updated to HPA version 13 <2014-11-15 Sat>

- Use BiocStyle for vignette <2014-11-15 Sat>

HTSFilter
---------

Changes in version 1.7.1:

- -- Minor bug fix for integration in edgeR pipeline -- Minor updates
  in vignette and documentation -- Removed IRanges and GenomicRanges
  from imports (hard-coded the mcols and colData functions for use in
  DESeq2 pipeline instead) -- Add BiocStyle to Suggests for vignette
  style

illuminaio
----------

Changes in version 0.9.1 (2015-02-25):

- Modified code for reading encrypted idat files to cope with VeraCode
  data.

Changes in version 0.9.0 (2014-10-13):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.1 for R (>= 3.2.0).

immunoClust
-----------

1.0.0: The first version caintains basically the functios and routines
        usesd to obtain the results of "Soerensen, T., Baumgart, S.,
        Durek, P., Gruetzkau, A. and Haeupl, T.  immunoClust - an
        automated analysis pipeline for the identification of
        immunophenotypic signatures in high-dimensional cytometric
        datasets.  Cytometry A (accepted)." CHANGES: * The code was
        cleaned and modified in the C-binding calls to make it runable
        on R 3.1.2.  * A bug in the cell.hclust function was fixed,
        which does not effect the general results but lead to minor
        differences in concrete numbers.

IMPCdata
--------

Changes in version 1.0.1:

FEATURES

- Package allows systematically explore the IMPC dataset's multiple
  dimensions until the correct combination of filters has been
  selected.

- Package helps to obtain datasets from IMPC database.

- The data retrieved from IMPCdata package can be directly used by
  PhenStat -- an R package that encapsulates the IMPC statistical
  pipeline, available at <URL:
  http://www.bioconductor.org/packages/release/bioc/html/PhenStat.html>.

InPAS
-----

Changes in version 0.99.7:

NEW FEATURES

- compensation by last CDS but not by ratio.

Changes in version 0.99.6:

BUG FIXES

- fix the bug in utr3UsageEstimation in calculating of PDUI

Changes in version 0.99.5:

BUG FIXES

- fix the bug in utr3Annotation when multiple transcripts sharing same
  transcript name by Gviz::GeneRegionTrack

Changes in version 0.99.4:

BUG FIXES

- fix the bug in documentation of utr3Annotation

Changes in version 0.99.3:

BUG FIXES

- fix the output ranges from short form to whole 3UTR

Changes in version 0.99.2:

BUG FIXES

- clear warnings when run check.  Warning: multiple methods tables
  found for 'species' Warning: multiple methods tables found for
  ‘organism’

Changes in version 0.99.1:

BUG FIXES

- errors when run BioCheck. import all the package should be import.

intansv
-------

Changes in version 1.7.1:

Notes

- Fixed a bug in svAnnotation.


IRanges
-------

Changes in version 2.2.0:

NEW FEATURES

- Add NCList() and NCLists() for preprocessing a Ranges or RangesList
  object into an NCList or NCLists object that can be used for fast
  overlap search with findOverlaps(). NCList() and NCLists() are
  replacements for IntervalTree() and IntervalForest() that use Nested
  Containment Lists instead of interval trees. For a one time use, it's
  not advised to explicitely preprocess the input. This is because
  findOverlaps() or countOverlaps() will take care of it and do a
  better job at it (that is, they preprocess only what's needed when
  it's needed and release memory as they go).

- Add coercion methods from Hits to CompressedIntegerList, to
  PartitioningByEnd, and to Partitioning.

SIGNIFICANT USER-VISIBLE CHANGES

- The code behind overlap-based operations like findOverlaps(),
  countOverlaps(), subsetByOverlaps(), summarizeOverlaps(), nearest(),
  etc... was refactored and improved. Some highlights on what has
  changed: - The underlying code used for finding/counting overlaps is
  now based on the Nested Containment List algorithm by Alexander V.
  Alekseyenko and Christopher J. Lee.  - The old algorithm based on
  interval trees is still available (but deprecated). The 'algorithm'
  argument was added to most overlap-based operations to let the user
  choose between the new (algorithm="nclist", the default) and the old
  (algorithm="intervaltree") algorithm.  - With the new algorithm, the
  hits returned by findOverlaps() are not fully ordered (i.e. ordered
  by queryHits and subject Hits) anymore, but only partially ordered
  (i.e. ordered by queryHits only). Other than that, and except for the
  3 particular situations mentioned below, choosing one or the other
  doesn't affect the output, only performance.  - Either the query or
  subject can be preprocessed with NCList() for a Ranges object
  (replacement for IntervalTree()), NCLists() for a RangesList object
  (replacement for IntervalForest()), and GNCList() for a GenomicRanges
  object (replacement for GIntervalTree()).  However, for a one time
  use, it's not advised to explicitely preprocess the input. This is
  because findOverlaps() or countOverlaps() will take care of it and do
  a better job at it (that is, they preprocess only what's needed when
  it's needed and release memory as they go).  - With the new
  algorithm, countOverlaps() on Ranges or GenomicRanges objects doesn't
  call findOverlaps() to collect all the hits in a growing Hits object
  and count them only at the end. Instead the counting happens at the C
  level and the hits are not kept. This reduces memory usage
  considerably when there is a lot of hits.  - When 'minoverlap=0',
  zero-width ranges are interpreted as insertion points and are
  considered to overlap with ranges that contain them.  This is the 1st
  situation where using 'algorithm="nclist"' or
  'algorithm="intervaltree"' produces different output.  - When using
  'select="arbitrary"', the new algorithm will generally not select the
  same hits as the old algorithm. This is the 2nd situation where using
  'algorithm="nclist"' or 'algorithm="intervaltree"' produces different
  output.  - When using the old interval tree algorithm, 'maxgap' has a
  special meaning if 'type' is "start", "end", or "within". This is not
  yet the case with the new algorithm. That feature seems somewhat
  useful though so maybe the new algorithm should also support it?
  Anyway, this is the 3rd situation where using 'algorithm="nclist"' or
  'algorithm="intervaltree"' produces different output.  - Objects
  preprocessed with NCList(), NCLists(), and GNCList() are
  serializable.

- The RleViewsList() constructor function now reorders its 'rleList'
  argument so that its names match the names on the 'rangesList'
  argument.

- Minor changes to breakInChunks(): - Add 'nchunk' arg.  - Now returns
  a PartitioningByEnd instead of a PartitioningByWidth object.  - Now
  accepts 'chunksize' of 0 if 'totalsize' is 0.

- 300x speedup or more when doing unique() on a CompressedRleList
  object.

- 20x speedup or more when doing unlist() on a SimpleRleList object.

- Moved the RleTricks.Rnw vignette to the S4Vectors package.

DEPRECATED AND DEFUNCT

- Deprecated mapCoords() and pmapCoords(). They're replaced by
  mapToTranscripts() and pmapToTranscripts() from the GenomicFeatures
  package and mapToAlignments() and pmapToAlignments() from the
  GenomicAlignments package.

- Deprecated IntervalTree and IntervalForest objects.

- seqapply(), seqby(), seqsplit(), etc are now defunct (were deprecated
  in IRanges 2.0.0).

- Removed map(), pmap(), and splitAsListReturnedClass() (were defunct
  in IRanges 2.0.0).

- Removed 'with.mapping' argunment from reduce() methods (was defunct
  in IRanges 2.0.0).

BUG FIXES

- findOverlaps,Vector,missing method now accepts extra arguments via
  ...  so for example one can specify 'ignore.strand=TRUE' when calling
  it on a GRanges object (before that, 'findOverlaps(gr,
  ignore.strand=TRUE)' would fail).

- PartitioningByEnd() and PartitioningByWidth() constructors now check
  that, when 'x' is an integer vector, it cannot contain NAs or
  negative values.

isobar
------

Changes in version 1.13.0:

- calcCumulativeProbXGreaterThanY() contributed by A. Stukalov, for
  calculating cummulative p-values with replicates

- Bioconductor version increment

kebabs
------

Changes in version 1.1.9:

- inclusion of dense LIBSVM 3.20 for dense kernel matrix support to
  provide a reliable way for training with kernel matrices

- new accessors folds and performance for CrossValidationResult

- removed fold performance from show of CV result

- adaptions for user defined sequence kernel with new export
  isUserDefined, example in inst/examples/UserDefinedKernel

- correction of errors with position offset for position specific
  kernels

- computation of AUC via trapezoidal rule

- changes for auto mode in CV, grid search, model selection

- check for non-negative mixing coefficients in spectrum and gappy pair
  kernel

- build warnings on Windows removed

- added definition of performance parameters for binary and multiclass
  classification to vignette

- update of citation file and reference section in help pages

Changes in version 1.1.8:

- new accessors selGridRow, selGridCol and fullModel for class
  ModelSelectionResult

- change of naming of feature weights because of change in LiblineaR
  1.94-2

- GCC warnings in Linux removed

Changes in version 1.1.7:

- change in LiblineaR - upgrade to LIBLINEAR 1.94 in function LiblineaR
  the parameter labels was renamed to target

- correction in model selection for performance parameters

- error correction of vector length overflow in sparse explicit
  representation for very large number of sequences in spectrum, gappy
  pair and motif kernel

- error correction for AUC in cross validation

- minor changes in help pages

- minor changes in vignette

Changes in version 1.1.6:

- error correction for training with position specific kernel and
  computation of feature weights

- error correction in coercion of kernel to character for distance
  weighting

- error correction in spectrum, gappy pair and motif kernel for kernel
  matrix - last feature was missing in kernel value in rare situations

- correction of Windows build problem in linearKernel

- build warnings on Windows removed

- minor changes in help pages

- minor changes in vignette

Changes in version 1.1.5:

- new method heatmap to display heatmap of prediction profiles

- extension of function linearKernel to optionally return a sparse
  kernel matrix

- correction of computation of feature weights for LiblineaR with more
  than 3 classes

- new accessor SVindex for class KBModel

- correction in subsetting of sparse explicit representation for head /
  tail

- error correction in subsetting of prediction profile

- error correction in mismatch kernel

- check uniqueness of motifs in motif kernel

- minor changes in help pages

- change name of vignette Rnw to lowercase

- minor changes in vignette

Changes in version 1.1.4:

- added two help pages

Changes in version 1.1.3:

- fix to adapt for changed Biostrings/S4Vectors API

Changes in version 1.1.2:

- minor C code changes for mismatch kernel

- correction of MCC

- new class ROCData and new function computeROCandAUC for binary
  classification added

- new plot function for ROCData to plot ROC for binary classes

- AUC as additional performance parameter in cross validation and as
  performance objective in grid search

Changes in version 1.1.1:

- correction for cross validation with factor label

- correction for storing prob model in kebabs model for kernlab

- removal of clang warnings for unused functions

Changes in version 1.1.0:

- first devel version created from first official release with version
  1.0.0

KEGGprofile
-----------

Changes in version 1.9.3:

- Fix a bug in example of convertId function, which was caused by
  update of BioMart database. Thanks Dan.

Changes in version 1.9.1:

- Fix a bug when pathway_min is set to 1. Thanks Norbert Auer.


limma
-----

Changes in version 3.24.0:

- Limma review article (Ritchie et al, Nucleic Acids Research, 2015) is
  now published, and becomes the primary citation for limma. References
  in User's Guide, help files and CITATION have been updated to refer
  to this article.

- New function plotMD() for mean-difference plots.  This will take over
  the functionality of the plotMA() function.

- mdplot() has new arguments 'columns' and 'main'.

- Arguments 'pch', 'cex' and 'col' of plotWithHighlights() renamed to
  'hi.pch', 'hi.cex' and 'hi.col'.  Improved help page for
  plotWithHighlights().

- plot() methods removed for EList, MAList and MArrayLM objects.  Use
  plotMD() instead.

- New function fry() provides a fast version of mroast() when there is
  little heteroscedasticity between genes.

- Minor change to mroast() output so that FDR values are never smaller
  than p-values, even when mid-p-values are used.

- camera(), roast(), mroast() and romer() now give a warning if any
  arguments found in ... are not used.

- romer() is now an S3 generic function. Speed improvements for romer()
  when the "mean" or "floormean" set statistics are used.

- topGO() can now sort gene ontology results by several p-value columns
  at once, using the minimum of the columns.  The default is now to
  sort results using all the p-value columns found in the results
  instead of just by the first column.

- goana() has new arguments 'plot' and 'covariate'. The plot argument
  allows users to view any estimated trend in the genewise significance
  values. The help files goana.Rd and goana.MArrayLM.Rd are
  consolidated into one file.

- plotDensities() has a new argument 'legend', allowing user to
  reposition the legend.

- voomWithQualityWeights() has a new argument 'col', allowing the
  barplot of sample weights to be colour-coded.

- Improvements and bug fixes to plotMDS(). plotMDS now stores
  'axislabel', depending on the method used to compute distances, and
  uses this to make more informative axis labels. plotMDS now traps
  problems when the distances are all zero or if one or more of the
  dimensions are degenerate. Also, a bug when gene.selection="common"
  and top=nrow(x) has been fixed.

- The arguments for predFCm() have been renamed.  Old argument 'VarRel'
  replaced by logical argument 'var.indep.of.fc'.  Old argument 'prop'
  replaced by 'all.de'.  New argument 'prop.true.null.method' passes
  through method to propTrueNull().

- vennDiagram() now fills circles with specified colors.

- plotMA() now has a method for EListRaw objects.

- voomWithQualityWeights() now stores the sample weights in the output
  EList object.  The values stored in $weights are unchanged, they are
  still the combined observation and sample weights.

- nec() and neqc() now remove the background Eb values if they exist in
  the EListRaw object.

- Acknowledgements edited in User's Guide.

- fix URLs for a number of references in help pages.

- Improvements to lowess C code.

- When loading suggested packages within function calls, usages of
  require() have been converted to requireNamespace() wherever
  possible.  Functions from suggested packages are called using `::'.

- topTable() now gives an information error message if eBayes() or
  treat() have not already been run on the fitted model object.

- genas() now gives an error message when applied to fitted model for
  which the coefficient correlations are gene-specific.

- Fix a bug in the MAList method for plotDensities that was introduced
  with the use of NextMethod() in version 3.21.16.

LowMACA
-------

Changes in version 0.99.4:

- custom perl command was added to the methods setup and alignsequences
  bug fixes:

- onLoad check added for: -> if "clustalo" is not in the PATH, a
  warning is raised -> if clustalo version is not 1.2.x, a warning is
  raised -> if XML:Simple and LWP modules for perl are not installed, a
  warning is raised

- check for exit status of clustalo system call added

Changes in version 0.99.3:

- alignSequences can now run in web mode. No clustalo is necessary but
  the limit is set to 2000 sequences and a email address is required.
  In addition, perl must be installed in the system

- the newLowMACA constructor is now able to check all possible gene
  symbols, even if they are not protein coding genes or if they are not
  mapped in LowMACAAnnotation package

- add a conversion from U (selenocystein) to A in amino acid sequences
  while creating the LowMACA object ('newLowMACA' function) in order to
  fix a problem with .Trident_Score function

- modified check validity function in order to accept 'auto' in the
  bandwidth

- 'protter' method is now compatible with Windows OS

- a check was added in newLowMACA constructor to prevent more than one
  pfam ID entry by the user

Changes in version 0.99.2:

- modified 'lfm' method in order to accept different thresholds of p
  value or q value and conservation score

- removed the output folder: -> now all plots can be evaluated on the
  graphical device -> clustalOmega output file is saved only if a path
  is provided either in 'alignSequences' method or 'setup' method

- in 'lmPlot' if the mode is 'gene' and only one sequence is analyzed,
  domains are represented within the plot and no logo plot and trident
  score are visualized

- the method 'paths' has been removed and all the options are now
  accessible through 'lmParams' method


MeSHDbi
-------

Changes in version 1.3.2:

- CITATION is added

Changes in version 1.3.1:

- dbBeginTransaction is changed to dbBegin caused by Hadley's
  specification of RSQLite

metagenomeSeq
-------------

Changes in version 1.9:

- Added flexibility in formula choice for fitTimeSeries

- Added readability in ssPermAnalysis

- Fixed default in plotClassTimeSeries (include = c("1",...))

- Added fitTimeSeries vignette

- Removed interactiveDisplay to namespace - moved to suggests

- Fixed ordering of MRtable,MRfulltable first four columns

- modified df estimated through responsibilities

- renamed fitMeta to fitLogNormal - a more appropriate name

metaMS
------

Changes in version 1.3.5:

- Include xset object in GCresults data file

Changes in version 1.3.4:

- Optimize alignmentLC function

Changes in version 1.3.3:

- various small changes in code (esthetic only)

Changes in version 1.3.2:

- Changed LC settings to explicitly include retcor arguments

- Corrected LC bug on the PeakTable export when data files are stored
  into multiple subfolders

Changes in version 1.3.1:

- Changed LC settings so that the retcor algorithm is passed correctly
  to the underlying XCMS function

- Changed LC pipeline so that all user-defined settings are actually
  used

metaseqR
--------

Changes in version 1.5.31 (2015-03-05):

NEW FEATURES

- Added support for both forward and reverse RNA-Seq library
  preparation protocol. Thanks to Ben Elsworth, Garvan Institute of
  Medical Research, Syndey, Australia.

- Added Reads per Gene Model length (rpgm) measurement, where gene
  model length is either the sum of exon lengths or the gene length
  (depending on count.type argument).

BUG FIXES

- Fixed a problem caused by a bug in baySeq 2.0.50 code which caused an
  analysis based in baySeq to crash.

Changes in version 1.5.21 (2015-01-09):

NEW FEATURES

- None

BUG FIXES

- Fixed small bug regarding annotation retrieval when directly using
  the get.annotation function for certain genomes (e.g. pig). Thanks to
  Rathi Ryan from Synthetic Genomics, Inc.

Changes in version 1.5.15 (2014-12-22):

NEW FEATURES

- Minor code improvements

BUG FIXES

- Introduced a (hopefully robust) fix for constant changes in Ensembl
  schema for human GRC37 (hg19). It increases annotation retrieval time
  though

- Fixed problem when input counts are a data.frame and count.type in
  exon, introduced during package upgrades. Thanks to Ben Elsworth,
  Garvan Institute of Medical Research, Austria.

- Replaced "==" with all.equal while checking if the sum of weights is
  one when meta.p="weight". Of course this is the correct way to
  compare. Thanks to Ben Elsworth, Garvan Institute of Medical
  Research, Syndey, Australia.

- Fixed problem related to a bug in the current version of edgeR GLM
  and the usage of ... in glmFit.

Changes in version 1.5.2 (2014-12-30):

NEW FEATURES

- Added the argument exclude.list, to provide the possibility of
  excluding samples from any analysis, e.g. in the case where after a
  preliminary metaseqr run, some samples are not of adequate quality.
  It is useful when a previous analysis is restored, so as not to
  repeat the construction of a gene model from exons for example.

BUG FIXES

- None

Changes in version 1.5.1 (2014-10-31):

NEW FEATURES

- Support for pig genome (Sus scrofa, susScr3 in UCSC)

BUG FIXES

- None

MethylAid
---------

Changes in version 1.1.10:

- added functionality for MethylAidData to show as a reference data set

Changes in version 1.1.9:

- coloring of individual samples add if number of is not >25

Changes in version 1.1.5:

- restructered ui.R

Changes in version 1.1.4:

BUG FIXES

- fixed outlier detection using all quality control plots using an
  initialization method

Changes in version 1.1.3:

- add option/functionality to show background data

Changes in version 1.1.2:

- moved filter thresholds from shiny application to R as argument to
  visualize

minfi
-----

Changes in version 1.13:

- read.450k.exp has support for argument base when targets is supplied.
  Thanks to Brent Pedersen <bpederse@gmail.com> for noticing this and
  providing an initial fix.

- changed the default behaviour of read.450k.exp.  If called using a
  targets argument created by read.450k.sheet, you should not also give
  it a base argument (which was always superfluous).

- Some NAMESPACE imports fixes.

- getGenomicRatioSetFromGEO added to read directly from GEO and create
  a GenomicsRatioSet. Thanks Tim Triche for writing the original
  function.

- makeGenomicRatioSetFromMatrix added. This function turns a matrix
  into a GenomicRatioSet. The 450K feature IDs need to be supplied or
  in the rownames of the matrix.

- makeGenomicRatioSetFromMatrix added to convert matrices to
  GenomicRatioSets. This can be useful for reading in files with beta
  values and turning into object that can be directly passed to
  bumphunter and blockFinder.

- readGEORawFile added to read raw intensity files provided as
  Supplementary Material on GEO. The files include the unmethylated and
  methylated signals. The new function returns a GenomicMethylSet which
  permits you to seamlessly apply minfi preprocessing functions.

- readTCGA is wrapper for makeGenomicRatioSetFromMatrix that reads in
  files in the TCGA format. The function is very specific to this
  format.

- Minor coding fixes including some NAMESPACE issues, missing pData<-
  methods, replace require() with requireNamespace().

- cpgCollapse now works for GenomicRatioSets since it no longer
  attempts to summarize CN data when passed a GenomicRatioSet.

- estimateCellCounts now works on only 2 cell types.

- Various NAMESPACE fixes.

- the gaphunter function by Shan Andrews has been added. We welcome
  Shan as a contributing author.

missMethyl
----------

Changes in version 1.1.10:

- This version of the package contains functions to perform SWAN
  normalisation, RUV normalisation and differential methylation
  analysis, differential variability analysis and gene ontology testing
  for the 450K array.

Changes in version 1.1.3:

- Added the gometh function to perform gene ontology analysis.

Changes in version 1.1.2:

- Added functions to perform RUV normalization and differential
  methylation analysis.

mogsa
-----

Changes in version 0.99.4:

NEW FEATURES

- update prepGraphite function according to the update of Graphite
  pacakge

monocle
-------

Changes in version 1.1.5:

- After grappling with Bioconductor build issues in HSMMSingleCell
  related to VGAM updates, we have moved to a different data layout in
  HSMMSingleCell, which has caused some changes in the vignette.

Changes in version 1.1.1:

- Fixed a bug in responseMatrix() that occurs when you don't have any
  genes that fail VGAM fitting

motifStack
----------

Changes in version 1.11.8:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- Fix the bugs if background values don't have a name.

Changes in version 1.11.7:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- Fix the bugs in readPWM

Changes in version 1.11.6:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- Fix the bugs in motifCircos when input r.pfm and r.pfm2

Changes in version 1.11.5:

NEW FEATURES

- export $ and $<- methods for classes

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.11.4:

NEW FEATURES

- add function getRankedUniqueMotifs and pfm2pwm

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.11.3:

NEW FEATURES

- add function motifPiles

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.11.2:

NEW FEATURES

- export function highlightCol

- add function motifCircos

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.11.1:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- Change the format of help file to avoid lines wider than 90
  characters NOTE

msa
---

Changes in version 1.0.0:

- first official release as part of Bioconductor 3.1

MSGFgui
-------

Changes in version 1.1.2:

- Update to Bootstrap 3 to work with shiny 0.11

MSnbase
-------

Changes in version 1.15.18:

- fix failing test_MSnExp::readMSData unit test on Windows i386
  (@.cache$size being different on that arch) [2015-04-14 Tue]

- merge @vladpetyuk PR #50 fix of combine features bug [2015-04-14 Tue]

Changes in version 1.15.17:

- add TMT10 paragraph and fig to demo vignette [2015-04-09 Thu]

Changes in version 1.15.16:

- support for TMT10 plex [2015-04-08 Wed]

Changes in version 1.15.15:

- using serial parallelisation during quantitation in unit tests
  [2015-04-02 Thu]

Changes in version 1.15.14:

- new msnset data, used in various examples instead of quantifying the
  itraqdata experiment over and over again [2015-04-01 Wed]

Changes in version 1.15.13:

- improve nbavg imputation description and add example [2015-03-22 Sun]

- reduce compareSpectra example timing [2015-03-30 Mon]

Changes in version 1.15.12:

- average neighbour imputation for ordered fractions along a gradient
  [2015-03-21 Sat]

Changes in version 1.15.11:

- update malformed Description [2015-03-20 Fri]

Changes in version 1.15.10:

- update compareSpectra man [2015-03-19 Thu]

- ExpressionSet <-> MSnSet coercion [2015-03-19 Thu]

Changes in version 1.15.9:

- use S4Vectors' isEmpty generic [2015-03-10 Tue]

Changes in version 1.15.8:

- Improved imputation section in demo vignette [2015-03-05 Thu]

Changes in version 1.15.7:

- Importing ProtGenerics [2015-02-28 Sat]

- selective import of MALDIquant [2015-03-02 Mon]

Changes in version 1.15.6:

- add intensity column to the calculateFragments output; closes #47;1
  [2015-02-02 Mon]

- add method argument to calculateFragments to allow the user choosing
  the highest/closest peak or all peaks in the tolerance range; closes
  #47;2 [2015-02-09 Mon]

- add neutralLoss argument to calculateFragments and calculate loss of
  water and ammonia; closes #47:3 [2015-02-19 Thu]

- new imputation methods via imputeLCMD and norm [2015-02-09 Mon]

- vignette updates [2015-02-09 Mon]

- imputation unit test [2015-02-24 Tue]

Changes in version 1.15.5:

- update dependency to mzID >= 1.5.2 [2015-01-28 Wed]

- rewrite addIdentificationData and change its signature [2015-01-28
  Wed]

- add methods addIdentificationData which work on MSnExp and MSnSets
  using filenames (characters), mzID objects and data.frames; see #42;
  closes #45; [2015-01-28 Wed]

- add a section about MSmaps in the vignette [2015-02-02 Mon]

- MSmap has a new zeroIsNA argument to set all 0 values of the map to
  NA. This simplifies the resulting plot3D figure. [2015-02-02 Mon]

Changes in version 1.15.4:

- partly rewrite readMgfData [2015-01-19 Mon]

- Typo in addIdentification data man [2015-01-20 Tue]

- use Biocstyle [2015-01-23 Fri]

- replace require with requireNamespace [2015-01-23 Fri]

Changes in version 1.15.3:

- plot,Spectrum2,character to add fragment ions based on peptide
  sequence [2014-11-20 Thu]

- update vignette with above [2014-11-20 Thu]

- comment parallel code in quantify man [2015-01-09 Fri]

Changes in version 1.15.2:

- merged sgibb's pull request exporting overloaded methods as well
  [2014-11-02 Sun]

- updated MSnSet validity, checking that exprs(.) is a matrix
  [2014-11-12 Wed]

- fix error (missing header extraction) is MSmap example [2014-11-18
  Tue]

Changes in version 1.15.1:

- Fixing error when id file has no spectrumFile info (see issue #39)
  and return a warning (instead of an error) when the file used to
  create the MSnExp/MSnSet and mzid file were different [2014-10-15
  Wed]

Changes in version 1.15.0:

- Devel version for Bioc 3.1

MSnID
-----

Changes in version 1.1.6:

- Added explicit mc.cores argument for "Grid" type of optimization

- Fix of trimming to basename on POSIX

Changes in version 1.1.5:

- Always trimming spectrumFile and databaseFile to basename after
  reading data from mzIdentML files.

Changes in version 1.1.4:

- Catching up with changes in mzID package version 1.5.3

Changes in version 1.1.3:

- Updated data, documentation and unit test to reflect the latest
  changes in mzID

Changes in version 1.1.2:

- id_quality and evaluate_filter may have multiple values for the level
  argument. By default it is all three: "PSM", "peptide" and
  "accession". The return value for both methods is switched to
  matrices.

- added Laurent Gatto as a contributor

Changes in version 1.1.1:

- fix mishandling PSM

- use mzR::psms generic

Changes in version 1.1.0:

- Bumping versions after creating 3.0 release branch.


mzID
----

Changes in version 1.5.3:

- Use "location" attribute instead of "name" to create the
  spectrumFile/databaseFile column in the flattened data.frame.

mzR
---

Changes in version 2.1.12:

- Remove superfluous BiocGenerics in Suggests [2015-03-01 Sun]

- Merged KK's Makefile to generate pwiz lib for windows [2015-03-01
  Sun]

Changes in version 2.1.11:

- Using generics from BiocGenerics (fileName, score) and ProtGenerics
  [2015-02-28 Sat]

Changes in version 2.1.10:

- Add instrumentInfo() and runInfo() for CDF backend (closes issue #22)

- Add precompiled libpwiz.a to reduce compile time on Windows (closes
  issue #21), thanks to KK

Changes in version 2.1.9:

- Fix a compiler warning on OSX in the pep_XML serialiser

Changes in version 2.1.8:

- Fix build on unix with netcdf.h in non-standard location

Changes in version 2.1.7:

- Fix score segfault (see issue #18) and remove usage of ListBuilder.

- Update unit test to reflext new score code.

Changes in version 2.1.6:

- adding an acquisitionNum column to psms (closes issue #17)
  [2015-02-05 Thu]

Changes in version 2.1.5:

- don't print '1' to the console when calling get3Dmap [2015-01-29 Thu]

- add sequence length in psms output (closes issue #19) [2015-02-04
  Wed]

Changes in version 2.1.4:

- documentation and vignette: mzIdentML version 1.1 support only
  [2015-01-23 Fri]

Changes in version 2.1.1:

- remove dependency on faahKO


nethet
------

0.99.5: First Bioconductor-devel release.

NetPathMiner
------------

Changes in version 1.3.1:

- Fixed a bug SBML2igraph (metabolic).

npGSEA
------

Changes in version 1.3.8:

- Adjusted email of Jessica L. Larson

Changes in version 1.3.7:

- Scaling of X to have unit variance is now optional (and is the
  current default)

Changes in version 1.3.6:

- Scaling of X and Y is now optional (and is the current default)

Changes in version 1.3.5:

- Set a minimum p-value for the Beta approximation as the default

Changes in version 1.3.4:

- Added some checks to the values of XG and Y

Changes in version 1.3.3:

- Made the epsilonBetaAdj=FALSE as the default option.

Changes in version 1.3.2:

- Set a minimum p-value for the Beta approximation as the default, but
  allowed it be changed with the 'epsilonBetaAdj' parameter in the main
  function.

Changes in version 1.3.1:

- Added a requirement that each treatment group have at least two
  observations.

- Set a minimum p-value for the Beta approximation so that it never
  returns a value of zero

oligo
-----

Changes in version 1.32:

USER VISIBLE CHANGES

- Fixed vignette

- Added initial support to Generic Arrays

- Using BiocGenerics definition of normalize

omicade4
--------

Changes in version 1.7.1:

NEW FEATURES

- mcia will return pairwise RV coefficient in the "$mcoa"

- updated reference

PAA
---

Changes in version 1.1.1 (2015-03-13):

BUG FIXES

- Bug in the function batchAdjust() removed. Because the argument "mod"
  was not passed as factor to sva's function ComBat() errors were
  thrown and execution stopped for the latest versions of sva. This has
  been corrected.

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

Pbase
-----

Changes in version 0.6.13:

- Updating Biomart query to reflect changes [2015-04-06 Mon]

Changes in version 0.6.12:

- Fix filling in tracks plots [2015-03-07 Sat]

Changes in version 0.6.11:

- Documentation for the plotting heler functions [2015-02-28 Sat]

Changes in version 0.6.10:

- plotting helper functions [2015-02-28 Sat]

- fixed peptide grouping [2015-02-28 Sat]

Changes in version 0.6.9:

- Fixing junction peptides splitting to support peptides reaching > 2
  exons [2015-02-27 Fri]

- supporting more protein coding biotypes [2015-02-27 Fri]

- rewrite and simplify proteinCoverage; closes #19 [2015-02-27 Fri]

Changes in version 0.6.8:

- fix splitting exon junction ranges [2015-02-26 Thu]

Changes in version 0.6.7:

- The genomic coordinates of peptides that overlap exon junctions are
  now correctly split instead of spanning over the exon and are grouped
  accordingly. [2015-02-26 Thu]

Changes in version 0.6.6:

- update mapping vignette and use alignment to identify correct
  transcript [2015-02-25 Wed]

Changes in version 0.6.5:

- see also mapping vignette in ?mapToGenome [2015-02-20 Fri]

- calculateHeavyLabels, isCleaved, proteotypic and proteinCoverage are
  is now a function [2015-02-20 Fri]

- added accessors to Pbase:::pplot() and update Pbase-data vignette
  accordingly [2015-02-21 Sat]

Changes in version 0.6.4:

- typos in man [2015-02-19 Thu]

- Added use.names argument ot etrid2gel [2015-02-19 Thu]

- mapToGenome,Proteins,GRangesList (remove GRanges) using names for
  1-to-many mapping and document in mapping vignette [2015-02-19 Thu]

Changes in version 0.6.3:

- new fig in Pbase-data vignette [2015-02-14 Sat]

- new section for pmapToGenome in mapping vignette [2015-02-14 Sat]

Changes in version 0.6.2:

- Add uniprot release to p metadata [2015-02-14 Sat]

- added a metadata update method [2015-02-14 Sat]

Changes in version 0.6.1:

- [p]mapToGenome methods [2015-02-14 Sat]

Changes in version 0.6.0:

- fixed typos in vigs [2015-02-12 Thu]

- new etrid2grl function to convert Ensembl transcript identifiers to
  GRangesList objects [2015-02-13 Fri]

- initial mapToGenome function [2015-02-13 Fri]

- regenerated p data using latest UniProt release and added
  acols(p)$ENST [2015-02-13 Fri]

pdInfoBuilder
-------------

Changes in version 1.32:

BUG FIXES

- Fixed dependency changes for RSQLite 1.0.0

USER VISIBLE CHANGES

- GenericArray is now available to accommodate oligo support for custom
  CDFs

PhenStat
--------

Changes in version 2.1.3:

NEW FUNCTIONALITY

- Box-Cox transformation is applied when needed for the continues data.
  Transformation values are shown in output functions. Genotype effects
  are converted back into original scale if transformation was applied.

- New statistical method is implemented for the categorical data -
  Logistic Regression.

- New function is created boxplotSexGenotypeBatchAdjusted() that allows
  to evaluate Batch effect.

- Major architectural changes: classes PhenList and PhenTestResult are
  S4 classes now.  Different new methods are provided for the classes.

- PhenTestResult object contains dataset analysedDataset() that was
  used for the analysis including original values, transformed values
  (if transformation has been applied) and values adjusted for batch
  effect.

COMPATIBILITY ISSUES

- The vectorOutput() has new element for transformation values which
  have increased its length.



plethy
------

Changes in version 1.5.10:

- Modified 'show' method of 'BuxcoDB' objects to utilize available
  metadata

Changes in version 1.5.7:

NEW FEATURES

- Added ability to generate SQLite packages similar to AnnotationDbi

polyester
---------

1.2.3: change URL of example chromosome 22 sequence dataset / do not
        test internet-connectivity-dependent stuff during check

1.2.2: bug fix in fragment start position simulation with positional
        bias

1.2.1: major refactor: - reduced number of arguments for wrapper
        functions (simulate_experiment and
        simulate_experiment_countmat) - shortened source code in
        wrapper functions by adding several internal helper functions

1.2.1: several new features added: - GC bias in expression - empirical
        fragment length distribution model available - positional bias
        in fragmentation step - empirical error models available -
        custom library size factors available -
        simulate_experiment_empirical function added as quick way to
        simulate an experiment with abundances derived from a data set

prebs
-----

Changes in version 1.7.1:

- Introduced a new parameter sum.method. It is used to choose the
  microarray summarization mode: RPA or RMA. The default value is RPA.

ProCoNA
-------

Changes in version 1.4.1:

- importing libraries in the namespace

pRoloc
------

Changes in version 1.7.13:

- use donttest instead of dontrun [2015-04-09 Thu]

Changes in version 1.7.12:

- don't run knntl example to reduce checking time and timeout on
  Windows [2015-04-08 Wed]

Changes in version 1.7.11:

- fix splitTh, closes issue #49 [2015-04-06 Mon]

Changes in version 1.7.10:

- Change in vignette to work on zin1: using 12 random th rows. See
  issue #49 for details. Fixed in version 1.7.11. [2015-04-03 Fri]

Changes in version 1.7.9:

- updated tl vignette [2015-04-02 Thu]

- depending on latest (1.5.8) pRolocdata [2015-04-02 Thu]

Changes in version 1.7.8:

- updating tl vig [2015-04-02 Thu]

- update getParams documentation [2015-04-02 Thu]

Changes in version 1.7.7:

- renaming theta.scores to knntl.scores [2015-03-24 Tue]

Changes in version 1.7.6:

- added the theta inductive transfer infrastructure [2015-02-05 Thu]

- theta vignette stub [2015-02-06 Fri]

- rename getClasses to getMarkerClasses [2015-02-06 Fri]

- added the infrastructure to create GO MSnSet [2015-02-07 Sat]

- Fixed ml vignette [2015-02-16 Mon]

- new filterZeroRows function [2015-03-10 Tue]

- hpa data section [2015-03-10 Tue]

- theta sections [2015-03-10 Tue]

- deprecate getRegulari[z|s]edParams [2015-03-11 Wed]

Changes in version 1.7.5:

- Fix vignettes: run bibtex and pdflatex twice and typo [2015-02-03
  Tue]

Changes in version 1.7.4:

- Use default Sweave call to build vignette [2015-01-24 Sat]

Changes in version 1.7.3:

- use Biocstyle [2015-01-23 Fri]

- replace library/require by requireNamespace [2015-01-23 Fri]

Changes in version 1.7.2:

- added t-SNE method to plot2D [2015-01-14 Wed]

- Updated NAMESPACE imports [2015-01-14 Wed]

Changes in version 1.7.1:

- updated vignettes with markers.orig [2014-10-30 Thu]

- updated ml tests [2014-10-30 Thu]

Changes in version 1.7.0:

- new devel version, Bioc 3.1

pRolocGUI
---------

Changes in version 1.1.5:

- updating failing unit test [2015-04-03 Fri]

Changes in version 1.1.4:

- fix R_HOME error [2015-02-26 Thu]

Changes in version 1.1.3:

- don't require GNU make [2015-02-11 Wed]

Changes in version 1.1.2:

- Fix bug with FoIs and multiple data sets (reported by Harriet
  Parsons) [2015-02-06 Fri]

Changes in version 1.1.1:

- Update README with uptodate installation instructions [2014-10-14
  Tue]

- handling and filtering missing value in input data [2014-10-21 Tue]

proteoQC
--------

Changes in version 1.3.2:

- Update description file

Changes in version 1.3.1:

- Update the function labelRatio

- Update the vignette

ProtGenerics
------------

Changes in version 0.99.3:

- update man with when and why ProtGenerics [2015-03-16 Mon]

Changes in version 0.99.2:

- updating DESCRIPTION as with BiocContributions::clean [2015-03-03
  Tue]

Changes in version 0.99.1:

- Adding NEWS file [2015-02-27 Fri]

Changes in version 0.99.0:

- adding ProtGenerics to Manifest

qcmetrics
---------

Changes in version 1.5.1:

- Using Biocstyle [2015-01-23 Fri]

QDNAseq
-------

Changes in version 1.2.4 (2015-01-21):

OTHER

- update package maintainer

Changes in version 1.2.3 (2015-01-20):

BUG FIXES

- createBins() now also works for species not supported by GenomeInfoDb
  (such as Canis lupus familiaris)

Changes in version 1.2.2 (2014-12-23):

IMPROVEMENTS

- applyFilters() now ignores the residual filter for chromosomes for
  which it does not exist (instead of always filtering out the entire
  chromosome)

Changes in version 1.2.1 (2014-11-01):

IMPROVEMENTS

- add an asterisk in the noise measure to clarify it's not a regular
  standard deviation or variance, but first scaled with the mean (so
  that the mean is 1.0 and the relationship holds between variance and
  1/N, where N is the average number of reads per bins)

- clarify the use of log2/sqrt transformation in segmentBins()

qpgraph
-------

Changes in version 2.20:

USER VISIBLE CHANGES

- Updated the vignette "Estimate eQTL networks using qpgraph". It
  includes more detailed simulations illustrating the steps involved in
  the estimation of eQTL networks with qpgraph.

BUG FIXES

- Bugfix on the display of eQTL networks with hive plots

QuartPAC
--------

Changes in version 0.99.0:

- First release of the QuartPAC package.

- Currently supports the iPAC, GraphPAC and SpacePAC packages when
  extending clustering results to quarternary data.

QuasR
-----

Changes in version 1.8.0:

PUBLICATION

- A paper describing the QuasR package has been published (see
  citation("QuasR")): Gaidatzis D. et al., Bioinformatics, 2014. doi:
  10.1093/bioinformatics/btu781

NEW FEATURES

- started implementing support for BiocParallel, allowing QuasR to also
  run on batch clusters (support is currently limited to qQCReport)

qvalue
------

Changes in version 1.99:

This update of the qvalue package includes

- removal of tcltk usage

- the use of ggplot2 in the plotting functions

- an updated vignette

- vectorizing sections of the code

- function name changes

- inclusion of local FDR estimation

R3CPET
------

2015.02.1:

r3Cseq
------

Changes in version 1.13.1 (2015-01-26):

- fixed the bugs in function "export3Cseq2bedGraph"

- fixed the bugs in function "getViewpoint"

- fixed the bugs in function "export3CseqRawReads2bedGraph"

- added "mm10" and "rn5" genomes into the package

Rbowtie
-------

Changes in version 1.7.9:

NEW FEATURES

- updated bowtie to version 1.1.1 (currently not using long index
  feature, renaming bowtie-align-s and bowtie-build-s to bowtie and
  bowtie-build)


ReactomePA
----------

Changes in version 1.11.9:

- use mapIds instead of mget in TERM2NAME. mapIds will not throw error
  when pathID is old/deprecated. <2015-04-10, Fri>

Changes in version 1.11.8:

- according to https://support.bioconductor.org/p/63024/#66438, modify
  gsePathway to only test pathways. <2015-04-09, Thu>

Changes in version 1.11.7:

- update use_internal_data parameter to ... according to the change of
  DOSE <2015-04-01, Wed>

Changes in version 1.11.6:

- update according to DOSE <2015-03-31, Tue> add use_internal_data
  parameter in internal functional of enrichPathway.

Changes in version 1.11.5:

- update according to DOSE <2015-03-01, Sun> add use_internal_data
  parameter in getGeneSet.

Changes in version 1.11.2:

- remove import org.Hs.eg.db <2015-01-22, Thu>

Changes in version 1.11.1:

- fix some pathways don't have a name <2014-11-17, Mon> ## >
  get("5493857", reactomePATHID2EXTID) ## [1] "510850" "523328"
  "282187" "282188" ...  ## > get("5493857", reactomePATHID2NAME) ##
  Error in .checkKeys(value, Lkeys(x), x@ifnotfound) : ## value for
  "5493857" not found --> This is not a bug in ReactomePA, but an issue
  of reactome.db, refer to: -->
  https://support.bioconductor.org/p/63024/

RedeR
-----

Changes in version 1.15.0:

- Improved app reliability.

- Implemented new compatibility checks for call-back functions.

- Fix a client-server issue related with some MS-Windows
  specifications.

regionReport
------------

Changes in version 1.1.9:

NEW FEATURES

- Introduced renderReport() which creates a simple exploratory report
  for any set of genomic regions. It allows the user to further
  customize the report by using a child file.

- You can now use the 'output_format' advanced parameter on both
  renderReport() and derfinderReport() to output a PDF file instead of
  an HTML file. The interactive tables are lost and only the top 20
  rows are shown.

Changes in version 1.1.8:

SIGNIFICANT USER-VISIBLE CHANGES

- Adapted to work with bumphunter >= 1.7.6

Changes in version 1.1.7:

NEW FEATURES

- Users can now control 'output_format' and 'clean' options from
  rmarkdown::render() when running derfinderReport()

Changes in version 1.1.3:

BUG FIXES

- Adapted derfinderReport() to derfinder 1.1.5

ReportingTools
--------------

Changes in version 2015-3-27:

- Updated email for Jessica L. Larson


rGREAT
------

Changes in version 0.99.5:

- because analysis depends on internet connection, jump to next version
  to enforce rebuilding of the vignette on bioconductor.

Changes in version 0.99.4:

- retrieve categories and corresponding ontologies from GREAT

rgsepd
------

Changes in version 0.99.15:

BUG FIXES

- bugfix in gene ID systems passed to goseq.

Changes in version 0.99.12:

NEW FEATURES

- expanding user manuals to clarify how figures and clustering is
  performed.

Changes in version 0.99.11:

BUG FIXES

- bugfix in ProcessAll cardinality message.

- Error catching around PCA.

- Heatmap genelist now respects Annote_Filter file when LIMIT$HARD.

- Cleaned up warnings when no rows for HMA.

- Bugfix in PCA regarding non-negative definite covariance from
  princomp's eigen: replaced all instances with prcomp().

Changes in version 0.99.10:

BUG FIXES

- bugfix in G$LIMIT$baseMean and explanations of naming conventions in
  the vignette.

Changes in version 0.99.4:

SIGNIFICANT USER-VISIBLE CHANGES

- replacing DESeq with DESeq2 Fall2014

rhdf5
-----

Changes in version 2.12.0:

NEW FEATURES

- Filenames are expanded with normalizePaths.

- New function h5set_extent implemented.

- New low level function H5Sset_extent_simple implemented.

BUG FIXES

- Segmentation fault while writing data type names for uncommitted data
  types.

Risa
----

1.9.1: 1. Added citation information. 2. The identification of
        investigation files now avoids editor backup files (ignoring
        files like "i_Investigation~") 3. Fixed issue where the
        investigation file wasn't fully read, due to the number of
        columns in higher rows being greater than the first five rows
        (as used by read.table) 4. Moved inst/doc to vignettes folder
        as required since R 3.1.0 5. Imported packages in Depends 6.
        Removed 'library' or 'require' calls to packages already
        attached by Depends. 7. Fixed no visible binding for global
        variable.

RNAprobR
--------

Changes in version 0.99.4:

NEW FEATURES

- BED2txDb is now using rtracklayer fucntion import() to load BED
  files.

RnaSeqSampleSize
----------------

Changes in version 0.99.8 (2015-04-12):

BUG FIXES

- Fix a bug in example of optimize_parameter function, which was caused
  by update of heatmap3 package. Thanks Dan.

Changes in version 0.99.7 (2015-04-06):

BUG FIXES

- Fix a bug in example of convertId function, which was caused by
  update of BioMart database. Thanks Dan.

Changes in version 0.99.5 (2015-04-06):

SIGNIFICANT USER-VISIBLE CHANGES

- Minor changes for genes with less than 1 read count;

Changes in version 0.99.4:

NEW FEATURES

- The user friendly web interface was improved

BUG FIXES

- Bug fixed: Max sample size for distribution based sample size
  estimation;

- Bug fixed: Keep consistent for all function names;

SIGNIFICANT USER-VISIBLE CHANGES

- Other improvement based on the comments from Bioconductor reviewer;

Changes in version 0.99.3 (2014-11-23):

SIGNIFICANT USER-VISIBLE CHANGES

- New parameter: countFilterInRawDistribution,
  selectedGeneFilterByCount;

- Other improvement based on the comments from Bioconductor reviewer;

Changes in version 0.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette was switched to a BiocStyle;

- Other improvement based on the comments from Bioconductor reviewer;

Changes in version 0.99.1 (2014-10-19):

SIGNIFICANT USER-VISIBLE CHANGES

- Parameter was used to determine the power of genes below minAveCount;

- Some recommendations from BiocCheck were improved;

Changes in version 0.99.0 (2014-10-16):

SIGNIFICANT USER-VISIBLE CHANGES

- Submit to Bioconductor.


RnBeads
-------

Changes in version 0.99.22:

- Incorporated bioconductor requirements for parallel processing.
  "parallel.disable" is now "parallel.teardown"

Changes in version 0.99.20:

- RnBeads no longer suggests IlluminaHumanMethylation450k.db.

Changes in version 0.99.19:

- Fixed a bug in the computation of surrogate variables (SVA) that did
  not take into account other adjustment variables in the model
  (compared to the null model)

- Updated locus profiles to work with the latest version of Gviz

- Differential methylation reports now contain information on group
  sizes

- Included plots of number of sites vs coverage percentiles for each
  sample in the dataset (rnb.plot.num.sites.covg)

- Bisulfite-seq QC reports now contain an additional section with
  summary statistics and plots when the coverage histogram plots are
  enabled

Changes in version 0.99.18:

- Implementation of class constructors (the S4 constructors are not
  available anymore)

- Analyses can now be started from an RnBSet which is stored on the
  hard drive (using save.rnb.set) via the rnb.set.dir data type

- New IDAT loading routine, based on the illuminaio package

- Multiple fixes in wrappers for wateRmelon normalization methods

- Fixed the RefFreeEWAS wrapper

- Updated the Houseman, 2012 reference-based method implementation

- RnBSets now know the version with which they were created

- Various other bugfixes and improvements to reports and documentation

Changes in version 0.99.17:

- Consistency issues with Bioconductor 3.0 fixed

- Fixed a bug relating to no coverage masking being conducted for
  sequencing data

- Various other bugfixes

Changes in version 0.99.16:

- Enhanced XML-based analysis though the tag preanalysis.script

- Optimized loading of the sequencing data sets

- Enhanced region profiles plots. Added region site distribution plots

- Locus profiles can now be generated in the exploratory report for
  regions specified in custom bed files and for genes listed as gene
  symbols

- Improved cell type heterogeneity inference (no cell types excluded
  anymore)

- Added more support for paired differential methylation with limma

- Some functions were renamed: merge.samples->mergeSamples,
  add.pheno->addPheno, add.region.subsegments->addRegionSubsegments,
  and all functions operating on Report

- Added helper classes for submitting RnBeads to a scientific compute
  cluster. Current implementation includes Sun Grid Engine

- The default pipeline now uses ff functionality and saves intermediate
  objects into the report directory

- Minor bugfixes

Changes in version 0.99.15:

- Restructuring of the pipeline modules: ** the loading module has been
  renamed to import ** the prefiltering,normalization and postfiltering
  modules have been summarized in a new module: preprocessing ** the
  batch and profiles modules have been summarized in a new module for
  exploratory analysis ** the export module has been renamed to tracks
  and tables ** corresponding option names have been changed ** see the
  overview figure on the website or the vignette for a quick overview
  on the new module structure

- Multiple RnBSets can now be concatenated with the add() function

- Multiple samples in an RnBSet can be merged using the merge.samples()
  function

- Gender prediction can be performed on Infinium 450k datasets

- Minor updates on covariate adjustment

- Updates on tissue heterogeneity estimation

- Calling differentially methylated sites with RefFreeEWAS now supports
  paired design

- Multiple minor bugfixes and performance improvements

- Vignette and documentation updates

Changes in version 0.99.13:

- Added new module on annotation inference

- Added correction for cell type heterogeneity

- Added SVA functionality

- Minor bugfixes

- Updates to the vignette and other documentation

Changes in version 0.99.12:

- Support Bismark coverage file loading

- Enhanced documentation and logging of loading steps

- Minor bugfixes

- Updates to the vignette and other documentation RnBeads 0.99.11

- Accommodate data packages for individual genomes

- Performance: disk space usage when using disk dumping, options for
  subsetting sites when computing distance matrices for clustering, PCA
  and MDS

- Region subsegmentation

- The default method for differential methylation p-values is limma

Changes in version 0.99.10:

- Performance improvements

- Option to keep big matrices on the hard drive rather than main memory

- Restructuring of the filtering modules. Parts of the filtering steps
  are executed before normalization others afterwards

- Improvements in normalization: more methods supported

- Cosmetic changes to some of the plots

Changes in version 0.99.9:

- Bugfixes in paired analysis

- Usability of bisulfite sequencing

Changes in version 0.99.8:

- New normalization methods integrated

- Improved arguments to rnb.run.analysis

- Improved parsing of the sample annotation table

- Paired analysis (testing stage)

- Webservice installed

- Multiple bugfixes

Changes in version 0.99.7:

- Support for background subtraction and BMIQ normalization of Infinium
  450k data

- Support for differential methylation analysis on all pairs of sample
  groups

Changes in version 0.99.6:

- Locus Profiles

- Support for parallel computing

Changes in version 0.99.0:

- Many additional features including bisulfite sequencing mode, the
  mouse genome, data export, ...



rpx
---

Changes in version 1.3.1:

- update vignette to reflect new files in px example, and avoid
  downloaded the raw data [2015-03-31 Tue]

Rqc
---

Changes in version 1.2:

NEW FEATURES

- Study design metadata system (group ID, strand, etc) added

- Statistics for trimming step added

- Per sample mean quality box plot added

- Per cycle mean quality Principal Component Analysis (PCA) plot added

- Overrepresented sequencing reads plot added

USER VISIBLE CHANGES

- RqcResultSet class stores data as frequency tables

- BPPARAM argument added to rqc and rqcQA functions

BUG FIXES

- Memory usage improved

- Input file validation system

- Rqc does not depend on X11 to generate plots

- Rqc's reporter can generate plots for more than 10 files

- Fixed bug related with reshape2 dependency

- Fixed bug related with table function

Rsamtools
---------

Changes in version 1.19:

SIGNIFICANT USER-VISIBLE CHANGES

- FaFile accepts a distinct index file

- Support for cigars > 32767 characters

- Mate pairs use pos and mpos values calculated modulo target length
  for pairing, facilitating some representations of mates on circular
  chromosomes.

- scanBam no longer translates mapq '255' to 'NA'

BUG FIXES

- segfault on file iteration, introduced in 1.19.35, fixed in 1.19.44

- scanBam correctly parses '=' and 'X'


RUVcorr
-------

Changes in version 0.99.1:

- Corrected spelling mistakes in manual.

Changes in version 0.99.0:

- Initial version submitted to Bioconductor.

RUVSeq
------

Changes in version 1.1:

- Corrected typo in vignette

- RUV* methods for matrix now support log data


SeqArray
--------

Changes in version 1.7.1-1.7.5:

- bug fix in getting genotypes if position > 2^31

- add an option 'ignore.chr.prefix' to the function 'seqVCF2GDS'

- 'seqVCF2GDS' ignores the INFO or FORMAT variables if they are not
  defined ahead

- a new action 'push+set' in the function 'seqSetFilter'

- bug fix if 'requireNamespace("SeqArray")' is called from other
  packages


SeqVarTools
-----------

Changes in version 1.5.1:

- Use existing isSNV generic from VariantAnnotation instead of
  redefining

- Use BiocStyle for vignette

SGSeq
-----

Changes in version 1.2.0:

- Renamed class TxVariants to SGVariants

- Renamed class TxVariantCounts to SGVariantCounts

- Renamed findTxVariants() to findSGVariants()

- Renamed getTxVariantCounts() to getSGVariantCounts()

- Parallelization is now controlled with a single cores argument

- Argument max_complexity for predictTxFeatures() controls skipping of
  problematic regions

- getBamInfo() is no longer run as part of analyzeFeatures()

- getSGVariantCounts() now supports obtaining counts from BAM files

- Bug fixes and other improvements

ShortRead
---------

Changes in version 1.25:

SIGNIFICANT USER-VISIBLE CHANGES

- srapply is defunct

- readAligned() for BAM files is deprecated; use
  GenomicAlignments::readGAligned instead.

BUG FIXES

- close opened files when parsing old bowtie, soap, and solexa export
  file formats.

- Don't allow R memory to be released prematurely when processing old
  bowtie file formats / creating external pointers.

- writeFastq,FastqFile obeys 'compress' argument; mode must be
  specified by the caller (typically mode="a")

SigCheck
--------

Changes in version 2.0.0:

- Major revamp. Survival analysis, S4 object, etc.

SIMAT
-----

Changes in version 0.99.3:

- exporting single functions instead of pattern [2015-03-12 Thu]

Changes in version 0.99.2:

- The package is updated based on the comments from Bioconductor
  reviewers [2015-03-01 Wed]

Changes in version 0.99.1:

- The package is updated based on the comments from Bioconductor
  reviewers [2015-02-11 Wed]

Changes in version 0.99.0:

- SIMAT package is created [2015-01-27 Tue]

SNPRelate
---------

Changes in version 1.1.0-1.1.11:

- fix a bug in snpgdsVCF2GDS when 'method="biallelic.only"'

- add 'snpgdsVCF2GDS_R' for the R implementation

- fix a bug in 'snpgdsBED2GDS' if 'family=TRUE'

- 'snpgdsGDS2BED' allows the file name of GDS

- improve 'snpgdsSlidingWindow'

- add an option 'ignore.chr.prefix' to the function 'snpgdsVCF2GDS'

- a new function 'snpgdsHWE'

- v1.1.5: add 'Fst estimation' to the vignette

- v1.1.6: bug fix if 'requireNamespace("SNPRelate")' is called from
  other packages

- v1.1.7: snpgdsPCA uses 'DSPEVX' to compute eigenvalues and
  eigenvectors instead of 'DSPEV' if top eigenvalues are required only
  (significant improvement on computing speed)

- v1.1.8: the original Rnw vignette is replaced by a R Markdown
  vignette

- v1.1.9: a new function 'snpgdsPED2GDS'


specL
-----

Changes in version 1.1.17:

USER UNVISIBLE CHANGES

- removed file argument in genSwathIonLib function

Changes in version 1.1.16:

USER UNVISIBLE CHANGES

- added unit test for genSwathIonLib

Changes in version 1.1.15:

USER UNVISIBLE CHANGES

- added circle plots to specLSet plot method

- added breaks argument in genSwathIonLib methode

Changes in version 1.1.14:

USER UNVISIBLE CHANGES

- LinkedTo Rcpp; added C++ STL lower bound function which is reqired
  for determining overlapping q1 and q3 SWATH windows

Changes in version 1.1.13:

USER VISIBLE CHANGES

- fixed man pages

Changes in version 1.1.12:

USER VISIBLE CHANGES

- impoved package vignette

Changes in version 1.1.11:

USER VISIBLE CHANGES

- modified default parameters of genSwathIonLib

- add content to vignette

Changes in version 1.1.10:

USER VISIBLE CHANGES

- added generate.consensus

Changes in version 1.1.9:

USER VISIBLE CHANGES

- new features in specLSet 'summary' plot

USER UNVISIBLE CHANGES

- refactored merge.specLSet; merge by group_id

- added unit test for merge.specLSet

Changes in version 1.1.8:

USER VISIBLE CHANGES

- renamed annotateProteinID to annotate.protein_id

- added graphics on plot.specLSet method

USER UNVISIBLE CHANGES

- refactored merge

Changes in version 1.1.7:

USER VISIBLE CHANGES

- introduce peakplot for bibliospec object

- introduce LCMS map for bibliospec object

- vignette cosmetics

Changes in version 1.1.6:

USER VISIBLE CHANGES

- introduce specL_bibliospec summary method

Changes in version 1.1.5:

USER VISIBLE CHANGES

- specLSet merge function

- work on specLSet summary method

- specLSet merge function

- work on specLSet summary method

Changes in version 1.1.4:

USER VISIBLE CHANGES

- summary method of specLSet class

- summary method of specLSet class

USER UNVISIBLE CHANGES

- unit test for data containing no iRT peptides

- unit test for data containing no iRT peptides

Changes in version 1.1.3:

USER VISIBLE CHANGES

- renamed write.Spectronaut to write.spectronaut

- write.spectronaut writes filename

- added benchmark section in package vignette

- renamed write.Spectronaut to write.spectronaut

- write.spectronaut writes filename

- added benchmark section in package vignette

Changes in version 1.1.2:

USER VISIBLE CHANGES

- uses modSeq in group_id iff existing

- uses modSeq in group_id iff existing

Changes in version 1.1.1:

USER VISIBLE CHANGES

- streamline modsequence, e.g., AAAMASATTM[+16.0]LTTK for compatibility
  with peakView V2.0

- streamline modsequence, e.g., AAAMASATTM[+16.0]LTTK for compatibility
  with peakView V2.0

SRAdb
-----

Changes in version 1.21.11 (2015-03-22):

- Fixed a bug with listSRAfile function.

Changes in version 1.21.8 (2014-12-12):

- Added fastq table to the database

- Re-wrote getFastqInfo function

- Fastq files can be downloaded by getFastqFile or getSRAfile

- Added <dir2> to the fastq ftp addresses for SRR accessions over
  SRR999999

supraHex
--------

Changes in version 1.5.1:

NEW FEATURES

- Fix the batch algorithm allowing for processing an input data
  containing many zero entries

synapter
--------

Changes in version 1.9.5:

- update cross reference to qvalue::plot.qvalue in Synapter man page;
  closes #86 [2015-03-31 Tue]

Changes in version 1.9.4:

- use biocstyle [2015-01-23 Fri]

- use requireNamespace [2015-01-23 Fri]

Changes in version 1.9.3:

- Filter entries in quantiation final peptides and Pep3D data that
  don't match in their intensity valus; see #42 for details.
  [2014-11-26 Wed]

Changes in version 1.9.2:

- Only suggesting tcltk and tcltk2, as gui is now deprecated
  [2014-11-21 Fri]

Changes in version 1.9.1:

- Deprecating synatperGUI [2014-11-10 Mon]

- Directing questions to the Bioc support site [2014-11-10 Mon]

Changes in version 1.9.0:

- new devel version for Bioc 3.1

TargetSearch
------------

Changes in version 1.24.0:

BUG FIXES

- Add clarification note for "Window" parameter in "RIcorrect"
  function.

TCC
---

Changes in version 1.7.13:

- fix bug in '.testBySamseq'.

- fix bug in 'simulateReadCounts'.

- rewrite the R code.

- 'fit0' and 'fit1' arguments were changed to 'full' and 'reduced',
  'DESeq.test' arguments were deleted.

- 'voom' in limma package wa added for DEG identification.

- the 'test.method' for the comparison of paired two-group can be
  choosed from 'edger', 'deseq', 'deseq2', 'bayseq' and 'voom'.

- change the mode of expression of tcc$simulation$trueDEG. '0'
  indicates a non-DEG and '1' indicates a DEG, only '0' and '1' are
  used in this field.

- 'groups' argument in 'plot.TCC' function has changed to 'group'.

TFBSTools
---------

Changes in version 1.5.0:

NEW FEATURES

- searchSeq now works on DNAStringSet.  The names in DNAStringSet will
  be used as seqname in SiteSet.

- More user-friendly output in data.frame or GRanges from SiteSet,
  SiteSetList, et al.

- reverseComplement works on XMatrix now.

- faster implmentation of searchPairBSgenome

TPP
---

Changes in version 1.0.0:

- First release of the TPP package for analyzing Thermal Proteome
  Profiling experiments.

trackViewer
-----------

Changes in version 1.3.4:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- improve efficiency of importScore

Changes in version 1.3.3:

NEW FEATURES

- could flip the x-axis if needed

BUG FIXES

- improve efficiency of interactiveViewer

Changes in version 1.3.2:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- fix bugs in importScore when format=="BigWig"

Changes in version 1.3.1:

NEW FEATURES

- add new function plotGRanges for quick view GRanges data

- add ranges when importScore

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

unifiedWMWqPCR
--------------

Changes in version 1.3.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated the citation information.

- Updated the help page of the package.

VariantAnnotation
-----------------

Changes in version 1.14.0:

NEW FEATURES

- gVCF support: - missing END header written out with writeVcf() -
  expand() handles <NON_REF> 'REF' value

- support 'Type=Character' in INFO header fields

- add 'row.names' argument to expand()

- add 'Efficient Usage' section to readVcf() man page

- efficiency improvements to info(..., row.names=) - anyDuplicated()
  less expensive than any(duplicated()) - use row.names=FALSE when not
  needed, e.g., show()

- add genotypeCodesToNucleotides()

- add support for gvcf in isSNP family of functions

- add VcfFile and VcfFileList classes

- support 'Type=Character' of unspecified length (.)

- add isDelins() from Robert Castelo

- add makeVRangesFromGRanges() from Thomas Sandmann

MODIFICATIONS

- VCFHeader support: - SAMPLE and PEDIGREE header fields are now parsed
  - meta(VCFHeader) returns DataFrameList instead of DataFrame -
  show(VCFHeader) displays the outer list names in meta -
  fixed(VCFHeader) returns 'ALT' and 'REF' if present

- 'ALT' in expandedVCF output is DNAStringSet, not *List

- remove .listCumsum() and .listCumsumShifted() helpers

- add multiple INFO field unit test from Julian Gehring

- add additional expand() unit tests

- modify readVccfAsVRanges() to use ScanVcfParam() as the 'param';
  deprecate VRangesScanVcfParam

- replace mapCoords() with mapToTranscripts()

- change 'CDSID' output from integer to IntegerList in locateVariants()
  and predictCoding()

- add readVcf,character,ANY,ANY; remove
  readVcf,character,ANY,ScanVcfParam

- replace rowData() accessor with rowRanges()

- replace 'rowData' argument with 'rowRanges' (construct SE, VCF
  classes)

- replace getTranscriptSeqs() with extractTranscripts()

BUG FIXES

- readVcf() properly handles Seqinfo class as 'genome'

- allow 'ignore.strand' to pass through mapCoords()

- writeVcf() no longer ignores rows with no genotype field

- expand() properly handles - less than all INFO fields are selected -
  VCF has only one row - only one INFO column

- don't call path() on non-*File objects

- split (relist) of VRanges now yields a CompressedVRangesList

- predictCoding() now ignores zero-width ranges

VariantFiltering
----------------

Changes in version 1.4:

USER VISIBLE CHANGES

- Genome sequence is not fixed to a specific human genome anymore and
  can be parametrized via a BSgenome package

- The calculation of splice site strength on variants is now
  parametrized (the user can feed his/her own splice site matrices),
  optional and, by default, disabled because is quite computationally
  intensive and the user must be sure he/she wants to calculate it.

- PhastConsDb and MafDb support has been updated. Older PhastConsDb and
  MafDb annotation packages will likely don't work with this newer
  version of VariantFilering. The user must upgrade also those
  annotation packages.

- Analysis functions unrelatedIndividuals(),
  autosomalRecessiveHomozygous(), etc. now can stream over the input
  VCF file to reduce the memory footprint.

- Multiple single-sample VCF files are no longer supported. The user
  must provide always a single multi-sample VCF file.

- Internally, the package uses now the VRanges-class to store variants
  and they can be now also retrieved in this container.

- Indels are now annotated separtely as Insertions and Deletions. Added
  the annotation of Delins (deletion followed by insertion).

- Added support for annotating non-SNVs using XtraSNPlocs.* packages.

- Added support for associating BAM files to the resulting
  VariantFilteringResults object.

- Added some basic variant visualization capabilities via de Gviz
  package and the support for associated BAM files.

- Added HGVS annotations at genomic, coding and protein level.

- Added annotations to Sequence Ontology terms.

- Changed the display of VariantFilteringResults objects and added a
  summary() method that returns a data.frame tallying variants to
  feature annotations, which can be either BioC/VariantAnnotation
  features or Sequence Ontology features.

- Updated the example VCF files of the CEU trio to include a smaller
  subset of variants but called with the latest version of the GATK
  pipeline. Added a subset of the corresponding BAM files with aligned
  reads around the 1Kb region of some of the called variants.

- Integration of the IRanges/FilterRules mechanism into the filtering
  functionality, enabling the addition of user-specific filters.

BUG FIXES

- Some bugfixes related to matching the sequence style between
  variants, annotations and genome sequence.

xcms
----

Changes in version 1.43.3:

BUG FIXES

- Give a more verbose error message when file not found

Changes in version 1.43.2:

BUG FIXES

- Use ProtGenerics, adapted xcms peaks()

Changes in version 1.43.1:

NEW FEATURE

- function plotQC() for plotting various QC plots on RT and m/z

xps
---

Changes in version 3.2:

VERSION xps-1.27.1

- add SystemRequirements GNU make



Packages removed since the last release.
=================================

The following packages are no longer in Bioconductor:

asmn, COPDSexualDimorphism, DNaseR, flowFlowJo, flowPhyto


