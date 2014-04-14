Bioconductors:

We are pleased to announce Bioconductor 2.14, consisting of 824
software packages, 200 experiment data packages, and more than 860
up-to-date annotation packages. 

There are 77 new software packages, and many updates and improvements
to existing packages; Bioconductor 2.14 is compatible with R 3.1.0,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release includes an updated Bioconductor [Amazon Machine Image]
(http://bioconductor.org/help/bioconductor-cloud-ami/).

Visit [http://bioconductor.org](http://bioconductor.org)
for details and downloads.

Contents
--------

* Getting Started with Bioconductor 2.14
* New Software Packages
* NEWS from new and existing packages
* Packages removed from the release

Getting Started with Bioconductor 2.14
======================================

To update to or install Bioconductor 2.14:

1. Install R 3.1.0.  Bioconductor 2.14 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 77 new packages in this release of Bioconductor.

ABSSeq - Inferring differential expression genes by absolute
expression differences between two groups, utilizing
generalized Poisson model to account for over-dispersion across
samples and heterogeneity of differential expression across
genes.

alsace - Alternating Least Squares (or Multivariate Curve
Resolution) for analytical chemical data, in particular
hyphenated data where the first direction is a retention time
axis, and the second a spectral axis. Package builds on the
basic als function from the ALS package and adds functionality
for high-throughput analysis, including definition of time
windows, clustering of profiles, retention time correction,
etcetera.

asmn - Performs all sample mean normalization using raw data
output from BeadStudio and MethyLumiM data.

ASSIGN - ASSIGN is a computational tool to evaluate the pathway
deregulation/activation status in individual patient samples.
ASSIGN employs a flexible Bayesian factor analysis approach
that adapts predetermined pathway signatures derived either
from knowledge-based literatures or from perturbation
experiments to the cell-/tissue-specific pathway signatures.
The deregulation/activation level of each context-specific
pathway is quantified to a score, which represents the extent
to which a patient sample encompasses the pathway
deregulation/activation signature.

AtlasRDF - Query the Gene Expression Atlas RDF data at the European
Bioinformatics Institute using genes, experimental factors (such as
disease, cell type, compound treatments), pathways and proteins. Also
contains a function to perform an enrichment of your gene list across
Experimental Factor Ontology (EFO) using the Atlas background set.

Basic4Cseq - Basic4Cseq is an R/Bioconductor package for basic
filtering, analysis and subsequent visualization of 4C-seq
data. Virtual fragment libraries can be created for any
BSGenome package, and filter functions for both reads and
fragments and basic quality controls are included. Fragment
data in the vicinity of the experiment's viewpoint can be
visualized as a coverage plot based on a running median
approach and a multi-scale contact profile.

BEAT - Model-based analysis of single-cell methylation data

BiocCheck - Bioconductor-specific package checks

biosvd - The biosvd package contains functions to reduce the input
data set from the feature x assay space to the reduced diagonalized
eigenfeature x eigenassay space, with the eigenfeatures and
eigenassays unique orthonormal superpositions of the features and
assays, respectively. Results of SVD applied to the data can
subsequently be inspected based on generated graphs, such as a heatmap
of the eigenfeature x assay matrix and a bar plot with the
eigenexpression fractions of all eigenfeatures. These graphs aid in
deciding which eigenfeatures and eigenassays to filter out (i.e.,
eigenfeatures representing steady state, noise, or experimental
artifacts; or when applied to the variance in the data, eigenfeatures
representing steady-scale variance). After possible removal of steady
state expression, steady-scale variance, noise and experimental
artifacts, and after re-applying SVD to the normalized data, a summary
html report of the eigensystem is generated, containing among others
polar plots of the assays and features, a table with the list of
features sortable according to their coordinates, radius and phase in
the polar plot, and a visualization of the data sorted according to
the two selected eigenfeatures and eigenassays with colored
feature/assay annotation information when provided. This gives a
global picture of the dynamics of expression/intensity levels, in
which individual features and assays are classified in groups of
similar regulation and function or similar cellular state and
biological phenotype.

CAFE - Detection and visualizations of gross chromosomal
aberrations using Affymetrix expression microarrays as input

ccrepe - The CCREPE (Compositionality Corrected by REnormalizaion and
PErmutation) package is designed to assess the significance of general
similarity measures in compositional datasets.  In microbial abundance
data, for example, the total abundances of all microbes sum to one;
CCREPE is designed to take this constraint into account when assigning
p-values to similarity measures between the microbes.  The package has
two functions: ccrepe: Calculates similarity measures, p-values and
q-values for relative abundances of bugs in one or two body sites
using bootstrap and permutation matrices of the data. nc.score:
Calculates species-level co-variation and co-exclusion patterns based
on an extension of the checkerboard score to ordinal data.

ChIPQC - Quality metrics for ChIPseq data

ChIPseeker - This package implements functions to retrieve  the
nearest genes around the peak, annotate genomic region of the peak.
Visualization functions are implemented to summarize genomic
annotation, distance to TSS, and overlap of peaks or genes.

Clomial - Clomial fits binomial distributions to counts obtained
from Next Gen Sequencing data of multiple samples of the same
tumor. The trained parameters can be interpreted to infer the
clonal structure of the tumor.

CNEr - Large-scale identification and advanced visualization of sets
of conserved noncoding elements.

COHCAP - This package provides a pipeline to analyze
single-nucleotide resolution methylation data (Illumina 450k
methylation array, targeted BS-Seq, etc.). It provides QC
metrics, differential methylation for CpG Sites, differential
methylation for CpG Islands, integration with gene expression
data, and visualization of methylation values.

COMPASS - COMPASS is a statistical framework that enables unbiased
analysis of antigen-specific T-cell subsets. COMPASS uses a Bayesian
hierarchical framework to model all observed cell-subsets and select
the most likely to be antigen-specific while regularizing the small
cell counts that often arise in multi-parameter space. The model
provides a posterior probability of specificity for each cell subset
and each sample, which can be used to profile a subject's immune
response to external stimuli such as infection or vaccination.

compcodeR - This package provides extensive functionality for
comparing results obtained by different methods for differential
expression analysis of RNAseq data. It also contains functions for
simulating count data and interfaces to several packages for
performing the differential expression analysis.

CompGO - This package contains functions to accomplish several
tasks. It is able to download full genome databases from UCSC,
import .bed files easily, annotate these .bed file regions with
genes (plus distance) from aforementioned database dumps,
interface with DAVID to create functional annotation and gene
ontology enrichment charts based on gene lists (such as those
generated from input .bed files) and finally visualise and
compare these enrichments using either directed acyclic graphs
or scatterplots.

COPDSexualDimorphism - Sexual dimoprhic and COPD differential  (SDCD)
analysis contrasts regression coefficients from two stratified
analysis. Stratification can be done in two ways: by COPD status or by
sex. For COPD-stratified analysis, SDCD analysis contrasts sexual
dimorphism between cases and controls, while sex-stratified SDCD
analsysis contrasts COPD differential expression pattern between males
and females. The package is meant to be used in conjunction with the
package limma.

CopyNumber450k - This package contains a set of functions that allow
CNV calling from Illumina 450k methylation microarrays.

CoverageView - This package provides a framework for the visualization
of genome coverage profiles. It can be used for ChIP-seq experiments,
but it can be also used for genome-wide nucleosome positioning
experiments or other experiment types where it is important to have a
framework in order to inspect how the coverage distributed across the
genome

CRISPRseek - The package includes functions to find potential guide
RNAs for input target sequences, optionally filter guide RNAs without
restriction enzyme cut site, or without paired guide RNAs, genome-wide
search for off-targets, score, rank, fetch flank sequence and indicate
whether the target and off-targets are located in exon region or not.
Potential guide RNAs are annotated with total score of the top5 and
topN off-targets, detailed topN mismatch sites, restriction enzyme cut
sites, and paired guide RNAs. This package leverages Biostrings and
BSgenome packages.

DMRcate - De novo identification and extraction of differentially
methylated regions (DMRs) in the human genome using Illumina Infinium
HumanMethylation450 BeadChip array data. Provides functionality for
filtering probes possibly confounded by SNPs and cross-hybridisation.
Includes bedGraph and plotting functions.

DMRforPairs - DMRforPairs allows researchers to compare n>=2 unique
samples with regard to their methylation profile. The (pairwise)
comparison of n unique single samples distinguishes DMRforPairs from
other existing pipelines as these often compare groups of samples in
either single CpG locus or region based analysis. DMRforPairs defines
regions of interest as genomic ranges with sufficient probes located
in close proximity to each other. Probes in one region are optionally
annotated to the same functional class(es). Differential methylation
is evaluated by comparing the methylation values within each region
between individual samples and (if the difference is sufficiently
large), testing this difference formally for statistical significance.

dualKS - This package implements a Kolmogorov Smirnov rank-sum based
algorithm for training (i.e. discriminant analysis--identification
of genes that discriminate between classes) and classification
of gene expression data sets.  One of the chief strengths of
this approach is that it is amenable to the "multiclass" problem.
That is, it can discriminate between more than 2 classes.

EDDA - EDDA can aid in the design of a range of common experiments
such as RNA-seq, Nanostring assays, RIP-seq and Metagenomic
sequencing, and enables researchers to comprehensively investigate the
impact of experimental decisions on the ability to detect differential
abundance.

ELBOW - Elbow an improved fold change test that uses cluster
analysis and pattern recognition to set cut off limits that are
derived directly from intrareplicate variance without assuming
a normal distribution for as few as 2 biological replicates.
Elbow also provides the same consistency as fold testing in
cross platform analysis. Elbow has lower false positive and
false negative rates than standard fold testing when both are
evaluated using T testing and Statistical Analysis of
Microarray using 12 replicates (six replicates each for initial
and final conditions). Elbow provides a null value based on
initial condition replicates and gives error bounds for results
to allow better evaluation of significance.

fastLiquidAssociation - This package extends the function of the
LiquidAssociation package for genome-wide application. It integrates a
screening method into the LA analysis to reduce the number of triplets
to be examined for a high LA value and provides code for use in
subsequent significance analyses.

flowBin - Software to combine flow cytometry data that has been
multiplexed into multiple tubes with common markers between them, by
establishing common bins across tubes in terms of the common markers,
then determining expression within each tube for each bin in terms of
the tube-specific markers.

flowCL - Semantic labelling of flow cytometric cell populations.

flowCyBar - A package to analyze flow cytometric data using gate
information to follow population/community dynamics

flowMatch - Matching cell populations and building meta-clusters and
templates from a collection of FC samples.

FRGEpistasis - A Tool for Epistasis Analysis Based on Functional
Regression Model

gaucho - Use genetic algorithms to determine the relationship
between clones in heterogenous populations such as cancer
sequencing samples

GeneOverlap - Test two sets of gene lists and visualize the results.

geneRxCluster - Detect Differential Clustering of Genomic Sites such
as gene therapy integrations.  The package provides some functions for
exploring genomic insertion sites originating from two different
sources. Possibly, the two sources are two different gene therapy
vectors.  Vectors are preferred that target sensitive regions less
frequently, motivating the search for localized clusters of insertions
and comparison of the clusters formed by integration of different
vectors.  Scan statistics allow the discovery of spatial differences
in clustering and calculation of False Discovery Rates (FDRs)
providing statistical methods for comparing retroviral vectors. A scan
statistic for comparing two vectors using multiple window widths to
detect clustering differentials and compute FDRs is implemented here.

GenomeInfoDb - Contains data and functions that define and allow
translation between different chromosome sequence naming conventions
(e.g., "chr1" versus "1"), including a function that attempts to place
sequence names in their natural, rather than lexicographic, order.

GenomicAlignments - Provides efficient containers for storing and
manipulating short genomic alignments (typically obtained by aligning
short reads to a reference genome). This includes read counting,
computing the coverage, junction detection, and working with the
nucleotide content of the alignments.

GenomicFiles - This package provides infrastructure for parallel
queries distributed 'by file' or 'by range'. User defined
map and reduce functions provide added flexibility for
data combination and manipulation.

GOTHiC - This is a Hi-C analysis package using a cumulative
binomial test to detect interactions between distal genomic
loci that have significantly more reads than expected by chance
in Hi-C experiments. It takes mapped paired NGS reads as input
and gives back the list of significant interactions for a given
bin size in the genome.

GSCA - GSCA takes as input several lists of activated and
repressed genes. GSCA then searches through a compendium of
publicly available gene expression profiles for biological
contexts that are enriched with a specified pattern of gene
expression. GSCA provides both traditional R functions and
interactive, user-friendly user interface.

iClusterPlus - Integrative clustering of multiple genomic data using a
joint latent variable model

INPower - An R package for computing the number of susceptibility
SNPs and power of future studies

massiR - Predicts the sex of samples in gene expression microarray
datasets

MeSHDbi - The package is unified implementation of MeSH.db,
MeSH.AOR.db, and MeSH.PCR.db and also is interface to construct Gene-
MeSH package (org.MeSH.XXX.db). loadMeSHDbiPkg import sqlite file and
generate org.MeSH.XXX.db.

meshr - A set of annotation maps describing the entire MeSH assembled
using data from MeSH

messina - Messina is a collection of algorithms for constructing
optimally robust single-gene classifiers, and for identifying
differential expression in the presence of outliers or unknown
sample subgroups.  The methods have application in identifying
lead features to develop into clinical tests (both diagnostic
and prognostic), and in identifying differential expression
when a fraction of samples show unusual patterns of expression.

metaMS - MS-based metabolomics data processing and compound
annotation pipeline.

metaseqR - Provides an interface to several normalization and
statistical testing packages for RNA-Seq gene expression data.
Additionally, it creates several diagnostic plots, performs
meta-analysis by combinining the results of several statistical
tests and reports the results in an interactive way.

MIMOSA - Modeling count data using Dirichlet-multinomial and beta-
binomial mixtures with applications to single-cell assays.

Mirsynergy - Detect synergistic miRNA regulatory modules by
overlapping neighbourhood expansion.

MLSeq - This package applies several machine learning methods, including
SVM, bagSVM, Random Forest and CART, to RNA-Seq data.

mmnet - This package gives the implementations microbiome
metabolic network constructing and analyzing. It introduces a
unique metagenomic systems biology approach, mapping
metagenomic data to the KEGG global metabolic pathway and
constructing a systems-level network. The system-level network
and the next topological analysis will be of great help to
analysis the various functional properties, including
regulation and metabolic functionality of the metagenome.

NetPathMiner - NetPathMiner is a general framework for network path
mining using genome-scale networks. It constructs networks from
KGML, SBML and BioPAX files, providing three network
representations, metabolic, reaction and gene representations.
NetPathMiner finds active paths and applies machine learning
methods to summarize found paths for easy interpretation. It
also provides static and interactive visualizations of networks
and paths to aid manual investigation.

nondetects - Methods to model and impute non-detects in the results of
qPCR experiments.

npGSEA - Current gene set enrichment methods rely upon permutations
for inference.  These approaches are computationally expensive and
have minimum achievable p-values based on the number of permutations,
not on the actual observed statistics.  We have derived three
parametric approximations to the permutation distributions of two gene
set enrichment test statistics.  We are able to reduce the
computational burden and granularity issues of permutation testing
with our method, which is implemented in this package. npGSEA
calculates gene set enrichment statistics and p-values without the
computational cost of permutations.  It is applicable in settings
where one or many gene sets are of interest.  There are also built-in
plotting functions to help users visualize results.

PECA - Calculates Probe-level Expression Change Averages (PECA) to
identify differential expression in Affymetrix gene expression
microarray studies or in proteomic studies using peptide-level
mesurements respectively.

PhenStat - Package contains methods for statistical analysis of
phenotypic data such as Mixed Models and Fisher Exact Test.

QDNAseq - Quantitative DNA sequencing for chromosomal aberrations.

Rariant - The 'Rariant' package identifies single nucleotide variants
from sequencing data based on the difference of binomially distributed
mismatch rates between matched samples.

Rcpi - The Rcpi package offers an R/Bioconductor package
emphasizing the comprehensive integration of bioinformatics and
chemoinformatics into a molecular informatics platform for drug
discovery.

RefNet - Molecular interactions with metadata, some archived, some
dynamically obtained

roar - Identify preferential usage of APA sites, comparing two
biological conditions, starting from known alternative sites
and alignments obtained from standard RNA-seq experiments.

rpx - This package implements an interface to proteomics
data submitted to the ProteomeXchange consortium.

sangerseqR - This package contains several tools for analyzing Sanger
Sequencing data files in R, including reading .scf and .ab1
files, making basecalls and plotting chromatograms.

sapFinder - sapFinder is developed to automate (1)
variation-associated database construction, (2) database searching,
(3) post-processing, (4) HTML-based report generation in shotgun
proteomics.

savR - Parse Illumina Sequence Analysis Viewer (SAV) files,
access data, and generate QC plots.

scsR - Corrects genome-wide siRNA screens for seed mediated
off-target effect. Suitable functions to identify the effective
seeds/miRNAs and to visualize their effect are also provided in
the package.

SomaticSignatures - The SomaticSignatures package identifies
mutational signatures of single nucleotide variants (SNVs).

Sushi - Flexible, quantitative, and integrative genomic
visualizations for publication-quality multi-panel figures

TitanCNA - Hidden Markov model to segment and predict regions of
subclonal copy number alterations (CNA) and loss of
heterozygosity (LOH), and estimate cellular prevalenece of
clonal clusters in tumour whole genome sequencing data.

trackViewer - plot ChIP-seq, RNA-seq, miRNA-seq, DNA-seq and etc NGS
sequence data, especially for big files.

UNDO - UNDO is an R package for unsupervised deconvolution of
tumor and stromal mixed expression data. It detects marker
genes and deconvolutes the mixing expression data without any
prior knowledge.

unifiedWMWqPCR - This packages implements the unified Wilcoxon-Mann-
Whitney Test for qPCR data. This modified test allows for testing
differential expression in qPCR data.

VariantFiltering - Filter genetic variants using different criteria
such as inheritance model, amino acid change consequence, minimum
allele frequencies across human populations, splice site strength,
conservation, etc.

viper - Inference of protein activity from gene expression data,
including the VIPER and msVIPER algorithms


NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages. The following package NEWS is available:


ADaCGH2
-------

Changes in version 2.3.10 (2013-12-27):

- Minor changes to main vignette with a new section "Why ADaCGH2 instead of a manual solution".

- Change in NAMESPACE to adapt to changes in ffbase or bit (we were getting warnings of "replacing
  previous import by ffbase::[.ff when loading ADaCGH2"")

Changes in version 2.3.9 (2013-11-28):

- Minor changes to "benchmarks.pdf": consistent usage of lty and pch for figures of reading and
  analysis.

Changes in version 2.3.8 (2013-11-26):

- Minor changes to "benchmarks.pdf".

Changes in version 2.3.7 (2013-11-26):

- Minor changes to "benchmarks.pdf".

Changes in version 2.3.6 (2013-11-24):

- Help for inputToADaCGH has a more verbose section on the need to use the right mc.cores.

- Made default for mc.cores to inputToADaCGH be half the number of cores.

- More changes to benchmarks.pdf, including the link to all the data.

Changes in version 2.3.5 (2013-11-11):

- Changes in vignette: unified all benchmarking in benchmark.pdf.

- Clarified help for pSegment on loadBalance and say explicitly not the default for HaarSeg.

- Note versions 2.3.3 and 2.3.4 were never placed in BioC repos.

Changes in version 2.3.4 (2013-11-09):

- Reorganized file location for additional files (benchmarks, long-vignette).

- Addedd loadBalance as an argument to the pSegment and pChromPlot functions, so that the user can
  choose to use load balancing-like with both forking and MPI. See additional benchmarking vignette
  for tables with results. Added code to the long vignette to exercise this.

Changes in version 2.3.3:

- Playing with clusterApplyLB.

Changes in version 2.3.2 (2013-10-22):

- Import aCGH and snapCGH again. Corresponding changes in Author, init.c, R code, and C code:
  basically, all the code taken from those packages is not here.

Changes in version 2.3.1 (2013-10-20):

- Added aCGH and snapCGH code: new R and C files, with minor modifications and updates (registration
  of routines, no partial matching of arguments, etc)

- NAMESPACE and Description reflect no longer dependency on snapCGH or aCGH.

- This version works, passes tests, is checked to give the same results as previous ones, etc. But
  never made it into BioC repos, as on 2013-10-23 I realized aCGH and snapCGH were again in BioC
  devel.

Changes in version 2.3.0 (2013-10-20):

- Version bump for new BioC devel.

affxparser
----------

Changes in version 1.36.0 (2014-04-11):

- The version number was bumped for the Bioconductor release version, which is now BioC v2.14 for R
  (>= 3.1.0).

Changes in version 1.35.3 (2014-02-28):

- Same updates as in release v1.34.2.

Changes in version 1.35.2 (2014-02-28):

- Patches to Fusion SDK based on clang v3.4.

Changes in version 1.35.1 (2014-02-27):

- Same updates as in release v1.34.1.

Changes in version 1.35.0 (2013-10-14):

- The version number was bumped for the Bioconductor devel version.

AllelicImbalance
----------------

Changes in version 1.2.0:

NEW FEATURES

- Faster getAlleleCounts()

BUG FIXES

- Fixed bug in the "realCigarPosition" function, that otherwise could generate an error.

alsace
------

Changes in version 0.2:

- implement overlapping segments! should be a real novelty _DONE_, Oct 24

- add quality criteria in summary _DONE_, Oct 20

- complete pipeline so that the output is an intensity matrix _DONE_

- add wrapper function for a complete analysis _DONE_ Nov 20

- add some parallellization support _IN PROGRESS_: some lapply statements have been changed into
  mclapply but the most time-consuming step, the ALS iterations, will remain slow. Different time
  windows are done in parallel, but some time windows require very few iterations, others require
  quite a lot. In addition, the ALS output is annoying. Also the warping can be slow with many
  components.

- add estimation of # of components

annmap
------

Changes in version 1.5.10:

- IMPROVED ANNMAP-121 Allow setting of WS base url for local testing

Changes in version 1.5.8:

- IMPROVED ANNMAP-119 as.vector=TRUE results are now named

- IMPROVED ANNMAP-120 genomeToTranscriptCoords could handle coding regions

Changes in version 1.5.7:

- BUG Check for uniqueness in transcriptToCodingExon, transcriptToCodingRange, transcriptToUtrExon
  and transcriptToUtrRange

Changes in version 1.5.6:

- NEW Added transcriptToUtrExon

- IMPROVED Sequence removed from transcriptToCodingExon

Changes in version 1.5.5:

- IMPROVED ANNMAP-115 Added method geneToGeneRegionTrack to generate a list of Gviz GeneRegionTrack
  elements that can be passed to plotTracks

Changes in version 1.5.4:

- BUG ANNMAP-118 transcriptToCodingRange( ..., end='3' ) returns whole transcript for untranslated
  transcripts.

- IMPROVED ANNMAP-114 Added additional logging when using webservice

- IMPROVED ANNMAP-116 Introduce function(s) to calculate the coding length of a gene(s) - See:
  ?nonIntronicGeneLength and ?nonIntronicTranscriptLength.

Changes in version 1.5.3:

- IMPROVED ANNMAP-114 Added additional logging when using webservice

- IMPROVED ANNMAP-113 Added documentation to cookbook for webservice

- IMPROVED webservice json to data.frame speed improvements

AnnotationDbi
-------------

Changes in version 1.26:

NEW FEATURES and API changes

- Support for new schemaless chip packages to allow more flexible suport of non-model organisms in
  chip packages.

- Support for Inparanoid8 objects which contain data from the most recent release of inparanoid.
  These objects will be available through the AnnotationHub.

BUG FIXES AND CODE MAINTENANCE

- PFAM and PROSITE mappings are going away!  This is because these mappings are tied to the now long
  defunct IPI ids.  Don't worry we have you covered.  These ids will still be available with
  select(), the mappings are just bad because they can't exist without IPI ids (which we can no
  longer get).  So for the next release if you try to use these mappings you will be warned about
  it.  After that, you need to have made the switch to select(). * * * 1.24.x SERIES NEWS * * *

AnnotationHub
-------------

Changes in version 1.4.0:

NEW FEATURES

- Add display,AnnotationHub-method using shiny

- Adds support for Inparanoid8 objects (which in turn support select, keys, columns and keytypes)

- Adds local metadata caching that dramatically improves performance

SIGNIFICANT USER VISIBLE CHANGES

- metadata argument 'columns' replaces deprecated 'cols'

BUG FIXES

- metadata() returns only information on current hub elements

aroma.light
-----------

Changes in version 2.0.0 (2014-04-11):

- The version number was bumped for the Bioconductor release version, which is now BioC v2.14 for R
  (>= 3.1.0).

Changes in version 1.99.3 (2014-03-31):

- Bumped the version such that the next release will be 2.0.0.

BaseSpaceR
----------

Changes in version 1.7:

NEW FEATURES

- Improved API design and a cleaner representation of the BaseSpace Data Model.

- New count<RESOURCE> methods. These methods should be used to get the number of instances of a
  particular resource which are visible under the current scope.  For example, the total number of
  Samples within a Project. Or the total number of Files for a given AppResults.

- Improved authentication process. One can now seamlessly used the R SDK from with Native Apps or
  Web Apps. New 'AppSessionAuth' class to handle the authentication in this cases. It extends the
  'AppAuth' class, by keeping track of the AppSession Id.  User-level constructors:
  'authWebClient()' and 'authNativeClient()' generator .AppSessionAuth().

- New mechanism for handling the BaseSpace server responses. 'ResponseStatus' class is now used for
  tracking the HTTP/server response status and integrated in the 'AppAuth' interface.
  'x$showResponse()' can now be used to print the (JSON) body of the last response from the server.

- Integration with some Bioconductor data structures. BAM data from BaseSpace can now be easily map
  to 'BamFile' objects.

SIGNIFICANT USER-VISIBLE CHANGES

- 'authenticateClient()' function remove and replaced by 'authWebClient()' and 'authNativeClient()'.
  'AppAuth' instances are used for desktop Apps only.

BUG FIXES

- Minor bug fixes with the authentication process.

- Method dispatchment now works as intended for 'listFiles' and 'Files'.

BiGGR
-----

Changes in version 1.1.4:

- Updated vignette. Equality constraints on GABA shunt and pentose phosphate pathway fluxes are now
  used in FBA and ensemble simulations.

Changes in version 1.1.3:

- Updated vignette and bibtex file. FBA and ensemble sampling are now done with integration of brain
  metabolite uptake measurement data.

- Added new function 'sampleFluxEnsemble' which can be used to generate a posterior distribution of
  flux vectors with respect to measuerement error in certain fluxes. An example of
  'sampleFluxEnsemble' has been added to the package vignette.

- Added LIM model file 'Glycolysis_TCA.LIM' to the examples directory. This file contains the model
  which is build in the examples.

- Fixed warnings in functions buildSBMLFromReactionIDs and createLIMFromSBML to seperate missing
  identifiers with a comma.

BiocInstaller
-------------

Changes in version 1.14.0:

NEW FEATURES

- biocUpdatePackages updates installed packages and their dependencies.

BiocParallel
------------

Changes in version 0.5.5:

NEW FEATURES

- multicoreWorkers() determines the number of workers based on operating system (Windows: 1), user
  preference (via the global option options(mc.cores=...)), or system capability (detectCores()).

- bpparam() selects a default BiocParallelParam, from global options or, if that fails, the most
  recently registered() back-end.

SIGNIFICANT USER-VISIBLE CHANGES

- Rename argument controlling resumption on error as BPRESUME

- Default to parallel back-end (multicore on non-Windows; snow on Windows).

BUG FIXES

- bpvec,ANY,MulticoreParam-method with fewer tasks than cores evaluates only the cores for which
  tasks are defined.

Changes in version 0.5.2:

NEW FEATURES

- mclapply(), pvec() require only length, [, and (for mclapply) [[.

- pvectorize() creates a parallel version of its vectorized function argument.

- MulticoreParam, SnowParam, DoparParam (foreach-derived), SerialParam to parameterize back-ends.

- bplapply, bpvec as parallel evaluation models.

- bpstart, bpstop, bpisup for back-end management.

- bpvec has a new argument AGGREGATE, a function to specify how results are to be combined.

- Support for BatchJobs back-end added, via GSOC Michel Lang.

SIGNIFICANT USER-VISIBLE CHANGES

- BPPARM is now used as the argument name for passing BiocParallelParam instances to functions.

- bplapply and bpvec now only dispatch on X and BPPARAM.

BiocStyle
---------

Changes in version 1.2.0:

USER VISIBLE CHANGES

- Remove dependency on 'helvet' LaTeX package to allow for the same font both in Sweave and knitr

- Improve package vignette by adding paragraphs about building vignettes and using bibliography

- knitr chunk_opts error=FALSE by default, so failures during vignette processing are signaled to R
  CMD build and R CMD check

- Name-mangle \comment mark-up to avoid conflicts with other LaTeX styles

- Introduce \bioctitle to allow for full and short (header) titles

- Add BiocStyle::latex option 'use.unsrturl=TRUE' to use the 'unsrturl' bibliography style by
  default

BUG FIXES

- Avoid use of 'titling' LaTeX package to circumvent the conflict w/ \footnote in \author

biomvRCNS
---------

Changes in version 1.3.3:

NEW FEATURES

- updates in package vignette

BUG FIXES

- remove limitation of maximum typecode in the plot, now using rainbow()

Changes in version 1.3.2:

BUG FIXES

- fix the returned posterior p

bsseq
-----

Changes in version 0.11:

- Converted to using Authors@R in the DESCRIPTIOn file.

- plotRegion did not respect the col/lty/lwd argument if given explicitely as opposed to through
  pData().  Reported by Karen Conneely <kconnee@emory.edu>.

- Fixed an issue introduced by the previous change (to plotRegion). Reported (with a fix) by Tim
  Triche Jr <tim.triche@gmail.com>.

- Fixed a serious bug in collapseBSseq reported by Julien Roux <jroux@uchicago.edu>: it would use
  the Meth values instead of Cov for the final object's Cov.  However, this will result in the
  return object having a methylation of 100 percent across all loci, so hopefully users will have
  seen this for themselves.

- Fixed a bug in combineList which made combineList use "slow" code even in a special case where a
  faster shortcut was possible.  This bug fix does not change the output of the function under any
  input, it only makes it faster.  Reported by Julien Roux <jroux@uchicago.edu>.

bumphunter
----------

Changes in version 1.3.5:

- nearestgene() now loads TT.rda into the local environment

- known_transcripts(), and consequently matchGenes(), was yielding mostly NA for the Refseq /
  'annotation' column (iff there were multiple refseq ids for the given Entrez/Gene id)

- Fixed data/TT.rda, as the new value of known_transcripts()

Changes in version 1.3:

- Fixed NAMESPACE issues.

CAFE
----

Changes in version 0.99.2:

- Code and man files are now limited to 80 character width per line

- Code now uniformly uses 4-space indentation

- Vignette now uses BiocStyle style

Changes in version 0.99.1:

- fixes a bug that gave an error when plotting with idiogram=TRUE.  This was due to an update in the
  ggbio package that broke CAFE

- CAFE can now be used with species other than human It is no longer limited to 22 autologous
  chromosomes

CAMERA
------

CHANGES IN VERSION 1.19.1

  * Fixed bug in the generateRules function, where in the negative mode
    rules with polycharged positve ions are not created


CellNOptR
---------

Changes in version 1.10.0 (2014-03-13):

- BUG FIX
  
  * typos leading to failure fixed in cutCNOlist
  
  * typos leading to failure in CNOlist if variance were provided

- CHANGES
  
  * Changed warning into errors in readMIDAS and checkSignals functions

- NEW FEATURES
  
  * self loop can be read from the SIF files
  
  * readSBMLQual do not need to add dummy nodes anymore

ChemmineOB
----------

Changes in version 1.2.0:

NEW FEATURES

- SMARTS search

ChemmineR
---------

Changes in version 2.16.0:

NEW FEATURES

- SMARTS Search availible through ChemmineOB

- Folding of FPset objects

- Support for SQL database updates by compound name

- Coordinate re-generation via ChemmineOB improves structure rendering

chimera
-------

Changes in version 1.5.2:

NEW FEATURES

- prettyPrint function was added to generate a table delimited file output form a list of fSet
  objects.

ChIPpeakAnno
------------

Changes in version 2.11.4:

NEW FEATURE

- Add selection of 'shortestDistance' to output parameter of annotatePeakInBatch.

Changes in version 2.11.3:

NEW FEATURE

- Improved efficiency of annotatePeakInBatch.

ChIPQC
------

Changes in version 1.0.0:

- Initial release.

cisPath
-------

Changes in version 1.3.3:

- improvements

Changes in version 1.3.1:

- documentation improvements

cleanUpdTSeq
------------

Changes in version 1.1.3:

NEW FEATURES

- No changes classified as "new features" (package under active development)

BUG FIXES

- Fix the bugs 'match' requires vector arguments when using %in%

Changes in version 1.1.2:

NEW FEATURES

- Add citation

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.1.1:

NEW FEATURES

- Re-name the inst/doc directory to vignettes

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.1.0:

NEW FEATURES

- New package released in bioconductor

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

cleaver
-------

Changes in version 1.1.8 (2014-03-26):

- Fix missedCleavages>1.

- Add argument "unique".

- Add methods, Biostrings, IRanges to NAMESPACE.

Changes in version 1.1.7 (2014-03-25):

- Typo in the manual page.

Changes in version 1.1.5 (2014-02-25):

- Using AAStringSetList constructor for list of characters instead of creating a lot of AAStringSets
  dramatically decreases running time for cleave,AAStringSet-method.

Changes in version 1.1.4 (2013-12-20):

- tests:

- move tests into tests/testthat to adapt to testthat 0.8 and new CRAN policy.

Clomial
-------

Changes in version 0.99.0 (2014-02-11):

- Under review by Bioconductor.

- Created.

clonotypeR
----------

Changes in version 1.2.0:

NEW FEATURES

- Added a “long” option to yassai_identifier(), where every amino acid is represented.  This solves
  the problem of ID collisions, where a single identifier could be produced by two different
  clonotypes.

- Added a new option to clonotype_table() for randomly sampling libraries.

- Added a new mode to common_clonotypes() for calculating the abundance relatively to one library.

- Unified the syntax of common_clonotypes() and unique_clonotypes().

BUG FIXES

- Removed infinite loop in yassai_identifier() when the germinal V sequence was completely absent
  from the CDR3.

- Corrected output bug where “integer(0)” was returned if the V–J boundary was a codon boundary.

- yassai_identifier(): properly encode the names of the V and J segments.

clusterProfiler
---------------

Changes in version 1.11.3:

- now plotting compareClusterResult support scale the dot size by geneRatio, and it is now setting
  as default. <2014-01-14, Tue>

- The original default parameter by="percentage" is now changing to "rowPercentage" <2014-01-14,
  Tue>

Changes in version 1.11.2:

- remove viewKEGG function <2013-11-12, Tue>

- In vignette, illustrates how to visualize KEGG pathway using use pathview package. <2013-11-12,
  Tue>

CNAnorm
-------

1.9.2: Added function CNAnormWorkflow Added window weighting (for segmentation) depending on
          dispersion. Fixed a bug that made plotGenome crash if the interval of interest had no
          valid values. Changed vignette to better reflect typical usage and defined basic and
          advanced use.

codelink
--------

Changes in version 1.32.0:

- Minor fixes to documentation.

- Dropped the old term "Bioarray" (originally used for the first Codelink platform) from the
  DESCRIPTION file.

COMPASS
-------

Changes in version 0.9.0:

- First release.

compcodeR
---------

0.99.2: Bug fixes

0.99.1: improved documentation of the compData class

0.2.0: more examples, checks and error messages

0.2.0: more options for providing user-specified values in the data simulation

0.2.0: changed the range of the color palette in the Spearman correlation plots to [-1,1]

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

DESeq2
------

Changes in version 1.4.0:

- *** USAGE NOTE *** Expanded model matrices are now used when betaPrior = TRUE (the default).
  Therefore, level comparison results should be extracted using the 'contrast' argument to the
  results() function. Expanded model matrices produce shrinkage of log fold changes that is
  independent of the choice of base level. Expanded model matrices are not used in the case of
  designs with an interaction term between factors with only 2 levels.

- The order of the arguments 'name' and 'contrast' to the results() function are swapped, to
  indicate that 'contrast' should be used for the standard comparisons of levels against each other.
  Calling results() with no arguments will still produce the same comparison: the fold change of the
  last level of the last design variable over the first level of the last design variable. See
  ?results for more details.

- The DESeq() function will automatically replace count outliers flagged by Cook's distance when
  there are 7 or more replicates. The DESeq() argument 'minReplicatesForReplace' (default 7) is used
  to decide which samples are eligible for automatic replacement. This default behavior helps to
  prevent filtering genes based on Cook's distance when there are many degrees of freedom.

Changes in version 1.3.58:

- Added a list() option to the 'contrast' argument of results(). See examples in ?results.

Changes in version 1.3.24:

- rlogTransformation() gains an argument 'fast', which switches to an approximation of the rlog
  transformation. Speed-up is ~ 2x.

- A more robust estimator for the beta prior variance is used: instead of taking the mean of squared
  MLE betas, the prior variance is found by matching an upper quantile of the absolute value of MLE
  betas with an upper quantile of a zero-centered Normal distribution.

Changes in version 1.3.17:

- It is possible to use a log2 fold change prior (beta prior) and obtain likelihood ratio test
  p-values, although the default for test="LRT" is still betaPrior=FALSE.

Changes in version 1.3.15:

- The DESeq() function will automatically replace count outliers flagged by Cook's distance when
  there are 7 or more replicates. The DESeq() argument 'minReplicatesForReplace' (default 7) is used
  to decide which samples are eligible for automatic replacement. This default behavior helps to
  prevent filtering genes based on Cook's distance when there are many degrees of freedom.

- The results() function produces an object of class 'DESeqResults' which is a simple subclass of
  'DataFrame'. This class allows for methods to be written specifically for DESeq2 results. For
  example, plotMA() can be called on a 'DESeqResults' object.

Changes in version 1.3.12:

- Added a check in nbinomWaldTest which ensures that priors on logarithmic fold changes are only
  estimated for interactions terms, in the case that interaction terms are present in the design
  formula.

Changes in version 1.3.6:

- Reduced the amount of filtering from Cook's cutoff: maximum no longer includes samples from
  experimental groups with only 2 samples, the default F quantile is raised to 0.99, and a robust
  estimate of dispersion is used to calculate Cook's distance instead of the fitted dispersion.

Changes in version 1.3.5:

- New arguments to results(), 'lfcThreshold' and 'alternativeHypothesis', allow for tests of log
  fold changes which are above or below a given threshold.

- plotMA() function now passes ellipses arguments to the results() function.

DEXSeq
------

Changes in version 2014-03-21:

- Major code revisions: the ExonCountSet object was deprecated and substituted by the DEXSeqDataSet
  class, a subclass of the DESeqDataSet.

- DEXseq now uses the DESeq2 package as internal engine and backbone for all analyses

- All functions and methods for the ExonCountSet object were replaced by new functions

- DEXSeq is now better integrated with other Bioconductor packages

- We now use knitr to build the vignette

DiffBind
--------

Changes in version 1.10.0:

- Counting
  
  * New: option to compute summits
  
  * New: option to center peaks with fixed width around summits
  
  * New: scores for summits (height, position) and CPM for TMM values
  
  * New: filter reads by mapping quality (mapQCth)
  
  * New: support for PE bam data using summarizeOverlaps
  
  * Remove: bCalledMask (now always TRUE)
  
  * Change: insertLength to fragmentSize
  
  * Add: fragmentSize can be a vector with a size for each sample
  
  * Change: fragmentSize default is 125 bp

- Plotting
  
  * Change: colors based on CRUK color scheme
  
  * PCA plots
  
  * New: legend
  
  * New: label parameter for adding text labels of points in 2D plot
  
  * Venn diagrams
  
  * New: plot overlaps of differentially bound sites by specifying contrasts, thresholds etc.
  
  * New: able to return overlapping peaksets as GRanges directly
  
  * New: able to generate new DBA object consisting of overlapping peaks
  
  * New: labelAttributes for controlling default labels
  
  * New: default main and sub titles
  
  * Heatmaps
  
  * Fix: don’t plot column vector for attributes where every sample has a different value

- General
  
  * New: add attribute value: DBA_ALL_ATTRIBUTES
  
  * Change: SN (signal/noise) to FRIP (fraction of reads in peaks)
  
  * Change: “Down” to “Loss” and Up” to “Gain”
  
  * Vignette
  
  * Change: vignette uses BiocStyles and dynamically generated figures
  
  * Change: example data based on hg19 instead of hg18
  
  * Change: example reads from bam files instead of bed files
  
  * New: section on using DiffBind and ChIPQC together
  
  * New configuration defaults options (DBA$config):
  
  * Metadata name strings: ID, Tissue, Factor, Condition, Treatment, Caller
  
  * th: significance threshold
  
  * bUsePval
  
  * fragmentSize
  
  * mapQCth: filter reads by mapping quality
  
  * fragments (for summarizeOverlaps)

- Bugs/Issues
  
  * Fix: bRemoveDuplicates had some unpredictable behaviour
  
  * Fix: chrN_random were being counted against chrN
  
  * Disable: tamoxifen_GEO.R doesn’t work after SRA changed format of archived data

easyRNASeq
----------

Changes in version 1.99.3:

- Replicated the bug fixing of version 1.8.8.

Changes in version 1.99.2:

- Added missing imports

Changes in version 1.99.1:

- Fixed a bug in the getBamFileList function that prevented reporting accurate error/warning
  messages.

- Updated the package dependencies.

Changes in version 1.99.0:

- Ported the final changes from the Git repository that make the simpleRNASeq function functional.

- Bumped the version number to 1.99 as these are major changes and hence next Bioc release will have
  easyRNASeq 2.0.0.

Changes in version 1.9.7:

- Same as release 1.8.7 change.

Changes in version 1.9.6:

- Ported changes from the Git repository including further class refinements and unit testing.

Changes in version 1.9.5:

- Automatic version number bump by Bioc

Changes in version 1.9.4:

- Mirrored changes in version 1.8.6

Changes in version 1.9.3:

- Mirrored changes in version 1.8.5

Changes in version 1.9.2:

- Ported version 1.8.4 changes to the development version

- Started to import new classes from the git repository as part of the package main function
  re-implementation

Changes in version 1.9.1:

- Change from Bioc to insert GenomicAlignments

Changes in version 1.9.0:

- No changes, new development version

EBImage
-------

Changes in version 4.6.0:

NEW FEATURES

- 'toRGB' function for convenient grayscale to RGB conversion

- 'luminance' mode in 'channel' for luminance-preserving RGB to grayscale conversion

SIGNIFICANT USER-VISIBLE CHANGES

- Performance improvements to: 'Image', 'is.Image', 'readImage', 'writeImage', 'show', 'normalize',
  'getFrame', 'selectChannel', 'rgbImage', 'colorLabels', 'flip'/'flop'

- Reduced memory footprint of 'readImage'

- When called on an 'Image' object, 'as.Image' returns its argument rather than the
  Grayscale-coerced copy

- 'displayRaster' sets 'par$usr' coordinates to image pixel coordinates easing up plot annotation

BUG FIXES

- 'getFrame', 'getNumberOfFrames' and 'colorLabels' support multi-dimensional images

- Proper handling of multi-dimensional character arrays by the 'Image' constructor

- Fixed 'getFrame' and 'combine' in case of single-channel Color Images

- Fixed color mode check in 'validImageObject'

- Proper 'fg.col' and 'bg.col' color handling in 'tile'

- Updates to documentation

EDASeq
------

Changes in version 1.9:

- Added method plotRLE for Relative Log Expression (RLE) plots.

EDDA
----

Changes in version 0.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Put "NEWS.Rd" in EDDA/.

- Change biocViews to DESCRIPTION file.

- Modify DESCRIPTION and NAMESPACE file.

- move the data to EDDA/inst/extdata.

Changes in version 0.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Put "EDDA.Rnw" in EDDA/vignettes.

- Remove "library("EDDA");" in file man/EDDA-package.Rd.

- Clean up pre-defined tempates in man files.

- Use TRUE instead of T and FALSE instead of F.

- Clean R code by breaking code less 80 characters, removing tabs and using indent by a multiple of
  4 spaces.

NEW FEATURES

- Add biocViews to DESCRIPTION file.

- Add KEYWORDS for each function man file.

- Add "NEWS.Rd" in EDDA/inst/.

Changes in version 0.99.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Submit EDDA

edgeR
-----

Changes in version 3.6.0:

- Improved treatment of fractional counts. Previously the classic edgeR pipeline permitted
  fractional counts but the glm pipeline did not. edgeR now permits fractional counts throughout.

- All glm-based functions in edgeR now accept quantitative observation-level weights. The glm
  fitting function mglmLS() and mglmSimple() are retired, and all glm fitting is now done by either
  mglmLevenberg() or mglmOneWay().

- New capabilities for robust estimation allowing for observation-level outliers. In particular, the
  new function estimateGLMRobustDisp() computes a robust dispersion estimate for each gene.

- More careful calculation of residual df in the presence of exactly zero fitted values for
  glmQLFTest() and estimateDisp(). The new code allows for deflation of residual df for more complex
  experimental designs.

- New function processHairpinReads() for analyzing data from shRNA-seq screens.

- New function sumTechReps() to collapse counts over technical replicate libraries.

- New functions nbinomDeviance() and nbinomUnitDeviance. Old function deviances.function() removed.

- New function validDGEList().

- rpkm() is now a generic function, and it now tries to find the gene lengths automatically if
  available from the annotation information in a DGEList object.

- Subsetting a DGEList object now has the option of resetting to the library sizes to the new column
  sums. Internally, the subsetting code for DGEList, DGEExact, DGEGLM, DGELRT and TopTags data
  objects has been simplified using the new utility function subsetListOfArrays in the limma
  package.

- To strengthen the interface and to strengthen the object-orientated nature of the functions, the
  DGEList methods for estimateDisp(), estimateGLMCommonDisp(), estimateGLMTrendedDisp() and
  estimateGLMTagwiseDisp no longer accept offset, weights or AveLogCPM as arguments. These
  quantities are now always taken from the DGEList object.

- The User's Guide has new sections on read alignment, producing a table of counts, and on how to
  translate scientific questions into contrasts when using a glm.

- camera.DGEList(), roast.DGEList() and mroast.DGEList() now include ... argument.

- The main computation of exactTestByDeviance() now implemented in C++ code.

- The big.count argument has been removed from functions exactTestByDeviance() and
  exactTestBySmallP().

- New default value for offset in dispCoxReid.

- More tolerant error checking for dispersion value when computing aveLogCPM().

- aveLogCPM() now returns a value even when all the counts are zero.

- The functions is.fullrank and nonEstimable are now imported from limma.

eiR
---

Changes in version 1.4.0:

NEW FEATURES

- eiCluster can now cluster subsets of database

- use new features in ChemmineR to store duplicate descriptors only once.

- embedded descriptors now stored in database and the matrix file is written out only as needed by
  LSH to create an index.

UPGRADING

- Database schema changes make this version incompatible with version 1.2 an earlier. Existing
  databases will need to be re-loaded.

ensemblVEP
----------

Changes in version 1.4.0:

NEW FEATURES

- VEPParam is now VIRTUAL with subclasses VECParam67 and VEPParam73

- Add helpers currentVEP() and supportedVEP()

- Add VEPParam75 class

MODIFICATIONS

- Modify VEPParam to support API versions 67 and 73

- Add slots 'scriptPath' and 'version' to VEPParam

- New man page for runtime options

- Update biocViews terms

FGNet
-----

Changes in version 2.0:

NEW FEATURES

- The networks can now be plotted for terms in addition to genes.

- In the report, clicking on the plots allows to see the plot at full size next to the terms table
  (if the screen resolution allows it).

- ToMatrix(): Added arguments 'key' (to choose either genes or terms) and 'removeFiltered' terms.
  Renamed main argument (geneTermSets) to 'results'.

- IntersectionNetwork() shows a warning if there is no intersection. Added argument 'plotAllMg' to
  allow choosing wether to plot unconnected metagroups or not.

- FunctionalNetwork() changed the two main arguments to a single one, which is the raw output from
  toMatrix() (a list with names: c("metagroupsMatrix", "gtSetsMatrix", "filteredOut")). Added
  arguments: 'eColor' to provide the edges color, and 'weighted' to draw the edge width according to
  the number of shared gene-term sets.

- Reports: Added argument 'downloadGOtree' to allow choosing wether to download the go term trees
  png (slower) or just provide the link to the web tool.

BUG FIXES

- functionalNetwork() now correctly writes either "cl" or "mg" in the legend

Changes in version 1.3.1:

NEW FEATURES

- adjMatrix() has been renamed to toMatrix()

- GO png trees are now automatically downloaded when generating the report

- functionalNetwork: metagroup/cluster legend order has been changed to alphabetical

BUG FIXES

- Minor bug fixes

flipflop
--------

Changes in version 1.1.6:

NEW FEATURES

- Add 'expected.counts' in the output objects.  It gives the expected raw count (ie number of mapped
  fragments) for each predicted transcripts.  This information is also available in the output gtf
  file.

Changes in version 1.1.4:

NEW FEATURES

- Add a 'OnlyPreprocess' option for performing only the pre-processing step of the input sam file.
  This step writes two files: one file named 'prefix.instance' and one other named
  'prefix.totalnumread', where 'prefix' is the prefix of the input sam file. The 'prefix.instance'
  file can then be given to the option 'preprocess.instance' and the total number of mapped reads is
  stored in the 'prefix.totalnumread' file.

- Allow to give 'NN' (total number of mapped fragments) even when using the sam file as input. This
  can be used to run flipflop in parallel on parsed sam files with the same NN constant.

USER-LEVEL CHANGES

- Handle '~' in input paths with path.expand function.

- Give a more detailed R output with the read count.

MINOR CHANGES

- Do not write individual Coverage into the pre-processing file anymore ie, the .instance file. (in
  practice comments lines in readgroup.cpp, part 'toStream')

Changes in version 1.1.2:

BUG FIXES

- Switch from GRangeList to a regular list as the number of metadata columns could not vary in
  GRangeList.

fmcsR
-----

Changes in version 1.6.0:

NEW FEATURES

- Added timeout option

- Allow parallel execution for batch computations

FRGEpistasis
------------

Changes in version 0.99.5:

- Added FDR control

Changes in version 0.99.4:

- Completed the man pages

- Added the NEWS file

- Updated the DESCRIPTION file

- Completed the vignette

- Updated the genotype file format

Changes in version 0.99.3:

- Corrected the biocViews

Changes in version 0.99.1:

- Added FRGEpistasis vignette

gage
----

Changes in version 2.12.1:

- remove pathview from the imports list for easier installation and loading.

gCMAP
-----

Changes in version 1.7.2:

- CHANGE: Helper functions KEGG2cmap, reactome2cmap, wiki2cmap and go2cmap have two additional
  parameters, min.size and max.size, restricting the output to gene sets containing members within
  the specified range.

- CHANGE: Method camera_score now outputs correlations in column 'effect'.

- BUGFIX: Correct numbers of annotated genes are now returned for each set.

Changes in version 1.7.1:

- CHANGE: Vignettes are now stored in vignettes directory.

gCMAPWeb
--------

Changes in version 1.3.2:

- BUGFIXMoved the 'Rook' package from 'Imports' to 'Depends'

Changes in version 1.3.1:

- CHANGEVignettes are now stored in vignettes directory.

geNetClassifier
---------------

Changes in version 1.3.1:

NEW FEATURES

- plotAssignments: Optimized speed and changed background color to blue. Added argument pointSize.

- plotExpressionProfiles: Added argument sampleColors.

BUG FIXES

- Added RUnit and BiocGenerics to "Suggests" to avoid R CHECK warnings

- Import methods, ebarrays and emfit in namespace to avoid R CHECK warnings

- plot.GenesNetwork is now completelly equivalent to plotNetwork (same default arguments)

GenomeInfoDb
------------

Changes in version 0.99.14:

SIGNIFICANT USER-VISIBLE CHANGES

- rename: package: Seqnames --> GenomeInfoDb supportedStyles -> genomeStyles makeSeqnameIds -->
  rankSeqlevels (add to export) seqnamesOrder --> orderSeqlevels extractSeqnameSet ->
  extractSeqlevels extractSeqnameSetByGroup -> extractSeqlevelsByGroup findSequenceRenamingMaps -->
  mapSeqlevels seqnamesInGroup --> seqlevelsInGroup seqnamesStyle --> seqlevelsStyle
  "seqnameStyle<-" --> "seqlevelsStyle<-"

Changes in version 0.99.7:

SIGNIFICANT USER-VISIBLE CHANGES

- rename: isSupportedSeqnames -> .isSupportedSeqnames supportedSeqnameStyles ->
  .supportedSeqnameStyles supportedSeqnameMappings -> .supportedSeqnameMappings
  isSupportedSeqnamesStyle -> .isSupportedSeqnamesStyle

Changes in version 0.99.6:

NEW FEATURES

- add new functions() seqnamesInGroup which will take a character vector of chromosomes and return
  the chromosomes specified by the group parameter supplied by the user. The user can also give the
  species and the style.  seqnamesOrder() internally calls Herve's function makeSeqnameIds()

- add seqnameStyles generic and method from GenomicRanges

SIGNIFICANT USER-VISIBLE CHANGES

- rename: testSeqnames -> isSupportedSeqnames

- move SeqnamesStyle generic from GenomicRanges and define a new method which works on a character
  vector.

DEPRECATED AND DEFUNCT

- deprecate listAllSupportedStylesBySpecies(), listAllSupportedSeqnameStyles(), supportedOrganisms()
  supportedSeqnameMappingsWithGroup()

DEPRECATED AND USED INTERNALLY(NOT EXPORTED

- deprecate supportedSeqnameMappings(), supportedSeqnameStyles(),
  isSupportedSeqnamesStyle(),issupportedSeqnames()

Changes in version 0.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- The Seqnames package will have functions which will be moved from AnnotationDbi , GenomicRanges

- List of 9 functions moved from AnnotationDbi supportedSeqnameMappings, findSequenceRenamingMaps,
  supportedSeqnameStyles, supportedSeqnames, extractSeqnameSet, testSeqnames,
  isSupportedSeqnamesStyle, listAllSupportedStylesBySpecies, listAllSupportedSeqnameStyles.

- makeSeqnameIds moved from GenomicRanges

- keepStandardChromosomes moved from GenomicRanges

- rename: keepStandardChromosomes -> keepChromosomes

NEW FEATURES

- added new functions: supportedOrganisms() supportedSeqnameMappingsWithGroup()
  extractSeqnameSetByGroup()

GenomicAlignments
-----------------

Changes in version 0.99:

NEW FEATURES

- coverage,BamFile-method uses yieldSize to iterate through large files.

- coverage,character-method calculates coverage from a BAM file.

GenomicFiles
------------

Changes in version 1.0.0:

NEW FEATURES

- First release of GenomicFiles package

GenomicRanges
-------------

Changes in version 1.16.0:

NEW FEATURES

- Add "subset" method for SummarizedExperiment objects.

- Allow DataFrame in SummarizedExperiment assays.

- Add 'use.mcols' arg (FALSE by default) to the granges(), grglist(), and rglist() generics (a.k.a.
  the range-squeezer generics).

- Add coercion method from GRangesList to RangesList.

- Add score() setter for GRangesList objects.

- findOverlaps(..., type="within") now works on circular chromosomes.

- Add 'ignore.strand' arg to "sort" method for GRanges objects.

- Support subsetting of a named list-like object *by* a GenomicRanges subscript.

- Support sort(granges, by = ~ score), i.e., a formula-based interface for sorting by the mcols.

SIGNIFICANT USER-LEVEL CHANGES

- Move many functionalities to the new GenomicAlignments package: - The GAlignments,
  GAlignmentPairs, and GAlignmentsList classes. - The qnarrow() generic and methods. - The "narrow"
  and "pintersect" methods for GAlignments and GAlignmentsList objects. - The low-level CIGAR
  utilities. - The "findOverlaps" methods for GAlignment* objects. - The summarizeOverlaps() generic
  and methods, and the "Counting reads with summarizeOverlaps" vignette. - findCompatibleOverlaps()
  and countCompatibleOverlaps(). - The findSpliceOverlaps() generic and methods. - The "overlap
  encodings" stuff i.e. the "encodeOverlaps" method for GRangesList objects, flipQuery(),
  selectEncodingWithCompatibleStrand(), isCompatibleWithSplicing(), isCompatibleWithSkippedExons(),
  extractSteppedExonRanks(), extractSpannedExonRanks(), extractSkippedExonRanks(), and
  extractQueryStartInTranscript(), and the "OverlapEncodings" vignette.

- Rename 'with.mapping' arg -> 'with.revmap' in "reduce" methods. The old arg name is still working
  but deprecated.

- Move makeSeqnameIds() function to the new GenomeInfoDb package and rename it rankSeqlevels(). The
  old name is still working but deprecated.

- The "strand" methods now perform stricter checking and are guaranteed to always return a factor
  (or factor-Rle) with the "standard strand levels" and no NAs. Or to fail.

BUG FIXES

- Tweaks and fixes to various "strand" methods: - Methods for character vectors and factors do not
  accept NAs anymore (they raise an error). - Methods for integer and logical vectors map NAs to *
  (instead of NA). - Method for Rle objects now also works on character-, factor-, and integer-Rle
  objects (in addition to logical-Rle objects).

genoset
-------

Changes in version 1.16.0:

NEW FEATURES

- cn2lr now has methods for vector, matrix, and DataFrame (of Rle) and allows you to center your
  log2ratio values on 2 copies, or on a specified 'ploidy' for your sample.

DEPRECATED AND DEFUNCT

- The BAFSet and CNSet classes have moved from deprecated to defunct. These classes only added
  getter/setter methods for baf/lrr/cn. Since these only cover some possible assayDataElements, it
  is better to use x[i,j,k], where k is the name of an assayDataElement.

- All RangedData-related things have now progressed to defunct. Please use GRanges for locData and
  everywhere else. Since the genoset package provides a common API for GenoSet, GRanges, and
  RangedData, I hope this will be a simple change for everyone.

- sampleNames and featureNames are now *un-deprecated*.  Feel free to use them. They just call
  colnames and rownames, respectively. I have defined my own rownames and colnames functions. The
  eSet ones seem to do a lot of extra work, and the getter versions read from the first
  assayDataElement, rather than pData and fData. I changed to the latter, so BigMatrix
  assayDataElements will remain untouched until you really mean to access them.

GGBase
------

Changes in version 3.25.1:

- plot_EvG will, for probeId() first argument, attempt to translate to gene symbol for y axis
  labeling

GGtools
-------

Changes in version 4.11:

- appraise and calfig are provided to foster evaluation of eQTL-based GWAS hit predictions

- cisAssoc has been added to obtain assay data from SummarizedExperiment and variants from VCF

- CisConfig instances keepMapCache now defaults to TRUE

- for naming symmetry we now have cisScores and transScores operating on CisConfig and TransConfig
  instances respectively

- Added transeqByCluster() for nested concurrency for trans searches

- Added TransConfig class to control trans searches, and modified transScores accordingly

- Added hmm878, chromatin map of GM12878 for assessment of chromatin state enrichments

- Added gffprocess(), of use when All.cis is used to generate chunk-specific gff3: gffprocess uses
  external sort/grep/tabix to unify the chunks into a single tabix-indexed gff3

- pifdr() has been changed to avoid approximation on a grid and to compute binning of permuted
  scores more rapidly using hist(); old behavior recoverable with legacy=TRUE

- ciseqByCluster() uses nested concurrency to perform cis searches

gmapR
-----

Changes in version 1.6.0:

NEW FEATURES

- Add median distance from nearest end (MDFNE) statistics to output of variantSummary.

- Updated GSNAP, which is orders of magnitude faster than the previous version, brings many fixes
  and offers many new features. One new feature is the clip_overlap argument, which clips
  overlapping ends of read pairs (important for variant calling).

- Updated bam_tally, which is faster and includes support for counting in soft-clipped regions.

USER-VISIBLE CHANGES

- Changes to tallyVariant statistics: drop the unique read position counts; renamed
  count.pos/count.neg to count.plus/count.minus (way better names)

- tallyVariants does a better job of carrying over the Seqinfo from the BAM file.

GOSemSim
--------

Changes in version 1.21.3:

- fixed minor bug in combineMethods <2013-12-16, Mon>

graphite
--------

Changes in version 1.9.4 (2014-04-04):

- Updated Biocarta, KEGG, NCI and Reactome data.

- New databases: HumanCyc and Panther.

Gviz
----

Changes in version 1.8.0:

NEW FEATURES

- The new HighlightTrack class to add a comon highlighting region for multiple tracks.

- The new OverlayTrack class to merge the panels of multiple tracks into a single panel.

- The reverseStrand display parameter lets you plot the data relative to the negative strand.

- The just.label display parameter adds control to the placement of group labls in AnnotationTrack
  and GeneRegionTrack objects.

- The box.legend display parameter to add a box around the legend in a DataTrack object.

- extend.right and extend.left now also take relative expansion factors (as values between -1 and
  1).

- A new shape type fixedArrow for AnnotationTrack and GeneRegionTrack objects, and the
  arrowHeadWidth and arrowHeadMaxWidth parameters to better control the arrow shapes.

- Display parameter schemes to persistently modify parameter settings.

- The new featureAnnotation and groupAnnotation parameters to better control the feature and group
  labels in AnnotationTracks.

- The new exonAnnotation and transcriptAnnotation parameters to better control the exon and
  transcript labels in GeneRegionTracks.

- The new AlignmentsTrack class to visualized aligned NGS reads in a BAM file.

BUG FIXES

- A number of significant fixes.

SIGNIFICANT USER-VISIBLE CHANGES

- Some display parameter names have been reworked, but the old ones should still work as aliases.

- Overplotting in AnnotationTrack and GeneRegionTrack objects has been minimized to be able to make
  better use of alpha blending. Also the way composite exons (e.g. part UTR, part CDS) are plotted
  has been changed. Those will now be merged into one feature as long as the exon identifier is
  identical and if they can be reduced into a single range wit min.gapwidth=1.

gwascat
-------

Changes in version 1.7.6:

USER VISIBLE CHANGES

- bindcadd_snv has been added to bind CADD scores of Kircher et al. 2014 for SNVs to any GRanges
  query

Changes in version 1.7.1:

USER VISIBLE CHANGES

- gwrngs is now serialized and explicitly made available in global workspace on attachment

- gwcat data.frame is not generated on attachment

CHANGES PRIOR TO VERSION 1.7.1

- May 30 2012 -- improve makeCurrentGwascat to handle chrX and p.Value suitably

GWASTools
---------

Changes in version 1.9.9:

- added block size support for GDS files stored with scan,snp dimensions

- gdsSubset and gdsSubsetCheck now operate on the fastest dimension of the GDS file

Changes in version 1.9.8:

- updates/bug fixes to gdsSubset/gdsSubsetCheck - different missing value attributes may be set if
  sub.storage type is different.

Changes in version 1.9.7:

- Added gdsSubset and gdsSubsetCheck functions to make a subset GDS file that includes only
  specified SNPs and scans from an existing GDS file

- Updated gdsImputedDosage and gdsCheckImputedDosage to account for IMPUTE2 gprobs files that have
  missing values (specified by three equal probability strings)

- Updated gdsCheckImputedDosage to produce optional logfile reporting any missing genotypes

Changes in version 1.9.6:

- Revised anomFilterBAF - fewer centromere spanning anomalies that aren't real, corrects some
  merging issues (previously it would merge sections that really were different split widths).
  Users should be aware that this will increase running time.

Changes in version 1.9.5:

- Remove defunct functions.

- Improve efficiency of gwasExactHW, mendelErr, assocTestRegression (reduce number of calls to
  rbind).

Changes in version 1.9.4:

- Bug fix in getChromosome method for SnpAnnotationDataFrame (proper behavior of unnamed "index"
  argument).

Changes in version 1.9.3:

- Added gdsImputedDosage function.

- GdsGenotypeReader can return transposed genotypes.

Changes in version 1.9.2:

- ScanAnnotationDataFrame and ScanAnnotationSQLite allow non-integer scanID.

- Fix getAlleleA and getAlleleB in GdsGenotypeReader to work with indels.

Changes in version 1.9.1:

- Documentation now located in vignettes/ folder.

- Added ibdAssignRelationshipsKing.

- Added support for genotype GDS files with scan x snp dimensions in GdsGenotypeReader.

HiTC
----

Changes in version 1.7.11:

NEW FEATURES

- Update NAMESPACE for BioC 2.14

BUG FIXES

- Fix a bug in the getRestrictionSitesPerChromosome and matchPattern function (fixed=TRUE). If
  FALSE, all 'N' are reported.

Changes in version 1.7.5:

NEW FEATURES

- Update for BioC 2.14

SIGNIFICANT USER-VISIBLE CHANGES

- Add fit.out parameter in the plotIntraDist to remove the outliers during the regression

Changes in version 1.7.4:

SIGNIFICANT USER-VISIBLE CHANGES

- New parameter for normICE - spars.filter, to filter out the more sparse bins before normlization

- Change the parameters of getAnnotatedRestrictionSites and setGenomicFeatures functions

- Update of the setGenomicAnnotation function to fit with the original HiCNorm method

BUG FIXES

- isSymmetrix - NA values

- Bug in definition of upstream and downstream flanking region for a restriction site in
  getAnnotatedRestrictionSites

Changes in version 1.7.3:

NEW FEATURES

- Add new normalization method for Hi-C data from Hu et al (HiCNorm). This method is based on linear
  regression model between interaction counts and sources of bias such as GC content, mappability,
  fragment length, etc. See normLGF(), setGenomicFeatures(), getAnnotatedRestrictionSites().

- Add new normalization method for Hi-C data from Imakaev et al.(ICE). The ICE procedure is an
  iterative normalization method used to remove any bias from HiC data.

- Add 'summary' method for HTCexp and HTClist objects

SIGNIFICANT USER-VISIBLE CHANGES

- Improve quality control methods based using sparse data

- Change method option for normPerTrans methods. The 'mean' method is in fact a 'max' method

- Update of the export.my5C list format

BUG FIXES

- max (na.rm=TRUE) in mapC function

HTSeqGenie
----------

Changes in version 3.13.13:

- removed the "no job lefts" tests in test.sclapply()

Changes in version 3.13.12:

- now uses TP53Which() (instead of which.rds) for variant calling tests

- now builds in BioC

Changes in version 3.13.11:

- now use TxDb.Hsapiens.UCSC.hg19.knownGene's version number in the TP53 genomic feature cache file,
  to avoid bugs when updating TxDb.Hsapiens.UCSC.hg19.knownGene

Changes in version 3.13.10:

- added an optional which config to limit VariantTools variant calling to certain regions

Changes in version 3.13.9:

- NA

Changes in version 3.13.8:

- fixed the "as.vector(x, mode) : invalid 'mode' argument" bug by adding the
  "importFrom(BiocGenerics,table)" NAMESPACE directive

Changes in version 3.13.7:

- enable indelRealigner in runPipeline controlled by config value of alignReads.do

Changes in version 3.13.6:

- set gatk.path in options if GATK_PATH env is set to an existing file. This is mostly for allowing
  the unit tests to know if and where a GATK is installed

Changes in version 3.13.5:

- require bamtally bugfixed version of gmapR

Changes in version 3.13.4:

- fixed bug in generateCountFeaturesPlots() when no genes left after filtering.  Practically only
  happens for bacterial or viral genomes when every gene is hit more than 200 times

Changes in version 3.13.3:

- now saving fragmentLength in summary_coverage.txt

- using Jinfeng's weighted mean to estimate fragmentLength

- removed the parallelization in saveCoverage (since saveCoverage is used in parallel)

Changes in version 3.13.2:

- activate read trimming for illumina 1.5 quality using the 'B' as qual indicator

Changes in version 3.13.1:

- deactived the unit tests (test.wrap.callVariants, test.wrap.callVariants.parallel,
  test.wrap.callVariants.rmsk_dbsnp) until VariantAnnotation is fixed, to be able to compile
  HTSeqGenie on BioC servers

Changes in version 3.13.0:

- change to new dev cycle version number

HTSFilter
---------

Changes in version 1.2.1:

- -- Minor bug fix for integeration in DESeq2 pipeline (which now has automatic independent
  filtering available) -- HTSFilter now accepts data with multi-factor designs for all input types

illuminaio
----------

Changes in version 0.6.0 (2014-04-11):

- The version number was bumped for the Bioconductor release version, which is now BioC v2.14 for R
  (>= 3.1.0).

Changes in version 0.5.6 (2014-03-10):

- Added a CITATION file, cf. citation("illuminaio").

Changes in version 0.5.5 (2014-01-19):

- Added support for encrypted IDAT files where not all data fields are of the same length e.g.
  HumanHap550 v1.

Changes in version 0.5.2 (2013-11-05):

- moved testData to IlluminaDataTestFiles.

- Added unitTests for encrypted IDATs.

Changes in version 0.5.1 (2013-11-04):

- Resshuffled internal code in readIDAT_nonenc so that seek() is always forward based.

Changes in version 0.5.0:

- The version number was bumped for the Bioconductor devel version.

- No updates.

inSilicoDb
----------

Changes in version 2.0.0:

NEW FEATURES

- Ability to access InSilico MySafe

- New methods make it now possible to check the accessibility of data

- Better integration of curated data

- Possibility to use SCAN and UPC normalizations

OTHER

- Changes in the API, check the documentation for more information

inSilicoMerging
---------------

Changes in version 1.8.0:

NEW FEATURES

- Possibility to use named lists instead of eSets, check the merge-method documentation for an
  example

BUG FIXES

- Measurements are now always sorted in the same order in phenoData and exprs

OTHER

- Previous "NONE" method only available when not specifying a batch effect removal method or by
  specifying NULL.

intansv
-------

Changes in version 1.2.1:

Notes

- Add support to SV prediction programs Lumpy and softSearch.

- Modify methodsMerge to make it more robust.

- Modify all the read functions to coordinate methodsMerge.

isobar
------

Changes in version 1.9.3:

- allow the use of a combination matrix in NoiseModel reporterTagNames

- added function getProteinInfoFromTheInternet which automatically recognizes Uniprot and Entrez
  ACs. Set now as default 'protein.info.f' in the report configuration.

- Overhauled isobar-analysis.Rnw

- use of ltxtables to allow optimal column widths (longer runtime, however)

- use [[ instead of $ for accession of lists

- better column and grouping description

- allow report generation without a proteinInfo object

- set report cache directory to 'cache' (instead of '.') by default

- 1.9.3.2:

- mascotParser2.pl: allow to skip modif-conversion with -no-modiconv

- mascotParser2.pl: set --lightXML as default

- report generation: set combn.method="versus.channel", which computes the ratios against the first
  channel

- various PDF report improvements

- 1.9.3.3:

- report tables are written into tables which are loaded with LTxtable report generation takes
  longer now

- fixes in correct.peptide.ratios

- use combined protein group for peptide-protein mapping

- use only reporter proteins for mapping

- fix in creation of protein groups from template, subset by peptide _and_ modif

Changes in version 1.9.2:

- fix issue of NA in 'n.spectra' when calculating summarized ratios

- various report improvements:

- use column of variable width to display class labels

- add attributes of quant table to summarized result table

- improved placement of tikz peptide group pictures

KEGGgraph
---------

Changes in version 1.21.1 (2013-10-17):

BUG FIXES

- mergeKEGGgraph can handle NULL objects now

- KEGGpathway2Graph eliminates duplicated edges so that parseKGML2DataFrame does not fail. (Thanks
  to Paul Shannon's bug report)

- Paul Shannon added a set of RUnit tests to KEGGgraph

limma
-----

Changes in version 3.20.0:

- New functions diffSplice(), topSplice() and plotSplice() provide functionality to analyse
  differential splicing using exon-level expression data from either microarrays or RNA-seq.

- New Pasilla case study added to User's Guide, demonstrating differential splicing analysis of
  RNA-seq data.

- new function weightedLowess() which fits a lowess curve with prior weights. Unlike previous
  implementations of lowess or loess, the weights are used in calculating which neighbouring points
  to include in each local regression as well as in the local regression itself.

- weightedLoess() now becomes the default method used by loessFit() to fit the loess curve when
  there are weights. The previous locfit and loess() methods are offered as options.

- linear model fit functions lm.series(), mrlm.series() and gls.series() no longer drop the
  dimensions of the components of the fitted object when there is just coefficient or just one gene.
  Previously this was done inconsistently in some cases but not others. Now the matrix components
  always keep dimensions.

- The functions lmFit(), eBayes() and tmixture.vector() now work even when there is just one gene
  (one row of data).

- New function subsetListOfArrays(), which is used to simplify the subsetting code for RGList,
  MAList, EList, EListRaw and MArrayLM objects.

- new function tricubeMovingAverage() for smoothing a time series.

- barcodeplot() has a new option to add enrichment worms to the plot, making use of
  tricubeMovingAverage().

- New plot() methods for RGList, MAList, EList and MArrayLM class objects.  In each case, this
  produces a similar result to plotMA(). When using plot() or plotMA() on an MArrayLM object, the
  column is now specified by the 'coef' argument instead of by 'array'.

- plotMA3by2() now works on single channel data objects as well as on MAList objects.

- New function read.idat() to read files from Illumina expression beadarrays in IDAT format.

- The ctrlpath argument of read.ilmn() now defaults to the same as path for regular probes. This
  means that only one path setting is required if the regular and control probe profiles are in the
  same directory.

- read.ilmn() now sets the same probe IDs as rownames for both the expression matrix E and the
  annotation data.frame genes, providing that the probe IDs are unique.

- beadCountWeights() can now work with either probe-wise standard errors or probe-wise standard
  deviations.

- treat() has new arguments robust and winsor.tail.p which are passed through to robust empirical
  Bayes estimation.

- topTreat() now includes ... argument which is passed to topTable().

- topTable() with confint=TRUE now produces confidence intervals based on the t-distribution instead
  of on the normal distribution. It also now accepts a numeric value for the confint argument to
  specify a confidence level other the default of 0.95.

- topTable() will now work on an MArrayLM fit object that is missing the lods component, for example
  as produced by treat().

- roast() and mroast() now permit array weights and observation weights to both be specified.

- camera(), roast() and mroast() now use getEAWP() to interpret the data object.  This means that
  they now work on any class of data object that lmFit() will.

- romer() now uses propTrueNull(method="lfdr") instead of convest(). This makes it substantially
  faster when the number of genes is large.

- genas() now uses fit$df.total from the MArrayLM object. This prevents df.total from exceeding the
  total pooled residual df for the dataset. The genas() results will change slightly for datasets
  for which df.prior was very lage.

- plotDensities() is now an S3 generic function with methods for RGList, MAList, EListRaw and EList
  objects.

- plotFB is now an S3 generic function with methods for RGList and EList data objects.

- New topic help pages 10GeneSetTests.Rd and 11RNAseq.Rd. The page 10Other.Rd is deleted. All topic
  help pages are now listed under 'See also' in the package introduction page accessible by ?limma.

- avereps() was never intended to be applied to RGList or EListRaw objects. It now gives an error
  when applied to these objects instead of returning a matrix of questionable value.

- Bug fix: fitFDistRobustly() was failing when there were missing values or zero df values and
  covariate was NULL.

- Bug fix: vennDiagram() wasn't passing extra arguments (...) to plot() when the number of sets was
  greater than 3.

- Bug fix to topTreat().  Rownames were incorrectly ordered when p<1.

- bug fix to genas(), which was not handling vector df.prior correctly when the fit object was
  generated using robust=TRUE.

- bug fix to squeezeVar(). Previously there was an error when robust=TRUE and trend=FALSE and some
  of the estimated df.prior were infinite.

- bug fix to topTable() and topTableF() when sorting by F-statistic combined with p-value or lfc
  cutoffs.

metagenomeSeq
-------------

Changes in version 1.5 (2014-04-17):

- Incorporating biom-format support with the biom2MRexperiment, MRexperiment2biom and load_biome
  function.

- Added uniqueFeatures, filterData, aggregateByTaxonomy / aggTax, plotFeature and
  calculateEffectiveSamples functions.

- Renamed MRfisher to fitPA (presence-absence fisher test).

- Added warnings for normalization

- Added fitDO (Discovery odds ratio test) and fitMeta (original metastats).

- Added match.call() info to fitZig output

- Fixed missing E-Step bounds

metaMS
------

Changes in version 0.99.8:

- fixed documentation for previous fix...

Changes in version 0.99.7:

- fixed incomplete GC settings

Changes in version 0.99.6:

- replaced paste(..., "/") with file.path() in runLC.Rnw

metaseqR
--------

Changes in version 1.0.0 (2014-04-11):

NEW FEATURES

- Function to calculate the F1-score (or harmonic mean of precision and recall)

- Mature and tested enough to go from 0.x.y version to 1.0.0

Changes in version 0.99.6 (2014-04-09):

NEW FEATURES

- Example on how to estimate statistical test weights in the vignette

BUG FIXES

- Bug when exporting flags when gene or exon filters are NULL

Changes in version 0.99.5 (2014-03-31):

NEW FEATURES

- Changed directory structure of the report to look more organized and pro

- Added the ability to save and retrieve the gene model counts list (counts for each exon) in an
  .RData file to be reused in another analysis later, as summarizing genes is one of the most
  time-consuming parts.

BUG FIXES

- Minor problems with the report

Changes in version 0.99.4 (2014-03-18):

BUG FIXES

- Removed all <<- assignments

- Changed the Depends of DESCRIPTION file to import less packages in the main environment

- Added qvalue package to Depends as there was a problem with NBPSeq

Changes in version 0.99.3 (2014-03-17):

NEW FEATURES

- Function to check and warn if a main metaseqr argument is invalid (may prevent crashes during run)

BUG FIXES

- Small bug fix in read2count function

Changes in version 0.99.2 (2014-03-14):

NEW FEATURES

- Ability to display only the top x% of significant genes to avoid excessively big reports

BUG FIXES

- More code and codestyle changes to comply with Bioconductor's guidelines

Changes in version 0.99.1 (2014-03-12):

NEW FEATURES

- Ability to export counts table (raw and normalized) when reading from BAM files or just normalized
  otherwise

- Valid examples for many functions

BUG FIXES

- More code and codestyle changes to comply with Bioconductor's guidelines

Changes in version 0.99.0 (2014-02-21):

NEW FEATURES

- Support for Pan troglodytes

- Improved the simulator to i) simulate length bias and ii) better length selection

BUG FIXES

- Code changes to comply with Bioconductor's guidelines

Changes in version 0.93.0 (2014-02-12):

NEW FEATURES

- Support for Arabidopsis thaliana

- Additional validation functions

BUG FIXES

- Bugs affecting multilevel factorial analysis (single factor, more than two conditions)

- Smaller bug fixes

Changes in version 0.92.0 (2014-01-25):

NEW FEATURES

- Functions to estimate weights for combining p-values based on simulated datasets from real data

BUG FIXES

- Smaller bug fixes

Changes in version 0.91.0 (2014-01-03):

NEW FEATURES

- Functions to create simulated datasets

- Functions to create false discovery and ROC curves

BUG FIXES

- A bug in p-value vector naming causing unordered p-values for limma and baySeq when combining
  methods

- Smaller bug fixes

Changes in version 0.9.1 (2013-11-27):

NEW FEATURES

- Permutation methods for combining p-values from multiple statistics

- More interactive and compact report

BUG FIXES

- A bug in exon filters resulting in removing genes/transcripts with only one exon

- A serious bug in exon filters causing wrong genes to be filtered out under circumstances

Changes in version 0.9.0 (2013-11-21):

NEW FEATURES

- First release

MIMOSA
------

Changes in version 0.99.2:

- Added unit tests

- Passing BiocCheck

- Remove dependency on `multicore`. Now `parallel` only.

Changes in version 0.9.12:

- MIMOSA now returns a MIMOSAResultList S3 class. A lightweight wrapper for a list, with methods
  defined to extract the `fdr` and `pData`, `getW`, and `getZ` to extract the component weights and
  posterior probabilities, respectively.

- `print` for `MIMOSAResult` and `MIMOSAResultList`

- extractors for `countsTable`, the table of counts used to fit the model. Takes an argument to
  return the proportions or the counts.

- MIMOSA checks throws a more informative warning if the data is filtered into oblivion when not
  properly paired (e.g. when the user doesn't aggregate over replicate negative controls, for
  example. )

- `volcanoPlot` implemented for MIMOSAResultList

- Cleaned up package warnings and notes.

- Cleaned up imports and depends.

Changes in version 0.9.9:

- Bug Fixes ** Filtered out empty categories on the conditioning variables. ** Support for parallel
  as well as multicore

- News ** multicore support (via the multicore package) will be phased out in favor of the parallel
  package.

Changes in version 0.8.0:

- Model fitting via the 'MIMOSA' function

- Data represented via an ExpressionSet

- Formula interface support

- Removed old code, including ICS class, and dependencies upon it.

- New documentation

Changes in version 0.7.0:

- First beta release of MIMOSA

- model fitting via MCMC through the .fitMCMC function

minfi
-----

Changes in version 1.9:

- Importing the changes from 1.8 into 1.9.

- Added the withColor argument to the getProbeType function, which allows the return of "IGrn",
  "IRed", "II", instead of only "I", "II".

- Added asList argument to getControlAddress to return result as a list.

- Moved reshape from Depends to Imports.

- Dramatic improvement in memory usage of preprocessRaw.

- Updated CITATION, the minif paper is in press.

- Fixed bug with mapToGenome(..., mergeManifest = TRUE) reported by Dale Watkins
  <dale.watkins@sahmri.com> and Allegra A. Petti <apetti@genome.wustl.edu>.

- Fixed bug with mapToGenome(rSet) with rSet being a RatioSet with the CN set to NULL reported
  byAllegra A. Petti <apetti@genome.wustl.edu>.

- Added preprocessFunnorm, a new preprocessing method.

- Improvements to the speed of getAnnotation by Martin Morgan <mtmorgan@fhcrc.org>.

MotifDb
-------

Changes in version 1.5.9:

NEW FEATURES

- JASPAR_2014 added, with 592 motifs

motifStack
----------

Changes in version 1.7.7:

NEW FEATURES

- No changes classified as 'new features' (package under active development)

- add ic.scale parameter to plotMotifLogo.

BUG FIXES

- Fix bugs for motifSignature when the distance is larger than maximal distance.

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.7.6:

NEW FEATURES

- Change the class of signature to include the color sets.

BUG FIXES

- Allow motifSignature to accept Non-Binary Trees.

Changes in version 1.7.5:

NEW FEATURES

- Change the output of motifStack for further steps.

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.7.3:

NEW FEATURES

- Change the output of readPCM from a matrix to a list of pcm objects.

BUG FIXES

- No changes classified as 'bug fixes' (package under active development)

Changes in version 1.7.2:

NEW FEATURES

- No changes classified as 'new features' (package under active development)

BUG FIXES

- fix the bugs of changing plot layout when using plotMotifLogoStackWithTree

Changes in version 1.7.1:

NEW FEATURES

- No changes classified as 'new features' (package under active development)

BUG FIXES

- fix the bugs of leave names and plot region

MSnbase
-------

Changes in version 1.11.14:

- update dependency to R >= 3.1 [2014-04-05 Sat]

Changes in version 1.11.13:

- Document [get|grep]Ecols in io vignette [2014-03-31 Mon]

- typo in readMSnSet man [2014-03-31 Mon]

- updated affiliation in vignettes [2014-03-31 Mon]

Changes in version 1.11.12:

- Fixed a bug in readMzTabData reported by Hendrik Weisser [2014-03-26 Wed]

Changes in version 1.11.11:

- precomputed msx test data now has id data [2014-03-25 Tue]

- quantitation unit tests [2014-03-25 Tue]

- quantify method now accepts label-free methods [2014-03-25 Tue]

Changes in version 1.11.10:

- import mzID [2014-03-20 Thu]

- dummyiTRAQ id extdata [2014-03-21 Fri]

- add addIdentificationData method for MSnExp and MSnSet [2014-03-21 Fri]

- using data from pRolocdata (was pRoloc) [2014-03-23 Sun]

- added removeNoId,MSnExp,MSnSet methods [2014-03-23 Sun]

- added idSummary,MSnExp,MSnSet methods [2014-03-24 Mon]

- Not generating different sample per file when reading raw data. pData(.) now has systematically 1
  row if not specifically provided by the user. A message is also reported by validity,pSet if
  row(pData(.)) > 1. [2014-03-24 Mon]

- NAnnotatedDataFrame now has default multiplex 1 [2014-03-25 Tue]

- NAnnotatedDataFrame unit tests [2014-03-25 Tue]

Changes in version 1.11.9:

- write.exprs can have fDataCol or fcol (for consistence) [2014-03-17 Mon]

- Fixing bug in combineFeatures(..., is.character(groupBy)) [2014-03-19 Wed]

- fixed combineFeatures [2014-03-20 Thu]

- added example, test and doc for combineFeatures with list [2014-03-20 Thu]

Changes in version 1.11.8:

- adding redundancy handling to combineFeatures (by vladpetyuk, pull request #18) [2014-03-14 Fri]

- updated combineFeatures signature to accomodate above changes [2014-03-14 Fri]

- updated unit tests for new testhat 0.8 [2014-03-14 Fri]

Changes in version 1.11.7:

- NA

Changes in version 1.11.6:

- add corresponding xcms functions to the chromatogram and xic manual page [2014-02-21 Fri]

- new bpca imputation methods [2014-02-27 Thu]

- replacing stop_on_error with option in vignette [2014-02-27 Thu]

Changes in version 1.11.5:

- typo in MSnSet droplevels man [2014-01-27 Mon]

- typo in MSnbase-demo vignette [2014-02-20 Thu]

- fix BPI legend in chromatogram [2014-02-20 Thu]

Changes in version 1.11.4:

- passing ... to sweep when normalising [2013-12-08 Sun]

- updated makeMTD to accomodate new MS ontology [2013-12-23 Mon]

Changes in version 1.11.3:

- updated mzTab example files to new url [2013-11-15 Fri]

- warning about mzTab versions [2013-11-15 Fri]

Changes in version 1.11.2:

- move inst/doc to vignettes [2013-10-19 Sat]

Changes in version 1.11.1:

- document na.rm in combineFeatures Rd [2013-10-18 Fri]

Changes in version 1.11.0:

- New devel version for Bioc 2.14

MSstats
-------

Changes in version 2.1.3:

- fix the groupComparison for label-free experiments.

- automatically generate progress report as .txt files

- add progress message for groupComparison and dataProcessPlots function.

Changes in version 2.1.1:

- fix the bug in Condition plot : 1. for label-based : match reference and endogenous 2. for
  label-free : when there is one observation in each group, SD=NA. make it zero.

- fix the bug in heatmap and comparison plots : remove NA result for plotting

- fix the bug for label-free groupComparison : how to get subject_nested parameter in
  make.contrast.free for unequal number per group

- fix the bug in group quantification : make.contrast.group.quantification fixed for subject_nested
  parameter

mzID
----

Changes in version 1.1.6:

- Added documentation

Changes in version 1.1.5:

- Now computes mzR compatible acquisitionNum for the scans (thanks to Sebastian Gibb)

- Checks for existence of local files (thanks to Sebastian Gibb)

- XML moved from depend to import

Changes in version 1.1.4:

- Now computes mzR compatible acquisitionNum for the scans (thanks to Sebastion Gibb)

- Various bug fixes

Changes in version 1.1.3:

- Introducing the mzIDCollection class to handle multiple mzID objects

Changes in version 1.1.2:

- Added the possibility to create the different 'sub'classes of mzID directly from a file without
  first having an internal xml representation and namespace

Changes in version 1.1.1:

- Fixed bug where multiple names in the modification rules would crash the parsing. Multiple names
  gets collapsed with '/'

mzR
---

Changes in version 1.9.8:

- Pointing to the relevant wiki page in the Rcpp compiler/linker warning [2014-04-03 Thu]

Changes in version 1.9.7:

- modify to new biocViews to DESCRIPTION file (s.arora)

Changes in version 1.9.6:

- import all of Rcpp to avoid warnings in reverse dependencies <2014-02-14 Fri>

Changes in version 1.9.3:

- fix a string in ramp.cpp to enable compilation on clang-3.4

Changes in version 1.9.2:

- version bump for Rcpp 0.10.6 <2013-10-30 Wed>

Changes in version 1.9.1:

- moved vignettes to /vignettes <2013-10-17 Thu>

nondetects
----------

Changes in version 0.99.2:

BUG FIXES

- Added informative NEWS file.

- Added version information to DESCRIPTION file.

Changes in version 0.99.1:

BUG FIXES

- Changes in coding style to make package consistent with Bioconductor.

- New biocViews added to package.

- Added importing of pData from Biobase package.

- Added unit tests.

Changes in version 0.99.0:

- Package released.

npGSEA
------

Changes in version 0.99.0:

- Package submitted to Bioconductor (March 31 2014)

pathview
--------

Changes in version 1.3.6:

- ajusted node x coordinate by +0.5 to better fit the color blocks in 2 layer native kegg views.

Changes in version 1.3.4:

- updated bods to included an extra column of id.type, the default gene ID type.

Changes in version 1.2.4:

- updated korg to included over 600 newly added species. Pathview can work with 2970 species now.

- Make returned values from pathview, keggview.native and keggview.graph functions invisible.

Changes in version 1.2.3:

- Fixed bug in node.map function, which produces 0 values when all multiple genes in a node are
  NA's.

Changes in version 1.2.2:

- Fixed bug in mol.sum function, which generates "incorrect number of dimension" or NA's when
  sum.method="median" etc.

Changes in version 1.2.1:

- Fixed bug in "missing red disease gene node labels" in diease pathways. To avoid interfering with
  node coloring, set all disease gene labels to black instead.

phyloseq
--------

Changes in version 1.7.24:

USER-VISIBLE CHANGES

- Added support for [Partial] Constrained Analysis of Principal Coordinates (CAP).

- A supported/documented option in `ordinate`, supported by `plot_ordination`.

- This solves [Issue 312](https://github.com/joey711/phyloseq/issues/312).

- The `ordinate` function now takes an explicit `formula` argument.

- This facilitates reliable contrained ordination calls for:

- CAP (this commit)

- RDA (partial redundancy analysis)

- CCA (constrained correspondence analysis)

Changes in version 1.7.23:

USER-VISIBLE CHANGES

- Refactor `plot_ordination` to be more stable, error less and give informative warnings.

- Expect no critical API changes. Some errors now informative warnings with useful auto-changes to
  parameters.

- The `type='biplot'` option no longer hard-specifies a discrete color scale. Available default
  pallette should work.

- For `type='biplot'`, the non-variable (Taxa or Sample) label will always appear first in a
  discrete legend.

Changes in version 1.7.22:

USER-VISIBLE CHANGES

- Revised psmelt to automatically modify data column-name conflicts, with warning

- Udpated `psmelt` doc to formally notify users of these potential conflicts.

- This solves Issue 307: https://github.com/joey711/phyloseq/issues/307

Changes in version 1.7.21:

USER-VISIBLE CHANGES

- Updated `import_qiime` doc to emphasize it is intended for legacy QIIME files.

- Much faster and mem-efficient import of legacy QIIME and usearch files.

- Uses data.table syntax to better manage import of large files.

- Entire HMPv35 now imports in about 1 minute, low risk of mem-swap.

- Added dependency to data.table

Changes in version 1.7.20:

USER-VISIBLE CHANGES

- No user-visible changes. All future compatibility changes.

BUG FIXES

- Unit test changes to work with upcoming R release and new testthat version.

Changes in version 1.7.19:

USER-VISIBLE CHANGES

- Documentation revisions. Faster examples, updated links.

Changes in version 1.7.18:

USER-VISIBLE CHANGES

- Fix minor bug that prohibited parallel execution of weighted UniFrac

Changes in version 1.7.17:

USER-VISIBLE CHANGES

- Added `import_usearch_uc` Added first-time support for usearch “.uc” style output table.
  Addresses Issue 286, importing from UPARSE.  https://github.com/joey711/phyloseq/issues/286
  Further feedback on performance, use-cases, should be posted there.

Changes in version 1.7.16:

USER-VISIBLE CHANGES

- Fixed minor bug affecting legend-order in `plot_network` Issue 288,
  https://github.com/joey711/phyloseq/issues/288

Changes in version 1.7.15:

USER-VISIBLE CHANGES

- `tip_glom` now uses standard R clustering tools, and takes their arguments documentation and tests
  updated to reflect the change much simpler, faster

- `merge_taxa` now uses abundance to determine the achetype by default. Previously arbitrary.

Changes in version 1.7.14:

USER-VISIBLE CHANGES

- Minor change in mixture model vignette, revised graphic

Changes in version 1.7.13:

USER-VISIBLE CHANGES

- Deprecated originalUniFrac() internal function old (original) unifrac algorithm no longer
  supported.  Addresses Issue 66: https://github.com/joey711/phyloseq/issues/66

Changes in version 1.7.12:

USER-VISIBLE CHANGES

- Formal deprecation of functions using .Deprecated Issue 269,
  https://github.com/joey711/phyloseq/issues/269

- Fixed bug in interface with vegan::fisher.alpha(..., se=TRUE).  vegan doc states that this returns
  a data.frame, but a data.frame is not returned in vegan version 1.7.10.  phyloseq no checks output
  dimensions before processing in `estimate_richness`

- Replaced deprecated functions in tests and documentation.

Changes in version 1.7.11:

USER-VISIBLE CHANGES

- Adds warning in make_network() and error in plot_network if empty graph encountered Issue 275,
  check/warning for empty igraph objects https://github.com/joey711/phyloseq/issues/275

- rarefy_even_depth() messages changed from cat() to messages(), and optional verbose argument added
  https://github.com/joey711/phyloseq/issues/263

Changes in version 1.7.10:

USER-VISIBLE CHANGES

- Fixes build-error originating from change in ade4 NAMESPACE in version 1.6.2

- Change minimum ade4 version to 1.6.2

- Uncommented examples now included in documentation for DPCoA function

Changes in version 1.7.9:

USER-VISIBLE CHANGES

- Fixed typo-derived bug in new vignette.

- These changes allow user to build from source without error.

Changes in version 1.7.8:

USER-VISIBLE CHANGES

- Package dependencies reduced/clarified: https://github.com/joey711/phyloseq/issues/259 Should
  reduce chances for collisions with other packages, and related issues.  Removed any dependencies
  on the picante package.

- Replaced picante::node.age() with a faster implementation, node_ages() Appears to be 3 times
  faster.  Speeds up UniFrac() and tip_glom() calculations.

Changes in version 1.7.7:

USER-VISIBLE CHANGES

- Tree fixes: https://github.com/joey711/phyloseq/issues/235
  https://github.com/joey711/phyloseq/issues/255

- If a tree has NA branch-length values, they are automatically set to 0.  This occurs within both
  phyloseq(), and read_tree().

- UniFrac calculations require a rooted tree. While a rooted tree is not required to be part of a
  phyloseq object, it is a helpful default behavior to select a random root when UniFrac is called
  and the tree is unrooted, flashing a notice to the user.

- Precise import from ape-package, rather than full-import.  Smaller chance for collisions.
  Precisely-defined dependencies listed in NAMESPACE

- As a result of the previous, phyloseq defines a placeholder "phylo" class, extended from "list".
  This seems to match the class from a full import of ape, and is necessary since ape does not
  export the "phylo" class.

Changes in version 1.7.6:

USER-VISIBLE CHANGES

- Merged branches 1.7.4 and 1.7.5

Changes in version 1.7.5:

NEW FEATURES

- User-specified axis ordering to plot_heatmap()

- User-specified axis edges to plot_heatmap()

- This addresses: [Issue 237](https://github.com/joey711/phyloseq/issues/237) [Issue
  230](https://github.com/joey711/phyloseq/issues/230)

USER-VISIBLE CHANGES

- New arguments to plot_heatmap(): `taxa.order`, `sample.order`, `first.sample`, `first.taxa`

Changes in version 1.7.4:

NEW FEATURES

- import_mothur now handles more formats

- Added documentation to discourage .group/.list formats

Changes in version 1.7.3:

NEW FEATURES

- Added phyloseq_to_deseq2() wrapper function and examples for computing multiple OTU tests using
  Negative Binomial model and GLM (DESeq2).

USER-VISIBLE CHANGES

- Also added new .Rmd vignette for using DESeq, with colorectal carcinoma data

Changes in version 1.7.2:

USER-VISIBLE CHANGES

- Reformat NEWS (this) file.

Changes in version 1.7.1:

USER-VISIBLE CHANGES

- Rmd/HTML-based vignettes. No more Rnw/Sweave/PDF

- Updated installer

piano
-----

Changes in version 1.4.0:

NEW FEATURES

- Added argument plot to consensusHeatmap() so that drawing the heatmap can be suppressed but the
  corresponding numerical matrix can be saved.

- Added argument cellnote to consensusHeatmap() so that the information inside each cell of the
  heatmap can be chosen to be either the consensus scores (as previously), the median p-values, the
  number of genes or empty.

- Added a matrix nGenesMat to the output of consensusHeatmap() containing the same information as
  printed in the heatmap if argument cellnote="nGenes".

- Added arguments columnnames, colorkey, colorgrad and cex to consensusHeatmap() for better control
  of the column labels, toggling of the colorkey, color selection and text size.

- Introduced the new function GSAheatmap(), which is similar to consensusHeatmap() but for only a
  single gene set result (gsaRes object).

BUG FIXES

- Fixed a bug which for some settings of runGSA would not output number of up- and downregulated
  genes in the gene sets.

DOCUMENTATION

- Updated consensusHeatmap() man page according to new changes in the function.

- Added man page for GSAheatmap().

plethy
------

Changes in version 1.1.3:

NEW FEATURES

- Added 'summaryMeasures' method.

- Added experimental heatmap-like plotting functionality 'plethy:::mvtsplot'

procoil
-------

Changes in version 1.13.2:

- cleared up package dependencies and namespace

- reference to Bioconductor Git-SVN bridge

- minor corrections to vignette

pRoloc
------

Changes in version 1.3.19:

- fixed error introduced with mclust 4.3 (that now returns the data in the Mclust output - see
  comment in pRoloc:::gmmOutliers for details) [2014-04-07 Mon]

Changes in version 1.3.18:

- getPredictions can take class-specific scores [2014-04-04 Fri]

Changes in version 1.3.17:

- fixed newly introduced bug (see 1.3.16) in pRoloc:::subsetAsDataFrame - thank you unit tests for
  saving me, again [2014-03-26 Wed]

Changes in version 1.3.16:

- pRoloc:::subsetAsDataFrame now preserved original sample/column names [2014-03-24 Mon]

- fixed wrong message when using col and pch in plot2D [2014-03-25 Tue]

Changes in version 1.3.15:

- updated pRolocmarkers("mmus") [2014-03-21 Fri]

- moved extdata/*csv to pRolocdata [2014-03-23 Sun]

- using *csv from pRolocdata [2014-03-23 Sun]

Changes in version 1.3.14:

- deleted tab character [2014-03-15 Sat]

- added support for GMM parametrisation to phenoDisco [2014-03-17 Mon]

- message instead of warning when using colour and pch [2014-03-21 Fri]

- remove 1 duplicated mouse marker [2014-03-21 Fri]

Changes in version 1.3.13:

- fixed a bug in addLegend [2014-03-14 Fri]

- updated testing to testthat 0.8 [2014-03-14 Fri]

- Fixing several warnings about symnbols being replaced upon pRoloc loading and note about usage of
  ::: [2014-03-15 Sat]

Changes in version 1.3.12:

- added phenoDisco2 for testing, allows choice of GMM parameters [2014-02-26 Wed]

- removed duplicated fly markers [2014-02-28 Fri]

- updated affiliations in vignettes [2014-03-10 Mon]

Changes in version 1.3.11:

- modify to new biocViews to DESCRIPTION file by s.arora [2014-03-04]

Changes in version 1.3.10:

- new phenoDisco ndims argument to use more than two two principal components as input for discovery
  analysis [2014-01-03 Mon]

- fixed and updated phenoDisco logging [2014-02-03 Mon]

- added support for parallel phenoDisco execution. See BPPARAM argument [2014-02-03 Mon]

- fixed issues with using PD and ndims [2014-02-10 Mon]

- fixed call to anyUnknown in PD code [2014-02-10 Mon]

- checking if duplicated markers in addMarkers [2014-02-14 Fri]

Changes in version 1.3.9:

- bump to force rebuild for new Rcpp

Changes in version 1.3.8:

- Removed trailing space in mmus nucleus markers [2014-01-21 Tue]

- using filterNA to remove features with missing values in plot2D [2014-01-21 Tue]

- fixed plot2D/addLegend [2014-01-23 Thu]

Changes in version 1.3.7:

- fixed addLegend to use correct colours (order) [2014-01-20 Mon]

- fix typo in addMarkers man [2014-01-20 Mon]

- re-arranged stockpch so that interleave full and empty plotting character [2014-01-20 Mon] 0
  Removed last stockcol (tomato), too cose to "#FF7F00" [2014-01-20 Mon]

Changes in version 1.3.6:

- updated human markers: keep new pd.markers phenotypes, validated by Lisa and remove singletons
  [2014-01-14 Tue]

- first stockpch is noe 19 [2014-01-16 Thu]

- plot2D(.., pch) now taken into account for labelled data [2014-01-16 Thu]

- removed alpha plot2D argument [2014-01-17 Fri]

- updated plot2D and addLegend function with support for more organelle groups than colours. The
  previous versions are available as plot2D_v1 and addLegend_v1. [2014-01-17 Fri]

Changes in version 1.3.5:

- pRoloc citation [2014-01-12 Sun]

Changes in version 1.3.4:

- new unknownSet function [2013-12-11 Wed]

- updated vignettes to account for tan2009r1 changes [2013-12-16 Mon]

Changes in version 1.3.3:

- update citation in phenoDisco.Rd [2013-11-22 Fri]

- typos in vignettes [2013-11-25 Mon]

Changes in version 1.3.2:

- new markers [2013-11-06 Wed]

- markers in vinette [2013-11-15 Fri]

- mouse markers [2013-11-18 Mon]

Changes in version 1.3.1:

- using combineFeatures(..., na.rm=TRUE) in lopims pipeine [2013-10-22 Tue]

- corrected spelling errors in phenoDisco doc [2013-10-29 Tue]

Changes in version 1.3.0:

- next devel version fot Bioc 2.14

PWMEnrich
---------

Changes in version 3.5:

- After further testing revert back to PWMEnrich 2.x group P-value algorithm

- Introduced group sorting by top motifs

Changes in version 3.1.4:

- New way of estimating P-value for groups of sequences. Note this will produce different P-values
  for groups of sequences than PWMEnrich 2.x !

QDNAseq
-------

Changes in version 1.0.0 (2014-04-14):

INITIAL RELEASE

- initial release as part of Bioconductor 2.14

qpgraph
-------

Changes in version 1.20:

USER VISIBLE CHANGES

- Functions qpRndHMGM() and qpSampleFromHMGM() which were defunct in the previous release, are now
  removed from the package.

- The execution of function qpGetCliques() can be now interrupted with CTRL+C and should allow to
  process GUI events.

- Function qpBoundary() with argument logscale.bdsize=TRUE now handles plotting of boundaries of
  size 0 by setting them to 1 (=log(1)=0) for plotting purposes only.

NEW FEATURES

- Function qpCItest() when assuming a mixed GMM returns the fraction of variance explained (partial
  eta-squared) as estimate. See manual page for further details.

BUG FIXES

- Several bugfixes related to the simulation of mixed GMMs and eQTL networks.

QuasR
-----

Changes in version 1.4.0:

NEW FEATURES

- new arguments mapqMin and mapqMaxm in qCount, qProfile, qMeth and qExportWig allow to select
  alignments based on mapping quality

- qAlign gains checkOnly argument which allows checking for existing alignments without triggering
  new alignments in case of missing ones

r3Cseq
------

Changes in version 1.9.2 (2013-11-11):

- fixed the bugs in function getFragmentsPerWindow for "overlapping window"

Changes in version 1.9.1 (2013-10-21):

- added the yLim input parameter of plotInteractionsNearViewpoint

- added the log2 cutoff ratio for the function plotInteractionsNearViewpoint

- fixed the shift of interaction in the viewpoint chromosome found in the function
  getFragmentsPerWindow

- fixed the missing input parameter for the function getBatchReadCountPerWindow

- fixed the class of r3CseqInBatch

RamiGO
------

Changes in version 1.9.5:

- RamiGO is now using AmiGO v2 Visualize as the default web server.

Changes in version 1.9.1:

- added citation information. Please cite
  http://bioinformatics.oxfordjournals.org/content/29/5/666.full

Rariant
-------

Changes in version 1.0.0:

- First stable release with Bioconductor 2.14 (2014-04)

rBiopaxParser
-------------

Changes in version 2.1.7:

- BUGFIX: fixed setkey for parsed biopax models

- BUGFIX: fixed possible uses of complete URLs for ids (e.g. in PathwayCommons)

- BUGFIX: fixed parameter checking failing due to data.table version 1.9.2. This will be fixed
  anyways when data.table 1.9.3. is out.

Changes in version 2.1.5:

- BUGFIX: Fixed several bugs relating to the default BioPAX version parameter biopaxlevel.

- The package requires R >= 3.0 due to data.table requiring Reshape2, requiring plyr, requiring
  Rcpp. :-(

Changes in version 2.1.3:

- BUGFIX: Switched from rep_len to rep in order to continue R<3.0 support.

Changes in version 2.1.1:

- Attention: This is a major update and might introduce backwards incompatibility to previously
  parsed data!

- Speed issues in previous versions were common with huge databases like Reactome or PathwayCommons,
  and with search & replace operations in the internal data.frame.

- The data.table package is now required and imported for rBiopaxParser.

- The internal data model has been switched from data.frame to data.table to speed up working with
  the rBiopaxParser.

- The "df" slot of the biopax object, a data.frame object, is gone. Use the "dt" slot, a data.table,
  to directly access to biopax data.

- If you work with the "df" slot of the biopax object in your current code, you have to make
  yourself comfortable with using the data.table syntax!

- All functions which return a subset of the internal data, return a data.table object now.

- A long standing bug of deeply nested XML instances has been fixed. See
  http://permalink.gmane.org/gmane.science.biology.informatics.conductor/48534

- The default parameter specifying the Biopax Level (usually called biopaxlevel) has been changed
  and defaults to 3 now.

- Functions selectInstances, getInstanceProperty, getReferencedIDs, splitComplex, getInstanceClass
  and many others now acceps either a biopax object or a compatible internal data.table.

Changes in version 1.3.3:

- fixed download links. Added links for the NCI databases with fixed XML encoding ( see discussion
  at http://sourceforge.net/mailarchive/message.php?msg_id=30106560)

- added TemplateReactionRegulation to be plotted in pathway2RegulatoryGraph

- future versions will likely switch to data.table support to deal with huge BioPAX databases from
  Reactome et al

RDAVIDWebService
----------------

Changes in version 1.1.4:

BUG FIXED

- `getFunctionalAnnotationChartFile` bug was fixed. Now threshold and count parameter works as
  supposed (Thanks to Ulrik Stervbo)

Changes in version 1.1.3:

MINOR CHANGES

- Added again into the repository.

Changes in version 1.1.2:

CODE CHANGES

- In setAs(from="data.frame", to="DAVIDFunctionalAnnotationChart") function PValue field force to
  numeric data type (thanks to Marc Carlson)

NEW FEATURES

- Apache Axis time out parameter modification set/getTimeOut

- Apache Axis http protocol version parameter modification set/getHttpProtocolVersion

DOCUMENTATION

- Vignette Trouble shooting section added with apache axis parameter examples.

Changes in version 1.1.1:

DOCUMENTATION

- CITATION file was included.

- Minor updates to DAVIDWebService-class and DAVIDWebService-package man pages and Roxygen2 quotes,
  in order to include the article citation.

- Vignette email parameter invocation corrected.

CODE CHANGES

- DAVIDWebService initialize email parameter order modification (email, ..., url)

DEPENDENCY CHANGES

- Category, GO.db, RBGL and rJava packages were moved to import.

- NAMESPACE file updated accordingly.

Changes in version 1.1.0:

MINOR CHANGES

- Bump y in version x.y.z to odd number in 2.14 devel

ReactomePA
----------

Changes in version 1.7.2:

- bug fixed for multi-organisms support <2013-12-09, Mon>

RedeR
-----

Changes in version 1.12.0:

- Improved interactive layout algorithm.

- User interface is simplified, some shortcuts have been added.

- Final release cycle of deprecated functions moved to the addGraph/getGraph related methods.

- Deprecated plugin builder methods have been removed.

RefNet
------

Changes in version 1.0.0:

NEW FEATURES

- Two sets of interactions, retrieved automatically from the AnnotationHub, in addition to the many
  provided by PSICQUIC: - gerstein-2012: "Architecture of the human regulatory network derived from
  ENCODE data", Gerstein et al, pmid 22955619 - hypoxiaSignaling-2006: "Hypoxia signalling in cancer
  and approaches to enforce tumour regression", Pouyssgur et al, pmid 16724055

- These interactions have been expressed in the standard RefNet 19-column format for easy querying.

- a.canonical

- b.canonical

- relation

- bidirectional

- detectionMethod

- pmid

- a.organism

- b.organism

- a.common

- a.canonicalIdType

- b.common

- b.canonicalIdType

- cellType

- a.modification

- a.cellularComponent

- b.modification

- b.cellularComponent

- provider

- comment


ReQON
-----

Changes in version 1.9.1:

- Minor changes were made to the vignette. September 4, 2012 Our manuscript describing ReQON has
  been published in BMC Bioinformatics. Refer to this paper for further description of the ReQON
  algorithm and a comparison with other quality score recalibration algorithms.  Please cite this
  paper if you use our package. Cabanski CR et al. (2012) ReQON: a Bioconductor package for
  recalibrating quality scores from next-generation sequencing data. BMC Bioinformatics 13(221).
  doi:10.1186/1471-2105-13-221.

Rgraphviz
---------

Changes in version 2.7:

- Fixed include statements in Graphviz code.

rhdf5
-----

Changes in version 2.8.0:

NEW FEATURES

- New function h5version implemented.

- New low level general library functions H5open, H5close, H5garbage_collect, H5get_libversion, and
  H5Dset_extent implemented.

USER VISIBLE CHANGES

- h5createDataset automatically uses chunking and compression.

- Added a warning if chunk size is equal to dimensions for large compressed datasets.

BUG FIXES

- C-stack overflow when reading large fixed-length strings.

- error in i/o with chunksize or blocksize parameters.

- compiling errors due to missing int return value.

roar
----

Changes in version 0.99.6:

NEW FEATURES

- No new features, just adopted the suggested format for the NEWS file

BUG FIXES

- No changes classified as 'bug fixes'

rols
----

Changes in version 1.5.2:

- pretty printing with strwrap(mtd[...]) in vignette and termMetadata S3 printing method for
  readable output (suggestion from Martin Morgan) [2014-02-18 Tue]

- Using BiocStyle vignette [2014-02-18 Tue]

Changes in version 1.5.1:

- split go terms (PRO:...) instead of ids (PR:...) in vignette tgnqueryShow chunk <2013-10-19 Sat>

Changes in version 1.5.0:

- new devel for Bioc 2.14

rpx
---

Changes in version 0.99.9:

- In pxref, dealing with cases without any or pending publications [2014-04-01 Tue]

Changes in version 0.99.8:

- fix leftover issues from renaming package (d.tenenbaum) [2014-03-30 Sun]

Changes in version 0.99.7:

- renamed to rpx due to CRAN name clash [2014-03-28]

Changes in version 0.99.6:

- more typos in vignette [2014-03-12 Wed]

Changes in version 0.99.5:

- Suggestions, biocView and typos reported by Nate Hayden [2014-03-12 Wed]

- added Runit test [2014-03-12 Wed]

- use ae fonts in vignette [2014-03-12 Wed]

Changes in version 0.99.4:

- Update reference in man and vignette [2014-03-11 Tue]

- Typos in vignette [2014-03-11 Tue]

Changes in version 0.99.3:

- reverting pxannounced to empty signature to fix warning [2014-03-05 Wed]

Changes in version 0.99.2:

- replacing download.file by getURL to download annoucement rss over https <2014-03-05 Wed>

Changes in version 0.99.1:

- fixing pxfiles ending by \r on Windows [2014-03-05 Wed]

Changes in version 0.99.0:

- added a vignette [2014-03-05 Wed]

- handle invalid identifiers [2014-03-05 Wed]

rqubic
------

Changes in version 1.8.1 (2011-08-02):

- Fix errors in documentations

Rsamtools
---------

Changes in version 1.15.0:

NEW FEATURES

- asSam converts BAM files to SAM files

- razip, bgzip re-compress directly from .gz files

- yieldReduce through a BAM or other file, applying a MAP function to each chunk and reducing the
  result to it's final representation

SIGNIFICANT USER-VISIBLE CHANGES

- bgzip default extension changed to '.bgz'

- seqinfo,BamFile-method attempts to return seqnames in 'natural' order, e.g., chr1, chr2, ...

- yieldSize now works on BAM files queried with ranges. Successive ranges are input until the total
  number of records first equals or exceeds yieldSize..

- scanFa supports DNA, RNA, and AAStringSet return objects

BUG FIXES

- scanFa returns correct sequence at the very end of files

- razip compresses small files

- applyPileups no longer crashes in the absence of an index file

Rsubread
--------

Changes in version 1.14.0:

NEW FEATURES

- featureCounts automatically re-orders paired end reads if they were sorted by chromosomal
  locations in the input. It can also deal with read pairs in which only one end was included in the
  input.

- featureCounts is more robust in processing different variants of GTF/GFF annotation files.

- Subjunc can detect exon-exon junctions which are located at the start or end of reads. It can
  detect non-canonical junctions ('reportAllJunctions' option) in addition to canonical junctions,
  and it can also be used to detect chimerism in both RNA-seq and gDNA-seq data.

- Both align (Subread aligner) and subjunc now take gzipped FASTQ files and outputs BAM files in
  their default settings. They accept multiple input files as well.

- Breakpoint locations are reported along with mapping location of each fusion read in SAM/BAM
  files, using tags including CC(chromosome name), CP(mapping position), CG(CIGAR string) and
  CT(strand).

- A full index (no gaps) can be built for a reference genome to further speed up read mapping.

- qualityScores() and propmapped() functions were rewritten.

- Bug fixes.

rTANDEM
-------

Changes in version 1.3.9:

- Capitalisation aliases added for most functions.

- The rtandem function now has an 'output.path' argument.  BUG FIXES

- Changed default parameters to prevent the risk of not generating output files.

- Fixed a bug with 'rtandem()' where the function would not take a complex rTTAXO object as
  argument.

Changes in version 1.3.8:

BUG FIXES

- Fixed a memory leak associated with PTMTreeSearch.

- Removed an unwanted side-effect from the function GetPeptides().

Changes in version 1.3.7:

NEW FEATURES

- rTANDEM now includes the PTMTreeSearch scoring algorithm.

- New function to set default values for PTMTreeSearch algorithm.

Changes in version 1.3.5:

NEW FEATURES

- New rTResult_s class support saving spectra.

- New function 'ms2.plot' to visualize spectra.

Changes in version 1.3.4:

NEW FEATURES

- rTANDEM now supports TPP's hrk-score function.

BUG FIXES

- Fixed a bug with k-score function.

Changes in version 1.3.3:

NEW FEATURES

- Functions have been implemented to create and manipulate parameter objects (S3 class rTParam).
  Those functions give default values for instrument-specific parameters for some types of mass
  spectrometer as well as default values for general parameters.

Changes in version 1.3.2:

NEW FEATURES

- rTANDEM now supports the k-score function. It can be used by setting the parameter 'scoring,
  algorithm' to "k-score".

RTN
---

Changes in version 1.2.0:

- Introduced variant set enrichment analysis for regulons.

rTRM
----

Changes in version 1.2:

- An arbitrary number of targets can be specified in findTRM(). This allows to use the same list of
  TFs for target and query, in order to identify TRMs using information from predicted TFBS only.

- BioGRID dataset updated to release 3.2.11 (April 2014)

- getBiogridData() retrieves by default the latest BioGRID release (www.thebiogrid.org)

SANTA
-----

Changes in version 1.3.10 (2014-03-03):

NOTES: 

- Fixed minor memory leak in the Knet function.

sapFinder
---------

Changes in version 0.99.4:

BUG FIXES

- *To fix a bug of "index.js" that may cause a confusing CSS layout of index html page.

Changes in version 0.99.3:

NEW FEATURES

- *Add some BiocViews items to packages. *Reduce the size of testing data. *Make import and export
  much more continent. *To make package consistent with bioconductor in the coding style. *Use
  BiocStyle for vignette.

BUG FIXES

- *To fix a bug of vignette document that may cause the package checking to crash.

Changes in version 0.99.2:

BUG FIXES

- *Resolve some errors and warnings.

Changes in version 0.99.0:

- *Package released.

SCAN.UPC
--------

Changes in version 2.5.8:

NEW FEATURES

- This package now provides "generic" functions for UPC normalizing any type of gene-expression
  data. These values can be input as vectors or ExpressionSet objects. This means that data from
  Illumina BeadChip microarrays or any other type of microarray can be UPC normalized with relative
  ease. It is also easy through integration with the GEOquery package to UPC normalize any
  preprocessed data set from GEO.

- This package now uses the ComBat function from the sva package to make it easy to adjust gene
  expression data for batch effects. This can be done for any type of expression data.

- RNA-Sequencing data can now be input as matrix files and/or where a header is present in the input
  data files.

NOTES

- The SCAN_TwoColor and UPC_TwoColor functions have changed. They now return ExpressionSet objects.

- The SCAN_TwoColor and UPC_TwoColor methods now add an integer suffix to each probe for which there
  is a duplicate name.

- The UPC_RNASeq function now also returns an ExpressionSet object.

- A bug was fixed in which UPC_RNASeq returned NA values if there was no entry in the annotation
  file.

SeqArray
--------

Changes in version 1.3.2:

- update test codes to avoid the conflict

Changes in version 1.3.1:

- update according to the new version of VariantAnnotation

shinyTANDEM
-----------

Changes in version 1.1.4:

BUG FIXES

- Better error handling, error messages and conditions for missing data.

Changes in version 1.1.3:

BUG FIXES

- Corrected problems arising from incompatibility between data.table and data.frame subsetting
  syntaxes.

Changes in version 1.1.2:

NEW FEATURES

- MS2 spectra can now be visualized.

- The 'Result Statistics' and 'Peptides view' tab are now ready.

BUG FIXES

- Corrected a bug in the selection of proteins by number of peptides.

ShortRead
---------

Changes in version 1.21:

NEW FEATURES

- writeFastq can write (and does so by default) gz-compressed files

SIGNIFICANT USER-VISIBLE CHANGES

- Use BiocParallel rather than srapply, mark srapply as 'Deprecated'

- qa,character-method defaults to type="fastq"

- Input of 'legacy' formats marked as such

- alphabetByCycle supports amino acid string sets

snm
---

Changes in version 1.11:

- Now uses lme4 R (>= 1.0)

- Default nbins=20

- New argument lmer.max.iter allows one to cap the number of lmer iterations

SomaticSignatures
-----------------

Changes in version 1.0.0:

- First stable release with Bioconductor 2.14 (2014-04)

supraHex
--------

Changes in version 1.1.17:

NEW FEATURES

- Add new functions for advanced heatmap visualisation and tree-based analysis of sample
  relationships

Changes in version 1.1.2:

NEW FEATURES

- Add a new function "sMapOverlay" to overlay additional data to the trained map

synapter
--------

Changes in version 1.5.7:

- Enable setting number of missed cleavages (closes #53) [2014-03-24 Mon]

Changes in version 1.5.6:

- change the default value of grid.ppm.from to 2 (was 5 before) [2014-03-05 Wed]

- change the separator in GridDetails' names to ":" (was "." before) [2014-03-07 Fri]

- test for corresponding Pep3D file (closes #42) [2014-03-19 Wed]

Changes in version 1.5.5:

- fix a bug in the calculation of non unique matches in gridSearch2 (part of searchGrid); results in
  a higher number of non unique matches (some of the reported -2 will now be reported as 2 in the
  grid details) [2014-03-05 Wed]

- partial rewrite of gridSearch2 (part of searchGrid) for faster grid calculation [2014-03-04 Tue]

Changes in version 1.5.4:

- biocViews update

Changes in version 1.5.3:

- modified findMSeEMRTs (part of searchGrid) for faster grid calculation [2014-02-28 Fri]

- typos in manual [2014-02-25 Tue]

- replace readFasta by Biostrings::readAAStringSet [2014-02-25 Tue]

- replace digest by cleaver::cleave [2014-02-25 Tue]

Changes in version 1.5.1:

- typo in vignette [2014-02-18 Tue]

Changes in version 1.5.0:

- new devel version for Bioc 2.14

TargetSearch
------------

Changes in version 1.20.0:

BUG FIXES

- Write the actual retention index value instead of the work 'RI' in `ProfileCleanUp`.

TCC
---

Changes in version 1.3.2:

- DESeq2 was implemented in TCC for identifying DEGs.

- EBSeq was removed from TCC.

- arguments of 'WAD' function were changed.

- fixed bug in '.testByDeseq' that missing to treat size factors.

- the strategies for DE analysis of paired two-group dataset were implemented.

- add the section for describing DE analysis of paired two-group dataset into vignette.

tigre
-----

Changes in version 1.16.1:

- More informative error messages

TransView
---------

Changes in version 1.7.4:

BUG FIXES

- Removed import error introduced with version 1.7.2

Changes in version 1.7.3:

BUG FIXES

- Minor bug fix

Changes in version 1.7.2:

NEW FEATURES

- annotatePeaks now assigns a gene to the peak center instead of the peak limits

- annotatePeaks takes the new argument 'reference' to enable alternative peak to gene body
  associations

Changes in version 1.7.1:

BUG FIXES

- Minor correction to avoid warning messages issued after GRanges to data.frame conversion

unifiedWMWqPCR
--------------

Changes in version 0.9.9:

SIGNIFICANT USER-VISIBLE CHANGES

- First submission of uWMWqPCR to BioConductor.

VariantAnnotation
-----------------

Changes in version 1.10.0:

NEW FEATURES

- add support for ##contig in VCF header

- add 'meta<-', 'info<-', 'geno<-' replacement methods for VCFHeader

- add 'header<-' replacement method for VCF

- add strand to output from locationVariants()

- add support for writeVcf() to process Rle data in geno matrix

- readVcf() now parses 'geno' fields with Number=G as ((#alleles + 1) * (#alleles + 2)) / 2

- writeVcf() now sorts the VCF when 'index=TRUE'

- add 'fixed<-,VCFHeader,DataFrameList' method

- add convenience functions for reading VCF into VRanges

- add Rplinkseq test script

- add 'isSNV', 'isInsertion', 'isDeletion', 'isIndel', 'isTransition', 'isPrecise', 'isSV' and
  'isSubstitution' generics

- add 'isSNV', 'isInsertion', 'isDeletion', 'isIndel' methods for VRanges and VCF classes

- add match methods between ExpandedVCF and VRanges

- add support for VRanges %in% TabixFile

MODIFICATIONS

- expand,VCF-method ignores 'AD' header of 'AD' geno is NULL

- add support for SIFT.Hsapiens.dbSNP137

- remove locateVariants() dependence on chr7-sub.vcf.gz

- modify expand() to handle 'AD' field where 'Number' is integer

- rename readVRangesFromVCF() to readVcfAsVRanges()

- remove check for circular chromosomes in locateVariants() and predictCoding() and
  refLocsToLocalLocs()

- modify filterVcf() to handle ranges in ScanVcfParam

- pass 'genetic.code' through predictCoding()

- change default to 'row.names=TRUE' for readGT(), readGeno(), and readInfo()

- fixed() on empty VCF now returns DataFrame with column names and data types vs an empty DataFrame

- update biocViews

- modify 'show,VCF' to represent empty values in XStringSet with '.'

- replace rtracklayer:::pasteCollapse with unstrsplit()

DEPRECATED and DEFUNCT

- remove defunct dbSNPFilter(), regionfilter() and MatrixToSnpMatrix()

- deprecate defunct readVcfLongForm()

BUG FIXES

- modify expand.geno() to handle case where header and geno don't match

- modify writeVcf() to write out rownames with ":" character instead of treating as missing

- fix how sample names were passed from 'ScanVcfParam' to scanVcf()

- fix bug in 'show,VCF' method

- fix bugs in VRanges -> VCF coercion methods

- fix bug in lightweight read* functions that were ignoring samples in ScanVcfParam

- fix bug in writeVcf() when no 'ALT' is present

VariantFiltering
----------------

Changes in version 0.99:

USER VISIBLE CHANGES

- Submission of the first version to the Bioconductor project (start date: 26 September, 2013)

VariantTools
------------

Changes in version 1.6.0:

NEW FEATURES

- Add pileupVariants function as an alternative to tallyVariants for computing nucleotide pileups
  using Rsamtools, instead of gmapR. In the future, this should allow VariantTools to become
  independent of gmapR, but the full variant statistics will only be available via tallyVariants.

USER-VISIBLE CHANGES

- See the NEWS for gmapR to learn about changes to tallyVariants output.

wateRmelon
----------

1.3.1: thanks to several users for bug reports!

xcms
----

Changes in version 1.39.6:

USER VISIBLE CHANGES

- Massifquant reports the maximum intensity for each isotope trace (peak). This is useful for
  interactive parameter optimization.

BUG FIXES

- Major memory reduction in parallel fillPeaks() thanks to Jan Stanstrup. Now using an environment
  to mirror gvals to each list item in the very large argList.

Changes in version 1.39.4:

BUG FIXES

- Fixed write.cdf(), which had an intensity offset of +1, added a unit test

Changes in version 1.39.3:

BUG FIXES

- New R-devel check unload better. Lingering ramp code removed, import from mzR. Cleaned up other
  errors in package check.

Changes in version 1.39.1:

BUG FIXES

- Updated doubleMatrix c code to allow for larger profile matrixes

REQUIRED CHANGES

- Moved inst/doc to vignettes


yaqcaffy
--------

Changes in version 1.23.1:

- moved inst/doc to vignettes [2013-10-19 Sat]



Packages removed from the release
=================================

The following packages are no longer in the release:

Agi4x44PreProcess, maDB, pgUtils
