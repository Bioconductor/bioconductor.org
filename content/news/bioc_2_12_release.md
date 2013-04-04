April 4, 2013

Bioconductors:

We are pleased to announce Bioconductor 2.12, consisting of 671
software packages and more than 675 up-to-date annotation packages.
There are 65 new software packages, and many updates and improvements
to existing packages; Bioconductor 2.12 is compatible with R 3.0,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release includes an updated Bioconductor
[Amazon Machine Image](http://bioconductor.org/help/bioconductor-cloud-ami/).
Visit [http://bioconductor.org](http://bioconductor.org)
for details and downloads.

Contents
--------

* Getting Started with Bioconductor 2.12
* New Software Packages
* NEWS from new and existing packages
* Packages removed from the release

Getting Started with Bioconductor 2.12
======================================

To install Bioconductor 2.12:

1. Install R 3.0.  Bioconductor 2.12 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 65 new packages in this release of Bioconductor.

- [AnnotationHub](http://bioconductor.org/packages/release/bioc/html/AnnotationHub.html):
A client for retrieving data from the Bioconductor AnnotationHub
online services.

- [antiProfiles](http://bioconductor.org/packages/release/bioc/html/antiProfiles.html):
Implements gene expression anti-profiles as described in Corrada Bravo
et al., BMC Bioinformatics 2012, 13:272 doi:10.1186/1471-2105-13-272.

- [ARRmNormalization](http://bioconductor.org/packages/release/bioc/html/ARRmNormalization.html):
Perform the Adaptive Robust Regression method (ARRm) for the
normalization of methylation data from the Illumina Infinium
HumanMethylation 450k assay.

- [BaseSpaceR](http://bioconductor.org/packages/release/bioc/html/BaseSpaceR.html):
A rich R interface to Illumina's BaseSpace cloud computing
environment, enabling the fast development of data analysis and
visualisation tools.

- [biomvRCNS](http://bioconductor.org/packages/release/bioc/html/biomvRCNS.html):
In this package, a Hidden Semi Markov Model and one homogeneous
segmentation model are designed and implemented for segmentation
genomic data, with the aim of assisting in transcripts detection using
high throughput technology like RNA-seq or tiling array, and copy
number analysis using aCGH or sequencing.

- [BiSeq](http://bioconductor.org/packages/release/bioc/html/BiSeq.html):
The BiSeq package provides useful classes and functions to handle and
analyze targeted bisulfite sequencing (BS) data such as reduced-
representation bisulfite sequencing (RRBS) data. In particular, it
implements an algorithm to detect differentially methylated regions
(DMRs).  The package takes already aligned BS data from one or
multiple samples.

- [bumphunter](http://bioconductor.org/packages/release/bioc/html/bumphunter.html):
Tools for finding bumps in genomic data

- [CAGEr](http://bioconductor.org/packages/release/bioc/html/CAGEr.html):
Preprocessing of CAGE sequencing data, identification and
normalization of transcription start sites and downstream analysis of
transcription start sites clusters (promoters).

- [casper](http://bioconductor.org/packages/release/bioc/html/casper.html):
Infer alternative splicing from paired-end RNA-seq data. The model is
based on counting paths across exons, rather than pairwise exon
connections, and estimates the fragment size and start distributions
non- parametrically, which improves estimation precision.

- [chimera](http://bioconductor.org/packages/release/bioc/html/chimera.html):
This package facilitates the characterisation of fusion products
events. It allows to import fusion data results from the following
fusion finders: bellerophontes, deFuse, FusionFinder, FusionHunter,
mapSplice, tophat-fusion, FusionMap

- [cisPath](http://bioconductor.org/packages/release/bioc/html/cisPath.html):
cisPath is an R package for identification and visualization of the
shortest functional paths between proteins in the protein-protein
interaction network.

- [clipper](http://bioconductor.org/packages/release/bioc/html/clipper.html):
clipper is a package for topological gene set analysis. It implements
a two-step empirical approach based on the exploitation of graph
decomposition into a junction tree to reconstruct the most relevant
signal path. In the first step clipper selects significant pathways
according to statistical tests on the means and the concentration
matrices of the graphs derived from pathway topologies. Then, it
"clips" the whole pathway identifying the signal paths having the
greatest association with a specific phenotype.

- [CNORfeeder](http://bioconductor.org/packages/release/bioc/html/CNORfeeder.html):
This package integrates literature-constrained and data-driven methods
to infer signalling networks from perturbation experiments. It permits
to extends a given network with links derived from the data via
various inference methods, and uses information on physical
interactions of proteins to guide and validate the integration of
links.

- [copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html):
Penalized least squares regression is applied to fit piecewise
constant curves to copy number data to locate genomic regions of
constant copy number. Procedures are available for individual
segmentation of each sample, joint segmentation of several samples and
joint segmentation of the two data tracks from SNP-arrays. Several
plotting functions are available for visualization of the data and the
segmentation results.

- [DASiR](http://bioconductor.org/packages/release/bioc/html/DASiR.html):
R package for programmatic retrieval of information from DAS servers

- [deltaGseg](http://bioconductor.org/packages/release/bioc/html/deltaGseg.html):
Identifying distinct subpopulations through multiscale time series
analysis

- [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html):
Estimate variance-mean dependence in count data from high-throughput
sequencing assays and test for differential expression based on a
model using the negative binomial distribution

- [dexus](http://bioconductor.org/packages/release/bioc/html/dexus.html):
DEXUS identifies differentially expressed genes in RNA-Seq data under
all possible study designs such as studies without replicates, without
sample groups, and with unknown conditions. DEXUS works also for known
conditions, for example for RNA-Seq data with two or multiple
conditions. RNA-Seq read count data can be provided both by the S4
class Count Data Set and by read count matrices. Differentially
expressed transcripts can be visualized by heatmaps, in which unknown
conditions, replicates, and samples groups are also indicated. This
software is fast since the core algorithm is written in C. For very
large data sets, a parallel version of DEXUS is provided in this
package. DEXUS is a statistical model that is selected in a Bayesian
framework by an EM algorithm. DEXUS does not need replicates to detect
differentially expressed transcripts, since the replicates (or
conditions) are estimated by the EM method for each transcript. The
method provides an informative/non-informative value to extract
differentially expressed transcripts at a desired significance level
or power.

- [DriverNet](http://bioconductor.org/packages/release/bioc/html/DriverNet.html):
DriverNet is a package to predict functional important driver genes in
cancer by integrating genome data (mutation and copy number variation
data) and transcriptome data (gene expression data). The different
kinds of data are combined by an influence graph, which is a gene-gene
interaction network deduced from pathway data. A greedy algorithm is
used to find the possible driver genes, which may mutated in a larger
number of patients and these mutations will push the gene expression
values of the connected genes to some extreme values.

- [DrugVsDisease](http://bioconductor.org/packages/release/bioc/html/DrugVsDisease.html):
This package generates ranked lists of differential gene expression
for either disease or drug profiles. Input data can be downloaded from
Array Express or GEO, or from local CEL files. Ranked lists of
differential expression and associated p-values are calculated using
Limma.  Enrichment scores (Subramanian et al. PNAS 2005) are
calculated to a reference set of default drug or disease profiles, or
a set of custom data supplied by the user. Network visualisation of
significant scores are output in Cytoscape format.

- [eiR](http://bioconductor.org/packages/release/bioc/html/eiR.html):
The eiR package provides utilities for accelerated structure
similarity searching of very large small molecule data sets using an
embedding and indexing approach.

- [ensemblVEP](http://bioconductor.org/packages/release/bioc/html/ensemblVEP.html):
Query the Ensembl Variant Effect Predictor via the perl API

- [epigenomix](http://bioconductor.org/packages/release/bioc/html/epigenomix.html):
A package for the integrative analysis of microarray based gene
expression and histone modification data obtained by ChIP-seq. The
package provides methods for data preprocessing and matching as well
as methods for fitting bayesian mixture models in order to detect
genes with differences in both data types.

- [gCMAPWeb](http://bioconductor.org/packages/release/bioc/html/gCMAPWeb.html):
The gCMAPWeb R package provides a graphical user interface for the
gCMAP package. gCMAPWeb uses the Rook package and can be used either
on a local machine, leveraging R's internal web server, or run on a
dedicated rApache web server installation.  gCMAPWeb allows users to
search their own data sources and instructions to generate reference
datasets from public repositories are included with the package. The
package supports three common types of analyses, specifically queries
with 1. one or two sets of query gene identifiers, whose members are
expected to show changes in gene expression in a consistent
direction. For example, an up-regulated gene set might contain genes
activated by a transcription factor, a down-regulated geneset targets
repressed by the same factor. 2. a single set of query gene
identifiers, whose members are expected to show divergent differential
expression (non-directional query). For example, members of a
particular signaling pathway, some of which may be up- some
down-regulated in response to a stimulus. 3. a query with the complete
results of a differential expression profiling experiment. For
example, gene identifiers and z-scores from a previous perturbation
experiment. gCMAPWeb accepts three types of identifiers: EntreIds,
gene Symbols and microarray probe ids and can be configured to work
with any species supported by Bioconductor. For each query submission,
significantly similar reference datasets will be identified and
reported in graphical and tabular form.

- [GENE.E](http://bioconductor.org/packages/release/bioc/html/GENE.E.html): Interactive exploration of matrices in GENE-E.

- [geNetClassifier](http://bioconductor.org/packages/release/bioc/html/geNetClassifier.html):
Comprehensive package to automatically train a multi- class SVM
classifier based on gene expression data. Provides transparent
selection of gene markers, their coexpression networks, and an
interface to query the classifier.

- [GraphPAC](http://bioconductor.org/packages/release/bioc/html/GraphPAC.html):
Identifies mutational clusters of amino acids in a protein while
utilizing the proteins tertiary structure via a graph theoretical
model.

- [HCsnip](http://bioconductor.org/packages/release/bioc/html/HCsnip.html):
Decompose given hierarchical clustering tree into non-overlapping
clusters in a semi-supervised way by using available patients
follow-up information as guidance. Contains functions for snipping HC
tree, various cluster quality evaluation criteria, assigning new
patients to one of the two given HC trees, testing the significance of
clusters with permutation argument and clusters visualization using
sample's molecular entropy.

- [HTSFilter](http://bioconductor.org/packages/release/bioc/html/HTSFilter.html):
This package implements a filtering procedure for replicated
transcriptome sequencing data based on a global Jaccard similarity
index in order to identify genes with low, constant levels of
expression across one or more experimental conditions.

- [iBMQ](http://bioconductor.org/packages/release/bioc/html/iBMQ.html):
integrated Bayesian Modeling of eQTL data

- [illuminaio](http://bioconductor.org/packages/release/bioc/html/illuminaio.html):
Tools for parsing Illumina's microarray output files, including IDAT.

- [jmosaics](http://bioconductor.org/packages/release/bioc/html/jmosaics.html):
jmosaics detects enriched regions of ChIP-seq data sets jointly.

- [KEGGREST](http://bioconductor.org/packages/release/bioc/html/KEGGREST.html):
A package that provides a client interface to the KEGG REST
server. Based on KEGGSOAP by J. Zhang, R. Gentleman, and Marc Carlson,
and KEGG (python package) by Aurelien Mazurie.

- [lpNet](http://bioconductor.org/packages/release/bioc/html/lpNet.html):
lpNet takes perturbation data as input and generates an LP model which
allows the inference of signaling networks. For parameter
identification either leave-one-out cross-validation or stratified
n-fold cross-validation can be used.

- [metagenomeSeq](http://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html):
metagenomeSeq is designed to determine features (be it Operational
Taxanomic Unit (OTU), species, etc.) that are differentially abundant
between two or more groups of multiple samples. metagenomeSeq is
designed to address the effects of both normalization and
under-sampling of microbial communities on disease association
detection and the testing of feature correlations.

- [MethylSeekR](http://bioconductor.org/packages/release/bioc/html/MethylSeekR.html):
This is a package for the discovery of regulatory regions from Bis-seq
data

- [MineICA](http://bioconductor.org/packages/release/bioc/html/MineICA.html):
The goal of MineICA is to make easier the interpretation of the
interpretation of a decomposition obtained by Independent Component
Analysis on transcriptomic data. It helps the biological
interpretation of the components by studying their association with
variables (e.g sample annotations) and gene sets, and enables the
comparison of components from different datasets using
correlation-based graph.

- [MMDiff](http://bioconductor.org/packages/release/bioc/html/MMDiff.html):
This package detects statistically significant difference between read
enrichment profiles in different ChIP-Seq samples. To take advantage
of shape differences it uses Kernel methods (Maximum Mean Discrepancy,
MMD).

- [PAPi](http://bioconductor.org/packages/release/bioc/html/PAPi.html):
The Pathway Activity Profiling - PAPi - is an R package for predicting
the activity of metabolic pathways based solely on a metabolomics data
set containing a list of metabolites identified and their respective
abundances in different biological samples. PAPi generates hypothesis
that improves the final biological interpretation. See Aggio, R.B.M;
Ruggiero, K.  and Villas-Boas, S.G. (2010) - Pathway Activity
Profiling (PAPi): from metabolite profile to metabolic pathway
activity. Bioinformatics.

- [PathNet](http://bioconductor.org/packages/release/bioc/html/PathNet.html):
PathNet uses topological information present in pathways and
differential expression levels of genes (obtained from microarray
experiment) to identify pathways that are 1) significantly enriched
and 2) associated with each other in the context of differential
expression. The algorithm is described in: PathNet: A tool for pathway
analysis using topological information. Dutta B, Wallqvist A, and
Reifman J. Source Code for Biology and Medicine 2012 Sep 24;7(1):10.

- [pathview](http://bioconductor.org/packages/release/bioc/html/pathview.html):
Pathview is a tool set for pathway based data integration and
visualization. It maps and renders a wide variety of biological data
on relevant pathway graphs. All users need is to supply their data and
specify the target pathway. Pathview automatically downloads the
pathway graph data, parses the data file, maps user data to the
pathway, and render pathway graph with the mapped data. In addition,
Pathview also seamlessly integrates with pathway and gene set analysis
tools for large-scale and fully automated analysis.

- [piano](http://bioconductor.org/packages/release/bioc/html/piano.html):
Piano performs gene set analysis using various statistical methods,
from different gene level statistics and a wide range of gene-set
collections . Futhermore, the Piano package contains functions for
combining the results of multiple runs of gene set analyses.

- [plrs](http://bioconductor.org/packages/release/bioc/html/plrs.html):
The present package implements a flexible framework for modeling the
relationship between DNA copy number and gene expression data using
Piecewise Linear Regression Splines (PLRS).

- [prebs](http://bioconductor.org/packages/release/bioc/html/prebs.html):
The prebs package aims at making RNA-sequencing (RNA-seq) data more
comparable to microarray data. The comparability is achieved by
summarizing sequencing-based expressions of probe regions using a
modified version of RMA algorithm. The pipeline takes mapped reads in
BAM format as an input and produces either gene expressions or
original microarray probe set expressions as an output.

- [pRoloc](http://bioconductor.org/packages/release/bioc/html/pRoloc.html):
This package implements pattern recognition techniques on
quantitiative mass spectrometry data to infer protein sub-cellular
localisation.

- [proteinProfiles](http://bioconductor.org/packages/release/bioc/html/proteinProfiles.html):
Significance assessment for distance measures of time- course protein
profiles

- [pvca](http://bioconductor.org/packages/release/bioc/html/pvca.html):
This package contains the function to assess the batch sources by
fitting all "sources" as random effects including two-way interaction
terms in the Mixed Model(depends on lme4 package) to selected
principal components, which were obtained from the original data
correlation matrix.  This package accompanies the book "Batch Effects
and Noise in Microarray Experiements, chapter 12.

- [QuasR](http://bioconductor.org/packages/release/bioc/html/QuasR.html):
This package provides a framework for the quantification and analysis
of Short Reads. It covers a complete workflow starting from raw
sequence reads, over creation of alignments and quality control plots,
to the quantification of genomic regions of interest.

- [rBiopaxParser](http://bioconductor.org/packages/release/bioc/html/rBiopaxParser.html):
Parses BioPAX files and represents them in R, at the moment BioPAX
level 2 and level 3 are supported.

- [Rbowtie](http://bioconductor.org/packages/release/bioc/html/Rbowtie.html):
This package provides an R wrapper around the popular bowtie short
read aligner and around SpliceMap, a de novo splice junction discovery
and alignment tool. The package is used by the QuasR bioconductor
package.  We recommend to use the QuasR package instead of using
Rbowtie directly.

- [RIPSeeker](http://bioconductor.org/packages/release/bioc/html/RIPSeeker.html):
Infer and discriminate RIP peaks from RIP-seq alignments using
two-state HMM with negative binomial emission probability. While
RIPSeeker is specifically tailored for RIP-seq data analysis, it also
provides a suite of bioinformatics tools integrated within this
self-contained software package comprehensively addressing issues
ranging from post-alignments processing to visualization and
annotation.

- [RNASeqPower](http://bioconductor.org/packages/release/bioc/html/RNASeqPower.html):
RNA-seq, sample size

- [ROntoTools](http://bioconductor.org/packages/release/bioc/html/ROntoTools.html):
Suite of tools for functional analysis

- [RSVSim](http://bioconductor.org/packages/release/bioc/html/RSVSim.html):
RSVSim is a package for the simulation of deletions, insertions,
inversion, tandem-duplications and translocations of various sizes in
any genome available as FASTA-file or BSgenome data package. SV
breakpoints can be placed uniformly accross the whole genome, with a
bias towards repeat regions and regions of high homology (for hg19) or
at user-supplied coordinates.

- [rTANDEM](http://bioconductor.org/packages/release/bioc/html/rTANDEM.html):
This package encapsulate X!Tandem in R. In its most basic
functionality, this package allows to call tandem(input) from R, just
as tandem.exe /path/to/input.xml would be used to run X!Tandem from
the command line. Classes are also provided for taxonomy and
parameters objects and methods are provided to convert xml files to R
objects and vice versa. This package is the first step in an attempt
to provide a reliable worflow for proteomics analysis in R.

- [SANTA](http://bioconductor.org/packages/release/bioc/html/SANTA.html):
This package provides methods for measuring the strength of
association between a network and a phenotype. It does this by
measuring clustering of the phenotype across the network. Vertices can
also be individually ranked by their strength of association with
high-weight vertices.

- [SeqArray](http://bioconductor.org/packages/release/bioc/html/SeqArray.html):
Big data management of genome-wide variants using the CoreArray
library, where genotypic data and annotations are stored in an array-
oriented manner, offering efficient access of genetic variants using
the R language.

- [SeqGSEA](http://bioconductor.org/packages/release/bioc/html/SeqGSEA.html):
Gene set enrichment analysis of high-throughput RNA-Seq data by
integrating differential expression and splicing. Using negative
binomial distribution to model read count data, which accounts for
sequencing biases and biological variation. Based on permutation,
significance analysis can also be done regarding each gene's
differential expression and splicing, respectively.

- [SNAGEE](http://bioconductor.org/packages/release/bioc/html/SNAGEE.html):
Signal-to-Noise applied to Gene Expression Experiments.
Signal-to-noise ratios can be used as a proxy for quality of gene
expression studies and samples. The SNRs can be calculated on any gene
expression data set as long as gene IDs are available, no access to
the raw data files is necessary. This allows to flag problematic
studies and samples in any public data set.

- [SomatiCA](http://bioconductor.org/packages/release/bioc/html/SomatiCA.html):
SomatiCA is a software suite that is capable of identifying,
characterizing, and quantifying somatic CNAs from cancer genome
sequencing. First, it uses read depths and lesser allele frequencies
(LAF) from mapped short sequence reads to segment the genome and
identify candidate CNAs. Second, SomatiCA estimates the admixture rate
from the relative copy-number profile of tumor-normal pair by a
Bayesian finite mixture model. Third, SomatiCA quantifies absolute
somatic copy-number and subclonality for each genomic segment to guide
its characterization. Results from SomatiCA can be further integrated
with single nucleotide variations (SNVs) to get a better understanding
of the tumor evolution.

- [SPEM](http://bioconductor.org/packages/release/bioc/html/SPEM.html):
This package can optimize the parameter in S-system models given time
series data

- [SplicingGraphs](http://bioconductor.org/packages/release/bioc/html/SplicingGraphs.html):
This package allows the user to create, manipulate, and visualize
splicing graphs and their bubbles based on a gene model for a given
organism. Additionally it allows the user to assign RNA- seq reads to
the edges of a set of splicing graphs, and to summarize them.

- [triplex](http://bioconductor.org/packages/release/bioc/html/triplex.html):
This package provides functions for identification and visualization
of potential intramolecular triplex patterns in DNA sequence.  The
main functionality is to detect the positions of subsequences capable
of folding into an intramolecular triplex (H-DNA) in a much larger
sequence.  The potential H-DNA (triplexes) should be made of as many
canonical nucleotide triplets as possible. The package includes
visualization showing the exact base-pairing in 1D, 2D or 3D.

- [UniProt.ws](http://bioconductor.org/packages/release/bioc/html/UniProt.ws.html):
A collection of functions for retrieving, processing and repackaging
the Uniprot web services.

- [wateRmelon](http://bioconductor.org/packages/release/bioc/html/wateRmelon.html):
15 flavours of betas and three performance metrics, with methods for
objects produced by methylumi, minfi and IMA packages.

NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages. The following package NEWS is available:

ADaCGH2
-------

Changes in version 2.1.3 (2010-10-06):

- Rmpi in "Enhances", not "Suggests", to allow for R CMD check in Mac
  and Windows.

Changes in version 2.1.2 (2010-10-04):

- rsprng no longer in depends; L'Ecuyer as default random number
  generator.

- Using R's registration mechanism for C routines.

- Decreased size of example to speed up R CMD check.

- SOCK is default cluster, and Rmpi not loaded by us.

Changes in version 2.1.1 (2010-09-30):

- Output to CGHregions and input from limma and snapCGH. Changes in
  functions, help, vignnette.

- Can also use rlecuyer.

- Works with R-2.11 (adapted to differences in "inherits").

Changes in version 2.1.0 (2010-09-23):

- First fully working version for BioC. Versioning changed!

Changes in version 1.9.0 (2013-03-23):

- Commented out unused code that gave warnings in checks.

- Removed partial matching (was giving notes).

- Note: version 2.0 is almost ready, so work on current version is now
  limited.

affxparser
----------

Changes in version 1.31.4 (2013-03-19):

- Made example(invertMap) a bit faster so 'R CMD check' won't complain.

Changes in version 1.31.3 (2013-03-18):

- CLEANUP: Internal isPackageLoaded() of findFiles() no longer uses
  defunct manglePackageName() function.

Changes in version 1.31.2 (2013-01-07):

- Same updates as in release v1.30.2.

Changes in version 1.31.1 (2012-10-18):

- Now compareCdfs() gives a more precise 'reason' attribute when there
  is a difference in (regular or QC) units.  It narrows down the first
  unit that differs and reports it unit number.

Changes in version 1.31.0 (2012-10-01):

- The version number was bumped for the Bioconductor devel version.

Changes in version 1.30.2 (2013-01-07):

- BUG FIX: writeCdf() did not encode unit types as decoded by
  readCdf().  Unit type 'unknown' was incorrectly encoded such that
  readCdf() would decode it as 'copynumber'.  Also, unit types
  'genotypingcontrol' and 'expressioncontrol' where not encoded at all.

affyPara
--------

Changes in version 1.19.2:

- more Sweave bugfix to get a clean build

Changes in version 1.19.1:

- bugfix to get a clean build

annmap
------

Changes in version 1.1.1:

- IMPROVED genomicPlot when passed a .genes data.frame used to call
  geneDetails on it to ensure the data was in the correct format. This
  has been removed, so that you can pass custom regions in and have
  them plotted on the graph.

AnnotationDbi
-------------

Changes in version 1.22:

NEW FEATURES

- There is now a convenience function for extracting data from the
  GO.db package as a graph object. The function is called: makeGOGraph.

- Adds supportedSeqnameMappings() and findSequenceRenamingMaps()
  utilities.

- Improved vignette for new users

BUG FIXES AND CODE MAINTENANCE

- Fixes a bug in select for users accessing reactome.db

- Fixes a bug in select for users requesting via a PROBEID

- Fixes a bug in select for users using rat chip packages

- Fixes a bug in mget, for Bimap objects (when ifnotfound=NA) * * *
  1.20.x SERIES NEWS * * *

aroma.light
-----------

Changes in version 1.29.0 (2012-10-01):

- The version number was bumped for the Bioconductor devel version.

ASEB
----

Changes in version 1.4.0:

- documentation improvements

Changes in version 1.3.2:

- documentation improvements

Changes in version 1.3.1:

- add release notes and documentation improvements

Changes in version 1.3.0:

- documentation improvements

BaseSpaceR
----------

Changes in version 1.0.0:

PACKAGE FEATURES

- BaseSpaceR is a R SDK for Illumina's BaseSpace cloud computing
  environment, enabling the rapid prototyping and development of
  production ready application for next-gen sequencing data.

PACKAGE FEATURES

- It provides a set of S4 classes and methods to interface with
  BaseSpace data model.

PACKAGE FEATURES

- It allows for persistent connections with BaseSpace REST server and
  offers support for the REST API query parameters.

PACKAGE FEATURES

- It allows for queries across multiple Projects, Samples, Files,
  AppResults, etc., using vectorized operations in line with the R
  semantic.

Biobase
-------

Changes in version 2.19:

USER VISIBLE CHANGES

- dimnames(), rownames(), colnames() and setters work on eSet-derived
  objects

BiocInstaller
-------------

Changes in version 1.10.0:

NEW FEATURES

- biocValid() checks that installed packages are consistent with those
  available via biocLite().

- biocVersion() returns the version of Bioconductor expected with this
  version of the BiocInstaller package.

USER-VISIBLE CHANGES

- biocLite() invoked with no arguments updates currently installed
  packages to their most-recent version.

biomvRCNS
---------

Changes in version 0.99.3:

NEW FEATURES

- add mvt and t to the emis.type and tmvtFit function

BUG FIXES

- fix shift index (-1) in pois and nbinom fitting

Changes in version 0.99.2:

NEW FEATURES

- add log-Vertibi in c

- add normalizing factor for estimation of emis$p and soj$d

- add avgFunc to be used in average calculation, median and mean with
  trim, default to median

BUG FIXES

- fix shift index (-1) in pois and nbinom fitting

Changes in version 0.99.1:

NEW FEATURES

- in Bioc-devel

BitSeq
------

Changes in version 1.3.11 (2013-04-02):

NEW FEATURES

- new functions for handling transcript information

NEW FEATURES

- getExpression returns also counts and transcript information

Changes in version 1.3.10 (2013-03-30):

NEW FEATURES

- proper handling of alignments from Bowtie2

Changes in version 1.3.6 (2013-03-19):

NEW FEATURES

- enabling prallelization in getExpression via OpenMP

- gene names can be extracted from Ensembl-like reference while
  estimating expression

- added seed option for estimateHyperPar and estimateDE

- improved output for getDE when using more than 2 conditions

BUG FIXES

- problem with long lines in reference sequence file in parseAlignment

Changes in version 1.2.2 (2013-01-29):

- occasional crash when using non-uniform read model in getExpression

Changes in version 1.2.1 (2012-11-06):

- parsing gapped alignments and half-alignments of paired end reads

BSgenome
--------

Changes in version 1.28.0:

NEW FEATURES

- Add seqnames() setter for BSgenome objects so users can rename the
  single sequences in those objects. This has been a popular user
  request for a while.

SIGNIFICANT USER-VISIBLE CHANGES

- Add sanity check to "getSeq" method for BSgenome objects, raising an
  error if the supplied BSgenome and GRanges objects are based on
  incompatible reference genomes.

MISCELLANEOUS

- Started the NEWS file (this file).

bsseq
-----

Changes in version 0.7:

- Removed the returnRaw argument to read.bismark() as it was
  unnecessary (Bismark output files does not have additional
  information beyond M and Cov and genomic positions, unlike BSmooth).

- Moved the Bismark example data from data to inst/extdata.

- combineList() now deals with the case where the list of BSseq objects
  have different genomic locations.  This speeds up read.bismark()
  substantially.

- Exposed combineList() as a faster alternative to Reduce(combine,
  list).

- Updated the code for the plotting routines (plotRegion).  This should
  not have an impact on user-visible code.

- Added read.bismark() function to parse output from the Bismark
  alignment suit [thanks to Pete Hickey].

- Refactorized plotting code.

CellNOptR
---------

Changes in version 1.6.0 (2013-03-12):

- BUG FIX
  
  * readMIDAS: DV, DA and TR can now be in the specy name
  
  * makeCNOlist: bug fixed when only one experiment was present
  
  * plotModel: figures where in B and W due to issue in Rgraphviz
  package
  
  * gaBinaryT family: (1) fix bug when only one model found within the
  tolerance. (2) a faster hash table used in the optimisation.

- CHANGES
  
  * writeSIF: can overwrite existing file with a parameter
  
  * remove cutAndPlotResultsT2: use cutAndPlotResultsTN instead or
  cutAndPlot
  
  * plotOptimResultsPan: NA are now rendered in gray and species are
  rectangle instead of ellipses. The black line connecting measurements
  is printed.  The ylim for the y axis can be set manually to overwrite
  the default behaviour.
  
  * readSIF can now read relation "A 1 B C D E" as excepted in SIF
  format.
  
  * makeCNOlist/CNOlist: variances available in the structure

- NEW FEATURES
  
  * CNOdata: a function to fetch data from www.cellnopt.org
  
  * plotModel
  
  * exhaustive function: a simple exhaustive function to perform
  optimisation for small models
  
  * cutCNOlist: cut a cnolist given a list of species
  
  * CNOlist: 3 new methods: length, randomize and plot accepts either
  one or 2 cnolist arguments.
  
  * randomizeCNOlist: function to perform different randomization of
  the data
  
  * model2igraph function
  
  * compatCNOlist to convert CNOlist back to old style; used by CNORode

CGEN
----

Changes in version 2.0:

USER VISIBLE CHANGES

- The function additive.test has been added

BUG FIXES

- The function snp.logisitic checks that the snp variable is used in
  main.vars or int.vars

PLANS

- The next release will have new features for snp.logistic and
  snp.matched, along with some new functions.

ChemmineR
---------

Changes in version 2.12.0:

NEW FEATURES

- Accelerated similarity searching of large small molecule data sets
  via new eiR add-on package

NEW FEATURES

- Jarvis-Patrick clustering of large small molecule data sets

NEW FEATURES

- SQLite support for small molecule management

cisPath
-------

Changes in version 1.0.0:

- documentation improvements

Changes in version 0.99.10:

- documentation improvements

Changes in version 0.99.9:

- documentation improvements

Changes in version 0.99.8:

- add the easyEditor method

Changes in version 0.99.7:

- unordered_map -> unordered_map (gcc > 4.6.0) or hash_map

Changes in version 0.99.6:

- Several improvements have been done.

Changes in version 0.99.5:

- Several improvements have been done.

Changes in version 0.99.4:

- A button has been added which can be used to freeze the graph.

Changes in version 0.99.3:

- hash_map -> unordered_map.

Changes in version 0.99.2:

- The size of output was reduced to about 30%.

Changes in version 0.99.1:

- A new method called networkView has been added to this new package.

Changes in version 0.99.0:

- Package released

clusterProfiler
---------------

Changes in version 1.7.1:

- support GO enrichment analysis for organism celegans <2013-01-22,
  Fri>

- add compress parameter in buildGOmap <2013-01-16, Wed>

- import ggtitle from ggplot2 <2012-09-07, Fri>

- update codes of plot functions accompaning with ggplot2 (version
  0.9.2) <2012-09-06, Thu>

- bug fixed of buildGOmap due to the empty GO annotation query from
  biomaRt <2012-07-18, Wed>

cqn
---

Changes in version 1.5:

- THe output object from cqn() now has an additional component:
  glm.offset which is an offset matrix which can directly be used in a
  GLM type model (specifically edgeR).  The usage is explained in the
  vignette secion on Import into edgeR.  Previously the vignette
  recommended using the offset component of the cqn output, which is
  wrong, due to a scaling issue.  The offset component of cqn is
  unchanged.  This bug was found by Mike Love <email:
  love@molgen.mpg.de> and fixed in CQN 1.5.1.

ddgraph
-------

Changes in version 1.3.1:

- Removed dependency on a broken infotheo package

Changes in version 1.3.0:

- Bioconductor 2.11 release version

DEGseq
------

Changes in version 1.14.0:

- documentation improvements

Changes in version 1.13.3:

- add release notes and documentation improvements

Changes in version 1.13.2:

- documentation improvements

Changes in version 1.13.1:

- documentation improvements

Changes in version 1.13.0:

- documentation improvements

DESeq2
------

Changes in version 1.0.0:

- Base class: SummarizedExperiment is used as the superclass for
  storing the data.

- Workflow: The wrapper function DESeq() performs all steps for a
  differential expression analysis. Individual steps are still
  accessible.

- Statistics: Incorporation of prior distributions into the estimation
  of dispersions and fold changes (empirical-Bayes shrinkage). A Wald
  test for significance is provided as the default inference method,
  with the likelihood ratio test of the previous version also
  available.

- Normalization: it is possible to provide a matrix of sample- *and*
  gene-specific normalization factors

DEXSeq
------

Changes in version 2013-02-27:

- A parameter -r was added to the python scripts, that allow the users
  either to ignore the exonic bins belonging to several genes and treat
  the genes separately, or merge the genes into an aggregate gene. The
  equivalent R implementations of the python scripts were finally
  added.

Changes in version 2012-11-28:

- The TRT method is implemented, for people with a big number of
  samples without completions or speed issues

DiffBind
--------

Changes in version 1.6.0:

- New: Low memory counting of bam files using Rsamtools and
  summarizeOverlaps (bLowMem in dba.count)

- New: Ability to read in externally derived counts (e.g. from htSeq)
  (dba.count)

- Improved: Features to deal with filtering intervals based on read
  scores (dba.count)
  
  * Change parameter name: maxFilter -> filter
  
  * Allow maxFilter to be a numerical vector to retrieve filtering rate
  
  * Add parameter: filterFun to control filtering method

- New: Support for SummarizedExperiment objects (dba and dba.report)
  
  * Add bSummarizedExperiment option to dba() to convert DBA object
  
  * Add DataType = DBA_DATA_SUMMARIZED_REPORT option to dba.report() to
  return SummarizedExperiment

- Documentation: Add section to vignette showing how to obtain full
  tamoxifen resistance dataset
  
  * Add section to vignette showing how to obtains full tamoxifen
  dataset
  
  * Add script (tamoxifen_GEO.R) and sample sheet (tamoxifen_GEO.csv)
  to extras for full tamoxifen dataset
  
  * Add examples to man page for dba.count to show filtering
  
  * Add examples to man pages for dba and dba.report to show retrieval
  of SummarizedExperiment objects
  
  * Update and cleanup vignette and man pages

- Various bugfixes and improved warnings

DOSE
----

Changes in version 1.5.1:

- bug fixed in enrich.internal, now return NA rather than throw error,
  if gene have no ontology annotation <2013-01-22, Fri>

DSS
---

Changes in version 1.2.0-1:

- Fixed a bug in computing local FDR.

- Change the way to deal with genes with all 0 counts. Now the FDRs for
  these genes are assigned as 0.

easyRNASeq
----------

Changes in version 1.5.1:

- Adapted the dependencies version to match the Bioconductor release
  version 2.11 BUG FIXES

- corrected an innapropriate function call to an internal function (as
  in stable version 1.4.2)

Changes in version 1.5.0:

- No changes, Bioconductor development version 2.12

Changes in version 1.4.2:

BUG FIXES

- corrected an innapropriate function call to an internal function

Changes in version 1.4.1:

- Adapted the dependencies version to match the Bioconductor release
  version 2.11

EBImage
-------

Changes in version 4.2.0:

NEW FEATURES

- 'localCurvature' function for computing local curvature along a line
  (J. Barry)

SIGNIFICANT USER-VISIBLE CHANGES

- the range of pixel coordinates displayed in the JavaScript viewer is
  now (1,1):(w,h) rather than (0,0):(w-1,h-1) and matches the indices
  of the corresponding Image array

BUG FIXES

- 'erode'/'dilate': fixed a bug introduced in the previous version
  (4.0.0)

- 'resize': new image width was calculated incorrectly when only height
  was provided (reported by B. Fischer)

- 'medianFilter': incorrect [0:1] <-> integer range conversion (thanks
  to K. Johnson)

edgeR
-----

Changes in version 3.2.0:

- The User's Guide has a new section on between and within subject
  designs and a new case study on RNA-seq profiling of unrelated
  Nigerian individuals. Section 2.9 (item 2) now gives a code example
  of how to pre-specify the dispersion value.

- New functions estimateDisp() and WLEB() to automate the estimation of
  common, trended and tagwise dispersions. The function estimateDisp()
  provides a simpler alternative pipeline and in principle replaces all
  the other dispersion estimation functions, for both glms and for
  classic edgeR. It can also incorporate automatic estimation of the
  prior degrees of freedom, and can do this in a robust fashion.

- glmLRT() now permits the contrast argument to be a matrix with
  multiple columns, making the treatment of this argument analogous to
  that of the coef argument.

- glmLRT() now has a new F-test option. This option takes into account
  the uncertainty with which the dispersion is estimated and is more
  conservative than the default chi-square test.

- glmQLFTest() has a number of important improvements. It now has a
  simpler alternative calling sequence: it can take either a fitted
  model object as before, or it can take a DGEList object and design
  matrix and do the model fit itself. If provided with a fitted model
  object, it now checks whether the dispersion is of a suitable type
  (common or trended). It now optionally produces a plot of the raw and
  shrunk residual mean deviances versus AveLogCPM. It now has the
  option of robustifying the empirical Bayes step. It now has a more
  careful calculation of residual df that takes special account of
  cases where all replicates in a group are identically zero.

- The gene set test functions roast(), mroast() and camera() now have
  methods defined for DGEList data objects. This facilitates gene set
  testing and pathway analysis of expression profiles within edgeR.

- The default method of plotMDS() for DGEList objects has changed. The
  new default forms log-counts-per-million and computes Euclidean
  distances. The old method based on BCV-distances is available by
  setting method="BCV". The annotation of the plot axes has been
  improved so that the distance method used is apparent from the plot.

- The argument prior.count.total used for shrinking log-fold-changes
  has been changed to prior.count in various functions throughout the
  package, and now refers to the average prior.count per observation
  rather than the total prior count across a transcript. The treatment
  of prior.counts has also been changed very slightly in cpm() when
  log=TRUE.

- New function aveLogCPM() to compute the average log count per million
  for each transcript across all libraries. This is now used by all
  functions in the package to set AveLogCPM, which is now the standard
  measure of abundance.  The value for AveLogCPM is now computed just
  once, and not updated when the dispersion is estimated or when a
  linear model is fitted. glmFit() now preserves the AveLogCPM vector
  found in the DGEList object rather than recomputing it. The use of
  the old abundance measure is being phased out.

- The glm dispersion estimation functions are now much faster.

- New function rpkm() to compute reads per kilobase per million (RPKM).

- New option method="none" for calcNormFactors().

- The default span used by dispBinTrend() has been reduced.

- Various improvements to internal C++ code.

- Functions binCMLDispersion() and bin.dispersion() have been removed
  as obsolete.

- Bug fix to subsetting for DGEGLM objects.

- Bug fix to plotMDS.DGEList to make consistent use of norm.factors.

eiR
---

Changes in version 1.0.0:

NEW FEATURES

- The eiR packages introduces efficient methods for accelerating
  structure similarity searches and clustering of very large compound
  datasets. The acceleration is achieved by applying embedding and
  indexing techniques to represent chemical compounds in a
  high-dimensional Euclidean space and to employ ultra-fast
  pre-screening of the compound dataset using the LSH-assisted nearest
  neighbor search in the embedding space. This method can drastically
  reduce the search time of large databases, by a factor of 40–200 fold
  when searching for the 100 closest compounds to a query.

ExiMiR
------

Changes in version 2.0.1:

- R/createAB.R: Adding support for non imagene data type read by limma

- R/NormiR.R: Fix median/mean normalization issue

exomeCopy
---------

Changes in version 1.6.0 (2013-03-05):

- switched from X*beta to exp(X*beta) for modeling the expected read
  depth

fmcsR
-----

Changes in version 1.1.3:

- Vignette updates

Changes in version 1.1.2:

- Help updates

Changes in version 1.1.1:

- Vignette updates

gage
----

Changes in version 2.9.4:

- add secondary vignette, "Gene set and data preparation", on data
  preparation.

Changes in version 2.9.2:

- suggests and connected to pathview package for results visualization.

Changes in version 2.9.1:

- removed dependency on multtest package for p-value FDR adjustment,
  use p.adjust function of the stat package instead.

- change "depends" to "imports" for graph package as we only need to
  import graphNEL class and connComp method.

- subset kegg.gs. From now, kegg.gs only include the subset of
  canonical signaling and metabolic pathways from KEGG pathway
  database, and kegg.gs.dise is the subset of disease pathways. And it
  is recommended to do KEGG pathway analysis with either kegg.gs or
  kegg.gs.dise seperately (rather than combined altogether) for better
  defined results. Note that kegg.gs and subsets are be defined
  slightly different in gageData package.

- In gage vignette, add citation section and an example of including
  all genes (rather than those selected using essGene function) in top
  gene set result check using geneData function.

GeneAnswers
-----------

Changes in version 2.0.0:

NEW FEATURES

- Reactome pathway analysis is based on local reactome.db. The previous
  xml queries had been disabled.

BUG FIXES

- Completely support igraph. igraph0 is not supported, which means the
  nodes index starts at 1, not 0 in igraph0 or less than 0.6.x version
  igraph.

- GO.CC and GO.MF level specified test had been fixed.

geNetClassifier
---------------

Changes in version 1.0:

NEW FEATURES

- Package released

genomeIntervals
---------------

Changes in version 1.15.3:

- genome intervals order now consisten with assumption that (start ==
  [start-1 and that stop) == stop-1]

Changes in version 1.15.2:

- Depends on intervals >=0.14.0 to fix a change in R's split behavior

- sort does not have byName argument any longer NEW FEATURES

- order, sort, rank and xtfrm consistently implemented.

Changes in version 1.15.1:

- Created that NEWS file to replace the CHANGES file and be compliant
  to the R standard package architecture NEW FEATURES

- introduced a coercion to data.frame

- introduced a writeGff3 function

- reverted the sort behavior to the default R behavior and added a
  method argument. Setting it to byName results in a more biologically
  relevant sorting of the object.

GenomicFeatures
---------------

Changes in version 1.12:

NEW FEATURES

- Support for new UCSC species

BUG FIXES

- Better support for GTF and GFF processing into TranscriptDb objects

NEW FEATURES

- Methods for making TranscriptDb objects from general sources have
  been made more useful

BUG FIXES

- Updates to allow continued access to ever changing services like UCSC

NEW FEATURES

- Corrections for seqnameStyle methods

BUG FIXES

- Over 10X performance gains for processing of GTF and GFF files

GenomicRanges
-------------

Changes in version 1.12.0:

NEW FEATURES

- Implement "seqnameStyle" replacement method for Seqinfo object.
  'seqnameStyle(x) <- style' works on any object with a "seqinfo"
  replacement method.

- Add trim,GenomicRanges-method to trim out of bound ranges.

- Add promoters,GenomicRanges and promoters,GRangesList methods.

- Add "overlapsAny" methods as a replacement for the deprecated "%in%"
  methods.

- Add 'ignore.strand' argument to match,GenomicRanges-method.

- Add 'with.mapping' argument to "reduce" method for GenomicRanges
  objects.

- Add "unname" method to remove dimnames from SummarizedExperiment.

- Add "cbind" and "rbind" methods for SummarizedExperiment.

- Add "seqselect", "seqselect<-" and "split" methods for
  SummarizedExperiment.

- Add GAlignmentsList class.

- Add readGAlignmentsList generic and methods.

SIGNIFICANT USER-LEVEL CHANGES

- resize,GenomicRanges method no longer checks that 'fix' is
  length-compatible with 'x' when 'x' is length zero. This allows for
  resize(x, w, fix = "end") without worrying about 'x' being
  zero-length.

- Change the behavior of "distance". Previously adjacent ranges had a
  distance of 1 and overlapping had a distance of 0. Now both adjacent
  AND overlapping have a distance of 0.

- shift,GenomicRanges-method no longer trims out of bound ranges.

- "distanceToNearest" no longer drops ranges that have no hit but
  returns 'NA' for 'subjectHits' and 'distance'.

- "genome" is no longer an invalid metadata colname for GenomicRanges
  objects.

- 4x-8x speedup for doing coverage() on a GRanges or GRangesList with
  many seqlevels.

- Remove ">=", "<", and ">" methods for GenomicRanges objects.

- Speedup "seqinfo" setters for GenomicRanges and GappedAlignments by
  avoiding validation when not necessary.

- readGappedAlignments can now pass a BamFile to
  readBamGappedAlignments.

- Remove unneeded "unique" and "sort" methods for GenomicRanges
  objects.

- Change behavior of "match" and "%in%" on GenomicRanges objects to use
  equality instead of overlap for comparing elements between
  GenomicRanges objects 'x' and 'table'.

- match,GenomicRanges-method gets the same 'method' argumnet as the
  "duplicated" method for these objects.

- Remove unneeded "countOverlaps" methods.

- "classNameForDisplay" shortens the name of data type when displayed.

- Add global options 'showHeadLines' and 'showTailLines' to control the
  number of head/tails lines displayed in show,GRanges and
  show,GappedAlignments methods.

- "distanceToNearest" now returns a Hits object instead of DataFrame.

DEPRECATED AND DEFUNCT

- Remove defunct countGenomicOverlaps(), grg(), and globalToQuery()

- Defunct previously deprecated '.ignoreElementMetadata' argmuent of
  c,GenomicRanges-method.

- Deprecate all "match" and "%in%" methods in the package except for
  those with the GenomicRanges,GenomicRanges signature.

- Deprecate "resolveHits" methods.

BUG FIXES

- Several bug fixes to "nearest".

- Output of "findSpliceOverlaps" now displays 'NA' for ranges with no
  hits.

GGtools
-------

Changes in version 4.7:

- The snplocsDefault() function has been added to simplify appropriate
  selection of SNPlocs.Hsapiens.dbSNP.* to a common value for all
  usages

- The sensanal() function now operates on a sensiCisInput instance to
  help provide an overview of sensitivity analysis for cis-eQTL
  searches

gmapR
-----

Changes in version 1.2.0:

NEW FEATURES

- New method getSeq,GmapGenome retrieves sequence from a GmapGenome
  index. This also supports a coercion to DNAStringSet and thus easy
  export to FASTA via rtracklayer.

- bam_tally gains an ignore_duplicates argument for ignoring BAM
  records flagged as PCR/optical duplicates.

- Read position mean and variance are now output by bam_tally.

USER-VISIBLE CHANGES

- GMAP has been updated to the July '12 version (yes, this is old).

- GSTRUCT (bamtally) updated to trunk as of 3/22/13.

BUG FIXES

- asBam,GsnapOutput now actually works.

GOSemSim
--------

Changes in version 1.17.1:

- update IC data for next release <2013-03-08, Fri>

- after removing NA row/col of similarity matrix, if only one row/col
  remains, R will turn it to be a vector, and combineScore function
  will not work properly. This bug was fixed <2013-01-11, Fri>

- bug fixed of infoContentMethod, now return NA when ID is not belong
  to the ontology <2012-10-11, Thu>

GraphPAC
--------

Changes in version 1.0.0:

- First release of the GraphPAC package.

- Two plotting types available.

- Insertion methods allowed are cheapest, nearest, farthest and random.

GSEABase
--------

Changes in version 1.21:

SIGNIFICANT USER-VISIBLE CHANGES

- GeneSetCollection,*,*,GOCollection-method respects evidenceCode and
  ontology

Gviz
----

Changes in version 1.4.0:

NEW FEATURES

- BiomartGeneRegionTracks will now make use of the available CDS
  information in Ensembl.

- The constructors to the AnnotationTrack, GeneRegionTrack, DataTrack
  and SequenceTrack classes now accept a character scalar that points
  to a file on the file system. A number of default parser functions
  have been implemented to read the standard file types. Alternatively,
  a user-defined import function can be provided. This feature also
  supports streaming from indexed file types like BAM or bigWig, in
  which case the data is fetched dynamically upon each plotting
  operation.

- The mart object in BiomartGeneRegionTrack objects is now cached in
  order to speed up subsequent queries to the same mart.

- When plotting DataTracks with type 'gradient' or 'heatmap', a color
  scale is plotted next to the regular y-axis to indicate the mapping
  of numeric values in the false color range. Thanks to Mark Heron for
  his code contribution.

- Sample names can now be shown in heatmap-type plots by setting the
  'showSampleNames' display parameter.

SIGNIFICANT USER-VISIBLE CHANGES

- Complete refactoring of the automatic font size adjustments to
  provide more reasonable defaults.

- Tick labels on the genomic axis are now show in between tick marks
  when zoomed in to single nucleotide level.

BUG FIXES

- Fixed a bug in IdeogramTracks where all bands in the rounded caps at
  the edges of the Ideogram were missing.

- The way genomic ranges are plotted is now according to the Lego block
  model suggested by Herve. This is only relevant when zooming in to
  the level of single nucleotides.

- Tick labels on the genome axis show only significant digits now.

- Sample ordering in heatmap plots is now correct.

- Numerous little fixes.

GWASTools
---------

Changes in version 1.5.9:

- assocTestRegression computes allele counts separately for each model.

- convertNcdfGds uses information from a SnpAnnotationDataFrame to
  store allele and chromosome codes in the GDS file.

Changes in version 1.5.8:

- Adding missing value support to GdsReader.

- Fixed bug in getAttribute method for GdsReader.

- Updated GdsReader for compatibility with gdsfmt 0.9.11 (no longer
  compatible with older versions).

Changes in version 1.5.7:

- Fixed bug in genotypeToCharacter that resulted in calls to
  getGenotype(char=TRUE) for a single SNP to return NA.

- Renamed minorAlleleSensitivitySpecificity to
  minorAlleleDetectionAccuracy and added additional output.

Changes in version 1.5.6:

- Added function minorAlleleSensitivitySpecificity.

Changes in version 1.5.5:

- Deprecated pedigreeClean and pedigreeFindDuplicates.  pedigreeCheck
  now encompasses all pedigree checks and should be used instead.

- Added pedigreeMaxUnrelated to find the maximum set of unrelated
  members of a pedigree.

- Added additional output column "MAF" to matrix returned by
  alleleFrequency.

Changes in version 1.5.4:

- Removed hard-coding of autosomes as 1:22; can now set a vector of
  integer codes corresponding to autosomes with "autosomeCode" argument
  at object creation and retrieve with "autosomeCode" methods.  This
  change makes GWASTools compatible with non-human organisms.

- Added option to duplicateDiscordanceAcrossDatasets to count missing
  data as discordance.

- Added option to start axes of genoClusterPlot at 0.

Changes in version 1.5.3:

- Removed "alleleA.col" and "alleleB.col" options from plink functions,
  as "alleleA" and "alleleB" are now standard names.

- Added "getAlleleA" and "getAlleleB" methods to GdsGenotypeReader.

- Added "getDimension" method to NcdfReader.

Changes in version 1.5.2:

- Added "getAlleleA" and "getAlleleB" methods to SnpAnnotation* and
  GenotypeData objects.

- Added genotypeToCharacter function to convert genotypes from number
  of A alleles to A/B format.

- getGenotype for GenotypeData has option char=TRUE to return character
  genotypes in A/B format.

Changes in version 1.5.1:

- Added option to duplicateDiscordanceAcrossDatasets to calculate minor
  allele discordance.

hapFabia
--------

Changes in version 1.2.0:

- bug fix vcftoFABIA

Changes in version 1.0.5:

- haplotype vcf files are now possible

Changes in version 1.0.4:

- rename in programs and manuals "haplotype cluster" in "IBD segment"

Changes in version 1.0.3:

- smaller improvements

Changes in version 1.0.2:

NEW FEATURES

- new parameter: distance between SNVs for computing bin sizes

NEW FEATURES

- new parameter: bin size in terms of SNVs directly

NEW FEATURES

- default are now genotype data

Changes in version 1.0.1:

NEW FEATURES

- chromosome number automatically extracted from annotation

NEW FEATURES

- annotation file now tab or blank separated (automatically checked)

HiTC
----

Changes in version 1.3.3:

NEW FEATURES

- Add PCA function on Hi-C interaction map as in Lieberman-Aiden et al.
  2009

BUG FIXES

- HTCexp. Error in constructor when the interaction map has a dim = 1

- mapC. Bug fixed in the visualization of two interaction maps of
  different sizes. The scale was updated, and the annotation tracks is
  always drawn from the larger map.

Changes in version 1.3.2:

NEW FEATURES

- New visualization for HTClist objects

- New methods for HTClist objects

- New HTClist class to manage list of HTCexp object (basically Hi-C
  data)

SIGNIFICANT USER-VISIBLE CHANGES

- Update of all man pages

- Update of the Nora_5C data. E14 and MEF are now HTClist objects

- Update of the importC/exportC function. The standard format is now
  matrix-based. This seems to be the most commonly used format.

- Update of all mapC methods. The view parameter is removed. The HTCexp
  object are now displayed in the triangle view, whereas the HTClist
  are displayed in the heatmap view

Changes in version 1.3.1:

NEW FEATURES

- MAJOR RELEASE : replace all GenomeIntervals objects by GRanges ones
  in order to improve the compatibility with other HT BioC packages

SIGNIFICANT USER-VISIBLE CHANGES

- The ExtractRegion method has a new MARGIN parameter. The idea is the
  same than for any apply function. If MARGIN is equal to 1 (resp. 2,
  resp. c(1,2)), the region is extracted from the 'x' (resp. 'y', resp.
  both) intervals

- Plot function. When two HTCexp objects are plot together, only the
  intersection of the 'x' and 'y' intervals are used.

- The 'range' method now returns a GRanges object

- Changes in the data windowing for the extreme bins

- mapC requires a HTCexp object only. Objects from the matrix class are
  no longer allowed

DEPRECATED AND DEFUNCT

- seq_name is now deprecated

- export and normPerZscore are now defunct

BUG FIXES

- exportC. Bug fixed in bin coordinates

- CQC. Bug fixed with NA values

- mapC. Bug fixed in the visualization of annotation features. Select
  the annotation in the same chromosome space before plotting.

- mapC. Bug fixed in the visualization of count values for
  interchromosomal data

hpar
----

Changes in version 1.1.2:

- Vignette update for knitr 1.0 compatibility, thanks Dan! <2013-01-15
  Tue>

Changes in version 1.1.1:

- fixing vignette <2012-10-02 Tue>

Changes in version 1.1.0:

- version bump for new devel <2012-10-01 Mon>

HTSeqGenie
----------

Changes in version 3.9.25:

- VariantTools "analyzeVariants.indels" are OK

Changes in version 3.9.24:

- aloow for new quality score range "GATK-rescaled" from 1-50 (33-83 in
  ASCII)

Changes in version 3.9.23:

- added the config parameter "analyzeVariants.indels"

Changes in version 3.9.22:

- removed the dependency towards the "logging" package

Changes in version 3.9.21:

- variant calling via GATK

Changes in version 3.9.20:

- preparation to BioC submission

- now using detectRRNA.do: FALSE in default-config.txt

Changes in version 3.9.19:

- removed mc.preschedule=FALSE from mergeBAMsAcrossDirs

- added the configuration parameter 'analyzeVariants.method' (GATK
  check has yet to be done)

Changes in version 3.9.18:

- now depends on VariantTools 1.1.13 that fixes the
  mclapply(mc.preschedule=FALSE) bug

Changes in version 3.9.17:

- added the config parameter 'alignReads.analyzedBam' to control how
  analyzed.bam are built

- removed the former config parameter 'alignReads.analyzed_bamregexp'
  that could not work on single ends

Changes in version 3.9.16:

- sessionInfo() is not called anymore in writePreprocessAlignReport()
  during generation of report, to prevent crash when PACKAGES have been
  updated while the pipeline is running

- sclapply() now uses a 'finally' cleanup procedure to kill all threads
  it has created

- added some unit tests to check that no leftover threads are present
  after sclapply() in different scenarios

Changes in version 3.9.15:

- use low lever variant calling interface from VariantTools.  This
  allows for access to the raw_variants as well as the filtered ones/

- variant calling now included in mergeLanes()

Changes in version 3.9.14:

- added the config parameter "alignReads.use_gmapR_gsnap" to control if
  gsnap should be called from gmapR or from the PATH

- default config parameter "alignReads.use_gmapR_gsnap" is now TRUE

- removed the duplicated default config parameters: path.picard_tools,
  markDuplicates.do

- added a check in checkConfig() to stop if some config paramters are
  duplicated

Changes in version 3.9.13:

- add variant calling using VariantTools (not yet parallelized yet)

Changes in version 3.9.12:

- include markDuplicates into runPipeline(), controlled by
  markDuplicates.do config

Changes in version 3.9.11:

- refactor setupTestFramework() to allow for injection of TP53 genome
  template

Changes in version 3.9.10:

- add function to mark duplicates via picard tools

Changes in version 3.9.9:

- fixed detectRRNA code, including bug in wrapGsnap

- add test for detectRRNA working on tp53 genome

Changes in version 3.9.8:

- the system command 'samtools' is no used anymore in the code

- removed unused functions: indexBAMFiles, filterBam,
  getReadLengthFromBam, getBamIndexStats

- (filterBam will be back in the xenograft module)

Changes in version 3.9.7:

- works with Biobase 2.18.0 (Bioconductor release 2.11)

- fixed the "x is not present in the PATH" bogus message

- old gmapR stuffs are now gone: parallelized_gsnap, consolidateSAM,
  consolidateGsnapOutput, consolidateBAm

- now use wrapGsnap, to facilitate the transition to the gsnap offered
  by gmapR

- now depends on gmapR (to load TP53Genome())

Changes in version 3.9.6:

- make remaining tests run with TP53 genome

- move detectRRNA tests to HTSeqGenie.gne as they depend on IGIS

Changes in version 3.9.5:

- remove runPipeline tests depending on IGIS. Instead use simple
  integration test based on TP53 genome. This requried additon of :
  R/runPipeline.R R/TP53GenomicFeatures.R copied from bioc branch and
  dependance on gmapR for the TP53Genome

Changes in version 3.9.4:

- converging with the BioC version: adding @internal keyword

- configuration parameter

Changes in version 3.9.3:

- minor comments (converging with the BioC version...)

- checks OK on module apps/ngs_pipeline/dev

Changes in version 3.9.2:

- removed everything related to calculateJunctionReads, junctionReads
  (due to the usage of an obsolete newCompressedList in BioC)

- checks OK on apps/ngs_pipeline/dev

Changes in version 3.9.1:

- removed everything related to SNVsOmuc, analyzeVariants,
  variantConcordance (due to gmapR conflict)

- renamed CHANGES into NEWS

Changes in version 3.9.0:

- strict copy from 3.8.0

htSeqTools
----------

Changes in version 1.5.1:

- Fixed bug in enrichedRegions method for RangedDataList objects

iPAC
----

Changes in version 1.1.2:

- Included get.Remapped.Order() function that displays the reshuffled
  amino acids after culling.

- Set the "AtomCount" column label to display "Can.Count" when the
  get.Positions function is called to be consistent with the
  get.AlignedPositions function.

- Included a "Plot.Protein.Linear" function that helps visualize the
  protein rearrangements.

- You can now specify the title of your choice to the
  "Plot.Protein.Linear" function.

IRanges
-------

Changes in version 1.18.0:

NEW FEATURES

- Add global options 'showHeadLines' and 'showTailLines' to control the
  number of head/tails lines displayed by "show" methods for Ranges,
  DataTable, and Hits objects.

- "subset" method for Vector objects now considers metadata columns.

- Add classNameForDisplay() generic and use it in all "show" methods
  defined in IRanges and GenomicRanges.

- as(x, "DataFrame") now works on *any* R object.

- Add findMatches(), an enhanced version of match() that returns all
  the matches between 'x' and 'table'. The hits are returned in a Hits
  object.  Also add countMatches() for counting the number of matches
  in 'table' for each element in 'x'.

- Add overlapsAny() as a replacement for %in% (now deprecated on
  range-based objects), and %over% and %within% as convenience wrappers
  for overlapsAny(). %over% is the replacement for %in%.

- Add 'with.mapping' arg to "reduce" methods for IRanges, Ranges,
  Views, RangesList, and CompressedIRangesList objects.

- Add "order" method for Rle objects.

- Add subsetByRanges() generic with methods for ANY, NULL, vector, and
  IRanges for now. This is work-in-progress and more methods will be
  added soon. The long term plan is to make this a replacement for
  seqselect(), but with a faster and cleaner implementation.

- Add promoters() generic with methods for Ranges, RangesList, Views,
  and CompressedIRangesList objects.

- elementLengths() now works on XVectorList objects (and thus works on
  DNAStringSet objects and family defined in the Biostrings package).
  Note that this is the first step towards having relist() work on
  XVector objects (e.g. DNAString objects) eventhough this is not ready
  yet.

- Add "mstack" method for DataFrame objects.

- Add 'name.var' argument to "stack" method for List objects for naming
  the optional column formed when the elements themselves have named
  elements.

SIGNIFICANT USER-VISIBLE CHANGES

- "distanceToNearest" methods now return a Hits instead of a DataFrame
  object.

- The behavior of distance() has changed. Adjacent and overlapping
  ranges now return a distance of 0L. See ?distance man page for
  details.  A temporary warning will be emitted by distance() until the
  release of Bioconductor 2.13.

- Change arg list of expand() generic: function(x, ...) instead of
  function(x, colnames, keepEmptyRows).

- Dramatic duplicated() and unique() speedups on CompressedAtomicList
  objects.

- Significant endoapply() speedup on XVectorList objects (this benefits
  DNAStringSet objects and family defined in the Biostrings package).

- 2x speedup to "c" method for CompressedList objects.

- classNameForDisplay() strips 'Simple' or 'Compressed', which affects
  all the "show" methods based on it. So now: > IntegerList(1:4, 2:-3)
  IntegerList of length 2 [[1]] 1 2 3 4 [[2]] 2 1 0 -1 -2 -3 instead
  of: > IntegerList(1:4, 2:-3) CompressedIntegerList of length 2 [[1]]
  1 2 3 4 [[2]] 2 1 0 -1 -2 -3

- Optimization of "[<-" method for Rle objects when no indices are
  selected (just return self).

- "stack" method for List objects now creates a factor for the optional
  name variable.

- Evaluating FilterRules now subsets by each filter individually,
  rather than subsetting by all at the end.

- Optimized which() on CompressedLogicalList objects.

- All the binary comparison operations (==, <=, etc...) on Ranges
  objects are now using compare() behind the scene. This makes them
  slightly faster and also slightly more memory efficient.

DEPRECATED AND DEFUNCT

- %in% is now deprecated on range-based objects. Please use %over%
  instead.  More precisely: - "match" and "%in%" methods that operate
  on Views, ViewsList, RangesList, or RangedData objects (20 methods in
  total) are now deprecated.  - Behavior of match() and %in% on Ranges
  objects was changed (and will issue a warning) to use equality
  instead of overlap for comparing elements between Ranges objects 'x'
  and 'table'. The old behavior is still available for match() via new
  'match.if.overlap' arg that is FALSE by default (the arg will be
  deprecated in BioC 2.13 and removed in BioC 2.14).

- tofactor() is now defunct.

- '.ignoreElementMetadata' argument of "c" method for IRanges objects
  is now defunct.

BUG FIXES

- Small fix to "unlist" method for CompressedList objects when
  'use.names' is TRUE and 'x' is a zero-length named List (the
  zero-length vector returned in that case was not named, now it is).

- "resize" method for Ranges objects now allows zero-length 'fix' when
  'x' is zero-length.

- Subsetting a Views object now subsets its metadata columns.

- Names on the vector-like columns of a DataFrame object are now
  preserved when calling DataFrame(), or when coercing to DataFrame, or
  when combining DataFrame objects with rbind().

- relist() now propagates the names on 'skeleton' when returning a
  SimpleList.

- Better argument checking in breakInChunks().

- Fix broken "showAsCell" method for ANY. Now tries to coerce
  uni-dimensional objects to vector instead of data.frame (which never
  worked anyway, due to a bug).

- Fix long standing bug in "flank" method for Ranges objects: it no
  longer returns an invalid object when NAs are passed thru the 'width'
  arg.  Now it's an error to try to do so.

- Fix issue with some of the "as.env" methods not being able to find
  the environment of the caller.

- Fix bug in "showAsCell" method for AtomicList objects: now returns
  character(0) instead of NULL on an object of length 0.

- sort() now drops NA's when 'na.last=NA' on an Rle object (consistent
  with base::sort).

- table() now handles NA's appropriately on an Rle object.

- table() now returns all the levels on a factor-Rle object.

- Fix sub-replacement of Rles when using Ranges as the index.

- Fix bug in [<- method for DataFrame objects. The fix corrects the way
  a new column created by a subset assignment is filled. Previously, if
  the first row was set, say, to '1', all values in the column were set
  to '1' when they needed to be set to NA (for consistency with
  data.frame).

- Fix bug in compare() (was not returning 0 when comparing a 0-width
  range to itself).

- Fix naming of column when passing an AsIs matrix to DataFrame() -- no
  more .X suffix.

- Fix "rbind" method for DataFrame objects when some columns are matrix
  objects.

isobar
------

Changes in version 1.5.2:

- added MSGF+ tsv import [one-line-per-psm format]

- refactored various parts of the code (proteinRatios, report-utils,
  isobar-import)

- PTM XLS report: report significance for protein ratio, and peptide
  ratio

Changes in version 1.5.1:

- added molecular weight correction to emPAI and dNSAF

- added property 'ratiodistr.class.labels': biological variability can
  be calculated in the report with other labels

- improved PDF analysis report: added number of proteins in each
  section

- added location scale family T distribution as biological variability
  distribution (distr class) and fitTlsd.

- better protein PDF analysis report layout.

Changes in version 1.5:

- Added modules for PTM validation and quantification

- Validation

- PhosphoRS XML import writers and outpout readers

- DeltaScore calculation when the data is provided

- Quantification

- All quantifications can be done now either on the protein level,
  peptide level, or modified peptide level. For modified peptide level,
  supply a matrix with a 'peptide' and 'modif' column to the
  appropriate functions.

- Correction of peptide ratios with protein ratios is possible. Also
  the variance can be adjusted (assuming no or full correlation)

- Report generation

- Import PhosphoSitePlus or neXtProt information on modification sites

KEGGREST
--------

Changes in version 1.0.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Package introduced.

NEW FEATURES

- Package introduced.

limma
-----

Changes in version 3.16.0:

- New section in User's Guide on time course experiments with many time
  points. The RNA-seq case study in User's Guide has also been revised.

- Improvements to various help pages including read.maimages.Rd,
  squeezeVar.Rd, fitFDist.Rd, trigammaInverse.Rd,
  normalizeRobustSpline.Rd, genas.Rd and roast.Rd. Previously the
  meaning of source="agilent" was mis-stated in read.maimages.Rd.

- New robust method for estimating the empirical Bayes prior, called by
  specifying robust=TRUE in the call to eBayes().  When this is TRUE
  the output df.prior is now a vector instead of a scalar.

- New function fitFDistRobustly() estimates the parameters of a scaled
  F-distribution robustly using Winsorized values.  Outlier
  observations receive smaller values for df.prior than non-outliers.
  This permits robust methods for squeezeVar(), ebayes() and eBayes(),
  all of which now have a new argument wins.tail.p to specify the tail
  proportions for Winsorizing.

- fitFDist() now permits infinite values for the covariate. It also
  gracefully handles cases where the covariate takes only a small
  number of distinct values. Similarly for eBayes() and squeezeVar()
  that call fitFDist().

- All the functions that perform gene set tests have been revised to
  make the input and output formats more consistent.
  
  roast(), mroast() and camera() are now S3 generic functions, with
  methods for EList and MAList objects.
  
  The order of arguments has been changed for roast(), mroast() and
  camera() so that the first argument is now y.
  
  All functions that perform gene sets now use the argument 'index' to
  specify which genes are included in the test set. Previously this
  argument was called 'iset' for roast() and romer() and 'indices' for
  camera().
  
  camera() and mroast() now produce a data.frames. Instead of separate
  up and down p-value columns, there is now a two-sided p-value and a
  column indicating direction of change.  There are new columns giving
  FDR values and the number of genes in each set. There is a new
  argument 'sort' to indicate whether output results should be sorted
  by p-value.
  
  mroast() has a new argument 'weights' for observational weights, to
  bring it into line with roast(),

- vennDiagram() can now plot up to five sets (previously limited to
  three).

- genas() now optionally draws a plot in which ellipses are used to
  represent the technical and biological components of correlation. It
  also now has the ability to automatically select which probes are
  used for the correlation analysis, and a new argument controls the
  method used for this selection.

- New options for the method argument of propTrueNull().

- New functions vooma() and voomaByGroup() for computing precision
  weights based on a mean-variance trend. vooma() is similar to voom()
  but for microarray data instead of RNA-Seq. voomaByGroup() allows
  different groups to have systematically different variances.

- New function predFCm() to compute predictive (shrunk) log fold
  changes.

- New function fitGammaIntercept() for estimating the intercept of a
  gamma glm with an offset. Used by genas().

- New function zscoreHyper() for computing z-score equivalents of
  deviates from a hypergeometric distribution.

- New function qqf() for qq-plots relative to an F-distribution.

- normalizeWithinArrays() with method="robustspline" now longer
  requires the layout argument to be set. The layout argument for
  normalizeRobustSpline() now defaults to a single print-tip group.

- fitFDist() now coerces degrees of freedom df1 below 1e-15 to zero.

- Due to changes in R, loessFit() no longer makes direct calls to
  foreign language code in the stats package, and instead calls R
  functions. Unfortunately, this makes loessFit() about 25-30% slower
  than previously when weights are used.

- Bug fix to read.maimages(), which was not accepting
  source="agilent.mean".

- Bug fix for contrasts.fit() when the covariance matrix of the
  coefficients (cov.coefficients) is not found in the fitted model
  object. This situation doesn't arise using any of the standard limma
  analysis pipelines.

- Bug fix to lmscFit() when the residual df = 1.

- Bug fix to readTargets() to avoid warning message when targets$Label
  is used to set row names but targets$Label contains duplicated
  entries.

MANOR
-----

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

metagenomeSeq
-------------

Changes in version 1.0.0 (2013-03-29):

- -- release!

minfi
-----

Changes in version 1.5:

- Added unit testing for the preprocessing algorithms.

- Improved the speed of SWAN for large datasets.

- Added the new class "GenomicRatioSet".  It is akin to
  "GenomicMethylSet" but instead of containing Meth and Unmeth it
  contains M and/or Beta and copy number.

- We now depend on illuminaio instead of crlmm in order to get
  readIDAT.

- Added unsrturl.bst to minimize dependences for running Sweave.

MLInterfaces
------------

Changes in version 1.39.5:

- ksvmI test/trainScores are have now rownames.  Resulted in error in
  predScores. <2013-03-12 Tue>

Changes in version 1.39.4:

- plsda prediction returns prob matrix instead of array of [, , 1] dims
  (lgatto) <2013-03-09 Sat>

Changes in version 1.39.3:

- fixed predScores (lgatto) <2013-03-01 Fri>

Changes in version 1.39.2:

- macroF1 fix ans with na.rm

Changes in version 1.39.1:

- not scaling kmeans' partition in partPlot (lgatto) <2012-11-10 Sat>

MotifDb
-------

Changes in version 1.1.3:

NEW FEATURES

- 683 new motifs derived from DGF (digital genome footprinting) from
  the stamlab (http://www.stamlab.org/) added.  Each is identified as
  either 'knownMotif' or 'novelMotif'

motifStack
----------

Changes in version 1.2.12:

NEW FEATURES

- New methods add to pcm and pfm class.

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.2.10:

NEW FEATURES

- New function motifStack

- motifCloud function is separated into two functions: motifSignature
  and motifCloud

- New class pcm and motifSig

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.2.9:

NEW FEATURES

- New layout for function motifCloud to draw motif cloud

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.2.8:

NEW FEATURES

- New function motifCloud to draw motif cloud

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.2.7:

NEW FEATURES

- Shorten the execution time

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.2.6:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- Fix bug object 'rayonWidth' not found when the pfms is not given for
  plotMotifStackWithRadialPhylog()

Changes in version 1.2.5:

NEW FEATURES

- Append a set of color prameter to control the behavior of
  plotMotifStackWithRadialPhylog()

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.2.4:

NEW FEATURES

- Append a prameter clockwise, init.angle and angle to control the
  behavior of plotMotifStackWithRadialPhylog()

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.2.3:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- fix the bug that plots x-axis repeatedly in plotMotifLogo().

Changes in version 1.2.2:

NEW FEATURES

- Append a prameter revcomp indicates whether the DNAmotifAlignment
  should use reverse-complement motif for alignment.

BUG FIXES

- Append a prameter rcpostfix, defaut is "(RC)", to the name slot of
  the pfm when DNAmotifAlignment using reverse-complement motif for
  alignment.

Changes in version 1.2.0:

NEW FEATURES

- draw aligned motif logo stacks with phylogenetic tree in radial style

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

MSnbase
-------

Changes in version 1.7.25:

- updated itraqdata to fix issue in vignette when combine(exp1, exp2)
  and different MIAPE versions <2013-03-22 Fri>

Changes in version 1.7.24:

- Mention scale in vignette <2013-03-02 Sat>

- exprsToRatio matrix method <2013-03-20 Wed>

Changes in version 1.7.23:

- new private nologging function <2013-02-21 Thu>

- adding total number of features on plotNA <2013-02-22 Fri>

- updated msnbase.r <2013-02-26 Tue>

Changes in version 1.7.22:

- msnbase.r na.rm arg <2013-02-20 Wed>

Changes in version 1.7.21:

- Added impute method <2013-02-19 Tue>

Changes in version 1.7.20:

- Explicitating that normalise and normalize are the same methods in
  the man. <2013-02-13 Wed>

Changes in version 1.7.19:

- adding MIAPE and pSet accessors: analyserDetails, analyzerDetails,
  ionSourceDetails, instrumentModel, instrumentManufacturer,
  instrumentCustomisations <2013-02-12 Tue>

- switching back to analyserDetails slot <2013-02-12 Tue>

Changes in version 1.7.18:

- readMgfData now supports comments and PEPMASS with precursor mz and
  intensity, requested by Thomas Taus <2013-02-11 Mon>

- improved and running read/writeMgfData example <2013-02-11 Mon>

- new scanIndex accessor method <2013-02-11 Mon>

Changes in version 1.7.17:

- added a analyzer accessors/slot to accomodate new mzTab files with
  more meta-data <2013-01-30 Wed>

Changes in version 1.7.16:

- fixing knitr 1.0 compatibility <2013-01-15 Tue>

Changes in version 1.7.15:

- new scale method <2013-01-11 Fri>

- renaming scale.mean and scale.median normalisation methods to
  center.mean and centre.median <2013-01-11 Fri>

Changes in version 1.7.14:

- new unexported readIspy15NData <2013-01-09 Wed>

- min.int readIspy[Silac|15N]Data arg <2013-01-11 Fri>

Changes in version 1.7.13:

- msnbase.r v0.1.1 with -h (help) arg <2013-01-08 Tue>

- msnbase.r coerce -b arg to numeric <2013-01-08 Tue>

- testing if any features left in readIspyData <2013-01-08 Tue>

Changes in version 1.7.12:

- updated makeImpuritiesMatrix to create matrix from csv file with
  correction factors <2012-12-23 Sun>

- makeImpuritiesMatrix test <2012-12-23 Sun>

- readIspyData: message instead of warning if NA in featureData
  <2012-12-24 Mon>

- Added msnbase.r script <2012-12-24 Mon>

Changes in version 1.7.11:

- new 'pattern' argument to filterNA <2012-12-15 Sat>

- vignette and man updates <2012-12-15 Sat>

- filterNA(, pattern) tests <2012-12-15 Sat>

Changes in version 1.7.10:

- new droplevels.MSnSet S3 method <2012-12-14 Fri>

- fixed errors in vignette and udpates <2012-12-14 Fri>

- vignette build stops in case of error <2012-12-14 Fri>

Changes in version 1.7.9:

- Updating processing data on readIspyData <2012-12-05 Wed>

- filterNA has a droplevels arg <2012-12-05 Wed>

- featureCV's default cv.norm is 'sum' now <2012-12-11 Tue>

- fixed featureCV for 1 sample <2012-12-11 Tue>

Changes in version 1.7.8:

- new featureCV function <2012-12-04 Tue>

- more MSnSet combineFeatures tests <2012-12-04 Tue>

- new TMT6 impurity matrix and fixed purityCorrect <2012-12-05 Wed>

- combineFeatures now automatically computes feature CVs (using
  featureCV) and collates this in featureData <2012-12-05 Wed>

- new exprsToRatios method (moved from pRoloc) <2012-12-05 Wed>

- initial implementation of impurity correction using Cramer's rule
  (see MSnbase:::cramer4) <2012-12-05 Wed>

Changes in version 1.7.7:

- added scale.mean and scale.median MSnSet normalisation method
  <2012-11-30 Fri>

- improvements to readMSData <2012-11-30 Fri>

- small updates to caching code, max level 2 <2012-11-30 Fri>

- readMSdata test <2012-12-01 Sat>

Changes in version 1.7.6:

- Fixed bug in readIspyData, reported by Claire Mulvey <2012-11-27 Tue>

Changes in version 1.7.5:

- dropping levels in readIspySilacData <2012-11-06 Tue>

- fixed plotNA <2012-11-09 Fri>

Changes in version 1.7.4:

- exporting log method <2012-11-02 Fri>

- private readIspySilacData function <2012-11-05 Mon>

- updating '['-MSnSet to log intial/final dims <2012-11-05 Mon>

Changes in version 1.7.3:

- updating readMzTabData to properly capture experiment description
  <2012-10-12 Fri>

Changes in version 1.7.2:

- fixed readMSData for MS1 <2012-10-08 Mon>

Changes in version 1.7.1:

- fixed vignettes <2012-10-02 Tue>

Changes in version 1.7.0:

- Version bump for next devel release <2012-10-01 Mon>

mzR
---

Changes in version 1.5.9:

- version bump for Rcpp 0.10.3

Changes in version 1.5.8:

- version bump for Rcpp 0.10.2

Changes in version 1.5.7:

- only load Rcpp modules after checking for Rcpp version conflict (DT)

Changes in version 1.5.6:

- Explicitely call utils::packageVersion() to avoid warning (SN)
  <2012-12-15 Sat>

Changes in version 1.5.5:

- ------------------------a o Added utils to Depends (SN) <2012-12-14
  Fri>

Changes in version 1.5.4:

- requiring Rcpp (>= 0.10.1) LG <2012-12-05 Wed>

- catching Rcpp build-time version (thanks to Dan for help!) LG
  <2012-12-05 Wed>

- checking Rcpp installed vs building versions and warn if these are
  different. LG <2012-12-05 Wed>

Changes in version 1.5.3:

- bumping version to force rebuild due to Rcpp change LG <2012-12-05
  Wed>

- added NEWS file <2012-12-05 Wed>

NarrowPeaks
-----------

Changes in version 1.3.2:

- The functionality of NarrowPeaks has been extended to multiple
  ChIP-seq sample comparison using FPCA, by implementing the function
  "narrowpeaksDiff.R".

- Vignette no. 2 has been added to the package.

- Package Title modified from "Functional Principal Component Analysis
  to Narrow Down Transcription Factor Binding Site Candidates" to
  "Analysis of Variation in ChIP-seq using Functional PCA Statistics"

ncdfFlow
--------

1.5.32: 1.add "file" argument to allow user to specify file path 2.add
        "rbind" method to allow combining more than two ncdfFlowSets
        once,

NOISeq
------

Changes in version 1.1.5 (2013-01-28):

- Some graphical improvements made

Changes in version 1.1.4 (2013-01-21):

- Fixed minor issues

Changes in version 1.1.3 (2013-01-18):

- Fixed the readData function so it can read the chromosome information
  if the chromosomes are not in numeric format.

- The NOISeq output includes now the biotype information, if provided
  to the readData function.

- A new exploratory plot for differential expression results has been
  added to the DE.plot function, in which the distribution of
  differentially expressed features across chromosomes or biotypes is
  shown.

Changes in version 1.1.2 (2012-11-30):

- Fixed some normalization issues

Changes in version 1.1.1:

- Fixed some problems with graphics

- Updated the vignette

oligo
-----

Changes in version 1.24:

BUG FIXES

- Removed dependency on RConverters.h

pathview
--------

Changes in version 1.0.0:

- Initial release with Bioconductor

- Main function: pathview

- Four functional modules: -Downloader: download.kegg; -Parser:
  node.info, combineKEGGnodes, reaction2edge; -Mapper: node.map, eg2id,
  id2eg, cpdkegg2name, cpdname2kegg, cpd2kegg, cpdidmap,
  kegg.species.code, mol.sum, sim.mol.data; -Viewer: keggview.native,
  keggview.graph, node.color, col.key, wordwrap, strfit

piano
-----

Changes in version 0.99-1:

- Added citation DOI to vignette.

- Fixed bug in geneSetSummary when no directions are available.

- Updated the man page for consensusScores, added correct output
  description.

- Fixed a bug in diffExp() regarding the result table, when gene names
  (annotation) are not available

- Fixed a bug in diffExp() so that the heatmap shows gene names if
  available, otherwise the probeset IDs

- Updated the Description field in the DESCRIPTION file.

- Removed man page for internal functions.

- Removed contrastName as output from runGSA and geneSetSummary,
  including man pages.

- Changed the man page for consensusHeatmap clarifying that the cutoff
  argument is consensus score (not rank)

- Updated the man page for loadGSC, clarifying the input.

- Reworked the vignette to fit Bioconductor, removed section on R
  introduction.

- Change name of folder for example data from exampleData to extdata,
  and updated man pages and vignette.

- Changed so that total number of gene-level statistics are printed
  during run, instead of total number of unique genes.

- Removed 'typical usages' section from man page of loadMAdata since
  this is covered in the vignette.

- Removed the arguments 'venn', 'heatmap' and 'polarPlot' from diffExp
  and replaced them with a new argument: 'plot'.

- Updated the examples for diffExp, networkPlot, consensusHeatmap and
  consensusScores to show how to handle the returned object.

- The consensusScores function now does not return its result
  invisibly.

- Added CITATION file.

- Added NEWS file.

Changes in version 0.99-0:

- Added more links to similar packages in runGSA help page.

- Updated the installation instructions in the Vignette to fit
  Bioconductor.

- Updated the loadMAdata function to use the justPlier function from
  package plier, instead of a modified version.

- Removed internal function justPlierSpec.

plgem
-----

Changes in version 1.31.2:

- new features: --added a parameter `prefix' to function `plgem.fit' so
  that multiple PLGEM fitting evaluation plots could be saved under
  different names. --added a file existence test in function
  `plgem.fit' to avoid overwriting of PLGEM fitting evaluation plot
  files. --added new parameter `gPar' to function `plgem.fit' to define
  plotting boundaries of PLGEM fitting evaluation plot. --added new
  function `setGpar' to facilitate passing graphical parameters to
  `plgem.fit'.

Changes in version 1.31.1:

- fixed bug: --corrected call to `png' in function `plgem.fit' to avoid
  partial argument match of `file' to `filename'.

- minor changes: --improved readability of progress report outputted by
  `plgem.resampledStn' when `verbose=TRUE'. --changed call to
  `packageStartupMessage' from an `.onLoad' to an `.onAttach' hook.
  --moved source of the vignette in new subfolder `vignettes'.

pRoloc
------

Changes in version 0.99.17:

- illustrating class.weights in the vignette <2013-03-24 Sun>

Changes in version 0.99.16:

- new addMarkers function <2013-03-22 Fri>

Changes in version 0.99.15:

- Fixing issues in vignette <2013-03-22 Fri>

Changes in version 0.99.14:

- Added scale in tutorial <2013-03-02 Sat>

- implemented viction <2013-03-09 Sat>

- New vignette section on phenoDisco follow up classification
  <2013-03-19 Tue>

- Using knitr engine <2013-03-19 Tue>

Changes in version 0.99.13:

- depending on MSnbase >= 1.7.23 as makeNaData needs
  MSnbase:::nologging <2013-02-27 Wed>

- updated phenoDisco parameters, added error messages when an
  insufficient number of markers per class and/or number of classes are
  used and updated phenoDisco help file <2012-02-28 Thu>

- added first belief diffscores function <2013-03-01 Fri>

Changes in version 0.99.12:

- plot2D has gained a plot argument <2013-02-19 Tue>

- new f1Count method <2013-02-20 Wed>

- new makeNaData function <2013-02-20 Wed>

- new makeNaData2 function <2013-02-21 Thu>

- new whichNAfunction <2013-02-20 Wed>

- updated tutorial vignette <2013-02-20 Wed>

- Now passing ... to predictor functions in xxxClassification (reported
  by Marianne Sandin) <2013-02-26 Tue>

Changes in version 0.99.11:

- setUnknowncol(NULL) and friedns reset to default values <2013-02-19
  Tue>

Changes in version 0.99.10:

- new default col/pch setters <2013-02-16 Sat>

Changes in version 0.99.9:

- Updates to phenoDisco: verbose param, fixed error in example, adding
  params to processingInfo <2013-02-11 Mon>

- Unexported getOtherParams method <2013-02-12 Tue>

- adding export param documentation <2013-02-12 Tue>

Changes in version 0.99.8:

- Integration of the perTurbo algorithms, contributed by Thomas Burger
  and Samuel Wieczorek <2013-01-18 Fri>

- summariseMatList now has na.rm = TRUE by default <2013-01-19 Sat>

- PerTurbo's inv/reg now as other hyperparams <2013-02-11 Mon>

Changes in version 0.99.7:

- knitr 1.0 compatibility <2013-01-15 Tue>

- Updated phenoDisco documentation and README <2013-01-10 Thu>

Changes in version 0.99.6:

- removed updateClass man <2013-01-03 Thu>

- removed old Rd files <2013-01-03 Thu>

- Deprecating *Regularisation and *Prediction function <2013-01-03 Thu>

- Updated new names in test_ml.R <2013-01-04 Fri>

Changes in version 0.99.5:

- more reg data into GenRegRes objects - cmMatrices (knn) <2012-11-15
  Thu> (other) <2012-11-30 Fri> - testPartitions (knn) <2012-11-17 Sat>
  (other) <2012-11-30 Fri>

- new minMarkers function <2012-11-16 Fri>

- renamed updateClass to minClassScore <2012-11-16 Fri>

- renamed xxxRegularisation to xxxOptimisation <2012-11-30 Fri>

- renamed xxxPrediction to xxxClassification <2012-11-30 Fri>

- renamed getRegularisedParams to getParams <2012-11-30 Fri>

- moved exprsToRatios to MSnbase <2012-12-05 Wed>

- updated phenoDisco help file <2012-12-07 Fri>

- updated phenoDisco code to cope with identical protein profiles
  <2012-12-07 Fri>

Changes in version 0.99.4:

- Adding README file describing Rd generation and suggesting roxygen2
  <2012-11-14 Wed>

- Added scol=NULL to ignore completely <2012-11-14 Wed>

Changes in version 0.99.3:

- pdres in extdata - updated vignette <2012-11-14 Wed>

- udpated pd's GS/times default <2012-11-14 Wed>

Changes in version 0.99.2:

- fixed MLearn("formula", "MSnSet", "clusteringSchema", "missing") -
  interface was wrong <2012-11-10 Sat>

- vignette updates <2012-11-10 Sat> <2012-11-11 Sun>

- nicer knn score names when scores = "all" <2012-11-11 Sun>

Changes in version 0.99.1:

- updated exprsToRatio when ncol(object) is 2 <2012-11-05 Mon>

- typos in vignette <2012-11-09 Fri> <2012-11-10 Sat>

- new MLearn method for signature c("formula", "MSnSet",
  "clusteringSchema", "missing") <2012-11-09 Fri>

- Several vignette udpates <2012-11-10 Sat>

Changes in version 0.99.0:

- Submission to Bioc <2012-11-04 Sun>

PWMEnrich
---------

Changes in version 2.0.0:

- General cleanup of the code with various small optimisations

- A FASTA file name is now also taken as input to motifEnrichment()

- The output of motifEnrichment() is now wrapped into a class
  MotifEnrichmentResults that provides a number of convenience methods
  for common tasks like ranking and plotting motifs

- Functions makeBackground() and getBackgroundFrequencies() can now
  take BSgenome objects as input. Thanks to Diego Diez for suggesting
  this and providing the code.

- Another version of motifScores() has been implemented that requires
  large amounts of memory, but is at least 2 times faster than the old
  motifScores() implementation. Use a new option
  useBigMemoryPWMEnrich() to switch to this implementation.

- PFMtoPWM now accepts a new parameter seq.count so that MotifDb motifs
  that are expressed as probabilities instead of frequencies can be
  easily used.

Changes in version 1.3.0:

- Bioconductor 2.11 release version

qpgraph
-------

Changes in version 1.16:

NEW FEATURES

- new procedures to simulate homogenous mixed graphical Markov models
  with a desired linear additive effect on the mixed linear
  associations.

- new object classes 'UGgmm', 'HMgmm' and methods 'rUGgmm()',
  'rHMgmm()', 'plot', etc. to create and simulate undirected Gaussian
  and homogeneous mixed graphical Markov models and data from them.

- new object class 'eQTLcross' and methods 'reQTLcross()', etc. to
  create and simulate expression quantitative trait loci (eQTL) models
  in experimental crosses and data from them in combination with the
  'qtl' package.

- 'qpNrr()' now also takes a qtl/cross object as input.

USER VISIBLE CHANGES

- 'qpRndHMGM()' and 'qpSampleFromHMGM()' have been deprecated in favor
  of the newer S4 classes and methods for simulation.

- new vignette on simulating molecular regulatory networks using
  qpgraph that illustrates the new collection of S4 object classes and
  methods to simulate data from graphical Markov models and from
  expression quantitative loci (eQTL) models in experimental crosses.

BUG FIXES

- correct calculation of SSD matrices and conditional independence
  tests when more than one discrete variable was involved in the test
  containing missing values using complete observations.

QuasR
-----

Changes in version 1.0.0:

INITIAL RELEASE

- ChIP-seq with support for single-end, paired-end and allele specific
  samples

INITIAL RELEASE

- RNA-seq with support for spliced alignment, single-end, paired-end
  and allele specific samples

INITIAL RELEASE

- Bis-seq with support for single-end, paired-end and allele specific
  samples

r3Cseq
------

Changes in version 1.5.0 (2012-10-25):

- added the new normalization method of 3C-seq analysis

- added the new statistical analysis for identification of 3C-seq
  interaction regions

- added the analysis for both restriction fragment and a user defined
  non-overlapping window

- added the analysis for 3C-seq data with replicates

- added the new visualization "domainograms"

- updated the existing plots

- update the exported methods

- added the options for counting reads per regions

- fixed the bug of counting reads per regions

- updated the vignette

R453Plus1Toolbox
----------------

Changes in version 1.9.2:

- Added function writeSFF to create a .sff file from a SFFContainer
  object.

- Added function qualityReportSFF to create a pdf quality report for
  SFF files.

- Fixed a bug in the function to subset SFFContainers.

rBiopaxParser
-------------

Changes in version 0.99.10:

- Changed function getInstanceProperty, this now returns all names for
  Biopaxl Level 3 instances: name, displayName and standardName by
  default. Use parameter to turn this behaivior off.

- getXrefAnnotations now correctly follows memberEntityReferences and
  memberPhysicalEntitys to obtain more annotations for physicalEntity
  instances.

- Fixed a problem with a too small buffer for generating BioPAX data in
  internal_propertyListToDF.

Changes in version 0.99.9:

- Added function removeNodes. This gracefully removes nodes from a
  regulatory graph: Connecting parents and children by their multiplied
  edge weights.

- Added function combineNodes. This gracefully combines nodes from a
  regulatory graph. This is basically a wrapper for
  graph::combineNodes(nodes, graph, newName, collapseFunction=max).

Changes in version 0.99.8:

- Speed up of many function which were beautifully programmed but
  became increasingly slow when using the whole Reactome database.

- Fixed a nasty bug where a complex without any pyhsical entities could
  cause many functions to fail badly.

- Fixed CITATION file.

Changes in version 0.99.7:

- Fixed getXrefAnnotations. This function now works with Biopax Level 3
  and can retrieve all Xref annotations of any given vector of instance
  ids.

- Added CITATION file! rBiopaxParser was pubished in Bioinformatics!

- rBiopaxParser - an R package to parse, modify and visualize BioPAX
  data. Kramer F, Bayerlova M, Klemm F, Bleckmann A, Beissbarth T.
  Bioinformatics (2013) 29(4): 520-522.

- rBiopaxParser has been accepted to Bioconductor!
  http://www.bioconductor.org/packages/devel/bioc/html/rBiopaxParser.html

- README.md has been updated accordingly with steps to install directly
  from Bioconductor.

- Function selectInstances correctly returns all referenced instances
  with parameter "includeReferencedInstances=TRUE" now.

Changes in version 0.99.6:

- Fixed pathway layouting. It broke simarily to:
  http://permalink.gmane.org/gmane.science.biology.informatics.conductor/44745

- Function createBiopax can create Biopax Level 3 models now!

- Function listInstances was very slow. It is alot faster now with huge
  databases like Reactome.

- Functions pathway2Geneset, pathway2RegulatoryGraph,
  pathway2AdjacencyMatrix now fully support nested Subpathways and
  PathwayOrder/PathwaySteps in Biopax Level 2 and Level 3.

- Started to integrate RUnit tests for some of my functions. More to
  come.

Changes in version 0.99.5:

- Submitted some fixes. colorGraphNodes is finally working (again)
  using graph and Rgraphviz packages.

Changes in version 0.99.3:

- We are finally supporting BioPAX level 3 now! All functions (except
  for getXrefAnnotations) have been updated and work correctly with
  BioPAX Level 2 and Level 3 objects.

Changes in version 0.99.0:

- Pushing version number to 0.99.0 according to Bioconductor checklist

Changes in version 0.23:

- Fixes to documentation. All function documentation is featuring
  examples now.

Changes in version 0.21:

- Attention: Major function renaming and restructuring in this version!
  This unfortunatly leads to incompatibilities with the previous
  versions due to the renaming.

- The variables "instancetype" and "instanceid" have been renamed to
  "class" and "id". This keeps everyones code shorter and matches
  general ontology syntax.

- Function names of the ever-growing number of convenience functions
  were a mess and hard to understand. I have taken the time to refactor
  and consolidate the collection of functions to increase usability and
  give function(name)s a common structure. This will also serve as a
  stepping stone to allow biopax level 3 integration.

- Functions with names starting "get" were previously being used for
  all selecting/accessing of data. These functions are now effectively
  being split up into functions with names "select" and "list".

- Functions with prefix "select" return a TRUE/FALSE index vector for
  the whole biopax internal data.frame. This allows easy combination of
  select statements via logical & and |.

- Functions with prefix "list" return a specific list of IDs or items
  depending on the function definition.

- Functions with prefix "get" return newly created objects or
  informations about the biopax data.frame or its instance(s).

Changes in version 0.15:

- Added function createBiopax which creates a new Biopax model from
  scratch.

- Added some more parameter checks.

- Added verbosity to writeBiopax since re-exporting the NCI Pathway
  Interaction Database took almost 15 minutes.

- Updated the vignette. Nice pictures, text and more!

Changes in version 0.14:

- Some more sanity checks added to fix sloppy BioPAX models

- Added functions getBiopaxIDsByName, which does as the name suggests,
  and getXrefAnnotations, which returns all Xref annontations of given
  Biopax IDs.

- Function pathway2RegulatoryGraph now has the optional parameter
  useIDasNodenames, this causes a regulatory graph to be created with
  the Biopax IDs used instead of the instance names.

- Function plotRegulatoryGraph can be passed logical argument
  layoutGraph, which can be used to stop the function from re-layouting
  the graph in case you modified it.

- Added a check for BioPAX Levels different from Level 2 and started to
  build a foundation to enable those too.

- Unzipping functionality of downloadBiopaxData only works using Linux.
  Added checks and instructions for Windows users.

Changes in version 0.13:

- Added functions removeInstance, removeProperties, addPhysicalEntity,
  addPhysicalEntityParticipant, addBiochemicalReaction, addControl.

- See new examples for modifying a biopax model with R in the short
  vignette!

Changes in version 0.12:

- Fixed some documentation and vignette pdfs.

Changes in version 0.11:

- Set packages RCurl and Rgraphviz from imports to suggests. This means
  you will not have to have these installed to run most of the
  commands. However layouting the graphs requires Rgraphviz and
  downloading data requires RCurl.

Changes in version 0.10:

- Detailed installation instructions are available now. rBiopaxParser
  depends on some packages which might be a hassle to install. Check
  out the README or github front page for instructions!

- Provided some more documentation.

- Fixed a very embarrassing typo in downloadBiopaxData.

Changes in version 0.09:

- Provided some more documentation.

- getNeighborhood: a function to build a neigborhood graph for
  molecules

Changes in version 0.08:

- First version on GitHub! Code still needs some cleaning up but
  everything works, longer examples comming next!

Changes in version 0.06:

- Visualization & Layouting added. Check out new functions

- layoutRegulatoryBiopaxGraph, plotRegulatoryBiopaxGraph

- intersectGraphs, uniteGraphs, diffGraphs

Changes in version 0.04:

- Building regulatory graphs from biopax model works now!

- pathway2adjacancyMatrix, pathway2RegulatoryGraph,

Changes in version 0.01:

- Parsing works! Alot of thing changing quickly now, still bugs to fix.

Rbowtie
-------

Changes in version 1.0.0:

INITIAL RELEASE

- bowtie version 0.12.8

- SpliceMap version 3.3.5.2 modified

BUG FIXES

- Several bug fixes and modifications in the source code of SpliceMap
  for integration in Rbowtie/QuasR.

RCytoscape
----------

Changes in version 1.9.8:

BUG FIXES

- Significant (3x) speedup.  A 5000-node, 6000-edge graph transmits to
  Cytoscape from R in about 20 seconds.

ReactomePA
----------

Changes in version 1.3.1:

- bug fixed of ALLEXTID. <2013-03-1, Fri>

ReportingTools
--------------

Changes in version 2013-4-1:

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

ReQON
-----

Changes in version 1.5.0:

BUG FIXES

- Some users were receiving errors from hist.default(), which depended
  on the training region specified.  This has been fixed.

Rgraphviz
---------

Changes in version 2.3:

- Fixing bug in Graphviz on OpenBSD regarding sincos (Thanks to Rainer
  Hurling).

- Now using path.expand to fix tilde-expansion.

- We no longer include <R_ext/RConverters.h>; which seemingly was not
  needed (Thanks to Brian D. Ripley).

- Fix handling of fontsize attrib for nodes and edges (seems it has
  been broken at least since 2006).

rhdf5
-----

Changes in version 2.4.0:

NEW FEATURES

- support for reading 64-bit integers added

- support for reading variable length strings added

- support for reading scalar objects added

USER VISIBLE CHANGES

- NEWS.Rd added

- display of chunksize.pdf as a vignette avoided

Risa
----

1.1.0: 1. The ISAtab-class was extended to include the definition of
        assay.names, factors, treatments and groups. 2. The previous
        processAssayXcmsSet method was split into
        processAssayXcmsSet.1factor and processAssayXcmsSet to consider
        the first factor only (as it was in the previous definition) or
        all the factors, respectively. 3. More methods for processing
        mass spectrometry assays were added: - the method
        getMSAssayFilenames retrieves a list with the mass spectrometry
        assay filenames - the method getRawDataFiles retrieves a list
        with all the files listed under the column 'Raw Spectral Data
        File' within the mass spectrometry assays. - the method is.ms
        receives an object from the ISAtab-class and an assay filename
        as parameter, and retrieves TRUE if the assay filename is a
        mass spectrometry assay, or FALSE otherwise. 4. A method called
        suggestBiocPackage was added, which retrieves a list of
        packages in Bioconductor that might be relevant to the ISAtab
        dataset, according to its assays' mesaurement and technology
        types. This method relies on the BiocViews annotations for each
        of the packages available in Bioconductor in a specific
        version. 5. Added AssayTab-class and subclasses for MS,
        Microarray, Seq and NMR. 6. Defined method
        getAssayRawDataFilenames and getRawDataFilenames 7. Defined
        method getExpressionSet for microarray assays relying on affy
        package.

RNASeqPower
-----------

Changes in version 0.99-2:

- Add the tests directory

- Repair a mistake, revealed in the tests, found previously, but not
  propogated between two copies of the code on different machines.

- Add this NEWS file.

rols
----

Changes in version 1.1.6:

- Using knitr as VignetteEngine and visual tweaking <2013-02-17 Sun>

Changes in version 1.1.5:

- update vignette for knitr 1.0 compatibility, based on Dan's updates
  in stable version <2013-01-15 Tue>

Changes in version 1.1.4:

- changed query to avoid failure (see issue #6) <2012-12-26 Wed>

- added missing space in `Empty query results after 3attempts`
  <2012-12-26 Wed>

Changes in version 1.1.3:

- updated olsQuery example <2012-12-25 Tue>

Changes in version 1.1.2:

- .rolsEnv now has emptyenv() as parent <2012-10-31 Wed>

Changes in version 1.1.1:

- fixing vignette <2012-10-02 Tue>

- new iface <2012-10-02 Tue>

Changes in version 1.1.0:

- version bump for next devel <2012-10-01 Mon>

RPA
---

Changes in version 1.15.34 (2013-03-07):

- rpa.online has been added and benchmarked

- vignette contents moved to online wiki

rqubic
------

Changes in version 1.8.1 (2011-08-02):

- Fix errors in documentations

Changes in version 1.8.0 (2011-07-11):

- Submission to Bioconductor

Rsamtools
---------

Changes in version 1.12.0:

NEW FEATURES

- BamSampler draws a random sample from BAM file records, obeying any
  restriction by ScanBamParam().

- Add argument 'obeyQname' to BamFile. Used with qname-sorted Bam files
  only.

- Add readBamGAlignmentsList function for reading qname-sorted Bam
  files into a GAlignmentsList object.

USER-VISIBLE CHANGES

- bamPath and bamIndicies applied to BamViews returns named vectors.

- 'yieldSize' argument in BamFile represents the number of unique
  qnames when 'obeyQname=TRUE'.

BUG FIXES

- completely free razip, bgzip files when done.

- sortBam, indexBam fail gracefully on non-BAM input.

- headerTabix on an open TabixFile no longer reads the first record

- scanBcfHeader provides informative error message when header line
  ('#CHROM POS ...') is missing

Rsubread
--------

Changes in version 1.10.0:

NEW FEATURES

- Rsubread package can now run on Mac OS X operating systems.

NEW FEATURES

- Major updates to the function 'featureCounts', a general-purpose read
  summarization function.

RTCA
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

rtracklayer
-----------

Changes in version 1.20:

NEW FEATURES

- Table query interface supports multiple query ranges.

- Files (RTLFile objects) can be directly uploaded to UCSC, via
  track<-.

SIGNIFICANT USER-VISIBLE CHANGES

- All methods with asRangedData argument now have it default to FALSE,
  instead of TRUE. A warning is issued if the argument is missing.
  Eventually, we will drop all support for RangedData import (export
  will still work via automatic coercion to GRanges).

BUG FIXES

- Chromosome list for a genome is now downloaded from the table
  browser, instead of the Genome Browser page. This supports genomes
  with more than 1000 contigs.

- BEDX+Y formats now work when a track line is present, and the
  extraCols argument is used for the column names.

- path.expand() is now called for paths passed off to the Kent library.

- Order of metadata columns upon GFF import no longer depends on
  LC_COLLATE.

SeqArray
--------

Changes in version 0.99.1:

- first submission

ShortRead
---------

Changes in version 1.17:

SIGNIFICANT USER-VISIBLE CHANGES

- FastqSampler can return records in the order encountered in the
  sampled file.

- Increase to 10000 the number of reads examined for determining Fastq
  quality type

- as(FastqQuality, "numeric") returns a vector of quality scores
  concatenated end to end (previously cycle to cycle), without padding
  to effective equal width

BUG FIXES

- trimTails, succesive=TRUE would return inconsistent results

- FastqStreamer, FastqSampler parse fastq files created with '\r'

SSPA
----

Changes in version 1.99.0:

NEW FEATURES

- Added new estimation method for the density of effect sizes and
  proportion of non differentially expressed genes

- power and sample size analysis can now also be performed for
  experimental designs that lead to F and \chi^2 statistics e.g. RNAseq
  data

BUG FIXES

- suppressWarings from dt or pt when ncp is large

synapter
--------

Changes in version 1.1.5:

- updated references <2013-03-22 Fri>

Changes in version 1.1.4:

- added citation <2013-03-21 Thu>

- vignette uses knitr engine and scrartcl class <2013-03-21 Thu>

Changes in version 1.1.3:

- knitr 1.0 compatibility, based on Dan's updates <2013-01-15 Tue>

Changes in version 1.1.2:

- added note about tcl BWidget package in synapterGUI man <2012-10-04
  Thu>

Changes in version 1.1.1:

- fixing vignette <2012-10-02 Tue>

Changes in version 1.1.0:

- new devel version bump <2012-10-01 Mon>

TargetSearch
------------

Changes in version 1.16.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Function 'fixRIcorrection' is now defunct (previously was
  deprecated). Use 'fixRI'.

BUG FIXES

- Removed references to deprecated R functions.

- Fixed bug in 'ImportLibrary.msp' that occurs if there is only one
  metabolite.

TEQC
----

Changes in version 2.7.2:

- bug fix in 'coverage.target' (chromosomes without any read but
  appearing in the BAM file header were a problem)

TransView
---------

Changes in version 1.2.11:

BUG FIXES

- Fix compilation problem with new Rsamtools version on windows.

Changes in version 1.2.9:

BUG FIXES

- Updated vignette

Changes in version 1.2.8:

NEW FEATURES

- meltPeaks() normalizes to Reads Per Million mapped reads after
  quality filtering according to the filtered_reads slot.

NEW FEATURES

- parseReads() now also returns the total and local base pairs covered
  which can be accessed by the slots gsize() and lsize() respectively.

Changes in version 1.2.7:

NEW FEATURES

- New convenience function meltPeaks(), which returns a data frame with
  normalized peak densities suitable for plotting with ggplot2.

BUG FIXES

- Minor help file corrections.

Changes in version 1.2.6:

BUG FIXES

- macs2gr() supports relative distances instead of absolute distances
  reported by older MACS versions.

Changes in version 1.2.5:

BUG FIXES

- annotatePeaks() did not accurately resolve ambiguities in some
  instances.

Changes in version 1.2.4:

BUG FIXES

- gtf2gr() can now handle gene and transcript IDs with white spaces.

Changes in version 1.2.3:

NEW FEATURES

- All slice methods can now return binned densities of predefined
  width.

NEW FEATURES

- The method to bin data can be specified by a new option in plotTV and
  the slice methods called bin_method. The amount of bins to be
  returned can be specified by nbins.

NEW FEATURES

- Default method to plot densities in plotTV is now mean instead of
  linear interpolation using approx. median or max are additional
  possibilities.

Changes in version 1.2.2:

BUG FIXES

- Minor bug fix in plotTV

Changes in version 1.2.1:

BUG FIXES

- Minor bug fix in plotTVData

Changes in version 1.2.0:

NEW FEATURES

- A new class TVResults is returned by plotTV containing all important
  results of the clustering and plotting.

- TVResults objects can be accessed by their slots and plotTVData.
  plotTVData returns a data frame with the summarized results.

- Added option showPlot to plotTV() to suppress plotting optionally.
  This is useful if only the clustering results are needed.

- Added an option name_width to plotTV() to enable customized widths of
  the row labels.

BUG FIXES

- Improved kmeans clustering

TurboNorm
---------

Changes in version 1.7.2:

BUG FIXES

- normalize.AffyBatch.pspline now correcly handles the weights

Changes in version 1.7.1:

NEW FEATURES

- faster turbotrend due to c implementation for the summation and
  couting by Paul Eilers

- add robustifying iterations using same approach as lowess code by
  Paul Eilers

BUG FIXES

- print array number when normalization of an array causes an error

VariantAnnotation
-----------------

Changes in version 1.6.0:

NEW FEATURES

- VCF is now VIRTUAL. Concrete subclasses are CollapsedVCF and
  ExpandedVCF.

- Add filterVcf() generic and methods for character and TabixFile.
  This method creates one VCF file from another, using FilterRules.

- Enhance show,VCF method with header information.

- Stephanie Gogarten added genotypeToSnpMatrix() generic and
  CollapsedVCF and matrix methods.

- Chris Wallace added snpSummary() generic and CollapsedVCF method.

- Add cbind and rbind for VCF objects.

MODIFICATIONS

- writeVcf,connection-method allows writing to console and appending.

- writeVcf,connection-method accepts connections with open="a", only
  adding a header if the file does not already exist.

- predictCoding and genotypeToSnpMatrix can now handle ALT as
  CharacterList. Structural variants are set to empty character ("").

- When no INFO data are present in a vcf file, the info() slot is now
  an empty DataFrame. Previously an empty column named 'INFO' was
  returned.

- Empty VCF class now has an empty VCFHeader

- expand,CollapsedVCF-method expands 'geno' data with Number=A.

- VCF class accessors "fixed", "info" now return DataFrame instead of
  GRanges. "rowData" returns fixed fileds as the mcols.

- Updates to the vignette.

DEPRECATED and DEFUNCT

- Deprecate dbSNPFilter() and regionFilter().

- Deprecate MatrixToSnpMatrix().

BUG FIXES

- Multiple bugs fixed in "locateVariants".

- Multiple bugs fixed in "writeVcf".

- Bug fixed in subsetting of VCF objects.

- Bug fixed in "predictCoding" related to QUERYID column not mapping
  back to original indices (rows).

VariantTools
------------

Changes in version 1.2.0:

NEW FEATURES

- Tally, call and export indels (using same algorithm as for SNVs).

- Add post-filter that discards variants that are clumped together on
  the chromosome (likely mapping errors).

- Add filter for masking regions like simple / low complexity repeats.

- Add a filter that performs a t-test on the alt vs. ref read
  positions.

- Add callWidtype() function for determining whether a position is
  variant, wildtype or uncallable, assuming the built-in variant
  calling filters. This is based on a power calculation that considers
  the coverage.

- Some functions for estimating concordance between samples have been
  added; these were developed for sample ID verification and should be
  considered experimental.

- matchVariants() utility for matching variants by pos and alt.

USER-VISIBLE CHANGES

- The VCF output is now always in expanded form (one alt per row). The
  AD (allele depth) geno tag contains the REF and ALT counts, while AP
  (allele present) indicates presence of the REF and/or ALT allele.
  Besides the DP tag, all other tags were removed. These changes bring
  VariantTools more in line with GATK.

- Control alt and total counts are returned from
  callSampleSpecificVariants.

BUG FIXES

- The power cutoff in the sample-specific algorithm was not considering
  the minimum alt read count filter.

wateRmelon
----------

0.99.16: seabird() no longer takes a vector for the stop argument

0.99.15: temporarily removed dependency on ROCR (and thus gdata)

0.99.14: normalized intensities as well as betas fixed a bug in
        as.methylumi

0.99.13:

0.99.12: Sentrix IDs should work (thanks to Elodie Portales-Casamar for
        bug report).

0.99.11:

0.99.10:

0.9.8: new MethyLumiSet method for pfilter RGtoy2 removed from
        data(minfitoy).  This is because minfi now expects the manifest
        to be a package.  swan function patched so that it can process
        subsets of array data (provided by Jovana Maksimovic).  Now
        works with melon data set.  as.methylumi generic function
        added. Gets fData from package:IlluminaHumanMethylation450k.db
        . The MethyLumiSet method is useful for ensuring methylumi
        objects contain this data.  Andrew Teschendorff's BMIQ function
        and a MethyLumiSet method added minfi 1.4.0 has a changed
        manifest structure resulting in SNP probes being left out of
        the MethylSet.  This breaks genki() for the minfi objects and
        also the minfitoy data() set.

xcms
----

Changes in version 1.35.7:

BUG FIXES

- fixed indexing bug in group.nearest, which under certain
  circumstances caused all peaks in the first sample to be ignored
  (reported by Tony Larson)

Changes in version 1.35.6:

BUG FIXES

- Obiwarp retention time alignment error-ed if scanrange was used as a
  parameter setting during xcmsSet/peak detection The method now tries
  to automatically find the set scanrange and uses this range for
  alignment.

Changes in version 1.35.4:

NEW FEATURES

- Introducing parallel fillPeaks

USER VISIBLE CHANGES

- Replace snow requirement with minimum R version 2.14.0

Changes in version 1.35.3:

BUG FIXES

- if group.density was used with very low minfrac settings (< 0.5) it
  did not return all feature groups, but only those that include
  features from at least 50% of samples in a group. This limitation was
  removed.

Changes in version 1.35.2:

UPDATED FEATURES

- Behind the scenes xcms now uses the xcmsSource class to read raw
  data.  This allows e.g. to write a class that pulls raw data from
  e.g. a database

BUG FIXES

- massifquant: simplified logic structure of Tracker::claimDataIdx
  resolved failure on new test case.

- massifquant: reporting features data structure compatible with
  multiple sample comparison within XCMS.

Changes in version 1.35.1:

UPDATED FEATURES

- The mzData export is now much faster and uses less memory

xps
---

Changes in version 2.16:

VERSION xps-1.19.10

- update XPSSchemes.cxx to correct possible memory error in TString
  array names[i]

VERSION xps-1.19.9

- update README

VERSION xps-1.19.8

- replace .path.package() with path.package()

VERSION xps-1.19.7

- update XPSSchemes.cxx to replace error with warning for missing
  annotation header '%genome-species'

VERSION xps-1.19.2 - 6

- update configure.win and Makefile.win for ROOT compiled with MinGW;
  xps works now on Windows 7

VERSION xps-1.19.1

- update Import..Annotation() methods in XPSchemes.cxx to protect
  against tabs in Affymetrix annotation files

- update script4schemes.R to include schemes with annotation na33

Packages removed from the release
=================================

The following packages are no longer in the release:
cosmo, cosmoGUI, gene2pathway
