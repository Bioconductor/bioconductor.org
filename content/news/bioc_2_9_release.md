Bioconductors:

We are pleased to announce Bioconductor 2.9, consisting of 517
software packages and more than 500 up-to-date annotation packages.
There are 60 new software packages, and many updates and improvements
to existing packages; 10 packages have been removed from this
release. Bioconductor 2.9 is compatible with R 2.14.0, and is
supported on Linux, 32- and 64-bit Windows, and Mac OS.  This release
includes an updated Bioconductor [Amazon Machine Image][ami].  Visit
[http://bioconductor.org][bioc] for details and downloads.

Contents
--------

* Getting Started with Bioconductor 2.9
* New Software Packages
* NEWS from new and existing packages
* Packages removed from the release

Getting Started with Bioconductor 2.9
=====================================

To install Bioconductor 2.9:

1. Install R 2.14.0.  Bioconductor 2.9 has been designed expressly for
this version of R.

2. Follow the instructions at [http://bioconductor.org/install/][install].

New Software Packages
=====================

There are 60 new packages in this release of Bioconductor.

- AGDEX: Agreement of Differential Expression Analysis

- BiocInstaller: Install/Update Bioconductor and CRAN Packages

- biovizBase: Basic graphic utilities for visualization of genomic
  data.

- CellNOptR: R version of CellNOpt, boolean features only

- CNAnorm: A normalization method for Copy Number Aberration in
  cancer samples.

- cn.mops: Mixture of Poisson for CNV detection in NGS data

- Cormotif: Correlation Motif Fit

- cqn: Conditional quantile normalization

- cummeRbund: cummeRbund: The finishing touch on your Tuxedo
  workflow. Analysis, exploration, manipulation, and visualization of
  Cufflinks HTS data.

- DECIPHER: Database Enabled Code for Ideal Probe Hybridization
  Employing R

- DEXSeq: Inference of differential exon usage in RNA-Seq

- DiffBind: Differential Binding Analysis of ChIP-Seq peak data

- dks: The double Kolmogorov-Smirnov package for evaluating multiple
  testing procedures.

- DOSE: Disease Ontology Semantic and Enrichment analysis

- DTA: Dynamic Transcriptome Analysis

- EDASeq: Exploratory Data Analysis and Normalization for RNA-Seq

- exomeCopy: Detection of CNV in exome/targeted sequencing data

- fastseg: fastseg - a fast segmentation algorithm

- flowType: Phenotyping Flow Cytometry Assays

- flowWorkspace: Import flowJo Workspaces into BioConductor and
  replicate flowJo gating with flowCore

- GeneExpressionSignature: Gene Expression Signature based Similarity
  Metric

- ggbio: Static visualization for genomic data.

- GOFunction: GO-function: deriving biologcially relevant functions
  from statistically significant functions

- graphite: GRAPH Interaction from pathway Topological Environment

- GRENITS: Gene Regulatory Network Inference Using Time Series

- GWASTools: GWASTools: Tools for Genome Wide Association Studies

- htSeqTools: Quality Control, Visualization and Processing for
  High-Throughput Sequencing data

- IdMappingRetrieval: ID Mapping Data Retrieval

- inSilicoDb: Access to the InSilico Database

- iontree: Data management and analysis of ion trees from ion-trap
  mass spectrometry

- isobar: Analysis and quantitation of isobarically tagged MSMS
  proteomics data

- minfi: Analyze Illumina's 450k methylation arrays

- MmPalateMiRNA: Murine Palate miRNA Expression Analysis

- mzR: parser for netCDF, mzXML, mzData and mzML files (mass
  spectrometry data)

- ncdfFlow: ncdfFlow: A package that provides ncdf based storage
  based flow cytometry data.

- NormqPCR: Functions for normalisation of RT-qPCR data

- nucleR: Nucleosome positioning package for R

- PAN: Posterior association networks and functional modules inferred
  from rich phenotypes of gene perturbations

- PREDA: Position RElated Data Anlysis

- predictionet: Inference for predictive networks designed for (but
  not limited to) genomic data

- qtbase: Interface between R and Qt

- qtpaint: Qt-based painting infrastructure

- r3Cseq: Analysis of Chromosome Conformation Capture and
  Next-generation Sequencing (3C-seq)

- RamiGO: AmiGO visualize R interface

- randPack: Randomization routines for Clinical Trials

- RCASPAR: A package for survival time prediction based on a
  piecewise baseline hazard Cox regression model.

- ReadqPCR: Read qPCR data

- RedeR: RedeR: bridging the gap between hierarchical network
  representation and functional analysis.

- REDseq: Analysis of high-throughput sequencing data processed by
  restriction enzyme digestion

- Repitools: Epigenomic tools

- ReQON: Recalibrating Quality Of Nucleotides

- rqubic: Qualitative biclustering algorithm for expression data
  analysis in R

- RTopper: This package is designed to perform Gene Set Analysis
  across multiple genomic platforms

- stepwiseCM: Stepwise Classification of Cancer Samples using
  Clinical and Molecular Data

- Streamer: Enabling stream processing of large files

- sva: Surrogate Variable Analysis

- TSSi: Transcription Start Site Identification

- tweeDEseq: RNA-seq data analysis using the Poisson-Tweedie family
  of distributions

- VariantAnnotation: Annotation of Genetic Variants

- zlibbioc: An R packaged zlib-1.2.5

NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages. The following package NEWS is available:

AnnotationDbi
-------------

Changes in version 1.16:

NEW FEATURES

- New Object types: Org.Db, ChipDb, GODb.  These are loaded when an
Annotation package is loaded.

- select method for these new objects to extract data.frames of
available annotations.  Users can specify keys, along with the
keytype, and the columns of data that they want extracted from the
annotation package.

- keys now will operate on new objects to expose most ID types as keys

- keytypes will show which kinds of IDs can be used as a key by select

- cols will display the kinds of data that can be extracted by select

- Users who need organism packages for organisms that are annotated at
NCBI can now use makeOrgPackageFromNCBI() to generate a basic
organism package on the fly

SIGNIFICANT USER-VISIBLE CHANGES

- New objects are always named after the package that loaded them. So
if you just loaded the hgu95av2.db package, your ChipDb object will
be called hgu95av2.db

- new fields are now required in the metadata table of annotation
packages to allow the on-the-fly creation of these new objects. If
you are constructing a custom package and want to extend this
infrastructure, you may need to add some fields to the metadata
table.

BUG FIXES

- makeProbePackage, and other generated Annotation packages have had
templates modified so that packages can load and also to enable
loading of new objects to allow new features in this release.

- updated GOFrame constructor to allow use of new types of evidence
codes.

- generated packages will create databases with the appropiate metadata
entries for the new objects.

clusterProfiler
---------------

Changes in version 1.1.22:

- update vignette. <2011-07-23 Thu>

Changes in version 1.1.21:

- update roxygen docs. <2011-07-23 Sat>

Changes in version 1.1.20:

- update roxygen docs. <2011-07-19 Wed>

Changes in version 1.1.19:

- change to using roxygen for generating Rd file. <2011-07-19 Tue>

Changes in version 1.1.18:

- clean up some codes. <2011-07-19 Tue>

Changes in version 1.1.17:

- Remove dependence of GOstats and Category. <2011-07-18 Mon>

Changes in version 1.1.16:

- Re-implement function enrichGO. <2011-07-18 Mon>

- add parameter readable, for mapping gene IDs to gene names.
<2011-07-18 Mon>

Changes in version 1.1.15:

- Modified plot method of groupGOResult, with colour coded. <2011-06-15
Wed>

DESeq
-----

Changes in version 1.6.0:

- changes to dispersion estimation scheme to improve handling of
outliers (see vignette for details)

- refactorized the workflow and overhauled and expanded the vignette
(see there for changes to user interface)

- various changes and bug fixes (see appendix in vignette and NEWS
file)

Changes in version 1.5.29 (2011-10-10):

- fixed a bug in parametricDispersionFit. (Thanks to Alejandro Reyes
for spotting and fixing this one.)

Changes in version 1.5.27 (2011-10-18):

- method="pooled" is now default in estimateDispersions

- added 'fitInfo' accessor function

Changes in version 1.5.24 (2011-08-10):

- fitNbinomGLM no longer stops when glm.fit throws an error

Changes in version 1.5.23 (2011-07-14):

- fixed bug in 'nbinomTest' again; now using a symmetric two-tailed
test (i.e., the one-tailed p value is doubled)

- added warning in case of method='blind' if not used together with
sharingMode='fit-only'

Changes in version 1.5.22 (2011-07-14):

- fixed bug in 'nbinomTest' that led to wrong p values for dispersions
larger than 1

Changes in version 1.5.21 (2011-07-14):

- minor changes to vignette

Changes in version 1.5.20 (2011-07-07):

- fixed bug in 'getVarianceStabilizedData'

- more informative error message if parametric fit fails

Changes in version 1.5.1:

refactorized the way how dispersion estimation is handled, overhauled the vignette

- vignette now uses data from 'pasilla' package as example

- GLMs explained in vignette

- 'estimateVarianceFunctions' renamed to 'estimateDispersions'

- estimateDispersions now stores intermediate results in 'fitInfo' slot

- new functionality for 'estimateDispersion': 'sharingMode' option with
new sharing modes 'maximum' and 'gene-est-only'; 'fitType' with new
fitting type 'parametric'; changed defaults for these; new method
'pooled-CR' to do edgeR-style CR estimation etc.

- several diagnostic functions removed

- lots of further small changes

DEXSeq
------

Changes in version 2011-10-03:

- Changes to the documentation and to the vignette.  New functions
added and more support for gene names and strange exonids.

Changes in version 2011-07-12:

- Parallelization possibility added to the code, with its description
in the vignette and single exon genes ignored properly.

Changes in version 2011-07-01:

- 'estimateSizeFactors' and 'estimatedispersions' function added as S4
methods, no more problems when loading DESeq and DEXSeq in the same R
session

DOSE
----

Changes in version 0.99.7:

- fixed bug in .SemVal_internal

Changes in version 0.99.6:

- Use qvalue instead of fdrtool to calculate qvalues.

edgeR
-----

Changes in version 2.4.0:

NEW FEATURES

- New function spliceVariants for detecting alternative exon usage from
exon-level count data.

- A choice of rejection regions is now implemented for exactTest, and
the default is changed from one based on small probabilities to one
based on doubling the smaller of the tail probabilities. This gives
better results than the original conditional test when the dispersion
is large (especially > 1).  A Beta distribution approximation to the
tail probability is also implemented when the counts are large,
allowing exactTest() to be much faster and use less memory.

- exactTest now uses tagwise.dispersion by default if found in the
object.

- estimateCRDisp is removed. It is now replaced by
estimateGLMCommonDisp, estimateGLMTrendedDisp and
estimateGLMTagwiseDisp.

- Changes to glmFit so that it automatically detects dispersion
estimates if in data object. It uses tagwise if available, then
trended, then common.

- Add getPriorN() to calculate the weight given to the common parameter
likelihood in order to smooth (or stabilize) the dispersion
estimates. Used as default for estimateTagwiseDisp and
estimateGLMTagwiseDisp.

- New function cutWithMinN used in binning methods.

- glmFit now S3 generic function, and glmFit has new method argument
specifying fitting algorithm.

- DGEGLM objects now subsettable.

- plotMDS.dge is retired, instead a DGEList method is now defined for
plotMDS in the limma package.  One advantage is that the plot can be
repeated with different graphical parameters without recomputing the
distances.  The MDS method is also now much faster.

- Add as.data.frame method for TopTags objects.

- New function cpm to calculate counts per million conveniently.

- Adding args to dispCoxReidInterpolateTagwise to give more access to
tuning parameters.

- estimateGLMTagwiseDisp now uses trended.dispersion by default if
trended.dispersion is found.

- Change to glmLRT to ensure character coefficient argument will work.

- Change to maPlot so that any really extreme logFCs are brought back
to a more reasonable scale.

- estimateGLMCommonDisp now returns NA when there are no residual df
rather than returning dispersion of zero.

- The trend computation of the local common likelihood in
dispCoxReidInterpolateTagwise is now based on moving averages rather
than lowess.

- Changes to binGLMDispersion to allow trended dispersion for data sets
with small numbers of genes, but with extra warnings.

BUG FIXES

- dispDeviance and dispPearson now give graceful estimates and messages
when the dispersion is outside the specified interval.

- Bug fix to mglmOneWay, which was confusing parametrizations when the
design matrix included negative values.

- mglmOneWay (and hence glmFit) no longer produces NA coefficients when
some of the fitted values were exactly zero.

- Changes to offset behaviour in estimateGLMCommonDisp,
estimateGLMTrendedDisp and estimateGLMTagwiseDisp to fix bug. Changes
to several other functions on the way to fixing bugs when computing
dispersions in data sets with genes that have all zero counts.

- Bug fix to mglmSimple with matrix offset.

- Bug fix to adjustedProfLik when there are fitted values exactly at
zero for one or more groups.

fabia
-----

Changes in version 2.0.0:

NEW FEATURES

- spfabia: fabia for a sparse data matrix (in sparse matrix format) and
sparse vector/matrix computations in the code to speed up
computations. spfabia applications: (a) detecting >identity by
descent< in next generation sequencing data with rare variants, (b)
detecting >shared haplotypes< in disease studies based on next
generation sequencing data with rare variants

- fabia for non-negative factorization (parameter: non_negative)

- changed to C and removed dependencies to Rcpp

- improved update for lambda (alpha should be smaller, e.g. 0.03)

- introduced maximal number of row elements (lL)

- introduced cycle bL when upper bounds nL or lL are effective

- reduced computational complexity

BUG FIXES

- update formula for lambda: tighter approximation

- corrected inverse of the conditional covariance matrix of z

GenomicFeatures
---------------

Changes in version 1.6:

NEW FEATURES

- TranscriptDbs are now available as standard packages.  Functions that
were made available before the last release allow users to create
these packages.

- TranscriptDb objects now can be used with select

- select method for TranscriptDb objects to extract data.frames of
available annotations.  Users can specify keys, along with the
keytype, and the columns of data that they want extracted from the
annotation package.

- keys now will operate on TranscriptDB objects to expose ID types as
potential keys

- keytypes will show which kinds of IDs can be used as a key by select

- cols will display the kinds of data that can be extracted by select

- isActiveSeq has been added to allow entire chromosomes to be toggled
active/inactive by the user.  By default, everything is exposed, but
if you wish you can now easily hide everything that you don't want to
see.  Subsequence to this, all your accessors will behave as if only
the "active" things are present in the database.

SIGNIFICANT USER-VISIBLE CHANGES

- saveDb and loadDb are here and will be replacing saveFeatures and
loadFeatures.  The reason for the name change is that they dispatch
on (and should work with a wider range of object types than just
trancriptDb objects (and their associated databases).

BUG FIXES

- ORDER BY clause has been added to SQL statements to enforce more
consistent ordering of returned rows.

- bug fixes to enable DB construction to still work even after changes
in schemas etc at UCSC, and ensembl sources.

- bug fixes to makeFeatureDbFromUCSC allow it to work more reliably (it
was being a little too optimistic about what UCSC would actually
supply data for)

GenomicRanges
-------------

Changes in version 1.6.0:

NEW FEATURES

- seqlevels() and seqinfo() setters have a new arg ('force', default is
FALSE) to force dropping sequence levels currently in use.

- Seqinfo objects now have a genome column that can be accessed with
genome() getter/setter.

- "pgap" method for c(x="GRanges", y="GRanges").

- Add comparison (==, <=, duplicated, unique, etc...) and ordering
(order, sort, rank) methods for GenomicRanges objects.

- Add "flank" method for GRangesList objects.

- Add "isDisjoint" and "restrict" methods for GRanges and GRangesList
objects.

- Add GRangesList constructor makeGRangesListFromFeatureFragments().

- Add "names" and "names<-" methods for GappedAlignments objects.

- Add 'ignore.strand' arg to a number of methods: -
findOverlaps,GRangesList,RangesList -
findOverlaps,GappedAlignments,ANY - findOverlaps,ANY,GappedAlignments

- 'shift' and 'weight' arguments of "coverage" method for GenomicRanges
objects now can be numeric vectors in addition to lists.

- Add "c" method for GappedAlignments objects.

SIGNIFICANT USER-VISIBLE CHANGES

- readGappedAlignments() supports 2 new arguments: (1) 'use.names'
(default is FALSE) for using the query template names (QNAME field in
a SAM/BAM file) to set the names of the returned object, and (2)
'param' (default is NULL, otherwise a ScanBamParam object) for
controlling what fields and which records are imported.
readGappedAlignments() doesn't support the 'which' arg anymore.

- The names of a GRanges/GRangesList/GappedAlignments object are not
required to be unique anymore.

- By default, the rownames are not set anymore on the DataFrame
returned by elementMetadata() on a
GRanges/GRangesList/GappedAlignments object.

- 'width' arg of "coverage" method for GenomicRanges objects now must
be NULL or numeric vector (instead of NULL or list).

DEPRECATED AND DEFUNCT

- Deprecate countGenomicOverlaps() in favor of summarizeOverlaps().

- Deprecate grg() in favor of granges().

BUG FIXES

- Fix bug in "pintersect" methods operating on GappedAlignments
objects.

GOSemSim
--------

Changes in version 1.11.2:

- Using methods implemented in DOSE for semantic similarity
calculation.

Heatplus
--------

Changes in version 2.0.0:

- Add both row- and column annotation via annHeatmap2

- New function regHeatmap as replacement for heatmap_2

- New function annHeatmap as replacement for heatmap_plus

- Support for heatmap_2, heatmap_plus will be dropped (deprecation
warning)

- Mutiple continous covariates for annotation possible (picketPlot)

- Simple legend also for annotated heatmaps (annHeatmap2)

- More detailed interface for specifiying graphical parameters
(annHeatmap2)

- Green-to-red palette as new default (annHeatmap2)

IRanges
-------

Changes in version 1.12.0:

NEW FEATURES

- Add "relist" method that works on a List skeleton.

- Add XDoubleViews class with support of most of the functionalities
available for XIntegerViews.

- c() now works on XInteger and XDouble objects (in addition to XRaw
objects).

- Add min, max, mean, sum, which.min, which.max methods as synonyms for
the view* functions.

SIGNIFICANT USER-VISIBLE CHANGES

- Views and RleViewsList classes don't derive from IRanges and
IRangesList classes anymore.

- When used on a List or a list, togroup() now returns an integer
vector (instead of a factor) for consistency with what it does on
other objects (e.g. on a Partitioning object).

- Move compact() generic from Biostrings to IRanges.

- Drop deprecated 'multiple' argument from "findOverlaps" methods.

- Drop deprecated 'start' and 'symmetric' arguments from "resize"
method for Ranges objects.

DEPRECATED AND DEFUNCT

- Using as.data.frame() and or as( , "data.frame") on an AtomicList
object is deprecated.

- Deprecate all coercion methods from AtomicList to atomic vectors.
Those methods were unlisting the object, which can still be done with
unlist().

- Deprecate the Binning class.

- Remove defunct overlap() and countOverlap().

BUG FIXES

- togroup() on a List or a list does not look at the names anymore to
infer the grouping, only at the shape of the list-like object.

- Fix 'relist(IRanges(), IRangesList())'.

- Fix 'rep.int(Rle(), integer(0))'.

- Fix some long-standing issues with the XIntegerViews code (better
handling of "out of limits" or empty views, overflows, NAs).

isobar
------

Changes in version 1.0.0:

- first Bioconductor version (version bump to 1.0.0)!

- slot name reporterMasses is renamed to reporterTagMasses to fix clash
with method reporterMasses which fetches
assayData(ibspectra)[["mass"]]

- slot name reporterNames is renamed to reporterTagNames to distinguish
from deprecated Biobase::reporterNames

- added option 'scan.lines' to readIBSpectra: read mgf files in parts
for large MGF files

- use a function for protein reporting: create.reports.R
properties.conf will be replaced by properties.R

Changes in version 0.2.5:

- MSnbase support: Added functions to coerce from MSnSet to IBSpectra
and vice versa. Added Msnbase to Suggests.

- support for multiple classes added

- updated Perl parsers: mascotParser2.pl and pidresParser2.pl instead
of isobarXParser.pl resulting XML files can be converted to id.csv
using psx2tab2.pl

- prob otion for readIBSpectra worked errornous - fixed (thanks to
Xavier Robin)

- added property use.na: Use NA values for quantification

- various Analysis Report beautifications (thanks to Xavier Robin)

- varous bug fixes

Changes in version 0.2.4:

- improved Vignette descriptions, added CITATION (still UTF-8 error)

- added possibility to revert Mascot escaped TITLEs

- if proteins are excluded w/ subsetIBSpectra, exclude all it's
peptides, not only reporter-specific ones

- fix error introduced in 0.2.3: When multiple MGFs were read, an false
error occured that not all id spectra could be matched

- add property 'author' for LaTeX reports

- section 'Significantly regulated proteins' not shown anymore by
default added property show.significant.proteins to reenable

- added properties isotope.impurities and fragment.outlier.prob

- bug fixes:

- naming not correct when class labels contain NAs

- numeric class labels are not handled correctly

- added naRegion to noise model

- data is now stored before normalization. Those values are then used
to normalize. (Thanks to observation of Chris Bielow)

Changes in version 0.2.3:

- specify combination matrix for proteinRatios and in properties.conf

- improved logging of IBSpectra creation and normalization

- fix: maplot crashed on all NA channels

- NA names in PDF report section 'Not quantified proteins' removed

- allow for NA class labels - they are ignored in the comuptation of
ratios

Changes in version 0.2.2:

- re-added ratio vs intensity plot in QC report

- issue warning when summarize property is incorrectly defined

- create cachedir if it does not exist

- estimateRatio.group_specific_proteins renamed to
quant.w.grouppeptides

- sanitize analysisname, uniprotlink, and subsection names for LaTeX

- use fancyhdr instead of fanctheadings

- added argument require.reporter.specific to reporterProteins

Changes in version 0.2.1:

- Bug fix: as.data.frame generated ions/mass colnames with a 'X' in
front

Changes in version 0.2:

- Published online with JPR Publication

lumi
----

- Adding supports for Infinium 450K methylation microarrays, which
includes updating the color bias correction and other related
functions.

MSnbase
-------

Changes in version 1.1.28:

- added pgfSweave to Suggests <2011-10-23 Sun>

Changes in version 1.1.27:

- fixed bug in readIspyData <2011-10-07 Fri>

- removed (unexported) ratio code <2011-10-12 Wed>

- added processing description when MSnSet['ing <2011-10-12 Wed>

- man typos corrected <2011-10-14 Fri>

- changed readMzXMLData to readMSData in tests <2011-10-14 Fri>

- expecting warning for readMzXMLData in test_io <2011-10-14 Fri>

Changes in version 1.1.26:

- extractSpectra deprecated <2011-10-05 Wed>

- small changes in cache.R functions <2011-10-05 Wed>

- changed [, extractSpectra and extractPrecSpectra to update the .cache
slot <2011-10-05 Wed>

- MSnExp show method uses .cache when level==1. <2011-10-05 Wed>

- Finished level 1 cache implementation: leads to an average 11.8 times
faster MSnExp show method. <2011-10-05 Wed>

Changes in version 1.1.25:

- added .cache slot to pSet class <2011-10-03 Mon>

- added pSet initialize method to set .cache and .cache$level as '0'.
<2011-10-03 Mon>

- Check that level is defined in .cache and env is locked in pSet
validity method. <2011-10-03 Mon>

- updated itraqdata.RData <2011-10-03 Mon>

- Deprecated readMzXMLData, added defunct.Rd <2011-10-03 Mon>

- changed readMzXMLData to readMSData in man <2011-10-03 Mon>

- updated show MSnExp for speed <2011-10-03 Mon>

- updated read*Data is support cache = [0|1]. <2011-10-03 Mon>

- updated precScanNum to use sapply(spectra(...), precScanNum) instead
of unlist(eapply(assayData(...), precScanNum)) to preserve splectra
order. <2011-10-03 Mon>

- new cache.R file with @.cache related code <2011-10-03 Mon>

Changes in version 1.1.24:

- in readMgfData, SCANS is not used to populate acquisitionNum anymore,
as several scans might have been combined upstreams <2011-09-26 Mon>

- added fillUp function in utils.R and exported <2011-09-27 Tue>

- new MSnbase-io vignette <2011-09-28 Wed>

Changes in version 1.1.23:

- fixed writeMgfData("MSnExp") <2011-09-21 Wed>

- readMgfData now works when mz and intensty are separated by a '\t'
(as exported by PRIDE Inspector) <2011-09-21 Wed>

- show("MSnExp") now works without retention time <2011-09-21 Wed>

- updated mgf2Spectrum2 to make it faster <2011-09-21 Wed>

- fixed missing fromFile slot in data created from readMgfData that
prevented calling header <2011-09-21 Wed>

- readMgfData now creates a fData based on the peak header <2011-09-21
Wed>

- modified getCurveWidth to work with centroided data <2011-09-21 Wed>

- fixed bug getCurveWidth <2011-09-21 Wed>

Changes in version 1.1.22:

- removed (internal) Mascot query link column in readIspyData to work
with latest ouput version <2011-09-12 Mon>

- removed the fillUp part in readIspyData <2011-09-14 Wed>

- exported readIspyData <2011-09-20 Tue>

- removed link to proteomics sig list <2011-09-20 Tue>

- removed url in DESCRIPTION <2011-09-20 Tue>

- Spectrum2 slot ms1scan renamed to precScanNum and populated using
mzR::header()$precursorScanNum. Updated affected methods/functions
and manual pages. New accessor method precScanNum is exported.  THIS
CHANGE IS NOT COMPATIBLE WITH PREVIOUSLY CREATED MSnSet OBJECTS!
<2011-09-20 Tue>

Changes in version 1.1.21:

- updated write.exprs to add fData columns <2011-09-08 Thu>

- added write.exprs and readMSnSet unit test <2011-09-08 Thu>

Changes in version 1.1.20:

- added a readMSnSet function <2011-09-08 Thu>

- added a write.MSnSet("MSnSet") method <2011-09-08 Thu>

- vignette/man updates - mainly data import section documenting
readMSData and readMgfData <2011-09-08 Thu>

- incorporating mzR io frame work <2011-09-05 Mon>

- use mzR's peaksCount and header generics <2011-09-05 Mon>

- added test to check that readMzXMLData and readMSData give same
output <2011-09-05 Mon>

- added readMSData.Rd doc file <2011-09-05 Mon>

- added Author@R field in DESCRIPTION <2011-09-07 Wed>

- compressed/resaved (using resaveRdaFiles) itraqdata.RData file to
<2011-09-07 Wed>

Changes in version 1.1.18:

- read support for mgf files contributed by Guangchuang Yu <2011-07-06
Wed>

- added as.ExpressionSet.MSnSet and setAs methods and updated MSnSet
doc <2011-07-09 Sat>

- fixed read/write mgf compatibility <2011-09-01 Thu>

- added mgf io test <2011-09-01 Thu>

- exported and document mfg read/write support <2011-09-01 Thu>

- added centroided parameter to rawToSpectrum[1|2] to set this directly
at object creation in readMzXmlData <2011-09-01 Thu>

- other minor changes in Rd files <2011-09-01 Thu>

- added warning checks for combineFeatures when verbose=TRUE in
test_MSnSet.R <2011-09-02 Fri>

Changes in version 1.1.17:

- updated itraqdata manual <2011-06-27 Mon>

- fixed bug in MSnSet validity method - msg was initialised to NULL
when testing Biobase:::isValidVersion and
Biobase::assayDataValidMembers <2011-07-02 Sat>

- added t.MSnSet method <2011-07-02 Sat>

- document t.MSnSet method <2011-07-05 Tue>

- test for t-MSnSet <2011-07-05 Tue>

- new "["-MSnSet method to properly subset qual slot <2011-07-06 Wed>

- updated MSnSet validity method to check qual dims <2011-07-06 Wed>

- combineFeatures now drops the spectrum-specific qual slot <2011-07-06
Wed>

- test for "["-MSnSet <2011-07-06 Wed>

Changes in version 1.1.16:

- added addVigs2WinMenu call in .onAttach <2011-06-26 Sun>

- added plotMzDelta paragraph in vignette <2011-06-26 Sun>

Changes in version 1.1.15:

- new names and description methods for ReporterIons <2011-06-16 Thu>

- added validObject() checks <2011-06-16 Thu>

- new removeReporters method <2011-06-16 Thu>

- plotMzDelta is exported and documented <2011-06-17 Fri>

Changes in version 1.1.14:

- Adding plotMzDelta QC plot (not yet exported), contributed by
Guangchuang Yu <2011-06-14 Tue>

- created a locked environment to store amino.acids dataframe
<2011-06-14 Tue>

- TODO plotMzDelta documentation and vignette section - man DONE in v
1.1.15 , vignette DONE in v 1.1.16

Changes in version 1.1.13:

- changed pSet [[ method <2011-06-13 Mon>

- additional updates pSet [[ method <2011-06-13 Mon>

- added MSnExp subsetting error tests <2011-06-13 Mon>

- created new plotting-dataframe.R file with former plotting utils
functions, now renamed plot*.header <2011-06-13 Mon>

Changes in version 1.1.12:

- added invisible(NULL) to all show methods <2011-06-08 Wed>

- added centroided argument to plot.MSnExp <2011-06-09 Thu>

Changes in version 1.1.11:

- harmonised MSnExp and Spectrum plot axes labels <2011-05-18 Wed>

- Added plotting customisation section in vignette <2011-05-18 Wed>

- updated signature of plot2d method to "MSnExp" only <2011-05-18 Wed>

- added/exported plotDensity methods <2011-05-18 Wed>

- started QC vignette section <2011-05-18 Wed>

- added preprocSelection and preprocSelectionTable functions
<2011-05-18 Wed>

- TODO document preprocSelection[Table] in vignette

- reduces plot2d-figure and plotDensity-figure sizes using png
<2011-05-19 Thu>

- added plotDensity doc <2011-05-19 Thu>

- added round param to preprocSelection[Table] <2011-05-19 Thu>

- preprocSelection[Table] documented and exported <2011-05-19 Thu>

- fixed bug in plot.Spectrum1 <2011-05-24 Tue>

- changed removePeaks setGeneric explicit signature <2011-05-26 Thu>

- added MassSpectrometry biocView <2011-05-27 Fri>

Changes in version 1.1.10:

- minor updates in demo vignette <2011-05-13 Fri>

- added plot argument to plot methods <2011-05-16 Mon>

- fix in makeImpuritiesMatrix <2011-05-16 Mon>

- added meanSdPlot MSnSet method <2011-05-17 Tue>

- minor cosmetic fix in purityCorrect error message <2011-05-17 Tue>

- added method="sum" to Spectrum/MSnExp normalisation <2011-05-17 Tue>

- typo in MSnSet-class.Rd corrected <2011-05-17 Tue>

Changes in version 1.1.9:

- exporting normali[s|z]e methods for MSnSet, Spectrum and MSnExp
<2011-05-12 Thu>

- added quantile normalisation <2011-05-12 Thu>

- added quantile.robust normalisation <2011-05-12 Thu>

- added vsn2 normalisation <2011-05-12 Thu>

- added normalise manual <2011-05-12 Thu>

- included normalisation section in vignette <2011-05-13 Fri>

- more vignette updates <2011-05-13 Fri>

- updated plotting methods to round MZ in title <2011-05-13 Fri>

Changes in version 1.1.8:

- added writeMgfData method from spectra and experiment <2011-05-09
Mon>

- added writeMgfData manual <2011-05-09 Mon>

- added itraqdata data set and updated vignette to use it <2011-05-11
Wed>

- updated man pages to use tiny dummyiTRAQ.mzXML <2011-05-11 Wed>

- added makeImpuritiesMatrix function <2011-05-11 Wed>

Changes in version 1.1.7:

- reporter ions purity correction method, man and test <2011-05-09 Mon>

- updated vignette <2011-05-09 Mon>

Changes in version 1.1.6:

- added incomplete dissociation and spectral counting sections to demo
vignette <2011-05-06 Fri>

- added bioc-sig-proteomics link to foreword <2011-05-06 Fri>

- type in foreword <2011-05-06 Fri>

Changes in version 1.1.5:

- added combineFeatures example in demo vignette <2011-05-05 Thu>

Changes in version 1.1.4:

- readIspyData updated to return updated factors <2011-05-03 Tue>

- added unexported/undocumented combineFeatures function for MSnSets
<2011-05-03 Tue>

- added basic tests for combineFeatures <2011-05-03 Tue>

- added combineFeatures manual <2011-05-04 Wed>

- exportig combineFeatures <2011-05-04 Wed>

Changes in version 1.1.3:

- as.data.frame.Spectrum columns now are (first) mz and (second)
intensity <2011-04-27 Wed>

- exporting as.data.frame.Spectrum and coerce <2011-04-27 Wed>

Changes in version 1.1.2:

- Simplified quantify generic signature - now only object argument
<2011-04-19 Tue>

- Added strict parameter to quantify method, man updated, added
relevant test

- Added illustrative plot for quantitation methods in MSnbase-demo
vignette <2011-04-19 Tue>

- Added illustrative plot for data pre-processing (removePeaks and
clean) in MSnbase-demo vignette <2011-04-20 Wed>

- No warnings are issued anymore when peaks expands outside of
mz(reporters) +/- width(reporters).  See ?quantify on how to check
this manually. <2011-04-19 Tue>

- No warnings are issued anymore when reporter peaks are missing.  See
?quantify on how to check this manually. <2011-04-20 Wed>

- pSet validity warns if length(unique(msLevel(object))) > 1, rather
than != 1.  The latter triggered a warning for a new("MSnExp").
<2011-04-20 Wed>

Changes in version 1.1.1:

- added setAs data.frame and as.data.frame methods for Spectrum objects
<2011-03-29 Tue>

- support for uncentroided MS2 spectra plots <2011-03-31 Thu>
<2011-04-02 Sat>

- support for uncentroided MS1 spectra plots <2011-04-02 Sat>

- minor modification to readIspyData <2011-04-04 Mon>

- removed centroided slot from MSnProcess and added to individial
Spectrum instances. Relevant for Velos HCD (profile)/CID
(uncentroided) data <2011-04-04 Mon>

- modified readMzXmlData accordingly <2011-04-04 Mon>

- added validObject(new(...)) tests for each class <2011-04-04 Mon>

- added centroided[<-] methods to Spectrum and pSet <2011-04-04 Mon>

- added 'keepAll' parameter to readIspyData <2011-04-11 Mon>

qpgraph
-------

Changes in version 1.10:

NEW FEATURES

- estimation of mixed interactions for genetical genomics data via
mixed graphical model theory.

- qpBoundary() to explore sparsity in graphs estimated from
non-rejection rates.

- MLE of covariance matrices via the Hastie-Tibshirani-Friedman (HTF)
algorithm (Hastie, Tibshirani and Friedman, 2009, pg. 634), which
enables a much faster simulation of these matrices via qpG2Sigma()
than with the previous version based on the IPF algorithm, and much
faster MLE of partial correlations via qpPAC().

- uniform sampling of d-regular graphs in qpRndGraph() with the
algorithm of Steger and Wormald (1990)

BUG FIXES

- proper handling of master node identification with snow so that
parallel computations work again with snow versions > 0.3-3

rqubic
------

Changes in version 1.8.1 (2011-08-02):

- Fix errors in documentations

Changes in version 1.8.0 (2011-07-11):

- Submission to Bioconductor

Changes in version 0.99.1 (2011-08-11):

- Rename quantile.discritize to quantileDiscritize in order to avoid
conflicts with the S3 method nomenclature

Changes in version 0.99.0 (2011-08-05):

- Re-starting version indexing due to submission to the Bioc repository

Rsamtools
---------

Changes in version 1.6.0:

NEW FEATURES

- TabixFile, indexTabix, scanTabix, yieldTabix index (sorted,
compressed) and parse tabix-indexed files

- readBamGappedReads(), bamFlagAsBitMatrix(), bamFlagAND()

- Add use.names and param args to readBamGappedAlignments(); dropped
which and ... args.

- PileupFiles, PileupParam, applyPileup for visiting several BAM files
and calculating pile-ups on each.

- Provide a zlib for Windows, as R does not currently do this

- BamFileList, BcfFileList, TabixFileList, FaFileList clases extend
IRanges::SimpleList, for managings lists of file references

- razfFa creates random access compressed fasta files.

- count and scanBam support input of larger numbers of records;
countBam nucleotide count is now numeric() and subject to rounding
error when large.

- Update to samtools 0.1.17

- asBcf and indexBcf coerces VCF files to BCF, and indexes BCF

- Update to samtools 0.1.18

- scanVcf parses VCF files; use scanVcf,connection,missing-method to
stream, scanVcf,TabixFile,*-method to select subsets. Use unpackVcf
to expand INFO and GENO fields.

SIGNIFICANT USER-VISIBLE CHANGES

- ScanBamParam argument 'what' defaults to character(0) (nothing)
rather than scanBamWhat() (everything)

- bamFlag returns a user-friendly description of flags by default

BUG FIXES

- scanBam (and readBamGappedAlignments) called with an invalid or
character(0) index no longer segfaults.

- scanBcfHeader parses values with embedded commas or =

- scanFa fails, rather than returns incorrect sequences, when file is
compressed and file positions are not accessed sequentially

- scanBcf parses VCF files with no genotype information.

- scanBam called with the first range having no reads returned invalid
results for subsequent ranges; introduced in svn r57138

- scanBamFlag isPrimaryRead changed to isNotPrimaryRead, correctly
reflecting the meaning of the flag.

Rsubread
--------

Changes in version 1.4.0:

NEW FEATURES

- Detect exon junctions (subjunc function).

- Make detection calls for each genes to tell whether or not they are
expressed (detectionCall function).

- Count number of reads falling into each gene or each exon
(featureCounts function).

RTCA
----

Changes in version 1.5.2:

- plateView function had a typo which is now fixed

Changes in version 1.5.1:

- Option "skipWells" of function "parseRTCA" is renamed as "maskWells".
It now accepts either well name(s) or a regular expression pattern as
character string.

- Internally the parseRTCA will not skip the first line automatically
as it did before. Instead it first detects whether there are lines
with duplicated time points; if so, it keeps the first appearance and
removes the rest items.

ShortRead
---------

Changes in version 1.11:

NEW FEATURES

- trimTails to trim low quality trailing nucleotides

- trimEnds to remove arbitrary (vectors of) letters from reads or
qualities

- FastqStreamer to iterate over a fastq file

- FastqFile, FastqFileList to represent fastq files

SIGNIFICANT USER-VISIBLE CHANGES

- writeFastq has argument full, default value FALSE, disabling printing
of identifier a second time in '+' line

- srapply requires that options(srapply_fapply="parallel") or
options(srapply_fapply="Rmpi") to enable parallel processing via
fapply

BUG FIXES

- SolexaRealign, SolexaAlign, and SolexaResult transposed strand
inforamtion

- FastqSampler segfaulted on some files

- writeFasta had a semi-documented argument mode; it is now documented
and as a consequence dis-allows argument 'append' that would
previously have been passed to underlying methods.

TEQC
----

Changes in version 2.0.0:

- new function 'TEQCreport' creates an html report with standard TEQC
analysis results

- besides bed files, now also BAM files can be used as input for
'get.reads'

- 'get.reads' and 'get.targets' now only read the columns required for
the analysis from the respective bed files

- when read IDs include '#0/1' and '#0/2' (to indicate read 1 and read
2 of a pair), those characters will be removed from the IDs within
'get.reads'. The reason is that in 'reads2pairs' the two IDs of a
read pair have to be identical.

- the package now depends on packages Rsamtools and hwriter

- 'chrom.barplot',' fraction.reads.target', 'insert.size.hist', and
'duplicates.barplot' now also can deal with with 'reads2pairs' output
having the two elements 'singleReads' and 'readpairs'

Changes in version 1.1.2:

- bug fix in 'coverage.target' and 'coverage.GC' (in very large
datasets, global coverage average and standard deviation were not
calculated)

Changes in version 1.1.0:

- fix in 'reads2pairs': when the two reads of a read pair map to
different chromosomes, they will be returned within the 'singleReads'
element of the output (before function gave en error in case of such
read pairs)

- added optional argument 'max.distance' to 'reads2pairs'; when the
reads of a read pair are further apart than 'max.dist' bases, they
will be added to the 'singleReads' element of the output instead to
the 'readpairs' element

VariantAnnotation
-----------------

Changes in version 1.0.0:

NEW FEATURES

- readVcf() for reading and parsing VCF files into a
SummarizedExperiment

- locateVariants() and predictCoding() for identifying amino acid
coding changes in nonsynonymous variants

- dbSNPFilter() and regionFilter() for filtering variants on membership
in dbSNP or on a particular location in the genome

- access to PolyPhen and SIFT predictions through keys() , cols() and
select() methods. See ?SIFT or ?PolyPhen.

BUG FIXES

- No changes classified as 'bug fixes' (package under active
development)

xps
---

Changes in version 1.12.0:

VERSION xps-1.11.12

- update functions: plotBoxplot(), plotAffyRNAdeg()

VERSION xps-1.11.11

- update functions: plotImage(), plotAffyRNAdeg()

VERSION xps-1.11.10

- add functions: plotImage(), plotBoxplot()

- update methods image()

VERSION xps-1.11.9

- add method: pcaplot()

- update method plotAffyRNAdeg()

VERSION xps-1.11.8

- add methods: corplot(), madplot()

VERSION xps-1.11.7

- update method plotAffyRNAdeg()

VERSION xps-1.11.6

- add methods for RNA degradation plots: xpsRNAdeg(), plotAffyRNAdeg()

- add man pages

VERSION xps-1.11.5

- correct minor bug in XPSProcessing.cxx: ExportExprTreeInfo()

- update method image()

- add man pages

VERSION xps-1.11.4

- add new class QualTreeSet to add quality control features

- add functions qualify() and fitQC() and derived functions

- add method image() for residual plots of IVT, Gene ST, Exon ST and
plate arrays

- add new plots coiplot() and borderplot()

- update boxplot(), callplot(), image() to be independent of slot
'data'

VERSION xps-1.11.3

- add function trma()

VERSION xps-1.11.1

- update DESCRIPTION to correct SystemRequirements to root_v5.27.04

- update function READ_WSTRING() to handle big endian for PPC
Packages removed from the release
=================================

The following packages are no longer in the release: GeneSpring,
GeneTraffic, makePlatformDesign, Rdbi, RdbiPgSQL, rflowcyt, Rredland,
RTools4TB, Ruuid, simulatorAPMS.

[bioc]: http://bioconductor.org
[install]: http://bioconductor.org/install/
[ami]: http://bioconductor.org/help/bioconductor-cloud-ami/
