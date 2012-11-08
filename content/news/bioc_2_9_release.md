November 1, 2011

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


a4
--

Changes in version 1.1.2:

- compact pdf vignette

a4Base
------

Changes in version 1.1.4:

- fix in limmaTwoLevels

Changes in version 1.1.3:

- update for revised glmnet package

Changes in version 1.1.1:

- on BioConductor

a4Classif
---------

Changes in version 1.1.3:

- adapt to new lognet objects from glmnet

Changes in version 1.1.1:

- on BioConductor

a4Core
------

Changes in version 1.1.3:

- remove duplicated code for glmnet type of objects; add topTable
  method for 'elnet' objects and accompanying print method
  print.topTableElnet

Changes in version 1.1.2:

- add topTable method for 'lognet' objects and accompanying print
  method print.topTableLognet

Changes in version 1.1.1:

- on BioConductor

a4Reporting
-----------

Changes in version 1.1.1:

- fixes for latest glmnet versions

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

Changes in version 1.3.3 (2011-10-19):

- Added more "Sys.sleep" to see if we can get it not to crash in
  Windoze.

Changes in version 1.3.2 (2011-10-18):

- Added unloading of rlecuyer, to prevent problems with cleanEx during
  R CMD check.

Changes in version 1.2.1 (2011-06-06):

- Changed "if" in stop.na.inf, as per Dundan Murdoch's suggestion.

affxparser
----------

Changes in version 1.25.1 (2011-09-27):

- Maintainer email was updated.

Changes in version 1.25.0 (2011-04-13):

- The version number was bumped for the Bioconductor devel version.

AnnotationFuncs
---------------

Changes in version 1.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Added function 'getOrthologs' and auxcillary functions that speeds up
  ortholog mappings in the Inparanoid databases.

- Added parameter 'simple' to 'translate' for improved usage with
  lapply.

aroma.light
-----------

Changes in version 1.21.2 (2011-10-10):

- Updated robustSmoothSpline() such that it works with the new
  "uniqueness" scheme of smooth.spline() in R v2.14.0 and newer. It is
  tricky, because robustSmoothSpline() is a reiterative algorithm which
  requires that the choosen "unique" x:s does not change in each
  iteration.  Previously, 'signif(x, 6)' was used to identify unique
  x:s, which gives the same set of values when called twice, whereas
  this is not true for the new choice with 'round((x - mean(x))/tol)'.

Changes in version 1.21.1 (2011-06-26):

- Added argument 'aspectRatio' to plotMvsA().  It can be used to adjust
  the range of the 'Mlim' argument relative to the 'Alim' argument.

Changes in version 1.21.0 (2011-04-13):

- The version number was bumped for the Bioconductor devel version.

ArrayExpressHTS
---------------

Changes in version 1.3.1:

NEW FEATURES

- Added test dataset, which is a cut version of E-GEOD-16190 experiment

- Added normalization of the reference folder

- Added automatic fixing of reference .fai index, namely duplicated
  chromosome entries

- Added automatic detection of quality scale where possible

- Added and option to specify quality scale for custom experiments

- Performed performance and memory optimization

- Removed options$reference_version. Since the pipeline allows
  processing of multi-species datasets, a single version cannot fit
  into multiple organism paradigm. The pipeline always runs on the
  latest version of reference.  Support for custom versions can be
  added upon request.

- Added memory monitoring

BUG FIXES

- Fixes to accommodate changes in Ensembl interfaces, required for
  preparation of reference and annotation data.

- Fixes to accommodate changes in ENA interfaces, required for querying
  metadata & downloading datasets for ArrayExpress experiments.

- Fixed examples in prepareReference help page

- Fixed ArrayExpressHTSFastQ scenario without .sdrf

- Fixed production of final alignment report

Changes in version 1.3.0:

- No changes. Version changed due to package propagation.

beadarray
---------

Changes in version 2.4:

MAJOR CHANGES

- The example data required to run examples and the vignette has been
  removed from the package. The experimental data package
  beadarrayExampleData has been created for this purpose

- Mapping of bead-level ArrayAddress IDs to Illumina IDs is now
  performed by loading a companion annotation package. e.g.
  illuminaHumanv3.db

- The addFeatureData function has been added to simplify the annotation
  of summary data objects

- The ggplot2 library is now used to produce boxplots of summary data
  with the option of including columns of featureData or phenoData as
  factors in the plot

- The combinedControlPlot has been added, which uses ggplot2 to
  consolidate the bead-level control plots into a single graphic

- BASH and HULK can now only operate on a single array at a time and
  accept an outlier-calling method as a parameter

OTHER CHANGES

- The imageplot function now uses ggplot2

- Accessor functions added for metrics, P95 value and signal-to-noise
  ratio

- The default background correction when running readIllumina with
  useImages = TRUE is now the medianBackground method

- getAnnotation replaced with annotation accessor function

- setAnnotation replaced with annotation <- replacement method

- checkPlatform replaced with suggestAnnotation function

- The platformSigs object was added and contains lists of
  ArrayAddressIDs for all know expression platforms

- The outlier calling method, squeezedVarOutlierMethod has been added,
  which shrinks the observed variance for a bead-type towards the
  predicted variance based on all bead-types on the array-section

- Errors in BASH that led to the einvasions and dinvasions arguments
  not being respected have been corrected

- lower-level BASH functions (e.g. BASHDiffuse) will no longer
  calculate a Neighbours matrix if one is not supplied (ergo one must
  be supplied)

- beadarray used to have problems if no outlier removal method were set
  in summarization (important if the outliers have already been removed
  by Illumina's scanner software). The function noOutlierMethod is
  provided to address this

- findAllOutliersE and findAllOutliersIgnore are to be deprecated, used
  internally by BASH have been removed

Changes in version 2.3:

NEW FEATURES

- checkRegistration() has been completely re-written to provide a far
  more meaningful result in a new class, beadRegistrationData

- beadRegistrationData objects can be passed to boxplot.

- The summarize function now has a default channel parameter and
  useSampleFac = FALSE

BUG FIXES in 2.3.6

- Fixed problem when providing arrays argument to checkRegistration()

- Fixed problem in readIllumina when useImages = TRUE

- Patched error with analyseDirectory() (thanks to Juerg Straubhaar for
  this)

- Fixed problem when providing arrays argument to checkRegistration()

- Fixed problem in readIllumina when useImages = TRUE

- Patched error with analyseDirectory() (thanks to Juerg Straubhaar for
  this)

BUG FIXES in 2.3.5

- Added support to readIllumina() for files that use something other
  than a period as the decimal point character.

BUG FIXES in 2.3.4

- Fixed issues when using a combination of iScan data and a locs file
  to generate neighbours, manifesting as problems with imageplot and
  BASH.

BUG FIXES in 2.3.3

- Stopped the presence of any other .txt files in a directory resulting
  in undesirable behaviour of readIllumina()

BUG FIXES in 2.3.2

- readBeadLevelTextFile will now take (and use) a seperator argument.

- locs information correctly taken from a bab file if required.  This
  effects functions such as generateNeighbours(), BASH(), HULK(),
  imageplot()

BUG FIXES in 2.3.1

- Image processing functions can now be passed matrices containing
  either integer or numeric values.

- Passed forceIScan argument to functions called by readIllumina()

BeadDataPackR
-------------

Changes in version 1.5:

BUG FIXES

- Beads with negative coordinates can now be included in the bab file

- If overlapping segments are found a warning is printed and the full
  locs index is used.

- Coordinate values returned by readCompressedData() are now rounded to
  match the precision of the original text file.

CAMERA
------

Changes in version 2011-28-04:

- Fix bug in getpspectra, which incorrect label isotope peaks

Changes in version 2011-24-05:

- Fix bug in findIsotopes, which could cause a crash, if maxiso was
  higher than 3

Changes in version 2011-23-09:

- Version 1.9.7

- Bugfix in calcIsotopes, causes groupCorr with calcIso to crash with
  (Error in rbind(resMat ...)

- Bugfix in groupCorr with given xcmsRaw

Changes in version 2011-22-08:

- Version 1.9.6

- Bugfix in plotEICs (Error in pks[, 1] : incorrect number of
  dimensions)

Changes in version 2011-20-25:

- Version 1.9.9

- added some "drop=FALSE" to fix "Error in isomatrix[, 1] : incorrect
  number of dimensions" error

Changes in version 2011-20-10:

- Version 1.9.8

- Correct rule table extended_adducts_pos.csv (typo in proton mass)

Changes in version 2011-12-05:

- plotPsSpectrum now accepts additional parameters for plot

Changes in version 2011-04-04:

- Fix bug in findAdducts, if pseudospectrum mass list has only NA
  values

Changes in version 2011-02-08:

- Version 1.9.5

- Bugfix for findIsotopes if ppm was very high it could occur that one
  peak is assigned as two or more isotope peaks

- Version 1.9.4

- Add parameter polarity to xsAnnotate constructor

- getPeaklist and getpspectra returns now correct annotation of
  negative charged ions [M]-

- Add snow as additonal possibilty for parallel processing

- Add function cleanParallel to clean up with spawned slave processes

ChIPpeakAnno
------------

Changes in version 2.2.0:

- Find the peaks with bi-directional promoters with summary statistics
  (peaksNearBDP).

- Summarize the occurrence of motifs in peaks
  (summarizePatternInPeaks).

- Add other IDs to annotated peaks or enrichedGO (addGeneIDs)The
  function makeVennDiagram supports 4-way venn diagram.

- Enriched Gene Ontology (GO) terms are annotated with a list of genes
  which can be traced back to peaks.

- FAQ: http://pgfe.umassmed.edu/ChIPpeakAnno/FAQ.html.

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

cqn
---

Changes in version 0.99:

- Version bump as well as added to Bioconductor devel branch.

- Wrote vignette.

- Added additional manpages as well as examples.

- Added data from Montgomery et al.

- Added NAMESPACE file.

- Added manpage for cqn.

- Renamed cqn, cqn2.

- Added NEWS.Rd file.

- Added vignette skeleton.

- Added CITATION file.

cummeRbund
----------

Changes in version 0.99.5:

- Significant speed improvements to readCufflinks() for large cuffdiff
  datasets.

- Tables written first then indexed.

- Added slot accessor methods to avoid using slots directly.

Changes in version 0.99.4:

- Second beta release and submission to Bioconductor

Changes in version 0.1.3 (2011-08-18):

- First Beta release of cummeRbund and submission to Bioconductor for
  review and hosting.

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

- refactorized the way how dispersion estimation is handled,
  overhauled the vignette

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

EDASeq
------

Changes in version 0.99.1:

- The read-level EDA is now based on the classes BamFileList and
  FastqFileList of Rsamtools. The class ShortReadSet doesn't exist
  anymore.

- Changed vignette

Changes in version 0.99.0:

- Version number changed for submission to Bioconductor

Changes in version 0.0.6:

- Fixed a bug in withinLaneNormalization to work with unnamed matrices

- Fixed a bug in biasPlot to work with a single lane experiment

- Added methods MDPlot and biasBoxplot

- Changed the vignette and added an help page for newSeqExpression
  function

Changes in version 0.0.5:

- Modified coerce from SeqExpressionSet to CountDataSet to treat the
  two class comparison in a separate way

Changes in version 0.0.4:

- Added Vignette

- Adjusted xlab in the boxplot method for signature FastqQuality

- Changed getOffset and setOffset<- methods in offst and offst<- for
  compatibility with edgeR

Changes in version 0.0.3:

- Added CITATION file

- Added Help pages

- Added data and extdata for the examples

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

flowFP
------

- We now require that length of sampleClasses be the same as the number of
instances in the flowFP object.  The documentation indicated that, but 
the code didn't enforce it.

- Fixed a bug regarding setting sampleClasses for a flowFPPlex.  It erroneously assigned classes to only the first member of the plex.

- Setting sampleClasses for a fp or a pled to NULL or to a zero-length factor now results in removing sampleClasses from the fp or plex.

- Cleaned up documentation for flowFP-class and flowFPPlex-class.

flowWorkspace
-------------

- 0.6.0 This version supports importing of flowJo XML workspaces
  generated using, up to version 9.2 of the Mac OS X version of
  flowJo. We do not yet support workpsaces generated using the windows
  version of flowJo.

FEATURES:

-import flowJo XML workspace from the Mac version. (Greg Finak)

-convert GatingHierarchy objects to workflows. (Mike Jiang)

-export workflows to flowJo XML workspace readable by Mac version of
 flowJo (Mose Andre)

-netCDF support for large data sets. Based on code contributed by
 Nisat Gopalakrishnan.

KNOWN ISSUES:
	
- export to flowJo - flowWorkspace will export the correct
  compensation matrices, but they will need to apply manually via the
  menu items in flowJo.

- Navigating between graphs of populations is buggy.

- export expects to receive a workflow that contains a compensation
  view and a transformation view as the first two actions in the
  workflow.

- import from flowJ - we may not support all gating
  configurations. The parsing code was written based on examples. If
  your workspace is not being imported correctly, contact the
  authors. If your workspace is not being imported correctly because
  it's generated by flowJo for windows, it is not supported. Are you
  sure you're trying to import the XML workspace and not the .jo or
  .wsp file?

- flowJo 8.2 for MAC. - XML has some differences and there may be
  issues importing workspaces into R. We have resolved those bugs that
  we've found but we have not had the opportunity for extensive
  testing with this version.

- speed - import can slow for very large workspaces with many gates.

- This version supports importing of flowJo XML workspaces generated
  using, up to version 9.2 of the Mac OS X version of flowJo. We do
  not yet support workpsaces generated using the windows version of
  flowJo.

geneplotter
-----------

Changes in version 1.31.4:

USER VISIBLE CHANGES

- Using NEWS.Rd

- The interpretation of the argument bw by the function multidensity
  has been changed with the aim to make the behaviour more robust when
  the data ranges or sample sizes of the groups are quite different.
  Also, to allow more control, bw can now also be a vector, with
  elements corresponding to the different groups.

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

genoset
-------

1.0.0 

- Object creation at user level changed so that GenoSet, BAFSet, and
CNSet functions all call an internal input checker and object
generator function initGenoSet rather than having GenoSet() create all
3 objects using a "type" arg.

- Segmentation functions now return DataFrame objects containing one
Rle vector per sample.  These give the same effect of having big
matrices with repeated data matching the dimensions of the
un-segmented data, but are considerably smaller.  For a dataset of ~1k
SNP6 arrays, the data in the lrr matrix is about 0.8% of the size
using this strategy. In all regards, these can be treated as just
another matrix in assayData, but in special cases it is useful to take
advantage of the data in its runValue and runLength form.

- New methods defined for DataFrame to allow usage as assayData
element: colMeans annotatedDataFrameFrom

- Methods that get genome information from a GenoSet now work on
RangedData too.  RangedData and GenoSet now share some API to make
this easier. RangedData gets chr and pos, GenoSet gets names, ranges,
and elementLengths. Now chrInfo, chrIndices, genoPos, etc. can work on
either object type.

- pos and genoPos now defined as floor of mean of start and end
  positions.

1.0.6

- Substantial speed improvements for boundingIndices() and rangeSampleMeans.

1.1.7

- Added loadGC function to get local GC content and version bump for
  bioC2.9

1.1.8

- genomeOrder becomes toGenomeOrder which now takes and returns a
  GenoSet or RangedData rather than just returning the index to
  reorder such objects.  locData<- now reorders the whole GenoSet and
  assures that all of the featureNames match.

1.1.9

- Added support for big.matrix objects from bigmemory as
  assayDataElements

1.1.10

- ***API Change*** list operator "[[" no longer used to subset by
  chromosome. It reverts back to extracting a column from pData like
  other eSets. chrIndices gains a "chr" argument that serves as a fast
  way to get the indices needed to subset rows by chromosome.

1.1.11

- ***API Addition*** "k" argument to "[" can be used to subset from a
  specific assayDataElement.  Numeric and character "k"s are allowed.
  assayDataElement(ds,k)[i,j] is the same as ds[i,j,k], but "i" can be
  a RangedData or RangesList.

1.1.12

- assayDataElements can be integer matrices with levels that will
  serve as factors for now.  Please see the help for
  convertToBigMatrix.

GOSemSim
--------

Changes in version 1.11.2:

- Using methods implemented in DOSE for semantic similarity
  calculation.

GRENITS
-------

Changes in version 0.99.1 (2011-08-07):

- Initial release!

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

HTqPCR
------

Changes in version 1.7:

NEW FEATURES

- Added "geometric.mean" as a normalisation method.

BUG FIXES

- Check for all-NA features in normalizeCtData rank-invariant methods.

- Modified plotCtCategory to handle cases with only 1 sample.

- Modified plotCtCategory to make colours more consistent.

- Fix in deltaCt for normalization with NAs

- Added warning in readCtData if there are NAs in input

htSeqTools
----------

Changes in version 0.99:

- Adapted enrichedChrRegions to adapt to change in IRanges behavior

- Fixed bug in PeakLocation and stdPeakLocation when strand information
  was stored in variables of type factor

- Included adjustments for number of reads in each sample for
  giniCoverage and ssdCoverage

- Improved filterDuplReads, using mixtures of NegBinom

- Adjusted behavior of giniCoverage

- Removed readAlignedBatch, export2Aligned, export2SGR, export2SAM

Changes in version 0.2.7:

- Fixed bug in filterDuplReads method for 'RangedDataList' which caused
  it not to pass on arguments to 'RangedData' method

- Changed default value of ppOverAmp=.95 to .99 in filterDuplReads

- Added function rowLogRegLRT

- Fixed bug in alignPeaks caused by change in IRanges behavior.

- Added function coverageDiff to compute difference in coverage between
  two objects

- Fixed bug regionsCoverage caused by convertion of class "by" objects

- Added plotChrRegions function

- Added enrichedChrRegions function

- Added countHitsWindow function

- stdPeakLocation now allows to set xlim

- plot method for cmdsFit objects now automatically puts both axis on
  the same scale

- closestGene now allows argument genome to be of class data.frame.
  This avoids downloading the refflat file at each call.

- Fixed bug in enrichedRegions which caused it to crash when sample1
  had more chromosomes than sample2

- Fixed bug in filterDuplReads which caused it to crash when no reads
  needed to be filtered

- Fixed bug in regionsCoverage which caused it to crash for empty
  chromosomes

- Removed functions genePlot, rangesPlot (moved to package casper)

- Added regionsCoverage, gridCoverage, stdGrid functions

- Added chipSeqQC function

- enrichedRegions was trying to load multicores library instead of
  multicore.

- Added ssdCoverage function

- Added enrichedRegions method for RangedDataList objects. This allows
  comparing proportion of reads in islands across >2 samples.

- Added islandCounts function

- Improved memory usage and speed in filterDuplReads, which no longer
  relies on countOverlaps from package IRanges

- Changed call to parallel to multicore::parallel to avoid confusion
  with other packages

Changes in version 0.2.5:

- Added arguments maxRepeats and ppOverAmp to filterDuplReads to allow
  it to eliminate reads appearing more than k times (before only k=1
  was allowed), or with a large posterior probability of being an
  artifact

- Added ppEnrichedCounts function

- Modified enrichedRegions to make it compatible with new IRanges
  version

- Added cmds function to perform MultiDimensional Scaling on sequencing
  data

- Exact P-value calculation in enrichedRegions is now based on Fisher's
  exact test, instead of simulations, to increase speed. Changed
  parameters simulate.p.value and B for exact.

- Chi-square P-value calculation in enrichedRegions is now
  substantially faster (reduced to about 25% of previous
  implementation)

- enrichedRegions now selects regions with pval <= pvalFilter, instead
  of pval < pvalFilter. This allows to obtain all regions by setting
  pvalFilter=1.

- enrichedRegions now returns reads per kilobase per million (RPKM), as
  well as raw counts

- Corrected enrichedRegions fold change calculation. Now it is
  (counts1/nsample1)/(counts2/nsample2) instead of
  (counts1/(nsample1-counts1))/(counts2/(nsample2-counts2))

- Reduced substantially the amount of memory required by rangesPlot

Changes in version 0.2.3:

- Fixed issue in enrichedRegions which caused it to return an object
  with too many spaces

- rangesPlot now produces a warning when no ranges are found in the
  specified region, but it still produces a plot

- Added genePlot method for signature 'RangedData'

- Fixed issue in genePlot which caused it to ignore the ... argument.

- Fixed issue in genePlot which sometimes caused it to plot additional
  splicing variants

- Fixed issue in 'filterDuplReads' which caused the returned object to
  have a too large number of spaces

- Added methods for signature(regions='RangedData',
  sample1='RangedData', sample2='missing') and

- Fixed bug in enrichedRegions when no hits were found. Also changed
  argument 'parallelComp' for 'mc.cores', so facilitate using mclapply
  instead of 'parallel' calls

- Added function tabDuplReads

- Fixed bug in filterDuplReads for RangedData signature which caused it
  to filter two reads from different chromosomes which had the same
  start and end positions

- Added extendRanges and filterDuplReads methods for RangedDataList
  objects. Added mc.cores argument to allow for parallel computing.

- Added export2SGR method for RangedDataList objects.

- export2SGR now uses 'coverage' function instead of 'pileup', as the
  latter has been deprecated

- Added export2aligned method for RangedDataList objects. Added
  argument dir.

- alignPeaks now returns an object of the same class as its input,
  instead of always returning IRangesList objects.

- Added alignPeaks method for RangedDataList objects, and parameter
  mc.cores to allow parallel computation.

- Major change in readAlignedBatch. Sequences from all files are now
  returned into a single RangedDataList object. Files are not directly
  saved to disk as in previous versions.

- Fixed bug in rangesPlot which caused it to plot reads longer than
  maxFragLength twice

- Fixed bug in genePlot which caused an error for genes with a single
  exon

- Added PeakLocation function. It produces the same plots as
  stdPeakLocation except that it doesn't standardize the coordinates.

- Added p-value adjustment in enrichedRegions. Defined methods for the
  case when sample2 is missing.

- readAlignedBatch now checks for the existence of the files, and
  returns a warning if they don't.

- Added listOverlap function, which allows to assess the overlap
  between lists e.g. two ChIP-seq experiments, or ChIP-Seq vs
  microarray results.

- Corrected bug in alignPeaks which could cause sequences to start at
  negative indexes.

Changes in version 0.2.1:

- ReadAlignedBatch now returns 'RangedData' objects. This required
  substantial structural changes, as new methods had to be defined and
  many existing ones needed to be adjusted.

- Added methods 'filterDuplReads' and 'mergeRegions' for signature
  'RangedData'

- Changed stdPeakLocation to a generic function and defined methods for
  signatures data.frame and RangedData.

- Added closestGene method for signature 'RangedData'.

- Added enrichedPeaks function to detect peaks in sequencing
  experiments.

- Added enrichedRegions function to detect significantly enriched
  regions in sequencing experiments.

- closestGene modified: After bioC upgrade parameter multiple stoped
  working due to findOverlaps function that deprecated Overlaps
  function (which was used in previous versions).

- closestGene: Some rbinds have been removed and do.call is used insted
  to use less memory.

- closestGene: Better control of downloaded refflat files has been
  added.

- closestGene: new parameter added (promoterInOther). If it's TRUE
  peaks that are in a promoter region and inside other gene will be
  assigned to both genes. Only works when multiple=TRUE.

Changes in version 0.2.0:

- Fixed bug in alignPeaks which caused it not to work when there were
  some chromosomes with no sequences. Also, now alignPeaks reports the
  estimated shift.

- readAlignedBatch deletes memory consuming objects after they have
  been saved to disk. This lowers memory requirements.

- readAlignedBatch forces the chromosome names to start by 'chr', to
  ensure compatibility with external softwares.

- Added vignette

- Added Rd file for alignPeaks

- Added export2aligned, export2SGR and rangesPlot methods for
  RangedData and RangedDataList objects.

- export2SGR now longer requires to specify chromosome length for
  IRangesList and RangedDataList objects. The range of positions where
  reads have been observed is used instead.

Changes in version 0.1.2:

- Added vignette

- Exported rangesPlot in NAMESPACE.

- Added function stdPeakLocation to plot peak location in ChIP-Seq and
  ChIP-chip experiments

- Corrected mismatching arguments in generic definition for genePlot.
  Also did minor changes to improve label appearance.

- closestGene did not work after upgrading to R 2.10.

Changes in version 0.1.1:

- closestGene now uses findOverlaps function instead of the deprecated
  'overlap'.

- Converted rangesPlot function to S4. Added methods for character,
  IRanges and IRangesList.

- Converted genePlot function to S4. Added methods for IRanges,
  IRangesList and CompressedIRangesList objects (additionally to the
  already existing one for character).

- Added alignPeaks function to align the '+' and '-' reads in ChipSeq
  experiments

- Now readAlignedBatch returns a RangeDataList object for single end
  experiments instead of an IRangesList. This allows to store the
  strand information, which could not be done with the IRangesList
  object.

- Added export2SAM function

- Added export2aligned method for 'IRanges' objects

- Added 'paired' option in readAlignedBatch

- Changed default type='Bowtie' in readAlignedBatch

- Added genePlot function

- closestGene had an error when looking if a region falls in a
  promoter.

- Distance to gene added to clsoestGene function.

Changes in version 0.1.0:

- Added readAlignedBatch function

Changes in version 0.0.1:

- extendRanges now outputs an IRangesList object instead of a simple
  list

- Added export2aligned method.

- Added export2SGR method for IRangesList objects

- Added closestGene function.

- Changed "c" for "rbind" in "enrichedRegions.r" line 92 (needed for
  IRanges version 1.6.0)

- Added require(multicore) to "enrichedRegions.r" and changed mclapply
  to lapply for non-parallel computation

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

KEGGgraph
---------

- 2011-08-23 KEGGgraph 1.9.2: Update getKGMLurl and retrieveKGML to
get KGML files from HTTP pages instead of from the FTP server, which
is from July 2011 only open to subscribed users. The KEGG database
needs support, and users of KEGGgraph could also support by
subscribing to the FTP service.

- 2011-08-23 KEGGgraph 1.9.1: Add parseKGML2DataFrame function. This
function is table to convert the KGML file into a data frame recording
nodes and edges in one step. For help use '?parseKGML2DataFrame' in R
command line.

limma
-----

Changes in version 3.10.0:

- New function voom() allows RNA-Seq experiments to be analysed using
  the standard limma pipeline. An RNA-Seq case study is added to User's
  Guide.

- treat(), roast() and mroast() can now estimate and work with a trend
  on the prior variance, bringing them into line with eBayes().

- barcodeplot() and barcodeplot2() merged into one function.

- removeBatchEffect() can now correct for two batch factors.

- plotMDS is now an S3 generic function.  This allows MDS plots to be
  redrawn with new labels without needing to repeat the distance or
  scaling calculations. New S4 class "MDS" to hold the multidimensional
  scaling information output from plotMDS.

- getEAWP() now gets probe annotation from the expression rownames of
  an EList object, if no other probe annotation is available.

- topRomer() now ranks gene sets by secondary columns as well the
  primary criterion specified, to give a more meaningful ranking when
  the p-values are tied.

- wilcoxGST() now accepts signed or unsigned test statistics. Change to
  p-value calculation in geneSetTest() when rank.only=FALSE to avoid
  zero p-values and follow recommendation of Phipson and Smyth (SAGMB,
  2010).

- plotMA() now recognizes ElistRaw and EList objects appropriately.

- Default span for normalizeCyclicLoess increased from 0.4 to 0.7.
  Speed improved when weights=NULL.

- Weaver case study (Section 11.5) in User's Guide is updated and
  rewritten. Data classes ElistRaw and Elist now described in the quick
  start section of the User's Guide.  Other minor updates to User's
  Guide.

- Bug fix for normalizeBetweenArrays() when object is an EListRaw and
  method="cyclicloess".  Previously this function was applying
  cyclicloess to the raw intensities, then logging.  Now it logs first,
  then applies cyclicloess.

- Bug fix to avereps() for EList objects x when x$other is not empty.

lumi
----

- Adding supports for Infinium 450K methylation microarrays,
  which\nincludes updating the color bias correction and other
  related\nfunctions.

minfi
-----

Changes in version 0.99:

- Initial release to Bioconductor.

- Added NEWS file.

- Bugfix to vignette.

- readIDAT is now exported by crlmm, implying that we can import this
  function through NAMESPACE.

MLP
---

Changes in version 1.1.1:

- update using new annotation packages

- apply changes to column permutations (contributed by Heather Turner)

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

Mulcom
------

Changes in version 1.3.0:

- Mulcom subroutines are substantially rewritten in C.

ncdfFlow
--------

0.1.0

FEATURES:

- netCDF support for large data sets.

- support all the related methods for flowSet
	 
KNOWN ISSUES:
	
- write meta data - ncdfFlow allows user to save the entire
  ncdfFlowSet object in ncdf file.  Currently the meta data is first
  serialized in R and stored as raw vector in cdf. It can fail when
  the meta data size exceeds the limit of serialization function.

netresponse
-----------

Changes in version 1.5.14 (2011-10-28):

- removed compiler package and functions as they are included in R-2.14

Changes in version 1.5.12 (2011-06-29):

- added options for read.network

Changes in version 1.5.11 (2011-06-28):

- compiler added

- corrected weight estimation from vdp.mixt

Changes in version 1.5.09 (2011-06-13):

- updated toydata

- compute.weights: now operating in log-domain to avoid occasional
  floating point errors

Changes in version 1.5.08 (2011-05-27):

- added plot.pca function to visualize subnets in 2D with ellipsoid
  confidence intervals and annotation colors

Changes in version 1.5.04 (2011-05-27):

BUG FIXES

- fixed discretization NAMESPACE issues

Changes in version 1.5.02 (2011-05-18):

SIGNIFICANT USER-VISIBLE CHANGES

- added speedup option speedup.max.edges in detect.responses. To
  consider only this many most similar neighborghs in merging,
  similarity is evaluated by empirical mutual information estimate with
  sqrt(n) bins where n is data sample size

- changed AIC to default information criterion (previously BIC). Both
  remain to be available through information.criterion argument

pdInfoBuilder
-------------

Changes in version 1.18:

USER VISIBLE CHANGES

- The annotation packages for SNP chips no longer have a
  'fragment_length*' column in the featureSet/featureSetCNV tables.
  Fragment length now has its own table, 'fragmentLength' (for SNP
  probes) and 'fragmentLengthCNV' (for CN probes).

- Rda files available through annotation packages created by
  pdInfoBuilder > 1.17 are compressed using 'xz'.

- Fixed License field on the template.

- Added an additional template for SNP arrays supported by CRLMM.

pint
----

Changes in version 1.5.35 (2011-10-13):

- Fixed bugs with GenomeModels accessors

Changes in version 1.5.34 (2011-08-15):

- changed the default match.probes to TRUE in screen.cgh.mrna; now also
  checks with match.probes = FALSE that the probes are matched

Changes in version 1.5.33 (2011-08-02):

- now artificially adding arm information if it is missing, to avoid
  halts

Changes in version 1.5.30 (2011-07-15):

- screen.cgh.mrna option segmented changed into match.probes

Changes in version 1.5.29 (2011-05-09):

BUG FIXES

- Fixed a bug in plot.ChromosomeModels

Changes in version 1.5.28 (2011-05-09):

BUG FIXES

- fixed a bug which prevented screeningn while missing values were
  present

Changes in version 1.5.06 (2011-04-16):

- added interpretation tools in report.R

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

RPA
---

Changes in version 1.9.19 (2011-10-29):

- removed compiler and corrected the broken compiler function
  references

Changes in version 1.9.17 (2011-06-27):

- added bg correction in rpa.online

- accurate and relatively fast variance hyperparameter update function
  added (see update.hyperparameters, update.s2)

Changes in version 1.9.11 (2011-06-02):

- speedups in hyperparameter estimation; alpha treated as scalar

- removed affinity.method

- added quantiles.online normalization method BUG FIXES

- alpha, beta updates fixed

- rpa.fit did not return updated alpha, beta: now corrected.

Changes in version 1.9.06 (2011-05-29):

- rpa.fit-class: added alpha, beta prior parameters for the inverse
  Gamma conjugate prior for the probe-specific variances

- modified RPA.iteration, RPA.update.sigma2 so as to directly utilize
  priors everywhere.

- added set.alpha, set.beta and update.hyperparameters, update.alpha,
  update.beta in internal functions

- added toydata generator function sample.probeset

- rpa.plot: added comparison of toydata and fitted data by adding the
  toydata.object argument

Changes in version 1.9.03 (2011-05-14):

- corrected AffyCompII result note

- polished plot functions in rpa.plot and plot-methods

Changes in version 1.9.02 (2011-04-21):

- added compiler to betahat in RPA.sigma2.update.R to speed up

- added NA/NaN check to RPA.iteration after a defected affybatch with
  NA values was found to cause crash

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

TargetSearch
------------

Changes in version 1.10.0:

NEW FEATURES

- New function "writeLibText" to save "tsLib" objects in a
  tab-delimited text file (which can consecuently be re-imported). It
  might be also used to convert a MSP file into a TEXT file.

BUG FIXES

- Added some 'drop=FALSE' in plotSpectra function.

- Calls to 'pdf' are now followed by 'on.exit(dev.off())' to make sure
  that the device is closed after unexpected errors or user
  interruption.

Changes in version 1.8.2:

BUG FIXES

- bug fixed error due to changes in TargetSearch object definitions in
  function writeMSP.

Changes in version 1.8.1:

BUG FIXES

- 'sampleRI' and 'Profile' functions would return NA if the minPairCor
  parameter was greater than the number of samples.

- Change GUI message to a more suitable one.

- 'FAMEoutliers' failed if only one sample was analyzed.

- 'ProfileCleanUP' accepts a minPairObs parameter like 'Profile'

- The parameter minPairObs is now checked so that cannot take a value
  lower than 5 in 'sampleRI', 'Profile', and 'ProfileCleanUP'
  functions.

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

tigre
-----

Changes in version 1.8.0:

NEW FEATURES

- Direct support for variances for Illumina Bead Array data from lumi.

- Prettier default plot style

BUG FIXES

- Fix to variance scaling in plots of models with TF mRNA measurements

Changes in version 1.6.1 (2011-04-26):

BUG FIXES

- Fix: target ranking now works with >= 10 known targets (thanks to
  Jake Michaelson for the report)

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

Changes in version 1.14.0:

VERSION xps-1.13.10

- update function xpsQAReport() to replace "_" with "\_" for *.Rnw

VERSION xps-1.13.9

- add quality report function xpsQAReport()

- update vignette xps.Rnw

VERSION xps-1.13.8

- add functions attachProbe(), removeProbe(), etc.

- add functions contentGC(), probeSequence()

- add functions inten2GCplot(), plotInten2GC()

- update vignette xps.Rnw

VERSION xps-1.13.7

- add function unitID2symbol()

- add function plotProbeset()

- update vignette xps.Rnw

VERSION xps-1.13.6

- update function indexUnits() for exon probesets

- add function symbol2unitID()

- add function probesetplot()

VERSION xps-1.13.5

- update script4schemes.R to include schemes with annotation na32

- update script4xps.R and script4exon.R

VERSION xps-1.13.4

- move to ROOT version 5.30/00, update README

- update vignettes xps.Rnw and xpsClasses.Rnw

VERSION xps-1.13.3

- correct problem in unitID2affyID() on Linux caused by sQuote()

VERSION xps-1.13.2

- add functions plotXXX()

- function hist() no longer requires to attachInten()

- add functions indexUnits(), pmindex(), mmindex()

- add functions probesetID2unitID(), unitID2probesetID(), etc

- add functions attachUnitNames(), removeUnitNames()

- add functions attachDataXY(), removeDataXY()

VERSION xps-1.13.1

- update function plotImage() to draw background images


Packages removed from the release
=================================

The following packages are no longer in the release: GeneSpring,
GeneTraffic, makePlatformDesign, Rdbi, RdbiPgSQL, rflowcyt, Rredland,
Ruuid, simulatorAPMS.

[bioc]: http://bioconductor.org
[install]: http://bioconductor.org/install/
[ami]: http://bioconductor.org/help/bioconductor-cloud-ami/
