April 2, 2012

Bioconductors:

We are pleased to announce Bioconductor 2.10, consisting of 554
software packages and more than 600 up-to-date annotation packages.
There are 45 new software packages, and many updates and improvements
to existing packages; 5 packages have been removed from this
release. Bioconductor 2.10 is compatible with R 2.15.0, and is
supported on Linux, 32- and 64-bit Windows, and Mac OS.  This release
includes an updated Bioconductor [Amazon Machine Image][ami].  Visit
[http://bioconductor.org][bioc] for details and downloads.

Contents
--------

* Getting Started with Bioconductor 2.10
* New Software Packages
* NEWS from new and existing packages
* Packages removed from the release

Getting Started with Bioconductor 2.10
======================================

To install Bioconductor 2.10:

1. Install R 2.15.0.  Bioconductor 2.10 has been designed expressly for
this version of R.

2. Follow the instructions at [http://bioconductor.org/install/][install].

New Software Packages
=====================

There are 45 new packages in this release of Bioconductor.

- AffyRNADegradation: Analyze and correct probe positional bias in
  microarray data due to RNA degradation

- ASEB: Predict Acetylated Lysine Sites

- BiocGenerics: Generic functions for Bioconductor

- birta: Bayesian Inference of Regulation of Transcriptional Activity

- BitSeq: Transcript expression inference and differential expression

- BRAIN: Baffling Recursive Algorithm for Isotope distributioN calculations

- BrainStars: query gene expression data and plots from BrainStars (B*)

- CancerMutationAnalysis: Cancer mutation analysis

- categoryCompare: Meta-analysis of high-throughput experiments using feature
  annotations

- cellGrowth: Fitting cell population growth models

- cnvGSA: Gene Set Analysis of (Rare) Copy Number Variants

- coGPS: cancer outlier Gene Profile Sets

- DART: Denoising Algorithm based on Relevance network Topology

- deepSNV: Test for subclonal SNVs in deep sequencing experiments.

- easyRNASeq: Count summarization and normalization for RNA-Seq data.

- EBcoexpress: EBcoexpress for Differential Co-Expression Analysis

- ffpe: Quality assessment and control for FFPE microarray expression

- GeneGroupAnalysis: Gene Functional Class Analysis

- GEWIST: Gene Environment Wide Interaction Search Threshold 

- gprege: Gaussian Process Ranking and Estimation of Gene Expression
  time-series

- Gviz: Plotting data and annotation information along genomic coordinates

- gwascat: representing and modeling data in the NHGRI GWAS catalog

- HiTC: High Throughput Chromosome Conformation Capture analysis

- HybridMTest: Hybrid Multiple Testing

- iASeq: iASeq: integrating multiple sequencing datasets for detecting
  allele-specific events

- iBBiG: Iterative Binary Biclustering of Genesets

- IdMappingAnalysis: ID Mapping Analysis

- inSilicoMerging: Collection of Merging Techniques for Gene Expression Data

- manta: Microbial Assemblage Normalized Transcript Analysis

- maskBAD: Masking probes with binding affinity differences

- MinimumDistance: A package for de novo CNV detection in case-parent trios

- motifRG: A package for discriminative motif discovery, designed for high
  throughput sequencing dataset

- NarrowPeaks: Functional Principal Component Analysis to Narrow Down
  Transcription Factor Binding Site Candidates
  
- pcaGoPromoter: pcaGoPromoter is used to analyze DNA micro array data

- phyloseq: Handling and analysis of high-throughput
  phylogenetic sequence data.
  
- PING: Probabilistic inference for Nucleosome Positioning with
  MNase-based or Sonicated Short-read Data
  
- QUALIFIER: Qualitiy Control of Gated Flow Cytometry Experiments

- RchyOptimyx: Optimyzed Cellular Hierarchies for Flow Cytometry

- ReactomePA: Reactome Pathway Analysis

- rhdf5: HDF5 interface to R

- sigaR: statistics for integrative genomics analyses in R

- spade: SPADE -- An analysis and visualization tool for Flow Cytometry

- ternarynet: Ternary Network Estimation

- VegaMC: VegaMC: A Package Implementing a Variational Piecewise Smooth
  Model for Identification of Driver Chromosomal Imbalances in Cancer

- virtualArray: Build virtual array from different microarray platforms

NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages. The following package NEWS is available:



a4Base
------

Changes in version 1.2.4:

- throw error in case 'group' variable contains missing values in
  spectalMap

Changes in version 1.2.3:

- circumvent build issues on certain versions of windows

Changes in version 1.2.2:

- more appropriate handling of the dots argument in plot1gene

Changes in version 1.2.1:

- fix in legend for spectralMap

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

affxparser
----------

Changes in version 1.28.0 (2012-03-30):

- The version number was bumped for the Bioconductor 2.10 release
  version.

Changes in version 1.27.5 (2012-03-19):

- ROBUSTNESS: Now the native code throws R errors, instead of printing
  an error message to stdout/stderr and then returning null, that is
  translated to errors at the R level.

Changes in version 1.27.4 (2012-03-05):

- BUG FIX: affxparser would not build on Windows with the new Rtools
  toolchain (Rtools 2.15.0.1915-1919).

Changes in version 1.27.3 (2011-11-18):

- ROBUSTNESS: Added sanity checks asserting that the internal readers
  did indeed read something and not just returned NULL. It should be
  very unlikely that this occurs, but there is still a small risk that
  after asserting that a file exists, but before the internal Fusion
  SDK parsers access the file, the file has been removed.

Changes in version 1.27.1 (2011-11-01):

- Same updates as in v1.26.1.

Changes in version 1.27.0 (2011-10-31):

- The version number was bumped for the Bioconductor devel version.

Changes in version 1.26.4 (2012-03-06):

- BUG FIX: affxparser would not build on Windows with the new Rtools
  toolchain (Rtools 2.15.0.1915-1919), which is for R (> 2.14.1), i.e.
  also for R v2.14.2 (but not v2.14.1). This is the same bug fix that
  was first done in v1.27.4.

Changes in version 1.26.2 (2011-11-16):

- The version number was bumped by Bioconductor to trigger a build.

Changes in version 1.26.1 (2011-11-01):

- FIX: Fixed warning on "In readBin(con, what = "integer", size = 4, n
  = 1, signed = FALSE, 'signed = FALSE' is only valid for integers of
  sizes 1 and 2" that some read methods would generated.

aroma.light
-----------

Changes in version 1.24.0 (2012-03-30):

- The version number was bumped for the Bioconductor release version.

Changes in version 1.23.0 (2011-10-31):

- The version number was bumped for the Bioconductor devel version.

ArrayExpressHTS
---------------

Changes in version 1.5.2:

UPDATES

- fixed trace names, cleanup

- updated the code downloading SDRFs to rename seq.srdf.txt to
  .sdrf.txt

- unified the functions providing default options

- exposed the options that are used in the process

- changed errors to warnings where process can still be completed

- fixed improper tools detection on the external R Cloud, added
  detection using elements in PATH.

- simplified setting of tools options, removed the need to invoke
  initEnvironmentVariables after the options are set, made it automatic

- updated help pages

- composed all options related pages into "package options" page and
  "processing options" page. Exposed all options.

- removed the initEnvironmentVariables

- removed outdated references

Changes in version 1.5.1:

UPDATES

- improved error handling within the ArrayExpressHTS

- added stop.on.warnings that allows detection of possible failures
  earlier. However if used it would narrow the amount of experiments
  that can be successfully processed.

- added log.error instead of log.info where necessary

- reorganised pipeline options, made them visible to users in the
  command syntax.

- reorganised R Cloud parameters into a single set of visible/usable
  options

- made the pipeline automatically detect "Organism" from SDRF and
  filter using a selected value.

- reworked reading from SDRF to cover most of the "flexible formatting"
  cases

- improved creation of R Cloud computation cluster, widened usage of
  retries.

- added proper functions to work with environment variables

- updated package help pages

- updated the annotation processing following Ensembl update.

- updated help pages, descriptions and references

- updated package vignette documentation

- slightly reorganised the structure of the document

- rewritten some sections

- added a few new sections

- added examples where needed

- reimplemented problematic searching and downloading of .seq.sdrf.txt
  files

- fixed problematic handling of "NA" or "missing" nominal length values
  from ENA

- made a number of reference preparation tasks be able to run in
  parallel and not affect each other.

- fixed the "count" method to properly read fields of bam files

biovizBase
----------

Changes in version 1.1.16:

NEW FEATURES

- addSteping add extra arugments group.selfish to control stepping
  method

- add estimateCoverage method for fast estimate coverage

SIGNIFICANT USER-LEVEL CHANGES

- addStepings rename to addStepping

CAMERA
------

Changes in version 2012-03-19:

- Version 1.11.10

- Bugfix in groupCorr, where the function throws an error if argument
  xraw is not null

- Version 1.11.9

- Bugfix in findAdducts, where adduct annotation filtering was to
  stringent, if ips score is higher 1.5

Changes in version 2012-02-29:

- Version 1.11.8

- Bugfix in findAdducts parallel mode, where wrong psgrp indices were
  stored in annoGrp

Changes in version 2012-02-17:

- Version 1.11.7

- Bugfix in findAdducts, where the ruleset is not saved within a
  xsAnnotate object with user defined rules

- Changed a general groupCorr behaviour. If no edge is above
  correlation threshold, all peaks are seperated in pspecs of size 1

Changes in version 2012-02-15:

- Version 1.11.6

- Fix error in groupCorr: If sample is set to something other than 1 or
  NA, it result in a crash

Changes in version 2012-02-10:

- Version 1.11.5

- Bugfix in groupCorr where parameter cor_eic_th doesn't influence
  correlation across samples

- Add new parameter for groupCorr cor_exp_th

Changes in version 2011-29-11:

- Version 1.11.3

- Add new function combinexsAnnos. It allows checking and reannotation
  of sample with a coressponding sample from the opposite ion mode

Changes in version 2011-28-04:

- Fix bug in getpspectra, which incorrect label isotope peaks

Changes in version 2011-24-11:

- Version 1.11.2

- Bugfix in annotateDiffreport: Since 1.7.7 mismatch of the peaklist
  from CAMERA and the diffreport function results in false ordered
  peaktable

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

Changes in version 2011-12-12:

- Version 1.11.4

- Add missing Rd page combinexsAnnos

Changes in version 2011-12-05:

- plotPsSpectrum now accepts additional parameters for plot

Changes in version 2011-10-11:

- Version 1.11.1

- add the possibility to extract multiple isotope intensity from
  different samples with getIsotopeCluster

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

Changes in version 2010-29-09:

- Add function findNeutralLoss and findNeutralLossSpecs

Changes in version 2010-23-11:

- Fix bug in findAdducts. Function failed, if run in parallel mode.

Changes in version 2010-20-09:

- Add function getIsotopeCluster for retrival of isotope cluster

Changes in version 2010-11-10:

- Rewrite Vignette

Changes in version 2010-11-06:

- Speed up in findAdducts

- Fix unit tests

Changes in version 2010-08-11:

- Add intval parameter to groupFWHM and findIsotopes

- Change getPeaklist to S4 Methods.

- Add getPeaklist parameter intval, where the intensity value can be
  selected

- Additional bugfix for findIsotopes. Could occur with 1.7.1, that all
  isotopes will be deleted.

Changes in version 2010-05-03:

- Fix bug in findIsotopes, crashed with dim error

Changes in version 2010-04-20:

- Fix bug in annotate function

Changes in version 2010-04-19:

- Last changes for the next BioC release

- Rewrite of the vignette

Changes in version 2010-03-17:

- Hugh changes in CAMERA for working with multiple sample

- Add check constrains in groupCorr for isotopes and primary adducts

Changes in version 2010-01-11:

- Fix bug in findIsotopes. Occurs if the first isotope peak, could be
  assigned from two different monoisotopic peaks.

- Reduce slightly the number of found isotopes

Changes in version 2009-11-24:

- Add experimental parallel mode with MPI

Changes in version 2009-11-20:

- Fix bug in setting ruletable

Changes in version 2009-08-26:

- add nSlaves as argument to findAdducts (for parallel annotation with
  MPI)

Changes in version 2009-08-24:

- small bugfixes in getPeaklist

- add neutral losses to rule set

Changes in version 2009-08-12:

- Fixed plotEICs for 1 peak groups

Changes in version 2009-08-04:

- Add ips Score to annoGrp

Changes in version 2009-06-03:

- Add adduct labels to plotPeaks()

- Add adduct labels to plotEICs()

Changes in version 2009-05-22:

- Add first visualisation function plotEICs() and plotPeaks()

Changes in version 2009-05-06:

- bump for devel 2.5

Changes in version 2009-03-30:

- add combine_xsanno

- small refactoring

- changes in scoring schemata

Changes in version 2009-03-18:

- after groupCorr every peak is now member of a group

- bugfix: remove xM-xH clones

- remove dependency on Hmisc

Changes in version 2009-03-10:

- speed-up of findAdducts

Changes in version 2009-02-24:

- refactoring findAdducts

- add methods for pos/neg polarity comparison

Changes in version 2009-02-20:

- replace na.omit with naOmit

Changes in version 2009-02-18:

- add neutral losses into ruleset

- other bugfixes

Changes in version 2009-01-19:

- Change isotope nomination

- add lists for fragments, ions and neutral losses

Changes in version 2008-10-15:

- Lot of bugfixes and speed up the correlation

Changes in version 2008-10-13:

- First build for version 0.1.1

Changes in version 2007-10-12:

- Inital release & collection of stuff

categoryCompare
---------------

Changes in version 0.99.0:

SIGNIFICANT USER-VISIBLE CHANGES

- First review in Bioconductor

Changes in version 0.7.8:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated the vignette. Removed long lines, added links to figures,
  fixed some errors

Changes in version 0.7.4:

SIGNIFICANT USER-VISIBLE CHANGES

- changed category in ccCompare to categoryName for
  GENccEnrichResults-class

Changes in version 0.7.3:

SIGNIFICANT USER-VISIBLE CHANGES

- First commit

cellGrowth
----------

Changes in version 0.99.0:

SIGNIFICANT USER-VISIBLE CHANGES

- no dontrun in examples

- use system.file() instead of installed.packages() in examples

- renamed the example_files directory to extdata

Changes in version 0.4.0:

SIGNIFICANT USER-VISIBLE CHANGES

- first submission

clusterProfiler
---------------

Changes in version 1.4.0:

- bump up version to 1.4.0 for 2.10 release

- add URL in DESCRIPTION

Changes in version 1.3.15:

- add support to organsims other than human, mouse, and yeast.
  <2012-03-22, Thu>

- add GFFparser.R, which provide function Gff2GeneTable to parse Gff
  files and build gene information table.<2012-03-22, Thu>

- add function buildGOmap for building GO mapping files for unsurported
  organisms, which can be analyzed now. <2012-03-22, Thu>

- add vignette explaining how to run GO analysis for unsupported
  organism. <2012-03-12, Thu>

Changes in version 1.3.14:

- re-implement enrichGO by import enrich.internal in DOSE implement S3
  methods for mapping IDs for GO analysis. <2012-03-18, Sun>

- re-implement enrichKEGG by import enrich.internal in DOSE implement
  S3 methods for mapping IDs for KEGG analysis. <2012-03-19, Mon>

- re-implement groupGOResult class by extended from enrichResult, and
  modified groupGO by using S3 methods designed for enrichGO.
  <2012-03-19,Mon>

- bug fixed for importing setReadable<- instead of setReadable.
  <2012-03-20, Tue>

- bug fixed of TERM2NAME.GO. <2012-03-21, Wed>

Changes in version 1.3.13:

- import plot summary from stats4, for BiocGenerics (version 0.1.10)
  removed them <2012-03-03, Sat>

Changes in version 1.3.12:

- fixed BibTeX database file .bib. month = , must be month = someMonth,
  or totally deleted, leave it blank will cause texi2dvi failed.
  <2012-03-01, Thu>

Changes in version 1.3.11:

- update vignette. <2012-02-28, Tue>

Changes in version 1.3.10:

- fixed warnings concerning documents of plot generics. <2012-02-27,
  Mon>

- import summary generic from BiocGenerics instead of stats4.
  <2012-02-27, Mon>

Changes in version 1.3.9:

- fixed build error in bioc repos <2012-02-26 Sun>

Changes in version 1.3.8:

- add visualization section in vignette <2012-02-16 Thu>

Changes in version 1.3.7:

- @exportMethod plot <2012-02-15 Wed>

- fix bug when calling summary method from plot, for summary defined in
  base is S3 method, instead import summary generic from stats4
  <2012-02-15 Wed>

Changes in version 1.3.6:

- remove generic definition of show, summary and plot, add NAMESPACE
  import show from methods and plot from graphics. summary need not to
  import, for is defined in the base package. <2012-02-13 Mon>

- rewrite vignette <2012-02-13 Mon>

Changes in version 1.3.5:

- update plot codes accompany with new version of ggplot2(>=0.9.0)
  <2012-02-07 Tue>

Changes in version 1.3.4:

- using apply instead of mdply to improve speed <2012-02-03 Fri>

- add citation of clusterProfiler <2012-02-03 Fri>

Changes in version 1.3.3:

- add showCategory parameter for plot functions <2011-01-11 Wed>

Changes in version 1.3.2:

- add qvalueCutoff parameter. <2011-01-05 Thu>

Changes in version 1.3.1:

- add plot.categoryNet for enrichGOResult. <2011-11-10 Thu>

CNAnorm
-------

Changes in version 1.1.8:

SIGNIFICANT USER-VISIBLE CHANGES

- Fixed plotGenome so that, when fixVaxes = TRUE horizontal lines are
  drawn even if there are no points in a certain area.

Changes in version 1.1.7:

SIGNIFICANT USER-VISIBLE CHANGES

- Removed arguments 'adjust' and 'n' from plotPeaks as it needs to
  retrieve the same values used in peakPloidy

- Added option fixVAxes to plotGenome to fix vertical axes to maxRatio
  and minRatio

- Fixed plotGenome so that the line from superimpose = "DNACopy" is
  actually a line and not a series of dots

- If plotting only one chromosome, instead of plotting chr name, it
  plots chromosome position.

Changes in version 1.1.6:

SIGNIFICANT USER-VISIBLE CHANGES

- Added method 'closest' to function peakPloidy for a "standard"
  normalisation

- Fixed the 'mode' method

- Internal tiding

Changes in version 1.1.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Added methods 'median' and 'mode' to to function peakPloidy for a
  "standard" normalisation

Changes in version 1.1.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Fixed an error in exportTable which was providing a meaningless value
  for SegMean

Changes in version 1.1.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Changed the provided data (LS041 and CN) to match what bam2windows.pl
  version 0.3.4

Changes in version 1.1.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Added an option to plot more colorful genome plots

Changes in version 1.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Corrected a bug in exportTable that crashed if smoothed ratio was not
  available.

codelink
--------

Changes in version 1.23.2:

- Improvements on how codPlotScatter() manages labels. Renamed generic
  scatter plot function plotma() to plotxy(), which better reflects the
  current status.

Changes in version 1.23.1:

- Implemented missing function codPlotScatter() for CodelinkSet objects
  (accessible via codPlot(..., what = "scatter")). The implementation
  reuses the plotma() function used by codPlotMA(). By default the
  first array is plotted against the median. Thanks to William Michels
  for reporting this and testing the CodelinkSet framework.

- Minor improvements to linear model section in CodelinkSet-vignette.
  Now the code is evaluated.

cqn
---

Changes in version 1.1:

- Same fixes as in 1.0.1

- Resaved the data files, so they take up less space.

- Added edgeR as a Suggests: since the vignette uses it.

Changes in version 1.0.1:

- The function alpha has been moved from ggplot2 to the new package
  scales.  Vignette and Suggests: fields have been changed accordingly.

- sizeFactors = NULL should now work for cqn and cq.fixedlength. Thanks
  to Maria Chikina for reporting this

cummeRbund
----------

Changes in version 1.1.5:

BUG FIXES

- Fixed minor bug in database setup that caused instability with
  cuffdiff --no-diff argument. oFixed bug in csDendro method for
  CuffData objects.

Changes in version 1.1.4:

NEW FEATURES

- Added MAplot() method to CuffData objects.

BUG FIXES

- Finished abrupt migration to reshape2. As a result fixed a bug in
  which 'cast' was still required for several functions and could not
  be found. Now appropriately using 'dcast' or 'acast'.

- Fixed minor bug in CuffFeature::fpkmMatrix

Changes in version 1.1.3:

NEW FEATURES

- getSig() has been split into two functions: getSig() now returns a
  vector of ids (no longer a list of vectors), and getSigTable()
  returns a 'testTable' of

- binary values indicating whether or not a gene was significant in a
  particular comparison.

- Added ability in getSig() to limit retrieval of significant genes to
  two provided conditions (arguments x & y). (reduces time for function
  call if you have a specific comparison in mind a priori)

- When you specify x & y with getSig(), q-values are recalculated from
  just those selected tests to reduce impact of multiple testing
  correction.

- If you do not specificy x & y getSig() will return a vector of
  tracking_ids for all comparisons (with appropriate MTC).

- You can now specify an 'alpha' for getSig() and getSigTable() [ 0.05
  by default to match cuffdiff default ] by which to filter the
  resulting significance calls.

- Added csSpecificity() method: This method returns a
  feature-X-condition matrix (same shape as fpkmMatrix) that provides a
  'condition-specificity' score defined as 1-(JSdist(p,q)) where p is
  is the density of expression (probability vector of log(FPKM+1)) of a
  given gene across all conditions, and q is the unit vector for that
  condition (ie. perfect expression in that particular condition)

- specificity = 1.0 if the feature is expressed exclusively in that
  condition

- Created csDendro() method: This method returns a object of class
  'dendrogram' (and plots using grid) of JS distances between
  conditions for all genes in a CuffData, CuffGeneSet, or
  CuffFeatureSet object.

- Useful for identifying relationships between conditions for subsets
  of features

- New visual cues in several plot types that indicates the
  quantification status ('quant_stat' field) of a particular
  gene:condition. This information is useful to indicate whether or not
  to trust the expression values for a given gene under a specific
  condition, and may provide insight into outlier expression values.
  This feature can be disabled by setting showStatus=F.

- csDensity() is now available for CuffFeatureSet and CuffGeneSet
  objects

BUG FIXES

- Fixed bug in getGenes that may have resulted in long query lag for
  retrieving promoter diffData. As a result all calls to getGenes
  should be significantly faster.

- CuffData fpkm argument 'features' now returns appropriate data.frame
  (includes previously un-reported data fields).

- Replaced all instances of 'ln_fold_change' with the actual
  'log2_fold_change'.  Values were previously log2 fold change but
  database headers were not updated to reflect this.

- Fixed bug that could cause readCufflinks() to die with error when
  using reshape2::melt instead of reshape::melt.

NOTES

- ***The structure of the underlying database has changed in this
  version.  As a consequence, you must rebuild you cuffData.db file to
  use new version. readCufflinks(rebuild=T)***

- Updated vignette

- 'fullnames' logical argument was added to fpkmMatrix. If True,
  rownames for fpkmMatrix will be a concatenation of gene_short_name
  and tracking_id. This has the added benefit of making row labels in
  csHeatmap easier to read, as well as preserving uniqueness.

- Slight speed improvements to JSdist (noticeable when using csCluster
  on large feature sets).

- 'testTable' argument to getSig() has been dropped in lieu of new
  getSigTable() method.

Changes in version 1.1.1:

BUG FIXES

- fixed issue in which there was no graceful error handling of missing
  CDS or TSS data in cuffdiff output.

- Fixed issue in which distribution test data (promoters, splicing,
  relCDS) were not appropriately added to objects on creation. oFixed
  bug that would sometimes cause csBoxplot() to throw an error when
  log-transforming fpkm data. Also added pseudocount argument.

- Fixed bug that would cause diffData() to return a filtered subset of
  results by default.

- Adjusted indexing of tables to improve performance on large datasets.

- Fixed bug that caused diffData method to not be registered with
  CuffFeature and CuffGene objects.

- Fixed bug that sometimes caused over-plotting of axis labels in
  csBarplots.

NEW FEATURES

- added getSig method to CuffSet class for rapid retrieval of
  significant features from all pairwise tests (as a list of IDs).

- By default the level is 'genes' but any feature level can be queried.

- csCluster now uses Jensen-Shannon distance by default (as opposed to
  Euclidean)

- Added 'xlimits' argument to csVolcano to constrain plot dimensions.

- Enforced requirement in csVolcano for x and y arguments (as sample
  names).

NOTES

- Changed dependency 'reshape' to 'reshape2'

- Changed the default orientation of expressionBarplot() for
  CuffFeatureSet objects.

- Changed output of csCluster to a list format that includes clustering
  information. As a result, I created the function csClusterPlot to
  replace the previous default drawing behavior of csCluster.  This
  allows for stable cluster analysis.

- For consistency, the 'testId' slot for CuffDist objects was renamed
  to 'idField'.  This brings the CuffDist class in line with the
  CuffData class.

- CuffGene and CuffGeneSet now include slots for promoter, splicing,
  and relCDS distribution test results.

deepSNV
-------

Changes in version 1.1.4:

BUG FIXES

- Fixed error in summary() when there were no significant SNVs.

- Some fixes if only a single column of the alignment is selected

Changes in version 1.0.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Added CITATION file

- Made NEWS (this file) R-readable

- Changed Vignette

Changes in version 0.99.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Changed plot to S3 method (to avoid warning in R-devel)

Changes in version 0.99.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Added small .bam example files test.bam, control.bam with 100
  positions.

- Modified man pages for bam2R()

- Modified man page for coordinates()

- Corrected example of consensusSequence()

- Compressed .RData files with tools::resaveRdaFiles

- Changed vignette to attach data, rather than load remotely.

- Argument "regions" of deepSNV can be a GRanges object.

Changes in version 0.99.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Added BiocViews field

- Added HIVmix data

- Added new examples

- Registered bam2R with R_registerRoutines

- New accessor functions "test", "control", "p.val", and "coordinates"

- Updated vignette

Changes in version 0.9.5:

SIGNIFICANT USER-VISIBLE CHANGES

- "summary" now reports additional columns from regions slot.

BUG FIXES

- drop=FALSE in subsetting and summary.

Changes in version 0.9.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Directly link to static samtools library provided by Rsamtools

- Load example .bam files over http

Changes in version 0.9.3:

SIGNIFICANT USER-VISIBLE CHANGES

- Added beta-binomial model

- Extended documentation

- Use summary instead of significantSNV

Changes in version 0.9.2:

BUG FIXES

- Minor bugfixes

Changes in version 0.9.1:

BUG FIXES

- Minor bugfixes

Changes in version 0.9.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Pre-release

DESeq
-----

Changes in version 1.7.10 (2011-03-26):

- Added pooled-CR method.

- Fixed error message for insufficient replication.

Changes in version 1.7.9 (2011-03-16):

- Fixed error in vignette.

- Silenced warnings from nbinomFitGLM about convergence; this
  information is given in the result anyway.

Changes in version 1.7.8 (2011-03-15):

- The VST is now roughly equal to log2 for local dispersion fits, too.

Changes in version 1.7.7 (2011-03-15):

- changed the variance stabilizing transformation once more; this time
  it is hopefully correct. The derivation of the formula is now
  documented in the new file vst.pdf.

Changes in version 1.7.6:

- fixed regression bug in getVarianceStabilizedData: the formula for
  parametric dispersion fit was wrong since I "fixed" the factor to get
  log2 asymptotic behaviour.

Changes in version 1.7.1:

- various changes to vignette by whuber

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

DiffBind
--------

Changes in version 1.2.0 (2012-03-30):

- GRanges is default class for peaksets and reports instead of
  RangedData, controlled by DataType parameter.

- Both analysis methods (edgeR and DESeq) use generalized linear models
  (GLMs) for two-group contrasts by default.

- Blocking factors (for two-factor analysis) can be specified flexibly
  such that arbitrary blocking factors can be used.  Section added to
  vignette showing an ananalysis using a blocking factor.

- Added new metadata type, DBA_TREATMENT.

- New DBA_SCORE_ options for specifying scoring method, including TMM
  normalized counts, and ability to change scoring method on the fly in
  dba.plotHeatmap and dba.plotPCA when plotting global binding matrix.

- bRemoveDuplicates parameter in dba.count allows duplicate reads to be
  discarded when computing counts

- More efficient use of memory when analyzing (controlled by
  bReduceObjects parameter in dba.analyze).

- various bugs fixed, man pages updated, and warning messages added.

DOSE
----

Changes in version 1.1.12:

- implement barplot for enrichResult <2012-03-18, Sun>

- bug fixed for setReadable method <2012-03-19, Mon>

- add logFC parameter for cnet plot, support color gene nodes by their
  expression value (log fold change) <2012-03-21, Wed>

- add mapping entrezgene ID and gene Name for organisms other than
  human, mouse and yeast. <2012-03-22, Thu>

- bug fixed for attempt to name logFC, when it is NULL. <2012-03-22,
  Thu>

- optimized readable method. <2012-03-26, Mon>

Changes in version 1.1.11:

- setReadable method for mapping gene ID to gene Symbol in enrichResult
  instance. <2012-03-12, Mon>

- export method show. <2012-03-12, Mon>

Changes in version 1.1.10:

- import plot summary from stats4, for BiocGenerics (version 0.1.10)
  removed them <2012-03-03, Sat>

- Add DO2ALLEG and EG2ALLDO, for mapping undirecte annotation.
  <2012-03-03, Sat>

- update vignette <2012-03-06, Tue>

Changes in version 1.1.9:

- fixed BibTeX database file .bib. month = , must be month = someMonth,
  or totally deleted, leave it blank will cause texi2dvi failed.
  <2012-03-01, Thu>

- update IC data and DO-EG mapping data. <2012-03-01, Thu>

Changes in version 1.1.8:

- update vignette, add semantic similarity algorithms' details.
  <2012-02-28, Tue>

Changes in version 1.1.7:

- fixed warnings concerning documents of plot generics. <2012-02-27,
  Mon>

- import summary generic from BiocGenerics instead of stats4.
  <2012-02-27, Mon>

Changes in version 1.1.6:

- defined S3 generic for ALLEXTID, EXTID2NAME, EXTID2TERMID, TERM2NAME,
  and TERMID2EXTID. <2012-02-26, Sun>

- update roxygen and regenerate man file. <2012-02-26, Sun>

- import S4 generics of plot from BiocGenerics.  <2012-02-26, Sun>

Changes in version 1.1.5:

- add S4 method of plot, which accept parameter type = "cnet", and call
  cnetplot.enrichResult method. <2012-02-23, Thu>

- add S3 method cnetplot.enrichResult for plotting enrichResult object.
  <2012-02-23, Thu>

- define cnetplot function for category-gene network visualization.
  <2012-02-23, Thu>

- remove generic definition of show and summary, import show from
  methods and summary from stats4 <2012-02-23, Thu>

- redefine functions as S3 methods for mapping ID among gene and Term.
  this will make enrich.internal which calling these mapping function
  more robust <2012-02-23, Thu>

Changes in version 1.1.4:

- add Enrichment Analysis session in vignette. <2012-02-22, Wed>

- optimize enrichDO, ten time faster. <2012-02-22, Wed>

- separate code of enrichDO to enrich.internal, make it more general,
  and can be applied to other ontology. <2012-02-22, Wed>

- rename enrichDOResult class to enrichResult and add slot
  geneInCategory. <2012-02-22, Wed>

- export infoContentMethod and wangMethod. <2012-02-22, Wed>

Changes in version 1.1.3:

- update infoContentMethod to make it consistent between DOSE and
  GOSemSim. <2011-12-31, Sat>

Changes in version 1.1.2:

- change to using roxygen for generating Rd docs

Changes in version 1.1.1:

- add function rebuildAnnoData

- update Disease-Gene Mapping data

easyRNASeq
----------

Changes in version 1.1.10:

NEW FEATURES

- Added a naPositionFilter extending ShortRead srFilters

BUG FIXES

- Worked on Wade Davis case with 3 different sets of chromosome names
  in the three different input (reads, annotation, chromosome sizes)

- Worked on smoother error handling when not using bam files. Again
  through Wade Davis example

- Ensured that chromosome names conversion occurs whether provided with
  a factor or a character vector

- extended the NAMESPACE

Changes in version 1.1.9:

BUG FIXES

- Added an \alias{RNAseq} to ease the class documentation access; an
  H.Pages suggestion

- Changed the DESCRIPTION file to make sure that the latest ShortRead
  (1.13.13) and BiocGenerics (0.1.11) package are required

Changes in version 1.1.8:

BUG FIXES

- Corrected the last occurence of the deprecated matchMatrix call

- Corrected an issue on windows raised by the parallel package. Thanks
  to Wade Davis for pointing that one out.

Changes in version 1.1.7:

BUG FIXES

- Fixed the NAMESPACE and the vignette generation

Changes in version 1.1.6:

BUG FIXES

- Thanks to Francesco Lescai, a bug was fixed. Namely I was not
  expecting the chromosome names in the bam files and in the chromosome
  name lists to be two different set with a common intersect. I always
  consider that one would be the subset of the other one. Now, when
  such situation occurs only the common set is kept and used for the
  calculations.

- Herve Pages changed the findOverlaps value. It is now an object of
  the Hits class that does not support the matchMatrix accessor
  anymore. The code was adapted to the new accessor queryHits.

- Corrected the package structure to add a vignettes sub-directory.
  Moved the relevant files there

Changes in version 1.1.5:

BUG FIXES

- Removed the fitInfo method extension to the DESeq package as it as
  been implemented in that package

- Reworked the plotDispersionEstimates and .normalizationDispatcher
  function to deal with the new fitInfo function (the information is
  stored in an environment rather than in a list)

Changes in version 1.1.4:

- Package introduced in Bioconductor

EDASeq
------

Changes in version 1.2.0:

- Fixed a bug in methods plotQuality and plotNtFrequency that now allow
  the strand information to be missing in the SAM/BAM file.

- Fixed a bug in method biasPlot. It doesn't require anymore that one
  of the column of pData is called "conditions".

- biasPlot method of SeqExpressionSet has now the argument "col" to
  specify the column of pData to use for color coding.

edgeR
-----

Changes in version 2.6.0:

- edgeR now depends on limma.

- Considerable work on the User's Guide. New case study added on
  Pathogen inoculated arabidopsis illustrating a two group comparison
  with batch effects. All the other case studies have been updated and
  streamlined. New section explaining why adjustments for GC content
  and mappability are not necessary in a differential expression
  context.

- New and more intuitive column headings for topTags() output. 'logFC'
  is now the first column. Log-concentration is now replaced by
  log-counts-per-million ('logCPM'). 'PValue' replaces 'P.Value'. These
  column headings are now inserted in the table of results by
  exactTest() and glmLRT() instead of being modified by the show method
  for the TopTags object generated by topTags(). This means that the
  column names will be correct even when users access the fitted model
  objects directing instead of using the show method.

- plotSmear() and plotMeanVar() now use logCPM instead of logConc.

- New function glmQLFTest() provides quasi-likelihood hypothesis
  testing using F-tests, as an alternative to likelihood ratio tests
  using the chisquare distribution.

- New functions normalizeChIPtoInput() and calcNormOffsetsforChIP() for
  normalization of ChIP-Seq counts relative to input control.

- New capabilities for formal shrinkage of the logFC. exactTest() now
  incorporates formal shrinkage of the logFC, controlled by argument
  'prior.count.total'. predFC() provides similar shrinkage capability
  for glms.

- estimateCommonDisp() and estimateGLMCommonDisp() now set the
  dispersion to NA when there is no replication, instead of setting the
  dispersion to zero. This means that users will need to set a
  dispersion value explicitly to use functions further down the
  analysis pipeline.

- New function estimateTrendedDisp() analogous to
  estimateGLMTrendedDisp() but for classic edgeR.

- The algorithms implemented in estimateTagwiseDisp() now uses fewer
  grid points but interpolates, similar to estimateGLMTagwiseDisp().

- The power trend fitted by dispCoxReidPowerTrend() now includes a
  positive asymptote. This greatly improves the fit on real data sets.
  This now becomes the default method for estimateGLMTrendedDisp() when
  the number of genes is less than 200.

- New user-friendly function plotBCV() displays estimated dispersions.

- New argument target.size for thinCounts().

- New utility functions getDispersion() and zscoreNBinom().

- dimnames() methods for DGEExact, DGELRT and TopTags classes.

- Function pooledVar() removed as no longer necessary.

- Minor fixes to various functions to ensure correct results in special
  cases.

ExiMiR
------

Changes in version 1.3.1:

- NEWS: add NEWS file for tracking changes.

- inst/doc/ExiMiR-vignette.Rnw: Fix vignette generation issue

exomeCopy
---------

Changes in version 1.2.0 (2012-02-26):

- support for multisample projects: distribute sample segmentation
  across workstations, compile results and plot CNVs across samples
  (see vignette for a reworked example)

- log odds ratios for CNV segments using the fitted distribution of
  counts for the predicted and normal state

- strand-aware read starts and duplicate read removal for counting
  reads from a BAM file in genomic ranges

- functions for generating background read depth and calculating GC
  content

flowCore
--------

Changes in version 1.21.5:

SIGNIFICANT USER-VISIBLE CHANGES

- add .readFCSdataRaw routine to read FCS containing bit-packed integer
  data (with odd-bitwidth like 9,11 instead of 8,16,32,64)

- Currently the bit-wise manipulation is done within R,it can be moved
  to C if speed issue becomes a problem in the future.

Changes in version 1.21.1:

SIGNIFICANT USER-VISIBLE CHANGES

- add TEXT segment parser in readFCStext function for FCS3 when the
  delimiter characters existing inside of keyword values. Note this
  parser require all keywords and their values to be non-empty, which
  conforms to FCS3 standard

GeneGroupAnalysis
-----------------

Changes in version 1.1:

- Considering on having more flexibility on the linear models that can
  be fitted by the Gibb's sampler.

GeneticsPed
-----------

Changes in version 2007-09-12:

- Added prune function for pruning/trimming the pedigree. It does not
  work on pedigree, but assumes a data.frame with defined structure.
  Adaption needed.

- We depend on genetics package since gpLong2Wide and hwp assume that
  input is of genotype class.

- Added two utility functions (gpLong2Wide and hwp) for work with
  gpi(). There is also a separate help page for them.

- Internal fixes in gpi - no need to transpose inputs for Fortran call
  anymore - check that there are no NA values in gp and hwp

Changes in version 2007-04-25:

- Version bump to follow BioC releases. # 0.1.4

Changes in version 2007-04-19:

- Added internal .get* functions for retrieving slot values and use of
  it.

Changes in version 2007-04-18:

- R CMD check should now fail also when R error (usually call to
  stop()) occurs in unit testing.

- Added unit tests for sort, examples in help page and clarified sort
  help page.

Changes in version 2007-04-06:

- New small dataset Falconer5.1.

- Added sex(), sex<-(), ascendantSex() and ascendantSex<-() functions.

- Added some more tests for relationshipAdditive, inverseAdditive and
  inbreeding. # 0.1.2

Changes in version 2007-04-01:

- Now we are more rigorous for value of ascendantSex argument in
  Pedigree(). It must accord to values in sex column, if that one is
  passed of course.

- MASS added to depends due to use of fractions() in many places in
  documnetation.

- Added vignette on quantitative genetic (animal) model and
  model.matrix.Pedigree() functions for educational purposes.

- Created data directory and added pedigree example Mrode2.1 and
  Mrode3.1.

- Added unit tests for genetic relationship matrices that should test
  all related functions - mainly against Mrode's book examples -
  runit.genRelMatrix.R.

- Added arguments sort and names to inverseAdditive(),
  relationshipAdditive() and inbreeding().

- Reworked core for relationshipAdditive() into vectorized form.

- Added geneFlowT(), geneFlowTinv(), geneFlowM() and
  mendelianSamplingD() functions.

Changes in version 2007-04:

- Implemented sort by pedigree information only. # 0.1.3

Changes in version 2007-03-02:

- Registration of native routines - src/register.cc.

- Added gpi() function. # 0.1.1

- Temporarily removed unknown funcs, due to planned move to BioC and
  change to S4.

- All depends are now in gdata --> removing ggmisc.

Changes in version 2006-03-29:

- codeUnit to get internal codes for subject and ascendants if factors
  are used.

- Handled unused levels in nlevels.pedigree and summary.pedigree now
  produces a simple summary

- Started to implement checks in check.pedigree and friends

- Proper factor handling

Changes in version 2006-03-16:

- Playing around with NA/unknown representation - it is very likely
  that I messed up some things with factors --> subject to changes. #
  Version 0.1

- Initial version # NEWS ends here

GenomicFeatures
---------------

Changes in version 1.8:

NEW FEATURES

- Added asBED and asGFF methods to convert a TranscriptDb to a GRanges
  that describes transcript structures according to either the BED or
  GFF format. This enables passing a TranscriptDb object to
  rtracklayer::export directly, when targeting GFF/BED.

GenomicRanges
-------------

Changes in version 1.8.0:

NEW FEATURES

- Add GappedAlignmentPairs class (with accessors first(), last(),
  left(), right(), seqnames(), strand(), isProperPair()), and
  readGappedAlignmentPairs() for dealing with paired-end reads. Most of
  the GappedAlignments functionalities (e.g. coercion to GRangesList,
  "findOverlaps" and related methods, "coverage", etc...) work on a
  GappedAlignmentPairs object.

- Add encodeOverlaps,GRangesList,GRangesList,missing and related
  utilities flipQuery(), selectEncodingWithCompatibleStrand(),
  isCompatibleWithSplicing(), isCompatibleWithSkippedExons() and
  extractSkippedExonRanks().

- Add 'order.as.in.query' arg to grglist() and rglist().

- SummarizedExperiment gains direct access to colData columns with $,
  $<-, [[, and [[<- methods

- Add map,GenomicRanges,GRangesList and
  map,GenomicRanges,GappedAlignments methods. These allow mapping from
  genome space to transcript space, and genome space to read space,
  respectively.

- Add seqinfo methods (and friends) for RangedData, RangesList, and
  other IRanges data structures. These use metadata(x)$seqinfo.

- Add disjointBins,GenomicRanges.

- Add score,GRangesList and score,GenomicRanges (gets the score column
  like for RangedData).

- Add RangedDataList -> GenomicRangesList coercion.

- Add RleViewsList -> GRanges coercion.

- Add pintersect,GRangesList,GRangesList

- Add stack,GenomicRangesList

- ignore.strand argument now more uniformly supported on set
  operations.

- Add Ops,GenomicRanges (from rtracklayer).

- Add strand,Rle (only logical-Rle is supported).

- Add compare,GenomicRanges

- Add 'drop.empty.ranges' arg (FALSE by default) to low-level cigar
  utilities cigarToIRanges(), cigarToIRangesListByAlignment(), and
  cigarToIRangesListByRName().

- Add 'reduce.ranges' arg to cigarToIRangesListByAlignment().

SIGNIFICANT USER-LEVEL CHANGES

- grglist,GappedAlignments now carries over element metadata.

- Names are no longer forced to be unique when unlisting a GRangesList
  with use.names=TRUE.

- seqnames() is now preferred over rname() on a GappedAlignments
  object.

- cigarToIRangesListByAlignment() now returns a CompressedIRangesList
  instead of CompressedNormalIRangesList.

- Low-level CIGAR utilities now ignore CIGAR operation P (instead of
  trowing an error).

- The 'weight' arg in "coverage" method for GenomicRanges objects now
  can also be a single string naming a column in elementMetadata(x).

- Ranges outside the sequences bounds of the underlying sequences are
  now accepted (with a warning) in
  GenomicRanges/GRangesList/GappedAlignments objects.

- When called with 'ignore.strand=TRUE', the "range" and "disjoin"
  methods for GenomicRanges objects now behave like if they set the
  strand of the input to "*" before they do any computation.

- When called with 'ignore.strand=TRUE', "reduce" method for
  GenomicRanges objects, and "union", "intersect" and "setdiff" methods
  for GRanges objects now set the strand of their arguments to "*"
  prior to any computation.

- No more mangling of the names when combining GRanges objects ("c"
  method for GRanges objects was trying to return unique names).

- Remove isCircularWithKnownLength() generic and methods (nobody knows,
  uses, or needs this).

BUG FIXES

- flank,GRangesList no longer forces 'use.names' to TRUE and 'both' to
  FALSE.

- range,GenomicRanges was broken when object had no ranges

- Fix integer overflow issue that can occur in cigarQNarrow() or
  cigarQNarrow() when the cigar vector is very long.

genoset
-------

Changes in version 1.4.19:

SIGNIFICANT USER-VISIBLE CHANGES

- *Minor API change* Rarely used (by me) method, orderedChrs, gone.
  It's just chrOrder(chrNames(x)) anyway. 'names' on GenoSet
  depricated.  chrNames gives universal way to get chromosome names
  (i.e. names, seqlevels) for GenoSet, RangedData, GRanges.

Changes in version 1.4.10:

SIGNIFICANT USER-VISIBLE CHANGES

- segTable for Rle now optionally takes chrIndices table, start and
  stop from locData for speed. segTable for DataFrame uses this trick.
  Much faster.  About 95% time reduction for a large dataset on a large
  chip.

Changes in version 1.4.9:

- *** API Changes *** segTable on a DataFrame of Rle now has a "stack"
  argument to rbind the resulting list of data.frames of per-sample
  segments into on giant data.frame.  A "Sample" column will be added
  to separate samples.  The list of individual data.frames no longer
  has an "ID" column.  Also, some refactoring to speed up this method.

ggbio
-----

Changes in version 1.1.8:

NEW FEATURES

- create lower level API and rewrite higher level API

- new geom: geom_alignment, geom_chevron, geom_arch, geom_arrow,
  geom_arrowrect

- redefined geom: geom_rect, geom_segment

- new stat: stat_aggregate, stat_coverage, stat_mismatch, stat_gene,
  stat_table, stat_stepping

- redefined stat: stat_identity

- new layout: layout_circle, layout_karyogram

- tracks function are more smart and with more accessors.

- themes provided.

- More supported object for autoplot: VCF, ExpressionSet,
  GenomicRangesList.

SIGNIFICANT USER-LEVEL CHANGES

- qplot changed to generic autoplot function.

- argument use only "facets", no alias "facet_gr" or "facet" accepted.

- plotMismatchSum will be replaced by stat_mismatch, or autoplot,
  BamFile.

Notes

- new website for ggbio: http://tengfei.github.com/ggbio hosting docs,
  tutorials and case study

- pdf version vignette will not longer supported or just provide a
  short form.

GOSemSim
--------

Changes in version 1.13.6:

- remove dependency of organism annotation packages. <2012-03-09, Fri>
  User not need to install all these annotation packages for using
  GOSemSim. User only need to install the specific organism annotation
  package they want to calculate.

- update IC data sets for 1.14 release. <2012-03-30, Fri>

Changes in version 1.13.5:

- fixed BibTeX database file .bib. month = , must be month = someMonth,
  leave it blank will cause texi2dvi failed. <2012-03-01, Thu>

Changes in version 1.13.4:

- update vignette. <2012-02-28, Tue>

Changes in version 1.13.3:

- bug fixed for multiple annotation. <2012-02-01, Fri>

Changes in version 1.13.2:

- update infoContentMethod to make it consistent between DOSE and
  GOSemSim. <2011-12-31, Sat>

Changes in version 1.13.1:

- remove dependency of DOSE

- remove Streptomyces coelicolor support, as the genome wide annotation
  package contributor no longer supports it.

gprege
------

Changes in version 0.99.1 (2011-07-26):

- Fixed error that exhaustivePlot gave on profiles with zero data
  variance. Now it simply ignores such profiles.

Changes in version 0.99.0 (2011-07-11):

- Initial version in Bioconductor.

GRENITS
-------

Changes in version 1.7.1 (2012-03-15):

- Corrected ggplot2 (0.9.0) problems

GSVA
----

Changes in version 1.4:

USER VISIBLE CHANGES

- removed the system-requirement dependency from the GNU Scientific
  Library

- added two additional gene-set expression summarization methods:
  single-sample GSEA from Barbie et al. (Nature, 2009) and a combined
  Z-score method similar to the one used by Lee et al. (PLos Comp Biol,
  2008) via a new 'method' argument in the 'gsva()' function

- added handling of RNA-seq expression data matrices by the GSVA method
  with a new 'rnaseq' argument in the 'gsva()' function

- added a method with signature(expr="matrix",
  gset.idx.list="GeneSetCollection", annotation="character") which did
  not exist before. Now gsva() accepts the following pairs of data
  structures storing expression data and gene sets:
  ExpressionSet-GeneSetCollection, ExpressionSet-list,
  matrix-GeneSetCollection and matrix-list

BUG FIXES

- matching of gene IDs from ExpressionSet objects to GeneSetCollection
  objects now also works with Entrez-based gene IDs in ExpressionSet
  objects (e.g., when annotation(eset) == "org.Hs.eg.db") by using
  GSEABase >= 1.17.4

GWASTools
---------

Changes in version 1.1.9:

- anomSegStats checks for SNPs in centromere gaps.

- anomStatsPlot has option to plot LRR/BAF individually (for greater
  flexibility in layout).

- Updates to arguments for plot titles in chromIntensityPlot,
  anomStatsPlot, and pseudoautoIntensityPlot for consistency.

- plinkCheck has map.alt argument to override default GenotypeData ->
  PLINK annotation conversion.

Changes in version 1.1.8:

- Updated positions of pseudoautosomal regions.

- Added plinkToNcdf to convert PLINK files to NetCDF for use in
  GWASTools.

Changes in version 1.1.7:

- chromIntensityPlot and pseudoautoIntensityPlot have cex=0.5 by
  default.

- chromIntensityPlot colors now match anomStatsPlot colors.

- plinkCheck has options to skip checking parents and sex.

- plinkCheck sorts alleles by character to avoid phase mismatches.

- plinkWrite and plinkCheck print progress messages if verbose=TRUE.

Changes in version 1.1.6:

- duplicateDiscordance and duplicateDiscordanceAcrossDatasets use only
  one pair of scans per subject by default.

- duplicateDiscordanceProbability sets small negative values to 0.

Changes in version 1.1.5:

- duplicateDiscordance has an option to compute correlation by SNP.

- Added scan.exclude argument to plinkCheck.

Changes in version 1.1.4:

- Added ncdfSetMissingGenotypes function.

- plinkCheck now writes a log file with all mismatches found.

- duplicateDiscordance excludes Y chrom SNPs for females.

- duplicateDiscordance has an option to consider only pairs involving
  the minor allele.

Changes in version 1.1.3:

- batchChisqTest and batchFisherTest now return n results for n batches
  even if n=2.

- batchFisherTest has return.by.snp=FALSE as default.

Changes in version 1.1.2:

- Added LR tests to assocTestRegression.

- Bug fix in calculation of mean odds ratio in batchFisherTest.

- Bug fix in missingGenotypeByScanChrom for data sets with only one
  female.

Changes in version 1.1.1:

- Added functions plinkWrite and plinkCheck for writing and checking
  PLINK ped and map files.

- Added pcaSnpFilters data set for identifying regions with high PC-SNP
  correlation.

HTqPCR
------

Changes in version 1.9:

SIGNIFICANT USER-VISIBLE CHANGES

- Altered plotCtCor to plot 1-correlation instead of correlation.

- Altered qPCRset object to inherit from eSet. This extends the range
  of (meta) data that can be included.

NEW FEATURES

- qPCRset now contains slots for phenoData, featureData, protocol
  experiment etc. inherited from eSet.

- readCtData has been expanded to include file formats from multiple
  qPCR detection systems and vendors.

BUG FIXES

- 

htSeqTools
----------

Changes in version 1.1.5:

- Added missing documentation

Changes in version 1.1.4:

- Fixed bug in cmds for method with empty chromosomes

Changes in version 1.1.3:

- Fixed bug in alignPeaks (no longer directly accessing matchMatrix
  slot of RangesMatchingList object)

- Added plotMeanCoverage

- Fixed bug in gini for method with no chromosome lengths

Changes in version 1.1.2:

- Add methods for objects of class GRanges and GRangesList.

- Add set.seed to vignette.

- Use pvec instead of mclapply in tabDuplReads and filterDuplReads.

- Correct rpkm form enrichedRegions and remove it from enrichedPeaks

- Fixed giniCoverage for ranges with low numbers of reads

Changes in version 1.1.1:

- Adjusted parallel computing in enrichedPeaks so that it no longer
  spans an uncontrolled number of child processes when mc.cores>1

- Added arguments "labels" and "cex.text" to the plot method for
  cmdsFit objects.

- Added monotonicity contraint to filterDuplReads and fdrEnrichedCounts
  to ensure that the estimated FDR decreases with the number of repeats

- Fixed overflow problem in enrichedChrRegions which occurred for long
  genomes (e.g. human)

- Adjusted behavior of stdPeakLocation so that it is consistent with
  PeakLocation

- Fixed bug in RPKM calculation by enrichedRegions

- Added option to compute Spearman correlations in cmds

IRanges
-------

Changes in version 1.14.0:

NEW FEATURES

- The map generic and RangesMapping class for mapping ranges between
  sequences according to some alignment. Some useful methods are
  implemented in GenomicRanges.

- The Hits class has experimental support for basic set operations,
  including setdiff, union and intersect.

- Added a number of data manipulation functions and methods, including
  mstack, multisplit, rename, unsplit for Vector.

- Added compare() generic for generalized range-wise comparison of 2
  range-based objects.

- Added OverlapEncodings class and encodeOverlaps() generic for dealing
  with "overlap encodings".

- subsetByOverlaps() should now work again on an RleViews object.

- DataFrame now supports storing an array (like a matrix) in a column.

- Added as.matrix,DataFrame method.

- Added merge,DataTable,DataTable method.

- Added disjointBins,RangesList method.

- Added ranges,Rle and ranges,RleList methods.

- Added which.max,Rle method.

- Added drop,AtomicList method.

- Added tofactor() wrapper around togroup().

- Added coercions from vector to any AtomicList subtype (compressed and
  uncompressed).

- Added AtomicList to Character/Numeric/Logical/Integer/Raw/ComplexList
  coercions.

- Added revElements() for reversing individual elements of a List
  object.

SIGNIFICANT USER-VISIBLE CHANGES

- RangesMatching has been renamed to Hits and extends Vector, so that
  it supports element metadata and other features.

- RangesMatchingList has been renamed to HitsList.

- The 2 columns of the matrix returned by the "as.matrix" method for
  Hits objects are now named queryHits/subjectHits instead of
  query/subject, for consistency with the queryHits() and subjectHits()
  getters.

- queryLength()/subjectLength() are recommended alternatives to
  dim,Hits.

- breakInChunks() returns a PartitioningByWidth object.

- The 'weight' arg in "coverage" methods for IRanges, Views and
  MaskCollection objects now can also be a single string naming a
  column in elementMetadata(x).

- "countOverlaps" methods now propagate the names of the query.

DEPRECATED AND DEFUNCT

- matchMatrix,Hits is deprecated.

- Moved the following deprecated features to defunct status: - use of
  as.data.frame() or as( , "data.frame") on an AtomicList object; - all
  coercion methods from AtomicList to atomic vectors; - subsetting an
  IRanges by Ranges; - subsetting a RangesList or RangedData by
  RangesList.

BUG FIXES

- within,RangedData/List now support replacing columns

- aggregate() override no longer breaks on . ~ x formulas

- "[", "c", "rep.int" and "seqselect" methods for Rle objects are now
  safer and will raise an error if the object to be returned has a
  length > .Machine$integer.max

- Avoid blowing up memory by not expanding 'logical' Rle's into logical
  vectors internally in "slice" method for RleList objects.

isobar
------

Changes in version 1.1.3:

- better matching of file patterns of peaklist and id in report

- tab2xls improvements:

- fix when there are cells with preceeding colons - would think they
  are cell properties

- fix row limitation of 65536 - add new worksheet with remaining lines

- re-added ibspiked_set2 dataset as the xz requirement allows for
  additional data

Changes in version 1.1.2:

- fixed handling of divergent identifications in one search engine

- fixed number of spectra in isobar-analysis report

- fixed recently introduced error when reading mgf file

- identifications tab in XLS report is now in concise format

- shared peptides are colored in gray

- added xls report format = wide

Changes in version 1.1.1:

NEW FEATURES

- normalization can now be performed on individual channels (and
  channel pairs)

- added semi-quantitative Quantitation with emPAI, dNSAF and spectral
  count

- proteinInfo can now be gathered from Uniprot directly

- added reporter intensity plot shpwing effect of normalization

- added linear regression as ratio estimator

- improved MA plot: added 'Infinity' on the axis

KEGGgraph
---------

Changes in version 1.11.1:

SIGNIFICANT USER-VISIBLE CHANGES

- The colorectalcancerSPIA dataset has been compressed with R CMD build
  --resave-data

limma
-----

Changes in version 3.12.0:

- read.maimages() with source="agilent" now reads median foreground
  estimates instead of mean foreground.  New option source=
  "agilent.mean" preserves earlier meaning of source="agilent".

- Agilent single-channel case study added to User's Guide.

- removeBatchEffect() now corrects for continuous covariates as well as
  qualitative factors.

- new function camera() performs competitive gene set tests while
  adjusting for inter-gene correlation.

- new function interGeneCorrelation() estimates the average intergene
  correlation for a set of genes.

- columns in output from roast() have been re-ordered.

- arguments 'selected' and 'selected2' renamed to 'index' and 'index2'
  in functions barcodeplot(), geneSetTest() and wilcoxGST().

- default labels for barcodeplot() are now somewhat more explicit.

- new function rankSumTestWithCorrelation extends the
  Wilcoxon-Mann-Whitney test to allow for correlation between cases in
  one of the groups.  geneSetTest() now calls this function instead of
  wilcox.test, with a consequence improvement in speed.

- The lfc (log-fold-change) cutoff argument of topTable() is now
  applied to the minimum absolute logFC when ranking by F-statistic.
  Previously lfc was only used when ranking by t-statistic.

- new methods "fast" and "affy" for normalizeCyclicLoess(), with "fast"
  becoming the default method. New argument 'cyclic.method' for
  normalizeBetweenArrays() gives access to the different cyclic loess
  methods.

- There were problems with using the argument gene.weights in mroast().
  This argument is now permitted to be of the same length as the number
  of probes in the data set.  It is then automatically subsetted for
  each gene set.

- mroast() now uses mid-p-values by default when adjusting for multiple
  testing.

- neqc(), nec() and normexp.fit.control() now give user-friendly error
  messages when no negative control probes or no regular probes are
  found.

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

minfi
-----

Changes in version 1.1:

- Changed NAMESPACE file

- Defined constructors for MethylSet, RGChannelSet,
  RGChannelSetExtended.

- Included a version number in the class definition for MethylSet and
  RGChannelSet.  Old objects can be updated by calls of the form
  updateObject(Mset).

- read.manifest (not exported) updated to include nCpGs.

- preprocessSwan was added.  Still work in progress.

- Changed background calculation in preprocessSwan.

- Added a section to the vignette describing preprocessSwan.

- Bug fix: ilogit2 is now in base (it used to be base e).  Thanks to
  Time Triche, Jr <tim.triche@gmail.com>.

- Added and dcoumented the IlluminaMethylationAnnotation class; still
  work in progess.

- Moved package vignette from inst/doc to vignettes.

MinimumDistance
---------------

Changes in version 1.0:

USER VISIBLE CHANGES

- Using NEWS.Rd

MLP
---

Changes in version 1.3.0:

- update using new annotation packages

- move to roxygen2

Changes in version 1.2.1:

- also export addGeneSetDescription

mosaics
-------

Changes in version 1.3.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Improve help documents for all classes and functions.

BUG FIXES

- mosaicsPeak(): Correct bin size calculation when binsize=NA.

Changes in version 1.3.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Simplify arguments of mosaicsRunAll(), constructBind(), and export().

- Add parallel argument in mosaicsFit().

- Extensive use of parallel processing/computing.

- Overall speed improvements in the package.

- Update the vignette.

- Use parallel package instead of multicore package.

Changes in version 1.2.5:

- Correct version number in DESCRIPTION and package?mosaics.

Changes in version 1.2.4:

SIGNIFICANT USER-VISIBLE CHANGES

- Add parallel argument in readBins().

- Add parallel argument in mosaicsRunAll().

BUG FIXES

- DESCRIPTION: 'multicore' package in 'Enhances' instead of 'Suggests'.

Changes in version 1.2.3:

NEW FEATURES

- New model for deeply sequenced ChIP-seq data.

- Genome-wide analysis of ChIP-seq data is now available.

- Supports more aligned read file formats: eland_result,
  eland_extended, eland_export, bowtie, SAM, BED, CSEM.

- Preprocessing of aligned read files can be done within the R
  environment using constructBins().

- Easier model fitting for the two sample analysis using
  mosaicsRunAll().

- Preprocessing and model fitting become much faster (Rcpp).

- Parallel processing/computing is now supported (multicore).

SIGNIFICANT USER-VISIBLE CHANGES

- Add constructBins(): Preprocess aligned read files to bin-level
  files.

- Add mosaicsRunAll(): Convenient two sample analysis.

- Add bgEst argument in mosaicsFit(): Choose background estimation
  approach.

- Add nCore argument in readBins(): Parallel processing.

- Vignettes is now extensively updated.

- Rcpp package is required and multicore package is suggested.

DEPRECATED AND DEFUNCT

- Drop chrID argument in export().

BUG FIXES

- Fix mosaicsPeak() for the case that no peak is called.

- Fix export() by removing unnecessary spaces in output text files.

MSnbase
-------

Changes in version 1.3.15:

- new updateFeatureNames function <2012-02-17 Fri>

- updated vignettes to illustrate vertical/horizontal combine
  <2012-02-17 Fri>

- typo in normalised.Rd <2012-02-19 Sun>

- TODO combine unit tests

Changes in version 1.3.14:

- new is.na.MSnSet <2012-02-16 Thu>

- updated vignette and NA related man pages with cross-links
  <2012-02-16 Thu>

Changes in version 1.3.13:

- new plotNA method + doc <2012-02-15 Wed>

- new filterNA method + doc + tests <2012-02-15 Wed>

- added a check on 'n' in topN <2012-02-16 Thu>

- created a .Rinstignore <2012-02-16 Thu>

- Update package Rd <2012-02-16 Thu>

Changes in version 1.3.12:

- new topN methods + doc + tests <2012-02-14 Tue>

Changes in version 1.3.11:

- changed explicit close(file) in writeMgf methods to
  on.exit(close(file)) <2012-02-12 Sun>

- typo in vignette <2012-02-13 Mon>

Changes in version 1.3.10:

- type in writeMgfData man <2012-02-07 Tue>

- updated TITLE in writeMgfData <2012-02-08 Wed>

Changes in version 1.3.9:

- updated demo vignette <2012-02-03 Fri>

- sorting numeric subsets in "[" pSet, as unsorted numerical indexes
  fails <2012-02-03 Fri>

- Added match.arg in combineFeatures so that a unique default value
  (the first) is used when no fun is specified <2012-02-03 Fri>

Changes in version 1.3.8:

- Modified trimMz warning to report acquisition number <2012-02-01 Wed>

- add 'experimentData(object, value) <- ' method for signature eSet and
  MIAPE <2012-02-02 Thu>

- combine methods for MIAPE instances <2012-02-02 Thu>

- combine methods for MSnProcess instances <2012-02-02 Thu>

- changed qual drop warning into message in combineFeatures, updated
  test_MSnSet accordingly <2012-02-02 Thu>

- new updateFvarLabels and updateSampleNames function <2012-02-03 Fri>

- combine method for MSnSets <2012-02-03 Fri>

- Updated demo vignette figure 8 <2012-02-03 Fri>

Changes in version 1.3.7:

- Speeded up writeMgfData <2012-01-28 Sat>

- fixes for ggplot2 0.9.0

- added import(grid) and import(reshape) <2012-01-30 Mon>

- importFrom(plyr, ...) instead of only llply <2012-01-31 Tue>

- loading reshape and grid in vignette <2012-01-31 Tue>

- fixed chunk 21 (label = quantitation-plot) <2012-01-31 Tue>

Changes in version 1.3.6:

- Updated NoteAboutSpeedAndMemory since parallel processing has been
  added. <2011-12-18 Sun>

- Added CITATION <2012-01-27 Fri>

- Added information to header output: acquisition number and precursor
  intensity <2012-01-27 Fri>

- Added a test in plot.Spectrum2 for empty dataframe <2012-01-27 Fri>

- moved foreach, doMC to enhances <2012-01-27 Fri>

Changes in version 1.3.5:

- added a gc() before mzR::close(msdata)... seems to help with Rcpp and
  ref classes issue. <2011-12-09 Fri>

- added a show parameter to getCacheEnv to define .cache should be
  printed out before being returned. <2011-12-09 Fri>

- added cache unit test <2011-12-09 Fri>

- readMzXMLData is now defunct and remove xcms from Imports <2011-12-16
  Fri>

Changes in version 1.3.4:

- fixed bug in show MSnExp method for MS1 experiments.  When loading
  MS1 spectra, cache is set to 0. Bug reported by Jesse Meyer.
  <2011-12-06 Tue>

- fixed another bug/typo in readMSData <2011-12-06 Tue>

- now running extractSprectum example again <2011-12-06 Tue>

- setting default cache to 0, as cache=1 introduces unstabel
  behavious... will investigate that <2011-12-06 Tue>

Changes in version 1.3.3:

- added parallel computation for MSnExp quantitation using foreach with
  llply(..., .parallel=TRUE) <2011-12-03 Sat>

- TODO document above in quantify-methods.Rd

- added foreach and doMC in Suggests <2011-12-03 Sat>

- added Spectrum removePeaks and clean'ing in readMSData <2011-12-05
  Mon>

Changes in version 1.3.2:

- \dontrun{} extractSpectrum example, as this seems to be a major
  offender producing the intermittent check 'Error in function (x) :
  attempt to apply non-function' error <2011-11-07 Mon>

- typo in Author@R <2011-11-14 Mon>

- modified utils.removePeaks and utils.clean to call sapply instead of
  IRanges:sapply <2011-12-01 Thu>

Changes in version 1.3.1:

- Herve added BioGenerics as a dependency and import statement in
  NAMESPACE <2011-11-29 Tue>

Changes in version 1.3.0:

- Version bump for Bioc 2.10 devel

ncdfFlow
--------

Changes in version 1.1.2:

SIGNIFICANT USER-VISIBLE CHANGES

- Using temporary directory instead of working directory to store cdf
  file in creating ncdfFlowSet from flowSet

- allow for user specified path in ncdfFlowSet_sync method to save the
  cdf in different location other than original one

- clone.ncdfFlowSet function: - change argument name to avoid
  confusion:sNewNcFile-->isNew ;newNcFile-->fileName -avoid copying the
  entire cdf repository when clone subsetted ncdfFlowSet -fix the bug
  of inconsistent dimensions (sample*colnames) when create the new cdf
  file

- check whether source file exist in read.ncdfFlowSet

- .writeSlice: -allow for either flowFrame or matrix to be added by
  -add sample name to the error message to help troubleshoot the
  problematic FCS file especially for loading large datasets

- add isNew=FALSE to split method to allow for splitting into multipe
  cdf files for the sake of parallel computing

- set compress=FALSE to disable compression mode of CDF

Changes in version 1.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- netCDF support for large data sets.

- centralized storage of flow data in 3-D matrix (sample*channel*event)

- fast data accessing,subsetting and splitting

- support all the related methods for flowSet

KNOWN ISSUES

- write meta data - ncdfFlow allows user to save the entire ncdfFlowSet
  object in ncdf file.

- Currently the meta data is first serialized in R and stored as raw
  vector in cdf. It can fail when the meta data size exceeds the limit
  of serialization function.

netresponse
-----------

Changes in version 1.7.38 (2012-03-28):

SIGNIFICANT USER-VISIBLE CHANGES

- added the option to ignore network in detect.responses (netw=NULL);
  then the methods assumes fully connected network. Without speedups
  the performance may be very slow.

- added (optional) initial mutual information based filtering of the
  network edges also in the first stage where pairwise similarities are
  calculated for all node pairs This can give considerable speedups
  with large networks

- merged netresponse.visualization package in the netresponse main
  package

- added mode = hard in sample2response function

- plot.associations -> plotAssociations

- plot.pca -> plotPCA

- get.gofz -> getqofz

BUG FIXES

- changed dependency RBGL into RColorBrewer

- changed dependencies multicore and doMC into parallel

- Modified mk.hp.posterior so as to accommodate 'max number of
  responses' (c.max) option. Validations in tests/vdpmixture.R ok

- now allowing max.subnet.size = 1 in detect.responses

- order.responses: fixed minor bug which occurs when no enrichments are
  detected

Changes in version 1.7.34 (2012-03-28):

- changed dependency RBGL into RColorBrewer

Changes in version 1.7.33 (2012-03-13):

- Modified mk.hp.posterior so as to accommodate 'max number of
  responses' (c.max) option. Validations in tests/vdpmixture.R ok

Changes in version 1.7.32 (2012-02-21):

- added the option to ignore network in detect.responses (netw=NULL);
  then the methods assumes fully connected network. Without speedups
  the performance may be very slow.

Changes in version 1.7.03 (2012-02-02):

- merged netresponse.visualization package in the netresponse main
  package

Changes in version 1.7.02 (2012-02-01):

- switched to R-2.14.1

- now allowing max.subnet.size = 1 in detect.responses

- order.responses: fixed minor bug which occurs when no enrichments are
  detected

- added mode = hard in sample2response function

- added (optional) initial mutual information based filtering of the
  network edges also in the first stage where pairwise similarities are
  calculated for all node pairs This can give considerable speedups
  with large networks

- plot.associations -> plotAssociations

- plot.pca -> plotPCA

- get.gofz -> getqofz

oligo
-----

Changes in version 1.20:

USER VISIBLE CHANGES

- New getProbeInfo() function added, to simplify probe selection
  without using SQL.

- New fitProbeLevelModel() function added. It allows probe-level models
  ('plm' and 'medianpolish'), which can be used for QC.

- fitPLM, coefs and resids are now Deprecated. Use fitProbeLevelModel,
  coef and residuals respectively. 'coef' and 'residuals' follow the
  standards used elsewhere in R.

- Now using foreach for parallelization.

BUG FIXES

- Addressed issue in which probes without chr info would be removed
  from pmChr, leading to results whose dimensions did not match PM
  matrix (on TilingFeatureSet).

oneChannelGUI
-------------

Changes in version 1.21.11:

SIGNIFICANT USER-VISIBLE CHANGES

- adding interface to EDASeq package

pdInfoBuilder
-------------

Changes in version 1.20:

USER VISIBLE CHANGES

- The annotation packages for Affymetrix Tiling arrays now contain
  strand information (suggestion and implementation: Kristof De Beuf)

phyloseq
--------

Changes in version 0.99:

- Bioconductor development release updates.

- BIOM format import: import_biom() function

- Parallel Fast UniFrac

- distance() wrapper for ecological distance calculations

- ordinate() wrapper, calculates many different ordination methods.

- plot_ordination() powerful, flexible ordination plotting using
  ggplot2

- make_sample_network(), plot_sample_network() - microbiome network
  visualization

- plot_richness_estimates() for easy, flexible summary of species
  richness

- Support for Double Principle Coordinate Analysis (DPCoA)

- Several published exampled datasets included

- General importer for all supported data formats: import()

- Lots of documentation updates.

- Lots and lots of fixes and improvements.

pint
----

Changes in version 1.7.03 (2012-03-22):

- fixed missing documentations

Changes in version 1.7.01 (2012-01-16):

- added remove.duplicates option to pint.match

procoil
-------

Changes in version 1.5.1:

- made arguments xlab and ylab accessible to users

qpgraph
-------

Changes in version 1.12:

NEW FEATURES

- Results of conditional independence tests with qpCItest() are now
  returned using the R standard htest class

- qpCItest() allows one to test mixed interactions involving phenotypic
  data variables and expression profiles in ExpressionSet objects, and
  phenotypic/genetic variables and expression profiles in smlSet
  objects

- qpEdgeNrr() estimates the non-rejection rate involving phenotypic
  data variables and expression profiles in ExpressionSet objects, and
  phenotypic/genetic variables and expression profiles in smlSet
  objects

QUALIFIER
---------

Changes in version 1.0.0:

FEATURES

- outlier detection for different population
  stats(counts,proportion,MFI,Spike) of gated/ungated flow data

- xyplot and boxplot on different dimensions

- scatter plot for gated population

- HTML report with svg support

- fuzzy match of population name to select the same population from the
  different gating path

KNOWN ISSUES

- extracting populations stats needs further optimization on speed

- extend the formula parser to allow more functions to be applied by
  formula

R453Plus1Toolbox
----------------

Changes in version 1.5.2:

- Added new graphical components, options and annotations to the
  plotVariants function.

RCytoscape
----------

Changes in version 1.6.0:

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

Changes in version 0.99.0:

- change package name to ReactomePA, for there is already an RPA
  package. <2012-03-02, Fri>

- Vignette issues: <2012-03-02, Fri> change image format from .eps to
  .pdf, make it easier to build.  remove the tolatex tag of
  sessionInfo(), make the output more readble.

- re-implement geneID2Name using select method. <2012-03-02, Fri>

- add examples in man pages. <2012-03-02, Fri>

- remove man pages of internal functions. <2012-03-02, Fri>

- import plot summary from stats4, for BiocGenerics (version 0.1.10)
  removed them <2012-03-03, Sat>

Changes in version 0.2.3:

- @exportMethod plot <2012-02-15 Wed>

- fix bug when calling summary method from plot, for summary defined in
  base is S3 method, instead import summary generic from stats4
  <2012-02-15 Wed>

Changes in version 0.2.2:

- remove generic definition of show, summary and plot, add NAMESPACE
  import show from methods and plot from graphics. summary need not to
  import, for is defined in the base package <2012-02-10 Fri>

Changes in version 0.2.1:

- update vignette <2012-02-09 Thu>

- add sample data (an example list of genes from ProfCom:
  http://webclu.bio.wzw.tum.de/profcom/gene_Lists/example1.txt the gene
  symbols were converted to entrezgene) <2012-02-09 Thu>

Changes in version 0.2.0:

- separate codes of mapping pathway ID to pathway Name to a new
  function pathID2Name <2012-02-09 Thu>

- implement geneID2Name function for mapping gene ID to gene Symbol
  <2012-02-09 Thu>

- add parameter readable in function enrichPathway <2012-02-09 Thu>

- implement cnetplot for plotting category net <2012-02-09 Thu>

- modify plot function for class enrichPathwayResult to use cnetplot
  <2012-02-09 Thu>

Changes in version 0.1.1:

- add vignette <2012-02-08 Wed>

- bug fixed for multiple mapping of pathway ID to pathway Name such as
  pathway 162906 can mapping to 1629061 and 1629062 when getting
  pathway name, remain the first one. <2012-02-08 Wed>

Changes in version 0.1.0:

- implement show, summary and plot method for enrichPathwayResult class
  <2012-02-08 Wed>

- define class enrichPathwayResult to store result of enrichPathway
  <2012-02-08 Wed>

- implement enrichPathway function for enrichment analysis. using
  hypergeometric model <2012-02-08 Wed>

- initial package skeleton <2012-02-08 Wed>

RedeR
-----

Changes in version 1.1.16:

NEW FEATURES

- Loading performance, xml serialization, data packing and post
  formats.

- Server/client connection.

- Control over the app main features from R.

- Node-container assignment options.

BUG FIXES

- Correction/fine-tune of dynamic layout function under
  zoom-in/zoom-out requests.

- Merge out-edge function is fixed to rescale with node size.

- Legend function is fixed for node/edge shape attributes.

- Legend color is fixed for two-palette option.

- The method 'addGraph' is fixed to load/convert directed graphs from
  igraph to RedeR.

ReQON
-----

Changes in version 1.2.0:

NEW FEATURES

- Output BAM file keeps input header and adds header line: "@CO Quality
  scores were recalibrated with ReQON."

- Allow threshold options (nerr and nrf) to remove positions from the
  training set that are likely to contain incorrect error calls.
  (e.g., novel variants and systematic mapping errors)

- diagnostic output now outputs flagged read positions ($FlagPos) and
  regression coefficients ($coeff)

RPA
---

Changes in version 1.11.13 (2012-02-25):

- modifications to accommodate single-probe probesets without errors

- sigma2.method default to "robust" in functions RPA.sigma2.update and
  rpa.fit

- changed defaults in set.alpha function

- added missing data imputation in rpa.fit

Changes in version 1.11.12 (2011-12-13):

- removed cind option from update.hyperparameters in RPA.online

Changes in version 1.11.05 (2011-11-11):

- online functions work properly

rqubic
------

Changes in version 1.8.1 (2011-08-02):

- Fix errors in documentations

Changes in version 1.8.0 (2011-07-11):

- Submission to Bioconductor

Changes in version 1.1-3 (2012-01-05):

- Minor bug fix in readBiclusterResults: bicluster files without any
  bicluster are recognized

- Minor bug fix in readBiclusterResults: bicluster files with only one
  bicluster and one row/column not report error anymore (previously due
  to matrix dropping)

- Minor bug fix in readBiclusterResults: feature/condition names are
  written in "features"/"conditions" items in the info list

- Feature improvement in readBiclusterResults: featureNames and
  sampleNames are written in Parameters, so as to be used by the coerce
  method in the eisa package, to coerce a Biclust object into an
  ISAModules object

- Methods features and conditions now have a two-step strategy: first
  try the info list, if failed, then try the matrix names

- length.Biclust returns 0 if the Biclust object contains no bicluster

- combineBiclusts support non-empty Biclust objects with empty Biclust
  objects, e.g. those without valid biclusters detected.

Changes in version 1.1-2 (2012-01-03):

- Add fcFilter for feature-condition filtering

Changes in version 1.1-1 (2011-12-22):

- Generalize S4 methods for QUBICBiclusterSet to Biclust: most of the
  S4 methods now can be applied to a Biclust object.

- Add the combineBiclusts method to combine multiple
  Biclust/QUBICBiclusterSet objects.

- Add readBiclusterResults function to complement the
  writeBiclusterResults function in the biclust package.

Rsamtools
---------

Changes in version 1.8.0:

NEW FEATURES

- Add readBamGappedAlignmentPairs() (plus related utilities
  findMateAlignment() and makeGappedAlignmentPairs()) to read a BAM
  file into a GappedAlignmentPairs object.

SIGNIFICANT USER-VISIBLE CHANGES

- update samtools to github commit
  dc27682f70713a70d4f31bca652cf78e00757da2

- Add 'bitnames' arg to bamFlagAsBitMatrix() utility.

- By default readBamGappedAlignments() and readBamGappedReads() don't
  drop PCR or optical duplicates anymore.

BUG FIXES

- readBamGappedAlignments handles empty 'tag' fields

- scanTabix would omit variants overlapping range ends

- scanFa would segfault on empty files or empty ids

Rsubread
--------

Changes in version 1.6.0:

NEW FEATURES

- Significant improvement on the exon-exon junction detection (subjunc
  function).

- Calling SNPs using a simple allele fraction approach (callSNPs
  function) .

- Removing duplicated reads (removeDupReads function).

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

Changes in version 1.15:

NEW FEATURES

- Integrated with tabix via Rsamtools. BED and GFF export methods now
  generate a tabix index, if index=TRUE. Most import() methods gain a
  'which' argument that takes advantage of tabix, when available.

- Added wigToBigWig() function for efficient conversion of WIG to
  BigWig.

- Added SeqinfoForBSGenome() and SeqinfoForUCSCGenome() for
  conveniently retrieving Seqinfo objects for a given genome.

- Added support for FASTA import/export via Biostrings.

- GTF and GVF files are now parsed as GFF.

SIGNIFICANT USER-VISIBLE CHANGES

- The import/export API is now based on RTLFile objects, which wrap a
  file path, URL or connection. There is an RTLFile subclass for every
  file format. This makes it easier to extend rtracklayer (export,
  import) with new file types. The existing API is still supported (and
  even encouraged for most uses).

- Handle CSV attributes in GFF3 using CharacterList columns.

- BED columns thickStart/thickEnd translate to an IRanges column named
  "thick". The blockStarts/Sizes/Count columns now map to a single
  RangesList "blocks" column.

BUG FIXES

- Numerous fixes in the import/export methods, as a result of
  implementing a full unit test suite. If something was not working for
  you in the past, please try again.

- Compression and connections should now work fairly uniformly across
  file types. (start date: 29 March, 2012)

ShortRead
---------

Changes in version 1.13:

SIGNIFICANT USER-VISIBLE CHANGES

- FastqSampler is considerably faster

- FastqSampler and FastqStreamer require explicit close() to avoid
  warnings about closing unused connections

BUG FIXES

- qa reports on very large lanes would overflow alphabetFrequency

- qa report scales adapaterContamination correctly

- FastqSampler would rarely sample fewer than requested reads

- FastqSampler supports outputs of >2^31 - 1 total nucleotides

- readFastq parses records with 0 width

TargetSearch
------------

Changes in version 1.12.0:

SIGNIFICANT USER-VISIBLE CHANGES

- New binary file format for the peak-list files, a.k.a. RI files. This
  speeds up metabolite searches by 5-10 fold. The old TEXT format is
  kept for compatibility. See method 'fileFormat'. Also see 'bin2text'
  and 'text2bin' funtions.

- New plot peak function 'plotPeak'. The old function was renamed as
  plotPeakSimple. The function show also the regions in which the
  searches were performed to provide better quality controls.

BUG FIXES

- Changed to check.names=TRUE in 'read.delim' call in functions
  TargetSearchGUI() and ImportLibrary(). This prevents empty column
  names that might produce errors in downstream functions.

VariantAnnotation
-----------------

Changes in version 1.2.0:

NEW FEATURES

- readVcf() has genome argument, can be subset on ranges or VCF
  elements with ScanVcfParam()

- scanVcfHeader() returns VCFHeader class with accessors fixed, info,
  geno, etc.

- writeVcf() writes out a VCF file from a VCF class

- locateVariants() - returns GRanges instead of DataFrame - output
  includes txID, geneID and cdsID - has cache argument for repeated
  calls over multiple vcf files

- predictCoding() - returns GRanges instead of DataFrame - output
  includes txID, geneID, cdsID, cDNA-based, cds-based and protein-based
  coordinates

xps
---

Changes in version 2.15.0:

VERSION xps-1.15.2

- move to ROOT version 5.32/01, update README

VERSION xps-1.15.1

- update zzz.R

Changes in version 2.14.0:

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

Changes in version 2.13.0:

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

Changes in version 2.12.0:

VERSION xps-1.9.9

- update method XExonChip::ProbesetLevel() for whole genome annotation
  na31 files

VERSION xps-1.9.8

- revert update of Makefile.win

VERSION xps-1.9.7

- update Makefile.win to clean xpsLinkDef.h

- update script4schemes.R for annotation na31 (updated release)

VERSION xps-1.9.6

- update information files for new ROOT Version 5.27/04 (root_v5.27/04)

- update script4xps.R, script4schemes.R for annotation na31

VERSION xps-1.9.5

- update import.data to use make.names(celnames) to protect against
  certain characters

VERSION xps-1.9.4

- update root.profile.R, macroDrawProfilePlot.C to allow selecting
  subset of trees

- in read.table() set stringsAsFactors=FALSE

VERSION xps-1.9.3

- update method XPreFilter::Calculate() to handle exon arrays correctly

VERSION xps-1.9.2

- update Makefile and Makefile.arch

VERSION xps-1.9.1

- add support in XPSAnalysis.cxx to export exon array probeset
  annotations (filters)

Changes in version 2.11.0:

VERSION xps-1.7.9

- update method validData() to handle slot data containing different
  column types

- update methods seExprTreeSet(), rleplot(), mvaplot(), nuseplot()

VERSION xps-1.7.8

- add method xpsFIRMA()

- add functions firma(), firma.expr(), firma.score()

VERSION xps-1.7.7

- update bgcorrect.R to warn from using tmpdir resulting in empty root
  file

- update normalize.R to warn from using tmpdir resulting in empty root
  file

VERSION xps-1.7.6

- minor change to allow computation with g++ 4.4.x

VERSION xps-1.7.4

- add NEWS

VERSION xps-1.7.3

- add ExprTreeSet methods validSE(), nuseplot(), rleplot()

- allow to export layout trees for incomplete *.CLF files

- update examples/updateAnnotation.R

VERSION xps-1.7.2

- add examples/updateAnnotation.R

- update script4xps.R

VERSION xps-1.7.1

- allow using mas5() and mas5.call() with plate arrays w/o MMs

- update script4xps.R

Changes in version 2.10.0:

VERSION xps-1.5.19

- update script4xps.R

VERSION xps-1.5.18

- add parameter bgcorrect.option to function mas5.call()

VERSION xps-1.5.17

- update README

VERSION xps-1.5.16

- allow handling of probesets w/o MMs on Citrus.CDF

VERSION xps-1.5.15

- validBgrd() implement 'which'

- add vignette xpsPreprocess.pdf

- add example/macro4xpsPreprocess.R

VERSION xps-1.5.14

- update export() to include read.table(..,comment.char='')

- update methods.DataTreeSet.R to allow probe-level lowess and supsmu
  normalization

VERSION xps-1.5.13

- update express() to allow setting bufsize for tree baskets

VERSION xps-1.5.12

- changes to allow import of miRNA-1_0.CDF

VERSION xps-1.5.9

- update validOption() to allow 'separate:none'

VERSION xps-1.5.8

- update rma() to allow improved ties handling as option like
  preprocessCore

VERSION xps-1.5.7

- update method validCall()

- add methods validExpr() and validPVal()

- update vignette APTvsXPS.pdf

- update examples script4xps2apt.R and script4bestmatch.R

VERSION xps-1.5.4

- update function exonLevel() to use affx=c(4,8,16,32)

- update function dataDataTreeSet() to return correct ids for mask

- add new internal function exonLevelIDs()

VERSION xps-1.5.3

- update help file exonLevel.Rd

VERSION xps-1.5.1

- update validData() to check for duplicate rownames

- allow reading of genetitan plate data

Changes in version 2.9.0:

VERSION xps-1.3.13

- update DESCRIPTION to mention root version

- update README

VERSION xps-1.3.12

- correct bug in implementation of FDR and Hochberg

VERSION xps-1.3.11

- make function exonLevel public and add exonLevel.Rd

VERSION xps-1.3.8

- update all initialize methods to prevent checkS3forClass warnings

- update bgcorrect.rma, bgcorrect.mas5

- update script4xps.R, script4exon.R

VERSION xps-1.3.6

- support using genome array as exon arrays

VERSION xps-1.3.5

- add code in function import.data() to replace dots, colons in
  celnames with underscores

VERSION xps-1.3.4

- correct bug in xpsPreprocess for add.data=FALSE

- correct sub(.root, .txt, x) to sub(\.root, .txt, x)

- update root.image() to get setname from setName()

VERSION xps-1.3.3

- change intensity<- to allow using slot data for further processing

VERSION xps-1.3.1

- add function root.merge.data()

Changes in version 2.8.0:

VERSION xps-1.1.9

- protect class XRMABackground against defect Affy chips, e.g. zero
  division

- protect root.data() etc against duplicate celnames or treenames

VERSION xps-1.1.8

- correct bug in class XINICall resulting in buffer overflow

VERSION xps-1.1.7

- add call method I/NI-call (Talloen et al)

VERSION xps-1.1.6

- prevent import of CEL-files with zero max intensity

- update functions returning ExprTreeSet to import results as option
  only

- update functions returning CallTreeSet to import results as option
  only

- update root.density etc to allow saving from R function

- add root.profile to use root graphics for boxplots

- add summarization method FARMS (Hochreiter et al)

- add summarization method DFW (Chen et al)

- update vignette xps.pdf

VERSION xps-1.1.5

- correct problem in validData for CEL-files starting with numbers

VERSION xps-1.1.4

- update vignette xps.pdf

- add new vignette APTvsXPS.pdf

- update examples script4exon.R

- add examples script4xps2apt.R and script4bestmatch.R

- update README

VERSION xps-1.1.3

- to allow CEL-names starting with a number, update to read.table(...,
  check.names=FALSE)

- update ExprTreeSet to set slot exprtype to correct type

- add functions root.expr() and root.call()

- need to change setname for dabg.call() from CallTreeSet to CallSet as
  for mas5.call()

VERSION xps-1.1.2

- allow different exonlevels for bgrd, normalization, summarization

- update replacement methods exprs, pvalData, presCall to allow
  subsetting

- add function metaProbesets to compute metacoreList.mps for apt

VERSION xps-1.1.1

- increase maximum root file size from 2GB to 2TB

- decrease computation time

- correct bug preventing export of exon probeset normalized data

Changes in version 2.7.0:

VERSION xps-0.99.11

- correct bug resulting in empty exon probeset column

VERSION xps-0.99.10

- update source code to handle tmpdir correctly on WinXP

- update examples script4xps.R and script4exon.R

VERSION xps-0.99.9

- update methods namePart, extenPart, validData, validBgrd to handle
  names with underscores

VERSION xps-0.99.8

- update functions root.xxx to work with WinXP

VERSION xps-0.99.3

- package can now be built for Windows XP

- added possibility to add current date and/or time to root filename

- added function existsROOTFile

- updated vignette xps.Snw

VERSION xps-0.4.3

- new method rawCELName() to get the names of the imported CEL-files

VERSION xps-0.4.2

- add support to import generic (calvin) CEL-files

- update method volcanoplot

VERSION xps-0.4.1

- change DESCRIPTION

- add method volcanoplot

- correct update bug in xpsUniFilter

VERSION xps-0.4.0

- import.data: import CEL-files from different directories

- update DESCRIPTION, NAMESPACE

- add possibility to apply non-sepcific filters and univariate filters

- add S4 classes Filter, PreFilter, UniFilter

- add S4 classes FilterTreeSet and AnalysisTreeSet

yaqcaffy
--------

Changes in version 1.15.1:

- remove url from DESCRIPTION <2011-12-23 Fri>


Packages removed from the release
=================================

The following packages are no longer in the release: 
edd, ontoTools, GeneR, RMAGEML, RTools4TB


[bioc]: http://bioconductor.org
[install]: http://bioconductor.org/install/
[ami]: http://bioconductor.org/help/bioconductor-cloud-ami/
