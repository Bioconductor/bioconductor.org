October 3, 2012

Bioconductors:

We are pleased to announce Bioconductor 2.11, consisting of 610
software packages and more than 650 up-to-date annotation packages.
There are 58 new software packages, and many updates and improvements
to existing packages; Bioconductor 2.11 is compatible with R 2.15.1,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS.  This
release includes an updated Bioconductor
[Amazon Machine Image](http://bioconductor.org/help/bioconductor-cloud-ami/).
Visit [http://bioconductor.org](http://bioconductor.org)
for details and downloads.

Contents
--------

* Getting Started with Bioconductor 2.11
* New Software Packages
* NEWS from new and existing packages
* Packages removed from the release

Getting Started with Bioconductor 2.11
======================================

To install Bioconductor 2.11:

1. Install R 2.15.1.  Bioconductor 2.11 has been designed expressly
   for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/).

3. If you have already been using Bioconductor 2.10 (with R-2.15), you
   can upgrade as follows:

    source("http://bioconductor.org/biocLite.R")
    
    biocLite("BiocUpgrade")




New Software Packages
=====================

There are 58 new packages in this release of Bioconductor.

- agilp: provides a pipeline for the low-level analysis of gene
  expression microarray data, primarily Agilent data

- annmap: annmap provides annotation mappings for Affymetrix exon
  arrays and coordinate based queries to support deep sequencing data
  analysis. Database access is hidden behind the API which provides a
  set of functions such as genesInRange(), geneToExon(), exonDetails(),
  etc. Functions to plot gene architecture and BAM file data are also
  provided. Underlying data are from Ensembl.

- AnnotationForge: Provides code for generating Annotation packages and
  their databases.  Packages produced are intended to be used with
  AnnotationDbi.

- bigmemoryExtras: This package defines a "BigMatrix" ReferenceClass
  which adds safety and convenience features to the
  filebacked.big.matrix class from the bigmemory package. BigMatrix
  protects against segfaults by monitoring and gracefully restoring the
  connection to on-disk data and it also protects against accidental
  data modification with a filesystem-based permissions system. We
  provide utilities for using BigMatrix-derived classes as assayData
  matrices within the Biobase package's eSet family of classes.
  BigMatrix provides some optimizations related to attaching to, and
  indexing into, file-backed matrices with dimnames. Additionally, the
  package provides a "BigMatrixFactor" class, a file-backed matrix with
  factor properties.

- bsseq: Tools for analyzing and visualizing bisulfite sequencing data

- cancerclass: The classification protocol starts with a feature
  selection step and continues with nearest-centroid classification.
  The accurarcy of the predictor can be evaluated using training and
  test set validation, leave-one-out cross-validation or in a multiple
  random validation protocol. Methods for calculation and visualization
  of continuous prediction scores allow to balance sensitivity and
  specificity and define a cutoff value according to clinical
  requirements.

- ChIPXpress: ChIPXpress takes as input predicted TF bound genes from
  ChIPx data and uses a corresponding database of gene expression
  profiles downloaded from NCBI GEO to rank the TF bound targets in
  order of which gene is most likely to be functional TF target.

- chroGPS: We provide intuitive maps to visualize the association
  between genetic elements, with emphasis on epigenetics. The approach
  is based on Multi-Dimensional Scaling. We provide several sensible
  distance metrics, and adjustment procedures to remove systematic
  biases typically observed when merging data obtained under different
  technologies or genetic backgrounds.

- CNORdt: This add-on to the package CellNOptR handles time-course
  data, as opposed to steady state data in CellNOptR. It scales the
  simulation step to allow comparison and model fitting for time-course
  data. Future versions will optimize delays and strengths for each
  edge.

- CNORfuzzy: This package is an extension to CellNOptR.  It contains
  additional functionality needed to simulate and train a prior
  knowledge network to experimental data using constrained fuzzy logic
  (cFL, rather than Boolean logic as is the case in CellNOptR).
  Additionally, this package will contain functions to use for the
  compilation of multiple optimization results (either Boolean or cFL).

- CNORode: ODE add-on to CellNOptR

- CorMut: CorMut provides functions for computing kaks for individual
  sites or specific amino acids and detecting correlated mutations
  among them. Two methods are provided for detecting correlated
  mutations ,including conditional selection pressure and mutual
  information. The computation consists of two steps: First, the
  positive selection sites are detected; Second, the mutation
  correlations are computed among the positive selection sites. Note
  that the first step is optional. Meanwhile, CorMut facilitates the
  comparison of the correlated mutations between two conditions by the
  means of correlated mutation network.

- DBChIP: DBChIP detects differentially bound sharp binding sites
  across multiple conditions, with or without matching control samples.

- ddgraph: Distinguish direct from indirect interactions in gene
  regulation and infer combinatorial code from highly correlated
  variables such as transcription factor binding profiles. The package
  implements the Neighbourhood Consistent PC algorithm (NCPC) and draws
  Direct Dependence Graphs to represent dependence structure around a
  target variable. The package also provides a unified interface to
  other Graphical Modelling (Bayesian Network) packages for
  distinguishing direct and indirect interactions.

- DeconRNASeq: DeconSeq is an R package for deconvolution of
  heterogeneous tissues based on mRNA-Seq data. It modeled expression
  levels from heterogeneous cell populations in mRNA-Seq as the
  weighted average of expression from different constituting cell types
  and predicted cell type proportions of single expression profiles.

- DirichletMultinomial: Dirichlet-multinomial mixture models can be
  used to describe variability in microbial metagenomic data. This
  package is an interface to code originally made available by Holmes,
  Harris, and Quince, 2012, PLoS ONE 7(2): 1-15, as discussed further
  in the man page for this package, ?DirichletMultinomial.

- DSS: DSS is an R library performing the differential expression
  analysis for RNA-seq count data. DSS implements a new dispersion
  shrinkage method to estimate the gene-specific biological variance.
  Extensive simulation results showed that DSS performs favorabily
  compared to DESeq and edgeR when the variation of biological
  variances is large.

- EasyqpcR: This package is based on the qBase algorithms published by
  Hellemans et al. in 2007. The EasyqpcR package allows you to import
  easily qPCR data files as described in the vignette. Thereafter, you
  can calculate amplification efficiencies, relative quantities and
  their standard errors, normalization factors based on the best
  reference genes choosen (using the SLqPCR package), and then the
  normalized relative quantities, the NRQs scaled to your control and
  their standard errors. This package has been created for
  low-throughput qPCR data analysis.

- flowPeaks: A fast and automatic clustering to classify the cells into
  subpopulations based on finding the peaks from the overall density
  function generated by K-means.

- flowQB: flowQB is a fully automated R Bioconductor package to
  calculate automatically the detector efficiency (Q), optical
  background (B), and electronic noise.

- fmcsR: The fmcsR package introduces an efficient maximum common
  substructure (MCS) algorithms combined with a novel matching strategy
  that allows for atom and/or bond mismatches in the substructures
  shared among two small molecules. The resulting flexible MCSs (FMCSs)
  are often larger than strict MCSs, resulting in the identification of
  more common features in their source structures, as well as a higher
  sensitivity in finding compounds with weak structural similarities.
  The fmcsR package provides several utilities to use the FMCS
  algorithm for pairwise compound comparisons, structure similarity
  searching and clustering.

- FunciSNP: FunciSNP integrates information from GWAS, 1000genomes and
  chromatin feature to identify functional SNP in coding or non-coding
  regions.

- gCMAP: The gCMAP package provides a toolkit for comparing
  differential gene expression profiles through gene set enrichment
  analysis. Starting from normalized microarray or RNA-seq gene
  expression values (stored in lists of ExpressionSet and CountDataSet
  objects) the package performs differential expression analysis using
  the limma or DESeq packages. Supplying a simple list of gene
  identifiers, global differential expression profiles or data from
  complete experiments as input, users can use a unified set of several
  well-known gene set enrichment analysis methods to retrieve
  experiments with similar changes in gene expression. To take into
  account the directionality of gene expression changes, gCMAPQuery
  introduces the SignedGeneSet class, directly extending GeneSet from
  the GSEABase package.  To increase performance of large queries,
  multiple gene sets are stored as sparse incidence matrices within
  CMAPCollection eSets. gCMAP offers implementations of 1.  Fisher's
  exact test (Fisher, J R Stat Soc, 1922) 2. The "connectivity map"
  method (Lamb et al, Science, 2006) 3. Parametric and non-parametric
  t-statistic summaries (Jiang & Gentleman, Bioinformatics, 2007) and
  4. Wilcoxon / Mann-Whitney rank sum statistics (Wilcoxon, Biometrics
  Bulletin, 1945) as well as wrappers for the 5. camera (Wu & Smyth,
  Nucleic Acid Res, 2012) 6. mroast and romer (Wu et al,
  Bioinformatics, 2010) functions from the limma package. All methods
  return CMAPResult objects, an S4 class inheriting from
  AnnotatedDataFrame, containing enrichment statistics as well as
  annotation data and providing simple high-level summary plots.

- GeneNetworkBuilder: Appliation for discovering direct or indirect
  targets of transcription factors using ChIP-chip or ChIP-seq, and
  microarray or RNA-seq gene expression data. Inputting a list of genes
  of potential targets of one TF from ChIP-chip or ChIP-seq, and the
  gene expression results, GeneNetworkBuilder generates a regulatory
  network of the TF.

- gmapR: GSNAP and GMAP are a pair of tools to align short-read data
  written by Tom Wu.  This package provides convenience methods to work
  with GMAP and GSNAP from within R. In addition, it provides methods
  to tally alignment results on a per-nucleotide basis using the
  bam_tally tool.

- hapFabia: A package to identify rare and short haplotype clusters in
  large sequencing data by FABIA biclustering. Individuals that
  inherited a particular DNA segment from the same founder constitute a
  haplotype cluster by sharing minor alleles of single nucleotide
  variants (SNVs) that tag/mark this segment. Knowledge of haplotype
  clusters are relevant for phasing of genotyping data, association
  studies, and for population genetics, where they shed light on the
  evolutionary history of humans. The package supports VCF formats, is
  based on sparse matrix operations, and provides visualization of
  haplotype clusters in different formats.

- HMMcopy: Corrects GC and mappability biases for readcounts (i.e.
  coverage) in non-overlapping windows of fixed length for single whole
  genome samples, yielding a rough estimate of copy number for furthur
  analysis.  Designed for rapid correction of high coverage whole
  genome tumour and normal samples.

- hpar: A simple interface to and data from the Human Protein Atlas
  project.

- HTSeqGenie: A software package to analyse high-throughput sequencing
  experiments

- iPAC: iPAC is a novel tool to identify somatic amino acid mutation
  clustering within proteins while taking into account protein
  structure.

- KEGGprofile: KEGGprofile is an annotation and visualization tool
  which integrated the expression profiles and the function annotation
  in KEGG pathway maps. The multi-types and multi-groups expression
  data can be visualized in one pathway map. KEGGprofile facilitated
  more detailed analysis about the specific function changes inner
  pathway or temporal correlations in different genes and samples.

- lmdme: linear ANOVA decomposition of Multivariate Designed
  Experiments implementation based on limma lmFit. Features: i)
  Flexible formula type interface, ii) Fast limma based implementation,
  iii) p and F values over deflacted coefficients and iv) ploting
  functions for PCA and PLS

- matchBox: The matchBox package enables comparing ranked vectors of
  features, merging multiple datasets, removing redundant features,
  using CAT-plots and Venn diagrams, and computing statistical
  significance.

- methyAnalysis: The methyAnalysis package aims for the DNA methylation
  data analysis and visualization. A new class is defined to keep the
  chromosome location information together with the data. The current
  version of the package mainly focus on analyzing the Illumina
  Infinium methylation array data, but most methods can be generalized
  to other methylation array or sequencing data.

- MiRaGE: The package contains functions for inferece of target gene
  regulation by miRNA, based on only target gene expression profile.

- MotifDb: More than 2000 annotated position frequency matrices from
  five public source, for multiple organisms

- motifStack: motifStack is a package that is able to draw amino acid
  sequence as easy as to draw DNA/RNA sequence. motifStack provides the
  flexibility for users to select the font type and symbol colors.
  motifStack is designed for graphical representation of multiple
  motifs.

- networkBMA: An extension of Bayesian Model Averaging (BMA) for
  network construction using time series gene expression data. Includes
  assessment functions and sample test data.

- NOISeq: Analysis of RNA-seq expression data or other similar kind of
  data. Exploratory plots to evualuate saturation, count distribution,
  expression per chromosome, type of detected features, features
  length, etc. Differential expression between two experimental
  conditions with no parametric assumptions.

- OrganismDbi: The package enables a simple unified interface to
  several annotation packages each of which has its own schema by
  taking advantage of the fact that each of these packages implements a
  select methods.

- OSAT: A sizable genomics study such as microarray often involves the
  use of multiple batches (groups) of experiment due to practical
  complication. To minimize batch effects, a careful experiment design
  should ensure the even distribution of biological groups and
  confounding factors across batches. OSAT (Optimal Sample Assignment
  Tool) is developed to facilitate the allocation of collected samples
  to different batches. With minimum steps, it produces setup that
  optimizes the even distribution of samples in groups of biological
  interest into different batches, reducing the confounding or
  correlation between batches and the biological variables of interest.
  It can also optimize the even distribution of confounding factors
  across batches. Our tool can handle challenging instances where
  incomplete and unbalanced sample collections are involved as well as
  ideal balanced RCBD. OSAT provides a number of predefined layout for
  some of the most commonly used genomics platform.

- PADOG: This package implements a general purpose gene set analysis
  method called PADOG that downplays the importance of genes that apear
  often accross the sets of genes to be analyzed. The package provides
  also a benchmark for gene set analysis methods in terms of
  sensitivity and ranking using 24 public datasets from
  KEGGdzPathwaysGEO package.

- PWMEnrich: Asses the enrichment of already known PWMs (e.g. from
  JASPAR) in DNA sequences. Motif hits in a sequence or DNA region are
  considered together and P-values derived for their joint pattern. The
  package implements multiple algorithms, including fixed-threshold
  (Z-score) and threshold-free (Lognormal normalization and Clover)
  methods. The main goal is to identify a set of transcription factors
  that most likely bind to a single sequence, group of sequences, or
  show significantly different binding affinity between two sets of
  sequences.

- Rcade: Rcade (which stands for "R-based analysis of ChIP-seq And
  Differential Expression") is a tool for integrating ChIP-seq data
  with differential expression summary data, through a Bayesian
  framework. A key application is in identifing the genes targeted by a
  transcription factor of interest - that is, we collect genes that are
  associated with a ChIP-seq peak, and differential expression under
  some perturbation related to that TF.

- ReportingTools: The ReportingTools software package enables users to
  easily display reports of analysis results generated from sources
  such as microarray and sequencing data.  The package allows users to
  create HTML pages that may be viewed on a web browser such as Safari,
  or in other formats readable by programs such as Excel.  Users can
  generate tables with sortable and filterable columns, make and
  display plots, and link table entries to other data sources such as
  NCBI or larger plots within the HTML page.  Using the package, users
  can also produce a table of contents page to link various reports
  together for a particular project that can be viewed in a web
  browser.

- RGalaxy: Given an R function and its manual page, make the documented
  function available in Galaxy.

- Risa: The Investigation / Study / Assay (ISA) tab-delimited format is
  a general purpose framework with which to collect and communicate
  complex metadata (i.e. sample characteristics, technologies used,
  type of measurements made) from experiments employing a combination
  of technologies, spanning from traditional approaches to
  high-throughput techniques. Risa allows to access metadata/data in
  ISA-Tab format and build Bioconductor data structures. Currently,
  data generated from microarray, flow cytometry and metabolomics-based
  (i.e. mass spectrometry) assays are supported.  The package is
  extendable and efforts are undergoing to support metadata associated
  to proteomics assays.

- RMassBank: Workflow to process tandem MS files and build MassBank
  records. Functions include automated extraction of tandem MS spectra,
  formula assignment to tandem MS fragments, recalibration of tandem MS
  spectra with assigned fragments, spectrum cleanup, automated
  retrieval of compound information from Internet databases, and export
  to MassBank records.

- rols: This package allows to query EBI's Ontology Lookup Service
  (OLS) using Simple Object Access Protocol (SOAP).

- rSFFreader: rSFFreader reads sequence, qualities and clip point
  values from sff files generated by Roche 454 and Life Sciences Ion
  Torrent sequencers. The plan is to also write out sff files and to
  read in flowgrams with some utils

- SCAN.UPC: SCAN is a microarray normalization method to facilitate
  personalized-medicine workflows. Rather than processing microarray
  samples as groups, which can introduce biases and present logistical
  challenges, SCAN normalizes each sample individually by modeling and
  removing probe- and array-specific background noise using only data
  from within each array. (The Universal Probability of expression
  Codes (UPC) method is an extension of SCAN and will be added to this
  package soon.)

- staRank: Detecting all relevant variables from a data set is
  challenging, especially when only few samples are available and data
  is noisy. Stability ranking provides improved variable rankings of
  increased robustness using resampling or subsampling.

- synapter: The synapter package provides functionality to reanalyse
  label-free proteomics data acquired on a Synapt G2 mass spectrometer.
  One or several runs, possibly processed with additional ion mobility
  separation to increase identification accuracy can be combined to
  other quantitation files to maximise identification and quantitation
  accuracy.

- TransView: This package provides efficient tools to generate, access
  and display read densities of sequencing based data sets such as from
  RNA-Seq and ChIP-Seq.

- triform: The Triform algorithm uses model-free statistics to identify
  peak-like distributions of TF ChIP sequencing reads, taking advantage
  of an improved peak definition in combination with known profile
  characteristics.

- VariantTools: Tools for Tools for detecting, filtering, calling,
  comparing and plotting variants.

- waveTiling: This package is designed to conduct transcriptome
  analysis for tiling arrays based on fast wavelet-based functional
  models.

NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages. The following package NEWS is available:

a4Base
------

Changes in version 1.4.1:

- plot1gene now coerces the 'groups' variable to a factor in plot1gene

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

Changes in version 1.29.13 (2012-09-26):

- Added argument 'cdf=FALSE' to createCel(). Note, the previous
  implementation corresponded to cdf=TRUE.

- ROBUSTNESS: Now createCel() validates/sets CEL header field 'total'
  based on 'cols' and 'rows'.

- ROBUSTNESS: Added a system test for validating that the package can
  write and read a CEL.  The test is spawning of another R process so
  that the test is robust against core dumps.

- Bumped up package dependencies.

Changes in version 1.29.12 (2012-09-12):

- DOCUMENTATION: Updated one Rd link.

Changes in version 1.29.11 (2012-09-01):

- Added argument 'aliases' to arrangeCelFilesByChipType(), e.g.
  arrangeCelFilesByChipType(..., aliases=c("Focus"="HG-Focus")).

- BUG FIX: arrangeCelFilesByChipType(pathnames) assumed 'pathnames'
  were files in the current directory.

Changes in version 1.29.10 (2012-08-29):

- Updated some internal files used solely for maintainance.

Changes in version 1.29.9 (2012-08-29):

- BUG FIX: The move to Fusion SDK 1.1.2 caused the package to not
  compile on Windows.

Changes in version 1.29.8 (2012-08-14):

- Upgraded to Fusion SDK 1.1.2.

Changes in version 1.29.7 (2012-08-14):

- Rearranged patchdir.

Changes in version 1.29.6 (2012-06-26):

- Same updates as in v1.28.1.

Changes in version 1.29.5 (2012-06-19):

- Added arrangeCelFilesByChipType() for moving CEL files to
  subdirectories named according to their chip types, which can be
  useful when for instance downloading GEO data sets.

Changes in version 1.29.4 (2012-06-14):

- readPgfEnv(..., indices=NULL) no longer gives a warning.

- Updated the error messages for the CLF and PGF parsers.

Changes in version 1.29.3 (2012-05-22):

- Now system test tests/testWriteAndReadEmptyCdf.R generates an error
  that is detected and reported by R CMD check.

Changes in version 1.29.2 (2012-05-22):

- GENERALIZATION: Now system tests that launch another R process no
  longer assumes R is on the OS's search path.

- ROBUSTNESS/CRAN POLICY: readCel() and readCelUnits() are no longer
  calling .Internal(qsort(...)).

Changes in version 1.29.1 (2012-05-18):

- Replaced several throw() with stop(), because the former assumes that
  R.methodsS3 is loaded, which it may not be.

- ROBUSTNESS: Added a system test for validating that the package can
  write and read a CDF.  The test is spawning of another R process so
  that the test is robust against core dumps.

Changes in version 1.29.0 (2012-03-30):

- The version number was bumped for the Bioconductor devel version.

Changes in version 1.28.1 (2012-06-26):

- COMPATIBILITY: Now package compile also with gcc/g++ 4.7. Thanks Dan
  Tenenbaum at the Bioconductor Core Team), Fred Hutchinson Cancer
  Research Center, USA for this.

affyPara
--------

Changes in version 1.17.1:

- bugfix for cdfname object in read.affypara, see BioC Mailinglist,
  thanks Tobias

annmap
------

Changes in version 0.99.3:

- BUG ANNMAP-110 XXXInRange( NULL ) used to throw an error.  It now
  returns NULL.

- BUG ANNMAP-111 Fix crash when translation data is incorrect and
  causes us to fall off the translation exon when finding the
  UTR/Coding range.

Changes in version 0.99.2:

- BUG ANNMAP-109 geneToSymbol when passed a vector of gene IDs,
  returned the symbols in the order of the genes start location rather
  than the order the function was called with.

Changes in version 0.99.1:

- BUG ANNMAP-101 ngsBridgePlot called with no genes in range had a huge
  border.

- BUG ANNMAP-106 Empty ngsBridgePlot trace should draw line at 0.

- BUG ANNMAP-107 transcriptToCodingRange was very slow.

Changes in version 0.99.0:

- NEW Initial Release to Bioconductor

- BUG ANNMAP-98 A call to transcriptToTranslatedprobes for transcripts
  in HPRT1 caused a crash

Changes in version 0.9.16:

- BUG ANNMAP-97 geneToExonProbeset threw an error when as.vector=TRUE

Changes in version 0.9.15:

- NEW ANNMAP-94 Added annmapAddConnection method for setting up
  databases.

Changes in version 0.9.14:

- NEW ANNMAP-90 Added generateBridgeData method for helping with BAM
  file imports.

- IMPROVED ANNMAP-93 Refactored out some NCBI->Ensembl seqname mapping
  functions.

- IMPROVED ANNMAP-84 Added ngsBridgePlot examples to the cookbook.

Changes in version 0.9.13:

- IMPROVED ANNMAP-78 Cookbook now appears when vignette(
  package='annmap' ) is called.

Changes in version 0.9.12:

- IMPROVED ANNMAP-74, ANNMAP-75, ANNMAP-76; Fixed warnings and notes
  when using R-2.15

Changes in version 0.9.11:

- IMPROVED ANNMAP-73 Methods generalisedNameToNCBI,
  generalisedNameToEnsembl, seqnameMapping, seqnamesToNCBI and
  seqnamesToEnsembl added to allow easy renaming of the seqnames column
  of GRanges data.

Changes in version 0.9.10:

- IMPROVED ANNMAP-71 Bug in geneToSymbol fixed

Changes in version 0.9.9:

- IMPROVED ANNMAP-63 We now use GenomicRanges instead of IRanges

Changes in version 0.9.8:

- IMPROVED ANNMAP-58 Added 4 pos translation documentation to cookbook

- BUG ANNMAP-60 Webservice tests were failing

- BUG ANNMAP-54 geneToExonProbeset didn't support as.vector

Changes in version 0.9.7:

- BUG ANNMAP-48 Error in proteinCoordsToGenome in annmap

- BUG ANNMAP-47 annmapGenePlot fails to show individual transcripts

- IMPROVED ANNMAP-45 Added a NEWS.Rd file to inst

Changes in version 0.9.6:

- BUG Fixed ANNMAP-44; genomicPlot crashes when called with very small
  overlapped areas

Changes in version 0.9.5:

- IMPROVED ANNMAP-41; All InRange methods are now S3 methods which work
  with character, data.frame or RangedData (see ?annmap.range).

- IMPROVED ANNMAP-42; genomicPlot highlights 'colour' column is now
  'col' as with all other methods.

- BUG ANNMAP-40; genomicPlot crashes if highlights data.frame has
  unexpected columns.

AnnotationForge
---------------

Changes in version 1.9.5:

NEW FEATURES

- An organism package has been added for Streptomyces coelicolor

- extensive overhaul of inparanoid packages means that inparanoid
  packages now match to 100 different organisms

- Extended support for ensembl mappings to yeast and flies.

NOTEWORTHY CHANGES BETWEEN THIS version and 1.3.11

- All chip packages now depend on org packages.  This simplifies the
  schema and also allows for more convenient updating of these packages
  and smaller downloads for users.

- chip package mappings that contain probes which map to multiple
  targets are now hidden by default, with the ability to be exposed
  when required.  See the use of the new toggleProbes() method. * * *
  1.3.11 SERIES NEWS * * *

aroma.light
-----------

Changes in version 1.27.1 (2012-09-12):

- ROBUSTNESS: Replaced an .Internal(psort(...)) call in medianPolish()
  with a call to matrixStats:::.psortKM().

Changes in version 1.27.0 (2012-08-30):

- CLEANUP: Removed weightedMedian(), which has been moved to the
  matrixStats package.

- BACKWARD COMPATIBILITY: Now package depends on the matrixStats (>=
  0.5.2) package, so that weightedMedian() is still available when
  loading this package.  In future releases, matrixStats will be
  downgraded to only be a suggested package.

Changes in version 1.26.1 (2012-08-30):

- BUG FIX: robustSmoothSpline() would not work with most recent R devel
  versions.

- Updated the package dependencies.

Changes in version 1.26.0 (2012-08-19):

- Changed the license of aroma.light to GPL (>= 2) from LGPL (>= 2),
  because some of the implementation was adopted from GPL (>= 2) code,
  i.e. robustSmoothSpline() uses code from stats::smooth.spline().

- R CMD check no longer warns about some examples depending on the
  R.basic package.

Changes in version 1.25.4 (2012-08-19):

- WORKAROUND: Now robustSmoothSpline() robustly locates the proper
  native R fit function for smooth splines, which vary with different
  releases of R.

Changes in version 1.25.3 (2012-04-16):

- Package no longer depends on R.methodsS3, only imports.

Changes in version 1.25.2 (2012-04-16):

- 'R CMD check' no longer complaints about .Internal() calls.

Changes in version 1.25.1 (2012-04-16):

- Added support for fitNaiveGenotypes(..., flavor="fixed").

- GENERALIZATION: Now fitNaiveGenotypes() returns also 'flavor' and
  'tau'.  The latter are the genotype threshholds used by the caller.

- CLEANUP: Dropped argument 'flavor' of callNaiveGenotypes(); it is
  instead passed to fitNaiveGenotypes() via '...'.

Changes in version 1.25.0 (2012-03-30):

- The version number was bumped for the Bioconductor devel version.

arrayQualityMetrics
-------------------

Changes in version 3.13.2:

USER VISIBLE CHANGES

- Using NEWS.Rd

- Added the arguments maxNumArrays and nrColumns to the functions
  aqm.maplot and aqm.spatial.

bigmemoryExtras
---------------

 Version Date Category Text
 0.1.8   <NA> <NA>

Biobase
-------

Changes in version 2.17:

USER VISIBLE CHANGES

- l2e(), previously deprecated, has been made defunct.

- All objects made defunct in previous release cycles have been
  removed.  This includes geneNames, getExpData, eList, reporterNames,
  getBiocRepos, read.exprSet, updateOldMiame, df2pD, read.pD,
  read.phenoData, exprData, exprList, and phenoData.

BitSeq
------

Version: 1.2.0 (27.9.2012)

- IMPORTANT change: the way samples-files are passed to getDE,
  estimateHyperPar, estimateDE changed.  -> instead of providing 2
  vectors of filenames for each condition, the files are passed as a
  list of vectors, each vector containing filenames for one condition
  (allowing use of more than 2 conditions)

- new internal structure (not visible to user)

- estimateExpression has a new convergence criterion which should
  result in producing fewer samples (faster and dropping the use of
  MCMC_scaleReduction and MCMC_samplesNmax flags) ->
  estimateExpressionLegacy uses the original convergence criterion

- library normalization option for getDE, estimateDE,
  estimateHyperPar, getMeanVariance (in form of providing the
  normalization constants, for getting the constants please use edgeR
  or similar)

CAMERA
------

Changes in version 2012-06-12:

- Bugfix in findIsotopes, clears 'vec' must be sorted non-decreasingly
  error

Changes in version 2012-06-11:

- Version 1.13.4

- Add ByteCompile: TRUE

- Bugfix for findIsotopes, clears subscript out of bound error

Changes in version 2012-05-25:

- Version 1.13.2

- First changes for improved isotope detection

Changes in version 2012-04-12:

- Version 1.13.1

- Bugfix for findIsotopes to fix not consecutive isotope label like
  [M]+,[M+2]+ without [M+1]+

ChemmineR
---------

Changes in version 2.10.0:

NEW FEATURES

- Streaming functionality for SDFs enables processing of millions of
  molecules on a laptop

- Fast and memory efficient fingerprint searches with atom pair
  fingerprints or PubChem fingerprints

- Flexible maximum common substructure (MCS) search support provided by
  new fmcs.R add-on package

clusterProfiler
---------------

Changes in version 1.5.1:

- import ggtitle from ggplot2 <2012-09-07, Fri>

- update codes of plot functions accompaning with ggplot2 (version
  0.9.2) <2012-09-06, Thu>

- bug fixed of buildGOmap due to the empty GO annotation query from
  biomaRt <2012-07-18, Wed>

Changes in version 1.5.0:

- bump up version to 1.5.0 for BioC 2.11 devel

- add URL in DESCRIPTION

- update citation and add biocView GeneSetEnrichment <2012-05-09, Wed>

- support zebrafish <2012-05-18, Mon>

cqn
---

Changes in version 1.3:

- Bugfix to the vignette; the two-panel (color) plot on page 6 used CQN
  corrected data as blue points in both panels.  Now the left plot
  shows standard RPKM in blue.  Thanks to Maria Keays <email:
  mkeays@ebi.ac.uk>.

- Small fix to the vignette in the edgeR example, caused by changes to
  edgeR.

- Updated the citations in the vignette and the CITATION file.

- cqn.fixedlength has been removed, using cqn(lengthMethod = "fixed")
  instead.

- A call slot has been added to cqn objects.

- cqn() now accepts count matrices with 1 column or vectors (although
  it makes little sense to use the function on such data).

- Added Questions and Answers to the vignette and moved vignette to
  vignettes dir.

cummeRbund
----------

v1.99.6
  Notes:
    - 'annotation' and "annotation<-" generics were moved to
          BiocGenerics 0.3.2. Now using appropriate generic function,
          but requiring BiocGenerics >= 0.3.2
v1.99.5
  Bugfixes:
    - Added replicates argument to csDistHeat to view distances
          between individual replicate samples.
    - Appropriately distinguish now between 'annotation' (external
          attributes) and features (gene-level sub-features).
    - csHeatmap now has 'method' argument to pass function for any
          dissimilarity metric you desire. You must pass a function
          that returns a 'dist' object applied to rows of a
          matrix. Default is still JS-distance.
    
v1.99.3
  New Features:
    - Added diffTable() method to return a table of differential
          results broken out by pairwise comparison. (more
          human-readable)
    - Added sigMatrix() method to CuffSet objects to draw heatmap
          showing number of significant genes by pairwise comparison
          at a given FDR.
    - A call to fpkm() now emits calculated (model-derived)
          standard deviation field as well.
    - Can now pass a GTF file as argument to readCufflinks() to
          integrate transcript model information into database backend
      * Added requirement for rtracklayer and GenomicFeatures
              packages.
      * You must also indicate which genome build the .gtf was
              created against by using the 'genome' argument to
              readCufflinks.
    - Integration with Gviz:
      * CuffGene objects now have a makeGeneRegionTrack()
              argument to create a GeneRegionTrack() from transcript
              model information
      * Can also make GRanges object
      * ONLY WORKS IF YOU READ .gtf FILE IN WITH readCufflinks()
    - Added csScatterMatrix() and csVolcanoMatrix() method to
          CuffData objects.
    - Added fpkmSCVPlot() as a CuffData method to visualize
          replicate-level coefficient of variation across fpkm range
          per condition.
    - Added PCAplot() and MDSplot() for dimensionality reduction
          visualizations (Principle components, and multi-dimensional
          scaling respectively)
    - Added csDistHeat() to create a heatmap of JS-distances
          between conditions.
    
  Bugfixes:
    - Fixed diffData 'features' argument so that it now does what
          it's supposed to do.
    - added DB() with signature(object="CuffSet") to NAMESPACE
    
  Notes:
    - Once again, there have been modifications to the underlying
          database schema so you will have to re-run
          readCufflinks(rebuild=T) to re-analyze existing datasets.
    - Importing 'defaults' from plyr instead of requiring entire
          package (keeps namespace cleaner).
    - Set pseudocount=0.0 as default for csDensity() and
          csScatter() methods (This prevents a visual bias for genes
          with FPKM <1 and ggplot2 handles removing true zero values).
    
v1.99.2
  Bugfixes:
    - Fixed bug in replicate table that did not apply
          make.db.names to match samples table.
    - Fixed bug for missing values in *.count_tracking files.
    - Now correctly applying make.db.names to
          *.read_group_tracking files.
    - Now correctly allows for empty *.count_tracking and
          *.read_group_tracking files
v1.99.1
  - This represents a major set of improvements and feature
      additions to cummeRbund.
  - cummeRbund now incorporates additional information emitted from
      cuffdiff 2.0 including:
    - run parameters and information.
    - sample-level information such as mass and scaling factors.
    - individual replicate fpkms and associated statistics for all
          features.
    - raw and normalized count tables and associated statistics
          all features.
  
  New Features:
    - Please see updated vignette for overview of new features.
    - New dispersionPlot() to visualize model fit (mean count vs
          dispersion) at all feature levels.
    - New runInfo() method returns cuffdiff run parameters.
    - New replicates() method returns a data.frame of
          replicate-level parameters and information.
    - getGene() and getGenes() can now take a list of any
      tracking_id or gene_short_name (not just gene_ids) to
      retrieve a gene or geneset.
    - Added getFeatures() method to retrieve a CuffFeatureSet
          independent of gene-level attributes.  This is ideal for
          looking at sets of features outside of the context of all
          other gene-related information (i.e. facilitates
          feature-level analysis)
    - Replicate-level fpkm data now available.
    - Condition-level raw and normalized count data now available.
    - repFpkm(), repFpkmMatrix, count(), and countMatrix are new
          accessor methods to CuffData, CuffFeatureSet, and
          CuffFeature objects.
    - All relevant plots now have a logical 'replicates' argument
          (default = F) that when set to TRUE will expose replicate
          FPKM values in appropriate ways.
    - MAPlot() now has 'useCount' argument to draw MA plots using
          count data as opposed to fpkm estimates.
  
  Notes:
    - Changed default csHeatmap colorscheme to the much more
          pleasing 'lightyellow' to 'darkred' through 'orange'.
    - SQLite journaling is no longer disabled by default (The
          benefits outweigh the moderate reduction in load times).
  
  Bugfixes:
    - Numerous random bug fixes to improve consistency and improve
          performance for large datasets.
v1.2.1
  Bugfixes:
    -Fixed bug in CuffFeatureSet::expressionBarplot to make
         compatible with ggplot2 v0.9.
  New Features:
    - Added 'distThresh' argument to findSimilar.  This allows you
          to retrieve all similar genes within a given JS distance as
          specified by distThresh.
    - Added 'returnGeneSet' argument to findSimilar.  [default =
          T] If true, findSimilar returns a CuffGeneSet of genes
          matching criteria (default). If false, a rank-ordered data
          frame of JS distance values is returned.
    - findSimilar can now take a 'sampleIdList' argument. This
          should be a vector of sample names across which the distance
          between genes should be evaluated.  This should be a subset
          of the output of samples(genes(cuff)).
  Notes:
    - Added requirement for 'fastcluster' package.  There is very
          little footprint, and it makes a significant improvement in
          speed for the clustering analyses.
DART
----

 Version Date Category                                  Text
 1.3.1   <NA> Dependency changed from igraph to igraph0

ddgraph
-------

Changes in version 1.0.3:

- Fixed bugs with multiple testing correction. Show actual
  multiple-corrected P-values on edges.

Changes in version 1.0.2:

- Change default representation of conditional independence in DDGraphs
  from arrows to dots

Changes in version 1.0.1:

- Add the toyExample dataset and extend the method introduction in
  vignette

Changes in version 1.0.0:

- Initial release of the package

deepSNV
-------

Changes in version 1.3.3 (2012-09-14):

Updates

- Changed CITATION
- Updated documentation
- Updated biocViews
- Using roxygen2 for generating man pages

Changes in version 1.3.2 (2012-04-10):

Updates

- New devel version, identical to release 1.2.3

Changes in version 1.2.3 (2012-04-10):

Bugfixes

- Fixed Vignette
- Jumped a few numbers due to automated bioc version numbering

DESeq
-----

Changes in version 1.9.12 (2012-09-05):

- added function newCountDataSetFromHTSeqCount

Changes in version 1.9.7 (2012-05-25):

- fixed handling of rank-deficient observations now also in
  estimateDispersions

Changes in version 1.9.6:

- fixed that varianceStabilizingTransformation was missing in exports

Changes in version 1.9.5:

- added function plotDispEsts

Changes in version 1.9.4 (2012-04-24):

- fixed handling of rank-deficient observations

Changes in version 1.9.2 (2012-04-05):

- fixed bug in pooled-CR dispersion estimation

Changes in version 1.9.1 (2012-04-01):

- Added a new function 'varianceStabilizingTransformation'

DEXSeq
------

Changes in version 2012-06-26:

- Now any function relies on the order of the levels of the factors.

Changes in version 2012-05-21:

- More options and flexibilty added to "estimatelog2FoldChanges"

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

Changes in version 1.4.0:

- Plotting

  * dba.plotMA

  * Smooth plots now default

  * Added fold parameter in addition to th (threshold)

  * dba.plotHeatmap

  * Side colorbars added

  * Add support for specifying sample mask to include any subset of
  samples in a contrast plot, including samples that were not in the
  original contrast

  * dba.plotVenn

  * Changed plotter from limma to T. Girke's overLapper

  * Added support for 4-way Venns (also in dba/overlap)

  * dba.plotPCA

  * Add support for specifying sample mask to include any subset of
  samples in a contrast plot, including samples that were not in the
  original contrast

- Peaksets (dba and dba.peakset)

  * Peakset formats

  * narrowPeaks format supported

  * Can override file format, score column, and score orientation
  defaults for supported peak callers

  * Consensus peaksets

  * Added ability to generate sets of consensus peaksets based on
  metadata attributes: for example create consensus peaksets for each
  tissue type and/or condition, or for all unique samples by taking the
  consensus of their replicate peaksets

- Read counting (dba.count)

  * Compute Signal-to-Noise ratio when counting

  * Added bScaleControl to down-scale control reads by default

  * Add option to specify a mask in peak parameter to limit which
  peaksets are used to for a consensus by overlap. Works with new
  consensus peakset options in dba.peakset

  * Remove references to support for SAM files

- Analysis (dba.analyze)

  * edgeR: updated calls to math change sin edgeR; updated vignette and
  references

  * DESeq: updated to work with current DESeq; use pooled-CR dispersion
  estimation method for blocking analysis; update vignette

- Various bug fixes; more informative warnings; update documentation
  including vignette, new examples and cross-referencing in man pages

Changes in version 1.2.3 (2012-09-01):

- more informative warnings and minor bug fixes.

DOSE
----

Changes in version 1.3.2:

- update codes of plot functions accompaning with ggplot2 (version
  0.9.2) <2012-09-06, Thu>

- update cnetplot corresponding to igraph version 0.6 <2012-07-11, Wed>

- parameter showCategory now support term ID vector <2012-07-11, Wed>

- import termSim and combineScores from GOSemSim. <2012-06,14, Thu>

- optimize setReadable <2012-05-09, Wed>

- bug fixed of setReadable method. For those unmapped genes, return the
  original gene ID. <2012-05-15, Tue>

- fill color in barplot for Ontology classification. <2012-05-22, Tue>

- update color scale of cnetplot <2012-06-18, Mon>

Changes in version 1.3.1:

- update color scheme of cnetplot <2012-04-10, Tue>

- add simplot for plotting semantic similarity matrix <2012-04-20, Fri>

- bug fixed of combineScore function for DO semantic similarity matrix
  which containing NA rows of NA coloumns <2012-04-20, Fri>

- export doSim and geneSim functions <2012-04-20, Fri>

DSS
---

Changes in version 0.99-1:

- Fixed a bug in newSeqCountSet in dealing with the input experimental
  designs.

- Fixed a bug in computing local FDR.

- Add options to model the relationship between dispersion and mean
  expression.

Changes in version 0.99:

NEW FEATURES

- Initial release.

EasyqpcR
--------

Changes in version 1.0.0:

- First version of this package.

easyRNASeq
----------

Changes in version 1.3.14:

NEW FEATURES

- easyRNASeq now returns a SummarizedExperiment in an effort to
  consolidate the objects used for Next Generation Sequencing in
  Bioconductor. This is the new default of the count function. The
  count function is a new function to supersed easyRNASeq in the coming
  development version (1.5.x) to consolidate the parameters and output
  of the easyRNASeq function.

BUG FIXES

- corrected a validity check that went permissive.

- changed the print method to display the read length range when
  dealing with variable read lengths rather than every single value.

Changes in version 1.3.13:

BUG FIXES

- Same correction as in the stable version 1.2.5, but for those already
  corrected in version 1.3.3.

Changes in version 1.3.12:

- Providing the 'outputFormat' argument is not necessary anymore, it
  defaults to matrix (i.e. count table).

- Relaxed the gtf file checking. If the gene_name is absent, the
  gene_id is used instead.

- Improved some reporting and remove a bottle-neck occuring when there
  are many sequences in the reference. BUG FIXES

- Ensure that only the matched ranges are returned when reading gapped
  alignments.

- The library size is more exactly calculated and is the number of
  aligned reads.

- Corrected a bug in the validity checking that prevented bam files
  created by different aligners using the same reference to be
  processed as the reference sequences were not ordered in the same
  fashion.

Changes in version 1.3.11:

BUG FIXES

- Fixed a bug in the gtf file handling reported by Mark Robinson.

Changes in version 1.3.10:

- Some vignette discrepancies have been corrected.  Thanks to Richard
  Friedman for spotting them.

- Providing the 'filesDirectory' argument is not necessary anymore, if
  the files to proceed are present in the current directory. Indeed,
  this parameter now defaults to the current directory as can be found
  out using 'getwd()'. BUG FIXES

- Fixed a bug introduced by a change in the IRanges coverage function
  return value.

Changes in version 1.3.9:

- Added the manuscript citation.

- Updated the package version dependencies. BUG FIXES

- Improved the support for reads of different lengths.

- A cosmetic change to report read lengths as well when read files with
  variable read length are processed.

- Corrected a bug and enhanced the loading of gtf annotation files.
  Thanks to Tomasz Kulinski for spotting the issue and providing the
  dataset to reproduce it.

Changes in version 1.3.8:

NEW FEATURES

- Now bam files can be processed in parallel (long time request from
  Wade Davis). If the easyRNASeq argument 'nbCore' is greater than 1 (1
  being the default), then that many core will be used to process the
  read files in parallel.  Pay attention not to use too many cores and
  have enough memory available. The memory load scales up linearly with
  the number of files processed.

Changes in version 1.3.7:

NEW FEATURES

- easyRNASeq now supports read of different lenghts. Thanks to Mark
  Robinson for the toy dataset.

- Added a function that lists existing organism conversion when
  applying the validity checks.

- Added a bp.coverage to the fetchCoverage function that defaults to
  FALSE. To allow for variable length reads, it now returns read
  coverage proportion per bp by default.

- Added additional checks in the .checkArguments internal function.

BUG FIXES

- Not a real bug, but more a consolidation. When an organism is unknown
  and no custom.map is provided, then the validity checks are turned
  off and a warning is emitted.

- Providing the chr.sizes as as list has been deprecated.  Only named
  numeric vector are supported.

- Removed a now useless warning in the .readGffGtf function.

- Modified the RPKM function generic to avoid using a 'protected' word
  as argument: i.e. 'unique' was replaced by 'simplify'

Changes in version 1.3.6:

NEW FEATURES

- It is now possible to pass arguments to list files through the three
  dots. I.e. setting recursive=TRUE is now possible.

BUG FIXES

- Corrected a bug in the .getArguments internal function.

Changes in version 1.3.5:

NEW FEATURES

- bam is now the default format for the easyRNASeq method.

- chromosome sizes are now extracted from the BAM header when the
  'chr.sizes' argument is set to "auto". Thanks to Simon Anders for
  pushing that off my TODO list and the nice implementation.

BUG FIXES

- Adapted to an API change of the edgeR package for estimating the
  tagwise dispersion.

Changes in version 1.3.4:

NEW FEATURES

- Added an additional validity check for chromosome names Thanks to
  Simon Anders for generating a reproducible use-case for that. Same
  change as in the stable version 1.2.3

- Ensure that gtf with non Ensembl ID are correctly parsed as well.

Changes in version 1.3.3:

- Converted the package to use Roxygen2, a Doxygen like in-source
  documentation system for generating the RD and NAMESPACE. The
  original man page were converted using the Rd2roxygen package and the
  resulting in-source documentation manually edited. NEW FEATURES

- Added a type accessor for Genome_intervals object

- Added a coercion to GRangesList from Genome_intervals object BUG
  FIXES

- Adapted to the new arguments of the edgeR estimateTagwiseDisp
  function

- Removed the dispersion.method argument from the plotMeanVar edgeR
  method call as this argument is defunct.

Changes in version 1.3.2:

BUG FIXES

- Corrected a bug that was considering a GTF file as a GFF file. Thanks
  to Simon Anders for spotting this.

Changes in version 1.3.1:

NEW FEATURES

- Added an enhanced read length check (same as stable 1.2.1 change)

Changes in version 1.3.0:

- New development version for Bioconductor 2.11

Changes in version 1.2.5:

BUG FIXES

- Corrected a bug in the condition file name checking.

- When using edgeR, it was not possible to de-activate the drawing of
  the quality assessment plots.

- Some edgeR changes to the API have been ported to the stable R
  version, should not have occured... The following are changes that
  adapt to that new API, changes ported from version the easyRNASeq
  development version 1.3.3...

- Adapted to the new arguments of the edgeR estimateTagwiseDisp
  function

- Removed the dispersion.method argument from the plotMeanVar edgeR
  method call as this argument is defunct.

Changes in version 1.2.4:

- Added the manuscript citation.

- Updated the package version dependencies.

Changes in version 1.2.3:

NEW FEATURES

- Added an additional validity check for chromosome names Thanks to
  Simon Anders for generating a reproducible use-case for that.

- Ensure that gtf with non Ensembl ID are correctly parsed as well.

Changes in version 1.2.2:

BUG FIXES

- Corrected a bug that was considering a GTF file as a GFF file. Thanks
  to Simon Anders for spotting this.

Changes in version 1.2.1:

NEW FEATURES

- Added an enhanced read length check

EBImage
-------

Changes in version 4.0.0:

NEW FEATURES

- 'transpose' function for transposing an image by swapping its spatial
  dimensions

- greyscale functions for computation of the self-complementary top-hat
  (I. Kats)

- a median filter based on Perreault's constant time median filter (J.
  Barry)

SIGNIFICANT USER-VISIBLE CHANGES

- removed all dependencies towards GTK+ and ImageMagick

- replaced the former GTK+ based 'display' function by a new one
  displaying images using either a JavaScript image viewer, or R's
  built-in raster graphics

- 'readImage' and 'writeImage' now rely on 'jpeg', 'png' and 'tiff'
  packages and do not depend on ImageMagick any more

- added support for images containing an alpha channel; both greyscale
  and color images with an alpha channel are stored as a 'colormode =
  Color' Image

- refactored the functions, not using ImageMagick any longer:
  'translate', 'affine', 'rotate', 'resize'

- deprecated: 'blur', 'equalize'

- deprecated: 'drawtext', 'drawfont'

- deprecated: 'getFeatures', 'hullFeatures', 'zernikeMoments',
  'edgeProfile', 'edgeFeatures'

- deprecated: 'haralickFeatures', 'haralickMatrix'

- deprecated: 'moments', 'smoments', 'rmoments', 'cmoments'

- removed 'animate'

- improved 'getFrame': better performance by reassigning array
  dimension only when needed

- modified 'as.raster'

- 'inst/images/lena.gif' is now 'inst/images/lena.png'

- overhauled the testing procedure in 'tests/test.R'

- added 'NEWS.Rd'

BUG FIXES

- 'erode' and 'dilate': incorrect range of loop indices caused memory
  reads from outside the kernel vector

EDASeq
------

Changes in version 1.3:

- Fixed a bug in biasPlot relative to the lwd, xlab and ylab arguments.

- Added a color_code option and changed the behavior of col in
  biasPlot.

- Updated CITATION file.

- Added an option to withinLaneNormalization and
  betweenLaneNormalization to return unrounded values.

- A new way to deal with zero counts by adding a small positive
  constant.

edgeR
-----

Changes in version 3.0.0:

- New chapter in the User's Guide covering a number of common types of
  experimental designs, including multiple groups, multiple factors and
  additive models.  Many other updates to the User's Guide and to the
  help pages.

- New function edgeRUsersGuide() to open the User's Guide in a pdf
  viewer.

- Many functions have made faster by rewriting the core computations in
  C++. This includes adjustedProfileLik(), mglmLevenberg(),
  maximizeInterpolant() and goodTuring().

- New argument verbose for estimateCommonDisp() and
  estimateGLMCommonDisp().

- The trended dispersion methods based on binning and interpolation
  have been rewritten to give more stable results when the number of
  genes is not large.

- The amount by which the tagwise dispersion estimates are squeezed
  towards the global value is now specified in estimateTagwiseDisp(),
  estimateGLMTagwiseDisp() and dispCoxReidInterpolateTagwise() by
  specifying the prior degrees of freedom prior.df instead of the prior
  number of samples prior.n.

- The weighted likelihood empirical Bayes code has been simplified or
  developed in a number of ways. The old functions weightedComLik() and
  weightedComLikMA() are now removed as no longer required.

- The functions estimateSmoothing() and approx.expected.info() have
  been removed as no longer recommended.

- The span used by estimateGLMTagwiseDisp() is now chosen by default as
  a decreasing function of the number of tags in the dataset.

- New method "loess" for the trend argument of estimateTagwiseDisp,
  with "tricube" now treated as a synonym.

- New functions loessByCol() and locfitByCol() for smoothing columns of
  matrix by non-robust loess curves. These functions are used in the
  weighted likelihood empirical Bayes procedures to compute local
  common likelihood.

- glmFit now shrinks the estimated fold-changes towards zero. The
  default shrinkage is as for exactTest().

- predFC output is now on the natural log scale instead of log2.

- mglmLevenberg() is now the default glm fitting algorithm, avoiding
  the occasional errors that occurred previously with mglmLS().

- The arguments of glmLRT() and glmQLFTest() have been simplified so
  that the argument y, previously the first argument of glmLRT, is no
  longer required.

- glmQLFTest() now ensures that no p-value is smaller than what would
  be obtained by treating the likelihood ratio test statistic as
  chisquare.

- glmQLFTest() now treats tags with all zero counts in replicate arrays
  as having zero residual df.

- gof() now optionally produces a qq-plot of the genewise goodness of
  fit statistics.

- Argument null.hypothesis removed from equalizeLibSizes().

- DGEList now longer outputs a component called all.zeros.

- goodTuring() now longer produces a plot. Instead there is a new
  function goodTuringPlot() for plotting log-probability versus
  log-frequency. goodTuring() has a new argument 'conf' giving the
  confidence factor for the linear regression approximation.

- Added plot.it argument to maPlot().

ExiMiR
------

Changes in version 1.99.1:

- NAMESPACE: removing the load of limma and affy that are in
  dependencies

Changes in version 1.99.0:

- inst/doc/fig0.png: new figure

- man/bg.correct.miR.Rd: new help file

- man/createAB.Rd: new help file

- man/NormiR.methods.Rd: new help file

- R/createAB.R: new feature implementation

- R/NormiR.methods: new file for giving availables methods

- All the others files have been updated according new features
  implementation

fabia
-----

Changes in version 2.3.1:

NEW FEATURES

- Getters and setters for class Factorization

flowCore
--------

- add new classes "filters", "filtersList" to allow flowViz to plot
  multiple filters/gates for one flowFrame
- add argument "emptyValue" to read.FCS API so that parser can still
  work correctly when either cases below occurs :
  1. there is double-delimiter in keyword values (sometime
     like\n\\\\c:\\\\path\\\\...)
  2. there is empty keyword
     value\n(\\\\keyword1\\\\value1\\\\keyword2\\\\\\\\keyword3)
- fix the bug that malformed spillover matrix in write.FCS


flowViz
-------

CHANGES IN VERSION 1.21.1

1.add modified lattice theme to flowViz and change the default color
  scheme for non-smoothed xyplot
2.add stat=TRUE to display population % in xyplot and add abs=FALSE
  and pos=0.5 to control the position of gate labels
3.made change to prepanel.xyplot.flowset so that it return an empty
  list instead of NULL value for empty panels.This was causing the
  error thrown by lattice:::limits.and.aspect which calculates the
  scales for each panel and expects non-null return value from
  prepanel function
4.remove the old src and R code for hexbin and add hexbin package
  based hexagon plot support within panel.xyplot.flowframe
5.add binTrans argument to xyplot that gets passed to hexin to
  transform the raw counts. sqrt is the default,NULL value means no
  transformation.
6.add new classes "filters", "filtersList" to allow flowViz to plot
  multiple filters/gates for one flowFrame

fmcsR
-----

Changes in version 0.99.0:

- fmcsR submitted to Bioconductor

FunciSNP
--------

Changes in version 0.99.0 (2012-05-25):

- Updated package to address issues and comments by BioC curator.

Changes in version 0.2.0 (2012-05-18):

- Initial Bioconductor release. New citation included.

Changes in version 0.1.0 (2012-02-12):

- Created. Initial build and release of FunciSNP.

GeneAnswers
-----------

Changes in version 1.14:

NEW FEATURES

- Multigroup concept-gene analysis html report supports interactive
  network with cytoscape web support as well as original fixed images

- caBIO pathway and REACTOME pathway are included with xml queries,
  therefore, internet access is required.

- The total number of pooled genes in Hypergeometric test can be set by
  the amount of genes in annotation library or the total annoted genes
  in the given species.

- Add gene annotation summarization functions.  Please check
  geneFunctionSummarize.pdf for more details.

GeneNetworkBuilder
------------------

Changes in version 1.0.0:

NEW FEATURES

- buildNetwork() to build the network by binding list and interaction
  map

- filterNetwork() to filter the network by expression data

- polishNetwork() to generate graphNEL object for visualization

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)


GenomicFeatures
---------------

Changes in version 1.10:

NEW FEATURES

- Add makeTranscriptDbFromGFF().  Users can now use GFF files to make
  TranscriptDb resources.

- Add *restricted* "seqinfo<-" method for TranscriptDb objects. It only
  supports replacement of the sequence names (for now), i.e., except
  for their sequence names, Seqinfo objects 'value' (supplied) and
  'seqinfo(x)' (current) must be identical.

- Add promoters() and getPromoterSeq().

- Add 'reassign.ids' arg (FALSE by default) to makeTranscriptDb().

SIGNIFICANT USER-VISIBLE CHANGES

- Updated vignette.

- Improve how makeTranscriptDbFromUCSC() and
  makeTranscriptDbFromBiomart() assign internal ids (see commit 65144
  for the details).

- 2.5x speedup of fiveUTRsByTranscript() and threeUTRsByTranscript().

DEPRECATED AND DEFUNCT

- Are now defunct: transcripts_deprecated(), exons_deprecated(), and
  introns_deprecated().

- Deprecate loadFeatures() and saveFeatures() in favor of loadDb() and
  saveDb(), respectively.

BUG FIXES

- Better handling of BioMart data anomalies.


GenomicRanges
-------------

Changes in version 1.10.0:

NEW FEATURES

- SummarizedExperiment gains direct GRanges / GRangesList interface to
  rowData.

- Add "distanceToNearest" method for GenomicRanges objects.

- SummarizedExperiment class can now be subset by row when there are no
  'columns', and by column when there are no 'rows'.

- Add 'drop.D.ranges' argument to coverage,GappedAlignments and
  coverage,GappedAlignmentPairs methods.

- findOverlaps() now supports 'select="last"' and 'select="arbitrary"'
  (in addition to 'select="all"' and 'select="first"') on GenomicRanges
  objects.

- summarizeOverlaps(..., mode="IntersectionStrict") now handles
  circular chromosomes. A warning is issued and circular chromosomes in
  'reads' are omitted from counting.

- Add disjoin,GRangesList method.

- Add findSpliceOverlaps() for identifyng ranges (reads) that are
  compatible with a specific transcript isoform (the non-compatible
  ranges are analyzed for the presence of novel splice events).

- Add ngap,GappedAlignmentPairs method.

- Add introns() generic with methods for GappedAlignments and
  GappedAlignmentPairs objects.

- No more arbitrary max of 3 gaps per read in
  isCompatibleWithSplicing() and isCompatibleWithSkippedExons().

- Add findCompatibleOverlaps() and countCompatibleOverlaps().

- Passing '...' down through as.data.frame(GRanges, ...) so user can
  tweak stringsAsFactors default for metadata columns.

- Add extractSteppedExonRanks(), extractSpannedExonRanks() and
  extractQueryStartInTranscript() utilities (work with single- and
  paired-end reads).

- Add 'flip.query.if.wrong.strand' arg (FALSE by default) to
  "encodeOverlaps" method for GRangesList objects.

- Add makeSeqnameIds() low-level utility.

SIGNIFICANT USER-LEVEL CHANGES

- SummarizedExperiment rowData and assays operations have significant
  performance improvements.

- mcols() is now the preferred way (over elementMetadata() or values())
  to access the metadata columns of a GenomicRanges, GRangesList,
  GappedAlignments, GappedAlignmentPairs, SummarizedExperiment object,
  or any Vector object. elementMetadata() and values() might go away at
  some point in the (not so close) future.

- Add "$" and "$<-" methods for GenomicRanges *only*. Provided as a
  convenience and as the result of strong popular demand. Note that
  those methods are not consistent with the other "$" and "$<-" methods
  in the IRanges/GenomicRanges infrastructure, and might confuse some
  users by making them believe that a GenomicRanges object can be
  manipulated as a data.frame-like object. It is therefore recommended
  to use them only interactively, and their use in scripts or packages
  is discouraged. For the latter, use 'mcols(x)$name' instead of
  'x$name'.

- No more warning when doing as(x, "GRanges") on a RangedData object
  with no "strand" column.

- Refactor "[" method for GenomicRanges objects. The new implementation
  always preserves the names of the selected elements instead of trying
  to return a GenomicRanges object with unique names. This new behavior
  is consistent with subsetting of ordinary vectors and other Vector
  objects defined in IRanges/GenomicRanges. Also modify "seqselect"
  method for GenomicRanges objects so it also preserves the names of
  the selected elements (and thus remains consistent with new behavior
  of "[" method for GenomicRanges objects).

- No more names on the integer vector returned by "ngap" method for
  GappedAlignments objects.

- Many improvements to the "Overlap encodings" vignette.

- Remove 'param' argument from summarizeOverlaps() generic.

DEPRECATED AND DEFUNCT

- Defunct previously deprecated grg() function.

- Defunct previously deprecated countGenomicOverlaps() generic and
  methods.

BUG FIXES

- Fix several issues with "precede", "follow", "nearest", and
  "distance" methods for GenomicRanges objects.

- Fix bug in summarizeOverlaps(..., ignore.strand=TRUE).

- 6x speedup (and a 6x memory footprint reduction) or more when using
  encodeOverlaps() on big GRangesList objects.

- Fix bug in renameSeqlevels() wrt order of rename vector.

- Fix bug in selectEncodingWithCompatibleStrand().


genoset
-------

1.9.8 GRanges everywhere!  GenoSet now supports GRanges in the locData
slot.  All functions that take RangedData now also take GRanges.  I
have unified the API for GRanges, RangedData, and GenoSet to the point
that GenoSet classes and the functions in the package are agnostic to
the type of range object.  I have not, however, fixed the contentious
issue of using the "$" operator with GRanges to access
elementMetadata.

1.9.10 Subsetting by location now only with GRanges and RangedData.
Dropped RangesList to avoid weird errors about the
RangedDataOrRangesListOrGRanges class union.  Apparently the
RangedDataOrGRanges class union is fine. I think RangesLists are not
used often anyway.

1.9.11 GenoSet creation and featureNames<- no longer do make.names.

1.9.12 GenoSet creation and sampleNames<- no longer do make.names.

GGBase
------

Changes in version 3.20.0:

- MAFfilter now treats lower, upper as a left-open interval, i.e.,
  retains x if lower < x <= upper; previously was closed at both ends

ggbio
-----

Changes in version 1.5.16:

NEW FEATURES

- ggplot generic method added.

- mold generic method added for molding object to data.frame.

- support ggplot(data) + stat_* style, original data being kept.

- tracks function updated, API is enhenced, utilities could control
  attributes of plots and trakcs.

- autoplot now support: VCF, SummarizedExperiments, matrix, where when
  genomic position is provided, options to visualize a heatmap sitting
  on the genomic position.

- theme could define track based themes.

- ideogram: support + xlim method, when embeded with tracks,
  automatically update zoomed region.

- pheno.plot added to SuumarizedExperiemnts and ExpressionSet.

SIGNIFICANT USER-LEVEL CHANGES

- align.plots is deprecated, alignPlots created

Notes

- updated website for ggbio: http://tengfei.github.com/ggbio manuals
  and vignettes and paper added

GGtools
-------

Changes in version 4.6:

- The primary tools for one-population analyses are best.cis.eQTLs and
  transScores.  Multipopulation analyses are handled with
  meta.best.cis.eQTLs and meta.transScores.

- High volume genotype data has been addressed by packaging
  ExpressionSet and chromosome-specific SnpMatrix instances; requests
  for expression plus genotype data are directed to packages mediated
  through GGBase::getSS.

- Two species of data filtering parameters that may be used jointly in
  the primary tools are exFilter, which operates on expression
  component prior to any analyses, and smFilter, which operates on the
  entire smlSet.  exFilter may be used to isolate samples of interest
  early in the workflow, for example when an expression plus genotype
  package includes samples from distinct tissues on the same
  individuals.

- june 2012: exFilter facility properly handled in best.cis.eQTLs.mchr

- may 2012: best.cis.eQTLs has getDFFITS option

GOSemSim
--------

Changes in version 1.15.3:

- remove all the S4 classes and methods <2012-09-12, Wed>

- add progress bar for mgeneSim <2012-09-12, Wed>

- re-implement calculating semantic values in Wang's method
  <2012-09-12, Wed>

- update IC data for next release <2012-09-12, Wed>

- bug fixed in getSV <2012-09-13, Thu>

Changes in version 1.15.2:

- re-implement gene2GO <2012-09-7, Fri>

- information content based methods implemented in c++ <2012-09-5, Wed>

Changes in version 1.15.1:

- export termSim, which can be used in other ontological semantic
  similarity measurement <2012-06-14, Thu>

- update vignette. <2012-06-14, Thu>

GSEABase
--------

Changes in version 1.19:

NEW FEATURES

- Added UniprotIdentifier class

Gviz
----

Changes in version 1.2.0:

NEW FEATURES

- A SequenceTrack class has been added to draw genomic sequence
  information on a Gviz plot. Possible inputs for the track are
  DNAStringSet objects or directly from BSgenome packages.

- GeneRegionTracks can now deal with coding and non-coding regions by
  means of the feature property in combination with the thinBoxFeature
  display parameter.

- StackedTracks now have a new display parameter 'reverseStacking'
  which reverts the horizontal ordering of stacked items. If set to
  TRUE, the lowest items are moved to the top of the stack, and vice
  versa.

SIGNIFICANT USER-VISIBLE CHANGES

- Updated the show methods for most tracks to give more meaningful and
  more compact information about the track's content. Availablability
  of data on other chromosomes than the currently active one should now
  be indicated.

- IdeogramTracks can now be constructed from a cytoband table via the
  new bands argument in the constructor.

- AnnotationTrack objects now by default draw connecting lines in a
  light gray color. This feature can be controlled via the col.line
  display parameter.

- Sliding window summarization can now deal with NA values.

- Exporting drawGD from the name space now to allow for sub-classing of
  GdObjects in other packages.

- When building GeneRegionTracks from TrasncriptDb objects, the
  information about UTRs and coding regions is now retained.

BUG FIXES

- When zooming into the emty space between two grouped features, the
  connecting line will now be plotted for all classes inheriting from
  AnnotationTrack.

- An error in calculating ylims when drawing AlignedReadTracks has been
  fixed.

- Numerous other little fixes that mainly aim at improving performance.

GWASTools
---------

Changes in version 1.3.16:

- Added convertVcfGds to extract bi-allelic SNPs from a VCF file.

- Added ncdfImputedDosage to convert output from common imputation
  programs to NetCDF.  assocTestRegression has an additional argument
  dosage=TRUE to be used with these files.

- Added vignette describing GWASTools data structures.

Changes in version 1.3.15:

- Bug fix in pedigreePairwiseRelatedness related to use of character
  identifiers.

Changes in version 1.3.14:

- assocTestRegression returns NA for snps where cases or controls are
  monomorphic, added assocTestFisherExact to use in that case.

- Added snp.exclude argument to pseudoautoIntensityPlot.

- Bug fix in messages reporting file read times when creating or
  checking netCDF files.

Changes in version 1.3.13:

- Added vignette on converting VCF to NetCDF with annotation.

- Prevent duplicateDiscordance from checking correlation by SNP in
  cases of no variation.

Changes in version 1.3.12:

- Added GdsReader and GdsGenotypeReader classes with dependency on
  gdsfmt.  GenotypeData objects can also be created with
  GdsGenotypeReader objects in the "data" slot.

Changes in version 1.3.11:

- Fixed bug in duplicateDiscordance when Y chromosome is not included.

Changes in version 1.3.10:

- Fixed bug in chromIntensityPlot so ideogram scales correctly if SNPs
  are excluded.

Changes in version 1.3.9:

- Fixed bug in assocTestCPH that could lead to false positives if
  additive model failed but GxE model did not.

- Allow multiple variables for stratified analysis in assocTestCPH.

Changes in version 1.3.8:

- Pedigree functions accept non-numeric identifiers and provide
  additional output.

Changes in version 1.3.7:

- In batchChisqTest, Yates correction cannot be bigger than the terms
  it corrects.  Changed to match bug fix to chisq.test in R 2.15.1.

Changes in version 1.3.6:

- Removed automatic subtitle from qqPlot.

- Allow selection of theoretical boundaries to draw in ibdPlot.

Changes in version 1.3.5:

- Added function asSnpMatrix to convert a GenotypeData object to a
  SnpMatrix object for use with snpStats.

Changes in version 1.3.4:

- Added chromosome ideograms to chromIntensityPlot and anomStatsPlot.
  anomStatsPlot has an option to put multiple anomalies on the same
  plot.

Changes in version 1.3.3:

- Updated vignette.

Changes in version 1.3.2:

- Use lazy loading of data.

- manhattanPlot and snpCorrelationPlot accept character vectors of
  chromosome; chrom.labels argument no longer used.

Changes in version 1.3.1:

- close method of NcdfReader returns invisibly.

HiTC
----

Changes in version 1.1.3:

NEW FEATURES

- Adding CITATION file

DEPRECATED AND DEFUNCT

- Update of getExpectedCount function to use the lowess() function
  (stats). The 'C' call in stats is now deprecated

Changes in version 1.1.2:

NEW FEATURES

- New package vignette

- Add new normalization method - normPerTrans

- Add importC and exportC function, to load and import csv file

SIGNIFICANT USER-VISIBLE CHANGES

- The CQC function now returns a matrix

- Simplify getExpectedCounts help page

- Update of import.my5C function. Simplify the import for matrix data

DEPRECATED AND DEFUNCT

- The export method is now replace by exportC. The standard csv format
  is now exported.

- The normPerZscore method is depracted. See normPerExpected instead

- Remove Bau et al. 5C dataset

BUG FIXES

- isBinned. Fix bug for interchromosomal interactions

- extracRegion. Add a chromosome parameter, and changes for
  interchromosomal data

- binningC. Changes for interchromosomal maps

- Sort xgi and ygi objects when the HTCexp constructor is called

- Force the xgi and ygi objects to have some ids

Changes in version 1.1.1:

NEW FEATURES

- Include the Nora et al (Nature 2012) 5C dataset (GSE35721).Two mouse
  samples are included in the package ; male undifferentiated ES cells
  (E14, GSM873935) and male embryonic fibroblasts (MEF, GSM873924).
  Only the cis interaction maps chrX vs chX are provided.

hpar
----

Changes in version 0.99.1:

- Added collate field <2012-09-14 Fri>

- Updated to HPA version 10 <2012-09-15 Sat>

- Updated installation part in section to use biocLite <2012-09-15 Sat>

Changes in version 0.99.0:

- Added vignette <2012-09-06 Thu>

Changes in version 0.1.0:

- Initial commit <2012-09-06 Thu>

htSeqTools
----------

Changes in version 1.3.1:

- Fixed bug in plotMeanCoverage

iPAC
----

Changes in version 0.3.0:

- First release of the iPAC package.

- ClusterFind method allows the use of MDS and linear mappers.

- Two beta methods available to reconcile data between the COSMIC and
  PDB databases.

IRanges
-------

Changes in version 1.16.0:

NEW FEATURES

- as( , "SimpleList"), as( , "CompressedList"), and as( , "List") now
  work on atomic vectors, and each element of the vector corresponds to
  an element of the returned List (this is consistent with as.list).

- Add as.list,Rle method.

- Add as.matrix,Views method. Each view corresponds to a row in the
  returned matrix. Rows corresponding to views shorter than the longest
  view are right-padded with NAs.

- Add FilterClosure closure class for functions placed into a
  FilterRules. Has methods for getting parameters and showing.

- Support 'na.rm' argument in "runsum", "runwtsum", "runq", and
  "runmean" methods for Rle and RleList objects.

- Add splitAsList() and splitAsListReturnedClass().

- Improve summary,FilterRules to support serial evaluation, discarded
  counts (instead of passed) and percentages.

- Make rename work on ordinary vector (in addition to Vector).

- Add coercion from RangedData to CompressedIRangesList, IRangesList,
  or RangesList. It propagates the data columns (aka values) of the
  RangedData object to the inner metadata columns of the RangesList
  object.

- Add 'NG' arg to PartitioningByEnd() and PartitioningByWidth()
  constructors.

- Make PartitioningByEnd() work on list-like objects (like
  PartitioningByWidth()).

- Fast disjoin() for moderate-sized CompressedIRangesList.

- Add countQueryHits() and countSubjectHits().

- coverage() now supports method="auto" and this is the new default.

- Add the flippedQuery(), levels(), ngap(), Lngap(), Rngap(),
  Lencoding(), and Rencoding() getters for OverlapEncodings objects.

- Add "encodeOverlaps" method for GRangesList objects.

- Enhance "[" methods for IRanges, XVector, XVectorList, and
  MaskCollection objects, as well as "[<-" method for IRanges objects,
  by supporting the following subscript types: NULL, Rle, numeric,
  logical, character, and factor. (All the methods listed above already
  supported some of those types but no method supported them all).

- Add remapHits() for remapping the query and subject hits of a Hits
  object.

- Add match,Hits method.

- Add %in%,Vector method.

- Add "compare", "==", "!=", "<=", ">=", "<", ">", "is.unsorted",
  "order", "rank", "match", and "duplicated" methods for XRawList
  objects. unique() and sort() also work on these objects via the
  "unique" and "sort" methods for Vector objects.

- Add expand() for expanding a DataFrame based on the contents of one
  or more designated columms.

- After being deprecated (in BioC 2.9) and defunct (in BioC 2.10), the
  "as.vector" method for AtomicList objects is back, but now it mimics
  what as.vector() does on an ordinary list i.e. it's equivalent to
  'as.vector(as.list(x), mode=mode)'. Also coercions from AtomicList to
  logical/integer/numeric/double/complex/character/raw are back and
  based on the "as.vector" method for AtomicList objects i.e. they work
  only on objects with top-level elements of length <= 1.

- DataFrame constructor now supports 'check.names' argument.

- Add revElements() generic with methods for List and CompressedList
  objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Splitting / relisting a Hits object now returns a HitsList instead of
  an ordinary list.

- Operations in the Ops group between a List and an atomic vector
  operand now coerce the atomic vector to List (SimpleList or
  CompressedList) before performing the operation. Also, operands are
  recycled and a better job is done returning zero length results of
  the correct type.

- Change the warning for 'Integer overflow ...' thrown by sum() on
  integer-Rle's

- DataFrame now coerces List/list value to DataFrame in [<-.

- Fix as.matrix,DataFrame for zero column DataFrames. Returns an
  nrow()x0 logical matrix.

- union,Hits method now sorts the returned hits first by query hit,
  then by subject hit.

- Add mcols() accessor as the preferred way (over elementMetadata() and
  values()) to access the metadata columns of a Vector object.

- By default, mcols(x) and elementMetadata(x) do NOT propagate the
  names of x as the row names of the returned DataTable anymore.
  However the user can still get the old behavior by doing mcols(x,
  use.names=TRUE).

- [<-,XVectorList now preserves the original names instead of
  propagating the names of the replacement value, which is consistent
  with how [<- operates on an ordinary vector/list.

- coverage() now returns a numeric-Rle when passed numeric weights.

- When called on a List object with use.names=TRUE, unlist() no longer
  tries to mimic the kind of non-sense name mangling that
  base::unlist() does (e.g. on list(a=1:3)) in a pointless effort to
  return a vector with unique names.

- Remove 'hits' argument from signature of encodeOverlaps() generic
  function.

- unique,Vector now drops the names for consistency with
  base::unique().

- Remove make.names() coercion in colnames<-,DataFrame for consistency
  with data.frame.

DEPRECATED AND DEFUNCT

- Deprecated tofactor().

- Remove RangesMatching, RangesMatchingList, and Binning classes.

- Change from deprecated to defunct: matchMatrix(), "dim" method for
  Hits objects, and RangesMatchingList().

BUG FIXES

- Fix bug in pintersect,IRanges,IRanges when input had empty ranges
  (broken since 2010-03-04).

- Avoid integer overflow in mean,Rle method by coercing integer-Rle to
  numeric-Rle internally.

- Change evaluation frame of with,List to parent.frame(), and get the
  enclosure correct in eval,List.

- Many fixes and improvements to coercion from RangesList to RangedData
  (see commit 68195 for the details).

- Fix "runValue" and "ranges" methods for CompressedRleList objects
  (broken for a very long time).

- shift,Ranges method now fails in case of integer overflow instead of
  returning an invalid Ranges object.

- mstack() now works on Vector objects with NULL metadata columns.

- In case of integer overflow, coverage() now puts NAs in the returned
  Rle and issues a warning.

- Fix bug in xvcopy,XRawList objects that prevented sequences from
  being removed from the cache of a BSgenome object. See commit 67171
  for the details.

- Fix issues related to duplicate column names in DataFrame (see commit
  67163 for the details).

- Fix a bunch of subsetting methods that were not subsetting the
  metadata columns: "[", "subseq", and "seqselect" methods for XVector
  objects, "seqselect" and "window" methods for XVectorList objects,
  and "[" method for MaskCollection objects.

- Fix empty replacement with [<-,Vector

- Make %in% robust on an empty 'table' argument when operating on Hits
  objects.


iSeq
----

Changes in version 1.9.0:

- Previously, peakreg fails when only one posterior probability is
  greater than ppcut or fdrcut, although this condition is really rare.
  This bug has been fixed.

- A tutorial of ChIP-seq data analysis using iSeq and an R script
  called iSeq.R that can be used as a command line program have been
  posted at
  https://sites.google.com/site/quincymobio/teaching-materials

- The document for iSeq has been updated.

isobar
------

 - Rockerbox import (just define the XXX.dat.peptides.csv as identifications)
 - possibility to define columns for XLS report in properties.R. e.g.
   xls.report.columns <- c("ratio","is.significant",
                           "ratio.minus.sd","ratio.plus.sd",
                           "p.value.ratio","p.value.sample",
                           "log2.ratio","log2.variance")


KEGGprofile
-----------

2012-09-19 KEGGprofile 0.99.2: Improvement in processing expression
matrix with only one time point

2012-06-29 KEGGprofile 0.99.1: Add the param lwd when type='lines'

limma
-----

Changes in version 3.14.0:

- limma license upgraded to GPL (>=2) instead of LGPL to match R
  itself.

- Many updates to the User's Guide.  Sections have been added on
  reading single channel Agilent and Illumina data. The chapter on
  experimental designs has been split into three chapters on
  single-channel, common reference and two-color designs respectively.
  The material on the fixed effect approach to technical replication
  has been deleted. There are new sections on nested interactions for
  factorial designs and on multi-level designs.

- The links to the Apoa1, Weaver and Bob1 datasets in the User's Guide
  have been updated to help users download the data themselves if they
  wish to repeat the case study analyses.

- The help page for camera() now cites published paper Wu and Smyth
  (NAR, 2012). In view of the results of this paper, the claim is no
  longer made on help page for geneSetTest() that genes might be
  treated as independent when the experimental units are genetically
  identical mice.

- Minor edits to CITATION file.

- New function propTrueNull() for fast estimation of the proportion of
  true null hypotheses from a vector of p-values.

- New function zscore() to compute z-score equivalents for deviates
  from any continuous distribution. Includes the functionality of the
  older functions zscoreGamma() and zscoreT() as special cases.

- roast() now accepts observation level weights, through a new argument
  'weights'.

- loessFit() now applies minimum and maximum bounds by default to avoid
  zero or infinite weights.  Equal weights are now treated as if the
  weights were NULL, even all zero weights, so that the lowess code is
  called instead of the loess code.

- When there are no weights, loessFit() now extracts residuals directly
  from the C code output instead of computing in R.

- fitFDist() now permits missing values for x or zero values for df1
  even when there is a covariate.  This means that squeezeVar() and
  eBayes() now work with trends even when not all the data values are
  informative.

- New argument 'file' for convest(), implementing edits contributed by
  Marcus Davy.  Arguments doplot and dereport renamed to 'plot' and
  'report'.

- Two improvements for plotMDS(). It now coerces labels to be
  character, and now makes extra room on the plot when the text labels
  are wide.

- plotMDS() no longer gives an error when the requested number of top
  genes is greater than the total number of rows of data.

- Code speed-up for alias2SymbolTable()

- any(duplicated()) replaced by anyDuplicated() in several functions.

- Fix to voom() so that it computes weights correctly even when the
  design matrix is not of full rank.

- Bug fix for roast() when the fitted model has only one coefficient.

lmdme
-----

Changes in version 0.99.1:

MINOR CHANGES

- Getters now accept a character vector in `term` parameter, in order
  to specify one than one term if required. In addition, `design` and
  `model` were added, and `pvalues` like `Fvalues` were changed to
  match slot names.  (Thanks to Valerie Obenchain)

- `lmdme` now works with NA presence in data matrix. This bug is
  related to lmFit intercept coefficient behavior, which breaks the
  data structure using drop (to numeric instead of keeping a matrix
  with one column) only if NA are present.

- `lmdme` is now the only constructor. Method `initialize` was erased
  due to different reasons as described in
  https://stat.ethz.ch/pipermail/bioc­devel/2012­August/003554.html

- `decomposition` example sections using `subset` parameter for
  simplicity (not all the data has to be decomposed in the example).

- Enhance of `biplot` and `screeplot` functions with `term` and `mfcol`
  to simplify the graphic output specification.

DOCUMENTATION

- `NEWS` file was added.

minfi
-----

Changes in version 1.3:

- Updated preprocessSwan to fix a bug when mSet was not set to the
  default value of NULL.  Specifically, now the "counts" tables is used
  to construct "subset".

- Changed the function manifestNew() to IlluminaMethylationManifest().

- Added IlluminaMethylationAnnotation().

- Added placeholders for unit testing based on RUnit.

- Introduced a new show method for MethylSet and RGChannelSet, derived
  from the eSet method in Biobase.

- The annotation slot of a MethylSet/RGChannelSet is now intended to
  _not_ be a scalar, but instead have length 2 with components 'array'
  and 'annotation'.  This foreshadows introdution of annotation
  packages for use with minfi.

- Reorganization of R files; rewriting of the man pages for MethylSet,
  RGChannelSet.

- getMeth, getUnmeth, getBeta, getM are now methods.

- bug fix to qcReport thanks to Tao Shi.

- Changes to getBeta / getM, both in terms of which arguments the
  methods take and how the values are computed.

- Changes to the manifest structure; it now has separate slots for
  genotype probes and these probes are no longer part of a MethylSet
  (using eg. preprocessRaw).  They can be accessed using
  getProbeInfo(rgSet, type) with type equal to "SnpI" or "SnpII".

- Introduction of mapToGenome, getLocations and the new class
  GenomicMethylSet.  man pages are reasonably complete, still need to
  add examples to the vignette.  This will be a standard part of an
  extended pipeline.

- Introduction of IlluminaHumanMethylation450lannotation.ilmn.v1.2
  which contains some new annotation needed for
  mapToGenome/getLocations.  This package will be split into several
  packages moving forward, in an attempt to harmonize efforts by us and
  Tim Triche.  getLocations/mapToGenome will stay the same.

- getControlTypes added (returns the different types of control
  probes).

- GenomicMethylSet now inherits a number of methods including
  granges(), start(), end() etc. from SummarizedExperiemnt.  They have
  therefore been deleted from minfi.

- Bugfix to getLocations(..., mergeManifest = TRUE).  It now longer
  throws an error.

- mapToGenome now returns a GenomicMethylSet ordered according to the
  chromosome name ordering chr1,..,chr22,chrX,chrY,unmapped, the last
  one not present if drop=TRUE (default).

MLInterfaces
------------

Changes in version 1.37.1:

- added NEWS file (lgatto) <2012-04-02 Mon>

- renamed 'predScores' to 'predScore' and new 'predScores' method that
  returns the full prediction score matrix (lgatto) <2012-04-02 Mon>

- remove ununsed t argument for predScore (lgatto) <2012-04-02 Mon>

- changed maxit to maxiter in RAB4es and decr to decreasing in fs.absT
  (partial argument matching note during checking) (lgatto) <2012-04-02
  Mon>

- fixed predScores <2012-04-03 Tue>

MLP
---

Changes in version 1.5.2:

- remove require(mouse4302mmentrezg.db) calls to keep R CMD check from
  giving inappropriate warnings

Changes in version 1.5.0:

- update exampleGeneSet.rda using GO.db for BioC 2.11

Changes in version 1.4.2:

- further improvements in getGeneSets and addGeneSetDescription

Changes in version 1.4.1:

- better checks on geneSetSource in getGeneSets

mosaics
-------

Changes in version 1.5.3:

BUG FIXES

- mosaicsRunAll(): Bug fix when byChr = TRUE.

Changes in version 1.5.2:

SIGNIFICANT USER-VISIBLE CHANGES

- constructBins(): Supports aligned read file formats for PET data
  (eland results and SAM formats).

- mosaicsRunAll(): Supports aligned read file formats for PET data.

- Add generateWig(): Constructs wiggle files for PET and SET data.

- Use tab separator instead of whitespaces for generateWig() and
  constructBins().

- Improve the vignette (case studies, example lines for input files,
  generateWig()).

BUG FIXES

- constructBins(): Bug fix for capping and excludeChr. Fix incorrect
  summary when byChr = TRUE.

- mosaicsRunAll(): Bug fix for excludeChr & handling the full path for
  chipFile and controlFile.

Changes in version 1.5.1:

SIGNIFICANT USER-VISIBLE CHANGES

- constructBins(): Chromosome information can now be specified.

- mosaicsRunAll(): Chromosome information can now be specified.

Changes in version 1.4.1:

BUG FIXES

- constructBins(): Bug fix for the "outfileLoc" argument.

- mosaicsFit(): Minor changes in two-signal-component model fitting.

- mosaicsPeak(): No warning with the updated IRanges package.

motifStack
----------

Changes in version 1.0.0:

NEW FEATURES

- draw DNA/RNA sequence logo

- draw Amino Acid sequence logo

- draw aligned motif logo stacks with phylogenetic tree

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

MSnbase
-------

Changes in version 1.5.25:

- fixed bug when quantifying exp of length 1 (reported by Colin
  Archer), added test <2012-09-26 Wed>

- fixed parallel default to FALSE <2012-09-26 Wed>

Changes in version 1.5.24:

- updated clean, removePeaks, combineFeatures, purityCorrect, trimMz,
  plot, MSnSet-class examples to not use readMSData <2012-09-25 Tue>

- fixed clean,MSnExp-method <2012-09-25 Tue>

- space in log message in extractPrecSpectra <2012-09-25 Tue>

Changes in version 1.5.23:

- updated quantify documentation <2012-09-24 Mon>

- changed all foo.class functions to foo_class <2012-09-24 Mon>

Changes in version 1.5.22:

- setting parallel default to FALSE <2012-09-22 Sat>

- enhanced parallel in DESCRIPTION and detectCores() passed to
  registerDoMC <2012-09-23 Sun>

Changes in version 1.5.21:

- removed registerDoMC no visible global function NOTE <2012-09-21 Fri>

- fixed NOTE about xvarname which was a bug <2012-09-21 Fri>

Changes in version 1.5.20:

- an immediate warning is thrown is any(centroided(object)) in
  quantify.MSnExp <2012-09-15 Sat>

- mzTab file and loading time is now recorded in processingData
  <2012-09-15 Sat>

- temporarily dontrun'ing' the mzTab read/write examples as OLS is down
  (and breaks rols) (note: modifed Rd files, not roxygen template)
  <2012-09-15 Sat>

Changes in version 1.5.19:

- updated all ReporterIons rda data <2012-09-13 Thu>

- Using && in testing parallel, require(foreach and doMC) <2012-09-14
  Fri>

Changes in version 1.5.18:

- new log method for MSnSet instances <2012-09-13 Thu>

- new MAplot methods using generic and ma.plot/mva.pairs from affy
  <2012-09-13 Thu>

Changes in version 1.5.17:

- changed TMT7[7] mass from 229.26 to 230.17 and ReporterIons
  descriptions <2012-09-12 Wed>

Changes in version 1.5.16:

- parallel quantify is now always set to FALSE on Windows, fixing
  example checking issues <2012-09-11 Tue>

- fixed types in plotMzDelta man <2012-09-11 Tue>

Changes in version 1.5.15:

- updating code to ggplot2 v0.9.2 <2012-09-07 Fri>

Changes in version 1.5.14:

- added platform test to use doMC (d.tenenbaum) <2012-08-31 Fri>

Changes in version 1.5.13:

- updated fillUp function <2012-08-15 Wed>

- added tikzDevice to suggests <2012-08-16 Thu>

- tikzDevice no longer on CRAN - removing from Suggests and using pdf
  as device in vignette <2012-08-16 Thu>

Changes in version 1.5.12:

- using knitr instead of pgfSweave and misc vignette updated<2012-08-13
  Mon>

- spectrum2 reporter plotting params updates <2012-08-14 Tue>

- added reporterNames to NAMESPACE <2012-08-14 Tue>

Changes in version 1.5.11:

- type in filterNA log messaging and also rounding pNA <2012-06-07 Thu>

- typo in demo vignette - Gb instead of Mb <2012-07-15 Sun>

Changes in version 1.5.10:

- Fixed readMSData instrumentInfo handling (reported by Gopuraja
  Dharmalingam) <2012-06-05 Tue>

- new multiple file loading test in test_io <2012-06-05 Tue>

Changes in version 1.5.9:

- topN now properly updates processingData <2012-06-01 Fri>

- combineFeatures updates featureNames based on the groupBy argument -
  updated demo vignette and man accordinlgy <2012-06-01 Fri>

- additional parameters were not passed when normalise using vsn
  <2012-06-01 Fri>

Changes in version 1.5.8:

- added aa data in environment data.frame <2012-05-22 Tue>

- fixed MSnbase:::subsetBy (used by topN) when ncol(object) == 1
  <2012-05-25 Fri>

- new nQuants function <2012-05-31 Thu>

- minor vignettes updates <2012-05-31 Thu>

- nQuants now return a matrix with col names, taken form
  sampleNames(object) <2012-05-31 Thu>

Changes in version 1.5.7:

- changed title method to exptitle to avoid conflict/confusion with
  graphics::title and consitency with expinfo method <2012-05-15 Tue>

- changed email to expemail <2012-05-15 Tue>

- Added rols to Suggests <2012-05-15 Tue>

Changes in version 1.5.6:

- new ionSource, analyser, detectorType, title accessor methods for
  MIAPE, pSet and MSnSet classes <2012-05-06 Sun>

- updated quantify example to use data(itraqdata) instead of reading
  dummyiTRAQ.mzXML <2012-05-09 Wed>

- initial mzTab write support <2012-05-10 Thu>

- mzTab read support <2012-05-10 Thu>

- updated demo and io vignettes with mzTab info <2012-05-10 Thu>

Changes in version 1.5.5:

- added email slot in MIAPE and accessor <2012-05-04 Fri>

- added expinfo methods to pSet and MSnSet <2012-05-04 Fri>

- updated itraqdata <2012-05-04 Fri>

- misc man typos fixed <2012-05-04 Fri>

- quantify now properly propagates processingData <2012-05-04 Fri>

Changes in version 1.5.4:

- automatically populating experiment data instrument info while
  reading data <2012-04-30 Mon>

- msInfo fixed and exported <2012-04-30 Mon>

- update demo vignette with 14 fractions analysis paragraph <2012-05-01
  Tue>

Changes in version 1.5.3:

- extractSpectra is now defunct <2012-04-20 Fri>

- caching full header in level 1; this is required when and MSnExp
  instance with *many* spectra (created from many raw files) is
  quantified - calling header(object) is a too big overhead compared to
  actual reporter quantification. <2012-04-20 Fri>

- The header() method now uses the cached dataframe if level >= 1; the
  (unexported) .header function can be used to generate the dataframe
  using the assayData slot data. <2012-04-20 Fri>

- Setting processingData in MSnSet initialisation. <2012-04-20 Fri>

- Dropping index column from header. <2012-04-20 Fri>

- new Spectrum class v0.2.0 has tic slot. <2012-04-21 Sat>

- *tic* method (data stored as a Spectrum slot) now returns _total ion
  current_ (as commonly used) and _total ion count_ is obtain using
  *ionCount*. <2012-04-21 Sat>

- fixed normalisation boxplot titles and other tic/ionCount changes in
  demo vignette. <2012-04-21 Sat>

- removed qual subetting in MSnSet's "[" method <2012-04-24 Tue>

Changes in version 1.5.2:

- transposing and MSnSet does _silently_ drop the protocolData now
  <2012-04-03 Tue>

- fixed MSnExp pData creation for multiple files, feature names have a
  .fileNumber extension now. <2012-04-19 Thu>

- testing for uniqueness of files (filenames) in readMSData <2012-04-19
  Thu>

- updated itraqdata.RData <2012-04-19 Thu>

- added a paralle argument to quantify and using paralle = FALSE in
  vignette, to avoid duplicated display of the command <2012-04-19 Thu>

- defined "reporterNames<-" generics <2012-04-19 Thu>

- fixed warning in readMSData where all not used for comparison of
  verctors of length > 1 <2012-04-20 Fri>

Changes in version 1.5.1:

- updated NA warning message in readIspyData <2012-03-05 Mon>

- fixed combineFeatures/combineMatrixFeatures for 1 sample <2012-03-20
  Tue>

- image method for MSnSet instances <2012-03-20 Tue>

- fixed bug in plotNA (first t was wring) <2012-03-21 Wed>

- import plot from stats4 <2012-03-21 Wed>

- using reshape2 <2012-03-30 Fri>

Changes in version 1.5.0:

- new devel version bump

ncdfFlow
--------

1.optimize [[ accessor by merging two routines into [[ method
and\nremoving some unnecessary check or subsetting

netresponse
-----------

Changes in version 1.9.17 (2012-06-28):

- added remove.negative.edges function

- added positive.edges argument to detect.responses

Changes in version 1.9.15 (2012-05-23):

- igraph dependency moved to igraph0

Changes in version 1.9.14 (2012-05-14):

- dmt added to dependencies

Changes in version 1.9.13 (2012-05-04):

USER-VISIBLE CHANGES

- added PlotMixtureMultivariate function

Changes in version 1.9.12 (2012-05-03):

USER-VISIBLE CHANGES

- added parallelization on update.model.pair step

BUG FIXES

- removed merge updates for potential subnet merges that would exceed
  max size

Changes in version 1.9.10 (2012-04-30):

SIGNIFICANT USER-VISIBLE CHANGES

- added bic.threshold option in detect responses and downstream
  functions

- latent.class.analysis function replaced with mixture.model

BUG FIXES

- Rd conflicts resolved

Changes in version 1.9.07 (2012-04-25):

SIGNIFICANT USER-VISIBLE CHANGES

- changed all documentation into Roxygen

- added bic.mixture.univariate function

- added argument mixture.method in detect.responses and to consecutive
  downstream functions

BUG FIXES

- added information.criterion option to edge.delta costind.ab
  calculation

- "c.max <- max.responses - 1"--> "c.max <- max.responses" in vdp.mixt

Changes in version 1.9.05 (2012-04-24):

- added the plot.mode = "pca" in plot.responses function

- added tests/mclust-mixture.R example

Changes in version 1.9.04 (2012-04-21):

- removed the redundant network.nodes parameter from independent.models
  input

oligo
-----

Changes in version 1.22:

USER VISIBLE CHANGES

- Improved results by getProbeInfo();

- fitPLM, coefs and resids are now Defunct. Use fitProbeLevelModel,
  coef and residuals respectively. 'coef' and 'residuals' follow the
  standards used elsewhere in R;

BUG FIXES

- Fixed problem caused by the fact that oligoClasses had its own
  annotation() method when BiocGenerics added a new one;

- Several fixes to probe selector, allowing 'target' to be used;

- PA Calls didn't know about target;

phyloseq
--------

* plot_ordination() powerful, flexible ordination plots built with ggplot2
* plot_heatmap() easy, flexible heat maps built with ggplot2
* plot_tree() easy, flexible annotated tree plots built with ggplot2
* make_sample_network(), plot_sample_network() - microbiome network
  visualization
* plot_richness() for easy, flexible summary of species richness
* Parallel Fast UniFrac
* distance() wrapper for ecological distance calculations
* ordinate() wrapper, calculates many different ordination methods.
* General importer for all supported data formats: import()
* BIOM format import: import_biom() function
* Support for Double Principle Coordinate Analysis (DPCoA)
* Several published exampled datasets included
* Bioconductor development release updates.
* Lots of documentation updates.
* Lots and lots of fixes and improvements.

phyloseq 1.1.50

- Fixed several compatibility issues to support latest version of
  ggplot2 (0.9.2).
- Also changes plot_richness_estimates() to plot_richness().

phyloseq 1.1.45

- Backward compatibility for import_qiime_sampleData, now superseded
  by import_qiime_sample_data
- Added a functioning example based on the GlobalPatterns example
  sample-map file included in the package extdata.

phyloseq 1.1.44

- Fixed minor bug in tax_glom function.  Thanks to Katie Shelef for
  the bug report. Bug only affected tax_glom behavior when the
  right-most rank was specified as the position for merging.

phyloseq 1.1.43

- fixed distance() issue from species/taxa replacement for type
  argument.

phyloseq 1.1.42

- fixed make_network/plot_network issue from species/taxa replacement
  for type argument.

phyloseq 1.1.41

- Fixed documentation for `prune_taxa` and `prune_samples`
- Updated `prune_samples` method to allow for logical vectors.
- Fixed `prune_taxa` so that it properly fails with a message if the
  taxa argument is a logical of wrong length.  There was some
  potential (and no warning) for unpredictable vector-recycling with
  short vectors in the old implementation.

phyloseq 1.1.40

## Huge Update and Renaming Event.
- Made all functions use an *underscore* for English word delimiter,
  if they were using an abbreviation.
- Replaced "species" in all function names with "taxa". 
- These changes are all backward compatible, for now, so your old code
  should work. Let me know if it doesn't and I will quickly make the
  adjustment. This will remain true through the next official release,
  but functional references to "species" will not be supported
  afterward, except in the occasions where you actually mean taxonomic
  species, like `tax_glom(x, "species")`.

phyloseq 1.1.33

* Revise taxglom() such that it handles phyloseq and taxonomyTable
  classes, throws warning otherwise. It should not take a
  manually-produced character vector, as this is roughly equivalent to
  functionality supported in other method, especially
  prune_species()/merge_species().
* Also added unit-tests and executable examples for taxglom().
* Got rid of taxglom.internal, incorporated directly into taxglom().
  taxglom() is no longer an S4-method, and doesn't need to be now that
  the character-vector argument option is omitted, with S4-class
  handling delegated to merge_species().
* Updated "taxTab<-" to be S4 assignment, clearer handling of taxonomy
  Table assignments, especially useful for taxglom.

phyloseq 1.1.29

* Add unit tests and example files for import_biom (as well as
  import("biom",...) ).

phyloseq 1.1.28

* Added rarefy_even_depth() function for random subsampling of
  microbiome samples to the same number of reads. Default uses the
  minimum total reads among the samples in the dataset. This is based
  on the core "sample" function, which can have its random number
  generator fixed by set.seed for reproducibility.

phyloseq 1.1.27

* Fix bug in plot_ordination that caused an error rather than produce
  unannotated plots when sampleData absent in the input.

phyloseq 1.1.23-26

* Added unit tests and bugfixes

phyloseq 1.1.19-22

* Improving import_qiime() importer to handle large datasets, like the
  HMPv35 dataset, for example, while also providing useful status
  messages during non-trivial imports that might take 10 minutes or
  more to complete.

phyloseq 1.1.18

* Added replicate labels as a "Sample" factor in the soilrep dataset.

phyloseq 1.1.17

* Fix possible bug that results from the latest version (0.6+) of
  igraph not being backward compatible. A stable igraph0 package is
  available on CRAN as a stop-gap, and so all igraph dependencies were
  migrated to "igraph0" until the phyloseq-source can be updated to
  match the igraph latest.

phyloseq 1.1.15

* plot_heatmap: Added default (but adjustable) threshold to omit
  taxa/sample labels

phyloseq 1.1.14

* Update import_qiime() function to import latest non-BIOM qiime
  output files.  Also added check for presence of taxonomy information
  (consensus lineage).

phyloseq 1.1.10

* Add plot_heatmap() function, for easy flexible heat maps built with
  ggplot2

phyloseq 1.1.8-9

* Fix bug for some variants of new BIOM format
* Add import_RDP_otu() import function for new RDP pipeline export
  format

phyloseq 1.1.7

* Removed the old plot_tree_phyloseq() function, in favor of the new
  ggplot2-based plot_tree()
* Uncommented / tested formal examples in documentation of
  plot-functions
* Updated variable names and doc for the plot_taxa_bar() function

phyloseq 1.1.6

* Update vignette with plot_tree() example, replacing the old
  base-graphics function, plot_tree_phyloseq().
* Fix bug in legend for trees with size mapped to abundance

phyloseq 1.1.5

* Add initial version of tree_plot(), built with ggplot2
* Adds several internal functions borrowed from devel version of ggphylo

phyloseq 1.1.4

* Add errorIfNULL option to auxiliary accessors
  (e.g. sample.variables(), rank.names())

phyloseq 1.1.1-3

* R version updated to match Bioconductor, R-2.15.0+
* ape-package version updated to 3.0+
* ape-package now import dependency
* ggplot2-package version updated to 0.9.0+
* ggplot2-package now import dependency

PWMEnrich
---------

Changes in version 1.0.2:

- Use core package parallel instead of doMC

Changes in version 1.0.1:

- Fix a typo in test cases and remove doMC as build dependency

Changes in version 1.0.0:

- Initial release of the package

qpgraph
-------

Changes in version 1.14:

NEW FEATURES

- new qpPlotMap() function to show associations between genetic markers
  and gene expression profiles using their positions along the genome

- conditional independence tests and non-rejection rate estimation with
  missing data via complete-case analysis and the EM algorithm

- new qpAllCItests() function to perform multiple conditional
  independence tests with a fixed conditioning set.

- new arguments fix.Q and restrict.Q for functions estimating the
  non-rejection rate in order to restrict and fix variables in the
  conditioning subsets of the independence tests.

R453Plus1Toolbox
----------------

Changes in version 1.7.4:

- Added ability to read in data from Roche's latest Amplicon Variant
  Analyzer (v2.7).

- New method ava2vcf to write the variants in an AVA-Set to a file in
  VCF format.

- Multiple bugfixes.

ReactomePA
----------

Changes in version 1.1.1:

- remove geneID2Name, instead import EXTID2NAME from
  DOSE. <2012-03-12, Mon>

- remove most of the codes in enrichPathway, instead import
  enrich.internal in DOSE.  implement some S3 method for mapping.
  pathID2Name was rename to TERM2NAME, which will called by
  enrich.internal. <2012-03-12, Tue>

RedeR
-----

 Version    Date Category Text
 2011-03-23 <NA> <NA>

ReQON
-----

Changes in version 1.3.4:

NEW FEATURES

- An additional option has been added to ReadPosErrorPlot.R.  The
  option "startpos" now allows users to designate the starting read
  position to be plotted.  The default start position is 1.

BUG FIXES

- The required version of R has increased to 2.15.  Major changes to
  the ReQON package were made in version 1.2.0, and older versions of R
  would download ReQON v. 1.0.0, which is incompatible with the current
  documentation.

- There were minor inconsistencies between the FWSE reported in the
  output plots and FWSE calculated from the recalibrated BAM file,
  which has now been fixed.

Changes in version 1.3.2:

NEW FEATURES

- ReQON no longer recalibrates 'N' bases.  It returns the original
  quality score for these bases in the output BAM file.

- Minor modifications have been made to FWSEplot.R.  Points are now
  shaded according to the relative frequency of bases assigned that
  quality score.  See the reference manual for more details.

Rgraphviz
---------

Changes in version 2.1:

- Rgraphviz now requires Graphviz >= 2.16.

- Added graphvizCapabilities() that reports the capabilities of
  Graphviz.  This requires Graphviz >= 2.28 and returns NULL if the
  Rgraphviz installation does not support it.

- Rgraphviz now comes bundled with Graphviz 2.28, greatly simplifying
  installation.

- Numerous bugfixes for bugs introduced in the 2.x.x versions.

Risa
----

This is the first release version of the package. It contains
functionality to parse ISAtab datasets into an R object from the
ISAtab class. It also provides functionality to save the ISA-tab
dataset, or each of its individual files. Additionally, it is also
possible to update assay files. Currently, metadata associated to
proteomics and metabolomics-based assays (i.e. mass spectrometry) can
be processed into an xcmsSet object (from the xcms R package

rols
----

Changes in version 0.99.12:

- updated README.md (for github) and .Rbuildignore <2012-09-19 Wed>

Changes in version 0.99.11:

- updated biocViews <2012-09-18 Tue>

Changes in version 0.99.10:

- using knitr instead of pgfSweave <2012-08-13 Mon>

Changes in version 0.99.9:

- fixed olsQuery issued warning when missing(ontologyNames)
  irrespective of extact <2012-05-22 Tue>

Changes in version 0.99.8:

- really using pdfSweave <2012-05-18 Fri>

Changes in version 0.99.7:

- using pdfSweave <2012-05-16 Wed>

- changed VignetteIndexEntry <2012-05-16 Wed>

- removed roxygen2 from Suggests <2012-05-16 Wed>

- added README file describinh Rd generation with roxygen2 <2012-05-16
  Wed>

Changes in version 0.99.6:

- implementing Mark Carlson's suggestions

- (Temporarily) removing pgf <2012-05-15 Tue>

- added Collate field <2012-05-15 Tue>

- moving generics and class definitions to AllGenerics.R and
  AllClasses.R <2012-05-15 Tue>

- olsQuery now repeats query 'n' times if reply is empty - see man page
  for a brief discussion <2012-05-15 Tue>

- updated vignette to discuss off/on-line data access <2012-05-15 Tue>

Changes in version 0.99.5:

- downgraded SSOAP dependency to (>= 0.8.0), as later SSOAP versions
  are not available through biocLite (due to check errors). rols works
  with SSOAP 0.8.0 and XMLSchema 0.7.2. <2012-05-11 Fri>

Changes in version 0.99.4:

- added url to DESCRIPTION <2012-05-11 Fri>

- released CVParam validity contrains: now CVParam name must match term
  or any synonym <2012-05-11 Fri>

- new CVParam constructor <2012-05-11 Fri>

Changes in version 0.99.3:

- minor vignette update <2012-05-04 Fri>

- new CVParam class <2012-05-08 Tue>

- More verbose CVParam validity error msg <2012-05-08 Tue>

- new 'rep' method for CVParam <2012-05-09 Wed>

Changes in version 0.99.2:

- olsQuery has a new 'exact' parameter <2012-05-04 Fri>

Changes in version 0.99.1:

- import(XMLSchema) to remove warning at startup <2012-05-04 Fri>

Changes in version 0.99.0:

- vignette update <2012-04-30 Mon>

- added GO.db to Suggests for vignette <2012-04-30 Mon>

- bumped version to 0.99 for Bioc submission <2012-04-30 Mon>

Changes in version 0.2.2:

- added a vignette <2011-12-18 Sun>

- added xtable in suggests for vignette <2011-12-18 Sun>

- added parents() and childrenRelations() query functions <2011-12-18
  Sun>

Changes in version 0.2.1:

- added specific SSOAP and XMLSchema version requirements <2011-12-18
  Sun>

- fixed a couple of typos <2011-12-18 Sun>

Changes in version 0.2.0:

- initial release

RPA
---

Changes in version 1.13.06 (2012-09-09):

- affinity estimates storing option added to rpa.online

- moved sigma2 -> tau2 to make notation compatible with publications

Changes in version 1.13.02 (2012-05-27):

NEW FEATURES

- added arguments in rpa.online function: save.batches.dir,
  keep.batch.files, unique.run.identifier

- changed the default for save.batches option in rpa.online into TRUE

BUG FIXES

- fixed load.batches bug in summarize.batches


Rsamtools
---------

Changes in version 1.10.0:

NEW FEATURES

- BamFile and TabixFile accept argument yieldSize; repeated calls to
  scanBam and scanTabix return successive yieldSize chunks of the file.
  readBamGappedAlignments, VariantAnnotation::readVcf automatically
  gain support for yield'ing through files.

- Add getDumpedAlignments(), countDumpedAlignments(), and
  flushDumpedAlignments() low-level utilities for manipulating
  alignments dumped by findMateAlignment().

- Add quickBamCounts() utility for classifying the records in a BAM
  file according to a set of predefined groups (based on the flag bits)
  and for counting the nb of records in each group.

SIGNIFICANT USER-VISIBLE CHANGES

- scanBamFlag isValidVendorRead deprecated in favor of
  isNotPassingQualityControls

- Rename makeGappedAlignmentPairs() arg 'keep.colnames' -> 'use.mcols'.

BUG FIXES

- close razip, bgzip files when done

- bamReverseComplement<- failed to return the updated object

- scanBcfHeader works on remote files

- allow asBam to work without warnings on header-only SAM files

- some bug fixes and and small performance improvements to
  findMateAlignment()

- fix bug in readBamGappedAlignmentPairs() where fields and tags
  specified by the user were not propagated to the returned
  GappedAlignmentPairs object

Rsubread
--------

Changes in version 1.8.0:

NEW FEATURES

- Fine-tuning of read alignment function (align) and exon-exon junction
  detection function (subjunc).

- Major update to SNP calling algorithm (now use Fisher' exact tests to
  call SNPs).

- More efficient implementation for removing duplicated reads (eg.
  supporting multi-threads and using less memory).

ShortRead
---------

SIGNIFICANT USER-VISIBLE CHANGES

- as(ShortReadQ, "matrix") now accepts ShortReadQ instances with
  heterogenous widths, returning a matrix x[i, j] with NA values in
  when j > width()[i].

NEW FEATURES

- FastqStreamer accepts IRanges for selecting input records

BUG FIXES

- readAligned, type="BAM" correctly adds requried 'what' elements

- FastqSampler would only randomize first read; introduced 1.13.9
  2011-12-02, fixed 1.15.4 2012-04-25

- report(qa, ...) no longer produces obviously confused base calls per
  cycle

- FastqFileList would fail to initialize correctly from a character
  vector

synapter
--------

Changes in version 0.99.15:

- updated ref in package man (2) <2012-09-14 Fri>

- updated Synapter show method to display short log <2012-09-14 Fri>

Changes in version 0.99.14:

- updated ref in package man <2012-09-11 Tue>

Changes in version 0.99.13:

- bumping version number to take all previous changes into account
  <2012-08-29 Wed>

Changes in version 0.99.12:

- fixed build problem on windows <2012-08-27 Mon>

- updated refs in vignette <2012-08-28 Tue>

Changes in version 0.99.11:

- using dev='pdf' in vignette and removing tikzDevice (not on CRAN
  anymore) from Suggests<2012-08-16 Thu>

Changes in version 0.99.10:

- added tikzDevice suggestion for vignette <2012-08-15 Wed>

Changes in version 0.99.9:

- using knitr instead of pgfSweave <2012-08-13 Mon>

- new verbose arg to loadIdentOnly taken into account by
  estimateMasterFdr and makeMaster <2012-08-13 Mon>

Changes in version 0.99.8:

- vignette update with summary table <2012-07-17 Tue>

- vignette update with additional PLGS slides <2012-07-17 Tue>

Changes in version 0.99.7:

- updated vignette <2012-07-12 Thu>

- removed README.org file <2012-07-12 Thu>

- added github page in DESCRIPTION <2012-07-12 Thu>

Changes in version 0.99.6:

- Masterpeptides has new fdr and method slots <2012-07-10 Tue>

- added method arg to makeMaster <2012-07-10 Tue>

- fixed non-match precuror.leIDs when using master <2012-07-11 Wed>

Changes in version 0.99.5:

- removed 'precursor.Mobility' columns from light merged and matched
  csv output, as these are not available when a MSe master is used
  <2012-07-09 Mon>

Changes in version 0.99.4:

- use log2 for t.test <2012-07-05 Thu>

- updates to vignette <2012-07-06 Fri>

Changes in version 0.99.3:

- git/svn integration messing up and testing <2012-07-03 Tue>

- late night testing <2012-07-03 Tue>

Changes in version 0.99.2:

- changed merged to total FDR in MasterFdrResults plot x label, as in
  manuscript <2012-06-13 Wed>

- mention multtest in vignette <2012-06-13 Wed>

- qvalue refs <2012-06-14 Thu>

- Synapter.Rd updates, to match manuscript nomenclature <2012-06-14
  Thu>

- Added 'Getting help' section in vignette and package man page
  <2012-06-14 Thu>

- Changed HDMSe to ident in plotGrid plot legend <2012-06-25 Mon>

Changes in version 0.99.1:

- Latex typo in vignette <2012-06-11 Mon>

- Dan's suggestions for <2012-06-13 Wed>

Changes in version 0.99.0:

- bump version for Bioconductor submission <2012-06-08 Fri>

TargetSearch
------------

Changes in version 1.14.0:

NEW FEATURES

- New function 'fixRI'. This function can be used to correct RI markers
  or to manually force their location to specific retention times if,
  for example, the RI markers were not co-injected with the biological
  samples. Replaces the now deprecated 'fixRIcorrection'.

- New function 'riMatrix'. This function searches RI markers in RI
  files instead of CDF files.

- Improvements in CDF import functions: - Automatic detection of m/z
  range. - Detection and correction of CDF files with non-integer m/z
  values.  These values are converted to nominal mass.

- ImportLibrary: - New parameter 'file.opt'. A list containing
  arguments to be passed to 'read.delim'. - It can take a data frame
  instead of a file to create a library object.

SIGNIFICANT USER-VISIBLE CHANGES

- Function 'fixRIcorrection' is deprecated. Use 'fixRI'.

BUG FIXES

- The .Call function in FindPeaks.R would incorrectly coerce RI limits
  to integers instead of double.

- ImportLibrary: - Fixed bug when reading one-metabolite libraries. -
  Check for unexpected quotation mark characters in input file or input
  data.frame. They will be removed.

- Profile: check for at least three top masses to calculate spectra
  similarity scores.

- quantMatrix: use character indeces instead of numeric indeces.

tigre
-----

Changes in version 1.12.0:

- Minor speed optimisation

- Minor documentation reorganisation and updates

TransView
---------

Changes in version 0.99.3:

NEW FEATURES

- Renamed all C files. Registering now takes place in a separate file.

- plotTV can now alternatively take a character vector of IDs matching
  'transcript_id' in the gtf instead of a GRanges object. This is
  useful to plot RNA-Seq data matching without regions.

BUG FIXES

- Added PACKAGE argument to .call

- readthrough_pairs argument to parseReads declared experimental.

- Removed the matching RNA-Seq data set and added pasillaBamSubset as
  an example for RNA-Seq visualisation.

Changes in version 0.99.2:

NEW FEATURES

- Added getter methods for all slots of DensityContainer and inherited
  classes for individual access. Added setter methods for spliced and
  ex_name

- removed dc.size() and added the slot 'size' instead. Also added env
  and data_pointer slots to class DensityContainer

- improved documentation of S4 methods with usage section

- gtf2df and macs2df renamed to gtf2gr and macs2gr respectively. They
  are now returning a GRanges object.

- annotatePeaks takes and returns only GRanges objects.

- id2tss renamed to peak2tss and now returning an updated GRanges
  object

- plotTV accepts GRanges objects rather than data.frames for gtf and
  peaks input

BUG FIXES

- Removed all direct slot accessions

- Removed unnecessary call to dyn.load and dyn.unload in parseReads()

Changes in version 0.99.1:

NEW FEATURES

- plotTV returns ordering to reproduce kmeans clustering

BUG FIXES

- sliceNT did not flip the exons on the negative strand

- Argument rowv to plotTV did not work with expression threshold
  remove_lowex

- Installation on 32bit Windows was not possible

xcms
----

CHANGES IN VERSION 1.33.16

USER VISIBLE CHANGES

    o diffreport and plotEIC have a new parameter mzdec, with is the
      number of decimal places of the m/z values in the EIC plot title

CHANGES IN VERSION 1.33.16

UPDATED FEATURES:

    o Lock mass gap filler now works with netCDF lock mass function
      file to find the exact times of the scans and works with the
      newer Waters MS instruments.

CHANGES IN VERSION 1.33.15

BUG FIXES
  
    o scanrage is now honoured in xcmsSet, also when in parallel mode

CHANGES IN VERSION 1.33.14

BUG FIXES
  
    o scanrage is now honoured in xcmsRaw, and consequently 
      also in xcmsSet(matchedFilter), where previously 
      it was ignored.

CHANGES IN VERSION 1.33.13

BUG FIXES
  
    o write.cdf() has been fixed to write files AMDIS can read

CHANGES IN VERSION 1.33.12

BUG FIXES
  
    o write.mzData adds Polarity to the file if available

CHANGES IN VERSION 1.33.11

USER VISIBLE CHANGES
  
    o centWave uses a new method to estimate local noise which
    improves detection of closely spaced peaks

NEW FEATURES

    o placeholder

BUG FIXES

    o group.mzClust was failing when result had one peak

For more details and all changes before May 2012 please see the (now
discontinued) CHANGELOG in the source package (inst/ folder).

CHANGED BEHAVIOUR since Version 1.32:

Other Changes since Version 1.32:

* improved mzData writing, now includes MSn spectra and less verbose.
* improved netCDF writing, but not yet good enough for AMDIS

yaqcaffy
--------

Changes in version 1.17.2:

- fixed NOTE about no visible binding for global variable 'latex'
  <2012-09-21 Fri>

- 'Warning: found methods to import for function 'as.list' but not the
  generic itself' comes from upstreams. <2012-09-22 Sat>

Changes in version 1.17.1:

- fixed boxplot alignment <2012-05-30 Wed>

- deleted vignette Makefile <2012-05-30 Wed>


Packages removed from the release
=================================

No packages have been removed from the release.
