Bioconductors:

We are pleased to announce Bioconductor 2.13, consisting of 749
software packages, 179 experiment data packages, and more than 690
up-to-date annotation packages. 

There are 84 new software packages, and many updates and improvements
to existing packages; Bioconductor 2.13 is compatible with R 3.0.2,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release includes an updated
[Bioconductor Amazon Machine Image](/help/bioconductor-cloud-ami/).

Visit [http://bioconductor.org](/)
for details and downloads.

Contents
--------

* Getting Started with Bioconductor 2.13
* New Software Packages
* NEWS from new and existing packages
* Packages removed from the release

Getting Started with Bioconductor 2.13
======================================

To update to or install Bioconductor 2.13:

1. Install R 3.0.2.  Bioconductor 2.13 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](/install/)

New Software Packages
=====================

There are 84 new packages in this release of Bioconductor.

AllelicImbalance: Provides a framework for allelic specific
expression investigation using RNA-seq data


ampliQueso: The package provides tools and reports for the analysis
of amplicon sequencing panels, such as AmpliSeq


ArrayTV: Wave correction for genotyping and copy number arrays


ASSET: An R package for subset-based analysis of heterogeneous traits
and subtypes


BADER: For RNA sequencing count data, BADER fits a Bayesian
hierarchical model. The algorithm returns the posterior
probability of differential expression for each gene between two
groups A and B. The joint posterior distribution of the variables
in the model can be returned in the form of posterior samples,
which can be used for further down-stream analyses such as gene
set enrichment.


BAGS: R package providing functions to perform geneset significance
analysis over simple cross-sectional data between 2 and 5
phenotypes of interest.


BiGGR: This package provides an interface to simulate metabolic
reconstruction from the BiGG database(http://bigg.ucsd.edu/) and
other metabolic reconstruction databases. The package aids in
performing flux balance analysis (FBA). Metabolic networks and
estimated fluxes can be visualized using hypergraphs.


bioassayR: bioassayR provides tools for statistical analysis of small
molecule bioactivity data


BiocParallel: This package provides modified versions and novel
implementation of functions for parallel evaluation, tailored to
use with Bioconductor objects.


BiocStyle: Provides standard formatting styles for Bioconductor
documents. The vignette illustrates use and functionality.


BiRewire: Fast functions for bipartite network rewiring through N
consecutive switching steps (See References) and for the
computation of the minimal number of switching steps to be
performed in order to maximise the dissimilarity with respect to
the original network. Includes function for the analysis of the
introduced randomness across the switching and several other
routines to analyse the resulting networks and their natural
projections. Extension to undirected networks (not bipartite) is
also provided.


CexoR: Strand specific peak-pair calling in ChIP-exo replicates. The
cumulative Skellam distribution function (package 'skellam') is
used to detect significant normalized count differences of
opposed sign at each DNA strand (peak-pairs). Irreproducible
discovery rate for overlapping peak-pairs across biological
replicates is estimated using the package 'idr'.


ChAMP: The package includes quality control metrics, a selection of
normalization methods and novel methods to identify
differentially methylated regions and to highlight copy number
aberrations.


ChemmineOB: ChemmineOB provides an R interface to a subset of
cheminformatics functionalities implemented by the OpelBabel C++
project. OpenBabel is an open source cheminformatics toolbox that
includes utilities for structure format interconversions,
descriptor calculations, compound similarity searching and mor.
ChemineOB aims to make a subset of these utilities available from
within R. For non-developers, ChemineOB is primarily intended to
be used from ChemmineR as an add-on package rather than used
directly.


chipenrich: ChIP-Enrich performs gene set enrichment testing using
peaks called from a ChIP-seq experiment. The method empirically
corrects for confounding factors such as the length of genes, and
the mappability of the sequence surrounding genes.


cleanUpdTSeq: This package uses the Naive Bayes classifier (from
e1071) to assign probability values to putative polyadenylation
sites (pA sites) based on training data from zebrafish. This will
allow the user to separate true, biologically relevant pA sites
from false, oligodT primed pA sites.


cleaver: In-silico cleavage of polypeptide sequences. The cleavage
rules are taken from:
http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html


clonotypeR: High throughput analysis of T cell antigen receptor
sequences The genes encoding T cell receptors are created by
somatic recombination, generating an immense combination of V,
(D) and J segments.  Additional processes during the
recombination create extra sequence diversity between the V an J
segments.  Collectively, this hyper-variable region is called the
CDR3 loop.

The purpose of this package is to process and quantitatively analyse
millions of V-CDR3-J combination, called clonotypes, from
multiple sequence libraries.

CNVrd2: CNVrd2 uses next-generation sequencing data to measure human
gene copy number for multiple samples, identify SNPs tagging copy number
variants and detect copy number polymorphic genomic regions.

cobindR: Finding and analysing co-occuring motifs of transcription
factor binding sites in groups of genes


CSSP: Power computation for ChIP-Seq data based on Bayesian
estimation for local poisson counting process.


customProDB: Generate customized protein sequence database from
RNA-Seq data for proteomics search


dagLogo: Visualize significant conserved amino acid sequence pattern
in groups based on probability theory


DNaseR: Strand-specific digital genomic footprinting in DNase-seq
data. The cumulative Skellam distribution function (package
'skellam') is used to detect significant normalized count
differences of opposed sign at each DNA strand. This is done in
order to determine the protein-binding footprint flanks.
Preprocessing of the mapped reads is recommended before running
DNaseR (e.g., quality checking and removal of sequence-specific
bias).


EBSeq: Differential Expression analysis at both gene and isoform
level using RNA-seq data


epivizr: This package provides Websocket communication to the epiviz
web app (http://epiviz.cbcb.umd.edu) for interactive
visualization of genomic data. Objects in R/bioc interactive
sessions can be displayed in genome browser tracks or plots to be
explored by navigation through genomic regions. Fundamental
Bioconductor data structures are supported (e.g., GenomicRanges
and SummarizedExperiment objects), while providing an easy
mechanism to support other data structures. Visualizations (using
d3.js) can be easily added to the web app as well.


exomePeak: The package is developed for the analysis of
affinity-based epitranscriptome shortgun sequencing data from
MeRIP-seq (maA-seq). It was built on the basis of the exomePeak
MATLAB package (Meng, Jia, et al. "Exome-based analysis for RNA
epigenome sequencing data." Bioinformatics 29.12 (2013):
1565-1567.) with new functions for differential analysis of two
experimental conditions to unveil the dynamics in
post-transcriptional regulation of the RNA methylome. The
exomePeak R-package accepts and statistically supports multiple
biological replicates, internally removes PCR artifacts and
multi-mapping reads, outputs exome-based binding sites (RNA
methylation sites) and detects differential post-transcriptional
RNA modification sites between two experimental conditions in
term of percentage rather the absolute amount. The package is
still under active development, and we welcome all biology and
computation scientist for all kinds of collaborations and
communications. Please feel free to contact Dr. Jia Meng
<jia.meng@hotmail.com> if you have any questions.


FGNet: Build and visualize functional gene networks from clustering
of enrichment analyses in multiple annotation spaces. The package
includes an interface to perform the analysis through David and
GeneTerm Linker.


flipflop: Flipflop discovers which isoforms of a gene are expressed
in a given sample together with their abundances, based on
RNA-Seq read data.


flowBeads: This package extends flowCore to provide functionality
specific to bead data. One of the goals of this package is to
automate analysis of bead data for the purpose of normalisation.


flowFit: This package estimate the proliferation of a cell population
in cell-tracking dye studies. The package uses an R
implementation of the Levenberg-Marquardt algorithm (minpack.lm)
to fit a set of peaks (corresponding to different generations of
cells) over the proliferation-tracking dye distribution in a FACS
experiment.


flowMap: This package provides an algorithm to compare and match cell
populations across multiple flow cytometry samples. The method is
based on the Friedman-Rafsky test, a nonparametric multivariate
statistical test, where two cell distributions match if they
occupy a similar feature space. The algorithm allows the users to
specify a reference sample for comparison or to construct a
reference sample from the available data. The output of the
algorithm is a set of text files where the cell population labels
are replaced by a metaset of population labels, generated from
the matching process.


GOSim: This package implements several functions useful for computing
similarities between GO terms and gene products based on their GO
annotation. Moreover it allows for computing a GO enrichment
analysis


h5vc: This package contains functions to interact with tally data
from NGS experiments that is stored in HDF5 files. For detail see
the webpage at http://www.ebi.ac.uk/~pyl/h5vc.


intansv: This package provides efficient tools to read and integrate
structural variations predicted by popular softwares. Annotation
and visulation of structural variations are also implemented in
the package.


interactiveDisplay: The interactiveDisplay package contains the
methods needed to generate interactive Shiny based display
methods for Bioconductor objects.


maPredictDSC: This package implements the classification pipeline of
the best overall team (Team221) in the IMPROVER Diagnostic
Signature Challenge. Additional functionality is added to compare
27 combinations of data preprocessing, feature selection and
classifier types.


metaSeq: The probabilities by one-sided NOISeq are combined by
Fisher's method or Stouffer's method


methylMnM: To give the exactly p-value and q-value of MeDIP-seq and
MRE-seq data for different samples comparation.


mitoODE: The package contains the methods to fit a cell-cycle model
on cell count data and the code to reproduce the results shown in
the paper "Dynamical modelling of phenotypes in a genome-wide
RNAi live-cell imaging assay" (submitted).


msmsEDA: Exploratory data analysis to assess the quality of a set of
LC-MS/MS experiments, and visualize de influence of the involved
factors.


msmsTests: Statistical tests for label-free LC-MS/MS data by spectral
counts, to discover differentially expressed proteins between two
biological conditions. Three tests are available: Poisson GLM
regression, quasi-likelihood GLM regression, and the negative
binomial of the edgeR package.The three models admit blocking
factors to control for nuissance variables.To assure a good level
of reproducibility a post-test filter is available, where we may
set the minimum effect size considered biologicaly relevant, and
the minimum expression of the most abundant condition.


MSstats: A set of tools for protein significance analysis in
label-free or LC-MS, SRM and DIA experiments.


mzID: A parser for mzIdentML files implemented using the XML package.
The parser tries to be general and able to handle all types of
mzIdentML files with the drawback of having less 'pretty' output
than a vendor specific parser. Please contact the maintainer with
any problems and supply an mzIdentML file so the problems can be
fixed quick.


neaGUI: neaGUI is an easy to use R package developed to perform the
network enrichment analysis (NEA) proposed by Alexeyenko et al.
(2012). The NEA method extends the overlap statistics in GSEA to
network links between genes in the experimental set and those in
the functional categories by exploiting biological information in
terms of gene interaction network. The neaGUI requires the
following R packages: tcltk, KEGG.db, GO.db, reactome.db,
org.Hs.eg.db, AnnotationDbi, and hwriter.


NetSAM: The NetSAM (Network Seriation and Modularization) package
takes an edge-list representation of a network as an input,
performs network seriation and modularization analysis, and
generates as files that can be used as an input for the
one-dimensional network visualization tool NetGestalt
(http://www.netgestalt.org) or other network analysis.


omicade4: Multiple co-inertia analysis of omics datasets


OmicCircos: OmicCircos is an R application and package for generating
high-quality circular maps for omic data


openCyto: This package is designed to facilitate the automated gating
methods in sequential way to mimic the manual gating strategy.


paircompviz: This package provides visualization of the results from
the multiple (i.e. pairwise) comparison tests such as
pairwise.t.test, pairwise.prop.test or pairwise.wilcox.test. The
groups being compared are visualized as nodes in Hasse diagram.
Such approach enables very clear and vivid depiction of which
group is significantly greater than which others, especially if
comparing a large number of groups.


pathifier: Pathifier is an algorithm that infers pathway deregulation
scores for each tumor sample on the basis of expression data.
This score is determined, in a context-specific manner, for every
particular dataset and type of cancer that is being investigated.
The algorithm transforms gene-level information into
pathway-level information, generating a compact and biologically
relevant representation of each sample.


plethy: This package provides the infrastructure and tools to import,
query and perform basic analysis of whole body plethysmography
and metabolism data.  Currently support is limited to data
derived from Buxco respirometry instruments as exported by their
FinePointe software.


ProCoNA: Protein co-expression network construction using peptide
level data, with statisical analysis. (Journal of Clinical
Bioinformatics 2013, 3:11 doi:10.1186/2043-9113-3-11)


prot2D: The purpose of this package is to analyze (i.e. Normalize and
select significant spots) data issued from 2D GEl experiments


PSICQUIC: PSICQUIC is a project within the HUPO Proteomics Standard
Initiative (HUPO-PSI).  It standardises programmatic access to
molecular interaction databases.


qcmetrics: The package provides a framework for generic quality
control of data. It permits to create, manage and visualise
individual or sets of quality control metrics and generate
quality control reports in various formats.


qusage: This package is an implementation the Quantitative Set
Analysis for Gene Expression (QuSAGE) method described in (Yaari
G. et al, Nucl Acids Res, 2013). This is a novel Gene Set
Enrichment-type test, which is designed to provide a faster, more
accurate, and easier to understand test for gene expression
studies. qusage accounts for inter-gene correlations using the
Variance Inflation Factor technique proposed by Wu et al.
(Nucleic Acids Res, 2012). In addition, rather than simply
evaluating the deviation from a null hypothesis with a single
number (a P value), qusage quantifies gene set activity with a
complete probability density function (PDF). From this PDF, P
values and confidence intervals can be easily extracted.
Preserving the PDF also allows for post-hoc analysis (e.g.,
pair-wise comparisons of gene set activity) while maintaining
statistical traceability. Finally, while qusage is compatible
with individual gene statistics from existing methods (e.g.,
LIMMA), a Welch-based method is implemented that is shown to
improve specificity. For questions, contact Chris Bolen
(cbolen1@gmail.com) or Steven Kleinstein
(steven.kleinstein@yale.edu)


Rchemcpp: The Rchemcpp package implements the marginalized graph
kernel and extensions, Tanimoto kernels, graph kernels,
pharmacophore and 3D kernels suggested for measuring the
similarity of molecules.


RDAVIDWebService: Tools for retrieving data from the Database for
Annotation, Visualization and Integrated Discovery (DAVID) using
Web Services into R objects. This package offers the main
functionalities of DAVID website including: i) user friendly
connectivity to upload gene/background list/s, change
gene/background position, select current specie/s, select
annotations, etc. ii) Reports of the submitted Gene List,
Annotation Category Summary, Gene/Term Clusters, Functional
Annotation Chart, Functional Annotation Table


rfPred: Based on external numerous data files where rfPred scores are
pre-calculated on all genomic positions of the human exome, the
package gives rfPred scores to missense variants identified by
the chromosome, the position (hg19 version), the referent and
alternative nucleotids and the uniprot identifier of the protein.
Note that for using the package, the user has to be connected on
the Internet or to download the TabixFile and index
(approximately 3.3 Go).


Roleswitch: Infer Probabilities of MiRNA-mRNA Interaction Signature
(ProMISe) using paired expression data from a single sample.
Roleswitch operates in two phases by inferring the probability of
mRNA (miRNA) being the targets ("targets") of miRNA (mRNA),
taking into account the expression of all of the mRNAs (miRNAs)
due to their potential competition for the same miRNA (mRNA). Due
to dynamic miRNA repression in the cell, Roleswitch assumes that
the total transcribed mRNA levels are higher than the observed
(equilibrium) mRNA levels and iteratively updates the total
transcription of each mRNA targets based on the above inference.


RRHO: The package is aimed at inference on the amount of agreement in
two sorted lists using the Rank-Rank Hypergeometric Overlap test.


RTN: This package provides classes and methods for transcriptional
network inference and analysis. Modulators of transcription
factor activity are assessed by conditional mutual information,
and master regulators are mapped to phenotypes using different
strategies, e.g., gene set enrichment, shadow and synergy
analyses.


rTRM: rTRM identifies transcriptional regulatory modules (TRMs) from
protein-protein interaction networks.


rTRMui: This package provides a web interface to compute
transcriptional regulatory modules with rTRM.


seqCNA: Copy number analysis of high-throughput sequencing cancer
data with fast summarization, extensive filtering and improved
normalization


SeqVarTools: An interface to the fast-access storage format for VCF
data provided in SeqArray, with tools for common operations and
analysis.


shinyTANDEM: This package provides a GUI interface for rTANDEM. The
GUI is primarily designed to visualize rTANDEM result object or
result xml files. But it also provides an interface for creating
parameter objects, launching searches or performing conversions
between R objects and xml files.


SigFuge: Algorithm for testing significance of clustering in RNA-seq
data.


SimBindProfiles: SimBindProfiles identifies common and unique binding
regions in genome tiling array data. This package does not rely
on peak calling, but directly compares binding profiles processed
on the same array platform. It implements a simple threshold
approach, thus allowing retrieval of commonly and differentially
bound regions between datasets as well as events of compensation
and increased binding.


SpacePAC: Identifies clustering of somatic mutations in proteins via
a simulation approach while considering the protein's tertiary
structure.


spliceR: An R package for classification of alternative splicing and
prediction of coding potential from RNA-seq data.


spliceSites: Align gap positions from RNA-seq data


sRAP: This package provides a pipeline for gene expression analysis
(primarily for RNA-Seq data).  The normalization function is
specific for RNA-Seq analysis, but all other functions (Quality
Control Figures, Differential Expression and Visualization, and
Functional Enrichment via BD-Func) will work with any type of
gene expression data.


sSeq: The purpose of this package is to discover the genes that are
differentially expressed between two conditions in RNA-seq
experiments. Gene expression is measured in counts of transcripts
and modeled with the Negative Binomial (NB) distribution using a
shrinkage approach for dispersion estimation. The method of
moment (MM) estimates for dispersion are shrunk towards an
estimated target, which minimizes the average squared difference
between the shrinkage estimates and the initial estimates. The
exact per-gene probability under the NB model is calculated, and
used to test the hypothesis that the expected expression of a
gene in two conditions identically follow a NB distribution.


STRINGdb: The STRINGdb package provides a user-friendly interface to
the STRING protein-protein interactions database (
http://www.string-db.org ).


supraHex: A supra-hexagonal map is a giant hexagon on a 2-dimensional
grid seamlessly consisting of smaller hexagons. It is supposed to
train, analyse and visualise a high-dimensional omics data. The
supraHex is able to carray out gene/meta-gene clustering and
sample correlation, plus intuitive visualisations to facilitate
exploratory analysis. Uniquely to this package, users can
simultaneously understand their own omics data in a
sample-specific fashion but without loss of information on large
genes.


SwimR: SwimR is an R-based suite that calculates, analyses, and plots
the frequency of C. elegans swimming behavior over time.  It
places a particular emphasis on identifying paralysis and
quantifying the kinetic elements of paralysis during swimming.
Data is input to SwipR from a custom built program that fits a 5
point morphometric spine to videos of single worms swimming in a
buffer called Worm Tracker.


TargetScore: Infer the posterior distributions of microRNA targets by
probabilistically modelling the likelihood
microRNA-overexpression fold-changes and sequence-based scores.
Variaitonal Bayesian Gaussian mixture model (VB-GMM) is applied
to log fold-changes and sequence scores to obtain the posteriors
of latent variable being the miRNA targets. The final targetScore
is computed as the sigmoid-transformed fold-change weighted by
the averaged posteriors of target components over all of the
features.


TCC: This package provides a series of functions for performing
differential expression analysis from RNA-seq count data using
robust normalization strategy (called DEGES). The basic idea of
DEGES is that potential differentially expressed genes or
transcripts (DEGs) among compared samples should be removed
before data normalization to obtain a well-ranked gene list where
true DEGs are top-ranked and non-DEGs are bottom ranked. This can
be done by performing a multi-step normalization strategy (called
DEGES for DEG elimination strategy). A major characteristic of
TCC is to provide the robust normalization methods for several
kinds of count data (two-group with or without replicates,
multi-group/multi-factor, and so on) by virtue of the use of
combinations of functions in other sophisticated packages
(especially edgeR, DESeq, and baySeq).


TFBSTools: Software package for TFBS.


tRanslatome: Detection of differentially expressed genes (DEGs) from
the comparison of two biological conditions (treated vs.
untreated, diseased vs. normal, mutant vs. wild-type) among
different levels of gene expression (transcriptome ,translatome,
proteome), using several statistical methods: Rank Product,
t-test, SAM, Limma, ANOTA, DESeq, edgeR. Possibility to plot the
results with scatterplots, histograms, MA plots, standard
deviation (SD) plots, coefficient of variation (CV) plots.
Detection of significantly enriched post-transcriptional
regulatory factors (RBPs, miRNAs, etc) and Gene Ontology terms in
the lists of DEGs previously identified for the two expression
levels. Comparison of GO terms enriched only in one of the levels
or in both. Calculation of the semantic similarity score between
the lists of enriched GO terms coming from the two expression
levels. Visual examination and comparison of the enriched terms
with heatmaps, radar plots and barplots.


trio: Testing SNPs and SNP interactions with a genotypic TDT. This
package furthermore contains functions for computing pairwise
values of LD measures and for identifying LD blocks, as well as
functions for setting up matched case pseudo-control genotype
data for case-parent trios in order to run trio logic regression,
for imputing missing genotypes in trios, for simulating
case-parent trios with disease risk dependent on SNP interaction,
and for power and sample size calculation in trio data.


vtpnet: variant-transcription factor-phenotype networks, inspired by
Maurano et al., Science (2012), PMID 22955828


XVector: Memory efficient S4 classes for storing sequences
"externally" (behind an R external pointer, or on disk).




NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages. The following package NEWS is available:



ADaCGH2
-------

Changes in version 2.1.6 (2013-09-23):

- Started testing against R beta 3.0.2. Fixed Imports and Depends.

- chromData: using "short", not "ushort", to catch more user errors.

- pSegmentGLAD: using a custom daglad that minimizes unneeded calls,
  specially sorting, that can be very expensive.

- Added "certain_noNA" to segmentation methods.

- Started testing against R-devel (to become R-3.1.0).

Changes in version 2.1.5 (2013-09-16):

- Fixed typos.

- Minimized usage of ":::", removed unused functions for Ansari, and
  some assignemts that no longer made sense (all packages now have
  namespaces).

- Minimize "Depends" and use "Suggests" and "Imports" in DESCRIPTION
  with "importFrom" in NAMESPACE.

- No longer using our own mergeLevels, since identical to the ones in
  aCGH package.

- GLAD uses now the recommended fastest option (smoothfunc=haarseg).

Changes in version 2.1.4 (2013-07-01):

- Fixed missing entry in bib of vignette.

Changes in version 2.1.3 (2013-06-20):

- Default merging of pSegmentDNAcopy changed to "MAD", to reflect our
  usage.

- Added more benchmarking results and recommendations to the vignette,
  and fixed some typos.

Changes in version 2.1.2 (2013-06-17):

- More changes in cutFile to try and get it to run under Mac.

- Fixed names in long examples that were leading to mistakenly
  reporting results as different.

- Added new benchmarking results.

Changes in version 2.1.1 (2013-06-16):

- Many small changes and adaptations in vignette and help to get it to
  work unded Win and Mac.

- Changes in cutFile to try and get it to run under Mac.

Changes in version 2.1.0 (2013-05-30):

- This is a major rewrite of a most of the code, has new functions,
  major changes in existing functions, new vignettes, etc.

- We no longer use snowfall.

- Major changes in parallelization, using forking.

- Reading of data: many more options, parallelized reading.

affxparser
----------

Changes in version 1.34.0 (2012-10-14):

- The version number was bumped for the Bioconductor release version,
  which is now Bioc v2.13 for R (>= 3.0.0).

Changes in version 1.33.4 (2013-09-23):

- SPEEDUP/CLEANUP: Package now uses which() instead of whichVector() of
  'R.utils'.  Before R (< 2.11.0), which() used to be 10x slower than
  whichVector(), but now it's 3x faster.

Changes in version 1.33.3 (2013-06-29):

- Same updates as in release v1.32.3.

Changes in version 1.33.2 (2013-05-25):

- Same updates as in release v1.32.2.

Changes in version 1.33.1 (2013-05-20):

- Same updates as in release v1.32.1.

Changes in version 1.33.0 (2013-04-03):

- The version number was bumped for the Bioconductor devel version.

- No updates.

Changes in version 1.32.3 (2013-06-29):

- BUG FIX: Since affxparser v1.30.2/1.31.2 (r72352; 2013-01-08),
  writeCdf() would incorrectly encode the unit types, iff the input
  'cdf' argument specified them as integers, e.g. as done by writeCdf()
  for AffyGenePDInfo in aroma.affymetrix.  More specifically, the unit
  type index would be off by one, e.g. an 'expression' unit (1) would
  be encoded as an 'unknown' unit (0) and so on.  On the other hand, if
  they were specified by their unit-type names (e.g. 'expression') the
  encoding should still be correct, e.g. if input is constructed from
  readCdf() of affxparser. Thanks to Guido Hooiveld at Wageningen UR
  (The Netherlands) for reporting on this.

- BUG FIX: Similarily, writeCdf() has "always", at least affxparser
  v1.7.4 since (r21888; 2007-01-09), encoded unit directions and QC
  unit types incorrectly, iff they were specified as integers.

Changes in version 1.32.2 (2013-05-25):

- SPEEDUP: Removed all remaining gc() calls.

- SPEEDUP: Replaced all rm() calls with NULL assignments.

Changes in version 1.32.1 (2013-05-20):

- CRAN POLICY: Now all Rd \usage{} lines are at most 90 characters
  long.

annmap
------

Changes in version 1.2.1:

- BUG ANNMAP-112 transcriptCoordsToGenome and proteinCoordsToGenome
  both failed when returning a GRanges object as I had strand as a
  numeric, rather than an integer.

AnnotationHub
-------------

Changes in version 1.2.0:

NEW FEATURES

- Support for TabixFile VCF files

BUG FIXES

- incorrect tab completion does not reset the completion line

aroma.light
-----------

Changes in version 1.32.0 (2012-10-14):

- The version number was bumped for the Bioconductor release version,
  which now is Bioc v2.13 for R (>= 3.0.0).

Changes in version 1.31.10 (2013-10-08):

- Added averageQuantile() for matrices in addition to lists.

- SPEEDUP: Now normalizeQuantileSpline(..., sortTarget=TRUE) sorts the
  target only once for lists of vectors just as done for matrices.

- DOCUMENTATION: Merged the documentation for normalizeQuantileSpline()
  for all data types into one help page.  Same for plotXYCurve().

- BUG FIX: Argument 'lwd' of plotXYCurve(X, ...) was ignored if 'X' was
  a matrix.

- Bumped up package dependencies.

Changes in version 1.31.9 (2013-10-07):

- Now library(aroma.light, quietly=TRUE) attaches the package
  completely silently without any messages.

- Now the 'aroma.light' Package object is also available when the
  package is only loaded (but not attached).

- DOCUMENTATION: Merged the documentation for normalizeQuantileRank()
  for numeric vectors and lists.

- DOCUMENTATION: Now documention S3 methods under their corresponding
  generic function.

Changes in version 1.31.8 (2013-10-02):

- DOCUMENTATION: More generic functions are now "aliased" under
  relevant corresponding methods.

Changes in version 1.31.7 (2013-09-27):

- SPEEDUP: Now all package functions utilizes 'matrixStats' functions
  where possible, e.g. anyMissing(), colMins(), and
  rowWeightedMedians().

- Bumped up package dependencies.

Changes in version 1.31.6 (2013-09-25):

- CLEANUP: Package no longer use a fallback attachment of the 'R.oo'
  package upon attachment.

Changes in version 1.31.5 (2013-09-23):

- ROBUSTNESS: Now properly declaring all S3 methods in the NAMESPACE
  file.

- SPEEDUP/CLEANUP: normalizeTumorBoost() now uses which() instead of
  whichVector() of 'R.utils'.  Before R (< 2.11.0), which() used to be
  10x slower than whichVector(), but now it's 3x faster.

- CLEANUP: Now only using 'Authors@R' in DESCRIPTION, which is possible
  since R (>= 2.14.0).  Hence the new requirement on the version of R.

- Bumped up package dependencies.

Changes in version 1.31.4 (2013-09-10):

- CLEANUP: Now package explicitly imports what it needs from
  matrixStats.

- Bumped up package dependencies.

Changes in version 1.31.3 (2013-05-25):

- SPEEDUP: Removed all remaining gc() calls, which were in
  normalizeQuantileSpline().

- SPEEDUP: Replaced all rm() calls with NULL assignments.

- Updated the package dependencies.

Changes in version 1.31.2 (2013-05-20):

- Same updates as in v1.30.2.

Changes in version 1.31.1 (2011-04-18):

- Same updates as in v1.30.1.

Changes in version 1.31.0 (2013-04-03):

- The version number was bumped for the Bioc devel version.

Changes in version 1.30.5 (2013-09-25):

- Backport from v1.31.5: Declaring all S3 methods in NAMESPACE.

- Backport from v1.31.5: normalizeTumorBoost() now uses which(), which
  also removes one dependency on 'R.utils'.

Changes in version 1.30.4 (2013-09-25):

- Backport from v1.31.4: Now package explicitly imports what it needs
  from matrixStats.

Changes in version 1.30.3 (2013-09-25):

- Backport from v1.31.3: Removal of all gc() calls and removal of
  variables is now faster.

- Removed one stray str() debug output in robustSmoothSpline().

Changes in version 1.30.2 (2013-05-20):

- CRAN POLICY: Now all Rd \usage{} lines are at most 90 characters
  long.

Changes in version 1.30.1 (2013-04-18):

- Now backtransformPrincipalCurve() preserves dimension names.

- BUG FIX: backtransformPrincipalCurve() gave an error if the pricipal
  curve was fitted using data with missing values.

- BUG FIX: fitPrincipalCurve() would not preserve dimension names if
  data contain missing values.

ASSET
-----

Changes in version 1.0:

USER VISIBLE CHANGES

- Initial release of the ASSET package.

BUG FIXES

- None.

PLANS

- The next release will have new features for h.traits and h.types.

BAGS
----

Changes in version 2.1:

- Considering on having a function for time series data.

BaseSpaceR
----------

Changes in version 1.1:

NEW FEATURES

- New high level function performOAuth() makes the App authentication
  easier.

- initializeAuth() now fires up a browser window and starts the OAuth
  v2 process. Users can use 'useBrowser' parameter to control this
  behaviour.

SIGNIFICANT USER-VISIBLE CHANGES

- requestAccessToken() now returns an integer giving the status of
  request: * 1 if the request was succesful and an access_token was
  received * -1 if the AppAuth handle already contains a feasible
  access token * An integer larger than one if the request fails. There
  is also a new 'verbose' parameter controling the messages shown to
  the user. By default 'verbose = TRUE'

BUG FIXES

- Added end of line character to some messages.

- fileItem$Size is now stored as a numeric, thus allowing for file
  larger than 2GB.

BiGGR
-----

Changes in version 0.99.3:

- Initial release of BiGGR in Bioconductor. Package was formerly
  available on CRAN.

bigmemoryExtras
---------------

1.3.1: # Thu 06-20-2013 - 11:29

1.4.0: # Wed 06-26-2013 - 19:56

1.4.1:

bioassayR
---------

Changes in version 1.0.0 (2013-10-15):

NEW FEATURES

- Provides a pre-built database containing PubChem Bioassay bioactivity
  data from hundreds of thousands of assays

NEW FEATURES

- Features for identifying target selective compounds

NEW FEATURES

- Features for identifying druggable protein targets

NEW FEATURES

- S4 classes and database structure for large bioactivity data

Biobase
-------

Changes in version 2.21:

USER VISIBLE CHANGES

- channelNames<-,NChannelSet,*-methods allow re-naming channels

USER VISIBLE CHANGES

- NChannelSet validity requires all assayDataElementNames() to be
  levels in varMetadata()$channel.

BiocParallel
------------

Changes in version 0.99.0:

NEW FEATURES

- mclapply(), pvec() require only length, [, and (for mclapply) [[.

- pvectorize() creates a parallel version of its vectorized function
  argument.

- MulticoreParam, SnowParam, DoparParam (foreach-derived), SerialParam
  to parameterize back-ends.

- bplapply, bpvec as parallel evaluation models.

- bpstart, bpstop, bpisup for back-end management.

- bpvec has a new argument AGGREGATE, a function to specify how results
  are to be combined.

SIGNIFICANT USER-VISIBLE CHANGES

- BPPARM is now used as the argument name for passing BiocParallelParam
  instances to functions.

- bplapply and bpvec now only dispatch on X and BPPARAM.

BiocStyle
---------

Changes in version 1.0.0:

USER VISIBLE CHANGES

- Rename \Rpkg{} as \CRANpkg{} to reflect functionality

BUG FIXES

- avoid option conflict with \usepackage{colors}

biomvRCNS
---------

Changes in version 1.1.9:

BUG FIXES

- NA bug when model fails

BUG FIXES

- using the length of the segment if shorter than maxk

Changes in version 1.1.8:

NEW FEATURES

- add error catching in biomvRhsmm

BUG FIXES

- updates in package vignette

NEW FEATURES

- fine tuning of memory usage for large count data structure

BUG FIXES

- minor bug fixes

Changes in version 1.1.7:

NEW FEATURES

- add support of Rle like sparse data in mcols of large count data
  structure.

NEW FEATURES

- add mclapply for parallel processing of different seqnames in
  biomvRhsmm

Changes in version 1.1.6:

BUG FIXES

- fix nbinom emis estimateSegCommonDisp -> nbinomCLLDD

Changes in version 1.1.5:

NEW FEATURES

- update biomvRGviz to have trackList returned when tofile=FALSE

BUG FIXES

- typo fix in function argument

Changes in version 1.1.4:

NEW FEATURES

- update biomvRGviz to plot multiple samples in one, and more param

- plot method now defaults to multiple samples in one, sampleInOne=T

BUG FIXES

- fix states rank when prior.m=cluster

Changes in version 1.1.3:

NEW FEATURES

- the x slot in the returning object for biomvRhsmm now also contains
  the probability of associated state for each position

NEW FEATURES

- the param slot in the returning object for biomvRhsmm now also
  contains the updated emission and sojourn parameters for each seqname
  and sample column

Changes in version 1.1.2:

NEW FEATURES

- add a new clustering method for emission prior estimation

- add a new example section of DMR detection in the vignette

BUG FIXES

- switch off several data checking warnings, use normal cat() instead

Changes in version 1.1.1:

NEW FEATURES

- add support of direct input of sojourn parameters as named list.

BUG FIXES

- suppress warnings generated when appending GRanges object for the
  output

BitSeq
------

Changes in version 1.5.5 (2013-10-12):

BUG FIXES

- fixed problem with processing Fasta reads

BUG FIXES

- fixed problem with processing mixed alignments in empirical
  estimation

BUG FIXES

- use proper includes for samtools

Changes in version 1.5.2 (2013-10-09):

BUG FIXES

- fixed force option in estimateHyperPar (it worked the wrong way
  around)

- fixed problem with effective lengths being too small

NEW FEATURES

- add estimateVBExpression function that uses Variational Bayes
  inference method to estimate expression levels of transcripts

Changes in version 1.4.3 (2013-07-08):

BUG FIXES

- fixed major bug in getGeneExpression and getWithinGeneExpression that
  caused wrong results for some genes if transcripts were not grouped
  by genes

NEW FEATURES

- fixed bug in parseAlignment which would occasionally cause underflow
  in effective length computation

BUG FIXES

- fixed taking square root of wrong column in getExpression result's
  $means field

NEW FEATURES

- much faster parseAlignment computation

BUG FIXES

- improved precision in getGeneExpression and getWithinGeneExpression

NEW FEATURES

- getWithinGeneExpression provides option to keep transcripts in
  original order even if they are not grouped by genes

bsseq
-----

Changes in version 0.9:

- Fixed a problem with "width" in the title of bsseq plots.

- plot.BSseqTstat now allows for BSseqTstat objects computed without
  correction.

- validObject(BSseq) has been extended to also check for sampleNames
  consistency.

- Fixed a bug related to validity checking.

- Increased maxk from 10,000 to 50,000 in calls to locfit, to allowing
  fitting the model on genomes with unusally long chromosomes (Thanks
  to Brian Herb for reporting).

- The class representation for class 'BSseq' has changed completely.
  The new class is build on 'SummarizedExperiment' from GenomicRanges
  instead of 'hasGRanges'.  Use 'x <- updateObject(x)' to update
  serialized (saved) objects.

- Fixing a problem in orderBSseq related to chromosome names.

- Allowed user specification of maxk, with a default of 10,000 in
  BSmooth.

- Many bugfixes made necessary by the new class representation.

- Better argument checking in BSmooth.tstat.

- A few undocumented functions are now documented.

- Rewrote orderBSseq

bumphunter
----------

Changes in version 1.1:

- Added NEWS file.

- Fixed a bug related to >= for numerics.

- Added smoothing using a gaussian kernal as implemented in the locfit
  package through the function locfitByCluster.

- Added closeSockets for cleanup for doParallel on Windows.

- More bugfixes for windows; now using foreachCleanup().

- Added a 'bumps' class and print method.

- annotateNearest / regionMatch now give NA annotations for queries
  with no nearest subject (perhaps because the seqname is missing from
  the subject). Previously this was taken to be mistaken input and a
  hard error raised.

- Speedup of fwer computations using foreach.

- Added boundedClusterMaker.

- bug fix to internal function .getModT (which are not used in the main
  bumphunter functions).  Now the t-statistics returned a correct.

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

1.5.1:

CellNOptR
---------

Changes in version 1.8.0 (2013-08-28):

- BUG FIX
  
  * readMIDAS: DV, DA and TR can now be in the specy name
  
  * prep4sim can read self loop
  
  * plotOptimResultsPan: fix special cases with one or no
  inhibitors/stimuli
  
  * fix bug in cutNONC: notMat was not populated
  
  * readSIF: can read and gates coded in big caps as well as small ones
  
  * writeMIDAS: manages absence of inhibitors

- CHANGES
  
  * CNOlist: subfield parameter has been removed. Subfield are
  automatically found from the header
  
  * expandAndGates are not limited to 4 inputs anymore
  
  * normaliseCNOlist: EC50 is set to 1 by default

- NEW FEATURES
  
  * readSBMLQual: a function to read prior knowledge network in
  SBMLqual format. !! This is a prototype. use with care for now.
  
  * cutCNOlist function can cut a MIDAS file over time

ChemmineOB
----------

Changes in version 1.0.0:

NEW FEATURES

- functions to convert between any two formats supported by OpenBabel

NEW FEATURES

- descriptor computation

NEW FEATURES

- fingerprint computation

chimera
-------

Changes in version 1.3.15:

NEW FEATURES

- GapFiller Nadalin et al BMC Bioinformatics was implemented as de novo
  validation tool for fusion break point.

Changes in version 1.3.12:

NEW FEATURES

- Picard-tools can be downloaded in chimara folder, picardInstallation.

NEW FEATURES

- The wrapper function validateSamFile uses picard tools to validate
  SAM/BAM files.

NEW FEATURES

- The wrapper function filterSamReads allows SAM and BAM filtering.

Changes in version 1.3.9:

NEW FEATURES

- Rsubread aligner is used instead of tophat.

Changes in version 1.3.3:

NEW FEATURES

- STAR fusion import function: importFusionData

NEW FEATURES

- STAR run function: starInstallation

NEW FEATURES

- STAR installation function: starRun

NEW FEATURES

- BiocParallel supported functions: STAR import fusionName
  supportingReads filterList

chipenrich
----------

Changes in version 1.0:

PKG FEATURES

- chipenrich performs gene set enrichment tests on peaks called from a
  ChIP-seq experiment

PKG FEATURES

- chipenrich empirically corrects for confounding factors such as the
  length of genes and mappability of sequence surrounding genes

PKG FEATURES

- Use multiple definitions of a gene "locus" when testing for
  enrichment, or provide your own definition

PKG FEATURES

- Test for enrichment using chipenrich or Fisher's exact test (should
  only be used for datasets where peaks are close to TSSs, see docs)

PKG FEATURES

- Test multiple sets of genesets (Gene Ontology, KEGG, Biocarta, OMIM,
  etc.)

PKG FEATURES

- Multiple plots to describe binding distance and likelihood of a peak
  as a function of gene length

PKG FEATURES

- Support for human (hg19), mouse (mm9), and rat (rn4) genomes

PKG FEATURES

- Many conveniences such as seeing which peaks were assigned to genes,
  their position relative to those genes and their TSS, etc.

PKG FEATURES

- See how many peaks were assigned to each gene along with the length
  and mappability of the gene

Changes in version 0.99.2:

USER-VISIBLE CHANGES

- Updated examples for various functions to be runnable (removed
  donttest)

- Updated DESCRIPTION to use Imports: rather than Depends:

- Updated license to GPL-3

- Updated NEWS file for bioconductor guidelines

BUG FIXES

- Added a correction for the case where a small gene set has a peak in
  every gene. This has the result of making a very few number of tests
  slightly conservative, at the benefit of actually being able to
  return a p-value for them.

Changes in version 0.99.1:

USER-VISIBLE CHANGES

- Minor updates to documentation for Bioconductor

Changes in version 0.99.0:

NEW FEATURES

- Initial submission to Bioconductor

Changes in version 0.9.6:

NEW FEATURES

- Added peaks per gene as a returned object / output file

Changes in version 0.9.5:

BUG FIXES

- Update to handle bioconductor/IRange's new "functionality" for
  distanceToNearest and distance

USER-VISIBLE CHANGES

- Changed sorting of results to put enriched terms first (sorted by
  p-value), then depleted (also sorted by p-value)

Changes in version 0.9.4:

USER-VISIBLE CHANGES

- Minor changes to vignette and documentation

Changes in version 0.9.3:

NEW FEATURES

- Addition of rat genome

BUG FIXES

- chipenrich() will correctly open both .bed and .bed.gz files now

Changes in version 0.9.2:

NEW FEATURES

- Added ability for user to input their own locus definition file (pass
  the full path to a file as the locusdef argument)

- Added a data frame to the results object that gives the
  arguments/values passed to chipenrich, also written to file
  *_opts.tab

- For FET and chipenrich methods, the outcome variable can be recoded
  to be >= 1 peak, 2 peaks, 3 peaks, etc. using the num_peak_threshold
  parameter

- Added a parameter to set the maximum size of gene set that should be
  tested (defaults to 2000)

USER-VISIBLE CHANGES

- Previously only peak midpoints were given in the peak --> gene
  assignments file, now the original peak start/ends are also given

- Updated help/man with new parameters and more information about the
  results

BUG FIXES

- Fixed an issue where status in results was not enriched if the odds
  ratio was infinite, and depleted if the odds ratio was exactly zero

Changes in version 0.9.1:

NEW FEATURES

- Added a QC plot for expected # of peaks and actual # of peaks vs.
  gene locus length. This will be automatically created if qc_plots is
  TRUE, or the plots can be created using the plot_expected_peaks
  function.

NEW FEATURES

- Distance to TSS is now signed for upstream (-) and downstream (+) of
  TSS

NEW FEATURES

- Column added to indicate whether the geneset is enriched or depleted

Changes in version 0.9:

NEW FEATURES

- Added support for reading BED files natively

BUG FIXES

- Fixed bug where invalid geneset in chipenrich() wasn't detected
  properly

Changes in version 0.8:

BUG FIXES

- Fixed crash when mappability contained an NA (will be removed from DB
  in future version)

Changes in version 0.7:

USER-VISIBLE CHANGES

- Updated binomial test to sum gene locus lengths to get genome length
  and remove genes that are not present in the set of genes being
  tested

- Updated spline fit plot to take into account mappability if requested
  (log mappable locus length plotted instead of simply log locus
  length)

- Removed SAMPLEABLE_GENOME* constants since they are no longer needed

- Updated help files to reflect changes to plot_spline_length and
  chipenrich functions

BUG FIXES

- Fixed bug where results for multiple gene set types (e.g. doing
  BioCarta and KEGG together) were not sorted by p-value

Changes in version 0.6:

BUG FIXES

- Fixed bug where 1kb/5kb locusdefs could fail if not all peaks were
  assigned to a gene

Changes in version 0.5:

USER-VISIBLE CHANGES

- Updated help to explain new mappability model

USER-VISIBLE CHANGES

- Changed how mappability is handled - now multiplies gene locus length
  by mappability, rather than adjusting as a spline term

ChIPpeakAnno
------------

Changes in version 2.9.6:

BUG FIXES

- Change Enhancer.Silencer to intergenic.Region in
  assignChromosomeRegion.

Changes in version 2.9.5:

BUG FIXES

- fixed bug error in negative widths are not allowed when call
  findOverlappingPeaks.

cisPath
-------

Changes in version 1.1.6:

- Fix a bug

Changes in version 1.1.5:

- documentation improvements

Changes in version 1.1.4:

- Fix a bug

Changes in version 1.1.3:

- documentation improvements

Changes in version 1.1.2:

- Fix a bug

Changes in version 1.1.1:

- Several methods have been added to format PPI data downloaded from
  PINA and STRING databases

cleaver
-------

Changes in version 0.99.5 (2013-07-24):

- vignette:

- remove duplicated sessionInfo entries.

Changes in version 0.99.4 (2013-07-23):

- vignette:

- use BiocStyle.

- add sessionInfo() and TOC.

Changes in version 0.99.3 (2013-07-08):

- vignette: add second BRAIN reference.

Changes in version 0.99.2 (2013-06-17):

- Replace own AAStringSetList constructor by
  Biostrings::AAStringSetList.

- man/cleaver-package.Rd: remove static date.

- vignette: add dynamic date and don't load Biostrings manually
  anymore.

- NAMESPACE: don't import from Biostrings and IRanges.

Changes in version 0.99.1 (2013-05-30):

- Add S4-methods for character, AAString, AAStringSet.

- man/cleave-methods.Rd: split table of cleavage rules to reduce table
  width.

- Extend vignette (add BRAIN and UniProt.ws based examples).

Changes in version 0.99.0 (2013-04-27):

- Initial release.

clonotypeR
----------

Changes in version 2013-10-07:

- Commited 0.99.6 to Bioconductor.

- Added unit tests for yassai_identifier().

Changes in version 2013-08-01:

- Unified the syntax of common_clonotypes and unique_clonotypes.

Changes in version 2013-05-07:

- Resubmitted 0.99.5 to Bioconductor.

- Distribute a copy of the clonotypes extracted from example_data.

- Execute all examples in the vignette.

Changes in version 2013-05-01:

- Resubmitted 0.99.4 to Bioconductor

- Moved the Markdown vignette to '/vignettes'.

- Added executable example with test data to the vignette.

Changes in version 2013-04-26:

- Resubmitted 0.99.3 to Bioconductor

- Moved extra data to 'inst/extdata'.

- Moved the Markdown vignette to 'inst/vignettes'.

Changes in version 2013-04-15:

- Resubmission after correcting yassai_identifier() and other warnings.

Changes in version 2013-01-08:

- Leaner Bioconductor package, without the wiki documentation.
  2012-.....  Charles Plessy <plessy@riken.jp>

- Many entries missing.

Changes in version 2012-10-10:

- Corrected Frenglish “improductive” with “unproductive”.  This can
  break backwards compatibility.

Changes in version 2012-09-06:

- Initial release

clusterProfiler
---------------

Changes in version 1.9.4:

- bug fixed of readGff, add default parameter fill=TRUE in read.table
  function <2013-07-09, Mon>

Changes in version 1.9.3:

- extend enrichGO to support 20 species <2013-07-09, Mon>

- update vignettes. <2013-07-09, Mon>

Changes in version 1.9.1:

- enrichGO and enrichKEGG support rat organism <2013-05-20, Mon>

- change some code according to DOSE <2013-03-27, Wed>

- modify enrichGO and enrichKEGG according to the change of
  enrich.internal, remove qvalueCutoff parameter, add pAdjustMethod,
  add universe paramter for user to specify background. <2013-05-29,
  Wed>

- add function viewKEGG for visualizing KEGG pathway and update
  vignette. <2013-06-14, Fri>

CNAnorm
-------

1.7.1:

CNORfuzzy
---------

Changes in version 1.4.0:

- cSimulator handles non integers values for the ihibitors and stimuli

- gaDiscreteT1.R: fix issue when only 1 model was returned within
  tolerance

- reduceFuzzy.R fix bug that causes seg faut (model was not cut
  properly)

- defaultParametersFuzzy.R: added nTF to set number of TF to arbitrary
  value (not tested)

- CNORwrapFuzzy.R: fixed pMutation argument that was not populated

- gaDiscrete functions return best score as well in the dataframe.

- output names of the fields returned by gaDiscrete are now using camel
  lower case so that plotFit from CellNoptR can be used

- add C simulator

codelink
--------

Changes in version 1.30.0:

- The most visible change is that the CodelinkSet interface has been
  adopted as the official supported system. Documentation of these
  topic has been improved, and a new vignette describing the
  CodelinkSet system is available (Codelink_Introduction).
  Documentation to the old Codelink class has moved to the
  Codelink_Legacy vignette.

- Before, readCodelink() would assign NA values to spots flagged as M
  (MSR spot), I (Irregular) and C (Contaminated). This could cause in
  large datasets that many spots would have at least one sample with NA
  at random, reducing drastically the number of de facto spots/probes
  used during normalization. Many thanks to Emmanuelle Dantony for
  spotting this problem and subsequent feedback. Because of this and
  other problems implementing and appropriate method to deal with this,
  automatic assigment of NA values is not performed anymore. The only
  exception are M flagged spots, which have intensity values of -9999,
  and hence do not represent any measure of intensity. Also, background
  and normalization methods are applied calling the appropriate
  functions in the limma package. Support for type- and flag-based
  weights has been included, and weights are automatically created by
  readCodelinkSet(). Weights can be used to modulate the contribution
  of some probes to normalization and during linear modeling more
  efficiently. Examples on how to use these approaches are documented
  in the vignette Codelink_Introduction.

- Added generic method normalize() for Codelink and CodelinkSet
  classes.

DART
----

Changes in version 1.7.1:

- Dependency changed from igraph0 to igraph (>=0.6.0) due to archiving
  of igraph0. In file PruneNet.R, the use of clusters() output was
  modified from zero-index to one-index accordingly.

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

Bugfixes

- Including shearwater vignette

Updates

- fixed issues with deletions in bf2Vcf()

Bugfixes

- makePrior() adds background on all sites

Changes in version 1.99.0 (2013-04-30):

Updates

- New shearwater algorithm

Updates

- Including VCF output through summary(deepSNV, value="VCF")

Changes in version 1.7.4 (2013-09-28):

Updates

- Only using the Dirichlet prior for control

Updates

- Change back version numbers from 1.99.x to 1.7.x (2.0 not quite there
  :)

deltaGseg
---------

Changes in version 1.1.3:

- corrected errors generated in denoiseSegments when segments are too
  small

- refined plotting with data that has many subpopulations (>10)

- more robust argument selection in functions

DESeq2
------

Changes in version 1.1.32:

- By default, use QR decomposition on the design matrix X. This
  stabilizes the GLM fitting. Can be turned off with the useQR argument
  of nbinomWaldTest() and nbinomLRT().

- Allow for "frozen" normalization of new samples using previous
  estimated parameters for the functions: estimateSizeFactors(),
  varianceStabilizingTransformation(), and rlogTransformation(). See
  manual pages for details and examples.

Changes in version 1.1.31:

- The adjustment of p-values and use of Cook's distance for outlier
  detection is moved to results() function instead of nbinomWaldTest(),
  nbinomLRT(), or DESeq(). This allows the user to change parameter
  settings without having to refit the model.

Changes in version 1.1.24:

- The results() function allows the user to specify a contrast of
  coefficients, either using the names of the factor and levels, or
  using a numeric contrast vector. Contrasts are only available for the
  Wald test differential analysis.

Changes in version 1.1.23:

- The results() function automatically performs independent filtering
  using the genefilter package and optimizing over the mean of
  normalized counts.

Changes in version 1.1.21:

- The regularized log transformation uses the fitted dispersions
  instead of the MAP dispersions. This prevents large, true log fold
  changes from being moderated due to a large dispersion estimate blind
  to the design formula. This behavior is also more consistent with the
  variance stabilizing transformation.

DEXSeq
------

Changes in version 2013-09-12:

- Major changes to vignette to provide more of an end-to-end
  description of the work flow.

- Major changes to function names to now make TRT rather than BM the
  default; changed vignette to reflect this.

- Added apepndix explaining TRT to vignette.

Changes in version 2013-02-27:

- A parameter -r was added to the python scripts, that allow the users
  either to ignore the exonic bins belonging to several genes and treat
  the genes separately, or merge the genes into an aggregate gene. The
  equivalent R implementations of the python scripts were finally
  added.

Changes in version 2012-11-28:

- The TRT method is implemented, for people with a big number of
  samples without completions or speed issues

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

Changes in version 1.8.0:

- Add support for DESeq2:
  
  * New: Add DBA_DESEQ2, DBA_ALL_METHODS and DBA_ALL_BLOCK method
  constants
  
  * Change: dba.analyze can analyze using DESeq2
  
  * Change: all reporting and plotting functions support DESeq2 results
  
  * Change: vignette includes comparison of edgeR, DESeq, and DESeq2

- Changes to counting using dba.count:
  
  * Change: optimize built-in counting code to use much less memory and
  run faster
  
  * Change: deprecate bLowMem, replaced with bUseSummarizeOverlaps
  
  * New: add readFormat parameter to specify read file type (instead of
  using file suffix)

- New: generation of result-based DBA object using dba.report (makes it
  easier to work with differentially bound peaksets)

- Changes to defaults:
  
  * Change: default score is now DBA_SCORE_TMM_MINUS_FULL instead of
  DBA_SCORE_TMM_MINUS_EFFECTIVE in dba.count
  
  * Change: default value for bFullLibrarySize is now TRUE in
  dba.analyze
  
  * New: add bCorPlot option to DBA$config to turn off CorPlot by
  default

- Various bugfixes, improved warnings, updated documentation

DOSE
----

Changes in version 1.99.6:

- bug fixed in EXTID2NAME. <2013-09-28, Sat>

Changes in version 1.99.5:

- fixed in calculating M when only one categroy presented, the object
  was matrix insted of list. <2013-09-16, Mon>

Changes in version 1.99.4:

- export gsea function <2013-07-10, Wed>

Changes in version 1.99.3:

- extend EXTID2NAME to support 20 species <2013-07-09, Tue>

- update vignette. <2013-07-09, Tue>

Changes in version 1.99.1:

- convert vignette file to knitr Sweave. <2013-06-24, Mon>

Changes in version 1.99.0:

- extent ggplot to support enrichResult by implementing fortify method.
  <2013-05-22, Wed>

- re-implement barplot.enrichResult. <2013-05-23, Thu>

- enrich.internal support user specifiy background by parameter
  universe. <2013-05-24, Fri>

- implement Gene Set Enrichment Analysis algorithm. <2013-05-29, Wed>

- change setReadable to support groupGO of clusterProfiler.
  <2013-05-29, Wed>

- fixed mclapply not support Windows platform issue. <2013-05-30, Fri>

- rename logFC parameter to foldChange. <2013-06-13, Thu>

Changes in version 1.7.1:

- use geom_bar(stat="identity") instead of geom_bar() in barplot for
  explicitly mapping y value. <2013-05-08, Wed>

- bug fixed when qvalue can't calculated. <2013-05-02, Thu>

- bug fixed of enrich.internal, drop those genes that without
  annotation when calculating sample gene number. <2013-05-02, Thu>

- change some code to satisfy ReactomePA <2013-03-27, Wed>

DSS
---

Changes in version 1.7.0:

- Added functionalities for multiple factor experimental design.
  Currently edgeR functions are used for glm fitting and hypothesis
  testing. DSS only provide functions for dispersion estimation.

Changes in version 1.6.0:

- Implemented methods for detecting differentially methylated loci
  (DML).

EBImage
-------

Changes in version 4.4.0:

NEW FEATURES

- New colorLabels function color-coding labels of object masks by a
  random permutation (Bernd Fisher)

- Additional argument inputRange to normalize allowing presetting a
  limited input intensity range

- Additional argument thick to paintObjects controlling the thickness
  of boundary contours

SIGNIFICANT USER-VISIBLE CHANGES

- normalize and combine use the generics from BiocGenerics

- removed the along argument from combine

- Re-introduced calculation of 's.radius.sd' (standard deviation of the
  mean radius) in cell features

BUG FIXES

- getFrame: XY dimensions equal 1 were dropped

edgeR
-----

Changes in version 3.3.8:

- predFC() with design=NULL now uses normalization factors correctly.
  However this use of predFC() to compute counts per million is being
  phased out in favour of cpm().

Changes in version 3.3.5:

- Refinement to cutWithMinN() to make the bin numbers more equal in the
  worst case.

- estimateDisp() now creates the design matrix correctly when the
  design matrix is not given as an argument and there is only one
  group.  Previously this case gave an error.

- Minor edit to glm.h code.

Changes in version 3.3.4:

- plotMDS.DGEList now gives a friendly error message when there are
  fewer than 3 data columns.

Changes in version 3.3.3:

- DGEList() accepts NULL as a possible value again for the group,
  lib.size and norm.factors arguments. It is treated the same way as a
  missing argument.

Changes in version 3.3.2:

- Update to cutWithMinN() so that it does not fail even when there are
  many repeated x values.

- Refinement to computation for nbins in dispBinTrend.  Now changes
  more smoothly with the number of genes.  trace argument is retired.

- Fixes to calcNormFactors with method="TMM" so that it takes account
  of lib.size and refCol if these are preset.

- Updates to help pages for the data classes.

Changes in version 3.3.1:

- Updates to DGEList() and DGEList-class documentation. Arguments
  lib.size, group and norm.factors are now set to their defaults in the
  function definition rather than set to NULL.

eiR
---

Changes in version 1.2.0:

NEW FEATURES

- speed improvements

NEW FEATURES

- ieInit now accepts a SNOW cluster for parallel inserts

NEW FEATURES

- allow compounds to be updated in-place by name

NEW FEATURES

- eiQuery can now return similarity values instead of distances

fabia
-----

Changes in version 2.6.1:

NEW FEATURES

- rescaling of lapla

NEW FEATURES

- extractPlot does not plot sorted matrices

FGNet
-----

Changes in version 1.0.0:

NEW FEATURES

- Package released

flowType
--------

Changes in version 2.0.0:

NEW FEATURES

- New flowType using an enhanced dynamic programming algorithm
  implemented in C++, allowing for searching only a subset of the
  phenotype space (to make very high-dimensional data like CyTOF and
  single cell RT-PCR feasible), and to allow the setting of multiple
  levels of expression.

NEW FEATURES

- Note: this version brings some minor but necessary changes to the way
  in which the flowType function is called, which may break backwards
  compatibility for some users. Also note that the return value from
  flowType is no longer compatible with version 1.x of RchyOptimyx;
  users wishing to use RchyOptimyx with flowType 2.x should use the new
  and improved version 2.x of RchyOptimyx, which is actually simpler
  and much more intuitive.

fmcsR
-----

Changes in version 1.4.0:

NEW FEATURES

- Package and fmcsR algorithm published in Bioinformatics (Sep 4, 2013,
  Epub ahead of print)

NEW FEATURES

- fmcsR outperforms most other virtual screening (VP) methods

gage
----

Changes in version 2.11.3:

- add secondary vignette, "RNA-Seq Data Pathway and Gene-set Analysis
  Workflows".

- add function kegg.gsets, which generates up-to-date KEGG pathway gene
  sets for any specified KEGG species.

gCMAP
-----

Changes in version 1.5.4:

- NEW: The bigmemory and bigmemoryExtras packages are not optional,
  enabling use of gCMAP on Windows

Changes in version 1.5.3:

- NEW: Updated titles of vignettes

Changes in version 1.5.2:

- BUGFIX: Updated eSet construction in mapNmerge function

Changes in version 1.5.1:

- BUGFIX: fixed incorrect varMetaData element definition in
  .process_counts

gCMAPWeb
--------

Changes in version 1.1.7:

- BUGFIXFixed memory leak and enabled html_table function to deal with
  multibyte strings.

Changes in version 1.1.6:

- UPDATEgCMAPWeb now issues an error when javascript is not available
  in the browser.

Changes in version 1.1.2:

- NEW: Added tutorial as a new vignette.

Changes in version 1.1.1:

- BUGFIX: Heatmaps are now displayed for results of directional query.

geNetClassifier
---------------

Changes in version 1.1.1:

BUG FIXES

- plotDiscriminantPower now accepts feature names starting by number
  (i.e. affy probes)

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

Changes in version 1.14:

NEW FEATURES

- keys method now has new arguments to allow for more sophisticated
  filtering.

- adds genes() extractor

- makeTranscriptDbFromGFF() now handles even more different kinds of
  GFF files.

BUG FIXES

- better argument checking for makeTranscriptDbFromGFF()

- cols arguments and methods will now be columns arguments and methods

GenomicRanges
-------------

Changes in version 1.14.0:

NEW FEATURES

- Add coercion from GenomicRangesList to RangedDataList.

- Add "c" method for GAlignmentsPairs objects.

- Add coercion from GAlignmentPairs to GAlignmentsList.

- Add 'inter.feature' and 'fragment' arguments to summarizeOverlaps().

- Add seqselect,GAlignments-method.

- Add CIGAR utilities: explodeCigarOps(), explodeCigarOpLengths()
  cigarRangesAlongReferenceSpace(), cigarRangesAlongQuerySpace()
  cigarRangesAlongPairwiseSpace(), extractAlignmentRangesOnReference()
  cigarWidthAlongReferenceSpace(), cigarWidthAlongQuerySpace()
  cigarWidthAlongPairwiseSpace().

- Add seqlevels0() and restoreSeqlevels().

- Add seqlevelsInUse() getter for GRanges, GRangesList, GAlignments
  GAlignmentPairs, GAlignmentsList and SummarizedExperiment objects.

- Add update,GAlignments method.

- Add GIntervalTree class.

- Add coercion from GAlignmentPairs to GAlignments.

- Add sortSeqlevels().

- Add tileGenome().

- Add makeGRangesFromDataFrame() and coercion from data.frame or
  DataFrame to GRanges.

SIGNIFICANT USER-LEVEL CHANGES

- Renaming (with aliases from old to new names): - classes
  GappedAlignments -> GAlignments GappedAlignmentPairs ->
  GAlignmentPairs - functions GappedAlignments() -> GAlignments()
  GappedAlignmentPairs() -> GAlignmentPairs() readGappedAlignments() ->
  readGAlignments() readGappedAlignmentPairs() -> readGAlignmentPairs()

- Remove 'asProperPairs' argument to readGAlignmentsList().

- Modify "show" method for Seqinfo object to honor showHeadLines and
  showTailLines global options.

- 50x speedup or more when merging 2 Seqinfo objects, 1 very small and
  1 very big.

- Add dependency on new XVector package.

- Enhanced examples for renaming seqlevels in seqlevels-utils.Rd.

- More efficient reference class constructor for 'assays' slot of
  SummarizedExperiment objects.

- 'colData' slot of SummarizedExperiment produced from call to
  summarizedOverlaps() now holds the class type and length of 'reads'.

- 4x speedup to cigarToRleList().

- Relax SummarizedExperiment class validity.

- Renaming (with aliases from old to new names): cigarToWidth() ->
  cigarWidthOnReferenceSpace(), and cigarToQWidth() ->
  cigarWidthOnQuerySpace().

- Improvements to summarizeOverlaps(): - mode 'Union': 1.5x to 2x
  speedup - mode 'IntersectionNotEmpty': 2x to 8x speedup + memory
  footprint reduced by ~ half

- Change default 'use.names' to FALSE for readGAlignmentsList().

- Implement 'type="equal"' for findOverlaps,SummarizedExperiment
  methods.

- Modify summarizeOverlaps() examples to use 'asMates=TRUE' instead of
  'obeyQname=TRUE'.

- Remove unneeded "window" method for GenomicRanges objects.

- Speed up seqinfo() getter and setter on SummarizedExperiment objects
  and derivatives (e.g. VCF) by using direct access to 'rowData' slot.

- coverage,GenomicRanges method now uses .Ranges.coverage() when using
  the defaults for 'shift' and 'width'.

- Remove restriction that metadata column names must be different on a
  GRangesList and the unlisted GRanges.

- GenomicRangesUseCases vignette has been redone and renamed to
  GenomicRangesHOWTOs

DEPRECATED AND DEFUNCT

- Defunct all "match" and "%in%" methods in the package except for
  those with the GenomicRanges,GenomicRanges signature.

- Deprecate GappedAlignment*: - GappedAlignments and
  GappedAlignmentPairs classes - GappedAlignments() and
  GappedAlignmentPairs() constructors - readGappedAlignments() and
  readGappedAlignmentPairs() functions

- Deprecate cigar util functions: cigarToWidth(), cigarToQWidth(),
  cigarToIRanges() splitCigar(), cigarToIRanges(),
  cigarToIRangesListByAlignment() cigarToIRangesListByRName(),
  cigarToWidth(), cigarToQWidth() cigarToCigarTable(),
  summarizeCigarTable()

- Deprecate seqselect().

BUG FIXES

- Fix bug in c,GAlignments for case when objects were unnamed.

- Fix bug in flank,GenomicRanges (when 'ignore.strand=TRUE' 'start' was
  being set to TRUE).

- Fix bug in behavior of summarizeOverlaps() count mode
  'IntersectionNotEmpty' when 'inter.features=FALSE'. Shared regions
  are now removed before counting.

- Fix bug in cigarToIRangesListByAlignment() when 'flag' is supplied
  and indicates some reads are unmapped.

- Fix bug in summarizeOverlaps(..., mode='IntersectionNotEmpty') when
  'features' has '-' and '+' elements and 'ignore.strand=TRUE'.

- match,GenomicRanges,GenomicRanges method now handles properly objects
  with seqlevels not in the same order.

GGtools
-------

Changes in version 4.10:

- Sharply revised vignette

Changes in version 4.9:

- All.cis accepts a CisConfig instance to define parameters of a cis
  search

gmapR
-----

Changes in version 1.4.0:

NEW FEATURES

- Add desired_read_group to BamTallyParam; will limit tallies to that
  specific read group (useful for multi-amplicon sequencing, like
  Fluidigm)

- Add keep_ref_rows argument to variantSummary() for keeping rows for
  positions where no alt is detected (the rows where ref == alt).

- gsnap() will now output a GsnapOutputList when there are multiple
  input files

- Support 'terminal_threshold' and 'gmap_mode' parameters in
  GsnapParam, and use different defaults for DNA vs. RNA. This means a
  big improvement in alignment quality for DNA.

- GmapGenome now accepts a direct path to the genome as its first
  argument

USER-VISIBLE CHANGES

- Renamed summarizeVariants to variantSummary

- The 'which' in GsnapParam is now a GenomicRanges instead of a
  RangesList

- Refactor bam_tally, so that bam_tally returns a TallyIIT object,
  which is then summarized via summarizeVariants; this allows computing
  tallies once and summarizing them in different ways (like maybe get
  the coverage). The summarizeVariants function yields a VRanges.

BUG FIXES

- fix minimum quality cutoff check to >=, instead of >

- fix asBam,GsnapOutput for when unique_only is TRUE

- package created by makeGmapGenomePackage now have a GmapGenome with
  the correct name

GOSemSim
--------

Changes in version 1.19.3:

- add getSupported_Org function for accessing all the names of
  supported organisms <2013-07-09, Mon>

Changes in version 1.19.2:

- export getDb and loadGOMap <2013-07-09, Mon>

- update vignettes <2013-07-9, Mon>

Changes in version 1.19.1:

- update vignettes <2013-06-13, Thu>

GSEABase
--------

Changes in version 1.23:

SIGNIFICANT USER-VISIBLE CHANGES

- Support 'c7' broad set

SIGNIFICANT USER-VISIBLE CHANGES

- Warn or remove duplicate geneIds when parsing GMT or Broad XML

GWASTools
---------

Changes in version 1.7.8:

- Added gdsSetMissingGenotypes, updated argument names in
  ncdfSetMissingGenotypes.

- Changed colorscheme in manhattanPlot.R.

- Bug fix in ibdPlot - diagonal orange bars are back.

- Bug fix in plinkWrite for writing just one sample.

- Bug fix in printing pedigreeCheck error message.

Changes in version 1.7.7:

- Changed handling of GxE interaction variables in assocTestRegression.

Changes in version 1.7.6:

- Updated vignette for SNPRelate 0.9.16.

Changes in version 1.7.5:

- gwasExactHW will run on all chromosomes except (Y,M), rather than
  (autosome,X,XY) only.

Changes in version 1.7.4:

- More informative error messages in anomDetectBAF and anomDetectLOH.

Changes in version 1.7.3:

- Changed labeling of IBD plots from "HS" to "Deg2" and "FC" to "Deg3."

- Bug fix in pedigreePairwiseRelatedness - no more warning about
  multiple values passed to stringsAsFactor.

- pedigreeClean and pedigreeFindDuplicates are now defunct.  Use
  pedigreeCheck instead.

hapFabia
--------

Changes in version 1.2.2:

- bug fix split_sparse_matrix

- plot functions with other arguments '...'

- plot arguments grid and pairs

- new function 'plotLarger' (add samples without IBD and borders)

- vcftoFABIA with command line options -s (SNVs_) and -o (output file)

- vcftoFABIA in R with output file name

HiTC
----

Changes in version 1.5.2:

NEW FEATURES

- Efficient memory matrix representation using the Matrix package. The
  memory usage for big sparse matrix is improved by a factor 7.
  However, some operation are much slower based on the Matrix
  implementation. Thus, for some task as the plotting function, the
  Matrix are converted in standard matrix based object

- 'show' and 'detail' method for HTClist object

- 'c(x, ...)' method for HTCexp and HTClist objects

SIGNIFICANT USER-VISIBLE CHANGES

- The option mask.data from the mapC function is deprecated

- Update of parallel computation for some functions (import, normalize)

BUG FIXES

- Bug Fixe in import.my5C for ChrM

Changes in version 1.5.1:

BUG FIXES

- Fix bug in HTCexp contructor for matrix of dim 1

BUG FIXES

- Fix bug in import.my5C for matrix of dim 1

BUG FIXES

- Fix bug for title display on HTCexp plots

BUG FIXES

- Fix bug in HiTClist plot due to chromosome order

BUG FIXES

- Correct errors in the man pages

hpar
----

Changes in version 1.3.1:

- Updated to HPA version 11 <2013-04-11 Thu>

Changes in version 1.3.0:

- Bioc 2.13 devel version bump

HTSeqGenie
----------

Changes in version 3.11.10:

- vignette updated to build on BioC dev

Changes in version 3.11.9:

- TP53GenomicFeatures() now uses the same temp dir as TP53Genome(), to
  fix a build issues on BioC servers

Changes in version 3.11.8:

- -choice of VT filters configurable via config

Changes in version 3.11.7:

- mergeBams(..., sort=FALSE) in alignReads.R: since bams are already
  sorted in chunks, merging preserves the order and no re-sorting
  should is necessary

Changes in version 3.11.6:

- saveCoverage now outputs a summary.coverage summary file

Changes in version 3.11.5:

- directly export coverage from Rle with rtracklayer::export (instead
  of to make a memory-costly convertsion to a RangedData object)

Changes in version 3.11.4:

- update variant calling code to work with VariantTools 1.3.6

- include Jens' minlength=1 modification to handle reads fully trimmed

Changes in version 3.11.3:

- exports loginfo, logdebug, logwarn, logerror

- uses TallyVariantsParam and as(,"VRanges") to fix a bug preventing
  compilation on BioC

- the number of threads used during processChunks is divided by 2 due
  to an erroneous extra mcparallel(...) step in sclapply/safeExecute

- use a maximum of 12 cores in preprocessReads

Changes in version 3.11.2:

- trim lowest quality tail of reads according to Illumina manual.  This
  new quality step precedes the regular quality filtering that is still
  done.

Changes in version 3.11.1:

- allow to use better filters for variant tools by default, uses
  standard filters (bad), but if repeat masked track and dbsnp are
  given in config it uses those.

Changes in version 3.11.0:

- starting the new development branch

Changes in version 3.10.1:

- release version

htSeqTools
----------

Changes in version 1.7.1:

- Added findPeakHeight function to estimate optimal minHeight value
  based on False Discovery Rate for using in enrichedPeaks

HTSFilter
---------

Changes in version 1.0.1:

- -- Updated citation information, minor updates to vignette -- DESeq2,
  IRanges, and GenomicRanges added to Imports field, and HTSFilter and
  HTSFilterBasic methods written for objects of with DESeqDataSet
  signature (from DESeq2 pipeline).  Examples illustrating the use of
  HTSFilter have also been added in the documentation and vignettes. --
  DESeq and edgeR moved to Imports field (thanks to Denis Laloë for the
  suggestion) -- Warning message added when gene ID's are included as a
  column in the data frame (thanks to Luc Jouneau for the suggestion)

hyperdraw
---------

Changes in version 1.13.1:

- Node labels were ignoring the 'cex' node data value (reported by
  Hannes Hettling)

- Added "start" and "both" options for "arrowLoc" graph attribute
  (which will draw arrowhead at both ends of hyper edges)

- Bug fix for converting Hypergraph to graphBPH so that hyperedge names
  are used as edge labels (reported by Hannes Hettling)

- Update the R version check in drawEdge() to cope properly with major
  versions greater than 2 (!)

illuminaio
----------

Changes in version 0.3.9:

- Added vignette giving examples of usage and demonstrating comparison
  with GenomeStudio output

Changes in version 0.3.8:

- Bug fixes, renaming of readIDAT_bin to readIDAT_nonenc.

Changes in version 0.3.6:

- Added one parent function readIDAT which checks file format and
  dispatches to the relevant subfunction.

Changes in version 0.3.5 (2013-08-02):

- Cleaned up internal code of readBPM().

- ROBUSTNESS: Added unit tests for readBPM().  Note that these are only
  run if environment variable '_R_CHECK_FULL_' is set, i.e. they will
  *not* be perfomed on the Bioconductor servers.

imageHTS
--------

Changes in version 1.12.0:

NEW FEATURES

- Additional arguments to highlightSegmentation controlling the color
  (col) and opacity (opac) of nuclei and cell contours, closing of
  contours at borders (border), and whether highligting of grayscale
  images should be done in color (toRGB).

SIGNIFICANT USER-VISIBLE CHANGES

- removed elongation of the nuclei mask along image borders in
  segmentATH

intansv
-------

Changes in version 0.99.3:

Notes

- Improvement of the vigenette to make it more comprehensive.

Notes

- Removal of useless example dataset.

iPAC
----

Changes in version 1.5.2:

- Updated file to be consistent with Bioconductor version numbering.
  Also fixed a rare bug that gave an error if every residue aligned
  between mutational and positional data.

IPPD
----

09-29-2010:

01-12-2011:

01-17-2012:

IRanges
-------

Changes in version 1.20.0:

NEW FEATURES

- Add IntervalForest class from Hector Corrada Bravo.

- Add a FilterMatrix class, for holding the results of multiple
  filters.

- Add selfmatch() as a faster equivalent of 'match(x, x)'.

- Add "c" method for Views objects (only combine objects with same
  subject).

- Add coercion from SimpleRangesList to SimpleIRangesList.

- Add an `%outside%` that is the opposite of `%over%`.

- Add validation of length() and names() of Vector objects.

- Add "duplicated" and "table" methods for Vector objects.

- Add some split methods that dispatch to splitAsList() even when only
  'f' is a Vector.

SIGNIFICANT USER-VISIBLE CHANGES

- All functionalities related to XVector objects have been moved to the
  new XVector package.

- Refine how isDisjoint() handles empty ranges.

- Remove 'keepLength' argument from "window<-" methods.

- unlist( , use.names=FALSE) on a CompressedSplitDataFrameList object
  now preserves the rownames of the list elements, which is more
  consistent with what unlist() does on other CompressedList objects.

- Splitting a list by a Vector just yields a list, not a List.

- The rbind,DataFrame method now handles the case where Rle and vector
  columns need to be combined (assuming an equivalence between Rle and
  vector). Also the way the result DataFrame is constructed was changed
  (avoids undesirable coercions and should be faster).

- as.data.frame.DataFrame now passes 'stringsAsFactors=FALSE' and
  'check.names=!optional' to the underlying data.frame() call.
  as(x,"DataFrame") sets 'optional=TRUE' when delegating. Most places
  where we called as.data.frame(), we now call 'as(x,"data.frame")'.

- The [<-,DataFrame method now coerces column sub-replacement value to
  class of column when the column already exists.

- DataFrame() now automatically derives rownames (from the first
  argument that has some). This is a fairly significant change in
  behavior, but it probably does better match user behavior.

- Make sure that SimpleList objects are coerced to a DataFrame with a
  single column. The automatic coecion methods created by the methods
  package were trying to create a DataFrame with one column per
  element, because DataFrame extends SimpleList.

- Change default to 'compress=TRUE' for RleList() constructor.

- tapply() now handles the case where only INDEX is a Vector (e.g.  an
  Rle object).

- Speedup coverage() in the "tiling case" (i.e. when 'x' is a tiling of
  the [1, width] interval). This makes it much faster to turn into an
  Rle a coverage loaded from a BigWig, WIG or BED as a GRanges object.

- Allow logical Rle return values from filter rules.

- FilterRules no longer requires its elements to be named.

- The select,Vector method now returns a DataFrame even when a single
  column is selected.

- Move is.unsorted() generic to BiocGenerics.

DEPRECATED AND DEFUNCT

- Deprecate seqselect() and subsetByRanges().

- Deprecate 'match.if.overlap' arg of "match" method for Ranges
  objects.

- "match" and "%in%" methods that operate on Views, ViewsList,
  RangesList, or RangedData objects (20 methods in total) are now
  defunct.

- Remove previously defunct tofactor().

BUG FIXES

- The subsetting code for Vector derivatives was substancially
  refactored.  As a consequence, it's now cleaner, simpler, and [ and
  [[ behave more consistently across Vector derivatives. Some obscure
  long-standing bugs have been eliminated and the code can be slightly
  faster in some circumstances.

- Fix bug in findOverlaps(); zero-width ranges in the query no longer
  produce hits ever (regardless of 'maxgap' and 'minoverlap' values).

- Correctly free memory allocated for linked list of results compiled
  for findOverlap(select="all").

- Various fixes for AsIs and DataFrames.

- Allow zero-row replacement values in [<-,DataFrame.

- Fix long standing segfault in "[" method for Rle objects (when doing
  Rle()[0]).

- "show" methods now display its most specific class when a column or
  slot is an S3 object for which class() returns more than one class.

- "show" methods now display properly cells that are arrays.

- Fix the [<-,DataFrame method for when a value DataFrame has matrix
  columns.

- Fix ifelse() for when one or more of the arguments are Rle objects.

- Fix coercion from SimpleList to CompressedList via AtomicList
  constructors.

- Make "show" methods robust to "showHeadLines" and "showTailLines"
  global options set to NA, Inf or non-integer values.

- Fix error condition in eval,FilterRules method.

- Corrected an error formatting in eval,FilterRules,ANY method.

isobar
------

Changes in version 1.7.6:

- added TMT 10plex (contribution from Florent Gluck)

- fixed bugs with system.file not working on R < 2.11 (contribution
  from Florent Gluck)

- fixed bug in isobar-qc which was not working without normalize=TRUE

- added writeHscoreData for usage with Hscorer.py (MM Savitski, 2010)

- shQuote commands correctly - should fix issues running report
  generation on Windows

- added calculations and XLS tab for peptides with unsure localization
  of the modification site

- updated scripts for creating multi-sample reports
  (create.meta.reports)

Changes in version 1.7.5:

- fixed critical bugs: Excel report output had wrong ordering, ie
  ratios did not correspond to the meta information [introduced in
  version 1.7.3].

- fix of real peptide names: Reexport I/L peptides in reports

Changes in version 1.7.4:

- improved MSGF+ search result import

- refactored report properties: all properties can now be defined in in
  the properties.R

- speed and memory usage improvements when creating wide XLS report

- ratio p-value adjustment now works per comparision instead of
  globally

Changes in version 1.7.3:

- fix wide XLS report format

- novel plot for ratio combinations in QC report showing individual
  ratio distributions and significant ratios

Changes in version 1.7.2:

- added TMT 6plex masses to phosphoRS export script

- fixed mascot parsers

- MzIdentML version 1.1.0 support implemented [not fully tested]

Changes in version 1.7.1:

- fixed import of MzIdentML files from Mascot and ProteomeDiscoverer

KEGGgraph
---------


KEGGprofile
-----------


lmdme
-----

Changes in version 1.3.1:

BUGS FIXED

- NA manipulation of lmFit for coefficient is now contemplated

Changes in version 1.3.0:

MINOR CHANGES

- Bump y in version x.y.z to odd number in devel

MANOR
-----

Changes in version 2013-07-03 (2013-07-03):

- Replaced calls to 'exit', 'fprintf' which now raise a WARNING upon
  check.

- Iteration index not printed anymore in 'nem' (raised an error when
  compiling the vignette).

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

Changes in version 1.2 (2013-08-20):

- Our paper got accepted and is available!

- Added methods for MRexperiment objects
  (colSums,colMeans,rowSums,rowMeans, usage is for example colSums(obj)
  or colSums(obj,norm=TRUE)) (09-25)

- Added two new functions, plotOrd and plotRare - a function to plot
  PCA/MDS coordinates and rarefaction effect (09-04,09-18)

- Updated MRfisher to include thresholding for presence-absence testing
  (08-19)

- Updated comments (roxygen2) style for all the functions using the
  Rd2roxygen package (07-13)

- Updated plotCorr and plotMRheatmap to allow various colors/not
  require trials(07-13)

- Rewrote vignette (and switched to knitr)

Changes in version 1.1 (2013-06-25):

- Rewrote load_meta and load_metaQ to be faster/use less memory

- Modified cumNormStat to remove NA samples from calculations (example
  would be samples without any counts)

- Re-added plotGenus' jitter

- Fixed uniqueNames call in the MR tables

- Changed thanks to Kasper Daniel Hansen's suggestions the following:
  plotOTU and plotGenus both have much better auto-generated axis
  MRtable, MRfulltable, MRcoefs have a sort by p-value option now
  MRtable, MRfulltable, MRcoefs now have an extra option to include
  unique numbers for OTU features (default would automatically add them
  previously) cumNorm.R - now returns the object as well - not just
  replacing the environment 0 Still need to turn the fitZig output to
  S3, consider subsetting function address low p-values

minfi
-----

Changes in version 1.7:

- Added getMethSignal(), a convenience function for programming.

- Changed the argument name of "type" to "what" for getMethSignal().

- Added the class "RatioSet", like "GenomicRatioSet" but without the
  genome information.

- Bugfixes to the "GenomicRatioSet()" constructor.

- Added the method ratioConvert(), for converting a "MethylSet" to a
  "RatioSet" or a "GenomicMethylSet" to a "GenomicRatioSet".

- Fixed an issue with GenomicMethylSet() and GenomicRatioSet() caused
  by a recent change to a non-exported function in the GenomicRanges
  package (Reported by Gustavo Fernandez Bayon <gbayon@gmail.com>).

- Added fixMethOutliers for thresholding extreme observations in the
  [un]methylation channels.

- Added getSex, addSex, plotSex for estimating sex of the samples.

- Added getQC, addQC, plotQC for a very simple quality control measure.

- Added minfiQC for a one-stop function for quality control measures.

- Changed some verbose=TRUE output in various functions.

- Added preprocessQuantile.

- Added bumphunter method for "GenomicRatioSet".

- Handling signed zero in minfi:::.digestMatrix which caused unit tests
  to fail on Windows.

- addSex and addQC lead to sampleNames() being dropped because of a
  likely bug in cbind(DataFrame, DataFrame).  Work-around has been
  implemented.

- Re-ran the test data generator.

- Fixed some Depends and Imports issues revealed by new features of R
  CMD check.

- Added blockFinder and cpgCollapse.

- (internal) added convenience functions for argument checking.

- Exposed and re-wrote getAnnotation().

- Changed getLocations() from being a method to a simple function.
  Arguments have been removed (for example, now the function always
  drops non-mapping loci).

- Implemented getIslandStatus(), getProbeType(), getSnpInfo() and
  addSnpInfo().  The two later functions retrieve pre-computed SNP
  overlaps, and the new annotation object includes SNPs based on dbSNP
  137, 135 and 132.

- Changed the IlluminaMethylatioAnnotation class to now include
  genomeBuild information as well as defaults.

- Added estimateCellCounts for deconvolution of cell types in whole
  blood.  Thanks to Andrew Jaffe and Andres Houseman.

mitoODE
-------

Changes in version 0.99.6:

- changed mc.cores to 1 to avoid Windows building error

Changes in version 0.99.5:

- changed mc.cores=12 by mc.cores=detectCores()

Changes in version 0.99.2:

- initial submission

motifStack
----------

Changes in version 1.5.4:

NEW FEATURES

- use BiocStyle to write vignette

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.5.3:

NEW FEATURES

- No changes classified as 'new features' (package under active
  development)

BUG FIXES

- fix the bug that motifStack can not plot color correctly for
  radialPhylog layout in motifStack().

Changes in version 1.5.2:

NEW FEATURES

- use log-likelihood ratio (ALLR) to do alignment

BUG FIXES

- No changes classified as 'bug fixes' (package under active
  development)

Changes in version 1.5.1:

NEW FEATURES

- plotMotifStackWithRadialPhylog can draw motif cloud in RadialPhylog
  style. Check vignette to see the details

BUG FIXES

- check mode of pcm and pfm matrix, which should be numeric

MSnbase
-------

Changes in version 1.9.12:

- fix MSnSet -> ExpressionSet <2013-10-13 Sun>

- MSnSet -> ExpressionSet unit test <2013-10-13 Sun>

Changes in version 1.9.11:

- MIAPE to MIAME conversion <2013-10-11 Fri>

- proper MIAME when MSnSet -> ExpressionSet <2013-10-11 Fri>

Changes in version 1.9.10:

- faster plotMzDelta <2013-09-28 Sat>

- faster plotMzDelta for mzRramp instances <2013-09-29 Sun>

- chromatogram method <2013-10-04 Fri>

- plotMzDelta has subset arg <2013-10-04 Fri>

- xic method <2013-10-04 Fri>

- suggesting msdata for chromatogram example <2013-10-04 Fri>

- renamed plotting arg 'centroided.' <2013-10-04 Fri>

Changes in version 1.9.9:

- typo in filterNA Rd <2013-09-18 Wed>

- writeMgfData now has a progress bar <2013-09-24 Tue>

- centroided(MSnExp) <- TRUE now allowed <2013-09-24 Tue>

Changes in version 1.9.8:

- using new.env(parent=emptyenv()) to get rid of enclosing env when
  creating new MSnExps <2013-09-17 Tue>

- new (private) MSnExp.size function <2013-09-17 Tue>

Changes in version 1.9.7:

- Passing ... to read.table in MSnbase:::readIspy[Silac|15N]Data
  <2013-09-16 Mon>

- QualityControl biocView <2013-09-16 Mon>

Changes in version 1.9.6:

- new as.data.frame.MSnSet method <2013-08-16 Fri>

- new ms2df function <2013-08-16 Fri>

- new getEcols and grepEcols helpers for readMSnSet2 <2013-08-16 Fri>

Changes in version 1.9.5:

- typo in Author[s]@R <2013-05-15 Wed>

Changes in version 1.9.4:

- new simple MSnSet constructor <2013-05-07 Tue>

Changes in version 1.9.3:

- Using knitr as VignetteEngine <2013-04-29 Mon>

- Remove LazyLoad from DESCRIPTION, which is default nowadays
  <2013-04-29 Mon>

- knitr dependency > 1.1.0 for VignetteBuilder <2013-04-29 Mon>

- Adding MSnSet creating sections in io vignette <2013-04-29 Mon>

- new readMSnSet2 function <2013-04-30 Tue>

Changes in version 1.9.2:

- clean has now a all param (default FALSE is retain original
  behavious) to remove all 0 intensity values <2013-04-17 Wed>

- using BiocGenerics::normalize <2013-04-25 Thu>

Changes in version 1.9.1:

- new logging utility function to update an MSnSet's
  processingData(object)$processing <2013-03-29 Fri>

- Proper logging in t.MSnSet <2013-03-29 Fri>

Changes in version 1.9.0:

- new Bioc 2.13 devel version

MSstats
-------

Changes in version 1.99.1:

- fixed several NOTES, added .Rbuildignore, compacted vignettes

- TODO: check remaining 'no visible binding for global variable' NOTES

- removed warn -1

- added validity check when returning MSnSet

- used inherits/is. for class testing

- TODO fix if conditions syntax

Changes in version 1.99.0:

- improve efficiency for computing groupComparison and quantification
  <2012-12-21>

- add .rnw <2012-12-03>

- update groupComparision for label-free time-course experiment with
  single Feature and with or without technical replicates <2013-04-08>

- add option for saving QQ plot and Residual plot in order to checkin
  the normality assumption in groupComparison function. <2013-04-08>

- use ggplot2 package for all plots. <2013-07-11>

- fix bug for volcano plot : different color labeling <2013-07-12>

- add power plot in sample size calculation plot <2013-07-12>

- add 'interference=TRUE or FALSE' in sample size calculation
  <2013-07-15>

- add 'legend.size=7'for size of feature name legends in
  dataProcessPlots <2013-07-23>

- add 'text.angle=0'for angle of condition labeling in dataProcessPlots
  <2013-07-23>

- fix bug for quantification : when there are missing values in
  endogenous intensities, but values in reference intensities.
  <2013-07-24>

- fix bug for groupComparison : when there are missing values in
  endogenous intensities, but values in reference intensities,
  .make.constast.based or .free sub function were changed. <2013-07-25>

- two function for transformation between required input for MSstats
  and MSnSet class <2013-09-04>

- flexibility for visualization : save as pdf files or show in window
  with selected proteins or all proteins. <2013-09-04>

- handle unequal variance for feature in groupComparison function with
  featureVar=TRUE <2013-09-04>

- Add 'missing.action' for impute missing values in group comparison
  stage. <2013-09-20>

mzID
----

Changes in version 0.3.1:

- Comment unused functions <2013-07-05 Fri>

- Minor typos in vignette <2013-07-05 Fri>

Changes in version 0.3.0:

NEW FEATURES AND FUNCTIONS

- Support for mzIdentML v1.0 (thanks Laurent Gatto)

Changes in version 0.2.1:

NEW FEATURES AND FUNCTIONS

- Version counter now in BioConductor style

NEW FEATURES AND FUNCTIONS

- Namespace now extracted from file instead of hardcoded in the
  constructor

NEW FEATURES AND FUNCTIONS

- mzID constructor now checks the mzIdentML version before parsing
  (only 1.1.0 supported atm)

Changes in version 0.1-1:

NEW FEATURES AND FUNCTIONS

- Added function `flatten` to create tabular representation of results

Changes in version 0.0-2:

NEW FEATURES AND FUNCTIONS

- The package is now fully documented

- created helper functions: `countChildren`, `attrExtract` and
  `attrExtractNameValuePair`

- Added NEWS file

- Added README.md file

- The parser now succesfully imports all test files on the HUPO
  mzIdentML page

PERFORMANCE

- Improvements in the mzIDpsm constructor makes it ~2-3x faster
  depending on the file size (thanks to `countChildren`, `attrExtract`
  and `attrExtractNameValuePair`)

Changes in version 0.0-1:

NEW FEATURES AND FUNCTIONS

- First build of the package

NEW FEATURES AND FUNCTIONS

- Introduction of the classes: `mzID`, `mzIDdatabase`, `mzIDevidence`,
  `mzIDparameters`, `mzIDpeptides` and `mzIDpsm` with related
  constructors

mzR
---

Changes in version 1.7.4:

- version bump for Rcpp 0.10.5 <2013-10-02 Wed>

Changes in version 1.7.3:

- Fix a compile error with the clang-3.3 compiler

Changes in version 1.7.2:

- updated Rcpp number mismatch warning to include versions <2013-08-01
  Thu>

Changes in version 1.7.1:

- version bump for Rcpp 0.10.4

NOISeq
------

Changes in version 2.2.0 (2013-10-14):

- New function to generate a Quality Control Report in PDF format
  including all the exploratory plots.

- Plot to evaluate RNA composition bias has been changed.

- Some bugs have been fixed.

OSAT
----

Changes in version 1.8.1:

BUG FIXES

- "exlude" arg in the gContainer construction now can be used when
  whole plates are excluded.

NEW FUNCTION

- "exlude<-" is defined as generic and can be used to exclude wells
  from an existing container object.

pathview
--------

Changes in version 1.1.7:

- Graphviz view can automatic choose different types of legends, either
  on nodes, edges or both depending on the specific pathways.

- Vigette has been reformatted: the "Quick start" section added

Changes in version 1.1.6:

- Pathview can plot/integrate/compare multiple states/samples in the
  same graph. Several functions and data objects are revised: including
  pathview, keggview.native, keggview.graph, render.kegg.node etc. New
  section on multiple state data with working exampleshas been added.

- Vigette has been reformatted: Data integration section splitted into
  multiple sub-sections.

Changes in version 1.1.5:

- Pathview works with species where default KEGG gene ID is not Entrez
  Gene. Several functions and data objects are revised: including
  pathview, node.map, sim.mol.data, kegg.species.code, korg. New
  section on KEGG species and Gene ID usage with working exampleshas
  been added.

Changes in version 1.1.3:

- Pathview paper published in Bioinformatics

phyloseq
--------


piano
-----

Changes in version 1.2.0:

NEW FEATURES

- Added argument ncpus to runGSA(). Enables parts of this function to
  run in parallel, thus decreasing runtime. Requires R package
  snowfall.

- Added function runGSAhyper() to perform gene set analysis using
  Fisher's exact test, as an alternative to runGSA.

- Added information about genes in each area of the Venn diagram in the
  output of diffExp().

- Added volcano plot as optional output of diffExp() and added argument
  volcanoFC.

- Added argument ncharLabel to networkPlot() and consensusHeatmap() to
  control the length of the labels in the plots and add the option to
  not truncate them.

- Added the yeast metabolic model iTO977 to be loaded with loadGSC(),
  for detecting reporter metabolites using gene set analysis. (See
  vignette.)

- Added support for running networkPlot() on objects returned by
  runGSAhyper().

USER-VISIBLE CHANGES

- Minor updates to the vignette.

- Updated diffExp() man page.

- Updated the main legend of the plot from consensusScores() to make it
  clearer.

- Updated error-messages in networkPlot().

- Fixed typo in PC variance plot produced by runQC().

- Restructered this NEWS file.

BUG FIXES

- Updated diffExp() to handle changes in lmFit() and topTable() from
  limma.

- Removed suppressWarnings(), as temporarily introduced in version
  1.0.1, in polarPlot() around the calls to radial.plot() since
  warnings are now fixed in package plotrix.

Changes in version 1.0.7:

USER-VISIBLE CHANGES

- Updated consensusScores() man page.

Changes in version 1.0.6:

USER-VISIBLE CHANGES

- Updated loadMAdata() man page.

Changes in version 1.0.5:

BUG FIXES

- Updated diffExp() to only use vennDiagram() of the limma package for
  Venn diagram plotting in order to correct a bug when plotting more
  than three circles. Also updated the corresponding man page.

USER-VISIBLE CHANGES

- Updated the SBML section in the loadGSC() man page.

Changes in version 1.0.4:

BUG FIXES

- Fixed bug in loadMAdata so that also compressed CEL-files (*.CEL.gz)
  can be loaded correctly.

Changes in version 1.0.3:

USER-VISIBLE CHANGES

- Updated references in vignette.

Changes in version 1.0.2:

USER-VISIBLE CHANGES

- Updated CITATION information.

Changes in version 1.0.1:

USER-VISIBLE CHANGES

- Updated CITATION information.

- Fixed typos in DESCRIPTION and piano-package.Rd.

- Updated the man file for loadGSC().

- Added URL in DESCRIPTION.

BUG FIXES

- Fixed bug in loadGSC() so that gmt-files are now loaded correctly.

- Temporarily added suppressWarnings() in polarPlot() around the calls
  to radial.plot() since warnings appeared for example("radial.plot")
  in plotrix v3.4-6.

plethy
------

Changes in version 0.99.4:

SIGNIFICANT USER-VISIBLE CHANGES

- The 'BuxcoR' package is now renamed 'plethy' as of version 0.99.3.

NEW FEATURES

- Added 'get.labels.by.sample', 'get.err.breaks' and 'adjust.labels' to
  assist with QA/QC of data.

plgem
-----

Changes in version 1.33.1:

- minor changes: --removed unnecessary calls to `which' within
  subsetting operations, to avoid potential unintended consequences
  pointed out by Wolfgang Huber:
  <https://stat.ethz.ch/pipermail/bioc-devel/2013-May/004353.html>

ProCoNA
-------

Changes in version 0.99.2:

- removed a random component in the unit tests that could cause
  failure.

Changes in version 0.99.1:

- Major update ... moved to S4 methods (getters, setters, setMethods,
  etc).

- Smaller sample data set for quick example runs in the man pages.

- Integration with the MSnbase package. See the vignette for an example
  of raw data processing in preparation for network building. The main
  buildProconaNetwork function now takes either matrices of peptide
  data, or MSnSet objects containing data.

- Man pages relating to the main proconaNet object have been combined.

pRoloc
------

Changes in version 1.1.8:

- semi-sup and sup comparison sections in ml vig <2013-10-09 Wed>

Changes in version 1.1.7:

- getMarkers has names arg <2013-10-01 Tue>

Changes in version 1.1.6:

- new outliers arg to plot2D <2013-09-23 Mon>

- cite addMarkers in vignette <2013-09-26 Thu>

- add code chunk in poster vig <2013-09-26 Thu>

Changes in version 1.1.5:

- added biocViews <2013-09-19 Thu>

- added knitr vig engine in ml <2013-09-19 Thu>

- import dependencies <2013-09-20 Fri>

Changes in version 1.1.4:

- building ml vignette in Makefile <2013-09-11 Wed>

Changes in version 1.1.3:

- new private anyUnknown function, used in phenoDisco, to check for the
  presence of unlabelled ('unknown') data <2013-08-27 Tue>

- added a note in vignette about "unknown" convention to define protein
  with unknown localisation <2013-08-28 Wed>

- Added HUPO 2011 poster <2013-09-09 Mon>

Changes in version 1.1.2:

- fixed Author[s]@R <2013-05-16 Thu>

- na.rm=TRUE in f1Count <2013-05-19 Sun>

- added CITATION file <2013-06-07 Fri>

- new testMarkers function <2013-06-29 Sat>

- error message in getBestParams suggests to use testMarkers
  <2013-06-29 Sat>

- nndist methods (see issue #23) <2013-07-01 Mon>

- remove unnecessary as.matrix in plot2D's cmdscale <2013-07-19 Fri>

- added plot2D(..., method = "MDS") example <2013-07-19 Fri>

- changed 'euclidian' to 'euclidean' in nndist_matrix <2013-07-26 Fri>

- fixed row ordering in phenoDisco, input and output rownames are now
  the same <2013-08-13 Tue>

- Using filterNA in phenoDisco <2013-08-13 Tue>

Changes in version 1.1.1:

- Update README.md to reflect availability in stable release and
  biocLite installation <2013-04-07 Sun>

- new MSe pipeline <2013-04-07 Sun>

- perTurbo using Rcpp <2013-04-16 Tue>

- initial clustering infrastructure (not exported) <2013-04-17 Wed>

- new markerSet function <2013-04-24 Wed>

- plsaOptim's ncomp is now 2:6 <2013-04-27 Sat>

- added forword to vignette <2013-04-29 Mon>

- default k is now seq(3, 15, 2) in knnOptim <2013-05-09 Thu>

Changes in version 1.1.0:

- Bioc 2.13 devel version bump

PWMEnrich
---------

Changes in version 2.4.4:

- Vignette update and fix naming of columns in the motif enrichment
  report

Changes in version 2.4.2:

- Improve promoter selection for human and mouse genomes (duplicates
  are now disregarded)

Changes in version 2.4.0:

- Major update with more functions and small bugfixes

- Added sequenceReport() and groupReport() for easier report generation

- Visualise motif scores along a sequence with plotMotifScores()

- Creation of empirical CDFs for motif scores

- Almost complete rewrite of the vignette to emphasize the main use
  cases

- Converted documenation to roxygen2

Changes in version 2.3.2:

- Subsetting functions for backgrounds from Diego Diez

Changes in version 2.3.1:

- Fix a bug with plotTopMotifsSequence() with calling an unknown
  function

- Implement group.only for all background, not only pval in
  motifEnrichment()

- New default to plotMultipleMotifs() so the margins are better

qcmetrics
---------

Changes in version 0.99.3:

- fixed template - toc after begin document <2013-10-01 Tue>

- updates to vignette <2013-10-01 Tue>

- re-reoxygenise rnadeg man <2013-10-01 Tue>

Changes in version 0.99.2:

- Implemented Andrzej's suggestions <2013-09-27 Fri>

Changes in version 0.99.1:

- Updated github README file <2013-09-18 Wed>

- Added Arnoud's suggestions to vig <2013-09-21 Sat>

- rnadeg wrapper function available in the package <2013-09-21 Sat>

- new QcMetadata class <2013-09-21 Sat>

- metadata and mdata synonym <2013-09-21 Sat>

- added metadata section in knitr reports <2013-09-21 Sat>

- selectively import pander::pander.table to fix warning uppon
  installation <2013-09-21 Sat>

- added n15qc wrapper <2013-09-21 Sat>

Changes in version 0.99.0:

- submission to Bioconductor <2013-09-17 Tue>

qpgraph
-------

Changes in version 1.18:

USER VISIBLE CHANGES

- Updated the formatting of the vignettes to adhere to the Bioconductor
  style

- Functions qpRndHMGM() and qpSampleFromHMGM() which were deprecated in
  the previous release, are now defunct.

NEW FEATURES

- qpCItest() takes now also R/qtl cross objects as input, i.e., the
  user can test for conditional independence directly on R/qtl cross
  objects

- added functions addGenes(), addeQTL(), and addGeneAssociation() to
  incrementally build and simulate eQTL networks

BUG FIXES

- Real or integer valued levels in discrete variables now prompt an
  error when they are not positive integers

- qpFunctionalCoherence() handles regulatory modules with just one
  target without giving an error

- reQTLcross() can now simulate an initial eQTL model with no genes

- fixed plotting for HMgmm objects

QuasR
-----

Changes in version 1.2.0:

NEW FEATURES

- new miRNA-seq sample data; vignette has now workflows for ChIP-seq,
  RNA-seq, smRNA-seq/miRNA-seq, Bis-seq and allele-specific analysis

NEW FEATURES

- qMeth can report methylation states for individual molecules
  (reportLevel="alignment")

NEW FEATURES

- qCount/qProfile gain argument to ignore spliced reads in counting

Rbowtie
-------

Changes in version 1.1.1:

NEW FEATURES

- updated bowtie to version 1.0.0

RchyOptimyx
-----------

2.0:

RDAVIDWebService
----------------

Changes in version 0.99.1:

MINOR CHANGES

- `dontrun` tags have been removed from plotting examples (Thanks to
  Valerie Obenchain).

Changes in version 0.99.0:

DOCUMENTATION

- `NEWS` file was added.

ReactomePA
----------

Changes in version 1.5.3:

- bug fixed in TERMID2EXTID <2013-09-18, Wed>

Changes in version 1.5.2:

- implement gseAnalyzer for GSEA analysis <2013-07-10, Wed>

- implement viewPathway for visualizing pathway <2013-07-10, Wed>

Changes in version 1.5.1:

- extend enrichPathway to support rat, mouse, celegans, zebrafish, and
  fly. <2013-03-27, Wed>

- modify enrichPathway according to the change of enrich.internal,
  remove qvalueCutoff parameter, add pAdjustMethod, add universe
  paramter for user to specify background. <2013-05-29, Wed>

RedeR
-----

Changes in version 2.0.0:

- Improved graph attribute settings and interface.

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

Rgraphviz
---------

Changes in version 2.5:

- Expose mai argument in plot,Ragraph-method.

- Changed the LICENSE to EPL.

- Fixed some wrong text in the DESCRIPTION file.

- Fixed pieGylph to adress the issue of warning messages about
  rep(NULL) (reported by Cristobal Fresno Rodríguez
  <Cristobalfresno@gmail.com>).

- Updated Imports, Depends, NAMESPACE as per R CMD check requirements.

- Fixing issue with the lt~obsolete.m4 file(s) in Graphviz; R CMD check
  was issueing a warning.

- Moved vignettes from inst/doc to vignettes directory.

rhdf5
-----

Changes in version 2.6.0:

NEW FEATURES

- support for logical added

- support for reading attributes added (use read.attributes=TRUE)

- enabeled compression for data.frame in h5write

USER VISIBLE CHANGES

- Use BiocStyles for package vignette

Risa
----

1.3.3: 1. Fixed definitions of assay.filenames.per.sample and factors.
        2. Fixed regulession of investigation file (i_) to be
        considered at the beginning of the string.  3. Added CITATION
        file.

RNASeqPower
-----------

Changes in version 1.0.1:

- Add a citation to the manual page, add CITATION file

rols
----

Changes in version 1.3.2:

- fixed bug in isIdObsolete, reported by Alex Thomas <2013-07-27 Sat>

Changes in version 1.3.1:

- other example in CVParam-class example, as PSI ontology is not
  available <2013-05-02 Thu>

Changes in version 1.3.0:

- Bioc 2.13 devel version bump

ROntoTools
----------

Changes in version 1.2.0:

- add the ability to analyze gene sets (pathways with no interaction)
  using only over-representation

- bug fixes: plot boundaries, etc.

rqubic
------

Changes in version 1.8.1 (2011-08-02):

- Fix errors in documentations

Changes in version 1.8.0 (2011-07-11):

- Submission to Bioconductor

Rsamtools
---------

Changes in version 1.14.0:

NEW FEATURES

- seqinfo(FaFile) returns available information on sequences and
  lengths on Fa (indexed fasta) files.

- filterBam accepts FilterRules call-backs for arbitrary filtering.

- add isIncomplete,BamFile-method to test for end-of-file

- add count.mapped.reads to summarizeOverlaps,*,BamFileList-method; set
  to TRUE to collect read and nucleotide counts via countBam.

- add summarizeOverlaps,*,character-method to count simple file paths

- add sequenceLayer() and stackStringsFromBam()

- add 'with.which_label' arg to readGAlignmentsFromBam(),
  readGappedReadsFromBam(), readGAlignmentPairsFromBam(), and
  readGAlignmentsListFromBam()

SIGNIFICANT USER-VISIBLE CHANGES

- rename: readBamGappedAlignments() -> readGAlignmentsFromBam()
  readBamGappedReads() -> readGappedReadsFromBam()
  readBamGappedAlignmentPairs() -> readGAlignmentPairsFromBam()
  readBamGAlignmentsList() -> readGAlignmentsListFromBam()
  makeGappedAlignmentPairs() -> makeGAlignmentPairs()

- speedup findMateAlignment()

DEPRECATED AND DEFUNCT

- deprecate readBamGappedAlignments(), readBamGappedReads(),
  readBamGappedAlignmentPairs(), readBamGAlignmentsList(), and
  makeGappedAlignmentPairs()

BUG FIXES

- scanVcfHeader tolerates records without ID fields, and with fields
  named similar to ID.

- close razip files only once.

- report tabix input errors

rSFFreader
----------


Rsubread
--------

Changes in version 1.12.0:

NEW FEATURES

- Added a number of new features to featureCounts read summarization
  function, including reordering of reads in BAM files to make reads
  from the same pair be adjacent to each other, support for chromosome
  aliases and output of complete annotation data for counting results
  from meta-feature level summarization.

NEW FEATURES

- It is described in more details in the Users Guide on how
  featureCounts program summarizes reads.

NEW FEATURES

- Improved short indel detection for both Subread (align function) and
  Subjunc aligners. This was achieved by building a consensus indel
  table and by realigning the reads. Discovered indels are reported in
  the output in addition to the read mapping results.

NEW FEATURES

- Support for detection of long indels (up to 200bp) was added in
  Subread. When the specified value of `indels' argument is greater
  than 16, Subread will automatically perform read assembly to detect
  long insertions and deletions.

NEW FEATURES

- Subread and Subjunc can now take FASTQ/FASTA, SAM and BAM files as
  input and output mapping results in both SAM and BAM formats.

NEW FEATURES

- Subjunc now directly operates on raw read data (it previously took
  Subread output as input), thus reducing running time by nearly half.

NEW FEATURES

- Subjunc can be instructed to output uniquely mapped reads. Hamming
  distance and mapping quality scores can be used to break ties when
  more than one best location was found.

NEW FEATURES

- More options were added to exactSNP program. Its documentation was
  also greatly improved.

NEW FEATURES

- A number of bug fixes.

rTANDEM
-------

Changes in version 1.1.6:

NEW FEATURES

- rTANDEM is now based on x!tandem version 13-09-01 'sledgehammer'

Changes in version 1.1.5:

BUG FIXES

- Fixed a memory management problen in Cpp

NEW FEATURES

- More parameters are supported.

Changes in version 1.1.4:

NEW FEATURES

- rTANDEM is now based on x!tandem version 13-06-15 'jackhammer'

BUG FIXES

- Corrected a bug that prevented the use of multi-threading

- Corrected a bug that could affect results when multiple fasta files
  are assigned to a same taxon.

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

RTN
---

Changes in version 1.0.0:

- 1st Bioconductor release of RTN [2013-10-15].

rtracklayer
-----------

Changes in version 1.22:

NEW FEATURES

- import,BigWigFile gains an asRle parameter that returns the data as
  an RleList (assuming it tiles the sequence); much faster than
  importing a GRanges and calling coverage() on it.

- add export,RleList,BigWigFile method for direct (and much more
  efficient) output of RleLists (like coverage) to BigWig files.

SIGNIFICANT USER-VISIBLE CHANGES

- UCSCData now extends GRanges instead of RangedData (thanks Herve)

BUG FIXES

- handle different Set-Cookie header field casing; often due to proxies
  (thanks to Altuna Akalin)

- attempt to more gracefully handle UCSC truncation of large data
  downloads

- handle re-directs due to UCSC mirroring (thanks Martin Morgan)

rTRM
----

Changes in version 1.0:

- First release in Bioconductor.

- Integration with Bioconductor packages: MotifDb and PWMEnrich, graph
  (limited).

- Human and mouse BioGRID datasets updated to version 3.2.105 (October
  2013).

rTRMui
------

Changes in version 1.0:

- First Bioconductor release.

- Provides a shiny interface to rTRM.

- All options available in rTRM can be accesses from rTRMui.

- Several utilities to control network appearance (node, edge and label
  size, etc.)

- Download results in PDF and text format.

SBMLR
-----

Changes in version 1.57.1:

NEW FEATURES

- Added new function S4toS3() for converting S4 SBML models of rsbml
  into S3 SBMLR models.

- Includes code from Vishak Venkateswaran's branch of SBMLR on github
  (July 2011). This is may allow it to read more models in
  http://www.ebi.ac.uk/biomodels-main/ I say may because I don't fully
  understand all of his codes, but add it in anyway.

- Problem of libsbml allowing multiple args to multiplication operator
  was solved by using prod(). Similarly for sum(). Note that "*"() is
  strictly a binary operator.

SIGNIFICANT USER-LEVEL CHANGES

- Call my model object call SBMLR now to let SBML refer to rsbml's SBML
  object. Similary my simulate() is now sim() to keep it clear of
  rsbml's simulate().

- SBMLR model files no longer need to have the reversible flag set. The
  default is FALSE, which is the opposite of Level 2: all of my
  reaction rate laws have always been positive, and a design objective
  I like is to keep SBMLR model files as lean as possible (subject to
  the constraint that they be valid R code).

- Simulate handles lsoda() event dataframes, see simulate help.

- curtoNatural.R (see Fig. 7 BMC Bioinformatics 2004, 5:190) is now in
  the models folder.

- The SOD model of Kowald et al, JTB 238 (2006) 828–840 is now also in
  this folder.

- Two pdfs of publications that are freely available were removed to
  make the package lighter.

- Similarly, only manual.doc remains: its redundant pdf is now out.

Notes

- SBMLR objects do need to have their reversible flag set in files, but
  do in SBMLR objects in R.

SCAN.UPC
--------

Changes in version 2.2.8:

NEW FEATURES

- Ability to download microarray data directly from Gene Expression
  Omnibus and normalize the files in a single command.

- Alternate functions (SCANfast and UPCfast) for performing SCAN and
  UPC normalization that use fewer probes and thus execute in ~60% less
  time.

- Ability to execute normalize multiple .CEL files in parallel either
  across multiple cores on a single computer or across multiple
  computers on a cluster.

- Ability to generate RNA-Seq annotation files to be used with the
  UPC_RNASeq function from GTF/FASTA source files.

- Ability to download and install BrainArray packages via an R
  function.

- Improved support for Affymetrix exon arrays.

- Improved support for Affymetrix HT_HG-U133A early access exon arrays.

OTHER

- Various improvements and clarifications to the docs.

SeqArray
--------

Changes in version 1.1.5:

- minor bug fix in asVCF

- update man page "SeqVarGDSClass-class.Rd" with new methods

- in DESCRIPTION, BiocGenerics listed in "Suggests" instead of
  "Imports" as suggested by R CMD check

- bug fix in seqDelete

- revise the function 'seqTranspose' according to the update of gdsfmt
  (v1.0.0)

Changes in version 1.1.4:

- add a new argument "action" to the function 'seqSetFilter'

- add a new function 'seqInfoNewVar' which allows adding new variables
  to the INFO fields

Changes in version 1.1.3:

- minor bug fix to avoid 'seqGetData' crashing when no value returned
  from a variable-length variable

- update documents

Changes in version 1.1.2:

- added methods 'qual', 'filt', 'asVCF'

- 'granges' method uses length of reference allele to set width

Changes in version 1.1.1:

- revise the argument 'var.index' in the function 'seqApply'

- basic supports of 'GRanges' and 'DNAStringSetList'

SeqGSEA
-------

Changes in version 1.1.5 (2013-09-01):

- Added an all-in function runSeqGSEA() for one step analysis

- Modified genpermuteMat() with invoking set.seed() for reproducibility

- Vignette updated

Changes in version 1.1.4 (2013-08-10):

- Updated functions to allow DE-only GSEA analysis.

- Fixed a few bugs.

- Vignette updated.

Changes in version 1.1.3 (2013-06-13):

- Sorted out an error due to biomaRt update (thanks to Martin Morgan).

Changes in version 1.1.2 (2013-05-01):

- Updated vignette.

Changes in version 1.1.1 (2013-04-23):

- The methodology paper of the SeqGSEA package published at BMC
  Bioinformatics (2013, 14(Suppl 5):S16).

- Added a function writeScores to generate DE/DS and gene scores
  output.

shinyTANDEM
-----------

Changes in version 0.99.1:

NEW FEATURES

- It is now possible to use a rTResult object from the R session as a
  parameter to the shinyTANDEM function.

Changes in version 0.99.0:

- First submission of shinyTANDEM to Bioconductor

ShortRead
---------

Changes in version 1.19:

SIGNIFICANT USER-VISIBLE CHANGES

- qa(..., type="fastq") uses a sample of n=1000000 reads by default,
  rather than then entire file; use sample=FALSE to revert to previous
  behavior.

NEW FEATURES

- encoding,FastqQuality and encoding,SFastqQuality provide a convenient
  map between letter encodings and their corresponding integer quality
  values.

- filterFastq transforms one fastq file to another, removing reads or
  nucleotides via a user-specified function. trimEnds,character-method
  & friends use this for an easy way to remove low-quality base.

BUG FIXES

- writeFastq successfully writes zero-length fastq files.

- FastqStreamer / FastqSampler warn on incomplete (corrupt) files

SpacePAC
--------

Changes in version 0.99.0:

- First release of the SpacePAC package.

- Two 3D clustering methods: Using a Simulation approach and using a
  Poisson approach.

- Allows a rudimentary plotting function to visualize the most
  significant cluster. Currently only one sphere can be plotted at a
  time.

spliceR
-------

Changes in version 1.0.0:

- Initial release.

supraHex
--------

Changes in version 1.0.0:

NEW FEATURES

- A pure R implementation of a self-organising learning algorithm as
  applied to the symmetric topology of the supra-hexagonal map

NEW FEATURES

- Dozens of functions for post-processing the trained map with multiple
  purposes, including: (i) visualizations of map indexes, hits and
  patterns; (ii) partitioning of the map into continuous clusters (i.e.
  gene meta-clusters) as they are different from gene clusters in an
  individual map node; and (iii) reordering of sample-specific map
  components (i.e. sample correlation)

NEW FEATURES

- Several omics datasets used to illustrate the functionalities
  supported in supraHex

synapter
--------

Changes in version 1.3.4:

- estimateMasterFdr now support list of vector as pepfiles <2013-09-27
  Fri>

- import(MSnbase) <2013-09-27 Fri>

Changes in version 1.3.3:

- Updated references. <2013-07-05 Fri>

Changes in version 1.3.2:

- fixed bug when using findEMRTs = "copy" <2013-05-30 Thu>

Changes in version 1.3.1:

- Reporting total number of peptides in dbUniquePeptideSet. Fixes issue
  #41. <2013-05-13 Mon>

- New mergedEMRTs arg in findEMRTs.  Closes issue #38. <2013-05-13 Mon>

- fixed synapterTiny$QuantPep3DData, which had the rownames as first
  column synapterTiny$QuantPep3DData$X. Detected thanks to new
  mergedEMRTs arg. <2013-05-13 Mon>

- added mergedEMRTs arg to synergise <2013-05-22 Wed>

- Synapter checks that one file per list element is passed <2013-05-28
  Tue>

- minor typo/fixes <2013-05-28 Tue>

- new idSource column when matching EMRTs <2013-05-29 Wed>

- new performance2 method that shows identification source and NA
  values contingency table <2013-05-29 Wed>

- new filterPeptideLength method <2013-05-29 Wed>

- added peplen argument to synergise to filter on peptide length
  <2013-05-29 Wed>

Changes in version 1.3.0:

- Bioc 2.13 devel version bump

TargetSearch
------------

Changes in version 1.18.0:

NEW FEATURES

- New options for function 'ProfileCleanUp' that allow fine tuning of
  the suggested metabolite in case of redundancy. This problem occurs
  when the reference library contains two or more metabolites in the
  same retention time window.

NEW FEATURES

- The function 'Write.Results' can create a quantification matrix based
  in one quantification mass. This mass can be selected automatically
  or specified by the user.

NEW FEATURES

- The above options have been incorporated in the GUI as well.

TCC
---

Changes in version 1.2.0:

NEW FEATURES

- This package was released as a Bioconductor package (previously
  CRAN).

- WAD method for identifying DEGs was added.

- ROKU method for identifying tissue-specific genes was added.

- 'increment' argument of 'calcNormFactor' function was added.

SIGNIFICANT USER-VISIBLE CHANGES

- 'replicates' field of TCC class was deleted.

Changes in version 1.1.3:

SIGNIFICANT USER-VISIBLE CHANGES

- 'generateSimulationData' function was renamed to 'simulateReadCount'.

SIGNIFICANT USER-VISIBLE CHANGES

- 'names' field of TCC class was changed to 'gene_id'.

SIGNIFICANT USER-VISIBLE CHANGES

- 'hypoData' was reduced to a smaller data set.

SIGNIFICANT USER-VISIBLE CHANGES

- 'hypoData_mg' was created. This is the simulation dataset which
  consists of 1,000 genes and 9 samples.

Changes in version 1.0.0:

SIGNIFICANT USER-VISIBLE CHANGES

- 'TCC' class was implemented as a R5 reference class. Wrapper
  functions with functional programming semantics were proviede.

TEQC
----

Changes in version 3.0.0:

- new function 'multiTEQCreport' collects results from 'TEQCreport'
  output from multiple samples and creates a combined quality report

Changes in version 2.9.2:

- bug fix in 'coverage.target': when getting reads data from a BAM
  file, output of average coverage values per target was wrong (in
  wrong order)

tigre
-----

Changes in version 1.14.1:

BUG FIXES

- Various fixes for non-cascade ODE model (i.e. model not using TF mRNA
  data) (thanks to Andrea Ocone for the report)

TransView
---------

Changes in version 1.5.9:

BUG FIXES

- Use of mcols in favor of values to access meta data in GRanges
  objects

Changes in version 1.5.8:

BUG FIXES

- Minor bug fix in annotatePeaks to accommodate recent changes in
  GRanges

Changes in version 1.5.7:

BUG FIXES

- Minor bug fix in annotatePeaks

Changes in version 1.5.6:

BUG FIXES

- Minor bug fix

Changes in version 1.5.5:

BUG FIXES

- Speed improvement and improved accuracy of annotatePeaks

Changes in version 1.5.4:

NEW FEATURES

- plotTV can now plot expression profiles including introns triggered
  by the new option pre_mRNA

BUG FIXES

- Fixed a bug which occurred if only two colors were passed to plotTV
  via colr

Changes in version 1.5.3:

NEW FEATURES

- meltPeak can return loess smoothed scores.

Changes in version 1.5.2:

BUG FIXES

- Fixed an issue with the gene names returned from gtf2gr.

Changes in version 1.5.1:

BUG FIXES

- The q-value cut off in macs2gr did not work as expected.

triplex
-------

Changes in version 1.2.0:

NEW FEATURES

- It's now possible to set custom scoring and isogroup tables through
  triplex.search interface.

NEW FEATURES

- New triplex.score.table and triplex.group.table functions for getting
  default scoring and isogroup tables to make the customization process
  almost effortless. OPTIMIZATIONS

NEW FEATURES

- Dynamic algorithm optimization technique implemented as a
  computational reduction based on minimal score option.

NEW FEATURES

- Optimized usage of processor data cache - the computation was divided
  into smaller pieces to prevent cache misses.

VariantAnnotation
-----------------

Changes in version 1.8.0:

NEW FEATURES

- Add 'upstream' and 'downstream' arguments to IntergenicVariants()
  constructor.

- Add 'samples' argument to ScanVcfParam().

- Add readGT(), readGeno() and readInfo().

- Add VRanges, VRangesList, SimpleVRangesList, and
  CompressedVRangesList classes.

- Add coercion VRanges -> VCF and VCF -> VRanges.

- Add methods for VRanges family: altDepth(), refDepth(), totalDepth(),
  altFraction() called(), hardFilters(), sampleNames(),
  softFilterMatrix() isIndel(), resetFilter().

- Add stackedSamples,VRangesList method.

MODIFICATIONS

- VCF validity method now requires the number of rows in info() to
  match the length of rowData().

- PRECEDEID and FOLLOWID from locateVariants() are now CharacterLists
  with all genes in 'upstream' and 'downstream' range.

- Modify rownames on rowData() GRanges to CHRAM:POS_REF/ALT for
  variants with no ID.

- readVcf() returns info() and geno() in the order specified in the
  ScanVcfParam.

- Work on scanVcf(): - free parse memory at first opportunity - define
  it_next in .c rather than .h - parse ALT "." in C - hash incoming
  strings - parse only param-requested 'fixed', 'info', 'geno' fields

- Add dimnames<-,VCF method to prevent 'fixed' fields from being copied
  into 'rowData' when new rownames or colnames were assigned.

- Support read/write for an emtpy VCF.

- readVcf(file=character, ...) method attempts coercion to TabixFile.

- Support for read/write an emtpy VCF.

- Add performance section to vignette; convert to BiocStyle.

- expand,CompressedVcf method expands geno() field 'AD' to length ALT +
  1. The expanded field is a (n x y x 2) array.

- 'genome' argument to readVcf() can be a character(1) or Seqinfo
  object.

DEPRECATED and DEFUNCT

- Defunct dbSNPFilter(), regionFilter() and MatrixToSnpMatrix().

- Deprecate readVcfLongForm().

BUG FIXES

- Fix bug in compatibility of read/writeVcf() when no INFO are columns
  present.

- Fix bug in locateVariants() when 'features' has no txid and cdsid.

- Fix bug in asVCF() when writing header lines.

- Fix bug in "expand" methods for VCF to handle multiple 'A' columns in
  info().

VariantTools
------------

Changes in version 1.4:

NEW FEATURES

- tallyVariants will now keep ref rows if variant_strand=0; this is
  useful for getting information when no alts are present (e.g., for
  making wildtype calls). Better have a big cluster to do this over the
  whole genome.

- add a keep_extra_stats param to TallyVariantsParam; setting this to
  FALSE will speed things up when the extra stats are not needed.

- idVerify now supports VCF input like that output by GATK.

- callableFraction() now supports GRangesList and TranscriptDb.

USER-VISIBLE CHANGES

- The API is now based on VRanges, a formal GRanges-derived class for
  representing variants; use of so-called "tally" or "variant" GRanges
  is deprecated.

- Disable proximity filter by default; we recommend this now only for
  whole genome calling.

- QA filtering is no longer a formal part of the calling pipeline; we
  recommend to apply QA filters "softly" via qaVariants() and use the
  results for diagnostics only.

- Use BiocParallel (BPPARAM argument) for tallyVariants

- VariantTallyParam deprecated; use TallyVariantsParam

BUG FIXES

- idVerify now correctly computes cliques instead of connected
  components

- use the total count, rather than the ref count when calculating the
  alt frequency

xcms
----

Changes in version 1.37.6:

NEW FEATURE

- Introducing write.mzQuantML(xcmsSet) to export the peak list and
  grouped matrix to the PSI format mzQuantML (see
  http://www.psidev.info/mzquantml)

USER VISIBLE CHANGES

- Add Brigham Young University to LICENSE file for copyright purposes.

- Add copyright information display when running
  findPeaks.massifquant() within xcmsRaw.R

- Clean and update documentation for findPeaks.massifquant-methods.Rd

BUG FIXES

- Remove unused parameters in findKalmanROIs() within xcmsRaw.R

Changes in version 1.37.5:

BUG FIXES

- fixed bug in retcor.obiwarp where the scanrange of the first sample
  would be checked instead of the center sample

Changes in version 1.37.4:

BUG FIXES

- Skip t-test in diffreport() if one class has less than 2 samples.

Changes in version 1.37.3:

BUG FIXES

- fixed bug in patternVsRowScore (group.nearest) that was introduced by
  the modifications in rev 65169 and caused features to be aligned that
  were far outside the given m/z and retention time windows.

Changes in version 1.37.1:

BUG FIXES

- fixed fillPeaks, which 1) dropped non-standard columns and 2) failed
  if nothing to do, based on patches by Tony Larson.

NEW FEATURES

- Introducing msn2xcmsRaw, to allow findPeaks() on MS2 and MSn data

xps
---

Changes in version 3.00:

VERSION xps-1.21.5

- add QualTreeSet methods NUSE() and RLE() to get stats and values

- update man export.Rd

VERSION xps-1.21.4

- update xpsQAReport.R for R-3.x

VERSION xps-1.21.3

- update XPSSchemes.cxx to replace error with warning for missing
  annotation header '%netaffx-annotation-'

VERSION xps-1.21.2

- update XPSNormalizer.cxx to correct for uninitialized variables

VERSION xps-1.21.1

- update README

- update Makefile to set include path (for ~/.R/Makevars)

- update XPSUtils.cxx to eliminate warning with clang

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

Changes in version 2.15:

VERSION xps-1.17.2

- update script4schemes.R to include schemes for HuGene-2_x-st arrays

- update script4xps.R with example 4a for HuGene-2_1-st arrays

- update README

VERSION xps-1.17.1

- remove warnings: partial argument match

- zzz.R: use .onAttach()

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


Packages removed from the release
=================================

The following packages are no longer in the release:

dualKS, externalVector, GeneGroupAnalysis, iFlow, KEGGSOAP, xmapcore
